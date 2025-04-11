use bevy::prelude::*;

// Core Components
#[derive(Component, Default, Debug, Clone, Copy)]
pub struct Position(pub Vec2);

#[derive(Component, Default, Debug, Clone, Copy)]
pub struct Velocity(pub Vec2);

#[derive(Component, Default, Debug, Clone, Copy)]
pub struct Mass(pub f32);

#[derive(Component, Default, Debug, Clone, Copy)]
pub struct AffineMomentum(pub Mat2);

#[derive(Component, Default, Debug, Clone)]
pub struct ParticleTag;

#[derive(Component, Default, Debug, Clone)]
pub struct CellMassMomentumContributions(pub [GridMassAndMomentumChange; 9]);

#[derive(Default, Debug, Clone, Copy)]
pub struct GridMassAndMomentumChange(pub usize, pub f32, pub Vec2);

// Grid Cell data
#[derive(Default, Debug, Clone, Copy)]
pub struct GridCell {
    pub mass: f32,
    pub velocity: Vec2,
}

// Grid Resource
#[derive(Resource, Default, Debug, Clone)]
pub struct Grid {
    pub width: usize,
    pub height: usize,
    pub cells: Vec<GridCell>,
}

impl Grid {
    pub fn new(width: usize, height: usize) -> Self {
        let num_cells = width * height;
        let mut cells = Vec::with_capacity(num_cells);
        cells.resize_with(num_cells, Default::default);

        Self {
            width,
            height,
            cells,
        }
    }

    pub fn index_at(&self, x: usize, y: usize) -> usize {
        // Clamp to grid boundaries
        let x = x.clamp(0, self.width - 1);
        let y = y.clamp(0, self.height - 1);
        y * self.width + x
    }

    pub fn update(&mut self, dt: f32, gravity: f32) {
        // Update grid velocities based on forces (like gravity)
        for cell in &mut self.cells {
            if cell.mass > 0.0 {
                // Apply gravity
                cell.velocity.y += gravity * dt;

                // Normalize velocity by mass (needed after P2G transfer)
                cell.velocity /= cell.mass;
            }
        }
    }

    pub fn reset(&mut self) {
        // Reset grid cells for the next simulation step
        for cell in &mut self.cells {
            cell.mass = 0.0;
            cell.velocity = Vec2::ZERO;
        }
    }
}

// World State Parameters
#[derive(Resource, Debug, Clone)]
pub struct WorldState {
    pub dt: f32,
    pub gravity: f32,
    pub gravity_enabled: bool,
    pub time: f32,
    pub iteration: usize,
}

impl Default for WorldState {
    fn default() -> Self {
        Self {
            dt: 1.0 / 60.0,
            gravity: -9.81,
            gravity_enabled: true,
            time: 0.0,
            iteration: 0,
        }
    }
}

impl WorldState {
    pub fn update(&mut self) {
        self.time += self.dt;
        self.iteration += 1;
    }
}

// Interpolation functions for MPM
pub fn quadratic_interpolation_weights(diff: Vec2) -> [Vec2; 3] {
    let mut weights = [Vec2::ZERO; 3];

    // Quadratic interpolation weights - using quadratic B-splines
    // For x-direction
    let x = diff.x;
    weights[0].x = 0.5 * (0.5 - x) * (0.5 - x);
    weights[1].x = 0.75 - x * x;
    weights[2].x = 0.5 * (0.5 + x) * (0.5 + x);

    // For y-direction
    let y = diff.y;
    weights[0].y = 0.5 * (0.5 - y) * (0.5 - y);
    weights[1].y = 0.75 - y * y;
    weights[2].y = 0.5 * (0.5 + y) * (0.5 + y);

    weights
}

// Helper for APIC calculation
pub fn weighted_velocity_and_cell_dist_to_term(weighted_velocity: Vec2, cell_dist: Vec2) -> Mat2 {
    Mat2::from_cols(
        weighted_velocity * cell_dist.x,
        weighted_velocity * cell_dist.y,
    )
}

// MPM Algorithm systems

// P2G (Particles to Grid) step
pub fn particles_to_grid(
    mut grid: ResMut<Grid>,
    particles: Query<(&Position, &Velocity, &Mass), With<ParticleTag>>,
) {
    // Reset grid before transferring particle data
    grid.reset();

    for (position, velocity, mass) in particles.iter() {
        let cell_x: u32 = position.0.x as u32;
        let cell_y: u32 = position.0.y as u32;
        let cell_diff = Vec2::new(
            position.0.x - cell_x as f32 - 0.5,
            position.0.y - cell_y as f32 - 0.5,
        );
        let weights = quadratic_interpolation_weights(cell_diff);

        // Transfer mass and momentum to surrounding 9 cells
        for gx in 0..3 {
            for gy in 0..3 {
                let weight = weights[gx].x * weights[gy].y;
                let cell_pos_x = (cell_x as i32 + gx as i32) - 1;
                let cell_pos_y = (cell_y as i32 + gy as i32) - 1;
                let cell_index = grid.index_at(cell_pos_x as usize, cell_pos_y as usize);

                // PIC transfer (mass and velocity)
                grid.cells[cell_index].mass += weight * mass.0;
                grid.cells[cell_index].velocity += weight * velocity.0 * mass.0;
            }
        }
    }
}

// Grid update step
pub fn update_grid(mut grid: ResMut<Grid>, mut world: ResMut<WorldState>) {
    world.update();
    grid.update(
        world.dt,
        if world.gravity_enabled {
            world.gravity
        } else {
            0.0
        },
    );
}

// G2P (Grid to Particles) step
pub fn grid_to_particles(
    world: Res<WorldState>,
    grid: Res<Grid>,
    mut particles: Query<(&mut Position, &mut Velocity), With<ParticleTag>>,
) {
    for (mut position, mut velocity) in particles.iter_mut() {
        // Reset particle velocity - will be recalculated from grid
        velocity.0 = Vec2::ZERO;

        let cell_x: u32 = position.0.x as u32;
        let cell_y: u32 = position.0.y as u32;
        let cell_diff = Vec2::new(
            position.0.x - cell_x as f32 - 0.5,
            position.0.y - cell_y as f32 - 0.5,
        );
        let weights = quadratic_interpolation_weights(cell_diff);

        // Gather velocity from surrounding cells
        for gx in 0..3 {
            for gy in 0..3 {
                let weight = weights[gx].x * weights[gy].y;
                let cell_pos_x = (cell_x as i32 + gx as i32) - 1;
                let cell_pos_y = (cell_y as i32 + gy as i32) - 1;
                let cell_index = grid.index_at(cell_pos_x as usize, cell_pos_y as usize);

                // Transfer velocity from grid to particle
                velocity.0 += grid.cells[cell_index].velocity * weight;
            }
        }

        // Advect particles using velocity
        position.0 += velocity.0 * world.dt;

        // Boundary conditions - keep particles inside grid
        position.0.x = position.0.x.clamp(1.0, (grid.width - 2) as f32);
        position.0.y = position.0.y.clamp(1.0, (grid.height - 2) as f32);
    }
}

// Spawn particles in a square region
pub fn spawn_particles(
    commands: &mut Commands,
    start_x: f32,
    start_y: f32,
    width: f32,
    height: f32,
    spacing: f32,
    mass: f32,
) {
    let cols = (width / spacing) as usize;
    let rows = (height / spacing) as usize;

    for i in 0..rows {
        for j in 0..cols {
            let position = Vec2::new(start_x + j as f32 * spacing, start_y + i as f32 * spacing);

            commands.spawn((
                Position(position),
                Velocity(Vec2::ZERO),
                Mass(mass),
                ParticleTag,
                CellMassMomentumContributions([GridMassAndMomentumChange(0, 0.0, Vec2::ZERO); 9]),
                // Add these components to fix hierarchy warnings
                Transform::default(),
                GlobalTransform::default(),
                Visibility::default(),
                InheritedVisibility::default(),
            ));
        }
    }
}

// Plugin for MPM simulation
pub struct MPMPlugin;

impl Plugin for MPMPlugin {
    fn build(&self, app: &mut App) {
        app.insert_resource(Grid::new(100, 100)) // Larger grid for more space
            .insert_resource(WorldState {
                dt: 1.0 / 120.0, // Higher framerate for stability
                gravity: -5.0,   // Less extreme gravity
                gravity_enabled: true,
                time: 0.0,
                iteration: 0,
            })
            .add_systems(
                Update,
                (
                    particles_to_grid,
                    update_grid.after(particles_to_grid),
                    grid_to_particles.after(update_grid),
                )
                    .chain(),
            );
    }
}

// Visual representation of particles
#[derive(Component)]
struct ParticleSprite;

// Setup system
fn setup(mut commands: Commands) {
    // Add a simple 2D camera with some zoom
    commands.spawn(Camera2d::default());

    // Create a dam break scenario - a tall column of particles on the left side
    spawn_particles(
        &mut commands,
        5.0,  // start_x - near the left edge
        5.0,  // start_y - from bottom
        15.0, // width - compact block
        50.0, // height - tall column
        0.8,  // spacing - slightly denser
        1.0,  // mass
    );
}

// Render particles
fn render_particles(
    mut commands: Commands,
    particles: Query<(Entity, &Position), With<ParticleTag>>,
    mut existing_sprites: Query<(&mut Transform, &mut Sprite, &ParticleSprite)>,
) {
    // Create sprites if we don't have enough
    let particle_count = particles.iter().count();
    let sprite_count = existing_sprites.iter().count();

    // Spawn more sprites if needed
    if sprite_count < particle_count {
        for _ in 0..(particle_count - sprite_count) {
            commands.spawn((
                Sprite {
                    color: Color::rgb(0.2, 0.6, 0.9), // Water blue color
                    custom_size: Some(Vec2::new(1.0, 1.0)),
                    ..Default::default()
                },
                Transform::default(),
                ParticleSprite,
            ));
        }
    }

    // Update sprite positions based on particles
    let mut sprites = existing_sprites.iter_mut();
    for (_, position) in &particles {
        if let Some((mut transform, _, _)) = sprites.next() {
            transform.translation.x = position.0.x;
            transform.translation.y = position.0.y;
            transform.translation.z = 1.0; // Make sure it's visible
        }
    }
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugins(MPMPlugin)
        .add_systems(Startup, setup)
        .add_systems(Update, render_particles)
        .run();
}
