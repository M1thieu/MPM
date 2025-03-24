use bevy::prelude::*;

// Constants
const PARTICLE_COUNT: usize = 1000;
const PARTICLE_RADIUS: f32 = 2.0;
const GRAVITY: Vec2 = Vec2::new(0.0, -98.0);
const RESTITUTION: f32 = 0.9;
const GRID_CELL_SIZE: f32 = 10.0;
const GRID_MARGIN: usize = 1;

// Component for particles
#[derive(Component)]
struct Particle {
    position: Vec2,
    velocity: Vec2,
    mass: f32,
}

// Resource for simulation parameters
#[derive(Resource)]
struct SimParams {
    gravity: Vec2,
    dt: f32,
    restitution: f32,
    show_grid: bool,
}

// Grid resource for MPM/PIC method
#[derive(Resource)]
struct Grid {
    cell_size: f32,
    dimensions: UVec2,
    world_bounds: Vec2,
    origin: Vec2,
    cell_mass: Vec<f32>,
    cell_velocity: Vec<Vec2>,
    cell_momentum: Vec<Vec2>,
}

impl Grid {
    fn world_to_grid(&self, world_pos: Vec2) -> UVec2 {
        let rel_pos = world_pos - self.origin;
        UVec2::new(
            (rel_pos.x / self.cell_size).floor() as u32,
            (rel_pos.y / self.cell_size).floor() as u32
        )
    }
    
    fn grid_to_world(&self, grid_pos: UVec2) -> Vec2 {
        self.origin + Vec2::new(
            (grid_pos.x as f32 + 0.5) * self.cell_size,
            (grid_pos.y as f32 + 0.5) * self.cell_size
        )
    }
    
    fn in_bounds(&self, grid_pos: UVec2) -> bool {
        grid_pos.x < self.dimensions.x && grid_pos.y < self.dimensions.y
    }
    
    fn get_index(&self, grid_pos: UVec2) -> usize {
        (grid_pos.y * self.dimensions.x + grid_pos.x) as usize
    }
    
    fn reset(&mut self) {
        for i in 0..self.cell_mass.len() {
            self.cell_mass[i] = 0.0;
            self.cell_momentum[i] = Vec2::ZERO;
            self.cell_velocity[i] = Vec2::ZERO;
        }
    }
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .insert_resource(SimParams {
            gravity: GRAVITY,
            dt: 1.0/60.0,
            restitution: RESTITUTION,
            show_grid: true,
        })
        .add_systems(Startup, (setup, initialize_grid, spawn_particles))
        .add_systems(Update, (
            update_particles, 
            render_particles,
            render_grid
        ))
        .run();
}

// Setup camera
fn setup(mut commands: Commands) {
    commands.spawn(Camera2d::default());
}

// Initialize the grid
fn initialize_grid(mut commands: Commands, windows: Query<&Window>) {
    let window = windows.single();
    let cells_x = (window.width() / GRID_CELL_SIZE).ceil() as u32 + (GRID_MARGIN as u32 * 2);
    let cells_y = (window.height() / GRID_CELL_SIZE).ceil() as u32 + (GRID_MARGIN as u32 * 2);
    let origin = Vec2::new(
        -window.width()/2.0 - (GRID_MARGIN as f32 * GRID_CELL_SIZE),
        -window.height()/2.0 - (GRID_MARGIN as f32 * GRID_CELL_SIZE)
    );
    let total_cells = (cells_x * cells_y) as usize;
    
    commands.insert_resource(Grid {
        cell_size: GRID_CELL_SIZE,
        dimensions: UVec2::new(cells_x, cells_y),
        world_bounds: Vec2::new(window.width(), window.height()),
        origin,
        cell_mass: vec![0.0; total_cells],
        cell_velocity: vec![Vec2::ZERO; total_cells],
        cell_momentum: vec![Vec2::ZERO; total_cells],
    });
}

// Spawn initial particles
fn spawn_particles(mut commands: Commands, windows: Query<&Window>) {
    let window = windows.single();
    let width = window.width();
    let height = window.height();
    let cols = (PARTICLE_COUNT as f32).sqrt().ceil() as usize;
    let rows = (PARTICLE_COUNT + cols - 1) / cols;
    let spacing = Vec2::new(width / (cols as f32 + 1.0), height / (rows as f32 + 1.0));
    
    for i in 0..PARTICLE_COUNT {
        let pos = Vec2::new(
            ((i % cols) as f32 + 0.5) * spacing.x - width / 2.0,
            ((i / cols) as f32 + 0.5) * spacing.y - height / 2.0);
        let dir = pos.normalize_or_zero();
        
        commands.spawn(Particle { 
            position: pos,
            velocity: Vec2::new(dir.y * 20.0, -dir.x * 20.0),
            mass: 1.0 
        });
    }
}

// Simple particle update with boundary collisions
fn update_particles(mut particles: Query<&mut Particle>, windows: Query<&Window>, params: Res<SimParams>) {
    let bounds = Vec2::new(windows.single().width() / 2.0, windows.single().height() / 2.0);
    
    for mut particle in &mut particles {
        // Apply gravity
        particle.velocity += params.gravity * params.dt;
        
        // Update position
        let velocity = particle.velocity; // Store locally to avoid borrow issues
        particle.position += velocity * params.dt;
        
        // Boundary collisions (X-axis)
        if particle.position.x < -bounds.x + PARTICLE_RADIUS {
            particle.position.x = -bounds.x + PARTICLE_RADIUS;
            particle.velocity.x *= -params.restitution;
        } else if particle.position.x > bounds.x - PARTICLE_RADIUS {
            particle.position.x = bounds.x - PARTICLE_RADIUS;
            particle.velocity.x *= -params.restitution;
        }
        
        // Boundary collisions (Y-axis)
        if particle.position.y < -bounds.y + PARTICLE_RADIUS {
            particle.position.y = -bounds.y + PARTICLE_RADIUS;
            particle.velocity.y *= -params.restitution;
        } else if particle.position.y > bounds.y - PARTICLE_RADIUS {
            particle.position.y = bounds.y - PARTICLE_RADIUS;
            particle.velocity.y *= -params.restitution;
        }
    }
}

// Render particles
fn render_particles(mut gizmos: Gizmos, particles: Query<&Particle>) {
    for particle in &particles {
        gizmos.circle_2d(particle.position, PARTICLE_RADIUS, Color::WHITE);
    }
}

// Render the grid
fn render_grid(mut gizmos: Gizmos, grid: Res<Grid>, params: Res<SimParams>) {
    if !params.show_grid { return; }
    
    let color = Color::rgba(0.3, 0.3, 0.8, 0.2);
    
    // Draw horizontal grid lines
    for y in 0..=grid.dimensions.y {
        let y_pos = grid.origin.y + y as f32 * grid.cell_size;
        gizmos.line_2d(
            Vec2::new(grid.origin.x, y_pos),
            Vec2::new(grid.origin.x + grid.dimensions.x as f32 * grid.cell_size, y_pos),
            color
        );
    }
    
    // Draw vertical grid lines
    for x in 0..=grid.dimensions.x {
        let x_pos = grid.origin.x + x as f32 * grid.cell_size;
        gizmos.line_2d(
            Vec2::new(x_pos, grid.origin.y),
            Vec2::new(x_pos, grid.origin.y + grid.dimensions.y as f32 * grid.cell_size),
            color
        );
    }
}