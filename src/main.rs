use bevy::prelude::*;

// Constants
const PARTICLE_COUNT: usize = 1000;
const PARTICLE_RADIUS: f32 = 2.0;
const GRAVITY: Vec2 = Vec2::new(0.0, -98.0);
const RESTITUTION: f32 = 0.9;

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
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .insert_resource(SimParams {
            gravity: GRAVITY,
            dt: 1.0/60.0,
            restitution: RESTITUTION,
        })
        .add_systems(Startup, (setup, spawn_particles))
        .add_systems(Update, (update_particles, render_particles))
        .run();
}

// Setup camera
fn setup(mut commands: Commands) {
    commands.spawn(Camera2d::default());
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