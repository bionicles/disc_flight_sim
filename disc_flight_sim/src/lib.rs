//! Simulate a disc golf flight in 3D with RK4 integration.
#![allow(clippy::multiple_crate_versions)]
use serde::{Deserialize, Serialize};
use uom::si::{
    acceleration::meter_per_second_squared,
    angle::radian,
    angular_acceleration::radian_per_second_squared,
    angular_velocity::radian_per_second,
    area::square_meter,
    f32::{
        Acceleration, Angle, AngularAcceleration, AngularVelocity, Area, Force, Length, Mass,
        MassDensity, Ratio, Time, Velocity,
    },
    force::newton,
    length::meter,
    mass::kilogram,
    mass_density::kilogram_per_cubic_meter,
    ratio::ratio,
    time::second,
    velocity::meter_per_second,
};
use wasm_bindgen::prelude::*;

/// Parameters describing a disc golf disc's physical and flight characteristics.
#[derive(Serialize, Deserialize, Clone, Copy, Debug)]
pub struct Disc {
    /// Initial launch speed of the disc.
    pub v0: Velocity, // m/s

    /// Launch angle of the disc above the horizontal plane.
    pub angle: Angle, // radians

    /// Initial angular velocity (spin rate) of the disc.
    pub spin: AngularVelocity, // rad/s

    /// Initial azimuthal angle (spin orientation).
    pub phi: Angle, // radians

    /// Mass of the disc.
    pub mass: Mass, // kg

    /// Cross-sectional area of the disc.
    pub area: Area, // m²

    /// Dimensionless drag coefficient (`Cd`).
    pub cd: Ratio, // unitless

    /// Dimensionless lift coefficient (`Cl`).
    pub cl: Ratio, // unitless

    /// Gyroscopic stability factor (dimensionless).
    pub gyro: Ratio, // unitless

    /// Total upper surface area (for advanced lift calculations).
    pub surface: Area, // m²
}

impl Default for Disc {
    fn default() -> Self {
        Self {
            v0: Velocity::new::<meter_per_second>(20.0),
            angle: Angle::new::<radian>(10.0_f32.to_radians()),
            spin: AngularVelocity::new::<radian_per_second>(40.0),
            phi: Angle::new::<radian>(0.0),
            mass: Mass::new::<kilogram>(0.175),
            area: Area::new::<square_meter>(0.056),
            cd: Ratio::new::<ratio>(0.6),
            cl: Ratio::new::<ratio>(0.7),
            gyro: Ratio::new::<ratio>(1.0),
            surface: Area::new::<square_meter>(0.057),
        }
    }
}

/// Environmental and simulation parameters for disc flight.
#[derive(Serialize, Deserialize, Clone, Copy, Debug)]
pub struct EnvParams {
    /// Air density.
    pub rho: MassDensity, // kg/m³

    /// Acceleration due to gravity.
    pub g: Acceleration, // m/s²

    /// Wind velocity along the X axis.
    pub wind_x: Velocity, // m/s

    /// Wind velocity along the Y axis.
    pub wind_y: Velocity, // m/s

    /// Wind velocity along the Z axis.
    pub wind_z: Velocity, // m/s

    /// Integration timestep.
    pub dt: Time, // s

    /// Maximum simulation time.
    pub max_t: Time, // s
}

impl Default for EnvParams {
    fn default() -> Self {
        Self {
            rho: MassDensity::new::<kilogram_per_cubic_meter>(1.225),
            g: Acceleration::new::<meter_per_second_squared>(9.81),
            wind_x: Velocity::new::<meter_per_second>(0.0),
            wind_y: Velocity::new::<meter_per_second>(0.0),
            wind_z: Velocity::new::<meter_per_second>(0.0),
            dt: Time::new::<second>(0.02),
            max_t: Time::new::<second>(10.0),
        }
    }
}

/// A single point in the simulated trajectory of a disc.
#[derive(Serialize, Deserialize, Clone, Copy, Debug)]
pub struct TrajectoryPoint {
    /// X coordinate (meters) at this time point.
    pub x: Length,
    /// Y coordinate (meters) at this time point.
    pub y: Length,
    /// Z coordinate (meters) (height above ground) at this time point.
    pub z: Length,
    /// Time (seconds) since the start of the simulation.
    pub t: Time,
}

/// Simulate a disc golf flight in 3D with RK4 integration.
///
/// # Errors
/// Returns a `JsValue` error if serialization or input parameters are invalid.
///
/// # Parameters
/// - `disc`: `JsValue` representing the disc parameters.
/// - `env`: `JsValue` representing the environment parameters.
///
/// # Returns
/// - `JsValue`: A JSON array of `TrajectoryPoints` on success.
#[must_use]
#[wasm_bindgen]
pub fn simulate_flight(disc: JsValue, env: JsValue) -> JsValue {
    let disc: Disc = match serde_wasm_bindgen::from_value(disc) {
        Ok(d) => d,
        Err(e) => return JsValue::from_str(&format!("invalid disc params: {e}")),
    };
    let env: EnvParams = match serde_wasm_bindgen::from_value(env) {
        Ok(e) => e,
        Err(e) => return JsValue::from_str(&format!("invalid env params: {e}")),
    };
    let traj = simulate(&disc, &env);
    match serde_wasm_bindgen::to_value(&traj) {
        Ok(val) => val,
        Err(e) => JsValue::from_str(&format!("serialize error: {e}")),
    }
}

/// The complete physical state of a disc in flight at a single instant in time.
///
/// This struct holds all the fundamental dynamic variables required to
/// fully describe the disc’s state for integration purposes.
///
/// # Fields
/// - `x`, `y`, `z`: Position coordinates in meters.
/// - `vx`, `vy`, `vz`: Velocity components in meters per second.
/// - `spin`: Spin rate (angular velocity) in radians per second.
#[derive(Copy, Clone, Debug)]
pub struct PhysicsState {
    /// X position (meters).
    pub x: Length,
    /// Y position (meters).
    pub y: Length,
    /// Z position (meters).
    pub z: Length,
    /// Velocity in X direction (meters/second).
    pub vx: Velocity,
    /// Velocity in Y direction (meters/second).
    pub vy: Velocity,
    /// Velocity in Z direction (meters/second).
    pub vz: Velocity,
    /// Spin rate (radians/second).
    pub spin: AngularVelocity,
}

/// The time derivative of a `PhysicsState`.
///
/// Represents the rate of change (first derivative) of every field in
/// `PhysicsState`. Used for ODE integration (e.g., RK4).
///
/// # Fields
/// - `dx`, `dy`, `dz`: Time derivatives of position (i.e., velocity).
/// - `dvx`, `dvy`, `dvz`: Time derivatives of velocity (i.e., acceleration).
/// - `dspin`: Time derivative of spin rate (i.e., angular acceleration).
pub struct PhysicsDerivative {
    /// Time derivative of X position (velocity in X, meters/second).
    pub dx: Velocity,
    /// Time derivative of Y position (velocity in Y, meters/second).
    pub dy: Velocity,
    /// Time derivative of Z position (velocity in Z, meters/second).
    pub dz: Velocity,
    /// Time derivative of X velocity (acceleration in X, meters/second²).
    pub dvx: Acceleration,
    /// Time derivative of Y velocity (acceleration in Y, meters/second²).
    pub dvy: Acceleration,
    /// Time derivative of Z velocity (acceleration in Z, meters/second²).
    pub dvz: Acceleration,
    /// Time derivative of spin (angular acceleration, radians/second²).
    pub dspin: AngularAcceleration,
}

/// Classic RK4 ODE solver for fixed step
pub fn rk4_step<F>(y: &PhysicsState, t: Time, dt: Time, f: F) -> PhysicsState
where
    F: Fn(&PhysicsState, Time) -> PhysicsDerivative,
{
    let k1 = f(y, t);
    let y2 = y.add_scaled(&k1, dt / 2.0);
    let k2 = f(&y2, t + dt / 2.0);
    let y3 = y.add_scaled(&k2, dt / 2.0);
    let k3 = f(&y3, t + dt / 2.0);
    let y4 = y.add_scaled(&k3, dt);
    let k4 = f(&y4, t + dt);

    y.linear_combination(&k1, &k2, &k3, &k4, dt)
}

impl PhysicsState {
    /// Returns a new `PhysicsState` computed as `self + scale * derivative` for each field.
    #[must_use]
    pub fn add_scaled(&self, k: &PhysicsDerivative, scale: Time) -> Self {
        Self {
            x: self.x + k.dx * scale,
            y: self.y + k.dy * scale,
            z: self.z + k.dz * scale,
            vx: self.vx + k.dvx * scale,
            vy: self.vy + k.dvy * scale,
            vz: self.vz + k.dvz * scale,
            spin: self.spin + AngularVelocity::new::<radian_per_second>((k.dspin * scale).value),
            //  (radians / second²) × (seconds) = radians / second      ^^^^^^^^^^^^^^^^^
        }
    }

    /// RK4 update: combine weighted derivatives to advance the state.
    #[must_use]
    pub fn linear_combination(
        &self,
        k1: &PhysicsDerivative,
        k2: &PhysicsDerivative,
        k3: &PhysicsDerivative,
        k4: &PhysicsDerivative,
        dt: Time,
    ) -> Self {
        let sixth = dt / 6.0;
        let spin_increment_value =
            (k1.dspin + 2.0 * k2.dspin + 2.0 * k3.dspin + k4.dspin).value * sixth.value;
        Self {
            x: self.x + sixth * (k1.dx + 2.0 * k2.dx + 2.0 * k3.dx + k4.dx),
            y: self.y + sixth * (k1.dy + 2.0 * k2.dy + 2.0 * k3.dy + k4.dy),
            z: self.z + sixth * (k1.dz + 2.0 * k2.dz + 2.0 * k3.dz + k4.dz),
            vx: self.vx + sixth * (k1.dvx + 2.0 * k2.dvx + 2.0 * k3.dvx + k4.dvx),
            vy: self.vy + sixth * (k1.dvy + 2.0 * k2.dvy + 2.0 * k3.dvy + k4.dvy),
            vz: self.vz + sixth * (k1.dvz + 2.0 * k2.dvz + 2.0 * k3.dvz + k4.dvz),
            spin: AngularVelocity::new::<radian_per_second>(
                self.spin.get::<radian_per_second>() + spin_increment_value,
            ),
        }
    }
}

/// Simulate the 3D flight of a disc golf disc under specified environmental conditions.
///
/// Uses a fixed-step fourth-order Runge-Kutta (RK4) integration method to numerically solve
/// the disc's equations of motion. The simulation proceeds in discrete time steps from the initial
/// conditions until either the maximum simulation time is reached or the disc hits the ground (z < 0).
///
/// # Parameters
/// - `disc`: Reference to a [`Disc`] struct containing the disc's initial state and aerodynamic parameters.
/// - `env`: Reference to an [`EnvParams`] struct specifying environmental parameters such as air density, gravity, wind, and integration step size.
///
/// # Returns
/// A vector of [`TrajectoryPoint`] structs representing the disc's position and time at each simulation step.
///
/// # Behavior
/// - The disc is launched from position (0, 0, 1.5) meters (1.5m above ground).
/// - The initial velocity and spin are set according to `disc`.
/// - The simulation ends when the disc's z-coordinate is less than 0 (touches ground) or `env.max_t` is reached.
///
/// # Panics
/// This function will not panic under normal circumstances. However, unphysical input (e.g., negative mass)
/// may yield invalid results.
///
/// # Example
/// ```rust
/// use disc_flight_sim::{Disc, EnvParams, simulate};
/// let disc = Disc::default();
/// let env = EnvParams::default();
/// let path = simulate(&disc, &env);
/// # assert!(!path.is_empty());
/// ```
#[must_use]
pub fn simulate(disc: &Disc, env: &EnvParams) -> Vec<TrajectoryPoint> {
    // State: [x, y, z, vx, vy, vz, spin]
    let mut t = Time::new::<second>(0.0);
    let mut state = PhysicsState {
        x: Length::new::<meter>(0.0),
        y: Length::new::<meter>(0.0),
        z: Length::new::<meter>(1.5),
        vx: disc.v0 * disc.angle.get::<radian>().cos(),
        vy: disc.v0 * disc.angle.get::<radian>().sin(),
        vz: disc.v0 * disc.phi.get::<radian>().sin(),
        spin: disc.spin,
    };
    let mut prev_state;
    let mut prev_t;
    let mut path = Vec::new();
    while t < env.max_t && state.z.get::<meter>() >= 0.0 {
        prev_state = state;
        prev_t = t;
        path.push(TrajectoryPoint {
            x: state.x,
            y: state.y,
            z: state.z,
            t,
        });
        state = rk4_step(
            &state,
            t,
            env.dt,
            // formatting comment
            |s, _t| dynamics(s, disc, env),
        );
        t += env.dt;

        // Check if we've crossed below zero
        if state.z.get::<meter>() < 0.0 {
            // Interpolate for landing point
            let z1 = prev_state.z.get::<meter>();
            let z2 = state.z.get::<meter>();
            let frac = z1 / (z1 - z2); // How far between prev and current?
            let interp = |a: f32, b: f32| frac.mul_add(b - a, a);

            path.push(TrajectoryPoint {
                x: Length::new::<meter>(interp(
                    prev_state.x.get::<meter>(),
                    state.x.get::<meter>(),
                )),
                y: Length::new::<meter>(interp(
                    prev_state.y.get::<meter>(),
                    state.y.get::<meter>(),
                )),
                z: Length::new::<meter>(0.0),
                t: prev_t + (t - prev_t) * frac, // estimate time at ground hit
            });
            break;
        }
    }
    path
}

/// Computes the time derivative of the disc's state vector for use in RK4 integration.
///
/// The state vector is `[x, y, z, vx, vy, vz, spin]`:
/// - `x, y, z`: Position (meters)
/// - `vx, vy, vz`: Velocity components (m/s)
/// - `spin`: Spin rate (rad/s)
///
/// The returned derivative is `[dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt, dspin/dt]`.
///
/// Forces modeled:
/// - **Drag** (opposes velocity)
/// - **Lift** (proportional to velocity, spin, and other disc factors)
/// - **Gravity**
/// - **Wind** (applies a constant force in each direction)
///
/// # Parameters
/// - `state`: Current state vector.
/// - `disc`: Disc parameters.
/// - `env`: Environmental parameters.
///
/// # Returns
/// `PhysicsDerivative` for each state variable.
#[allow(clippy::similar_names)]
fn dynamics(state: &PhysicsState, disc: &Disc, env: &EnvParams) -> PhysicsDerivative {
    // All physical values: use units!
    let mass = disc.mass;
    let area = disc.area;
    let cd = disc.cd;
    let cl = disc.cl;
    let gyro = disc.gyro;
    let surface = disc.surface;
    let rho = env.rho;
    let g = env.g;

    // Total forces: sum with distinct names for Clippy!
    let rel_vx = state.vx.get::<meter_per_second>() - env.wind_x.get::<meter_per_second>();
    let rel_vy = state.vy.get::<meter_per_second>() - env.wind_y.get::<meter_per_second>();
    let rel_vz = state.vz.get::<meter_per_second>() - env.wind_z.get::<meter_per_second>();
    let rel_v = rel_vz
        .mul_add(rel_vz, rel_vx.mul_add(rel_vx, rel_vy * rel_vy))
        .sqrt();

    let dir = if rel_v > 1e-6 {
        [rel_vx / rel_v, rel_vy / rel_v, rel_vz / rel_v]
    } else {
        [0.0, 0.0, 0.0]
    };

    // SAFETY: -0.5 is unitless.
    // rho.value: [kg / m^3] (air density, scalar value)
    // cd.value: [unitless]   (drag coefficient, scalar value)
    // area.value: [m^2]      (cross-sectional area, scalar value)
    // rel_v: [m / s]         (relative velocity, scalar value)
    //
    // [kg/m^3] * [unitless] * [m^2] * [m/s] * [m/s]
    // = [kg/m^3] * [m^2] * [m^2/s^2]
    // = [kg] * [m^2] / [m^3] * [m^2] / [s^2]
    // = [kg] * [m^4] / [m^3 * s^2]
    // = [kg] * [m] / [s^2]
    // = [kg·m/s^2] => newton
    let drag_force: Force =
        Force::new::<newton>(-0.5 * rho.value * cd.value * area.value * rel_v * rel_v);
    let force_drag_x = drag_force * dir[0];
    let force_drag_y = drag_force * dir[1];
    let force_drag_z: Force = drag_force * dir[2];

    // Gravity force in Newtons
    let force_gravity_z = -mass * g; // this is a force (N)

    // Lift (if you model it as force)
    let spin_rps = state.spin.get::<radian_per_second>();

    // SAFETY: 0.5 is unitless.
    // cl.value: [unitless]   (lift coefficient, scalar value)
    // area.value: [m^2]      (cross-sectional area, scalar value)
    // rho.value: [kg/m^3]    (air density, scalar value)
    // rel_v: [m/s]           (relative velocity, scalar value)
    // spin_rps: [1/s]        (spin rate in rotations per second, scalar value)
    // gyro.value: [unitless] (gyroscopic stability, scalar value)
    // surface.value: [m^2]   (surface area, scalar value)
    //
    // [unitless] * [m^2] * [kg/m^3] * [m/s] * [1/s] * [unitless] * [m^2]
    // = [kg/m^3] * [m^2] * [m/s] * [1/s] * [m^2]
    // = [kg/m^3] * [m^2] * [m/s^2] * [m^2]
    // = [kg] * [m^4] / [m^3 * s^2]
    // = [kg] * [m] / [s^2]
    // = [kg·m/s^2] => newton
    let force_lift_z: Force = Force::new::<newton>(
        0.5 * cl.value * area.value * rho.value * rel_v * spin_rps * gyro.value * surface.value,
    );

    // Total force in each direction
    let total_force_x = force_drag_x;
    let total_force_y = force_drag_y;
    let total_force_z = force_drag_z + force_gravity_z + force_lift_z;

    // Accelerations (dV/dt): force / mass
    let acc_x = total_force_x / mass;
    let acc_y = total_force_y / mass;
    let acc_z = total_force_z / mass;

    PhysicsDerivative {
        dx: state.vx,
        dy: state.vy,
        dz: state.vz,
        dvx: acc_x,
        dvy: acc_y,
        dvz: acc_z,
        dspin: AngularAcceleration::new::<radian_per_second_squared>(0.0), // No spin acceleration yet
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::anyhow;
    use uom::si::length::meter;

    #[test]
    fn zero_velocity_disc_drops() -> anyhow::Result<()> {
        let v0 = Velocity::new::<meter_per_second>(0.0);
        let disc = Disc {
            v0,
            ..Default::default()
        };
        let env = EnvParams::default();
        let path = simulate(&disc, &env);
        assert!(!path.is_empty());
        assert!(path.last().ok_or_else(|| anyhow!("no path"))?.z <= Length::new::<meter>(0.0));
        Ok(())
    }

    #[test]
    fn high_spin_gives_lift() {
        let mut disc = Disc {
            spin: AngularVelocity::new::<radian_per_second>(100.0),
            ..Default::default()
        };
        let env = EnvParams::default();
        let path_spin = simulate(&disc, &env);

        disc.spin = AngularVelocity::new::<radian_per_second>(0.0);
        let path_no_spin = simulate(&disc, &env);

        // Expect high spin disc to stay in the air longer
        assert!(path_spin.len() > path_no_spin.len());
    }

    #[test]
    fn negative_mass_is_nonsense() {
        let disc = Disc {
            mass: Mass::new::<kilogram>(-1.0),
            ..Default::default()
        };
        let env = EnvParams::default();
        let path = simulate(&disc, &env);
        // Should not panic, but path will do weird things
        assert!(!path.is_empty());
    }

    #[test]
    fn high_angle_doesnt_crash() {
        let disc = Disc {
            angle: Angle::new::<radian>(80.0_f32.to_radians()),
            ..Default::default()
        };
        let env = EnvParams::default();
        let path = simulate(&disc, &env);
        assert!(!path.is_empty());
    }
}
