#![no_main]

use libfuzzer_sys::fuzz_target;
use disc_flight_sim::{simulate, Disc, EnvParams};
use serde::{Serialize, Deserialize};
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

// Helper macro to create UOM types from f32 values
macro_rules! uom_f32 {
    ($type:ident, $unit:ident, $value:expr) => {
        $type::new::<$unit>($value)
    };
}

fuzz_target!(|data: &[u8]| {
    if data.len() < std::mem::size_of::<DiscFuzzInput>() {
        return;
    }

    let mut disc_input = DiscFuzzInput::default();
    let mut env_input = EnvFuzzInput::default();

    // Use arbitrary_derive when it's available for a more robust way to create these inputs
    // For now, we'll manually create them from slices of data
    let mut offset = 0;
    let disc_input_size = std::mem::size_of::<DiscFuzzInput>();
    let env_input_size = std::mem::size_of::<EnvFuzzInput>();

    if data.len() >= disc_input_size {
        let (disc_bytes, rest) = data.split_at(disc_input_size);
        if let Ok(input) = postcard::from_bytes::<DiscFuzzInput>(disc_bytes) {
            disc_input = input;
        } else {
            // If deserialization fails, just use the default or skip
            return;
        }
        if rest.len() >= env_input_size {
            let (env_bytes, _) = rest.split_at(env_input_size);
             if let Ok(input) = postcard::from_bytes::<EnvFuzzInput>(env_bytes) {
                env_input = input;
            } else {
                // If deserialization fails, just use the default or skip
                return;
            }
        }
    } else {
        return; // Not enough data
    }

    let disc = Disc {
        v0: uom_f32!(Velocity, meter_per_second, disc_input.v0),
        angle: uom_f32!(Angle, radian, disc_input.angle),
        spin: uom_f32!(AngularVelocity, radian_per_second, disc_input.spin),
        phi: uom_f32!(Angle, radian, disc_input.phi),
        mass: uom_f32!(Mass, kilogram, disc_input.mass),
        area: uom_f32!(Area, square_meter, disc_input.area),
        cd: uom_f32!(Ratio, ratio, disc_input.cd),
        cl: uom_f32!(Ratio, ratio, disc_input.cl),
        gyro: uom_f32!(Ratio, ratio, disc_input.gyro),
        surface: uom_f32!(Area, square_meter, disc_input.surface),
    };

    let env = EnvParams {
        rho: uom_f32!(MassDensity, kilogram_per_cubic_meter, env_input.rho),
        g: uom_f32!(Acceleration, meter_per_second_squared, env_input.g),
        wind_x: uom_f32!(Velocity, meter_per_second, env_input.wind_x),
        wind_y: uom_f32!(Velocity, meter_per_second, env_input.wind_y),
        wind_z: uom_f32!(Velocity, meter_per_second, env_input.wind_z),
        dt: uom_f32!(Time, second, env_input.dt),
        max_t: uom_f32!(Time, second, env_input.max_t),
    };

    // Call the simulation function
    let _ = simulate(&disc, &env);
});

#[derive(Serialize, Deserialize, Debug, Default, Clone, Copy)]
struct DiscFuzzInput {
    v0: f32,
    angle: f32,
    spin: f32,
    phi: f32,
    mass: f32,
    area: f32,
    cd: f32,
    cl: f32,
    gyro: f32,
    surface: f32,
}

#[derive(Serialize, Deserialize, Debug, Default, Clone, Copy)]
struct EnvFuzzInput {
    rho: f32,
    g: f32,
    wind_x: f32,
    wind_y: f32,
    wind_z: f32,
    dt: f32,
    max_t: f32,
}
