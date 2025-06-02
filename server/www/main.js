import init, { simulate_flight } from './pkg/disc_flight_sim.js';

async function run() {
    console.log("hello from main.js::run")
    await init();
    console.log("we did init")
    const disc = {
        v0: 20,
        angle: Math.PI/12,
        spin: 40,
        phi: 0,
        mass: 0.175,
        area: 0.056,
        cd: 0.6,
        cl: 0.7,
        gyro: 1.0,
        surface: 0.057
    };
    const env = {
        rho: 1.225,
        g: 9.81,
        wind_x: 0.0,
        wind_y: 0.0,
        wind_z: 0.0,
        dt: 0.02,
        max_t: 10.0
    };
    console.log("simulate", {disc, env})
    const path = simulate_flight(disc, env);

    if (path && path.length > 0) {
        console.log("First trajectory point from Rust:", JSON.stringify(path[0]));
        // console.log("Sample of path data from Rust:", JSON.stringify(path.slice(0, 3)));
    } else {
        console.log("Path data from Rust is empty or undefined.");
    }
    // The console.log("worked!", { path }); can be kept if desired, or removed.
    // For now, let's keep it as it shows the raw object structure if path is not empty.
    console.log("worked!", { path });

    // Data Processing
    const times = path.map(p => p.t);
    const displacements = path.map(p => Math.sqrt(p.x**2 + p.y**2 + p.z**2));
    const velocities = path.map(p => Math.sqrt(p.vx**2 + p.vy**2 + p.vz**2));
    const accelerations = path.map(p => Math.sqrt(p.ax**2 + p.ay**2 + p.az**2));

    console.log("Processed data for charts:", { times, displacements, velocities, accelerations });

    // Chart Creation
    const displacementCtx = document.getElementById('displacementChart').getContext('2d');
    const velocityCtx = document.getElementById('velocityChart').getContext('2d');
    const accelerationCtx = document.getElementById('accelerationChart').getContext('2d');

    new Chart(displacementCtx, {
        type: 'line',
        data: {
            labels: times,
            datasets: [{
                label: 'Position Magnitude',
                data: displacements,
                borderColor: 'rgb(75, 192, 192)',
                tension: 0.1
            }]
        },
        options: {
            scales: {
                y: {
                    beginAtZero: true // Magnitude starts at or above zero
                },
                x: {
                    title: {
                        display: true,
                        text: 'Time (s)'
                    }
                }
            }
        }
    });

    new Chart(velocityCtx, {
        type: 'line',
        data: {
            labels: times,
            datasets: [{
                label: 'Velocity Magnitude',
                data: velocities,
                borderColor: 'rgb(255, 99, 132)',
                tension: 0.1
            }]
        },
        options: {
            scales: {
                y: {
                    beginAtZero: true
                },
                x: {
                    title: {
                        display: true,
                        text: 'Time (s)'
                    }
                }
            }
        }
    });

    new Chart(accelerationCtx, {
        type: 'line',
        data: {
            labels: times,
            datasets: [{
                label: 'Acceleration Magnitude',
                data: accelerations,
                borderColor: 'rgb(54, 162, 235)',
                tension: 0.1
            }]
        },
        options: {
            scales: {
                y: {
                    beginAtZero: true
                },
                x: {
                    title: {
                        display: true,
                        text: 'Time (s)'
                    }
                }
            }
        }
    });
}

run();
