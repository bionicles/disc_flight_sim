import init, { simulate_flight } from './pkg/disc_flight_sim.js';

// Declare Three.js variables in a broader scope
let scene, camera, renderer, line = null;

// Environment parameters (assuming these are fixed for now)
const env = {
    rho: 1.225,
    g: 9.81,
    wind_x: 0.0,
    wind_y: 0.0,
    wind_z: 0.0,
    dt: 0.02,
    max_t: 10.0 // Max simulation time
};

function getDiscParametersFromUI() {
    // Default values
    const defaults = {
        v0: 20,
        angleDegrees: 15, // Math.PI/12 radians
        spin: 40,
        phiDegrees: 0,
        mass: 0.175,
        area: 0.056,
        cd: 0.6,
        cl: 0.7,
        gyro: 1.0,
        surface: 0.057
    };

    let v0 = parseFloat(document.getElementById('v0').value);
    if (isNaN(v0)) v0 = defaults.v0;

    let angleDegrees = parseFloat(document.getElementById('angle').value);
    if (isNaN(angleDegrees)) angleDegrees = defaults.angleDegrees;

    let spin = parseFloat(document.getElementById('spin').value);
    if (isNaN(spin)) spin = defaults.spin;

    let phiDegrees = parseFloat(document.getElementById('phi').value);
    if (isNaN(phiDegrees)) phiDegrees = defaults.phiDegrees;

    let mass = parseFloat(document.getElementById('mass').value);
    if (isNaN(mass)) mass = defaults.mass;

    let area = parseFloat(document.getElementById('area').value);
    if (isNaN(area)) area = defaults.area;

    let cd = parseFloat(document.getElementById('cd').value);
    if (isNaN(cd)) cd = defaults.cd;

    let cl = parseFloat(document.getElementById('cl').value);
    if (isNaN(cl)) cl = defaults.cl;

    let gyro = parseFloat(document.getElementById('gyro').value);
    if (isNaN(gyro)) gyro = defaults.gyro;

    let surface = parseFloat(document.getElementById('surface').value);
    if (isNaN(surface)) surface = defaults.surface;

    // Convert angles from degrees to radians for the simulation
    const angleRadians = angleDegrees * Math.PI / 180;
    const phiRadians = phiDegrees * Math.PI / 180;

    return {
        v0,
        angle: angleRadians,
        spin,
        phi: phiRadians,
        mass,
        area,
        cd,
        cl,
        gyro,
        surface
    };
}

function updateSimulation(discParams, environmentParams) {
    console.log("Updating simulation with:", { discParams, environmentParams });

    // Clear previous line from scene if it exists
    if (line) {
        scene.remove(line);
        line.geometry.dispose(); // Dispose of the old geometry to free up resources
        line.material.dispose(); // Dispose of the old material
        line = null;
    }

    const path = simulate_flight(discParams, environmentParams);
    console.log("Simulation worked!", { path });

    if (path && path.length > 0) {
        const points = path.map(p => new THREE.Vector3(p.x, p.y, p.z));
        const geometry = new THREE.BufferGeometry().setFromPoints(points);
        const material = new THREE.LineBasicMaterial({ color: 0xff00ff }); // Changed color for visibility of update
        line = new THREE.Line(geometry, material);
        scene.add(line);

        // Optional: Adjust camera to fit the new path, or keep it static
        // For now, keeping camera static as per original setup
    } else {
        console.warn("Simulation returned no path or an empty path.");
    }
}

async function run() {
    console.log("hello from main.js::run");
    await init();
    console.log("WASM module initialized");

    // Initialize Three.js scene, camera, and renderer
    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
    renderer = new THREE.WebGLRenderer();
    renderer.setSize(window.innerWidth, window.innerHeight);

    // Check if renderer.domElement is already a child of body
    // This can happen if run() is called multiple times, though current structure doesn't do that.
    // For safety, especially during development, ensure it's not added multiple times.
    if (!document.body.contains(renderer.domElement)) {
        const container = document.getElementById('container');
        if (container) {
            container.appendChild(renderer.domElement);
        } else {
            console.error("Container element not found for renderer!");
            document.body.appendChild(renderer.domElement); // Fallback to body
        }
    }

    // Initial camera setup
    camera.position.set(10, -10, 20); // Initial camera position
    camera.lookAt(new THREE.Vector3(0, 0, 0)); // Look at origin

    // Animation loop
    function animate() {
        requestAnimationFrame(animate);
        renderer.render(scene, camera);
    }
    animate();

    // Get initial parameters from UI and run simulation
    const initialDiscParams = getDiscParametersFromUI();
    updateSimulation(initialDiscParams, env);

    // Add event listener to the "Simulate" button
    const simulateButton = document.getElementById('simulate-button');
    if (simulateButton) {
        simulateButton.addEventListener('click', () => {
            console.log("Simulate button clicked");
            const newDiscParams = getDiscParametersFromUI();
            updateSimulation(newDiscParams, env);
        });
    } else {
        console.error("Simulate button not found!");
    }
}

run().catch(console.error);
