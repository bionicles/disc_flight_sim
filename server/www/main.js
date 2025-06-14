import init, { simulate_flight } from './pkg/disc_flight_sim.js';
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const MIN_DRAG_COEFFICIENT = 0.001;

function showToast(message) {
    const container = document.getElementById('toast-container');
    if (!container) {
        console.error('Toast container not found!');
        return;
    }

    const toast = document.createElement('div');
    toast.className = 'toast';
    toast.textContent = message;

    container.appendChild(toast);

    // Trigger reflow to enable animation
    requestAnimationFrame(() => {
        toast.classList.add('show');
    });

    setTimeout(() => {
        toast.classList.remove('show');
        // Remove the toast from DOM after animation
        setTimeout(() => {
            if (toast.parentNode === container) { // Check if still child before removing
                container.removeChild(toast);
            }
        }, 500); // Match CSS transition duration
    }, 3000); // Toast visible for 3 seconds
}

// Declare Three.js variables in a broader scope
let scene, camera, renderer, line = null, controls = null;

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

    // Add event listener for drag coefficient input clamping
    const cdInput = document.getElementById('cd');
    if (cdInput) {
        cdInput.addEventListener('input', () => {
            let val = parseFloat(cdInput.value);
            // Check if the value is a number and less than MIN_DRAG_COEFFICIENT
            if (!isNaN(val) && val < MIN_DRAG_COEFFICIENT) {
                cdInput.value = MIN_DRAG_COEFFICIENT.toString();
                showToast(`Drag coefficient was adjusted to the minimum allowed value of ${MIN_DRAG_COEFFICIENT}.`);
            } else if (isNaN(val) && cdInput.value.trim() !== '' && cdInput.value.trim() !== '-') {
                // Handle cases where input is not a valid number but not empty (e.g. "abc")
                // Optionally, clear it or set to default, or just let getDiscParametersFromUI handle it
                // For now, we'll let getDiscParametersFromUI handle NaN from non-numeric strings
            }
        });
    } else {
        console.error("Drag coefficient input field ('cd') not found!");
    }

    // Initialize Three.js scene, camera, and renderer
    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
    renderer = new THREE.WebGLRenderer();
    renderer.setSize(window.innerWidth, window.innerHeight);
    controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true; // an optional feature that adds inertia to the controls
    controls.dampingFactor = 0.25;
    controls.enableZoom = true;

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
        if (controls) {
            controls.update();
        }
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
