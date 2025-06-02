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
    console.log("worked!", { path })
    // Draw path with Three.js, basic line geometry
    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(75, window.innerWidth/window.innerHeight, 0.1, 1000);
    const renderer = new THREE.WebGLRenderer();
    renderer.setSize(window.innerWidth, window.innerHeight);
    document.body.appendChild(renderer.domElement);

    const points = path.map(p => new THREE.Vector3(p.x, p.y, p.z));
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    const material = new THREE.LineBasicMaterial({ color: 0xff0000 });
    const line = new THREE.Line(geometry, material);
    scene.add(line);

    camera.position.z = 20;
    camera.position.y = -10;
    camera.position.x = 10;
    camera.lookAt(new THREE.Vector3(0, 0, 0));

    function animate() {
        requestAnimationFrame(animate);
        renderer.render(scene, camera);
    }
    animate();
}

run();
