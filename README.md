# disc_flight_sim

**A Rust + WASM disc golf flight physics simulator using 3D Runge-Kutta integration.**

## Features

- Real 3D physics (RK4)
- WASM-ready: run in browser with [wasm-pack](https://rustwasm.github.io/wasm-pack/)
- Fully documented equations
- Tests for edge cases, happy paths
- Three.js integration for 3D rendering (see `/www`)

---

## Quick Start

### Prerequisites

- Rust (`rustup`)
- [wasm-pack](https://rustwasm.github.io/wasm-pack/): `cargo install wasm-pack`
- Node.js & npm (for the Three.js demo)

### Build & Test

```sh
git clone https://github.com/YOURNAME/disc_flight_sim.git
cd disc_flight_sim
cargo test
```

# Build WASM and Serve Demo
```sh
wasm-pack build --target web
```
Open http://localhost:8080 and toss some discs.

# Usage (Rust)
```rust
use disc_flight_sim::{Disc, EnvParams, simulate_flight};

let disc = Disc::default();
let env = EnvParams::default();
let path = simulate_flight(&disc, &env);
```

# Usage (JS/Three.js)

```js
import init, { simulate_flight } from '../pkg/disc_flight_sim.js';
// call simulate_flight and use output points with Three.js
```

# TODO:

- [ ] fuzz testing

# Contributing
PRs and bug reports are welcome. If your disc does a "loop-de-loop," tell us.