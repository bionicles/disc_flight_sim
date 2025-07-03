# disc_flight_sim

**A Rust + WASM disc golf flight physics simulator using 3D Runge-Kutta integration.**

## Features

- Real 3D physics (RK4)
- WASM-ready: run in browser with [wasm-pack](https://rustwasm.github.io/wasm-pack/)
- Fully documented equations
- Tests for edge cases, happy paths
- [WIP] Three.js integration for 3D rendering (see `/www`)

---

## Quick Start

### Prerequisites

- [Rust](https://www.rust-lang.org/tools/install): `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`
- [Just](https://github.com/casey/just): `cargo install just`
- [wasm-pack](https://rustwasm.github.io/wasm-pack/): `cargo install wasm-pack`

### Build & Test

```sh
git clone https://github.com/bionicles/disc_flight_sim.git
cargo test
```

# Build WASM and Serve Demo
```sh
just run
```
Open http://localhost:3000 and check console log.

# Usage (Rust)
```rust
use disc_flight_sim::{Disc, EnvParams, simulate_flight};

let disc = Disc::default();
let env = EnvParams::default();
let path = simulate_flight(&disc, &env);
```

# Usage (JS/Three.js)

you need to copy the `pkg` directory built by `wasm_pack` into the folder with your javascript
see `disc_flight_sim/Justfile/build` for details
```js
import init, { simulate_flight } from '../pkg/disc_flight_sim.js';
// call simulate_flight and use output points with Three.js
```

# TODO:

- [ ] render 3d trajectory from paths
- [ ] calculate s(t), v(t), a(t) on the rust side
- [ ] plot s(t), v(t), a(t) with `Chart.js` on the js side
- [ ] fuzz test the simulation and see if it crashes
- [x] add nice ui to input disc parameters
- [ ] omit negative drag coefficient
- [ ] simulate other things like flying saucers 🛸

# Contributing
PRs and bug reports are welcome. If your disc does a "loop-de-loop," tell us.
