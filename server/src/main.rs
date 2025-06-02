//! Minimal Axum server for `disc_flight_sim` (Rust + WASM + REST)
//!
//! # Quickstart
//! 1. cargo run
//! 2. curl `http://localhost:3000/`
//! 3. curl -X POST `http://localhost:3000/simulate_flight` -d '{"disc":{...},"env":{...}}' -H "Content-Type: application/json"
#![allow(clippy::multiple_crate_versions)]

use anyhow::anyhow;
use axum::{
    Router,
    extract::Json,
    response::{Html, IntoResponse},
    routing::{get, post},
};
use disc_flight_sim::{Disc, EnvParams, TrajectoryPoint, simulate};
use serde::{Deserialize, Serialize};
use tokio::{net::TcpListener, signal};
use tower_http::services::ServeDir;

#[derive(Deserialize)]
struct SimRequest {
    disc: Disc,
    env: EnvParams,
}

#[derive(Serialize)]
struct SimResponse {
    path: Vec<TrajectoryPoint>,
}

async fn root() -> Html<&'static str> {
    Html(include_str!("../www/index.html"))
}

async fn simulate_flight(Json(req): Json<SimRequest>) -> impl IntoResponse {
    let path = simulate(&req.disc, &req.env);
    axum::Json(SimResponse { path })
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    // build the app
    let static_dir = std::env::current_dir()?.join("www");
    let app = Router::new()
        .route("/", get(root))
        .route("/simulate_flight", post(simulate_flight))
        .fallback_service(ServeDir::new(static_dir));
    println!("app {app:#?}");

    // build the listener
    let listener = TcpListener::bind("127.0.0.1:3000").await?;
    println!("listening on {}", listener.local_addr()?);

    // handle interrupt signals
    let ctrl_c = async {
        match signal::ctrl_c().await {
            Ok(()) => Ok(()),
            Err(e) => Err(anyhow!("Failed to install Ctrl+C handler ({e})")),
        }
    };

    // run the server
    let server = axum::serve(listener, app);

    tokio::select! {
        _ = server => {
            println!("Server stopped");
        }
        _ = ctrl_c => {
            println!("Ctrl+C received, stopping server");
        }
    }
    Ok(())
}
