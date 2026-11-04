# Constraint Lab

A small C++ 2D circle-physics engine and a no-install WebAssembly playground. Semi-implicit Euler at 120 Hz; Jacobian distance constraints, unilateral circle/wall contacts and Coulomb friction solved by sequential impulses.

Live: https://aditya2375.github.io/constraint-lab/

## Run the checks

```sh
g++ -std=c++20 -O2 -Wall -Wextra -Werror -DTEST_ENGINE engine.cpp -o engine-test
./engine-test
```

Tests check free fall, impulse transfer, anchor stability, joint drift, four 50-second stress runs, deterministic replay, storage/work bounds and dragging. Address/undefined sanitizers are also supported. These are numerical regression tests, not a proof of physical accuracy for every scene.

## Browser build

The committed demo is in `docs/`. No account, install, API key or network call is needed after its static assets load. Build from source with Emscripten 6.0.13 and `sh build.sh`. The demo has a 16 MB fixed WASM heap. The generated runtime is the compiler's loader; physics is in engine.cpp.

## Solver

For a center-to-center distance constraint C = |x_b - x_a| - L, the Jacobian is [-n, n]. We solve lambda = -(Jv + beta C/h) / (invMass_a + invMass_b), then apply M^-1 J^T lambda. Contacts accumulate and clamp normal impulse to >=0; friction clamps tangent impulse to +/- mu times accumulated normal impulse. Position drift correction adds energy, so this is not an energy-conserving integrator.

Reference: Erin Catto, [Modeling and Solving Constraints, GDC 2009](https://box2d.org/files/ErinCatto_ModelingAndSolvingConstraints_GDC2009.pdf). This is a separate educational implementation, not Box2D.

## Limits

Circles and center joints only: no angular dynamics, polygons, continuous collision detection, warm starting or sleeping. O(n²) collision broad phase and a 160-body limit. Discrete fast bodies can tunnel; speed is capped at 100 m/s for stability, which alters physics. The renderer caps catch-up and drops excess wall-clock time after slow frames. Hidden tabs stop simulation. Reduced motion starts paused. Drag uses a bounded velocity target, not a physical spring model. The HUD CPU time measures simulation work only, not total frame time.

Static frontend only, no backend/accounts/payments/telemetry. The repo's source is available for inspection; no third-party code or artwork copied into the solver.
