# Fluid Dynamics
### or: Fluids are fun. But man, are they hard to simulate

This is a simple repo for simple men. Always wanted to write a decent fluid simulation, this is currently my best shot at it.
Of course, I will be nowhere near this poit were it not for the invaluable "Fluid Simulation for Computer Graphics" by Rober Bridson.
But some other smart guys (unwittingly) helped me along my days of debugging: the GitHub repo 
[`incremental-fluids`](https://github.com/tunabrain/incremental-fluids) by *tunabrain* stands out (with the small downside of being a 
fluid-solid solver, no space for air, but hey, I got you covered for that).

What is this program so far: a PIC/FLIP 2D Fluid Simulation, absolutely unoptimized.

### Features:

- It does not crash (at least not immediatly);
- It seems to provide reasonable output;
- Preconditioned-Conjugate-Gradient Pressure Solver;
- PIC/FLIP for improved particle-based advection;
- Level-Set tracking of FLUID-AIR-SOLID geometry;
- `lodepng.h` for saving simulation frame to `.png` images;
- Improved Curved Boundaries (not tested);
- Improved Free Surfaces control;

### What still need to be implemented:

- Extensive optimization. Seriously, a man's patience has limits.
- ~~Better free surfaces (Bridson: Cap 8.2)~~ Done! ✓;
- Fluid sources/sinks;
- Better output graphics;
- Variable densities;
- Vorticity confinement;
- ...3D (what a dream!);

### Show me that it works!

Marvel at this magnificent melting earth (rendered with online GifMaker tool [egzif.com](https://ezgif.com/)).

![Earth is melting. Do something!](earth.gif)
