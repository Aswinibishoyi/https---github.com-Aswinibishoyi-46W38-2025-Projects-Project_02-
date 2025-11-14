The project takes time stamped wind speed from multiple files(which are parked in input file, however as per the project requirement it is in multiple files/three sub-folders).
Using input wind speed, interpolated thrust coefficient is calculated.
To derive Turbie’s mass, stiffness and damping matrices, we make a few assumptions about our dynamical system:
The turbine can only move in the fore-aft direction; The 3 blade deflections in the fore-aft direction are syncronized, i.e., only collective flapwise deflections.
With these assumptions, the Turbie is reduced to a 2DOF mass-spring-damper system.
