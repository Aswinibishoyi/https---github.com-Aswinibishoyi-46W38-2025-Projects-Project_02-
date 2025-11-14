The project takes time stamped wind speed from multiple files(which are parked in input file, however as per the project requirement it is in multiple files/three sub-folders).
Using input wind speed, interpolated thrust coefficient is calculated.
To derive Turbie’s mass, stiffness and damping matrices, we make a few assumptions about our dynamical system:
The turbine can only move in the fore-aft direction; The 3 blade deflections in the fore-aft direction are syncronized, i.e., only collective flapwise deflections.
With these assumptions, the Turbie is reduced to a 2DOF mass-spring-damper system.
With the above assumption a 2DOF system was create with the given matrix(M,C and K) in the input file.
aerodynamic forcing on the blades is formulated as per question, where where A is the rotor area and u(t) is the wind speed at time t. The thrust coefficient Ct is determined from the mean wind speed U=u(t) using the look-up table defined in C_T.txt.
