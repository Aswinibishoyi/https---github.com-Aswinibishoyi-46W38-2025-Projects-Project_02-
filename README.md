The project takes time stamped wind speed from multiple files(which are parked in input file, however as per the project requirement it is in multiple files/three sub-folders).
Using input wind speed, interpolated thrust coefficient is calculated.
To derive Turbie’s mass, stiffness and damping matrices, we make a few assumptions about our dynamical system:
The turbine can only move in the fore-aft direction; The 3 blade deflections in the fore-aft direction are syncronized, i.e., only collective flapwise deflections.
With these assumptions, the Turbie is reduced to a 2DOF mass-spring-damper system.
With the above assumption a 2DOF system was create with the given matrix(M,C and K) in the input file.
aerodynamic forcing on the blades is formulated as per question, where where A is the rotor area and u(t) is the wind speed at time t. The thrust coefficient Ct is determined from the mean wind speed U=u(t) using the look-up table defined in C_T.txt.
The full dynamical equations for Turbie formulated with the state vector, and the forcing vector. By simulating Turbie's response to wind loads, the derivatives of the blade and tower deflections and their velocities at each time step determined.
With a second Python script, these functions are called to solve Turbie for each wind speed case.
y¯′(t) is passed to the scipy.integrate.solve_ivp function (along with initial conditions for y) which output the blade and tower displacements and velocities of Turbie at each time step.
For each wind speed cases, time-marching variation of wind and the blade and tower displacements are plotted and saved in the  results folder.
The means and standard deviations of the blade and tower displacements for each wind speed calculated and saved within a separate text file for each TI category in the results folder.
Discussion : 
The mean deflections primarily driven by the average wind speed, not the turbulence intensity. The mean does not change much across different TIs if the average wind speed is kept the same.
TheStandard deviation is a measure of variation. Higher Turbulence Intensity (TI) at the same average wind speed should lead to larger fluctuations in aerodynamic force and thus larger standard deviations in the blade and tower displacements.
