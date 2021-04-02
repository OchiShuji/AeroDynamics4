### AeroDynamics4

#### Numerical Analysis of scalar 1D linear hyperbolic PDE
main_1.m : Solve scalar 1D linear hyperbolic PDE by the following three methods. <p> 1) 1st order up-wind differential scheme <p> 2) Lax-Wendroff Scheme <p>3)Yee's Symetric-TVD Scheme
<p>When adapting TVD scheme to the program, two types of limiter function are used; minmod, and superbee.

main_2.m : Solve scalar PDE system(2D) by Godunov's method. 1st-order up-wind scheem is used for the caluculation of numerical flux.

main_3.m : Solve shock tube problem by FDS(Flux Difference Spliiting). Roe's approximate method is used for Riemann Solver. 
Symmetric TVD is used for the caluculation of numerical flux. Minimod is used as a limiter function.