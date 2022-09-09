# Unsteady Laminar Flow Solver

This code models a low speed 2D flow around a square cylinder using a finite volume approach with explicit Pressure-Correction and explicit Forward Euler time stepping.The analysis focuses on the effects of the Reynolds number, as determined by the input velocity parameters, on the flow around the square.

- Written in C++ using the Eigen LIbrary
- Use of the SIMPLE Algorithm on a staggered grid
- Valid only for orthogonal mesh
- Post Processing in Paraview

The code has been used to also solve simple orthogonal geometries, such as Lid Driven Cavity problems, Channel Flows and Flow over Square cylinder problems. Bibliographical data have been used to validate the results.

[1] Atul Sharma & V. Eswaran (2004): HEAT AND FLUID FLOW ACROSS A SQUARE CYLINDER IN THE TWO-DIMENSIONAL LAMINAR FLOW REGIME, Numerical Heat Transfer, Part A: Applications: An International Journal of Computation and Methodology, 45:3, 247-269.

[2] A. Sohankar, C. Norberg, and L. Davidson, Low-Reynolds Number Flow around a Square Cylinder at Incidence: Study of Blockage, Onset of Vortex Shedding and Outlet Boundary Condition, Int. J. Numer. Meth. Fluids, vol. 26, pp. 39–56, 1998.

[3] M. Provansal, C. Mathis, and L. Boyer, Be ́ nard-von Ka ́rma ́ n Instability: Transient and Forced Regimes, J. Fluid Mech., vol. 182, pp. 1–22, 1987.

[4] Calhoun, D., A Cartesian Grid Method for Solving the Two-Dimensional Streamfunction-Vorticity Equa- tions in Irregular Regions, Journal of Computational Physics, 176 (2002), 2, pp. 231-275
