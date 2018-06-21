# Optimal Shape Design for Poisson Equation

We study the shape design problem through the minimization of the cost functional

<p align="center">
  <img src="https://latex.codecogs.com/gif.latex?%5Cmathcal%7BJ%7D%20%5Cleft%28%20%5Ctheta%2C%20y%20%5Cright%29%20%3D%20%5Cfrac%7B1%7D%7B2%7D%20%5Cint_%7B%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%7D%20%5Cleft%28%20y%20-%20y_d%20%5Cright%29%20%5E2%20%5C%2C%20%5Cmathrm%7Bd%7D%20%5COmega%2C">
</p>

where <img src="https://latex.codecogs.com/gif.latex?y%20%5Cin%20L%5E2%20%5Cleft%28%20%5COmega%20%5Cright%29"> is the state variable, <img src="https://latex.codecogs.com/gif.latex?y_d%20%5Cin%20L%5E2%20%5Cleft%28%20%5COmega%20%5Cright%29"> is a target function, and <img src="https://latex.codecogs.com/gif.latex?%5Ctheta%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29"> the normal displacement to a reference boundary,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5CGamma%20%5Cleft%28%20%5Ctheta%20%5Cright%29%20%3D%20%5Cleft%5C%7B%20%5Cmathbf%7Bx%7D%20&plus;%20%5Ctheta%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29%20%5Cmathbf%7Bn%7D%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29%20%3A%20%5Cmathbf%7Bx%7D%20%5Cin%20%5CGamma_%7B0%7D%20%5Cright%5C%7D.">
</p>

The problem is subject to the following elliptic PDE in the domain <img src="https://latex.codecogs.com/gif.latex?%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%20%5Csubset%20%5Cmathbb%7BR%7D%5Ed"> with Dirichlet boundary conditions on <img src="https://latex.codecogs.com/gif.latex?%5CGamma%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%20%3D%20%5CGamma_w%20%5Ccup%20%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29">,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20-%5CDelta%20y%20%3D%20f%20%26%20%5Ctext%7Bin%20%7D%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%2C%5C%5C%20y%20%3D%20y_%7Bw%7D%20%26%20%5Ctext%7Bin%20%7D%20%5CGamma_w%2C%5C%5C%20y%20%3D%20y_%7Bs%7D%20%26%20%5Ctext%7Bin%20%7D%20%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%2C%20%5Cend%7Bcases%7D">
</p>

and <img src="https://latex.codecogs.com/gif.latex?f%20%5Cin%20L%5E2%20%5Cleft%28%20%5COmega%20%5Cright%29">.

We consider the set of admissible domains whose measure is fixed,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cmathcal%7BU%7D_%7Bad%7D%20%3D%20%5Cleft%5C%7B%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%20%3A%20%5Clvert%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%20%5Crvert%20%3D%20%5COmega_0%20%5Cright%5C%7D%2C">
</p>

and aim at finding the optimal shape that minimizes the cost function. In order to do so, we use the steepest descent method with descent direction given by

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cdelta%20%5Ctheta%20%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20%3D%20-%20%5Cleft%5B%20%5Cfrac%7B1%7D%7B2%7D%20%5Cleft%28%20y%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20-%20y_d%20%5Cright%29%20%5E2%20&plus;%20%5Cfrac%7B%5Cpartial%20v%7D%7B%5Cpartial%20n%7D%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20%5Cfrac%7B%5Cpartial%20y%7D%7B%5Cpartial%20n%7D%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20&plus;%20q%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20%5Cright%5D%2C">
</p>

where <img src="https://latex.codecogs.com/gif.latex?v"> is solution of the adjoint problem

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20-%20%5CDelta%20v%20%3D%20y%20-%20y_d%20%26%20%5Ctext%7Bin%20%7D%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%2C%20%5C%5C%20v%20%3D%200%20%26%20%5Ctext%7Bon%20%7D%20%5CGamma%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%20%3D%20%5CGamma_w%20%5Ccup%20%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29.%20%5Cend%7Bcases%7D">
</p>

The normal displacement field is updated at every iteration,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Ctheta%5E%7B%5Cleft%28%20n%20&plus;%201%20%5Cright%29%7D%20%3D%20%5Ctheta%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20&plus;%20%5Cepsilon%20%5Cdelta%7B%5Ctheta%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%7D%2C">
</p>

with <img src="https://latex.codecogs.com/gif.latex?\epsilon"> sufficiently small. The Lagrange multiplier is computed in order to ensure that the volume contraint is fulfilled, thus

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?q%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20%3D%20-%20%5Cfrac%7B1%7D%7B%5Cleft%7C%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%5Cright%7C%7D%20%5Cint_%7B%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%7D%20%5Cleft%5B%20%5Cfrac%7B1%7D%7B2%7D%20%5Cleft%28%20y%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20-%20y_d%20%5Cright%29%20%5E2%20&plus;%20%5Cfrac%7B%5Cpartial%20v%7D%7B%5Cpartial%20n%7D%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20%5Cfrac%7B%5Cpartial%20y%7D%7B%5Cpartial%20n%7D%5E%7B%5Cleft%28%20n%20%5Cright%29%7D%20%5Cright%5D%20%5Cmathrm%7Bd%7D%20%5CGamma.">
</p>

## Mesh Motion Solver

The solutions to the primal and adjoint problems are commonly approximated by means of numerical methods, such as the finite element method (FEM) or the finite volume method (FVM). In order to apply these techniques, the domain under study must be tessellated with a mesh. Nevertheless, applying the displacement directly to the controlled boundary will deteriorate the surrounding elements after a few iterations and the computation will crash if the interior nodes of the domain are not reallocated. In order to avoid this, the domain can be re-meshed after a number of iterations. However, this can be very expensive as a completely new mesh must be generated. A commonly used alternative is to move the interior nodes of the mesh according to the displacements prescribed on the boundary.  By doing this the number of elements and the nodes connectivities remain the same, only the mesh nodes positions are updated. The Solid Body Rotation Stress method and the Laplacian smoothing included in the OpenFOAM library have been used.

### Linear Elasticity

Mesh motion can be achieved by treating the mesh as an elastic body and solving the equations of Linear Elasticity for solids with prescribed displacements on the domain boundary,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20%5Cpartial_j%20%5Cleft%28%20G%20%5Cpartial_j%20u_i%20%5Cright%29%20&plus;%20%5Cpartial_j%20%5Cleft%28%20G%20%5Cpartial_i%20u_j%20%5Cright%29%20&plus;%20%5Cpartial_i%20%5Cleft%28%20%5Clambda%20%5Cpartial_j%20u_j%20%5Cright%29%20%3D%200%20%26%20%5Ctext%7Bin%20%7D%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%2C%20%5C%5C%20u_i%20%3D%200%20%26%20%5Ctext%7Bon%20%7D%20%5CGamma_w%2C%20%5C%5C%20u_i%20%3D%20%5Cdelta%20%5Ctheta%20n_i%20%26%20%5Ctext%7Bon%20%7D%20%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29.%20%5Cend%7Bcases%7D">
</p>

### Solid Body Rotation (SBR) Stress

The Linear Elasticity model fails when dealing with rotating meshes. This can be mitigated by selecting the material properties in a proper manner,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%202%20%5Cpartial_j%20%5Cleft%28%20%5Cgamma%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29%20%5Cpartial_j%20u_i%20%5Cright%29%20&plus;%20%5Cpartial_j%20%5Cleft%5B%20%5Cgamma%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29%20%5Cleft%28%20%5Cpartial_i%20u_j%20-%20%5Cpartial_j%20u_i%20-%20%5Cdelta_%7Bij%7D%20%5Cpartial_k%20u_k%20%5Cright%29%20%5Cright%5D%20%3D%200%2C%20%26%20%5Ctext%7Bin%20%7D%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%2C%20%5C%5C%20u_i%20%3D%200%20%26%20%5Ctext%7Bon%20%7D%20%5CGamma_w%2C%20%5C%5C%20u_i%20%3D%20%5Cdelta%20%5Ctheta%20n_i%20%26%20%5Ctext%7Bon%20%7D%20%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29.%20%5Cend%7Bcases%7D">
</p>

### Laplacian Smoothing

A simple and widely used practice is to solve a Laplace equation with the prescribed displacements as boundary conditions,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20%5CDelta%20u_i%20%3D%200%20%26%20%5Ctext%7Bin%20%7D%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29%2C%20%5C%5C%20u_i%20%3D%200%20%26%20%5Ctext%7Bon%20%7D%20%5CGamma_w%2C%20%5C%5C%20u_i%20%3D%20%5Cdelta%20%5Ctheta%20n_i%20%26%20%5Ctext%7Bon%20%7D%20%5CGamma_s%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29.%20%5Cend%7Bcases%7D">
</p>

The behavior of the Laplacian smoothing can be improved by adding a non-uniform diffusivity term,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cpartial_k%20%5Cleft%28%20%5Cgamma%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29%20%5Cpartial_k%20u_i%20%5Cright%29%20%3D%200%20%5Cquad%20%5Ctext%7Bin%20%7D%20%5COmega%20%5Cleft%28%20%5Ctheta%20%5Cright%20%29.">
</p>

The diffusion field <img src="https://latex.codecogs.com/gif.latex?%5Cgamma%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29%3A%20%5Cmathbb%7BR%7D%5Ed%20%5Crightarrow%20%5Cmathbb%7BR%7D"> decreases with the distance to the controlled boundary. A common choice is to make it depend on the inverse of the distance as

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cgamma%20%5Cleft%28%20%5Cmathbf%7Bx%7D%20%5Cright%29%20%3D%20%5Cfrac%7B1%7D%7Bd%5Em%7D%2C">
</p>

so that the nodes next to the deforming boundaries move with similar displacements as those on the boundary.

## Getting Started

The solver must be compiled in the terminal. It is advisable to first clean previous compilations with

```
wclean
```

and then use

```
wmake
```

### Dynamic Mesh

The mesh motion solver is specified in the dictionary _constant/dynamicMeshDict_.

* For the Laplacian solver:

```C++
dynamicFvMesh   dynamicMotionSolverFvMesh;

motionSolverLibs ( "libfvMotionSolvers.so" );

solver          displacementLaplacian;

displacementLaplacianCoeffs
{
    //diffusivity  	uniform;
    //diffusivity     	inversePointDistance (deformedWall);
    diffusivity     	quadratic inversePointDistance (deformedWall);
}
```

* For SBR Stress method:

```C++
motionSolver 	displacementSBRStress;

displacementSBRStressCoeffs
{
    // diffusivity  	uniform;
    // diffusivity  	directional (1 200 0);
    // diffusivity  	motionDirectional (1 1000 0);
    // diffusivity  	file motionDiffusivity;
    // diffusivity  	quadratic inverseFaceDistance (deformedWall);
    // diffusivity  	quadratic inverseDistance (deformedWall);
    diffusivity  	quadratic inversePointDistance (deformedWall);
    // diffusivity	inversePointDistance (deformedWall);
}
```

### Prerequisites

OpenFOAM C++ library must be installed in order to compile the code.

The OpenFOAM distribution provided by the [OpenFOAM Foundation](https://openfoam.org/) was used.

## Running a Case

In order to run the solver move to the case folder _poissonOptShapeFoamCase_ and type in the command line

```
./Allprepare

poissonOptShapeFoam
```

The shape optimization method described in the previous section has been tested with a simple example. The Poisson equation is posed in a two-dimensional circular domain with boundary

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5CGamma_w%20%3D%20%5Cleft%5C%7B%20%5Cleft%28%20x%2C%20y%20%5Cright%29%20%5Cin%20%5Cmathbb%7BR%7D%5E2%3A%20%5Csqrt%7Bx%5E2%20&plus;%20y%5E2%7D%20%3D%20R_w%20%5Cright%5C%7D.">
</p>

The reference geometry to be optimized is an inner hole with boundary given by

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5CGamma_s%20%5Cleft%28%200%5Cright%29%20%3D%20%5Cleft%5C%7B%20%5Cleft%28%20x%2C%20y%20%5Cright%29%20%5Cin%20%5Cmathbb%7BR%7D%5E2%20%3A%20%5Csqrt%7B%5Cleft%28%20x%20-%20c_x%20%5Cright%29%5E2%20&plus;%20%5Cleft%28%20y%20-%20c_y%20%5Cright%29%5E2%7D%20%3D%20R_s%20%5Cright%5C%7D.">
</p>

<p align="center">
  <img src="poissonOptShapeFoamCase/figs/fig1.png" width="400" height="340">
</p>

The target function <img src="https://latex.codecogs.com/gif.latex?u_d"> will be the analytical solution of the Poisson equation with the inner hole centered in the origin,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?u_%7Banalytic%7D%5Cleft%28%20r%20%5Cright%29%20%3D%20-%5Cfrac%7Bf%7D%7B4%7Dr%5E2%20&plus;%20c_1%20%5Clog%20r%20&plus;%20c_2%2C">
</p>

where the integration constants are obtained from the boundary conditions,

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Balign*%7D%20c_1%20%3D%20%26%20%5Cfrac%7Bu_w%20-%20u_s%20&plus;%20%5Cdisplaystyle%20%5Cfrac%7Bf%7D%7B4%7D%20%5Cleft%28%20R_w%5E2%20-%20R_s%5E2%20%5Cright%29%7D%7B%5Clog%5Cleft%28%20%5Cdisplaystyle%20%5Cfrac%7BR_w%7D%7BR_s%7D%20%5Cright%29%7D%2C%20%5C%5C%20c_2%20%3D%20%26%20%5Cfrac%7Bu_w%20&plus;%20u_s%7D%7B2%7D%20&plus;%20%5Cfrac%7Bf%7D%7B8%7D%5Cleft%28%20R_w%5E2%20&plus;%20R_s%5E2%20%5Cright%29%20-%20%5Cfrac%7Bc_1%7D%7B2%7D%20%5Clog%5Cleft%28%20R_w%20R_s%20%5Cright%29.%20%5Cend%7Balign*%7D">
</p>

When the hole center is displaced from the origin, the target state must be extended in the region <img src="https://latex.codecogs.com/gif.latex?r%20%5Cleq%20R_1">, thus

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?u_d%5Cleft%28%20r%20%5Cright%29%3D%20%5Cbegin%7Bcases%7D%20u_s%2C%20%5Cquad%20%26%20r%20%5Cleq%20R_s%2C%20%5C%5C%20u_%7Banalytic%7D%2C%20%5Cquad%20%26%20R_s%20%3C%20r%20%5Cleq%20R_w.%20%5Cend%7Bcases%7D">
</p>

It is clear that for

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5CGamma%5E%7B*%7D_s%20%3D%20%5Cleft%5C%7B%20%5Cleft%28%20x%2C%20y%20%5Cright%29%20%5Cin%20%5Cmathbb%7BR%7D%5E2%20%3A%20%5Csqrt%7Bx%5E2%20&plus;%20y%5E2%7D%20%3D%20R_s%20%5Cright%5C%7D">
</p>

the cost function equals zero, thus it is an optimal solution. The steepest descent algorithm has been coded in the OpenFOAM solver _poissonOptShapeFoam_ with <img src="https://latex.codecogs.com/gif.latex?%5Cepsilon%20%3D%2010%5E%7B-3%7D">.

<p align="center">
  <img src="poissonOptShapeFoamCase/figs/fig2.png" width="400" height="300">
  <img src="poissonOptShapeFoamCase/figs/fig3.png" width="400" height="300">
</p>

<p align="center">
  <img src="poissonOptShapeFoamCase/figs/laplace_uniform.gif" width="600" height="350">
</p>

### Warning

It might be needed to use 

```
sed -i -e 's/\r$//' filename
```

and

```
chmod +x filename
```

in order to be able to execute 

```
./filename
```

## References

* O. Pironneau. _Optimal shape design for elliptic systems_. Springer Science & Business Media, 2012.
* J Simon. _Diferenciación de problemas de contorno respecto del dominio_. Technical report, Universidad de Sevilla, Facultad de Matemáticas, Departamento de Análisis Matemático, 1989.
* [The OpenFOAM Foundation](http://openfoam.org).
* [CFD Direct](http://cfd.direct/openfoam/user-guide/).
* Moving boundary problem based on calculated data, [CFDonline](http://www.cfd-online.com/Forums/openfoam-programming-development/122557-moving-boundary-problem-based-calculated-data.html), 2013.

Linear Elasticity mesh motion method:
* Andrew A. Johnson and Tayfun E. Tezduyar. Mesh update strategies in parallel finite element computations of flow problems with moving boundaries and interfaces. _Computer methods in applied mechanics and engineering_, 119(1-2):73–94, 1994.
* Thomas D. Economon, Francisco Palacios, and Juan J. Alonso. Unsteady continuous adjoint approach for aerodynamic design on dynamic meshes. _AIAA Journal_, 53(9):2437–2453, 2015.
* George S. Eleftheriou and Guillaume Pierrot. Rigid motion mesh morpher: A novel approach for mesh deformation. In _OPT-i, An International Conference on Engineering and Applied Sciences Optimization_, Kos, Greece, pages 4–6, 2014.

Solid Body Rotation (SBR) Stress method:
* Richard P. Dwight. Robust mesh deformation using the linear elasticity equations. In _Computational fluid dynamics 2006_, pages 401–406. Springer, 2009.

Laplacian Smoothing:
* Peter Hansbo. Generalized laplacian smoothing of unstructured grids. _International Journal for Numerical Methods in Biomedical Engineering_, 11(5):455–464, 1995.
* Rainald Löhner and Chi Yang. Improved ale mesh velocities for moving bodies. _Communications in numerical methods in engineering_, 12(10):599–608, 1996.
* Hrvoje Jasak and Zeljko Tukovic. Automatic mesh motion for the unstructured finite volume method. _Transactions of FAMENA_, 30(2):1–20, 2006.
