# A NITSCHE-BASED METHOD FOR SOLVING A PROBLEM IN CONTACT MECHANICS - Parallel version
## S. Galati
### A.A. 2024-2025

The project aims to solve a contact mechanics problem using a Nitsche-based method for imposing both the non-penetrability and the Coulomb friction condition at a fault cutting a 3-dimensional bulk domain.
The FEM framework is given by GetFEM++, a very powerful library for solving many FEM problems. In order to compile and execute the code, the user will need the following libraries installed:
- getfem
- mumps (parallel version): for the linear algebra procedures
- qhull: for the mesh generation of getfem

Before installing getfem:
- be sure you have installed blas, lapack (usually comes automatically with mumps), qhull, metis, openblas and the parallel version of mumps, scalapack
- configure it with the following options:
``` ! bash
./configure --enable-shared --enable-mumps --enable-parallel=2 --enable-qhull --enable-metis --enable-openmp
```

The mesh can be built either internally to getfem (using qhull) or by using gmsh and giving the option `-m` at runtime (in such a case, gmsh needs to be installed). The user has just to modify the data.pot file properly.

The `src` and `include` directories contain the code, `external` just contains a copy of muparserx (which is used to read functions form the .pot file), while the `build` directoy is used to store object and dependency files built when the code gets compiled.

### Source code
The main implemented classes are:
`Params`:
A struct holding all the parameters contained in `data.pot`.
- It implements a function for reading the parameters from a GetPot datafile, including some options passed at runtime.
- Overloads the streaming operator.

`Mesh` and `MeshBuilderStrategy`: the first is just a wrapper for the `getfem::mesh` object, which is built using the Builder in `MeshBuilder.hpp`. It allows direct access to the underlying regions definied within the mesh (Bulks, Fault and Boundary regions).

`FEMManager`: to manage the `getfem::mesh_fem` objects.

`BCHandler` and `BC`: to handle the boundary conditions. The supported BC types are Dirichlet, Neumann and normal Dirichlet. In the data.pot:
- region< type > = [1 4 5]: Dirichlet is imposed on boundary regions with tags 1,4 and 5
- bd< type >< i > = {x[1]+sqrt(x[0]),x[2],3.5*t^2}: the < type > to be imposed on the i-th region given in region< type >

`ContactProblem`: the main class handling the problem. It implements methods to:
- initialize the problem: set BCs, set FEM, set integration method
- prepare the assemble: use GFWL of getfem to instruct the solver on how to assemble the matrices at each iteration
- solve the problem: implements the time-advancing scheme and calls `getfem::standard_solve(...)`, which has already Netwon-like solver implemented and optimally chooses the best solver depending on the information given in the `assemble()` method.
- export the solution in `vtk` format.


## Notes for compiling the code
In the Makefile, if you have the aforementioned libraries installed in non-standard paths, you need to add them in `CPPFLAGS`, then simply run `make`.

## Optional runtime flags
Input data can be changed at runtime modifying `data.pot`. To run the executable type 
```bash
mpirun -n <nprocs> ./main
```
with possible options:
  - `-v`: be verbose
  - `-m`: generate the mesh with gmsh
  - `-f <filename>`: to give a different filename (`data.pot` by default).

Notice that if the mesh isn't generated using gmsh, the angle parameter will be ignored
