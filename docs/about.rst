About StagBL
============

StagBL is a C library providing a discretization (grid) and solver layer for a focused set of target applications:
regional and global geodynamics codes, on orthogonal, regular finite volume meshes of quadrilaterals or hexahedra.
In particular, we focus on the applications `StagYY <https://doi.org/10.1016/j.pepi.2008.08.005>`__, I3ELVIS \cite{GeryaYuen2009},
and `LaMEM <https://bitbucket.org/bkaus/lamem>`__.

The scalability bottleneck in these applications is the solution of large, coupled systems of conservation equations, in particular the Stokes equations with heterogeneous viscosity.
The development of the library also produced contributions to `PETSc <https://mcs.anl.gov/petsc>`__, notably the DMStag object for working with distributed, regular, staggered grids. This allows StagBL to remain focused on its applications, but for reusable components to be more widely used (hence maintained).
 Scalable solvers for these systems have been a subject of much research, but making these advanced solvers available to practitioners requires composable software tools. The use of a unified layer also allows operations to optimized for new computational architectures.
