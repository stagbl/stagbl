====================
The Design of StagBL
====================

Here, some motivating design principles for StagBL are discussed, along
with some of their consequences.

Focus
-----

StagBL maintains a sharp focus on its target applications: finite volume
solvers for simulation of lithospheric and mantle dynamics, using orthogonal,
regular finite volume discretizations, often coupled to particle systems.

Maintainability
---------------

StagBL strives to have a small code base.

A key technique to achieve this aim is extensive leverage of tools from the
PETSc library.

* StagBL can take advantage of PETSc's robust configuration system, which takes
  the place of a tool like CMake or GNU Autotools, as well as assuming much of
  the functionality of a package manager, for example making it simple to
  automatically download and install direct solver packages.
* StagBL uses PETSc's error handling.
* StagBL uses data types corresponding to those from the PETSc library, for
  example ``PetscInt``, ``PetscErrorCode``, and ``PetscScalar``.

PETSc's permissive BSD-2 license exposes the useful possibility to contribute
functionality to the PETSc library. This has been done with the DMStag object
and related code, implemented as part of the StagBL PASC Project. This allows
some of our components to be used (and hence tested and improved) beyond the scope
of our geodynamics applications, and keeps the StagBL library's native code
more focused.

StagBL's documentation is by example as much as possible.
