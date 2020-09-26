====================
The Design of StagBL
====================

Here, some motivating design principles for StagBL are discussed, along
with some of their consequences.

Focus
=====

StagBL maintains a sharp focus on its target applications: finite volume
solvers for simulation of lithospheric and mantle dynamics, using orthogonal,
regular finite volume discretizations, often coupled to particle systems.

Maintainability
===============

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

Core Objects
============

StagBL presents a small set of main classes.

These are implemented with a minimalist object-oriented C design, inspired by PETSc.
We define classes which may have several implementations defined by a "type". This type
may be determined based both on a required backend library and/or specifics of the
implementation.

Whenever possible, objects are used to create compatible objects. For instance,
a grid can create a compatible array.

Classes here have default "PETSc" implementations, which includes/wraps an
associated PETSc class or classes. This possibility is a good rule of thumb as
to whether something is "core" or not. We note that StagBL's scope is smaller
than that of the more generalized objects being wrapped (e.g. we support Stokes
and temperature operators useful for geodynamics, not the much larger class of
linear operators that PETSc's Mat object represents).  This might cause
existing users of PETSc to be left without tools and flexibility that they've
become accustomed to. To alleviate this, we always expose an "escape hatch" to
access internal PETSc objects and to interact with them via command-line
options.

StagBLGrid
----------

An object representing a parallel domain comprised of one or more
regular staggered grids.

The default implementation uses a ``DM`` object from PETSc, typically
a DMStag object or composite thereof.

StagBLArray
-----------
A distributed array, plus an immutable reference to a StagBLGrid.

The default implmentation uses a Petsc ``Vec`` object.

StagBLSystem
------------
Information on how to compute a residual on a grid,
and how to "solve" to compute a zero/minimum of this residual.
The default implementation uses a Petsc ``SNES`` object.

Contains an immutable reference to a StagBLGrid.
