# StagBL Core Components Proposal

## Design

A minimalist object-oriented C design, inspired by PETSc. We define classes
which may have several implementations defined by a "type". This type is
determined based both on a required backend library and/or specifics of the
implementation.

We avoid a deep class hierarchy, but code duplication is kept to a minimum 
by reusing implementations of API functions when possible.

Whenever possible, objects are used to create compatible objects. For instance,
a grid can create a compatible array.

Data structures, numberings, function signatures, and other conventions should
correspond exactly to those used in PETSc, whenever possible.

When appropriate, classes here will have a "PETSc" implementation, which
includes/wraps an associated PETSc class or classes. This possibility is a good 
rule of thumb as to whether something is "core" or not. We note that StagBL's scope is
smaller than that of the more generalized objects being wrapped (e.g. we support
Stokes and temperature operators useful for geodynamics, not the much larger class
of linear operators that PETSc's Mat object represents).  This comes at the cost
that existing users of PETSc are left without tools and flexibility that they've
become accustomed to. To alleviate this, we always attempt to allow an "escape
hatch" to access internal PETSc objects and to interact with them via command-line
options.

## Core Classes

### StagBLGrid
Implementations:
 - Native (data layouts identical to DMStag)
 - PETSc (DMStag)
 - Multi (several compatible grids)
 - YY (two grids)
 - Composite (general combination of grids)

Note: the composite implementations use DMComposite when all members also use PETSc

### StagBLArray
A distributed array, plus an immutable reference to a StagBLGrid.

At any given time, it may store explicit representations of local and/or global
data, and maintains state as to which are "up to date".

It features routines to obtain read-write or write-only access to the array data.

Implementations:
 - Default (C arrays)
 - PETSc (Petsc Vecs, holding refernce to DMStag)
 - UTOPIA
 - Multi
 - YY
 - Composite

### StagBLSystem
This knows how to compute a residual on a grid, from values on zero or more compatible grids.

Contains an immutable reference to a StagBLGrid.

Implementations:
 - Native Assembled
 - Native Matrix-Free
 - PETSc (Mat + Vec)
 - UTOPIA

### StagBLSolver
Holds an immutable reference to a StagBLSystem.

This knows how to take an initial guess (could be implicit zero) of an array,
and compute a lower-residual update with respect to a StagBLSystem.

Implementations:
 - Direct interface(s) to some direct solver package(s)
 - PETSc (KSP)
 - Native MG using StagBL components only
 - UTOPIA
 - Native Krylov using StagBL components and ILUPACK
