# StagBL Core Components Proposal

## Design

A minimalist object-oriented C design, inspired by PETSc. We define classes
which may have several implementations defined by a "type".

Data structures, numberings, and other conventions should correspond exactly
to those used in PETSc, whenever possible.

When appropriate, classes here will have a "PETSc" implementation, which just
wraps the associated PETSc class. This possibility is a good rule of thumb
as to whether something is "core" or not.

## Core Classes

### StagBLGrid
- Global elementwise sizes
- Local element sizes

Implementations:
 - Default
 - PETSc (DMStag, but probably not using the section directly)

### StagBLArray
A multi-dimensional array, endowed with global
and local numberings.

Implementations:
 - Default (C arrays)
 - Implicit (Compute using a user-defined ctx and function)
 - Uniform (special case of implicit, optimized for speed)
 - Interpolated (special case of implicit, common use case) (?)
 - PETSc (Vec)
 - UTOPIA

### StagBLField

Little more than pointers to a StagBLGrid and a StagBLArray.

### StagBLScatter

An object which maps entries in one StagBLArray to another.

Implementations:
 - Default (MPI-1)
 - PETSc (VecScatter)
 - PETSc (PetscSF) (?)

### StagBLBoundaryIterator

### StagBLOperator

- assembled and matrix-free implementations
- Stencil information(?)
- Needs to support re-using the same array (same non-zero pattern) and changing the values.

### StagBLLinearSolver

Implementations:
 - Direct interface(s) to some direct solver package(s) (MUMPs or SuperLU or maybe even Elemental ..)
 - PETSc (KSP)
 - Native Krylov using StagBL components and ILUPACK (not a huge priority..)
 - Native MG using StagBL components only (again, aim to make data between StagBLGrid and DMStag data-compatible, so we can provide our StagBL components to PCMG if desired)
 - UTOPIA

### StagBLNonlinearSolver (?)
