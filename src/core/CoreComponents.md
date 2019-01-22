# StagBL Core Components Proposal

## Design

A minimalist object-oriented C design, inspired by PETSc. We define classes
which may have several implementations defined by a "type". This type is
determined based both on a required backend library and/or specifics of the
implementation.

Data structures, numberings, function signatures, and other conventions should
correspond exactly to those used in PETSc, whenever possible.

When appropriate, classes here will have a "PETSc" implementation, which
includes/wraps an associated PETSc class. This possibility is a good rule of
thumb as to whether something is "core" or not. We note that StagBL's scope is
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

### StagBLGridYY
Implementations :
 - Composite (2 StagBLGrid members)
 - PETSc (DMComposite object holding 2 DMStag members)

### StagBLArray
A multi-dimensional array, plus a reference to a StagBLGrid.

At any given time, it may store explicit representations of local and/or global
data.

Implementations:
 - Default (C arrays)
 - PETSc (Petsc Vecs)
 - UTOPIA

### StagBLOperator
- Stencil information(?)
- Should support re-using the same array (same non-zero pattern) w/ new values.

Implementations:
 - Native Assembled
 - Native Matrix-Free
 - PETSc (Mat)
 - UTOPIA

### StagBLLinearSolver

Implementations:
 - Direct interface(s) to some direct solver package(s)
 - PETSc (KSP)
 - Native MG using StagBL components only
 - UTOPIA
 - Native Krylov using StagBL components and ILUPACK
