This document records design choices an guidelines for StagBL.
Each point should include a rationale [in brackets].

# Code style

General : [Patrick's (in progress) style guide](https://bitbucket.org/psanan/pdsstyle)

PETSc : [PETSc Developer's Manual](http://www.mcs.anl.gov/petsc/petsc-current/docs/developers.pdf)

# Conventions 

## Grids

Grid sizes are specified in terms of elements (top-level cells, e.g quads in 2D and hexes in 3D). [The discretization is based on these control volumes. It is arguably easier to reason about. It might make more sense with ghost points].


