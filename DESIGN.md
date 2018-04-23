This document records design choices and guidelines for StagBL.
Each point should include a rationale [in brackets].

# Code style

[Patrick's (in progress) style guide](https://bitbucket.org/psanan/pdsstyle)

# Documentation

Use .rst format [Human-readable, can go on readthedocs.org]

# Conventions

## Grids

Grid sizes are specified in terms of elements (top-level cells, e.g quads in 2D and hexes in 3D). [The discretization is based on these control volumes. It is arguably easier to reason about. It might make more sense with ghost points].


