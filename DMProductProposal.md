DMProduct: Proposal
-------------------

# Motivation
For applications relying on a regular, orthogonal grid, coordinates
are often representable in a highly compressed format. For a grid in N dimensions,
a grid point (i,j,k) has coordinates (f(i),g(i),h(i)).

PETSc's DM object holds a second DM object which is assumed to be "compatible"
and which represents coordinate fields on the DM.

A motivating case is that of modelling the mantle of the Earth. Important structure is known
a priori in the z (radial) direction, motivating nonuniform coordinates to give,
for example, more resolution near the surface, the core-mantle boundary, or the
660 km discontinuity. In the other two dimensions, however, a regular grid is often
used. Thus, one can efficiently represent the coefficients as a product of functions,
some of which may be implmented as table lookups or lower-dimensional DMs.

Another important example arised in atmospheric modelling. A domain may be represented with an unstructured
mesh in the plan, times a uniform vertical mesh.

It is inefficient to use the same type of DM object to represent the coordinates
as is used to represent solution or coefficient fields on the domain. PETSc does
not provide a native object to represent these products, so applications with this
structure must implement these lookups themselves.

This proposal describes a new DM implementation, DMProduct, intended to represent
this sort of structured grid.

Key is the notion of "compatiblity", meaning that a section on one DM can be interpreted
as a section on another.

# Design

For dimension D, DMProduct holds a D parameters of type DMProductSlotType. These
determines how coordinates in this dimension are computed:
- by a given coordinate of a sub-DM
- by a built-in function (likely only uniform spacing supported)
- by a given coordinate of a user-specified function of the global index

(Partially-filled) arrays of size D hold pointers to DMs and custom functions.
