DMStag - Component Proposal
---------------------------

`DMStag` is an implementation of PETSc's `DM` abstract class.

It represents a topologically Cartesian cell complex based on a grid of line segments, quadrilaterals, or hexahedra.
It provides routines to manipulate interacting fields on these cells and all lower-dimensional cells in the complex.

The intended use is for finite volume or DEC schemes, on regular grids.

It is intended to be intermediate between `DMDA` an `DMPlex`; it represents a structured grid like `DMDA`, but supports different strata like `DMPlex`.

## Terminology

- Point   : a cell of any dimension
- Face    : 2-cell
- Edge    : 1-cell
- Vertex  : 0-cell
- Element : 3-cell in 3D, 2-cell in 2D
- Stratum : set of all k-cells for a given k

### Working Definition of DM
We think of a DM as a combination of up to four other concepts:

1. Topology (required), including parallel decomposition
2. An embedding/immersion of this topology (coordinates)
3. A default "Section" - what the primary fields living on the DM are. 


The fields in 3. are assumed to be strongly interacting, implying that it's natural to store them interleaved. 
`DM` object themselves are lightweight - heavy data are pointed to by references. Thus, if non-interleaved fields are desired,
it is natural to use a `DMComposite` object comprised of several `DM`s which may refer to the same topology and coordinates, with different sections/fields.

In the case of `DMPlex`, 1. (topology) is in 1-1 correspondence with the `DMPlex`-specific information in `DM_Plex` (pointed to by the `data` field in `_p_DM`). This data in turn contains pointers to heavy topology data and a reference count. This means that two `DMPlex` objects with the same topology can both simply point to the same `DM_Plex` object.

This is not the case for `DMDA`: while all of the topology (the grid sizes) are indeed captured by `DM_DA`, additional information is also included, in particular the DOF specifications. The DOF specification would more naturally be included with the section information (3.). The stencil information might be associated with a particular operator, 
though might naturally considered as part of the topology specification, if the concept of topology is extended to include a parallel decomposition.

## Design Considerations

### PETSc design
- Maintain parallels to `DMDA` and `DMPlex`
- Obeys the working definition above
- Provide interfaces that "make MPI invisible" as much as possible.

### Scope and generality/optimization
- Allow for highly-efficient operations for the applications of primary interest (narrow-stencil Stokes operators)
- Reduce code duplication where possible, but focus on 1-,2, and 3-dimensional cases.

### Ease of use
- Provide interfaces that minimize the number of indexing errors possible from the user

### Extensibility
- Allow for a pathway to intelligent cache-blocking. (Determinining the size of the blocks is out of scope here, but maybe we can provide some of the index-twiddling to make this easy enough for people to use.)
- Allow for future extension using as a base `DM` for `DMForest`
- Allow for future custom stencil types (higher-order FVM schemes?)

## Indexing Convention

We refer to all points by the x,y[,z] indices of an element, and "half index" offsets in each of the three dimensions.
These offsets can be positive or negative, giving non-unique representations of non-primal cells (one for each primal
cell intersected).

```
         +-----------+-----------+
         |           |           |
         |           |           |  B: (i,j|0,+) or (i,j+1|0,-)
         | i,j+1|0,0 |           D  
         |           |           |  C: (i,j|+,+),(i,j+1|+,-),(i+1,j|-,+) or (i+1,j+1,-,-)
         |           |           |
         +-----B-----C-----------+  D: (i+1,j+1|+,0), or (i+2,j+1|-,0)
         |           |           |
         |           |           |
 x       |  i,j|0,0  |           |
 ^       |           |           |
 |       |           |           |
 +--> y  +-----------+-----------+
```

## DOF ordering
Each point is associated with a cell by choosing it's "-" representation. 
Global "Natural" ordering within each element is then done in "x fast" ordering, with - < 0.
```
0: (0,0|-,-)
1: (0,0|0,-)
2: (0,0|-,0)
3: (0,0|0,0)
...
```

Note that on the right/top/back boundaries, numbering is respect to a fictitious cell.
This ordering is used in an extension to `MatStencil`. This allows the user to simply
think in terms of which element (in global numbering) and which boundary, as opposed to having to remember
a convention for index numbering of face/edge/vertex degrees of freedom.

This also corresponds to a general convention of ordering the lowest-dimensional object first
(i.e. vertices, edges, faces, 3-cells).

### 2D example
```
20 - 21 - 22 - 23 - 24
|         |         |
12   13   16   17   19
|         |         |
10 - 11 --14 - 15 - 18
|         |         |
2    3    6    7    9
|         |         |
0 -- 1 -- 4 -- 5 -- 8
```
Note that not all strata must be included in a `DMStag` object (that is,
the domain on Sections need not include all of these).
For instance, our standard Stokes solve might use two `DMStag` objects, one for velocity and pressure fields,
and one to hold density and viscosity fields.

Unknowns (pressure and velocity)

```
    +--- 14 --+--- 15 --+
    |         |         |
    8    9    11   12   13
    |         |         |
    +--- 7 ---+--- 10 --+
    |         |         |
    1    2    4    5    6
    |         |         |
    +--- 0 ---+--- 3 ---+
```


Vertex- and Cell-based material parameters
```
    10 ------ 11 ------ 12
    |         |         |
    |    6    |    8    |
    |         |         |
    5 ------- 7 ------- 9
    |         |         |
    |    1    |    3    |
    |         |         |
    0 --------2 ------- 4
```

### Parallel Decomposition
In parallel, partitioning is always done by cell. As with `DMDA`, the decomposition must be 
a "product decomposition". Fictious extra cells are associated with the last rank in each direction.

## Compatibility

We deem two (or more) `DMStag` objects "compatible" if 
1. they are different "views" of the same abstract staggered grid (the same grid dimension and sizes).
2. they have identical parallel decompositions of the underlying cells

This amounts to checking that the dimensions on each rank are identical.
```
PetscErrorCode DMStagCheckCompatibility(dm dmstag,dm dmstag2,...); /* varargs won't work in Fortran */
```
We might also consider introducing a function for this into the general `DM` API.

This allows us to formalize the process of iterating over several DMs together
```
PetscErrorCode DMStagGetCornersCompatible(DM dmstag0,PetscInt *x,PetscInt *y,...,PetscInt *n,PescInt *m,...,DM dmstag1,...);  
PetscErrorCode DMStagVecGetArraysCompatible(DM dmstag0,**PetscScalar arr0,PetscBool readonly0,DM dmstag1,...);
```

Note that we don't inlude anything about the coordinates in the above definition of compatibility.
In most cases, one would only consider two `DMStag` objects to be compatible if their coordinates
aligned on the overlapping strata.

## Stencil and Points Description

`DMStag` has a flag for each stratum to determine if its points are included.

`DMStag` accepts a number of degrees of freedom associated with points in each stratum.
This may be 0, but typically one would rather clone the DM and turn the flag
off for that stratum entirely. Thus, in our creation routines below we
simply set these flags by checking if the corresponding DOF count is zero.
The user could turn them back on with the `DMStag` API before SetUp.

We define stencil types, which may apply to point in one, some, or all strata.

- `DMSTAG_STENCIL_NONE` (placeholder when stratum is not active, or dof is 0)
- `DMSTAG_STENCIL_SELF` (update from the same point)
- `DMSTAG_STENCIL_SELF_STAR` (star in the same stratum)
- `DMSTAG_STENCIL_SELF_BOX` (box in the same stratum)
- `DMSTAG_STENCIL_BOUNDARY` (neighbors in the next lower-dimensional stratum)
- `DMDSTAG_STENCIL_EXTERIOR` (neighbors in the next higher-dimensional stratum)
- `DMSTAG_STENCIL_MOMENTUM`  (9-point simple Stokes FV momentum stencil, bad name)
- `DMDSTAG_STENCIL_MOMENTUM_COEFFICIENT` (stencil required for simple Stokes FV coefficients on an edge, bad name)

## Ghosts and Boundaries

As with `DMPlex`, `DMLabel` objects will be used to defined subset of points.

## Creation

Creation Routines mirror those for `DMDA`.

```
PetscErrorCode DMStagCreate2d(
  MPI_Comm comm,
  DMBoundaryType bx, DMBoundaryType by,                   // boundary types
  DMStagStencilType stencilCell,                          // Stencil information
  DMStagStencilType stencilEdge,                          
  DMStagStencilType stencilVertex,                        
  PetscInt M,PetscInt N,                                  // global sizes
  PetscInt m,PetscInt n,                                  // local sizes
  PetscInt dofVertex, PetscInt dofEdge, PetscInt dofElement, // DOF information (0 turns off a point type)
  DM *dm
)
```

## Conversion

Creating derivative `DMStag` objects will be done with the Clone/Set/SetUp paradigm,
or with `DMGetSubDM()`.

Useful potential conversion routines include

- Conversion from `DMStag` to `DMPlex`
- Extraction of a single-point-type subgrid as a `DMDA`

## Helpers

It might be useful to provide a function to interpolate values of field from one stratum to another (even though this is typically not the best thing to do). This could be useful for quickly propagating fields to a collocated grid for plotting, or for moving available parameter fields around.

## Members

As with `DMDA`, we keep extra, perhaps-unused, fields corresponding to the maximum
dimensionality of the grid (likely to remain at 3).

```
#define MAX_GRID_DIM 3
#define MAX_STRATA MAX_GRID_DIM+1
typedef struct {
  PetscInt N[MAX_GRID_DIM];            /* Global dimensions              */
  PetscInt n[MAX_GRID_DIM];            /* Local dimensions               */
  PetscInt dof[MAX_STRATA];            /* dof per point for each stratum */
  PetscBool stratumActive[MAX_STRATA]; /* which strata are active        */
  DMStagStencilType[MAX_STRATA];       /* Stencil type for each stratum  */
  ...
} DMStag;
#undef MAX_GRID_DIM
#undef MAX_STRATA

```

