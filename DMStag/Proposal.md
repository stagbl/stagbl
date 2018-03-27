DMStag - Component Proposal
---------------------------

`DMStag` is a proposed implementation of PETSc's `DM` abstract class.

It represents a topologically Cartesian cell complex based on a rectangular grid of line segments, quadrilaterals, or hexahedra.
It provides routines to manipulate interacting fields with respect to these top-level element ands lower-dimensional cells in the complex.

The intended use is for finite volume or DEC schemes, on regular grids.
It is intended to be intermediate between `DMDA` and `DMPlex`; it represents a structured grid like `DMDA`, but supports different strata like `DMPlex`.

It will also support a notion of "compressed coordinates," allowing for efficiency associated with knowing that coordinates are regularly-spaced in one or more dimensions, or describable as a product of 1- or 2-dimensional coordinate grids.

## Comparison with DMDA 

It is illuminating to consider how the functionality of DMStag can be achieved using DMDA.

Most staggered-grid codes in PETSc thus far have opted to either use 
1. A single DMDA with "dummy" points
2. Multiple DMDA objects

The first approach is practical but saddles the user with remembering a convention about node numbering, and that certain points on the boundary are unused "dummy" points.

The second requires care in ensuring that parallel decompositions and numberings are used correctly, and does not allow for interlaced storage of data.

It should be noted that one will often work with multiple DMStag objects. A DM being the combination of a (parallel decompositon of a) topology, an embedding, and a section on a subset of that topology, we introduce the notion of "compatible" DMStag objects being those wherein the topology and its parallel decomposition is identical, and where embedding on all active points is equal (though this check may be omitted for performance reasons).

Since the cell complex represented by DMDA is included in those which can be represented by a DMPlex object, it may of course be represented this way. However, one expects to find some added performance and simplicity in restricting to a much-more-structured case.

## Terminology

- Point   : a cell of any dimension
- Face    : 2-cell
- Edge    : 1-cell
- Vertex  : 0-cell
- Element : 3-cell in 3D, 2-cell in 2D, 1-cell in 1D
- Stratum : set of all k-cells for a given k
- Entry   : a single entry in a vector representing one or more fields on the complex
- Ghost   : describes an extra point or degree of freedom corresponding to one stored on a neighboring rank (which may be defined by periodicity) in the global decomposition
- Dummy   : describes an extra point or degree of freedom in the local representation, used on the top/right/front of the domain, to give a representation with an equal number of points of each type

### Working Definition of DM
We think of a `DM` as a combination of up to three other concepts:

1. Topology (required), including parallel decomposition
2. An embedding/immersion of this topology (coordinates)
3. A default "Section" - what the primary fields living on the DM are. 

"Topology" includes both "local" and "global" topology and the maps between them. In the case of `DMDA`, we consider the stencil description to be part of the topology, as it is used to define the local spaces.

The fields in the default section are assumed to be strongly interacting, implying that it's natural to store them interleaved. 
`DM` objects are lightweight - heavy data are pointed to by references. Thus, if non-interleaved fields are desired,
one may use a `DMComposite` object comprised of several `DM`s which may refer to the same topology and coordinates, with different sections/fields.

In the case of `DMPlex`, Topology is the only thing defined by the `DMPlex`-specific information in `DM_Plex` (pointed to by the `data` field in `_p_DM`). This data in turn contains a reference count and pointers to heavy topology data. This means that two `DMPlex` objects with the same topology can share a pointer to the same `DM_Plex` object.
This is not the case for `DMDA`: while all of the topology (the grid sizes) are indeed captured by `DM_DA`, additional information is also included, in particular the DOF specifications. The DOF specification would more naturally be included with the default section information.  

## Design Considerations

### PETSc design

- Maintain parallels to `DMDA` and `DMPlex`
- Obeys the working definition above
- Provide interfaces that "make MPI invisible" as much as possible.
- Prioritize simple, efficient, easily-debuggable implementation

### Scope and generality/optimization

- Allow for highly-efficient operations for the applications of primary interest (narrow-stencil Stokes operators)
- Reduce code duplication where possible, but focus on 1-,2, and 3-dimensional cases within a single object. Prioritize efficiency in the 3D case.

### Ease of use

- Provide interfaces that reduce the number of indexing errors possible from the user

### Extensibility

- Allow for a pathway to intelligent cache-blocking. (Determining the size of the blocks is out of scope here, but maybe we can provide some of the index-twiddling to make this easy enough for people to use.)
- Allow for future extension using as a base `DM` for `DMForest`
- Allow for future custom stencil types (in particular, higher-order FVM schemes)
- Allow for interaction with `DMSwarm` or other PIC/MIC schemes
- Allow for the introduction of `PetscSection`/`PetscSF` to define fields, as in `DMPlex` (and, undocumented it seems, in `DMDA`).

## Indexing Convention

We refer to all points by the x,y[,z] indices of an element, and "half index" offsets in each of the three dimensions.
These offsets can be positive or negative, giving non-unique representations of non-primal cells (one for each primal
cell (element) intersected).

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
 y       |  i,j|0,0  |           |
 ^       |           |           |
 |       |           |           |
 +--> x  +-----------+-----------+
```

## DOF ordering
Each point is associated with an element by choosing it's "-" representation. 
Global "Natural" ordering within each element is then done in "x fast" ordering, with - < 0.
```
0: (0,0|-,-)
1: (0,0|0,-)
2: (0,0|-,0)
3: (0,0|0,0)
...
```

Note that on the right/top/front boundaries, numbering is respect to a fictitious element.
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


Vertex- and Element-based material parameters
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
In parallel, partitioning is always done by element, associating each lower-dimensional cell with the element it is below/left of/behind. As with `DMDA`, the decomposition must be a "product decomposition". Partial elements with extra dummy points are associated with the last rank in each direction.
Global ("PETSc") ordering of degrees of freedom is defined by ordering ranks in "x-fast" order starting from the bottom back left,
ordering points locally by the same scheme as above, ordering dof sequentially within each point. For example, consider a 2D 3x3 grid, with 1 dof on each vertex and two dof on each element, decomposed on 4 ranks in a 2x2 grid.

Global numbering of dof:
```
 [rank 2]                           :  [rank 3]
                                    :
    26 ------ 27 ------             :     32 ------ 33
    |         |                     :     |         | 
    |  21,22  |  24,25              :     |  29,30  | 
    |         |                     :     |         | 
    20------- 23-------             :     28------- 31
                                    :
....................................:..............................
                                    :
 [rank 0]                           :   [rank 1]
                                    :
    |         |                     :     |         |
    |  7,8    |  10,11              :     |  17,18  |
    |         |                     :     |         |
    6 ------- 9 -------             :     16------- 19
    |         |                     :     |         |
    |  1,2    |  4,5                :     |  13,14  |
    |         |                     :     |         |
    0 ------- 3 -------             :     12 ------ 15
```

Local numbering of dof:
```

 [rank 2]                           :  [rank 3]
                                    :
    |         |         |           :     |         |        
    |  10,11  |  13,14  |  16,17    :     |  7,8    |  10,11 
    |         |         |           :     |         |        
    9 ------- 12 ------ 15 -----    :     6 ------- 9 -------
    |         |         |           :     |         |        
    |  1,2    |  4,5    |  7,8      :     |  1,2    |  4,5   
    |         |         |           :     |         |        
    0 ------- 3 ------- 6 ------    :     0 ------- 3 -------
                                    :
....................................:..............................
                                    :
  [rank 0]                          :  [rank 1]
                                    :
    |         |         |           :     |         |        
    |  19,20  |  22,23  |  25,26    :     |  13,14  |  16,17 
    |         |         |           :     |         |        
    18 ------ 21 ------ 24 -----    :     12 ------ 15 ------
    |         |         |           :     |         |        
    |  10,11  |  13,14  |  16,17    :     |  7,8    |  10,11 
    |         |         |           :     |         |        
    9 ------- 12 ------ 15 -----    :     6 ------- 9 -------
    |         |         |           :     |         |        
    |  1,2    |  4,5    |  7,8      :     |  1,2    |  4,5   
    |         |         |           :     |         |        
    0 ------- 3 ------- 6 ------    :     0 ------- 3 -------

```

As with `DMDA`, "natural" ordering is defined as the global ordering in the 1-rank case.

## Compatibility

We deem two (or more) `DMStag` objects "compatible" if 

1. they are different "views" of the same abstract staggered grid (the same grid dimension and sizes).
2. they have identical local representations (the same parallel decomposition and ghost region sizes).
3. If coordinates are assigned, they agree. That is, two DMStag objects are deemed incompatible if they assign different coordinates to the same point.

The last point can be further clarified by stating
- some condition for "different coordinates" must be specified. The default might reasonably be bitwise equality, as coordinates are likely to be copied or generated by an identical process. 
- this check might be skipped by default, for performance purposes, except in the cases where such a check can be performed in O(1) time for regular coordinates.

```
PetscErrorCode DMStagCheckCompatibility(dm dmstag,PetscInt ndms,DM dm[],PetscBool *isCompatible);
```
We might also consider introducing a function for this into the general `DM` API.

This allows us to formalize the process of iterating over several DMs together
```
PetscErrorCode DMStagGetCornersCompatible(DM dmstag0,PetscInt *x,PetscInt *y,...,PetscInt *n,PescInt *m,...,DM dmstag1,...);  
PetscErrorCode DMStagVecGetArraysCompatible(DM dmstag0,**PetscScalar arr0,PetscBool readonly0,DM dmstag1,...);
```

Note that we don't include anything about the coordinates or section in the above definition of compatibility.

## Compressed coordinate description for orthogonal grids

TODO: further formalize
Inspired by LaMEM, introduce an optional, compressed way to represent coordinates for an
orthogonal grid. In each dimension, coordinates may be declared as uniform between a given minimum and maximum, 
or by a 1D array. 

## Intuitive indexing

TODO: further formalize
An extension to `MatStencil` and `MatSetValuesStencil` should be introduced,
to allow a natural way to work values corresponding to lower-dimensional cells.
For getting and setting vector entries, this may entail introducing a second,
perhaps less-efficient method to index with the help of an intermediate array
of stencil values.

## Stencil and Points Description

`DMStag` has a flag for each stratum to determine if its points are included.

`DMStag` accepts a number of degrees of freedom associated with points in each stratum.
There might be some argument for allowing this to be zero, but as in most cases one would 
simply create a new `DMStag` object.
Thus, in our creation routines below we
simply set these flags by checking if the corresponding dof count is zero.
The user could turn them back on with the `DMStag` API before SetUp.

The creation routines below accept a single element-wise stencil argument
which determines the sizes of the ghost regions in the local representations.

Operations should eventually be expanded to add more specific stencil information which can be used
to create operators and perhaps to reduce halo exchange information.
We define stencil types, which may apply to point in one, some, or all strata.

- `DMSTAG_STENCIL_NONE` (placeholder when stratum is not active, or dof is 0)
- `DMSTAG_STENCIL_SELF` (update from the same point)
- `DMSTAG_STENCIL_SELF_STAR` (star in the same stratum)
- `DMSTAG_STENCIL_SELF_BOX` (box in the same stratum)
- `DMSTAG_STENCIL_BOUNDARY` (neighbors in the next lower-dimensional stratum)
- `DMDSTAG_STENCIL_EXTERIOR` (neighbors in the next higher-dimensional stratum)
- `DMSTAG_STENCIL_MOMENTUM`  (11-point simple Stokes FV momentum stencil, bad name)
- `DMDSTAG_STENCIL_MOMENTUM_COEFFICIENT` (stencil required for simple Stokes FV coefficients on an edge, bad name)

## Creation

Creation Routines mirror those for `DMDA`.

```
PetscErrorCode DMStagCreate2d(
  MPI_Comm comm,
  DMBoundaryType bx, DMBoundaryType by,                      // boundary types
  DMStagGhostStencilType stencil,                            // element-wise stencil type
  PetscInt stencilWidth,                                     // element-wise stencil width
  PetscInt M,PetscInt N,                                     // global sizes
  PetscInt m,PetscInt n,                                     // local sizes
  PetscInt dofVertex, PetscInt dofEdge, PetscInt dofElement, // DOF information (0 turns off a point type)
  DM *dm
)
```

## Conversion

Creating derivative `DMStag` objects will be done with the Clone/Set/SetUp paradigm,
or with `DMCreateSubDM()`.

## Members

As with `DMDA`, we keep extra, perhaps-unused, fields corresponding to the maximum
dimensionality of the grid (likely to remain at 3).

Members for a current prototype include the following:
```
#define DMSTAG_MAX_DIM 3
#define DMSTAG_MAX_STRATA MAX_DIM+1
typedef struct {
  PetscInt               dim;
  PetscInt               N[DMSTAG_MAX_DIM];                /* Global dimensions (elements)    */
  PetscInt               n[DMSTAG_MAX_DIM];                /* Local dimensions (elements)     */
  PetscInt               nghost[DMSTAG_MAX_DIM];           /* Local dimensions (with ghosts)  */
  PetscInt               start[DMSTAG_MAX_DIM];            /* First element number            */
  PetscInt               startGhost[DMSTAG_MAX_DIM];       /* First element number (ghosted)  */
  PetscMPIInt            proc[DMSTAG_MAX_DIM];             /* Location in processor grid      */
  PetscBool              lastproc[DMSTAG_MAX_DIM];         /* Last proc in this dim?          */ 
  PetscBool              firstproc[DMSTAG_MAX_DIM];        /* First proc in this dim?         */ 
  PetscMPIInt            nprocs[DMSTAG_MAX_DIM];           /* Procs in each direction         */
  PetscInt               dof[DMSTAG_MAX_STRATA];           /* dof per point for each stratum  */
  PetscBool              stratumActive[DMSTAG_MAX_STRATA]; 
  DMStagGhostStencilType ghostStencil;                     /* element-wise ghost stencil      */
  PetscInt               ghostStencilWidth;                /* elementwise ghost width         */
  DMBoundaryType         boundaryType[DMSTAG_MAX_DIM];
  VecScatter             gton;                             /* Global  --> Natural             */
  VecScatter             gtol;                             /* Global  --> Local               */
  // Potential : additional, more-specific stencil information
} DM_Stag;

```

## Examples

A small set of examples/tests is key to the PETSc contribution. A small set of examples should cover the new functionality (and in the future we can confirm this with the automated gcov information in the nightlies).

1. An example for the basic operations: 
 - creation
 - viewing, 
 - generating a vector 
 - iterating over that vector (AoS)
 - creating compatible DMStags with different active strata and/or default sections
 - checking for compatibility
 - simultaneous iteration (SoA)
 - destroying
2. An example for using multiple grids and checking compatibility.
3. A simple operator construction (isoviscous Stokes, eta=1, using MMS), solution, and converge rate check.
4. A linear Stokes solver, modelled on KSP tutorial example 70, solving a common benchmark (or several), e.g. solck, solkx, Rayleigh-Benard (compute Nusselt number), viscous sinker, etc.
5. Additional tests of individual operations, as needed for coverage or for better unit testing
