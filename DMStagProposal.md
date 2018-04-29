DMStag - Component Proposal
---------------------------

(See psanan/dmstag branch at https://bitbucket.org/psanan/petsc)

`DMStag` is a proposed implementation of PETSc's `DM` abstract class.

It represents a topologically Cartesian cell complex based on a rectangular grid of line segments, quadrilaterals, or hexahedra.
It provides routines to manipulate interacting fields with respect to these top-level element ands lower-dimensional cells in the complex.

The intended use is for finite volume or DEC schemes, on regular grids.
It is intended to be intermediate between `DMDA` and `DMPlex`; it represents a structured grid like `DMDA`, but supports different strata like `DMPlex`.

It will also support a notion of "structured coordinates," allowing for efficiency associated with knowing that coordinates are regularly-spaced in one or more dimensions, or describable as a product of 1- or 2-dimensional coordinate grids.

## Comparison with DMDA

It is illuminating to consider how the functionality of DMStag can be achieved using DMDA.

Most staggered-grid codes in PETSc thus far have opted to either use
1. A single DMDA with "dummy" points
2. Multiple DMDA objects

The first approach is practical but saddles the user with remembering a convention about node numbering, and that certain points on the boundary are unused "dummy" points.

The second requires care in ensuring that parallel decompositions and numberings are used correctly, and does not allow for interlaced storage of data.


Since the cell complex represented by DMDA is included in those which can be represented by a DMPlex object, it may of course be represented this way. However, one expects to find some added performance and simplicity in restricting to a much-more-structured case.

## Terminology

- Point        : a cell of any dimension
- Vertex       : 0-cell
- Element      : 3-cell in 3D, 2-cell in 2D, 1-cell in 1D
- Stratum      : set of all k-cells for a given k
- Entry        : a single entry in a vector representing one or more fields on the complex
- Ghost        : describes an extra point or degree of freedom used locally, in addition to those which correspond to global degrees of freedom stored on the current rank
- Dummy point  : a ghost point used for padding purposes, which does not participate in global<-->local mappings. These can either on the right/top/front of the physical domain, or in the corners of interior subdomains when using a "star" stencil.

### Working Definition of DM
We think of a `DM` as a combination of up to three other concepts:

1. Topology (required), including parallel decomposition
2. An embedding/immersion of this topology (coordinates)
3. A "Section" - what the primary fields living on the DM are.

"Topology" includes both "local" and "global" topology and the maps between them. In the case of `DMDA`, we consider the stencil description to be part of the topology, as it is used to define the local spaces.

The fields in the section are assumed to be strongly interacting, implying that it's natural to store them interleaved.
`DM` objects are lightweight - heavy data are pointed to by references. Thus, if non-interleaved fields are desired,
one may use a `DMComposite` object comprised of several `DM`s which may refer to the same topology and coordinates, with different sections/fields.

In the case of `DMPlex`, Topology is the only thing defined by the `DMPlex`-specific information in `DM_Plex` (pointed to by the `data` field in `_p_DM`). This data in turn contains a reference count and pointers to heavy topology data. This means that two `DMPlex` objects with the same topology can share a pointer to the same `DM_Plex` object.
This is not the case for `DMDA`: while all of the topology (the grid sizes) are indeed captured by `DM_DA`, additional information is also included, in particular the DOF specifications. The DOF specification would more naturally be included with the section information.

It should be noted that one will often work with multiple DM objects corresponding to the same simulation domain. Thus, we introduce the notion of "compatible" DM objects as those with identical topology (including parallel decomposition).

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

Points are referred to by specifying an element with i,j,k indices, and then by specifying up to three "half element" offsets. These can be positive or negative in each of the three dimensions (LEFT, RIGHT, UP, DOWN, FRONT, BACK). Note that this gives a non-unique representation of all points except those centered in elements.

```
         +-----------+-----------+
         |           |           |
         |           |           |  B: i,j,UP or i,j+1,DOWN
         |   i,j+1   |  i+1,j+1  |
         |           |           |  C: i,j,UP,RIGHT or i,j+1,DOWN,RIGHT
         |           |           |       or i+1,j,LEFT,UP or i+1,j+1,LEFT,DOWN
         +-----B-----C-----------+
         |           |           |
         |           |           |
 y       |    i,j    |   i+1,j   |
 ^       |           |           |
 |       |           |           |
 +--> x  +-----------+-----------+
```

## DOF ordering
Elements are ordered "x fast", e.g. (0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(0,1,1),(1,1,1) for a 3d 2x2x2 mesh.
Each point is associated with an element by choosing its representation involving only LEFT,DOWN,and BACK
Global "Natural" ordering within each element is then done in "x fast" ordering, e.g. for the first element
```
0: 0,0,BACK,DOWN,LEFT
1: 0,0,BACK,DOWN
2: 0,0,BACK,LEFT
3: 0,0,BACK
4: 0,0,DOWN,LEFT
5: 0,0,DOWN
6: 0,0,LEFT
7: 0,0
...
```

Note that on the right/top/front boundaries, numbering is respect to a dummy element.
This ordering is used in an extension to `MatStencil`. This allows the user to simply
think in terms of which element (in global numbering) and which boundary, as opposed to having to remember
a convention for index numbering of face/edge/vertex degrees of freedom.

(Note: more-verbose names like `PETSC_STENCIL_LEFT` will likely be used to avoid clashes in the global namespace).


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

We deem two (or more) `DMStag` objects "compatible" if they have the same topology (including parallel decomposition). That is, they have the same sizes and live on the same communicator. This in particular asserts that they can be safely iterated over together, using element numbering from either of them.

## Structured coordinate description for orthogonal grids

DMStag features optimized additional data structures for grid coordinates, beyond those available with DM.
There, a single coordinate DM is used for all coordinate information. With DMStag (as opposed to the other available option of defining special DM implementations for this purpose), we store a flag signifying whether to use "structured" coordinates, along with a coordinate type for each dimension. This includes "uniform" and "DMStag 1D". Data is also stored giving a (global) minimum and maximum real value for each dimension, a spacing (redundant), and a pointer to a DM. Functions in the DMStag API can be used to efficiently access these coordinates. In the

## Intuitive indexing

"Raw" access to vectors is of course still available, for expert users.

In addition, two other access methods are provided, one which is likely more performant, but requires
more assumptions from the user, and another which is less error-prone but likely will involve more overhead.

### Method 1

Here, the user must accept the convention about which element degress of freedom are attached to,
"under the hood". Namely, that a degree of freedom which is not element-centered is grouped
with the element above, right, or in front of it, and that ordering within an element follows the same "S shape". This, and the fact that local
vectors are padded with ghosts, allows access to vectors much like a DMDA with multiple degrees of freedom.
This simply uses e.g. `VecGetValues2d()` internally.

### Method 2

The method to assemble matrices with stencils (`MatStencil`,`MatSetValuesStencil`) is generalized
to allow set and get to all entries bordering a given element, and to allow usage with `Vec`s
as well as `Mat`s.

The approach is to create auxilary arrays `PetscScalar vals[]` and `StagStencil[] entries`,
which are populated with values and locations, e.g.
```
vals[idx]        = 1.0;
entries[idx].i   = 3;
entries[idx].j   = 5;
entries[idx].loc = LEFT;
...
ierr = DMStagVecSetValuesStencil(v,vals,INSERT_VALUES,entries);CHKERRQ(ierr);
...
ierr = VecAssemblyBegin(v);CHKERRQ(ierr);
ierr = VecAssemblyEnd(v);CHKERRQ(ierr);
```

It's important to note a potential source of errors for the user here, due to a simplification in our implementation. There is (currently) an asymmetry about which dof are actually available for given stencil width. With a width of 1 and an "edge" dof (a vertex, in 1d), we know that adjacent "interior" (and element in 1d) dof will be available, but we do NOT know that the point 1 element away will be available. On the left it will be, but on the right it will not. This is something of an artifact, and a more advanced implementation would NOT correctly fill the value on the left, though it would allocate a space for it so we could iterate in a blocked fashion.
```
x - | x - x - x - x - | x -
                    a     b
                    something is probably wrong if your stencil has b depending on a...
```

## Stencil and Points Description

`DMStag` accepts a number of degrees of freedom associated with points in each stratum.
0 points implies that the stratum is inactive.

The creation routines accept a single element-wise stencil argument
which determines the sizes of the ghost regions in the local representations.

Operations should eventually be expanded to add more specific stencil information which can be used
to create operators and perhaps to reduce halo exchange information.

This will not be implemented until an actual performance penalty is demonstrated
resulting from sendin complete elements-worth of data in ghost regions.

## Creation

Creation Routines mirror those for `DMDA`.

(See current working implementation)

## Conversion

Creating derivative `DMStag` objects is accomplished with a new function `DMStagCreateReducedDMStag()` which (despite the name) creates a compatible `DMStag` object with different numbers of dof on each point type, ignoring or zeroing values present in only one of the source or destination `DMStag`.

## Members

As with `DMDA`, we keep extra, perhaps-unused, fields corresponding to the maximum
dimensionality of the grid (likely to remain at 3).

(See current working implemetation)

## Examples

A small set of examples/tests is key to the PETSc contribution. A small set of examples should cover the new functionality.

* `stag_ex1` : for unit tests. Don't do anything useful but simply exercise the API, controlling with command line options
* `stag_ex2` : a 1-D toy problem (u''(x) = f(x) in mixed form)
* `stag_ex3` : a 2-D problem (isoviscous Stokes, MMS)
* `stag_ex4` : a 3-D problem, similar to the previous example
* `stag_ex5` : a 2-D linear Stokes solver, analogous to KSP tutorial ex70
