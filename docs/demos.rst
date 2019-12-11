StagBL Demos
------------

StagBLDemo2d and StagBLDemo3d are two mini-apps designed for several purposes:

1. To provide general examples of usage of StagBL. Whenever possible, these demos will be used
   as examples/tutorials.
2. "The path to performance," To demonstrate how programs similar to those in Gerya's textbook may be written with StagBL
3. To reproduce simple benchmarks with StagBL
4. To provide application testing for StagBL
5. To provide platforms for performance testing, particular for emerging architectures and HPC systems.

These codes are PETSc C codes, sharing some code between 2d and 3d demos, as
this allows for more concise demos, but StagBL can be called from any C/C++
application code, including those with their own MPI-based logic. (Usage from
Fortran is also possible with the use of some wrappers which are not currently
included with StagBL itself. Please contact the authors if this is of
interest.)

.. image:: resources/3d_sinker_box_20.png

DMStag Demos
------------

In addition to these demos, there are several tutorial examples included with
PETSc, that demonstrate direct usage of DMStag, the object providing the main
implementation of the ``StagBLGrid`` object. See ``src/dm/impls/stag/examples/tutorials``
in the PETSc source tree.

This object can be used for tasks beyond the scope of StagBL, as it provides
a generic interface for working with 1- to 3-dimensional, logically regular
cell complexes, in particular those corresponding to orthogonal grids.

For example, `DMStag tutorial ex6 <https://bitbucket.org/psanan/petsc/src/6f35e31b9f2989e6fe59ddc38ff726d76adaefc9/src/dm/impls/stag/examples/tutorials/ex6.c?at=psanan%2Fstagbl-working-base>`__ shows how to use DMStag to simulate seismic waves.

.. code-block:: bash

      rm -f *.vtr && ./ex6 -nsteps 1000 -stag_grid_x 1000 -stag_grid_y 1000 -dt 4e-4

The resulting `.vtr` file for frame 485 can be visualized in Paraview. This is the y velocity:

.. image:: resources/dmstag_ex6_yvel_485.png



  
