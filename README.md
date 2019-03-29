![StagBL](documentation/resources/logo/logo_half.pdf)

## About

StagBL is a C library designed to allow optimized, massively-parallel Stokes solvers
for geodynamics application codes based on finite-volume methods on regular,
orthogonal grids, usually coupled to a particle-based advection scheme.

It aims to be as lightweight as possible while still providing the flexibility
and extensibility required for scientific application codes. This accomplished
with careful design and interfaces to powerful external libraries like
[PETSc](https://www.mcs.anl.gov/petsc).

## Dependencies

* MPI is required.
* PETSc is currently required (but will not always be).
* Python is required for the configuration process.

## Quickstart

First, you need a working branch of PETSc's master branch, configured with
`--download-suitesparse`.

    git clone -b master https://bitbucket.org/petsc/petsc petsc-master

Then, from this directory,

    export PETSC_DIR=yyy
    export PETSC_ARCH=xxx
    ./configure.py         # follow instructions to make
    cd demos/2d
    make
    ./stagbldemo2 -pc_type lu -pc_factor_mat_solver_type umfpack
    paraview &             # open out_element.vtr

![stagbl2ddemo quickstart](documentation/resources/stagbldemo2d_quickstart.pdf)

## Support
Development of StagBL is supported by the [Platform for Advanced Scientific Computing](https://www.pasc-ch.org).
