![StagBL Logo](documentation/resources/logo/logo_half.pdf)
# StagBL

StagBL is a C library designed to allow optimized, massively-parallel solvers
for geodynamics application codes based on finite-volume methods on regular,
orthogonal grids, usually coupled to a particle-based advection scheme.

It aims to be as lightweight as possible while still providing access to modern
solvers.

StagBL is in the initial heavy development phase. Everything and anything may change.
See the Trello Board on the Bitbucket repository.

## Dependencies

* MPI is required.
* PETSc is currently required (but will not always be).
* Python is required for the configuration process.

## Quickstart

First, you need a working branch of PETSc's master branch

    git clone -b master https://bitbucket.org/petsc/petsc petsc-master

After it is configured, with `--download-suitesparse`,

    export PETSC_DIR=yyy
    export PETSC_ARCH=xxx
    ./configure.py         # follow instructions to make
    cd demos/2d
    make
    ./stagbldemo2 -pc_type lu -pc_factor_mat_solver_type umfpack
    paraview &             # open out_element.vtr

![stagbl2ddemo quickstart](documentation/resources/stagbldemo2d_quickstart.pdf)
