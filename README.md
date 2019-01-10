![StagBL Logo](documentation/resources/logo/logo.pdf)
# StagBL

StagBL is in the initial heavy development phase. Everything and anything may change.
See the Trello Board on the Bitbucket repository.

## Dependencies

* MPI is required.
* PETSc is currently required (but will not always be).
* Python is required for the configuration process.

## Quickstart

First, you need a working branch of PETSc which includes `DMStag`, currently "master":

    git clone -b master https://bitbucket.org/petsc/petsc petsc-master

After it is configured,

    export PETSC_DIR=yyy
    export PETSC_ARCH=xxx
    ./configure.py
    cd demos/2d
    make
    ./stagbldemo2d                     # -pc_type lu -pc_factor_mat_solver_type umfpack
    paraview & # open out_element.vtr

![stagbl2ddemo quickstart](documentation/resources/stagbldemo2d_quickstart.pdf)
