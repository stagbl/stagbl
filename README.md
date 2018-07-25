![StagBL Logo](documentation/resources/logo/logo.pdf)
# StagBL

StagBL is in the initial heavy development phase. Everything and anything may change.
See the Trello Board on the Bitbucket repository.

## Dependencies

* MPI is required.
* PETSc is currently required (but will not always be).
* Python is required for the configuration process.

## Quickstart

First, you need a working branch of PETSc which includes `DMStag`, obtainable with

    git clone -b psanan/dmstag https://bitbucket.org/psanan/petsc

After it is configured,

    export PETSC_DIR=yyy
    export PETSC_ARCH=xxx
    ./configure.py
    cd demos/StagBLDemo2d
    make
    ./stagbldemo2d
    paraview & # open out_element.vtr

![stagbl2ddemo quickstart](documentation/resources/stagbldemo2d_quickstart.pdf)
