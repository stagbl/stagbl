![StagBL Logo](documentation/resources/logo/logo.pdf)
# StagBL

StagBL is in the initial heavy development phase. Everything and anything may change.
See the Trello Board on the Bitbucket repository.

## Quickstart

    export PETSC_DIR=yyy
    export PETSC_ARCH=xxx
    ./configure.py
    cd demos/StagBLDemo2d
    make
    ./stagbldemo2d
    paraview &

![stagbl2ddemo quickstart](documentation/resources/stagbl2demo_quickstart.png)

## Dependencies

* MPI is required.
* PETSc is currently required (but will not always be).
* Python is required for the configuration process.
