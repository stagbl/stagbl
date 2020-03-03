.. image:: docs/resources/logo/logo_half.png
   :alt: StagBL

.. image:: https://travis-ci.com/stagbl/stagbl.svg?branch=master
    :target: https://travis-ci.com/stagbl/stagbl
    :alt: CI Status

.. image:: https://readthedocs.org/projects/stagbl/badge/?version=latest
    :target: https://stagbl.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

About
=====

StagBL is a C library designed to allow optimized, massively-parallel
Stokes solvers for geodynamics application codes based on finite-volume
methods on regular, orthogonal grids, usually coupled to a
particle-based advection scheme.

It aims to be as lightweight as possible while still providing the
flexibility and extensibility required for scientific application codes.
This accomplished with careful design and interfaces to powerful
external libraries. In particular, its parallel staggered-grid data structure
leverages the `DMStag component <https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DMSTAG/index.html>`__
within `PETSc <https://www.mcs.anl.gov/petsc>`__.

Development of StagBL is supported by the `Platform for Advanced
Scientific Computing <https://www.pasc-ch.org>`__.

Dependencies
============

- `PETSc <https://www.mcs.anl.gov/petsc>`__
-  Python (for configuration)

Documentation
=============
See `StagBL on Read the Docs <https://stagbl.rtfd.io>`__ for additional information.

Quickstart
==========

Step 1/2: PETSc
---------------

Clone a `custom branch <https://bitbucket.org/psanan/petsc/branch/psanan/stagbl-working-base>`__ of PETSc

.. code-block:: bash

    git clone https://bitbucket.org/psanan/petsc -b psanan/stagbl-working-base petsc-stagbl

Configure PETSc with some direct solver packages (SuiteSparse, SuperLU_dist, MUMPS), build, and check. See
`the PETSc website <https://www.mcs.anl.gov/petsc/documentation/installation.html>`__
for full information.

An example configuration command for a local GNU/Linux system is

.. code-block:: bash

    cd petsc-stagbl
    ./configure --download-mpich --with-debugging=0 --download-suitesparse --download-superlu_dist --download-mumps --download-metis --download-parmetis --download-scalapack --FOPTFLAGS='-g -O3' --COPTFLAGS='-g -O3' --CXXOPTFLAGS='-g -O3'
    # Build and check as instructed
    export PETSC_DIR=$PWD # this directory
    export PETSC_ARCH=xxx # replace 'xxx' with the value shown during configuration.
    cd ..

On an OS X system, first `install XCode <https://guide.macports.org/chunked/installing.html#installing.xcode>`__, and check that your compilers work

.. code-block:: bash

    printf '#include<stdio.h>\nint main(){printf("OK!\\n");}' > t.c && /usr/bin/gcc t.c && ./a.out && rm t.c a.out
    printf '#include<iostream>\nint main(){std::cout<<"OK!"<<std::endl;}' > t.cpp && /usr/bin/g++ t.cpp && ./a.out && rm t.cpp a.out

Next, install a Fortran compiler like `gfortran`. If you use MacPorts, Homebrew, Nix, etc. you can use those (recommended),
or you can install it using `these instructions <http://hpc.sourceforge.net>`__. Test your compiler works

.. code-block:: bash

    printf 'program t\nprint"(a)","OK!"\nend program' > t.f90 && gfortran t.f90 && ./a.out && rm t.f90 a.out

Then, configure as

.. code-block:: bash

    cd petsc-stagbl
     ./configure --with-fc=gfortran --with-cc=/usr/bin/gcc --with-cxx=/usr/bin/g++ --download-mpich --download-hdf5 --download-metis --download-parmetis --download-scalapack --download-mumps --download-suitesparse --download-superlu_dist --with-debugging=no --FOPTFLAGS='-g -O3' --COPTFLAGS='-g -O3' --CXXOPTFLAGS='-g -O3' --download-cmake
    # Build and check as instructed
    export PETSC_DIR=$PWD # this directory
    export PETSC_ARCH=xxx # replace 'xxx' with the value shown during configuration.
    cd ..

In either case, note the values of ``PETSC_ARCH`` and ``PETSC_DIR`` defined during configuration.
You may want to define this in a login file, or otherwise record them for next time you log in.

Step 2/2: StagBL
----------------

Clone this repository, including submodules (using git 2.13 or later)

.. code-block:: bash

    git clone --recurse-submodules https://github.com/stagbl/stagbl

Configure and build StagBL, making sure `PETSC_ARCH` and `PETSC_DIR` are defined,
as above. (If you forget these values, `PETSC_DIR` is where you configured PETSc from,
and `PETSC_ARCH` is the name of the directory, e.g. `arch-linux-c-opt` or `arch-darwin-c-opt`,
that was created there during the configuration process).

.. code-block:: bash

    cd stagbl
    ./configure.py         # follow instructions to make
    cd demos
    make 2d
    ./stagbldemo2d
    paraview out_element_0000.vtr &

.. image:: docs/resources/stagbldemo2d_quickstart.png
   :alt: stagbl2ddemo quickstart

In parallel, try

.. code-block:: bash

    $PETSC_DIR/$PETSC_ARCH/bin/mpiexec -np 4 ./stagbldemo2d -mode sinker -stag_grid_x 30 -stag_grid_y 50
    paraview out_element_0000.vtr &

.. image:: docs/resources/stagbldemo2d_quickstart2.png
   :alt: stagbl2ddemo quickstart 2
