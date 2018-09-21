StagBL Tests
------------

StagBL's tests are organized here. This directory includes

1. Unit and integration tests
2. Scripts for running the test suite

The test suite relies on executables defined in main library source as well as on StagBLDemo[2d,3d].

## Design

The tests will be defined with python scripts leveraging https://bitbucket.org/dmay/pythontestharness
The major rationale for this, over other, more-standard testing systems, is that this framework
is dedicated to providing easy testing of MPI-based code on parallel clusters with batch/queue systems.

Inspired with our experience with the pTatin3D tests, we allow the directory structure here
to define the test names. As such, moving or renaming a test should be as simple as changing
directory names.

## Misc.

### DMStag Debugging Routines
Note that some useful checks to confirm consistency of local->global and global->local maps in DMStag were removed,
but are archived on this branch:
https://bitbucket.org/psanan/petsc-private/branch/psanan/dmstag-l2g-debugging
