StagBL Tests
------------

StagBL's test suite resides here. This includes

1. Unit and integration tests
2. Scripts for running the test suite

The test suite itself includes the unit and integration tests here as well as runs of StagBLDemo[2d,3d].

The tests will be defined with python scripts leveraging https://bitbucket.org/dmay/pythontestharness
The major rationale for this, over other, more-standard testing systems, is that this framework
is dedicated to providing easy testing of MPI-based code on parallel clusters with batch/queue systems.
