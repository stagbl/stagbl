StagBL Documentation
--------------------

As much as possible, documentation is by example.

The demos in `../demos` are central in this, demonstrating functionality in
terms of reproducing standard examples and benchmarks.

The DMStag tutorials in `$PETSC_DIR/src/dm/impls/examples/tutorials` may also
be useful.

Some unit/integration tests in `../tests` might also be of use to some users.

Any additional documentation is generated using Sphinx. From this directory,

  make html

And see `_build/html/index.html`.

Prefer to document functions, macros, etc. with javadoc-style comment blocks,
e.g. `/** ... */`, as these are recognized by Hawkmoth and Doxygen.
