Toy Code 1
----------

A toy 2D Stokes code, intended to solve problems 7.1 and 7.2 in Gerya's textbook.

Uses multiple PETSc DMDA objects.

Does not run in parallel, due to lazy/naive implementation of the operator assembly.

7.1 : 
  ./runme -isoviscous
7.2 : 
  ./runme

Output: see paraviewMe.xmf (use xdmf3 Paraview reader)

Additional debugging output includes other .xmf files and PETSc binary data in .petscbin files. See loadData.m to examine in MATLAB or OCTAVE.
