% Barebones script to load data into MATLAB/OCTAVE for examination

% Update this to point to your PETSc installation
addpath('~/petsc-maint/share/petsc/matlab')

clear A x B
A = PetscBinaryRead('A.pbin');
b = PetscBinaryRead('b.pbin');
x = PescBinaryRead('x.pbin');
