% Barebones script to load data into MATLAB/OCTAVE for examination

% Update this to point to your PETSc installation
addpath('/Users/patrick/petsc-maint/share/petsc/matlab')

clear A x B
A = pbinaryRead('A.pbin');
b = pbinaryRead('b.pbin');
x = pbinaryRead('x.pbin');
