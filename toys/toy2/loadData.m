% Barebones script to load data into MATLAB/OCTAVE for examination

% Update this to point to your PETSc installation
addpath('/Users/patrick/petsc-maint/share/petsc/matlab')

clear A x B
A = PetscBinaryRead('A.petscbin');
b = PetscBinaryRead('b.petscbin');
x = PetscBinaryRead('x.petscbin');
