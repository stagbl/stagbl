#include "stagbl.h"
#include <stdio.h>

StagBLErrorCode StagBLInitialize(int argc,char** argv,MPI_Comm comm)
{
  PetscErrorCode ierr;

  if (comm) PETSC_COMM_WORLD = comm;
  ierr = PetscInitialize(&argc,&argv,(char*)0,(void*)0);CHKERRQ(ierr);
  return 0;
}

StagBLErrorCode StagBLFinalize()
{
  PetscErrorCode ierr;

  ierr = PetscFinalize();
  return ierr;
}
