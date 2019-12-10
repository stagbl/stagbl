#include "stagbl.h"
#include <stdio.h>

PetscErrorCode StagBLInitialize(int argc,char** argv,const char *help,MPI_Comm comm)
{
  PetscErrorCode ierr;

  if (comm && comm != MPI_COMM_NULL) PETSC_COMM_WORLD = comm;
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr);
  return 0;
}

PetscErrorCode StagBLFinalize()
{
  PetscErrorCode ierr;

  ierr = PetscFinalize();
  return ierr;
}
