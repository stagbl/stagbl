#include "stagbl.h"
#include <string.h>
#include <stdio.h>

PetscBool StagBLCheckType(const char *type1, const char *type2)
{
  return (PetscBool) strcmp(type1, type2) == 0;
}

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
