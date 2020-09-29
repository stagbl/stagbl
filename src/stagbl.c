#include "stagbl.h"
#include <string.h>
#include <stdio.h>

PetscBool StagBLCheckType(const char *type1, const char *type2)
{
  return (PetscBool) strcmp(type1, type2) == 0;
}

PetscErrorCode StagBLEnforceType(const char *type, const char *type_required)
{
  // FIXME - after rebase to PETSc master, use PetscDefined(USE_DEBUG)
#if defined(PETSC_USE_DEBUG)
  PetscFunctionBegin;
  if (!StagBLCheckType(type,type_required)) StagBLError2(PETSC_COMM_WORLD,"Type %s is required, but %s used",type_required,type);
  PetscFunctionReturn(0);
#else
  STAGBL_UNUSED(type);
  STAGBL_UNUSED(type_required);
  return 0;
#endif
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
