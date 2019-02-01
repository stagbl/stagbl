/* Use PETSc tools for convenience is accepting command-line arguments */
#include "args.h"
#include <petscsys.h>

StagBLErrorCode GetIntArg(const char *flg,StagBLInt def,StagBLInt *dest)
{
  PetscErrorCode ierr;

  *dest = def;
  ierr = PetscOptionsGetInt(NULL,NULL,flg,dest,NULL);CHKERRQ(ierr);
  return 0;
}
