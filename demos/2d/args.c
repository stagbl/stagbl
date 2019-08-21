/* Use PETSc tools for convenience is accepting command-line arguments */
#include "args.h"
#include <petscsys.h>

PetscErrorCode GetStringArg(const char *flg,const char *def,size_t len,char *dest)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscStrcpy(dest,def);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,flg,dest,len,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
