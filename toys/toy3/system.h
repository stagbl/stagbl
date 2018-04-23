#if !defined(SYSTEM_H_)
#define SYSTEM_H_

#include <petsc.h>
#include "ctx.h"

PetscErrorCode CreateSystem(Ctx,Mat*,Vec*);
PetscReal getRho(Ctx ctx,PetscReal,PetscReal);
PetscReal getEta(Ctx ctx,PetscReal,PetscReal);

#endif
