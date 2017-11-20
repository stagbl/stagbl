#ifndef SYSTEM_H__
#define SYSTEM_H__

#include <petscmat.h>
#include "ctx.h"

PetscErrorCode CreateSystem(Mat *pA,Vec *pb, Ctx ctx);
PetscErrorCode ApplyStokes(Mat A,Vec X,Vec Y);

#endif
