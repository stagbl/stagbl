#ifndef DUMP_H_
#define DUMP_H_

#include "ctx.h"

PetscErrorCode DumpParticles(Ctx,PetscInt);
PetscErrorCode DumpStokes(Ctx,PetscInt);
PetscErrorCode DumpTemperature(Ctx,PetscInt);

#endif
