#ifndef SYSTEM_H_
#define SYSTEM_H_
#include "ctx.h"
#include <petscdmstag.h>
PetscErrorCode CreateSystem(const Ctx,Mat*,Vec*);
#endif
