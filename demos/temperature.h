#ifndef TEMPERATURE_H_
#define TEMPERATURE_H_
#include "ctx.h"
PetscErrorCode InitializeTemperature(Ctx);
PetscErrorCode PopulateTemperatureSystem(Ctx);
PetscErrorCode UpdateTemperature(Ctx);
#endif
