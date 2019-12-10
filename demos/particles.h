#ifndef PARTICLES_H_
#define PARTICLES_H_

#include "ctx.h"

PetscErrorCode CreateParticleSystem(Ctx);
PetscErrorCode InterpolateTemperatureToParticles(Ctx);
PetscErrorCode MaterialPoint_AdvectRK1(Ctx,Vec,PetscReal);

#endif
