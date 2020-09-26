#if !defined(STAGBLSYSTEMSIMPLEIMPL_H_)
#define STAGBLSYSTEMSIMPLEIMPL_H_

#include "stagbl.h"
#include <petsc.h>

typedef struct {
  PetscScalar *mat;
  StagBLArray rhs;
  PetscInt    system_size;
} StagBLSystem_Simple;

#endif
