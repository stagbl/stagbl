#if !defined(STAGBLSYSTEMPETSCIMPL_H_)
#define STAGBLSYSTEMPETSCIMPL_H_

#include "stagbl.h"
#include <petsc.h>

typedef struct {
  Mat  mat;
  Vec  rhs;
  SNES snes;
  PetscErrorCode (*residual_function) (SNES,Vec,Vec,void*);
  PetscErrorCode (*jacobian_function) (SNES,Vec,Mat,Mat,void*);
} StagBLSystem_PETSc;

#endif
