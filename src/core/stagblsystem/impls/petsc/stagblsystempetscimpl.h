#if !defined(STAGBLSYSTEMPETSCIMPL_H_)
#define STAGBLSYSTEMPETSCIMPL_H_

#include "stagbl.h"
#include <petsc.h>

typedef struct {
  Mat mat;
  Vec rhs;
} StagBLSystem_PETSc;

#endif
