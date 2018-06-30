#if !defined(STAGBLLINEARSOLVERPETSCIMPL_H_)
#define STAGBLLINEARSOLVERPETSCIMPL_H_

#include "stagbl.h"
#include <petsc.h>

typedef struct {
  KSP ksp;
} StagBLLinearSolver_PETSc;

#endif
