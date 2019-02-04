#if !defined(STAGBLSOLVERPETSCIMPL_H_)
#define STAGBLSOLVERPETSCIMPL_H_

#include "stagbl.h"
#include <petsc.h>

typedef struct {
  KSP ksp;
} StagBLSolver_PETSc;

#endif
