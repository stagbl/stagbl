#if !defined(STAGBLSOLVERPETSCIMPL_H_)
#define STAGBLSOLVERPETSCIMPL_H_

#include "stagbl.h"
#include <petsc.h>

typedef struct {
  SNES snes;
} StagBLSolver_PETSc;

#endif
