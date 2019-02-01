#if !defined(STAGBLOPERATORPETSCIMPL_H_)
#define STAGBLOPERATORPETSCIMPL_H_

#include "stagbl.h"
#include <petsc.h>

typedef struct {
  Mat mat;
} StagBLOperator_PETSc;

#endif
