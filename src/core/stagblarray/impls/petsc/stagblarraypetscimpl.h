#if !defined(STAGBLARRAYPETSCIMPL_H_)
#define STAGBLARRAYPETSCIMPL_H_

#include "stagbl.h"
#include <petsc.h>

typedef struct {
  Vec local, global;
} StagBLArray_PETSc;

#endif
