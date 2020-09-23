#if !defined(STAGBLARRAYSIMPLEIMPL_H_)
#define STAGBLARRAYSIMPLEIMPL_H_

#include "stagbl.h"

typedef struct {
  PetscScalar *local, *global;
} StagBLArray_Simple;

#endif
