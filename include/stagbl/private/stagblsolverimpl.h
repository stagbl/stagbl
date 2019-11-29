#if !defined(STAGBLSOLVERIMPL_H_)
#define STAGBLSOLVERIMPL_H_

#include "stagbl.h"

struct data_StagBLSolverOps {
  PetscErrorCode (*create)(StagBLSolver);
  PetscErrorCode (*destroy)(StagBLSolver);
  PetscErrorCode (*solve)(StagBLSolver,StagBLArray);
};

typedef struct data_StagBLSolverOps *StagBLSolverOps;

struct data_StagBLSolver
{
  StagBLSolverOps ops;
  const char    *type;
  void          *data;
  StagBLSystem  system;
};

#endif
