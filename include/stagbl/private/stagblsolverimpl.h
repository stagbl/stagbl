#if !defined(STAGBLSOLVERIMPL_H_)
#define STAGBLSOLVERIMPL_H_

#include "stagbl.h"

struct _p_StagBLSolverOps {
  StagBLErrorCode (*create)(StagBLSolver);
  StagBLErrorCode (*destroy)(StagBLSolver);
};

typedef struct _p_StagBLSolverOps *StagBLSolverOps;

struct _p_StagBLSolver
{
  StagBLSolverOps ops;
  const char    *type;
  void          *data;
};

#endif
