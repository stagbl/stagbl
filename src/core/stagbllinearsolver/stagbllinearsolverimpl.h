#if !defined(STAGBLLINEARSOLVERIMPL_H_)
#define STAGBLLINEARSOLVERIMPL_H_

#include "stagbl.h"

struct _p_StagBLLinearSolverOps {
  StagBLErrorCode (*create)(StagBLLinearSolver);
  StagBLErrorCode (*destroy)(StagBLLinearSolver);
};

typedef struct _p_StagBLLinearSolverOps *StagBLLinearSolverOps;

struct _p_StagBLLinearSolver
{
  StagBLLinearSolverOps ops;
  const char    *type;
  void          *data;
};

#endif
