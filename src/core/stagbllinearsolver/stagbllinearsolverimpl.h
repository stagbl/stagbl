#if !defined(STAGBLLINEARSOLVERIMPL_H_)
#define STAGBLLINEARSOLVERIMPL_H_

#include "stagbl.h"

struct _p_StagBLLinearSolverOps {
  void (*create)(StagBLLinearSolver);
  void (*destroy)(StagBLLinearSolver);
};

typedef struct _p_StagBLLinearSolverOps *StagBLLinearSolverOps;

struct _p_StagBLLinearSolver
{
  StagBLLinearSolverOps ops;
  const char    *type;
  void          *data;
};

#endif
