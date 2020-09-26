#if !defined(STAGBLSYSTEMIMPL_H_)
#define STAGBLSYSTEMIMPL_H_

#include "stagbl.h"

struct data_StagBLSystemOps {
  PetscErrorCode (*create)(StagBLSystem);
  PetscErrorCode (*destroy)(StagBLSystem);
  PetscErrorCode (*operatorsetvaluesstencil)(StagBLSystem,PetscInt,const DMStagStencil*,PetscInt,const DMStagStencil*,const PetscScalar*);
  PetscErrorCode (*rhssetconstant)(StagBLSystem,PetscScalar);
  PetscErrorCode (*rhssetvaluesstencil)(StagBLSystem,PetscInt,const DMStagStencil*,const PetscScalar*);
};

typedef struct data_StagBLSystemOps *StagBLSystemOps;

struct data_StagBLSystem
{
  StagBLSystemOps  ops;
  const char       *type;
  void             *data;
  StagBLGrid       grid;
  StagBLSolverType solver_type;
};

#endif
