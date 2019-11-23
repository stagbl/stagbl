#if !defined(STAGBLGRIDIMPL_H_)
#define STAGBLGRIDIMPL_H_

#include "stagbl.h"

struct _p_StagBLGridOps {
  PetscErrorCode (*create)(StagBLGrid);
  PetscErrorCode (*createcompatiblestagblgrid)(StagBLGrid,PetscInt,PetscInt,PetscInt,PetscInt,StagBLGrid*);
  PetscErrorCode (*createstagblarray)(StagBLGrid,StagBLArray*);
  PetscErrorCode (*destroy)(StagBLGrid);
};

typedef struct _p_StagBLGridOps *StagBLGridOps;

struct _p_StagBLGrid
{
  StagBLGridOps ops;
  const char    *type;
  void          *data;
};

#define STAGBLGRIDPETSC "petsc"
PetscErrorCode StagBLGridCreate_PETSc(StagBLGrid);

#endif
