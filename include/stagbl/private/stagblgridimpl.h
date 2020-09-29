#if !defined(STAGBLGRIDIMPL_H_)
#define STAGBLGRIDIMPL_H_

#include "stagbl.h"

struct data_StagBLGridOps {
  PetscErrorCode (*create)(StagBLGrid);
  PetscErrorCode (*createcompatiblestagblgrid)(StagBLGrid,PetscInt,PetscInt,PetscInt,PetscInt,StagBLGrid*);
  PetscErrorCode (*createstagblarray)(StagBLGrid,StagBLArray*);
  PetscErrorCode (*destroy)(StagBLGrid);
};

typedef struct data_StagBLGridOps *StagBLGridOps;

struct data_StagBLGrid
{
  StagBLGridOps ops;
  const char       *type;
  StagBLArrayType  array_type;
  StagBLSystemType system_type;
  void             *data;
};

#define STAGBLGRIDPETSC "petsc"
PetscErrorCode StagBLGridCreate_PETSc(StagBLGrid);

#endif
