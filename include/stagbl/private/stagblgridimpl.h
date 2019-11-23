#if !defined(STAGBLGRIDIMPL_H_)
#define STAGBLGRIDIMPL_H_

#include "stagbl.h"

struct _p_StagBLGridOps {
  StagBLErrorCode (*create)(StagBLGrid);
  StagBLErrorCode (*createcompatiblestagblgrid)(StagBLGrid,StagBLInt,StagBLInt,StagBLInt,StagBLInt,StagBLGrid*);
  StagBLErrorCode (*createstagblarray)(StagBLGrid,StagBLArray*);
  StagBLErrorCode (*destroy)(StagBLGrid);
};

typedef struct _p_StagBLGridOps *StagBLGridOps;

struct _p_StagBLGrid
{
  StagBLGridOps ops;
  const char    *type;
  void          *data;
};

#define STAGBLGRIDPETSC "petsc"
StagBLErrorCode StagBLGridCreate_PETSc(StagBLGrid);

#endif
