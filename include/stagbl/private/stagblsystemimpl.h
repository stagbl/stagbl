#if !defined(STAGBLSYSTEMIMPL_H_)
#define STAGBLSYSTEMIMPL_H_

#include "stagbl.h"

struct _p_StagBLSystemOps {
  PetscErrorCode (*create)(StagBLSystem);
  PetscErrorCode (*destroy)(StagBLSystem);
};

typedef struct _p_StagBLSystemOps *StagBLSystemOps;

struct _p_StagBLSystem
{
  StagBLSystemOps ops;
  const char    *type;
  void          *data;
  StagBLGrid    grid;
};

#endif
