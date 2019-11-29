#if !defined(STAGBLSYSTEMIMPL_H_)
#define STAGBLSYSTEMIMPL_H_

#include "stagbl.h"

struct data_StagBLSystemOps {
  PetscErrorCode (*create)(StagBLSystem);
  PetscErrorCode (*destroy)(StagBLSystem);
};

typedef struct data_StagBLSystemOps *StagBLSystemOps;

struct data_StagBLSystem
{
  StagBLSystemOps ops;
  const char    *type;
  void          *data;
  StagBLGrid    grid;
};

#endif
