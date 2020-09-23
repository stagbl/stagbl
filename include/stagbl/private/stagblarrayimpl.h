#if !defined(STAGBLARRAYIMPL_H_)
#define STAGBLARRAYIMPL_H_

#include "stagbl.h"

struct data_StagBLArrayOps {
  PetscErrorCode (*create)(StagBLArray);
  PetscErrorCode (*destroy)(StagBLArray);
  PetscErrorCode (*print)(StagBLArray);
  PetscErrorCode (*setlocalconstant)(StagBLArray,PetscScalar);
};

typedef struct data_StagBLArrayOps *StagBLArrayOps;

struct data_StagBLArray
{
  StagBLArrayOps   ops;
  StagBLGrid       grid;
  StagBLArrayType  type;
  PetscBool        current_local, current_global;
  void             *data;
};

#endif
