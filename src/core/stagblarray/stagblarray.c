#include "stagblarrayimpl.h"
#include <stdlib.h>

void StagBLArrayCreate(StagBLArray *stagblarray)
{
  *stagblarray = malloc(sizeof(struct _p_StagBLArray));
  (*stagblarray)->ops = calloc(1,sizeof(struct _p_StagBLArrayOps));

  // Setting Type and calling creation routine hard-coded for now
  (*stagblarray)->type = STAGBLARRAYPETSC;
  (*stagblarray)->ops->create = StagBLArrayCreate_PETSc; // Sets other ops
  ((*stagblarray)->ops->create)(*stagblarray);
}

void StagBLArrayDestroy(StagBLArray *stagblarray)
{
  if ((*stagblarray)->ops->destroy) {
    ((*stagblarray)->ops->destroy)(*stagblarray);
  }
  free((*stagblarray)->ops);
  free(*stagblarray);
  *stagblarray = NULL;
}
