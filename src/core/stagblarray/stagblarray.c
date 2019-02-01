#include "stagbl/private/stagblarrayimpl.h"
#include <stdlib.h>

StagBLErrorCode StagBLArrayCreate(StagBLGrid grid, StagBLArray *stagblarray)
{
  *stagblarray = malloc(sizeof(struct _p_StagBLArray));
  (*stagblarray)->ops = calloc(1,sizeof(struct _p_StagBLArrayOps));

  (*stagblarray)->grid = grid;

  // Setting Type and calling creation routine hard-coded for now
  (*stagblarray)->type = STAGBLARRAYPETSC;
  (*stagblarray)->ops->create = StagBLArrayCreate_PETSc; // Sets other ops
  ((*stagblarray)->ops->create)(*stagblarray);

  return 0;
}

StagBLErrorCode StagBLArrayDestroy(StagBLArray *stagblarray)
{
  if ((*stagblarray)->ops->destroy) {
    ((*stagblarray)->ops->destroy)(*stagblarray);
  }
  free((*stagblarray)->ops);
  free(*stagblarray);
  *stagblarray = NULL;

  return 0;
}
