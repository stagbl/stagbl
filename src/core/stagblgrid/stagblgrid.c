#include "stagblgridimpl.h"
#include <stdlib.h>

void StagBLGridCreate(StagBLGrid *stagblgrid)
{
  *stagblgrid = malloc(sizeof(struct _p_StagBLGrid));
  (*stagblgrid)->ops = calloc(1,sizeof(struct _p_StagBLGridOps));

  // Setting Type and calling creation routine hard-coded for now
  (*stagblgrid)->type = STAGBLGRIDPETSC;
  (*stagblgrid)->ops->create = StagBLGridCreate_PETSc; // Sets other ops
  ((*stagblgrid)->ops->create)(*stagblgrid);
}

void StagBLGridDestroy(StagBLGrid *stagblgrid)
{
  if ((*stagblgrid)->ops->destroy) {
    ((*stagblgrid)->ops->destroy)(*stagblgrid);
  }
  free((*stagblgrid)->ops);
  free(*stagblgrid);
  *stagblgrid = NULL;
}
