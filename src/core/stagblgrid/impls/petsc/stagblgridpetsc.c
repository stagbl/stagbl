#include "../../stagblgridimpl.h"
#include "stagblgridpetscimpl.h"
#include <stdlib.h>

StagBLErrorCode StagBLGridDestroy_PETSc(StagBLGrid stagblgrid)
{
  StagBLGrid_PETSc *data = (StagBLGrid_PETSc*) stagblgrid->data;
  if (data->dm) {
    DMDestroy(&data->dm);
  }
  free(stagblgrid->data);
  stagblgrid->data = NULL;

  return 0;
}

StagBLErrorCode StagBLGridCreate_PETSc(StagBLGrid stagblgrid)
{
  StagBLGrid_PETSc *data;
  stagblgrid->data = (void*) malloc(sizeof(StagBLGrid_PETSc));
  data = (StagBLGrid_PETSc*) stagblgrid->data;
  data->dm = NULL;
  stagblgrid->ops->destroy = StagBLGridDestroy_PETSc;

  return 0;
}

StagBLErrorCode StagBLGridPETScGetDMPointer(StagBLGrid stagblgrid,DM **dm)
{
  StagBLGrid_PETSc * const data = (StagBLGrid_PETSc*) stagblgrid->data;
  *dm = &(data->dm);

  return 0;
}
