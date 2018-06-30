#include "../../stagbloperatorimpl.h"
#include "stagbloperatorpetscimpl.h"
#include <stdlib.h>

void StagBLOperatorDestroy_PETSc(StagBLOperator stagbloperator)
{
  StagBLOperator_PETSc *data = (StagBLOperator_PETSc*) stagbloperator->data;

  if (data->mat) {
    MatDestroy(&data->mat);
  }
  free(stagbloperator->data);
  stagbloperator->data = NULL;
}

void StagBLOperatorCreate_PETSc(StagBLOperator stagbloperator)
{
  StagBLOperator_PETSc *data;

  stagbloperator->data = (void*) malloc(sizeof(StagBLOperator_PETSc));
  data = (StagBLOperator_PETSc*) stagbloperator->data;
  data->mat = NULL;
  stagbloperator->ops->destroy = StagBLOperatorDestroy_PETSc;
}

void StagBLOperatorPETScGetMatPointer(StagBLOperator stagbloperator,Mat **mat)
{
  StagBLOperator_PETSc * const data = (StagBLOperator_PETSc*) stagbloperator->data;
  *mat = &(data->mat);
}
