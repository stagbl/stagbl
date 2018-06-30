#include "stagbloperatorimpl.h"
#include <stdlib.h>

void StagBLOperatorCreate(StagBLOperator *stagbloperator)
{
  *stagbloperator = malloc(sizeof(struct _p_StagBLOperator));
  (*stagbloperator)->ops = calloc(1,sizeof(struct _p_StagBLOperatorOps));

  // Setting Type and calling creation routine hard-coded for now
  (*stagbloperator)->type = STAGBLOPERATORPETSC;
  (*stagbloperator)->ops->create = StagBLOperatorCreate_PETSc; // Sets other ops
  ((*stagbloperator)->ops->create)(*stagbloperator);
}

void StagBLOperatorDestroy(StagBLOperator *stagbloperator)
{
  if ((*stagbloperator)->ops->destroy) {
    ((*stagbloperator)->ops->destroy)(*stagbloperator);
  }
  free((*stagbloperator)->ops);
  free(*stagbloperator);
  *stagbloperator = NULL;
}
