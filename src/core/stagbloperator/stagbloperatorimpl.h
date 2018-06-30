#if !defined(STAGBLOPERATORIMPL_H_)
#define STAGBLOPERATORIMPL_H_

#include "stagbl.h"

struct _p_StagBLOperatorOps {
  void (*create)(StagBLOperator);
  void (*destroy)(StagBLOperator);
};

typedef struct _p_StagBLOperatorOps *StagBLOperatorOps;

struct _p_StagBLOperator
{
  StagBLOperatorOps ops;
  const char    *type;
  void          *data;
};

#endif
