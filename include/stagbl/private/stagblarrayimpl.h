#if !defined(STAGBLARRAYIMPL_H_)
#define STAGBLARRAYIMPL_H_

#include "stagbl.h"

struct _p_StagBLArrayOps {
  StagBLErrorCode (*create)(StagBLArray);
  StagBLErrorCode (*destroy)(StagBLArray);
};

typedef struct _p_StagBLArrayOps *StagBLArrayOps;

struct _p_StagBLArray
{
  StagBLArrayOps ops;
  const char    *type;
  void          *data;
};

#endif
