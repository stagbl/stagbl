#include "stagbl.h"
#include <string.h>

PetscBool StagBLCheckType(const char *type1, const char *type2)
{
  return (PetscBool) strcmp(type1, type2) == 0;
}
