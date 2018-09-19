#include "stagbl.h"
#include <stdio.h>

void StagBLInitialize(int argc,char** argv,MPI_Comm comm)
{
  printf("Hello from StagBL\n");
#if defined(STAGBL_WITH_PETSC)
  {
    if (comm) PETSC_COMM_WORLD = comm;
    PetscInitialize(&argc,&argv,(char*)0,(void*)0);
    PetscPrintf(PETSC_COMM_WORLD,"PETSc active\n");
  }

#endif
}

void StagBLFinalize()
{
#if defined(STAGBL_WITH_PETSC)
  PetscFinalize();
#endif
}
