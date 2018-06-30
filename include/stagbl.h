#if !defined(STAGBL_H_)
#define STAGBL_H_
#include <mpi.h>

void StagBLInitialize(int,char**,MPI_Comm);
void StagBLFinalize();

#endif
