#include "mpi.h"
#include "stagbl.h"

/* Note that this is usually called from a macro StagBLError(), which invokes
   the PETSc error handler if available */
void StagBLErrorFileLine(MPI_Comm comm, const char* msg,const char* file, long int line)
{
  int rank;
  MPI_Comm_rank(comm,&rank);
  if (rank == 0) {
    fprintf(stderr,"StagBL Error at %s:%ld: %s\n",file,line,msg);  
  }
  fflush(stderr);
  MPI_Barrier(comm);
  MPI_Abort(comm,MPI_ERR_UNKNOWN);
}
