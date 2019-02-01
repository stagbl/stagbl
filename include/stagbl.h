#if !defined(STAGBL_H_)
#define STAGBL_H_

// Hard-coding for now (later put in configure-generated header)
#define STAGBL_WITH_PETSC

#if defined(STAGBL_WITH_PETSC)
#include <petsc.h>
#endif
#include <mpi.h>

// Data Types (obviously needs to precede all the other declarations here)
#if defined(STAGBL_WITH_PETSC)
#define StagBLErrorCode PetscErrorCode
#define StagBLInt       PetscInt
#define StagBLReal      PetscReal
#else
#define StagBLErrorCode int
#define StagBLInt       int
#define StagBLReal      double
#endif

StagBLErrorCode StagBLInitialize(int,char**,MPI_Comm);
StagBLErrorCode StagBLFinalize();

// Errors
void StagBLErrorFileLine(MPI_Comm,const char*,const char*,long int);
#if defined(STAGBL_WITH_PETSC)
#define StagBLError(comm,msg) SETERRQ1(comm,PETSC_ERR_LIB,"StagBL Error: %s",msg);
#else
#define StagBLError(comm,msg) StagBLErrorFileLine(comm,msg,__FILE__,__LINE__)
// Note: this check could be disabled in opt mode (though almost certainly not a real performance concern, one comparison with zero)
#define CHKERRQ(ierr) do {if(ierr) {StagBLError(MPI_COMM_SELF,"StagBL Error Check failed");}} while(0)
#endif

// StagBLGrid Data
struct _p_StagBLGrid;
typedef struct _p_StagBLGrid *StagBLGrid;

// StagBLArray Data
struct _p_StagBLArray;
typedef struct _p_StagBLArray *StagBLArray;

// StagBLSystem Data
struct _p_StagBLSystem;
typedef struct _p_StagBLSystem *StagBLSystem;

// StagBLSolver Data
struct _p_StagBLSolver;
typedef struct _p_StagBLSolver *StagBLSolver;

// StagBLGrid Functions
StagBLErrorCode StagBLGridCreate(StagBLGrid*);
StagBLErrorCode StagBLGridCreateCompatibleStagBLGrid(StagBLGrid,StagBLInt,StagBLInt,StagBLInt,StagBLInt,StagBLGrid*);
StagBLErrorCode StagBLGridCreateStagBLArray(StagBLGrid,StagBLArray*);
StagBLErrorCode StagBLGridCreateStagBLSystem(StagBLGrid,StagBLSystem*);
StagBLErrorCode StagBLGridDestroy(StagBLGrid*);

// StagBLGrid impls
#if defined (STAGBL_WITH_PETSC)
#define STAGBLGRIDPETSC "petsc"
StagBLErrorCode StagBLGridCreate_PETSc(StagBLGrid); // TODO this and similar functions are private, so need to be moved to a private header
StagBLErrorCode StagBLGridPETScGetDM(StagBLGrid,DM*);
StagBLErrorCode StagBLGridPETScGetDMPointer(StagBLGrid,DM**); // TODO this may not be needed anymore
#endif

// StagBLArray Functions
StagBLErrorCode StagBLArrayCreate(StagBLGrid,StagBLArray*); // TODO this needn't be public (or even exist?) because we can only create StagBLArrays from StagBLGrids
StagBLErrorCode StagBLArrayDestroy(StagBLArray*);

// StagBLArray impls
#if defined (STAGBL_WITH_PETSC)
#define STAGBLARRAYPETSC "petsc"
StagBLErrorCode StagBLArrayCreate_PETSc(StagBLArray);
StagBLErrorCode StagBLArrayPETScGetLocalVec(StagBLArray,Vec*);
StagBLErrorCode StagBLArrayPETScGetGlobalVec(StagBLArray,Vec*);
StagBLErrorCode StagBLArrayPETScGetLocalVecPointer(StagBLArray,Vec**);
StagBLErrorCode StagBLArrayPETScGetGlobalVecPointer(StagBLArray,Vec**);
#endif

// StagBLSystem Functions
StagBLErrorCode StagBLSystemCreate(StagBLSystem*);
StagBLErrorCode StagBLSystemCreateStagBLSolver(StagBLSystem,StagBLSolver*);
StagBLErrorCode StagBLSystemDestroy(StagBLSystem*);

// StagBLSystem impls
#if defined (STAGBL_WITH_PETSC)
#define STAGBLSYSTEMPETSC "petsc"
StagBLErrorCode StagBLSystemCreate_PETSc(StagBLSystem);
StagBLErrorCode StagBLSystemPETScGetMatPointer(StagBLSystem,Mat**);
StagBLErrorCode StagBLSystemPETScGetVecPointer(StagBLSystem,Vec**);
#endif

// StagBLSolver Functions
StagBLErrorCode StagBLSolverCreate(StagBLSolver*);
StagBLErrorCode StagBLSolverDestroy(StagBLSolver*);

// StagBLSolver impls
#if defined (STAGBL_WITH_PETSC)
#define STAGBLSOLVERPETSC "petsc"
StagBLErrorCode StagBLSolverCreate_PETSc(StagBLSolver);
StagBLErrorCode StagBLSolverPETScGetKSPPointer(StagBLSolver,KSP**);
#endif

// Stokes
StagBLErrorCode StagBLGridCreateStokes2DBox(MPI_Comm,StagBLInt,StagBLInt,StagBLReal,StagBLReal,StagBLReal,StagBLReal,StagBLGrid*);

#endif
