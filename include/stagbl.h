#if !defined(STAGBL_H_)
#define STAGBL_H_

#include <petsc.h>

// Data Types (obviously needs to precede all the other declarations here)
#define StagBLErrorCode PetscErrorCode
#define StagBLInt       PetscInt
#define StagBLReal      PetscReal

StagBLErrorCode StagBLInitialize(int,char**,MPI_Comm);
StagBLErrorCode StagBLFinalize();

// Errors
void StagBLErrorFileLine(MPI_Comm,const char*,const char*,long int);
#define StagBLError(comm,msg) SETERRQ1(comm,PETSC_ERR_LIB,"StagBL Error: %s",msg);

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
#define STAGBLGRIDPETSC "petsc"
StagBLErrorCode StagBLGridCreate_PETSc(StagBLGrid); // TODO this and similar functions are private, so need to be moved to a private header
StagBLErrorCode StagBLGridPETScGetDM(StagBLGrid,DM*);
StagBLErrorCode StagBLGridPETScGetDMPointer(StagBLGrid,DM**); // TODO this may not be needed anymore

// StagBLArray Functions
StagBLErrorCode StagBLArrayCreate(StagBLGrid,StagBLArray*); // TODO this needn't be public (or even exist?) because we can only create StagBLArrays from StagBLGrids
StagBLErrorCode StagBLArrayDestroy(StagBLArray*);
StagBLErrorCode StagBLArrayGetStagBLGrid(StagBLArray,StagBLGrid*);

// StagBLArray impls
#define STAGBLARRAYPETSC "petsc"
StagBLErrorCode StagBLArrayCreate_PETSc(StagBLArray);
StagBLErrorCode StagBLArrayPETScGetLocalVec(StagBLArray,Vec*);
StagBLErrorCode StagBLArrayPETScGetGlobalVec(StagBLArray,Vec*);
StagBLErrorCode StagBLArrayPETScGetLocalVecPointer(StagBLArray,Vec**);
StagBLErrorCode StagBLArrayPETScGetGlobalVecPointer(StagBLArray,Vec**);

// StagBLSystem Functions
StagBLErrorCode StagBLSystemCreate(StagBLGrid,StagBLSystem*);
StagBLErrorCode StagBLSystemCreateStagBLSolver(StagBLSystem,StagBLSolver*);
StagBLErrorCode StagBLSystemDestroy(StagBLSystem*);

// StagBLSystem impls
#define STAGBLSYSTEMPETSC "petsc"
StagBLErrorCode StagBLSystemCreate_PETSc(StagBLSystem);
StagBLErrorCode StagBLSystemGetGrid(StagBLSystem,StagBLGrid*);
StagBLErrorCode StagBLSystemPETScGetMat(StagBLSystem,Mat*);
StagBLErrorCode StagBLSystemPETScGetMatPointer(StagBLSystem,Mat**);
StagBLErrorCode StagBLSystemPETScGetVec(StagBLSystem,Vec*);
StagBLErrorCode StagBLSystemPETScGetVecPointer(StagBLSystem,Vec**);

// StagBLSolver Functions
StagBLErrorCode StagBLSolverCreate(StagBLSystem,StagBLSolver*);
StagBLErrorCode StagBLSolverDestroy(StagBLSolver*);
StagBLErrorCode StagBLSolverGetSystem(StagBLSolver,StagBLSystem*);
StagBLErrorCode StagBLSolverSolve(StagBLSolver,StagBLArray);

// StagBLSolver impls
#define STAGBLSOLVERPETSC "petsc"
StagBLErrorCode StagBLSolverCreate_PETSc(StagBLSolver);
StagBLErrorCode StagBLSolverPETScGetKSPPointer(StagBLSolver,KSP**);

// Stokes
StagBLErrorCode StagBLGridCreateStokes2DBox(MPI_Comm,StagBLInt,StagBLInt,StagBLReal,StagBLReal,StagBLReal,StagBLReal,StagBLGrid*);

#endif
