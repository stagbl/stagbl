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

// StagBLGrid
struct _p_StagBLGrid;
typedef struct _p_StagBLGrid *StagBLGrid;
StagBLErrorCode StagBLGridCreate(StagBLGrid*);
StagBLErrorCode StagBLGridDestroy(StagBLGrid*);

// StagBLGrid impls
#if defined (STAGBL_WITH_PETSC)
#define STAGBLGRIDPETSC "petsc"
StagBLErrorCode StagBLGridCreate_PETSc(StagBLGrid); // TODO this and similar functions are private, so need to be moved to a private header
StagBLErrorCode StagBLGridPETScGetDMPointer(StagBLGrid,DM**);
#endif

// StagBLArray
struct _p_StagBLArray;
typedef struct _p_StagBLArray *StagBLArray;
StagBLErrorCode StagBLArrayCreate(StagBLArray*);
StagBLErrorCode StagBLArrayDestroy(StagBLArray*);

// StagBLArray impls
#if defined (STAGBL_WITH_PETSC)
#define STAGBLARRAYPETSC "petsc"
StagBLErrorCode StagBLArrayCreate_PETSc(StagBLArray);
StagBLErrorCode StagBLArrayPETScGetVecPointer(StagBLArray,Vec**);
#endif

// StagBLOperator
struct _p_StagBLOperator;
typedef struct _p_StagBLOperator *StagBLOperator;
StagBLErrorCode StagBLOperatorCreate(StagBLOperator*);
StagBLErrorCode StagBLOperatorDestroy(StagBLOperator*);

// StagBLOperator impls
#if defined (STAGBL_WITH_PETSC)
#define STAGBLOPERATORPETSC "petsc"
StagBLErrorCode StagBLOperatorCreate_PETSc(StagBLOperator);
StagBLErrorCode StagBLOperatorPETScGetMatPointer(StagBLOperator,Mat**);
#endif

// StagBLLinearSolver
struct _p_StagBLLinearSolver;
typedef struct _p_StagBLLinearSolver *StagBLLinearSolver;
StagBLErrorCode StagBLLinearSolverCreate(StagBLLinearSolver*);
StagBLErrorCode StagBLLinearSolverDestroy(StagBLLinearSolver*);

// StagBLLinearSolver impls
#if defined (STAGBL_WITH_PETSC)
#define STAGBLLINEARSOLVERPETSC "petsc"
StagBLErrorCode StagBLLinearSolverCreate_PETSc(StagBLLinearSolver);
StagBLErrorCode StagBLLinearSolverPETScGetKSPPointer(StagBLLinearSolver,KSP**);
#endif

// Stokes
StagBLErrorCode StagBLGridCreateStokes2DBox(MPI_Comm,StagBLInt,StagBLInt,StagBLReal,StagBLReal,StagBLReal,StagBLReal,StagBLGrid*);

#endif
