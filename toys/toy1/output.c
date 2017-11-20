#include "output.h"
#include <petscdmda.h>

PetscErrorCode OutputVecBinary(Vec v,const char *filename,PetscBool headless)
{
  PetscErrorCode ierr;
  MPI_Comm       comm;
  PetscViewer    viewer;

  PetscFunctionBeginUser;
  ierr = PetscObjectGetComm((PetscObject)v,&comm);CHKERRQ(ierr);
  ierr = PetscViewerCreate(comm,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERBINARY);CHKERRQ(ierr);
  if (headless) {
    ierr = PetscViewerBinarySetSkipHeader(viewer,PETSC_TRUE);CHKERRQ(ierr);
  }
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr); 
  ierr = VecView(v,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OutputMatBinary(Mat A,const char *filename)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;
  MPI_Comm       comm;

  PetscFunctionBeginUser;
  ierr = PetscObjectGetComm((PetscObject)A,&comm);CHKERRQ(ierr);
  ierr = PetscViewerCreate(comm,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer,PETSCVIEWERBINARY);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer,filename);CHKERRQ(ierr); 
  ierr = MatView(A,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OutputDMCoordsBinary(DM da,const char* filename,PetscBool headless)
{
  PetscErrorCode ierr;
  Vec            coords;

  PetscFunctionBeginUser;
  ierr = DMGetCoordinates(da,&coords);CHKERRQ(ierr);
  ierr = OutputVecBinary(coords,filename,headless);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMDA2dXDMFStart(DM da,const char* filename,const char* coordsBinFilename,PetscViewer *pViewer)
{
  PetscErrorCode ierr;
  MPI_Comm       comm;
  PetscInt       M,N,P;
  PetscViewer    viewer;

  PetscFunctionBeginUser;

  ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(comm,filename,pViewer);CHKERRQ(ierr);
  viewer = *pViewer;
  ierr = PetscViewerASCIIPrintf(viewer,"<?xml version=\"1.0\" ?>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"3.3\">\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"<Domain>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"<Grid GridType=\"Uniform\">\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
  // Topology
  ierr = PetscViewerASCIIPrintf(viewer,"<Topology TopologyType=\"2DSMesh\" Dimensions=\"%D %D\"/>\n",N,M);CHKERRQ(ierr);
  // Geometry
  ierr = PetscViewerASCIIPrintf(viewer,"<Geometry GeometryType=\"XY\">\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"<DataItem Dimensions=\"%D %D 2\"\n",N,M);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer," NumberType=\"Float\" Precision=\"%d\"\n",sizeof(PetscScalar));CHKERRQ(ierr); // 8 for PetscScalar = double
  ierr = PetscViewerASCIIPrintf(viewer," Format=\"Binary\" Endian=\"Big\">\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"   %s\n",coordsBinFilename);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"</DataItem>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);CHKERRQ(ierr);
  PetscViewerASCIIPrintf(viewer,"</Geometry>\n");
  ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMDA2dXDMFAddAttribute(DM da,const char* binFilename,const char* attributeName, PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscInt       M,N,P;
  PetscInt       dof;
  const char     *attributeType;

  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0,&dof,0,0,0,0,0);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
  attributeType = dof == 1 ? "Scalar" : "Vector";
  ierr = PetscViewerASCIIPrintf(viewer,"<Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Node\">\n",attributeName,attributeType);CHKERRQ(ierr); // "Node" centered here, even though these are cells in our simulation. Requires hdf3 reader in Paraview, apparently
  ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
  PetscViewerASCIIPrintf(viewer,"<DataItem Dimensions=\"%D %D %D\" NumberType=\"Float\" Precision=\"%d\" Format=\"Binary\" Endian=\"Big\">\n",N,M,dof,sizeof(PetscScalar));
  PetscViewerASCIIPrintf(viewer,"  %s\n",binFilename);
  PetscViewerASCIIPrintf(viewer,"</DataItem>\n");
  ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  PetscViewerASCIIPrintf(viewer,"</Attribute>\n");
  ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMDA2dXDMFFinish(PetscViewer *pViewer) 
{
  PetscErrorCode ierr;
  PetscViewer viewer;

  PetscFunctionBeginUser;
  viewer = *pViewer;
  ierr = PetscViewerASCIIPrintf(viewer,"</Grid>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"</Domain>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"</Xdmf>\n");CHKERRQ(ierr);
  ierr = PetscViewerDestroy(pViewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
