#include "nusselt.h"

/**
  * Compute the Nusselt number, as described in BlankenBachEtAl1989, 2.2 (i)
  * We approximate the temperature gradient at the surface with a lowest-order
  * finite difference.
  * We assume a constant bottom temperature.
  *
  */
PetscErrorCode ComputeNusseltNumber(Ctx ctx,PetscReal *nusselt_number)
{
  PetscErrorCode ierr;
  DM             dm_temperature;
  Vec            vec_temperature,vec_temperature_local;
  PetscInt       startx,starty,nx,ny,nextrax,ex;
  PetscBool      is_top_rank;
  PetscScalar    ***arr_temperature;
  PetscScalar    **arr_coordinates_x,**arr_coordinates_y;
  PetscInt       slot_vertex,slot_prev,slot_next;
  PetscReal      sum_local,sum_global;

  PetscFunctionBeginUser;

  /* Get access to PETSc DM and Vec */
  ierr = StagBLGridPETScGetDM(ctx->temperature_grid,&dm_temperature);CHKERRQ(ierr);
  ierr = StagBLArrayPETScGetGlobalVec(ctx->temperature_array,&vec_temperature);CHKERRQ(ierr);

  /* Iterate over the top row of elements, if we're on the top */
  ierr = DMStagGetCorners(dm_temperature,&startx,&starty,NULL,&nx,&ny,NULL,&nextrax,NULL,NULL);CHKERRQ(ierr);
  ierr = DMStagGetIsLastRank(dm_temperature,NULL,&is_top_rank,NULL);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm_temperature,&vec_temperature_local);CHKERRQ(ierr);
  ierr = DMGlobalToLocal(dm_temperature,vec_temperature,INSERT_VALUES,vec_temperature_local);CHKERRQ(ierr);
  ierr = DMStagVecGetArrayRead(dm_temperature,vec_temperature_local,&arr_temperature);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dm_temperature,DMSTAG_UP_LEFT,0,&slot_vertex);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateLocationSlot(dm_temperature,DMSTAG_LEFT,&slot_prev);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateLocationSlot(dm_temperature,DMSTAG_RIGHT,&slot_next);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateArraysRead(dm_temperature,&arr_coordinates_x,&arr_coordinates_y,NULL);CHKERRQ(ierr);
  sum_local = 0.0;
  if (is_top_rank) {
    const PetscInt ey = starty+ny-1; /* The top element */
    for (ex=0; ex<startx+nx+nextrax; ++ex) { /* Include "extra" vertex on the right */
      PetscScalar hx = arr_coordinates_x[ex][slot_next]-arr_coordinates_x[ex][slot_prev];
      PetscScalar hy = arr_coordinates_y[ey][slot_next]-arr_coordinates_y[ey][slot_prev];
      const PetscScalar dT = arr_temperature[ey][ex][slot_vertex] - arr_temperature[ey-1][ex][slot_vertex];
      sum_local += hx * dT / hy;
    }
  }
  ierr = DMRestoreLocalVector(dm_temperature,&vec_temperature_local);CHKERRQ(ierr);

  /* Perform a reduction */
  ierr = MPI_Allreduce(&sum_local,&sum_global,1,MPIU_SCALAR,MPI_SUM,PetscObjectComm((PetscObject)dm_temperature));CHKERRQ(ierr);

  /* Compute the result (on all ranks, but one could just compute it on rank 0) */
  *nusselt_number = -1.0 * ((ctx->ymax-ctx->ymin) / ((ctx->xmax - ctx->xmin) * ctx->temperature_bottom)) * sum_global;

  PetscFunctionReturn(0);
}
