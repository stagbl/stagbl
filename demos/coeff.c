#include "coeff.h"

/* Coefficient/forcing Functions */

/* Constant */
static PetscScalar getRho_constant(void *ptr,PetscScalar x, PetscScalar y, PetscScalar z) {
  Ctx ctx = (Ctx) ptr;

  STAGBL_UNUSED(x);
  STAGBL_UNUSED(y);
  STAGBL_UNUSED(z);
  return ctx->rho1;
}

static PetscScalar getEta_constant(void *ptr,PetscScalar x, PetscScalar y, PetscScalar z) {
  Ctx ctx = (Ctx) ptr;

  STAGBL_UNUSED(x);
  STAGBL_UNUSED(y);
  STAGBL_UNUSED(z);
  return ctx->eta1;
}

/* Sinker */
static PetscScalar getRho_sinker2(void *ptr,PetscScalar x, PetscScalar y, PetscScalar z) {
  Ctx ctx = (Ctx) ptr;

  const PetscScalar d = ctx->xmax-ctx->xmin;
  const PetscScalar xx = x/d - 0.5;
  const PetscScalar yy = y/d - 0.5;
  STAGBL_UNUSED(z);
  return (xx*xx + yy*yy) > 0.3*0.3 ? ctx->rho1 : ctx->rho2;
}

static PetscScalar getEta_sinker2(void *ptr,PetscScalar x, PetscScalar y, PetscScalar z) {
  Ctx ctx = (Ctx) ptr;

  const PetscScalar d = ctx->xmax-ctx->xmin;
  const PetscScalar xx = x/d - 0.5;
  const PetscScalar yy = y/d - 0.5;
  STAGBL_UNUSED(z);
  return (xx*xx + yy*yy) > 0.3*0.3 ? ctx->eta1 : ctx->eta2;
}
static PetscScalar getRho_sinker3(void *ptr,PetscScalar x, PetscScalar y, PetscScalar z) {
  Ctx ctx = (Ctx) ptr;

  const PetscScalar d = ctx->xmax-ctx->xmin;
  const PetscScalar xx = x/d - 0.5;
  const PetscScalar yy = y/d - 0.5;
  const PetscScalar zz = z/d - 0.5;
  return (xx*xx + yy*yy + zz*zz) > 0.3*0.3 ? ctx->rho1 : ctx->rho2;
}

static PetscScalar getEta_sinker3(void *ptr,PetscScalar x, PetscScalar y, PetscScalar z) {
  Ctx ctx = (Ctx) ptr;

  const PetscScalar d = ctx->xmax-ctx->xmin;
  const PetscScalar xx = x/d - 0.5;
  const PetscScalar yy = y/d - 0.5;
  const PetscScalar zz = z/d - 0.5;
  return (xx*xx + yy*yy + zz*zz) > 0.3*0.3 ? ctx->eta1 : ctx->eta2;
}

static PetscScalar getRho_sinker_box2(void *ptr,PetscScalar x, PetscScalar y, PetscScalar z) {
  Ctx ctx = (Ctx) ptr;

  STAGBL_UNUSED(z);
  const PetscScalar d = ctx->xmax-ctx->xmin;
  const PetscScalar xx = x/d - 0.5;
  const PetscScalar yy = y/d - 0.5;
  return (xx*xx > 0.3*0.3 || yy*yy > 0.3*0.3) ? ctx->rho1 : ctx->rho2;
}

static PetscScalar getEta_sinker_box2(void *ptr,PetscScalar x, PetscScalar y, PetscScalar z) {
  Ctx ctx = (Ctx) ptr;

  STAGBL_UNUSED(z);
  const PetscScalar d = ctx->xmax-ctx->xmin;
  const PetscScalar xx = x/d - 0.5;
  const PetscScalar yy = y/d - 0.5;
  return (xx*xx > 0.3*0.3 || yy*yy > 0.3*0.3) ? ctx->eta1 : ctx->eta2;
}

static PetscScalar getRho_sinker_box3(void *ptr,PetscScalar x, PetscScalar y, PetscScalar z) {
  Ctx ctx = (Ctx) ptr;

  const PetscScalar d = ctx->xmax-ctx->xmin;
  const PetscScalar xx = x/d - 0.5;
  const PetscScalar yy = y/d - 0.5;
  const PetscScalar zz = z/d - 0.5;
  return (xx*xx > 0.3*0.3 || yy*yy > 0.3*0.3 || zz*zz > 0.3*0.3) ? ctx->rho1 : ctx->rho2;
}

static PetscScalar getEta_sinker_box3(void *ptr,PetscScalar x, PetscScalar y, PetscScalar z) {
  Ctx ctx = (Ctx) ptr;

  const PetscScalar d = ctx->xmax-ctx->xmin;
  const PetscScalar xx = x/d - 0.5;
  const PetscScalar yy = y/d - 0.5;
  const PetscScalar zz = z/d - 0.5;
  return (xx*xx > 0.3*0.3 || yy*yy > 0.3*0.3 || zz*zz > 0.3*0.3) ? ctx->eta1 : ctx->eta2;
}

/* Vertical layers */
static PetscScalar getRho_gerya72(void *ptr,PetscScalar x,PetscScalar y,PetscScalar z)
{
  Ctx ctx = (Ctx) ptr;

  STAGBL_UNUSED(y);
  STAGBL_UNUSED(z);
  if (x < (ctx->xmax-ctx->xmin)/2.0) {
    return ctx->rho1;
  } else {
    return ctx->rho2;
  }
}

static PetscScalar getEta_gerya72(void *ptr,PetscScalar x,PetscScalar y, PetscScalar z)
{
  Ctx ctx = (Ctx) ptr;

  STAGBL_UNUSED(y);
  STAGBL_UNUSED(z);
  if (x < (ctx->xmax-ctx->xmin)/2.0) {
    return ctx->eta1;
  } else {
    return ctx->eta2;
  }
}

PetscErrorCode PopulateCoefficientData(Ctx ctx,const char* mode)
{
  PetscErrorCode  ierr;
  PetscInt        dim;
  PetscInt        N[3];
  PetscInt        ex,ey,ez,startx,starty,startz,nx,ny,nz;
  PetscInt        slot_prev,slot_center;
  PetscInt        slot_rho_downleft,slot_rho_backleft,slot_rho_backdown,slot_eta_element,slot_eta_downleft,slot_eta_backleft,slot_eta_backdown;
  DM              dm_coefficients;
  Vec             *p_coeff_local;
  Vec             coeff_local;
  PetscReal       **arr_coordinates_x,**arr_coordinates_y,**arr_coordinates_z;
  PetscBool       flg;
  StagBLArrayType array_type;

  PetscFunctionBeginUser;

  /* Pull out DM object */
  ierr = StagBLGridPETScGetDM(ctx->coefficient_grid,&dm_coefficients);CHKERRQ(ierr);
  ierr = DMGetDimension(dm_coefficients,&dim);CHKERRQ(ierr);

  /* Set coefficient evaluation functions from mode */
  flg = PETSC_FALSE;
  ierr = PetscStrcmp(mode,"gerya72",&flg);CHKERRQ(ierr);
  if (flg) {
    ctx->getEta = getEta_gerya72;
    ctx->getRho = getRho_gerya72;
  }
  if (!flg) {
    ierr = PetscStrcmp(mode,"sinker",&flg);CHKERRQ(ierr);
    if (flg) {
      switch (dim) {
        case 2:
          ctx->getEta = getEta_sinker2;
          ctx->getRho = getRho_sinker2;
          break;
        case 3:
          ctx->getEta = getEta_sinker3;
          ctx->getRho = getRho_sinker3;
          break;
        default: SETERRQ1(ctx->comm,PETSC_ERR_SUP,"Unsupported dimension %D",dim);
      }
    }
  }
  if (!flg) {
    ierr = PetscStrcmp(mode,"sinker_box",&flg);CHKERRQ(ierr);
    if (flg) {
      switch (dim) {
        case 2:
          ctx->getEta = getEta_sinker_box2;
          ctx->getRho = getRho_sinker_box2;
          break;
        case 3:
          ctx->getEta = getEta_sinker_box3;
          ctx->getRho = getRho_sinker_box3;
          break;
        default: SETERRQ1(ctx->comm,PETSC_ERR_SUP,"Unsupported dimension %D",dim);
      }
    }
  }
  if (!flg) {
    ierr = PetscStrcmp(mode,"blankenbach",&flg);CHKERRQ(ierr);
    if (flg) {
      ctx->getEta = getEta_constant;
      ctx->getRho = getRho_constant;
    }
  }
  if (!flg) {
    SETERRQ1(ctx->comm,PETSC_ERR_ARG_OUTOFRANGE,"Unrecognized mode %s",mode);
  }

  // FIXME: this logic is a kludge which may need to be improved.
  /* If array doesn't exist, create it and pull out a local Vec. Otherwise, get the local Vec */
  if (!ctx->coefficient_array) {
    ierr = StagBLGridCreateStagBLArray(ctx->coefficient_grid,&ctx->coefficient_array);CHKERRQ(ierr);
    ierr = StagBLArrayGetType(ctx->coefficient_array,&array_type);CHKERRQ(ierr);
    if (StagBLCheckType(array_type,STAGBLARRAYPETSC)) {
      ierr = StagBLArrayPETScGetLocalVecPointer(ctx->coefficient_array,&p_coeff_local);CHKERRQ(ierr);
      ierr = DMCreateLocalVector(dm_coefficients,p_coeff_local);CHKERRQ(ierr);
    }
  } else {
    ierr = StagBLArrayGetType(ctx->coefficient_array,&array_type);CHKERRQ(ierr);
  }


  if (StagBLCheckType(array_type,STAGBLARRAYPETSC)) {
    ierr = StagBLArrayPETScGetLocalVec(ctx->coefficient_array,&coeff_local);CHKERRQ(ierr);
  } else if (StagBLCheckType(array_type,STAGBLARRAYSIMPLE)) {
    ierr = DMGetLocalVector(dm_coefficients,&coeff_local);
  } else StagBLError1(PetscObjectComm((PetscObject)dm_coefficients),"Unsupported array type %s",array_type);

  ierr = DMStagGetGhostCorners(dm_coefficients,&startx,&starty,&startz,&nx,&ny,&nz);CHKERRQ(ierr); /* Iterate over all local elements */
  ierr = DMStagGetGlobalSizes(dm_coefficients,&N[0],&N[1],&N[2]);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateArraysRead(dm_coefficients,&arr_coordinates_x,&arr_coordinates_y,&arr_coordinates_z);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateLocationSlot(dm_coefficients,DMSTAG_ELEMENT,&slot_center);CHKERRQ(ierr);
  ierr = DMStagGetProductCoordinateLocationSlot(dm_coefficients,DMSTAG_LEFT,&slot_prev);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dm_coefficients,DMSTAG_ELEMENT,  0,&slot_eta_element);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dm_coefficients,DMSTAG_DOWN_LEFT,0,&slot_eta_downleft);CHKERRQ(ierr);
  ierr = DMStagGetLocationSlot(dm_coefficients,DMSTAG_DOWN_LEFT,1,&slot_rho_downleft);CHKERRQ(ierr);
  if (dim == 2) {
    PetscScalar ***arr_coefficients;

    ierr = DMStagVecGetArray(dm_coefficients,coeff_local,&arr_coefficients);CHKERRQ(ierr);
    for (ey = starty; ey<starty+ny; ++ey) {
      for (ex = startx; ex<startx+nx; ++ex) {
        arr_coefficients[ey][ex][slot_eta_element]  = ctx->getEta(ctx,arr_coordinates_x[ex][slot_center],arr_coordinates_y[ey][slot_center],0.0);
        arr_coefficients[ey][ex][slot_eta_downleft] = ctx->getEta(ctx,arr_coordinates_x[ex][slot_prev],  arr_coordinates_y[ey][slot_prev],  0.0);
        arr_coefficients[ey][ex][slot_rho_downleft] = ctx->getRho(ctx,arr_coordinates_x[ex][slot_prev],  arr_coordinates_y[ey][slot_prev],  0.0);
      }
    }
    ierr = DMStagVecRestoreArray(dm_coefficients,coeff_local,&arr_coefficients);CHKERRQ(ierr);
  } else if (dim == 3) {
    PetscScalar ****arr_coefficients;

    ierr = DMStagGetLocationSlot(dm_coefficients,DMSTAG_BACK_LEFT,0,&slot_eta_backleft);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dm_coefficients,DMSTAG_BACK_LEFT,1,&slot_rho_backleft);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dm_coefficients,DMSTAG_BACK_DOWN,0,&slot_eta_backdown);CHKERRQ(ierr);
    ierr = DMStagGetLocationSlot(dm_coefficients,DMSTAG_BACK_DOWN,1,&slot_rho_backdown);CHKERRQ(ierr);
    ierr = DMStagVecGetArray(dm_coefficients,coeff_local,&arr_coefficients);CHKERRQ(ierr);
    for (ez = startz; ez<startz+nz; ++ez) {
      for (ey = starty; ey<starty+ny; ++ey) {
        for (ex = startx; ex<startx+nx; ++ex) {
          const PetscScalar x_prev = arr_coordinates_x[ex][slot_prev];
          const PetscScalar y_prev = arr_coordinates_y[ey][slot_prev];
          const PetscScalar z_prev = arr_coordinates_z[ez][slot_prev];
          const PetscScalar x_center = arr_coordinates_x[ex][slot_center];
          const PetscScalar y_center = arr_coordinates_y[ey][slot_center];
          const PetscScalar z_center = arr_coordinates_z[ez][slot_center];

          arr_coefficients[ez][ey][ex][slot_eta_element]  = ctx->getEta(ctx,x_center,y_center,z_center);
          arr_coefficients[ez][ey][ex][slot_eta_downleft] = ctx->getEta(ctx,x_prev,  y_prev,  z_center);
          arr_coefficients[ez][ey][ex][slot_rho_downleft] = ctx->getRho(ctx,x_prev,  y_prev,  z_center);
          arr_coefficients[ez][ey][ex][slot_eta_backleft] = ctx->getEta(ctx,x_prev,  y_center,z_prev  );
          arr_coefficients[ez][ey][ex][slot_rho_backleft] = ctx->getRho(ctx,x_prev,  y_center,z_prev  );
          arr_coefficients[ez][ey][ex][slot_eta_backdown] = ctx->getEta(ctx,x_center,y_prev,  z_prev  );
          arr_coefficients[ez][ey][ex][slot_rho_backdown] = ctx->getRho(ctx,x_center,y_prev,  z_prev  );
        }
      }
    }
    ierr = DMStagVecRestoreArray(dm_coefficients,coeff_local,&arr_coefficients);CHKERRQ(ierr);
  } else SETERRQ1(PetscObjectComm((PetscObject)dm_coefficients),PETSC_ERR_SUP,"Unsupported dimension %d",dim);
  ierr = DMStagRestoreProductCoordinateArraysRead(dm_coefficients,&arr_coordinates_x,&arr_coordinates_y,&arr_coordinates_z);CHKERRQ(ierr);

  // FIXME kludge for simple case
  if (StagBLCheckType(array_type,STAGBLARRAYSIMPLE)) {
    PetscScalar       *local_raw;
    const PetscScalar *coeff_local_array;
    PetscInt          n;

    ierr = StagBLArraySimpleGetLocalRaw(ctx->coefficient_array,&local_raw);CHKERRQ(ierr);
    ierr = VecGetLocalSize(coeff_local,&n);CHKERRQ(ierr);
    ierr = VecGetArrayRead(coeff_local,&coeff_local_array);CHKERRQ(ierr);
    for (PetscInt i=0; i<n; ++i) local_raw[i] = coeff_local_array[i];
    ierr = VecRestoreArrayRead(coeff_local,&coeff_local_array);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm_coefficients,&coeff_local);
  }

  /* Declare that, through our custom operations, the local representation is now current and global is not */
  ierr = StagBLArraySetLocalCurrent(ctx->coefficient_array,PETSC_TRUE);CHKERRQ(ierr);
  ierr = StagBLArraySetGlobalCurrent(ctx->coefficient_array,PETSC_FALSE);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
