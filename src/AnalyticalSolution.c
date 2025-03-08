/**
 * @file AnalyticalSolution.c  //  Particle In Cell main.
 * @brief Provides the methods required to generate analytical solutions that can be used to test the particle tracking mechanism.
 * 
**/

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscdmcomposite.h>

#include "AnalyticalSolution.h"

/**
 * @brief Sets the local Cartesian vector field based on input coordinates.
 *
 * This function computes the vector components by applying the sine function to the
 * input coordinate values along the x, y, and z directions.
 *
 * @param[in]     fieldName Pointer to a string representing the field name (for logging purposes).
 * @param[in,out] vecField  Pointer to the `Cmpnts` structure where the computed vector field will be stored.
 * @param[in]     coor      Pointer to the `Cmpnts` structure containing the input coordinate values.
 *
 * @return PetscErrorCode Returns 0 on successful execution, non-zero on failure.
 */
PetscErrorCode SetLocalCartesianField_Vector(const char *fieldName, Cmpnts *vecField, Cmpnts *coor) {
    // Log input coordinate values for debugging purposes.
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
        "SetLocalCartesianField_Vector: Input Coordinates - x: %f, y: %f, z: %f\n",
        coor->x, coor->y, coor->z);

    // Compute vector components as the sine of the coordinate values.
    vecField->x = sin(coor->x);
    vecField->y = sin(coor->y);
    vecField->z = sin(coor->z);

    // Log computed vector field values along with the field name.
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
        "SetLocalCartesianField_Vector: Computed Vector %s - x: %f, y: %f, z: %f\n",
        fieldName, vecField->x, vecField->y, vecField->z);

    return 0;
}

/**
 * @brief Sets the local Cartesian scalar field based on input coordinates.
 *
 * This function computes the scalar field value by combining the sine of the input
 * coordinate values. In this example, the scalar field is computed as the sum of the
 * sine functions of the x, y, and z coordinates.
 *
 * @param[in]     fieldName   Pointer to a string representing the field name (for logging purposes).
 * @param[in,out] scalarField Pointer to the PetscReal where the computed scalar field value will be stored.
 * @param[in]     coor        Pointer to the `Cmpnts` structure containing the input coordinate values.
 *
 * @return PetscErrorCode Returns 0 on successful execution, non-zero on failure.
 */
PetscErrorCode SetLocalCartesianField_Scalar(const char *fieldName, PetscReal *scalarField, Cmpnts *coor) {
    // Log input coordinate values for debugging purposes.
    LOG_ALLOW(LOCAL, LOG_DEBUG,
        "SetLocalCartesianField_Scalar: Input Coordinates - x: %f, y: %f, z: %f\n",
        coor->x, coor->y, coor->z);

    // Compute scalar field value as the sum of the sine functions of the coordinates.
    *scalarField = sin(coor->x) + sin(coor->y) + sin(coor->z);

    // Log computed scalar field value along with the field name.
    LOG_ALLOW(LOCAL, LOG_DEBUG,
        "SetLocalCartesianField_Scalar: Computed Scalar %s - value: %f\n",
        fieldName, *scalarField);

    return 0;
}

/**
 * @brief Sets an analytical Cartesian field (scalar or vector) for cell centers based on a field name.
 *
 * This function looks up the field within the user context by comparing the provided field name
 * against the supported names:
 *   - Vector fields: "Ucat" and "Ucont"
 *   - Scalar fields: "P" and "nvert"
 *
 * If the field is found, the function verifies that the DM block size matches the expected value
 * (3 for vector fields and 1 for scalar fields), retrieves local DM information for both the coordinate DM (da)
 * and the cell-centered DM (fda), interpolates the corner-based coordinates (from da) to cell centers,
 * and then updates the field using the generic helper macro SetLocalCartesianField. Interior cells are
 * updated with the interpolated coordinates, and boundary cells are updated using the original coordinate data.
 *
 * If the field name is not found in the user context, the function throws an error.
 *
 * @param[in]  user      Pointer to the UserCtx structure containing:
 *                         - da: DM for coordinate (corner) data.
 *                         - fda: DM for cell-centered data.
 *                         - Ucat: vector field (Cmpnts ***)
 *                         - Ucont: vector field (Cmpnts ***)
 *                         - P: scalar field (PetscReal ***)
 *                         - nvert: scalar field (PetscReal ***)
 * @param[in]  fieldName Name of the field to update.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode SetAnalyticalCartesianField(UserCtx *user, const char *fieldName)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,
       "SetAnalyticalCartesianField: Starting for field '%s' on rank %d.\n", fieldName, rank);

  /* Look up the field Vec in user context. */
  Vec fieldVec = NULL;
  int fieldIsVector = -1;
  if (strcmp(fieldName, "Ucat") == 0) {
    fieldVec = user->Ucat;
    fieldIsVector = 1;
  } else if (strcmp(fieldName, "Ucont") == 0) {
    fieldVec = user->Ucont;
    fieldIsVector = 1;
  } else if (strcmp(fieldName, "P") == 0) {
    fieldVec = user->P;
    fieldIsVector = 0;
  } else if (strcmp(fieldName, "Nvert") == 0) {
    fieldVec = user->Nvert;
    fieldIsVector = 0;
  } else {
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
             "Field '%s' not found in user context", fieldName);
  }

  /* Verify DM block size. */
  PetscInt expected_bs = fieldIsVector ? 3 : 1;
  PetscInt bs;
  ierr = DMGetBlockSize(user->fda, &bs); CHKERRQ(ierr);
  if (bs != expected_bs) {
    SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
             "Expected block size %d for field '%s', got %d", expected_bs, fieldName, bs);
  }

  /* Get local DM info for da. */
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(user->da,&info);

  /* Get coordinate array from da (read-only). */
  Vec Coor;
  ierr = DMGetCoordinatesLocal(user->da, &Coor); CHKERRQ(ierr);
  Cmpnts ***coor;
  ierr = DMDAVecGetArrayRead(user->fda, Coor, &coor); CHKERRQ(ierr);

  /* Allocate temporary array for interpolated cell-center coordinates.*/  

  /* Define the physical region of da. */
  PetscInt xs = info.xs, xe = info.xs + info.xm;
  PetscInt ys = info.ys, ye = info.ys + info.ym;
  PetscInt zs = info.zs, ze = info.zs + info.zm;
  PetscInt mx = info.mx, my = info.my, mz = info.mz;

  PetscInt lxs = xs,lxe = xe;
  PetscInt lys = ys,lye = ye;
  PetscInt lzs = zs,lze = ze;
 
  /*  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  */

  Cmpnts ***centcoor = NULL;
  ierr = Allocate3DArray(&centcoor, info.zm, info.ym, info.xm); CHKERRQ(ierr);

  /* Interpolate corner coordinates (from da) to cell centers. */
  ierr = InterpolateFieldFromCornerToCenter(coor, centcoor, user); CHKERRQ(ierr);

  if (fieldIsVector) {
    Cmpnts ***vecField;
    ierr = DMDAVecGetArray(user->fda, fieldVec, &vecField); CHKERRQ(ierr);
    /*--- Update interior cells ---*/
    for (PetscInt k = lzs; k < lze ; k++) {
      for (PetscInt j = lys; j < lye; j++) {
        for (PetscInt i = lxs; i < lxe; i++) {
          ierr = SetLocalCartesianField(fieldName,
                    &vecField[k][j][i],
                    &centcoor[k - lzs][j - lys][i - lxs]); CHKERRQ(ierr);
        }
      }
    }
    /*--- Update ghost cells ---*/
    

    ierr = DMDAVecRestoreArray(user->fda, fieldVec, &vecField); CHKERRQ(ierr);
  } else {
    PetscReal ***scalarField;
    ierr = DMDAVecGetArray(user->fda, fieldVec, &scalarField); CHKERRQ(ierr);
    for (PetscInt k = lzs; k < lze - 1; k++) {
      for (PetscInt j = lys; j < lye - 1; j++) {
        for (PetscInt i = lxs; i < lxe - 1; i++) {
          ierr = SetLocalCartesianField(fieldName,
                    &scalarField[k][j][i],
                    &centcoor[k - zs][j - ys][i - xs]); CHKERRQ(ierr);
        }
      }
    }
    /*-- Update ghost cells --*/

    ierr = DMDAVecRestoreArray(user->fda, fieldVec, &scalarField); CHKERRQ(ierr);
  }

  ierr = Deallocate3DArray(centcoor, info.zm, info.ym); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(user->fda, Coor, &coor); CHKERRQ(ierr);

  LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,
       "SetAnalyticalCartesianField: Completed for field '%s' on rank %d.\n", fieldName, rank);
  return 0;
}

/**
 * @brief Applies the analytical solution to the position vector.
 *
 * This function updates each entry in the provided PETSc vector by computing its sine,
 * thereby replacing each position with sin(position).
 *
 * @param tempVec The PETSc Vec containing particle positions which will be used to store velocities.
 * @return PetscErrorCode Returns 0 on success.
 */
PetscErrorCode SetAnalyticalSolution(Vec tempVec)
 {
     PetscErrorCode ierr;
     PetscInt nParticles;
     PetscReal *vels;
 
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "SetAnalyticalSolution - Starting analytical solution computation.\n");
 
     ierr = VecGetLocalSize(tempVec, &nParticles); CHKERRQ(ierr);
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "SetAnalyticalSolution - Number of local particles: %d.\n", nParticles);
 
     ierr = VecGetArray(tempVec, &vels); CHKERRQ(ierr);
     for (PetscInt i = 0; i < nParticles; i++) {
         vels[i] = sin(vels[i]);
     }
     ierr = VecRestoreArray(tempVec, &vels); CHKERRQ(ierr);
 
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "SetAnalyticalSolution - Completed analytical solution computation.\n");
     return 0;
 }
