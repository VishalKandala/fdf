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
 * @brief Sets the local Cartesian velocity components based on coordinates.
 *
 * This function computes the velocity components by applying the sine function to the 
 * input coordinate values along the x, y, and z directions.
 *
 * @param[in,out] ucont Pointer to the `Cmpnts` structure where the velocity components will be stored.
 * @param[in]     coor  Pointer to the `Cmpnts` structure containing the input coordinate values.
 *
 * @return PetscErrorCode Returns 0 on successful execution, non-zero on failure.
 */
 static inline PetscErrorCode SetLocalCartesianVelocity(Cmpnts *ucat, Cmpnts *coor) {
    // Log input coordinate values for debugging purposes
    LOG_ALLOW(LOCAL, LOG_DEBUG, "SetLocalCartesianVelocity: Input Coordinates - x: %f, y: %f, z: %f \n", coor->x, coor->y, coor->z);

    // Compute velocity components as the sine of the coordinates
    ucat->x = sin(coor->x);
    ucat->y = sin(coor->y);
    ucat->z = sin(coor->z);

    // Log computed velocity components
    LOG_ALLOW(LOCAL, LOG_DEBUG, "SetLocalCartesianVelocity: Computed Velocity - x: %f, y: %f, z: %f \n", ucat->x, ucat->y, ucat->z);

    return 0;
}

/**
 * @brief Updates the local Cartesian velocity field based on interpolated coordinates.
 *
 * This function performs the following operations:
 * 1. Retrieves the local grid information from the DMDA associated with the user context.
 * 2. Accesses the local coordinate vector and retrieves a read-only array representation.
 * 3. Allocates memory for a temporary 3D array (`centcoor`) to store interpolated coordinate values at cell centers.
 * 4. Accesses the Cartesian velocity vector and obtains a writable array representation.
 * 5. Interpolates the coordinate values from cell corners to cell centers using `InterpolateFieldFromCornerToCenter`.
 * 6. Updates the Cartesian velocity values at each cell center based on the interpolated coordinates using `SetLocalCartesianVelocity`.
 * 7. Deallocates the temporary array (`centcoor`) and restores the arrays to maintain PETSc's internal state.
 *
 * @param[in,out] user Pointer to a `UserCtx` structure containing:
 *                     - `user->da`: DMDA for the grid.
 *                     - `user->fda`: DMDA for the Cartesian velocity field.
 *                     - `user->Ucat`: Cartesian velocity vector.
 *                     - `user->Coor`: Coordinate vector.
 *                     - `user->info`: Local DMDA grid information.
 *
 * @return PetscErrorCode 
 *         - Returns 0 on successful execution.
 *         - Non-zero error code on failure.
 *
 * @note 
 * - Assumes that the coordinate vector (`user->Coor`) has been properly initialized and populated before this function is called.
 * - The function is grid-aware and handles memory allocation and interpolation for local subdomains only.
 */

PetscErrorCode UpdateCartesianVelocity(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo info = user->info;
    Vec Coor;
    Cmpnts ***coor = NULL, ***centcoor = NULL, ***ucat = NULL;
    PetscInt rank;

    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;
 
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "UpdateCartesianVelocity: Starting on rank : %d \n",rank);

    /* 1) Access local coords (Coor) & velocity array (Ucat) */
    ierr = DMGetCoordinatesLocal(user->da, &Coor); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, Coor, &coor); CHKERRQ(ierr);

    ierr = DMDAVecGetArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);

    /* 2) Allocate scratch array 'centcoor' the same local size as [zs..ze, ys..ye, xs..xe]. */
    ierr = Allocate3DArray(&centcoor, ze, ye, xe); CHKERRQ(ierr);

    /* 3) Interpolate corner->center coords into centcoor */
    ierr = InterpolateFieldFromCornerToCenter(coor, centcoor, &info); CHKERRQ(ierr);

    /* 4) Fill interior cell velocity using 'SetLocalCartesianVelocity' with the center coords. 
       (Stop at xe-1, ye-1, ze-1 so i+1 doesn't go out-of-range) */
    for (PetscInt k = zs; k < ze - 1; k++) {
      for (PetscInt j = ys; j < ye - 1; j++) {
        for (PetscInt i = xs; i < xe - 1; i++) {
          ierr = SetLocalCartesianVelocity(&ucat[k][j][i], &centcoor[k][j][i]); CHKERRQ(ierr);
        }
      }
    }

    /* 5) Set boundary nodes to sine-of-corner-coords 
       i.e. if i==xe-1 or i==xs or similarly for j,k. 
       That ensures no boundary node is left zero. */
    for (PetscInt k = zs; k < ze; k++) {
      for (PetscInt j = ys; j < ye; j++) {
        for (PetscInt i = xs; i < xe; i++) {
          PetscBool isBoundary = PETSC_FALSE;
          if ((i == xe - 1) || (i == xs) ||
              (j == ye - 1) || (j == ys)  ||
              (k == ze - 1) || (k == zs)) {
            isBoundary = PETSC_TRUE;
          }
          if (isBoundary) {
            ierr = SetLocalCartesianVelocity(&ucat[k][j][i], &coor[k][j][i]); CHKERRQ(ierr);
          }
        }
      }
    }

    /* 6) Clean up */
    ierr = Deallocate3DArray(centcoor, ze, ye); CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, Coor, &coor); CHKERRQ(ierr);

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "UpdateCartesianVelocity: Completed on rank : %d \n",rank);
    return 0;
}
