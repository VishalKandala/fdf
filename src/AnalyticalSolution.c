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
 * @brief Interpolates a vector field from cell corners to cell centers using simple averaging.
 *
 * This version loops only up to `xe-1, ye-1, ze-1`, ensuring `i+1, j+1, k+1` are in range.
 * The result is stored in `centfield[k][j][i]` at the same (k,j,i) offsets
 * (assuming the caller knows to shift or match indexing carefully).
 *
 * @param[in]  field     A 3D array of `Cmpnts` at cell corners (size at least [info->ze][info->ye][info->xe]).
 * @param[out] centfield A 3D array of `Cmpnts` into which cell-center data is written.
 *                       Must be allocated by the caller with same local extents or an offset approach.
 * @param[in]  info      DMDALocalInfo with local domain indices [xs..xe), etc.
 *
 * @return PetscErrorCode
 */

PetscErrorCode InterpolateFieldFromCornerToCenter(Cmpnts ***field,
                                                  Cmpnts ***centfield,
                                                  DMDALocalInfo *info)
{
    PetscInt i, j, k;
    PetscInt xs = info->xs, xe = info->xs + info->xm;
    PetscInt ys = info->ys, ye = info->ys + info->ym;
    PetscInt zs = info->zs, ze = info->zs + info->zm;

    /* we must stop at xe-1 so i+1 is valid local indexing */
    PetscInt iend = xe - 1, jend = ye - 1, kend = ze - 1;

    for (k = zs; k < kend; k++) {
        for (j = ys; j < jend; j++) {
            for (i = xs; i < iend; i++) {
                centfield[k][j][i].x =
                   ( field[k][j][i].x + field[k][j][i+1].x +
                     field[k][j+1][i].x + field[k][j+1][i+1].x +
                     field[k+1][j][i].x + field[k+1][j][i+1].x +
                     field[k+1][j+1][i].x + field[k+1][j+1][i+1].x ) / 8.0;

                centfield[k][j][i].y =
                   ( field[k][j][i].y + field[k][j][i+1].y +
                     field[k][j+1][i].y + field[k][j+1][i+1].y +
                     field[k+1][j][i].y + field[k+1][j][i+1].y +
                     field[k+1][j+1][i].y + field[k+1][j+1][i+1].y ) / 8.0;

                centfield[k][j][i].z =
                   ( field[k][j][i].z + field[k][j][i+1].z +
                     field[k][j+1][i].z + field[k][j+1][i+1].z +
                     field[k+1][j][i].z + field[k+1][j][i+1].z +
                     field[k+1][j+1][i].z + field[k+1][j+1][i+1].z ) / 8.0;
            }
        }
    }
    return 0;
}

/**
 * @brief Allocates a 3D array of `Cmpnts` structures.
 *
 * This function dynamically allocates memory for a 3D array of `Cmpnts` structures.
 * Each component of the array (x, y, z) is initialized to 0.
 *
 * @param[out] array Pointer to the 3D array to be allocated.
 * @param[in]  nz    Number of layers in the z-direction.
 * @param[in]  ny    Number of rows in the y-direction.
 * @param[in]  nx    Number of columns in the x-direction.
 *
 * @return PetscErrorCode Returns 0 on successful allocation, non-zero on failure.
 */
PetscErrorCode Allocate3DArray(Cmpnts ****array, PetscInt nz, PetscInt ny, PetscInt nx) {
    PetscErrorCode ierr;

    LOG_ALLOW(LOCAL, LOG_INFO, "Allocate3DArray: Allocating 3D array of size (%d x %d x %d).\n", nz, ny, nx);

    // Allocate memory for each dimension
    ierr = PetscMalloc1(nz, array); CHKERRQ(ierr);
    for (PetscInt k = 0; k < nz; k++) {
        ierr = PetscMalloc1(ny, &(*array)[k]); CHKERRQ(ierr);
        for (PetscInt j = 0; j < ny; j++) {
            ierr = PetscMalloc1(nx, &(*array)[k][j]); CHKERRQ(ierr);
            for (PetscInt i = 0; i < nx; i++) {
                // Initialize components to zero
                (*array)[k][j][i].x = 0.0;
                (*array)[k][j][i].y = 0.0;
                (*array)[k][j][i].z = 0.0;
            }
        }
    }

    LOG_ALLOW(LOCAL, LOG_INFO, "Allocate3DArray: Successfully allocated 3D array.\n");
    return 0;
}

/**
 * @brief Deallocates a 3D array of `Cmpnts` structures.
 *
 * This function frees the memory allocated for a 3D array of `Cmpnts` structures.
 *
 * @param[in]  array Pointer to the 3D array to be deallocated.
 * @param[in]  nz    Number of layers in the z-direction.
 * @param[in]  ny    Number of rows in the y-direction.
 *
 * @return PetscErrorCode Returns 0 on successful deallocation, non-zero on failure.
 */
PetscErrorCode Deallocate3DArray(Cmpnts ***array, PetscInt nz, PetscInt ny) {
    PetscErrorCode ierr;

    LOG_ALLOW(LOCAL, LOG_INFO, "Deallocate3DArray: Deallocating 3D array of size (%d x %d).\n", nz, ny);

    // Free memory layer by layer
    for (PetscInt k = 0; k < nz; k++) {
        for (PetscInt j = 0; j < ny; j++) {
            ierr = PetscFree(array[k][j]); CHKERRQ(ierr);
        }
        ierr = PetscFree(array[k]); CHKERRQ(ierr);
    }
    ierr = PetscFree(array); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_INFO, "Deallocate3DArray: Successfully deallocated 3D array.\n");
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
