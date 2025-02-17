#ifndef ANALYTICALSOLUTION_H
#define ANALYTICALSOLUTION_H

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscdmcomposite.h>

// Include additional headers
#include "common.h"         // Shared type definitions
#include "ParticleSwarm.h"  // Particle swarm functions
#include "walkingsearch.h"  // Particle location functions
#include "grid.h"           // Grid functions
#include "logging.h"        // Logging macros
#include "io.h"             // Data Input and Output functions

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
    DMDALocalInfo *info);


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
 PetscErrorCode Allocate3DArray(Cmpnts ****array, PetscInt nz, PetscInt ny, PetscInt nx);


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
PetscErrorCode Deallocate3DArray(Cmpnts ***array, PetscInt nz, PetscInt ny);

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
 **/

PetscErrorCode UpdateCartesianVelocity(UserCtx *user);

#endif // ANALYTICALSOLUTION_H