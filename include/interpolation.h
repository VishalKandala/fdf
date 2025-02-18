#ifndef INTERPOLATION_H
#define INTERPOLATION_H

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

// Macros and constants
#define NUM_WEIGHTS 8 // Number of weights in trilinear interpolation

#define InterpolateFieldFromCornerToCenter(field, centfield, info)          \
  _Generic((field),                                                        \
    PetscReal ***: InterpolateFieldFromCornerToCenter_Scalar,                \
    Cmpnts ***:    InterpolateFieldFromCornerToCenter_Vector                  \
  )(field, centfield, info)


#define InterpolateFieldFromCenterToCorner(centfield, field, info)         \
  _Generic((centfield),                                                    \
    PetscReal ***: InterpolateFieldFromCenterToCorner_Scalar,                \
    Cmpnts ***:    InterpolateFieldFromCenterToCorner_Vector                  \
  )(centfield, field, info)

// Function declarations

/**
 * @brief Interpolates particle velocities using trilinear interpolation.
 *
 * @param[in] user Pointer to the user-defined context containing grid and swarm information.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InterpolateParticleVelocities(UserCtx *user);

/**
 * @brief Safely interpolate a scalar field from corners to centers, including boundary checks.
 *
 * For each cell center at (i,j,k), we consider up to eight corners around it:
 * (i-1..i, j-1..j, k-1..k). We check each candidate corner to ensure
 * 0 <= ci < mx, 0 <= cj < my, 0 <= ck < mz, where mx,my,mz are the global array dimensions
 * obtained from the DMDA. If a corner index is out of range, we skip it.
 *
 * This prevents out-of-bounds reads when xs=0 or near domain edges.
 * On boundaries, fewer corners are used, resulting in a partial average.
 *
 * @param[in]  field     3D array (PetscReal ***) containing corner values
 *                       (allocated with dimension [mz][my][mx])
 * @param[out] centfield 3D array (PetscReal ***) to store the interpolated cell-center values
 *                       (allocated with dimension [mz][my][mx])
 * @param[in]  info      DMDALocalInfo for the cell-center DMDA, containing local indices xs..xe, ys..ye, zs..ze
 *
 * @return PetscErrorCode
 */
PetscErrorCode InterpolateFieldFromCornerToCenter_Scalar(PetscReal ***field,
                                                  PetscReal ***centfield,
                                                  DMDALocalInfo *info);

/**
 * @brief Safely interpolate a vector field (Cmpnts) from corners to centers, including boundary checks.
 *
 * Similar to the scalar version, but each corner contributes three components (x,y,z).
 * On boundaries, fewer than eight corners may be valid. The code checks each candidate
 * corner index before reading field[ck][cj][ci].
 *
 * @param[in]  field     3D array (Cmpnts ***) of corner vectors
 * @param[out] centfield 3D array (Cmpnts ***) for interpolated cell-center vectors
 * @param[in]  info      DMDALocalInfo for the cell-center DMDA
 *
 * @return PetscErrorCode
 */
PetscErrorCode InterpolateFieldFromCornerToCenter_Vector(Cmpnts ***field,
                                                  Cmpnts ***centfield,
                                                  DMDALocalInfo *info);


/**
 * @brief Safely interpolate a scalar field from cell centers to corners, including boundary checks.
 *
 * For each corner at (i,j,k), we consider up to eight adjacent cell centers
 * at indices (i..i-1, j..j-1, k..k-1). If a center is out of valid [xs..xe), [ys..ye), etc.,
 * we skip it. This ensures no out-of-bounds read when i=0 or i=mx-1, etc.
 *
 * @param[in]  centfield 3D array (PetscReal ***) of cell-center values
 * @param[out] field     3D array (PetscReal ***) to store interpolated corner values
 * @param[in]  info      DMDALocalInfo for the cell-center DMDA
 *
 * @return PetscErrorCode
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Scalar(PetscReal ***field,
                                                  PetscReal ***centfield,
                                                  DMDALocalInfo *info);

/**
 * @brief Safely interpolate a vector field (Cmpnts) from cell centers to corners, including boundary checks.
 *
 * We loop over each corner (i,j,k) and consider up to eight cell centers
 * (i+di, j+dj, k+dk) in {0,1}. Each valid center is averaged.
 * On boundaries, partial stencils are used (skipping out-of-range indices).
 *
 * @param[in]  centfield 3D array (Cmpnts ***) of cell-center vectors
 * @param[out] field     3D array (Cmpnts ***) for corner vectors
 * @param[in]  info      DMDALocalInfo for the cell-center DMDA
 *
 * @return PetscErrorCode
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Vector(Cmpnts ***field,
                                                  Cmpnts ***centfield,
                                                  DMDALocalInfo *info);



#endif // INTERPOLATION_H
