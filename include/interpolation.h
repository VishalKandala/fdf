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
#include "setup.h"          // Setup functions that are used across the codebase

// Macros and constants
#define NUM_WEIGHTS 8 // Number of weights in trilinear interpolation

/**
 * @brief Generic macro to call the appropriate interpolation function based on the field type.
 *
 * This macro will select either the scalar or vector interpolation function based on
 * the type of the 'field' pointer. It also performs a compile-time check that 'centfield'
 * is of the same type as 'field'. If the types are not the same, a compile-time error is produced.
 *
 * Usage:
 *   InterpolateFieldFromCornerToCenter(field, centfield, user);
 */
#define InterpolateFieldFromCornerToCenter(field, centfield, user)                         \
  ( (void)sizeof(char[1 - 2*!!(!__builtin_types_compatible_p(typeof(field), typeof(centfield)))]),\
    _Generic((field),                                                                       \
      PetscReal ***: InterpolateFieldFromCornerToCenter_Scalar,                             \
      Cmpnts ***:    InterpolateFieldFromCornerToCenter_Vector                                \
    )(field, centfield, user) )


#define InterpolateFieldFromCenterToCorner(centfield, field, info)         \
  _Generic((centfield),                                                    \
    PetscReal ***: InterpolateFieldFromCenterToCorner_Scalar,                \
    Cmpnts ***:    InterpolateFieldFromCenterToCorner_Vector                  \
  )(centfield, field, info)


/**
 * @brief Macro that calls either the scalar or vector piecewise interpolation function
 *        based on the type of the `fieldPtr` parameter (3D array).
 *
 * Usage example:
 *
 *   // For scalar:
 *   PetscReal ***fieldScal;
 *   PetscReal outVal;
 *   PieceWiseLinearInterpolation(fieldName, fieldScal, i, j, k, &outVal);
 *
 *   // For vector:
 *   Cmpnts ***fieldVec;
 *   Cmpnts vec;
 *   PieceWiseLinearInterpolation(fieldName, fieldVec, i, j, k, &vec);
 */
#define PieceWiseLinearInterpolation(fieldName, fieldPtr, i, j, k, outPtr) \
  _Generic((fieldPtr),                                                     \
    PetscReal ***: PieceWiseLinearInterpolation_Scalar,                    \
    Cmpnts      ***: PieceWiseLinearInterpolation_Vector                   \
  )(fieldName, fieldPtr, i, j, k, outPtr)

/**
 * @brief Macro that calls either the scalar or vector trilinear interpolation function
 *        based on the type of the `fieldPtr` parameter (3D array).
 *
 * Usage example:
 * 
 *   PetscReal result;
 *   Cmpnts vec;
 *   
 *   // For scalars:
 *   TrilinearInterpolation(fieldName, fieldScal, i, j, k, a1, a2, a3, &result);
 *
 *   // For vectors:
 *   TrilinearInterpolation(fieldName, fieldVec, i, j, k, a1, a2, a3, &vec);
 */
#define TrilinearInterpolation(fieldName, fieldPtr, i, j, k, a1, a2, a3, outPtr)           \
  _Generic((fieldPtr),                                                                     \
    PetscReal ***: TrilinearInterpolation_Scalar,                                          \
    Cmpnts      ***: TrilinearInterpolation_Vector                                         \
  )(fieldName, fieldPtr, i, j, k, a1, a2, a3, outPtr)

// Function declarations

/**
 * @brief Computes the trilinear interpolated scalar at a given point.
 *
 * @param[in]  fieldName A string representing the name of the scalar field (e.g., "temperature").
 * @param[in]  fieldScal 3D array of the field from a DMDA (indexed as [k][j][i]),
 *                      each cell a PetscReal.
 * @param[in]  i, j, k   Integral cell indices (the "lower" corner in each dimension).
 * @param[in]  a1, a2, a3 Normalized coordinates within the cell ([0,1] range).
 * @param[out] val       Pointer to a PetscReal that will store the interpolated scalar.
 *
 * This function uses the standard 8-corner trilinear formula via `ComputeTrilinearWeights()`.
 * If a different scheme is desired, implement a new function with the same interface.
 */
PetscErrorCode TrilinearInterpolation_Scalar(
    const char   *fieldName,
    PetscReal  ***fieldScal,
    PetscInt      i,
    PetscInt      j,
    PetscInt      k,
    PetscReal     a1,
    PetscReal     a2,
    PetscReal     a3,
    PetscReal    *val);

/**
 * @brief Computes the trilinear interpolated vector (e.g., velocity) at a given point.
 *
 * @param[in]  fieldName  A string representing the name of the vector field (e.g., "velocity").
 * @param[in]  fieldVec   3D array of the field from a DMDA (indexed as [k][j][i]),
 *                        each cell of type Cmpnts.
 * @param[in]  i, j, k    Integral cell indices (the "lower" corner in each dimension).
 * @param[in]  a1, a2, a3 Normalized coordinates within the cell ([0,1] range).
 * @param[out] vec        Pointer to a Cmpnts struct that will store the interpolated vector (x, y, z).
 *
 * This function uses the standard 8-corner trilinear formula via `ComputeTrilinearWeights()`.
 * If a different scheme is desired, implement a new function with the same interface.
 */
PetscErrorCode TrilinearInterpolation_Vector(
    const char   *fieldName,
    Cmpnts     ***fieldVec,
    PetscInt      i,
    PetscInt      j,
    PetscInt      k,
    PetscReal     a1,
    PetscReal     a2,
    PetscReal     a3,
    Cmpnts       *vec);

/**
 * @brief High-level function to interpolate one field (either scalar or vector) onto all local particles.
 *
 * Steps:
 *  1) Get blockSize from fieldGlobal (1 => scalar, 3 => vector).
 *  2) Retrieve local pointer from DMDAVecGetArray().
 *  3) Retrieve "DMSwarm_CellID", "weight", and swarmOutFieldName from the DMSwarm.
 *  4) Loop over local particles. For each:
 *     - read iCell,jCell,kCell
 *     - clamp or skip if out of range
 *     - read a1,a2,a3
 *     - call InterpolateSingleFieldForParticle(...) with the typed pointer
 *  5) Restore arrays and fields
 */
PetscErrorCode InterpolateEulerFieldToSwarm(
    UserCtx    *user,
    Vec         fieldGlobal,       /* DMDA Vec for the field */
    const char *fieldName,         /* e.g., "velocity", "temp" */
    const char *swarmOutFieldName); /* DMSwarm field for storing results */

/**
 * @brief Interpolates all relevant fields from the DMDA to the DMSwarm.
 *
 * Currently, it interpolates:
 *   - user->Ucat (vector field) into the DMSwarm field "swarmVelocity".
 *
 * To add more fields, duplicate the call to InterpolateOneFieldOverSwarm and provide:
 *   - The global Vec for that field (e.g. user->Tcat for temperature),
 *   - A human-readable field name (for logging),
 *   - A DMSwarm output field name (e.g. "swarmTemperature").
 *
 * @param[in,out] user Pointer to a UserCtx containing:
 *                     - user->da (DM for the grid),
 *                     - user->swarm (DMSwarm for particles),
 *                     - user->Ucat (Vec for the vector field),
 *                     - possibly more fields like user->Tcat, user->Pcat, etc.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InterpolateAllFieldsToSwarm(UserCtx *user);

/**
 * @brief Interpolates particle velocities using trilinear interpolation.
 *
 * @param[in] user Pointer to the user-defined context containing grid and swarm information.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */

PetscErrorCode InterpolateParticleVelocities(UserCtx *user);

/**
 * @brief Safely interpolate a scalar field from corner nodes (from the coordinate DM)
 *        to cell centers (from the cell-centered DM) using the provided UserCtx.
 *
 * For each cell center in the physical region of the cell-centered DM (fda), this function
 * averages the up to 8 surrounding scalar values from the coordinate DM (da). On boundaries,
 * where fewer corners are available, a partial average is computed.
 *
 * The coordinate DM (da) is built on corners (IM+1 x JM+1 x KM+1) while the cell-centered DM (fda)
 * covers the physical cells (IM x JM x KM). Index offsets are adjusted via DMDAGetLocalInfo.
 *
 * @param[in]  field     3D array of corner-based scalar data (from user->da).
 * @param[out] centfield 3D array for the interpolated cell-center scalar data (for user->fda).
 * @param[in]  user      User context containing:
 *                       - da  : DM for the coordinate (corner) data.
 *                       - fda : DM for the cell-centered data.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCornerToCenter_Scalar(
    PetscReal ***field,
    PetscReal ***centfield,
    UserCtx *user);

/**
 * @brief Safely interpolate a vector field from corner nodes (from the coordinate DM)
 *        to cell centers (from the cell-centered DM) using the provided UserCtx.
 *
 * For each cell center in the physical region of the cell-centered DM (fda), this function
 * averages the 8 surrounding corner values from the coordinate DM (da). The coordinate DM
 * (da) is built on corners (IM+1 x JM+1 x KM+1) while the cell-centered DM (fda) covers
 * the physical cells (IM x JM x KM). Index offsets are adjusted using DMDAGetLocalInfo.
 *
 * @param[in]  field     3D array of corner-based vector data (from user->da).
 * @param[out] centfield 3D array for the interpolated cell-center vector data (for user->fda).
 * @param[in]  user      User context containing:
 *                       - da  : DM for the coordinate (corner) data.
 *                       - fda : DM for the cell-centered data.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCornerToCenter_Vector(
    Cmpnts ***field,
    Cmpnts ***centfield,
    UserCtx *user);

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
