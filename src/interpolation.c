/**
 * @file interpolation.c  //  Particle In Cell main.
 * @brief Main program for DMSwarm interpolation using the fdf-curvIB method.
 *
 * Initializes a particle swarm, reads velocity fields, and performs particle-grid interpolation.
 */

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscdmcomposite.h>

// Include the updated headers
//#include "common.h"         // Shared type definitions
//#include "ParticleSwarm.h"  // Particle swarm functions
//#include "walkingsearch.h"  // Particle location functions
//#include "grid.h"           // Grid functions
//#include "logging.h"        // Logging macros
//#include "io.h"             // Data Input and Output functions
#include "interpolation.h"


#define NUM_WEIGHTS 8 // Number of weights in trilinear interpolation

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
                                                                  DMDALocalInfo *info)
{
  PetscErrorCode ierr;
  PetscInt       mx, my, mz;             /* Global grid sizes in x,y,z   */
  PetscInt       xs, xe, ys, ye, zs, ze; /* Local cell-center range      */
  PetscInt       i, j, k, di, dj, dk;

  xs = info->xs; xe = info->xs + info->xm; /* i in [xs, xe) */
  ys = info->ys; ye = info->ys + info->ym; /* j in [ys, ye) */
  zs = info->zs; ze = info->zs + info->zm; /* k in [zs, ze) */
  mx = info->mx; 
  my = info->my; 
  mz = info->mz;

  /* Loop over all cell centers in the local patch. */
  for (k = zs; k < ze; k++) {
    for (j = ys; j < ye; j++) {
      for (i = xs; i < xe; i++) {

        /* We'll average up to eight corners around (i,j,k). */
        PetscReal sum   = 0.0;
        PetscInt  count = 0;

        /* Each corner is offset by di, dj, dk in {0,1}, but we subtract 1 from i, j, k first. */
        for (dk = 0; dk < 2; dk++) {
          for (dj = 0; dj < 2; dj++) {
            for (di = 0; di < 2; di++) {
              PetscInt ci = i - 1 + di; /* corner index in x */
              PetscInt cj = j - 1 + dj; /* corner index in y */
              PetscInt ck = k - 1 + dk; /* corner index in z */

              /* Check if this corner index is within valid bounds: 0 <= ci < mx, etc. */
              if (ci >= 0 && ci < mx &&
                  cj >= 0 && cj < my &&
                  ck >= 0 && ck < mz)
              {
                sum += field[ck][cj][ci];
                count++;
              }
            }
          }
        }

        /* Average the sum over however many corners were valid. */
        if (count > 0) {
          centfield[k][j][i] = sum / (PetscReal)count;
        } else {
          /* For corner cases with no valid neighbor, set 0 or apply a boundary condition. */
          centfield[k][j][i] = 0.0;
        }
      }
    }
  }

  return 0;
}

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
                                                         DMDALocalInfo *info)
{
  PetscErrorCode ierr;
  PetscInt       mx, my, mz;
  PetscInt       xs, xe, ys, ye, zs, ze;
  PetscInt       i, j, k, di, dj, dk;

  xs = info->xs; xe = info->xs + info->xm;
  ys = info->ys; ye = info->ys + info->ym;
  zs = info->zs; ze = info->zs + info->zm;
  mx = info->mx; my = info->my; mz = info->mz;

  for (k = zs; k < ze; k++) {
    for (j = ys; j < ye; j++) {
      for (i = xs; i < xe; i++) {

        Cmpnts sum = {0.0, 0.0, 0.0};
        PetscInt count = 0;

        for (dk = 0; dk < 2; dk++) {
          for (dj = 0; dj < 2; dj++) {
            for (di = 0; di < 2; di++) {
              PetscInt ci = i - 1 + di;
              PetscInt cj = j - 1 + dj;
              PetscInt ck = k - 1 + dk;
              if (ci >= 0 && ci < mx &&
                  cj >= 0 && cj < my &&
                  ck >= 0 && ck < mz)
              {
                sum.x += field[ck][cj][ci].x;
                sum.y += field[ck][cj][ci].y;
                sum.z += field[ck][cj][ci].z;
                count++;
              }
            }
          }
        }

        if (count > 0) {
          centfield[k][j][i].x = sum.x / (PetscReal)count;
          centfield[k][j][i].y = sum.y / (PetscReal)count;
          centfield[k][j][i].z = sum.z / (PetscReal)count;
        } else {
          centfield[k][j][i].x = 0.0;
          centfield[k][j][i].y = 0.0;
          centfield[k][j][i].z = 0.0;
        }
      }
    }
  }

  return 0;
}

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
PetscErrorCode InterpolateFieldFromCenterToCorner_Scalar(PetscReal ***centfield,
                                                                  PetscReal ***field,
                                                                  DMDALocalInfo *info)
{
  PetscErrorCode ierr;
  PetscInt       mx, my, mz;
  PetscInt       xs, xe, ys, ye, zs, ze;
  PetscInt       i, j, k, di, dj, dk;

  xs = info->xs; xe = info->xs + info->xm;
  ys = info->ys; ye = info->ys + info->ym;
  zs = info->zs; ze = info->zs + info->zm;
  mx = info->mx; my = info->my; mz = info->mz;

  /*
    The corner array is typically 1 cell larger in each dimension.
    We loop over corner indices.
    Here, we'll assume the local corner domain can be [xs-1..xe], etc.
    We'll just do a naive approach: let i go from (xs-1) to xe (not inclusive).
    We'll still do boundary checks inside the loop.
  */

  PetscInt xs_corner = xs - 1;
  PetscInt ys_corner = ys - 1;
  PetscInt zs_corner = zs - 1;
  PetscInt xe_corner = xe; 
  PetscInt ye_corner = ye; 
  PetscInt ze_corner = ze; 

  for (k = zs_corner; k < ze_corner; k++) {
    for (j = ys_corner; j < ye_corner; j++) {
      for (i = xs_corner; i < xe_corner; i++) {

        PetscReal sum   = 0.0;
        PetscInt  count = 0;

        /* The adjacent cell centers are (i + di, j + dj, k + dk) for di,dj,dk in {0,1}. */
        for (dk = 0; dk < 2; dk++) {
          for (dj = 0; dj < 2; dj++) {
            for (di = 0; di < 2; di++) {
              PetscInt ic = i + di;
              PetscInt jc = j + dj;
              PetscInt kc = k + dk;
              /* Check if (ic, jc, kc) is within the local cell-center domain. */
              if (ic >= xs && ic < xe &&
                  jc >= ys && jc < ye &&
                  kc >= zs && kc < ze)
              {
                sum += centfield[kc][jc][ic];
                count++;
              }
            }
          }
        }

        if (count > 0) {
          field[k][j][i] = sum / (PetscReal)count;
        } else {
          /* If no cell center is valid, set to 0 or apply boundary condition. */
          field[k][j][i] = 0.0;
        }
      }
    }
  }

  return 0;
}

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
PetscErrorCode InterpolateFieldFromCenterToCorner_Vector(Cmpnts ***centfield,
                                                         Cmpnts ***field,
                                                         DMDALocalInfo *info)
{
  PetscErrorCode ierr;
  PetscInt       mx, my, mz;
  PetscInt       xs, xe, ys, ye, zs, ze;
  PetscInt       i, j, k, di, dj, dk;

  xs = info->xs; xe = info->xs + info->xm;
  ys = info->ys; ye = info->ys + info->ym;
  zs = info->zs; ze = info->zs + info->zm;
  mx = info->mx; my = info->my; mz = info->mz;

  /* We assume the corner domain is slightly bigger: [xs-1..xe], etc., with boundary checks. */
  PetscInt xs_corner = xs - 1;
  PetscInt xe_corner = xe;
  PetscInt ys_corner = ys - 1;
  PetscInt ye_corner = ye;
  PetscInt zs_corner = zs - 1;
  PetscInt ze_corner = ze;

  for (k = zs_corner; k < ze_corner; k++) {
    for (j = ys_corner; j < ye_corner; j++) {
      for (i = xs_corner; i < xe_corner; i++) {

        Cmpnts sum  = {0.0, 0.0, 0.0};
        PetscInt count = 0;

        for (dk = 0; dk < 2; dk++) {
          for (dj = 0; dj < 2; dj++) {
            for (di = 0; di < 2; di++) {
              PetscInt ic = i + di;
              PetscInt jc = j + dj;
              PetscInt kc = k + dk;
              if (ic >= xs && ic < xe &&
                  jc >= ys && jc < ye &&
                  kc >= zs && kc < ze)
              {
                sum.x += centfield[kc][jc][ic].x;
                sum.y += centfield[kc][jc][ic].y;
                sum.z += centfield[kc][jc][ic].z;
                count++;
              }
            }
          }
        }

        if (count > 0) {
          field[k][j][i].x = sum.x / (PetscReal)count;
          field[k][j][i].y = sum.y / (PetscReal)count;
          field[k][j][i].z = sum.z / (PetscReal)count;
        } else {
          field[k][j][i].x = 0.0;
          field[k][j][i].y = 0.0;
          field[k][j][i].z = 0.0;
        }
      }
    }
  }

  return 0;
}


/**
 * @brief Computes the trilinear interpolation weights from the interpolation coefficients.
 *
 * This function computes the weights for trilinear interpolation at the eight corners of a cell
 * using the interpolation coefficients provided along the x, y, and z directions.
 *
 * @param[in]  a1 Interpolation coefficient along the x-direction (normalized coordinate within the cell).
 * @param[in]  a2 Interpolation coefficient along the y-direction (normalized coordinate within the cell).
 * @param[in]  a3 Interpolation coefficient along the z-direction (normalized coordinate within the cell).
 * @param[out] w  Array of 8 weights, each corresponding to one corner of the cell.
 *
 * @note
 * - The coefficients `a1`, `a2`, and `a3` should be in the range [0, 1].
 * - The order of weights corresponds to the eight corners of a hexahedral cell.
 */
static inline void ComputeTrilinearWeights(PetscReal a1, PetscReal a2, PetscReal a3, PetscReal *w) {
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ComputeTrilinearWeights: Computing weights for a1=%f, a2=%f, a3=%f.\n", a1, a2, a3);

    // Ensure a1, a2, a3 are within [0,1]
    a1 = PetscMax(0.0, PetscMin(1.0, a1));
    a2 = PetscMax(0.0, PetscMin(1.0, a2));
    a3 = PetscMax(0.0, PetscMin(1.0, a3));

    // Compute complementary coefficients (1 - a)
    PetscReal oa1 = 1.0 - a1;  
    PetscReal oa2 = 1.0 - a2;  
    PetscReal oa3 = 1.0 - a3;

    // Compute weights for each corner of the cell
    w[0] = a1  * a2  * a3;   // Front-bottom-left
    w[1] = a1  * oa2 * oa3;  // Front-top-right
    w[2] = oa1 * a2  * oa3;  // Back-bottom-right
    w[3] = a1  * a2  * oa3;  // Front-bottom-right
    w[4] = oa1 * oa2 * a3;   // Back-top-left
    w[5] = oa1 * oa2 * oa3;  // Back-top-right
    w[6] = oa1 * a2  * a3;   // Back-bottom-left
    w[7] = a1  * oa2 * a3;   // Front-top-left

    // Log the computed weights for debugging
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ComputeTrilinearWeights: Weights computed - "
        "w0=%f, w1=%f, w2=%f, w3=%f, w4=%f, w5=%f, w6=%f, w7=%f. \n",
        w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7]);
}

/**
 * @brief Computes the trilinear interpolated velocity at a given point.
 *
 * @param[in]  ucat   3D array of velocity field from a DMDA (indexed as [k][j][i]), each cell of type Cmpnts.
 * @param[in]  i,j,k  Integral cell indices (the "lower" corner in each dimension).
 * @param[in]  a1,a2,a3  Normalized coordinates within the cell ([0,1] range).
 * @param[out] vel    Pointer to a Cmpnts struct that will hold the interpolated velocity (x, y, z).
 *
 * The function uses the 8-corner trilinear formula with `ComputeTrilinearWeights()`.
 * If a different scheme is desired, implement a new function with the same interface.
 */
static inline PetscErrorCode InterpolateTrilinearVelocity(
    Cmpnts ***ucat,
    PetscInt i, PetscInt j, PetscInt k,
    PetscReal a1, PetscReal a2, PetscReal a3,
    Cmpnts *vel)
{
    PetscFunctionBegin; // PETSc macro for error/stack tracing

    // Compute the 8 corner weights
    PetscReal wcorner[8];
    ComputeTrilinearWeights(a1, a2, a3, wcorner);

    // Offsets for cell corners
    PetscInt i1 = i + 1;
    PetscInt j1 = j + 1;
    PetscInt k1 = k + 1;

    // Initialize velocity
    vel->x = 0.0;
    vel->y = 0.0;
    vel->z = 0.0;

    // Corner 0 => (i, j, k)
    vel->x += wcorner[0] * ucat[k ][j ][i ].x;
    vel->y += wcorner[0] * ucat[k ][j ][i ].y;
    vel->z += wcorner[0] * ucat[k ][j ][i ].z;

    // Corner 1 => (i+1, j, k)
    vel->x += wcorner[1] * ucat[k ][j ][i1].x;
    vel->y += wcorner[1] * ucat[k ][j ][i1].y;
    vel->z += wcorner[1] * ucat[k ][j ][i1].z;

    // Corner 2 => (i, j+1, k)
    vel->x += wcorner[2] * ucat[k ][j1][i ].x;
    vel->y += wcorner[2] * ucat[k ][j1][i ].y;
    vel->z += wcorner[2] * ucat[k ][j1][i ].z;

    // Corner 3 => (i+1, j+1, k)
    vel->x += wcorner[3] * ucat[k ][j1][i1].x;
    vel->y += wcorner[3] * ucat[k ][j1][i1].y;
    vel->z += wcorner[3] * ucat[k ][j1][i1].z;

    // Corner 4 => (i, j, k+1)
    vel->x += wcorner[4] * ucat[k1][j ][i ].x;
    vel->y += wcorner[4] * ucat[k1][j ][i ].y;
    vel->z += wcorner[4] * ucat[k1][j ][i ].z;

    // Corner 5 => (i+1, j, k+1)
    vel->x += wcorner[5] * ucat[k1][j ][i1].x;
    vel->y += wcorner[5] * ucat[k1][j ][i1].y;
    vel->z += wcorner[5] * ucat[k1][j ][i1].z;

    // Corner 6 => (i, j+1, k+1)
    vel->x += wcorner[6] * ucat[k1][j1][i ].x;
    vel->y += wcorner[6] * ucat[k1][j1][i ].y;
    vel->z += wcorner[6] * ucat[k1][j1][i ].z;

    // Corner 7 => (i+1, j+1, k+1)
    vel->x += wcorner[7] * ucat[k1][j1][i1].x;
    vel->y += wcorner[7] * ucat[k1][j1][i1].y;
    vel->z += wcorner[7] * ucat[k1][j1][i1].z;


    LOG_ALLOW(LOCAL, LOG_DEBUG, 
    "InterpolateTrilinearVelocity: Interpolated velocity at (i=%d, j=%d, k=%d) with a1=%.6f, a2=%.6f, a3=%.6f -> (x=%.6f, y=%.6f, z=%.6f).\n",
    i, j, k, a1, a2, a3, vel->x, vel->y, vel->z);


    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Completed particle velocity interpolation across all ranks.\n");


    PetscFunctionReturn(0);
}

/**
 * @brief Performs trilinear interpolation of velocities from grid to particles using interpolation coefficients.
 *
 * This function interpolates velocities for particles based on the velocity field defined on the computational grid.
 * It retrieves cell indices for each particle from the "DMSwarm_CellID" field and assigns the interpolated velocity
 * using trilinear interpolation from the surrounding grid cells. The function handles boundary checks and ensures
 * that particles outside the valid grid range are appropriately managed.
 *
 * Key Steps:
 * 1. Retrieve the number of local particles and their associated data fields (cell indices and velocities).
 * 2. Map the global velocity field (`Ucat`) to the local portion for efficient access.
 * 3. Retrieve interpolation coefficients (`a1`, `a2`, `a3`) for each particle.
 * 4. Compute trilinear interpolation weights and interpolate velocities from the surrounding grid cells.
 * 5. Handle edge cases where particles may lie outside the valid grid range.
 * 6. Restore all data fields and ensure PETSc arrays and vectors are correctly finalized.
 *
 * @param[in] user Pointer to the `UserCtx` structure containing:
 *                 - `user->da`: DMDA for the grid.
 *                 - `user->swarm`: DMSwarm for particles.
 *                 - `user->Ucat`: Global velocity vector on the grid.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InterpolateParticleVelocities(UserCtx *user) {
    PetscErrorCode ierr;
    PetscInt n_local;             // Number of local particles
    PetscInt64 *cellIDs = NULL;     // Array to store cell indices for each particle
    PetscReal *velocities = NULL; // Array to store particle velocities
    PetscReal *weights = NULL;    // Array to store interpolation weights
    Cmpnts ***ucat;               // 3D array to map local grid velocities
    Cmpnts uinterp;               // Temporary variable to store interpolated velocities.
    PetscInt i,j,k;
    PetscReal a1,a2,a3;           // The weights of a particles are stored here for trilinear coefficient calculation.
    DM fda = user->fda;           // Field DA 
    DM swarm = user->swarm;       // DMSwarm for the particles


    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Starting particle velocity interpolation.\n");

    // Verify global velocity field
    PetscReal max_val;
    ierr = VecNorm(user->Ucat, NORM_INFINITY, &max_val); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Global velocity field Ucat maximum magnitude: %f \n", max_val);

    // Retrieve the number of local particles
    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Found %d local particles.\n", n_local);

    // Retrieve particle data fields
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr); // Ensure 'weight' field exists

    // Access the local portion of the global velocity vector (Ucat) using 'fda'
    ierr = DMDAVecGetArrayRead(user->fda,user->Ucat,&ucat);CHKERRQ(ierr);


    LOG_ALLOW(GLOBAL, LOG_DEBUG, "InterpolateParticleVelocities: Starting velocity assignment for particles.\n");

    // Log grid dimensions
    LOG_ALLOW(GLOBAL, LOG_INFO, "Grid dimensions: mx=%d, my=%d, mz=%d \n", user->info.mx, user->info.my, user->info.mz);

    // Loop over all local particles
    for (PetscInt p = 0; p < n_local; p++) {
        // Retrieve cell indices for the particle
        i = cellIDs[3 * p];
        j = cellIDs[3 * p + 1];
        k = cellIDs[3 * p + 2];

	LOG_LOOP_ALLOW(GLOBAL, LOG_DEBUG, p, 10, "Particle %d: Host Cell = (%d, %d, %d)\n", p, i, j, k);

	// Clamp i, j, k to [0..mx-2], [0..my-2], [0..mz-2]
	if (i >= user->info.mx - 1) i = user->info.mx - 2;
	if (j >= user->info.my - 1) j = user->info.my - 2;
	if (k >= user->info.mz - 1) k = user->info.mz - 2;

        // Validate cell indices (boundary check)
        if (i < 0 || j < 0 || k < 0 || i >= user->info.mx || j >= user->info.my || k >= user->info.mz) {
            LOG_ALLOW(GLOBAL, LOG_WARNING, "Particle %d has invalid cell indices (%d, %d, %d)\n. Skipping interpolation.\n", p, i, j, k);
            velocities[3 * p    ] = 0.0;
            velocities[3 * p + 1] = 0.0;
            velocities[3 * p + 2] = 0.0;
            continue;
        }

        // Retrieve a1, a2, a3 from the 'weights' field (if that's where you're storing them)
        a1 = weights[3*p + 0];
        a2 = weights[3*p + 1];
        a3 = weights[3*p + 2];

        
	// Apply interpolation method to obtain velocity at the particle location
	//  ierr = InterpolateTrilinearVelocity(ucat,i,j,k,a1,a2,a3,&uinterp);
	
        // Assign interpolated velocity to the particle

        // zeroth order Interpolation
        
        velocities[3 * p]     = ucat[k+1][j+1][i+1].x; // u-component
        velocities[3 * p + 1] = ucat[k+1][j+1][i+1].y; // v-component
        velocities[3 * p + 2] = ucat[k+1][j+1][i+1].z; // w-component
       
	LOG_LOOP_ALLOW(GLOBAL, LOG_DEBUG, p, 10, "Particle %d: Interpolated velocity: (velocity.x=%f, velocity.y=%f, velocity.z=%f).\n",
                        p, velocities[3 * p], velocities[3 * p + 1], velocities[3 * p + 2]);	
    }

    // Restore the local velocity array
    ierr = DMDAVecRestoreArrayRead(fda,user->Ucat, &ucat); CHKERRQ(ierr);

    // Restore particle data fields
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Particle velocity interpolation completed.\n");

    // Ensure all ranks finish interpolation before proceeding
    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: All ranks completed particle interpolation.\n");

    return 0;
}
