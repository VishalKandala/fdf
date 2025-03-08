/**
 * @file interpolation.c
 * @brief Main program for DMSwarm interpolation using the fdf-curvIB method.
 *
 * Provides routines for interpolation between corner-based and center-based
 * fields in the cell-centered DM (fda), plus partial usage examples for
 * DMSwarm-based field sampling.
 */

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscdmcomposite.h>

#include "interpolation.h"

// Number of weights used in certain trilinear interpolation examples
#define NUM_WEIGHTS 8 
// -----------------------------------------------------------------------------
// Interpolation: Corner -> Center (vector)
// -----------------------------------------------------------------------------
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
    Cmpnts ***field,         /* coordinate DM array (da, ghosted) */
    Cmpnts ***centfield,     /* output: interpolated interior cell–center values for fda */
    UserCtx *user)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
    "InterpolateFieldFromCornerToCenter_Vector - Rank %d starting interpolation.\n", rank);

  /* Get local info for da (owned cells) */
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

  /* For da, get ghost region info instead of owned-only info */
  PetscInt gxs, gys, gzs, gxm, gym, gzm;
  ierr = DMDAGetGhostCorners(user->da, &gxs, &gys, &gzs, &gxm, &gym, &gzm); CHKERRQ(ierr);
  LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
    "Rank %d -> DM-da Ghost Corners: gxs=%d, gxm=%d, gys=%d, gym=%d, gzs=%d, gzm=%d\n",
    rank, gxs, gxm, gys, gym, gzs, gzm);

  /* Compute physical region for fda */
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

  for (PetscInt k = lzs; k < lze; k++) {
    for (PetscInt j = lys; j < lye; j++) {
      for (PetscInt i = lxs; i < lxe; i++) {

        Cmpnts sum = {0.0, 0.0, 0.0};
        PetscInt count = 0;

        /* For cell center (i,j,k) in fda, its 8 surrounding corners in da have global indices
           from (i-1,j-1,k-1) to (i,j,k). Use ghost region information from da.
         */
        for (PetscInt dk = 0; dk < 2; dk++) {
          for (PetscInt dj = 0; dj < 2; dj++) {
            for (PetscInt di = 0; di < 2; di++) {
             
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

	        LOG_LOOP_ALLOW(GLOBAL,LOG_DEBUG,i*j*k,1000," Rank %d| i,j,k - %d,%d,%d |ci,cj,ck - %d,%d,%d| \n",rank,i,j,k,ci,cj,ck);
              }
	    }
	  }
	}
	
        if (count > 0) {
          // As centfield is a local array with index 0 at xs,xy,xz.
          centfield[k-info.zs][j-info.ys][i-info.xs].x = sum.x / (PetscReal)count;
          centfield[k-info.zs][j-info.ys][i-info.xs].y = sum.y / (PetscReal)count;
          centfield[k-info.zs][j-info.ys][i-info.xs].z = sum.z / (PetscReal)count;
        } else {
	  //  centfield[k-zs][j-ys][i-xs].x = 0.0;
	  //  centfield[k-zs][j-ys][i-xs].y = 0.0;
	  //  centfield[k-zs][j-ys][i-xs].z = 0.0;
          LOG_ALLOW(GLOBAL, LOG_DEBUG,
                    "Rank %d: Cell (i=%d,j=%d,k=%d) got no valid corner data.\n", rank, i, j, k);
        }
	/*
        if (count > 0) {
          centfield[k][j][i].x = sum.x / (PetscReal)count;
          centfield[k][j][i].y = sum.y / (PetscReal)count;
          centfield[k][j][i].z = sum.z / (PetscReal)count;
        } else {
          centfield[k][j][i].x = 0.0;
          centfield[k][j][i].y = 0.0;
          centfield[k][j][i].z = 0.0;
        }
        */
      }
    }
  }

  LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
            "InterpolateFieldFromCornerToCenter_Vector_Interior - Rank %d completed interpolation.\n", rank);
  return 0;
}

// -----------------------------------------------------------------------------
// Interpolation: Corner -> Center (scalar)
// -----------------------------------------------------------------------------
/**
 * @brief Safely interpolate a scalar field from corner nodes (from the coordinate DM)
 *        to cell centers (from the cell-centered DM) using the provided UserCtx.
 *
 * For each cell center in the physical region of the cell-centered DM (fda), this function
 * averages the 8 surrounding corner scalar values from the coordinate DM (da). The coordinate DM
 * (da) is built on corners while the cell-centered DM (fda) covers the physical cells.
 * Index offsets are adjusted using DMDAGetLocalInfo.
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
    UserCtx *user)
{
    PetscErrorCode ierr;
    PetscMPIInt rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "InterpolateFieldFromCornerToCenter_Scalar - Rank %d starting interpolation.\n", rank);

    /* Obtain local DM information from the cell-centered DM (fda) and the coordinate DM (da). */
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(user->da,  &info);  CHKERRQ(ierr);

    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;
    PetscInt mx = info.mx, my = info.my, mz = info.mz;

    PetscInt lxs = xs, lxe = xe;
    PetscInt lys = ys, lye = ye;
    PetscInt lzs = zs, lze = ze;
  
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
  
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;


    /* Loop over each cell center. */
    for (PetscInt k = lzs; k < lze; k++) {
        for (PetscInt j = lys; j < lye; j++) {
            for (PetscInt i = lxs; i < lxe; i++) {

                PetscReal sum = 0.0;
                PetscInt count = 0;

                /* For a cell center (i,j,k) in fda, the 8 surrounding corners in da have
                   global indices ranging from (i-1,j-1,k-1) to (i,j,k). */
                for (PetscInt dk = 0; dk < 2; dk++) {
                    for (PetscInt dj = 0; dj < 2; dj++) {
                        for (PetscInt di = 0; di < 2; di++) {
                            PetscInt ci = i - 1 + di;  // global index in x for corner
                            PetscInt cj = j - 1 + dj;  // global index in y for corner
                            PetscInt ck = k - 1 + dk;  // global index in z for corner

                                sum += field[ck][cj][ci];
                                count++;   
                        }
                    }
                }

                /* Compute and store the average into the cell-centered array. */
                if (count > 0) {
                    centfield[k-zs][j-ys][i-xs] = sum / count;
                } else {
                    centfield[k-zs][j-ys][i-xs] = 0.0;
                }
            }
        }
    }

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "InterpolateFieldFromCornerToCenter_Scalar - Rank %d completed interpolation.\n", rank);
    return 0;
}


// -----------------------------------------------------------------------------
// Interpolation: Center -> Corner (scalar)
// -----------------------------------------------------------------------------
/**
 * @brief Safely interpolate a scalar field from cell centers to corners.
 *
 * For each corner at (i,j,k), we look up to 8 adjacent cell centers in
 * [i..i-1, j..j-1, k..k-1]. If a center is out of valid bounds, we skip it.
 *
 * @param[in]  centfield 3D array [mz][my][mx] of cell-center data
 * @param[out] field     3D array [mz][my][mx] for corner data
 * @param[in]  info      DMDALocalInfo from the cell-centered DM (fda).
 *
 * @return PetscErrorCode 0 on success
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Scalar(
    PetscReal     ***centfield,
    PetscReal     ***field,
    DMDALocalInfo *info)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL, LOG_DEBUG,
    "InterpolateFieldFromCenterToCorner_Scalar - Rank %d start.\n", rank);

  PetscInt mx = info->mx; 
  PetscInt my = info->my; 
  PetscInt mz = info->mz; 
  PetscInt xs = info->xs; 
  PetscInt xe = xs + info->xm;
  PetscInt ys = info->ys; 
  PetscInt ye = ys + info->ym; 
  PetscInt zs = info->zs; 
  PetscInt ze = zs + info->zm;

  // We'll define corner indices in [xs-1..xe, ys-1..ye, zs-1..ze].
  PetscInt xs_corner = xs - 1;
  PetscInt ys_corner = ys - 1;
  PetscInt zs_corner = zs - 1;
  PetscInt xe_corner = xe;
  PetscInt ye_corner = ye;
  PetscInt ze_corner = ze;

  for (PetscInt k = zs_corner; k < ze_corner; k++) {
    for (PetscInt j = ys_corner; j < ye_corner; j++) {
      for (PetscInt i = xs_corner; i < xe_corner; i++) {
        PetscReal sum   = 0.0;
        PetscInt  count = 0;

        // The adjacent cell centers are [i..i+1], [j..j+1], [k..k+1]
        for (PetscInt dk = 0; dk < 2; dk++) {
          for (PetscInt dj = 0; dj < 2; dj++) {
            for (PetscInt di = 0; di < 2; di++) {
              PetscInt ic = i + di;
              PetscInt jc = j + dj;
              PetscInt kc = k + dk;
              if ((ic >= xs && ic < xe) &&
                  (jc >= ys && jc < ye) &&
                  (kc >= zs && kc < ze))
              {
                sum += centfield[kc][jc][ic];
                count++;
              }
            }
          }
        }

        field[k][j][i] = (count > 0) ? (sum / count) : 0.0;
      }
    }
  }

  LOG_ALLOW(GLOBAL, LOG_DEBUG,
    "InterpolateFieldFromCenterToCorner_Scalar - Rank %d done.\n", rank);
  return 0;
}

// -----------------------------------------------------------------------------
// Interpolation: Center -> Corner (vector)
// -----------------------------------------------------------------------------
/**
 * @brief Safely interpolate a vector field (Cmpnts) from cell centers to corners.
 *
 * Same logic, but for vector data.
 *
 * @param[in]  centfield 3D array [mz][my][mx] of cell-centered vectors
 * @param[out] field     3D array [mz][my][mx] to store corner-based vectors
 * @param[in]  info      Local info from the cell-centered DM (fda).
 *
 * @return PetscErrorCode 0 on success
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Vector(
    Cmpnts       ***centfield,
    Cmpnts       ***field,
    DMDALocalInfo *info)
{
  PetscErrorCode ierr;
  PetscMPIInt rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL, LOG_DEBUG,
    "InterpolateFieldFromCenterToCorner_Vector - Rank %d start.\n", rank);

  PetscInt mx = info->mx; 
  PetscInt my = info->my; 
  PetscInt mz = info->mz; 
  PetscInt xs = info->xs; 
  PetscInt xe = xs + info->xm; 
  PetscInt ys = info->ys; 
  PetscInt ye = ys + info->ym; 
  PetscInt zs = info->zs; 
  PetscInt ze = zs + info->zm;

  LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG," Rank, %d, bounds: x - (%d,%d),y - (%d,%d), z - (%d,%d) \n",rank,xs,xe,ys,ye,zs,ze);

  PetscInt lxs = xs,lxe = xe;
  PetscInt lys = ys,lye = ye;
  PetscInt lzs = zs,lze = ze;
  
 
  lxs = xs+1;
  lys = ys+1;
  lzs = zs+1; 
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  

    /* Loop over each corner. */
    for (PetscInt k = lzs; k < lze; k++) {
        for (PetscInt j = lys; j < lye; j++) {
            for (PetscInt i = lxs; i < lxe; i++) {

        Cmpnts sum = {0.0, 0.0, 0.0};
        PetscInt count = 0;

        for (PetscInt dk = 0; dk < 2; dk++) {
          for (PetscInt dj = 0; dj < 2; dj++) {
            for (PetscInt di = 0; di < 2; di++) {
              PetscInt ic = i + di - 1;
              PetscInt jc = j + dj - 1;
              PetscInt kc = k + dk - 1;

              if ((ic >= xs && ic < xe) &&
                  (jc >= ys && jc < ye) &&
                  (kc >= zs && kc < ze))
              {
                sum.x += centfield[kc-zs][jc-ys][ic-xs].x;
                sum.y += centfield[kc-zs][jc-ys][ic-xs].y;
                sum.z += centfield[kc-zs][jc-ys][ic-xs].z;
                count++;

		LOG_LOOP_ALLOW(GLOBAL,LOG_DEBUG,i*j*k,500," Rank %d| i,j,k - %d,%d,%d |ci,cj,ck - %d,%d,%d| \n",rank,i,j,k,ic,jc,kc);
              }
            }
          }
        }

        if (count > 0) {
          field[k-zs][j-ys][i-xs].x = sum.x / count;
          field[k-zs][j-ys][i-xs].y = sum.y / count;
          field[k-zs][j-ys][i-xs].z = sum.z / count;
        } else {
          field[k-zs][j-ys][i-xs].x = 0.0;
          field[k-zs][j-ys][i-xs].y = 0.0;
          field[k-zs][j-ys][i-xs].z = 0.0;
          LOG_ALLOW(GLOBAL, LOG_DEBUG,
                    "Rank %d: Cell (i=%d,j=%d,k=%d) got no valid corner data.\n", rank, i, j, k);

        }
      }
    }
  }

  LOG_ALLOW(GLOBAL, LOG_DEBUG,
    "InterpolateFieldFromCenterToCorner_Vector - Rank %d done.\n", rank);
  return 0;
}

/**
 * @brief Retrieves the scalar value at the cell (iCell, jCell, kCell).
 *
 * This function does a simple piecewise (zeroth‐order) interpolation
 * by returning fieldScal[kCell][jCell][iCell], ignoring any fractional coordinates.
 *
 * @param[in]  fieldName  A string identifying the field (e.g. "temperature").
 * @param[in]  fieldScal  3D array of PetscReal (indexed as [k][j][i]).
 * @param[in]  iCell, jCell, kCell  Integral cell indices to sample.
 * @param[out] val        Pointer to a PetscReal that receives the sampled scalar.
 *
 * @return PetscErrorCode  0 on success
 */
PetscErrorCode PieceWiseLinearInterpolation_Scalar(
    const char   *fieldName,
    PetscReal  ***fieldScal,
    PetscInt      iCell,
    PetscInt      jCell,
    PetscInt      kCell,
    PetscReal    *val)
{
  PetscFunctionBegin;
  *val = fieldScal[kCell][jCell][iCell];
  
  // Optional logging
  LOG_ALLOW(LOCAL, LOG_DEBUG,
      "PieceWiseLinearInterpolation_Scalar: Field '%s' at (i=%d, j=%d, k=%d) => val=%.6f\n",
      fieldName, iCell, jCell, kCell, *val);

  PetscFunctionReturn(0);
}


/**
 * @brief Retrieves the vector (Cmpnts) at the cell (iCell, jCell, kCell).
 *
 * This function simply sets:
 *   vec->x = fieldVec[kCell][jCell][iCell].x
 *   vec->y = fieldVec[kCell][jCell][iCell].y
 *   vec->z = fieldVec[kCell][jCell][iCell].z
 * effectively a nearest‐cell or piecewise approach.
 *
 * @param[in]  fieldName  String identifying the field (e.g. "velocity").
 * @param[in]  fieldVec   3D array of Cmpnts (indexed as [k][j][i]).
 * @param[in]  iCell, jCell, kCell  Integral cell indices to sample.
 * @param[out] vec        Pointer to a Cmpnts struct that receives the sampled vector.
 *
 * @return PetscErrorCode  0 on success
 */
PetscErrorCode PieceWiseLinearInterpolation_Vector(
    const char   *fieldName,
    Cmpnts     ***fieldVec,
    PetscInt      iCell,
    PetscInt      jCell,
    PetscInt      kCell,
    Cmpnts       *vec)
{
  PetscFunctionBegin;
  vec->x = fieldVec[kCell][jCell][iCell].x;
  vec->y = fieldVec[kCell][jCell][iCell].y;
  vec->z = fieldVec[kCell][jCell][iCell].z;

  // Optional logging
  LOG_ALLOW(LOCAL, LOG_DEBUG,
      "PieceWiseLinearInterpolation_Vector: Field '%s' at (i=%d, j=%d, k=%d) => (x=%.6f, y=%.6f, z=%.6f)\n",
      fieldName, iCell, jCell, kCell, vec->x, vec->y, vec->z);

  PetscFunctionReturn(0);
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

    const PetscReal oa1 = 1.0 - a1;
    const PetscReal oa2 = 1.0 - a2;
    const PetscReal oa3 = 1.0 - a3;

    w[0] = oa1 * oa2 * oa3;  /* cornerOffsets[0] => (0,0,0) */
    w[1] = a1  * oa2 * oa3;  /* cornerOffsets[1] => (1,0,0) */
    w[2] = oa1 * a2  * oa3;  /* cornerOffsets[2] => (0,1,0) */
    w[3] = a1  * a2  * oa3;  /* cornerOffsets[3] => (1,1,0) */
    w[4] = oa1 * oa2 * a3;   /* cornerOffsets[4] => (0,0,1) */
    w[5] = a1  * oa2 * a3;   /* cornerOffsets[5] => (1,0,1) */
    w[6] = oa1 * a2  * a3;   /* cornerOffsets[6] => (0,1,1) */
    w[7] = a1  * a2  * a3;   /* cornerOffsets[7] => (1,1,1) */

    // Log the computed weights for debugging
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ComputeTrilinearWeights: Weights computed - "
        "w0=%f, w1=%f, w2=%f, w3=%f, w4=%f, w5=%f, w6=%f, w7=%f. \n",
        w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7]);
}

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
    PetscReal    *val)
{
    PetscFunctionBegin; // PETSc macro for error/stack tracing

    // Compute the 8 corner weights
    PetscReal wcorner[8];
    ComputeTrilinearWeights(a1, a2, a3, wcorner);

    // Offsets for cell corners
    PetscInt i1 = i + 1;
    PetscInt j1 = j + 1;
    PetscInt k1 = k + 1;

    // Initialize the output scalar
    PetscReal sum = 0.0;

    // Corner 0 => (i, j, k)
    sum += wcorner[0] * fieldScal[k ][j ][i ];
    // Corner 1 => (i+1, j, k)
    sum += wcorner[1] * fieldScal[k ][j ][i1];
    // Corner 2 => (i, j+1, k)
    sum += wcorner[2] * fieldScal[k ][j1][i ];
    // Corner 3 => (i+1, j+1, k)
    sum += wcorner[3] * fieldScal[k ][j1][i1];
    // Corner 4 => (i, j, k+1)
    sum += wcorner[4] * fieldScal[k1][j ][i ];
    // Corner 5 => (i+1, j, k+1)
    sum += wcorner[5] * fieldScal[k1][j ][i1];
    // Corner 6 => (i, j+1, k+1)
    sum += wcorner[6] * fieldScal[k1][j1][i ];
    // Corner 7 => (i+1, j+1, k+1)
    sum += wcorner[7] * fieldScal[k1][j1][i1];

    *val = sum;

    // Logging (optional)
    LOG_ALLOW(LOCAL, LOG_DEBUG,
        "TrilinearInterpolation_Scalar: Field '%s' at (i=%d, j=%d, k=%d), "
        "a1=%.6f, a2=%.6f, a3=%.6f -> val=%.6f.\n",
        fieldName, i, j, k, a1, a2, a3, *val);

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,
        "TrilinearInterpolation_Scalar: Completed interpolation for field '%s' across local cells.\n",
        fieldName);

    PetscFunctionReturn(0);
}


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
/**
 * @brief Computes the trilinear interpolated vector at a given point, with partial weighting if corners go out of range.
 *
 * If any of the 8 corners (i or i+1, j or j+1, k or k+1) is out of [0..mx), [0..my), [0..mz),
 * that corner is skipped and the corresponding weight is omitted. The total is normalized
 * by the sum of used weights, so we get a partial interpolation near boundaries.
 */
PetscErrorCode TrilinearInterpolation_Vector(
    const char   *fieldName,
    Cmpnts     ***fieldVec,  /* 3D array [k][j][i], dimension [mz][my][mx] */
    PetscInt      i,
    PetscInt      j,
    PetscInt      k,
    PetscReal     a1,
    PetscReal     a2,
    PetscReal     a3,
    Cmpnts       *vec)
{
  PetscFunctionBegin; // PETSc macro for error/stack tracing

  // Compute the 8 corner weights
  PetscReal wcorner[8];
  ComputeTrilinearWeights(a1, a2, a3, wcorner);

  // For partial interpolation, we'll keep track of how many corners are valid
  // and how much sum of weights is used. Then we do a final normalization.
  PetscReal sumW = 0.0;
  Cmpnts    accum = {0.0, 0.0, 0.0};

  // The eight corner indices, with their weights:
  // corners: (i,j,k), (i+1,j,k), (i,j+1,k), (i+1,j+1,k), etc.
  // We store them in an array to iterate cleanly.
  const PetscInt cornerOffsets[8][3] = {
    {0, 0, 0},
    {1, 0, 0},
    {0, 1, 0},
    {1, 1, 0},
    {0, 0, 1},
    {1, 0, 1},
    {0, 1, 1},
    {1, 1, 1}
  };

  // Weighted partial sum
  for (int c = 0; c < 8; c++) {
    const PetscInt di = cornerOffsets[c][0];
    const PetscInt dj = cornerOffsets[c][1];
    const PetscInt dk = cornerOffsets[c][2];
    PetscInt iC = i + di;
    PetscInt jC = j + dj;
    PetscInt kC = k + dk;

    /*
    // skip if out of domain
    // (Assuming you know global domain is [0..mx), [0..my), [0..mz).)
    if (iC < 0 || iC >= (PetscInt)userGlobalMx ||
        jC < 0 || jC >= (PetscInt)userGlobalMy ||
        kC < 0 || kC >= (PetscInt)userGlobalMz)
    {
      // skip this corner
      continue;
    }

    */

    // Otherwise, accumulate
    accum.x += wcorner[c] * fieldVec[kC][jC][iC].x;
    accum.y += wcorner[c] * fieldVec[kC][jC][iC].y;
    accum.z += wcorner[c] * fieldVec[kC][jC][iC].z;
    sumW += wcorner[c];
  }

  // If sumW=0 => out-of-range or zero weighting => set (0,0,0)
  if (sumW > 1.0e-14) {
    vec->x = accum.x / sumW;
    vec->y = accum.y / sumW;
    vec->z = accum.z / sumW;
  } else {
    vec->x = 0.0;  vec->y = 0.0;  vec->z = 0.0;
  }

  PetscFunctionReturn(0);
}


/**
 * @brief Interpolates a single field (scalar or vector) for one particle, storing the result in a swarm array.
 *
 * This helper function is used inside a larger loop over local particles.
 *
 * @param[in]  fieldName   A string identifying the field (e.g. "velocity", "temperature", etc.).
 * @param[in]  fieldPtr    A pointer to the local DMDA array for the field:
 *                         - (PetscReal ***) if scalar (blockSize = 1)
 *                         - (Cmpnts    ***) if vector (blockSize = 3)
 * @param[in]  iCell, jCell, kCell  The cell indices for this particle (already clamped if needed).
 * @param[in]  a1, a2, a3   Interpolation coefficients in [0..1].
 *                          If using PiecewiseLinearInterpolation, these are currently ignored.
 * @param[out] swarmOut     A pointer to the DMSwarm output array:
 *                            - (PetscReal*) if scalar dof=1
 *                            - (PetscReal*) if vector dof=3 (storing x,y,z in consecutive reals)
 *                            - or (Cmpnts*) if you store the result as a Cmpnts struct
 * @param[in]  p            The local particle index (used to compute the correct offset in swarmOut).
 * @param[in]  blockSize    The number of degrees of freedom (1 => scalar, 3 => vector).
 *
 * This routine demonstrates a switch between:
 *  - PiecewiseLinearInterpolation (zeroth order / nearest cell)
 *  - TrilinearInterpolation (8-corner weighted).
 *
 * By default, PiecewiseLinearInterpolation is active, while the TrilinearInterpolation calls are commented out.
 * To switch to trilinear, simply comment/uncomment appropriately.
 *
 * Logging calls are provided (LOG_ALLOW) that you can adapt to your existing logging system.
 *
 * @return PetscErrorCode 0 on success, non-zero on error.
 */
static inline PetscErrorCode InterpolateEulerFieldToSwarmForParticle(
    const char  *fieldName,
    void        *fieldPtr,   /* typed pointer => either (PetscReal***) or (Cmpnts***) */
    PetscInt     iCell,
    PetscInt     jCell,
    PetscInt     kCell,
    PetscReal    a1,
    PetscReal    a2,
    PetscReal    a3,
    void        *swarmOut,   /* typed pointer => (PetscReal*) or (Cmpnts*) or dof=3 array */
    PetscInt     p,          /* particle index */
    PetscInt     blockSize)  /* dof=1 => scalar, dof=3 => vector */
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // Optional logging at start
  LOG_ALLOW(LOCAL, LOG_DEBUG,
    "InterpolateEulerFieldToSwarmForParticle: field='%s', blockSize=%d, "
    "cell IDs=(%d,%d,%d), weights=(%.4f,%.4f,%.4f)\n",
    fieldName, blockSize, iCell, jCell, kCell, a1, a2, a3);

  /*
     If blockSize=1, we interpret the fieldPtr as a 3D array of PetscReal (scalar).
     If blockSize=3, we interpret the fieldPtr as a 3D array of Cmpnts (vector).
   */
  if (blockSize == 1) {
    /* Scalar field: Cast fieldPtr to (PetscReal ***). */
    PetscReal ***fieldScal = (PetscReal ***) fieldPtr;
    PetscReal val;

    // Currently using trilinear.
    ierr = TrilinearInterpolation(fieldName, fieldScal,
                                   iCell, jCell, kCell,
                                   a1, a2, a3,
                                   &val);
     CHKERRQ(ierr);

    // Alternative (commented) call to PiecewiseLinearInterpolation (zeroth order) :
    //   PetscErrorCode ierr = PiecewiseLinearInterpolation(fieldName,
    //                                                   fieldScal,
    //                                                   iCell, jCell, kCell,
    //                                                   &val);
    // CHKERRQ(ierr);

    // Write the scalar result to the swarm output at index [p].
    ((PetscReal*)swarmOut)[p] = val;

    LOG_ALLOW(LOCAL, LOG_DEBUG,
      "InterpolateEulerFieldToSwarmForParticle [Scalar]: field='%s', result=%.6f "
      "stored at swarmOut index p=%d.\n", fieldName, val, (int)p);
  }
  else if (blockSize == 3) {
    /* Vector field: Cast fieldPtr to (Cmpnts ***). */
    Cmpnts ***fieldVec = (Cmpnts ***) fieldPtr;
    Cmpnts vec;

    // Piecewise interpolation (zeroth order).
    //  PetscErrorCode ierr = PieceWiseLinearInterpolation(fieldName,
    //                                                   fieldVec,
    //                                                   iCell, jCell, kCell,
    //                                                   &vec);
    // CHKERRQ(ierr);

    // Alternative (commented) call to trilinear:
     ierr = TrilinearInterpolation(fieldName, fieldVec,
                                   iCell, jCell, kCell,
                                   a1, a2, a3,
                                   &vec);
     CHKERRQ(ierr);

    // If swarmOut is an array of 3 reals per particle:
    ((PetscReal*)swarmOut)[3*p + 0] = vec.x;
    ((PetscReal*)swarmOut)[3*p + 1] = vec.y;
    ((PetscReal*)swarmOut)[3*p + 2] = vec.z;

    LOG_ALLOW(LOCAL, LOG_DEBUG,
      "InterpolateEulerFieldToSwarmForParticle [Vector]: field='%s', result=(%.6f,%.6f,%.6f) "
      "stored at swarmOut[3p..3p+2], p=%d.\n",
      fieldName, vec.x, vec.y, vec.z, (int)p);

    /*
       If you store the vector result as a Cmpnts in the swarm, do instead:
         ((Cmpnts*)swarmOut)[p] = vec;
       but ensure your DMSwarm field is sized for a Cmpnts struct.
    */
  }
  else {
    /* If blockSize isn't 1 or 3, we raise an error. */
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP,
      "InterpolateEulerFieldToSwarmForParticle: only blockSize=1 or 3 supported, got %d.",
      (int)blockSize);
  }

  PetscFunctionReturn(0);
}


/**
 * @brief Interpolates a cell-centered field (scalar or vector) onto DMSwarm particles,
 *        converting the cell-center data to corner data first, then looping over particles.
 *
 * Steps:
 *   1) Check that the Vec has blockSize=1 or 3 (scalar vs. vector).
 *   2) Map the cell-centered Vec to a local array (fieldGlobal -> localPtr).
 *   3) Allocate a corner array (cornerPtr) via Allocate3DArray(...), sized (zm+1, ym+1, xm+1).
 *   4) Convert from cell-centers to corners via InterpolateFieldFromCenterToCorner(...).
 *   5) Restore the cell-centered array.
 *   6) Retrieve DMSwarm fields: "DMSwarm_CellID", "weight", and swarmOutFieldName.
 *   7) Loop over local particles, clamp i/j/k, skip or zero out if out of range, read (a1,a2,a3).
 *   8) Call InterpolateEulerFieldToSwarmForParticle(...) with cornerPtr to do final interpolation.
 *   9) Restore swarm fields, free the corner array.
 *
 * @param[in]  user              User context with:
 *                                - user->da     (cell-centered DMDA),
 *                                - user->swarm  (DMSwarm).
 * @param[in]  fieldGlobal       Vec with blockSize=1 or 3, storing the cell-centered field.
 * @param[in]  fieldName         Human-readable field name for logging (e.g. "velocity").
 * @param[in]  swarmOutFieldName Name of the DMSwarm field where interpolation results go.
 *
 * @return PetscErrorCode  0 on success, non-zero on error.
 */
PetscErrorCode InterpolateEulerFieldToSwarm(
    UserCtx    *user,
    Vec         fieldGlobal,       /* DMDA Vec containing cell‐center data */
    const char *fieldName,         /* e.g., "Ucat" */
    const char *swarmOutFieldName) /* Name of the output DMSwarm field */
{
  PetscErrorCode ierr;
  DM             fda    = user->fda;      /* DM for cell‐center field data */
  DM             da     = user->da;       /* DM for grid information (local indices) */
  DM             swarm  = user->swarm;    /* DMSwarm for particles */
  PetscInt       bs;                  /* Block size: 1 (scalar) or 3 (vector) */
  DMDALocalInfo  info;                /* Local grid info */
  void          *localPtr   = NULL;     /* Pointer to cell‐center data from fda */
  void          *cornerPtr  = NULL;     /* Will hold the typed corner array */
  void          *swarmOut   = NULL;     /* Pointer to the swarm output field */
  PetscInt64    *cellIDs    = NULL;     /* Particle cell indices from swarm */
  PetscReal     *weights    = NULL;     /* Interpolation coefficients from swarm */
  PetscInt       nLocal;

  PetscFunctionBegin;

  /* (A) Check block size and get local domain info */
  ierr = VecGetBlockSize(fieldGlobal, &bs); CHKERRQ(ierr);
  if (bs != 1 && bs != 3) {
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP,
             "InterpolateEulerFieldToSwarm: blockSize must be 1 or 3, got %d.", (int)bs);
  }
  ierr = DMDAGetLocalInfo(fda, &info); CHKERRQ(ierr);
  LOG_ALLOW(GLOBAL, LOG_INFO,
    "InterpolateEulerFieldToSwarm: Starting with field='%s', blockSize=%d, local domain: (%d x %d x %d)\n",
    fieldName, bs, info.xm, info.ym, info.zm);

  /* (B) Map the cell-centered Vec to a local array using the DM attached to fieldGlobal */
  ierr = DMDAVecGetArray(fda, fieldGlobal, &localPtr); CHKERRQ(ierr);

    // Corner domain is larger than cell center domain 
    PetscInt nz = info.zm + 1; 
    PetscInt ny = info.ym + 1;
    PetscInt nx = info.xm + 1;
    LOG_ALLOW(LOCAL, LOG_DEBUG,
      "InterpolateEulerFieldToSwarm: Allocating corner array of size (%d x %d x %d), blockSize = %d\n",
      nx, ny, nz, bs);
    if (bs == 1) {
      /* Declare a typed pointer for scalar corners */
      PetscReal ***cornerScal = NULL;
      ierr = Allocate3DArray(&cornerScal, nz, ny, nx); CHKERRQ(ierr);
      if (!cornerScal) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "InterpolateEulerFieldToSwarm: 3D Array Allocation for scalar corners failed.\n");
      }
      /* Save typed pointer into cornerPtr so later code can cast appropriately */
      cornerPtr = (void*) cornerScal;
      /* (D) Convert cell-center data to corners for scalar field */
      ierr = InterpolateFieldFromCenterToCorner( (PetscReal ***) localPtr, cornerScal, &info); CHKERRQ(ierr);
    } else {
      /* For vector fields */
      Cmpnts ***cornerVec = NULL;
      ierr = Allocate3DArray(&cornerVec, nz, ny, nx); CHKERRQ(ierr);
      if (!cornerVec) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "InterpolateEulerFieldToSwarm: 3D Array Allocation for vector corners failed.\n");
      }
      cornerPtr = (void*) cornerVec;
      ierr = InterpolateFieldFromCenterToCorner( (Cmpnts ***) localPtr, cornerVec, &info); CHKERRQ(ierr);
    }
    LOG_ALLOW(LOCAL, LOG_INFO,
      "InterpolateEulerFieldToSwarm: Completed center-to-corner interpolation for field='%s'.\n",
      fieldName);

  /* (E) Restore the cell-centered array since we now have corner data */
  ierr = DMDAVecRestoreArray(fda, fieldGlobal, &localPtr); CHKERRQ(ierr);

  /* (F) Retrieve swarm fields: cell IDs, weights, and the output field */
  ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "DMSwarm_CellID",  NULL, NULL, (void**)&cellIDs);  CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "weight",          NULL, NULL, (void**)&weights);  CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, swarmOutFieldName, NULL, NULL, &swarmOut);         CHKERRQ(ierr);
  LOG_ALLOW(GLOBAL, LOG_INFO,
    "InterpolateEulerFieldToSwarm: Found %d local particles to process for field='%s'.\n",
    nLocal, fieldName);

  /* (G) Loop over each local particle and perform the final interpolation from corners */
  for (PetscInt p = 0; p < nLocal; p++) {
    PetscInt iCell = (PetscInt)cellIDs[3*p + 0];
    PetscInt jCell = (PetscInt)cellIDs[3*p + 1];
    PetscInt kCell = (PetscInt)cellIDs[3*p + 2];

    /* Boundary clamp: adjust indices to be within [0, mx), [0, my), [0, mz) */
    if (iCell >= user->info.mx) iCell = user->info.mx - 1;
    if (jCell >= user->info.my) jCell = user->info.my - 1;
    if (kCell >= user->info.mz) kCell = user->info.mz - 1;
    if (iCell < 0 || jCell < 0 || kCell < 0 ||
        iCell >= user->info.mx ||
        jCell >= user->info.my ||
        kCell >= user->info.mz)
    {
      /* Out-of-range: set output to zero */
      if (bs == 1) {
        ((PetscReal*)swarmOut)[p] = 0.0;
      } else {
        ((PetscReal*)swarmOut)[3*p + 0] = 0.0;
        ((PetscReal*)swarmOut)[3*p + 1] = 0.0;
        ((PetscReal*)swarmOut)[3*p + 2] = 0.0;
      }
      continue;
    }

    /* Retrieve interpolation coefficients (a1, a2, a3) for this particle */
    PetscReal alpha1 = weights[3*p + 0];
    PetscReal alpha2 = weights[3*p + 1];
    PetscReal alpha3 = weights[3*p + 2];

    /* (Optional) If your final interpolation expects a corner offset, adjust here.
       For example, if the cell center corresponds to the average of corners at (iCell,jCell,kCell)
       and (iCell+1, jCell+1, kCell+1), you might add +1. For now, we pass indices as is.
    */
    PetscInt iUse = iCell;  // + 1 if required
    PetscInt jUse = jCell;  // + 1 if required
    PetscInt kUse = kCell;  // + 1 if required

    /* (H) Call the per-particle interpolation function.
       This function will use your _Generic macro (TrilinearInterpolation or PiecewiseLinearInterpolation)
       on the corner data.
    */
    ierr = InterpolateEulerFieldToSwarmForParticle(
              fieldName,    /* e.g., "Ucat" */
              cornerPtr,    /* typed pointer: (PetscReal***) or (Cmpnts***) */
              iUse, jUse, kUse,
              alpha1, alpha2, alpha3,
              swarmOut,     /* pointer to swarm output array */
              p,            /* particle index */
              bs);          /* block size: 1 or 3 */
    CHKERRQ(ierr);
  }

  /* (I) Restore swarm fields */
  ierr = DMSwarmRestoreField(swarm, swarmOutFieldName, NULL, NULL, &swarmOut);        CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "weight",          NULL, NULL, (void**)&weights); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID",  NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);

  /* (J) Deallocate the corner array using the generic deallocation macro */
  if (bs == 1) {
    PetscReal ***cornerScal = (PetscReal ***) cornerPtr;
    ierr = Deallocate3DArray(cornerScal, info.zm, info.ym); CHKERRQ(ierr);
  } else {
    Cmpnts ***cornerVec = (Cmpnts ***) cornerPtr;
    ierr = Deallocate3DArray(cornerVec, info.zm, info.ym); CHKERRQ(ierr);
  }

  LOG_ALLOW(GLOBAL, LOG_INFO,
    "InterpolateEulerFieldToSwarm: Completed interpolation of field='%s' for %d local particles.\n",
    fieldName, nLocal);

  PetscFunctionReturn(0);
}

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
PetscErrorCode InterpolateAllFieldsToSwarm(UserCtx *user)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /* 
     1) Interpolate the 'velocity' field (user->Ucat) into DMSwarm's "swarmVelocity".
        - The function InterpolateOneFieldOverSwarm() loops over all local particles,
          retrieves iCell, jCell, kCell, (a1,a2,a3) from the DMSwarm,
          and calls the appropriate trilinear interpolation routine.

        - "velocity" => just a label used in logging or debugging.
        - "swarmVelocity" => the DMSwarm field name where we store the result
          (assume dof=3 if user->Ucat has blockSize=3).
  */

  LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
    "interpolateAllFieldsToSwarm: Interpolation of ucat to velocity begins.\n");
  ierr = InterpolateEulerFieldToSwarm(user, user->Ucat, 
                                      "Ucat", 
                                      "velocity"); CHKERRQ(ierr);

  /* 
     2) (OPTIONAL) If you have more fields, you can interpolate them similarly.

     For example, if user->Tcat is a scalar Vec for "temperature" and the DMSwarm
     has a field "Temperature" (with dof=1), you could do:

         ierr = InterpolateOneFieldOverSwarm(user, user->Tcat,
                                             "temperature",
                                             "swarmTemperature"); CHKERRQ(ierr);

     For pressure:
         ierr = InterpolateOneFieldOverSwarm(user, user->Pcat,
                                             "pressure",
                                             "swarmPressure"); CHKERRQ(ierr);
     
     ...and so forth.
   */

  /* 
     3) Optionally, synchronize or log that all fields are done 
  */
  ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);
  LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,
    "interpolateAllFieldsToSwarm: Completed interpolating all fields to the swarm.\n");

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
	if (i >= user->info.mx) i = user->info.mx - 1;
	if (j >= user->info.my) j = user->info.my - 1;
	if (k >= user->info.mz) k = user->info.mz - 1;

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
