
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
 * @param[in]     fieldName PoPetscInter to a string representing the field name (for logging purposes).
 * @param[in,out] vecField  PoPetscInter to the `Cmpnts` structure where the computed vector field will be stored.
 * @param[in]     coor      PoPetscInter to the `Cmpnts` structure containing the input coordinate values.
 *
 * @return PetscErrorCode Returns 0 on successful execution, non-zero on failure.
 */
PetscErrorCode SetLocalCartesianField_Vector(const char *fieldName, Cmpnts *vecField, Cmpnts *coor,PetscInt FieldInitialization,PetscReal Constant) {
    // Log input coordinate values for debugging purposes.
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
        "SetLocalCartesianField_Vector: Input Coordinates - x: %f, y: %f, z: %f\n",
        coor->x, coor->y, coor->z);

    // Compute vector components as the sine of the coordinate values.
    if(FieldInitialization == 1){
      vecField->x = sin(coor->x);
      vecField->y = sin(coor->y);
      vecField->z = sin(coor->z);
    }
    else if(FieldInitialization == 2){
      vecField->x = cos(coor->x)*sin(coor->y)*sin(coor->z);
      vecField->y = sin(coor->x)*cos(coor->y)*sin(coor->z);
      vecField->z = sin(coor->x)*sin(coor->y)*cos(coor->z);
    }
    else if(FieldInitialization == 0){
      vecField->x = 0.0;
      vecField->y = 0.0;
      vecField->z = Constant;
    }
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
 * @param[in]     fieldName   PoPetscInter to a string representing the field name (for logging purposes).
 * @param[in,out] scalarField PoPetscInter to the PetscReal where the computed scalar field value will be stored.
 * @param[in]     coor        PoPetscInter to the `Cmpnts` structure containing the input coordinate values.
 *
 * @return PetscErrorCode Returns 0 on successful execution, non-zero on failure.
 */
    PetscErrorCode SetLocalCartesianField_Scalar(const char *fieldName, PetscReal *scalarField, Cmpnts *coor,PetscInt FieldInitialization, PetscReal Constant) {
    // Log input coordinate values for debugging purposes.
    LOG_ALLOW(LOCAL, LOG_DEBUG,
        "SetLocalCartesianField_Scalar: Input Coordinates - x: %f, y: %f, z: %f\n",
        coor->x, coor->y, coor->z);

    // Compute scalar field value as the sum of the sine functions of the coordinates.
    if(FieldInitialization == 1){
    *scalarField = sin(coor->x) + sin(coor->y) + sin(coor->z);
    }
    else if(FieldInitialization == 0){
      *scalarField = Constant;
	}
    // Log computed scalar field value along with the field name.
    LOG_ALLOW(LOCAL, LOG_DEBUG,
        "SetLocalCartesianField_Scalar: Computed Scalar %s - value: %f\n",
        fieldName, *scalarField);

    return 0;
}


/**
 * @brief Sets an analytical Cartesian field (scalar or vector) based on a field name.
 *
 * This function calculates analytical values for specified fields.
 * For scalar fields (P, Nvert on user->da), values are computed at cell centers and
 * stored at the DMDA index corresponding to the cell's origin node.
 * For vector fields (Ucat, Ucont on user->fda), values representing the center of
 * cell C(i,j,k) (origin N(i,j,k)) are stored at the NODE index [k][j][i] of user->fda,
 * but ONLY for globally interior nodes. Boundary values are typically set by BC routines.
 *
 * @param[in]  user      Pointer to the UserCtx structure.
 * @param[in]  fieldName Name of the field to update ("Ucat", "Ucont", "P", or "Nvert").
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

    /* --- 1. Identify Target Vec and Type --- */
    Vec fieldVec = NULL;
    DM  targetDM = NULL; // DM associated with the target field's Vec
    DM  accessDM = NULL; // DM used for DMDAVecGet/RestoreArray (usually same as targetDM)
    DM  infoDM_for_cell_loop = NULL; // DM whose DMDALocalInfo will define the primary cell/node loop
    PetscInt fieldIsVector = -1;
    PetscInt expected_bs;
    PetscReal ConstantValue;

    if (strcmp(fieldName, "Ucat") == 0) {
        fieldVec = user->Ucat;
        targetDM = user->fda;
        accessDM = user->fda;
        infoDM_for_cell_loop = user->fda; // Loop over nodes of fda, treat as cell origins
        fieldIsVector = 1;
        expected_bs = 3;
        ConstantValue = user->ConstantVelocity;
    } else if (strcmp(fieldName, "Ucont") == 0) {
        fieldVec = user->Ucont;
        targetDM = user->fda;
        accessDM = user->fda;
        infoDM_for_cell_loop = user->fda; // Loop over nodes of fda, treat as cell origins
        fieldIsVector = 1;
        expected_bs = 3;
        ConstantValue = user->ConstantContra;
    } else if (strcmp(fieldName, "P") == 0) {
        fieldVec = user->P;
        targetDM = user->da;
        accessDM = user->da;
        infoDM_for_cell_loop = user->da;  // Loop over "cells" defined by da's node origins
        fieldIsVector = 0;
        expected_bs = 1;
        ConstantValue = user->ConstantPressure;
    } else if (strcmp(fieldName, "Nvert") == 0) {
        fieldVec = user->Nvert;
        targetDM = user->da;
        accessDM = user->da;
        infoDM_for_cell_loop = user->da; // Loop over "cells" defined by da's node origins
        fieldIsVector = 0;
        expected_bs = 1;
        ConstantValue = user->ConstantNvert;
    } else {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Field '%s' not supported", fieldName);
    }

    /* --- 2. Verify Target DM Block Size --- */
    PetscInt bs;
    ierr = DMGetBlockSize(targetDM, &bs); CHKERRQ(ierr);
    if (bs != expected_bs) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                 "Expected block size %d for field '%s' based on its DM, got %d", expected_bs, fieldName, bs);
    }

    /* --- 3. Get Coordinates and Prepare for Cell Center Calculation --- */
    Vec localCoor;
    ierr = DMGetCoordinatesLocal(user->da, &localCoor); CHKERRQ(ierr); // Coords are on da, layout fda

    Vec globalCoor;
    ierr = DMGetCoordinates(user->da, &globalCoor); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->fda, globalCoor, INSERT_VALUES, localCoor); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->fda, globalCoor, INSERT_VALUES, localCoor); CHKERRQ(ierr);

    Cmpnts ***coor_arr; // Node coordinates, layout by fda
    ierr = DMDAVecGetArrayRead(user->fda, localCoor, &coor_arr); CHKERRQ(ierr);

    // Get DMDALocalInfo for the DM that defines our primary loop (cell origins)
    DMDALocalInfo info_loop_dm;
    ierr = DMDAGetLocalInfo(infoDM_for_cell_loop, &info_loop_dm); CHKERRQ(ierr);

    // Determine owned CELL ranges based on the info_loop_dm
    // These are cells whose origin node is owned by info_loop_dm's partitioning
    PetscInt xs_cell_global_i, xm_cell_local_i, xe_cell_global_i_excl;
    PetscInt ys_cell_global_j, ym_cell_local_j, ye_cell_global_j_excl;
    PetscInt zs_cell_global_k, zm_cell_local_k, ze_cell_global_k_excl;

    
    ierr = GetOwnedCellRange(&info_loop_dm, 0, &xs_cell_global_i, &xm_cell_local_i); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_loop_dm, 1, &ys_cell_global_j, &ym_cell_local_j); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_loop_dm, 2, &zs_cell_global_k, &zm_cell_local_k); CHKERRQ(ierr);

    xe_cell_global_i_excl = xs_cell_global_i + xm_cell_local_i;
    ye_cell_global_j_excl = ys_cell_global_j + ym_cell_local_j;
    ze_cell_global_k_excl = zs_cell_global_k + zm_cell_local_k;

    // Allocate temporary local array for cell-center coordinates
    Cmpnts ***centcoor = NULL; // Will be indexed [k_local][j_local][i_local]
    if (xm_cell_local_i > 0 && ym_cell_local_j > 0 && zm_cell_local_k > 0) { // Only allocate if there are cells
        ierr = Allocate3DArray(&centcoor, zm_cell_local_k, ym_cell_local_j, xm_cell_local_i); CHKERRQ(ierr);
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Allocated centcoor[%d][%d][%d] for field %s.\n", rank, zm_cell_local_k, ym_cell_local_j, xm_cell_local_i, fieldName);

        // Interpolate node coordinates to cell centers
        // This function needs to be aware that coor_arr is global-indexed,
        // and centcoor is local-indexed. It needs info_loop_dm to map.
        // Or, more simply, InterpolateFieldFromCornerToCenter itself loops over owned cells
        // using GetOwnedCellRange internally or takes xs_cell_global, xm_cell_local as input.
        // Assuming the provided InterpolateFieldFromCornerToCenter handles this correctly
        // based on its own GetOwnedCellRange call (using user->fda for info_nodes).
        // Let's pass the necessary info if Interpolate is to be generic.
        // For now, we assume InterpolateFieldFromCornerToCenter_Vector (called by generic)
        // has been fixed to use the corrected GetOwnedCellRange with user->fda info.
        ierr = InterpolateFieldFromCornerToCenter(coor_arr, centcoor, user); CHKERRQ(ierr);
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Interpolated coordinates to cell centers for field %s.\n", rank, fieldName);
    }


    /* --- 4. Populate Target Field Vec --- */
    if (fieldIsVector) { // Handle Ucat or Ucont
        Cmpnts ***vecField_arr;
        ierr = DMDAVecGetArray(accessDM, fieldVec, &vecField_arr); CHKERRQ(ierr);

        // Loop over the GLOBAL indices of the cell ORIGIN NODES owned by this rank
        // (as determined by info_loop_dm and GetOwnedCellRange)
        for (PetscInt k_glob_cell_origin = zs_cell_global_k; k_glob_cell_origin < ze_cell_global_k_excl; k_glob_cell_origin++) {
            for (PetscInt j_glob_cell_origin = ys_cell_global_j; j_glob_cell_origin < ye_cell_global_j_excl; j_glob_cell_origin++) {
                for (PetscInt i_glob_cell_origin = xs_cell_global_i; i_glob_cell_origin < xe_cell_global_i_excl; i_glob_cell_origin++) {

                    PetscInt k_local = k_glob_cell_origin - zs_cell_global_k;
                    PetscInt j_local = j_glob_cell_origin - ys_cell_global_j;
                    PetscInt i_local = i_glob_cell_origin - xs_cell_global_i;

                    // Defensive check for centcoor access (should be guaranteed by loop bounds if allocation was correct)
                    if (!centcoor || k_local < 0 || k_local >= zm_cell_local_k ||
                        j_local < 0 || j_local >= ym_cell_local_j ||
                        i_local < 0 || i_local >= xm_cell_local_i) {
                        // This should not happen if xm_cell_local > 0 checks were done before allocation/loop
                        LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: Skipping cell origin (%d,%d,%d) due to invalid local index for centcoor for field %s.",
                                  rank, i_glob_cell_origin, j_glob_cell_origin, k_glob_cell_origin, fieldName);
                        continue;
                    }
                    Cmpnts* cell_center_coord_ptr = &centcoor[k_local][j_local][i_local];

                    // The value for cell C(i_glob_cell_origin, ...) is stored at NODE index [i_glob_cell_origin][...]
                    PetscInt node_k_target = k_glob_cell_origin;
                    PetscInt node_j_target = j_glob_cell_origin;
                    PetscInt node_i_target = i_glob_cell_origin;

                    // Correct interior check: GLOBAL node index vs GLOBAL domain cell counts
                    // A node N(i,j,k) is interior if 0 < i < IM, 0 < j < JM, 0 < k < KM
                    // (where IM, JM, KM are global CELL counts, so nodes are 0..IM, 0..JM, 0..KM)
                    // Node i is interior if global_node_i > 0 AND global_node_i < (user->IM)
                    // (i.e., not N0 and not N_IM)
                    PetscBool is_globally_interior_node =
                        (node_i_target > 0 && node_i_target < user->IM &&  // Not node 0 or node IM (max cell index is IM-1)
                         node_j_target > 0 && node_j_target < user->JM &&
                         node_k_target > 0 && node_k_target < user->KM);

                    if (is_globally_interior_node) {
                        // vecField_arr expects global node indices
                        ierr = SetLocalCartesianField(fieldName, &vecField_arr[node_k_target][node_j_target][node_i_target],
                                                      cell_center_coord_ptr, user->FieldInitialization, ConstantValue); CHKERRQ(ierr);
                    }
                }
            }
        }
        ierr = DMDAVecRestoreArray(accessDM, fieldVec, &vecField_arr); CHKERRQ(ierr);
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Populated interior node indices of %s.\n", rank, fieldName);

    } else { // Handle P or Nvert (Scalar Fields)
        PetscReal ***scalarField_arr;
        ierr = DMDAVecGetArray(accessDM, fieldVec, &scalarField_arr); CHKERRQ(ierr);

        // Loop over the GLOBAL indices of the cell ORIGIN NODES owned by this rank
        for (PetscInt k_glob_cell_origin = zs_cell_global_k; k_glob_cell_origin < ze_cell_global_k_excl; k_glob_cell_origin++) {
            for (PetscInt j_glob_cell_origin = ys_cell_global_j; j_glob_cell_origin < ye_cell_global_j_excl; j_glob_cell_origin++) {
                for (PetscInt i_glob_cell_origin = xs_cell_global_i; i_glob_cell_origin < xe_cell_global_i_excl; i_glob_cell_origin++) {

                    PetscInt k_local = k_glob_cell_origin - zs_cell_global_k;
                    PetscInt j_local = j_glob_cell_origin - ys_cell_global_j;
                    PetscInt i_local = i_glob_cell_origin - xs_cell_global_i;

                    if (!centcoor || k_local < 0 || k_local >= zm_cell_local_k ||
                        j_local < 0 || j_local >= ym_cell_local_j ||
                        i_local < 0 || i_local >= xm_cell_local_i) {
                        LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: Skipping cell origin (%d,%d,%d) due to invalid local index for centcoor for field %s.",
                                  rank, i_glob_cell_origin, j_glob_cell_origin, k_glob_cell_origin, fieldName);
                        continue;
                    }
                    Cmpnts* cell_center_coord_ptr = &centcoor[k_local][j_local][i_local];

                    // scalarField_arr (from user->da) expects global indices corresponding to cell origin nodes
                    ierr = SetLocalCartesianField(fieldName, &scalarField_arr[k_glob_cell_origin][j_glob_cell_origin][i_glob_cell_origin],
                                                  cell_center_coord_ptr, user->FieldInitialization, ConstantValue); CHKERRQ(ierr);
                }
            }
        }
        ierr = DMDAVecRestoreArray(accessDM, fieldVec, &scalarField_arr); CHKERRQ(ierr);
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Populated scalar field %s.\n", rank, fieldName);
    }

    /* --- 5. Cleanup --- */
    if (centcoor) { // Deallocate only if it was allocated
        ierr = Deallocate3DArray(centcoor, zm_cell_local_k, ym_cell_local_j); CHKERRQ(ierr);
    }
    ierr = DMDAVecRestoreArrayRead(user->fda, localCoor, &coor_arr); CHKERRQ(ierr);

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "SetAnalyticalCartesianField: Completed for field '%s' on rank %d.\n", fieldName, rank);
    MPI_Barrier(PETSC_COMM_WORLD);
    return 0;
}
  

// Forward declaration if needed, or place ApplyAnalyticalBC_Vector before ApplyAnalyticalBC
static PetscErrorCode ApplyAnalyticalBC_Vector(UserCtx *user, const char *fieldName, DM vectorDM, Vec coordVec, Vec fieldVec);

// Static function to handle vector fields like Ucat/Ucont
static PetscErrorCode ApplyAnalyticalBC_Vector(
    UserCtx    *user,
    const char *fieldName,
    DM          vectorDM,   // e.g., user->fda
    Vec         coordVec,   // e.g., user->coordinates
    Vec         fieldVec    // e.g., user->Ucat or user->Ucont
)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    Cmpnts       ***field_arr;
    Vec            localCoords;
    Cmpnts       ***coor_arr;
    DMDALocalInfo  info;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    // No top-level log here, the dispatcher will log

    // Get Local Info & Global Dimensions from the specific vectorDM
    ierr = DMDAGetLocalInfo(vectorDM, &info); CHKERRQ(ierr);
    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;
    PetscInt mx = info.mx;
    PetscInt my = info.my;
    PetscInt mz = info.mz;

    // Get Access to Vectors
    ierr = DMDAVecGetArray(vectorDM, fieldVec, &field_arr); CHKERRQ(ierr);

    /*
    // Get Coordinates (using the same DM as the field for consistency)
    ierr = DMGetLocalVector(vectorDM, &localCoords); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(vectorDM, coordVec, INSERT_VALUES, localCoords); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(vectorDM, coordVec, INSERT_VALUES, localCoords); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(vectorDM, localCoords, &coor_arr); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Obtained coord array pointer %p\n", rank, (void*)coor_arr);
    */

    DMDAVecGetArrayRead(vectorDM,coordVec,&coor_arr);

    /*
    // --- BEGIN DEBUG: Test Coordinate Ghost Reads ---
    if (!coor_arr) {
      LOG_ALLOW(LOCAL, LOG_ERROR, "Rank %d: ERROR - coor_arr pointer is NULL after GetArray!\n", rank);
    } else {
      // Test reading a point in the ghost region (if applicable)
      // Example: Read the point corresponding to the global index one layer below owned start (zs-1)
      // Requires stencil width s >= 1
      if (info.zs > 0 && info.dim >= 3) { // Check if there's a layer below
        PetscInt k_ghost_global = info.zs - 1;
        PetscInt j_ghost_global = info.ys; // Just use the start of owned y
        PetscInt i_ghost_global = info.xs; // Just use the start of owned x

        // Convert global ghost index to local index relative to GHOSTED array start
        PetscInt k_ghost_local = k_ghost_global - info.gzs;
        PetscInt j_ghost_local = j_ghost_global - info.gys;
        PetscInt i_ghost_local = i_ghost_global - info.gxs;

        // Check if calculated local index is within the ghosted array bounds
        if (k_ghost_local >= 0 && k_ghost_local < info.gzm &&
            j_ghost_local >= 0 && j_ghost_local < info.gym &&
            i_ghost_local >= 0 && i_ghost_local < info.gxm)
	  {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: DEBUG Attempting coord ghost read at GLOBAL (%d,%d,%d) -> LOCAL (%d,%d,%d)\n",
                      rank, k_ghost_global, j_ghost_global, i_ghost_global, k_ghost_local, j_ghost_local, i_ghost_local);
            Cmpnts test_coord_ghost = coor_arr[k_ghost_local][j_ghost_local][i_ghost_local];
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: DEBUG SUCCESS reading coord ghost: x=%g\n", rank, test_coord_ghost.x);
	  } else {
	  LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: DEBUG Calculated coord ghost local index (%d,%d,%d) is out of bounds [%d,%d,%d]. Skipping read.\n",
		    rank, k_ghost_local, j_ghost_local, i_ghost_local, info.gzm, info.gym, info.gxm);
        }
      } else {
	LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: DEBUG Skipping coord ghost read test (zs=0 or dim<3).\n", rank);
      }
    
      // Test reading an owned point (e.g., the first owned point)
      PetscInt k_owned_global = info.zs;
      PetscInt j_owned_global = info.ys;
      PetscInt i_owned_global = info.xs;
      PetscInt k_owned_local = k_owned_global - info.gzs;
      PetscInt j_owned_local = j_owned_global - info.gys;
      PetscInt i_owned_local = i_owned_global - info.gxs;
      if (k_owned_local >= 0 && k_owned_local < info.gzm &&
	  j_owned_local >= 0 && j_owned_local < info.gym &&
	  i_owned_local >= 0 && i_owned_local < info.gxm)
	{
	  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: DEBUG Attempting coord owned read at GLOBAL (%d,%d,%d) -> LOCAL (%d,%d,%d)\n",
		    rank, k_owned_global, j_owned_global, i_owned_global, k_owned_local, j_owned_local, i_owned_local);
	  Cmpnts test_coord_owned = coor_arr[k_owned_local][j_owned_local][i_owned_local];
	  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: DEBUG SUCCESS reading coord owned: x=%g\n", rank, test_coord_owned.x);
	} else {
        LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: DEBUG Calculated coord owned local index (%d,%d,%d) is out of bounds [%d,%d,%d]. Skipping read.\n",
                  rank, k_owned_local, j_owned_local, i_owned_local, info.gzm, info.gym, info.gxm);
      }
    }
    // --- END DEBUG ---
    */
    
    // Loop over OWNED nodes
    for (PetscInt k = zs; k < ze; k++) {
        for (PetscInt j = ys; j < ye; j++) {
            for (PetscInt i = xs; i < xe; i++) {
                // Check if node (i,j,k) is on a GLOBAL boundary
                PetscBool isBoundaryNode = (i == 0 || i == mx - 1 ||
                                            j == 0 || j == my - 1 ||
                                            k == 0 || k == mz - 1);

                if (isBoundaryNode) {
                    PetscInt i_local = i - info.gxs;
                    PetscInt j_local = j - info.gys;
                    PetscInt k_local = k - info.gzs;

		    // --- BEGIN DEBUG: Print indices before access ---
                    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: DEBUG Pre-Access Node(G:%d,%d,%d) -> CoordIndex(L:%d,%d,%d) ArrayBounds(L:%d,%d,%d)\n",
                              rank, i, j, k, i_local, j_local, k_local, info.gxm, info.gym, info.gzm);
                    // --- END DEBUG ---

		    /*
                    // --- Defensive Check ---
                    if (k_local < 0 || k_local >= info.gzm || // Use gzm/gym/gxm for array bounds
                        j_local < 0 || j_local >= info.gym ||
                        i_local < 0 || i_local >= info.gxm) {
                         // Restore arrays before erroring
                         ierr = DMDAVecRestoreArrayRead(vectorDM, localCoords, &coor_arr); CHKERRQ(ierr);
                         ierr = DMRestoreLocalVector(vectorDM, &localCoords); CHKERRQ(ierr);
                         ierr = DMDAVecRestoreArray(vectorDM, fieldVec, &field_arr); CHKERRQ(ierr);
                         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Local coordinate index [%d][%d][%d] out of bounds [%d][%d][%d] during BC application!", k_local, j_local, i_local, info.gzm, info.gym, info.gxm);
                    }
                    // --- End check ---
		    */

		    Cmpnts node_coord = coor_arr[k][j][i];
		    // Cmpnts node_coord = coor_arr[k_local][j_local][i_local];
		    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: DEBUG Accessed Coords OK for local index (%d,%d,%d)\n", rank, k_local, j_local, i_local);
		    
                    // Use existing helper to set field = sin(coord)
                    ierr = SetLocalCartesianField_Vector(fieldName, &field_arr[k][j][i], &node_coord,user->FieldInitialization,user->ConstantVelocity); CHKERRQ(ierr);

                    // TRACE logging can be intensive, maybe DEBUG or remove later
                    LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG,i+j+k,100, "Rank %d: Applied BC Vector at global node (%d,%d,%d) -> val (%.2f,%.2f,%.2f)\n",
                              rank, i, j, k, field_arr[k][j][i].x, field_arr[k][j][i].y, field_arr[k][j][i].z);
                }
            }
        }
    }

    // Restore Access

    ierr = DMDAVecRestoreArrayRead(vectorDM,coordVec,&coor_arr);
    /*
    ierr = DMDAVecRestoreArrayRead(vectorDM, localCoords, &coor_arr); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(vectorDM, &localCoords); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(vectorDM, fieldVec, &field_arr); CHKERRQ(ierr);
    */
    
    PetscFunctionReturn(0);
}

// Placeholder for scalar BCs
static PetscErrorCode ApplyAnalyticalBC_Scalar(
    UserCtx    *user,
    const char *fieldName,
    DM          scalarDM,  // e.g., user->da
    Vec         coordVec,  // e.g., user->coordinates
    Vec         fieldVec   // e.g., user->P or user->Nvert
)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: ApplyAnalyticalBC_Scalar for field '%s' not yet implemented.\n", rank, fieldName);
    // TODO: Implement scalar boundary condition logic
    // - Get scalar field array (e.g., PetscReal ***) using scalarDM and fieldVec
    // - Get coordinate array (careful: coords usually use fda, scalar field uses da - might need separate coord access or interpolation)
    // - Loop over owned SCALAR grid points (cells or nodes depending on convention for P/Nvert)
    // - Check if the point is on a global boundary
    // - If yes, get coordinates (might need interpolation if coords are nodal and field is cell-centered)
    // - Apply analytical scalar function (e.g., using SetLocalCartesianField_Scalar)
    // - Restore arrays
    PetscFunctionReturn(0);
}

/**
 * @brief Applies analytical boundary conditions to a specified global vector or scalar field.
 *
 * This function acts as a dispatcher based on the fieldName. It determines the
 * appropriate DM and Vec from the UserCtx and calls field-specific helper routines
 * (e.g., ApplyAnalyticalBC_Vector, ApplyAnalyticalBC_Scalar) to set boundary values
 * according to a predefined analytical function (e.g., sin(coord)).
 *
 * This should be called AFTER setting interior values and BEFORE updating local ghosts.
 *
 * @param[in] user      Pointer to the UserCtx containing DMs, coordinates, and field Vecs.
 * @param[in] fieldName Name of the field to apply BCs to ("Ucat", "Ucont", "P", "Nvert").
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode ApplyAnalyticalBC(UserCtx *user, const char *fieldName)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    Vec            fieldVec = NULL;
    DM             fieldDM  = NULL;
    Vec            coordVec = NULL;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "Rank %d: Dispatching analytical BC application for field '%s'.\n", rank, fieldName);

    // --- Get Global Coordinates ---
    // Assume coordinates are associated with the 'da' but layout matches 'fda'
    ierr = DMGetCoordinates(user->da, &coordVec); CHKERRQ(ierr);
    if (!coordVec) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "Global coordinates not found/set on user->da.");
    }

    // --- Dispatch based on fieldName ---
    if (strcmp(fieldName, "Ucat") == 0) {
        fieldVec = user->Ucat;
        fieldDM  = user->fda; // Ucat uses fda
        ierr = ApplyAnalyticalBC_Vector(user, fieldName, fieldDM, coordVec, fieldVec); CHKERRQ(ierr);
    } else if (strcmp(fieldName, "Ucont") == 0) {
        fieldVec = user->Ucont;
        fieldDM  = user->fda; // Ucont uses fda
        ierr = ApplyAnalyticalBC_Vector(user, fieldName, fieldDM, coordVec, fieldVec); CHKERRQ(ierr);
    } else if (strcmp(fieldName, "P") == 0) {
        fieldVec = user->P;
        fieldDM  = user->da;  // P uses da
        ierr = ApplyAnalyticalBC_Scalar(user, fieldName, fieldDM, coordVec, fieldVec); CHKERRQ(ierr);
    } else if (strcmp(fieldName, "Nvert") == 0) {
        fieldVec = user->Nvert;
        fieldDM  = user->da;  // Nvert uses da
        ierr = ApplyAnalyticalBC_Scalar(user, fieldName, fieldDM, coordVec, fieldVec); CHKERRQ(ierr);
    } else {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Field '%s' not recognized for ApplyAnalyticalBC.", fieldName);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Finished applying analytical BCs for field '%s'.\n", rank, fieldName);
    PetscFunctionReturn(0);
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
PetscErrorCode SetAnalyticalSolution(Vec tempVec,PetscInt FieldInitialization)
 {
     PetscErrorCode ierr;
     PetscInt nParticles;
     PetscReal *vels;
 
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "SetAnalyticalSolution - Starting analytical solution computation.\n");
 
     ierr = VecGetLocalSize(tempVec, &nParticles); CHKERRQ(ierr);
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "SetAnalyticalSolution - Number of local particles: %d.\n", nParticles);
 
     ierr = VecGetArray(tempVec, &vels); CHKERRQ(ierr);

     for (PetscInt i = 0; i < nParticles; i++) {
         if(FieldInitialization==1){
	   vels[i] = sin(vels[i]);
	 }
	 else if(FieldInitialization==0){
	   if ((i % 3) == 2) { // Check if remainder is 2 when divided by 3
            vels[i] = 1;
	   } else {
            vels[i] = 0;
	   }
	 }
     }
     ierr = VecRestoreArray(tempVec, &vels); CHKERRQ(ierr);
 
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "SetAnalyticalSolution - Completed analytical solution computation.\n");
     return 0;
 }

