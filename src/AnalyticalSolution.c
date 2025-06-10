
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
	//    ConstantValue = user->ConstantVelocity;
    } else if (strcmp(fieldName, "P") == 0) {
        fieldVec = user->P;
        targetDM = user->da;
        accessDM = user->da;
        infoDM_for_cell_loop = user->da;  // Loop over "cells" defined by da's node origins
        fieldIsVector = 0;
        expected_bs = 1;
	//  ConstantValue = user->ConstantPressure;
    } else if (strcmp(fieldName, "Nvert") == 0) {
        fieldVec = user->Nvert;
        targetDM = user->da;
        accessDM = user->da;
        infoDM_for_cell_loop = user->da; // Loop over "cells" defined by da's node origins
        fieldIsVector = 0;
        expected_bs = 1;
	//  ConstantValue = user->ConstantNvert;
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

    //  Vec globalCoor;
    //  ierr = DMGetCoordinates(user->da, &globalCoor); CHKERRQ(ierr);
    //  ierr = DMGlobalToLocalBegin(user->fda, globalCoor, INSERT_VALUES, localCoor); CHKERRQ(ierr);
    //  ierr = DMGlobalToLocalEnd(user->fda, globalCoor, INSERT_VALUES, localCoor); CHKERRQ(ierr);

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
  

/**
 * @brief Sets an analytical contravariant field (Ucont) at nodes, where each component
 *        is evaluated using the physical coordinates of its respective face center.
 *
 * This function:
 * 1. Interpolates nodal physical coordinates to the centers of i-faces, j-faces, and k-faces
 *    for the cells relevant to the current rank, storing them in temporary local arrays.
 * 2. For each node (i,j,k) owned by the current rank (where Ucont components are stored):
 *    a. Retrieves the pre-calculated physical coordinate of the i-face center.
 *    b. Evaluates the analytical function for U^1 using this i-face center coordinate.
 *    c. Retrieves the pre-calculated physical coordinate of the j-face center.
 *    d. Evaluates the analytical function for U^2 using this j-face center coordinate.
 *    e. Retrieves the pre-calculated physical coordinate of the k-face center.
 *    f. Evaluates the analytical function for U^3 using this k-face center coordinate.
 *    g. Stores these (U^1, U^2, U^3) into user->Ucont[k][j][i].
 *
 * @param[in]  user      Pointer to the UserCtx structure.
 * @param[in]  fieldName Name of the field to update (currently only "Ucont").
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
/*
PetscErrorCode SetAnalyticalContravariantField(UserCtx *user, const char *fieldName)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    if (strcmp(fieldName, "Ucont") != 0) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Field '%s' not supported by this version of SetAnalyticalContravariantField. Only 'Ucont' is supported.", fieldName);
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,
                   "Starting for field '%s' (face-centered logic) on rank %d.\n", fieldName, rank);

    DM             fda = user->fda;
    DMDALocalInfo  info_fda;
    ierr = DMDAGetLocalInfo(fda, &info_fda); CHKERRQ(ierr);

    // Get array for nodal physical coordinates (local, ghosted)
    Cmpnts ***nodal_phy_coords_arr;
    Vec localCoor;

    DMGetCoordinatesLocal(user->da, &localCoor);
    ierr = DMDAVecGetArrayRead(fda,localCoor, &nodal_phy_coords_arr); CHKERRQ(ierr);

    // Determine owned CELL ranges (xs, xm, etc.) which dictates the size of face arrays
    PetscInt xs_cell, xm_cell, ys_cell, ym_cell, zs_cell, zm_cell;
    ierr = GetOwnedCellRange(&info_fda, 0, &xs_cell, &xm_cell); CHKERRQ(ierr); // info_fda because faces are tied to fda nodes
    ierr = GetOwnedCellRange(&info_fda, 1, &ys_cell, &ym_cell); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_fda, 2, &zs_cell, &zm_cell); CHKERRQ(ierr);

    // Allocate temporary local C-arrays for face-centered physical coordinates
    Cmpnts ***local_iface_phy_coords, ***local_jface_phy_coords, ***local_kface_phy_coords;
    // Dimensions based on your InterpolateCornerToFaceCenter_Vector output:
    // faceX_arr: [zm_cell][ym_cell][xm_cell+1]
    // faceY_arr: [zm_cell][ym_cell+1][xm_cell]
    // faceZ_arr: [zm_cell+1][ym_cell][xm_cell]
    if (xm_cell > 0 || ym_cell > 0 || zm_cell > 0) { // Only allocate if there are cells
        ierr = Allocate3DArray(&local_iface_phy_coords, zm_cell, ym_cell, xm_cell + 1); CHKERRQ(ierr);
        ierr = Allocate3DArray(&local_jface_phy_coords, zm_cell, ym_cell + 1, xm_cell); CHKERRQ(ierr);
        ierr = Allocate3DArray(&local_kface_phy_coords, zm_cell + 1, ym_cell, xm_cell); CHKERRQ(ierr);
    } else { // No cells on this rank, so no faces to compute coords for.
        local_iface_phy_coords = NULL;
        local_jface_phy_coords = NULL;
        local_kface_phy_coords = NULL;
    }


    // Interpolate nodal physical coordinates to face centers (populates the local_***_phy_coords arrays)
    if (local_iface_phy_coords) { // Check if allocation happened
        ierr = InterpolateCornerToFaceCenter(nodal_phy_coords_arr,
                                                    local_iface_phy_coords,
                                                    local_jface_phy_coords,
                                                    local_kface_phy_coords,
                                                    user); CHKERRQ(ierr);
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Physical coordinates interpolated to face centers.\n", rank);
    }


    // Get WRITE array for Ucont (Global Vector user->Ucont)
    Cmpnts ***ucont_arr_global;
    ierr = DMDAVecGetArray(fda, user->Ucont, &ucont_arr_global); CHKERRQ(ierr);
    ierr = VecSet(user->Ucont, 0.0);CHKERRQ(ierr); // Initialize

    // Loop over OWNED NODES of user->fda (where Ucont components are stored)
    PetscInt xs_node_own = info_fda.xs, xe_node_own = info_fda.xs + info_fda.xm;
    PetscInt ys_node_own = info_fda.ys, ye_node_own = info_fda.ys + info_fda.ym;
    PetscInt zs_node_own = info_fda.zs, ze_node_own = info_fda.zs + info_fda.zm;

    // Cmpnts ConstantValue = user->ConstantContra; // From UserCtx

    for (PetscInt k_glob_node = zs_node_own; k_glob_node < ze_node_own; k_glob_node++) {
        for (PetscInt j_glob_node = ys_node_own; j_glob_node < ye_node_own; j_glob_node++) {
            for (PetscInt i_glob_node = xs_node_own; i_glob_node < xe_node_own; i_glob_node++) {
                Cmpnts ucont_val_at_node = {0.0, 0.0, 0.0}; // Initialize
                Cmpnts face_coord;

                // --- U^1 component (conceptually on i-face at this node's i-coordinate) ---
                // i_glob_node is the global index of the i-face.
                // local_iface_phy_coords is indexed [k_loc_cell][j_loc_cell][i_loc_faceX]
                // i_loc_faceX = i_glob_node - xs_cell (ranges 0 to xm_cell for owned faces)
                PetscInt k_loc_cell_for_iface = k_glob_node - zs_cell;
                PetscInt j_loc_cell_for_iface = j_glob_node - ys_cell;
                PetscInt i_loc_face           = i_glob_node - xs_cell;

                if (local_iface_phy_coords &&
                    k_loc_cell_for_iface >= 0 && k_loc_cell_for_iface < zm_cell &&
                    j_loc_cell_for_iface >= 0 && j_loc_cell_for_iface < ym_cell &&
                    i_loc_face >= 0 && i_loc_face <= xm_cell) { // Note: i_loc_face goes up to xm_cell (xm_cell+1 values)
                    face_coord = local_iface_phy_coords[k_loc_cell_for_iface][j_loc_cell_for_iface][i_loc_face];

                    if(user->FieldInitialization == 1) ucont_val_at_node.x = sin(face_coord.x);
                    else if(user->FieldInitialization == 2) ucont_val_at_node.x = cos(face_coord.x)*sin(face_coord.y)*sin(face_coord.z);
                    else if(user->FieldInitialization == 0) ucont_val_at_node.x = ConstantValue.x;
                } else if (local_iface_phy_coords) {
                     // This case means the node (i_glob_node, j_glob_node, k_glob_node) is outside the range
                     // of cells for which local_iface_phy_coords were computed (e.g. node on max boundary of domain)
                     // Or, xm_cell/ym_cell/zm_cell is 0 for this rank.
                     // Set to default, or could use nodal coordinate as fallback.
                     // This requires careful thought on boundary conditions for Ucont.
                     // For now, if not found, it will be zero from initialization.
                }


                // --- U^2 component (conceptually on j-face at this node's j-coordinate) ---
                // j_glob_node is the global index of the j-face.
                // local_jface_phy_coords is indexed [k_loc_cell][j_loc_faceY][i_loc_cell]
                PetscInt k_loc_cell_for_jface = k_glob_node - zs_cell;
                PetscInt j_loc_face           = j_glob_node - ys_cell;
                PetscInt i_loc_cell_for_jface = i_glob_node - xs_cell;

                if (local_jface_phy_coords &&
                    k_loc_cell_for_jface >= 0 && k_loc_cell_for_jface < zm_cell &&
                    j_loc_face >= 0 && j_loc_face <= ym_cell &&
                    i_loc_cell_for_jface >= 0 && i_loc_cell_for_jface < xm_cell) {
                    face_coord = local_jface_phy_coords[k_loc_cell_for_jface][j_loc_face][i_loc_cell_for_jface];

                    if(user->FieldInitialization == 1) ucont_val_at_node.y = sin(face_coord.y);
                    else if(user->FieldInitialization == 2) ucont_val_at_node.y = sin(face_coord.x)*cos(face_coord.y)*sin(face_coord.z);
                    else if(user->FieldInitialization == 0) ucont_val_at_node.y = ConstantValue.y;
                } // else fallback to zero

                // --- U^3 component (conceptually on k-face at this node's k-coordinate) ---
                // k_glob_node is the global index of the k-face.
                // local_kface_phy_coords is indexed [k_loc_faceZ][j_loc_cell][i_loc_cell]
                PetscInt k_loc_face           = k_glob_node - zs_cell;
                PetscInt j_loc_cell_for_kface = j_glob_node - ys_cell;
                PetscInt i_loc_cell_for_kface = i_glob_node - xs_cell;

                if (local_kface_phy_coords &&
                    k_loc_face >= 0 && k_loc_face <= zm_cell &&
                    j_loc_cell_for_kface >= 0 && j_loc_cell_for_kface < ym_cell &&
                    i_loc_cell_for_kface >= 0 && i_loc_cell_for_kface < xm_cell) {
                    face_coord = local_kface_phy_coords[k_loc_face][j_loc_cell_for_kface][i_loc_cell_for_kface];

                    if(user->FieldInitialization == 1) ucont_val_at_node.z = sin(face_coord.z);
                    else if(user->FieldInitialization == 2) ucont_val_at_node.z = sin(face_coord.x)*sin(face_coord.y)*cos(face_coord.z);
                    else if(user->FieldInitialization == 0) ucont_val_at_node.z = ConstantValue.z;
                } // else fallback to zero


                // Store computed Ucont components into the global Ucont array at the current node
                ucont_arr_global[k_glob_node][j_glob_node][i_glob_node] = ucont_val_at_node;
            }
        }
    }

    // Restore and deallocate
    ierr = DMDAVecRestoreArrayRead(fda,localCoor, &nodal_phy_coords_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(fda, user->Ucont, &ucont_arr_global); CHKERRQ(ierr);

    if (local_iface_phy_coords) { ierr = Deallocate3DArray(local_iface_phy_coords, zm_cell, ym_cell); CHKERRQ(ierr); }
    if (local_jface_phy_coords) { ierr = Deallocate3DArray(local_jface_phy_coords, zm_cell, ym_cell + 1); CHKERRQ(ierr); }
    if (local_kface_phy_coords) { ierr = Deallocate3DArray(local_kface_phy_coords, zm_cell + 1, ym_cell); CHKERRQ(ierr); }

    ierr = VecAssemblyBegin(user->Ucont); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(user->Ucont); CHKERRQ(ierr);

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "SetAnalyticalContravariantField: Completed for field '%s' on rank %d.\n", fieldName, rank);
    MPI_Barrier(PETSC_COMM_WORLD); // Good for debugging field setup
    return 0;
}


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

    
    // Get Coordinates (using the same DM as the field for consistency)
   // ierr = DMGetLocalVector(vectorDM, &localCoords); CHKERRQ(ierr);
   // ierr = DMGlobalToLocalBegin(vectorDM, coordVec, INSERT_VALUES, localCoords); CHKERRQ(ierr);
   // ierr = DMGlobalToLocalEnd(vectorDM, coordVec, INSERT_VALUES, localCoords); CHKERRQ(ierr);
   // ierr = DMDAVecGetArrayRead(vectorDM, localCoords, &coor_arr); CHKERRQ(ierr);
   // LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Obtained coord array pointer %p\n", rank, (void*)coor_arr);
    

    DMDAVecGetArrayRead(vectorDM,coordVec,&coor_arr);
    */
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
    /*
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
		    */
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
                    /*
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
    
    //ierr = DMDAVecRestoreArrayRead(vectorDM, localCoords, &coor_arr); CHKERRQ(ierr);
    //ierr = DMRestoreLocalVector(vectorDM, &localCoords); CHKERRQ(ierr);
    //ierr = DMDAVecRestoreArray(vectorDM, fieldVec, &field_arr); CHKERRQ(ierr);
    
    
    PetscFunctionReturn(0);
}
*/

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
 * @brief Applies analytical boundary conditions to a contravariant vector field (e.g., Ucont).
 *
 * For each component of the contravariant vector field stored at a boundary node,
 * this function determines the physical coordinates of the center of the corresponding
 * face (i-face for .x component, j-face for .y, k-face for .z). It then evaluates
 * the analytical function using these face-center coordinates to set the boundary value
 * for that specific component at that node.
 *
 * This function is self-contained: it fetches nodal coordinates, interpolates them
 * to face centers, applies the logic, and cleans up temporary resources.
 *
 * @param user         User context, providing parameters like FieldInitialization and ConstantContra.
 * @param fieldName    Name of the contravariant field (e.g., "Ucont").
 * @param dm_field     DMDA for the contravariant field (e.g., user->fda).
 * @param fieldVec_global      Global PETSc Vec for the contravariant field being set (e.g., user->Ucont).
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
/*
static PetscErrorCode ApplyAnalyticalBC_Contra_Vector(
    UserCtx    *user,
    const char *fieldName,
    DM          dm_field,
    Vec         globalCoor,
    Vec         fieldVec_global)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    Cmpnts       ***field_arr_global;
    Cmpnts       ***nodal_phy_coords_arr;
    Vec            localCoords_with_ghosts; // This will be our temporary local vector

    Cmpnts ***local_iface_phy_coords = NULL;
    Cmpnts ***local_jface_phy_coords = NULL;
    Cmpnts ***local_kface_phy_coords = NULL;

    DMDALocalInfo  info_field_dm;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Applying Face-Centered BCs for Contravariant field '%s'.\n", rank, fieldName);

    ierr = DMDAGetLocalInfo(dm_field, &info_field_dm); CHKERRQ(ierr);

    // --- 1. Get Nodal Physical Coordinates (Local Ghosted Array) ---
    // THIS IS THE CORRECTED PATTERN:
    // a. Create an empty temporary local vector compatible with the DM.
    ierr = DMGetLocalVector(dm_field, &localCoords_with_ghosts); CHKERRQ(ierr);
    if (!localCoords_with_ghosts) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_MEM, "Failed to get local vector for coordinates.");

    // b. Scatter data from the global coordinate vector into the new local vector.
    ierr = DMGlobalToLocalBegin(dm_field, globalCoor, INSERT_VALUES, localCoords_with_ghosts); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dm_field, globalCoor, INSERT_VALUES, localCoords_with_ghosts); CHKERRQ(ierr);
    
    ierr = DMDAVecGetArrayRead(dm_field, localCoords_with_ghosts, &nodal_phy_coords_arr); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d:Local Coordinates read .\n", rank);
    
    // --- 2. Interpolate Nodal Coordinates to Face Centers ---
    // Determine owned CELL ranges to size temp face coord arrays.
    PetscInt xs_cell, xm_cell, ys_cell, ym_cell, zs_cell, zm_cell;
    ierr = GetOwnedCellRange(&info_field_dm, 0, &xs_cell, &xm_cell); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_field_dm, 1, &ys_cell, &ym_cell); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_field_dm, 2, &zs_cell, &zm_cell); CHKERRQ(ierr);

    if (xm_cell > 0 || ym_cell > 0 || zm_cell > 0) { // Only if this rank has cells to define faces
        ierr = Allocate3DArray(&local_iface_phy_coords, zm_cell, ym_cell, xm_cell + 1); CHKERRQ(ierr);
        ierr = Allocate3DArray(&local_jface_phy_coords, zm_cell, ym_cell + 1, xm_cell); CHKERRQ(ierr);
        ierr = Allocate3DArray(&local_kface_phy_coords, zm_cell + 1, ym_cell, xm_cell); CHKERRQ(ierr);

        // Interpolate using the generic macro dispatcher
        ierr = InterpolateCornerToFaceCenter(nodal_phy_coords_arr,
                                             local_iface_phy_coords,
                                             local_jface_phy_coords,
                                             local_kface_phy_coords,
                                             user); CHKERRQ(ierr);
    }

    // --- 3. Apply Boundary Conditions ---
    // Get the array for the contravariant field vector (e.g., Ucont)
    ierr = DMDAVecGetArray(dm_field, fieldVec_global, &field_arr_global); CHKERRQ(ierr);

    // Loop over OWNED nodes of dm_field
    PetscInt xs_node_own = info_field_dm.xs, xe_node_own = info_field_dm.xs + info_field_dm.xm;
    PetscInt ys_node_own = info_field_dm.ys, ye_node_own = info_field_dm.ys + info_field_dm.ym;
    PetscInt zs_node_own = info_field_dm.zs, ze_node_own = info_field_dm.zs + info_field_dm.zm;

    // Global domain dimensions (number of NODES on dm_field)
    PetscInt mx_global_nodes = info_field_dm.mx;
    PetscInt my_global_nodes = info_field_dm.my;
    PetscInt mz_global_nodes = info_field_dm.mz;

    // Cmpnts ConstantValue = user->ConstantContra; // Correctly use the Cmpnts struct

    for (PetscInt k_glob_node = zs_node_own; k_glob_node < ze_node_own; k_glob_node++) {
        for (PetscInt j_glob_node = ys_node_own; j_glob_node < ye_node_own; j_glob_node++) {
            for (PetscInt i_glob_node = xs_node_own; i_glob_node < xe_node_own; i_glob_node++) {

                PetscBool isBoundaryNode = (i_glob_node == 0 || i_glob_node == mx_global_nodes - 1 ||
                                            j_glob_node == 0 || j_glob_node == my_global_nodes - 1 ||
                                            k_glob_node == 0 || k_glob_node == mz_global_nodes - 1);

                if (isBoundaryNode) {
                    Cmpnts field_bc_val; // Will be fully set component-wise
                    Cmpnts face_coord;

                    // --- .x component (for I-Face i_glob_node) ---
                    PetscInt k_loc_c_iface = k_glob_node - zs_cell;
                    PetscInt j_loc_c_iface = j_glob_node - ys_cell;
                    PetscInt i_loc_iface   = i_glob_node - xs_cell;
                    if (local_iface_phy_coords && // Check if allocated (rank has cells)
                        k_loc_c_iface >= 0 && k_loc_c_iface < zm_cell &&
                        j_loc_c_iface >= 0 && j_loc_c_iface < ym_cell &&
                        i_loc_iface >= 0 && i_loc_iface <= xm_cell) {
                        face_coord = local_iface_phy_coords[k_loc_c_iface][j_loc_c_iface][i_loc_iface];
                    } else { // Fallback if face coords not available
                        face_coord = nodal_phy_coords_arr[k_glob_node][j_glob_node][i_glob_node];
                        LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Contravariant BC for .x comp at node (%d,%d,%d) using fallback (nodal coord).\n", rank, i_glob_node, j_glob_node, k_glob_node);
                    }
                    if(user->FieldInitialization == 1)      field_bc_val.x = sin(face_coord.x);
                    else if(user->FieldInitialization == 2) field_bc_val.x = cos(face_coord.x)*sin(face_coord.y)*sin(face_coord.z);
                    else if(user->FieldInitialization == 0) field_bc_val.x = ConstantValue.x;


                    // --- .y component (for J-Face j_glob_node) ---
                    PetscInt k_loc_c_jface = k_glob_node - zs_cell;
                    PetscInt j_loc_jface   = j_glob_node - ys_cell;
                    PetscInt i_loc_c_jface = i_glob_node - xs_cell;
                    if (local_jface_phy_coords &&
                        k_loc_c_jface >= 0 && k_loc_c_jface < zm_cell &&
                        j_loc_jface >= 0 && j_loc_jface <= ym_cell &&
                        i_loc_c_jface >= 0 && i_loc_c_jface < xm_cell) {
                        face_coord = local_jface_phy_coords[k_loc_c_jface][j_loc_jface][i_loc_c_jface];
                    }  else { // Fallback
                        face_coord = nodal_phy_coords_arr[k_glob_node][j_glob_node][i_glob_node];
                        LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Contravariant BC for .y comp at node (%d,%d,%d) using fallback.\n", rank, i_glob_node, j_glob_node, k_glob_node);
                    }
                    if(user->FieldInitialization == 1)      field_bc_val.y = sin(face_coord.y);
                    else if(user->FieldInitialization == 2) field_bc_val.y = sin(face_coord.x)*cos(face_coord.y)*sin(face_coord.z);
                    else if(user->FieldInitialization == 0) field_bc_val.y = ConstantValue.y;


                    // --- .z component (for K-Face k_glob_node) ---
                    PetscInt k_loc_kface   = k_glob_node - zs_cell;
                    PetscInt j_loc_c_kface = j_glob_node - ys_cell;
                    PetscInt i_loc_c_kface = i_glob_node - xs_cell;
                    if (local_kface_phy_coords &&
                        k_loc_kface >= 0 && k_loc_kface <= zm_cell &&
                        j_loc_c_kface >= 0 && j_loc_c_kface < ym_cell &&
                        i_loc_c_kface >= 0 && i_loc_c_kface < xm_cell) {
                        face_coord = local_kface_phy_coords[k_loc_kface][j_loc_c_kface][i_loc_c_kface];
                    } else { // Fallback
                        face_coord = nodal_phy_coords_arr[k_glob_node][j_glob_node][i_glob_node];
                        LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Contravariant BC for .z comp at node (%d,%d,%d) using fallback.\n", rank, i_glob_node, j_glob_node, k_glob_node);
                    }
                    if(user->FieldInitialization == 1)      field_bc_val.z = sin(face_coord.z);
                    else if(user->FieldInitialization == 2) field_bc_val.z = sin(face_coord.x)*sin(face_coord.y)*cos(face_coord.z);
                    else if(user->FieldInitialization == 0) field_bc_val.z = ConstantValue.z;

                    field_arr_global[k_glob_node][j_glob_node][i_glob_node] = field_bc_val;
                } // end if isBoundaryNode
            } // end i_glob_node
        } // end j_glob_node
    } // end k_glob_node

    ierr = DMDAVecRestoreArray(dm_field, fieldVec_global, &field_arr_global); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(dm_field, localCoords_with_ghosts, &nodal_phy_coords_arr); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(dm_field, &localCoords_with_ghosts); CHKERRQ(ierr); // MUST restore the temp local Vec
    
    if (local_iface_phy_coords) { ierr = Deallocate3DArray(local_iface_phy_coords, zm_cell, ym_cell); CHKERRQ(ierr); }
    if (local_jface_phy_coords) { ierr = Deallocate3DArray(local_jface_phy_coords, zm_cell, ym_cell + 1); CHKERRQ(ierr); }
    if (local_kface_phy_coords) { ierr = Deallocate3DArray(local_kface_phy_coords, zm_cell + 1, ym_cell); CHKERRQ(ierr); }

    PetscFunctionReturn(0);
}
*/

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
/*
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
        ierr = ApplyAnalyticalBC_Contra_Vector(user, fieldName, fieldDM, coordVec, fieldVec); CHKERRQ(ierr);
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
*/
////////// FOR ERROR CALCULATION //////////

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


/**
 * @brief Sets the initial values for the INTERIOR of a specified Eulerian field.
 *
 * This function initializes the interior nodes of `Ucont` based on a profile selected
 * by `user->FieldInitialization`. It explicitly skips any node that lies on a global
 * boundary, as those values are set by the Boundary System's `Initialize` methods.
 *
 * The initialization is directional, aligned with the primary INLET face that was
 * identified by the parser. This ensures the initial flow is physically meaningful.
 *
 * Supported `user->FieldInitialization` profiles for "Ucont":
 *  - 0: Zero Velocity. All interior components of Ucont are set to 0.
 *  - 1: Constant Normal Velocity. The contravariant velocity component normal to the
 *       inlet direction is set such that the physical velocity normal to those grid
 *       planes is a constant `uin`. Other contravariant components are zero.
 *  - 2: Poiseuille Normal Velocity. The contravariant component normal to the
 *       inlet direction is set with a parabolic profile.
 *
 * @param user      The main UserCtx struct, containing all simulation data and configuration.
 * @param fieldName A string ("Ucont" or "P") identifying which field to initialize.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetInitialInteriorField(UserCtx *user, const char *fieldName)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    LOG_ALLOW(GLOBAL, LOG_INFO, "Setting initial INTERIOR field for '%s' with profile %d.\n", fieldName, user->FieldInitialization);

    // This function currently only implements logic for Ucont.
    if (strcmp(fieldName, "Ucont") != 0) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Skipping SetInitialInteriorField for non-Ucont field '%s'.\n", fieldName);
        PetscFunctionReturn(0);
    }

    // --- 1. Get references to required data and PETSc arrays ---
    DM            fieldDM = user->fda;
    Vec           fieldVec = user->Ucont;
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(fieldDM, &info); CHKERRQ(ierr);

    Vec      localCoor;
    Cmpnts ***coor_arr;
    ierr = DMGetCoordinatesLocal(user->da, &localCoor); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, localCoor, &coor_arr); CHKERRQ(ierr);
    
    Cmpnts ***csi_arr, ***eta_arr, ***zet_arr;
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, &csi_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, &eta_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, &zet_arr); CHKERRQ(ierr);
   
    // --- 2. Compute Cell-Center Coordinates (only if needed by the selected profile) ---
    Cmpnts ***cent_coor = NULL;
    PetscInt xs_cell=0, xm_cell=0, ys_cell=0, ym_cell=0, zs_cell=0, zm_cell=0;
    
    if (user->FieldInitialization == 2) { // Profile 2 (Poiseuille) requires cell centers.
        ierr = GetOwnedCellRange(&info, 0, &xs_cell, &xm_cell); CHKERRQ(ierr);
        ierr = GetOwnedCellRange(&info, 1, &ys_cell, &ym_cell); CHKERRQ(ierr);
        ierr = GetOwnedCellRange(&info, 2, &zs_cell, &zm_cell); CHKERRQ(ierr);

        if (xm_cell > 0 && ym_cell > 0 && zm_cell > 0) {
            ierr = Allocate3DArray(&cent_coor, zm_cell, ym_cell, xm_cell); CHKERRQ(ierr);
            ierr = InterpolateFieldFromCornerToCenter(coor_arr, cent_coor, user); CHKERRQ(ierr);
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Computed temporary cell-center coordinates for Poiseuille profile.\n");
        }
    }
    
    // --- 3. Loop Over Owned Nodes and Apply Initial Condition to Interior ---
    Cmpnts ***ucont_arr;
    ierr = DMDAVecGetArray(fieldDM, fieldVec, &ucont_arr); CHKERRQ(ierr);
    
    PetscInt i, j, k;
    const PetscInt mx = info.mx, my = info.my, mz = info.mz; // Global node dimensions
    const PetscInt xs = info.xs, xe = info.xs + info.xm;
    const PetscInt ys = info.ys, ye = info.ys + info.ym;
    const PetscInt zs = info.zs, ze = info.zs + info.zm;
    
    const Cmpnts uin = user->InitialConstantContra; // Max/average velocity from user options
    // Flow into a negative face (e.g., -Zeta at k=0) is in the positive physical direction (+z).
    // Flow into a positive face (e.g., +Zeta at k=mz-1) is in the negative physical direction (-z).

    LOG_ALLOW(GLOBAL,LOG_DEBUG,"InitialConstantContra = (%.3f, %.3f, %.3f)\n",(double)uin.x, (double)uin.y, (double)uin.z);
    
    const PetscReal flow_direction_sign = (user->identifiedInletBCFace % 2 == 0) ? -1.0 : 1.0;
        
    for (k = zs; k < ze; k++) {
        for (j = ys; j < ye; j++) {
            for (i = xs; i < xe; i++) {
                
                // The crucial check to ensure we only modify interior nodes.
                const PetscBool is_interior = (i > 0 && i < mx - 1 &&
                                               j > 0 && j < my - 1 &&
                                               k > 0 && k < mz - 1);

                if (is_interior) {
                    Cmpnts ucont_val = {0.0, 0.0, 0.0}; // Default to zero velocity
                    PetscReal normal_velocity_mag = 0.0;

                    // Step A: Determine the magnitude of the desired physical normal velocity.
                    switch (user->FieldInitialization) {
                        case 0: // Zero initial velocity
                            normal_velocity_mag = 0.0;
                            break;
 		        case 1: /* Constant Normal Velocity */
		            if (user->identifiedInletBCFace == BC_FACE_NEG_X ||
			    user->identifiedInletBCFace == BC_FACE_POS_X) {
			    normal_velocity_mag = uin.x;
		            } else if (user->identifiedInletBCFace == BC_FACE_NEG_Y || user->identifiedInletBCFace == BC_FACE_POS_Y) {
			   normal_velocity_mag = uin.y;
		            } else {
			   normal_velocity_mag = uin.z;
		            }
		      break; 
                        case 2: // Poiseuille Normal Velocity
                            {
                                PetscReal r_sq = 0.0;
                                const PetscInt i_local = i - xs_cell, j_local = j - ys_cell, k_local = k - zs_cell;
                                if (cent_coor && i_local >= 0 && i_local < xm_cell && j_local >= 0 && j_local < ym_cell && k_local >= 0 && k_local < zm_cell) {
                                    const Cmpnts* center = &cent_coor[k_local][j_local][i_local];
                                    if (user->identifiedInletBCFace <= BC_FACE_POS_X) r_sq = center->y * center->y + center->z * center->z;
                                    else if (user->identifiedInletBCFace <= BC_FACE_POS_Y) r_sq = center->x * center->x + center->z * center->z;
                                    else r_sq = center->x * center->x + center->y * center->y;
				    /* pick the correct contravariant component for center‐line speed */
				    PetscReal u0;
				    if (user->identifiedInletBCFace == BC_FACE_NEG_X ||
					user->identifiedInletBCFace == BC_FACE_POS_X) {
				      u0 = uin.x;
				    } else if (user->identifiedInletBCFace == BC_FACE_NEG_Y ||
					       user->identifiedInletBCFace == BC_FACE_POS_Y) {
				      u0 = uin.y;
				    } else {
				      u0 = uin.z;
				    }
				    /* now form the parabolic profile as before */
                                    normal_velocity_mag = 2.0 * u0 * (1.0 - 4.0 * r_sq);
                                }
                            }
                            break;
                        default:
                            LOG_ALLOW(LOCAL, LOG_WARNING, "Unrecognized FieldInitialization profile %d. Defaulting to zero.\n", user->FieldInitialization);
                            normal_velocity_mag = 0.0;
                            break;
                    }

                    // Step B: Apply direction sign and set the correct contravariant component.
                    // The contravariant component U^n = v_n * Area_n, where v_n is the physical normal velocity.
                    if (normal_velocity_mag != 0.0) {
                        const PetscReal signed_normal_vel = normal_velocity_mag * flow_direction_sign;
                        
                        if (user->identifiedInletBCFace == BC_FACE_NEG_X || user->identifiedInletBCFace == BC_FACE_POS_X) {
                            const PetscReal area_i = sqrt(csi_arr[k][j][i].x * csi_arr[k][j][i].x + csi_arr[k][j][i].y * csi_arr[k][j][i].y + csi_arr[k][j][i].z * csi_arr[k][j][i].z);
			    
                            ucont_val.x = signed_normal_vel * area_i;

			    LOG_LOOP_ALLOW(GLOBAL,LOG_DEBUG,k*j,500," ucont_val.x = %.6f (signed_normal_vel=%.3f × area=%.4f)\n",ucont_val.x, signed_normal_vel, area_i);
                        } 
                        else if (user->identifiedInletBCFace == BC_FACE_NEG_Y || user->identifiedInletBCFace == BC_FACE_POS_Y) {
                            const PetscReal area_j = sqrt(eta_arr[k][j][i].x * eta_arr[k][j][i].x + eta_arr[k][j][i].y * eta_arr[k][j][i].y + eta_arr[k][j][i].z * eta_arr[k][j][i].z);
			     
                            ucont_val.y = signed_normal_vel * area_j;

			    LOG_LOOP_ALLOW(GLOBAL,LOG_DEBUG,k*i,500," ucont_val.y = %.6f (signed_normal_vel=%.3f × area=%.4f)\n",ucont_val.y, signed_normal_vel, area_j);
                        } 
                        else { // Z-inlet
                            const PetscReal area_k = sqrt(zet_arr[k][j][i].x * zet_arr[k][j][i].x + zet_arr[k][j][i].y * zet_arr[k][j][i].y + zet_arr[k][j][i].z * zet_arr[k][j][i].z);

                            ucont_val.z = signed_normal_vel * area_k;

			    LOG_LOOP_ALLOW(GLOBAL,LOG_DEBUG,i*j,500," i,j,k,ucont_val.z = %d, %d, %d, %.6f (signed_normal_vel=%.3f × area=%.4f)\n",i,j,k,ucont_val.z, signed_normal_vel, area_k);
                        }
                    }
                    ucont_arr[k][j][i] = ucont_val;
                } // end if(is_interior)
            }
        }
    }
    ierr = DMDAVecRestoreArray(fieldDM, fieldVec, &ucont_arr); CHKERRQ(ierr);

    // --- 5. Cleanup: Restore arrays and free temporary memory ---
    ierr = DMDAVecRestoreArrayRead(user->fda, localCoor, &coor_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, &csi_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, &eta_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, &zet_arr); CHKERRQ(ierr);
    
    if (cent_coor) {
        ierr = Deallocate3DArray(cent_coor, zm_cell, ym_cell); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
