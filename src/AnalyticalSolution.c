
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
 * This function looks up the field within the user context by comparing the provided field name
 * against supported names ("Ucat", "Ucont", "P", "Nvert").
 *
 * It calculates the analytical value based on the coordinate of the corresponding cell center.
 * For scalar fields (P, Nvert), the value is stored directly at the cell index.
 * For vector fields (Ucat, Ucont), adhering to the original code's convention, the value
 * representing the center of cell C(i,j,k) is stored at the NODE index [k][j][i], but ONLY
 * for interior nodes. Boundary node indices are left untouched by this function and must
 * be populated separately by boundary condition routines (like FormBCS).
 *
 * @param[in]  user      Pointer to the UserCtx structure containing DMs and field Vecs.
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
  DM  targetDM = NULL; // DM associated with the target field
  DM  accessDM = NULL; // DM used for Get/Restore Array (matches targetDM)
  PetscInt fieldIsVector = -1;
  PetscInt expected_bs;
  PetscReal ConstantValue;

  if (strcmp(fieldName, "Ucat") == 0) {
    fieldVec = user->Ucat;
    targetDM = user->fda; // Ucat uses fda (DOF=3, node-indexed)
    accessDM = user->fda;
    fieldIsVector = 1;
    expected_bs = 3;
    ConstantValue = user->ConstantVelocity;
  } else if (strcmp(fieldName, "Ucont") == 0) {
    fieldVec = user->Ucont;
    targetDM = user->fda; // Ucont uses fda (DOF=3, node-indexed)
    accessDM = user->fda;
    fieldIsVector = 1;
    expected_bs = 3;
     ConstantValue = user->ConstantContra;
  } else if (strcmp(fieldName, "P") == 0) {
    fieldVec = user->P;
    targetDM = user->da;  // P uses da (DOF=1, conceptually cell-based)
    accessDM = user->da;
    fieldIsVector = 0;
    expected_bs = 1;
    ConstantValue = user->ConstantPressure;
  } else if (strcmp(fieldName, "Nvert") == 0) {
    fieldVec = user->Nvert;
    targetDM = user->da;  // Nvert uses da (DOF=1, conceptually cell-based)
    accessDM = user->da;
    fieldIsVector = 0;
    expected_bs = 1;
    ConstantValue = user->ConstantNvert;
  } else {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
             "Field '%s' not found in user context", fieldName);
  }

  /* --- 2. Verify Target DM Block Size --- */
  PetscInt bs;
  ierr = DMGetBlockSize(targetDM, &bs); CHKERRQ(ierr);
  if (bs != expected_bs) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
             "Expected block size %d for field '%s' based on its DM, got %d", expected_bs, fieldName, bs);
  }
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Target DM block size verified for %s on rank %d \n", fieldName,rank);

  /* --- 3. Get Coordinates and Calculate Cell Center Coordinates --- */
  
  // Get local, ghosted coordinate vector (associated with da, layout fda)
  Vec localCoor;
  ierr = DMGetCoordinatesLocal(user->da, &localCoor); CHKERRQ(ierr);
  // !!! IMPORTANT: Ensure localCoor ghosts are up-to-date before Interpolate... call
  // Example: If user->coordinates is the global coordinate vector:
  Vec globalCoor;
  ierr = DMGetCoordinates(user->da, &globalCoor); CHKERRQ(ierr);
  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Before Coord G2L: fda=%p, globalCoor=%p, localCoor=%p\n", rank, (void*)user->fda, (void*)globalCoor, (void*)localCoor);
  // For debugging, to ensure the globalCoor is created properly.
  // ierr = VecView(globalCoor, PETSC_VIEWER_STDOUT_WORLD);
  ierr = DMGlobalToLocalBegin(user->fda, globalCoor, INSERT_VALUES, localCoor); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->fda, globalCoor, INSERT_VALUES, localCoor); CHKERRQ(ierr);
  // For debugging, to ensure the localCoor is created properly.
  // ierr = VecView(localCoor, PETSC_VIEWER_STDOUT_WORLD);
  
  MPI_Barrier(PETSC_COMM_WORLD);

  // Get read access to coordinate array using fda layout (node indices)
  Cmpnts ***coor_arr;
  ierr = DMDAVecGetArrayRead(user->fda, localCoor, &coor_arr); CHKERRQ(ierr);
  LOG_ALLOW(LOCAL, LOG_DEBUG, "Obtained read access to local coordinate array on rank %d.\n",rank);

  // Get local node info from fda (needed for cell range calculation and interior check)
  DMDALocalInfo info_nodes;
  ierr = DMDAGetLocalInfo(user->fda, &info_nodes); CHKERRQ(ierr);

  // Determine owned CELL ranges using helper
  PetscInt xs_cell, xm_cell, xe_cell;
  PetscInt ys_cell, ym_cell, ye_cell;
  PetscInt zs_cell, zm_cell, ze_cell;
  ierr = GetOwnedCellRange(&info_nodes, 0, &xs_cell, &xm_cell); CHKERRQ(ierr);
  ierr = GetOwnedCellRange(&info_nodes, 1, &ys_cell, &ym_cell); CHKERRQ(ierr);
  ierr = GetOwnedCellRange(&info_nodes, 2, &zs_cell, &zm_cell); CHKERRQ(ierr);
  xe_cell = xs_cell + xm_cell;
  ye_cell = ys_cell + ym_cell;
  ze_cell = zs_cell + zm_cell;

  // Allocate temporary local array for cell-center coordinates using the generic macro
  // Size based on the number of owned cells (xm_cell, ym_cell, zm_cell)
  Cmpnts ***centcoor = NULL;
  // Pass the address of the pointer (centcoor) which has type Cmpnts****
  ierr = Allocate3DArray(&centcoor, zm_cell, ym_cell, xm_cell); CHKERRQ(ierr);
  LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "Allocated temporary array centcoor[%d][%d][%d] on rank %d.\n", zm_cell, ym_cell, xm_cell,rank);

  // Interpolate node coordinates to cell centers using the generic macro
  // This will dispatch to InterpolateFieldFromCornerToCenter_Vector
  ierr = InterpolateFieldFromCornerToCenter(coor_arr, centcoor, user); CHKERRQ(ierr);
  LOG_ALLOW(LOCAL, LOG_DEBUG, "Interpolated coordinates to cell centers on rank %d.\n",rank);
  

  // DEBUG -----------------------------------------------------
  // Get local node info from fda (still needed for interior check)
  /*
  DMDALocalInfo info_nodes;
  ierr = DMDAGetLocalInfo(user->fda, &info_nodes); CHKERRQ(ierr);

  // Get owned NODE ranges (we will loop over nodes directly now)
  PetscInt xs_node = info_nodes.xs;
  PetscInt xm_node = info_nodes.xm;
  PetscInt xe_node = xs_node + xm_node;
  PetscInt ys_node = info_nodes.ys;
  PetscInt ym_node = info_nodes.ym;
  PetscInt ye_node = ys_node + ym_node;
  PetscInt zs_node = info_nodes.zs;
  PetscInt zm_node = info_nodes.zm;
  PetscInt ze_node = zs_node + zm_node;
  */
  /* --- 4. Populate Target Field Vec --- */
  // --- DEBUG ---------------------------------------------------
  
  if (fieldIsVector) { // Handle Ucat or Ucont
    Cmpnts ***vecField_arr;
    
    ierr = DMDAVecGetArray(accessDM, fieldVec, &vecField_arr); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Obtained write access to vector field %s.\n", fieldName);
    
    // Loop over the GLOBAL indices of the OWNED CELLS
    for (PetscInt k = zs_cell; k < ze_cell; k++) { // Global CELL index k
      for (PetscInt j = ys_cell; j < ye_cell; j++) { // Global CELL index j
        for (PetscInt i = xs_cell; i < xe_cell; i++) { // Global CELL index i

            PetscInt k_local = k - zs_cell;
            PetscInt j_local = j - ys_cell;
            PetscInt i_local = i - xs_cell;

            // Check local bounds (defensive)
            if (k_local < 0 || k_local >= zm_cell || j_local < 0 || j_local >= ym_cell || i_local < 0 || i_local >= xm_cell) {
                  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Calculated local index out of bounds for centcoor!");
            }
            Cmpnts* cell_center_coord_ptr = &centcoor[k_local][j_local][i_local];

            // Determine the GLOBAL NODE index where this cell's value should be stored
            PetscInt node_k = k;
            PetscInt node_j = j;
            PetscInt node_i = i;

             // Only write if it's an INTERIOR node index
             if (node_i >= 1 && node_i < info_nodes.mx - 1 &&
                 node_j >= 1 && node_j < info_nodes.my - 1 &&
                 node_k >= 1 && node_k < info_nodes.mz - 1)
             {
                // Use the generic macro - dispatches to SetLocalCartesianField_Vector
	       ierr = SetLocalCartesianField(fieldName, &vecField_arr[node_k][node_j][node_i], cell_center_coord_ptr,user->FieldInitialization,ConstantValue); CHKERRQ(ierr);
             }
             // Else: Boundary node indices left untouched for FormBCS
        }
      }
    }
    
    /*  ---- DEBUG ------------------
    // --- MODIFIED LOOP: Loop over OWNED NODES ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "[DEBUG] Populating interior nodes of %s with constant values.\n", fieldName);
    for (PetscInt k = zs_node; k < ze_node; k++) { // Global NODE index k
      for (PetscInt j = ys_node; j < ye_node; j++) { // Global NODE index j
        for (PetscInt i = xs_node; i < xe_node; i++) { // Global NODE index i

             // Only write if it's an INTERIOR node index
             if (i >= 1 && i < info_nodes.mx - 1 &&
                 j >= 1 && j < info_nodes.my - 1 &&
                 k >= 1 && k < info_nodes.mz - 1)
             {
                 // Set a constant value, bypassing analytical calculation and centcoor
                 vecField_arr[k][j][i].x = 1.0; // Example constant
                 vecField_arr[k][j][i].y = 0.5;
                 vecField_arr[k][j][i].z = 0.1;
             }
             // Else: Boundary node indices left untouched
        }
      }
    }
    // --- END MODIFIED LOOP ---
    -------------- DEBUG ------------------
    */
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Populated interior node indices of %s with analytical values on rank %d.\n", fieldName,rank);
    ierr = DMDAVecRestoreArray(accessDM, fieldVec, &vecField_arr); CHKERRQ(ierr); // MOVE to BOTTOM if error!

  } else { // Handle P or Nvert (Scalar Fields)
    PetscReal ***scalarField_arr;
    // Use the correct DM (da) for scalar fields
    ierr = DMDAVecGetArray(accessDM, fieldVec, &scalarField_arr); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Obtained write access to scalar field %s.\n", fieldName);

    // Loop over the GLOBAL indices of the OWNED CELLS
    for (PetscInt k = zs_cell; k < ze_cell; k++) { // Global CELL index k
      for (PetscInt j = ys_cell; j < ye_cell; j++) { // Global CELL index j
        for (PetscInt i = xs_cell; i < xe_cell; i++) { // Global CELL index i

            PetscInt k_local = k - zs_cell;
            PetscInt j_local = j - ys_cell;
            PetscInt i_local = i - xs_cell;

             // Check local bounds (defensive)
            if (k_local < 0 || k_local >= zm_cell || j_local < 0 || j_local >= ym_cell || i_local < 0 || i_local >= xm_cell) {
                  SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Calculated local index out of bounds for centcoor!");
            }
            Cmpnts* cell_center_coord_ptr = &centcoor[k_local][j_local][i_local];

            // Use the generic macro - dispatches to SetLocalCartesianField_Scalar
            // Store using the GLOBAL CELL index [k][j][i]
            ierr = SetLocalCartesianField(fieldName, &scalarField_arr[k][j][i], cell_center_coord_ptr,user->FieldInitialization,ConstantValue); CHKERRQ(ierr);
        }
      }
    }
    

    /*  --------- DEBUG ---------------------
        // --- MODIFIED LOOP: Loop over OWNED CELLS (using node ranges for now, simpler) ---
    // For scalars, we actually want cell indices. Get them properly.
    PetscInt xs_cell, xm_cell, xe_cell;
    PetscInt ys_cell, ym_cell, ye_cell;
    PetscInt zs_cell, zm_cell, ze_cell;
    ierr = GetOwnedCellRange(&info_nodes, 0, &xs_cell, &xm_cell); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_nodes, 1, &ys_cell, &ym_cell); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_nodes, 2, &zs_cell, &zm_cell); CHKERRQ(ierr);
    xe_cell = xs_cell + xm_cell;
    ye_cell = ys_cell + ym_cell;
    ze_cell = zs_cell + zm_cell;

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "[DEBUG] Populating owned cells of %s with constant values.\n", fieldName);
    for (PetscInt k = zs_cell; k < ze_cell; k++) { // Global CELL index k
      for (PetscInt j = ys_cell; j < ye_cell; j++) { // Global CELL index j
        for (PetscInt i = xs_cell; i < xe_cell; i++) { // Global CELL index i
            // Set a constant value, bypassing analytical calculation and centcoor
            scalarField_arr[k][j][i] = 10.0; // Example constant
        }
      }
    }
    // --- END MODIFIED LOOP ---
    ---------- DEBUG ----------------------
    */
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Populated interior node indices of %s with analytical values.\n", fieldName);
    ierr = DMDAVecRestoreArray(accessDM, fieldVec, &scalarField_arr); CHKERRQ(ierr); // MOVE to BOTTOM if error!
  }

  /* --- 5. Cleanup --- */
  // Use the generic macro for deallocation - pass the pointer centcoor (type Cmpnts***)
  ierr = Deallocate3DArray(centcoor, zm_cell, ym_cell); CHKERRQ(ierr); // Use correct dimensions from GetOwnedCellRange
  ierr = DMDAVecRestoreArrayRead(user->fda, localCoor, &coor_arr); CHKERRQ(ierr); 

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Analytical Solution Setup Complete for Interior.\n");
  
  // LOG_ALLOW(GLOBAL, LOG_DEBUG, "[DEBUG] Skipped cleanup for bypassed coordinate steps.\n")
   MPI_Barrier(PETSC_COMM_WORLD); // Ensure rank 1 finishes test before proceeding
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

