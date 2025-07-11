// Boundaries.c

#include "Boundaries.h"

/**
 * @brief Determines if the current MPI rank owns any part of the globally defined inlet face,
 *        making it responsible for placing particles on that portion of the surface.
 *
 * The determination is based on the rank's owned nodes (from `DMDALocalInfo`) and
 * the global node counts, in conjunction with the `user->identifiedInletBCFace`.
 * A rank can service an inlet face if it owns the cells adjacent to that global boundary
 * and has a non-zero extent (owns cells) in the tangential dimensions of that face.
 *
 * @param user Pointer to the UserCtx structure, containing `identifiedInletBCFace`.
 * @param info Pointer to the DMDALocalInfo for the current rank's DA (node-based).
 * @param IM_nodes_global Global number of nodes in the I-direction (e.g., user->IM + 1 if user->IM is cell count).
 * @param JM_nodes_global Global number of nodes in the J-direction.
 * @param KM_nodes_global Global number of nodes in the K-direction.
 * @param[out] can_service_inlet_out Pointer to a PetscBool; set to PETSC_TRUE if the rank
 *                                   services (part of) the inlet, PETSC_FALSE otherwise.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode CanRankServiceInletFace(UserCtx *user, const DMDALocalInfo *info,
                                              PetscInt IM_nodes_global, PetscInt JM_nodes_global, PetscInt KM_nodes_global,
                                              PetscBool *can_service_inlet_out)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank_for_logging; // For detailed debugging logs
    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank_for_logging); CHKERRQ(ierr);

    *can_service_inlet_out = PETSC_FALSE; // Default to no service

    if (!user->inletFaceDefined) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d]: Inlet face not defined in user context. Cannot service.\n", rank_for_logging);
        PetscFunctionReturn(0);
    }

    // Get the range of cells owned by this rank in each dimension
    PetscInt owned_start_cell_i, num_owned_cells_on_rank_i;
    PetscInt owned_start_cell_j, num_owned_cells_on_rank_j;
    PetscInt owned_start_cell_k, num_owned_cells_on_rank_k;

    ierr = GetOwnedCellRange(info, 0, &owned_start_cell_i, &num_owned_cells_on_rank_i); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 1, &owned_start_cell_j, &num_owned_cells_on_rank_j); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 2, &owned_start_cell_k, &num_owned_cells_on_rank_k); CHKERRQ(ierr);

    // Determine the global index of the last cell (0-indexed) in each direction.
    // Example: If IM_nodes_global = 11 (nodes 0-10), there are 10 cells (0-9). Last cell index is 9.
    // Formula: global_nodes - 1 (num cells) - 1 (0-indexed) = global_nodes - 2.
    PetscInt last_global_cell_idx_i = (IM_nodes_global > 1) ? (IM_nodes_global - 2) : -1; // -1 if 0 or 1 node (i.e., 0 cells)
    PetscInt last_global_cell_idx_j = (JM_nodes_global > 1) ? (JM_nodes_global - 2) : -1;
    PetscInt last_global_cell_idx_k = (KM_nodes_global > 1) ? (KM_nodes_global - 2) : -1;

    switch (user->identifiedInletBCFace) {
        case BC_FACE_NEG_X: // Inlet on the global I-minimum face (face of cell C_i=0)
            // Rank services if its first owned node is global node 0 (info->xs == 0),
            // and it owns cells in I, J, and K directions.
            if (info->xs == 0 && num_owned_cells_on_rank_i > 0 &&
                num_owned_cells_on_rank_j > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_inlet_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_X: // Inlet on the global I-maximum face (face of cell C_i=last_global_cell_idx_i)
            // Rank services if it owns the last cell in I-direction,
            // and has extent in J and K.
            if (last_global_cell_idx_i >= 0 && /* Check for valid global domain */
                (owned_start_cell_i + num_owned_cells_on_rank_i - 1) == last_global_cell_idx_i && /* Rank's last cell is the global last cell */
                num_owned_cells_on_rank_j > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_inlet_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_NEG_Y:
            if (info->ys == 0 && num_owned_cells_on_rank_j > 0 &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_inlet_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_Y:
            if (last_global_cell_idx_j >= 0 &&
                (owned_start_cell_j + num_owned_cells_on_rank_j - 1) == last_global_cell_idx_j &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_inlet_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_NEG_Z:
            if (info->zs == 0 && num_owned_cells_on_rank_k > 0 &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_j > 0) {
                *can_service_inlet_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_Z:
            if (last_global_cell_idx_k >= 0 &&
                (owned_start_cell_k + num_owned_cells_on_rank_k - 1) == last_global_cell_idx_k &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_j > 0) {
                *can_service_inlet_out = PETSC_TRUE;
            }
            break;
        default:
             LOG_ALLOW(LOCAL, LOG_WARNING, "[Rank %d]: Unknown inlet face enum %d.\n", rank_for_logging, user->identifiedInletBCFace);
            break;
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] Inlet face enum %d. Owns cells (i,j,k):(%d,%d,%d) starting at cell (%d,%d,%d). Global nodes(I,J,K):(%d,%d,%d). ==> Can service: %s.\n",
        rank_for_logging, user->identifiedInletBCFace,
        num_owned_cells_on_rank_i, num_owned_cells_on_rank_j, num_owned_cells_on_rank_k,
        owned_start_cell_i, owned_start_cell_j, owned_start_cell_k,
        IM_nodes_global, JM_nodes_global, KM_nodes_global,
        (*can_service_inlet_out) ? "TRUE" : "FALSE");

    PetscFunctionReturn(0);
}

/**
 * @brief Determines if the current MPI rank owns any part of a specified global face.
 *
 * This function is a general utility for parallel boundary operations. It checks if the
 * local domain of the current MPI rank is adjacent to a specified global boundary face.
 * A rank "services" a face if it owns the cells adjacent to that global boundary and has
 * a non-zero extent (i.e., owns at least one cell) in the tangential dimensions of that face.
 *
 * @param info              Pointer to the DMDALocalInfo for the current rank's DA.
 * @param face_id           The specific global face (e.g., BC_FACE_NEG_Z) to check.
 * @param[out] can_service_out Pointer to a PetscBool; set to PETSC_TRUE if the rank
 *                           services the face, PETSC_FALSE otherwise.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode CanRankServiceFace(const DMDALocalInfo *info, BCFace face_id, PetscBool *can_service_out)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank_for_logging;
    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank_for_logging); CHKERRQ(ierr);

    *can_service_out = PETSC_FALSE; // Default to no service

    // Get the global dimensions (total number of nodes) from the DMDALocalInfo
    PetscInt IM_nodes_global = info->mx;
    PetscInt JM_nodes_global = info->my;
    PetscInt KM_nodes_global = info->mz;

    // Get the range of cells owned by this rank
    PetscInt owned_start_cell_i, num_owned_cells_on_rank_i;
    PetscInt owned_start_cell_j, num_owned_cells_on_rank_j;
    PetscInt owned_start_cell_k, num_owned_cells_on_rank_k;
    ierr = GetOwnedCellRange(info, 0, &owned_start_cell_i, &num_owned_cells_on_rank_i); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 1, &owned_start_cell_j, &num_owned_cells_on_rank_j); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 2, &owned_start_cell_k, &num_owned_cells_on_rank_k); CHKERRQ(ierr);

    // Determine the global index of the last cell (0-indexed) in each direction.
    PetscInt last_global_cell_idx_i = (IM_nodes_global > 1) ? (IM_nodes_global - 2) : -1;
    PetscInt last_global_cell_idx_j = (JM_nodes_global > 1) ? (JM_nodes_global - 2) : -1;
    PetscInt last_global_cell_idx_k = (KM_nodes_global > 1) ? (KM_nodes_global - 2) : -1;

    switch (face_id) {
        case BC_FACE_NEG_X:
            if (info->xs == 0 && num_owned_cells_on_rank_i > 0 &&
                num_owned_cells_on_rank_j > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_X:
            if (last_global_cell_idx_i >= 0 &&
                (owned_start_cell_i + num_owned_cells_on_rank_i - 1) == last_global_cell_idx_i &&
                num_owned_cells_on_rank_j > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_NEG_Y:
            if (info->ys == 0 && num_owned_cells_on_rank_j > 0 &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_Y:
            if (last_global_cell_idx_j >= 0 &&
                (owned_start_cell_j + num_owned_cells_on_rank_j - 1) == last_global_cell_idx_j &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_NEG_Z:
            if (info->zs == 0 && num_owned_cells_on_rank_k > 0 &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_j > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_Z:
            if (last_global_cell_idx_k >= 0 &&
                (owned_start_cell_k + num_owned_cells_on_rank_k - 1) == last_global_cell_idx_k &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_j > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        default:
             LOG_ALLOW(LOCAL, LOG_WARNING, "CanRankServiceFace - Rank %d: Unknown face enum %d. \n", rank_for_logging, face_id);
            break;
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "CanRankServiceFace - Rank %d check for face %d: Result=%s. \n",
        rank_for_logging, face_id, (*can_service_out ? "TRUE" : "FALSE"));

    PetscFunctionReturn(0);
}

/**
 * @brief Assuming the current rank services the inlet face, this function selects a random
 *        cell (owned by this rank on that face) and random logical coordinates within that cell,
 *        suitable for placing a particle on the inlet surface.
 *
 * It is the caller's responsibility to ensure CanRankServiceInletFace returned true.
 *
 * @param user Pointer to UserCtx.
 * @param info Pointer to DMDALocalInfo for the current rank (node-based).
 * @param xs_gnode, ys_gnode, zs_gnode Local starting node indices (incl. ghosts) for the rank's DA.
 * @param IM_nodes_global, JM_nodes_global, KM_nodes_global Global node counts.
 * @param rand_logic_i_ptr, rand_logic_j_ptr, rand_logic_k_ptr Pointers to RNGs for logical coords.
 * @param[out] ci_metric_lnode_out, cj_metric_lnode_out, ck_metric_lnode_out Local node indices of the selected cell's origin (these are local to the rank's DA including ghosts).
 * @param[out] xi_metric_logic_out, eta_metric_logic_out, zta_metric_logic_out Logical coords [0,1] within the cell.
 * @return PetscErrorCode
 */
PetscErrorCode GetRandomCellAndLogicOnInletFace(
    UserCtx *user, const DMDALocalInfo *info,
    PetscInt xs_gnode_rank, PetscInt ys_gnode_rank, PetscInt zs_gnode_rank, // Local starting node index (with ghosts) of the rank's DA patch
    PetscInt IM_nodes_global, PetscInt JM_nodes_global, PetscInt KM_nodes_global,
    PetscRandom *rand_logic_i_ptr, PetscRandom *rand_logic_j_ptr, PetscRandom *rand_logic_k_ptr,
    PetscInt *ci_metric_lnode_out, PetscInt *cj_metric_lnode_out, PetscInt *ck_metric_lnode_out,
    PetscReal *xi_metric_logic_out, PetscReal *eta_metric_logic_out, PetscReal *zta_metric_logic_out)
{
    PetscErrorCode ierr = 0;
    PetscReal r_val_i_sel, r_val_j_sel, r_val_k_sel;
    PetscInt local_cell_idx_on_face_dim1 = 0; // 0-indexed relative to owned cells on face
    PetscInt local_cell_idx_on_face_dim2 = 0;
    PetscMPIInt rank_for_logging; 

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank_for_logging); CHKERRQ(ierr); 

    // Defaults for cell origin node (local index for the rank's DA patch, including ghosts)
    *ci_metric_lnode_out = xs_gnode_rank; *cj_metric_lnode_out = ys_gnode_rank; *ck_metric_lnode_out = zs_gnode_rank;
    // Defaults for logical coordinates
    *xi_metric_logic_out = 0.5; *eta_metric_logic_out = 0.5; *zta_metric_logic_out = 0.5;

    // Get number of cells this rank owns in each dimension (tangential to the face mainly)
    PetscInt owned_start_cell_i, num_owned_cells_on_rank_i;
    PetscInt owned_start_cell_j, num_owned_cells_on_rank_j;
    PetscInt owned_start_cell_k, num_owned_cells_on_rank_k;

    ierr = GetOwnedCellRange(info, 0, &owned_start_cell_i, &num_owned_cells_on_rank_i); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 1, &owned_start_cell_j, &num_owned_cells_on_rank_j); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 2, &owned_start_cell_k, &num_owned_cells_on_rank_k); CHKERRQ(ierr);

    // Index of the last cell (0-indexed) in each global direction
    PetscInt last_global_cell_idx_i = (IM_nodes_global > 1) ? (IM_nodes_global - 2) : -1;
    PetscInt last_global_cell_idx_j = (JM_nodes_global > 1) ? (JM_nodes_global - 2) : -1;
    PetscInt last_global_cell_idx_k = (KM_nodes_global > 1) ? (KM_nodes_global - 2) : -1;
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Inlet face %d. Owned cells(i,j,k):(%d,%d,%d). GlobNodes(I,J,K):(%d,%d,%d). Rank's DA node starts (xs_g,ys_g,zs_g): (%d,%d,%d).\n",
        rank_for_logging, user->identifiedInletBCFace, num_owned_cells_on_rank_i,num_owned_cells_on_rank_j,num_owned_cells_on_rank_k,
        IM_nodes_global,JM_nodes_global,KM_nodes_global, xs_gnode_rank,ys_gnode_rank,zs_gnode_rank);


    switch (user->identifiedInletBCFace) {
        case BC_FACE_NEG_X: // Particle on -X face of cell C_0 (origin node N_0)
            // Cell origin node is the first owned node in I by this rank (global index info->xs).
            // Its local index within the rank's DA (incl ghosts) is xs_gnode_rank.
            *ci_metric_lnode_out = xs_gnode_rank;
            *xi_metric_logic_out = 1.0e-6;

            // Tangential dimensions are J and K. Select an owned cell randomly on this face.
            // num_owned_cells_on_rank_j/k must be > 0 (checked by CanRankServiceInletFace)
            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val_j_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim1 = (PetscInt)(r_val_j_sel * num_owned_cells_on_rank_j); // Index among owned J-cells
            local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_on_rank_j - 1);
            *cj_metric_lnode_out = ys_gnode_rank + local_cell_idx_on_face_dim1; // Offset from start of rank's J-nodes

            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val_k_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim2 = (PetscInt)(r_val_k_sel * num_owned_cells_on_rank_k);
            local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_on_rank_k - 1);
            *ck_metric_lnode_out = zs_gnode_rank + local_cell_idx_on_face_dim2;

            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
            break;

        case BC_FACE_POS_X: // Particle on +X face of cell C_last_I (origin node N_last_I_origin)
            // Origin node of the last I-cell is global_node_idx = last_global_cell_idx_i.
            // Its local index in rank's DA: (last_global_cell_idx_i - info->xs) + xs_gnode_rank
            *ci_metric_lnode_out = xs_gnode_rank + (last_global_cell_idx_i - info->xs);
            *xi_metric_logic_out = 1.0 - 1.0e-6;

            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val_j_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim1 = (PetscInt)(r_val_j_sel * num_owned_cells_on_rank_j);
            local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_on_rank_j - 1);
            *cj_metric_lnode_out = ys_gnode_rank + local_cell_idx_on_face_dim1;

            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val_k_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim2 = (PetscInt)(r_val_k_sel * num_owned_cells_on_rank_k);
            local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_on_rank_k - 1);
            *ck_metric_lnode_out = zs_gnode_rank + local_cell_idx_on_face_dim2;

            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
            break;
        // ... (Cases for Y and Z faces, following the same pattern) ...
        case BC_FACE_NEG_Y:
            *cj_metric_lnode_out = ys_gnode_rank;
            *eta_metric_logic_out = 1.0e-6;
            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val_i_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim1 = (PetscInt)(r_val_i_sel * num_owned_cells_on_rank_i);
            local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_on_rank_i - 1);
            *ci_metric_lnode_out = xs_gnode_rank + local_cell_idx_on_face_dim1;
            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val_k_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim2 = (PetscInt)(r_val_k_sel * num_owned_cells_on_rank_k);
            local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_on_rank_k - 1);
            *ck_metric_lnode_out = zs_gnode_rank + local_cell_idx_on_face_dim2;
            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr);
            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
            break;
        case BC_FACE_POS_Y:
            *cj_metric_lnode_out = ys_gnode_rank + (last_global_cell_idx_j - info->ys);
            *eta_metric_logic_out = 1.0 - 1.0e-6;
            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val_i_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim1 = (PetscInt)(r_val_i_sel * num_owned_cells_on_rank_i);
            local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_on_rank_i - 1);
            *ci_metric_lnode_out = xs_gnode_rank + local_cell_idx_on_face_dim1;
            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val_k_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim2 = (PetscInt)(r_val_k_sel * num_owned_cells_on_rank_k);
            local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_on_rank_k - 1);
            *ck_metric_lnode_out = zs_gnode_rank + local_cell_idx_on_face_dim2;
            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr);
            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
            break;
        case BC_FACE_NEG_Z: // Your example case
            *ck_metric_lnode_out = zs_gnode_rank; // Cell origin is the first owned node in K by this rank
            *zta_metric_logic_out = 1.0e-6;      // Place particle slightly inside this cell from its -Z face
            // Tangential dimensions are I and J
            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val_i_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim1 = (PetscInt)(r_val_i_sel * num_owned_cells_on_rank_i);
            local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_on_rank_i - 1);
            *ci_metric_lnode_out = xs_gnode_rank + local_cell_idx_on_face_dim1;

            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val_j_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim2 = (PetscInt)(r_val_j_sel * num_owned_cells_on_rank_j);
            local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_on_rank_j - 1);
            *cj_metric_lnode_out = ys_gnode_rank + local_cell_idx_on_face_dim2;

            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr); // Intra-cell logical for I
            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr); // Intra-cell logical for J
            break;
        case BC_FACE_POS_Z:
            *ck_metric_lnode_out = zs_gnode_rank + (last_global_cell_idx_k - info->zs);
            *zta_metric_logic_out = 1.0 - 1.0e-6;
            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val_i_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim1 = (PetscInt)(r_val_i_sel * num_owned_cells_on_rank_i);
            local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_on_rank_i - 1);
            *ci_metric_lnode_out = xs_gnode_rank + local_cell_idx_on_face_dim1;
            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val_j_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim2 = (PetscInt)(r_val_j_sel * num_owned_cells_on_rank_j);
            local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_on_rank_j - 1);
            *cj_metric_lnode_out = ys_gnode_rank + local_cell_idx_on_face_dim2;
            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr);
            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
            break;
        default:
             SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "GetRandomCellAndLogicOnInletFace: Invalid user->identifiedInletBCFace %d. \n", user->identifiedInletBCFace);
    }

    PetscReal eps = 1.0e-7;
    if (user->identifiedInletBCFace == BC_FACE_NEG_X || user->identifiedInletBCFace == BC_FACE_POS_X) {
        *eta_metric_logic_out = PetscMin(PetscMax(0.0, *eta_metric_logic_out), 1.0 - eps);
        *zta_metric_logic_out = PetscMin(PetscMax(0.0, *zta_metric_logic_out), 1.0 - eps);
    } else if (user->identifiedInletBCFace == BC_FACE_NEG_Y || user->identifiedInletBCFace == BC_FACE_POS_Y) {
        *xi_metric_logic_out  = PetscMin(PetscMax(0.0, *xi_metric_logic_out),  1.0 - eps);
        *zta_metric_logic_out = PetscMin(PetscMax(0.0, *zta_metric_logic_out), 1.0 - eps);
    } else { 
        *xi_metric_logic_out  = PetscMin(PetscMax(0.0, *xi_metric_logic_out),  1.0 - eps);
        *eta_metric_logic_out = PetscMin(PetscMax(0.0, *eta_metric_logic_out), 1.0 - eps);
    }
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Target CellNode(loc lnode idx)=(%d,%d,%d). Logic(xi,et,zt)=(%.2e,%.2f,%.2f). \n",
        rank_for_logging, *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out,
        *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);

    PetscFunctionReturn(0);
}

/**
 * @brief (Private) Creates and configures a specific BoundaryCondition handler object.
 *
 * This function acts as a factory. Based on the requested handler_type, it allocates
 * a BoundaryCondition object and populates it with the correct set of function
 * pointers corresponding to that specific behavior.
 *
 * @param handler_type The specific handler to create (e.g., BC_HANDLER_WALL_NOSLIP).
 * @param[out] new_bc_ptr  A pointer to where the newly created BoundaryCondition
 *                         object's address will be stored.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode BoundaryCondition_Create(BCHandlerType handler_type, BoundaryCondition **new_bc_ptr)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    const char* handler_name = BCHandlerTypeToString(handler_type);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Factory called for handler type %d (%s). \n", handler_type, handler_name);

    ierr = PetscMalloc1(1, new_bc_ptr); CHKERRQ(ierr);
    BoundaryCondition *bc = *new_bc_ptr;

    bc->type        = handler_type;
    bc->data        = NULL;
    bc->Initialize  = NULL;
    bc->PreStep     = NULL;
    bc->Apply       = NULL;
    bc->PlaceSource = NULL;
    bc->Destroy     = NULL;

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Allocated generic handler object at address %p.\n", (void*)bc);

    switch (handler_type) {

        case BC_HANDLER_NOGRAD_COPY_GHOST:
	     LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_UniformFlowCopyGhost().\n");
             ierr = Create_NogradCopyGhost(bc); CHKERRQ(ierr);
             break;
	     
        case BC_HANDLER_WALL_NOSLIP:
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_WallNoSlip().\n");
            ierr = Create_WallNoSlip(bc); CHKERRQ(ierr);
            break;

        case BC_HANDLER_INLET_CONSTANT_VELOCITY:
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_InletConstantVelocity().\n");
	    ierr = Create_InletConstantVelocity(bc); CHKERRQ(ierr);
            break;

    case BC_HANDLER_INLET_PARABOLIC:
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_InletParabolicProfile().\n");
	    ierr = Create_InletParabolicProfile(bc); CHKERRQ(ierr);
            break;
        /* Add cases for other handlers here in future phases */
        
        default:
            LOG_ALLOW(GLOBAL, LOG_ERROR, "Handler type %d (%s) is not recognized or implemented in the factory.\n", handler_type, handler_name);
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Boundary handler type %d (%s) not recognized in factory.\n", handler_type, handler_name);
    }
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "BoundaryCondition_Create: Successfully created and configured handler for %s.\n", handler_name);
    PetscFunctionReturn(0);
}


//================================================================================
//
//                       PUBLIC MASTER SETUP FUNCTION
//
//================================================================================

/**
 * @brief Initializes the entire boundary system based on a configuration file.
 *
 * (Full Doxygen from header file goes here)
 */
PetscErrorCode BoundarySystem_Create(UserCtx *user, const char *bcs_filename)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting creation and initialization of all boundary handlers.\n");

    // Step 1: Parse the configuration file to determine user intent.
    // This function, defined in io.c, populates the configuration enums and parameter
    // lists within the user->boundary_faces array on all MPI ranks.
    ierr = ParseAllBoundaryConditions(user, bcs_filename); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Configuration file '%s' parsed successfully.\n", bcs_filename);

    // Step 2: Create and Initialize the handler object for each of the 6 faces.
    for (int i = 0; i < 6; i++) {
        BoundaryFaceConfig *face_cfg = &user->boundary_faces[i];
        
        const char *face_name = BCFaceToString(face_cfg->face_id);
        const char *handler_name = BCHandlerTypeToString(face_cfg->handler_type);

        LOG_ALLOW(LOCAL, LOG_DEBUG, "Creating handler for Face %d (%s) with handler type '%s'.\n", i, face_name, handler_name);

        // Use the private factory to construct the correct handler object based on the parsed type.
        // The factory returns a pointer to the new handler object, which we store in the config struct.
        ierr = BoundaryCondition_Create(face_cfg->handler_type, &face_cfg->handler); CHKERRQ(ierr);

        // Step 3: Call the specific Initialize() method for the newly created handler.
        // This allows the handler to perform its own setup, like reading parameters from the
        // face_cfg->params list and setting the initial field values on its face.
        if (face_cfg->handler && face_cfg->handler->Initialize) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Calling Initialize() method for handler on Face %d (%s).\n", i, face_name);
            
            // Prepare the context needed by the Initialize() function.
            BCContext ctx = {
                .user = user,
                .face_id = face_cfg->face_id,
                .global_inflow_sum = NULL,  // Global flux sums are not relevant during initialization.
                .global_outflow_sum = NULL
            };
            
            ierr = face_cfg->handler->Initialize(face_cfg->handler, &ctx); CHKERRQ(ierr);
        } else {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Handler for Face %d (%s) has no Initialize() method, skipping.\n", i, face_name);
        }
    }
    // ====================================================================================
    // --- NEW: Step 4: Synchronize Vectors After Initialization ---
    // This is the CRITICAL fix. The Initialize() calls have modified local vector
    // arrays on some ranks but not others. We must now update the global vector state
    // and then update all local ghost regions to be consistent.
    // ====================================================================================
     
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Committing global boundary initializations to local vectors.\n");

    // Commit changes from the global vectors (Ucat, Ucont) to the local vectors (lUcat, lUcont)
    // NOTE: The Apply functions modified Ucat and Ucont via GetArray, which works on the global
    // representation.
    ierr = DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat); CHKERRQ(ierr);
    
    ierr = DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont); CHKERRQ(ierr);
    /*
     // Now, update all local vectors (including ghost cells) from the newly consistent global vectors

    ierr = DMLocalToGlobalBegin(user->fda, user->lUcat, INSERT_VALUES, user->Ucat); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(user->fda, user->lUcat, INSERT_VALUES, user->Ucat); CHKERRQ(ierr);
    
    ierr = DMLocalToGlobalBegin(user->fda, user->lUcont, INSERT_VALUES, user->Ucont); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(user->fda, user->lUcont, INSERT_VALUES, user->Ucont); CHKERRQ(ierr);
    */
    

    LOG_ALLOW(GLOBAL, LOG_INFO, "All boundary handlers created and initialized successfully.\n");
    PetscFunctionReturn(0);
}

//================================================================================
//
//                      PUBLIC MASTER TIME-STEP FUNCTION
//
//================================================================================

/**
 * @brief Executes one full boundary condition update cycle for a time step.
 *
 * This function is the main entry point from the time-stepping loop. It orchestrates
 * the three-phase process required for robustly applying boundary conditions in parallel:
 *   1. PreStep: Each handler is called to measure its local contribution to the
 *      domain's overall flux balance (e.g., target inflow, measured outflow).
 *   2. Reduction: A parallel sum (`MPI_Allreduce`) is performed to get the global,
 *      domain-wide totals for inflow and outflow.
 *   3. Apply: Each handler is called again. Now, they have access to the global
 *      flux data and can apply their conditions correctly (e.g., an outlet can
 *      now apply a correction to enforce mass conservation).
 *
 * @param user The main UserCtx struct, containing the live boundary handlers.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode BoundarySystem_ExecuteStep(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    
    // Time-step-scoped variables to hold the flux calculations for this step.
    PetscReal local_inflow_sum = 0.0;
    PetscReal local_outflow_sum = 0.0;
    PetscReal global_inflow_sum = 0.0;
    PetscReal global_outflow_sum = 0.0;

    LOG_ALLOW(LOCAL, LOG_DEBUG, "BoundarySystem: Executing PreStep phase for all handlers.\n");
    
    // --- PHASE 1: PRE-STEP (Measure local fluxes) ---
    // Each rank calculates its contribution without knowledge of other ranks.
    for (int i = 0; i < 6; i++) {
        BoundaryFaceConfig *face_cfg = &user->boundary_faces[i];
        if (face_cfg->handler && face_cfg->handler->PreStep) {
            // Prepare a minimal context for the PreStep phase.
            BCContext ctx = { .user = user, .face_id = face_cfg->face_id, .global_inflow_sum = NULL, .global_outflow_sum = NULL };
            ierr = face_cfg->handler->PreStep(face_cfg->handler, &ctx, &local_inflow_sum, &local_outflow_sum); CHKERRQ(ierr);
        }
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "PreStep complete. Local contributions: Inflow=%.4e, Outflow=%.4e \n", local_inflow_sum, local_outflow_sum);

    // --- PHASE 2: PARALLEL REDUCTION (Synchronize global fluxes) ---
    // This is the only point of global communication in the BC step. It ensures
    // every process has the same, correct totals.
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Performing MPI_Allreduce to synchronize fluxes.\n");
    ierr = MPI_Allreduce(&local_inflow_sum, &global_inflow_sum, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Allreduce(&local_outflow_sum, &global_outflow_sum, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "MPI_Allreduce complete. Global sums: Inflow=%.4e, Outflow=%.4e \n", global_inflow_sum, global_outflow_sum);

    // --- PHASE 3: APPLY (Update boundary values using global information) ---
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Executing Apply phase for all handlers.\n");
    
    // Now, prepare the context with pointers to the synchronized global data.
    // Any handler that needs to enforce conservation can now access these values.
    BCContext apply_ctx = {
        .user = user,
        .global_inflow_sum = &global_inflow_sum,
        .global_outflow_sum = &global_outflow_sum
    };

    for (int i = 0; i < 6; i++) {
        BoundaryFaceConfig *face_cfg = &user->boundary_faces[i];
        if (face_cfg->handler && face_cfg->handler->Apply) {
            apply_ctx.face_id = face_cfg->face_id; // Update the face_id for each call
            ierr = face_cfg->handler->Apply(face_cfg->handler, &apply_ctx); CHKERRQ(ierr);
        }
    }
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Apply phase complete. \n");

    // =========================================================================
    // ---PHASE 4: COMMIT GLOBAL CHANGES TO LOCAL VECTORS ---
    // =========================================================================

    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Committing global boundary changes to local vectors.\n");
      
    ierr = DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat); CHKERRQ(ierr);
    
    ierr = DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont); CHKERRQ(ierr);

    /*     
    // Commit changes from lUcat to Ucat
     ierr = DMLocalToGlobalBegin(user->fda, user->lUcat, INSERT_VALUES, user->Ucat); CHKERRQ(ierr);
     ierr = DMLocalToGlobalEnd(user->fda, user->lUcat, INSERT_VALUES, user->Ucat); CHKERRQ(ierr);
    
    // Commit changes from lUcont to Ucont
     ierr = DMLocalToGlobalBegin(user->fda, user->lUcont, INSERT_VALUES, user->Ucont); CHKERRQ(ierr);
     ierr = DMLocalToGlobalEnd(user->fda, user->lUcont, INSERT_VALUES, user->Ucont); CHKERRQ(ierr);

     LOG_ALLOW(GLOBAL, LOG_INFO, "changes for Ucat and Ucont committed to global and local states.\n");
    // =========================================================================
    */
    PetscFunctionReturn(0);
}


//================================================================================
//
//                         PUBLIC MASTER CLEANUP FUNCTION
//
//================================================================================

/**
 * @brief Cleans up and destroys all resources allocated by the boundary system.
 *
 * This function should be called once at the end of the simulation. It iterates
 * through all created handlers and calls their respective Destroy methods to free
 * any privately allocated data (like parameter lists or handler-specific data),
 * and then frees the handler object itself. This prevents memory leaks.
 *
 * @param user The main UserCtx struct containing the boundary system to be destroyed.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode BoundarySystem_Destroy(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting destruction of all boundary handlers. \n");

    for (int i = 0; i < 6; i++) {
        BoundaryFaceConfig *face_cfg = &user->boundary_faces[i];
        const char *face_name = BCFaceToString(face_cfg->face_id);

        // --- Step 1: Free the parameter linked list associated with this face ---
        if (face_cfg->params) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "  Freeing parameter list for Face %d (%s). \n", i, face_name);
            FreeBC_ParamList(face_cfg->params);
            face_cfg->params = NULL; // Good practice to nullify dangling pointers
        }

        // --- Step 2: Destroy the handler object itself ---
        if (face_cfg->handler) {
            const char *handler_name = BCHandlerTypeToString(face_cfg->handler->type);
            LOG_ALLOW(LOCAL, LOG_DEBUG, "  Destroying handler '%s' on Face %d (%s).\n", handler_name, i, face_name);
            
            // Call the handler's specific cleanup function first, if it exists.
            // This will free any memory stored in the handler's private `data` pointer.
            if (face_cfg->handler->Destroy) {
                ierr = face_cfg->handler->Destroy(face_cfg->handler); CHKERRQ(ierr);
            }

            // Finally, free the generic BoundaryCondition object itself.
            ierr = PetscFree(face_cfg->handler); CHKERRQ(ierr);
            face_cfg->handler = NULL;
        }
    }
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Destruction complete.\n");
    PetscFunctionReturn(0);
}
