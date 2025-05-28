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
        LOG_ALLOW(LOCAL, LOG_DEBUG, "CanRankServiceInletFace - Rank %d: Inlet face not defined in user context. Cannot service.", rank_for_logging);
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
        case BC_FACE_NEG_Z: // Your example case
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
             LOG_ALLOW(LOCAL, LOG_WARNING, "CanRankServiceInletFace - Rank %d: Unknown inlet face enum %d.", rank_for_logging, user->identifiedInletBCFace);
            break;
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "CanRankServiceInletFace - Rank %d: Inlet face enum %d. Owns cells (i,j,k):(%d,%d,%d) starting at cell (%d,%d,%d). Global nodes(I,J,K):(%d,%d,%d). ==> Can service: %s.",
        rank_for_logging, user->identifiedInletBCFace,
        num_owned_cells_on_rank_i, num_owned_cells_on_rank_j, num_owned_cells_on_rank_k,
        owned_start_cell_i, owned_start_cell_j, owned_start_cell_k,
        IM_nodes_global, JM_nodes_global, KM_nodes_global,
        (*can_service_inlet_out) ? "TRUE" : "FALSE");

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
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Inlet face %d. Owned cells(i,j,k):(%d,%d,%d). GlobNodes(I,J,K):(%d,%d,%d). Rank's DA node starts (xs_g,ys_g,zs_g): (%d,%d,%d).",
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
             SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "GetRandomCellAndLogicOnInletFace: Invalid user->identifiedInletBCFace %d", user->identifiedInletBCFace);
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
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Target CellNode(loc lnode idx)=(%d,%d,%d). Logic(xi,et,zt)=(%.2e,%.2f,%.2f)",
        rank_for_logging, *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out,
        *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);

    PetscFunctionReturn(0);
}
