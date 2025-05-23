// ParticleSwarm.c

#include "ParticleSwarm.h"

#define INTERPOLATION_DISTANCE_TOLERANCE 1.0e-14
/**
 * @brief Initializes the DMSwarm object within the UserCtx structure.
 *
 * This function creates the DMSwarm, sets its type and dimension, and configures basic swarm properties.
 *
 * @param[in,out] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeSwarm(UserCtx* user) {
    PetscErrorCode ierr;  // Error code for PETSc functions

    // Create the DMSwarm object for particle management
    ierr = DMCreate(PETSC_COMM_WORLD, &user->swarm); CHKERRQ(ierr);
    ierr = DMSetType(user->swarm, DMSWARM); CHKERRQ(ierr);
    ierr = DMSetDimension(user->swarm, 3); CHKERRQ(ierr);
    ierr = DMSwarmSetType(user->swarm, DMSWARM_BASIC); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO, "InitializeSwarm - DMSwarm created and configured.\n");

    return 0;
}

/**
 * @brief Registers a swarm field without finalizing registration.
 *
 * This function calls DMSwarmRegisterPetscDatatypeField for the given field,
 * but does not finalize the registration. The finalization is deferred until
 * all fields have been registered.
 *
 * @param swarm      [in]  The DMSwarm object.
 * @param fieldName  [in]  Name of the field to register.
 * @param fieldDim   [in]  Dimension of the field (1 for scalar, 3 for vector, etc.).
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
static PetscErrorCode RegisterSwarmField(DM swarm, const char *fieldName, PetscInt fieldDim)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    
    ierr = DMSwarmRegisterPetscDatatypeField(swarm, fieldName, fieldDim, PETSC_REAL); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"RegisterSwarmField - Registered field '%s' with dimension=%d.\n",
              fieldName, fieldDim);
    
    PetscFunctionReturn(0);
}


/**
 * @brief Registers necessary particle fields within the DMSwarm.
 *
 * This function registers fields such as position, velocity, CellID, and weight for each particle.
 *
 * @param[in,out] swarm   The DMSwarm object managing the particle swarm.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */

PetscErrorCode RegisterParticleFields(DM swarm)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    
    // Register each field using the helper function
    ierr = RegisterSwarmField(swarm, "position", 3); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG, "RegisterParticleFields - Registered field 'position'.\n");
    
    ierr = RegisterSwarmField(swarm, "velocity", 3); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"RegisterParticleFields - Registered field 'velocity'.\n");
    
    ierr = RegisterSwarmField(swarm, "DMSwarm_CellID", 3); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"RegisterParticleFields - Registered field 'DMSwarm_CellID'.\n");
    
    ierr = RegisterSwarmField(swarm, "weight", 3); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"RegisterParticleFields - Registered field 'weight'.\n");

    ierr = RegisterSwarmField(swarm,"P", 1); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"RegisterParticleFields - Registered field 'P' - Pressure.\n");

    ierr = RegisterSwarmField(swarm,"pos_phy", 3); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"RegisterParticleFields - Registered field 'pos_phy' - Physical Position.\n");
    
    // Finalize the field registration after all fields have been added
    ierr = DMSwarmFinalizeFieldRegister(swarm); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO,"RegisterParticleFields - Finalized field registration.\n");
    
    PetscFunctionReturn(0);
}


/**
 * @brief Initializes random number generators for assigning particle properties.
 *
 * This function creates and configures separate PETSc random number generators for the x, y, and z coordinates.
 *
 * @param[in,out] user    Pointer to the UserCtx structure containing simulation context.
 * @param[out]    randx    Pointer to store the RNG for the x-coordinate.
 * @param[out]    randy    Pointer to store the RNG for the y-coordinate.
 * @param[out]    randz    Pointer to store the RNG for the z-coordinate.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeRandomGenerators(UserCtx* user, PetscRandom *randx, PetscRandom *randy, PetscRandom *randz) {
    PetscErrorCode ierr;  // Error code for PETSc functions
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Initialize RNG for x-coordinate
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, randx); CHKERRQ(ierr);
    ierr = PetscRandomSetType((*randx), PETSCRAND48); CHKERRQ(ierr);
    ierr = PetscRandomSetInterval(*randx, user->bbox.min_coords.x, user->bbox.max_coords.x); CHKERRQ(ierr);
    ierr = PetscRandomSetSeed(*randx, rank + 12345); CHKERRQ(ierr);  // Unique seed per rank
    ierr = PetscRandomSeed(*randx); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOCAL,LOG_DEBUG, "InitializeRandomGenerators - Initialized RNG for X-axis.\n");

    // Initialize RNG for y-coordinate
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, randy); CHKERRQ(ierr);
    ierr = PetscRandomSetType((*randy), PETSCRAND48); CHKERRQ(ierr);
    ierr = PetscRandomSetInterval(*randy, user->bbox.min_coords.y, user->bbox.max_coords.y); CHKERRQ(ierr);
    ierr = PetscRandomSetSeed(*randy, rank + 67890); CHKERRQ(ierr);  // Unique seed per rank
    ierr = PetscRandomSeed(*randy); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOCAL,LOG_DEBUG, "InitializeRandomGenerators - Initialized RNG for Y-axis.\n");

    // Initialize RNG for z-coordinate
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, randz); CHKERRQ(ierr);
    ierr = PetscRandomSetType((*randz), PETSCRAND48); CHKERRQ(ierr);
    ierr = PetscRandomSetInterval(*randz, user->bbox.min_coords.z, user->bbox.max_coords.z); CHKERRQ(ierr);
    ierr = PetscRandomSetSeed(*randz, rank + 54321); CHKERRQ(ierr);  // Unique seed per rank
    ierr = PetscRandomSeed(*randz); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOCAL,LOG_DEBUG, "InitializeRandomGenerators - Initialized RNG for Z-axis.\n");

    return 0;
}

/**
 * @brief Initializes random number generators for logical space operations [0.0, 1.0).
 *
 * This function creates and configures three separate PETSc random number generators,
 * one for each logical dimension (i, j, k or xi, eta, zeta equivalent).
 * Each RNG is configured to produce uniformly distributed real numbers in the interval [0.0, 1.0).
 * These are typically used for selecting owned cells or generating intra-cell logical coordinates.
 *
 * @param[out]   rand_logic_i Pointer to store the RNG for the i-logical dimension.
 * @param[out]   rand_logic_j Pointer to store the RNG for the j-logical dimension.
 * @param[out]   rand_logic_k Pointer to store the RNG for the k-logical dimension.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeLogicalSpaceRNGs(PetscRandom *rand_logic_i, PetscRandom *rand_logic_j, PetscRandom *rand_logic_k) {
    PetscErrorCode ierr;
    PetscMPIInt rank;
    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // --- RNG for i-logical dimension ---
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, rand_logic_i); CHKERRQ(ierr);
    ierr = PetscRandomSetType((*rand_logic_i), PETSCRAND48); CHKERRQ(ierr);
    ierr = PetscRandomSetInterval(*rand_logic_i, 0.0, 1.0); CHKERRQ(ierr); // Key change: [0,1)
    ierr = PetscRandomSetSeed(*rand_logic_i, rank + 202401); CHKERRQ(ierr); // Unique seed
    ierr = PetscRandomSeed(*rand_logic_i); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOCAL,LOG_DEBUG, "InitializeLogicalSpaceRNGs - Initialized RNG for i-logical dimension [0,1).\n");

    // --- RNG for j-logical dimension ---
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, rand_logic_j); CHKERRQ(ierr);
    ierr = PetscRandomSetType((*rand_logic_j), PETSCRAND48); CHKERRQ(ierr);
    ierr = PetscRandomSetInterval(*rand_logic_j, 0.0, 1.0); CHKERRQ(ierr); // Key change: [0,1)
    ierr = PetscRandomSetSeed(*rand_logic_j, rank + 202402); CHKERRQ(ierr);
    ierr = PetscRandomSeed(*rand_logic_j); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOCAL,LOG_DEBUG, "InitializeLogicalSpaceRNGs - Initialized RNG for j-logical dimension [0,1).\n");

    // --- RNG for k-logical dimension ---
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, rand_logic_k); CHKERRQ(ierr);
    ierr = PetscRandomSetType((*rand_logic_k), PETSCRAND48); CHKERRQ(ierr);
    ierr = PetscRandomSetInterval(*rand_logic_k, 0.0, 1.0); CHKERRQ(ierr); // Key change: [0,1)
    ierr = PetscRandomSetSeed(*rand_logic_k, rank + 202403); CHKERRQ(ierr);
    ierr = PetscRandomSeed(*rand_logic_k); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOCAL,LOG_DEBUG, "InitializeLogicalSpaceRNGs - Initialized RNG for k-logical dimension [0,1).\n");

    PetscFunctionReturn(0);
}

/**
 * @brief Determines cell selection and intra-cell logical coordinates for surface initialization (Mode 0).
 *
 * This function is called when `user->ParticleInitialization == 0`. Based on the
 * `user->identifiedInletBCFace` (which is determined by parsing "bcs.dat"), this function
 * calculates parameters for placing a particle on the specified global inlet face.
 *
 * The process involves:
 * 1. Checking if the current MPI rank owns any portion of the designated global inlet face.
 * 2. If it does, it randomly selects an *owned cell* that lies on this rank's portion of the inlet face.
 *    The selection uses the number of owned cells in the two dimensions tangential to the face.
 * 3. It sets the intra-cell logical coordinate normal to the inlet face to a very small value
 *    (e.g., 1e-6 for a MIN face like BC_FACE_NEG_X, or 1.0 - 1.0e-6 for a MAX face like BC_FACE_POS_X),
 *    ensuring the particle is just inside the chosen cell, on its "inlet-facing" logical surface.
 * 4. The other two intra-cell logical coordinates (tangential to the inlet face) are chosen randomly
 *    within the `[0,1)` range.
 * 5. The outputs are the local node indices (`ci/cj/ck_metric_lnode_out`) of the origin of the
 *    selected cell (for use with `MetricLogicalToPhysical`) and the calculated intra-cell
 *    logical coordinates (`xi/eta/zta_metric_logic_out`).
 * 6. A flag (`can_place_on_surface_out`) indicates if a valid placement could be determined by this rank.
 *
 * **Important Note on `DMDALocalInfo info` members:**
 *   - `info->xs, info->ys, info->zs`: Global starting indices of *owned cells*.
 *   - `info->mx, info->my, info->mz`: Number of *grid points (nodes)* in each local dimension on this process.
 *     Therefore, the number of *owned cells* in a dimension is `info->mX - 1` (if `info->mX > 0`).
 *
 * @param[in]  user Pointer to `UserCtx` (contains `identifiedInletBCFace`).
 * @param[in]  info Pointer to `DMDALocalInfo` for the current rank's grid portion.
 * @param[in]  xs_gnode, ys_gnode, zs_gnode Local indices (in the ghosted array) of the first *owned node*.
 * @param[in]  IM_gcells, JM_gcells, KM_gcells Total number of cells in the global domain in I, J, K.
 * @param[in]  rand_logic_i_ptr Pointer to the RNG for i-dimension tasks [0,1).
 * @param[in]  rand_logic_j_ptr Pointer to the RNG for j-dimension tasks [0,1).
 * @param[in]  rand_logic_k_ptr Pointer to the RNG for k-dimension tasks [0,1).
 * @param[out] ci_metric_lnode_out Pointer to store the local i-node index of the selected cell's origin.
 * @param[out] cj_metric_lnode_out Pointer to store the local j-node index of the selected cell's origin.
 * @param[out] ck_metric_lnode_out Pointer to store the local k-node index of the selected cell's origin.
 * @param[out] xi_metric_logic_out Pointer to store the intra-cell logical xi-coordinate [0,1).
 * @param[out] eta_metric_logic_out Pointer to store the intra-cell logical eta-coordinate [0,1).
 * @param[out] zta_metric_logic_out Pointer to store the intra-cell logical zeta-coordinate [0,1).
 * @param[out] can_place_on_surface_out PETSC_TRUE if placement parameters were successfully determined, PETSC_FALSE otherwise.
 * @return PetscErrorCode 0 on success, or a PETSc error code.
 */
static PetscErrorCode DetermineSurfaceInitializationParameters(
    UserCtx *user, DMDALocalInfo *info,
    PetscInt xs_gnode, PetscInt ys_gnode, PetscInt zs_gnode,
    PetscInt IM_gcells, PetscInt JM_gcells, PetscInt KM_gcells,
    PetscRandom *rand_logic_i_ptr, PetscRandom *rand_logic_j_ptr, PetscRandom *rand_logic_k_ptr, /* Pointers to RNGs */
    PetscInt *ci_metric_lnode_out, PetscInt *cj_metric_lnode_out, PetscInt *ck_metric_lnode_out,
    PetscReal *xi_metric_logic_out, PetscReal *eta_metric_logic_out, PetscReal *zta_metric_logic_out,
    PetscBool *can_place_on_surface_out)
{
    PetscErrorCode ierr = 0;
    PetscReal r_val; // For storing random numbers from [0,1) RNGs
    PetscInt local_owned_cell_idx_face_dim1 = 0; // Local owned cell index in the first tangential dimension of the face
    PetscInt local_owned_cell_idx_face_dim2 = 0; // Local owned cell index in the second tangential dimension of the face
    PetscMPIInt rank_for_logging; // For more informative logs if needed

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank_for_logging); CHKERRQ(ierr); // Get rank for detailed logging

    *can_place_on_surface_out = PETSC_FALSE; // Default to: cannot place on surface from this rank

    // Check the flag set by ParseBCSFileForInlet
    if (user->inletFaceDefined == PETSC_FALSE) { // Or simply: if (!user->inletFaceDefined)
        LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP - Rank %d: Inlet face is not defined (from bcs.dat). Cannot place particle on surface.", rank_for_logging);
        PetscFunctionReturn(0); // Exit early if no inlet face was defined
    }
    
    // Default intra-cell logicals to cell center; will be overridden for the face-normal direction
    // and the other two will be randomized if placement occurs.
    *xi_metric_logic_out = 0.5; *eta_metric_logic_out = 0.5; *zta_metric_logic_out = 0.5;
    // Default cell node indices to the start of the rank's owned region
    *ci_metric_lnode_out = xs_gnode; *cj_metric_lnode_out = ys_gnode; *ck_metric_lnode_out = zs_gnode;

    // A rank must own some 3D cells to be able to define a face for particle placement.
    // info->mx,my,mz are number of NODES in local processor including ghosts if any that are not part of global domain.
    // So, number of CELLS is (info->mX - 1), (info->mY - 1), (info->mZ - 1) for non-zero node counts.
    PetscInt num_owned_cells_i = (info->mx > 1) ? info->mx - 1 : 0; // Num cells = Num nodes - 1
    PetscInt num_owned_cells_j = (info->my > 1) ? info->my - 1 : 0;
    PetscInt num_owned_cells_k = (info->mz > 1) ? info->mz - 1 : 0;

    if (num_owned_cells_i == 0 || num_owned_cells_j == 0 || num_owned_cells_k == 0) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP - Rank %d: Has zero cells in at least one dimension (owned cells i,j,k: %d,%d,%d). Cannot place particle on surface.",
                  rank_for_logging, num_owned_cells_i, num_owned_cells_j, num_owned_cells_k);
        PetscFunctionReturn(0);
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP - Rank %d: Processing inlet face %d. Owned cell counts (i,j,k): (%d,%d,%d). Global cell counts (I,J,K): (%d,%d,%d). Ghosted node starts (xs_g,ys_g,zs_g): (%d,%d,%d). Owned cell starts (info.xs/ys/ks): (%d,%d,%d)",
        rank_for_logging, user->identifiedInletBCFace, num_owned_cells_i, num_owned_cells_j, num_owned_cells_k, IM_gcells, JM_gcells, KM_gcells, xs_gnode, ys_gnode, zs_gnode, info->xs, info->ys, info->zs);

    switch (user->identifiedInletBCFace) {
        case BC_FACE_NEG_X: // Global I-MIN face (DMDA 'x'-min)
            // Check if this rank owns cells at the global I=0 boundary
            if (info->xs == 0) {
                *can_place_on_surface_out = PETSC_TRUE;
                *ci_metric_lnode_out = xs_gnode;         // Cell origin is at the first layer of owned nodes in i
                *xi_metric_logic_out = 1.0e-6;           // Intra-cell logical xi is near the "left" face of this cell

                // Tangential dimensions for an I-face are J and K. Select a cell randomly on this face.
                // Select local owned cell index in J-direction
                ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val); CHKERRQ(ierr); // Dereference RNG pointer
                local_owned_cell_idx_face_dim1 = (PetscInt)(r_val * num_owned_cells_j);
                local_owned_cell_idx_face_dim1 = PetscMin(PetscMax(0, local_owned_cell_idx_face_dim1), num_owned_cells_j - 1);
                *cj_metric_lnode_out = ys_gnode + local_owned_cell_idx_face_dim1;

                // Select local owned cell index in K-direction
                ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val); CHKERRQ(ierr); // Dereference RNG pointer
                local_owned_cell_idx_face_dim2 = (PetscInt)(r_val * num_owned_cells_k);
                local_owned_cell_idx_face_dim2 = PetscMin(PetscMax(0, local_owned_cell_idx_face_dim2), num_owned_cells_k - 1);
                *ck_metric_lnode_out = zs_gnode + local_owned_cell_idx_face_dim2;

                // Random intra-cell logical coordinates for eta and zeta (tangential to I-face)
                ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
                ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
                LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/NEG_Xi: Target. CellNode(i,j,k)=(%d,%d,%d). Logic(xi,et,zt)=(%.2e,%.2f,%.2f)", *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out, *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);
            } else { LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/NEG_Xi: Rank %d not on this face (info->xs=%d).", rank_for_logging, info->xs); }
            break;

        case BC_FACE_POS_X: // Global I-MAX face
            if (info->xs + num_owned_cells_i == IM_gcells) {
                *can_place_on_surface_out = PETSC_TRUE;
                *ci_metric_lnode_out = xs_gnode + (num_owned_cells_i - 1); // Node starting the last layer of owned cells in i
                *xi_metric_logic_out = 1.0 - 1.0e-6;            // Intra-cell logical xi is near the "right" face
                
                ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val); CHKERRQ(ierr);
                local_owned_cell_idx_face_dim1 = (PetscInt)(r_val * num_owned_cells_j);
                local_owned_cell_idx_face_dim1 = PetscMin(PetscMax(0, local_owned_cell_idx_face_dim1), num_owned_cells_j - 1);
                *cj_metric_lnode_out = ys_gnode + local_owned_cell_idx_face_dim1;

                ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val); CHKERRQ(ierr);
                local_owned_cell_idx_face_dim2 = (PetscInt)(r_val * num_owned_cells_k);
                local_owned_cell_idx_face_dim2 = PetscMin(PetscMax(0, local_owned_cell_idx_face_dim2), num_owned_cells_k - 1);
                *ck_metric_lnode_out = zs_gnode + local_owned_cell_idx_face_dim2;

                ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
                ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
                LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/POS_Xi: Target. CellNode(i,j,k)=(%d,%d,%d). Logic(xi,et,zt)=(%.2e,%.2f,%.2f)", *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out, *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);
            } else { LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/POS_Xi: Rank %d not on this face (info->xs=%d, num_i_cells=%d, IM_gcells=%d).", rank_for_logging, info->xs, num_owned_cells_i, IM_gcells); }
            break;

        case BC_FACE_NEG_Y: // Global J-MIN face
             if (info->ys == 0) {
                *can_place_on_surface_out = PETSC_TRUE;
                *cj_metric_lnode_out = ys_gnode;
                *eta_metric_logic_out = 1.0e-6;

                ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val); CHKERRQ(ierr);
                local_owned_cell_idx_face_dim1 = (PetscInt)(r_val * num_owned_cells_i);
                local_owned_cell_idx_face_dim1 = PetscMin(PetscMax(0, local_owned_cell_idx_face_dim1), num_owned_cells_i - 1);
                *ci_metric_lnode_out = xs_gnode + local_owned_cell_idx_face_dim1;
                
                ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val); CHKERRQ(ierr);
                local_owned_cell_idx_face_dim2 = (PetscInt)(r_val * num_owned_cells_k);
                local_owned_cell_idx_face_dim2 = PetscMin(PetscMax(0, local_owned_cell_idx_face_dim2), num_owned_cells_k - 1);
                *ck_metric_lnode_out = zs_gnode + local_owned_cell_idx_face_dim2;
                
                ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr);
                ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
                LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/NEG_Eta: Target. CellNode(i,j,k)=(%d,%d,%d). Logic(xi,et,zt)=(%.2f,%.2e,%.2f)", *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out, *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);
            } else { LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/NEG_Eta: Rank %d not on this face (info->ys=%d).", rank_for_logging, info->ys); }
            break;

        case BC_FACE_POS_Y: // Global J-MAX face
            if (info->ys + num_owned_cells_j == JM_gcells) {
                *can_place_on_surface_out = PETSC_TRUE;
                *cj_metric_lnode_out = ys_gnode + (num_owned_cells_j - 1);
                *eta_metric_logic_out = 1.0 - 1.0e-6;

                ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val); CHKERRQ(ierr);
                local_owned_cell_idx_face_dim1 = (PetscInt)(r_val * num_owned_cells_i);
                local_owned_cell_idx_face_dim1 = PetscMin(PetscMax(0, local_owned_cell_idx_face_dim1), num_owned_cells_i - 1);
                *ci_metric_lnode_out = xs_gnode + local_owned_cell_idx_face_dim1;

                ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val); CHKERRQ(ierr);
                local_owned_cell_idx_face_dim2 = (PetscInt)(r_val * num_owned_cells_k);
                local_owned_cell_idx_face_dim2 = PetscMin(PetscMax(0, local_owned_cell_idx_face_dim2), num_owned_cells_k - 1);
                *ck_metric_lnode_out = zs_gnode + local_owned_cell_idx_face_dim2;

                ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr);
                ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
                LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/POS_Eta: Target. CellNode(i,j,k)=(%d,%d,%d). Logic(xi,et,zt)=(%.2f,%.2e,%.2f)", *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out, *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);
            } else { LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/POS_Eta: Rank %d not on this face (info->ys=%d, num_j_cells=%d, JM_gcells=%d).", rank_for_logging, info->ys, num_owned_cells_j, JM_gcells); }
            break;

        case BC_FACE_NEG_Z: // Global K-MIN face
            if (info->zs == 0) {
                *can_place_on_surface_out = PETSC_TRUE;
                *ck_metric_lnode_out = zs_gnode;
                *zta_metric_logic_out = 1.0e-6;

                ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val); CHKERRQ(ierr);
                local_owned_cell_idx_face_dim1 = (PetscInt)(r_val * num_owned_cells_i); // Scale by num_owned_cells_i
                *ci_metric_lnode_out = xs_gnode + PetscMin(PetscMax(0, local_owned_cell_idx_face_dim1), num_owned_cells_i - 1); // Clamp with num_owned_cells_i -1

                ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val); CHKERRQ(ierr);
                local_owned_cell_idx_face_dim2 = (PetscInt)(r_val * num_owned_cells_j); // Scale by num_owned_cells_j
                *cj_metric_lnode_out = ys_gnode + PetscMin(PetscMax(0, local_owned_cell_idx_face_dim2), num_owned_cells_j - 1); // Clamp with num_owned_cells_j - 1
                
                LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/NEG_Zeta Pre-Metric: r_val_i=%.10g, num_i_cells=%d, unclamp_idx_i=%d, clamp_idx_i=%d (rel. to owned), xs_gnode=%d, ci_node=%d",
                          (double)r_val, num_owned_cells_i, (PetscInt)(r_val*num_owned_cells_i), local_owned_cell_idx_face_dim1, xs_gnode, *ci_metric_lnode_out); // r_val here is last one used (for j)
                LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/NEG_Zeta Pre-Metric: r_val_j=%.10g, num_j_cells=%d, unclamp_idx_j=%d, clamp_idx_j=%d (rel. to owned), ys_gnode=%d, cj_node=%d",
                          (double)r_val, num_owned_cells_j, (PetscInt)(r_val*num_owned_cells_j), local_owned_cell_idx_face_dim2, ys_gnode, *cj_metric_lnode_out); // r_val here is last one used (for j)


                ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr);
                ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
                LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/NEG_Zeta: Target. CellNode(i,j,k)=(%d,%d,%d). Logic(xi,et,zt)=(%.2f,%.2f,%.2e)", *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out, *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);
            } else { LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/NEG_Zeta: Rank %d not on this face (info->zs=%d).", rank_for_logging, info->zs); }
            break;

        case BC_FACE_POS_Z: // Global K-MAX face
            if (info->zs + num_owned_cells_k == KM_gcells) {
                *can_place_on_surface_out = PETSC_TRUE;
                *ck_metric_lnode_out = zs_gnode + (num_owned_cells_k - 1);
                *zta_metric_logic_out = 1.0 - 1.0e-6;

                ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val); CHKERRQ(ierr);
                local_owned_cell_idx_face_dim1 = (PetscInt)(r_val * num_owned_cells_i);
                *ci_metric_lnode_out = xs_gnode + PetscMin(PetscMax(0, local_owned_cell_idx_face_dim1), num_owned_cells_i - 1);

                ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val); CHKERRQ(ierr);
                local_owned_cell_idx_face_dim2 = (PetscInt)(r_val * num_owned_cells_j);
                *cj_metric_lnode_out = ys_gnode + PetscMin(PetscMax(0, local_owned_cell_idx_face_dim2), num_owned_cells_j - 1);

                ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr);
                ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
                LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/POS_Zeta: Target. CellNode(i,j,k)=(%d,%d,%d). Logic(xi,et,zt)=(%.2f,%.2f,%.2e)", *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out, *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);
            } else { LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/POS_Zeta: Rank %d not on this face (info->zs=%d, num_k_cells=%d, KM_gcells=%d).", rank_for_logging, info->zs, num_owned_cells_k, KM_gcells); }
            break;
        default:
            LOG_ALLOW(LOCAL, LOG_ERROR, "DSP - Rank %d: Invalid user->identifiedInletBCFace value: %d\n", rank_for_logging, user->identifiedInletBCFace);
            *can_place_on_surface_out = PETSC_FALSE; // Should already be false
            break;
    }

    if (*can_place_on_surface_out) {
        // Clamp the two random intra-cell logical coordinates (tangential to the face)
        // The one normal to the face is already set (e.g., xi_metric_logic_out for I-faces)
        if (user->identifiedInletBCFace == BC_FACE_NEG_X || user->identifiedInletBCFace == BC_FACE_POS_X) {
            *eta_metric_logic_out = PetscMin(*eta_metric_logic_out, 1.0 - 1.0e-7);
            *zta_metric_logic_out = PetscMin(*zta_metric_logic_out, 1.0 - 1.0e-7);
        } else if (user->identifiedInletBCFace == BC_FACE_NEG_Y || user->identifiedInletBCFace == BC_FACE_POS_Y) {
            *xi_metric_logic_out  = PetscMin(*xi_metric_logic_out,  1.0 - 1.0e-7);
            *zta_metric_logic_out = PetscMin(*zta_metric_logic_out, 1.0 - 1.0e-7);
        } else { // Z-faces (NEG_Z or POS_Z)
            *xi_metric_logic_out  = PetscMin(*xi_metric_logic_out,  1.0 - 1.0e-7);
            *eta_metric_logic_out = PetscMin(*eta_metric_logic_out, 1.0 - 1.0e-7);
        }
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Determines cell selection and intra-cell logical coordinates for volumetric initialization (Mode 1).
 *
 * This function is called when `user->ParticleInitialization == 1`. It randomly selects
 * an *owned cell* on the current MPI rank and then generates random intra-cell logical
 * coordinates `[0,1)^3` within that chosen cell.
 *
 * The process involves:
 * 1. Checking if the current MPI rank owns any 3D cells.
 * 2. If it does, it randomly selects an *owned cell index* in each logical direction (i, j, k)
 *    by scaling a `[0,1)` random number with the number of owned cells in that direction.
 * 3. These local owned cell indices are then converted to the *local node indices*
 *    (`ci/cj/ck_metric_lnode_out`) corresponding to the origin of the selected cell,
 *    for use with `MetricLogicalToPhysical`. This conversion uses `xs/ys/zs_gnode`.
 * 4. All three intra-cell logical coordinates (`xi/eta/zta_metric_logic_out`) for
 *    `MetricLogicalToPhysical` are chosen randomly within the `[0,1)` range.
 * 5. A flag (`can_place_in_volume_out`) indicates if a valid placement could be determined.
 *
 * **Important Note on `DMDALocalInfo info` members (same as for surface init):**
 *   - `info->xs, info->ys, info->ks`: Global starting indices of *owned cells*.
 *   - `info->mx, info->my, info->mz`: Number of *grid points (nodes)* in each local dimension on this process.
 *     Therefore, the number of *owned cells* in a dimension is `info->mX - 1` (if `info->mX > 0`).
 *
 * @param[in]  user Pointer to `UserCtx`. (Currently not used in this specific helper, but kept for API consistency).
 * @param[in]  info Pointer to `DMDALocalInfo` for the current rank's grid portion.
 * @param[in]  xs_gnode, ys_gnode, zs_gnode Local indices (in the ghosted array) of the first *owned node*.
 * @param[in]  rand_logic_i_ptr Pointer to the RNG for i-dimension tasks [0,1).
 * @param[in]  rand_logic_j_ptr Pointer to the RNG for j-dimension tasks [0,1).
 * @param[in]  rand_logic_k_ptr Pointer to the RNG for k-dimension tasks [0,1).
 * @param[out] ci_metric_lnode_out Pointer to store the local i-node index of the selected cell's origin.
 * @param[out] cj_metric_lnode_out Pointer to store the local j-node index of the selected cell's origin.
 * @param[out] ck_metric_lnode_out Pointer to store the local k-node index of the selected cell's origin.
 * @param[out] xi_metric_logic_out Pointer to store the intra-cell logical xi-coordinate [0,1).
 * @param[out] eta_metric_logic_out Pointer to store the intra-cell logical eta-coordinate [0,1).
 * @param[out] zta_metric_logic_out Pointer to store the intra-cell logical zeta-coordinate [0,1).
 * @param[out] can_place_in_volume_out PETSC_TRUE if placement parameters were successfully determined, PETSC_FALSE otherwise.
 * @return PetscErrorCode 0 on success, or a PETSc error code.
 */
static PetscErrorCode DetermineVolumetricInitializationParameters(
    UserCtx *user, DMDALocalInfo *info,
    PetscInt xs_gnode, PetscInt ys_gnode, PetscInt zs_gnode,
    PetscRandom *rand_logic_i_ptr, PetscRandom *rand_logic_j_ptr, PetscRandom *rand_logic_k_ptr, /* Pointers to RNGs */
    PetscInt *ci_metric_lnode_out, PetscInt *cj_metric_lnode_out, PetscInt *ck_metric_lnode_out,
    PetscReal *xi_metric_logic_out, PetscReal *eta_metric_logic_out, PetscReal *zta_metric_logic_out,
    PetscBool *can_place_in_volume_out)
{
    PetscErrorCode ierr = 0;
    PetscReal      r_val; // Temporary for random numbers from [0,1) RNGs
    PetscInt       local_owned_cell_idx_i, local_owned_cell_idx_j, local_owned_cell_idx_k;
    PetscMPIInt    rank_for_logging; // For logging if needed

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank_for_logging); CHKERRQ(ierr);

    *can_place_in_volume_out = PETSC_FALSE; // Default to: cannot place

    // Default intra-cell logicals and cell node indices (e.g. if placement fails)
    *xi_metric_logic_out = 0.5; *eta_metric_logic_out = 0.5; *zta_metric_logic_out = 0.5;
    *ci_metric_lnode_out = xs_gnode; *cj_metric_lnode_out = ys_gnode; *ck_metric_lnode_out = zs_gnode;

    // Calculate number of owned cells in each direction from node counts in info
    // Assumes info->mx, info->my, info->mz are node counts on this process for each dimension.
    // Number of cells = Number of nodes - 1 (if > 0 nodes).
    PetscInt num_owned_cells_i = (info->mx > 1) ? info->mx - 1 : 0;
    PetscInt num_owned_cells_j = (info->my > 1) ? info->my - 1 : 0;
    PetscInt num_owned_cells_k = (info->mz > 1) ? info->mz - 1 : 0;

    if (num_owned_cells_i > 0 && num_owned_cells_j > 0 && num_owned_cells_k > 0) { // If rank owns any 3D cells
        *can_place_in_volume_out = PETSC_TRUE;

        // --- 1. Select a Random Owned Cell ---
        // The selected index will be a 0-based index relative to the start of this rank's owned cells.

        // Select random local owned cell index in I-direction
        ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val); CHKERRQ(ierr); // Dereference RNG pointer
        local_owned_cell_idx_i = (PetscInt)(r_val * num_owned_cells_i);
        // Clamp to be safe: local_owned_cell_idx_i should be in [0, num_owned_cells_i - 1]
        local_owned_cell_idx_i = PetscMin(PetscMax(0, local_owned_cell_idx_i), num_owned_cells_i - 1);
        *ci_metric_lnode_out = xs_gnode + local_owned_cell_idx_i; // Convert to local node index for cell origin

        // Select random local owned cell index in J-direction
        ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val); CHKERRQ(ierr); // Dereference RNG pointer
        local_owned_cell_idx_j = (PetscInt)(r_val * num_owned_cells_j);
        local_owned_cell_idx_j = PetscMin(PetscMax(0, local_owned_cell_idx_j), num_owned_cells_j - 1);
        *cj_metric_lnode_out = ys_gnode + local_owned_cell_idx_j;

        // Select random local owned cell index in K-direction
        ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val); CHKERRQ(ierr); // Dereference RNG pointer
        local_owned_cell_idx_k = (PetscInt)(r_val * num_owned_cells_k);
        local_owned_cell_idx_k = PetscMin(PetscMax(0, local_owned_cell_idx_k), num_owned_cells_k - 1);
        *ck_metric_lnode_out = zs_gnode + local_owned_cell_idx_k;

        LOG_ALLOW(LOCAL, LOG_DEBUG, "DVP - Rank %d: Selected Cell (Owned Idx: %d,%d,%d -> LNodeStart: %d,%d,%d). OwnedCells(i,j,k): (%d,%d,%d). GhostNodeStarts(xs,ys,zs): (%d,%d,%d)",
                  rank_for_logging, local_owned_cell_idx_i, local_owned_cell_idx_j, local_owned_cell_idx_k,
                  *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out,
                  num_owned_cells_i, num_owned_cells_j, num_owned_cells_k,
                  xs_gnode, ys_gnode, zs_gnode);


        // --- 2. Generate Random Intra-Cell Logical Coordinates [0,1) for MetricLogicalToPhysical ---
        ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out);  CHKERRQ(ierr); // Re-use RNGs
        ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
        ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);

        // Ensure logical coordinates are strictly within [0,1) for robustness with MetricLogicalToPhysical
        *xi_metric_logic_out  = PetscMin(*xi_metric_logic_out,  1.0 - 1.0e-7);
        *eta_metric_logic_out = PetscMin(*eta_metric_logic_out, 1.0 - 1.0e-7);
        *zta_metric_logic_out = PetscMin(*zta_metric_logic_out, 1.0 - 1.0e-7);
        // Ensure they are not negative either (though [0,1) RNGs shouldn't produce this)
        *xi_metric_logic_out  = PetscMax(*xi_metric_logic_out,  0.0);
        *eta_metric_logic_out = PetscMax(*eta_metric_logic_out, 0.0);
        *zta_metric_logic_out = PetscMax(*zta_metric_logic_out, 0.0);

    } else {
        // This rank does not own any 3D cells (e.g., in a 1D or 2D decomposition,
        // or if the global domain itself is not 3D in terms of cells).
        // *can_place_in_volume_out remains PETSC_FALSE.
        LOG_ALLOW(LOCAL, LOG_WARNING, "DVP - Rank %d: Cannot place particle volumetrically. Rank has zero owned cells in at least one dimension (owned cells i,j,k: %d,%d,%d).\n",
                  rank_for_logging, num_owned_cells_i, num_owned_cells_j, num_owned_cells_k);
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Initializes basic particle properties: physical position, particle ID, and cell ID placeholder.
 *
 * This function orchestrates the initialization of particles based on `user->ParticleInitialization`.
 * It retrieves necessary grid and swarm information, then loops through particles, calling
 * helper subroutines to determine cell selection and intra-cell logical coordinates.
 * It then performs the logical-to-physical transformation and sets standard particle fields.
 *
 * @param[in,out] user               Pointer to the `UserCtx` structure.
 * @param[in]     particlesPerProcess Number of particles to initialize on this MPI process.
 * @param[in]     rand_logic_i       RNG for i-dimension tasks, range `[0.0, 1.0)`.
 * @param[in]     rand_logic_j       RNG for j-dimension tasks, range `[0.0, 1.0)`.
 * @param[in]     rand_logic_k       RNG for k-dimension tasks, range `[0.0, 1.0)`.
 * @param[in]     bboxlist           (Currently unused).
 * @return PetscErrorCode Returns `0` on success, or a PETSc error code on failure.
 */
static PetscErrorCode InitializeParticleBasicProperties(UserCtx *user,
                                                   PetscInt particlesPerProcess,
                                                   PetscRandom *rand_logic_i,
                                                   PetscRandom *rand_logic_j,
                                                   PetscRandom *rand_logic_k,
                                                   BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    DM             swarm;
    PetscReal      *positions_field = NULL, *pos_phy_field = NULL;
    PetscInt64     *particleIDs = NULL, *cellIDs_petsc = NULL;
    PetscMPIInt    rank;
    const Cmpnts   ***coor_nodes_local_array;
    Vec            Coor_local;
    DMDALocalInfo  info;
    PetscInt       xs_gnode, ys_gnode, zs_gnode;
    PetscInt       xm_onode, ym_onode, zm_onode; // Unused in this main func but fetched
    PetscInt       IM_gcells, JM_gcells, KM_gcells;

    PetscFunctionBeginUser;

    // --- 1. Input Validation and Basic Setup ---
    if (!user || !rand_logic_i || !rand_logic_j || !rand_logic_k) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "InitializeParticleBasicProperties - Null user or RNG pointer.");
    }
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    swarm = user->swarm;

    ierr = DMGetCoordinatesLocal(user->da, &Coor_local); CHKERRQ(ierr);
    if (!Coor_local) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "DMGetCoordinatesLocal for user->da returned NULL.");
    ierr = DMDAVecGetArrayRead(user->fda, Coor_local, (void*)&coor_nodes_local_array); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_INFO, "InitializeParticleBasicProperties - Rank %d: Initializing %d particles. Mode: %d. Identified Inlet (if Mode 0): %d\n",
              rank, particlesPerProcess, user->ParticleInitialization, (user->ParticleInitialization == 0) ? user->identifiedInletBCFace : -1);

    ierr = DMSwarmGetField(swarm, "position",       NULL, NULL, (void**)&positions_field); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "pos_phy",        NULL, NULL, (void**)&pos_phy_field);   CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid",    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs_petsc);  CHKERRQ(ierr);

    // --- 2. Get Grid Information (once for all particles on this rank) ---
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = DMDAGetGhostCorners(user->da, &xs_gnode, &ys_gnode, &zs_gnode, &xm_onode, &ym_onode, &zm_onode); CHKERRQ(ierr);
    ierr = DMDAGetInfo(user->da, NULL, &IM_gcells, &JM_gcells, &KM_gcells, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL); CHKERRQ(ierr);

    // --- 3. Loop Over Particles ---
    for (PetscInt p = 0; p < particlesPerProcess; p++) {
        PetscInt  ci_metric_lnode = 0, cj_metric_lnode = 0, ck_metric_lnode = 0; // Cell origin local node index
        PetscReal xi_metric_logic = 0.5, eta_metric_logic = 0.5, zta_metric_logic = 0.5; // Intra-cell logical [0,1]
        Cmpnts    phys_coords = {0.0, 0.0, 0.0}; // Default physical coordinates
        PetscBool can_place_particle = PETSC_FALSE;    // Flag: can this particle be properly placed by this rank

        // --- 3.a. Determine Cell and Intra-Cell Logical Coordinates ---
        if (user->ParticleInitialization == 0) {
            ierr = DetermineSurfaceInitializationParameters(user, &info, xs_gnode, ys_gnode, zs_gnode,
                                                            IM_gcells, JM_gcells, KM_gcells,
                                                            rand_logic_i, rand_logic_j, rand_logic_k,
                                                            &ci_metric_lnode, &cj_metric_lnode, &ck_metric_lnode,
                                                            &xi_metric_logic, &eta_metric_logic, &zta_metric_logic,
                                                            &can_place_particle); CHKERRQ(ierr);
        } else { // user->ParticleInitialization == 1
            ierr = DetermineVolumetricInitializationParameters(user, &info, xs_gnode, ys_gnode, zs_gnode,
                                                               rand_logic_i, rand_logic_j, rand_logic_k,
                                                               &ci_metric_lnode, &cj_metric_lnode, &ck_metric_lnode,
                                                               &xi_metric_logic, &eta_metric_logic, &zta_metric_logic,
                                                               &can_place_particle); CHKERRQ(ierr);
        }

        // --- 3.b. Perform Logical to Physical Transformation if placement is possible ---
        if (can_place_particle) {
            ierr = MetricLogicalToPhysical(user, coor_nodes_local_array,
                                           ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                                           xi_metric_logic, eta_metric_logic, zta_metric_logic,
                                           &phys_coords); CHKERRQ(ierr);
        } else {
            // phys_coords remains (0,0,0) or whatever default you prefer for unplaced particles
            LOG_LOOP_ALLOW(LOCAL, LOG_WARNING, p, 1,
                "InitProps - Rank %d: Particle %d - Placement parameters not determined (e.g., not on inlet face, or no owned cells). Phys default: (%.2f,%.2f,%.2f).\n",
                rank, p, phys_coords.x, phys_coords.y, phys_coords.z);
        }

        // --- 3.c. Store Common Particle Properties ---
        positions_field[3*p+0] = phys_coords.x; positions_field[3*p+1] = phys_coords.y; positions_field[3*p+2] = phys_coords.z;
        pos_phy_field[3*p+0]   = phys_coords.x; pos_phy_field[3*p+1]   = phys_coords.y; pos_phy_field[3*p+2]   = phys_coords.z;
        particleIDs[p]         = (PetscInt64)rank * particlesPerProcess + p;
        cellIDs_petsc[3*p+0]   = -1; cellIDs_petsc[3*p+1] = -1; cellIDs_petsc[3*p+2] = -1;

        // --- 3.d. Logging for this particle ---
        if (can_place_particle) {
            LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, p, (particlesPerProcess > 20 ? particlesPerProcess/10 : 1),
                "InitProps - Rank %d: PID %lld (idx %d) PLACED. Mode %d. Cell(LNodeStart):(%d,%d,%d). Logic(Metric): (%.2e,%.2f,%.2f). Phys: (%.6f,%.6f,%.6f).\n",
                rank, (long long)particleIDs[p], p, user->ParticleInitialization, ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                xi_metric_logic, eta_metric_logic, zta_metric_logic, phys_coords.x, phys_coords.y, phys_coords.z);
        } // Warning for unplaced particles already logged in step 3.b
    } // --- End of per-particle loop ---

    // --- 4. Restore Pointers ---
    ierr = DMSwarmRestoreField(swarm, "position",       NULL, NULL, (void**)&positions_field); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "pos_phy",        NULL, NULL, (void**)&pos_phy_field);   CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid",    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs_petsc);  CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, Coor_local, (void*)&coor_nodes_local_array); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_INFO, "InitializeParticleBasicProperties - Rank %d: Completed processing for %d particles.\n",
              rank, particlesPerProcess);
    PetscFunctionReturn(0);
}

/**
 * @brief Helper function to update a given particle’s field value.
 *
 * This function performs conditional, point-level updates for a swarm field based on its name.
 * For example, you might want to initialize the "velocity" field to 0.0, but the "temperature"
 * field to a nonzero default (e.g., 300.0). This function can be extended for other fields.
 *
 * @param[in] fieldName  Name of the swarm field.
 * @param[in] p          Particle index.
 * @param[in] fieldDim   Dimension of the field.
 * @param[out] fieldData Pointer to the field’s data array.
 *
 * @return PetscErrorCode Returns 0 on success.
 */
static PetscErrorCode UpdateSwarmFieldValue(const char *fieldName, PetscInt p, PetscInt fieldDim, PetscReal *fieldData)
{
  PetscFunctionBeginUser;
  if (strcmp(fieldName, "velocity") == 0) {
    // For velocity, initialize all components to zero
    for (PetscInt d = 0; d < fieldDim; d++) {
      fieldData[fieldDim * p + d] = 0.0;
    }
  } else if (strcmp(fieldName, "temperature") == 0) {
    // For temperature, for example, initialize to a default value (e.g., 300.0)
    for (PetscInt d = 0; d < fieldDim; d++) {
      fieldData[fieldDim * p + d] = 300.0;
    }
  } else if (strcmp(fieldName, "P") == 0) {
    // For pressure, initialize to a default value (e.g., 101325.0)
    for (PetscInt d = 0; d < fieldDim; d++) {
      fieldData[fieldDim * p + d] = 101325.0;
    }
  } else {
    // Default: initialize all components to zero
    for (PetscInt d = 0; d < fieldDim; d++) {
      fieldData[fieldDim * p + d] = 0.0;
    }
  }
  PetscFunctionReturn(0);
}

/**
 * @brief Initializes a generic swarm field with point-level updates.
 *
 * This field-agnostic function retrieves the specified swarm field (which may be
 * scalar or multi-component) and initializes each particle's entry using a helper
 * that performs conditional updates based on the field name.
 *
 * @param[in,out] user       Pointer to the UserCtx structure containing the swarm.
 * @param[in]     fieldName  Name of the swarm field to initialize.
 * @param[in]     fieldDim   Dimension of the field (e.g., 1 for scalar, 3 for vector).
 *
 * @return PetscErrorCode    Returns 0 on success, non-zero on failure.
 */
static PetscErrorCode AssignInitialFieldToSwarm(UserCtx *user, const char *fieldName, PetscInt fieldDim)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscReal     *fieldData = NULL;
    PetscInt       nLocal;

    PetscFunctionBeginUser;
    
    // Get the number of local particles
    ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO, "AssignInitialFieldToSwarm - %d local particles found.\n", nLocal);

    // Retrieve the swarm field pointer for the specified fieldName
    ierr = DMSwarmGetField(swarm, fieldName, NULL, NULL, (void**)&fieldData); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG, "AssignInitialFieldToSwarm - Retrieved field '%s'.\n", fieldName);

    // Loop over all particles and update the field using the helper function
    for (PetscInt p = 0; p < nLocal; p++) {
        ierr = UpdateSwarmFieldValue(fieldName, p, fieldDim, fieldData); CHKERRQ(ierr);
        LOG_LOOP_ALLOW(LOCAL,LOG_DEBUG,p, 100,
            "AssignInitialFieldToSwarm - Particle %d: %s = [", p, fieldName);
        for (PetscInt d = 0; d < fieldDim; d++) {
	  LOG_ALLOW(LOCAL,LOG_DEBUG,"%.6f ", fieldData[fieldDim * p + d]);
        }
        LOG_ALLOW(LOCAL,LOG_DEBUG,"]\n");
    }
    
    // Restore the swarm field pointer
    ierr = DMSwarmRestoreField(swarm, fieldName, NULL, NULL, (void**)&fieldData); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO, "AssignInitialFieldToSwarm - Initialization of field '%s' complete.\n", fieldName);
    
    PetscFunctionReturn(0);
}

/**
 * @brief Initializes all particle properties in the swarm.
 *
 * This function orchestrates the initialization of particle properties.
 * It first determines the inlet face if surface initialization (Mode 0) is selected
 * by parsing "bcs.dat".
 * Then, it initializes basic particle properties (physical position, Particle ID,
 * and placeholder Cell IDs) by calling `InitializeParticleBasicProperties`. This call
 * uses the provided `rand_logic_i/j/k` RNGs, which must be pre-initialized for [0,1).
 * The `rand_phys_x/y/z` RNGs (physically bounded) are passed but may not be used by
 * `InitializeParticleBasicProperties` for position setting if all initialization paths
 * use logical-to-physical mapping.
 * Finally, it calls helper functions to initialize other registered swarm fields
 * like "velocity", "weight", and "P" (pressure) to default values.
 *
 * @param[in,out] user               Pointer to the `UserCtx` structure.
 * @param[in]     particlesPerProcess Number of particles assigned to this MPI process.
 * @param[in]     rand_phys_x        RNG for physical x-coordinates (from `InitializeRandomGenerators`).
 * @param[in]     rand_phys_y        RNG for physical y-coordinates (from `InitializeRandomGenerators`).
 * @param[in]     rand_phys_z        RNG for physical z-coordinates (from `InitializeRandomGenerators`).
 * @param[in]     rand_logic_i       RNG for i-logical dimension tasks [0,1) (from `InitializeLogicalSpaceRNGs`).
 * @param[in]     rand_logic_j       RNG for j-logical dimension tasks [0,1) (from `InitializeLogicalSpaceRNGs`).
 * @param[in]     rand_logic_k       RNG for k-logical dimension tasks [0,1) (from `InitializeLogicalSpaceRNGs`).
 * @param[in]     bboxlist           Array of BoundingBox structures (potentially unused by IPBP).
 *
 * @return PetscErrorCode            Returns 0 on success, non-zero on failure.
 */
PetscErrorCode AssignInitialPropertiesToSwarm(UserCtx* user,
                                              PetscInt particlesPerProcess,
                                              PetscRandom *rand_phys_x, // RNG from original InitializeRandomGenerators
                                              PetscRandom *rand_phys_y, // RNG from original InitializeRandomGenerators
                                              PetscRandom *rand_phys_z, // RNG from original InitializeRandomGenerators
                                              PetscRandom *rand_logic_i, // RNG from InitializeLogicalSpaceRNGs
                                              PetscRandom *rand_logic_j, // RNG from InitializeLogicalSpaceRNGs
                                              PetscRandom *rand_logic_k, // RNG from InitializeLogicalSpaceRNGs
                                              BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    // --- 0. Input Validation ---
    if (!user || !bboxlist || !rand_logic_i || !rand_logic_j || !rand_logic_k || !rand_phys_x || !rand_phys_y || !rand_phys_z) {
        // Check all RNGs now as they are passed in
        LOG_ALLOW(GLOBAL, LOG_ERROR, "AssignInitialPropertiesToSwarm - Null user, bboxlist, or RNG pointer.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "AssignInitialPropertiesToSwarm - Null input detected.");
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "AssignInitialPropertiesToSwarm - Initializing swarm with %d particles per process. Mode: %d.\n",
              particlesPerProcess, user->ParticleInitialization);

    // --- 1. Parse BCS File for Inlet Information (if Mode 0) ---
    if (user->ParticleInitialization == 0) {
      if(user->inletFaceDefined == PETSC_FALSE){
	LOG_ALLOW(GLOBAL, LOG_ERROR, "AssignInitialPropertiesToSwarm - ParticleInitialization Mode 0 (Surface Init) selected, but no INLET face was identified from bcs.dat. Cannot proceed.\n");
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "ParticleInitialization Mode 0 requires an INLET face to be defined in bcs.dat.");
      }else{
	LOG_ALLOW(GLOBAL, LOG_INFO, "AssignInitialPropertiesToSwarm - After ParseBCSFileForInlet, identifiedInletBCFace = %d\n", user->identifiedInletBCFace);
      }
    }

    // --- 2. Initialize Basic Particle Properties (Position, PID, Cell IDs placeholder) ---
    // The rand_logic_i/j/k are now passed directly.
    // The rand_phys_x/y/z are passed but InitializeParticleBasicProperties (refactored version)
    // will not use them for setting positions if all its paths use logical-to-physical mapping.
    LOG_ALLOW(LOCAL, LOG_DEBUG, "AssignInitialPropertiesToSwarm - Calling InitializeParticleBasicProperties.\n");
    ierr = InitializeParticleBasicProperties(user, particlesPerProcess,
                                             rand_logic_i, rand_logic_j, rand_logic_k,
                                             bboxlist); // bboxlist passed along
    CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "AssignInitialPropertiesToSwarm - Successfully initialized basic particle properties.\n");

    // Note: The logical RNGs (rand_logic_i/j/k) are NOT destroyed here.
    // They were created externally (e.g., by InitializeLogicalSpaceRNGs) and
    // should be destroyed externally (e.g., in FinalizeSwarmSetup).
    // Same for rand_phys_x/y/z.

    // --- 3. Initialize Other Swarm Fields (Velocity, Weight, Pressure, etc.) ---
    LOG_ALLOW(LOCAL, LOG_DEBUG, "AssignInitialPropertiesToSwarm - Initializing 'velocity' field.\n");
    ierr = AssignInitialFieldToSwarm(user, "velocity", 3); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "AssignInitialPropertiesToSwarm - 'velocity' field initialization complete.\n");

    LOG_ALLOW(LOCAL, LOG_DEBUG, "AssignInitialPropertiesToSwarm - Initializing 'weight' field.\n");
    ierr = AssignInitialFieldToSwarm(user, "weight", 3); CHKERRQ(ierr); // Assuming weight is vec3
    LOG_ALLOW(LOCAL, LOG_INFO, "AssignInitialPropertiesToSwarm - 'weight' field initialization complete.\n");

    LOG_ALLOW(LOCAL, LOG_DEBUG, "AssignInitialPropertiesToSwarm - Initializing 'P' (Pressure) field.\n");
    ierr = AssignInitialFieldToSwarm(user, "P", 1); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "AssignInitialPropertiesToSwarm - 'P' field initialization complete.\n");
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "AssignInitialPropertiesToSwarm - Successfully completed all swarm property initialization.\n");

    PetscFunctionReturn(0);
}

/**
 * @brief Distributes particles evenly across MPI processes, handling any remainders.
 *
 * This function calculates the number of particles each MPI process should handle,
 * distributing the remainder particles to the first few ranks if necessary.
 *
 * @param[in]     numParticles       Total number of particles to create across all MPI processes.
 * @param[in]     rank               MPI rank of the current process.
 * @param[in]     size               Total number of MPI processes.
 * @param[out]    particlesPerProcess Number of particles assigned to the current MPI process.
 * @param[out]    remainder           Remainder particles when dividing numParticles by size.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode DistributeParticles(PetscInt numParticles, PetscMPIInt rank, PetscMPIInt size, PetscInt* particlesPerProcess, PetscInt* remainder) {

    // Calculate the base number of particles per process
    *particlesPerProcess = numParticles / size;
    *remainder = numParticles % size;

    // Distribute the remainder particles to the first 'remainder' ranks
    if (rank < *remainder) {
        *particlesPerProcess += 1;
        LOG_ALLOW_SYNC(GLOBAL,LOG_INFO,"DistributeParticles - Rank %d receives an extra particle. Total: %d\n", rank, *particlesPerProcess);
    } else {
        LOG_ALLOW_SYNC(GLOBAL,LOG_INFO, "DistributeParticles - Rank %d receives %d particles.\n", rank, *particlesPerProcess);
    }

    return 0;
}

/**
 * @brief Finalizes the swarm setup by destroying random generators and logging completion.
 *
 * This function cleans up resources by destroying random number generators and LOG_ALLOWs the completion of swarm setup.
 *
 * @param[in]     randx             Random number generator for the x-coordinate.
 * @param[in]     randy             Random number generator for the y-coordinate.
 * @param[in]     randz             Random number generator for the z-coordinate.
 * @param[in]     rand_logic_i      Random number generator for the xi-coordinate.
 * @param[in]     rand_logic_j      Random number generator for the eta-coordinate.
 * @param[in]     rand_logic_k      Random number generator for the zeta-coordinate.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode FinalizeSwarmSetup(PetscRandom *randx, PetscRandom *randy, PetscRandom *randz, PetscRandom *rand_logic_i, PetscRandom *rand_logic_j, PetscRandom *rand_logic_k) {
    PetscErrorCode ierr;  // Error code for PETSc functions
    PetscInt  ParticleInitialization; 

    ierr = PetscOptionsGetInt(NULL, NULL, "-pinit", &ParticleInitialization, NULL); CHKERRQ(ierr);
 
    if(ParticleInitialization==1){

      // Destroy random number generators to free resources
      // Physical space
      ierr = PetscRandomDestroy(randx); CHKERRQ(ierr);
      ierr = PetscRandomDestroy(randy); CHKERRQ(ierr);
      ierr = PetscRandomDestroy(randz); CHKERRQ(ierr);

      // Logical space
      ierr = PetscRandomDestroy(rand_logic_i); CHKERRQ(ierr);
      ierr = PetscRandomDestroy(rand_logic_j); CHKERRQ(ierr);
      ierr = PetscRandomDestroy(rand_logic_k); CHKERRQ(ierr);      
      
      LOG_ALLOW(LOCAL,LOG_DEBUG,"FinalizeSwarmSetup - Destroyed all random number generators.\n");
    }else if(ParticleInitialization==0){
      LOG_ALLOW(LOCAL,LOG_DEBUG,"FinalizeSwarmSetup - Not a Random Initialization of Particles.\n");
    }

    return 0;
}

/**
 * @brief Creates and initializes a Particle Swarm.
 *
 * This function sets up a DMSwarm within the provided UserCtx structure, initializes
 * particle fields, and distributes particles across MPI processes. It ensures that
 * the number of particles is evenly divided among the available MPI ranks. If the total
 * number of particles isn't divisible by the number of processes, the remainder is distributed
 * to the first few ranks.
 *.
 *
 * @param[in,out] user          Pointer to the UserCtx structure containing the simulation context.
 * @param[in]     numParticles  Total number of particles to create across all MPI processes.
 * @param[in]     bboxlist      Pointer to an array of BoundingBox structures, one per rank.
 *
 * @param[in]     particlesPerProcess   
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 *
 * @note
 * - Ensure that `numParticles` is a positive integer.
 * - The `control.dat` file should contain necessary PETSc options.
 * - The `bboxlist` array should be properly populated before calling this function.
 */
PetscErrorCode CreateParticleSwarm(UserCtx *user, PetscInt numParticles, PetscInt *particlesPerProcess, BoundingBox *bboxlist) {
    PetscErrorCode ierr;                      // PETSc error handling variable
    PetscMPIInt rank, size;                   // Variables to store MPI rank and size
    PetscInt remainder = 0;                   // Remainder of particles after division
    
    // Validate input parameters
    if (numParticles <= 0) {
      LOG_ALLOW(GLOBAL,LOG_DEBUG, "CreateParticleSwarm - Number of particles must be positive. Given: %d\n", numParticles);
        return PETSC_ERR_ARG_OUTOFRANGE;
    }

    // Insert PETSc options from "control.dat" into the PETSc options database
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, "control.dat", PETSC_TRUE); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG,"CreateParticleSwarm - Inserted options from control.dat\n");
    
    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG, "CreateParticleSwarm - Domain dimensions: xMin=%.2f, xMax=%.2f,yMin=%.2f, yMax=%.2f,zMin=%.2f, zMax=%.2f \n", 
		   user->xMin,user->xMax,user->yMin,user->yMax, user->zMin,user->zMax);
    
    // Retrieve MPI rank and size
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_INFO, "CreateParticleSwarm - Rank %d out of %d processes.\n", rank, size);

    // Distribute particles among MPI processes
    ierr = DistributeParticles(numParticles, rank, size, particlesPerProcess, &remainder); CHKERRQ(ierr);

    // Initialize the DMSwarm - creates the swarm, sets the type and dimension
    ierr = InitializeSwarm(user); CHKERRQ(ierr);

    // Register particle fields (position, velocity, CellID, weight, etc.)
    ierr = RegisterParticleFields(user->swarm); CHKERRQ(ierr);

    // Set the local number of particles for this rank and additional buffer for particle migration
    ierr = DMSwarmSetLocalSizes(user->swarm, *particlesPerProcess, 4); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_INFO, "CreateParticleSwarm - Set local swarm size: %d particles.\n", *particlesPerProcess);

    // Optionally, LOG_ALLOW detailed DM info in debug mode
    if (get_log_level() == LOG_DEBUG && is_function_allowed(__func__)) {
      LOG_ALLOW(LOCAL,LOG_DEBUG,"CreateParticleSwarm - Viewing DMSwarm:\n");
        ierr = DMView(user->swarm, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    }

    LOG_ALLOW(LOCAL,LOG_INFO, "CreateParticleSwarm - Particle swarm creation and initialization complete.\n");

    return 0;
}

// This particle struct implementation is specifically for search and may not be used for later work! but can be extended!
// -----------------------------------------------------------------------------------------------------------------------
/**
 * @brief Initializes a Particle struct with data from DMSwarm fields.
 *
 * This helper function populates a Particle structure using data retrieved from DMSwarm fields.
 *
 * @param[in]     i            Index of the particle in the DMSwarm.
 * @param[in]     PIDs         Pointer to the array of particle IDs.
 * @param[in]     weights      Pointer to the array of particle weights.
 * @param[in]     positions    Pointer to the array of particle positions.
 * @param[in]     cellIndices  Pointer to the array of particle cell indices.
 * @param[out]    particle     Pointer to the Particle struct to initialize.
 *
 * @return PetscErrorCode  Returns `0` on success, non-zero on failure.
 */
PetscErrorCode InitializeParticle(PetscInt i, const PetscInt64 *PIDs, const PetscReal *weights,
                                         const PetscReal *positions, const PetscInt64 *cellIndices,
                                         Particle *particle) {
    PetscFunctionBeginUser;
    
    if (particle == NULL) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "InitializeParticle - Output Particle pointer is NULL. \n");
    }
    
    // logging the start of particle initialization
    LOG_ALLOW_SYNC(GLOBAL,LOG_INFO, "InitializeParticle - Initializing Particle [%d] with PID: %ld.\n", i, PIDs[i]);
    
    // Initialize PID
    particle->PID = PIDs[i];
    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG, "InitializeParticle - Particle [%d] PID set to: %ld.\n", i, particle->PID);
    
    // Initialize weights
    particle->weights.x = weights[3 * i];
    particle->weights.y = weights[3 * i + 1];
    particle->weights.z = weights[3 * i + 2];
    LOG_ALLOW(LOCAL,LOG_DEBUG, "InitializeParticle - Particle [%d] weights set to: (%.6f, %.6f, %.6f).\n", 
        i, particle->weights.x, particle->weights.y, particle->weights.z);
    
    // Initialize locations
    particle->loc.x = positions[3 * i];
    particle->loc.y = positions[3 * i + 1];
    particle->loc.z = positions[3 * i + 2];
    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG, "InitializeParticle - Particle [%d] location set to: (%.6f, %.6f, %.6f).\n", 
        i, particle->loc.x, particle->loc.y, particle->loc.z);
    
    // Initialize velocities (assuming default zero; modify if necessary)
    particle->vel.x = 0.0;
    particle->vel.y = 0.0;
    particle->vel.z = 0.0;
    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG,"InitializeParticle - Particle [%d] velocities initialized to zero.\n", i);
    
    // Initialize cell indices
    particle->cell[0] = cellIndices[3 * i];
    particle->cell[1] = cellIndices[3 * i + 1];
    particle->cell[2] = cellIndices[3 * i + 2];
    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG,"InitializeParticle - Particle [%d] cell indices set to: [%ld, %ld, %ld].\n", 
        i, particle->cell[0], particle->cell[1], particle->cell[2]);
    
    // logging the completion of particle initialization
    LOG_ALLOW_SYNC(GLOBAL,LOG_INFO,"InitializeParticle - Completed initialization of Particle [%d]. \n", i);
    
    PetscFunctionReturn(0);
}

/**
 * @brief Updates DMSwarm fields with data from a Particle struct.
 *
 * This helper function writes back the modified Particle data to the corresponding DMSwarm fields.
 *
 * @param[in] i            Index of the particle in the DMSwarm.
 * @param[in] particle     Pointer to the Particle struct containing updated data.
 * @param[in,out] weights  Pointer to the array of particle weights.
 * @param[in,out] cellIndices Pointer to the array of particle cell indices.
 *
 * @return PetscErrorCode  Returns `0` on success, non-zero on failure.
 */
PetscErrorCode UpdateSwarmFields(PetscInt i, const Particle *particle,
                                        PetscReal *weights, PetscInt64 *cellIndices) {
    PetscFunctionBeginUser;
    
    if (particle == NULL) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UpdateSwarmFields - Input Particle pointer is NULL.\n");
    }
    
    // logging the start of swarm fields update
    LOG_ALLOW_SYNC(GLOBAL,LOG_INFO,"Updating DMSwarm fields for Particle [%d].\n", i);
    
    // Update weights
    weights[3 * i]     = particle->weights.x;
    weights[3 * i + 1] = particle->weights.y;
    weights[3 * i + 2] = particle->weights.z;
    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG,"UpdateSwarmFields - Updated weights for Particle [%d]: (%.6f, %.6f, %.6f).\n", 
        i, weights[3 * i], weights[3 * i + 1], weights[3 * i + 2]);
    
    // Update cell indices
    cellIndices[3 * i]     = particle->cell[0];
    cellIndices[3 * i + 1] = particle->cell[1];
    cellIndices[3 * i + 2] = particle->cell[2];
    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG, "UpdateSwarmFields -  Updated cell indices for Particle [%d]: [%ld, %ld, %ld].\n", 
        i, cellIndices[3 * i], cellIndices[3 * i + 1], cellIndices[3 * i + 2]);
    
    // logging the completion of swarm fields update
    LOG_ALLOW_SYNC(GLOBAL,LOG_INFO,"UpdateSwarmFields  - Completed updating DMSwarm fields for Particle [%d].\n", i);
    
    PetscFunctionReturn(0);
}

/**
 * @brief Locates all particles within the grid and calculates their interpolation weights.
 * @ingroup ParticleLocation
 *
 * This function iterates through all particles currently local to this MPI rank.
 * For each particle, it first checks if the particle is within the rank's
 * pre-calculated bounding box (`user->bbox`). If it is, it calls the
 * `LocateParticleInGrid` function to perform the walking search.
 *
 * `LocateParticleInGrid` is responsible for finding the containing cell `(i,j,k)`
 * and calculating the corresponding interpolation weights `(w1,w2,w3)`. It updates
 * the `particle->cell` and `particle->weights` fields directly upon success.
 * If the search fails (particle not found within MAX_TRAVERSAL, goes out of bounds,
 * or gets stuck without resolution), `LocateParticleInGrid` sets the particle's
 * `cell` to `{-1,-1,-1}` and `weights` to `{0.0, 0.0, 0.0}`.
 *
 * After attempting location, this function updates the corresponding entries in the
 * DMSwarm's "DMSwarm_CellID" and "weight" fields using the potentially modified
 * data from the `particle` struct.
 *
 * @param[in] user Pointer to the UserCtx structure containing grid, swarm, and bounding box info.
 *
 * @return PetscErrorCode Returns `0` on success, non-zero on failure (e.g., errors accessing DMSwarm fields).
 *
 * @note Assumes `user->bbox` is correctly initialized for the local rank.
 * @note Assumes `InitializeParticle` correctly populates the temporary `particle` struct.
 * @note Assumes `UpdateSwarmFields` correctly writes data back to the DMSwarm.
 */
PetscErrorCode LocateAllParticlesInGrid(UserCtx *user) {
    PetscErrorCode ierr;
    PetscMPIInt rank, size;
    PetscInt localNumParticles;
    PetscReal *positions = NULL, *weights = NULL; // Pointers to DMSwarm data arrays
    PetscInt64 *cellIndices = NULL, *PIDs = NULL; // Pointers to DMSwarm data arrays
    DM swarm = user->swarm;                 // Convenience pointer to the swarm DM
    Particle particle;                      // Reusable temporary Particle struct for processing

    PetscFunctionBeginUser;
    LOG_FUNC_TIMER_BEGIN_EVENT(EVENT_walkingsearch, LOCAL);

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "LocateAllParticlesInGrid - Start on Rank %d/%d.\n", rank, size);

    // Optional barrier for debugging synchronization
    // ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

    // --- Access DMSwarm Data Arrays ---
    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "LocateAllParticlesInGrid - Number of local particles: %d.\n", localNumParticles);

    // Get direct pointers to the underlying data arrays for efficiency
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr); // Array to write weights back to
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIndices); CHKERRQ(ierr); // Array to write cell indices back to
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&PIDs); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "LocateAllParticlesInGrid - DMSwarm fields accessed successfully.\n");

    // --- Iterate over each local particle ---
    for (PetscInt i = 0; i < localNumParticles; ++i) {
        // Load current particle data into the temporary struct
        ierr = InitializeParticle(i, PIDs, weights, positions, cellIndices, &particle); CHKERRQ(ierr);

        LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, i, 10, "LocateAllParticlesInGrid - Processing Particle [%d]: PID=%lld.\n", i, particle.PID);

        // --- Coarse Check: Is particle within this rank's bounding box? ---
        // This is a quick check; particle could still be in a ghost cell managed by this rank.
        PetscBool particle_detected = IsParticleInsideBoundingBox(&(user->bbox), &particle);
        LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, i, 10, "LocateAllParticlesInGrid - Particle [%d] (PID %lld) inside local bbox: %s.\n",
                       i, particle.PID, particle_detected ? "YES" : "NO");

        if (particle_detected) {
            // --- Perform Detailed Location Search ---
            LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, i, 10, "LocateAllParticlesInGrid - Locating Particle [%d] (PID %lld) in grid...\n", i, particle.PID);

            // Call the walking search. This function will update particle.cell and particle.weights
            // internally if successful, or set them to -1 / 0 if it fails.
            ierr = LocateParticleInGrid(user, &particle); CHKERRQ(ierr); // Pass only user and particle struct

            // Log the outcome of the search for this particle
            if (particle.cell[0] >= 0) {
                 LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, i, 10,
                                "LocateAllParticlesInGrid - Particle [%d] (PID %lld) located/assigned to cell [%ld, %ld, %ld].\n",
                                i, particle.PID, particle.cell[0], particle.cell[1], particle.cell[2]);
            } else {
                 LOG_LOOP_ALLOW(LOCAL, LOG_WARNING, i, 1, // Log all failures
                                "LocateAllParticlesInGrid - Particle [%d] (PID %lld) FAILED TO LOCATE (CellID = -1).\n",
                                i, particle.PID);
            }
            // --- Weight calculation is now handled inside LocateParticleInGrid ---

        } else {
            // Particle was outside the local bounding box - mark as invalid for this rank
            LOG_LOOP_ALLOW(LOCAL, LOG_WARNING, i, 1,
                           "LocateAllParticlesInGrid - Particle [%d] (PID %lld) outside local bbox. Marking invalid (CellID = -1).\n",
                           i, particle.PID);
            particle.cell[0] = -1;
            particle.cell[1] = -1;
            particle.cell[2] = -1;
        } // end if (particle_detected)

        // --- Update DMSwarm Data ---
        // Write the potentially modified cell index and weights from the 'particle' struct
        // back into the main DMSwarm data arrays.
        ierr = UpdateSwarmFields(i, &particle, weights, cellIndices); CHKERRQ(ierr);

    } // --- End particle loop ---

    // --- Restore DMSwarm Data Arrays ---
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIndices); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&PIDs); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOCAL, LOG_INFO, "LocateAllParticlesInGrid - DMSwarm fields restored successfully on Rank %d.\n", rank);

    LOG_ALLOW(LOCAL, LOG_DEBUG, "LocateAllParticlesInGrid - Completed function on Rank %d.\n", rank);
    LOG_FUNC_TIMER_END_EVENT(EVENT_walkingsearch, LOCAL);
    PetscFunctionReturn(0);
}


/**
 * @brief Checks if a particle's location is within a specified bounding box.
 *
 * This function determines whether the given particle's location lies inside the provided bounding box.
 * It performs an axis-aligned bounding box (AABB) check by comparing the particle's coordinates to the
 * minimum and maximum coordinates of the bounding box in each dimension (x, y, z).
 *
 * logging statements are included to provide detailed information about the function's execution.
 *
 * @param[in]  bbox     Pointer to the BoundingBox structure containing minimum and maximum coordinates.
 * @param[in]  particle Pointer to the Particle structure containing the particle's location and identifier.
 *
 * @return PetscBool    Returns `PETSC_TRUE` if the particle is inside the bounding box, `PETSC_FALSE` otherwise.
 *
 * @note
 * - The function assumes that the `bbox` and `particle` pointers are valid and non-NULL.
 * - The function includes logging statements that start with the function name.
 * - The `LOG_ALLOW_SCOPE` variable is used to distinguish between `GLOBAL` and `LOCAL` LOG_ALLOW outputs.
 * - Be cautious when logging in performance-critical code sections, especially if the function is called frequently.
 */
PetscBool IsParticleInsideBoundingBox(const BoundingBox *bbox, const Particle *particle)
{
    // Function name for logging purposes
    const char *funcName = "IsParticleInsideBoundingBox";

    // Validate input pointers
    if (!bbox) {
        // LOG_ALLOW error message and return PETSC_FALSE
      LOG_ALLOW(LOCAL,LOG_ERROR, "%s: Error - 'bbox' pointer is NULL.", funcName);
        return PETSC_FALSE;
    }
    if (!particle) {
      LOG_ALLOW(LOCAL,LOG_ERROR,"%s: Error - 'particle' pointer is NULL.", funcName);
        return PETSC_FALSE;
    }

    // Extract particle location and bounding box coordinates
    const Cmpnts loc = particle->loc;
    const Cmpnts min_coords = bbox->min_coords;
    const Cmpnts max_coords = bbox->max_coords;

    // LOG_ALLOW the particle location and bounding box coordinates for debugging
    LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "%s: Particle PID %ld location: (%.6f, %.6f, %.6f).\n", funcName, particle->PID, loc.x, loc.y, loc.z);
    LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "%s: BoundingBox min_coords: (%.6f, %.6f, %.6f), max_coords: (%.6f, %.6f, %.6f).\n",
        funcName, min_coords.x, min_coords.y, min_coords.z, max_coords.x, max_coords.y, max_coords.z);

    // Check if the particle's location is within the bounding box
    if ((loc.x >= min_coords.x && loc.x <= max_coords.x) &&
        (loc.y >= min_coords.y && loc.y <= max_coords.y) &&
        (loc.z >= min_coords.z && loc.z <= max_coords.z)) {
        // Particle is inside the bounding box
        LOG_ALLOW_SYNC(LOCAL,LOG_DEBUG, "%s: Particle PID %ld is inside the bounding box.\n", funcName, particle->PID);
        return PETSC_TRUE;
    }

    // Particle is outside the bounding box
    LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG,"%s: Particle PID %ld is outside the bounding box.\n", funcName, particle->PID);
    return PETSC_FALSE;
}

/**
 * @brief Updates a particle's interpolation weights based on distances to cell faces.
 *
 * This function computes interpolation weights using distances to the six
 * cell faces (`d`) and updates the `weight` field of the provided particle.
 *
 * @param[in]  d        Pointer to an array of distances to the six cell faces.
 * @param[out] particle Pointer to the Particle structure whose weights are to be updated.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode UpdateParticleWeights(PetscReal *d, Particle *particle) {

    // Validate input pointers
    if (!d || !particle) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
                "UpdateParticleWeights - Null pointer argument (d or particle).");
    }


    // Validate distances
    for (PetscInt i = LEFT; i < NUM_FACES; i++) {
        if (d[i] <= INTERPOLATION_DISTANCE_TOLERANCE) {
            LOG_ALLOW_SYNC(GLOBAL, LOG_WARNING,
                "UpdateParticleWeights - face distance d[%d] = %f <= %f; "
                "clamping to 1e-14 to avoid zero/negative.\n",
                i, (double)d[i], INTERPOLATION_DISTANCE_TOLERANCE);
            d[i] = INTERPOLATION_DISTANCE_TOLERANCE;
        }
    }

    // LOG_ALLOW the input distances
    LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG,
        "UpdateParticleWeights - Calculating weights with distances: "
        "[LEFT=%f, RIGHT=%f, BOTTOM=%f, TOP=%f, FRONT=%f, BACK=%f].\n",
        d[LEFT], d[RIGHT], d[BOTTOM], d[TOP], d[FRONT], d[BACK]);

    // Compute and update the particle's weights
    particle->weights.x = d[LEFT] / (d[LEFT] + d[RIGHT]);
    particle->weights.y = d[BOTTOM] / (d[BOTTOM] + d[TOP]);
    particle->weights.z = d[FRONT] / (d[FRONT] + d[BACK]);

    // LOG_ALLOW the updated weights
    LOG_ALLOW_SYNC(LOCAL,LOG_DEBUG,
        "UpdateParticleWeights - Updated particle weights: x=%f, y=%f, z=%f.\n",
        particle->weights.x, particle->weights.y, particle->weights.z);

    return 0;
}

/**
 * @brief Perform particle swarm initialization, particle-grid interaction, and related operations.
 *
 * This function handles the following tasks:
 * 1. Initializes the particle swarm using the provided bounding box list (bboxlist) to determine initial placement
 *    if ParticleInitialization is 0.
 * 2. Locates particles within the computational grid.
 * 3. Updates particle positions based on grid interactions (if such logic exists elsewhere in the code).
 * 4. Interpolates particle velocities from grid points using trilinear interpolation.
 *
 * @param[in,out] user     Pointer to the UserCtx structure containing grid and particle swarm information.
 * @param[in]     np       Number of particles to initialize in the swarm.
 * @param[in]     bboxlist Pointer to an array of BoundingBox structures, one per MPI rank.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 *
 * @note
 * - Ensure that `np` (number of particles) is positive.
 * - The `bboxlist` array must be correctly computed and passed in before calling this function.
 * - If ParticleInitialization == 0, particles will be placed at the midpoint of the local bounding box.
 */
 PetscErrorCode InitializeParticleSwarm(UserCtx *user, PetscInt np, BoundingBox *bboxlist) {
    PetscErrorCode ierr;
    PetscInt particlesPerProcess = 0;         // Number of particles assigned to the local MPI process.
    PetscRandom randx,randy,randz;     // Random number generators[x,y,z]. (used if ParticleInitialization==1).
    PetscRandom rand_logic_i, rand_logic_j,rand_logic_k; // RNGs for Logical space. 
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting particle swarm Initialization with %d particles.\n", np);

    // Step 1: Create and initialize the particle swarm
    // Here we pass in the bboxlist, which will be used by CreateParticleSwarm() and subsequently
    // by AssignInitialProperties() to position particles if ParticleInitialization == 0.
    LOG_ALLOW(GLOBAL, LOG_INFO, "Creating particle swarm.\n");
    ierr = CreateParticleSwarm(user, np, &particlesPerProcess,bboxlist); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Particle swarm created successfully.\n");

    // Create the RNGs
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Initializing physical space RNGs.\n");
    ierr = InitializeRandomGenerators(user, &randx, &randy, &randz); CHKERRQ(ierr);
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Initializing logical space RNGs [0,1).\n");
    ierr = InitializeLogicalSpaceRNGs(&rand_logic_i, &rand_logic_j, &rand_logic_k); CHKERRQ(ierr);
    // Assign initial properties to particles
    // The bboxlist array is passed here so that if ParticleInitialization == 0,
    // particles can be placed at the midpoint of the local bounding box corresponding to this rank.
    ierr = AssignInitialPropertiesToSwarm(user, particlesPerProcess, &randx, &randy, &randz, &rand_logic_i,&rand_logic_j,&rand_logic_k,bboxlist); CHKERRQ(ierr);
    // Finalize swarm setup by destroying RNGs if ParticleInitialization == 1
    ierr = FinalizeSwarmSetup(&randx, &randy, &randz, &rand_logic_i, &rand_logic_j, &rand_logic_k); CHKERRQ(ierr);

    // Ensure all ranks complete before proceeding
    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, " Particles generated & initialized.\n");

    return 0;
}
