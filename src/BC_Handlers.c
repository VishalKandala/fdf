#include "BC_Handlers.h"     // The header that declares this file's "constructor" functions

//================================================================================
//
//               HANDLER IMPLEMENTATION: NO-SLIP WALL
//               (Corresponds to BC_HANDLER_WALL_NOSLIP)
//
// This handler implements a stationary, impenetrable wall where the fluid
// velocity is zero (no-slip condition). It is one of the simplest but most
// common boundary conditions.
//
//================================================================================

/**
 * @brief (Handler Action) Applies the no-slip wall condition to a specified face.
 *
 * This function is the core implementation for the no-slip wall. It is called by the
 * BoundarySystem during the Apply phase of each time step. Its responsibilities are:
 *   1. Check if the current MPI rank owns any part of the face to be processed.
 *   2. If so, iterate over the local portion of that face.
 *   3. For each boundary cell face, set the normal contravariant velocity (flux) to zero.
 *   4. Set the Cartesian ghost cell velocity to enforce the no-slip condition
 *      (e.g., u_ghost = -u_interior), which is a common second-order accurate
 *      approximation for a wall located exactly on the cell face.
 *
 * @param self A pointer to the BoundaryCondition object. This simple handler does not
 *             use this parameter, but it is required by the standard interface.
 * @param ctx  A pointer to the BCContext, which provides access to the UserCtx (containing
 *             the vectors to be modified) and the ID of the face to be processed.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode Apply_WallNoSlip(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    PetscBool      can_service;
    PetscMPIInt    rank;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    
    // This unused-variable pragma silences compiler warnings for simple handlers
    // that don't need their own 'self' pointer.
    (void)self;
    
    PetscFunctionBeginUser;

    // Step 1: Check if this MPI rank has any work to do for this face.
    // This is a critical efficiency step to avoid unnecessary work and memory access.
    // It uses a utility function that understands the parallel grid decomposition.
    ierr = CanRankServiceFace(&user->info, face_id, &can_service); CHKERRQ(ierr);
    if (!can_service) {
        PetscFunctionReturn(0); // This rank does not own this boundary face. Exit silently.
    }

    // If we proceed, this rank has work to do. Log the action.
    const char *face_name = BCFaceToString(face_id);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Apply_WallNoSlip: Rank %d applying condition to Face %d (%s).\n",rank, face_id, face_name);

    // Step 2: Get safe access to the local PETSc vector arrays.
    // We get the local vectors (Ucat, Ucont) because ghost cell data is only
    // guaranteed to be correct on the local representation after a ghost update.
    DMDALocalInfo  *info = &user->info;
    Cmpnts       ***ucat, ***ucont;
    ierr = DMDAVecGetArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);

    // Step 3: Apply the no-slip condition based on which face is being processed.
    // The loop bounds are for the *owned* nodes on this rank (xs to xe, etc.).
    // The logic inside correctly accesses the ghost nodes or face-normal fluxes.
    PetscInt i, j, k;
    PetscInt xs = info->xs, xe = info->xs + info->xm;
    PetscInt ys = info->ys, ye = info->ys + info->ym;
    PetscInt zs = info->zs, ze = info->zs + info->zm;
    PetscInt mx = info->mx, my = info->my, mz = info->mz; // Global dimensions

    switch (face_id) {
        case BC_FACE_NEG_X: // -Xi face at global index i=0
            i = 0;
            for (k = zs; k < ze; k++) {
                for (j = ys; j < ye; j++) {
                    ucont[k][j][i].x = 0.0; // Set normal contravariant flux to zero.
                    // Set Cartesian ghost cell velocity for no-slip: u_ghost(i) = -u_interior(i+1)
                    ucat[k][j][i].x = -ucat[k][j][i+1].x;
                    ucat[k][j][i].y = -ucat[k][j][i+1].y;
                    ucat[k][j][i].z = -ucat[k][j][i+1].z;
                }
            }
            break;

        case BC_FACE_POS_X: // +Xi face at global index i = mx-1
            i = mx - 1;
            for (k = zs; k < ze; k++) {
                for (j = ys; j < ye; j++) {
                    ucont[k][j][i-1].x = 0.0;
                    // u_ghost(i) = -u_interior(i-1)
                    ucat[k][j][i].x = -ucat[k][j][i-1].x;
                    ucat[k][j][i].y = -ucat[k][j][i-1].y;
                    ucat[k][j][i].z = -ucat[k][j][i-1].z;
                }
            }
            break;

        case BC_FACE_NEG_Y: // -Eta face at global index j=0
            j = 0;
            for (k = zs; k < ze; k++) {
                for (i = xs; i < xe; i++) {
                    ucont[k][j][i].y = 0.0;
                    ucat[k][j][i].x = -ucat[k][j+1][i].x;
                    ucat[k][j][i].y = -ucat[k][j+1][i].y;
                    ucat[k][j][i].z = -ucat[k][j+1][i].z;
                }
            }
            break;

        case BC_FACE_POS_Y: // +Eta face at global index j = my-1
            j = my - 1;
            for (k = zs; k < ze; k++) {
                for (i = xs; i < xe; i++) {
                    ucont[k][j-1][i].y = 0.0;
                    ucat[k][j][i].x = -ucat[k][j-1][i].x;
                    ucat[k][j][i].y = -ucat[k][j-1][i].y;
                    ucat[k][j][i].z = -ucat[k][j-1][i].z;
                }
            }
            break;
            
        case BC_FACE_NEG_Z: // -Zeta face at global index k=0
            k = 0;
            for (j = ys; j < ye; j++) {
                for (i = xs; i < xe; i++) {
                    ucont[k][j][i].z = 0.0;
                    ucat[k][j][i].x = -ucat[k+1][j][i].x;
                    ucat[k][j][i].y = -ucat[k+1][j][i].y;
                    ucat[k][j][i].z = -ucat[k+1][j][i].z;
                }
            }
            break;

        case BC_FACE_POS_Z: // +Zeta face at global index k = mz-1
            k = mz - 1;
            for (j = ys; j < ye; j++) {
                for (i = xs; i < xe; i++) {
                    ucont[k-1][j][i].z = 0.0;
                    ucat[k][j][i].x = -ucat[k-1][j][i].x;
                    ucat[k][j][i].y = -ucat[k-1][j][i].y;
                    ucat[k][j][i].z = -ucat[k-1][j][i].z;
                }
            }
            break;
    }

    // Step 4: Restore safe access to the PETSc vector arrays.
    ierr = DMDAVecRestoreArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/**
 * @brief (Handler Constructor) Populates a BoundaryCondition object with No-Slip Wall behavior.
 *
 * This function is called by the BoundarySystem factory (`BoundaryCondition_Create` in
 * `Boundaries.c`). It "constructs" a no-slip wall handler by setting the function
 * pointers in the provided BoundaryCondition struct to point to the static functions
 * defined in this file.
 *
 * A no-slip wall is simple and requires only the `Apply` method:
 *  - It needs no special initialization (`Initialize` is NULL).
 *  - It does not contribute to the global mass balance (`PreStep` is NULL).
 *  - It has a specific `Apply` function to enforce zero velocity.
 *  - It allocates no private data, so `Destroy` is NULL.
 *
 * @param bc A pointer to the generic BoundaryCondition object to be configured.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode Create_WallNoSlip(BoundaryCondition *bc)
{
    PetscFunctionBeginUser;
    if (!bc) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input BoundaryCondition object is NULL in Create_WallNoSlip");

    // Assign the appropriate function pointers for this handler type.
    // This is the essence of the "Strategy" pattern in C.
    bc->Initialize = NULL;
    bc->PreStep    = NULL;
    bc->Apply      = Apply_WallNoSlip; // The ONLY action this handler needs to perform.
    bc->Destroy    = NULL;
    
    // No private data struct is needed for this handler, so bc->data remains NULL.

    PetscFunctionReturn(0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////

//================================================================================
//
//          HANDLER IMPLEMENTATION: CONSTANT VELOCITY INLET
//          (Corresponds to BC_HANDLER_INLET_CONSTANT_VELOCITY)
//
// This handler implements an inlet with a prescribed, uniform Cartesian velocity.
//
//================================================================================

// --- 1. FORWARD DECLARATIONS for this handler's static methods ---
// This tells the compiler that these functions exist before they are used.
// It resolves the "declared implicitly" warning.
static PetscErrorCode Initialize_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode PreStep_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx, PetscReal *in, PetscReal *out);
static PetscErrorCode Apply_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx);
static PetscErrorCode Destroy_InletConstantVelocity(BoundaryCondition *self);

/**
 * @brief Private data structure for the Constant Velocity Inlet handler.
 *
 * This struct holds the specific parameters needed for this handler, which are
 * parsed from the bcs.dat file during initialization.
 */
typedef struct {
    Cmpnts specified_velocity; // The desired Cartesian velocity (vx, vy, vz)
} InletConstantData;


/**
 * @brief (Handler Action) Sets the initial state of the boundary face.
 *
 * This function is called once by BoundarySystem_Create. It performs two tasks:
 * 1. Parses the `params` list (from bcs.dat) to find the 'vx', 'vy', and 'vz'
 *    parameters and stores them in its private `data` struct.
 * 2. Sets the initial Ucont and Ucat values on the boundary face to reflect
 *    this constant velocity.
 *
 * @param self The BoundaryCondition object for this handler.
 * @param ctx  The context, providing access to UserCtx and face_id.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode Initialize_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    InletConstantData *data = (InletConstantData*)self->data;
    
    PetscFunctionBeginUser;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Initialize_InletConstantVelocity: Initializing handler for Face %d. \n", face_id);

    // 1. Parse parameters from the linked list stored in the face configuration.
    data->specified_velocity = (Cmpnts){0.0, 0.0, 0.0}; // Default to zero
    for (BC_Param *p = user->boundary_faces[face_id].params; p; p = p->next) {
        if (strcasecmp(p->key, "vx") == 0) data->specified_velocity.x = atof(p->value);
        else if (strcasecmp(p->key, "vy") == 0) data->specified_velocity.y = atof(p->value);
        else if (strcasecmp(p->key, "vz") == 0) data->specified_velocity.z = atof(p->value);
    }
    LOG_ALLOW(LOCAL, LOG_INFO, "  Inlet Face %d configured with velocity (vx,vy,vz) = (%.2f, %.2f, %.2f) \n",
              face_id, data->specified_velocity.x, data->specified_velocity.y, data->specified_velocity.z);

    // 2. Set the initial boundary state. We can simply call the Apply function to do this,
    //    as the logic is the same for the initial state and subsequent steps for this handler.
    ierr = Apply_InletConstantVelocity(self, ctx); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/**
 * @brief (Handler Action) Calculates the target inflow flux for the PreStep phase.
 *
 * This function calculates the total volumetric flux that *should* be entering
 * through the portion of the inlet face owned by this MPI rank. It does this by
 * dotting the specified Cartesian velocity with the face-normal area vectors.
 *
 * @param self The BoundaryCondition object for this handler.
 * @param ctx  The context, providing UserCtx and face_id.
 * @param[out] local_inflow_contribution  The calculated inflow flux is added to this value.
 * @param[out] local_outflow_contribution This handler does not produce outflow, so this is unused.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode PreStep_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx, PetscReal *local_inflow_contribution, PetscReal *local_outflow_contribution)
{
    // ... (This function would calculate the flux: V_spec · Area_vector) ...
    // For Phase 2, we can leave this as a placeholder, as it's only needed
    // when we have a mass-conserving outlet that needs to know the target inflow.
    (void)self; (void)ctx; (void)local_inflow_contribution; (void)local_outflow_contribution;
    PetscFunctionBeginUser;
    PetscFunctionReturn(0);
}

/**
 * @brief (Handler Action) Applies the constant velocity inlet condition to a specified face.
 *
 * This function is the core "workhorse" for the constant velocity inlet. It is called
 * by the BoundarySystem during the Apply phase of each time step. Its responsibility
 * is to enforce the pre-configured velocity on its assigned face. It does this by:
 *
 * 1.  Setting the ghost-cell Cartesian velocity (`Ucat`) directly to the specified
 *     velocity vector (e.g., `ucat[k][j][0] = {vx, vy, vz}`). This provides a
 *     Dirichlet condition for any calculations using Cartesian velocity.
 *
 * 2.  Calculating the corresponding contravariant velocity components (`Ucont`). This
 *     is done by projecting the specified Cartesian velocity vector onto the
 *     contravariant basis vectors at each point on the face. For example, the first
 *     contravariant component is U¹ = v_cartesian ⋅ g¹ = (vx,vy,vz) ⋅ (csi.x,csi.y,csi.z).
 *     This ensures that the volumetric flux through each face is consistent with the
 *     specified Cartesian velocity.
 *
 * @param self The BoundaryCondition object for this handler, which contains the
 *             private `data` struct holding the configured velocity.
 * @param ctx  The context, providing access to UserCtx (for PETSc vectors and grid info)
 *             and the ID of the face to be processed.
 * @return PetscErrorCode 0 on success.
 */
static PetscErrorCode Apply_InletConstantVelocity(BoundaryCondition *self, BCContext *ctx)
{
    PetscErrorCode ierr;
    UserCtx*       user = ctx->user;
    BCFace         face_id = ctx->face_id;
    InletConstantData *data = (InletConstantData*)self->data; // Cast private data
    PetscBool      can_service;
    PetscMPIInt    rank;
    
    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // Step 1: Check if this MPI rank owns any part of the face to be processed.
    ierr = CanRankServiceFace(&user->info, face_id, &can_service); CHKERRQ(ierr);
    if (!can_service) {
        PetscFunctionReturn(0); // This rank has no work to do for this face.
    }

    const char* face_name = BCFaceToString(face_id);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Apply_InletConstantVelocity: Rank %d applying condition to Face %d (%s).\n", rank, face_id, face_name);

    // Step 2: Get access to the necessary PETSc vector arrays.
    DMDALocalInfo  *info = &user->info;
    Cmpnts       ***ucat, ***ucont;
    Cmpnts       ***l_csi, ***l_eta, ***l_zet;
    ierr = DMDAVecGetArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, &l_csi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, &l_eta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, &l_zet); CHKERRQ(ierr);

    const Cmpnts *v_spec = &data->specified_velocity; // Create a convenient shortcut

    // Step 3: Loop over the owned nodes on the specified face and apply the BC.
    PetscInt i, j, k;
    const PetscInt xs = info->xs, xe = info->xs + info->xm;
    const PetscInt ys = info->ys, ye = info->ys + info->ym;
    const PetscInt zs = info->zs, ze = info->zs + info->zm;
    const PetscInt mx = info->mx, my = info->my, mz = info->mz;

    // The loops iterate over the entire local ownership range. The `if` conditions
    // inside ensure that the logic is only applied to the nodes on the correct boundary face.
    for (k = zs; k < ze; k++) {
        for (j = ys; j < ye; j++) {
            for (i = xs; i < xe; i++) {

                if (i == 0 && face_id == BC_FACE_NEG_X) {
                    ucat[k][j][i] = *v_spec;
                    ucont[k][j][i].x = v_spec->x * l_csi[k][j][i].x + v_spec->y * l_csi[k][j][i].y + v_spec->z * l_csi[k][j][i].z;
                }
                else if (i == mx - 1 && face_id == BC_FACE_POS_X) {
                    ucat[k][j][i] = *v_spec;
                    // Note: Ucont.x is on the i-face, so for the max face, it's at index i-1
                    ucont[k][j][i-1].x = v_spec->x * l_csi[k][j][i-1].x + v_spec->y * l_csi[k][j][i-1].y + v_spec->z * l_csi[k][j][i-1].z;
                }
                else if (j == 0 && face_id == BC_FACE_NEG_Y) {
                    ucat[k][j][i] = *v_spec;
                    ucont[k][j][i].y = v_spec->x * l_eta[k][j][i].x + v_spec->y * l_eta[k][j][i].y + v_spec->z * l_eta[k][j][i].z;
                }
                else if (j == my - 1 && face_id == BC_FACE_POS_Y) {
                    ucat[k][j][i] = *v_spec;
                    ucont[k][j-1][i].y = v_spec->x * l_eta[k][j-1][i].x + v_spec->y * l_eta[k][j-1][i].y + v_spec->z * l_eta[k][j-1][i].z;
                }
                else if (k == 0 && face_id == BC_FACE_NEG_Z) {
                    ucat[k][j][i] = *v_spec;
                    ucont[k][j][i].z = v_spec->x * l_zet[k][j][i].x + v_spec->y * l_zet[k][j][i].y + v_spec->z * l_zet[k][j][i].z;
                }
                else if (k == mz - 1 && face_id == BC_FACE_POS_Z) {
                    ucat[k][j][i] = *v_spec;
                    ucont[k-1][j][i].z = v_spec->x * l_zet[k-1][j][i].x + v_spec->y * l_zet[k-1][j][i].y + v_spec->z * l_zet[k-1][j][i].z;
                }
            }
        }
    }
    
    // Step 4: Restore safe access to the PETSc vector arrays.
    ierr = DMDAVecRestoreArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Ucont, &ucont); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, &l_csi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, &l_eta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, &l_zet); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

/**
 * @brief (Handler Destructor) Frees memory allocated by the Constant Velocity Inlet handler.
 */
static PetscErrorCode Destroy_InletConstantVelocity(BoundaryCondition *self)
{
    PetscFunctionBeginUser;
    if (self && self->data) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Destroy_InletConstantVelocity: Freeing private data struct.");
        PetscFree(self->data);
        self->data = NULL;
    }
    PetscFunctionReturn(0);
}


/**
 * @brief (Handler Constructor) Populates a BoundaryCondition object with Constant Velocity Inlet behavior.
 */
PetscErrorCode Create_InletConstantVelocity(BoundaryCondition *bc)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    if (!bc) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input BoundaryCondition object is NULL");

    // Allocate the private data structure for this handler
    InletConstantData *data = NULL;
    ierr = PetscMalloc1(1, &data); CHKERRQ(ierr);
    bc->data = (void*)data;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Create_InletConstantVelocity: Allocated InletConstantData (%zu bytes) at %p.\n",
              sizeof(*data), (void*)data);
    
    // Assign the function pointers for this handler type
    bc->Initialize = Initialize_InletConstantVelocity;
    bc->PreStep    = PreStep_InletConstantVelocity;
    bc->Apply      = Apply_InletConstantVelocity;
    bc->Destroy    = Destroy_InletConstantVelocity;
    
    PetscFunctionReturn(0);
}

/*****************************************************************************
 *  NOGRAD – copy-ghost handler
 *
 *  Behaviour
 *  ---------
 *    • Works only on the face(s) where it is prescribed in bcs.dat.
 *    • Copies the entire *first interior* node/face layer onto the ghost
 *      layer for both Cartesian velocity  Ucat  and contravariant velocity
 *      Ucont (normal component only – i.e. the one that lives on that face).
 *    • No Initialise() or PreStep() are required.
 *
 *  Integration
 *  -----------
 *    Called by BoundarySystem_ExecuteStep() after contra2cart().
 *****************************************************************************/

/* --- no private data is needed, keep an empty struct for future growth --- */
typedef struct { int dummy; } NgData;

/* ------------------------------------------------------------------------- */
static PetscErrorCode Apply_NogradCopyGhost(BoundaryCondition *self,
                                            BCContext         *ctx)
{
    PetscErrorCode ierr;
    UserCtx       *u   = ctx->user;
    DMDALocalInfo *inf = &u->info;

    /* Arrays */
    Cmpnts ***ucat, ***ucont;
    ierr = DMDAVecGetArray(u->fda, u->Ucat , &ucat ); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(u->fda, u->Ucont, &ucont); CHKERRQ(ierr);

    const PetscInt xs = inf->xs, xe = inf->xs + inf->xm;
    const PetscInt ys = inf->ys, ye = inf->ys + inf->ym;
    const PetscInt zs = inf->zs, ze = inf->zs + inf->zm;
    const PetscInt mx = inf->mx, my = inf->my, mz = inf->mz;

    PetscInt i,j,k;

    PetscInt lxs  = (xs == 0) ? xs + 1 : xs;
    PetscInt lxe   = (xe == mx) ? xe - 1 : xe;

    PetscInt lys  = (ys == 0) ? ys + 1 : ys;
    PetscInt lye   = (ye == my) ? ye - 1 : ye;

    PetscInt lzs  = (zs == 0) ? zs + 1 : zs;
    PetscInt lze   = (ze == mz) ? ze - 1 : ze;

    
    switch (ctx->face_id)
    {
    /* ------------------------------------------------------------------ */
    case BC_FACE_NEG_X:   /* i = 0 copies from i = 1 */
        if (xs == 0)
        {
            i = 0;
            for (k=lzs;k<lze;k++)
              for (j=lys;j<lye;j++)
              {
                  ucat [k][j][i]   = ucat [k][j][i+1];
                  ucont[k][j][i].x = ucont[k][j][i+1].x;    /* flux normal to face */
              }
        }
        break;

    case BC_FACE_POS_X:   /* i = mx-1 copies from i = mx-2 */
        if (xe == mx)
        {
            i = mx-1;
            for (k=lzs;k<lze;k++)
              for (j=lys;j<lye;j++)
              {
                  ucat [k][j][i]     = ucat [k][j][i-1];
		  // Add a guard to prevent out-of-bounds read on thin domains
		  if (mx >= 3) {
		    ucont[k][j][i].x = ucont[k][j][i-1].x;
		  }
              }
        }
        break;

    /* ------------------------------------------------------------------ */
    case BC_FACE_NEG_Y:   /* j = 0 from j = 1 */
        if (ys == 0)
        {
            j = 0;
            for (k=lzs;k<lze;k++)
              for (i=lxs;i<lxe;i++)
              {
                  ucat [k][j][i]   = ucat [k][j+1][i];
                  ucont[k][j][i].y = ucont[k][j+1][i].y;
              }
        }
        break;

    case BC_FACE_POS_Y:   /* j = my-1 from j = my-2 */
        if (ye == my)
        {
            j = my-1;
            for (k=zs;k<ze;k++)
              for (i=xs;i<xe;i++)
              {
                  ucat [k][j][i]     = ucat [k][j-1][i];
		  // Add a guard to prevent out-of-bounds read on thin domains
		  if (my >= 3) {
		    ucont[k][j][i].y = ucont[k][j-1][i].y;
		  }
	      }
        }
        break;

    /* ------------------------------------------------------------------ */
    case BC_FACE_NEG_Z:   /* k = 0 from k = 1 */
        if (zs == 0)
        {
            k = 0;
            for (j=lys;j<lye;j++)
              for (i=lxs;i<lxe;i++)
              {
                  ucat [k][j][i]   = ucat [k+1][j][i];
                  ucont[k][j][i].z = ucont[k+1][j][i].z;
              }
        }
        break;

    case BC_FACE_POS_Z:   /* k = mz-1 from k = mz-2 */
        if (ze == mz)
        {
            k = mz-1;
            for (j=lys;j<lye;j++)
              for (i=lxs;i<lxe;i++)
              {
                  ucat [k][j][i]     = ucat [k-1][j][i];
		  // Add a guard to prevent out-of-bounds read on thin domains
		  if (mz >= 3) {
		    ucont[k][j][i].z = ucont[k-1][j][i].z;
		  }
              }
        }
        break;
    }

    ierr = DMDAVecRestoreArray(u->fda, u->Ucat , &ucat ); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(u->fda, u->Ucont, &ucont); CHKERRQ(ierr);
    return 0;
}

/* ------------------------------------------------------------------------- */
static PetscErrorCode Destroy_NogradCopyGhost(BoundaryCondition *self)
{
    PetscFree(self->data);      /* nothing in it yet, but keep it symmetric  */
    self->data = NULL;
    return 0;
}

/* ------------------------------------------------------------------------- */
PetscErrorCode Create_NogradCopyGhost(BoundaryCondition *bc)
{
    PetscErrorCode ierr;
    ierr = PetscMalloc1(1,&bc->data); CHKERRQ(ierr);   /* allocate empty struct */

    bc->Initialize = NULL;       /* no parameters to parse */
    bc->PreStep    = NULL;       /* no flux bookkeeping    */
    bc->Apply      = Apply_NogradCopyGhost;
    bc->Destroy    = Destroy_NogradCopyGhost;
    return 0;
}
