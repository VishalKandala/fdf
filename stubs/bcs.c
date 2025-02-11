/**
 *  @file bcs.c
 *  @brief Describes a general boundary condition application framework and several specific boundary conditions 
 */
///// Boundary Condition Application Framework //////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * @brief Retrieves the PETSc Vec and DM for a given field from the user context.
 *
 * This function is field-agnostic. Given a field name (e.g., "ucat" or "coor"),
 * it looks up the corresponding PETSc Vec and DM stored in the UserCtx structure.
 *
 * @param[in]  user      Pointer to the UserCtx structure.
 * @param[in]  fieldName Name of the field.
 * @param[out] vec       Pointer to the PETSc Vec associated with the field.
 * @param[out] dm        Pointer to the DM associated with the field.
 *
 * @return PetscErrorCode Returns 0 on success, nonzero on failure.
 */
PetscErrorCode GetFieldVec(UserCtx *user, const char *fieldName, Vec *vec, DM *dm)
{
    PetscFunctionBegin;
    if (strcmp(fieldName, "ucat") == 0) {
        if (!user->ucat)
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Field 'ucat' not found in UserCtx.");
        *vec = user->ucat;
        *dm = user->fda;
    } else if (strcmp(fieldName, "coor") == 0) {
        if (!user->coor)
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Field 'coor' not found in UserCtx.");
        *vec = user->coor;
        *dm = user->cda;
    } else {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Field '%s' not found in UserCtx.", fieldName);
    }
    PetscFunctionReturn(0);
}


/**
 * @brief Creates a new boundary condition registry.
 *
 * This function allocates and initializes a BCRegistry structure, which will hold
 * a linked list of BC objects representing different boundary conditions.
 *
 * @param[out] registry Pointer to the newly created BCRegistry.
 *
 * @return PetscErrorCode Returns 0 on success, nonzero on failure.
 */
PetscErrorCode BCRegistryCreate(BCRegistry **registry)
{
    PetscFunctionBegin;
    *registry = (BCRegistry*) malloc(sizeof(BCRegistry));
    if (!*registry)
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_MEM, "Unable to allocate BCRegistry");
    (*registry)->head = NULL;
    PetscFunctionReturn(0);
}


/**
 * @brief Adds a boundary condition (BC) object to the registry.
 *
 * This function inserts the given BC object at the beginning of the registry's
 * linked list.
 *
 * @param[in,out] registry Pointer to the BCRegistry.
 * @param[in]     bcObj    Pointer to the BC object to add.
 *
 * @return PetscErrorCode Returns 0 on success, nonzero on failure.
 */
PetscErrorCode BCRegistryAdd(BCRegistry *registry, BC *bcObj)
{
    PetscFunctionBegin;
    bcObj->next = registry->head;
    registry->head = bcObj;
    PetscFunctionReturn(0);
}


/**
 * @brief Destroys the boundary condition registry and frees all associated memory.
 *
 * This function iterates over all BC objects in the registry, calls each object's
 * destroy method (if provided), and frees the memory allocated for the registry.
 *
 * @param[in,out] registry Pointer to the BCRegistry to destroy.
 *
 * @return PetscErrorCode Returns 0 on success, nonzero on failure.
 */
PetscErrorCode BCRegistryDestroy(BCRegistry *registry)
{
    BC *current = registry->head;
    while (current) {
        BC *temp = current;
        current = current->next;
        if (temp->destroy)
            temp->destroy(temp);
    }
    free(registry);
    return 0;
}

/**
 * @brief Applies all registered boundary condition objects to a specified field.
 *
 * This function iterates over the BCRegistry linked list and calls the apply method
 * of each BC object for the given field name. The function is field-agnostic.
 *
 * @param[in] user      Pointer to the UserCtx structure.
 * @param[in] fieldName Name of the field to which the BCs will be applied.
 * @param[in] registry  Pointer to the BCRegistry containing the BC objects.
 *
 * @return PetscErrorCode Returns 0 on success, nonzero on failure.
 */
PetscErrorCode ApplyAllBCs(UserCtx *user, const char *fieldName, BCRegistry *registry)
{
    PetscErrorCode ierr;
    BC *current = registry->head;
    PetscFunctionBegin;
    while (current) {
        ierr = current->apply(user, fieldName, current); CHKERRQ(ierr);
        current = current->next;
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Sets boundary conditions for a given field by applying all registered BCs.
 *
 * This field-agnostic wrapper function creates a BC registry, registers the desired
 * BC objects, applies them to the specified field (by field name), and cleans up
 * the registry. This single function can be used for any field.
 *
 * @param[in] user      Pointer to the UserCtx structure.
 * @param[in] fieldName Name of the field on which to apply boundary conditions.
 *
 * @return PetscErrorCode Returns 0 on success, nonzero on failure.
 */
PetscErrorCode SetBoundaryConditions(UserCtx *user, const char *fieldName)
{
    PetscErrorCode ierr;
    BCRegistry *registry;
    PetscFunctionBegin;
    ierr = BCRegistryCreate(&registry); CHKERRQ(ierr);

    /* 
     * Registration of BC objects goes here.
     * For example, one can register BC_DirichletZero, BC_DirichletSin, BC_NeumannExtrapolate,
     * or any new BC types. Each BC object is created with its target face, type, pointers
     * to its apply and destroy functions, and its context. See individual BC implementations.
     */


/* Example: Creating a new BC using the template functions
   This example shows how to instantiate a BC that uses the template apply
   and destroy functions. You should replace the TODO comment with your
   specific computation.

   // Allocate memory for the new BC object
   BC *bc_new = (BC*) malloc(sizeof(BC));
   if (!bc_new) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_MEM, "Unable to allocate bc_new");
   
   // Set the target boundary face (e.g., BC_FACE_POS_Z)
   bc_new->face = BC_FACE_POS_Z;
   
   // Set the BC type (e.g., BC_TYPE_DIRICHLET or BC_TYPE_NEUMANN)
   bc_new->type = BC_TYPE_DIRICHLET;
   
   // Assign the apply method pointer using your template function (or a custom function)
   bc_new->apply = BC_Template_Apply;
   
   // Assign the destroy method pointer using your template destroy function
   bc_new->destroy = BC_Template_Destroy;
   
   // Allocate and initialize the context structure for your new BC.
   // For example, if your new BC needs a scaling factor, you might have:
   // typedef struct { int comp; PetscReal factor; } BC_ScaleCtx;
   // BC_ScaleCtx *scaleCtx = (BC_ScaleCtx*) malloc(sizeof(BC_ScaleCtx));
   // scaleCtx->comp = 0; // (e.g., x-component)
   // scaleCtx->factor = 0.5; // (e.g., scale by 0.5)
   // bc_new->ctx = scaleCtx;
   // For the template, you can simply set a dummy context if no parameters are needed:
   bc_new->ctx = NULL;  // Or allocate and assign a specific context as needed.
   
   // Add the new BC object to the registry
   ierr = BCRegistryAdd(registry, bc_new); CHKERRQ(ierr);
   
   // End of example for new BC creation using template functions.
*/



    ierr = ApplyAllBCs(user, fieldName, registry); CHKERRQ(ierr);
    ierr = BCRegistryDestroy(registry); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}



///////////// Boundary Conditions ///////////////////////////////////////

/**
 * @brief Template boundary condition apply function.
 *
 * This is a generic template for a BC apply function. Use this as a
 * starting point for implementing a specific boundary condition.
 *
 * The function performs the following steps:
 * 1. Retrieves the field vector and DM associated with the given field name.
 * 2. Obtains the local array from the DMDA.
 * 3. Checks if the current process's local subdomain touches the global boundary
 *    for the BC object's specified face.
 * 4. If so, fixes the appropriate index (e.g., i = info.xs for BC_FACE_NEG_X)
 *    and loops over the remaining indices.
 * 5. Computes the new boundary value using parameters stored in the BC object's
 *    context (self->ctx) and updates the array.
 * 6. Restores the local array and updates ghost cells.
 *
 * @param[in]  user      Pointer to the UserCtx structure.
 * @param[in]  fieldName Name of the field to which the BC is to be applied.
 * @param[in]  self      Pointer to the BC object containing BC parameters.
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
/*
PetscErrorCode BC_Template_Apply(UserCtx *user, const char *fieldName, BC *self)
{
    PetscErrorCode ierr;
    DMDALocalInfo info = user->info;
    Vec fieldVec;
    DM  fieldDM;
    Cmpnts ***array = NULL;
    PetscInt i, j, k;
    PetscBool applyBC = PETSC_FALSE;
    
    PetscFunctionBegin;
    // 1. Retrieve the field Vec and DM by the given field name 
    ierr = GetFieldVec(user, fieldName, &fieldVec, &fieldDM); CHKERRQ(ierr);
    
    // 2. Get the local array representation 
    ierr = DMDAVecGetArray(fieldDM, fieldVec, &array); CHKERRQ(ierr);
    
    // 3. Determine if this process's local subdomain touches the boundary
       for the BC object's specified face.  Adjust the fixed index accordingly.
       (Add cases for all faces as needed.)
    
    switch (self->face) {
        case BC_FACE_NEG_X:
            if (info.xs == 0) { applyBC = PETSC_TRUE; i = info.xs; }
            break;
        case BC_FACE_POS_X:
            if (info.xs + info.xm == info.mx) { applyBC = PETSC_TRUE; i = info.xs + info.xm - 1; }
            break;
        case BC_FACE_NEG_Y:
            if (info.ys == 0) { applyBC = PETSC_TRUE; j = info.ys; }
            break;
        case BC_FACE_POS_Y:
            if (info.ys + info.ym == info.my) { applyBC = PETSC_TRUE; j = info.ys + info.ym - 1; }
            break;
        case BC_FACE_NEG_Z:
            if (info.zs == 0) { applyBC = PETSC_TRUE; k = info.zs; }
            break;
        case BC_FACE_POS_Z:
            if (info.zs + info.zm == info.mz) { applyBC = PETSC_TRUE; k = info.zs + info.zm - 1; }
            break;
        default:
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "BC_Template_Apply: Unknown face.");
    }
    
    // 4. If the current process touches the boundary, loop over the remaining indices.
       For example, for BC_FACE_NEG_X, fix i and loop over j and k.
       Insert your computation to update the boundary node value below.
    
    if (applyBC) {
        if (self->face == BC_FACE_NEG_X || self->face == BC_FACE_POS_X) {
            for (k = info.zs; k < info.zs + info.zm; k++) {
                for (j = info.ys; j < info.ys + info.ym; j++) {
                    // Example: retrieve the current value, compute the new value,
                    // and update the field.
                    // Use self->ctx to access BC-specific parameters.
                    
                    Cmpnts newVal = array[k][j][i];
                    // TODO: Insert computation for newVal based on self->ctx.
                    array[k][j][i] = newVal;
                }
            }
        }
        // Repeat similar loops for Y and Z faces as needed.
    }
    
    // 5. Restore the local array 
    ierr = DMDAVecRestoreArray(fieldDM, fieldVec, &array); CHKERRQ(ierr);
    
    // 6. Update ghost cells to ensure all processes have consistent boundary values 
    ierr = DMGlobalToLocalBegin(fieldDM, fieldVec, INSERT_VALUES, fieldVec); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(fieldDM, fieldVec, INSERT_VALUES, fieldVec); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}
*/

/**
 * @brief Template boundary condition destroy function.
 *
 * This function is a template for cleaning up a BC object.
 * It should free any dynamically allocated memory in the BC object's context,
 * and then free the BC object itself.
 *
 * @param[in] self Pointer to the BC object to destroy.
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
/*
PetscErrorCode BC_Template_Destroy(BC *self)
{
    PetscFunctionBegin;
    // Free any memory allocated in the context, if applicable 
    if (self->ctx) {
        free(self->ctx);
    }
   // Free the BC object itself 
    free(self);
    PetscFunctionReturn(0);
}
*/

/**
 * @brief Applies a Dirichlet sine boundary condition.
 *
 * This function sets a specified component of the field at a boundary to
 * sin(coordinate), where the coordinate is obtained from the coordinate field.
 *
 * The function works only on a designated boundary (e.g., the positive y-face).
 * It retrieves the target field (by fieldName) and its DM, accesses the local array,
 * and for nodes on the boundary, reads the coordinate value (using the coordinate DM/Vec)
 * and sets the corresponding field component to sin(coordinate).
 *
 * @param[in]  user      Pointer to the UserCtx structure.
 * @param[in]  fieldName Name of the field to update (e.g., "ucat").
 * @param[in]  self      Pointer to the BC object containing BC parameters.
 *
 * @return PetscErrorCode Returns 0 on success, nonzero on failure.
 */
PetscErrorCode BC_DirichletSin_Apply(UserCtx *user, const char *fieldName, BC *self)
{
    PetscErrorCode ierr;
    DMDALocalInfo info = user->info;
    Vec fieldVec;
    DM  fieldDM;
    Cmpnts ***array = NULL;
    /* The context holds the component index and pointers to the coordinate DM and Vec */
    BC_SinCtx *ctx = (BC_SinCtx*) self->ctx;
    PetscInt i, j, k;
    PetscBool applyBC = PETSC_FALSE;
    
    PetscFunctionBegin;
    /* Retrieve the target field vector and its DM by name */
    ierr = GetFieldVec(user, fieldName, &fieldVec, &fieldDM); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(fieldDM, fieldVec, &array); CHKERRQ(ierr);
    
    /* For this sine BC, we assume it applies on a positive face.
       For example, if self->face == BC_FACE_POS_Y, check that the local subdomain touches the global positive y boundary.
    */
    if (self->face == BC_FACE_POS_Y) {
        if (info.ys + info.ym == info.my) { applyBC = PETSC_TRUE; j = info.ys + info.ym - 1; }
    }
    /* Additional cases for other faces could be added similarly */
    
    if (!applyBC) {
        ierr = DMDAVecRestoreArray(fieldDM, fieldVec, &array); CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }
    
    /* Loop over the boundary nodes on the fixed face.
       In this example, for the positive y-face, we loop over i and k while j is fixed.
    */
    for (k = info.zs; k < info.zs + info.zm; k++) {
        for (i = info.xs; i < info.xs + info.xm; i++) {
            /* Obtain the corresponding coordinate from the coordinate field */
            Vec coorVec;
            DM  coorDM;
            Cmpnts ***carray;
            ierr = GetFieldVec(user, "coor", &coorVec, &coorDM); CHKERRQ(ierr);
            ierr = DMDAVecGetArray(coorDM, coorVec, &carray); CHKERRQ(ierr);
            /* Here we assume that the coordinate of interest is the y-component */
            PetscReal coord_val = carray[k][j][i].y;
            ierr = DMDAVecRestoreArray(coorDM, coorVec, &carray); CHKERRQ(ierr);
            
            /* Compute the new value: set the designated component to sin(coordinate) */
            Cmpnts newVal = array[k][j][i];  /* start with the current value */
            if (ctx->comp == 0) newVal.x = sin(coord_val);
            else if (ctx->comp == 1) newVal.y = sin(coord_val);
            else if (ctx->comp == 2) newVal.z = sin(coord_val);
            array[k][j][i] = newVal;
        }
    }
    
    /* Restore the local array and update ghost cells */
    ierr = DMDAVecRestoreArray(fieldDM, fieldVec, &array); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(fieldDM, fieldVec, INSERT_VALUES, fieldVec); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(fieldDM, fieldVec, INSERT_VALUES, fieldVec); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}
