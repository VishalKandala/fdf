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
 * @brief Initializes basic particle properties: position, particle ID, and cell ID.
 *
 * This function sets the position field using either random values (if ParticleInitialization==1)
 * or at the z=0 plane. It also assigns a unique particle ID and initializes
 * cell IDs to -1.
 *
 * @param[in,out] user               Pointer to the UserCtx structure containing the swarm.
 * @param[in]     particlesPerProcess Number of particles assigned to this MPI process.
 * @param[in]     randx              Random number generator for the x-coordinate.
 * @param[in]     randy              Random number generator for the y-coordinate.
 * @param[in]     randz              Random number generator for the z-coordinate.
 * @param[in]     bboxlist           Array of BoundingBox structures (one per MPI rank).
 *
 * @return PetscErrorCode            Returns 0 on success, non-zero on failure.
 */
static PetscErrorCode InitializeParticleBasicProperties(UserCtx *user,
                                                   PetscInt particlesPerProcess,
                                                   PetscRandom *randx, PetscRandom *randy, PetscRandom *randz,
                                                   BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    DM swarm;
    PetscReal *positions = NULL;
    PetscInt64 *particleIDs = NULL, *cellIDs = NULL;
    PetscMPIInt rank;
    const Cmpnts ***coor;
    Vec Coor;
    PetscReal *pos_phys = NULL;
    PetscFunctionBeginUser;

    // Validate input pointers
    if (!user || !randx || !randy || !randz || !bboxlist) {
        LOG_ALLOW(GLOBAL, LOG_ERROR, "InitializeParticleBasicProperties - One or more input pointers are NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "InitializeParticleBasicProperties - Null input detected.");
    }

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    swarm = user->swarm;

    ierr = DMGetCoordinatesLocal(user->da,&Coor);
    ierr = DMDAVecGetArrayRead(user->fda,Coor,&coor);
    
    LOG_ALLOW(LOCAL, LOG_INFO, "InitializeParticleBasicProperties - Initializing %d particles on rank %d.\n", 
              particlesPerProcess, rank);

    /*--- Retrieve the swarm fields for position, particle IDs, and cell IDs ---*/
    LOG_ALLOW(LOCAL, LOG_DEBUG, "InitializeParticleBasicProperties - Retrieving swarm fields.\n");
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&particleIDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "pos_phy",NULL,NULL,(void**)&pos_phys); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_DEBUG, "InitializeParticleBasicProperties - Particle Initialization = %d\n",user->ParticleInitialization);

    // Initialize random number generators
    ierr = InitializeRandomGenerators(user, randx, randy, randz); CHKERRQ(ierr);

    /// EXPERIMENTAL CODE ! LOGICAL to PHYSICAL SPACE///
    /* choose cell indices owned by this rank */
    PetscInt is,js,ks,nx,ny,nz;
    DMDAGetCorners(user->da,&is,&js,&ks,&nx,&ny,&nz);    
    /*--- Loop over each particle and initialize the basic properties ---*/
    for (PetscInt p = 0; p < particlesPerProcess; p++) {
      PetscReal r;                       /* temporary */
      PetscInt ci,cj,ck;
      PetscRandomGetValueReal(*randx,&r);   /* correct signature */
      ci = is + (PetscInt)((nx-3) * r);

      PetscRandomGetValueReal(*randy,&r);
      cj = js + (PetscInt)((ny-3) * r);

      PetscRandomGetValueReal(*randz,&r);
      ck = ks + (PetscInt)((nz-3) * r);

      /* draw logical coordinates */
      PetscReal xi,eta,zta;
      PetscRandomGetValueReal(*randx,&xi);
      PetscRandomGetValueReal(*randy,&eta);
      PetscRandomGetValueReal(*randz,&zta);
        if (user->ParticleInitialization == 1) {
	  /* logical storage (unchanged) */
	  positions[3*p+0] = xi;
	  positions[3*p+1] = eta;
	  positions[3*p+2] = zta;	   
	  /*
            ierr = PetscRandomGetValue(*randx, &positions[3 * p + 0]); CHKERRQ(ierr);
            ierr = PetscRandomGetValue(*randy, &positions[3 * p + 1]); CHKERRQ(ierr);
            ierr = PetscRandomGetValue(*randz, &positions[3 * p + 2]); CHKERRQ(ierr);
	  */
        } else if (user->ParticleInitialization == 0){
	  /*
	    ierr = PetscRandomGetValue(*randx, &positions[3 * p + 0]); CHKERRQ(ierr);
            ierr = PetscRandomGetValue(*randy, &positions[3 * p + 1]); CHKERRQ(ierr);
	  */
	  positions[3*p+0] = xi;
	  positions[3*p+1] = eta;
	  positions[3 * p + 2] = user->zMin;
	}

            LOG_LOOP_ALLOW(LOCAL,LOG_DEBUG,p,10,"InitializeParticleBasicProperties - Particle %d initialized at (%.6f, %.6f, %.6f).\n",
                      p, positions[3 * p + 0], positions[3 * p + 1], positions[3 * p + 2]);
        /* Assign unique particle IDs */
        particleIDs[p] = rank * particlesPerProcess + p;
        LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG,p,10, "InitializeParticleBasicProperties - Assigned particle ID %lld to particle %d.\n",
                  (long long)particleIDs[p], p);

        /* Initialize cell IDs to -1 */
        cellIDs[3 * p + 0] = -1;
        cellIDs[3 * p + 1] = -1;
        cellIDs[3 * p + 2] = -1;

	/* NEW : physical mapping */
	Cmpnts phys;
	ierr = MetricLogicalToPhysical(user,coor,ci,cj,ck,xi,eta,zta,&phys); CHKERRQ(ierr);
	pos_phys[3*p+0] = phys.x;
	pos_phys[3*p+1] = phys.y;
	pos_phys[3*p+2] = phys.z;
    }

    /*--- Restore the swarm fields ---*/
    LOG_ALLOW(LOCAL, LOG_DEBUG, "InitializeParticleBasicProperties - Restoring swarm fields.\n");
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&particleIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm,"pos_phy",NULL,NULL,(void**)&pos_phys); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda,Coor,&coor);CHKERRQ(ierr);


    LOG_ALLOW(LOCAL, LOG_INFO, "InitializeParticleBasicProperties - Successfully initialized %d particles on rank %d.\n",
              particlesPerProcess, rank);

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
 * @brief Initializes all particle properties in the swarm, using a field-agnostic pipeline.
 *
 * This function first initializes the basic particle properties (position, particle ID, cell ID)
 * by calling InitializeParticleBasicProperties. Then it calls the field-agnostic
 * initializer for additional fields such as "velocity", "weight", and "temperature".
 *
 * @param[in,out] user               Pointer to the UserCtx structure containing simulation context.
 * @param[in]     particlesPerProcess Number of particles assigned to this MPI process.
 * @param[in]     randx              Random number generator for the x-coordinate (if ParticleInitialization==1).
 * @param[in]     randy              Random number generator for the y-coordinate.
 * @param[in]     randz              Random number generator for the z-coordinate.
 * @param[in]     bboxlist           Array of BoundingBox structures, one per MPI rank.
 *
 * @return PetscErrorCode            Returns 0 on success, non-zero on failure.
 */
PetscErrorCode AssignInitialPropertiesToSwarm(UserCtx* user,
                                              PetscInt particlesPerProcess,
                                              PetscRandom *randx, PetscRandom *randy, PetscRandom *randz,
                                              BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    // Validate input pointers
    if (!user || !randx || !randy || !randz || !bboxlist) {
        LOG_ALLOW(GLOBAL, LOG_ERROR, "AssignInitialPropertiesToSwarm - One or more input pointers are NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "AssignInitialPropertiesToSwarm - Null input detected.");
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "AssignInitialPropertiesToSwarm - Initializing swarm with %d particles per process.\n",
              particlesPerProcess);

    /*--- Step 1: Initialize basic particle properties (position, PID, cell IDs) ---*/
    LOG_ALLOW(LOCAL, LOG_DEBUG, "AssignInitialPropertiesToSwarm - Initializing basic particle properties.\n");
    ierr = InitializeParticleBasicProperties(user, particlesPerProcess, randx, randy, randz, bboxlist);
    CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "AssignInitialPropertiesToSwarm - Successfully initialized basic particle properties.\n");

    /*--- Step 2: Initialize fields using the generic pipeline ---*/
    /* Initialize velocity (3 components) */
    LOG_ALLOW(LOCAL, LOG_DEBUG, "AssignInitialPropertiesToSwarm - Initializing velocity field.\n");
    ierr = AssignInitialFieldToSwarm(user, "velocity", 3);
    CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "AssignInitialPropertiesToSwarm - Velocity field initialization complete.\n");

    /* Initialize interpolation weights (3 components) */
    LOG_ALLOW(LOCAL, LOG_DEBUG, "AssignInitialPropertiesToSwarm - Initializing weight field.\n");
    ierr = AssignInitialFieldToSwarm(user, "weight", 3);
    CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "AssignInitialPropertiesToSwarm - Weight field initialization complete.\n");

    /* Initialize interpolation Pressure (1 component) */
    LOG_ALLOW(LOCAL, LOG_DEBUG, "AssignInitialPropertiesToSwarm - Initializing Pressure field.\n");
    ierr = AssignInitialFieldToSwarm(user, "P", 1);
    CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "AssignInitialPropertiesToSwarm - Pressure field initialization complete.\n");
    
    // Future extensions can be added here (ensure corresponding fields are registered)

    LOG_ALLOW(GLOBAL, LOG_INFO, "AssignInitialPropertiesToSwarm - Successfully completed swarm initialization.\n");

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
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode FinalizeSwarmSetup(PetscRandom *randx, PetscRandom *randy, PetscRandom *randz) {
    PetscErrorCode ierr;  // Error code for PETSc functions
    PetscInt  ParticleInitialization; 

    ierr = PetscOptionsGetInt(NULL, NULL, "-pinit", &ParticleInitialization, NULL); CHKERRQ(ierr);
 
    if(ParticleInitialization==1){

      // Destroy random number generators to free resources
      ierr = PetscRandomDestroy(randx); CHKERRQ(ierr);
      ierr = PetscRandomDestroy(randy); CHKERRQ(ierr);
      ierr = PetscRandomDestroy(randz); CHKERRQ(ierr);
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
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting particle swarm Initialization with %d particles.\n", np);

    // Step 1: Create and initialize the particle swarm
    // Here we pass in the bboxlist, which will be used by CreateParticleSwarm() and subsequently
    // by AssignInitialProperties() to position particles if ParticleInitialization == 0.
    LOG_ALLOW(GLOBAL, LOG_INFO, "Creating particle swarm.\n");
    ierr = CreateParticleSwarm(user, np, &particlesPerProcess,bboxlist); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Particle swarm created successfully.\n");

    // Assign initial properties to particles
    // The bboxlist array is passed here so that if ParticleInitialization == 0,
    // particles can be placed at the midpoint of the local bounding box corresponding to this rank.
    ierr = AssignInitialPropertiesToSwarm(user, particlesPerProcess, &randx, &randy, &randz, bboxlist); CHKERRQ(ierr);
    // Finalize swarm setup by destroying RNGs if ParticleInitialization == 1
    ierr = FinalizeSwarmSetup(&randx, &randy, &randz); CHKERRQ(ierr);

    // Ensure all ranks complete before proceeding
    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, " Particles generated & initialized.\n");

    return 0;
}
