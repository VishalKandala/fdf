/**
 * @file setup.c  //  Setup code for running any simulation 
 * @brief Test program for DMSwarm interpolation using the fdf-curvIB method.
 * Provides the setup to start any simulation with DMSwarm and DMDAs.
 **/

 #include "setup.h"

 /**
    * @brief Register events for logging.
    * 
    * Registers events for logging purposes using PetscLogEventRegister.
    *  
 */
 PetscErrorCode registerEvents(void) {
    PetscErrorCode ierr;
    // Register the event with a descriptive name
    ierr = PetscLogEventRegister("walkingsearch", PETSC_OBJECT_CLASSID, &EVENT_walkingsearch);
    ierr = PetscLogEventRegister("Individualwalkingsearch", PETSC_OBJECT_CLASSID, &EVENT_Individualwalkingsearch);
    CHKERRQ(ierr);
    return 0;
}

/**
 * @brief Initialize the simulation context, read options, and parse boundary conditions.
 *
 * - Initializes PETSc and MPI rank/size.
 * - Allocates the UserCtx structure.
 * - Reads runtime options from "control.dat" and command line for simulation parameters
 *   (number of particles, timesteps, initialization modes, etc.).
 * - Calls `ParseAllBoundaryConditions` to read "bcs.dat", determine the BCType for all
 *   6 global faces, and identify the first INLET face for particle initialization (if Mode 0).
 *   This information is stored in the UserCtx.
 * - Initializes PETSc's logging.
 *
 * @param[out] user_ptr    Pointer to the UserCtx structure pointer (will be allocated here).
 * @param[out] rank_out    MPI rank of the process.
 * @param[out] size_out    Number of MPI processes.
 * @param[out] np_out      Number of particles.
 * @param[out] StartStep_out Simulation Starting Timestep.
 * @param[out] StepsToRun_out Number of Timesteps to run simulation.
 * @param[out] StartTime_out Time of start of simulation.
 * @param[out] nblk_out    Number of grid blocks (typically 1 for this setup).
 * @param[out] outputFreq_out Frequency at which data should be output.
 * @param[out] readFields_out Flag to decide if Eulerian fields are read or generated.
 * @param[out] allowedFuncs_out List of functions allowed to show log output.
 * @param[out] nAllowed_out Number of functions allowed to show log output.
 * @param[in,out] allowedFile_inout Path to the config file for allowed log functions.
 * @param[out] useCfg_out Flag indicating if a config file for logging was used.
 * @param[out] OnlySetup Flag indicating if only setup should be run (for debugging) without advancing the simulation in time.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode InitializeSimulation(UserCtx **user_ptr, /* Changed name for clarity */
                                    PetscMPIInt *rank_out, PetscMPIInt *size_out,
                                    PetscInt *np_out, PetscInt *StartStep_out,
                                    PetscInt *StepsToRun_out, PetscReal *StartTime_out,
                                    PetscInt *nblk_out, PetscInt *outputFreq_out,
                                    PetscBool *readFields_out,
                                    char ***allowedFuncs_out, PetscInt *nAllowed_out,
                                    char *allowedFile_inout, PetscBool *useCfg_out,PetscBool* OnlySetup)
{
    PetscErrorCode ierr;
    UserCtx*       local_user; // Temporary local pointer to the user context
    PetscBool      bcsFileOptionFound;
    char           bcs_filename_buffer[PETSC_MAX_PATH_LEN]; // Local buffer for the filename    

    PetscFunctionBeginUser;

    // Get MPI rank and size
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, rank_out); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, size_out); CHKERRQ(ierr);

    // --- 0. PETSc Options and MPI Setup ---
    // Attempt to insert options from "control.dat"
    // This should ideally be done after PetscInitialize, but before most other PETSc calls
    // if options in control.dat are meant to influence PETSc's own setup.
    // However, for user-specific options, it's fine here.
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, "control.dat", PETSC_TRUE);
    if (ierr == PETSC_ERR_FILE_OPEN) {
        // This is a warning if control.dat is optional, an error if it's mandatory.
        LOG_ALLOW(GLOBAL, LOG_WARNING, "InitializeSimulation - Could not open 'control.dat'. Proceeding with command-line options and defaults.\n");
        // SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "InitializeSimulation - Could not open 'control.dat'. Please ensure it exists.");
    } else {
        CHKERRQ(ierr); // Handle other errors from PetscOptionsInsertFile
    }
    
    // --- 1. Logging Configuration First (as it might affect subsequent PETSc/MPI calls) ---
    ierr = PetscOptionsGetString(NULL,NULL,
				 "-func_config_file",
				 allowedFile_inout, // Use the passed-in buffer
				 PETSC_MAX_PATH_LEN,
				 useCfg_out);CHKERRQ(ierr);

    if (*useCfg_out) {
      // Attempt to load from the specified file
      ierr = LoadAllowedFunctionsFromFile(allowedFile_inout, allowedFuncs_out, nAllowed_out);
      if (ierr) { // If loading fails (e.g., file not found, parse error)
	LOG_ALLOW(GLOBAL, LOG_WARNING, "InitializeSimulation - Failed to load allowed functions from '%s' (Error %d). Using default allowed functions (main, AdvanceSimulation).\n", allowedFile_inout, ierr);
	*useCfg_out = PETSC_FALSE; // Mark as not using the config file successfully
	ierr = 0; // Clear the error from LoadAllowedFunctionsFromFile if we are handling it by defaulting
      }
    }

    if (!(*useCfg_out)) { // If config file was not specified, or if loading it failed
      // Default to allowing only "main" and "AdvanceSimulation"
      // Free any previously allocated memory for allowedFuncs_out if LoadAllowedFunctionsFromFile failed partway
      if (*allowedFuncs_out) {
	for (PetscInt i = 0; i < *nAllowed_out; ++i) {
	  ierr = PetscFree((*allowedFuncs_out)[i]); CHKERRQ(ierr);
	}
	ierr = PetscFree(*allowedFuncs_out); CHKERRQ(ierr);
	*allowedFuncs_out = NULL;
      }

      *nAllowed_out = 2; // We are defining 2 default functions
      ierr = PetscMalloc1(*nAllowed_out, allowedFuncs_out); CHKERRQ(ierr);
      ierr = PetscStrallocpy("main", &((*allowedFuncs_out)[0])); CHKERRQ(ierr);
      ierr = PetscStrallocpy("AdvanceSimulation", &((*allowedFuncs_out)[1])); CHKERRQ(ierr);
      LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeSimulation - Using default allowed log functions: main, AdvanceSimulation.\n");
    }
    
    // Register the determined set of allowed functions with the logger
    set_allowed_functions((const char **)*allowedFuncs_out, (size_t)*nAllowed_out);
    
    print_log_level(); // Print current log level early

    // Allocate user context
    ierr = PetscCalloc1(1, &local_user); CHKERRQ(ierr);
    *user_ptr = local_user; // Assign to output parameter

    // Initialize user context flags (good practice to initialize all members)
    local_user->averaging = PETSC_FALSE;
    local_user->les = PETSC_FALSE;
    local_user->rans = PETSC_FALSE;
    local_user->inletFaceDefined = PETSC_FALSE; // Initialize new BC members
    local_user->identifiedInletBCFace = BC_FACE_NEG_X; // Default, will be overwritten
    for (int i = 0; i < 6; ++i) {
        local_user->face_bc_types[i] = WALL; // Default all faces to WALL
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeSimulation - Initialized on rank %d out of %d processes.\n", *rank_out, *size_out);

   // Set default BCs filename
    strcpy(bcs_filename_buffer, "bcs.dat");
    // Attempt to override with command-line/control.dat option
    ierr = PetscOptionsGetString(NULL, NULL, "-bcs_file", bcs_filename_buffer, PETSC_MAX_PATH_LEN, &bcsFileOptionFound); CHKERRQ(ierr);
    if (bcsFileOptionFound) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeSimulation - Using BCs file specified by option: %s\n", bcs_filename_buffer);
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeSimulation - No -bcs_file option found, using default: %s\n", bcs_filename_buffer);
    }
    
    // --- 2. Read Runtime Options for Simulation Parameters ---
    ierr = PetscOptionsGetInt(NULL, NULL, "-numParticles", np_out, NULL); CHKERRQ(ierr);
    local_user->NumberofParticles = *np_out;
    ierr = PetscOptionsGetInt(NULL, NULL, "-rstart", StartStep_out, NULL); CHKERRQ(ierr);
    if((*StartStep_out) > 0) { // Only get StartTime if restarting
        ierr = PetscOptionsGetReal(NULL, NULL, "-ti", StartTime_out, NULL); CHKERRQ(ierr);
    } else {
        *StartTime_out = 0.0; // Ensure StartTime is 0 if not restarting
    }
    ierr = PetscOptionsGetInt(NULL,NULL, "-totalsteps", StepsToRun_out, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-nblk", nblk_out, NULL); CHKERRQ(ierr); // Typically nblk=1
    ierr = PetscOptionsGetInt(NULL, NULL, "-pinit", &(local_user->ParticleInitialization), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-finit", &(local_user->FieldInitialization), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-tio", outputFreq_out, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-logfreq", &(local_user->LoggingFrequency), NULL); CHKERRQ(ierr); // Ensure LoggingFrequency is in UserCtx

    ierr = PetscOptionsGetReal(NULL, NULL, "-dt", &(local_user->dt), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-uin", &(local_user->ConstantVelocity), NULL); CHKERRQ(ierr);

    ierr = PetscOptionsGetBool(NULL, NULL, "-read_fields", readFields_out, NULL); CHKERRQ(ierr);

    ierr = PetscOptionsGetBool(NULL, NULL, "-Setup_Only", OnlySetup, NULL); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL,LOG_DEBUG," -- Console Output Functions [Total : %d] : --\n",*nAllowed_out);
    for (PetscInt i = 0; i < *nAllowed_out; ++i) {
      LOG_ALLOW(GLOBAL,LOG_DEBUG,"   [%2d] «%s»\n", i, (*allowedFuncs_out)[i]);
    }

    // --- 3. Parse Boundary Conditions File ---
    // This is called by all ranks. Rank 0 reads, then broadcasts to all.
    // It populates user->face_bc_types[], user->inletFaceDefined, and user->identifiedInletBCFace.
    LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeSimulation - Parsing boundary conditions from %s.\n", bcs_filename_buffer);
    ierr = ParseAllBoundaryConditions(local_user,bcs_filename_buffer); CHKERRQ(ierr);
    // After this call, user->identifiedInletBCFace (and other BC info) is set on all ranks.

    // --- 4. PETSc Logging Setup ---
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);
    ierr = registerEvents(); CHKERRQ(ierr); // Assuming registerEvents is defined elsewhere

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "InitializeSimulation - Completed successfully.\n");
    PetscFunctionReturn(0);
}

/** 
 * @brief Setup grid and vectors for the simulation.
 *
 * @param[in,out] user Pointer to the UserCtx structure.
 * @param[in] block_number The number of grid blocks in the domain.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode SetupGridAndVectors(UserCtx *user, PetscInt block_number) {
    PetscErrorCode ierr;

    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx pointer is null.");
    if (block_number <= 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "block_number must be > 0");

    // Define grid coordinates
    ierr = DefineGridCoordinates(user); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "SetupGridAndVectors - Grid setup complete.\n");

    ierr = DMDAGetLocalInfo(user->da, &user->info); CHKERRQ(ierr);
   
    // Create global vectors for each block
      for (PetscInt bi = 0; bi < block_number; bi++) {
        // Check that the DMs are valid
        if (!user[bi].fda || !user[bi].da) {
           SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "DMs not initialized \n");
     }

        LOG_ALLOW(GLOBAL,LOG_DEBUG,"Creating vectors for block %d: fda=%p, da=%p\n", bi, (void*)user[bi].fda, (void*)user[bi].da);

        // Create global vectors (Destroyed in FinalizeSimulation)
	ierr = DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user[bi].fda, &user[bi].Ucont); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user[bi].da, &user[bi].P); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user[bi].da, &user[bi].Nvert); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user[bi].da, &user[bi].Nvert_o); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(user[bi].da, &user[bi].ParticleCount); CHKERRQ(ierr);

	// Create local vectors (Destroyed in FinalizeSimulation)
        ierr = DMCreateLocalVector(user[bi].fda, &user[bi].lUcat); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(user[bi].fda, &user[bi].lUcont); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(user[bi].da, &user[bi].lP); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(user[bi].da, &user[bi].lNvert); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(user[bi].da, &user[bi].lNvert_o); CHKERRQ(ierr);

	// Set Initial values for all vectors
	ierr = VecSet(user[bi].Ucat,0.0);
	ierr = VecSet(user[bi].lUcat,0.0);
	
	ierr = VecSet(user[bi].Ucont,0.0);	
	ierr = VecSet(user[bi].lUcont,0.0);

	ierr = VecSet(user[bi].P,0.0);
	ierr = VecSet(user[bi].lP,0.0);

	ierr = VecSet(user[bi].Nvert,0.0);
	ierr = VecSet(user[bi].lNvert,0.0);

	ierr = VecSet(user[bi].Nvert_o,0.0);
	ierr = VecSet(user[bi].lNvert_o,0.0);

	ierr = VecSet(user[bi].ParticleCount,0.0);
      } //bi 

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "SetupGridAndVectors - Grid and vectors setup completed on all ranks. \n");
    return 0;
}

/**
 * @brief Finalize the simulation and free resources.
 *
 * @param[in,out] user Pointer to the UserCtx structure.
 * @param[in] block_number The number of grid blocks in the domain.
 * @param[in,out] bboxlist  Pointer to the array of BoundingBoxes.
 * @param[in,out] allowedFuncs list of functions that are allowed to show output.
 * @param[in] nAllowed No.of functions allowed to show output.
 * @param[in] PetSc viewer data type for logging 
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode FinalizeSimulation(UserCtx *user, PetscInt block_number, BoundingBox *bboxlist, char **allowedFuncs, PetscInt nAllowed,PetscViewer *logviewer) {
    PetscErrorCode ierr;
    PetscMPIInt rank;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr); // Corrected: use &rank


    // Create an ASCII viewer to write log output to file
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "simulationlog.txt", logviewer);

    LOG_ALLOW(GLOBAL,LOG_INFO," PETSC Logs written \n");
    
    // Print PETSc logging results at the end
    ierr = PetscLogView(*logviewer); CHKERRQ(ierr);
    
    // Destroy DM and vectors for each block
    for (PetscInt bi = 0; bi < block_number; bi++) {
        ierr = VecDestroy(&(user[bi].Ucat)); CHKERRQ(ierr);
	ierr = VecDestroy(&(user[bi].lUcat)); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," Ucat Destroyed \n");
        ierr = VecDestroy(&(user[bi].Ucont)); CHKERRQ(ierr);
        ierr = VecDestroy(&(user[bi].lUcont)); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," Ucont Destroyed \n");
        ierr = VecDestroy(&(user[bi].P)); CHKERRQ(ierr);
	ierr = VecDestroy(&(user[bi].lP)); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," P Destroyed \n");
        ierr = VecDestroy(&(user[bi].Nvert)); CHKERRQ(ierr);
	ierr = VecDestroy(&(user[bi].lNvert)); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," Nvert Destroyed \n");
        ierr = VecDestroy(&(user[bi].Nvert_o)); CHKERRQ(ierr);
        ierr = VecDestroy(&(user[bi].lNvert_o)); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," Nvert_o Destroyed \n");
	ierr = VecDestroy(&(user[bi].ParticleCount)); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," ParticleCount Destroyed \n");
	
	//  ierr = DMDestroy(&(user[bi].fda)); CHKERRQ(ierr);
	//	LOG_ALLOW(GLOBAL,LOG_DEBUG," fda Destroyed \n");
	
        ierr = DMDestroy(&(user[bi].da)); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," da Destroyed \n");
        ierr = DMDestroy(&(user[bi].swarm)); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," swarm Destroyed \n");
    }

    // Free user context
    ierr = PetscFree(user); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_DEBUG," user Destroyed \n");     
    // Now free bboxlist on all ranks since all allocated their own copy
    if (bboxlist) {
        free(bboxlist);
        bboxlist = NULL;
    }

    if (allowedFuncs && nAllowed > 0) {
        ierr = FreeAllowedFunctions(allowedFuncs, nAllowed); CHKERRQ(ierr);
    }
    
    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "Simulation finalized and resources cleaned up on all ranks. \n");

    // Ensure all MPI ranks reach this point before finalizing PETSc
    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

    // Finalize PETSc
    ierr = PetscFinalize(); CHKERRQ(ierr);
    return 0;
}

/**
 * @brief Allocates a 3D array of PetscReal values using PetscCalloc.
 *
 * This function dynamically allocates memory for a 3D array of PetscReal values
 * with dimensions nz (layers) x ny (rows) x nx (columns). It uses PetscCalloc1
 * to ensure the memory is zero-initialized.
 *
 * The allocation is done in three steps:
 *  1. Allocate an array of nz pointers (one for each layer).
 *  2. Allocate a contiguous block for nz*ny row pointers and assign each layer’s row pointers.
 *  3. Allocate a contiguous block for all nz*ny*nx PetscReal values.
 *
 * This setup allows the array to be accessed as array[k][j][i], and the memory
 * for the data is contiguous, which improves cache efficiency.
 *
 * @param[out] array Pointer to the 3D array to be allocated.
 * @param[in]  nz    Number of layers (z-direction).
 * @param[in]  ny    Number of rows (y-direction).
 * @param[in]  nx    Number of columns (x-direction).
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
PetscErrorCode Allocate3DArrayScalar(PetscReal ****array, PetscInt nz, PetscInt ny, PetscInt nx)
{
  PetscErrorCode ierr;
  PetscReal      ***data;
  PetscReal      *dataContiguous;
  PetscInt       k, j;

  PetscFunctionBegin;
  /* Step 1: Allocate memory for an array of nz layer pointers (zero-initialized) */
  ierr = PetscCalloc1(nz, &data); CHKERRQ(ierr);

  /* Step 2: Allocate memory for all row pointers (nz * ny pointers) */
  ierr = PetscCalloc1(nz * ny, &data[0]); CHKERRQ(ierr);
  for (k = 1; k < nz; k++) {
    data[k] = data[0] + k * ny;
  }

  /* Step 3: Allocate one contiguous block for all data elements (nz*ny*nx) */
  ierr = PetscCalloc1(nz * ny * nx, &dataContiguous); CHKERRQ(ierr);

  /* Build the 3D pointer structure: each row pointer gets the correct segment of data */
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      data[k][j] = dataContiguous + (k * ny + j) * nx;
      /* Memory is already zeroed by PetscCalloc1, so no manual initialization is needed */
    }
  }
  *array = data;
  PetscFunctionReturn(0);
}

/**
 * @brief Deallocates a 3D array of PetscReal values allocated by Allocate3DArrayScalar.
 *
 * This function frees the memory allocated for a 3D array of PetscReal values.
 * It assumes the memory was allocated using Allocate3DArrayScalar, which allocated
 * three separate memory blocks: one for the contiguous data, one for the row pointers,
 * and one for the layer pointers.
 *
 * @param[in] array Pointer to the 3D array to be deallocated.
 * @param[in] nz    Number of layers (z-direction).
 * @param[in] ny    Number of rows (y-direction).
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
PetscErrorCode Deallocate3DArrayScalar(PetscReal ***array, PetscInt nz, PetscInt ny)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
   if (!array || !array[0] || !array[0][0] ) { // Added more robust check
      LOG_ALLOW(GLOBAL, LOG_WARNING, "Deallocate3DArrayScalar called with potentially unallocated or NULL array.\n");
       if (array) {
           if (array[0]) { // Check if row pointers might exist
               // Cannot safely access array[0][0] if array[0] might be invalid/freed
               // Standard deallocation below assumes valid pointers.
                ierr = PetscFree(array[0]); CHKERRQ(ierr); // Free row pointers if they exist
           }
            ierr = PetscFree(array); CHKERRQ(ierr); // Free layer pointers if they exist
       }
       PetscFunctionReturn(0);
   }

  // --- Standard Deallocation (assuming valid allocation) ---

  /* 1. Free the contiguous block of PetscReal values.
     The starting address was stored in array[0][0]. */
  ierr = PetscFree(array[0][0]); CHKERRQ(ierr); // Free the ACTUAL DATA

  /* 2. Free the contiguous block of row pointers.
     The starting address was stored in array[0]. */
  ierr = PetscFree(array[0]); CHKERRQ(ierr); // Free the ROW POINTERS

  /* 3. Free the layer pointer array.
     The starting address is 'array' itself. */
  ierr = PetscFree(array); CHKERRQ(ierr); // Free the LAYER POINTERS

  PetscFunctionReturn(0);
}

/**
 * @brief Allocates a 3D array of Cmpnts structures using PetscCalloc.
 *
 * This function dynamically allocates memory for a 3D array of Cmpnts (vector) structures
 * with dimensions nz (layers) x ny (rows) x nx (columns). It uses PetscCalloc1 to ensure
 * that all allocated memory is zero-initialized.
 *
 * The allocation procedure is similar to Allocate3DArrayScalar:
 *  1. Allocate an array of nz pointers (one for each layer).
 *  2. Allocate a contiguous block for nz*ny row pointers.
 *  3. Allocate one contiguous block for nz*ny*nx Cmpnts structures.
 *
 * After allocation, the array can be accessed as array[k][j][i] and each element
 * (a Cmpnts structure) will have its x, y, and z fields initialized to 0.0.
 *
 * @param[out] array Pointer to the 3D array to be allocated.
 * @param[in]  nz    Number of layers in the z-direction.
 * @param[in]  ny    Number of rows in the y-direction.
 * @param[in]  nx    Number of columns in the x-direction.
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
PetscErrorCode Allocate3DArrayVector(Cmpnts ****array, PetscInt nz, PetscInt ny, PetscInt nx)
{
  PetscErrorCode ierr;
  Cmpnts         ***data;
  Cmpnts         *dataContiguous;
  PetscInt       k, j;

  PetscFunctionBegin;
  /* Step 1: Allocate memory for nz layer pointers (zeroed) */
  ierr = PetscCalloc1(nz, &data); CHKERRQ(ierr);

  /* Step 2: Allocate memory for all row pointers (nz * ny pointers) */
  ierr = PetscCalloc1(nz * ny, &data[0]); CHKERRQ(ierr);
  for (k = 1; k < nz; k++) {
    data[k] = data[0] + k * ny;
  }

  /* Step 3: Allocate one contiguous block for nz*ny*nx Cmpnts structures (zeroed) */
  ierr = PetscCalloc1(nz * ny * nx, &dataContiguous); CHKERRQ(ierr);

  /* Build the 3D pointer structure for vector data */
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      data[k][j] = dataContiguous + (k * ny + j) * nx;
      /* The PetscCalloc1 call has already initialized each Cmpnts to zero. */
    }
  }
  *array = data;
  PetscFunctionReturn(0);
}

/**
 * @brief Deallocates a 3D array of Cmpnts structures allocated by Allocate3DArrayVector.
 *
 * This function frees the memory allocated for a 3D array of Cmpnts structures.
 * It assumes the memory was allocated using Allocate3DArrayVector, which created three
 * separate memory blocks: one for the contiguous vector data, one for the row pointers,
 * and one for the layer pointers.
 *
 * @param[in] array Pointer to the 3D array to be deallocated.
 * @param[in] nz    Number of layers in the z-direction.
 * @param[in] ny    Number of rows in the y-direction.
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
 PetscErrorCode Deallocate3DArrayVector(Cmpnts ***array, PetscInt nz, PetscInt ny)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  // If array is NULL or hasn't been allocated properly, just return.
  if (!array || !array[0] || !array[0][0] ) {
      LOG_ALLOW(GLOBAL, LOG_WARNING, "Deallocate3DArrayVector called with potentially unallocated or NULL array.\n");
      // Attempt to free what might exist, but be cautious
      if (array) {
          if (array[0]) { // Check if row pointers were allocated
             // We don't have a direct pointer to the contiguous data block
             // saved separately in this allocation scheme. The allocation relies
             // on array[0][0] pointing to it. If array[0] was freed first,
             // accessing array[0][0] is unsafe.
             // The allocation scheme where the contiguous data block is not
             // stored separately makes safe deallocation tricky if freeing
             // happens out of order or if parts are NULL.

             // A SAFER ALLOCATION/DEALLOCATION would store the data pointer separately.
             // Given the current allocation scheme, the order MUST be:
             // 1. Free the data block (pointed to by array[0][0])
             // 2. Free the row pointer block (pointed to by array[0])
             // 3. Free the layer pointer block (pointed to by array)

             // Let's assume the allocation was successful and pointers are valid.
             // Get pointer to the contiguous data block *before* freeing row pointers
             Cmpnts *dataContiguous = array[0][0];
             ierr = PetscFree(dataContiguous); CHKERRQ(ierr); // Free data block

             // Now free the row pointers block
             ierr = PetscFree(array[0]); CHKERRQ(ierr); // Free row pointers

          }
          // Finally, free the array of layer pointers
          ierr = PetscFree(array); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0); // Return gracefully if input was NULL initially
  }


  // --- Standard Deallocation (assuming valid allocation) ---

  /* 1. Free the contiguous block of Cmpnts structures.
     The starting address was stored in array[0][0] by Allocate3DArrayVector. */
  ierr = PetscFree(array[0][0]); CHKERRQ(ierr); // Free the ACTUAL DATA

  /* 2. Free the contiguous block of row pointers.
     The starting address was stored in array[0]. */
  ierr = PetscFree(array[0]); CHKERRQ(ierr); // Free the ROW POINTERS

  /* 3. Free the layer pointer array.
     The starting address is 'array' itself. */
  ierr = PetscFree(array); CHKERRQ(ierr); // Free the LAYER POINTERS

  PetscFunctionReturn(0);
}

/**
 * @brief Gets the global starting index of cells owned by this rank and the number of such cells.
 *
 * A cell's global index is considered the same as its origin node's global index.
 * This function assumes a node-centered DMDA where `info_nodes` provides all necessary
 * information:
 *  - `info_nodes->xs, ys, zs`: Global starting index of the first node owned by this rank (excluding ghosts).
 *  - `info_nodes->xm, ym, zm`: Number of nodes owned by this rank in each dimension (excluding ghosts).
 *  - `info_nodes->mx, my, mz`: Total number of global nodes in each dimension for the entire domain.
 *
 * A cell `C_k` (0-indexed) is defined by its origin node `N_k` and extends to node `N_{k+1}`.
 * Thus, the last node in the global domain cannot be an origin for a cell. The last possible
 * cell origin node index is `GlobalNodesInDim - 2`.
 *
 * @param[in] info_nodes Pointer to the DMDALocalInfo struct for the current rank.
 *                       This struct contains local ownership information (xs, xm, etc.)
 *                       and global domain dimensions (mx, my, mz for nodes).
 * @param[in] dim        The dimension for which to get the cell range (0 for i/x, 1 for j/y, 2 for k/z).
 * @param[out] xs_cell_global_out Pointer to store the global index of the first cell whose origin node
 *                                is owned by this rank. If the rank owns no valid cell origins in this
 *                                dimension, this will be the rank's starting node index, but
 *                                `xm_cell_local_out` will be 0.
 * @param[out] xm_cell_local_out  Pointer to store the number of cells for which this rank owns the
 *                                origin node AND that origin node is a valid cell origin within the
 *                                global domain.
 *
 * @return PetscErrorCode 0 on success, or an error code on failure.
 *
 * @note Example: If GlobalNodesInDim = 11 (nodes N0 to N10), there are 10 cells (C0 to C9).
 *       The last cell, C9, has its origin at node N9. So, N9 (index 9) is the last valid
 *       cell origin (GlobalNodesInDim - 2 = 11 - 2 = 9).
 *       If a rank owns nodes N8, N9, N10 (xs=8, xm=3):
 *         - First potential origin on rank = N8.
 *         - Last potential origin on rank (node that is not the last owned node) = N9.
 *         - Actual last origin this rank can form = min(N9, GlobalMaxOrigin=N9) = N9.
 *         - Number of cells = (N9 - N8 + 1) = 2 cells (C8, C9).
 *       If a rank owns only node N10 (xs=10, xm=1):
 *         - First potential origin on rank = N10.
 *         - Actual last origin rank can form = min(N9, GlobalMaxOrigin=N9) (since N10-1=N9).
 *         - first_potential_origin_on_rank (N10) > actual_last_origin_this_rank_can_form (N9) => 0 cells.
 */
PetscErrorCode GetOwnedCellRange(const DMDALocalInfo *info_nodes,
                                 PetscInt dim,
                                 PetscInt *xs_cell_global_out,
                                 PetscInt *xm_cell_local_out)
{
    PetscErrorCode ierr = 0; // Standard PETSc error code, not explicitly set here but good practice.
    PetscInt xs_node_global_rank;   // Global index of the first node owned by this rank in the specified dimension.
    PetscInt num_nodes_owned_rank;  // Number of nodes owned by this rank in this dimension (local count, excluding ghosts).
    PetscInt GlobalNodesInDim_from_info; // Total number of global nodes in this dimension, from DMDALocalInfo.

    PetscFunctionBeginUser;

    // --- 1. Input Validation ---
    if (!info_nodes || !xs_cell_global_out || !xm_cell_local_out) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null pointer passed to GetOwnedCellRange.");
    }

    // --- 2. Extract Node Ownership and Global Dimension Information from DMDALocalInfo ---
    if (dim == 0) { // I-direction
        xs_node_global_rank = info_nodes->xs;       // Starting owned node index (global)
        num_nodes_owned_rank  = info_nodes->xm;     // Number of nodes owned by this rank (local count)
        GlobalNodesInDim_from_info = info_nodes->mx; // Total global nodes in this dimension
    } else if (dim == 1) { // J-direction
        xs_node_global_rank = info_nodes->ys;
        num_nodes_owned_rank  = info_nodes->ym;
        GlobalNodesInDim_from_info = info_nodes->my;
    } else if (dim == 2) { // K-direction
        xs_node_global_rank = info_nodes->zs;
        num_nodes_owned_rank  = info_nodes->zm;
        GlobalNodesInDim_from_info = info_nodes->mz;
    } else {
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Invalid dimension %d in GetOwnedCellRange. Must be 0, 1, or 2.", dim);
    }

    // --- 3. Handle Edge Cases for Global Domain Size ---
    // If the global domain has 0 or 1 node plane, no cells can be formed.
    if (GlobalNodesInDim_from_info <= 1) {
        *xs_cell_global_out = xs_node_global_rank; // Still report the rank's starting node
        *xm_cell_local_out = 0;                    // But 0 cells
        PetscFunctionReturn(0);
    }
    // Negative global dimension is an error (should be caught by DMDA setup, but defensive)
    if (GlobalNodesInDim_from_info < 0 ) {
         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "GlobalNodesInDim %d from DMDALocalInfo must be non-negative for dimension %d.", GlobalNodesInDim_from_info, dim);
    }

    // --- 4. Determine Cell Ownership Based on Node Ownership ---
    // The first cell this rank *could* define has its origin at the first node this rank owns.
    *xs_cell_global_out = xs_node_global_rank;

    // If the rank owns no nodes in this dimension, it can't form any cell origins.
    if (num_nodes_owned_rank == 0) {
        *xm_cell_local_out = 0;
    } else {
        // Calculate the global index of the last possible node that can serve as a cell origin.
        // If GlobalNodesInDim = N (nodes 0 to N-1), cells are C_0 to C_{N-2}.
        // The origin of cell C_{N-2} is node N_{N-2}.
        // So, the last valid cell origin node index is (GlobalNodesInDim - 2).
        PetscInt last_possible_origin_global_idx = GlobalNodesInDim_from_info - 2;

        // Determine the range of nodes owned by this rank that could *potentially* be cell origins.
        // The first node owned by the rank is a potential origin.
        PetscInt first_potential_origin_on_rank = xs_node_global_rank;

        // A node can be an origin if there's at least one node after it to form the cell.
        // So, the last node owned by the rank that could *potentially* be an origin is
        // the second-to-last node it owns: (xs_node_global_rank + num_nodes_owned_rank - 1) - 1
        // which simplifies to: xs_node_global_rank + num_nodes_owned_rank - 2.
        PetscInt last_potential_origin_on_rank = xs_node_global_rank + num_nodes_owned_rank - 2;

        // The actual last origin this rank can provide is capped by the global domain limit.
        PetscInt actual_last_origin_this_rank_can_form = PetscMin(last_potential_origin_on_rank, last_possible_origin_global_idx);

        // If the first potential origin this rank owns is already beyond the actual last origin it can form,
        // then this rank forms no valid cell origins. This happens if:
        //  - num_nodes_owned_rank is 1 (so last_potential_origin_on_rank = first_potential_origin_on_rank - 1).
        //  - The rank only owns nodes at the very end of the global domain (e.g., only the last global node).
        if (first_potential_origin_on_rank > actual_last_origin_this_rank_can_form) {
            *xm_cell_local_out = 0;
        } else {
            // The number of cells is the count of valid origins this rank owns.
            // (Count = Last Index - First Index + 1)
            *xm_cell_local_out = actual_last_origin_this_rank_can_form - first_potential_origin_on_rank + 1;
        }
    }
    PetscFunctionReturn(ierr);
}

/**
 * @brief Updates the local vector (including ghost points) from its corresponding global vector.
 *
 * This function identifies the correct global vector, local vector, and DM based on the
 * provided fieldName and performs the standard PETSc DMGlobalToLocalBegin/End sequence.
 * Includes optional debugging output (max norms before/after).
 *
 * @param user       The UserCtx structure containing the vectors and DMs.
 * @param fieldName  The name of the field to update ("Ucat", "Ucont", "P", "Nvert", etc.).
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 *
 * @note This function assumes the global vector associated with fieldName has already
 *       been populated with the desired data (including any boundary conditions).
 */
PetscErrorCode UpdateLocalGhosts(UserCtx* user, const char *fieldName)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    Vec            globalVec = NULL;
    Vec            localVec = NULL;
    DM             dm = NULL; // The DM associated with this field pair

    PetscFunctionBeginUser; // Use User version for application code
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Starting ghost update for field '%s'.\n", rank, fieldName);

    // --- 1. Identify the correct Vectors and DM ---
    if (strcmp(fieldName, "Ucat") == 0) {
        globalVec = user->Ucat;
        localVec  = user->lUcat;
        dm        = user->fda;
    } else if (strcmp(fieldName, "Ucont") == 0) {
        globalVec = user->Ucont;
        localVec  = user->lUcont;
        dm        = user->fda;
    } else if (strcmp(fieldName, "P") == 0) {
        globalVec = user->P;
        localVec  = user->lP;
        dm        = user->da;
    } else if (strcmp(fieldName, "Nvert") == 0) {
        globalVec = user->Nvert;
        localVec  = user->lNvert;
        dm        = user->da;
    } // Add other fields as needed
    else {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Field '%s' not recognized for ghost update.", fieldName);
    }

    // --- 2. Check if components were found ---
    if (!globalVec) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Global vector for field '%s' is NULL.", fieldName);
    }
    if (!localVec) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Local vector for field '%s' is NULL.", fieldName);
    }
    if (!dm) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "DM for field '%s' is NULL.", fieldName);
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Identified components for '%s': DM=%p, GlobalVec=%p, LocalVec=%p.\n",
              rank, fieldName, (void*)dm, (void*)globalVec, (void*)localVec);

    // --- 3. Optional Debugging: Norm Before Update ---
    // Use your logging convention check
    // if (get_log_level() >= LOG_LEVEL_DEBUG && is_function_allowed("UpdateLocalGhosts")) { // Example check
    if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)){
        PetscReal norm_global_before;
        ierr = VecNorm(globalVec, NORM_INFINITY, &norm_global_before); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_INFO,"Max norm '%s' (Global) BEFORE Ghost Update: %g\n", fieldName, norm_global_before);
        // Optional: Norm of local vector before update (might contain old ghost values)
        // PetscReal norm_local_before;
        // ierr = VecNorm(localVec, NORM_INFINITY, &norm_local_before); CHKERRQ(ierr);
        // LOG_ALLOW(GLOBAL, LOG_DEBUG,"Max norm '%s' (Local) BEFORE Ghost Update: %g\n", fieldName, norm_local_before);
    }

    // --- 4. Perform the Global-to-Local Transfer (Ghost Update) ---
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Calling DMGlobalToLocalBegin/End for '%s'.\n", rank, fieldName);
    ierr = DMGlobalToLocalBegin(dm, globalVec, INSERT_VALUES, localVec); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dm, globalVec, INSERT_VALUES, localVec); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Completed DMGlobalToLocalBegin/End for '%s'.\n", rank, fieldName);

    // --- 5. Optional Debugging: Norm After Update ---
    // Use your logging convention check
    // if (get_log_level() >= LOG_LEVEL_DEBUG && is_function_allowed("UpdateLocalGhosts")) { // Example check
    if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)){ // Using your specific check
        PetscReal norm_local_after;
        ierr = VecNorm(localVec, NORM_INFINITY, &norm_local_after); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_INFO,"Max norm '%s' (Local) AFTER Ghost Update: %g\n", fieldName, norm_local_after);

        // --- 6. Optional Debugging: Specific Point Checks (Example for Ucat on Rank 0/1) ---
        //    (Keep this conditional if it's only for specific debug scenarios)
        if (strcmp(fieldName, "Ucat") == 0) { // Only do detailed checks for Ucat for now
           PetscMPIInt rank_test;
           MPI_Comm_rank(PETSC_COMM_WORLD, &rank_test);

           // Get Local Info needed for indexing checks
           DMDALocalInfo info_check;
           ierr = DMDAGetLocalInfo(dm, &info_check); CHKERRQ(ierr); // Use the correct dm

           // Buffer for array pointer
           Cmpnts ***lUcat_arr_test = NULL;
           PetscErrorCode ierr_test = 0;

           LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Testing '%s' access immediately after ghost update...\n", rank_test, fieldName);
           ierr_test = DMDAVecGetArrayDOFRead(dm, localVec, &lUcat_arr_test); // Use correct dm and localVec

           if (ierr_test) {
               LOG_ALLOW(LOCAL, LOG_ERROR, "Rank %d: ERROR %d getting '%s' array after ghost update!\n", rank_test, ierr_test, fieldName);
           } else if (!lUcat_arr_test) {
                LOG_ALLOW(LOCAL, LOG_ERROR, "Rank %d: ERROR NULL pointer getting '%s' array after ghost update!\n", rank_test, fieldName);
           }
           else {
               // Check owned interior point (e.g., first interior point)
               PetscInt k_int = info_check.zs + (info_check.zm > 1 ? 1 : 0); // Global k index (at least zs+1 if possible)
               PetscInt j_int = info_check.ys + (info_check.ym > 1 ? 1 : 0); // Global j index
               PetscInt i_int = info_check.xs + (info_check.xm > 1 ? 1 : 0); // Global i index
                // Ensure indices are within global bounds if domain is very small
               if (k_int >= info_check.mz-1) k_int = info_check.mz-2; if (k_int < 1) k_int = 1;
               if (j_int >= info_check.my-1) j_int = info_check.my-2; if (j_int < 1) j_int = 1;
               if (i_int >= info_check.mx-1) i_int = info_check.mx-2; if (i_int < 1) i_int = 1;

               // Only attempt read if indices are actually owned (relevant for multi-rank)
               if (k_int >= info_check.zs && k_int < info_check.zs + info_check.zm &&
                   j_int >= info_check.ys && j_int < info_check.ys + info_check.ym &&
                   i_int >= info_check.xs && i_int < info_check.xs + info_check.xm)
               {
                  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Attempting test read OWNED INTERIOR [%d][%d][%d] (Global)\n", rank_test, k_int, j_int, i_int);
                  Cmpnts test_val_owned_interior = lUcat_arr_test[k_int][j_int][i_int];
                  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: SUCCESS reading owned interior: x=%g\n", rank_test, test_val_owned_interior.x);
               } else {
                  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Skipping interior test read for non-owned index [%d][%d][%d].\n", rank_test, k_int, j_int, i_int);
               }


               // Check owned boundary point (e.g., first owned point)
               PetscInt k_bnd = info_check.zs; // Global k index
               PetscInt j_bnd = info_check.ys; // Global j index
               PetscInt i_bnd = info_check.xs; // Global i index
               LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Attempting test read OWNED BOUNDARY [%d][%d][%d] (Global)\n", rank_test, k_bnd, j_bnd, i_bnd);
               Cmpnts test_val_owned_boundary = lUcat_arr_test[k_bnd][j_bnd][i_bnd];
               LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: SUCCESS reading owned boundary: x=%g\n", rank_test, test_val_owned_boundary.x);


               // Check ghost point (e.g., one layer below in k, if applicable)
               if (info_check.zs > 0) { // Only if there's a rank below
                   PetscInt k_ghost = info_check.zs - 1;
                   PetscInt j_ghost = info_check.ys; // Use start of owned y, simple example
                   PetscInt i_ghost = info_check.xs; // Use start of owned x, simple example
                   LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Attempting test read GHOST [%d][%d][%d] (Global)\n", rank_test, k_ghost, j_ghost, i_ghost);
                   Cmpnts test_val_ghost = lUcat_arr_test[k_ghost][j_ghost][i_ghost];
                   LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: SUCCESS reading ghost: x=%g\n", rank_test, test_val_ghost.x);
               } else {
                   LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Skipping ghost test read (zs=0).\n", rank_test);
               }

               // Restore the array
               ierr_test = DMDAVecRestoreArrayDOFRead(dm, localVec, &lUcat_arr_test);
               if(ierr_test){ LOG_ALLOW(LOCAL, LOG_ERROR, "Rank %d: ERROR %d restoring '%s' array after test read!\n", rank_test, ierr_test, fieldName); }
               LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Finished testing '%s' access.\n", rank_test, fieldName);
           }
        } // end if Ucat
    } // end debug logging check

    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Completed ghost update for field '%s'.\n", rank, fieldName);
    PetscFunctionReturn(0);
}

/**
 * @brief Performs the complete initial setup for the particle simulation at time t=0.
 *
 * This includes:
 * 1. Initial locating of particles (based on their potentially arbitrary initial assignment).
 * 2. A preliminary migration cycle to ensure particles are on the MPI rank that owns
 *    their initial physical region.
 * 3. If `user->ParticleInitialization == 0` (Surface Init), re-initializes particles on the
 *    designated inlet surface. This ensures particles migrated to an inlet-owning rank
 *    are correctly distributed on that surface.
 * 4. A final locating of all particles to get their correct cell indices and interpolation weights.
 * 5. Interpolation of initial Eulerian fields to the particles.
 * 6. Scattering of particle data to Eulerian fields (if applicable).
 * 7. Outputting initial data if requested.
 *
 * @param user Pointer to the UserCtx structure.
 * @param currentTime The current simulation time (should be StartTime, typically 0.0).
 * @param step The current simulation step (should be StartStep, typically 0).
 * @param readFields Flag indicating if Eulerian fields were read from file (influences output).
 * @param bboxlist Array of BoundingBox structures for domain decomposition.
 * @param OutputFreq Frequency for writing output files.
 * @param StepsToRun Total number of simulation steps planned (used for output logic on setup-only runs).
 * @param StartStep The starting step of the simulation (used for output logic).
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformInitialSetup(UserCtx *user, PetscReal currentTime, PetscInt step,
                                   PetscBool readFields, const BoundingBox *bboxlist,
                                   PetscInt OutputFreq, PetscInt StepsToRun, PetscInt StartStep)
{
    PetscErrorCode ierr;
    MigrationInfo  *migrationList_for_initial_sort = NULL; // Managed locally for this initial sort
    PetscInt       migrationCount_initial_sort = 0;
    PetscInt       migrationListCapacity_initial_sort = 0;
    PetscInt       globalMigrationCount_initial_sort;
    PetscMPIInt    rank; // For logging

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Performing initial particle setup procedures.\n", currentTime, step);

    // --- 1. Initial Locate (based on positions from InitializeParticleBasicProperties) ---
    // This gets a first guess of where particles are, which might be (0,0,0) for many if not placed by their initial rank.
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Initial Locate: Determining particle cells before preliminary migration.\n", currentTime, step);
    ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL,LOG_INFO," Particle layout (pre-preliminary migration): \n");
    ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);

    // --- 2. Preliminary Spatial Sorting Migration ---
    // This moves particles (e.g., those at (0,0,0)) to the rank that actually owns their physical region.
    ierr = PerformSingleParticleMigrationCycle(user, bboxlist,
                                               &migrationList_for_initial_sort, &migrationCount_initial_sort, &migrationListCapacity_initial_sort,
                                               currentTime, step, "Preliminary Sort", &globalMigrationCount_initial_sort); CHKERRQ(ierr);
    // migrationList_for_initial_sort is allocated/reallocated within PerformSingleParticleMigrationCycle if needed.

    ierr = PetscBarrier(NULL);
    
    // --- 3. Re-initialize Particles on Inlet Surface (if applicable) ---
    // This is crucial for ParticleInitialization == 0 (Surface Init). After particles have been
    // migrated to the rank(s) owning the inlet surface, this step ensures they are properly
    // distributed *on* that surface, rather than remaining at their migrated position (e.g., (0,0,0)).
    if (user->ParticleInitialization == 0 && user->inletFaceDefined) {
        ierr = ReinitializeParticlesOnInletSurface(user, currentTime, step); CHKERRQ(ierr);
    }

    // --- 4. Final Locate and Interpolate ---
    // Now that particles are on their correct ranks and (if surface init) correctly positioned on the inlet,
    // perform the definitive location and interpolation for t=0.
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Final Initial Locate & Interpolate (post-migration & re-placement).\n", currentTime, step);
    ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);
    ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr);

    // --- 5. Scatter Particle Data to Eulerian Grid (if part of initialization) ---
    ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);

    // --- 6. Initial Output ---
    // Output if OutputFreq > 0 OR if it's a setup-only run (StepsToRun==0 and StartStep==0).
    if (OutputFreq > 0 || (StepsToRun == 0 && StartStep == 0)) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Logging initial data and (if applicable) interpolation error.\n", currentTime, step);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Printing Initial particle fields at t=%.4f (step %d completed)\n", currentTime, step);
        ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);

        if(user->FieldInitialization==1) { // If analytic fields are available for comparison
            ierr = LOG_INTERPOLATION_ERROR(user); CHKERRQ(ierr);
        }

        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Writing initial particle data (for step %d).\n", currentTime, step, step);
        ierr = WriteSwarmField(user, "position", step, "dat"); CHKERRQ(ierr);
        ierr = WriteSwarmField(user, "velocity", step, "dat"); CHKERRQ(ierr);
        // ierr = WriteSwarmField(user, "pos_phy",  step, "dat"); CHKERRQ(ierr); // If needed
        if (!readFields) { // Only write grid fields if they were generated, not read
            ierr = WriteSimulationFields(user); CHKERRQ(ierr);
        }
    }

    // --- 7. Cleanup Memory ---
    // Free the migration list used specifically for this initial sort.
    ierr = PetscFree(migrationList_for_initial_sort); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/**
 * @brief Executes the main time-marching loop for the particle simulation.
 *
 * This function performs the following steps repeatedly for `StepsToRun`:
 * 1. Updates/Sets the background fluid velocity field (Ucat) for the current step.
 * 2. If it's the very first step (StartStep=0), calls `PerformInitialSetup` for
 *    preliminary migration, location, interpolation, and output.
 * 3. Updates particle positions using velocity from the *previous* step's interpolation.
 * 4. Performs boundary checks and removes out-of-bounds particles.
 * 5. Migrates particles between MPI processors using `PerformSingleParticleMigrationCycle`.
 * 6. Advances simulation time and step counter.
 * 7. Locates particles in the grid based on their *new* positions.
 * 8. Interpolates the fluid velocity (from the *current* Ucat) to the new particle locations
 *    to get velocities for the *next* advection step.
 * 9. Scatters particle data back to Eulerian fields.
 *10. Logs errors and outputs data at specified `OutputFreq` intervals.
 *
 * @param user         Pointer to the UserCtx structure.
 * @param StartStep    The initial step number (e.g., 0 for a new run, >0 for restart).
 * @param StartTime    The simulation time corresponding to StartStep.
 * @param StepsToRun   The number of steps to execute in this run. If 0 and StartStep is 0,
 *                     only `PerformInitialSetup` is executed.
 * @param OutputFreq   Frequency (in number of steps) at which to output data and log errors.
 * @param readFields   Flag indicating whether to read initial fields (used by `SetEulerianFields`
 *                     and `PerformInitialSetup`).
 * @param bboxlist     Array of BoundingBox structures for domain decomposition, used for migration.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode AdvanceSimulation(UserCtx *user, PetscInt StartStep, PetscReal StartTime,
                                 PetscInt StepsToRun, PetscInt OutputFreq, PetscBool readFields,
                                 const BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;            // MPI rank of the current process
    PetscReal      dt = user->dt;   // Timestep size
    PetscInt       step_loop_counter; // Loop counter for simulation steps
    PetscReal      currentTime;     // Current simulation time
    PetscInt       removed_local, removed_global; // Counters for out-of-bounds particles

    // Variables for particle migration within the main loop
    MigrationInfo  *migrationList_main_loop = NULL;    // Managed by AdvanceSimulation for the loop
    PetscInt       migrationCount_main_loop = 0;
    PetscInt       migrationListCapacity_main_loop = 0;
    PetscInt       globalMigrationCount_in_loop_cycle; // Output from migration cycle

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting simulation run: %d steps from step %d (t=%.4f), dt=%.4f\n",
              StepsToRun, StartStep, StartTime, dt);

    currentTime = StartTime; // Initialize current simulation time

    // --- Handle Initial Setup (if StartStep is 0) ---
    if (StartStep == 0) {
        // Set Eulerian fields for t=0 before particle setup
        ierr = SetEulerianFields(user, StartStep, StartStep, currentTime, readFields); CHKERRQ(ierr);
        // Perform comprehensive initial particle setup
        ierr = PerformInitialSetup(user, currentTime, StartStep, readFields, bboxlist, OutputFreq, StepsToRun, StartStep); CHKERRQ(ierr);

        // If only initial setup was requested (StepsToRun == 0), exit now.
        if (StepsToRun == 0) {
            LOG_ALLOW(GLOBAL, LOG_INFO, "Initial setup completed as StepsToRun is 0. Time t=%.4f. Exiting AdvanceSimulation.\n", currentTime);
            // The migrationList_main_loop is not yet used, but free it for consistency.
            ierr = PetscFree(migrationList_main_loop); CHKERRQ(ierr);
            PetscFunctionReturn(0);
        }
    }

    // --- Time Marching Loop ---
    // Loops from the StartStep up to (but not including) StartStep + StepsToRun
    for (step_loop_counter = StartStep; step_loop_counter < StartStep + StepsToRun; step_loop_counter++)
    {
      LOG_ALLOW(GLOBAL, LOG_INFO, "Starting step %d, t=%.4f\n", step_loop_counter, currentTime);

      // Set/Update Eulerian Fields for the current time step
      // If StartStep was 0, fields for t=0 were already set.
      // If restarting (StartStep > 0), or for subsequent steps (step_loop_counter > StartStep),
      // set/update fields for the current `currentTime`.
      if (step_loop_counter > StartStep || (step_loop_counter == StartStep && StartStep > 0)) {
          ierr = SetEulerianFields(user, step_loop_counter, StartStep, currentTime, readFields); CHKERRQ(ierr);
      }

      // Step 3: Update Particle Positions
      // Moves particles from P(currentTime) to P(currentTime + dt) using velocity_particle(currentTime).
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f -> T=%.4f, Step=%d] Updating particle positions (Rank: %d).\n", currentTime, currentTime+dt, step_loop_counter, rank);
      ierr = UpdateAllParticlePositions(user); CHKERRQ(ierr);

      // Step 4: Check Boundaries and Remove Out-Of-Bounds Particles
      ierr = CheckAndRemoveOutOfBoundsParticles(user, &removed_local, &removed_global); CHKERRQ(ierr);
      if (removed_global > 0) {
          LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Removed %d out-of-bounds particles globally.\n", currentTime, step_loop_counter, removed_global);
      }

      // Step 5 & 6: Migrate Particles between processors
      // Migration is based on new positions P(currentTime + dt).
      // The time logged for migration is the target time of the new positions.
      ierr = PerformSingleParticleMigrationCycle(user, bboxlist,
                                                 &migrationList_main_loop, &migrationCount_main_loop, &migrationListCapacity_main_loop,
                                                 currentTime + dt, step_loop_counter, "Main Loop", &globalMigrationCount_in_loop_cycle); CHKERRQ(ierr);
      // migrationCount_main_loop is reset inside PerformSingleParticleMigrationCycle for the *next* call within this loop.

      // Step 7: Advance Time
      currentTime += dt; // currentTime now represents the time at the *end* of the step just completed.

      // Step 8: Update global step counter in user context.
      // This should reflect the step number *completed*.
      user->step = step_loop_counter + 1;

      // Step 9: Locate Particles at New Positions (t = currentTime)
      // This is for particles now on the current rank after migration, using their P(currentTime) positions.
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d completed] Locating particles at new positions (Rank: %d).\n", currentTime, step_loop_counter, rank);
      ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);

      // Step 10: Interpolate Eulerian Field to New Particle Positions
      // Eulerian field Ucat is typically from the beginning of the current step (time effectively currentTime - dt).
      // Interpolate this Ucat to P(currentTime) to get V_particle(currentTime) for the *next* advection.
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d completed] Interpolating field (from Ucat at t=%.4f) to new particle positions (Rank: %d).\n", currentTime, step_loop_counter, currentTime-dt, rank);
      ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr);

      // Step 11: Scatter Particle Fields back to Eulerian Grid (if applicable)
      ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d completed] Scattering particle fields to grid (Rank: %d).\n", currentTime, step_loop_counter, rank);

      // Step 12: Output and Error Logging (based on the step *just completed*, using user->step)
      if (OutputFreq > 0 && (user->step % OutputFreq == 0) ) {
          LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d completed (Output for step %d)] Logging interpolation error.\n", currentTime, step_loop_counter, user->step);
          if(user->FieldInitialization==1) { // If analytic fields are available for comparison
              ierr = LOG_INTERPOLATION_ERROR(user); CHKERRQ(ierr);
          }
          LOG_ALLOW(GLOBAL, LOG_DEBUG, "Printing particle fields at t=%.4f (step %d completed, output for step %d)\n", currentTime, step_loop_counter, user->step);
          ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);

          LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d completed] Writing particle data (for step %d).\n", currentTime, step_loop_counter, user->step);
          ierr = WriteSwarmField(user, "position", user->step, "dat"); CHKERRQ(ierr); // Write P(t+dt)
          ierr = WriteSwarmField(user, "velocity", user->step, "dat"); CHKERRQ(ierr); // Write V_interp @ P(t+dt)
          // ierr = WriteSwarmField(user, "pos_phy",  user->step, "dat"); CHKERRQ(ierr);
          ierr = WriteSimulationFields(user); CHKERRQ(ierr); // Optional grid field write
      }
      LOG_ALLOW(GLOBAL, LOG_INFO, "Finished step %d, current time t=%.4f (user->step is now %d)\n", step_loop_counter, currentTime, user->step);
    } // End of time marching loop

    PetscReal finalTimeRun = StartTime + StepsToRun * dt; // Target final time if all steps ran
    LOG_ALLOW(GLOBAL, LOG_INFO, "Time marching completed. %d steps run from StartStep %d. Current time t=%.4f (target final: %.4f).\n", StepsToRun, StartStep, currentTime, finalTimeRun);

    // --- Final Particle Field Log (if not already done by OutputFreq and if steps were run) ---
    if (StepsToRun > 0 && !(OutputFreq > 0 && (user->step % OutputFreq == 0))) {
        // Avoid double logging if the last step was an output step.
        // Only log if actual simulation steps were performed.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Printing final particle fields at end of run (t=%.4f):\n", currentTime);
        ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);
    }

    // --- Cleanup migration list used in the main loop ---
    ierr = PetscFree(migrationList_main_loop); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/**
 * @brief Initializes or updates all necessary simulation fields for a given timestep.
 *
 * This function handles the logic for either:
 * A) Initializing fields analytically (for the first step or if not reading):
 *    - Sets interior values using SetAnalyticalCartesianField.
 *    - Applies boundary conditions using ApplyAnalyticalBC.
 *    - Updates local vectors with ghosts using UpdateLocalGhosts.
 *    - Optionally writes the initial fields.
 * B) Reading fields from a file for a specific timestep index:
 *    - Reads global vectors using ReadSimulationFields.
 *    - Updates local vectors with ghosts using UpdateLocalGhosts.
 * C) Updating fields using a fluid solver (Placeholder for future integration):
 *    - Calls a placeholder function SolveFluidEquations.
 *    - Applies boundary conditions using ApplyAnalyticalBC.
 *    - Updates local vectors with ghosts using UpdateLocalGhosts.
 *
 * @param user        Pointer to the UserCtx structure.
 * @param step        The current timestep number (0 for initial step).
 * @param time        The current simulation time.
 * @param readFields  Flag indicating whether to read fields from file.
 * @param fieldSource Source for field data (e.g., ANALYTICAL, FILE, SOLVER).
 *                    (Here using readFields bool for simplicity based on original code)
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode SetEulerianFields(UserCtx *user, PetscInt step, PetscInt StartStep, PetscReal time, PetscBool readFields)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Updating simulation fields.\n", time, step);

    if (step == StartStep){
      if(!readFields) { 
        // --- Initial Analytical Setup (Step 0, Not Reading) ---
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Setting analytical interior fields.\n", time, step);
        ierr = SetAnalyticalCartesianField(user,"Ucat"); CHKERRQ(ierr);
        // Add other fields (P, Nvert, Ucont...) if needed
        // ierr = SetAnalyticalCartesianField(user,"P"); CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Applying analytical boundary conditions.\n", time, step);
        ierr = ApplyAnalyticalBC(user, "Ucat"); CHKERRQ(ierr);
        // Add other fields
        // ierr = ApplyAnalyticalBC(user, "P"); CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Updating local ghosts after initial setup.\n", time, step);
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
        // Add other fields
        // ierr = UpdateLocalGhosts(user, "P"); CHKERRQ(ierr);

        // Optional: Write the very initial fields (at step 0)
        // This might be better done once in main after this call returns
        // LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Writing initial simulation fields.\n", time, step);
        // ierr = WriteSimulationFields(user, step); CHKERRQ(ierr); // Pass step for filename

    } else if (readFields) {
        // -- Initial Field Read Setup ----------
	
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Reading fields from file (timestep index %d).\n", time, step, step); // Assuming step is timestep index 'ti'
        ierr = ReadSimulationFields(user, step); CHKERRQ(ierr); // Pass step as timestep index 'ti'

        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Updating local ghosts after reading fields.\n", time, step);
        ierr = UpdateLocalGhosts(user,"Ucat");CHKERRQ(ierr);
        // Update ghosts for all fields read
        // ierr = UpdateLocalGhosts(user, "P"); CHKERRQ(ierr);

       }
    }else {
        // -- Could add another if/else loop if we want to read fields for every timestep!
      if(!readFields){
        // --- Update Fields During Time Marching (step > 0, Not Reading) ---
        // This is where you would call your fluid solver if Ucat wasn't static.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating fields via solver (Placeholder).\n", time, step);

        // Placeholder: In a real simulation, update Ucat, P etc. here
        // ierr = SolveFluidEquations(user, user->dt); CHKERRQ(ierr);

        // For the analytical case, Ucat is constant in time, but we might
        // re-apply BCs if they were time-dependent. If BCs are also static,
        // we might not strictly need to do anything here for Ucat.
        // However, it's good practice to update ghosts if the solver *could* change things.
        // Re-applying static analytical BCs doesn't hurt.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Re-applying boundary conditions.\n", time, step);
        ierr = ApplyAnalyticalBC(user, "Ucat"); CHKERRQ(ierr);
        // ierr = ApplyAnalyticalBC(user, "P"); CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghosts after field update/BCs.\n", time, step);      }
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
        // ierr = UpdateLocalGhosts(user, "P"); CHKERRQ(ierr);
      }

    PetscFunctionReturn(0);
}

/**
 * @brief Computes and stores the Cartesian neighbor ranks for the DMDA decomposition.
 *
 * This function retrieves the neighbor information from the primary DMDA (user->da)
 * and stores the face neighbors (xm, xp, ym, yp, zm, zp) in the user->neighbors structure.
 * It assumes a standard PETSc ordering for the neighbors array returned by DMDAGetNeighbors.
 * Logs warnings if the assumed indices seem incorrect (e.g., center rank mismatch).
 *
 * @param[in,out] user Pointer to the UserCtx structure where neighbor info will be stored.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode ComputeAndStoreNeighborRanks(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    const PetscMPIInt *neighbor_ranks_ptr; // Use const pointer from PETSc

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Computing DMDA neighbor ranks.\n", rank);

    if (!user || !user->da) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx or user->da is NULL in ComputeAndStoreNeighborRanks.");
    }

    // Get the neighbor information from the DMDA
    ierr = DMDAGetNeighbors(user->da, &neighbor_ranks_ptr); CHKERRQ(ierr);

    // --- Extract face neighbors assuming standard PETSc 3D ordering ---
    // Center index in 3x3x3 stencil is 13 (0-based)
    // Z varies slowest, then Y, then X.
    // Index = k*9 + j*3 + i (where i,j,k are relative offsets -1, 0, +1 mapped to 0, 1, 2)
    // Example: (i,j,k)=(0,0,0) relative => index = 1*9 + 1*3 + 1 = 13 (Center)
    // Example: (i,j,k)=(-1,0,0) relative => index = 1*9 + 1*3 + 0 = 12 (xm)
    // Example: (i,j,k)=(+1,0,0) relative => index = 1*9 + 1*3 + 2 = 14 (xp)
    // Example: (i,j,k)=(0,-1,0) relative => index = 1*9 + 0*3 + 1 = 10 (ym)
    // Example: (i,j,k)=(0,+1,0) relative => index = 1*9 + 2*3 + 1 = 16 (yp)
    // Example: (i,j,k)=(0,0,-1) relative => index = 0*9 + 1*3 + 1 = 4  (zm)
    // Example: (i,j,k)=(0,0,+1) relative => index = 2*9 + 1*3 + 1 = 22 (zp)

    if (neighbor_ranks_ptr[13] != rank) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "Rank %d: DMDAGetNeighbors center index mismatch (expected %d, got %d). Neighbor indexing might be incorrect.\n",
                  rank, rank, neighbor_ranks_ptr[13]);
    }

    user->neighbors.rank_xm = neighbor_ranks_ptr[12];
    user->neighbors.rank_xp = neighbor_ranks_ptr[14];
    user->neighbors.rank_ym = neighbor_ranks_ptr[10];
    user->neighbors.rank_yp = neighbor_ranks_ptr[16];
    user->neighbors.rank_zm = neighbor_ranks_ptr[4];
    user->neighbors.rank_zp = neighbor_ranks_ptr[22];

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Neighbors: xm=%d, xp=%d, ym=%d, yp=%d, zm=%d, zp=%d\n", rank,
              user->neighbors.rank_xm, user->neighbors.rank_xp,
              user->neighbors.rank_ym, user->neighbors.rank_yp,
              user->neighbors.rank_zm, user->neighbors.rank_zp);

    // Note: neighbor_ranks_ptr memory is managed by PETSc, do not free it.

    PetscFunctionReturn(0);
}

/**
 * @brief Sets the processor layout for a given DMDA based on PETSc options.
 *
 * Reads the desired number of processors in x, y, and z directions using
 * PETSc options (e.g., -dm_processors_x, -dm_processors_y, -dm_processors_z).
 * If an option is not provided for a direction, PETSC_DECIDE is used for that direction.
 * Applies the layout using DMDASetNumProcs.
 *
 * Also stores the retrieved/decided values in user->procs_x/y/z if user context is provided.
 *
 * @param dm   The DMDA object to configure the layout for.
 * @param user Pointer to the UserCtx structure (optional, used to store layout values).
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode SetDMDAProcLayout(DM dm, UserCtx *user)
{
    PetscErrorCode ierr;
    PetscMPIInt    size, rank;
    PetscInt       px = PETSC_DECIDE, py = PETSC_DECIDE, pz = PETSC_DECIDE;
    PetscBool      px_set = PETSC_FALSE, py_set = PETSC_FALSE, pz_set = PETSC_FALSE;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)dm), &size); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)dm), &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Configuring DMDA processor layout for %d total processes.\n", rank, size);

    // --- Read desired layout from options ---
    // Use names like "-dm_processors_x" which are somewhat standard in PETSc examples
    ierr = PetscOptionsGetInt(NULL, NULL, "-dm_processors_x", &px, &px_set); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dm_processors_y", &py, &py_set); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-dm_processors_z", &pz, &pz_set); CHKERRQ(ierr);

    // --- Validate User Input (Optional but Recommended) ---
    // Check if specified processor counts multiply to the total MPI size
    if (px_set && py_set && pz_set) {
        if (px * py * pz != size) {
             SETERRQ(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_INCOMP,
                     "Specified processor layout %d x %d x %d = %d does not match MPI size %d",
                     px, py, pz, px * py * pz, size);
        }
         LOG_ALLOW(GLOBAL, LOG_INFO, "Using specified processor layout: %d x %d x %d\n", px, py, pz);
    } else if (px_set || py_set || pz_set) {
         // If only some are set, PETSC_DECIDE will be used for others
         LOG_ALLOW(GLOBAL, LOG_INFO, "Using partially specified processor layout: %d x %d x %d (PETSC_DECIDE for unspecified)\n", px, py, pz);
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Using fully automatic processor layout (PETSC_DECIDE x PETSC_DECIDE x PETSC_DECIDE)\n");
    }
     // Additional checks: Ensure px, py, pz are positive if set
     if ((px_set && px <= 0) || (py_set && py <= 0) || (pz_set && pz <= 0)) {
         SETERRQ(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_OUTOFRANGE, "Specified processor counts must be positive.");
     }


    // --- Apply the layout to the DMDA ---
    ierr = DMDASetNumProcs(dm, px, py, pz); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Rank %d: DMDASetNumProcs called with px=%d, py=%d, pz=%d.\n", rank, px, py, pz);


    // --- Store the values in UserCtx (Optional) ---
    // Note: If PETSC_DECIDE was used, PETSc calculates the actual values during DMSetUp.
    // We store the *requested* values here. To get the *actual* values used,
    // you would need to call DMDAGetInfo after DMSetUp.
    /*
    if (user) {
        user->procs_x = px;
        user->procs_y = py;
        user->procs_z = pz;
    }
    */
    
    PetscFunctionReturn(0);
}

/**
 * @brief Sets up the full rank communication infrastructure, including neighbor ranks and bounding box exchange.
 *
 * This function orchestrates the following steps:
 * 1. Compute and store the neighbor ranks in the user context.
 * 2. Gather all local bounding boxes to rank 0.
 * 3. Broadcast the complete bounding box list to all ranks.
 *
 * The final result is that each rank has access to its immediate neighbors and the bounding box information of all ranks.
 *
 * @param[in,out] user      Pointer to the UserCtx structure (must be initialized).
 * @param[in,out] bboxlist  Pointer to BoundingBox array pointer; after this call, it will point to the broadcasted list.
 *
 * @return PetscErrorCode Returns 0 on success or non-zero PETSc error code.
 */
PetscErrorCode SetupDomainRankInfo(UserCtx *user, BoundingBox **bboxlist)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    LOG_ALLOW(GLOBAL, LOG_INFO, "SetupDomainRankCommunications: Starting full rank communication setup.\n");

    // Step 1: Compute and store neighbor ranks
    ierr = ComputeAndStoreNeighborRanks(user); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "SetupDomainRankCommunications: Neighbor ranks computed and stored.\n");

    // Step 2: Gather all local bounding boxes on rank 0
    ierr = GatherAllBoundingBoxes(user, bboxlist); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "SetupDomainRankCommunications: Bounding boxes gathered on rank 0.\n");

    // Step 3: Broadcast bounding box list to all ranks
    ierr = BroadcastAllBoundingBoxes(user, bboxlist); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "SetupDomainRankCommunications: Bounding boxes broadcasted to all ranks.\n");

    LOG_ALLOW(GLOBAL, LOG_INFO, "SetupDomainRankCommunications: Completed successfully.\n");

    PetscFunctionReturn(0);
}
