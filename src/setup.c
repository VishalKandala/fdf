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
    UserCtx*       local_user; // Temporary local pointer to the user context;    

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
    
    ierr = PetscOptionsGetReal(NULL, NULL, "-ucont_x", &(local_user->InitialConstantContra.x), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-ucont_y", &(local_user->InitialConstantContra.y), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-ucont_z", &(local_user->InitialConstantContra.z), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(NULL, NULL, "-read_fields", readFields_out, NULL); CHKERRQ(ierr);

    ierr = PetscOptionsGetBool(NULL, NULL, "-Setup_Only", OnlySetup, NULL); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL,LOG_DEBUG," -- Console Output Functions [Total : %d] : --\n",*nAllowed_out);
    for (PetscInt i = 0; i < *nAllowed_out; ++i) {
      LOG_ALLOW(GLOBAL,LOG_DEBUG,"   [%2d] «%s»\n", i, (*allowedFuncs_out)[i]);
    }

    // --- 3. PETSc Logging Setup ---
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);
    ierr = registerEvents(); CHKERRQ(ierr); // Assuming registerEvents is defined elsewhere

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "InitializeSimulation - Completed successfully.\n");
    PetscFunctionReturn(0);
}

/**
 * @brief Setup grid DMs, field vectors, and compute grid metrics for the simulation block.
 *
 * This function orchestrates the setup for the simulation grid block. It:
 * 1. Defines the grid geometry and physical coordinates using DefineGridCoordinates.
 * 2. Creates the primary global and local PETSc Vecs for flow variables (Ucat, P, etc.).
 * 3. Creates the Vecs required for storing grid metrics.
 * 4. Initializes all newly created vectors to zero.
 * 5. Calls subroutines to compute and populate the grid metrics (Csi, Eta, Zet, Aj).
 *    This includes interior calculation and boundary extrapolation.
 *
 * @param[in,out] user           Pointer to the single UserCtx structure for the simulation.
 * @param[in]     block_number   The number of grid blocks (expected to be 1 for this setup).
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode SetupGridAndVectors(UserCtx *user, PetscInt block_number) {
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx pointer is null.");
    // This function assumes it's handling a single block context passed from main.
    if (block_number != 1) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "SetupGridAndVectors is optimized for block_number=1 but was called with %d. The loop will run, but ensure `user` is an array.", block_number);
    }

    // --- 1. Define Grid DMs and Assign Physical Coordinates ---
    // DefineGridCoordinates is assumed to work on the single 'user' context passed in.
    // If it were truly multi-block, its signature would be `DefineGridCoordinates(user, block_number)`.
    // Let's assume it sets up user->da, user->fda, and coordinates correctly.
    ierr = DefineGridCoordinates(user); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Grid DMs and coordinates defined.\n");

    // --- 2. Get local info and store it in the context ---
    ierr = DMDAGetLocalInfo(user->da, &user->info); CHKERRQ(ierr);

    // Check that the DMs are valid
    if (!user->fda || !user->da) {
       SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "DMs not initialized after DefineGridCoordinates");
    }

    // --- 3. Create and Initialize All Vectors ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Creating and initializing all simulation vectors.");
    // Primary Field Vecs
    ierr = DMCreateGlobalVector(user->fda, &user->Ucat); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(user->fda, &user->Ucont); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(user->da,  &user->P); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(user->da,  &user->Nvert); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(user->da,  &user->Nvert_o); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(user->da,  &user->ParticleCount); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(user->fda, &user->lUcat); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(user->fda, &user->lUcont); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(user->da,  &user->lP); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(user->da,  &user->lNvert); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(user->da,  &user->lNvert_o); CHKERRQ(ierr);

    // Metric Vecs
    ierr = DMCreateGlobalVector(user->fda, &user->Csi); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(user->fda, &user->Eta); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(user->fda, &user->Zet); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(user->da,  &user->Aj); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(user->fda, &user->lCsi); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(user->fda, &user->lEta); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(user->fda, &user->lZet); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(user->da,  &user->lAj); CHKERRQ(ierr);
    // (Create ICsi, IAj, etc. here if needed)

    // Zero out global vectors
    ierr = VecSet(user->Ucat, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->Ucont, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->P, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->Nvert, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->Nvert_o, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->ParticleCount, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->Csi, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->Eta, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->Zet, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->Aj, 0.0); CHKERRQ(ierr);

    // Zero out local vectors
    ierr = VecSet(user->lUcat, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->lUcont, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->lP, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->lNvert, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->lNvert_o, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->lCsi, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->lEta, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->lZet, 0.0); CHKERRQ(ierr);
    ierr = VecSet(user->lAj, 0.0); CHKERRQ(ierr);    

    // =========================================================================
    // --- 4. CALL THE METRIC COMPUTATION SUBROUTINES ---
    // =========================================================================
    // These functions now operate directly on the 'user' context, handling
    // all internal steps including getting coordinates, calculation, extrapolation,
    // assembly, and updating the local ghosted vectors.

    ierr = ComputeFaceMetrics(user); CHKERRQ(ierr);
    ierr = ComputeCellCenteredJacobianInverse(user); CHKERRQ(ierr);
    ierr = CheckAndFixGridOrientation(user); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Primary metrics (Csi, Eta, Zet, Aj) computed and prepared.");

    // --- (Optional) Compute face-centered metrics ---
    // ierr = ComputeFaceCenteredMetrics(user); CHKERRQ(ierr);

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "Grid and vectors setup completed on all ranks.\n");
    PetscFunctionReturn(0);
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
	ierr = VecDestroy(&(user[bi].Csi)); CHKERRQ(ierr);
        ierr = VecDestroy(&(user[bi].lCsi)); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," Csi Destroyed \n");
	ierr = VecDestroy(&(user[bi].Eta)); CHKERRQ(ierr);
        ierr = VecDestroy(&(user[bi].lEta)); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," Eta Destroyed \n");
	ierr = VecDestroy(&(user[bi].Zet)); CHKERRQ(ierr);
        ierr = VecDestroy(&(user[bi].lZet)); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," Zet Destroyed \n");
	ierr = VecDestroy(&(user[bi].Aj)); CHKERRQ(ierr);
        ierr = VecDestroy(&(user[bi].lAj)); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," Ucont Destroyed \n");
	
	//  ierr = DMDestroy(&(user[bi].fda)); CHKERRQ(ierr);
	//	LOG_ALLOW(GLOBAL,LOG_DEBUG," fda Destroyed \n");

	LOG_ALLOW(GLOBAL, LOG_INFO, "Finalizing simulation, destroying boundary system.\n");
	ierr = BoundarySystem_Destroy(&user[bi]); CHKERRQ(ierr);
	
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
    } else if (strcmp(fieldName, "Csi") == 0) {
        globalVec = user->Csi;
        localVec  = user->lCsi;
        dm        = user->fda;
    } else if (strcmp(fieldName, "Eta") == 0) {
        globalVec = user->Eta;
        localVec  = user->lEta;
        dm        = user->fda;
    }  else if (strcmp(fieldName, "Zet") == 0) {
        globalVec = user->Zet;
        localVec  = user->lZet;
        dm        = user->fda;
    }else if (strcmp(fieldName, "Nvert") == 0) {
        globalVec = user->Nvert;
        localVec  = user->lNvert;
        dm        = user->da;
     // Add other fields as needed
    } else if (strcmp(fieldName, "Aj") == 0) {
        globalVec = user->Aj;
        localVec  = user->lAj;
        dm        = user->da;
    }else {
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
 * @brief Computes and stores the Cartesian neighbor ranks for the DMDA decomposition.
 *
 * This function retrieves the neighbor information from the primary DMDA (user->da)
 * and stores the face neighbors (xm, xp, ym, yp, zm, zp) in the user->neighbors structure.
 * It assumes a standard PETSc ordering for the neighbors array returned by DMDAGetNeighbors.
 * If DMDAGetNeighbors returns a negative rank that is not MPI_PROC_NULL (which can happen
 * in some PETSc/MPI configurations for non-periodic boundaries if not fully standard),
 * this function will sanitize it to MPI_PROC_NULL to prevent issues.
 *
 * @param[in,out] user Pointer to the UserCtx structure where neighbor info will be stored.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode ComputeAndStoreNeighborRanks(UserCtx *user)
{
    PetscErrorCode    ierr;
    PetscMPIInt       rank;
    PetscMPIInt       size; // MPI communicator size
    const PetscMPIInt *neighbor_ranks_ptr; // Pointer to raw neighbor data from PETSc

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr); // Get MPI size for validation

    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Computing DMDA neighbor ranks.\n", rank);

    if (!user || !user->da) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx or user->da is NULL in ComputeAndStoreNeighborRanks.");
    }

    // Get the neighbor information from the DMDA
    // neighbor_ranks_ptr will point to an internal PETSc array of 27 ranks.
    ierr = DMDAGetNeighbors(user->da, &neighbor_ranks_ptr); CHKERRQ(ierr);

    // Log the raw values from DMDAGetNeighbors for boundary-relevant directions for debugging
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "[CASNR - Rank %d] Raw DMDAGetNeighbors: xm_raw=%d, xp_raw=%d, ym_raw=%d, yp_raw=%d, zm_raw=%d, zp_raw=%d. MPI_PROC_NULL is %d.",
                   rank,
                   neighbor_ranks_ptr[12], neighbor_ranks_ptr[14],
                   neighbor_ranks_ptr[10], neighbor_ranks_ptr[16],
                   neighbor_ranks_ptr[4],  neighbor_ranks_ptr[22],
                   (int)MPI_PROC_NULL);

    // PETSc standard indices for 3D face neighbors from the 27-point stencil:
    // Index = k_offset*9 + j_offset*3 + i_offset (where offsets -1,0,1 map to 0,1,2)
    // Center: (i_off=1, j_off=1, k_off=1) => 1*9 + 1*3 + 1 = 13
    // X-min:  (i_off=0, j_off=1, k_off=1) => 1*9 + 1*3 + 0 = 12
    // X-plus: (i_off=2, j_off=1, k_off=1) => 1*9 + 1*3 + 2 = 14
    // Y-min:  (i_off=1, j_off=0, k_off=1) => 1*9 + 0*3 + 1 = 10
    // Y-plus: (i_off=1, j_off=2, k_off=1) => 1*9 + 2*3 + 1 = 16
    // Z-min:  (i_off=1, j_off=1, k_off=0) => 0*9 + 1*3 + 1 = 4
    // Z-plus: (i_off=1, j_off=1, k_off=2) => 2*9 + 1*3 + 1 = 22

    if (neighbor_ranks_ptr[13] != rank) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "Rank %d: DMDAGetNeighbors center index (13) is %d, expected current rank %d. Neighbor indexing might be non-standard or DMDA small.\n",
                  rank, neighbor_ranks_ptr[13], rank);
        // This warning is important. If the center isn't the current rank, the offsets are likely wrong.
        // However, PETSc should ensure this unless the DM is too small for a 3x3x3 stencil.
    }

    // Assign and sanitize each neighbor rank
    PetscMPIInt temp_neighbor;

    temp_neighbor = neighbor_ranks_ptr[12]; // xm
    if ((temp_neighbor < 0 && temp_neighbor != MPI_PROC_NULL) || temp_neighbor >= size) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "[CASNR - Rank %d] Correcting invalid xm neighbor %d to MPI_PROC_NULL (%d).", rank, temp_neighbor, (int)MPI_PROC_NULL);
        user->neighbors.rank_xm = MPI_PROC_NULL;
    } else {
        user->neighbors.rank_xm = temp_neighbor;
    }

    temp_neighbor = neighbor_ranks_ptr[14]; // xp
    if ((temp_neighbor < 0 && temp_neighbor != MPI_PROC_NULL) || temp_neighbor >= size) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "[CASNR - Rank %d] Correcting invalid xp neighbor %d to MPI_PROC_NULL (%d).", rank, temp_neighbor, (int)MPI_PROC_NULL);
        user->neighbors.rank_xp = MPI_PROC_NULL;
    } else {
        user->neighbors.rank_xp = temp_neighbor;
    }

    temp_neighbor = neighbor_ranks_ptr[10]; // ym
    if ((temp_neighbor < 0 && temp_neighbor != MPI_PROC_NULL) || temp_neighbor >= size) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "[CASNR - Rank %d] Correcting invalid ym neighbor %d to MPI_PROC_NULL (%d).", rank, temp_neighbor, (int)MPI_PROC_NULL);
        user->neighbors.rank_ym = MPI_PROC_NULL;
    } else {
        user->neighbors.rank_ym = temp_neighbor;
    }

    temp_neighbor = neighbor_ranks_ptr[16]; // yp
    if ((temp_neighbor < 0 && temp_neighbor != MPI_PROC_NULL) || temp_neighbor >= size) {
        // The log for index 16 was "zm" in your output, should be yp
        LOG_ALLOW(GLOBAL, LOG_WARNING, "[CASNR - Rank %d] Correcting invalid yp neighbor (raw index 16) %d to MPI_PROC_NULL (%d).", rank, temp_neighbor, (int)MPI_PROC_NULL);
        user->neighbors.rank_yp = MPI_PROC_NULL;
    } else {
        user->neighbors.rank_yp = temp_neighbor;
    }

    temp_neighbor = neighbor_ranks_ptr[4]; // zm
    if ((temp_neighbor < 0 && temp_neighbor != MPI_PROC_NULL) || temp_neighbor >= size) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "[CASNR - Rank %d] Correcting invalid zm neighbor %d to MPI_PROC_NULL (%d).", rank, temp_neighbor, (int)MPI_PROC_NULL);
        user->neighbors.rank_zm = MPI_PROC_NULL;
    } else {
        user->neighbors.rank_zm = temp_neighbor;
    }

    temp_neighbor = neighbor_ranks_ptr[22]; // zp
    if ((temp_neighbor < 0 && temp_neighbor != MPI_PROC_NULL) || temp_neighbor >= size) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "[CASNR - Rank %d] Correcting invalid zp neighbor %d to MPI_PROC_NULL (%d).", rank, temp_neighbor, (int)MPI_PROC_NULL);
        user->neighbors.rank_zp = MPI_PROC_NULL;
    } else {
        user->neighbors.rank_zp = temp_neighbor;
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "[CASNR - Rank %d] Stored user->neighbors: xm=%d, xp=%d, ym=%d, yp=%d, zm=%d, zp=%d\n", rank,
              user->neighbors.rank_xm, user->neighbors.rank_xp,
              user->neighbors.rank_ym, user->neighbors.rank_yp,
              user->neighbors.rank_zm, user->neighbors.rank_zp);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); // Ensure logs are flushed

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting full rank communication setup.\n");

    // Step 1: Compute and store neighbor ranks
    ierr = ComputeAndStoreNeighborRanks(user); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Neighbor ranks computed and stored.\n");

    // Step 2: Gather all local bounding boxes on rank 0
    ierr = GatherAllBoundingBoxes(user, bboxlist); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Bounding boxes gathered on rank 0.\n");

    // Step 3: Broadcast bounding box list to all ranks
    ierr = BroadcastAllBoundingBoxes(user, bboxlist); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Bounding boxes broadcasted to all ranks.\n");

    // Step 4: Setup Domain Cell Composition
    ierr = SetupDomainCellDecompositionMap(user);
    LOG_ALLOW(GLOBAL,LOG_INFO, "Domain Cell Composition set and broadcasted to all ranks. \n");

    LOG_ALLOW(GLOBAL, LOG_INFO, "SetupDomainRankCommunications: Completed successfully.\n");

    PetscFunctionReturn(0);
}

/**
 * @brief Reconstructs Cartesian velocity (Ucat) at cell centers from contravariant
 *        velocity (Ucont) defined on cell faces.
 *
 * This function performs the transformation from a contravariant velocity representation
 * (which is natural on a curvilinear grid) to a Cartesian (x,y,z) representation.
 * For each interior computational cell owned by the rank, it performs the following:
 *
 * 1.  It averages the contravariant velocity components (U¹, U², U³) from the
 *     surrounding faces to get an estimate of the contravariant velocity at the cell center.
 * 2.  It averages the metric vectors (Csi, Eta, Zet) from the surrounding faces
 *     to get an estimate of the metric tensor at the cell center. This tensor forms
 *     the transformation matrix.
 * 3.  It solves the linear system `[MetricTensor] * [ucat] = [ucont]` for the
 *     Cartesian velocity vector `ucat = (u,v,w)` using Cramer's rule.
 * 4.  The computed Cartesian velocity is stored in the global `user->Ucat` vector.
 *
 * The function operates on local, ghosted versions of the input vectors (`user->lUcont`,
 * `user->lCsi`, etc.) to ensure stencils are valid across processor boundaries.
 *
 * @param[in,out] user      Pointer to the UserCtx structure. The function reads from
 *                          `user->lUcont`, `user->lCsi`, `user->lEta`, `user->lZet`, `user->lNvert`
 *                          and writes to the global `user->Ucat` vector.
 *
 * @return PetscErrorCode 0 on success.
 *
 * @note
 *  - This function should be called AFTER `user->lUcont` and all local metric vectors
 *    (`user->lCsi`, etc.) have been populated with up-to-date ghost values via `UpdateLocalGhosts`.
 *  - It only computes `Ucat` for interior cells (not on physical boundaries) and for
 *    cells not marked as solid/blanked by `user->lNvert`.
 *  - The caller is responsible for subsequently applying boundary conditions to `user->Ucat`
 *    and calling `UpdateLocalGhosts(user, "Ucat")` to populate `user->lUcat`.
 */
PetscErrorCode Contra2Cart(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Cmpnts       ***lcsi_arr, ***leta_arr, ***lzet_arr; // Local metric arrays
    Cmpnts       ***lucont_arr;                       // Local contravariant velocity array
    Cmpnts       ***gucat_arr;                        // Global Cartesian velocity array
    PetscReal    ***lnvert_arr;                       // Local Nvert array
    PetscReal    ***laj_arr;                          // Local Jacobian Determinant inverse array

    PetscFunctionBeginUser;
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Starting Contravariant-to-Cartesian velocity transformation.\n");

    // --- 1. Get DMDA Info and Check for Valid Inputs ---
    // All inputs (lUcont, lCsi, etc.) and outputs (Ucat) are on DMs from the UserCtx.
    // We get local info from fda, which governs the layout of most arrays here.
    ierr = DMDAGetLocalInfo(user->fda, &info); CHKERRQ(ierr);
    if (!user->lUcont || !user->lCsi || !user->lEta || !user->lZet || !user->lNvert || !user->Ucat) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Contra2Cart requires lUcont, lCsi/Eta/Zet, lNvert, and Ucat to be non-NULL.");
    }


    // --- 2. Get Read-Only Array Access to Local Input Vectors (with ghosts) ---
    ierr = DMDAVecGetArrayRead(user->fda, user->lUcont, &lucont_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi,   &lcsi_arr);   CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta,   &leta_arr);   CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet,   &lzet_arr);   CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da,  user->lNvert, &lnvert_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da,  user->lAj, &laj_arr); CHKERRQ(ierr);

    // --- 3. Get Write-Only Array Access to the Global Output Vector ---
    // We compute for local owned cells and write into the global vector.
    // PETSc handles mapping the global indices to the correct local memory locations.
    ierr = DMDAVecGetArray(user->fda, user->Ucat, &gucat_arr); CHKERRQ(ierr);


    // --- 4. Define Loop Bounds for INTERIOR Cells ---
    // We use adjusted bounds to avoid calculating Ucat on the physical domain boundaries,
    // as these are typically set explicitly by boundary condition functions.
    // The stencils use indices like i-1, j-1, k-1, so we must start loops at least at index 1.
    PetscInt i_start = (info.xs == 0) ? info.xs + 1 : info.xs;
    PetscInt i_end   = (info.xs + info.xm == info.mx) ? info.xs + info.xm - 1 : info.xs + info.xm;

    PetscInt j_start = (info.ys == 0) ? info.ys + 1 : info.ys;
    PetscInt j_end   = (info.ys + info.ym == info.my) ? info.ys + info.ym - 1 : info.ys + info.ym;

    PetscInt k_start = (info.zs == 0) ? info.zs + 1 : info.zs;
    PetscInt k_end   = (info.zs + info.zm == info.mz) ? info.zs + info.zm - 1 : info.zs + info.zm;

    // --- 5. Main Computation Loop ---
    // Loops over the GLOBAL indices of interior cells owned by this rank.
    for (PetscInt k_cell = k_start; k_cell < k_end; ++k_cell) {
        for (PetscInt j_cell = j_start; j_cell < j_end; ++j_cell) {
            for (PetscInt i_cell = i_start; i_cell < i_end; ++i_cell) {

                // Check if the cell is a fluid cell (not solid/blanked)
	      //    if (lnvert_arr[k_cell][j_cell][i_cell] > 0.1) continue; // Skip solid/blanked cells

                // Transformation matrix [mat] is the metric tensor at the cell center,
                // estimated by averaging metrics from adjacent faces.
                PetscReal mat[3][3];

		PetscReal aj_center = laj_arr[k_cell+1][j_cell+1][i_cell+1];
		
                mat[0][0] = 0.5 * (lcsi_arr[k_cell][j_cell][i_cell-1].x + lcsi_arr[k_cell][j_cell][i_cell].x); //* aj_center;
                mat[0][1] = 0.5 * (lcsi_arr[k_cell][j_cell][i_cell-1].y + lcsi_arr[k_cell][j_cell][i_cell].y); //* aj_center;
                mat[0][2] = 0.5 * (lcsi_arr[k_cell][j_cell][i_cell-1].z + lcsi_arr[k_cell][j_cell][i_cell].z); //* aj_center;

                mat[1][0] = 0.5 * (leta_arr[k_cell][j_cell-1][i_cell].x + leta_arr[k_cell][j_cell][i_cell].x); //* aj_center;
                mat[1][1] = 0.5 * (leta_arr[k_cell][j_cell-1][i_cell].y + leta_arr[k_cell][j_cell][i_cell].y); //* aj_center;
                mat[1][2] = 0.5 * (leta_arr[k_cell][j_cell-1][i_cell].z + leta_arr[k_cell][j_cell][i_cell].z); //* aj_center;

                mat[2][0] = 0.5 * (lzet_arr[k_cell-1][j_cell][i_cell].x + lzet_arr[k_cell][j_cell][i_cell].x); //* aj_center;
                mat[2][1] = 0.5 * (lzet_arr[k_cell-1][j_cell][i_cell].y + lzet_arr[k_cell][j_cell][i_cell].y); //* aj_center;
                mat[2][2] = 0.5 * (lzet_arr[k_cell-1][j_cell][i_cell].z + lzet_arr[k_cell][j_cell][i_cell].z); //* aj_center;
	
                // Contravariant velocity vector `q` at the cell center,
                // estimated by averaging face-based contravariant velocities.
                PetscReal q[3];
                q[0] = 0.5 * (lucont_arr[k_cell][j_cell][i_cell-1].x + lucont_arr[k_cell][j_cell][i_cell].x); // U¹ at cell center
                q[1] = 0.5 * (lucont_arr[k_cell][j_cell-1][i_cell].y + lucont_arr[k_cell][j_cell][i_cell].y); // U² at cell center
                q[2] = 0.5 * (lucont_arr[k_cell-1][j_cell][i_cell].z + lucont_arr[k_cell][j_cell][i_cell].z); // U³ at cell center

                // Solve the 3x3 system `mat * ucat = q` using Cramer's rule.
                PetscReal det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
                                mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
                                mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

		  if (PetscAbsReal(det) < 1.0e-18) {
		      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FLOP_COUNT, "Transformation matrix determinant is near zero at cell (%d,%d,%d) \n", i_cell, j_cell, k_cell);
		                }

                PetscReal det_inv = 1.0 / det;

                PetscReal det0 = q[0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
                                 q[1] * (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]) +
                                 q[2] * (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]);

                PetscReal det1 = -q[0] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
                                  q[1] * (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) -
                                  q[2] * (mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0]);

                PetscReal det2 = q[0] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) -
                                 q[1] * (mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0]) +
                                 q[2] * (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);

                // Store computed Cartesian velocity in the GLOBAL Ucat array at the
                // array index corresponding to the cell's origin node.
                gucat_arr[k_cell][j_cell][i_cell].x = det0 * det_inv;
                gucat_arr[k_cell][j_cell][i_cell].y = det1 * det_inv;
                gucat_arr[k_cell][j_cell][i_cell].z = det2 * det_inv;
            }
        }
    }

    // --- 6. Restore Array Access ---
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lUcont, &lucont_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi,   &lcsi_arr);   CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta,   &leta_arr);   CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet,   &lzet_arr);   CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da,  user->lNvert, &lnvert_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da,  user->lAj, &laj_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Ucat, &gucat_arr); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Completed Contravariant-to-Cartesian velocity transformation. \n");
    PetscFunctionReturn(0);
}

/**
 * @brief Sets up the entire boundary condition system for the simulation.
 *
 * This function is the main entry point for all boundary condition setup. It performs
 * two main tasks:
 *   1. It determines the name of the boundary condition file, using "bcs.dat" as a
 *      default but allowing it to be overridden by the PETSc option `-bcs_file`.
 *   2. It then calls the core `BoundarySystem_Create` function, which reads the file,
 *      creates all the necessary handler objects, and calls their Initialize() methods
 *      to set the initial state of the boundary fields.
 *
 * This function should be called in `main()` AFTER `SetupGridAndVectors()` has completed
 * to ensure that the grid DMs and DMDALocalInfo are valid.
 *
 * @param user The main UserCtx struct.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetupBoundaryConditions(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscBool      bcsFileOptionFound = PETSC_FALSE;
    char           bcs_filename_buffer[PETSC_MAX_PATH_LEN];
    PetscFunctionBeginUser;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting boundary condition system setup...");

    // --- Step 1: Determine the BCs configuration file name ---
    // Set a default name first.
    ierr = PetscStrcpy(bcs_filename_buffer, "bcs.dat"); CHKERRQ(ierr);
    
    // Check for a user-provided override from the command line or a control file.
    ierr = PetscOptionsGetString(NULL, NULL, "-bcs_file", bcs_filename_buffer, sizeof(bcs_filename_buffer), &bcsFileOptionFound); CHKERRQ(ierr);
    
    if (bcsFileOptionFound) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Using user-specified boundary conditions file: '%s'", bcs_filename_buffer);
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO, "No -bcs_file option found. Using default: '%s'", bcs_filename_buffer);
    }

    // --- Step 2: Call the main creator for the boundary system ---
    // This single call will parse the file, create all handlers, and initialize them.
    ierr = BoundarySystem_Create(user, bcs_filename_buffer); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Boundary condition system setup complete.");
    PetscFunctionReturn(0);
}

/**
 * @brief Creates and distributes a map of the domain's cell decomposition to all ranks.
 * @ingroup DomainInfo
 *
 * This function is a critical part of the simulation setup. It determines the global
 * cell ownership for each MPI rank and makes this information available to all
 * other ranks. This "decomposition map" is essential for the robust "Walk and Handoff"
 * particle migration strategy, allowing any rank to quickly identify the owner of a
 * target cell.
 *
 * The process involves:
 * 1. Each rank gets its own node ownership information from the DMDA.
 * 2. It converts this node information into cell ownership ranges using the
 *    `GetOwnedCellRange` helper function.
 * 3. It participates in an `MPI_Allgather` collective operation to build a complete
 *    array (`user->RankCellInfoMap`) containing the ownership information for every rank.
 *
 * This function should be called once during initialization after the primary DMDA
 * (user->da) has been set up.
 *
 * @param[in,out] user Pointer to the UserCtx structure. The function will allocate and
 *                     populate `user->RankCellInfoMap` and set `user->num_ranks`.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 *         Errors can occur if input pointers are NULL or if MPI communication fails.
 */
PetscErrorCode SetupDomainCellDecompositionMap(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  local_node_info;
    RankCellInfo   my_cell_info;
    PetscMPIInt    rank, size;

    PetscFunctionBeginUser;

    // --- 1. Input Validation and MPI Info ---
    if (!user) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx pointer is NULL in SetupDomainCellDecompositionMap.");
    }
    if (!user->da) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "user->da is not initialized in SetupDomainCellDecompositionMap.");
    }

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Setting up domain cell decomposition map for %d ranks.\n", size);

    // --- 2. Determine Local Cell Ownership ---
    // Get the local node ownership information from the primary DMDA.
    ierr = DMDAGetLocalInfo(user->da, &local_node_info); CHKERRQ(ierr);

    // Use the robust helper function to convert node ownership to cell ownership.
    // A cell's index is defined by its origin node.
    ierr = GetOwnedCellRange(&local_node_info, 0, &my_cell_info.xs_cell, &my_cell_info.xm_cell); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&local_node_info, 1, &my_cell_info.ys_cell, &my_cell_info.ym_cell); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&local_node_info, 2, &my_cell_info.zs_cell, &my_cell_info.zm_cell); CHKERRQ(ierr);

    // Log the calculated local ownership for debugging purposes.
    LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] Owns cells: i[%d, %d), j[%d, %d), k[%d, %d)\n",
              rank, my_cell_info.xs_cell, my_cell_info.xs_cell + my_cell_info.xm_cell,
              my_cell_info.ys_cell, my_cell_info.ys_cell + my_cell_info.ym_cell,
              my_cell_info.zs_cell, my_cell_info.zs_cell + my_cell_info.zm_cell);

    // --- 3. Allocate and Distribute the Global Map ---
    // Allocate memory for the global map that will hold information from all ranks.
    ierr = PetscMalloc1(size, &user->RankCellInfoMap); CHKERRQ(ierr);

    // Perform the collective communication to gather the `RankCellInfo` struct from every rank.
    // Each rank sends its `my_cell_info` and receives the complete array in `user->RankCellInfoMap`.
    // We use MPI_BYTE to ensure portability across different systems and struct padding.
    ierr = MPI_Allgather(&my_cell_info, sizeof(RankCellInfo), MPI_BYTE,
                         user->RankCellInfoMap, sizeof(RankCellInfo), MPI_BYTE,
                         PETSC_COMM_WORLD); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Domain cell decomposition map created and distributed successfully.\n");

    PetscFunctionReturn(0);
}

/**
 * @brief Performs a binary search for a key in a sorted array of PetscInt64.
 *
 * This is a standard binary search algorithm implemented as a PETSc-style helper function.
 * It efficiently determines if a given `key` exists within a `sorted` array.
 *
 * @param[in]  n      The number of elements in the array.
 * @param[in]  arr    A pointer to the sorted array of PetscInt64 values to be searched.
 * @param[in]  key    The PetscInt64 value to search for.
 * @param[out] found  A pointer to a PetscBool that will be set to PETSC_TRUE if the key
 *                    is found, and PETSC_FALSE otherwise.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 *
 * @note The input array `arr` **must** be sorted in ascending order for the algorithm
 *       to work correctly.
 */
PetscErrorCode BinarySearchInt64(PetscInt n, const PetscInt64 arr[], PetscInt64 key, PetscBool *found)
{
    PetscInt low = 0, high = n - 1;

    PetscFunctionBeginUser;

    // --- 1. Input Validation ---
    if (!found) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output pointer 'found' is NULL in PetscBinarySearchInt64.");
    }
    if (n > 0 && !arr) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input array 'arr' is NULL for n > 0.");
    }
    
    // Initialize output
    *found = PETSC_FALSE;

    // --- 2. Binary Search Algorithm ---
    while (low <= high) {
        // Use this form to prevent potential integer overflow on very large arrays
        PetscInt mid = low + (high - low) / 2;

        if (arr[mid] == key) {
            *found = PETSC_TRUE; // Key found!
            break;               // Exit the loop
        }
        
        if (arr[mid] < key) {
            low = mid + 1; // Search in the right half
        } else {
            high = mid - 1; // Search in the left half
        }
    }

    PetscFunctionReturn(0);
}
