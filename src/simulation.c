/**
 * @file simulation.c  // code for simulation loop 
 * @brief Test program for DMSwarm interpolation using the fdf-curvIB method.
 * Provides the setup to start any simulation with DMSwarm and DMDAs.
 **/

#include "simulation.h"

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
    PetscReal      uMax = 0.0, uMaxcat = 0.0;

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
    
    // --- 3. Re-initialize Particles on Inlet Surface (if applicable) ---
    // This is crucial for ParticleInitialization == 0 (Surface Init). After particles have been
    // migrated to the rank(s) owning the inlet surface, this step ensures they are properly
    // distributed *on* that surface, rather than remaining at their migrated position (e.g., (0,0,0)).
 
    if (user->ParticleInitialization == 0 && user->inletFaceDefined) {
      ierr = ReinitializeParticlesOnInletSurface(user, currentTime, step); CHKERRQ(ierr);
      // ---   DEBUG
      // if (rank == 0) { // Or the rank that does the re-init
      //  PetscReal *coords_check;
      //  PetscInt nlocal_check;
      //  ierr = DMSwarmGetLocalSize(user->swarm, &nlocal_check); CHKERRQ(ierr);
      //  if (nlocal_check > 0) {
      //     ierr = DMSwarmGetField(user->swarm, "position", NULL, NULL, (void**)&coords_check); CHKERRQ(ierr);
      //	  LOG_ALLOW(LOCAL,LOG_DEBUG, "[Rank %d DEBUG] After Reinit, first particle pos: (%.2f, %.2f, %.2f)\n",
      //	      rank, coords_check[0], coords_check[1], coords_check[2]);
	  // Check for NaNs/Infs in a few particles if suspicious
      //  for (int p_chk = 0; p_chk < PetscMin(5, nlocal_check); ++p_chk) {
      //	    if (PetscIsInfOrNanReal(coords_check[3*p_chk+0]) ||
      //	PetscIsInfOrNanReal(coords_check[3*p_chk+1]) ||
      //	PetscIsInfOrNanReal(coords_check[3*p_chk+2])) {
      //      LOG_ALLOW(LOCAL,LOG_DEBUG, "[Rank %d DEBUG ERROR] Bad coord for particle %d after reinit!\n", rank, p_chk);
      //    }
      //  }
      //  ierr = DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void**)&coords_check); CHKERRQ(ierr);
      //  }
      //  fflush(stdout);
	// }
    // --------- DEBUG	
    }

    LOG_ALLOW(LOCAL,LOG_DEBUG," [ Rank %d] ReinitializeParticlesOnInletSurface completed, heading into PetscBarrier.\n",rank);
    
    // Ensure all ranks wait till Reinitialization is done to proceed.
    ierr = PetscBarrier(NULL);
    
    // --- 4. Final Locate and Interpolate ---
    // Now that particles are on their correct ranks and (if surface init) correctly positioned on the inlet,
    // perform the definitive location and interpolation for t=0.
    LOG_ALLOW(LOCAL, LOG_INFO, "[T=%.4f, Step=%d , Rank=%d] Second Initial Locate & Interpolate (post-migration & re-placement).\n", currentTime, step,rank);
    ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);
    ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr);

    VecNorm(user->Ucont, NORM_INFINITY, &uMax);
    LOG_ALLOW(GLOBAL, LOG_INFO,"[Step %d] Initial  max(|Ucont|) before Scatter  = %.6e\n",step , uMax);
    uMax = 0.0;

    VecNorm(user->Ucat, NORM_INFINITY, &uMaxcat);
    LOG_ALLOW(GLOBAL, LOG_INFO,"[Step %d]  Initial max(|Ucat|) before Scatter  = %.6e\n", step, uMaxcat);
    uMaxcat = 0.0;
    
    // --- 5. Scatter Particle Data to Eulerian Grid (if part of initialization) ---
    ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);

    VecNorm(user->Ucont, NORM_INFINITY, &uMax);
    LOG_ALLOW(GLOBAL, LOG_INFO,"[Step %d] Initial  max(|Ucont|) after Scatter  = %.6e\n",step , uMax);
    uMax = 0.0;

    VecNorm(user->Ucat, NORM_INFINITY, &uMaxcat);
    LOG_ALLOW(GLOBAL, LOG_INFO,"[Step %d]  Initial max(|Ucat|) after Scatter  = %.6e\n", step, uMaxcat);
    uMaxcat = 0.0;
    
    // --- 6. Initial Output ---
    // Output if OutputFreq > 0 OR if it's a setup-only run (StepsToRun==0 and StartStep==0).
    if (OutputFreq > 0 || (StepsToRun == 0 && StartStep == 0)) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Logging initial data and (if applicable) interpolation error.\n", currentTime, step);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Printing Initial particle fields at t=%.4f (step %d completed)\n", currentTime, step);
       
	//  LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] ABOUT TO CALL LOG_PARTICLE_FIELDS.\n", rank);
	ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);
	//  LOG_ALLOW(LOCAL, LOG_INFO, "[Rank %d] RETURNED FROM LOG_PARTICLE_FIELDS.\n", rank); 

	// DEPRECIATED - Testing against sinusoidal (now FieldInitialization = 1 has changed meaning.
	// if(user->FieldInitialization==1) { // If analytic fields are available for comparison
	//    ierr = LOG_INTERPOLATION_ERROR(user); CHKERRQ(ierr);
        //}

        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Writing initial particle data (for step %d).\n", currentTime, step, step);
	//	LOG_ALLOW(LOCAL, LOG_INFO, "[%d] About to call WriteSwarmField for position.\n", rank);
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
 * @brief Initializes or updates the complete, consistent state of all Eulerian fields for a given timestep.
 *
 * This function is a high-level wrapper that orchestrates the entire process of preparing
 * the fluid fields for a single time step. It follows the standard procedure for a
 * curvilinear solver: first resolving contravariant velocities (`Ucont`) and then
 * converting them to Cartesian (`Ucat`).
 *
 * Its sequential operations are:
 * 1.  Update the INTERIOR of the domain:
 *     - For the initial step, it calls `SetInitialInteriorField` to generate values.
 *     - For subsequent steps, it calls the main fluid solver.
 *     - If restarting from a file, it reads the data, overwriting the whole field.
 *
 * 2.  Apply Boundary Conditions:
 *     - It then calls the modular `BoundarySystem_ExecuteStep` to enforce all configured
 *       boundary conditions on the domain edges.
 *
 * 3.  Convert to Cartesian and Finalize:
 *     - It calls `Contra2Cart` to compute `Ucat` from `Ucont`.
 *     - It calls `UpdateLocalGhosts` to ensure all parallel data is synchronized.
 *
 * @param user        Pointer to the UserCtx structure, containing all simulation data.
 * @param step        The current timestep number being processed.
 * @param StartStep   The initial timestep number of the simulation.
 * @param time        The current simulation time.
 * @param readFields  A boolean flag. If true, the simulation attempts to read fields
 *                    from files at the StartStep instead of generating them.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetEulerianFields(UserCtx *user, PetscInt step, PetscInt StartStep, PetscReal time, PetscBool readFields)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Preparing complete Eulerian state.\n", time, step);

    PetscReal umax=0.0;
    PetscReal ucont_max=0.0;
    PetscReal umin=0.0;
    // ==============================================================================
    // --- STEP 1: Update the INTERIOR of the domain ---
    // ==============================================================================

    if (step == StartStep && readFields) {
        // --- RESTART from file ---
        // This case reads the full, previously saved state, overwriting everything.
        LOG_ALLOW(GLOBAL, LOG_INFO, "RESTART condition: Reading all fields from file for step %d.\n", step);
        ierr = ReadSimulationFields(user, step); CHKERRQ(ierr);

	// Even after reading, we must update local ghosts to ensure consistency for subsequent steps.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghost regions for all fields after reading.\n", time, step);
        ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
	
    } else {
        // --- This block handles both initial setup and time-advancement ---
        
          if (step == StartStep) {
            // --- Initial Field Setup ---
            // The boundaries have already been initialized by BoundarySystem_Create.
            // We now fill the domain interior using the newly named function.
	    // ierr = VecZeroEntries(user->Ucont);CHKERRQ(ierr);
	    // ----DEBUG TEST
	    // ierr = VecSet(user->Ucont,1.0);
	    LOG_ALLOW(GLOBAL, LOG_INFO, "INITIAL start: Generating INTERIOR fields for initial step %d.\n", step);
            ierr = SetInitialInteriorField(user, "Ucont"); CHKERRQ(ierr);

	    ierr = VecNorm(user->Ucont, NORM_INFINITY, &ucont_max); CHKERRQ(ierr);
	    LOG_ALLOW(GLOBAL,LOG_INFO,"[DEBUG] max(|Ucont|) after SetInitialInteriorField = %.6e\n", ucont_max);
	    ucont_max=0.0;
	    
	         } else { // Advancing the simulation (step > StartStep)
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "ADVANCE step: Updating INTERIOR fields for step %d.\n", step);


	    ierr = VecNorm(user->Ucont, NORM_INFINITY, &ucont_max); CHKERRQ(ierr);
	    LOG_ALLOW(GLOBAL,LOG_INFO,"[DEBUG] max(|Ucont|) (timestep = %d)  = %.6e\n",step,ucont_max);
	    ucont_max=0.0;
	    // This is the hook for the actual fluid dynamics solver, which would update Ucont.
            // ierr = YourNavierStokesSolver(user, user->dt); CHKERRQ(ierr);
	        }

        // ==============================================================================
        // --- STEP 2: APPLY BOUNDARY CONDITIONS ---
        // ==============================================================================
        // The boundary system applies conditions. For a wall, this sets the normal
        // component of Ucont to 0 and also sets the ghost-cell Ucat values.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Executing boundary condition system.\n", time, step);
	ierr = BoundarySystem_ExecuteStep(user); CHKERRQ(ierr);
	
	// ==============================================================================
        // --- STEP 3: SYNCHRONIZE Ucont BEFORE CONVERSION ---
        // ==============================================================================
        // THIS IS THE CRITICAL FIX: Update lUcont with correct global data and ghost cells
        // BEFORE calling Contra2Cart, which relies on lUcont for its stencil.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghost regions for Ucont.\n", time, step);
        ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);

        // ==============================================================================
        // --- STEP 4: CONVERT CONTRAVARIANT TO CARTESIAN ---
        // ==============================================================================
        // With Ucont fully defined (interior + boundaries), compute the Cartesian velocity.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Converting Ucont to Ucat.\n", time, step);
	//	ierr = VecZeroEntries(user->Ucat); CHKERRQ(ierr);
	// DEBUG--------

        ierr = Contra2Cart(user); CHKERRQ(ierr);
	
	ierr = VecNorm(user->Ucat, NORM_INFINITY, &umax); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_INFO,"[DEBUG] max(|Ucat|) after Contra2Cart = %.6e\n", umax);
	umax = 0.0;

        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Executing boundary condition system.\n", time, step);
	ierr = BoundarySystem_ExecuteStep(user); CHKERRQ(ierr);	
	
	ierr = VecMin(user->Ucat, NULL, &umin);
	LOG_ALLOW(GLOBAL,LOG_DEBUG, "[DEBUG] min(Ucat) after Contra2Cart = %.6e\n", umin);
	// ==============================================================================
        // --- STEP 5: SYNCHRONIZE Ucat AFTER CONVERSION ---
        // ==============================================================================
        // Finally, update the local lUcat with the newly computed global Ucat data.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghost regions for Ucat.\n", time, step);
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
    }
    
    // ==============================================================================
    // --- STEP 4: UPDATE GHOST CELLS (FINALIZATION) ---
    // ==============================================================================
    // This final step synchronizes the ghost cell layers across all MPI ranks for
    // all relevant fields, ensuring a consistent state before they are used.
    //  LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghost regions for all fields.\n", time, step);
    //  ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);
    // ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Complete Eulerian state is now finalized and consistent.\n", time, step);
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
    PetscReal      uMax = 0.0, uMaxcat = 0.0;

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
      //    if (step_loop_counter > StartStep || (step_loop_counter == StartStep && StartStep > 0)) {
      //   ierr = SetEulerianFields(user, step_loop_counter, StartStep, currentTime, readFields); CHKERRQ(ierr);
	  //	    }

      // Step 3: Update Particle Positions
      // Moves particles from P(currentTime) to P(currentTime + dt) using velocity_particle(currentTime).
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f -> T=%.4f, Step=%d] Updating particle positions (Rank: %d).\n", currentTime, currentTime+dt, step_loop_counter, rank);
      ierr = UpdateAllParticlePositions(user); CHKERRQ(ierr);

      // Step 4: Check Boundaries and Remove Out-Of-Bounds Particles
      ierr = CheckAndRemoveOutOfBoundsParticles(user, &removed_local, &removed_global, bboxlist); CHKERRQ(ierr);
      if (removed_global > 0) {
          LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Removed %d out-of-bounds particles globally.\n", currentTime, step_loop_counter, removed_global);
      }

      LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG,"[Rank %d] ENTERING MIGRATION AFTER UPDATE.\n",rank);
      
      // Step 5 & 6: Migrate Particles between processors
      // Migration is based on new positions P(currentTime + dt).
      // The time logged for migration is the target time of the new positions.
      ierr = PerformSingleParticleMigrationCycle(user, bboxlist,
                                                 &migrationList_main_loop, &migrationCount_main_loop, &migrationListCapacity_main_loop,
                                                 currentTime + dt, step_loop_counter, "Main Loop", &globalMigrationCount_in_loop_cycle); CHKERRQ(ierr);
      // migrationCount_main_loop is reset inside PerformSingleParticleMigrationCycle for the *next* call within this loop.

      // Step 7: Advance Time
      currentTime += dt; // currentTime now represents the time at the *end* of the step just completed.

      // ierr = SetEulerianFields(user, step_loop_counter+1, StartStep, currentTime, readFields); CHKERRQ(ierr);

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

      // Step 8: Update global step counter in user context.
      // This should reflect the step number *completed*.
      user->step = step_loop_counter + 1;
      
      // Step 12: Output and Error Logging (based on the step *just completed*, using user->step)
      if (OutputFreq > 0 && (user->step % OutputFreq == 0) ) {
          LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d completed (Output for step %d)] Logging interpolation error.\n", currentTime, step_loop_counter, user->step);

	// DEPRECIATED - Testing against sinusoidal (now FieldInitialization = 1 has changed meaning.
	//  if(user->FieldInitialization==1) { // If analytic fields are available for comparison
        //      ierr = LOG_INTERPOLATION_ERROR(user); CHKERRQ(ierr);
        //  }

	  
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
