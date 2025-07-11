/**
 * @file simulation.c  // code for simulation loop 
 * @brief Test program for DMSwarm interpolation using the fdf-curvIB method.
 * Provides the setup to start any simulation with DMSwarm and DMDAs.
 **/

#include "simulation.h"


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
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetEulerianFields(UserCtx *user, PetscInt step, PetscInt StartStep, PetscReal time)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Preparing complete Eulerian state.\n", time, step);

    PetscReal umax=0.0;
    PetscReal ucont_max=0.0;
    PetscReal umin=0.0;

    // ==============================================================================
    // --- STEP 1: Update the INTERIOR of the domain based on the simulation phase ---
    // ==============================================================================

    if (step == StartStep && StartStep > 0) {
        // --- PATH 1: RESTART from file ---
        // This is the first time this function is called in a restarted run.
        LOG_ALLOW(GLOBAL, LOG_INFO, "RESTART condition: Reading all grid fields from file for step %d.\n", step);
        ierr = ReadSimulationFields(user, step); CHKERRQ(ierr); // Assumes this function reads Ucat, Ucont, etc.

        // After loading, we MUST update local ghosts to ensure consistency for any subsequent calculations.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghost regions for all fields after reading.\n", time, step);
        ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);

    } else {
        // --- PATH 2 & 3: FRESH START or TIME ADVANCEMENT ---
        // This block handles both generating initial fields and advancing the solver in time.

        if (step == 0) { // Condition is now simply step == 0 for a fresh start
            // --- PATH 2: Initial Field Setup (t=0) ---
            LOG_ALLOW(GLOBAL, LOG_INFO, "FRESH START: Generating INTERIOR fields for initial step 0.\n");
            ierr = SetInitialInteriorField(user, "Ucont"); CHKERRQ(ierr);
        } else {
            // --- PATH 3: Advancing the simulation (step > 0) ---
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "TIME ADVANCE: Updating INTERIOR fields for step %d.\n", step);
            // This is the hook for the actual fluid dynamics solver.
            // ierr = YourNavierStokesSolver(user, user->dt); CHKERRQ(ierr);
	    // ierr = SetInitialInteriorField(user, "Ucont"); CHKERRQ(ierr);
        }

        // The following logic is common to both fresh starts and time advancement,
        // but not to a file-based restart (which loads the final Ucat directly).

        // STEP 2: APPLY BOUNDARY CONDITIONS
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Executing boundary condition system.\n", time, step);
        ierr = BoundarySystem_ExecuteStep(user); CHKERRQ(ierr);

        // STEP 3: SYNCHRONIZE Ucont BEFORE CONVERSION
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghost regions for Ucont.\n", time, step);
        ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);

        // STEP 4: CONVERT CONTRAVARIANT TO CARTESIAN
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Converting Ucont to Ucat.\n", time, step);
        ierr = Contra2Cart(user); CHKERRQ(ierr);

        // STEP 5: Re-apply BCs and SYNCHRONIZE Ucat
        // It's often necessary to apply BCs again to ensure Ucat is correct at boundaries.
        ierr = BoundarySystem_ExecuteStep(user); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghost regions for Ucat.\n", time, step);
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Complete Eulerian state is now finalized and consistent.\n", time, step);
    PetscFunctionReturn(0);
}

/**
 * @brief Performs the complete initial setup for the particle simulation at time t=0. [TEST VERSION]
 *
 * This version uses the new, integrated `LocateAllParticlesInGrid_TEST` orchestrator,
 * which handles both location and migration in a single, robust, iterative process.
 *
 * Its sequential operations are:
 * 1. A single, comprehensive call to `LocateAllParticlesInGrid_TEST` to sort all particles
 *    to their correct owner ranks and find their initial host cells.
 * 2. If `user->ParticleInitialization == 0` (Surface Init), it re-initializes particles on the
 *    designated inlet surface, now that they are on the correct MPI ranks.
 * 3. A second call to `LocateAllParticlesInGrid_TEST` is needed after re-initialization to
 *    find the new, correct host cells for the surface-placed particles.
 * 4. Interpolates initial Eulerian fields to the settled particles.
 * 5. Scatters particle data to Eulerian fields (if applicable).
 * 6. Outputs initial data if requested.
 *
 * @param user Pointer to the UserCtx structure.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformInitialSetup_TEST(UserCtx *user, PetscReal currentTime, PetscInt step,
                                        PetscInt OutputFreq,
                                        PetscInt StepsToRun, PetscInt StartStep,
					BoundingBox *bboxlist)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Performing initial particle setup procedures [TEST].\n", currentTime, step);

    // --- 1. Initial Particle Settlement (Location and Migration) ---
    // This single call replaces the old sequence of Locate -> Migrate. The new
    // orchestrator handles the iterative process internally until all particles are settled.
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Initial Settlement: Locating and migrating all particles to their correct ranks and cells.\n", currentTime, step);
    ierr = LocateAllParticlesInGrid_TEST(user,bboxlist); CHKERRQ(ierr);

    // --- 2. Re-initialize Particles on Inlet Surface (if applicable) ---
    if (user->ParticleInitialization == 0 && user->inletFaceDefined) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Re-initializing particles on inlet surface now that they are on correct ranks.\n", currentTime, step);
        ierr = ReinitializeParticlesOnInletSurface(user, currentTime, step); CHKERRQ(ierr);

        // --- CRITICAL: After re-placing particles, we MUST locate them again. ---
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Post-Reinitialization Settlement: Finding host cells for newly placed inlet particles.\n", currentTime, step);
        ierr = LocateAllParticlesInGrid_TEST(user,bboxlist); CHKERRQ(ierr);
    }
    
    // --- 3. Finalize State for t=0 ---
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Interpolating initial fields to settled particles.\n", currentTime, step);
    ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr);
    ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);

    // --- 4. Initial Output ---
    if (OutputFreq > 0 || (StepsToRun == 0 && StartStep == 0)) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Writing initial simulation data.\n", currentTime, step);
        ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);
        ierr = WriteSwarmField(user, "position", step, "dat"); CHKERRQ(ierr);
        ierr = WriteSwarmField(user, "velocity", step, "dat"); CHKERRQ(ierr);
	//  if (!readFields) {
        ierr = WriteSimulationFields(user); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


/**
 * @brief Executes the main time-marching loop for the particle simulation. [TEST VERSION]
 *
 * This version uses the new, integrated `LocateAllParticlesInGrid_TEST` orchestrator
 * and the `ResetAllParticleStatuses` helper for a clean, robust, and understandable workflow.
 *
 * For each timestep, it performs:
 *  1. Sets the background fluid velocity field (Ucat) for the current step.
 *  2. Updates particle positions using velocity from the *previous* step's interpolation.
 *  3. Removes any particles that have left the global domain.
 *  4. A single call to `LocateAllParticlesInGrid_TEST`, which handles all
 *     particle location and migration until the swarm is fully settled.
 *  5. Interpolates the current fluid velocity to the newly settled particle locations.
 *  6. Scatters particle data back to Eulerian fields.
 *  7. Outputs data at specified intervals.
 *
 * @param user       Pointer to the UserCtx structure.
 * @param StartStep  Index of the first step to execute (0-based).
 * @param StartTime  Physical time at StartStep.
 * @param StepsToRun Number of timesteps to advance.
 * @param OutputFreq Frequency (in steps) at which to write output (0 = no output).
 * @param bboxlist   Array of bounding boxes for out-of-bounds checking.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode AdvanceSimulation_TEST(UserCtx *user,
                                      PetscInt  StartStep,
                                      PetscReal StartTime,
                                      PetscInt  StepsToRun,
                                      PetscInt  OutputFreq,
                                      BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    const PetscReal dt          = user->dt;
    PetscReal       currentTime = StartTime;
    PetscInt        removed_local_ob, removed_global_ob;
    PetscInt        removed_local_lost, removed_global_lost;
    PetscInt        removed_local,removed_global;
    PetscInt        output_step;

    PetscFunctionBeginUser;
    LOG_ALLOW(GLOBAL, LOG_INFO,
              "Starting simulation run [PRODUCTION]: %d steps from step %d (t=%.4f), dt=%.4f\n",
              StepsToRun, StartStep, StartTime, dt);

    // --- Handle Initial Setup (t = StartTime, step = StartStep) ---
    if (StartStep == 0) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "--- Preparing state at t=%.4f (Step 0) ---\n", currentTime);
        user->step = 0;

        ierr = SetEulerianFields(user, 0, StartStep, currentTime); CHKERRQ(ierr);
        ierr = PerformInitialSetup_TEST(user, currentTime, 0, OutputFreq, StepsToRun, StartStep, bboxlist);
        CHKERRQ(ierr);

        if (StepsToRun == 0) {
            LOG_ALLOW(GLOBAL, LOG_INFO, "Initial setup completed. No steps to run. Exiting.\n");
            PetscFunctionReturn(0);
        }
        // NOTE: do not advance currentTime here
        LOG_ALLOW(GLOBAL, LOG_INFO, "--- Initial setup complete. Beginning time marching. ---\n\n");
    }

    // --- Time Marching Loop ---
    for (PetscInt step = StartStep; step < StartStep + StepsToRun; ++step) {
        output_step = step + 1;
        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "--- Advancing from Step %d (t=%.4f) to Step %d (t=%.4f) ---\n",
                  step, currentTime, output_step, currentTime + dt);

        // 1) Reset statuses
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Resetting particle statuses for next timestep.\n");
        ierr = ResetAllParticleStatuses(user); CHKERRQ(ierr);

        // 2) Update Eulerian fields for time t
        ierr = SetEulerianFields(user, step, StartStep, currentTime); CHKERRQ(ierr);

        // 3) Advect particles P(t→t+dt)
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f] Updating particle positions to T=%.4f.\n",
                  currentTime, currentTime + dt);
        ierr = UpdateAllParticlePositions(user); CHKERRQ(ierr);

	// 4) Settle particles (location + migration)
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f] Settling all particles.\n", currentTime + dt);
        ierr = LocateAllParticlesInGrid_TEST(user, bboxlist); CHKERRQ(ierr);
	
        // 5a) Remove Lost Particles
	ierr = CheckAndRemoveLostParticles(user,&removed_local_lost,&removed_global_lost);CHKERRQ(ierr);
	// 5b) Remove Out-of-Bounds Particles
        ierr = CheckAndRemoveOutOfBoundsParticles(user, &removed_local_ob, &removed_global_ob, bboxlist);
	// 5c) Accumulate all removed particles
	removed_local = removed_local_lost + removed_local_ob;
	removed_global = removed_global_lost + removed_global_ob; 
        CHKERRQ(ierr);
        if (removed_global > 0) {
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "[T=%.4f] Removed %d out-of-bounds particles globally.\n",
                      currentTime + dt, removed_global);
        }


        // 6) Interpolate & scatter
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f] Interpolating to settled particles.\n", currentTime + dt);
        ierr = InterpolateAllFieldsToSwarm(user);            CHKERRQ(ierr);
        ierr = ScatterAllParticleFieldsToEulerFields(user);  CHKERRQ(ierr);

        // 7) Advance time and step count
        currentTime += dt;
        user->step = output_step;

        // 8) Output if requested
        if (OutputFreq > 0 && (user->step % OutputFreq) == 0) {
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "[T=%.4f] Writing output at Step %d.\n",
                      currentTime, user->step);
            ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);
            ierr = WriteSwarmField(user, "position", user->step, "dat");          CHKERRQ(ierr);
            ierr = WriteSwarmField(user, "velocity", user->step, "dat");          CHKERRQ(ierr);
            ierr = WriteSimulationFields(user);                                   CHKERRQ(ierr);
        }

        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "--- Completed Step %d at t=%.4f ---\n\n", step, currentTime);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Time marching completed. Final time t=%.4f.\n", currentTime);
    PetscFunctionReturn(0);
}
