// ParticleMotion.c

#include "ParticleMotion.h"

// Define a buffer size for error messages if not already available
#ifndef ERROR_MSG_BUFFER_SIZE
#define ERROR_MSG_BUFFER_SIZE 256 // Or use PETSC_MAX_PATH_LEN if appropriate
#endif

/**
 * @brief Updates a particle's position based on its velocity and the timestep dt (stored in user->dt).
 *
 * @param[in]     user     Pointer to your UserCtx (must contain user->dt).
 * @param[in,out] position Pointer to the particle's current position (Cmpnts).
 * @param[in]     velocity Pointer to the particle's velocity (Cmpnts).
 *
 * @return PetscErrorCode  Returns 0 on success, or an error code on failure.
 */
PetscErrorCode UpdateParticlePosition(UserCtx *user, Cmpnts *position, const Cmpnts *velocity)
{
  PetscFunctionBeginUser; // PETSc macro for error/stack tracing

  /* Update the position with velocity * dt */
  position->x += velocity->x * user->dt;
  position->y += velocity->y * user->dt;
  position->z += velocity->z * user->dt;

  PetscFunctionReturn(0);
}

/**
 * @brief Loops over all local particles in the DMSwarm, updating their positions
 *        based on velocity and the global timestep user->dt.
 *
 * @param[in,out] user    Pointer to UserCtx (must contain dt).
 *
 * @return PetscErrorCode Returns 0 on success, or an error code on failure.
 */
PetscErrorCode UpdateAllParticlePositions(UserCtx *user)
{
  PetscErrorCode ierr;
  DM swarm = user->swarm;
  PetscInt       nLocal, p;
  Cmpnts        *pos = NULL;
  Cmpnts        *vel = NULL;
  PetscMPIInt rank;

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscFunctionBeginUser;  // PETSc macro for error/stack tracing

  // 1) Get the number of local particles
  ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
  if (nLocal == 0) PetscFunctionReturn(0);   /* nothing to do, no fields held */

  // 2) Access the "position" and "velocity" fields
  ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&pos); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&vel); CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL,LOG_DEBUG," [Rank %d] No.of Particles to update: %d.\n",rank,nLocal);

  // 3) Loop over all local particles, updating each position by velocity * dt
  for (p = 0; p < nLocal; p++) {
    ierr = UpdateParticlePosition(user, &pos[p], &vel[p]); CHKERRQ(ierr);
  }

  // 4) Restore the fields
  ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&pos); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&vel); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @brief Checks for particles outside the physical domain boundaries and removes them
 *        using DMSwarmRemovePointAtIndex.
 *
 * This function iterates through all particles local to the current MPI rank.
 * It checks if a particle's position (x, y, or z) is outside the specified
 * physical domain boundaries [xMin, xMax], [yMin, yMax], [zMin, zMax].
 *
 * If a particle is found out of bounds, it is removed using DMSwarmRemovePointAtIndex.
 * NOTE: Removing points changes the indices of subsequent points in the iteration.
 *       Therefore, it's crucial to iterate BACKWARDS or carefully manage indices
 *       after a removal. Iterating backwards is generally safer.
 *
 * @param user    Pointer to the UserCtx structure.
 * @param[out] removedCountLocal Pointer to store the number of particles removed *on this rank*.
 * @param[out] removedCountGlobal Pointer to store the total number of particles removed *across all ranks*.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode CheckAndRemoveOutOfBoundsParticles(UserCtx *user,
                                              PetscInt *removedCountLocal,
                                              PetscInt *removedCountGlobal)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscInt       nLocalInitial, p;
    Cmpnts        *pos = NULL;
    // You only *need* the position field to check bounds
    PetscInt       local_removed_count = 0;
    PetscInt       global_removed_count = 0;
    PetscMPIInt    rank;
    PetscReal xMin,xMax,yMin,yMax,zMin,zMax;

    xMin = user->xMin;
    xMax = user->xMax;

    yMin = user->yMin;
    yMax = user->yMax;

    zMin = user->zMin;
    zMax = user->zMax;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Checking for out-of-bounds particles (using RemovePointAtIndex)...\n", rank);

    ierr = DMSwarmGetLocalSize(swarm, &nLocalInitial); CHKERRQ(ierr);
    if (nLocalInitial == 0) {
        *removedCountLocal = 0;
        *removedCountGlobal = 0; // Will be summed later
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: No local particles to check for removal.\n", rank);
        PetscFunctionReturn(0);
    }

    // Get read-only access to positions
    // Note: We get it before the loop because removal might invalidate pointers
    // if the swarm reallocates internally, although Get/Restore is generally safe.
    // For maximum safety with removal, access inside the loop might be better,
    // but less efficient. Let's try getting it once first.
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
     if (!pos) {
         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Could not access position field.");
     }

    // --- Iterate BACKWARDS to handle index changes during removal ---
    local_removed_count = 0;
    for (p = nLocalInitial - 1; p >= 0; p--) { // Loop from n-1 down to 0
        PetscBool isOutOfBounds = PETSC_FALSE;
        if (pos[p].x < xMin || pos[p].x > xMax ||
            pos[p].y < yMin || pos[p].y > yMax ||
            pos[p].z < zMin || pos[p].z > zMax)
        {
            isOutOfBounds = PETSC_TRUE;
        }

        if (isOutOfBounds) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Removing particle at local index %d at (%g,%g,%g).\n",
                      rank, p, pos[p].x, pos[p].y, pos[p].z);
            // Restore position field BEFORE removing point, as removal modifies the swarm state
            ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
            pos = NULL; // Invalidate pointer

            // --- Remove the particle at the current local index 'p' ---
            ierr = DMSwarmRemovePointAtIndex(swarm, p); CHKERRQ(ierr);
            local_removed_count++;

            // --- Re-get fields if we continue looping (or break) ---
            // Since we removed a point, indices might have shifted if implementation
            // involves compaction. Re-getting fields is safest if loop continues.
            // However, since we iterate backwards, indices p-1, p-2 etc. remain valid
            // relative to the *current* state after removal.
            // Let's try without re-getting inside the loop first for efficiency.
            // We need to re-get 'pos' for the next iteration check.
            if (p > 0) { // Only re-get if there are more iterations
               ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
               if (!pos) {
                   SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Could not re-access position field after removal.");
               }
            }
        }
    } // End backwards loop

    // Restore position field if it was gotten last time inside loop or if loop didn't run/remove
    if (pos) {
         ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
    }


    // Get the *final* local size after removals
    PetscInt nLocalFinal;
    ierr = DMSwarmGetLocalSize(swarm, &nLocalFinal); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Finished removing %d particles. Final local size: %d.\n", rank, local_removed_count, nLocalFinal);


    // Calculate global removed count
    *removedCountLocal = local_removed_count;
    ierr = MPI_Allreduce(&local_removed_count, &global_removed_count, 1, MPIU_INT, MPI_SUM, PetscObjectComm((PetscObject)swarm)); CHKERRQ(ierr);
    *removedCountGlobal = global_removed_count;

    LOG_ALLOW(GLOBAL, LOG_INFO, "CheckAndRemoveOutOfBoundsParticles: Removed %d particles globally.\n", global_removed_count);

    PetscFunctionReturn(0);
}

/**
 * @brief Defines the basic migration pattern for particles within the swarm.
 *
 * This function establishes the migration pattern that dictates how particles
 * move between different MPI ranks in the simulation. It initializes a migration
 * list where each particle is assigned a target rank based on predefined conditions.
 * The migration pattern can be customized to implement various migration behaviors.
 *
 * @param[in,out] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode DefineBasicMigrationPattern(UserCtx* user) {
    DM swarm = user->swarm;           // DMSwarm object managing the particle swarm
    PetscErrorCode ierr;              // Error code for PETSc functions
    PetscMPIInt *miglist;             // Migration list indicating target MPI ranks for particles
    PetscInt localNumParticles;       // Number of particles managed by the local MPI process
    PetscMPIInt rank, size;           // MPI rank of the current process and total number of processes

    // Retrieve the MPI rank of the current process
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    // Retrieve the total number of MPI processes
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL,LOG_INFO,"DefineBasicMigrationPattern - Rank %d out of %d processes.\n", rank, size);

    // Get the number of particles managed by the local MPI process
    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL,LOG_INFO,"DefineBasicMigrationPattern - Rank %d handling %d particles.\n", rank, localNumParticles);

    // Allocate memory for the migration list
    ierr = PetscCalloc1(localNumParticles, &miglist); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG,"DefineBasicMigrationPattern - Allocated migration list for %d particles.\n", localNumParticles);

    // Initialize the migration list: assign each particle to migrate to the current rank by default
    for (PetscInt p = 0; p < localNumParticles; p++) {
        miglist[p] = rank;
    }
    LOG_ALLOW(LOG_DEBUG, LOCAL,"DefineBasicMigrationPattern - Initialized migration list with default rank assignments.\n");

    // Define custom migration conditions based on the number of MPI processes
    if (size > 1) {
        // Example condition: Assign the first particle in rank 0 to migrate to rank 2
        if (rank == 0 && localNumParticles > 0) {
            miglist[0] = 2;
            LOG_ALLOW(LOG_INFO,LOCAL,"DefineBasicMigrationPattern - Rank 0, Particle 0 assigned to migrate to Rank 2.\n");
        }

        // Additional custom conditions can be added here for other ranks
        // Example:
        // if(rank == 1 && localNumParticles > 1){
        //     miglist[1] = 3;
        //     LOG_ALLOW_SYNC(LOG_INFO, "DefineBasicMigrationPattern - Rank 1, Particle 1 assigned to migrate to Rank 3.\n");
        // }

        // ... add more custom conditions as needed ...
    }

    // Assign the migration list to the user context for later use
    user->miglist = miglist;
    LOG_ALLOW(LOG_DEBUG, LOCAL,"DefineBasicMigrationPattern - Migration list assigned to user context.\n");

    /*
    // Optional: Debugging output to verify migration assignments
    for(PetscInt p = 0; p < localNumParticles; p++) {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,
            "DefineBasicMigrationPattern - Rank %d - miglist[%ld] = %ld\n",
            rank, p, miglist[p]);
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "***********************\n");
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    */

    return 0;
}


/**
 * @brief Performs the basic migration of particles based on the defined migration pattern.
 *
 * This function updates the positions of particles within the swarm by migrating them
 * to target MPI ranks as specified in the migration list. It handles the migration process
 * by setting the 'DMSwarm_rank' field for each particle and invokes the DMSwarm migration
 * mechanism to relocate particles across MPI processes. After migration, it cleans up
 * allocated resources and ensures synchronization across all MPI ranks.
 *
 * @param[in,out] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode PerformBasicMigration(UserCtx* user) {
    DM swarm = user->swarm;                // DMSwarm object managing the particle swarm
    PetscErrorCode ierr;                   // Error code for PETSc functions
    PetscMPIInt *miglist;                  // Migration list indicating target MPI ranks for particles
    PetscMPIInt *rankval;                  // Array to store current MPI rank of each particle
    PetscInt localNumParticles;            // Number of particles managed by the local MPI process
    PetscMPIInt rank;                      // MPI rank of the current process
    PetscBool removePoints = PETSC_TRUE;   // Flag indicating whether to remove migrated particles from the local swarm

    // Retrieve the MPI rank of the current process
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL,LOG_INFO,"PerformBasicMigration - Rank %d is initiating migration.\n", rank);

    // Execute the migration pattern to define target ranks for particles
    ierr = DefineBasicMigrationPattern(user); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO,"PerformBasicMigration - Migration pattern defined.\n");

    // Get the number of particles in the local swarm
    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG,"PerformBasicMigration - Rank %d handling %d particles.\n", rank, localNumParticles);

    // Retrieve the migration list from the user context
    miglist = user->miglist;
    LOG_ALLOW(LOG_DEBUG,LOCAL,"PerformBasicMigration - Retrieved migration list from user context.\n");

    /*
    // Optional: Debugging output to verify migration assignments before migration
    for(PetscInt p = 0; p < localNumParticles; p++) {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,
            "PerformBasicMigration - Rank %d - miglist[%ld] = %d; user->miglist[p] = %d\n",
            rank, p, miglist[p], user->miglist[p]);
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "***********************\n");
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    */

    // Access the 'DMSwarm_rank' field from the DMSwarm to update particle ranks
    ierr = DMSwarmGetField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&rankval); CHKERRQ(ierr);
    LOG_ALLOW_SYNC( LOCAL,LOG_DEBUG,"PerformBasicMigration - Retrieved 'DMSwarm_rank' field.\n");

    // Update the 'DMSwarm_rank' field based on the migration list
    for (PetscInt p = 0; p < localNumParticles; p++) {
        rankval[p] = miglist[p];
        LOG_ALLOW_SYNC( LOCAL,LOG_DEBUG,"PerformBasicMigration - Particle %d assigned to Rank %d.\n", p, rankval[p]);
    }

    /*
    // Optional: Debugging output to verify migration assignments after rank updates
    for(PetscInt p = 0; p < localNumParticles; p++) {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,
            "PerformBasicMigration - After change - Rank %ld - rankval[%ld] = %ld; user->miglist[p] = %ld\n",
            rank, p, rankval[p], user->miglist[p]);
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "***********************\n");
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    */

    // Restore the 'DMSwarm_rank' field after modification
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&rankval); CHKERRQ(ierr);
    LOG_ALLOW_SYNC( LOCAL,LOG_DEBUG,"PerformBasicMigration - Restored 'DMSwarm_rank' field.\n");

    // Invoke the DMSwarm migration process to relocate particles based on updated ranks
    ierr = DMSwarmMigrate(swarm, removePoints); CHKERRQ(ierr);
    LOG_ALLOW_SYNC( LOCAL,LOG_INFO,"PerformBasicMigration - DMSwarm migration executed.\n");

    // Free the allocated migration list to prevent memory leaks
    ierr = PetscFree(user->miglist); CHKERRQ(ierr);
    LOG_ALLOW_SYNC( LOCAL,LOG_DEBUG,"PerformBasicMigration - Freed migration list memory.\n");

    // Synchronize all MPI processes to ensure migration completion before proceeding
    ierr = PetscBarrier(NULL); CHKERRQ(ierr);
    LOG_ALLOW_SYNC( LOCAL,LOG_INFO, "PerformBasicMigration - Migration synchronization completed.\n");

    return 0;
}

/**
 * @brief Checks if a particle position is within the bounds of a given bounding box.
 *
 * @param bbox Pointer to the BoundingBox structure.
 * @param pos  Pointer to the particle's position (Cmpnts).
 *
 * @return PetscBool PETSC_TRUE if the particle is inside or on the boundary, PETSC_FALSE otherwise.
 */
static inline PetscBool IsParticleInBox(const BoundingBox *bbox, const Cmpnts *pos) {
    return (pos->x >= bbox->min_coords.x && pos->x <= bbox->max_coords.x &&
            pos->y >= bbox->min_coords.y && pos->y <= bbox->max_coords.y &&
            pos->z >= bbox->min_coords.z && pos->z <= bbox->max_coords.z);
}

/**
 * @brief Identifies particles leaving the local bounding box and finds their target neighbor rank.
 *
 * Iterates local particles, checks against local bounding box `user->bbox`. If outside, checks
 * the pre-computed immediate neighbors (`user->neighbors`) using the global `bboxlist`
 * to see if the particle landed in one of them, prioritizing the exit direction.
 * Populates the `migrationList` with particles found in a neighbor's box.
 * Assumes particles leaving the global domain were handled previously.
 *
 * @param user           Pointer to the UserCtx (contains local bbox and neighbors).
 * @param bboxlist       Array of BoundingBox structs for all ranks (for checking neighbor boxes).
 * @param migrationList  Pointer to an array of MigrationInfo structs (output, allocated/reallocated by this func).
 * @param migrationCount Pointer to the number of particles marked for migration (output).
 * @param listCapacity   Pointer to the current allocated capacity of migrationList (in/out).
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode IdentifyMigratingParticles(UserCtx *user,
                                        const BoundingBox *bboxlist,
                                        MigrationInfo **migrationList,
                                        PetscInt *migrationCount,
                                        PetscInt *listCapacity)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscInt       nLocal, p;
    Cmpnts        *pos = NULL;
    PetscMPIInt    rank;
    BoundingBox    localBBox = user->bbox;
    RankNeighbors  neighbors = user->neighbors; // Use stored neighbors
    PetscInt       currentMigrationCount = 0;
    PetscInt       currentListCapacity = *listCapacity;
    MigrationInfo *currentMigrationList = *migrationList;
    // Add PID pointer if logging PIDs
    // PetscInt64    *pids = NULL;


    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
    if (nLocal == 0) {
        *migrationCount = 0;
        // Ensure output pointers are consistent even if no allocation happened
        *migrationList = currentMigrationList;
        *listCapacity = currentListCapacity;
        PetscFunctionReturn(0);
    }

    // Get read-only access to position
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
    // ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pids); CHKERRQ(ierr); // If logging PIDs
    if (!pos /*|| !pids*/) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Could not access required DMSwarm fields.");


    currentMigrationCount = 0; // Reset count for this call
    for (p = 0; p < nLocal; p++) {
        // Check if particle is OUTSIDE the local bounding box
        if (!IsParticleInBox(&localBBox, &pos[p]))
        {
            PetscInt targetRank = -1; // Target rank not yet found

            // Determine likely exit direction(s) to prioritize neighbor check
            PetscBool exit_xm = pos[p].x < localBBox.min_coords.x;
            PetscBool exit_xp = pos[p].x > localBBox.max_coords.x;
            PetscBool exit_ym = pos[p].y < localBBox.min_coords.y;
            PetscBool exit_yp = pos[p].y > localBBox.max_coords.y;
            PetscBool exit_zm = pos[p].z < localBBox.min_coords.z;
            PetscBool exit_zp = pos[p].z > localBBox.max_coords.z;

            LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Particle %d at (%g,%g,%g) left local bbox. Checking neighbors...\n",
                      rank, p, pos[p].x, pos[p].y, pos[p].z);

            // Check neighbors preferentially
            if (exit_xm && neighbors.rank_xm != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_xm], &pos[p])) targetRank = neighbors.rank_xm;
            else if (exit_xp && neighbors.rank_xp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_xp], &pos[p])) targetRank = neighbors.rank_xp;
            else if (exit_ym && neighbors.rank_ym != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_ym], &pos[p])) targetRank = neighbors.rank_ym;
            else if (exit_yp && neighbors.rank_yp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_yp], &pos[p])) targetRank = neighbors.rank_yp;
            else if (exit_zm && neighbors.rank_zm != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_zm], &pos[p])) targetRank = neighbors.rank_zm;
            else if (exit_zp && neighbors.rank_zp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_zp], &pos[p])) targetRank = neighbors.rank_zp;
            // Add checks for edge/corner neighbors if needed and if they were stored

            // --- Optional Fallback (if strict Cartesian neighbors aren't enough) ---
            // if (targetRank == -1) { /* ... loop through all ranks in bboxlist ... */ }

            if (targetRank != -1) {
                 // Resize list if needed (using PetscRealloc for safety)
                if (currentMigrationCount >= currentListCapacity) {
                    PetscInt newCapacity = (currentListCapacity == 0) ? 16 : currentListCapacity * 2;

		    ierr = PetscRealloc((size_t)newCapacity * sizeof(MigrationInfo),
					(void**)&currentMigrationList); CHKERRQ(ierr);
		    currentListCapacity = newCapacity;
		    
                    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Reallocated migrationList capacity to %d\n", rank, newCapacity);
                 }
                // Add to migration list
                currentMigrationList[currentMigrationCount].local_index = p;
                currentMigrationList[currentMigrationCount].target_rank = targetRank;
                // currentMigrationList[currentMigrationCount].pid = pids[p]; // If storing PID
                currentMigrationCount++;
                LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Particle %d marked for migration to rank %d.\n", rank, p, targetRank);
            } else {
                // Particle left local box but was not found in any *checked* neighbor box.
                // Since CheckAndRemove should have run first, this might indicate a particle
                // moved more than one cell width into a diagonal neighbor's domain not checked here,
                // or there's an issue with BBox overlap/gaps.
                 LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: Particle %d at (%g,%g,%g) left local bbox but target neighbor rank not found! (May be lost or need wider neighbor check).\n",
                           rank, p, pos[p].x, pos[p].y, pos[p].z);
                 // Consider marking for removal here if this case should not happen:
                 // Maybe add to a separate 'removalList' or use the marking technique
                 // from the CheckAndRemove function if it's readily available.
            }
        } // end if isOutsideLocal
    } // end particle loop

    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
    // ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pids); CHKERRQ(ierr);

    *migrationList = currentMigrationList;
    *migrationCount = currentMigrationCount;
    *listCapacity = currentListCapacity;

    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Identified %d particles for potential migration.\n", rank, currentMigrationCount);

    PetscFunctionReturn(0);
}

/**
 * @brief Sets the target rank field (DMSwarmPICField_rank) for particles scheduled for migration.
 *
 * @param user           Pointer to UserCtx (contains swarm).
 * @param migrationList  Array of MigrationInfo structs containing local indices and target ranks.
 * @param migrationCount Number of particles in the migrationList.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode SetMigrationRanks(UserCtx* user, const MigrationInfo *migrationList, PetscInt migrationCount)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscInt       p_idx;
    PetscInt      *rankField = NULL; // Field storing target rank

    PetscFunctionBeginUser;

    // Ensure the migration rank field exists
    ierr = DMSwarmGetField(swarm, "DMSwarm_rank", NULL, NULL, (void **)&rankField); CHKERRQ(ierr);

    // Set the target rank for migrating particles
    for(p_idx = 0; p_idx < migrationCount; ++p_idx) {
        rankField[migrationList[p_idx].local_index] = migrationList[p_idx].target_rank;
    }

    ierr = DMSwarmRestoreField(swarm, "DMSwarm_rank", NULL, NULL, (void **)&rankField); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

/**
 * @brief Performs particle migration based on the pre-populated DMSwarmPICField_rank field.
 *
 * Assumes SetMigrationRanks has already been called to mark particles with their target ranks.
 * Calls DMSwarmMigrate to execute the communication and removal of un-migrated particles.
 *
 * @param user Pointer to the UserCtx structure containing the swarm.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformMigration(UserCtx *user)
{
    PetscErrorCode ierr;
    DM swarm = user->swarm;
    PetscMPIInt rank;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Starting DMSwarmMigrate...\n", rank);

    // Perform the migration - PETSC_TRUE removes particles that fail to land
    // in a valid cell on the target rank (or were marked with an invalid rank).
    ierr = DMSwarmMigrate(swarm, PETSC_TRUE); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Migration complete.\n", rank);
    PetscFunctionReturn(0);
}

//-----------------------------------------------------------------------------
// MODULE (COUNT): Calculates Particle Count Per Cell - REVISED FOR DMDAVecGetArray
//-----------------------------------------------------------------------------
/**
 * @brief Counts particles in each cell of the DMDA 'da' and stores the result in user->ParticleCount.
 * @ingroup scatter_module
 *
 * Zeros the user->ParticleCount vector, then iterates through local particles.
 * Reads the **GLOBAL** cell index (I, J, K) stored in the "DMSwarm_CellID" field.
 * Uses DMDAVecGetArray to access the local portion of the count vector and increments
 * the count at the global index (I, J, K) if it belongs to the local patch (including ghosts).
 *
 * @param[in,out] user Pointer to the UserCtx structure containing da, swarm, and ParticleCount.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.

 */
PetscErrorCode CalculateParticleCountPerCell(UserCtx *user) {
    PetscErrorCode ierr;
    DM             da = user->da;
    DM             swarm = user->swarm;
    Vec            countVec = user->ParticleCount;
    PetscInt       nlocal, p;
    PetscInt64       *global_cell_id_arr; // Read GLOBAL cell IDs
    PetscScalar    ***count_arr_3d;     // Use 3D accessor
    PetscInt64       *PID_arr;
    PetscMPIInt    rank;
    char           msg[ERROR_MSG_BUFFER_SIZE];
    PetscInt       particles_counted_locally = 0;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // --- Input Validation ---
    if (!da) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx->da is NULL.");
    if (!swarm) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx->swarm is NULL.");
    if (!countVec) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx->ParticleCount is NULL.");
    // Check DOF of da
    PetscInt count_dof;
    ierr = DMDAGetInfo(da, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &count_dof, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
     if (count_dof != 1) { PetscSNPrintf(msg, sizeof(msg), "countDM must have DOF=1, got %d.", count_dof); SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, msg); }

    // --- Zero the count vector ---
    ierr = VecSet(countVec, 0.0); CHKERRQ(ierr);

    // --- Get Particle Data ---
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "CalculateParticleCountPerCell: Accessing particle data.\n");
    ierr = DMSwarmGetLocalSize(swarm, &nlocal); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm,"DMSwarm_CellID", NULL, NULL, (void **)&global_cell_id_arr); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm,"DMSwarm_pid",NULL,NULL,(void **)&PID_arr);CHKERRQ(ierr);

    // --- Get Grid Vector Array using DMDA accessor ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "CalculateParticleCountPerCell: Accessing ParticleCount vector array (using DMDAVecGetArray).\n");
    ierr = DMDAVecGetArray(da, countVec, &count_arr_3d); CHKERRQ(ierr);

    // Get local owned range for validation/logging if needed, but not for indexing with DMDAVecGetArray
    PetscInt xs, ys, zs, xm, ym, zm;
    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);

    // --- Accumulate Counts Locally ---
    LOG_ALLOW(LOCAL, LOG_DEBUG, "CalculateParticleCountPerCell (Rank %d): Processing %d local particles using GLOBAL CellIDs.\n",rank,nlocal);
    for (p = 0; p < nlocal; p++) {
        // Read the GLOBAL indices stored for this particle
        PetscInt64 i = global_cell_id_arr[p * 3 + 0]; // Global i index
        PetscInt64 j = global_cell_id_arr[p * 3 + 1]; // Global j index
        PetscInt64 k = global_cell_id_arr[p * 3 + 2]; // Global k index

        // *** Bounds check is implicitly handled by DMDAVecGetArray for owned+ghost region ***
        // However, accessing outside this region using global indices WILL cause an error.
        // A preliminary check might still be wise if global IDs could be wild.
        // We rely on LocateAllParticles to provide valid global indices [0..IM-1] etc.
    // *** ADD PRINTF HERE ***
	LOG_LOOP_ALLOW(LOCAL,LOG_DEBUG,p,100,"[Rank %d CalcCount] Read CellID for p=%d, PID = %ld: (%ld, %ld, %ld)\n", rank, p,PID_arr[p],i, j, k);
    // *** END PRINTF ***
        // *** Access the local array using GLOBAL indices ***
        // DMDAVecGetArray allows this, mapping global (I,J,K) to the correct
        // location within the local ghosted array segment.
        // This assumes (I,J,K) corresponds to a cell within the local owned+ghost region.
        // If (I,J,K) is for a cell owned by *another* rank (and not in the ghost region
        // of this rank), accessing count_arr_3d[K][J][I] will likely lead to a crash
        // or incorrect results. This highlights why storing LOCAL indices is preferred
        // for parallel runs. But proceeding with GLOBAL as requested:
        if (i >= xs && i < xs + xm  && 
            j >= ys && j < ys + ym  && // Adjust based on actual ghost width
            k >= zs && k < zs + zm  )   // This check prevents definite crashes but doesn't guarantee ownership
        {
             // Increment count at the location corresponding to GLOBAL index (I,J,K)
	  //  LOG_ALLOW(LOCAL, LOG_DEBUG, "CalculateParticleCountPerCell (Rank %d): Particle %d with global CellID (%d, %d, %d) incremented with a particle.\n",rank, p, i, j, k);
             count_arr_3d[k][j][i] += 1.0;
             particles_counted_locally++;
         } else {
              // This particle's global ID is likely outside the range this rank handles (even ghosts)
              LOG_ALLOW(LOCAL, LOG_DEBUG, "CalculateParticleCountPerCell (Rank %d): Skipping particle %ld with global CellID (%ld, %ld, %ld) - likely outside local+ghost range.\n",rank, PID_arr[p] , i, j, k);
	}
    }
    LOG_ALLOW(LOCAL, LOG_DEBUG, "CalculateParticleCountPerCell (Rank %d): Local counting finished. Processed %d particles locally.\n", rank, particles_counted_locally);

    // --- Restore Access ---
    ierr = DMDAVecRestoreArray(da, countVec, &count_arr_3d); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm,"DMSwarm_CellID", NULL, NULL, (void **)&global_cell_id_arr); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm,"DMSwarm_pid",NULL,NULL,(void **)&PID_arr);CHKERRQ(ierr);

    // --- Assemble Global Vector ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "CalculateParticleCountPerCell: Assembling global ParticleCount vector.\n");
    ierr = VecAssemblyBegin(countVec); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(countVec); CHKERRQ(ierr);

    // --- Verification Logging ---
    PetscReal total_counted_particles = 0.0, max_count_in_cell = 0.0;
    ierr = VecSum(countVec, &total_counted_particles); CHKERRQ(ierr);
    PetscInt max_idx_global = -1;
    ierr = VecMax(countVec, &max_idx_global, &max_count_in_cell); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "CalculateParticleCountPerCell: Total counted globally = %.0f, Max count in cell = %.0f\n",
              total_counted_particles, max_count_in_cell);

    // --- ADD THIS DEBUGGING BLOCK ---
    if (max_idx_global >= 0) { // Check if VecMax found a location
         // Need to convert the flat global index back to 3D global index (I, J, K)
         // Get global grid dimensions (Nodes, NOT Cells IM/JM/KM)
         PetscInt M, N, P;
         ierr = DMDAGetInfo(da, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
         // Note: Assuming DOF=1 for countVec, index mapping uses node dimensions M,N,P from DMDA creation (IM+1, etc)
         // Re-check if your DMDA uses cell counts (IM) or node counts (IM+1) for Vec layout. Let's assume Node counts M,N,P.
         PetscInt Kmax = max_idx_global / (M * N);
         PetscInt Jmax = (max_idx_global % (M * N)) / M;
         PetscInt Imax = max_idx_global % M;
         LOG_ALLOW(GLOBAL, LOG_INFO, "  -> Max count located at global index (I,J,K) = (%d, %d, %d) [Flat index: %d]\n",
                   (int)Imax, (int)Jmax, (int)Kmax, (int)max_idx_global);

        // Also, let's explicitly check the count at (0,0,0)
        PetscScalar count_at_origin = 0.0;
        PetscScalar ***count_arr_for_check;
        ierr = DMDAVecGetArrayRead(da, countVec, &count_arr_for_check); CHKERRQ(ierr);
        // Check bounds before accessing - crucial if using global indices
        PetscInt xs, ys, zs, xm, ym, zm;
        ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);
        if (0 >= xs && 0 < xs+xm && 0 >= ys && 0 < ys+ym && 0 >= zs && 0 < zs+zm) {
             count_at_origin = count_arr_for_check[0][0][0]; // Access using global index (0,0,0)
        } else {
            // Origin is not on this rank (relevant for parallel, but check anyway)
            count_at_origin = -999.0; // Indicate it wasn't accessible locally
        }
        ierr = DMDAVecRestoreArrayRead(da, countVec, &count_arr_for_check); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_INFO, "  -> Count at global index (0,0,0) = %.1f\n", count_at_origin);

    } else {
         LOG_ALLOW(GLOBAL, LOG_WARNING, "  -> VecMax did not return a location for the maximum value.\n");
    }
    // --- END DEBUGGING BLOCK ---
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "CalculateParticleCountPerCell: Particle counting complete.\n");

    PetscFunctionReturn(0);
}

// NOTE: AccumulateParticleField also needs the same modification if it is to handle global CellIDs
// and use DMDAVecGetArray/DMDAVecGetArrayDOF

// --- Helper function to resize swarm globally (add or remove) ---
// This assumes removing excess particles means removing the globally last ones.
PetscErrorCode ResizeSwarmGlobally(DM swarm, PetscInt N_target)
{
    PetscErrorCode ierr;
    PetscInt       N_current, nlocal_current;
    PetscMPIInt    rank;
    MPI_Comm       comm;

    PetscFunctionBeginUser;
    ierr = PetscObjectGetComm((PetscObject)swarm, &comm); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);
    ierr = DMSwarmGetSize(swarm, &N_current); CHKERRQ(ierr);
    ierr = DMSwarmGetLocalSize(swarm, &nlocal_current); CHKERRQ(ierr);

    PetscInt delta = N_target - N_current;

    if (delta == 0) {
        PetscFunctionReturn(0); // Nothing to do
    }

    if (delta < 0) { // Remove particles
        PetscInt num_to_remove_global = -delta;
        LOG_ALLOW(GLOBAL, LOG_INFO, "ResizeSwarmGlobally: Current size %d > target size %d. Removing %d particles globally.\n", N_current, N_target, num_to_remove_global);

        // --- Strategy: Remove the globally last 'num_to_remove_global' particles ---
        // Each rank needs to determine how many of its *local* particles fall
        // within the range of global indices [N_target, N_current - 1].

        PetscInt rstart = 0;
	PetscInt  rend;
	// Global range owned by this rank [rstart, rend)

	ierr = MPI_Exscan(&nlocal_current, &rstart, 1, MPIU_INT, MPI_SUM, comm); CHKERRMPI(ierr); // Use CHKERRMPI for MPI calls

	rend = rstart + nlocal_current;
	

        // Calculate how many local particles have global indices >= N_target
        PetscInt nlocal_remove_count = 0;
        if (rend > N_target) { // If this rank owns any particles slated for removal
            PetscInt start_remove_local_idx = (N_target > rstart) ? (N_target - rstart) : 0;
            nlocal_remove_count = nlocal_current - start_remove_local_idx;
        }

        if (nlocal_remove_count < 0) nlocal_remove_count = 0; // Sanity check
        if (nlocal_remove_count > nlocal_current) nlocal_remove_count = nlocal_current; // Sanity check

        LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Global range [%d, %d). Target size %d. Need to remove %d local particles (from end).\n", rank, rstart, rend, N_target, nlocal_remove_count);

        // Remove the last 'nlocal_remove_count' particles *locally* by iterating backwards
        PetscInt removal_ops_done = 0;
        for (PetscInt p = nlocal_current - 1; p >= 0 && removal_ops_done < nlocal_remove_count; --p) {
            ierr = DMSwarmRemovePointAtIndex(swarm, p); CHKERRQ(ierr);
            removal_ops_done++;
        }

        if (removal_ops_done != nlocal_remove_count) {
             SETERRQ(comm, PETSC_ERR_PLIB, "Rank %d: Failed to remove the expected number of local particles (%d != %d)", rank, removal_ops_done, nlocal_remove_count);
        }
         LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Removed %d local particles.\n", rank, removal_ops_done);
 
	// Barrier to ensure all removals are done before size check
        ierr = MPI_Barrier(comm); CHKERRMPI(ierr);

    } else { // delta > 0: Add particles
        PetscInt num_to_add_global = delta;
        LOG_ALLOW(GLOBAL, LOG_INFO, "ResizeSwarmGlobally: Current size %d < target size %d. Adding %d particles globally.\n", N_current, N_target, num_to_add_global);
        ierr = DMSwarmAddNPoints(swarm, num_to_add_global); CHKERRQ(ierr);
        // Note: Added particles will have uninitialized field data. Reading will overwrite.
    }

    // Verify final size
    PetscInt N_final;
    ierr = DMSwarmGetSize(swarm, &N_final); CHKERRQ(ierr);
    if (N_final != N_target) {
        SETERRQ(comm, PETSC_ERR_PLIB, "Failed to resize swarm: expected %d particles, got %d", N_target, N_final);
    }
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ResizeSwarmGlobally: Swarm successfully resized to %d particles.\n", N_final);
    PetscFunctionReturn(0);
}

/**
 * @brief Checks particle count in the position file and resizes the swarm if needed (PETSc 3.18 compatible).
 *
 * Always uses "position" as the reference field and assumes a block size of 3.
 * Reads the position file (expected at results/position<#####>_0.ext format)
 * into a temporary Vec to determine the number of particles (`N_file`).
 * Compares `N_file` with the current swarm size (`N_current`). If they differ,
 * resizes the swarm globally.
 * If the position file is not found, returns PETSC_ERR_FILE_OPEN.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the DMSwarm.
 * @param[in]     ti   Time index for constructing the file name.
 * @param[in]     ext  File extension (e.g., "dat").
 *
 * @return PetscErrorCode 0 on success, non-zero on failure (including PETSC_ERR_FILE_OPEN).
 */
PetscErrorCode PreCheckAndResizeSwarm(UserCtx *user,
                                      PetscInt ti,
                                      const char *ext)
{
    PetscErrorCode ierr;
    char           filename[PETSC_MAX_PATH_LEN];
    PetscViewer    viewer;
    Vec            tmpVec = NULL;
    PetscInt       N_file = 0;
    PetscInt       N_current = 0;
    MPI_Comm       comm;
    PetscBool      fileExists = PETSC_FALSE;
    const char    *refFieldName = "position"; // Hardcoded reference field
    const PetscInt bs = 3;                   // Hardcoded block size for position
    const int      placeholder_int = 0;      // Placeholder integer for filename format

    PetscFunctionBeginUser;
    ierr = PetscObjectGetComm((PetscObject)user->swarm, &comm); CHKERRQ(ierr);

    // --- Construct filename using the specified format ---
    // results/%s%05<PetscInt_FMT>_%d.%s
    ierr = PetscSNPrintf(filename, sizeof(filename), "results/%s%05" PetscInt_FMT "_%d.%s",
                         refFieldName, ti, placeholder_int, ext); CHKERRQ(ierr);
    // Note: Make sure the "results" directory exists or handle directory creation elsewhere.

    // Log message indicating the check using the constructed filename
    LOG_ALLOW(GLOBAL, LOG_INFO, "PreCheckAndResizeSwarm: Checking particle count for timestep %d using reference file '%s' (bs=%d assumed).\n", ti, filename, bs);

    // --- Check existence of the position file ---
    ierr = PetscTestFile(filename, 'r', &fileExists); CHKERRQ(ierr);

    if (!fileExists) {
        SETERRQ(comm, PETSC_ERR_FILE_OPEN, "Mandatory reference file '%s' not found for timestep %d.", filename, ti);
    }

    // LOG_ALLOW(GLOBAL, LOG_DEBUG, "PreCheckAndResizeSwarm: Found position file %s. Determining particle count.\n", filename);

    // --- Load position file into temporary Vec to get size ---
    ierr = PetscViewerBinaryOpen(comm, filename, FILE_MODE_READ, &viewer); CHKERRQ(ierr);
    ierr = VecCreate(comm, &tmpVec); CHKERRQ(ierr);

    ierr = VecLoad(tmpVec, viewer);
    if (ierr) {
        PetscErrorCode ierr_destroy;
        ierr_destroy = PetscViewerDestroy(&viewer); CHKERRQ(ierr_destroy);
        ierr_destroy = VecDestroy(&tmpVec); CHKERRQ(ierr_destroy);
        // LOG_ALLOW(GLOBAL, LOG_ERROR, "PreCheckAndResizeSwarm: Failed to load reference file %s into temporary vector. Error %d\n", filename, ierr);
        PetscFunctionReturn(ierr); // Return the VecLoad error
    }

    // --- Calculate particle count using hardcoded block size ---
    PetscInt vecSize;
    ierr = VecGetSize(tmpVec, &vecSize); CHKERRQ(ierr);
    if (vecSize % bs != 0) {
        PetscErrorCode ierr_destroy;
        ierr_destroy = PetscViewerDestroy(&viewer); CHKERRQ(ierr_destroy);
        ierr_destroy = VecDestroy(&tmpVec); CHKERRQ(ierr_destroy);
        SETERRQ(comm, PETSC_ERR_FILE_READ, "Temporary vector size %d from file %s is not divisible by assumed block size %d for '%s'", vecSize, filename, bs, refFieldName);
    }
    N_file = vecSize / bs; // Calculate particle count

    // LOG_ALLOW(GLOBAL, LOG_DEBUG, "PreCheckAndResizeSwarm: File %s contains %d particles (Vec size %d, BS %d).\n", filename, N_file, vecSize, bs);

    // --- Clean up temporary objects ---
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    ierr = VecDestroy(&tmpVec); CHKERRQ(ierr);

    // --- Compare with current swarm size and resize if needed ---
    ierr = DMSwarmGetSize(user->swarm, &N_current); CHKERRQ(ierr);

    if (N_file == N_current) {
        // LOG_ALLOW(GLOBAL, LOG_DEBUG, "PreCheckAndResizeSwarm: Swarm size (%d) already matches file size (%d). No resize needed for timestep %d.\n", N_current, N_file, ti);
    } else {
        // LOG_ALLOW(GLOBAL, LOG_INFO, "PreCheckAndResizeSwarm: Swarm size %d differs from file size %d for timestep %d. Resizing swarm.\n", N_current, N_file, ti);
        ierr = ResizeSwarmGlobally(user->swarm, N_file); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0); // Return 0 only if everything succeeded
}
