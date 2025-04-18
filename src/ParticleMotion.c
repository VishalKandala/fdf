// ParticleMotion.c

#include "ParticleMotion.h"

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

  PetscFunctionBeginUser;  // PETSc macro for error/stack tracing

  // 1) Get the number of local particles
  ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);

  // 2) Access the "position" and "velocity" fields
  ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&pos); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&vel); CHKERRQ(ierr);

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
		    
                    ierr = PetscRealloc(newCapacity, &currentMigrationList); CHKERRQ(ierr);
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
