/**
 * @file setup.c  //  Setup code for running any simulation 
 * @brief Test program for DMSwarm interpolation using the fdf-curvIB method.
 * Provides the setup to start any simulation with DMSwarm and DMDAs.
 **/

 #include <petscpf.h>
 #include <petscdmswarm.h>
 #include <stdlib.h>
 #include <time.h>
 #include <math.h>
 #include <petsctime.h>
 #include <petscdmcomposite.h>
 
 // Include the updated headers
 //#include "common.h"         // Shared type definitions
 //#include "ParticleSwarm.h"  // Particle swarm functions
 //#include "walkingsearch.h"  // Particle location functions
 //#include "grid.h"           // Grid functions
 //#include "logging.h"        // Logging macros
 //#include "io.h"             // Data Input and Output functions

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
 * @brief Initialize the simulation context.
 *
 * Checks for the presence of "control.dat" file, reads runtime options, and sets up the user context.
 *
 * @param[out] user    Pointer to the allocated UserCtx structure.
 * @param[out] rank    MPI rank of the process.
 * @param[out] size    Number of MPI processes.
 * @param[out] np      Number of particles.
 * @param[out] rstart  Flag to restart (1) or start from t=0 (0).
 * @param[out] ti      The timestep to start from if restarting.
 * @param[out] nblk    Number of grid blocks.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
 PetscErrorCode InitializeSimulation(UserCtx **user, PetscInt *rank, PetscInt *size, PetscInt *np, PetscInt *rstart, PetscInt *ti, PetscInt *nblk) {
    PetscErrorCode ierr;

    // Attempt to insert options from "control.dat"
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, "control.dat", PETSC_TRUE);
    if (ierr == PETSC_ERR_FILE_OPEN) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
                "InitializeSimulation - Could not open 'control.dat'. Please ensure it exists in the current directory.");
    } else {
        CHKERRQ(ierr);
    }

    // Allocate user context
    ierr = PetscCalloc1(1, user); CHKERRQ(ierr);

    // Get MPI rank and size
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, size); CHKERRQ(ierr);

    // Initialize user context flags
    (*user)->averaging = PETSC_FALSE;
    (*user)->les = PETSC_FALSE;
    (*user)->rans = PETSC_FALSE;

    LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeSimulation - Initialized on rank %d out of %d processes.\n", *rank, *size);

    // Read runtime options
    ierr = PetscOptionsGetInt(NULL, NULL, "-numParticles", np, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-rstart", rstart, NULL); CHKERRQ(ierr);
    if((*rstart) == 1) {
        ierr = PetscOptionsGetInt(NULL, NULL, "-ti", ti, NULL); CHKERRQ(ierr);
    }
    ierr = PetscOptionsGetInt(NULL, NULL, "-nblk", nblk, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-pinit", &((*user)->ParticleInitialization), NULL);
    ierr = PetscOptionsGetReal(NULL, NULL, "-dt", &((*user)->dt), NULL);

    LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeSimulation Runtime Options:\n");
    LOG_ALLOW(GLOBAL, LOG_INFO, "rstart = %d\n", *rstart);
    if((*rstart) == 1) LOG_ALLOW(GLOBAL, LOG_INFO, "Restarting from time: %d\n", *ti);
    LOG_ALLOW(GLOBAL, LOG_INFO, "No. of Particles: %d\n", *np);
    LOG_ALLOW(GLOBAL, LOG_INFO, "No. of Grid blocks: %d\n", *nblk);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Time-step size : %f\n", (*user)->dt);

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "InitializeSimulation: Completed setup across all ranks.\n");


    return 0;
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
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode FinalizeSimulation(UserCtx *user, PetscInt block_number, BoundingBox *bboxlist) {
    PetscErrorCode ierr;
    PetscMPIInt rank;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr); // Corrected: use &rank

    // Destroy DM and vectors for each block
    for (PetscInt bi = 0; bi < block_number; bi++) {
        ierr = VecDestroy(&(user[bi].Ucat)); CHKERRQ(ierr);
        ierr = VecDestroy(&(user[bi].Ucont)); CHKERRQ(ierr);
        ierr = VecDestroy(&(user[bi].P)); CHKERRQ(ierr);
        ierr = VecDestroy(&(user[bi].Nvert)); CHKERRQ(ierr);
        ierr = VecDestroy(&(user[bi].Nvert_o)); CHKERRQ(ierr);

        ierr = DMDestroy(&(user[bi].fda)); CHKERRQ(ierr);
        ierr = DMDestroy(&(user[bi].da)); CHKERRQ(ierr);
        ierr = DMDestroy(&(user[bi].swarm)); CHKERRQ(ierr);
    }

    // Free user context
    ierr = PetscFree(user); CHKERRQ(ierr);

    // Now free bboxlist on all ranks since all allocated their own copy
    if (bboxlist) {
        free(bboxlist);
        bboxlist = NULL;
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "Simulation finalized and resources cleaned up on all ranks. \n");

    // Ensure all MPI ranks reach this point before finalizing PETSc
    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

    // Finalize PETSc
    ierr = PetscFinalize(); CHKERRQ(ierr);
    return 0;
}