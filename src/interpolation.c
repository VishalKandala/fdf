/**
 * @file interpolation.c  //  Particle In Cell main.
 * @brief Main program for DMSwarm interpolation using the fdf-curvIB method.
 *
 * Initializes a particle swarm, reads velocity fields, and performs particle-grid interpolation.
 */

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
#include "interpolation.h"


#define NUM_WEIGHTS 8 // Number of weights in trilinear interpolation

/**
 * @brief Computes the trilinear interpolation weights from the interpolation coefficients.
 *
 * This function computes the weights for trilinear interpolation at the eight corners of a cell
 * using the interpolation coefficients provided along the x, y, and z directions.
 *
 * @param[in]  a1 Interpolation coefficient along the x-direction (normalized coordinate within the cell).
 * @param[in]  a2 Interpolation coefficient along the y-direction (normalized coordinate within the cell).
 * @param[in]  a3 Interpolation coefficient along the z-direction (normalized coordinate within the cell).
 * @param[out] w  Array of 8 weights, each corresponding to one corner of the cell.
 *
 * @note
 * - The coefficients `a1`, `a2`, and `a3` should be in the range [0, 1].
 * - The order of weights corresponds to the eight corners of a hexahedral cell.
 */
static inline void ComputeTrilinearWeights(PetscReal a1, PetscReal a2, PetscReal a3, PetscReal *w) {
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ComputeTrilinearWeights: Computing weights for a1=%f, a2=%f, a3=%f.\n", a1, a2, a3);

    // Compute complementary coefficients (1 - a)
    PetscReal oa1 = a1 - 1.0; 
    PetscReal oa2 = a2 - 1.0;
    PetscReal oa3 = a3 - 1.0;

    // Compute weights for each corner of the cell
    w[0] = a1  * a2  * a3;   // Front-bottom-left
    w[1] = a1  * oa2 * oa3;  // Front-top-right
    w[2] = oa1 * a2  * oa3;  // Back-bottom-right
    w[3] = a1  * a2  * oa3;  // Front-bottom-right
    w[4] = oa1 * oa2 * a3;   // Back-top-left
    w[5] = oa1 * oa2 * oa3;  // Back-top-right
    w[6] = oa1 * a2  * a3;   // Back-bottom-left
    w[7] = a1  * oa2 * a3;   // Front-top-left

    // Log the computed weights for debugging
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ComputeTrilinearWeights: Weights computed - "
        "w0=%f, w1=%f, w2=%f, w3=%f, w4=%f, w5=%f, w6=%f, w7=%f. \n",
        w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7]);
}

/**
 * @brief Performs trilinear interpolation of velocities from grid to particles using interpolation coefficients.
 *
 * This function interpolates velocities for particles based on the velocity field defined on the computational grid.
 * It retrieves cell indices for each particle from the "DMSwarm_CellID" field and assigns the interpolated velocity
 * using trilinear interpolation from the surrounding grid cells. The function handles boundary checks and ensures
 * that particles outside the valid grid range are appropriately managed.
 *
 * Key Steps:
 * 1. Retrieve the number of local particles and their associated data fields (cell indices and velocities).
 * 2. Map the global velocity field (`Ucat`) to the local portion for efficient access.
 * 3. Retrieve interpolation coefficients (`a1`, `a2`, `a3`) for each particle.
 * 4. Compute trilinear interpolation weights and interpolate velocities from the surrounding grid cells.
 * 5. Handle edge cases where particles may lie outside the valid grid range.
 * 6. Restore all data fields and ensure PETSc arrays and vectors are correctly finalized.
 *
 * @param[in] user Pointer to the `UserCtx` structure containing:
 *                 - `user->da`: DMDA for the grid.
 *                 - `user->swarm`: DMSwarm for particles.
 *                 - `user->Ucat`: Global velocity vector on the grid.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InterpolateParticleVelocities(UserCtx *user) {
    PetscErrorCode ierr;
    PetscInt n_local;             // Number of local particles
    PetscInt64 *cellIDs = NULL;     // Array to store cell indices for each particle
    PetscReal *velocities = NULL; // Array to store particle velocities
    PetscReal *weights = NULL;    // Array to store interpolation weights
    Cmpnts ***ucat;               // 3D array to map local grid velocities
    PetscInt i,j,k;
    DM fda = user->fda;           // Field DA 
    DM swarm = user->swarm;       // DMSwarm for the particles


    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Starting particle velocity interpolation.\n");

    // Verify global velocity field
    PetscReal max_val;
    ierr = VecNorm(user->Ucat, NORM_INFINITY, &max_val); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Global velocity field Ucat maximum magnitude: %f \n", max_val);

    // Retrieve the number of local particles
    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Found %d local particles.\n", n_local);

    // Retrieve particle data fields
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr); // Ensure 'weight' field exists

    // Access the local portion of the global velocity vector (Ucat) using 'fda'
    ierr = DMDAVecGetArrayRead(user->fda,user->Ucat,&ucat);


    LOG_ALLOW(GLOBAL, LOG_DEBUG, "InterpolateParticleVelocities: Starting velocity assignment for particles.\n");

    // Log grid dimensions
    LOG_ALLOW(GLOBAL, LOG_INFO, "Grid dimensions: mx=%d, my=%d, mz=%d \n", user->info.mx, user->info.my, user->info.mz);

    // Loop over all local particles
    for (PetscInt p = 0; p < n_local; p++) {
        // Retrieve cell indices for the particle
        i = cellIDs[3 * p];
        j = cellIDs[3 * p + 1];
        k = cellIDs[3 * p + 2];

        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Particle %d: Host Cell = (%d, %d, %d)\n", p, i, j, k);

        // Validate cell indices (boundary check)
        if (i < 0 || j < 0 || k < 0 || i >= user->info.mx || j >= user->info.my || k >= user->info.mz) {
            LOG_ALLOW(GLOBAL, LOG_WARNING, "Particle %d has invalid cell indices (%d, %d, %d)\n. Skipping interpolation.\n", p, i, j, k);
            velocities[3 * p    ] = 0.0;
            velocities[3 * p + 1] = 0.0;
            velocities[3 * p + 2] = 0.0;
            continue;
        }

        // Assign interpolated velocity to the particle

        velocities[3 * p]     = ucat[k][j][i].x; // u-component
        velocities[3 * p + 1] = ucat[k][j][i].y; // v-component
        velocities[3 * p + 2] = ucat[k][j][i].z; // w-component

        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Particle %d: Interpolated velocity from (x=%f, y=%f, z=%f) to (x=%f, y=%f, z=%f).\n",
            p, ucat[k][j][i].x, ucat[k][j][i].y, ucat[k][j][i].z,velocities[3 * p], velocities[3 * p + 1], velocities[3 * p + 2]);
    }

    // Restore the local velocity array
    ierr = DMDAVecRestoreArrayRead(fda,user->Ucat, &ucat); CHKERRQ(ierr);

    // Restore particle data fields
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Particle velocity interpolation completed.\n");
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

    LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeSimulation Runtime Options:\n");
    LOG_ALLOW(GLOBAL, LOG_INFO, "rstart = %d\n", *rstart);
    if((*rstart) == 1) LOG_ALLOW(GLOBAL, LOG_INFO, "Restarting from time: %d\n", *ti);
    LOG_ALLOW(GLOBAL, LOG_INFO, "No. of Particles: %d\n", *np);
    LOG_ALLOW(GLOBAL, LOG_INFO, "No. of Grid blocks: %d\n", *nblk);

    return 0;
}


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

    LOG_ALLOW(GLOBAL, LOG_INFO, "SetupGridAndVectors - Grid and vectors setup completed.\n");
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
PetscErrorCode PerformParticleSwarmOperations(UserCtx *user, PetscInt np, BoundingBox *bboxlist) {
    PetscErrorCode ierr;
    PetscInt particlesPerProcess = 0;         // Number of particles assigned to the loca  l MPI process.
    PetscRandom randx,randy,randz;     // Random number generators[x,y,z]. (used if ParticleInitialization==1).       

    LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Starting particle swarm operations with %d particles.\n", np);

    // Step 1: Create and initialize the particle swarm
    // Here we pass in the bboxlist, which will be used by CreateParticleSwarm() and subsequently
    // by AssignInitialProperties() to position particles if ParticleInitialization == 0.
    LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Initializing particle swarm.\n");
    ierr = CreateParticleSwarm(user, np, &particlesPerProcess,bboxlist); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Particle swarm initialized successfully.\n");

    // Assign initial properties to particles
    // The bboxlist array is passed here so that if ParticleInitialization == 0,
    // particles can be placed at the midpoint of the local bounding box corresponding to this rank.
    ierr = AssignInitialPropertiesToSwarm(user, particlesPerProcess, &randx, &randy, &randz, bboxlist); CHKERRQ(ierr);

    // Finalize swarm setup by destroying RNGs if ParticleInitialization == 1
    ierr = FinalizeSwarmSetup(&randx, &randy, &randz); CHKERRQ(ierr);

    // Step 2: Print particle fields for debugging (optional)
 //   LOG_ALLOW(GLOBAL, LOG_DEBUG, "PerformParticleSwarmOperations - Printing initial particle fields (optional).\n");
  //  ierr = PrintParticleFields(user); CHKERRQ(ierr);

    // Step 3: Locate particles within the computational grid
    // This updates each particle's cell indices and interpolation weights as needed.
    LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Locating particles within the grid.\n");
    ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Particle positions updated after grid search.\n");

    // Print particle fields again after location update (optional)
  //  LOG_ALLOW(GLOBAL, LOG_DEBUG, "PerformParticleSwarmOperations - Printing particle fields after location update.\n");
 //   ierr = PrintParticleFields(user); CHKERRQ(ierr);

    // Step 4: Interpolate particle velocities using trilinear interpolation
    // This requires that particles have valid cell indices and weights from the previous step.
    LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Interpolating particle velocities using trilinear interpolation.\n");
    ierr = InterpolateParticleVelocities(user); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Particle velocities interpolated successfully.\n");

    // Print particle fields again after velocity interpolation (optional)
   LOG_ALLOW(GLOBAL, LOG_DEBUG, "PerformParticleSwarmOperations - Printing particle fields after velocity interpolation.\n");
   ierr = PrintParticleFields(user); CHKERRQ(ierr);

   // Write the particle positions to file.
   LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Writing the particle positions to file.\n");
   ierr = WriteSwarmField(user,"position",0,"dat");

   // Write the interpolated velocity to file.
   LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Writing the interpolated velocities to file.\n");
   ierr = WriteSwarmField(user,"velocity",0,"dat");

    LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Particle swarm operations completed.\n");

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "Simulation finalized and resources cleaned up.\n");

    // Finalize PETSc
    ierr = PetscFinalize(); CHKERRQ(ierr);
    return 0;
}
