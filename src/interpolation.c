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
#include "common.h"         // Shared type definitions
#include "ParticleSwarm.h"  // Particle swarm functions
#include "walkingsearch.h"  // Particle location functions
#include "grid.h"           // Grid functions
#include "logging.h"        // Logging macros
#include "io.h"             // Data Input and Output functions

static char help[] = "DMSwarm Interpolation - fdf-curvIB ";

/**
 * @brief Initialize PETSc and the simulation context.
 *
 * @param[out] user Pointer to the allocated UserCtx structure.
 * @param[out] rank MPI rank of the process.
 * @param[out] size Number of MPI processes.
 * @param[out] np Number of particles.
 * @param[out] rstart Flag to restart(1) or start from t = 0 (0).
 * @param[out] ti The timestep to start from if restarting.
 * @param[out] nblk Number of grid blocks.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeSimulation(UserCtx **user, PetscInt *rank, PetscInt *size, PetscInt *np, PetscInt *rstart, PetscInt *ti, PetscInt *nblk) {
    PetscErrorCode ierr;

    // Initialize PETSc
    ierr = PetscInitialize(NULL, NULL, NULL, help); CHKERRQ(ierr);

    // Read options from control.dat
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, "control.dat", PETSC_TRUE); CHKERRQ(ierr);

    // Allocate user context
    ierr = PetscCalloc1(1, user); CHKERRQ(ierr);

    // Get MPI rank and size
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, size); CHKERRQ(ierr);

    // Initialize user context
    (*user)->averaging = PETSC_FALSE;
    (*user)->les = PETSC_FALSE;
    (*user)->rans = PETSC_FALSE;

    LOG(GLOBAL, LOG_INFO, "InitializeSimulation - Initialized on rank %d out of %d processes.\n", *rank, *size);

   // Read runtime options
    ierr = PetscOptionsGetInt(NULL, NULL, "-numParticles", np, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-rstart", rstart, NULL); CHKERRQ(ierr);
    if((*rstart)==1){
      ierr = PetscOptionsGetInt(NULL,NULL,"-ti",ti,NULL); CHKERRQ(ierr);
    }
    ierr = PetscOptionsGetInt(NULL, NULL, "-nblk",nblk, NULL); CHKERRQ(ierr);

    LOG(GLOBAL, LOG_INFO, "InitializeSimulation Runtime Options: \n");
    LOG(GLOBAL, LOG_INFO, "rstart = %d \n",*rstart);
    if(rstart) LOG(GLOBAL, LOG_INFO, "Restarting from time: %d \n",*ti);
    LOG(GLOBAL, LOG_INFO, "No.of Particles: %d \n", *np);
    LOG(GLOBAL, LOG_INFO, "No.of Grid blocks: %d\n",*nblk);

    return 0;
}

/**
 * @brief Set up the simulation grid and vectors.
 *
 * @param[in,out] user Pointer to the UserCtx structure.
 * @param[in] block_number The number of grid blocks in the domain.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode SetupGridAndVectors(UserCtx *user,PetscInt block_number) {
    PetscErrorCode ierr;

    // Define grid coordinates
    ierr = DefineGridCoordinates(user); CHKERRQ(ierr);

    // Create global vectors for each block
    for (PetscInt bi = 0; bi < block_number; bi++) {
        ierr = DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat); CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(user[bi].fda, &user[bi].Ucont); CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(user[bi].da, &user[bi].P); CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(user[bi].da, &user[bi].Nvert); CHKERRQ(ierr);
        ierr = DMCreateGlobalVector(user[bi].da, &user[bi].Nvert_o); CHKERRQ(ierr);
    }

    LOG(GLOBAL, LOG_INFO, "Grid and vectors setup completed.\n");
    return 0;
}

/**
 * @brief Perform particle swarm initialization, particle-grid interaction, and related operations.
 *
 * @param[in,out] user Pointer to the UserCtx structure.
 * @param[in] np Number of particles.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode PerformParticleSwarmOperations(UserCtx *user, PetscInt np) {
    PetscErrorCode ierr;

    // Create and initialize particle swarm
    ierr = CreateParticleSwarm(user, np); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_INFO, "Particle swarm initialized.\n");

    // Locate particles in the grid
    ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_INFO, "Particle positions updated after search.\n");

    // Print particle fields (optional debugging)
    ierr = PrintParticleFields(user); CHKERRQ(ierr);

    return 0;
}

/**
 * @brief Finalize the simulation and free resources.
 *
 * @param[in,out] user Pointer to the UserCtx structure.
 * @param[in] block_number The number of grid blocks in the domain.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode FinalizeSimulation(UserCtx *user, PetscInt block_number) {
    PetscErrorCode ierr;

    // Destroy DM and vectors for each block
    for (PetscInt bi = 0; bi < block_number; bi++) {
        ierr = DMDestroy(&(user[bi].fda)); CHKERRQ(ierr);
        ierr = DMDestroy(&(user[bi].da)); CHKERRQ(ierr);
        ierr = DMDestroy(&(user[bi].swarm)); CHKERRQ(ierr);
    }

    // Free user context
    ierr = PetscFree(user); CHKERRQ(ierr);

    LOG(GLOBAL, LOG_INFO, "Simulation finalized and resources cleaned up.\n");

    // Finalize PETSc
    ierr = PetscFinalize(); CHKERRQ(ierr);
    return 0;
}

/**
 * @brief Main function for DMSwarm interpolation.
 *
 * Initializes the grid, reads the velocity field, creates particles, and performs particle-grid interpolation.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return int Returns 0 on success, non-zero on failure.
 */
int main(int argc, char **argv) {
    PetscErrorCode ierr;
    UserCtx *user = NULL;
    PetscInt rank, size;
    PetscInt np = 0,ti = 0;
    PetscReal umax;
    PetscInt block_number = 1, rstart = 0;
    BoundingBox *bboxlist;

    // Initialize simulation
    ierr = InitializeSimulation(&user, &rank, &size,&np,&rstart,&ti,&block_number); CHKERRQ(ierr);
    // Setup grid and vectors
    ierr = SetupGridAndVectors(user,block_number); CHKERRQ(ierr);

    // Read simulation fields
    ierr = ReadSimulationFields(user, ti); CHKERRQ(ierr);

    // Compute and log max velocity
    ierr = VecNorm(user->Ucat, NORM_INFINITY, &umax); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_INFO, "Maximum velocity magnitude: %f\n", umax);

    // Create Bounding Boxes for each rank
    ierr = GatherAllBoundingBoxes(user, &bboxlist); CHKERRQ(ierr);

    // Perform particle swarm operations
    ierr = PerformParticleSwarmOperations(user, np); CHKERRQ(ierr);

    // Finalize simulation
    ierr = FinalizeSimulation(user,block_number); CHKERRQ(ierr);

    return 0;
}
