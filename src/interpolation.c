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


#define NUM_WEIGHTS 8 // Number of weights in trilinear interpolation




/**
 * @brief Computes the trilinear interpolation weights from the interpolation coefficients.
 *
 * @param[in]  a1 Interpolation coefficient along the x-direction.
 * @param[in]  a2 Interpolation coefficient along the y-direction.
 * @param[in]  a3 Interpolation coefficient along the z-direction.
 * @param[out] w  Array of 8 weights to be computed.
 */
static inline void ComputeTrilinearWeights(PetscReal a1, PetscReal a2, PetscReal a3, PetscReal *w)
{
    PetscReal oa1 = 1.0 - a1;
    PetscReal oa2 = 1.0 - a2;
    PetscReal oa3 = 1.0 - a3;

    // Compute weights for each corner of the cell
    w[0] = oa1 * oa2 * oa3;
    w[1] = a1  * oa2 * oa3;
    w[2] = oa1 * a2  * oa3;
    w[3] = a1  * a2  * oa3;
    w[4] = oa1 * oa2 * a3;
    w[5] = a1  * oa2 * a3;
    w[6] = oa1 * a2  * a3;
    w[7] = a1  * a2  * a3;
}

/**
 * @brief Performs trilinear interpolation of velocities from grid to particles using interpolation coefficients.
 *
 * This function interpolates the velocities at particle positions using trilinear interpolation.
 * It extracts the interpolation coefficients `a1`, `a2`, and `a3` from the `"weights"` field
 * of each particle and computes the interpolation weights on-the-fly.
 *
 * @param[in] user Pointer to the user-defined context containing grid and swarm information.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InterpolateParticleVelocities(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscInt n_local;
    PetscReal *positions = NULL;
    PetscReal *weights = NULL;      // Contains a1, a2, a3 for each particle
    PetscInt *cellIDs;
    PetscReal *velocities = NULL;

    Vec Ucat_local;
    PetscReal ****u_array = NULL; // 4D array: [z][y][x][dof]

    DM da = user->da;          // DMDA for the grid
    DM swarm = user->swarm;    // DMSwarm for the particles

    PetscFunctionBeginUser;

    // Access DMSwarm fields
    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);

    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);

    // Access the local portion of Ucat
    ierr = DMGetLocalVector(da, &Ucat_local); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(da, user->Ucat, INSERT_VALUES, Ucat_local); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, user->Ucat, INSERT_VALUES, Ucat_local); CHKERRQ(ierr);

    // Access grid velocities
    ierr = DMDAVecGetArrayDOFRead(da, Ucat_local, &u_array); CHKERRQ(ierr);

    // Loop over local particles
    for (PetscInt p = 0; p < n_local; ++p) {
        // Retrieve cell indices for the particle
        PetscInt i = cellIDs[3 * p];
        PetscInt j = cellIDs[3 * p + 1];
        PetscInt k = cellIDs[3 * p + 2];

        // Retrieve interpolation coefficients from the "weights" field
        PetscReal a1 = weights[3 * p];       // a1 for particle p
        PetscReal a2 = weights[3 * p + 1];   // a2 for particle p
        PetscReal a3 = weights[3 * p + 2];   // a3 for particle p

        // Compute interpolation weights from coefficients
        PetscReal w[NUM_WEIGHTS];
        ComputeTrilinearWeights(a1, a2, a3, w);

        // Initialize interpolated velocity components
        PetscReal vel_p[3] = {0.0, 0.0, 0.0};

        // Loop over the 8 corners of the cell
        for (PetscInt corner = 0; corner < NUM_WEIGHTS; ++corner) {
            // Determine offsets based on corner index
            PetscInt di = (corner & 1) ? 1 : 0;
            PetscInt dj = (corner & 2) ? 1 : 0;
            PetscInt dk = (corner & 4) ? 1 : 0;

            PetscInt ii = i + di;
            PetscInt jj = j + dj;
            PetscInt kk = k + dk;

            // Access grid velocities at the grid point (ii, jj, kk)
            PetscReal *grid_vel = u_array[kk][jj][ii]; // u_array[z][y][x][dof]

            // Accumulate weighted velocities
            vel_p[0] += w[corner] * grid_vel[0]; // u-component
            vel_p[1] += w[corner] * grid_vel[1]; // v-component
            vel_p[2] += w[corner] * grid_vel[2]; // w-component
        }

        // Store interpolated velocities back into DMSwarm field
        velocities[3 * p]     = vel_p[0];
        velocities[3 * p + 1] = vel_p[1];
        velocities[3 * p + 2] = vel_p[2];
    }

    // Restore grid velocities
    ierr = DMDAVecRestoreArrayDOFRead(da, Ucat_local, &u_array); CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(da, &Ucat_local); CHKERRQ(ierr);

    // Restore DMSwarm fields
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////

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

    PetscPrintf(PETSC_COMM_WORLD," Before Locating Host \n");

    // Print particle fields (optional debugging)
    ierr = PrintParticleFields(user); CHKERRQ(ierr);

    // Locate particles in the grid
    ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_INFO, "Particle positions updated after search.\n");

    PetscPrintf(PETSC_COMM_WORLD," After Locating Host \n");

    // Print particle fields (optional debugging)
    ierr = PrintParticleFields(user); CHKERRQ(ierr);

    // Interpolate particle velocity from grid points (Trilinear)
    ierr = InterpolateParticleVelocities(user); CHKERRQ(ierr);

   PetscPrintf(PETSC_COMM_WORLD," After Interpolating Velocity \n");

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
