/**
 * @file inttest.c  //  Particle In Cell main.
 * @brief Test program for DMSwarm interpolation using the fdf-curvIB method.
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
#include "setup.h"
#include "AnalyticalSolution.h"
#include "ParticleMotion.h"

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Starting particle swarm operations with %ld particles.\n", np);

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
  //  ierr = LOG_PARTICLE_FIELDS(user,100); CHKERRQ(ierr);

    // Step 3: Locate particles within the computational grid
    // This updates each particle's cell indices and interpolation weights as needed.
    LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Locating particles within the grid.\n");
    ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Particle positions updated after grid search.\n");

    // Interpolate particle velocities using trilinear interpolation
    // This requires that particles have valid cell indices and weights from the previous step.
    LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Interpolating particle velocities using trilinear interpolation.\n");
    //  ierr = InterpolateParticleVelocities(user); CHKERRQ(ierr);
    ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Particle velocities interpolated successfully.\n");

    // Print particle fields again after velocity interpolation (optional)
   LOG_ALLOW(GLOBAL, LOG_DEBUG, "PerformParticleSwarmOperations - Printing particle fields after velocity interpolation.\n");
    ierr = LOG_PARTICLE_FIELDS(user,10); CHKERRQ(ierr);

    // Print out the interpolation error
    ierr = LOG_INTERPOLATION_ERROR(user); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Particle positions updated after velocity interpolation.\n");

   // Write the particle positions to file.
    //  LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Writing the particle positions to file.\n");
    //  ierr = WriteSwarmField(user,"position",0,"dat");

   // Write the interpolated velocity to file.
   // LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Writing the interpolated velocities to file.\n");
   // ierr = WriteSwarmField(user,"velocity",0,"dat");

    // Update the particle positions based on the interpolated velocities.
    //  ierr = UpdateAllParticlePositions(user); CHKERRQ(ierr);
    //  LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Particle positions updated after velocity interpolation.\n");
   
    // Step 3: Locate particles within the computational grid
    // This updates each particle's cell indices and interpolation weights as needed.
    //  LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Locating particles within the grid.\n");
    //  ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);
    //   LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Particle positions updated after grid search.\n");

    // Interpolate particle velocities using trilinear interpolation
    // This requires that particles have valid cell indices and weights from the previous step.
    //  LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Interpolating particle velocities using trilinear interpolation.\n");
    //  ierr = InterpolateParticleVelocities(user); CHKERRQ(ierr);
    //  LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Particle velocities interpolated successfully.\n");

    // Print particle fields again after velocity interpolation (optional)
    //  LOG_ALLOW(GLOBAL, LOG_DEBUG, "PerformParticleSwarmOperations - Printing particle fields after velocity interpolation.\n");
    //   ierr = LOG_PARTICLE_FIELDS(user,10); CHKERRQ(ierr);

    // Print out the interpolation error
    // ierr = LOG_INTERPOLATION_ERROR(user); CHKERRQ(ierr);

    // Write the particle positions to file.
    //  LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Writing the particle positions to file.\n");
    //  ierr = WriteSwarmField(user,"position",1,"dat");

   // Write the interpolated velocity to file.
    //  LOG_ALLOW(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Writing the interpolated velocities to file.\n");
    // ierr = WriteSwarmField(user,"velocity",1,"dat");   

    // Ensure all ranks complete before proceeding
    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "PerformParticleSwarmOperations - Particle data writing completed across all ranks.\n");

    return 0;
}

#undef _FUNCT_
#define __FUNCT__ "main"

/**
 * @brief Main function for DMSwarm interpolation.
 *
 * Initializes the grid, reads/writes the velocity field depending on command line options,
 * creates particles, and performs particle-grid interpolation.
 *
 * Introduces a command-line option `-read_fields` to toggle between:
 * - Updating and writing the Cartesian velocity fields (default)
 * - Reading the Cartesian velocity fields from file
 *
 * Additionally, after gathering the bounding boxes on rank 0, it now calls a separate function,
 * `BroadcastAllBoundingBoxes()`, to distribute bboxlist to all ranks.
 *
 * Usage example:
 * - To run with update and write steps:
 *   `mpirun -np X ./your_executable`
 * - To run with reading fields:
 *   `mpirun -np X ./your_executable -read_fields`
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 *
 * @return int Returns 0 on success, non-zero on failure.
 */
int main(int argc, char **argv) {
    UserCtx *user = NULL;     // User context
    PetscErrorCode ierr;      // PETSc error handling
    PetscInt block_number = 1;
    PetscInt rstart = 0;
    PetscInt np = 0;
    PetscInt ti = 0;
    PetscMPIInt rank, size;
    BoundingBox *bboxlist;    // Array of bounding boxes
    PetscReal umax;
    PetscBool readFields = PETSC_FALSE;
    static char help[] = " Test for interpolation - swarm-curvIB";
    PetscViewer logviewer;

    // -------------------- 1. PETSc Initialization --------------------
    ierr = PetscInitialize(&argc, &argv, (char *)0, help); CHKERRQ(ierr);

    // -------------------- 2. Setup Logging Allow-List ----------------
    // Only these function names will produce LOG_ALLOW (or LOG_ALLOW_SYNC) output.
    // You can add more as needed (e.g., "InitializeSimulation", "PerformParticleSwarmOperations", etc.).
    const char *allowedFuncs[] = {
        "main",
       // "InitializeSimulation",                 
       // "SetupGridAndVectors",
       // "UpdateCartesianVelocity",
      // "InterpolateFieldFromCornerToCenter",
       // "WriteSimulationFields",
       // "ReadSimulationFields",
       //  "GatherAllBoundingBoxes",
        // "BroadcastAllBoundingBoxes",
       "PerformParticleSwarmOperations"
      // "CreateParticleSwarm",
      // "AssignInitialPropertiesToSwarm",
      // "InitializeParticleBasicProperties",
      // "FinalizeSwarmSetup",
      //  "LocateAllParticlesInGrid",
      //"InterpolateParticleVelocities",
           //"ComputeTrilinearWeights",
       //"FinalizeSimulation"
    };
    set_allowed_functions(allowedFuncs, 2);

    // Enable PETSc default logging
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);

    registerEvents();   

    print_log_level();

    // Check if user requested to read fields instead of updating
    ierr = PetscOptionsGetBool(NULL, NULL, "-read_fields", &readFields, NULL); CHKERRQ(ierr);

    // Initialize simulation: user context, MPI rank/size, np, etc.
    ierr = InitializeSimulation(&user, &rank, &size, &np, &rstart, &ti, &block_number); CHKERRQ(ierr);

    // Another demonstration of LOG_ALLOW
    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,
              "main: readFields = %s, rank = %d, size = %d\n",
              readFields ? "true" : "false", rank, size);

    // Setup the computational grid
    ierr = SetupGridAndVectors(user, block_number); CHKERRQ(ierr);

    // Either update and write fields, or read fields from file
    if (!readFields) {
        ierr = SetAnalyticalCartesianField(user,"Ucat"); CHKERRQ(ierr);
        ierr = WriteSimulationFields(user); CHKERRQ(ierr);
    } else {
        ierr = ReadSimulationFields(user, ti); CHKERRQ(ierr);
    }

    // Compute and print maximum velocity magnitude
    ierr = VecNorm(user->Ucat, NORM_INFINITY, &umax); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_INFO,"Maximum velocity magnitude: %f\n", umax);

    // Gather bounding boxes on rank 0
    ierr = GatherAllBoundingBoxes(user, &bboxlist); CHKERRQ(ierr);

    // Broadcast bboxlist to all ranks
    ierr = BroadcastAllBoundingBoxes(user, &bboxlist); CHKERRQ(ierr);

    // Perform particle swarm operations with bboxlist knowledge on all ranks
    ierr = PerformParticleSwarmOperations(user, np, bboxlist); CHKERRQ(ierr);
    
    // Create an ASCII viewer to write log output to file
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "petsc_log.txt", &logviewer);

    // Print PETSc logging results at the end
    ierr = PetscLogView(logviewer); CHKERRQ(ierr);

    // Finalize simulation
    ierr = FinalizeSimulation(user, block_number, bboxlist); CHKERRQ(ierr);

    return 0;
}


