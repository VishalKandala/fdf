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

#include "setup.h"

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
    PetscInt StartStep = 0;
    PetscInt StepsToRun;
    PetscInt np = 0;
    PetscInt OutputFreq;
    PetscReal StartTime = 0.0;
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
       //"DefineGridCoordinates",
       //"ParseGridInputs",
       // "DetermineGridSizes",
       //"InitializeGridDM",
       //"AssignGridCoordinates",
       //"FinalizeGridSetup",
       //"SetAnalyticalCartesianField",
       //"ApplyAnalyticalBC",
       // "ApplyAnalyticalBC_Vector",
       //"InterpolateFieldFromCornerToCenter_Vector",
       //"WriteSimulationFields",
       //"ReadSimulationFields",
       //"GatherAllBoundingBoxes",
       //"BroadcastAllBoundingBoxes",
       //"InitializeParticleSwarm",
       //"CreateParticleSwarm",
       //"AssignInitialPropertiesToSwarm",
       //"InitializeParticleBasicProperties",
       //"FinalizeSwarmSetup",
       "LocateAllParticlesInGrid",
       "InterpolateAllFieldsToSwarm",
       //"InterpolateEulerFieldToSwarm",
       //"InterpolateFieldFromCenterToCorner_Vector",
       //"ComputeTrilinearWeights",
       //"InterpolateEulerFieldToSwarmForParticle",
       //"TrilinearInterpolation_Vector",
       //"FinalizeSimulation"
       "AdvanceSimulation",
       //"SetEulerianFields"
       //"CheckAndRemoveOutOfBoundsParticles",
       //"IdentifyMigratingParticles",
       //"SetMigrationRanks",
       //"PerformMigration",
       //"SetDMDAProcLayout",
       //"ComputeAndStoreNeighborsRanks",
       "CalculateParticleCountPerCell",
	"ScatterAllParticleFieldsToEulerFields",
       "ScatterParticleFieldToEulerField",
       "ScatterParticleFieldToEulerField_Internal",
	"NormalizeGridVectorByCount",
	"AccumulateParticleField",
	"GetScatterTargetInfo"
    };
    set_allowed_functions(allowedFuncs, sizeof(allowedFuncs) / sizeof(allowedFuncs[0])); // Use sizeof for count

    // Enable PETSc default logging
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);

    registerEvents();   

    print_log_level();

    // Check if user requested to read fields instead of updating
    ierr = PetscOptionsGetBool(NULL, NULL, "-read_fields", &readFields, NULL); CHKERRQ(ierr);

    // Initialize simulation: user context, MPI rank/size, np, etc.
    ierr = InitializeSimulation(&user, &rank, &size, &np, &StartStep, &StepsToRun, &StartTime, &block_number, &OutputFreq); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL,LOG_INFO," Simulation Initialized \n");

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,
              "readFields = %s, size = %d, rank = %d\n",
		   readFields ? "true" : "false", size,rank);


    
    // Setup the computational grid
    ierr = SetupGridAndVectors(user, block_number); CHKERRQ(ierr);

    // Compute and Store the neighboring ranks for each rank.
    ierr = ComputeAndStoreNeighborRanks(user); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL,LOG_INFO," Grid & Fields Setup! \n");

    PetscPrintf(PETSC_COMM_WORLD," ---------------------------------------- \n");
    PetscPrintf(PETSC_COMM_WORLD," Grid: %d,%d,%d \n",user->IM,user->JM,user->KM);
    PetscPrintf(PETSC_COMM_WORLD," StartTime: %0.4f | Timestep: %0.4f \n",StartTime,user->dt);
    PetscPrintf(PETSC_COMM_WORLD," StartStep: %d | Simulation Length (Timesteps): %d \n",StartStep,StepsToRun);
    PetscPrintf(PETSC_COMM_WORLD," No.of Processors: %d \n",size);
    PetscPrintf(PETSC_COMM_WORLD," No.of Particles: %d \n",np);
    PetscPrintf(PETSC_COMM_WORLD," Particle Initialization Mode: %d \n",user->ParticleInitialization);
    PetscPrintf(PETSC_COMM_WORLD," Field Initialization Mode: %d \n",user->FieldInitialization);
    if(user->FieldInitialization == 0){
    PetscPrintf(PETSC_COMM_WORLD," Constant Velocity: %.4f \n",user->ConstantVelocity); 
    }
    PetscPrintf(PETSC_COMM_WORLD," ---------------------------------------- \n");

    LOG_ALLOW(GLOBAL,LOG_INFO," Simulation Fields %s \n",readFields ? "read":"generated");    

    if(get_log_level() == LOG_INFO && is_function_allowed(__func__)){
    // Compute and print maximum velocity magnitude
    ierr = VecNorm(user->Ucat, NORM_INFINITY, &umax); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_INFO,"Maximum velocity magnitude:%f \n", umax);
    }
    
    // Gather bounding boxes on rank 0
    ierr = GatherAllBoundingBoxes(user, &bboxlist); CHKERRQ(ierr);

    // Broadcast bboxlist to all ranks
    ierr = BroadcastAllBoundingBoxes(user, &bboxlist); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL,LOG_INFO," Bounding Boxes setup \n");

    // Initialize particle swarm with bboxlist knowledge on all ranks
    ierr = InitializeParticleSwarm(user, np, bboxlist); CHKERRQ(ierr);

    // Advance the Lagrangian Particle Simulation
    ierr = AdvanceSimulation(user,StartStep,StartTime,StepsToRun,OutputFreq,readFields,bboxlist);
 
    // Create an ASCII viewer to write log output to file
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "simulationlog.txt", &logviewer);

    LOG_ALLOW(GLOBAL,LOG_INFO," PETSC Logs written \n");
    
    // Print PETSc logging results at the end
    ierr = PetscLogView(logviewer); CHKERRQ(ierr);
    
    // Finalize simulation
    ierr = FinalizeSimulation(user, block_number, bboxlist); CHKERRQ(ierr);
 
    return 0;
}


