/**
 * @file inttest.c  //  Particle In Cell main.
 * @brief Test program for DMSwarm interpolation using the fdf-curvIB method.
 *
 * Initializes a particle swarm, reads velocity fields, and performs particle-grid interpolation.
 */

#include "simulation.h"

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
    PetscInt ActualStepsToRun;
    PetscInt np = 0;
    PetscInt OutputFreq;
    PetscReal StartTime = 0.0;
    PetscMPIInt rank, size;
    BoundingBox *bboxlist = NULL;    // Array of bounding boxes
    PetscReal umax;
    PetscBool readFields = PETSC_FALSE;
    static char help[] = " Test for interpolation - swarm-curvIB";
    PetscViewer logviewer;
    char allowedFile[PETSC_MAX_PATH_LEN]  = "config.dat";
    PetscBool useCfg = PETSC_FALSE;
    char **allowedFuncs  = NULL;
    PetscInt nAllowed   = 0;
    PetscBool OnlySetup = PETSC_FALSE;

    // -------------------- 1. PETSc Initialization --------------------
    ierr = PetscInitialize(&argc, &argv, (char *)0, help); CHKERRQ(ierr);

    // Initialize simulation: user context, MPI rank/size, np, etc.
    ierr = InitializeSimulation(&user, &rank, &size, &np, &StartStep, &StepsToRun, &StartTime, &block_number, &OutputFreq,&readFields,&allowedFuncs,&nAllowed,allowedFile,&useCfg,&OnlySetup); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL,LOG_INFO," Simulation Initialized \n");

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,
              "readFields = %s, size = %d, rank = %d\n",
               readFields ? "true" : "false", size,rank);

    // Setup the computational grid
    ierr = SetupGridAndVectors(user, block_number); CHKERRQ(ierr);

    LOG_ALLOW_SYNC(GLOBAL,LOG_INFO," Grid & Fields Setup on rank %d! \n",rank);

    LOG_ALLOW(GLOBAL,LOG_INFO," Simulation Fields %s \n",readFields ? "read":"generated");    

    // Setup the Domain Rank Information.
    ierr = SetupDomainRankInfo(user, &bboxlist);
    
    LOG_ALLOW_SYNC(GLOBAL,LOG_INFO,"Domain Decomposition Information setup on rank %d! \n",rank);

    ierr = SetupBoundaryConditions(user);

    LOG_ALLOW_SYNC(GLOBAL,LOG_INFO,"Boundary Condition system setup on rank %d! \n",rank);
    

    //if(get_log_level() == LOG_INFO && is_function_allowed(__func__)){
    //   print_log_level();
       // Compute and print maximum velocity magnitude
       //   ierr = VecNorm(user->Ucat, NORM_INFINITY, &umax); CHKERRQ(ierr);
    //   LOG_ALLOW(GLOBAL,LOG_INFO,"Maximum velocity magnitude:%f \n", umax);
    //  }
    
    
    // Initialize particle swarm with bboxlist knowledge on all ranks
    ierr = InitializeParticleSwarm(user, np, bboxlist); CHKERRQ(ierr);

    // Display Banner for simulation
    ierr = DisplayBanner(user, StartTime, StartStep, StepsToRun, size, np, bboxlist); CHKERRQ(ierr);


    // Setup Only Condition
     ActualStepsToRun=StepsToRun;
    if (OnlySetup) {
      LOG_ALLOW(GLOBAL, LOG_INFO, "SETUP ONLY MODE: Forcing StepsToRun to 0 for AdvanceSimulation call.\n");
      ActualStepsToRun = 0; // This will trigger the setup-only path in AdvanceSimulation
     }
    // Advance the Lagrangian Particle Simulation
    //  ierr = AdvanceSimulation(user,StartStep,StartTime,ActualStepsToRun,OutputFreq,readFields,bboxlist);
    //ierr = AdvanceSimulation_TEST(user,StartStep,StartTime,ActualStepsToRun,OutputFreq,bboxlist);
    
    // Finalize simulation
    ierr = FinalizeSimulation(user, block_number, bboxlist,allowedFuncs,nAllowed,&logviewer); CHKERRQ(ierr);
 
    return 0;
}
