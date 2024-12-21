// ParticleSwarm.c

#include "ParticleSwarm.h"
#include "walkingsearch.h"
#include "logging.h"
#include <petsc.h>
#include <petscdmswarm.h>
#include <stdbool.h>
#include <math.h>

/**
 * @brief Initializes the DMSwarm object within the UserCtx structure.
 *
 * This function creates the DMSwarm, sets its type and dimension, and configures basic swarm properties.
 *
 * @param[in,out] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeSwarm(UserCtx* user) {
    PetscErrorCode ierr;  // Error code for PETSc functions

    // Create the DMSwarm object for particle management
    ierr = DMCreate(PETSC_COMM_WORLD, &user->swarm); CHKERRQ(ierr);
    ierr = DMSetType(user->swarm, DMSWARM); CHKERRQ(ierr);
    ierr = DMSetDimension(user->swarm, 3); CHKERRQ(ierr);
    ierr = DMSwarmSetType(user->swarm, DMSWARM_BASIC); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_INFO, "InitializeSwarm - DMSwarm created and configured.\n");

    return 0;
}


/**
 * @brief Registers necessary particle fields within the DMSwarm.
 *
 * This function registers fields such as position, velocity, CellID, and weight for each particle.
 *
 * @param[in,out] swarm   The DMSwarm object managing the particle swarm.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode RegisterParticleFields(DM swarm) {
    PetscErrorCode ierr;  // Error code for PETSc functions

    // Register particle fields
    ierr = DMSwarmRegisterPetscDatatypeField(swarm, "position", 3, PETSC_REAL); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "RegisterParticleFields - Registered field 'position'.\n");

    ierr = DMSwarmRegisterPetscDatatypeField(swarm, "velocity", 3, PETSC_REAL); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "RegisterParticleFields - Registered field 'velocity'.\n");

    ierr = DMSwarmRegisterPetscDatatypeField(swarm, "DMSwarm_CellID", 3, PETSC_INT64); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "RegisterParticleFields - Registered field 'DMSwarm_CellID'.\n");

    ierr = DMSwarmRegisterPetscDatatypeField(swarm, "weight", 3, PETSC_REAL); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "RegisterParticleFields - Registered field 'weight'.\n");

    // Finalize field registration
    ierr = DMSwarmFinalizeFieldRegister(swarm); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_INFO, "RegisterParticleFields - Finalized field registration.\n");

    return 0;
}

/**
 * @brief Initializes random number generators for assigning particle properties.
 *
 * This function creates and configures separate PETSc random number generators for the x, y, and z coordinates.
 *
 * @param[in,out] user    Pointer to the UserCtx structure containing simulation context.
 * @param[out]    randx    Pointer to store the RNG for the x-coordinate.
 * @param[out]    randy    Pointer to store the RNG for the y-coordinate.
 * @param[out]    randz    Pointer to store the RNG for the z-coordinate.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeRandomGenerators(UserCtx* user, PetscRandom *randx, PetscRandom *randy, PetscRandom *randz) {
    PetscErrorCode ierr;  // Error code for PETSc functions

    // Initialize RNG for x-coordinate
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, randx); CHKERRQ(ierr);
    ierr = PetscRandomSetType((*randx), PETSCRAND48); CHKERRQ(ierr);
    ierr = PetscRandomSetInterval(*randx, user->bbox.min_coords.x, user->bbox.max_coords.x); CHKERRQ(ierr);
    ierr = PetscRandomSeed(*randx); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "InitializeRandomGenerators - Initialized RNG for X-axis.\n");

    // Initialize RNG for y-coordinate
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, randy); CHKERRQ(ierr);
    ierr = PetscRandomSetType((*randy), PETSCRAND48); CHKERRQ(ierr);
    ierr = PetscRandomSetInterval(*randy, user->bbox.min_coords.y, user->bbox.max_coords.y); CHKERRQ(ierr);
    ierr = PetscRandomSeed(*randy); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "InitializeRandomGenerators - Initialized RNG for Y-axis.\n");

    // Initialize RNG for z-coordinate
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, randz); CHKERRQ(ierr);
    ierr = PetscRandomSetType((*randz), PETSCRAND48); CHKERRQ(ierr);
    ierr = PetscRandomSetInterval(*randz, user->bbox.min_coords.z, user->bbox.max_coords.z); CHKERRQ(ierr);
    ierr = PetscRandomSeed(*randz); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "InitializeRandomGenerators - Initialized RNG for Z-axis.\n");

    return 0;
}

/**
 * @brief Assigns initial positions, velocities, IDs, CellIDs, and weights to particles.
 *
 * This function populates the particle fields with initial data, including random positions within the domain,
 * zero velocities, unique particle IDs, default CellIDs, and zero weights.
 *
 * @param[in,out] user               Pointer to the UserCtx structure containing simulation context.
 * @param[in]     particlesPerProcess Number of particles assigned to the local MPI process.
 *
 * @param[in]     randx             Random number generator for the x-coordinate.
 * @param[in]     randy             Random number generator for the y-coordinate.
 * @param[in]     randz             Random number generator for the z-coordinate.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode AssignInitialProperties(UserCtx* user, PetscInt particlesPerProcess, PetscRandom *randx, PetscRandom *randy, PetscRandom *randz) {
    PetscErrorCode ierr;                   // Error code for PETSc functions
    DM swarm = user->swarm;                // DMSwarm object managing particles
    PetscReal *positions, *velocities, *weights; // Pointers to particle data fields
    PetscInt64 *particleIDs, *cellIDs;     // Pointers to particle ID and CellID fields
    PetscInt rank;           // MPI rank of the current process
   

    PetscPrintf(PETSC_COMM_WORLD," Initialization : %d \n",user->ParticleInitialization);

    if(user->ParticleInitialization==1){
      // Initialize random number generators
      ierr = InitializeRandomGenerators(user, randx, randy, randz); CHKERRQ(ierr);
    }
  
    // Determine the current process rank
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    // Access the 'position' field from the DMSwarm
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "AssignInitialProperties - Accessed field 'position'.\n");

    // Access the 'velocity' field from the DMSwarm
    ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "AssignInitialProperties - Accessed field 'velocity'.\n");

    // Access the 'DMSwarm_pid' field from the DMSwarm
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&particleIDs); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "AssignInitialProperties - Accessed field 'DMSwarm_pid'.\n");

    // Access the 'DMSwarm_CellID' field from the DMSwarm
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "AssignInitialProperties - Accessed field 'DMSwarm_CellID'.\n");

    // Access the 'weight' field from the DMSwarm
    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "AssignInitialProperties - Accessed field 'weight'.\n");

    // Initialize particle properties
    for (PetscInt p = 0; p < particlesPerProcess; p++) {

      if(user->ParticleInitialization==1){
          // Generate random positions within the simulation domain
          ierr = PetscRandomGetValue(*randx, &positions[p * 3 + 0]); CHKERRQ(ierr);
          ierr = PetscRandomGetValue(*randy, &positions[p * 3 + 1]); CHKERRQ(ierr);
          ierr = PetscRandomGetValue(*randz, &positions[p * 3 + 2]); CHKERRQ(ierr);
          LOG_DEFAULT(LOG_DEBUG, "AssignInitialProperties - Particle %d position set to (%.4f, %.4f, %.4f).\n",p, positions[p * 3 + 0], positions[p * 3 + 1], positions[p * 3 + 2]);
      }else {
   
          positions[3 * p + 0] = 0.6;
          positions[3 * p + 1] = 0.6;
          positions[3 * p + 2] = 0.6; 
      }
      
        // Initialize velocities to zero
        velocities[3 * p + 0] = 0.0;
        velocities[3 * p + 1] = 0.0;
        velocities[3 * p + 2] = 0.0;
        LOG_DEFAULT(LOG_DEBUG, "AssignInitialProperties - Particle %d velocity initialized to (0.0, 0.0, 0.0).\n", p);

        // Assign unique IDs to particles
        particleIDs[p] = rank * particlesPerProcess + p;
        LOG_DEFAULT(LOG_DEBUG, "AssignInitialProperties - Particle %d assigned ID %" PetscInt64_FMT ".\n", p, particleIDs[p]);

        // Initialize CellIDs to -1 (front, bottom, left corner grid node indices)
        cellIDs[3 * p + 0] = -1;
        cellIDs[3 * p + 1] = -1;
        cellIDs[3 * p + 2] = -1;
        LOG_DEFAULT(LOG_DEBUG, "AssignInitialProperties - Particle %d CellIDs initialized to (0, 0, 0).\n", p);

        // Initialize interpolation weights to zero
        weights[3 * p + 0] = 0.0;
        weights[3 * p + 1] = 0.0;
        weights[3 * p + 2] = 0.0;
        LOG_DEFAULT(LOG_DEBUG, "AssignInitialProperties - Particle %d weights initialized to (0.0, 0.0, 0.0).\n", p);
    }

    // Restore particle data fields
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&particleIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "AssignInitialProperties - Restored all particle fields.\n");

    return 0;
}


/**
 * @brief Distributes particles evenly across MPI processes, handling any remainders.
 *
 * This function calculates the number of particles each MPI process should handle,
 * distributing the remainder particles to the first few ranks if necessary.
 *
 * @param[in]     numParticles       Total number of particles to create across all MPI processes.
 * @param[in]     rank               MPI rank of the current process.
 * @param[in]     size               Total number of MPI processes.
 * @param[out]    particlesPerProcess Number of particles assigned to the current MPI process.
 * @param[out]    remainder           Remainder particles when dividing numParticles by size.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode DistributeParticles(PetscInt numParticles, PetscMPIInt rank, PetscMPIInt size, PetscInt* particlesPerProcess, PetscInt* remainder) {

    // Calculate the base number of particles per process
    *particlesPerProcess = numParticles / size;
    *remainder = numParticles % size;

    // Distribute the remainder particles to the first 'remainder' ranks
    if (rank < *remainder) {
        *particlesPerProcess += 1;
        LOG_DEFAULT(LOG_INFO, "DistributeParticles - Rank %d receives an extra particle. Total: %d\n", rank, *particlesPerProcess);
    } else {
        LOG_DEFAULT(LOG_INFO, "DistributeParticles - Rank %d receives %d particles.\n", rank, *particlesPerProcess);
    }

    return 0;
}

/**
 * @brief Finalizes the swarm setup by destroying random generators and logging completion.
 *
 * This function cleans up resources by destroying random number generators and logs the completion of swarm setup.
 *
 * @param[in]     randx             Random number generator for the x-coordinate.
 * @param[in]     randy             Random number generator for the y-coordinate.
 * @param[in]     randz             Random number generator for the z-coordinate.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode FinalizeSwarmSetup(PetscRandom *randx, PetscRandom *randy, PetscRandom *randz) {
    PetscErrorCode ierr;  // Error code for PETSc functions
    PetscInt  ParticleInitialization; 

    ierr = PetscOptionsGetInt(NULL, NULL, "-pinit", &ParticleInitialization, NULL); CHKERRQ(ierr);
 
    if(ParticleInitialization==1){

      // Destroy random number generators to free resources
      ierr = PetscRandomDestroy(randx); CHKERRQ(ierr);
      ierr = PetscRandomDestroy(randy); CHKERRQ(ierr);
      ierr = PetscRandomDestroy(randz); CHKERRQ(ierr);
      LOG_DEFAULT(LOG_DEBUG, "FinalizeSwarmSetup - Destroyed all random number generators.\n");
    }else if(ParticleInitialization==0){
      LOG_DEFAULT(LOG_DEBUG, "FinalizeSwarmSetup - Not a Random Initialization of Particles.\n");
    }

    return 0;
}

/**
 * @brief Creates and initializes a Particle Swarm.
 *
 * This function sets up a DMSwarm within the provided UserCtx structure, initializes
 * particle fields, and distributes particles across MPI processes. It ensures that
 * the number of particles is evenly divided among the available MPI ranks. If the total
 * number of particles isn't divisible by the number of processes, the remainder is distributed
 * to the first few ranks.
 *
 * @param[in,out] user          Pointer to the UserCtx structure containing simulation context.
 * @param[in]     numParticles  Total number of particles to create across all MPI processes.
 *
 * @return PetscErrorCode         Returns 0 on success, non-zero on failure.
 *
 * @note
 * - Ensure that `numParticles` is a positive integer.
 * - The `control.dat` file should contain necessary PETSc options.
 */
PetscErrorCode CreateParticleSwarm(UserCtx *user, PetscInt numParticles) {
    PetscErrorCode ierr;                      // Error code for PETSc functions
    PetscMPIInt rank, size;                   // MPI rank and size
    PetscInt particlesPerProcess = 0;         // Number of particles per MPI process
    PetscInt remainder = 0;                   // Remainder particles
    PetscReal domainLengthY, domainLengthZ;   // Simulation domain dimensions in y and z
    PetscRandom randx,randy,randz;            // Random Number Generators if assigning particles random locations.
 
    // Validate input parameters
    if (numParticles <= 0) {
        LOG_DEFAULT(LOG_ERROR, "CreateParticleSwarm - Number of particles must be positive. Given: %d\n", numParticles);
        return PETSC_ERR_ARG_OUTOFRANGE;  
    }

    // Insert PETSc options from "control.dat" into the option database
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, "control.dat", PETSC_TRUE); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "CreateParticleSwarm - Inserted options from control.dat\n");

    // Retrieve simulation domain dimensions from PETSc options
    // Uncomment and set Lx if needed
    // ierr = PetscOptionsGetReal(NULL, NULL, "-L_x", &domainLengthX, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-L_y", &domainLengthY, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-L_z", &domainLengthZ, NULL); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "CreateParticleSwarm - Domain dimensions: Lx=%.2f, Ly=%.2f, Lz=%.2f\n", 
                1.0, domainLengthY, domainLengthZ); // Assuming Lx=1.0 if not retrieved

    // Initialize MPI rank and size
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_INFO, "CreateParticleSwarm - Rank %d out of %d processes.\n", rank, size);

    // Distribute particles across MPI processes
   // ierr = DistributeParticles(numParticles, rank, size, &particlesPerProcess, &remainder); CHKERRQ(ierr);

    // Initialize the DMSwarm
    ierr = InitializeSwarm(user); CHKERRQ(ierr);

    // Register particle fields
    ierr = RegisterParticleFields(user->swarm); CHKERRQ(ierr);

    // Set the local number of particles and number of fields per particle (4 particles buffer is added for migration).
    ierr = DMSwarmSetLocalSizes(user->swarm, particlesPerProcess, 4); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_INFO, "CreateParticleSwarm - Set local swarm size: %d particles.\n", particlesPerProcess);

    // Log detailed information about the swarm for debug.
    if(get_log_level()==LOG_DEBUG){
      
      LOG(GLOBAL, LOG_INFO, "CreateParticleSwarm - Viewing DMSwarm:\n");
      ierr = DMView(user->swarm, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    }

      // Assign initial properties to particles
    ierr = AssignInitialProperties(user, particlesPerProcess,&randx,&randy,&randz); CHKERRQ(ierr);

    // Finalize swarm setup by destroying RNGs
    ierr = FinalizeSwarmSetup(&randx, &randy, &randz); CHKERRQ(ierr);

    LOG_DEFAULT(LOG_INFO, "CreateParticleSwarm - Particle swarm creation and initialization complete.\n");

    return 0;
}

/**
 * @brief Prints the coordinates of all particles in the swarm.
 *
 * This function retrieves the local number of particles and their coordinates
 * from the DMSwarm associated with the provided UserCtx. It then prints out
 * the coordinates of each particle in a synchronized manner across all MPI processes.
 * The output includes the MPI rank, global particle ID, local particle index, and
 * the (x, y, z) coordinates.
 *
 * @param[in] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode PrintParticleCoordinates(UserCtx* user) {
    DM swarm = user->swarm;            // DMSwarm object containing particles
    PetscErrorCode ierr;               // Error code for PETSc functions
    PetscInt localNumParticles;        // Number of particles on the local MPI process
    PetscReal *coordinates;            // Array to store particle coordinates
    PetscMPIInt rank;                  // MPI rank of the current process

    // Retrieve the MPI rank of the current process
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_INFO, "PrintParticleCoordinates - Rank %d is retrieving particle coordinates.\n", rank);

    // Get the number of particles in the local swarm
    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PrintParticleCoordinates - Rank %d has %d particles.\n", rank, localNumParticles);

    // Access the 'position' field from the DMSwarm
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&coordinates); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PrintParticleCoordinates - Retrieved 'position' field.\n");

    // Iterate over each local particle and print its coordinates
    for (PetscInt i = 0; i < localNumParticles; i++) {
        // Calculate the global particle ID (assuming particles are evenly distributed)
        PetscInt64 globalParticleID = rank * localNumParticles + i + 1;

        // Synchronized printing to ensure orderly output across MPI processes
        ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,
            "Rank %d - Global Particle %" PetscInt64_FMT " - Local Particle %d : Coordinates = (%.6f, %.6f, %.6f)\n",
            rank, globalParticleID, i + 1,
            coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2]); CHKERRQ(ierr);
    }

    // Flush the synchronized output to ensure all messages are printed
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PrintParticleCoordinates - Completed printing coordinates on Rank %d.\n", rank);

    // Restore the 'position' field to clean up
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&coordinates); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PrintParticleCoordinates - Restored 'position' field.\n");

    return 0;
}

/**
 * @brief Prints the positions and associated metadata of all particles in the swarm.
 *
 * This function retrieves the local number of particles, their positions,
 * unique identifiers, and the MPI rank from the DMSwarm associated with the provided UserCtx.
 * It then prints out the positions of each particle along with their IDs and ranks
 * in a synchronized manner across all MPI processes. The output includes the MPI rank,
 * global particle ID, local particle index, position coordinates, and associated metadata.
 *
 * @param[in] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode PrintParticleFields(UserCtx* user) {
    DM swarm = user->swarm;                // DMSwarm object containing particles
    PetscErrorCode ierr;                   // Error code for PETSc functions
    PetscInt localNumParticles;            // Number of particles on the local MPI process
    PetscReal *positions;                  // Array to store particle positions.
    PetscInt64 *particleIDs;               // Array to store particle unique IDs.
    PetscMPIInt *particleRanks;            // Array to store particle MPI ranks.
    PetscMPIInt rank;                      // MPI rank of the current process.
    PetscInt64  *cellIDs;                  // Array to store (host)cell IDs of particles.
    PetscReal *weights;                    // Array to store particle weights.
    PetscReal *velocities;                    // Array to store particle velocities.

    // Retrieve the MPI rank of the current process
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_INFO, "PrintParticleFields - Rank %d is retrieving particle positions.\n", rank);

    // Get the number of particles in the local swarm
    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PrintParticleFields - Rank %d has %d particles.\n", rank, localNumParticles);

    // Access the 'position', 'DMSwarm_pid','DMSwarm_rank','DMSwarm_CellID' and 'weights' fields from the DMSwarm
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PrintParticleFields - Retrieved 'position' field.\n");

    ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&particleIDs); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PrintParticleFields - Retrieved 'DMSwarm_pid' field.\n");

    ierr = DMSwarmGetField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&particleRanks); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PrintParticleFields - Retrieved 'DMSwarm_rank' field.\n");

    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PrintParticleFields - Retrieved 'DMSwarm_CellID' field.\n");

    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PrintParticleFields - Retrieved 'weight' field.\n");

    ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PrintParticleFields - Retrieved 'velocity' field.\n");

    ierr = PetscPrintf(PETSC_COMM_WORLD,"___________________________________________________________________________________________________________________________________________________________________\n");
    ierr = PetscPrintf(PETSC_COMM_WORLD,"|Rank | PID | Host IDs: i,j,k |        Position : x,y,z                   |           Velocity: x,y,z                    |           Weights: a1,a2,a3            | \n");  
    // Iterate over each local particle and print its position and metadata
    for (PetscInt i = 0; i < localNumParticles; i++) {
        // Synchronized printing to ensure orderly output across MPI processes
         ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"___________________________________________________________________________________________________________________________________________________________________\n");
         ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,
				       "|  %d  |  %d  |   %d,  %d,  %d   | %.6f, %.6f, %.6f | %.6f, %.6f, %.6f |  %.6f, %.6f, %.6f  |\n",particleRanks[i],particleIDs[i],cellIDs[3 * i],cellIDs[3 * i + 1], cellIDs[3 * i + 2],positions[3 * i], positions[3 * i + 1], positions[3 * i + 2],velocities[3 * i],velocities[3 * i + 1], velocities[3 * i + 2],weights[3 * i],weights[3 * i + 1], weights[3 * i + 2]); CHKERRQ(ierr);
    }
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"___________________________________________________________________________________________________________________________________________________________________\n");
    // Add a blank line after each rank's output
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);     

    // Flush the synchronized output to ensure all messages are printed
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PrintParticleFields - Completed printing positions on Rank %d.\n", rank);

    // Restore the 'position', 'DMSwarm_pid','DMSwarm_rank','DMSwarm_CellID' fields to clean up
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&particleIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&particleRanks); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);

    LOG_DEFAULT(LOG_DEBUG, "PrintParticleFields - Restored all particle fields.\n");

    return 0;
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
    LOG_DEFAULT(LOG_INFO, "DefineBasicMigrationPattern - Rank %d out of %d processes.\n", rank, size);

    // Get the number of particles managed by the local MPI process
    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_INFO, "DefineBasicMigrationPattern - Rank %d handling %d particles.\n", rank, localNumParticles);

    // Allocate memory for the migration list
    ierr = PetscCalloc1(localNumParticles, &miglist); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "DefineBasicMigrationPattern - Allocated migration list for %d particles.\n", localNumParticles);

    // Initialize the migration list: assign each particle to migrate to the current rank by default
    for (PetscInt p = 0; p < localNumParticles; p++) {
        miglist[p] = rank;
    }
    LOG_DEFAULT(LOG_DEBUG, "DefineBasicMigrationPattern - Initialized migration list with default rank assignments.\n");

    // Define custom migration conditions based on the number of MPI processes
    if (size > 1) {
        // Example condition: Assign the first particle in rank 0 to migrate to rank 2
        if (rank == 0 && localNumParticles > 0) {
            miglist[0] = 2;
            LOG_DEFAULT(LOG_INFO, "DefineBasicMigrationPattern - Rank 0, Particle 0 assigned to migrate to Rank 2.\n");
        }

        // Additional custom conditions can be added here for other ranks
        // Example:
        // if(rank == 1 && localNumParticles > 1){
        //     miglist[1] = 3;
        //     LOG_DEFAULT(LOG_INFO, "DefineBasicMigrationPattern - Rank 1, Particle 1 assigned to migrate to Rank 3.\n");
        // }

        // ... add more custom conditions as needed ...
    }

    // Assign the migration list to the user context for later use
    user->miglist = miglist;
    LOG_DEFAULT(LOG_DEBUG, "DefineBasicMigrationPattern - Migration list assigned to user context.\n");

    /*
    // Optional: Debugging output to verify migration assignments
    for(PetscInt p = 0; p < localNumParticles; p++) {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,
            "DefineBasicMigrationPattern - Rank %d - miglist[%d] = %d\n",
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
    LOG_DEFAULT(LOG_INFO, "PerformBasicMigration - Rank %d is initiating migration.\n", rank);

    // Execute the migration pattern to define target ranks for particles
    ierr = DefineBasicMigrationPattern(user); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_INFO, "PerformBasicMigration - Migration pattern defined.\n");

    // Get the number of particles in the local swarm
    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PerformBasicMigration - Rank %d handling %d particles.\n", rank, localNumParticles);

    // Retrieve the migration list from the user context
    miglist = user->miglist;
    LOG_DEFAULT(LOG_DEBUG, "PerformBasicMigration - Retrieved migration list from user context.\n");

    /*
    // Optional: Debugging output to verify migration assignments before migration
    for(PetscInt p = 0; p < localNumParticles; p++) {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,
            "PerformBasicMigration - Rank %d - miglist[%d] = %d; user->miglist[p] = %d\n",
            rank, p, miglist[p], user->miglist[p]);
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "***********************\n");
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    */

    // Access the 'DMSwarm_rank' field from the DMSwarm to update particle ranks
    ierr = DMSwarmGetField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&rankval); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PerformBasicMigration - Retrieved 'DMSwarm_rank' field.\n");

    // Update the 'DMSwarm_rank' field based on the migration list
    for (PetscInt p = 0; p < localNumParticles; p++) {
        rankval[p] = miglist[p];
        LOG_DEFAULT(LOG_DEBUG, "PerformBasicMigration - Particle %d assigned to Rank %d.\n", p, rankval[p]);
    }

    /*
    // Optional: Debugging output to verify migration assignments after rank updates
    for(PetscInt p = 0; p < localNumParticles; p++) {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,
            "PerformBasicMigration - After change - Rank %d - rankval[%d] = %d; user->miglist[p] = %d\n",
            rank, p, rankval[p], user->miglist[p]);
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "***********************\n");
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    */

    // Restore the 'DMSwarm_rank' field after modification
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&rankval); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PerformBasicMigration - Restored 'DMSwarm_rank' field.\n");

    // Invoke the DMSwarm migration process to relocate particles based on updated ranks
    ierr = DMSwarmMigrate(swarm, removePoints); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_INFO, "PerformBasicMigration - DMSwarm migration executed.\n");

    // Free the allocated migration list to prevent memory leaks
    ierr = PetscFree(user->miglist); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_DEBUG, "PerformBasicMigration - Freed migration list memory.\n");

    // Synchronize all MPI processes to ensure migration completion before proceeding
    ierr = PetscBarrier(NULL); CHKERRQ(ierr);
    LOG_DEFAULT(LOG_INFO, "PerformBasicMigration - Migration synchronization completed.\n");

    return 0;
}

/**
 * @brief Initializes a Particle struct with data from DMSwarm fields.
 *
 * This helper function populates a Particle structure using data retrieved from DMSwarm fields.
 *
 * @param[in]     i            Index of the particle in the DMSwarm.
 * @param[in]     PIDs         Pointer to the array of particle IDs.
 * @param[in]     weights      Pointer to the array of particle weights.
 * @param[in]     positions    Pointer to the array of particle positions.
 * @param[in]     cellIndices  Pointer to the array of particle cell indices.
 * @param[out]    particle     Pointer to the Particle struct to initialize.
 *
 * @return PetscErrorCode  Returns `0` on success, non-zero on failure.
 */
PetscErrorCode InitializeParticle(PetscInt i, const PetscInt64 *PIDs, const PetscReal *weights,
                                         const PetscReal *positions, const PetscInt64 *cellIndices,
                                         Particle *particle) {
    PetscFunctionBeginUser;
    
    if (particle == NULL) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "InitializeParticle - Output Particle pointer is NULL. \n");
    }
    
    // Logging the start of particle initialization
    LOG(GLOBAL,LOG_INFO, "InitializeParticle - Initializing Particle [%D] with PID: %lld.\n", i, PIDs[i]);
    
    // Initialize PID
    particle->PID = PIDs[i];
    LOG(GLOBAL,LOG_DEBUG, "InitializeParticle - Particle [%D] PID set to: %lld.\n", i, particle->PID);
    
    // Initialize weights
    particle->weights.x = weights[3 * i];
    particle->weights.y = weights[3 * i + 1];
    particle->weights.z = weights[3 * i + 2];
    LOG(GLOBAL,LOG_DEBUG, "InitializeParticle - Particle [%D] weights set to: (%.6f, %.6f, %.6f).\n", 
        i, particle->weights.x, particle->weights.y, particle->weights.z);
    
    // Initialize locations
    particle->loc.x = positions[3 * i];
    particle->loc.y = positions[3 * i + 1];
    particle->loc.z = positions[3 * i + 2];
    LOG(GLOBAL,LOG_DEBUG, "InitializeParticle - Particle [%D] location set to: (%.6f, %.6f, %.6f).\n", 
        i, particle->loc.x, particle->loc.y, particle->loc.z);
    
    // Initialize velocities (assuming default zero; modify if necessary)
    particle->vel.x = 0.0;
    particle->vel.y = 0.0;
    particle->vel.z = 0.0;
    LOG(GLOBAL,LOG_DEBUG, "InitializeParticle - Particle [%D] velocities initialized to zero.\n", i);
    
    // Initialize cell indices
    particle->cell[0] = cellIndices[3 * i];
    particle->cell[1] = cellIndices[3 * i + 1];
    particle->cell[2] = cellIndices[3 * i + 2];
    LOG(GLOBAL,LOG_DEBUG, "InitializeParticle - Particle [%D] cell indices set to: [%lld, %lld, %lld].\n", 
        i, particle->cell[0], particle->cell[1], particle->cell[2]);
    
    // Logging the completion of particle initialization
    LOG(GLOBAL,LOG_INFO, "InitializeParticle - Completed initialization of Particle [%D].\n", i);
    
    PetscFunctionReturn(0);
}

/**
 * @brief Updates DMSwarm fields with data from a Particle struct.
 *
 * This helper function writes back the modified Particle data to the corresponding DMSwarm fields.
 *
 * @param[in] i            Index of the particle in the DMSwarm.
 * @param[in] particle     Pointer to the Particle struct containing updated data.
 * @param[in,out] weights  Pointer to the array of particle weights.
 * @param[in,out] cellIndices Pointer to the array of particle cell indices.
 *
 * @return PetscErrorCode  Returns `0` on success, non-zero on failure.
 */
PetscErrorCode UpdateSwarmFields(PetscInt i, const Particle *particle,
                                        PetscReal *weights, PetscInt64 *cellIndices) {
    PetscFunctionBeginUser;
    
    if (particle == NULL) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UpdateSwarmFields - Input Particle pointer is NULL.\n");
    }
    
    // Logging the start of swarm fields update
    LOG(GLOBAL,LOG_INFO, "Updating DMSwarm fields for Particle [%D].\n", i);
    
    // Update weights
    weights[3 * i]     = particle->weights.x;
    weights[3 * i + 1] = particle->weights.y;
    weights[3 * i + 2] = particle->weights.z;
    LOG(GLOBAL,LOG_DEBUG, "UpdateSwarmFields - Updated weights for Particle [%D]: (%.6f, %.6f, %.6f).\n", 
        i, weights[3 * i], weights[3 * i + 1], weights[3 * i + 2]);
    
    // Update cell indices
    cellIndices[3 * i]     = particle->cell[0];
    cellIndices[3 * i + 1] = particle->cell[1];
    cellIndices[3 * i + 2] = particle->cell[2];
    LOG(GLOBAL,LOG_DEBUG, "UpdateSwarmFields -  Updated cell indices for Particle [%D]: [%lld, %lld, %lld].\n", 
        i, cellIndices[3 * i], cellIndices[3 * i + 1], cellIndices[3 * i + 2]);
    
    // Logging the completion of swarm fields update
    LOG(GLOBAL,LOG_INFO, "UpdateSwarmFields  - Completed updating DMSwarm fields for Particle [%D].\n", i);
    
    PetscFunctionReturn(0);
}

/**
 * @brief Locates all particles within the grid and calculates their interpolation weights.
 *
 * This function iterates through all local particles, checks if they intersect the bounding box,
 * locates them within the grid using `LocateParticleInGrid`, and updates their interpolation weights
 * and cell indices accordingly.
 *
 * @param[in] user Pointer to the user-defined context containing grid and swarm information.
 *
 * @return PetscErrorCode Returns `0` on success, non-zero on failure.
 */
PetscErrorCode LocateAllParticlesInGrid(UserCtx *user) {
    PetscErrorCode ierr;
    PetscMPIInt rank, size;
    PetscInt localNumParticles;
    PetscReal *positions = NULL, *weights = NULL;
    PetscInt64 *cellIndices = NULL, *PIDs = NULL;
    PetscReal d[NUM_FACES];
    DM swarm = user->swarm;
    Particle particle;  // Reusable Particle struct
    
    PetscFunctionBeginUser;
    
    // Retrieve MPI rank and size
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    LOG(GLOBAL,LOG_INFO, "LocateAllParticlesInGrid - MPI Rank: %d, Size: %d.\n", rank, size);
    
    // Synchronize all processes to ensure all have reached this point
    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);
    LOG(GLOBAL,LOG_INFO, "LocateAllParticlesInGrid - All processes synchronized at barrier.\n");
    
    // Access DMSwarm fields
    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG(GLOBAL,LOG_INFO, "LocateAllParticlesInGrid - Number of local particles: %d.\n", localNumParticles);
    
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIndices); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&PIDs); CHKERRQ(ierr);
    LOG(GLOBAL,LOG_INFO, "LocateAllParticlesInGrid - DMSwarm fields accessed successfully.\n");
    
    // Iterate over each local particle
    for (PetscInt i = 0; i < localNumParticles; ++i) {
        // Initialize Particle struct with data from DMSwarm fields
        ierr = InitializeParticle(i, PIDs, weights, positions, cellIndices, &particle); CHKERRQ(ierr);
        
        // Log particle initialization
        LOG(GLOBAL,LOG_DEBUG, "LocateAllParticlesInGrid - Processing Particle [%D]: PID=%lld.\n", 
            i, particle.PID);
        
        // Check if the particle intersects the bounding box
        PetscBool particle_detected = IsParticleInsideBoundingBox(&(user->bbox), &particle);
        LOG(GLOBAL,LOG_DEBUG, "LocateAllParticlesInGrid - Particle [%D] intersected bounding box: %s.\n", 
            i, particle_detected ? "YES" : "NO");
        
        if (particle_detected) {
            // Locate the particle within the grid
            LOG(GLOBAL,LOG_DEBUG, "LocateAllParticlesInGrid - Locating Particle [%D] in grid. \n", i);
            ierr = LocateParticleInGrid(user, &particle, d); CHKERRQ(ierr);
            LOG(GLOBAL,LOG_DEBUG, "LocateAllParticlesInGrid - Particle [%D] located in cell [%lld, %lld, %lld].\n", 
                i, particle.cell[0], particle.cell[1], particle.cell[2]);

            // Update the weights of the particle for interpolation.
          ierr = UpdateParticleWeights(d,&particle);
	} // particle_detected
        
        // Update DMSwarm fields with possibly modified Particle data
        ierr = UpdateSwarmFields(i, &particle, weights, cellIndices); CHKERRQ(ierr);
    }
    
    // Restore DMSwarm fields to release resources
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIndices); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&PIDs); CHKERRQ(ierr);
    LOG(GLOBAL,LOG_INFO, "LocateAllParticlesInGrid - DMSwarm fields restored successfully.\n");
    
    // Logging the completion of particle location
    LOG(GLOBAL,LOG_INFO, "LocateAllParticlesInGrid - Completed locating all particles in grid.\n");
    
    PetscFunctionReturn(0);
}

/**
 * @brief Checks if a particle's location is within a specified bounding box.
 *
 * This function determines whether the given particle's location lies inside the provided bounding box.
 * It performs an axis-aligned bounding box (AABB) check by comparing the particle's coordinates to the
 * minimum and maximum coordinates of the bounding box in each dimension (x, y, z).
 *
 * Logging statements are included to provide detailed information about the function's execution.
 *
 * @param[in]  bbox     Pointer to the BoundingBox structure containing minimum and maximum coordinates.
 * @param[in]  particle Pointer to the Particle structure containing the particle's location and identifier.
 *
 * @return PetscBool    Returns `PETSC_TRUE` if the particle is inside the bounding box, `PETSC_FALSE` otherwise.
 *
 * @note
 * - The function assumes that the `bbox` and `particle` pointers are valid and non-NULL.
 * - The function includes logging statements that start with the function name.
 * - The `LOG_SCOPE` variable is used to distinguish between `GLOBAL` and `LOCAL` log outputs.
 * - Be cautious when logging in performance-critical code sections, especially if the function is called frequently.
 */
PetscBool IsParticleInsideBoundingBox(const BoundingBox *bbox, const Particle *particle)
{
    // Function name for logging purposes
    const char *funcName = "IsParticleInsideBoundingBox";

    // Validate input pointers
    if (!bbox) {
        // Log error message and return PETSC_FALSE
        LOG(LOCAL, LOG_ERROR, "%s: Error - 'bbox' pointer is NULL.", funcName);
        return PETSC_FALSE;
    }
    if (!particle) {
        LOG(LOCAL, LOG_ERROR, "%s: Error - 'particle' pointer is NULL.", funcName);
        return PETSC_FALSE;
    }

    // Extract particle location and bounding box coordinates
    const Cmpnts loc = particle->loc;
    const Cmpnts min_coords = bbox->min_coords;
    const Cmpnts max_coords = bbox->max_coords;

    // Log the particle location and bounding box coordinates for debugging
    LOG(LOCAL, LOG_DEBUG, "%s: Particle PID %lld location: (%.6f, %.6f, %.6f).\n", funcName, particle->PID, loc.x, loc.y, loc.z);
    LOG(LOCAL, LOG_DEBUG, "%s: BoundingBox min_coords: (%.6f, %.6f, %.6f), max_coords: (%.6f, %.6f, %.6f).\n",
        funcName, min_coords.x, min_coords.y, min_coords.z, max_coords.x, max_coords.y, max_coords.z);

    // Check if the particle's location is within the bounding box
    if ((loc.x >= min_coords.x && loc.x <= max_coords.x) &&
        (loc.y >= min_coords.y && loc.y <= max_coords.y) &&
        (loc.z >= min_coords.z && loc.z <= max_coords.z)) {
        // Particle is inside the bounding box
        LOG(LOCAL, LOG_DEBUG, "%s: Particle PID %lld is inside the bounding box.\n", funcName, particle->PID);
        return PETSC_TRUE;
    }

    // Particle is outside the bounding box
    LOG(LOCAL, LOG_DEBUG, "%s: Particle PID %lld is outside the bounding box.\n", funcName, particle->PID);
    return PETSC_FALSE;
}

/**
 * @brief Updates a particle's interpolation weights based on distances to cell faces.
 *
 * This function computes interpolation weights using distances to the six
 * cell faces (`d`) and updates the `weight` field of the provided particle.
 *
 * @param[in]  d        Pointer to an array of distances to the six cell faces.
 * @param[out] particle Pointer to the Particle structure whose weights are to be updated.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode UpdateParticleWeights(PetscReal *d, Particle *particle) {

    // Validate input pointers
    if (!d || !particle) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
                "UpdateParticleWeights - Null pointer argument (d or particle).");
    }

    // Validate distances
    for (PetscInt i = LEFT; i < NUM_FACES; i++) {
        if (d[i] <= 0) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                    "UpdateParticleWeights - Distances must be positive values.");
        }
    }

    // Log the input distances
    LOG(GLOBAL, LOG_DEBUG,
        "UpdateParticleWeights - Calculating weights with distances: "
        "[LEFT=%f, RIGHT=%f, BOTTOM=%f, TOP=%f, FRONT=%f, BACK=%f].\n",
        d[LEFT], d[RIGHT], d[BOTTOM], d[TOP], d[FRONT], d[BACK]);

    // Compute and update the particle's weights
    particle->weights.x = d[LEFT] / (d[LEFT] + d[RIGHT]);
    particle->weights.y = d[BOTTOM] / (d[BOTTOM] + d[TOP]);
    particle->weights.z = d[FRONT] / (d[FRONT] + d[BACK]);

    // Log the updated weights
    LOG(GLOBAL, LOG_DEBUG,
        "UpdateParticleWeights - Updated particle weights: x=%f, y=%f, z=%f.\n",
        particle->weights.x, particle->weights.y, particle->weights.z);

    return 0;
}

