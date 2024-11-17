// ParticleSwarm.h

#ifndef PARTICLE_SWARM_H
#define PARTICLE_SWARM_H


#include <petsc.h>        // PETSc library header
#include <stdbool.h>
#include <math.h>
#include <petscdmswarm.h>
#include "petsctime.h"
#include "variables.h"


/**
 * @file ParticleSwarm.h
 * @brief Header file for Particle Swarm management functions.
 *
 * This file contains declarations of functions responsible for creating, managing,
 * initializing, migrating, and printing particle swarms within a simulation using PETSc's DMSwarm.
 */

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
PetscErrorCode CreateParticleSwarm(UserCtx *user, PetscInt numParticles);

/**
 * @brief Initializes the DMSwarm object within the UserCtx structure.
 *
 * This function creates the DMSwarm, sets its type and dimension, and configures basic swarm properties.
 *
 * @param[in,out] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeSwarm(UserCtx* user);

/**
 * @brief Registers necessary particle fields within the DMSwarm.
 *
 * This function registers fields such as position, velocity, CellID, and weight for each particle.
 *
 * @param[in,out] swarm   The DMSwarm object managing the particle swarm.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode RegisterParticleFields(DM swarm);

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
PetscErrorCode InitializeRandomGenerators(UserCtx *user, PetscRandom *randx, PetscRandom *randy, PetscRandom *randz);

/**
 * @brief Assigns initial positions, velocities, IDs, CellIDs, and weights to particles.
 *
 * This function populates the particle fields with initial data, including random positions within the domain,
 * zero velocities, unique particle IDs, default CellIDs, and zero weights.
 *
 * @param[in,out] user               Pointer to the UserCtx structure containing simulation context.
 * @param[in]     particlesPerProcess Number of particles assigned to the local MPI process.
 * @param[in]     randx             Random number generator for the x-coordinate.
 * @param[in]     randy             Random number generator for the y-coordinate.
 * @param[in]     randz             Random number generator for the z-coordinate.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode AssignInitialProperties(UserCtx *user, PetscInt particlesPerProcess, PetscRandom randx, PetscRandom randy, PetscRandom randz);

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
PetscErrorCode DistributeParticles(PetscInt numParticles, PetscMPIInt rank, PetscMPIInt size, PetscInt* particlesPerProcess, PetscInt* remainder);

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
PetscErrorCode FinalizeSwarmSetup(PetscRandom randx, PetscRandom randy, PetscRandom randz);

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
PetscErrorCode PrintParticleCoordinates(UserCtx* user);

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
PetscErrorCode PrintParticlePositions(UserCtx* user);

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
PetscErrorCode DefineBasicMigrationPattern(UserCtx* user);

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
PetscErrorCode PerformBasicMigration(UserCtx* user);

#endif // PARTICLE_SWARM_H
