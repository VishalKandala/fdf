/**
 * @file ParticleSwarm.h
 * @brief Header file for Particle Swarm management functions.
 *
 * This file contains declarations of functions responsible for creating, managing,
 * initializing, migrating, and printing particle swarms within a simulation using PETSc's DMSwarm.
 */

#ifndef PARTICLE_SWARM_H
#define PARTICLE_SWARM_H

// Include necessary headers
#include <petsc.h>        // PETSc library header
#include <petscdmswarm.h> // PETSc DMSwarm header
#include <stdbool.h>
#include <math.h>
#include "common.h"       // Common type definitions
#include "logging.h"      // Logging macros and definitions
#include "walkingsearch.h"
#include "Metric.h"
#include "io.h"
// --------------------- Function Declarations ---------------------

/**
 * @brief Creates and initializes a Particle Swarm.
 *
 * This function sets up a DMSwarm within the provided UserCtx structure, initializes
 * particle fields, and distributes particles across MPI processes. It ensures that
 * the number of particles is evenly divided among the available MPI ranks. If the total
 * number of particles isn't divisible by the number of processes, the remainder is distributed
 * to the first few ranks.
 *
 * Additionally, it now takes a 'bboxlist' array as an input parameter and passes it on to
 * AssignInitialProperties(), enabling particle initialization at the midpoint of each rank's
 * bounding box if ParticleInitialization is set to 0.
 *
 * @param[in,out] user          Pointer to the UserCtx structure containing the simulation context.
 * @param[in]     numParticles  Total number of particles to create across all MPI processes.
 * @param[in]     bboxlist      Pointer to an array of BoundingBox structures, one per rank.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 *
 * @note
 * - Ensure that `numParticles` is a positive integer.
 * - The `control.dat` file should contain necessary PETSc options.
 * - The `bboxlist` array should be properly populated before calling this function.
 */
PetscErrorCode CreateParticleSwarm(UserCtx *user, PetscInt numParticles, PetscInt *particlesPerProcess, BoundingBox *bboxlist);

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
 * @param[out]    randx   Pointer to store the RNG for the x-coordinate.
 * @param[out]    randy   Pointer to store the RNG for the y-coordinate.
 * @param[out]    randz   Pointer to store the RNG for the z-coordinate.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeRandomGenerators(UserCtx *user, PetscRandom *randx, PetscRandom *randy, PetscRandom *randz);

/**
 * @brief Initializes random number generators for logical space operations [0.0, 1.0).
 *
 * This function creates and configures three separate PETSc random number generators,
 * one for each logical dimension (i, j, k or xi, eta, zeta equivalent).
 * Each RNG is configured to produce uniformly distributed real numbers in the interval [0.0, 1.0).
 * These are typically used for selecting owned cells or generating intra-cell logical coordinates.
 *
 * @param[out]   rand_logic_i Pointer to store the RNG for the i-logical dimension.
 * @param[out]   rand_logic_j Pointer to store the RNG for the j-logical dimension.
 * @param[out]   rand_logic_k Pointer to store the RNG for the k-logical dimension.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeLogicalSpaceRNGs(PetscRandom *rand_logic_i, PetscRandom *rand_logic_j, PetscRandom *rand_logic_k);

/**
 * @brief Initializes all particle properties in the swarm.
 *
 * This function orchestrates the initialization of particle properties.
 * It first determines the inlet face if surface initialization (Mode 0) is selected
 * by parsing "bcs.dat".
 * Then, it initializes basic particle properties (physical position, Particle ID,
 * and placeholder Cell IDs) by calling `InitializeParticleBasicProperties`. This call
 * uses the provided `rand_logic_i/j/k` RNGs, which must be pre-initialized for [0,1).
 * The `rand_phys_x/y/z` RNGs (physically bounded) are passed but may not be used by
 * `InitializeParticleBasicProperties` for position setting if all initialization paths
 * use logical-to-physical mapping.
 * Finally, it calls helper functions to initialize other registered swarm fields
 * like "velocity", "weight", and "P" (pressure) to default values.
 *
 * @param[in,out] user               Pointer to the `UserCtx` structure.
 * @param[in]     particlesPerProcess Number of particles assigned to this MPI process.
 * @param[in]     rand_phys_x        RNG for physical x-coordinates (from `InitializeRandomGenerators`).
 * @param[in]     rand_phys_y        RNG for physical y-coordinates (from `InitializeRandomGenerators`).
 * @param[in]     rand_phys_z        RNG for physical z-coordinates (from `InitializeRandomGenerators`).
 * @param[in]     rand_logic_i       RNG for i-logical dimension tasks [0,1) (from `InitializeLogicalSpaceRNGs`).
 * @param[in]     rand_logic_j       RNG for j-logical dimension tasks [0,1) (from `InitializeLogicalSpaceRNGs`).
 * @param[in]     rand_logic_k       RNG for k-logical dimension tasks [0,1) (from `InitializeLogicalSpaceRNGs`).
 * @param[in]     bboxlist           Array of BoundingBox structures (potentially unused by IPBP).
 *
 * @return PetscErrorCode            Returns 0 on success, non-zero on failure.
 */
PetscErrorCode AssignInitialPropertiesToSwarm(UserCtx* user,
                                              PetscInt particlesPerProcess,
                                              PetscRandom *rand_phys_x, // RNG from original InitializeRandomGenerators
                                              PetscRandom *rand_phys_y, // RNG from original InitializeRandomGenerators
                                              PetscRandom *rand_phys_z, // RNG from original InitializeRandomGenerators
                                              PetscRandom *rand_logic_i, // RNG from InitializeLogicalSpaceRNGs
                                              PetscRandom *rand_logic_j, // RNG from InitializeLogicalSpaceRNGs
                                              PetscRandom *rand_logic_k, // RNG from InitializeLogicalSpaceRNGs
                                              BoundingBox *bboxlist);

/**
 * @brief Distributes particles evenly across MPI processes, handling any remainders.
 *
 * This function calculates the number of particles each MPI process should handle,
 * distributing the remainder particles to the first few ranks if necessary.
 *
 * @param[in]     numParticles        Total number of particles to create across all MPI processes.
 * @param[in]     rank                MPI rank of the current process.
 * @param[in]     size                Total number of MPI processes.
 * @param[out]    particlesPerProcess Number of particles assigned to the current MPI process.
 * @param[out]    remainder           Remainder particles when dividing numParticles by size.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode DistributeParticles(PetscInt numParticles, PetscMPIInt rank, PetscMPIInt size, PetscInt* particlesPerProcess, PetscInt* remainder);

/**
 * @brief Finalizes the swarm setup by destroying random generators and logging completion.
 *
 * This function cleans up resources by destroying random number generators and LOG_ALLOWs the completion of swarm setup.
 *
 * @param[in]     randx             Random number generator for the x-coordinate.
 * @param[in]     randy             Random number generator for the y-coordinate.
 * @param[in]     randz             Random number generator for the z-coordinate.
 * @param[in]     rand_logic_i      Random number generator for the xi-coordinate.
 * @param[in]     rand_logic_j      Random number generator for the eta-coordinate.
 * @param[in]     rand_logic_k      Random number generator for the zeta-coordinate.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode FinalizeSwarmSetup(PetscRandom *randx, PetscRandom *randy, PetscRandom *randz, PetscRandom *rand_logic_i, PetscRandom *rand_logic_j, PetscRandom *rand_logic_k);

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
                                         const PetscReal *positions, const PetscInt *cellIndices,
                                         Particle *particle);

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
                                        PetscReal *weights, PetscInt *cellIndices);


/**
 * @brief Locates all particles within the grid and calculates their interpolation weights.
 * @ingroup ParticleLocation
 *
 * This function iterates through all particles currently local to this MPI rank.
 * For each particle, it first checks if the particle is within the rank's
 * pre-calculated bounding box (`user->bbox`). If it is, it calls the
 * `LocateParticleInGrid` function to perform the walking search.
 *
 * `LocateParticleInGrid` is responsible for finding the containing cell `(i,j,k)`
 * and calculating the corresponding interpolation weights `(w1,w2,w3)`. It updates
 * the `particle->cell` and `particle->weights` fields directly upon success.
 * If the search fails (particle not found within MAX_TRAVERSAL, goes out of bounds,
 * or gets stuck without resolution), `LocateParticleInGrid` sets the particle's
 * `cell` to `{-1,-1,-1}` and `weights` to `{0.0, 0.0, 0.0}`.
 *
 * After attempting location, this function updates the corresponding entries in the
 * DMSwarm's "DMSwarm_CellID" and "weight" fields using the potentially modified
 * data from the `particle` struct.
 *
 * @param[in] user Pointer to the UserCtx structure containing grid, swarm, and bounding box info.
 *
 * @return PetscErrorCode Returns `0` on success, non-zero on failure (e.g., errors accessing DMSwarm fields).
 *
 * @note Assumes `user->bbox` is correctly initialized for the local rank.
 * @note Assumes `InitializeParticle` correctly populates the temporary `particle` struct.
 * @note Assumes `UpdateSwarmFields` correctly writes data back to the DMSwarm.
 */
PetscErrorCode LocateAllParticlesInGrid(UserCtx *user);


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
 * - Be cautious when logging in performance-critical code sections, especially if the function is called frequently.
 */
PetscBool IsParticleInsideBoundingBox(const BoundingBox *bbox, const Particle *particle);

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
PetscErrorCode UpdateParticleWeights(PetscReal *d, Particle *particle);

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
PetscErrorCode InitializeParticleSwarm(UserCtx *user, PetscInt np, BoundingBox *bboxlist);

#endif // PARTICLE_SWARM_H
