#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscdmcomposite.h>

// Include additional headers
#include "common.h"         // Shared type definitions
#include "ParticleSwarm.h"  // Particle swarm functions
#include "walkingsearch.h"  // Particle location functions
#include "grid.h"           // Grid functions
#include "logging.h"        // Logging macros
#include "io.h"             // Data Input and Output functions

// Macros and constants
#define NUM_WEIGHTS 8 // Number of weights in trilinear interpolation

// Function declarations

/**
 * @brief Computes the trilinear interpolation weights from the interpolation coefficients.
 *
 * @param[in]  a1 Interpolation coefficient along the x-direction.
 * @param[in]  a2 Interpolation coefficient along the y-direction.
 * @param[in]  a3 Interpolation coefficient along the z-direction.
 * @param[out] w  Array of 8 weights to be computed.
 * void ComputeTrilinearWeights(PetscReal a1, PetscReal a2, PetscReal a3, PetscReal *w);
 */

/**
 * @brief Interpolates particle velocities using trilinear interpolation.
 *
 * @param[in] user Pointer to the user-defined context containing grid and swarm information.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InterpolateParticleVelocities(UserCtx *user);

/**
 * @brief Initializes the simulation context and reads runtime options.
 *
 * @param[out] user Pointer to the allocated UserCtx structure.
 * @param[out] rank MPI rank of the process.
 * @param[out] size Number of MPI processes.
 * @param[out] np Number of particles.
 * @param[out] rstart Flag to restart(1) or start from t = 0 (0).
 * @param[out] ti The timestep to start from if restarting.
 * @param[out] nblk Number of grid blocks.
 * @param[in] help message required by 'PetscInitialize'
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeSimulation(UserCtx **user, PetscInt *rank, PetscInt *size, PetscInt *np, PetscInt *rstart, PetscInt *ti, PetscInt *nblk);

/**
 * @brief Sets up the simulation grid and initializes vectors.
 *
 * @param[in,out] user Pointer to the UserCtx structure.
 * @param[in] block_number Number of grid blocks.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode SetupGridAndVectors(UserCtx *user, PetscInt block_number);

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
PetscErrorCode PerformParticleSwarmOperations(UserCtx *user, PetscInt np, BoundingBox *bboxlist);


/**
 * @brief Cleans up resources and finalizes the simulation.
 *
 * @param[in,out] user Pointer to the UserCtx structure.
 * @param[in] block_number Number of grid blocks.
 * @param[in]     bboxlist Pointer to an array of BoundingBox structures, one per MPI rank.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode FinalizeSimulation(UserCtx *user, PetscInt block_number, BoundingBox* bboxlist);

#endif // INTERPOLATION_H
