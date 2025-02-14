/**
 * @file ParticleSwarm.h
 * @brief Header file for Particle Motion and migration related functions.
 *
 * This file contains declarations of functions responsible for moving and migrating particle swarms within a simulation using PETSc's DMSwarm.
 */

 #ifndef PARTICLE_MOTION_H
 #define PARTICLE_MOTION_H

// Include necessary headers
#include <petsc.h>        // PETSc library header
#include <petscdmswarm.h> // PETSc DMSwarm header
#include <stdbool.h>
#include <math.h>
#include "common.h"       // Common type definitions
#include "logging.h"      // Logging macros and definitions
#include "walkingsearch.h"  // Walking search function for particle migration   
/**
 * @brief Updates a particle's position based on its velocity and the timestep dt (stored in user->dt).
 *
 * @param[in]     user     Pointer to your UserCtx (must contain user->dt).
 * @param[in,out] position Pointer to the particle's current position (Cmpnts).
 * @param[in]     velocity Pointer to the particle's velocity (Cmpnts).
 *
 * @return PetscErrorCode  Returns 0 on success, or an error code on failure.
 */
 PetscErrorCode UpdateParticlePosition(UserCtx *user, Cmpnts *position, const Cmpnts *velocity);

/**
 * @brief Loops over all local particles in the DMSwarm, updating their positions
 *        based on velocity and the global timestep user->dt.
 * @param[in,out] user    Pointer to UserCtx (must contain dt).
 *
 * @return PetscErrorCode Returns 0 on success, or an error code on failure.
 */
 PetscErrorCode UpdateAllParticlePositions(UserCtx *user);


 #endif // PARTICLE_MOTION_H