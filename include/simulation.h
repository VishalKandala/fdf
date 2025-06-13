#ifndef SIMULATION_H
#define SIMULATION_H

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscsys.h>
#include <petscdmcomposite.h>
#include <petscsystypes.h>

// Include additional headers
#include "common.h"         // Shared type definitions
#include "ParticleSwarm.h"  // Particle swarm functions
#include "walkingsearch.h"  // Particle location functions
#include "grid.h"           // Grid functions
#include "logging.h"        // Logging macros
#include "io.h"             // Data Input and Output functions
#include "interpolation.h"  // Interpolation routines
#include "AnalyticalSolution.h" // Analytical Solution for testing
#include "ParticleMotion.h" // Functions related to motion of particles
#include "Boundaries.h"     //  Functions related to Boundary conditions
#include "setup.h"          // Functions  related to setup

/**
 * @brief Performs the complete initial setup for the particle simulation at time t=0.
 *
 * This includes:
 * 1. Initial locating of particles (based on their potentially arbitrary initial assignment).
 * 2. A preliminary migration cycle to ensure particles are on the MPI rank that owns
 *    their initial physical region.
 * 3. If `user->ParticleInitialization == 0` (Surface Init), re-initializes particles on the
 *    designated inlet surface. This ensures particles migrated to an inlet-owning rank
 *    are correctly distributed on that surface.
 * 4. A final locating of all particles to get their correct cell indices and interpolation weights.
 * 5. Interpolation of initial Eulerian fields to the particles.
 * 6. Scattering of particle data to Eulerian fields (if applicable).
 * 7. Outputting initial data if requested.
 *
 * @param user Pointer to the UserCtx structure.
 * @param currentTime The current simulation time (should be StartTime, typically 0.0).
 * @param step The current simulation step (should be StartStep, typically 0).
 * @param readFields Flag indicating if Eulerian fields were read from file (influences output).
 * @param bboxlist Array of BoundingBox structures for domain decomposition.
 * @param OutputFreq Frequency for writing output files.
 * @param StepsToRun Total number of simulation steps planned (used for output logic on setup-only runs).
 * @param StartStep The starting step of the simulation (used for output logic).
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformInitialSetup(UserCtx *user, PetscReal currentTime, PetscInt step,
                                   PetscBool readFields, const BoundingBox *bboxlist,
                                   PetscInt OutputFreq, PetscInt StepsToRun, PetscInt StartStep);

/**
 * @brief Initializes or updates the complete, consistent state of all Eulerian fields for a given timestep.
 *
 * This function is a high-level wrapper that orchestrates the entire process of preparing
 * the fluid fields for a single time step. It follows the standard procedure for a
 * curvilinear solver: first resolving contravariant velocities (`Ucont`) and then
 * converting them to Cartesian (`Ucat`).
 *
 * Its sequential operations are:
 * 1.  Update the INTERIOR of the domain:
 *     - For the initial step, it calls `SetInitialInteriorField` to generate values.
 *     - For subsequent steps, it calls the main fluid solver.
 *     - If restarting from a file, it reads the data, overwriting the whole field.
 *
 * 2.  Apply Boundary Conditions:
 *     - It then calls the modular `BoundarySystem_ExecuteStep` to enforce all configured
 *       boundary conditions on the domain edges.
 *
 * 3.  Convert to Cartesian and Finalize:
 *     - It calls `Contra2Cart` to compute `Ucat` from `Ucont`.
 *     - It calls `UpdateLocalGhosts` to ensure all parallel data is synchronized.
 *
 * @param user        Pointer to the UserCtx structure, containing all simulation data.
 * @param step        The current timestep number being processed.
 * @param StartStep   The initial timestep number of the simulation.
 * @param time        The current simulation time.
 * @param readFields  A boolean flag. If true, the simulation attempts to read fields
 *                    from files at the StartStep instead of generating them.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetEulerianFields(UserCtx *user, PetscInt step, PetscInt StartStep, PetscReal time, PetscBool readFields);

/**
 * @brief Executes the main time-marching loop for the particle simulation.
 *
 * This function performs the following steps repeatedly for `StepsToRun`:
 * 1. Updates/Sets the background fluid velocity field (Ucat) for the current step.
 * 2. If it's the very first step (StartStep=0), calls `PerformInitialSetup` for
 *    preliminary migration, location, interpolation, and output.
 * 3. Updates particle positions using velocity from the *previous* step's interpolation.
 * 4. Performs boundary checks and removes out-of-bounds particles.
 * 5. Migrates particles between MPI processors using `PerformSingleParticleMigrationCycle`.
 * 6. Advances simulation time and step counter.
 * 7. Locates particles in the grid based on their *new* positions.
 * 8. Interpolates the fluid velocity (from the *current* Ucat) to the new particle locations
 *    to get velocities for the *next* advection step.
 * 9. Scatters particle data back to Eulerian fields.
 *10. Logs errors and outputs data at specified `OutputFreq` intervals.
 *
 * @param user         Pointer to the UserCtx structure.
 * @param StartStep    The initial step number (e.g., 0 for a new run, >0 for restart).
 * @param StartTime    The simulation time corresponding to StartStep.
 * @param StepsToRun   The number of steps to execute in this run. If 0 and StartStep is 0,
 *                     only `PerformInitialSetup` is executed.
 * @param OutputFreq   Frequency (in number of steps) at which to output data and log errors.
 * @param readFields   Flag indicating whether to read initial fields (used by `SetEulerianFields`
 *                     and `PerformInitialSetup`).
 * @param bboxlist     Array of BoundingBox structures for domain decomposition, used for migration.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode AdvanceSimulation(UserCtx *user, PetscInt StartStep, PetscReal StartTime,
                                 PetscInt StepsToRun, PetscInt OutputFreq, PetscBool readFields,
                                 const BoundingBox *bboxlist);
/**
 * @brief Executes the main time-marching loop for the particle simulation. [TEST VERSION]
 *
 * This version uses the new, integrated `LocateAllParticlesInGrid_TEST` orchestrator
 * and the `ResetAllParticleStatuses` helper for a clean, robust, and understandable workflow.
 *
 * For each timestep, it performs:
 * 1. Sets the background fluid velocity field (Ucat) for the current step.
 * 2. Updates particle positions using velocity from the *previous* step's interpolation.
 * 3. Removes any particles that have left the global domain.
 * 4. **A single call to `LocateAllParticlesInGrid_TEST`**, which robustly handles all
 *    particle location and migration until the swarm is fully settled for the current timestep.
 * 5. Interpolates the current fluid velocity to the newly settled particle locations.
 * 6. Scatters particle data back to Eulerian fields.
 * 7. **Resets particle statuses to `NEEDS_LOCATION`** to prepare for the next timestep.
 * 8. Outputs data at specified intervals.
 *
 * @param user Pointer to the UserCtx structure.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode AdvanceSimulation_TEST(UserCtx *user, PetscInt StartStep, PetscReal StartTime,
                                      PetscInt StepsToRun, PetscInt OutputFreq, PetscBool readFields, BoundingBox *bboxlist);

/**
 * @brief Performs the complete initial setup for the particle simulation at time t=0. [TEST VERSION]
 *
 * This version uses the new, integrated `LocateAllParticlesInGrid_TEST` orchestrator,
 * which handles both location and migration in a single, robust, iterative process.
 *
 * Its sequential operations are:
 * 1. A single, comprehensive call to `LocateAllParticlesInGrid_TEST` to sort all particles
 *    to their correct owner ranks and find their initial host cells.
 * 2. If `user->ParticleInitialization == 0` (Surface Init), it re-initializes particles on the
 *    designated inlet surface, now that they are on the correct MPI ranks.
 * 3. A second call to `LocateAllParticlesInGrid_TEST` is needed after re-initialization to
 *    find the new, correct host cells for the surface-placed particles.
 * 4. Interpolates initial Eulerian fields to the settled particles.
 * 5. Scatters particle data to Eulerian fields (if applicable).
 * 6. Outputs initial data if requested.
 *
 * @param user Pointer to the UserCtx structure.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformInitialSetup_TEST(UserCtx *user, PetscReal currentTime, PetscInt step,
                                        PetscBool readFields, PetscInt OutputFreq,
                                        PetscInt StepsToRun, PetscInt StartStep,
					BoundingBox *bboxlist);
#endif  // SIMULATION_H
