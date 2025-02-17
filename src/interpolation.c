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
//#include "common.h"         // Shared type definitions
//#include "ParticleSwarm.h"  // Particle swarm functions
//#include "walkingsearch.h"  // Particle location functions
//#include "grid.h"           // Grid functions
//#include "logging.h"        // Logging macros
//#include "io.h"             // Data Input and Output functions
#include "interpolation.h"


#define NUM_WEIGHTS 8 // Number of weights in trilinear interpolation

/**
 * @brief Computes the trilinear interpolation weights from the interpolation coefficients.
 *
 * This function computes the weights for trilinear interpolation at the eight corners of a cell
 * using the interpolation coefficients provided along the x, y, and z directions.
 *
 * @param[in]  a1 Interpolation coefficient along the x-direction (normalized coordinate within the cell).
 * @param[in]  a2 Interpolation coefficient along the y-direction (normalized coordinate within the cell).
 * @param[in]  a3 Interpolation coefficient along the z-direction (normalized coordinate within the cell).
 * @param[out] w  Array of 8 weights, each corresponding to one corner of the cell.
 *
 * @note
 * - The coefficients `a1`, `a2`, and `a3` should be in the range [0, 1].
 * - The order of weights corresponds to the eight corners of a hexahedral cell.
 */
static inline void ComputeTrilinearWeights(PetscReal a1, PetscReal a2, PetscReal a3, PetscReal *w) {
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ComputeTrilinearWeights: Computing weights for a1=%f, a2=%f, a3=%f.\n", a1, a2, a3);

    // Ensure a1, a2, a3 are within [0,1]
    a1 = PetscMax(0.0, PetscMin(1.0, a1));
    a2 = PetscMax(0.0, PetscMin(1.0, a2));
    a3 = PetscMax(0.0, PetscMin(1.0, a3));

    // Compute complementary coefficients (1 - a)
    PetscReal oa1 = 1.0 - a1;  
    PetscReal oa2 = 1.0 - a2;  
    PetscReal oa3 = 1.0 - a3;

    // Compute weights for each corner of the cell
    w[0] = a1  * a2  * a3;   // Front-bottom-left
    w[1] = a1  * oa2 * oa3;  // Front-top-right
    w[2] = oa1 * a2  * oa3;  // Back-bottom-right
    w[3] = a1  * a2  * oa3;  // Front-bottom-right
    w[4] = oa1 * oa2 * a3;   // Back-top-left
    w[5] = oa1 * oa2 * oa3;  // Back-top-right
    w[6] = oa1 * a2  * a3;   // Back-bottom-left
    w[7] = a1  * oa2 * a3;   // Front-top-left

    // Log the computed weights for debugging
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ComputeTrilinearWeights: Weights computed - "
        "w0=%f, w1=%f, w2=%f, w3=%f, w4=%f, w5=%f, w6=%f, w7=%f. \n",
        w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7]);
}

/**
 * @brief Computes the trilinear interpolated velocity at a given point.
 *
 * @param[in]  ucat   3D array of velocity field from a DMDA (indexed as [k][j][i]), each cell of type Cmpnts.
 * @param[in]  i,j,k  Integral cell indices (the "lower" corner in each dimension).
 * @param[in]  a1,a2,a3  Normalized coordinates within the cell ([0,1] range).
 * @param[out] vel    Pointer to a Cmpnts struct that will hold the interpolated velocity (x, y, z).
 *
 * The function uses the 8-corner trilinear formula with `ComputeTrilinearWeights()`.
 * If a different scheme is desired, implement a new function with the same interface.
 */
static inline PetscErrorCode InterpolateTrilinearVelocity(
    Cmpnts ***ucat,
    PetscInt i, PetscInt j, PetscInt k,
    PetscReal a1, PetscReal a2, PetscReal a3,
    Cmpnts *vel)
{
    PetscFunctionBegin; // PETSc macro for error/stack tracing

    // Compute the 8 corner weights
    PetscReal wcorner[8];
    ComputeTrilinearWeights(a1, a2, a3, wcorner);

    // Offsets for cell corners
    PetscInt i1 = i + 1;
    PetscInt j1 = j + 1;
    PetscInt k1 = k + 1;

    // Initialize velocity
    vel->x = 0.0;
    vel->y = 0.0;
    vel->z = 0.0;

    // Corner 0 => (i, j, k)
    vel->x += wcorner[0] * ucat[k ][j ][i ].x;
    vel->y += wcorner[0] * ucat[k ][j ][i ].y;
    vel->z += wcorner[0] * ucat[k ][j ][i ].z;

    // Corner 1 => (i+1, j, k)
    vel->x += wcorner[1] * ucat[k ][j ][i1].x;
    vel->y += wcorner[1] * ucat[k ][j ][i1].y;
    vel->z += wcorner[1] * ucat[k ][j ][i1].z;

    // Corner 2 => (i, j+1, k)
    vel->x += wcorner[2] * ucat[k ][j1][i ].x;
    vel->y += wcorner[2] * ucat[k ][j1][i ].y;
    vel->z += wcorner[2] * ucat[k ][j1][i ].z;

    // Corner 3 => (i+1, j+1, k)
    vel->x += wcorner[3] * ucat[k ][j1][i1].x;
    vel->y += wcorner[3] * ucat[k ][j1][i1].y;
    vel->z += wcorner[3] * ucat[k ][j1][i1].z;

    // Corner 4 => (i, j, k+1)
    vel->x += wcorner[4] * ucat[k1][j ][i ].x;
    vel->y += wcorner[4] * ucat[k1][j ][i ].y;
    vel->z += wcorner[4] * ucat[k1][j ][i ].z;

    // Corner 5 => (i+1, j, k+1)
    vel->x += wcorner[5] * ucat[k1][j ][i1].x;
    vel->y += wcorner[5] * ucat[k1][j ][i1].y;
    vel->z += wcorner[5] * ucat[k1][j ][i1].z;

    // Corner 6 => (i, j+1, k+1)
    vel->x += wcorner[6] * ucat[k1][j1][i ].x;
    vel->y += wcorner[6] * ucat[k1][j1][i ].y;
    vel->z += wcorner[6] * ucat[k1][j1][i ].z;

    // Corner 7 => (i+1, j+1, k+1)
    vel->x += wcorner[7] * ucat[k1][j1][i1].x;
    vel->y += wcorner[7] * ucat[k1][j1][i1].y;
    vel->z += wcorner[7] * ucat[k1][j1][i1].z;


    LOG_ALLOW(LOCAL, LOG_DEBUG, 
    "InterpolateTrilinearVelocity: Interpolated velocity at (i=%d, j=%d, k=%d) with a1=%.6f, a2=%.6f, a3=%.6f -> (x=%.6f, y=%.6f, z=%.6f).\n",
    i, j, k, a1, a2, a3, vel->x, vel->y, vel->z);


    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Completed particle velocity interpolation across all ranks.\n");


    PetscFunctionReturn(0);
}

/**
 * @brief Performs trilinear interpolation of velocities from grid to particles using interpolation coefficients.
 *
 * This function interpolates velocities for particles based on the velocity field defined on the computational grid.
 * It retrieves cell indices for each particle from the "DMSwarm_CellID" field and assigns the interpolated velocity
 * using trilinear interpolation from the surrounding grid cells. The function handles boundary checks and ensures
 * that particles outside the valid grid range are appropriately managed.
 *
 * Key Steps:
 * 1. Retrieve the number of local particles and their associated data fields (cell indices and velocities).
 * 2. Map the global velocity field (`Ucat`) to the local portion for efficient access.
 * 3. Retrieve interpolation coefficients (`a1`, `a2`, `a3`) for each particle.
 * 4. Compute trilinear interpolation weights and interpolate velocities from the surrounding grid cells.
 * 5. Handle edge cases where particles may lie outside the valid grid range.
 * 6. Restore all data fields and ensure PETSc arrays and vectors are correctly finalized.
 *
 * @param[in] user Pointer to the `UserCtx` structure containing:
 *                 - `user->da`: DMDA for the grid.
 *                 - `user->swarm`: DMSwarm for particles.
 *                 - `user->Ucat`: Global velocity vector on the grid.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InterpolateParticleVelocities(UserCtx *user) {
    PetscErrorCode ierr;
    PetscInt n_local;             // Number of local particles
    PetscInt64 *cellIDs = NULL;     // Array to store cell indices for each particle
    PetscReal *velocities = NULL; // Array to store particle velocities
    PetscReal *weights = NULL;    // Array to store interpolation weights
    Cmpnts ***ucat;               // 3D array to map local grid velocities
    Cmpnts uinterp;               // Temporary variable to store interpolated velocities.
    PetscInt i,j,k;
    PetscReal a1,a2,a3;           // The weights of a particles are stored here for trilinear coefficient calculation.
    DM fda = user->fda;           // Field DA 
    DM swarm = user->swarm;       // DMSwarm for the particles


    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Starting particle velocity interpolation.\n");

    // Verify global velocity field
    PetscReal max_val;
    ierr = VecNorm(user->Ucat, NORM_INFINITY, &max_val); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Global velocity field Ucat maximum magnitude: %f \n", max_val);

    // Retrieve the number of local particles
    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Found %d local particles.\n", n_local);

    // Retrieve particle data fields
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr); // Ensure 'weight' field exists

    // Access the local portion of the global velocity vector (Ucat) using 'fda'
    ierr = DMDAVecGetArrayRead(user->fda,user->Ucat,&ucat);CHKERRQ(ierr);


    LOG_ALLOW(GLOBAL, LOG_DEBUG, "InterpolateParticleVelocities: Starting velocity assignment for particles.\n");

    // Log grid dimensions
    LOG_ALLOW(GLOBAL, LOG_INFO, "Grid dimensions: mx=%d, my=%d, mz=%d \n", user->info.mx, user->info.my, user->info.mz);

    // Loop over all local particles
    for (PetscInt p = 0; p < n_local; p++) {
        // Retrieve cell indices for the particle
        i = cellIDs[3 * p];
        j = cellIDs[3 * p + 1];
        k = cellIDs[3 * p + 2];

	LOG_LOOP_ALLOW(GLOBAL, LOG_DEBUG, p, 10, "Particle %d: Host Cell = (%d, %d, %d)\n", p, i, j, k);

	// Clamp i, j, k to [0..mx-2], [0..my-2], [0..mz-2]
	if (i >= user->info.mx - 1) i = user->info.mx - 2;
	if (j >= user->info.my - 1) j = user->info.my - 2;
	if (k >= user->info.mz - 1) k = user->info.mz - 2;

        // Validate cell indices (boundary check)
        if (i < 0 || j < 0 || k < 0 || i >= user->info.mx || j >= user->info.my || k >= user->info.mz) {
            LOG_ALLOW(GLOBAL, LOG_WARNING, "Particle %d has invalid cell indices (%d, %d, %d)\n. Skipping interpolation.\n", p, i, j, k);
            velocities[3 * p    ] = 0.0;
            velocities[3 * p + 1] = 0.0;
            velocities[3 * p + 2] = 0.0;
            continue;
        }

        // Retrieve a1, a2, a3 from the 'weights' field (if that's where you're storing them)
        a1 = weights[3*p + 0];
        a2 = weights[3*p + 1];
        a3 = weights[3*p + 2];

        
	// Apply interpolation method to obtain velocity at the particle location
	//  ierr = InterpolateTrilinearVelocity(ucat,i,j,k,a1,a2,a3,&uinterp);
	
        // Assign interpolated velocity to the particle

        // zeroth order Interpolation
        
        velocities[3 * p]     = ucat[k+1][j+1][i+1].x; // u-component
        velocities[3 * p + 1] = ucat[k+1][j+1][i+1].y; // v-component
        velocities[3 * p + 2] = ucat[k+1][j+1][i+1].z; // w-component
       
	LOG_LOOP_ALLOW(GLOBAL, LOG_DEBUG, p, 10, "Particle %d: Interpolated velocity: (velocity.x=%f, velocity.y=%f, velocity.z=%f).\n",
                        p, velocities[3 * p], velocities[3 * p + 1], velocities[3 * p + 2]);	
    }

    // Restore the local velocity array
    ierr = DMDAVecRestoreArrayRead(fda,user->Ucat, &ucat); CHKERRQ(ierr);

    // Restore particle data fields
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Particle velocity interpolation completed.\n");

    // Ensure all ranks finish interpolation before proceeding
    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: All ranks completed particle interpolation.\n");

    return 0;
}
