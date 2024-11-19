/**
 * @file swarm_interp.c
 * @brief Main program for DMSwarm interpolation using the fdf-curvIB method.
 *
 * This program initializes a particle swarm, reads velocity fields, and performs particle-grid interpolation.
 */

static char help[] = "DMSwarm Interpolation - fdf-curvIB ";

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscdmcomposite.h>

// Include the updated headers
#include "common.h"         // Shared type definitions
#include "ParticleSwarm.h"  // Particle swarm functions
#include "walkingsearch.h"  // Particle location functions
#include "grid.h"           // Grid functions
#include "logging.h"        // Logging macros

// Global variables
PetscInt np = 0;        // Number of particles
PetscInt ti = 0;        // Time index
PetscReal L_dim = 1.0;  // Domain length (unused in this snippet)
PetscInt block_number = 1; // Number of blocks (assuming 1 for simplicity)
PetscInt visflg = 0;    // Visualization flag

// --------------------- Function Implementations ---------------------

/**
 * @brief Reads the velocity field (Ucat) from a binary file.
 *
 * This function reads the velocity field from a binary file and loads it into the UserCtx's Ucat vector.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode Ucat_Binary_Input(UserCtx *user)
{
    PetscViewer viewer;
    char filen[90];
    PetscInt ctr = 7;
    PetscErrorCode ierr;

    sprintf(filen, "results/ufield%5.5d_%1.1d.dat", ti, user->_this);

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer); CHKERRQ(ierr);

    PetscInt N;
    ierr = VecGetSize(user->Ucat, &N); CHKERRQ(ierr);
    if (visflg == ctr)
        PetscPrintf(PETSC_COMM_WORLD, "Ucat_Binary_Input - user, SizeOf(Ucat) - %p, %d \n", (void*)user, N);
    ierr = VecLoad(user->Ucat, viewer); CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    ierr = PetscBarrier(NULL); CHKERRQ(ierr);

    return 0;
}

/**
 * @brief Checks if a particle's location intersects with the given bounding box.
 *
 * @param[in] bbox     Pointer to the BoundingBox structure.
 * @param[in] particle Pointer to the Particle structure.
 *
 * @return PetscBool Returns PETSC_TRUE if the particle intersects with the bounding box, PETSC_FALSE otherwise.
 */
PetscBool CPUPointIntersectCheck(BoundingBox *bbox, Particle *particle)
{
    Cmpnts loc = particle->loc;
    Cmpnts min_coords = bbox->min_coords;
    Cmpnts max_coords = bbox->max_coords;
    PetscBool Intersects = PETSC_FALSE;

    if ((loc.x >= min_coords.x && loc.x <= max_coords.x) &&
        (loc.y >= min_coords.y && loc.y <= max_coords.y) &&
        (loc.z >= min_coords.z && loc.z <= max_coords.z)) {
        Intersects = PETSC_TRUE;
    }

    return Intersects;
}

/**
 * @brief Calculates interpolation weights based on distances to cell faces.
 *
 * @param[out] a Pointer to a Cmpnts structure to store the weights.
 * @param[in]  d Array of distances to the six cell faces.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InterpolationWeightsCalculate(Cmpnts *a, PetscReal d[6])
{
    a->x = d[0] / (d[0] + d[1]);
    a->y = d[2] / (d[2] + d[3]);
    a->z = d[4] / (d[4] + d[5]);
    return 0;
}

// --------------------- Main Function ---------------------

#undef __FUNCT__
#define __FUNCT__ "main"

/**
 * @brief Main function for DMSwarm interpolation.
 *
 * Initializes the grid, reads the velocity field, creates particles, and performs particle-grid interpolation.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 *
 * @return int Returns 0 on success, non-zero on failure.
 */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    UserCtx *user;
    PetscInt ctr = 0; // Output counter
    PetscInt rank, size, bi;
    PetscReal umax;
    BoundingBox *bboxlist;

    ierr = PetscInitialize(&argc, &argv, (char *)0, help); if (ierr) return ierr;

    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, "control.dat", PETSC_TRUE); CHKERRQ(ierr);

    ierr = PetscOptionsGetInt(NULL, NULL, "-visflg", &visflg, NULL); CHKERRQ(ierr);

    // Allocate memory for user context
    ierr = PetscMalloc1(block_number, &user); CHKERRQ(ierr);

    // Identify all processors
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    // Read number of particles and time index from control file
    ierr = PetscOptionsGetInt(NULL, NULL, "-ti", &ti, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-numParticles", &np, NULL); CHKERRQ(ierr);

    if (visflg == ctr) {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD, "main - user, rank - %p, %d  \n", (void*)user, rank);
        PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

        PetscPrintf(PETSC_COMM_WORLD, "main - ti, number of particles %d \n", ti, np);
    }

    // Define grid coordinates
    ierr = DefineGridCoordinates(user); CHKERRQ(ierr);

    // Create Vector for Ucat
    for (bi = 0; bi < block_number; bi++) {
        ierr = DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat); CHKERRQ(ierr);
    }

    // Read Ucat from ufield
    ierr = Ucat_Binary_Input(user); CHKERRQ(ierr);

    ierr = VecNorm(user->Ucat, NORM_INFINITY, &umax); CHKERRQ(ierr);

    if (visflg == ctr)
        PetscPrintf(PETSC_COMM_WORLD, "main - max(Ucat) - %f \n", umax);

    // Create Bounding Boxes for each rank
    ierr = GatherAllBoundingBoxes(user, &bboxlist); CHKERRQ(ierr);

    // Create and Initialize Particle Swarm
    ierr = CreateParticleSwarm(user, np); CHKERRQ(ierr);

    if (visflg == ctr) {
        PetscPrintf(PETSC_COMM_WORLD, "main - ParticleSwarm: Pre-Search. \n");
        PetscPrintf(PETSC_COMM_WORLD, "\n");
    }
    ierr = PrintParticlePositions(user); CHKERRQ(ierr);

    // Locate all particles in the grid
    ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);

    if (visflg == ctr) {
        PetscPrintf(PETSC_COMM_WORLD, "main - ParticleSwarm: Post-Search. \n");
        PetscPrintf(PETSC_COMM_WORLD, "\n");
    }
    ierr = PrintParticlePositions(user); CHKERRQ(ierr);

    // Finalize and clean up
    for (bi = 0; bi < block_number; bi++) {
        ierr = DMDestroy(&(user[bi].fda)); CHKERRQ(ierr);
        ierr = DMDestroy(&(user[bi].da)); CHKERRQ(ierr);
        ierr = DMDestroy(&(user[bi].swarm)); CHKERRQ(ierr);
    }

    ierr = PetscFree(user); CHKERRQ(ierr);

    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}
