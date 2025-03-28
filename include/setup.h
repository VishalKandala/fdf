#ifndef SETUP_H
#define SETUP_H

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


/* Macro to automatically select the correct allocation function */
#define Allocate3DArray(array, nz, ny, nx)  \
  _Generic((array),                          \
    PetscReal ****: Allocate3DArrayScalar,       \
    Cmpnts    ****: Allocate3DArrayVector            \
  )(array, nz, ny, nx)

/* Macro for deallocation */
#define Deallocate3DArray(array, nz, ny)  \
  _Generic((array),                        \
    PetscReal ***: Deallocate3DArrayScalar,     \
    Cmpnts    ***: Deallocate3DArrayVector          \
  )(array, nz, ny)


 /**
    * @brief Register events for logging.
    * 
    * Registers events for logging purposes using PetscLogEventRegister.
    *  
 */
PetscErrorCode registerEvents(void);

/**
 * @brief Initialize the simulation context.
 *
 * Checks for the presence of "control.dat" file, reads runtime options, and sets up the user context.
 *
 * @param[out] user    Pointer to the allocated UserCtx structure.
 * @param[out] rank    MPI rank of the process.
 * @param[out] size    Number of MPI processes.
 * @param[out] np      Number of particles.
 * @param[out] rstart  Flag to restart (1) or start from t=0 (0).
 * @param[out] ti      The timestep to start from if restarting.
 * @param[out] nblk    Number of grid blocks.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
 PetscErrorCode InitializeSimulation(UserCtx **user, PetscMPIInt *rank, PetscMPIInt *size, PetscInt *np, PetscInt *rstart, PetscInt *ti, PetscInt *nblk);

/** 
 * @brief Setup grid and vectors for the simulation.
 *
 * @param[in,out] user Pointer to the UserCtx structure.
 * @param[in] block_number The number of grid blocks in the domain.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode SetupGridAndVectors(UserCtx *user, PetscInt block_number);

/**
 * @brief Finalize the simulation and free resources.
 *
 * @param[in,out] user Pointer to the UserCtx structure.
 * @param[in] block_number The number of grid blocks in the domain.
 * @param[in,out] bboxlist  Pointer to the array of BoundingBoxes.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
 PetscErrorCode FinalizeSimulation(UserCtx *user, PetscInt block_number, BoundingBox *bboxlist);

/**
 * @brief Allocates a 3D array of PetscReal values using PetscCalloc.
 *
 * This function dynamically allocates memory for a 3D array of PetscReal values
 * with dimensions nz (layers) x ny (rows) x nx (columns). It uses PetscCalloc1
 * to ensure the memory is zero-initialized.
 *
 * The allocation is done in three steps:
 *  1. Allocate an array of nz pointers (one for each layer).
 *  2. Allocate a contiguous block for nz*ny row pointers and assign each layer’s row pointers.
 *  3. Allocate a contiguous block for all nz*ny*nx PetscReal values.
 *
 * This setup allows the array to be accessed as array[k][j][i], and the memory
 * for the data is contiguous, which improves cache efficiency.
 *
 * @param[out] array Pointer to the 3D array to be allocated.
 * @param[in]  nz    Number of layers (z-direction).
 * @param[in]  ny    Number of rows (y-direction).
 * @param[in]  nx    Number of columns (x-direction).
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
 PetscErrorCode Allocate3DArrayScalar(PetscReal ****array, PetscInt nz, PetscInt ny, PetscInt nx);

/**
 * @brief Deallocates a 3D array of PetscReal values allocated by Allocate3DArrayScalar.
 *
 * This function frees the memory allocated for a 3D array of PetscReal values.
 * It assumes the memory was allocated using Allocate3DArrayScalar, which allocated
 * three separate memory blocks: one for the contiguous data, one for the row pointers,
 * and one for the layer pointers.
 *
 * @param[in] array Pointer to the 3D array to be deallocated.
 * @param[in] nz    Number of layers (z-direction).
 * @param[in] ny    Number of rows (y-direction).
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
PetscErrorCode Deallocate3DArrayScalar(PetscReal ***array, PetscInt nz, PetscInt ny);

/**
 * @brief Deallocates a 3D array of Cmpnts structures allocated by Allocate3DArrayVector.
 *
 * This function frees the memory allocated for a 3D array of Cmpnts structures.
 * It assumes the memory was allocated using Allocate3DArrayVector, which created three
 * separate memory blocks: one for the contiguous vector data, one for the row pointers,
 * and one for the layer pointers.
 *
 * @param[in] array Pointer to the 3D array to be deallocated.
 * @param[in] nz    Number of layers in the z-direction.
 * @param[in] ny    Number of rows in the y-direction.
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
 PetscErrorCode Allocate3DArrayVector(Cmpnts ****array, PetscInt nz, PetscInt ny, PetscInt nx);

/**
 * @brief Deallocates a 3D array of Cmpnts structures allocated by Allocate3DArrayVector.
 *
 * This function frees the memory allocated for a 3D array of Cmpnts structures.
 * It assumes the memory was allocated using Allocate3DArrayVector, which created three
 * separate memory blocks: one for the contiguous vector data, one for the row pointers,
 * and one for the layer pointers.
 *
 * @param[in] array Pointer to the 3D array to be deallocated.
 * @param[in] nz    Number of layers in the z-direction.
 * @param[in] ny    Number of rows in the y-direction.
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
PetscErrorCode Deallocate3DArrayVector(Cmpnts ***array, PetscInt nz, PetscInt ny);

 #endif // SETUP_H
