/**
 * @file grid.h
 * @brief Header file for grid management and coordinate definitions.
 *
 * This file contains declarations of functions responsible for parsing grid inputs,
 * determining grid sizes, initializing DM structures, and assigning coordinates to grid points.
 */

#ifndef GRID_H
#define GRID_H

// Include necessary headers
#include <petsc.h>
#include "common.h"   // Common type definitions
#include "logging.h"  // Logging macros and definitions
#include <stdlib.h>

// --------------------- Function Declarations ---------------------

/**
 * @brief Parses grid inputs from options or files.
 *
 * @param[in,out] user           Pointer to UserCtx structure.
 * @param[out]    generate_grid  Flag indicating whether the grid is programmatically generated (1) or read from a file (0).
 * @param[out]    grid1d         Flag indicating whether the grid is 1D (1) or 3D (0).
 * @param[out]    L_x, L_y, L_z  Domain lengths in x, y, and z directions.
 * @param[out]    imm, jmm, kmm  Grid sizes for blocks in x, y, and z directions.
 * @param[out]    nblk           Number of blocks.
 * @param[in,out] fd             File pointer for grid input file.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ParseGridInputs(UserCtx *user, PetscInt *generate_grid, PetscInt *grid1d,
                               PetscReal *L_x, PetscReal *L_y, PetscReal *L_z,
                               PetscInt *imm, PetscInt *jmm, PetscInt *kmm,
                               PetscInt *nblk, FILE *fd);

/**
 * @brief Determines grid sizes for a specific block.
 *
 * @param[in]     bi            Block index.
 * @param[in,out] user          Pointer to UserCtx structure.
 * @param[out]    IM, JM, KM    Grid dimensions for the block in x, y, and z directions.
 * @param[in,out] fd            File pointer for grid input file.
 * @param[in]     generate_grid Flag indicating whether the grid is programmatically generated (1) or read from a file (0).
 * @param[in]     imm, jmm, kmm Grid sizes for blocks in x, y, and z directions.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode DetermineGridSizes(PetscInt bi, UserCtx *user, PetscInt *IM, PetscInt *JM, PetscInt *KM,
                                  FILE *fd, PetscInt generate_grid, PetscInt *imm, PetscInt *jmm, PetscInt *kmm);

/**
 * @brief Initializes the DM structure for a grid.
 *
 * @param[in,out] user Pointer to UserCtx structure.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeGridDM(UserCtx *user);

/**
 * @brief Assigns coordinates to the grid points.
 *
 * @param[in,out] user          Pointer to UserCtx structure.
 * @param[in]     generate_grid Flag indicating whether the grid is programmatically generated (1) or read from a file (0).
 * @param[in]     grid1d        Flag indicating whether the grid is 1D (1) or 3D (0).
 * @param[in]     IM, JM, KM    Grid dimensions in x, y, and z directions.
 * @param[in]     L_x, L_y, L_z Domain lengths in x, y, and z directions.
 * @param[in]     fd            File pointer for grid input file.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode AssignGridCoordinates(UserCtx *user, PetscInt generate_grid, PetscInt grid1d,
                                     PetscInt IM, PetscInt JM, PetscInt KM,
                                     PetscReal L_x, PetscReal L_y, PetscReal L_z, FILE *fd);

/**
 * @brief Finalizes the grid setup by closing the input file.
 *
 * @param[in] generate_grid Flag indicating whether the grid is programmatically generated (1) or read from a file (0).
 * @param[in] fd            File pointer for grid input file.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode FinalizeGridSetup(PetscInt generate_grid, FILE *fd);

/**
 * @brief Defines grid coordinates for the computational domain.
 *
 * This function orchestrates the grid setup process for all blocks, including parsing inputs,
 * setting up the DM structures, and assigning coordinates.
 *
 * @param[in,out] user Pointer to an array of UserCtx structures, one for each block.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode DefineGridCoordinates(UserCtx *user);

/**
 * @brief Computes the local bounding box of the grid on the current process.
 *
 * This function calculates the minimum and maximum coordinates of the local grid points owned
 * by the current MPI process and stores the computed bounding box in the provided structure.
 *
 * @param[in]  user      Pointer to the user-defined context containing grid information.
 * @param[out] localBBox Pointer to the BoundingBox structure to store the computed bounding box.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ComputeLocalBoundingBox(UserCtx *user, BoundingBox *localBBox);

/**
 * @brief Gathers local bounding boxes from all MPI processes to rank 0.
 *
 * This function computes the local bounding box on each process, then collects all local
 * bounding boxes on the root process (rank 0) using MPI. The result is stored in an array
 * of BoundingBox structures on rank 0.
 *
 * @param[in]  user       Pointer to the user-defined context containing grid information.
 * @param[out] allBBoxes  Pointer to a pointer where the array of gathered bounding boxes
 *                        will be stored on rank 0. The caller on rank 0 must free this array.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode GatherAllBoundingBoxes(UserCtx *user, BoundingBox **allBBoxes);

#endif // GRID_H
