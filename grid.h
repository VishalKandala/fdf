#ifndef GRID_H
#define GRID_H

#include <petsc.h>
#include "logging.h"
#include "variables.h"
#include <stdlib.h>


// Function Declarations

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
 * @param[in,out] user         Pointer to UserCtx structure.
 * @param[in]     generate_grid Flag indicating whether the grid is programmatically generated (1) or read from a file (0).
 * @param[in]     grid1d       Flag indicating whether the grid is 1D (1) or 3D (0).
 * @param[in]     IM, JM, KM   Grid dimensions in x, y, and z directions.
 * @param[in]     L_x, L_y, L_z Domain lengths in x, y, and z directions.
 * @param[in]     fd           File pointer for grid input file.
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

#endif // GRID_H
