/**
 * @file io.h
 * @brief Header file for input/output operations in the simulation.
 *
 * This file declares functions for reading and writing simulation fields, grid generation parameters,
 * and optional statistical, LES, and RANS data.
 */

#ifndef IO_H
#define IO_H

#include <petsc.h> // System dependency, required.
#include "common.h" // For shared definitions like UserCtx.

// ------------------------ Function Declarations ------------------------

/**
 * @brief Reads grid generation parameters from PETSc options or configuration files.
 *
 * Reads domain dimensions and block parameters for programmatically generating a grid.
 *
 * @param[out] grid1d       Pointer to flag indicating if the grid is 1D (1) or 3D (0).
 * @param[out] L_x          Pointer to domain length in the x-direction.
 * @param[out] L_y          Pointer to domain length in the y-direction.
 * @param[out] L_z          Pointer to domain length in the z-direction.
 * @param[out] imm          Pointer to array storing the i-dimensions of each block.
 * @param[out] jmm          Pointer to array storing the j-dimensions of each block.
 * @param[out] kmm          Pointer to array storing the k-dimensions of each block.
 * @param[out] nblk         Pointer to number of blocks.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadGridGenerationInputs(PetscInt *grid1d, PetscReal *L_x, PetscReal *L_y, PetscReal *L_z,
                                        PetscInt **imm, PetscInt **jmm, PetscInt **kmm, PetscInt *nblk);
/**
 * @brief Reads grid parameters from a file.
 *
 * This function retrieves grid block dimensions, number of blocks, and grid type (1D or 3D)
 * from a file and broadcasts the data to all MPI processes.
 *
 * @param[in]  filename  Name of the file containing grid data.
 * @param[out] nblk      Pointer to store the number of blocks.
 * @param[out] imm       Pointer to Array to store the i-dimensions for each block.
 * @param[out] jmm       Pointer to Array to store the j-dimensions for each block.
 * @param[out] kmm       Pointer to Array to store the k-dimensions for each block.
 * @param[out] grid1d    Pointer to store the grid type (1 for 1D, 0 for 3D).
 * @param[in]  comm      MPI communicator for broadcasting data.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode ReadGridFile(const char *filename, PetscInt *nblk, PetscInt **imm, PetscInt **jmm, 
                            PetscInt **kmm, PetscInt *grid1d, MPI_Comm comm);

/**
 * @brief Reads binary field data for velocity, pressure, and other required vectors.
 *
 * Reads contravariant velocity (`Ucont`) from `vfield`, Cartesian velocity (`Ucat`) from `ufield`,
 * pressure (`P`), node state (`Nvert_o`), and optionally statistical quantities, LES, and RANS fields
 * from binary files. Logs missing files but continues execution.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the simulation context.
 * @param[in]     The timestep at which the simulation field data needs to be read.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadSimulationFields(UserCtx *user, PetscInt ti);


/**
 * @brief Reads data for a specific field from a file into the provided vector.
 *
 * This function uses the field name to construct the file path and reads the data
 * from the corresponding file into the provided PETSc vector.
 *
 * @param[in]     user        Pointer to the UserCtx structure containing simulation context.
 * @param[in]     field_name  Name of the field (e.g., "ufield", "vfield", "pfield").
 * @param[out]    field_vec   PETSc vector to store the field data.
 * @param[in]     ti          Time index for constructing the file name.
 * @param[in]     ext         File extension (e.g., "dat").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadFieldData(UserCtx *user, const char *field_name, Vec field_vec, PetscInt ti, const char *ext);


/**
 * @brief Reads statistical fields used for time-averaged simulations.
 *
 * Reads statistical quantities such as velocity sums and pressure sums. Logs missing files and initializes
 * fields to zero if files are unavailable.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the simulation context.
 * @param[in]     The timestep at which the simulation field data needs to be read.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadStatisticalFields(UserCtx *user,PetscInt ti);

/**
 * @brief Reads LES-related fields used in turbulence modeling.
 *
 * Reads the Smagorinsky constant (`Cs`) and transfers data to local vectors. Logs missing files and
 * initializes fields to zero if files are unavailable.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the simulation context.
 * @param[in]     ti          Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadLESFields(UserCtx *user, PetscInt ti);

/**
 * @brief Reads RANS-related fields for turbulence modeling.
 *
 * Reads `K_Omega` fields (used in RANS modeling) and initializes them if files are unavailable.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the simulation context.
 * @param[in]     ti          Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadRANSFields(UserCtx *user, PetscInt ti);

/**
 * @brief Writes data from a specific PETSc vector to a file.
 *
 * This function uses the field name to construct the file path and writes the data
 * from the provided PETSc vector to the corresponding file.
 *
 * @param[in] user       Pointer to the UserCtx structure containing simulation context.
 * @param[in] field_name Name of the field (e.g., "ufield", "vfield", "pfield").
 * @param[in] field_vec  PETSc vector containing the field data to write.
 * @param[in] ti         Time index for constructing the file name.
 * @param[in] ext        File extension (e.g., "dat").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteFieldData(UserCtx *user, const char *field_name, Vec field_vec, PetscInt ti, const char *ext);


/**
 * @brief Writes simulation fields to files.
 *
 * This function writes contravariant velocity, Cartesian velocity, pressure, and node state
 * fields to their respective binary files. It also conditionally writes LES, RANS, and
 * statistical fields if they are enabled.
 *
 * @param[in] user Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteSimulationFields(UserCtx *user);


/**
 * @brief Writes statistical fields for averaging purposes.
 *
 * This function writes data for fields such as Ucat_sum, Ucat_cross_sum, Ucat_square_sum,
 * and P_sum to their respective binary files.
 *
 * @param[in] user Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteStatisticalFields(UserCtx *user);

/**
 * @brief Writes LES-related fields.
 *
 * This function writes LES-related fields such as Cs (Smagorinsky constant)
 * to their respective binary files.
 *
 * @param[in] user Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteLESFields(UserCtx *user);

/**
 * @brief Writes RANS-related fields.
 *
 * This function writes RANS-related fields such as K_Omega to their respective
 * binary files.
 *
 * @param[in] user Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteRANSFields(UserCtx *user);

#endif // IO_H
