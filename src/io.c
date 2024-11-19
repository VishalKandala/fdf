/**
 * @file io.c
 * @brief Handles input/output operations for velocity, pressure, grid, and other simulation fields.
 *
 * This file provides functions for reading and writing simulation fields to binary files,
 * including optional handling for statistical, LES, and RANS data.
 */

#include <petsc.h>          // System dependency.
#include "common.h"         // For UserCtx and shared definitions.
#include "logging.h"        // For logging macros.
#include "grid.h"           // Only if grid-related functions are directly used.
#include "io.h"             

// ------------------------ Function Definitions ------------------------

/**
 * @brief Reads grid generation parameters from input options or configuration.
 *
 * This function reads domain dimensions and block parameters for programmatically generating a grid.
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
                                        PetscInt **imm, PetscInt **jmm, PetscInt **kmm, PetscInt *nblk) {
    PetscErrorCode ierr;

    LOG(GLOBAL, LOG_INFO, "ReadGridGenerationInputs - Reading grid generation parameters.\n");

    // Read the flag indicating if grid is 1D(1) or 3D(0)
    ierr = PetscOptionsGetInt(NULL, NULL, "-grid1d", grid1d, NULL); CHKERRQ(ierr);

    // Read grid dimensions & number of blocks
    ierr = PetscOptionsGetInt(NULL, NULL, "-nblk", nblk, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-L_x", L_x, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-L_y", L_y, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-L_z", L_z, NULL); CHKERRQ(ierr);
    
    LOG(GLOBAL, LOG_INFO, "ReadGridGenerationInputs - Number of grid blocks: %d.\n", *nblk);
  
    // Ensure if 'nblk' is not mentioned in the control.dat file, the default value is set to 1.
    if(*nblk==0) *nblk = 1;

    // Allocate memory for imm, jmm, kmm
    ierr = PetscCalloc1(*nblk, imm); CHKERRQ(ierr);
    ierr = PetscCalloc1(*nblk, jmm); CHKERRQ(ierr);
    ierr = PetscCalloc1(*nblk, kmm); CHKERRQ(ierr);

    // Read the number of grid points for each block
    ierr = PetscOptionsGetIntArray(NULL, NULL, "-im", *imm, nblk, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetIntArray(NULL, NULL, "-jm", *jmm, nblk, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetIntArray(NULL, NULL, "-km", *kmm, nblk, NULL); CHKERRQ(ierr);

    LOG(GLOBAL, LOG_INFO, "ReadGridGenerationInputs - Grid parameters read successfully.\n");

    return 0;
}


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
                            PetscInt **kmm, PetscInt *grid1d, MPI_Comm comm) {
    PetscErrorCode ierr;
    PetscInt rank;

    
    LOG(GLOBAL, LOG_INFO, "ReadGridFile - Reading grid data from file: %s.\n", filename);

    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);


    if (rank == 0) {
        FILE *fd = fopen(filename, "r");
        if (!fd) {
            SETERRQ1(comm, PETSC_ERR_FILE_OPEN, "Cannot open file: %s", filename);
        }

        // Read number of blocks
        fscanf(fd, "%d\n", nblk);

	LOG(GLOBAL, LOG_INFO, "ReadGridGridFile - Number of grid blocks: %d.\n", *nblk);
  
	// Ensure if 'nblk' is not mentioned in the control.dat file, the default value is set to 1.
	if(*nblk==0) *nblk = 1;

	// Allocate memory for imm, jmm, kmm
	ierr = PetscCalloc1(*nblk, imm); CHKERRQ(ierr);
	ierr = PetscCalloc1(*nblk, jmm); CHKERRQ(ierr);
	ierr = PetscCalloc1(*nblk, kmm); CHKERRQ(ierr);

        // Read block dimensions
        for (PetscInt i = 0; i < *nblk; i++) {
            fscanf(fd, "%d %d %d\n", &imm[i], &jmm[i], &kmm[i]);
        }

        // Determine if it's a 1D grid (all dimensions > 1)
        *grid1d = ((*nblk == 1) && (jmm[0] == 1) && (kmm[0] == 1)) ? 1 : 0;

        fclose(fd);
    }

    // Broadcast data to all processes
    ierr = MPI_Bcast(nblk, 1, MPI_INT, 0, comm); CHKERRQ(ierr);
    ierr = MPI_Bcast(imm, *nblk, MPI_INT, 0, comm); CHKERRQ(ierr);
    ierr = MPI_Bcast(jmm, *nblk, MPI_INT, 0, comm); CHKERRQ(ierr);
    ierr = MPI_Bcast(kmm, *nblk, MPI_INT, 0, comm); CHKERRQ(ierr);
    ierr = MPI_Bcast(grid1d, 1, MPI_INT, 0, comm); CHKERRQ(ierr);

    LOG(GLOBAL, LOG_DEBUG, "ReadGridFile - Blocks: nblk=%d, grid1d=%d.\n", *nblk, *grid1d);
    LOG(GLOBAL, LOG_INFO, "ReadGridFile - Grid file data retrieved successfully.\n");

    return 0;
}

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
PetscErrorCode ReadFieldData(UserCtx *user, const char *field_name, Vec field_vec, PetscInt ti, const char *ext)
{
    PetscErrorCode ierr;   // PETSc error handling
    PetscViewer viewer;    // PETSc binary viewer for reading files
    char filename[256];    // File name buffer
    PetscBool fileExists;  // Flag to check file existence
    PetscInt rank;         // MPI rank

    // Validate inputs
    if (!user || !field_name || !field_vec || !ext) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, 
                "ReadFieldData - Null argument provided (user, field, vec, or extension).");
    }

    // Retrieve MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Construct the file name
    user->_this = 0; // Bad hack for now! until the use of _this is figured out! 
    snprintf(filename, sizeof(filename), "results/%s%05d_%d.%s", field_name, ti, user->_this, ext);

    LOG(GLOBAL, LOG_DEBUG, "ReadFieldData - Attempting to read file: %s\n", filename);

    // Check if the file exists before attempting to open it
    ierr = PetscTestFile(filename, FILE_MODE_READ, &fileExists); CHKERRQ(ierr);
    if (!fileExists) {
        LOG(GLOBAL, LOG_WARNING, "ReadFieldData - File '%s' does not exist. Skipping.\n", filename);
        return PETSC_ERR_FILE_OPEN; // File does not exist
    }

    // Attempt to open the file
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer);
    if (ierr) {
        LOG(GLOBAL, LOG_WARNING, "ReadFieldData - Failed to open file: %s. Skipping this field.\n", filename);
        return 0; // Continue execution despite missing file
    }

    // Load data into the vector
    ierr = VecLoad(field_vec, viewer); CHKERRQ(ierr);

    // Close the file viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    LOG(GLOBAL, LOG_INFO, "ReadFieldData - Successfully loaded data for field: %s\n", field_name);

    return 0;
}

/**
 * @brief Reads simulation fields from files into their respective PETSc vectors.
 *
 * This function reads contravariant velocity, Cartesian velocity, pressure, and node state
 * fields from their respective binary files. It also conditionally reads LES, RANS, and
 * statistical fields if they are enabled.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing simulation context.
 * @param[in]     ti          Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadSimulationFields(UserCtx *user,PetscInt ti)
{
    PetscErrorCode ierr;

    LOG(GLOBAL, LOG_INFO, "ReadSimulationFields - Starting to read simulation fields.\n");

    // Read Cartesian velocity field
    ierr = ReadFieldData(user, "ufield", user->Ucat, ti, "dat"); CHKERRQ(ierr);

    // Read contravariant velocity field
    ierr = ReadFieldData(user, "vfield", user->Ucont, ti, "dat"); CHKERRQ(ierr);

    // Read pressure field
    ierr = ReadFieldData(user, "pfield", user->P, ti, "dat"); CHKERRQ(ierr);

    // Read node state field (nvert)
    ierr = ReadFieldData(user, "nvfield", user->Nvert_o, ti, "dat"); CHKERRQ(ierr);

    // Process LES fields if enabled
    if (user->les) {
      ierr = ReadLESFields(user,ti); CHKERRQ(ierr);
    }

    // Process RANS fields if enabled
    if (user->rans) {
      ierr = ReadRANSFields(user,ti); CHKERRQ(ierr);
    }

    // Process statistical fields if averaging is enabled
    if (user->averaging) {
      ierr = ReadStatisticalFields(user,ti); CHKERRQ(ierr);
    }

    LOG(GLOBAL, LOG_INFO, "ReadSimulationFields - Finished reading simulation fields.\n");

    return 0;
}

/**
 * @brief Reads statistical fields for averaging purposes.
 *
 * This function reads data for fields such as Ucat_sum, Ucat_cross_sum, Ucat_square_sum,
 * and P_sum, used for statistical analysis during simulation.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing simulation context.
 * @param[in]     ti          Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadStatisticalFields(UserCtx *user,PetscInt ti)
{
    PetscErrorCode ierr;

    LOG(GLOBAL, LOG_INFO, "ReadStatisticalFields - Starting to read statistical fields.\n");

    ierr = ReadFieldData(user, "su0", user->Ucat_sum, ti, "dat"); CHKERRQ(ierr);
    ierr = ReadFieldData(user, "su1", user->Ucat_cross_sum, ti, "dat"); CHKERRQ(ierr);
    ierr = ReadFieldData(user, "su2", user->Ucat_square_sum, ti, "dat"); CHKERRQ(ierr);
    ierr = ReadFieldData(user, "sp", user->P_sum, ti, "dat"); CHKERRQ(ierr);

    LOG(GLOBAL, LOG_INFO, "ReadStatisticalFields - Finished reading statistical fields.\n");

    return 0;
}

/**
 * @brief Reads LES-related fields.
 *
 * This function reads LES-related fields such as Cs (Smagorinsky constant)
 * into their respective PETSc vectors.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing simulation context.
 * @param[in]     ti          Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadLESFields(UserCtx *user,PetscInt ti)
{
    PetscErrorCode ierr;
    Vec Cs;

    LOG(GLOBAL, LOG_INFO, "ReadLESFields - Starting to read LES fields.\n");

    VecDuplicate(user->P, &Cs);
    ierr = ReadFieldData(user, "cs", Cs, ti, "dat"); CHKERRQ(ierr);
    DMGlobalToLocalBegin(user->da, Cs, INSERT_VALUES, user->lCs);
    DMGlobalToLocalEnd(user->da, Cs, INSERT_VALUES, user->lCs);
    VecDestroy(&Cs);

    LOG(GLOBAL, LOG_INFO, "ReadLESFields - Finished reading LES fields.\n");

    return 0;
}

/**
 * @brief Reads RANS-related fields.
 *
 * This function reads RANS-related fields such as K_Omega into their respective
 * PETSc vectors.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing simulation context.
 * @param[in]     ti          Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadRANSFields(UserCtx *user,PetscInt ti)
{
    PetscErrorCode ierr;

    LOG(GLOBAL, LOG_INFO, "ReadRANSFields - Starting to read RANS fields.\n");

    ierr = ReadFieldData(user, "kfield", user->K_Omega, ti, "dat"); CHKERRQ(ierr);
    VecCopy(user->K_Omega, user->K_Omega_o);

    DMGlobalToLocalBegin(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);

    DMGlobalToLocalBegin(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);

    LOG(GLOBAL, LOG_INFO, "ReadRANSFields - Finished reading RANS fields.\n");

    return 0;
}

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
PetscErrorCode WriteFieldData(UserCtx *user, const char *field_name, Vec field_vec, PetscInt ti, const char *ext)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    char filen[128];

    // Construct the file name
    snprintf(filen, sizeof(filen), "results/%s%05d_%d.%s", field_name, ti, user->_this, ext);

    LOG(GLOBAL, LOG_DEBUG, "WriteFieldData - Attempting to write file: %s\n", filen);

    // Open the file for writing
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);

    // Write data from the vector
    ierr = VecView(field_vec, viewer); CHKERRQ(ierr);

    // Close the file viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    LOG(GLOBAL, LOG_INFO, "WriteFieldData - Successfully wrote data for field: %s\n", field_name);

    return 0;
}

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
PetscErrorCode WriteSimulationFields(UserCtx *user)
{
    PetscErrorCode ierr;

    LOG(GLOBAL, LOG_INFO, "WriteSimulationFields - Starting to write simulation fields.\n");

    // Write contravariant velocity field
    ierr = WriteFieldData(user, "vfield", user->Ucont, user->ti, "dat"); CHKERRQ(ierr);

    // Write Cartesian velocity field
    ierr = WriteFieldData(user, "ufield", user->Ucat, user->ti, "dat"); CHKERRQ(ierr);

    // Write pressure field
    ierr = WriteFieldData(user, "pfield", user->P, user->ti, "dat"); CHKERRQ(ierr);

    // Write node state field (nvert)
    ierr = WriteFieldData(user, "nvfield", user->Nvert, user->ti, "dat"); CHKERRQ(ierr);

    // Write LES fields if enabled
    if (user->les) {
        ierr = WriteLESFields(user); CHKERRQ(ierr);
    }

    // Write RANS fields if enabled
    if (user->rans) {
        ierr = WriteRANSFields(user); CHKERRQ(ierr);
    }

    // Write statistical fields if averaging is enabled
    if (user->averaging) {
        ierr = WriteStatisticalFields(user); CHKERRQ(ierr);
    }

    LOG(GLOBAL, LOG_INFO, "WriteSimulationFields - Finished writing simulation fields.\n");

    return 0;
}

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
PetscErrorCode WriteStatisticalFields(UserCtx *user)
{
    PetscErrorCode ierr;

    LOG(GLOBAL, LOG_INFO, "WriteStatisticalFields - Starting to write statistical fields.\n");

    ierr = WriteFieldData(user, "su0", user->Ucat_sum, user->ti, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "su1", user->Ucat_cross_sum, user->ti, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "su2", user->Ucat_square_sum, user->ti, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "sp", user->P_sum, user->ti, "dat"); CHKERRQ(ierr);

    LOG(GLOBAL, LOG_INFO, "WriteStatisticalFields - Finished writing statistical fields.\n");

    return 0;
}

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
PetscErrorCode WriteLESFields(UserCtx *user)
{
    PetscErrorCode ierr;
    Vec Cs;

    LOG(GLOBAL, LOG_INFO, "WriteLESFields - Starting to write LES fields.\n");

    VecDuplicate(user->P, &Cs);
    DMLocalToGlobalBegin(user->da, user->lCs, INSERT_VALUES, Cs);
    DMLocalToGlobalEnd(user->da, user->lCs, INSERT_VALUES, Cs);
    ierr = WriteFieldData(user, "cs", Cs, user->ti, "dat"); CHKERRQ(ierr);
    VecDestroy(&Cs);

    LOG(GLOBAL, LOG_INFO, "WriteLESFields - Finished writing LES fields.\n");

    return 0;
}

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
PetscErrorCode WriteRANSFields(UserCtx *user)
{
    PetscErrorCode ierr;

    LOG(GLOBAL, LOG_INFO, "WriteRANSFields - Starting to write RANS fields.\n");

    ierr = WriteFieldData(user, "kfield", user->K_Omega, user->ti, "dat"); CHKERRQ(ierr);

    LOG(GLOBAL, LOG_INFO, "WriteRANSFields - Finished writing RANS fields.\n");

    return 0;
}






