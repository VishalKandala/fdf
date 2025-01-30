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
 * @param[in,out] user      Pointer to the UserCtx structure containing grid details.
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
PetscErrorCode ReadGridGenerationInputs(UserCtx *user, PetscInt *grid1d, PetscReal *L_x, PetscReal *L_y, PetscReal *L_z,
                                        PetscInt **imm, PetscInt **jmm, PetscInt **kmm, PetscInt *nblk) {
    PetscErrorCode ierr;

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadGridGenerationInputs - Reading grid generation parameters.\n");

    // Read the flag indicating if grid is 1D(1) or 3D(0)
    ierr = PetscOptionsGetInt(NULL, NULL, "-grid1d", grid1d, NULL); CHKERRQ(ierr);

    // Set the default stretching ratios
    user->rx = 1.0;
    user->ry = 1.0;
    user->rz = 1.0;

    // Read grid dimensions, number of blocks and stretching ratios.
    ierr = PetscOptionsGetInt(NULL, NULL, "-nblk", nblk, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-L_x", L_x, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-L_y", L_y, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-L_z", L_z, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-r_x", &user->rx, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-r_y", &user->ry, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-r_z", &user->rz, NULL); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadGridGenerationInputs - Number of grid blocks: %d.\n", *nblk);
  
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

    
    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadGridGenerationInputs - Grid parameters read successfully.\n");

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

    
    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadGridFile - Reading grid data from file: %s.\n", filename);

    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);


    if (rank == 0) {
        FILE *fd = fopen(filename, "r");
        if (!fd) {
            SETERRQ1(comm, PETSC_ERR_FILE_OPEN, "Cannot open file: %s", filename);
        }

        // Read number of blocks
        fscanf(fd, "%d\n", nblk);

	LOG_ALLOW(GLOBAL, LOG_INFO, "ReadGridGridFile - Number of grid blocks: %d.\n", *nblk);
  
	// Ensure if 'nblk' is not mentioned in the control.dat file, the default value is set to 1.
	if(*nblk==0) *nblk = 1;

	// Allocate memory for imm, jmm, kmm
	ierr = PetscCalloc1(*nblk, imm); CHKERRQ(ierr);
	ierr = PetscCalloc1(*nblk, jmm); CHKERRQ(ierr);
	ierr = PetscCalloc1(*nblk, kmm); CHKERRQ(ierr);

        // Read block dimensions
        for (PetscInt i = 0; i < *nblk; i++) {
            fscanf(fd, "%" PetscInt_FMT " %" PetscInt_FMT " %" PetscInt_FMT "\n", imm[i], jmm[i], kmm[i]);

	    //	    fscanf(fd, "%" PetscInt_FMT " %" PetscInt_FMT " %" PetscInt_FMT "\n", &imm[i], &jmm[i], &kmm[i]);
	    //  fscanf(fd, "%d %d %d\n", &imm[i], &jmm[i], &kmm[i]);
        }

        // Determine if it's a 1D grid (all dimensions > 1)
	*grid1d = ((*nblk == 1) && (*jmm[0] == 1) && (*kmm[0] == 1)) ? 1 : 0;

	//  *grid1d = ((*nblk == 1) && (jmm[0] == 1) && (kmm[0] == 1)) ? 1 : 0;

        fclose(fd);
    }

    // Broadcast data to all processes
    ierr = MPI_Bcast(nblk, 1, MPI_INT, 0, comm); CHKERRQ(ierr);
    ierr = MPI_Bcast(imm, *nblk, MPI_INT, 0, comm); CHKERRQ(ierr);
    ierr = MPI_Bcast(jmm, *nblk, MPI_INT, 0, comm); CHKERRQ(ierr);
    ierr = MPI_Bcast(kmm, *nblk, MPI_INT, 0, comm); CHKERRQ(ierr);
    ierr = MPI_Bcast(grid1d, 1, MPI_INT, 0, comm); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadGridFile - Blocks: nblk=%d, grid1d=%d.\n", *nblk, *grid1d);
    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadGridFile - Grid file data retrieved successfully.\n");

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

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadFieldData - Attempting to read file: %s\n", filename);

    // Check if the file exists before attempting to open it
    ierr = PetscTestFile(filename, FILE_MODE_READ, &fileExists); CHKERRQ(ierr);
    if (!fileExists) {
    LOG_ALLOW(GLOBAL, LOG_WARNING, "ReadFieldData - File '%s' does not exist.\n", filename);
    char err_msg[512];
    PetscSNPrintf(err_msg, sizeof(err_msg), 
              "ReadFieldData - Could not open file '%s' for reading. File does not exist.", 
              filename);
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, err_msg);
}  
    // Attempt to open the file
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer);
    if (ierr) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "ReadFieldData - Failed to open file: %s. Skipping this field.\n", filename);
        return 0; // Continue execution despite missing file
    }

    // Load data into the vector
    ierr = VecLoad(field_vec, viewer); CHKERRQ(ierr);

    // Close the file viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadFieldData - Successfully loaded data for field: %s\n", field_name);

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadSimulationFields - Starting to read simulation fields.\n");

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadSimulationFields - Finished reading simulation fields.\n");

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadStatisticalFields - Starting to read statistical fields.\n");

    ierr = ReadFieldData(user, "su0", user->Ucat_sum, ti, "dat"); CHKERRQ(ierr);
    ierr = ReadFieldData(user, "su1", user->Ucat_cross_sum, ti, "dat"); CHKERRQ(ierr);
    ierr = ReadFieldData(user, "su2", user->Ucat_square_sum, ti, "dat"); CHKERRQ(ierr);
    ierr = ReadFieldData(user, "sp", user->P_sum, ti, "dat"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadStatisticalFields - Finished reading statistical fields.\n");

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadLESFields - Starting to read LES fields.\n");

    VecDuplicate(user->P, &Cs);
    ierr = ReadFieldData(user, "cs", Cs, ti, "dat"); CHKERRQ(ierr);
    DMGlobalToLocalBegin(user->da, Cs, INSERT_VALUES, user->lCs);
    DMGlobalToLocalEnd(user->da, Cs, INSERT_VALUES, user->lCs);
    VecDestroy(&Cs);

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadLESFields - Finished reading LES fields.\n");

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadRANSFields - Starting to read RANS fields.\n");

    ierr = ReadFieldData(user, "kfield", user->K_Omega, ti, "dat"); CHKERRQ(ierr);
    VecCopy(user->K_Omega, user->K_Omega_o);

    DMGlobalToLocalBegin(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);

    DMGlobalToLocalBegin(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadRANSFields - Finished reading RANS fields.\n");

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

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteFieldData - Attempting to write file: %s\n", filen);

    // Open the file for writing
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);

    // Write data from the vector
    ierr = VecView(field_vec, viewer); CHKERRQ(ierr);

    // Close the file viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteFieldData - Successfully wrote data for field: %s\n", field_name);

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteSimulationFields - Starting to write simulation fields.\n");

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteSimulationFields - Finished writing simulation fields.\n");

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteStatisticalFields - Starting to write statistical fields.\n");

    ierr = WriteFieldData(user, "su0", user->Ucat_sum, user->ti, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "su1", user->Ucat_cross_sum, user->ti, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "su2", user->Ucat_square_sum, user->ti, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "sp", user->P_sum, user->ti, "dat"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteStatisticalFields - Finished writing statistical fields.\n");

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteLESFields - Starting to write LES fields.\n");

    VecDuplicate(user->P, &Cs);
    DMLocalToGlobalBegin(user->da, user->lCs, INSERT_VALUES, Cs);
    DMLocalToGlobalEnd(user->da, user->lCs, INSERT_VALUES, Cs);
    ierr = WriteFieldData(user, "cs", Cs, user->ti, "dat"); CHKERRQ(ierr);
    VecDestroy(&Cs);

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteLESFields - Finished writing LES fields.\n");

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteRANSFields - Starting to write RANS fields.\n");

    ierr = WriteFieldData(user, "kfield", user->K_Omega, user->ti, "dat"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "WriteRANSFields - Finished writing RANS fields.\n");

    return 0;
}

/**
 * @brief Writes data from a specific field in a PETSc Swarm to a file.
 *
 * This function retrieves the Swarm from the UserCtx (i.e., `user->swarm`) and
 * creates a global PETSc vector from the specified Swarm field. It then calls
 * the existing WriteFieldData() function to handle the actual I/O operation.
 * After writing the data, the function destroys the temporary global vector 
 * to avoid memory leaks.
 *
 * @param[in] user       Pointer to the UserCtx structure containing simulation context
 *                       and the PetscSwarm (as `user->swarm`).
 * @param[in] field_name Name of the Swarm field to be written (e.g., "my_field").
 * @param[in] ti         Time index used to construct the output file name.
 * @param[in] ext        File extension (e.g., "dat", "bin").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 *
 * @note Compatible with PETSc 3.14.4.
 */
PetscErrorCode WriteSwarmField(UserCtx *user, const char *field_name, PetscInt ti, const char *ext)
{
  PetscErrorCode ierr;
  Vec            fieldVec;
  DM     swarm;  

  PetscFunctionBegin; /* PETSc macro indicating start of function */

  /* 
   * 1) Retrieve the PetscSwarm from the user context.
   *    Ensure user->swarm is initialized and not NULL.
   */
  swarm = user->swarm;      

  /* 
   * 2) Create a global vector from the specified swarm field.
   *    This function is available in PETSc 3.14.4.
   *    It provides a read/write "view" of the swarm field as a global Vec.
   */
  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "WriteSwarmField - Attempting to create global vector from field: %s\n",
            field_name);
  ierr = DMSwarmCreateGlobalVectorFromField(swarm, field_name, &fieldVec);CHKERRQ(ierr);

  /*
   * 3) Use your existing WriteFieldData() to write the global vector to a file.
   *    The field name, time index, and extension are passed along for naming.
   */
  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "WriteSwarmField - Calling WriteFieldData for field: %s\n",
            field_name);
  ierr = WriteFieldData(user, field_name, fieldVec, ti, ext);CHKERRQ(ierr);

  /*
   * 4) Destroy the global vector once the data is successfully written.
   *    This step is crucial for avoiding memory leaks. 
   *    DMSwarmDestroyGlobalVectorFromField() is also available in PETSc 3.14.4.
   */
  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "WriteSwarmField - Destroying the global vector for field: %s\n",
            field_name);
  ierr = DMSwarmDestroyGlobalVectorFromField(swarm, field_name, &fieldVec);CHKERRQ(ierr);

  /* Log and return success. */
  LOG_ALLOW(GLOBAL, LOG_INFO,
            "WriteSwarmField - Successfully wrote swarm data for field: %s\n",
            field_name);

  PetscFunctionReturn(0); /* PETSc macro indicating end of function */
}

/* 
   -------------- Prototypes of Internal-Use-Only (static) ---------------
   These do not appear in io.h since they are "private" to io.c:
     - WriteVTKAppendedBlock
     - WriteVTSXMLHeader, WriteVTPXMLHeader
     - WriteVTSXMLFooter, WriteVTPXMLFooter
     - WriteVTKFileHeader, WriteVTKFileFooter

static int WriteVTKAppendedBlock(FILE *fp,
                                 const void *buf,
                                 int nvals,
                                 int dsize);

static int WriteVTSXMLHeader(FILE       *fp,
                             int         mx,
                             int         my,
                             int         mz,
                             int         nnodes,
                             const char *fieldName,
                             int         numComponents,
                             int         boffset,
                             int        *boffsetOut);

static int WriteVTSXMLFooter(FILE *fp);

static int WriteVTPXMLHeader(FILE       *fp,
                             int         npoints,
                             const char *fieldName,
                             int         boffset,
                             int        *boffsetOut);

static int WriteVTPXMLFooter(FILE *fp);

static int WriteVTKFileHeader(FILE *fp,
                              const VTKMetaData *meta,
                              int boffset,
                              int *boffsetOut);

static int WriteVTKFileFooter(FILE *fp,
                              const VTKMetaData *meta);


=====================================================================
   1) ReadDataFileToArray

   See the function-level comments in io.h for a summary.
   This reads one value per line from an ASCII file. Rank 0 does I/O,
   broadcasts the data to all ranks.
===================================================================== */
int ReadDataFileToArray(const char   *filename,
                        double      **data_out,
                        int          *Nout,
                        MPI_Comm      comm)
{
    /* STEP 0: Prepare local variables & log function entry */
    int    rank, size;
    PetscErrorCode ierr;
    FILE  *fp = NULL;
    int    N   = 0;            /* number of lines/values read on rank 0 */
    double *array = NULL;      /* pointer to local array on each rank */
    int    fileExistsFlag = 0; /* 0 = doesn't exist, 1 = does exist */

    PetscPrintf(comm, "=== STEP 0: Entering ReadDataFileToArray => file: %s ===\n", filename);
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "ReadDataFileToArray - Start reading from file: %s\n",
              filename);

    /* Basic error checking: data_out, Nout must be non-null. */
    if (!filename || !data_out || !Nout) {
        PetscPrintf(comm,
            "[ERROR] ReadDataFileToArray - Null pointer argument (filename/data_out/Nout).\n");
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "ReadDataFileToArray - Null pointer argument provided.\n");
        return 1;
    }

    /* Determine rank/size for coordinating I/O. */
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* STEP 1: On rank 0, check if file can be opened. */
    PetscPrintf(comm, "STEP 1: Checking file existence on rank %d...\n", rank);
    if (!rank) {
        fp = fopen(filename, "r");
        if (fp) {
            fileExistsFlag = 1;
            fclose(fp);
        }
    }

    /* STEP 2: Broadcast file existence to all ranks. */
    PetscPrintf(comm, "STEP 2: Broadcasting fileExistsFlag (%d) to all %d ranks...\n",
                fileExistsFlag, size);
    // In ReadDataFileToArray:
    ierr = MPI_Bcast(&fileExistsFlag, 1, MPI_INT, 0, comm); CHKERRQ(ierr);

    if (!fileExistsFlag) {
        /* If file does not exist, log & return. */
        if (!rank) {
            PetscPrintf(comm,
                "[ERROR] ReadDataFileToArray - File '%s' not found.\n", filename);
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ReadDataFileToArray - File '%s' not found.\n",
                      filename);
        }
        return 2;
    }

    /* STEP 3: Rank 0 re-opens and reads the file, counting lines, etc. */
    PetscPrintf(comm, "STEP 3: If rank=0, open and parse file: %s\n", filename);
    if (!rank) {
        fp = fopen(filename, "r");
        if (!fp) {
            PetscPrintf(comm,
                "[ERROR] ReadDataFileToArray - Could not reopen file '%s'.\n", filename);
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ReadDataFileToArray - File '%s' could not be opened for reading.\n",
                      filename);
            return 3;
        }

        /* (3a) Count lines first. */
        {
            char line[256];
            while (fgets(line, sizeof(line), fp)) {
                N++;
            }
        }
        PetscPrintf(comm,"ReadDataFileToArray - File '%s' has %d lines.\n",filename, N);

        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "ReadDataFileToArray - File '%s' has %d lines.\n",
                  filename, N);

        /* (3b) Allocate array on rank 0. */
        array = (double*)malloc(N * sizeof(double));
        if (!array) {
            fclose(fp);
            PetscPrintf(comm,
                "[ERROR] ReadDataFileToArray - malloc failed for array.\n");
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ReadDataFileToArray - malloc failed for array.\n");
            return 4;
        }

        /* (3c) Rewind & read values into array. */
        rewind(fp);
        {
            int i = 0;
            char line[256];
            while (fgets(line, sizeof(line), fp)) {
                double val;
                if (sscanf(line, "%lf", &val) == 1) {
                    array[i++] = val;
                }
            }
        }
        fclose(fp);

        PetscPrintf(comm,
            "   => Rank 0 read %d values from '%s'.\n", N, filename);
        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "ReadDataFileToArray - Successfully read %d values from '%s'.\n",
                  N, filename);
    }

    /* STEP 4: Broadcast the integer N to all ranks. */
    PetscPrintf(comm,
        "STEP 4: Broadcasting total count N=%d to all ranks...\n", N);
    ierr = MPI_Bcast(&N, 1, MPI_INT, 0, comm); CHKERRQ(ierr);

    /* STEP 5: Each rank allocates an array to receive the broadcast if rank>0. */
    PetscPrintf(comm,
        "STEP 5: Each rank (including rank 0) ensures array allocated => length=%d.\n", N);
    if (rank) {
        array = (double*)malloc(N * sizeof(double));
        if (!array) {
            PetscPrintf(comm,
                "[ERROR] ReadDataFileToArray - malloc failed on rank %d.\n", rank);
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ReadDataFileToArray - malloc failed on rank %d.\n",
                      rank);
            return 5;
        }
    }

    /* STEP 6: Broadcast the actual data from rank 0 to all. */
    PetscPrintf(comm,
        "STEP 6: Broadcasting data array of length=%d from rank 0 to others...\n", N);
    ierr = MPI_Bcast(array, N, MPI_DOUBLE, 0, comm); CHKERRQ(ierr);

    /* STEP 7: Assign outputs on all ranks. */
    PetscPrintf(comm,
        "STEP 7: Storing array pointer & count in data_out/Nout.\n");
    *data_out = array;
    *Nout     = N;

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "ReadDataFileToArray - Done. Provided array of length=%d to all ranks.\n",
              N);

    PetscPrintf(comm, "=== Exiting ReadDataFileToArray => success, array length=%d ===\n", N);
    return 0; /* success */
}

/**
 * @brief Writes a data block in appended format for a VTK file, including a 4-byte size prefix.
 *
 * This function writes the size of the data block (in bytes) as a 4-byte integer, 
 * followed immediately by the raw data bytes. It returns 0 on success, or 1 on error.
 *
 * @param fp            File pointer to write to (must be open for binary write).
 * @param data          Pointer to the data to be written.
 * @param num_elements  Number of elements in the data array.
 * @param element_size  Size (in bytes) of each element.
 *
 * @return int  Returns 0 if successful, non-zero otherwise.
 */
int WriteVTKAppendedBlock(FILE *fp, const void *data, int num_elements, size_t element_size) {

    // Log the function call with parameters
    LOG_ALLOW_SYNC(LOG_INFO, "WriteVTKAppendedBlock - Called with %d elements, each of size %zu bytes.\n",
                   num_elements, element_size);

    // Calculate the block size
    int block_size = num_elements * (int)element_size;
    LOG_ALLOW_SYNC(LOG_DEBUG, "WriteVTKAppendedBlock - Calculated block size: %d bytes.\n", block_size);

    // Write the block size as a 4-byte integer
    if (fwrite(&block_size, sizeof(int), 1, fp) != 1) {
        fprintf(stderr, "[ERROR] Failed to write block size.\n");
        LOG_ALLOW_SYNC(LOG_ERROR, "WriteVTKAppendedBlock - Error writing block size.\n");
        return 1;
    }

    // Write the actual data
    if (fwrite(data, element_size, num_elements, fp) != (size_t)num_elements) {
        fprintf(stderr, "[ERROR] Failed to write data block.\n");
        LOG_ALLOW_SYNC(LOG_ERROR, "WriteVTKAppendedBlock - Error writing data block.\n");
        return 1;
    }

    // Log success
    LOG_ALLOW_SYNC(LOG_INFO, "WriteVTKAppendedBlock - Successfully wrote block of %d bytes.\n", block_size);

    return 0; // Success
}

/**
 * @brief Writes the XML header portion of a .vts file, defining the structured grid extents and data arrays.
 *
 * This function prints XML tags for the grid, including whole and piece extents,
 * configures the file's byte order and floating-point precision, and sets up
 * appended data sections (for points and a single scalar field). The offsets in
 * the appended data are updated accordingly.
 *
 * @param[in]  fp         File pointer (already open) for writing the .vts header.
 * @param[in]  mx         Number of cells in the x-direction (plus ghost cells).
 * @param[in]  my         Number of cells in the y-direction (plus ghost cells).
 * @param[in]  mz         Number of cells in the z-direction (plus ghost cells).
 * @param[in]  nnodes     Total number of grid nodes.
 * @param[in]  fieldName  Name of the scalar field written to PointData.
 * @param[in]  boffset    Current byte offset in the appended data section.
 * @param[out] boffsetOut Updated byte offset after writing the header information.
 *
 * @return 0 on success, non-zero on failure (though no explicit failure cases are handled here).
 */
static int WriteVTSXMLHeader(FILE       *fp,
                             int         mx,
                             int         my,
                             int         mz,
                             int         nnodes,
                             const char *fieldName,
                             int         boffset,
                             int        *boffsetOut)
{
    // Log the entry into this function and the parameters
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLHeader - Writing .vts header (mx=%d, my=%d, mz=%d, nnodes=%d).\n",
              mx, my, mz, nnodes);

    // Set XML configuration strings
    const char *byte_order = "LittleEndian";
    const char *precision  = "Float64";

    // Begin XML header
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",
            byte_order);
    fprintf(fp, "  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n",
            0, mx-2, 0, my-2, 0, mz-2);
    fprintf(fp, "    <Piece Extent=\"%d %d %d %d %d %d\">\n",
            0, mx-2, 0, my-2, 0, mz-2);

    // Points section
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" "
                "format=\"appended\" offset=\"%d\" />\n",
            precision, boffset);
    boffset += (int)sizeof(int) + 3 * nnodes * (int)sizeof(double);
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLHeader - Updated offset to %d after Points.\n", boffset);
    fprintf(fp, "      </Points>\n");

    // Single scalar field
    fprintf(fp, "      <PointData Scalars=\"%s\">\n", fieldName);
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"1\" "
                "format=\"appended\" offset=\"%d\" />\n",
            precision, fieldName, boffset);
    boffset += (int)sizeof(int) + nnodes * (int)sizeof(double);
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLHeader - Updated offset to %d after PointData.\n", boffset);
    fprintf(fp, "      </PointData>\n");

    // Close the piece and structured grid tags
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </StructuredGrid>\n");

    // Open the appended data section
    fprintf(fp, "  <AppendedData encoding=\"raw\">\n");
    fprintf(fp, "_");

    // Return new offset
    *boffsetOut = boffset;

    // Log successful completion
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLHeader - Completed writing .vts header. New boffset=%d.\n", boffset);

    return 0;
}

/**
 * @brief Finalizes the .vts file by closing the appended data section and the VTKFile XML tag.
 *
 * This function writes the closing XML tags for the appended data section and the 
 * overall VTK file, ensuring a well-formed .vts file structure.
 *
 * @param[in] fp  File pointer (already open) for writing the .vts footer.
 *
 * @return int  Returns 0 on success.
 */
static int WriteVTSXMLFooter(FILE *fp)
{
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLFooter - Closing .vts XML.\n");

    fprintf(fp, "\n </AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");

    // Log completion
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLFooter - Successfully wrote .vts footer.\n");

    return 0;
}


/**
 * @brief Writes the XML header portion of a .vtp file for a point-cloud representation.
 *
 * This function sets up the appended data sections for point coordinates, a named field
 * (e.g., velocity) with a specified number of components, and the connectivity data for
 * each point as a vertex. The byte offset for the appended data is updated after each
 * section, ensuring correct placement of raw binary data in the final output file.
 *
 * @param[in]  fp            File pointer for writing the XML header (already open).
 * @param[in]  npoints       Number of points in the point-cloud.
 * @param[in]  fieldName     Name of the point data field (e.g., "velocity").
 * @param[in]  numComponents Number of components in the fieldName array (e.g., 3 for velocity).
 * @param[in]  boffset       Current byte offset into the appended data section.
 * @param[out] boffsetOut    Updated byte offset after writing the header information.
 *
 * @return int  Returns 0 on success.
 */
static int WriteVTPXMLHeader(FILE *fp,
                             int npoints,
                             const char *fieldName,
                             int numComponents,
                             int boffset,
                             int *boffsetOut)
{
    // Log function entry with parameters
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Called with npoints=%d, fieldName=%s, numComponents=%d, boffset=%d.\n",
              npoints, fieldName, numComponents, boffset);

    // We'll assume double precision => Float64
    const char *precision = "Float64";

    // 1) Standard XML + <VTKFile> opening
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "  <PolyData>\n");

    // For simple point-cloud => #Verts = npoints
    fprintf(fp,
            "    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"%d\" NumberOfLines=\"0\" NumberOfPolys=\"0\" NumberOfStrips=\"0\">\n",
            npoints, npoints);

    // 2) Points => offset=boffset
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" "
                "format=\"appended\" offset=\"%d\"/>\n",
                precision, boffset);
    boffset += sizeof(int) + 3 * npoints * sizeof(double);
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Updated offset to %d after Points.\n", boffset);
    fprintf(fp, "      </Points>\n");

    // 3) PointData => e.g. "velocity" => next offset
    fprintf(fp, "      <PointData>\n");
    fprintf(fp,
            "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%d\" "
            "format=\"appended\" offset=\"%d\"/>\n",
            precision, fieldName, numComponents, boffset);
    boffset += sizeof(int) + numComponents * npoints * sizeof(double);
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Updated offset to %d after PointData.\n", boffset);
    fprintf(fp, "      </PointData>\n");

    // 4) Verts => connectivity + offsets
    fprintf(fp, "      <Verts>\n");
    fprintf(fp,
            "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%d\"/>\n",
            boffset);
    boffset += sizeof(int) + npoints * sizeof(int);
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Updated offset to %d after connectivity.\n", boffset);

    fprintf(fp,
            "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%d\"/>\n",
            boffset);
    boffset += sizeof(int) + npoints * sizeof(int);
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Updated offset to %d after offsets.\n", boffset);
    fprintf(fp, "      </Verts>\n");

    // 5) Close out the piece + PolyData
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </PolyData>\n");

    // 6) <AppendedData encoding="raw"> + underscore
    //    => all binary must follow AFTER this line
    fprintf(fp, "  <AppendedData encoding=\"raw\">\n");
    fprintf(fp, "_");

    *boffsetOut = boffset;

    // Log function completion
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Completed writing VTP XML header. Updated boffset=%d.\n",
              boffset);

    return 0;
}

/**
 * @brief Finalizes the .vtp file by closing the appended data section and the VTKFile tag.
 *
 * This function writes the closing XML tags for the appended data section 
 * and completes the overall VTK structure, ensuring a well-formed .vtp file.
 *
 * @param[in] fp  File pointer to the .vtp file (open for writing).
 *
 * @return int  Returns 0 on success.
 */
static int WriteVTPXMLFooter(FILE *fp)
{
    // Log the function entry
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTPXMLFooter - Closing VTP XML.\n");

    // no stray prints or data -> directly close
    fprintf(fp, "\n  </AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");

    // Log successful completion
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTPXMLFooter - Completed writing VTP XML footer.\n");

    return 0;
}

/**
 * @brief Writes the initial XML header for a VTK file based on the provided metadata.
 *
 * This function currently handles only \c VTK_POLYDATA type. For structured grids or other types,
 * you would add the necessary conditionals and calls. It determines whether to write a scalar or
 * vector field (based on the metadata) and calls \c WriteVTPXMLHeader to construct the appropriate
 * XML header.
 *
 * @param[in]  fp          File pointer (open for writing) to which the header will be written.
 * @param[in]  meta        Pointer to a \c VTKMetaData structure containing fileType, npoints, 
 *                         field names, and other necessary information.
 * @param[in]  boffset     Current byte offset in the appended data section.
 * @param[out] boffsetOut  Updated byte offset after writing any header-related data.
 *
 * @return int  Returns 0 on success, or -1 if the fileType is not handled (e.g., structured grids).
 */
static int WriteVTKFileHeader(FILE *fp, const VTKMetaData *meta, int boffset, int *boffsetOut)
{
    // Log the entry into this function
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTKFileHeader - Entered with fileType=%d, npoints=%d, scalarField=%s, vectorField=%s, numVectorFields=%d.\n",
              meta->fileType, meta->npoints, meta->scalarFieldName, meta->vectorFieldName, meta->numVectorFields);

    // For brevity, only handle VTK_POLYDATA here:
    if (meta->fileType == VTK_POLYDATA) {
        // Decide field name / number of components
        const char *fieldName = NULL;
        int numComponents = 1;

        if (meta->numVectorFields > 0 && meta->vectorFieldName) {
            fieldName = meta->vectorFieldName;
            numComponents = 3; // vector
            LOG_ALLOW(GLOBAL, LOG_DEBUG,
                      "WriteVTKFileHeader - Using vector field '%s' with 3 components.\n", fieldName);
        } else {
            fieldName = meta->scalarFieldName;
            numComponents = 1; // scalar
            LOG_ALLOW(GLOBAL, LOG_DEBUG,
                      "WriteVTKFileHeader - Using scalar field '%s' with 1 component.\n", fieldName);
        }

        // Call the function to write the VTP XML header
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "WriteVTKFileHeader - Calling WriteVTPXMLHeader with boffset=%d.\n", boffset);
        return WriteVTPXMLHeader(fp,
                                 meta->npoints,
                                 fieldName,
                                 numComponents,
                                 boffset,
                                 boffsetOut);
    } else {
        // implement structured later
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "WriteVTKFileHeader - fileType=%d not implemented yet.\n", meta->fileType);
        return -1;
    }
}

/**
 * @brief Writes the XML footer for a VTK file based on the provided metadata.
 *
 * This function currently handles only \c VTK_POLYDATA. If the file type is 
 * \c VTK_POLYDATA, it delegates to \c WriteVTPXMLFooter. Otherwise, 
 * it logs a warning and returns -1 (indicating that other file types are not yet implemented).
 *
 * @param[in] fp    File pointer (open for writing) to which the footer will be written.
 * @param[in] meta  Pointer to a \c VTKMetaData structure containing the fileType 
 *                  (e.g., \c VTK_POLYDATA or \c VTK_STRUCTURED).
 *
 * @return int  Returns 0 or the value from \c WriteVTPXMLFooter on success, 
 *              or -1 if the file type is not handled.
 */
static int WriteVTKFileFooter(FILE *fp, const VTKMetaData *meta)
{
    // Log the entry into this function
    LOG_ALLOW(GLOBAL, LOG_DEBUG, 
              "WriteVTKFileFooter - Called with fileType=%d.\n", meta->fileType);

    if (meta->fileType == VTK_POLYDATA) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTKFileFooter - Writing VTP XML footer.\n");
        return WriteVTPXMLFooter(fp);
    } else {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "WriteVTKFileFooter - fileType=%d not implemented, returning -1.\n", meta->fileType);
        fprintf(stderr, "VTK_STRUCTURED path not shown\n");
        return -1;
    }
}

/**
 * @brief Creates and writes a VTK file (either .vts or .vtp) based on the provided metadata.
 *
 * This function gathers the necessary information from the \c VTKMetaData structure (e.g., \c coords,
 * \c scalarField, \c vectorField, and connectivity or offsets if it is \c VTK_POLYDATA). It writes
 * the XML header, followed by the appended binary data blocks in the correct order, and finally the
 * XML footer. The I/O occurs only on \c rank 0 of the provided \c MPI_Comm.
 *
 * @param[in]  filename  The output file name (e.g., "output.vtp" or "output.vts").
 * @param[in]  meta      Pointer to a \c VTKMetaData structure containing all necessary fields.
 * @param[in]  comm      The MPI communicator used for parallel execution.
 *
 * @return int  Returns 0 on success, or 1 on file-opening failure. Always returns 0 on non-\c rank 0 processes.
 *
 * @note
 * - Currently, only \c VTK_POLYDATA is demonstrated. 
 * - For structured grids (\c VTK_STRUCTURED), further implementation is needed in the header and footer logic.
 */
int CreateVTKFileFromMetadata(const char *filename,
                              const VTKMetaData *meta,
                              MPI_Comm comm)
{
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Entry point. Filename: %s\n", filename);

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - rank=%d, size=%d.\n", rank, size);

    // 1) Only rank 0 writes the file
    if (!rank) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "CreateVTKFileFromMetadata - Rank 0 writing file '%s'.\n", filename);

        FILE *fp = fopen(filename, "wb");
        if (!fp) {
            PetscPrintf(comm, "[ERROR] Could not open '%s' for writing.\n", filename);
            LOG_ALLOW(GLOBAL, LOG_ERROR, "CreateVTKFileFromMetadata - fopen failed for %s.\n", filename);
            return 1;
        }
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Successfully opened file: %s\n", filename);

        int boffset = 0;
        // 2) Write the XML header => sets <DataArray ... offset="...">
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing header (initial boffset=%d).\n", boffset);
        WriteVTKFileHeader(fp, meta, boffset, &boffset);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Header written (updated boffset=%d).\n", boffset);

        // 3) Immediately write appended data blocks in EXACT order:
        //    (a) coords (3*npoints doubles)
        if (meta->coords) {
            int ncoords = 3 * meta->npoints;
            if (ncoords > 0) {
                LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing coords block: %d doubles.\n", ncoords);
                WriteVTKAppendedBlock(fp, meta->coords, ncoords, sizeof(double));
                boffset += sizeof(int) + ncoords * sizeof(double);
            }
        }

        //    (b) scalar or vector field
        if (meta->numScalarFields > 0 && meta->scalarField) {
            int nvals = meta->npoints; // 1 component
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing scalar field block: %d doubles.\n", nvals);
            WriteVTKAppendedBlock(fp, meta->scalarField, nvals, sizeof(double));
            boffset += sizeof(int) + nvals * sizeof(double);
        }
        else if (meta->numVectorFields > 0 && meta->vectorField) {
            int nvals = 3 * meta->npoints; // 3 comps
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing vector field block: %d doubles.\n", nvals);
            WriteVTKAppendedBlock(fp, meta->vectorField, nvals, sizeof(double));
            boffset += sizeof(int) + nvals * sizeof(double);
        }

        //    (c) connectivity (if polydata)
        if (meta->fileType == VTK_POLYDATA && meta->connectivity) {
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing connectivity block: %d ints.\n", meta->npoints);
            WriteVTKAppendedBlock(fp, meta->connectivity, meta->npoints, sizeof(int));
            boffset += sizeof(int) + meta->npoints * sizeof(int);
        }

        //    (d) offsets (if polydata)
        if (meta->fileType == VTK_POLYDATA && meta->offsets) {
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing offsets block: %d ints.\n", meta->npoints);
            WriteVTKAppendedBlock(fp, meta->offsets, meta->npoints, sizeof(int));
            boffset += sizeof(int) + meta->npoints * sizeof(int);
        }

        // 4) Write footer => </AppendedData> & </VTKFile>
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing footer.\n");
        WriteVTKFileFooter(fp, meta);

        fclose(fp);
        LOG_ALLOW(GLOBAL, LOG_INFO, "CreateVTKFileFromMetadata - File '%s' closed.\n", filename);
    }

    return 0;
}
