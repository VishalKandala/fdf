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
*/

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


/* =====================================================================
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
    int         rank, size;
    FILE       *fp = NULL;
    int         N   = 0;          /* number of lines/values read on rank 0 */
    double     *array = NULL;     /* pointer to local array (allocated on each rank) */
    int         fileExistsFlag = 0; /* 0 = doesn't exist, 1 = does exist */

    /* Log function entry (optional) */
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "ReadDataFileToArray - Start reading from file: %s\n",
              filename);

    /* Basic error checking: data_out, Nout must be non-null. */
    if (!filename || !data_out || !Nout) {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "ReadDataFileToArray - Null pointer argument provided.\n");
        return 1; /* or use a PETSc error code, e.g. SETERRQ(...) if integrated with PETSc */
    }

    /* Get rank and size for coordinating I/O. */
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* 1) On rank 0, try to open the file to see if it exists. */
    if (!rank) {
        fp = fopen(filename, "r");
        if (fp) {
            fileExistsFlag = 1;
            fclose(fp);
        }
    }

    /* 2) Broadcast existence to all ranks. */
    MPI_Bcast(&fileExistsFlag, 1, MPI_INT, 0, comm);

    if (!fileExistsFlag) {
        /* If file does not exist, log a warning and return an error. */
        if (!rank) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ReadDataFileToArray - File '%s' not found.\n",
                      filename);
        }
        return 2; /* some error code. If integrated with PETSc, you might do SETERRQ(...) */
    }

    /* 3) Now rank 0 re-opens the file and reads all lines. */
    if (!rank) {
        fp = fopen(filename, "r");
        if (!fp) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ReadDataFileToArray - File '%s' could not be opened for reading.\n",
                      filename);
            return 3;
        }

        /* a) Count lines first. */
        {
            char line[256];
            while (fgets(line, sizeof(line), fp)) {
                N++;
            }
        }
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "ReadDataFileToArray - File '%s' has %d lines.\n",
                  filename, N);

        /* b) Allocate array. We use double here; if you prefer "PetscScalar",
              adapt accordingly. */
        /* We'll broadcast the array to all ranks, so rank 0 keeps it. */
        array = (double*)malloc(N * sizeof(double));
        if (!array) {
            fclose(fp);
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ReadDataFileToArray - malloc failed for array.\n");
            return 4;
        }

        /* c) Rewind file and parse each line into array. */
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

        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "ReadDataFileToArray - Successfully read %d values from '%s'.\n",
                  N, filename);
    }

    /* 4) Broadcast the integer N to all ranks so they can allocate the array too. */
    MPI_Bcast(&N, 1, MPI_INT, 0, comm);

    /* 5) Each rank allocates an array to receive the broadcast. */
    if (rank) {
        array = (double*)malloc(N * sizeof(double));
        if (!array) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ReadDataFileToArray - malloc failed on rank %d.\n",
                      rank);
            return 5;
        }
    }

    /* 6) Broadcast the actual data from rank 0 to others. */
    MPI_Bcast(array, N, MPI_DOUBLE, 0, comm);

    /* 7) Assign outputs. */
    *data_out = array;
    *Nout     = N;

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "ReadDataFileToArray - Done. Provided array of length=%d to all ranks.\n",
              N);

    return 0; /* success */
}

/* =====================================================================
   2) VTK Writing: Private Helper Functions

   We define a set of static (file-scope) functions for writing:
   - block data (size + raw bytes)
   - .vts headers/footers
   - .vtp headers/footers
   - combined wrappers (header, footer)

   These are used by CreateVTKFileFromMetadata (the public function).
===================================================================== */

/* -- WriteVTKAppendedBlock -----------------------------------------
   Writes one "appended data" block:
     1) An integer: blockSize = (nvals * dsize) bytes
     2) The raw data bytes themselves (buf).

   Example usage:
     int blockSize = nvals * dsize;
     fwrite(&blockSize, sizeof(int), 1, fp);
     fwrite(buf, dsize, nvals, fp);
-------------------------------------------------------------------- */
static int WriteVTKAppendedBlock(FILE       *fp,
                                 const void *buf,
                                 int         nvals,
                                 int         dsize)
{
    int blockSize = nvals * dsize;

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTKAppendedBlock - Writing block with %d elems, blockSize=%d bytes.\n",
              nvals, blockSize);

    /* 1) Write the blockSize as a 4-byte integer. */
    fwrite(&blockSize, sizeof(int), 1, fp);

    /* 2) Write the raw data array. */
    fwrite(buf, dsize, (size_t)nvals, fp);

    return 0; /* success */
}

/* -- WriteVTSXMLHeader ---------------------------------------------
   Writes the XML tags for a .vts file:
     <VTKFile type="StructuredGrid" ...>
       <StructuredGrid WholeExtent="...">
         <Piece Extent="...">
           <Points>
             <DataArray ... offset="..." />
           </Points>
           <PointData>
             <DataArray ... offset="..." />
           </PointData>
         </Piece>
       </StructuredGrid>
       <AppendedData encoding="raw">
       _
   Also updates boffset to reflect the appended data blocks that follow.
-------------------------------------------------------------------- */
static int WriteVTSXMLHeader(FILE       *fp,
                             int         mx,
                             int         my,
                             int         mz,
                             int         nnodes,
                             const char *fieldName,
                             int         boffset,
                             int        *boffsetOut)
{
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLHeader - Writing .vts header (mx=%d, my=%d, mz=%d, nnodes=%d).\n",
              mx, my, mz, nnodes);

    const char *byte_order = "LittleEndian";
    const char *precision  = "Float64";

    /* XML header part */
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",
            byte_order);
    fprintf(fp, "  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n",
            0, mx-2, 0, my-2, 0, mz-2);
    fprintf(fp, "    <Piece Extent=\"%d %d %d %d %d %d\">\n",
            0, mx-2, 0, my-2, 0, mz-2);

    /* Points => offset1 */
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" "
                "format=\"appended\" offset=\"%d\" />\n",
            precision, boffset);
    boffset += (int)sizeof(int) + 3 * nnodes * (int)sizeof(double);
    fprintf(fp, "      </Points>\n");

    /* Single scalar field => offset2 */
    fprintf(fp, "      <PointData Scalars=\"%s\">\n", fieldName);
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"1\" "
                "format=\"appended\" offset=\"%d\" />\n",
            precision, fieldName, boffset);
    boffset += (int)sizeof(int) + nnodes * (int)sizeof(double);
    fprintf(fp, "      </PointData>\n");

    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </StructuredGrid>\n");
    fprintf(fp, "  <AppendedData encoding=\"raw\">\n");
    fprintf(fp, "_");

    /* Return new offset */
    *boffsetOut = boffset;
    return 0;
}

/* -- WriteVTSXMLFooter ---------------------------------------------
   Closes the XML for .vts:
     </AppendedData>
     </VTKFile>
-------------------------------------------------------------------- */
static int WriteVTSXMLFooter(FILE *fp)
{
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLFooter - Closing .vts XML.\n");

    fprintf(fp, "\n </AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");

    return 0;
}

/* -- WriteVTPXMLHeader ---------------------------------------------
   Writes the XML tags for a .vtp (PolyData) file:
     <VTKFile type="PolyData" ...>
       <PolyData>
         <Piece NumberOfPoints="npoints" NumberOfVerts="npoints">
           <Points>
             <DataArray ... offset="..." />
           </Points>
           <PointData>
             <DataArray ... offset="..." />
           </PointData>
           <Verts>
             <DataArray connectivity offset="..." />
             <DataArray offsets offset="..." />
           </Verts>
         </Piece>
       </PolyData>
       <AppendedData encoding="raw">
       _
   We update boffset for coords, scalar array, connectivity, offsets.
-------------------------------------------------------------------- */
static int WriteVTPXMLHeader(FILE       *fp,
                             int         npoints,
                             const char *fieldName,
                             int         boffset,
                             int        *boffsetOut)
{
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Writing .vtp header for %d points.\n", npoints);

    const char *byte_order = "LittleEndian";
    const char *precision  = "Float64";

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"%s\">\n",
            byte_order);
    fprintf(fp, "  <PolyData>\n");
    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"%d\">\n",
            npoints, npoints);

    /* Points => 3*npoints => offset1 */
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" "
                "format=\"appended\" offset=\"%d\" />\n",
            precision, boffset);
    boffset += (int)sizeof(int) + 3 * npoints * (int)sizeof(double);
    fprintf(fp, "      </Points>\n");

    /* Single scalar field => offset2 */
    fprintf(fp, "      <PointData>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"1\" "
                "format=\"appended\" offset=\"%d\" />\n",
            precision, fieldName, boffset);
    boffset += (int)sizeof(int) + npoints * (int)sizeof(double);
    fprintf(fp, "      </PointData>\n");

    /* Verts => connectivity + offsets => offset3, offset4 */
    fprintf(fp, "      <Verts>\n");
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" "
                "offset=\"%d\" />\n", boffset);
    boffset += (int)sizeof(int) + npoints * (int)sizeof(int);

    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" "
                "offset=\"%d\" />\n", boffset);
    boffset += (int)sizeof(int) + npoints * (int)sizeof(int);
    fprintf(fp, "      </Verts>\n");

    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </PolyData>\n");
    fprintf(fp, "  <AppendedData encoding=\"raw\">\n");
    fprintf(fp, "_");

    *boffsetOut = boffset;
    return 0;
}

/* -- WriteVTPXMLFooter ---------------------------------------------
   Closes the XML for .vtp:
     </AppendedData>
     </VTKFile>
-------------------------------------------------------------------- */
static int WriteVTPXMLFooter(FILE *fp)
{
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLFooter - Closing .vtp XML.\n");

    fprintf(fp, "\n </AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");

    return 0;
}

/* -- WriteVTKFileHeader --------------------------------------------
   A wrapper that picks either WriteVTSXMLHeader or WriteVTPXMLHeader
   based on meta->fileType. 
   We maintain a 'boffset' for appended data offset tracking.
-------------------------------------------------------------------- */
static int WriteVTKFileHeader(FILE *fp,
                              const VTKMetaData *meta,
                              int boffset,
                              int *boffsetOut)
{
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTKFileHeader - Deciding .vts or .vtp from fileType=%d\n",
              (int)meta->fileType);

    if (meta->fileType == VTK_STRUCTURED) {
        return WriteVTSXMLHeader(fp,
                                 meta->mx, meta->my, meta->mz,
                                 meta->nnodes,
                                 meta->scalarFieldName,
                                 boffset, boffsetOut);
    } else {
        return WriteVTPXMLHeader(fp,
                                 meta->npoints,
                                 meta->scalarFieldName,
                                 boffset, boffsetOut);
    }
}

/* -- WriteVTKFileFooter --------------------------------------------
   A wrapper that picks either .vts or .vtp footer.
-------------------------------------------------------------------- */
static int WriteVTKFileFooter(FILE *fp,
                              const VTKMetaData *meta)
{
    if (meta->fileType == VTK_STRUCTURED) {
        return WriteVTSXMLFooter(fp);
    } else {
        return WriteVTPXMLFooter(fp);
    }
}

/* =====================================================================
   3) CreateVTKFileFromMetadata

   The "public" function that writes either a .vts or .vtp file
   depending on meta->fileType. It uses the private functions above.

   Sequence:
     1) On rank 0, open the file in "wb".
     2) Write header (structured or polydata).
     3) Write appended data blocks:
        - coords
        - scalar field
        - (if polydata) connectivity & offsets
     4) Write the footer
     5) Close file
   All other ranks do nothing.
===================================================================== */
int CreateVTKFileFromMetadata(const char       *filename,
                              const VTKMetaData *meta,
                              MPI_Comm          comm)
{
    int rank, size;
    int boffset = 0;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    LOG_ALLOW(GLOBAL, LOG_INFO,
              "CreateVTKFileFromMetadata - Writing file '%s' using fileType=%s\n",
              filename,
              (meta->fileType == VTK_STRUCTURED ? "VTK_STRUCTURED" : "VTK_POLYDATA"));

    /* Only rank 0 does the actual file writing. */
    if (!rank) {
        FILE *fp = fopen(filename, "wb");
        if (!fp) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "CreateVTKFileFromMetadata - Could not open '%s' for writing.\n",
                      filename);
            return 11;
        }

        /* 1) Write the file header => sets up appended data offsets. */
        WriteVTKFileHeader(fp, meta, boffset, &boffset);

        /* 2) Write the coordinate array. 
              - If .vts => length=3*nnodes
              - If .vtp => length=3*npoints
        */
        int ncoords = (meta->fileType == VTK_STRUCTURED)
                      ? 3 * meta->nnodes
                      : 3 * meta->npoints;
        if (ncoords > 0 && meta->coords) {
            WriteVTKAppendedBlock(fp, meta->coords, ncoords, sizeof(double));
        }

        /* 3) If we have a scalarField, write it. 
              - If .vts => length=nnodes
              - If .vtp => length=npoints
        */
        if (meta->scalarField && meta->numScalarFields > 0) {
            int nvals = (meta->fileType == VTK_STRUCTURED)
                        ? meta->nnodes
                        : meta->npoints;
            WriteVTKAppendedBlock(fp, meta->scalarField, nvals, sizeof(double));
        }

        /* 4) If .vtp => also write connectivity and offsets. */
        if (meta->fileType == VTK_POLYDATA) {
            WriteVTKAppendedBlock(fp, meta->connectivity, meta->npoints, sizeof(int));
            WriteVTKAppendedBlock(fp, meta->offsets, meta->npoints, sizeof(int));
        }

        /* 5) Write the footer and close. */
        WriteVTKFileFooter(fp, meta);
        fclose(fp);

        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "CreateVTKFileFromMetadata - Successfully wrote file '%s'.\n",
                  filename);
    }

    return 0; /* success */
}


