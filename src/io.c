/**
 * @file io.c
 * @brief Handles input/output operations for velocity, pressure, grid, and other simulation fields.
 *
 * This file provides functions for reading and writing simulation fields to binary files,
 * including optional handling for statistical, LES, and RANS data.
 */

//#include <petsc.h>          // System dependency.
//#include "common.h"         // For UserCtx and shared definitions.
#include "io.h"             

// ------------------------ Function Definitions ------------------------

/**
 * @brief Reads grid generation parameters from input options or configuration.
 *
 * This function reads domain dimensions and block parameters for programmatically generating a grid.
 *
 * @param[in,out] user      Pointer to the UserCtx structure containing grid details.
 * @param[out] grid1d       Pointer to flag indicating if the grid is 1D (1) or 3D (0).
 * @param[out]    xMax          Pointer to store the upper bound in x-direction.
 * @param[out]    yMin          Pointer to store the lower bound in y-direction.
 * @param[out]    yMax          Pointer to store the upper bound in y-direction.
 * @param[out]    zMin          Pointer to store the lower bound in z-direction.
 * @param[out]    zMax          Pointer to store the upper bound in z-direction.
 * @param[out] imm          Pointer to array storing the i-dimensions of each block.
 * @param[out] jmm          Pointer to array storing the j-dimensions of each block.
 * @param[out] kmm          Pointer to array storing the k-dimensions of each block.
 * @param[out] nblk         Pointer to number of blocks.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadGridGenerationInputs(UserCtx *user, PetscInt *grid1d, PetscReal *xMin, PetscReal *xMax,
                               PetscReal *yMin, PetscReal *yMax, PetscReal *zMin, PetscReal *zMax,  PetscInt **imm, PetscInt **jmm, PetscInt **kmm,
                               PetscInt *nblk ){
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
    ierr = PetscOptionsGetReal(NULL, NULL, "-xMin", xMin, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-yMin", yMin, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-zMin", zMin, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-xMax", xMax, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-yMax", yMax, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-zMax", zMax, NULL); CHKERRQ(ierr);
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
 * @brief Reads grid parameters from a file, ignoring an optional "FDFGRID" header line.
 *
 * The first non‑blank line may be the literal string "FDFGRID". If present, it is skipped.
 * The next line must contain the number of blocks (nblk).  Subsequent lines list the
 * I-, J-, and K‑dimensions for each block.
 *
 * All data are broadcast to every MPI rank.
 *
 * @param[in]  filename  Name of the file containing grid data.
 * @param[out] nblk      Pointer to store the number of blocks.
 * @param[out] imm       Pointer to array of I‑dimensions (length nblk).
 * @param[out] jmm       Pointer to array of J‑dimensions (length nblk).
 * @param[out] kmm       Pointer to array of K‑dimensions (length nblk).
 * @param[out] grid1d    Pointer to store the grid type (1 = 1‑D, 0 = 3‑D).
 * @param[in]  comm      MPI communicator for broadcasting data.
 *
 * @return PetscErrorCode 0 on success, otherwise a PETSc error code.
 */
PetscErrorCode ReadGridFile(const char *filename, PetscInt *nblk,
                            PetscInt **imm, PetscInt **jmm, PetscInt **kmm,
                            PetscInt *grid1d, MPI_Comm comm)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadGridFile ‑ reading grid data from file: %s\n", filename);

    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);

    /* ---------------- Rank‑0 reads the file ---------------- */
    if (rank == 0) {
        FILE *fd = fopen(filename, "r");
        if (!fd) SETERRQ(comm, PETSC_ERR_FILE_OPEN, "Cannot open file: %s", filename);

        /* ---- Read first token; decide whether it is the header or nblk ---- */
        char firstTok[32] = {0};
        if (fscanf(fd, "%31s", firstTok) != 1)
            SETERRQ(comm, PETSC_ERR_FILE_READ, "Empty grid file: %s", filename);

        if (strcmp(firstTok, "FDFGRID") == 0) {
            /* Header is present – read nblk from the next line */
            if (fscanf(fd, "%" PetscInt_FMT, nblk) != 1)
                SETERRQ(comm, PETSC_ERR_FILE_READ,
                        "Expected number of blocks after \"FDFGRID\" in %s", filename);
        } else {
            /* No header – the token we just read is actually nblk */
            *nblk = (PetscInt)strtol(firstTok, NULL, 10);
        }

        LOG_ALLOW(GLOBAL, LOG_INFO, "Number of grid blocks: %d\n", *nblk);
        if (*nblk == 0) *nblk = 1;               /* default to one block */

        /* Allocate arrays */
        ierr = PetscCalloc1(*nblk, imm); CHKERRQ(ierr);
        ierr = PetscCalloc1(*nblk, jmm); CHKERRQ(ierr);
        ierr = PetscCalloc1(*nblk, kmm); CHKERRQ(ierr);

        /* Read block sizes */
        for (PetscInt i = 0; i < *nblk; ++i) {
            if (fscanf(fd, "%" PetscInt_FMT " %" PetscInt_FMT " %" PetscInt_FMT,
                       &(*imm)[i], &(*jmm)[i], &(*kmm)[i]) != 3)
                SETERRQ(comm, PETSC_ERR_FILE_READ,
                        "Expected 3 integers for block %d in %s", (int)i, filename);
        }

        /* 1‑D grid if there is only one block and J=K=1 */
        *grid1d = ((*nblk == 1) && ((*jmm)[0] == 1) && ((*kmm)[0] == 1)) ? 1 : 0;

        fclose(fd);
    }

    /* ------------- Broadcast metadata and arrays to every rank ------------- */
    ierr = MPI_Bcast(nblk, 1, MPI_INT, 0, comm);                CHKERRQ(ierr);

    if (rank != 0) { /* allocate on the other ranks before receiving */
        ierr = PetscCalloc1(*nblk, imm); CHKERRQ(ierr);
        ierr = PetscCalloc1(*nblk, jmm); CHKERRQ(ierr);
        ierr = PetscCalloc1(*nblk, kmm); CHKERRQ(ierr);
    }
    ierr = MPI_Bcast(*imm,  *nblk, MPI_INT, 0, comm);           CHKERRQ(ierr);
    ierr = MPI_Bcast(*jmm,  *nblk, MPI_INT, 0, comm);           CHKERRQ(ierr);
    ierr = MPI_Bcast(*kmm,  *nblk, MPI_INT, 0, comm);           CHKERRQ(ierr);
    ierr = MPI_Bcast(grid1d, 1,     MPI_INT, 0, comm);          CHKERRQ(ierr);

    /* ------------- Verbose diagnostic output ------------- */
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Grid is %s\n", *grid1d ? "1‑D" : "3‑D");
    for (PetscInt i = 0; i < *nblk; ++i)
        LOG_ALLOW(LOCAL, LOG_DEBUG,
                  "Block %d: IM=%d  JM=%d  KM=%d\n",
                  (int)i, (int)(*imm)[i], (int)(*jmm)[i], (int)(*kmm)[i]);

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadGridFile – grid file data retrieved successfully\n");
    return 0;
}

// Helper function to construct the filename and test existence
// This avoids repeating the snprintf and PetscTestFile logic for each field.
static PetscErrorCode CheckDataFile(PetscInt ti, const char *fieldName, const char *ext, PetscBool *fileExists)
{
    PetscErrorCode ierr;
    char           filename[PETSC_MAX_PATH_LEN];
    PetscMPIInt    rank; // Use PetscMPIInt for MPI rank

    PetscFunctionBeginUser;
    *fileExists = PETSC_FALSE; // Default to not existing

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Construct filename: results/FIELDNAME<#####>_RANK.EXT
    // Ensure the format string correctly uses PetscInt_FMT for ti if it's PetscInt
    // and %d for rank if it's int/PetscMPIInt.
    // Based on your logs: "results/ufield00000_0.dat"
    // fieldName, ti (5 digits), rank (integer), ext
    ierr = PetscSNPrintf(filename, sizeof(filename), "results/%s%05" PetscInt_FMT "_%d.%s",
                         fieldName, ti, (int)rank, ext); CHKERRQ(ierr);
    
    ierr = PetscTestFile(filename, 'r', fileExists); CHKERRQ(ierr);

    if (!(*fileExists)) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "CheckDataFile - Optional data file '%s' for field '%s' at ti=%d not found.\n", filename, fieldName, ti);
    } else {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "CheckDataFile - Data file '%s' for field '%s' at ti=%d found.\n", filename, fieldName, ti);
    }
    PetscFunctionReturn(0);
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
    PetscErrorCode ierr, ierr_load; // Use a separate variable for VecLoad's potential error
    PetscViewer    viewer;
    char           filename[PETSC_MAX_PATH_LEN]; // Use PETSc constant
    PetscBool      fileExists;
    PetscMPIInt    rank;

    PetscFunctionBeginUser; // Use User version
    // Validate inputs
    if (!user || !field_name || !field_vec || !ext) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "ReadFieldData - Null argument provided.");
    }
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr); // Use CHKERRQ for MPI errors

    // Construct filename
    // Consider making results directory path configurable
    ierr = PetscSNPrintf(filename, sizeof(filename), "results/%s%05" PetscInt_FMT "_%d.%s", field_name, ti, 0 /*user->_this placeholder*/, ext); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadFieldData - Attempting to read file: %s\n", filename);

    // Check file existence
    ierr = PetscTestFile(filename, FILE_MODE_READ, &fileExists); CHKERRQ(ierr);
    if (!fileExists) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "ReadFieldData - File '%s' does not exist. Skipping field read.\n", filename);
        // Return a specific error code that the caller (ReadAllSwarmFields) can check
        PetscFunctionReturn(PETSC_ERR_FILE_OPEN);
    }

    // Open the viewer
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &viewer);
    // Handle viewer open error separately - might be permissions etc.
    if (ierr) {
         LOG_ALLOW(GLOBAL, LOG_WARNING, "ReadFieldData - PetscViewerBinaryOpen failed for '%s' (Error code %d). Skipping field read.\n", filename, ierr);
         // We still want to return PETSC_ERR_FILE_OPEN to the caller
         PetscFunctionReturn(PETSC_ERR_FILE_OPEN);
    }

    // Attempt to load data, store return code separately
    ierr_load = VecLoad(field_vec, viewer);

    // Always destroy the viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr); // Use CHKERRQ for critical errors like viewer destroy failure

    // Check the return code from VecLoad *after* destroying the viewer
    if (ierr_load) {
        // An error occurred during VecLoad (e.g., size mismatch, corrupted data)
        LOG_ALLOW(GLOBAL, LOG_WARNING,"ReadFieldData - VecLoad failed for file '%s' (Error code %d). Vector data may be incorrect/incomplete. Continuing post-processing.\n", filename, ierr_load);
        // Decide how to proceed. Returning 0 allows post-processing to continue,
        // but the vector might be in an inconsistent state (e.g., wrong size
        // if the resize failed internally before erroring).
        // Returning the error might be safer if subsequent steps rely on valid data.
        // Let's return 0 for now, assuming the caller checks vector size.
        PetscFunctionReturn(0);
         // Alternatively, return the specific error:
         // PetscFunctionReturn(ierr_load);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "ReadFieldData - Successfully loaded data for field: %s from %s\n", field_name, filename);

    PetscFunctionReturn(0);
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

    // Read Particle Count field
    ierr = ReadFieldData(user, "ParticleCount", user->ParticleCount, ti, "dat"); CHKERRQ(ierr);

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

    PetscReal vmin, vmax;
    ierr = VecMin(field_vec, NULL, &vmin); CHKERRQ(ierr);
    ierr = VecMax(field_vec, NULL, &vmax); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"% s step %d  min=%.6e  max=%.6e\n",field_name, ti, (double)vmin, (double)vmax);

    
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
    ierr = WriteFieldData(user, "vfield", user->Ucont, user->step, "dat"); CHKERRQ(ierr);

    // Write Cartesian velocity field
    ierr = WriteFieldData(user, "ufield", user->Ucat, user->step, "dat"); CHKERRQ(ierr);

    // Write pressure field
    ierr = WriteFieldData(user, "pfield", user->P, user->step, "dat"); CHKERRQ(ierr);

    // Write node state field (nvert)
    ierr = WriteFieldData(user, "nvfield", user->Nvert, user->step, "dat"); CHKERRQ(ierr);

    ierr = WriteFieldData(user, "ParticleCount",user->ParticleCount,user->step,"dat"); CHKERRQ(ierr);

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

    ierr = WriteFieldData(user, "su0", user->Ucat_sum, user->step, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "su1", user->Ucat_cross_sum, user->step, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "su2", user->Ucat_square_sum, user->step, "dat"); CHKERRQ(ierr);
    ierr = WriteFieldData(user, "sp", user->P_sum, user->step, "dat"); CHKERRQ(ierr);

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
    ierr = WriteFieldData(user, "cs", Cs, user->step, "dat"); CHKERRQ(ierr);
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

    ierr = WriteFieldData(user, "kfield", user->K_Omega, user->step, "dat"); CHKERRQ(ierr);

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

static PetscInt WriteVTKAppendedBlock(FILE *fp,
                                 const void *buf,
                                 PetscInt nvals,
                                 PetscInt dsize);

static PetscInt WriteVTSXMLHeader(FILE       *fp,
                             PetscInt         mx,
                             PetscInt         my,
                             PetscInt         mz,
                             PetscInt         nnodes,
                             const char *fieldName,
                             PetscInt         numComponents,
                             PetscInt         boffset,
                             PetscInt        *boffsetOut);

static PetscInt WriteVTSXMLFooter(FILE *fp);

static PetscInt WriteVTPXMLHeader(FILE       *fp,
                             PetscInt         npoints,
                             const char *fieldName,
                             PetscInt         boffset,
                             PetscInt        *boffsetOut);

static PetscInt WriteVTPXMLFooter(FILE *fp);

static PetscInt WriteVTKFileHeader(FILE *fp,
                              const VTKMetaData *meta,
                              PetscInt boffset,
                              PetscInt *boffsetOut);

static PetscInt WriteVTKFileFooter(FILE *fp,
                              const VTKMetaData *meta);


=====================================================================
   1) ReadDataFileToArray

   See the function-level comments in io.h for a summary.
   This reads one value per line from an ASCII file. Rank 0 does I/O,
   broadcasts the data to all ranks.
===================================================================== */
PetscInt ReadDataFileToArray(const char   *filename,
                        double      **data_out,
                        PetscInt          *Nout,
                        MPI_Comm      comm)
{
    /* STEP 0: Prepare local variables & log function entry */
    PetscMPIInt    rank, size;
    PetscErrorCode ierr;
    FILE  *fp = NULL;
    PetscInt    N   = 0;            /* number of lines/values read on rank 0 */
    double *array = NULL;      /* pointer to local array on each rank */
    PetscInt    fileExistsFlag = 0; /* 0 = doesn't exist, 1 = does exist */

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "ReadDataFileToArray - Start reading from file: %s\n",
              filename);

    /* Basic error checking: data_out, Nout must be non-null. */
    if (!filename || !data_out || !Nout) {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "ReadDataFileToArray - Null pointer argument provided.\n");
        return 1;
    }

    /* Determine rank/size for coordinating I/O. */
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* STEP 1: On rank 0, check if file can be opened. */
    if (!rank) {
        fp = fopen(filename, "r");
        if (fp) {
            fileExistsFlag = 1;
            fclose(fp);
        }
    }

    /* STEP 2: Broadcast file existence to all ranks. */
    // In ReadDataFileToArray:
    ierr = MPI_Bcast(&fileExistsFlag, 1, MPI_INT, 0, comm); CHKERRQ(ierr);

    if (!fileExistsFlag) {
        /* If file does not exist, log & return. */
        if (!rank) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ReadDataFileToArray - File '%s' not found.\n",
                      filename);
        }
        return 2;
    }

    /* STEP 3: Rank 0 re-opens and reads the file, counting lines, etc. */
    if (!rank) {
        fp = fopen(filename, "r");
        if (!fp) {
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

        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "ReadDataFileToArray - File '%s' has %d lines.\n",
                  filename, N);

        /* (3b) Allocate array on rank 0. */
        array = (double*)malloc(N * sizeof(double));
        if (!array) {
            fclose(fp);
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ReadDataFileToArray - malloc failed for array.\n");
            return 4;
        }

        /* (3c) Rewind & read values into array. */
        rewind(fp);
        {
            PetscInt i = 0;
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

    /* STEP 4: Broadcast the integer N to all ranks. */
    ierr = MPI_Bcast(&N, 1, MPI_INT, 0, comm); CHKERRQ(ierr);

    /* STEP 5: Each rank allocates an array to receive the broadcast if rank>0. */
    if (rank) {
        array = (double*)malloc(N * sizeof(double));
        if (!array) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ReadDataFileToArray - malloc failed on rank %d.\n",
                      rank);
            return 5;
        }
    }

    /* STEP 6: Broadcast the actual data from rank 0 to all. */
    ierr = MPI_Bcast(array, N, MPI_DOUBLE, 0, comm); CHKERRQ(ierr);

    /* STEP 7: Assign outputs on all ranks. */
    *data_out = array;
    *Nout     = N;

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "ReadDataFileToArray - Done. Provided array of length=%d to all ranks.\n",
              N);
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
 * @return PetscInt  Returns 0 if successful, non-zero otherwise.
 */
PetscInt WriteVTKAppendedBlock(FILE *fp, const void *data, PetscInt num_elements, size_t element_size) {

    // Log the function call with parameters
  LOG_ALLOW_SYNC(LOCAL,LOG_INFO,"WriteVTKAppendedBlock - Called with %d elements, each of size %zu bytes.\n",
                   num_elements, element_size);

    // Calculate the block size
    PetscInt block_size = num_elements * (PetscInt)element_size;
    LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "WriteVTKAppendedBlock - Calculated block size: %d bytes.\n", block_size);

    // Write the block size as a 4-byte integer
    if (fwrite(&block_size, sizeof(PetscInt), 1, fp) != 1) {
        fprintf(stderr, "[ERROR] Failed to write block size.\n");
        LOG_ALLOW_SYNC(LOCAL, LOG_ERROR, "WriteVTKAppendedBlock - Error writing block size.\n");
        return 1;
    }

    // Write the actual data
    if (fwrite(data, element_size, num_elements, fp) != (size_t)num_elements) {
        fprintf(stderr, "[ERROR] Failed to write data block.\n");
        LOG_ALLOW_SYNC(LOCAL, LOG_ERROR, "WriteVTKAppendedBlock - Error writing data block.\n");
        return 1;
    }

    // Log success
    LOG_ALLOW_SYNC(LOCAL, LOG_INFO, "WriteVTKAppendedBlock - Successfully wrote block of %d bytes.\n", block_size);

    return 0; // Success
}

/**
 * @brief Writes the XML header portion of a .vts file, defining the structured grid extents and data arrays.
 *
 * This function prPetscInts XML tags for the grid, including whole and piece extents,
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
 * @param[in]  numComponents Number of components in the fieldName array (e.g., 3 for velocity).
 * @param[in]  boffset    Current byte offset in the appended data section.
 * @param[out] boffsetOut Updated byte offset after writing the header information.
 *
 * @return 0 on success, non-zero on failure (though no explicit failure cases are handled here).
 */
static PetscInt WriteVTSXMLHeader(FILE       *fp,
                             PetscInt         mx,
                             PetscInt         my,
                             PetscInt         mz,
                             PetscInt         nnodes,
                             const char *fieldName,
                             PetscInt         numComponents,
                             PetscInt         boffset,
                             PetscInt        *boffsetOut)
{
    // Log the entry into this function and the parameters
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLHeader - Writing .vts header (fieldName=%s, numComponents=%d,boffset=%d,mx=%d, my=%d, mz=%d, nnodes=%d).\n",fieldName,numComponents,boffset,
              mx, my, mz, nnodes);

    // Set XML configuration strings
    const char *byte_order = "LittleEndian";
    const char *precision  = "Float64";

    // Begin XML header
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",
            byte_order);
    fprintf(fp, "  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n",
            0, mx-1, 0, my-1, 0, mz-1);
    fprintf(fp, "    <Piece Extent=\"%d %d %d %d %d %d\">\n",
            0, mx-1, 0, my-1, 0, mz-1);

    // Points section
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" "
                "format=\"appended\" offset=\"%d\" />\n",
            precision, boffset);
    boffset += (PetscInt)sizeof(PetscInt) + 3 * nnodes * (PetscInt)sizeof(double);
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLHeader - Updated offset to %d after Points.\n", boffset);
    fprintf(fp, "      </Points>\n");

    // PointData => "ufield,vfield,p
    fprintf(fp, "      <PointData Scalars=\"%s\">\n", fieldName);
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%d\" "
                "format=\"appended\" offset=\"%d\" />\n",
            precision, fieldName, numComponents,boffset);
    boffset += (PetscInt)sizeof(PetscInt) + nnodes * numComponents * (PetscInt)sizeof(double);
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
 * @return PetscInt  Returns 0 on success.
 */
static PetscInt WriteVTSXMLFooter(FILE *fp)
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
 * @return PetscInt  Returns 0 on success.
 */
/**
 * @brief Writes the XML header portion of a .vtp file for a point-cloud representation.
 *
 * This function sets up the appended data sections for point coordinates, a named field
 * (e.g., velocity) with a specified number of components, and the connectivity/offsets data for
 * each point as a vertex. The byte offset for the appended data is updated after each
 * section. The data type for connectivity and offsets ("Int32" or "Int64") is now
 * determined dynamically based on the size of PetscInt during compilation.
 *
 * @param[in]  fp            File pointer for writing the XML header (already open).
 * @param[in]  npoints       Number of points in the point-cloud.
 * @param[in]  fieldName     Name of the point data field (e.g., "velocity").
 * @param[in]  numComponents Number of components in the fieldName array (e.g., 3 for velocity).
 * @param[in]  boffset       Current byte offset into the appended data section.
 * @param[out] boffsetOut    Updated byte offset after writing the header information.
 *
 * @return PetscInt  Returns 0 on success.
 */
static PetscInt WriteVTPXMLHeader(FILE       *fp,
                                  PetscInt    npoints,
                                  const char *fieldName,
                                  PetscInt    numComponents,
                                  PetscInt    boffset,
                                  PetscInt   *boffsetOut)
{
    // Log function entry with parameters
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Called with npoints=%d, fieldName=%s, numComponents=%d, boffset=%d.\n",
              npoints, fieldName, numComponents, boffset);

    // Determine data types based on compilation settings
    const char *precision     = "Float64"; // Assuming coordinates and field data are doubles
    const char *int_type_str  = (sizeof(PetscInt) == 8) ? "Int64" : "Int32"; // Dynamic check for PetscInt size

    // Pre-calculate byte sizes for data blocks
    PetscInt points_bytes       = npoints * 3 * sizeof(double); // 3 components, double precision
    PetscInt field_data_bytes   = npoints * numComponents * sizeof(double); // numComponents, double precision
    PetscInt connectivity_bytes = npoints * sizeof(PetscInt); // Based on actual PetscInt size
    PetscInt offsets_bytes      = npoints * sizeof(PetscInt); // Based on actual PetscInt size

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Determined types: precision=%s, int_type=%s (sizeof(PetscInt)=%zu)\n",
              precision, int_type_str, sizeof(PetscInt));

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
    boffset += sizeof(PetscInt) + points_bytes; // Increment by size prefix + actual data bytes
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Updated offset to %d after Points (%d bytes data).\n", boffset, points_bytes);
    fprintf(fp, "      </Points>\n");

    // 3) PointData => e.g. "velocity" => next offset
    fprintf(fp, "      <PointData>\n");
    fprintf(fp,
            "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%d\" "
            "format=\"appended\" offset=\"%d\"/>\n",
            precision, fieldName, numComponents, boffset);
    boffset += sizeof(PetscInt) + field_data_bytes; // Increment by size prefix + actual data bytes
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Updated offset to %d after PointData ('%s', %d bytes data).\n", boffset, fieldName, field_data_bytes);
    fprintf(fp, "      </PointData>\n");

    // 4) Verts => connectivity + offsets
    fprintf(fp, "      <Verts>\n");
    // Use the dynamically determined integer type string
    fprintf(fp,
            "        <DataArray type=\"%s\" Name=\"connectivity\" format=\"appended\" offset=\"%d\"/>\n",
            int_type_str, boffset);
    boffset += sizeof(PetscInt) + connectivity_bytes; // Increment by size prefix + actual data bytes
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Updated offset to %d after connectivity (%s, %d bytes data).\n", boffset, int_type_str, connectivity_bytes);

    // Use the dynamically determined integer type string
    fprintf(fp,
            "        <DataArray type=\"%s\" Name=\"offsets\" format=\"appended\" offset=\"%d\"/>\n",
            int_type_str, boffset);
    boffset += sizeof(PetscInt) + offsets_bytes; // Increment by size prefix + actual data bytes
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Updated offset to %d after offsets (%s, %d bytes data).\n", boffset, int_type_str, offsets_bytes);
    fprintf(fp, "      </Verts>\n");

    // 5) Close out the piece + PolyData
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </PolyData>\n");

    // 6) <AppendedData encoding="raw"> + underscore
    //    => all binary must follow AFTER this line
    fprintf(fp, "  <AppendedData encoding=\"raw\">\n");
    fprintf(fp, "_"); // This underscore MUST be present before the binary data blocks

    // Assign the final calculated offset to the output parameter
    *boffsetOut = boffset;

    // Log function completion
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTPXMLHeader - Completed writing VTP XML header. Final calculated boffset (before actual data write)=%d.\n",
              boffset);

    return 0; // Success
}


/**
 * @brief Finalizes the .vtp file by closing the appended data section and the VTKFile tag.
 *
 * This function writes the closing XML tags for the appended data section 
 * and completes the overall VTK structure, ensuring a well-formed .vtp file.
 *
 * @param[in] fp  File pointer to the .vtp file (open for writing).
 *
 * @return PetscInt  Returns 0 on success.
 */
static PetscInt WriteVTPXMLFooter(FILE *fp)
{
    // Log the function entry
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTPXMLFooter - Closing VTP XML.\n");

    // no stray prPetscInts or data -> directly close
    fprintf(fp, "\n  </AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");

    // Log successful completion
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTPXMLFooter - Completed writing VTP XML footer.\n");

    return 0;
}

/**
 * @brief Writes the initial XML header for a VTK file based on the provided metadata.
 *
 * This function writes the XML header for VTK output. It now supports both polydata
 * (VTK_POLYDATA) and structured grid (VTK_STRUCTURED) file types. For each type, the
 * function determines whether to write a scalar or vector field header by inspecting
 * the metadata. For polydata, it calls \c WriteVTPXMLHeader and for structured grids,
 * it calls \c WriteVTSXMLHeader.
 *
 * @param[in]  fp          File pointer (open for writing) to which the header will be written.
 * @param[in]  meta        Pointer to a \c VTKMetaData structure containing fileType, number of points,
 *                         grid dimensions, field names, and other necessary information.
 * @param[in]  boffset     Current byte offset in the appended data section.
 * @param[out] boffsetOut  Updated byte offset after writing any header-related data.
 *
 * @return PetscInt Returns 0 on success, or -1 if the fileType is not handled.
 *
 * @note The function uses PETSc-style logging (LOG_ALLOW) to trace its execution.
 *       For structured grid output, it assumes that grid dimensions (mx, my, mz) and
 *       the number of nodes (nnodes) have been set in the metadata.
 */
static PetscInt WriteVTKFileHeader(FILE *fp, const VTKMetaData *meta, PetscInt boffset, PetscInt *boffsetOut)
{
    const char *fieldName = NULL;
    PetscInt numComponents = 1;

    /*-----------------------------------------------------------------------
     * Log the entry into this function with the key metadata parameters.
     *-----------------------------------------------------------------------*/
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTKFileHeader - Entered: fileType=%d, npoints=%d, scalarField=%s, vectorField=%s, "
              "numVectorFields=%d, mx=%d, my=%d, mz=%d, nnodes=%d.\n",
              meta->fileType, meta->npoints, meta->scalarFieldName, meta->vectorFieldName,
              meta->numVectorFields, meta->mx, meta->my, meta->mz, meta->nnodes);

    if (meta->fileType == VTK_POLYDATA) {
        /*===================================================================
         * Polydata Branch: VTK_POLYDATA
         *===================================================================*/
        LOG_ALLOW(GLOBAL, LOG_INFO, "WriteVTKFileHeader - Processing VTK_POLYDATA output.\n");

        /* Determine the field type for polydata:
         * If a vector field is available, use it (with 3 components);
         * otherwise, use the scalar field (with 1 component).
         */
        if (meta->numVectorFields > 0 && meta->vectorFieldName) {
            fieldName = meta->vectorFieldName;
            numComponents = 3;
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "WriteVTKFileHeader - Using vector field '%s' with 3 components for polydata.\n",
                      fieldName);
        } else {
            fieldName = meta->scalarFieldName;
            numComponents = 1;
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "WriteVTKFileHeader - Using scalar field '%s' with 1 component for polydata.\n",
                      fieldName);
        }

        /* Log the call to the polydata header-writing function */
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "WriteVTKFileHeader - Calling WriteVTPXMLHeader with boffset=%d.\n", boffset);
        return WriteVTPXMLHeader(fp,
                                 meta->npoints,
                                 fieldName,
                                 numComponents,
                                 boffset,
                                 boffsetOut);
    } else if (meta->fileType == VTK_STRUCTURED) {
        /*===================================================================
         * Structured Grid Branch: VTK_STRUCTURED
         *===================================================================*/
        LOG_ALLOW(GLOBAL, LOG_INFO, "WriteVTKFileHeader - Processing VTK_STRUCTURED output.\n");

        /* Determine the field type for structured grid:
         * If a vector field is available, use it (with 3 components);
         * otherwise, use the scalar field (with 1 component).
         */
        if (meta->numVectorFields > 0 && meta->vectorFieldName) {
            fieldName = meta->vectorFieldName;
            numComponents = 3;
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "WriteVTKFileHeader - Using vector field '%s' with 3 components for structured grid.\n",
                      fieldName);
        } else {
            fieldName = meta->scalarFieldName;
            numComponents = 1;
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "WriteVTKFileHeader - Using scalar field '%s' with 1 component for structured grid.\n",
                      fieldName);
        }

        /* Log the grid dimensions and number of nodes, then call the structured grid header-writing function */
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "WriteVTKFileHeader - Calling WriteVTSXMLHeader with boffset=%d, grid dimensions (%d, %d, %d) and nnodes=%d.\n",
                  boffset, meta->mx, meta->my, meta->mz, meta->nnodes);
        return WriteVTSXMLHeader(fp,
                                 meta->mx,
                                 meta->my,
                                 meta->mz,
                                 meta->nnodes,
                                 fieldName,
                                 numComponents,
                                 boffset,
                                 boffsetOut);
    } else {
        /*-------------------------------------------------------------------
         * Unsupported File Type: Log a warning and return an error.
         *-------------------------------------------------------------------*/
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "WriteVTKFileHeader - Unsupported file type %d encountered.\n", meta->fileType);
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
 * @return PetscInt  Returns 0 or the value from \c WriteVTPXMLFooter on success, 
 *              or -1 if the file type is not handled.
 */
static PetscInt WriteVTKFileFooter(FILE *fp, const VTKMetaData *meta)
{
    // Log the entry into this function
    LOG_ALLOW(GLOBAL, LOG_DEBUG, 
              "WriteVTKFileFooter - Called with fileType=%d.\n", meta->fileType);

    if (meta->fileType == VTK_POLYDATA) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTKFileFooter - Writing VTP XML footer.\n");
        return WriteVTPXMLFooter(fp);
    } else {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTKFileFooter - Writing VTS XML footer.\n");
	return WriteVTSXMLFooter(fp);
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
 * @return PetscInt  Returns 0 on success, or 1 on file-opening failure. Always returns 0 on non-\c rank 0 processes.
 */
PetscInt CreateVTKFileFromMetadata(const char *filename,
                                   const VTKMetaData *meta,
                                   MPI_Comm comm)
{
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Entry point. Filename: %s\n", filename);

    PetscMPIInt rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - rank=%d, size=%d.\n", rank, size);

    // 1) Only rank 0 writes the file
    if (!rank) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "CreateVTKFileFromMetadata - Rank 0 writing file '%s'.\n", filename);

        FILE *fp = fopen(filename, "wb");
        if (!fp) {
            LOG_ALLOW(GLOBAL, LOG_ERROR, "CreateVTKFileFromMetadata - fopen failed for %s.\n", filename);
            return 1;
        }
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Successfully opened file: %s\n", filename);

        PetscInt boffset = 0;
        // 2) Write the XML header => sets <DataArray ... offset="...">
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing header (initial boffset=%d).\n", boffset);
        WriteVTKFileHeader(fp, meta, boffset, &boffset);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Header written (updated boffset=%d).\n", boffset);

        // 3) Immediately write appended data blocks in EXACT order:
        //    (a) coords (3*npoints doubles)
        if (meta->coords) {
	  PetscInt ncoords = 3 * (meta->fileType == VTK_STRUCTURED ? meta->nnodes : meta->npoints);
            if (ncoords > 0) {
                LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing coords block: %d doubles.\n", ncoords);
                WriteVTKAppendedBlock(fp, meta->coords, ncoords, sizeof(double));
                boffset += sizeof(PetscInt) + ncoords * sizeof(double);
            }
        }

        //    (b) scalar or vector field
        if (meta->numScalarFields > 0 && meta->scalarField) {
	  PetscInt nvals = (meta->fileType == VTK_STRUCTURED ? meta->nnodes : meta->npoints);
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing scalar field block: %d doubles.\n", nvals);
            WriteVTKAppendedBlock(fp, meta->scalarField, nvals, sizeof(double));
            boffset += sizeof(PetscInt) + nvals * sizeof(double);
        }
        else if (meta->numVectorFields > 0 && meta->vectorField) {
	  PetscInt nvals = 3 * (meta->fileType == VTK_STRUCTURED ? meta->nnodes : meta->npoints); // 3 components
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing vector field block: %d doubles.\n", nvals);
            WriteVTKAppendedBlock(fp, meta->vectorField, nvals, sizeof(double));
            boffset += sizeof(PetscInt) + nvals * sizeof(double);
        }

        //    (c) connectivity (if polydata)
        if (meta->fileType == VTK_POLYDATA && meta->connectivity) {
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing connectivity block: %d PetscInts.\n", meta->npoints);
            WriteVTKAppendedBlock(fp, meta->connectivity, meta->npoints, sizeof(PetscInt));
            boffset += sizeof(PetscInt) + meta->npoints * sizeof(PetscInt);
        }

        //    (d) offsets (if polydata)
        if (meta->fileType == VTK_POLYDATA && meta->offsets) {
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing offsets block: %d PetscInts.\n", meta->npoints);
            WriteVTKAppendedBlock(fp, meta->offsets, meta->npoints, sizeof(PetscInt));
            boffset += sizeof(PetscInt) + meta->npoints * sizeof(PetscInt);
        }

        // 4) Write footer => </AppendedData> & </VTKFile>
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateVTKFileFromMetadata - Writing footer.\n");
        WriteVTKFileFooter(fp, meta);

        fclose(fp);
        LOG_ALLOW(GLOBAL, LOG_INFO, "CreateVTKFileFromMetadata - File '%s' closed.\n", filename);
    }

    return 0;
}


/**
 * @brief Gathers a PETSc vector onto rank 0 as a contiguous array of doubles.
 *
 * This function retrieves the local portions of the input vector \p inVec from
 * all MPI ranks via \c MPI_Gatherv and assembles them into a single array on rank 0.
 * The global size of the vector is stored in \p N, and a pointer to the newly
 * allocated array is returned in \p arrayOut (valid only on rank 0).
 *
 * @param[in]  inVec      The PETSc vector to gather.
 * @param[out] N          The global size of the vector (output).
 * @param[out] arrayOut   On rank 0, points to the newly allocated array of size \p N.
 *                        On other ranks, it is set to NULL.
 *
 * @return PetscErrorCode  Returns 0 on success, or a non-zero PETSc error code.
 */
PetscErrorCode VecToArrayOnRank0(Vec inVec, PetscInt *N, double **arrayOut)
{
    PetscErrorCode    ierr;
    MPI_Comm          comm;
    PetscMPIInt       rank, size;
    PetscInt          globalSize, localSize;
    const PetscScalar *localArr = NULL;

    // Log entry into the function
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "VecToArrayOnRank0 - Start gathering vector onto rank 0.\n");

    /* Get MPI comm, rank, size */
    ierr = PetscObjectGetComm((PetscObject)inVec, &comm);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm, &size);CHKERRQ(ierr);

    /* Get global size (for the entire Vec) */
    ierr = VecGetSize(inVec, &globalSize);CHKERRQ(ierr);
    *N = globalSize;

    /* Get local size (portion on this rank) */
    ierr = VecGetLocalSize(inVec, &localSize);CHKERRQ(ierr);

    // Log vector sizes and process info
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "VecToArrayOnRank0 - rank=%d of %d, globalSize=%d, localSize=%d.\n",
              rank, size, globalSize, localSize);

    /* Access the local array data */
    ierr = VecGetArrayRead(inVec, &localArr);CHKERRQ(ierr);

    /*
       We'll gather the local chunks via MPI_Gatherv:
       - First, gather all local sizes into recvcounts[] on rank 0.
       - Then set up a displacement array (displs[]) to place each chunk in the correct spot.
       - Finally, gather the actual data.
    */
    PetscMPIInt *recvcounts = NULL;
    PetscMPIInt *displs     = NULL;

    if (!rank) {
        recvcounts = (PetscMPIInt *) malloc(size * sizeof(PetscMPIInt));
        displs     = (PetscMPIInt *) malloc(size * sizeof(PetscMPIInt));
    }

    /* Convert localSize (PetscInt) to PetscMPIInt for MPI calls */
    PetscMPIInt localSizeMPI = (PetscMPIInt)localSize;

    /* Gather local sizes to rank 0 */
    ierr = MPI_Gather(&localSizeMPI, 1, MPI_INT,
                      recvcounts,      1, MPI_INT,
                      0, comm);CHKERRQ(ierr);

    /* On rank 0, build displacements and allocate the big array */
    if (!rank) {
        displs[0] = 0;
        for (PetscMPIInt i = 1; i < size; i++) {
            displs[i] = displs[i - 1] + recvcounts[i - 1];
        }

        /* Allocate a buffer for the entire (global) array */
        *arrayOut = (double *) malloc(globalSize * sizeof(double));
        if (!(*arrayOut)) SETERRQ(comm, PETSC_ERR_MEM, "Failed to allocate array on rank 0.");
    } else {
        /* On other ranks, we do not allocate anything */
        *arrayOut = NULL;
    }

    /* Gather the actual data (assuming real scalars => MPI_DOUBLE) */
    ierr = MPI_Gatherv((void *) localArr,         /* sendbuf on this rank */
                       localSizeMPI, MPI_DOUBLE,  /* how many, and type */
                       (rank == 0 ? *arrayOut : NULL),  /* recvbuf on rank 0 */
                       (rank == 0 ? recvcounts : NULL),
                       (rank == 0 ? displs    : NULL),
                       MPI_DOUBLE, 0, comm);CHKERRQ(ierr);

    /* Restore local array (cleanup) */
    ierr = VecRestoreArrayRead(inVec, &localArr);CHKERRQ(ierr);

    if (!rank) {
        free(recvcounts);
        free(displs);
    }

    // Log successful completion
    LOG_ALLOW(GLOBAL, LOG_INFO, "VecToArrayOnRank0 - Successfully gathered data on rank 0.\n");

    PetscFunctionReturn(0);
}

/**
 * @brief Reads data from a file into a specified field of a PETSc DMSwarm.
 *
 * This function is the counterpart to WriteSwarmField(). It creates a global PETSc vector 
 * that references the specified DMSwarm field, uses ReadFieldData() to read the data from 
 * a file, and then destroys the global vector reference.
 *
 * @param[in]  user       Pointer to the UserCtx structure (containing `user->swarm`).
 * @param[in]  field_name Name of the DMSwarm field to read into (must be previously declared/allocated).
 * @param[in]  ti         Time index used to construct the input file name.
 * @param[in]  ext        File extension (e.g., "dat" or "bin").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 *
 * @note Compatible with PETSc 3.14.x.
 */
PetscErrorCode ReadSwarmField(UserCtx *user, const char *field_name, PetscInt ti, const char *ext)
{
  PetscErrorCode ierr;
  DM             swarm;
  Vec            fieldVec;

  PetscFunctionBegin;

  swarm = user->swarm;

  LOG_ALLOW(GLOBAL,LOG_DEBUG," ReadSwarmField Begins \n");
 
  /* 2) Create a global vector that references the specified Swarm field. */
  ierr = DMSwarmCreateGlobalVectorFromField(swarm, field_name, &fieldVec);CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL,LOG_DEBUG," Vector created from Field \n");

  /* 3) Use the provided ReadFieldData() function to read data into fieldVec. */
  ierr = ReadFieldData(user, field_name, fieldVec, ti, ext);CHKERRQ(ierr);

  /* 4) Destroy the global vector reference. */
  ierr = DMSwarmDestroyGlobalVectorFromField(swarm, field_name, &fieldVec);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/**
 * @brief Reads multiple fields (positions, velocity, CellID, and weight) into a DMSwarm.
 *
 * This function is analogous to ReadSimulationFields() but targets a DMSwarm. 
 * Each Swarm field is read from a separate file using ReadSwarmField().
 * 
 * @param[in,out] user Pointer to the UserCtx structure containing the DMSwarm (user->swarm).
 * @param[in]     ti   Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
/**
 * @brief Reads multiple fields (positions, velocity, etc.) into a DMSwarm.
 *
 * ASSUMES the swarm has already been resized (if necessary) by a function
 * like PreCheckAndResizeSwarm to match the particle count in the input files
 * for the current timestep.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the DMSwarm (user->swarm).
 * @param[in]     ti   Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure. Handles
 *         missing optional files as warnings.
 */
PetscErrorCode ReadAllSwarmFields(UserCtx *user, PetscInt ti)
{
  PetscErrorCode ierr;
  PetscInt nGlobal;
  PetscBool fileExists;

  PetscFunctionBegin;
  ierr = DMSwarmGetSize(user->swarm, &nGlobal); CHKERRQ(ierr); // Get final size
  LOG_ALLOW(GLOBAL, LOG_INFO, "ReadAllSwarmFields - Reading DMSwarm fields for timestep %d (expected size %d).\n", ti, nGlobal);

  // If the swarm is empty after potential resize, we might still technically proceed
  // but ReadSwarmField might do nothing or error depending on implementation.
  // It's safer to just return if empty.
  if (nGlobal == 0) {
      LOG_ALLOW(GLOBAL, LOG_INFO, "ReadAllSwarmFields - Swarm is empty for timestep %d. Nothing to read.\n", ti);
      PetscFunctionReturn(0);
  }

  // 1) Read positions (required, should succeed now)
  ierr = ReadSwarmField(user, "position", ti, "dat");
  if (ierr == PETSC_ERR_FILE_OPEN) {
      // This shouldn't happen if PreCheckAndResizeSwarm ran successfully, but handle defensively.
      LOG_ALLOW(GLOBAL, LOG_ERROR, "ReadAllSwarmFields - Position file missing for step %d, but PreCheck succeeded? Inconsistent state.\n", ti);
      PetscFunctionReturn(PETSC_ERR_PLIB);
  } else if (ierr) {
      LOG_ALLOW(GLOBAL, LOG_ERROR, "ReadAllSwarmFields - Failed to read position field for step %d (Error %d), even though size should match.\n", ti, ierr);
      CHKERRQ(ierr); // Treat other errors as fatal
  } else {
       LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadAllSwarmFields - Successfully read position field for step %d.\n", ti);
  }

  // 2) Read velocity
  ierr = CheckDataFile(ti, "velocity", "dat", &fileExists);
  if(fileExists) ierr = ReadSwarmField(user, "velocity", ti, "dat");
   if (ierr == PETSC_ERR_FILE_OPEN) { // Handle potentially missing file as warning
        LOG_ALLOW(GLOBAL, LOG_WARNING, "ReadAllSwarmFields - Velocity file for step %d not found. Velocity data may be invalid/uninitialized.\n", ti);
   } else if (ierr) {
       LOG_ALLOW(GLOBAL, LOG_ERROR, "ReadAllSwarmFields - Failed to read velocity field for step %d (Error %d). File might be corrupt or size incorrect despite pre-check.\n", ti, ierr);
       CHKERRQ(ierr); // Treat other errors as fatal
   } else {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadAllSwarmFields - Successfully read velocity field for step %d.\n", ti);
   }

  // 3) Read pos_phy
  ierr = CheckDataFile(ti, "pos_phy", "dat", &fileExists); 
  if(fileExists) ierr = ReadSwarmField(user, "pos_phy", ti, "dat");
   if (ierr == PETSC_ERR_FILE_OPEN) { // Handle potentially missing file as warning
        LOG_ALLOW(GLOBAL, LOG_WARNING, "ReadAllSwarmFields - Physical position file for step %d not found. physical position data may be invalid/uninitialized.\n", ti);
   } else if (ierr) {
       LOG_ALLOW(GLOBAL, LOG_ERROR, "ReadAllSwarmFields - Failed to read physical position field for step %d (Error %d). File might be corrupt or size incorrect despite pre-check.\n", ti, ierr);
       CHKERRQ(ierr); // Treat other errors as fatal
   } else {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadAllSwarmFields - Successfully read physical position field for step %d.\n", ti);
   }
   
  // 3) Read CellID (Optional)
  // ierr = CheckDataFile(ti, "velocity", ext, &fileExists); 
  // if(fileExists)ierr = ReadSwarmField(user, "DMSwarm_CellID", ti, "dat");
  // if (ierr == PETSC_ERR_FILE_OPEN) { LOG_ALLOW(GLOBAL, LOG_WARNING, "ReadAllSwarmFields - CellID file for step %d not found.\n", ti); }
  // else CHKERRQ(ierr); // Treat other errors as fatal unless CellID read can also fail gracefully

  // 4) Read weight (Optional)
  // ierr = ReadSwarmField(user, "weight", ti, "dat");
  // if (ierr == PETSC_ERR_FILE_OPEN) { LOG_ALLOW(GLOBAL, LOG_WARNING, "ReadAllSwarmFields - Weight file for step %d not found.\n", ti); }
  // else CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL, LOG_INFO, "ReadAllSwarmFields - Finished reading DMSwarm fields for timestep %d.\n", ti);
  PetscFunctionReturn(0);
}

/**
 * @brief Reads coordinate data (for particles)  from file into a PETSc Vec, then gathers it to rank 0.
 *
 * This function uses \c ReadFieldData to fill a PETSc Vec with coordinate data,
 * then leverages \c VecToArrayOnRank0 to gather that data into a contiguous array
 * (valid on rank 0 only).
 *
 * @param[in]  timeIndex    The time index used to construct file names.
 * @param[in]  user         Pointer to the user context.
 * @param[out] coordsArray  On rank 0, will point to a newly allocated array holding the coordinates.
 * @param[out] Ncoords      On rank 0, the length of \p coordsArray. On other ranks, 0.
 *
 * @return PetscErrorCode  Returns 0 on success, or non-zero on failures.
 */
PetscErrorCode ReadPositionsFromFile(PetscInt timeIndex,
                                      UserCtx *user,
                                      double **coordsArray,
                                      PetscInt *Ncoords)
{
  PetscFunctionBeginUser;

  PetscErrorCode ierr;
  Vec            coordsVec;

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadPositionsFromFile - Creating coords Vec.\n");
  ierr = VecCreate(PETSC_COMM_WORLD, &coordsVec);CHKERRQ(ierr);
  ierr = VecSetFromOptions(coordsVec);CHKERRQ(ierr);

  // For example: "position" is the name of the coordinate data
  ierr = ReadFieldData(user, "position", coordsVec, timeIndex, "dat");
  if (ierr) {
    LOG_ALLOW(GLOBAL, LOG_ERROR,
              "ReadPositionsFromFile - Error reading position data (ti=%d).\n",
              timeIndex);
    PetscFunctionReturn(ierr);
  }

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadPositions - Gathering coords Vec to rank 0.\n");
  ierr = VecToArrayOnRank0(coordsVec, Ncoords, coordsArray);CHKERRQ(ierr);

  ierr = VecDestroy(&coordsVec);CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "ReadPositionsFromFile - Successfully gathered coordinates. Ncoords=%d.\n", *Ncoords);
  PetscFunctionReturn(0);
}


/**
 * @brief Reads a named field from file into a PETSc Vec, then gathers it to rank 0.
 *
 * This function wraps \c ReadFieldData and \c VecToArrayOnRank0 into a single step.
 * The gathered data is stored in \p scalarArray on rank 0, with its length in \p Nscalars.
 *
 * @param[in]  timeIndex     The time index used to construct file names.
 * @param[in]  fieldName     Name of the field to be read (e.g., "velocity").
 * @param[in]  user          Pointer to the user context.
 * @param[out] scalarArray   On rank 0, a newly allocated array holding the field data.
 * @param[out] Nscalars      On rank 0, length of \p scalarArray. On other ranks, 0.
 *
 * @return PetscErrorCode  Returns 0 on success, or non-zero on failures.
 */
PetscErrorCode ReadFieldDataToRank0(PetscInt timeIndex,
                                           const char *fieldName,
                                           UserCtx *user,
                                           double **scalarArray,
                                           PetscInt *Nscalars)
{
  PetscFunctionBeginUser;

  PetscErrorCode ierr;
  Vec            fieldVec;

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadFieldDataWrapper - Creating field Vec.\n");
  ierr = VecCreate(PETSC_COMM_WORLD, &fieldVec);CHKERRQ(ierr);
  ierr = VecSetFromOptions(fieldVec);CHKERRQ(ierr);

  ierr = ReadFieldData(user, fieldName, fieldVec, timeIndex, "dat");
  if (ierr) {
    LOG_ALLOW(GLOBAL, LOG_ERROR,
              "ReadFieldDataWrapper - Error reading field '%s' (ti=%d).\n",
              fieldName, timeIndex);
    PetscFunctionReturn(ierr);
  }

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ReadFieldDataWrapper - Gathering field Vec to rank 0.\n");
  ierr = VecToArrayOnRank0(fieldVec, Nscalars, scalarArray);CHKERRQ(ierr);

  ierr = VecDestroy(&fieldVec);CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "ReadFieldDataWrapper - Successfully gathered field '%s'. Nscalars=%d.\n",
            fieldName, *Nscalars);
  PetscFunctionReturn(0);
}

/**
 * @brief Displays a structured banner summarizing the simulation configuration.
 *
 * This function prints key simulation parameters to standard output.
 * It is intended to be called ONLY by MPI rank 0.
 * It retrieves global domain bounds from `user->global_domain_bbox` and boundary
 * conditions for all faces from `user->face_bc_types`.
 *
 * @param[in] user         Pointer to `UserCtx` structure.
 * @param[in] StartTime    Initial simulation time.
 * @param[in] StartStep    Starting timestep index.
 * @param[in] StepsToRun   Total number of timesteps to run.
 * @param[in] num_mpi_procs Total number of MPI processes.
 * @param[in] total_num_particles Total number of particles.
 * @param[in] bboxlist     (If rank 0 needed to compute global_domain_bbox here, otherwise NULL)
 *
 * @return PetscErrorCode  Returns `0` on success.
 */
PetscErrorCode DisplayBanner(UserCtx *user,
                             PetscReal StartTime,
                             PetscInt StartStep,
                             PetscInt StepsToRun,
                             PetscMPIInt num_mpi_procs,
                             PetscInt total_num_particles,
                             BoundingBox *bboxlist_on_rank0) // bboxlist is only valid on rank 0
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    Cmpnts         global_min_coords, global_max_coords;

    PetscFunctionBeginUser;

    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "DisplayBanner - UserCtx pointer is NULL.");
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    if (rank == 0) {
        // If global_domain_bbox is not pre-populated in UserCtx, compute it here from bboxlist_on_rank0
        // This assumes bboxlist_on_rank0 is valid and contains all local bounding boxes on rank 0.
        if (bboxlist_on_rank0 && num_mpi_procs > 0) {
            global_min_coords = bboxlist_on_rank0[0].min_coords;
            global_max_coords = bboxlist_on_rank0[0].max_coords;
            for (PetscMPIInt p = 1; p < num_mpi_procs; ++p) {
                global_min_coords.x = PetscMin(global_min_coords.x, bboxlist_on_rank0[p].min_coords.x);
                global_min_coords.y = PetscMin(global_min_coords.y, bboxlist_on_rank0[p].min_coords.y);
                global_min_coords.z = PetscMin(global_min_coords.z, bboxlist_on_rank0[p].min_coords.z);
                global_max_coords.x = PetscMax(global_max_coords.x, bboxlist_on_rank0[p].max_coords.x);
                global_max_coords.y = PetscMax(global_max_coords.y, bboxlist_on_rank0[p].max_coords.y);
                global_max_coords.z = PetscMax(global_max_coords.z, bboxlist_on_rank0[p].max_coords.z);
            }
            // Optionally store this in user->global_domain_bbox if it's useful elsewhere
            // user->global_domain_bbox.min_coords = global_min_coords;
            // user->global_domain_bbox.max_coords = global_max_coords;
        } else {
            // Fallback or warning if bboxlist is not available for global calculation
            LOG_ALLOW(PETSC_COMM_SELF, LOG_WARNING, "DisplayBanner (Rank 0) - bboxlist not provided or num_mpi_procs <=0; using user->bbox for domain bounds.\n");
            global_min_coords = user->bbox.min_coords; // Use local bbox of rank 0 as fallback
            global_max_coords = user->bbox.max_coords;
        }


        ierr = PetscPrintf(PETSC_COMM_SELF, "\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "=============================================================\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "                          CASE SUMMARY                       \n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "=============================================================\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Grid Dimensions (Cells)     : %d X %d X %d\n", user->IM, user->JM, user->KM); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Global Domain Bounds (X)    : %.6f to %.6f\n", (double)global_min_coords.x, (double)global_max_coords.x); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Global Domain Bounds (Y)    : %.6f to %.6f\n", (double)global_min_coords.y, (double)global_max_coords.y); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Global Domain Bounds (Z)    : %.6f to %.6f\n", (double)global_min_coords.z, (double)global_max_coords.z); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "-------------------- Boundary Conditions --------------------\n"); CHKERRQ(ierr);
	const int face_name_width = 17; // Adjusted for longer names (Zeta,Eta,Xi)
        for (PetscInt i_face = 0; i_face < 6; ++i_face) {
	  BCFace current_face = (BCFace)i_face;
	  // The BCFaceToString will now return the Xi, Eta, Zeta versions
	  const char* face_str = BCFaceToString(current_face); 
	  const char* bc_type_str = BCTypeToString(user->boundary_faces[current_face].mathematical_type);
            
	  ierr = PetscPrintf(PETSC_COMM_SELF, " Face %-*s : %s\n", 
			     face_name_width, face_str, bc_type_str); CHKERRQ(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "-------------------------------------------------------------\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Start Time                  : %.4f\n", (double)StartTime); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Timestep Size               : %.4f\n", (double)user->dt); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Starting Step               : %d\n", StartStep); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Total Steps to Run          : %d\n", StepsToRun); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Number of MPI Processes     : %d\n", num_mpi_procs); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD," Number of Particles         : %d\n", total_num_particles); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, " Particle Initialization Mode: %d\n", user->ParticleInitialization); CHKERRQ(ierr);
        if (user->ParticleInitialization == 0) {
            if (user->inletFaceDefined) {
                ierr = PetscPrintf(PETSC_COMM_SELF, " Particles Initialized At    : %s (Enum Val: %d)\n", BCFaceToString(user->identifiedInletBCFace), user->identifiedInletBCFace); CHKERRQ(ierr);
            } else {
                ierr = PetscPrintf(PETSC_COMM_SELF, " Particles Initialized At    : --- (No INLET face identified)\n"); CHKERRQ(ierr);
            }
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, " Field Initialization Mode   : %d\n", user->FieldInitialization); CHKERRQ(ierr);
        if (user->FieldInitialization == 1) {
	  ierr = PetscPrintf(PETSC_COMM_SELF, " Constant Velocity           : x - %.4f, y - %.4f, z - %.4f \n", (double)user->InitialConstantContra.x,(double)user->InitialConstantContra.y,(double)user->InitialConstantContra.z ); CHKERRQ(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "=============================================================\n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "\n"); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Parses the specified boundary conditions file to identify BC types for all 6 global faces
 *        and specifically notes the first "INLET" face found.
 *
 * Reads the file specified by `bcs_input_filename`. The file is expected to specify
 * boundary conditions for the six logical faces of the global computational domain
 * (-Xi, +Xi, -Eta, +Eta, -Zeta, +Zeta).
 * - Populates `user->face_bc_types[BCFace]` with the `BCType` for each face.
 * - Sets `user->inletFaceDefined` to `PETSC_TRUE` if at least one "INLET" is found.
 * - Sets `user->identifiedInletBCFace` to the `BCFace` of the *first* "INLET" found.
 *
 * @param[in,out] user                 Pointer to the UserCtx structure where parsed BC info will be stored.
 * @param[in]     bcs_input_filename   The name/path of the boundary conditions file to parse.
 * @return PetscErrorCode 0 on success, error code on failure.
 */
/*
PetscErrorCode ParseAllBoundaryConditions(UserCtx *user, const char *bcs_input_filename) { // Added const char *bcs_input_filename
    PetscErrorCode ierr = 0;
    FILE           *file;
    char           line[256];
    char           face_specifier_str[10];
    char           bc_type_str[100];
    PetscMPIInt    rank;

    BCType    parsed_face_bc_types_rank0[6];
    BCFace    parsed_inlet_face_rank0 = BC_FACE_NEG_X; 
    PetscBool inlet_defined_on_rank0 = PETSC_FALSE;
    PetscInt  i;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Initialize UserCtx fields (these will be set by rank 0 and then broadcasted)
    // These initializations on non-root ranks are temporary until broadcast.
    for (i = 0; i < 6; ++i) { user->face_bc_types[i] = WALL; }
    user->inletFaceDefined = PETSC_FALSE;
    user->identifiedInletBCFace = BC_FACE_NEG_X; // A default valid BCFace

    if (rank == 0) {
        // Initialize rank 0's temporary parsing buffers
        for (i = 0; i < 6; ++i) { parsed_face_bc_types_rank0[i] = WALL; }
        // inlet_defined_on_rank0 and parsed_inlet_face_rank0 already have defaults.

        LOG_ALLOW(GLOBAL, LOG_INFO, "ParseAllBoundaryConditions (Rank 0) - Attempting to open BCs file: %s\n", bcs_input_filename);
        file = fopen(bcs_input_filename, "r"); // Use filename passed as argument
        if (!file) {
            LOG_ALLOW(GLOBAL, LOG_WARNING, "ParseAllBoundaryConditions (Rank 0) - Could not open BCs file '%s'. All faces will default to WALL, no INLET identified.\n", bcs_input_filename);
        } else {
            LOG_ALLOW(GLOBAL, LOG_INFO, "ParseAllBoundaryConditions (Rank 0) - Reading BCs file: %s\n", bcs_input_filename);
            if (fgets(line, sizeof(line), file) != NULL) {
                line[strcspn(line, "\r\n")] = 0;
                if (strcmp(line, "BCS") != 0) { fseek(file, 0, SEEK_SET); }
            }

            PetscInt inlet_count_temp = 0; // To track multiple INLET definitions

            while (fgets(line, sizeof(line), file) != NULL) {
                if (sscanf(line, "%9s %99s", face_specifier_str, bc_type_str) == 2) {
                    BCFace current_face_enum_temp; 
                    PetscBool valid_face_specifier = PETSC_TRUE;

                    if      (strcmp(face_specifier_str, "-Xi") == 0) current_face_enum_temp = BC_FACE_NEG_X;
                    else if (strcmp(face_specifier_str, "+Xi") == 0) current_face_enum_temp = BC_FACE_POS_X;
                    else if (strcmp(face_specifier_str, "-Eta") == 0) current_face_enum_temp = BC_FACE_NEG_Y;
                    else if (strcmp(face_specifier_str, "+Eta") == 0) current_face_enum_temp = BC_FACE_POS_Y;
                    else if (strcmp(face_specifier_str, "-Zeta") == 0) current_face_enum_temp = BC_FACE_NEG_Z;
                    else if (strcmp(face_specifier_str, "+Zeta") == 0) current_face_enum_temp = BC_FACE_POS_Z;
                    else {
                        LOG_ALLOW(GLOBAL, LOG_WARNING, "ParseAllBoundaryConditions (Rank 0) - Unknown face specifier '%s' in %s. Line ignored.\n", face_specifier_str, bcs_input_filename);
                        valid_face_specifier = PETSC_FALSE;
                    }

                    if (valid_face_specifier) {
                        BCType current_bc_type_enum_temp = WALL; // Default for this line
                        for (i = 0; bc_type_str[i]; i++) { bc_type_str[i] = toupper(bc_type_str[i]); }

                        if      (strcmp(bc_type_str, "DIRICHLET") == 0) current_bc_type_enum_temp = DIRICHLET;
                        else if (strcmp(bc_type_str, "NEUMANN")   == 0) current_bc_type_enum_temp = NEUMANN;
                        else if (strcmp(bc_type_str, "WALL")      == 0) current_bc_type_enum_temp = WALL;
                        else if (strcmp(bc_type_str, "INLET")     == 0) current_bc_type_enum_temp = INLET;
                        else if (strcmp(bc_type_str, "OUTLET")    == 0) current_bc_type_enum_temp = OUTLET;
                        // else if (strcmp(bc_type_str, "CUSTOM") == 0) current_bc_type_enum_temp = CUSTOM;
                        else {
                            LOG_ALLOW(GLOBAL, LOG_WARNING, "ParseAllBoundaryConditions (Rank 0) - Unknown BC type '%s' for face '%s' in %s. Defaulting to WALL.\n", bc_type_str, face_specifier_str, bcs_input_filename);
                            // current_bc_type_enum_temp remains WALL
                        }
                        
                        parsed_face_bc_types_rank0[current_face_enum_temp] = current_bc_type_enum_temp;

                        if (current_bc_type_enum_temp == INLET) {
                            inlet_count_temp++;
                            if (!inlet_defined_on_rank0) { // Store the first INLET found
                                inlet_defined_on_rank0 = PETSC_TRUE;
                                parsed_inlet_face_rank0 = current_face_enum_temp;
                                LOG_ALLOW(GLOBAL, LOG_INFO, "ParseAllBoundaryConditions (Rank 0) - Identified FIRST INLET face: %s (BCFace Enum: %d) from %s\n", face_specifier_str, parsed_inlet_face_rank0, bcs_input_filename);
                            } else if (inlet_count_temp > 1) {
                                LOG_ALLOW(GLOBAL, LOG_WARNING, "ParseAllBoundaryConditions (Rank 0) - Additional INLET face '%s' found. Only the first one ('%d') will be used for particle init.\n", face_specifier_str, parsed_inlet_face_rank0);
                            }
                        }
                    }
                }
            }
            fclose(file);
            if (!inlet_defined_on_rank0 && inlet_count_temp == 0) { // Check inlet_count_temp as well
                 LOG_ALLOW(GLOBAL, LOG_WARNING, "ParseAllBoundaryConditions (Rank 0) - No INLET face defined in %s.\n", bcs_input_filename);
            }
        } // end if file opened

        // Update UserCtx directly on rank 0 with the parsed values
        for(i=0; i<6; ++i) {
            user->face_bc_types[i] = parsed_face_bc_types_rank0[i];
        }
        user->inletFaceDefined = inlet_defined_on_rank0;
        if (inlet_defined_on_rank0) {
            user->identifiedInletBCFace = parsed_inlet_face_rank0;
        }
    } // end if rank == 0

    // Broadcast the UserCtx fields that were populated by rank 0
    ierr = MPI_Bcast(user->face_bc_types, 6 * sizeof(BCType), MPI_BYTE, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    
    PetscMPIInt temp_inlet_defined_bcast = (PetscMPIInt)user->inletFaceDefined; // Use user struct directly on rank 0
    PetscInt    temp_inlet_face_bcast    = (PetscInt)user->identifiedInletBCFace;  // Use user struct directly on rank 0

    ierr = MPI_Bcast(&temp_inlet_defined_bcast, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Bcast(&temp_inlet_face_bcast, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

    if (rank != 0) { // Non-root ranks update from broadcasted values (face_bc_types already updated by Bcast above)
        user->inletFaceDefined = (PetscBool)temp_inlet_defined_bcast;
        if (user->inletFaceDefined) {
            user->identifiedInletBCFace = (BCFace)temp_inlet_face_bcast;
        }
    }

    // Log final status on all ranks
    for(i=0; i<6; ++i) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "ParseAllBoundaryConditions - Rank %d: Face %d (%s) has BC Type %d (%s)\n",
                  rank, i, BCFaceToString((BCFace)i), user->face_bc_types[i], BCTypeToString(user->face_bc_types[i]));
    }
    if (user->inletFaceDefined) {
        LOG_ALLOW(LOCAL, LOG_INFO, "ParseAllBoundaryConditions - Rank %d confirmed INLET face for particle init: %d (%s) from %s.\n",
                  rank, user->identifiedInletBCFace, BCFaceToString(user->identifiedInletBCFace), bcs_input_filename);
    } else if (user->ParticleInitialization == 0) { // Only warn if Mode 0 and no inlet was found
        LOG_ALLOW(LOCAL, LOG_WARNING, "ParseAllBoundaryConditions - Rank %d: Mode 0 particle init, but NO INLET face identified from %s.\n",
                  rank, bcs_input_filename);
    }

    PetscFunctionReturn(0);
}
*/

//================================================================================
//
//                        PRIVATE HELPER FUNCTIONS
//
//================================================================================

/**
 * @brief Frees the memory allocated for a linked list of BC_Param structs.
 * @param head A pointer to the head of the linked list to be freed.
 */
void FreeBC_ParamList(BC_Param *head) {
    BC_Param *current = head;
    while (current != NULL) {
        BC_Param *next = current->next;
        PetscFree(current->key);
        PetscFree(current->value);
        PetscFree(current);
        current = next;
    }
}

/**
 * @brief Converts a string representation of a face to a BCFace enum.
 * @param str The input string (e.g., "-Xi", "+Zeta"). Case-insensitive.
 * @param[out] face_out The resulting BCFace enum.
 * @return 0 on success.
 */
PetscErrorCode StringToBCFace(const char* str, BCFace* face_out) {
    if      (strcasecmp(str, "-Xi")   == 0) *face_out = BC_FACE_NEG_X;
    else if (strcasecmp(str, "+Xi")   == 0) *face_out = BC_FACE_POS_X;
    else if (strcasecmp(str, "-Eta")  == 0) *face_out = BC_FACE_NEG_Y;
    else if (strcasecmp(str, "+Eta")  == 0) *face_out = BC_FACE_POS_Y;
    else if (strcasecmp(str, "-Zeta") == 0) *face_out = BC_FACE_NEG_Z;
    else if (strcasecmp(str, "+Zeta") == 0) *face_out = BC_FACE_POS_Z;
    else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Unknown face specifier: %s", str);
    return 0;
}

/**
 * @brief Converts a string representation of a BC type to a BCType enum.
 * @param str The input string (e.g., "WALL", "INLET"). Case-insensitive.
 * @param[out] type_out The resulting BCType enum.
 * @return 0 on success.
 */
PetscErrorCode StringToBCType(const char* str, BCType* type_out) {
    if      (strcasecmp(str, "WALL")      == 0) *type_out = WALL;
    else if (strcasecmp(str, "SYMMETRY")  == 0) *type_out = SYMMETRY;
    else if (strcasecmp(str, "INLET")     == 0) *type_out = INLET;
    else if (strcasecmp(str, "OUTLET")    == 0) *type_out = OUTLET;
    else if (strcasecmp(str, "NOGRAD")    == 0) *type_out = NOGRAD;
    // ... add other BCTypes here ...
    else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Unknown BC Type string: %s", str);
    return 0;
}

/**
 * @brief Converts a string representation of a handler to a BCHandlerType enum.
 * @param str The input string (e.g., "noslip", "constant_velocity"). Case-insensitive.
 * @param[out] handler_out The resulting BCHandlerType enum.
 * @return 0 on success.
 */
PetscErrorCode StringToBCHandlerType(const char* str, BCHandlerType* handler_out) {
    if      (strcasecmp(str, "noslip")              == 0) *handler_out = BC_HANDLER_WALL_NOSLIP;
    else if (strcasecmp(str, "constant_velocity")   == 0) *handler_out = BC_HANDLER_INLET_CONSTANT_VELOCITY;
    else if (strcasecmp(str, "conservation")        == 0) *handler_out = BC_HANDLER_OUTLET_CONSERVATION;
    else if (strcasecmp(str, "allcopy")             == 0) *handler_out = BC_HANDLER_NOGRAD_COPY_GHOST;
    // ... add other BCHandlerTypes here ...
    else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Unknown BC Handler string: %s", str);
    return 0;
}

/**
 * @brief Validates that a specific handler is compatible with a general BC type.
 * @param type The general BCType.
 * @param handler The specific BCHandlerType.
 * @return 0 if compatible, error code otherwise.
 */
PetscErrorCode ValidateBCHandlerForBCType(BCType type, BCHandlerType handler) {
    switch (type) {
        case NOGRAD:
	    if(handler != BC_HANDLER_NOGRAD_COPY_GHOST) return PETSC_ERR_ARG_WRONG;
	    break;
        case WALL:
            if (handler != BC_HANDLER_WALL_NOSLIP && handler != BC_HANDLER_WALL_MOVING) return PETSC_ERR_ARG_WRONG;
            break;
        case INLET:
            if (handler != BC_HANDLER_INLET_CONSTANT_VELOCITY && handler != BC_HANDLER_INLET_PULSANTILE_FLUX) return PETSC_ERR_ARG_WRONG;
            break;
        // ... add other validation cases here ...
        default: break;
    }
    return 0; // Combination is valid
}

//================================================================================
//
//                        PUBLIC PARSING FUNCTION
//
//================================================================================

/**
 * @brief Parses the boundary conditions file to configure the type, handler, and
 *        any associated parameters for all 6 global faces of the domain.
 *
 * This function performs the following steps:
 * 1.  On MPI rank 0, it reads the specified configuration file line-by-line.
 * 2.  It parses each line for `<Face> <Type> <Handler> [param=value]...` format.
 * 3.  It validates the parsed strings and stores the configuration, including a
 *     linked list of parameters, in a temporary array.
 * 4.  It then serializes this configuration and broadcasts it to all other MPI ranks.
 * 5.  All ranks (including rank 0) then deserialize the broadcasted data to populate
 *     their local `user->boundary_faces` array identically.
 * 6.  It also sets legacy fields in UserCtx for compatibility with other modules.
 *
 * @param[in,out] user               The main UserCtx struct where the final configuration
 *                                   for all ranks will be stored.
 * @param[in]     bcs_input_filename The path to the boundary conditions configuration file.
 * @return PetscErrorCode 0 on success, error code on failure.
 */
PetscErrorCode ParseAllBoundaryConditions(UserCtx *user, const char *bcs_input_filename)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;

    // Temporary storage for rank 0 to build the configuration before broadcasting.
    BoundaryFaceConfig configs_rank0[6];
    
    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    if (rank == 0) {
        FILE *file;
        char line_buffer[1024];

        // Initialize the temporary config array with safe defaults on rank 0.
        for (int i = 0; i < 6; i++) {
            configs_rank0[i].face_id = (BCFace)i;
            configs_rank0[i].mathematical_type = WALL;
            configs_rank0[i].handler_type = BC_HANDLER_WALL_NOSLIP;
            configs_rank0[i].params = NULL;
            configs_rank0[i].handler = NULL; // Handler object is not created here.
        }

        LOG_ALLOW(GLOBAL, LOG_INFO, "Parsing BC configuration from '%s' on rank 0... \n", bcs_input_filename);
        file = fopen(bcs_input_filename, "r");
        if (!file) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Could not open BCs file '%s'.", bcs_input_filename);

        while (fgets(line_buffer, sizeof(line_buffer), file)) {
            char *current_pos = line_buffer;
            while (isspace((unsigned char)*current_pos)) current_pos++; // Skip leading whitespace
            if (*current_pos == '#' || *current_pos == '\0' || *current_pos == '\n' || *current_pos == '\r') continue;

            char *face_str = strtok(current_pos, " \t\n\r");
            char *type_str = strtok(NULL, " \t\n\r");
            char *handler_str = strtok(NULL, " \t\n\r");

            if (!face_str || !type_str || !handler_str) {
                LOG_ALLOW(GLOBAL, LOG_WARNING, "Malformed line in bcs.dat, skipping: %s \n", line_buffer);
                continue;
            }

            BCFace      face_enum;
            BCType      type_enum;
            BCHandlerType handler_enum;
            const char* handler_name_for_log;

            // --- Convert strings to enums and validate ---
            ierr = StringToBCFace(face_str, &face_enum); CHKERRQ(ierr);
            ierr = StringToBCType(type_str, &type_enum); CHKERRQ(ierr);
            ierr = StringToBCHandlerType(handler_str, &handler_enum); CHKERRQ(ierr);
            ierr = ValidateBCHandlerForBCType(type_enum, handler_enum);
            if (ierr) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Validation failed: Handler '%s' is not valid for Type '%s' on Face '%s'.\n", handler_str, type_str, face_str);

            // Store the core types for the corresponding face
            configs_rank0[face_enum].mathematical_type = type_enum;
            configs_rank0[face_enum].handler_type = handler_enum;
            handler_name_for_log = BCHandlerTypeToString(handler_enum); // Assumes this utility exists
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "  Parsed Face '%s': Type=%s, Handler=%s \n", face_str, type_str, handler_name_for_log);

            // --- Parse optional key=value parameters for this face ---
            FreeBC_ParamList(configs_rank0[face_enum].params); // Clear any previous (default) params
            configs_rank0[face_enum].params = NULL;
            BC_Param **param_next_ptr = &configs_rank0[face_enum].params; // Pointer to the 'next' pointer to build the list

            char* token;
            while ((token = strtok(NULL, " \t\n\r")) != NULL) {
                char* equals_ptr = strchr(token, '=');
                if (!equals_ptr) {
                    LOG_ALLOW(GLOBAL, LOG_WARNING, "Malformed parameter '%s' on face '%s', skipping. \n", token, face_str);
                    continue;
                }
                
                *equals_ptr = '\0'; // Temporarily split the string at '=' to separate key and value
                char* key_str = token;
                char* value_str = equals_ptr + 1;

                BC_Param *new_param;
                ierr = PetscMalloc1(1, &new_param); CHKERRQ(ierr);
                ierr = PetscStrallocpy(key_str, &new_param->key); CHKERRQ(ierr);
                ierr = PetscStrallocpy(value_str, &new_param->value); CHKERRQ(ierr);
                new_param->next = NULL;

                *param_next_ptr = new_param;
                param_next_ptr = &new_param->next;
                LOG_ALLOW(GLOBAL, LOG_DEBUG, "    - Found param: [%s] = [%s] \n", new_param->key, new_param->value);
            }
        }
        fclose(file);
    }

    // =========================================================================
    //               BROADCASTING THE CONFIGURATION FROM RANK 0
    // =========================================================================
    // This is a critical step to ensure all processes have the same configuration.
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d broadcasting/receiving BC configuration.\n", rank);

    for (int i = 0; i < 6; i++) {
        // --- Broadcast simple enums ---
        if (rank == 0) {
            user->boundary_faces[i] = configs_rank0[i]; // Rank 0 populates its final struct
        }
        MPI_Bcast(&user->boundary_faces[i].mathematical_type, 1, MPI_INT, 0, PETSC_COMM_WORLD);
        MPI_Bcast(&user->boundary_faces[i].handler_type, 1, MPI_INT, 0, PETSC_COMM_WORLD);
        
        // --- Serialize and Broadcast the parameter linked list ---
        PetscInt n_params = 0;
        if (rank == 0) { // On rank 0, count the number of parameters to send
            for (BC_Param *p = user->boundary_faces[i].params; p; p = p->next) n_params++;
        }
        MPI_Bcast(&n_params, 1, MPI_INT, 0, PETSC_COMM_WORLD);
        
        if (rank != 0) { // Non-root ranks need to receive and build the list
            FreeBC_ParamList(user->boundary_faces[i].params); // Ensure list is empty before building
            user->boundary_faces[i].params = NULL;
        }

        BC_Param **param_next_ptr = &user->boundary_faces[i].params;

        for (int j = 0; j < n_params; j++) {
            char key_buf[256] = {0}, val_buf[256] = {0};
            if (rank == 0) {
                // On rank 0, navigate to the j-th param and copy its data to buffers
                BC_Param *p = user->boundary_faces[i].params;
                for (int k = 0; k < j; k++) p = p->next;
                strncpy(key_buf, p->key, 255);
                strncpy(val_buf, p->value, 255);
            }

            MPI_Bcast(key_buf, 256, MPI_CHAR, 0, PETSC_COMM_WORLD);
            MPI_Bcast(val_buf, 256, MPI_CHAR, 0, PETSC_COMM_WORLD);

            if (rank != 0) {
                // On non-root ranks, deserialize: create a new node and append it
                BC_Param *new_param;
                ierr = PetscMalloc1(1, &new_param); CHKERRQ(ierr);
                ierr = PetscStrallocpy(key_buf, &new_param->key); CHKERRQ(ierr);
                ierr = PetscStrallocpy(val_buf, &new_param->value); CHKERRQ(ierr);
                new_param->next = NULL;
                *param_next_ptr = new_param;
                param_next_ptr = &new_param->next;
            } else {
                 // On rank 0, just advance the pointer for the next iteration
                 param_next_ptr = &((*param_next_ptr)->next);
            }
        }
        user->boundary_faces[i].face_id = (BCFace)i; // Ensure face_id is set on all ranks
    }
    
    // --- Set legacy fields for compatibility with particle system ---
    user->inletFaceDefined = PETSC_FALSE;
    for (int i=0; i<6; i++) {
        if (user->boundary_faces[i].mathematical_type == INLET) {
            user->inletFaceDefined = PETSC_TRUE;
            user->identifiedInletBCFace = (BCFace)i;
            LOG_ALLOW(GLOBAL, LOG_INFO, "Inlet face for particle initialization identified as Face %d.\n", i);
            break; // Found the first one, stop looking
        }
    }

    if (rank == 0) {
        // Rank 0 can now free the linked lists it created, as they have been broadcast.
        // Or, if user->boundary_faces was used directly, this is not needed.
    }
    
    PetscFunctionReturn(0);
}
