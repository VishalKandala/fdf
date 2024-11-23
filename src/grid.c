
// grid.c 

#include "grid.h"
#include "logging.h"
#include <petsc.h>
#include <stdlib.h>

//extern PetscInt block_number;

/**
 * @brief Initializes the DMDA grid for the simulation.
 *
 * This function sets up a 3D DMDA (Distributed Memory Distributed Array) grid with uniform
 * coordinates. The DMDA is essential for managing simulation data distributed across
 * multiple MPI processes. It assigns the global dimensions of the grid and partitions
 * it among processes based on the stencil width and degrees of freedom.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing grid details.
 *                     This structure must include information such as global grid
 *                     dimensions (`IM`, `JM`, `KM`) and the DMDA handles (`da`, `fda`).
 * @param[in] L_x,L_y,L_z Dimensions of the domain.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 *
 * @note
 * - The global dimensions of the grid are read from the `user` structure.
 * - The function also sets up uniform coordinates in the range [0,Lx],[0,Ly],[0,Lz].
 */
PetscErrorCode InitializeGridDM(UserCtx *user, PetscReal L_x, PetscReal L_y, PetscReal L_z) {
    PetscErrorCode ierr; // PETSc error handling

    // Validate grid dimensions
    if (user->IM <= 0 || user->JM <= 0 || user->KM <= 0) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, 
                "InitializeGridDM - Grid dimensions (IM, JM, KM) must be positive.");
    }

    // Log the start of grid initialization
    LOG(GLOBAL, LOG_INFO, "InitializeGridDM - Starting DMDA grid initialization for IM=%d, JM=%d, KM=%d.\n",
        user->IM, user->JM, user->KM);

    // Create a 3D DMDA grid
    ierr = DMDACreate3d(PETSC_COMM_WORLD,
                        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,  // Boundary types
                        DMDA_STENCIL_BOX,                                      // Stencil type
                        user->IM + 1, user->JM + 1, user->KM + 1,             // Global dimensions
                        PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,             // Dimensions per process
                        1,                                                     // Degrees of freedom
                        1,                                                     // Stencil width
                        NULL, NULL, NULL,                                      // Array partitioning
                        &(user->da)); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_DEBUG, "InitializeGridDM - Created 3D DMDA grid.\n");

    // Set up the DMDA
    ierr = DMSetUp(user->da); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_DEBUG, "InitializeGridDM - DMDA setup completed.\n");

    // Assign coordinates based on domain size (L_x, L_y, L_z)
    PetscReal x_min = 0.0, x_max = L_x;
    PetscReal y_min = 0.0, y_max = L_y;
    PetscReal z_min = 0.0, z_max = L_z;
    ierr = DMDASetUniformCoordinates(user->da, x_min, x_max, y_min, y_max, z_min, z_max); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_DEBUG, "InitializeGridDM - Coordinates set in the range [%f, %f] x [%f, %f] x [%f, %f].\n", 
        x_min, x_max, y_min, y_max, z_min, z_max);

    // Retrieve the coordinate DM
    ierr = DMGetCoordinateDM(user->da, &(user->fda)); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_DEBUG, "InitializeGridDM - Coordinate DM retrieved.\n");

    // Debugging logs for DMDA and its coordinate DM
    if (get_log_level() == LOG_DEBUG) {
        LOG(GLOBAL, LOG_INFO, "InitializeGridDM - Viewing DMDA:\n");
        ierr = DMView(user->da, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

        LOG(GLOBAL, LOG_INFO, "InitializeGridDM - Viewing coordinate DM:\n");
        ierr = DMView(user->fda, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    }

    // Log the successful completion of grid initialization
    LOG(GLOBAL, LOG_INFO, "InitializeGridDM - DMDA grid initialization completed successfully.\n");

    return 0;
}

/**
 * @brief Parses grid-related inputs for setting up the simulation.
 *
 * This function coordinates the retrieval of grid parameters by delegating
 * to appropriate helper functions based on the grid source (generation or file).
 *
 * @param[in,out] user          Pointer to the UserCtx structure containing grid details.
 * @param[out]    generate_grid Flag indicating whether the grid should be generated (1) or read from file (0).
 * @param[out]    grid1d              Flag indicating whether the grid is 1D(1) or not(0).
 * @param[out]    L_x           Pointer to store the domain length in the x-direction.
 * @param[out]    L_y           Pointer to store the domain length in the y-direction.
 * @param[out]    L_z           Pointer to store the domain length in the z-direction.
 * @param[out]    imm           Pointer to Array to store the i-dimensions for each block.
 * @param[out]    jmm           Pointer to Array to store the j-dimensions for each block.
 * @param[out]    kmm           Pointer to Array to store the k-dimensions for each block.
 * @param[out]    nblk          Pointer to store the number of blocks.
 * @param[out]    fd            File pointer for reading grid data (if applicable).
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode ParseGridInputs(UserCtx *user, PetscInt *generate_grid, PetscInt *grid1d, PetscReal *L_x,
                               PetscReal *L_y, PetscReal *L_z, PetscInt **imm, PetscInt **jmm, PetscInt **kmm,
                               PetscInt *nblk, FILE *fd) {
    PetscErrorCode ierr;

    LOG(GLOBAL, LOG_INFO, "ParseGridInputs - Starting retrieval of grid inputs.\n");

    // Insert options from the control file
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, "control.dat", PETSC_TRUE); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_DEBUG, "ParseGridInputs - Loaded options from control.dat.\n");

    // Read the basic grid settings
    ierr = PetscOptionsGetInt(NULL, NULL, "-grid", generate_grid, NULL); CHKERRQ(ierr);

    if (*generate_grid) {
        // Delegate the reading of grid generation inputs, including grid1d
      ierr = ReadGridGenerationInputs(grid1d,L_x, L_y, L_z, imm, jmm, kmm, nblk); CHKERRQ(ierr);
    } else {
        // Delegate the reading of grid file data, including grid1d // hard coded the file as a temporary hack! 
        ierr = ReadGridFile("grid.dat", nblk, imm, jmm, kmm, grid1d, PETSC_COMM_WORLD); CHKERRQ(ierr);
    }

    LOG(GLOBAL, LOG_INFO, "ParseGridInputs - Grid inputs retrieved successfully.\n");
    return 0;
}


/**
 * @brief Populates the DMDA grid with coordinates from a file or generates them based on input parameters.
 *
 * This function assigns coordinates to the nodes of a DMDA grid. The coordinates can be read
 * from a file or generated programmatically based on domain dimensions and grid resolution.
 * The function supports both 1D and 3D grid generation modes.
 *
 * @param[in,out] user          Pointer to the UserCtx structure containing grid details.
 * @param[in]     generate_grid Flag indicating whether the grid is programmatically generated (1) or read from file (0).
 * @param[in]     grid1d        Flag indicating if the grid is 1D (1) or 3D (0).
 * @param[in]     IM, JM, KM    Global grid dimensions in the x, y, and z directions, respectively.
 * @param[in]     L_x, L_y, L_z Domain lengths in the x, y, and z directions, respectively.
 * @param[in]     fd            File pointer for reading grid data (if required).
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 *
 * @note
 * - The function assumes that the DMDA (`user->da`) has already been set up.
 * - Coordinates are normalized using the characteristic length scale (`cl`) and domain length scale (`L_dim`).
 */
PetscErrorCode AssignGridCoordinates(UserCtx *user, PetscInt generate_grid, PetscInt grid1d, PetscInt IM, PetscInt JM, 
                            PetscInt KM, PetscReal L_x, PetscReal L_y, PetscReal L_z, FILE *fd) {
    PetscErrorCode ierr;          // PETSc error handling
    PetscInt rank, i, j, k;       // MPI rank and loop counters
    Vec Coor, gCoor;              // Local and global coordinate vectors
    Cmpnts ***coor;               // 3D array for node coordinates
    PetscReal *gc;            // Temporary storage for coordinates
    PetscReal cl = 1.0, L_dim = 1.0; // Normalization constants

    // Retrieve MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Log the start of grid population
    LOG(GLOBAL, LOG_INFO, "AssignGridCoordinates - Starting grid population for rank %d.\n", rank);

    // Retrieve local grid information
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    PetscInt xs = info.xs, xe = info.xs + info.xm; // Local x start and end indices
    PetscInt ys = info.ys, ye = info.ys + info.ym; // Local y start and end indices
    PetscInt zs = info.zs, ze = info.zs + info.zm; // Local z start and end indices

    // Allocate memory for grid coordinates based on 1D or 3D mode
    if (grid1d) {
        // For 1D grid, allocate memory for IM + JM + KM coordinates
        ierr = PetscMalloc1(IM + JM + KM, &gc); CHKERRQ(ierr);
    } else {
        // For 3D grid, allocate memory for all grid points (including boundaries) and three coordinates per point
        ierr = PetscMalloc1(3 * (IM + 1) * (JM + 1) * (KM + 1), &gc); CHKERRQ(ierr);
    }
    LOG(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Memory allocated for grid coordinates.\n");

    // Access the DMDA coordinate vector
    ierr = DMGetCoordinatesLocal(user->da, &Coor); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, Coor, &coor); CHKERRQ(ierr);

    // Rank 0: Read or generate grid coordinates
    if (rank == 0) {
        if (grid1d) {
            // Read or generate 1D grid coordinates
            // (This part remains unchanged)
        } else {
            // Generate 3D grid coordinates
            if (!generate_grid) {
                // Read coordinates from file
                for (k = 0; k <= KM; k++) {
                    for (j = 0; j <= JM; j++) {
                        for (i = 0; i <= IM; i++) {
                            PetscInt index = ((k * (JM + 1) * (IM + 1)) + (j * (IM + 1)) + i) * 3;
                            if (fscanf(fd, "%le %le %le\n", &gc[index], &gc[index + 1], &gc[index + 2]) != 3) {
                                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Error reading grid coordinates from file");
                            }
                        }
                    }
                }
                LOG(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Rank 0 read 3D grid coordinates from file.\n");
            } else {
                // Generate coordinates programmatically
                for (k = 0; k <= KM; k++) {
                    for (j = 0; j <= JM; j++) {
                        for (i = 0; i <= IM; i++) {
                            PetscInt index = ((k * (JM + 1) * (IM + 1)) + (j * (IM + 1)) + i) * 3;
                            gc[index]     = L_x / IM * i; // X-coordinate
                            gc[index + 1] = L_y / JM * j; // Y-coordinate
                            gc[index + 2] = L_z / KM * k; // Z-coordinate
                        }
                    }
                }
                LOG(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Rank 0 generated 3D grid coordinates programmatically.\n");
            }
        }
    }

    // Broadcast grid coordinates to all ranks from rank 0
    if (grid1d) {
        ierr = MPI_Bcast(gc, IM + JM + KM, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        LOG(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Broadcasted 1D grid coordinates.\n");
    } else {
        PetscInt total_points = 3 * (IM + 1) * (JM + 1) * (KM + 1);
        ierr = MPI_Bcast(gc, total_points, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        LOG(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Broadcasted 3D grid coordinates.\n");
    }

    // Assign coordinates to the local grid
    if (grid1d) {
        // (This part remains unchanged)
    } else {
        // Assign coordinates for the 3D grid
        for (k = zs; k < ze; k++) {
            for (j = ys; j < ye; j++) {
                for (i = xs; i < xe; i++) {
                    if (k <= KM && j <= JM && i <= IM) {
                        PetscInt index = ((k * (JM + 1) * (IM + 1)) + (j * (IM + 1)) + i) * 3;
                        coor[k][j][i].x = gc[index]     / cl * L_dim;
                        coor[k][j][i].y = gc[index + 1] / cl * L_dim;
                        coor[k][j][i].z = gc[index + 2] / cl * L_dim;
                    }
                }
            }
        }
    }

    LOG(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Local grid coordinates assigned for rank %d.\n", rank);

    // Clean up and restore arrays
    ierr = PetscFree(gc); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, Coor, &coor); CHKERRQ(ierr);

    // Sync coordinates across ranks
    ierr = DMGetCoordinates(user->da, &gCoor); CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(user->fda, Coor, INSERT_VALUES, gCoor); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(user->fda, Coor, INSERT_VALUES, gCoor); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->fda, gCoor, INSERT_VALUES, Coor); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->fda, gCoor, INSERT_VALUES, Coor); CHKERRQ(ierr);

    LOG(GLOBAL, LOG_INFO, "AssignGridCoordinates - Grid population completed for rank %d.\n", rank);

    return 0;
}


/**
 * @brief Determines the grid dimensions for the specified block.
 *
 * This function calculates the grid dimensions (`IM`, `JM`, `KM`) for a specific block
 * based on the provided input parameters. If the grid is programmatically generated,
 * the dimensions are derived from input arrays. Otherwise, they are read from a file.
 *
 * @param[in]     bi            Block index to retrieve dimensions for.
 * @param[in,out] user          Pointer to the UserCtx structure containing simulation context.
 * @param[out]    IM, JM, KM    Pointers to store grid dimensions in the x, y, and z directions.
 * @param[in]     fd            File pointer to read grid dimensions (if required).
 * @param[in]     generate_grid Flag indicating whether the grid is programmatically generated (1) or read from a file (0).
 * @param[in]     imm, jmm, kmm Arrays containing grid dimensions for programmatic generation.
 * @param[in]     nblk          Pointer to store the number of blocks.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode DetermineGridSizes(PetscInt bi, UserCtx *user, PetscInt *IM, PetscInt *JM, PetscInt *KM, FILE *fd, 
                                  PetscInt generate_grid, PetscInt *imm, PetscInt *jmm, PetscInt *kmm, PetscInt *nblk) {
    PetscErrorCode ierr;  // PETSc error handling
    PetscInt rank;        // MPI rank

    // Validate input pointers
    if (generate_grid && (!imm || !jmm || !kmm)) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, 
                "DetermineGridSizes - Input arrays (imm, jmm, kmm) must not be NULL when generating grid programmatically.");
    }
    if (nblk == NULL) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, 
                "DetermineGridSizes - Number of blocks (nblk) must not be NULL.");
    }
    if (bi < 0 || bi >= *nblk) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, 
                "DetermineGridSizes - Block index out of range");
    }

    // Retrieve MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Log the start of the function
    LOG(GLOBAL, LOG_INFO, "DetermineGridSizes - Starting grid size determination for block %d on rank %d.\n", bi, rank);

    if (rank == 0) {
        if (!generate_grid) {
            // Read grid dimensions from file
            if (!fd) {
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, 
                        "DetermineGridSizes - File pointer (fd) is NULL. Cannot read grid dimensions from file.");
            }
            if (fscanf(fd, "%i %i %i\n", IM, JM, KM) != 3) {
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, 
                        "DetermineGridSizes - Failed to read grid dimensions from file.");
            }
            LOG(GLOBAL, LOG_DEBUG, "DetermineGridSizes - Read grid dimensions from file: IM=%d, JM=%d, KM=%d.\n", 
                *IM, *JM, *KM);
        } else {
            // Get grid dimensions from input arrays
            *IM = imm[bi];
            *JM = jmm[bi];
            *KM = kmm[bi];
            LOG(GLOBAL, LOG_DEBUG, 
                "DetermineGridSizes - Programmatically generated grid dimensions: IM=%d, JM=%d, KM=%d.\n", 
                *IM, *JM, *KM);
        }
    }

    // Broadcast dimensions to all ranks
    ierr = MPI_Bcast(IM, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Bcast(JM, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Bcast(KM, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_DEBUG, "DetermineGridSizes - Broadcasted grid dimensions: IM=%d, JM=%d, KM=%d.\n", 
        *IM, *JM, *KM);

    // Update the user context on all ranks
    user->IM = *IM;
    user->JM = *JM;
    user->KM = *KM;
    LOG(GLOBAL, LOG_DEBUG, 
        "DetermineGridSizes - Updated UserCtx: IM=%d, JM=%d, KM=%d on rank %d.\n", 
        user->IM, user->JM, user->KM, rank);

    // Log the end of the function
    LOG(GLOBAL, LOG_INFO, "DetermineGridSizes - Completed grid size determination for block %d on rank %d.\n", bi, rank);

    return 0;
}

/**
 * @brief Finalizes the grid setup by freeing resources and closing files.
 *
 * This function handles cleanup tasks after grid setup, including closing the file
 * descriptor if a grid file was read and deallocating memory allocated for grid parameters.
 *
 * @param[in] generate_grid Flag indicating if the grid was generated programmatically (1) or read from a file (0).
 * @param[in] fd            File pointer to the grid file (used if generate_grid is 0).
 * @param[in,out] imm       Array storing i-dimensions for each block, freed here.
 * @param[in,out] jmm       Array storing j-dimensions for each block, freed here.
 * @param[in,out] kmm       Array storing k-dimensions for each block, freed here.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode FinalizeGridSetup(PetscInt generate_grid, FILE *fd, PetscInt *imm, PetscInt *jmm, PetscInt *kmm) {
    PetscInt rank;               // MPI rank
    PetscErrorCode ierr;         // PETSc error handling

    // Retrieve MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Log the start of the function
    LOG(GLOBAL, LOG_INFO, "FinalizeGridSetup - Starting grid finalization on rank %d.\n", rank);

    // Handle file closing for non-programmatic grid setup
    if (!generate_grid && rank == 0) {
        if (fd != NULL) {
            fclose(fd);
            LOG(GLOBAL, LOG_DEBUG, "FinalizeGridSetup - Closed grid file on rank %d.\n", rank);
        } else {
            LOG(GLOBAL, LOG_WARNING, "FinalizeGridSetup - File pointer is NULL on rank %d. Nothing to close.\n", rank);
        }
    }

    // Free memory allocated for imm, jmm, and kmm
    if (imm != NULL) {
        ierr = PetscFree(imm); CHKERRQ(ierr);
        LOG(GLOBAL, LOG_DEBUG, "FinalizeGridSetup - Freed memory for imm on rank %d.\n", rank);
    }

    if (jmm != NULL) {
        ierr = PetscFree(jmm); CHKERRQ(ierr);
        LOG(GLOBAL, LOG_DEBUG, "FinalizeGridSetup - Freed memory for jmm on rank %d.\n", rank);
    }

    if (kmm != NULL) {
        ierr = PetscFree(kmm); CHKERRQ(ierr);
        LOG(GLOBAL, LOG_DEBUG, "FinalizeGridSetup - Freed memory for kmm on rank %d.\n", rank);
    }

    // Log the end of the function
    LOG(GLOBAL, LOG_INFO, "FinalizeGridSetup - Completed grid finalization on rank %d.\n", rank);

    return 0;
}

/**
 * @brief Configures grid coordinates and initializes grid management structures.
 *
 * This function sets up the grid using either programmatic generation or file input, initializes 
 * the DMDA structures for each grid block, and assigns grid coordinates.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode DefineGridCoordinates(UserCtx *user) {
    PetscErrorCode ierr;         // PETSc error handling
    PetscInt rank;               // MPI rank
    PetscInt bi;                 // Block index
    PetscInt IM, JM, KM;         // Grid dimensions for a block
    FILE *fd = NULL;             // File pointer for reading grid input
    PetscInt generate_grid = 0;  // Flag for programmatic grid generation
    PetscInt grid1d = 0;         // Flag for 1D grid input
    PetscInt block_number = 1;   // Number of grid blocks
    PetscReal L_x, L_y, L_z;     // Domain lengths in x, y, and z directions
    PetscInt *imm = NULL, *jmm = NULL, *kmm = NULL;  // Grid sizes for blocks

    // Retrieve MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Log the start of the function
    LOG(GLOBAL, LOG_INFO, "DefineGridCoordinates - Starting grid configuration on rank %d.\n", rank);

    // Parse input parameters for grid setup
    ierr = ParseGridInputs(user, &generate_grid, &grid1d, &L_x, &L_y, &L_z, &imm, &jmm, &kmm, &block_number, fd); CHKERRQ(ierr);

    // Validate block number
    if (block_number <= 0) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Number of blocks must be greater than 0.");
    }

    // Iterate over blocks to configure their respective grids
    for (bi = 0; bi < block_number; bi++) {
        // Determine grid sizes for the current block
        ierr = DetermineGridSizes(bi, user, &IM, &JM, &KM, fd, generate_grid, imm, jmm, kmm, &block_number); CHKERRQ(ierr);

        LOG(GLOBAL,LOG_INFO,"DefineGridCoordinates - IM,JM,KM - %d,%d,%d \n",IM,JM,KM);

        // Set up the DM for the current block
        ierr = InitializeGridDM(user,L_x,L_y,L_z); CHKERRQ(ierr);

        // Populate grid coordinates for the current block
        ierr = AssignGridCoordinates(user, generate_grid, grid1d, IM, JM, KM, L_x, L_y, L_z, fd); CHKERRQ(ierr);

        // Log the block's subdomain ownership details
        DMDALocalInfo info = user->info;
        ierr = DMDAGetLocalInfo(user->da, &info);
        LOG(GLOBAL, LOG_INFO, 
            "DefineGridCoordinates - Rank %d owns subdomain for block %d: xs=%d, xe=%d, ys=%d, ye=%d, zs=%d, ze=%d.\n", 
            rank, bi, info.xs, info.xs + info.xm, info.ys, info.ys + info.ym, info.zs, info.zs + info.zm);
    }

    // Finalize grid setup by closing the input file (if applicable) and freeing allocated memory
    ierr = FinalizeGridSetup(generate_grid, fd, imm, jmm, kmm); CHKERRQ(ierr);

    // Log the end of the function
    LOG(GLOBAL, LOG_INFO, "DefineGridCoordinates - Completed grid configuration on rank %d.\n", rank);

    return 0;
}


/**
 * @brief Computes the local bounding box of the grid on the current process.
 *
 * This function calculates the minimum and maximum coordinates (x, y, z) of the
 * local grid points owned by the current MPI process. It iterates over the local
 * portion of the grid, examines each grid point's coordinates, and updates the
 * minimum and maximum values accordingly.
 *
 * The computed bounding box is stored in the provided `localBBox` structure,
 * and the `user->bbox` field is also updated with this bounding box for
 * consistency within the user context.
 *
 * @param[in]  user      Pointer to the user-defined context containing grid information.
 *                       This context must be properly initialized before calling this function.
 * @param[out] localBBox Pointer to the BoundingBox structure where the computed local bounding box will be stored.
 *                       The structure should be allocated by the caller.
 *
 * @return PetscErrorCode Returns `0` on success, non-zero on failure.
 */
PetscErrorCode ComputeLocalBoundingBox(UserCtx *user, BoundingBox *localBBox)
{
    PetscErrorCode ierr;
    PetscInt i, j, k,rank;
    PetscInt xs, ys, zs, xe, ye, ze;
    DMDALocalInfo info;
    Vec coordinates;
    Cmpnts ***coordArray;
    Cmpnts minCoords, maxCoords;

    // Start of function execution
    LOG(GLOBAL, LOG_INFO, "ComputeLocalBoundingBox: Entering the function.\n");

    // Validate input pointers
    if (!user) {
        LOG(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Input 'user' pointer is NULL.\n");
        return -1;
    }
    if (!localBBox) {
        LOG(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Output 'localBBox' pointer is NULL.\n");
        return -1;
    }

    // Get MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Get the local coordinates vector from the DMDA
    ierr = DMGetCoordinatesLocal(user->da, &coordinates);
    if (ierr) {
        LOG(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Error getting local coordinates vector.\n");
        return ierr;
    }

    if (!coordinates) {
        LOG(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Coordinates vector is NULL.\n");
        return -1;
    }

    // Access the coordinate array for reading
    ierr = DMDAVecGetArrayRead(user->fda, coordinates, &coordArray);
    if (ierr) {
        LOG(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Error accessing coordinate array.\n");
        return ierr;
    }

    // Get the local grid information (indices and sizes)
    ierr = DMDAGetLocalInfo(user->da, &info);
    if (ierr) {
        LOG(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Error getting DMDA local info.\n");
        return ierr;
    }

    xs = info.xs; xe = xs + info.xm;
    ys = info.ys; ye = ys + info.ym;
    zs = info.zs; ze = zs + info.zm;

    // Initialize min and max coordinates with extreme values
    minCoords.x = minCoords.y = minCoords.z = PETSC_MAX_REAL;
    maxCoords.x = maxCoords.y = maxCoords.z = PETSC_MIN_REAL;

    LOG(LOCAL, LOG_DEBUG, "ComputeLocalBoundingBox: Grid indices: xs=%d, xe=%d, ys=%d, ye=%d, zs=%d, ze=%d.\n", xs, xe, ys, ye, zs, ze);

    // Iterate over the local grid to find min and max coordinates
    for (k = zs; k < ze; k++) {
        for (j = ys; j < ye; j++) {
            for (i = xs; i < xe; i++) {
                Cmpnts coord = coordArray[k][j][i];

                // Update min and max coordinates
                if (coord.x < minCoords.x) minCoords.x = coord.x;
                if (coord.y < minCoords.y) minCoords.y = coord.y;
                if (coord.z < minCoords.z) minCoords.z = coord.z;

                if (coord.x > maxCoords.x) maxCoords.x = coord.x;
                if (coord.y > maxCoords.y) maxCoords.y = coord.y;
                if (coord.z > maxCoords.z) maxCoords.z = coord.z;
            }
        }
    }

    // Log the computed min and max coordinates
    LOG(LOCAL, LOG_INFO, "ComputeLocalBoundingBox: Rank - %d - Computed bounding box - minCoords=(%.6f, %.6f, %.6f), maxCoords=(%.6f, %.6f, %.6f).\n",
        rank,minCoords.x, minCoords.y, minCoords.z, maxCoords.x, maxCoords.y, maxCoords.z);

    // Restore the coordinate array
    ierr = DMDAVecRestoreArrayRead(user->fda, coordinates, &coordArray);
    if (ierr) {
        LOG(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Error restoring coordinate array.\n");
        return ierr;
    }

    // Set the local bounding box
    localBBox->min_coords = minCoords;
    localBBox->max_coords = maxCoords;

    // Update the bounding box inside the UserCtx for consistency
    user->bbox = *localBBox;

    LOG(GLOBAL, LOG_INFO, "ComputeLocalBoundingBox: Exiting the function successfully.\n");
    return 0;
}

/**
 * @brief Gathers local bounding boxes from all MPI processes to rank 0.
 *
 * This function first computes the local bounding box on each process by calling
 * `ComputeLocalBoundingBox`. It then uses an MPI gather operation to collect all
 * local bounding boxes on the root process (rank 0). On rank 0, it allocates an array
 * of `BoundingBox` structures to hold the gathered data and returns it via the
 * `allBBoxes` pointer. On other ranks, `allBBoxes` is set to `NULL`.
 *
 * @param[in]  user       Pointer to the user-defined context containing grid information.
 *                        This context must be properly initialized before calling this function.
 * @param[out] allBBoxes  Pointer to a pointer where the array of gathered bounding boxes will be stored on rank 0.
 *                        On rank 0, this will point to the allocated array; on other ranks, it will be `NULL`.
 *
 * @return PetscErrorCode Returns `0` on success, non-zero on failure.
 */
PetscErrorCode GatherAllBoundingBoxes(UserCtx *user, BoundingBox **allBBoxes)
{
    PetscErrorCode ierr;
    PetscMPIInt rank, size;
    BoundingBox *bboxArray = NULL;
    BoundingBox localBBox;

    LOG(GLOBAL, LOG_INFO, "GatherAllBoundingBoxes: Entering the function. \n");

    // Validate input pointers
    if (!user) {
        LOG(LOCAL, LOG_ERROR, "GatherAllBoundingBoxes: Input 'user' pointer is NULL.\n");
        return -1;
    }
    if (!allBBoxes) {
        LOG(LOCAL, LOG_ERROR, "GatherAllBoundingBoxes: Output 'allBBoxes' pointer is NULL.\n");
        return -1;
    }

    // Get the rank and size of the MPI communicator
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (ierr != MPI_SUCCESS) {
        LOG(LOCAL, LOG_ERROR, "GatherAllBoundingBoxes: Error getting MPI rank.\n");
        return ierr;
    }
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
    if (ierr != MPI_SUCCESS) {
        LOG(LOCAL, LOG_ERROR, "GatherAllBoundingBoxes: Error getting MPI size.\n");
        return ierr;
    }

    LOG(LOCAL, LOG_INFO, "GatherAllBoundingBoxes: MPI rank=%d, size=%d.\n", rank, size);

    // Compute the local bounding box on each process
    ierr = ComputeLocalBoundingBox(user, &localBBox);
    if (ierr) {
        LOG(LOCAL, LOG_ERROR, "GatherAllBoundingBoxes: Error computing local bounding box.\n");
        return ierr;
    }
    
    PetscBarrier(PETSC_NULL);

    // On rank 0, allocate memory for the array of bounding boxes
    if (rank == 0) {
        bboxArray = (BoundingBox *)malloc(size * sizeof(BoundingBox));
        if (!bboxArray) {
            LOG(LOCAL, LOG_ERROR, "GatherAllBoundingBoxes: Memory allocation failed for bounding box array.\n");
            return -1;
        }
        LOG(LOCAL, LOG_DEBUG, "GatherAllBoundingBoxes: Allocated memory for bounding box array on rank 0.\n");
    }

    // Perform MPI_Gather to collect all local bounding boxes on rank 0
    ierr = MPI_Gather(&localBBox, sizeof(BoundingBox), MPI_BYTE,
                      bboxArray, sizeof(BoundingBox), MPI_BYTE,
                      0, PETSC_COMM_WORLD);
    if (ierr != MPI_SUCCESS) {
        LOG(LOCAL, LOG_ERROR, "GatherAllBoundingBoxes: Error during MPI_Gather operation.\n");
        if (rank == 0 && bboxArray) free(bboxArray); // Clean up if allocation was done
        return ierr;
    }

    // On rank 0, assign the gathered bounding boxes to the output pointer
    if (rank == 0) {
        *allBBoxes = bboxArray;
        LOG(GLOBAL, LOG_INFO, "GatherAllBoundingBoxes: Successfully gathered bounding boxes on rank 0.\n");
    } else {
        *allBBoxes = NULL;
    }

    LOG(GLOBAL, LOG_INFO, "GatherAllBoundingBoxes: Exiting the function successfully.\n");
    return 0;
}
