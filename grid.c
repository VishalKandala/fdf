
// grid.c 

#include "grid.h"
#include "logging.h"

extern PetscInt block_number;

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
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 *
 * @note
 * - The global dimensions of the grid are read from the `user` structure.
 * - The function also sets up uniform coordinates in the range [0, 1] along all axes.
 */
PetscErrorCode InitializeGridDM(UserCtx *user) {
    PetscErrorCode ierr; // PETSc error handling

    // Log the start of grid initialization
    LOG(GLOBAL, LOG_INFO, "InitializeGridDM - Starting DMDA grid initialization.\n");

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

    // Assign uniform coordinates [0, 1] along all axes
    ierr = DMDASetUniformCoordinates(user->da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_DEBUG, "InitializeGridDM - Uniform coordinates set in the range [0, 1].\n");

    // Retrieve the coordinate DM
    ierr = DMGetCoordinateDM(user->da, &(user->fda)); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_DEBUG, "InitializeGridDM - Coordinate DM retrieved.\n");

    // Log detailed information about DMDA and its coordinate DM
    if(get_log_level()==LOG_DEBUG){
      
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
 * @brief Retrieves grid-related inputs for setting up the simulation.
 *
 * This function reads grid parameters from the PETSc options database or, if specified,
 * from a file. The inputs include global dimensions, domain lengths, and block structure
 * for the grid. It also handles broadcasting relevant parameters across MPI ranks.
 *
 * @param[in,out] user          Pointer to the UserCtx structure containing grid details.
 * @param[out]    generate_grid Flag indicating whether the grid should be generated (1) or read from file (0).
 * @param[out]    grid1d        Flag indicating if the grid is 1D (1) or 3D (0).
 * @param[out]    L_x           Pointer to store the domain length in the x-direction.
 * @param[out]    L_y           Pointer to store the domain length in the y-direction.
 * @param[out]    L_z           Pointer to store the domain length in the z-direction.
 * @param[out]    imm           Array to store the i-dimensions for each block.
 * @param[out]    jmm           Array to store the j-dimensions for each block.
 * @param[out]    kmm           Array to store the k-dimensions for each block.
 * @param[out]    nblk          Pointer to store the number of blocks.
 * @param[out]    fd            File pointer for reading grid data (if needed).
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 *
 * @note
 * - This function requires `control.dat` to be present for reading PETSc options.
 * - The number of blocks (`nblk`) and dimensions (`imm`, `jmm`, `kmm`) are broadcast
 *   across all MPI ranks to ensure consistency.
 */
PetscErrorCode ParseGridInputs(UserCtx *user, PetscInt *generate_grid, PetscInt *grid1d, PetscReal *L_x, 
                                  PetscReal *L_y, PetscReal *L_z, PetscInt *imm, PetscInt *jmm, PetscInt *kmm, 
                                  PetscInt *nblk, FILE *fd) {
    PetscErrorCode ierr;  // PETSc error handling
    PetscInt rank;        // MPI rank
    PetscReal cl = 1.0;   // Characteristic length scale (default)

    // Log the start of grid input retrieval
    LOG(GLOBAL, LOG_INFO, "ParseGridInputs - Starting retrieval of grid inputs.\n");

    // Insert options from the "control.dat" file into the PETSc options database
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, "control.dat", PETSC_TRUE); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_DEBUG, "ParseGridInputs - Loaded options from control.dat.\n");

    // Retrieve MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Read grid options from the PETSc options database
    ierr = PetscOptionsGetInt(NULL, NULL, "-grid", generate_grid, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-grid1d", grid1d, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-chact_leng", &cl, NULL); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_DEBUG, "ParseGridInputs - Grid options retrieved: generate_grid=%d, grid1d=%d, cl=%.2f.\n", 
        *generate_grid, *grid1d, cl);

    // If generating a grid, retrieve the domain dimensions
    if (*generate_grid) {
        ierr = PetscOptionsGetReal(NULL, NULL, "-L_x", L_x, NULL); CHKERRQ(ierr);
        ierr = PetscOptionsGetReal(NULL, NULL, "-L_y", L_y, NULL); CHKERRQ(ierr);
        ierr = PetscOptionsGetReal(NULL, NULL, "-L_z", L_z, NULL); CHKERRQ(ierr);
        LOG(GLOBAL, LOG_DEBUG, "ParseGridInputs - Domain dimensions: L_x=%.2f, L_y=%.2f, L_z=%.2f.\n", 
            *L_x, *L_y, *L_z);

        ierr = PetscOptionsGetIntArray(NULL, NULL, "-im", imm, nblk, NULL); CHKERRQ(ierr);
        ierr = PetscOptionsGetIntArray(NULL, NULL, "-jm", jmm, nblk, NULL); CHKERRQ(ierr);
        ierr = PetscOptionsGetIntArray(NULL, NULL, "-km", kmm, nblk, NULL); CHKERRQ(ierr);
        LOG(GLOBAL, LOG_DEBUG, "ParseGridInputs - Block dimensions retrieved: imm[0]=%d, jmm[0]=%d, kmm[0]=%d.\n", 
            imm[0], jmm[0], kmm[0]);
    } else {
        // Open the grid data file on rank 0 and broadcast the number of blocks
        if (rank == 0) {
            fd = fopen("grid.dat", "r");
            fscanf(fd, "%i\n", nblk);
            LOG(GLOBAL, LOG_DEBUG, "ParseGridInputs - Opened grid.dat. Number of blocks: %d.\n", *nblk);
        }
        ierr = MPI_Bcast(nblk, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        LOG(GLOBAL, LOG_DEBUG, "ParseGridInputs - Broadcast number of blocks: %d.\n", *nblk);
    }

    // Log the successful completion of grid input retrieval
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
    PetscReal xx, *gc;            // Temporary storage for coordinates
    PetscReal cl = 1.0, L_dim = 1.0; // Normalization constants

    // Retrieve MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Log the start of grid population
    LOG(GLOBAL, LOG_INFO, "AssignGridCoordinates - Starting grid population for rank %d.\n", rank);

    // Retrieve local grid information
    DMDALocalInfo info = user->info;
    PetscInt xs = info.xs, xe = info.xs + info.xm; // Local x start and end indices
    PetscInt ys = info.ys, ye = info.ys + info.ym; // Local y start and end indices
    PetscInt zs = info.zs, ze = info.zs + info.zm; // Local z start and end indices

    // Allocate memory for grid coordinates based on 1D or 3D mode
    if (grid1d) {
        ierr = PetscMalloc1(IM + JM + KM, &gc); CHKERRQ(ierr);
    } else {
        ierr = PetscMalloc1(3 * IM * JM * KM, &gc); CHKERRQ(ierr);
    }
    LOG(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Memory allocated for grid coordinates.\n");

    // Access the DMDA coordinate vector
    ierr = DMGetCoordinatesLocal(user->da, &Coor); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, Coor, &coor); CHKERRQ(ierr);

    // Rank 0: Read or generate grid coordinates
    if (rank == 0) {
      if (grid1d) {
            // Read or generate 1D grid coordinates
            for (i = 0; i < IM; i++) fscanf(fd, "%le %le %le\n", &gc[i], &xx, &xx);
            for (j = 0; j < JM; j++) fscanf(fd, "%le %le %le\n", &xx, &gc[IM + j], &xx);
            for (k = 0; k < KM; k++) fscanf(fd, "%le %le %le\n", &xx, &xx, &gc[IM + JM + k]);

            LOG(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Rank 0 read 1D grid coordinates from file.\n");

            // Broadcast coordinates to all ranks
            ierr = MPI_Bcast(gc, IM + JM + KM, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
            LOG(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Broadcasted 1D grid coordinates.\n");
      } else {
            // Read or generate 3D grid coordinates
            for (k = 0; k < KM; k++) {
                for (j = 0; j < JM; j++) {
                    for (i = 0; i < IM; i++) {
                        if (!generate_grid) fscanf(fd, "%le", gc + (k * JM * IM + j * IM + i) * 3);
                        else *(gc + (k * JM * IM + j * IM + i) * 3) = L_x / IM * i;
                    }
                }
            }

            LOG(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Rank 0 read/generated 3D grid coordinates.\n");

            // Broadcast coordinates to all ranks
            ierr = MPI_Bcast(gc, 3 * IM * JM * KM, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
            LOG(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Broadcasted 3D grid coordinates.\n");
        }
    } // rank 0.


    // Broadcast call from all other processors to recieve from rank 0
    if (grid1d) {
    ierr = MPI_Bcast(gc, IM + JM + KM, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    } else {
    ierr = MPI_Bcast(gc, 3 * IM * JM * KM, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    }


    // Assign coordinates to the local grid
    if (grid1d) {
        for (k = zs; k < ze; k++) {
            for (j = ys; j < ye; j++) {
                for (i = xs; i < xe; i++) {
                    if (k < KM && j < JM && i < IM) {
                        coor[k][j][i].x = gc[i] / cl * L_dim;
                        coor[k][j][i].y = gc[IM + j] / cl * L_dim;
                        coor[k][j][i].z = gc[IM + JM + k] / cl * L_dim;
                    }
                }
            }
        }
    } else {
        for (k = zs; k < ze; k++) {
            for (j = ys; j < ye; j++) {
                for (i = xs; i < xe; i++) {
                    if (k < KM && j < JM && i < IM) {
                        coor[k][j][i].x = gc[(k * JM * IM + j * IM + i) * 3] / cl;
                        coor[k][j][i].y = gc[(k * JM * IM + j * IM + i) * 3 + 1] / cl;
                        coor[k][j][i].z = gc[(k * JM * IM + j * IM + i) * 3 + 2] / cl;
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
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 *
 * @note
 * - This function performs MPI communication to synchronize grid dimensions across all processes.
 */
PetscErrorCode DetermineGridSizes(PetscInt bi, UserCtx *user, PetscInt *IM, PetscInt *JM, PetscInt *KM, FILE *fd, 
                            PetscInt generate_grid, PetscInt *imm, PetscInt *jmm, PetscInt *kmm) {
    PetscErrorCode ierr;         // PETSc error handling
    PetscInt rank;               // MPI rank

    // Retrieve MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Log the start of the function
    LOG(GLOBAL, LOG_INFO, "DetermineGridSizes - Starting grid size determination for block %d on rank %d.\n", bi, rank);

    // Rank 0: Read dimensions from file or input arrays
    if (rank == 0) {
        if (!generate_grid) {
            // Read grid dimensions from file
            fscanf(fd, "%i %i %i\n", &(user->IM), &(user->JM), &(user->KM));
            LOG(GLOBAL, LOG_DEBUG, "DetermineGridSizes - Read grid dimensions from file: IM=%d, JM=%d, KM=%d.\n", 
                user->IM, user->JM, user->KM);
        } else {
            // Get grid dimensions from input arrays
            user->IM = imm[bi];
            user->JM = jmm[bi];
            user->KM = kmm[bi];
            LOG(GLOBAL, LOG_DEBUG, "DetermineGridSizes - Programmatically generated grid dimensions: IM=%d, JM=%d, KM=%d.\n", 
                user->IM, user->JM, user->KM);
        }

        // Update local pointers
        *IM = user->IM;
        *JM = user->JM;
        *KM = user->KM;

        // Broadcast dimensions to all ranks
        ierr = MPI_Bcast(&(user->IM), 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Bcast(&(user->JM), 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Bcast(&(user->KM), 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        LOG(GLOBAL, LOG_DEBUG, "DetermineGridSizes - Broadcasted grid dimensions to all ranks.\n");
    } else {
        // Non-rank 0: Receive dimensions from rank 0
        ierr = MPI_Bcast(&(user->IM), 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Bcast(&(user->JM), 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Bcast(&(user->KM), 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

        // Update local pointers
        *IM = user->IM;
        *JM = user->JM;
        *KM = user->KM;

        LOG(GLOBAL, LOG_DEBUG, "DetermineGridSizes - Received grid dimensions: IM=%d, JM=%d, KM=%d.\n", 
            *IM, *JM, *KM);
    }

    // Log the end of the function
    LOG(GLOBAL, LOG_INFO, "DetermineGridSizes - Completed grid size determination for block %d on rank %d.\n", bi, rank);

    return 0;
}

/**
 * @brief Finalizes grid operations by closing the grid file (if required).
 *
 * This function ensures proper cleanup after grid setup by closing any opened file
 * associated with the grid generation or reading process.
 *
 * @param[in] generate_grid Flag indicating whether the grid was programmatically generated (1) or read from a file (0).
 * @param[in] fd            File pointer for the grid file. It is closed if `generate_grid` is 0.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 *
 * @note
 * - This function only closes the file on rank 0.
 */
PetscErrorCode FinalizeGridSetup(PetscInt generate_grid, FILE *fd) {
    PetscInt rank;               // MPI rank
    PetscErrorCode ierr;         // PETSc error handling

    // Retrieve MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Log the start of the function
    LOG(GLOBAL, LOG_INFO, "FinalizeGridSetup - Starting grid finalization on rank %d.\n", rank);

    // Rank 0: Close the file if not programmatically generated
    if (rank == 0 && !generate_grid) {
        if (fd != NULL) {
            fclose(fd);
            LOG(GLOBAL, LOG_DEBUG, "FinalizeGridSetup - Closed grid file on rank %d.\n", rank);
        } else {
            LOG(GLOBAL, LOG_WARNING, "FinalizeGridSetup - File pointer is NULL on rank %d. Nothing to close.\n", rank);
        }
    }

    // Log the end of the function
    LOG(GLOBAL, LOG_INFO, "FinalizeGridSetup - Completed grid finalization on rank %d.\n", rank);

    return 0;
}

/**
 * @brief Defines the grid coordinates for the computational domain.
 *
 * This function orchestrates the process of defining grid coordinates by parsing input parameters,
 * determining grid sizes, setting up the DM structures, and assigning grid coordinates.
 * It supports both programmatically generated grids and grids read from an input file.
 *
 * @param[in,out] user Pointer to an array of UserCtx structures, one for each block.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 *
 * @note
 * - The function supports both 1D and 3D grids.
 * - The grid file is closed after use if the grid is read from a file.
 * - This function handles multi-block configurations.
 */
PetscErrorCode DefineGridCoordinates(UserCtx *user) {
    PetscErrorCode ierr;         // PETSc error handling
    PetscInt rank;               // MPI rank
    PetscInt bi;                 // Block index
    PetscInt IM, JM, KM;         // Grid dimensions for a block
    FILE *fd = NULL;             // File pointer for reading grid input
    PetscReal cl = 1.0;          // Characteristic length for normalization
    PetscInt generate_grid = 0;  // Flag for programmatic grid generation
    PetscInt grid1d = 0;         // Flag for 1D grid input
    PetscInt nblk = block_number; // Number of blocks
    PetscReal L_x, L_y, L_z;     // Domain lengths in x, y, and z directions
    PetscInt imm[block_number], jmm[block_number], kmm[block_number]; // Grid sizes for blocks

    // Get MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Log the start of the function
    LOG(GLOBAL, LOG_INFO, "DefineGridCoordinates - Starting grid configuration on rank %d.\n", rank);

    // Parse input parameters for grid setup
    ierr = ParseGridInputs(user, &generate_grid, &grid1d, &L_x, &L_y, &L_z, imm, jmm, kmm, &nblk, fd); CHKERRQ(ierr);

    // Iterate over blocks to configure their respective grids
    for (bi = 0; bi < block_number; bi++) {
        // Determine grid sizes for the current block
        ierr = DetermineGridSizes(bi, user, &IM, &JM, &KM, fd, generate_grid, imm, jmm, kmm); CHKERRQ(ierr);

        // Set up the DM for the current block
        ierr = InitializeGridDM(&user[bi]); CHKERRQ(ierr);

        // Populate grid coordinates for the current block
        ierr = AssignGridCoordinates(&user[bi], generate_grid, grid1d, IM, JM, KM, L_x, L_y, L_z, fd); CHKERRQ(ierr);

        // Log the block's subdomain ownership details
        DMDALocalInfo info = user[bi].info;
        LOG(GLOBAL, LOG_INFO, 
            "DefineGridCoordinates - Rank %d owns subdomain for block %d: xs=%d, xe=%d, ys=%d, ye=%d, zs=%d, ze=%d.\n", 
            rank, bi, info.xs, info.xs + info.xm, info.ys, info.ys + info.ym, info.zs, info.zs + info.zm);
    }

    // Finalize grid setup by closing the input file (if applicable)
    ierr = FinalizeGridSetup(generate_grid, fd); CHKERRQ(ierr);

    // Log the end of the function
    LOG(GLOBAL, LOG_INFO, "DefineGridCoordinates - Completed grid configuration on rank %d.\n", rank);

    return 0;
}
