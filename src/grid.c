 
// grid.c 

#include "grid.h"

#define BBOX_TOLERANCE 1e-9

//extern PetscInt block_number;

// DM da;  // For DOF=1, cell-centered scalars (Pressure, Nvert)
// DM fda; // For DOF=3, node-based, used for Coor, Ucat, Ucont storage

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
 * @param[in] generate_grid Flag indicating if the grid is being generated (1) or read from file (0).
 *                         Used to control logging verbosity.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 *
 * @note
 * - The global dimensions of the grid are read from the `user` structure.
 * - The function also sets up uniform coordinates in the range [0,Lx],[0,Ly],[0,Lz].
 * - When generate_grid=0 (reading from file), certain log messages about temporary coordinates are suppressed.
 */
PetscErrorCode InitializeGridDM(UserCtx *user, PetscInt *generate_grid) {
    PetscErrorCode ierr; // PETSc error handling
    PetscInt M,N,P; // No.of grid points/vector length
    PetscInt stencil_width = 2; //
    PetscBool file_grid_mode = (generate_grid && *generate_grid == 0) ? PETSC_TRUE : PETSC_FALSE;
    
    // Validate grid dimensions
    if (user->IM <= 0 || user->JM <= 0 || user->KM <= 0) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                "InitializeGridDM - Grid dimensions (IM, JM, KM) must be positive.");
    }

    // Log the start of grid initialization (always shown)
    LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeGridDM - Starting DMDA grid initialization for IM=%d, JM=%d, KM=%d.\n",
        user->IM, user->JM, user->KM);

    M = user->IM+1;
    N = user->JM+1;
    P = user->KM+1;

    // --- Option 1: Standard Approach (da=DOF1, fda=Derived DOF3) ---
    if (!file_grid_mode) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "InitializeGridDM - Creating DMDA 'da' (DOF=1, s=%d).\n", stencil_width);
    }
    
    ierr = DMDACreate3d(PETSC_COMM_WORLD,
                        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,  // Boundary types
                        DMDA_STENCIL_BOX,                                      // Stencil type
                        M,N,P,                                                 // Global dimensions (Nodes)
                        PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,             // Dimensions per process
                        1,                                                     // Degrees of freedom (Scalar)
                        stencil_width,                                         // Stencil width (s=2)
                        NULL, NULL, NULL,                                      // Array partitioning
                        &(user->da)); CHKERRQ(ierr);

    // --- Explicitly set the process grid layout ---

    // --- Set the Processor Layout from Options ---
    ierr = SetDMDAProcLayout(user->da, user); CHKERRQ(ierr); // Call the new function
    
    ierr = DMSetUp(user->da); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)user->da, "Scalar_Cell_DM"); CHKERRQ(ierr);
    
    // Retrieve the coordinate DM (PETSc creates/associates this)
    // This fda will have DOF=3 and compatible stencil width s=2
    ierr = DMGetCoordinateDM(user->da, &(user->fda)); CHKERRQ(ierr);
    ierr = DMSetUp(user->fda); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)user->fda, "Vector_Coord_DM_Derived"); CHKERRQ(ierr);
    
    if (!file_grid_mode) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "InitializeGridDM - Coordinate DM 'fda' retrieved/created by PETSc.\n");
    }

    // --- [DEBUGGING STEP 1: Explicitly Set Coordinate DM] ---
    // Although DMGetCoordinateDM likely does this, adding an explicit call
    // can sometimes resolve subtle association issues.
    ierr = DMSetCoordinateDM(user->da, user->fda); CHKERRQ(ierr);
    
    if (!file_grid_mode) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "InitializeGridDM - EXPLICITLY associated 'fda' as coordinate DM for 'da'.\n");
    }
    
    // Assign coordinates based on domain size
    PetscReal x_min = user->xMin, x_max = user->xMax;
    PetscReal y_min = user->yMin, y_max = user->yMax;
    PetscReal z_min = user->zMin, z_max = user->zMax;

    // Set coordinates using the DM that holds the coordinate vector.
    DM dm_for_coords = (user->da) ? user->da : user->fda; // Use da if it exists, else fda
    
    if (!file_grid_mode) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "InitializeGridDM - Setting uniform coordinates using DM %s...\n", 
                 (user->da)?"da":"fda");
    }
    
    ierr = DMDASetUniformCoordinates(dm_for_coords, x_min, x_max, y_min, y_max, z_min, z_max); CHKERRQ(ierr);
    
    if (!file_grid_mode) {
        LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "InitializeGridDM - Coordinates set in the range [%f, %f] x [%f, %f] x [%f, %f].\n",
            x_min, x_max, y_min, y_max, z_min, z_max);
    } else {
        // For file grid mode, we'll add a more informative message
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "InitializeGridDM - Creating initial uniform grid as placeholder (will be replaced by file data).\n");
    }

    // Verify fda is correctly retrieved/set if we used Option 1
    if (user->da) {
         ierr = DMGetCoordinateDM(user->da, &(user->fda)); CHKERRQ(ierr); // Re-get to be sure
         if (!user->fda) {
             SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Failed to retrieve Coordinate DM 'fda' after setting coordinates.");
         }
         
         if (!file_grid_mode) {
             LOG_ALLOW(GLOBAL, LOG_DEBUG, "InitializeGridDM - Verified Coordinate DM 'fda' retrieval after setting coords.\n");
         }
    } else {
         if (!file_grid_mode) {
             LOG_ALLOW(GLOBAL, LOG_DEBUG, "InitializeGridDM - Running in simplified mode (only fda).\n");
         }
    }

    // Debugging logs for DMDA and its coordinate DM - only when not in file grid mode
    if (!file_grid_mode && get_log_level() == LOG_DEBUG && is_function_allowed(__func__)) {
        if (user->da) {
             LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeGridDM - Viewing DMDA ('da'):\n");
             ierr = DMView(user->da, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        } else {
             LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeGridDM - 'da' is NULL (Simplified Setup).\n");
        }
        if (user->fda) {
             LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeGridDM - Viewing Coordinate DM ('fda'):\n");
             ierr = DMView(user->fda, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        } else {
            // This should not happen if Option 1 ran correctly
             LOG_ALLOW(GLOBAL, LOG_ERROR, "InitializeGridDM - 'fda' is NULL!\n");
        }
    }

    // Log the successful completion of grid initialization
    if (!file_grid_mode) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeGridDM - DMDA grid initialization completed successfully.\n");
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeGridDM - DMDA grid structure initialized (coordinates will be updated from file).\n");
    }

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
 * @param[out]    grid1d        Flag indicating whether the grid is 1D(1) or not(0).
 * @param[out]    xMin          Pointer to store the lower bound in x-direction.
 * @param[out]    xMax          Pointer to store the upper bound in x-direction.
 * @param[out]    yMin          Pointer to store the lower bound in y-direction.
 * @param[out]    yMax          Pointer to store the upper bound in y-direction.
 * @param[out]    zMin          Pointer to store the lower bound in z-direction.
 * @param[out]    zMax          Pointer to store the upper bound in z-direction.
 * @param[out]    imm           Pointer to Array to store the i-dimensions for each block.
 * @param[out]    jmm           Pointer to Array to store the j-dimensions for each block.
 * @param[out]    kmm           Pointer to Array to store the k-dimensions for each block.
 * @param[out]    nblk          Pointer to store the number of blocks.
 * @param[out]    fd            File Pointer for reading grid data (if applicable).
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode ParseGridInputs(UserCtx *user, PetscInt *generate_grid, PetscInt *grid1d, PetscReal *xMin, PetscReal *xMax,
                               PetscReal *yMin, PetscReal *yMax, PetscReal *zMin, PetscReal *zMax,  PetscInt **imm, PetscInt **jmm, PetscInt **kmm,
                               PetscInt *nblk, FILE **fd) {
    PetscErrorCode ierr;

    LOG_ALLOW(GLOBAL, LOG_INFO, "ParseGridInputs - Starting retrieval of grid inputs.\n");

    // Insert options from the control file
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, "control.dat", PETSC_TRUE); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ParseGridInputs - Loaded options from control.dat.\n");

    // Read the basic grid settings
    ierr = PetscOptionsGetInt(NULL, NULL, "-grid", generate_grid, NULL); CHKERRQ(ierr);

    if (*generate_grid) {
        // Delegate the reading of grid generation inputs
      ierr = ReadGridGenerationInputs(user, grid1d, xMin, xMax, yMin, yMax, zMin, zMax, imm, jmm, kmm, nblk); CHKERRQ(ierr);
    } else {
      // We are reading grid data from a file, check if 'grid.dat' is accessible
      char grid_filename[PETSC_MAX_PATH_LEN] = "grid.dat";
      ierr = PetscOptionsGetString(NULL, NULL, "-grid_file", grid_filename, 
				   PETSC_MAX_PATH_LEN, NULL); CHKERRQ(ierr);
    
      *fd = fopen(grid_filename, "r");
      if (!*fd) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
                "ParseGridInputs - Could not open grid file '%s'...", grid_filename);
      }

      LOG_ALLOW(GLOBAL, LOG_INFO, "ParseGridInputs - Reading grid from file: %s\n", grid_filename);
    
      // Validate file format with a magic number or header check
      char header[256];
      if (fgets(header, sizeof(header), *fd) == NULL || 
	  strncmp(header, "FDFGRID", 7) != 0) {
        fclose(*fd);
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
                "ParseGridInputs - Invalid grid file format. Expected FDFGRID header.");
      }
      rewind(*fd); // Reset file position for normal reading
      
      // If file is accessible, delegate reading of grid file data, including grid1d
      ierr = ReadGridFile(grid_filename, nblk, imm, jmm, kmm, grid1d, PETSC_COMM_WORLD); CHKERRQ(ierr);
      // HACK needs to be fixed later in a better way!!
      *xMax = *yMax = *zMax = 1.0; 
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "ParseGridInputs - Grid inputs retrieved successfully.\n");
    return 0;
}

/**
 * @brief Computes a stretched coordinate along one dimension.
 *
 * This function computes a coordinate based on a geometric stretching ratio.
 * If the ratio (r) is 1.0, a uniform distribution is used:
 *     x(i) = L * (i/N)
 *
 * If r != 1.0, a geometric stretching is applied:
 *     x(i) = L * [ (r^(i/N) - 1 ) / (r - 1) ]
 *
 * Here:
 * - i   : The current index along the dimension.
 * - N   : The total number of divisions along that dimension.
 * - L   : The length of the domain along that dimension.
 * - r   : The stretching ratio. r > 1.0 stretches the grid in a geometric fashion
 *         increasing spacing away from the start, whereas 0 < r < 1.0 would
 *         cluster points near the start.
 *
 * @param[in] i Current index (0 <= i <= N).
 * @param[in] N Number of segments along the dimension.
 * @param[in] L Total length of the domain.
 * @param[in] r Stretching ratio.
 *
 * @return PetscReal The computed coordinate at index i.
 *
 * @note This function does not return a PetscErrorCode because it
 *       does not allocate memory or call PETSc routines that can fail.
 *       It is just a helper calculation function.
 */
static inline PetscReal ComputeStretchedCoord(PetscInt i, PetscInt N, PetscReal L, PetscReal r) {
    // Handle the case where N=0 to avoid division by zero (though it should never happen in normal usage)
    if (N == 0) return 0.0;

    PetscReal fraction = (PetscReal)i / (PetscReal)N;

    if (r == 1.0) {
        // No stretching (uniform distribution)
        return L * fraction;
    } else {
        // Geometric stretching
        // r^(fraction) grows (or decays) geometrically depending on r.
        // When i=0, r^(0)=1, so x(0)=0; when i=N, r^(1)=r, so x(N)=L.
        // This distributes poPetscInts between 0 and L non-linearly.
        PetscReal numerator = PetscPowReal(r, fraction) - 1.0;
        PetscReal denominator = r - 1.0;
        return L * (numerator / denominator);
    }
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
 * @param[in]     fd            File Pointer for reading grid data (if required).
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 *
 * @note
 * - The function assumes that the DMDA (`user->da`) has already been set up.
 * - Coordinates are normalized using the characteristic length scale (`cl`) and domain length scale (`L_dim`).
 */
PetscErrorCode AssignGridCoordinates(UserCtx *user, PetscInt generate_grid, PetscInt grid1d, PetscInt IM, PetscInt JM, 
                            PetscInt KM,  FILE *fd) {
    PetscErrorCode ierr;          // PETSc error handling
    PetscInt i, j, k;             // loop counters
    PetscMPIInt rank;             // MPI rank
    Vec Coor, gCoor;              // Local and global coordinate vectors
    Cmpnts ***coor;               // 3D array for node coordinates
    PetscReal *gc;            // Temporary storage for coordinates
    PetscReal cl = 1.0, L_dim = 1.0; // Normalization constants
    PetscReal L_x,L_y,L_z;

    // Length of domain in each dimension. 
    L_x = user->xMax - user->xMin;
    L_y = user->yMax - user->yMin;
    L_z = user->zMax - user->zMin;

   
    // Retrieve MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Log the start of grid population
    LOG_ALLOW(GLOBAL, LOG_INFO, "AssignGridCoordinates - Starting grid population for rank %d.\n", rank);

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
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Memory allocated for grid coordinates.\n");

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
	      // Skipping header
	      /* we are in AssignGridCoordinates(), rank==0, !generate_grid, grid1d==0 */
	      if (rank == 0 && !generate_grid && !grid1d) {
		/* skip the header line, nblk line + one IM‑JM‑KM line per block */
		PetscInt headerLines = user->nblk + 2;          /* e.g. 1 + 1 = 2 */
		char dummy[256];
		for (PetscInt s = 0; s < headerLines; ++s) {
		  if (!fgets(dummy, sizeof(dummy), fd))
		    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
			    "Unexpected EOF while skipping grid header");
		}
		LOG_ALLOW(GLOBAL,LOG_DEBUG,
			  "AssignGridCoordinates - skipped %d header lines, now at coordinates.\n",
			  headerLines);
	      }
                // Read coordinates from file
                for (k = 0; k <= KM; k++) {
                    for (j = 0; j <= JM; j++) {
                        for (i = 0; i <= IM; i++) {
			  PetscInt index = (((k * (JM + 1) * (IM + 1)) + (j * (IM + 1)) + i) * 3); 
                            if (fscanf(fd, "%le %le %le\n", &gc[index], &gc[index + 1], &gc[index + 2]) != 3) {
                                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Error reading grid coordinates from file");
                            }
                        }
                    }
                }
                LOG_ALLOW(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Rank 0 read 3D grid coordinates from file.\n");
            } else {
                // Generate coordinates programmatically
                for (k = 0; k <= KM; k++) {
                    for (j = 0; j <= JM; j++) {
                        for (i = 0; i <= IM; i++) {
                            PetscInt index = ((k * (JM + 1) * (IM + 1)) + (j * (IM + 1)) + i) * 3;
                            gc[index]     = user->xMin + ComputeStretchedCoord(i, IM, L_x, user->rx);
                            gc[index + 1] = user->yMin + ComputeStretchedCoord(j, JM, L_y, user->ry);
                            gc[index + 2] = user->zMin + ComputeStretchedCoord(k, KM, L_z, user->rz);
                        }
                    }
                }
                LOG_ALLOW(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Rank 0 generated 3D grid coordinates programmatically.\n");
            }
        }
    }

    // Broadcast grid coordinates to all ranks from rank 0
    if (grid1d) {
        ierr = MPI_Bcast(gc, IM + JM + KM, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Broadcasted 1D grid coordinates.\n");
    } else {
        PetscInt total_points = 3 * (IM + 1) * (JM + 1) * (KM + 1);
        ierr = MPI_Bcast(gc, total_points, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
        LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Broadcasted 3D grid coordinates.\n");
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

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Local grid coordinates assigned for rank %d.\n", rank);

    // Clean up and restore arrays
    ierr = PetscFree(gc); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, Coor, &coor); CHKERRQ(ierr);

    // Sync coordinates across ranks
    ierr = DMGetCoordinates(user->da, &gCoor); CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(user->fda, Coor, INSERT_VALUES, gCoor); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(user->fda, Coor, INSERT_VALUES, gCoor); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL,LOG_DEBUG," Coordinates Local to Global completed on rank %d. \n",rank);
    
    ierr = DMGlobalToLocalBegin(user->fda, gCoor, INSERT_VALUES, Coor); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->fda, gCoor, INSERT_VALUES, Coor); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "AssignGridCoordinates - Grid population completed for rank %d.\n", rank);

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
 * @param[in]     fd            File Pointer to read grid dimensions (if required).
 * @param[in]     generate_grid Flag indicating whether the grid is programmatically generated (1) or read from a file (0).
 * @param[in]     imm, jmm, kmm Arrays containing grid dimensions for programmatic generation.
 * @param[in]     nblk          Pointer to store the number of blocks.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode DetermineGridSizes(PetscInt bi, UserCtx *user, PetscInt *IM, PetscInt *JM, PetscInt *KM, FILE *fd, 
                                  PetscInt generate_grid, PetscInt *imm, PetscInt *jmm, PetscInt *kmm, PetscInt *nblk) {
    PetscErrorCode ierr;  // PETSc error handling
    PetscMPIInt rank;        // MPI rank

    // Validate input Pointers
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
    LOG_ALLOW(GLOBAL, LOG_INFO, "DetermineGridSizes - Starting grid size determination for block %d on rank %d.\n", bi, rank);

    if (rank == 0) {
            // Get grid dimensions from input arrays
            *IM = imm[bi];
            *JM = jmm[bi];
            *KM = kmm[bi];
            LOG_ALLOW(GLOBAL, LOG_DEBUG, 
                "grid dimensions(No.of cells) : IM=%d, JM=%d, KM=%d.\n", 
                *IM, *JM, *KM);
      // }
    }

    // Broadcast dimensions to all ranks
    ierr = MPI_Bcast(IM, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Bcast(JM, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Bcast(KM, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "DetermineGridSizes - Broadcasted grid dimensions: IM=%d, JM=%d, KM=%d.\n", 
        *IM, *JM, *KM);

    // Update the user context on all ranks
    user->IM = *IM;
    user->JM = *JM;
    user->KM = *KM;
    LOG_ALLOW(GLOBAL, LOG_DEBUG, 
        "DetermineGridSizes - Updated UserCtx: IM=%d, JM=%d, KM=%d on rank %d.\n", 
        user->IM, user->JM, user->KM, rank);

    // Log the end of the function
    LOG_ALLOW(GLOBAL, LOG_INFO, "DetermineGridSizes - Completed grid size determination for block %d on rank %d.\n", bi, rank);

    return 0;
}

/**
 * @brief Finalizes the grid setup by freeing resources and closing files.
 *
 * This function handles cleanup tasks after grid setup, including closing the file
 * descriptor if a grid file was read and deallocating memory allocated for grid parameters.
 *
 * @param[in] generate_grid Flag indicating if the grid was generated programmatically (1) or read from a file (0).
 * @param[in] fd            File Pointer to the grid file (used if generate_grid is 0).
 * @param[in,out] imm       Array storing i-dimensions for each block, freed here.
 * @param[in,out] jmm       Array storing j-dimensions for each block, freed here.
 * @param[in,out] kmm       Array storing k-dimensions for each block, freed here.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode FinalizeGridSetup(PetscInt generate_grid, FILE *fd, PetscInt *imm, PetscInt *jmm, PetscInt *kmm) {
    PetscMPIInt rank;               // MPI rank
    PetscErrorCode ierr;         // PETSc error handling

    // Retrieve MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Log the start of the function
    LOG_ALLOW(GLOBAL, LOG_INFO, "FinalizeGridSetup - Starting grid finalization on rank %d.\n", rank);

    // Handle file closing for non-programmatic grid setup
    if (!generate_grid && rank == 0) {
        if (fd != NULL) {
            fclose(fd);
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "FinalizeGridSetup - Closed grid file on rank %d.\n", rank);
        } else {
            LOG_ALLOW(GLOBAL, LOG_WARNING, "FinalizeGridSetup - File Pointer is NULL on rank %d. Nothing to close.\n", rank);
        }
    }

    // Free memory allocated for imm, jmm, and kmm
    if (imm != NULL) {
        ierr = PetscFree(imm); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "FinalizeGridSetup - Freed memory for imm on rank %d.\n", rank);
    }

    if (jmm != NULL) {
        ierr = PetscFree(jmm); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "FinalizeGridSetup - Freed memory for jmm on rank %d.\n", rank);
    }

    if (kmm != NULL) {
        ierr = PetscFree(kmm); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "FinalizeGridSetup - Freed memory for kmm on rank %d.\n", rank);
    }

    // Log the end of the function
    LOG_ALLOW(GLOBAL, LOG_INFO, "FinalizeGridSetup - Completed grid finalization on rank %d.\n", rank);

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
/*
PetscErrorCode DefineGridCoordinates(UserCtx *user) {
    PetscErrorCode ierr;         // PETSc error handling
    PetscMPIInt rank;               // MPI rank
    PetscInt bi;                 // Block index
    PetscInt IM, JM, KM;         // Grid dimensions for a block
    FILE *fd = NULL;             // File Pointer for reading grid input
    PetscInt generate_grid = 0;  // Flag for programmatic grid generation
    PetscInt grid1d = 0;         // Flag for 1D grid input
    PetscInt block_number = 1;   // Number of grid blocks
    PetscReal xMin,yMin,zMin;    // Lower bound of grid in all dimensions.
    PetscReal xMax,yMax,zMax;    // Upper bound of grid in all dimensions.
    PetscInt *imm = NULL, *jmm = NULL, *kmm = NULL;  // Grid sizes for blocks

    // Retrieve MPI rank 
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Log the start of the function
    LOG_ALLOW(GLOBAL, LOG_INFO, "DefineGridCoordinates - Starting grid configuration on rank %d.\n", rank);

    // Parse input parameters for grid setup
    ierr = ParseGridInputs(user, &generate_grid, &grid1d, &xMin, &xMax, &yMin, &yMax, &zMin, &zMax, &imm, &jmm, &kmm, &block_number, &fd); CHKERRQ(ierr);

    // Update grid boundaries & number of blocks in user.
    user->xMin = xMin;
    user->xMax = xMax;
    user->yMin = yMin;
    user->yMax = yMax;
    user->zMin = zMin;
    user->zMax = zMax;
    user->nblk  = block_number;
    
    // Validate block number
    if (block_number <= 0) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Number of blocks must be greater than 0.");
    }

    // Iterate over blocks to configure their respective grids
    for (bi = 0; bi < block_number; bi++) {
        // Determine grid sizes for the current block
        ierr = DetermineGridSizes(bi, user, &IM, &JM, &KM, fd, generate_grid, imm, jmm, kmm, &block_number); CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL,LOG_INFO,"DefineGridCoordinates - IM,JM,KM - %d,%d,%d \n",IM,JM,KM);

        // Set up the DM for the current block
        ierr = InitializeGridDM(user,&generate_grid   ); CHKERRQ(ierr);

        // Populate grid coordinates for the current block
        ierr = AssignGridCoordinates(user, generate_grid, grid1d, IM, JM, KM,fd); CHKERRQ(ierr);

        // Log the block's subdomain ownership details
        DMDALocalInfo info = user->info;
        ierr = DMDAGetLocalInfo(user->da, &info);
        LOG_ALLOW(GLOBAL, LOG_INFO, 
            "DefineGridCoordinates - Rank %d owns subdomain for block %d: xs=%d, xe=%d, ys=%d, ye=%d, zs=%d, ze=%d.\n", 
            rank, bi, info.xs, info.xs + info.xm, info.ys, info.ys + info.ym, info.zs, info.zs + info.zm);
    }

    // Finalize grid setup by closing the input file (if applicable) and freeing allocated memory
    ierr = FinalizeGridSetup(generate_grid, fd, imm, jmm, kmm); CHKERRQ(ierr);

    // Log the end of the function
    LOG_ALLOW(GLOBAL, LOG_INFO, "DefineGridCoordinates - Completed grid configuration on rank %d.\n", rank);

    return 0;
}
*/

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
    PetscMPIInt rank;               // MPI rank
    PetscInt bi;                 // Block index
    PetscInt IM, JM, KM;         // Grid dimensions for a block
    FILE *fd = NULL;             // File Pointer for reading grid input
    PetscInt generate_grid = 0;  // Flag for programmatic grid generation
    PetscInt grid1d = 0;         // Flag for 1D grid input
    PetscInt block_number = 1;   // Number of grid blocks
    PetscReal xMin,yMin,zMin;    // Lower bound of grid in all dimensions.
    PetscReal xMax,yMax,zMax;    // Upper bound of grid in all dimensions.
    PetscInt *imm = NULL, *jmm = NULL, *kmm = NULL;  // Grid sizes for blocks

    // Retrieve MPI rank 
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    
    // Log the start of the function
    LOG_ALLOW(GLOBAL, LOG_INFO, "DefineGridCoordinates - Starting grid configuration on rank %d.\n", rank);
    
    // Parse input parameters for grid setup
    ierr = ParseGridInputs(user, &generate_grid, &grid1d, &xMin, &xMax, &yMin, &yMax, &zMin, &zMax, &imm, &jmm, &kmm, &block_number, &fd); CHKERRQ(ierr);
    
    // Update grid boundaries & number of blocks in user.
    user->xMin = xMin;
    user->xMax = xMax;
    user->yMin = yMin;
    user->yMax = yMax;
    user->zMin = zMin;
    user->zMax = zMax;
    user->nblk  = block_number;
    
    // Validate block number
    if (block_number <= 0) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Number of blocks must be greater than 0.");
    }
    
    // Iterate over blocks to configure their respective grids
    for (bi = 0; bi < block_number; bi++) {
        // Determine grid sizes for the current block
        ierr = DetermineGridSizes(bi, user, &IM, &JM, &KM, fd, generate_grid, imm, jmm, kmm, &block_number); CHKERRQ(ierr);
        
        // Only log grid dimensions when actually generating a grid or in debug mode
        if (generate_grid || get_log_level() >= LOG_DEBUG) {
            LOG_ALLOW(GLOBAL, LOG_INFO, "DefineGridCoordinates - IM,JM,KM - %d,%d,%d \n", IM, JM, KM);
        }
        
        // Set up the DM for the current block - pass generate_grid to control logging
        ierr = InitializeGridDM(user, &generate_grid); CHKERRQ(ierr);
        
        // Populate grid coordinates for the current block
        ierr = AssignGridCoordinates(user, generate_grid, grid1d, IM, JM, KM, fd); CHKERRQ(ierr);
        
        // Log the block's subdomain ownership details - always useful regardless of grid source
        DMDALocalInfo info = user->info;
        ierr = DMDAGetLocalInfo(user->da, &info);
        LOG_ALLOW(GLOBAL, LOG_INFO, 
            "DefineGridCoordinates - Rank %d owns subdomain for block %d: xs=%d, xe=%d, ys=%d, ye=%d, zs=%d, ze=%d.\n", 
            rank, bi, info.xs, info.xs + info.xm, info.ys, info.ys + info.ym, info.zs, info.zs + info.zm);
    }
    
    // Finalize grid setup by closing the input file (if applicable) and freeing allocated memory
    ierr = FinalizeGridSetup(generate_grid, fd, imm, jmm, kmm); CHKERRQ(ierr);
    
    // Log the end of the function
    if (generate_grid) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "DefineGridCoordinates - Completed grid generation on rank %d.\n", rank);
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO, "DefineGridCoordinates - Completed grid loading from file on rank %d.\n", rank);
    }
    
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
    PetscInt i, j, k;
    PetscMPIInt rank;
    PetscInt xs, ys, zs, xe, ye, ze;
    DMDALocalInfo info;
    Vec coordinates;
    Cmpnts ***coordArray;
    Cmpnts minCoords, maxCoords;

    // Start of function execution
    LOG_ALLOW(GLOBAL, LOG_INFO, "ComputeLocalBoundingBox: Entering the function.\n");

    // Validate input Pointers
    if (!user) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Input 'user' Pointer is NULL.\n");
        return PETSC_ERR_ARG_NULL;
    }
    if (!localBBox) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Output 'localBBox' Pointer is NULL.\n");
        return PETSC_ERR_ARG_NULL;
    }

    // Get MPI rank
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Get the local coordinates vector from the DMDA
    ierr = DMGetCoordinatesLocal(user->da, &coordinates);
    if (ierr) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Error getting local coordinates vector.\n");
        return ierr;
    }

    if (!coordinates) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Coordinates vector is NULL.\n");
        return PETSC_ERR_ARG_NULL;
    }

    // Access the coordinate array for reading
    ierr = DMDAVecGetArrayRead(user->fda, coordinates, &coordArray);
    if (ierr) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Error accessing coordinate array.\n");
        return ierr;
    }

    // Get the local grid information (indices and sizes)
    ierr = DMDAGetLocalInfo(user->da, &info);
    if (ierr) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Error getting DMDA local info.\n");
        return ierr;
    }

    xs = info.xs; xe = xs + info.xm;
    ys = info.ys; ye = ys + info.ym;
    zs = info.zs; ze = zs + info.zm;

    // Initialize min and max coordinates with extreme values
    minCoords.x = minCoords.y = minCoords.z = PETSC_MAX_REAL;
    maxCoords.x = maxCoords.y = maxCoords.z = PETSC_MIN_REAL;

    LOG_ALLOW(LOCAL, LOG_DEBUG, "ComputeLocalBoundingBox: Grid indices: xs=%d, xe=%d, ys=%d, ye=%d, zs=%d, ze=%d.\n", xs, xe, ys, ye, zs, ze);

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


    // Add tolerance to bboxes.
    minCoords.x =  minCoords.x - BBOX_TOLERANCE;
    minCoords.y =  minCoords.y - BBOX_TOLERANCE;
    minCoords.z =  minCoords.z - BBOX_TOLERANCE;

    maxCoords.x =  maxCoords.x + BBOX_TOLERANCE;
    maxCoords.y =  maxCoords.y + BBOX_TOLERANCE;
    maxCoords.z =  maxCoords.z + BBOX_TOLERANCE;

    LOG_ALLOW(LOCAL,LOG_INFO," Tolerance added to the limits: %.8e .\n",(PetscReal)BBOX_TOLERANCE);
       
    // Log the computed min and max coordinates
     LOG_ALLOW(LOCAL, LOG_INFO,"Rank - %d - Computed bounding box - minCoords=(%.6f, %.6f, %.6f), maxCoords=(%.6f, %.6f, %.6f).\n",rank,minCoords.x, minCoords.y, minCoords.z, maxCoords.x, maxCoords.y, maxCoords.z);


    
    // Restore the coordinate array
    ierr = DMDAVecRestoreArrayRead(user->fda, coordinates, &coordArray);
    if (ierr) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Error restoring coordinate array.\n");
        return ierr;
    }

    // Set the local bounding box
    localBBox->min_coords = minCoords;
    localBBox->max_coords = maxCoords;

    // Update the bounding box inside the UserCtx for consistency
    user->bbox = *localBBox;

    LOG_ALLOW(GLOBAL, LOG_INFO, "ComputeLocalBoundingBox: Exiting the function successfully.\n");
    return 0;
}

/**
 * @brief Gathers local bounding boxes from all MPI processes to rank 0.
 *
 * This function first computes the local bounding box on each process by calling
 * `ComputeLocalBoundingBox`. It then uses an MPI gather operation to collect all
 * local bounding boxes on the root process (rank 0). On rank 0, it allocates an array
 * of `BoundingBox` structures to hold the gathered data and returns it via the
 * `allBBoxes` Pointer. On other ranks, `allBBoxes` is set to `NULL`.
 *
 * @param[in]  user       Pointer to the user-defined context containing grid information.
 *                        This context must be properly initialized before calling this function.
 * @param[out] allBBoxes  Pointer to a Pointer where the array of gathered bounding boxes will be stored on rank 0.
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

    LOG_ALLOW(GLOBAL, LOG_INFO, "Entering the function. \n");

    // Validate input Pointers
    if (!user) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Input 'user' Pointer is NULL.\n");
        return PETSC_ERR_ARG_NULL;
    }
    if (!allBBoxes) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Output 'allBBoxes' Pointer is NULL.\n");
        return PETSC_ERR_ARG_NULL;
    }

    // Get the rank and size of the MPI communicator
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (ierr != MPI_SUCCESS) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Error getting MPI rank.\n");
        return ierr;
    }
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
    if (ierr != MPI_SUCCESS) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Error getting MPI size.\n");
        return ierr;
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "MPI rank=%d, size=%d.\n", rank, size);

    // Compute the local bounding box on each process
    ierr = ComputeLocalBoundingBox(user, &localBBox);
    if (ierr) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Error computing local bounding box.\n");
        return ierr;
    }
    
    PetscBarrier(PETSC_NULLPTR);

    // On rank 0, allocate memory for the array of bounding boxes
    if (rank == 0) {
        bboxArray = (BoundingBox *)malloc(size * sizeof(BoundingBox));
        if (!bboxArray) {
            LOG_ALLOW(LOCAL, LOG_ERROR, "GatherAllBoundingBoxes: Memory allocation failed for bounding box array.\n");
            return PETSC_ERR_MEM;
        }
        LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAllBoundingBoxes: Allocated memory for bounding box array on rank 0.\n");
    }

    // Perform MPI_Gather to collect all local bounding boxes on rank 0
    // Corrected MPI_Gather call
ierr = MPI_Gather(&localBBox, sizeof(BoundingBox), MPI_BYTE,
                  (rank == 0) ? bboxArray : NULL,  // Explicitly NULL on non-roots
                  sizeof(BoundingBox), MPI_BYTE,   // Recv count is ignored on non-roots
                  0, PETSC_COMM_WORLD); CHKERRMPI(ierr);
 
//   ierr = MPI_Gather(&localBBox, sizeof(BoundingBox), MPI_BYTE,
// bboxArray, sizeof(BoundingBox), MPI_BYTE,
//                      0, PETSC_COMM_WORLD);
    if (ierr != MPI_SUCCESS) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "GatherAllBoundingBoxes: Error during MPI_Gather operation.\n");
        if (rank == 0 && bboxArray) free(bboxArray); // Clean up if allocation was done
        return ierr;
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "[Rank %d] Successfully gathered bounding boxes on rank 0.\n",rank);    
    
    // On rank 0, assign the gathered bounding boxes to the output Pointer
    if (rank == 0) {
        *allBBoxes = bboxArray;
    } else {
        *allBBoxes = NULL;
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Exiting the function successfully.\n");
    return 0;
}

/**
 * @brief Broadcasts the bounding box information collected on rank 0 to all other ranks.
 *
 * This function assumes that `GatherAllBoundingBoxes()` was previously called, so `bboxlist`
 * is allocated and populated on rank 0. All other ranks will allocate memory for `bboxlist`,
 * and this function will use MPI_Bcast to distribute the bounding box data to them.
 *
 * @param[in]     user      Pointer to the UserCtx structure. (Currently unused in this function, but kept for consistency.)
 * @param[in,out] bboxlist  Pointer to the array of BoundingBoxes. On rank 0, this should point to
 *                          a valid array of size 'size' (where size is the number of MPI ranks).
 *                          On non-root ranks, this function will allocate memory for `bboxlist`.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on MPI or PETSc-related errors.
 */
PetscErrorCode BroadcastAllBoundingBoxes(UserCtx *user, BoundingBox **bboxlist) {
    PetscErrorCode ierr;
    PetscMPIInt rank, size;

    // Get MPI rank and size
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    // On non-root ranks, allocate memory for bboxlist before receiving the broadcast
    if (rank != 0) {
        *bboxlist = (BoundingBox *)malloc(size * sizeof(BoundingBox));
        if (!*bboxlist) SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_MEM, "Failed to allocate memory for bboxlist on non-root ranks.");
    }

    LOG_ALLOW(LOCAL, LOG_INFO, "Broadcasting bounding box information from rank 0.\n");

    // Broadcast bboxlist from rank 0 to all other ranks
    ierr = MPI_Bcast(*bboxlist, (PetscInt)(size * sizeof(BoundingBox)), MPI_BYTE, 0, PETSC_COMM_WORLD);
    if (ierr != MPI_SUCCESS) {
        SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_LIB, "MPI_Bcast failed for bboxlist.");
    }

    LOG_ALLOW(LOCAL, LOG_INFO, "Broadcasted bounding box information from rank 0.\n");    

    return 0;
}
