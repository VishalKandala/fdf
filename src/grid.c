// ============================================================================
// grid.c - Always use da for coordinate data and fda for cell-centered solution data
//
//   user->fda : cell-centered DM with boundary ghost cells (IM+2 x JM+2 x KM+2).
//   user->da  : coordinate DM with physical node layout (IM+1 x JM+1 x KM+1).
//
// ============================================================================

#include "grid.h"

/**
 * @brief Initializes the DMDA grid for the simulation.
 *
 * This function sets up two 3D DMDAs:
 * 1) A "cell-centered" DM (user->fda) with ghost boundaries and an extra layer on each side.
 * 2) A separate "coordinate" DM (user->da) that holds the physical node coordinates
 *    with no extra boundary cells.
 *
 * @param[in,out] user   Pointer to the UserCtx structure containing grid details.
 * @param[in]     L_x,L_y,L_z  Physical dimensions of the domain.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 *
 * @note
 * - The cell-centered DM is created with size (IM+2, JM+2, KM+2) and DM_BOUNDARY_GHOSTED.
 * - The coordinate DM is created with size (IM+1, JM+1, KM+1) and DM_BOUNDARY_NONE.
 * - The coordinate DM is then attached to the cell-centered DM via DMSetCoordinateDM.
 */
PetscErrorCode InitializeGridDM(UserCtx *user, PetscReal L_x, PetscReal L_y, PetscReal L_z)
{
    PetscErrorCode ierr;
    PetscInt       Nx, Ny, Nz;
    PetscInt       NxCoord, NyCoord, NzCoord;

    if (L_x <= 0.0 || L_y <= 0.0 || L_z <= 0.0) {
        LOG_ALLOW(GLOBAL, LOG_ERROR,
            "InitializeGridDM - Invalid domain lengths: L_x=%.3f, L_y=%.3f, L_z=%.3f\n",
            (double)L_x, (double)L_y, (double)L_z);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
            "InitializeGridDM - Domain lengths must be positive in all directions.");
    }
    if (user->IM <= 0 || user->JM <= 0 || user->KM <= 0) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                "InitializeGridDM - Grid dimensions (IM, JM, KM) must be positive.");
    }

    LOG_ALLOW(GLOBAL, LOG_INFO,
              "InitializeGridDM - Initializing DM for IM=%d, JM=%d, KM=%d\n",
              user->IM, user->JM, user->KM);

    // -------------------------------------------------------------------------
    // 1) Create the cell-centered DM (fda) with ghost boundaries and +2 in each dimension.
    // -------------------------------------------------------------------------
    Nx = user->IM + 2;
    Ny = user->JM + 2;
    Nz = user->KM + 2;

    ierr = DMDACreate3d(
                PETSC_COMM_WORLD,
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                DMDA_STENCIL_BOX,
                Nx, Ny, Nz,
                PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                3, // 3 dof
                1, // stencil width
                NULL, NULL, NULL,
                &user->fda);
    CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "InitializeGridDM - Created cell-centered DM (fda) with ghost boundaries.\n");

    ierr = DMSetUp(user->fda); CHKERRQ(ierr);

    // -------------------------------------------------------------------------
    // 2) Create a separate DM (da) for coordinates, sized (IM+1, JM+1, KM+1).
    // -------------------------------------------------------------------------
    NxCoord = user->IM + 1;
    NyCoord = user->JM + 1;
    NzCoord = user->KM + 1;

    ierr = DMDACreate3d(
                PETSC_COMM_WORLD,
                DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
                DMDA_STENCIL_BOX,
                NxCoord, NyCoord, NzCoord,
                PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                3, // 3 dof for x,y,z
                0, // no extra stencil
                NULL, NULL, NULL,
                &user->da);
    CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "InitializeGridDM - Created coordinate DM (da) with ghost boundaries.\n");

    ierr = DMSetUp(user->da); CHKERRQ(ierr);

    // -------------------------------------------------------------------------
    // 3) Set uniform coordinates in [0,L_x], [0,L_y], [0,L_z] for the coordinate DM.
    // -------------------------------------------------------------------------
    ierr = DMDASetUniformCoordinates(user->da, 0.0, L_x, 0.0, L_y, 0.0, L_z);
    CHKERRQ(ierr);

    // 4) Attach the coordinate DM to the cell-centered DM
    ierr = DMSetCoordinateDM(user->fda, user->da); CHKERRQ(ierr);

    if (get_log_level() == LOG_DEBUG && is_function_allowed(__func__)) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeGridDM - Viewing fda (cell-centered DM):\n");
        ierr = DMView(user->fda, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeGridDM - Viewing da (coordinate DM):\n");
        ierr = DMView(user->da, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeGridDM - Completed DM setup.\n");
    return 0;
}

/**
 * @brief Parses grid-related inputs for setting up the simulation.
 *
 * This function coordinates the retrieval of grid parameters by delegating
 * to appropriate helper functions based on the grid source (generation or file).
 *
 * @param[in,out] user   Pointer to the UserCtx structure containing grid details.
 * @param[out]    generate_grid 1 if the grid should be generated programmatically, 0 if from file.
 * @param[out]    grid1d        1 if the grid is 1D, 0 if not.
 * @param[out]    L_x,L_y,L_z   Domain extents in x, y, z.
 * @param[out]    imm,jmm,kmm   Arrays holding per-block dimensions if generating programmatically.
 * @param[out]    nblk          Number of blocks.
 * @param[out]    fd            File pointer for reading grid data (if needed).
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode ParseGridInputs(UserCtx *user,
                               PetscInt *generate_grid, PetscInt *grid1d,
                               PetscReal *L_x, PetscReal *L_y, PetscReal *L_z,
                               PetscInt **imm, PetscInt **jmm, PetscInt **kmm,
                               PetscInt *nblk, FILE *fd)
{
    PetscErrorCode ierr;

    LOG_ALLOW(GLOBAL, LOG_INFO, "ParseGridInputs - Retrieving grid inputs.\n");

    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, "control.dat", PETSC_TRUE); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ParseGridInputs - Loaded options from control.dat.\n");

    ierr = PetscOptionsGetInt(NULL, NULL, "-grid", generate_grid, NULL); CHKERRQ(ierr);

    if (*generate_grid) {
        ierr = ReadGridGenerationInputs(user, grid1d, L_x, L_y, L_z, imm, jmm, kmm, nblk); CHKERRQ(ierr);
    } else {
        FILE *testFile = fopen("grid.dat", "r");
        if (!testFile) {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
                    "ParseGridInputs - Could not open 'grid.dat'.");
        } else {
            fclose(testFile);
        }
        ierr = ReadGridFile("grid.dat", nblk, imm, jmm, kmm, grid1d, PETSC_COMM_WORLD); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "ParseGridInputs - Finished reading grid inputs.\n");
    return 0;
}

/**
 * @brief Computes a stretched coordinate along one dimension.
 *
 * If r=1.0, it uses uniform spacing; otherwise a geometric stretch.
 */
static inline PetscReal ComputeStretchedCoord(PetscInt i, PetscInt N, PetscReal L, PetscReal r)
{
    if (N == 0) return 0.0;
    PetscReal fraction = (PetscReal)i / (PetscReal)N;

    if (r == 1.0) {
        return L * fraction;
    } else {
        PetscReal numerator   = PetscPowReal(r, fraction) - 1.0;
        PetscReal denominator = r - 1.0;
        return L * (numerator / denominator);
    }
}

/**
 * @brief Populates the coordinate DM (da) with coordinates from file or generates them.
 *
 * This function assigns coordinate values to the vector stored in user->da. By default,
 * user->da has (IM+1, JM+1, KM+1) global nodes in x,y,z. The user may have indicated
 * generation (stretching) or reading from a file.
 *
 * @param[in,out] user        Pointer to the UserCtx structure.
 * @param[in]     generate_grid 1 if generating programmatically, 0 if reading from file.
 * @param[in]     grid1d      1 if grid is 1D, 0 otherwise.
 * @param[in]     IM,JM,KM    Physical dimension counts for x,y,z cells (so # of nodes is IM+1, etc.).
 * @param[in]     L_x,L_y,L_z Domain lengths in x,y,z.
 * @param[in]     fd          File pointer if reading from file.
 */
PetscErrorCode AssignGridCoordinates(UserCtx *user,
                                     PetscInt generate_grid, PetscInt grid1d,
                                     PetscInt IM, PetscInt JM, PetscInt KM,
                                     PetscReal L_x, PetscReal L_y, PetscReal L_z,
                                     FILE *fd)
{
    PetscErrorCode ierr;
    PetscInt       rank, i, j, k;
    Vec            Coor;
    Cmpnts       ***coorArray;
    PetscReal     *gc;
    PetscReal      cl = 1.0, L_dim = 1.0;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "AssignGridCoordinates - Start on rank %d.\n", rank);

    // CHANGED: We now get the local info from the coordinate DM (user->da):
    DMDALocalInfo cinfo; 
    ierr = DMDAGetLocalInfo(user->da, &cinfo); CHKERRQ(ierr);
    PetscInt xs = cinfo.xs, xe = xs + cinfo.xm;
    PetscInt ys = cinfo.ys, ye = ys + cinfo.ym;
    PetscInt zs = cinfo.zs, ze = zs + cinfo.zm;

    // Allocate a buffer for reading or generating node coordinates
    if (grid1d) {
        ierr = PetscMalloc1(IM + JM + KM, &gc); CHKERRQ(ierr);
    } else {
        ierr = PetscMalloc1(3 * (IM + 1) * (JM + 1) * (KM + 1), &gc); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Allocated coordinate buffer.\n");

    // CHANGED: Retrieve the local coordinate vector from user->da directly
    ierr = DMGetCoordinatesLocal(user->da, &Coor); CHKERRQ(ierr);
    // CHANGED: Then get the array from user->da
    ierr = DMDAVecGetArray(user->da, Coor, &coorArray); CHKERRQ(ierr);

    // On rank 0, read/generate
    if (rank == 0) {
        if (grid1d) {
            // Possibly read or generate 1D coords...
        } else {
            if (!generate_grid) {
                // read from file
                for (k = 0; k <= KM; k++) {
                    for (j = 0; j <= JM; j++) {
                        for (i = 0; i <= IM; i++) {
                            PetscInt idx = ((k*(JM+1)*(IM+1)) + (j*(IM+1)) + i)*3;
                            if (fscanf(fd, "%le %le %le\n",
                                       &gc[idx], &gc[idx+1], &gc[idx+2]) != 3) {
                                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
                                        "Error reading grid coords");
                            }
                        }
                    }
                }
            } else {
                // generate coords
                for (k = 0; k <= KM; k++) {
                    for (j = 0; j <= JM; j++) {
                        for (i = 0; i <= IM; i++) {
                            PetscInt idx = ((k*(JM+1)*(IM+1)) + (j*(IM+1)) + i)*3;
                            gc[idx]   = ComputeStretchedCoord(i, IM, L_x, user->rx);
                            gc[idx+1] = ComputeStretchedCoord(j, JM, L_y, user->ry);
                            gc[idx+2] = ComputeStretchedCoord(k, KM, L_z, user->rz);
                        }
                    }
                }
            }
        }
    }

    // Broadcast from rank 0 to all
    if (grid1d) {
        ierr = MPI_Bcast(gc, IM + JM + KM, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    } else {
        PetscInt total_pts = 3 * (IM+1)*(JM+1)*(KM+1);
        ierr = MPI_Bcast(gc, total_pts, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    }

    // Fill the local portion of coorArray
    if (!grid1d) {
        for (k = zs; k < ze; k++) {
            for (j = ys; j < ye; j++) {
                for (i = xs; i < xe; i++) {
                    if (k <= KM && j <= JM && i <= IM) {
                        PetscInt idx = ((k*(JM+1)*(IM+1)) + (j*(IM+1)) + i)*3;
                        coorArray[k][j][i].x = gc[idx]     / cl * L_dim;
                        coorArray[k][j][i].y = gc[idx + 1] / cl * L_dim;
                        coorArray[k][j][i].z = gc[idx + 2] / cl * L_dim;
                    }
                }
            }
        }
    }

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Assigned local coords rank %d.\n", rank);

    ierr = PetscFree(gc); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, Coor, &coorArray); CHKERRQ(ierr);

    // Sync local->global
    Vec globalCoor, localCoor;
    ierr = DMGetCoordinates(user->da, &globalCoor); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(user->da, &localCoor); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(user->da, globalCoor, INSERT_VALUES, localCoor); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->da, globalCoor, INSERT_VALUES, localCoor); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "AssignGridCoordinates - Done on rank %d.\n", rank);
    return 0;
}

/**
 * @brief Determines the grid dimensions for the specified block.
 *
 * (Unchanged logic)
 */
PetscErrorCode DetermineGridSizes(PetscInt bi, UserCtx *user,
                                  PetscInt *IM, PetscInt *JM, PetscInt *KM,
                                  FILE *fd,
                                  PetscInt generate_grid,
                                  PetscInt *imm, PetscInt *jmm, PetscInt *kmm,
                                  PetscInt *nblk)
{
    PetscErrorCode ierr;
    PetscInt       rank;

    if (generate_grid && (!imm || !jmm || !kmm)) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
                "DetermineGridSizes - (imm, jmm, kmm) must not be NULL when generating.");
    }
    if (!nblk) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
                "DetermineGridSizes - nblk must not be NULL.");
    }
    if (bi < 0 || bi >= *nblk) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                "DetermineGridSizes - Block index out of range");
    }

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO,
        "DetermineGridSizes - Block %d on rank %d.\n", bi, rank);

    if (rank == 0) {
        if (!generate_grid) {
            if (!fd) {
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
                        "DetermineGridSizes - No file pointer for reading dims.");
            }
            if (fscanf(fd, "%i %i %i\n", IM, JM, KM) != 3) {
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
                        "DetermineGridSizes - Failed reading grid dims.");
            }
            LOG_ALLOW(GLOBAL, LOG_DEBUG,
                "DetermineGridSizes - Read IM=%d, JM=%d, KM=%d from file.\n",
                *IM, *JM, *KM);
        } else {
            *IM = imm[bi];
            *JM = jmm[bi];
            *KM = kmm[bi];
            LOG_ALLOW(GLOBAL, LOG_DEBUG,
                "DetermineGridSizes - Programmatic IM=%d, JM=%d, KM=%d.\n",
                *IM, *JM, *KM);
        }
    }

    ierr = MPI_Bcast(IM, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Bcast(JM, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Bcast(KM, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

    user->IM = *IM;
    user->JM = *JM;
    user->KM = *KM;

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
        "DetermineGridSizes - user->IM=%d, JM=%d, KM=%d on rank %d.\n",
        user->IM, user->JM, user->KM, rank);

    return 0;
}

/**
 * @brief Finalizes the grid setup by freeing resources and closing files.
 *
 * (Unchanged)
 */
PetscErrorCode FinalizeGridSetup(PetscInt generate_grid, FILE *fd,
                                 PetscInt *imm, PetscInt *jmm, PetscInt *kmm)
{
    PetscErrorCode ierr;
    PetscInt       rank;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "FinalizeGridSetup - rank %d.\n", rank);

    if (!generate_grid && rank == 0) {
        if (fd != NULL) {
            fclose(fd);
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "FinalizeGridSetup - Closed grid file rank %d.\n", rank);
        }
    }

    if (imm != NULL) {
        ierr = PetscFree(imm); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "FinalizeGridSetup - Freed imm rank %d.\n", rank);
    }
    if (jmm != NULL) {
        ierr = PetscFree(jmm); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "FinalizeGridSetup - Freed jmm rank %d.\n", rank);
    }
    if (kmm != NULL) {
        ierr = PetscFree(kmm); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "FinalizeGridSetup - Freed kmm rank %d.\n", rank);
    }

    return 0;
}

/**
 * @brief Configures grid coordinates and initializes grid management structures.
 *
 * Sets up the cell-centered DM (fda) with ghost boundaries and the coordinate DM (da).
 * Then assigns coordinates either from a file or via programmatic generation.
 */
PetscErrorCode DefineGridCoordinates(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscInt       rank;
    PetscInt       bi;
    PetscInt       IM, JM, KM;
    FILE          *fd = NULL;
    PetscInt       generate_grid = 0;
    PetscInt       grid1d        = 0;
    PetscInt       block_number  = 1;
    PetscReal      L_x, L_y, L_z;
    PetscInt      *imm = NULL, *jmm = NULL, *kmm = NULL;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "DefineGridCoordinates - rank %d.\n", rank);

    ierr = ParseGridInputs(user,
                           &generate_grid, &grid1d,
                           &L_x, &L_y, &L_z,
                           &imm, &jmm, &kmm,
                           &block_number,
                           fd);
    CHKERRQ(ierr);

    if (block_number <= 0) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                "DefineGridCoordinates - block_number must be > 0.");
    }

    for (bi = 0; bi < block_number; bi++) {
        ierr = DetermineGridSizes(bi, user, &IM, &JM, &KM, fd,
                                  generate_grid, imm, jmm, kmm,
                                  &block_number);
        CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL, LOG_INFO, "DefineGridCoordinates - IM,JM,KM = %d,%d,%d\n", IM, JM, KM);

        // 1) Create the cell-centered DM (fda) + coordinate DM (da)
        ierr = InitializeGridDM(user, L_x, L_y, L_z); CHKERRQ(ierr);

        // 2) Assign coordinates into da
        ierr = AssignGridCoordinates(user,
                                     generate_grid, grid1d,
                                     IM, JM, KM,
                                     L_x, L_y, L_z,
                                     fd);
        CHKERRQ(ierr);

        // Log subdomain ownership for the cell-centered DM
        DMDALocalInfo info;
        ierr = DMDAGetLocalInfo(user->fda, &info); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_INFO,
            "DefineGridCoordinates - rank %d block %d: xs=%d, xe=%d, ys=%d, ye=%d, zs=%d, ze=%d.\n",
            rank, bi, info.xs, info.xs+info.xm, info.ys, info.ys+info.ym, info.zs, info.zs+info.zm);
    }

    ierr = FinalizeGridSetup(generate_grid, fd, imm, jmm, kmm); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "DefineGridCoordinates - Completed on rank %d.\n", rank);
    return 0;
}

/**
 * @brief Computes the local bounding box of the grid on the current process.
 *
 * Iterates over the local coordinate DM (user->da) to find min/max x,y,z.
 *
 * @param[in]  user      Pointer to the user-defined context containing grid info.
 * @param[out] localBBox Stores the local bounding box.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ComputeLocalBoundingBox(UserCtx *user, BoundingBox *localBBox)
{
    PetscErrorCode ierr;
    PetscInt       i, j, k, rank;
    DMDALocalInfo  info;
    Vec            coordinates;
    Cmpnts       ***/**coordArray*/coordArray;
    Cmpnts         minCoords, maxCoords;

    LOG_ALLOW(GLOBAL, LOG_INFO, "ComputeLocalBoundingBox: Enter.\n");

    if (!user || !localBBox) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: Null pointers.\n");
        return PETSC_ERR_ARG_NULL;
    }

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // CHANGED: Now use user->da to get local info and coordinate vector
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    PetscInt xs = info.xs, xe = xs + info.xm;
    PetscInt ys = info.ys, ye = ys + info.ym;
    PetscInt zs = info.zs, ze = zs + info.zm;

    ierr = DMGetCoordinatesLocal(user->da, &coordinates); CHKERRQ(ierr);
    if (!coordinates) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "ComputeLocalBoundingBox: No coordinate vector.\n");
        return PETSC_ERR_ARG_NULL;
    }

    ierr = DMDAVecGetArrayRead(user->da, coordinates, &coordArray); CHKERRQ(ierr);

    minCoords.x = minCoords.y = minCoords.z = PETSC_MAX_REAL;
    maxCoords.x = maxCoords.y = maxCoords.z = PETSC_MIN_REAL;

    for (k = zs; k < ze; k++) {
        for (j = ys; j < ye; j++) {
            for (i = xs; i < xe; i++) {
                Cmpnts c = coordArray[k][j][i];
                if (c.x < minCoords.x) minCoords.x = c.x;
                if (c.y < minCoords.y) minCoords.y = c.y;
                if (c.z < minCoords.z) minCoords.z = c.z;

                if (c.x > maxCoords.x) maxCoords.x = c.x;
                if (c.y > maxCoords.y) maxCoords.y = c.y;
                if (c.z > maxCoords.z) maxCoords.z = c.z;
            }
        }
    }

    ierr = DMDAVecRestoreArrayRead(user->da, coordinates, &coordArray); CHKERRQ(ierr);

    localBBox->min_coords = minCoords;
    localBBox->max_coords = maxCoords;
    user->bbox = *localBBox;

    LOG_ALLOW(GLOBAL, LOG_INFO,
        "ComputeLocalBoundingBox: rank %d local bounding box min=(%.6f,%.6f,%.6f), max=(%.6f,%.6f,%.6f).\n",
        rank, minCoords.x, minCoords.y, minCoords.z,
              maxCoords.x, maxCoords.y, maxCoords.z);

    return 0;
}

/**
 * @brief Gathers local bounding boxes from all MPI processes to rank 0.
 *
 * This calls ComputeLocalBoundingBox on each rank, then gathers them.
 */
PetscErrorCode GatherAllBoundingBoxes(UserCtx *user, BoundingBox **allBBoxes)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank, size;
    BoundingBox   *bboxArray = NULL;
    BoundingBox    localBBox;

    LOG_ALLOW(GLOBAL, LOG_INFO, "GatherAllBoundingBoxes: Enter.\n");
    if (!user || !allBBoxes) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "GatherAllBoundingBoxes: Null pointers.\n");
        return PETSC_ERR_ARG_NULL;
    }

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    ierr = ComputeLocalBoundingBox(user, &localBBox); CHKERRQ(ierr);

    if (rank == 0) {
        bboxArray = (BoundingBox *)malloc(size * sizeof(BoundingBox));
        if (!bboxArray) {
            LOG_ALLOW(LOCAL, LOG_ERROR, "GatherAllBoundingBoxes: malloc failed.\n");
            return PETSC_ERR_MEM;
        }
    }

    ierr = MPI_Gather(&localBBox, sizeof(BoundingBox), MPI_BYTE,
                      bboxArray, sizeof(BoundingBox), MPI_BYTE,
                      0, PETSC_COMM_WORLD);
    if (ierr != MPI_SUCCESS) {
        if (rank == 0 && bboxArray) free(bboxArray);
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "MPI_Gather failed in GatherAllBoundingBoxes.");
    }

    if (rank == 0) {
        *allBBoxes = bboxArray;
    } else {
        *allBBoxes = NULL;
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "GatherAllBoundingBoxes: Exit.\n");
    return 0;
}

/**
 * @brief Broadcasts the bounding box info from rank 0 to others.
 */
PetscErrorCode BroadcastAllBoundingBoxes(UserCtx *user, BoundingBox **bboxlist)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank, size;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    if (rank != 0) {
        *bboxlist = (BoundingBox*)malloc(size * sizeof(BoundingBox));
        if (!*bboxlist) {
            SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_MEM, "Failed to allocate bboxlist on non-root ranks.");
        }
    }

    ierr = MPI_Bcast(*bboxlist, (int)(size * sizeof(BoundingBox)), MPI_BYTE, 0, PETSC_COMM_WORLD);
    if (ierr != MPI_SUCCESS) {
        SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_LIB, "MPI_Bcast failed in BroadcastAllBoundingBoxes.");
    }

    return 0;
}
