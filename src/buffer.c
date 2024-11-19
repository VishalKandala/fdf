
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
    ierr = DMDAGetLocalInfo(user->da, &info);
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
      } else { // grid1d
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

        }
    } // rank 0.

    // Broadcast grid coordinates to all ranks from rank 0 (collective operation that needs to be called on the entire communicator)
    if (grid1d) {
    ierr = MPI_Bcast(gc, IM + JM + KM, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Broadcasted 1D grid coordinates.\n");
    } else {
    ierr = MPI_Bcast(gc, 3 * IM * JM * KM, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_DEBUG, "AssignGridCoordinates - Broadcasted 3D grid coordinates.\n");
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
    } else {  // grid1d
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
