
/**
 * @brief Performs trilinear interpolation of velocities from grid to particles using interpolation coefficients.
 *
 * This function interpolates velocities for particles based on the velocity field defined on the computational grid.
 * It retrieves cell indices for each particle from the "DMSwarm_CellID" field and assigns the interpolated velocity
 * using trilinear interpolation from the surrounding grid cells. The function handles boundary checks and ensures
 * that particles outside the valid grid range are appropriately managed.
 *
 * Key Steps:
 * 1. Retrieve the number of local particles and their associated data fields (cell indices and velocities).
 * 2. Map the global velocity field (`Ucat`) to the local portion for efficient access.
 * 3. Retrieve interpolation coefficients (`a1`, `a2`, `a3`) for each particle.
 * 4. Compute trilinear interpolation weights and interpolate velocities from the surrounding grid cells.
 * 5. Handle edge cases where particles may lie outside the valid grid range.
 * 6. Restore all data fields and ensure PETSc arrays and vectors are correctly finalized.
 *
 * @param[in] user Pointer to the `UserCtx` structure containing:
 *                 - `user->da`: DMDA for the grid.
 *                 - `user->swarm`: DMSwarm for particles.
 *                 - `user->Ucat`: Global velocity vector on the grid.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InterpolateParticleVelocities(UserCtx *user) {
    PetscErrorCode ierr;
    PetscInt n_local;             // Number of local particles
    PetscInt64 *cellIDs = NULL;     // Array to store cell indices for each particle
    PetscReal *velocities = NULL; // Array to store particle velocities
    PetscReal *weights = NULL;    // Array to store interpolation weights
    Cmpnts ***ucat;               // 3D array to map local grid velocities
    PetscInt i,j,k;
    DM fda = user->fda;           // Field DA 
    DM swarm = user->swarm;       // DMSwarm for the particles


    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Starting particle velocity interpolation.\n");

    // Verify global velocity field
    PetscReal max_val;
    ierr = VecNorm(user->Ucat, NORM_INFINITY, &max_val); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Global velocity field Ucat maximum magnitude: %f \n", max_val);

    // Retrieve the number of local particles
    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Found %d local particles.\n", n_local);

    // Retrieve particle data fields
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr); // Ensure 'weight' field exists

    // Access the local portion of the global velocity vector (Ucat) using 'fda'
    ierr = DMDAVecGetArrayRead(user->fda,user->Ucat,&ucat);


    LOG_ALLOW(GLOBAL, LOG_DEBUG, "InterpolateParticleVelocities: Starting velocity assignment for particles.\n");

    // Log grid dimensions
    LOG_ALLOW(GLOBAL, LOG_INFO, "Grid dimensions: mx=%d, my=%d, mz=%d \n", user->info.mx, user->info.my, user->info.mz);

    // Loop over all local particles
    for (PetscInt p = 0; p < n_local; p++) {
        // Retrieve cell indices for the particle
        i = cellIDs[3 * p];
        j = cellIDs[3 * p + 1];
        k = cellIDs[3 * p + 2];

        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Particle %d: Host Cell = (%d, %d, %d)\n", p, i, j, k);

        // Validate cell indices (boundary check)
        if (i < 0 || j < 0 || k < 0 || i >= user->info.mx || j >= user->info.my || k >= user->info.mz) {
            LOG_ALLOW(GLOBAL, LOG_WARNING, "Particle %d has invalid cell indices (%d, %d, %d)\n. Skipping interpolation.\n", p, i, j, k);
            velocities[3 * p    ] = 0.0;
            velocities[3 * p + 1] = 0.0;
            velocities[3 * p + 2] = 0.0;
            continue;
        }

        // Assign interpolated velocity to the particle

        velocities[3 * p]     = ucat[k][j][i].x; // u-component
        velocities[3 * p + 1] = ucat[k][j][i].y; // v-component
        velocities[3 * p + 2] = ucat[k][j][i].z; // w-component

        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Particle %d: Interpolated velocity from (x=%f, y=%f, z=%f) to (x=%f, y=%f, z=%f).\n",
            p, ucat[k][j][i].x, ucat[k][j][i].y, ucat[k][j][i].z,velocities[3 * p], velocities[3 * p + 1], velocities[3 * p + 2]);
    }

    // Restore the local velocity array
    ierr = DMDAVecRestoreArrayRead(fda,user->Ucat, &ucat); CHKERRQ(ierr);

    // Restore particle data fields
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Particle velocity interpolation completed.\n");
    return 0;
}


PetscErrorCode UpdateCartesianVelocity(UserCtx *user) {
    PetscErrorCode ierr;

    // Declare required variables
    Vec Coor;
    Cmpnts ***ucat = NULL, ***coor = NULL, ***centcoor = NULL;
    PetscInt i, j, k;

    LOG_ALLOW(GLOBAL, LOG_INFO, "UpdateCartesianVelocity - Starting velocity update process.\n");

    // Retrieve local DMDA grid information
    DMDALocalInfo info = user->info;

    PetscInt mx = info.mx, my = info.my, mz = info.mz;

    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;

    PetscInt nz = ze, ny = ye, nx = xe; // or the local sizes

    LOG_ALLOW(GLOBAL, LOG_INFO, "UpdateCartesianVelocity - Local subdomain ranges - xs: %d, xe: %d, ys: %d, ye: %d, zs: %d, ze: %d.\n", xs, xe, ys, ye, zs,ze);

    // Access the DMDA coordinate vector
    ierr = DMGetCoordinatesLocal(user->da, &Coor); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, Coor, &coor); CHKERRQ(ierr);

    // Allocate memory for the temporary 3D array to store interpolated coordinates
    ierr = Allocate3DArray(&centcoor, nz, ny, nx); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "UpdateCartesianVelocity - Allocated centcoor for interpolated coordinates.\n");

    // Access the Cartesian velocity vector
    ierr = DMDAVecGetArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);

    // Interpolate coordinate values from cell corners to cell centers
    LOG_ALLOW(GLOBAL, LOG_INFO, "UpdateCartesianVelocity - Interpolating coordinates from corners to centers.\n");
    ierr = InterpolateFieldFromCornerToCenter(coor, centcoor, &info); CHKERRQ(ierr);

    // For each valid "cell-center" index, set velocity of centcoor */
      {
	PetscInt iend = xe - 1; // up to the last interior cell
	PetscInt jend = ye - 1;
	PetscInt kend = ze - 1;


    // Update the Cartesian velocity at each cell center
    LOG_ALLOW(GLOBAL, LOG_INFO, "UpdateCartesianVelocity - Updating velocity values at cell centers.\n");
    for (k = lzs; k < lze; k++) {
        for (j = lys; j < lye; j++) {
            for (i = lxs; i < lxe; i++) {
	        ierr = SetLocalCartesianVelocity(&ucat[k][j][i], &centcoor[k-lzs][j-lys][i-lxs]); CHKERRQ(ierr);
                LOG_ALLOW(GLOBAL, LOG_DEBUG, "UpdateCartesianVelocity - Updated velocity at (%d, %d, %d): x=%f, y=%f, z=%f \n",
                    k, j, i, ucat[k][j][i].x, ucat[k][j][i].y, ucat[k][j][i].z);
            }
        }
    }

    // Deallocate memory for the temporary interpolated array
    ierr = Deallocate3DArray(centcoor,lze-lzs,lye-lys); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "UpdateCartesianVelocity: Deallocated centcoor.\n");

    // Restore coordinate and velocity arrays
    ierr = DMDAVecRestoreArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, Coor, &coor); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "UpdateCartesianVelocity: Completed velocity update process.\n");
    return 0;
}

PetscErrorCode InterpolateFieldFromCornerToCenter(Cmpnts ***field,
                                                  Cmpnts ***centfield,
                                                  DMDALocalInfo *info)
{
  PetscInt i, j, k;
  PetscInt xs = info->xs, xe = info->xs + info->xm;
  PetscInt ys = info->ys, ye = info->ys + info->ym;
  PetscInt zs = info->zs, ze = info->zs + info->zm;

  /* We must stop one short of the top boundary so that (i+1) etc. are in range */
  PetscInt iend = PetscMin(xe - 1, info->mx - 1); 
  PetscInt jend = PetscMin(ye - 1, info->my - 1);
  PetscInt kend = PetscMin(ze - 1, info->mz - 1);

  for (k = zs; k < kend; k++) {
    for (j = ys; j < jend; j++) {
      for (i = xs; i < iend; i++) {
        /* Average corner values for 8 corners of the cell:
             (k,j,i), (k,j,i+1), (k,j+1,i), (k,j+1,i+1),
             (k+1,j,i), etc.
         */
        centfield[k][j][i].x = (field[k][j][i].x       + field[k][j][i+1].x +
                                field[k][j+1][i].x     + field[k][j+1][i+1].x +
                                field[k+1][j][i].x     + field[k+1][j][i+1].x +
                                field[k+1][j+1][i].x   + field[k+1][j+1][i+1].x ) / 8.0;

        centfield[k][j][i].y = (field[k][j][i].y       + field[k][j][i+1].y +
                                field[k][j+1][i].y     + field[k][j+1][i+1].y +
                                field[k+1][j][i].y     + field[k+1][j][i+1].y +
                                field[k+1][j+1][i].y   + field[k+1][j+1][i+1].y ) / 8.0;

        centfield[k][j][i].z = (field[k][j][i].z       + field[k][j][i+1].z +
                                field[k][j+1][i].z     + field[k][j+1][i+1].z +
                                field[k+1][j][i].z     + field[k+1][j][i+1].z +
                                field[k+1][j+1][i].z   + field[k+1][j+1][i+1].z ) / 8.0;
      }
    }
  }

  return 0; // success
}



///////////////////////////// Single VTP Implementation ////////

/**
 * @brief Parses command-line options using PETSc's option database.
 *
 * This function initializes PETSc (via \c PetscInitialize) and extracts the
 * user-specified time index, field name, and output extension from the
 * command line. Defaults are assigned if these are not specified.
 *
 * @param[in]     argc           Number of command-line arguments.
 * @param[in]     argv           The command-line argument strings.
 * @param[out]    timeIndex      The time index (default = 0).
 * @param[out]    fieldName      Buffer to store the parsed field name (default "velocity").
 * @param[out]    outExt         Buffer to store the file extension (default "vtp").
 * @param[in]     fieldNameSize  The maximum size of \p fieldName buffer.
 * @param[in]     outExtSize     The maximum size of \p outExt buffer.
 *
 * @return PetscErrorCode  Returns 0 on success, otherwise an error code from PETSc.
 */
/*
static PetscErrorCode ParseCommandLineOptions(int argc, char **argv,
                                              int *timeIndex,
                                              char *fieldName,
                                              char *outExt,
                                              size_t fieldNameSize,
                                              size_t outExtSize)
{
  PetscFunctionBeginUser;

  PetscErrorCode ierr;

  *timeIndex = 0; // default
  PetscStrcpy(fieldName, "velocity");
  PetscStrcpy(outExt,    "vtp");

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ParseCommandLineOptions - Parsing PETSc options.\n");
  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);

  // parse
  PetscOptionsGetInt(NULL, NULL, "-ti", timeIndex, NULL);
  PetscOptionsGetString(NULL, NULL, "-field", fieldName, fieldNameSize, NULL);
  PetscOptionsGetString(NULL, NULL, "-out_ext", outExt, outExtSize, NULL);

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ParseCommandLineOptions - Completed parsing.\n");
  PetscFunctionReturn(0);
}
*/



 /*
 * @brief Main entry point for the PETSc-based VTK post-processing tool.
 *
 * This function reads coordinate and field data from distributed PETSc \c Vecs,
 * gathers them onto rank 0, and writes them as a VTK file (either \c .vts or \c .vtp).
 * The file name, time index, and field name are configured via command-line options.
 *
 * Usage:
 * \code
 *   mpirun -n X ./this_executable -ti 10 -field velocity -out_ext vtp
 * \endcode
 *
 * @param[in] argc Number of command-line arguments.
 * @param[in] argv Array of command-line argument strings.
 *
 * @return int Returns 0 on success, or an error code upon failures.
 */

/*
int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  int            rank, size;
  int            ti;
  char           field_name[64];
  char           out_ext[16];
  double        *coordsArray = NULL;
  double        *scalarArray = NULL;
  PetscInt       Ncoords     = 0;
  PetscInt       Nscalars    = 0;

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "main - Starting PETSc-based VTK post-processing.\n");

  // 1) Parse command-line options and initialize PETSc 
  ierr = ParseCommandLineOptions(argc, argv, &ti,
                                 field_name, out_ext,
                                 sizeof(field_name), sizeof(out_ext));CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
  LOG_ALLOW(GLOBAL, LOG_INFO, "main - PETSc initialized. rank=%d of %d.\n", rank, size);

  // 2) Build a user context 
  UserCtx user;

  ierr = BuildUserContext(&user);CHKERRQ(ierr);

  // 3) Read coordinate data into coordsArray 
  ierr = ReadPositions(ti, &user, &coordsArray, &Ncoords);CHKERRQ(ierr);

  // 4) Read the field data into scalarArray 
  ierr = ReadFieldDataWrapper(ti, field_name, &user, &scalarArray, &Nscalars);CHKERRQ(ierr);

  // 5) Prepare the VTKMetaData struct 
  VTKMetaData meta;
  ierr = PrepareVTKMetaData(coordsArray, Ncoords,
                            scalarArray, Nscalars,
                            field_name, out_ext,
                            &meta, rank);CHKERRQ(ierr);

  // 6) Construct the output file name 
  char outFile[256];
  ierr = ConstructOutputFilename(field_name, ti, out_ext, outFile, sizeof(outFile));CHKERRQ(ierr);

  // 7) Create the VTK file 
  if (!rank) {
    LOG_ALLOW(GLOBAL, LOG_INFO, "main - Creating VTK file '%s'.\n", outFile);
  }
  PetscErrorCode errorCode = 0;
  errorCode = CreateAndWriteVTKFile(outFile, &meta, PETSC_COMM_WORLD);
  if (errorCode) {
    PetscPrintf(PETSC_COMM_WORLD,
                "[ERROR] CreateVTKFileFromMetadata returned %d.\n", errorCode);
    LOG_ALLOW(GLOBAL, LOG_ERROR,
              "main - CreateVTKFileFromMetadata failed with code=%d.\n", errorCode);
    goto finalize;
  }
  if (!rank) {
    PetscPrintf(PETSC_COMM_SELF, "[postprocess] Wrote file: %s\n", outFile);
    LOG_ALLOW(GLOBAL, LOG_INFO, "main - Successfully wrote file: %s\n", outFile);
  }

finalize:
  // Cleanup rank-0 memory 
  if (!rank) {
    free(coordsArray);
    free(scalarArray);
    if (meta.fileType == VTK_POLYDATA) {
      free(meta.connectivity);
      free(meta.offsets);
    }
  }

  // Finalize PETSc //
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "main - Finalizing PETSc.\n");
  PetscFinalize();
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "main - Program completed.\n");
  return 0;
}


// @brief Constructs an output file path of the form "results/<fieldName><timeIndex>.<outExt>".
// 
//  @param[in]  fieldName   The name of the field (e.g., "velocity").
//  @param[in]  timeIndex   The time index to be appended (zero-padded).
//  @param[in]  outExt      The desired file extension (e.g., "vtp" or "vts").
//  @param[out] outFile     Buffer to store the resulting file path.
//  @param[in]  outFileSize The size of \p outFile buffer.
// 
//  @return PetscErrorCode  Returns 0 on success, or error from \c PetscSNPrintf.
 
static PetscErrorCode ConstructOutputFilename(const char *fieldName,
                                              PetscInt    timeIndex,
                                              const char *outExt,
                                              char       *outFile,
                                              size_t      outFileSize)
{
  PetscFunctionBeginUser;
  PetscErrorCode ierr;
  ierr = PetscSNPrintf(outFile, outFileSize,
                       "results/%s%05d.%s", fieldName, timeIndex, outExt);CHKERRQ(ierr);
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ConstructOutputFilename - Created output file '%s'.\n", outFile);
  PetscFunctionReturn(0);
}

*/ 


///////////////////////////////////////////////////


/**
 * @brief Prints the coordinates of all particles in the swarm.
 *
 * This function retrieves the local number of particles and their coordinates
 * from the DMSwarm associated with the provided UserCtx. It then prints out
 * the coordinates of each particle in a synchronized manner across all MPI processes.
 * The output includes the MPI rank, global particle ID, local particle index, and
 * the (x, y, z) coordinates.
 *
 * @param[in] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode PrintParticleCoordinates(UserCtx* user) {
    DM swarm = user->swarm;            // DMSwarm object containing particles
    PetscErrorCode ierr;               // Error code for PETSc functions
    PetscInt localNumParticles;        // Number of particles on the local MPI process
    PetscReal *coordinates;            // Array to store particle coordinates
    PetscMPIInt rank;                  // MPI rank of the current process

    // Retrieve the MPI rank of the current process
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOG_INFO,LOCAL, "PrintParticleCoordinates - Rank %d is retrieving particle coordinates.\n", rank);

    // Get the number of particles in the local swarm
    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOG_DEBUG,LOCAL, "PrintParticleCoordinates - Rank %d has %d particles.\n", rank, localNumParticles);

    // Access the 'position' field from the DMSwarm
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&coordinates); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOG_DEBUG,LOCAL, "PrintParticleCoordinates - Retrieved 'position' field.\n");

    // Iterate over each local particle and print its coordinates
    for (PetscInt i = 0; i < localNumParticles; i++) {
        // Calculate the global particle ID (assuming particles are evenly distributed)
        PetscInt64 globalParticleID = rank * localNumParticles + i + 1;

        // Synchronized printing to ensure orderly output across MPI processes
        ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,
            "Rank %d - Global Particle %" PetscInt64_FMT " - Local Particle %d : Coordinates = (%.6f, %.6f, %.6f)\n",
            rank, globalParticleID, i + 1,
            coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2]); CHKERRQ(ierr);
    }

    // Flush the synchronized output to ensure all messages are printed
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOG_DEBUG,LOCAL, "PrintParticleCoordinates - Completed printing coordinates on Rank %d.\n", rank);

    // Restore the 'position' field to clean up
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&coordinates); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOG_DEBUG, LOCAL,"PrintParticleCoordinates - Restored 'position' field.\n");

    return 0;
}


//////////////////////////////////


/**
 * @brief Prints the positions and associated metadata of all particles in the swarm.
 *
 * This function retrieves the local number of particles, their positions,
 * unique identifiers, and the MPI rank from the DMSwarm associated with the provided UserCtx.
 * It then prints out the positions of each particle along with their IDs and ranks
 * in a synchronized manner across all MPI processes. The output includes the MPI rank,
 * global particle ID, local particle index, position coordinates, and associated metadata.
 *
 * @param[in] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode PrintParticleFields(UserCtx* user) {
    DM swarm = user->swarm;                // DMSwarm object containing particles
    PetscErrorCode ierr;                   // Error code for PETSc functions
    PetscInt localNumParticles;            // Number of particles on the local MPI process
    PetscReal *positions;                  // Array to store particle positions.
    PetscInt64 *particleIDs;               // Array to store particle unique IDs.
    PetscMPIInt *particleRanks;            // Array to store particle MPI ranks.
    PetscMPIInt rank;                      // MPI rank of the current process.
    PetscInt64  *cellIDs;                  // Array to store (host)cell IDs of particles.
    PetscReal *weights;                    // Array to store particle weights.
    PetscReal *velocities;                    // Array to store particle velocities.

    // Retrieve the MPI rank of the current process
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(LOG_INFO, GLOBAL,"PrintParticleFields - Rank %d is retrieving particle positions.\n", rank);

    // Get the number of particles in the local swarm
    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_ALLOW(LOG_DEBUG, GLOBAL,"PrintParticleFields - Rank %d has %d particles.\n", rank, localNumParticles);

    // Access the 'position', 'DMSwarm_pid','DMSwarm_rank','DMSwarm_CellID' and 'weights' fields from the DMSwarm
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    LOG_ALLOW(LOG_DEBUG, LOCAL,"PrintParticleFields - Retrieved 'position' field.\n");

    ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&particleIDs); CHKERRQ(ierr);
    LOG_ALLOW(LOG_DEBUG, LOCAL,"PrintParticleFields - Retrieved 'DMSwarm_pid' field.\n");

    ierr = DMSwarmGetField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&particleRanks); CHKERRQ(ierr);
    LOG_ALLOW(LOG_DEBUG, LOCAL,"PrintParticleFields - Retrieved 'DMSwarm_rank' field.\n");

    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    LOG_ALLOW(LOG_DEBUG, LOCAL,"PrintParticleFields - Retrieved 'DMSwarm_CellID' field.\n");

    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    LOG_ALLOW(LOG_DEBUG, LOCAL,"PrintParticleFields - Retrieved 'weight' field.\n");

    ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    LOG_ALLOW(LOG_DEBUG, LOCAL,"PrintParticleFields - Retrieved 'velocity' field.\n");

    ierr = PetscPrintf(PETSC_COMM_WORLD,"___________________________________________________________________________________________________________________________________________________________________\n");
    ierr = PetscPrintf(PETSC_COMM_WORLD,"|Rank | PID | Host IDs: i,j,k |        Position : x,y,z                   |           Velocity: x,y,z                    |           Weights: a1,a2,a3            | \n");  
    // Iterate over each local particle and print its position and metadata
    for (PetscInt i = 0; i < localNumParticles; i++) {
        // Synchronized printing to ensure orderly output across MPI processes
         ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"___________________________________________________________________________________________________________________________________________________________________\n");
         ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,
				       "|  %d  |  %d  |   %d,  %d,  %d   | %.6f, %.6f, %.6f | %.6f, %.6f, %.6f |  %.6f, %.6f, %.6f  |\n",particleRanks[i],particleIDs[i],cellIDs[3 * i],cellIDs[3 * i + 1], cellIDs[3 * i + 2],positions[3 * i], positions[3 * i + 1], positions[3 * i + 2],velocities[3 * i],velocities[3 * i + 1], velocities[3 * i + 2],weights[3 * i],weights[3 * i + 1], weights[3 * i + 2]); CHKERRQ(ierr);
    }
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"___________________________________________________________________________________________________________________________________________________________________\n");
    // Add a blank line after each rank's output
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);     

    // Flush the synchronized output to ensure all messages are printed
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOG_DEBUG,GLOBAL,"PrintParticleFields - Completed printing positions on Rank %d.\n", rank);

    // Restore the 'position', 'DMSwarm_pid','DMSwarm_rank','DMSwarm_CellID' fields to clean up
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&particleIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&particleRanks); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);

    LOG_ALLOW(LOG_DEBUG, LOCAL,"PrintParticleFields - Restored all particle fields.\n");

    return 0;
}

///////////////////////////////

// -----------------------------------------------------------------------------
// Interpolation: Corner -> Center (vector)
// -----------------------------------------------------------------------------
/**
 * @brief Safely Interpolatete a vector field from corner nodes (from the coordinate DM)
 *        to cell centers (from the cell-centered DM) using the provided UserCtx.
 *
 * For each cell center in the physical region of the cell-centered DM (fda), this function
 * averages the 8 surrounding corner values from the coordinate DM (da). The coordinate DM
 * (da) is built on corners (IM+1 x JM+1 x KM+1) while the cell-centered DM (fda) covers
 * the physical cells (IM x JM x KM). Index offsets are adjusted using DMDAGetLocalInfo.
 *
 * @param[in]  field     3D array of corner-based vector data (from user->da).
 * @param[out] centfield 3D array for the Interpolateted cell-center vector data (for user->fda).
 * @param[in]  user      User context containing:
 *                       - da  : DM for the coordinate (corner) data.
 *                       - fda : DM for the cell-centered data.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCornerToCenter_Vector(
    Cmpnts ***field,         /* coordinate DM array (da, ghosted) */
    Cmpnts ***centfield,     /* output: Interpolateted Interior cell–center values for fda */
    UserCtx *user)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
    "InterpolateFieldFromCornerToCenter_Vector - Rank %d starting Interpolatetion.\n", rank);

  /* Get local info for da (owned cells) */
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

  /* For da, get ghost region info instead of owned-only info */
  PetscInt gxs, gys, gzs, gxm, gym, gzm;
  ierr = DMDAGetGhostCorners(user->da, &gxs, &gys, &gzs, &gxm, &gym, &gzm); CHKERRQ(ierr);
  LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
    "Rank %d -> DM-da Ghost Corners: gxs=%d, gxm=%d, gys=%d, gym=%d, gzs=%d, gzm=%d\n",
    rank, gxs, gxm, gys, gym, gzs, gzm);

  /* Compute physical region for fda */
  PetscInt xs = info.xs, xe = info.xs + info.xm;
  PetscInt ys = info.ys, ye = info.ys + info.ym;
  PetscInt zs = info.zs, ze = info.zs + info.zm;
  PetscInt mx = info.mx, my = info.my, mz = info.mz;
  
  PetscInt lxs = xs,lxe = xe;
  PetscInt lys = ys,lye = ye;
  PetscInt lzs = zs,lze = ze;
  
  /*
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  */

  for (PetscInt k = lzs; k < lze; k++) {
    for (PetscInt j = lys; j < lye; j++) {
      for (PetscInt i = lxs; i < lxe; i++) {

        Cmpnts sum = {0.0, 0.0, 0.0};
        PetscInt count = 0;

        /* For cell center (i,j,k) in fda, its 8 surrounding corners in da have global indices
           from (i-1,j-1,k-1) to (i,j,k). Use ghost region information from da.
         */
        for (PetscInt dk = 0; dk < 2; dk++) {
          for (PetscInt dj = 0; dj < 2; dj++) {
            for (PetscInt di = 0; di < 2; di++) {
             
              PetscInt ci = i - 1 + di;
              PetscInt cj = j - 1 + dj;
              PetscInt ck = k - 1 + dk;
              
              if (ci >= 0 && ci < mx &&
                  cj >= 0 && cj < my &&
                  ck >= 0 && ck < mz)
              {
                sum.x += field[ck][cj][ci].x;
                sum.y += field[ck][cj][ci].y;
                sum.z += field[ck][cj][ci].z;
                count++;

	        LOG_LOOP_ALLOW(GLOBAL,LOG_DEBUG,i*j*k,1000," Rank %d| i,j,k - %d,%d,%d |ci,cj,ck - %d,%d,%d| \n",rank,i,j,k,ci,cj,ck);
              }
	    }
	  }
	}
	
        if (count > 0) {
          // As centfield is a local array with index 0 at xs,xy,xz.
          centfield[k-info.zs][j-info.ys][i-info.xs].x = sum.x / (PetscReal)count;
          centfield[k-info.zs][j-info.ys][i-info.xs].y = sum.y / (PetscReal)count;
          centfield[k-info.zs][j-info.ys][i-info.xs].z = sum.z / (PetscReal)count;
        } else {
	  //  centfield[k-zs][j-ys][i-xs].x = 0.0;
	  //  centfield[k-zs][j-ys][i-xs].y = 0.0;
	  //  centfield[k-zs][j-ys][i-xs].z = 0.0;
          LOG_ALLOW(GLOBAL, LOG_DEBUG,
                    "Rank %d: Cell (i=%d,j=%d,k=%d) got no valid corner data.\n", rank, i, j, k);
        }
	/*
        if (count > 0) {
          centfield[k][j][i].x = sum.x / (PetscReal)count;
          centfield[k][j][i].y = sum.y / (PetscReal)count;
          centfield[k][j][i].z = sum.z / (PetscReal)count;
        } else {
          centfield[k][j][i].x = 0.0;
          centfield[k][j][i].y = 0.0;
          centfield[k][j][i].z = 0.0;
        }
        */
      }
    }
  }

  LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
            "InterpolateFieldFromCornerToCenter_Vector_Interior - Rank %d completed Interpolatetion.\n", rank);
  return 0;
}

/**
 * @brief Performs trilinear Interpolatetion of velocities from grid to particles using Interpolatetion coefficients.
 *
 * This function Interpolatetes velocities for particles based on the velocity field defined on the computational grid.
 * It retrieves cell indices for each particle from the "DMSwarm_CellID" field and assigns the Interpolateted velocity
 * using trilinear Interpolatetion from the surrounding grid cells. The function handles boundary checks and ensures
 * that particles outside the valid grid range are appropriately managed.
 *
 * Key Steps:
 * 1. Retrieve the number of local particles and their associated data fields (cell indices and velocities).
 * 2. Map the global velocity field (`Ucat`) to the local portion for efficient access.
 * 3. Retrieve Interpolatetion coefficients (`a1`, `a2`, `a3`) for each particle.
 * 4. Compute trilinear Interpolatetion weights and Interpolatete velocities from the surrounding grid cells.
 * 5. Handle edge cases where particles may lie outside the valid grid range.
 * 6. Restore all data fields and ensure PETSc arrays and vectors are correctly finalized.
 *
 * @param[in] user Pointer to the `UserCtx` structure containing:
 *                 - `user->da`: DMDA for the grid.
 *                 - `user->swarm`: DMSwarm for particles.
 *                 - `user->Ucat`: Global velocity vector on the grid.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InterpolateParticleVelocities(UserCtx *user) {
    PetscErrorCode ierr;
    PetscInt n_local;             // Number of local particles
    PetscInt64 *cellIDs = NULL;     // Array to store cell indices for each particle
    PetscReal *velocities = NULL; // Array to store particle velocities
    PetscReal *weights = NULL;    // Array to store Interpolatetion weights
    Cmpnts ***ucat;               // 3D array to map local grid velocities
    Cmpnts uPetscInterp;               // Temporary variable to store Interpolateted velocities.
    PetscInt i,j,k;
    PetscReal a1,a2,a3;           // The weights of a particles are stored here for trilinear coefficient calculation.
    DM fda = user->fda;           // Field DA 
    DM swarm = user->swarm;       // DMSwarm for the particles


    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Starting particle velocity Interpolatetion.\n");

    // Verify global velocity field
    PetscReal max_val;
    ierr = VecNorm(user->Ucat, NORM_INFINITY, &max_val); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Global velocity field Ucat maximum magnitude: %f \n", max_val);

    // Retrieve the number of local particles
    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Found %d local particles.\n", n_local);

    // Retrieve particle data fields
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr); // Ensure 'weight' field exists

    // Access the local portion of the global velocity vector (Ucat) using 'fda'
    ierr = DMDAVecGetArrayRead(user->fda,user->Ucat,&ucat);CHKERRQ(ierr);


    LOG_ALLOW(GLOBAL, LOG_DEBUG, "InterpolateParticleVelocities: Starting velocity assignment for particles.\n");

    // Log grid dimensions
    LOG_ALLOW(GLOBAL, LOG_INFO, "Grid dimensions: mx=%d, my=%d, mz=%d \n", user->info.mx, user->info.my, user->info.mz);

    // Loop over all local particles
    for (PetscInt p = 0; p < n_local; p++) {
        // Retrieve cell indices for the particle
        i = cellIDs[3 * p];
        j = cellIDs[3 * p + 1];
        k = cellIDs[3 * p + 2];

	LOG_LOOP_ALLOW(GLOBAL, LOG_DEBUG, p, 10, "Particle %d: Host Cell = (%d, %d, %d)\n", p, i, j, k);

	// Clamp i, j, k to [0..mx-2], [0..my-2], [0..mz-2]
	if (i >= user->info.mx) i = user->info.mx - 1;
	if (j >= user->info.my) j = user->info.my - 1;
	if (k >= user->info.mz) k = user->info.mz - 1;

        // Validate cell indices (boundary check)
        if (i < 0 || j < 0 || k < 0 || i >= user->info.mx || j >= user->info.my || k >= user->info.mz) {
            LOG_ALLOW(GLOBAL, LOG_WARNING, "Particle %d has invalid cell indices (%d, %d, %d)\n. Skipping Interpolatetion.\n", p, i, j, k);
            velocities[3 * p    ] = 0.0;
            velocities[3 * p + 1] = 0.0;
            velocities[3 * p + 2] = 0.0;
            continue;
        }

        // Retrieve a1, a2, a3 from the 'weights' field (if that's where you're storing them)
        a1 = weights[3*p + 0];
        a2 = weights[3*p + 1];
        a3 = weights[3*p + 2];

        
	// Apply Interpolatetion method to obtain velocity at the particle location
	//  ierr = InterpolateTrilinearVelocity(ucat,i,j,k,a1,a2,a3,&uPetscInterp);
	
        // Assign Interpolateted velocity to the particle

        // zeroth order Interpolation
        
        velocities[3 * p]     = ucat[k+1][j+1][i+1].x; // u-component
        velocities[3 * p + 1] = ucat[k+1][j+1][i+1].y; // v-component
        velocities[3 * p + 2] = ucat[k+1][j+1][i+1].z; // w-component
       
	LOG_LOOP_ALLOW(GLOBAL, LOG_DEBUG, p, 10, "Particle %d: Interpolated velocity: (velocity.x=%f, velocity.y=%f, velocity.z=%f).\n",
                        p, velocities[3 * p], velocities[3 * p + 1], velocities[3 * p + 2]);	
    }

    // Restore the local velocity array
    ierr = DMDAVecRestoreArrayRead(fda,user->Ucat, &ucat); CHKERRQ(ierr);

    // Restore particle data fields
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: Particle velocity Interpolatetion completed.\n");

    // Ensure all ranks finish Interpolatetion before proceeding
    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "InterpolateParticleVelocities: All ranks completed particle Interpolatetion.\n");

    return 0;
}


/////////////////////////////////////


/**
 * @brief Initialize the simulation context.
 *
 * Checks for the presence of "control.dat" file, reads runtime options, and sets up the user context.
 *
 * @param[out] user    Pointer to the allocated UserCtx structure.
 * @param[out] rank    MPI rank of the process.
 * @param[out] size    Number of MPI processes.
 * @param[out] np      Number of particles.
 * @param[out] StartStep Simulation Starting Timestep
 * @param[out] StepsToRun No.of Timesteps to run simulation.
 * @param[out] StartTIme Time of start of simulation.
 * @param[out] ti      The timestep to start from if restarting.
 * @param[out] nblk    Number of grid blocks.
 * @param[out] outputFreq The Frequency at which data should be output from the simulation.
 * @param[out] readFields The flag to decide if eulerian fields are read or generated.
 * @param[out] allowedFuncs list of functions that are allowed to show output
 * @param[out] nAllowed No.of functions allowed to show output
 * @param[out] allowedFile indicates the file  that contains the list of allowed functions.
 * @param[out] useCfg Flag for whether a config file is prescribed or to use default.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode InitializeSimulation(UserCtx **user, PetscMPIInt *rank, PetscMPIInt *size, PetscInt *np, PetscInt *StartStep, PetscInt *StepsToRun,PetscReal *StartTime, PetscInt *nblk, PetscInt *outputFreq, PetscBool *readFields,char ***allowedFuncs, PetscInt *nAllowed, char *allowedFile, PetscBool *useCfg) {
    PetscErrorCode ierr;

    ierr = PetscOptionsGetString(NULL,NULL,
				 "-allowed_funcs_file",
				 (char *)allowedFile,
				 PETSC_MAX_PATH_LEN,
				 useCfg);CHKERRQ(ierr);

    /*  Read the allowed functions list from disk */
    ierr = LoadAllowedFunctionsFromFile(allowedFile,allowedFuncs,nAllowed);CHKERRQ(ierr);

        
    /* 2c  Register Allowed functions with logger. */
    set_allowed_functions((const char **)*allowedFuncs,(size_t)*nAllowed);
    
    // Attempt to insert options from "control.dat"
    ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, "control.dat", PETSC_TRUE);
    if (ierr == PETSC_ERR_FILE_OPEN) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
                "InitializeSimulation - Could not open 'control.dat'. Please ensure it exists in the current directory.");
    } else {
        CHKERRQ(ierr);
    }

    // Allocate user context
    ierr = PetscCalloc1(1, user); CHKERRQ(ierr);

    // Get MPI rank and size
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, size); CHKERRQ(ierr);

    // Initialize user context flags
    (*user)->averaging = PETSC_FALSE;
    (*user)->les = PETSC_FALSE;
    (*user)->rans = PETSC_FALSE;

    LOG_ALLOW(GLOBAL, LOG_INFO, "InitializeSimulation - Initialized on rank %d out of %d processes.\n", *rank, *size);

    // Read runtime options
    ierr = PetscOptionsGetInt(NULL, NULL, "-numParticles", np, NULL); CHKERRQ(ierr);
    (*user)->NumberofParticles = *np; 
    ierr = PetscOptionsGetInt(NULL, NULL, "-rstart", StartStep, NULL); CHKERRQ(ierr);
    if((*StartStep)>0) {
        ierr = PetscOptionsGetReal(NULL, NULL, "-ti", StartTime, NULL); CHKERRQ(ierr);
    }
    ierr = PetscOptionsGetInt(NULL,NULL, "-totalsteps",StepsToRun,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-nblk", nblk, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-pinit", &((*user)->ParticleInitialization), NULL);
    ierr = PetscOptionsGetInt(NULL, NULL, "-finit", &((*user)->FieldInitialization), NULL);
    ierr = PetscOptionsGetInt(NULL, NULL, "-tio",outputFreq,NULL);
    ierr = PetscOptionsGetInt(NULL, NULL, "-logfreq", &((*user)->LoggingFrequency), NULL);

    ierr = PetscOptionsGetReal(NULL, NULL, "-dt", &((*user)->dt), NULL);
    ierr = PetscOptionsGetReal(NULL, NULL, "-uin",&((*user)->ConstantVelocity), NULL);

    // Check if user requested to read fields instead of updating
    ierr = PetscOptionsGetBool(NULL, NULL, "-read_fields", readFields, NULL); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL,LOG_DEBUG," -- Console Output Functions[Total : %d] : --\n",*nAllowed);
    
    for (int i = 0; i < *nAllowed; ++i)
      LOG_ALLOW(GLOBAL,LOG_DEBUG,"   [%2d] «%s»\n", i, (*allowedFuncs)[i]);

    if (*user->ParticleInitialization == 0) {
         if(rank == 0) { // Only rank 0 parses, then bcasts
            ierr = ParseBCSFileForInlet(user); CHKERRQ(ierr);
         } else { // Non-root ranks receive broadcasted info
            PetscInt temp_inlet_face_int_bcast;
            PetscMPIInt temp_found_flag_int;
            ierr = MPI_Bcast(&temp_found_flag_int, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
            ierr = MPI_Bcast(&temp_inlet_face_int_bcast, 1, MPI_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
            user->inletFaceDefined = (PetscBool)temp_found_flag_int;
            if (user->inletFaceDefined) {
                user->identifiedInletBCFace = (BCFace)temp_inlet_face_int_bcast;
            }
         }
    }
    
    
    // Enable PETSc default logging
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);
    
    registerEvents();   

    print_log_level();
    
    return 0;
}


/////

/**
 * @brief Performs the initial setup of particles: preliminary location, migration (if needed),
 *        final location, interpolation, and initial output.
 *
 * This function is called only at the very beginning of a simulation run (step = StartStep = 0).
 * It ensures particles are spatially sorted to their correct ranks before the first
 * velocity interpolation for advection.
 *
 * @param user         Pointer to the UserCtx structure.
 * @param currentTime  The current simulation time (should be StartTime).
 * @param step         The current step number (should be StartStep).
 * @param readFields   Flag indicating whether to read initial Eulerian fields.
 * @param bboxlist     Array of bounding boxes for all ranks.
 * @param OutputFreq   Frequency for outputting data.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformInitialSetup(UserCtx *user, PetscReal currentTime, PetscInt step,
                                   PetscBool readFields, const BoundingBox *bboxlist,
                                   PetscInt OutputFreq)
{
    PetscErrorCode ierr;
    MigrationInfo  *migrationList = NULL;
    PetscInt       migrationCount = 0;
    PetscInt       migrationListCapacity = 0; // Initial capacity for migrationList
    PetscMPIInt    rank;
    PetscFunctionBeginUser;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Performing initial particle setup procedures.\n", currentTime, step);

    // --- Optional: Preliminary Spatial Sorting / Migration at T=0 ---
    // This step is to ensure particles are on the rank that owns their initial
    // physical region BEFORE the first velocity interpolation.
    // This is particularly important if InitializeParticleSwarm distributes PIDs
    // without immediate regard to spatial decomposition (e.g., round-robin PID assignment).
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Preliminary: Locating particles for initial spatial sorting.\n", currentTime, step);
    ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr); // Locates based on initial physical positions
 
    LOG_ALLOW(GLOBAL,LOG_INFO," Th intial layout of particles: \n");   
    ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Preliminary: Identifying migrating particles for spatial sorting.\n", currentTime, step);
    ierr = IdentifyMigratingParticles(user, bboxlist, &migrationList, &migrationCount, &migrationListCapacity); CHKERRQ(ierr);

    PetscInt globalMigrationCount = 0;
    ierr = MPI_Allreduce(&migrationCount, &globalMigrationCount, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
    
    if (globalMigrationCount > 0) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Preliminary: Performing spatial sort migration (%d particles).\n", currentTime, step, migrationCount);
        ierr = SetMigrationRanks(user, migrationList, migrationCount); CHKERRQ(ierr);
	LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG, "Rank %d | migrationCount %d\n",rank,migrationCount);
        ierr = PerformMigration(user); CHKERRQ(ierr);
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Preliminary: No particles identified globally for initial spatial sort migration.\n", currentTime, step);
    }
    // Reset for subsequent use in the main loop if AdvanceSimulation continues
    migrationCount = 0;
    // Note: migrationList will be freed at the end of AdvanceSimulation or if it's re-scoped.
    // If migrationList is malloc'd in IdentifyMigratingParticles, ensure it's freed if not PETSC_NULL.
    // For simplicity here, assuming IdentifyMigratingParticles handles its memory or it's managed by PetscFree later.

    // --- End of Preliminary Migration ---

    // Now, particles should be on the correct rank spatially for their initial (t=0) positions.
    // Proceed with the standard initial location (to get weights) and interpolation.
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Performing main initial locate & interpolate (post-preliminary migration).\n", currentTime, step);
    ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr); // Locate again, now they should be on the right rank and get correct cell/weights
    ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr); // This is where Vz=1.0 MUST BE OBTAINED for z=0 particles

    ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);

    // --- Initial Output ---
    // Force output if OutputFreq > 0, or consider a specific flag if needed for setup-only runs
    // that might otherwise have OutputFreq = 0.
    if (OutputFreq > 0) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Logging initial interpolation error.\n", currentTime, step);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Printing Initial particle fields at start: %.4f \n", currentTime);
        ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);
        if(user->FieldInitialization==1) {
            ierr = LOG_INTERPOLATION_ERROR(user); CHKERRQ(ierr);
        }
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Writing initial particle data (step %d).\n", currentTime, step, step);
        ierr = WriteSwarmField(user, "position", step, "dat"); CHKERRQ(ierr);
        ierr = WriteSwarmField(user, "velocity", step, "dat"); CHKERRQ(ierr);
        // ierr = WriteSwarmField(user, "pos_phy",  step, "dat"); CHKERRQ(ierr);
        if (!readFields) {
            ierr = WriteSimulationFields(user); CHKERRQ(ierr);
        }
    }

    // Cleanup memory allocated within this scope if it's not needed later
    // If migrationList was allocated here and not by IdentifyMigratingParticles for external use.
    // However, the main AdvanceSimulation expects to free it, so we leave it.
    // If IdentifyMigratingParticles reuses/reallocs migrationList, this is fine.

    PetscFunctionReturn(0);
}

/**
 * @brief Executes the main time-marching loop for the particle simulation.
 *
 * This function performs the following steps repeatedly:
 * 1. Updates/Sets the background fluid velocity field (Ucat) for the current step.
 * 2. Updates particle positions using velocity from the *previous* step's interpolation.
 *    (Note: For the very first step (step=StartStep), the velocity used might be zero
 *     or an initial guess if not handled carefully).
 * 3. Locates particles in the grid based on their *new* positions.
 * 4. Interpolates the fluid velocity (from the *current* Ucat) to the new particle locations.
 * 5. Logs errors and outputs data at specified intervals.
 *
 * @param user         Pointer to the UserCtx structure.
 * @param StartStep    The initial step number (e.g., 0 for a new run, >0 for restart).
 * @param StartTime    The simulation time corresponding to StartStep.
 * @param StepsToRun   The number of steps to execute in this run.
 * @param OutputFreq   Frequency (in number of steps) at which to output data and log errors.
 * @param readFields   Flag indicating whether to read initial fields (only at StartStep).
 * @param bboxlist     A list that contains the bounding boxes of all the ranks.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
/*
PetscErrorCode AdvanceSimulation(UserCtx *user, PetscInt StartStep, PetscReal StartTime, PetscInt StepsToRun, PetscInt OutputFreq, PetscBool readFields, const BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscReal      dt = user->dt;
    PetscInt       step; // Loop counter starting from StartStep
    PetscReal      currentTime;
    PetscInt       removed_local, removed_global;
    MigrationInfo  *migrationList = NULL;
    PetscInt       migrationCount = 0;
    PetscInt       migrationListCapacity = 0;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting simulation run: %d steps from step %d (t=%.4f), dt=%.4f\n",
              StepsToRun, StartStep, StartTime, dt);

    currentTime = StartTime; // Initialize time

    // --- Time Marching Loop ---
    // Loop from the starting step up to (but not including) the final step number
    for (step = StartStep; step < StartStep + StepsToRun; step++)
    {
      //---------------------------------------------
      PetscPrintf(PETSC_COMM_WORLD, " Current Time: %.4f | Start Time: %.4f |Current Time Step: %d | Start Time Step: %d \n",currentTime,StartTime,step,StartStep);
      LOG_ALLOW(GLOBAL, LOG_INFO, "Starting step %d, t=%.4f\n", step, currentTime);
      // Update Current Time step in user.
        // --- Step 1: Set/Update Eulerian Fields for the START of this step (time = currentTime) ---
        // This handles initialization on step==StartStep OR updates for subsequent steps
        ierr = SetEulerianFields(user, step, StartStep, currentTime, readFields); CHKERRQ(ierr);

        // --- Step 2: Initial Locate & Interpolate IF it's the very first step ---
        // This is needed *only* if step == StartStep to get the velocity BEFORE the first position update.
        // If restarting (StartStep > 0), assume velocity was read or is already correct.
        if (step == StartStep) {
             LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Performing initial locate & interpolate.\n", currentTime, step);

	     // --- << NEW: PRELIMINARY SPATIAL SORTING / MIGRATION AT T=0 (Step 0) >> ---
	     // Particles are at their initial physical positions (e.g., z=0.0),
	     // but might be on a rank that doesn't own that spatial region based on Y-decomposition.
	     // This first sequence ensures particles are on the correct rank *before* initial velocity interpolation.
	     LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Preliminary: Locating particles for initial migration.\n", currentTime, step);
	     ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr); // Locates based on initial physical positions

	     LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Preliminary: Identifying migrating particles.\n", currentTime, step);
	     ierr = IdentifyMigratingParticles(user, bboxlist, &migrationList, &migrationCount, &migrationListCapacity); CHKERRQ(ierr);
          
	     //if (migrationCount > 0) { // It's good practice to check, but for initial sort, there will likely be migrations
	     LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Preliminary: Performing migration (%d particles).\n", currentTime, step, migrationCount);
	     ierr = SetMigrationRanks(user, migrationList, migrationCount); CHKERRQ(ierr);
	     ierr = PerformMigration(user); CHKERRQ(ierr);
	     migrationCount = 0; // Reset for the main loop's migration
	     //}
	     // --- << END OF PRELIMINARY MIGRATION >> ---
	     LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Performing initial locate & interpolate (post-preliminary migration).\n", currentTime, step);
	     
             ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);
             ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr);

	     // Scatter Particle Fields to Euler Field.
	     ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);
	     
             // --- Initial Output ---
             if (OutputFreq > 0) { // Always output initial state if output is on
                 LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Logging initial interpolation error.\n", currentTime, step);
		 LOG_ALLOW(GLOBAL, LOG_DEBUG, "Printing Initial particle fields at start: %.4f \n", StartTime);
		 ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency);
		 if(user->FieldInitialization==1) ierr = LOG_INTERPOLATION_ERROR(user); CHKERRQ(ierr);
                 LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Writing initial particle data (step %d).\n", currentTime, step, step);
		 
                 ierr = WriteSwarmField(user, "position", step, "dat"); CHKERRQ(ierr);
                 ierr = WriteSwarmField(user, "velocity", step, "dat"); CHKERRQ(ierr);
		 // ierr = WriteSwarmField(user, "pos_phy",  step, "dat"); CHKERRQ(ierr);
                 if (!readFields) { ierr = WriteSimulationFields(user); CHKERRQ(ierr); } // Optional initial grid write
             }
        }
     
	// --- Step 3: Update Particle Positions ---
        // Moves particle from currentTime to currentTime + dt using velocity interpolated at currentTime
        LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f -> T=%.4f, Step=%d] Updating particle positions(Rank: %d).\n", currentTime, currentTime+dt, step,rank);
        ierr = UpdateAllParticlePositions(user); CHKERRQ(ierr); // P(t+dt) = P(t) + V(t)*dt
	
        // --- Step 4: Check Boundaries and Remove Out-Of-Bounds Particles ---
	ierr = CheckAndRemoveOutOfBoundsParticles(user,&removed_local,&removed_global);
	if (removed_global > 0) {
            LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Removed %d out-of-bounds particles.\n", currentTime, step, removed_global);
	}

	// --- Step 5: Migrate Particles between processors (if any identified) --
	ierr = IdentifyMigratingParticles(user, bboxlist, &migrationList, &migrationCount, &migrationListCapacity); CHKERRQ(ierr);

        // --- Step 6: Migrate Particles Between Processors (if any identified) ---
	//   if (migrationCount > 0) {
            ierr = SetMigrationRanks(user, migrationList, migrationCount); CHKERRQ(ierr);
            ierr = PerformMigration(user); CHKERRQ(ierr);
	    migrationCount = 0; // Reset count
	    // }
	
        // --- Step 7: Advance Time ---
        // Now currentTime represents the time at the *end* of the step completed
        currentTime += dt;
	
        // --- Step 8: Locate Particles at New Positions (t = currentTime) ---
        LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d completed] Locating particles at new positions.(Rank: %d)\n", currentTime, step,rank);
        ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr); // Finds cell/weights for P(t+dt)

        // --- Step 9: Interpolate Field at New Positions (for the *next* step's update) ---
        LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d completed] Interpolating field at new particle positions.(Rank: %d)\n", currentTime, step,rank);
        ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr); // Calculates V_interp @ P(t+dt)

	// --- Step 10: Scatter Fields from particles to grid
	ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);
	LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d completed] Scattering particle fields to grid at new particle positions.(Rank: %d)\n", currentTime, step,rank);

	user->step = step + 1;
        // --- Step 10: Output and Error Logging (based on the step *just completed*) ---
        // Output frequency check uses 'step+1' if we want output *after* completing step 'n*OutputFreq'
        // Or use 'step' if step 0 counts and we want output at 0, F, 2F, ...
        // Let's output *after* completing steps that are multiples of OutputFreq (and not step 0 again if already done)
        if (OutputFreq > 0 && (step + 1) % OutputFreq == 0) {
            LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d completed] Logging interpolation error.\n", currentTime, step);
            if(user->FieldInitialization==1)ierr = LOG_INTERPOLATION_ERROR(user); CHKERRQ(ierr); // Compares V_interp @ P(t+dt) with V_analytic @ P(t+dt)
	    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Printing particle fields at t: %.4f \n", currentTime);
	    ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency);

            LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d completed] Writing particle data (step %d).\n", currentTime, step, step+1); // Save state for start of next step
            ierr = WriteSwarmField(user, "position", step + 1, "dat"); CHKERRQ(ierr); // Write P(t+dt)
            ierr = WriteSwarmField(user, "velocity", step + 1, "dat"); CHKERRQ(ierr); // Write V_interp @ P(t+dt)
	    // ierr = WriteSwarmField(user, "pos_phy",  step + 1, "dat"); CHKERRQ(ierr);
            ierr = WriteSimulationFields(user); CHKERRQ(ierr); // Optional grid field write
        }

        LOG_ALLOW(GLOBAL, LOG_INFO, "Finished step %d, t=%.4f\n", step, currentTime);
	

    } // End for loop

    PetscReal finalTime = StartTime + StepsToRun * dt; // Calculate actual final time
    LOG_ALLOW(GLOBAL, LOG_INFO, "Time marching completed. %d steps run from step %d. Final time t=%.4f\n", StepsToRun, StartStep, finalTime);

    // --- Final Particle Field Log ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Printing final particle fields (t=%.4f):\n", finalTime);
    ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency);

    // --- Cleanup migration list
    ierr = PetscFree(migrationList); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
*/

/**
 * @brief Executes the main time-marching loop for the particle simulation.
 *
 * This function performs the following steps repeatedly:
 * 1. Updates/Sets the background fluid velocity field (Ucat) for the current step.
 * 2. If it's the very first step (StartStep=0), calls PerformInitialSetup for
 *    preliminary migration, location, interpolation, and output.
 * 3. Updates particle positions using velocity from the *previous* step's interpolation.
 * 4. Performs boundary checks and particle migration.
 * 5. Locates particles in the grid based on their *new* positions.
 * 6. Interpolates the fluid velocity (from the *current* Ucat) to the new particle locations
 *    to get velocities for the *next* advection step.
 * 7. Logs errors and outputs data at specified intervals.
 *
 * @param user         Pointer to the UserCtx structure.
 * @param StartStep    The initial step number (e.g., 0 for a new run, >0 for restart).
 * @param StartTime    The simulation time corresponding to StartStep.
 * @param StepsToRun   The number of steps to execute in this run. If 0 and StartStep is 0,
 *                     only PerformInitialSetup is executed.
 * @param OutputFreq   Frequency (in number of steps) at which to output data and log errors.
 * @param readFields   Flag indicating whether to read initial fields (only at StartStep).
 * @param bboxlist     Array of bounding boxes for all ranks.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode AdvanceSimulation(UserCtx *user, PetscInt StartStep, PetscReal StartTime,
                                 PetscInt StepsToRun, PetscInt OutputFreq, PetscBool readFields,
                                 const BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscReal      dt = user->dt;
    PetscInt       step; // Loop counter starting from StartStep
    PetscReal      currentTime;
    PetscInt       removed_local, removed_global;
    MigrationInfo  *migrationList = NULL;
    PetscInt       migrationCount = 0;
    PetscInt       migrationListCapacity = 0;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting simulation run: %d steps from step %d (t=%.4f), dt=%.4f\n",
              StepsToRun, StartStep, StartTime, dt);

    currentTime = StartTime; // Initialize time

    // --- Handle Initial Setup Separately ---
    if (StartStep == 0) { // Only do this if it's a fresh start (not a restart from step > 0)
        ierr = SetEulerianFields(user, StartStep, StartStep, currentTime, readFields); CHKERRQ(ierr);
        ierr = PerformInitialSetup(user, currentTime, StartStep, readFields, bboxlist, OutputFreq); CHKERRQ(ierr);

        if (StepsToRun == 0) { // If only initial setup was requested
            LOG_ALLOW(GLOBAL, LOG_INFO, "Initial setup completed as StepsToRun is 0. Time t=%.4f. Exiting AdvanceSimulation.\n", currentTime);
            // Note: currentTime has not been advanced by dt yet.
            PetscReal finalTime = currentTime; // Should be StartTime
            LOG_ALLOW(GLOBAL, LOG_INFO, "Time marching completed. 0 full steps run from step %d. Final time t=%.4f\n", StartStep, finalTime);
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "Printing final particle fields (t=%.4f):\n", finalTime);
            ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);
            ierr = PetscFree(migrationList); CHKERRQ(ierr); // migrationList from PerformInitialSetup if it was used
            PetscFunctionReturn(0);
        }
    }

    // --- Time Marching Loop ---
    // Loop from the starting step up to (but not including) the final step number
    for (step = StartStep; step < StartStep + StepsToRun; step++)
    {
      PetscPrintf(PETSC_COMM_WORLD, " Current Time: %.4f | Start Time: %.4f |Current Time Step: %d | Start Time Step: %d \n",currentTime,StartTime,step,StartStep);
      LOG_ALLOW(GLOBAL, LOG_INFO, "Starting step %d, t=%.4f\n", step, currentTime);

      // Set/Update Eulerian Fields:
      // If StartStep was 0, fields for t=0 were already set before the loop.
      // If StartStep > 0 (restart), or for step > StartStep, set/update fields for currentTime.
      if (step > StartStep || (step == StartStep && StartStep > 0)) {
          ierr = SetEulerianFields(user, step, StartStep, currentTime, readFields); CHKERRQ(ierr);
          // If restarting (StartStep > 0 and step == StartStep), and velocities were not read from file,
          // an interpolation might be needed here if particle positions were read but velocities were not.
          // This example assumes either velocities are read, or StartStep=0 handles initial V.
      }

      // --- Step 3: Update Particle Positions ---
      // Moves particle from currentTime to currentTime + dt using velocity interpolated at currentTime
      // (or from previous step's end for step > StartStep)
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f -> T=%.4f, Step=%d] Updating particle positions(Rank: %d).\n", currentTime, currentTime+dt, step,rank);
      ierr = UpdateAllParticlePositions(user); CHKERRQ(ierr); // P(t+dt) = P(t) + V_particle(t)*dt

      // --- Step 4: Check Boundaries and Remove Out-Of-Bounds Particles ---
      ierr = CheckAndRemoveOutOfBoundsParticles(user,&removed_local,&removed_global); CHKERRQ(ierr);
      if (removed_global > 0) {
          LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Removed %d out-of-bounds particles.\n", currentTime, step, removed_global);
      }

      // --- Step 5 & 6: Migrate Particles between processors ---
      // Migration is based on P(t+dt)
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d] Main migration phase: Identifying particles at new positions.\n", currentTime + dt, step);
      ierr = IdentifyMigratingParticles(user, bboxlist, &migrationList, &migrationCount, &migrationListCapacity); CHKERRQ(ierr);
      // if (migrationCount > 0) { // Standard PETSc practice is to call even if count is 0.
          LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d] Main migration phase: Performing migration (%d particles).\n", currentTime + dt, step, migrationCount);
          ierr = SetMigrationRanks(user, migrationList, migrationCount); CHKERRQ(ierr);
          ierr = PerformMigration(user); CHKERRQ(ierr);
          migrationCount = 0; // Reset count
      // }

      // --- Step 7: Advance Time ---
      currentTime += dt; // currentTime now represents the time at the *end* of the step completed

      // --- Step 8: Locate Particles at New Positions (t = currentTime) ---
      // This is for particles that are now on the current rank after migration.
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d completed] Locating particles at new positions.(Rank: %d)\n", currentTime, step,rank);
      ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr); // Finds cell/weights for P(t+dt)

      // --- Step 9: Interpolate Field at New Positions (for the *next* step's update) ---
      // Eulerian field Ucat is still from the beginning of the current step (time = currentTime - dt, effectively).
      // We interpolate this Ucat to P(t+dt) to get V_particle(t+dt) which will be used for advection in the *next* step.
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d completed] Interpolating field (from t=%.4f) at new particle positions.(Rank: %d)\n", currentTime, step, currentTime-dt, rank);
      ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr); // Calculates V_interp @ P(t+dt) using Ucat(t_current_step_start)

      // --- Step 10: Scatter Fields from particles to grid ---
      ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d completed] Scattering particle fields to grid at new particle positions.(Rank: %d)\n", currentTime, step,rank);

      user->step = step + 1; // Update global step counter in user context

      // --- Output and Error Logging (based on the step *just completed*) ---
      if (OutputFreq > 0 && (user->step % OutputFreq == 0) ) { // Use user->step for consistency
          LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d completed] Logging interpolation error.\n", currentTime, step);
          if(user->FieldInitialization==1) {
              ierr = LOG_INTERPOLATION_ERROR(user); CHKERRQ(ierr); // Compares V_interp @ P(t+dt) with V_analytic @ P(t+dt)
          }
          LOG_ALLOW(GLOBAL, LOG_DEBUG, "Printing particle fields at t: %.4f (step %d completed)\n", currentTime, step);
          ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);

          LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d completed] Writing particle data (for step %d).\n", currentTime, step, user->step);
          ierr = WriteSwarmField(user, "position", user->step, "dat"); CHKERRQ(ierr); // Write P(t+dt)
          ierr = WriteSwarmField(user, "velocity", user->step, "dat"); CHKERRQ(ierr); // Write V_interp @ P(t+dt)
          // ierr = WriteSwarmField(user, "pos_phy",  user->step, "dat"); CHKERRQ(ierr);
          ierr = WriteSimulationFields(user); CHKERRQ(ierr); // Optional grid field write
      }
      LOG_ALLOW(GLOBAL, LOG_INFO, "Finished step %d, t=%.4f\n", step, currentTime);
    } // End for loop

    PetscReal finalTimeRun = StartTime + StepsToRun * dt; // Target final time if all steps run
    LOG_ALLOW(GLOBAL, LOG_INFO, "Time marching completed. %d steps run from step %d. Current time t=%.4f (target final: %.4f)\n", StepsToRun, StartStep, currentTime, finalTimeRun);

    // --- Final Particle Field Log (if not already done by OutputFreq) ---
    if (!(OutputFreq > 0 && (user->step % OutputFreq == 0)) || StepsToRun == 0) { // Avoid double logging if last step was an output step
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Printing final particle fields (t=%.4f):\n", currentTime);
        ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);
    }

    // --- Cleanup migration list ---
    ierr = PetscFree(migrationList); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/////////////

/**
 * @brief Initializes basic particle properties: physical position, particle ID, and cell ID placeholder.
 *
 * This function orchestrates the initialization of particles based on `user->ParticleInitialization`.
 * It retrieves necessary grid and swarm information, then loops through particles, calling
 * helper subroutines to determine cell selection and intra-cell logical coordinates.
 * It then performs the logical-to-physical transformation and sets standard particle fields.
 *
 * @param[in,out] user               Pointer to the `UserCtx` structure.
 * @param[in]     particlesPerProcess Number of particles to initialize on this MPI process.
 * @param[in]     rand_logic_i       RNG for i-dimension tasks, range `[0.0, 1.0)`.
 * @param[in]     rand_logic_j       RNG for j-dimension tasks, range `[0.0, 1.0)`.
 * @param[in]     rand_logic_k       RNG for k-dimension tasks, range `[0.0, 1.0)`.
 * @param[in]     bboxlist           (Currently unused).
 * @return PetscErrorCode Returns `0` on success, or a PETSc error code on failure.
 */
static PetscErrorCode InitializeParticleBasicProperties(UserCtx *user,
                                                   PetscInt particlesPerProcess,
                                                   PetscRandom *rand_logic_i,
                                                   PetscRandom *rand_logic_j,
                                                   PetscRandom *rand_logic_k,
                                                   BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    DM             swarm;
    PetscReal      *positions_field = NULL, *pos_phy_field = NULL;
    PetscInt64     *particleIDs = NULL, *cellIDs_petsc = NULL;
    PetscMPIInt    rank;
    const Cmpnts   ***coor_nodes_local_array;
    Vec            Coor_local;
    DMDALocalInfo  info;
    PetscInt       xs_gnode, ys_gnode, zs_gnode;
    PetscInt       xm_onode, ym_onode, zm_onode; // Unused in this main func but fetched
    PetscInt       IM_gcells, JM_gcells, KM_gcells;

    PetscFunctionBeginUser;

    // --- 1. Input Validation and Basic Setup ---
    if (!user || !rand_logic_i || !rand_logic_j || !rand_logic_k) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "InitializeParticleBasicProperties - Null user or RNG pointer.");
    }
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    swarm = user->swarm;

    ierr = DMGetCoordinatesLocal(user->da, &Coor_local); CHKERRQ(ierr);
    if (!Coor_local) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "DMGetCoordinatesLocal for user->da returned NULL.");
    ierr = DMDAVecGetArrayRead(user->fda, Coor_local, (void*)&coor_nodes_local_array); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_INFO, "InitializeParticleBasicProperties - Rank %d: Initializing %d particles. Mode: %d. Identified Inlet (if Mode 0): %d\n",
              rank, particlesPerProcess, user->ParticleInitialization, (user->ParticleInitialization == 0) ? user->identifiedInletBCFace : -1);

    ierr = DMSwarmGetField(swarm, "position",       NULL, NULL, (void**)&positions_field); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "pos_phy",        NULL, NULL, (void**)&pos_phy_field);   CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid",    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs_petsc);  CHKERRQ(ierr);

    // --- 2. Get Grid Information (once for all particles on this rank) ---
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = DMDAGetGhostCorners(user->da, &xs_gnode, &ys_gnode, &zs_gnode, &xm_onode, &ym_onode, &zm_onode); CHKERRQ(ierr);
    ierr = DMDAGetInfo(user->da, NULL, &IM_gcells, &JM_gcells, &KM_gcells, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL); CHKERRQ(ierr);

    // --- 3. Loop Over Particles ---
    for (PetscInt p = 0; p < particlesPerProcess; p++) {
        PetscInt  ci_metric_lnode = 0, cj_metric_lnode = 0, ck_metric_lnode = 0; // Cell origin local node index
        PetscReal xi_metric_logic = 0.5, eta_metric_logic = 0.5, zta_metric_logic = 0.5; // Intra-cell logical [0,1]
        Cmpnts    phys_coords = {0.0, 0.0, 0.0}; // Default physical coordinates
        PetscBool can_place_particle = PETSC_FALSE;    // Flag: can this particle be properly placed by this rank

        // --- 3.a. Determine Cell and Intra-Cell Logical Coordinates ---
        if (user->ParticleInitialization == 0) {
            ierr = DetermineSurfaceInitializationParameters(user, &info, xs_gnode, ys_gnode, zs_gnode,
                                                            IM_gcells, JM_gcells, KM_gcells,
                                                            rand_logic_i, rand_logic_j, rand_logic_k,
                                                            &ci_metric_lnode, &cj_metric_lnode, &ck_metric_lnode,
                                                            &xi_metric_logic, &eta_metric_logic, &zta_metric_logic,
                                                            &can_place_particle); CHKERRQ(ierr);
        } else { // user->ParticleInitialization == 1
            ierr = DetermineVolumetricInitializationParameters(user, &info, xs_gnode, ys_gnode, zs_gnode,
                                                               rand_logic_i, rand_logic_j, rand_logic_k,
                                                               &ci_metric_lnode, &cj_metric_lnode, &ck_metric_lnode,
                                                               &xi_metric_logic, &eta_metric_logic, &zta_metric_logic,
                                                               &can_place_particle); CHKERRQ(ierr);
        }

        // --- 3.b. Perform Logical to Physical Transformation if placement is possible ---
        if (can_place_particle) {
            ierr = MetricLogicalToPhysical(user, coor_nodes_local_array,
                                           ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                                           xi_metric_logic, eta_metric_logic, zta_metric_logic,
                                           &phys_coords); CHKERRQ(ierr);
        } else {
            // phys_coords remains (0,0,0) or whatever default you prefer for unplaced particles
            LOG_LOOP_ALLOW(LOCAL, LOG_WARNING, p, 1,
                "InitProps - Rank %d: Particle %d - Placement parameters not determined (e.g., not on inlet face, or no owned cells). Phys default: (%.2f,%.2f,%.2f).\n",
                rank, p, phys_coords.x, phys_coords.y, phys_coords.z);
        }

        // --- 3.c. Store Common Particle Properties ---
        positions_field[3*p+0] = phys_coords.x; positions_field[3*p+1] = phys_coords.y; positions_field[3*p+2] = phys_coords.z;
        pos_phy_field[3*p+0]   = phys_coords.x; pos_phy_field[3*p+1]   = phys_coords.y; pos_phy_field[3*p+2]   = phys_coords.z;
        particleIDs[p]         = (PetscInt64)rank * particlesPerProcess + p;
        cellIDs_petsc[3*p+0]   = -1; cellIDs_petsc[3*p+1] = -1; cellIDs_petsc[3*p+2] = -1;

        // --- 3.d. Logging for this particle ---
        if (can_place_particle) {
            LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, p, (particlesPerProcess > 20 ? particlesPerProcess/10 : 1),
                "InitProps - Rank %d: PID %lld (idx %d) PLACED. Mode %d. Cell(LNodeStart):(%d,%d,%d). Logic(Metric): (%.2e,%.2f,%.2f). Phys: (%.6f,%.6f,%.6f).\n",
                rank, (long long)particleIDs[p], p, user->ParticleInitialization, ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                xi_metric_logic, eta_metric_logic, zta_metric_logic, phys_coords.x, phys_coords.y, phys_coords.z);
        } // Warning for unplaced particles already logged in step 3.b
    } // --- End of per-particle loop ---

    // --- 4. Restore Pointers ---
    ierr = DMSwarmRestoreField(swarm, "position",       NULL, NULL, (void**)&positions_field); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "pos_phy",        NULL, NULL, (void**)&pos_phy_field);   CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid",    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs_petsc);  CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, Coor_local, (void*)&coor_nodes_local_array); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_INFO, "InitializeParticleBasicProperties - Rank %d: Completed processing for %d particles.\n",
              rank, particlesPerProcess);
    PetscFunctionReturn(0);
}
