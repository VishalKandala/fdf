
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

/////////////////

/**
 * @brief Determines cell selection and intra-cell logical coordinates for surface initialization (Mode 0).
 *
 * (Keep existing detailed Doxygen comments, just ensure parameters align)
 *
 * @param[in]  user Pointer to `UserCtx` (contains `identifiedInletBCFace`, IM, JM, KM).
 * @param[in]  info Pointer to `DMDALocalInfo` for the current rank's grid portion (from user->da).
 * @param[in]  xs_gnode, ys_gnode, zs_gnode Local indices (in the ghosted array) of the first *owned node*.
 * @param[in]  IM_gcells_global, JM_gcells_global, KM_gcells_global Total number of CELLS in the global domain in I, J, K. (user->IM, user->JM, user->KM)
 * @param[in]  rand_logic_i_ptr Pointer to the RNG for i-dimension tasks [0,1).
 * @param[in]  rand_logic_j_ptr Pointer to the RNG for j-dimension tasks [0,1).
 * @param[in]  rand_logic_k_ptr Pointer to the RNG for k-dimension tasks [0,1).
 * @param[out] ci_metric_lnode_out Pointer to store the local i-node index of the selected cell's origin.
 * @param[out] cj_metric_lnode_out Pointer to store the local j-node index of the selected cell's origin.
 * @param[out] ck_metric_lnode_out Pointer to store the local k-node index of the selected cell's origin.
 * @param[out] xi_metric_logic_out Pointer to store the intra-cell logical xi-coordinate [0,1).
 * @param[out] eta_metric_logic_out Pointer to store the intra-cell logical eta-coordinate [0,1).
 * @param[out] zta_metric_logic_out Pointer to store the intra-cell logical zeta-coordinate [0,1).
 * @param[out] can_place_on_surface_out PETSC_TRUE if placement parameters were successfully determined, PETSC_FALSE otherwise.
 * @return PetscErrorCode 0 on success, or a PETSc error code.
 *
 */
static PetscErrorCode DetermineSurfaceInitializationParameters(
    UserCtx *user, DMDALocalInfo *info, /* DMDALocalInfo from user->da */
    PetscInt xs_gnode, PetscInt ys_gnode, PetscInt zs_gnode,
    PetscInt IM_gcells_global, PetscInt JM_gcells_global, PetscInt KM_gcells_global, /* Pass user->IM, user->JM, user->KM */
    PetscRandom *rand_logic_i_ptr, PetscRandom *rand_logic_j_ptr, PetscRandom *rand_logic_k_ptr,
    PetscInt *ci_metric_lnode_out, PetscInt *cj_metric_lnode_out, PetscInt *ck_metric_lnode_out,
    PetscReal *xi_metric_logic_out, PetscReal *eta_metric_logic_out, PetscReal *zta_metric_logic_out,
    PetscBool *can_place_on_surface_out)
{
    PetscErrorCode ierr = 0;
    PetscReal r_val_i_sel, r_val_j_sel, r_val_k_sel; // For storing random numbers from [0,1) RNGs
    PetscInt local_cell_idx_on_face_dim1 = 0; // Local owned cell index in the first tangential dimension of the face (0-based relative to owned cells on face)
    PetscInt local_cell_idx_on_face_dim2 = 0; // Local owned cell index in the second tangential dimension of the face (0-based)
    PetscMPIInt rank_for_logging;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank_for_logging); CHKERRQ(ierr);

    *can_place_on_surface_out = PETSC_FALSE;

    if (user->inletFaceDefined == PETSC_FALSE) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP - Rank %d: Inlet face is not defined. Cannot place particle on surface. \n", rank_for_logging);
        PetscFunctionReturn(0);
    }

    *xi_metric_logic_out = 0.5; *eta_metric_logic_out = 0.5; *zta_metric_logic_out = 0.5;
    *ci_metric_lnode_out = xs_gnode; *cj_metric_lnode_out = ys_gnode; *ck_metric_lnode_out = zs_gnode; // Default

    // Get the number of cells this rank owns in each dimension using the new GetOwnedCellRange
    // Note: GlobalNodesInDim = GlobalCellsInDim + 1
    PetscInt xs_cell_global_i, num_owned_cells_i;
    PetscInt xs_cell_global_j, num_owned_cells_j;
    PetscInt xs_cell_global_k, num_owned_cells_k;

    ierr = GetOwnedCellRange(info, 0, &xs_cell_global_i, &num_owned_cells_i); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 1, &xs_cell_global_j, &num_owned_cells_j); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 2, &xs_cell_global_k, &num_owned_cells_k); CHKERRQ(ierr);

    // A rank must own some 3D cells to be able to define a face for particle placement.
    // Check if it owns at least one cell in each dimension.
    if (num_owned_cells_i == 0 || num_owned_cells_j == 0 || num_owned_cells_k == 0) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP - Rank %d: Has zero owned cells in at least one dimension (owned cells i,j,k: %d,%d,%d). Cannot place particle on surface.\n",
                  rank_for_logging, num_owned_cells_i, num_owned_cells_j, num_owned_cells_k);
        PetscFunctionReturn(0);
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP - Rank %d: Processing inlet face %d. Owned cell counts (i,j,k): (%d,%d,%d). Global CELL counts (I,J,K): (%d,%d,%d). Ghosted node starts (xs_g,ys_g,zs_g): (%d,%d,%d). Owned cell global starts (xs_c,ys_c,zs_c): (%d,%d,%d) \n",
        rank_for_logging, user->identifiedInletBCFace, num_owned_cells_i, num_owned_cells_j, num_owned_cells_k,
        IM_gcells_global, JM_gcells_global, KM_gcells_global,
        xs_gnode, ys_gnode, zs_gnode, xs_cell_global_i, xs_cell_global_j, xs_cell_global_k);

    switch (user->identifiedInletBCFace) {
        case BC_FACE_NEG_X: // Global I-MIN face (global cell index i=0)
            // Check if this rank owns cells at the global I=0 boundary
            // info->xs is the global starting NODE index. If it's 0, this rank owns node N0, so it owns cell C0.
            if (info->xs == 0 && num_owned_cells_i > 0) { // Rank owns cell C0
                *can_place_on_surface_out = PETSC_TRUE;
                *ci_metric_lnode_out = xs_gnode;         // Cell C0 origin is the first owned node in i (local index xs_gnode)
                *xi_metric_logic_out = 1.0e-6;

                // Tangential dimensions are J and K. Select an owned cell randomly on this face.
                if (num_owned_cells_j > 0) {
                    ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val_j_sel); CHKERRQ(ierr);
                    local_cell_idx_on_face_dim1 = (PetscInt)(r_val_j_sel * num_owned_cells_j);
                    local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_j - 1);
                    *cj_metric_lnode_out = ys_gnode + local_cell_idx_on_face_dim1; // Offset from start of owned nodes in J
                } else { *can_place_on_surface_out = PETSC_FALSE; }

                if (*can_place_on_surface_out && num_owned_cells_k > 0) {
                    ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val_k_sel); CHKERRQ(ierr);
                    local_cell_idx_on_face_dim2 = (PetscInt)(r_val_k_sel * num_owned_cells_k);
                    local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_k - 1);
                    *ck_metric_lnode_out = zs_gnode + local_cell_idx_on_face_dim2; // Offset from start of owned nodes in K
                } else if (*can_place_on_surface_out) { *can_place_on_surface_out = PETSC_FALSE; }

                if (*can_place_on_surface_out) {
                    ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
                    ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
                    LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/NEG_Xi: Target. CellNode(loc lnode idx)=(%d,%d,%d). Logic(xi,et,zt)=(%.2e,%.2f,%.2f) \n", *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out, *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);
                }
            } else { LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/NEG_Xi: Rank %d not on this face (info->xs=%d or num_owned_i=%d is 0). \n", rank_for_logging, info->xs, num_owned_cells_i); }
            break;

        case BC_FACE_POS_X: // Global I-MAX face (last cell C_{IM_gcells_global - 1})
            // Check if this rank owns the last cell in I.
            // The last cell C_{IM-1} has origin N_{IM-1}. Global node indices are 0..IM.
            // Global cell indices are 0..IM-1.
            // Last cell origin global index = IM_gcells_global - 1.
            // This rank owns this cell if xs_cell_global_i <= (IM_gcells_global - 1) AND
            // (xs_cell_global_i + num_owned_cells_i -1) >= (IM_gcells_global -1)
            // Simpler: if it owns node N_{IM_gcells_global - 1} (which is the origin of the last cell)
            if ( (xs_cell_global_i + num_owned_cells_i -1 >= IM_gcells_global - 1) && (xs_cell_global_i <= IM_gcells_global -1) && num_owned_cells_i > 0) {
                 // This rank's owned cell range includes the last global cell.
                *can_place_on_surface_out = PETSC_TRUE;
                // The origin node of the last cell (global index IM_gcells_global - 1)
                // Its local node index is (IM_gcells_global - 1) - info->xs + xs_gnode
                // Or, more directly, it's the (num_owned_cells_i - 1)-th owned cell origin.
                *ci_metric_lnode_out = xs_gnode + ((IM_gcells_global - 1) - xs_cell_global_i) ;
                *xi_metric_logic_out = 1.0 - 1.0e-6;

                if (num_owned_cells_j > 0) {
                    ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val_j_sel); CHKERRQ(ierr);
                    local_cell_idx_on_face_dim1 = (PetscInt)(r_val_j_sel * num_owned_cells_j);
                    local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_j - 1);
                    *cj_metric_lnode_out = ys_gnode + local_cell_idx_on_face_dim1;
                } else { *can_place_on_surface_out = PETSC_FALSE; }

                if (*can_place_on_surface_out && num_owned_cells_k > 0) {
                    ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val_k_sel); CHKERRQ(ierr);
                    local_cell_idx_on_face_dim2 = (PetscInt)(r_val_k_sel * num_owned_cells_k);
                    local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_k - 1);
                    *ck_metric_lnode_out = zs_gnode + local_cell_idx_on_face_dim2;
                } else if (*can_place_on_surface_out) { *can_place_on_surface_out = PETSC_FALSE; }

                if (*can_place_on_surface_out) {
                    ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
                    ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
                    LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/POS_Xi: Target. CellNode(loc lnode idx)=(%d,%d,%d). Logic(xi,et,zt)=(%.2e,%.2f,%.2f) \n", *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out, *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);
                }
            } else { LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/POS_Xi: Rank %d not on this face (xs_c_i=%d, num_i=%d, IM_g=%d). \n", rank_for_logging, xs_cell_global_i, num_owned_cells_i, IM_gcells_global); }
            break;

        // --- Cases for BC_FACE_NEG_Y, BC_FACE_POS_Y, BC_FACE_NEG_Z, BC_FACE_POS_Z ---
        // --- Follow the same pattern as NEG_X and POS_X, swapping dimensions accordingly ---

        case BC_FACE_NEG_Y: // Global J-MIN face (global cell index j=0)
            if (info->ys == 0 && num_owned_cells_j > 0) { // Rank owns cell C(i,0,k) for some i,k
                *can_place_on_surface_out = PETSC_TRUE;
                *cj_metric_lnode_out = ys_gnode;
                *eta_metric_logic_out = 1.0e-6;

                if (num_owned_cells_i > 0) { // Tangential I
                    ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val_i_sel); CHKERRQ(ierr);
                    local_cell_idx_on_face_dim1 = (PetscInt)(r_val_i_sel * num_owned_cells_i);
                    local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_i - 1);
                    *ci_metric_lnode_out = xs_gnode + local_cell_idx_on_face_dim1;
                } else { *can_place_on_surface_out = PETSC_FALSE; }

                if (*can_place_on_surface_out && num_owned_cells_k > 0) { // Tangential K
                    ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val_k_sel); CHKERRQ(ierr);
                    local_cell_idx_on_face_dim2 = (PetscInt)(r_val_k_sel * num_owned_cells_k);
                    local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_k - 1);
                    *ck_metric_lnode_out = zs_gnode + local_cell_idx_on_face_dim2;
                } else if (*can_place_on_surface_out) { *can_place_on_surface_out = PETSC_FALSE; }

                if (*can_place_on_surface_out) {
                    ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr);
                    ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
                    LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/NEG_Eta: Target. CellNode(loc lnode idx)=(%d,%d,%d). Logic(xi,et,zt)=(%.2f,%.2e,%.2f) \n", *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out, *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);
                }
            } else { LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/NEG_Eta: Rank %d not on this face (info->ys=%d or num_owned_j=%d is 0).\n", rank_for_logging, info->ys, num_owned_cells_j); }
            break;

        case BC_FACE_POS_Y: // Global J-MAX face
            if ( (xs_cell_global_j + num_owned_cells_j -1 >= JM_gcells_global - 1) && (xs_cell_global_j <= JM_gcells_global -1) && num_owned_cells_j > 0) {
                *can_place_on_surface_out = PETSC_TRUE;
                *cj_metric_lnode_out = ys_gnode + ((JM_gcells_global - 1) - xs_cell_global_j);
                *eta_metric_logic_out = 1.0 - 1.0e-6;

                if (num_owned_cells_i > 0) {
                    ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val_i_sel); CHKERRQ(ierr);
                    local_cell_idx_on_face_dim1 = (PetscInt)(r_val_i_sel * num_owned_cells_i);
                    local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_i - 1);
                    *ci_metric_lnode_out = xs_gnode + local_cell_idx_on_face_dim1;
                } else { *can_place_on_surface_out = PETSC_FALSE; }

                if (*can_place_on_surface_out && num_owned_cells_k > 0) {
                    ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val_k_sel); CHKERRQ(ierr);
                    local_cell_idx_on_face_dim2 = (PetscInt)(r_val_k_sel * num_owned_cells_k);
                    local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_k - 1);
                    *ck_metric_lnode_out = zs_gnode + local_cell_idx_on_face_dim2;
                } else if (*can_place_on_surface_out) { *can_place_on_surface_out = PETSC_FALSE; }

                if (*can_place_on_surface_out) {
                    ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr);
                    ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
                    LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/POS_Eta: Target. CellNode(loc lnode idx)=(%d,%d,%d). Logic(xi,et,zt)=(%.2f,%.2e,%.2f) \n", *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out, *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);
                }
            } else { LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/POS_Eta: Rank %d not on this face (xs_c_j=%d, num_j=%d, JM_g=%d). \n", rank_for_logging, xs_cell_global_j, num_owned_cells_j, JM_gcells_global); }
            break;

        case BC_FACE_NEG_Z: // Global K-MIN face (global cell index k=0)
            if (info->zs == 0 && num_owned_cells_k > 0) { // Rank owns cell C(i,j,0) for some i,j
                *can_place_on_surface_out = PETSC_TRUE;
                *ck_metric_lnode_out = zs_gnode; // Cell C(i,j,0) origin is the first owned node in k
                *zta_metric_logic_out = 1.0e-6;

                if (num_owned_cells_i > 0) { // Tangential I
                    ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val_i_sel); CHKERRQ(ierr);
                    local_cell_idx_on_face_dim1 = (PetscInt)(r_val_i_sel * num_owned_cells_i);
                    local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_i - 1);
                    *ci_metric_lnode_out = xs_gnode + local_cell_idx_on_face_dim1;
                } else { *can_place_on_surface_out = PETSC_FALSE; }

                if (*can_place_on_surface_out && num_owned_cells_j > 0) { // Tangential J
                    ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val_j_sel); CHKERRQ(ierr);
                    local_cell_idx_on_face_dim2 = (PetscInt)(r_val_j_sel * num_owned_cells_j);
                    local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_j - 1);
                    *cj_metric_lnode_out = ys_gnode + local_cell_idx_on_face_dim2;
                } else if (*can_place_on_surface_out) { *can_place_on_surface_out = PETSC_FALSE; }

                if (*can_place_on_surface_out) {
                    ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr);
                    ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
                    LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/NEG_Zeta: Target. CellNode(loc lnode idx)=(%d,%d,%d). Logic(xi,et,zt)=(%.2f,%.2f,%.2e) \n", *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out, *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);
                }
            } else { LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/NEG_Zeta: Rank %d not on this face (info->zs=%d or num_owned_k=%d is 0). \n", rank_for_logging, info->zs, num_owned_cells_k); }
            break;

        case BC_FACE_POS_Z: // Global K-MAX face
            if ( (xs_cell_global_k + num_owned_cells_k -1 >= KM_gcells_global - 1) && (xs_cell_global_k <= KM_gcells_global -1) && num_owned_cells_k > 0) {
                *can_place_on_surface_out = PETSC_TRUE;
                *ck_metric_lnode_out = zs_gnode + ((KM_gcells_global - 1) - xs_cell_global_k);
                *zta_metric_logic_out = 1.0 - 1.0e-6;

                if (num_owned_cells_i > 0) {
                    ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val_i_sel); CHKERRQ(ierr);
                    local_cell_idx_on_face_dim1 = (PetscInt)(r_val_i_sel * num_owned_cells_i);
                    local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_i - 1);
                    *ci_metric_lnode_out = xs_gnode + local_cell_idx_on_face_dim1;
                } else { *can_place_on_surface_out = PETSC_FALSE; }

                if (*can_place_on_surface_out && num_owned_cells_j > 0) {
                    ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val_j_sel); CHKERRQ(ierr);
                    local_cell_idx_on_face_dim2 = (PetscInt)(r_val_j_sel * num_owned_cells_j);
                    local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_j - 1);
                    *cj_metric_lnode_out = ys_gnode + local_cell_idx_on_face_dim2;
                } else if (*can_place_on_surface_out) { *can_place_on_surface_out = PETSC_FALSE; }

                if (*can_place_on_surface_out) {
                    ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr);
                    ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
                    LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/POS_Zeta: Target. CellNode(loc lnode idx)=(%d,%d,%d). Logic(xi,et,zt)=(%.2f,%.2f,%.2e) \n", *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out, *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);
                }
            } else { LOG_ALLOW(LOCAL, LOG_DEBUG, "DSP/POS_Zeta: Rank %d not on this face (xs_c_k=%d, num_k=%d, KM_g=%d).\n", rank_for_logging, xs_cell_global_k, num_owned_cells_k, KM_gcells_global); }
            break;

    default:
            LOG_ALLOW(LOCAL, LOG_ERROR, "DSP - Rank %d: Invalid user->identifiedInletBCFace value: %d\n", rank_for_logging, user->identifiedInletBCFace);
            *can_place_on_surface_out = PETSC_FALSE;
            break;
    }

    if (*can_place_on_surface_out) {
        if (user->identifiedInletBCFace == BC_FACE_NEG_X || user->identifiedInletBCFace == BC_FACE_POS_X) {
            *eta_metric_logic_out = PetscMin(PetscMax(0.0, *eta_metric_logic_out), 1.0 - 1.0e-7); // Clamp [0, 1-eps]
            *zta_metric_logic_out = PetscMin(PetscMax(0.0, *zta_metric_logic_out), 1.0 - 1.0e-7);
        } else if (user->identifiedInletBCFace == BC_FACE_NEG_Y || user->identifiedInletBCFace == BC_FACE_POS_Y) {
            *xi_metric_logic_out  = PetscMin(PetscMax(0.0, *xi_metric_logic_out),  1.0 - 1.0e-7);
            *zta_metric_logic_out = PetscMin(PetscMax(0.0, *zta_metric_logic_out), 1.0 - 1.0e-7);
        } else { // Z-faces (NEG_Z or POS_Z)
            *xi_metric_logic_out  = PetscMin(PetscMax(0.0, *xi_metric_logic_out),  1.0 - 1.0e-7);
            *eta_metric_logic_out = PetscMin(PetscMax(0.0, *eta_metric_logic_out), 1.0 - 1.0e-7);
        }
    }
    PetscFunctionReturn(0);
}

/////////////////////////////////////////////////////////////////////////////

/**
 * @brief Defines the basic migration pattern for particles within the swarm.
 *
 * This function establishes the migration pattern that dictates how particles
 * move between different MPI ranks in the simulation. It initializes a migration
 * list where each particle is assigned a target rank based on predefined conditions.
 * The migration pattern can be customized to implement various migration behaviors.
 *
 * @param[in,out] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode DefineBasicMigrationPattern(UserCtx* user) {
    DM swarm = user->swarm;           // DMSwarm object managing the particle swarm
    PetscErrorCode ierr;              // Error code for PETSc functions
    PetscMPIInt *miglist;             // Migration list indicating target MPI ranks for particles
    PetscInt localNumParticles;       // Number of particles managed by the local MPI process
    PetscMPIInt rank, size;           // MPI rank of the current process and total number of processes

    // Retrieve the MPI rank of the current process
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    // Retrieve the total number of MPI processes
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL,LOG_INFO,"DefineBasicMigrationPattern - Rank %d out of %d processes.\n", rank, size);

    // Get the number of particles managed by the local MPI process
    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL,LOG_INFO,"DefineBasicMigrationPattern - Rank %d handling %d particles.\n", rank, localNumParticles);

    // Allocate memory for the migration list
    ierr = PetscCalloc1(localNumParticles, &miglist); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG,"DefineBasicMigrationPattern - Allocated migration list for %d particles.\n", localNumParticles);

    // Initialize the migration list: assign each particle to migrate to the current rank by default
    for (PetscInt p = 0; p < localNumParticles; p++) {
        miglist[p] = rank;
    }
    LOG_ALLOW(LOG_DEBUG, LOCAL,"DefineBasicMigrationPattern - Initialized migration list with default rank assignments.\n");

    // Define custom migration conditions based on the number of MPI processes
    if (size > 1) {
        // Example condition: Assign the first particle in rank 0 to migrate to rank 2
        if (rank == 0 && localNumParticles > 0) {
            miglist[0] = 2;
            LOG_ALLOW(LOG_INFO,LOCAL,"DefineBasicMigrationPattern - Rank 0, Particle 0 assigned to migrate to Rank 2.\n");
        }

        // Additional custom conditions can be added here for other ranks
        // Example:
        // if(rank == 1 && localNumParticles > 1){
        //     miglist[1] = 3;
        //     LOG_ALLOW_SYNC(LOG_INFO, "DefineBasicMigrationPattern - Rank 1, Particle 1 assigned to migrate to Rank 3.\n");
        // }

        // ... add more custom conditions as needed ...
    }

    // Assign the migration list to the user context for later use
    user->miglist = miglist;
    LOG_ALLOW(LOG_DEBUG, LOCAL,"DefineBasicMigrationPattern - Migration list assigned to user context.\n");

    /*
    // Optional: Debugging output to verify migration assignments
    for(PetscInt p = 0; p < localNumParticles; p++) {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,
            "DefineBasicMigrationPattern - Rank %d - miglist[%ld] = %ld\n",
            rank, p, miglist[p]);
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "***********************\n");
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    */

    return 0;
}


/**
 * @brief Performs the basic migration of particles based on the defined migration pattern.
 *
 * This function updates the positions of particles within the swarm by migrating them
 * to target MPI ranks as specified in the migration list. It handles the migration process
 * by setting the 'DMSwarm_rank' field for each particle and invokes the DMSwarm migration
 * mechanism to relocate particles across MPI processes. After migration, it cleans up
 * allocated resources and ensures synchronization across all MPI ranks.
 *
 * @param[in,out] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode PerformBasicMigration(UserCtx* user) {
    DM swarm = user->swarm;                // DMSwarm object managing the particle swarm
    PetscErrorCode ierr;                   // Error code for PETSc functions
    PetscMPIInt *miglist;                  // Migration list indicating target MPI ranks for particles
    PetscMPIInt *rankval;                  // Array to store current MPI rank of each particle
    PetscInt localNumParticles;            // Number of particles managed by the local MPI process
    PetscMPIInt rank;                      // MPI rank of the current process
    PetscBool removePoints = PETSC_TRUE;   // Flag indicating whether to remove migrated particles from the local swarm

    // Retrieve the MPI rank of the current process
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL,LOG_INFO,"PerformBasicMigration - Rank %d is initiating migration.\n", rank);

    // Execute the migration pattern to define target ranks for particles
    ierr = DefineBasicMigrationPattern(user); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_INFO,"PerformBasicMigration - Migration pattern defined.\n");

    // Get the number of particles in the local swarm
    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG,"PerformBasicMigration - Rank %d handling %d particles.\n", rank, localNumParticles);

    // Retrieve the migration list from the user context
    miglist = user->miglist;
    LOG_ALLOW(LOG_DEBUG,LOCAL,"PerformBasicMigration - Retrieved migration list from user context.\n");

    /*
    // Optional: Debugging output to verify migration assignments before migration
    for(PetscInt p = 0; p < localNumParticles; p++) {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,
            "PerformBasicMigration - Rank %d - miglist[%ld] = %d; user->miglist[p] = %d\n",
            rank, p, miglist[p], user->miglist[p]);
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "***********************\n");
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    */

    // Access the 'DMSwarm_rank' field from the DMSwarm to update particle ranks
    ierr = DMSwarmGetField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&rankval); CHKERRQ(ierr);
    LOG_ALLOW_SYNC( LOCAL,LOG_DEBUG,"PerformBasicMigration - Retrieved 'DMSwarm_rank' field.\n");

    // Update the 'DMSwarm_rank' field based on the migration list
    for (PetscInt p = 0; p < localNumParticles; p++) {
        rankval[p] = miglist[p];
        LOG_ALLOW_SYNC( LOCAL,LOG_DEBUG,"PerformBasicMigration - Particle %d assigned to Rank %d.\n", p, rankval[p]);
    }

    /*
    // Optional: Debugging output to verify migration assignments after rank updates
    for(PetscInt p = 0; p < localNumParticles; p++) {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,
            "PerformBasicMigration - After change - Rank %ld - rankval[%ld] = %ld; user->miglist[p] = %ld\n",
            rank, p, rankval[p], user->miglist[p]);
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, "***********************\n");
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    */

    // Restore the 'DMSwarm_rank' field after modification
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&rankval); CHKERRQ(ierr);
    LOG_ALLOW_SYNC( LOCAL,LOG_DEBUG,"PerformBasicMigration - Restored 'DMSwarm_rank' field.\n");

    // Invoke the DMSwarm migration process to relocate particles based on updated ranks
    ierr = DMSwarmMigrate(swarm, removePoints); CHKERRQ(ierr);
    LOG_ALLOW_SYNC( LOCAL,LOG_INFO,"PerformBasicMigration - DMSwarm migration executed.\n");

    // Free the allocated migration list to prevent memory leaks
    ierr = PetscFree(user->miglist); CHKERRQ(ierr);
    LOG_ALLOW_SYNC( LOCAL,LOG_DEBUG,"PerformBasicMigration - Freed migration list memory.\n");

    // Synchronize all MPI processes to ensure migration completion before proceeding
    ierr = PetscBarrier(NULL); CHKERRQ(ierr);
    LOG_ALLOW_SYNC( LOCAL,LOG_INFO, "PerformBasicMigration - Migration synchronization completed.\n");

    return 0;
}


/**
 * @brief Locates the grid cell containing a particle using a walking search.
 * @ingroup ParticleLocation
 *
 * This function implements a walking search algorithm to find the specific cell
 * (identified by global indices i, j, k) that encloses the particle's physical
 * location (`particle->loc`).
 *
 * The search starts from an initial guess cell (either the particle's previously known
 * cell or the corner of the local process domain) and iteratively steps to adjacent
 * cells based on the signed distances from the particle to the faces of the current
 * cell.
 *
 * Upon successful location (`position == 0`), the particle's `cell` field is updated
 * with the found indices (i,j,k), and the corresponding interpolation `weights`
 * are calculated and stored using the distances (`d`) relative to the final cell.
 *
 * Handles particles exactly on cell boundaries by attempting a tie-breaker if the
 * search gets stuck oscillating between adjacent cells due to the boundary condition.
 *
 * @param[in]  user     Pointer to the UserCtx structure containing grid information (DMDA, coordinates)
 *                      and domain boundaries.
 * @param[in,out] particle Pointer to the Particle structure. Its `loc` field provides the
 *                      target position. On successful return, its `cell` and `weights`
 *                      fields are updated. On failure, `cell` is set to `{-1, -1, -1}`
 *                      and `weights` to `{0.0, 0.0, 0.0}`.
 * @param[out]    status_out   The final status of the particle after the search, using the official enum.
 *
 * @return PetscErrorCode 0 on success. Non-zero error codes may indicate issues during
 *                        coordinate access, distance calculation, or other internal errors.
 *                        A return code of 0 does not guarantee the particle was found;
 *                        check `particle->cell[0] >= 0` afterward.
 *
 * @note Relies on helper functions like `InitializeTraversalParameters`, `CheckCellWithinLocalGrid`,
 *       `RetrieveCurrentCell`, `EvaluateParticlePosition`, `UpdateCellIndicesBasedOnDistances`,
 *       `UpdateParticleWeights`, and `FinalizeTraversal`.
 * @warning Ensure `particle->loc` is set correctly before calling.
 * @warning The function may fail to find the particle if it lies outside the domain accessible
 *          by the current process (including ghost cells) or if `MAX_TRAVERSAL` steps are exceeded.
 */
/*
PetscErrorCode LocateParticleInGridTEST(UserCtx *user, Particle *particle)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;

    
    PetscInt idx, idy, idz;           // Current search cell indices
    PetscInt traversal_steps;         // Counter for search steps
    PetscBool search_concluded = PETSC_FALSE;// Flag indicating end of search

    
    //const Cmpnts p = particle->loc;   // Particle position (local copy)
    Cell current_cell;                // Geometry of the cell being checked
    //const PetscReal threshold = DISTANCE_THRESHOLD; // Base threshold for distance checks
    // DMDALocalInfo info;               // Local grid info (for bounds)
    
    PetscInt repeatedIndexCount = 0;  // Counter for consecutive iterations at the same index
    PetscInt prevIdx = PETSC_MIN_INT; // Previous iteration's index i (for repeat check)
    PetscInt prevIdy = PETSC_MIN_INT; // Previous iteration's index j (for repeat check)
    PetscInt prevIdz = PETSC_MIN_INT; // Previous iteration's index k (for repeat check)
    PetscInt last_position = -999;    // Position result (-1, 0, >=1) from *previous* iteration

    PetscFunctionBeginUser;
     ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_FUNC_TIMER_BEGIN_EVENT(EVENT_Individualwalkingsearch, LOCAL);

    // Get local grid information (needed for bounds and index updates)
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

    // Determine starting cell indices (previous cell or local corner)
    ierr = InitializeTraversalParameters(user, particle, &idx, &idy, &idz, &traversal_steps); CHKERRQ(ierr);

    // --- Main Walking Search Loop ---
    while (!cell_found && traversal_steps < MAX_TRAVERSAL) {
        traversal_steps++;

        // Log current state for debugging
        LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, traversal_steps, 1,
                       "LocateParticleInGrid [PID %lld, Step %d]: Checking cell (%d, %d, %d). Pos (%.6f, %.6f, %.6f)\n",
                       (long long)particle->PID, traversal_steps, idx, idy, idz, p.x, p.y, p.z);

        // --- Check for Stuck Loop / Oscillation ---
        if (idx == prevIdx && idy == prevIdy && idz == prevIdz) {
            repeatedIndexCount++;
            LOG_ALLOW(LOCAL, LOG_DEBUG,"LocateParticleInGrid [PID %lld, Step %d]: Repeated index (%d,%d,%d) count = %d\n",
                        (long long)particle->PID, traversal_steps, idx, idy, idz, repeatedIndexCount);

            if (repeatedIndexCount > REPEAT_COUNT_THRESHOLD) {
                // --- Option B Tie-Breaker Logic ---
                if (last_position >= 1) { // Was the last step result 'on boundary'? Likely oscillation.
                    LOG_ALLOW(LOCAL, LOG_WARNING,
                              "LocateParticleInGrid [PID %lld]: Stuck at index (%d,%d,%d) after being on boundary for >%d steps. Applying boundary tie-break.\n",
                              (long long)particle->PID, idx, idy, idz, REPEAT_COUNT_THRESHOLD);

                    // Attempt to assign to the cell where we are stuck.
                    // Need to ensure the cell is valid and re-get distances for weights.
                    Cell final_cell;
                    PetscReal final_d[NUM_FACES];
                    PetscInt final_position; // Result not used, just need distances
                    PetscBool is_stuck_cell_valid;

                    ierr = CheckCellWithinLocalGrid(user, idx, idy, idz, &is_stuck_cell_valid); CHKERRQ(ierr);
                    if (is_stuck_cell_valid) {
                         ierr = RetrieveCurrentCell(user, idx, idy, idz, &final_cell); CHKERRQ(ierr);
                         // Re-evaluate to get distances 'final_d' relative to this specific cell
                         ierr = EvaluateParticlePosition(&final_cell, final_d, p, &final_position, threshold);
                         if (ierr == 0) { // Check if evaluation succeeded
                            ierr = UpdateParticleWeights(final_d, particle); CHKERRQ(ierr); // Calculate weights
                            particle->cell[0] = idx; // Assign the stuck cell index
                            particle->cell[1] = idy;
                            particle->cell[2] = idz;
                            cell_found = PETSC_TRUE; // Mark as found via tie-break
                            LOG_ALLOW(LOCAL, LOG_INFO,
                                      "LocateParticleInGrid [PID %lld]: Assigned via tie-break to cell (%d,%d,%d).\n",
                                      (long long)particle->PID, idx, idy, idz);
                         } else {
                             // Error during final evaluation
                             LOG_ALLOW(LOCAL, LOG_ERROR,
                                       "LocateParticleInGrid [PID %lld]: Error %d during final evaluation for tie-break at cell (%d,%d,%d). Search fails.\n",
                                       (long long)particle->PID, ierr, idx, idy, idz);
                             cell_found = PETSC_FALSE; // Ensure failure
                         }
                    } else {
                         LOG_ALLOW(LOCAL, LOG_WARNING,
                                   "LocateParticleInGrid [PID %lld]: Stuck at invalid cell index (%d,%d,%d) during tie-break attempt. Search fails.\n",
                                   (long long)particle->PID, idx, idy, idz);
                         cell_found = PETSC_FALSE; // Mark as failed
                    }
                    break; // Exit loop after attempting tie-break
                } else {
                    // Stuck for > threshold steps, but last known state wasn't 'on boundary'.
                    LOG_ALLOW(LOCAL, LOG_WARNING,
                              "LocateParticleInGrid [PID %lld]: Stuck at index (%d,%d,%d) for >%d steps (last_pos=%d). Breaking search (failed).\n",
                              (long long)particle->PID, idx, idy, idz, REPEAT_COUNT_THRESHOLD, last_position);
                    cell_found = PETSC_FALSE; // Mark as failure
                    break; // Exit loop (failed)
                }
            } // end if (repeatedIndexCount > threshold)
        } else { // Index changed from previous iteration
            repeatedIndexCount = 0;
        }
        // --- End Stuck Loop Check ---

        // --- Update previous index state for the *next* iteration's check ---
        prevIdx = idx;
        prevIdy = idy;
        prevIdz = idz;

        // --- Check if current cell index is within accessible grid bounds ---
        PetscBool is_within;
        ierr = CheckCellWithinLocalGrid(user, idx, idy, idz, &is_within); CHKERRQ(ierr);
        if (!is_within) {
            LOG_ALLOW(LOCAL, LOG_WARNING,
                      "LocateParticleInGrid [PID %lld]: Search moved outside local ghosted grid boundaries at cell (%d, %d, %d). Breaking search.\n",
                      (long long)particle->PID, idx, idy, idz);
            cell_found = PETSC_FALSE; // Mark as failed
            break; // Exit loop
        }

        // --- Retrieve geometry and evaluate particle position ---
        ierr = RetrieveCurrentCell(user, idx, idy, idz, &current_cell); CHKERRQ(ierr);
        PetscReal current_d[NUM_FACES]; // Holds distances relative to current_cell
        PetscInt position;              // Holds result: -1 (outside), 0 (inside), >=1 (boundary)
        ierr = EvaluateParticlePosition(&current_cell, current_d, p, &position, threshold); CHKERRQ(ierr);

        // Store the position result for the *next* iteration's stuck check
        last_position = position;

        LOG_ALLOW(LOCAL, LOG_DEBUG,
                   "LocateParticleInGrid [PID %lld, Step %d]: Evaluated pos in cell (%d,%d,%d): position=%d\n",
                   (long long)particle->PID, traversal_steps, idx, idy, idz, position);

        // --- Main Logic: Handle position result ---
        if (position == 0) { // Strictly INSIDE the cell
            cell_found = PETSC_TRUE;
            particle->cell[0] = idx; // Assign the indices of the cell we just evaluated
            particle->cell[1] = idy;
            particle->cell[2] = idz;
            // Calculate weights using the distances 'current_d' from this successful evaluation
            ierr = UpdateParticleWeights(current_d, particle); CHKERRQ(ierr);
            LOG_ALLOW(LOCAL, LOG_INFO,
                      "LocateParticleInGrid [PID %lld]: Found INSIDE cell (%d, %d, %d) after %d steps. Weights calculated.\n",
                      (long long)particle->PID, idx, idy, idz, traversal_steps);
            break; // Found the correct cell, exit loop

        } else { // Particle is OUTSIDE (position == -1) or ON BOUNDARY (position >= 1)
            LOG_ALLOW(LOCAL, LOG_DEBUG,
                      "LocateParticleInGrid [PID %lld, Step %d]: Particle OUTSIDE or ON BOUNDARY of cell (%d, %d, %d), position = %d. Updating indices.\n",
                      (long long)particle->PID, traversal_steps, idx, idy, idz, position);
            // Both cases require stepping to the next potential cell.
            // Update search indices (idx, idy, idz) based on the distances 'current_d'
            ierr = UpdateCellIndicesBasedOnDistances(current_d, &idx, &idy, &idz, &info); CHKERRQ(ierr);
            // Loop continues to the next iteration with potentially updated idx, idy, idz
        }
    } // --- End while loop ---

    // --- Handle cases where the loop terminated WITHOUT finding the cell ---
    // (Could be MAX_TRAVERSAL, out of bounds, or non-boundary stuck state)
    if (!cell_found) {
        LOG_ALLOW(LOCAL, LOG_WARNING,
                  "LocateParticleInGrid [PID %lld]: Search FAILED after %d steps. Final check was at cell (%d,%d,%d). Resetting cell/weights.\n",
                  (long long)particle->PID, traversal_steps, idx, idy, idz);
        // Ensure particle state reflects failure explicitly
        particle->cell[0] = -1;
        particle->cell[1] = -1;
        particle->cell[2] = -1;
    }

    // Finalize traversal by reporting the results (logs success/failure message)
    // Note: idx, idy, idz passed here are the *last checked* indices for logging purposes.
    // The actual found cell index is stored in particle->cell if cell_found is true.
    ierr = FinalizeTraversal(user, particle, traversal_steps, cell_found, idx, idy, idz); CHKERRQ(ierr);

    LOG_FUNC_TIMER_END_EVENT(EVENT_Individualwalkingsearch, LOCAL);
    PetscFunctionReturn(0);
}
*/

//////////////////////////////


/**
 * @brief Identifies particles that have left the local MPI rank's domain and determines their target rank.
 *
 * This function iterates through all particles currently local to this MPI rank.
 * For each particle, it first checks if the particle is still within the rank's
 * primary bounding box (defined by `user->bbox`).
 *
 * If a particle is outside `user->bbox`:
 * 1. It preferentially checks the 6 immediate Cartesian neighbors (defined in `user->neighbors`).
 *    The bounding boxes for these neighbors are looked up in the global `bboxlist`.
 * 2. If the particle is not found in any of the immediate Cartesian neighbors, a fallback
 *    search is initiated. This fallback search iterates through *all other* MPI ranks
 *    (excluding the current rank and already checked immediate neighbors if optimized)
 *    using the `bboxlist` to find a rank whose bounding box contains the particle.
 *
 * Particles successfully assigned a `targetRank` are added to the `migrationList`.
 * If a particle leaves the local box but is not found in any other rank's bounding box
 * (even after the fallback), a warning is logged, as this particle might be lost from
 * the global computational domain (assuming `bboxlist` covers the entire domain).
 * Such "lost" particles should ideally be handled by a separate global boundary condition
 * check (e.g., `CheckAndRemoveOutOfBoundsParticles`) prior to calling this function.
 *
 * The `migrationList` is dynamically reallocated if its current capacity is exceeded.
 *
 * @param[in]      user           Pointer to the UserCtx structure. It must contain:
 *                                - `swarm`: The DMSwarm object.
 *                                - `bbox`: The BoundingBox of the current MPI rank.
 *                                - `neighbors`: A RankNeighbors struct with the ranks of the 6 Cartesian neighbors.
 *                                - `size`: The total number of MPI ranks in PETSC_COMM_WORLD (implicitly via MPI_Comm_size).
 * @param[in]      bboxlist       An array of BoundingBox structures for ALL MPI ranks, indexed 0 to (size-1).
 *                                This array must be up-to-date and available on all ranks.
 * @param[in,out]  migrationList  Pointer to an array of MigrationInfo structures. This array will be
 *                                populated with particles to be migrated. The function may reallocate
 *                                this array if more space is needed. The caller is responsible for
 *                                eventually freeing this memory if it's not NULL.
 * @param[out]     migrationCount Pointer to a PetscInt that will be set to the number of particles
 *                                identified for migration (i.e., the number of valid entries in `*migrationList`).
 * @param[in,out]  listCapacity   Pointer to a PetscInt representing the current allocated capacity of
 *                                `*migrationList` (in terms of number of MigrationInfo structs).
 *                                This function will update it if `*migrationList` is reallocated.
 *
 * @return PetscErrorCode 0 on success, non-zero on PETSc/MPI errors or if essential fields are missing.
 *
 * @note It is assumed that `user->neighbors` contains valid rank identifiers or `MPI_PROC_NULL`.
 * @note It is assumed that `bboxlist` is correctly populated and broadcast to all ranks before calling this.
 * @note For the fallback search to be effective, `bboxlist` should represent a tiling of the
 *       global domain. Particles not found in any box in `bboxlist` are effectively outside this
 *       tiled global domain.
 */
/*
PetscErrorCode IdentifyMigratingParticles(UserCtx *user,
					  const BoundingBox *bboxlist,
					  MigrationInfo **migrationList,
					  PetscInt *migrationCount,
					  PetscInt *listCapacity)
{
  PetscErrorCode ierr;
  DM             swarm = user->swarm;
  PetscInt       nLocal, p;
  Cmpnts        *pos = NULL;
  PetscMPIInt    rank,size;
  BoundingBox    localBBox = user->bbox;
  RankNeighbors  neighbors = user->neighbors; // Use stored neighbors
  
  //PetscInt       currentMigrationCount = 0;
  //PetscInt       currentListCapacity = *listCapacity;
  //MigrationInfo *currentMigrationList = *migrationList;
  // Add PID pointer if logging PIDs
  PetscInt64    *pids = NULL;
  PetscFunctionBeginUser;

  // --- Input Validation and Initialization ---
  if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx 'user' is NULL.");
  if (!bboxlist) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Global 'bboxlist' is NULL.");
  if (!migrationList) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "'migrationList' output pointer is NULL.");
  if (!migrationCount) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "'migrationCount' output pointer is NULL.");
  if (!listCapacity) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "'listCapacity' output pointer is NULL.");

  if (!swarm) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "UserCtx->swarm is NULL.");
  
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
  
  ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
  if (nLocal == 0) {
    *migrationCount = 0;
    // Ensure output pointers are consistent even if no allocation happened
    // *migrationList = currentMigrationList;
    // *listCapacity = currentListCapacity;
    PetscFunctionReturn(0);
  }
  // Get read-only access to position
  ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pids); CHKERRQ(ierr); // If logging PIDs
  if (!pos || !pids) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Could not access required DMSwarm fields.");
  /*
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
    "[IdentifyMigratingParticles - Rank %d INCOMING user->neighbors] xm=%d, xp=%d, ym=%d, yp=%d, zm=%d, zp=%d. MPI_PROC_NULL is %d. \n",
    rank, user->neighbors.rank_xm, user->neighbors.rank_xp,
    user->neighbors.rank_ym, user->neighbors.rank_yp,
    user->neighbors.rank_zm, user->neighbors.rank_zp, (int)MPI_PROC_NULL);
  */
  /*
  *migrationCount = 0;
  //currentMigrationCount = 0; // Reset count for this call
  for (p = 0; p < nLocal; p++) {
    // Check if particle is OUTSIDE the local bounding box
    if (!IsParticleInBox(&localBBox, &pos[p]))
      {
	PetscInt targetRank = MPI_PROC_NULL; // Target rank not yet found

	// Determine likely exit direction(s) to prioritize neighbor check
	PetscBool exit_xm = pos[p].x < localBBox.min_coords.x;
	PetscBool exit_xp = pos[p].x > localBBox.max_coords.x;
	PetscBool exit_ym = pos[p].y < localBBox.min_coords.y;
	PetscBool exit_yp = pos[p].y > localBBox.max_coords.y;
	PetscBool exit_zm = pos[p].z < localBBox.min_coords.z;
	PetscBool exit_zp = pos[p].z > localBBox.max_coords.z;

	// DEBUG ------------
	/*
	  if (rank == 1 && p < 5) { // Log for first few particles on Rank 1
	  LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
	  "[Identify - Rank 1, p=%d] Particle Pos(%.2f,%.2f,%.2f). Exits(xm%d,xp%d,ym%d,yp%d,zm%d,zp%d). Neighbors(xm%d,xp%d,ym%d,yp%d,zm%d,zp%d)",
	  p, pos[p].x, pos[p].y, pos[p].z,
	  exit_xm, exit_xp, exit_ym, exit_yp, exit_zm, exit_zp,
	  neighbors.rank_xm, neighbors.rank_xp, neighbors.rank_ym, neighbors.rank_yp, neighbors.rank_zm, neighbors.rank_zp);
	  }
	*/
	// DEBUG -----------
        /*
	
	LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Particle %d[PID = %ld] at (%g,%g,%g) left local bbox. Checking neighbors...\n",rank, p,pids[p], pos[p].x, pos[p].y, pos[p].z);

	//1.  Check neighbors preferentially
	if (exit_xm && neighbors.rank_xm != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_xm], &pos[p])) targetRank = neighbors.rank_xm;
	else if (exit_xp && neighbors.rank_xp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_xp], &pos[p])) targetRank = neighbors.rank_xp;
	else if (exit_ym && neighbors.rank_ym != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_ym], &pos[p])) targetRank = neighbors.rank_ym;
	else if (exit_yp && neighbors.rank_yp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_yp], &pos[p])) targetRank = neighbors.rank_yp;
	else if (exit_zm && neighbors.rank_zm != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_zm], &pos[p])) targetRank = neighbors.rank_zm;
	else if (exit_zp && neighbors.rank_zp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_zp], &pos[p])) targetRank = neighbors.rank_zp;
	// Add checks for edge/corner neighbors if needed and if they were stored

	
	// 2.--- Fallback (if strict Cartesian neighbors aren't enough) ---
	
	if (targetRank == MPI_PROC_NULL){
	  /* ... loop through all ranks in bboxlist ... */
        /*
	  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Particle %" PetscInt64_FMT " not in immediate neighbors. Fallback: checking all %d other ranks...",rank, pids[p],size);

	  for(PetscMPIInt r = 0; r < size; r++){

	    if(IsParticleInBox(&bboxlist[r], &pos[p])){

	      targetRank = r;

	      LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] Particle %ld FOUND in bboxlist[%d] during fallback search. \n",rank, pids[p],r);
	      break;
	    }
	    
	  }
	  
	}
	// 3. ---- Add to migration list if a target rank was found.
	if (targetRank != MPI_PROC_NULL){

	  /* OLD LOGIC
          // Resize list if needed (using PetscRealloc for safety)
	  if (currentMigrationCount >= currentListCapacity) {
	    PetscInt OldCapacity = currentListCapacity;
	    
	      PetscInt newCapacity = (currentListCapacity == 0) ? 16 : currentListCapacity * 2;

	    ierr = PetscRealloc((size_t)newCapacity * sizeof(MigrationInfo),
				migrationList); CHKERRQ(ierr);

	    currentMigrationList = *migrationList;
	    
	    
	    currentListCapacity = newCapacity;

	    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Reallocated migrationList capacity from %d to %d\n", rank,OldCapacity, newCapacity);
	   }
	    
	  // Add to migration list
	  currentMigrationList[currentMigrationCount].local_index = p;
	  currentMigrationList[currentMigrationCount].target_rank = targetRank;
	  // currentMigrationList[currentMigrationCount].pid = pids[p]; // If storing PID
	  currentMigrationCount++;
	  */

	  // NEW LOGIC (TEST)
          /*
	  ierr = AddToMigrationList(migrationList,    // Pointer to the list pointer
				    listCapacity,     // Pointer to the capacity variable
				    migrationCount,   // Pointer to the count variable
				    p,                // The particle's local index
				    targetRank);      // The destination rank
	  CHKERRQ(ierr);
	  
	  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Particle %d marked for migration to rank %d.\n", rank, p, targetRank);
	  
	}
	  else {
	  // Particle left local box but was not found in any *checked* neighbor box.
	  // Since CheckAndRemove should have run first, this might indicate a particle
	  // moved more than one cell width into a diagonal neighbor's domain not checked here,
	  // or there's an issue with BBox overlap/gaps.
          LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: Particle %d[PID = %ld] at (%g,%g,%g) left local bbox but target neighbor rank not found! (May be lost or need wider neighbor check).\n",
                    rank, p, pids[p],pos[p].x, pos[p].y, pos[p].z);
          // Consider marking for removal here if this case should not happen:
          // Maybe add to a separate 'removalList' or use the marking technique
          // from the CheckAndRemove function if it's readily available.
	  }
      } // end if isOutsideLocal
  } // end particle loop
  ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pids); CHKERRQ(ierr);
  //  *migrationList = currentMigrationList;
  // *migrationCount = currentMigrationCount;
  //  *listCapacity = currentListCapacity;
  
  //LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Identified %d particles for potential migration.\n", rank, currentMigrationCount);
  LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Identified %d particles for potential migration.\n", rank, *migrationCount);
  PetscFunctionReturn(0);
}
*/

//////////////////////////

/*
PetscErrorCode CheckCellWithinLocalGrid(UserCtx *user, PetscInt idx, PetscInt idy, PetscInt idz, PetscBool *is_within)
{
    PetscErrorCode ierr;
    DMDALocalInfo info_nodes;
    PetscInt stencil_width;

    // Validate inputs
    if (user == NULL || is_within == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "Input pointer is NULL.\n");
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input pointer is NULL.");
    }

    // Get node info from fda (to derive cell info)
    ierr = DMDAGetLocalInfo(user->fda, &info_nodes); CHKERRQ(ierr);
    // Get stencil width used for fda (and implicitly for da's ghosts)
    ierr = DMDAGetInfo(user->fda, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &stencil_width, NULL, NULL, NULL, NULL); CHKERRQ(ierr);


    // Get owned cell ranges using the helper
    PetscInt xs_cell, xm_cell, xe_cell;
    PetscInt ys_cell, ym_cell, ye_cell;
    PetscInt zs_cell, zm_cell, ze_cell;
    ierr = GetOwnedCellRange(&info_nodes, 0, &xs_cell, &xm_cell); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_nodes, 1, &ys_cell, &ym_cell); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_nodes, 2, &zs_cell, &zm_cell); CHKERRQ(ierr);
    xe_cell = xs_cell + xm_cell; // Exclusive end of owned cells
    ye_cell = ys_cell + ym_cell;
    ze_cell = zs_cell + zm_cell;

    // Define ghosted cell ranges (inclusive start, exclusive end)
    PetscInt gxs_cell = xs_cell - stencil_width;
    PetscInt gxe_cell = xe_cell + stencil_width;
    PetscInt gys_cell = ys_cell - stencil_width;
    PetscInt gye_cell = ye_cell + stencil_width;
    PetscInt gzs_cell = zs_cell - stencil_width;
    PetscInt gze_cell = ze_cell + stencil_width;

    // Check if cell index (idx, idy, idz) is within the ghosted range
    // Ensure we also check against physical domain bounds (0 to IM-1 etc.)
    // Although typically ghost indices outside might be valid temporarily before migration.
    // Let's keep the check relative to the ghosted region accessible by this process.
    if (idx >= gxs_cell && idx < gxe_cell &&
        idy >= gys_cell && idy < gye_cell &&
        idz >= gzs_cell && idz < gze_cell) {
        *is_within = PETSC_TRUE; // It's within the area this proc can access (owned + ghost)
    }
    else {
        *is_within = PETSC_FALSE; // It's outside the accessible ghosted region
    }

    LOG_ALLOW(LOCAL,LOG_DEBUG, "Cell (%d, %d, %d) is %s the ghosted local grid (x:%d..%d, y:%d..%d, z:%d..%d).\n",
        idx, idy, idz, (*is_within) ? "within" : "outside",
        gxs_cell, gxe_cell, gys_cell, gye_cell, gzs_cell, gze_cell);

    return 0;
}
*/


/**
 * @brief Checks if the current CELL indices are within the LOCAL GHOSTED grid boundaries.
 *
 * This function determines if the provided global cell indices fall within the
 * range accessible by the current process, including its ghost cell layers.
 *
 * @param[in]  user       Pointer to the user-defined context containing grid information (needs fda for node info).
 * @param[in]  idx        The global i-index of the current cell.
 * @param[in]  idy        The global j-index of the current cell.
 * @param[in]  idz        The global k-index of the current cell.
 * @param[out] is_within  Pointer to a PetscBool that will be set to PETSC_TRUE if within ghosted bounds, else PETSC_FALSE.
 *
 * @return PetscErrorCode  Returns 0 on success, non-zero on failure.
 */


/**
 * @brief Updates the cell indices based on the signed distances to each face.
 *
 * This function modifies the cell indices (`idx`, `idy`, `idz`) to move towards the direction
 * where the particle is likely to be located, based on positive distances indicating
 * that the particle is outside in that particular direction.
 *
 * @param[in]  d    An array of six `PetscReal` values representing the signed distances to each face:
 *                  - d[LEFT]: Left Face
 *                  - d[RIGHT]: Right Face
 *                  - d[BOTTOM]: Bottom Face
 *                  - d[TOP]: Top Face
 *                  - d[FRONT]: Front Face
 *                  - d[BACK]: Back Face
 * @param[out] idx  Pointer to the i-index of the cell to be updated.
 * @param[out] idy  Pointer to the j-index of the cell to be updated.
 * @param[out] idz  Pointer to the k-index of the cell to be updated.
 * @param[in]  info DMDALocalInfo structure that holds local & global domain bounds.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode UpdateCellIndicesBasedOnDistances( PetscReal d[NUM_FACES], PetscInt *idx, PetscInt *idy, PetscInt *idz, DMDALocalInfo *info)
{
    PetscInt cxm,cxs;  // maximum & minimum cell ID in x
    PetscInt cym,cys;  // maximum & minimum cell ID in y
    PetscInt czm,czs;  // maximum & minimum cell ID in z

    cxs = info->xs; cxm = cxs + info->xm - 2;
    cys = info->ys; cym = cys + info->ym - 2;
    czs = info->zs; czm = czs + info->zm - 2; 

    LOG_ALLOW(LOCAL, LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Received d: "
	      "d[LEFT=%d]=%.3e, d[RIGHT=%d]=%.3e, d[BOTTOM=%d]=%.3e, "
	      "d[TOP=%d]=%.3e, d[FRONT=%d]=%.3e, d[BACK=%d]=%.3e\n",
	      LEFT, d[LEFT], RIGHT, d[RIGHT], BOTTOM, d[BOTTOM],
	      TOP, d[TOP], FRONT, d[FRONT], BACK, d[BACK]);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Raw d: "
	      "[%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]\n",
	      d[0], d[1], d[2], d[3], d[4], d[5]);
    
    // Validate input pointers
    if (d == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "UpdateCellIndicesBasedOnDistances - 'd' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UpdateCellIndicesBasedOnDistances - Input array 'd' is NULL.");
    }
    if (idx == NULL || idy == NULL || idz == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "UpdateCellIndicesBasedOnDistances - One or more index pointers are NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UpdateCellIndicesBasedOnDistances - One or more index pointers are NULL.");
    }

    // Debug: Print current face distances
    LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Current Face Distances:\n");
    if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)) LOG_FACE_DISTANCES(d);

    // Update k-direction based on FRONT and BACK distances
    if (d[FRONT] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Condition met: d[FRONT] < 0.0, incrementing idz.\n");
        (*idz) += 1;
    }
    else if(d[BACK] < 0.0){
      LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Condition met: d[BACK] < 0.0, decrementing idz.\n");
        (*idz) -= 1;
    }

    // Update i-direction based on LEFT and RIGHT distances
    if (d[LEFT] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Condition met: d[LEFT] < 0.0, decrementing idx.\n");
        (*idx) -= 1;
    }
    else if (d[RIGHT] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Condition met: d[RIGHT] < 0.0, incrementing idx.\n");
        (*idx) += 1;
    }

    // Update j-direction based on BOTTOM and TOP distances
    if (d[BOTTOM] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Condition met: d[BOTTOM] < 0.0, decrementing idy.\n");
        (*idy) -= 1;
    }
    else if (d[TOP] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Condition met: d[TOP] < 0.0, incrementing idy.\n");
        (*idy) += 1;
    }

    // The 'cell' corners you can reference go from [xs .. xs+xm-1], but
    // to form a valid cell in x, you need (idx+1) in range, so max is (xs+xm-2).
    *idx = PetscMax(cxs,               PetscMin(*idx, cxm));
    *idy = PetscMax(cys,               PetscMin(*idy, cym));
    *idz = PetscMax(czs,               PetscMin(*idz, czm));

    LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Updated Indices after clamping (inside domain bounds)  - idx, idy, idz: %d, %d, %d\n", *idx, *idy, *idz);

    return 0; // Indicate successful execution
}


/**
 * @brief Locates the grid cell containing a particle using a walking search.
 * @ingroup ParticleLocation
 *
 * This function implements a walking search algorithm to find the specific cell
 * (identified by global indices i, j, k) that encloses the particle's physical
 * location (`particle->loc`).
 *
 * The search starts from an initial guess cell (either the particle's previously known
 * cell or the corner of the local process domain) and iteratively steps to adjacent
 * cells based on the signed distances from the particle to the faces of the current
 * cell.
 *
 * Upon successful location (`position == 0`), the particle's `cell` field is updated
 * with the found indices (i,j,k), and the corresponding interpolation `weights`
 * are calculated and stored using the distances (`d`) relative to the final cell.
 *
 * Handles particles exactly on cell boundaries by attempting a tie-breaker if the
 * search gets stuck oscillating between adjacent cells due to the boundary condition.
 *
 * @param[in]  user     Pointer to the UserCtx structure containing grid information (DMDA, coordinates)
 *                      and domain boundaries.
 * @param[in,out] particle Pointer to the Particle structure. Its `loc` field provides the
 *                      target position. On successful return, its `cell` and `weights`
 *                      fields are updated. On failure, `cell` is set to `{-1, -1, -1}`
 *                      and `weights` to `{0.0, 0.0, 0.0}`.
 *
 * @return PetscErrorCode 0 on success. Non-zero error codes may indicate issues during
 *                        coordinate access, distance calculation, or other internal errors.
 *                        A return code of 0 does not guarantee the particle was found;
 *                        check `particle->cell[0] >= 0` afterward.
 *
 * @note Relies on helper functions like `InitializeTraversalParameters`, `CheckCellWithinLocalGrid`,
 *       `RetrieveCurrentCell`, `EvaluateParticlePosition`, `UpdateCellIndicesBasedOnDistances`,
 *       `UpdateParticleWeights`, and `FinalizeTraversal`.
 * @warning Ensure `particle->loc` is set correctly before calling.
 * @warning The function may fail to find the particle if it lies outside the domain accessible
 *          by the current process (including ghost cells) or if `MAX_TRAVERSAL` steps are exceeded.
 */
PetscErrorCode LocateParticleInGrid(UserCtx *user, Particle *particle)
{
    PetscErrorCode ierr;
    PetscInt idx, idy, idz;           // Current search cell indices
    PetscInt traversal_steps;         // Counter for search steps
    PetscBool cell_found = PETSC_FALSE;// Flag indicating success
    const Cmpnts p = particle->loc;   // Particle position (local copy)
    Cell current_cell;                // Geometry of the cell being checked
    const PetscReal threshold = DISTANCE_THRESHOLD; // Base threshold for distance checks
    DMDALocalInfo info;               // Local grid info (for bounds)
    PetscInt repeatedIndexCount = 0;  // Counter for consecutive iterations at the same index
    PetscInt prevIdx = PETSC_MIN_INT; // Previous iteration's index i (for repeat check)
    PetscInt prevIdy = PETSC_MIN_INT; // Previous iteration's index j (for repeat check)
    PetscInt prevIdz = PETSC_MIN_INT; // Previous iteration's index k (for repeat check)
    PetscInt last_position = -999;    // Position result (-1, 0, >=1) from *previous* iteration

    PetscFunctionBeginUser;
    LOG_FUNC_TIMER_BEGIN_EVENT(EVENT_Individualwalkingsearch, LOCAL);

    // Get local grid information (needed for bounds and index updates)
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

    // Determine starting cell indices (previous cell or local corner)
    ierr = InitializeTraversalParameters(user, particle, &idx, &idy, &idz, &traversal_steps); CHKERRQ(ierr);

    // --- Main Walking Search Loop ---
    while (!cell_found && traversal_steps < MAX_TRAVERSAL) {
        traversal_steps++;

        // Log current state for debugging
        LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, traversal_steps, 1,
                       "LocateParticleInGrid [PID %lld, Step %d]: Checking cell (%d, %d, %d). Pos (%.6f, %.6f, %.6f)\n",
                       (long long)particle->PID, traversal_steps, idx, idy, idz, p.x, p.y, p.z);

        // --- Check for Stuck Loop / Oscillation ---
        if (idx == prevIdx && idy == prevIdy && idz == prevIdz) {
            repeatedIndexCount++;
            LOG_ALLOW(LOCAL, LOG_DEBUG,"LocateParticleInGrid [PID %lld, Step %d]: Repeated index (%d,%d,%d) count = %d\n",
                        (long long)particle->PID, traversal_steps, idx, idy, idz, repeatedIndexCount);

            if (repeatedIndexCount > REPEAT_COUNT_THRESHOLD) {
                // --- Option B Tie-Breaker Logic ---
                if (last_position >= 1) { // Was the last step result 'on boundary'? Likely oscillation.
                    LOG_ALLOW(LOCAL, LOG_WARNING,
                              "LocateParticleInGrid [PID %lld]: Stuck at index (%d,%d,%d) after being on boundary for >%d steps. Applying boundary tie-break.\n",
                              (long long)particle->PID, idx, idy, idz, REPEAT_COUNT_THRESHOLD);

                    // Attempt to assign to the cell where we are stuck.
                    // Need to ensure the cell is valid and re-get distances for weights.
                    Cell final_cell;
                    PetscReal final_d[NUM_FACES];
                    PetscInt final_position; // Result not used, just need distances
                    PetscBool is_stuck_cell_valid;

                    ierr = CheckCellWithinLocalGrid(user, idx, idy, idz, &is_stuck_cell_valid); CHKERRQ(ierr);
                    if (is_stuck_cell_valid) {
                         ierr = RetrieveCurrentCell(user, idx, idy, idz, &final_cell); CHKERRQ(ierr);
                         // Re-evaluate to get distances 'final_d' relative to this specific cell
                         ierr = EvaluateParticlePosition(&final_cell, final_d, p, &final_position, threshold);
                         if (ierr == 0) { // Check if evaluation succeeded
                            ierr = UpdateParticleWeights(final_d, particle); CHKERRQ(ierr); // Calculate weights
                            particle->cell[0] = idx; // Assign the stuck cell index
                            particle->cell[1] = idy;
                            particle->cell[2] = idz;
                            cell_found = PETSC_TRUE; // Mark as found via tie-break
                            LOG_ALLOW(LOCAL, LOG_INFO,
                                      "LocateParticleInGrid [PID %lld]: Assigned via tie-break to cell (%d,%d,%d).\n",
                                      (long long)particle->PID, idx, idy, idz);
                         } else {
                             // Error during final evaluation
                             LOG_ALLOW(LOCAL, LOG_ERROR,
                                       "LocateParticleInGrid [PID %lld]: Error %d during final evaluation for tie-break at cell (%d,%d,%d). Search fails.\n",
                                       (long long)particle->PID, ierr, idx, idy, idz);
                             cell_found = PETSC_FALSE; // Ensure failure
                         }
                    } else {
                         LOG_ALLOW(LOCAL, LOG_WARNING,
                                   "LocateParticleInGrid [PID %lld]: Stuck at invalid cell index (%d,%d,%d) during tie-break attempt. Search fails.\n",
                                   (long long)particle->PID, idx, idy, idz);
                         cell_found = PETSC_FALSE; // Mark as failed
                    }
                    break; // Exit loop after attempting tie-break
                } else {
                    // Stuck for > threshold steps, but last known state wasn't 'on boundary'.
                    LOG_ALLOW(LOCAL, LOG_WARNING,
                              "LocateParticleInGrid [PID %lld]: Stuck at index (%d,%d,%d) for >%d steps (last_pos=%d). Breaking search (failed).\n",
                              (long long)particle->PID, idx, idy, idz, REPEAT_COUNT_THRESHOLD, last_position);
                    cell_found = PETSC_FALSE; // Mark as failure
                    break; // Exit loop (failed)
                }
            } // end if (repeatedIndexCount > threshold)
        } else { // Index changed from previous iteration
            repeatedIndexCount = 0;
        }
        // --- End Stuck Loop Check ---

        // --- Update previous index state for the *next* iteration's check ---
        prevIdx = idx;
        prevIdy = idy;
        prevIdz = idz;

        // --- Check if current cell index is within accessible grid bounds ---
        PetscBool is_within;
        ierr = CheckCellWithinLocalGrid(user, idx, idy, idz, &is_within); CHKERRQ(ierr);
        if (!is_within) {
            LOG_ALLOW(LOCAL, LOG_WARNING,
                      "LocateParticleInGrid [PID %lld]: Search moved outside local ghosted grid boundaries at cell (%d, %d, %d). Breaking search.\n",
                      (long long)particle->PID, idx, idy, idz);
            cell_found = PETSC_FALSE; // Mark as failed
            break; // Exit loop
        }

        // --- Retrieve geometry and evaluate particle position ---
        ierr = RetrieveCurrentCell(user, idx, idy, idz, &current_cell); CHKERRQ(ierr);
        PetscReal current_d[NUM_FACES]; // Holds distances relative to current_cell
        PetscInt position;              // Holds result: -1 (outside), 0 (inside), >=1 (boundary)
        ierr = EvaluateParticlePosition(&current_cell, current_d, p, &position, threshold); CHKERRQ(ierr);

        // Store the position result for the *next* iteration's stuck check
        last_position = position;

        LOG_ALLOW(LOCAL, LOG_DEBUG,
                   "LocateParticleInGrid [PID %lld, Step %d]: Evaluated pos in cell (%d,%d,%d): position=%d\n",
                   (long long)particle->PID, traversal_steps, idx, idy, idz, position);

        // --- Main Logic: Handle position result ---
        if (position == 0) { // Strictly INSIDE the cell
            cell_found = PETSC_TRUE;
            particle->cell[0] = idx; // Assign the indices of the cell we just evaluated
            particle->cell[1] = idy;
            particle->cell[2] = idz;
            // Calculate weights using the distances 'current_d' from this successful evaluation
            ierr = UpdateParticleWeights(current_d, particle); CHKERRQ(ierr);
            LOG_ALLOW(LOCAL, LOG_INFO,
                      "LocateParticleInGrid [PID %lld]: Found INSIDE cell (%d, %d, %d) after %d steps. Weights calculated.\n",
                      (long long)particle->PID, idx, idy, idz, traversal_steps);
            break; // Found the correct cell, exit loop

        } else { // Particle is OUTSIDE (position == -1) or ON BOUNDARY (position >= 1)
            LOG_ALLOW(LOCAL, LOG_DEBUG,
                      "LocateParticleInGrid [PID %lld, Step %d]: Particle OUTSIDE or ON BOUNDARY of cell (%d, %d, %d), position = %d. Updating indices.\n",
                      (long long)particle->PID, traversal_steps, idx, idy, idz, position);
            // Both cases require stepping to the next potential cell.
            // Update search indices (idx, idy, idz) based on the distances 'current_d'
            ierr = UpdateCellIndicesBasedOnDistances(current_d, &idx, &idy, &idz, &info); CHKERRQ(ierr);
            // Loop continues to the next iteration with potentially updated idx, idy, idz
        }
    } // --- End while loop ---

    // --- Handle cases where the loop terminated WITHOUT finding the cell ---
    // (Could be MAX_TRAVERSAL, out of bounds, or non-boundary stuck state)
    if (!cell_found) {
        LOG_ALLOW(LOCAL, LOG_WARNING,
                  "LocateParticleInGrid [PID %lld]: Search FAILED after %d steps. Final check was at cell (%d,%d,%d). Resetting cell/weights.\n",
                  (long long)particle->PID, traversal_steps, idx, idy, idz);
        // Ensure particle state reflects failure explicitly
        particle->cell[0] = -1;
        particle->cell[1] = -1;
        particle->cell[2] = -1;
    }

    // Finalize traversal by reporting the results (logs success/failure message)
    // Note: idx, idy, idz passed here are the *last checked* indices for logging purposes.
    // The actual found cell index is stored in particle->cell if cell_found is true.
    ierr = FinalizeTraversal(user, particle, traversal_steps, cell_found, idx, idy, idz); CHKERRQ(ierr);

    LOG_FUNC_TIMER_END_EVENT(EVENT_Individualwalkingsearch, LOCAL);
    PetscFunctionReturn(0);
}
//////////////////////////


/**
 * @brief Locates all particles within the grid and calculates their interpolation weights.
 * @ingroup ParticleLocation
 *
 * This function iterates through all particles currently local to this MPI rank.
 * For each particle, it first checks if the particle is within the rank's
 * pre-calculated bounding box (`user->bbox`). If it is, it calls the
 * `LocateParticleInGrid` function to perform the walking search.
 *
 * `LocateParticleInGrid` is responsible for finding the containing cell `(i,j,k)`
 * and calculating the corresponding interpolation weights `(w1,w2,w3)`. It updates
 * the `particle->cell` and `particle->weights` fields directly upon success.
 * If the search fails (particle not found within MAX_TRAVERSAL, goes out of bounds,
 * or gets stuck without resolution), `LocateParticleInGrid` sets the particle's
 * `cell` to `{-1,-1,-1}` and `weights` to `{0.0, 0.0, 0.0}`.
 *
 * After attempting location, this function updates the corresponding entries in the
 * DMSwarm's "DMSwarm_CellID" and "weight" fields using the potentially modified
 * data from the `particle` struct.
 *
 * @param[in] user Pointer to the UserCtx structure containing grid, swarm, and bounding box info.
 *
 * @return PetscErrorCode Returns `0` on success, non-zero on failure (e.g., errors accessing DMSwarm fields).
 *
 * @note Assumes `user->bbox` is correctly initialized for the local rank.
 * @note Assumes `InitializeParticle` correctly populates the temporary `particle` struct.
 * @note Assumes `UpdateSwarmFields` correctly writes data back to the DMSwarm.
 */
PetscErrorCode LocateAllParticlesInGrid(UserCtx *user) {
    PetscErrorCode ierr;
    PetscMPIInt rank, size;
    PetscInt localNumParticles;
    PetscReal *positions = NULL, *weights = NULL, *velocities = NULL; // Pointers to DMSwarm data arrays
    PetscInt *cellIndices = NULL;
    PetscInt *LocStatus = NULL;
    PetscInt64 *PIDs = NULL;      // Pointers to DMSwarm data arrays
    DM swarm = user->swarm;                 // Convenience pointer to the swarm DM
    Particle particle;                      // Reusable temporary Particle struct for processing

    PetscFunctionBeginUser;
    LOG_FUNC_TIMER_BEGIN_EVENT(EVENT_walkingsearch, LOCAL);

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "LocateAllParticlesInGrid - Start on Rank %d/%d.\n", rank, size);

    

    // Optional barrier for debugging synchronization
    // ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

    // --- Access DMSwarm Data Arrays ---
    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "LocateAllParticlesInGrid - Number of local particles: %d.\n", localNumParticles);

    // Get direct pointers to the underlying data arrays for efficiency
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr); // Array to write weights back to
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIndices); CHKERRQ(ierr); // Array to write cell indices back to
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&PIDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_location_status", NULL, NULL, (void**)&LocStatus); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "LocateAllParticlesInGrid - DMSwarm fields accessed successfully.\n");

    //---  TEST
    /*
    PetscInt64 *pid_snapshot = NULL;
    ierr = GetLocalPIDSnapshot(PIDs,localNumParticles,&pid_snapshot);

    for(PetscInt p = 0; p < localNumParticles;++p){
      LOG_LOOP_ALLOW(LOCAL,LOG_DEBUG,p,10," S.No = %d | Snapshot PID = %ld | PID = %ld .\n",p,pid_snapshot[p],PIDs[p]); 
    }

    PetscFree(pid_snapshot);
    */
    // ------------
    
    // --- Iterate over each local particle ---
    for (PetscInt i = 0; i < localNumParticles; ++i) {
        // Load current particle data into the temporary struct
      ierr = UnpackSwarmFields(i,PIDs, weights, positions, cellIndices, velocities, LocStatus, &particle); CHKERRQ(ierr);

        LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, i, 10, "LocateAllParticlesInGrid - Processing Particle [%d]: PID=%ld.\n", i, particle.PID);

        // --- Coarse Check: Is particle within this rank's bounding box? ---
        // This is a quick check; particle could still be in a ghost cell managed by this rank.
        PetscBool particle_detected = IsParticleInsideBoundingBox(&(user->bbox), &particle);
        LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, i, 10, "LocateAllParticlesInGrid - Particle [%d] (PID %ld) inside local bbox: %s.\n",
                       i, particle.PID, particle_detected ? "YES" : "NO");

        if (particle_detected) {
            // --- Perform Detailed Location Search ---
            LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, i, 10, "LocateAllParticlesInGrid - Locating Particle [%d] (PID %ld) in grid...\n", i, particle.PID);

            // Call the walking search. This function will update particle.cell and particle.weights
            // internally if successful, or set them to -1 / 0 if it fails.
            ierr = LocateParticleInGrid(user, &particle); CHKERRQ(ierr); // Pass only user and particle struct

            // Log the outcome of the search for this particle
            if (particle.cell[0] >= 0) {
                 LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, i, 10,
                                "LocateAllParticlesInGrid - Particle [%d] (PID %ld) located/assigned to cell [%d, %d, %d].\n",
                                i, particle.PID, particle.cell[0], particle.cell[1], particle.cell[2]);
            } else {
                 LOG_LOOP_ALLOW(LOCAL, LOG_WARNING, i, 1, // Log all failures
                                "LocateAllParticlesInGrid - Particle [%d] (PID %ld) FAILED TO LOCATE (CellID = -1).\n",
                                i, particle.PID);
            }
            // --- Weight calculation is now handled inside LocateParticleInGrid ---

        } else {
            // Particle was outside the local bounding box - mark as invalid for this rank
            LOG_LOOP_ALLOW(LOCAL, LOG_WARNING, i, 1,
                           "LocateAllParticlesInGrid - Particle [%d] (PID %ld) outside local bbox. Marking invalid (CellID = -1).\n",
                           i, particle.PID);
            particle.cell[0] = -1;
            particle.cell[1] = -1;
            particle.cell[2] = -1;
        } // end if (particle_detected)

        // --- Update DMSwarm Data ---
        // Write the potentially modified cell index and weights from the 'particle' struct
        // back into the main DMSwarm data arrays.
        ierr = UpdateSwarmFields(i, &particle, weights, cellIndices, LocStatus); CHKERRQ(ierr);

    } // --- End particle loop ---

    // --- Restore DMSwarm Data Arrays ---
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIndices); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&PIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_location_status", NULL, NULL, (void**)&LocStatus); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(LOCAL, LOG_INFO, "LocateAllParticlesInGrid - DMSwarm fields restored successfully on Rank %d.\n", rank);

    LOG_ALLOW(LOCAL, LOG_DEBUG, "LocateAllParticlesInGrid - Completed function on Rank %d.\n", rank);
    LOG_FUNC_TIMER_END_EVENT(EVENT_walkingsearch, LOCAL);
    PetscFunctionReturn(0);
}



/**
 * @brief Performs the complete initial setup for the particle simulation at time t=0.
 *
 * This includes:
 * 1. Initial locating of particles (based on their potentially arbitrary initial assignment).
 * 2. A preliminary migration cycle to ensure particles are on the MPI rank that owns
 *    their initial physical region.
 * 3. If `user->ParticleInitialization == 0` (Surface Init), re-initializes particles on the
 *    designated inlet surface. This ensures particles migrated to an inlet-owning rank
 *    are correctly distributed on that surface.
 * 4. A final locating of all particles to get their correct cell indices and interpolation weights.
 * 5. Interpolation of initial Eulerian fields to the particles.
 * 6. Scattering of particle data to Eulerian fields (if applicable).
 * 7. Outputting initial data if requested.
 *
 * @param user Pointer to the UserCtx structure.
 * @param currentTime The current simulation time (should be StartTime, typically 0.0).
 * @param step The current simulation step (should be StartStep, typically 0).
 * @param readFields Flag indicating if Eulerian fields were read from file (influences output).
 * @param bboxlist Array of BoundingBox structures for domain decomposition.
 * @param OutputFreq Frequency for writing output files.
 * @param StepsToRun Total number of simulation steps planned (used for output logic on setup-only runs).
 * @param StartStep The starting step of the simulation (used for output logic).
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformInitialSetup(UserCtx *user, PetscReal currentTime, PetscInt step,
                                   PetscBool readFields, const BoundingBox *bboxlist,
                                   PetscInt OutputFreq, PetscInt StepsToRun, PetscInt StartStep)
{
    PetscErrorCode ierr;
    MigrationInfo  *migrationList_for_initial_sort = NULL; // Managed locally for this initial sort
    PetscInt       migrationCount_initial_sort = 0;
    PetscInt       migrationListCapacity_initial_sort = 0;
    PetscInt       globalMigrationCount_initial_sort;
    PetscMPIInt    rank; // For logging
    PetscReal      uMax = 0.0, uMaxcat = 0.0;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Performing initial particle setup procedures.\n", currentTime, step);

    // --- 1. Initial Locate (based on positions from InitializeParticleBasicProperties) ---
    // This gets a first guess of where particles are, which might be (0,0,0) for many if not placed by their initial rank.
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Initial Locate: Determining particle cells before preliminary migration.\n", currentTime, step);
    ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL,LOG_INFO," Particle layout (pre-preliminary migration): \n");
    ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);

    // --- 2. Preliminary Spatial Sorting Migration ---
    // This moves particles (e.g., those at (0,0,0)) to the rank that actually owns their physical region.
    ierr = PerformSingleParticleMigrationCycle(user, bboxlist,
                                               &migrationList_for_initial_sort, &migrationCount_initial_sort, &migrationListCapacity_initial_sort,
                                               currentTime, step, "Preliminary Sort", &globalMigrationCount_initial_sort); CHKERRQ(ierr);
    // migrationList_for_initial_sort is allocated/reallocated within PerformSingleParticleMigrationCycle if needed.
    
    // --- 3. Re-initialize Particles on Inlet Surface (if applicable) ---
    // This is crucial for ParticleInitialization == 0 (Surface Init). After particles have been
    // migrated to the rank(s) owning the inlet surface, this step ensures they are properly
    // distributed *on* that surface, rather than remaining at their migrated position (e.g., (0,0,0)).
 
    if (user->ParticleInitialization == 0 && user->inletFaceDefined) {
      ierr = ReinitializeParticlesOnInletSurface(user, currentTime, step); CHKERRQ(ierr);
      // ---   DEBUG
      // if (rank == 0) { // Or the rank that does the re-init
      //  PetscReal *coords_check;
      //  PetscInt nlocal_check;
      //  ierr = DMSwarmGetLocalSize(user->swarm, &nlocal_check); CHKERRQ(ierr);
      //  if (nlocal_check > 0) {
      //     ierr = DMSwarmGetField(user->swarm, "position", NULL, NULL, (void**)&coords_check); CHKERRQ(ierr);
      //	  LOG_ALLOW(LOCAL,LOG_DEBUG, "[Rank %d DEBUG] After Reinit, first particle pos: (%.2f, %.2f, %.2f)\n",
      //	      rank, coords_check[0], coords_check[1], coords_check[2]);
	  // Check for NaNs/Infs in a few particles if suspicious
      //  for (int p_chk = 0; p_chk < PetscMin(5, nlocal_check); ++p_chk) {
      //	    if (PetscIsInfOrNanReal(coords_check[3*p_chk+0]) ||
      //	PetscIsInfOrNanReal(coords_check[3*p_chk+1]) ||
      //	PetscIsInfOrNanReal(coords_check[3*p_chk+2])) {
      //      LOG_ALLOW(LOCAL,LOG_DEBUG, "[Rank %d DEBUG ERROR] Bad coord for particle %d after reinit!\n", rank, p_chk);
      //    }
      //  }
      //  ierr = DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void**)&coords_check); CHKERRQ(ierr);
      //  }
      //  fflush(stdout);
	// }
    // --------- DEBUG	
    }

    LOG_ALLOW(LOCAL,LOG_DEBUG," [ Rank %d] ReinitializeParticlesOnInletSurface completed, heading into PetscBarrier.\n",rank);
    
    // Ensure all ranks wait till Reinitialization is done to proceed.
    ierr = PetscBarrier(NULL);
    
    // --- 4. Final Locate and Interpolate ---
    // Now that particles are on their correct ranks and (if surface init) correctly positioned on the inlet,
    // perform the definitive location and interpolation for t=0.
    LOG_ALLOW(LOCAL, LOG_INFO, "[T=%.4f, Step=%d , Rank=%d] Second Initial Locate & Interpolate (post-migration & re-placement).\n", currentTime, step,rank);
    ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);
    ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr);

    VecNorm(user->Ucont, NORM_INFINITY, &uMax);
    LOG_ALLOW(GLOBAL, LOG_INFO,"[Step %d] Initial  max(|Ucont|) before Scatter  = %.6e\n",step , uMax);
    uMax = 0.0;

    VecNorm(user->Ucat, NORM_INFINITY, &uMaxcat);
    LOG_ALLOW(GLOBAL, LOG_INFO,"[Step %d]  Initial max(|Ucat|) before Scatter  = %.6e\n", step, uMaxcat);
    uMaxcat = 0.0;
    
    // --- 5. Scatter Particle Data to Eulerian Grid (if part of initialization) ---
    ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);

    VecNorm(user->Ucont, NORM_INFINITY, &uMax);
    LOG_ALLOW(GLOBAL, LOG_INFO,"[Step %d] Initial  max(|Ucont|) after Scatter  = %.6e\n",step , uMax);
    uMax = 0.0;

    VecNorm(user->Ucat, NORM_INFINITY, &uMaxcat);
    LOG_ALLOW(GLOBAL, LOG_INFO,"[Step %d]  Initial max(|Ucat|) after Scatter  = %.6e\n", step, uMaxcat);
    uMaxcat = 0.0;
    
    // --- 6. Initial Output ---
    // Output if OutputFreq > 0 OR if it's a setup-only run (StepsToRun==0 and StartStep==0).
    if (OutputFreq > 0 || (StepsToRun == 0 && StartStep == 0)) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Logging initial data and (if applicable) interpolation error.\n", currentTime, step);
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Printing Initial particle fields at t=%.4f (step %d completed)\n", currentTime, step);
       
	//  LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] ABOUT TO CALL LOG_PARTICLE_FIELDS.\n", rank);
	ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);
	//  LOG_ALLOW(LOCAL, LOG_INFO, "[Rank %d] RETURNED FROM LOG_PARTICLE_FIELDS.\n", rank); 

	// DEPRECIATED - Testing against sinusoidal (now FieldInitialization = 1 has changed meaning.
	// if(user->FieldInitialization==1) { // If analytic fields are available for comparison
	//    ierr = LOG_INTERPOLATION_ERROR(user); CHKERRQ(ierr);
        //}

        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Writing initial particle data (for step %d).\n", currentTime, step, step);
	//	LOG_ALLOW(LOCAL, LOG_INFO, "[%d] About to call WriteSwarmField for position.\n", rank);
        ierr = WriteSwarmField(user, "position", step, "dat"); CHKERRQ(ierr);
        ierr = WriteSwarmField(user, "velocity", step, "dat"); CHKERRQ(ierr);
        // ierr = WriteSwarmField(user, "pos_phy",  step, "dat"); CHKERRQ(ierr); // If needed
        if (!readFields) { // Only write grid fields if they were generated, not read
            ierr = WriteSimulationFields(user); CHKERRQ(ierr);
        }
    }

    // --- 7. Cleanup Memory ---
    // Free the migration list used specifically for this initial sort.
    ierr = PetscFree(migrationList_for_initial_sort); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


    /*
    // ==============================================================================
    // --- STEP 1: Update the INTERIOR of the domain ---
    // ==============================================================================

    if (step == StartStep && readFields) {
        // --- RESTART from file ---
        // This case reads the full, previously saved state, overwriting everything.
        LOG_ALLOW(GLOBAL, LOG_INFO, "RESTART condition: Reading all fields from file for step %d.\n", step);
        ierr = ReadSimulationFields(user, step); CHKERRQ(ierr);

	// Even after reading, we must update local ghosts to ensure consistency for subsequent steps.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghost regions for all fields after reading.\n", time, step);
        ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
	
    } else {
        // --- This block handles both initial setup and time-advancement ---
        
          if (step == StartStep) {
            // --- Initial Field Setup ---
            // The boundaries have already been initialized by BoundarySystem_Create.
            // We now fill the domain interior using the newly named function.
	    // ierr = VecZeroEntries(user->Ucont);CHKERRQ(ierr);
	    // ----DEBUG TEST
	    // ierr = VecSet(user->Ucont,1.0);
	    LOG_ALLOW(GLOBAL, LOG_INFO, "INITIAL start: Generating INTERIOR fields for initial step %d.\n", step);
            ierr = SetInitialInteriorField(user, "Ucont"); CHKERRQ(ierr);

	    ierr = VecNorm(user->Ucont, NORM_INFINITY, &ucont_max); CHKERRQ(ierr);
	    LOG_ALLOW(GLOBAL,LOG_INFO,"[DEBUG] max(|Ucont|) after SetInitialInteriorField = %.6e\n", ucont_max);
	    ucont_max=0.0;
	    
	         } else { // Advancing the simulation (step > StartStep)
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "ADVANCE step: Updating INTERIOR fields for step %d.\n", step);


	    ierr = VecNorm(user->Ucont, NORM_INFINITY, &ucont_max); CHKERRQ(ierr);
	    LOG_ALLOW(GLOBAL,LOG_INFO,"[DEBUG] max(|Ucont|) (timestep = %d)  = %.6e\n",step,ucont_max);
	    ucont_max=0.0;
	    // This is the hook for the actual fluid dynamics solver, which would update Ucont.
            // ierr = YourNavierStokesSolver(user, user->dt); CHKERRQ(ierr);
	        }

        // ==============================================================================
        // --- STEP 2: APPLY BOUNDARY CONDITIONS ---
        // ==============================================================================
        // The boundary system applies conditions. For a wall, this sets the normal
        // component of Ucont to 0 and also sets the ghost-cell Ucat values.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Executing boundary condition system.\n", time, step);
	ierr = BoundarySystem_ExecuteStep(user); CHKERRQ(ierr);
	
	// ==============================================================================
        // --- STEP 3: SYNCHRONIZE Ucont BEFORE CONVERSION ---
        // ==============================================================================
        // THIS IS THE CRITICAL FIX: Update lUcont with correct global data and ghost cells
        // BEFORE calling Contra2Cart, which relies on lUcont for its stencil.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghost regions for Ucont.\n", time, step);
        ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);

        // ==============================================================================
        // --- STEP 4: CONVERT CONTRAVARIANT TO CARTESIAN ---
        // ==============================================================================
        // With Ucont fully defined (interior + boundaries), compute the Cartesian velocity.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Converting Ucont to Ucat.\n", time, step);
	//	ierr = VecZeroEntries(user->Ucat); CHKERRQ(ierr);
	// DEBUG--------

        ierr = Contra2Cart(user); CHKERRQ(ierr);
	
	ierr = VecNorm(user->Ucat, NORM_INFINITY, &umax); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_INFO,"[DEBUG] max(|Ucat|) after Contra2Cart = %.6e\n", umax);
	umax = 0.0;

        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Executing boundary condition system.\n", time, step);
	ierr = BoundarySystem_ExecuteStep(user); CHKERRQ(ierr);	
	
	ierr = VecMin(user->Ucat, NULL, &umin);
	LOG_ALLOW(GLOBAL,LOG_DEBUG, "[DEBUG] min(Ucat) after Contra2Cart = %.6e\n", umin);
	// ==============================================================================
        // --- STEP 5: SYNCHRONIZE Ucat AFTER CONVERSION ---
        // ==============================================================================
        // Finally, update the local lUcat with the newly computed global Ucat data.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghost regions for Ucat.\n", time, step);
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
    }
    
    // ==============================================================================
    // --- STEP 4: UPDATE GHOST CELLS (FINALIZATION) ---
    // ==============================================================================
    // This final step synchronizes the ghost cell layers across all MPI ranks for
    // all relevant fields, ensuring a consistent state before they are used.
    //  LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Updating local ghost regions for all fields.\n", time, step);
    //  ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);
    // ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Complete Eulerian state is now finalized and consistent.\n", time, step);
    PetscFunctionReturn(0);
}
    */

/**
 * @brief Executes the main time-marching loop for the particle simulation.
 *
 * This function performs the following steps repeatedly for `StepsToRun`:
 * 1. Updates/Sets the background fluid velocity field (Ucat) for the current step.
 * 2. If it's the very first step (StartStep=0), calls `PerformInitialSetup` for
 *    preliminary migration, location, interpolation, and output.
 * 3. Updates particle positions using velocity from the *previous* step's interpolation.
 * 4. Performs boundary checks and removes out-of-bounds particles.
 * 5. Migrates particles between MPI processors using `PerformSingleParticleMigrationCycle`.
 * 6. Advances simulation time and step counter.
 * 7. Locates particles in the grid based on their *new* positions.
 * 8. Interpolates the fluid velocity (from the *current* Ucat) to the new particle locations
 *    to get velocities for the *next* advection step.
 * 9. Scatters particle data back to Eulerian fields.
 *10. Logs errors and outputs data at specified `OutputFreq` intervals.
 *
 * @param user         Pointer to the UserCtx structure.
 * @param StartStep    The initial step number (e.g., 0 for a new run, >0 for restart).
 * @param StartTime    The simulation time corresponding to StartStep.
 * @param StepsToRun   The number of steps to execute in this run. If 0 and StartStep is 0,
 *                     only `PerformInitialSetup` is executed.
 * @param OutputFreq   Frequency (in number of steps) at which to output data and log errors.
 * @param readFields   Flag indicating whether to read initial fields (used by `SetEulerianFields`
 *                     and `PerformInitialSetup`).
 * @param bboxlist     Array of BoundingBox structures for domain decomposition, used for migration.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode AdvanceSimulation(UserCtx *user, PetscInt StartStep, PetscReal StartTime,
                                 PetscInt StepsToRun, PetscInt OutputFreq, PetscBool readFields,
                                 const BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;            // MPI rank of the current process
    PetscReal      dt = user->dt;   // Timestep size
    PetscInt       step_loop_counter; // Loop counter for simulation steps
    PetscReal      currentTime;     // Current simulation time
    PetscInt       removed_local, removed_global; // Counters for out-of-bounds particles
    PetscReal      uMax = 0.0, uMaxcat = 0.0;

    // Variables for particle migration within the main loop
    MigrationInfo  *migrationList_main_loop = NULL;    // Managed by AdvanceSimulation for the loop
    PetscInt       migrationCount_main_loop = 0;
    PetscInt       migrationListCapacity_main_loop = 0;
    PetscInt       globalMigrationCount_in_loop_cycle; // Output from migration cycle

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting simulation run: %d steps from step %d (t=%.4f), dt=%.4f\n",
              StepsToRun, StartStep, StartTime, dt);

    currentTime = StartTime; // Initialize current simulation time

    // --- Handle Initial Setup (if StartStep is 0) ---
    if (StartStep == 0) {
        // Set Eulerian fields for t=0 before particle setup
        ierr = SetEulerianFields(user, StartStep, StartStep, currentTime); CHKERRQ(ierr);
        // Perform comprehensive initial particle setup
        ierr = PerformInitialSetup(user, currentTime, StartStep, readFields, bboxlist, OutputFreq, StepsToRun, StartStep); CHKERRQ(ierr);

        // If only initial setup was requested (StepsToRun == 0), exit now.
        if (StepsToRun == 0) {
            LOG_ALLOW(GLOBAL, LOG_INFO, "Initial setup completed as StepsToRun is 0. Time t=%.4f. Exiting AdvanceSimulation.\n", currentTime);
            // The migrationList_main_loop is not yet used, but free it for consistency.
            ierr = PetscFree(migrationList_main_loop); CHKERRQ(ierr);
            PetscFunctionReturn(0);
        }
    }
    
    // --- Time Marching Loop ---
    // Loops from the StartStep up to (but not including) StartStep + StepsToRun
    for (step_loop_counter = StartStep; step_loop_counter < StartStep + StepsToRun; step_loop_counter++)
    {
      LOG_ALLOW(GLOBAL, LOG_INFO, "Starting step %d, t=%.4f\n", step_loop_counter, currentTime);

      // Set/Update Eulerian Fields for the current time step
      // If StartStep was 0, fields for t=0 were already set.
      // If restarting (StartStep > 0), or for subsequent steps (step_loop_counter > StartStep),
      // set/update fields for the current `currentTime`.
      //    if (step_loop_counter > StartStep || (step_loop_counter == StartStep && StartStep > 0)) {
      //   ierr = SetEulerianFields(user, step_loop_counter, StartStep, currentTime, readFields); CHKERRQ(ierr);
	  //	    }

      // Step 3: Update Particle Positions
      // Moves particles from P(currentTime) to P(currentTime + dt) using velocity_particle(currentTime).
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f -> T=%.4f, Step=%d] Updating particle positions (Rank: %d).\n", currentTime, currentTime+dt, step_loop_counter, rank);
      ierr = UpdateAllParticlePositions(user); CHKERRQ(ierr);

      // Step 4: Check Boundaries and Remove Out-Of-Bounds Particles
      ierr = CheckAndRemoveOutOfBoundsParticles(user, &removed_local, &removed_global, bboxlist); CHKERRQ(ierr);
      if (removed_global > 0) {
          LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Removed %d out-of-bounds particles globally.\n", currentTime, step_loop_counter, removed_global);
      }

      LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG,"[Rank %d] ENTERING MIGRATION AFTER UPDATE.\n",rank);
      
      // Step 5 & 6: Migrate Particles between processors
      // Migration is based on new positions P(currentTime + dt).
      // The time logged for migration is the target time of the new positions.
      ierr = PerformSingleParticleMigrationCycle(user, bboxlist,
                                                 &migrationList_main_loop, &migrationCount_main_loop, &migrationListCapacity_main_loop,
                                                 currentTime + dt, step_loop_counter, "Main Loop", &globalMigrationCount_in_loop_cycle); CHKERRQ(ierr);
      // migrationCount_main_loop is reset inside PerformSingleParticleMigrationCycle for the *next* call within this loop.

      // Step 7: Advance Time
      currentTime += dt; // currentTime now represents the time at the *end* of the step just completed.

      // ierr = SetEulerianFields(user, step_loop_counter+1, StartStep, currentTime, readFields); CHKERRQ(ierr);

      // Step 9: Locate Particles at New Positions (t = currentTime)
      // This is for particles now on the current rank after migration, using their P(currentTime) positions.
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d completed] Locating particles at new positions (Rank: %d).\n", currentTime, step_loop_counter, rank);
      ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);
      // After migration,since the cellID is from previous rank, first locate fails(sets CellIDs to -1)
      // Second locate will successfully locate the particles in new rank.
      ierr = LocateAllParticlesInGrid(user); CHKERRQ(ierr);

      // Step 10: Interpolate Eulerian Field to New Particle Positions
      // Eulerian field Ucat is typically from the beginning of the current step (time effectively currentTime - dt).
      // Interpolate this Ucat to P(currentTime) to get V_particle(currentTime) for the *next* advection.
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d completed] Interpolating field (from Ucat at t=%.4f) to new particle positions (Rank: %d).\n", currentTime, step_loop_counter, currentTime-dt, rank);
      ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr);
      
      // Step 11: Scatter Particle Fields back to Eulerian Grid (if applicable)
      ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d completed] Scattering particle fields to grid (Rank: %d).\n", currentTime, step_loop_counter, rank);

      // Step 8: Update global step counter in user context.
      // This should reflect the step number *completed*.
      user->step = step_loop_counter + 1;
      
      // Step 12: Output and Error Logging (based on the step *just completed*, using user->step)
      if (OutputFreq > 0 && (user->step % OutputFreq == 0) ) {
          LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d completed (Output for step %d)] Logging interpolation error.\n", currentTime, step_loop_counter, user->step);

	// DEPRECIATED - Testing against sinusoidal (now FieldInitialization = 1 has changed meaning.
	//  if(user->FieldInitialization==1) { // If analytic fields are available for comparison
        //      ierr = LOG_INTERPOLATION_ERROR(user); CHKERRQ(ierr);
        //  }

	  
          LOG_ALLOW(GLOBAL, LOG_DEBUG, "Printing particle fields at t=%.4f (step %d completed, output for step %d)\n", currentTime, step_loop_counter, user->step);
          ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);

          LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d completed] Writing particle data (for step %d).\n", currentTime, step_loop_counter, user->step);
          ierr = WriteSwarmField(user, "position", user->step, "dat"); CHKERRQ(ierr); // Write P(t+dt)
          ierr = WriteSwarmField(user, "velocity", user->step, "dat"); CHKERRQ(ierr); // Write V_interp @ P(t+dt)
          // ierr = WriteSwarmField(user, "pos_phy",  user->step, "dat"); CHKERRQ(ierr);
          ierr = WriteSimulationFields(user); CHKERRQ(ierr); // Optional grid field write
      }
      LOG_ALLOW(GLOBAL, LOG_INFO, "Finished step %d, current time t=%.4f (user->step is now %d)\n", step_loop_counter, currentTime, user->step);
    } // End of time marching loop

    PetscReal finalTimeRun = StartTime + StepsToRun * dt; // Target final time if all steps ran
    LOG_ALLOW(GLOBAL, LOG_INFO, "Time marching completed. %d steps run from StartStep %d. Current time t=%.4f (target final: %.4f).\n", StepsToRun, StartStep, currentTime, finalTimeRun);

    // --- Final Particle Field Log (if not already done by OutputFreq and if steps were run) ---
    if (StepsToRun > 0 && !(OutputFreq > 0 && (user->step % OutputFreq == 0))) {
        // Avoid double logging if the last step was an output step.
        // Only log if actual simulation steps were performed.
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Printing final particle fields at end of run (t=%.4f):\n", currentTime);
        ierr = LOG_PARTICLE_FIELDS(user, user->LoggingFrequency); CHKERRQ(ierr);
    }

    // --- Cleanup migration list used in the main loop ---
    ierr = PetscFree(migrationList_main_loop); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


////////////////////////////////////////////////



/**
 * @brief Checks for particles outside the physical domain boundaries and removes them
 *        using DMSwarmRemovePointAtIndex.
 *
 * This function iterates through all particles local to the current MPI rank.
 * It checks if a particle's position (x, y, or z) is outside the specified
 * physical domain boundaries [xMin, xMax], [yMin, yMax], [zMin, zMax].
 *
 * If a particle is found out of bounds, it is removed using DMSwarmRemovePointAtIndex.
 * NOTE: Removing points changes the indices of subsequent points in the iteration.
 *       Therefore, it's crucial to iterate BACKWARDS or carefully manage indices
 *       after a removal. Iterating backwards is generally safer.
 *
 * @param user    Pointer to the UserCtx structure.
 * @param[out] removedCountLocal Pointer to store the number of particles removed *on this rank*.
 * @param[out] removedCountGlobal Pointer to store the total number of particles removed *across all ranks*.
 * @param[in]      bboxlist       An array of BoundingBox structures for ALL MPI ranks, indexed 0 to (size-1).
 *                                This array must be up-to-date and available on all ranks.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
/*
PetscErrorCode CheckAndRemoveOutOfBoundsParticles(UserCtx *user,
                                              PetscInt *removedCountLocal,
					      PetscInt *removedCountGlobal,
					      const BoundingBox *bboxlist )
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscInt       nLocalInitial, p;
    PetscReal        *pos = NULL;
    //   PetscInt64       *PIDs = NULL;
    // You only *need* the position field to check bounds
    PetscInt       local_removed_count = 0;
    PetscMPIInt       global_removed_count = 0;
    PetscMPIInt    rank,size;
    
    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "[Rank %d]Checking for out-of-bounds particles (using RemovePointAtIndex)...\n", rank);

    *removedCountLocal = 0;
    *removedCountGlobal = 0; // Will be summed later
    
    ierr = DMSwarmGetLocalSize(swarm, &nLocalInitial); CHKERRQ(ierr);

    // Get read-only access to positions
    // Note: We get it before the loop because removal might invalidate pointers
    // if the swarm reallocates internally, although Get/Restore is generally safe.
    // For maximum safety with removal, access inside the loop might be better,
    // but less efficient. Let's try getting it once first.

    if(nLocalInitial == 0){
      LOG_ALLOW(LOCAL,LOG_DEBUG," [Rank %d] No  Particles in this rank!.\n",rank);
    }
    
    else{
      
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
    //   ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&PIDs);CHKERRQ(ierr);
     if (!pos) {
         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Could not access position field.");
     }
    // --- Iterate BACKWARDS to handle index changes during removal ---
    local_removed_count = 0;
    
    for (p = nLocalInitial - 1; p >= 0; p--) { // Loop from n-1 down to 0

        Particle temp_particle;
        temp_particle.loc.x = pos[3*p + 0];
        temp_particle.loc.y = pos[3*p + 1];
        temp_particle.loc.z = pos[3*p + 2];
        temp_particle.PID   = -1; // PID isn't needed for this check
      
        PetscBool isOutOfBounds = PETSC_FALSE;

	PetscInt  IsParticleInLocalBboxCount = 0;
	
        //if (pos[p].x < xMin || pos[p].x > xMax ||
        //    pos[p].y < yMin || pos[p].y > yMax ||
        //    pos[p].z < zMin || pos[p].z > zMax)

	// Going through the bbox of each rank, if a particle is found to be in it's bbox, we increment IsParticleInLocalBboxCount
	// IsParticleInLocalBboxCount should be always 1 (as only one rank will contain the particle) or 0 (no rank contains the particle, so it is not incremented).
	for (PetscMPIInt proc = 0; proc < size; proc++){
	  if(IsParticleInsideBoundingBox(&bboxlist[proc],&temp_particle)){
	    IsParticleInLocalBboxCount++;
	  }
	}

	if(IsParticleInLocalBboxCount == 0){
	  isOutOfBounds = PETSC_TRUE;
	}

        if (isOutOfBounds) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Removing at local index %d at (%g,%g,%g).\n",
                      rank, p, temp_particle.loc.x, temp_particle.loc.y,temp_particle.loc.z);
            // Restore position field BEFORE removing point, as removal modifies the swarm state
            ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
	    //  ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&PIDs); CHKERRQ(ierr);
            pos = NULL; // Invalidate pointer.
	    
            // --- Remove the particle at the current local index 'p' ---
            ierr = DMSwarmRemovePointAtIndex(swarm, p); CHKERRQ(ierr);
            local_removed_count++;

            // --- Re-get fields if we continue looping (or break) ---
            // Since we removed a point, indices might have shifted if implementation
            // involves compaction. Re-getting fields is safest if loop continues.
            // However, since we iterate backwards, indices p-1, p-2 etc. remain valid
            // relative to the *current* state after removal.
            // Let's try without re-getting inside the loop first for efficiency.
            // We need to re-get 'pos' for the next iteration check.
            if (p > 0) { // Only re-get if there are more iterations
	      PetscInt nLocalAfterRemoval;
	      ierr = DMSwarmGetLocalSize(swarm,&nLocalAfterRemoval);
	      if(nLocalAfterRemoval>0){
		ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
               if (!pos) {
                   SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Could not re-access position field after removal.");
                   }
	        }else{
		  // All remaining particles were removed, or this was the last one
                  // No need to get 'pos' again as the loop will terminate or p will be 0.
		  break;
	          }
	    }
	  }
        } // End backwards loop

       // Restore position field if it was gotten last time inside loop or if loop didn't run/remove
       if (pos) {
         ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
       }
    } // nLocalInitial > 0

    // Get the *final* local size after removals
    PetscInt nLocalFinal;
    ierr = DMSwarmGetLocalSize(swarm, &nLocalFinal); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "[Rank %d] Finished removing %d particles. Final local size: %d.\n", rank, local_removed_count, nLocalFinal);


    // Calculate global removed count
    *removedCountLocal = local_removed_count;
    ierr = MPI_Allreduce(&local_removed_count, &global_removed_count, 1, MPIU_INT, MPI_SUM, PetscObjectComm((PetscObject)swarm)); CHKERRQ(ierr);
    *removedCountGlobal = global_removed_count;

    LOG_ALLOW(LOCAL, LOG_INFO, "[Rank %d]Removed %d particles globally.\n", rank,global_removed_count);

    PetscFunctionReturn(0);
}
*/


/**
 * @brief Removes particles that have been definitively flagged as LOST by the location algorithm.
 *
 * This function is the designated cleanup utility. It should be called after the
 * `LocateAllParticlesInGrid` orchestrator has run and every particle's status
 * has been definitively determined.
 *
 * It iterates through all locally owned particles and checks their `DMSwarm_location_status`
 * field. If a particle's status is `LOST`, it is permanently removed from the simulation
 * using `DMSwarmRemovePointAtIndex`.
 *
 * This approach centralizes the removal logic, making the `DMSwarm_location_status`
 * the single source of truth for a particle's validity, which is more robust than
 * relying on secondary geometric checks (like bounding boxes).
 *
 * @param[in,out]  user              Pointer to the UserCtx structure containing the swarm.
 * @param[out]     removedCountLocal Pointer to store the number of particles removed on this rank.
 * @param[out]     removedCountGlobal Pointer to store the total number of particles removed across all ranks.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
/*
PetscErrorCode CheckAndRemoveLostParticles(UserCtx *user,
                                           PetscInt *removedCountLocal,
                                           PetscInt *removedCountGlobal)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscInt       nLocalInitial;
    PetscInt       *status_p = NULL;
    PetscInt64     *pid_p = NULL; // For better logging
    PetscReal      *pos_p = NULL; // For better logging
    PetscInt       local_removed_count = 0;
    PetscMPIInt    global_removed_count_mpi = 0;
    PetscMPIInt    rank;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Checking for and removing LOST particles...\n", rank);

    *removedCountLocal = 0;
    if (removedCountGlobal) *removedCountGlobal = 0;

    ierr = DMSwarmGetLocalSize(swarm, &nLocalInitial); CHKERRQ(ierr);
    if (nLocalInitial == 0) {
        PetscFunctionReturn(0); // Nothing to do
    }

    // We only need the status field for the check, but getting PID and position makes for better logs.
    ierr = DMSwarmGetField(swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status_p); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid",             NULL, NULL, (void **)&pid_p);    CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "position",                NULL, NULL, (void **)&pos_p);    CHKERRQ(ierr);

    // --- Iterate BACKWARDS to handle index changes safely during removal ---
    for (PetscInt p = nLocalInitial - 1; p >= 0; p--) {
        if (status_p[p] == LOST) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Removing LOST particle [PID %lld] at local index %d. Position: (%g, %g, %g).\n",
                      rank, (long long)pid_p[p], p, pos_p[3*p], pos_p[3*p+1], pos_p[3*p+2]);

            // Restore all fields BEFORE modifying the swarm structure.
            // This is the safest pattern when removing points.
            ierr = DMSwarmRestoreField(swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status_p); CHKERRQ(ierr);
            ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid",             NULL, NULL, (void **)&pid_p);    CHKERRQ(ierr);
            ierr = DMSwarmRestoreField(swarm, "position",                NULL, NULL, (void **)&pos_p);    CHKERRQ(ierr);

            // --- Remove the particle at the current local index 'p' ---
            ierr = DMSwarmRemovePointAtIndex(swarm, p); CHKERRQ(ierr);
            local_removed_count++;

            // --- After removal, the swarm is modified. We MUST re-acquire pointers. ---
            PetscInt nLocalCurrent;
            ierr = DMSwarmGetLocalSize(swarm, &nLocalCurrent); CHKERRQ(ierr);
            if (nLocalCurrent > 0 && p > 0) { // If there are still particles left to check
                ierr = DMSwarmGetField(swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status_p); CHKERRQ(ierr);
                ierr = DMSwarmGetField(swarm, "DMSwarm_pid",             NULL, NULL, (void **)&pid_p);    CHKERRQ(ierr);
                ierr = DMSwarmGetField(swarm, "position",                NULL, NULL, (void **)&pos_p);    CHKERRQ(ierr);
            } else {
                // All remaining particles were removed, or this was the last one.
                // Invalidate pointers and break the loop since there's nothing left to check.
                status_p = NULL;
                pid_p = NULL;
                pos_p = NULL;
                break;
            }
        }
    } // End of backwards loop

    // Restore fields if the loop completed without a final removal (i.e., pointers are still valid)
    if (status_p) {
        ierr = DMSwarmRestoreField(swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status_p); CHKERRQ(ierr);
    }
    if (pid_p) {
        ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid",             NULL, NULL, (void **)&pid_p);    CHKERRQ(ierr);
    }
    if (pos_p) {
        ierr = DMSwarmRestoreField(swarm, "position",                NULL, NULL, (void **)&pos_p);    CHKERRQ(ierr);
    }

    // Get the *final* local size after removals for logging
    PetscInt nLocalFinal;
    ierr = DMSwarmGetLocalSize(swarm, &nLocalFinal); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Finished removing %d LOST particles. Final local size: %d.\n", rank, local_removed_count, nLocalFinal);

    // --- Synchronize counts across all ranks ---
    *removedCountLocal = local_removed_count;
    if (removedCountGlobal) {
        ierr = MPI_Allreduce(&local_removed_count, &global_removed_count_mpi, 1, MPI_INT, MPI_SUM, PetscObjectComm((PetscObject)swarm)); CHKERRQ(ierr);
        *removedCountGlobal = global_removed_count_mpi;
        LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "[Rank %d] Removed %d particles globally.\n",rank, *removedCountGlobal);
    }

    PetscFunctionReturn(0);
}
*/
