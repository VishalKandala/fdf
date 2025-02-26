
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



/**
 * @brief Safely interpolate a vector field (Cmpnts) from cell centers to corners, including boundary checks.
 *
 * We loop over each corner (i,j,k) and consider up to eight cell centers
 * (i+di, j+dj, k+dk) in {0,1}. Each valid center is averaged.
 * On boundaries, partial stencils are used (skipping out-of-range indices).
 *
 * @param[in]  centfield 3D array (Cmpnts ***) of cell-center vectors
 * @param[out] field     3D array (Cmpnts ***) for corner vectors
 * @param[in]  info      DMDALocalInfo for the cell-center DMDA
 *
 * @return PetscErrorCode
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Vector(Cmpnts ***centfield,
                                                         Cmpnts ***field,
                                                         DMDALocalInfo *info)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "InterpolateFieldFromCenterToCorner_Vector - Rank %d starting vector center-to-corner interpolation.\n", rank);

  PetscInt       mx, my, mz;
  PetscInt       xs, xe, ys, ye, zs, ze;
  PetscInt       i, j, k, di, dj, dk;

  xs = info->xs; xe = info->xs + info->xm;
  ys = info->ys; ye = info->ys + info->ym;
  zs = info->zs; ze = info->zs + info->zm;
  mx = info->mx; my = info->my; mz = info->mz;

  /* We assume the corner domain is slightly bigger: [xs-1..xe], etc., with boundary checks. */
  PetscInt xs_corner = xs - 1;
  PetscInt xe_corner = xe;
  PetscInt ys_corner = ys - 1;
  PetscInt ye_corner = ye;
  PetscInt zs_corner = zs - 1;
  PetscInt ze_corner = ze;

  for (k = zs_corner; k < ze_corner; k++) {
    for (j = ys_corner; j < ye_corner; j++) {
      for (i = xs_corner; i < xe_corner; i++) {

        Cmpnts sum  = {0.0, 0.0, 0.0};
        PetscInt count = 0;

        for (dk = 0; dk < 2; dk++) {
          for (dj = 0; dj < 2; dj++) {
            for (di = 0; di < 2; di++) {
              PetscInt ic = i + di;
              PetscInt jc = j + dj;
              PetscInt kc = k + dk;
              if (ic >= xs && ic < xe &&
                  jc >= ys && jc < ye &&
                  kc >= zs && kc < ze)
              {
                sum.x += centfield[kc][jc][ic].x;
                sum.y += centfield[kc][jc][ic].y;
                sum.z += centfield[kc][jc][ic].z;
                count++;
              }
            }
          }
        }

        if (count > 0) {
          field[k][j][i].x = sum.x / (PetscReal)count;
          field[k][j][i].y = sum.y / (PetscReal)count;
          field[k][j][i].z = sum.z / (PetscReal)count;
        } else {
          field[k][j][i].x = 0.0;
          field[k][j][i].y = 0.0;
          field[k][j][i].z = 0.0;
        }
      }
    }
  }
  
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "InterpolateFieldFromCenterToCorner_Vector - Rank %d completed vector center-to-corner interpolation.\n", rank);
  return 0;
}



PetscErrorCode InterpolateFieldFromCornerToCenter_Vector(
    Cmpnts ***field,         /* coordinate DM array (da, including ghost cells) */
    Cmpnts ***centfield,     /* output: interpolated interior cell–center values for fda */
    UserCtx *user)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
    "InterpolateFieldFromCornerToCenter_Vector - Rank %d starting interpolation.\n", rank);

  /* Obtain local info for fda and da */
  DMDALocalInfo info_fda, info_da;
  ierr = DMDAGetLocalInfo(user->fda, &info_fda); CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(user->da, &info_da); CHKERRQ(ierr);

  /* Obtain global domain dimensions from fda */
  PetscInt gNx, gNy, gNz;
  ierr = DMDAGetInfo(user->fda, NULL, &gNx, &gNy, &gNz, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);

  /* Compute interior bounds for fda:
     If the local start equals the global minimum then skip that boundary cell,
     and similarly for the upper boundary.
  */
  PetscInt lxs = (info_fda.xs == 0) ? info_fda.xs + 1 : info_fda.xs;
  PetscInt lxe = ((info_fda.xs + info_fda.xm) == gNx) ? info_fda.xs + info_fda.xm - 1 : info_fda.xs + info_fda.xm;
  PetscInt lys = (info_fda.ys == 0) ? info_fda.ys + 1 : info_fda.ys;
  PetscInt lye = ((info_fda.ys + info_fda.ym) == gNy) ? info_fda.ys + info_fda.ym - 1 : info_fda.ys + info_fda.ym;
  PetscInt lzs = (info_fda.zs == 0) ? info_fda.zs + 1 : info_fda.zs;
  PetscInt lze = ((info_fda.zs + info_fda.zm) == gNz) ? info_fda.zs + info_fda.zm - 1 : info_fda.zs + info_fda.zm;

  LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
            "Rank %d: fda interior bounds: lxs=%d, lxe=%d, lys=%d, lye=%d, lzs=%d, lze=%d (global dims: %d,%d,%d)\n",
            rank, lxs, lxe, lys, lye, lzs, lze, gNx, gNy, gNz);
  
  /* Loop over the interior cells of fda */
  for (PetscInt k = lzs; k < lze; k++) {
    for (PetscInt j = lys; j < lye; j++) {
      for (PetscInt i = lxs; i < lxe; i++) {

        Cmpnts sum = {0.0, 0.0, 0.0};
        PetscInt count = 0;
        PetscInt iterVar = 0;

        /* For cell center (i,j,k) in fda, its 8 surrounding corners in da have global indices
           from (i-1,j-1,k-1) to (i,j,k). Here we assume that the global indices match between fda and da.
           We then subtract info_da.xs, etc., to index into da's local array.
         */
        for (PetscInt dk = 0; dk < 2; dk++) {
          for (PetscInt dj = 0; dj < 2; dj++) {
            for (PetscInt di = 0; di < 2; di++) {
              PetscInt ci = i - 1 + di;
              PetscInt cj = j - 1 + dj;
              PetscInt ck = k - 1 + dk;
              /* No clamping needed here if ghost cells exist; assume indices are valid.
                 (You may add an assertion if desired.) */
              if ((ci >= info_da.xs) && (ci < info_da.xs + info_da.xm) &&
                  (cj >= info_da.ys) && (cj < info_da.ys + info_da.ym) &&
                  (ck >= info_da.zs) && (ck < info_da.zs + info_da.zm))
              {
                sum.x += field[ck - info_da.zs][cj - info_da.ys][ci - info_da.xs].x;
                sum.y += field[ck - info_da.zs][cj - info_da.ys][ci - info_da.xs].y;
                sum.z += field[ck - info_da.zs][cj - info_da.ys][ci - info_da.xs].z;
                count++;
                iterVar = ci*cj*ck;
		  LOG_LOOP_ALLOW(GLOBAL, LOG_DEBUG,iterVar,10,
                    "Rank %d: Interior cell (i=%d,j=%d,k=%d):corner (ci=%d,cj=%d,ck=%d) da range [%d,%d) [%d,%d) [%d,%d)\n",
                    rank, i, j, k, ci, cj, ck, info_da.xs, info_da.xs+info_da.xm, info_da.ys, info_da.ys+info_da.ym, info_da.zs, info_da.zs+info_da.zm);
              } else {
                iterVar = ci*cj*ck;
		  LOG_LOOP_ALLOW(GLOBAL, LOG_DEBUG,iterVar,1,
                    "Rank %d: Interior cell (i=%d,j=%d,k=%d): skipping corner (ci=%d,cj=%d,ck=%d) outside da range [%d,%d) [%d,%d) [%d,%d)\n",
                    rank, i, j, k, ci, cj, ck, info_da.xs, info_da.xs+info_da.xm, info_da.ys, info_da.ys+info_da.ym, info_da.zs, info_da.zs+info_da.zm);
              }
            }
          }
        }
        if (count > 0) {
          centfield[k - info_fda.zs][j - info_fda.ys][i - info_fda.xs].x = sum.x / count;
          centfield[k - info_fda.zs][j - info_fda.ys][i - info_fda.xs].y = sum.y / count;
          centfield[k - info_fda.zs][j - info_fda.ys][i - info_fda.xs].z = sum.z / count;
        } else {
          centfield[k - info_fda.zs][j - info_fda.ys][i - info_fda.xs].x = 0.0;
          centfield[k - info_fda.zs][j - info_fda.ys][i - info_fda.xs].y = 0.0;
          centfield[k - info_fda.zs][j - info_fda.ys][i - info_fda.xs].z = 0.0;
          LOG_ALLOW(GLOBAL, LOG_DEBUG,
                    "Rank %d: Interior cell (i=%d,j=%d,k=%d) got no valid corner data.\n",
                    rank, i, j, k);
        }
      }
    }
  }

  LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
            "InterpolateFieldFromCornerToCenter_Vector - Rank %d completed interpolation.\n", rank);
  return 0;
}
