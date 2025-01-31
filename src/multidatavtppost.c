/**
 * @brief Gathers particle positions from user->swarm, then calls the same VTK pipeline.
 */
PetscErrorCode WriteParticleVTK(UserCtx *user, PetscInt timeIndex, const char *outExt)
{
    PetscFunctionBeginUser;
    MPI_Comm    comm = PETSC_COMM_WORLD;
    PetscMPIInt rank;
    MPI_Comm_rank(comm, &rank);

    // a) Gather particle positions into an array (x,y,z).
    double   *partCoords = NULL; 
    PetscInt  Nparticles = 0;
    // e.g. "GatherParticlePositions(user->swarm, &Nparticles, &partCoords)"

    // b) For particles, we might treat them as polydata (vtp).
    //    Optionally gather a scalar field (like particle temperature).
    double *particleScalars = NULL;
    PetscInt Nscalars       = 0; 
    // gather if needed

    // c) Prepare metadata
    VTKMetaData meta;
    PrepareVTKMetaData(partCoords, 3*Nparticles,
                       particleScalars, Nscalars,
                       "Particles", outExt,
                       &meta, rank);

    // d) Construct filename
    char outFile[256];
    ConstructOutputFilename("particles", timeIndex, outExt, outFile, sizeof(outFile));

    // e) Write
    CreateAndWriteVTKFile(outFile, &meta, comm);

    // f) Cleanup on rank 0
    if (!rank) {
        if (meta.fileType == 2 /* VTK_POLYDATA */) {
            free(meta.connectivity);
            free(meta.offsets);
        }
        // free(partCoords);
        // free(particleScalars);
    }

    PetscFunctionReturn(0);
}


#include <petscsys.h>
#include <string.h>
#include <stdio.h>

/**
 * @brief Prepares a \c VTKMetaData object for .vtp output with:
 *        - \p coordsArray as <Points>
 *        - \p velArray as a single vector field in <PointData>
 *
 * Future fields can be added similarly by incrementing meta->numFields.
 *
 * @param[in]  coordsArray  3*Npoints array of positions (valid on rank 0).
 * @param[in]  Ncoords      Length of \p coordsArray (should be 3*Npoints).
 * @param[in]  velArray     3*Npoints array of velocity components (valid on rank 0).
 * @param[in]  Nvel         Length of \p velArray (should be 3*Npoints).
 * @param[in]  rank         The MPI rank (only rank 0 populates the struct).
 * @param[out] meta         The VTKMetaData struct to fill for .vtp writing.
 *
 * @return PetscErrorCode (0 = success).
 */
PetscErrorCode PrepareVTPWithPositionsAndVelocity(const double *coordsArray,
                                                  PetscInt Ncoords,
                                                  const double *velArray,
                                                  PetscInt Nvel,
                                                  int rank,
                                                  VTKMetaData *meta)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr = 0;

    /* Initialize the struct to zeros */
    memset(meta, 0, sizeof(*meta));

    meta->fileType = 2; /* VTK_POLYDATA (assuming you define 2 = polydata) */

    if (!rank) {
        /* Check coords array => must be multiple of 3 */
        if (Ncoords % 3 != 0) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP,
                    "coordsArray length must be multiple of 3 for .vtp polydata.\n");
        }
        meta->npoints = (int)(Ncoords / 3);
        meta->coords  = (double*)coordsArray;  /* We store pointer (no copy) */

        /* Create connectivity & offsets for 'npoints' vertices */
        meta->connectivity = (int*)calloc(meta->npoints, sizeof(int));
        meta->offsets      = (int*)calloc(meta->npoints, sizeof(int));
        for (int i = 0; i < meta->npoints; i++) {
            meta->connectivity[i] = i;
            meta->offsets[i]      = i + 1;
        }

        /* Now handle velocity (3*Npoints) => 3 comps => vector field */
        if (Nvel != 3 * meta->npoints) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP,
                    "velArray length must be 3*Npoints for a vector field.\n");
        }

        meta->numFields = 1; /* We'll store exactly 1 field in <PointData> for now */

        strcpy(meta->pointDataFields[0].name, "velocity");
        meta->pointDataFields[0].numComponents = 3;
        meta->pointDataFields[0].data          = (double*)velArray;
    }

    PetscFunctionReturn(ierr);
}

/**
 * @brief Gathers position and velocity arrays (3*N each) onto rank 0, 
 *        prepares a single VTKMetaData for .vtp, and writes one file.
 *
 * @param[in]  coordsVec   PETSc Vec of length 3*N for positions (x,y,z).
 * @param[in]  velVec      PETSc Vec of length 3*N for velocity (vx,vy,vz).
 * @param[in]  npoints     On rank 0, total # of points (N). On other ranks unused.
 * @param[in]  timeIndex   Timestep index for the filename.
 * @param[in]  outExt      Typically "vtp".
 *
 * @return PetscErrorCode  0 on success, non-zero on errors.
 */
PetscErrorCode GatherPosVelWriteSingleVTP(Vec coordsVec, Vec velVec,
                                          PetscInt npoints,
                                          PetscInt timeIndex,
                                          const char *outExt)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    MPI_Comm       comm = PETSC_COMM_WORLD;
    PetscMPIInt    rank;
    ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);

    /* 1) Gather coords => coordsArray (3*N total) */
    PetscInt  Ncoords = 0;
    double   *coordsArray = NULL;
    if (coordsVec) {
        ierr = VecToArrayOnRank0(coordsVec, &Ncoords, &coordsArray);CHKERRQ(ierr);
    }

    /* 2) Gather velocity => velArray (3*N total) */
    PetscInt  Nvel = 0;
    double   *velArray = NULL;
    if (velVec) {
        ierr = VecToArrayOnRank0(velVec, &Nvel, &velArray);CHKERRQ(ierr);
    }

    /* 3) Prepare the VTKMetaData for a single .vtp with coords + velocity */
    VTKMetaData meta;
    ierr = PrepareVTPWithPositionsAndVelocity(coordsArray, Ncoords,
                                              velArray, Nvel,
                                              rank, &meta);CHKERRQ(ierr);

    /* 4) Construct output filename, e.g. "particles_00010.vtp" */
    char outFile[256];
    sprintf(outFile, "particles_%05d.%s", (int)timeIndex, outExt);

    /* 5) Actually write the file (rank 0 does the I/O) */
    int errCode = CreateVTKFileFromMetadata(outFile, &meta, comm);
    if (errCode != 0) {
        PetscPrintf(comm, "[ERROR] CreateVTKFileFromMetadata returned %d for %s.\n", errCode, outFile);
    }

    /* 6) Cleanup arrays on rank 0 */
    if (!rank) {
        free(meta.connectivity);
        free(meta.offsets);
        free(coordsArray);
        free(velArray);
    }

    PetscFunctionReturn(0);
}
