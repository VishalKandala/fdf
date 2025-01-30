/*****************************************************************************/
/*   postprocess.c

   This file implements a simple post-processing executable that:
    1) Initializes PETSc and MPI.
    2) Reads command-line options for input data (time index, field name, etc.).
    3) Reads .dat file(s) into arrays.
    4) Decides whether to write a .vts or .vtp file based on '-out_ext'.
    5) Constructs a VTKMetaData object and calls CreateVTKFileFromMetadata().
    6) Finalizes PETSc.

/*****************************************************************************/

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscdmcomposite.h>

#include "postprocess.h"
#include "logging.h"
#include "common.h"
// #include "io.h"
// #include "grid.h"
// #include "interpolation.h"
// #include "ParticleSwarm.h"

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
              "VecToArrayOnRank0 - rank=%d of %d, globalSize=%D, localSize=%D.\n",
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
 * @brief Main entry point for the PETSc-based VTK post-processing tool.
 *
 * This function demonstrates how to read distributed data (coordinates and field values)
 * from PETSc Vecs, gather them to rank 0, and write them into VTK file formats (.vts or .vtp)
 * based on the command-line options. It leverages the functions defined elsewhere for reading
 * data (ReadFieldData), gathering vectors (VecToArrayOnRank0), and creating VTK files
 * (CreateVTKFileFromMetadata). 
 *
 * Command-line options:
 * - \c -ti <time_index> : Time index for naming convention.
 * - \c -field <field_name> : Name of the field to read and write (default "velocity").
 * - \c -out_ext <extension> : Desired VTK file extension, either "vts" or "vtp" (default "vtp").
 *
 * @param[in] argc Number of command-line arguments.
 * @param[in] argv Array of command-line argument strings.
 *
 * @return int Returns 0 on success, or a non-zero value on errors (e.g., file I/O failures).
 */
int main(int argc, char **argv)
{
    // Log the start of the main function
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "main - Starting PETSc-based VTK post-processing.\n");

    PetscErrorCode ierr;
    int            rank, size;
    int            ti           = 0;              /* time index */
    char           field_name[64] = "velocity";   /* default field name */
    char           out_ext[16]    = "vtp";        /* default => .vtp */
    double        *coordsArray    = NULL;         /* rank-0 coords array */
    double        *scalarArray    = NULL;         /* rank-0 field array */
    PetscInt       Ncoords        = 0;            /* total # of coords */
    PetscInt       Nscalars       = 0;            /* total # of scalars */

    /* 1) Initialize PETSc */
    ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "main - PETSc initialized. rank=%d of %d.\n", rank, size);

    /* 2) Parse command-line options. Adjust as needed in your code. */
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "main - Parsing command-line options.\n");
    PetscOptionsGetInt(NULL, NULL, "-ti", &ti, NULL);
    PetscOptionsGetString(NULL, NULL, "-field", field_name, sizeof(field_name), NULL);
    PetscOptionsGetString(NULL, NULL, "-out_ext", out_ext, sizeof(out_ext), NULL);

    /* 3) Build a user context. (In your code, you might do more here.) */
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "main - Building user context.\n");
    UserCtx user;
    user._this = 0;  /* e.g., partition index => appended to file name. */

    /* 4) Read coordinate data into a PETSc Vec, then gather into coordsArray */
    {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "main - Reading and gathering coordinate data.\n");
        Vec coordsVec;
        ierr = VecCreate(PETSC_COMM_WORLD, &coordsVec);CHKERRQ(ierr);
        ierr = VecSetFromOptions(coordsVec);CHKERRQ(ierr);

        /* We read from "results/position%05d_%d.dat", etc. */
        ierr = ReadFieldData(&user, "position", coordsVec, ti, "dat");
        if (ierr) {
            LOG_ALLOW(GLOBAL, LOG_ERROR, 
                      "main - Error reading position data (ti=%d). Aborting.\n", ti);
            goto finalize;
        }

        /* Gather distributed Vec -> rank-0 array */
        ierr = VecToArrayOnRank0(coordsVec, &Ncoords, &coordsArray);CHKERRQ(ierr);

        ierr = VecDestroy(&coordsVec);CHKERRQ(ierr);
    }

    /* 5) Read the field data (velocity) into a PETSc Vec, gather into scalarArray. */
    {
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "main - Reading and gathering field data '%s'.\n", field_name);
        Vec fieldVec;
        ierr = VecCreate(PETSC_COMM_WORLD, &fieldVec);CHKERRQ(ierr);
        ierr = VecSetFromOptions(fieldVec);CHKERRQ(ierr);

        ierr = ReadFieldData(&user, field_name, fieldVec, ti, "dat");
        if (ierr) {
            LOG_ALLOW(GLOBAL, LOG_ERROR, 
                      "main - Error reading field data '%s' (ti=%d). Aborting.\n",
                      field_name, ti);
            goto finalize;
        }

        ierr = VecToArrayOnRank0(fieldVec, &Nscalars, &scalarArray);CHKERRQ(ierr);

        ierr = VecDestroy(&fieldVec);CHKERRQ(ierr);
    }

    /* 6) Prepare a VTKMetaData struct to describe how to interpret coords & field. */
    LOG_ALLOW(GLOBAL, LOG_DEBUG, 
              "main - Preparing VTKMetaData (field=%s, out_ext=%s).\n", field_name, out_ext);
    VTKMetaData meta;
    memset(&meta, 0, sizeof(meta)); /* zero out */

    /* Decide if we want .vts (structured) or .vtp (polydata) by comparing out_ext */
    PetscBool isVTP;
    {
        PetscBool match;
        ierr = PetscStrcmp(out_ext, "vtp", &match);CHKERRQ(ierr);
        isVTP = match ? PETSC_TRUE : PETSC_FALSE;
    }

    if (!isVTP) {
        /* => VTK_STRUCTURED => .vts */
        meta.fileType = VTK_STRUCTURED;
        meta.mx = 4; meta.my = 4; meta.mz = 4;
        meta.nnodes = (meta.mx - 1) * (meta.my - 1) * (meta.mz - 1);

        if (!rank) {
            meta.coords = coordsArray;           
            meta.scalarField = scalarArray;
            meta.scalarFieldName = field_name;
            meta.numScalarFields = 1;
        }
    } else {
        /* => VTK_POLYDATA => .vtp */
        meta.fileType = VTK_POLYDATA;
        if (!rank) {
            meta.coords = coordsArray; 

            /* Check coords size is multiple of 3. */
            if (Ncoords % 3 != 0) {
                PetscPrintf(PETSC_COMM_SELF, 
                            "Error: coords array length %d not multiple of 3.\n", (int)Ncoords);
                LOG_ALLOW(GLOBAL, LOG_ERROR, 
                          "main - coords array length %D not multiple of 3, abort.\n", Ncoords);
                goto finalize;
            }
            meta.npoints = Ncoords / 3;

            /* If field is 3*Npoints, assume vector field. Otherwise assume scalar. */
            if (Nscalars == meta.npoints * 3) {
                meta.vectorField = scalarArray;
                meta.vectorFieldName = field_name;
                meta.numVectorFields = 1;
                meta.numScalarFields = 0;
            } else {
                meta.scalarField = scalarArray;
                meta.scalarFieldName = field_name;
                meta.numScalarFields = 1;
                meta.numVectorFields = 0;
            }

            /* Create connectivity & offsets => 1 cell (vertex) per point. */
            meta.connectivity = (int *)calloc(meta.npoints, sizeof(int));
            meta.offsets      = (int *)calloc(meta.npoints, sizeof(int));
            for (int i = 0; i < meta.npoints; i++) {
                meta.connectivity[i] = i;     
                meta.offsets[i]      = i + 1; 
            }
        }
    }

    /* 7) Construct output file name. E.g., "results/velocity00000.vtp" */
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "main - Constructing output file name.\n");
    char outFile[256];
    PetscSNPrintf(outFile, sizeof(outFile), "results/%s%05d.%s", field_name, ti, out_ext);

    /* 8) Call the file writer */
    {
        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "main - Creating VTK file '%s'.\n", outFile);
        int err = CreateVTKFileFromMetadata(outFile, &meta, PETSC_COMM_WORLD);
        if (err) {
            PetscPrintf(PETSC_COMM_WORLD,
                        "[ERROR] CreateVTKFileFromMetadata returned %d.\n", err);
            LOG_ALLOW(GLOBAL, LOG_ERROR, 
                      "main - CreateVTKFileFromMetadata failed with code=%d.\n", err);
            goto finalize;
        }
        if (!rank) {
            PetscPrintf(PETSC_COMM_SELF, "[postprocess] Wrote file: %s\n", outFile);
            LOG_ALLOW(GLOBAL, LOG_INFO, "main - Successfully wrote file: %s\n", outFile);
        }
    }

finalize:
    /* 9) Clean up memory on rank 0. */
    if (!rank) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "main - Cleaning up memory on rank 0.\n");
        free(coordsArray);
        free(scalarArray);
        if (meta.fileType == VTK_POLYDATA) {
            free(meta.connectivity);
            free(meta.offsets);
        }
    }

    /* 10) Finalize PETSc. */
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "main - Finalizing PETSc.\n");
    PetscFinalize();
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "main - Program completed.\n");
    return 0;
}
