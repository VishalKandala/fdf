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
 * @brief Gathers the contents of a distributed PETSc Vec into a single array on rank 0.
 *
 * @param[in]  inVec       The input (possibly distributed) Vec.
 * @param[out] N           The global size of the vector.
 * @param[out] arrayOut    On rank 0, points to the newly allocated array holding all data.
 *                         On other ranks, it is set to NULL.
 *
 * @return PetscErrorCode  Return 0 on success, nonzero on failure.
 */
PetscErrorCode VecToArrayOnRank0(Vec inVec, PetscInt *N, double **arrayOut)
{
    PetscErrorCode    ierr;
    MPI_Comm          comm;
    PetscMPIInt       rank, size;
    PetscInt          globalSize, localSize;
    const PetscScalar *localArr = NULL;
    
    /* Get MPI comm, rank, size */
    ierr = PetscObjectGetComm((PetscObject)inVec, &comm);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm, &size);CHKERRQ(ierr);

    /* Get global size (for the entire Vec) */
    ierr = VecGetSize(inVec, &globalSize);CHKERRQ(ierr);
    *N = globalSize;

    /* Get local size (portion on this rank) */
    ierr = VecGetLocalSize(inVec, &localSize);CHKERRQ(ierr);

    /* Access the local array data */
    ierr = VecGetArrayRead(inVec, &localArr);CHKERRQ(ierr);

    /*
       We'll gather the local chunks via MPI_Gatherv.
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

    PetscFunctionReturn(0);
}


/**************************************************************
 * MAIN POSTPROCESSING PROGRAM
 *
 * Reads two PETSc vectors:
 *  1) position (coordinates)
 *  2) velocity (scalar or vector)
 *
 * Gathers them into rank-0 arrays (coordsArray, scalarArray),
 * then calls CreateVTKFileFromMetadata() to produce .vts or .vtp.
 *
 * KEY CORRECTION:
 *  - For .vtp, we must set meta->coords to coordsArray
 *    so that the coordinate block is actually written.
 **************************************************************/

int main(int argc, char **argv)
{
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

    /* 2) Parse command-line options. Adjust as needed in your code. */
    PetscOptionsGetInt(NULL, NULL, "-ti", &ti, NULL);
    PetscOptionsGetString(NULL, NULL, "-field", field_name, sizeof(field_name), NULL);
    PetscOptionsGetString(NULL, NULL, "-out_ext", out_ext, sizeof(out_ext), NULL);

    /* 3) Build a user context. (In your code, you might do more here.) */
    UserCtx user;
    user._this = 0;  /* e.g., partition index => appended to file name. */

    /* 4) Read coordinate data into a PETSc Vec, then gather into coordsArray */
    {
        Vec coordsVec;
        ierr = VecCreate(PETSC_COMM_WORLD, &coordsVec);CHKERRQ(ierr);
        ierr = VecSetFromOptions(coordsVec);CHKERRQ(ierr);

        /* We read from "results/position%05d_%d.dat", etc. */
        ierr = ReadFieldData(&user, "position", coordsVec, ti, "dat");
        if (ierr) {
            /* Could handle the error or goto finalize */
            goto finalize;
        }

        /* Gather distributed Vec -> rank-0 array */
        ierr = VecToArrayOnRank0(coordsVec, &Ncoords, &coordsArray);CHKERRQ(ierr);

        ierr = VecDestroy(&coordsVec);CHKERRQ(ierr);
    }

    /* 5) Read the field data (velocity) into a PETSc Vec, gather into scalarArray. */
    {
        Vec fieldVec;
        ierr = VecCreate(PETSC_COMM_WORLD, &fieldVec);CHKERRQ(ierr);
        ierr = VecSetFromOptions(fieldVec);CHKERRQ(ierr);

        ierr = ReadFieldData(&user, field_name, fieldVec, ti, "dat");
        if (ierr) {
            /* Handle error or skip. */
            goto finalize;
        }

        ierr = VecToArrayOnRank0(fieldVec, &Nscalars, &scalarArray);CHKERRQ(ierr);

        ierr = VecDestroy(&fieldVec);CHKERRQ(ierr);
    }

    /* 6) Prepare a VTKMetaData struct to describe how to interpret coords & field. */
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
        /* Suppose you know real domain dims: (mx, my, mz).
           For demonstration we pick 4,4,4. Adjust as needed. */
        meta.mx = 4; meta.my = 4; meta.mz = 4;
        meta.nnodes = (meta.mx - 1) * (meta.my - 1) * (meta.mz - 1);

        if (!rank) {
            meta.coords = coordsArray;           /* must not be NULL now */
            meta.scalarField = scalarArray;
            meta.scalarFieldName = field_name;
            meta.numScalarFields = 1;
        }
    } else {
        /* => VTK_POLYDATA => .vtp */
        meta.fileType = VTK_POLYDATA;
        if (!rank) {
            /* We must set meta->coords here! Otherwise no coords are written. */
            meta.coords = coordsArray; /* CRUCIAL FIX. */

            /* Check coords size is multiple of 3. */
            if (Ncoords % 3 != 0) {
                PetscPrintf(PETSC_COMM_SELF, "Error: coords array length %d not multiple of 3.\n", (int)Ncoords);
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
                meta.connectivity[i] = i;      /* each cell has 1 index: i */
                meta.offsets[i]      = i + 1;  /* cell boundary at i+1 */
            }
        }
    }

    /* 7) Construct output file name. E.g., "results/velocity00000.vtp" */
    char outFile[256];
    PetscSNPrintf(outFile, sizeof(outFile), "results/%s%05d.%s", field_name, ti, out_ext);

    /* 8) Call the file writer */
    {
        int err = CreateVTKFileFromMetadata(outFile, &meta, PETSC_COMM_WORLD);
        if (err) {
            PetscPrintf(PETSC_COMM_WORLD,
                        "[ERROR] CreateVTKFileFromMetadata returned %d.\n", err);
            goto finalize;
        }
        if (!rank) {
            PetscPrintf(PETSC_COMM_SELF, "[postprocess] Wrote file: %s\n", outFile);
        }
    }

finalize:
    /* 9) Clean up memory on rank 0. */
    if (!rank) {
        free(coordsArray);
        free(scalarArray);
        if (meta.fileType == VTK_POLYDATA) {
            free(meta.connectivity);
            free(meta.offsets);
        }
    }

    /* 10) Finalize PETSc. */
    PetscFinalize();
    return 0;
}
