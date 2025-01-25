/*****************************************************************************/
/* postprocess.c                                                              */
/*                                                                            */
/* This file contains routines for writing data to VTK files (XML/.vts        */
/* format), as well as other post-processing and output-related functions.    */
/*                                                                            */
/* Example includes:                                                          */
/*  - GatherParallelVecToSeq() for gathering PETSc Vecs to rank 0             */
/*  - WriteVTKXMLHeader(), WriteVTKXMLFooter() for generating VTK XML markup  */
/*  - WriteVTKAppendedBlock() for appended binary data in VTK                 */
/*  - CreateVTKFileFromData() for producing a .vts file from a scalar Vec     */
/*  - Additional placeholders for extended post-processing routines           */
/*    (e.g., multi-field outputs, vector fields, etc.).                       */
/*                                                                            */
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

/* 
   postprocess.c

   This file implements a simple post-processing executable that:
    1) Initializes PETSc and MPI.
    2) Reads command-line options for input data (time index, field name, etc.).
    3) Reads .dat file(s) into arrays.
    4) Decides whether to write a .vts or .vtp file based on '-out_ext'.
    5) Constructs a VTKMetaData object and calls CreateVTKFileFromMetadata().
    6) Finalizes PETSc.
*/

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

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    int            rank, size;
    int            ti           = 0;             
    char           field_name[64] = "velocity";  
    char           out_ext[16]    = "vtp";       
    double        *coordsArray  = NULL;  // Global pointer for coordinates
    double        *scalarArray  = NULL;  // Global pointer for scalar field
    PetscInt       Ncoords = 0;          // Global size of coordsArray
    PetscInt       Nscalars = 0;         // Global size of scalarArray

    // Initialize PETSc
    ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);

    UserCtx user; 
    user._this = 0; // e.g., partition index

    // Parse command line options
    PetscOptionsGetInt(NULL, NULL, "-ti", &ti, NULL);
    PetscOptionsGetString(NULL, NULL, "-field", field_name, sizeof(field_name), NULL);
    PetscOptionsGetString(NULL, NULL, "-out_ext", out_ext, sizeof(out_ext), NULL);

    // STEP: Read coords
    {
        Vec coordsVec;
        ierr = VecCreate(PETSC_COMM_WORLD, &coordsVec);CHKERRQ(ierr);
        ierr = VecSetFromOptions(coordsVec);CHKERRQ(ierr);

        // e.g. "position" data
        ierr = ReadFieldData(&user, "position", coordsVec, ti, "dat");
        if (ierr) goto finalize;

        // Gather to rank 0
        ierr = VecToArrayOnRank0(coordsVec, &Ncoords, &coordsArray);CHKERRQ(ierr);
        ierr = VecDestroy(&coordsVec);CHKERRQ(ierr);

        PetscPrintf(PETSC_COMM_WORLD,
                    "   ... Successfully read %d coordinate values.\n", (int)Ncoords);
    }

    // STEP: Read scalar field
    {
        Vec fieldVec;
        ierr = VecCreate(PETSC_COMM_WORLD, &fieldVec);CHKERRQ(ierr);
        ierr = VecSetFromOptions(fieldVec);CHKERRQ(ierr);

        ierr = ReadFieldData(&user, field_name, fieldVec, ti, "dat");
        if (ierr) goto finalize;

        // Gather to rank 0
        ierr = VecToArrayOnRank0(fieldVec, &Nscalars, &scalarArray);CHKERRQ(ierr);
        ierr = VecDestroy(&fieldVec);CHKERRQ(ierr);

        PetscPrintf(PETSC_COMM_WORLD,
                    "   ... Successfully read %d scalar values.\n", (int)Nscalars);
    }

    // After reading coordinates and field data:
    PetscPrintf(PETSC_COMM_WORLD, "Ncoords = %d, Nscalars = %d\n",
		(int)Ncoords, (int)Nscalars);

    // Decide .vts or .vtp
    PetscBool isVTP;
    {
        PetscBool match;
        ierr = PetscStrcmp(out_ext, "vtp", &match);CHKERRQ(ierr);
        isVTP = match ? PETSC_TRUE : PETSC_FALSE;
    }

    // Create VTK meta struct
    VTKMetaData meta;
    memset(&meta, 0, sizeof(meta));

    // For .vts or .vtp...
    if (!isVTP) {
        // structured grid...
      // these toy numbers must be augmented to get IM,JM,KM.
        meta.fileType = VTK_STRUCTURED;
        meta.mx = 4;  meta.my = 4;  meta.mz = 4;
        meta.nnodes   = (meta.mx - 1) * (meta.my - 1) * (meta.mz - 1);

        if (!rank) {
            meta.coords          = coordsArray;
            meta.scalarField     = scalarArray;
            meta.scalarFieldName = field_name;
            meta.numScalarFields = 1;
        }
    } else {
        // polydata...
        meta.fileType = VTK_POLYDATA;
        if (!rank) {
            if (Ncoords % 3 != 0) {
	      PetscPrintf(PETSC_COMM_SELF, "Error: Coordinates array size invalid.\n");
                goto finalize;
            }

            meta.npoints         = Ncoords / 3;
	    // After setting metadata npoints:
	    PetscPrintf(PETSC_COMM_WORLD, "Ncoords = %d, npoints = %d, Nscalars = %d\n",
			(int)Ncoords, (int)meta.npoints, (int)Nscalars);
            /* Detect if the field is a 3D vector (3 components per point) */
            if (Nscalars == meta.npoints * 3) {
                meta.vectorField = scalarArray;          // Assign to vector field
                meta.vectorFieldName = field_name;
                meta.numVectorFields = 1;
                meta.numScalarFields = 0;                // Disable scalar fields
            } else {
            meta.scalarField     = scalarArray;
            meta.scalarFieldName = field_name;
            meta.numScalarFields = 1;
	    meta.numVectorFields = 0;
	    }

            // create connectivity, offsets
            meta.connectivity = (int*)calloc(meta.npoints, sizeof(int));
            meta.offsets      = (int*)calloc(meta.npoints, sizeof(int));
            for (int i = 0; i < meta.npoints; i++) {
                meta.connectivity[i] = i;
                meta.offsets[i]      = i + 1;
            }
        } // rank
    } // vtp or vts 

    // Construct output file name
    char outFile[256];
    PetscSNPrintf(outFile, sizeof(outFile),
                  "results/%s%05d.%s", field_name, ti, out_ext);

    // Write VTK file
    ierr = CreateVTKFileFromMetadata(outFile, &meta, PETSC_COMM_WORLD);
    if (ierr) { /* handle error */ goto finalize; }

    if (!rank) {
        PetscPrintf(PETSC_COMM_SELF,
                    "[postprocess] Successfully wrote file: %s\n", outFile);
    }

finalize:
    // Cleanup
    if (!rank) {
        free(coordsArray);
        free(scalarArray);
        if (meta.fileType == VTK_POLYDATA) {
            free(meta.connectivity);
            free(meta.offsets);
        }
    }

    PetscFinalize();
    return 0;
}

