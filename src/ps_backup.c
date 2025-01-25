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

   Logging macro usage:
     LOG_ALLOW(GLOBAL, LOG_DEBUG,   "Debug message...\n");
     LOG_ALLOW(GLOBAL, LOG_INFO,    "Info message...\n");
     LOG_ALLOW(GLOBAL, LOG_WARNING, "Warning message...\n");
   
   Adjust logging levels, arguments, etc. to match your project’s logging setup.
*/

/* 
   We assume you have includes in your build environment for:
   - <petscsys.h>
   - <mpi.h>
   - <stdio.h>
   - <string.h>
   - "postprocess.h"
   - "io.h"
   etc.
*/

/* Forward declaration (if needed):
   int PostprocessMain(int argc, char **argv);
   int main(int argc, char **argv);
*/

int PostprocessMain(int argc, char **argv)
{
    PetscErrorCode ierr; /* or int if not using PETSc error codes */
    int            rank, size;
    int            ti           = 0;            /* time index */
    char           field_name[64] = "velocity"; /* default field */
    char           out_ext[16]    = "vts";      /* default extension => .vts */
    char           coordsFile[256], fieldFile[256];
    double        *coordsArray  = NULL; /* will store coordinate data from .dat */
    double        *scalarArray  = NULL; /* will store field data from .dat */
    int            Ncoords = 0;         /* length of coordsArray */
    int            Nscalars = 0;        /* length of scalarArray */

    /* Example of how your logging system can be told to allow logs from certain functions: */
    const char *allowedFuncs[] = {
        "PostProcessMain" /* We'll allow logging from this main function */
        /* You can add more function names if you have logs there */
    };
    set_allowed_functions(allowedFuncs, 1); /* set_allowed_functions from your logging framework */

    /* Step 0: Starting main postprocessing routine */
    PetscPrintf(PETSC_COMM_WORLD, "=== STEP 0: Entering PostprocessMain ===\n");
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "PostprocessMain - Starting.\n");

    /* 1) Initialize PETSc (and MPI). */
    PetscPrintf(PETSC_COMM_WORLD, "STEP 1: Initializing PETSc...\n");
    ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO,
              "PostprocessMain - Running on %d MPI ranks.\n", size);
    PetscPrintf(PETSC_COMM_WORLD, "  ... PETSc initialized. Running on %d ranks.\n", size);

    /* 2) Get command line options (-ti, -field, -out_ext). 
          If they aren't provided, we stick to the defaults. */
    PetscPrintf(PETSC_COMM_WORLD, "STEP 2: Parsing command-line options...\n");
    PetscOptionsGetInt(NULL, NULL, "-ti", &ti, NULL);
    PetscOptionsGetString(NULL, NULL, "-field", field_name, sizeof(field_name), NULL);
    PetscOptionsGetString(NULL, NULL, "-out_ext", out_ext, sizeof(out_ext), NULL);

    if (!rank) {
        PetscPrintf(PETSC_COMM_SELF,
                    "\n[postprocess] Using:\n"
                    "  time index (ti) = %d\n"
                    "  field name      = %s\n"
                    "  out extension   = %s\n\n",
                    ti, field_name, out_ext);
    }

    /* 3) Build filenames for reading coordinate data and the field data.
          For demonstration:
            - coords => "results/position<ti>_0.dat"
            - field  => "results/<field_name><ti>_0.dat"
       (assuming some naming conventions in your code)
    */
    PetscPrintf(PETSC_COMM_WORLD, "STEP 3: Constructing input file names...\n");
    PetscSNPrintf(coordsFile, sizeof(coordsFile), "results/position%05d_0.dat", ti);
    PetscSNPrintf(fieldFile,  sizeof(fieldFile),  "results/%s%05d_0.dat", field_name, ti);

    PetscPrintf(PETSC_COMM_WORLD,"   coordsFile = '%s'\n", coordsFile);
    PetscPrintf(PETSC_COMM_WORLD,"   fieldFile  = '%s'\n", fieldFile);

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "PostprocessMain - coords file: '%s'\n", coordsFile);
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "PostprocessMain - field file: '%s'\n",  fieldFile);

    /* 4) Read coords data into coordsArray. 
          readDataFile => length => Ncoords. */
    PetscPrintf(PETSC_COMM_WORLD, "STEP 4: Reading coordinates from file...\n");
    ierr = ReadDataFileToArray(coordsFile, &coordsArray, &Ncoords, PETSC_COMM_WORLD);
    if (ierr) {
        PetscPrintf(PETSC_COMM_WORLD,"[ERROR] PostprocessMain - Failed to read coords from '%s'.\n", coordsFile);
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "PostprocessMain - Failed to read coords from '%s'.\n", coordsFile);
        goto finalize;
    }
    PetscPrintf(PETSC_COMM_WORLD, "   ... Successfully read %d coordinate values.\n", Ncoords);

    /* Read scalar data => scalarArray => length => Nscalars. */
    PetscPrintf(PETSC_COMM_WORLD, "STEP 4.1: Reading scalar field from file...\n");
    ierr = ReadDataFileToArray(fieldFile, &scalarArray, &Nscalars, PETSC_COMM_WORLD);
    if (ierr) {
        PetscPrintf(PETSC_COMM_WORLD,"[ERROR] PostprocessMain - Failed to read field from '%s'.\n", fieldFile);
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "PostprocessMain - Failed to read field from '%s'.\n", fieldFile);
        goto finalize;
    }
    PetscPrintf(PETSC_COMM_WORLD, "   ... Successfully read %d scalar values.\n", Nscalars);

    /* 5) Decide if we want .vts (structured) or .vtp (polydata) 
          based on out_ext. 
          If out_ext = "vtp", we set meta.fileType=VTK_POLYDATA, else VTK_STRUCTURED. */
    PetscPrintf(PETSC_COMM_WORLD, "STEP 5: Determining output file type from '%s'...\n", out_ext);
    PetscBool isVTP = PETSC_FALSE;
    {
        PetscBool match;
        ierr = PetscStrcmp(out_ext, "vtp", &match);CHKERRQ(ierr);
        if (match) isVTP = PETSC_TRUE;
    }
    PetscPrintf(PETSC_COMM_WORLD, "   isVTP? => %s\n", isVTP ? "YES (will write .vtp)" : "NO (will write .vts)");

    /* 6) Populate a VTKMetaData struct. 
          - If .vts => domain is 4x4x4 => 3*3*3 = 27 nodes (example).
          - If .vtp => interpret coords as (3*N) for N points.
    */
    PetscPrintf(PETSC_COMM_WORLD, "STEP 6: Filling VTKMetaData...\n");
    VTKMetaData meta;
    memset(&meta, 0, sizeof(meta));

    if (!isVTP) {
        /* => VTK_STRUCTURED => .vts */
        meta.fileType          = VTK_STRUCTURED;
        meta.mx                = 4; /* example dimension */
        meta.my                = 4;
        meta.mz                = 4;
        meta.nnodes            = (meta.mx - 1)*(meta.my - 1)*(meta.mz - 1);

        /* coords => length=3*nnodes, scalar => nnodes */
        meta.coords            = coordsArray;
        meta.scalarField       = scalarArray;
        meta.scalarFieldName   = field_name;
        meta.numScalarFields   = 1;

        if (!rank) {
            PetscPrintf(PETSC_COMM_SELF,
                        "[postprocess] Interpreting data as structured grid => .vts\n"
                        "  Domain = %dx%dx%d => nnodes = %d\n",
                        meta.mx, meta.my, meta.mz, meta.nnodes);
        }
    }
    else {
        /* => VTK_POLYDATA => .vtp */
        meta.fileType          = VTK_POLYDATA;

        if (Ncoords % 3 != 0) {
            PetscPrintf(PETSC_COMM_WORLD, 
                        "[ERROR] coords array length %d is not multiple of 3 for .vtp.\n", 
                        Ncoords);
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "PostprocessMain - coords array length %d not multiple of 3.\n",
                      Ncoords);
            goto finalize;
        }
        meta.npoints           = Ncoords / 3;
        meta.coords            = coordsArray;
        meta.scalarField       = scalarArray;
        meta.scalarFieldName   = field_name;
        meta.numScalarFields   = 1;

        if (!rank) {
            PetscPrintf(PETSC_COMM_SELF,
                        "[postprocess] Interpreting data as polydata => .vtp\n"
                        "  npoints=%d\n", meta.npoints);
        }

        /* Create connectivity & offsets => 1 cell per point. */
        meta.connectivity = (int*)calloc(meta.npoints, sizeof(int));
        meta.offsets      = (int*)calloc(meta.npoints, sizeof(int));
        for (int i = 0; i < meta.npoints; i++) {
            meta.connectivity[i] = i;
            meta.offsets[i]      = i+1;
        }
    }

    /* 7) Construct the output filename, e.g. "results/velocity00010.vts" or ".vtp". */
    PetscPrintf(PETSC_COMM_WORLD, "STEP 7: Constructing output filename...\n");
    char outFile[256];
    PetscSNPrintf(outFile, sizeof(outFile),
                  "results/%s%05d.%s",
                  field_name, ti, out_ext);

    LOG_ALLOW(GLOBAL, LOG_INFO,
              "PostprocessMain - Will write to '%s'.\n", outFile);
    PetscPrintf(PETSC_COMM_WORLD, "   Will write output to file: %s\n", outFile);

    /* 8) Call CreateVTKFileFromMetadata to generate the file. */
    PetscPrintf(PETSC_COMM_WORLD, "STEP 8: Creating VTK file from metadata...\n");
    ierr = CreateVTKFileFromMetadata(outFile, &meta, PETSC_COMM_WORLD);
    if (ierr) {
        PetscPrintf(PETSC_COMM_WORLD,
                    "[ERROR] CreateVTKFileFromMetadata returned error %d for '%s'.\n",
                    (int)ierr, outFile);
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "PostprocessMain - CreateVTKFileFromMetadata failed (err=%d).\n",
                  (int)ierr);
        goto finalize;
    }

    if (!rank) {
        PetscPrintf(PETSC_COMM_SELF,
                    "[postprocess] Successfully wrote file: %s\n", outFile);
    }

finalize:
    /* 9) Cleanup
          free the arrays allocated by ReadDataFileToArray
          plus connectivity/offsets if we did polydata
    */
    PetscPrintf(PETSC_COMM_WORLD, "STEP 9: Cleaning up allocated resources...\n");
    free(coordsArray);
    free(scalarArray);
    if (meta.fileType == VTK_POLYDATA) {
        free(meta.connectivity);
        free(meta.offsets);
    }

    /* 10) Finalize PETSc */
    PetscPrintf(PETSC_COMM_WORLD, "STEP 10: Finalizing PETSc...\n");
    ierr = PetscFinalize();
    if (ierr) {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "PostprocessMain - PetscFinalize returned error=%d\n",
                  (int)ierr);
        PetscPrintf(PETSC_COMM_WORLD,
                    "[ERROR] PetscFinalize returned error=%d\n", (int)ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "PostprocessMain - Done.\n");
    PetscPrintf(PETSC_COMM_WORLD, "=== Exiting PostprocessMain ===\n");

    return 0; /* Return success code. */
}

int main(int argc, char **argv)
{
    /* If you want a separate main: simply call PostprocessMain */
    return PostprocessMain(argc, argv);
}

/* End of postprocess.c */
