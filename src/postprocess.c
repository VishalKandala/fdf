/*****************************************************************************/
/* postprocess.c                                                          */
/*                                                                           */
/* This file contains routines for writing data to VTK files (XML/.vts       */
/* format), as well as other post-processing and output-related functions.   */
/*                                                                           */
/* Example includes:                                                         */
/*  - GatherParallelVecToSeq() for gathering PETSc Vecs to rank 0            */
/*  - WriteVTKXMLHeader(), WriteVTKXMLFooter() for generating VTK XML markup */
/*  - WriteVTKAppendedBlock() for appended binary data in VTK                */
/*  - CreateVTKFileFromData() for producing a .vts file from a scalar Vec    */
/*  - Additional placeholders for extended post-processing routines          */
/*    (e.g., multi-field outputs, vector fields, etc.).                      */
/*                                                                           */
/*****************************************************************************/

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscdmcomposite.h>

#include "postprocess.h"
//#include "logging.h"
//#include "common.h"
//#include "io.h"
//#include "grid.h"
//#include "interpolation.h"
//  #include "ParticleSwarm.h"
/* ========================== Utility Functions ========================== */

/**
 * @brief Gathers a parallel PETSc Vec into a sequential Vec on rank 0.
 *
 * This function checks if MPI size > 1. If so, it creates a new Seq Vec on rank 0
 * and scatters `in` into that sequential Vec. Otherwise, it just references `in`.
 *
 * @param[in]  in    The parallel (or possibly serial) input Vec.
 * @param[out] out   The sequential Vec on rank 0 (created in this function or references in).
 * @param[in]  comm  The MPI communicator for the Vec.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode GatherParallelVecToSeq(Vec in, Vec *out, MPI_Comm comm)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank, size;
    PetscInt       Nglobal;

    PetscFunctionBegin;
    ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm, &size);CHKERRQ(ierr);
    ierr = VecGetSize(in, &Nglobal);CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
        "GatherParallelVecToSeq - Gathering Vec of global size %d to rank 0.\n",
        Nglobal);

    if (size > 1) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "GatherParallelVecToSeq - MPI size > 1, creating seq Vec on rank 0.\n");

        ierr = VecCreateSeq(PETSC_COMM_SELF, Nglobal, out);CHKERRQ(ierr);

        IS from, to;
        VecScatter scatter;
        ierr = ISCreateStride(comm, Nglobal, 0, 1, &from);CHKERRQ(ierr);
        ierr = ISCreateStride(PETSC_COMM_SELF, Nglobal, 0, 1, &to);CHKERRQ(ierr);

        ierr = VecScatterCreate(in, from, *out, to, &scatter);CHKERRQ(ierr);
        ierr = VecScatterBegin(scatter, in, *out, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        ierr = VecScatterEnd(  scatter, in, *out, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);

        ierr = VecScatterDestroy(&scatter);CHKERRQ(ierr);
        ierr = ISDestroy(&from);CHKERRQ(ierr);
        ierr = ISDestroy(&to);CHKERRQ(ierr);

    } else {
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "GatherParallelVecToSeq - Running in serial, reusing the same Vec.\n");
        *out = in;
    }

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
        "GatherParallelVecToSeq - Finished gathering onto rank 0.\n");

    PetscFunctionReturn(0);
}

/**
 * @brief Writes the XML header (up to AppendedData) for a VTK .vts file.
 *
 * @param[in]  fp          The opened file pointer (on rank 0).
 * @param[in]  mx,my,mz    Global grid dimensions in x, y, z.
 * @param[in]  nnodes      Number of interior nodes (e.g. (mx-1)*(my-1)*(mz-1)).
 * @param[in]  boffset     Current byte offset for appended data blocks.
 * @param[out] boffsetOut  Updated byte offset after writing the tags.
 * @param[in]  field_name  Name of the scalar field.
 *
 * This writes `<VTKFile>`, `<StructuredGrid>`, `<Points>`, `<PointData>` tags,
 * plus the DataArray placeholders with `offset`.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteVTKXMLHeader(FILE *fp,
                                 PetscInt mx, PetscInt my, PetscInt mz,
                                 PetscInt nnodes,
                                 PetscInt boffset,
                                 const char *field_name,
                                 PetscInt *boffsetOut)
{
    PetscFunctionBegin;

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
        "WriteVTKXMLHeader - Writing header for field '%s', offset=%d\n",
        field_name, boffset);

    const char *byte_order = "LittleEndian";
    const char precision[] = "Float64";

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n", byte_order);
    fprintf(fp, "  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0,mx-2, 0,my-2, 0,mz-2);
    fprintf(fp, "    <Piece Extent=\"%d %d %d %d %d %d\">\n", 0,mx-2, 0,my-2, 0,mz-2);

    /* Points */
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\" />\n",
            precision, boffset);
    boffset += (PetscInt)sizeof(int) + 3*nnodes*(PetscInt)sizeof(PetscScalar);
    fprintf(fp, "      </Points>\n");

    /* Single scalar field in <PointData> */
    fprintf(fp, "      <PointData Scalars=\"%s\">\n", field_name);
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\" />\n",
            precision, field_name, boffset);
    boffset += (PetscInt)sizeof(int) + nnodes*(PetscInt)sizeof(PetscScalar);
    fprintf(fp, "      </PointData>\n");

    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </StructuredGrid>\n");
    fprintf(fp, "  <AppendedData encoding=\"raw\">\n");
    fprintf(fp, "_");

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
        "WriteVTKXMLHeader - New boffset after header: %d\n", boffset);

    *boffsetOut = boffset;
    PetscFunctionReturn(0);
}

/**
 * @brief Writes a single block of appended binary data to the file.
 *
 * @param[in] fp      The file pointer (rank 0).
 * @param[in] buf     Pointer to the data (already in correct type layout).
 * @param[in] nvals   Number of values in the block (e.g., 3*nnodes for coords).
 * @param[in] dsize   Size in bytes per value (e.g., size of double).
 *
 * Writes:
 *   1) An integer "blockSize" = nvals*dsize
 *   2) The raw data.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteVTKAppendedBlock(FILE *fp,
                                     const PetscScalar *buf,
                                     PetscInt nvals,
                                     PetscInt dsize)
{
    PetscFunctionBegin;
    PetscInt blockSize = nvals * dsize;

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
       "WriteVTKAppendedBlock - Writing block of %d values, blockSize=%d bytes.\n",
       nvals, blockSize);

    fwrite(&blockSize, sizeof(int), 1, fp);
    fwrite(buf, dsize, (size_t)nvals, fp);

    PetscFunctionReturn(0);
}

/**
 * @brief Writes the closing tags for the VTK XML file.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteVTKXMLFooter(FILE *fp)
{
    PetscFunctionBegin;
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
        "WriteVTKXMLFooter - Writing closing XML tags.\n");

    fprintf(fp, "\n </AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");

    PetscFunctionReturn(0);
}

/* ===================== Main Refactored Function ===================== */

/**
 * @brief Creates a VTK .vts file from an existing data file (refactored + logging).
 *
 * 1) Reads scalar field from file (via ReadFieldData).
 * 2) Gathers field & coordinates on rank 0.
 * 3) Writes a .vts file with appended binary data.
 *
 * @param[in] user       The UserCtx (must have user->da for domain & coords).
 * @param[in] field_name Name of the field (e.g. "myField").
 * @param[in] ti         Time index for reading/writing files.
 * @param[in] read_ext   File extension for the input data (e.g., "dat").
 * @param[in] out_ext    File extension for the output (e.g., "vts").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode CreateVTKFileFromData(UserCtx    *user,
                                     const char *field_name,
                                     PetscInt    ti,
                                     const char *read_ext,
                                     const char *out_ext)
{
    PetscErrorCode ierr;
    MPI_Comm       comm;
    PetscMPIInt    rank, size;
    PetscInt       mx, my, mz, nnodes;
    PetscInt       dsize;
    char           filen[256];
    FILE          *fp           = NULL;
    Vec            fieldVec      = NULL;
    Vec            fieldVecSeq   = NULL; /* seq on rank 0 */
    Vec            coordsGlobal  = NULL; /* parallel coords */
    Vec            coordsSeq     = NULL; /* seq coords on rank 0 */
    PetscScalar   *coords3d      = NULL; /* pointer to store 3*(mx-1)*(my-1)*(mz-1) coords */
    PetscScalar   *fieldArray    = NULL; /* pointer to store (mx-1)*(my-1)*(mz-1) for the field. */

    PetscFunctionBegin;
    comm = PetscObjectComm((PetscObject)user->da);
    ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm, &size);CHKERRQ(ierr);
    ierr = MPI_Type_size(MPI_DOUBLE, &dsize);CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO,
        "CreateVTKFileFromData - Start creating VTK for field '%s' at time %d.\n",
        field_name, ti);

    /* 1) Domain size from the DMDA. */
    {
        DMDALocalInfo info;
        ierr = DMDAGetLocalInfo(user->da, &info);CHKERRQ(ierr);
        mx = info.mx;  my = info.my;  mz = info.mz;
    }
    nnodes = (mx-1)*(my-1)*(mz-1);

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
        "CreateVTKFileFromData - DM global size: mx=%d my=%d mz=%d => nnodes=%d.\n",
        mx, my, mz, nnodes);

    /* 2) Create an empty parallel Vec, read from file into it. */
    ierr = VecCreate(comm, &fieldVec);CHKERRQ(ierr);
    ierr = VecSetFromOptions(fieldVec);CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO,
      "CreateVTKFileFromData - Reading field '%s' from time %d with ext='%s'.\n",
      field_name, ti, read_ext);

    /* 3) Read from file into fieldVec (via your existing routine). */
    ierr = ReadFieldData(user, field_name, fieldVec, ti, read_ext);CHKERRQ(ierr);

    /* 4) Gather that Vec onto rank 0. */
    ierr = GatherParallelVecToSeq(fieldVec, &fieldVecSeq, comm);CHKERRQ(ierr);

    /* 5) Gather coordinates similarly. */
    ierr = DMGetCoordinates(user->da, &coordsGlobal);CHKERRQ(ierr);
    ierr = GatherParallelVecToSeq(coordsGlobal, &coordsSeq, comm);CHKERRQ(ierr);

    /* 6) Allocate arrays on rank 0 for the interior data. */
    if (!rank) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "CreateVTKFileFromData - Allocating memory for coords3d & fieldArray.\n");
        ierr = PetscMalloc1(3*nnodes, &coords3d);CHKERRQ(ierr);
        ierr = PetscMalloc1(nnodes,   &fieldArray);CHKERRQ(ierr);
    }

    /* 7) Copy coordinate data (seq) into coords3d for the interior. */
    if (!rank) {
        PetscScalar *seqA;
        PetscInt    i, j, k, idx=0;
        PetscInt    ncoor;

        ierr = VecGetSize(coordsSeq, &ncoor);CHKERRQ(ierr);
        ierr = VecGetArray(coordsSeq, &seqA);CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "CreateVTKFileFromData - Copying coordinate data (size=%d) to interior array.\n",
            ncoor);

        idx=0;
        for (k=0; k<mz-1; k++) {
            for (j=0; j<my-1; j++) {
                for (i=0; i<mx-1; i++) {
                    PetscInt nloc = ((k*my)+j)*mx + i;
                    coords3d[idx*3 + 0] = seqA[nloc*3 + 0];
                    coords3d[idx*3 + 1] = seqA[nloc*3 + 1];
                    coords3d[idx*3 + 2] = seqA[nloc*3 + 2];
                    idx++;
                }
            }
        }
        ierr = VecRestoreArray(coordsSeq, &seqA);CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "CreateVTKFileFromData - Finished copying coordinates.\n");
    }

    /* 8) Copy field data (seq) into fieldArray for the interior. */
    if (!rank) {
        PetscScalar *seqF;
        PetscInt    i, j, k, idx=0;
        PetscInt    Nf;

        ierr = VecGetSize(fieldVecSeq, &Nf);CHKERRQ(ierr);
        ierr = VecGetArray(fieldVecSeq, &seqF);CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "CreateVTKFileFromData - Copying field data (size=%d) to interior array.\n",
            Nf);

        idx=0;
        for (k=0; k<mz-1; k++) {
            for (j=0; j<my-1; j++) {
                for (i=0; i<mx-1; i++) {
                    PetscInt nloc = ((k*my)+j)*mx + i;
                    fieldArray[idx++] = seqF[nloc];
                }
            }
        }
        ierr = VecRestoreArray(fieldVecSeq, &seqF);CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "CreateVTKFileFromData - Finished copying scalar field.\n");
    }

    /* 9) Write the .vts file on rank 0. */
    if (!rank) {
        PetscInt boffset=0;
        PetscSNPrintf(filen, sizeof(filen), "results/%s%05d.%s", field_name, ti, out_ext);

        LOG_ALLOW(GLOBAL, LOG_INFO,
          "CreateVTKFileFromData - Opening file '%s' for writing.\n", filen);

        FILE *fpOut = fopen(filen, "wb");
        if (!fpOut) {
            char err_msg[512];
            PetscSNPrintf(err_msg, sizeof(err_msg),
                "CreateVTKFileFromData - Could not open file '%s' for writing.", filen);
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, err_msg);
        }

        /* Write XML header. */
        ierr = WriteVTKXMLHeader(fpOut, mx, my, mz, nnodes, boffset, field_name, &boffset);CHKERRQ(ierr);

        /* coords3d => 3*nnodes */
        ierr = WriteVTKAppendedBlock(fpOut, coords3d, 3*nnodes, dsize);CHKERRQ(ierr);

        /* fieldArray => nnodes */
        ierr = WriteVTKAppendedBlock(fpOut, fieldArray, nnodes, dsize);CHKERRQ(ierr);

        /* Finish the file. */
        ierr = WriteVTKXMLFooter(fpOut);CHKERRQ(ierr);
        fclose(fpOut);

        PetscFree(coords3d);
        PetscFree(fieldArray);

        LOG_ALLOW(GLOBAL, LOG_INFO,
          "CreateVTKFileFromData - Successfully wrote file: %s\n", filen);
    }

    /* 10) Cleanup: destroy everything we created. */
    if (size>1 && fieldVecSeq!=fieldVec) {
        ierr = VecDestroy(&fieldVecSeq);CHKERRQ(ierr);
    }
    if (size>1 && coordsSeq!=coordsGlobal) {
        ierr = VecDestroy(&coordsSeq);CHKERRQ(ierr);
    }
    ierr = VecDestroy(&fieldVec);CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO,
      "CreateVTKFileFromData - Done.\n");

    PetscFunctionReturn(0);
}


/* =================== Additional Post-Processing Routines =================== */

/**
 * @brief Example placeholder for future expansions, e.g. writing multi-field .vts.
 */
PetscErrorCode CreateVTKFileMultiField(UserCtx *user, PetscInt ti, const char *ext)
{
    /* Implementation of a multi-field approach:
       - Possibly gather velocity (3-component) + pressure + etc.
       - Write multiple <DataArray> entries in <PointData>.
    */
    PetscFunctionBegin;
    LOG_ALLOW(GLOBAL, LOG_INFO,
      "CreateVTKFileMultiField - Not yet implemented.\n");
    PetscFunctionReturn(0);
}

/**
 * @brief Another placeholder for specialized post-processing, e.g. Q-criterion or shear output.
 */
PetscErrorCode PostProcessCustom(UserCtx *user)
{
    /* Implementation for advanced post-processing:
       - e.g., computing Q-criterion, or doing time averaging, etc.
    */
    PetscFunctionBegin;
    LOG_ALLOW(GLOBAL, LOG_INFO,
      "PostProcessCustom - Not yet implemented.\n");
    PetscFunctionReturn(0);
}

/**
 * @brief Creates a VTK .vtp file from swarm particle data (positions and velocities).
 *
 * 1) Retrieves particle data (positions and velocities) from the swarm.
 * 2) Gathers particle information globally and writes a .vtp file on rank 0.
 * 3) Includes particle positions, velocities, and vertex connectivity in the .vtp file.
 *
 * @param[in] user       The UserCtx containing the particle swarm (user->swarm).
 * @param[in] field_name Name of the field for naming the output file (e.g., "particles").
 * @param[in] ti         Time index for the file (used to create unique filenames).
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode CreateVTPFileFromSwarm(UserCtx *user, const char *field_name, PetscInt ti)
{
    PetscErrorCode ierr;        // Error code returned by PETSc functions
    MPI_Comm comm;              // MPI communicator
    PetscMPIInt rank, size;     // MPI rank and size
    PetscInt n_local, n_global; // Local and global number of particles
    DM swarm = user->swarm;     // Particle swarm from user context
    FILE *fp = NULL;            // File pointer for output file
    char filen[256];            // Buffer for the output filename
    PetscReal *positions = NULL; // Pointer to particle positions
    PetscReal *velocities = NULL; // Pointer to particle velocities

    PetscFunctionBegin;

    // Get MPI communicator
    comm = PETSC_COMM_WORLD;

    // Get the rank of the current process (used for determining file output responsibility)
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);

    // Get the total number of processes in the communicator
    ierr = MPI_Comm_size(comm, &size); CHKERRQ(ierr);

    // Get the number of particles managed by the current process
    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);

    // Get the total number of particles across all processes
    ierr = DMSwarmGetSize(swarm, &n_global); CHKERRQ(ierr);

    // Access the "position" field from the swarm (each particle has 3 components: x, y, z)
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);

    // Access the "velocity" field from the swarm (each particle has 3 components: vx, vy, vz)
    ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);

    // Only the process with rank 0 writes the output file
    if (rank == 0) {
        // Generate the output filename using the field name and time index
        snprintf(filen, sizeof(filen), "results/%s%05d.vtp", field_name, ti);

        // Open the file for writing
        fp = fopen(filen, "w");
        if (!fp) SETERRQ(comm, PETSC_ERR_FILE_OPEN, "Failed to open VTP file for writing");

        // Write the VTP XML header
        fprintf(fp, "<?xml version=\"1.0\"?>\n");
        fprintf(fp, "<VTKFile type=\"PolyData\" version=\"0.1\">\n");
        fprintf(fp, "  <PolyData>\n");
        fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"%d\">\n", 
                n_global, n_global);

        // Write particle positions (Points section)
        fprintf(fp, "      <Points>\n");
        fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for (PetscInt i = 0; i < n_local; i++) {
            // Each line corresponds to one particle's position: x, y, z
            fprintf(fp, "          %g %g %g\n", 
                    positions[3*i], positions[3*i+1], positions[3*i+2]);
        }
        fprintf(fp, "        </DataArray>\n");
        fprintf(fp, "      </Points>\n");

        // Write particle velocities as point data (PointData section)
        fprintf(fp, "      <PointData Scalars=\"Velocity\">\n");
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for (PetscInt i = 0; i < n_local; i++) {
            // Each line corresponds to one particle's velocity: vx, vy, vz
            fprintf(fp, "          %g %g %g\n",
                    velocities[3*i], velocities[3*i+1], velocities[3*i+2]);
        }
        fprintf(fp, "        </DataArray>\n");
        fprintf(fp, "      </PointData>\n");

        // Write vertices (Verts section, where each vertex corresponds to a particle)
        fprintf(fp, "      <Verts>\n");
        fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
        for (PetscInt i = 0; i < n_global; i++) {
            // Each particle is connected to its index
            fprintf(fp, "          %d\n", i);
        }
        fprintf(fp, "        </DataArray>\n");
        fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
        for (PetscInt i = 0; i < n_global; i++) {
            // Offset for each vertex
            fprintf(fp, "          %d\n", i+1);
        }
        fprintf(fp, "        </DataArray>\n");
        fprintf(fp, "      </Verts>\n");

        // Close the XML tags for the VTK file
        fprintf(fp, "    </Piece>\n");
        fprintf(fp, "  </PolyData>\n");
        fprintf(fp, "</VTKFile>\n");

        // Close the file
        fclose(fp);
    }

    // Restore the "position" field after use (cleanup)
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);

    // Restore the "velocity" field after use (cleanup)
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);

    // Mark the end of the function and return a success code
    PetscFunctionReturn(0);
}


static const char help[] = "Reads a .dat field and writes a .vts file for post-processing.\n";

/* Main function for 'postprocess' executable */
int main(int argc, char **argv)
{
    UserCtx     *user = NULL;     // User context 
    PetscMPIInt rank, size;       // MPI rank and size
    PetscInt    ti = 0;           // Time index
    PetscInt    np = 0;           // Number of particles
    PetscInt    block_number = 1; // Number of blocks
    PetscBool   rstart = PETSC_FALSE;
    PetscBool   is_particle_data = PETSC_FALSE;
    PetscErrorCode ierr;
    BoundingBox *bboxlist;    // Array of bounding boxes


    /* Default options */
    char field_name[64] = "velocity"; 
    char read_ext[16] = "dat";    
    char out_ext[16] = "vts";     

    /* 1) Initialize PETSc */
    ierr = PetscInitialize(&argc, &argv, NULL, help);CHKERRQ(ierr);

    /* 2) Get command line options */
    ierr = PetscOptionsGetInt(NULL, NULL, "-ti", &ti, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL, NULL, "-field", field_name, sizeof(field_name), NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL, NULL, "-read_ext", read_ext, sizeof(read_ext), NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL, NULL, "-out_ext", out_ext, sizeof(out_ext), NULL);CHKERRQ(ierr);
    
    // Check if this is particle data
    ierr = PetscStrcmp(out_ext, "vtp", &is_particle_data);CHKERRQ(ierr);

    /* 3) Initialize simulation context */
    ierr = InitializeSimulation(&user, &rank, &size, &np, &rstart, &ti, &block_number);CHKERRQ(ierr);
    
    // Print the value of is_particle_data
    ierr = PetscPrintf(PETSC_COMM_WORLD, "is_particle_data: %s\n", is_particle_data ? "PETSC_TRUE" : "PETSC_FALSE");CHKERRQ(ierr);

   // Setup the computational grid
    ierr = SetupGridAndVectors(user, block_number); CHKERRQ(ierr);

    // Gather bounding boxes on rank 0
    ierr = GatherAllBoundingBoxes(user, &bboxlist); CHKERRQ(ierr);

    // Broadcast bboxlist to all ranks
    ierr = BroadcastAllBoundingBoxes(user, &bboxlist); CHKERRQ(ierr);

    // Perform particle swarm operations with bboxlist knowledge on all ranks
    ierr = CreateParticleSwarm(user, np, bboxlist); CHKERRQ(ierr);

    /* 4) Setup grid/vectors or particle swarm based on data type */
    if (is_particle_data) {
        // Initialize swarm for particle data
        InitializeSwarm(user);
        // Setup particle fields (position, velocity)

        ierr = DMSwarmRegisterPetscDatatypeField(user->swarm, "position", 3, PETSC_REAL);CHKERRQ(ierr);
        ierr = DMSwarmRegisterPetscDatatypeField(user->swarm, "velocity", 3, PETSC_REAL);CHKERRQ(ierr);
        ierr = DMSwarmFinalizeFieldRegister(user->swarm);CHKERRQ(ierr);
        
        LOG_ALLOW(GLOBAL, LOG_INFO, "postprocess: Creating VTP file for particle data...\n");
        ierr = CreateVTPFileFromSwarm(user, field_name, ti);CHKERRQ(ierr);
    } else {
        // Regular grid data
        ierr = SetupGridAndVectors(user, block_number);CHKERRQ(ierr);
        
        LOG_ALLOW(GLOBAL, LOG_INFO, "postprocess: Creating VTK file from grid data...\n");
        ierr = CreateVTKFileFromData(user, field_name, ti, read_ext, out_ext);CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "postprocess: Successfully wrote output file.\n");

    /* 5) Cleanup */
    if (is_particle_data) {
        ierr = DMDestroy(&user->swarm);CHKERRQ(ierr);
    }
    ierr = PetscFinalize();CHKERRQ(ierr);
    return 0;
}



/* End of postprocessing.c */
