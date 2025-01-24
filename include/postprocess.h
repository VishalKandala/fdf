/******************************************************************************/
/* postprocess.h                                                              */
/*                                                                            */
/* Header file for postprocessing.c. Contains prototypes for VTK output       */
/* routines, helper functions, and any additional post-processing hooks.      */
/*                                                                            */
/******************************************************************************/

#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscdmcomposite.h>

#include "logging.h"
#include "common.h"
#include "io.h"
#include "grid.h"
#include "interpolation.h"
#include "ParticleSwarm.h"
#include <petscsys.h>   /* For PETSc types, e.g. PetscErrorCode, etc. */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Gathers a parallel PETSc Vec into a sequential Vec on rank 0.
 * 
 * @param[in]  in    The parallel (or possibly serial) input Vec.
 * @param[out] out   The sequential Vec on rank 0 (created in this function).
 * @param[in]  comm  The MPI communicator for the Vec.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode GatherParallelVecToSeq(Vec in, Vec *out, MPI_Comm comm);

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
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode WriteVTKXMLHeader(FILE *fp,
                                 PetscInt mx, PetscInt my, PetscInt mz,
                                 PetscInt nnodes,
                                 PetscInt boffset,
                                 const char *field_name,
                                 PetscInt *boffsetOut);

/**
 * @brief Writes a single block of appended binary data to the .vts file.
 *
 * @param[in] fp      The file pointer (rank 0).
 * @param[in] buf     Pointer to the data (already in correct type layout).
 * @param[in] nvals   Number of values in the block (e.g., 3*nnodes).
 * @param[in] dsize   Size in bytes per value (e.g., sizeof(double)).
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode WriteVTKAppendedBlock(FILE *fp,
                                     const PetscScalar *buf,
                                     PetscInt nvals,
                                     PetscInt dsize);

/**
 * @brief Writes the closing tags for the VTK XML file.
 *
 * @param[in] fp  The file pointer (rank 0).
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode WriteVTKXMLFooter(FILE *fp);

/**
 * @brief Creates a VTK .vts file from an existing data file.
 *
 * 1) Reads scalar field (via ReadFieldData).
 * 2) Gathers field & coordinates onto rank 0.
 * 3) Writes a .vts file in appended binary format.
 *
 * @param[in] user       Pointer to UserCtx (must hold a valid DM).
 * @param[in] field_name Field name (e.g., "myField").
 * @param[in] ti         Time index for constructing filenames.
 * @param[in] read_ext   Input file extension (e.g., "dat").
 * @param[in] out_ext    Output file extension (e.g., "vts").
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode CreateVTKFileFromData(UserCtx    *user,
                                     const char *field_name,
                                     PetscInt    ti,
                                     const char *read_ext,
                                     const char *out_ext);

/**
 * @brief Placeholder function for creating a VTK .vts file with multiple fields.
 *
 * @param[in] user  Pointer to UserCtx.
 * @param[in] ti    Time index.
 * @param[in] ext   File extension or format specifier.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode CreateVTKFileMultiField(UserCtx *user, PetscInt ti, const char *ext);

/**
 * @brief Another placeholder for specialized post-processing, e.g. Q-criterion or shear output.
 *
 * @param[in] user  Pointer to UserCtx.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PostProcessCustom(UserCtx *user);

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
PetscErrorCode CreateVTPFileFromSwarm(UserCtx *user, const char *field_name, PetscInt ti);

#ifdef __cplusplus
}
#endif

#endif /* POSTPROCESS_H */
