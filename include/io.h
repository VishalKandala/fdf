/**
 * @file io.h
 * @brief Header file for input/output operations in the simulation.
 *
 * This file declares functions for reading and writing simulation fields, grid generation parameters,
 * and optional statistical, LES, and RANS data.
 */

#ifndef IO_H
#define IO_H

#include <petsc.h> // System dependency, required.
#include "common.h" // For shared definitions like UserCtx.
#include "logging.h"        // For logging macros.
#include "grid.h"           // Only if grid-related functions are directly used


// ------------------------ Function Declarations ------------------------

/**
 * @brief Reads grid generation parameters from input options or configuration.
 *
 * This function reads domain dimensions and block parameters for programmatically generating a grid.
 *
 * @param[in,out] user      Pointer to the UserCtx structure containing grid details.
 * @param[out] grid1d       Pointer to flag indicating if the grid is 1D (1) or 3D (0).
 * @param[out]    xMax          Pointer to store the upper bound in x-direction.
 * @param[out]    yMin          Pointer to store the lower bound in y-direction.
 * @param[out]    yMax          Pointer to store the upper bound in y-direction.
 * @param[out]    zMin          Pointer to store the lower bound in z-direction.
 * @param[out]    zMax          Pointer to store the upper bound in z-direction.
 * @param[out] imm          Pointer to array storing the i-dimensions of each block.
 * @param[out] jmm          Pointer to array storing the j-dimensions of each block.
 * @param[out] kmm          Pointer to array storing the k-dimensions of each block.
 * @param[out] nblk         Pointer to number of blocks.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadGridGenerationInputs(UserCtx *user, PetscInt *grid1d, PetscReal *xMin, PetscReal *xMax,
                               PetscReal *yMin, PetscReal *yMax, PetscReal *zMin, PetscReal *zMax,  PetscInt **imm, PetscInt **jmm, PetscInt **kmm,
					PetscInt *nblk);
/**
 * @brief Reads grid parameters from a file.
 *
 * This function retrieves grid block dimensions, number of blocks, and grid type (1D or 3D)
 * from a file and broadcasts the data to all MPI processes.
 *
 * @param[in]  filename  Name of the file containing grid data.
 * @param[out] nblk      Pointer to store the number of blocks.
 * @param[out] imm       Pointer to Array to store the i-dimensions for each block.
 * @param[out] jmm       Pointer to Array to store the j-dimensions for each block.
 * @param[out] kmm       Pointer to Array to store the k-dimensions for each block.
 * @param[out] grid1d    Pointer to store the grid type (1 for 1D, 0 for 3D).
 * @param[in]  comm      MPI communicator for broadcasting data.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode ReadGridFile(const char *filename, PetscInt *nblk, PetscInt **imm, PetscInt **jmm, 
                            PetscInt **kmm, PetscInt *grid1d, MPI_Comm comm);

/**
 * @brief Reads binary field data for velocity, pressure, and other required vectors.
 *
 * Reads contravariant velocity (`Ucont`) from `vfield`, Cartesian velocity (`Ucat`) from `ufield`,
 * pressure (`P`), node state (`Nvert_o`), and optionally statistical quantities, LES, and RANS fields
 * from binary files. Logs missing files but continues execution.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the simulation context.
 * @param[in]     The timestep at which the simulation field data needs to be read.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadSimulationFields(UserCtx *user, PetscInt ti);


/**
 * @brief Reads data for a specific field from a file into the provided vector.
 *
 * This function uses the field name to construct the file path and reads the data
 * from the corresponding file into the provided PETSc vector.
 *
 * @param[in]     user        Pointer to the UserCtx structure containing simulation context.
 * @param[in]     field_name  Name of the field (e.g., "ufield", "vfield", "pfield").
 * @param[out]    field_vec   PETSc vector to store the field data.
 * @param[in]     ti          Time index for constructing the file name.
 * @param[in]     ext         File extension (e.g., "dat").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadFieldData(UserCtx *user, const char *field_name, Vec field_vec, PetscInt ti, const char *ext);


/**
 * @brief Reads statistical fields used for time-averaged simulations.
 *
 * Reads statistical quantities such as velocity sums and pressure sums. Logs missing files and initializes
 * fields to zero if files are unavailable.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the simulation context.
 * @param[in]     The timestep at which the simulation field data needs to be read.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadStatisticalFields(UserCtx *user,PetscInt ti);

/**
 * @brief Reads LES-related fields used in turbulence modeling.
 *
 * Reads the Smagorinsky constant (`Cs`) and transfers data to local vectors. Logs missing files and
 * initializes fields to zero if files are unavailable.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the simulation context.
 * @param[in]     ti          Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadLESFields(UserCtx *user, PetscInt ti);

/**
 * @brief Reads RANS-related fields for turbulence modeling.
 *
 * Reads `K_Omega` fields (used in RANS modeling) and initializes them if files are unavailable.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the simulation context.
 * @param[in]     ti          Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadRANSFields(UserCtx *user, PetscInt ti);

/**
 * @brief Writes data from a specific PETSc vector to a file.
 *
 * This function uses the field name to construct the file path and writes the data
 * from the provided PETSc vector to the corresponding file.
 *
 * @param[in] user       Pointer to the UserCtx structure containing simulation context.
 * @param[in] field_name Name of the field (e.g., "ufield", "vfield", "pfield").
 * @param[in] field_vec  PETSc vector containing the field data to write.
 * @param[in] ti         Time index for constructing the file name.
 * @param[in] ext        File extension (e.g., "dat").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteFieldData(UserCtx *user, const char *field_name, Vec field_vec, PetscInt ti, const char *ext);


/**
 * @brief Writes simulation fields to files.
 *
 * This function writes contravariant velocity, Cartesian velocity, pressure, and node state
 * fields to their respective binary files. It also conditionally writes LES, RANS, and
 * statistical fields if they are enabled.
 *
 * @param[in] user Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteSimulationFields(UserCtx *user);


/**
 * @brief Writes statistical fields for averaging purposes.
 *
 * This function writes data for fields such as Ucat_sum, Ucat_cross_sum, Ucat_square_sum,
 * and P_sum to their respective binary files.
 *
 * @param[in] user Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteStatisticalFields(UserCtx *user);

/**
 * @brief Writes LES-related fields.
 *
 * This function writes LES-related fields such as Cs (Smagorinsky constant)
 * to their respective binary files.
 *
 * @param[in] user Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteLESFields(UserCtx *user);

/**
 * @brief Writes RANS-related fields.
 *
 * This function writes RANS-related fields such as K_Omega to their respective
 * binary files.
 *
 * @param[in] user Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode WriteRANSFields(UserCtx *user);

/**
 * @brief Writes data from a specific field in a PETSc Swarm to a file.
 *
 * This function retrieves the Swarm from the UserCtx (i.e., `user->swarm`) and
 * creates a global PETSc vector from the specified Swarm field. It then calls
 * the existing WriteFieldData() function to handle the actual I/O operation.
 * After writing the data, the function destroys the temporary global vector 
 * to avoid memory leaks.
 *
 * @param[in] user       Pointer to the UserCtx structure containing simulation context
 *                       and the PetscSwarm (as `user->swarm`).
 * @param[in] field_name Name of the Swarm field to be written (e.g., "my_field").
 * @param[in] ti         Time index used to construct the output file name.
 * @param[in] ext        File extension (e.g., "dat", "bin").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 *
 * @note Compatible with PETSc 3.14.4.
 */
PetscErrorCode WriteSwarmField(UserCtx *user, const char *field_name, PetscInt ti, const char *ext);

/* --------------------------------------------------------------------
   ReadDataFileToArray

   Reads a simple ASCII .dat file containing one numeric value per line.

   PARAMETERS:
     filename   - (input)  path to the .dat file
     data_out   - (output) pointer to newly allocated array (rank 0 
                          broadcasts to all ranks so each rank gets a copy)
     Nout       - (output) number of values read; same on all ranks
     comm       - (input)  MPI communicator (we use rank 0 for I/O)

   RETURN:
     PetscErrorCode (0 = success, or error code on failure)

   NOTES:
     1) This function checks file existence on rank 0 and broadcasts
        the result. If the file doesn't exist, it sets an error.
     2) The array is allocated on all ranks, and the data is
        broadcast so each rank has a local copy.
     3) If your file format is more complex (e.g. multiple columns,
        binary, etc.), you can adapt the reading logic here.
-------------------------------------------------------------------- */
PetscInt ReadDataFileToArray(const char   *filename,
                        double      **data_out,
                        PetscInt          *Nout,
                        MPI_Comm      comm);

/* --------------------------------------------------------------------
   CreateVTKFileFromMetadata

   Creates a .vts (if fileType=VTK_STRUCTURED) or .vtp (if fileType=VTK_POLYDATA),
   writing out the data contained in the VTKMetaData structure.

   PARAMETERS:
     filename   - (input) path to the output file, e.g. "myfile.vts" or "myfile.vtp"
     meta       - (input) pointer to the VTKMetaData struct with all geometry/fields
     comm       - (input) MPI communicator (only rank 0 writes the file)

   RETURN:
     PetscErrorCode (0 = success, or error code on failure)

   NOTES:
     1) The function writes an XML header, the appended binary data blocks
        (coordinates, scalar field, connectivity if needed), and the closing tags.
     2) If meta->fileType = VTK_STRUCTURED, it calls .vts-specific logic.
        If meta->fileType = VTK_POLYDATA, it calls .vtp logic.
     3) Coordinates and field arrays must be properly sized. For instance,
        if fileType=VTK_STRUCTURED, coords => length=3*nnodes, field=>nnodes.
        If fileType=VTK_POLYDATA, coords => length=3*npoints, field=>npoints,
        plus connectivity/offsets => length=npoints each.
-------------------------------------------------------------------- */
PetscInt CreateVTKFileFromMetadata(const char       *filename,
                              const VTKMetaData *meta,
                              MPI_Comm          comm);

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
PetscErrorCode VecToArrayOnRank0(Vec inVec, PetscInt *N, double **arrayOut);

/**
 * @brief Reads data from a file into a specified field of a PETSc DMSwarm.
 *
 * This function is the counterpart to WriteSwarmField(). It creates a global PETSc vector 
 * that references the specified DMSwarm field, uses ReadFieldData() to read the data from 
 * a file, and then destroys the global vector reference.
 *
 * @param[in]  user       Pointer to the UserCtx structure (containing `user->swarm`).
 * @param[in]  field_name Name of the DMSwarm field to read into (must be previously declared/allocated).
 * @param[in]  ti         Time index used to construct the input file name.
 * @param[in]  ext        File extension (e.g., "dat" or "bin").
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 *
 * @note Compatible with PETSc 3.14.x.
 */
PetscErrorCode ReadSwarmField(UserCtx *user, const char *field_name, PetscInt ti, const char *ext);


/**
 * @brief Reads multiple fields (positions, velocity, CellID, and weight) into a DMSwarm.
 *
 * This function is analogous to ReadSimulationFields() but targets a DMSwarm. 
 * Each Swarm field is read from a separate file using ReadSwarmField().
 * 
 * @param[in,out] user Pointer to the UserCtx structure containing the DMSwarm (user->swarm).
 * @param[in]     ti   Time index for constructing the file name.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ReadAllSwarmFields(UserCtx *user, PetscInt ti);

/**
 * @brief Reads coordinate data (for particles)  from file into a PETSc Vec, then gathers it to rank 0.
 *
 * This function uses \c ReadFieldData to fill a PETSc Vec with coordinate data,
 * then leverages \c VecToArrayOnRank0 to gather that data into a contiguous array
 * (valid on rank 0 only).
 *
 * @param[in]  timeIndex    The time index used to construct file names.
 * @param[in]  user         Pointer to the user context.
 * @param[out] coordsArray  On rank 0, will point to a newly allocated array holding the coordinates.
 * @param[out] Ncoords      On rank 0, the length of \p coordsArray. On other ranks, 0.
 *
 * @return PetscErrorCode  Returns 0 on success, or non-zero on failures.
 */
PetscErrorCode ReadPositionsFromFile(PetscInt timeIndex,
                                      UserCtx *user,
                                      double **coordsArray,
				     PetscInt *Ncoords);


/**
 * @brief Reads a named field from file into a PETSc Vec, then gathers it to rank 0.
 *
 * This function wraps \c ReadFieldData and \c VecToArrayOnRank0 into a single step.
 * The gathered data is stored in \p scalarArray on rank 0, with its length in \p Nscalars.
 *
 * @param[in]  timeIndex     The time index used to construct file names.
 * @param[in]  fieldName     Name of the field to be read (e.g., "velocity").
 * @param[in]  user          Pointer to the user context.
 * @param[out] scalarArray   On rank 0, a newly allocated array holding the field data.
 * @param[out] Nscalars      On rank 0, length of \p scalarArray. On other ranks, 0.
 *
 * @return PetscErrorCode  Returns 0 on success, or non-zero on failures.
 */
PetscErrorCode ReadFieldDataToRank0(PetscInt timeIndex,
                                           const char *fieldName,
                                           UserCtx *user,
                                           double **scalarArray,
				    PetscInt *Nscalars);


#endif // IO_H
