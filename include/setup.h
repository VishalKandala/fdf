#ifndef SETUP_H
#define SETUP_H

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscsys.h>
#include <petscdmcomposite.h>
#include <petscsystypes.h>

// Include additional headers
#include "common.h"         // Shared type definitions
#include "ParticleSwarm.h"  // Particle swarm functions
#include "walkingsearch.h"  // Particle location functions
#include "grid.h"           // Grid functions
#include "logging.h"        // Logging macros
#include "io.h"             // Data Input and Output functions
#include "interpolation.h"  // Interpolation routines
#include "AnalyticalSolution.h" // Analytical Solution for testing
#include "ParticleMotion.h" // Functions related to motion of particles
#include "Boundaries.h"     //  Functions related to Boundary conditions



/* Macro to automatically select the correct allocation function */
#define Allocate3DArray(array, nz, ny, nx)  \
  _Generic((array),                          \
    PetscReal ****: Allocate3DArrayScalar,       \
    Cmpnts    ****: Allocate3DArrayVector            \
  )(array, nz, ny, nx)

/* Macro for deallocation */
#define Deallocate3DArray(array, nz, ny)  \
  _Generic((array),                        \
    PetscReal ***: Deallocate3DArrayScalar,     \
    Cmpnts    ***: Deallocate3DArrayVector          \
  )(array, nz, ny)


 /**
    * @brief Register events for logging.
    * 
    * Registers events for logging purposes using PetscLogEventRegister.
    *  
 */
PetscErrorCode registerEvents(void);

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
 * @param[out] OnlySetup Flag indicating if only setup should be run (for debugging) without advancing the simulation in time.
 *
 * @return PetscErrorCode Returns 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode InitializeSimulation(UserCtx **user, PetscMPIInt *rank, PetscMPIInt *size, PetscInt *np, PetscInt *StartStep, PetscInt *StepsToRun,PetscReal *StartTime, PetscInt *nblk, PetscInt *outputFreq, PetscBool *readFields, char ***allowedFuncs, PetscInt *nAllowed, char *allowedFile, PetscBool *useCfg, PetscBool *OnlySetup);

/** 
 * @brief Setup grid and vectors for the simulation.
 *
 * @param[in,out] user Pointer to the UserCtx structure.
 * @param[in] block_number The number of grid blocks in the domain.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode SetupGridAndVectors(UserCtx *user, PetscInt block_number);

/**
 * @brief Finalize the simulation and free resources.
 *
 * @param[in,out] user Pointer to the UserCtx structure.
 * @param[in] block_number The number of grid blocks in the domain.
 * @param[in,out] bboxlist  Pointer to the array of BoundingBoxes.
 * @param[in,out] allowedFuncs list of functions that are allowed to show output.
 * @param[in] nAllowed No.of functions allowed to show output.
 * @param[in] PetSc viewer data type for logging 
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode FinalizeSimulation(UserCtx *user, PetscInt block_number, BoundingBox *bboxlist, char **allowedFuncs, PetscInt nAllowed,PetscViewer *logviewer);

/**
 * @brief Allocates a 3D array of PetscReal values using PetscCalloc.
 *
 * This function dynamically allocates memory for a 3D array of PetscReal values
 * with dimensions nz (layers) x ny (rows) x nx (columns). It uses PetscCalloc1
 * to ensure the memory is zero-initialized.
 *
 * The allocation is done in three steps:
 *  1. Allocate an array of nz pointers (one for each layer).
 *  2. Allocate a contiguous block for nz*ny row pointers and assign each layer’s row pointers.
 *  3. Allocate a contiguous block for all nz*ny*nx PetscReal values.
 *
 * This setup allows the array to be accessed as array[k][j][i], and the memory
 * for the data is contiguous, which improves cache efficiency.
 *
 * @param[out] array Pointer to the 3D array to be allocated.
 * @param[in]  nz    Number of layers (z-direction).
 * @param[in]  ny    Number of rows (y-direction).
 * @param[in]  nx    Number of columns (x-direction).
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
 PetscErrorCode Allocate3DArrayScalar(PetscReal ****array, PetscInt nz, PetscInt ny, PetscInt nx);

/**
 * @brief Deallocates a 3D array of PetscReal values allocated by Allocate3DArrayScalar.
 *
 * This function frees the memory allocated for a 3D array of PetscReal values.
 * It assumes the memory was allocated using Allocate3DArrayScalar, which allocated
 * three separate memory blocks: one for the contiguous data, one for the row pointers,
 * and one for the layer pointers.
 *
 * @param[in] array Pointer to the 3D array to be deallocated.
 * @param[in] nz    Number of layers (z-direction).
 * @param[in] ny    Number of rows (y-direction).
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
PetscErrorCode Deallocate3DArrayScalar(PetscReal ***array, PetscInt nz, PetscInt ny);

/**
 * @brief Deallocates a 3D array of Cmpnts structures allocated by Allocate3DArrayVector.
 *
 * This function frees the memory allocated for a 3D array of Cmpnts structures.
 * It assumes the memory was allocated using Allocate3DArrayVector, which created three
 * separate memory blocks: one for the contiguous vector data, one for the row pointers,
 * and one for the layer pointers.
 *
 * @param[in] array Pointer to the 3D array to be deallocated.
 * @param[in] nz    Number of layers in the z-direction.
 * @param[in] ny    Number of rows in the y-direction.
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
 PetscErrorCode Allocate3DArrayVector(Cmpnts ****array, PetscInt nz, PetscInt ny, PetscInt nx);

/**
 * @brief Deallocates a 3D array of Cmpnts structures allocated by Allocate3DArrayVector.
 *
 * This function frees the memory allocated for a 3D array of Cmpnts structures.
 * It assumes the memory was allocated using Allocate3DArrayVector, which created three
 * separate memory blocks: one for the contiguous vector data, one for the row pointers,
 * and one for the layer pointers.
 *
 * @param[in] array Pointer to the 3D array to be deallocated.
 * @param[in] nz    Number of layers in the z-direction.
 * @param[in] ny    Number of rows in the y-direction.
 *
 * @return PetscErrorCode 0 on success, nonzero on failure.
 */
PetscErrorCode Deallocate3DArrayVector(Cmpnts ***array, PetscInt nz, PetscInt ny);

/**
 * @brief Determines the global starting index and number of CELLS owned by the
 *        current processor in a specified dimension. Ownership is defined by the
 *        rank owning the cell's origin node (min i,j,k corner).
 *
 * @param[in]  info_nodes     Pointer to the DMDALocalInfo struct obtained from the NODE-based DMDA
 *                            (e.g., user->da or user->fda, assuming they have consistent nodal partitioning
 *                            for defining cell origins).
 * @param[in]  dim            The dimension to compute the range for (0 for x/i, 1 for y/j, 2 for z/k).
 * @param[out] xs_cell_global Pointer to store the starting GLOBAL CELL index owned by this process.
 *                            A cell C(i) is defined by nodes N(i) and N(i+1). Its global index is i.
 * @param[out] xm_cell_local  Pointer to store the NUMBER of CELLs owned by this process in this dimension.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode GetOwnedCellRange(const DMDALocalInfo *info_nodes,
                                 PetscInt dim,
                                 PetscInt *xs_cell_global,
                                 PetscInt *xm_cell_local);

/**
 * @brief Gets the global cell range for a rank, including boundary cells.
 *
 * This function first calls GetOwnedCellRange to get the conservative range of
 * fully-contained cells. It then extends this range by applying the
 * "Lower-Rank-Owns-Boundary" principle. A rank claims ownership of the
 * boundary cells it shares with neighbors in the positive (+x, +y, +z)
 * directions.
 *
 * This results in a final cell range that is gap-free and suitable for building
 * the definitive particle ownership map.
 *
 * @param[in]  info_nodes       Pointer to the DMDALocalInfo struct.
 * @param[in]  neighbors        Pointer to the RankNeighbors struct containing neighbor info.
 * @param[in]  dim              The dimension (0 for i/x, 1 for j/y, 2 for k/z).
 * @param[out] xs_cell_global_out Pointer to store the final starting cell index.
 * @param[out] xm_cell_local_out  Pointer to store the final number of cells.
 *
 * @return PetscErrorCode 0 on success, or an error code on failure.
 */
PetscErrorCode GetGhostedCellRange(const DMDALocalInfo *info_nodes,
                                   const RankNeighbors *neighbors,
                                   PetscInt dim,
                                   PetscInt *xs_cell_global_out,
                                   PetscInt *xm_cell_local_out);

/**
 * @brief Updates the local vector (including ghost points) from its corresponding global vector.
 *
 * This function identifies the correct global vector, local vector, and DM based on the
 * provided fieldName and performs the standard PETSc DMGlobalToLocalBegin/End sequence.
 * Includes optional debugging output (max norms before/after).
 *
 * @param user       The UserCtx structure containing the vectors and DMs.
 * @param fieldName  The name of the field to update ("Ucat", "Ucont", "P", "Nvert", etc.).
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 *
 * @note This function assumes the global vector associated with fieldName has already
 *       been populated with the desired data (including any boundary conditions).
 */
PetscErrorCode UpdateLocalGhosts(UserCtx* user, const char *fieldName);

/**
 * @brief Computes and stores the Cartesian neighbor ranks for the DMDA decomposition.
 *
 * This function retrieves the neighbor information from the primary DMDA (user->da)
 * and stores the face neighbors (xm, xp, ym, yp, zm, zp) in the user->neighbors structure.
 * It assumes a standard PETSc ordering for the neighbors array returned by DMDAGetNeighbors.
 * Logs warnings if the assumed indices seem incorrect (e.g., center rank mismatch).
 *
 * @param[in,out] user Pointer to the UserCtx structure where neighbor info will be stored.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode ComputeAndStoreNeighborRanks(UserCtx *user);

/**
 * @brief Sets the processor layout for a given DMDA based on PETSc options.
 *
 * Reads the desired number of processors in x, y, and z directions using
 * PETSc options (e.g., -dm_processors_x, -dm_processors_y, -dm_processors_z).
 * If an option is not provided for a direction, PETSC_DECIDE is used for that direction.
 * Applies the layout using DMDASetNumProcs.
 *
 * Also stores the retrieved/decided values in user->procs_x/y/z if user context is provided.
 *
 * @param dm   The DMDA object to configure the layout for.
 * @param user Pointer to the UserCtx structure (optional, used to store layout values).
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode SetDMDAProcLayout(DM dm, UserCtx *user);

/**
 * @brief Sets up the full rank communication infrastructure, including neighbor ranks and bounding box exchange.
 *
 * This function orchestrates the following steps:
 * 1. Compute and store the neighbor ranks in the user context.
 * 2. Gather all local bounding boxes to rank 0.
 * 3. Broadcast the complete bounding box list to all ranks.
 *
 * The final result is that each rank has access to its immediate neighbors and the bounding box information of all ranks.
 *
 * @param[in,out] user      Pointer to the UserCtx structure (must be initialized).
 * @param[in,out] bboxlist  Pointer to BoundingBox array pointer; after this call, it will point to the broadcasted list.
 *
 * @return PetscErrorCode Returns 0 on success or non-zero PETSc error code.
 */
PetscErrorCode SetupDomainRankInfo(UserCtx *user, BoundingBox **bboxlist);

/**
 * @brief Reconstructs Cartesian velocity (Ucat) at cell centers from contravariant
 *        velocity (Ucont) defined on cell faces.
 *
 * This function performs the transformation from a contravariant velocity representation
 * (which is natural on a curvilinear grid) to a Cartesian (x,y,z) representation.
 * For each interior computational cell owned by the rank, it performs the following:
 *
 * 1.  It averages the contravariant velocity components (U¹, U², U³) from the
 *     surrounding faces to get an estimate of the contravariant velocity at the cell center.
 * 2.  It averages the metric vectors (Csi, Eta, Zet) from the surrounding faces
 *     to get an estimate of the metric tensor at the cell center. This tensor forms
 *     the transformation matrix.
 * 3.  It solves the linear system `[MetricTensor] * [ucat] = [ucont]` for the
 *     Cartesian velocity vector `ucat = (u,v,w)` using Cramer's rule.
 * 4.  The computed Cartesian velocity is stored in the global `user->Ucat` vector.
 *
 * The function operates on local, ghosted versions of the input vectors (`user->lUcont`,
 * `user->lCsi`, etc.) to ensure stencils are valid across processor boundaries.
 *
 * @param[in,out] user      Pointer to the UserCtx structure. The function reads from
 *                          `user->lUcont`, `user->lCsi`, `user->lEta`, `user->lZet`, `user->lNvert`
 *                          and writes to the global `user->Ucat` vector.
 *
 * @return PetscErrorCode 0 on success.
 *
 * @note
 *  - This function should be called AFTER `user->lUcont` and all local metric vectors
 *    (`user->lCsi`, etc.) have been populated with up-to-date ghost values via `UpdateLocalGhosts`.
 *  - It only computes `Ucat` for interior cells (not on physical boundaries) and for
 *    cells not marked as solid/blanked by `user->lNvert`.
 *  - The caller is responsible for subsequently applying boundary conditions to `user->Ucat`
 *    and calling `UpdateLocalGhosts(user, "Ucat")` to populate `user->lUcat`.
 */
PetscErrorCode Contra2Cart(UserCtx *user);

/**
 * @brief Sets up the entire boundary condition system for the simulation.
 *
 * This function is the main entry point for all boundary condition setup. It performs
 * two main tasks:
 *   1. It determines the name of the boundary condition file, using "bcs.dat" as a
 *      default but allowing it to be overridden by the PETSc option `-bcs_file`.
 *   2. It then calls the core `BoundarySystem_Create` function, which reads the file,
 *      creates all the necessary handler objects, and calls their Initialize() methods
 *      to set the initial state of the boundary fields.
 *
 * This function should be called in `main()` AFTER `SetupGridAndVectors()` has completed
 * to ensure that the grid DMs and DMDALocalInfo are valid.
 *
 * @param user The main UserCtx struct.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetupBoundaryConditions(UserCtx *user);

/**
 * @brief Creates and distributes a map of the domain's cell decomposition to all ranks.
 * @ingroup DomainInfo
 *
 * This function is a critical part of the simulation setup. It determines the global
 * cell ownership for each MPI rank and makes this information available to all
 * other ranks. This "decomposition map" is essential for the robust "Walk and Handoff"
 * particle migration strategy, allowing any rank to quickly identify the owner of a
 * target cell.
 *
 * The process involves:
 * 1. Each rank gets its own node ownership information from the DMDA.
 * 2. It converts this node information into cell ownership ranges using the
 *    `GetOwnedCellRange` helper function.
 * 3. It participates in an `MPI_Allgather` collective operation to build a complete
 *    array (`user->RankCellInfoMap`) containing the ownership information for every rank.
 *
 * This function should be called once during initialization after the primary DMDA
 * (user->da) has been set up.
 *
 * @param[in,out] user Pointer to the UserCtx structure. The function will allocate and
 *                     populate `user->RankCellInfoMap` and set `user->num_ranks`.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 *         Errors can occur if input pointers are NULL or if MPI communication fails.
 */
PetscErrorCode SetupDomainCellDecompositionMap(UserCtx *user);

/**
 * @brief Performs a binary search for a key in a sorted array of PetscInt64.
 *
 * This is a standard binary search algorithm implemented as a PETSc-style helper function.
 * It efficiently determines if a given `key` exists within a `sorted` array.
 *
 * @param[in]  n      The number of elements in the array.
 * @param[in]  arr    A pointer to the sorted array of PetscInt64 values to be searched.
 * @param[in]  key    The PetscInt64 value to search for.
 * @param[out] found  A pointer to a PetscBool that will be set to PETSC_TRUE if the key
 *                    is found, and PETSC_FALSE otherwise.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 *
 * @note The input array `arr` **must** be sorted in ascending order for the algorithm
 *       to work correctly.
 */
PetscErrorCode BinarySearchInt64(PetscInt n, const PetscInt64 arr[], PetscInt64 key, PetscBool *found);

 #endif // SETUP_H
