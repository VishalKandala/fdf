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
 * @brief Initializes or updates all necessary simulation fields for a given timestep.
 *
 * This function handles the logic for either:
 * A) Initializing fields analytically (for the first step or if not reading):
 *    - Sets interior values using SetAnalyticalCartesianField.
 *    - Applies boundary conditions using ApplyAnalyticalBC.
 *    - Updates local vectors with ghosts using UpdateLocalGhosts.
 *    - Optionally writes the initial fields.
 * B) Reading fields from a file for a specific timestep index:
 *    - Reads global vectors using ReadSimulationFields.
 *    - Updates local vectors with ghosts using UpdateLocalGhosts.
 * C) Updating fields using a fluid solver (Placeholder for future integration):
 *    - Calls a placeholder function SolveFluidEquations.
 *    - Applies boundary conditions using ApplyAnalyticalBC.
 *    - Updates local vectors with ghosts using UpdateLocalGhosts.
 *
 * @param user        Pointer to the UserCtx structure.
 * @param step        The current timestep number (0 for initial step).
 * @param time        The current simulation time.
 * @param readFields  Flag indicating whether to read fields from file.
 * @param fieldSource Source for field data (e.g., ANALYTICAL, FILE, SOLVER).
 *                    (Here using readFields bool for simplicity based on original code)
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode SetEulerianFields(UserCtx *user, PetscInt step, PetscInt StartStep, PetscReal time, PetscBool readFields);


/**
 * @brief Performs the complete initial setup for the particle simulation at time t=0.
 *
 * This includes:
 * 1. Initial locating of particles (based on their potentially arbitrary initial assignment).
 * 2. A preliminary migration cycle to ensure particles are on the MPI rank that owns
 *    their initial physical region.
 * 3. If `user->ParticleInitialization == 0` (Surface Init), re-initializes particles on the
 *    designated inlet surface. This ensures particles migrated to an inlet-owning rank
 *    are correctly distributed on that surface.
 * 4. A final locating of all particles to get their correct cell indices and interpolation weights.
 * 5. Interpolation of initial Eulerian fields to the particles.
 * 6. Scattering of particle data to Eulerian fields (if applicable).
 * 7. Outputting initial data if requested.
 *
 * @param user Pointer to the UserCtx structure.
 * @param currentTime The current simulation time (should be StartTime, typically 0.0).
 * @param step The current simulation step (should be StartStep, typically 0).
 * @param readFields Flag indicating if Eulerian fields were read from file (influences output).
 * @param bboxlist Array of BoundingBox structures for domain decomposition.
 * @param OutputFreq Frequency for writing output files.
 * @param StepsToRun Total number of simulation steps planned (used for output logic on setup-only runs).
 * @param StartStep The starting step of the simulation (used for output logic).
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformInitialSetup(UserCtx *user, PetscReal currentTime, PetscInt step,
                                   PetscBool readFields, const BoundingBox *bboxlist,
                                   PetscInt OutputFreq, PetscInt StepsToRun, PetscInt StartStep);

/**
 * @brief Executes the main time-marching loop for the particle simulation.
 *
 * This function performs the following steps repeatedly:
 * 1. Updates/Sets the background fluid velocity field (Ucat) for the current step.
 * 2. Updates particle positions using velocity from the *previous* step's interpolation.
 *    (Note: For the very first step (step=StartStep), the velocity used might be zero
 *     or an initial guess if not handled carefully).
 * 3. Locates particles in the grid based on their *new* positions.
 * 4. Interpolates the fluid velocity (from the *current* Ucat) to the new particle locations.
 * 5. Logs errors and outputs data at specified intervals.
 *
 * @param user         Pointer to the UserCtx structure.
 * @param StartStep    The initial step number (e.g., 0 for a new run, >0 for restart).
 * @param StartTime    The simulation time corresponding to StartStep.
 * @param StepsToRun   The number of steps to execute in this run.
 * @param OutputFreq   Frequency (in number of steps) at which to output data and log errors.
 * @param readFields   Flag indicating whether to read initial fields (only at StartStep).
 * @param bboxlist     A list that contains the bounding boxes of all the ranks.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode AdvanceSimulation(UserCtx *user, PetscInt StartStep, PetscReal StartTime, PetscInt StepsToRun, PetscInt OutputFreq, PetscBool readFields, const BoundingBox *bboxlist);

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

 #endif // SETUP_H
