/**
 * @file inttest.c  //  Particle In Cell main.
 * @brief Test program for DMSwarm interpolation using the fdf-curvIB method.
 *
 * Initializes a particle swarm, reads velocity fields, and performs particle-grid interpolation.
 */

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscdmcomposite.h>

// Include the updated headers
//#include "common.h"         // Shared type definitions
//#include "ParticleSwarm.h"  // Particle swarm functions
//#include "walkingsearch.h"  // Particle location functions
//#include "grid.h"           // Grid functions
//#include "logging.h"        // Logging macros
//#include "io.h"             // Data Input and Output functions

#include "interpolation.h"

PetscErrorCode registerEvents(void) {
    PetscErrorCode ierr;
    // Register the event with a descriptive name
    ierr = PetscLogEventRegister("walkingsearch", PETSC_OBJECT_CLASSID, &EVENT_walkingsearch);
    ierr = PetscLogEventRegister("Individualwalkingsearch", PETSC_OBJECT_CLASSID, &EVENT_Individualwalkingsearch);
    CHKERRQ(ierr);
    return 0;
}


/**
 * @brief Sets the local Cartesian velocity components based on coordinates.
 *
 * This function computes the velocity components by applying the sine function to the 
 * input coordinate values along the x, y, and z directions.
 *
 * @param[in,out] ucont Pointer to the `Cmpnts` structure where the velocity components will be stored.
 * @param[in]     coor  Pointer to the `Cmpnts` structure containing the input coordinate values.
 *
 * @return PetscErrorCode Returns 0 on successful execution, non-zero on failure.
 */
static inline PetscErrorCode SetLocalCartesianVelocity(Cmpnts *ucat, Cmpnts *coor) {
    // Log input coordinate values for debugging purposes
    LOG_ALLOW(LOCAL, LOG_DEBUG, "SetLocalCartesianVelocity: Input Coordinates - x: %f, y: %f, z: %f \n", coor->x, coor->y, coor->z);

    // Compute velocity components as the sine of the coordinates
    ucat->x = sin(coor->x);
    ucat->y = sin(coor->y);
    ucat->z = sin(coor->z);

    // Log computed velocity components
    LOG_ALLOW(LOCAL, LOG_DEBUG, "SetLocalCartesianVelocity: Computed Velocity - x: %f, y: %f, z: %f \n", ucat->x, ucat->y, ucat->z);

    return 0;
}

/**
 * @brief Interpolates a vector field from cell corners to cell centers using simple averaging.
 *
 * This version loops only up to `xe-1, ye-1, ze-1`, ensuring `i+1, j+1, k+1` are in range.
 * The result is stored in `centfield[k][j][i]` at the same (k,j,i) offsets
 * (assuming the caller knows to shift or match indexing carefully).
 *
 * @param[in]  field     A 3D array of `Cmpnts` at cell corners (size at least [info->ze][info->ye][info->xe]).
 * @param[out] centfield A 3D array of `Cmpnts` into which cell-center data is written.
 *                       Must be allocated by the caller with same local extents or an offset approach.
 * @param[in]  info      DMDALocalInfo with local domain indices [xs..xe), etc.
 *
 * @return PetscErrorCode
 */

PetscErrorCode InterpolateFieldFromCornerToCenter(Cmpnts ***field,
                                                  Cmpnts ***centfield,
                                                  DMDALocalInfo *info)
{
    PetscInt i, j, k;
    PetscInt xs = info->xs, xe = info->xs + info->xm;
    PetscInt ys = info->ys, ye = info->ys + info->ym;
    PetscInt zs = info->zs, ze = info->zs + info->zm;

    /* we must stop at xe-1 so i+1 is valid local indexing */
    PetscInt iend = xe - 1, jend = ye - 1, kend = ze - 1;

    for (k = zs; k < kend; k++) {
        for (j = ys; j < jend; j++) {
            for (i = xs; i < iend; i++) {
                centfield[k][j][i].x =
                   ( field[k][j][i].x + field[k][j][i+1].x +
                     field[k][j+1][i].x + field[k][j+1][i+1].x +
                     field[k+1][j][i].x + field[k+1][j][i+1].x +
                     field[k+1][j+1][i].x + field[k+1][j+1][i+1].x ) / 8.0;

                centfield[k][j][i].y =
                   ( field[k][j][i].y + field[k][j][i+1].y +
                     field[k][j+1][i].y + field[k][j+1][i+1].y +
                     field[k+1][j][i].y + field[k+1][j][i+1].y +
                     field[k+1][j+1][i].y + field[k+1][j+1][i+1].y ) / 8.0;

                centfield[k][j][i].z =
                   ( field[k][j][i].z + field[k][j][i+1].z +
                     field[k][j+1][i].z + field[k][j+1][i+1].z +
                     field[k+1][j][i].z + field[k+1][j][i+1].z +
                     field[k+1][j+1][i].z + field[k+1][j+1][i+1].z ) / 8.0;
            }
        }
    }
    return 0;
}

/**
 * @brief Allocates a 3D array of `Cmpnts` structures.
 *
 * This function dynamically allocates memory for a 3D array of `Cmpnts` structures.
 * Each component of the array (x, y, z) is initialized to 0.
 *
 * @param[out] array Pointer to the 3D array to be allocated.
 * @param[in]  nz    Number of layers in the z-direction.
 * @param[in]  ny    Number of rows in the y-direction.
 * @param[in]  nx    Number of columns in the x-direction.
 *
 * @return PetscErrorCode Returns 0 on successful allocation, non-zero on failure.
 */
PetscErrorCode Allocate3DArray(Cmpnts ****array, PetscInt nz, PetscInt ny, PetscInt nx) {
    PetscErrorCode ierr;

    LOG_ALLOW(LOCAL, LOG_INFO, "Allocate3DArray: Allocating 3D array of size (%d x %d x %d).\n", nz, ny, nx);

    // Allocate memory for each dimension
    ierr = PetscMalloc1(nz, array); CHKERRQ(ierr);
    for (PetscInt k = 0; k < nz; k++) {
        ierr = PetscMalloc1(ny, &(*array)[k]); CHKERRQ(ierr);
        for (PetscInt j = 0; j < ny; j++) {
            ierr = PetscMalloc1(nx, &(*array)[k][j]); CHKERRQ(ierr);
            for (PetscInt i = 0; i < nx; i++) {
                // Initialize components to zero
                (*array)[k][j][i].x = 0.0;
                (*array)[k][j][i].y = 0.0;
                (*array)[k][j][i].z = 0.0;
            }
        }
    }

    LOG_ALLOW(LOCAL, LOG_INFO, "Allocate3DArray: Successfully allocated 3D array.\n");
    return 0;
}

/**
 * @brief Deallocates a 3D array of `Cmpnts` structures.
 *
 * This function frees the memory allocated for a 3D array of `Cmpnts` structures.
 *
 * @param[in]  array Pointer to the 3D array to be deallocated.
 * @param[in]  nz    Number of layers in the z-direction.
 * @param[in]  ny    Number of rows in the y-direction.
 *
 * @return PetscErrorCode Returns 0 on successful deallocation, non-zero on failure.
 */
PetscErrorCode Deallocate3DArray(Cmpnts ***array, PetscInt nz, PetscInt ny) {
    PetscErrorCode ierr;

    LOG_ALLOW(LOCAL, LOG_INFO, "Deallocate3DArray: Deallocating 3D array of size (%d x %d).\n", nz, ny);

    // Free memory layer by layer
    for (PetscInt k = 0; k < nz; k++) {
        for (PetscInt j = 0; j < ny; j++) {
            ierr = PetscFree(array[k][j]); CHKERRQ(ierr);
        }
        ierr = PetscFree(array[k]); CHKERRQ(ierr);
    }
    ierr = PetscFree(array); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_INFO, "Deallocate3DArray: Successfully deallocated 3D array.\n");
    return 0;
}

/**
 * @brief Updates the local Cartesian velocity field based on interpolated coordinates.
 *
 * This function performs the following operations:
 * 1. Retrieves the local grid information from the DMDA associated with the user context.
 * 2. Accesses the local coordinate vector and retrieves a read-only array representation.
 * 3. Allocates memory for a temporary 3D array (`centcoor`) to store interpolated coordinate values at cell centers.
 * 4. Accesses the Cartesian velocity vector and obtains a writable array representation.
 * 5. Interpolates the coordinate values from cell corners to cell centers using `InterpolateFieldFromCornerToCenter`.
 * 6. Updates the Cartesian velocity values at each cell center based on the interpolated coordinates using `SetLocalCartesianVelocity`.
 * 7. Deallocates the temporary array (`centcoor`) and restores the arrays to maintain PETSc's internal state.
 *
 * @param[in,out] user Pointer to a `UserCtx` structure containing:
 *                     - `user->da`: DMDA for the grid.
 *                     - `user->fda`: DMDA for the Cartesian velocity field.
 *                     - `user->Ucat`: Cartesian velocity vector.
 *                     - `user->Coor`: Coordinate vector.
 *                     - `user->info`: Local DMDA grid information.
 *
 * @return PetscErrorCode 
 *         - Returns 0 on successful execution.
 *         - Non-zero error code on failure.
 *
 * @note 
 * - Assumes that the coordinate vector (`user->Coor`) has been properly initialized and populated before this function is called.
 * - The function is grid-aware and handles memory allocation and interpolation for local subdomains only.
 */

PetscErrorCode UpdateCartesianVelocity(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo info = user->info;
    Vec Coor;
    Cmpnts ***coor = NULL, ***centcoor = NULL, ***ucat = NULL;
    PetscInt rank;

    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;
 
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "UpdateCartesianVelocity: Starting on rank : %d \n",rank);

    /* 1) Access local coords (Coor) & velocity array (Ucat) */
    ierr = DMGetCoordinatesLocal(user->da, &Coor); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, Coor, &coor); CHKERRQ(ierr);

    ierr = DMDAVecGetArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);

    /* 2) Allocate scratch array 'centcoor' the same local size as [zs..ze, ys..ye, xs..xe]. */
    ierr = Allocate3DArray(&centcoor, ze, ye, xe); CHKERRQ(ierr);

    /* 3) Interpolate corner->center coords into centcoor */
    ierr = InterpolateFieldFromCornerToCenter(coor, centcoor, &info); CHKERRQ(ierr);

    /* 4) Fill interior cell velocity using 'SetLocalCartesianVelocity' with the center coords. 
       (Stop at xe-1, ye-1, ze-1 so i+1 doesn't go out-of-range) */
    for (PetscInt k = zs; k < ze - 1; k++) {
      for (PetscInt j = ys; j < ye - 1; j++) {
        for (PetscInt i = xs; i < xe - 1; i++) {
          ierr = SetLocalCartesianVelocity(&ucat[k][j][i], &centcoor[k][j][i]); CHKERRQ(ierr);
        }
      }
    }

    /* 5) Set boundary nodes to sine-of-corner-coords 
       i.e. if i==xe-1 or i==xs or similarly for j,k. 
       That ensures no boundary node is left zero. */
    for (PetscInt k = zs; k < ze; k++) {
      for (PetscInt j = ys; j < ye; j++) {
        for (PetscInt i = xs; i < xe; i++) {
          PetscBool isBoundary = PETSC_FALSE;
          if ((i == xe - 1) || (i == xs) ||
              (j == ye - 1) || (j == ys)  ||
              (k == ze - 1) || (k == zs)) {
            isBoundary = PETSC_TRUE;
          }
          if (isBoundary) {
            ierr = SetLocalCartesianVelocity(&ucat[k][j][i], &coor[k][j][i]); CHKERRQ(ierr);
          }
        }
      }
    }

    /* 6) Clean up */
    ierr = Deallocate3DArray(centcoor, ze, ye); CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, Coor, &coor); CHKERRQ(ierr);

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "UpdateCartesianVelocity: Completed on rank : %d \n",rank);
    return 0;
}



#undef _FUNCT_
#define __FUNCT__ "main"

/**
 * @brief Main function for DMSwarm interpolation.
 *
 * Initializes the grid, reads/writes the velocity field depending on command line options,
 * creates particles, and performs particle-grid interpolation.
 *
 * Introduces a command-line option `-read_fields` to toggle between:
 * - Updating and writing the Cartesian velocity fields (default)
 * - Reading the Cartesian velocity fields from file
 *
 * Additionally, after gathering the bounding boxes on rank 0, it now calls a separate function,
 * `BroadcastAllBoundingBoxes()`, to distribute bboxlist to all ranks.
 *
 * Usage example:
 * - To run with update and write steps:
 *   `mpirun -np X ./your_executable`
 * - To run with reading fields:
 *   `mpirun -np X ./your_executable -read_fields`
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 *
 * @return int Returns 0 on success, non-zero on failure.
 */
int main(int argc, char **argv) {
    UserCtx *user = NULL;     // User context
    PetscErrorCode ierr;      // PETSc error handling
    PetscInt block_number = 1;
    PetscInt rstart = 0;
    PetscInt np = 0;
    PetscInt ti = 0;
    PetscInt rank, size;
    BoundingBox *bboxlist;    // Array of bounding boxes
    PetscReal umax;
    PetscBool readFields = PETSC_FALSE;
    static char help[] = " Test for interpolation - swarm-curvIB";
    PetscViewer logviewer;

    // -------------------- 1. PETSc Initialization --------------------
    ierr = PetscInitialize(&argc, &argv, (char *)0, help); CHKERRQ(ierr);

    // -------------------- 2. Setup Logging Allow-List ----------------
    // Only these function names will produce LOG_ALLOW (or LOG_ALLOW_SYNC) output.
    // You can add more as needed (e.g., "InitializeSimulation", "PerformParticleSwarmOperations", etc.).
    const char *allowedFuncs[] = {
      //  "main",
      //  "InitializeSimulation",                 
      //  "SetupGridAndVectors",
      //  "UpdateCartesianVelocity",
      // "InterpolateFieldFromCornerToCenter",
      //  "WriteSimulationFields",
      //  "ReadSimulationFields",
      //   "GatherAllBoundingBoxes",
      //   "BroadcastAllBoundingBoxes",
      // "PerformParticleSwarmOperations",
      // "CreateParticleSwarm",
      // "AssignInitialPropertiesToSwarm",
      // "InitializeParticleBasicProperties",
      // "FinalizeSwarmSetup",
      //  "LocateAllParticlesInGrid",
      //"InterpolateParticleVelocities",
           //"ComputeTrilinearWeights",
      // "FinalizeSimulation"
    };
    set_allowed_functions(allowedFuncs, 0);

    // Enable PETSc default logging
    ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);

    registerEvents();   

    print_log_level();

    // Check if user requested to read fields instead of updating
    ierr = PetscOptionsGetBool(NULL, NULL, "-read_fields", &readFields, NULL); CHKERRQ(ierr);

    // Initialize simulation: user context, MPI rank/size, np, etc.
    ierr = InitializeSimulation(&user, &rank, &size, &np, &rstart, &ti, &block_number); CHKERRQ(ierr);

    // Another demonstration of LOG_ALLOW
    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,
              "main: readFields = %s, rank = %d, size = %d\n",
              readFields ? "true" : "false", rank, size);

    // Setup the computational grid
    ierr = SetupGridAndVectors(user, block_number); CHKERRQ(ierr);

    // Either update and write fields, or read fields from file
    if (!readFields) {
        ierr = UpdateCartesianVelocity(user); CHKERRQ(ierr);
        ierr = WriteSimulationFields(user); CHKERRQ(ierr);
    } else {
        ierr = ReadSimulationFields(user, ti); CHKERRQ(ierr);
    }

    // Compute and print maximum velocity magnitude
    ierr = VecNorm(user->Ucat, NORM_INFINITY, &umax); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_INFO,"Maximum velocity magnitude: %f\n", umax);

    // Gather bounding boxes on rank 0
    ierr = GatherAllBoundingBoxes(user, &bboxlist); CHKERRQ(ierr);

    // Broadcast bboxlist to all ranks
    ierr = BroadcastAllBoundingBoxes(user, &bboxlist); CHKERRQ(ierr);

    // Perform particle swarm operations with bboxlist knowledge on all ranks
    ierr = PerformParticleSwarmOperations(user, np, bboxlist); CHKERRQ(ierr);
    
    // Create an ASCII viewer to write log output to file
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "petsc_log.txt", &logviewer);

    // Print PETSc logging results at the end
    ierr = PetscLogView(logviewer); CHKERRQ(ierr);

    // Finalize simulation
    ierr = FinalizeSimulation(user, block_number, bboxlist); CHKERRQ(ierr);

    return 0;
}


