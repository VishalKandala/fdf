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
    LOG(LOCAL, LOG_DEBUG, "SetLocalCartesianVelocity: Input Coordinates - x: %f, y: %f, z: %f \n", coor->x, coor->y, coor->z);

    // Compute velocity components as the sine of the coordinates
    ucat->x = sin(coor->x);
    ucat->y = sin(coor->y);
    ucat->z = sin(coor->z);

    // Log computed velocity components
    LOG(LOCAL, LOG_DEBUG, "SetLocalCartesianVelocity: Computed Velocity - x: %f, y: %f, z: %f \n", ucat->x, ucat->y, ucat->z);

    return 0;
}

/**
 * @brief Interpolates a vector field from cell corners to cell centers using averaging.
 *
 * This function computes the field value at the center of each cell by averaging the field 
 * values at the eight surrounding corners of the cell. This is often used to obtain smoother
 * data for computations requiring cell-centered quantities.
 *
 * @param[in]  field       A 3D array of `Cmpnts` structures representing the field vectors at cell corners.
 * @param[out] centfield   A 3D array of `Cmpnts` structures where interpolated field vectors at cell centers will be stored.
 * @param[in]  info        Pointer to a `DMDALocalInfo` structure containing grid information.
 *
 * @return PetscErrorCode Returns 0 on successful execution, non-zero on failure.
 */
PetscErrorCode InterpolateFieldFromCornerToCenter(Cmpnts ***field, Cmpnts ***centfield, DMDALocalInfo *info) {
    PetscInt i, j, k, ic, jc, kc;
    PetscInt lxs,lys,lzs,lxe,lye,lze;

    PetscInt mx = info->mx, my = info->my, mz = info->mz;

    PetscInt xs = info->xs, xe = info->xs + info->xm;
    PetscInt ys = info->ys, ye = info->ys + info->ym;
    PetscInt zs = info->zs, ze = info->zs + info->zm;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
  
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
  
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

    LOG(LOCAL, LOG_INFO, "InterpolateFieldFromCornerToCenter: Starting interpolation with local ranges xs=%d, xe=%d, ys=%d, ye=%d, zs=%d, ze=%d. \n", xs, xe, ys, ye, zs, ze);

    // Loop over all cell centers in the local grid
    for (k = lzs; k < lze; k++) {
        for (j = lys; j < lye; j++) {
            for (i = lxs; i < lxe; i++) {
                // Compute local center indices
                ic = i - lxs;
                jc = j - lys;
                kc = k - lzs;

                // Average the surrounding corner values to compute the center value
                centfield[kc][jc][ic].x = (field[k][j][i].x + field[k][j][i+1].x +
                                           field[k][j+1][i].x + field[k][j+1][i+1].x +
                                           field[k+1][j][i].x + field[k+1][j][i+1].x +
                                           field[k+1][j+1][i].x + field[k+1][j+1][i+1].x) / 8.0;

                centfield[kc][jc][ic].y = (field[k][j][i].y + field[k][j][i+1].y +
                                           field[k][j+1][i].y + field[k][j+1][i+1].y +
                                           field[k+1][j][i].y + field[k+1][j][i+1].y +
                                           field[k+1][j+1][i].y + field[k+1][j+1][i+1].y) / 8.0;

                centfield[kc][jc][ic].z = (field[k][j][i].z + field[k][j][i+1].z +
                                           field[k][j+1][i].z + field[k][j+1][i+1].z +
                                           field[k+1][j][i].z + field[k+1][j][i+1].z +
                                           field[k+1][j+1][i].z + field[k+1][j+1][i+1].z) / 8.0;

                LOG(LOCAL, LOG_DEBUG, "InterpolateFieldFromCornerToCenter: Center (%d, %d, %d) - x: %f, y: %f, z: %f \n",
                    ic, jc, kc, centfield[kc][jc][ic].x, centfield[kc][jc][ic].y, centfield[kc][jc][ic].z);
            }
        }
    }

    LOG(LOCAL, LOG_INFO, "InterpolateFieldFromCornerToCenter: Completed interpolation.\n");
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

    LOG(LOCAL, LOG_INFO, "Allocate3DArray: Allocating 3D array of size (%d x %d x %d).\n", nz, ny, nx);

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

    LOG(LOCAL, LOG_INFO, "Allocate3DArray: Successfully allocated 3D array.\n");
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

    LOG(LOCAL, LOG_INFO, "Deallocate3DArray: Deallocating 3D array of size (%d x %d).\n", nz, ny);

    // Free memory layer by layer
    for (PetscInt k = 0; k < nz; k++) {
        for (PetscInt j = 0; j < ny; j++) {
            ierr = PetscFree(array[k][j]); CHKERRQ(ierr);
        }
        ierr = PetscFree(array[k]); CHKERRQ(ierr);
    }
    ierr = PetscFree(array); CHKERRQ(ierr);

    LOG(LOCAL, LOG_INFO, "Deallocate3DArray: Successfully deallocated 3D array.\n");
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
PetscErrorCode UpdateCartesianVelocity(UserCtx *user) {
    PetscErrorCode ierr;

    // Declare required variables
    Vec Coor;
    Cmpnts ***ucat = NULL, ***coor = NULL, ***centcoor = NULL;
    PetscInt i, j, k;
    PetscInt lxs,lys,lzs,lxe,lye,lze;

    LOG(GLOBAL, LOG_INFO, "UpdateCartesianVelocity: Starting velocity update process.\n");

    // Retrieve local DMDA grid information
    DMDALocalInfo info = user->info;

    PetscInt mx = info.mx, my = info.my, mz = info.mz;

    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
  
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
  
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

    LOG(GLOBAL, LOG_INFO, "UpdateCartesianVelocity: Local subdomain ranges - lxs: %d, lxe: %d, lys: %d, lye: %d, lzs: %d, lze: %d.\n", lxs, lxe, lys, lye, lzs,lze);

    // Access the DMDA coordinate vector
    ierr = DMGetCoordinatesLocal(user->da, &Coor); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, Coor, &coor); CHKERRQ(ierr);

    // Allocate memory for the temporary 3D array to store interpolated coordinates
    ierr = Allocate3DArray(&centcoor, lze-lzs, lye-lys,lxe-lxs); CHKERRQ(ierr);

    LOG(GLOBAL, LOG_DEBUG, "UpdateCartesianVelocity: Allocated centcoor for interpolated coordinates.\n");

    // Access the Cartesian velocity vector
    ierr = DMDAVecGetArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);

    // Interpolate coordinate values from cell corners to cell centers
    LOG(GLOBAL, LOG_INFO, "UpdateCartesianVelocity: Interpolating coordinates from corners to centers.\n");
    ierr = InterpolateFieldFromCornerToCenter(coor, centcoor, &info); CHKERRQ(ierr);

    // Update the Cartesian velocity at each cell center
    LOG(GLOBAL, LOG_INFO, "UpdateCartesianVelocity: Updating velocity values at cell centers.\n");
    for (k = lzs; k < lze; k++) {
        for (j = lys; j < lye; j++) {
            for (i = lxs; i < lxe; i++) {
                ierr = SetLocalCartesianVelocity(&ucat[k][j][i], &centcoor[k-lzs][j-lys][i-lxs]); CHKERRQ(ierr);
                LOG(GLOBAL, LOG_DEBUG, "Updated velocity at (%d, %d, %d): x=%f, y=%f, z=%f \n",
                    k, j, i, ucat[k][j][i].x, ucat[k][j][i].y, ucat[k][j][i].z);
            }
        }
    }

    // Deallocate memory for the temporary interpolated array
    ierr = Deallocate3DArray(centcoor,lze-lzs,lye-lys); CHKERRQ(ierr);
    LOG(GLOBAL, LOG_DEBUG, "UpdateCartesianVelocity: Deallocated centcoor.\n");

    // Restore coordinate and velocity arrays
    ierr = DMDAVecRestoreArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, Coor, &coor); CHKERRQ(ierr);

    LOG(GLOBAL, LOG_INFO, "UpdateCartesianVelocity: Completed velocity update process.\n");
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
 * This version introduces a command-line option `-read_fields` to toggle between:
 * - Updating and writing the Cartesian velocity fields (default)
 * - Reading the Cartesian velocity fields from file
 *
 * Usage example:
 * - To run the code with update and write steps:
 *   `mpirun -np X ./your_executable`
 * - To run the code with reading fields from file:
 *   `mpirun -np X ./your_executable -read_fields`
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 *
 * @return int Returns 0 on success, non-zero on failure.
 */
int main(int argc, char **argv) {
    UserCtx *user = NULL;     // User-defined context holding simulation parameters and data structures
    PetscErrorCode ierr;      // Error handling variable for PETSc routines
    PetscInt block_number = 1; // Number of grid blocks (default to 1, may be updated by initialization)
    PetscInt rstart = 0;      // Starting index for particles or subdomain ranges
    PetscInt np = 0;          // Number of particles
    PetscInt ti = 0;          // Time index or iteration counter
    PetscInt rank, size;      // MPI rank and size variables
    BoundingBox *bboxlist;    // Array of bounding boxes collected from all ranks
    PetscReal umax;           // Variable to store the maximum velocity magnitude
    PetscBool readFields = PETSC_FALSE; // Flag indicating whether to read fields from file
    static char help[] = " Test for interpolation - swarm-curvIB";

    // Initialize PETSc
    ierr = PetscInitialize(&argc, &argv, (char *)0, help); CHKERRQ(ierr);

    // Check if the user requested reading fields from files instead of updating them
    ierr = PetscOptionsGetBool(NULL, NULL, "-read_fields", &readFields, NULL); CHKERRQ(ierr);

    // Initialize simulation data structures, MPI, and other parameters
    ierr = InitializeSimulation(&user, &rank, &size, &np, &rstart, &ti, &block_number); CHKERRQ(ierr);

    // Setup the computational grid, including DMDAs and solution vectors
    ierr = SetupGridAndVectors(user, block_number); CHKERRQ(ierr);

    // If we are NOT reading the fields, we update and write them out
    if (!readFields) {
        // Update the Cartesian Velocity Vector based on current conditions
        ierr = UpdateCartesianVelocity(user); CHKERRQ(ierr);

        // Write the current Cartesian Velocity Vector to file for future reference
        ierr = WriteSimulationFields(user); CHKERRQ(ierr);
    } else {
        // If readFields == PETSC_TRUE, we skip the update and write steps, and read from file instead
        // Read the Simulation Fields (velocity, etc.) from the file for the given time index ti
        ierr = ReadSimulationFields(user, ti); CHKERRQ(ierr);
    }

    // Compute and display the maximum velocity magnitude in the global Ucat vector
    ierr = VecNorm(user->Ucat, NORM_INFINITY, &umax); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "Maximum velocity magnitude: %f\n", umax);

    // Gather the bounding boxes of subdomains from all ranks
    ierr = GatherAllBoundingBoxes(user, &bboxlist); CHKERRQ(ierr);

    // Perform operations related to the particle swarm, such as initialization, movement, and interpolation
    ierr = PerformParticleSwarmOperations(user, np); CHKERRQ(ierr);

    // Finalize the simulation by cleaning up data structures, closing files, and freeing memory
    ierr = FinalizeSimulation(user, block_number); CHKERRQ(ierr);

    return 0;
}
