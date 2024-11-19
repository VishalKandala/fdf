/**
 * @file common.h
 * @brief Common type definitions and structures shared across multiple modules.
 *
 * This header file contains the definitions of shared types such as `UserCtx`, `Particle`,
 * `BoundingBox`, `Cmpnts`, and other structures used throughout the simulation.
 * Including this file in other headers and source files ensures consistent type definitions
 * and avoids circular dependencies.
 */

#ifndef COMMON_H
#define COMMON_H

// Include PETSc library header
#include <petsc.h>
#include <petscdmda.h>
#include <petscdmswarm.h>

// Standard C library headers
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

// --------------------- Type Definitions ---------------------

/**
 * @brief Represents a 3D point or vector with x, y, z components.
 */
typedef struct {
    PetscScalar x, y, z;
} Cmpnts;

/**
 * @brief Represents a 2D point or vector with x, y components.
 */
typedef struct {
    PetscScalar x, y;
} Cmpnts2;

/**
 * @brief Represents a flow wave with time and frequency components.
 */
typedef struct {
    PetscReal t, f;
} FlowWave;

/**
 * @brief Represents a bounding box in 3D space defined by minimum and maximum coordinates.
 */
typedef struct {
    Cmpnts min_coords; /**< Minimum x, y, z coordinates of the bounding box. */
    Cmpnts max_coords; /**< Maximum x, y, z coordinates of the bounding box. */
} BoundingBox;

/**
 * @brief Represents a particle with its properties for use in particle simulations.
 */
typedef struct {
    PetscInt64 PID;     /**< Unique Particle ID. */
    PetscInt64 cell[3]; /**< Indices of the cell containing the particle (i, j, k). */
    Cmpnts loc;         /**< Location of the particle in 3D space. */
    Cmpnts vel;         /**< Velocity of the particle in 3D space. */
    Cmpnts weights;     /**< Weights associated with the particle (e.g., for interpolation). */
} Particle;

/**
 * @brief Represents a cell in the computational grid with its eight vertices.
 */
typedef struct {
    Cmpnts vertices[8]; /**< Coordinates of the eight vertices of the cell. */
} Cell;

/**
 * @brief User-defined context containing simulation data and configurations.
 *
 * This structure holds various simulation parameters, PETSc data structures, and
 * other information needed throughout the simulation.
 */
typedef struct UserCtx {
    // PETSc Data Structures
    DM da;              /**< Data structure for scalars (grid geometry information). */
    DM fda;             /**< Data structure for vectors. */
    DM dmcell;          /**< Shell data structure for point location (for PIC methods). */
    DM swarm;           /**< Data structure for particles (DMSwarm). */
    DMDALocalInfo info; /**< Local information about the DMDA grid. */

    // Particle Migration
    PetscInt *miglist;  /**< List of ranks to migrate to during particle migration. */

    // Bounding Box
    BoundingBox bbox;   /**< Local bounding box of the grid on the current process. */

    // Grid and Flow Variables
    Vec Cent;           /**< Coordinates of cell centers. */
    Vec Ucat;           /**< Cartesian velocity components (u, v, w). */
    Vec P;              /**< Pressure field. */
    // Add other Vecs and variables as needed

    // Simulation Parameters
    PetscReal dt;       /**< Time step size. */
    PetscReal ren;      /**< Reynolds number. */

    // Random Number Generators
    PetscRandom randx;  /**< Random number generator for x-coordinate. */
    PetscRandom randy;  /**< Random number generator for y-coordinate. */
    PetscRandom randz;  /**< Random number generator for z-coordinate. */

    // Logging Level
    PetscInt log_level; /**< Logging level for controlling output verbosity. */

    PetscInt	IM, JM, KM; // dimensions of grid
    
    PetscInt _this; // Add other fields as needed

} UserCtx;

// Add other shared types (e.g., IBMNodes, IBMInfo) if needed

#endif // COMMON_H
