/**
 * @file walkingsearch.h
 * @brief Header file for particle location functions using the walking search algorithm.
 *
 * This file contains declarations of functions that implement the walking search algorithm
 * to locate particles within a computational grid. It includes functions for computing
 * distances to cell faces, determining particle positions relative to cells, and updating
 * traversal parameters.
 */

#ifndef WALKINGSEARCH_H
#define WALKINGSEARCH_H

// Include necessary headers
#include <petsc.h>
#include <stdbool.h>
#include <math.h>
#include "common.h"   // Common type definitions
#include "logging.h"  // Logging macros and definitions

// --------------------- Enumerations ---------------------

/**
 * @brief Enumerates the six faces of a cubic cell for distance calculations.
 */
typedef enum {
    LEFT = 0,    /**< Left face (x-) */
    RIGHT,       /**< Right face (x+) */
    BOTTOM,      /**< Bottom face (y-) */
    TOP,         /**< Top face (y+) */
    FRONT,       /**< Front face (z-) */
    BACK,        /**< Back face (z+) */
    NUM_FACES    /**< Total number of faces */
} Face;

// --------------------- Function Declarations ---------------------

/**
 * @brief Computes the signed distance from a point to the plane defined by four other points.
 *
 * This function calculates the signed distance from a given point `p` to the plane defined by points
 * `p1`, `p2`, `p3`, and `p4`. The distance is projected along the normal vector of the plane.
 *
 * @param[in]  p1        First point defining the plane.
 * @param[in]  p2        Second point defining the plane.
 * @param[in]  p3        Third point defining the plane.
 * @param[in]  p4        Fourth point defining the plane (should lie on the same plane as p1, p2, p3).
 * @param[in]  p         The point from which the distance to the plane is calculated.
 * @param[out] d         Pointer to store the computed signed distance.
 * @param[in]  threshold The threshold below which the distance is considered zero.
 *
 * @return void
 */
void ComputeSignedDistanceToPlane(const Cmpnts p1, const Cmpnts p2, const Cmpnts p3, const Cmpnts p4, const Cmpnts p, PetscReal *d, const PetscReal threshold);

/**
 * @brief Computes the signed distances from a point to each face of a cubic cell.
 *
 * This function calculates the signed distances from a specified point `p` to each of the six
 * faces of a cubic cell. The cell is defined by its eight vertices, and the distances are
 * stored in the memory location pointed to by `d`.
 *
 * @param[in]  p         The target point in 3D space.
 * @param[in]  cell      Pointer to a `Cell` structure that defines the cubic cell via its vertices.
 * @param[out] d         A pointer to an array of six `PetscReal` values where the computed signed distances will be stored.
 * @param[in]  threshold A `PetscReal` value that specifies the minimum distance below which a computed distance is considered to be zero.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode CalculateDistancesToCellFaces(const Cmpnts p, const Cell *cell, PetscReal *d, const PetscReal threshold);

/**
 * @brief Determines whether a point is inside, on the boundary, or outside a cell.
 *
 * This function calculates the signed distances from a point `p` to each face of the cell
 * and determines the point's position relative to the cell.
 *
 * @param[in]  p         The point in 3D space to be tested.
 * @param[in]  cell      Pointer to a `Cell` structure representing the cell.
 * @param[out] result    A pointer to an integer that will be set based on the point's position:
 *                        - `0`: Inside the cell
 *                        - `1`: On the boundary of the cell
 *                        - `-1`: Outside the cell
 * @param[in]  threshold The threshold below which the distance is considered zero.
 *
 * @return PetscErrorCode  Returns 0 on success, non-zero on failure.
 */
PetscErrorCode DeterminePointPosition(const Cmpnts p, const Cell *cell, int *result, const PetscReal threshold);

/**
 * @brief Retrieves the coordinates of the eight vertices of a cell based on grid indices.
 *
 * This function populates the `cell` structure with the coordinates of the eight vertices
 * of a hexahedral cell given its grid indices `(idx, idy, idz)`.
 *
 * @param[in]  coor   Pointer to a 3D array of `Cmpnts` structures representing the grid coordinates.
 * @param[in]  idx    The i-index of the cell.
 * @param[in]  idy    The j-index of the cell.
 * @param[in]  idz    The k-index of the cell.
 * @param[out] cell   Pointer to a `Cell` structure to store the cell's vertices.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode GetCellVerticesFromGrid(Cmpnts ***coor, PetscInt idx, PetscInt idy, PetscInt idz, Cell *cell);

/**
 * @brief Initializes traversal parameters for locating a particle.
 *
 * This function sets the initial cell indices based on the local grid boundaries
 * and initializes the traversal step counter.
 *
 * @param[in]  user            Pointer to the user-defined context containing grid information.
 * @param[in]  particle        Pointer to the Particle structure containing its location.
 * @param[out] idx             Pointer to store the initial i-index of the cell.
 * @param[out] idy             Pointer to store the initial j-index of the cell.
 * @param[out] idz             Pointer to store the initial k-index of the cell.
 * @param[out] traversal_steps Pointer to store the initial traversal step count.
 *
 * @return PetscErrorCode     Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InitializeTraversalParameters(UserCtx *user, Particle *particle, PetscInt *idx, PetscInt *idy, PetscInt *idz, PetscInt *traversal_steps);

/**
 * @brief Checks if the current cell indices are within the local grid boundaries.
 *
 * @param[in]  user       Pointer to the user-defined context containing grid information.
 * @param[in]  idx        The i-index of the current cell.
 * @param[in]  idy        The j-index of the current cell.
 * @param[in]  idz        The k-index of the current cell.
 * @param[out] is_within  Pointer to a PetscBool that will be set to PETSC_TRUE if within bounds, else PETSC_FALSE.
 *
 * @return PetscErrorCode  Returns 0 on success, non-zero on failure.
 */
PetscErrorCode CheckCellWithinLocalGrid(UserCtx *user, PetscInt idx, PetscInt idy, PetscInt idz, PetscBool *is_within);

/**
 * @brief Retrieves the coordinates of the eight vertices of the current cell.
 *
 * @param[in]  user  Pointer to the user-defined context containing grid information.
 * @param[in]  idx   The i-index of the current cell.
 * @param[in]  idy   The j-index of the current cell.
 * @param[in]  idz   The k-index of the current cell.
 * @param[out] cell  Pointer to a Cell structure to store the cell's vertices.
 *
 * @return PetscErrorCode  Returns 0 on success, non-zero on failure.
 */
PetscErrorCode RetrieveCurrentCell(UserCtx *user, PetscInt idx, PetscInt idy, PetscInt idz, Cell *cell);

/**
 * @brief Locates the cell within the grid that contains the given particle.
 *
 * This function navigates through cells in a 3D grid to find the cell that contains the specified particle.
 * It uses the signed distances from the particle to the cell faces to determine the direction to move.
 *
 * @param[in]  user     Pointer to the user-defined context containing grid information (DMDA, etc.).
 * @param[in]  particle Pointer to the Particle structure containing its location and identifiers.
 * @param[in]  d         A pointer to an array of six `PetscReal` values that store the
 *                       signed distances from the particle to each face of the cell.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 *
 * @note
 * - Ensure that the `user` and `particle` pointers are not `NULL` before calling this function.
 * - The function assumes that the grid is properly partitioned and that each process has access to its local grid.
 * - The `Particle` structure should have its `loc` field accurately set before calling this function.
 */
PetscErrorCode LocateParticleInGrid(UserCtx *user, Particle *particle, PetscReal* d);

#endif // WALKINGSEARCH_H
