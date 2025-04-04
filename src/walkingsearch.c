// walkingsearch.c

#include "walkingsearch.h"
#include "logging.h"
#include <petsc.h>
#include <stdbool.h>
#include <math.h>

// Define maximum traversal steps to prevent infinite loops
#define MAX_TRAVERSAL 1000
#define DISTANCE_THRESHOLD 1e-14

/**
 * @brief Estimates a characteristic length of the cell for threshold scaling.
 *
 * For a hexahedral cell with vertices cell->vertices[0..7], we approximate
 * the cell size by some measure—e.g. average edge length or diagonal.
 * @param[in]  cell A pointer to the Cell structure
 * @param[out] cellSize A pointer to a PetscReal where the characteristic size is stored.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode GetCellCharacteristicSize(const Cell *cell, PetscReal *cellSize)
{
    PetscErrorCode ierr = 0;
    if (!cell || !cellSize) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, 
        "GetCellCharacteristicSize: Null pointer(s).");

    // A simple approach: compute the average of the distances between
    // all pairs of adjacent vertices that share an edge. That gives
    // a typical measure of cell dimension. For a uniform grid, you
    // could do something simpler.

    PetscInt edges[12][2] = {
        {0,1},{1,2},{2,3},{3,0},  // bottom face edges
        {4,5},{5,6},{6,7},{7,4},  // top face edges
        {0,7},{1,6},{2,5},{3,4}   // vertical edges
    };

    PetscReal totalEdgeLen = 0.0;
    for (PetscInt i=0; i<12; i++) {
        PetscInt vA = edges[i][0];
        PetscInt vB = edges[i][1];
        Cmpnts A = cell->vertices[vA];
        Cmpnts B = cell->vertices[vB];
        PetscReal dx = B.x - A.x;
        PetscReal dy = B.y - A.y;
        PetscReal dz = B.z - A.z;
        PetscReal edgeLen = sqrt(dx*dx + dy*dy + dz*dz);
        totalEdgeLen += edgeLen;
    }

    // Average edge length
    *cellSize = totalEdgeLen / 12.0;

    return ierr;
}

/**
 * @brief Computes the signed distance from a point to the plane defined by four other points.
 *
 * This function calculates the signed distance from a given point `p` to the plane defined by points
 * `p1`, `p2`, `p3`, and `p4`. The distance is projected along the normal vector of the plane.
 * A positive distance indicates that the point lies in the direction of the normal vector,
 * while a negative distance indicates the opposite. If the absolute distance is less than the
 * specified threshold, it is considered to be zero to account for floating-point inaccuracies.
 *
 * The plane is determined by three non-colinear points among `p1`, `p2`, `p3`, and `p4`. The fourth
 * point is used to ensure that all four points lie on the same plane.
 *
 * @param[in]  p1        First point defining the plane.
 * @param[in]  p2        Second point defining the plane.
 * @param[in]  p3        Third point defining the plane.
 * @param[in]  p4        Fourth point defining the plane (should lie on the same plane as p1, p2, p3).
 * @param[in]  p         The point from which the distance to the plane is calculated.
 * @param[out] d         Pointer to store the computed signed distance.
 * @param[in]  threshold The threshold below which the distance is considered zero.
 *
 * @note
 * - Ensure that the four points (`p1`, `p2`, `p3`, `p4`) are not colinear and indeed define a valid plane.
 * - The `threshold` parameter allows flexibility in handling floating-point precision based on application needs.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ComputeSignedDistanceToPlane(const Cmpnts p1, const Cmpnts p2, const Cmpnts p3, const Cmpnts p4, const Cmpnts p, PetscReal *d, const PetscReal threshold)
{
  
   PetscMPIInt rank;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    // Validate output pointer
    if (d == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "ComputeSignedDistanceToPlane - Output pointer 'd' is NULL on rank %d \n",rank);
        return PETSC_ERR_ARG_NULL;
    }

    // Debug: Print points defining the plane and the point p
    LOG_ALLOW(LOCAL,LOG_DEBUG, "ComputeSignedDistanceToPlane - Computing distance for point (%.3f, %.3f, %.3f) to plane defined by:\n", p.x, p.y, p.z);
    LOG_ALLOW(LOCAL,LOG_DEBUG, "  p1: (%.3f, %.3f, %.3f)\n", p1.x, p1.y, p1.z);
    LOG_ALLOW(LOCAL,LOG_DEBUG, "  p2: (%.3f, %.3f, %.3f)\n", p2.x, p2.y, p2.z);
    LOG_ALLOW(LOCAL,LOG_DEBUG, "  p3: (%.3f, %.3f, %.3f)\n", p3.x, p3.y, p3.z);
    LOG_ALLOW(LOCAL,LOG_DEBUG, "  p4: (%.3f, %.3f, %.3f)\n", p4.x, p4.y, p4.z);
    
    // Calculate vectors in the plane
    PetscReal vector1_x = p3.x - p1.x;
    PetscReal vector1_y = p3.y - p1.y;
    PetscReal vector1_z = p3.z - p1.z;

    PetscReal vector2_x = p4.x - p2.x;
    PetscReal vector2_y = p4.y - p2.y;
    PetscReal vector2_z = p4.z - p2.z;

    // Compute the normal vector via cross product
    PetscReal normal_x = vector1_y * vector2_z - vector1_z * vector2_y;
    PetscReal normal_y = -(vector1_x * vector2_z - vector1_z * vector2_x);
    PetscReal normal_z = vector1_x * vector2_y - vector1_y * vector2_x;

    // Compute the magnitude of the normal vector
    PetscReal normal_magnitude = sqrt(normal_x * normal_x + normal_y * normal_y + normal_z * normal_z);

    // Check for degenerate plane
    if (normal_magnitude == 0.0) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "ComputeSignedDistanceToPlane - Degenerate plane detected (zero normal vector) on rank %d \n",rank);

        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,
                "Degenerate plane detected (normal vector is zero).");
        return PETSC_ERR_USER;
    }

    // Normalize the normal vector
    normal_x /= normal_magnitude;
    normal_y /= normal_magnitude;
    normal_z /= normal_magnitude;

    // Calculate the centroid of the four defining points
    PetscReal centroid_x = 0.25 * (p1.x + p2.x + p3.x + p4.x);
    PetscReal centroid_y = 0.25 * (p1.y + p2.y + p3.y + p4.y);
    PetscReal centroid_z = 0.25 * (p1.z + p2.z + p3.z + p4.z);

    // Compute the signed distance from point p to the plane
    *d = (centroid_x - p.x) * normal_x + (centroid_y - p.y) * normal_y + (centroid_z - p.z) * normal_z;

    // Apply threshold to handle floating-point precision
    if (fabs(*d) < threshold) {
        *d = 0.0;
    }

    // Debug: Print the computed distance
    LOG_ALLOW(LOCAL,LOG_DEBUG, "ComputeSignedDistanceToPlane - Signed distance: %.6f\n", *d);

    return 0;
}

/**
 * @brief Computes the signed distances from a point to each face of a cubic cell.
 *
 * This function calculates the signed distances from a specified point `p` to each of the six
 * faces of a cubic cell. The cell is defined by its eight vertices, and the distances are
 * stored in the memory location pointed to by `d`. Each distance corresponds to a specific
 * face of the cell, as enumerated by the `Face` enumeration. A threshold value is used to
 * determine when a distance should be considered effectively zero, accounting for floating-point
 * precision limitations.
 *
 * @param[in]  p         The target point in 3D space for which distances to the cell's faces
 *                       are to be calculated. Represented by the `Cmpnts` structure.
 * @param[in]  cell      A pointer to a `Cell` structure that defines the cubic cell via its
 *                       eight vertices. The vertices must be ordered consistently to
 *                       ensure correct normal directions for each face.
 * @param[out] d         A pointer to an array of six `PetscReal` values where the computed
 *                       signed distances will be stored. Each index in the array corresponds
 *                       to a specific face as defined by the `Face` enumeration:
 *                       - d[LEFT]: Distance to the Left Face of the cell.
 *                       - d[RIGHT]: Distance to the Right Face of the cell.
 *                       - d[BOTTOM]: Distance to the Bottom Face of the cell.
 *                       - d[TOP]: Distance to the Top Face of the cell.
 *                       - d[FRONT]: Distance to the Front Face of the cell.
 *                       - d[BACK]: Distance to the Back Face of the cell.
 * @param[in] threshold A `PetscReal` value that specifies the minimum distance below which
 *                       a computed distance is considered to be zero. This helps in
 *                       mitigating issues arising from floating-point arithmetic inaccuracies.
 *
 * @return PetscErrorCode Returns 0 if the function executes successfully. If an error occurs,
 *                        a non-zero error code is returned, indicating the type of failure.
 *
 * @note
 * - The vertices of the cell must be ordered in a counter-clockwise manner when viewed from
 *   outside the cell. This ordering ensures that the normal vectors computed for each face
 *   point outward, which is essential for correctly determining the sign of the distances.
 * - All four vertices defining a face must lie on a non-degenerate (i.e., non-flat) plane.
 *   Degenerate planes can lead to undefined normal vectors and incorrect distance calculations.
 * - The `threshold` parameter provides flexibility in handling cases where the point `p` is
 *   extremely close to a face, allowing the user to define what constitutes "close enough" to
 *   be treated as zero distance based on the specific requirements of their application.
 */
PetscErrorCode CalculateDistancesToCellFaces(const Cmpnts p, const Cell *cell, PetscReal *d, const PetscReal threshold)
{
    PetscErrorCode ierr;
    // Validate that the 'cell' pointer is not NULL to prevent dereferencing a null pointer.
    if (cell == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "CalculateDistancesToCellFaces - 'cell' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "CalculateDistancesToCellFaces - 'cell' is NULL.");
    }

    // Validate that the 'd' pointer is not NULL to ensure there is memory allocated for distance storage.
    if (d == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "CalculateDistancesToCellFaces - 'd' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "CalculateDistancesToCellFaces - 'd' is NULL.");
    }

    // Compute the signed distance from point 'p' to the BACK face of the cell.
    // The BACK face is defined by vertices 0, 3, 2, and 1, with its normal vector pointing in the -z direction.
    ierr = ComputeSignedDistanceToPlane(
        cell->vertices[0], // Vertex 0
        cell->vertices[3], // Vertex 3
        cell->vertices[2], // Vertex 2
        cell->vertices[1], // Vertex 1
        p,                  // Target point
        &d[BACK],           // Storage location for the BACK face distance
        threshold           // Threshold for zero distance
    );  CHKERRQ(ierr);

    // Compute the signed distance from point 'p' to the FRONT face of the cell.
    // The FRONT face is defined by vertices 4, 7, 6, and 5, with its normal vector pointing in the +z direction.
    ierr = ComputeSignedDistanceToPlane(
        cell->vertices[4], // Vertex 4
        cell->vertices[7], // Vertex 7
        cell->vertices[6], // Vertex 6
        cell->vertices[5], // Vertex 5
        p,                  // Target point
        &d[FRONT],          // Storage location for the FRONT face distance
        threshold           // Threshold for zero distance
    );  CHKERRQ(ierr);

    // Compute the signed distance from point 'p' to the BOTTOM face of the cell.
    // The BOTTOM face is defined by vertices 0, 1, 6, and 7, with its normal vector pointing in the -y direction.
    ierr = ComputeSignedDistanceToPlane(
        cell->vertices[0], // Vertex 0
        cell->vertices[1], // Vertex 1
        cell->vertices[6], // Vertex 6
        cell->vertices[7], // Vertex 7
        p,                  // Target point
        &d[BOTTOM],         // Storage location for the BOTTOM face distance
        threshold           // Threshold for zero distance
    );  CHKERRQ(ierr);

    // Compute the signed distance from point 'p' to the TOP face of the cell.
    // The TOP face is defined by vertices 3, 4, 5, and 2, with its normal vector pointing in the +y direction.
    ierr = ComputeSignedDistanceToPlane(
        cell->vertices[3], // Vertex 3
        cell->vertices[4], // Vertex 4
        cell->vertices[5], // Vertex 5
        cell->vertices[2], // Vertex 2
        p,                  // Target point
        &d[TOP],            // Storage location for the TOP face distance
        threshold           // Threshold for zero distance
    );  CHKERRQ(ierr);

    // Compute the signed distance from point 'p' to the LEFT face of the cell.
    // The LEFT face is defined by vertices 0, 7, 4, and 3, with its normal vector pointing in the -x direction.
    ierr = ComputeSignedDistanceToPlane(
        cell->vertices[0], // Vertex 0
        cell->vertices[7], // Vertex 7
        cell->vertices[4], // Vertex 4
        cell->vertices[3], // Vertex 3
        p,                  // Target point
        &d[LEFT],           // Storage location for the LEFT face distance
        threshold           // Threshold for zero distance
    );  CHKERRQ(ierr);

    // Compute the signed distance from point 'p' to the RIGHT face of the cell.
    // The RIGHT face is defined by vertices 1, 2, 5, and 6, with its normal vector pointing in the +x direction.
    ierr = ComputeSignedDistanceToPlane(
        cell->vertices[1], // Vertex 1
        cell->vertices[2], // Vertex 2
        cell->vertices[5], // Vertex 5
        cell->vertices[6], // Vertex 6
        p,                  // Target point
        &d[RIGHT],          // Storage location for the RIGHT face distance
        threshold           // Threshold for zero distance
    );  CHKERRQ(ierr);

    return 0; // Indicate successful execution of the function.
}

/**
 * @brief Classifies a point based on precomputed face distances.
 *
 * Given an array of six distances d[NUM_FACES] from the point to each face
 * of a hexahedral cell, this function determines:
 *    0 => inside
 *    1 => on a single face
 *    2 => on an edge (2 faces = 0)
 *    3 => on a corner (3+ faces = 0)
 *   -1 => outside
 *
 * @param[in]  d        Array of six face distances.
 * @param[out] result   Pointer to an integer classification: {0,1,2,3,-1}.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode DeterminePointPosition(PetscReal *d, PetscInt *result)
{
    

    // Validate input pointers
    if (d == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "DeterminePointPosition - 'd' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "DeterminePointPosition - Input parameter 'd' is NULL.");
    }
    if (result == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "DeterminePointPosition - 'result' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "DeterminePointPosition - Output parameter 'result' is NULL.");
    }

    // Initialize flags
    PetscBool isInside = PETSC_TRUE;
    PetscBool isOnBoundary = PETSC_FALSE;
    PetscInt IntersectionCount = 0; // Counts the number of intersections of the point with various planes of the cell.

    // Analyze distances to determine position
    for(int i = 0; i < NUM_FACES; i++) {
        if(d[i] <  0.0) {
            isInside = PETSC_FALSE; // Point is outside in at least one direction
        }
        if(d[i] == 0.0) {
            isOnBoundary = PETSC_TRUE; // Point is exactly at least one face
            IntersectionCount++; 
        }
    }

    // Set the result based on flags
    if(isInside) {
        *result = 0; // Inside the cell
        LOG_ALLOW(LOCAL,LOG_DEBUG, "DeterminePointPosition - Particle is inside the cell.\n");
    }
    else if(isOnBoundary) {
      if(IntersectionCount == 1){ 
          *result = 1; // on face
          LOG_ALLOW(LOCAL,LOG_DEBUG, "DeterminePointPosition - Particle is on a face of the cell.\n");
      }
      else if(IntersectionCount == 2){ 
          *result = 2; // on edge
          LOG_ALLOW(LOCAL,LOG_DEBUG, "DeterminePointPosition - Particle is on an edge of the cell.\n");
      }
      else if(IntersectionCount >= 3){ 
          *result = 3; // on corner
          LOG_ALLOW(LOCAL,LOG_DEBUG, "DeterminePointPosition - Particle is on a corner of the cell.\n");
      }
    }
    else {
        *result = -1; // Outside the cell
        LOG_ALLOW(LOCAL,LOG_DEBUG, "DeterminePointPosition - Particle is outside the cell.\n");
    }

    return 0; // Indicate successful execution
}

/**
 * @brief Retrieves the coordinates of the eight vertices of a cell based on grid indices.
 *
 * This function populates the `cell` structure with the coordinates of the eight vertices
 * of a hexahedral cell given its grid indices `(idx, idy, idz)`. The vertex numbering follows
 * the convention depicted in the provided figure:
 *
 * ```
 *               (i,j+1,k)
 *              3---------------------------2(i+1,j+1,k)
 *             /|                          /|
 *            / |                         / |
 *           /  |                        /  |
 *          /   |                       /   |
 *         /    |                      /    |
 *        /     |                     /     |
 *       /      |                    /      |
 *      4-------|------------------5        |
 *      |       |                   |       |
 *      |(i,j+1,k+1)                |       |       
 *      |       |                   |       |
 *      |       0-------------------|-------1 (i+1,j,k)
 *      |       /                   |      /
 *      |      /                    |     /
 *      |     /                     |    /
 *      |    /                      |   /
 *      |   /                       |  /
 *      |  /                        | /
 *      | /                         |/
 *      7---------------------------6
 *      (i,j,k+1)                  (i+1,j,k+1)
 * ```
 *
 * **Vertex Numbering and Positions:**
 * - **Vertex 0**: `(i,   j,   k)`       → `cell->vertices[0]`
 * - **Vertex 1**: `(i+1, j,   k)`       → `cell->vertices[1]`
 * - **Vertex 2**: `(i+1, j+1, k)`       → `cell->vertices[2]`
 * - **Vertex 3**: `(i,   j+1, k)`       → `cell->vertices[3]`
 * - **Vertex 4**: `(i,   j+1, k+1)`     → `cell->vertices[4]`
 * - **Vertex 5**: `(i+1, j+1, k+1)`     → `cell->vertices[5]`
 * - **Vertex 6**: `(i+1, j,   k+1)`     → `cell->vertices[6]`
 * - **Vertex 7**: `(i,   j,   k+1)`     → `cell->vertices[7]`
 *
 * @param[in]  coor   Pointer to a 3D array of `Cmpnts` structures representing the grid coordinates.
 * @param[in]  idx    The i-index of the cell.
 * @param[in]  idy    The j-index of the cell.
 * @param[in]  idz    The k-index of the cell.
 * @param[out] cell   Pointer to a `Cell` structure to store the cell's vertices.
 *
 * @return PetscErrorCode Returns 0 to indicate successful execution. Non-zero on failure.
 *
 * @note
 * - It is assumed that the grid indices `(idx, idy, idz)` are within valid bounds.
 * - The `coor` array should be properly initialized with accurate coordinates for all grid points.
 * - The function assumes unit spacing between grid points. Modify the coordinate assignments if your grid spacing differs.
 * - The boundary checks have been removed as they are performed before calling this function.
 */
PetscErrorCode GetCellVerticesFromGrid(Cmpnts ***coor, PetscInt idx, PetscInt idy, PetscInt idz,
                                       Cell *cell)
{

    // Validate input pointers
    if (coor == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "GetCellVerticesFromGrid - 'coor' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "GetCellVerticesFromGrid - Input array 'coor' is NULL.");
    }
    if (cell == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "GetCellVerticesFromGrid - 'cell' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "GetCellVerticesFromGrid - Output parameter 'cell' is NULL.");
    }

    // Assign vertices based on grid indices matching the figure
    // Vertex numbering follows the convention depicted in the figure
    cell->vertices[0] = coor[idz][idy][idx];         // Vertex 0: (i,   j,   k)
    cell->vertices[1] = coor[idz][idy][idx+1];       // Vertex 1: (i+1, j,   k)
    cell->vertices[2] = coor[idz][idy+1][idx+1];     // Vertex 2: (i+1, j+1, k)
    cell->vertices[3] = coor[idz][idy+1][idx];       // Vertex 3: (i,   j+1, k)
    cell->vertices[4] = coor[idz+1][idy+1][idx];     // Vertex 4: (i,   j+1, k+1)
    cell->vertices[5] = coor[idz+1][idy+1][idx+1];   // Vertex 5: (i+1, j+1, k+1)
    cell->vertices[6] = coor[idz+1][idy][idx+1];     // Vertex 6: (i+1, j,   k+1)
    cell->vertices[7] = coor[idz+1][idy][idx];       // Vertex 7: (i,   j,   k+1)

    LOG_ALLOW(LOCAL,LOG_DEBUG, "GetCellVerticesFromGrid - Retrieved vertices for cell (%ld, %ld, %ld).\n", idx, idy, idz);

    return 0; // Indicate successful execution
}

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
PetscErrorCode InitializeTraversalParameters(UserCtx *user, Particle *particle, PetscInt *idx, PetscInt *idy, PetscInt *idz, PetscInt *traversal_steps)
{
    PetscErrorCode ierr;
    DMDALocalInfo info;

    // Validate input pointers
    if (user == NULL || particle == NULL || idx == NULL || idy == NULL || idz == NULL || traversal_steps == NULL) {
      LOG_ALLOW(GLOBAL,LOG_ERROR, "InitializeTraversalParameters - One or more input pointers are NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "InitializeTraversalParameters - One or more input pointers are NULL.");
    }

    // Get grid information
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

    // Option 1: Start from the first cell in the local grid
    *idx = info.xs;
    *idy = info.ys;
    *idz = info.zs;

    // Option 2: Compute initial indices based on particle location (optional)
    // Assuming unit grid spacing
    /*
    *idx = (PetscInt)floor(particle->loc.x);
    *idy = (PetscInt)floor(particle->loc.y);
    *idz = (PetscInt)floor(particle->loc.z);
    */

    // Initialize traversal step counter
    *traversal_steps = 0;

    LOG_ALLOW(LOCAL,LOG_INFO, "InitializeTraversalParameters - Starting traversal at cell (%ld, %ld, %ld).\n", *idx, *idy, *idz);

    return 0;
}

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
PetscErrorCode CheckCellWithinLocalGrid(UserCtx *user, PetscInt idx, PetscInt idy, PetscInt idz, PetscBool *is_within)
{
    PetscErrorCode ierr;
    DMDALocalInfo info;

    // Validate input pointers
    if (user == NULL || is_within == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "CheckCellWithinLocalGrid - One or more input pointers are NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "CheckCellWithinLocalGrid - One or more input pointers are NULL.");
    }

    // Get grid information
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

    // Check if indices are within bounds
    if (idx >= info.xs && idx < (info.xs + info.xm - 1) &&
        idy >= info.ys && idy < (info.ys + info.ym - 1) &&
        idz >= info.zs && idz < (info.zs + info.zm - 1)) {
        *is_within = PETSC_TRUE;
    }
    else {
        *is_within = PETSC_FALSE;
    }

    LOG_ALLOW(LOCAL,LOG_DEBUG, "CheckCellWithinLocalGrid - Cell (%ld, %ld, %ld) is %s the local grid.\n",
        idx, idy, idz, (*is_within) ? "within" : "outside");

    return 0;
}

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
PetscErrorCode RetrieveCurrentCell(UserCtx *user, PetscInt idx, PetscInt idy, PetscInt idz, Cell *cell)
{
    PetscErrorCode ierr;
    Vec Coor;
    Cmpnts ***coor_array;
    PetscMPIInt rank;

    // Validate input pointers
    if (user == NULL || cell == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "RetrieveCurrentCell - One or more input pointers are NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "RetrieveCurrentCell - One or more input pointers are NULL.");
    }

    // Get local coordinates
    ierr = DMGetCoordinatesLocal(user->da, &Coor); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, Coor, &coor_array); CHKERRQ(ierr);

    // Get the current cell's vertices
    ierr = GetCellVerticesFromGrid(coor_array, idx, idy, idz, cell); CHKERRQ(ierr);

    // Restore array
    ierr = DMDAVecRestoreArrayRead(user->fda, Coor, &coor_array); CHKERRQ(ierr);

    // Get MPI rank for debugging
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Debug: Print cell vertices
    LOG_ALLOW(LOCAL,LOG_DEBUG, "RetrieveCurrentCell - Cell (%ld, %ld, %ld) vertices \n", idx, idy, idz);
    ierr = LOG_CELL_VERTICES(cell, rank); CHKERRQ(ierr);

    return 0;
}

/**
 * @brief Determines the spatial relationship of a particle relative to a cubic cell.
 *
 * This function evaluates whether a particle located at a specific point `p` in 3D space
 * is positioned inside the cell, on the boundary of the cell, or outside the cell. The
 * determination is based on the signed distances from the particle to each of the six
 * faces of the cell. The function utilizes the computed distances to ascertain the particle's
 * position with respect to the cell boundaries, considering a threshold to account for
 * floating-point precision.
 *
 * @param[in]  cell      A pointer to a `Cell` structure that defines the cubic cell via its
 *                       vertices. The cell's geometry is essential for accurately computing
 *                       the distances to each face.
 * @param[in]  d         A pointer to an array of six `PetscReal` values that store the
 *                       signed distances from the particle to each face of the cell. These
 *                       distances are typically computed using the `CalculateDistancesToCellFaces`
 *                       function.
 * @param[in]  p         The location of the particle in 3D space, represented by the `Cmpnts`
 *                       structure. This point is the reference for distance calculations to the
 *                       cell's faces.
 * @param[out] position  A pointer to an integer that will be set based on the particle's position
 *                       relative to the cell:
 *                       - `0`: The particle is inside the cell.
 *                       - `1`: The particle is on the boundary of the cell.
 *                       - `-1`: The particle is outside the cell.
 * @param[in]  threshold A `PetscReal` value that defines the minimum distance below which a
 *                       computed distance is considered to be zero. This threshold helps in
 *                       mitigating inaccuracies due to floating-point arithmetic, especially
 *                       when determining if the particle lies exactly on the boundary.
 *
 * @return PetscErrorCode Returns `0` if the function executes successfully. If an error occurs,
 *                        a non-zero error code is returned, indicating the type of failure.
 *
 * @note
 * - It is assumed that the `d` array has been properly allocated and contains valid distance
 *   measurements before calling this function.
 * - The function relies on `CalculateDistancesToCellFaces` to accurately compute the signed
 *   distances to each face. Any inaccuracies in distance calculations can affect the
 *   determination of the particle's position.
 * - The `threshold` parameter should be chosen based on the specific precision requirements
 *   of the application to balance between sensitivity and robustness against floating-point
 *   errors.
 * - The function includes a debug statement that prints the face distances, which can be useful
 *   for verifying the correctness of distance computations during development or troubleshooting.
 */
PetscErrorCode EvaluateParticlePosition(const Cell *cell, PetscReal *d, const Cmpnts p, PetscInt *position, const PetscReal threshold)
{
    PetscErrorCode ierr;
    PetscReal cellSize;
    PetscReal cellThreshold;

    // Validate input pointers to ensure they are not NULL, preventing potential segmentation faults.
    if (cell == NULL || position == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "EvaluateParticlePosition - One or more input pointers are NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "EvaluateParticlePosition - One or more input pointers are NULL.");
    }
  
    // Compute a local cell size
    ierr = GetCellCharacteristicSize(cell, &cellSize); CHKERRQ(ierr);

    // scale base threshold by cell size (ex: if threshold = 1e-6 and cell size = 0.01, cellThreshold = 1e-8)
    cellThreshold = threshold*cellSize; 

    // Invoke the function to calculate signed distances from the particle to each face of the cell.
    // The distances are stored in the array pointed to by 'd'.
    ierr = CalculateDistancesToCellFaces(p, cell, d, cellThreshold); CHKERRQ(ierr);
    CHKERRQ(ierr); // Check for errors in distance calculation.


    // Catch degenerate-plane error manually:
    if (ierr == PETSC_ERR_USER) {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "EvaluateParticlePosition - Skipping cell due to degenerate face.\n");
        // We can set *position = -1 here
        *position = -1; // treat as outside
        return 0; // not a fatal error, just skip
    } else {
        CHKERRQ(ierr);
    }


    // Debugging output: Print the computed distances to each face for verification purposes.
    LOG_ALLOW(LOCAL,LOG_DEBUG, "EvaluateParticlePosition - Face Distances:\n");
    ierr = LOG_FACE_DISTANCES(d); 
    CHKERRQ(ierr); // Check for errors in printing distances.

    // Determine the particle's position relative to the cell based on the computed distances.
    // The function sets the value pointed to by 'position' accordingly:
    // 0 for inside, 1 (2 or 3) for on the boundary, and -1 for outside.
    ierr = DeterminePointPosition(d,position); 
    CHKERRQ(ierr); // Check for errors in position determination.

    return 0; // Indicate successful execution of the function.
}

/**
 * @brief Updates the cell indices based on the signed distances to each face.
 *
 * This function modifies the cell indices (`idx`, `idy`, `idz`) to move towards the direction
 * where the particle is likely to be located, based on positive distances indicating
 * that the particle is outside in that particular direction.
 *
 * @param[in]  d    An array of six `PetscReal` values representing the signed distances to each face:
 *                  - d[LEFT]: Left Face
 *                  - d[RIGHT]: Right Face
 *                  - d[BOTTOM]: Bottom Face
 *                  - d[TOP]: Top Face
 *                  - d[FRONT]: Front Face
 *                  - d[BACK]: Back Face
 * @param[out] idx  Pointer to the i-index of the cell to be updated.
 * @param[out] idy  Pointer to the j-index of the cell to be updated.
 * @param[out] idz  Pointer to the k-index of the cell to be updated.
 * @param[in]  info DMDALocalInfo structure that holds local & global domain bounds.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode UpdateCellIndicesBasedOnDistances( PetscReal d[NUM_FACES], PetscInt *idx, PetscInt *idy, PetscInt *idz, DMDALocalInfo *info)
{
    PetscInt cxm,cxs;  // maximum & minimum cell ID in x
    PetscInt cym,cys;  // maximum & minimum cell ID in y
    PetscInt czm,czs;  // maximum & minimum cell ID in z

    cxs = info->xs; cxm = cxs + info->xm - 2;
    cys = info->ys; cym = cys + info->ym - 2;
    czs = info->zs; czm = czs + info->zm - 2; 
    
    // Validate input pointers
    if (d == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "UpdateCellIndicesBasedOnDistances - 'd' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UpdateCellIndicesBasedOnDistances - Input array 'd' is NULL.");
    }
    if (idx == NULL || idy == NULL || idz == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "UpdateCellIndicesBasedOnDistances - One or more index pointers are NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UpdateCellIndicesBasedOnDistances - One or more index pointers are NULL.");
    }

    // Debug: Print current face distances
    LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Current Face Distances:\n");
    LOG_FACE_DISTANCES(d);

    // Update k-direction based on FRONT and BACK distances
    if (d[FRONT] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Condition met: d[FRONT] < 0.0, incrementing idz.\n");
        (*idz) += 1;
    }
    else if(d[BACK] < 0.0){
      LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Condition met: d[BACK] < 0.0, decrementing idz.\n");
        (*idz) -= 1;
    }

    // Update i-direction based on LEFT and RIGHT distances
    if (d[LEFT] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Condition met: d[LEFT] < 0.0, decrementing idx.\n");
        (*idx) -= 1;
    }
    else if (d[RIGHT] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Condition met: d[RIGHT] < 0.0, incrementing idx.\n");
        (*idx) += 1;
    }

    // Update j-direction based on BOTTOM and TOP distances
    if (d[BOTTOM] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Condition met: d[BOTTOM] < 0.0, decrementing idy.\n");
        (*idy) -= 1;
    }
    else if (d[TOP] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Condition met: d[TOP] < 0.0, incrementing idy.\n");
        (*idy) += 1;
    }

    // The 'cell' corners you can reference go from [xs .. xs+xm-1], but
    // to form a valid cell in x, you need (idx+1) in range, so max is (xs+xm-2).
    *idx = PetscMax(cxs,               PetscMin(*idx, cxm));
    *idy = PetscMax(cys,               PetscMin(*idy, cym));
    *idz = PetscMax(czs,               PetscMin(*idz, czm));

    LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Updated Indices after clamping (inside domain bounds)  - idx, idy, idz: %ld, %ld, %ld\n", *idx, *idy, *idz);

    return 0; // Indicate successful execution
}

/**
 * @brief Finalizes the traversal by reporting the results.
 *
 * This function prints the outcome of the traversal, indicating whether the particle
 * was found within a cell or not, and updates the particle's cell indices accordingly.
 *
 * @param[in]  user           Pointer to the user-defined context containing grid information.
 * @param[out] particle       Pointer to the Particle structure to update with cell indices.
 * @param[in]  traversal_steps The number of traversal steps taken.
 * @param[in]  cell_found      Flag indicating whether the particle was found within a cell.
 * @param[in]  idx             The i-index of the found cell.
 * @param[in]  idy             The j-index of the found cell.
 * @param[in]  idz             The k-index of the found cell.
 *
 * @return PetscErrorCode     Returns 0 on success, non-zero on failure.
 */
PetscErrorCode FinalizeTraversal(UserCtx *user, Particle *particle, PetscInt traversal_steps, PetscBool cell_found, PetscInt idx, PetscInt idy, PetscInt idz)
{
    // Validate input pointers
    if (user == NULL || particle == NULL) {
      LOG_ALLOW(GLOBAL,LOG_ERROR, "FinalizeTraversal - One or more input pointers are NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "FinalizeTraversal - One or more input pointers are NULL.");
    }

    if (cell_found) {
      LOG_ALLOW(LOCAL,LOG_INFO, "FinalizeTraversal - Particle located in cell (%ld, %ld, %ld) after %ld traversal steps.\n",
            idx, idy, idz, traversal_steps);
    }
    else {
      LOG_ALLOW(LOCAL,LOG_WARNING, "FinalizeTraversal - Particle could not be located within the grid after %ld traversal steps.\n", (PetscInt)traversal_steps);
        particle->cell[0] = -1;
        particle->cell[1] = -1;
        particle->cell[2] = -1;
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "FinalizeTraversal - Completed final traversal sync across all ranks.\n");


    return 0;
}

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
PetscErrorCode LocateParticleInGrid(UserCtx *user, Particle *particle, PetscReal *d)
{
    PetscErrorCode ierr;
    PetscInt idx, idy, idz;
    PetscInt traversal_steps;
    PetscBool Cell_found = PETSC_FALSE;
    Cmpnts p = particle->loc;
    Cell current_cell;
    const PetscReal threshold = DISTANCE_THRESHOLD ;
    DMDALocalInfo info;
    PetscInt repeatedIndexCount = 0;
    PetscInt prevIdx = PETSC_MIN_INT, prevIdy = PETSC_MIN_INT, prevIdz = PETSC_MIN_INT;
    
    //   LOG_FUNC_TIMER_BEGIN_EVENT(EVENT_Individualwalkingsearch,LOCAL);    

    // Retrieve local DMDA info (which holds xs, xm, ys, ym, zs, zm)
    ierr = DMDAGetLocalInfo(user->da,&info); CHKERRQ(ierr);
    
    // Initialize traversal parameters
    ierr = InitializeTraversalParameters(user, particle, &idx, &idy, &idz, &traversal_steps); CHKERRQ(ierr);

    // Traverse the grid to locate the particle
    while (!Cell_found && traversal_steps < MAX_TRAVERSAL) {
        traversal_steps++;

        // Detect if we haven't changed indices from the last iteration
        if (idx == prevIdx && idy == prevIdy && idz == prevIdz) {
            repeatedIndexCount++;
            if (repeatedIndexCount > 3) {
                // We toggled or got stuck in the same cell multiple times
                LOG_ALLOW(LOCAL, LOG_WARNING,
                "LocateParticleInGrid - Toggling or repeated index detected at cell (%ld,%ld,%ld); breaking.\n",
                    idx, idy, idz);
                break;
            }
        } else {
            // We moved to a new cell index, so reset the counter
            repeatedIndexCount = 0;
            prevIdx = idx;
            prevIdy = idy;
            prevIdz = idz;
        }

        // Check if current cell is within the local grid
        PetscBool is_within;
        ierr = CheckCellWithinLocalGrid(user, idx, idy, idz, &is_within); CHKERRQ(ierr);
        if (!is_within) {
	        LOG_ALLOW(LOCAL,LOG_WARNING, "LocateParticleInGrid - Particle is outside the local grid boundaries at cell (%ld, %ld, %ld).\n", idx, idy, idz);
            break;
        }

        // Retrieve the current cell's vertices
        ierr = RetrieveCurrentCell(user, idx, idy, idz, &current_cell); CHKERRQ(ierr);

        // Evaluate the particle's position relative to the current cell
        PetscInt position;
        ierr = EvaluateParticlePosition(&current_cell, d, p, &position, threshold); CHKERRQ(ierr);
	
	// Log every 10th iteration of evaluation
	LOG_LOOP_ALLOW(GLOBAL, LOG_DEBUG, traversal_steps, 10,
		       "LocateParticleInGrid - At traversal step %ld, evaluated particle position relative to cell (%ld, %ld, %ld): position=%ld.\n",
		       traversal_steps, idx, idy, idz, position);


        if (position == 0) { // Inside the cell
            Cell_found = PETSC_TRUE;
            particle->cell[0] = idx;
            particle->cell[1] = idy;
            particle->cell[2] = idz;
            LOG_ALLOW(LOCAL,LOG_INFO, "LocateParticleInGrid - Particle found in cell (%ld, %ld, %ld).\n", idx, idy, idz);
            break;
        }
        else if (position >= 1) { // On boundary (face,edge or corner) [ can be expanded for specific cases if necessary by having conditions position==1,2 or 3]
           // Depending on application, decide whether to consider it inside or check neighbors
           // Here, we treat it as inside
           Cell_found = PETSC_TRUE;
           particle->cell[0] = idx;
           particle->cell[1] = idy;
           particle->cell[2] = idz;
           LOG_ALLOW(LOCAL,LOG_INFO, "LocateParticleInGrid - Particle is on the boundary of cell (%ld, %ld, %ld).\n", idx, idy, idz);
           break;
        }
        else { // Outside the cell
            // Update cell indices based on positive distances
	  ierr = UpdateCellIndicesBasedOnDistances(d, &idx, &idy, &idz,&info); CHKERRQ(ierr);
        }
    } // !cell_found 

    // Finalize traversal by reporting the results
    ierr = FinalizeTraversal(user, particle, traversal_steps, Cell_found, idx, idy, idz); CHKERRQ(ierr);

    LOG_ALLOW_SYNC(GLOBAL,LOG_INFO, "LocateParticleInGrid - Finalized particle search across all ranks.\n");


    //  LOG_FUNC_TIMER_END_EVENT(EVENT_Individualwalkingsearch,LOCAL);

    return 0;
}


