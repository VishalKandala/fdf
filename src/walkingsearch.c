 // walkingsearch.c

#include "walkingsearch.h"
#include "logging.h"
#include <petsc.h>
#include <stdbool.h>
#include <math.h>

// Define maximum traversal steps to prevent infinite loops
#define MAX_TRAVERSAL 1000
#define DISTANCE_THRESHOLD 1e-14
#define REPEAT_COUNT_THRESHOLD 5 

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
 * @brief Computes the signed distance from a point to the plane approximating a quadrilateral face.
 *
 * This function calculates the signed distance from a given point `p_target` to the plane
 * approximating a quadrilateral face defined by points `v1`, `v2`, `v3`, and `v4`.
 * The plane's initial normal is determined using two edge vectors emanating from `v1`
 * (i.e., `v2-v1` and `v4-v1`).
 *
 * **Normal Orientation:**
 * The initial normal's orientation is checked against the cell's interior.
 * A vector from the face centroid to the cell centroid (`vec_face_to_cell_centroid`) is computed.
 * If the dot product of the initial normal and `vec_face_to_cell_centroid` is positive,
 * it means the initial normal is pointing towards the cell's interior (it's an inward normal
 * relative to the cell). In this case, the normal is flipped to ensure it points outward.
 *
 * **Signed Distance Convention:**
 * The signed distance is the projection of the vector from `p_target` to the face's centroid
 * onto the (now guaranteed) outward-pointing plane's normal vector.
 *   - `d_signed < 0`: `p_target` is "outside" and "beyond" the face, in the direction of the outward normal.
 *                     This is the primary case indicating the particle has crossed this face.
 *   - `d_signed > 0`: `p_target` is "outside" but "behind" the face (on the same side as the cell interior
 *                     relative to this face plane).
 *   - `d_signed = 0`: `p_target` is on the face plane (within threshold).
 *
 * @param[in]  v1            First vertex defining the face (used as origin for initial normal calculation).
 * @param[in]  v2            Second vertex defining the face (used for first edge vector: v2-v1).
 * @param[in]  v3            Third vertex defining the face (used for face centroid calculation).
 * @param[in]  v4            Fourth vertex defining the face (used for second edge vector: v4-v1).
 * @param[in]  cell_centroid The centroid of the 8-vertex cell to which this face belongs.
 * @param[in]  p_target      The point from which the distance to the plane is calculated.
 * @param[out] d_signed      Pointer to store the computed signed distance.
 * @param[in]  threshold     The threshold below which the absolute distance is considered zero.
 *
 * @note
 * - The order of vertices `v1, v2, v3, v4` is used for calculating the face centroid.
 *   The order of `v1, v2, v4` is used for the initial normal calculation.
 * - The `cell_centroid` is crucial for robustly determining the outward direction of the normal.
 *
 * @return PetscErrorCode Returns 0 on success, PETSC_ERR_ARG_NULL if `d_signed` is NULL,
 *                        or PETSC_ERR_USER if a degenerate plane (near-zero normal) is detected.
 */
PetscErrorCode ComputeSignedDistanceToPlane(const Cmpnts v1, const Cmpnts v2, const Cmpnts v3, const Cmpnts v4,
                                            const Cmpnts cell_centroid,const Cmpnts p_target,
					    PetscReal *d_signed, const PetscReal threshold)
{
    PetscErrorCode ierr = 0;
    PetscMPIInt rank;

    PetscFunctionBeginUser;

    if (d_signed == NULL) {
        ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank); CHKERRQ(ierr);
        LOG_ALLOW(LOCAL, LOG_ERROR, "ComputeSignedDistanceToPlane - Output pointer 'd_signed' is NULL on rank %d.\n", rank);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output pointer 'd_signed' must not be NULL.");
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "ComputeSignedDistanceToPlane - Target point: (%.6e, %.6e, %.6e)\n", p_target.x, p_target.y, p_target.z);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Cell Centroid: (%.6e, %.6e, %.6e)\n", cell_centroid.x, cell_centroid.y, cell_centroid.z);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Face Vertices:\n");
    LOG_ALLOW(LOCAL, LOG_DEBUG, "    v1: (%.6e, %.6e, %.6e)\n", v1.x, v1.y, v1.z);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "    v2: (%.6e, %.6e, %.6e)\n", v2.x, v2.y, v2.z);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "    v3: (%.6e, %.6e, %.6e)\n", v3.x, v3.y, v3.z);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "    v4: (%.6e, %.6e, %.6e)\n", v4.x, v4.y, v4.z);

    // --- Calculate Edge Vectors for Initial Normal Computation ---
    PetscReal edge1_x = v2.x - v1.x;
    PetscReal edge1_y = v2.y - v1.y;
    PetscReal edge1_z = v2.z - v1.z;

    PetscReal edge2_x = v4.x - v1.x;
    PetscReal edge2_y = v4.y - v1.y;
    PetscReal edge2_z = v4.z - v1.z;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Edge1 (v2-v1): (%.6e, %.6e, %.6e)\n", edge1_x, edge1_y, edge1_z);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Edge2 (v4-v1): (%.6e, %.6e, %.6e)\n", edge2_x, edge2_y, edge2_z);

    // --- Compute Initial Normal Vector (Cross Product: edge1 x edge2) ---
    PetscReal normal_x_initial = edge1_y * edge2_z - edge1_z * edge2_y;
    PetscReal normal_y_initial = edge1_z * edge2_x - edge1_x * edge2_z;
    PetscReal normal_z_initial = edge1_x * edge2_y - edge1_y * edge2_x;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Initial Raw Normal (edge1 x edge2): (%.6e, %.6e, %.6e)\n", normal_x_initial, normal_y_initial, normal_z_initial);

    PetscReal normal_magnitude = sqrt(normal_x_initial * normal_x_initial +
                                   normal_y_initial * normal_y_initial +
                                   normal_z_initial * normal_z_initial);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Initial Normal Magnitude: %.6e\n", normal_magnitude);

    if (normal_magnitude < 1.0e-12) {
        ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank); CHKERRQ(ierr);
        LOG_ALLOW(LOCAL, LOG_ERROR,
                  "ComputeSignedDistanceToPlane - Degenerate plane detected on rank %d. Normal magnitude (%.3e) is too small.\n",
                  rank, normal_magnitude);
        LOG_ALLOW(LOCAL, LOG_ERROR, "  Offending vertices for normal: v1(%.3e,%.3e,%.3e), v2(%.3e,%.3e,%.3e), v4(%.3e,%.3e,%.3e)\n",
                  v1.x,v1.y,v1.z, v2.x,v2.y,v2.z, v4.x,v4.y,v4.z);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Degenerate plane detected (normal vector is near zero).");
    }

    // --- Calculate the Centroid of the Four Face Vertices (v1, v2, v3, v4) ---
    PetscReal face_centroid_x = 0.25 * (v1.x + v2.x + v3.x + v4.x);
    PetscReal face_centroid_y = 0.25 * (v1.y + v2.y + v3.y + v4.y);
    PetscReal face_centroid_z = 0.25 * (v1.z + v2.z + v3.z + v4.z);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Face Centroid: (%.6e, %.6e, %.6e)\n", face_centroid_x, face_centroid_y, face_centroid_z);

    // --- Orient the Normal to Point Outward from the Cell ---
    // Vector from face centroid to cell centroid
    PetscReal vec_fc_to_cc_x = cell_centroid.x - face_centroid_x;
    PetscReal vec_fc_to_cc_y = cell_centroid.y - face_centroid_y;
    PetscReal vec_fc_to_cc_z = cell_centroid.z - face_centroid_z;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Vec (FaceCentroid -> CellCentroid): (%.6e, %.6e, %.6e)\n", vec_fc_to_cc_x, vec_fc_to_cc_y, vec_fc_to_cc_z);

    // Dot product of initial normal with vector from face centroid to cell centroid
    PetscReal dot_prod_orientation = normal_x_initial * vec_fc_to_cc_x +
                                     normal_y_initial * vec_fc_to_cc_y +
                                     normal_z_initial * vec_fc_to_cc_z;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Dot Product for Orientation (N_initial . Vec_FC_to_CC): %.6e\n", dot_prod_orientation);

    PetscReal normal_x = normal_x_initial;
    PetscReal normal_y = normal_y_initial;
    PetscReal normal_z = normal_z_initial;

    // If dot_prod_orientation > 0, initial normal points towards cell interior, so flip it.
    if (dot_prod_orientation > 1.0e-9) { // Use a small epsilon to avoid issues if dot product is extremely close to zero
        normal_x = -normal_x_initial;
        normal_y = -normal_y_initial;
        normal_z = -normal_z_initial;
        LOG_ALLOW(LOCAL, LOG_DEBUG, "  Initial normal was inward (dot_prod > 0). Flipped normal.\n");
    } else if (dot_prod_orientation == 0.0 && normal_magnitude > 1e-12) {
        // This case is ambiguous or face plane contains cell centroid.
        // This might happen for highly symmetric cells or if face_centroid IS cell_centroid (e.g. 2D cell).
        // For now, we keep the original normal direction based on v1,v2,v4 ordering.
        // A more robust solution for this edge case might be needed if it occurs often.
        ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank); CHKERRQ(ierr);
         LOG_ALLOW(LOCAL, LOG_WARNING, "ComputeSignedDistanceToPlane - Rank %d: Dot product for normal orientation is zero. Normal direction might be ambiguous. Keeping initial normal direction from (v2-v1)x(v4-v1).\n", rank);
    }
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Oriented Raw Normal: (%.6e, %.6e, %.6e)\n", normal_x, normal_y, normal_z);


    // --- Normalize the (Now Outward-Pointing) Normal Vector ---
    // Note: normal_magnitude was calculated from initial normal.
    // If we flipped, magnitude is the same.
    normal_x /= normal_magnitude;
    normal_y /= normal_magnitude;
    normal_z /= normal_magnitude;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Normalized Outward Normal: (%.6f, %.6f, %.6f)\n", normal_x, normal_y, normal_z);


    // --- Compute Vector from Target Point to Face Centroid ---
    PetscReal vec_p_to_fc_x = face_centroid_x - p_target.x;
    PetscReal vec_p_to_fc_y = face_centroid_y - p_target.y;
    PetscReal vec_p_to_fc_z = face_centroid_z - p_target.z;

    // --- Compute Signed Distance ---
    *d_signed = vec_p_to_fc_x * normal_x +
                vec_p_to_fc_y * normal_y +
                vec_p_to_fc_z * normal_z;

    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Raw Signed Distance (using outward normal): %.15e\n", *d_signed);

    // --- Apply Threshold ---
    if (fabs(*d_signed) < threshold) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "  Distance %.15e is less than threshold %.1e. Setting to 0.0.\n", *d_signed, threshold);
        *d_signed = 0.0;
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "ComputeSignedDistanceToPlane - Final Signed Distance: %.15e\n", *d_signed);

    PetscFunctionReturn(0);
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

    // Compute the centroid of the entire cell
    Cmpnts cell_centroid = {0.0, 0.0, 0.0};
    for (int i = 0; i < 8; ++i) {
      cell_centroid.x += cell->vertices[i].x;
      cell_centroid.y += cell->vertices[i].y;
      cell_centroid.z += cell->vertices[i].z;
    }
    cell_centroid.x /= 8.0;
    cell_centroid.y /= 8.0;
    cell_centroid.z /= 8.0;

    LOG_ALLOW(LOCAL,LOG_DEBUG, "CalculateDistancesToCellFaces - Cell Centroid: (%.6e, %.6e, %.6e)\n",
	      cell_centroid.x, cell_centroid.y, cell_centroid.z);    

    
    // Compute the signed distance from point 'p' to the BACK face of the cell.
    // The BACK face is defined by vertices 0, 3, 2, and 1, with its normal vector pointing in the -z direction.
    ierr = ComputeSignedDistanceToPlane(
        cell->vertices[0], // Vertex 0
        cell->vertices[3], // Vertex 3
        cell->vertices[2], // Vertex 2
        cell->vertices[1], // Vertex 1
	cell_centroid,     // Cell centroid 
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
	cell_centroid,     // Cell centroid 
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
	cell_centroid,     // Cell centroid 
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
	cell_centroid,     // Cell centroid 
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
	cell_centroid,     // Cell centroid 
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
	cell_centroid,     // Cell centroid 
        p,                  // Target point
        &d[RIGHT],          // Storage location for the RIGHT face distance
        threshold           // Threshold for zero distance
    );  CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_DEBUG, "CalculateDistancesToCellFaces - Populated d: "
	      "d[LEFT=%d]=%.3e, d[RIGHT=%d]=%.3e, d[BOTTOM=%d]=%.3e, "
	      "d[TOP=%d]=%.3e, d[FRONT=%d]=%.3e, d[BACK=%d]=%.3e\n",
	      LEFT, d[LEFT], RIGHT, d[RIGHT], BOTTOM, d[BOTTOM],
	      TOP, d[TOP], FRONT, d[FRONT], BACK, d[BACK]);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "CalculateDistancesToCellFaces - Raw d: "
	      "[%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]\n",
	      d[0], d[1], d[2], d[3], d[4], d[5]);
    
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

    LOG_ALLOW(LOCAL,LOG_DEBUG, "GetCellVerticesFromGrid - Retrieved vertices for cell (%d, %d, %d).\n", idx, idy, idz);

    return 0; // Indicate successful execution
}

/**
 * @brief Initializes traversal parameters for locating a particle.
 *
 * This function sets the initial cell indices for the grid search.
 * - If the particle has valid previous cell indices (not -1,-1,-1), the search
 *   starts from that previous cell.
 * - Otherwise (e.g., first time locating the particle), the search starts
 *   from the beginning corner (xs, ys, zs) of the local grid domain owned
 *   by the current process.
 * It also initializes the traversal step counter.
 *
 * @param[in]  user            Pointer to the user-defined context containing grid information.
 * @param[in]  particle        Pointer to the Particle structure containing its location and previous cell (if any).
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

    // Get grid information (needed for fallback start and potentially validation)
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

    // --- Check if the particle has a valid previous cell ID ---
    // Assuming particle->cell stores the *global* indices
    if (particle->cell[0] >= 0 && particle->cell[1] >= 0 && particle->cell[2] >= 0) {
        // Previous cell ID is valid, start search from there.
        *idx = particle->cell[0];
        *idy = particle->cell[1];
        *idz = particle->cell[2];

        // Optional: Add a check/clamp to ensure the previous cell ID is at least
        // within reasonable global bounds, though LocateParticleInGrid should handle
        // out-of-local-bounds cases.
        // Example clamp (adjust based on your global grid indexing, 0..IM-1):
        // *idx = PetscMax(0, PetscMin(info.mx - 2, *idx));
        // *idy = PetscMax(0, PetscMin(info.my - 2, *idy));
        // *idz = PetscMax(0, PetscMin(info.mz - 2, *idz));

        LOG_ALLOW(LOCAL,LOG_DEBUG, "InitializeTraversalParameters - Particle %lld has previous cell (%d, %d, %d). Starting search there.\n",
                  (long long)particle->PID, *idx, *idy, *idz); // Cast id to long long for printing PetscInt64
    } else {
        // No valid previous cell ID (e.g., -1,-1,-1 or first time).
        // Start from the first cell in the local owned grid domain.
        *idx = info.xs;
        *idy = info.ys;
        *idz = info.zs;
        LOG_ALLOW(LOCAL,LOG_DEBUG, "InitializeTraversalParameters - Particle %lld has no valid previous cell. Starting search at local corner (%d, %d, %d).\n",
                  (long long)particle->PID, *idx, *idy, *idz);
    }

    // Initialize traversal step counter
    *traversal_steps = 0;

    // Log the chosen starting point
    LOG_ALLOW(LOCAL,LOG_INFO, "InitializeTraversalParameters - Traversal for particle %lld initialized to start at cell (%d, %d, %d).\n",
               (long long)particle->PID, *idx, *idy, *idz);

    return 0;
}

/**
 * @brief Checks if the current CELL indices are within the LOCAL GHOSTED grid boundaries.
 *
 * This function determines if the provided global cell indices fall within the
 * range accessible by the current process, including its ghost cell layers.
 *
 * @param[in]  user       Pointer to the user-defined context containing grid information (needs fda for node info).
 * @param[in]  idx        The global i-index of the current cell.
 * @param[in]  idy        The global j-index of the current cell.
 * @param[in]  idz        The global k-index of the current cell.
 * @param[out] is_within  Pointer to a PetscBool that will be set to PETSC_TRUE if within ghosted bounds, else PETSC_FALSE.
 *
 * @return PetscErrorCode  Returns 0 on success, non-zero on failure.
 */
/*
PetscErrorCode CheckCellWithinLocalGrid(UserCtx *user, PetscInt idx, PetscInt idy, PetscInt idz, PetscBool *is_within)
{
    PetscErrorCode ierr;
    DMDALocalInfo info_nodes;
    PetscInt stencil_width;

    // Validate inputs
    if (user == NULL || is_within == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "Input pointer is NULL.\n");
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input pointer is NULL.");
    }

    // Get node info from fda (to derive cell info)
    ierr = DMDAGetLocalInfo(user->fda, &info_nodes); CHKERRQ(ierr);
    // Get stencil width used for fda (and implicitly for da's ghosts)
    ierr = DMDAGetInfo(user->fda, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &stencil_width, NULL, NULL, NULL, NULL); CHKERRQ(ierr);


    // Get owned cell ranges using the helper
    PetscInt xs_cell, xm_cell, xe_cell;
    PetscInt ys_cell, ym_cell, ye_cell;
    PetscInt zs_cell, zm_cell, ze_cell;
    ierr = GetOwnedCellRange(&info_nodes, 0, &xs_cell, &xm_cell); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_nodes, 1, &ys_cell, &ym_cell); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_nodes, 2, &zs_cell, &zm_cell); CHKERRQ(ierr);
    xe_cell = xs_cell + xm_cell; // Exclusive end of owned cells
    ye_cell = ys_cell + ym_cell;
    ze_cell = zs_cell + zm_cell;

    // Define ghosted cell ranges (inclusive start, exclusive end)
    PetscInt gxs_cell = xs_cell - stencil_width;
    PetscInt gxe_cell = xe_cell + stencil_width;
    PetscInt gys_cell = ys_cell - stencil_width;
    PetscInt gye_cell = ye_cell + stencil_width;
    PetscInt gzs_cell = zs_cell - stencil_width;
    PetscInt gze_cell = ze_cell + stencil_width;

    // Check if cell index (idx, idy, idz) is within the ghosted range
    // Ensure we also check against physical domain bounds (0 to IM-1 etc.)
    // Although typically ghost indices outside might be valid temporarily before migration.
    // Let's keep the check relative to the ghosted region accessible by this process.
    if (idx >= gxs_cell && idx < gxe_cell &&
        idy >= gys_cell && idy < gye_cell &&
        idz >= gzs_cell && idz < gze_cell) {
        *is_within = PETSC_TRUE; // It's within the area this proc can access (owned + ghost)
    }
    else {
        *is_within = PETSC_FALSE; // It's outside the accessible ghosted region
    }

    LOG_ALLOW(LOCAL,LOG_DEBUG, "Cell (%d, %d, %d) is %s the ghosted local grid (x:%d..%d, y:%d..%d, z:%d..%d).\n",
        idx, idy, idz, (*is_within) ? "within" : "outside",
        gxs_cell, gxe_cell, gys_cell, gye_cell, gzs_cell, gze_cell);

    return 0;
}
*/

/**
 * @brief Checks if the current GLOBAL CELL indices are within the LOCAL GHOSTED grid boundaries
 *        accessible by this MPI process.
 *
 * This function determines if the provided global cell indices (idx, idy, idz) fall within the
 * range of cells covered by the current process's owned and ghost NODES.
 * A cell C(i,j,k) (origin N(i,j,k)) is considered within this ghosted region if its origin node N(i,j,k)
 * is within the rank's ghosted nodal region, AND node N(i+1,j+1,k+1) is also within it (to ensure
 * the entire cell extent is covered by available node data). More simply, we check if the cell's
 * origin node is within the range of nodes that can form the start of a ghosted cell.
 *
 * @param[in]  user       Pointer to the user-defined context (needs user->fda for node info).
 * @param[in]  idx        The global i-index of the cell's origin node.
 * @param[in]  idy        The global j-index of the cell's origin node.
 * @param[in]  idz        The global k-index of the cell's origin node.
 * @param[out] is_within  Pointer to a PetscBool that will be set to PETSC_TRUE if within ghosted bounds, else PETSC_FALSE.
 *
 * @return PetscErrorCode  Returns 0 on success, non-zero on failure.
 */
PetscErrorCode CheckCellWithinLocalGrid(UserCtx *user, PetscInt idx, PetscInt idy, PetscInt idz, PetscBool *is_within)
{
    PetscErrorCode ierr;
    DMDALocalInfo info_nodes; // Node information from the DMDA that defines ghost regions (user->fda)

    PetscFunctionBeginUser; // Assuming this is part of your PETSc style

    // Validate inputs
    if (user == NULL || is_within == NULL) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Input pointer is NULL in CheckCellWithinLocalGrid.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input pointer is NULL in CheckCellWithinLocalGrid.");
    }
    if (user->fda == NULL) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "user->fda is NULL in CheckCellWithinLocalGrid. Cannot get ghost info.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "user->fda is NULL. Cannot get ghost info.");
    }

    // Get node info from user->fda (this DMDA has the ghost layer information for nodes)
    ierr = DMDAGetLocalInfo(user->fda, &info_nodes); CHKERRQ(ierr);

    // Determine the range of GLOBAL CELL INDICES that are covered by this rank's ghosted NODAL region.
    // A cell C(i,j,k) has origin node N(i,j,k).
    // The ghosted nodal region starts at global node index info_nodes.gxs and has info_nodes.gxm nodes.

    // Global starting index of the first cell whose origin node is within the ghosted nodal region.
    PetscInt gxs_cell_global_start = info_nodes.gxs;
    PetscInt gys_cell_global_start = info_nodes.gys;
    PetscInt gzs_cell_global_start = info_nodes.gzs;

    // Number of cells that can be formed starting from nodes within the ghosted nodal region.
    // If there are N ghosted nodes (info_nodes.gxm), they can be origins for N-1 cells.
    PetscInt gxm_cell_local_count = (info_nodes.gxm > 0) ? info_nodes.gxm - 1 : 0;
    PetscInt gym_cell_local_count = (info_nodes.gym > 0) ? info_nodes.gym - 1 : 0;
    PetscInt gzm_cell_local_count = (info_nodes.gzm > 0) ? info_nodes.gzm - 1 : 0;

    // Global exclusive end index for cells whose origins are in the ghosted nodal region.
    PetscInt gxe_cell_global_end_exclusive = gxs_cell_global_start + gxm_cell_local_count;
    PetscInt gye_cell_global_end_exclusive = gys_cell_global_start + gym_cell_local_count;
    PetscInt gze_cell_global_end_exclusive = gzs_cell_global_start + gzm_cell_local_count;

    // Check if the given global cell index (idx, idy, idz) falls within this range.
    // This means the origin node of cell (idx,idy,idz) is within the rank's accessible ghosted node region,
    // and that node can indeed serve as a cell origin (i.e., it's not the very last node in the ghosted region).
    if (idx >= gxs_cell_global_start && idx < gxe_cell_global_end_exclusive &&
        idy >= gys_cell_global_start && idy < gye_cell_global_end_exclusive &&
        idz >= gzs_cell_global_start && idz < gze_cell_global_end_exclusive) {
        *is_within = PETSC_TRUE;
    } else {
        *is_within = PETSC_FALSE;
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Cell (origin node glob idx) (%d, %d, %d) is %s the ghosted local grid (covered cell origins x:[%d..%d), y:[%d..%d), z:[%d..%d)).\n",
              idx, idy, idz, (*is_within) ? "within" : "outside",
              gxs_cell_global_start, gxe_cell_global_end_exclusive,
              gys_cell_global_start, gye_cell_global_end_exclusive,
              gzs_cell_global_start, gze_cell_global_end_exclusive);

    PetscFunctionReturn(0);
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
    LOG_ALLOW(LOCAL,LOG_DEBUG, "RetrieveCurrentCell - Cell (%d, %d, %d) vertices \n", idx, idy, idz);
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
    if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)) LOG_FACE_DISTANCES(d); 
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

    LOG_ALLOW(LOCAL, LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Received d: "
	      "d[LEFT=%d]=%.3e, d[RIGHT=%d]=%.3e, d[BOTTOM=%d]=%.3e, "
	      "d[TOP=%d]=%.3e, d[FRONT=%d]=%.3e, d[BACK=%d]=%.3e\n",
	      LEFT, d[LEFT], RIGHT, d[RIGHT], BOTTOM, d[BOTTOM],
	      TOP, d[TOP], FRONT, d[FRONT], BACK, d[BACK]);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Raw d: "
	      "[%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]\n",
	      d[0], d[1], d[2], d[3], d[4], d[5]);
    
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
    if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)) LOG_FACE_DISTANCES(d);

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

    LOG_ALLOW(LOCAL,LOG_DEBUG, "UpdateCellIndicesBasedOnDistances - Updated Indices after clamping (inside domain bounds)  - idx, idy, idz: %d, %d, %d\n", *idx, *idy, *idz);

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
      LOG_ALLOW(LOCAL,LOG_INFO, "FinalizeTraversal - Particle located in cell (%d, %d, %d) after %d traversal steps.\n",
            idx, idy, idz, traversal_steps);
    }
    else {
      LOG_ALLOW(LOCAL,LOG_WARNING, "FinalizeTraversal - Particle could not be located within the grid after %d traversal steps.\n", (PetscInt)traversal_steps);
        particle->cell[0] = -1;
        particle->cell[1] = -1;
        particle->cell[2] = -1;
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "FinalizeTraversal - Completed final traversal sync across all ranks.\n");


    return 0;
}


/**
 * @brief Locates the grid cell containing a particle using a walking search.
 * @ingroup ParticleLocation
 *
 * This function implements a walking search algorithm to find the specific cell
 * (identified by global indices i, j, k) that encloses the particle's physical
 * location (`particle->loc`).
 *
 * The search starts from an initial guess cell (either the particle's previously known
 * cell or the corner of the local process domain) and iteratively steps to adjacent
 * cells based on the signed distances from the particle to the faces of the current
 * cell.
 *
 * Upon successful location (`position == 0`), the particle's `cell` field is updated
 * with the found indices (i,j,k), and the corresponding interpolation `weights`
 * are calculated and stored using the distances (`d`) relative to the final cell.
 *
 * Handles particles exactly on cell boundaries by attempting a tie-breaker if the
 * search gets stuck oscillating between adjacent cells due to the boundary condition.
 *
 * @param[in]  user     Pointer to the UserCtx structure containing grid information (DMDA, coordinates)
 *                      and domain boundaries.
 * @param[in,out] particle Pointer to the Particle structure. Its `loc` field provides the
 *                      target position. On successful return, its `cell` and `weights`
 *                      fields are updated. On failure, `cell` is set to `{-1, -1, -1}`
 *                      and `weights` to `{0.0, 0.0, 0.0}`.
 *
 * @return PetscErrorCode 0 on success. Non-zero error codes may indicate issues during
 *                        coordinate access, distance calculation, or other internal errors.
 *                        A return code of 0 does not guarantee the particle was found;
 *                        check `particle->cell[0] >= 0` afterward.
 *
 * @note Relies on helper functions like `InitializeTraversalParameters`, `CheckCellWithinLocalGrid`,
 *       `RetrieveCurrentCell`, `EvaluateParticlePosition`, `UpdateCellIndicesBasedOnDistances`,
 *       `UpdateParticleWeights`, and `FinalizeTraversal`.
 * @warning Ensure `particle->loc` is set correctly before calling.
 * @warning The function may fail to find the particle if it lies outside the domain accessible
 *          by the current process (including ghost cells) or if `MAX_TRAVERSAL` steps are exceeded.
 */
PetscErrorCode LocateParticleInGrid(UserCtx *user, Particle *particle)
{
    PetscErrorCode ierr;
    PetscInt idx, idy, idz;           // Current search cell indices
    PetscInt traversal_steps;         // Counter for search steps
    PetscBool cell_found = PETSC_FALSE;// Flag indicating success
    const Cmpnts p = particle->loc;   // Particle position (local copy)
    Cell current_cell;                // Geometry of the cell being checked
    const PetscReal threshold = DISTANCE_THRESHOLD; // Base threshold for distance checks
    DMDALocalInfo info;               // Local grid info (for bounds)
    PetscInt repeatedIndexCount = 0;  // Counter for consecutive iterations at the same index
    PetscInt prevIdx = PETSC_MIN_INT; // Previous iteration's index i (for repeat check)
    PetscInt prevIdy = PETSC_MIN_INT; // Previous iteration's index j (for repeat check)
    PetscInt prevIdz = PETSC_MIN_INT; // Previous iteration's index k (for repeat check)
    PetscInt last_position = -999;    // Position result (-1, 0, >=1) from *previous* iteration

    PetscFunctionBeginUser;
    LOG_FUNC_TIMER_BEGIN_EVENT(EVENT_Individualwalkingsearch, LOCAL);

    // Get local grid information (needed for bounds and index updates)
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

    // Determine starting cell indices (previous cell or local corner)
    ierr = InitializeTraversalParameters(user, particle, &idx, &idy, &idz, &traversal_steps); CHKERRQ(ierr);

    // --- Main Walking Search Loop ---
    while (!cell_found && traversal_steps < MAX_TRAVERSAL) {
        traversal_steps++;

        // Log current state for debugging
        LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, traversal_steps, 1,
                       "LocateParticleInGrid [PID %lld, Step %d]: Checking cell (%d, %d, %d). Pos (%.6f, %.6f, %.6f)\n",
                       (long long)particle->PID, traversal_steps, idx, idy, idz, p.x, p.y, p.z);

        // --- Check for Stuck Loop / Oscillation ---
        if (idx == prevIdx && idy == prevIdy && idz == prevIdz) {
            repeatedIndexCount++;
            LOG_ALLOW(LOCAL, LOG_DEBUG,"LocateParticleInGrid [PID %lld, Step %d]: Repeated index (%d,%d,%d) count = %d\n",
                        (long long)particle->PID, traversal_steps, idx, idy, idz, repeatedIndexCount);

            if (repeatedIndexCount > REPEAT_COUNT_THRESHOLD) {
                // --- Option B Tie-Breaker Logic ---
                if (last_position >= 1) { // Was the last step result 'on boundary'? Likely oscillation.
                    LOG_ALLOW(LOCAL, LOG_WARNING,
                              "LocateParticleInGrid [PID %lld]: Stuck at index (%d,%d,%d) after being on boundary for >%d steps. Applying boundary tie-break.\n",
                              (long long)particle->PID, idx, idy, idz, REPEAT_COUNT_THRESHOLD);

                    // Attempt to assign to the cell where we are stuck.
                    // Need to ensure the cell is valid and re-get distances for weights.
                    Cell final_cell;
                    PetscReal final_d[NUM_FACES];
                    PetscInt final_position; // Result not used, just need distances
                    PetscBool is_stuck_cell_valid;

                    ierr = CheckCellWithinLocalGrid(user, idx, idy, idz, &is_stuck_cell_valid); CHKERRQ(ierr);
                    if (is_stuck_cell_valid) {
                         ierr = RetrieveCurrentCell(user, idx, idy, idz, &final_cell); CHKERRQ(ierr);
                         // Re-evaluate to get distances 'final_d' relative to this specific cell
                         ierr = EvaluateParticlePosition(&final_cell, final_d, p, &final_position, threshold);
                         if (ierr == 0) { // Check if evaluation succeeded
                            ierr = UpdateParticleWeights(final_d, particle); CHKERRQ(ierr); // Calculate weights
                            particle->cell[0] = idx; // Assign the stuck cell index
                            particle->cell[1] = idy;
                            particle->cell[2] = idz;
                            cell_found = PETSC_TRUE; // Mark as found via tie-break
                            LOG_ALLOW(LOCAL, LOG_INFO,
                                      "LocateParticleInGrid [PID %lld]: Assigned via tie-break to cell (%d,%d,%d).\n",
                                      (long long)particle->PID, idx, idy, idz);
                         } else {
                             // Error during final evaluation
                             LOG_ALLOW(LOCAL, LOG_ERROR,
                                       "LocateParticleInGrid [PID %lld]: Error %d during final evaluation for tie-break at cell (%d,%d,%d). Search fails.\n",
                                       (long long)particle->PID, ierr, idx, idy, idz);
                             cell_found = PETSC_FALSE; // Ensure failure
                         }
                    } else {
                         LOG_ALLOW(LOCAL, LOG_WARNING,
                                   "LocateParticleInGrid [PID %lld]: Stuck at invalid cell index (%d,%d,%d) during tie-break attempt. Search fails.\n",
                                   (long long)particle->PID, idx, idy, idz);
                         cell_found = PETSC_FALSE; // Mark as failed
                    }
                    break; // Exit loop after attempting tie-break
                } else {
                    // Stuck for > threshold steps, but last known state wasn't 'on boundary'.
                    LOG_ALLOW(LOCAL, LOG_WARNING,
                              "LocateParticleInGrid [PID %lld]: Stuck at index (%d,%d,%d) for >%d steps (last_pos=%d). Breaking search (failed).\n",
                              (long long)particle->PID, idx, idy, idz, REPEAT_COUNT_THRESHOLD, last_position);
                    cell_found = PETSC_FALSE; // Mark as failure
                    break; // Exit loop (failed)
                }
            } // end if (repeatedIndexCount > threshold)
        } else { // Index changed from previous iteration
            repeatedIndexCount = 0;
        }
        // --- End Stuck Loop Check ---

        // --- Update previous index state for the *next* iteration's check ---
        prevIdx = idx;
        prevIdy = idy;
        prevIdz = idz;

        // --- Check if current cell index is within accessible grid bounds ---
        PetscBool is_within;
        ierr = CheckCellWithinLocalGrid(user, idx, idy, idz, &is_within); CHKERRQ(ierr);
        if (!is_within) {
            LOG_ALLOW(LOCAL, LOG_WARNING,
                      "LocateParticleInGrid [PID %lld]: Search moved outside local ghosted grid boundaries at cell (%d, %d, %d). Breaking search.\n",
                      (long long)particle->PID, idx, idy, idz);
            cell_found = PETSC_FALSE; // Mark as failed
            break; // Exit loop
        }

        // --- Retrieve geometry and evaluate particle position ---
        ierr = RetrieveCurrentCell(user, idx, idy, idz, &current_cell); CHKERRQ(ierr);
        PetscReal current_d[NUM_FACES]; // Holds distances relative to current_cell
        PetscInt position;              // Holds result: -1 (outside), 0 (inside), >=1 (boundary)
        ierr = EvaluateParticlePosition(&current_cell, current_d, p, &position, threshold); CHKERRQ(ierr);

        // Store the position result for the *next* iteration's stuck check
        last_position = position;

        LOG_ALLOW(LOCAL, LOG_DEBUG,
                   "LocateParticleInGrid [PID %lld, Step %d]: Evaluated pos in cell (%d,%d,%d): position=%d\n",
                   (long long)particle->PID, traversal_steps, idx, idy, idz, position);

        // --- Main Logic: Handle position result ---
        if (position == 0) { // Strictly INSIDE the cell
            cell_found = PETSC_TRUE;
            particle->cell[0] = idx; // Assign the indices of the cell we just evaluated
            particle->cell[1] = idy;
            particle->cell[2] = idz;
            // Calculate weights using the distances 'current_d' from this successful evaluation
            ierr = UpdateParticleWeights(current_d, particle); CHKERRQ(ierr);
            LOG_ALLOW(LOCAL, LOG_INFO,
                      "LocateParticleInGrid [PID %lld]: Found INSIDE cell (%d, %d, %d) after %d steps. Weights calculated.\n",
                      (long long)particle->PID, idx, idy, idz, traversal_steps);
            break; // Found the correct cell, exit loop

        } else { // Particle is OUTSIDE (position == -1) or ON BOUNDARY (position >= 1)
            LOG_ALLOW(LOCAL, LOG_DEBUG,
                      "LocateParticleInGrid [PID %lld, Step %d]: Particle OUTSIDE or ON BOUNDARY of cell (%d, %d, %d), position = %d. Updating indices.\n",
                      (long long)particle->PID, traversal_steps, idx, idy, idz, position);
            // Both cases require stepping to the next potential cell.
            // Update search indices (idx, idy, idz) based on the distances 'current_d'
            ierr = UpdateCellIndicesBasedOnDistances(current_d, &idx, &idy, &idz, &info); CHKERRQ(ierr);
            // Loop continues to the next iteration with potentially updated idx, idy, idz
        }
    } // --- End while loop ---

    // --- Handle cases where the loop terminated WITHOUT finding the cell ---
    // (Could be MAX_TRAVERSAL, out of bounds, or non-boundary stuck state)
    if (!cell_found) {
        LOG_ALLOW(LOCAL, LOG_WARNING,
                  "LocateParticleInGrid [PID %lld]: Search FAILED after %d steps. Final check was at cell (%d,%d,%d). Resetting cell/weights.\n",
                  (long long)particle->PID, traversal_steps, idx, idy, idz);
        // Ensure particle state reflects failure explicitly
        particle->cell[0] = -1;
        particle->cell[1] = -1;
        particle->cell[2] = -1;
    }

    // Finalize traversal by reporting the results (logs success/failure message)
    // Note: idx, idy, idz passed here are the *last checked* indices for logging purposes.
    // The actual found cell index is stored in particle->cell if cell_found is true.
    ierr = FinalizeTraversal(user, particle, traversal_steps, cell_found, idx, idy, idz); CHKERRQ(ierr);

    LOG_FUNC_TIMER_END_EVENT(EVENT_Individualwalkingsearch, LOCAL);
    PetscFunctionReturn(0);
}

/**
 * @brief Finds the MPI rank that owns a given global cell index.
 * @ingroup DomainInfo
 *
 * This function performs a linear search through the pre-computed decomposition map
 * (`user->RankCellInfoMap`) to determine which process is responsible for the cell
 * with global indices (i, j, k). It is the definitive method for resolving cell
 * ownership in the "Walk and Handoff" migration algorithm.
 *
 * If the provided indices are outside the range of any rank (e.g., negative or
 * beyond the global domain), the function will not find an owner and `owner_rank`
 * will be set to -1.
 *
 * @param[in]  user       Pointer to the UserCtx structure, which must contain the
 *                        initialized `RankCellInfoMap` and `num_ranks`.
 * @param[in]  i          Global i-index of the cell to find.
 * @param[in]  j          Global j-index of the cell to find.
 * @param[in]  k          Global k-index of the cell to find.
 * @param[out] owner_rank Pointer to a `PetscMPIInt` where the resulting owner rank will
 *                        be stored. It is set to -1 if no owner is found.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
PetscErrorCode FindOwnerOfCell(UserCtx *user, PetscInt i, PetscInt j, PetscInt k, PetscMPIInt *owner_rank)
{
    PetscErrorCode ierr;
    PetscMPIInt size;

    PetscFunctionBeginUser;

    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

    // --- 1. Input Validation ---
    if (!user || !user->RankCellInfoMap) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "UserCtx or RankCellInfoMap is not initialized in FindOwnerOfCell.");
    }
    if (!owner_rank) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output pointer owner_rank is NULL in FindOwnerOfCell.");
    }

    // --- 2. Linear Search through the Decomposition Map ---
    // Initialize to a "not found" state.
    *owner_rank = -1;

    // Loop through the map, which contains the ownership info for every rank 'r'.
    for (PetscMPIInt r = 0; r < size; ++r) {
        const RankCellInfo *info = &user->RankCellInfoMap[r];

        // A rank owns a cell if the cell's index is within its start (inclusive)
        // and end (exclusive) range for all three dimensions.
        if ((i >= info->xs_cell && i < info->xs_cell + info->xm_cell) &&
            (j >= info->ys_cell && j < info->ys_cell + info->ym_cell) &&
            (k >= info->zs_cell && k < info->zs_cell + info->zm_cell))
        {
            *owner_rank = r; // We found the owner.
            break;           // The search is over, exit the loop.
        }
    }

    // --- 3. Logging for Diagnostics ---
    // This is extremely useful for debugging particle migration issues.
    if (*owner_rank == -1) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "FindOwnerOfCell: No owner found for global cell (%d, %d, %d). It is likely outside the domain.\n", i, j, k);
    } else {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "FindOwnerOfCell: Owner of cell (%d, %d, %d) is Rank %d.\n", i, j, k, *owner_rank);
    }

    PetscFunctionReturn(0);
}
