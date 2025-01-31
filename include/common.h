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
typedef struct {
    // Grid-related fields
    DM da;                  ///< Data structure for scalar fields.
    DM fda;                 ///< Data structure for vector fields.
    DM fda2;                  ///< Data structure for RANS fields.
    PetscInt IM, JM, KM;    ///< Global grid dimensions in x, y, z directions.
    BoundingBox bbox;       ///< Bounding box for the local grid domain.
    DMDALocalInfo info;     ///< Local information about the DMDA.
    PetscReal rx;           // Stretching ratio in x-direction
    PetscReal ry;           // Stretching ratio in y-direction
    PetscReal rz;           // Stretching ratio in z-direction

    // Simulation fields
    Vec Ucont;              ///< Contravariant velocity field.
    Vec Ucat;               ///< Cartesian velocity field.
    Vec P;                  ///< Pressure field.
    Vec Nvert;              ///< Node state field (fluid, solid, etc.).
    Vec Nvert_o;            ///< Node state field in the previous timestep.
    Vec lUcont, lNvert;     ///< Local versions of Ucont and Nvert.

    // Statistical fields
    Vec Ucat_sum;           ///< Sum of Cartesian velocity for averaging.
    Vec Ucat_cross_sum;     ///< Cross-product sum of Cartesian velocities.
    Vec Ucat_square_sum;    ///< Squared velocity sum for RMS calculations.
    Vec P_sum;              ///< Sum of pressure values.

    // LES-specific fields
    Vec lCs;                ///< Eddy viscosity field for LES.
    Vec Cs;                 ///< Global version of the LES constant field.

    // RANS-specific fields
    Vec K_Omega;            ///< Turbulent kinetic energy (K) and specific dissipation rate (Omega).
    Vec K_Omega_o;          ///< Old values of K_Omega for time-stepping.
  Vec lK_Omega,lK_Omega_o;  ///< Local K_Omega and old K_Omega for each process.

    // Particle-related fields
    DM swarm;               ///< Particle data structure using DMSwarm.
    PetscInt *miglist;      ///< List of ranks for particle migration.
    PetscInt ParticleInitialization;

    // Simulation parameters
    PetscReal dt;           ///< Time step.
    PetscReal ren;          ///< Reynolds number.

    // Flags for simulation modes
    PetscBool averaging;    ///< Flag to indicate whether statistical averaging is enabled.
    PetscBool les;          ///< Flag to indicate if LES is active.
    PetscBool rans;         ///< Flag to indicate if RANS is active.

    // Miscellaneous fields
    PetscInt _this;         ///< Current block index.
    PetscInt ti;            ///< Current timestep index.
} UserCtx;

/* --------------------------------------------------------------------
   VTKFileType

   This enum indicates whether we should write a .vts (StructuredGrid)
   or a .vtp (PolyData) output file when we finalize the VTK data.

   - VTK_STRUCTURED => .vts
   - VTK_POLYDATA   => .vtp

   Example usage:
     VTKFileType fileType = VTK_STRUCTURED; // for structured data
     VTKFileType fileType = VTK_POLYDATA;   // for particle data
-------------------------------------------------------------------- */
typedef struct {
    PetscInt  startTime;       /* start of time loop */
    PetscInt  endTime;         /* end of time loop */
    PetscInt  timeStep;        /* increment of time loop */

    char      eulerianExt[8];  /* "vts", "vtk", etc. for grid data */
    char      particleExt[8];  /* "vtp", "vtk", etc. for particle data */

    PetscBool outputParticles; /* whether to write particle data or not */
} PostProcessParams;

typedef enum {
    VTK_STRUCTURED,
    VTK_POLYDATA
} VTKFileType;

/* You can define a max # of fields, or make this dynamic. */
#define MAX_POINT_DATA_FIELDS 10

/**
 * @brief Metadata structure for VTK output (.vts or .vtp files).
 *
 * This struct is used for both Eulerian grid data and Lagrangian (particle) data.
 * It stores coordinate arrays, scalar/vector fields, and connectivity details
 * required to write valid VTK files.
 */
typedef struct _n_VTKMetaData {
    /* 1) File type: VTK_STRUCTURED (Eulerian grid) or VTK_POLYDATA (Particle data) */
    VTKFileType  fileType;
 
    /* 2) For VTK_STRUCTURED => .vts */
    int          mx, my, mz;   /* Grid dimensions */
    int          nnodes;       /* Total number of structured grid nodes */

    /* 3) For VTK_POLYDATA => .vtp */
    int          npoints;      /* Number of particle points */

    /* 4) Coordinate and field data (shared) */
    double      *coords;       /* Coordinates array (3 components per point/node) */

    /* Scalar field metadata */
    double      *scalarField;
    const char  *scalarFieldName;
    int          numScalarFields;

    /* Vector field metadata */
    double      *vectorField;       /* Interleaved x,y,z vector field */
    const char  *vectorFieldName;   /* Name of the vector field (e.g., "Velocity") */
    int          numVectorFields;   /* Number of vector fields (0 or 1) */

    /* 5) For polydata connectivity */
    int         *connectivity; 
    int         *offsets;
} VTKMetaData;


// -------- MultiField Particle vtp file implementation ------------

/*
// A small struct to hold one field's metadata 
typedef struct {
    char    name[64];       // e.g. "velocity", "pressure", "temp", etc.
    int     numComponents;  // e.g. 1 = scalar, 3 = vector, etc. 
    double *data;           // pointer to the actual data array 
} PointDataField;

// --------------------------------------------------------------------
//   VTKMetaData

//   Works for both .vts (structured) and .vtp (polydata).
//   Now uses an array of 'pointDataFields' for all fields in <PointData>.
//-------------------------------------------------------------------- 
typedef struct _n_VTKMetaData {
    // 1) File type: VTK_STRUCTURED or VTK_POLYDATA 
    VTKFileType  fileType;
 
    // 2) For VTK_STRUCTURED => .vts 
    int          mx, my, mz;   // domain dimensions 
    int          nnodes;       // e.g., (mx-1)*(my-1)*(mz-1) or actual node count 

    // 3) For VTK_POLYDATA => .vtp 
    int          npoints;     //  number of point "particles" 

 // 4) Shared coordinate array 
    double      *coords;       // (3 * npoints) or (3 * nnodes) if you store them explicitly 

    // 5) Array-based approach for multiple fields //
    int           numFields;   // how many fields in <PointData>? 
    PointDataField pointDataFields[MAX_POINT_DATA_FIELDS];

    // 6) For polydata connectivity (only used if fileType == VTK_POLYDATA) 
    int         *connectivity; 
    int         *offsets;
} VTKMetaData;
*/

// Add other shared types (e.g., IBMNodes, IBMInfo) if needed

#endif // COMMON_H
