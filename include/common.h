/**
 * @file common.h
 * @brief Foundational type definitions and structures for the entire simulation.
 *
 * This header is the central location for all shared data structures and enumerations.
 * It is structured to avoid circular dependencies by using forward declarations for
 * complex types before their full definition. It should be included by most other
 * files in the project to ensure access to a consistent set of types like UserCtx,
 * Cmpnts, and the boundary condition system structs.
 */

#ifndef COMMON_H
#define COMMON_H

// --- Primary Library Includes ---
#include <petsc.h>
#include <petscdmda.h>
#include <petscdmswarm.h>

// --- Standard C Library Includes ---
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h> 

//================================================================================
//
//                 1. FORWARD DECLARATIONS & BASIC TYPES
//
//================================================================================

// --- Forward Declarations ---
// These declarations inform the compiler that these struct types exist, allowing
// pointers to them to be used before their full content is defined. This is
// essential for breaking circular dependencies between complex structs.
typedef struct BC_Param_s BC_Param;
typedef struct BoundaryCondition_s BoundaryCondition;
typedef struct BoundaryFaceConfig_s BoundaryFaceConfig;
typedef struct UserCtx_s UserCtx;

// --- Foundational Geometric and Data Types ---

/** @brief A 3D point or vector with PetscScalar components. */
typedef struct {
    PetscScalar x, y, z;
} Cmpnts;

/** @brief A 2D point or vector with PetscScalar components. */
typedef struct {
    PetscScalar x, y;
} Cmpnts2;

/** @brief Represents a single point in a time-varying flow waveform. */
typedef struct {
    PetscReal t, f;
} FlowWave;

/** @brief Defines a 3D axis-aligned bounding box. */
typedef struct {
    Cmpnts min_coords; ///< Minimum x, y, z coordinates of the bounding box.
    Cmpnts max_coords; ///< Maximum x, y, z coordinates of the bounding box.
} BoundingBox;

/** @brief Defines the vertices of a single hexahedral grid cell. */
typedef struct {
    Cmpnts vertices[8]; ///< Coordinates of the eight vertices of the cell.
} Cell;

/** @brief Defines a particle's core properties for Lagrangian tracking. */
typedef struct {
    PetscInt64 PID;     ///< Unique Particle ID.
    PetscInt cell[3];   ///< Computational indices (i, j, k) of the cell containing the particle.
    Cmpnts loc;         ///< Physical location (x,y,z) of the particle.
    Cmpnts vel;         ///< Physical velocity (vx,vy,vz) of the particle.
    Cmpnts weights;     ///< Interpolation weights within its host cell.
} Particle;

/** @brief Stores the MPI ranks of neighboring subdomains. */
typedef struct {
    PetscMPIInt rank_xm, rank_xp; // Neighbors at -x, +x
    PetscMPIInt rank_ym, rank_yp; // Neighbors at -y, +y
    PetscMPIInt rank_zm, rank_zp; // Neighbors at -z, +z
} RankNeighbors;


//================================================================================
//
//                 2. BOUNDARY CONDITION SYSTEM ENUMS
//
//================================================================================

/** @brief Identifies the six logical faces of a structured computational block. */
typedef enum {
    BC_FACE_NEG_X = 0, BC_FACE_POS_X = 1,
    BC_FACE_NEG_Y = 2, BC_FACE_POS_Y = 3,
    BC_FACE_NEG_Z = 4, BC_FACE_POS_Z = 5
} BCFace;

/** @brief Defines the general mathematical/physical category of a boundary. */
typedef enum {
    UNDEFINED = 0, NOGRAD, WALL, SYMMETRY, INLET,
    OUTLET, FARFIELD, PERIODIC, INTERFACE
} BCType;

/** @brief Defines the specific computational "strategy" for a boundary handler. */
typedef enum {
    BC_HANDLER_UNDEFINED = 0,BC_HANDLER_NOGRAD_COPY_GHOST,
    BC_HANDLER_WALL_NOSLIP, BC_HANDLER_WALL_MOVING,
    BC_HANDLER_SYMMETRY_PLANE,
    BC_HANDLER_INLET_CONSTANT_VELOCITY, BC_HANDLER_INLET_PULSANTILE_FLUX, BC_HANDLER_INLET_DEVELOPED_PROFILE,
    BC_HANDLER_OUTLET_CONSERVATION, BC_HANDLER_OUTLET_PRESSURE,
    BC_HANDLER_FARFIELD_NONREFLECTING,
    BC_HANDLER_PERIODIC, BC_HANDLER_INTERFACE_OVERSET
} BCHandlerType;


//================================================================================
//
//               3. BOUNDARY CONDITION SYSTEM STRUCTS
//
//================================================================================

/** @brief A node in a linked list for storing key-value parameters from the bcs.dat file. */
struct BC_Param_s {
    char *key;
    char *value;
    struct BC_Param_s *next;
};

/** @brief Provides execution context for a boundary condition handler. */
typedef struct {
    UserCtx  *user;             ///< Access to all global simulation data.
    BCFace    face_id;          ///< The geometric face (0-5) being processed.
    const PetscReal *global_inflow_sum;  ///< Pointer to total domain inflow for the current step.
    const PetscReal *global_outflow_sum; ///< Pointer to total measured domain outflow for the current step.
} BCContext;

/** @brief The "virtual table" struct for a boundary condition handler object. */
struct BoundaryCondition_s {
    BCHandlerType type;
    void         *data;
    PetscErrorCode (*Initialize)(BoundaryCondition *self, BCContext *ctx);
    PetscErrorCode (*PreStep)(BoundaryCondition *self, BCContext *ctx, PetscReal *local_inflow, PetscReal *local_outflow);
    PetscErrorCode (*Apply)(BoundaryCondition *self, BCContext *ctx);
    PetscErrorCode (*PlaceSource)(BoundaryCondition *self, BCContext *ctx, ...);
    PetscErrorCode (*Destroy)(BoundaryCondition *self);
};

/** @brief Holds the complete configuration for one of the six boundary faces. */
struct BoundaryFaceConfig_s {
    BCFace             face_id;
    BCType             mathematical_type;
    BCHandlerType      handler_type;
    BC_Param           *params;
    BoundaryCondition *handler;
};


//================================================================================
//
//                4. MAIN USER CONTEXT (GOD) STRUCT
//
//================================================================================

/**
 * @brief User-defined context containing all simulation data and configurations.
 */
struct UserCtx_s {
    // --- Grid & Parallelization ---
    DM da;                      ///< DMDA for scalar fields (P, Nvert).
    DM fda;                     ///< DMDA for primary vector fields (Ucat, Ucont, Metrics).
    DM fda2;                    ///< DMDA for secondary vector fields (e.g., RANS).
    DMDALocalInfo info;         ///< Cached local grid info for the current rank.
    PetscInt IM, JM, KM;        ///< Global grid dimensions (number of cells in i,j,k).
    PetscReal xMin,yMin,zMin;   ///< Physical minimum bounds of the grid.
    PetscReal xMax,yMax,zMax;   ///< Physical maximum bounds of the grid.
    PetscReal rx, ry, rz;       ///< Grid stretching ratios.
    PetscInt nblk;              ///< Number of grid blocks in the simulation.
    BoundingBox bbox;           ///< Bounding box for the local processor's grid domain.
    BoundingBox global_domain_bbox; ///< Bounding box for the entire global domain.
    RankNeighbors neighbors;    ///< MPI ranks of neighboring subdomains.
    PetscInt GridOrientation; ///< Integerer that determines Whether grid is right-handed (+1) or left-handed (-1) i.e normals are facing outward or inward.
    // --- Boundary Condition System ---
    BoundaryFaceConfig boundary_faces[6]; ///< The new, primary BC configuration array.
    BCType    face_bc_types[6];
    PetscBool inletFaceDefined;           ///< Legacy flag for particle system compatibility.
    BCFace    identifiedInletBCFace;      ///< Legacy field for particle system compatibility.

    // --- Primary Simulation Fields ---
    Vec Ucont, lUcont;          ///< Global and local Contravariant velocity.
    Vec Ucat, lUcat;            ///< Global and local Cartesian velocity.
    Vec P, lP;                  ///< Global and local Pressure.
    Vec Nvert, lNvert;          ///< Global and local Node state field (fluid, solid, etc.).
    Vec Nvert_o, lNvert_o;      ///< Global and local Node state from previous timestep.

    // --- Curvilinear Grid Metrics ---
    Vec Csi, Eta, Zet;          ///< (DEPRECATED: Use local) Global metric vectors.
    Vec ICsi, IEta, IZet;       ///< (DEPRECATED)
    Vec JCsi, JEta, JZet;       ///< (DEPRECATED)
    Vec KCsi, KEta, KZet;       ///< (DEPRECATED)
    Vec lCsi, lEta, lZet;       ///< Local metric vectors (primary).
    Vec lICsi, lIEta, lIZet;    ///< Local inverse metric vectors.
    Vec lJCsi, lJEta, lJZet;    ///< Local inverse metric vectors.
    Vec lKCsi, lKEta, lKZet;    ///< Local inverse metric vectors.
    Vec Aj, IAj, JAj, KAj;      ///< (DEPRECATED) Global Jacobian vectors.
    Vec lAj, lIAj, lJAj, lKAj;  ///< Local Jacobian vectors.

    // --- Particle System ---
    DM swarm;                   ///< DMSwarm object for particle data.
    PetscMPIInt *miglist;      ///< List of ranks for particle migration.
    PetscInt NumberofParticles; ///< Total number of particles in the simulation.
    Vec ParticleCount;          ///< Eulerian field to count particles per cell.
    PetscInt ParticleInitialization; ///< Flag controlling how particles are initially placed.

    // --- Simulation Parameters & State ---
    PetscReal dt;               ///< Time step size.
    PetscReal ren;              ///< Reynolds number.
    PetscReal ti;               ///< Current simulation time.
    PetscInt step;              ///< Current timestep index.
    PetscInt FieldInitialization; ///< Flag controlling how fields are initially set.
    PetscInt LoggingFrequency;  ///< Frequency for detailed logging.
    
    // --- Initial Condition Parameters ---.
    Cmpnts    InitialConstantContra;   ///< A constant contravariant velocity vector.
    PetscReal InitialConstantPressure; ///< A constant pressure value.
    PetscReal InitialConstantNvert;    ///< A constant nvert value.
    
    // --- Simulation Mode Flags ---
    PetscBool averaging;        ///< Flag to enable statistical averaging.
    PetscBool les;              ///< Flag to enable Large Eddy Simulation.
    PetscBool rans;             ///< Flag to enable Reynolds-Averaged Navier-Stokes.

    // --- Statistical Averaging Fields ---
    Vec Ucat_sum;
    Vec Ucat_cross_sum;
    Vec Ucat_square_sum;
    Vec P_sum;

    // --- Turbulence Model Fields (LES/RANS) ---
    Vec lCs, Cs;                // LES specific.
    Vec K_Omega, K_Omega_o;     // RANS specific.
    Vec lK_Omega, lK_Omega_o;   // RANS specific.
    
    // --- Miscellaneous ---
    PetscInt _this;             ///< Legacy block index, likely for multi-block contexts.
};

//================================================================================
//
//                     5. OTHER MISC. DEFINITIONS
//
//================================================================================
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
    char      eulerianPrefix[20]; /* directory to save euler fields in */
    char      particlePrefix[20]; /* directory to save particle fields in */
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
    PetscInt          mx, my, mz;   /* Grid dimensions */
    PetscInt          nnodes;       /* Total number of structured grid nodes */

    /* 3) For VTK_POLYDATA => .vtp */
    PetscInt          npoints;      /* Number of particle points */

    /* 4) Coordinate and field data (shared) */
    PetscScalar      *coords;       /* Coordinates array (3 components per point/node) */

    /* Scalar field metadata */
    PetscScalar      *scalarField;
    const char  *scalarFieldName;
    PetscInt          numScalarFields;

    /* Vector field metadata */
    PetscScalar      *vectorField;       /* Interleaved x,y,z vector field */
    const char  *vectorFieldName;   /* Name of the vector field (e.g., "Velocity") */
    PetscInt          numVectorFields;   /* Number of vector fields (0 or 1) */

    /* 5) For polydata connectivity */
    PetscInt         *connectivity; 
    PetscInt         *offsets;
} VTKMetaData;

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

/** @brief Information needed to migrate a single particle between MPI ranks. */
typedef struct {
    PetscInt local_index;
    PetscInt target_rank;
} MigrationInfo;

#endif // COMMON_H
