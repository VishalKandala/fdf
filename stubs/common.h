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
#include <ctype.h> 

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
    PetscInt cell[3]; /**< Indices of the cell containing the particle (i, j, k). */
    Cmpnts loc;         /**< Location of the particle in 3D space. */
    Cmpnts vel;         /**< Velocity of the particle in 3D space. */
    Cmpnts weights;     /**< Weights associated with the particle (e.g., for interpolation). */
  //PetscReal P;          // Pressure associated with the particle //
} Particle;

// Structure to hold neighbor ranks (using common directions)
typedef struct {
    PetscMPIInt rank_xm; // Neighbor at x-
    PetscMPIInt rank_xp; // Neighbor at x+
    PetscMPIInt rank_ym; // Neighbor at y-
    PetscMPIInt rank_yp; // Neighbor at y+
    PetscMPIInt rank_zm; // Neighbor at z-
    PetscMPIInt rank_zp; // Neighbor at z+
    // Add edge/corner neighbors if needed later
} RankNeighbors;

/**
 * @brief Represents a cell in the computational grid with its eight vertices.
 */
typedef struct {
    Cmpnts vertices[8]; /**< Coordinates of the eight vertices of the cell. */
} Cell;

//////////// BC Structs //////////////////////////////

/* Domain face identifiers */
typedef enum {
    BC_FACE_NEG_X = 0,
    BC_FACE_POS_X = 1,
    BC_FACE_NEG_Y = 2,
    BC_FACE_POS_Y = 3,
    BC_FACE_NEG_Z = 4,
    BC_FACE_POS_Z = 5
} BCFace;

/* BC type identifiers */
typedef enum {
    UNDEFINED = 0,
    WALL,
    SYMMETRY,
    INLET,
    OUTLET,
    FARFIELD,
    PERIODIC,
    INTERFACE
    /* Extend as needed (e.g., Robin conditions) */
} BCType;

/**
 * @brief Defines the specific computational "strategy" or implementation for a boundary.
 *
 * This enum allows for multiple behaviors for the same general BCType. For example,
 * a BC_TYPE_INLET can be implemented by a constant velocity handler or a pulsatile one.
 */
typedef enum {
    BC_HANDLER_UNDEFINED = 0,

    // --- Wall & Symmetry Handlers ---
    BC_HANDLER_WALL_NOSLIP,             ///< Behavior: Stationary, no-slip (zero velocity) wall.
    BC_HANDLER_WALL_MOVING,             ///< Behavior: Wall with a prescribed, non-zero velocity.
    BC_HANDLER_SYMMETRY_PLANE,          ///< Behavior: Zero normal velocity, zero gradient for tangential components.

    // --- Inlet Handlers ---
    BC_HANDLER_INLET_CONSTANT_VELOCITY, ///< Behavior: Prescribed, uniform velocity vector.
    BC_HANDLER_INLET_PULSANTILE_FLUX,   ///< Behavior: Time-varying flux profile read from a file.
    BC_HANDLER_INLET_DEVELOPED_PROFILE, ///< Behavior: Spatially-varying velocity (e.g., parabolic).

    // --- Outlet Handlers ---
    BC_HANDLER_OUTLET_CONSERVATION,     ///< Behavior: Extrapolated velocity, corrected for global mass conservation.
    BC_HANDLER_OUTLET_PRESSURE,         ///< Behavior: Prescribes a static pressure, allowing velocity to develop.

    // --- Other Physical Handlers ---
    BC_HANDLER_FARFIELD_NONREFLECTING,  ///< Behavior: Characteristic-based non-reflecting condition.

    // --- Multi-Block / Interface Handlers ---
    BC_HANDLER_PERIODIC,                ///< Behavior: Copies ghost cell data from the opposite domain face.
    BC_HANDLER_INTERFACE_OVERSET        ///< Behavior: Interpolates data from an overlapping grid block.

} BCHandlerType;

/**
 * @brief Provides the necessary context for a boundary condition handler to perform its work.
 *
 * This struct is a clean, temporary container populated by the master controller
 * each time step and passed to the handler functions. It avoids the use of global
 * state variables for time-dependent data.
 */
typedef struct {
    UserCtx  *user;                       ///< Pointer to all global user data (DMs, Vecs, parameters).
    BCFace    face_id;                    ///< The geometric face (0-5) currently being processed.

    // --- Time-Step Scoped Data ---
    // These pointers are valid during the Apply phase, pointing to the results
    // of the PreStep and parallel reduction phases.
    const PetscReal *global_inflow_sum;   ///< Pointer to the total domain inflow for this step.
    const PetscReal *global_outflow_sum;  ///< Pointer to the total measured domain outflow for this step.

} BCContext;

/**
 * @brief Represents a "live" boundary condition handler object.
 *
 * This struct encapsulates the specific behavior and data for a boundary face
 * through a standardized set of function pointers (a "virtual table"). An instance
 * of this struct is a concrete implementation of a boundary condition strategy.
 */
struct BoundaryCondition_s {
    BCHandlerType type;   ///< The specific handler type of this object.
    void         *data;   ///< Pointer to a private, handler-specific data struct (e.g., for waveform data). Must be cast to the correct type inside handler functions.

    /**
     * @brief (Phase 1) Called once at setup to initialize the handler.
     *
     * This function is responsible for allocating and populating the private `data`
     * struct, reading necessary parameters from options/files, and setting the
     * initial field values (e.g., Ucat, Ucont) on its assigned boundary face.
     * @param self Pointer to this BoundaryCondition object.
     * @param ctx  The context for initialization.
     * @return 0 on success.
     */
    PetscErrorCode (*Initialize)(BoundaryCondition *self, BCContext *ctx);

    /**
     * @brief (Phase 2) Called every time step to measure local flux contributions.
     *
     * This function calculates the inflow or outflow flux on the portion of the
     * boundary face owned by the current MPI rank.
     * @param self Pointer to this BoundaryCondition object.
     * @param ctx  The context for the current step.
     * @param[out] local_inflow_contribution  The calculated inflow flux on this rank.
     * @param[out] local_outflow_contribution The calculated outflow flux on this rank.
     * @return 0 on success.
     */
    PetscErrorCode (*PreStep)(BoundaryCondition *self, BCContext *ctx, PetscReal *local_inflow_contribution, PetscReal *local_outflow_contribution);

    /**
     * @brief (Phase 3) Called every time step to apply the boundary condition.
     *
     * This function updates the contravariant fluxes (ucont) and ghost-cell
     * Cartesian velocities (ucat) on its assigned face. For corrective BCs like
     * outlets, it will use the global flux data available in the `ctx`.
     * @param self Pointer to this BoundaryCondition object.
     * @param ctx  The context for the current step, now including global flux data.
     * @return 0 on success.
     */
    PetscErrorCode (*Apply)(BoundaryCondition *self, BCContext *ctx);

    /**
     * @brief (Optional) Called to place a particle or other source on the boundary.
     * @param self Pointer to this BoundaryCondition object.
     * @param ctx The context for the current step.
     * @param ... (variable arguments) e.g., pointers to store output coordinates.
     * @return 0 on success.
     */
    PetscErrorCode (*PlaceSource)(BoundaryCondition *self, BCContext *ctx, ...);

    /**
     * @brief Called once at program exit to free memory allocated by the handler.
     *
     * This function is responsible for freeing the private `data` struct.
     * @param self Pointer to this BoundaryCondition object.
     * @return 0 on success.
     */
    PetscErrorCode (*Destroy)(BoundaryCondition *self);
};

/**
 * @brief A simple struct to hold a key-value string pair parsed from the bcs.dat file.
 */
typedef struct BC_Param_s {
    char *key;
    char *value;
    struct BC_Param_s *next;
} BC_Param;

/**
 * @brief Holds the complete, static configuration for one of the six boundary faces of a block.
 *
 * An array of these structs will live inside the main UserCtx, providing the primary
 * access point for the boundary system.
 */
typedef struct {
    BCFace             face_id;             ///< Geometric ID (e.g., BC_FACE_NEG_Z).
    BCType             mathematical_type;   ///< General type from .dat file (e.g., INLET).
    BCHandlerType      handler_type;        ///< Specific handler from .dat file (e.g., BC_HANDLER_INLET_CONSTANT_VELOCITY).

  BC_Param             *params;               // Head of a linked list of parameters for this face.
    /**
     * @brief Pointer to the live handler object that implements the behavior.
     * This is populated by the BoundarySystem_Create "factory" function based on handler_type.
     */
    BoundaryCondition *handler;

} BoundaryFaceConfig;


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
    PetscReal xMin,yMin,zMin; /// Minimum bounds of the grid.
    PetscReal xMax,yMax,zMax; /// Maximum bounds of the grid.
    PetscInt IM, JM, KM;    ///< Global grid dimensions in x, y, z directions.
    PetscInt nblk;          /// No.of blocks in the grid.  
    BoundingBox bbox;       ///< Bounding box for the local grid domain.
    DMDALocalInfo info;     ///< Local information about the DMDA.
    PetscReal rx;           // Stretching ratio in x-direction
    PetscReal ry;           // Stretching ratio in y-direction
    PetscReal rz;           // Stretching ratio in z-direction

   // Boundary Condition Related fields
  BoundaryFaceConfig boundary_faces[6]; 
    BCType    face_bc_types[6];     //TO BE DEPRECIATED (LEGACY) can be indexed directly using BCFace enum-ex: face_bc_types[BC_FACE_NEG_X]. 
    PetscBool inletFaceDefined;     // Flag: PETSC_TRUE if an INLET face was found in bcs.dat
    BCFace    identifiedInletBCFace; // Stores the BCFace enum corresponding to the INLET

    // Simulation fields
    Vec Ucont,lUcont;              ///< Contravariant velocity field.
    Vec Ucat, lUcat;               ///< Cartesian velocity field.
    Vec P, lP;                  ///< Pressure field.
    Vec Nvert, lNvert;              ///< Node state field (fluid, solid, etc.).
    Vec Nvert_o, lNvert_o;            ///< Node state field in the previous timestep.
    PetscReal ConstantVelocity;
    Cmpnts    ConstantContra;
    PetscReal ConstantPressure;
    PetscReal ConstantNvert;

    // Curvilinear Transformation fields.
   Vec Csi,Eta,Zet;  // Populated in Metrics.c and then DESTROYED.
   Vec ICsi,IEta,IZet;  // Populated in Metrics.c and then DESTROYED.
   Vec JCsi,JEta,JZet;  // Populated in Metrics.c and then DESTROYED.
   Vec KCsi,KEta,KZet;  // Populated in Metrics.c and then DESTROYED.
   Vec lCsi,lEta,lZet;  // VALID & used everywhere.
   Vec lICsi,lIEta,lIZet;  // VALID & used everywhere.
   Vec lJCsi,lJEta,lJZet;  // VALID & used everywhere.
   Vec lKCsi,lKEta,lKZet;  // VALID & used everywhere.

   Vec Aj, IAj, JAj,KAj;
   Vec lAj, lIAj, lJAj, lKAj; 

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
    PetscMPIInt *miglist;      ///< List of ranks for particle migration.
    PetscInt ParticleInitialization;
    PetscInt NumberofParticles;   /// No.of particles in the domain.
    Vec ParticleCount;        /// Count of number of particles in each cell.

    // Simulation parameters
    PetscReal dt;           ///< Time step.
    PetscReal ren;          ///< Reynolds number.
    PetscReal ti;            ///< Current time.
    PetscInt step;            /// Current Timestep Index.
    PetscInt  FieldInitialization;
    PetscInt  VelocityInitializationMethod;
    PetscInt  LoggingFrequency;  // Logging frequency for particle data logging (LOG_PARTICLE_FIELDS)
  
    // Flags for simulation modes
    PetscBool averaging;    ///< Flag to indicate whether statistical averaging is enabled.
    PetscBool les;          ///< Flag to indicate if LES is active.
    PetscBool rans;         ///< Flag to indicate if RANS is active.

    // Parallelization Parameters 
    RankNeighbors neighbors;  // Neighbor ranks
    BoundingBox global_domain_bbox;

    // Miscellaneous fields
    PetscInt _this;         ///< Current block index.
    
  
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

// Structure to hold migration information for a particle
typedef struct {
    PetscInt local_index; // Original local index of the particle on this rank
    PetscInt target_rank; // Rank this particle should migrate to
    // Add PID if needed for debugging: PetscInt64 pid;
} MigrationInfo;

#endif // COMMON_H
