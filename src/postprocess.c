/*****************************************************************************/
/*   postprocess.c

   This file implements a simple post-processing executable that:
    1) Initializes PETSc and MPI.
    2) Reads command-line options for input data (time index, field name, etc.).
    3) Reads .dat file(s) into arrays.
    4) Decides whether to write a .vts or .vtp file based on '-out_ext'.
    5) Constructs a VTKMetaData object and calls CreateVTKFileFromMetadata().
    6) Finalizes PETSc.
*/
/*****************************************************************************/

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscdmcomposite.h>

#include "postprocess.h"

// #include "io.h"
// #include "grid.h"
// #include "interpolation.h"
// #include "ParticleSwarm.h"

#define MAX_FILENAME_LENGTH 256

/**
 * @brief Prepares a VTKMetaData structure for VTK output.
 *
 * This function fills in the VTKMetaData structure for either structured grid
 * output (VTK_STRUCTURED) or polydata output (VTK_POLYDATA). For structured grids,
 * the grid dimensions (mx, my, mz) are obtained from the PETSc options database,
 * and the function verifies that the coordinate array length matches the expected
 * grid size. The field array is then checked to determine if it contains scalar
 * data (one value per grid point) or vector data (three values per grid point).
 *
 * For polydata, the function verifies that the coordinate array length is a multiple
 * of 3 (each point having x, y, and z coordinates), computes the number of points,
 * and determines whether the field array contains scalar or vector data. Additionally,
 * connectivity and offsets arrays are created to define the connectivity of the points.
 *
 * @param[in]  coordsArray   Pointer to the array of coordinate values.
 * @param[in]  Ncoords       Total number of coordinate values in coordsArray.
 * @param[in]  fieldArray    Pointer to the field data (can be scalar or vector data).
 * @param[in]  Nscalars      Total number of values in the fieldArray.
 * @param[in]  fieldName     Name of the field (e.g., "velocity").
 * @param[in]  outExt        The output file extension. "vtp" indicates polydata;
 *                           any other value is interpreted as structured grid.
 * @param[out] meta          Pointer to the VTKMetaData structure that will be populated.
 * @param[in]  rank          MPI rank of the current process (only rank 0 performs full setup).
 *
 * @return PetscErrorCode Returns 0 on success, or an error code on failure.
 *
 * @note For structured grid output, grid dimensions are provided via PETSc options:
 *       -grid_x, -grid_y, -grid_z. The field is interpreted as:
 *         - a vector if Nscalars equals (mx * my * mz * 3), or
 *         - a scalar if Nscalars equals (mx * my * mz).
 */
static PetscErrorCode PrepareVTKMetaData(const double *coordsArray,
                                         PetscInt      Ncoords,
                                         double       *fieldArray,
                                         PetscInt      Nscalars,
                                         const char   *fieldName,
                                         const char   *outExt,
                                         VTKMetaData  *meta,
                                         int           rank)
{
  PetscFunctionBeginUser;
  PetscErrorCode ierr;
  PetscBool      match, isVTP;

  /*-----------------------------------------------------------------------
   * Initialize the VTKMetaData structure by zeroing out all fields.
   *-----------------------------------------------------------------------*/
  memset(meta, 0, sizeof(VTKMetaData));
  LOG_ALLOW(LOCAL, LOG_DEBUG, "VTKMetaData structure initialized.\n");

  /*-----------------------------------------------------------------------
   * Determine the output file type by comparing outExt with "vtp".
   * If outExt is "vtp", then configure for polydata output; otherwise,
   * assume structured grid output.
   *-----------------------------------------------------------------------*/
  ierr = PetscStrcmp(outExt, "vtp", &match);CHKERRQ(ierr);
  isVTP = match ? PETSC_TRUE : PETSC_FALSE;
  LOG_ALLOW(LOCAL, LOG_INFO, "Output file extension '%s' detected; isVTP = %d.\n", outExt, isVTP);

  if (!isVTP) {
    /*=====================================================================
     * Structured Grid Output: VTK_STRUCTURED
     *=====================================================================*/
    LOG_ALLOW(LOCAL, LOG_INFO, "Configuring for VTK_STRUCTURED output.\n");
    meta->fileType = VTK_STRUCTURED;

    /*---------------------------------------------------------------------
     * Retrieve grid dimensions from PETSc options.
     * These dimensions can be provided via:
     *    -grid_x <value> -grid_y <value> -grid_z <value>
     * The grid dimensions must be at least 2 in each direction.
     *---------------------------------------------------------------------*/
    PetscInt mx, my, mz;
    ierr = PetscOptionsGetInt(NULL, NULL, "-im", &mx, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-jm", &my, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-km", &mz, NULL);CHKERRQ(ierr);

    if (mx < 2 || my < 2 || mz < 2) {
      LOG_ALLOW(LOCAL, LOG_WARNING, "Invalid grid dimensions: mx=%d, my=%d, mz=%d. All dimensions must be >= 2.\n", mx, my, mz);
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
              "Grid dimensions must be at least 2 in each direction.");
    }
    LOG_ALLOW(LOCAL, LOG_INFO, "Grid dimensions obtained from options: mx=%d, my=%d, mz=%d.\n", mx, my, mz);

    // These are number of cells , hence we add one for each dimension to get number of grid points

    meta->mx = mx + 1;
    meta->my = my + 1; 
    meta->mz = mz + 1;
    /* Compute the number of  cells (nodes) for the structured grid.*/
    
    meta->nnodes = (mx+1) * (my+1) * (mz+1);
    LOG_ALLOW(LOCAL, LOG_INFO, "Computed number of nodes (cells): %d.\n", meta->nnodes);

    /*---------------------------------------------------------------------
     * Populate the metadata only on rank 0.
     * Validate that the coordsArray length matches the expected number of values.
     * Expected length = (mx+1) *( my+1) * (mz+1) points * 3 coordinates per point.
     *---------------------------------------------------------------------*/
    if (!rank) {
      PetscInt npoints = (mx+1) * (my+1) * (mz+1);
           if (Ncoords != npoints * 3) {
	     LOG_ALLOW(LOCAL, LOG_WARNING, "Coordinates array length mismatch: expected %d, got %d.\n", npoints * 3, Ncoords);
	     SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                 "Coordinates array length does not match grid dimensions.");
        }
      LOG_ALLOW(LOCAL, LOG_INFO, "Coordinates array validated: %d points (%d values).\n", npoints, Ncoords);
      meta->coords = (double*)coordsArray;

      /*-------------------------------------------------------------------
       * Determine the field type based on the field array length.
       * If Nscalars equals npoints * 3, the field is interpreted as a vector field;
       * if Nscalars equals npoints, it is treated as a scalar field.
       *-------------------------------------------------------------------*/
      if (Nscalars == npoints * 3) {
        meta->vectorField      = fieldArray;
        meta->vectorFieldName  = fieldName;
        meta->numVectorFields  = 1;
        meta->numScalarFields  = 0;
        LOG_ALLOW(LOCAL, LOG_INFO, "Field '%s' interpreted as a vector for structured grid.\n", fieldName);
      } else if (Nscalars == npoints) {
        meta->scalarField      = fieldArray;
        meta->scalarFieldName  = fieldName;
        meta->numScalarFields  = 1;
        meta->numVectorFields  = 0;
        LOG_ALLOW(LOCAL, LOG_INFO, "Field '%s' interpreted as a scalar for structured grid.\n", fieldName);
      } else {
        LOG_ALLOW(LOCAL, LOG_WARNING, "Field array length %d does not match expected sizes for grid points (%d or %d).\n", Nscalars, npoints, npoints * 3);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                "Field array length does not match number of grid points for scalar or vector.");
      }
      LOG_ALLOW(LOCAL, LOG_DEBUG, "Structured grid metadata setup complete on rank 0.\n");
    }
  } else {
    /*=====================================================================
     * Polydata Output: VTK_POLYDATA
     *=====================================================================*/
    LOG_ALLOW(LOCAL, LOG_INFO, "Configuring for VTK_POLYDATA output.\n");
    meta->fileType = VTK_POLYDATA;
    if (!rank) {
      /*-------------------------------------------------------------------
       * For polydata, assign the coordinate array pointer and validate its length.
       * The coordinate array length must be a multiple of 3 (x, y, z per point).
       *-------------------------------------------------------------------*/
      meta->coords = (double*)coordsArray;
      if (Ncoords % 3 != 0) {
        LOG_ALLOW(LOCAL, LOG_WARNING, "Coordinates array length %d is not a multiple of 3.\n", Ncoords);
        LOG_ALLOW(LOCAL, LOG_INFO, "Coordinates array: %d values provided.\n", Ncoords);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                "Coordinates length must be multiple of 3.");
      }
      /* Compute the number of points. */
      meta->npoints = Ncoords / 3;
      LOG_ALLOW(LOCAL, LOG_INFO, "Number of points computed from coordinates: %d.\n", meta->npoints);

      /*-------------------------------------------------------------------
       * Identify the field type.
       * If the field array length equals npoints * 3, the field is interpreted as
       * a vector; otherwise, it is treated as a scalar field.
       *-------------------------------------------------------------------*/
      if (Nscalars == meta->npoints * 3) {
        meta->vectorField      = fieldArray;
        meta->vectorFieldName  = fieldName;
        meta->numVectorFields  = 1;
        meta->numScalarFields  = 0;
        LOG_ALLOW(LOCAL, LOG_INFO, "Field '%s' interpreted as a vector for polydata.\n", fieldName);
      } else {
        meta->scalarField      = fieldArray;
        meta->scalarFieldName  = fieldName;
        meta->numScalarFields  = 1;
        meta->numVectorFields  = 0;
        LOG_ALLOW(LOCAL, LOG_INFO, "Field '%s' interpreted as a scalar for polydata.\n", fieldName);
      }

      /*-------------------------------------------------------------------
       * Create connectivity and offsets arrays for polydata.
       * Each point is treated as an individual entity.
       *-------------------------------------------------------------------*/
      meta->connectivity = (int*)calloc(meta->npoints, sizeof(int));
      meta->offsets      = (int*)calloc(meta->npoints, sizeof(int));
      for (int i = 0; i < meta->npoints; i++) {
        meta->connectivity[i] = i;
        meta->offsets[i]      = i + 1;
      }
      LOG_ALLOW(LOCAL, LOG_INFO, "Connectivity and offsets arrays created for %d points.\n", meta->npoints);
    }
  }

  LOG_ALLOW(LOCAL, LOG_DEBUG, "PrepareVTKMetaData - Setup complete.\n");
  PetscFunctionReturn(0);
}

static PetscErrorCode CreateAndWriteVTKFile(const char   *outFile,
                                            VTKMetaData  *meta,
                                            MPI_Comm      comm)
{
  PetscFunctionBeginUser;
  LOG_ALLOW(LOCAL, LOG_DEBUG,
            "CreateAndWriteVTKFile - Calling CreateVTKFileFromMetadata for '%s'.\n", outFile);

  int err = CreateVTKFileFromMetadata(outFile, meta, comm);
  if (err) {
    LOG_ALLOW_SYNC(LOCAL, LOG_ERROR,
              "CreateAndWriteVTKFile - CreateVTKFileFromMetadata returned code=%d.\n", err);
    SETERRQ(comm, PETSC_ERR_FILE_WRITE,
            "CreateVTKFileFromMetadata returned an error. \n");
  }
  LOG_ALLOW(LOCAL, LOG_DEBUG, "CreateAndWriteVTKFile - Successfully wrote VTK file '%s'.\n", outFile);
  PetscFunctionReturn(0);
}


/**
 * @brief Parses post-processing parameters from command-line options.
 *
 * This function sets the start time, end time, and time step for batch post-processing.
 * If only `-startTime` is provided, the function sets:
 *   - `endTime = startTime`
 *   - `timeStep = 0`  (indicating a single timestep)
 *
 * @param[out] pps  Pointer to `PostProcessSettings` structure to populate.
 *
 * @return PetscErrorCode  Returns 0 on success, or a PETSc error code on failure.
 */
static PetscErrorCode ParsePostProcessingParams(PostProcessParams *pps)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;
    PetscBool startTimeSet = PETSC_FALSE, endTimeSet = PETSC_FALSE, timeStepSet = PETSC_FALSE;

    // 1) Initialize defaults
    pps->startTime = 0;
    pps->endTime   = 1;
    pps->timeStep  = 1;

    // 2) Parse command-line options
    ierr = PetscOptionsGetInt(NULL, NULL, "-startTime", &pps->startTime, &startTimeSet);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-endTime", &pps->endTime, &endTimeSet);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-timeStep", &pps->timeStep, &timeStepSet);CHKERRQ(ierr);

    // 3) Handle missing values
    if (startTimeSet && !endTimeSet) {
        pps->endTime = pps->startTime;  // Single timestep mode
    }

    if (startTimeSet && !timeStepSet) {
        pps->timeStep = 0;  // Single timestep mode
    }

    // 4) Set default processing fields and file extension
    pps->outputParticles = PETSC_TRUE;     // Flag to Output(or not) Particle Fields.
    strcpy(pps->eulerianExt, "vts");  // Eulerian File extension
    strcpy(pps->particleExt, "vtp");  // Lagrangian File extension
    strcpy(pps->eulerianPrefix,"viz"); // Prefix(directory to save eulerian fields in.
    strcpy(pps->particlePrefix,"viz"); // Prefix(directory to save eulerian fields in.                         

    // 5) Log parsed values
    LOG_ALLOW(LOCAL,LOG_INFO, "Parsed PostProcessingParams: startTime=%d, endTime=%d, timeStep=%d\n",
                pps->startTime, pps->endTime, pps->timeStep);

    PetscFunctionReturn(0);
}

/**
 * @brief Gathers a PETSc vector or a DMSwarm field, prepares VTKMetaData, and writes to a VTK file.
 *
 * This function dynamically handles **both Eulerian (.vts) and Lagrangian (.vtp) data**:
 *  - **If the field exists in DMSwarm**, it **extracts it into a PETSc Vec**.
 *  - **If already a Vec**, it directly processes it.
 *  - **For `.vtp` files, it also gathers particle positions** (from the swarm field "position").
 *  - **For `.vts` files, the coordinate vector ("coor") is gathered** and passed to the metadata.
 *
 * @param[in]  user       Pointer to user context (contains swarm and Eulerian fields).
 * @param[in]  fieldName  Name of the field (e.g., "Ucat" for vts, "Velocity" for vtp).
 * @param[in]  timeIndex  Timestep index for filename (e.g., "Ucat_00010.vts").
 * @param[in]  outExt     File extension ("vts" for structured, "vtp" for particles).
 * @param[in]  prefix     File prefix (directory to store file in).
 * @param[in]  comm       MPI communicator (usually PETSC_COMM_WORLD).
 *
 * @return PetscErrorCode  0 on success, or an error code on failure.
 */
PetscErrorCode GatherAndWriteField(UserCtx *user,
                                   const char *fieldName,
                                   PetscInt timeIndex,
                                   const char *outExt,
                                   const char *prefix,
                                   MPI_Comm comm)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "GatherAndWriteField - Rank %d processing field '%s' at time index %d with extension '%s'.\n",
              rank, fieldName, (int)timeIndex, outExt);

    Vec fieldVec = NULL;
    PetscBool isSwarmField = PETSC_FALSE;

    /* 1) Detect if the field exists in DMSwarm (for ".vtp" cases) */
    if (strcmp(outExt, "vtp") == 0) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAndWriteField - Attempting to create global vector from swarm field '%s'.\n", fieldName);
        ierr = DMSwarmCreateGlobalVectorFromField(user->swarm, fieldName, &fieldVec);
        if (ierr) {
            /* An error occurred so assume the field does not exist in the swarm.
               Clear the error code and mark isSwarmField as false. */
            LOG_ALLOW(LOCAL, LOG_WARNING, "GatherAndWriteField - Field '%s' not found in swarm, treating as Eulerian field.\n", fieldName);
            isSwarmField = PETSC_FALSE;
            ierr = 0;  /* Reset the error code so that we can continue */
        } else {
            isSwarmField = PETSC_TRUE;
            LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAndWriteField - Successfully created swarm vector for field '%s'.\n", fieldName);
        }
    }

    /* 2) If not a swarm field, assume it's a normal Eulerian field */
    if (!isSwarmField) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAndWriteField - Using Eulerian field for '%s'.\n", fieldName);
        if (strcmp(fieldName, "Ucat") == 0) fieldVec = user->Ucat;
        else if (strcmp(fieldName, "P") == 0) fieldVec = user->P;
        else if (strcmp(fieldName, "Ucont") == 0) fieldVec = user->Ucont;
        else if (strcmp(fieldName, "Nvert") == 0) fieldVec = user->Nvert;
        else SETERRQ(comm, PETSC_ERR_ARG_WRONG, "Unknown field requested.");
    }

    /* 3) Gather the field data (Eulerian OR swarm field) */
    PetscInt  Nfield = 0;
    double   *fieldArray = NULL;
    if (fieldVec) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAndWriteField - Gathering field data for '%s'.\n", fieldName);
        ierr = VecToArrayOnRank0(fieldVec, &Nfield, &fieldArray);CHKERRQ(ierr);
        LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAndWriteField - Field data gathered: Nfield=%D.\n", Nfield);
    }

    /* 4) Gather coordinate data based on file type.
     * For ".vtp", gather particle positions from the swarm ("position" field).
     * For ".vts", gather the Eulerian coordinate field ("coor") from the user context.
     */
    PetscInt  Ncoords = 0;
    double   *coordsArray = NULL;
    Vec coordsVec;
    if (strcmp(outExt, "vtp") == 0) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAndWriteField - Gathering particle positions from swarm for .vtp output.\n");
        ierr = DMSwarmCreateGlobalVectorFromField(user->swarm, "position", &coordsVec);CHKERRQ(ierr);
        ierr = VecToArrayOnRank0(coordsVec, &Ncoords, &coordsArray);CHKERRQ(ierr);
        ierr = DMSwarmDestroyGlobalVectorFromField(user->swarm, "position", &coordsVec);CHKERRQ(ierr);
        LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAndWriteField - Particle positions gathered: Ncoords=%D.\n", Ncoords);
    } else if (strcmp(outExt, "vts") == 0) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAndWriteField - Gathering Eulerian coordinates for .vts output.\n");
         /* For structured grid, get the coordinates vector from DMDA */
         ierr = DMGetCoordinatesLocal(user->da, &coordsVec);CHKERRQ(ierr);
         ierr = VecToArrayOnRank0(coordsVec, &Ncoords, &coordsArray);CHKERRQ(ierr);
         LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAndWriteField - Eulerian coordinates gathered: Ncoords=%D.\n", Ncoords);
    }

    /* 5) Prepare VTK metadata.
     * For ".vts" files, coordsArray now contains the Eulerian grid coordinates.
     */
    VTKMetaData meta;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAndWriteField - Preparing VTK metadata.\n");
    ierr = PrepareVTKMetaData(coordsArray, Ncoords,  /* coordsArray will be non-NULL for "vts" */
                              fieldArray, Nfield,
                              fieldName, outExt,
                              &meta, rank);CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAndWriteField - VTK metadata prepared.\n");

    /* 6) Construct the filename (e.g., "Ucat_00010.vts" OR "particleVelocity_00010.vtp") */
    char outFile[256];
    snprintf(outFile, MAX_FILENAME_LENGTH, "%s/%s_%05d.%s", prefix, fieldName, (int)timeIndex, outExt);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAndWriteField - Constructed output filename: %s\n", outFile);

    /* 7) Write the VTK file */
    LOG_ALLOW(LOCAL, LOG_INFO, "GatherAndWriteField - Writing VTK file '%s'.\n", outFile);
    ierr = CreateAndWriteVTKFile(outFile, &meta, comm);CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "GatherAndWriteField - Successfully wrote VTK file '%s'.\n", outFile);

    /* 8) Cleanup on rank 0 */
    if (rank == 0) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAndWriteField - Cleaning up rank 0 memory.\n");
        if (meta.fileType == VTK_POLYDATA) {
            free(meta.connectivity);
            free(meta.offsets);
        }
        free(coordsArray);
        free(fieldArray);
    }

    /* 9) Cleanup DMSwarm Vec reference (if created) */
    if (isSwarmField) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "GatherAndWriteField - Cleaning up swarm field vector for '%s'.\n", fieldName);
        ierr = DMSwarmDestroyGlobalVectorFromField(user->swarm, fieldName, &fieldVec);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


/**
 * @brief Writes up to four Eulerian fields (Ucat, P, Ucont, Nvert) using the
 *        GatherAndWriteField() wrapper. Each field becomes a separate file.
 *
 * @param[in]  user       Pointer to the UserCtx holding PETSc vectors.
 * @param[in]  timeIndex  Timestep index for filenames.
 * @param[in]  outExt     File extension (e.g., "vts" or "vtp").
 * @param[in]  prefix     File prefix(directory to store file in).
 *
 * @return PetscErrorCode  Returns 0 on success, error code on failure.
 */
PetscErrorCode WriteEulerianVTK(UserCtx *user, PetscInt timeIndex, const char *outExt, const char *prefix)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    /*
       1) For each field, call GatherAndWriteField() if it's not NULL.
       2) This approach is concise and avoids repeated code.
    */

    // Field #1: Ucat
    if (user->Ucat) {
      ierr = GatherAndWriteField(user, "Ucat", timeIndex, outExt, prefix,PETSC_COMM_WORLD);CHKERRQ(ierr);
    }

    // Field #2: P
    if (user->P) {
      ierr = GatherAndWriteField(user,    "P",    timeIndex, outExt, prefix,PETSC_COMM_WORLD);CHKERRQ(ierr);
    }

    // Field #3: Ucont
    if (user->Ucont) {
      ierr = GatherAndWriteField(user, "Ucont", timeIndex, outExt, prefix,PETSC_COMM_WORLD);CHKERRQ(ierr);
    }

    // Field #4: Nvert
    if (user->Nvert) {
      //  ierr = GatherAndWriteField(user, "Nvert", timeIndex, outExt, prefix,PETSC_COMM_WORLD);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}


/**
 * @brief Writes particle (Lagrangian) data fields to separate .vtp files.
 *
 * This function loops over key DMSwarm fields (positions, velocity, etc.) and writes 
 * each as a separate .vtp file, ensuring modularity and flexibility.
 *
 * @param[in] user       Pointer to the user context (contains the particle swarm).
 * @param[in] timeIndex  Current timestep index (used in filenames).
 * @param[in] outExt     File extension ("vtp" for particles).
 * @param[in]  prefix     File prefix(directory to store file in).
 *
 * @return PetscErrorCode  0 on success, nonzero on failure.
 */
PetscErrorCode WriteParticleVTK(UserCtx *user, PetscInt timeIndex, const char *outExt,const char *prefix)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    // 1) Velocity field (assumes particles have a "Velocity" field)
    ierr = GatherAndWriteField(user, "velocity", timeIndex, outExt, prefix,PETSC_COMM_WORLD);CHKERRQ(ierr);

    // 2) Additional fields, e.g., Temperature (optional, add more as needed)
    //   ierr = GatherAndWriteSwarmField(user, "Temperature",
    //                                "particleTemperature",
    //                                timeIndex, outExt,
    //                                PETSC_COMM_WORLD);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


int main(int argc, char **argv)
{

    UserCtx *user = NULL;     // User context
    PetscErrorCode ierr;      // PETSc error handling
    PetscInt block_number = 1;
    PetscInt rstart = 0;
    PetscInt np = 0;
    PetscInt particlesPerProcess = 0;
    PetscInt ti = 0;
    PetscInt rank, size;
    BoundingBox *bboxlist = NULL;   // Array of bounding boxes
    PostProcessParams pps;
    static char help[] = " Postprocessing Tool - swarm-curvIB";

    //    PetscBool ParticleOut = PETSC_FALSE;

    // -------------------- PETSc Initialization --------------------
    ierr = PetscInitialize(&argc, &argv, (char *)0, help); CHKERRQ(ierr);

    // -------------------- Setup Logging Allow-List ----------------
    // Only these function names will produce LOG_ALLOW (or LOG_ALLOW_SYNC) output.
    // You can add more as needed (e.g., "InitializeSimulation", "PerformParticleSwarmOperations", etc.).
    const char *allowedFuncs[] = {
      // "SetupGridAndVectors",
      // "DefineGridCoordinates",
      // "AssignGridCoordinates"
      // "WriteEulerianVTK",
      "GatherAndWriteField",
      "VecToArrayOnRank0",
      "PrepareVTKMetaData",
     "CreateAndWriteVTKFile",
      "CreateVTKFileFromMetadata",
      // "InitializeSimulation",
      // "ReadSwarmField",
      // "ReadFieldData",
      // "ReadAllSwarmFields",
      // "InitializeSwarm"
       "main"                 
      // "InitializeSimulation",
      // "CreateParticleSwarm" 
      // "ParsePostProcessingParams",
      // "ReadSimulationFields",
    };
    set_allowed_functions(allowedFuncs, 6);

    print_log_level();

    // Parse the Post Processing parameters.
    ierr = ParsePostProcessingParams(&pps);

    // Initialize simulation: user context, MPI rank/size, np, etc.
    ierr = InitializeSimulation(&user, &rank, &size, &np, &rstart, &ti, &block_number); CHKERRQ(ierr);

    // Setup the computational grid
    ierr = SetupGridAndVectors(user, block_number); CHKERRQ(ierr);

    if(pps.outputParticles){
      // Create the Particle Swarm
      ierr = CreateParticleSwarm(user,np,&particlesPerProcess,bboxlist);
    }

    // Loop over timesteps from startTime to endTime, step by pps.timeStep. 
    for (PetscInt ti = pps.startTime; ti < pps.endTime; ti += pps.timeStep) {
      // We store the current timestep in user->ti if needed. 

         user->ti = (PetscInt) ti;
	 LOG_ALLOW(LOCAL,LOG_INFO, "  Output Timestep : %d \n",ti);

	 // Read the eulerian fields from data files.
	 ierr = ReadSimulationFields(user,ti); CHKERRQ(ierr);
	 LOG_ALLOW(LOCAL,LOG_INFO, " Eulerian Fields Read \n",ti);

         // Write Eulerian field data to vts files
         ierr =  WriteEulerianVTK(user, ti, pps.eulerianExt,pps.eulerianPrefix);

	 if(pps.outputParticles){
	 // Read lagrangian fields from data files.
	   ierr = ReadAllSwarmFields(user,ti);	 
	   LOG_ALLOW(LOCAL,LOG_INFO, " Particle Fields Read \n",ti);
    
	 // Write Lagrangian/particle data to vtp file/files.
	   ierr = WriteParticleVTK(user, ti, pps.particleExt,pps.particlePrefix); 
	   } // Particles Output loop.

	 } // timestep loop.
	 
        ierr = FinalizeSimulation(user,block_number,bboxlist);

    return 0;
}


