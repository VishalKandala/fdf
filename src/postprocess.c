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
#include "logging.h"
#include "common.h"
// #include "io.h"
// #include "grid.h"
// #include "interpolation.h"
// #include "ParticleSwarm.h"

/**
 * @brief Fills a VTKMetaData structure for either structured or polydata output.
 *
 * This function inspects \p outExt ("vts" or "vtp") to decide whether to create
 * \c VTK_STRUCTURED or \c VTK_POLYDATA. In the latter case, it also sets up
 * connectivity and offsets arrays for each point.
 *
 * @param[in]  coordsArray  The gathered array of coordinates (valid on rank 0).
 * @param[in]  Ncoords      The length of the coordinate array on rank 0.
 * @param[in]  fieldArray   The gathered array of field data (scalar or vector).
 * @param[in]  Nscalars     The length of the field array on rank 0.
 * @param[in]  fieldName    Name of the field (e.g., "velocity").
 * @param[in]  outExt       The desired file extension (e.g., "vtp" or "vts").
 * @param[out] meta         The \c VTKMetaData structure to fill.
 * @param[in]  rank         The current MPI rank (only rank 0 populates \p meta fully).
 *
 * @return PetscErrorCode  Returns 0 on success, or an error if the coords array
 *                         is not a multiple of 3 (polydata) or if memory fails.
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

  memset(meta, 0, sizeof(VTKMetaData));

  // Check if outExt is "vtp"
  ierr = PetscStrcmp(outExt, "vtp", &match);CHKERRQ(ierr);
  isVTP = match ? PETSC_TRUE : PETSC_FALSE;

  // Log decision for structured vs polydata
  if (!isVTP) {
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "PrepareVTKMetaData - Configuring for VTK_STRUCTURED output.\n");
    meta->fileType = VTK_STRUCTURED;
    // Example sizes for structured data
    meta->mx = 4; meta->my = 4; meta->mz = 4;
    meta->nnodes = (meta->mx - 1) * (meta->my - 1) * (meta->mz - 1);

    if (!rank) {
      meta->coords            = (double*)coordsArray;
      meta->scalarField       = fieldArray;
      meta->scalarFieldName   = fieldName;
      meta->numScalarFields   = 1;
      meta->numVectorFields   = 0;
    }
  } else {
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "PrepareVTKMetaData - Configuring for VTK_POLYDATA output.\n");
    meta->fileType = VTK_POLYDATA;
    if (!rank) {
      meta->coords = (double*)coordsArray;
      if (Ncoords % 3 != 0) {
        LOG_ALLOW(GLOBAL, LOG_ERROR,
                  "PrepareVTKMetaData - coords array length %D not multiple of 3.\n",
                  Ncoords);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                "Coordinates length must be multiple of 3.");
      }
      meta->npoints = Ncoords / 3;

      // Check if the field data is vector or scalar
      if (Nscalars == meta->npoints * 3) {
        meta->vectorField      = fieldArray;
        meta->vectorFieldName  = fieldName;
        meta->numVectorFields  = 1;
        meta->numScalarFields  = 0;
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "PrepareVTKMetaData - Interpreting field '%s' as a vector.\n", fieldName);
      } else {
        meta->scalarField      = fieldArray;
        meta->scalarFieldName  = fieldName;
        meta->numScalarFields  = 1;
        meta->numVectorFields  = 0;
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "PrepareVTKMetaData - Interpreting field '%s' as a scalar.\n", fieldName);
      }

      // Create connectivity & offsets for polydata
      meta->connectivity = (int*)calloc(meta->npoints, sizeof(int));
      meta->offsets      = (int*)calloc(meta->npoints, sizeof(int));
      for (int i = 0; i < meta->npoints; i++) {
        meta->connectivity[i] = i;
        meta->offsets[i]      = i + 1;
      }
      LOG_ALLOW(GLOBAL, LOG_DEBUG,
                "PrepareVTKMetaData - Created connectivity & offsets for %d points.\n",
                meta->npoints);
    }
  }

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "PrepareVTKMetaData - Setup complete.\n");
  PetscFunctionReturn(0);
}

/**
 * @brief Constructs an output file path of the form "results/<fieldName><timeIndex>.<outExt>".
 *
 * @param[in]  fieldName   The name of the field (e.g., "velocity").
 * @param[in]  timeIndex   The time index to be appended (zero-padded).
 * @param[in]  outExt      The desired file extension (e.g., "vtp" or "vts").
 * @param[out] outFile     Buffer to store the resulting file path.
 * @param[in]  outFileSize The size of \p outFile buffer.
 *
 * @return PetscErrorCode  Returns 0 on success, or error from \c PetscSNPrintf.
 */
static PetscErrorCode ConstructOutputFilename(const char *fieldName,
                                              PetscInt    timeIndex,
                                              const char *outExt,
                                              char       *outFile,
                                              size_t      outFileSize)
{
  PetscFunctionBeginUser;
  PetscErrorCode ierr;
  ierr = PetscSNPrintf(outFile, outFileSize,
                       "results/%s%05d.%s", fieldName, timeIndex, outExt);CHKERRQ(ierr);
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ConstructOutputFilename - Created output file '%s'.\n", outFile);
  PetscFunctionReturn(0);
}

/**
 * @brief Invokes CreateVTKFileFromMetadata to write a VTK file, then checks for errors.
 *
 * @param[in] outFile  The file name/path for the VTK file to create.
 * @param[in] meta     Pointer to the initialized \c VTKMetaData describing the data layout.
 * @param[in] comm     The MPI communicator.
 *
 * @return PetscErrorCode  Returns 0 on success, or a PETSc error if the underlying
 *                         function call returns a non-zero integer.
 */
static PetscErrorCode CreateAndWriteVTKFile(const char   *outFile,
                                            VTKMetaData  *meta,
                                            MPI_Comm      comm)
{
  PetscFunctionBeginUser;
  LOG_ALLOW(GLOBAL, LOG_DEBUG,
            "CreateAndWriteVTKFile - Calling CreateVTKFileFromMetadata for '%s'.\n", outFile);

  int err = CreateVTKFileFromMetadata(outFile, meta, comm);
  if (err) {
    LOG_ALLOW(GLOBAL, LOG_ERROR,
              "CreateAndWriteVTKFile - CreateVTKFileFromMetadata returned code=%d.\n", err);
    SETERRQ(comm, PETSC_ERR_FILE_WRITE,
            "CreateVTKFileFromMetadata returned an error. \n");
  }
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "CreateAndWriteVTKFile - Successfully wrote VTK file '%s'.\n", outFile);
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
static PetscErrorCode ParsePostProcessingParams(PostProcessSettings *pps)
{
    PetscFunctionBeginUser;

    PetscErrorCode ierr;
    PetscBool startTimeSet = PETSC_FALSE, endTimeSet = PETSC_FALSE, timeStepSet = PETSC_FALSE;

    // 1) Initialize defaults
    pps->startTime = 0;
    pps->endTime   = 0;
    pps->timeStep  = 0;

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
    pps->numFields = 2; 
    strcpy(pps->fieldNames[0], "Ucat");   // e.g., velocity field
    strcpy(pps->fieldNames[1], "P");      // e.g., pressure field
    strcpy(pps->ext, "dat");              // File extension

    // 5) Log parsed values
    PetscPrintf(PETSC_COMM_WORLD,
                "[INFO] Parsed PostProcessingParams: startTime=%d, endTime=%d, timeStep=%d\n",
                pps->startTime, pps->endTime, pps->timeStep);

    PetscFunctionReturn(0);
}

/**
 * @brief Gathers a PETSc vector or a DMSwarm field, prepares VTKMetaData, and writes to a VTK file.
 *
 * This function dynamically handles **both Eulerian (.vts) and Lagrangian (.vtp) data**:
 *  - **If the field exists in `DMSwarm`**, it **extracts it into a PETSc Vec**.
 *  - **If already a Vec**, it directly processes it.
 *  - **For `.vtp` files, it also gathers particle positions** (`DMSwarmPICField_coor`).
 *  - **For `.vts`, only the field is gathered** (structured grid).
 *
 * @param[in]  user       Pointer to user context (contains swarm and Eulerian fields).
 * @param[in]  fieldName  Name of the field (e.g., "Ucat" for vts, "Velocity" for vtp).
 * @param[in]  timeIndex  Timestep index for filename (e.g., "Ucat_00010.vts").
 * @param[in]  outExt     File extension ("vts" for structured, "vtp" for particles).
 * @param[in]  comm       MPI communicator (usually PETSC_COMM_WORLD).
 *
 * @return PetscErrorCode  0 on success, or an error code on failure.
 */
PetscErrorCode GatherAndWriteField(UserCtx *user,
                                   const char *fieldName,
                                   PetscInt timeIndex,
                                   const char *outExt,
                                   MPI_Comm comm)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);

    Vec fieldVec = NULL;
    PetscBool isSwarmField = PETSC_FALSE;

    // 1) Detect if the field exists in DMSwarm (for ".vtp" cases)
    if (strcmp(outExt, "vtp") == 0) {
        ierr = DMSwarmHasField(user->swarm, fieldName, &isSwarmField);CHKERRQ(ierr);
        if (isSwarmField) {
            // Create Vec from Swarm field
            ierr = DMSwarmCreateGlobalVectorFromField(user->swarm, fieldName, &fieldVec);CHKERRQ(ierr);
        }
    }

    // 2) If not a swarm field, assume it's a normal Eulerian field
    if (!isSwarmField) {
        if (strcmp(fieldName, "Ucat") == 0) fieldVec = user->Ucat;
        else if (strcmp(fieldName, "P") == 0) fieldVec = user->P;
        else if (strcmp(fieldName, "Ucont") == 0) fieldVec = user->Ucont;
        else if (strcmp(fieldName, "Nvert") == 0) fieldVec = user->Nvert;
        else SETERRQ(comm, PETSC_ERR_ARG_UNKNOWN, "Unknown field requested.");
    }

    // 3) Gather the field data (Eulerian OR swarm field)
    PetscInt  Nfield = 0;
    double   *fieldArray  = NULL;
    if (fieldVec) {
        ierr = VecToArrayOnRank0(fieldVec, &Nfield, &fieldArray);CHKERRQ(ierr);
    }

    // 4) If ".vtp", also gather particle positions from the swarm
    PetscInt  Ncoords = 0;
    double   *coordsArray = NULL;
    if (strcmp(outExt, "vtp") == 0) {
        Vec coordsVec;
        ierr = DMSwarmCreateGlobalVectorFromField(user->swarm, "DMSwarmPICField_coor", &coordsVec);CHKERRQ(ierr);
        ierr = VecToArrayOnRank0(coordsVec, &Ncoords, &coordsArray);CHKERRQ(ierr);
        ierr = DMSwarmDestroyGlobalVectorFromField(user->swarm, "DMSwarmPICField_coor", &coordsVec);CHKERRQ(ierr);
    }

    // 5) Prepare VTK metadata
    VTKMetaData meta;
    ierr = PrepareVTKMetaData(coordsArray, Ncoords,  // If "vts", coordsArray=NULL
                              fieldArray, Nfield,
                              fieldName, outExt,
                              &meta, rank);CHKERRQ(ierr);

    // 6) Construct the filename (e.g., "Ucat_00010.vts" OR "particleVelocity_00010.vtp")
    char outFile[256];
    sprintf(outFile, "%s_%05d.%s", fieldName, (int)timeIndex, outExt);

    // 7) Write the VTK file
    ierr = CreateAndWriteVTKFile(outFile, &meta, comm);CHKERRQ(ierr);

    // 8) Cleanup on rank 0
    if (!rank) {
        if (meta.fileType == VTK_POLYDATA) {
            free(meta.connectivity);
            free(meta.offsets);
        }
        free(coordsArray);
        free(fieldArray);
    }

    // 9) Cleanup DMSwarm Vec reference (if created)
    if (isSwarmField) {
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
 *
 * @return PetscErrorCode  Returns 0 on success, error code on failure.
 */
PetscErrorCode WriteEulerianVTK(UserCtx *user, PetscInt timeIndex, const char *outExt)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    /*
       1) For each field, call GatherAndWriteField() if it's not NULL.
       2) This approach is concise and avoids repeated code.
    */

    // Field #1: Ucat
    if (user->Ucat) {
      ierr = GatherAndWriteField(user, "Ucat", timeIndex, outExt, PETSC_COMM_WORLD);CHKERRQ(ierr);
    }

    // Field #2: P
    if (user->P) {
      ierr = GatherAndWriteField(user,    "P",    timeIndex, outExt, PETSC_COMM_WORLD);CHKERRQ(ierr);
    }

    // Field #3: Ucont
    if (user->Ucont) {
      ierr = GatherAndWriteField(user, "Ucont", timeIndex, outExt, PETSC_COMM_WORLD);CHKERRQ(ierr);
    }

    // Field #4: Nvert
    if (user->Nvert) {
      ierr = GatherAndWriteField(user, "Nvert", timeIndex, outExt, PETSC_COMM_WORLD);CHKERRQ(ierr);
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
 *
 * @return PetscErrorCode  0 on success, nonzero on failure.
 */
PetscErrorCode WriteParticleVTK(UserCtx *user, PetscInt timeIndex, const char *outExt)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    // 1) Positions (DMSwarm typically stores this under "DMSwarmPICField_coor")
    ierr = GatherAndWriteField(user, "position", timeIndex, outExt, PETSC_COMM_WORLD);CHKERRQ(ierr);

    // 2) Velocity field (assumes particles have a "Velocity" field)
    ierr = GatherAndWriteField(user, "velocity", timeIndex, outExt, PETSC_COMM_WORLD);CHKERRQ(ierr);

    // 3) Additional fields, e.g., Temperature (optional, add more as needed)
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
    PetscInt ti = 0;
    PetscInt rank, size;
    BoundingBox *bboxlist;    // Array of bounding boxes
    PostProcessParams pps;
    PetscBool ParticleOut = PETSC_FALSE;
    static char help[] = " Postprocessing Tool - swarm-curvIB";



    // -------------------- PETSc Initialization --------------------
    ierr = PetscInitialize(&argc, &argv, (char *)0, help); CHKERRQ(ierr);





    // -------------------- Setup Logging Allow-List ----------------
    // Only these function names will produce LOG_ALLOW (or LOG_ALLOW_SYNC) output.
    // You can add more as needed (e.g., "InitializeSimulation", "PerformParticleSwarmOperations", etc.).
    const char *allowedFuncs[] = {
        "main",                 // We'll allow logging from this main function
        "SetupGridAndVectors",
        "InitializeSimulation", // Example: also allow logs from InitializeSimulation
        // "BroadcastAllBoundingBoxes",    // Uncomment to allow logs from that function, etc.
        // "PerformParticleSwarmOperations"
    };
    set_allowed_functions(allowedFuncs, 2);




    // -------------------- Demonstrate LOG_ALLOW in main -----------
    // This message will only be printed if "main" is in the allow-list
    // AND if the LOG_LEVEL environment variable is high enough (e.g., INFO or DEBUG).
    LOG_ALLOW(GLOBAL, LOG_INFO, "==> Starting main with function-based logging...\n");



    // Check if user requested to process particle data or grid data.
    ierr = PetscOptionsGetBool(NULL, NULL, "-p", &ParticleOut, NULL); CHKERRQ(ierr);


    // Initialize simulation: user context, MPI rank/size, np, etc.
    ierr = InitializeSimulation(&user, &rank, &size, &np, &rstart, &ti, &block_number); CHKERRQ(ierr);


    // Setup the computational grid
    ierr = SetupGridAndVectors(user, block_number); CHKERRQ(ierr);


    // Create the Particle Swarm
    ierr = CreateParticleSwarm(user,np,bboxlist);

    
    ierr = ParsePostProcessingParams(pps);

    // Loop over timesteps from startTime to endTime, step by pps.timeStep. 
    for (PetscInt ti = pps.startTime; ti <= pps.endTime; ti += pps.timeStep) {
      // We store the current timestep in user->ti if needed. 

         user->ti = (PetscInt) ti;
	 
	 // Read the eulerian fields from data files.
	 ierr = ReadSimulationFields(user,ti); CHKERRQ(ierr);

	 // Read lagrangian fields from data files.
	 ierr = ReadAllSwarmFields(user,ti);	 
	 
         // Write Eulerian field data to vts files
         ierr =  WriteEulerianVTK(user, ti, "vts"); 
    
	 // Write Lagrangian/particle data to vtp file/files.
         ierr = WriteParticleVTK(user, ti, "vtp"); 
	 
	 } // timestep loop.
	 
        ierr = FinalizeSimulation(user,block_number,bboxlist);

    PetscFinalize();
    return 0;
}


///////////////////////////// Single VTP Implementation ////////

/**
 * @brief Parses command-line options using PETSc's option database.
 *
 * This function initializes PETSc (via \c PetscInitialize) and extracts the
 * user-specified time index, field name, and output extension from the
 * command line. Defaults are assigned if these are not specified.
 *
 * @param[in]     argc           Number of command-line arguments.
 * @param[in]     argv           The command-line argument strings.
 * @param[out]    timeIndex      The time index (default = 0).
 * @param[out]    fieldName      Buffer to store the parsed field name (default "velocity").
 * @param[out]    outExt         Buffer to store the file extension (default "vtp").
 * @param[in]     fieldNameSize  The maximum size of \p fieldName buffer.
 * @param[in]     outExtSize     The maximum size of \p outExt buffer.
 *
 * @return PetscErrorCode  Returns 0 on success, otherwise an error code from PETSc.
 */
/*
static PetscErrorCode ParseCommandLineOptions(int argc, char **argv,
                                              int *timeIndex,
                                              char *fieldName,
                                              char *outExt,
                                              size_t fieldNameSize,
                                              size_t outExtSize)
{
  PetscFunctionBeginUser;

  PetscErrorCode ierr;

  *timeIndex = 0; // default
  PetscStrcpy(fieldName, "velocity");
  PetscStrcpy(outExt,    "vtp");

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ParseCommandLineOptions - Parsing PETSc options.\n");
  ierr = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(ierr);

  // parse
  PetscOptionsGetInt(NULL, NULL, "-ti", timeIndex, NULL);
  PetscOptionsGetString(NULL, NULL, "-field", fieldName, fieldNameSize, NULL);
  PetscOptionsGetString(NULL, NULL, "-out_ext", outExt, outExtSize, NULL);

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "ParseCommandLineOptions - Completed parsing.\n");
  PetscFunctionReturn(0);
}
*/



 /*
 * @brief Main entry point for the PETSc-based VTK post-processing tool.
 *
 * This function reads coordinate and field data from distributed PETSc \c Vecs,
 * gathers them onto rank 0, and writes them as a VTK file (either \c .vts or \c .vtp).
 * The file name, time index, and field name are configured via command-line options.
 *
 * Usage:
 * \code
 *   mpirun -n X ./this_executable -ti 10 -field velocity -out_ext vtp
 * \endcode
 *
 * @param[in] argc Number of command-line arguments.
 * @param[in] argv Array of command-line argument strings.
 *
 * @return int Returns 0 on success, or an error code upon failures.
 */
/*
int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  int            rank, size;
  int            ti;
  char           field_name[64];
  char           out_ext[16];
  double        *coordsArray = NULL;
  double        *scalarArray = NULL;
  PetscInt       Ncoords     = 0;
  PetscInt       Nscalars    = 0;

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "main - Starting PETSc-based VTK post-processing.\n");

  // 1) Parse command-line options and initialize PETSc 
  ierr = ParseCommandLineOptions(argc, argv, &ti,
                                 field_name, out_ext,
                                 sizeof(field_name), sizeof(out_ext));CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
  LOG_ALLOW(GLOBAL, LOG_INFO, "main - PETSc initialized. rank=%d of %d.\n", rank, size);

  // 2) Build a user context 
  UserCtx user;

  ierr = BuildUserContext(&user);CHKERRQ(ierr);

  // 3) Read coordinate data into coordsArray 
  ierr = ReadPositions(ti, &user, &coordsArray, &Ncoords);CHKERRQ(ierr);

  // 4) Read the field data into scalarArray 
  ierr = ReadFieldDataWrapper(ti, field_name, &user, &scalarArray, &Nscalars);CHKERRQ(ierr);

  // 5) Prepare the VTKMetaData struct 
  VTKMetaData meta;
  ierr = PrepareVTKMetaData(coordsArray, Ncoords,
                            scalarArray, Nscalars,
                            field_name, out_ext,
                            &meta, rank);CHKERRQ(ierr);

  // 6) Construct the output file name 
  char outFile[256];
  ierr = ConstructOutputFilename(field_name, ti, out_ext, outFile, sizeof(outFile));CHKERRQ(ierr);

  // 7) Create the VTK file 
  if (!rank) {
    LOG_ALLOW(GLOBAL, LOG_INFO, "main - Creating VTK file '%s'.\n", outFile);
  }
  PetscErrorCode errorCode = 0;
  errorCode = CreateAndWriteVTKFile(outFile, &meta, PETSC_COMM_WORLD);
  if (errorCode) {
    PetscPrintf(PETSC_COMM_WORLD,
                "[ERROR] CreateVTKFileFromMetadata returned %d.\n", errorCode);
    LOG_ALLOW(GLOBAL, LOG_ERROR,
              "main - CreateVTKFileFromMetadata failed with code=%d.\n", errorCode);
    goto finalize;
  }
  if (!rank) {
    PetscPrintf(PETSC_COMM_SELF, "[postprocess] Wrote file: %s\n", outFile);
    LOG_ALLOW(GLOBAL, LOG_INFO, "main - Successfully wrote file: %s\n", outFile);
  }

finalize:
  // Cleanup rank-0 memory 
  if (!rank) {
    free(coordsArray);
    free(scalarArray);
    if (meta.fileType == VTK_POLYDATA) {
      free(meta.connectivity);
      free(meta.offsets);
    }
  }

  // Finalize PETSc //
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "main - Finalizing PETSc.\n");
  PetscFinalize();
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "main - Program completed.\n");
  return 0;
}
*/ 
