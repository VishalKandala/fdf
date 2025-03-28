// logging.c
#include "logging.h"

/* Maximum temporary buffer size for converting numbers to strings */
#define TMP_BUF_SIZE 128

// --------------------- Static Variable for Log Level ---------------------

/**
 * @brief Static variable to cache the current logging level.
 *
 * Initialized to -1 to indicate that the log level has not been set yet.
 */
static LogLevel current_log_level = -1;

// --------------------- Static Variables for Allow-List -------------------

/**
 * @brief Global/static array of function names allowed to log.
 */
static char** gAllowedFunctions = NULL;

/**
 * @brief Number of entries in the gAllowedFunctions array.
 */
static int gNumAllowed = 0;




//-----------------------------  Logging Events definition (PetscLog) ---------------------

PetscLogEvent EVENT_Individualwalkingsearch = 0 ; //Individual walking search in (walkingsearch.c/LocateParticleInGrid()
PetscLogEvent EVENT_walkingsearch = 0 ; // Total walking search in (ParticleSwarm.c/LocateAllParticlesInGrid()

// --------------------- Function Implementations ---------------------

/**
 * @brief Retrieves the current logging level from the environment variable `LOG_LEVEL`.
 *
 * The function checks the `LOG_LEVEL` environment variable and sets the logging level accordingly.
 * Supported levels are "DEBUG", "INFO", "WARNING", and defaults to "ERROR" if not set or unrecognized.
 * The log level is cached after the first call to avoid repeated environment variable checks.
 *
 * @return LogLevel The current logging level.
 */
LogLevel get_log_level() {
    if (current_log_level == -1) { // Log level not set yet
        const char *env = getenv("LOG_LEVEL");
        if (!env) {
            current_log_level = LOG_ERROR; // Default level
        }
        else if (strcmp(env, "DEBUG") == 0) {
            current_log_level = LOG_DEBUG;
        }
        else if (strcmp(env, "INFO") == 0) {
            current_log_level = LOG_INFO;
        }
        else if (strcmp(env, "WARNING") == 0) {
            current_log_level = LOG_WARNING;
        }
        else if (strcmp(env, "PROFILE") == 0) {  // <-- New profile level
            current_log_level = LOG_PROFILE;
        }
        else {
            current_log_level = LOG_ERROR; // Default if unrecognized
        }
    }
    return current_log_level;
}

/**
 * @brief Prints the current logging level to the console.
 *
 * This function retrieves the log level using `get_log_level()` and prints 
 * the corresponding log level name. It helps verify the logging configuration 
 * at runtime.
 *
 * The log levels supported are:
 * - `LOG_PROFILE` (0) : Logs performance profiling details.
 * - `LOG_ERROR`   (1) : Logs only critical errors.
 * - `LOG_WARNING` (2) : Logs warnings and errors.
 * - `LOG_INFO`    (3) : Logs general information, warnings, and errors.
 * - `LOG_DEBUG`   (4) : Logs debugging information, info, warnings, and errors.
 *
 * @note The log level is determined from the `LOG_LEVEL` environment variable.
 * If `LOG_LEVEL` is not set, it defaults to `LOG_INFO`.
 *
 * @see get_log_level()
 */
void print_log_level() {
    int level = get_log_level();
    const char *level_name = (level == LOG_ERROR)   ? "ERROR" :
                             (level == LOG_WARNING) ? "WARNING" :
                             (level == LOG_INFO)    ? "INFO" :
                             (level == LOG_DEBUG)   ? "DEBUG" :
                             (level == LOG_PROFILE) ? "PROFILE" : "UNKNOWN";
    
    LOG(GLOBAL, LOG_INFO, "Current log level: %s (%d)\n", level_name, level);
}




/**
 * @brief Sets the global list of function names that are allowed to log.
 *
 * @param functionList An array of function name strings (e.g., {"foo", "bar"}).
 * @param count The number of entries in the array.
 *
 * The existing allow-list is cleared and replaced by the new one.
 * If you pass an empty list (count = 0), then no function will be allowed
 * unless you change it later.
 */
void set_allowed_functions(const char** functionList, int count)
{
    // 1. Free any existing entries
    if (gAllowedFunctions) {
        for (int i = 0; i < gNumAllowed; ++i) {
            free(gAllowedFunctions[i]); // each was strdup'ed
        }
        free(gAllowedFunctions);
        gAllowedFunctions = NULL;
        gNumAllowed = 0;
    }

    // 2. Allocate new array
    if (count > 0) {
        gAllowedFunctions = (char**)malloc(sizeof(char*) * count);
    }

    // 3. Copy the new entries
    for (int i = 0; i < count; ++i) {
        // strdup is a POSIX function. If not available, implement your own string copy.
        gAllowedFunctions[i] = strdup(functionList[i]);
    }
    gNumAllowed = count;
}

/**
 * @brief Checks if the given function name is in the allow-list.
 *
 * @param functionName The name of the function to check.
 * @return PETSC_TRUE if functionName is allowed, otherwise PETSC_FALSE.
 *
 * If no functions are in the list, nothing is allowed by default.
 * You can reverse this logic if you prefer to allow everything
 * unless specified otherwise.
 */
PetscBool is_function_allowed(const char* functionName)
{
    // If no allow-list entries, default to disallow all
    for (int i = 0; i < gNumAllowed; ++i) {
        if (strcmp(gAllowedFunctions[i], functionName) == 0) {
            return PETSC_TRUE;
        }
    }
    return PETSC_FALSE;
}

/**
 * @brief Prints the coordinates of a cell's vertices.
 *
 * This function iterates through the eight vertices of a given cell and prints their
 * coordinates. It is primarily used for debugging purposes to verify the correctness
 * of cell vertex assignments.
 *
 * @param[in]  cell     Pointer to a `Cell` structure representing the cell, containing its vertices.
 * @param[in]  rank     MPI rank for identification (useful in parallel environments).
 * @return PetscErrorCode Returns 0 to indicate successful execution. Non-zero on failure.
 *
 * @note
 * - Ensure that the `cell` pointer is not `NULL` before calling this function..
 */
PetscErrorCode LOG_CELL_VERTICES(const Cell *cell, PetscMPIInt rank)
{

    // Validate input pointers
    if (cell == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "LOG_CELL_VERTICES - 'cell' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "LOG_CELL_VERTICES - Input parameter 'cell' is NULL.");
    }

      LOG_ALLOW(LOCAL,LOG_DEBUG, "LOG_CELL_VERTICES - Rank %d, Cell Vertices:\n", rank);
        for(int i = 0; i < 8; i++){ 
	  LOG_ALLOW(LOCAL,LOG_DEBUG, "  Vertex[%d]: (%.2f, %.2f, %.2f)\n", 
                       i, cell->vertices[i].x, cell->vertices[i].y, cell->vertices[i].z);
        }

    return 0; // Indicate successful execution
}


/**
 * @brief Prints the signed distances to each face of the cell.
 *
 * This function iterates through the six signed distances from a point to each face of a given cell
 * and prints their values. It is primarily used for debugging purposes to verify the correctness
 * of distance calculations.
 *
 * @param[in]  d        An array of six `PetscReal` values representing the signed distances.
 *                      The indices correspond to:
 *                      - d[LEFT]: Left Face
 *                      - d[RIGHT]: Right Face
 *                      - d[BOTTOM]: Bottom Face
 *                      - d[TOP]: Top Face
 *                      - d[FRONT]: Front Face
 *                      - d[BACK]: Back Face
 *
 * @return PetscErrorCode Returns 0 to indicate successful execution. Non-zero on failure.
 *
 * @note
 * - Ensure that the `d` array is correctly populated with signed distances before calling this function.
 */
PetscErrorCode LOG_FACE_DISTANCES(PetscReal* d)
{

    // Validate input array
    if (d == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, " LOG_FACE_DISTANCES - 'd' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, " LOG_FACE_DISTANCES - Input array 'd' is NULL.");
    }

      LOG_ALLOW(LOCAL,LOG_DEBUG, " LOG_FACE_DISTANCES - Face Distances:\n");
      LOG_ALLOW(LOCAL,LOG_DEBUG, "  LEFT(%d):   %.15f\n", LEFT, d[LEFT]);
      LOG_ALLOW(LOCAL,LOG_DEBUG, "  RIGHT(%d):  %.15f\n", RIGHT, d[RIGHT]);
      LOG_ALLOW(LOCAL,LOG_DEBUG, "  BOTTOM(%d): %.15f\n", BOTTOM, d[BOTTOM]);
      LOG_ALLOW(LOCAL,LOG_DEBUG, "  TOP(%d):    %.15f\n", TOP, d[TOP]);
      LOG_ALLOW(LOCAL,LOG_DEBUG, "  FRONT(%d):  %.15f\n", FRONT, d[FRONT]);
      LOG_ALLOW(LOCAL,LOG_DEBUG, "  BACK(%d):   %.15f\n", BACK, d[BACK]);

    return 0; // Indicate successful execution
}

/*
 * Helper function: Converts an integer (of type int) to a string.
 */
static void IntToStr(int value, char *buf, size_t bufsize)
{
    snprintf(buf, bufsize, "%d", value);
}

/*
 * Helper function: Converts a 64‐bit integer to a string.
 */
static void Int64ToStr(PetscInt64 value, char *buf, size_t bufsize)
{
    snprintf(buf, bufsize, "%ld", value);
}

/*
 * Helper function: Converts three integers into a formatted string "(i, j, k)".
 */
static void CellToStr(const PetscInt64 *cell, char *buf, size_t bufsize)
{
    snprintf(buf, bufsize, "(%ld, %ld, %ld)", cell[0], cell[1], cell[2]);
}

/*
 * Helper function: Converts three PetscReal values into a formatted string "(x, y, z)".
 */
static void TripleRealToStr(const PetscReal *arr, char *buf, size_t bufsize)
{
    snprintf(buf, bufsize, "(%.4f, %.4f, %.4f)", arr[0], arr[1], arr[2]);
}

/*
 * Helper function: Computes the maximum string length for each column (across all particles).
 *
 * The function examines every particle (from 0 to nParticles-1) and converts the value to a
 * string using the helper functions above. The maximum length is stored in the pointers provided.
 *
 * @param nParticles       Number of particles.
 * @param ranks            Array of particle MPI ranks.
 * @param pids             Array of particle IDs.
 * @param cellIDs          Array of cell IDs (stored consecutively, 3 per particle).
 * @param positions        Array of positions (3 per particle).
 * @param velocities       Array of velocities (3 per particle).
 * @param weights          Array of weights (3 per particle).
 * @param wRank            [out] Maximum width for Rank column.
 * @param wPID             [out] Maximum width for PID column.
 * @param wCell            [out] Maximum width for Cell column.
 * @param wPos             [out] Maximum width for Position column.
 * @param wVel             [out] Maximum width for Velocity column.
 * @param wWt              [out] Maximum width for Weights column.
 */
static PetscErrorCode ComputeMaxColumnWidths(PetscInt nParticles,
                                               const PetscMPIInt *ranks,
                                               const PetscInt64  *pids,
                                               const PetscInt64  *cellIDs,
                                               const PetscReal   *positions,
                                               const PetscReal   *velocities,
                                               const PetscReal   *weights,
                                               int *wRank, int *wPID, int *wCell,
                                               int *wPos, int *wVel, int *wWt)
{
    char tmp[TMP_BUF_SIZE];

    *wRank = strlen("Rank");  /* Start with the header label lengths */
    *wPID  = strlen("PID");
    *wCell = strlen("Cell (i,j,k)");
    *wPos  = strlen("Position (x,y,z)");
    *wVel  = strlen("Velocity (x,y,z)");
    *wWt   = strlen("Weights (a1,a2,a3)");

    for (PetscInt i = 0; i < nParticles; i++) {
        /* Rank */
        IntToStr(ranks[i], tmp, TMP_BUF_SIZE);
        if ((int)strlen(tmp) > *wRank) *wRank = (int)strlen(tmp);

        /* PID */
        Int64ToStr(pids[i], tmp, TMP_BUF_SIZE);
        if ((int)strlen(tmp) > *wPID) *wPID = (int)strlen(tmp);

        /* Cell: use the three consecutive values */
        CellToStr(&cellIDs[3 * i], tmp, TMP_BUF_SIZE);
        if ((int)strlen(tmp) > *wCell) *wCell = (int)strlen(tmp);

        /* Position */
        TripleRealToStr(&positions[3 * i], tmp, TMP_BUF_SIZE);
        if ((int)strlen(tmp) > *wPos) *wPos = (int)strlen(tmp);

        /* Velocity */
        TripleRealToStr(&velocities[3 * i], tmp, TMP_BUF_SIZE);
        if ((int)strlen(tmp) > *wVel) *wVel = (int)strlen(tmp);

        /* Weights */
        TripleRealToStr(&weights[3 * i], tmp, TMP_BUF_SIZE);
        if ((int)strlen(tmp) > *wWt) *wWt = (int)strlen(tmp);
    }
    return 0;
}

/*
 * Helper function: Builds a format string for a table row.
 *
 * The format string will include proper width specifiers for each column.
 * For example, it might create something like:
 *
 * "| %-6s | %-8s | %-20s | %-25s | %-25s | %-25s |\n"
 *
 * @param wRank   Maximum width for the Rank column.
 * @param wPID    Maximum width for the PID column.
 * @param wCell   Maximum width for the Cell column.
 * @param wPos    Maximum width for the Position column.
 * @param wVel    Maximum width for the Velocity column.
 * @param wWt     Maximum width for the Weights column.
 * @param fmtStr  Buffer in which to build the format string.
 * @param bufSize Size of fmtStr.
 */
static void BuildRowFormatString(PetscMPIInt wRank, PetscInt wPID, PetscInt wCell, PetscInt wPos, PetscInt wVel, PetscInt wWt, char *fmtStr, size_t bufSize)
{
    // Build a format string using snprintf.
    // We assume that the Rank is an int (%d), PID is a 64-bit int (%ld)
    // and the remaining columns are strings (which have been formatted already).
    snprintf(fmtStr, bufSize,
             "| %%-%dd | %%-%ldld | %%-%lds | %%-%lds | %%-%lds | %%-%lds |\n",
             wRank, wPID, wCell, wPos, wVel, wWt);
}

/*
 * Helper function: Builds a header string for the table using column titles.
 */
static void BuildHeaderString(char *headerStr, size_t bufSize, PetscMPIInt wRank, PetscInt wPID, PetscInt wCell, PetscInt wPos, PetscInt wVel, PetscInt wWt)
{
    snprintf(headerStr, bufSize,
             "| %-*s | %-*s | %-*s | %-*s | %-*s | %-*s |\n",
             (int)wRank, "Rank",
             (int)wPID, "PID",
             (int)wCell, "Cell (i,j,k)",
             (int)wPos, "Position (x,y,z)",
             (int)wVel, "Velocity (x,y,z)",
             (int)wWt, "Weights (a1,a2,a3)");
}

/**
 * @brief Prints particle fields in a table that automatically adjusts its column widths.
 *
 * This function retrieves data from the particle swarm and prints a table where the
 * width of each column is determined by the maximum width needed to display the data.
 * Only every 'printInterval'-th particle is printed.
 *
 * @param[in] user           Pointer to the UserCtx structure.
 * @param[in] printInterval  Only every printInterval‑th particle is printed.
 *
 * @return PetscErrorCode Returns 0 on success.
 */
PetscErrorCode LOG_PARTICLE_FIELDS(UserCtx* user, PetscInt printInterval)
{
    DM swarm = user->swarm;
    PetscErrorCode ierr;
    PetscInt localNumParticles;
    PetscReal *positions = NULL;
    PetscInt64 *particleIDs = NULL;
    PetscMPIInt *particleRanks = NULL;
    PetscInt64 *cellIDs = NULL;
    PetscReal *weights = NULL;
    PetscReal *velocities = NULL;
    PetscMPIInt rank;
    
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_INFO, "PrintParticleFields - Rank %d is retrieving particle data.\n", rank);

    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"PrintParticleFields - Rank %d has %ld particles.\n", rank, localNumParticles);

    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&particleIDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&particleRanks); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);

    /* Compute maximum column widths. */
    int wRank, wPID, wCell, wPos, wVel, wWt;
    wRank = wPID = wCell = wPos = wVel = wWt = 0;
    ierr = ComputeMaxColumnWidths(localNumParticles, particleRanks, particleIDs, cellIDs,
                                  positions, velocities, weights,
                                  &wRank, &wPID, &wCell, &wPos, &wVel, &wWt); CHKERRQ(ierr);

    /* Build a header string and a row format string. */
    char headerFmt[256];
    char rowFmt[256];
    BuildHeaderString(headerFmt, sizeof(headerFmt), wRank, wPID, wCell, wPos, wVel, wWt);
    BuildRowFormatString(wRank, wPID, wCell, wPos, wVel, wWt, rowFmt, sizeof(rowFmt));

    /* Print header (using synchronized printing for parallel output). */
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------------------------------------------\n"); CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s", headerFmt); CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------------------------------------------\n"); CHKERRQ(ierr);

    /* Loop over particles and print every printInterval-th row. */
    char rowStr[256];
    for (PetscInt i = 0; i < localNumParticles; i++) {
        if (i % printInterval == 0) {
            /* Format the row by converting each field to a string first.
             * We use temporary buffers and then build the row string.
             */
            char cellStr[TMP_BUF_SIZE], posStr[TMP_BUF_SIZE], velStr[TMP_BUF_SIZE], wtStr[TMP_BUF_SIZE];
            CellToStr(&cellIDs[3*i], cellStr, TMP_BUF_SIZE);
            TripleRealToStr(&positions[3*i], posStr, TMP_BUF_SIZE);
            TripleRealToStr(&velocities[3*i], velStr, TMP_BUF_SIZE);
            TripleRealToStr(&weights[3*i], wtStr, TMP_BUF_SIZE);
            
            /* Build the row string. Note that for the integer fields we can use the row format string. */
            snprintf(rowStr, sizeof(rowStr), rowFmt,
                     particleRanks[i],
                     particleIDs[i],
                     cellStr,
                     posStr,
                     velStr,
                     wtStr);
            ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%s", rowStr); CHKERRQ(ierr);
        }
    }
    
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "--------------------------------------------------------------------------------------------------------------\n"); CHKERRQ(ierr);
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n"); CHKERRQ(ierr);
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);

    LOG_ALLOW_SYNC(GLOBAL,LOG_DEBUG,"PrintParticleFields - Completed printing on Rank %d.\n", rank);

    /* Restore fields */
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&particleIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_rank", NULL, NULL, (void**)&particleRanks); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL,LOG_DEBUG, "PrintParticleFields - Restored all particle fields.\n");
    return 0;
}
 
 /**
  * @brief Logs the interpolation error between the analytical and computed solutions.
  *
  * This function creates global vectors for the "position" and "velocity" fields from the DMSwarm,
  * applies the analytical solution to the position vector, and then computes the L2 norm of the difference
  * between the analytical and computed solutions. The resulting interpolation error is logged.
  *
  * @param user The user context.
  * @return PetscErrorCode Returns 0 on success.
  */
 PetscErrorCode LOG_INTERPOLATION_ERROR(UserCtx *user)
 {
     PetscErrorCode ierr;
     DM swarm = user->swarm;
     Vec positionVec,analyticalvelocityVec,velocityVec, errorVec;
     PetscReal Interpolation_error = 0.0;
     PetscReal AnalyticalSolution_magnitude = 0.0;
     PetscReal ErrorPercentage = 0.0;
 
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "LOG_INTERPOLATION_ERROR - Creating global vectors for 'position' and 'velocity'.\n");
     ierr = DMSwarmCreateGlobalVectorFromField(swarm, "position", &positionVec); CHKERRQ(ierr);
     ierr = DMSwarmCreateGlobalVectorFromField(swarm, "velocity", &velocityVec); CHKERRQ(ierr);
 
     ierr = VecDuplicate(positionVec, &analyticalvelocityVec); CHKERRQ(ierr);
     ierr = VecCopy(positionVec, analyticalvelocityVec); CHKERRQ(ierr);

     LOG_ALLOW(GLOBAL, LOG_DEBUG, "LOG_INTERPOLATION_ERROR - Applying analytical solution to position vector.\n");
     ierr = SetAnalyticalSolution(analyticalvelocityVec); CHKERRQ(ierr);

     ierr = VecDuplicate(analyticalvelocityVec, &errorVec); CHKERRQ(ierr);
     ierr = VecCopy(analyticalvelocityVec, errorVec); CHKERRQ(ierr);
     
     // Get the magnitude of analytical solution
     ierr = VecNorm(analyticalvelocityVec,NORM_2,&AnalyticalSolution_magnitude); CHKERRQ(ierr);
 
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "LOG_INTERPOLATION_ERROR - Computing difference between velocity and analytical solution.\n");
     ierr = VecAXPY(errorVec, -1.0, velocityVec); CHKERRQ(ierr);
 
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "LOG_INTERPOLATION_ERROR - Calculating L2 norm of the difference.\n");
     ierr = VecNorm(errorVec, NORM_2, &Interpolation_error); CHKERRQ(ierr);

     ErrorPercentage = Interpolation_error/AnalyticalSolution_magnitude;
     
     ErrorPercentage = ErrorPercentage*100;
     
     LOG_ALLOW(GLOBAL, LOG_INFO, "LOG_INTERPOLATION_ERROR - Interpolation error (%%): %g\n", ErrorPercentage);
 
     PetscPrintf(PETSC_COMM_WORLD, "Interpolation error (%%): %g\n", ErrorPercentage);

     ierr = VecDestroy(&errorVec); CHKERRQ(ierr);
     ierr = DMSwarmDestroyGlobalVectorFromField(swarm, "position", &positionVec); CHKERRQ(ierr);
     ierr = DMSwarmDestroyGlobalVectorFromField(swarm, "velocity", &velocityVec); CHKERRQ(ierr);
 
     LOG_ALLOW(GLOBAL, LOG_DEBUG, "LOG_INTERPOLATION_ERROR - Completed logging interpolation error.\n");
     return 0;
 }
 
