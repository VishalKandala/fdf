/**
 * @file logging.h
 * @brief Logging utilities and macros for PETSc-based applications.
 *
 * This header defines logging levels, scopes, and macros for consistent logging throughout the application.
 * It provides functions to retrieve the current logging level and macros to simplify logging with scope control.
 */

#ifndef LOGGING_H
#define LOGGING_H

// Include necessary headers
#include <petsc.h>   // PETSc library header
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "AnalyticalSolution.h"
// --------------------- Logging Levels Definition ---------------------

/**
 * @brief Enumeration of logging levels.
 *
 * Defines various severity levels for logging messages.
 */

typedef enum {
    LOG_ERROR = 0,   /**< Critical errors that may halt the program */
    LOG_WARNING,     /**< Non-critical issues that warrant attention */
    LOG_PROFILE,      /**< Exclusive log level for performance timing and profiling */
    LOG_INFO,        /**< Informational messages about program execution */
    LOG_DEBUG,       /**< Detailed debugging information */
} LogLevel;

// -------------------- Logging Scope Definitions ------------------

/**
 * @brief Logging scope definitions for controlling message output.
 *
 * - LOCAL: Logs on the current process using MPI_COMM_SELF.
 * - GLOBAL: Logs across all processes using MPI_COMM_WORLD.
 */
#define LOCAL  0  ///< Scope for local logging on the current process.
#define GLOBAL 1  ///< Scope for global logging across all processes.


// ---------------------Logging Events declarations ----------------

extern PetscLogEvent EVENT_Individualwalkingsearch;
extern PetscLogEvent EVENT_walkingsearch;


// --------------------- Logging Macros ---------------------

/**
 * @brief Logging macro for PETSc-based applications with scope control.
 *
 * This macro provides a convenient way to log messages with different scopes
 * (LOCAL or GLOBAL) and severity levels. It utilizes PETSc's `PetscPrintf` 
 * function for message output.
 *
 * @param scope Specifies the logging scope:
 *              - LOCAL: Logs on the current process using MPI_COMM_SELF.
 *              - GLOBAL: Logs on all processes using MPI_COMM_WORLD.
 * @param level The severity level of the message (e.g., LOG_INFO, LOG_ERROR).
 * @param fmt   The format string for the message (similar to printf).
 * @param ...   Additional arguments for the format string (optional).
 *
 * Example usage:
 *     LOG(LOCAL, LOG_ERROR, "An error occurred at index %d.\n", idx);
 *     LOG(GLOBAL, LOG_INFO, "Grid size: %d x %d x %d.\n", nx, ny, nz);
 */
#define LOG(scope, level, fmt, ...) \
    do { \
        /* Determine the MPI communicator based on the scope */ \
        MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
        /* Check if the log level is within the allowed range */ \
        if ((int)(level) <= (int)get_log_level()) { \
            /* Print the message to the specified communicator */ \
            PetscPrintf(comm, fmt, ##__VA_ARGS__); \
        } \
    } while (0)

/**
 * @brief Default logging macro for PETSc-based applications.
 *
 * This macro simplifies logging by defaulting the scope to GLOBAL
 * (i.e., `MPI_COMM_WORLD`) and providing a convenient interface for
 * common logging needs.
 *
 * @param level The severity level of the message (e.g., LOG_ERROR, LOG_INFO).
 * @param fmt   The format string for the log message (similar to printf).
 * @param ...   Additional arguments for the format string (optional).
 *
 * Example usage:
 *     LOG_DEFAULT(LOG_ERROR, "Error occurred at index %d.\n", idx);
 *     LOG_DEFAULT(LOG_INFO, "Grid size: %d x %d x %d.\n", nx, ny, nz);
 *
 * @note
 * - By default, this macro logs across all MPI processes using `MPI_COMM_WORLD`.
 * - If finer control (e.g., local logging) is required, use the more general `LOG` macro.
 * - The log level is filtered based on the value returned by `get_log_level()`.
 */
#define LOG_DEFAULT(level, fmt, ...) \
    do { \
        /* Set the communicator to global (MPI_COMM_WORLD) by default */ \
        MPI_Comm comm = MPI_COMM_WORLD; \
        /* Check if the log level is within the allowed range */ \
        if ((int)(level) <= (int)get_log_level()) { \
            /* Print the message using PetscPrintf with the global communicator */ \
            PetscPrintf(comm, fmt, ##__VA_ARGS__); \
        } \
    } while (0)

/**
 * @brief Logging macro for PETSc-based applications with scope control, 
 *        using synchronized output across processes.
 *
 * This macro uses `PetscSynchronizedPrintf` and `PetscSynchronizedFlush` to 
 * ensure messages from different ranks are printed in a synchronized (rank-by-rank) 
 * manner, preventing interleaved outputs.
 *
 * @param scope Specifies the logging scope:
 *              - LOCAL:  Logs on the current process using MPI_COMM_SELF.
 *              - GLOBAL: Logs on all processes using MPI_COMM_WORLD.
 * @param level The severity level of the message (e.g., LOG_INFO, LOG_ERROR).
 * @param fmt   The format string for the message (similar to printf).
 * @param ...   Additional arguments for the format string (optional).
 *
 * Example usage:
 *     LOG_SYNC(LOCAL, LOG_ERROR, "An error occurred at index %d.\n", idx);
 *     LOG_SYNC(GLOBAL, LOG_INFO, "Synchronized info: rank = %d.\n", rank);
 */
#define LOG_SYNC(scope, level, fmt, ...) \
    do { \
        /* Determine the MPI communicator based on the scope */ \
        MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
        /* Check if the log level is within the allowed range */ \
        if ((int)(level) <= (int)get_log_level()) { \
            /* Synchronized print (collective) on the specified communicator */ \
            PetscSynchronizedPrintf(comm, fmt, ##__VA_ARGS__); \
            /* Ensure all ranks have finished printing before continuing */ \
            PetscSynchronizedFlush(comm, PETSC_STDOUT); \
        } \
    } while (0)

/**
 * @brief Default synchronized logging macro for PETSc-based applications.
 *
 * This macro simplifies logging by defaulting the scope to GLOBAL 
 * (i.e., `MPI_COMM_WORLD`) and provides synchronized output across 
 * all processes.
 *
 * @param level The severity level of the message (e.g., LOG_ERROR, LOG_INFO).
 * @param fmt   The format string for the log message (similar to printf).
 * @param ...   Additional arguments for the format string (optional).
 *
 * Example usage:
 *     LOG_SYNC_DEFAULT(LOG_ERROR, "Error at index %d.\n", idx);
 *     LOG_SYNC_DEFAULT(LOG_INFO,  "Process rank: %d.\n", rank);
 *
 * @note
 * - By default, this macro logs across all MPI processes using `MPI_COMM_WORLD`.
 * - If local (per-process) logging is required, use the more general `LOG_SYNC` macro.
 * - The log level is filtered based on the value returned by `get_log_level()`.
 */
#define LOG_SYNC_DEFAULT(level, fmt, ...) \
    do { \
        if ((int)(level) <= (int)get_log_level()) { \
            PetscSynchronizedPrintf(MPI_COMM_WORLD, fmt, ##__VA_ARGS__); \
            PetscSynchronizedFlush(MPI_COMM_WORLD, PETSC_STDOUT); \
        } \
    } while (0)



/**
 * @brief Logging macro that checks both the log level and whether the calling function
 *        is in the allowed-function list before printing. Useful for selective, per-function logging.
 *
 * @param scope Specifies the logging scope (LOCAL or GLOBAL).
 * @param level The severity level of the message (e.g., LOG_INFO, LOG_ERROR).
 * @param fmt   The format string for the message (similar to printf).
 * @param ...   Additional arguments for the format string (optional).
 *
 * Example usage:
 *     LOG_ALLOW(LOCAL, LOG_DEBUG, "Debugging info in function: %s\n", __func__);
 */
#define LOG_ALLOW(scope, level, fmt, ...) \
    do { \
        MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
        if ((int)(level) <= (int)get_log_level() && is_function_allowed(__func__)) { \
            PetscPrintf(comm, "[%s] " fmt, __func__, ##__VA_ARGS__); \
        } \
    } while (0)


/**
 * @brief Synchronized logging macro that checks both the log level 
 *        and whether the calling function is in the allow-list.
 *
 * This macro uses `PetscSynchronizedPrintf` and `PetscSynchronizedFlush` to
 * ensure messages from different ranks are printed in a rank-ordered fashion
 * (i.e., to avoid interleaving). It also filters out messages if the current
 * function is not in the allow-list (`is_function_allowed(__func__)`) or the
 * requested log level is higher than `get_log_level()`.
 *
 * @param scope  Either LOCAL (MPI_COMM_SELF) or GLOBAL (MPI_COMM_WORLD).
 * @param level  One of LOG_ERROR, LOG_WARNING, LOG_INFO, LOG_DEBUG.
 * @param fmt    A `printf`-style format string (e.g., "Message: %d\n").
 * @param ...    Variadic arguments to fill in `fmt`.
 *
 * Example usage:
 *
 * \code{.c}
 * LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Debug info: rank = %d\n", rank);
 * LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,  "Synchronized info in %s\n", __func__);
 * \endcode
 */
#define LOG_ALLOW_SYNC(scope,level, fmt, ...)                                 \
    do {                                                                       \
        if ((scope != LOCAL && scope != GLOBAL)) {                             \
            fprintf(stderr, "LOG_ALLOW_SYNC ERROR: Invalid scope at %s:%d\n",  \
                    __FILE__, __LINE__);                                       \
        } else if (is_function_allowed(__func__) &&                            \
                  (int)(level) <= (int)get_log_level()) {                      \
            MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
            PetscSynchronizedPrintf(comm, "[%s] " fmt, __func__, ##__VA_ARGS__); \
            PetscSynchronizedFlush(comm, PETSC_STDOUT);                        \
        }                                                                      \
    } while (0)

/**
 * @brief Logs a message inside a loop, but only every `interval` iterations.
 *
 * @param scope     LOCAL or GLOBAL.
 * @param level     LOG_* level.
 * @param iterVar   The loop variable (e.g., i).
 * @param interval  Only log when (iterVar % interval == 0).
 * @param fmt       printf-style format string.
 * @param ...       Variadic arguments to include in the formatted message.
 *
 * Example:
 *    for (int i = 0; i < 100; i++) {
 *        LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, i, 10, "Value of i=%d\n", i);
 *    }
 */
#define LOG_LOOP_ALLOW(scope,level, iterVar, interval, fmt, ...)              \
    do {                                                                       \
        if (is_function_allowed(__func__) && (int)(level) <= (int)get_log_level()) { \
            if ((iterVar) % (interval) == 0) {                                 \
                MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
                PetscPrintf(comm, "[%s] [Iter=%d] " fmt,                       \
                            __func__, (iterVar), ##__VA_ARGS__);               \
            }                                                                  \
        }                                                                      \
    } while (0)

/**
 * @brief Logs a single element of an array, given an index.
 *
 * @param scope   Either LOCAL or GLOBAL.
 * @param level   LOG_ERROR, LOG_WARNING, LOG_INFO, or LOG_DEBUG.
 * @param arr     Pointer to the array to log from.
 * @param length  The length of the array (to prevent out-of-bounds).
 * @param idx     The index of the element to print.
 * @param fmt     The printf-style format specifier (e.g. "%g", "%f", etc.).
 *
 * This macro only logs if:
 *  1) The current function is in the allow-list (`is_function_allowed(__func__)`).
 *  2) The requested logging `level` <= the current global `get_log_level()`.
 *  3) The index `idx` is valid (0 <= idx < length).
 */
#define LOG_ARRAY_ELEMENT_ALLOW(scope,level, arr, length, idx, fmt)          \
    do {                                                                      \
        if (is_function_allowed(__func__) && (int)(level) <= (int)get_log_level()) { \
            if ((idx) >= 0 && (idx) < (length)) {                             \
                MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
                PetscPrintf(comm, "[%s] arr[%d] = " fmt "\n",                 \
                            __func__, (idx), (arr)[idx]);                     \
            }                                                                 \
        }                                                                     \
    } while (0)

/**
 * @brief Logs a consecutive subrange of an array.
 *
 * @param scope   Either LOCAL or GLOBAL.
 * @param level   LOG_ERROR, LOG_WARNING, LOG_INFO, or LOG_DEBUG.
 * @param arr     Pointer to the array to log from.
 * @param length  Total length of the array.
 * @param start   Starting index of the subrange.
 * @param end     Ending index of the subrange (inclusive).
 * @param fmt     The printf-style format specifier (e.g., "%g", "%f").
 *
 * This macro prints each element arr[i] for i in [start, end], bounded by [0, length-1].
 */
#define LOG_ARRAY_SUBRANGE_ALLOW(scope,level, arr, length, start, end, fmt)  \
    do {                                                                      \
        if (is_function_allowed(__func__) && (int)(level) <= (int)get_log_level()) { \
            MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
            PetscInt _start = (start) < 0 ? 0 : (start);                      \
            PetscInt _end   = (end) >= (length) ? (length) - 1 : (end);       \
            for (PetscInt i = _start; i <= _end; i++) {                       \
                PetscPrintf(comm, "[%s] arr[%d] = " fmt "\n", __func__, i, (arr)[i]); \
            }                                                                 \
        }                                                                     \
    } while (0)


/**
 * @brief Begins timing a function by:
 *        1. Starting a wall-clock timer (PetscTime).
 *        2. Beginning a PETSc log event.
 *
 * @param eventID   A previously registered PetscLogEvent (e.g., EVENT_EvalPosition).
 * @param scope     LOCAL or GLOBAL (for potential later usage in logging).
 * @param level     LOG_ERROR, LOG_WARNING, LOG_INFO, LOG_DEBUG.
 *
 * Example usage:
 * \code{.c}
 * LOG_FUNC_TIMER_BEGIN_EVENT(EVENT_EvalPosition, LOCAL, LOG_DEBUG);
 * // ... function body ...
 * LOG_FUNC_TIMER_END_EVENT(EVENT_EvalPosition, LOCAL, LOG_DEBUG);
 * \endcode
 */
#define LOG_FUNC_TIMER_BEGIN_EVENT(eventID, scope)			\
    double __funcTimerStart = 0.0;                                              \
    PetscBool __funcTimerActive = PETSC_FALSE;                                  \
    do {                                                                        \
        if (is_function_allowed(__func__) && (int)(LOG_PROFILE) <= (int)get_log_level()) { \
            PetscLogDouble _timeStamp = 0.0;                                    \
            PetscTime(&_timeStamp);                                             \
            __funcTimerStart = (double)_timeStamp;                              \
            __funcTimerActive = PETSC_TRUE;                                     \
            /* Start the PETSc log event (rank-wide). */                        \
            PetscLogEventBegin(eventID, 0, 0, 0, 0);                            \
        }                                                                       \
    } while (0)

/**
 * @brief Ends timing a function by:
 *        1. Ending a PETSc log event.
 *        2. Logging the wall-clock elapsed time, if active.
 *
 * @param eventID  The same PetscLogEvent handle passed to LOG_FUNC_TIMER_BEGIN_EVENT.
 * @param scope    LOCAL or GLOBAL for possible MPI_Comm usage.
 * @param level    The log level at which to print the timing message.
 *
 * The log message is only printed if the function is in the allow-list
 * and the current log level is >= the requested `level`.
 */
#define LOG_FUNC_TIMER_END_EVENT(eventID, scope)                         \
    do {                                                                        \
        if (__funcTimerActive == PETSC_TRUE) {                                  \
            /* End the PETSc log event */                                       \
            PetscLogEventEnd(eventID, 0, 0, 0, 0);                              \
            /* Log the wall-clock elapsed time */                               \
            if (is_function_allowed(__func__) && (int)(LOG_PROFILE) <= (int)get_log_level()) { \
                PetscLogDouble _timeEnd = 0.0;                                  \
                PetscTime(&_timeEnd);                                           \
                double elapsed = (double)_timeEnd - __funcTimerStart;           \
                MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
                PetscPrintf(comm, "[%s] Elapsed Time: %f seconds\n", __func__, elapsed); \
            }                                                                   \
        }                                                                       \
    } while (0)


#define LOG_PROFILE_MSG(scope, fmt, ...)                                     \
    do {                                                                     \
        if ((int)(LOG_PROFILE) <= (int)get_log_level()) {                    \
            MPI_Comm comm = (scope == LOCAL) ? MPI_COMM_SELF : MPI_COMM_WORLD; \
            PetscPrintf(comm, "[PROFILE] " fmt, ##__VA_ARGS__);              \
        }                                                                    \
    } while (0)

// --------------------- Function Declarations ---------------------

/**
 * @brief Retrieves the current logging level from the environment variable `LOG_LEVEL`.
 *
 * The function checks the `LOG_LEVEL` environment variable and sets the logging level accordingly.
 * Supported levels are "DEBUG", "INFO", "WARNING", and defaults to "ERROR" if not set or unrecognized.
 *
 * @return LogLevel The current logging level.
 */
LogLevel get_log_level();

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
void print_log_level();

/**
 * @brief Sets the global list of function names that are allowed to log.
 *
 * You can replace the entire list of allowed function names at runtime.
 */
void set_allowed_functions(const char** functionList, int count);

/**
 * @brief Checks if a given function is in the allow-list.
 *
 * This helper is used internally by the LOG_ALLOW macro.
 */
PetscBool is_function_allowed(const char* functionName);

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
PetscErrorCode LOG_CELL_VERTICES(const Cell *cell, PetscInt rank);

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
PetscErrorCode LOG_FACE_DISTANCES(PetscReal* d);

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
PetscErrorCode LOG_PARTICLE_FIELDS(UserCtx* user, PetscInt printInterval);

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
  PetscErrorCode LOG_INTERPOLATION_ERROR(UserCtx *user);

#endif // LOGGING_H
