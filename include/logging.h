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

// --------------------- Logging Levels Definition ---------------------

/**
 * @brief Enumeration of logging levels.
 *
 * Defines various severity levels for logging messages.
 */
typedef enum {
    LOG_ERROR = 0,   /**< Critical errors that may halt the program */
    LOG_WARNING,     /**< Non-critical issues that warrant attention */
    LOG_INFO,        /**< Informational messages about program execution */
    LOG_DEBUG        /**< Detailed debugging information */
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

/* ------------------------------------------------------------------- */
/*  Per-function allow-list for logging                                */
/* ------------------------------------------------------------------- */

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
 * @brief Logging macro that checks both the log level and whether the calling function
 *        is in the allowed-function list, using synchronized output across ranks.
 *
 * This macro uses `PetscSynchronizedPrintf` and `PetscSynchronizedFlush` to ensure
 * messages from different ranks are printed in a rank-ordered fashion. It also filters
 * out messages from functions not present in the allow-list (set via `set_allowed_functions()`).
 *
 * @param scope Specifies the logging scope:
 *              - LOCAL:  Logs on the current process using MPI_COMM_SELF.
 *              - GLOBAL: Logs on all processes using MPI_COMM_WORLD.
 * @param level The severity level of the message (e.g., LOG_INFO, LOG_ERROR).
 * @param fmt   The format string for the message (similar to printf).
 * @param ...   Additional arguments for the format string (optional).
 *
 * Example usage:
 *     LOG_ALLOW_SYNC(LOCAL,  LOG_DEBUG, "Debug info: rank = %d\n", rank);
 *     LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,  "Synchronized info in %s\n", __func__);
 *
 * Notes:
 * - A function only produces output if it appears in the allow-list (`is_function_allowed(__func__)`).
 * - The log level is filtered based on the value returned by `get_log_level()`.
 * - If the function is allowed, printing is rank-synchronized (i.e., each rank’s output appears in order).
 */
#define LOG_ALLOW_SYNC(level, fmt, ...) \
    do { \
        /* Check both the logging level and the function allow-list */ \
        if ((int)(level) <= (int)get_log_level() && is_function_allowed(__func__)) { \
            /* Synchronized print (collective) on the specified communicator */ \
            PetscSynchronizedPrintf(MPI_COMM_WORLD, "[%s] " fmt, __func__, ##__VA_ARGS__); \
            /* Flush to ensure messages appear in rank order (0,1,2,...) */ \
            PetscSynchronizedFlush(MPI_COMM_WORLD, PETSC_STDOUT); \
        } \
    } while (0)


#endif // LOGGING_H
