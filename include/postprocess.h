#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include "io.h"
#include "common.h"
#include "logging.h" 
/* --------------------------------------------------------------------
   postprocess.h

   This header declares the interface for the post-processing executable
   or library. Typically, you'd have a function that runs your main 
   post-processing routine, or you might declare other helper functions.

   Here, we declare a single function: PostprocessMain, 
   which could be your main entry point if you want to call it from
   another place (or you might just put main() in postprocess.c).
-------------------------------------------------------------------- */
/**
 * @brief Gathers the contents of a distributed PETSc Vec into a single array on rank 0.
 *
 * @param[in]  inVec       The input (possibly distributed) Vec.
 * @param[out] N           The global size of the vector.
 * @param[out] arrayOut    On rank 0, points to the newly allocated array holding all data.
 *                         On other ranks, it is set to NULL.
 *
 * @return PetscErrorCode  Return 0 on success, nonzero on failure.
 */
PetscErrorCode VecToArrayOnRank0(Vec inVec, PetscInt *N, double **arrayOut);

//int PostprocessMain(int argc, char **argv);

#endif /* POSTPROCESS_H */
