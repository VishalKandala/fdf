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

int PostprocessMain(int argc, char **argv);

#endif /* POSTPROCESS_H */
