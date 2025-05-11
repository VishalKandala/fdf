#ifndef METRIC_H
#define METRIC_H

// Include necessary headers
#include <petsc.h>
#include "common.h"   // Common type definitions
#include "logging.h"  // Logging macros and definitions
#include <stdlib.h>
#include "io.h"
#include "setup.h"    // For SetDMDAProcLayout

PetscErrorCode MetricLogicalToPhysical(UserCtx *user, const Cmpnts ***X,
                                       PetscInt i,PetscInt j,PetscInt k,
                                       PetscReal xi,PetscReal eta,PetscReal zta,
                                       Cmpnts *Xp);

PetscErrorCode MetricGetCellVertices(UserCtx *user,const Cmpnts ***X,
                                     PetscInt i,PetscInt j,PetscInt k,
                                     Cmpnts V[8]);

PetscErrorCode MetricJacobian(UserCtx *user,const Cmpnts ***X,
                              PetscInt i,PetscInt j,PetscInt k,
                              PetscReal xi,PetscReal eta,PetscReal zta,
                              PetscReal J[3][3],PetscReal *detJ);

PetscErrorCode MetricVelocityContravariant(const PetscReal J[3][3],
                                           PetscReal detJ,
                                           const PetscReal u[3],PetscReal uc[3]);

#endif /* METRIC_H */
