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

/**
 * @brief Computes the primary face metric components (Csi, Eta, Zet), including
 *        boundary extrapolation, and stores them in the corresponding global Vec
 *        members of the UserCtx structure (user->Csi, user->Eta, user->Zet).
 *
 *        This is a self-contained routine that performs the following steps:
 *        1. Obtains local ghosted nodal coordinates using DMGetCoordinatesLocal.
 *        2. Calculates metrics for INTERIOR faces where finite difference stencils are valid.
 *        3. EXTRAPOLATES metrics for faces on the physical domain boundaries by copying
 *           from the nearest computed interior face.
 *        4. Assembles the global `user->Csi`, `user->Eta`, `user->Zet` Vecs.
 *        5. Updates the local ghosted `user->lCsi`, `user->lEta`, `user->lZet` Vecs.
 *
 * @param[in,out] user             Pointer to the UserCtx structure.
 *
 * @return PetscErrorCode 0 on success.
 *
 * @note
 *  - This function is a complete "compute and make ready" unit for Csi, Eta, and Zet.
 *  - It's recommended to call `VecZeroEntries` on user->Csi, Eta, Zet before this
 *    if they might contain old data.
 */
PetscErrorCode ComputeFaceMetrics(UserCtx *user);

/**
 * @brief Calculates the cell-centered inverse Jacobian determinant (1/J) for INTERIOR cells
 *        and stores it in `user->Aj`. This version includes boundary extrapolation.
 *
 *        Nodal coordinates are obtained internally.
 *        Refer to previous Doxygen comments for details on physical locations and
 *        storage convention (`aj_arr[k_n][j_n][i_n]` for cell `C(i_n-1,j_n-1,k_n-1)`).
 *
 * @param[in,out] user             Pointer to the UserCtx structure.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeCellCenteredJacobianInverse(UserCtx *user);

/**
 * @brief Ensure a **right-handed** metric basis (`Csi`, `Eta`, `Zet`) and a
 *        **positive Jacobian** (`Aj`) over the whole domain.
 *
 * The metric-generation kernels are completely algebraic, so they will happily
 * deliver a *left-handed* basis if the mesh file enumerates nodes in the
 * opposite ζ-direction.  
 * This routine makes the orientation explicit and—if needed—repairs it
 * **once per run**:
 *
 * | Step | Action |
 * |------|--------|
 * | 1 | Compute global `Aj_min`, `Aj_max`.                          |
 * | 2 | **Mixed signs** (`Aj_min < 0 && Aj_max > 0`) &rarr; abort: the mesh is topologically inconsistent. |
 * | 3 | **All negative** (`Aj_max < 0`) &rarr; flip <br>`Csi`, `Eta`, `Zet`, `Aj` & update local ghosts. |
 * | 4 | Store `user->orientation = ±1` so BC / IC routines can apply sign-aware logic if they care about inlet direction. |
 *
 * @param[in,out] user  Fully initialised #UserCtx that already contains  
 *                      `Csi`, `Eta`, `Zet`, `Aj`, their **local** ghosts, and
 *                      valid distributed DMs.
 *
 * @return `0` on success or a PETSc error code on failure.
 *
 * @note  Call **immediately after** `ComputeCellCenteredJacobianInverse()` and
 *        before any routine that differentiates or applies BCs.
 *
 * @author <your name>
 */
PetscErrorCode CheckAndFixGridOrientation(UserCtx *user);
#endif /* METRIC_H */
