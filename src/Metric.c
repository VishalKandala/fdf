/* Metric.c ------------------------------------------------------------------
 * Utility routines for curvilinear-grid metric operations used by the
 * particle-swarm module.
 *
 *  –  Logical (xi,eta,zta) → physical (x,y,z) mapping via trilinear blend.
 *  –  Jacobian matrix and determinant for contravariant velocity conversion.
 *
 * The only data this file needs from the application is the DMDA that stores
 * vertex coordinates in the usual PETSc coordinate DM (user->da) and the
 * coordinate array type `Cmpnts` (three-component struct {x,y,z}).
 * ---------------------------------------------------------------------------*/

#include <petsc.h>
#include "Metric.h"          /* forward declarations + Cmpnts + UserCtx */

/* ------------------------------------------------------------------------- */
/**
 * @brief Extract the eight vertex coordinates of the hexahedral cell (i,j,k).
 *
 * Vertices are returned in the standard trilinear ordering: bit 0 → x-corner,
 * bit 1 → y-corner, bit 2 → z-corner.  (000 = origin of the cell, 111 = far
 * corner.)
 */
PetscErrorCode MetricGetCellVertices(UserCtx *user,
                                     const Cmpnts ***X,   /* coord array */
                                     PetscInt i,PetscInt j,PetscInt k,
                                     Cmpnts V[8])
{
  PetscFunctionBeginUser;
  for (PetscInt c = 0; c < 8; ++c) {
    PetscInt ii = i + ((c & 1) ? 1 : 0);
    PetscInt jj = j + ((c & 2) ? 1 : 0);
    PetscInt kk = k + ((c & 4) ? 1 : 0);
    LOG_LOOP_ALLOW(GLOBAL, LOG_DEBUG,i+j+k,10," ii: %d,jj:%d,kk:%d - Retrieved.\n",ii,jj,kk);
    V[c] = X[kk][jj][ii];
  }
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------------- */
/**
 * @brief Map logical coordinates to physical space using trilinear basis.
 *
 * @param[in]   V   Array of the eight vertex coordinates (MetricGetCellVertices).
 * @param[in]   xi,eta,zta  Logical coordinates in [0,1].
 * @param[out]  Xp  Physical coordinate.
 */
static inline void TrilinearBlend(const Cmpnts V[8],
                                  PetscReal xi,PetscReal eta,PetscReal zta,
                                  Cmpnts *Xp)
{
  PetscReal x=0,y=0,z=0;
  for (PetscInt c=0;c<8;++c) {
    PetscReal N = ((c&1)?xi : 1.0-xi ) *
                  ((c&2)?eta: 1.0-eta) *
                  ((c&4)?zta: 1.0-zta);
    x += N * V[c].x;
    y += N * V[c].y;
    z += N * V[c].z;
  }
  Xp->x = x; Xp->y = y; Xp->z = z;
}

/* ------------------------------------------------------------------------- */
/**
 * @brief Public wrapper: map (cell index, ξ,η,ζ) to (x,y,z).
 */
PetscErrorCode MetricLogicalToPhysical(UserCtx  *user,
                                       const Cmpnts ***X,
                                       PetscInt i,PetscInt j,PetscInt k,
                                       PetscReal xi,PetscReal eta,PetscReal zta,
                                       Cmpnts *Xp)
{
  PetscErrorCode ierr;
  Cmpnts V[8];
  PetscFunctionBeginUser;
  ierr = MetricGetCellVertices(user,X,i,j,k,V); CHKERRQ(ierr);
  TrilinearBlend(V,xi,eta,zta,Xp);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------------- */
/**
 * @brief Compute Jacobian matrix and its determinant at (xi,eta,zta).
 *
 *        J = [ x_ξ  x_η  x_ζ ]
 *            [ y_ξ  y_η  y_ζ ]
 *            [ z_ξ  z_η  z_ζ ]
 *
 * This is handy for converting physical velocities (u,v,w) into contravariant
 * components and for volume weighting.
 */
PetscErrorCode MetricJacobian(UserCtx *user,
                              const Cmpnts ***X,
                              PetscInt i,PetscInt j,PetscInt k,
                              PetscReal xi,PetscReal eta,PetscReal zta,
                              PetscReal J[3][3], PetscReal *detJ)
{
  PetscErrorCode ierr;
  Cmpnts V[8];
  PetscFunctionBeginUser;
  ierr = MetricGetCellVertices(user,X,i,j,k,V); CHKERRQ(ierr);

  /* derivatives of trilinear shape functions */
  PetscReal dN_dXi[8], dN_dEta[8], dN_dZta[8];
  for (PetscInt c=0;c<8;++c) {
    PetscReal sx = (c & 1) ?  1.0 : -1.0;
    PetscReal sy = (c & 2) ?  1.0 : -1.0;
    PetscReal sz = (c & 4) ?  1.0 : -1.0;
    dN_dXi [c] = 0.125 * sx * ( (c&2?eta:1-eta) ) * ( (c&4?zta:1-zta) );
    dN_dEta[c] = 0.125 * sy * ( (c&1?xi :1-xi ) ) * ( (c&4?zta:1-zta) );
    dN_dZta[c] = 0.125 * sz * ( (c&1?xi :1-xi ) ) * ( (c&2?eta:1-eta) );
  }

  /* assemble Jacobian */
  PetscReal x_xi=0,y_xi=0,z_xi=0,
            x_eta=0,y_eta=0,z_eta=0,
            x_zta=0,y_zta=0,z_zta=0;
  for (PetscInt c=0;c<8;++c) {
    x_xi  += dN_dXi [c]*V[c].x;  y_xi  += dN_dXi [c]*V[c].y;  z_xi  += dN_dXi [c]*V[c].z;
    x_eta += dN_dEta[c]*V[c].x;  y_eta += dN_dEta[c]*V[c].y;  z_eta += dN_dEta[c]*V[c].z;
    x_zta += dN_dZta[c]*V[c].x;  y_zta += dN_dZta[c]*V[c].y;  z_zta += dN_dZta[c]*V[c].z;
  }

  J[0][0]=x_xi;  J[0][1]=x_eta;  J[0][2]=x_zta;
  J[1][0]=y_xi;  J[1][1]=y_eta;  J[1][2]=y_zta;
  J[2][0]=z_xi;  J[2][1]=z_eta;  J[2][2]=z_zta;

  if (detJ) {
    *detJ = x_xi*(y_eta*z_zta - y_zta*z_eta)
          - x_eta*(y_xi*z_zta - y_zta*z_xi)
          + x_zta*(y_xi*z_eta - y_eta*z_xi);
  }
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------------- */
/**
 * @brief Convert physical velocity (u,v,w) to contravariant components (u^xi, u^eta, u^zta).
 */
PetscErrorCode MetricVelocityContravariant(const PetscReal J[3][3], PetscReal detJ,
                                           const PetscReal u[3], PetscReal uc[3])
{
  PetscFunctionBeginUser;
  /* contravariant basis vectors (row of adjugate(J)) divided by detJ */
  PetscReal gxi[3]  = {  J[1][1]*J[2][2]-J[1][2]*J[2][1],
                        -J[0][1]*J[2][2]+J[0][2]*J[2][1],
                         J[0][1]*J[1][2]-J[0][2]*J[1][1] };
  PetscReal geta[3] = { -J[1][0]*J[2][2]+J[1][2]*J[2][0],
                         J[0][0]*J[2][2]-J[0][2]*J[2][0],
                        -J[0][0]*J[1][2]+J[0][2]*J[1][0] };
  PetscReal gzta[3] = {  J[1][0]*J[2][1]-J[1][1]*J[2][0],
                        -J[0][0]*J[2][1]+J[0][1]*J[2][0],
                         J[0][0]*J[1][1]-J[0][1]*J[1][0] };

  PetscReal invDet = 1.0 / detJ;
  for (int d=0; d<3; ++d) { gxi[d]  *= invDet; geta[d] *= invDet; gzta[d] *= invDet; }

  uc[0] = gxi [0]*u[0] + gxi [1]*u[1] + gxi [2]*u[2];
  uc[1] = geta[0]*u[0] + geta[1]*u[1] + geta[2]*u[2];
  uc[2] = gzta[0]*u[0] + gzta[1]*u[1] + gzta[2]*u[2];
  PetscFunctionReturn(0);
}

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
PetscErrorCode ComputeFaceMetrics(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Cmpnts       ***csi_arr, ***eta_arr, ***zet_arr;
    Cmpnts       ***nodal_coords_arr;
    Vec            localCoords_from_dm;

    PetscFunctionBeginUser;
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting calculation and update for Csi, Eta, Zet.\n");

    ierr = DMDAGetLocalInfo(user->fda, &info); CHKERRQ(ierr);

    // --- 1. Get Nodal Physical Coordinates (Local Ghosted Array directly) ---
    ierr = DMGetCoordinatesLocal(user->da, &localCoords_from_dm); CHKERRQ(ierr);
    if (!localCoords_from_dm) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "DMGetCoordinatesLocal failed to return a coordinate vector. \n");
    ierr = DMDAVecGetArrayRead(user->fda, localCoords_from_dm, &nodal_coords_arr); CHKERRQ(ierr);

    // --- 2. Get arrays for output global Vecs from UserCtx ---
    ierr = DMDAVecGetArray(user->fda, user->Csi, &csi_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Eta, &eta_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Zet, &zet_arr); CHKERRQ(ierr);

    // Define owned node ranges (global indices)
    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;

    // Global domain dimensions (total number of nodes)
    PetscInt mx = info.mx;
    PetscInt my = info.my;
    PetscInt mz = info.mz;

    // --- 3. Calculate Csi, Eta, Zet for INTERIOR Stencils ---
    // Start loops from 1 if at global boundary 0 to ensure k_node-1 etc. are valid.
    PetscInt k_loop_start = (zs == 0) ? zs + 1 : zs;
    PetscInt j_loop_start = (ys == 0) ? ys + 1 : ys;
    PetscInt i_loop_start = (xs == 0) ? xs + 1 : xs;

    // Calculate Csi
    for (PetscInt k_node = k_loop_start; k_node < ze; ++k_node) {
        for (PetscInt j_node = j_loop_start; j_node < ye; ++j_node) {
            for (PetscInt i_node = xs; i_node < xe; ++i_node) {
	      
                PetscReal dx_deta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node-1][j_node][i_node].x - nodal_coords_arr[k_node][j_node-1][i_node].x - nodal_coords_arr[k_node-1][j_node-1][i_node].x);
                PetscReal dy_deta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node-1][j_node][i_node].y - nodal_coords_arr[k_node][j_node-1][i_node].y - nodal_coords_arr[k_node-1][j_node-1][i_node].y);
                PetscReal dz_deta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node-1][j_node][i_node].z - nodal_coords_arr[k_node][j_node-1][i_node].z - nodal_coords_arr[k_node-1][j_node-1][i_node].z);
                PetscReal dx_dzeta = 0.5 * (nodal_coords_arr[k_node][j_node-1][i_node].x + nodal_coords_arr[k_node][j_node][i_node].x - nodal_coords_arr[k_node-1][j_node-1][i_node].x - nodal_coords_arr[k_node-1][j_node][i_node].x);
                PetscReal dy_dzeta = 0.5 * (nodal_coords_arr[k_node][j_node-1][i_node].y + nodal_coords_arr[k_node][j_node][i_node].y - nodal_coords_arr[k_node-1][j_node-1][i_node].y - nodal_coords_arr[k_node-1][j_node][i_node].y);
                PetscReal dz_dzeta = 0.5 * (nodal_coords_arr[k_node][j_node-1][i_node].z + nodal_coords_arr[k_node][j_node][i_node].z - nodal_coords_arr[k_node-1][j_node-1][i_node].z - nodal_coords_arr[k_node-1][j_node][i_node].z);

		csi_arr[k_node][j_node][i_node].x = dy_deta * dz_dzeta - dz_deta * dy_dzeta;
                csi_arr[k_node][j_node][i_node].y = dz_deta * dx_dzeta - dx_deta * dz_dzeta;
                csi_arr[k_node][j_node][i_node].z = dx_deta * dy_dzeta - dy_deta * dx_dzeta;
            }
        }
    }

    // Calculate Eta
    for (PetscInt k_node = k_loop_start; k_node < ze; ++k_node) {
        for (PetscInt j_node = ys; j_node < ye; ++j_node) {
            for (PetscInt i_node = i_loop_start; i_node < xe; ++i_node) {
	      
                PetscReal dx_dxi = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node-1][j_node][i_node].x - nodal_coords_arr[k_node][j_node][i_node-1].x - nodal_coords_arr[k_node-1][j_node][i_node-1].x);
                PetscReal dy_dxi = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node-1][j_node][i_node].y - nodal_coords_arr[k_node][j_node][i_node-1].y - nodal_coords_arr[k_node-1][j_node][i_node-1].y);
                PetscReal dz_dxi = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node-1][j_node][i_node].z - nodal_coords_arr[k_node][j_node][i_node-1].z - nodal_coords_arr[k_node-1][j_node][i_node-1].z);
                PetscReal dx_dzeta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node][j_node][i_node-1].x - nodal_coords_arr[k_node-1][j_node][i_node].x - nodal_coords_arr[k_node-1][j_node][i_node-1].x);
                PetscReal dy_dzeta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node][j_node][i_node-1].y - nodal_coords_arr[k_node-1][j_node][i_node].y - nodal_coords_arr[k_node-1][j_node][i_node-1].y);
                PetscReal dz_dzeta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node][j_node][i_node-1].z - nodal_coords_arr[k_node-1][j_node][i_node].z - nodal_coords_arr[k_node-1][j_node][i_node-1].z);

		eta_arr[k_node][j_node][i_node].x = dy_dzeta * dz_dxi  - dz_dzeta * dy_dxi;
                eta_arr[k_node][j_node][i_node].y = dz_dzeta * dx_dxi  - dx_dzeta * dz_dxi;
                eta_arr[k_node][j_node][i_node].z = dx_dzeta * dy_dxi  - dy_dzeta * dx_dxi;
            }
        }
    }

    // Calculate Zet
    for (PetscInt k_node = zs; k_node < ze; ++k_node) {
        for (PetscInt j_node = j_loop_start; j_node < ye; ++j_node) {
            for (PetscInt i_node = i_loop_start; i_node < xe; ++i_node) {
	      
                PetscReal dx_dxi = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node][j_node-1][i_node].x - nodal_coords_arr[k_node][j_node][i_node-1].x - nodal_coords_arr[k_node][j_node-1][i_node-1].x);
                PetscReal dy_dxi = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node][j_node-1][i_node].y - nodal_coords_arr[k_node][j_node][i_node-1].y - nodal_coords_arr[k_node][j_node-1][i_node-1].y);
                PetscReal dz_dxi = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node][j_node-1][i_node].z - nodal_coords_arr[k_node][j_node][i_node-1].z - nodal_coords_arr[k_node][j_node-1][i_node-1].z);
                PetscReal dx_deta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node][j_node][i_node-1].x - nodal_coords_arr[k_node][j_node-1][i_node].x - nodal_coords_arr[k_node][j_node-1][i_node-1].x);
                PetscReal dy_deta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node][j_node][i_node-1].y - nodal_coords_arr[k_node][j_node-1][i_node].y - nodal_coords_arr[k_node][j_node-1][i_node-1].y);
                PetscReal dz_deta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node][j_node][i_node-1].z - nodal_coords_arr[k_node][j_node-1][i_node].z - nodal_coords_arr[k_node][j_node-1][i_node-1].z);

		zet_arr[k_node][j_node][i_node].x = dy_dxi * dz_deta - dz_dxi * dy_deta;
                zet_arr[k_node][j_node][i_node].y = dz_dxi * dx_deta - dx_dxi * dz_deta;
                zet_arr[k_node][j_node][i_node].z = dx_dxi * dy_deta - dy_dxi * dx_deta;
            }
        }
    }

    // --- 4. Boundary Extrapolation ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Extrapolating boundary values for Csi, Eta, Zet.\n");
    PetscInt i_bnd, j_bnd, k_bnd;

    if (xs == 0) { // If this rank owns the global i=0 boundary
        i_bnd = 0;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
                if (i_bnd + 1 < mx) {
                    eta_arr[k_bnd][j_bnd][i_bnd] = eta_arr[k_bnd][j_bnd][i_bnd+1];
                    zet_arr[k_bnd][j_bnd][i_bnd] = zet_arr[k_bnd][j_bnd][i_bnd+1];
                }
            }
        }
    }
    if (xe == mx) { // If this rank owns the global i=mx-1 boundary
        i_bnd = mx - 1;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
                if (i_bnd - 1 >= 0) {
                    eta_arr[k_bnd][j_bnd][i_bnd] = eta_arr[k_bnd][j_bnd][i_bnd-1];
                    zet_arr[k_bnd][j_bnd][i_bnd] = zet_arr[k_bnd][j_bnd][i_bnd-1];
                }
            }
        }
    }
    if (ys == 0) {
        j_bnd = 0;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                if (j_bnd + 1 < my) {
                    csi_arr[k_bnd][j_bnd][i_bnd] = csi_arr[k_bnd][j_bnd+1][i_bnd];
                    zet_arr[k_bnd][j_bnd][i_bnd] = zet_arr[k_bnd][j_bnd+1][i_bnd];
                }
            }
        }
    }
    if (ye == my) {
        j_bnd = my - 1;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                if (j_bnd - 1 >= 0) {
                    csi_arr[k_bnd][j_bnd][i_bnd] = csi_arr[k_bnd][j_bnd-1][i_bnd];
                    zet_arr[k_bnd][j_bnd][i_bnd] = zet_arr[k_bnd][j_bnd-1][i_bnd];
                }
            }
        }
    }
    if (zs == 0) {
        k_bnd = 0;
        for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                if (k_bnd + 1 < mz) {
                    csi_arr[k_bnd][j_bnd][i_bnd] = csi_arr[k_bnd+1][j_bnd][i_bnd];
                    eta_arr[k_bnd][j_bnd][i_bnd] = eta_arr[k_bnd+1][j_bnd][i_bnd];
                }
            }
        }
    }
    if (ze == mz) {
        k_bnd = mz - 1;
        for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                if (k_bnd - 1 >= 0) {
                    csi_arr[k_bnd][j_bnd][i_bnd] = csi_arr[k_bnd-1][j_bnd][i_bnd];
                    eta_arr[k_bnd][j_bnd][i_bnd] = eta_arr[k_bnd-1][j_bnd][i_bnd];
                }
            }
        }
    }

    // --- 5. Restore all arrays ---
    ierr = DMDAVecRestoreArrayRead(user->fda, localCoords_from_dm, &nodal_coords_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Csi, &csi_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Eta, &eta_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Zet, &zet_arr); CHKERRQ(ierr);

    // --- 6. Assemble Global Vectors ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Assembling global Csi, Eta, Zet.\n");
    ierr = VecAssemblyBegin(user->Csi); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->Csi); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->Eta); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->Eta); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->Zet); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->Zet); CHKERRQ(ierr);

    // --- 7. Update Local Ghosted Versions ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Updating local lCsi, lEta, lZet.\n");
    ierr = UpdateLocalGhosts(user, "Csi"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "Eta"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "Zet"); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Completed calculation, extrapolation, and update for Csi, Eta, Zet.\n");
    PetscFunctionReturn(0);
}

/**
 * @brief Calculates the cell-centered inverse Jacobian determinant (1/J), including
 *        boundary extrapolation, stores it in `user->Aj`, assembles `user->Aj`, and
 *        updates `user->lAj`.
 *
 * @param[in,out] user             Pointer to the UserCtx structure.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeCellCenteredJacobianInverse(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    PetscScalar  ***aj_arr;
    Cmpnts       ***nodal_coords_arr;
    Vec            localCoords_from_dm;

    PetscFunctionBeginUser;
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting calculation, extrapolation, and update for Aj.\n");

    // --- 1. Get Nodal Coordinates and Output Array ---
    ierr = DMGetCoordinatesLocal(user->da, &localCoords_from_dm); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, localCoords_from_dm, &nodal_coords_arr); CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->Aj, &aj_arr); CHKERRQ(ierr);

    // Define owned node ranges (global indices)
    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;

    // Global domain dimensions (total number of nodes)
    PetscInt mx = info.mx;
    PetscInt my = info.my;
    PetscInt mz = info.mz;

    // --- 2. Calculate Aj for INTERIOR Stencils ---
    // Start loops from 1 if at global boundary 0 because stencil accesses i_node-1 etc.
    PetscInt k_start_node = (zs == 0) ? zs + 1 : zs;
    PetscInt j_start_node = (ys == 0) ? ys + 1 : ys;
    PetscInt i_start_node = (xs == 0) ? xs + 1 : xs;

    for (PetscInt k_node = k_start_node; k_node < ze; ++k_node) {
        for (PetscInt j_node = j_start_node; j_node < ye; ++j_node) {
            for (PetscInt i_node = i_start_node; i_node < xe; ++i_node) {
	      
                PetscReal dx_dxi = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node][j_node-1][i_node].x + nodal_coords_arr[k_node-1][j_node][i_node].x + nodal_coords_arr[k_node-1][j_node-1][i_node].x) - (nodal_coords_arr[k_node][j_node][i_node-1].x + nodal_coords_arr[k_node][j_node-1][i_node-1].x + nodal_coords_arr[k_node-1][j_node][i_node-1].x + nodal_coords_arr[k_node-1][j_node-1][i_node-1].x) );

		PetscReal dy_dxi = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node][j_node-1][i_node].y + nodal_coords_arr[k_node-1][j_node][i_node].y + nodal_coords_arr[k_node-1][j_node-1][i_node].y) - (nodal_coords_arr[k_node][j_node][i_node-1].y + nodal_coords_arr[k_node][j_node-1][i_node-1].y + nodal_coords_arr[k_node-1][j_node][i_node-1].y + nodal_coords_arr[k_node-1][j_node-1][i_node-1].y) );

		PetscReal dz_dxi = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node][j_node-1][i_node].z + nodal_coords_arr[k_node-1][j_node][i_node].z + nodal_coords_arr[k_node-1][j_node-1][i_node].z) - (nodal_coords_arr[k_node][j_node][i_node-1].z + nodal_coords_arr[k_node][j_node-1][i_node-1].z + nodal_coords_arr[k_node-1][j_node][i_node-1].z + nodal_coords_arr[k_node-1][j_node-1][i_node-1].z) );

		PetscReal dx_deta = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node][j_node][i_node-1].x + nodal_coords_arr[k_node-1][j_node][i_node].x + nodal_coords_arr[k_node-1][j_node][i_node-1].x) - (nodal_coords_arr[k_node][j_node-1][i_node].x + nodal_coords_arr[k_node][j_node-1][i_node-1].x + nodal_coords_arr[k_node-1][j_node-1][i_node].x + nodal_coords_arr[k_node-1][j_node-1][i_node-1].x) );

		PetscReal dy_deta = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node][j_node][i_node-1].y + nodal_coords_arr[k_node-1][j_node][i_node].y + nodal_coords_arr[k_node-1][j_node][i_node-1].y) - (nodal_coords_arr[k_node][j_node-1][i_node].y + nodal_coords_arr[k_node][j_node-1][i_node-1].y + nodal_coords_arr[k_node-1][j_node-1][i_node].y + nodal_coords_arr[k_node-1][j_node-1][i_node-1].y) );

		PetscReal dz_deta = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node][j_node][i_node-1].z + nodal_coords_arr[k_node-1][j_node][i_node].z + nodal_coords_arr[k_node-1][j_node][i_node-1].z) - (nodal_coords_arr[k_node][j_node-1][i_node].z + nodal_coords_arr[k_node][j_node-1][i_node-1].z + nodal_coords_arr[k_node-1][j_node-1][i_node].z + nodal_coords_arr[k_node-1][j_node-1][i_node-1].z) );

		PetscReal dx_dzeta = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node][j_node-1][i_node].x + nodal_coords_arr[k_node][j_node][i_node-1].x + nodal_coords_arr[k_node][j_node-1][i_node-1].x) - (nodal_coords_arr[k_node-1][j_node][i_node].x + nodal_coords_arr[k_node-1][j_node-1][i_node].x + nodal_coords_arr[k_node-1][j_node][i_node-1].x + nodal_coords_arr[k_node-1][j_node-1][i_node-1].x) );

		PetscReal dy_dzeta = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node][j_node-1][i_node].y + nodal_coords_arr[k_node][j_node][i_node-1].y + nodal_coords_arr[k_node][j_node-1][i_node-1].y) - (nodal_coords_arr[k_node-1][j_node][i_node].y + nodal_coords_arr[k_node-1][j_node-1][i_node].y + nodal_coords_arr[k_node-1][j_node][i_node-1].y + nodal_coords_arr[k_node-1][j_node-1][i_node-1].y) );

		PetscReal dz_dzeta = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node][j_node-1][i_node].z + nodal_coords_arr[k_node][j_node][i_node-1].z + nodal_coords_arr[k_node][j_node-1][i_node-1].z) - (nodal_coords_arr[k_node-1][j_node][i_node].z + nodal_coords_arr[k_node-1][j_node-1][i_node].z + nodal_coords_arr[k_node-1][j_node][i_node-1].z + nodal_coords_arr[k_node-1][j_node-1][i_node-1].z) );

		PetscReal jacobian_det = dx_dxi * (dy_deta * dz_dzeta - dz_deta * dy_dzeta) - dy_dxi * (dx_deta * dz_dzeta - dz_deta * dx_dzeta) + dz_dxi * (dx_deta * dy_dzeta - dy_deta * dx_dzeta);
                if (PetscAbsReal(jacobian_det) < 1.0e-12) { SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FLOP_COUNT, "Jacobian is near zero..."); }
                aj_arr[k_node][j_node][i_node] = 1.0 / jacobian_det;
            }
        }
    }

    // --- 4. Boundary Extrapolation for Aj ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Extrapolating boundary values for Aj. \n");
    PetscInt i_bnd, j_bnd, k_bnd;

    if (xs == 0) {
        i_bnd = 0;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
                if (i_bnd + 1 < mx) aj_arr[k_bnd][j_bnd][i_bnd] = aj_arr[k_bnd][j_bnd][i_bnd+1];
            }
        }
    }
    if (xe == mx) {
        i_bnd = mx - 1;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
                if (i_bnd - 1 >= 0) aj_arr[k_bnd][j_bnd][i_bnd] = aj_arr[k_bnd][j_bnd][i_bnd-1];
            }
        }
    }
    // (Similar extrapolation blocks for Y and Z boundaries for aj_arr)
    if (ys == 0) {
        j_bnd = 0;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                 if (j_bnd + 1 < my) aj_arr[k_bnd][j_bnd][i_bnd] = aj_arr[k_bnd][j_bnd+1][i_bnd];
            }
        }
    }
    if (ye == my) {
        j_bnd = my - 1;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                 if (j_bnd - 1 >= 0) aj_arr[k_bnd][j_bnd][i_bnd] = aj_arr[k_bnd][j_bnd-1][i_bnd];
            }
        }
    }
    if (zs == 0) {
        k_bnd = 0;
        for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                if (k_bnd + 1 < mz) aj_arr[k_bnd][j_bnd][i_bnd] = aj_arr[k_bnd+1][j_bnd][i_bnd];
            }
        }
    }
    if (ze == mz) {
        k_bnd = mz - 1;
        for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                 if (k_bnd - 1 >= 0) aj_arr[k_bnd][j_bnd][i_bnd] = aj_arr[k_bnd-1][j_bnd][i_bnd];
            }
        }
    }

    // --- 5. Restore arrays ---
    ierr = DMDAVecRestoreArrayRead(user->fda, localCoords_from_dm, &nodal_coords_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->Aj, &aj_arr); CHKERRQ(ierr);

    // --- 6. Assemble Global Vector ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Assembling global Aj.\n");
    ierr = VecAssemblyBegin(user->Aj); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(user->Aj); CHKERRQ(ierr);

    // --- 7. Update Local Ghosted Version ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Updating local lAj.\n");
    ierr = UpdateLocalGhosts(user, "Aj"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Completed calculation, extrapolation, and update for Aj.\n");
    PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------------- */
/* End of Metric.c */
