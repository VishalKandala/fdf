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

/* ------------------------------------------------------------------------- */
/* End of Metric.c */
