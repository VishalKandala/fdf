// ParticleMotion.c

#include "ParticleMotion.h"

/**
 * @brief Updates a particle's position based on its velocity and the timestep dt (stored in user->dt).
 *
 * @param[in]     user     Pointer to your UserCtx (must contain user->dt).
 * @param[in,out] position Pointer to the particle's current position (Cmpnts).
 * @param[in]     velocity Pointer to the particle's velocity (Cmpnts).
 *
 * @return PetscErrorCode  Returns 0 on success, or an error code on failure.
 */
PetscErrorCode UpdateParticlePosition(UserCtx *user, Cmpnts *position, const Cmpnts *velocity)
{
  PetscFunctionBeginUser; // PETSc macro for error/stack tracing

  /* Update the position with velocity * dt */
  position->x += velocity->x * user->dt;
  position->y += velocity->y * user->dt;
  position->z += velocity->z * user->dt;

  PetscFunctionReturn(0);
}

/**
 * @brief Loops over all local particles in the DMSwarm, updating their positions
 *        based on velocity and the global timestep user->dt.
 *
 * @param[in,out] user    Pointer to UserCtx (must contain dt).
 *
 * @return PetscErrorCode Returns 0 on success, or an error code on failure.
 */
PetscErrorCode UpdateAllParticlePositions(UserCtx *user)
{
  PetscErrorCode ierr;
  DM swarm = user->swarm;
  PetscInt       nLocal, p;
  Cmpnts        *pos = NULL;
  Cmpnts        *vel = NULL;

  PetscFunctionBeginUser;  // PETSc macro for error/stack tracing

  // 1) Get the number of local particles
  ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);

  // 2) Access the "position" and "velocity" fields
  ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&pos); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&vel); CHKERRQ(ierr);

  // 3) Loop over all local particles, updating each position by velocity * dt
  for (p = 0; p < nLocal; p++) {
    ierr = UpdateParticlePosition(user, &pos[p], &vel[p]); CHKERRQ(ierr);
  }

  // 4) Restore the fields
  ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&pos); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&vel); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
