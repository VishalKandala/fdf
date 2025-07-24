#include <petscsys.h>
#include <petscdm.h>
#include <petscdmswarm.h>
#include <petscdmda.h>

#define CHKERR(func) do { PetscErrorCode ierr = func; if (ierr) { PetscError(PETSC_COMM_SELF, __LINE__, PETSC_FUNCTION_NAME, __FILE__, ierr, PETSC_ERROR_REPEAT, " "); } } while (0)

int main(int argc, char **argv) {
  DM             dm;              // DM object for DMSwarm
  DM             da;              // DM object for DMDA (background grid)
  PetscInt       dim = 2;         // Dimension
  PetscInt       particles_per_cell = 5;  // Number of particles per cell
  PetscMPIInt    rank;
  Vec            coor;
  PetscReal      domain_min[2] = {0.0, 0.0};
  PetscReal      domain_max[2] = {1.0, 1.0};
  PetscInt       M = 4, N = 4;    // Number of grid cells in x and y direction

  PetscFunctionBeginUser;
  CHKERR(PetscInitialize(&argc, &argv, NULL, NULL));
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // Create a DMDA to represent the background domain
  CHKERR(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX, M, N, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &da));
  CHKERR(DMSetFromOptions(da));
  CHKERR(DMSetUp(da));

  // Create a DMSwarm to represent the particle system
  CHKERR(DMSwarmCreate(PETSC_COMM_WORLD, &dm));
  CHKERR(DMSetDimension(dm, dim));
  CHKERR(DMSwarmSetType(dm, DMSWARM_PIC));
  CHKERR(DMSwarmSetCellDM(dm, da));

  // Register particle fields if needed (e.g., velocities, IDs)
  CHKERR(DMSwarmSetLocalSizes(dm, particles_per_cell * M * N, 4));

  // Set up the swarm
  CHKERR(DMSwarmSetFromOptions(dm));
  CHKERR(DMSwarmSetUp(dm));
  CHKERR(DMSwarmInitializeCoordinates(dm));

  // Migrate particles according to domain decomposition
  CHKERR(DMSwarmMigrate(dm, PETSC_TRUE));

  // Get and print the particle coordinates to verify the result
  CHKERR(DMSwarmCreateGlobalVectorFromField(dm, "DMSwarmPIC_coor", &coor));
  CHKERR(VecView(coor, PETSC_VIEWER_STDOUT_WORLD));
  CHKERR(DMSwarmDestroyGlobalVectorFromField(dm, "DMSwarmPIC_coor", &coor));

  // Clean up
  CHKERR(DMDestroy(&da));
  CHKERR(DMDestroy(&dm));
  CHKERR(PetscFinalize());
  return 0;
}
