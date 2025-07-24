#include <petsc.h>
#include <stdio.h>
#include <math.h>
#include "walkingsearch.h"  // or your path to the prototypes

// Simple macros for pass/fail
#define TEST_ASSERT(cond, desc) do {                     \
  if (!(cond)) {                                         \
    PetscPrintf(PETSC_COMM_WORLD, "[TEST] FAIL: %s\n", desc); \
    fails++;                                             \
  } else {                                               \
    PetscPrintf(PETSC_COMM_WORLD, "[TEST] PASS: %s\n", desc); \
    passes++;                                            \
  }                                                      \
} while(0)

static int passes = 0;
static int fails  = 0;

// Helper to create a minimal mock "UserCtx" and DM for testing local cells
// For actual usage, you may already have a more robust user->da creation.
PetscErrorCode CreateMockUserCtx(UserCtx **userOut, PetscInt IM, PetscInt JM, PetscInt KM)
{
  PetscErrorCode ierr;
  UserCtx *user;
  ierr = PetscMalloc1(1, &user); CHKERRQ(ierr);
  // Suppose your user->da is created as a 3D DMDA:
  // If you want 2D or 1D, pass JM=2 or KM=2, etc.
  
  // Initialize user
  user->IM = IM; user->JM = JM; user->KM = KM;
  // (some other fields to zero or defaults)
  user->da = NULL;
  user->fda = NULL;

  // Create DMDA
  // For a minimal test, we do single rank and small block:
  ierr = DMDACreate3d(PETSC_COMM_WORLD,
                      DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                      DMDA_STENCIL_BOX,
                      IM+1, JM+1, KM+1,  // Global dims
                      1, 1, 1,          // decomposition
                      1,                // dof
                      1,                // stencil width
                      NULL, NULL, NULL, // partition
                      &user->da); CHKERRQ(ierr);

  ierr = DMSetUp(user->da); CHKERRQ(ierr);
  // For coordinate DM
  ierr = DMGetCoordinateDM(user->da, &user->fda); CHKERRQ(ierr);

  // Assign uniform coords [0..1] in x,y,z
  ierr = DMDASetUniformCoordinates(user->da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0); CHKERRQ(ierr);

  *userOut = user;
  return 0;
}

// Minimal function to tear down the DM we created
PetscErrorCode DestroyMockUserCtx(UserCtx *user)
{
  if (!user) return 0;
  if (user->da)  { DMDestroy(&user->da); }
  // user->fda is destroyed together with user->da typically
  PetscFree(user);
  return 0;
}

// ============ Tests begin here ============

// Test 1: Degenerate Plane (Corner Case #1)
PetscErrorCode TestDegeneratePlane()
{
  PetscErrorCode ierr;
  // We'll build a "Cell" that duplicates a vertex so the plane normal is zero
  Cell degenerateCell;
  // Simple: all corners the same => definitely degenerate
  for (int i=0; i<8; i++) {
    degenerateCell.vertices[i].x = 0.0;
    degenerateCell.vertices[i].y = 0.0;
    degenerateCell.vertices[i].z = 0.0;
  }

  PetscReal d[NUM_FACES];
  Cmpnts p = {0.0,0.0,0.0};

  // If we call CalculateDistancesToCellFaces, we expect PETSC_ERR_PLANE_DEGENERATE
  ierr = CalculateDistancesToCellFaces(p, &degenerateCell, d, 1e-14);
  if (ierr == PETSC_ERR_PLANE_DEGENERATE) {
    TEST_ASSERT(1, "Degenerate plane => correct error code");
    ierr = 0;
  } else if (ierr) {
    TEST_ASSERT(0, "Degenerate plane => unexpected error code");
  } else {
    TEST_ASSERT(0, "Degenerate plane => expected failure but got success");
  }

  return ierr;
}

// Test 2: Particle near shared boundary toggling => we skip
// We'll do a smaller mock scenario
PetscErrorCode TestBoundaryToggle()
{
  // We'll create a 2x2x2 domain => 1 cell in x,y,z
  // Then place a particle near boundary
  UserCtx *user=NULL;
  PetscErrorCode ierr;
  ierr = CreateMockUserCtx(&user, 1,1,1); CHKERRQ(ierr);

  // We have exactly one cell => i=0..1, j=0..1, k=0..1
  // Place a particle near x=1.0 => if distance is borderline negative or zero
  Particle part;
  part.loc.x = 0.9999999999999; // near boundary
  part.loc.y = 0.5;
  part.loc.z = 0.5;

  PetscReal d[NUM_FACES];
  ierr = LocateParticleInGrid(user, &part, d); // We rely on your code
  // If it toggled or messed up, it might not find a cell or might break loop
  // Expect it to find boundary => cell[0]=0,0,0
  TEST_ASSERT( (part.cell[0] >= 0), 
               "Boundary toggle test => found a valid cell (treated boundary as inside)" );

  // Cleanup
  ierr = DestroyMockUserCtx(user); CHKERRQ(ierr);
  return ierr;
}

// More tests: Edge/corner, toggling, 1D/2D usage, etc. You can add them similarly.

// Main test harness
int main(int argc, char **argv)
{
  PetscInitialize(&argc, &argv, NULL, NULL);
  passes=0; fails=0;

  // We do each corner case test in a function:

  // 1) Degenerate plane
  TestDegeneratePlane();

  // 2) Boundary toggle test
  TestBoundaryToggle();

  // 3) Edge/corner => we can do a small domain, place a particle at corner exactly.
  // 4) ...
  // etc. (Add your other corner-case tests similarly.)

  PetscPrintf(PETSC_COMM_WORLD, 
    "\n=====================================\n"
    "Test Summary: passes=%d, fails=%d\n"
    "=====================================\n", passes, fails);

  PetscFinalize();
  return (fails>0) ? 1 : 0;
}
