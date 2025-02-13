#include <mpi.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <string.h>

/*
   Suppose `UserCtx` stores:
   - an array of bounding boxes, one per rank: user->bboxlist[rankID]
   - swarm pointer: user->swarm
   - ...
*/
typedef struct {
  DM             swarm;        // Your DMSwarm
  PetscInt       dim;          // Dimension: 2 or 3, etc.
  PetscInt       size;         // # of ranks in MPI
  BoundingBox   *bboxlist;     // Array of bounding boxes for each rank
  /* ... other fields ... */
} UserCtx;

/* A simple bounding box struct for each rank */
typedef struct {
  PetscReal min_coords[3];
  PetscReal max_coords[3];
} BoundingBox;

/* 
   MigrateParticles_Manual:
   Moves any out-of-domain particles to the correct rank by manual MPI exchange,
   assuming that you have an array of bounding boxes `bboxlist` (size = #ranks),
   where rank r covers the region [bboxlist[r].min_coords..bboxlist[r].max_coords].
*/
PetscErrorCode MigrateParticles_Manual(UserCtx *user)
{
  PetscErrorCode ierr;
  MPI_Comm       comm;
  PetscMPIInt    rank, size;
  DM             swarm = user->swarm;

  PetscReal     *pos    = NULL; 
  PetscInt       nLocal, p;
  
  /*--- 1) Basic MPI info ---*/
  ierr = PetscObjectGetComm((PetscObject)swarm, &comm); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm, &size); CHKERRQ(ierr);

  /*--- 2) Access the particle positions (and anything else needed) ---*/
  ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&pos); CHKERRQ(ierr);

  /*
     We'll build a "send-list" of particles that left our domain. For each outgoing
     particle, we need:
       - which rank it should go to
       - the (x,y,z) coordinates, plus velocity or other fields if needed
  */

  /* data needed for each outgoing particle */
  typedef struct {
    PetscInt  destRank;             // rank to which this particle should go
    PetscReal x[3];                // position
    // If you need velocity, cellID, etc., add them here
  } OutgoingParticle;

  /* We'll store them in a temporary array `outgoing`. Overallocate by nLocal in worst case. */
  OutgoingParticle *outgoing = (OutgoingParticle*) malloc(nLocal * sizeof(OutgoingParticle));
  PetscInt nOut = 0;

  /*--- 3) Identify particles which left local bounding box. ---*/
  BoundingBox mybbox  = user->bboxlist[rank]; // bounding box of *this* rank

  for (p = 0; p < nLocal; p++) {
    PetscReal *pp = &pos[p * user->dim]; /* pointer to i-th particle’s coords */
    PetscBool  inLocalDomain = PETSC_TRUE;

    /* Check if inside local bounding box */
    for (PetscInt d = 0; d < user->dim; d++) {
      if (pp[d] < mybbox.min_coords[d] || pp[d] >= mybbox.max_coords[d]) {
        inLocalDomain = PETSC_FALSE;
        break;
      }
    }

    if (!inLocalDomain) {
      /* (a) find which rank's bounding box this belongs to */
      PetscInt dest = -1;
      for (PetscInt r = 0; r < size; r++) {
        BoundingBox bb = user->bboxlist[r];
        /* check if (pp[0],pp[1],pp[2]) lies in [bb.min_coords .. bb.max_coords] */
        PetscBool inside = PETSC_TRUE;
        for (PetscInt d=0; d<user->dim; d++) {
          if (pp[d] < bb.min_coords[d] || pp[d] >= bb.max_coords[d]) {
            inside = PETSC_FALSE;
            break;
          }
        }
        if (inside) {
          dest = r;
          break;
        }
      }
      /* (b) Add to outgoing array if we found a destination rank */
      if (dest >= 0 && dest != rank) {
        outgoing[nOut].destRank = dest;
        outgoing[nOut].x[0]     = pp[0];
        if (user->dim>1) outgoing[nOut].x[1] = pp[1];
        if (user->dim>2) outgoing[nOut].x[2] = pp[2];
        // copy velocity or other fields if you want
        nOut++;
      }
    }
  }

  /*--- 4) We now must remove these out-of-domain particles from local storage. ---*/
  /* A simple approach is to build a compacted array of "kept" particles. */
  /* Alternatively, we can just DMSwarmRemovePoint() in a loop (slower if many). */

  /* For demonstration, let's call DMSwarmRemovePoint(); */
  for (PetscInt o = 0; o < nOut; o++) {
    /* We need to find the matching position in the swarm & remove it.
       Because we haven't removed them yet, the index is stable 
       only if we remove from last to first. Or we do a second pass. */
  }
  /* A robust approach: we do a second pass from the end. */
  for (p = nLocal-1; p >= 0; p--) {
    PetscReal *pp = &pos[p * user->dim];
    PetscBool  inLocal = PETSC_TRUE;
    for (PetscInt d = 0; d < user->dim; d++) {
      if (pp[d] < mybbox.min_coords[d] || pp[d] >= mybbox.max_coords[d]) {
        inLocal = PETSC_FALSE;
        break;
      }
    }
    if (!inLocal) {
      ierr = DMSwarmRemovePointAtIndex(swarm, p); CHKERRQ(ierr);
    }
  }

  /* After removal, restore field array and reacquire if you plan to do more with it. */
  ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&pos); CHKERRQ(ierr);

  /*--- 5) Now we do an MPI exchange to send outgoing[] to the correct ranks. ---*/

  /* Build array "sendcounts" so that each rank knows how many we will send it. */
  PetscInt *sendcounts = (PetscInt*) calloc(size, sizeof(PetscInt));
  for (PetscInt o = 0; o < nOut; o++) {
    sendcounts[outgoing[o].destRank]++;
  }

  /* Build displacements in 'sdispls' */
  PetscInt *sdispls = (PetscInt*) malloc(size * sizeof(PetscInt));
  sdispls[0] = 0;
  for (PetscInt r = 1; r < size; r++) {
    sdispls[r] = sdispls[r-1] + sendcounts[r-1];
  }

  /* Allocate a "sortedOut" array in which we group all particles by destRank. */
  OutgoingParticle *sortedOut = (OutgoingParticle*) malloc(nOut * sizeof(OutgoingParticle));
  /* fill it */
  int *curPos = (int*) calloc(size, sizeof(int));
  for (PetscInt o=0; o<nOut; o++) {
    PetscInt r = outgoing[o].destRank;
    PetscInt idx = sdispls[r] + curPos[r];
    sortedOut[idx] = outgoing[o];
    curPos[r]++;
  }
  free(curPos);

  /* We'll do an AllToAll to let each rank know how many we get from each rank. */
  PetscInt *recvcounts = (PetscInt*) calloc(size, sizeof(PetscInt));
  ierr = MPI_Alltoall(sendcounts, 1, MPI_INT,    /* send 1 integer to each rank */
                      recvcounts, 1, MPI_INT, comm);  /* receive 1 integer from each rank */
  CHKERRMPI(ierr);

  /* Build 'rdispls' for the incoming data */
  PetscInt *rdispls = (PetscInt*) malloc(size * sizeof(PetscInt));
  rdispls[0] = 0;
  PetscInt nIn = recvcounts[0];
  for (PetscInt r=1; r<size; r++){
    rdispls[r] = rdispls[r-1] + recvcounts[r-1];
    nIn += recvcounts[r];
  }

  /* We'll receive nIn OutgoingParticles from various ranks. */
  OutgoingParticle *inBuf = (OutgoingParticle*) malloc(nIn * sizeof(OutgoingParticle));

  /* Exchange the actual data via MPI_Alltoallv */
  ierr = MPI_Alltoallv(sortedOut, sendcounts, sdispls, MPI_CHAR,
                       inBuf,      recvcounts, rdispls, MPI_CHAR,
                       comm); CHKERRMPI(ierr);

  /* free the intermediate arrays we no longer need */
  free(sortedOut);
  free(sendcounts);
  free(sdispls);
  /* keep recvcounts, rdispls for now until we have inBuf */

  /*--- 6) Add the newly incoming particles to our local DMSwarm. ---*/
  ierr = DMSwarmAddNPoints(swarm, DMSWARM_INSERT_RANDOM, nIn); CHKERRQ(ierr);

  /* Reacquire the position array after adding extra slots. */
  ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&pos); CHKERRQ(ierr);

  /* The “end” of the local array is the region for new points, 
     i.e. old nLocal..nLocal+(nIn-1). */
  PetscInt newStart = 0;
  ierr = DMSwarmGetLocalSize(swarm, &newStart); CHKERRQ(ierr);
  newStart -= nIn; /* the first new one is at [newStart]. */

  for (PetscInt i=0; i<nIn; i++) {
    PetscInt idx = newStart + i;
    PetscReal *pp = &pos[idx * user->dim];
    pp[0] = inBuf[i].x[0];
    if (user->dim>1) pp[1] = inBuf[i].x[1];
    if (user->dim>2) pp[2] = inBuf[i].x[2];
    /* If you stored velocity, copy that in as well, etc. */
  }

  /* Clean up & restore */
  free(inBuf);
  free(recvcounts);
  free(rdispls);

  ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&pos); CHKERRQ(ierr);

  free(outgoing);

  PetscFunctionReturn(0);
}
