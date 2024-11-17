static char help[] = "DMSwarm Interpolation - fdf-curvIB ";

#include <petscpf.h>
#include <petscdmswarm.h>
#include "variables.h"
#include <stdlib.h>
#include "time.h"
#include <math.h>
//#include "parallelComm.h"
#include "petsctime.h"
#include "petscdmcomposite.h"

PetscInt np = 0; // No.of Particles
PetscInt ti = 0;
PetscReal L_dim = 1.0,cl = 1.0;
PetscInt block_number = 1;
PetscInt visflg = 0;

//////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


// Read Velocity field - ufield //

PetscErrorCode Ucat_Binary_Input(UserCtx *user)
{
  PetscViewer viewer;
  char filen[90];
  PetscInt bi=user->_this;
  PetscInt ctr = 7;
  
  sprintf(filen, "results/ufield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);

  PetscInt N;

  VecGetSize(user->Ucat, &N);
  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "Ucat_Binary_Input - user,SizeOf(ucat) - %p,%d \n", user,N);
  VecLoad((user->Ucat),viewer);
 
  PetscViewerDestroy(&viewer);

  PetscBarrier(PETSC_NULL);

  return 0;

}

/////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


// Function to create bounding boxes for processors and gather them on rank 0

PetscErrorCode  WriteBoundingBox(UserCtx *user, BoundingBox *local_bbox) {

    PetscErrorCode ierr;
    PetscInt i, j, k, rank, size;
    PetscInt xs, ys, zs, xe, ye, ze;
    DMDALocalInfo info;
    Vec Coor;
    Cmpnts ***coor, min_coords, max_coords;
    PetscInt ctr = 0;

    // Get the rank and size of the communicator
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, " ***************************** \n");

    if (visflg == ctr) PetscSynchronizedPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - user - %p \n", user);
    PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

    DM da = user->da;
    DM fda = user->fda;

    //   PetscBarrier(PETSC_NULL);

    ierr = DMGetCoordinatesLocal(user->da, &Coor); CHKERRQ(ierr);

    if (Coor == NULL) {
        if (visflg == ctr) PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.34F - Error: Coor vector is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Coordinates vector is NULL");
    }

    ierr = DMDAVecGetArrayRead(user->fda, Coor, &coor); CHKERRQ(ierr);

    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

    xs = info.xs; xe = xs + info.xm;
    ys = info.ys; ye = ys + info.ym;
    zs = info.zs; ze = zs + info.zm;



    //    if (visflg == ctr) PetscSynchronizedPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.36 - rank - %d; xs - %d; xe - %d; ys - %d; ye - %d; zs - %d; ze - %d \n",
    //                                   rank, xs, xe, ys, ye, zs, ze);
    //  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);




    // Initialize min and max coordinates
    min_coords.x = min_coords.y = min_coords.z = PETSC_MAX_REAL;
    max_coords.x = max_coords.y = max_coords.z = PETSC_MIN_REAL;

    // Find min and max coordinates for the local grid
    for (k = zs; k < ze; k++) {
        for (j = ys; j < ye; j++) {
            for (i = xs; i < xe; i++) {
                // Update min and max coordinates using the coor array
                min_coords.x = PetscMin(min_coords.x, coor[k][j][i].x);
                min_coords.y = PetscMin(min_coords.y, coor[k][j][i].y);
                min_coords.z = PetscMin(min_coords.z, coor[k][j][i].z);

                max_coords.x = PetscMax(max_coords.x, coor[k][j][i].x);
                max_coords.y = PetscMax(max_coords.y, coor[k][j][i].y);
                max_coords.z = PetscMax(max_coords.z, coor[k][j][i].z);
            }
        }
    }

    ierr = DMDAVecRestoreArrayRead(fda, Coor, &coor); CHKERRQ(ierr);

    // Create local bounding box
    local_bbox->min_coords = min_coords;
    local_bbox->max_coords = max_coords;

    // if(visflg==ctr){
    //     PetscSynchronizedPrintf(PETSC_COMM_WORLD, "WriteBoundingBox - - rank - %d; coords: x=[%f,%f] y=[%f,%f] z=[%f,%f] \n",rank,bbox_array[i].min_coords.x, bbox_array[i].max_coords.x,bbox_array[i].min_coords.y,bbox_array[i].max_coords.y,bbox_array[i].min_coords.z, bbox_array[i].max_coords.z);
    //    PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
    //  }// visflg

    user->bbox = *local_bbox; // update bounding box inside local UserCtx;
    // Allocate memory on rank 0 for all bounding boxes1

    PetscBarrier(PETSC_NULL);

   if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, " ***************************** \n");

   return 0;
}


GatherBoundingBoxes(UserCtx* user,BoundingBox **all_bboxes){
  PetscInt rank,size,ctr=0;
  PetscErrorCode ierr;
  BoundingBox* bbox_array = NULL;
  BoundingBox local_bbox;

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
  
  if (rank == 0) {
        bbox_array = (BoundingBox *)malloc(size * sizeof(BoundingBox));
        if (!bbox_array) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_MEM, "WriteBoundingBox -- Memory allocation failed for bounding box array on rank 0");
        }
       if (visflg == ctr) PetscPrintf(PETSC_COMM_SELF, "GatherBoundingBoxes - rank - %d - bbox_array %p \n",rank,bbox_array);

    }// rank 0

  WriteBoundingBox(user,&local_bbox); // Write Local Bounding Box for each rank.
    
    // Gather local bounding boxes from all processors to rank 0
    ierr = MPI_Gather(&local_bbox, sizeof(BoundingBox), MPI_BYTE,
                      bbox_array, sizeof(BoundingBox), MPI_BYTE,
                      0, PETSC_COMM_WORLD); CHKERRQ(ierr);

    // On rank 0, assign gathered bounding boxes to output pointer
    if (rank == 0) {
        *all_bboxes = bbox_array;

        // Print all gathered bounding boxes for verification
        for (PetscInt i = 0; i < size; i++) {
	  if(visflg==ctr) { PetscPrintf(PETSC_COMM_WORLD, "GatherBoundingBoxes - Gathered Data - Host Rank %d - Receiver Rank - %d : coords: x=[%f,%f] y=[%f,%f] z=[%f,%f] \n", i+1,i,rank,bbox_array[i].min_coords.x, bbox_array[i].max_coords.x,bbox_array[i].min_coords.y,bbox_array[i].max_coords.y,bbox_array[i].min_coords.z, bbox_array[i].max_coords.z);
	  }
        }
    }else *all_bboxes = NULL;  // rank 0 
  
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

PetscBool  CPUPointIntersectCheck(BoundingBox *bbox, Particle *particle){
   
   Cmpnts loc,min_coords,max_coords;
   PetscErrorCode ierr;
   PetscInt ctr = 0;
   PetscInt rank;
   loc = particle->loc;
   PetscBool Intersects = PETSC_FALSE;
   min_coords = bbox->min_coords;
   max_coords = bbox->max_coords;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
   
    // if(visflg==ctr) PetscPrintf(PETSC_COMM_SELF, "CPUPointIntersectCheck - rank - %d particle - %d ; loc [%f,%f,%f]|  bbox:x [%f,%f] | y - [%f,%f] |z - [%f,%f] \n",rank,particle->PID,loc.x,loc.y,loc.z,min_coords.x,max_coords.x,min_coords.y,max_coords.y,min_coords.z,max_coords.z);

    if ((loc.x >= min_coords.x && loc.x <= max_coords.x) &&
        (loc.y >= min_coords.y && loc.y <= max_coords.y) &&
        (loc.z >= min_coords.z && loc.z <= max_coords.z)) {
        Intersects = PETSC_TRUE;
    }
   return Intersects;
}

///////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// Walking Search within a rank //////////////////////////////////////////

PetscErrorCode  InterpolationWeightsCalculate(Cmpnts a, PetscReal d[6])
{
  a.x = d[0]/(d[0]+d[1]);
  a.y = d[2]/(d[2]+d[3]);
  a.z = d[4]/(d[4]+d[5]);
  return(0);
}

// Function to create bounding boxes for processors, and locate particles in the grid and calculate interpolation weights for each particle.


PetscErrorCode ParticlesLocate(UserCtx *user) {
    PetscErrorCode ierr;
    PetscInt rank;
    PetscMPIInt size;
    PetscInt ctr = 0;
    PetscBool Particle_Detected = PETSC_FALSE;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    //  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - user - %p \n",user);

    PetscBarrier(PETSC_NULL);

    //  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - bbox complete %p \n",&bboxlist);

    // Access DMSwarm fields

    PetscInt localNumParticles;
    PetscReal *positions;
    PetscReal *weights;
    PetscInt64  *cellIndices,*PIDs;

    DM swarm = user->swarm;

    ierr = DMSwarmGetLocalSize(swarm, &localNumParticles); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIndices); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&PIDs); CHKERRQ(ierr);
    // Access other fields as needed

    //    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - swarm fields accessed \n");

    Particle particle;

    // Loop through each particle
    for (PetscInt i = 0; i < localNumParticles; i++) {

        // Populate the Particle struct with data from DMSwarm fields

        // Initialize PID
        particle.PID = PIDs[i];
        // Initialize weights
        particle.weights.x = weights[3*i + 0];
	particle.weights.y = weights[3*i + 1];
	particle.weights.z = weights[3*i + 2];
        // Initialize locations
        particle.loc.x = positions[3*i + 0];
	particle.loc.y = positions[3*i + 1];
	particle.loc.z = positions[3*i + 2];
        // Initialize velocities
        particle.vel.x = 0.0;
	particle.vel.y = 0.0;
	particle.vel.z = 0.0;
	//Initialize cell indices

	particle.cell[0] = cellIndices[3*i + 0];
	particle.cell[1] = cellIndices[3*i + 1];
	particle.cell[2] = cellIndices[3*i + 2];

	//  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate -  particle - [%d] - %p ; particle.loc.x - %f, positions.x - %f \n",i,&particle,particle.loc.x,positions[3*i+0]);
	//   if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate  - particle - [%d] - %p ; particle.cell[0,1,2] - [%d,%d,%d] - cellIndices[0,1,2] - [%d,%d,%d]  \n",i,&particle,particle.cell[0],particle.cell[1],particle.cell[2],cellIndices[3*i + 0],cellIndices[3*i + 1],cellIndices[3*i + 2]);

        //Check if the particle intersects the bounding box
  
	Particle_Detected = CPUPointIntersectCheck(&(user->bbox), &particle);
	
	if (Particle_Detected == PETSC_TRUE) {
            ierr = LocateParticleInGrid(user, &particle); CHKERRQ(ierr);
	  }

	// Restore data into the swarm
        weights[3*i + 0] = particle.weights.x;
        weights[3*i + 1] = particle.weights.y;
        weights[3*i + 2] = particle.weights.z;

        cellIndices[3*i + 0] = particle.cell[0];
        cellIndices[3*i + 1] = particle.cell[1];
        cellIndices[3*i + 2] = particle.cell[2];

	//     if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - cp10.7A%d - particle - [%d] - %p ; particle.loc.x - %f,\n",i+1,i,&particle,particle.loc.x);
	//  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - cp10.8A%d - particle - [%d] - %p ; particle.cell[0,1,2] - [%d,%d,%d] \n",i+1,i,&particle,particle.cell[0],particle.cell[1],particle.cell[2]);

    } // particle loop

    // Restore DMSwarm fields
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIndices); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&PIDs); CHKERRQ(ierr);

    //  free(bboxlist);

    return 0;
}

//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **argv){

PetscErrorCode ierr;
  
 DM da,fda;

UserCtx *user;

PetscInt ctr = 0; // Output counter is at 0 if only main outputs are necessary.

PetscInt rank, size,bi, ibi,i;

PetscReal umax;

BoundingBox *bboxlist;
  
ierr = PetscInitialize(&argc, &argv, (char *)0, help); if(ierr) return ierr;

PetscOptionsInsertFile(PETSC_COMM_WORLD, PETSC_NULL, "control.dat", PETSC_TRUE);

PetscOptionsGetInt(PETSC_NULL,PETSC_NULL, "-visflg", &visflg, PETSC_NULL);
//if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, " ***************************** \n");

// Allocate memory for user

PetscMalloc1(block_number, &user); // PetscMalloc1 automatically uses the sizeof(UserCtx) to allocate memory.

PetscBarrier(PETSC_NULL);

// Identify all processors 

//////////////////////////////////////////////////////////////////////////////

MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

MPI_Comm_size(PETSC_COMM_WORLD, &size);

PetscBarrier(PETSC_NULL);

//////////////////////////////////////////////////////////////////////////////

// Read no.of particles from control file 

PetscOptionsGetInt(PETSC_NULL,PETSC_NULL, "-ti", &ti, PETSC_NULL);

PetscOptionsGetInt(PETSC_NULL,PETSC_NULL, "-numParticles", &np, PETSC_NULL);

if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, " ***************************** \n");

//////////////////////////////////////////////////////////////////////////////

// Output the inputs & Allocations.

 if(visflg==ctr) {
PetscSynchronizedPrintf(PETSC_COMM_WORLD, "main - user,rank - %p,%d  \n",user,rank);
PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

PetscPrintf(PETSC_COMM_WORLD, "main -  ti,no.of particles %d \n", ti,np);
 }

/////////////////////////////////////////////////////////////////////////////////

// Define coordinates

 DefineGridCoordinates(user);


 ///////////////////////////////////////////////////////////////////////////////////

// Create Vector for ucat

for(bi = 0; bi < block_number; bi ++ ){

  ierr = DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat); CHKERRQ(ierr);
}

PetscBarrier(PETSC_NULL);

/////////////////////////////////////////////////////////////////////////////////

// Read Ucat from ufield

Ucat_Binary_Input(user);

 ierr = VecNorm(user->Ucat,NORM_INFINITY,&umax); CHKERRQ(ierr);

if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "main - max(ucat) - %f \n",umax);

 PetscBarrier(PETSC_NULL);

 ///////////////////////////////////////////////////////////////////////////////

// Create Bounding Boxes for each rank,store them in user as well as in rank 0.

// ierr =  GatherBoundingBoxes(user,&bboxlist); CHKERRQ(ierr);

// Create and Initialize Particle //

// ierr = CreateParticleSwarm(user,np);

/*
 if (visflg == ctr){
 PetscPrintf(PETSC_COMM_WORLD, "main - ParticleSwarm: Pre-Search. \n");
 PetscPrintf(PETSC_COMM_WORLD, "\n");
 }
 ierr = PrintParticlePositions(user);

 // Particles are located in the grid

 ierr =  ParticlesLocate(user);

 if (visflg == ctr){
  PetscPrintf(PETSC_COMM_WORLD, "main - ParticleSwarm: Post-Search. \n");
  PetscPrintf(PETSC_COMM_WORLD, "\n");
 }

 ierr = PrintParticlePositions(user);

// Finalize 

*/

 for(bi = 0; bi < block_number; bi ++){
   
   DMDestroy(&(user[bi].fda));

   DMDestroy(&(user[bi].da)); 

   DMDestroy(&(user[bi].swarm));
 }

 //if(rank == 0) free(bboxlist);

PetscFree(user);

PetscFinalize();

return 0;
}
