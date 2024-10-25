static char help[] = "Interpolation - fdf-curvIB ";

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

//extern PetscInt numParticles;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// Creation of Particle Vector ///////////////////////////////////////////////////////////

// Function to create a PETSc vector for particles

PetscErrorCode ParticleVectorCreate(PetscInt numParticles, UserCtx *user) {
  PetscErrorCode ierr;
  PetscInt ctr = 8;
  const PetscInt particleSize = (sizeof(Particle) / sizeof(PetscScalar));

  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.1 - user - %p  \n",user);

  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.2 - particleSize - %d  \n",particleSize);  

    if (!user) {
        if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.3F - Error: user is NULL\n");
        return PETSC_ERR_ARG_NULL;
    }

      if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorCreate - cp8.3S - ParticleVec before creation: %p \n", &user->ParticleVec);
      if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ParticelVectorCrate - cp8.4 - Local ParticleVec before creation: %p \n", &user->lParticleVec);    

    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorCreate - cp8.5 \n"); 
    
    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorCreate - cp8.6 - No.of Particles - %d \n",numParticles); 

    if (user->ParticleVec == NULL) {
        if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.7F - ParticleVec does not exist - Creating ParticleVec...\n");
	ierr = VecCreate(PETSC_COMM_WORLD,&user->ParticleVec); CHKERRQ(ierr);  // Check the PETSc error code and return if failed 
	if (ierr) {
	  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.71F - Error in DMCreateGlobalVector: %d\n", ierr);

	} else {
	  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.71S - Global vector created successfully.\n");
	}

        ierr = VecSetSizes(user->ParticleVec, PETSC_DECIDE, particleSize*numParticles); CHKERRQ(ierr);
        if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.72 -  vector size set \n"); 
   
        ierr = VecSetFromOptions(user->ParticleVec); CHKERRQ(ierr);
        if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.73 -  vector initialized to be set from options \n"); 

        ierr = VecSetType(user->ParticleVec,VECMPI); CHKERRQ(ierr);
        if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.74 -  vector initialized to be set from options \n"); 


	ierr = VecSet(user->ParticleVec,0.0); CHKERRQ(ierr);
    	if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.75 -  vector initialized to zeros \n");    

    } else {
        if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.7S - ParticleVec already exists.\n");
    }

    /// Diagnostics

    PetscInt N;
    PetscInt low,high;
    PetscInt rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
   
    VecGetSize(user->ParticleVec, &N);

    // Get the ownership range for this rank
    ierr = VecGetOwnershipRange(user->ParticleVec, &low, &high); CHKERRQ(ierr);

    // Synchronized output to display ownership information
    if(visflg==ctr){
    ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,
				   "ParticleVectorCreate - cp8.8%d - Rank %d owns indices from %d to %d\n", rank+1, rank, (int)low, (int)(high - 1)); CHKERRQ(ierr);
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
    }
    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorCreate - cp8.9 - sizeof(ParticeVec) - %d \n",N);    
    
    return(0);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// Function to initialize particle vector with random locations, zero velocity, and zero weight

PetscErrorCode ParticleVectorInsert(UserCtx *user, PetscInt index, Particle *particle)
{
  PetscInt ctr = 9;
  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "  ParticleVectorInsert - cp9.6%d1 - user - %p \n", index,user);
  
  const PetscInt particleSize = sizeof(Particle) / sizeof(PetscScalar);

  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "  ParticleVectorInsert - cp9.6%d2 - particleSize - %d \n", index+1,particleSize);

 if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "  ParticleVectorInsert - cp9.6%d3 - particle - %p \n", index+1,particle);

  PetscInt idx[particleSize];
    
  PetscScalar values[particleSize];
  
  PetscErrorCode ierr;
 
  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "  ParticleVectorInsert - cp9.6%d4 - idx - %p - values - %p \n", index+1,&idx,&values);

  values[0] = particle->cell[0];
  values[1] = particle->cell[1];
  values[2] = particle->cell[2];
  values[3] = particle->loc.x;
  values[4] = particle->loc.y;
  values[5] = particle->loc.z;
  values[6] = particle->vel.x;
  values[7] = particle->vel.y;
  values[8] = particle->vel.z;
  values[9] = particle->weights.x;
  values[10] = particle->weights.y;
  values[11] = particle->weights.z;
  
  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "  ParticleVectorInsert - cp9.6%d4 - values updated - sample - values[3]- %f == loc.x - %f \n", index+1, values[3],particle->loc.x);

  for(PetscInt j = 0; j< particleSize;j++){
     
    idx[j] = index*particleSize+j;
  
    ierr = VecSetValues(user->ParticleVec,1,&idx[j],&values[j],INSERT_VALUES); CHKERRQ(ierr);

  }

  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "  ParticleVectorInsert - cp9.6%d5 - Values set into ParticleVec - %p \n", index+1,&user->ParticleVec);
  
  return(0);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


PetscErrorCode ParticleVectorInitialize(UserCtx *user, PetscInt numParticles) {
  PetscErrorCode ierr;
  PetscRandom rand;
  PetscInt ctr = 9;
   
  PetscInt particleSize = sizeof(Particle) / sizeof(PetscScalar);
  
    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.1 - user - %p \n",user);
    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.2 - No.of Particles = %d \n",numParticles); 
    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.3 - ParticleSize = %d \n",particleSize);
    // Seed the random number generator
    srand(time(NULL));

    // Check if ParticleVec is NULL
      if (!user->ParticleVec) {
        if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorInitialize - cp9.4F - Error: user->ParticleVec is NULL\n");
        return PETSC_ERR_ARG_NULL;
     }

    // Create a PETSc random number generator
    PetscRandomCreate(PETSC_COMM_WORLD, &rand);
    PetscRandomSetType(rand, PETSCRAND48);
    PetscRandomSetInterval(rand,0.0 ,1.0 );
    PetscRandomSeed(rand);

    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.4S - Random Seed - %d \n",(PetscInt)rand); 

    // Initialize particle data
    Particle *particleData = (Particle *)calloc(numParticles,sizeof(Particle));

if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.5 - particleData  - %p \n",particleData);  

    for (PetscInt i = 0; i < numParticles; i++) {
        
      if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.5%dB - particleData[i]: %p - i: %d - loc.x - %f, loc.y - %f, loc.z - %f \n",i+1,particleData+i,i,particleData[i].loc.x, particleData[i].loc.y, particleData[i].loc.z); 
        
        PetscRandomGetValue(rand, &particleData[i].loc.x);
        PetscRandomGetValue(rand, &particleData[i].loc.y);
        PetscRandomGetValue(rand, &particleData[i].loc.z);

        particleData[i].vel.x = 0.0;
        particleData[i].vel.y = 0.0;
        particleData[i].vel.z = 0.0;
        particleData[i].weights.x = 0.0;
        particleData[i].weights.y = 0.0;
        particleData[i].weights.z = 0.0;
        particleData[i].cell[0] = -1;
        particleData[i].cell[1] = -1;
        particleData[i].cell[2] = -1;
        

      if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.51%dA - i  %d : loc.x - %f, loc.y - %f, loc.z - %f \n",i+1,i,particleData[i].loc.x, particleData[i].loc.y, particleData[i].loc.z);         
    }

   if (!user) {
    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorInitialize - cp9.6F - Error: user is NULL\n");
    return PETSC_ERR_ARG_NULL;
    }

    // Set values in the PETSc vector
    for (PetscInt i = 0; i < numParticles; ++i) {

      if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD," ParticleVectorInitialize - cp9.6S%d - user: %p - particleData[i]: %p - i: %d \n",i+1,user,&particleData[i],i);
       
      
      ParticleVectorInsert(user,i, &particleData[i]);
    }


   if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorInitialize - cp9.7 - ParticleVector Randomized \n");

    // Assemble the vector
    VecAssemblyBegin(user->ParticleVec);
    VecAssemblyEnd(user->ParticleVec);

 if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorInitialize - cp9.8 - ParticleVector Assembled \n");

    // Clean up
    free(particleData);
    PetscRandomDestroy(&rand);

    // Diagnostics 
    PetscInt low,high;
    PetscInt rank;
    Particle sample;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Get the ownership range for this rank
    ierr = VecGetOwnershipRange(user->ParticleVec, &low, &high); CHKERRQ(ierr);

    // Synchronized output to display ownership information
    if(visflg==ctr){
      ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,
				   "ParticleVectorInitialize - cp9.9%d - Rank %d owns indices from %d to %d\n", rank+1,rank, (int)low, (int)(high - 1)); CHKERRQ(ierr);
      ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
    }
    // if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorInitialize - cp9.9 - particleData and rand freed \n");

    // Output sample location
    // PetscInt range = high - low;
    // PetscInt PIdx = low
    // PetscInt locIdx = PIdx + 3; 
    // ierr = VecGetValues(
    // Synchronized output to display ownership information
    // if(visflg==ctr){
    //  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.9%d - Rank %d owns indices from %d to %d\n", rank+1,rank, (int)low, (int)(high - 1)); CHKERRQ(ierr);
// ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
//  }
    return(0);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*

//Quick test to determine the processor of a particle based on it's location //////

PetscErrorCode  WriteBoundingBox(UserCtx *user, BoundingBox *bbox){

  PetscInt i,j,k,rank;
  PetscErrorCode ierr;
  PetscInt xs,ys,zs,xe,ye,ze;
  DMDALocalInfo info;
  Vec Coor;
  Cmpnts ***coor,min_coords,max_coords;
  PetscInt ctr = 10;

  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.31 - user - %p \n",user);

  DM da = user->da;
  DM fda = user->fda; 

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

  
  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.32 - da - %p - user->da - %p \n",da,user->da);

  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.33 - fda - %p - user->fda - %p \n",fda,user->fda);
  
  ierr = DMGetCoordinatesLocal(user->da, &Coor); CHKERRQ(ierr);

if (Coor == NULL) {
    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp 10.34F - Error: Coor vector is NULL.\n");
    return 1; 
}

  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.34S - coor - %p - Coor - %p \n",Coor,coor);
 
  ierr = DMDAVecGetArrayRead(user->fda,Coor,&coor); CHKERRQ(ierr);

  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.35 - fda - %p - Coor - %p - coor - %p \n",user->fda,Coor,coor);

  ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  if(visflg==ctr) PetscPrintf(PETSC_COMM_SELF, " WriteBoundingBox - cp10.36 - rank - %d; xs - %d; xe - %d; ys - %d; ye - %d; zs - %d; ze - %d \n",rank,xs,xe,ys,ye,zs,ze);
   // Initialize min and max coordinates
  min_coords.x = min_coords.y = min_coords.z =  PETSC_MAX_REAL;
  max_coords.x = max_coords.y = max_coords.z =  PETSC_MIN_REAL;

  // Find min and max coordinates for the local grid
  for (k = zs; k < ze; k++) {
     for (j = ys; j < ye ; j++) {
         for (i = xs; i < xe; i++) {

                // Access numerical values from the coor array using . operator

                min_coords.x = PetscMin(min_coords.x,coor[k][j][i].x);
                min_coords.y = PetscMin(min_coords.y,coor[k][j][i].y);
                min_coords.z = PetscMin(min_coords.z,coor[k][j][i].z);                

                max_coords.x = PetscMax(max_coords.x,coor[k][j][i].x);
                max_coords.y = PetscMax(max_coords.y,coor[k][j][i].y);
                max_coords.z = PetscMax(max_coords.z,coor[k][j][i].z);

		//	if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.36%d - i,j,k - %d,%d,%d; max_coords - x,y,z - %f, %f, %f; min_coords - x,y,z - %f,%f,%f \n",(i+1)*(j+1)*(k+1),i,j,k,max_coords.x,max_coords.y,max_coords.z,min_coords.x,min_coords.y,min_coords.z);
            }
        }
    }
    ierr = DMDAVecRestoreArrayRead(fda,Coor,&coor); CHKERRQ(ierr);
    
    if(visflg==ctr) PetscPrintf(PETSC_COMM_SELF, " WriteBoundingBox - cp10.37%d - rank - %d - max_coords - x,y,z - %f, %f, %f; min_coords - x,y,z - %f,%f,%f \n",rank+1,rank,max_coords.x,max_coords.y,max_coords.z,min_coords.x,min_coords.y,min_coords.z);

    
    bbox->max_coords = max_coords;
    bbox->min_coords = min_coords;

    PetscBarrier(PETSC_NULL);

    return 0;

}
*/


// Function to create bounding boxes for processors and gather them on rank 0

PetscErrorCode WriteBoundingBox(UserCtx *user, BoundingBox **all_bboxes) {
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

    if (visflg == ctr) PetscSynchronizedPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.31 - user - %p \n", user);
    PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

    DM da = user->da;
    DM fda = user->fda;

    //  if (visflg == ctr) PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.32 - rank - da - %p - user->da - %p \n", da, user->da);
    //   if (visflg == ctr) PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.33 - fda - %p - user->fda - %p \n", fda, user->fda);

    PetscBarrier(PETSC_NULL);

    ierr = DMGetCoordinatesLocal(user->da, &Coor); CHKERRQ(ierr);

    if (Coor == NULL) {
        if (visflg == ctr) PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.34F - Error: Coor vector is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Coordinates vector is NULL");
    }

    if (visflg == ctr) PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.34S - Coor - %p \n", Coor);

    ierr = DMDAVecGetArrayRead(user->fda, Coor, &coor); CHKERRQ(ierr);

    if (visflg == ctr) PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.35 - fda - %p - Coor - %p - coor - %p \n", user->fda, Coor, coor);

    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

    xs = info.xs; xe = xs + info.xm;
    ys = info.ys; ye = ys + info.ym;
    zs = info.zs; ze = zs + info.zm;

    if (visflg == ctr) PetscSynchronizedPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.36 - rank - %d; xs - %d; xe - %d; ys - %d; ye - %d; zs - %d; ze - %d \n",
                                   rank, xs, xe, ys, ye, zs, ze);
    PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

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

    // if (visflg == ctr) PetscPrintf(PETSC_COMM_SELF, " WriteBoundingBox - cp10.37 - rank - %d - max_coords: x=%f, y=%f, z=%f; min_coords: x=%f, y=%f, z=%f \n",
    //    rank, max_coords.x, max_coords.y, max_coords.z, min_coords.x, min_coords.y, min_coords.z);

    // Create local bounding box
    BoundingBox local_bbox;
    local_bbox.min_coords = min_coords;
    local_bbox.max_coords = max_coords;

    if(visflg==ctr){
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "WriteBoundingBox - cp10.38 - rank - %d; bbox: min.x - %f, min.y - %f,min.z - %f; max.x - %f, max.y - %f, max.z - %f; \n",rank,local_bbox.min_coords.x,local_bbox.min_coords.y,local_bbox.min_coords.z,local_bbox.max_coords.x,local_bbox.max_coords.y,local_bbox.max_coords.z);
      PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
    }// visflg

    user->bbox = local_bbox; // update bounding box inside local UserCtx;
    // Allocate memory on rank 0 for all bounding boxes1
    BoundingBox *bbox_array = NULL;
    if (rank == 0) {
        bbox_array = (BoundingBox *)malloc(size * sizeof(BoundingBox));
        if (!bbox_array) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_MEM, "WriteBoundingBox -cp10.39F - Memory allocation failed for bounding box array on rank 0");
        }
       if (visflg == ctr) PetscPrintf(PETSC_COMM_SELF, " WriteBoundingBox - cp10.39S - rank - %d - bbox_array %p \n",rank,bbox_array);

    }// rank 0

    // Gather local bounding boxes from all processors to rank 0
    ierr = MPI_Gather(&local_bbox, sizeof(BoundingBox), MPI_BYTE,
                      bbox_array, sizeof(BoundingBox), MPI_BYTE,
                      0, PETSC_COMM_WORLD); CHKERRQ(ierr);

    // On rank 0, assign gathered bounding boxes to output pointer
    if (rank == 0) {
        *all_bboxes = bbox_array;

        // Print all gathered bounding boxes for verification
        for (PetscInt i = 0; i < size; i++) {
	  if(visflg==ctr) { PetscPrintf(PETSC_COMM_WORLD, "WriteBoundingBox - 10.40%d - Gathered Data - Host Rank %d - Receiver Rank - %d : coords: x=[%f,%f] y=[%f,%f] z=[%f,%f] \n", i+1,i,rank,bbox_array[i].min_coords.x, bbox_array[i].max_coords.x,bbox_array[i].min_coords.y,bbox_array[i].max_coords.y,bbox_array[i].min_coords.z, bbox_array[i].max_coords.z);
	  }
        }
    }else *all_bboxes = NULL;  // rank 0 
   



    PetscBarrier(PETSC_NULL);

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

PetscBool  CPUPointIntersectCheck(BoundingBox *bbox, Particle *particle){
   
   Cmpnts loc,min_coords,max_coords;
   PetscErrorCode ierr;
   PetscInt ctr = 10;
   PetscInt rank;
   loc = particle->loc;
   PetscBool Intersects = PETSC_FALSE;
   min_coords = bbox->min_coords;
   max_coords = bbox->max_coords;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
   
    if(visflg==ctr) PetscPrintf(PETSC_COMM_SELF, "CPUPointIntersectCheck - cp10.7V1 - rank - %d particle - %p ; loc.x - %f, bbox: min_x - %f | min_y - %f | min_z - %f |  max_x - %f | max_y - %f | max_z - %f \n",rank,particle,loc.x,min_coords.x,min_coords.y,min_coords.z,max_coords.x,max_coords.y,max_coords.z);

    if ((loc.x >= min_coords.x && loc.x <= max_coords.x) &&
        (loc.y >= min_coords.y && loc.y <= max_coords.y) &&
        (loc.z >= min_coords.z && loc.z <= max_coords.z)) {
        Intersects = PETSC_TRUE;
    }
   return Intersects;
}

///////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// Walking Search within a rank //////////////////////////////////////////

PetscReal distance_search(Cmpnts p1, Cmpnts p2, Cmpnts p3, Cmpnts p4, Cmpnts p, PetscReal *d)
{
  PetscReal xn1, yn1, zn1;
  PetscReal xc, yc, zc;
  
  PetscReal dx1, dy1, dz1, dx2, dy2, dz2, r;

  dx1 = p3.x - p1.x;
  dy1 = p3.y - p1.y;
  dz1 = p3.z - p1.z;

  dx2 = p4.x - p2.x;
  dy2 = p4.y - p2.y;
  dz2 = p4.z - p2.z;

  xn1 = dy1 * dz2 - dz1 * dy2;
  yn1 = - (dx1 * dz2 - dz1 * dx2);
  zn1 = dx1 * dy2 - dy1 * dx2;

  r = sqrt(xn1 * xn1 + yn1 * yn1 + zn1 * zn1);
  xn1 /= r; yn1 /= r; zn1 /= r;

  xc = 0.25 * (p1.x + p2.x + p3.x + p4.x);
  yc = 0.25 * (p1.y + p2.y + p3.y + p4.y);
  zc = 0.25 * (p1.z + p2.z + p3.z + p4.z);

  *d = (p.x - xc) * xn1 + (p.y - yc) * yn1 + (p.z - zc) * zn1;
  if (PetscAbsReal(*d)<1.e-6) *d=0.;
  return (0);
}

PetscErrorCode CellPointInterceptCalculate(Cmpnts p, Cmpnts cell[8], PetscReal d[6])
{ 
  // d is an empty array of zeros that will be populated with intercepts from cell faces to determine whether the        particle is inside the cel or not, and if not, what direction to go to for the next cell.

  // k direction
  distance_search(cell[0], cell[1], cell[2], cell[3], p, &(d[4]));
  distance_search(cell[4], cell[7], cell[6], cell[5], p, &(d[5]));

  // j direction
  distance_search(cell[0], cell[4], cell[5], cell[1], p, &(d[2]));
  distance_search(cell[3], cell[2], cell[6], cell[7], p, &(d[3]));

  // i direction
  distance_search(cell[0], cell[3], cell[7], cell[4], p, &(d[0]));
  distance_search(cell[1], cell[5], cell[6], cell[2], p, &(d[1]));
  
  return(0);
}

PetscErrorCode CellValuesGet(Cmpnts ***coor, PetscInt idx, PetscInt idy, PetscInt idz, Cmpnts cell[8])
{
  cell[0] = coor[idz][idy][idx];
  cell[1] = coor[idz][idy][idx+1];
  cell[2] = coor[idz+1][idy][idx+1];
  cell[3] = coor[idz+1][idy][idx];
  cell[4] = coor[idz][idy+1][idx];
  cell[5] = coor[idz][idy+1][idx+1];
  cell[6] = coor[idz+1][idy+1][idx+1];
  cell[7] = coor[idz+1][idy+1][idx];

  return(0);
}

PetscBool CellIntersectCheck(PetscReal d[6])
{
  PetscBool Intersects = PETSC_FALSE;
  for(int i = 0; i<6; i++){
    if(d[i] <= 0.0) {
       Intersects = PETSC_TRUE;
       break;
    }
  }
  return Intersects;
}


PetscErrorCode  InterpolationWeightsCalculate(Cmpnts a, PetscReal d[6])
{
  a.x = d[0]/(d[0]+d[1]);
  a.y = d[2]/(d[2]+d[3]);
  a.z = d[4]/(d[4]+d[5]);
  return(0);
}

PetscErrorCode CellIndexUpdate(PetscReal *d, PetscInt *idx, PetscInt *idy, PetscInt *idz)
{
  // i direction
  if(d[0] || d[1] <0){
    if(d[0]<0) idx--;
    else idx++;
    }
 
  // j direction
  if(d[2] || d[3] <0){
    if(d[2]<0) idy--;
    else idy++;
    }

  // k direction 
  if(d[4] || d[5] <0){
    if(d[4]<0) idz--;
    else idz++;  
    }
  return(0);
}



PetscErrorCode WalkingSearch(UserCtx *user, Particle *particle)
{
  PetscInt i,j,k;
  PetscInt idx,idy,idz;
  PetscInt xs,ys,zs,xe,ye,ze;
  PetscInt mz, my, mx;
  PetscInt lxs, lxe, lys, lye, lzs, lze;
  DM da = user->da;
  DM fda = user->fda; 
  DMDALocalInfo info;
  Vec Coor;
  Cmpnts ***coor, p, cell[8];
  PetscReal d[6];
  PetscBool Cell_found = PETSC_FALSE;  
  
  DMGetCoordinatesLocal(da, &Coor); 
  DMDAVecGetArrayRead(fda,Coor,&coor);

  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  
  lxs = xs-1; lxe = xe+1;
  lys = ys-1; lye = ye+1;
  lzs = zs-1; lze = ze+1;

  if (xs==0) lxs = xs;
  if (ys==0) lys = ys;
  if (zs==0) lzs = zs;

  if (xe==mx) lxe=xe;
  if (ye==my) lye=ye;
  if (ze==mz) lze=ze;
  
  idx = lxs;
  idy = lys;
  idz = lzs;
  p = particle->loc;

  //  cell = (Cmpnts *)malloc(8 * sizeof(Cmpnts));

  while (Cell_found== PETSC_FALSE){
    CellValuesGet(coor,idx,idy,idz,cell);  // Get the coordinates for each corner of the cell 
    CellPointInterceptCalculate(p,cell,d); // Calculate the intercepts/distances of the particle from the faces of cell  
    Cell_found = CellIntersectCheck(d); // Check if the cell intersects the particle, if yes, update the flag cell_found.
    if(!Cell_found)  CellIndexUpdate(&d,&idx,&idy,&idz);  // Update the IDs of the cell to check for intersection (Walking update)
  }
  
  InterpolationWeightsCalculate(particle->weights,d); // If the cell enclosing the particle is found, update the interpolation weights of the particle.
  
  


  // particle->cell[0] = idx; particle->cell[1] = idy; particle->cell[2] = idz; // Update the particle location by identifying global indices of the bottom front left corner vertex of the cell that houses the particle.
   
  DMDAVecRestoreArrayRead(fda, Coor, &coor);
  // free(cell);
  return (0);
}

// Function to create bounding boxes for processors, and locate particles in the grid and calculate interpolation weights for each particle.

PetscErrorCode ParticlesLocate(UserCtx *user, PetscInt numParticles) {

    PetscReal *particledata;
    PetscInt ctr = 10;
    PetscInt rank;
    PetscErrorCode ierr;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    
    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - cp10.1 - user - %p \n",user);

    PetscInt particleSize = sizeof(Particle) / sizeof(PetscScalar);

    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - cp10.2 - particleSize - %d - numParticles - %d \n",particleSize,numParticles);

    BoundingBox *bboxlist; // Array of bboxes for each processor.
    
    // Create bounding boxes for all processors

    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - cp10.3 - bboxlist %p  \n",&bboxlist);

    WriteBoundingBox(user,&bboxlist);
    
    PetscBarrier(PETSC_NULL);

 ierr = ParticleSwarmCreate(user, np); CHKERRQ(ierr);

 if (visflg == ctr) PetscPrintf(PETSC_COMM_WORLD, "main - Particle Swarm Created and Initialized - cp9 \n");

  ParticleSwarmViewPositions(user);

 

    /* // Get the array of particles from the vector */
    /* ierr = VecGetArray(user->ParticleVec, &particledata); CHKERRQ(ierr); */

    /* if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - cp10.5 - particledata array - %p - size - %\n",particledata); */


    /* PetscReal *pdata; */
    /* Particle particle; */
     
    /* // Loop through each particle */
    /* for (PetscInt i = 0; i < numParticles; i++) { */
   
    /*   pdata = &particledata[i * particleSize]; // this is a pointer to a specific part of particledata - associated with the current particle */

    /*   // Safely copy data into the Particle structure */
    /*   memcpy(&particle, pdata, sizeof(Particle)); */
      
    /*   if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - cp10.7B%d - particle - [%d] - %p ; particle.loc.x - %f,\n",i+1,i,&particle,particle.loc.x); */
    /*   if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - cp10.8B%d - particle - [%d] - %p ; particle.cell[0,1,2] - [%f,%f,%f] \n",i+1,i,&particle,particle.cell[0],particle.cell[1],particle.cell[2]); */
        
    /*     if (CPUPointIntersectCheck(bboxlist, &particle)) { */
    /*         WalkingSearch(user, &particle); */
    /*     } */

    /*  // Copy back into pdata */
    /*  memcpy(pdata, &particle, sizeof(Particle)); */


    /*   if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - cp10.7A%d - particle - [%d] - %p ; particle.loc.x - %f,\n",i+1,i,&particle,particle.loc.x); */
    /*   if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - cp10.8A%d - particle - [%d] - %p ; particle.cell[0,1,2] - [%f,%f,%f] \n",i+1,i,&particle,particle.cell[0],particle.cell[1],particle.cell[2]); */
      
    /* } */

    /* // Restore the array */
    /* VecRestoreArray(user->ParticleVec,&particledata); */

    if(rank==0) free(bboxlist);

    return(0);
}


PetscErrorCode InterpolationCoefficientsCalculate(Cmpnts a,PetscReal coeffs[8])
{
  PetscReal a1,a2,a3;
  a1 = a.x;
  a2 = a.y;
  a3 = a.z;
  
  coeffs[0] = a1*a2*a3;
  coeffs[1] = (a1-1.0)*a2*a3;
  coeffs[2] = a1*(a2-1.0)*a3;
  coeffs[3] = (a1-1.0)*(a2-1.0)*a3;
  coeffs[4] = a1*a2*(a3-1.0);
  coeffs[5] = (a1-1.0)*a2*(a3-1.0);
  coeffs[6] = a1*(a2-1.0)*a3;
  coeffs[7] = a1*(a2-1.0)*(a3-1.0);
  coeffs[8] = (a1-1.0)*(a2-1.0)*(a3-1.0);
  
  return(0);
}


PetscErrorCode InterpolationMatrixInitialize(PetscInt numParticles)
{
  
  PetscInt rows, cols;
  Mat InterpMat;
  
  PetscReal coeffs[8];
 
  rows = numParticles;
  cols =  8 * numParticles;
  
  // Create Matrix
  MatCreate(PETSC_COMM_WORLD, InterpMat);
  MatSetSizes(InterpMat, PETSC_DECIDE,PETSC_DECIDE,rows,cols);
  MatSetType(InterpMat, MATMPIAIJ);
  MatSetUp(InterpMat);
  MatZeroEntries(InterpMat);  

  return(0);
}


PetscErrorCode InterpolationMatrixPopulate(Mat InterpMat, Vec particleVec,PetscInt numParticles)
{
  Particle *particles;
  PetscReal coeffs[8];
  PetscInt row;
  PetscInt colIndices[8];
  
  // Get the array of particles from the vector
  VecGetArray(particleVec, (PetscScalar**)&particles);
  
  for(int p = 0; p < numParticles; p++){
      
      row = p;  // Each particle is associated with a row in the matrix. 
   
      for(int q = 0; q<8; q++){
         
          colIndices[q] = p*8 + q; // Diagonal block matrix 
      
      }
      
      InterpolationCoefficientsCalculate(particles[p].weights,coeffs); // Calculate the interpolation coefficients for the particle using it's weights.
   
      MatSetValues(InterpMat,1, &row, 8, colIndices,coeffs,INSERT_VALUES); 
     
  }

  VecRestoreArray(particleVec, (PetscScalar**)&particles);

  // Assembly of the matrix.
  MatAssemblyBegin(InterpMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(InterpMat, MAT_FINAL_ASSEMBLY);
 
  return(0);
}


PetscErrorCode InterpolationVectorCreate(UserCtx* user, Vec particleVec, PetscInt numParticles)
{
  Particle *particles; 

  PetscInt Vecsize = numParticles*8, indices[8],startid;
  
  PetscInt CmpntsSize = sizeof(Cmpnts) / sizeof(PetscReal);
  
  Vec Ucat = user->Ucat, InterpVec;

  DM  da = user->da, fda = user->fda;
  DMDALocalInfo info = user->info;

  Cmpnts ***ucat,host_vel[8];
  
  VecCreate(MPI_COMM_WORLD, &InterpVec);
  VecSetSizes(InterpVec, PETSC_DECIDE,Vecsize*CmpntsSize);
  VecSetFromOptions(InterpVec);
  
  VecSet(InterpVec, 0.0);
  

  DMDAVecGetArray(fda,Ucat,&ucat);

  // Get the array of particles from the vector
  VecGetArray(particleVec, (PetscScalar**)&particles);

  for(int i = 0 ; i < numParticles; i++){
    
    startid = i;
      
    for( int j = 0; j < 8; j++){

        indices[j] = startid + j;
            
    }
 
      CellValuesGet(ucat,particles[i].cell[0], particles[i].cell[1], particles[i].cell[2],host_vel);
      VecSetValues(InterpVec,8,indices, host_vel,INSERT_VALUES);
  }
  
  
  DMDAVecRestoreArray(fda, Ucat, &ucat);
  
  VecRestoreArray(particleVec, (PetscScalar**)&particles);
  
  return(0);
}

PetscErrorCode Interpolate(Mat InterpMat, Vec InterpVec, PetscInt numParticles)
{

  Vec R;
  PetscInt Vecsize = numParticles;
  
  PetscInt CmpntsSize = sizeof(Cmpnts)/sizeof(PetscReal);
  
  
  VecCreate(MPI_COMM_WORLD,&R);
  VecSetSizes(R,PETSC_DECIDE,Vecsize*CmpntsSize);
  VecSetFromOptions(R);
  VecSet(R,0.0);
 
  MatMult(InterpMat, InterpVec, R);
  
  return(0);

}

PetscErrorCode ParticleVelocityUpdate(Vec ParticleVec, Vec R, PetscInt numParticles)
{
 
  Particle *particles;
  Cmpnts *resultant;
  PetscInt particleSize = sizeof(Particle) / sizeof(PetscScalar);

  VecGetArray(ParticleVec, (PetscScalar**)&particles);
  VecGetArray(R, (PetscScalar**)&resultant);
  
  for(PetscInt i = 0; i<numParticles; i++){
  
    particles[i].vel.x = resultant[i].x;
    particles[i].vel.y = resultant[i].y;
    particles[i].vel.z = resultant[i].z;
  }
  
  
  VecRestoreArray(ParticleVec, (PetscScalar**)&particles);
  VecRestoreArray(R, (PetscScalar**)&resultant);
    
  return(0);
}

PetscErrorCode InterpolationMatrixDestroy(Mat InterpMat)
{
  
  MatDestroy(&InterpMat);

  return(0);
}


PetscErrorCode InterpolationVectorDestroy(Vec InterpVec, Vec R)
{

  VecDestroy(&InterpVec);
  VecDestroy(&R);
  
  return(0);

}

// Read Velocity field - ufield //

PetscErrorCode Ucat_Binary_Input(UserCtx *user)
{
  PetscViewer viewer;
  char filen[90];
  PetscInt bi=user->_this;
  PetscInt ctr = 7;
  
  sprintf(filen, "results/ufield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);

  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "Ucat_Binary_Input - cp7.2 - user - %p  \n",user);

  PetscInt N;

  VecGetSize(user->Ucat, &N);
  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "Ucat_Binary_Input - cp7.3  SizeOf(ucat) - %d \n", N);
  VecLoad((user->Ucat),viewer);
 
  PetscViewerDestroy(&viewer);

  PetscBarrier(PETSC_NULL);

  return 0;

}

PetscErrorCode ReadCoordinates(UserCtx *user) {

  Cmpnts ***coor;
  PetscInt ctr = 0;
  Vec Coor;
  PetscInt bi, i, j, k, rank, IM, JM, KM;
  PetscReal *gc;
  FILE *fd;
  PetscReal	d0 = 1.;
  PetscInt    generate_grid=0, grid1d=0, nblk=block_number;

  PetscOptionsInsertFile(PETSC_COMM_WORLD,PETSC_NULL, "control.dat", PETSC_TRUE); 

  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ReadCoordinates - cp5.1 - user - %p \n", user);

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ReadCoordinates - cp5.2 - rank - %d \n", rank);

  PetscOptionsGetInt(PETSC_NULL,PETSC_NULL, "-grid", &generate_grid, PETSC_NULL);

  PetscReal	cl = 1.;
  PetscOptionsGetReal(PETSC_NULL,PETSC_NULL, "-chact_leng", &cl, PETSC_NULL);

  PetscReal L_x,L_y,L_z;


  if (generate_grid) {
    PetscOptionsGetReal(PETSC_NULL,PETSC_NULL, "-L_x", &L_x, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL,PETSC_NULL, "-L_y", &L_y, PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL,PETSC_NULL, "-L_z", &L_z, PETSC_NULL);

    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.3 - ksi eta zeta -  %le %le %le \n",L_x,L_y,L_z);
    //  block_number=1;
  } else {
    if (!rank) {
      fd = fopen("grid.dat", "r");
      fscanf(fd, "%i\n", &block_number);
      MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }
  }

  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ReadCoordinates - cp5.4 - rank - %d \n", rank);
 
  PetscInt imm[block_number], kmm[block_number], jmm[block_number];
  if (generate_grid) {
    PetscOptionsGetIntArray(PETSC_NULL,PETSC_NULL, "-im", imm, &nblk, PETSC_NULL);
    PetscOptionsGetIntArray(PETSC_NULL,PETSC_NULL, "-jm", jmm, &nblk, PETSC_NULL);
    PetscOptionsGetIntArray(PETSC_NULL,PETSC_NULL, "-km", kmm, &nblk, PETSC_NULL);
  }

  if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "ReadCoordinates - cp5.5 -  imm[0] - %d; jmm[0] - %d; kmm[0] - %d  \n",imm[0],jmm[0],kmm[0]);  

  for (bi=0; bi<block_number; bi++) {

    //if (bi==0)  cl = 1.;
    //if (bi==1)  cl = 0.85;

    if (!rank) {
      if (!generate_grid)
	fscanf(fd, "%i %i %i\n", &(user[bi].IM), &(user[bi].JM), &(user[bi].KM));
      
      else {
	user[bi].IM=imm[bi];
	user[bi].JM=jmm[bi];
	user[bi].KM=kmm[bi];
      }	
      IM = user[bi].IM; JM = user[bi].JM; KM = user[bi].KM;
      
      MPI_Bcast(&(user[bi].IM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].JM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].KM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(&(user[bi].IM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].JM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].KM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      
      IM = user[bi].IM; JM = user[bi].JM; KM = user[bi].KM;
    }
    
    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.6 - IM - %d; JM - %d; KM - %d \n",IM,JM,KM);


    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.61 - user[bi].da - %p;\n",(user[bi].da));

    PetscErrorCode ierr;

    ierr = DMDACreate3d(PETSC_COMM_WORLD,
		    DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,  // Boundary types 
		    DMDA_STENCIL_BOX,                                      // Stencil type
		    user[bi].IM + 1, user[bi].JM + 1, user[bi].KM + 1,     // Global Dimensions (M,N,P)
                    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,              // Dimensions per process (m, n, p)
                    1,                                                     // dof
                    1,                                                     // s (stencil width)
                    NULL, NULL, NULL,                                      // lx[], ly[], lz[]
			&(user[bi].da)); CHKERRQ(ierr);                    // &da 
  

    ierr = DMSetUp(user[bi].da); CHKERRQ(ierr);
   
    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.7 - user[bi].da - %p;\n",user[bi].da);

    
    ierr = DMDASetUniformCoordinates(user[bi].da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0); CHKERRQ(ierr);

    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.71 - user[bi].fda - %p;\n",(user[bi].fda));

    if(user[bi].fda == NULL) if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.72F - fda is NULL before DMGetCoordinateDM \n");

    ierr = DMGetCoordinateDM(user[bi].da, &(user[bi].fda)); CHKERRQ(ierr);

    if(user[bi].fda == NULL) if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.73F - fda is NULL after DMGetCoordinateDM  \n");

    if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.8 - user[bi].fda - %p;\n",(user[bi].fda));

    ierr = DMView(user[bi].fda, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    DMDAGetLocalInfo(user[bi].da, &(user[bi].info));
    
    DMDALocalInfo	info = user[bi].info;

    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    //  PetscInt	gmx = info.gmx, gmy = info.gmy, gmz = info.gmz;
  
    PetscOptionsGetInt(PETSC_NULL,PETSC_NULL, "-grid", &generate_grid, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL,PETSC_NULL, "-grid1d", &grid1d, PETSC_NULL);

    if (grid1d) PetscMalloc((IM+JM+KM)*sizeof(PetscReal), &gc);

    else        PetscMalloc(3*(IM*JM*KM)*sizeof(PetscReal), &gc);

    DMGetCoordinatesLocal(user[bi].da, &Coor);
    DMDAVecGetArray(user[bi].fda, Coor, &coor);

    if (!rank) {
      if (grid1d) {
	PetscReal xx;
	// read i
	for (i=0; i<IM; i++) 
	  fscanf(fd, "%le %le %le\n",&gc[i],&xx,&xx);
	// read j
	for (j=0; j<JM; j++) 
	  fscanf(fd, "%le %le %le\n",&xx,&gc[IM+j],&xx);
	// read k
	for (i=0; i<KM; i++) 
	  fscanf(fd, "%le %le %le\n",&xx,&xx,&gc[IM+JM+i]);
	
	MPI_Bcast(gc, (IM+JM+KM), MPIU_REAL, 0, PETSC_COMM_WORLD);
	
	for (k=zs; k<=ze; k++) {
	  for (j=ys; j<=ye; j++) {
	    for (i=xs; i<=xe; i++) {
	      if (k<KM && j<JM && i<IM) {
		coor[k][j][i].x = *(gc + i)/cl*L_dim;
		coor[k][j][i].y = *(gc + IM + j)/cl*L_dim;
		coor[k][j][i].z = *(gc + IM + JM + k)/cl*L_dim;
	      }
	    }
	  }
	} // grid1d
	
      } else { // if 3d gridgen file
	for (k=0; k<KM; k++) {
	  for (j=0; j<JM; j++) {
	    for (i=0; i<IM; i++) {
	      if (!generate_grid)
		fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3);
	      else
		*(gc+(k*JM*IM + j*IM + i)*3) = L_x/(IM) * i;
	    }
	  }
	}
	
	for (k=0; k<KM; k++) {
	  for (j=0; j<JM; j++) {
	    for (i=0; i<IM; i++) {
	      if (!generate_grid)
		fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 1);
	      else
		*(gc+(k*JM*IM + j*IM + i)*3+1) = L_y/(JM) * j;
	    }
	  }
	}
	
	for (k=0; k<KM; k++) {
	  for (j=0; j<JM; j++) {
	    for (i=0; i<IM; i++) {
	      if (!generate_grid)
		fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 2);
	    else
	      *(gc+(k*JM*IM + j*IM + i)*3+2) = L_z/(KM) * k;
	    }
	  }
	}	
	
	MPI_Bcast(gc, 3*(IM*JM*KM), MPIU_REAL, 0, PETSC_COMM_WORLD);
      
	for (k=zs; k<=ze; k++) {
	  for (j=ys; j<=ye; j++) {
	    for (i=xs; i<=xe; i++) {
	      if (k<KM && j<JM && i<IM) {
		coor[k][j][i].x = *(gc + (k * (IM*JM) + j * IM + i) * 3  )/cl;
		coor[k][j][i].y = *(gc + (k * (IM*JM) + j * IM + i) * 3+1)/cl;
		coor[k][j][i].z = *(gc + (k * (IM*JM) + j * IM + i) * 3+2)/cl;
	      }
	    }
	  }
	}
      }  
    }
    else {
      if (grid1d) {
	MPI_Bcast(gc, (IM+JM+KM), MPIU_REAL, 0, PETSC_COMM_WORLD);
	
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    for (i=xs; i<xe; i++) {
	      if (k<KM && j<JM && i<IM) {
		coor[k][j][i].x = *(gc + i)/cl*L_dim;
		coor[k][j][i].y = *(gc + IM + j)/cl*L_dim;
		coor[k][j][i].z = *(gc + IM + JM + k)/cl*L_dim;
	      }
	    }
	  }
	}

      } else { // if 3d gridgen file

	MPI_Bcast(gc, 3*(IM*JM*KM), MPIU_REAL, 0, PETSC_COMM_WORLD);
	
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    for (i=xs; i<xe; i++) {
	      if (k<KM && j<JM && i<IM) {
		coor[k][j][i].x = *(gc + (k * (IM*JM) + j * IM + i) * 3  )/cl;
		coor[k][j][i].y = *(gc + (k * (IM*JM) + j * IM + i) * 3+1)/cl;
		coor[k][j][i].z = *(gc + (k * (IM*JM) + j * IM + i) * 3+2)/cl;
	      }
	    }
	  }
	}
      }
    }
    PetscFree(gc);
    DMDAVecRestoreArray(user[bi].fda, Coor, &coor);
    
    Vec	gCoor;
    DMGetCoordinates(user[bi].da, &gCoor);
    DMLocalToGlobalBegin(user[bi].fda, Coor, INSERT_VALUES, gCoor);
    DMLocalToGlobalEnd(user[bi].fda, Coor, INSERT_VALUES, gCoor);
    
    DMGlobalToLocalBegin(user[bi].fda, gCoor, INSERT_VALUES, Coor);
    DMGlobalToLocalEnd(user[bi].fda, gCoor, INSERT_VALUES, Coor);
    
    //   VecDestroy(&gCoor);
  }
  
  if (!rank) {
    if(!generate_grid)
      fclose(fd);
  }
  
  for (bi=0; bi<block_number; bi++) {
    user[bi]._this = bi;
    if(visflg==ctr) PetscSynchronizedPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.9%d - rank - %d, user[bi] - %p \n",bi+1,rank,&user[bi]);
  }
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
  // VecDestroy(&Coor);
  
  return(0);
}

////////////////////////////////////////////////


/*
PetscErrorCode InitializeParticles(UserCtx* user, PetscInt numParticles){
  
  PetscErrorCode ierr;
  PetscReal *coors;
  PetscInt localnumParticles = 0;
  PetscInt rank,size;
  PetscRandom randx,randy,randz;
  PetscReal Lx,Ly,Lz;
  
  PetscOptionsInsertFile(PETSC_COMM_WORLD, PETSC_NULL, "control.dat", PETSC_TRUE);

  PetscOptionsGetReal(PETSC_NULL,PETSC_NULL, "-L_x", &Lx, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL,PETSC_NULL, "-L_y", &Ly, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL,PETSC_NULL, "-L_z", &Lz, PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD," Lx - %f; Ly - %f; Lz - %f \n",Lx,Ly,Lz);
  
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
  
  localnumParticles = numParticles/size; 
 
  ierr = DMSwarmGetField(user->swarm, "DMSwarmPIC_coor",NULL,NULL,(void**)&coors); CHKERRQ(ierr);
 
  // Initialize random number generators
  // x 
  ierr = PetscRandomCreate(PETSC_COMM_WORLD, &randx); CHKERRQ(ierr);
  ierr = PetscRandomSetType(randx, PETSCRAND48); CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(randx, 0.0, Lx); CHKERRQ(ierr);
  ierr = PetscRandomSeed(randx); CHKERRQ(ierr);
  // y
  ierr = PetscRandomCreate(PETSC_COMM_WORLD, &randy); CHKERRQ(ierr);
  ierr = PetscRandomSetType(randy, PETSCRAND48); CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(randy, 0.0, Ly); CHKERRQ(ierr);
  ierr = PetscRandomSeed(randy); CHKERRQ(ierr);
  // z
  ierr = PetscRandomCreate(PETSC_COMM_WORLD, &randz); CHKERRQ(ierr);
  ierr = PetscRandomSetType(randz, PETSCRAND48); CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(randz, 0.0, Lz); CHKERRQ(ierr);
  ierr = PetscRandomSeed(randz); CHKERRQ(ierr);

  for(PetscInt p = 0; p < localnumParticles; p ++ ){
         
    //    coors[p*3 + 0] = (0.5+rank)*Lx;
    //        coors[p*3 + 1] = (0.5*rank)*Ly;            
    //        coors[p*3 + 2] = (0.5*rank)*Lz;
    
            ierr = PetscRandomGetValue(randx, &coors[p * 3 + 0]); CHKERRQ(ierr);
	    ierr = PetscRandomGetValue(randy, &coors[p * 3 + 1]); CHKERRQ(ierr);
	    ierr = PetscRandomGetValue(randz, &coors[p * 3 + 2]); CHKERRQ(ierr);
    
   }
        ierr = DMSwarmRestoreField(user->swarm, "DMSwarmPIC_coor",NULL,NULL,(void**)&coors); CHKERRQ(ierr);
       
        ierr = PetscRandomDestroy(&randx); CHKERRQ(ierr);
        ierr = PetscRandomDestroy(&randy); CHKERRQ(ierr);
        ierr = PetscRandomDestroy(&randz); CHKERRQ(ierr);
    return 0;
}
*/

PetscErrorCode ParticleSwarmCreate(UserCtx *user, PetscInt numParticles) {
    PetscErrorCode ierr;
    PetscInt rank;
    PetscInt size;
    PetscRandom randx,randy,randz;
    PetscInt localnumParticles = 0;
    PetscReal Lx,Ly,Lz;
    DM dmcell;

    dmcell = user->dmcell;

    PetscOptionsInsertFile(PETSC_COMM_WORLD, PETSC_NULL, "control.dat", PETSC_TRUE);

    //  PetscOptionsGetReal(PETSC_NULL,PETSC_NULL, "-L_x", &Lx, PETSC_NULL);
      PetscOptionsGetReal(PETSC_NULL,PETSC_NULL, "-L_y", &Ly, PETSC_NULL);
      PetscOptionsGetReal(PETSC_NULL,PETSC_NULL, "-L_z", &Lz, PETSC_NULL);

    //    PetscPrintf(PETSC_COMM_WORLD," Lx - %f; Ly - %f; Lz - %f \n",Lx,Ly,Lz);

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD," rank - %d; size - %d \n",rank,size);

    localnumParticles = numParticles/size; // ***** MAKE SURE NO.OF PARTICLES IS A MULTIPLE OF NO.OF PROCESSORS ******************  VIK - 10/2024  
  
  /* Create a DMShell for point location purposes */
    //   ierr = DMShellCreate(PETSC_COMM_WORLD,dmcell);CHKERRQ(ierr);
    //  ierr = DMSetApplicationContext(dmcell,(void*)user->da);CHKERRQ(ierr);
    //  dmcell->ops->locatepoints = DMLocatePoints_DMDARegular;
    // dmcell->ops->getneighbors = DMGetNeighbors_DMDARegular;

    // Create the DMSwarm
    ierr = DMCreate(PETSC_COMM_WORLD, &user->swarm); CHKERRQ(ierr);
    ierr = DMSetType(user->swarm, DMSWARM); CHKERRQ(ierr);
    ierr = DMSetDimension(user->swarm,3);
    ierr = DMSwarmSetType(user->swarm, DMSWARM_BASIC); CHKERRQ(ierr);
    //  ierr = DMSwarmSetCellDM(user->swarm, user->da); CHKERRQ(ierr);  // dmcell

    // Register particle fields
      ierr = DMSwarmRegisterPetscDatatypeField(user->swarm, "position", 3, PETSC_REAL); CHKERRQ(ierr);
      ierr = DMSwarmRegisterPetscDatatypeField(user->swarm, "velocity", 3, PETSC_REAL); CHKERRQ(ierr);
    // Add other fields as needed
    // ..
    // ..
    // Finalize DMSwarm
    ierr = DMSwarmFinalizeFieldRegister(user->swarm); CHKERRQ(ierr);

    //  ierr = SwarmViewGP(dms,"step0");CHKERRQ(ierr);

    // Only initialize particles on rank 0
    // if (rank == 0) {
    ierr = DMSwarmSetLocalSizes(user->swarm, localnumParticles,4); CHKERRQ(ierr);

    ierr = DMView(user->swarm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        // Access particle data
        PetscReal *positions, *velocities,*coor;
        PetscInt64 *IDs;

	  ierr = DMSwarmGetField(user->swarm, "position",NULL, NULL, (void**)&positions); CHKERRQ(ierr);
	  ierr = DMSwarmGetField(user->swarm,"DMSwarm_pid",NULL,NULL,(void**)&IDs);
	  ierr = DMSwarmGetField(user->swarm, "velocity",NULL,NULL, (void**)&velocities); CHKERRQ(ierr);

	// Initialize random number generators
        // x 
	ierr = PetscRandomCreate(PETSC_COMM_WORLD, &randx); CHKERRQ(ierr);
	ierr = PetscRandomSetType(randx, PETSCRAND48); CHKERRQ(ierr);
	ierr = PetscRandomSetInterval(randx, user->bbox.min_coords.x, user->bbox.max_coords.x); CHKERRQ(ierr);
	ierr = PetscRandomSeed(randx); CHKERRQ(ierr);
        // y
	ierr = PetscRandomCreate(PETSC_COMM_WORLD, &randy); CHKERRQ(ierr);
	ierr = PetscRandomSetType(randy, PETSCRAND48); CHKERRQ(ierr);
ierr = PetscRandomSetInterval(randy, user->bbox.min_coords.y, user->bbox.max_coords.y); CHKERRQ(ierr);
	//	ierr = PetscRandomSetInterval(randy, 0.0, Ly); CHKERRQ(ierr);
	ierr = PetscRandomSeed(randy); CHKERRQ(ierr);
        // z
	ierr = PetscRandomCreate(PETSC_COMM_WORLD, &randz); CHKERRQ(ierr);
	ierr = PetscRandomSetType(randz, PETSCRAND48); CHKERRQ(ierr);
	ierr = PetscRandomSetInterval(randz, user->bbox.min_coords.z, user->bbox.max_coords.z); CHKERRQ(ierr);
	//	ierr = PetscRandomSetInterval(randz, 0.0, Lz); CHKERRQ(ierr);
	ierr = PetscRandomSeed(randz); CHKERRQ(ierr);

        // Initialize particle positions and velocities
        for (PetscInt p = 0; p < localnumParticles; p++) {
            // Set particle positions (e.g., randomly within the domain)
	  // for(PetscInt d = 0; d < 3; d++){
            ierr = PetscRandomGetValue(randx, &positions[p * 3 + 0]); CHKERRQ(ierr);
	    ierr = PetscRandomGetValue(randy, &positions[p * 3 + 1]); CHKERRQ(ierr);
	    ierr = PetscRandomGetValue(randz, &positions[p * 3 + 2]); CHKERRQ(ierr);
           
            velocities[3*p + 0] = 0.0;
            velocities[3*p + 1] = 0.0;
            velocities[3*p + 2] = 0.0;

	    IDs[p] = rank*localnumParticles + p; 
	    // } 
            
        }

        // Restore fields
	  ierr = DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
	  ierr = DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void**)&velocities); CHKERRQ(ierr);
	  ierr = DMSwarmRestoreField(user->swarm,"DMSwarm_pid",NULL,NULL,(void**)&IDs);
     
        ierr = PetscRandomDestroy(&randx); CHKERRQ(ierr);
        ierr = PetscRandomDestroy(&randy); CHKERRQ(ierr);
        ierr = PetscRandomDestroy(&randz); CHKERRQ(ierr);
        
	//   } else {

        // Ensure other ranks have zero local particles
	// ierr = DMSwarmSetLocalSizes(user->swarm, 0, 0); CHKERRQ(ierr);
        //  }

    return 0;
}

PetscErrorCode ParticleSwarmViewCoors(UserCtx* user) {
    DM swarm = user->swarm;
    PetscErrorCode ierr;
    PetscInt localnumParticles;
    PetscReal *coors;
    PetscMPIInt rank;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
 
    PetscPrintf(PETSC_COMM_SELF,"ParticleSwarmViewCoors - rank %d \n",rank);
    
    // Get the number of particles in the swarm
    ierr = DMSwarmGetLocalSize(swarm, &localnumParticles); CHKERRQ(ierr);

    // Access the particle positions
    ierr = DMSwarmGetField(swarm, "DMSwarmPIC_coor", NULL, NULL, (void**)&coors); CHKERRQ(ierr);

    // Print out the positions
    for (PetscInt i = 0; i < localnumParticles; i++) {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "rank - %d - Global Particle %d - Local Particle %d : Coors  = (%f, %f, %f)\n", rank,rank*localnumParticles + i+1,i+1,
                    coors[3 * i],coors[3 * i + 1], coors[3 * i + 2]);
    }
    
    PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
    // Restore the field
    ierr = DMSwarmRestoreField(swarm, "DMSwarmPIC_coor", NULL, NULL, (void**)&coors); CHKERRQ(ierr);

    return 0;
}

PetscErrorCode ParticleSwarmViewPositions(UserCtx* user) {
    DM swarm = user->swarm;
    PetscErrorCode ierr;
    PetscInt localnumParticles;
    PetscReal *positions;
    PetscInt64* IDs;
    PetscMPIInt *ranks;
    PetscMPIInt rank;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Get the number of particles in the swarm
    ierr = DMSwarmGetLocalSize(swarm, &localnumParticles); CHKERRQ(ierr);

    // Access the particle positions
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm,"DMSwarm_pid",NULL,NULL,(void**)&IDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm,"DMSwarm_rank",NULL,NULL,(void**)&ranks); CHKERRQ(ierr);

    // Print out the positions
    for (PetscInt i = 0; i < localnumParticles; i++) {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "rank - %d - DMSwarm_rank - %d - Global Particle %" PRId64 " - Local Particle %d : Position = (%f, %f, %f)\n",rank, ranks[i], IDs[i], i + 1, positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**************** \n");
    
    PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
    // Restore the field
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&positions); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm,"DMSwarm_pid",NULL,NULL,(void**)&IDs); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm,"DMSwarm_rank",NULL,NULL,(void**)&ranks); CHKERRQ(ierr);
    
    return 0;
}

PetscErrorCode ParticleSwarmBasicMigrationPattern(UserCtx* user) {
    DM swarm = user->swarm;
    PetscErrorCode ierr;
    PetscMPIInt * miglist;   
    PetscInt localnumParticles;
    PetscMPIInt rank,size;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
   
 
    DMSwarmGetLocalSize(swarm,&localnumParticles);
   
    PetscCalloc1(localnumParticles, &miglist);
  
    for(PetscInt p = 0; p < localnumParticles; p ++) miglist[p] = rank;

    // set conditions, can be modified for different cases //
    
    if(size>1){

    if(rank == 0){
      miglist[0] = 2;
      }

    //   if(rank == 1){

    //     }

    //   if(rank == 2){
      //   miglist[0] = 0;

    // ... add custom conditions .. //
    // ... //
    // ... // 
    }
 
    (user->miglist) = miglist;

    //  for(PetscInt p = 0; p < localnumParticles; p++){
    //     PetscSynchronizedPrintf(PETSC_COMM_WORLD,"ParticleSwarmBasicMigrationPattern -  rank %d - miglist[%d] - %d; user->miglist[p] - %d \n",rank,p,miglist[p],user->miglist[p]);
    //   }
   
    //   PetscSynchronizedPrintf(PETSC_COMM_WORLD," *********************** \n");
    //    PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

    return 0;

}

PetscErrorCode ParticleSwarmMigrateBasic(UserCtx* user) {
 
  DM swarm = user->swarm;
  
  PetscErrorCode ierr;
  
  PetscMPIInt* miglist;
  PetscMPIInt* rankval;
  
  PetscInt localnumParticles;
  PetscMPIInt rank;

  PetscBool removePoints = PETSC_TRUE;

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  
  ierr = ParticleSwarmBasicMigrationPattern(user); CHKERRQ(ierr);

  ierr = DMSwarmGetLocalSize(user->swarm, &localnumParticles);
   
  miglist = user->miglist;

  // for(PetscInt p = 0; p < localnumParticles; p++){
    //    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"ParticleSwarmMigrateBasic -  rank %d - miglist[%d] - %d; user->miglist[p] - %d \n",rank,p,miglist[p],user->miglist[p]);
    // }
    
  // PetscSynchronizedPrintf(PETSC_COMM_WORLD," *********************** \n");
  // PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
 
  ierr = DMSwarmGetField(user->swarm,"DMSwarm_rank",NULL,NULL,(void**)&rankval);CHKERRQ(ierr);


  for(int p = 0; p < localnumParticles; p++){
    rankval[p] = miglist[p];
  }


  //  for(PetscInt p = 0; p < localnumParticles; p++){
  //    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"ParticleSwarmMigrateBasic - After change -  rank %d - rankval[%d] - %d; user->miglist[p] - %d \n",rank,p,rankval[p],user->miglist[p]);
  // }
    
  // PetscSynchronizedPrintf(PETSC_COMM_WORLD," *********************** \n");
  // PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);


  ierr = DMSwarmRestoreField(user->swarm,"DMSwarm_rank",NULL,NULL,(void**)&rankval);CHKERRQ(ierr);
  
  // if (visflg == ctr) PetscPrintf(PETSC_COMM_WORLD, "main - Particle Swarm Initialized  - cp12.5 \n");
 
  // ierr = DMSwarmSetMigrateType(user->swarm,DMSWARM_MIGRATE_BASIC);

  ierr = DMSwarmMigrate(user->swarm,removePoints); CHKERRQ(ierr);

  ierr = PetscFree(user->miglist);
  
  PetscBarrier(NULL);


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
  
DM da,fda,pda;

UserCtx *user;

PetscInt ctr = 0; // Output counter is at 0 if only main outputs are necessary.

 PetscInt rank, size,bi, ibi,i;

PetscReal umax;

// BoundingBox *bboxlist;
  
ierr = PetscInitialize(&argc, &argv, (char *)0, help); if(ierr) return ierr;

PetscOptionsInsertFile(PETSC_COMM_WORLD, PETSC_NULL, "control.dat", PETSC_TRUE);

PetscOptionsGetInt(PETSC_NULL,PETSC_NULL, "-visflg", &visflg, PETSC_NULL);
  
if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "main - cp1 - control.dat read \n");

// Allocate memory for user

 PetscMalloc1(block_number, &user); // PetscMalloc1 automatically uses the sizeof(UserCtx) to allocate memory.

if(visflg==ctr) PetscPrintf(PETSC_COMM_SELF, "main - cp2 - user - %p  \n",user);

PetscBarrier(PETSC_NULL);

// Identify all processors 

MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

MPI_Comm_size(PETSC_COMM_WORLD, &size);

if(visflg==ctr) PetscPrintf(PETSC_COMM_SELF, "main - cp3 - rank - %d  \n", rank);

PetscBarrier(PETSC_NULL);

// Read no.of particles from control file 

//PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);

PetscOptionsGetInt(PETSC_NULL,PETSC_NULL, "-ti", &ti, PETSC_NULL);

if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "main - cp4 - ti - %d \n", ti);

PetscOptionsGetInt(PETSC_NULL,PETSC_NULL, "-numParticles", &np, PETSC_NULL);

if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "main - cp5 - No.of Particles - %d  \n",np);

// Read coordinates and create fda

ReadCoordinates(user);

if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "main - cp6 - Coordinates Read/Generated \n");

// Create Vector for ucat

PetscBarrier(PETSC_NULL);

 DMDALocalInfo info;
 ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
 PetscSynchronizedPrintf(PETSC_COMM_WORLD, "user %p - Rank %d owns the subdomain: xs=%d, xe=%d, ys=%d, ye=%d, zs=%d, ze=%d\n",user,
	     rank, info.xs, info.xs + info.xm, info.ys, info.ys + info.ym, info.zs, info.zs + info.zm);

 PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

for(bi = 0; bi < block_number; bi ++ ){

  ierr = DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat); CHKERRQ(ierr);

  if(visflg==ctr) PetscPrintf(PETSC_COMM_SELF, "main - cp7.1%d - rank - %d; ucat - %p \n",bi+1,rank,user[bi].Ucat);

}

PetscBarrier(PETSC_NULL);

// Read Ucat from ufield

Ucat_Binary_Input(user);

 ierr = VecNorm(user->Ucat,NORM_INFINITY,&umax); CHKERRQ(ierr);

if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "main - cp8 - max(ucat) - %f \n",umax);

PetscBarrier(PETSC_NULL);

ParticlesLocate(user,np);

// Create Bounding Boxes for each rank and store them locally.

// WriteBoundingBox(user,&bboxlist);

// Create and Initialize Particle //
               
//  ierr = ParticleSwarmCreate(user, np); CHKERRQ(ierr);

// if (visflg == ctr) PetscPrintf(PETSC_COMM_WORLD, "main - Particle Swarm Created and Initialized - cp9 \n");

 //   ParticleSwarmViewPositions(user);

// Migrate Particles i.e move them between ranks  
 
//  ParticleSwarmMigrateBasic(user);

// if (visflg == ctr) PetscPrintf(PETSC_COMM_WORLD, "main - Particle Swarm Migrated - cp10 \n");

 //  ParticleSwarmViewPositions(user);
/*

if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "main - cp 8A - np - %d  \n",np);

if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "main - cp 8B - user - %p  \n",user);

ParticleVectorCreate(np, user);

if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "main - Particle Vector Created - cp9 \n");

// Initialize Particles

ParticleVectorInitialize(user,np);

if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD, "main - Particle Vector Initialized - cp10 \n");

// Locate Particles in grid 


if(visflg==ctr) PetscPrintf(PETSC_COMM_WORLD,"main - cp11 \n");

*/

// Finalize 

 

 for(bi = 0; bi < block_number; bi ++){
   
   DMDestroy(&(user[bi].fda));

   DMDestroy(&(user[bi].da)); 

   DMDestroy(&(user[bi].swarm));
 }

 // if(rank == 0) free(bboxlist);

PetscFree(user);

PetscFinalize();

return 0;
}
