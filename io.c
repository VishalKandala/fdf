#include "variables.h"
extern PetscInt ti;
extern PetscInt averaging;

PetscErrorCode Ucont_P_Binary_Input(UserCtx *user)
{
  
  PetscViewer	viewer;
  
  char filen[90];
  
  
  sprintf(filen, "results/vfield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  PetscInt N;

  VecGetSize(user->Ucont, &N);
  PetscPrintf(PETSC_COMM_WORLD, "PPP %d\n", N);
  VecLoad((user->Ucont),viewer);
  
  PetscViewerDestroy(&viewer);

  PetscBarrier(PETSC_NULL);

  PetscOptionsClearValue("-vecload_block_size");

  sprintf(filen, "results/pfield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad((user->P),viewer);
  PetscViewerDestroy(&viewer);

  sprintf(filen, "results/nvfield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad((user->Nvert_o),viewer);
  PetscViewerDestroy(&viewer);

  DMGlobalToLocalBegin(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o);
  DMGlobalToLocalEnd(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o);

  if(averaging) {	// Seokkoo Kang
    sprintf(filen, "results/su0_%06d_%1d.dat",ti, user->_this);
    FILE *fp=fopen(filen, "r");
    
    VecSet(user->Ucat_sum, 0);
    VecSet(user->Ucat_cross_sum, 0);
    VecSet(user->Ucat_square_sum, 0);
    VecSet(user->P_sum, 0);
    
    if(fp==NULL) {
      PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, setting the statistical quantities to zero and continues the computation ... ***\n\n", filen);
    }
    else {
      fclose(fp);
      PetscBarrier(PETSC_NULL);
      sprintf(filen, "su0_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->Ucat_sum,viewer);
      PetscViewerDestroy(&viewer);
      
      sprintf(filen, "su1_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad(user->Ucat_cross_sum,viewer);
      PetscViewerDestroy(&viewer);
      
      sprintf(filen, "su2_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad((user->Ucat_square_sum),viewer);
      PetscViewerDestroy(&viewer);
      
      sprintf(filen, "sp_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->P_sum,viewer);
      PetscViewerDestroy(&viewer);
      
      PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Read %s, continuing averaging ... ***\n\n", filen);
    }
  }
  
  if(les) {
    Vec Cs;
    
    VecDuplicate(user->P, &Cs);
    
    sprintf(filen, "cs_%06d_%1d.dat", ti, user->_this);
    FILE *fp=fopen(filen, "r");
    
    if(fp==NULL) {
      PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, setting Cs to 0 and contiues the computation ... ***\n\n", filen);
      VecSet(Cs, 0);
    }
    else {
      fclose(fp);
      
      PetscBarrier(PETSC_NULL);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( Cs,viewer);
      PetscViewerDestroy(&viewer);
    }
    
    DMGlobalToLocalBegin(user->da, Cs, INSERT_VALUES, user->lCs);
    DMGlobalToLocalEnd(user->da, Cs, INSERT_VALUES, user->lCs);
    
    VecDestroy(&Cs);
  }

  
  if(rans) {
    // K-Omega
    sprintf(filen, "kfield%06d_%1d.dat", ti, user->_this);
    FILE *fp=fopen(filen, "r");
    
    if(fp!=NULL) {
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad(user->K_Omega,viewer);
      PetscViewerDestroy(&viewer);
    }
    else {
      K_Omega_IC(user);
    }
    
    VecCopy(user->K_Omega, user->K_Omega_o);
    
    DMGlobalToLocalBegin(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
    
    DMGlobalToLocalBegin(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);    
  }

  return 0;
}

PetscErrorCode Ucont_P_Binary_Output(UserCtx *user, PetscInt bi)
{
  PetscViewer	viewer;
  char filen[80];

  int rank;
  
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  sprintf(filen, "results/vfield%5.5d_%1.1d.dat", ti, user->_this);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  
  VecView(user->Ucont, viewer);
  
  PetscViewerDestroy(&viewer);

  
  sprintf(filen, "results/ufield%5.5d_%1.1d.dat", ti, user->_this);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  
  VecView(user->Ucat, viewer);
  
  PetscViewerDestroy(&viewer);
  
  sprintf(filen, "results/pfield%5.5d_%1.1d.dat", ti, user->_this);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  
  VecView(user->P, viewer);
  PetscViewerDestroy(&viewer);

  sprintf(filen, "results/nvfield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  
  VecView(user->Nvert, viewer);
  PetscViewerDestroy(&viewer);
  
  if(averaging && ti!=0) {	// Seokkoo Kang
    sprintf(filen, "su0_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat_sum, viewer);
    PetscViewerDestroy(&viewer);
    sprintf(filen, "su0_%06d_%1d.dat.info", ti, user->_this);	if(!rank) unlink(filen);
    
    PetscBarrier(PETSC_NULL);
    
    sprintf(filen, "su1_%06d_%1d.dat",  ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat_cross_sum, viewer);
    PetscViewerDestroy(&viewer);  
    sprintf(filen, "su1_%06d_%1d.dat.info", ti, user->_this);	if(!rank) unlink(filen);
    
    PetscBarrier(PETSC_NULL);
    
    sprintf(filen, "su2_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat_square_sum, viewer);
    PetscViewerDestroy(&viewer);
    sprintf(filen, "su2_%06d_%1d.dat.info",ti, user->_this);	if(!rank) unlink(filen);
    
    PetscBarrier(PETSC_NULL);
    
    sprintf(filen, "sp_%06d_%1d.dat",ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->P_sum, viewer);
    PetscViewerDestroy(&viewer);
    sprintf(filen, "sp_%06d_%1d.dat.info",ti, user->_this);	if(!rank) unlink(filen);
    
    PetscBarrier(PETSC_NULL);
  }
  
  if(les) {
    Vec Cs;
    
    VecDuplicate(user->P, &Cs);
    DMLocalToGlobalBegin(user->da, user->lCs, INSERT_VALUES, Cs);
    DMLocalToGlobalEnd(user->da, user->lCs, INSERT_VALUES, Cs);
    
    sprintf(filen, "cs_%06d_%1d.dat",  ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(Cs, viewer);
    PetscViewerDestroy(&viewer);
    
    PetscBarrier(PETSC_NULL);
    VecDestroy(&Cs);
  }
  
  if(rans) {
    sprintf(filen, "kfield%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->K_Omega, viewer);
    PetscViewerDestroy(&viewer);
    
    PetscBarrier(PETSC_NULL);
  }
  
  return 0;
}
