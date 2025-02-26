#include "variables.h"
PetscErrorCode ReadParameters(char* ctrl_file){
  PetscPrintf(PETSC_COMM_WORLD,"Parameters being read from the file : %s \n",ctrl_file);
// ------- INPUT PARAMETERS ----------------------
  PetscOptionsInsertFile(PETSC_COMM_WORLD,ctrl_file, PETSC_TRUE);
  PetscOptionsGetInt(PETSC_NULL, "-vis_flg", &visflg, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-tio", &tiout, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-imm", &immersed, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-catcorr", &catcorr, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-inv", &invicid, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rstart", &tistart, &rstart_flg);
  PetscOptionsGetInt(PETSC_NULL,"-rstart_fem",&rstart_fem,PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-imp", &implicit, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-fdf", &fdf, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-imp_type", &implicit_type, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-imp_MAX_IT", &imp_MAX_IT, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-fsi", &movefsi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rfsi", &rotatefsi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL,"-move_ibm",&moveibm,PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-radi", &radi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-inlet", &inletprofile, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-str", &STRONG_COUPLING, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rs_fsi", &rstart_fsi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-r_grid", &rotate_grid, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-cop", &cop, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-fish", &fish, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-aneurysm", &aneurysm, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-aneu_dom_bn", &aneu_dom_bn, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rheology", &rheology, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-turbine", &turbine, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-fishcyl", &fishcyl, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-eel", &eel, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-fish_c", &fish_c, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-wing", &wing, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-sediment", &sediment, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-mhv", &MHV, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-hydro", &hydro, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-lv", &LV, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-lvad", &LVAD, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-reg", &regime, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-twoD", &TwoD, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-thin", &thin, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_z", &dgf_z, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_y", &dgf_y, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_x", &dgf_x, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_az", &dgf_az, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_ay", &dgf_ay, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_ax", &dgf_ax, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-body", &NumberOfBodies, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-mframe", &moveframe, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rframe", &rotateframe, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-blk", &blank, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-init1", &InitialGuessOne, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-CMx_c", &(CMx_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-CMy_c", &(CMy_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-CMz_c", &(CMz_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-imp_atol", &(imp_atol), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-imp_rtol", &(imp_rtol), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-imp_stol", &(imp_stol), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-poisson_tol", &poisson_tol, PETSC_NULL);		// Seokkoo Kang: tolerance of implicit matrix free solver. 1.e-4 is enough for most cases.
  PetscOptionsGetInt(PETSC_NULL, "-les", &les, PETSC_NULL);				// Seokkoo Kang: if 1 Smagorinsky with Cs=0.1, if 2 Dynamic model
  PetscOptionsGetInt(PETSC_NULL, "-rans", &rans, PETSC_NULL);			// Seokkoo Kang
  PetscOptionsGetInt(PETSC_NULL, "-wallfunction", &wallfunction, PETSC_NULL);	// Seokkoo Kang: 1 or 2
  PetscOptionsGetInt(PETSC_NULL, "-mixed", &mixed, PETSC_NULL);			// Seokkoo Kang: mixed model option for LES
  PetscOptionsGetInt(PETSC_NULL, "-clark", &clark, PETSC_NULL);			// Seokkoo Kang: mixed model option for LES
  PetscOptionsGetInt(PETSC_NULL, "-dynamic_freq", &dynamic_freq, PETSC_NULL);		// Seokkoo Kang: LES dynamic compute frequency 
  PetscOptionsGetInt(PETSC_NULL, "-averaging", &averaging, PETSC_NULL);	// Seokkoo Kang: if 1 do averaging; always begin with -rstart 0
  PetscOptionsGetInt(PETSC_NULL, "-grid1d", &grid1d, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-i_periodic", &i_periodic, PETSC_NULL);	
  PetscOptionsGetInt(PETSC_NULL, "-j_periodic", &j_periodic, PETSC_NULL);	
  PetscOptionsGetInt(PETSC_NULL, "-k_periodic", &k_periodic, PETSC_NULL);	
  PetscOptionsGetInt(PETSC_NULL, "-pbc_domain", &blkpbc, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-testfilter_ik", &testfilter_ik, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-testfilter_1d", &testfilter_1d, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-poisson", &poisson, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-i_homo_filter", &i_homo_filter, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-j_homo_filter", &j_homo_filter, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-k_homo_filter", &k_homo_filter, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-max_cs", &max_cs, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-St_exp", &(St_exp), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-wlngth", &(wavelength), PETSC_NULL);\
  PetscOptionsGetString(PETSC_NULL, "-orient", orient,sizeof(orient),PETSC_NULL);
  PetscOptionsGetString(PETSC_NULL, "-g_orient", gridrotorient,sizeof(orient),PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-totalsteps", &tisteps, &flg); // Get input for real timesteps
  if(visflg==5) PetscPrintf(PETSC_COMM_WORLD,"Parameters setup done! \n");  
  return 0;
}

PetscErrorCode CaseInitialize(){
 
  // Case specific parameters are adjusted, this function can be extended as more cases are added to the code. 
    L_dim = 1.;
  
  if(MHV) max_angle = -54.*PI/180.; 

  if(fdf){
    PetscPrintf(PETSC_COMM_WORLD,"fdf-test \n");
    numParticles = 1;
  }  

  if(LVAD){  // For LVAD; Set Characteristic Length to be 1.0, number of immersed boundaries to be 1 and also set the center of moment of the immersed boundary to be at 0,0,0.
    L_dim=1.; 
    
} 
 

  return 0;
}

PetscErrorCode IBMInitialize(IBMNodes *ibm, FSInfo *fsi){
  
  PetscInt ibi,bi,i;
  // The immersed boundary is initialized by allocating memory and calling the FSI Initialization also. 
 
  if (immersed) {  // If an immersed boundary exists, allocate memory for IBMNodes and FSInfo (data-types), then call the FsiInitialize() function.
    PetscMalloc(NumberOfBodies*sizeof(IBMNodes), &ibm);
    PetscMalloc(NumberOfBodies*sizeof(FSInfo), &fsi);
    for (ibi=0;ibi<NumberOfBodies;ibi++)
         FsiInitialize(0, &fsi[ibi], ibi);     
  }
  if(visflg==5) PetscPrintf(PETSC_COMM_WORLD,"memory allocated for ibmnodes \n");
   
  return 0;
}


PetscErrorCode CaseSetup(UserMG *usermg,IBMNodes *ibm,FSInfo *fsi){
 
  UserCtx *user;  // local user pointer that is used to point to data in usermg.
  PetscInt i,bi,ibi;
  PetscInt level;
  MGCtx *mgctx;

  mgctx = usermg->mgctx;
  
 
   // Case specific parameters are adjusted, this function can be extended as more cases are added to the code. 
   // if(fdf){
   
   //   for(bi=0; bi<block_number;bi++){

   //   PetscPrintf(PETSC_COMM_WORLD," in-main - numparticle %d \n",numParticles);
                
	    //	    ParticleVectorCreate(numParticles,user[bi]);
	   
   //       ParticleVectorInitialize(&user[bi],numParticles);
   //   }
   // }

  if (immersed) {
    level = usermg->mglevels-1;
    user = mgctx[level].user;
    for (bi=0; bi<block_number; bi++) {
      
      PetscMalloc(NumberOfBodies*sizeof(IBMList), &(user[bi].ibmlist)); //For each block, allocate memory for IBMList and also initialize it.
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	InitIBMList(&(user[bi].ibmlist[ibi]));
      }
    }

    ibi = 0;
    ibm_read_Icem(&ibm[ibi],ibi);  // This function reads the grid either from grid.dat or from the cartesian grid setup in the control file. 
    if(visflg==5) PetscPrintf(PETSC_COMM_WORLD,"IBM Read complete \n"); 
    if(rotatefsi)  Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi],user[bi].dt,ibi); // Rotate the immersed boundary if rotation is turned  on.
    ibm_surface_VTKOut(&ibm[ibi],ibi,0); // This function generates the VTK file that can be used to visualize the immersed boundary surface on ParaVi
    if(visflg==5)PetscPrintf(PETSC_COMM_WORLD,"IBM Surface Out complete \n");   
    PetscBarrier(PETSC_NULL);   

  } // immersed 

  return 0;
}


PetscErrorCode HandleRestart(UserMG *usermg, IBMNodes *ibm, FSInfo *fsi){

  if (tistart==0) tisteps ++; // if the simulation is starting from zero,run until the final timestep(hence increase timesteps by 1, so loop can be less than timesteps)

  // If restarting from a certain timestep, the data is setup in this function.

  UserCtx *user;
  PetscInt i,bi,ibi;
  PetscInt level;
  MGCtx *mgctx;
  PetscReal rx,ry,rz;
 
 //--------------- RESTART setup: if Restarting -----------------------------------
  level = usermg->mglevels-1;
  mgctx = usermg->mgctx;
  user = mgctx[level].user;
 
  
   
  if (rstart_flg) {
    ti = tistart; 
    tistart++;
    
    for (bi=0; bi<block_number; bi++) {
      Ucont_P_Binary_Input(&(user[bi])); // Read ucont from binary file stored at  tistart.
      
      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont,INSERT_VALUES, user[bi].lUcont);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont,INSERT_VALUES, user[bi].lUcont);

      DMGlobalToLocalBegin(user[bi].da, user[bi].P,INSERT_VALUES, user[bi].lP);
      DMGlobalToLocalEnd(user[bi].da, user[bi].P,INSERT_VALUES, user[bi].lP);
      
      Contra2Cart(&(user[bi]));
     
      // FSI Restart Options 
      if (rstart_fsi) {
	for (ibi=1;ibi<NumberOfBodies;ibi++) {

	  FSI_DATA_Input(&fsi[ibi],ibi);

	  if (movefsi) {
	    if (!moveframe)
	      Elmt_Move_FSI_TRANS(&fsi[ibi], &ibm[ibi],ibi);
	    for (i=0;i<6;i++){
	      fsi[ibi].S_realm1[i]=fsi[ibi].S_real[i];
	      fsi[ibi].S_real[i]=fsi[ibi].S_new[i];
	    }
	    for (i=0; i<ibm[ibi].n_v; i++) {
	      ibm[ibi].uold[i].x = fsi[ibi].S_real[1];
	      ibm[ibi].uold[i].y = fsi[ibi].S_real[3];
	      ibm[ibi].uold[i].z = fsi[ibi].S_real[5];
	    }
	    for (i=0; i<ibm[ibi].n_v; i++) {
	      ibm[ibi].urm1[i].x = fsi[ibi].S_realm1[1];
	      ibm[ibi].urm1[i].y = fsi[ibi].S_realm1[3];
	      ibm[ibi].urm1[i].z = fsi[ibi].S_realm1[5];
	    }
	  } // for movefsi

	  if ((rotatefsi &&  !rotateframe)|| MHV) {
	    //rotate_ibm(&ibm[ibi],&fsi[ibi]);
	    // calc_ibm_normal(&ibm[ibi]);
	    /* // */
	    Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi],user[bi].dt,ibi);
	    /* // */
	    ibm_surface_VTKOut(&ibm[ibi],ibi,0);
	    // if read ti, then will start for ti+1
	    for (i=0;i<6;i++){
	      fsi[ibi].S_ang_rm1[i]=fsi[ibi].S_ang_r[i];
	      fsi[ibi].S_ang_r[i]=fsi[ibi].S_ang_n[i];
	    }

	    fsi[ibi].F_x_real=fsi[ibi].F_x;
	    fsi[ibi].F_y_real=fsi[ibi].F_y;
	    fsi[ibi].F_z_real=fsi[ibi].F_z;
 
	    fsi[ibi].M_x_rm3=fsi[ibi].M_x;
	    fsi[ibi].M_y_rm3=fsi[ibi].M_y;
	    fsi[ibi].M_z_rm3=fsi[ibi].M_z;
	    
	    fsi[ibi].M_x_rm2=fsi[ibi].M_x;
	    fsi[ibi].M_y_rm2=fsi[ibi].M_y;
	    fsi[ibi].M_z_rm2=fsi[ibi].M_z;

	    fsi[ibi].M_x_real=fsi[ibi].M_x;
	    fsi[ibi].M_y_real=fsi[ibi].M_y;
	    fsi[ibi].M_z_real=fsi[ibi].M_z;



	    for (i=0; i<ibm[ibi].n_v; i++) {
	      rx = ibm[ibi].x_bp[i]-fsi[ibi].x_c;
	      ry = ibm[ibi].y_bp[i]-fsi[ibi].y_c;
	      rz = ibm[ibi].z_bp[i]-fsi[ibi].z_c;

	      ibm[ibi].u[i].x =   0.0  ;
	      ibm[ibi].u[i].y =-( -fsi[ibi].S_ang_n[1]*rz );
	      ibm[ibi].u[i].z =   -fsi[ibi].S_ang_n[1]*ry  ;

	      ibm[ibi].uold[i].x =  0.0  ;
	      ibm[ibi].uold[i].y =-( -fsi[ibi].S_ang_r[1]*rz );
	      ibm[ibi].uold[i].z =   -fsi[ibi].S_ang_r[1]*ry  ;
	      
	      ibm[ibi].urm1[i].x =   0.0  ;
	      ibm[ibi].urm1[i].y =-( -fsi[ibi].S_ang_rm1[1]*rz );
	      ibm[ibi].urm1[i].z =   -fsi[ibi].S_ang_rm1[1]*ry  ;
      
	    } // for nv
	  } // for rotatefsi and not rotateframe

	  if (rotatefsi &&  rotateframe) {
	    // Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi], user[bi].dt);
	    fsi[ibi].S_ang_n[5]=St_exp;
	    fsi[ibi].S_ang_r[5]=St_exp;

	    // if read ti, then will start for ti+1
	    for (i=0;i<6;i++){
	      fsi[ibi].S_ang_rm1[i]=fsi[ibi].S_ang_r[i];
	      fsi[ibi].S_ang_r[i]=fsi[ibi].S_ang_n[i];
	    }

	    fsi[ibi].F_x_real=fsi[ibi].F_x;
	    fsi[ibi].F_y_real=fsi[ibi].F_y;
	    fsi[ibi].F_z_real=fsi[ibi].F_z;
 
	    fsi[ibi].M_x_rm3=fsi[ibi].M_x;
	    fsi[ibi].M_y_rm3=fsi[ibi].M_y;
	    fsi[ibi].M_z_rm3=fsi[ibi].M_z;
	    
	    fsi[ibi].M_x_rm2=fsi[ibi].M_x;
	    fsi[ibi].M_y_rm2=fsi[ibi].M_y;
	    fsi[ibi].M_z_rm2=fsi[ibi].M_z;

	    fsi[ibi].M_x_real=fsi[ibi].M_x;
	    fsi[ibi].M_y_real=fsi[ibi].M_y;
	    fsi[ibi].M_z_real=fsi[ibi].M_z;

	    for (i=0; i<ibm[ibi].n_v; i++) {
	      rx = ibm[ibi].x_bp[i]-fsi[ibi].x_c;
	      ry = ibm[ibi].y_bp[i]-fsi[ibi].y_c;
	      rz = ibm[ibi].z_bp[i]-fsi[ibi].z_c;

	      ibm[ibi].u[i].x =-( ry*fsi[ibi].S_ang_n[5]-fsi[ibi].S_ang_n[3]*rz );
	      ibm[ibi].u[i].y = ( rx*fsi[ibi].S_ang_n[5]-fsi[ibi].S_ang_n[1]*rz );
	      ibm[ibi].u[i].z =-( rx*fsi[ibi].S_ang_n[3]-fsi[ibi].S_ang_n[1]*ry );

	      ibm[ibi].uold[i].x =-( ry*fsi[ibi].S_ang_r[5]-fsi[ibi].S_ang_r[3]*rz );
	      ibm[ibi].uold[i].y = ( rx*fsi[ibi/*  ibm->nf_x[i] =  -ibm->nf_x[i]; */
/*     ibm->nf_y[i] =  -ibm->nf_y[i]; */
/*     ibm->nf_z[i] =  -ibm->nf_z[i]; */].S_ang_r[5]-fsi[ibi].S_ang_r[1]*rz );
	      ibm[ibi].uold[i].z =-( rx*fsi[ibi].S_ang_r[3]-fsi[ibi].S_ang_r[1]*ry );

	      ibm[ibi].urm1[i].x =-( ry*fsi[ibi].S_ang_rm1[5]-fsi[ibi].S_ang_rm1[3]*rz );
	      ibm[ibi].urm1[i].y = ( rx*fsi[ibi].S_ang_rm1[5]-fsi[ibi].S_ang_rm1[1]*rz );
	      ibm[ibi].urm1[i].z =-( rx*fsi[ibi].S_ang_rm1[3]-fsi[ibi].S_ang_rm1[1]*ry );
      
	    } // for n_v
	  } // for rotateframe and rotatefsi

	}//ibi
      } // if rstart fsi

    }// bi

  } // if rstart

  return 0;

}

PetscErrorCode SearchInitialize(UserMG *usermg, IBMNodes *ibm){

  // do the search once if elmt is not moving!
  PetscInt bis=0,bie=block_number;
  PetscInt bi,i,ibi;
  PetscInt level;
  UserCtx *user;
  MGCtx *mgctx;
  
  
  if (immersed) {
    for (level = usermg->mglevels-1; level>=usermg->mglevels-1; level--) {
      
      mgctx = usermg->mgctx;
      user = mgctx[level].user;

      for (bi=bis; bi<bie; bi++) {  // chimera fluid grid blocks.
 
	user[bi].ibm=ibm;

	for (ibi=0;ibi<NumberOfBodies;ibi++) {  // immersed bodies
	    if(LVAD){
	      ibm_search_advanced_rev(&(user[bi]), &ibm[ibi], ibi);  
  	      if(visflg)  PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA REV LVAD  ibi %d bi %d\n", ibi,bi);        
            } 
            else {
	      ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
             if(visflg)   PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA General ibi %d bi %d\n", ibi, bi); 
	    }
	} //ibi

	PetscBarrier(PETSC_NULL);
	for (ibi=0;ibi<NumberOfBodies;ibi++) {
	  if(visflg)  PetscPrintf(PETSC_COMM_WORLD, "IBM_INTP\n");
	  ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi, 1);
	  if(visflg)  PetscPrintf(PETSC_COMM_WORLD, "IBM_INTP End\n");

	} //ibi
      }// bi
    } //mglevels
  } // immersed
  return 0;
}


PetscErrorCode BC_ICSetup(UserMG *usermg){
      
      PetscInt i,bi;
      PetscInt level;
      UserCtx *user;
      MGCtx *mgctx;
      PetscReal normZ, normX, normY;
      
      level = usermg->mglevels-1;
      mgctx = usermg->mgctx;
      user = mgctx[level].user;
         
      for (bi=0; bi<block_number; bi++) {
      // Initial Condition (if t = 0; if not, it is set by that last file used to restart.)
        if(ti = 0){
	  if(InitialGuessOne){
            SetInitialGuessToOne(&(user[bi]));   // Setup the initial condition for flow
	  } // IC 
        } // ti = 0

     // Boundary Conditions
        InflowFlux(&(user[bi]));  // Setup Inflow boundary condition
        OutflowFlux(&(user[bi])); // Setup Outflow boundary condition
        FormBCS(&(user[bi]));     // Setup all boundary conditions with info from bcs.dat
        if (wallfunction) {       // If wall function is used, setup Inflow and BCS again.
	 InflowFlux(&(user[bi]));
	 FormBCS(&(user[bi]));
        }
 
      } // chimera blocks 
      
      if (block_number>1) {
	Block_Interface_U(user);   // If more than one blocks(meshes) are present, setup interface conditions
      } 
     

    /* Print initial max Cartesian velocities */
    // This part of the code is useful for verifying the initial setup of velocity fields.
 
      if(visflg){ 
	for (bi=0; bi<block_number; bi++) {
	  VecStrideMax(user[bi].Ucat, 0, PETSC_NULL, &normX);
	  VecStrideMax(user[bi].Ucat, 1, PETSC_NULL, &normY);
	  VecStrideMax(user[bi].Ucat, 2, PETSC_NULL, &normZ);
	  PetscPrintf(PETSC_COMM_WORLD, "Initial max cartesian velocities along: x -  %le; y- %le; z- %le\n",normX, normY, normZ);
	} // for bi
      } // visflg

  return 0;
}

PetscErrorCode SaveHistoryFields(UserMG *usermg){
  
        
  PetscInt i,bi;
  PetscInt level;
  UserCtx *user;
  MGCtx *mgctx;      
  level = usermg->mglevels-1;
  mgctx = usermg->mgctx;
  user = mgctx[level].user;
  
  for (bi=0; bi<block_number; bi++) {

    if (immersed) {
    	VecCopy(user[bi].Nvert, user[bi].Nvert_o);   // Saving nvert at ti; storing in a local array.
        DMGlobalToLocalBegin(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
    	DMGlobalToLocalEnd(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
    }
 
    VecCopy(user[bi].Ucont, user[bi].Ucont_o);      // Copy Ucont to Ucont_o for the finest level
    VecCopy(user[bi].Ucont, user[bi].Ucont_rm1);    // Copy Ucont to Ucont_rm1 for the finest level
    VecCopy(user[bi].Ucat, user[bi].Ucat_o);        // Copy Ucat to Ucat_o for the finest level
    VecCopy(user[bi].P, user[bi].P_o);              // Copy P to P_o for the finest level

    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); // Global to Local for ucont,ucont_o and ucont_rm.
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);

    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
  } // for bi

  

  return 0;
}


PetscErrorCode CoupledSolve(UserMG *usermg, IBMNodes *ibm, FSInfo *fsi, FE *fem){

    PetscBool DoSCLoop= PETSC_TRUE ; // if TRUE- Strong Coupling, if FALSE-weak coupling.
    PetscInt itr_sc = 0; // Strong coupling iteration count.
    Cstart     cstart; // used for fish case.

    while (DoSCLoop) {
      itr_sc++;
      if(visflg == 5) PetscPrintf(PETSC_COMM_WORLD, "SC LOOP itr # %d\n", itr_sc);
      if (immersed){
      	 Struc_Solver(usermg, ibm, fsi, &cstart, itr_sc, tistart, &DoSCLoop, fem);        //Structural Solver!       
      } else  DoSCLoop = PETSC_FALSE;
      Flow_Solver(usermg, ibm, fsi);     //Flow Solver!
    }// End of while SC loop
  return 0;
}

PetscErrorCode CalculateMaxError(UserMG *usermg,PetscReal * MaxError){
   
  /*------Calculate Max Error------- */

  PetscInt bi;
  PetscInt level;
  UserCtx *user;
  MGCtx *mgctx;
  Vec ErrorVec;
  PetscReal Error;
      
  level = usermg->mglevels-1;
  mgctx = usermg->mgctx;
  user = mgctx[level].user;
  
  *MaxError=0.0;
    
  for (bi=0; bi<block_number; bi++) {
    Error = 0.0;   
    VecDuplicate(user[bi].Ucat,&ErrorVec); 
    VecWAXPY(ErrorVec,-1.0,user[bi].Ucat,user[bi].Ucat_o);
    VecNorm(ErrorVec,NORM_INFINITY,&Error);
    *MaxError = PetscMax(*MaxError,Error);
    VecDestroy(&ErrorVec);
  } // bi
  
  return 0;
}

PetscErrorCode SaveIBM_FSIHistory(IBMNodes *ibm, FSInfo* fsi){
  
  PetscInt ibi,i;
  
  for (ibi=0;ibi<NumberOfBodies;ibi++) {

    for (i=0; i<ibm[ibi].n_v; i++) {   // loop through every immersed boundary node
      ibm[ibi].x_bp_o[i] = ibm[ibi].x_bp[i]; // copy x,y,z co-ordinates of each node to old co-ordinate arrays.
      ibm[ibi].y_bp_o[i] = ibm[ibi].y_bp[i];
      ibm[ibi].z_bp_o[i] = ibm[ibi].z_bp[i];

      ibm[ibi].urm1[i].x = ibm[ibi].uold[i].x;   // copy x,y,z velocities(? what is ibm.u??) to old arrays.
      ibm[ibi].urm1[i].y = ibm[ibi].uold[i].y;
      ibm[ibi].urm1[i].z = ibm[ibi].uold[i].z;

      ibm[ibi].uold[i].x = ibm[ibi].u[i].x;
      ibm[ibi].uold[i].y = ibm[ibi].u[i].y;
      ibm[ibi].uold[i].z = ibm[ibi].u[i].z;
    }
    // FSI Related data to be stored for next iteration (elaborated later)
    for (i=0;i<6;i++){
      fsi[ibi].S_realm1[i]=fsi[ibi].S_real[i];
      fsi[ibi].S_real[i]=fsi[ibi].S_new[i];

      fsi[ibi].S_ang_rm1[i]=fsi[ibi].S_ang_r[i];
      fsi[ibi].S_ang_r[i]=fsi[ibi].S_ang_n[i];
    }
      

    fsi[ibi].F_x_real=fsi[ibi].F_x;
    fsi[ibi].F_y_real=fsi[ibi].F_y;
    fsi[ibi].F_z_real=fsi[ibi].F_z;
      
    fsi[ibi].M_x_rm3=fsi[ibi].M_x_rm2;
    fsi[ibi].M_y_rm3=fsi[ibi].M_y_rm2;
    fsi[ibi].M_z_rm3=fsi[ibi].M_z_rm2;

    fsi[ibi].M_x_rm2=fsi[ibi].M_x_real;
    fsi[ibi].M_y_rm2=fsi[ibi].M_y_real;
    fsi[ibi].M_z_rm2=fsi[ibi].M_z_real;

    fsi[ibi].M_x_real=fsi[ibi].M_x;
    fsi[ibi].M_y_real=fsi[ibi].M_y;
    fsi[ibi].M_z_real=fsi[ibi].M_z;

  } //ibi  

  return 0;
}
