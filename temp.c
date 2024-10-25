//--------------- RESTART setup: if Restarting -----------------------------------
  level = usermg.mglevels-1;
  user = usermg.mgctx[level].user;
  if (rstart_flg) {
    ti = tistart; tistart++;
    
    for (bi=0; bi<block_number; bi++) {
      Ucont_P_Binary_Input(&(user[bi])); // Read ucont from binary file stored at  tistart.
      
      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont,INSERT_VALUES, user[bi].lUcont);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont,INSERT_VALUES, user[bi].lUcont);

      DMGlobalToLocalBegin(user[bi].da, user[bi].P,INSERT_VALUES, user[bi].lP);
      DMGlobalToLocalEnd(user[bi].da, user[bi].P,INSERT_VALUES, user[bi].lP);
      
      Contra2Cart(&(user[bi]));
     

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
       
	    Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi],user[bi].dt,ibi);
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


	    PetscReal rx,ry,rz;
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

	    PetscReal rx,ry,rz;
	    for (i=0; i<ibm[ibi].n_v; i++) {
	      rx = ibm[ibi].x_bp[i]-fsi[ibi].x_c;
	      ry = ibm[ibi].y_bp[i]-fsi[ibi].y_c;
	      rz = ibm[ibi].z_bp[i]-fsi[ibi].z_c;

	      ibm[ibi].u[i].x =-( ry*fsi[ibi].S_ang_n[5]-fsi[ibi].S_ang_n[3]*rz );
	      ibm[ibi].u[i].y = ( rx*fsi[ibi].S_ang_n[5]-fsi[ibi].S_ang_n[1]*rz );
	      ibm[ibi].u[i].z =-( rx*fsi[ibi].S_ang_n[3]-fsi[ibi].S_ang_n[1]*ry );

	      ibm[ibi].uold[i].x =-( ry*fsi[ibi].S_ang_r[5]-fsi[ibi].S_ang_r[3]*rz );
	      ibm[ibi].uold[i].y = ( rx*fsi[ibi].S_ang_r[5]-fsi[ibi].S_ang_r[1]*rz );
	      ibm[ibi].uold[i].z =-( rx*fsi[ibi].S_ang_r[3]-fsi[ibi].S_ang_r[1]*ry );

              ibm->nf_x[i] =  -ibm->nf_x[i];
              ibm->nf_y[i] =  -ibm->nf_y[i];
              ibm->nf_z[i] =  -ibm->nf_z[i];

	      ibm[ibi].urm1[i].x =-( ry*fsi[ibi].S_ang_rm1[5]-fsi[ibi].S_ang_rm1[3]*rz );
	      ibm[ibi].urm1[i].y = ( rx*fsi[ibi].S_ang_rm1[5]-fsi[ibi].S_ang_rm1[1]*rz );
	      ibm[ibi].urm1[i].z =-( rx*fsi[ibi].S_ang_rm1[3]-fsi[ibi].S_ang_rm1[1]*ry );
      
	    } // for n_v
	  } // for rotateframe and rotatefsi

	}//ibi
      } // if rstart fsi

    }// bi

  } // if rstart

 
