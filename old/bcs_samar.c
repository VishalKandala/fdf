#include "variables.h"
#include "petscvec.h" 

extern PetscInt ti, moveframe, blank;
extern PetscReal FluxInSum, FluxOutSum;
extern PetscReal Flux_in, angle,CMy_c, CMx_c,CMz_c;
extern PetscInt block_number;
extern PetscInt inletprofile;
extern PetscReal flux_flow[8700];
PetscReal U_bc;

PetscErrorCode Contra2Cart(UserCtx *user);
PetscErrorCode VTKOut(UserCtx *user);
void wall_function (UserCtx *user, double sc, double sb, 
		    Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, 
		    PetscReal *ustar, double nx, double ny, double nz);
void Calculate_normal(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], 
		      double nj[3], double nk[3]);
PetscErrorCode GhostNodeVelocity(UserCtx *user);


#define INLET 5
#define OUTLET 4
#define SOLIDWALL 1
#define SYMMETRIC 3
#define FARFIELD 6

PetscErrorCode InletRead(UserCtx *user)
{
  PetscInt	i, j, rank;

  PetscPrintf(PETSC_COMM_WORLD, "Inlet Read\n");

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (!rank) {
    FILE *fd;
    fd = fopen("inlet.dat", "r");
/*   for (i=0; i<100; i++) { */
/*     fscanf(fd, "%le", &(user->r[i])); */
/*   } */
/*   for (i=0; i<100; i++) { */
/*     fscanf(fd, "%le", &(user->tin[i])); */
/*   } */

    for (i=0; i<101; i++) {
      user->r[i] = i * 0.005;
    }
    for (j=0; j<1001; j++) {
      for (i=0; i<101; i++) {
	fscanf(fd, "%le", &(user->uinr[i][j]));
      }
    }
    
    fclose(fd);
    PetscPrintf(PETSC_COMM_WORLD, "Inlet Read End\n");
    MPI_Bcast(user->r, 101, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(user->uinr, 1001*101, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  else {
    MPI_Bcast(user->r, 101, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(user->uinr, 1001*101, MPIU_REAL, 0, PETSC_COMM_WORLD);
  }
  return 0;
}

PetscReal InletInterpolation(PetscReal r, UserCtx *user)
{
  PetscInt i;
  PetscReal temp;
  PetscInt tstep, ts_p_cycle=2000;

  tstep = ti/2 - ((ti / ts_p_cycle) * ts_p_cycle);

  if (inletprofile == 8) {
    Flux_in=sin(2*3.14159*ti/200.);
    return(Flux_in);
  }
  
  if (inletprofile==3 || inletprofile==6) 
    return(Flux_in);
  
  if (r>1.) {
    temp = user->uinr[99][tstep];
    return(temp);
  }
  for (i=0; i<100; i++) {
    if (r>= (user->r[i]) && r< (user->r[i+1])) {
      temp = user->uinr[i][tstep] + (user->uinr[i+1][tstep] - user->uinr[i][tstep]) *
	(r-user->r[i]) / (user->r[i+1]-user->r[i]);
    
      return (temp);
    }
  }
  return 0;   
}

PetscErrorCode InflowFlux(UserCtx *user) 
{
  //  PetscInt bi;
  PetscInt   i, j, k;
  PetscReal  FluxIn, r, uin, xc, yc,zc;
  Vec        Coor;
  Cmpnts     ***ucont, ***ubcs, ***ucat, ***coor, ***csi, ***eta, ***zet;
  PetscReal  H=4.1, Umax=1.5;
  PetscReal  lAreaIn, lAreaOut, AreaSumIn, AreaSumOut;
  PetscReal	mat[3][3], det, det0, det1, det2;
  PetscReal	q[3]; //local working array

  DM            da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscReal	***nvert; //local working array

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  // DMDAGetGhostedCoordinates(da, &Coor);

  DMGetCoordinatesLocal(da, &Coor);

  DMDAVecGetArray(fda, Coor, &coor);
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecGetArray(fda, user->Ucat,  &ucat);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);

  if (user->bctype[4] == FARFIELD) FluxIn=0.;
  if (user->bctype[4] == SOLIDWALL) FluxIn=0.;
  if (user->bctype[4] == SYMMETRIC) FluxIn=0.;

  if (user->bctype[4] == INLET) {
  
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  
	  xc = (coor[k][j][i  ].x + coor[k][j-1][i  ].x +
		coor[k][j][i-1].x + coor[k][j-1][i-1].x) * 0.25-CMx_c;
	  yc = (coor[k][j][i  ].y + coor[k][j-1][i  ].y +
		coor[k][j][i-1].y + coor[k][j-1][i-1].y) * 0.25-CMy_c;
	  zc = (coor[k][j][i  ].z + coor[k][j-1][i  ].z +
		coor[k][j][i-1].z + coor[k][j-1][i-1].z) * 0.25-CMz_c;
	  
	  r = sqrt(zc * zc + yc * yc);
	  
	  if (inletprofile == 0) {
	    uin = InletInterpolation(r, user);
	  } else if (inletprofile == 1) {
	    uin=1.;
	  } else if (inletprofile == -1) {
	    uin = -1.;
	  } else if (inletprofile == 2) {
	    uin = 4.*Umax*yc*(H-yc)/(H*H);//InletInterpolation(r, user);
	  } else if (inletprofile == 3 || inletprofile == 6 || inletprofile == 8) {
	    uin = InletInterpolation(r, user);
	  } else if (inletprofile == 4) { //fully-developed pipe flow
	    uin = 2.*(1.-4.*r*r);
	    if (r>0.5) uin=0.; 
	   
	  } else if (inletprofile == 5) {
	    uin = -InletInterpolation(r, user);
	  } else if (inletprofile == 7) {
	    uin = 1.5*(1.-4.*yc*yc);
	  } else if (inletprofile == 10) {
	    double _y = (yc-1.5)*2., _x=xc-0.5;
			double A=1.;
			#ifndef M_PI 
			#define M_PI 3.14159265358979323846264338327950288
			#endif
			uin = 1- _y*_y;
			int n;
			for(n=0; n<20; n++) {	// exact solution for fully developed flows
				uin -= 4.*pow(2./M_PI,3) *  
				  pow(-1., n) / pow(2*n+1.,3.) * 
				  cosh((2*n+1)*M_PI*_x/2.) * 
				  cos((2.*n+1)*M_PI*_y/2.) / 
				  cosh((2*n+1)*M_PI*A/2.);
			}
	  } else if (inletprofile == 11) {
	    uin = 0.185;//0.2654;//	  			 
	  } else {
	    PetscPrintf(PETSC_COMM_SELF, "WRONG INLET PROFILE TYPE!!!! U_in = 0\n");
	    uin = 0.;
	  }

	  if (nvert[k+1][j][i]<0.1) {
	
	  ubcs[k][j][i].y = uin*zet[k][j][i].y/sqrt(zet[k][j][i].z*zet[k][j][i].z +
						    zet[k][j][i].y*zet[k][j][i].y +
						    zet[k][j][i].x*zet[k][j][i].x);
	  ubcs[k][j][i].z = uin*zet[k][j][i].z/sqrt(zet[k][j][i].z*zet[k][j][i].z +
						    zet[k][j][i].y*zet[k][j][i].y +
						    zet[k][j][i].x*zet[k][j][i].x);
	  ubcs[k][j][i].x = uin*zet[k][j][i].x/sqrt(zet[k][j][i].z*zet[k][j][i].z +
						    zet[k][j][i].y*zet[k][j][i].y +
						    zet[k][j][i].x*zet[k][j][i].x);
	  
	  ucont[k][j][i].z = uin*sqrt(zet[k][j][i].z*zet[k][j][i].z +
				      zet[k][j][i].y*zet[k][j][i].y +
				      zet[k][j][i].x*zet[k][j][i].x);
	  }
	}
      }
    }
    
    DMDAVecRestoreArray(fda, user->Ucont, &ucont);
    
    PetscBarrier(PETSC_NULL);
    
    
    DMDAVecGetArray(fda, user->Ucont, &ucont);
    FluxIn = 0;
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k+1][j][i]<0.1)
	  FluxIn += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxIn = 0;
    }

    MPI_Allreduce(&FluxIn, &FluxInSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      

  } else if (user->bctype[5] == INLET) {

    /*     calc lAreaOut */
    if (zs == 0) {
      lAreaOut=0.;
      k = zs;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k+1][j][i]<0.1)
	  lAreaOut += sqrt(zet[k][j][i].z*zet[k][j][i].z+
			   zet[k][j][i].y*zet[k][j][i].y+
			   zet[k][j][i].x*zet[k][j][i].x);

	}
      }
    }
    else {
      lAreaOut = 0.;
    }
    MPI_Allreduce(&lAreaOut,&AreaSumOut,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
   

    /*     calc lAreaIn */
    if (ze == mz) {
      lAreaIn=0.;
      k = ze-2;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i]<0.1)
	  lAreaIn += sqrt(zet[k][j][i].z*zet[k][j][i].z+
			  zet[k][j][i].y*zet[k][j][i].y+
			  zet[k][j][i].x*zet[k][j][i].x);

	}
      }
    }
    else {
      lAreaIn = 0.;
    }
    MPI_Allreduce(&lAreaIn,&AreaSumIn,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
   

    if (ze==mz) {
      k = mz-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  
	  xc = (coor[k-1][j][i].x + coor[k-1][j-1][i].x +
		coor[k-1][j][i-1].x + coor[k-1][j-1][i-1].x) * 0.25;
	  yc = (coor[k-1][j][i].y + coor[k-1][j-1][i].y +
		coor[k-1][j][i-1].y + coor[k-1][j-1][i-1].y) * 0.25;
	  r = sqrt(xc * xc + yc * yc);
	  
	  if (inletprofile == 0) {
	    uin = InletInterpolation(r, user);
	  } else if (inletprofile == 1) {
	    uin=1.;
	  } else if (inletprofile == -1) {
	    uin = -1.;
	  } else if (inletprofile == 2) {
	    uin = 4.*Umax*yc*(H-yc)/(H*H);//InletInterpolation(r, user);
	  } else if (inletprofile == 3) {
	    uin = InletInterpolation(r, user)*AreaSumOut/AreaSumIn;
	  } else {
	    PetscPrintf(PETSC_COMM_SELF, "WRONG INLET PROFILE TYPE!!!! U_in = 0\n");
	    uin = 0.;
	  }
	  
	  ucont[k-1][j][i].z = uin*sqrt(zet[k-1][j][i].z*zet[k-1][j][i].z+
					zet[k-1][j][i].y*zet[k-1][j][i].y+
					zet[k-1][j][i].x*zet[k-1][j][i].x);

	  q[0] = 0.;
	  q[1] = 0.;
	  q[2] = ucont[k-1][j][i].z;

	  mat[0][0] = 0.5 * (csi[k-1][j][i-1].x + csi[k-1][j][i].x);
	  mat[0][1] = 0.5 * (csi[k-1][j][i-1].y + csi[k-1][j][i].y);
	  mat[0][2] = 0.5 * (csi[k-1][j][i-1].z + csi[k-1][j][i].z);
	  
	  mat[1][0] = 0.5 * (eta[k-1][j-1][i].x + eta[k-1][j][i].x);
	  mat[1][1] = 0.5 * (eta[k-1][j-1][i].y + eta[k-1][j][i].y);
	  mat[1][2] = 0.5 * (eta[k-1][j-1][i].z + eta[k-1][j][i].z);
	  
	  mat[2][0] = zet[k-1][j][i].x;
	  mat[2][1] = zet[k-1][j][i].y;
	  mat[2][2] = zet[k-1][j][i].z;

	  det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
	    mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
	    mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

	  det0 = q[0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
	    q[1] * (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]) +
	    q[2] * (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]);

	  det1 = -q[0] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
	    q[1] * (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) -
	    q[2] * (mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0]);

	  det2 = q[0] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) -
	    q[1] * (mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0]) +
	    q[2] * (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);

	  ubcs[k][j][i].x = det0 / det;
	  ubcs[k][j][i].y = det1 / det;
	  ubcs[k][j][i].z = det2 / det;

	  ucat[k-1][j][i].x = det0 / det;
	  ucat[k-1][j][i].y = det1 / det;
	  ucat[k-1][j][i].z = det2 / det;
	  
	}
      }
    }
    
    DMDAVecRestoreArray(fda, user->Ucont, &ucont);

    PetscBarrier(PETSC_NULL);


    DMDAVecGetArray(fda, user->Ucont, &ucont);    
    if (ze==mz) {
      FluxIn = 0;
      k = mz-2;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  FluxIn += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxIn = 0;
    }
    MPI_Allreduce(&FluxIn,&FluxInSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  
    
    PetscPrintf(PETSC_COMM_WORLD,"Inflow Area in out %le %le %le\n", AreaSumIn, AreaSumOut, FluxInSum);
  } // if bstype == INLET
        
  user->FluxInSum = FluxInSum;    


  
  DMDAVecRestoreArray(fda, Coor, &coor);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  
  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  return 0; 
}

PetscErrorCode OutflowFlux(UserCtx *user) {
  
  PetscInt i, j, k;
  PetscReal FluxOut;
 
  Cmpnts	***ucont, ***ucat, ***zet;

  DM fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  

  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Ucat,  &ucat);
  DMDAVecGetArray(fda, user->lZet,  &zet);
  
  FluxOut = 0;
  
  if (user->bctype[5] == 4 || user->bctype[5] == 0 || user->bctype[5] == 14) {    
    if (ze==mz) {
      k = mz-2;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  FluxOut += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxOut = 0;
    }

  } else if (user->bctype[4] == 4 || user->bctype[4] == 0) {    
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  FluxOut += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxOut = 0;
    }
  }
  MPI_Allreduce(&FluxOut,&FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 
  user->FluxOutSum = FluxOutSum;


  DMDAVecRestoreArray(fda, user->Ucont, &ucont);

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  return 0;
}

PetscErrorCode Blank_Interface(UserCtx *user) {
  PetscInt sb;
  //  PetscInt ci, cj, ck;
  //  PetscInt hi, hj, hk, hb;
  PetscInt i, j, k;
  //  Cmpnts ***itfc;
  PetscReal ***nvert, counter=0.,counterSum;

  PetscInt ip,jp,kp;
  PetscInt dispx, dispy, dispz, dispnn;

/*   for (bi=0; bi<block_number; bi++) { */
    DM		da =user->da;
    DMDALocalInfo	info = user->info;

    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    PetscInt	lxs, lxe, lys, lye, lzs, lze;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
    
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
    
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;
    
    DMDAVecGetArray(da, user->Nvert, &nvert);

/*     for (i=0; i<user->itfcptsnumber; i++) { */
/*       ci = user->itfcI[i]; */
/*       cj = user->itfcJ[i]; */
/*       ck = user->itfcK[i]; */

/*       if (ci>xs && ci<lxe && */
/* 	  cj>lys && cj<lye && */
/* 	  ck>lzs && ck<lze) { */
/* 	counter++; */
/* 	nvert[ck][cj][ci]=1.5; */

/*       } */

/*     } */
    for (sb=1; sb<block_number; sb++) { 
      ip=user->ip[sb]; jp=user->jp[sb]; kp=user->kp[sb];
      dispx=user->dispx[sb]; dispy=user->dispy[sb]; dispz=user->dispz[sb];
      dispnn=user->dispnn[sb]; 

      for(i=ip-dispx+1; i<=ip+dispx; i++){
	for (j=jp-dispy+1; j<=jp+dispy; j++){
	  for (k=kp-dispz+1; k<=kp+dispz; k++){
	    if (i>=lxs && i<lxe && j>=lys && j<lye && k>=lzs && k<lze){
	      nvert[k][j][i]=sb*10.;
	      counter+=1.;
	    }
	  }
	}
      }  
    }
    DMDAVecRestoreArray(da, user->Nvert, &nvert);
  
    MPI_Allreduce(&counter,&counterSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
   
    DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
    PetscPrintf(PETSC_COMM_WORLD, "Interface pts blanked!!!! %i ip,etc %d %d %d %d %d %d\n", counterSum,ip,jp,kp,dispx,dispy,dispz);

    return(0);
}

/* 
 INTERFACE TYPE inttype
   Normal       0
   Mixed Inlet  1
   Mixed Outlet 2
   Side no flow 3

inttype_ratio is set by inttype[6]
   0      constant velocity correction
   1      proportional (needs FluxOutSum>0)
   2      proportional to flux
   3      proportional to normal velocity (flux/area)
*/

PetscErrorCode Calc_Correction(PetscReal lFluxOutSum[6], PetscReal lAreaSum[6], 
			       PetscReal lFluxOutSum_abs[6],
			       PetscReal *ratio,PetscReal *FluxIn, PetscReal *FluxOutSum,
			       PetscReal *AreaSum, PetscInt inttype, 
			       PetscInt inttype_ratio, PetscInt surf, UserCtx *user)
{
  PetscInt      bi,i;
/*   PetscScalar	FluxIn, FluxOutSum; */
/*   PetscScalar   AreaSum; */
  PetscReal     FluxOutSum_abs=0.;
  PetscScalar   epsilon=1.e-10, sign_surf;

  if (surf%2) 
    sign_surf=1;
  else 
    sign_surf=-1.;
  bi=user->_this;

  *FluxIn = FluxInSum;
  
  if (inttype==0) {    
    *FluxIn = sign_surf*FluxInSum;
    *FluxOutSum = lFluxOutSum[surf];
    *AreaSum = lAreaSum[surf];
    FluxOutSum_abs=lFluxOutSum_abs[surf];
  } else  if (inttype==1) {
    *FluxIn = -FluxInSum;
    *FluxOutSum=0.;*AreaSum=0.;
    for (i=0; i<6; i++) {
      *FluxOutSum += lFluxOutSum[i];
      *AreaSum += lAreaSum[i];
      FluxOutSum_abs +=lFluxOutSum_abs[i];
    }
  } else if (inttype==2) {
    *FluxOutSum=0.;*AreaSum=0.;
    for (i=0; i<6; i++) {
      *FluxOutSum += lFluxOutSum[i];
      *AreaSum += lAreaSum[i];
      FluxOutSum_abs +=lFluxOutSum_abs[i];
    }
  } else if (inttype==3) {
    *FluxIn = 0.;
    
    *FluxOutSum=0.;*AreaSum=0.;
    for (i=0; i<6; i++) {
      *FluxOutSum += lFluxOutSum[i];
      *AreaSum += lAreaSum[i];
      FluxOutSum_abs +=lFluxOutSum_abs[i];
    }
  } else if (inttype==4) {
    *FluxIn = user[bi].FluxInSum;
    *FluxOutSum=0.;*AreaSum=0.;
    for (i=0; i<6; i++) {
      *FluxOutSum += lFluxOutSum[i];
      *AreaSum += lAreaSum[i];
      FluxOutSum_abs +=lFluxOutSum_abs[i];
    }
  }
    
  if (*AreaSum<epsilon) *AreaSum=1.;

  if (inttype_ratio>1) {
    if (fabs(*FluxOutSum)< 1.e-8)
      *ratio=0.;
    else
      *ratio = (*FluxIn - *FluxOutSum) / FluxOutSum_abs;
  } else if (inttype_ratio)
    *ratio = *FluxIn/ *FluxOutSum;
  else
    *ratio = (*FluxIn - *FluxOutSum) / *AreaSum;
  
  return(0);
}

/* INTERFACE TYPE inttype
   Normal       0
   Mixed Inlet  1
   Mixed Outlet 2
   Side no flow 3

   flg : assigns the correction type
   flg == 0 normal corection proportional to area (const vel)
   flg == 2 correction porportional to flux
   flg == 3 correction porportional to normal velocity
*/

PetscErrorCode Block_Blank_Correction_adv(UserCtx *user, //Vec lUcor, 					  
					  PetscInt flg) 
{
  DM	da = user->da, fda = user->fda;

  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt      i, j, k, sb;
  PetscInt      lxs, lys, lzs, lxe, lye, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  
  PetscReal epsilon=-1.e-8;
  PetscReal ***nvert, ibmval,ibm_Flux,ibm_Area;
  Cmpnts ***ucor, ***csi, ***eta, ***zet;
  DMDAVecGetArray(fda, user->Ucont, &ucor);
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  PetscReal libm_Flux, libm_area, libm_Flux_abs, ibm_Flux_abs;

  for (sb=1; sb<block_number; sb++) {
  ibmval=sb*10.-1.;

  libm_Flux = 0;
  libm_area = 0;
  libm_Flux_abs=0.;

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > ibmval-0.5 && nvert[k][j][i+1] < ibmval+0.5 && i < mx-2) {
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	    libm_Flux += ucor[k][j][i].x;
	    if (flg==3) 
	      libm_Flux_abs += fabs(ucor[k][j][i].x)/sqrt(csi[k][j][i].x * csi[k][j][i].x +
							  csi[k][j][i].y * csi[k][j][i].y +
							  csi[k][j][i].z * csi[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].x);
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
	    } else 
	      ucor[k][j][i].x=0.;
			  
	  }
	  if (nvert[k][j+1][i] > ibmval-0.5 && nvert[k][j+1][i] < ibmval+0.5 && j < my-2) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    libm_Flux += ucor[k][j][i].y;
	    if (flg==3) 
	      libm_Flux_abs += fabs(ucor[k][j][i].y)/sqrt(eta[k][j][i].x * eta[k][j][i].x +
							  eta[k][j][i].y * eta[k][j][i].y +
							  eta[k][j][i].z * eta[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].y);
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	    } else 
	      ucor[k][j][i].y=0.;
	  }
	  if (nvert[k+1][j][i] > ibmval-0.5 && nvert[k+1][j][i] < ibmval+0.5 && k < mz-2) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    libm_Flux += ucor[k][j][i].z;
	    if (flg==3)
	      libm_Flux_abs += fabs(ucor[k][j][i].z)/sqrt(zet[k][j][i].x * zet[k][j][i].x +
							  zet[k][j][i].y * zet[k][j][i].y +
							  zet[k][j][i].z * zet[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].z);
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	    }else 
	      ucor[k][j][i].z=0.;
	  }
	}

	if (nvert[k][j][i] > ibmval-0.5 && nvert[k][j][i] < ibmval+0.5) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	    libm_Flux -= ucor[k][j][i].x;
	    if (flg==3)
	    libm_Flux_abs += fabs(ucor[k][j][i].x)/sqrt(csi[k][j][i].x * csi[k][j][i].x +
							csi[k][j][i].y * csi[k][j][i].y +
							csi[k][j][i].z * csi[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].x);
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	    }else 
	      ucor[k][j][i].x=0.;
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    libm_Flux -= ucor[k][j][i].y;
	    if (flg==3)
	      libm_Flux_abs += fabs(ucor[k][j][i].y)/ sqrt(eta[k][j][i].x * eta[k][j][i].x +
							   eta[k][j][i].y * eta[k][j][i].y +
							   eta[k][j][i].z * eta[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].y);
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	    }else 
	      ucor[k][j][i].y=0.;
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    libm_Flux -= ucor[k][j][i].z;
	    if (flg==3)
	      libm_Flux_abs += fabs(ucor[k][j][i].z)/sqrt(zet[k][j][i].x * zet[k][j][i].x +
							  zet[k][j][i].y * zet[k][j][i].y +
							  zet[k][j][i].z * zet[k][j][i].z);
	    else
	      libm_Flux_abs += fabs(ucor[k][j][i].z);
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	    }else 
	      ucor[k][j][i].z=0.;
	  }
	}

      }
    }
  }

  MPI_Allreduce(&libm_Flux,&ibm_Flux,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_Flux_abs,&ibm_Flux_abs,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_area,&ibm_Area,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 /*  PetscGlobalSum(PETSC_COMM_WORLD,&libm_Flux, &ibm_Flux); */
/*   PetscGlobalSum(PETSC_COMM_WORLD,&libm_Flux_abs, &ibm_Flux_abs); */
/*   PetscGlobalSum(PETSC_COMM_WORLD,&libm_area, &ibm_Area); */

  PetscPrintf(PETSC_COMM_WORLD, "BLANKFlux %le %le\n", ibm_Flux, ibm_Area);

  PetscReal correction;

  if (ibm_Area > 1.e-15) {
    if (flg>1) 
      correction = ibm_Flux / ibm_Flux_abs;
    else if (flg)
      correction = (ibm_Flux + FluxInSum) / ibm_Area;
    else
      correction = ibm_Flux / ibm_Area;
  }
  else {
    correction = 0;
  }

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > ibmval-0.5 && nvert[k][j][i+1] < ibmval+0.5 && i < mx-2) {
	    if (fabs(ucor[k][j][i].x)>epsilon){
	    if (flg==3) 
	      ucor[k][j][i].x -=correction*fabs(ucor[k][j][i].x)/
		sqrt(csi[k][j][i].x * csi[k][j][i].x +
		     csi[k][j][i].y * csi[k][j][i].y +
		     csi[k][j][i].z * csi[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].x -=correction*fabs(ucor[k][j][i].x);
	    else
	    ucor[k][j][i].x -= sqrt(csi[k][j][i].x * csi[k][j][i].x +
				    csi[k][j][i].y * csi[k][j][i].y +
				    csi[k][j][i].z * csi[k][j][i].z) *
				    correction;
	    }
			  
	  }
	  if (nvert[k][j+1][i] > ibmval-0.5 && nvert[k][j+1][i] < ibmval+0.5 && j < my-2) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].y -=correction*fabs(ucor[k][j][i].y)/
		sqrt(eta[k][j][i].x * eta[k][j][i].x + 
		     eta[k][j][i].y * eta[k][j][i].y +
		     eta[k][j][i].z * eta[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].y -=correction*fabs(ucor[k][j][i].y);
	    else
	    ucor[k][j][i].y -= sqrt(eta[k][j][i].x * eta[k][j][i].x + 
				    eta[k][j][i].y * eta[k][j][i].y +
				    eta[k][j][i].z * eta[k][j][i].z) *
				    correction;
	    }
	  }
	  if (nvert[k+1][j][i] > ibmval-0.5 && nvert[k+1][j][i] < ibmval+0.5 && k < mz-2) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].z -= correction*fabs(ucor[k][j][i].z)/
		sqrt(zet[k][j][i].x * zet[k][j][i].x +
		     zet[k][j][i].y * zet[k][j][i].y +
		     zet[k][j][i].z * zet[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].z -= correction*fabs(ucor[k][j][i].z);
	    else
	    ucor[k][j][i].z -= sqrt(zet[k][j][i].x * zet[k][j][i].x +
				    zet[k][j][i].y * zet[k][j][i].y +
				    zet[k][j][i].z * zet[k][j][i].z) *
				    correction;
	    }
	  }
	}

	if (nvert[k][j][i] > ibmval-0.5 && nvert[k][j][i] < ibmval+0.5) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].x += correction*fabs(ucor[k][j][i].x)/
		sqrt(csi[k][j][i].x * csi[k][j][i].x +
		     csi[k][j][i].y * csi[k][j][i].y +
		     csi[k][j][i].z * csi[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].x += correction*fabs(ucor[k][j][i].x);
	    else
	    ucor[k][j][i].x += sqrt(csi[k][j][i].x * csi[k][j][i].x +
				    csi[k][j][i].y * csi[k][j][i].y +
				    csi[k][j][i].z * csi[k][j][i].z) *
				    correction;
	    }
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].y +=correction*fabs(ucor[k][j][i].y)/
		sqrt(eta[k][j][i].x * eta[k][j][i].x +
		     eta[k][j][i].y * eta[k][j][i].y +
		     eta[k][j][i].z * eta[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].y +=correction*fabs(ucor[k][j][i].y);
	    else
	    ucor[k][j][i].y += sqrt(eta[k][j][i].x * eta[k][j][i].x +
				    eta[k][j][i].y * eta[k][j][i].y +
				    eta[k][j][i].z * eta[k][j][i].z) *
				    correction;
	    }
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].z += correction*fabs(ucor[k][j][i].z)/
		sqrt(zet[k][j][i].x * zet[k][j][i].x +
		     zet[k][j][i].y * zet[k][j][i].y +
		     zet[k][j][i].z * zet[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].z += correction*fabs(ucor[k][j][i].z);
	    else
	    ucor[k][j][i].z += sqrt(zet[k][j][i].x * zet[k][j][i].x +
				    zet[k][j][i].y * zet[k][j][i].y +
				    zet[k][j][i].z * zet[k][j][i].z) *
				    correction;
	    }
	  }
	}

      }
    }
  }
  
  libm_Flux = 0;
  libm_area = 0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > ibmval-0.5 && nvert[k][j][i+1] < ibmval+0.5 && i < mx-2) {
	    libm_Flux += ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] > ibmval-0.5 && nvert[k][j+1][i] < ibmval+0.5 && j < my-2) {
	    libm_Flux += ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] > ibmval-0.5 && nvert[k+1][j][i] < ibmval+0.5 && k < mz-2) {
	    libm_Flux += ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

	if (nvert[k][j][i] > ibmval-0.5 && nvert[k][j][i] < ibmval+0.5) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    libm_Flux -= ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    libm_Flux -= ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    libm_Flux -= ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

      }
    }
  }

  MPI_Allreduce(&libm_Flux,&ibm_Flux,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_area,&ibm_Area,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 /*  PetscGlobalSum(PETSC_COMM_WORLD,&libm_Flux, &ibm_Flux); */
/*   PetscGlobalSum(PETSC_COMM_WORLD,&libm_area, &ibm_Area); */
  PetscPrintf(PETSC_COMM_WORLD, "BLANK Flux22 %le %le\n", ibm_Flux, ibm_Area);

  //  FormBCS(user);        
  } //sb

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(fda, user->Ucont, &ucor);

  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

  //DALocalToGlobal(user->fda, user->lUcont,INSERT_VALUES,user->Ucont);
  return 0;
}

/*
 intcontrol
 0 single inlet or outlet
 1 multiple inlet
 2 multiple outlet
 3 summation zero
*/
PetscErrorCode Block_Blank_Correction(UserCtx *user) {
  PetscInt      bi,sb;
  DM            da, fda;
  PetscInt      i, j, k;

  DMDALocalInfo	info ;
  PetscInt	xs, xe, ys, ye, zs, ze;
  PetscInt  	mx,my,mz;	
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts	***ucont,***csi, ***eta, ***zet;
  PetscScalar	FluxIn,ratio;
  PetscScalar	lFluxOut[6];//FluxOutKmax,FluxOutKmin,FluxOutJmax,FluxOutJmin,FluxOutImax,FluxOutImin;
  PetscScalar	lArea[6];//lAreaKmax,lAreaKmin,lAreaJmax,lAreaJmin,lAreaImax,lAreaImin;
  PetscScalar	lFluxOutSum[6]; //FluxOutKmaxSum=0.,FluxOutKminSum=0.,FluxOutJmaxSum=0.,FluxOutJminSum=0.,FluxOutImaxSum=0.,FluxOutIminSum=0.;
  PetscScalar	lAreaSum[6];//AreaKmaxSum=0.,AreaKminSum=0.,AreaJmaxSum=0.,AreaJminSum=0.,AreaImaxSum=0.,AreaIminSum=0.;
  PetscScalar   Area, AreaSum;
  PetscScalar   epsilon=-1.e-8;
 

  PetscInt ip,jp,kp;
  PetscInt dispx, dispy, dispz, dispnn;
  /**************************************************************************************************************************/
  /* Calculate the fluxIn*/
  /**************************************************************************************************************************/
  //  for (bi=0; bi<block_number; bi++) { 
  bi=user->_this;
  InflowFlux(&(user[bi]));
  //  }

  /**************************************************************************************************************************/
  /* Calculate fluxes and Correct fluxes at Interface
   ----currently only at k direction. easily extendable to other
   ----directions */
  /**************************************************************************************************************************/
 
    //  for (bi=0; bi<block_number; bi++) { 

  da = user[bi].da;
  fda = user[bi].fda;
  info = user[bi].info;
  
  xs = info.xs; xe = info.xs + info.xm;
  ys = info.ys; ye = info.ys + info.ym;
  zs = info.zs; ze = info.zs + info.zm;
  mx = info.mx; my = info.my; mz = info.mz;
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;
  
  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  
  
  DMDAVecGetArray(fda, user[bi].lCsi,  &csi);
  DMDAVecGetArray(fda, user[bi].lEta,  &eta);
  DMDAVecGetArray(fda, user[bi].lZet,  &zet);
  
  //  DMDAVecGetArray(da, user[bi].Nvert,  &nvert);
  
  for (sb=1; sb<block_number; sb++) { 
    ip=user[bi].ip[sb]; jp=user[bi].jp[sb]; kp=user[bi].kp[sb];
    dispx=user[bi].dispx[sb]; dispy=user[bi].dispy[sb]; dispz=user[bi].dispz[sb];
    dispnn=user[bi].dispnn[sb]; 
  /**************************************************************************************************************************/
  /* Calculate the fluxOut and Area on all interfaces*/
  /**************************************************************************************************************************/
    for (i=0; i<6; i++) {
      lFluxOutSum[i]=0.;
      lAreaSum[i]=0.;
    }

    /**************************************************************************************************************************/
    /* Interface at Kmax */
    /**************************************************************************************************************************/

    DMDAVecGetArray(fda, user[bi].Ucont, &ucont);
    //    DMDAVecGetArray(fda, user[bi].Bcs.Ubcs, &ubcs);
      
    lArea[5]=0.;
    lFluxOut[5] = 0.;
    k=kp+dispz;
    if (k>=lzs && k<lze){ 
      // for (i=lxs; i<lxe; i++){ 
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (j=jp-dispy+1; j<=jp+dispy; j++){ 
	  if (i>=lxs && i<lxe && j>=lys && j<lye){
	    Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			 (zet[k][j][i].y) * (zet[k][j][i].y) +
			 (zet[k][j][i].z) * (zet[k][j][i].z));
	    if (fabs(ucont[k][j][i].z/Area)>epsilon){
	      lFluxOut[5] -= ucont[k][j][i].z;
	      lArea[5] += Area;
	    }
	  }
	}
      }
    }
    
    MPI_Allreduce(&lFluxOut[5],&lFluxOutSum[5],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea[5],&lAreaSum[5],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

    //  PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[5], &lFluxOutSum[5]);
    //  PetscGlobalSum(PETSC_COMM_WORLD,&lArea[5], &lAreaSum[5]);


    /**************************************************************************************************************************/
    /* Interface at Kmin */
    /**************************************************************************************************************************/
    lArea[4]=0.;
    lFluxOut[4] = 0.;
    k=kp-dispz;
    if (k>=lzs && k<lze){
      //      for (i=lxs; i<lxe; i++){
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (j=jp-dispy+1; j<=jp+dispy; j++){
	  if (i>=lxs && i<lxe && j>=lys && j<lye){	
            Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
                         (zet[k][j][i].y) * (zet[k][j][i].y) +
                         (zet[k][j][i].z) * (zet[k][j][i].z));
            if (fabs(ucont[k][j][i].z/Area)>epsilon){
              lFluxOut[4] += ucont[k][j][i].z; 
              lArea[4] += Area; 
	    }  
	  }      
	} 
      } 
    } 

    MPI_Allreduce(&lFluxOut[4],&lFluxOutSum[4],1,MPIU_REAL,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea[4],&lAreaSum[4],1,MPIU_REAL,MPI_SUM,PETSC_COMM_WORLD);
 
    // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[4], &lFluxOutSum[4]);
    // PetscGlobalSum(PETSC_COMM_WORLD,&lArea[4], &lAreaSum[4]);
    
    /**************************************************************************************************************************/
    /* Interface at Jmax */
    /**************************************************************************************************************************/
    lArea[3]=0.;
    lFluxOut[3] = 0.;
    j=jp+dispy;
    if (j>=lys && j<lye){
      //for (i=lxs; i<lxe; i++){
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){	  
	  if (i>=lxs && i<lxe && k>=lzs && k<lze){
            Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
                         (eta[k][j][i].y) * (eta[k][j][i].y) +
                         (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (user[bi].inttype[3]==5)
	      ucont[k][j][i].y = 0.;
            if (fabs(ucont[k][j][i].y/Area)>epsilon){
              lFluxOut[3] -= ucont[k][j][i].y;
              lArea[3] += Area;
	    }	
	  }
	}
      }
    }

    MPI_Allreduce(&lFluxOut[3],&lFluxOutSum[3],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea[3],&lAreaSum[3],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

    // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[3], &lFluxOutSum[3]);
    // PetscGlobalSum(PETSC_COMM_WORLD,&lArea[3], &lAreaSum[3]);


    /**************************************************************************************************************************/
    /* Interface at Jmin */
    /**************************************************************************************************************************/
    lArea[2]=0.;
    lFluxOut[2] = 0.;
    j=jp-dispy;
    if (j>=lys && j<lye){ 
      //for (i=lxs; i<lxe; i++){
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){
	  if (i>=lxs && i<lxe && k>=lzs && k<lze){
            Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
                         (eta[k][j][i].y) * (eta[k][j][i].y) +
                         (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (user[bi].inttype[2]==5)
	      ucont[k][j][i].y = 0.;
            if (fabs(ucont[k][j][i].y/Area)>epsilon){
              lFluxOut[2] += ucont[k][j][i].y;
              lArea[2] += Area;
	    }
	  }
	} 
      }
    }

    MPI_Allreduce(&lFluxOut[2],&lFluxOutSum[2],1,MPIU_REAL,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea[2],&lAreaSum[2],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    //  PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[2], &lFluxOutSum[2]);
    //  PetscGlobalSum(PETSC_COMM_WORLD,&lArea[2], &lAreaSum[2]);
    
    /**************************************************************************************************************************/
    /* Interface at Imax */
    /**************************************************************************************************************************/
    lArea[1]=0.;
    lFluxOut[1] = 0.;
    i=ip+dispx;
    if (i>=lxs && i<lxe){
      for (j=jp-dispy+1; j<=jp+dispy; j++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){
	  if (j>=lys && j<lye && k>=lzs && k<lze){
            Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
                         (csi[k][j][i].y) * (csi[k][j][i].y) +
                         (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (user[bi].inttype[1]==5)
	      ucont[k][j][i].x = 0.;
            if (fabs(ucont[k][j][i].x/Area)>epsilon){
              lFluxOut[1] -= ucont[k][j][i].x;
              lArea[1] += Area;
	    }
	  }
	}
      }
    }

    MPI_Allreduce(&lFluxOut[1],&lFluxOutSum[1],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea[1],&lAreaSum[1],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

    //   PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[1], &lFluxOutSum[1]);
    //   PetscGlobalSum(PETSC_COMM_WORLD,&lArea[1], &lAreaSum[1]);

    /**************************************************************************************************************************/
    /* Interface at Imin */
    /**************************************************************************************************************************/
    lArea[0]=0.;
    lFluxOut[0] = 0.;
    i=ip-dispx;
    if (i>=lxs && i<lxe){
      for (j=jp-dispy+1; j<=jp+dispy; j++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){
	  if (j>=lys && j<lye && k>=lzs && k<lze){
	    Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			 (csi[k][j][i].y) * (csi[k][j][i].y) +
			 (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (user[bi].inttype[0]==5)
	      ucont[k][j][i].x = 0.;
	    if (fabs(ucont[k][j][i].x/Area)>epsilon){
	      lFluxOut[0] += ucont[k][j][i].x;
	      lArea[0] += Area;
	    }
	  }
	}
      }
    }

    MPI_Allreduce(&lFluxOut[0],&lFluxOutSum[0],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea[0],&lAreaSum[0],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[0], &lFluxOutSum[0]);
    // PetscGlobalSum(PETSC_COMM_WORLD,&lArea[0], &lAreaSum[0]);

  /**************************************************************************************************************************/
  /* Correct fluxes at Interface
   ----Extended to all directions from the previous version
   ---- */
  /**************************************************************************************************************************/
      
    FluxIn = FluxInSum;
    
    if (user[bi].inttype[5]==0) {
      FluxOutSum = lFluxOutSum[5];
      AreaSum = lAreaSum[5];
    } else  if (user[bi].inttype[5]==1) {
      FluxIn = -FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[5]==2) {
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[5]==3) {
      FluxIn = 0.;
      
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[5]==4) {
      FluxIn = user[bi].FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }
    
    if (AreaSum<epsilon) AreaSum=1.;
    if (user[bi].inttype[6])
      ratio = FluxIn/FluxOutSum;
    else 
      ratio = (FluxIn - FluxOutSum) / AreaSum;
    PetscPrintf(PETSC_COMM_WORLD, "Blank Ratio Kmax Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[5], lAreaSum[5]);

    k=kp+dispz;
    if (k>=lzs && k<lze){
      //      for (i=lxs; i<lxe; i++){
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (j=jp-dispy+1; j<=jp+dispy; j++){
	  if (i>=lxs && i<lxe && j>=lys && j<lye){
	    Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			 (zet[k][j][i].y) * (zet[k][j][i].y) +
			 (zet[k][j][i].z) * (zet[k][j][i].z));
	    if (fabs(ucont[k][j][i].z/Area)>epsilon){
	      ucont[k][j][i].z -= ratio* Area;
	    } else
	      ucont[k][j][i].z = 0.;
	  }
	}
      }
    }
    //////////// Kmax End
    
    //////////// Kmin begin    
    FluxIn = FluxInSum;

    if (user[bi].inttype[4]==0) {
      FluxIn = -FluxInSum;
      FluxOutSum = lFluxOutSum[4];
      AreaSum = lAreaSum[4];
    } else  if (user[bi].inttype[4]==1) {
      FluxIn = -FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[4]==2){
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[4]==3) {
      FluxIn = 0.;

      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }else if (user[bi].inttype[4]==4) {
      FluxIn = user[bi].FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }


    if (AreaSum<epsilon) AreaSum=1.;
    if (user[bi].inttype[6])
      ratio = FluxIn/FluxOutSum;
    else      
      ratio = (FluxIn - FluxOutSum) / AreaSum;
    PetscPrintf(PETSC_COMM_WORLD, "Blank Ratio Kmin Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[4], lAreaSum[4]);
    k=kp-dispz;
    if (k>=lzs && k<lze){
      //      for (i=lxs; i<lxe; i++){
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (j=jp-dispy+1; j<=jp+dispy; j++){
	  if (i>=lxs && i<lxe && j>=lys && j<lye){
            Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
                         (zet[k][j][i].y) * (zet[k][j][i].y) +
                         (zet[k][j][i].z) * (zet[k][j][i].z));
            if (fabs(ucont[k][j][i].z/Area)>epsilon){
              ucont[k][j][i].z += ratio* Area;
            } else
              ucont[k][j][i].z = 0.;
	  }
	}
      } 
    }
      
    // kmin 

    // jmax
    FluxIn = FluxInSum;
    
    if (user[bi].inttype[3]==0) {
      FluxOutSum = lFluxOutSum[3];
      AreaSum = lAreaSum[3];
    } else  if (user[bi].inttype[3]==1) {
      FluxIn = -FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[3]==2){
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[3]==3) {
      FluxIn = 0.;
      
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[3]==4) {
      FluxIn = user[bi].FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } 


    if (AreaSum<epsilon) AreaSum=1.;
    if (user[bi].inttype[6])
      ratio = FluxIn/FluxOutSum;
    else 
      ratio = (FluxIn - FluxOutSum) / AreaSum;
    PetscPrintf(PETSC_COMM_WORLD, "Blank Ratio Jmax Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[3], lAreaSum[3]);
    j=jp+dispy;
    if (j>=lys && j<lye){
      //      for (i=lxs; i<lxe; i++){
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){
	  if (i>=lxs && i<lxe && k>=lzs && k<lze){
            Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
                         (eta[k][j][i].y) * (eta[k][j][i].y) +
                         (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (user[bi].inttype[3]==5){
	      ucont[k][j][i].y = 0.;
	    } else if (fabs(ucont[k][j][i].y/Area)>epsilon){
              ucont[k][j][i].y -= ratio * Area;
            } else
              ucont[k][j][i].y = 0.;    
	  }
	}
      }
    }
    // jmax

    // jmin
    FluxIn = FluxInSum;
    
    if (user[bi].inttype[2]==0) {
      FluxIn = -FluxInSum;
      FluxOutSum = lFluxOutSum[2];
      AreaSum = lAreaSum[2];
    } else  if (user[bi].inttype[2]==1) {
      FluxIn = -FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }else if (user[bi].inttype[2]==2){
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }else if (user[bi].inttype[2]==3) {
      FluxIn = 0.;
      
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[2]==4) {
      FluxIn = user[bi].FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }
    
    if (AreaSum<epsilon) AreaSum=1.;
    if (user[bi].inttype[6])
      ratio = FluxIn/FluxOutSum;
    else 
      ratio = (FluxIn - FluxOutSum) / AreaSum;
    PetscPrintf(PETSC_COMM_WORLD, "Blank Ratio Jmin Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[2], lAreaSum[2]);
    j=jp-dispy;
    if (j>=lys && j<lye){
      //      for (i=lxs; i<lxe; i++){
      for (i=ip-dispx+1; i<=ip+dispx; i++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){
	  if (i>=lxs && i<lxe && k>=lzs && k<lze){
            Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
                         (eta[k][j][i].y) * (eta[k][j][i].y) +
                         (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (user[bi].inttype[2]==5){
	      ucont[k][j][i].y = 0.;
            } else if (fabs(ucont[k][j][i].y/Area)>epsilon){
              ucont[k][j][i].y += ratio * Area;
            } else
              ucont[k][j][i].y = 0.;
	  }
	}
      }
    }    
    // jmin

    // imax
    FluxIn = FluxInSum;
    
    if (user[bi].inttype[1]==0) {
      FluxOutSum = lFluxOutSum[1];
      AreaSum = lAreaSum[1];
    } else  if (user[bi].inttype[1]==1) {
      FluxIn = -FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[1]==2){
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }else if (user[bi].inttype[1]==3) {
      FluxIn = 0.;
      
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[1]==4) {
      FluxIn = user[bi].FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    }
    
    
    if (AreaSum<epsilon) AreaSum=1.;

    if (user[bi].inttype[6])
      ratio = FluxIn/FluxOutSum;
    else 
      ratio = (FluxIn - FluxOutSum) / AreaSum;
 
    PetscPrintf(PETSC_COMM_WORLD, "Blank Ratio Imax Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[1], lAreaSum[1]);
    i=ip+dispx;
    if (i>=lxs && i<lxe){ 
      for (j=jp-dispy+1; j<=jp+dispy; j++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){
	  if (j>=lys && j<lye && k>=lzs && k<lze){
            Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
                         (csi[k][j][i].y) * (csi[k][j][i].y) +
                         (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (user[bi].inttype[1]==5){
	      ucont[k][j][i].x = 0.;
/* 	    }else if (ratio*lFluxOutSum[1]< 0 && */
/* 		      fabs(ratio*lAreaSum[1])>fabs(lFluxOutSum[1])){ */
/* 	      ucont[k][j][i].x = 0.; */
	    }
            else if (fabs(ucont[k][j][i].x/Area)>epsilon){
              ucont[k][j][i].x -= ratio * Area;
	    } else
	      ucont[k][j][i].x = 0.;
	  }
	}
      }
    }    
    // imax 

    // imin
    FluxIn = FluxInSum;
    
    if (user[bi].inttype[0]==0) {
      FluxIn = -FluxInSum;
      FluxOutSum = lFluxOutSum[0];
      AreaSum = lAreaSum[0];
    } else  if (user[bi].inttype[0]==1) {
      FluxIn = -FluxInSum;
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[0]==2) {
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[0]==3) {
      FluxIn = 0.;
      
      FluxOutSum=0.;AreaSum=0.;
      for (i=0; i<6; i++) {
	FluxOutSum += lFluxOutSum[i];
	AreaSum += lAreaSum[i];
      }
    } else if (user[bi].inttype[0]==4) {
      FluxIn = user[bi].FluxInSum;
	FluxOutSum=0.;AreaSum=0.;
	for (i=0; i<6; i++) {
	  FluxOutSum += lFluxOutSum[i];
	  AreaSum += lAreaSum[i];
	}
    }
    
    
    if (AreaSum<epsilon) AreaSum=1.;
    if (user[bi].inttype[6])
      ratio = FluxIn/FluxOutSum;
    else
      ratio = (FluxIn - FluxOutSum) / AreaSum;
    PetscPrintf(PETSC_COMM_WORLD, "Blank Ratio Imin Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[0], lAreaSum[0]);
    
    i=ip-dispx;
    if (i>=lxs && i<lxe){
      for (j=jp-dispy+1; j<=jp+dispy; j++){
	for (k=kp-dispz+1; k<=kp+dispz; k++){
	  if (j>=lys && j<lye && k>=lzs && k<lze){
            Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
                         (csi[k][j][i].y) * (csi[k][j][i].y) +
                         (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (user[bi].inttype[0]==5){
	      ucont[k][j][i].x = 0.;
/* 	    } else if (ratio*lFluxOutSum[0]< 0 && */
/* 		      fabs(ratio*lAreaSum[0])>fabs(lFluxOutSum[0])){ */
/* 	      ucont[k][j][i].x = 0.; */
	    }
	    else if (fabs(ucont[k][j][i].x/Area)>epsilon){
              ucont[k][j][i].x += ratio * Area;
	    } else
	      ucont[k][j][i].x = 0.;
	  }
	}
        }
    }
    
    // imin 

    /**************************************************************************************************************************/
    /* Interface at  */
    /**************************************************************************************************************************/

    DMDAVecRestoreArray(fda, user[bi].Ucont, &ucont);
    //    DMDAVecRestoreArray(fda, user[bi].Bcs.Ubcs, &ubcs);
/*     DALocalToGlobalBegin(fda, user[bi].lUcont,  user[bi].Ucont); */
/*     DALocalToGlobalEnd(fda, user[bi].lUcont,  user[bi].Ucont); */
    DMGlobalToLocalBegin(fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
    DMGlobalToLocalEnd(fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
    
    } //for sb

    DMDAVecRestoreArray(fda, user[bi].lCsi,  &csi);
    DMDAVecRestoreArray(fda, user[bi].lEta,  &eta);
    DMDAVecRestoreArray(fda, user[bi].lZet,  &zet);

    //    DMDAVecRestoreArray(da, user[bi].Nvert,  &nvert);

    Contra2Cart(&(user[bi]));

    //    GhostNodeVelocity(&user[bi]);    

    //    FormBCS(&(user[bi]));        
    // } //bi
  return(0);
}

PetscErrorCode Block_Interface_Correction(UserCtx *user) {
  PetscInt      bi;
  DM            da, fda;
  PetscInt      i, j, k;

  DMDALocalInfo	info ;
  PetscInt	xs, xe, ys, ye, zs, ze;
  PetscInt  	mx,my,mz;	
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts	***ucont,***csi, ***eta, ***zet,***ucat,***ubcs;
  PetscReal     ***nvert;
  PetscScalar	FluxIn,ratio;
  PetscScalar	lFluxOut[6],lFluxOut_abs[6];//FluxOutKmax,FluxOutKmin,FluxOutJmax,FluxOutJmin,FluxOutImax,FluxOutImin;
  PetscScalar	lArea[6];//lAreaKmax,lAreaKmin,lAreaJmax,lAreaJmin,lAreaImax,lAreaImin;
  PetscScalar	lFluxOutSum[6], lFluxOutSum_abs[6]; //FluxOutKmaxSum=0.,FluxOutKminSum=0.,FluxOutJmaxSum=0.,FluxOutJminSum=0.,FluxOutImaxSum=0.,FluxOutIminSum=0.;
  PetscScalar	lAreaSum[6];//AreaKmaxSum=0.,AreaKminSum=0.,AreaJmaxSum=0.,AreaJminSum=0.,AreaImaxSum=0.,AreaIminSum=0.;
  PetscScalar   Area, AreaSum;
  PetscScalar   epsilon=-1.e-8;
 
  /**************************************************************************************************************************/
  /* Calculate the fluxIn*/
  /**************************************************************************************************************************/
  for (bi=0; bi<block_number; bi++) { 
    InflowFlux(&(user[bi]));
  }

  /**************************************************************************************************************************/
  /* Calculate fluxes and Correct fluxes at Interface
   ----currently only at k direction. easily extendable to other
   ----directions */
  /**************************************************************************************************************************/
 
  for (bi=0; bi<block_number; bi++) { 

    da = user[bi].da;
    fda = user[bi].fda;
    info = user[bi].info;

    xs = info.xs; xe = info.xs + info.xm;
    ys = info.ys; ye = info.ys + info.ym;
    zs = info.zs; ze = info.zs + info.zm;
    mx = info.mx; my = info.my; mz = info.mz;
    
    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
    
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
    
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

 
    DMDAVecGetArray(fda, user[bi].lCsi,  &csi);
    DMDAVecGetArray(fda, user[bi].lEta,  &eta);
    DMDAVecGetArray(fda, user[bi].lZet,  &zet);

    DMDAVecGetArray(da, user[bi].Nvert,  &nvert);
  /**************************************************************************************************************************/
  /* Calculate the fluxOut and Area on all interfaces*/
  /**************************************************************************************************************************/
    for (i=0; i<6; i++) {
      lFluxOutSum[i]=0.;
      lFluxOutSum_abs[i]=0.;
      lAreaSum[i]=0.;
      lArea[i]=0.; 
      lFluxOut[i]=0.;
      lFluxOut_abs[i]=0.;
    }

    /**************************************************************************************************************************/
    /* Interface at Kmax */
    /**************************************************************************************************************************/
    if (user[bi].bctype[5]==14) {

      DMDAVecGetArray(fda, user[bi].Ucont, &ucont);
      DMDAVecGetArray(fda, user[bi].Bcs.Ubcs, &ubcs);
      DMDAVecGetArray(fda, user[bi].lUcat, &ucat);
      
      lArea[5]=0.;
      lFluxOut[5] = 0.;
      if (ze == mz) {
	k=ze-1;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (zet[k-1][j][i].x) * (zet[k-1][j][i].x) +
			 (zet[k-1][j][i].y) * (zet[k-1][j][i].y) +
			 (zet[k-1][j][i].z) * (zet[k-1][j][i].z));

	    ubcs[k][j][i].x = ucat[k-1][j][i].x;//+ratio;
	    ubcs[k][j][i].y = ucat[k-1][j][i].y;
	    ubcs[k][j][i].z = ucat[k-1][j][i].z;// + ratio;//*n_z;

	    ucont[k-1][j][i].z = (ubcs[k][j][i].x * (zet[k-1][j][i].x ) +
				  ubcs[k][j][i].y * (zet[k-1][j][i].y ) +
				  ubcs[k][j][i].z * (zet[k-1][j][i].z ));

	    lFluxOut[5] += ucont[k-1][j][i].z;
	    lArea[5] += Area;
	  }
	}
      }
      
      DMDAVecRestoreArray(fda, user[bi].Ucont, &ucont);
      DMDAVecRestoreArray(fda, user[bi].Bcs.Ubcs, &ubcs);
      DMDAVecRestoreArray(fda, user[bi].lUcat, &ucat);
      
      MPI_Allreduce(&lFluxOut[5],&lFluxOutSum[5],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lArea[5],&lAreaSum[5],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);      
      // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[5], &lFluxOutSum[5]);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lArea[5], &lAreaSum[5]);
    }

    DMDAVecGetArray(fda, user[bi].Ucont, &ucont);
    //    DMDAVecGetArray(fda, user[bi].Bcs.Ubcs, &ubcs);

    if (user[bi].bctype[5]==0) {
      lArea[5]=0.;
      lFluxOut[5] = 0.;
      if (ze == mz) {
	k=ze-2;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			 (zet[k][j][i].y) * (zet[k][j][i].y) +
			 (zet[k][j][i].z) * (zet[k][j][i].z));
	    if (fabs(ucont[k][j][i].z/Area)>epsilon 
		&& nvert[k][j][i]<0.1){
	      lFluxOut[5] += ucont[k][j][i].z;
	      if (user[bi].inttype[6]==2) 
		lFluxOut_abs[5] += fabs(ucont[k][j][i].z);
	      else if (user[bi].inttype[6]==3)
		lFluxOut_abs[5] += fabs(ucont[k][j][i].z)/Area;
	      lArea[5] += Area;
	    } 
	  }
	}
      }

      MPI_Allreduce(&lFluxOut[5],&lFluxOutSum[5],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lFluxOut_abs[5],&lFluxOutSum_abs[5],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lArea[5],&lAreaSum[5],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[5], &lFluxOutSum[5]);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut_abs[5], &lFluxOutSum_abs[5]);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lArea[5], &lAreaSum[5]);
    }

    /**************************************************************************************************************************/
    /* Interface at Kmin */
    /**************************************************************************************************************************/

    if (user[bi].bctype[4]==0) {
      if (zs == 0) {	
	k=0;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			 (zet[k][j][i].y) * (zet[k][j][i].y) +
			 (zet[k][j][i].z) * (zet[k][j][i].z));
	    if (fabs(ucont[k][j][i].z/Area)>epsilon
		&& nvert[k+1][j][i]<0.1){
	      lFluxOut[4] -= ucont[k][j][i].z;
	      if (user[bi].inttype[6]==2) 
		lFluxOut_abs[4] += fabs(ucont[k][j][i].z);
	      else if (user[bi].inttype[6]==3)
		lFluxOut_abs[4] += fabs(ucont[k][j][i].z)/Area;
	      lArea[4] += Area;
	    } 
	  }
	}
      }    
      MPI_Allreduce(&lFluxOut[4],&lFluxOutSum[4],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lFluxOut_abs[4],&lFluxOutSum_abs[4],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lArea[4],&lAreaSum[4],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);  
      // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[4], &lFluxOutSum[4]);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut_abs[4], &lFluxOutSum_abs[4]);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lArea[4], &lAreaSum[4]);
    }

    /**************************************************************************************************************************/
    /* Interface at Jmax */
    /**************************************************************************************************************************/

    if (user[bi].bctype[3]==0) {
      lArea[3]=0.;
      lFluxOut[3] = 0.;
      if (ye == my) {
	j=ye-2;
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			 (eta[k][j][i].y) * (eta[k][j][i].y) +
			 (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (fabs(ucont[k][j][i].y/Area)>epsilon
		&& nvert[k][j][i]<0.1){
	      lFluxOut[3] += ucont[k][j][i].y;
	      if (user[bi].inttype[6]==2) 
		lFluxOut_abs[3] += fabs(ucont[k][j][i].y);
	      else if (user[bi].inttype[6]==3)
		lFluxOut_abs[3] += fabs(ucont[k][j][i].y)/Area;
	      lArea[3] += Area;
	    }
	  }
	}
      }
      MPI_Allreduce(&lFluxOut[3],&lFluxOutSum[3],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lFluxOut_abs[3],&lFluxOutSum_abs[3],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lArea[3],&lAreaSum[3],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[3], &lFluxOutSum[3]);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut_abs[3], &lFluxOutSum_abs[3]);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lArea[3], &lAreaSum[3]);
    }

    /**************************************************************************************************************************/
    /* Interface at Jmin */
    /**************************************************************************************************************************/

    if (user[bi].bctype[2]==0) {
      lArea[2]=0.;
      lFluxOut[2] = 0.;
      if (ys == 0) {
	j=0;
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			 (eta[k][j][i].y) * (eta[k][j][i].y) +
			 (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (fabs(ucont[k][j][i].y/Area)>epsilon
		&& nvert[k][j+1][i]<0.1){
	      lFluxOut[2] -= ucont[k][j][i].y;
	      if (user[bi].inttype[6]==2) 
		lFluxOut_abs[2] += fabs(ucont[k][j][i].y);
	      else if (user[bi].inttype[6]==3)
		lFluxOut_abs[2] += fabs(ucont[k][j][i].y)/Area;
	      lArea[2] += Area;
	    } 
	  }
	}
      }
      MPI_Allreduce(&lFluxOut[2],&lFluxOutSum[2],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lFluxOut_abs[2],&lFluxOutSum_abs[2],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lArea[2],&lAreaSum[2],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[2], &lFluxOutSum[2]);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut_abs[2], &lFluxOutSum_abs[2]);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lArea[2], &lAreaSum[2]);
    }

    /**************************************************************************************************************************/
    /* Interface at Imax */
    /**************************************************************************************************************************/

    if (user[bi].bctype[1]==0) {
      lArea[1]=0.;
      lFluxOut[1] = 0.;
      if (xe == mx) {
	i=xe-2;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			 (csi[k][j][i].y) * (csi[k][j][i].y) +
			 (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (fabs(ucont[k][j][i].x/Area)>epsilon
		&& nvert[k][j][i]<0.1){
	      lFluxOut[1] += ucont[k][j][i].x;
	      if (user[bi].inttype[6]==2) 
		lFluxOut_abs[1] += fabs(ucont[k][j][i].x);
	      else if (user[bi].inttype[6]==3)
		lFluxOut_abs[1] += fabs(ucont[k][j][i].x)/Area;
	      lArea[1] += Area;
	    } 
	  }
	}
      }
      MPI_Allreduce(&lFluxOut[1],&lFluxOutSum[1],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lFluxOut_abs[1],&lFluxOutSum_abs[1],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lArea[1],&lAreaSum[1],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[1], &lFluxOutSum[1]);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut_abs[1], &lFluxOutSum_abs[1]);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lArea[1], &lAreaSum[1]);
    }

    /**************************************************************************************************************************/
    /* Interface at Imin */
    /**************************************************************************************************************************/

    if (user[bi].bctype[0]==0) {
      lArea[0]=0.;
      lFluxOut[0] = 0.;
      if (xs == 0) {
	i=0;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			 (csi[k][j][i].y) * (csi[k][j][i].y) +
			 (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (fabs(ucont[k][j][i].x/Area)>epsilon
		&& nvert[k][j][i+1]<0.1){
	      lFluxOut[0] -= ucont[k][j][i].x;
	      if (user[bi].inttype[6]==2) 
		lFluxOut_abs[0] += fabs(ucont[k][j][i].x);
	      else if (user[bi].inttype[6]==3)
		lFluxOut_abs[0] += fabs(ucont[k][j][i].x)/Area;
	      lArea[0] += Area;
	    }
	  }
	}
      }

      MPI_Allreduce(&lFluxOut[0],&lFluxOutSum[0],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lFluxOut_abs[0],&lFluxOutSum_abs[0],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lArea[0],&lAreaSum[0],1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut[0], &lFluxOutSum[0]);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut_abs[0], &lFluxOutSum_abs[0]);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lArea[0], &lAreaSum[0]);
    }

  /**************************************************************************************************************************/
  /* Correct fluxes at Interface
   ----Extended to all directions from the previous version
   ---- */
  /**************************************************************************************************************************/
    if (user[bi].bctype[5]==0 || user[bi].bctype[5]==14) {
      Calc_Correction(lFluxOutSum, lAreaSum, lFluxOutSum_abs, &ratio, &FluxIn, &FluxOutSum, &AreaSum, 
		      user[bi].inttype[5], user[bi].inttype[6], 5, &user[bi]);

      PetscPrintf(PETSC_COMM_WORLD, "Ratio Kmax Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[5], lAreaSum[5]);
      //user[bi].FluxIntfcSum+=ratio*user[bi].AreaIntfcSum;

      if (ze==mz) {
	k = ze-2;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			 (zet[k][j][i].y) * (zet[k][j][i].y) +
			 (zet[k][j][i].z) * (zet[k][j][i].z));
	    if (fabs(ucont[k][j][i].z/Area)>epsilon
		&& nvert[k][j][i]<0.1){
	      if (user[bi].inttype[6]==3)
		ucont[k][j][i].z += ratio * fabs(ucont[k][j][i].z)/Area;
	      else if (user[bi].inttype[6]==2)
		ucont[k][j][i].z += ratio * fabs(ucont[k][j][i].z);
	      else if (user[bi].inttype[6])
		ucont[k][j][i].z *= ratio; 
	      else
		ucont[k][j][i].z += ratio* Area; 
/* 	      ubcs[k+1][j][i].x += ratio*zet[k][j][i].x/Area;  */
/* 	      ubcs[k+1][j][i].y += ratio*zet[k][j][i].y/Area;  */
/* 	      ubcs[k+1][j][i].z += ratio*zet[k][j][i].z/Area;  */
	    } else 
	      ucont[k][j][i].z = 0.;
	  }
	}
      }

      if (user[bi].bctype[5]==14) {
	FluxOutSum=0.;AreaSum=0.;
	for (i=0; i<6; i++) {
	  if (user[bi].bctype[i]==0) {
	    if( ( i % 2 ) == 0) 
	      FluxOutSum += lFluxOutSum[i]-ratio*lAreaSum[i];
	    else 
	      FluxOutSum += lFluxOutSum[i]+ratio*lAreaSum[i];
	    AreaSum += lAreaSum[i];
	  }
	}
	user[bi].FluxIntfcSum=FluxOutSum;
	user[bi].AreaIntfcSum=AreaSum;      
      }
      
    } // kmax bctype[5]== INTERFACE

    if (user[bi].bctype[4]==0) {
      Calc_Correction(lFluxOutSum, lAreaSum, lFluxOutSum_abs, &ratio, &FluxIn, &FluxOutSum, &AreaSum, 
		      user[bi].inttype[4], 
		      user[bi].inttype[6],4, &user[bi]);

      PetscPrintf(PETSC_COMM_WORLD, "Ratio Kmin Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[4], lAreaSum[4]);
      //      user[bi].FluxIntfcSum+=ratio*user[bi].AreaIntfcSum;
      
      if (zs == 0) {
	k = 0;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			 (zet[k][j][i].y) * (zet[k][j][i].y) +
			 (zet[k][j][i].z) * (zet[k][j][i].z));
	    if (fabs(ucont[k][j][i].z/Area)>epsilon
		&& nvert[k+1][j][i]<0.1){
	      if (user[bi].inttype[6]==3)
		ucont[k][j][i].z -= ratio * fabs(ucont[k][j][i].z)/Area;
	      else if (user[bi].inttype[6]==2)
		ucont[k][j][i].z -= ratio * fabs(ucont[k][j][i].z);
	      else if (user[bi].inttype[6])
		ucont[k][j][i].z *= ratio; 
	      else
		ucont[k][j][i].z -= ratio*Area;
	      //ucont[k][j][i].z *= ratio;
/* 	      ubcs[k][j][i].x -= ratio*zet[k][j][i].x/Area;  */
/* 	      ubcs[k][j][i].y -= ratio*zet[k][j][i].y/Area;  */
/* 	      ubcs[k][j][i].z -= ratio*zet[k][j][i].z/Area;  */
	    } else 
	      ucont[k][j][i].z = 0.;	 
	  }
	}
      }
      
    } // kmin bctype[4]== INTERFACE

    if (user[bi].bctype[3]==0) {
      Calc_Correction(lFluxOutSum, lAreaSum, lFluxOutSum_abs, &ratio, &FluxIn, &FluxOutSum, &AreaSum,
		      user[bi].inttype[3], 
		      user[bi].inttype[6],3, &user[bi]);
      PetscPrintf(PETSC_COMM_WORLD, "Ratio Jmax Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[3], lAreaSum[3]);
      //      user[bi].FluxIntfcSum+=ratio*user[bi].AreaIntfcSum;

      if (ye==my) {
	j = ye-2;
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			 (eta[k][j][i].y) * (eta[k][j][i].y) +
			 (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (fabs(ucont[k][j][i].y/Area)>epsilon
		&& nvert[k][j][i]<0.1){
	      if (user[bi].inttype[6]==3)
		ucont[k][j][i].y += ratio * fabs(ucont[k][j][i].y)/Area;
	      else if (user[bi].inttype[6]==2)
		ucont[k][j][i].y += ratio * fabs(ucont[k][j][i].y);
	      else if (user[bi].inttype[6])
		ucont[k][j][i].y *= ratio; 
	      else
	      ucont[k][j][i].y += ratio * Area;
/* 	      ubcs[k][j+1][i].x += ratio*eta[k][j][i].x/Area;  */
/* 	      ubcs[k][j+1][i].y += ratio*eta[k][j][i].y/Area;  */
/* 	      ubcs[k][j+1][i].z += ratio*eta[k][j][i].z/Area;  */
	    } else 
	      ucont[k][j][i].y = 0.;
	  }
	}
      }      
    } // jmax bctype[3]== INTERFACE

    if (user[bi].bctype[2]==0) {      
      Calc_Correction(lFluxOutSum, lAreaSum, lFluxOutSum_abs,&ratio, &FluxIn, &FluxOutSum, &AreaSum,
		      user[bi].inttype[2], 
		      user[bi].inttype[6],2, &user[bi]);
      PetscPrintf(PETSC_COMM_WORLD, "Ratio Jmin Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[2], lAreaSum[2]);
      
      //      user[bi].FluxIntfcSum+=ratio*user[bi].AreaIntfcSum;

      if (ys == 0) {
	j = 0;
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    Area = sqrt( (eta[k][j][i].x) * (eta[k][j][i].x) +
			 (eta[k][j][i].y) * (eta[k][j][i].y) +
			 (eta[k][j][i].z) * (eta[k][j][i].z));
	    if (fabs(ucont[k][j][i].y/Area)>epsilon
		&& nvert[k][j+1][i]<0.1){
	      if (user[bi].inttype[6]==3)
		ucont[k][j][i].y -= ratio * fabs(ucont[k][j][i].y)/Area;
	      else if (user[bi].inttype[6]==2)
		ucont[k][j][i].y -= ratio * fabs(ucont[k][j][i].y);
	      else if (user[bi].inttype[6])
		ucont[k][j][i].y *= ratio; 
	      else
		ucont[k][j][i].y -= ratio * Area;
/* 	      ubcs[k][j][i].x -= ratio*eta[k][j][i].x/Area;  */
/* 	      ubcs[k][j][i].y -= ratio*eta[k][j][i].y/Area;  */
/* 	      ubcs[k][j][i].z -= ratio*eta[k][j][i].z/Area;  */
	    } else
	      ucont[k][j][i].y = 0.;
	  }
	}
      }      
    } // jmin bctype[2]== INTERFACE


    if (user[bi].bctype[1]==0) {
      Calc_Correction(lFluxOutSum, lAreaSum,lFluxOutSum_abs, &ratio, &FluxIn, &FluxOutSum, &AreaSum,
		      user[bi].inttype[1], 
		      user[bi].inttype[6],1, &user[bi]);
      PetscPrintf(PETSC_COMM_WORLD, "Ratio Imax Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[1], lAreaSum[1]);

      //      user[bi].FluxIntfcSum+=ratio*user[bi].AreaIntfcSum;
      
      if (xe == mx) {
	i=xe-2;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			 (csi[k][j][i].y) * (csi[k][j][i].y) +
			 (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (fabs(ucont[k][j][i].x/Area)>epsilon
		&& nvert[k][j][i]<0.1){
	      if (user[bi].inttype[6]==3)
		ucont[k][j][i].x += ratio * fabs(ucont[k][j][i].x)/Area;
	      else if (user[bi].inttype[6]==2)
		ucont[k][j][i].x += ratio * fabs(ucont[k][j][i].x);
	      else if (user[bi].inttype[6])
		ucont[k][j][i].x *= ratio; 
	      else
		ucont[k][j][i].x += ratio * Area;
/* 	      ubcs[k][j][i+1].x += ratio*csi[k][j][i].x/Area;  */
/* 	      ubcs[k][j][i+1].y += ratio*csi[k][j][i].y/Area;  */
/* 	      ubcs[k][j][i+1].z += ratio*csi[k][j][i].z/Area;  */
	    } else 
	      ucont[k][j][i].x = 0.;
	  }
	}
      }      
    } // imax bctype[1]== INTERFACE

    if (user[bi].bctype[0]==0) {
      Calc_Correction(lFluxOutSum, lAreaSum,lFluxOutSum_abs, &ratio, &FluxIn, &FluxOutSum, &AreaSum,
		      user[bi].inttype[0], 
		      user[bi].inttype[6],0, &user[bi]);
      PetscPrintf(PETSC_COMM_WORLD, "Ratio Imin Interface %d %d %le %le %le %le local %le %le\n", bi, ti, ratio, FluxIn, FluxOutSum, AreaSum, lFluxOutSum[0], lAreaSum[0]);

      //      user[bi].FluxIntfcSum+=ratio*user[bi].AreaIntfcSum;

      if (xs == 0) {
	i=0;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    Area = sqrt( (csi[k][j][i].x) * (csi[k][j][i].x) +
			 (csi[k][j][i].y) * (csi[k][j][i].y) +
			 (csi[k][j][i].z) * (csi[k][j][i].z));
	    if (fabs(ucont[k][j][i].x/Area)>epsilon
		&& nvert[k][j][i]<0.1){
	      if (user[bi].inttype[6]==3)
		ucont[k][j][i].x -= ratio * fabs(ucont[k][j][i].x)/Area;
	      else if (user[bi].inttype[6]==2)
		ucont[k][j][i].x -= ratio * fabs(ucont[k][j][i].x);
	      else if (user[bi].inttype[6])
		ucont[k][j][i].x *= ratio; 
	      else 
		ucont[k][j][i].x -= ratio * Area;
/* 	      ubcs[k][j][i].x -= ratio*csi[k][j][i].x/Area;  */
/* 	      ubcs[k][j][i].y -= ratio*csi[k][j][i].y/Area;  */
/* 	      ubcs[k][j][i].z -= ratio*csi[k][j][i].z/Area;  */
	    } else 
	      ucont[k][j][i].x = 0.;
	  }
	}
      }      
    } // imin bctype[0]== INTERFACE


    /**************************************************************************************************************************/
    /* Interface at  */
    /**************************************************************************************************************************/

    DMDAVecRestoreArray(fda, user[bi].Ucont, &ucont);
    //    DMDAVecRestoreArray(fda, user[bi].Bcs.Ubcs, &ubcs);
    DMDAVecRestoreArray(fda, user[bi].lCsi,  &csi);
    DMDAVecRestoreArray(fda, user[bi].lEta,  &eta);
    DMDAVecRestoreArray(fda, user[bi].lZet,  &zet);

    DMDAVecRestoreArray(da, user[bi].Nvert,  &nvert);

    DMGlobalToLocalBegin(fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
    DMGlobalToLocalEnd(fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

    Contra2Cart(&(user[bi]));

    GhostNodeVelocity(&user[bi]);    

    //    FormBCS(&(user[bi]));        
  }//bi
  return(0);
}

PetscErrorCode Block_Interface_U(UserCtx *user) {
  PetscInt bi,sb;
  PetscInt ci, cj, ck;
  PetscInt hi, hj, hk, hb;
  PetscInt i, j, k;
  PetscReal x, y, z;
  Vec	hostU;
  Cmpnts ***itfc,***ubcs, ucent;
  PetscReal *hostu,***nvert;

  VecScatter tolocalall;
  PetscErrorCode ierr;
    
  /* ucat is at cell centers while litfc is now on the cell corners! */
  
  for (bi=0; bi<block_number; bi++) {
    /* hostU is a parallel PETSc vector that will hold vector values 
       in the natural numbering, rather than in the PETSc parallel 
       numbering associated with the DA */
    ierr = DMDACreateNaturalVector(user[bi].fda, &hostU);
    

    // put the center node velocities in the natural ordering in nhostU 
    DMDAGlobalToNaturalBegin(user[bi].fda, user[bi].Ucat, INSERT_VALUES, hostU);
    DMDAGlobalToNaturalEnd(user[bi].fda, user[bi].Ucat, INSERT_VALUES, hostU);

    /* the velocities at cell centers are saved in user[bi].nhostU
       which is a sequential vector*/
   
    VecScatterCreateToAll(hostU, &tolocalall, &(user[bi].nhostU));

    VecScatterBegin(tolocalall, hostU, user[bi].nhostU, INSERT_VALUES,
		    SCATTER_FORWARD);
    VecScatterEnd(tolocalall, hostU, user[bi].nhostU, INSERT_VALUES,
		  SCATTER_FORWARD);
   
    VecScatterDestroy(&tolocalall);
    VecDestroy(&hostU);
   
  }

  for (bi=0; bi<block_number; bi++) {
    DMDALocalInfo	info = user[bi].info;

    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    
    PetscInt    lmx, lmy, lmz;
    /**************************************************************************************************************************/
    /* Interpolate the Velocity on the interface nodes
       from the host nodes.
       itfc is the velocity at the cell corners
       hostU is the velocity at the cell centers of the host block */
    /**************************************************************************************************************************/
    DMDAVecGetArray(user[bi].fda, user[bi].lItfc, &itfc);
   
    for (i=0; i<user[bi].itfcptsnumber; i++) {
      ci = user[bi].itfcI[i];
      cj = user[bi].itfcJ[i];
      ck = user[bi].itfcK[i];

      hi = user[bi].itfchostI[i];
      hj = user[bi].itfchostJ[i];
      hk = user[bi].itfchostK[i];
      hb = user[bi].itfchostB[i];

      x = user[bi].itfchostx[i];
      y = user[bi].itfchosty[i];
      z = user[bi].itfchostz[i];

      if ((x+y+z<-1.) &&
	  ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze   ) {
	itfc[ck][cj][ci].x = 0.;
	itfc[ck][cj][ci].y = 0.;
	itfc[ck][cj][ci].z = 0.;

      } else
      if (ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze   ) {

	VecGetArray(user[hb].nhostU, &hostu);
	lmx = user[hb].info.mx; lmy = user[hb].info.my; lmz = user[hb].info.mz;
	itfc[ck][cj][ci].x = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  )) * 3]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1)) * 3]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  )) * 3]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1)) * 3]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  )) * 3]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1)) * 3]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  )) * 3]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1)) * 3]
			      * x     * y     * z); //i+1,j+1,k+1

	itfc[ck][cj][ci].y = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))*3+1]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))*3+1]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))*3+1]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))*3+1]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))*3+1]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))*3+1]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))*3+1]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))*3+1]
			      * x     * y     * z); //i+1,j+1,k+1

	itfc[ck][cj][ci].z = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))*3+2]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))*3+2]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))*3+2]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))*3+2]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))*3+2]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))*3+2]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))*3+2]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))*3+2]
			      * x     * y     * z); //i+1,j+1,k+1

	VecRestoreArray(user[hb].nhostU, &hostu);
      }
    } // for itfcnumber
  
    DMDAVecRestoreArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMLocalToLocalBegin(user[bi].fda, user[bi].lItfc, INSERT_VALUES,
			user[bi].lItfc);
    DMLocalToLocalEnd(user[bi].fda, user[bi].lItfc, INSERT_VALUES,
		      user[bi].lItfc);
  } // for bi
  //  PetscBarrier(PETSC_NULL);
  //  PetscPrintf(PETSC_COMM_WORLD, "Interface Interpolation DDDD\n");

  for (bi=0; bi<block_number; bi++) {

    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    PetscInt	lxs, lxe, lys, lye, lzs, lze;

    lxs = xs; lxe = xe;
    lys = ys; lye = ye;
    lzs = zs; lze = ze;
    
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
    
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

    Cmpnts ***ucont, ***kzet, ***jeta,***icsi;
   
    //    DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);
    DMDAVecGetArray(user[bi].fda, user[bi].Bcs.Ubcs, &ubcs);
    DMDAVecGetArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMDAVecGetArray(user[bi].fda, user[bi].Ucont, &ucont);
    DMDAVecGetArray(user[bi].fda, user[bi].lZet, &kzet);
    DMDAVecGetArray(user[bi].fda, user[bi].lEta, &jeta);
    DMDAVecGetArray(user[bi].fda, user[bi].lCsi, &icsi);

    /**************************************************************************************************************************/
    /* Create boundary condition for flux (ucont) and cell surface
       center velocity (ubcs) from the interpolated velocity (itfc).
       
       itfc is the velocity at the cell surface centers */
    /**************************************************************************************************************************/
    if (user[bi].bctype[0] == 0 && xs==0) {
      i=0;//1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
/* 	  ubcs[k][j][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j][i].z = itfc[k][j][i].z; */
	  ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				    itfc[k  ][j-1][i  ].x +
				    itfc[k  ][j  ][i  ].x +
				    itfc[k-1][j-1][i  ].x) ;
	  ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				    itfc[k  ][j-1][i  ].y +
				    itfc[k  ][j  ][i  ].y +
				    itfc[k-1][j-1][i  ].y);
	  ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				    itfc[k  ][j-1][i  ].z +
				    itfc[k  ][j  ][i  ].z +
				    itfc[k-1][j-1][i  ].z);
	  ucont[k][j][i].x = (ubcs[k][j][i].x*
			      icsi[k][j][i].x +
			      ubcs[k][j][i].y*
			      icsi[k][j][i].y +
			      ubcs[k][j][i].z*
			      icsi[k][j][i].z);

	}
      }
    }

    if (user[bi].bctype[1] == 0 && xe==mx) {
      i=mx-2;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
/* 	  ubcs[k][j][i+1].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j][i+1].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j][i+1].z = itfc[k][j][i].z; */
	  ubcs[k][j][i+1].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				      itfc[k  ][j-1][i  ].x +
				      itfc[k  ][j  ][i  ].x +
				      itfc[k-1][j-1][i  ].x);
	  ubcs[k][j][i+1].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				      itfc[k  ][j-1][i  ].y +
				      itfc[k  ][j  ][i  ].y +
				      itfc[k-1][j-1][i  ].y);
	  ubcs[k][j][i+1].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				      itfc[k  ][j-1][i  ].z +
				      itfc[k  ][j  ][i  ].z +
				      itfc[k-1][j-1][i  ].z);
	  ucont[k][j][i].x = (ubcs[k][j][i+1].x  *
			      icsi[k][j][i].x +
			      ubcs[k][j][i+1].y  *
			      icsi[k][j][i].y +
			      ubcs[k][j][i+1].z *
			      icsi[k][j][i].z);

	}
      }
    }

    if (user[bi].bctype[2] == 0 && ys==0) {
      j=0;//1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
/* 	  ubcs[k][j][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j][i].z = itfc[k][j][i].z; */
	  ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				    itfc[k  ][j  ][i-1].x +
				    itfc[k  ][j  ][i  ].x +
				    itfc[k-1][j  ][i-1].x);
	  ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				    itfc[k  ][j  ][i-1].y +
				    itfc[k  ][j  ][i  ].y +
				    itfc[k-1][j  ][i-1].y);
	  ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				    itfc[k  ][j  ][i-1].z +
				    itfc[k  ][j  ][i  ].z +
				    itfc[k-1][j  ][i-1].z);
	  ucont[k][j][i].y = ( ubcs[k][j][i].x*
			       jeta[k][j][i].x +
			       ubcs[k][j][i].y *
			       jeta[k][j][i].y +
			       ubcs[k][j][i].z *
			       jeta[k][j][i].z);
	}
      }
    }

    if (user[bi].bctype[3] == 0 && ye==my) {
      j=my-2;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
/* 	  ubcs[k][j+1][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j+1][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j+1][i].z = itfc[k][j][i].z; */
	  ubcs[k][j+1][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				      itfc[k  ][j  ][i-1].x +
				      itfc[k  ][j  ][i  ].x +
				      itfc[k-1][j  ][i-1].x);
	  ubcs[k][j+1][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				      itfc[k  ][j  ][i-1].y +
				      itfc[k  ][j  ][i  ].y +
				      itfc[k-1][j  ][i-1].y);
	  ubcs[k][j+1][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				      itfc[k  ][j  ][i-1].z +
				      itfc[k  ][j  ][i  ].z +
				      itfc[k-1][j  ][i-1].z);
	  ucont[k][j][i].y = ( ubcs[k][j+1][i].x*
			       jeta[k][j  ][i].x +
			       ubcs[k][j+1][i].y*
			       jeta[k][j  ][i].y +
			       ubcs[k][j+1][i].z*
			       jeta[k][j  ][i].z);
	 /*  if (i==1 && k==50)  PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d ubcs.x is %le ubcs.y is %le ubcs.z is %le \n",i,j,k,ubcs[k][j+1][i].x,ubcs[k][j+1][i].y,ubcs[k][j+1][i].z); */
/* 	  if (i==1 && k==100)  PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d ubcs.x is %le ubcs.y is %le ubcs.z is %le \n",i,j,k,ubcs[k][j+1][i].x,ubcs[k][j+1][i].y,ubcs[k][j+1][i].z); */
/* 	  if (i==1 && k==150)  PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d ubcs.x is %le ubcs.y is %le ubcs.z is %le \n",i,j,k,ubcs[k][j+1][i].x,ubcs[k][j+1][i].y,ubcs[k][j+1][i].z); */
/* 	  if (i==1 && k==200)  PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d ubcs.x is %le ubcs.y is %le ubcs.z is %le \n",i,j,k,ubcs[k][j+1][i].x,ubcs[k][j+1][i].y,ubcs[k][j+1][i].z); */
/* 	  if (i==1 && k==50)  PetscPrintf(PETSC_COMM_SELF, "Ucont.y @ i=%d j=%d k=%d is %le \n",i,j,k,ucont[k][j][i].y); */
/* 	  if (i==1 && k==100)  PetscPrintf(PETSC_COMM_SELF, "Ucont.y @ i=%d j=%d k=%d is %le \n",i,j,k,ucont[k][j][i].y); */
/* 	  if (i==1 && k==150)  PetscPrintf(PETSC_COMM_SELF, "Ucont.y @ i=%d j=%d k=%d is %le \n",i,j,k,ucont[k][j][i].y); */
/* 	  if (i==1 && k==200)  PetscPrintf(PETSC_COMM_SELF, "Ucont.y @ i=%d j=%d k=%d is %le \n",i,j,k,ucont[k][j][i].y); */
	}
      }
    }


    if (user[bi].bctype[4] == 0 && zs==0) {
      k=0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
/* 	  ubcs[k][j][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k][j][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k][j][i].z = itfc[k][j][i].z; */
	  ubcs[k][j][i].x = 0.25 * (itfc[k  ][j  ][i  ].x +
				    itfc[k  ][j  ][i-1].x +
				    itfc[k  ][j-1][i  ].x +
				    itfc[k  ][j-1][i-1].x);
	  ubcs[k][j][i].y = 0.25 * (itfc[k  ][j  ][i  ].y +
				    itfc[k  ][j  ][i-1].y +
				    itfc[k  ][j-1][i  ].y +
				    itfc[k  ][j-1][i-1].y);
	  ubcs[k][j][i].z = 0.25 * (itfc[k  ][j  ][i  ].z +
				    itfc[k  ][j  ][i-1].z +
				    itfc[k  ][j-1][i  ].z +
				    itfc[k  ][j-1][i-1].z);
	  ucont[k][j][i].z = ( ubcs[k][j][i].x*
			       kzet[k][j][i].x +
			       ubcs[k][j][i].y*
			       kzet[k][j][i].y +
			       ubcs[k][j][i].z *
			       kzet[k][j][i].z);
	  
	}
      }
    }

    if (user[bi].bctype[5] == 0 && ze==mz) {
      k=mz-2;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
/* 	  ubcs[k+1][j][i].x = itfc[k][j][i].x; */
/* 	  ubcs[k+1][j][i].y = itfc[k][j][i].y; */
/* 	  ubcs[k+1][j][i].z = itfc[k][j][i].z; */
	  ubcs[k+1][j][i].x = 0.25 * (itfc[k  ][j  ][i  ].x +
				      itfc[k  ][j  ][i-1].x +
				      itfc[k  ][j-1][i  ].x +
				      itfc[k  ][j-1][i-1].x);
	  ubcs[k+1][j][i].y = 0.25 * (itfc[k  ][j  ][i  ].y +
				      itfc[k  ][j  ][i-1].y +
				      itfc[k  ][j-1][i  ].y +
				      itfc[k  ][j-1][i-1].y);
	  ubcs[k+1][j][i].z = 0.25 * (itfc[k  ][j  ][i  ].z +
				      itfc[k  ][j  ][i-1].z +
				      itfc[k  ][j-1][i  ].z +
				      itfc[k  ][j-1][i-1].z);
	  ucont[k][j][i].z = ( ubcs[k+1][j][i].x*
			       kzet[k  ][j][i].x +
			       ubcs[k+1][j][i].y*
			       kzet[k  ][j][i].y +
			       ubcs[k+1][j][i].z *
			       kzet[k  ][j][i].z);

	}
      }
    }

    DMDAVecRestoreArray(user[bi].fda, user[bi].Ucont, &ucont);

    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

    // This part is for blanking
    if(blank && bi==0){
           
      DMDAVecGetArray(user[bi].da, user[bi].lNvert, &nvert);
      DMDAVecGetArray(user[bi].fda, user[bi].lUcont, &ucont);
      for (sb=1; sb<block_number; sb++) {
	PetscInt ip, im, jp, jm, kp, km;
	PetscInt ii, jj, kk;
	  
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    for (i=lxs; i<lxe; i++) {
	      ip = (i<mx-2?(i+1):(i));
	      im = (i>1   ?(i-1):(i));
	      
	      jp = (j<my-2?(j+1):(j));
	      jm = (j>1   ?(j-1):(j));
	      
	      kp = (k<mz-2?(k+1):(k));
	      km = (k>1   ?(k-1):(k));
	      
	      if (((int)(nvert[k][j][i]+0.5) < sb*10) &&
		  ((int)(nvert[k][j][i]+0.5) > sb*10-3) ) {
		// flux in x direction
		kk=k;	jj=j;
		for (ii=im; ii<ip; ii++) {
		  ucent.x = 0.25 * (itfc[kk-1][jj  ][ii  ].x +
				    itfc[kk  ][jj-1][ii  ].x +
				    itfc[kk  ][jj  ][ii  ].x +
				    itfc[kk-1][jj-1][ii  ].x) ;
		  ucent.y = 0.25 * (itfc[kk-1][jj  ][ii  ].y +
				    itfc[kk  ][jj-1][ii  ].y +
				    itfc[kk  ][jj  ][ii  ].y +
				    itfc[kk-1][jj-1][ii  ].y);
		  ucent.z = 0.25 * (itfc[kk-1][jj  ][ii  ].z +
				    itfc[kk  ][jj-1][ii  ].z +
				    itfc[kk  ][jj  ][ii  ].z +
				    itfc[kk-1][jj-1][ii  ].z);
		  ucont[kk][jj][ii].x = (ucent.x*
					 icsi[kk][jj][ii].x +
					 ucent.y*
					 icsi[kk][jj][ii].y +
					 ucent.z*
					 icsi[kk][jj][ii].z);
		}
		// flux in y direction
		kk=k;   ii=i;
		for (jj=jm; jj<jp; jj++) {
		  ucent.x = 0.25 * (itfc[kk-1][jj  ][ii  ].x +
				    itfc[kk  ][jj  ][ii-1].x +
				    itfc[kk  ][jj  ][ii  ].x +
				    itfc[kk-1][jj  ][ii-1].x);
		  ucent.y = 0.25 * (itfc[kk-1][jj  ][ii  ].y +
				    itfc[kk  ][jj  ][ii-1].y +
				    itfc[kk  ][jj  ][ii  ].y +
				    itfc[kk-1][jj  ][ii-1].y);
		  ucent.z = 0.25 * (itfc[kk-1][jj  ][ii  ].z +
				    itfc[kk  ][jj  ][ii-1].z +
				    itfc[kk  ][jj  ][ii  ].z +
				    itfc[kk-1][jj  ][ii-1].z);
		  ucont[kk][jj][ii].y = ( ucent.x*
					  jeta[kk][jj][ii].x +
					  ucent.y *
					  jeta[kk][jj][ii].y +
					  ucent.z *
					  jeta[kk][jj][ii].z);
		}
		// flux in z direction
		jj=j;  ii=i;
		for (kk=km; kk<kp; kk++) {
		  ucent.x = 0.25 * (itfc[kk  ][jj  ][ii  ].x +
				    itfc[kk  ][jj  ][ii-1].x +
				    itfc[kk  ][jj-1][ii  ].x +
				    itfc[kk  ][jj-1][ii-1].x);
		  ucent.y = 0.25 * (itfc[kk  ][jj  ][ii  ].y +
				    itfc[kk  ][jj  ][ii-1].y +
				    itfc[kk  ][jj-1][ii  ].y +
				    itfc[kk  ][jj-1][ii-1].y);
		  ucent.z = 0.25 * (itfc[kk  ][jj  ][ii  ].z +
				    itfc[kk  ][jj  ][ii-1].z +
				    itfc[kk  ][jj-1][ii  ].z +
				    itfc[kk  ][jj-1][ii-1].z);
		  ucont[kk][jj][ii].z = ( ucent.x*
					  kzet[kk][jj][ii].x +
					  ucent.y*
					  kzet[kk][jj][ii].y +
					  ucent.z *
					  kzet[kk][jj][ii].z);
		}
	      }// if (nvert)
	    }
	  }
	}
      } //for sb
      DMDAVecRestoreArray(user[bi].da, user[bi].lNvert, &nvert);
      DMDAVecRestoreArray(user[bi].fda, user[bi].lUcont, &ucont);
      //      PetscPrintf(PETSC_COMM_WORLD, "Local to global lUcont _U ");

      DMLocalToGlobalBegin(user[bi].fda, user[bi].lUcont,INSERT_VALUES,user[bi].Ucont);
      DMLocalToGlobalEnd(user[bi].fda, user[bi].lUcont,INSERT_VALUES,user[bi].Ucont);

      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

    } // if blank
   
    DMDAVecRestoreArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMDAVecRestoreArray(user[bi].fda, user[bi].Bcs.Ubcs, &ubcs);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lZet, &kzet);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lEta, &jeta);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lCsi, &icsi);
 
  } // bi

  for (bi=0; bi<block_number; bi++) {
    VecDestroy(&user[bi].nhostU);
    Contra2Cart(&(user[bi]));
    GhostNodeVelocity(&user[bi]);
  }

  return(0);      
}

PetscErrorCode Block_Interface_U_old(UserCtx *user) {
  PetscInt bi,sb;
  PetscInt ci, cj, ck;
  PetscInt hi, hj, hk, hb;
  PetscInt i, j, k;
  PetscReal x, y, z;
  Vec	 nhostU;
  Cmpnts ***itfc, ***ucat,***ubcs;
  PetscReal *hostu;

  // First calculate Phi components at grid nodes

  /* ucat is at cell centers while itfc is cell surface nodes! */
  for (bi=0; bi<block_number; bi++) {
    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    DMDAVecGetArray(user[bi].fda, user[bi].Itfc, &itfc);
    DMDAVecGetArray(user[bi].fda, user[bi].lUcat, &ucat);

    DMDACreateNaturalVector(user[bi].fda, &nhostU);
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  if (k<mz-1 && j<my-1 && i<mx-1) {
	    itfc[k][j][i].x = 0.125 * (ucat[k  ][j  ][i  ].x +
				       ucat[k  ][j  ][i+1].x +
				       ucat[k  ][j+1][i  ].x +
				       ucat[k  ][j+1][i+1].x +
				       ucat[k+1][j  ][i  ].x +
				       ucat[k+1][j  ][i+1].x +
				       ucat[k+1][j+1][i  ].x +
				       ucat[k+1][j+1][i+1].x);

	    itfc[k][j][i].y = 0.125 * (ucat[k  ][j  ][i  ].y +
				       ucat[k  ][j  ][i+1].y +
				       ucat[k  ][j+1][i  ].y +
				       ucat[k  ][j+1][i+1].y +
				       ucat[k+1][j  ][i  ].y +
				       ucat[k+1][j  ][i+1].y +
				       ucat[k+1][j+1][i  ].y +
				       ucat[k+1][j+1][i+1].y);
	    itfc[k][j][i].z = 0.125 * (ucat[k  ][j  ][i  ].z +
				       ucat[k  ][j  ][i+1].z +
				       ucat[k  ][j+1][i  ].z +
				       ucat[k  ][j+1][i+1].z +
				       ucat[k+1][j  ][i  ].z +
				       ucat[k+1][j  ][i+1].z +
				       ucat[k+1][j+1][i  ].z +
				       ucat[k+1][j+1][i+1].z);
	  }
	}
      }
    }
    DMDAVecRestoreArray(user[bi].fda, user[bi].Itfc, &itfc);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lUcat, &ucat);

    /* the velocities at cell nodes are saved in nhostU*/
    DMDAGlobalToNaturalBegin(user[bi].fda, user[bi].Itfc, INSERT_VALUES, nhostU);
    DMDAGlobalToNaturalEnd(user[bi].fda, user[bi].Itfc, INSERT_VALUES, nhostU);
    //    DALocalToGlobal(user[bi].fda, user[bi].lItfc, INSERT_VALUES, user[bi].Itfc);
    VecScatter tolocalall;
    VecScatterCreateToAll(nhostU, &tolocalall, &(user[bi].nhostU));
    VecScatterBegin(tolocalall, nhostU, user[bi].nhostU, INSERT_VALUES,
		    SCATTER_FORWARD);

    VecScatterEnd(tolocalall, nhostU, user[bi].nhostU, INSERT_VALUES,
		  SCATTER_FORWARD);

    VecScatterDestroy(&tolocalall);
    VecDestroy(&nhostU);
  }


  //  VecCreateSeq(PETSC_COMM_SELF, 24, &hostU);
  for (bi=0; bi<block_number; bi++) {
    DMDALocalInfo	info = user[bi].info;

    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    
    PetscInt    lmx, lmy, lmz;
    /**************************************************************************************************************************/
    /* Interpolate the Velocity on the interface nodes 
       from the host nodes.
       itfc is the velocity at the cell nodes 
       hostU is the velocity at the cell nodes of the host block */
    /**************************************************************************************************************************/
    DMDAVecGetArray(user[bi].fda, user[bi].lItfc, &itfc);
    for (i=0; i<user[bi].itfcptsnumber; i++) {
      ci = user[bi].itfcI[i];
      cj = user[bi].itfcJ[i];
      ck = user[bi].itfcK[i];

      hi = user[bi].itfchostI[i];
      hj = user[bi].itfchostJ[i];
      hk = user[bi].itfchostK[i];
      hb = user[bi].itfchostB[i];

      x = user[bi].itfchostx[i];
      y = user[bi].itfchosty[i];
      z = user[bi].itfchostz[i];

      if ((x+y+z<-1.) &&
	  ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze   ) {
	itfc[ck][cj][ci].x = 0.;
	itfc[ck][cj][ci].y = 0.;
	itfc[ck][cj][ci].z = 0.;
      } else
      if (ci>=xs && ci<xe &&
	  cj>=ys && cj<ye &&
	  ck>=zs && ck<ze   ) {

	VecGetArray(user[hb].nhostU, &hostu);
	lmx = user[hb].info.mx; lmy = user[hb].info.my; lmz = user[hb].info.mz;
	itfc[ck][cj][ci].x = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  )) * 3]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1)) * 3]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  )) * 3]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1)) * 3]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  )) * 3]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1)) * 3]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  )) * 3]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1)) * 3]
			      * x     * y     * z); //i+1,j+1,k+1

	itfc[ck][cj][ci].y = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))*3+1]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))*3+1]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))*3+1]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))*3+1]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))*3+1]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))*3+1]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))*3+1]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))*3+1]
			      * x     * y     * z); //i+1,j+1,k+1

	itfc[ck][cj][ci].z = (hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi  ))*3+2]
			      * (1-x) * (1-y) * (1-z) + //i,j,k
			      hostu[((hk  )*lmx*lmy + (hj  )*lmx + (hi+1))*3+2]
			      * x     * (1-y) * (1-z) + //i+1,j,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi  ))*3+2]
			      * (1-x) * y     * (1-z) + //i,j+1,k
 			      hostu[((hk  )*lmx*lmy + (hj+1)*lmx + (hi+1))*3+2]
			      * x     * y     * (1-z) + //i+1,j+1,k
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi  ))*3+2]
			      * (1-x) * (1-y) * z     + //i,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj  )*lmx + (hi+1))*3+2]
			      * x     * (1-y) * z     + //i+1,j,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi  ))*3+2]
			      * (1-x) * y     * z     + //i,j+1,k+1
 			      hostu[((hk+1)*lmx*lmy + (hj+1)*lmx + (hi+1))*3+2]
			      * x     * y     * z); //i+1,j+1,k+1

	VecRestoreArray(user[hb].nhostU, &hostu);
      }
      // Is the point a local point?
      // Get the host cell information from the host CPU
      // Update
    }
    PetscBarrier(PETSC_NULL);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMLocalToLocalBegin(user[bi].fda, user[bi].lItfc, INSERT_VALUES,
			user[bi].lItfc);
    DMLocalToLocalEnd(user[bi].fda, user[bi].lItfc, INSERT_VALUES,
		      user[bi].lItfc);
  }
  //  VecDestroy(&hostU);
  
  for (bi=0; bi<block_number; bi++) {

    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    PetscInt	lxs, lxe, lys, lye, lzs, lze;

    lxs = xs; lxe = xe;

    lys = ys; lye = ye;
    lzs = zs; lze = ze;
    
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
    
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;

    Cmpnts ***ucont, ***kzet, ***jeta,***icsi;

    //    DMDAVecGetArray(user[bi].fda, user[bi].Ucat, &ucat);
    DMDAVecGetArray(user[bi].fda, user[bi].Bcs.Ubcs, &ubcs);
    DMDAVecGetArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMDAVecGetArray(user[bi].fda, user[bi].Ucont, &ucont);
    DMDAVecGetArray(user[bi].fda, user[bi].lZet, &kzet);
    DMDAVecGetArray(user[bi].fda, user[bi].lEta, &jeta);
    DMDAVecGetArray(user[bi].fda, user[bi].lCsi, &icsi);

    /**************************************************************************************************************************/
    /* Create boundary condition for flux (ucont) and cell center
       velocity (ucat) from the interpolated velocity (itfc).
       itfc is the velocity at the cell nodes */
    /**************************************************************************************************************************/
    if (user[bi].bctype[0] == 0 && xs==0) {
      i=0;//1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				    itfc[k  ][j-1][i  ].x +
				    itfc[k  ][j  ][i  ].x +
				    itfc[k-1][j-1][i  ].x) ;
	  ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				    itfc[k  ][j-1][i  ].y +
				    itfc[k  ][j  ][i  ].y +
				    itfc[k-1][j-1][i  ].y);
	  ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				    itfc[k  ][j-1][i  ].z +
				    itfc[k  ][j  ][i  ].z +
				    itfc[k-1][j-1][i  ].z);

	  ucont[k][j][i].x = (ubcs[k][j][i].x*
			      icsi[k][j][i].x +
			      ubcs[k][j][i].y*
			      icsi[k][j][i].y +
			      ubcs[k][j][i].z*
			      icsi[k][j][i].z);

	}
      }
    }

    if (user[bi].bctype[1] == 0 && xe==mx) {
      i=mx-2;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i+1].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				      itfc[k  ][j-1][i  ].x +
				      itfc[k  ][j  ][i  ].x +
				      itfc[k-1][j-1][i  ].x);
	  ubcs[k][j][i+1].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				      itfc[k  ][j-1][i  ].y +
				      itfc[k  ][j  ][i  ].y +
				      itfc[k-1][j-1][i  ].y);
	  ubcs[k][j][i+1].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				      itfc[k  ][j-1][i  ].z +
				      itfc[k  ][j  ][i  ].z +
				      itfc[k-1][j-1][i  ].z);
	  ucont[k][j][i].x = (ubcs[k][j][i+1].x  *
			      icsi[k][j][i].x +
			      ubcs[k][j][i+1].y  *
			      icsi[k][j][i].y +
			      ubcs[k][j][i+1].z *
			      icsi[k][j][i].z);

	}
      }
    }

    if (user[bi].bctype[2] == 0 && ys==0) {
      j=0;//1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				    itfc[k  ][j  ][i-1].x +
				    itfc[k  ][j  ][i  ].x +
				    itfc[k-1][j  ][i-1].x);
	  ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				    itfc[k  ][j  ][i-1].y +
				    itfc[k  ][j  ][i  ].y +
				    itfc[k-1][j  ][i-1].y);
	  ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				    itfc[k  ][j  ][i-1].z +
				    itfc[k  ][j  ][i  ].z +
				    itfc[k-1][j  ][i-1].z);
	  ucont[k][j][i].y = ( ubcs[k][j][i].x*
			       jeta[k][j][i].x +
			       ubcs[k][j][i].y *
			       jeta[k][j][i].y +
			       ubcs[k][j][i].z *
			       jeta[k][j][i].z);
	}
      }
    }

    if (user[bi].bctype[3] == 0 && ye==my) {
      j=my-2;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j+1][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				      itfc[k  ][j  ][i-1].x +
				      itfc[k  ][j  ][i  ].x +
				      itfc[k-1][j  ][i-1].x);
	  ubcs[k][j+1][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				      itfc[k  ][j  ][i-1].y +
				      itfc[k  ][j  ][i  ].y +
				      itfc[k-1][j  ][i-1].y);
	  ubcs[k][j+1][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				      itfc[k  ][j  ][i-1].z +
				      itfc[k  ][j  ][i  ].z +
				      itfc[k-1][j  ][i-1].z);
	    
	  ucont[k][j][i].y = ( ubcs[k][j+1][i].x*
			       jeta[k][j  ][i].x +
			       ubcs[k][j+1][i].y*
			       jeta[k][j  ][i].y +
			       ubcs[k][j+1][i].z*
			       jeta[k][j  ][i].z);
	}
      }
    }


    if (user[bi].bctype[4] == 0 && zs==0) {
      k=0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.25 * (itfc[k  ][j  ][i  ].x +
				    itfc[k  ][j  ][i-1].x +
				    itfc[k  ][j-1][i  ].x +
				    itfc[k  ][j-1][i-1].x);
	  ubcs[k][j][i].y = 0.25 * (itfc[k  ][j  ][i  ].y +
				    itfc[k  ][j  ][i-1].y +
				    itfc[k  ][j-1][i  ].y +
				    itfc[k  ][j-1][i-1].y);
	  ubcs[k][j][i].z = 0.25 * (itfc[k  ][j  ][i  ].z +
				    itfc[k  ][j  ][i-1].z +
				    itfc[k  ][j-1][i  ].z +
				    itfc[k  ][j-1][i-1].z);
	  ucont[k][j][i].z = ( ubcs[k][j][i].x*
			       kzet[k][j][i].x +
			       ubcs[k][j][i].y*
			       kzet[k][j][i].y +
			       ubcs[k][j][i].z *
			       kzet[k][j][i].z);
	  
	}
      }
    }

    if (user[bi].bctype[5] == 0 && ze==mz) {
      k=mz-2;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k+1][j][i].x = 0.25 * (itfc[k  ][j  ][i  ].x +
				      itfc[k  ][j  ][i-1].x +
				      itfc[k  ][j-1][i  ].x +
				      itfc[k  ][j-1][i-1].x);
	  ubcs[k+1][j][i].y = 0.25 * (itfc[k  ][j  ][i  ].y +
				      itfc[k  ][j  ][i-1].y +
				      itfc[k  ][j-1][i  ].y +
				      itfc[k  ][j-1][i-1].y);
	  ubcs[k+1][j][i].z = 0.25 * (itfc[k  ][j  ][i  ].z +
				      itfc[k  ][j  ][i-1].z +
				      itfc[k  ][j-1][i  ].z +
				      itfc[k  ][j-1][i-1].z);
	  ucont[k][j][i].z = ( ubcs[k+1][j][i].x*
			       kzet[k  ][j][i].x +
			       ubcs[k+1][j][i].y*
			       kzet[k  ][j][i].y +
			       ubcs[k+1][j][i].z *
			       kzet[k  ][j][i].z);

	}
      }
    }

    // This part is for blanking
    if(blank && bi==0){
      PetscInt ip,jp,kp;
      PetscInt dispx, dispy, dispz, dispnn;

      for (sb=1; sb<block_number; sb++) { 
	ip=user[bi].ip[sb]; jp=user[bi].jp[sb]; kp=user[bi].kp[sb];
	dispx=user[bi].dispx[sb]; dispy=user[bi].dispy[sb]; dispz=user[bi].dispz[sb];
	dispnn=user[bi].dispnn[sb]; 
       
	//   i=ip-disp;
	for(i=ip-dispx; i<=ip-dispx+dispnn; i++) {
          if (i>=lxs && i<lxe){
	    for (j=jp-dispy+1; j<=jp+dispy; j++){
	      for (k=kp-dispz+1; k<=kp+dispz; k++){
                if (j>=lys && j<lye && k>=lzs && k<lze){
		  ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
					    itfc[k  ][j-1][i  ].x +
					    itfc[k  ][j  ][i  ].x +
					    itfc[k-1][j-1][i  ].x) ;
		  ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
					    itfc[k  ][j-1][i  ].y +
					    itfc[k  ][j  ][i  ].y +
					    itfc[k-1][j-1][i  ].y);
		  ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
					    itfc[k  ][j-1][i  ].z +
					    itfc[k  ][j  ][i  ].z +
					    itfc[k-1][j-1][i  ].z);
		  ucont[k][j][i].x = (ubcs[k][j][i].x*
				      icsi[k][j][i].x +
				      ubcs[k][j][i].y*
				      icsi[k][j][i].y +
				      ubcs[k][j][i].z*
				      icsi[k][j][i].z);
		  //	if (i==ip-dispx)  PetscPrintf(PETSC_COMM_SELF, "Blank Imin ijk %d %d %d  ucont %le!\n", i,j,k, ucont[k][j][i].x);
		  if (i>ip-dispx) {
		    ubcs[k][j][i].x = 0.25 * (itfc[k  ][j  ][i  ].x +
					      itfc[k  ][j  ][i-1].x +
					      itfc[k  ][j-1][i  ].x +
					      itfc[k  ][j-1][i-1].x);
		    ubcs[k][j][i].y = 0.25 * (itfc[k  ][j  ][i  ].y +
					      itfc[k  ][j  ][i-1].y +
					      itfc[k  ][j-1][i  ].y +
					      itfc[k  ][j-1][i-1].y);
		    ubcs[k][j][i].z = 0.25 * (itfc[k  ][j  ][i  ].z +
					      itfc[k  ][j  ][i-1].z +
					      itfc[k  ][j-1][i  ].z +
					      itfc[k  ][j-1][i-1].z); //1.5*(1-(ucat[k][j][i].y-0.5)*(ucat[k][j][i].y-0.5)/0.25); //itfc[k][j][i].z;
		    ucont[k][j][i].z = ( ubcs[k][j][i].x*
					 kzet[k][j][i].x +
					 ubcs[k][j][i].y*
					 kzet[k][j][i].y +
					 ubcs[k][j][i].z *
					 kzet[k][j][i].z);
		    ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
					      itfc[k  ][j  ][i-1].x +
					      itfc[k  ][j  ][i  ].x +
					      itfc[k-1][j  ][i-1].x);
		    ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
					      itfc[k  ][j  ][i-1].y +
					      itfc[k  ][j  ][i  ].y +
					      itfc[k-1][j  ][i-1].y);
		    ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
					      itfc[k  ][j  ][i-1].z +
					      itfc[k  ][j  ][i  ].z +
					      itfc[k-1][j  ][i-1].z);
		    ucont[k][j][i].y = ( ubcs[k][j][i].x*
					 jeta[k][j][i].x +
					 ubcs[k][j][i].y *
					 jeta[k][j][i].y +
					 ubcs[k][j][i].z *
					 jeta[k][j][i].z);
		  }
		}
	      }
	    }
	  }
	}
	
      
         // i=ip+disp;
	for(i=ip+dispx; i>=ip+dispx-dispnn; i--){
          if (i>=lxs && i<lxe){
	    for (j=jp-dispy+1; j<=jp+dispy; j++){
	      for (k=kp-dispz+1; k<=kp+dispz; k++){
		if (j>=lys && j<lye && k>=lzs && k<lze){
		  
/*           ubcs[k][j][i].x = itfc[k][j][i].x; */
/*           ubcs[k][j][i].y = itfc[k][j][i].y; */
/*           ubcs[k][j][i].z = itfc[k][j][i].z; */
		  ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
					    itfc[k  ][j-1][i  ].x +
					    itfc[k  ][j  ][i  ].x +
					    itfc[k-1][j-1][i  ].x) ;
		  ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
					    itfc[k  ][j-1][i  ].y +
					    itfc[k  ][j  ][i  ].y +
					    itfc[k-1][j-1][i  ].y);
		  ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
					    itfc[k  ][j-1][i  ].z +
					    itfc[k  ][j  ][i  ].z +
					    itfc[k-1][j-1][i  ].z);
		  ucont[k][j][i].x = (ubcs[k][j][i].x*
				      icsi[k][j][i].x +
				      ubcs[k][j][i].y*
				      icsi[k][j][i].y +
				      ubcs[k][j][i].z*
				      icsi[k][j][i].z);
		  //		  if (i==ip+dispx)  PetscPrintf(PETSC_COMM_SELF, "Blank Imax ijk %d %d %d  ucont %le!\n", i,j,k, ucont[k][j][i].x);
		  if (i>ip+dispx-dispnn) {
		    ubcs[k][j][i].x = 0.25 * (itfc[k  ][j  ][i  ].x +
					      itfc[k  ][j  ][i-1].x +
					      itfc[k  ][j-1][i  ].x +
					      itfc[k  ][j-1][i-1].x);
		    ubcs[k][j][i].y = 0.25 * (itfc[k  ][j  ][i  ].y +
					      itfc[k  ][j  ][i-1].y +
					      itfc[k  ][j-1][i  ].y +
					      itfc[k  ][j-1][i-1].y);
		    ubcs[k][j][i].z = 0.25 * (itfc[k  ][j  ][i  ].z +
					      itfc[k  ][j  ][i-1].z +
					      itfc[k  ][j-1][i  ].z +
					      itfc[k  ][j-1][i-1].z); //1.5*(1-(ucat[k][j][i].y-0.5)*(ucat[k][j][i].y-0.5)/0.25); //itfc[k][j][i].z;
		    ucont[k][j][i].z = ( ubcs[k][j][i].x*
					 kzet[k][j][i].x +
					 ubcs[k][j][i].y*
					 kzet[k][j][i].y +
					 ubcs[k][j][i].z *
					 kzet[k][j][i].z);
		    ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
					      itfc[k  ][j  ][i-1].x +
					      itfc[k  ][j  ][i  ].x +
					      itfc[k-1][j  ][i-1].x);
		    ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
					      itfc[k  ][j  ][i-1].y +
					      itfc[k  ][j  ][i  ].y +
					      itfc[k-1][j  ][i-1].y);
		    ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
					      itfc[k  ][j  ][i-1].z +
					      itfc[k  ][j  ][i  ].z +
					      itfc[k-1][j  ][i-1].z);
		    ucont[k][j][i].y = ( ubcs[k][j][i].x*
					 jeta[k][j][i].x +
					 ubcs[k][j][i].y *
					 jeta[k][j][i].y +
					 ubcs[k][j][i].z *
					 jeta[k][j][i].z);
		  }
		}
	      }

	    }
	  }
	}
      
	//      PetscPrintf(PETSC_COMM_WORLD, "Interface Interpolation DDDD\n");

	for (j=jp-dispy; j<=jp-dispy+dispnn; j++){
	  if (j>=lys && j<lye){
	    for (i=ip-dispx+1; i<=ip+dispx; i++){
	      //          for (i=lxs; i<lxe; i++){
	      for (k=kp-dispz+1; k<=kp+dispz; k++){
		if (i>=lxs && i<lxe && k>=lzs && k<lze){
		  /* 		ubcs[k][j][i].x = itfc[k][j][i].x; */
		  /* 		ubcs[k][j][i].y = itfc[k][j][i].y; */
		  /* 		ubcs[k][j][i].z = itfc[k][j][i].z; */
		  ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
					    itfc[k  ][j  ][i-1].x +
					    itfc[k  ][j  ][i  ].x +
					    itfc[k-1][j  ][i-1].x);
		  ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
					    itfc[k  ][j  ][i-1].y +
					    itfc[k  ][j  ][i  ].y +
					    itfc[k-1][j  ][i-1].y);
		  ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
					    itfc[k  ][j  ][i-1].z +
					    itfc[k  ][j  ][i  ].z +
					    itfc[k-1][j  ][i-1].z);
		  //	  ucont[k][j][i].y = ucont[k][j-1][i].y;
		  
		  ucont[k][j][i].y = ( ubcs[k][j][i].x*
				       jeta[k][j][i].x +
				       ubcs[k][j][i].y *
				       jeta[k][j][i].y +
				       ubcs[k][j][i].z *
				       jeta[k][j][i].z);
		  //		  if (j==jp-dispy)PetscPrintf(PETSC_COMM_SELF, "Blank Jmin ijk %d %d %d  ucont %le!\n", i,j,k, ucont[k][j][i].y);
		  if (j>jp-dispy) {
		    ubcs[k][j][i].x = 0.25 * (itfc[k  ][j  ][i  ].x +
					      itfc[k  ][j  ][i-1].x +
					      itfc[k  ][j-1][i  ].x +
					      itfc[k  ][j-1][i-1].x);
		    ubcs[k][j][i].y = 0.25 * (itfc[k  ][j  ][i  ].y +
					      itfc[k  ][j  ][i-1].y +
					      itfc[k  ][j-1][i  ].y +
					      itfc[k  ][j-1][i-1].y);
		    ubcs[k][j][i].z = 0.25 * (itfc[k  ][j  ][i  ].z +
					      itfc[k  ][j  ][i-1].z +
					      itfc[k  ][j-1][i  ].z +
					      itfc[k  ][j-1][i-1].z); //1.5*(1-(ucat[k][j][i].y-0.5)*(ucat[k][j][i].y-0.5)/0.25); //itfc[k][j][i].z;
		    ucont[k][j][i].z = ( ubcs[k][j][i].x*
					 kzet[k][j][i].x +
					 ubcs[k][j][i].y*
					 kzet[k][j][i].y +
					 ubcs[k][j][i].z *
					 kzet[k][j][i].z);
		    ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
					      itfc[k  ][j-1][i  ].x +
					      itfc[k  ][j  ][i  ].x +
					      itfc[k-1][j-1][i  ].x) ;
		    ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
					      itfc[k  ][j-1][i  ].y +
					      itfc[k  ][j  ][i  ].y +
					      itfc[k-1][j-1][i  ].y);
		    ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
					      itfc[k  ][j-1][i  ].z +
					      itfc[k  ][j  ][i  ].z +
					      itfc[k-1][j-1][i  ].z);
		    ucont[k][j][i].x = (ubcs[k][j][i].x*
					icsi[k][j][i].x +
					ubcs[k][j][i].y*
					icsi[k][j][i].y +
					ubcs[k][j][i].z*
					icsi[k][j][i].z);
		  }
		}
	      }
	    }
	  }
	}
	
	//    PetscPrintf(PETSC_COMM_WORLD, "Interface Interpolation DDDD\n");
	for (j=jp+dispy; j>=jp+dispy-dispnn; j--){
	  if (j>=lys && j<lye){
	    for (i=ip-dispx+1; i<=ip+dispx; i++){
	      //          for (i=lxs; i<lxe; i++){
	      for (k=kp-dispz+1; k<=kp+dispz; k++){
		if (i>=lxs && i<lxe && k>=lzs && k<lze){
		  /* 		ubcs[k][j+1][i].x = itfc[k][j][i].x; */
		  /* 		ubcs[k][j+1][i].y = itfc[k][j][i].y; */
		  /* 		ubcs[k][j+1][i].z = itfc[k][j][i].z; */
		  ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
					    itfc[k  ][j  ][i-1].x +
					    itfc[k  ][j  ][i  ].x +
					    itfc[k-1][j  ][i-1].x);
		  ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
					    itfc[k  ][j  ][i-1].y +
					    itfc[k  ][j  ][i  ].y +
					    itfc[k-1][j  ][i-1].y);
		  ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
					    itfc[k  ][j  ][i-1].z +
					    itfc[k  ][j  ][i  ].z +
					    itfc[k-1][j  ][i-1].z);
		  //          ucont[k][j][i].y = ucont[k][j+1][i].y;
		  //ucont[k][j][i].y = ucat[k-12][j-12][i].y;
		  ucont[k][j][i].y = ( ubcs[k][j][i].x*
				       jeta[k][j][i].x +
				       ubcs[k][j][i].y *
				       jeta[k][j][i].y +
				       ubcs[k][j][i].z *
				       jeta[k][j][i].z);
		  
		  //if (j==jp+dispy )PetscPrintf(PETSC_COMM_SELF, "Blank Jmax ijk %d %d %d  ucont %le!\n", i,j,k, ucont[k][j][i].y);
		  if (j>jp+dispy-dispnn) {
		    ubcs[k][j][i].x = 0.25 * (itfc[k  ][j  ][i  ].x +
					      itfc[k  ][j  ][i-1].x +
					      itfc[k  ][j-1][i  ].x +
					      itfc[k  ][j-1][i-1].x);
		    ubcs[k][j][i].y = 0.25 * (itfc[k  ][j  ][i  ].y +
					      itfc[k  ][j  ][i-1].y +
					      itfc[k  ][j-1][i  ].y +
					      itfc[k  ][j-1][i-1].y);
		    ubcs[k][j][i].z = 0.25 * (itfc[k  ][j  ][i  ].z +
					      itfc[k  ][j  ][i-1].z +
					      itfc[k  ][j-1][i  ].z +
					      itfc[k  ][j-1][i-1].z); //1.5*(1-(ucat[k][j][i].y-0.5)*(ucat[k][j][i].y-0.5)/0.25); //itfc[k][j][i].z;
		    ucont[k][j][i].z = ( ubcs[k][j][i].x*
					 kzet[k][j][i].x +
					 ubcs[k][j][i].y*
					 kzet[k][j][i].y +
					 ubcs[k][j][i].z *
					 kzet[k][j][i].z);
		    ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
					      itfc[k  ][j-1][i  ].x +
					      itfc[k  ][j  ][i  ].x +
					      itfc[k-1][j-1][i  ].x) ;
		    ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
					      itfc[k  ][j-1][i  ].y +
					      itfc[k  ][j  ][i  ].y +
					      itfc[k-1][j-1][i  ].y);
		    ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
					      itfc[k  ][j-1][i  ].z +
					      itfc[k  ][j  ][i  ].z +
					      itfc[k-1][j-1][i  ].z);
		    ucont[k][j][i].x = (ubcs[k][j][i].x*
					icsi[k][j][i].x +
					ubcs[k][j][i].y*
					icsi[k][j][i].y +
					ubcs[k][j][i].z*
					icsi[k][j][i].z);
		  }
		  
		}
	      }
	    }
	  }
	}
	
	//PetscPrintf(PETSC_COMM_WORLD, "Interface Interpolation DDDD\n");
	for (k=kp-dispz; k<=kp-dispz+dispnn; k++){
	  if (k>=lzs && k<lze){
	    for (i=ip-dispx+1; i<=ip+dispx; i++){
	      //for (i=lxs; i<lxe; i++){
	      for (j=jp-dispy+1; j<=jp+dispy; j++){
		if (i>=lxs && i<lxe && j>=lys && j<lye){
		  ubcs[k][j][i].x = 0.25 * (itfc[k  ][j  ][i  ].x +
					    itfc[k  ][j  ][i-1].x +
					    itfc[k  ][j-1][i  ].x +
					    itfc[k  ][j-1][i-1].x);
		  ubcs[k][j][i].y = 0.25 * (itfc[k  ][j  ][i  ].y +
					    itfc[k  ][j  ][i-1].y +
					    itfc[k  ][j-1][i  ].y +
					    itfc[k  ][j-1][i-1].y);
		  ubcs[k][j][i].z = 0.25 * (itfc[k  ][j  ][i  ].z +
					    itfc[k  ][j  ][i-1].z +
					    itfc[k  ][j-1][i  ].z +
					    itfc[k  ][j-1][i-1].z); //1.5*(1-(ucat[k][j][i].y-0.5)*(ucat[k][j][i].y-0.5)/0.25); //itfc[k][j][i].z;
		  //	if(ti==11)PetscPrintf(PETSC_COMM_WORLD,"1. ucont[k][j][i].z = %le\n",ucont[k][j][i].z); 
		  
		  //	ucont[k][j][i].z = ucont[k-1][j][i].z; 
		  ucont[k][j][i].z = ( ubcs[k][j][i].x*
				       kzet[k][j][i].x +
				       ubcs[k][j][i].y*
				       kzet[k][j][i].y +
				       ubcs[k][j][i].z *
				       kzet[k][j][i].z);
		  
		  if (k>kp-dispz) {
		    ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
					      itfc[k  ][j  ][i-1].x +
					      itfc[k  ][j  ][i  ].x +
					      itfc[k-1][j  ][i-1].x);
		    ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
					      itfc[k  ][j  ][i-1].y +
					      itfc[k  ][j  ][i  ].y +
					      itfc[k-1][j  ][i-1].y);
		    ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
					      itfc[k  ][j  ][i-1].z +
					      itfc[k  ][j  ][i  ].z +
					      itfc[k-1][j  ][i-1].z);
		    ucont[k][j][i].y = ( ubcs[k][j][i].x*
					 jeta[k][j][i].x +
					 ubcs[k][j][i].y *
					 jeta[k][j][i].y +
					 ubcs[k][j][i].z *
					 jeta[k][j][i].z);
		    ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
					      itfc[k  ][j-1][i  ].x +
					      itfc[k  ][j  ][i  ].x +
					      itfc[k-1][j-1][i  ].x) ;
		    ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
					      itfc[k  ][j-1][i  ].y +
					      itfc[k  ][j  ][i  ].y +
					      itfc[k-1][j-1][i  ].y);
		    ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
					      itfc[k  ][j-1][i  ].z +
					      itfc[k  ][j  ][i  ].z +
					      itfc[k-1][j-1][i  ].z);
		    ucont[k][j][i].x = (ubcs[k][j][i].x*
					icsi[k][j][i].x +
					ubcs[k][j][i].y*
					icsi[k][j][i].y +
					ubcs[k][j][i].z*
					icsi[k][j][i].z);
		  }   
		  
		}
	      }
	    }
	  }
	}
	

	for (k=kp+dispz; k>=kp+dispz-dispnn; k--){
	  if (k>=lzs && k<lze){
	    for (i=ip-dispx+1; i<=ip+dispx; i++){
	      //          for (i=lxs; i<lxe; i++){
	      for (j=jp-dispy+1; j<=jp+dispy; j++){
		if (i>=lxs && i<lxe && j>=lys && j<lye){
		  ubcs[k][j][i].x = 0.25 * (itfc[k  ][j  ][i  ].x +
					    itfc[k  ][j  ][i-1].x +
					    itfc[k  ][j-1][i  ].x +
					    itfc[k  ][j-1][i-1].x);
		  ubcs[k][j][i].y = 0.25 * (itfc[k  ][j  ][i  ].y +
					    itfc[k  ][j  ][i-1].y +
					    itfc[k  ][j-1][i  ].y +
					    itfc[k  ][j-1][i-1].y);
		  ubcs[k][j][i].z = 0.25 * (itfc[k  ][j  ][i  ].z +
					    itfc[k  ][j  ][i-1].z +
					    itfc[k  ][j-1][i  ].z +
					    itfc[k  ][j-1][i-1].z);//1.5*(1-(ucat[k][j][i].y-0.5)*(ucat[k][j][i].y-0.5)/0.25); //itfc[k][j][i].z;
		  
		  //          ucont[k][j][i].z = ucont[k+1][j][i].z; 
		  //            ucont[k][j][i].z =ucat[k-12][j-12][i].z; 
		  ucont[k][j][i].z = ( ubcs[k][j][i].x*
				       kzet[k][j][i].x +
				       ubcs[k][j][i].y*
				       kzet[k][j][i].y +
				       ubcs[k][j][i].z *
				       kzet[k][j][i].z);
		  
		  if (k>kp+dispz-dispnn) {
		    ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
					      itfc[k  ][j  ][i-1].x +
					      itfc[k  ][j  ][i  ].x +
					      itfc[k-1][j  ][i-1].x);
		    ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
					      itfc[k  ][j  ][i-1].y +
					      itfc[k  ][j  ][i  ].y +
					      itfc[k-1][j  ][i-1].y);
		    ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
					      itfc[k  ][j  ][i-1].z +
					      itfc[k  ][j  ][i  ].z +
					      itfc[k-1][j  ][i-1].z);
		    ucont[k][j][i].y = ( ubcs[k][j][i].x*
					 jeta[k][j][i].x +
					 ubcs[k][j][i].y *
					 jeta[k][j][i].y +
					 ubcs[k][j][i].z *
					 jeta[k][j][i].z);
		    ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
					      itfc[k  ][j-1][i  ].x +
					      itfc[k  ][j  ][i  ].x +
					      itfc[k-1][j-1][i  ].x) ;
		    ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
					      itfc[k  ][j-1][i  ].y +
					      itfc[k  ][j  ][i  ].y +
					      itfc[k-1][j-1][i  ].y);
		    ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
					      itfc[k  ][j-1][i  ].z +
					      itfc[k  ][j  ][i  ].z +
					      itfc[k-1][j-1][i  ].z);
		    ucont[k][j][i].x = (ubcs[k][j][i].x*
					icsi[k][j][i].x +
					ubcs[k][j][i].y*
					icsi[k][j][i].y +
					ubcs[k][j][i].z*
					icsi[k][j][i].z);
		  }  
		  
		}
	      }
	    }
	  }
	}
	
      } //for sb
      //   PetscPrintf(PETSC_COMM_WORLD, "Interface Interpolation DDDD k\n");   
    } // if blank

    DMDAVecRestoreArray(user[bi].fda, user[bi].lItfc, &itfc);
  
    DMDAVecRestoreArray(user[bi].fda, user[bi].Bcs.Ubcs, &ubcs);
    DMDAVecRestoreArray(user[bi].fda, user[bi].Ucont, &ucont);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lZet, &kzet);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lEta, &jeta);
    
  }

  for (bi=0; bi<block_number; bi++) {
    VecDestroy(&user[bi].nhostU);
    Contra2Cart(&(user[bi]));
  }

  return 0;
    
}

      
PetscErrorCode GhostNodeVelocity(UserCtx *user)
{

  DM            da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	i, j, k;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts        ***ubcs, ***ucat,***lucat, ***csi, ***eta, ***zet,***cent;


  PetscReal Un, nx,ny,nz,A;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  PetscInt	gxs, gxe, gys, gye, gzs, gze;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);
 
  
/* ==================================================================================             */
/*   WALL FUNCTION */
/* ==================================================================================             */
  extern PetscInt wallfunction;
  if (wallfunction) {
    PetscReal ***ustar, ***aj, ***nvert;
  DMDAVecGetArray(fda, user->Ucat, &ucat);
  DMDAVecGetArray(da, user->lUstar, &ustar);
  DMDAVecGetArray(da, user->lAj,  &aj);
  DMDAVecGetArray(da, user->lNvert,  &nvert);

  // wall function for boundary
  //Dec 2015
  for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++) {
	if( nvert[k][j][i]<1.1 && ( (user->bctype[0]==-1 && i==1) || (user->bctype[1]==-1 &&  i==mx-2) ) ) {
	  double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
	  double sb, sc; 
	  double ni[3], nj[3], nk[3];
	  Cmpnts Uc, Ua;
	  
	  Ua.x = Ua.y = Ua.z = 0;
	  sb = 0.5/aj[k][j][i]/area;
	  
	  if(i==1) {
	    sc = 2*sb + 0.5/aj[k][j][i+1]/area;
	    Uc = ucat[k][j][i+1];
	  }
	  else {
	    sc = 2*sb + 0.5/aj[k][j][i-1]/area;
	    Uc = ucat[k][j][i-1];
	  }
	  
	  Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
	  if(i==mx-2) ni[0]*=-1, ni[1]*=-1, ni[2]*=-1;
	  
	  //if(i==1) printf("%f %f, %f %f %f, %e %e %e, %e\n", sc, sb, Ua.x, Ua.y, Ua.z, Uc.x, Uc.y, Uc.z, ucat[k][j][i+2].z);
	  wall_function (user, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], ni[0], ni[1], ni[2]);
	  nvert[k][j][i]=1;	/* set nvert to 1 to exclude from rhs */
	  //	  ucat[k][j][i] = lucat[k][j][i];	  
	}
	
	if( nvert[k][j][i]<1.1 && ( (user->bctype[2]==-1 && j==1) || (user->bctype[3]==-1 &&  j==my-2) ) ) {
	  double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
	  double sb, sc; 
	  double ni[3], nj[3], nk[3];
	  Cmpnts Uc, Ua;
	  
	  Ua.x = Ua.y = Ua.z = 0;
	  sb = 0.5/aj[k][j][i]/area;
	  
	  if(j==1) {
	    sc = 2*sb + 0.5/aj[k][j+1][i]/area;
	    Uc = ucat[k][j+1][i];
	  }
	  else {
	    sc = 2*sb + 0.5/aj[k][j-1][i]/area;
	    Uc = ucat[k][j-1][i];
	  }
	  
	  Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
	  if(j==my-2) nj[0]*=-1, nj[1]*=-1, nj[2]*=-1;
	  
	  wall_function (user, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], nj[0], nj[1], nj[2]);
	  nvert[k][j][i]=1;	/* set nvert to 1 to exclude from rhs */
	  //	  ucat[k][j][i] = lucat[k][j][i];
	}
	
      }
  DMDAVecRestoreArray(da, user->lAj,  &aj);
  DMDAVecRestoreArray(da, user->lNvert,  &nvert);
  DMDAVecRestoreArray(da, user->lUstar, &ustar);   
  DMDAVecRestoreArray(fda, user->Ucat, &ucat);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
 /*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

  } // if wallfunction

/* ==================================================================================             */
/*   Driven Cavity  */
/* ==================================================================================             */

  DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecGetArray(fda, user->Ucat,  &ucat);

  /* Designed for driven cavity problem (top(k=kmax) wall moving)
   u_x = 1 at k==kmax */
  if (user->bctype[5]==2) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;// - ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = 1.;//- ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = 0.;//- ucat[k-1][j][i].z;
	}
      }
    }
  }
  /*  ==================================================================================== */
  /*     Cylinder O-grid */
  /*  ==================================================================================== */
  if (user->bctype[3]==12) {
  /* Designed to test O-grid for flow over a cylinder at jmax velocity is 1 (horizontal) 
   u_x = 1 at k==kmax */
    if (ye==my) {
      j = ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 1.;
	}
      }
    }
  }
/*  ==================================================================================== */
/*     Annulus */
/*  ==================================================================================== */
/* designed to test periodic boundary condition for c grid j=0 rotates */
 
  if (user->bctype[2]==11) {
    DMDAVecGetArray(fda, user->Cent, &cent); 
   if (ys==0){
      j=0;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x =cent[k][j+1][i].y/sqrt(cent[k][j+1][i].x*cent[k][j+1][i].x+cent[k][j+1][i].y*cent[k][j+1][i].y);;
	  // ubcs[k][j][i].x=0.0;
	  ubcs[k][j][i].y =-cent[k][j+1][i].x/sqrt(cent[k][j+1][i].x*cent[k][j+1][i].x+cent[k][j+1][i].y*cent[k][j+1][i].y);
	  // ubcs[k][j][i].y = -cent[k][j+1][i].z/sqrt(cent[k][j+1][i].z*cent[k][j+1][i].z+cent[k][j+1][i].y*cent[k][j+1][i].y);
	  ubcs[k][j][i].z = 0.0;
	  // ubcs[k][j][i].z =cent[k][j+1][i].y/sqrt(cent[k][j+1][i].z*cent[k][j+1][i].z+cent[k][j+1][i].y*cent[k][j+1][i].y);
	}
      }
    }
   DMDAVecRestoreArray(fda, user->Cent,  &cent);
  }

/* ==================================================================================             */
/*   SYMMETRY BC */
/* ==================================================================================             */

  if (user->bctype[0]==3) {
    if (xs==0) {
    i= xs;

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	A=sqrt(csi[k][j][i].z*csi[k][j][i].z +
	       csi[k][j][i].y*csi[k][j][i].y +
	       csi[k][j][i].x*csi[k][j][i].x);
	nx=csi[k][j][i].x/A;
	ny=csi[k][j][i].y/A;
	nz=csi[k][j][i].z/A;
	Un=ucat[k][j][i+1].x*nx+ucat[k][j][i+1].y*ny+ucat[k][j][i+1].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i+1].x-Un*nx;//-V_frame.x;
	ubcs[k][j][i].y = ucat[k][j][i+1].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i+1].z-Un*nz;
      }
    }
    }
  }

  if (user->bctype[1]==3) {
    if (xe==mx) {
    i= xe-1;

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	A=sqrt(csi[k][j][i-1].z*csi[k][j][i-1].z +
	       csi[k][j][i-1].y*csi[k][j][i-1].y +
	       csi[k][j][i-1].x*csi[k][j][i-1].x);
	nx=csi[k][j][i-1].x/A;
	ny=csi[k][j][i-1].y/A;
	nz=csi[k][j][i-1].z/A;
	Un=ucat[k][j][i-1].x*nx+ucat[k][j][i-1].y*ny+ucat[k][j][i-1].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i-1].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j][i-1].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i-1].z-Un*nz;
      }
    }
    }
  }

  if (user->bctype[2]==3) {
    if (ys==0) {
    j= ys;

    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	A=sqrt(eta[k][j][i].z*eta[k][j][i].z +
	       eta[k][j][i].y*eta[k][j][i].y +
	       eta[k][j][i].x*eta[k][j][i].x);
	nx=eta[k][j][i].x/A;
	ny=eta[k][j][i].y/A;
	nz=eta[k][j][i].z/A;
	Un=ucat[k][j+1][i].x*nx+ucat[k][j+1][i].y*ny+ucat[k][j+1][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j+1][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j+1][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j+1][i].z-Un*nz;
      }
    }
    }
  }

  if (user->bctype[3]==3) {
    if (ye==my) {
    j=ye-1;

    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	A=sqrt(eta[k][j-1][i].z*eta[k][j-1][i].z +
	       eta[k][j-1][i].y*eta[k][j-1][i].y +
	       eta[k][j-1][i].x*eta[k][j-1][i].x);
	nx=eta[k][j-1][i].x/A;
	ny=eta[k][j-1][i].y/A;
	nz=eta[k][j-1][i].z/A;
	Un=ucat[k][j-1][i].x*nx+ucat[k][j-1][i].y*ny+ucat[k][j-1][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j-1][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j-1][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j-1][i].z-Un*nz;
      }
    }
    }
  }


  if (user->bctype[4]==3) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	A=sqrt(zet[k][j][i].z*zet[k][j][i].z +
	       zet[k][j][i].y*zet[k][j][i].y +
	       zet[k][j][i].x*zet[k][j][i].x);
	nx=zet[k][j][i].x/A;
	ny=zet[k][j][i].y/A;
	nz=zet[k][j][i].z/A;
	Un=ucat[k][j][i].x*nx+ucat[k][j][i].y*ny+ucat[k][j][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i].z-Un*nz;
	}
      }
    }
  }

  if (user->bctype[5]==3) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {  
	A=sqrt(zet[k-1][j][i].z*zet[k-1][j][i].z +
	       zet[k-1][j][i].y*zet[k-1][j][i].y +
	       zet[k-1][j][i].x*zet[k-1][j][i].x);
	nx=zet[k-1][j][i].x/A;
	ny=zet[k-1][j][i].y/A;
	nz=zet[k-1][j][i].z/A;
	Un=ucat[k-1][j][i].x*nx+ucat[k-1][j][i].y*ny+ucat[k-1][j][i].z*nz;
	ubcs[k][j][i].x = ucat[k-1][j][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k-1][j][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k-1][j][i].z-Un*nz;
	}
      }
    }
  }
/*  ==================================================================================== */
  /*     Rheology */
  /*  ==================================================================================== */
 
  // PetscPrintf(PETSC_COMM_WORLD, "moving plate velocity for rheology setup is %le \n",U_bc);

  if (user->bctype[2]==13){
    if (ys==0){
      j=0;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z =-U_bc;
	}
      }
    }
  }
  if (user->bctype[3]==13){
    if (ye==my){
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = U_bc;
	}
      }
    }
  }
  if (user->bctype[4]==13){
    if (zs==0){
      k=0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x =-U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 0.;
	}
      }
    }
  }
  if (user->bctype[5]==13){
    if (ze==mz){
      k=ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 0.;
	}
      }
    }
  }
/* ==================================================================================             */
  // boundary conditions on ghost nodes
/* ==================================================================================             */
//Mohsen Aug 2012
  if (xs==0 && user->bctype[0]!=7) {
    i = xs;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i+1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i+1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i+1].z;
      }
    }
  }

  if (xe==mx && user->bctype[0]!=7) {
    i = xe-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i-1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i-1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i-1].z;
      }
    }
  }

  if (ys==0 && user->bctype[2]!=7) {
    j = ys;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j+1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j+1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j+1][i].z;
      }
    }
  }

  if (ye==my && user->bctype[2]!=7) {
    j = ye-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j-1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j-1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j-1][i].z;
      }
    }
  }

  if (zs==0 && user->bctype[4]!=7) {
    k = zs;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k+1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k+1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k+1][j][i].z;
      }
    }
  }

  if (ze==mz && user->bctype[4]!=7) {
    k = ze-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k-1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k-1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k-1][j][i].z;
      }
    }
  }
  DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecRestoreArray(fda, user->Ucat,&ucat);
/* ==================================================================================             */
/*   Periodic BC Mohsen Aug 2012
/* ==================================================================================             */

  if (user->bctype[0]==7 || user->bctype[2]==7 ||  user->bctype[4]==7 ){

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
 
    DMDAVecGetArray(fda, user->lUcat, &lucat);
    DMDAVecGetArray(fda, user->Ucat, &ucat);
   
    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xs==0){
	i=xs;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    if(k>0 && k<user->KM && j>0 && j<user->JM){
	      ucat[k][j][i]=lucat[k][j][i-2];
	    }
	  }
	}
      }
    }
    if (user->bctype[2]==7 || user->bctype[3]==7){
      if (ys==0){
	j=ys;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    if(k>0 && k<user->KM && i>0 && i<user->IM){
	      ucat[k][j][i]=lucat[k][j-2][i];
	    }
	  }
	}
      }
    }
    if (user->bctype[4]==7 || user->bctype[5]==7){
      if (zs==0){
	k=zs;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if(j>0 && j<user->JM && i>0 && i<user->IM){
	      ucat[k][j][i]=lucat[k-2][j][i];
	    }
	  }
	}
      }
    }
    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xe==mx){
	i=mx-1;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    if(k>0 && k<user->KM && j>0 && j<user->JM){
	      ucat[k][j][i]=lucat[k][j][i+2];
	    }
	  }
	}
      }
    }
    if (user->bctype[2]==7 || user->bctype[3]==7){
      if (ye==my){
	j=my-1;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    if(k>0 && k<user->KM && i>0 && i<user->IM){
	    ucat[k][j][i]=lucat[k][j+2][i];
	    }
	  }
	}
      }
    }
    if (user->bctype[4]==7 || user->bctype[5]==7){
      if (ze==mz){
	k=mz-1;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if(j>0 && j<user->JM && i>0 && i<user->IM){
	      ucat[k][j][i].x=lucat[k+2][j][i].x;
	    }
	  }
	}
      }
    }
    
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
 /*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
  }

  // velocity on the corner point
  DMDAVecGetArray(fda, user->Ucat,  &ucat);

  if (zs==0 ) {
    k=0;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j][i+1].z);
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j][i-1].z);
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j+1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j+1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j+1][i].z);
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j-1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j-1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j-1][i].z);
      }
    }
  }

  if (ze==mz) {
    k=mz-1;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j][i+1].z);
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j][i-1].z);
      }
    }
    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j+1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j+1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j+1][i].z);
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j-1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j-1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j-1][i].z);
      }
    }
  }

  if (ys==0) {
    j=0;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j+1][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j+1][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j+1][i].z+ucat[k][j][i+1].z);
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j+1][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j+1][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j+1][i].z+ucat[k][j][i-1].z);
      }
    }
  }

  if (ye==my) {
    j=my-1;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j-1][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j-1][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j-1][i].z+ucat[k][j][i+1].z);
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j-1][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j-1][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j-1][i].z+ucat[k][j][i-1].z);
      }
    }
  }

  DMDAVecGetArray(fda, user->Ucat,  &ucat);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  if (user->bctype[0]==7 || user->bctype[2]==7 ||  user->bctype[4]==7 ){
    //i-direction
    DMDAVecGetArray(fda, user->lUcat, &lucat);
    DMDAVecGetArray(fda, user->Ucat, &ucat);
 
    if (user->bctype[0]==7){
      if (xs==0){
	i=xs;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    ucat[k][j][i]=lucat[k][j][i-2];
	  }
	}
      }
    }
    if (user->bctype[1]==7){
      if (xe==mx){
	i=xe-1;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    ucat[k][j][i]=lucat[k][j][i+2];
	  }
	}
      }
    }
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);

 /*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  //j-direction
  
    DMDAVecGetArray(fda, user->lUcat, &lucat);
    DMDAVecGetArray(fda, user->Ucat, &ucat);
 
    if (user->bctype[2]==7){
      if (ys==0){
	j=ys;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i]=lucat[k][j-2][i];
	  }
	}
      }
    }
  
    if (user->bctype[3]==7){
      if (ye==my){
	j=my-1;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	  ucat[k][j][i]=lucat[k][j+2][i];
	  }
	}
      }
    }

    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
  
 /*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  //k-direction
  
    DMDAVecGetArray(fda, user->lUcat, &lucat);
    DMDAVecGetArray(fda, user->Ucat, &ucat);
 
    if (user->bctype[4]==7){
      if (zs==0){
	k=zs;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	  ucat[k][j][i]=lucat[k-2][j][i];
	  }
	}
      }
    }
    if (user->bctype[5]==7){
      if (ze==mz){
	k=mz-1;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i].x=lucat[k+2][j][i].x;
	  }
	}
      }
    }
   
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);

 /*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  }
 
  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  return(0);
  
}

/* Boundary condition defination (array user->bctype[0-5]):
   0:	interpolation/interface
   -1:  wallfunction
   1:	solid wall (not moving)
   2:	moving solid wall (U=1)
   3:   slip wall/symmetry
   5:	Inlet
   4:	Outlet
   6:   farfield
   7:   periodic
   8:   Characteristic BC
   9:   Analytical Vortex
   10:  Oulet Junction
   11:  Annulus
   12:  Ogrid
   13:  Rheology
   14:  Outlet with Interface
  
*/

PetscErrorCode FormBCS(UserCtx *user)
{
  DM            da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	i, j, k;

  PetscReal	***nvert,***lnvert; //local working array

  PetscReal	***p,***lp;
  Cmpnts	***ucont, ***ubcs, ***ucat,***lucat, ***csi, ***eta, ***zet;
  Cmpnts	***cent,***centx,***centy,***centz;
  PetscScalar	FluxIn, FluxOut,ratio;
  PetscScalar   lArea, AreaSum;
 
  PetscScalar   FarFluxIn=0., FarFluxOut=0., FarFluxInSum, FarFluxOutSum;
  PetscScalar   FarAreaIn=0., FarAreaOut=0., FarAreaInSum, FarAreaOutSum;
  PetscScalar   FluxDiff, VelDiffIn, VelDiffOut;
  Cmpnts        V_frame;
 
  PetscReal Un, nx,ny,nz,A;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  PetscInt	gxs, gxe, gys, gye, gzs, gze;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx ) lxe = xe-1;
  if (ye==my ) lye = ye-1;
  if (ze==mz ) lze = ze-1;


  DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);

  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);

  PetscInt ttemp;
  for (ttemp=0; ttemp<3; ttemp++) {
    DMDAVecGetArray(da, user->Nvert, &nvert); 
    DMDAVecGetArray(fda, user->lUcat,  &ucat);
    DMDAVecGetArray(fda, user->Ucont, &ucont);
/* ==================================================================================             */
/*   FAR-FIELD BC */
/* ==================================================================================             */
 
   // reset FAR FLUXES
    FarFluxIn = 0.; FarFluxOut=0.;
    FarAreaIn = 0.; FarAreaOut=0.;
    
    V_frame.x=0.;
    V_frame.y=0.;
    V_frame.z=0.;
    
    PetscReal lFlux_abs=0.0,FluxSum_abs=0.0,ratio=0.0;

   if (user->bctype[0]==6) {
    if (xs == 0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i+1].x;
	  ubcs[k][j][i].y = ucat[k][j][i+1].y;
	  ubcs[k][j][i].z = ucat[k][j][i+1].z;
	  ucont[k][j][i].x = ubcs[k][j][i].x * csi[k][j][i].x;
	  FarFluxIn += ucont[k][j][i].x;
	  lFlux_abs += fabs(ucont[k][j][i].x);
	  FarAreaIn += csi[k][j][i].x;
	 
	}
      }
    }
  }
    
  if (user->bctype[1]==6) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i-1].x;
	  ubcs[k][j][i].y = ucat[k][j][i-1].y;
	  ubcs[k][j][i].z = ucat[k][j][i-1].z;
	  ucont[k][j][i-1].x = ubcs[k][j][i].x * csi[k][j][i-1].x;
	  FarFluxOut += ucont[k][j][i-1].x;
	  lFlux_abs  += fabs(ucont[k][j][i-1].x);
	  FarAreaOut += csi[k][j][i-1].x;
	}
      }
    }
  }

  if (user->bctype[2]==6) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k][j+1][i].x;
	  ubcs[k][j][i].y = ucat[k][j+1][i].y;
	  ubcs[k][j][i].z = ucat[k][j+1][i].z;
	  ucont[k][j][i].y = ubcs[k][j][i].y * eta[k][j][i].y;
	  FarFluxIn += ucont[k][j][i].y;
	  lFlux_abs += fabs(ucont[k][j][i].y);
	  FarAreaIn += eta[k][j][i].y;
	}
      }
    }
  }
  
  if (user->bctype[3]==6) {
    if (ye==my) {
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k][j-1][i].x;
	  ubcs[k][j][i].y = ucat[k][j-1][i].y;
	  ubcs[k][j][i].z = ucat[k][j-1][i].z;
	  ucont[k][j-1][i].y = ubcs[k][j][i].y * eta[k][j-1][i].y;
	  FarFluxOut += ucont[k][j-1][i].y;
	  lFlux_abs  += fabs(ucont[k][j-1][i].y);
	  FarAreaOut += eta[k][j-1][i].y;
	}
      }
    }
  }

  if (user->bctype[4]==6 || user->bctype[4]==10) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k+1][j][i].x;
	  ubcs[k][j][i].y = ucat[k+1][j][i].y;
	  ubcs[k][j][i].z = ucat[k+1][j][i].z;
	  ucont[k][j][i].z = ubcs[k][j][i].z * zet[k][j][i].z;
	  FarFluxIn += ucont[k][j][i].z;
	  lFlux_abs += fabs(ucont[k][j][i].z);
	  FarAreaIn += zet[k][j][i].z;
	}
      }
    }
  }

  if (user->bctype[5]==6 || user->bctype[5]==10) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = ucat[k-1][j][i].z;
	  ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	  FarFluxOut += ucont[k-1][j][i].z;
	  lFlux_abs  += fabs(ucont[k-1][j][i].z); 
	  FarAreaOut += zet[k-1][j][i].z;
	}
      }
    }
  }
  
  MPI_Allreduce(&FarFluxIn,&FarFluxInSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&FarFluxOut,&FarFluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&lFlux_abs,&FluxSum_abs,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&FarAreaIn,&FarAreaInSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&FarAreaOut,&FarAreaOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 
  if (user->bctype[5]==6 || user->bctype[3]==6 || user->bctype[1]==6) {
    FluxDiff = 0.5*(FarFluxInSum - FarFluxOutSum) ; //Mohsen 
    ratio=(FarFluxInSum - FarFluxOutSum)/FluxSum_abs;
    if (fabs(FluxSum_abs) <1.e-10) ratio = 0.;

    PetscPrintf(PETSC_COMM_WORLD, "/FluxSum_abs %le ratio %le \n", FluxSum_abs,ratio);

    VelDiffIn  = FluxDiff / FarAreaInSum ;
   
    if (fabs(FarAreaInSum) <1.e-6) VelDiffIn = 0.; 
    // VelDiffIn=0.;
 
    VelDiffOut  = FluxDiff / FarAreaOutSum ;
   

    if (fabs(FarAreaOutSum) <1.e-6) VelDiffOut = 0.;

    PetscPrintf(PETSC_COMM_WORLD, "Far Flux Diff %d %le %le %le %le %le %le %le\n", ti, FarFluxInSum, FarFluxOutSum, FluxDiff, FarAreaInSum, FarAreaOutSum, VelDiffIn, VelDiffOut);
    
  }
  
  if (user->bctype[5]==10) {
    FluxDiff = FluxInSum -( FarFluxOutSum -FarFluxInSum) ;
    VelDiffIn  = 1/3.*FluxDiff / (FarAreaInSum);// +  FarAreaOutSum);
    if (fabs(FarAreaInSum) <1.e-6) VelDiffIn = 0.;

    VelDiffOut  = 2./3.* FluxDiff / (FarAreaOutSum) ;//(FarAreaInSum +  FarAreaOutSum) ;
   
    if (fabs(FarAreaOutSum) <1.e-6) VelDiffOut = 0.;

    PetscPrintf(PETSC_COMM_WORLD, "Far Flux Diff %d %le %le %le %le %le %le %le\n", ti, FarFluxInSum, FarFluxOutSum, FluxDiff, FarAreaInSum, FarAreaOutSum, VelDiffIn, VelDiffOut);
           
  }
  
  
  // scale global mass conservation

  if (user->bctype[5]==6 || user->bctype[5]==10) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ucont[k-1][j][i].z  += ratio*fabs(ucont[k-1][j][i].z);
	  ubcs[k][j][i].z = ucont[k-1][j][i].z/zet[k-1][j][i].z;
	  // ubcs[k][j][i].z = ucat[k-1][j][i].z + VelDiffOut ;//+ V_frame.z;
	  //ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;


	}
      }
    }
  }

  if (user->bctype[3]==6) {
    if (ye==my) {
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ucont[k][j-1][i].y +=ratio*fabs(ucont[k][j-1][i].y);
	  ubcs[k][j][i].y = ucont[k][j-1][i].y /eta[k][j-1][i].y;
	  //  ubcs[k][j][i].y = ucat[k][j-1][i].y + VelDiffOut;// + V_frame.y;
	  // ucont[k][j-1][i].y = ubcs[k][j][i].y * eta[k][j-1][i].y;

	}
      }
    }
  }
    
  if (user->bctype[1]==6) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ucont[k][j][i-1].x +=ratio*fabs(ucont[k][j][i-1].x);
	  ubcs[k][j][i].x = ucont[k][j][i-1].x / csi[k][j][i-1].x ;
	  //  ubcs[k][j][i].x = ucat[k][j][i-1].x + VelDiffOut;// + V_frame.x;
	  // ucont[k][j][i-1].x = ubcs[k][j][i].x * csi[k][j][i-1].x;

	}
      }
    }
  }


  if (user->bctype[0]==6) {
    if (xs == 0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ucont[k][j][i].x  -=ratio*fabs(ucont[k][j][i].x);
	  ubcs[k][j][i].x = ucont[k][j][i].x / csi[k][j][i].x;
	  // ubcs[k][j][i].x = ucat[k][j][i+1].x - VelDiffIn;// + V_frame.x;
	  // ucont[k][j][i].x = ubcs[k][j][i].x * csi[k][j][i].x;

	}
      }
    }
  }
  

  if (user->bctype[2]==6) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ucont[k][j][i].y -=ratio*fabs(ucont[k][j][i].y);
	  ubcs[k][j][i].y = ucont[k][j][i].y / eta[k][j][i].y;
	  // ubcs[k][j][i].y = ucat[k][j+1][i].y - VelDiffIn;// + V_frame.y;
	   // ucont[k][j][i].y = ubcs[k][j][i].y * eta[k][j][i].y;

	}
      }
    }
  }
  
  
  if (user->bctype[4]==6 || user->bctype[5]==10) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ucont[k][j][i].z -=ratio*fabs(ucont[k][j][i].z);
	  ubcs[k][j][i].z =ucont[k][j][i].z / zet[k][j][i].z;
	  // ubcs[k][j][i].z = ucat[k+1][j][i].z - VelDiffIn;// + V_frame.z;
	  // ucont[k][j][i].z = ubcs[k][j][i].z * zet[k][j][i].z;

	}
      }
    }
  }
 
/* ==================================================================================             */
/*   SYMMETRY BC */
/* ==================================================================================             */

  if (user->bctype[0]==3) {
    if (xs==0) {
    i= xs;

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	A=sqrt(csi[k][j][i].z*csi[k][j][i].z +
	       csi[k][j][i].y*csi[k][j][i].y +
	       csi[k][j][i].x*csi[k][j][i].x);
	nx=csi[k][j][i].x/A;
	ny=csi[k][j][i].y/A;
	nz=csi[k][j][i].z/A;
	Un=ucat[k][j][i+1].x*nx+ucat[k][j][i+1].y*ny+ucat[k][j][i+1].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i+1].x-Un*nx;//-V_frame.x;
	ubcs[k][j][i].y = ucat[k][j][i+1].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i+1].z-Un*nz;
	ucont[k][j][i].x = 0.;
      }
    }
    }
  }

  if (user->bctype[1]==3) {
    if (xe==mx) {
    i= xe-1;

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	A=sqrt(csi[k][j][i-1].z*csi[k][j][i-1].z +
	       csi[k][j][i-1].y*csi[k][j][i-1].y +
	       csi[k][j][i-1].x*csi[k][j][i-1].x);
	nx=csi[k][j][i-1].x/A;
	ny=csi[k][j][i-1].y/A;
	nz=csi[k][j][i-1].z/A;
	Un=ucat[k][j][i-1].x*nx+ucat[k][j][i-1].y*ny+ucat[k][j][i-1].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i-1].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j][i-1].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i-1].z-Un*nz;
	ucont[k][j][i-1].x = 0.;
      }
    }
    }
  }

  if (user->bctype[2]==3) {
    if (ys==0) {
    j= ys;

    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	A=sqrt(eta[k][j][i].z*eta[k][j][i].z +
	       eta[k][j][i].y*eta[k][j][i].y +
	       eta[k][j][i].x*eta[k][j][i].x);
	nx=eta[k][j][i].x/A;
	ny=eta[k][j][i].y/A;
	nz=eta[k][j][i].z/A;
	Un=ucat[k][j+1][i].x*nx+ucat[k][j+1][i].y*ny+ucat[k][j+1][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j+1][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j+1][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j+1][i].z-Un*nz;
	ucont[k][j][i].y = 0.;
      }
    }
    }
  }

  if (user->bctype[3]==3) {
    if (ye==my) {
    j=ye-1;

    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	A=sqrt(eta[k][j-1][i].z*eta[k][j-1][i].z +
	       eta[k][j-1][i].y*eta[k][j-1][i].y +
	       eta[k][j-1][i].x*eta[k][j-1][i].x);
	nx=eta[k][j-1][i].x/A;
	ny=eta[k][j-1][i].y/A;
	nz=eta[k][j-1][i].z/A;
	Un=ucat[k][j-1][i].x*nx+ucat[k][j-1][i].y*ny+ucat[k][j-1][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j-1][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j-1][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j-1][i].z-Un*nz;
	ucont[k][j-1][i].y = 0.;
      }
    }
    }
  }
  

  if (user->bctype[4]==3) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	A=sqrt(zet[k][j][i].z*zet[k][j][i].z +
	       zet[k][j][i].y*zet[k][j][i].y +
	       zet[k][j][i].x*zet[k][j][i].x);
	nx=zet[k][j][i].x/A;
	ny=zet[k][j][i].y/A;
	nz=zet[k][j][i].z/A;
	Un=ucat[k+1][j][i].x*nx+ucat[k+1][j][i].y*ny+ucat[k+1][j][i].z*nz;
	ubcs[k][j][i].x = ucat[k+1][j][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k+1][j][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k+1][j][i].z-Un*nz;
	ucont[k][j][i].z = 0.;
	}
      }
    }
  }

  if (user->bctype[5]==3) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	A=sqrt(zet[k-1][j][i].z*zet[k-1][j][i].z +
	       zet[k-1][j][i].y*zet[k-1][j][i].y +
	       zet[k-1][j][i].x*zet[k-1][j][i].x);
	nx=zet[k-1][j][i].x/A;
	ny=zet[k-1][j][i].y/A;
	nz=zet[k-1][j][i].z/A;
	Un=ucat[k-1][j][i].x*nx+ucat[k-1][j][i].y*ny+ucat[k-1][j][i].z*nz;
	ubcs[k][j][i].x = ucat[k-1][j][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k-1][j][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k-1][j][i].z-Un*nz;
	ucont[k-1][j][i].z = 0.;
	}
      }
    }
  }
 
/* ==================================================================================             */
/*     CHARACTERISTIC OUTLET BC :8 */
/* ==================================================================================             */

  if (user->bctype[5]==8) {
    if (ze == mz) {
      k = ze-2;
      FluxOut = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  FluxOut += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxOut = 0.;
    }
    
    FluxIn = FluxInSum + FarFluxInSum;

    MPI_Allreduce(&FluxOut,&FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    // PetscGlobalSum(PETSC_COMM_WORLD,&FluxOut, &FluxOutSum);

    //ratio = FluxInSum / FluxOutSum;
    ratio = FluxIn / FluxOutSum;
    if (fabs(FluxOutSum) < 1.e-6) ratio = 1.;
    //if (fabs(FluxInSum) <1.e-6) ratio = 0.;
    if (fabs(FluxIn) <1.e-6) ratio = 0.;
    PetscPrintf(PETSC_COMM_WORLD, "Char Ratio %d %le %le %le %le %d %d\n", ti, ratio, FluxIn, FluxOutSum, FarFluxInSum,zs, ze);

    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  if (ti==0 || ti==1)
	    if (inletprofile<0)
	      ubcs[k][j][i].z = -1.;
	    else if (user->bctype[4]==6)
	      ubcs[k][j][i].z = 0.;
	    else
	      ubcs[k][j][i].z = 1.;//ubcs[0][j][i].z;//-1.;//1.;
	  
	  else
	    ucont[k-1][j][i].z = ucont[k-1][j][i].z*ratio;
	  
	  ubcs[k][j][i].z = ucont[k-1][j][i].z / zet[k-1][j][i].z;
	}
      }
    }
  }

  
/* ==================================================================================             */
/*     OUTLET BC :4 */
/* ==================================================================================             */

  
  if (user->bctype[5]==4 || user->bctype[5]==14 || user->bctype[5]==20) {
    lArea=0.;
    if (ze == mz) {
      //    k = ze-3;
      k=ze-1;
      FluxOut = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {

	  if ((nvert[k-1][j][i])<0.1) {
	  FluxOut += (ucat[k-1][j][i].x * (zet[k-1][j][i].x) +
		      ucat[k-1][j][i].y * (zet[k-1][j][i].y) +
		      ucat[k-1][j][i].z * (zet[k-1][j][i].z));

	  lArea += sqrt( (zet[k-1][j][i].x) * (zet[k-1][j][i].x) +
			 (zet[k-1][j][i].y) * (zet[k-1][j][i].y) +
			 (zet[k-1][j][i].z) * (zet[k-1][j][i].z));
	  }

	}
      }
    }
    else {
      FluxOut = 0.;
    }
    
    MPI_Allreduce(&FluxOut,&FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
   
    user->FluxOutSum = FluxOutSum;
    user->AreaOutSum = AreaSum;

    if (block_number>1 && user->bctype[5]==14) {
      FluxOutSum += user->FluxIntfcSum;
      //      AreaSum    += user->AreaIntfcSum;
    }

    FluxIn = FluxInSum + FarFluxInSum + user->FluxIntpSum;
    if (user->bctype[5]==20)
      ratio = (FluxIn / FluxOutSum);
    else
      ratio = (FluxIn - FluxOutSum) / AreaSum;
   

  /*   user->FluxOutSum += ratio*user->AreaOutSum; */
    user->FluxOutSum =0.0;
    FluxOut=0.0;
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if ((nvert[k-1][j][i])<0.1) {
	 
	    ubcs[k][j][i].x = ucat[k-1][j][i].x;//+ratio;
	    ubcs[k][j][i].y = ucat[k-1][j][i].y;
	    ubcs[k][j][i].z = ucat[k-1][j][i].z;// + ratio;//*n_z;

	    //  ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	    if (user->bctype[5]==20)
	      ucont[k-1][j][i].z = (ubcs[k][j][i].x * (zet[k-1][j][i].x ) +
				    ubcs[k][j][i].y * (zet[k-1][j][i].y ) +
				    ubcs[k][j][i].z * (zet[k-1][j][i].z ))*ratio;
	    
	    else{
	      ucont[k-1][j][i].z = (ubcs[k][j][i].x * (zet[k-1][j][i].x ) +
				    ubcs[k][j][i].y * (zet[k-1][j][i].y ) +
				    ubcs[k][j][i].z * (zet[k-1][j][i].z ))
		+ ratio * sqrt( (zet[k-1][j][i].x) * (zet[k-1][j][i].x) +
				(zet[k-1][j][i].y) * (zet[k-1][j][i].y) +
				(zet[k-1][j][i].z) * (zet[k-1][j][i].z)); 

	      FluxOut += ucont[k-1][j][i].z;
	    }
	  }//if
	}
      }
    }
    
    MPI_Allreduce(&FluxOut,&FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    user->FluxOutSum = FluxOutSum;
    // PetscPrintf(PETSC_COMM_WORLD, "Ratio %d %le %.15le %.15le %le %le %le\n", ti, ratio, FluxInSum, FluxOutSum, user->FluxIntfcSum,user->FluxIntpSum, AreaSum);

  } else if (user->bctype[5]==2) {
  /* Designed for driven cavity problem (top(k=kmax) wall moving)
   u_x = 1 at k==kmax */
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;// - ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = 1;//sin(2*3.14*ti*user->dt);//1.;//- ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = 0.;//- ucat[k-1][j][i].z;
	}
      }
    }
  }
  

  /*   OUTLET at k==0 */
  if (user->bctype[4]==4) {
    lArea=0.;
    if (zs == 0) {
      k = zs;
      //      k= zs + 1;
      FluxOut = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {


	  FluxOut += ucat[k+1][j][i].z * zet[k][j][i].z ;

	  lArea += zet[k][j][i].z;



	}
      }
    }
    else {
      FluxOut = 0.;
    }
    
    FluxIn = FluxInSum + FarFluxInSum;

    MPI_Allreduce(&FluxOut,&FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  

    ratio = (FluxInSum - FluxOutSum) / AreaSum;
    PetscPrintf(PETSC_COMM_WORLD, "Ratio b %d  %le %le %le %le %d %d\n", ti, ratio, FluxInSum, FluxOutSum, AreaSum,zs, ze);
    
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	 
	  ubcs[k][j][i].x = ucat[k+1][j][i].x;
	  ubcs[k][j][i].y = ucat[k+1][j][i].y;
	  ubcs[k][j][i].z = ucat[k+1][j][i].z;
	  ucont[k][j][i].z = (ubcs[k][j][i].z+ratio) * zet[k][j][i].z;
	}
      }
    }
  }
  /*  ==================================================================================== */
  /*     Cylinder O-grid */
  /*  ==================================================================================== */
  if (user->bctype[3]==12) {
  /* Designed to test O-grid for flow over a cylinder at jmax velocity is 1 (horizontal) 
   u_x = 1 at k==kmax */
    if (ye==my) {
      j = ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 1.;
	}
      }
    }
  }
  /*  ==================================================================================== */
  /*     Annulus */
  /*  ==================================================================================== */
  /* designed to test periodic boundary condition for O-grid j=0 rotates */
  DMDAVecGetArray(fda, user->Cent, &cent);
  if (user->bctype[2]==11) {
    if (ys==0){
      j=0;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	
	/*   ubcs[k][j][i].x=0.0; */
	 
/* 	  ubcs[k][j][i].y = -cent[k][j+1][i].z/sqrt(cent[k][j+1][i].z*cent[k][j+1][i].z+cent[k][j+1][i].y*cent[k][j+1][i].y); */
	 
/* 	  ubcs[k][j][i].z =cent[k][j+1][i].y/sqrt(cent[k][j+1][i].z*cent[k][j+1][i].z+cent[k][j+1][i].y*cent[k][j+1][i].y); */
	  ubcs[k][j][i].x = cent[k][j+1][i].y/sqrt(cent[k][j+1][i].x*cent[k][j+1][i].x+cent[k][j+1][i].y*cent[k][j+1][i].y);;
	  ubcs[k][j][i].y =-cent[k][j+1][i].x/sqrt(cent[k][j+1][i].x*cent[k][j+1][i].x+cent[k][j+1][i].y*cent[k][j+1][i].y);
	  ubcs[k][j][i].z =0.0;
	  //  if(k==1)  PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d ubcs.y is %le\n",i,j,k,ubcs[k][j][i].y);
	}
      }
    }
  }
  /*  ==================================================================================== */
  /*     Rheology */
  /*  ==================================================================================== */
 

  PetscOptionsGetReal(PETSC_NULL, "-U_bc", &(U_bc), PETSC_NULL);
  //  PetscPrintf(PETSC_COMM_WORLD, "moving plate velocity for rheology setup is %le \n",U_bc);

  if (user->bctype[2]==13){
    if (ys==0){
      j=0;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = -U_bc;
	}
      }
    }
  }
  if (user->bctype[3]==13){
    if (ye==my){
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = U_bc;
	}
      }
    }
  }
 if (user->bctype[4]==13){
    if (zs==0){
      k=0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x =-U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 0.;
	}
      }
    }
  }
  if (user->bctype[5]==13){
    if (ze==mz){
      k=ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 0.;
	}
      }
    }
  }
  DMDAVecRestoreArray(fda, user->Cent, &cent);
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(da, user->Nvert, &nvert);
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

 
 
  Contra2Cart(user); // it also does global to local for Ucat
 
/* ==================================================================================             */
/*   WALL FUNCTION */
/* ==================================================================================             */
  extern PetscInt wallfunction;
  if (wallfunction) {
  PetscReal ***ustar, ***aj;
  //Mohsen Dec 2015
  DMDAVecGetArray(fda, user->Ucat, &ucat);
  DMDAVecGetArray(da, user->lUstar, &ustar);
  DMDAVecGetArray(da, user->lAj,  &aj);
  DMDAVecGetArray(da, user->Nvert, &nvert);

  // wall function for boundary
  for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++) {
	if( nvert[k][j][i]<1.1 && ( (user->bctype[0]==-1 && i==1) ||
				    (user->bctype[1]==-1 && i==mx-2) ) ) {
	  double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
	  double sb, sc;
	  double ni[3], nj[3], nk[3];
	  Cmpnts Uc, Ua;
	  
	  Ua.x = Ua.y = Ua.z = 0;
	  sb = 0.5/aj[k][j][i]/area;
	  
	  if(i==1) {
	    sc = 2*sb + 0.5/aj[k][j][i+1]/area;
	    Uc = ucat[k][j][i+1];
	  }
	  else {
	    sc = 2*sb + 0.5/aj[k][j][i-1]/area;
	    Uc = ucat[k][j][i-1];
	  }
	  
	  Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
	  if(i==mx-2) ni[0]*=-1, ni[1]*=-1, ni[2]*=-1;
	  
	  wall_function (user, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], ni[0], ni[1], ni[2]);
	  nvert[k][j][i]=1;	/* set nvert to 1 to exclude from rhs */
	}
	
	if( nvert[k][j][i]<1.1 && ( (user->bctype[2]==-1 && j==1) ||
				    (user->bctype[3]==-1 && j==my-2) ) ) {
	  double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
	  double sb, sc;
	  double ni[3], nj[3], nk[3];
	  Cmpnts Uc, Ua;
	  
	  Ua.x = Ua.y = Ua.z = 0;
	  sb = 0.5/aj[k][j][i]/area;
	  
	  if(j==1) {
	    sc = 2*sb + 0.5/aj[k][j+1][i]/area;
	    Uc = ucat[k][j+1][i];
	  }
	  else {
	    sc = 2*sb + 0.5/aj[k][j-1][i]/area;
	    Uc = ucat[k][j-1][i];
	  }
	  
	  Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
	  if(j==my-2) nj[0]*=-1, nj[1]*=-1, nj[2]*=-1;
	  
	  wall_function (user, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], nj[0], nj[1], nj[2]);
	  nvert[k][j][i]=1;	/* set nvert to 1 to exclude from rhs */
	  
	}
      }

  DMDAVecRestoreArray(da, user->lAj,  &aj);
  DMDAVecRestoreArray(da, user->lUstar, &ustar);
  DMDAVecRestoreArray(fda, user->Ucat, &ucat);
  DMDAVecRestoreArray(da, user->Nvert, &nvert);
  DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
  DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
 /*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

 
  }

/* ==================================================================================             */

  DMDAVecGetArray(fda, user->Ucat, &ucat);
 
/* ==================================================================================             */
 
  // boundary conditions on ghost nodes
  if (xs==0 && user->bctype[0]!=7) {
    i = xs;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i+1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i+1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i+1].z;
      }
    }
  }

  if (xe==mx && user->bctype[0]!=7) {
    i = xe-1;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i-1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i-1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i-1].z;
      }
    }
  }


  if (ys==0 && user->bctype[2]!=7) {
    j = ys;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j+1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j+1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j+1][i].z;
      }
    }
  }

  if (ye==my && user->bctype[2]!=7) {
    j = ye-1;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j-1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j-1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j-1][i].z;
      }
    }
  }
 
  if (zs==0 && user->bctype[4]!=7) {
    k = zs;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k+1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k+1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k+1][j][i].z;
      }
    }
  }

  if (ze==mz && user->bctype[4]!=7) {
    k = ze-1;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k-1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k-1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k-1][j][i].z;
      }
    }
  }

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
 /* ==================================================================================             */
  /*   Periodic BC *///Mohsen
/* ==================================================================================             */
  if (user->bctype[0]==7 || user->bctype[2]==7 || user->bctype[4]==7){
    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
/* /\*   DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert); *\/ */
/* /\*   DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert); *\/ */
  //Mohsen Dec 2015
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lNvert, &lnvert);
    DMDAVecGetArray(fda, user->lUcat, &lucat);
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(da, user->Nvert, &nvert);
    DMDAVecGetArray(fda, user->Ucat, &ucat);
   
    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xs==0){
	i=xs;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    if(j>0 && k>0 && j<user->JM && k<user->KM){
	      ucat[k][j][i]=lucat[k][j][i-2];
	      p[k][j][i]=lp[k][j][i-2];
	      nvert[k][j][i]=lnvert[k][j][i-2];
	    }
	  }
	}
      }
    }
    if (user->bctype[2]==7 || user->bctype[3]==7){
      if (ys==0){
	j=ys;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    if(i>0 && k>0 && i<user->IM && k<user->KM){
	      ucat[k][j][i]=lucat[k][j-2][i];
	      p[k][j][i]=lp[k][j-2][i];
	      nvert[k][j][i]=lnvert[k][j-2][i];
	    }
	  }
	}
      }
    }
    if (user->bctype[4]==7 || user->bctype[5]==7){
      if (zs==0){
	k=zs;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if(i>0 && j>0 && i<user->IM && j<user->JM){
	      ucat[k][j][i]=lucat[k-2][j][i];
	      p[k][j][i]=lp[k-2][j][i];
	      nvert[k][j][i]=lnvert[k-2][j][i];
	    }
	  }
	}
      }
    }
    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xe==mx){
	i=mx-1;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    if(j>0 && k>0 && j<user->JM && k<user->KM){
	      ucat[k][j][i]=lucat[k][j][i+2];
	      p[k][j][i]=lp[k][j][i+2];
	      nvert[k][j][i]=lnvert[k][j][i+2];
	    }
	  }
	}
      }
    }
    if (user->bctype[2]==7 || user->bctype[3]==7){
      if (ye==my){
	j=my-1;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    if(i>0 && k>0 && i<user->IM && k<user->KM){
	      ucat[k][j][i]=lucat[k][j+2][i];
	      p[k][j][i]=lp[k][j+2][i];
	      nvert[k][j][i]=lnvert[k][j+2][i];
	    }
	  }
	}
      }
    }
    if (user->bctype[4]==7 || user->bctype[5]==7){
      if (ze==mz){
	k=mz-1;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if(i>0 && j>0 && i<user->IM && j<user->JM){
	      ucat[k][j][i]=lucat[k+2][j][i];
	      p[k][j][i]=lp[k+2][j][i];
	      nvert[k][j][i]=lnvert[k+2][j][i];
	    }
	  }
	}
      }
    }
 
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lNvert, &lnvert);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->Nvert, &nvert);

/*  /\*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); *\/ */
/* /\*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); *\/ */

/* /\*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); *\/ */
/* /\*   DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); *\/ */

/* /\*   DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert); *\/ */
/* /\*   DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert); *\/ */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
  }
 // 0 velocity on the corner point

  DMDAVecGetArray(fda, user->Ucat, &ucat);
  DMDAVecGetArray(da, user->P, &p);
  
  if (zs==0) {
    k=0;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j][i+1].z);
	p[k][j][i]= 0.5*(p[k+1][j][i]+p[k][j][i+1]);
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k+1][j][i]+p[k][j][i-1]);
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j+1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j+1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j+1][i].z);
	p[k][j][i] = 0.5*(p[k+1][j][i]+p[k][j+1][i]);
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j-1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j-1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j-1][i].z);
	p[k][j][i] = 0.5*(p[k+1][j][i]+p[k][j-1][i]);
      }
    }
  }
 
  if (ze==mz) {
    k=mz-1;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j][i+1].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j][i+1]);
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j][i-1]);
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j+1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j+1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j+1][i].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j+1][i]);
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j-1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j-1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j-1][i].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j-1][i]);
      }
    }
  }
 
  if (ys==0) {
    j=0;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j+1][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j+1][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j+1][i].z+ucat[k][j][i+1].z);
	p[k][j][i]= 0.5*(p[k][j+1][i]+p[k][j][i+1]);
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j+1][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j+1][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j+1][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k][j+1][i]+p[k][j][i-1]);
      }
    }
  }
 
  if (ye==my) {
    j=my-1;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j-1][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j-1][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j-1][i].z+ucat[k][j][i+1].z);
	p[k][j][i] = 0.5*(p[k][j-1][i]+p[k][j][i+1]);
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j-1][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j-1][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j-1][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k][j-1][i]+p[k][j][i-1]);
      }
    }
  }

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  DMDAVecRestoreArray(da, user->P, &p);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

  //Mohsen Nov 2012
  //Velocity and Presurre at corners for Periodic BC's

  if (user->bctype[0]==7 || user->bctype[2]==7 || user->bctype[4]==7){
  //i-direction

    DMDAVecGetArray(fda, user->lUcat,  &lucat);
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lNvert, &lnvert);
    DMDAVecGetArray(fda, user->Ucat,  &ucat);
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(da, user->Nvert, &nvert);
    
    if (user->bctype[0]==7){
      if (xs==0){
	i=xs;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    ucat[k][j][i]=lucat[k][j][i-2];
	    p[k][j][i]=lp[k][j][i-2];
	    nvert[k][j][i]=lnvert[k][j][i-2];
	  }
	}
      }
    }
    if (user->bctype[1]==7){
      if (xe==mx){
	i=xe-1;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    ucat[k][j][i].x=lucat[k][j][i+2].x;
	    p[k][j][i]=lp[k][j][i+2];
	    nvert[k][j][i]=lnvert[k][j][i+2];
	  }
	}
      }
    }
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lNvert, &lnvert);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->Nvert, &nvert);

 /*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

/*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*   DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */

/*   DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert); */
/*   DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert); */
  
    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
    
    //j-direction
    DMDAVecGetArray(fda, user->lUcat,  &lucat);
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lNvert, &lnvert);
    DMDAVecGetArray(fda, user->Ucat,  &ucat);
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(da, user->Nvert, &nvert);
    
    if (user->bctype[2]==7){
      if (ys==0){
	j=ys;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i]=lucat[k][j-2][i];
	    p[k][j][i]=lp[k][j-2][i];
	    nvert[k][j][i]=lnvert[k][j-2][i];
	  }
	}
      }
    }
    
    if (user->bctype[3]==7){
      if (ye==my){
	j=my-1;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i]=lucat[k][j+2][i];
	    p[k][j][i]=lp[k][j+2][i];
	    nvert[k][j][i]=lnvert[k][j+2][i];
	  }
	}
      }
    }
    
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lNvert, &lnvert);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->Nvert, &nvert);
    
/*   DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

/*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*   DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */

/*   DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert); */
/*   DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert); */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
    
    //k-direction
    DMDAVecGetArray(fda, user->lUcat,  &lucat);
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lNvert, &lnvert);
    DMDAVecGetArray(fda, user->Ucat,  &ucat);
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(da, user->Nvert, &nvert);
    
    if (user->bctype[4]==7){
      if (zs==0){
	k=zs;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i]=lucat[k-2][j][i];
	    p[k][j][i]=lp[k-2][j][i];
	    nvert[k][j][i]=lnvert[k-2][j][i];
	  }
	}
      }
    }
    if (user->bctype[5]==7){
      if (ze==mz){
	k=mz-1;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i]=lucat[k+2][j][i];
	    p[k][j][i]=lp[k+2][j][i];
	    nvert[k][j][i]=lnvert[k+2][j][i];
	  }
	}
      }
    }
    
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lNvert, &lnvert);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->Nvert, &nvert);
    
/*   DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

/*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*   DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */

/*   DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert); */
/*   DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert); */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
  }
  /* ==================================================================================             */
/*   Analytical Vortex BC */
/* ==================================================================================             */
 
  DMDAVecGetArray(fda, user->Ucat, &ucat);
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Cent, &cent); 
  DMDAVecGetArray(fda, user->Centx, &centx);
  DMDAVecGetArray(fda, user->Centy, &centy);
  DMDAVecGetArray(fda, user->Centz, &centz);

  if (user->bctype[0]==9) {
    if (xs==0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ucat[k][j][i].x=-cos(cent[k][j][i+1].x)*sin(cent[k][j][i+1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=-sin(cent[k][j][i+1].x)*cos(cent[k][j][i+1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;
	 
	  ucont[k][j][i].x =-(cos(centx[k][j][i].x)*sin(centx[k][j][i].y)*csi[k][j][i].x)*exp(-2.0*user->dt*(ti+1)/user->ren);
	}
      }
      if (ys==0) {
	j=ys;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j+1][i+1].x)*sin(cent[k][j+1][i+1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=-sin(cent[k][j+1][i+1].x)*cos(cent[k][j+1][i+1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
      if (ye==my) {
	j=ye-1;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j-1][i+1].x)*sin(cent[k][j-1][i+1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=-sin(cent[k][j-1][i+1].x)*cos(cent[k][j-1][i+1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
    }
  }
  if (user->bctype[1]==9) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ucat[k][j][i].x=-cos(cent[k][j][i-1].x)*sin(cent[k][j][i-1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=-sin(cent[k][j][i-1].x)*cos(cent[k][j][i-1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;
	
	  ucont[k][j][i-1].x =(-cos(centx[k][j][i-1].x)*sin(centx[k][j][i-1].y)*csi[k][j][i-1].x)*exp(-2.0*user->dt*(ti+1)/user->ren);
	}
      }
      if (ys==0) {
	j=ys;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j+1][i-1].x)*sin(cent[k][j+1][i-1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=-sin(cent[k][j+1][i-1].x)*cos(cent[k][j+1][i-1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
      if (ye==my) {
	j=ye-1;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j-1][i-1].x)*sin(cent[k][j-1][i-1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=-sin(cent[k][j-1][i-1].x)*cos(cent[k][j-1][i-1].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
    }
  }

  if (user->bctype[2]==9) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ucat[k][j][i].x=cos(cent[k][j+1][i].x)*sin(cent[k][j+1][i].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=sin(cent[k][j+1][i].x)*cos(cent[k][j+1][i].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;

	  ucont[k][j][i].y=(sin(centy[k][j][i].x)*cos(centy[k][j][i].y)*eta[k][j][i].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	}
      }
    }
  }
 
 
  if (user->bctype[3]==9) {
    if (ye==my) {
      j= ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ucat[k][j][i].x=cos(cent[k][j-1][i].x)*sin(cent[k][j-1][i].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].y=sin(cent[k][j-1][i].x)*cos(cent[k][j-1][i].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	  ucat[k][j][i].z =0.0;
	
	  ucont[k][j-1][i].y=(sin(centy[k][j-1][i].x)*cos(centy[k][j-1][i].y)*eta[k][j-1][i].y)*exp(-2.0*user->dt*(ti+1)/user->ren);
	}
      }
    }
  }
  if (user->bctype[4]==9) {
    if (zs==0) {
      k= zs;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  ucat[k][j][i].x=ucat[k+1][j][i].x;
	  ucat[k][j][i].y=ucat[k+1][j][i].y;
	  ucat[k][j][i].z=ucat[k+1][j][i].z;

	  ucont[k][j][i].z=0.0;
	}
      }
    }
  }
  if (user->bctype[5]==9) {
    if (ze==mz) {
      k= ze-1;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  ucat[k][j][i].x=ucat[k-1][j][i].x;
	  ucat[k][j][i].y=ucat[k-1][j][i].y;
	  ucat[k][j][i].z=ucat[k-1][j][i].z;

	  ucont[k-1][j][i].z=0.0;
	}
      }
    }
  }

  DMDAVecRestoreArray(fda, user->Cent, &cent);
  DMDAVecRestoreArray(fda, user->Centx, &centx);
  DMDAVecRestoreArray(fda, user->Centy, &centy);
  DMDAVecRestoreArray(fda, user->Centz, &centz);

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  DMDAVecRestoreArray(fda, user->Ucont,  &ucont);
 
  } // ttemp

 
  DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont); 

  return(0);
  }

PetscErrorCode fluxin(UserCtx *user)
{
 
  PetscInt  iRotate;

  PetscInt ts_p_cycle;
  PetscInt opening, closing;
  PetscInt open_steps, close_steps;

/*   ts_p_cycle = 500; */
/*   opening = 0; */
/*   open_steps = 50; */

/*   closing = 225; */
/*   close_steps = 50; */

  //ts_p_cycle = 1290;
  ts_p_cycle = 2500;//10000;//5000;

  opening = 10;
  open_steps = 100;

  closing = 580;
  close_steps = 80;

  PetscReal t_rel;

  iRotate = ti - ((ti / ts_p_cycle) * ts_p_cycle);

  // if (angle>.0 && iRotate>1058) iRotate-=angle;

  t_rel = iRotate * (1. / ts_p_cycle) * 860 + 6.8  ; //+7.15;
/*   t_rel = (iRotate-940) * (1. / ts_p_cycle) * 860/2 + */
/*                    940. * (1. / 2500.     ) * 860 + 6.8;     */

  PetscInt i;
  PetscBool interpolated = PETSC_FALSE;
  //PetscPrintf(PETSC_COMM_WORLD, "Inflow00 Rate %d %e %e %i\n",ti, Flux_in, t_rel, user->number_flowwave);
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    for (i=0; i<user->number_flowwave-1; i++) {
      /*       PetscPrintf(PETSC_COMM_WORLD, "Inflow Rate %e %e\n", Flux_in, user->inflow[i].t); */
      if (t_rel >= user->inflow[i].t && t_rel <= user->inflow[i+1].t) {
	Flux_in = user->inflow[i].f + (user->inflow[i+1].f - user->inflow[i].f) /
	  (user->inflow[i+1].t - user->inflow[i].t) *
	  (t_rel - user->inflow[i].t);
	/* 	PetscPrintf(PETSC_COMM_SELF, "Inflow Rate %i %e %e %e %e %e %e\n", i, Flux_in, t_rel, user->inflow[i].f, user->inflow[i].t, user->inflow[i+1].f, user->inflow[i+1].t); */
	interpolated = PETSC_TRUE;
      }
      if (interpolated) break;
    }  
    //if (t_rel > 350) Flux_in = 0.;
    
    //Flux_in=Flux_in;
    
    MPI_Bcast(&Flux_in, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  }
  else {
    MPI_Bcast(&Flux_in, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  }
  
  if (user->bctype[4]>1) {
    if (Flux_in<0.) {
      user->bctype[5]= 5;
      user->bctype[4]= 4;
      PetscPrintf(PETSC_COMM_WORLD, "IINNFLOW change inlet!!!!!!!%e\n", Flux_in);
    } else {
      user->bctype[5]= 4;
      user->bctype[4]= 5;
    }
  }

  PetscPrintf(PETSC_COMM_WORLD, "Angle %d %le %le flux-in %le intp%d\n",ti, t_rel, angle, Flux_in, interpolated);

  return 0;
}

PetscErrorCode OutflowVelocity(UserCtx *user, Vec Ucont)
{
  DM            fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal	lFluxOut = 0., ratio;
  Cmpnts ***ucont, ***zet, ***ucat, ***ubcs;
  
  PetscInt i, j, k;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  
  if (user->bctype[5] == 4) {

    //    OutflowFlux(user);

    Contra2Cart(user);
    DMDAVecGetArray(fda, Ucont, &ucont);
    DMDAVecGetArray(fda, user->lZet, &zet);
    DMDAVecGetArray(fda, user->Ucat, &ucat);
    DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
    /* Inflow flux at previous time step is 0, under this condition, it's assumed
       the flux difference at two time steps is uniformly distributed 
       to all outflow boundary nodes.*/
    if (ti==1) {//Flux_in_old < 1.e-6) {
      PetscReal lArea = 0., AreaSum;
      
      if (ze == mz) {
	k = ze-2;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    lArea += zet[k][j][i].z;
	  }
	}
      }
      MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      // PetscGlobalSum(PETSC_COMM_WORLD,&lArea, &AreaSum);

      PetscReal vd;
      vd = (user->FluxInSum) / AreaSum;
      PetscPrintf(PETSC_COMM_SELF, "FluxOOOO %e %e %e\n", vd, Flux_in, Flux_in);
      
      if (ze==mz) {
	k = ze-1;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    ubcs[k][j][i].z = vd;
	    ucont[k-1][j][i].z = vd * zet[k-1][j][i].z;
	  }
	}
      }
    }
    /* Scale the outflow flux to ensure global flux conservation */
    else {
      lFluxOut = 0.;
      if (ze==mz) {
	k = ze-2;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    lFluxOut += (ucat[k][j][i].x * (zet[k][j][i].x + zet[k-1][j][i].x) +
			 ucat[k][j][i].y * (zet[k][j][i].y + zet[k-1][j][i].y) +
			 ucat[k][j][i].z * (zet[k][j][i].z + zet[k-1][j][i].z))
	      * 0.5;
	  }
	}
      }

      MPI_Allreduce(&lFluxOut,&FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      //  PetscGlobalSum(PETSC_COMM_WORLD,&lFluxOut, &FluxOutSum);
      PetscBarrier(PETSC_NULL);
      ratio = user->FluxInSum / FluxOutSum;

      PetscPrintf(PETSC_COMM_WORLD, "Ratio %e %e\n", ratio, FluxOutSum);
      if (ze==mz) {
	k = ze-1;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    ubcs[k][j][i].x = ucat[k-1][j][i].x;
	    ubcs[k][j][i].y = ucat[k-1][j][i].y;
	    ubcs[k][j][i].z = ucat[k-1][j][i].z * ratio;
	    ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	  }
	}
      }
      
    }
    DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
    DMDAVecRestoreArray(fda, Ucont, &ucont);
    DMDAVecRestoreArray(fda, user->lZet, &zet);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);

/*     DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont); */
/*     DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont); */

/*     Contra2Cart(user, user->lUcont, user->Ucat); */
  }
  return 0;
}

PetscErrorCode SetInitialGuessToOne(UserCtx *user)
{
  DM da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  Vec           Coor;
  PetscReal	***nvert; 

  
  struct comp1{
    PetscScalar x, y, z;
  };


  struct comp1  ***ucont, ***zet,***coor, ***csi, ***eta, ***cent,***centx,***centy,***centz;
  
  PetscReal   gammadot=2.0;
  PetscOptionsGetReal(PETSC_NULL, "-gammadot", &gammadot, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-U_bc", &U_bc, PETSC_NULL);
  //  PetscPrintf(PETSC_COMM_WORLD, "shear rate for rheology setup is %le \n",gammadot);

  PetscReal xc,yc,r,uin,***p;
  PetscInt i, j, k;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

 
  //  DMDAGetGhostedCoordinates(da, &Coor);
  DMGetCoordinatesLocal(da, &Coor);
  DMDAVecGetArray(fda, Coor, &coor);

  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Centz, &centz);
  DMDAVecGetArray(fda, user->lZet,  &zet);
  DMDAVecGetArray(fda, user->lCsi,  &csi);  
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  
  extern PetscInt InitialGuessOne;

  for (k=zs ; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {	

	xc=centz[k][j][i].x-CMx_c;
	yc=centz[k][j][i].y-CMy_c;
	/*   xc = (coor[k][j][i  ].x + coor[k][j-1][i  ].x + */
/* 		coor[k][j][i-1].x + coor[k][j-1][i-1].x) * 0.25-CMx_c; */
/* 	  yc = (coor[k][j][i  ].y + coor[k][j-1][i  ].y + */
/* 		coor[k][j][i-1].y + coor[k][j-1][i-1].y) * 0.25-CMy_c; */
/* 	  zc = (coor[k][j][i  ].z + coor[k][j-1][i  ].z + */
/* 		coor[k][j][i-1].z + coor[k][j-1][i-1].z) * 0.25-CMz_c; */
	  //	  r = sqrt(xc * xc + yc * yc);
	r = sqrt(xc * xc + yc * yc);
	  
	  if (inletprofile == 4) { //fully-developed pipe flow
	    uin = 2.0*(1.0-4.0*r*r);
	  } else if (inletprofile == 7) {
	    uin = 1.5*(1.-4.*yc*yc);
	  } else if (inletprofile == 11) {
	    uin = 0.185;//0.2654;//
	  } else {
	    uin=1.;
	  }
	 
	  if (InitialGuessOne==2) {
	    ucont[k][j][i].x  = csi[k][j][i].x;
	    ucont[k][j][i].y  = eta[k][j][i].y;
	    ucont[k][j][i].z  = zet[k][j][i].z;
	  } else if (InitialGuessOne==3) {
	    // if (nvert[k][j][i]+nvert[k][j][i+1]<0.1) 
	    ucont[k][j][i].x  = uin*csi[k][j][i].z;
	    // if (nvert[k][j][i]+nvert[k][j+1][i]<0.1) 
	    ucont[k][j][i].y  = uin*eta[k][j][i].z;
	    // if (nvert[k][j][i]+nvert[k+1][j][i]<0.1)
	    ucont[k][j][i].z  = uin*zet[k][j][i].z;	    
	  } else if (InitialGuessOne==1) {
	    ucont[k][j][i].x  = 0.;
	    ucont[k][j][i].y  = 0.;
	    //if (nvert[k][j][i]+nvert[k+1][j][i]<0.1) 
	    ucont[k][j][i].z  	    
	      = uin* sqrt( (zet[k][j][i].x) * (zet[k][j][i].x) +
			   (zet[k][j][i].y) * (zet[k][j][i].y) +
			   (zet[k][j][i].z) * (zet[k][j][i].z));
	  }
      }
    }
  }
 
  if (InitialGuessOne==5) {
    for (k=zs ; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i]<9.5){
	    ucont[k][j][i].x  = csi[k][j][i].z;
	    ucont[k][j][i].y  = eta[k][j][i].z;
	    ucont[k][j][i].z  = zet[k][j][i].z;
	  }   
	}
      }
    }
  }
  
  if (InitialGuessOne==6) {
  
    for (k=zs ; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  
	  ucont[k][j][i].x  = 0.0;
	  ucont[k][j][i].y  = 0.0;
	  ucont[k][j][i].z  = (-U_bc+gammadot*centz[k][j][i].y)*zet[k][j][i].z;
	  if (j==my-2 && i==1 && k==1) PetscPrintf(PETSC_COMM_SELF, "w is  %le\n",(-U_bc+gammadot*centz[k][j][i].y));
	}
      }
    }
  
  }

  if (InitialGuessOne==4) {//Mohsen March 2012//
    
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(fda, user->Cent, &cent);
    DMDAVecGetArray(fda, user->Centx, &centx);
    DMDAVecGetArray(fda, user->Centy, &centy);
  /*   DMDAVecGetArray(fda, user->Centz, &centz); */
    
    for (k=lzs ; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=xs; i<lxe; i++) {
	  ucont[k][j][i].x=-(cos(centx[k][j][i].x)*sin(centx[k][j][i].y))*csi[k][j][i].x;
	}
      }
    }
    for (k=lzs ; k<lze; k++) {
      for (j=ys; j<lye; j++)  {
	for (i=lxs; i<lxe; i++) {
	  ucont[k][j][i].y=(sin(centy[k][j][i].x)*cos(centy[k][j][i].y))*eta[k][j][i].y;
	}
      }
    }
    for (k=zs ; k<lze; k++) {
      for (j=lys; j<lye; j++)  {
	for (i=lxs; i<lxe; i++) {
	  ucont[k][j][i].z=0.0;
	}
      }
    }
    
    for (k=lzs ; k<lze; k++) {
      for (j=lys; j<lye; j++)  {
	for (i=lxs; i<lxe; i++) {
	  p[k][j][i]=-0.25*(cos(2*cent[k][j][i].x)+cos(2*cent[k][j][i].y));
	}
      }
    }
    DMDAVecRestoreArray(fda, user->Centx, &centx);
    DMDAVecRestoreArray(fda, user->Centy, &centy);
  /*   DMDAVecRestoreArray(fda, user->Centz, &centz); */
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(fda, user->Cent, &cent);
    
  }//Mohsen

  DMDAVecRestoreArray(fda, Coor, &coor);
  DMDAVecRestoreArray(fda, user->Centz, &centz);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);
  DMDAVecRestoreArray(fda, user->lCsi,  &csi);  
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  // VecDestroy(&Coor);

 
  DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);

  DMGlobalToLocalBegin(user->da, user->P, INSERT_VALUES, user->lP);	
  DMGlobalToLocalEnd(user->da, user->P, INSERT_VALUES, user->lP);

  Contra2Cart(user);

  PetscReal normZ, normX, normY;
  VecStrideMax(user->Ucat, 0, PETSC_NULL, &normX);
  VecStrideMax(user->Ucat, 1, PETSC_NULL, &normY);
  VecStrideMax(user->Ucat, 2, PETSC_NULL, &normZ);  
  PetscPrintf(PETSC_COMM_WORLD, "Initial Eq 11111111111! %le %le %le\n",normX, normY, normZ);

  //Mohsen Sep 2012
  if ( InitialGuessOne==4){
 
    Divergence(user);
 
    Vec       Ucat_Exact,Ucat_Err;
    Cmpnts    ***ucat_exact;
    
    VecDuplicate(user->Ucat,&Ucat_Exact);
    VecDuplicate(user->Ucat,&Ucat_Err);
    VecSet(Ucat_Exact,0.0);
    VecSet(Ucat_Err,0.0);
    
    DMDAVecGetArray(fda, user->Cent, &cent);
    DMDAVecGetArray(fda,Ucat_Exact, &ucat_exact);
    
    for(k=lzs; k<lze ; k++){
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ucat_exact[k][j][i].x=-cos(cent[k][j][i].x)*sin(cent[k][j][i].y);
	  ucat_exact[k][j][i].y=sin(cent[k][j][i].x)*cos(cent[k][j][i].y);
	}
      }
    }
    
    DMDAVecRestoreArray(fda, user->Cent, &cent);
    DMDAVecRestoreArray(fda,Ucat_Exact, &ucat_exact);
    
    VecWAXPY(Ucat_Err,-1.0,user->Ucat,Ucat_Exact);
    
    PetscReal ucat_err,max_ucat_err,min_ucat_err,max_abs_ucat;
    PetscInt  mi,min_i,max_i;
    VecMax(Ucat_Err, &max_i, &max_ucat_err);
    VecMin(Ucat_Err, &min_i, &min_ucat_err);
    
    if (fabs(max_ucat_err)>fabs(min_ucat_err)) {
      max_abs_ucat=fabs(max_ucat_err);
      mi=max_i;
    } else {
      max_abs_ucat=fabs(min_ucat_err);
      mi=min_i;
    }
    
    VecNorm(Ucat_Err,NORM_INFINITY, & ucat_err);
    PetscPrintf(PETSC_COMM_WORLD, "Max initial ucat error is %le \n",ucat_err);
    VecDestroy(&Ucat_Exact);
    VecDestroy(&Ucat_Err);
  }
 
  return(0);
}

PetscErrorCode MomentumJet(UserCtx *user, PetscInt Kplane)
{
  DM            fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts ***ucont, ***ucat;   
  PetscInt i, j, k;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->lUcat , &ucat);
 
  k=Kplane;
  PetscReal lMomentum =0., Momentum=0.;
  if (k>zs && k<ze) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {	 
	lMomentum += fabs(ucont[k][j][i].z) * 0.5 *
	  (ucat[k][j][i].z + ucat[k+1][j][i].z);
      }
    }
  }

  MPI_Allreduce(&lMomentum,&Momentum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  //  PetscGlobalSum(PETSC_COMM_WORLD,&lMomentum, &Momentum);
  
  PetscReal lMomentum_in =0., Momentum_in=0.;
  if (zs==0) {
    k=0;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {	 
	lMomentum_in += fabs(ucont[k][j][i].z) * 0.5 *
	  (ucat[k][j][i].z + ucat[k+1][j][i].z);
      }
    }
  }

  MPI_Allreduce(&lMomentum_in,&Momentum_in,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  //  PetscGlobalSum(PETSC_COMM_WORLD,&lMomentum_in, &Momentum_in);

  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(fda, user->lUcat , &ucat);

  PetscPrintf(PETSC_COMM_WORLD, "Jet Momentum %d %d %le %le %le fluxin %le\n", ti, k, Momentum, Momentum_in, Flux_in*fabs(user->FluxInSum), Flux_in);

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {	
    char filen[80];
    sprintf(filen, "JetMomentum_%2.2d.dat", user->_this);
    FILE *f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d\t%.7e %le %le %le\n", ti, Momentum, Momentum_in, Flux_in*fabs(user->FluxInSum), Flux_in);
    fclose(f);
  }

  return(0);
}

PetscErrorCode fluxin_novostia(UserCtx *user)
{
  //iRotate = ti - ((ti / ts_p_cycle) * ts_p_cycle);
  PetscInt ts_p_cycle=8700;
  PetscInt i;
  PetscBool interpolated = PETSC_FALSE;
  PetscInt rank;
  FILE *fd, *fd1;

  if(!ti){
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (!rank) {
      fd = fopen("flow00", "r");
      for (i=0; i<ts_p_cycle; i++) 
	fscanf(fd, "%le\n",&flux_flow[i]);
      //
      MPI_Bcast(&flux_flow, 8700, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(&flux_flow, 8700, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    }
  }
  
  if (user->bctype[4]>1) {
    if (Flux_in<0.) {
      user->bctype[5]= 5;
      user->bctype[4]= 4;
      PetscPrintf(PETSC_COMM_WORLD, "IINNFLOW change inlet!!!!!!!%e\n", Flux_in);
    } else {
      user->bctype[5]= 4;
      user->bctype[4]= 5;
    }
  }

  Flux_in=flux_flow[ti];
  PetscPrintf(PETSC_COMM_WORLD, "ti: %d flux-in: %le\n",ti, Flux_in);
  
  return 0;
}
