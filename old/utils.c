#include <petscpf.h>
#include "variables.h"
#include <stdlib.h>
#include "time.h"
#include <math.h>
//#include "parallelComm.h"
#include "petsctime.h"
#include "petscdmcomposite.h"

extern void kaiser_wrap_(double a[3][3],int nrows,int n,double eigenv[3],double *trace,double *sume,int *ier);

extern void eigen_decomposition( double V[3][3], double d[3]);

extern PetscInt block_number, NumberOfBodies,blank,TwoD,rotateframe,tistart,NumberOfBlank;
extern PetscInt dIM,dJM,dKM;
extern PetscReal L_dim;
extern int nsend,nrecv,*sndMap,*snd_int,*rcvMap,**sndMapAll,**rcvMapAl,size;
extern PetscReal theta_x,theta_x_p;
//////////////////////////////
Mat Int_matrix;
DM  int_packer;
Vec U_int,Umult_int;
extern int cpu_size;

PetscInt lidxLocal1_matrix(PetscInt i, PetscInt j, PetscInt k, UserCtx *user,PetscInt blk);
PetscInt lidxLocal_matrix(PetscInt i, PetscInt j, PetscInt k,PetscInt blk,PetscInt cpu,PetscInt **xs,PetscInt **xm,PetscInt **ys,PetscInt **ym,PetscInt **zs,PetscInt **zm);

PetscInt find_cpu(PetscInt i,PetscInt j,PetscInt k,PetscInt hb,PetscInt size,PetscInt **xs,PetscInt **xm,PetscInt **ys,PetscInt **ym,PetscInt **zs,PetscInt **zm);



PetscInt *dg_start,**dl_start_blk;
//
PetscInt **cpu_xs;PetscInt **cpu_xm;
PetscInt **cpu_ys;PetscInt **cpu_ym;
PetscInt **cpu_zs;PetscInt **cpu_zm;
/////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// Quick test to determine the processor of a particle based on it's location //////

PetscBool  CPUPointIntersectCheck(UserCtx *user, Particle *particle){
  
  PetscInt i,j,k;
  Cmpnts loc;
  PetscInt xs,ys,zs,xe,ye,ze,mx,my,mz;
  PetscInt gxs,gys,gzs,gxe,gye,gze;
  PetscInt lxs,lys,lzs,lxe,lye,lze;
  DM da = user->da;
  DM fda = user->fda; 
  DMDALocalInfo info;
  Vec Coor;
  Cmpnts ***coor;
  PetscBool Intersects = PETSC_FALSE;
  
  DMGetCoordinatesLocal(user.da, &Coor); 
  DMDAVecGetArrayRead(fda,Coor,&coor);

  loc = particle->loc;

  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  lxs = xs-1; lxe = xe+1;
  lys = ys-1; lye = ye+1;
  lzs = zs-1; lze = ze+1;

  if (xs==0) lxs = xs;
  if (ys==0) lys = ys;
  if (zs==0) lzs = zs;

  if (xe==mx) lxe=xe;
  if (ye==my) lye=ye;
  if (ze==mz) lze=ze;
  
  for (k=gxs; k<gze; k++){
    for (j=gys; j<gye; j++){
      for (i=gxs; i<gxe; i++){
     
        if ((loc.x >= PetscRealPart(coor[k][j][i*3 + 0]) &&
             loc.x <= PetscRealPart(coor[k][j][(i+1)*3 + 0]) &&
             loc.y >= PetscRealPart(coor[k][j][i*3 + 1]) &&
             loc.y <= PetscRealPart(coor[k][j][(i+1)*3 + 1]) &&
	     loc.y >= PetscRealPart(coor[k][j][i*3 + 2]) &&
             loc.z <= PetscRealPart(coor[k][j][(i+1)*3 + 2]))) {
                    Intersects = PETSC_TRUE;
                    break;    
         } // if 
       } // i
       
       if(Intersects) break;
     } // j 
    
     if(Intersects) break;
   } // k

   DMDAVecRestoreArrayRead(fda,Coor,&coor)
   
   return Intersects;
} 

void findOBB(UserCtx *user,int bi,int size,int rank)
{
  int i,j,k,m,i3;
  double aa[3][3];
  double eigenv[3];
  double trace,sume;
  int ier;
  double xd[3];
  double xmin[3],xmax[3];
  int nrows,ncols;
  double xc[3],vec[3][3],dxc[3];

  Vec		coords;
  Cmpnts	***coor;
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info;
  PetscInt	xs, ys, zs, xe, ye, ze;
  PetscInt	mz, my, mx;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	gxs, gxe, gys, gye, gzs, gze;
  PetscInt      nnodes=0;

  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  lxs = xs-1; lxe = xe+1;
  lys = ys-1; lye = ye+1;
  lzs = zs-1; lze = ze+1;

  if (xs==0) lxs = xs;
  if (ys==0) lys = ys;
  if (zs==0) lzs = zs;

  if (xe==mx) lxe=xe;
  if (ye==my) lye=ye;
  if (ze==mz) lze=ze;

  

  //
  xc[0]=xc[1]=xc[2]=0;
  //
  // find centroid coordinates (not true centroid)
  //
  //DMDAGetGhostedCoordinates(da, &coords);
  DMDAVecGetArray(fda, user->cent_search, &coor);

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	if(i==lxs||i==(lxe-1)){
	xc[0]+=coor[k][j][i].x;
	xc[1]+=coor[k][j][i].y;
	xc[2]+=coor[k][j][i].z;
	nnodes++;
	//if (rank==0) PetscPrintf(PETSC_COMM_SELF,"bi:%d--- rank:%d---x(%d)=%le y(%d)=%le z(%d)=%le \n",bi,rank,i,coor[k][j][i].x,j,coor[k][j][i].y,k,coor[k][j][i].z);
	}
      }
    }
  }

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	if(j==lys||j==(lye-1)){
	xc[0]+=coor[k][j][i].x;
	xc[1]+=coor[k][j][i].y;
	xc[2]+=coor[k][j][i].z;
	nnodes++;
	}
      }
    }
  }

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	if(k==lzs||k==(lze-1)){
	xc[0]+=coor[k][j][i].x;
	xc[1]+=coor[k][j][i].y;
	xc[2]+=coor[k][j][i].z;
	nnodes++;
	}
      }
    }
  }

  //
  //PetscPrintf(PETSC_COMM_WORLD, "bi:%d nnodes:%d!\n",bi,nnodes);
  xc[0]/=nnodes;
  xc[1]/=nnodes;
  xc[2]/=nnodes;
  //PetscPrintf(PETSC_COMM_SELF, "bi:%d,rank:%d xs:%d ys:%d zs:%d lxe:%d lye:%d lze:%d %le %le %le!\n", bi,rank,xs,ys,zs,lxe,lye,lze, xc[0],xc[1],xc[2]);
  //PetscPrintf(PETSC_COMM_SELF, "bi:%d,lxe:%d lye:%d lze:%d %le %le %le!\n", bi,lze,zs,lze, xc[0],xc[1],xc[2]);
  //
  // find co-variance matrix
  // aa = [I11 I12 I13;I21 I22 I23;I31 I32 I33]
  //
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      aa[i][j]=0;
  //
  
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	if(i==lxs||i==(lxe-1)){
	aa[0][0]+=((coor[k][j][i].x-xc[0])*(coor[k][j][i].x-xc[0]));
	aa[1][1]+=((coor[k][j][i].y-xc[1])*(coor[k][j][i].y-xc[1]));
	aa[2][2]+=((coor[k][j][i].z-xc[2])*(coor[k][j][i].z-xc[2]));
	aa[1][0]+=((coor[k][j][i].x-xc[0])*(coor[k][j][i].y-xc[1]));
	aa[2][0]+=((coor[k][j][i].x-xc[0])*(coor[k][j][i].z-xc[2]));
	aa[2][1]+=((coor[k][j][i].y-xc[1])*(coor[k][j][i].z-xc[2]));
	}
      }
    }
  }

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	if(j==lys||j==(lye-1)){
	aa[0][0]+=((coor[k][j][i].x-xc[0])*(coor[k][j][i].x-xc[0]));
	aa[1][1]+=((coor[k][j][i].y-xc[1])*(coor[k][j][i].y-xc[1]));
	aa[2][2]+=((coor[k][j][i].z-xc[2])*(coor[k][j][i].z-xc[2]));
	aa[1][0]+=((coor[k][j][i].x-xc[0])*(coor[k][j][i].y-xc[1]));
	aa[2][0]+=((coor[k][j][i].x-xc[0])*(coor[k][j][i].z-xc[2]));
	aa[2][1]+=((coor[k][j][i].y-xc[1])*(coor[k][j][i].z-xc[2]));
	}
      }
    }
  }

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	if(k==lzs||k==(lze-1)){
	aa[0][0]+=((coor[k][j][i].x-xc[0])*(coor[k][j][i].x-xc[0]));
	aa[1][1]+=((coor[k][j][i].y-xc[1])*(coor[k][j][i].y-xc[1]));
	aa[2][2]+=((coor[k][j][i].z-xc[2])*(coor[k][j][i].z-xc[2]));
	aa[1][0]+=((coor[k][j][i].x-xc[0])*(coor[k][j][i].y-xc[1]));
	aa[2][0]+=((coor[k][j][i].x-xc[0])*(coor[k][j][i].z-xc[2]));
	aa[2][1]+=((coor[k][j][i].y-xc[1])*(coor[k][j][i].z-xc[2]));
	}
      }
    }
  }

  aa[0][1]=aa[1][0];
  aa[0][2]=aa[2][0];
  aa[1][2]=aa[2][1];
  
  // use kaisers method to estimate
  //eigen values and vectors of the covariance matrix
  
  nrows=3;
  ncols=3;

   eigen_decomposition(aa,eigenv);
  //kaiser_wrap_(aa,nrows,ncols,eigenv,&trace,&sume,&ier);
  //
  // copy the eigen vector basis on to vec
  //
  m=0;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      {
	vec[i][j]=aa[j][i];
      }
  //
  // find min and max bounds in the bounding box
  // vector basis
  //
  for(j=0;j<3;j++)
    {
      xmax[j]=-BIGVALUE;
      xmin[j]=BIGVALUE;
    }

  PetscInt kk;

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	if(i==lxs||i==(lxe-1)){

	for(kk=0;kk<3;kk++) xd[kk]=0;
	//
	for(kk=0;kk<3;kk++)
	  xd[kk]=(coor[k][j][i].x-xc[0])*vec[kk][0]+(coor[k][j][i].y-xc[1])*vec[kk][1]+(coor[k][j][i].z-xc[2])*vec[kk][2];
	//
	for(kk=0;kk<3;kk++)
	  {
	    xmax[kk]=max(xmax[kk],xd[kk]);
	    xmin[kk]=min(xmin[kk],xd[kk]);
	  }
	}
      }
    }
  }

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	if(j==lys||j==(lye-1)){

	for(kk=0;kk<3;kk++) xd[kk]=0;
	//
	for(kk=0;kk<3;kk++)
	  xd[kk]=(coor[k][j][i].x-xc[0])*vec[kk][0]+(coor[k][j][i].y-xc[1])*vec[kk][1]+(coor[k][j][i].z-xc[2])*vec[kk][2];
	//
	for(kk=0;kk<3;kk++)
	  {
	    xmax[kk]=max(xmax[kk],xd[kk]);
	    xmin[kk]=min(xmin[kk],xd[kk]);
	  }
	}
      }
    }
  }

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	if(k==lzs||k==(lze-1)){

	for(kk=0;kk<3;kk++) xd[kk]=0;
	//
	for(kk=0;kk<3;kk++)
	  xd[kk]=(coor[k][j][i].x-xc[0])*vec[kk][0]+(coor[k][j][i].y-xc[1])*vec[kk][1]+(coor[k][j][i].z-xc[2])*vec[kk][2];
	//
	for(kk=0;kk<3;kk++)
	  {
	    xmax[kk]=max(xmax[kk],xd[kk]);
	    xmin[kk]=min(xmin[kk],xd[kk]);
	  }
	}
      }
    }
  }
  //
  // find the extents of the box
  // and coordinates of the center w.r.t. xc
  // increase extents by 1% for tolerance
  //
  for(j=0;j<3;j++)
    {
      dxc[j]=1.01*(xmax[j]-xmin[j])*0.5;
      xd[j]=(xmax[j]+xmin[j])*0.5;
    }

  //
  // find the center of the box in
  // actual cartesian coordinates
  //
  
  for(j=0;j<3;j++)
    {
      for(k=0;k<3;k++)
	xc[j]+=(xd[k]*vec[k][j]);
    }
  
  //
  // PetscPrintf(PETSC_COMM_WORLD, "bi:%d, %le %le %le!\n", bi, xc[0],xc[1],xc[2]);

  ///////////////////////////////store in user
  for(k=0;k<3;k++)
    {
      user->lOBB_dxc[rank*3+k]=dxc[k];
      user->lOBB_xc[rank*3+k]=xc[k];
      for(j=0;j<3;j++){
	user->lOBB_ori_vec[rank*9+3*k+j]=vec[j][k];
      }
    }
  
  ///////////////////////////////
  
  DMDAVecRestoreArray(fda, user->cent_search, &coor);
  
}


int obbIntersectCheck(UserCtx *user, int bi, int sb, int cpu1, int rank1)
{
  int iflag;
  int i,j,k;
  int i1,i2,j1,j2;
  double r,r0,r1;
  double d1,d2;
  double eps=1e-12;
  double D[3];
  double c[3][3];
  double vA[3][3],xA[3],dxA[3],vB[3][3],xB[3],dxB[3];
  

  for(i=0;i<3;i++){
    // PetscPrintf(PETSC_COMM_SELF, "XA: %le\n",user[bi].OBB_xc[rank1*3+i]);
    xA[i]=user[bi].OBB_xc[rank1*3+i];
    xB[i]=user[sb].OBB_xc[cpu1*3+i];
    dxA[i]=user[bi].OBB_dxc[rank1*3+i];
    dxB[i]=user[sb].OBB_dxc[cpu1*3+i];
  }

  
  for(j=0;j<3;j++){
    for(i=0;i<3;i++){
      vA[i][j]=user[bi].OBB_ori_vec[rank1*9+3*j+i];
      vB[i][j]=user[sb].OBB_ori_vec[cpu1*9+3*j+i];
    }
  }
  

  //
  // D=distance between centers
  // C=scalar product of axes
  //
  for(i=0;i<3;i++) D[i]=xB[i]-xA[i];
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      {
	c[i][j]=0;
	for(k=0;k<3;k++)
	  c[i][j]=c[i][j]+vA[i][k]*vB[j][k];
      }
  //
  // separating axes based on the faces of box A
  //
  for(i=0;i<3;i++)
    {
      r0=dxA[i];
      r1=0;
      r=0;
      for(j=0;j<3;j++) 
	{
	  r1+=dxB[j]*fabs(c[i][j]);
	  r+=fabs(vA[i][j])*D[j];
	}
      if (r > (r0+r1+eps)) return 0;
    }
  //
  // separating axes based on the faces of box B
  //
  for(i=0;i<3;i++)
    {
      r1=dxB[i];
      r0=0;
      r=0;
      for(j=0;j<3;j++) 
	{
	  r0+=dxA[j]*fabs(c[j][i]);
	  r+=fabs(vB[i][j])*D[j];
	}
      if (r > (r0+r1+eps)) return 0;
    }
  //
  // cross products
  //
  for(i=0;i<3;i++)
    {
      i1=(i+1)%3;
      i2=(i+2)%3;
      for(j=0;j<3;j++)
	{
	  j1=(j+1)%3;
	  j2=(j+2)%3;
	  
	  r0=dxA[i1]*fabs(c[i2][j])+dxA[i2]*fabs(c[i1][j]);
	  r1=dxB[j1]*fabs(c[i][j2])+dxB[j2]*fabs(c[i][j1]);
	  
	  d2=0;
	  d1=0;
	  for(k=0;k<3;k++)
	    {
	      d2+=vA[i2][k]*D[k];
	      d1+=vA[i1][k]*D[k];
	    }
	  
	  r=fabs(c[i1][j]*d2-c[i2][j]*d1);
	  
	  if (r > (r0+r1+eps)) {
	    return 0;
	  }
	}
    }
  //
  // return zero if no separation can be found
  //
  return 1;
}
    

void writebbox(UserCtx *user,int bi,int rank,int ti)
{
  FILE *fp;
  char intstring[7];
  char fname[80];
  int l,i,k,j,m,il,ik,ij;
  double xx[3];
  double vec[3][3],dxc[3],xc[3];
  
    
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      vec[i][j]=user->lOBB_ori_vec[rank*9+3*j+i];

   for(i=0;i<3;i++){
     dxc[i]=user->lOBB_dxc[rank*3+i];
     xc[i]=user->lOBB_xc[rank*3+i];
   }
   
   
   sprintf(fname,"qbox%d_%d_%d.dat",bi,rank,ti);
   fp=fopen(fname,"w");
   fprintf(fp,"TITLE =\"Box file\"\n");
   fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\"\n");
   fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEPOINT\n",8,1);
   
   for(l=0;l<2;l++)
     {
       il=2*(l%2)-1;
       for(k=0;k<2;k++)
	 {
	   ik=2*(k%2)-1;
	   for(j=0;j<2;j++)
	     {
	       ij=2*(j%2)-1;
	       xx[0]=xx[1]=xx[2]=0;
	       for(m=0;m<3;m++)
		 xx[m]=xc[m]+ij*vec[0][m]*dxc[0]
		   +ik*vec[1][m]*dxc[1]
		   +il*vec[2][m]*dxc[2];	      
	      fprintf(fp,"%f %f %f\n",xx[0],xx[1],xx[2]);
	     }
	 }
     }
   fprintf(fp,"1 2 4 3 5 6 8 7\n");
   fclose(fp);
}



int cpuOBBPointIntersectCheck(UserCtx *user,double x[3],int rank)
{
  int i,j,k,l,m,n,p,i3;
  double xd[3];
  int cell_count=0;
  
  for(j=0;j<3;j++) xd[j]=0;
  //
  for(j=0;j<3;j++)
    xd[j]=(x[0]-user->OBB_xc[rank*3+0])*user->OBB_ori_vec[rank*9+3*0+j]+(x[1]-user->OBB_xc[rank*3+1])*user->OBB_ori_vec[rank*9+3*1+j]+(x[2]-user->OBB_xc[rank*3+2])*user->OBB_ori_vec[rank*9+3*2+j];
  
  //PetscPrintf(PETSC_COMM_SELF, "xd[0]=%le,xd[1]=%le,xd[2]=%le,!\n",xd[0],xd[1],xd[2]);
  //PetscPrintf(PETSC_COMM_SELF, "xd[0]=%le,xd[1]=%le,xd[2]=%le,!\n",(user->OBB_dxc[rank*3+0]),(user->OBB_dxc[rank*3+1]),(user->OBB_dxc[rank*3+2]));
  
  if ((fabs(xd[0]) <= (user->OBB_dxc[rank*3+0])) &&
      (fabs(xd[1]) <= (user->OBB_dxc[rank*3+1])) &&
      (fabs(xd[2]) <= (user->OBB_dxc[rank*3+2]))) 
    {
      cell_count++;
    }
  
  return cell_count;
}


PetscErrorCode distance_search(Cmpnts p1, Cmpnts p2, Cmpnts p3, Cmpnts p4, Cmpnts p, PetscReal *d)
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



PetscBool ISInsideCell_search(Cmpnts p, Cmpnts cell[8], PetscReal d[6])
{
  // k direction
  distance_search(cell[0], cell[1], cell[2], cell[3], p, &(d[4]));
  distance_search(cell[4], cell[7], cell[6], cell[5], p, &(d[5]));

  // j direction
  distance_search(cell[0], cell[4], cell[5], cell[1], p, &(d[2]));
  distance_search(cell[3], cell[2], cell[6], cell[7], p, &(d[3]));

  // i direction
  distance_search(cell[0], cell[3], cell[7], cell[4], p, &(d[0]));
  distance_search(cell[1], cell[5], cell[6], cell[2], p, &(d[1]));


  if ((d[0]<0)||(d[1]<0)||(d[2]<0)||(d[3]<0)||(d[4]<0)||(d[5]<0)) return(PETSC_FALSE);
  return(PETSC_TRUE);
}



void getQueryPoints(UserCtx *user, int *bnum, int *nints,int *nreals, double **realData,int **index,int **intData,int **sb_index,int rank,int cpu,int bi,int sb,int kk)
{
  int i,j,k,l,m,n,p,i3,jj;
  double xd[3];
  int cell_count=0;
  int *inode,*jnode,*knode;
  int sb_number=0;
  int idx=0;
  
  Cmpnts	***coor;
  DM		da = user[bi].da, fda = user[bi].fda;
  DMDALocalInfo	info;
  PetscInt	xs, ys, zs, xe, ye, ze;
  PetscInt	mz, my, mx;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	gxs, gxe, gys, gye, gzs, gze;
  PetscInt      nnodes=*nints;
  PetscReal	***bcs;


  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  
 
  DMDAVecGetArray(fda, user[bi].Coor, &coor);/// should change to sb
  DMDAVecGetArray(user[bi].da, user[bi].Interface, &bcs);

  for (k=zs; k<ze; k++){
    for (j=ys; j<ye; j++){
      for (i=xs; i<xe; i++){
  	if(bcs[k][j][i]<1.e-6){
  	  nnodes++;
  	}
      }
    }
  }

  //PetscPrintf(PETSC_COMM_SELF, "bi:%d nnodes:%d!\n",bi,nnodes);
  
  inode=(int *)malloc(sizeof(int)*nnodes);
  jnode=(int *)malloc(sizeof(int)*nnodes);
  knode=(int *)malloc(sizeof(int)*nnodes);
  

  ////////////////////////////////////////////////////////// boundary nodes
  for (k=zs; k<ze; k++){
    for (j=ys; j<ye; j++){
      for (i=xs; i<xe; i++){
	if(bcs[k][j][i]<1.e-6 ){
	  
	  for(jj=0;jj<3;jj++) xd[jj]=0;
	  //
	  for(jj=0;jj<3;jj++)
	    xd[jj]=(coor[k][j][i].x-user[sb].OBB_xc[cpu*3+0])*user[sb].OBB_ori_vec[cpu*9+3*0+jj]+(coor[k][j][i].y-user[sb].OBB_xc[cpu*3+1])*user[sb].OBB_ori_vec[cpu*9+3*1+jj]+(coor[k][j][i].z-user[sb].OBB_xc[cpu*3+2])*user[sb].OBB_ori_vec[cpu*9+3*2+jj];
	  

	  if ((fabs(xd[0]) <= (user[sb].OBB_dxc[cpu*3+0])) &&
	      (fabs(xd[1]) <= (user[sb].OBB_dxc[cpu*3+1])) &&
	      (fabs(xd[2]) <= (user[sb].OBB_dxc[cpu*3+2])))
	    {
	      inode[*nints]=i;
	      jnode[*nints]=j;
	      knode[*nints]=k;
	      (*nints)++;
	      (*nreals)+=3;
	      sb_number++;
	    }
	
	}
      }
    }
  }

  //
  //(*intData)=(int *)malloc(sizeof(int)*(*nints));
  
  if(kk==1){
    (*realData)=(double *)malloc(sizeof(double)*(*nreals));/////////reallocate the current array
    (*index)=(int *)malloc(sizeof(int)*(block_number));
    (*intData)=(int *)malloc(sizeof(int)*(*nreals));
    (*sb_index)=(int *)malloc(sizeof(int)*(*nints));
    //
    for(m=0;m<block_number;m++){
      (*index)[m]=0;
    }
    //
  }
  else if(kk!=1 && sb_number){
    (*realData)=(double*) realloc ((*realData), sizeof(double)*(*nreals));
    (*intData)=  (int*)realloc ((*intData), sizeof(int)*(*nreals));
    (*sb_index)= (int*) realloc ((*sb_index), sizeof(int)*(*nints));
  }
  //
  m=0;

  if(sb_number>0)
    (*index)[sb]= sb_number;

  for(m=0;m<sb;m++){
    idx+=((*index)[m]);
  }

  
  for(p=idx;p<idx+sb_number;p++)////// problem in multi-block
    {
      i=inode[p];
      j=jnode[p];
      k=knode[p];
      //
      (*intData)[3*p]=i;
      (*intData)[3*p+1]=j;
      (*intData)[3*p+2]=k;
      //
      (*realData)[3*p]=coor[k][j][i].x;
      (*realData)[3*p+1]=coor[k][j][i].y;
      (*realData)[3*p+2]=coor[k][j][i].z;
      //
      (*sb_index)[p]=sb;
    }
  //
  
  *bnum=1;
  
  DMDAVecRestoreArray(fda, user[bi].Coor, &coor);
  DMDAVecRestoreArray(user[bi].da, user[bi].Interface, &bcs);
  free(inode);
  free(jnode);
  free(knode);
  
}




void sendRecvPackets(PACKET *sndPack,PACKET *rcvPack,int rank,int bi)
{
  int i;
  int *scount,*rcount;
  int tag,irnum;
  MPI_Request *request;
  MPI_Status *status;
  //
  scount=(int *)malloc(2*sizeof(int)*size);
  rcount=(int *) malloc(2*sizeof(int)*size);
  request=(MPI_Request *) malloc(sizeof(MPI_Request)*2*(size+size));
  status=(MPI_Status *) malloc(sizeof(MPI_Status)*2*(size+size));
  //
  for(i=0;i<size;i++){
    scount[2*i]=sndPack[i].nints;			
    scount[2*i+1]=sndPack[i].nreals;
    //
    rcount[2*i]=0;
    rcount[2*i+1]=0;
  }
  //
  irnum=0;
  tag=1;
  //
  for(i=0;i<size;i++){
    if(i!=rank)
      MPI_Irecv(&(rcount[2*i]),2,MPI_INT,rcvMap[i],tag,MPI_COMM_WORLD,&request[irnum++]);
  }
  //
  for(i=0;i<size;i++){
    if(i!=rank)
    MPI_Isend(&(scount[2*i]),2,MPI_INT,sndMap[i],tag,MPI_COMM_WORLD,&request[irnum++]);
  }
  //
  MPI_Waitall(irnum,request,status);
  for(i=0;i<size;i++)
    {
      rcvPack[i].nints=rcount[2*i];
      rcvPack[i].nreals=rcount[2*i+1];
    }
  //
  irnum=0;
  for(i=0;i<size;i++)
    {
      if (rcvPack[i].nints > 0) {
	tag=1;
	rcvPack[i].index=(int *) malloc(sizeof(int)*block_number);
	MPI_Irecv(rcvPack[i].index,block_number,
		  MPI_INT,rcvMap[i],
		  tag,MPI_COMM_WORLD,&request[irnum++]);
      }
      if (rcvPack[i].nreals > 0) {
	tag=2;
	rcvPack[i].realData=(REAL *) malloc(sizeof(REAL)*rcvPack[i].nreals);
	MPI_Irecv(rcvPack[i].realData,rcvPack[i].nreals,
		  MPI_DOUBLE,rcvMap[i],
		  tag,MPI_COMM_WORLD,&request[irnum++]);
      }
    }
  //
  for(i=0;i<size;i++)
    {
      if (sndPack[i].nints > 0){
	tag=1;
	MPI_Isend(sndPack[i].index,block_number,
		  MPI_INT,sndMap[i],
		  tag,MPI_COMM_WORLD,&request[irnum++]);
      }
      if (sndPack[i].nreals > 0){
	tag=2;
	MPI_Isend(sndPack[i].realData,sndPack[i].nreals,
		  MPI_DOUBLE,sndMap[i],
		  tag,MPI_COMM_WORLD,&request[irnum++]);
      }
    }
  MPI_Waitall(irnum,request,status);
  //

  free(scount);
  free(rcount);
  free(request);
  free(status);
}



void sendRecvPacketsRev(PACKET *sndPack,PACKET *rcvPack,int rank,int bi)
{
  int i;
  int *scount,*rcount;
  int tag,irnum;
  MPI_Request *request;
  MPI_Status *status;
  //
  scount=(int *)malloc(2*sizeof(int)*size);
  rcount=(int *) malloc(2*sizeof(int)*size);
  request=(MPI_Request *) malloc(sizeof(MPI_Request)*2*(size+size));
  status=(MPI_Status *) malloc(sizeof(MPI_Status)*2*(size+size));
  //
  for(i=0;i<size;i++){
    scount[2*i]=sndPack[i].nints;			
    scount[2*i+1]=sndPack[i].nreals;
    //
    rcount[2*i]=0;
    rcount[2*i+1]=0;
  }
  //
  irnum=0;
  tag=1;
  //
  for(i=0;i<size;i++){
    if(i!=rank)
      MPI_Irecv(&(rcount[2*i]),2,MPI_INT,rcvMap[i],tag,MPI_COMM_WORLD,&request[irnum++]);
  }
  //
  for(i=0;i<size;i++){
    if(i!=rank)
      MPI_Isend(&(scount[2*i]),2,MPI_INT,sndMap[i],tag,MPI_COMM_WORLD,&request[irnum++]);
  }
  //
  MPI_Waitall(irnum,request,status);
  for(i=0;i<size;i++)
    {
      rcvPack[i].nints=rcount[2*i];
      rcvPack[i].nreals=rcount[2*i+1];
    }
  //
  irnum=0;
  for(i=0;i<size;i++)
    {
      if (rcvPack[i].nreals > 0) {
	tag=1;
	rcvPack[i].intData=(int *) malloc(sizeof(int)*rcvPack[i].nreals);
	MPI_Irecv(rcvPack[i].intData,rcvPack[i].nreals,
		  MPI_INT,rcvMap[i],
		  tag,MPI_COMM_WORLD,&request[irnum++]);
      }

      if (rcvPack[i].nreals > 0) {
	tag=2;
	rcvPack[i].realData=(REAL *) malloc(sizeof(REAL)*rcvPack[i].nreals);
	MPI_Irecv(rcvPack[i].realData,rcvPack[i].nreals,
		  MPI_DOUBLE,rcvMap[i],
		  tag,MPI_COMM_WORLD,&request[irnum++]);
      }
      
    }
  //
  for(i=0;i<size;i++)
    {
      if (sndPack[i].nreals > 0){
	tag=1;
	MPI_Isend(sndPack[i].intData,sndPack[i].nreals,
		  MPI_INT,sndMap[i],
		  tag,MPI_COMM_WORLD,&request[irnum++]);
      }

      if (sndPack[i].nreals > 0){
	tag=2;
	MPI_Isend(sndPack[i].realData,sndPack[i].nreals,
		  MPI_DOUBLE,sndMap[i],
		  tag,MPI_COMM_WORLD,&request[irnum++]);
      }
      
    }
  MPI_Waitall(irnum,request,status);
  //

  free(scount);
  free(rcount);
  free(request);
  free(status);
}


/* void setMap(int ns,int nr, int *snd,int *rcv) */
/* { */
/*   int i; */
/*   // */
/*   if (sndMap) free(sndMap); sndMap=NULL; */
/*   if (rcvMap) free(rcvMap); rcvMap=NULL; */
/*   // */
/*   nsend=ns; */
/*   nrecv=nr; */
/*   sndMap=(int *) malloc(sizeof(int)*nsend); */
/*   rcvMap=(int *) malloc(sizeof(int)*nrecv); */
/*   // */
/*   for(i=0;i<nsend;i++) sndMap[i]=snd[i]; */
/*   for(i=0;i<nrecv;i++) rcvMap[i]=rcv[i]; */
/* } */


void setMap(int ns,int nr, int *snd,int *rcv)
{
  int i;
  int *sndMap1,*rcvMap1;
  //
  int n = sizeof(sndMap)/sizeof(int);
  sndMap1=(int *) malloc(sizeof(int)*n);
  rcvMap1=(int *) malloc(sizeof(int)*n);
  
  for(i=0;i<n;i++){
    sndMap1[i]=sndMap[i];
    rcvMap1[i]=rcvMap[i];
  }
  
  //////////////////////////////////////
  free(sndMap); sndMap=NULL;
  free(rcvMap); rcvMap=NULL;
  //
  nsend=ns;
  nrecv=nr;
  sndMap=(int *) malloc(sizeof(int)*nsend);
  rcvMap=(int *) malloc(sizeof(int)*nrecv);
  //PetscPrintf(PETSC_COMM_SELF,"size:%i int:%i\n",sizeof(rcvMap),sizeof(int));
  //
  for(i=0;i<nsend;i++) sndMap[i]=sndMap1[i];
  for(i=0;i<nrecv;i++) rcvMap[i]=rcvMap1[i];
  
  // PetscPrintf(PETSC_COMM_SELF,"size:%i int:%i\n",sizeof(rcvMap),sizeof(int));
  
  free(sndMap1);
  free(rcvMap1);
}



void getMap(int *ns, int *nr, int **snd,int **rcv)
{
  *ns=nsend;
  *nr=nrecv;
  
  *snd=sndMap;
  *rcv=rcvMap;
  return;
}


void clearPackets(PACKET *sndPack, PACKET *rcvPack)
{
  int i;
  //
  // free Send and recv data
  //
  for(i=0;i<size;i++)
    {
      if (sndPack[i].nints > 0) free(sndPack[i].index);
      if (sndPack[i].nreals > 0) free(sndPack[i].realData);
      if (sndPack[i].nreals > 0) free(sndPack[i].intData);
      sndPack[i].index=NULL;
      sndPack[i].realData=NULL;
      sndPack[i].intData=NULL;
      sndPack[i].nints=sndPack[i].nreals=0;
    }
  for(i=0;i<size;i++)
    {
      if (rcvPack[i].nints > 0) free(rcvPack[i].index);
      if (rcvPack[i].nreals > 0) free(rcvPack[i].realData);
      if (rcvPack[i].nreals > 0) free(rcvPack[i].intData);
      rcvPack[i].index=NULL;
      rcvPack[i].realData=NULL;
      rcvPack[i].nints=rcvPack[i].nreals=0;
    }
  //
}


void initPackets(PACKET *sndPack, PACKET *rcvPack)
{
  int i;
  //
  for(i=0;i<size;i++)     
    {
      sndPack[i].nints=sndPack[i].nreals=0;
       sndPack[i].index=NULL;
      sndPack[i].realData=NULL;
      sndPack[i].intData=NULL;
    }
  //
  for(i=0;i<size;i++)     
    {
      rcvPack[i].nints=rcvPack[i].nreals=0;
      rcvPack[i].index=NULL;
      rcvPack[i].realData=NULL;
      rcvPack[i].intData=NULL;
    }
  //
}

void clearPacketsAll(PACKET *sndPackAll[block_number], PACKET *rcvPackAll[block_number])
{
  int i,bi;
  //
  // free Send and recv data
  //
  for (bi=0;bi<block_number;bi++){
    for(i=0;i<size;i++)
      {
	if (sndPackAll[bi][i].nints > 0) free(sndPackAll[bi][i].index);
	if (sndPackAll[bi][i].nreals > 0) free(sndPackAll[bi][i].realData);
	if (sndPackAll[bi][i].nreals > 0) free(sndPackAll[bi][i].intData);
	sndPackAll[bi][i].index=NULL;
	sndPackAll[bi][i].realData=NULL;
	sndPackAll[bi][i].intData=NULL;
	sndPackAll[bi][i].nints=sndPackAll[bi][i].nreals=0;
      }
    for(i=0;i<size;i++)
      {
	if (rcvPackAll[bi][i].nints > 0) free(rcvPackAll[bi][i].index);
	if (rcvPackAll[bi][i].nreals > 0) free(rcvPackAll[bi][i].realData);
	if (rcvPackAll[bi][i].nreals > 0) free(rcvPackAll[bi][i].intData);
	rcvPackAll[bi][i].index=NULL;
	rcvPackAll[bi][i].realData=NULL;
	rcvPackAll[bi][i].intData=NULL;
	rcvPackAll[bi][i].nints=rcvPackAll[bi][i].nreals=0;
      }
  }
  //
}


void initPacketsAll(PACKET *sndPackAll[block_number], PACKET *rcvPackAll[block_number])
{
  int i,bi;
  //
  for (bi=0;bi<block_number;bi++){
    for(i=0;i<size;i++)     
      {
	sndPackAll[bi][i].nints=sndPackAll[bi][i].nreals=0;
	sndPackAll[bi][i].index=NULL;
	sndPackAll[bi][i].realData=NULL;
	sndPackAll[bi][i].intData=NULL;
      }
    //
    for(i=0;i<size;i++)     
      {
	rcvPackAll[bi][i].nints=rcvPackAll[bi][i].nreals=0;
	rcvPackAll[bi][i].index=NULL;
	rcvPackAll[bi][i].realData=NULL;
	rcvPackAll[bi][i].intData=NULL;
      }
    //
  }
}


void copyPacketAll(PACKET *sndPackAll[block_number], PACKET *rcvPackAll[block_number],int **sndMapAll,int **rcvMapAll,PACKET *sndPack, PACKET *rcvPack,int *sndMap,int *rcvMap,int bi,int rank)
{
  int i,cpu,block;
  //
  //////////////////////////////////////////////////////////// Allocation
  
  for(i=0;i<size;i++) {
    sndPackAll[bi][i].nreals=sndPack[i].nreals;
    sndPackAll[bi][i].nints=sndPack[i].nints;
    //
    rcvPackAll[bi][i].nreals=rcvPack[i].nreals;
    rcvPackAll[bi][i].nints=rcvPack[i].nints;
    
    sndPackAll[bi][i].realData=(REAL *) malloc(sizeof(REAL)*sndPackAll[bi][i].nreals);
    //sndPackAll[bi][i].intData=(int *) malloc(sizeof(int)*sndPackAll[bi][i].nreals);
    sndPackAll[bi][i].index=(int *)malloc(sizeof(int)*(block_number));
    
    rcvPackAll[bi][i].realData=(REAL *) malloc(sizeof(REAL)*rcvPackAll[bi][i].nreals);
    rcvPackAll[bi][i].intData=(int *) malloc(sizeof(int)*rcvPackAll[bi][i].nreals);
    rcvPackAll[bi][i].index=(int *)malloc(sizeof(int)*(block_number));
    
	
  }
  
    
  ////////////////////////////////////////////////////////////
  for(cpu=0;cpu<size;cpu++) {
    for(i=0;i<sndPackAll[bi][cpu].nreals;i++) {
      sndPackAll[bi][cpu].realData[i]=sndPack[cpu].realData[i];
      //sndPackAll[bi][cpu].intData[i]=sndPack[cpu].intData[i];
    }
  }

  for(cpu=0;cpu<size;cpu++) {
    for(i=0;i<rcvPackAll[bi][cpu].nreals;i++) {
      rcvPackAll[bi][cpu].realData[i]=rcvPack[cpu].realData[i];
    }
  }
  
  //
  for(cpu=0;cpu<size;cpu++) {
    for(i=0;i<block_number;i++) {
      if(sndPackAll[bi][cpu].nreals)
      	sndPackAll[bi][cpu].index[i]=sndPack[cpu].index[i];
      if(rcvPackAll[bi][cpu].nreals)
      	rcvPackAll[bi][cpu].index[i]=rcvPack[cpu].index[i];
    }
  }

  for(cpu=0;cpu<size;cpu++) {
    for(i=0;i<block_number;i++) {
      sndMapAll[bi][cpu]=sndMap[cpu];
      rcvMapAll[bi][cpu]=rcvMap[cpu];
    }
  }
  
}



void copyPacketAllRev(PACKET *rcvPackAllRev[block_number],PACKET *rcvPack,int bi,int rank)
{
  int i,cpu,block;
  //
  //////////////////////////////////////////////////////////// Allocation
  
  for(i=0;i<size;i++) {
    //
    rcvPackAllRev[bi][i].nreals=rcvPack[i].nreals;
    rcvPackAllRev[bi][i].nints=rcvPack[i].nints;
    
    rcvPackAllRev[bi][i].realData=(REAL *) malloc(sizeof(REAL)*rcvPackAllRev[bi][i].nreals);
    rcvPackAllRev[bi][i].intData=(int *) malloc(sizeof(int)*rcvPackAllRev[bi][i].nreals);
  }
  
  ////////////////////////////////////////////////////////////

  for(cpu=0;cpu<size;cpu++) {
    for(i=0;i<rcvPackAllRev[bi][cpu].nreals;i++) {
      rcvPackAllRev[bi][cpu].realData[i]=rcvPack[cpu].realData[i];
      rcvPackAllRev[bi][cpu].intData[i]=rcvPack[cpu].intData[i];
    }
  }
  
  
}


void ReturnCoeff(PACKET *rcvPackAll[block_number],PACKET *sndPack,int bi,int rank)
{
  int i,cpu,block;

  for(cpu=0;cpu<size;cpu++) {
    if(rcvPackAll[bi][cpu].nreals>0)
      for(i=0;i<rcvPackAll[bi][cpu].nreals;i++) {
	sndPack[cpu].realData[i]=rcvPackAll[bi][cpu].realData[i];
	sndPack[cpu].intData[i]=rcvPackAll[bi][cpu].intData[i];
      }
  }
  
  //
  for(cpu=0;cpu<size;cpu++) {
    sndPack[cpu].nreals=rcvPackAll[bi][cpu].nreals;
    sndPack[cpu].nints=rcvPackAll[bi][cpu].nints;
  }
  
}



void find_sb(int **index,int *sb,int p){
  int i=-1;
  int idx=0;  
  
  do {
    i++;
    idx+=((*index)[i]);
    if(i>block_number){ 
      PetscPrintf(PETSC_COMM_SELF,"ERROR in finding sb -----------------------------------------------------------------------\n ERROR in finding sb -----------------------------------------------------------------------\n");
      break;
    }
  }while(p>(idx-1));
  
  *sb=i;
  
}

void find_host(UserCtx *user,Cmpnts p,int sb,int rank,int num,double **realData,int **intData,PetscInt ids[dKM][dJM][dIM][10],PetscInt jds[dKM][dJM][dIM][10],PetscInt kds[dKM][dJM][dIM][10],PetscInt ide[dKM][dJM][dIM][10],PetscInt jde[dKM][dJM][dIM][10],PetscInt kde[dKM][dJM][dIM][10]){
  
  PetscReal d[6];
  Cmpnts cell[8];
  PetscInt bh, si, sj, sk;
  Cmpnts ***host;
  PetscReal epss= 1e-8;
  PetscReal ddx, ddy, ddz;
  PetscInt dI, dJ, dK;
  PetscInt out_control=0;
  int count=0,count_tot=0,total_num=0;
  int si_old=-1,sj_old=-1,sk_old=-1;
  PetscBool found=PETSC_FALSE;
  PetscReal ***nvert;
  PetscReal x,y,z;
  PetscReal ibmval=5.0;
  
  DMDAGetLocalInfo(user[sb].da, &(user[sb].info));
  DMDALocalInfo	info = user[sb].info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  
  DMDAVecGetArray(user[sb].fda, user[sb].cent_search, &(host));
  DMDAVecGetArray(user[sb].da, user[sb].lNvert, &nvert);

  
  ddx = (user[sb].Max_X-user[sb].Min_X)/(double)(dIM);
  ddy = (user[sb].Max_Y-user[sb].Min_Y)/(double)(dJM);
  ddz = (user[sb].Max_Z-user[sb].Min_Z)/(double)(dKM);
  
  //PetscPrintf(PETSC_COMM_SELF,"p.x=%le ,p.x=%le,p.x=%le \n",p.x,p.y,p.z);
  
  dI = floor((p.x - user[sb].Min_X) / ddx);
  dJ = floor((p.y - user[sb].Min_Y) / ddy);
  dK = floor((p.z - user[sb].Min_Z) / ddz);
  
  if((dI<0||dI>=dIM)||(dJ<0||dJ>=dJM)||(dK<0||dK>=dKM)) goto nextp;
  if((kds[dK][dJ][dI][sb]>kde[dK][dJ][dI][sb])||(jds[dK][dJ][dI][sb]>jde[dK][dJ][dI][sb])||(ids[dK][dJ][dI][sb]>ide[dK][dJ][dI][sb])) goto nextp;
  
  sk=kds[dK][dJ][dI][sb];sj=jds[dK][dJ][dI][sb];si=ids[dK][dJ][dI][sb];
  
  
  do{
    
    if(d[0]>epss && d[1]<-epss && si<(ide[dK][dJ][dI][sb]-1))si++;else if(d[0]<-epss && d[1]>epss && si>ids[dK][dJ][dI][sb])si--;
    if(d[2]>epss && d[3]<-epss && sj<(jde[dK][dJ][dI][sb]-1))sj++;else if(d[2]<-epss && d[3]>epss && sj>jds[dK][dJ][dI][sb])sj--;
    if(d[4]>epss && d[5]<-epss && sk<(kde[dK][dJ][dI][sb])-1)sk++;else if(d[4]<-epss && d[5]>epss && sk>kds[dK][dJ][dI][sb])sk--;
    
    cell[0] = host[sk  ][sj  ][si  ];
    cell[1] = host[sk  ][sj  ][si+1];
    cell[2] = host[sk  ][sj+1][si+1];
    cell[3] = host[sk  ][sj+1][si  ];
    
    cell[4] = host[sk+1][sj  ][si  ];
    cell[5] = host[sk+1][sj  ][si+1];
    cell[6] = host[sk+1][sj+1][si+1];
    cell[7] = host[sk+1][sj+1][si  ];
    
    if((si==si_old) && (sj==sj_old) && (sk==sk_old)){
      //PetscPrintf(PETSC_COMM_SELF,"p.x=%le ,p.x=%le,p.x=%le \n",p.x,p.y,p.z);
      //PetscPrintf(PETSC_COMM_SELF, "rank=%d si %i %i %i \n", rank, si, sj, sk);
      //PetscPrintf(PETSC_COMM_SELF, "rank=%d si %i %i %i %i %i %i \n", rank, ids[dK][dJ][dI][sb], jds[dK][dJ][dI][sb], kds[dK][dJ][dI][sb], ide[dK][dJ][dI][sb], jde[dK][dJ][dI][sb], kde[dK][dJ][dI][sb]);
      // PetscPrintf(PETSC_COMM_SELF, "rank=%d si %le %le %le \n", rank, user[sb].Min_X, user[sb].Min_Y,user[sb].Min_Z);
      count++;
      count_tot++;
    }
    /////////////////////////////////////////
    if(si<ids[dK][dJ][dI][sb] || sj<jds[dK][dJ][dI][sb] || sk<kds[dK][dJ][dI][sb]|| si>(ide[dK][dJ][dI][sb]-1) || sj>(jde[dK][dJ][dI][sb]-1) || sk>(kde[dK][dJ][dI][sb]-1)) out_control=1; 
    if((count_tot>2)||out_control||total_num>300) goto Sctrl;
    /////////////////////////////////////////
    si_old=si;
    sj_old=sj;
    sk_old=sk;
    
    if(count){
      
      if(d[0]<-epss && d[1]<-epss){
	if(d[0]<-epss && si>xs) si--;;
      }
      if(d[2]<-epss && d[3]<-epss){
	if(d[2]<-epss && sj>ys) sj--;
      }
      if(d[4]<-epss && d[5]<-epss){
	if(d[4]<-epss && sk>zs) sk--;
      }
      
      count=0;
    }
    
    total_num++;
    
  }while(ISInsideCell_search(p, cell, d)==PETSC_FALSE);
  

  if(ISInsideCell_search(p, cell, d) &&
     nvert[sk][sj][si]<ibmval && nvert[sk][sj+1][si]<ibmval && nvert[sk][sj][si+1]<ibmval && nvert[sk][sj+1][si+1]<ibmval &&
     nvert[sk+1][sj][si]<ibmval && nvert[sk+1][sj+1][si]<ibmval && nvert[sk+1][sj][si+1]<ibmval && nvert[sk+1][sj+1][si+1]<ibmval) {
    
    found = PETSC_TRUE;
    //////////////////////////////   change realData for returning the packet
    
    x = d[0] / (d[0] + d[1]);
    y = d[2] / (d[2] + d[3]);
    z = d[4] / (d[4] + d[5]);


    (*realData)[3*num]=x;
    (*realData)[3*num+1]=y;
    (*realData)[3*num+2]=z;
    //
    (*intData)[3*num]=si;
    (*intData)[3*num+1]=sj;
    (*intData)[3*num+2]=sk;
    //////////////////////////////
    
  }
  
  
  //PetscPrintf(PETSC_COMM_SELF, "sk %i %i %i rank:%i d[0]=%le,d[1]=%le,d[2]=%le,d[3]=%le,d[4]=%le,d[5]=%le\n", si, sj, sk,rank,d[0],d[1],d[2],d[3],d[4],d[5]);
  
  //////////////////////////////////////////////////////////////// Search control cell
 Sctrl: if(count_tot>2||out_control||total_num>300) {
    
    for (sk=kds[dK][dJ][dI][sb]; sk<kde[dK][dJ][dI][sb]; sk++) {
      for (sj=jds[dK][dJ][dI][sb]; sj<jde[dK][dJ][dI][sb]; sj++) {
	for (si=ids[dK][dJ][dI][sb]; si<ide[dK][dJ][dI][sb]; si++) {
	  
	  cell[0] = host[sk  ][sj  ][si  ];
	  cell[1] = host[sk  ][sj  ][si+1];
	  cell[2] = host[sk  ][sj+1][si+1];
	  cell[3] = host[sk  ][sj+1][si  ];
	  
	  cell[4] = host[sk+1][sj  ][si  ];
	  cell[5] = host[sk+1][sj  ][si+1];
	  cell[6] = host[sk+1][sj+1][si+1];
	  cell[7] = host[sk+1][sj+1][si  ]; 
	  
	  if(ISInsideCell_search(p, cell, d) &&
	     nvert[sk][sj][si]<ibmval && nvert[sk][sj+1][si]<ibmval && nvert[sk][sj][si+1]<ibmval && nvert[sk][sj+1][si+1]<ibmval &&
	     nvert[sk+1][sj][si]<ibmval && nvert[sk+1][sj+1][si]<ibmval && nvert[sk+1][sj][si+1]<ibmval && nvert[sk+1][sj+1][si+1]<ibmval) {
	    found = PETSC_TRUE;
	    //PetscPrintf(PETSC_COMM_SELF, "sk1 %i %i %i rank:%i d[0]=%le,d[1]=%le,d[2]=%le,d[3]=%le,d[4]=%le,d[5]=%le\n", si, sj, sk,rank,d[0],d[1],d[2],d[3],d[4],d[5]);

	    //////////////////////////////   change realData for returning the packet

	    x = d[0] / (d[0] + d[1]);
	    y = d[2] / (d[2] + d[3]);
	    z = d[4] / (d[4] + d[5]);
	    
	    
	    (*realData)[3*num]=x;
	    (*realData)[3*num+1]=y;
	    (*realData)[3*num+2]=z;
	    //
	    (*intData)[3*num]=si;
	    (*intData)[3*num+1]=sj;
	    (*intData)[3*num+2]=sk;
	    //////////////////////////////
	    goto nextp;
	    
	  }
	  
	}
      }
    }
    
  }
  //////////////////////////////////////////////////////////////// 
 nextp:if (!found) { 
    // PetscPrintf(PETSC_COMM_SELF, "rank=%d si %le %le %le \n", rank, user[sb].Min_X, user[sb].Min_Y,user[sb].Min_Z);
    (*realData)[3*num]=-200.0;//-200.00
    (*realData)[3*num+1]=0.0;
    (*realData)[3*num+2]=0.0;
    //
    (*intData)[3*num]=-10;
    (*intData)[3*num+1]=-10;
    (*intData)[3*num+2]=-10;
    //PetscPrintf(PETSC_COMM_SELF, "NOT FOUND! \n");
  }
  
  DMDAVecRestoreArray(user[sb].fda, user[sb].cent_search, &(host));
  DMDAVecRestoreArray(user[sb].da, user[sb].lNvert, &nvert);
  
}


PetscErrorCode CenterNodeCoor(UserCtx *user)
{
  
  DM		da = user->da, fda = user->fda;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe=xe-1;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;

  PetscInt       i, j, k, IM, JM, KM;
  Cmpnts ***coor, ***cent;


  DMDAVecGetArray(user->fda, user->Coor, &coor);
  DMDAVecGetArray(user->fda, user->gcent_search, &cent);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	cent[k][j][i].x = 0.125 *
	  (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
	   coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x +
	   coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x +
	   coor[k-1][j  ][i-1].x + coor[k-1][j-1][i-1].x);
	cent[k][j][i].y = 0.125 *
	  (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
	   coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y +
	   coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y +
	   coor[k-1][j  ][i-1].y + coor[k-1][j-1][i-1].y);
	cent[k][j][i].z = 0.125 *
	  (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
	   coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z +
	   coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z +
	   coor[k-1][j  ][i-1].z + coor[k-1][j-1][i-1].z);

      }
    }
  }
  
  // Ghost nodes

  for (k=zs; k<ze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	
	if(k==0){
	  cent[k][j][i].x = 2.*0.25 * (coor[k  ][j  ][i  ].x +
				       coor[k  ][j  ][i-1].x +
				       coor[k  ][j-1][i  ].x +
				       coor[k  ][j-1][i-1].x )- cent[k+1][j][i].x;
	  
	  cent[k][j][i].y = 2.*0.25 * (coor[k  ][j  ][i  ].y +
				       coor[k  ][j  ][i-1].y +
				       coor[k  ][j-1][i  ].y +
				       coor[k  ][j-1][i-1].y )- cent[k+1][j][i].y;
	  
	  cent[k][j][i].z = 2.*0.25 * (coor[k  ][j  ][i  ].z +
				       coor[k  ][j  ][i-1].z +
				       coor[k  ][j-1][i  ].z +
				       coor[k  ][j-1][i-1].z )- cent[k+1][j][i].z;
	}
	
	if(k==user->KM){
	  cent[k][j][i].x = 2.*0.25 * (coor[k-1][j  ][i  ].x +
				       coor[k-1][j  ][i-1].x +
				       coor[k-1][j-1][i  ].x +
				       coor[k-1][j-1][i-1].x )- cent[k-1][j][i].x;
	  
	  cent[k][j][i].y = 2.*0.25 * (coor[k-1][j  ][i  ].y +
				       coor[k-1][j  ][i-1].y +
				       coor[k-1][j-1][i  ].y +
				       coor[k-1][j-1][i-1].y )- cent[k-1][j][i].y;
	  
	  cent[k][j][i].z = 2.*0.25 * (coor[k-1][j  ][i  ].z +
				       coor[k-1][j  ][i-1].z +
				       coor[k-1][j-1][i  ].z +
				       coor[k-1][j-1][i-1].z )- cent[k-1][j][i].z;
	}
	
      }
    }
  }
  

  for (k=lzs; k<lze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=lxs; i<lxe; i++) {
	
	if(j==0){
	  cent[k][j][i].x = 2.*0.25 * (coor[k  ][j  ][i  ].x +
				       coor[k  ][j  ][i-1].x +
				       coor[k-1][j  ][i  ].x +
				       coor[k-1][j  ][i-1].x)- cent[k][j+1][i].x;
	  
	  cent[k][j][i].y = 2.*0.25 * (coor[k  ][j  ][i  ].y +
				       coor[k  ][j  ][i-1].y +
				       coor[k-1][j  ][i  ].y +
				       coor[k-1][j  ][i-1].y)- cent[k][j+1][i].y;
	  
	  cent[k][j][i].z = 2.*0.25 * (coor[k  ][j  ][i  ].z +
				       coor[k  ][j  ][i-1].z +
				       coor[k-1][j  ][i  ].z +
				       coor[k-1][j  ][i-1].z)- cent[k][j+1][i].z;
	}
	
	if(j==user->JM){
	  cent[k][j][i].x = 2.*0.25 * (coor[k  ][j-1  ][i  ].x +
				       coor[k  ][j-1  ][i-1].x +
				       coor[k-1][j-1  ][i  ].x +
				       coor[k-1][j-1  ][i-1].x)- cent[k][j-1][i].x;
	  
	  cent[k][j][i].y = 2.*0.25 * (coor[k  ][j-1  ][i  ].y +
				       coor[k  ][j-1  ][i-1].y +
				       coor[k-1][j-1  ][i  ].y +
				       coor[k-1][j-1  ][i-1].y)- cent[k][j-1][i].y;
	  
	  cent[k][j][i].z = 2.*0.25 * (coor[k  ][j-1  ][i  ].z +
				       coor[k  ][j-1  ][i-1].z +
				       coor[k-1][j-1  ][i  ].z +
				       coor[k-1][j-1  ][i-1].z)- cent[k][j-1][i].z;
	}
	
      }
    }
  }


  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<xe; i++) {

	if(i==0){
	  cent[k][j][i].x = 2.*0.25 * (coor[k  ][j  ][i  ].x +
				       coor[k  ][j-1][i  ].x +
				       coor[k-1][j  ][i  ].x +
				       coor[k-1][j-1][i  ].x)- cent[k][j][i+1].x;
	  
	  cent[k][j][i].y = 2.*0.25 * (coor[k  ][j  ][i  ].y +
				       coor[k  ][j-1][i  ].y +
				       coor[k-1][j  ][i  ].y +
				       coor[k-1][j-1][i  ].y)- cent[k][j][i+1].y;
	  
	  
	  cent[k][j][i].z = 2.*0.25 * (coor[k  ][j  ][i  ].z +
				       coor[k  ][j-1][i  ].z +
				       coor[k-1][j  ][i  ].z +
				       coor[k-1][j-1][i  ].z)- cent[k][j][i+1].z;
	}
	
	if(i==user->IM){
	  cent[k][j][i].x = 2.*0.25 * (coor[k  ][j  ][i-1  ].x +
				       coor[k  ][j-1][i-1  ].x +
				       coor[k-1][j  ][i-1  ].x +
				       coor[k-1][j-1][i-1  ].x)- cent[k][j][i-1].x;
	  
	  cent[k][j][i].y = 2.*0.25 * (coor[k  ][j  ][i-1  ].y +
				       coor[k  ][j-1][i-1  ].y +
				       coor[k-1][j  ][i-1  ].y +
				       coor[k-1][j-1][i-1  ].y)- cent[k][j][i-1].y;
	  
	  
	  cent[k][j][i].z = 2.*0.25 * (coor[k  ][j  ][i-1  ].z +
				       coor[k  ][j-1][i-1  ].z +
				       coor[k-1][j  ][i-1  ].z +
				       coor[k-1][j-1][i-1  ].z)- cent[k][j][i-1].z;
	}
	
      }
    }
  }
  
  


  // lines

  // kmin plane lines
  for (k= zs; k< ze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=lxs; i<lxe; i++) {
	if(k==0 && j==0){
	  
	  cent[k][j][i].x = 2.*0.5*(coor[k][j][i].x+coor[k][j][i-1].x) - cent[k+1][j+1][i].x;
	  cent[k][j][i].y = 2.*0.5*(coor[k][j][i].y+coor[k][j][i-1].y) - cent[k+1][j+1][i].y;
	  cent[k][j][i].z = 2.*0.5*(coor[k][j][i].z+coor[k][j][i-1].z) - cent[k+1][j+1][i].z;
	  
	}
      }
    }
  }
    


  for (k= zs; k< ze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=lxs; i<lxe; i++) {
	if(k==0 && j== user->JM ){

	  cent[k][j][i].x = 2.*0.5*(coor[k][j-1][i].x+coor[k][j-1][i-1].x) - cent[k+1][j-1][i].x;
	  cent[k][j][i].y = 2.*0.5*(coor[k][j-1][i].y+coor[k][j-1][i-1].y) - cent[k+1][j-1][i].y;
	  cent[k][j][i].z = 2.*0.5*(coor[k][j-1][i].z+coor[k][j-1][i-1].z) - cent[k+1][j-1][i].z;

	}
      }
    }
  }
  
  for (k=zs; k<ze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=xs; i<xe; i++) {
	if(k==0 && i== 0){
	  
	  cent[k][j][i].x = 2.*0.5*(coor[k][j][i].x+coor[k][j-1][i].x) - cent[k+1][j][i+1].x;
	  cent[k][j][i].y = 2.*0.5*(coor[k][j][i].y+coor[k][j-1][i].y) - cent[k+1][j][i+1].y;
	  cent[k][j][i].z = 2.*0.5*(coor[k][j][i].z+coor[k][j-1][i].z) - cent[k+1][j][i+1].z;
	  	  
	}
      }
    }
  }
  
  for (k= zs; k< ze; k++) {
    for (j= lys; j< lye; j++) {
      for (i=xs; i<xe; i++) {
	if(k==0 && i== user->IM ){

	  cent[k][j][i].x = 2.*0.5*(coor[k][j][i-1].x+coor[k][j-1][i-1].x) - cent[k+1][j][i-1].x;
	  cent[k][j][i].y = 2.*0.5*(coor[k][j][i-1].y+coor[k][j-1][i-1].y) - cent[k+1][j][i-1].y;
	  cent[k][j][i].z = 2.*0.5*(coor[k][j][i-1].z+coor[k][j-1][i-1].z) - cent[k+1][j][i-1].z;


	}
      }
    }
  }
  
/*   // kmax plane lines */


  for (k= zs; k< ze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=lxs; i<lxe; i++) {
	if(k==user->KM && j==0 ){


	  cent[k][j][i].x = 2.*0.5*(coor[k-1][j][i-1].x+coor[k-1][j][i].x) - cent[k-1][j+1][i].x;
	  cent[k][j][i].y = 2.*0.5*(coor[k-1][j][i-1].y+coor[k-1][j][i].y) - cent[k-1][j+1][i].y;
	  cent[k][j][i].z = 2.*0.5*(coor[k-1][j][i-1].z+coor[k-1][j][i].z) - cent[k-1][j+1][i].z;


	  
	}
      }
    }
  }
  

  for (k= zs; k< ze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=lxs; i<lxe; i++) {
	if(k==user->KM && j==user->JM){

	  cent[k][j][i].x = 2.*0.5*(coor[k-1][j-1][i-1].x+coor[k-1][j-1][i].x) - cent[k-1][j-1][i].x;
	  cent[k][j][i].y = 2.*0.5*(coor[k-1][j-1][i-1].y+coor[k-1][j-1][i].y) - cent[k-1][j-1][i].y;
	  cent[k][j][i].z = 2.*0.5*(coor[k-1][j-1][i-1].z+coor[k-1][j-1][i].z) - cent[k-1][j-1][i].z;


	}
      }
    }
  }


  for (k= zs; k< ze; k++) {
    for (j= lys; j< lye; j++) {
      for (i=xs; i<xe; i++) {
	if(k==user->KM && i==0){
	
	  cent[k][j][i].x = 2.*0.5*(coor[k-1][j-1][i].x+coor[k-1][j][i].x) - cent[k-1][j][i+1].x;
	  cent[k][j][i].y = 2.*0.5*(coor[k-1][j-1][i].y+coor[k-1][j][i].y) - cent[k-1][j][i+1].y;
	  cent[k][j][i].z = 2.*0.5*(coor[k-1][j-1][i].z+coor[k-1][j][i].z) - cent[k-1][j][i+1].z;

	}
      }
    }
  }

  for (k= zs; k< ze; k++) {
    for (j= lys; j< lye; j++) {
      for (i=xs; i<xe; i++) {
	if(k==user->KM && i==user->IM ){

	  cent[k][j][i].x = 2.*0.5*(coor[k-1][j-1][i-1].x+coor[k-1][j][i-1].x) - cent[k-1][j][i-1].x;
	  cent[k][j][i].y = 2.*0.5*(coor[k-1][j-1][i-1].y+coor[k-1][j][i-1].y) - cent[k-1][j][i-1].y;
	  cent[k][j][i].z = 2.*0.5*(coor[k-1][j-1][i-1].z+coor[k-1][j][i-1].z) - cent[k-1][j][i-1].z;


	}
      }
    }
  }

  // jmin lines

  for (k= lzs; k< lze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=xs; i<xe; i++) {
	if(j==0 && i==0 ){

	  cent[k][j][i].x = 2.*0.5*(coor[k-1][j][i].x+coor[k][j][i].x) - cent[k][j+1][i+1].x;
	  cent[k][j][i].y = 2.*0.5*(coor[k-1][j][i].y+coor[k][j][i].y) - cent[k][j+1][i+1].y;
	  cent[k][j][i].z = 2.*0.5*(coor[k-1][j][i].z+coor[k][j][i].z) - cent[k][j+1][i+1].z;


	}
      }
    }
  }

  for (k= lzs; k< lze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=xs; i<xe; i++) {
	if(j==0 && i== user->IM ){

	  cent[k][j][i].x = 2.*0.5*(coor[k-1][j][i-1].x+coor[k][j][i-1].x) - cent[k][j+1][i-1].x;
	  cent[k][j][i].y = 2.*0.5*(coor[k-1][j][i-1].y+coor[k][j][i-1].y) - cent[k][j+1][i-1].y;
	  cent[k][j][i].z = 2.*0.5*(coor[k-1][j][i-1].z+coor[k][j][i-1].z) - cent[k][j+1][i-1].z;


	}
      }
    }
  }

  // jmax lines

  for (k= lzs; k< lze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=xs; i<xe; i++) {
	if(j==user->JM  && i==0 ){

	  cent[k][j][i].x = 2.*0.5*(coor[k-1][j-1][i].x+coor[k][j-1][i].x) - cent[k][j-1][i+1].x;
	  cent[k][j][i].y = 2.*0.5*(coor[k-1][j-1][i].y+coor[k][j-1][i].y) - cent[k][j-1][i+1].y;
	  cent[k][j][i].z = 2.*0.5*(coor[k-1][j-1][i].z+coor[k][j-1][i].z) - cent[k][j-1][i+1].z;


	}
      }
    }
  }

  
  for (k= lzs; k< lze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=xs; i<xe; i++) {
	if(j==user->JM  && i== user->IM ){

	  cent[k][j][i].x = 2.*0.5*(coor[k-1][j-1][i-1].x+coor[k][j-1][i-1].x) - cent[k][j-1][i-1].x;
	  cent[k][j][i].y = 2.*0.5*(coor[k-1][j-1][i-1].y+coor[k][j-1][i-1].y) - cent[k][j-1][i-1].y;
	  cent[k][j][i].z = 2.*0.5*(coor[k-1][j-1][i-1].z+coor[k][j-1][i-1].z) - cent[k][j-1][i-1].z;

	}
      }
    }
  }



/*   // Corners */

  for (k= zs; k< ze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=xs; i<xe; i++) {
	if(j==0 && i== 0 && k==0){
	  cent[k][j][i].x = 2.*coor[k][j][i].x - cent[k+1][j+1][i+1].x;
	  cent[k][j][i].y = 2.*coor[k][j][i].y - cent[k+1][j+1][i+1].y;
	  cent[k][j][i].z = 2.*coor[k][j][i].z - cent[k+1][j+1][i+1].z;
	}
      }
    }
  }


  for (k= zs; k< ze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=xs; i<xe; i++) {
	if(j==0 && i== 0 && k== user->KM){
	  cent[k][j][i].x = 2.*coor[k-1][j][i].x - cent[k-1][j+1][i+1].x;
	  cent[k][j][i].y = 2.*coor[k-1][j][i].y - cent[k-1][j+1][i+1].y;
	  cent[k][j][i].z = 2.*coor[k-1][j][i].z - cent[k-1][j+1][i+1].z;
	}
      }
    }
  }
  

  for (k= zs; k< ze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=xs; i<xe; i++) {
	if(j== user->JM &&i==0 && k==0){
	  cent[k][j][i].x = 2.*coor[k][j-1][i].x - cent[k+1][j-1][i+1].x;
	  cent[k][j][i].y = 2.*coor[k][j-1][i].y - cent[k+1][j-1][i+1].y;
	  cent[k][j][i].z = 2.*coor[k][j-1][i].z - cent[k+1][j-1][i+1].z;
	}
      }
    }
  }

  for (k= zs; k< ze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=xs; i<xe; i++) {
	if(j== user->JM &&i==0 && k== user->KM){
	  cent[k][j][i].x = 2.*coor[k-1][j-1][i].x - cent[k-1][j-1][i+1].x;
	  cent[k][j][i].y = 2.*coor[k-1][j-1][i].y - cent[k-1][j-1][i+1].y;
	  cent[k][j][i].z = 2.*coor[k-1][j-1][i].z - cent[k-1][j-1][i+1].z;
	}
      }
    }
  }


  for (k= zs; k< ze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=xs; i<xe; i++) {
	if(j== 0 &&i==user->IM && k== 0){
	  cent[k][j][i].x = 2.*coor[k][j][i-1].x - cent[k+1][j+1][i-1].x;
	  cent[k][j][i].y = 2.*coor[k][j][i-1].y - cent[k+1][j+1][i-1].y;
	  cent[k][j][i].z = 2.*coor[k][j][i-1].z - cent[k+1][j+1][i-1].z;
	}
      }
    }
  }
  

  for (k= zs; k< ze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=xs; i<xe; i++) {
	if(j== 0 &&i==user->IM && k== user->KM){
	  cent[k][j][i].x = 2.*coor[k-1][j][i-1].x - cent[k-1][j+1][i-1].x;
	  cent[k][j][i].y = 2.*coor[k-1][j][i-1].y - cent[k-1][j+1][i-1].y;
	  cent[k][j][i].z = 2.*coor[k-1][j][i-1].z - cent[k-1][j+1][i-1].z;
	}
      }
    }
  }

 
  for (k= zs; k< ze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=xs; i<xe; i++) {
	if(j==user->JM &&i==user->IM && k==0){
	  cent[k][j][i].x = 2.*coor[k][j-1][i-1].x - cent[k+1][j-1][i-1].x;
	  cent[k][j][i].y = 2.*coor[k][j-1][i-1].y - cent[k+1][j-1][i-1].y;
	  cent[k][j][i].z = 2.*coor[k][j-1][i-1].z - cent[k+1][j-1][i-1].z;
	}
      }
    }
  }
  

  for (k= zs; k< ze; k++) {
    for (j= ys; j< ye; j++) {
      for (i=xs; i<xe; i++) {
	if(j==user->JM &&i==user->IM && k==user->KM){
	  cent[k][j][i].x = 2.*coor[k-1][j-1][i-1].x - cent[k-1][j-1][i-1].x;
	  cent[k][j][i].y = 2.*coor[k-1][j-1][i-1].y - cent[k-1][j-1][i-1].y;
	  cent[k][j][i].z = 2.*coor[k-1][j-1][i-1].z - cent[k-1][j-1][i-1].z;
	}
      }
    }
  }


  DMDAVecRestoreArray(user->fda, user->gcent_search, &cent);
  DMDAVecRestoreArray(user->fda, user->Coor, &coor);
  /////////////
  DMGlobalToLocalBegin(fda, user->gcent_search, INSERT_VALUES,user->cent_search);
  DMGlobalToLocalEnd(fda, user->gcent_search, INSERT_VALUES, user->cent_search);
  
  
}


void getQueryPoints_rank(UserCtx *user, int *nints,int *nreals, double **realData,int **intData_rank_index,int **intData_rank,int rank,int cpu,int bi,int sb)
{
  int i,j,k,l,m,n,p,i3,jj;
  double xd[3];
  int cell_count=0;
  int *inode,*jnode,*knode;
  int sb_number=*nints;
  
  Cmpnts	***coor;
  DM		da = user[bi].da, fda = user[bi].fda;
  DMDALocalInfo	info;
  PetscInt	xs, ys, zs, xe, ye, ze;
  PetscInt	mz, my, mx;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	gxs, gxe, gys, gye, gzs, gze;
  PetscInt      nnodes=0;
  PetscReal     ***bcs;

  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  DMDAVecGetArray(fda, user[bi].Coor, &coor);/// should change to sb
  DMDAVecGetArray(user[bi].da, user[bi].Interface, &bcs); 
  
  for (k=zs; k<ze; k++){
    for (j=ys; j<ye; j++){
      for (i=xs; i<xe; i++){
  	if(bcs[k][j][i]<1.e-6){
  	  nnodes++;
  	}
      }
    }
  }
  
  inode=(int *)malloc(sizeof(int)*nnodes);
  jnode=(int *)malloc(sizeof(int)*nnodes);
  knode=(int *)malloc(sizeof(int)*nnodes);
  
  
  for (k=zs; k<ze; k++){
    for (j=ys; j<ye; j++){
      for (i=xs; i<xe; i++){
	if(bcs[k][j][i]<1.e-6){
	  
	  for(jj=0;jj<3;jj++) xd[jj]=0;
	  //
	  for(jj=0;jj<3;jj++)
	    xd[jj]=(coor[k][j][i].x-user[sb].OBB_xc[cpu*3+0])*user[sb].OBB_ori_vec[cpu*9+3*0+jj]+(coor[k][j][i].y-user[sb].OBB_xc[cpu*3+1])*user[sb].OBB_ori_vec[cpu*9+3*1+jj]+(coor[k][j][i].z-user[sb].OBB_xc[cpu*3+2])*user[sb].OBB_ori_vec[cpu*9+3*2+jj];
	  
	  
  
	  if ((fabs(xd[0]) <= (user[sb].OBB_dxc[cpu*3+0])) &&
	      (fabs(xd[1]) <= (user[sb].OBB_dxc[cpu*3+1])) &&
	      (fabs(xd[2]) <= (user[sb].OBB_dxc[cpu*3+2])))
	    {
	      inode[*nints]=i;
	      jnode[*nints]=j;
	      knode[*nints]=k;
	      (*nints)++;
	      (*nreals)+=3;
	      sb_number++;
	      //PetscPrintf(PETSC_COMM_SELF," i=%i j=%i k=%i ,x=%le , y=%le , z=%le \n",i, j, k, coor[k][j][i].x,coor[k][j][i].y,coor[k][j][i].z);
	    }
	
	}
      }
    }
  }

 

  //
  
    (*realData)=(double *)malloc(sizeof(double)*(*nreals));/////////reallocate the current array
    (*intData_rank_index)=(int *)malloc(sizeof(int)*(*nreals));
    (*intData_rank)=(int *)malloc(sizeof(int)*(*nreals));
  //
  m=0;
  
  for(p=0;p<*nints;p++)
    {
      i=inode[p];
      j=jnode[p];
      k=knode[p];
      //
      (*realData)[m++]=coor[k][j][i].x;
      (*realData)[m++]=coor[k][j][i].y;
      (*realData)[m++]=coor[k][j][i].z;
      //
      (*intData_rank_index)[3*p]=i;
      (*intData_rank_index)[3*p+1]=j;
      (*intData_rank_index)[3*p+2]=k;
    }
  //

    
   DMDAVecRestoreArray(fda, user[bi].Coor, &coor);
   DMDAVecRestoreArray(user[bi].da, user[bi].Interface, &bcs);
   free(inode);
   free(jnode);
   free(knode);
  
}

void move_grid(UserCtx *user,FSInfo *fsi,int ti){
  
  int i,j,k,l,m,n,p,i3,jj;
  Cmpnts	***coor,***coor_init;
  
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs, ys, zs, xe, ye, ze;
  PetscInt	mz, my, mx;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	gxs, gxe, gys, gye, gzs, gze;
  
 
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  
  DMDAVecGetArray(fda, user->Coor, &coor);
  DMDAVecGetArray(fda, user->Coor_init, &coor_init);
  
  /////////////////////////////////////////////// Translation
  PetscReal ax=0.0,ay=0.0;
      
  for (k=gzs; k<gze; k++){
    for (j=gys; j<gye; j++){
      for (i=gxs; i<gxe; i++){

	coor[k][j][i].x=(coor_init[k][j][i].x-ax)*cos(theta_x)-(coor_init[k][j][i].y-ay)*sin(theta_x)+ax;
	coor[k][j][i].y=(coor_init[k][j][i].x-ax)*sin(theta_x)+(coor_init[k][j][i].y-ay)*cos(theta_x)+ay;
	coor[k][j][i].z=coor_init[k][j][i].z;
      }
    }
  }
  PetscPrintf(PETSC_COMM_WORLD, "rotated angle=%f\n", theta_x);

  
  DMDAVecRestoreArray(fda, user->Coor, &coor);
  DMDAVecRestoreArray(fda, user->Coor_init, &coor_init);

  PetscBarrier(NULL);
}



void move_grid_ibm(UserCtx *user,FSInfo *fsi,int ti){
  
  int i,j,k,l,m,n,p,i3,jj;
  Cmpnts	***coor,***coor_init;
  
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs, ys, zs, xe, ye, ze;
  PetscInt	mz, my, mx;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	gxs, gxe, gys, gye, gzs, gze;
  
 
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  
  DMDAVecGetArray(fda, user->Cent, &coor);
  DMDAVecGetArray(fda, user->Cent_init, &coor_init);
  
  /////////////////////////////////////////////// Translation
  PetscReal ax=0.0,ay=0.0;
      
  for (k=zs; k<ze; k++){
    for (j=ys; j<ye; j++){
      for (i=xs; i<xe; i++){
	
	coor[k][j][i].x=(coor_init[k][j][i].x-ax)*cos(theta_x)-(coor_init[k][j][i].y-ay)*sin(theta_x)+ax;
	coor[k][j][i].y=(coor_init[k][j][i].x-ax)*sin(theta_x)+(coor_init[k][j][i].y-ay)*cos(theta_x)+ay;
	coor[k][j][i].z=coor_init[k][j][i].z;
      }
    }
  }

  
  DMDAVecRestoreArray(fda, user->Cent, &coor);
  DMDAVecRestoreArray(fda, user->Cent_init, &coor_init);
  //
  DMGlobalToLocalBegin(fda, user->Cent, INSERT_VALUES, user->lCent);
  DMGlobalToLocalEnd(fda, user->Cent, INSERT_VALUES, user->lCent);

  PetscBarrier(NULL);
}


void grid_out(UserCtx *user, int ti){

  Cmpnts	***coor,***gcoor;
  DMDALocalInfo	info = user->info;
  PetscInt	xs, ys, zs, xe, ye, ze;
  PetscInt	mz, my, mx;
  int i,j,k;
  
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  
  Vec      gCoor;
  DMCreateGlobalVector(user->fda, &(gCoor));
  VecSet(gCoor,0.0);
  //
  DMDAVecGetArray(user->fda, user->Coor, &coor);
  DMDAVecGetArray(user->fda, gCoor, &gcoor);

  for (k=zs; k<ze; k++){
    for (j=ys; j<ye; j++){
      for (i=xs; i<xe; i++){
	
	gcoor[k][j][i].x=coor[k][j][i].x;
	gcoor[k][j][i].y=coor[k][j][i].y;
	gcoor[k][j][i].z=coor[k][j][i].z;
      }
    }
  }
  
  DMDAVecRestoreArray(user->fda, user->Coor, &coor);
  DMDAVecRestoreArray(user->fda, gCoor, &gcoor);
  
  PetscViewer	viewer;
  char filen[90];
  PetscInt N;
  
  VecGetSize(gCoor, &N);
  PetscPrintf(PETSC_COMM_WORLD, "grid out time=%d\n", ti);
  
  sprintf(filen, "grid%5.5d_%1.1d.dat", ti, user->_this);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  VecView( gCoor, viewer);
  PetscViewerDestroy(&viewer);
    
  VecDestroy(&gCoor);
  
}


PetscErrorCode Interpolation_matrix(UserCtx *user,PACKET *rcvPackAllRev[block_number],int ***intData_index,int ***sb_index) // for two block
{
  
  PetscErrorCode       ierr;
  PetscInt             i,j,k;
  char                 filen[80];
  PetscInt             rank,size,bi,level,cpu;
  //DM                   packer;
  DMDALocalInfo        info1,info2,info3;
  PetscViewer          viewer;
  PetscReal           ***Check,***bcs;
  
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);


   /////////////////////////////////////----------------------------///////////////////////////////////////////////// reading the coeficients and host
  
  PetscLogDouble v1,v2,elapsed_time;

  PetscTime(&v1);
  
  PetscInt  hb=0;
  
  for (bi=0; bi<block_number; bi++) {
    
    DMDALocalInfo	info = user[bi].info;
    
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;


    PetscInt d,dg,cll=0;
    PetscInt itfc_number=0;
    PetscReal lval[8];
    PetscInt	row=0,lcol[8];
    PetscInt col[8];
    PetscReal val[8];
    PetscReal x, y, z;
    PetscInt  itfnumber=0,itfnumber1=0;
    PetscInt d_freedom=0;
    PetscInt coeff[block_number];
    PetscInt p;
    PetscInt cpu_host=0;
    PetscInt blk=0;

    DMDAVecGetArray(user[bi].da, user[bi].check, &(Check));
    DMDAVecGetArray(user[bi].da, user[bi].lInterface, &bcs);
    /**************************************************************************************************************************/
    /* Interpolate the Velocity on the interface nodes
       from the host nodes.
       itfc is the velocity at the cell corners
       hostU is the velocity at the cell centers of the host block */
    /**************************************************************************************************************************/
    for(cpu=0;cpu<size;cpu++){/////////////////
      
      for (p=0; p<rcvPackAllRev[bi][cpu].nints; p++) {
	
	x = rcvPackAllRev[bi][cpu].realData[3*p];
	y = rcvPackAllRev[bi][cpu].realData[3*p+1];
	z = rcvPackAllRev[bi][cpu].realData[3*p+2];
	
	
	dg=0;d=0;cll=0;
	
	for (i=0;i<8;i++){
	  lcol[i]=0;
	  col[i]=0;
	  val[i]=0.0;
	}
	
	//itfnumber=0;
	PetscInt ii,jj,kk;  
	
	ii=intData_index[bi][cpu][3*p];
	jj=intData_index[bi][cpu][3*p+1];
	kk=intData_index[bi][cpu][3*p+2];

	////////////////////////////////////// each blank just interpolate from related block
	PetscBool blk_block= PETSC_TRUE;
	hb = sb_index[bi][cpu][p];

	if((bcs[kk][jj][ii]<(-10*hb+1.)&& bcs[kk][jj][ii]>(-10*hb-1.))|| (bcs[kk][jj][ii]< (-1000)) ||(bcs[kk][jj][ii]<0 && bcs[kk][jj][ii]>-7.0)){
	  blk_block=PETSC_TRUE;
	}else{
	  blk_block=PETSC_FALSE;
	}
	  
	  //}	
	//////////////////////////////////////
      
	if(kk>=zs && kk<ze){
	  if(jj>=ys && jj<ye){
	    if(ii>=xs && ii<xe){
		if(rcvPackAllRev[bi][cpu].realData[3*p]>-1.25 && Check[kk][jj][ii]<0 && blk_block)//////////////////////////
		  {
		    
		    hb = sb_index[bi][cpu][p];
			    
	  
		    k=rcvPackAllRev[bi][cpu].intData[3*p+2];j=rcvPackAllRev[bi][cpu].intData[3*p+1];i=rcvPackAllRev[bi][cpu].intData[3*p]; //i,j,k
		    cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    //
		    dg=dg_start[cpu_host]+3*d;
		    for (blk=0; blk<hb; blk++)
		      dg+=dl_start_blk[blk][cpu_host];
		    //
		    lcol[0]=dg;
		    val[0]=(1-x) * (1-y) * (1-z);
		    itfnumber++;
		    
		    
		    k=rcvPackAllRev[bi][cpu].intData[3*p+2];j=rcvPackAllRev[bi][cpu].intData[3*p+1];i=rcvPackAllRev[bi][cpu].intData[3*p]+1; //i+1,j,k
		    cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    //
		    dg=dg_start[cpu_host]+3*d;
		    for (blk=0; blk<hb; blk++)
		      dg+=dl_start_blk[blk][cpu_host];
		    //
		    lcol[1]=dg;
		    val[1]=x * (1-y) * (1-z);
		    itfnumber++;
		    
		    
		    k=rcvPackAllRev[bi][cpu].intData[3*p+2];j=rcvPackAllRev[bi][cpu].intData[3*p+1]+1;i=rcvPackAllRev[bi][cpu].intData[3*p]; //i,j+1,k
		    cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    //
		    dg=dg_start[cpu_host]+3*d;
		    for (blk=0; blk<hb; blk++)
		      dg+=dl_start_blk[blk][cpu_host];
		    //
		    lcol[2]=dg;
		    val[2]=(1-x) * y * (1-z);
		    itfnumber++;
		    
		    
		    k=rcvPackAllRev[bi][cpu].intData[3*p+2]+1;j=rcvPackAllRev[bi][cpu].intData[3*p+1];i=rcvPackAllRev[bi][cpu].intData[3*p]; //i,j,k+1
		    cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    //
		    dg=dg_start[cpu_host]+3*d;
		    for (blk=0; blk<hb; blk++)
		      dg+=dl_start_blk[blk][cpu_host];
		    //
		    lcol[3]=dg;
		    val[3]=(1-x) * (1-y) * z;
		    itfnumber++;
		    
		    
		    k=rcvPackAllRev[bi][cpu].intData[3*p+2];j=rcvPackAllRev[bi][cpu].intData[3*p+1]+1;i=rcvPackAllRev[bi][cpu].intData[3*p]+1; //i+1,j+1,k
		    cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    //
		    dg=dg_start[cpu_host]+3*d;
		    for (blk=0; blk<hb; blk++)
		      dg+=dl_start_blk[blk][cpu_host];
		    //
		    lcol[4]=dg;
		    val[4]= x * y * (1-z);
		    itfnumber++;
		    
		    
		    k=rcvPackAllRev[bi][cpu].intData[3*p+2]+1;j=rcvPackAllRev[bi][cpu].intData[3*p+1];i=rcvPackAllRev[bi][cpu].intData[3*p]+1; //i+1,j,k+1
		    cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    //
		    dg=dg_start[cpu_host]+3*d;
		    for (blk=0; blk<hb; blk++)
		      dg+=dl_start_blk[blk][cpu_host];
		    //
		    lcol[5]=dg;
		    val[5]=x * (1-y) * z;
		    itfnumber++;
		    
		    
		    k=rcvPackAllRev[bi][cpu].intData[3*p+2]+1;j=rcvPackAllRev[bi][cpu].intData[3*p+1]+1;i=rcvPackAllRev[bi][cpu].intData[3*p]; //i,j+1,k+1
		    cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    //
		    dg=dg_start[cpu_host]+3*d;
		    for (blk=0; blk<hb; blk++)
		      dg+=dl_start_blk[blk][cpu_host];
		    //
		    lcol[6]=dg;
		    val[6]=(1-x) * y * z;
		    itfnumber++;
			  
		    
		    k=rcvPackAllRev[bi][cpu].intData[3*p+2]+1;j=rcvPackAllRev[bi][cpu].intData[3*p+1]+1;i=rcvPackAllRev[bi][cpu].intData[3*p]+1; //i+1,j+1,k+1
		    cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
		    //
		    dg=dg_start[cpu_host]+3*d;
		    for (blk=0; blk<hb; blk++)
		      dg+=dl_start_blk[blk][cpu_host];
		    //
		    lcol[7]=dg;
		    val[7]=x * y * z;
		    itfnumber++;
		    
		    
		    ///////////////////////////////////////
		    
		    /* for (find=0; find<block_number; find++) */
		    /*   { */
		    /* 	if(bi>find) */
		    /* 	  coeff[find]=1; */
		    /* 	else */
		    /* 	  coeff[find]=0; */
		    /*   } */
		    
		    /////////////////////////////////////////////////////////////////////////////////////////
		    
		    k=intData_index[bi][cpu][3*p+2];j=intData_index[bi][cpu][3*p+1];i=intData_index[bi][cpu][3*p];
		    
		    d=lidxLocal1_matrix(i, j, k, &user[bi],bi);
		    //
		    dg=dg_start[rank]+3*d;
		    for (blk=0; blk<bi; blk++)
		      dg+=dl_start_blk[blk][rank];
		    //
		    
		    for(d_freedom=0;d_freedom<3;d_freedom++){ 
		      row=dg+d_freedom;//lidxLocal(8, 8, 8, &user,0);
		      itfnumber1++;
		    
		      PetscInt             iii;
		      
		      for(iii=0;iii<8;iii++){
			col[iii]=lcol[iii]+d_freedom;
		      }
		      // if(k==33&&j==39)
		      //PetscPrintf(PETSC_COMM_SELF, "row:%d col:%d %d %d %d %d %d %d %d!\n",row,col[0],col[1],col[2],col[3],col[4],col[5],col[6],col[7]);
		      //PetscPrintf(PETSC_COMM_SELF, "itfsearch:%d , col:%d %d %d %d %d %d %d %d -------row:%d!\n",itfc_number,col[0],col[1],col[2],col[3],col[4],col[5],col[6],col[7],row);
		      
		      MatSetValues(Int_matrix,1,&row,8,col,val,INSERT_VALUES);
		      
		    }
		    
		    
		    Check[kk][jj][ii]=20;
		  }
		
		/* 	if(intData_index[bi][cpu][3*p+2]==19) */
		/* 	  PetscPrintf(PETSC_COMM_SELF, "point:%d %d %d  host: %d %d %d!\n",intData_index[bi][cpu][3*p],intData_index[bi][cpu][3*p+1],intData_index[bi][cpu][3*p+2],rcvPackAllRev[bi][cpu].intData[3*p],rcvPackAllRev[bi][cpu].intData[3*p+1],rcvPackAllRev[bi][cpu].intData[3*p+2]); */
	      }
	  }
	}
      }
      
      //PetscPrintf(PETSC_COMM_SELF, "Done matrix!\n");
      
      
    }
    DMDAVecRestoreArray(user[bi].da, user[bi].check, &(Check)); 
    DMDAVecRestoreArray(user[bi].da, user[bi].lInterface, &bcs);
  }


  // if(ass){
  //MatAssemblyBegin(Int_matrix,MAT_FINAL_ASSEMBLY);
  //MatAssemblyEnd(Int_matrix,MAT_FINAL_ASSEMBLY);
    // }

  PetscTime(&v2);
  elapsed_time = v2 - v1;


  // sprintf(filen, "Matfield.dat"); 
/*   PetscViewerASCIIOpen(PETSC_COMM_WORLD, filen, &viewer); */
/*   MatView(Int_matrix, viewer); */
  
  
  //PetscViewerDestroy(&viewer)


  //PetscPrintf(PETSC_COMM_SELF, "Done!\n");
  
  return 0;
}


PetscInt find_cpu(PetscInt i,PetscInt j,PetscInt k,PetscInt hb,PetscInt size,PetscInt **xs,PetscInt **xm,PetscInt **ys,PetscInt **ym,PetscInt **zs,PetscInt **zm){
 
  PetscInt cpu,host=-10;
  
  for (cpu=0;cpu<size;cpu++){
    if(k>=zs[hb][cpu] && k<(zs[hb][cpu]+zm[hb][cpu]))
      if(j>=ys[hb][cpu] && j<(ys[hb][cpu]+ym[hb][cpu]))
	if(i>=xs[hb][cpu] && i<(xs[hb][cpu]+xm[hb][cpu]))
	  host=cpu;
  }
  
  return host;
}





PetscInt lidxLocal_matrix(PetscInt i, PetscInt j, PetscInt k,PetscInt blk,PetscInt cpu,PetscInt **xs,PetscInt **xm,PetscInt **ys,PetscInt **ym,PetscInt **zs,PetscInt **zm)
{
  
  return ((k-zs[blk][cpu]) * (xm[blk][cpu]*ym[blk][cpu]) + (j-ys[blk][cpu])*(xm[blk][cpu]) + (i-xs[blk][cpu]));
}



PetscInt lidxLocal1_matrix(PetscInt i, PetscInt j, PetscInt k, UserCtx *user,PetscInt blk)
{
  DMDALocalInfo	info;
  
  DMDAGetLocalInfo(user->da,&info); 
  
  PetscInt	xs, xe, ys, ye, zs, ze;
  
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
    
  return ((k-zs) * (info.xm*info.ym) + (j-ys)*(info.xm) + (i-xs));
}


PetscErrorCode Interpolation_matrix_rank(UserCtx *user,int *intData_rank_index,double *realData_rank,int *intData_rank,PetscInt bi,PetscInt sb,PetscInt int_num) // for two block
{
   PetscErrorCode       ierr;
  PetscInt             i,j,k;
  char                 filen[80];
  PetscInt             rank,size,blk,level;
  //DM                   packer;
  DMDALocalInfo        info1,info2,info3;
  PetscViewer          viewer;
  PetscReal           ***Check,***bcs;
  
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);


  /////////////////////////////////////----------------------------///////////////////////////////////////////////// reading the coeficients and host
  
  PetscLogDouble v1,v2,elapsed_time;

  PetscTime(&v1);
  
  PetscInt  hb=0;
  
  
  DMDALocalInfo	info = user[bi].info;
  
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  
  
  PetscInt d,dg,cll=0;
  PetscInt itfc_number=0;
  PetscReal lval[8];
  PetscInt	row=0,lcol[8];
  PetscInt col[8];
  PetscReal val[8];
  PetscReal x, y, z;
  PetscInt  itfnumber=0,itfnumber1=0;
  PetscInt d_freedom=0;
  PetscInt coeff[block_number];
  PetscInt p;
  PetscInt cpu_host=0;
  
  
  /**************************************************************************************************************************/
  /* Interpolate the Velocity on the interface nodes
     from the host nodes.
     itfc is the velocity at the cell corners
     hostU is the velocity at the cell centers of the host block */
  /**************************************************************************************************************************/
  
  DMDAVecGetArray(user[bi].da, user[bi].check, &(Check));
  DMDAVecGetArray(user[bi].da, user[bi].lInterface, &bcs);
  for (p=0; p<int_num; p++) {
    
    x = realData_rank[3*p];
    y = realData_rank[3*p+1];
    z = realData_rank[3*p+2];
    
    
    dg=0;d=0;cll=0;
    
    for (i=0;i<8;i++){
      lcol[i]=0;
      col[i]=0;
      lval[i]=0.0;
      val[i]=0.0;
    }
    
    //itfnumber=0;
    PetscInt ii,jj,kk;  
    
    ii=intData_rank_index[3*p];
    jj=intData_rank_index[3*p+1];
    kk=intData_rank_index[3*p+2];

    //////////////////////////////////////each blank just interpolate from related block
    PetscBool blk_block= PETSC_TRUE;
    hb = sb;
    
    if((bcs[kk][jj][ii]<(-10*hb+1.)&& bcs[kk][jj][ii]>(-10*hb-1.))|| (bcs[kk][jj][ii]< (-1000)) ||(bcs[kk][jj][ii]<0 && bcs[kk][jj][ii]>-7.0)){
      blk_block=PETSC_TRUE;
    }else{
      blk_block=PETSC_FALSE;
    }
      
    
    //////////////////////////////////////
    
    if(realData_rank[3*p]>-1.25 && Check[kk][jj][ii]<0 && blk_block)//////////////////////////
      {
	
	hb = sb;
	
	
	k=intData_rank[3*p+2];j=intData_rank[3*p+1];i=intData_rank[3*p]; //i,j,k
	cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	//
	dg=dg_start[cpu_host]+3*d;
	for (blk=0; blk<hb; blk++)
	  dg+=dl_start_blk[blk][cpu_host];
	//
	lcol[0]=dg;
	val[0]=(1-x) * (1-y) * (1-z);
	itfnumber++;
	
	
	k=intData_rank[3*p+2];j=intData_rank[3*p+1];i=intData_rank[3*p]+1; //i+1,j,k
	cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	//
	dg=dg_start[cpu_host]+3*d;
	for (blk=0; blk<hb; blk++)
	  dg+=dl_start_blk[blk][cpu_host];
	//
	lcol[1]=dg;
	val[1]=x * (1-y) * (1-z);
	itfnumber++;
	
	
	k=intData_rank[3*p+2];j=intData_rank[3*p+1]+1;i=intData_rank[3*p]; //i,j+1,k
	cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	//
	dg=dg_start[cpu_host]+3*d;
	for (blk=0; blk<hb; blk++)
	  dg+=dl_start_blk[blk][cpu_host];
	//
	lcol[2]=dg;
	val[2]=(1-x) * y * (1-z);
	itfnumber++;
	
	
	k=intData_rank[3*p+2]+1;j=intData_rank[3*p+1];i=intData_rank[3*p]; //i,j,k+1
	cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	//
	dg=dg_start[cpu_host]+3*d;
	for (blk=0; blk<hb; blk++)
	  dg+=dl_start_blk[blk][cpu_host];
	//
	lcol[3]=dg;
	val[3]=(1-x) * (1-y) * z;
	itfnumber++;
	
	
	k=intData_rank[3*p+2];j=intData_rank[3*p+1]+1;i=intData_rank[3*p]+1; //i+1,j+1,k
	cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	//
	dg=dg_start[cpu_host]+3*d;
	for (blk=0; blk<hb; blk++)
	  dg+=dl_start_blk[blk][cpu_host];
	//
	lcol[4]=dg;
	val[4]= x * y * (1-z);
	itfnumber++;
	
	
	k=intData_rank[3*p+2]+1;j=intData_rank[3*p+1];i=intData_rank[3*p]+1; //i+1,j,k+1
	cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	//
	dg=dg_start[cpu_host]+3*d;
	for (blk=0; blk<hb; blk++)
	  dg+=dl_start_blk[blk][cpu_host];
	//
	lcol[5]=dg;
	val[5]=x * (1-y) * z;
	itfnumber++;
	
	  
	k=intData_rank[3*p+2]+1;j=intData_rank[3*p+1]+1;i=intData_rank[3*p]; //i,j+1,k+1
	cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	//
	dg=dg_start[cpu_host]+3*d;
	for (blk=0; blk<hb; blk++)
	  dg+=dl_start_blk[blk][cpu_host];
	//
	lcol[6]=dg;
	val[6]=(1-x) * y * z;
	itfnumber++;
	
	
	k=intData_rank[3*p+2]+1;j=intData_rank[3*p+1]+1;i=intData_rank[3*p]+1; //i+1,j+1,k+1
	cpu_host=find_cpu(i,j,k,hb,size,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	d=lidxLocal_matrix(i, j, k,hb,cpu_host,cpu_xs,cpu_xm,cpu_ys,cpu_ym,cpu_zs,cpu_zm);
	//
	dg=dg_start[cpu_host]+3*d;
	for (blk=0; blk<hb; blk++)
	  dg+=dl_start_blk[blk][cpu_host];
	//
	lcol[7]=dg;
	val[7]=x * y * z;
	itfnumber++;
	
	
	///////////////////////////////////////
	
	/* for (find=0; find<block_number; find++) */
	/*   { */
	/*     if(bi>find) */
	/*       coeff[find]=1; */
	/*     else */
	/*       coeff[find]=0; */
	/*   } */
	
	/////////////////////////////////////////////////////////////////////////////////////////
	
	k=intData_rank_index[3*p+2];j=intData_rank_index[3*p+1];i=intData_rank_index[3*p];
	if(k>=zs && k<ze)
	  if(j>=ys && j<ye)
	    if(i>=xs && i<xe)
	      {
		
		
		d=lidxLocal1_matrix(i, j, k, &user[bi],bi);
		//
		dg=dg_start[rank]+3*d;
		for (blk=0; blk<bi; blk++)
		  dg+=dl_start_blk[blk][rank];
		//
		//dg=dg_start[rank]+coeff[0]*dl_start_blk[0][rank]+3*d;
		
		for(d_freedom=0;d_freedom<3;d_freedom++){ 
		  row=dg+d_freedom;//lidxLocal(8, 8, 8, &user,0);
		  itfnumber1++;
		  
		  PetscInt             iii;
		  
		  for(iii=0;iii<8;iii++){
		    col[iii]=lcol[iii]+d_freedom;
		  }
		  // if(k==33&&j==39)
		  //PetscPrintf(PETSC_COMM_SELF, "row:%d col:%d %d %d %d %d %d %d %d!\n",row,col[0],col[1],col[2],col[3],col[4],col[5],col[6],col[7]);
		  //PetscPrintf(PETSC_COMM_SELF, "itfsearch:%d , col:%d %d %d %d %d %d %d %d -------row:%d!\n",itfc_number,col[0],col[1],col[2],col[3],col[4],col[5],col[6],col[7],row);
		  
		  MatSetValues(Int_matrix,1,&row,8,col,val,INSERT_VALUES);
		    
		}
		}
	
	Check[kk][jj][ii]=20;
      }
    
  }
  
  //PetscPrintf(PETSC_COMM_SELF, "bi:%d number: %d %d!\n",bi,itfnumber,itfnumber1);
  DMDAVecRestoreArray(user[bi].da, user[bi].check, &(Check));
  DMDAVecRestoreArray(user[bi].da, user[bi].lInterface, &bcs);
  
  
  
  //if(ass){
  //MatAssemblyBegin(Int_matrix,MAT_FINAL_ASSEMBLY);
  //MatAssemblyEnd(Int_matrix,MAT_FINAL_ASSEMBLY);
    //}
  
  PetscTime(&v2);
  elapsed_time = v2 - v1;
  
  
  /*   sprintf(filen, "Matfield.dat"); */
  /*   PetscViewerASCIIOpen(PETSC_COMM_WORLD, filen, &viewer); */
  /*   MatView(Int_matrix, viewer); */
  
  
  //PetscViewerDestroy(&viewer)
  
  
  //PetscPrintf(PETSC_COMM_SELF, "Done!\n");
  
  
  return 0;
}


PetscErrorCode create_matrix (UserCtx *user){

  PetscErrorCode       ierr;
  PetscInt             i,j,k;
  char                 filen[80];
  PetscInt             rank,size,bi,level;
  //DM                   packer;
   PetscViewer          viewer;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  
  ////////////////////
  ierr = DMCompositeCreate(PETSC_COMM_WORLD,&int_packer);CHKERRQ(ierr);
  
  for(bi=0;bi<block_number;bi++){
    ierr = DMCompositeAddDM(int_packer,user[bi].fda);CHKERRQ(ierr);
  }
  
  DMSetUp(int_packer);
  
  ierr = DMCreateGlobalVector(int_packer,&U_int);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(int_packer,&Umult_int);CHKERRQ(ierr);
  ////////////////////
  
  PetscInt N11=0,N1_local=0,NN=0,NN_local=0;
  
  VecGetSize(U_int,&N11);
  NN+=N11;
  VecGetLocalSize(U_int,&N1_local);
  NN_local+=N1_local;

  /////////////////////////////////////////////////////////////////// create and allocate a matrix
  
  MatCreateAIJ(PETSC_COMM_WORLD, NN_local, NN_local, NN, NN,8, PETSC_NULL,8, PETSC_NULL,&Int_matrix);
  MatSetOption(Int_matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  
  ///////////////////////////////////////////////////////////////////

}


PetscErrorCode init_matrix (){
  
  MatZeroEntries(Int_matrix);
  VecZeroEntries(U_int);
  VecZeroEntries(Umult_int);

}



PetscErrorCode da_info (UserCtx *user){

  PetscErrorCode       ierr;
  PetscInt             i,j,k;
  char                 filen[80];
  PetscInt             rank,size,bi,level;
  DMDALocalInfo        info[block_number];
  PetscViewer          viewer;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);

  ////////////////////
  for (i=0; i<block_number; i++)
    DMDAGetLocalInfo(user[i].da,&info[i]);
  
  /////////////////////////////////////////////////////////////// obtaining global indeces(starting index for each cpu)

  PetscMalloc(size*sizeof(PetscInt), &(dg_start));

  ////////////
  dl_start_blk=(int **)calloc(block_number, sizeof(int *));
  for (bi=0; bi<block_number; bi++)
    dl_start_blk[bi] = (int *)calloc(size , sizeof(int));

  ///////////
  cpu_xs=(int **)calloc(block_number, sizeof(int *));
  for (bi=0; bi<block_number; bi++)
    cpu_xs[bi] = (int *)calloc(size , sizeof(int));
  //
  cpu_xm=(int **)calloc(block_number, sizeof(int *));
  for (bi=0; bi<block_number; bi++)
    cpu_xm[bi] = (int *)calloc(size , sizeof(int));
  //
  cpu_ys=(int **)calloc(block_number, sizeof(int *));
  for (bi=0; bi<block_number; bi++)
    cpu_ys[bi] = (int *)calloc(size , sizeof(int));
  //
  cpu_ym=(int **)calloc(block_number, sizeof(int *));
  for (bi=0; bi<block_number; bi++)
    cpu_ym[bi] = (int *)calloc(size , sizeof(int));
  //
  cpu_zs=(int **)calloc(block_number, sizeof(int *));
  for (bi=0; bi<block_number; bi++)
    cpu_zs[bi] = (int *)calloc(size , sizeof(int));
  //
  cpu_zm=(int **)calloc(block_number, sizeof(int *));
  for (bi=0; bi<block_number; bi++)
    cpu_zm[bi] = (int *)calloc(size , sizeof(int));

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////// sending distribution info to all cpus

  ///------------------------///
  for (bi=0; bi<block_number; bi++){
    cpu_xs[bi][rank]=info[bi].xs;cpu_ys[bi][rank]=info[bi].ys;cpu_zs[bi][rank]=info[bi].zs;
    cpu_xm[bi][rank]=info[bi].xm;cpu_ym[bi][rank]=info[bi].ym;cpu_zm[bi][rank]=info[bi].zm;
  }


  PetscBarrier(NULL);
  
  for (bi=0; bi<block_number; bi++)
    for (i=0;i<size;i++){
      MPI_Allreduce(MPI_IN_PLACE,&cpu_xs[bi][i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE,&cpu_ys[bi][i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE,&cpu_zs[bi][i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE,&cpu_xm[bi][i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE,&cpu_ym[bi][i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE,&cpu_zm[bi][i],1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);

    }
  
  PetscBarrier(NULL);
  //PetscPrintf(PETSC_COMM_SELF, "rank:%d, lx: %d %d!\n",rank,cpu_ys[0][0],cpu_ys[0][1]);

  ////////////////////////////////////////////////////////////////////starting index of each cpu and each block
  PetscInt cpu=0;

  for(cpu=0;cpu<size;cpu++){
    dg_start[cpu]=0;
  }
  
  ////
  for (bi=0; bi<block_number; bi++)
    for(cpu=0;cpu<size;cpu++){
      dl_start_blk[bi][cpu]=(3*(cpu_xm[bi][cpu]*cpu_ym[bi][cpu]*cpu_zm[bi][cpu]));
    }
  ////

  ////
  for (i=0; i<size; i++){
    cpu=0;
    while(cpu<(i)){
      for (bi=0; bi<block_number; bi++){
	dg_start[i]+=(3*(cpu_xm[bi][cpu]*cpu_ym[bi][cpu]*cpu_zm[bi][cpu]));
      }
      cpu++;
    }
  }
  

  //PetscBarrier(NULL);

  //////////////////////////////////////////////////////////////////////
  return 0;
}


PetscErrorCode destroy_matrix (){
  MatDestroy(&Int_matrix);
  return 0;
}

PetscErrorCode destroy_Info (){
  PetscInt             rank,size,bi,i;

  for(i = 0; i < block_number; i++){
    free( cpu_xs[i]);
    free( cpu_xm[i]);
    free( cpu_ys[i]);
    free( cpu_ym[i]);
    free( cpu_zs[i]);    free( cpu_zm[i]);
  }
  free(cpu_xs);  free(cpu_xm);  free(cpu_ys);  free(cpu_ym);  free(cpu_zs);  free(cpu_zm);
  //
  PetscFree(dg_start);//PetscFree(dl_start_blk0); PetscFree(dl_start_blk1);

  for(i = 0; i < block_number; i++){
    free( dl_start_blk[i]);
  }
  free(dl_start_blk);

 return 0;
}


PetscErrorCode destroy_Interpolation (){
  
  MatDestroy(&Int_matrix);
  VecDestroy(&U_int);
  VecDestroy(&Umult_int);
  //
  DMDestroy(&int_packer);
 return 0;
}


PetscErrorCode Block_Interface_U(UserCtx *user) {
  PetscInt bi,sb;
  PetscInt ci, cj, ck;
  PetscInt hi, hj, hk, hb;
  PetscInt i, j, k;
  PetscReal x, y, z;
  Vec	hostU;
  Cmpnts ***itfc,***ubcs;
  PetscReal *hostu,***nvert;
  
  VecScatter tolocalall;
  PetscErrorCode ierr;
  Vec Ub[block_number];
  PetscLogDouble v1,v2,elapsed_time;
  
  PetscTime(&v1);
  
  /////////////////////////////////////////////////////////////////////////////
  
  DMCompositeGetAccess(int_packer,U_int,&Ub[0],&Ub[1],&Ub[2]);
  
  for (bi=0; bi<block_number; bi++) {
    VecCopy(user[bi].Ucat, Ub[bi]);
  }
  
  DMCompositeRestoreAccess(int_packer,U_int,&Ub[0],&Ub[1],&Ub[2]);

  //////////////////////

  MatMult(Int_matrix,U_int,Umult_int);
  
  //////////////////////
  PetscBarrier(PETSC_NULL);

  DMCompositeGetAccess(int_packer,Umult_int,&Ub[0],&Ub[1],&Ub[2]);
  
  for (bi=0; bi<block_number; bi++) {
    VecCopy(Ub[bi],user[bi].Itfc);
  }
  
  DMCompositeRestoreAccess(int_packer,Umult_int,&Ub[0],&Ub[1],&Ub[2]);
  
  for (bi=0; bi<block_number; bi++) {
  
  DMGlobalToLocalBegin(user[bi].fda, user[bi].Itfc, INSERT_VALUES,user[bi].lItfc);
  DMGlobalToLocalEnd(user[bi].fda, user[bi].Itfc, INSERT_VALUES,user[bi].lItfc);
  }

  for (bi=0; bi<block_number; bi++) {
    ierr = VecDestroy(&Ub[bi]);CHKERRQ(ierr);
  }
  
  PetscBarrier(PETSC_NULL);
  ///////////////////////////////////////////////
  
  if(rotateframe){
    for (bi=0; bi<block_number; bi++) {
      Interface_Rotation_matrix(&(user[bi]),bi);
    }
  }

  PetscBarrier(PETSC_NULL);
  ///////-------------------------------------------------------------  interface timing
  PetscTime(&v2);
  elapsed_time = v2 - v1;
  PetscPrintf(PETSC_COMM_WORLD,"elapsed_time:%le \n", elapsed_time);
  ///////-------------------------------------------------------------


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
   if (ys==0 && bi==1){ 
	j=0;
     for (k=zs; k<lze; k++) {
	for (i=xs; i<lxe; i++) {

	itfc[k][j][i].x=0.0; 
	itfc[k][j][i].y=0.0;
	itfc[k][j][i].z=0.0;
	 }
	}
    }  

    if (user[bi].bctype[0] == 0 && xs==0) {
      i=0;//1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	 
 
	  ubcs[k][j][i].x = 0.25 * (itfc[k-1][j  ][i  ].x +
				    itfc[k  ][j-1][i  ].x +
				    itfc[k  ][j  ][i  ].x +
				    itfc[k-1][j-1][i  ].x);
	  ubcs[k][j][i].y = 0.25 * (itfc[k-1][j  ][i  ].y +
				    itfc[k  ][j-1][i  ].y +
				    itfc[k  ][j  ][i  ].y +
				    itfc[k-1][j-1][i  ].y);
	  ubcs[k][j][i].z = 0.25 * (itfc[k-1][j  ][i  ].z +
				    itfc[k  ][j-1][i  ].z +
				    itfc[k  ][j  ][i  ].z +
				    itfc[k-1][j-1][i  ].z);



	  if (bi==1 && j==1){


	    ubcs[k][j][i].x = 0.5 * (itfc[k-1][j  ][i  ].x +				    
				     itfc[k  ][j  ][i  ].x); 
				    
	    ubcs[k][j][i].y = 0.5 * (itfc[k-1][j  ][i  ].y +
				     itfc[k  ][j  ][i  ].y );

	    ubcs[k][j][i].z = 0.5 * (itfc[k-1][j  ][i  ].z +
				     itfc[k  ][j  ][i  ].z);


	  }
	

	  
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



	  if (bi==1 && j==1){


	    ubcs[k][j][i].x = 0.5 * (itfc[k-1][j  ][i  ].x +				    
				     itfc[k  ][j  ][i  ].x); 
				    
	    ubcs[k][j][i].y = 0.5 * (itfc[k-1][j  ][i  ].y +
				     itfc[k  ][j  ][i  ].y );

	    ubcs[k][j][i].z = 0.5 * (itfc[k-1][j  ][i  ].z +
				     itfc[k  ][j  ][i  ].z);


	  }
	


	  
	  ucont[k][j][i].x = (ubcs[k][j][i+1].x  *
			      icsi[k][j][i].x +
			      ubcs[k][j][i+1].y  *
			      icsi[k][j][i].y +
			      ubcs[k][j][i+1].z  *
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
	  ucont[k][j][i].y = ( ubcs[k][j][i].x *
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
	  if (bi==1 && j==1){

	  ubcs[k][j][i].x = 0.5 * (itfc[k  ][j  ][i  ].x +
				    itfc[k  ][j  ][i-1].x);// +
				   // itfc[k  ][j-1][i  ].x +
				   // itfc[k  ][j-1][i-1].x);
	  ubcs[k][j][i].y = 0.5 * (itfc[k  ][j  ][i  ].y +
				    itfc[k  ][j  ][i-1].y );//+
				 //   itfc[k  ][j-1][i  ].y +
				  //  itfc[k  ][j-1][i-1].y);
	  ubcs[k][j][i].z = 0.5 * (itfc[k  ][j  ][i  ].z +
				    itfc[k  ][j  ][i-1].z); 
				 //   itfc[k  ][j-1][i  ].z +
				  //  itfc[k  ][j-1][i-1].z);

		
	  ucont[k][j][i].z = ( ubcs[k][j][i].x*
			       kzet[k][j][i].x +
			       ubcs[k][j][i].y*
			       kzet[k][j][i].y +
			       ubcs[k][j][i].z *
			       kzet[k][j][i].z);
	
	}
	  
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
	  if(bi==1 && j==1){
	    ubcs[k+1][j][i].x = 0.5 * (itfc[k  ][j  ][i  ].x +
				     itfc[k  ][j  ][i-1].x) ;
				//      itfc[k  ][j-1][i  ].x +
	    //  itfc[k  ][j-1][i-1].x);
	    ubcs[k+1][j][i].y = 0.5 * (itfc[k  ][j  ][i  ].y +
				      itfc[k  ][j  ][i-1].y); 
				   //   itfc[k  ][j-1][i  ].y +
				    //  itfc[k  ][j-1][i-1].y);
	    ubcs[k+1][j][i].z = 0.5 * (itfc[k  ][j  ][i  ].z +
				      itfc[k  ][j  ][i-1].z);



	  //	    if (ubcs[k+1][j][i].z<.9){
	  //	    PetscPrintf(PETSC_COMM_SELF, "i%d  j%d  k%d    ubcs.z   %le   itfc[i].z     %le   itfc[i-1].z    %le !\n",i,j,k,ubcs[k+1][j][i].z,itfc[k][j][i].z,itfc[k][j][i-1].z);   
	  //	  }
				   //   itfc[k  ][j-1][i  ].z +
				   //   itfc[k  ][j-1][i-1].z);
	  ucont[k][j][i].z = ( ubcs[k+1][j][i].x*
			       kzet[k  ][j][i].x +
			       ubcs[k+1][j][i].y*
			       kzet[k  ][j][i].y +
			       ubcs[k+1][j][i].z *
			       kzet[k  ][j][i].z);

 }
	}
      }
    }

    DMDAVecRestoreArray(user[bi].fda, user[bi].Ucont, &ucont);

    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

      //This part is for blanking
    if(blank){
      PetscInt	llxs, llxe, llys, llye, llzs, llze;
      llxs = xs-1; llxe = xe+1;
      llys = ys-1; llye = ye+1;
      llzs = zs-1; llze = ze+1;
      
      if (xs==0) llxs = xs+1;
      if (ys==0) llys = ys+1;
      if (zs==0) llzs = zs+1;
      
      if (xe==mx) llxe = xe-1;
      if (ye==my) llye = ye-1;
      if (ze==mz) llze = ze-1;
      
      DMDAVecGetArray(user[bi].da, user[bi].lNvert, &nvert);
      DMDAVecGetArray(user[bi].fda, user[bi].Ucont, &ucont);
      
      for (sb=0; sb<NumberOfBlank; sb++) {
	if(blank_id[bi][sb]){
	  PetscInt ip, im, jp, jm, kp, km;
	  PetscInt ii, jj, kk;
	  
	  for (k=llzs; k<llze; k++) {
	    for (j=llys; j<llye; j++) {
	      for (i=llxs; i<llxe; i++) {
		
		ip = (i<mx-2?(i+1):(i));
	      im = (i>1   ?(i-1):(i));
	      
	      jp = (j<my-2?(j+1):(j));
	      jm = (j>1   ?(j-1):(j));
	      
	      kp = (k<mz-2?(k+1):(k));
	      km = (k>1   ?(k-1):(k));

	      if (((int)(nvert[k][j][i]+0.1) < (blank_id[bi][sb]-1)*100+10*(sb+1)) && /// for all blank nverts are 10 or 1010
		  ((int)(nvert[k][j][i]+0.1) > (blank_id[bi][sb]-1)*100+10*(sb+1)-3)) {
		// flux in x direction
		kk=k;	jj=j;
		for (ii=im; ii<ip; ii++) {
		  if(kk<lze && kk>=zs){
		    if(jj<lye && jj>=ys){
		      if(ii<lxe && ii>=xs){
			ubcs[kk][jj][ii].x = 0.25 * (itfc[kk-1][jj  ][ii  ].x +
						     itfc[kk  ][jj-1][ii  ].x +
						     itfc[kk  ][jj  ][ii  ].x +
						     itfc[kk-1][jj-1][ii  ].x) ;
			ubcs[kk][jj][ii].y = 0.25 * (itfc[kk-1][jj  ][ii  ].y +
						     itfc[kk  ][jj-1][ii  ].y +
						     itfc[kk  ][jj  ][ii  ].y +
						     itfc[kk-1][jj-1][ii  ].y);
			ubcs[kk][jj][ii].z = 0.25 * (itfc[kk-1][jj  ][ii  ].z +
						     itfc[kk  ][jj-1][ii  ].z +
						     itfc[kk  ][jj  ][ii  ].z +
						     itfc[kk-1][jj-1][ii  ].z);



			if (bi==0 && (jj==1 || jj==2)) {		  
			  
			  ubcs[kk][jj][ii].x = 0.5 * (itfc[kk  ][jj+1  ][ii  ].x +
						    itfc[kk-1  ][jj+1  ][ii].x);
			}

			
			
			ucont[kk][jj][ii].x = (ubcs[kk][jj][ii].x*
					       icsi[kk][jj][ii].x +
					       ubcs[kk][jj][ii].y*
					       icsi[kk][jj][ii].y +
					       ubcs[kk][jj][ii].z*
					       icsi[kk][jj][ii].z);
		      }
		    }
		  }
		}
		// flux in y direction
		kk=k;   ii=i;
		for (jj=jm; jj<jp; jj++) {
		  if(kk<lze && kk>=zs){
		    if(jj<lye && jj>=ys){
		      if(ii<lxe && ii>=xs){
			ubcs[kk][jj][ii].x = 0.25 * (itfc[kk-1][jj  ][ii  ].x +
						     itfc[kk  ][jj  ][ii-1].x +
						     itfc[kk  ][jj  ][ii  ].x +
						     itfc[kk-1][jj  ][ii-1].x);
			ubcs[kk][jj][ii].y = 0.25 * (itfc[kk-1][jj  ][ii  ].y +
						     itfc[kk  ][jj  ][ii-1].y +
						     itfc[kk  ][jj  ][ii  ].y +
						     itfc[kk-1][jj  ][ii-1].y);
			ubcs[kk][jj][ii].z = 0.25 * (itfc[kk-1][jj  ][ii  ].z +
						     itfc[kk  ][jj  ][ii-1].z +
						     itfc[kk  ][jj  ][ii  ].z +
						     itfc[kk-1][jj  ][ii-1].z);
			ucont[kk][jj][ii].y =  ( ubcs[kk][jj][ii].x*
						 jeta[kk][jj][ii].x +
						 ubcs[kk][jj][ii].y *
						 jeta[kk][jj][ii].y +
						 ubcs[kk][jj][ii].z *
						 jeta[kk][jj][ii].z);
			if (bi==0 && jj==0){
			  
			  ubcs[kk][jj][ii].y=0.0;
			  
			}
		      }
		    }
		  }
		}
		// flux in z direction
		jj=j;  ii=i;
		for (kk=km; kk<kp; kk++) {
		  if(kk<lze && kk>=zs){
		    if(jj<lye && jj>=ys){
		      if(ii<lxe && ii>=xs){
			ubcs[kk][jj][ii].x = 0.25 * (itfc[kk  ][jj  ][ii  ].x +
						  itfc[kk  ][jj  ][ii-1].x +
						  itfc[kk  ][jj-1][ii  ].x +
						  itfc[kk  ][jj-1][ii-1].x);
			ubcs[kk][jj][ii].y = 0.25 * (itfc[kk  ][jj  ][ii  ].y +
						  itfc[kk  ][jj  ][ii-1].y +
						  itfc[kk  ][jj-1][ii  ].y +
						  itfc[kk  ][jj-1][ii-1].y);
			ubcs[kk][jj][ii].z = 0.25 * (itfc[kk  ][jj  ][ii  ].z +
						  itfc[kk  ][jj  ][ii-1].z +
						  itfc[kk  ][jj-1][ii  ].z +
						  itfc[kk  ][jj-1][ii-1].z);

			
			if (bi==0 && (jj==1 || jj==2)) {		  
			  //	  PetscPrintf(PETSC_COMM_SELF, "i%d  j%d  k%d    ubcs.z   %le   itfc[i].z     %le   itfc[i-1].z    %le !\n",ii,jj,kk,ubcs[kk][jj][ii].z,itfc[kk][jj][ii].z,itfc[kk][jj+1][ii].z);   
			
			  ubcs[kk][jj][ii].z = 0.5 * (itfc[kk  ][jj+1  ][ii  ].z +
						    itfc[kk  ][jj+1  ][ii-1].z);
			  // ubcs[kk][jj][ii].z=1.0;			  
			}
			  
			  ucont[kk][jj][ii].z =  ( ubcs[kk][jj][ii].x*
						   kzet[kk][jj][ii].x +
						   ubcs[kk][jj][ii].y*
						   kzet[kk][jj][ii].y +
						   ubcs[kk][jj][ii].z *
						   kzet[kk][jj][ii].z);
			}
		      }
		    }
		  }
		}
	      }// if (nvert)
	    }
	  }
	  }
      }
     //for sb
      DMDAVecRestoreArray(user[bi].da, user[bi].lNvert, &nvert);
      DMDAVecRestoreArray(user[bi].fda, user[bi].Ucont, &ucont);
      //      PetscPrintf(PETSC_COMM_WORLD, "Local to global lUcont _U ");

      //DMLocalToGlobalBegin(user[bi].fda, user[bi].lUcont,INSERT_VALUES,user[bi].Ucont);
      //DMLocalToGlobalEnd(user[bi].fda, user[bi].lUcont,INSERT_VALUES,user[bi].Ucont);

      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);

      PetscPrintf(PETSC_COMM_WORLD, "blank flux done!\n");   

    } // if blank
   
  

    DMDAVecRestoreArray(user[bi].fda, user[bi].lItfc, &itfc);
    DMDAVecRestoreArray(user[bi].fda, user[bi].Bcs.Ubcs, &ubcs);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lZet, &kzet);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lEta, &jeta);
    DMDAVecRestoreArray(user[bi].fda, user[bi].lCsi, &icsi);
 
  } // bi



  for (bi=0; bi<block_number; bi++) {
    Contra2Cart(&(user[bi]));
    GhostNodeVelocity(&user[bi]);
  }

  return(0);      
}


PetscErrorCode Interface_Rotation_matrix(UserCtx *user,PetscInt bi) {

  PetscReal *hostu,***nvert,***inter;
  Cmpnts ***itfc,***litfc;
  //PetscReal theta_x,theta_y,theta_z;
  PetscInt i, j, k,sb;
    

/*   PetscReal t; */
/*   t=ti*user->dt; */
/*   theta_x=(PI/1.)*t; */
  /////////////////////////////////////////////////////////////////// Rotation 
    
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  
  lxs = xs-1; lxe = xe+1;
  lys = ys-1; lye = ye+1;
  lzs = zs-1; lze = ze+1;
  
  if (xs==0) lxs = xs;
  if (ys==0) lys = ys;
  if (zs==0) lzs = zs;
  
  if (xe==mx) lxe = xe;
  if (ye==my) lye = ye;
  if (ze==mz) lze = ze;
  
  DMDAVecGetArray(user->fda, user->Itfc, &itfc);
  DMDAVecGetArray(user->fda, user->lItfc, &litfc);
  /**************************************************************************************************************************/
  /* Create boundary condition for flux (ucont) and cell surface
     center velocity (ubcs) from the interpolated velocity (itfc).
     
     itfc is the velocity at the cell surface centers */
  /**************************************************************************************************************************/
  if(bi==2){
    if (user->bctype[0] == 0 && xs==0) {
      i=0;//1;
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  ///---------------------------------------------------------------------------///x-rotation
	  itfc[k][j][i].x = litfc[k][j][i].x*cos(-theta_x)-litfc[k][j][i].y*sin(-theta_x);
	  
	  itfc[k][j][i].y = litfc[k][j][i].x*sin(-theta_x)+litfc[k][j][i].y*cos(-theta_x);
	  
	  itfc[k][j][i].z = litfc[k][j][i].z;

	}
      }
    }
    
    if (user->bctype[1] == 0 && xe==mx) {
      i=mx-2;
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  ///---------------------------------------------------------------------------///x-rotation
	  itfc[k][j][i].x = litfc[k][j][i].x*cos(-theta_x)-litfc[k][j][i].y*sin(-theta_x);
	  
	  itfc[k][j][i].y = litfc[k][j][i].x*sin(-theta_x)+litfc[k][j][i].y*cos(-theta_x);
	  
	  itfc[k][j][i].z = litfc[k][j][i].z;
	  
	}
      }
    }
    
    if (user->bctype[2] == 0 && ys==0) {
      j=0;//1;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  ///---------------------------------------------------------------------------///x-rotation
	  itfc[k][j][i].x = litfc[k][j][i].x*cos(-theta_x)-litfc[k][j][i].y*sin(-theta_x);
	  
	  itfc[k][j][i].y = litfc[k][j][i].x*sin(-theta_x)+litfc[k][j][i].y*cos(-theta_x);
	  
	  itfc[k][j][i].z = litfc[k][j][i].z;
	  
	}
      }
    }
    
    if (user->bctype[3] == 0 && ye==my) {
      j=my-2;
      for (k=zs; k<ze; k++) {
	for (i=xs; i<xe; i++) {
	  ///---------------------------------------------------------------------------///x-rotation
	  itfc[k][j][i].x = litfc[k][j][i].x*cos(-theta_x)-litfc[k][j][i].y*sin(-theta_x);
	  
	  itfc[k][j][i].y = litfc[k][j][i].x*sin(-theta_x)+litfc[k][j][i].y*cos(-theta_x);
	  
	  itfc[k][j][i].z = litfc[k][j][i].z;
	  
	}
      }
    }
    
    
    if (user->bctype[4] == 0 && zs==0) {
      k=0;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  ///---------------------------------------------------------------------------///x-rotation
	  itfc[k][j][i].x = litfc[k][j][i].x*cos(-theta_x)-litfc[k][j][i].y*sin(-theta_x);
	  
	  itfc[k][j][i].y = litfc[k][j][i].x*sin(-theta_x)+litfc[k][j][i].y*cos(-theta_x);
	  
	  itfc[k][j][i].z = litfc[k][j][i].z;
	  
	}
      }
    }
    
    if (user->bctype[5] == 0 && ze==mz) {
      k=mz-2;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  ///---------------------------------------------------------------------------///x-rotation
	  itfc[k][j][i].x = litfc[k][j][i].x*cos(-theta_x)-litfc[k][j][i].y*sin(-theta_x);
	  
	  itfc[k][j][i].y = litfc[k][j][i].x*sin(-theta_x)+litfc[k][j][i].y*cos(-theta_x);
	  
	  itfc[k][j][i].z = litfc[k][j][i].z;
	  
	}
      }
    }
    
  }
  
  //This part is for blanking
  if(blank && bi==1){
    
    DMDAVecGetArray(user->da, user->Interface, &inter);
    
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  
	  
	  if (inter[k][j][i] < -1.0 ) { //-1 to identify blank interface (-10)
	    
	    itfc[k][j][i].x = litfc[k][j][i].x*cos(theta_x_p)-litfc[k][j][i].y*sin(theta_x_p);
	    
	    itfc[k][j][i].y = litfc[k][j][i].x*sin(theta_x_p)+litfc[k][j][i].y*cos(theta_x_p);
	    
	    itfc[k][j][i].z = litfc[k][j][i].z;
	  }
	  
	  
	}
      }
    }
    //for sb
    DMDAVecRestoreArray(user->da, user->Interface, &inter);  
  } // if blank
  
  
  DMDAVecRestoreArray(user->fda, user->Itfc, &itfc);
  DMDAVecRestoreArray(user->fda, user->lItfc, &litfc);

  ///////////////////////////////////////////////////////////////////
  DMGlobalToLocalBegin(user->fda, user->Itfc, INSERT_VALUES,user->lItfc);
  DMGlobalToLocalEnd(user->fda, user->Itfc, INSERT_VALUES,user->lItfc);



return(0); 
}



/* PetscErrorCode Interpolation_mult(UserCtx *user) // for two block */
/* { */
/*   PetscLogDouble v1,v2,elapsed_time; */
/*   PetscInt bi; */
/*   Vec Ub[block_number]; */
  
  
/*   DMCompositeGetAccess(int_packer,U_int,&Ub[0],&Ub[1]); */
  
/*   for (bi=0; bi<block_number; bi++) { */
/*     VecCopy(user[bi].Q, Ub[bi]); */
/*   } */
  
/*   DMCompositeRestoreAccess(int_packer,U_int,&Ub[0],&Ub[1]); */

/*   ////////////////////// */

/*   MatMult(Int_matrix,U_int,Umult_int); */
  
/*   ////////////////////// */
/*   PetscBarrier(PETSC_NULL); */

/*   DMCompositeGetAccess(int_packer,Umult_int,&Ub[0],&Ub[1]); */
  
/*   for (bi=0; bi<block_number; bi++) { */
/*     VecCopy(Ub[bi],user[bi].Q); */
/*   } */
  
/*   DMCompositeRestoreAccess(int_packer,Umult_int,&Ub[0],&Ub[1]); */
  

/*   for (bi=0; bi<block_number; bi++) { */
/*     VecDestroy(&Ub[bi]); */
/*   } */
   
/*  return 0;  */
/* } */



PetscErrorCode mat_ass (){
  
  MatAssemblyBegin(Int_matrix,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Int_matrix,MAT_FINAL_ASSEMBLY);

  //////////////////////////
  /* char                 filen[80]; */
  /* PetscViewer          viewer; */

  /* sprintf(filen, "Matfield.dat"); */
  /* PetscViewerASCIIOpen(PETSC_COMM_WORLD, filen, &viewer); */
  /* MatView(Int_matrix, viewer); */
  /* PetscViewerDestroy(&viewer); */
  
  return 0;
}


////////////////////////////////////////////////////////////////////////// blank

PetscErrorCode blank_randomdirection(Cmpnts p, PetscInt ip, PetscInt jp,
			       PetscReal xbp_min, PetscReal ybp_min,
			       PetscReal zbp_max, PetscReal dcx, PetscReal dcy,
			       PetscReal dir[3],PetscInt seed)
{
  Cmpnts endpoint;
  PetscReal s;

  PetscReal xpc, ypc; 
  
  xpc = dcx * (ip+0.5) + xbp_min;
  ypc = dcy * (jp+0.5) + ybp_min;
    
  // init rand()
  //  srand(time(NULL)+seed);
  srand(seed);
  // Generate a random number [-0.5, 0.5)
  s = rand() / ((double)RAND_MAX + 1) - 0.5;
  endpoint.x = xpc + s * dcx;
  endpoint.y = ypc + s * dcy;
  endpoint.z = zbp_max + 0.2;

  dir[0] = endpoint.x - p.x;
  dir[1] = endpoint.y - p.y;
  dir[2] = endpoint.z - p.z;

  s = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
  dir[0] /= s;
  dir[1] /= s;
  dir[2] /= s;
  return 0;
}

PetscInt blank_point_cell_advanced(Cmpnts p, PetscInt ip, PetscInt jp, PetscInt kp,
			     IBMNodesblank *ibm, PetscInt ncx, PetscInt ncy,
			     PetscInt ncz, PetscReal dcx, PetscReal dcy,
			     PetscReal xbp_min, PetscReal ybp_min,
			     PetscReal zbp_max, List *cell_trg,
			     PetscInt flg)
{
  PetscInt	*nv1 = ibm->nv1, *nv2 = ibm->nv2, *nv3 = ibm->nv3;
  //  PetscReal	*nf_x = ibm->nf_x, *nf_y = ibm->nf_y, *nf_z = ibm->nf_z;
  PetscReal	*x_bp = ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp;

  PetscInt	i, j, k, ln_v, n1e, n2e, n3e, nintp;
  
  PetscInt	nvert_l;
  PetscReal	dt[1000], ndotn=0., dirdotn;
  //  Cmpnts        dnn[1000],nn;

  PetscReal	epsilon = 1.e-8;
  PetscReal     eps_tangent=1.e-10;

  PetscBool	*Element_Searched;
  j = jp; i = ip;

  PetscBool NotDecided = PETSC_TRUE, Singularity = PETSC_FALSE;
  PetscReal t, u, v;
  PetscReal orig[3], dir[3], vert0[3], vert1[3], vert2[3];

  node *current;
  PetscInt searchtimes=0;
  PetscMalloc(ibm->n_elmt*sizeof(PetscBool), &Element_Searched);
  if (flg) 
    PetscPrintf(PETSC_COMM_SELF, " serch itr\n");
  
  while (NotDecided) {
    
    searchtimes++;
    nintp = 0 ;
    blank_randomdirection(p, ip, jp, xbp_min, ybp_min, zbp_max, dcx, dcy, dir, searchtimes);
    Singularity = PETSC_FALSE;
    if (flg) 
      PetscPrintf(PETSC_COMM_SELF, " serch itr, dir %d %le %le %le\n", searchtimes,dir[0],dir[1],dir[2]);
    
    for (ln_v=0; ln_v<ibm->n_elmt; ln_v++) {
      Element_Searched[ln_v] = PETSC_FALSE;
    }
    
    for (k=kp; k<ncz; k++) {
      current = cell_trg[k*ncx*ncy+j*ncx+i].head;
      while (current) {
	ln_v = current->Node;
	if (!Element_Searched[ln_v]) {
	  Element_Searched[ln_v] = PETSC_TRUE;
	  n1e = nv1[ln_v]; n2e = nv2[ln_v]; n3e = nv3[ln_v];
	 
	  orig[0] = p.x; orig[1] = p.y, orig[2] = p.z;

	  vert0[0] = x_bp[n1e]; vert0[1] = y_bp[n1e]; vert0[2] = z_bp[n1e];
	  vert1[0] = x_bp[n2e]; vert1[1] = y_bp[n2e]; vert1[2] = z_bp[n2e];
	  vert2[0] = x_bp[n3e]; vert2[1] = y_bp[n3e]; vert2[2] = z_bp[n3e];
            
	  //	  dirdotn=dir[0]*nn.x+dir[1]*nn.y+dir[2]*nn.z;

	  nvert_l = intsect_triangle(orig, dir, vert0, vert1, vert2, &t, &u, &v);
	  
	  if (flg) 
	    PetscPrintf(PETSC_COMM_SELF, "elm, %d %d %le %le %le %d %d %d %le\n",ln_v,nvert_l,t,u,v,n1e,n2e,n3e,dirdotn);
	  
	  if (nvert_l > 0 && t>0) {
	    dt[nintp] = t;
	    //    dnn[nintp].x=nn.x;dnn[nintp].y=nn.y;dnn[nintp].z=nn.z;

	    nintp ++;
	    PetscInt temp;
	    for (temp = 0; temp < nintp-1; temp++) {
	      // Two interception points are the same, this leads to huge
	      // trouble for crossing number test
	      // Rather to program for all cases, we use a new line to
	      // repeat the test
	      //	      ndotn=dnn[temp].x*nn.x+dnn[temp].y*nn.y+dnn[temp].z*nn.z;	      
	      
	      if ((fabs(t-dt[temp]) < epsilon && ndotn>-0.97)){ 
		
		Singularity = PETSC_TRUE;
	      }
	    }
	    if (Singularity) break;
	  }
	}
	if (Singularity) {
	  break;
	}
	else {
	  current = current->next;
	}
	} // Search through the list
      if (Singularity) {
	break;
      }
    } // for k
    if (flg) 
      PetscPrintf(PETSC_COMM_SELF, " serch itr, %d %le \n",nintp,dirdotn);
    
    if (!Singularity) {
      NotDecided = PETSC_TRUE;
      if (nintp%2) { // The interception point number is odd, inside body
	PetscFree(Element_Searched);
	return 4;
      }
      else {
	PetscFree(Element_Searched);
	return 0;
      }
    }
  }
  PetscFree(Element_Searched);
  return 0;
}

PetscErrorCode blank_change_nvert(UserCtx *user){// change all the blank nvert to 10//10-1//10-2 for the cases with collision
  
  DM	da = user->da, fda = user->fda;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscReal	***nvert,***lnvert;
  PetscInt      sb=0;
  PetscInt	i, j, k;
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs;
  if (ys==0) lys = ys;
  if (zs==0) lzs = zs;

  if (xe==mx) lxe = xe;
  if (ye==my) lye = ye;
  if (ze==mz) lze = ze;

  DMDAVecGetArray(da, user->Nvert, &nvert);

  for(sb=0;sb<NumberOfBlank;sb++){

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if(blank_id[user->_this][sb]==1){
	    if ((int)(nvert[k][j][i]+0.5) == 10*(sb+1)) nvert[k][j][i]=10.;
	    if ((int)(nvert[k][j][i]+0.5) == 10*(sb+1)-1) nvert[k][j][i]=10-1.;
	    if ((int)(nvert[k][j][i]+0.5) == 10*(sb+1)-2) nvert[k][j][i]=10-2.;
	  }else if(blank_id[user->_this][sb]==11){ //for moving immersed boundaries 
	    if ((int)(nvert[k][j][i]+0.5) == (blank_id[user->_this][sb]-1)*100+10*(sb+1)) nvert[k][j][i]=1010.0;
	    if ((int)(nvert[k][j][i]+0.5) == (blank_id[user->_this][sb]-1)*100+10*(sb+1)-1) nvert[k][j][i]=1010-1.;
	    if ((int)(nvert[k][j][i]+0.5) == (blank_id[user->_this][sb]-1)*100+10*(sb+1)-2) nvert[k][j][i]=1010-2.;
	  }
	}
      }
    }
  }
  
  DMDAVecRestoreArray(da, user->Nvert, &nvert);
  //
  DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
  DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
  
  return 0;
}


PetscErrorCode blank_reset_nvert(UserCtx *user){// reset all the nverts before blank search
  
  DM	da = user->da, fda = user->fda;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscReal	***nvert,***lnvert;
  PetscInt	i, j, k;
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs;
  if (ys==0) lys = ys;
  if (zs==0) lzs = zs;

  if (xe==mx) lxe = xe;
  if (ye==my) lye = ye;
  if (ze==mz) lze = ze;

  DMDAVecGetArray(da, user->Nvert, &nvert);
  
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i]>(7.0)){
	  nvert[k][j][i] = 0.0;
	}
      }
    }
  }
  
  
  DMDAVecRestoreArray(da, user->Nvert, &nvert);
  //
  DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
  DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);

  return 0;
}



PetscErrorCode blank_search_advanced(UserCtx *user, IBMNodesblank *ibm, 
				     PetscInt ibi)

/*      Note : Always go from ibi (immersed body number) 0 -> NumberOfBodies  */
/*             Nvert should be set to zero before any new search and this */
/*             happens if ibi==0--not anymore! set nvert=0 manually!*/
{
  DM	da = user->da, fda = user->fda;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt	ncx = 10, ncy = 10, ncz = 10;
  List          *cell_trg;
  PetscReal	xbp_min, ybp_min, zbp_min, xbp_max, ybp_max, zbp_max;
  PetscReal	*x_bp = ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp;

  PetscInt 	ln_v, n_v = ibm->n_v;

  PetscInt	i, j, k;

  PetscReal	dcx, dcy, dcz;
  PetscInt	n1e, n2e, n3e;
  PetscReal	xv_min, yv_min, zv_min, xv_max, yv_max, zv_max;
  PetscInt	iv_min, iv_max, jv_min, jv_max, kv_min, kv_max;
  PetscReal	***nvert,***lnvert;
  PetscInt	ic, jc, kc;
  PetscInt      sb=0;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs;
  if (ys==0) lys = ys;
  if (zs==0) lzs = zs;

  if (xe==mx) lxe = xe;
  if (ye==my) lye = ye;
  if (ze==mz) lze = ze;

  xbp_min = 1.e23; xbp_max = -1.e23;
  ybp_min = 1.e23; ybp_max = -1.e23;
  zbp_min = 1.e23; zbp_max = -1.e23;

  for(i=0; i<n_v; i++) {
    
    xbp_min = PetscMin(xbp_min, x_bp[i]);
    xbp_max = PetscMax(xbp_max, x_bp[i]);
    
    ybp_min = PetscMin(ybp_min, y_bp[i]);
    ybp_max = PetscMax(ybp_max, y_bp[i]);

    zbp_min = PetscMin(zbp_min, z_bp[i]);
    zbp_max = PetscMax(zbp_max, z_bp[i]);
  }

  xbp_min -= 0.01; xbp_max += 0.01;
  ybp_min -= 0.01; ybp_max += 0.01;
  zbp_min -= 0.01; zbp_max += 0.01;
  
  dcx = (xbp_max - xbp_min) / (ncx - 1.);
  dcy = (ybp_max - ybp_min) / (ncy - 1.);
  dcz = (zbp_max - zbp_min) / (ncz - 1.);

  PetscMalloc(ncz * ncy * ncx * sizeof(List), &cell_trg);
 
  for (k=0; k<ncz; k++) {
    for (j=0; j<ncy; j++) {
      for (i=0; i<ncx; i++) {
	initlist(&cell_trg[k*ncx*ncy + j*ncx + i]);
      }
    }
  }
    
  for (ln_v=0; ln_v < ibm->n_elmt; ln_v++) {

    n1e = ibm->nv1[ln_v]; n2e = ibm->nv2[ln_v]; n3e = ibm->nv3[ln_v];
    
    xv_min = PetscMin(PetscMin(x_bp[n1e], x_bp[n2e]), x_bp[n3e]);
    xv_max = PetscMax(PetscMax(x_bp[n1e], x_bp[n2e]), x_bp[n3e]);

    yv_min = PetscMin(PetscMin(y_bp[n1e], y_bp[n2e]), y_bp[n3e]);
    yv_max = PetscMax(PetscMax(y_bp[n1e], y_bp[n2e]), y_bp[n3e]);
    
    zv_min = PetscMin(PetscMin(z_bp[n1e], z_bp[n2e]), z_bp[n3e]);
    zv_max = PetscMax(PetscMax(z_bp[n1e], z_bp[n2e]), z_bp[n3e]);
    
    iv_min = floor((xv_min - xbp_min) / dcx); //  +1???
    iv_max = floor((xv_max - xbp_min) / dcx) +1;

    jv_min = floor((yv_min - ybp_min) / dcy); //  +1???
    jv_max = floor((yv_max - ybp_min) / dcy) +1;

    kv_min = floor((zv_min - zbp_min) / dcz); //  +1???
    kv_max = floor((zv_max - zbp_min) / dcz) +1;

    iv_min = (iv_min<0) ? 0:iv_min;
    iv_max = (iv_max>ncx) ? ncx:iv_max;

    jv_min = (jv_min<0) ? 0:jv_min;
    jv_max = (jv_max>ncx) ? ncy:jv_max;

    kv_min = (kv_min<0) ? 0:kv_min;
    kv_max = (kv_max>ncz) ? ncz:kv_max;

    // Insert IBM node information into a list
    for (k=kv_min; k<kv_max; k++) {
      for (j=jv_min; j<jv_max; j++) {
	for (i=iv_min; i<iv_max; i++) {
	  insertnode(&(cell_trg[k *ncx*ncy + j*ncx +i]), ln_v);
	}
      }
    }
  }

  PetscInt rank, flg=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  Cmpnts ***coor;
  DMDAVecGetArray(fda, user->gcent_search, &coor);
  DMDAVecGetArray(da, user->Nvert, &nvert);

  PetscBarrier(PETSC_NULL);
  // for this body nvert 4 is inside, 2 is near bndry
  // for previous bodies nvert sb*10 inside
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	
	if (coor[k][j][i].x > xbp_min && coor[k][j][i].x < xbp_max &&
	    coor[k][j][i].y > ybp_min && coor[k][j][i].y < ybp_max &&
	    coor[k][j][i].z > zbp_min && coor[k][j][i].z < zbp_max) {

	  ic = floor((coor[k][j][i].x - xbp_min )/ dcx);
	  jc = floor((coor[k][j][i].y - ybp_min )/ dcy);
	  kc = floor((coor[k][j][i].z - zbp_min )/ dcz);
	  
	  nvert[k][j][i] =  PetscMax(nvert[k][j][i],
				     blank_point_cell_advanced(coor[k][j][i], ic, jc, kc, ibm, ncx, ncy, ncz, dcx, dcy, xbp_min, ybp_min, zbp_max, cell_trg, flg));
	  if (nvert[k][j][i] < 0) nvert[k][j][i] = 0;
	}
      }
    }
  }
  
  DMDAVecRestoreArray(fda, user->gcent_search, &coor);


  DMDAVecRestoreArray(da, user->Nvert, &nvert);

  DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
  DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);


  DMDAVecGetArray(da, user->Nvert, &nvert);
  
  // Back to the old nvert 3 and 1 
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((int)(nvert[k][j][i]+0.5) == 4) nvert[k][j][i]=(blank_id[user->_this][ibi]-1)*100+10*(ibi+1);
      }
    }
  }

  DMDAVecRestoreArray(da, user->Nvert, &nvert);
  DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
  DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
  
  
  VecCopy(user->Nvert, user->Nvert_o);
  VecCopy(user->lNvert, user->lNvert_o);
  
  
  for (k=0; k<ncz; k++) {
    for (j=0; j<ncy; j++) {
      for (i=0; i<ncx; i++) {
	destroy(&cell_trg[k*ncx*ncy+j*ncx+i]);
      }
    }
  }

  PetscFree(cell_trg);
  
  PetscViewer	viewer;
  char filen[80];
  sprintf(filen, "nvfield_blank.dat");
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  VecView(user->Nvert, viewer);
  PetscViewerDestroy(&viewer);

  return 0;
}




PetscErrorCode blank_boundary_nvert(UserCtx *user, PetscInt sb){// blank boundary nverts forming 10*sb-1/10*sb-2
  
  DM	da = user->da, fda = user->fda;
  
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscReal	***nvert,***lnvert;
  PetscInt	i, j, k;
  
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;
  
  if (xs==0) lxs = xs;
  if (ys==0) lys = ys;
  if (zs==0) lzs = zs;
  
  if (xe==mx) lxe = xe;
  if (ye==my) lye = ye;
  if (ze==mz) lze = ze;
  
  DMDAVecGetArray(da, user->Nvert, &nvert);
  DMDAVecGetArray(da, user->lNvert, &lnvert);
  
  PetscInt ip, im, jp, jm, kp, km;
  PetscInt ii, jj, kk;
  // Near Boundary?
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	ip = (i<mx-2?(i+1):(i));
	im = (i>1   ?(i-1):(i));
	
	jp = (j<my-2?(j+1):(j));
	jm = (j>1   ?(j-1):(j));
	
	kp = (k<mz-2?(k+1):(k));
	km = (k>1   ?(k-1):(k));
	
	if ((int)(nvert[k][j][i]+0.5) == (blank_id[user->_this][sb]-1)*100+10*(sb+1)) {
	  for (kk=km; kk<kp+1; kk++) {
	    for (jj=jm; jj<jp+1; jj++) {
	      for (ii=im; ii<ip+1; ii++) {
		if ((int)(lnvert[kk][jj][ii] +0.5) == 0) {
		  nvert[k][j][i] = PetscMin((blank_id[user->_this][sb]-1)*100+10*(sb+1)-1, nvert[k][j][i]);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  DMDAVecRestoreArray(da, user->Nvert, &nvert);
  DMDAVecRestoreArray(da, user->lNvert, &lnvert);
  DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
  DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
  
  DMDAVecGetArray(da, user->Nvert, &nvert);
  DMDAVecGetArray(da, user->lNvert, &lnvert);
  
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	ip = (i<mx-2?(i+1):(i));
	im = (i>1   ?(i-1):(i));
	
	jp = (j<my-2?(j+1):(j));
	jm = (j>1   ?(j-1):(j));
	
	kp = (k<mz-2?(k+1):(k));
	km = (k>1   ?(k-1):(k));
	
	if ((int)(nvert[k][j][i]+0.5) == (blank_id[user->_this][sb]-1)*100+10*(sb+1)) {
	  for (kk=km; kk<kp+1; kk++) {
	    for (jj=jm; jj<jp+1; jj++) {
	      for (ii=im; ii<ip+1; ii++) {
		if ((int)(lnvert[kk][jj][ii] +0.5) == (blank_id[user->_this][sb]-1)*100+10*(sb+1)-1) {
		  nvert[k][j][i] = PetscMin((blank_id[user->_this][sb]-1)*100+10*(sb+1)-2, nvert[k][j][i]);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  
  DMDAVecRestoreArray(da, user->Nvert, &nvert);
  DMDAVecRestoreArray(da, user->lNvert, &lnvert);
  
  DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
  DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);


  VecCopy(user->Nvert, user->Nvert_o);
  VecCopy(user->lNvert, user->lNvert_o);
  return 0;
}



PetscErrorCode blank_ibm_read_ucd(IBMNodesblank *ibm, PetscInt ibi)
{
   PetscInt	rank;
  PetscInt	n_v , n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscInt	i,ii;
  PetscInt	n1e, n2e, n3e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal     dr;
  //Added 4/1/06 iman
  PetscReal     *dA ;//area
  PetscReal	*nt_x, *nt_y, *nt_z;
  PetscReal	*ns_x, *ns_y, *ns_z;
  PetscReal     cl=L_dim;
  char   ss[20];
  //double xt;
  char string[128];
  
 

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) { // root processor read in the data
    FILE *fd;
    PetscPrintf(PETSC_COMM_SELF, "READ ibmdata\n");
    char filen[80];  
    sprintf(filen,"blank_grid%2.2d.dat" , ibi);
 
    fd = fopen(filen, "r"); if (!fd) SETERRQ(PETSC_COMM_WORLD,1, "Cannot open IBM node file");
    n_v =0;

    if (fd) {
      fgets(string, 128, fd);
      fgets(string, 128, fd);
      fgets(string, 128, fd);
      
      fscanf(fd, "%i %i %i %i %i",&n_v,&n_elmt,&ii,&ii,&ii);
      PetscPrintf(PETSC_COMM_SELF, "number of nodes & elements %d %d\n",n_v, n_elmt);
      
      ibm->n_v = n_v;
      ibm->n_elmt = n_elmt;      
      
      MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      //    PetscPrintf(PETSC_COMM_WORLD, "nv, %d %e \n", n_v, xt);
      PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &z_bp);
      
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
      
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
      
      for (i=0; i<n_v; i++) {
	fscanf(fd, "%i %le %le %le", &ii, &x_bp[i], &y_bp[i], &z_bp[i]);//, &t, &t, &t);
	x_bp[i] = x_bp[i]/cl;
	y_bp[i] = y_bp[i]/cl;
	z_bp[i] = z_bp[i]/cl;
	
	if (ibi==0){
	  x_bp[i] = x_bp[i]*9.0/6.0;	
	  y_bp[i] = y_bp[i]*9/6.0-4.5;
          z_bp[i] = 1.3*z_bp[i]-3.0;
	
	}

	if (ibi==1){
	 // x_bp[i] = x_bp[i]*1.1;
	  y_bp[i] = y_bp[i]*0.8;
          z_bp[i] = .2*z_bp[i]-.20;
	
	}

	//+ ibi*2.;//5.;//2.;//8.;//6.;//2.   ;// 24.;
   // if (ibi==1)  x_bp[i] = x_bp[i]*1.047;
	
	ibm->x_bp[i] = x_bp[i];
	ibm->y_bp[i] = y_bp[i];
	ibm->z_bp[i] = z_bp[i];
	
	ibm->x_bp0[i] = x_bp[i];
	ibm->y_bp0[i] = y_bp[i];
	ibm->z_bp0[i] = z_bp[i];
      }
      
      i=0;
      PetscPrintf(PETSC_COMM_WORLD, "xyz_bp %le %le %le\n", x_bp[i], y_bp[i], z_bp[i]);



      MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);



      MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);
      
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
      
      // end added

      for (i=0; i<n_elmt; i++) {

	fscanf(fd, "%i %i %s %i %i %i\n", &ii,&ii, ss, nv1+i, nv2+i, nv3+i);
	nv1[i] = nv1[i] - 1; nv2[i] = nv2[i]-1; nv3[i] = nv3[i] - 1;

      }
      ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

      i=0;
      PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);

      fclose(fd);
    }
      
    for (i=0; i<n_elmt; i++) {
      
      n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i];
      dx12 = x_bp[n2e] - x_bp[n1e];
      dy12 = y_bp[n2e] - y_bp[n1e];
      dz12 = z_bp[n2e] - z_bp[n1e];
      
      dx13 = x_bp[n3e] - x_bp[n1e];
      dy13 = y_bp[n3e] - y_bp[n1e];
      dz13 = z_bp[n3e] - z_bp[n1e];
      
      nf_x[i] = dy12 * dz13 - dz12 * dy13;
      nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      nf_z[i] = dx12 * dy13 - dy12 * dx13;
      
      dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);
      
      nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;
      	
    }
    
    
    ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;
    
    //Added 4/1/06 iman
    
    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);

  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_v = n_v;
    
    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));


    MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);



      
    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_elmt = n_elmt;

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;

    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
  }
  PetscBarrier(NULL);

	for(i=0;i<ibm->n_v;i++)
	{
	   PetscPrintf(PETSC_COMM_SELF, "READ ibmdata  ibm%le    %le\n",ibm->x_bp0[i],ibm->x_bp[i]);
  
	}
  return(0);
}



//------------------------------------------------------------------
//ibm_read_Ansys is updated for icem Hafez 1/20/14
//
//
//






PetscErrorCode blank_ibm_read_Icem(IBMNodesblank *ibm, PetscInt ibi)
{
  PetscInt	rank;
  PetscInt	n_v , n_elmt ;
  PetscReal	*x_bp , *y_bp , *z_bp ;
  PetscInt	*nv1 , *nv2 , *nv3 ;
  PetscReal	*nf_x, *nf_y, *nf_z;
  PetscInt	i,ii;
  PetscInt	n1e, n2e, n3e;
  PetscReal	dx12, dy12, dz12, dx13, dy13, dz13;
  PetscReal     dr;
  //Added 4/1/06 iman
  PetscReal     *dA ;//area
  PetscReal	*nt_x, *nt_y, *nt_z;
  PetscReal	*ns_x, *ns_y, *ns_z;
  PetscInt      itr;
  PetscReal     tt;
  PetscReal     cl=L_dim;//168. for copepod;

  PetscOptionsGetReal(PETSC_NULL, "-char_length_ibm", &cl, PETSC_NULL);     
  
  char string[128];

 

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) { // root processor read in the data
    FILE *fd;
    PetscPrintf(PETSC_COMM_SELF, "READ nlist, %le\n", 1./cl);
    char filen[80];  
    sprintf(filen,"blank_nlist%2.2d" , ibi);


    fd = fopen(filen, "r"); if (!fd) SETERRQ(PETSC_COMM_WORLD,1, "Cannot open IBM node file");
 
  

    if (fd) {
      
      fscanf(fd, "%i",&n_v);
      PetscPrintf(PETSC_COMM_SELF, "number of nodes of list %d %d \n",ibi, n_v);
      
      ibm->n_v = n_v;
          
      
      MPI_Bcast(&(ibm->n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    
      PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
      PetscMalloc(n_v*sizeof(PetscReal), &z_bp);
      
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));
      
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));
      
      i=-1;
         fgets(string,128, fd);// skip line one

    	while (i+1<n_v) {

         i++;
	 fscanf(fd, "  %d %le %le %le\n",&ii, &(x_bp[i]), &(y_bp[i]), &(z_bp[i]));
	
	 x_bp[i] = x_bp[i]/cl ;//0.25 ;// 24.;
	 y_bp[i] = y_bp[i]/cl ;//2.;//8.;//6.;//2.   ;// 24.;
	 z_bp[i] = z_bp[i]/cl ;//2.;//8.;//15.;//2.   ;// 24.;
	 
	 ibm->x_bp[i]=x_bp[i];
	 ibm->y_bp[i]=y_bp[i];
	 ibm->z_bp[i]=z_bp[i];
	
	 ibm->x_bp0[i]=x_bp[i];
	 ibm->y_bp0[i]=y_bp[i];
	 ibm->z_bp0[i]=z_bp[i];
	
	}



      MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
 
      MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
      MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
      fclose(fd);
    }

    //Reading elements list
    PetscPrintf(PETSC_COMM_SELF, "READ elist\n");

    sprintf(filen,"blank_elist%2.2d" , ibi);
    fd = fopen(filen, "r"); if (!fd) SETERRQ(PETSC_COMM_SELF,1, "Cannot open IBM node file");

    if (fd) {
     
      
      fscanf(fd, "%i",&n_elmt);
      PetscPrintf(PETSC_COMM_SELF, "number of element of list %d %d \n",ibi, n_elmt);
      
      ibm->n_elmt = n_elmt;      
     
      
      MPI_Bcast(&(ibm->n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);
      
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);
      
      
      i=0;
      fgets(string, 128, fd);//skip one line 
      //    while(fgets(string, 128, fd)) {
	itr = 0;	

	while (i<n_elmt) {
	i++;

	fscanf(fd, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", &ii,&ii,&ii,&ii,&ii,&ii,&ii,&ii,&ii,&ii,&ii,&nv1[i-1], &nv2[i-1], &nv3[i-1],&ii);
	nv1[i-1] = nv1[i-1] - 1; nv2[i-1] = nv2[i-1]-1; nv3[i-1] = nv3[i-1] - 1;

	}
	//  } closing of first while loop
      ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;

      i=0;
      PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);
      i=30;
      PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);
      i=n_elmt-1;
      PetscPrintf(PETSC_COMM_WORLD, "nv %d %d %d\n", nv1[i], nv2[i], nv3[i]);

      fclose(fd);
    }

    x_bp=ibm->x_bp;  y_bp=ibm->y_bp ; z_bp = ibm->z_bp ;
    for (i=0; i<n_elmt; i++) {
      //PetscPrintf(PETSC_COMM_WORLD, "cop nf %d !\n",i);       
      n1e = nv1[i]; n2e =nv2[i]; n3e = nv3[i];
      dx12 = x_bp[n2e] - x_bp[n1e];
      dy12 = y_bp[n2e] - y_bp[n1e];
      dz12 = z_bp[n2e] - z_bp[n1e];
      
      dx13 = x_bp[n3e] - x_bp[n1e];
      dy13 = y_bp[n3e] - y_bp[n1e];
      dz13 = z_bp[n3e] - z_bp[n1e];
      
      nf_x[i] = dy12 * dz13 - dz12 * dy13;
      nf_y[i] = -dx12 * dz13 + dz12 * dx13;
      nf_z[i] = dx12 * dy13 - dy12 * dx13;
      
      dr = sqrt(nf_x[i]*nf_x[i] + nf_y[i]*nf_y[i] + nf_z[i]*nf_z[i]);
      
      nf_x[i] /=dr; nf_y[i]/=dr; nf_z[i]/=dr;
      
 
    }
          

    ibm->nf_x = nf_x; ibm->nf_y = nf_y;  ibm->nf_z = nf_z;

    
    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    
    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
  

  }
  else if (rank) {
    MPI_Bcast(&(n_v), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_v = n_v;
    
    PetscMalloc(n_v*sizeof(PetscReal), &x_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &y_bp);
    PetscMalloc(n_v*sizeof(PetscReal), &z_bp);
    
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
    PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));



    MPI_Bcast(ibm->x_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

 
    MPI_Bcast(ibm->x_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->y_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->z_bp0, n_v, MPIU_REAL, 0, PETSC_COMM_WORLD);

   
    MPI_Bcast(&(n_elmt), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    ibm->n_elmt = n_elmt;

    PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
    PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);

    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_x);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_y);
    PetscMalloc(n_elmt*sizeof(PetscReal), &nf_z);

    ibm->nv1 = nv1; ibm->nv2 = nv2; ibm->nv3 = nv3;
    ibm->nf_x = nf_x; ibm->nf_y = nf_y; ibm->nf_z = nf_z;


    MPI_Bcast(ibm->nv1, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv2, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nv3, n_elmt, MPI_INT, 0, PETSC_COMM_WORLD);

    MPI_Bcast(ibm->nf_x, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_y, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    MPI_Bcast(ibm->nf_z, n_elmt, MPIU_REAL, 0, PETSC_COMM_WORLD);
    
  }
	for(i=0;i<ibm->n_v;i++)
	{
	   PetscPrintf(PETSC_COMM_SELF, "READ ibmdata  ibm%le    %le\n",ibm->x_bp0[i],ibm->x_bp[i]);
  
	}

  return(0);
}
//---------------------------------------------------------------------- blank pressure correction


PetscErrorCode Blank_Pressure_Correction(UserCtx *user, PetscInt flg) 
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
  
  PetscReal ***nvert,***p,***lp, ibmval;
  Cmpnts	***ucont,***ucont_o;
  
  DMDAVecGetArray(da, user->P, &p);
  DMDAVecGetArray(da, user->lP, &lp);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  
  DMDAVecGetArray(fda, user->lUcont, &ucont);
  DMDAVecGetArray(fda, user->lUcont_o, &ucont_o);

  for (sb=0; sb<NumberOfBlank; sb++) {
    ibmval = (blank_id[user->_this][sb]-1)*100+10*(sb+1)-1;
    
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {

	if (nvert[k][j][i] > ibmval-0.5 && nvert[k][j][i] < ibmval+0.5) {
	  if (nvert[k][j][i-1] < 0.1 && i> 1) {
	    p[k][j][i]=lp[k][j][i-1];
	    
	  }
	  ///
	  if (nvert[k][j-1][i] < 0.1 && j>1) {
	    p[k][j][i]=lp[k][j-1][i];
	  }
	  ///
	  if (nvert[k-1][j][i] < 0.1 && k>1) {
	    p[k][j][i]=lp[k-1][j][i];
	  }
	}
	
	/////////
	
	if (nvert[k][j][i] > ibmval-0.5 && nvert[k][j][i] < ibmval+0.5) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    p[k][j][i]=lp[k][j][i+1];
	  }
	  ///
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    p[k][j][i]=lp[k][j+1][i];
	  }
	  ///
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    p[k][j][i]=lp[k+1][j][i];
	  }
	}

      }
    }
  }

 
  } //sb

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da, user->P, &p);
  DMDAVecRestoreArray(da, user->lP, &lp);
  
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);
  DMDAVecRestoreArray(fda, user->lUcont_o, &ucont_o);
  
  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

  
 return 0;
}


PetscErrorCode Copy_Coor(UserCtx *user) {
  Vec           Coor;
  PetscInt      bi;
  Cmpnts	***coords,***coor;
  PetscInt	xs, ys, zs, xe, ye, ze;
  PetscInt	i, j, k;

  for (bi=0; bi<block_number; bi++) {
    DM		da = user[bi].da, fda = user[bi].fda;
    
    DMDALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
    
    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coords);
    DMDAVecGetArray(fda, user[bi].g_Coor, &coor);
    
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  coor[k][j][i].x= coords[k][j][i].x;
	  coor[k][j][i].y= coords[k][j][i].y;
	  coor[k][j][i].z= coords[k][j][i].z;
	}
      }
    }
    

    DMDAVecRestoreArray(fda, Coor, &coords);
    DMDAVecRestoreArray(fda, user[bi].g_Coor, &coor);

    DMGlobalToLocalBegin(fda, user[bi].g_Coor, INSERT_VALUES, user[bi].Coor);
    DMGlobalToLocalEnd(fda, user[bi].g_Coor, INSERT_VALUES, user[bi].Coor);
    //
    VecCopy( user[bi].Coor,user[bi].Coor_init);
  }
  
  

  return 0;
}


PetscErrorCode Blank_flux_calc(UserCtx *user, PetscInt flg) //calculate flux on nvert=1010
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
  PetscReal ***nvert, ibmval=0.0,ibm_Flux=0.0,ibm_Area=0.0;
  Cmpnts ***ucor, ***csi, ***eta, ***zet;
  DMDAVecGetArray(fda, user->Ucont, &ucor);
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  PetscReal libm_Flux=0.0, libm_area=0.0, libm_Flux_abs=0.0, ibm_Flux_abs=0.0;

  for (sb=0; sb<NumberOfBlank; sb++) {
    
    if(blank_id[user->_this][sb]==11){
      
      ibmval=(blank_id[user->_this][sb]-1)*100+10-1.;
      
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
      
      PetscPrintf(PETSC_COMM_WORLD, "BLANK-FLUX1000: %le %le\n", ibm_Flux, ibm_Area);
    }        
  }//sb
  
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(fda, user->Ucont, &ucor);
  
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  ////////////////
  user->blank_FluxIntpSum=ibm_Flux;
  ////////////////
 
 return 0;
}
