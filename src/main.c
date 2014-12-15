#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>
#include<time.h>
#include<string.h>

#include "splinesurf.h"
#include "jacobi.h"
#include "Jflux.h"
#include "readf.h"
#include "morse_vib.h"
#include "fourier.h"

/*
Program: Raman

Authors: Freddy Fernandes Guimarães
         Vinicius Vaz da Cruz - viniciusvcruz@gmail.com

History: Based on FluxVE written by Vinicius

Purpose: Compute raman cross section from a eSPec wavepacket propagation

Method: Perform FFT from |f(R,r,t)> to |f(R,r,E)>, then choose a selected energy wavepacket to propagate in the final state

Goiania, 08th of december of 2014
 */

#define NAV   6.0221415E+23  //number of entities in 1 mol
#define Me    9.10938188E-31 //mass of electron in Kg
#define FSAU(x) (x)*41.3413745758e+0 ///0.024189E+0 //fs to au units of time conversion                        mtrxdiag_ (char *dim, int *il, int *iu, int *info, int *mxdct, int *n, double *abstol, int *iwork, int *np, double *eigvl, double *shm, double *vpot, double *work, double *wk, double *eigvc);

int rdinput(FILE *arq,char *dim,int *np,char *file,int *stfil,int *endfil, double *ti, double *tf, double *pti, double *ptf, double *pstept, double *m,char *potfile, int *nf,int *twopow, double *width, double *Ef, int *type, double *crosst);

//void readallspl_(char *file,int *iflth,double *X,double *Y,double *T,double *bcoefRe,double *bcoefIm,double *xknot,double *yknot,double *tknot,int *nx, int *np[0], int *nf, int *kx, int *nxr, int *np[0], int *fpr, double *norm, int *prtRMSwhole, int *noreg, int *wholegrid, int *stfil);

int main(){
  int i,j,k,l,np[3],nf,nxr,nfr,xs,ys,ks,ldf,mdf,nw,nz,nwr,nzr;
  int fp[10],fpr[10],stfil,endfil,iflth,nrmse,type,noreg;
  int wholegrid,twopow,prtwpE,prtRMSwhole;
  double norm;
  double ti,tf,pti,ptf,pstept,xi,xf,yi,yf,stept,sh[3],width,potshift,maxF[2],maxaE[2],stepw,stepz;
  double *X,*Y,*T,rmse,xwork,*W,*Z,mem;
  double m[3],x,y,pk,pki,pkf,steppk,ansk,val[5];
  char axis,file[30],potfile[30],wp_Enam[20],num[20],dim[6];
  FILE *arq=fopen("raman.inp","r");
  FILE *deb=fopen("debug.dat","w");
  //--- spline variables
  int NTG,fpspl[4],ftinfo,kx;
  double *bcoefre,*bcoefim,*xknot,*yknot,*tknot;
  double stepxspl,stepyspl,steptspl,Xspl[3],Yspl[3],Respl[3],Imspl[3];
  //-- initial wavepacket parameters
  double crosst;
  //--- fourier variables
  int nE;
  double *WPERe,*WPEIm,*E,Ei,aE,stepE,*eknot,Ef;
  fftw_complex *workin,*workout,*workk; //*aE;
  double kmin,kmax,Emin,Emax;
  fftw_plan p;
  clock_t begint,endt;
  double timediff;
  FILE *wp_E;


  //--- Default values
  type=0;
  noreg=1;
  np[1] = 1;
  wholegrid=1;
  stepxspl=1.0E-3;
  stepyspl=1.0E-3;
  twopow=9; //11;
  width = 1.0e-9;
  kx=6.0e+0;
  prtwpE=1; // =0 to print transformed wp into files -> NEED TO ADD TO rdinput
  prtRMSwhole=1; // = 0 to compute RMSE in relation to whole data -> NEED TO ADD TO rdinput

  //-------------------------------------------------//
  printf("\n =========================================\n");
  printf("||                                         ||\n");
  printf("||                  Raman                  ||\n");
  printf("||                                         ||\n");
  printf(" ==========================================\n");
  printf("Vinicius V. Cruz and Freddy F. Guimaraes\n\n");
  printf("Reading input parameters...\n\n");

  //--- Read Input File
  rdinput(arq,dim,np,file,&stfil,&endfil,&ti,&tf,&pti,&ptf,&pstept,m,potfile,&nf,&twopow,&width,&Ef,&type,&crosst);
  fclose(arq);
  //-----

  //--- Units Conversion
  //m1=(m1*1.0E-3)/(NAV * Me);
  //m2=(m2*1.0E-3)/(NAV * Me);
  //conversion below used in eSPec
 
  sh[0]=(yf -yi)/(np[0] -1.0);
  sh[1]=(xf -xi)/(np[1] -1.0);
  stept=(tf - ti)/(nf -1.0);
  for(i=0;i<3;i++)m[i] = 1.836152667539e+3*m[i];

  //NTG = pow(2,twopow+1);
  NTG = pow(2,twopow);
  //steptspl=(2.0*tf)/(NTG); //-1.0e+0);
  steptspl=(tf - ti)/(NTG);
  //--- spline parameters determination 


  //--- Input Parameters -----------------------------------------------------------------------//
  printf("\n<< Starting input section >>\n\n");

  /*  printf("> Run type:\n");
  if(type==0){
    printf("This is a non-reactive calculation\n");
  }else if(type==1){
    printf("This is a reactive flux calculation\n");
    printf("a Jacobi coordinate transformation will be performed \n\n");
     wholegrid = 0;
     }*/

  printf(">Grid parameters:\n");
  printf("dimension: %s \n", dim);
  printf("npoints: %d ",np[0]);
  if(strncasecmp(dim,".2D",3)==0)printf("%d \n", np[1]);
  printf("\n");
  //printf("npoints x %d \n", np[1]);
  //printf("npoints y %d \n", np[0]);
  //printf("xrange: %lf %lf\n", xi, xf);
  //printf("yrange %lf %lf\n", yi, yf);
  //printf("step X %.7lf Bohr\n",stepx);
  //printf("step Y %.7lf Bohr\n", stepy);
  
  printf("\n>Wavepacket files:\n");
  printf("File name: %s\n",file);
  printf("ti %lf tf %lf step %lf fs nfile: %d \n",ti,tf,stept,nf);
  printf("masses m1 = %lf a.u.\n", m[0]);
  if(strncasecmp(dim,".2D",3)==0) printf("       m2 = %lf a.u.\n", m[1]);
  if(strncasecmp(dim,".2DCT",5)==0) printf("       m3 = %lf a.u.\n", m[2]);  

  printf("\n>Final state propagation:\n");
  printf("final state potential file: %s \n",potfile);
  printf("propagation time: %lf %lf, step: %lf \n",pti,ptf,pstept);
  printf("Energy: %lf \n",Ef);

  printf("\n>Fourier transformparameters:\n");
  printf("grid is 2^%d, time step = %E \n",twopow,steptspl);
  stepE = 2*M_PI/(NTG*FSAU(steptspl));
  Ei = -NTG*stepE/(2.0E+0);
  nE = NTG;
  printf("Energy step = %E a.u., Ef = %E a.u. \n\n",stepE,-Ei);


  mem = (4.0*(np[0]*np[1]*nf*sizeof(double))/(1024*1024*1024) + 2*(np[0]+np[1]+nf+NTG)*sizeof(double) + 2.0*(np[1]*np[0]*NTG*sizeof(double)))/(1024*1024*1024);
  printf("Memory requirement estimation: %E GB\n",mem);
  
  
  printf("\nFinished input section!\n");
  printf("------------------------------------------------------\n\n\n");
  

  /*

  //--- Allocate arrays
  bcoefre = malloc(np[0]*np[1]*nf*sizeof(double));
  bcoefim = malloc(np[0]*np[1]*nf*sizeof(double));
  Y = malloc(np[0]*sizeof(double));
  X = malloc(np[1]*sizeof(double));
  T = malloc(nf*sizeof(double));
  xknot = malloc((np[1]+kx)*sizeof(double));
  yknot = malloc((np[0]+kx)*sizeof(double));
  tknot = malloc((nf+kx)*sizeof(double));
  workk = fftw_malloc(sizeof(fftw_complex) * NTG);
  WPERe = malloc(np[0]*np[1]*(NTG)*sizeof(double)); 
  WPEIm = malloc(np[0]*np[1]*(NTG)*sizeof(double)); 

  printf("\n<< Starting reading and spline section >>\n");
  
  // eSPec output files have unitary euclidean norm, in order to renormalize it to the integral norm:
  norm = 1.0e+0/sqrt(stepx*stepy);
  //norm = 1.0;

  // convert time variables from fs to au
  ti = FSAU(ti);
  tf = FSAU(tf);
  stept = FSAU(stept);
  steptspl = FSAU(steptspl);

  //--- Read all wavepackets and generate spline coefficient files
  iflth = strlen(file);
  readallspl_(file,&iflth,X,Y,T,bcoefre,bcoefim,xknot,yknot,tknot,&np[1],&np[0],&nf,&kx,&np[1],&np[0],fpr,&norm,&prtRMSwhole,&noreg,&wholegrid,&stfil);

  printf("\n<< Starting Fourier section >>\n");

  begint = clock();
  ftinfo=run_all_fft(bcoefre,bcoefim,X,Y,xknot,yknot,tknot,kx,ti,tf,steptspl,np[1],np[0],nf,NTG,fpr,width,WPERe,WPEIm,maxF[0]);
  endt = clock();
  timediff = (double)(endt - begint)/CLOCKS_PER_SEC;

  stepE = 2*M_PI/(NTG*steptspl);
  Ei = -NTG*stepE/(2.0E+0);
  nE = NTG;
  
  if(ftinfo!=0){
    printf("\n An error occured while computing the wavepacket's Fourier transform: The program will stop!\n");
    return 666;
  }
  printf("\n Wavepacket has been trasformed to the Energy domain!\n");
  printf("\n it took %lf seconds \n",timediff);


   //checkmatrix_ (WPERe,WPEIm,&np[1],&np[0],&nE,X,Y,E);

  printf("\nGenerating 3D spline representation of the fourier transformed wavepacket\n");
  free(bcoefre);
  free(bcoefim);
  //free(T); can I free T ???
  free(xknot);
  free(yknot);
  free(tknot);
  xknot = malloc((np[1]+kx)*sizeof(double));
  yknot = malloc((np[0]+kx)*sizeof(double));
  eknot = malloc((nE +kx)*sizeof(double));
  E = malloc(nE*sizeof(double));
  bcoefre = malloc(np[0]*np[1]*nE*sizeof(double));
  bcoefim = malloc(np[0]*np[1]*nE*sizeof(double));

  for(k=0;k<nE;k++){
    E[k] = Ei + k*stepE;
  }

  //---
  begint = clock();
  dbsnak_ (&np[0], Y, &kx, yknot);
  dbsnak_ (&np[1], X, &kx, xknot);
  dbsnak_ (&nE, E, &kx, eknot);
  //----------
  dbs3in_ (&np[0],Y,&np[1],X,&nE,E,WPERe,&np[0],&np[1],&kx,&kx,&kx,yknot,xknot,eknot,bcoefre);
  dbs3in_ (&np[0],Y,&np[1],X,&nE,E,WPEIm,&np[0],&np[1],&kx,&kx,&kx,yknot,xknot,eknot,bcoefim);
  endt = clock();
  timediff = (double)(endt - begint)/CLOCKS_PER_SEC;
  //----------

  //free(WPERe);
  //free(WPEIm);

  printf("3D spline matrices calculated!\n");
  printf("it took %lf seconds\n",timediff);

  prtwpE=1;
  //-----print energy wavepackets ----------------------
  if(prtwpE==0){
    l=0;
    for(k=nE/2;k<nE;k=k + 3){
      strcpy(wp_Enam,"wpE_");
      sprintf(num,"%d.dat",l+1);
      strcat(wp_Enam,num);
      wp_E = fopen(wp_Enam,"w");
      val[2] = Ei + k*stepE;
      fprintf(wp_E,"# Energy = %E \n",val[2]);

      for(i=0;i<np[1];i++){
	for(j=0;j<np[0];j++){
	  val[0] = dbs3vl_ (&Y[j],&X[i],&val[2],&kx,&kx,&kx,yknot,xknot,eknot,&np[0],&np[1],&nE,bcoefre);
	  val[1] = dbs3vl_ (&Y[j],&X[i],&val[2],&kx,&kx,&kx,yknot,xknot,eknot,&np[0],&np[1],&nE,bcoefim);
	  fprintf(wp_E,"%E %E %E %E\n",X[i],Y[j],val[0],val[1]);
	}
	fprintf(wp_E,"\n");
      }
      fclose(wp_E);
      l=l+1;
    }
  }


  printf("\n Finished Fourier section\n");
  printf("------------------------------------------------------\n\n\n");
  */
  
  printf("\n# End of Calculation!\n");
  printf("# Raman terminated successfully!\n");
  return 0;
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
