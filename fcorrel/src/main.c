#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>
#include<time.h>
#include<string.h>
#include <unistd.h>

#include "splinesurf.h"
#include "fourier.h"

/*
Program: 2D+1D RIXS

Authors: Vinicius Vaz da Cruz - viniciusvcruz@gmail.com 
         Freddy Fernandes Guimarães
	 Rafael Couto
	 Victor Kimberg
	 Faris Gel'mukhanov
	
History: Based on Raman written by Vinicius

Purpose: Compute the final RIXS cross section for a 2D+1D model

Goiania, 01st of February of 2015
 */

#define NAV   6.0221415E+23  //number of entities in 1 mol
#define Me    9.10938188E-31 //mass of electron in Kg
#define FSAU(x) (x)*41.3413745758e+0 ///0.024189E+0 //fs to au units of time conversion                        mtrxdiag_ (char *dim, int *il, int *iu, int *info, int *mxdct, int *n, double *abstol, int *iwork, int *np, double *eigvl, double *shm, double *vpot, double *work, double *wk, double *eigvc);

int rdinput(FILE *arq,char *dim,int *np,char *file,int *stfil,int *endfil, double *ti, double *tf, double *pti, double *ptf, double *pstept, double *m,char *potfile, int *nf,int *twopow, double *width, int *nEf,double *Ef, int *type, double *crosst,char *windtype,char *jobnam, int *nfunc, char *funam,int *nvc,int *nvf,double *Evf,double *Evc,char *fcnam, char *fcornam, double *Delta,double *shift, char *rexfnam, char *xasfnam, int *nxas, int *nrexs, double *omega);

void readallspl(char *file,int *iflth,double *X,double *Y,double *T,double *bcoefRe,double *bcoefIm,double *xknot,double *yknot,double *tknot, int *np, int *nf, int *kx, int *stfil, char *dim);
  


int main(){
  int i,j,k,l,ii,np[3],nf,nxr,nfr,xs,ys,ks,ldf,mdf,nw,nz,nwr,nzr,rdbcf,nfunc,nvc,nvf;
  int fp[10],fpr[10],stfil,endfil,iflth,nrmse,type,nEf,nshift;
  int wholegrid,twopow,prtwpE,prtRMSwhole,qop;
  double norm,t,s;
  double ti,tf,pti,ptf,pstept,xi,xf,yi,yf,stept,sh[3],width,potshift,stepw,stepz,Delta;
  double *X,*Y,*T,rmse,xwork,*W,*Z,mem,Ereso,*FCGVc,**FCVcVf,*intens;
  double m[3],x,y,pk,pki,pkf,steppk,ansk,val[5],FC,Evf[20],Evc[20];
  char axis,file[30],potfile[30],wp_Enam[20],num[20],dim[6],windtype[10],fnam[20],fnam2[20],funam[20],jobnam[50];
  char fcnam[20],fcornam[20];
  FILE *arq=fopen("correl.inp","r");
  FILE *deb=fopen("debug.dat","w");
  //--- spline variables
  int NTG,fpspl[4],ftinfo,kx;
  double *bcoefre,*bcoefim,*xknot,*yknot,*tknot;
  double stepxspl,stepyspl,steptspl,Xspl[3],Yspl[3],Respl[3],Imspl[3];
  //-- initial wavepacket parameters
  double crosst;
  //--- fourier variables
  int nE;
  double *E,Ei,aE,stepE,*eknot,Ef[200],shift,***inWPRe,***inWPIm,**fcorrelRe,**fcorrelIm,*tfcorrelRe,*tfcorrelIm;
  fftw_complex *workin,*workout,*workk; //*aE;
  double kmin,kmax,Emin,Emax;
  fftw_plan p;
  clock_t begint,endt;
  double timediff,window,taux;
  FILE *wpin,*fcf,*fcorr,*deb2;
  //---- self absorption
  FILE *rexsfile,*xasfile;
  char rexfnam[20],xasfnam[20];
  int nxas,nrexs;
  double *xas_cross,*xas_bcoef,*xas_omega, *xas_knot,xas_cross_om,xas_cross_omp;
  double *rexs_cross_sa,*rexs_cross,*rexs_omegap,omega;


  //--- Default values
  sprintf(jobnam,"default_jobname");
  type=0;
  np[1] = 1;
  nshift=0;
  rdbcf=1;
  wholegrid=1;
  stepxspl=1.0E-3;
  stepyspl=1.0E-3;
  twopow=9; //11;
  width = 1.0e-9;
  kx=6.0e+0;
  strcpy(windtype,"no input");
  qop = 1;
  prtwpE=1;       // =0 to print transformed wp into files -> NEED TO ADD TO rdinput
  prtRMSwhole=1; // = 0 to compute RMSE in relation to whole data -> NEED TO ADD TO rdinput

  //-------------------------------------------------//
  printf("\n =========================================\n");
  printf("||                                         ||\n");
  printf("||              2D + 1D RIXS               ||\n");
  printf("||                                         ||\n");
  printf(" ==========================================\n");
  printf("Vinicius V. Cruz and Freddy F. Guimaraes\n\n");
  printf("Reading input parameters...\n\n");

  //--- Read Input File
  rdinput(arq,dim,np,file,&stfil,&endfil,&ti,&tf,&pti,&ptf,&pstept,m,potfile,&nf,&twopow,&width,&nEf,Ef,&type,&crosst,windtype,jobnam,&nfunc,funam,&nvc,&nvf,Evf,Evc,fcnam,fcornam,&Delta,&shift,rexfnam,xasfnam,&nxas,&nrexs,&omega);
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
  printf("> Job name: %s \n",jobnam);
  if(type != 3) printf("number of points in the correlation functions: %d \n",nf);

  if(type == 0){
    printf(">Grid parameters:\n");
    printf("dimension: %s \n", dim);
    printf("npoints: %d ",np[0]);
    if(strncasecmp(dim,".2D",3)==0)printf("%d \n", np[1]);
    printf("\n");

    printf("\n>Wavepacket files:\n");
    printf("File name: %s\n",file);
    printf("ti %lf tf %lf step %lf fs nfile: %d \n",ti,tf,stept,nf);
    printf("masses m1 = %lf a.u.\n", m[0]);
    if(strncasecmp(dim,".2D",3)==0) printf("       m2 = %lf a.u.\n", m[1]);
    if(strncasecmp(dim,".2DCT",5)==0) printf("       m3 = %lf a.u.\n", m[2]);  
 
    printf("\n>Initial wavefunctions:\n");
    printf("the correlation function will be calculated in relation to %d initial wavefunctions",nfunc);
    printf("the wavefunctions will ve read from file %s",funam);

  }

  printf("\nFinished input section!\n");
  printf("------------------------------------------------------\n\n\n");
 

 

  if(type==0){
    //--------------------- select run type
    
    // convert time variables from fs to au
    ti = FSAU(ti);
    tf = FSAU(tf);
    stept = FSAU(stept);
    steptspl = FSAU(steptspl);

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
    inWPRe =  malloc(nfunc*sizeof(double)); 
    inWPIm =  malloc(nfunc*sizeof(double)); 
    for(i=0;i<nfunc;i++){
      inWPRe[i] =  malloc(np[1]*sizeof(double)); 
      inWPIm[i] =  malloc(np[1]*sizeof(double)); 
      for(j=0;j<np[1];j++){
	inWPRe[i][j] =  malloc(np[0]*sizeof(double)); 
	inWPIm[i][j] =  malloc(np[0]*sizeof(double)); 
      }
    }

    fcorrelRe=malloc(nfunc*sizeof(double));
    fcorrelIm=malloc(nfunc*sizeof(double));
    for(i=0;i<nfunc;i++){
      fcorrelRe[i] = malloc(nf*sizeof(double));
      fcorrelIm[i] = malloc(nf*sizeof(double));
    }



    printf("\n<< Starting reading and spline section >>\n");

    //------- read initial wavefucntions
    printf("reading initial wavefunctions \n");
    wpin = fopen(funam,"r");
    for(k=0;k<nfunc;k++){
      for(i=0;i<np[1];i++){
	for(j=0;j<np[0];j++){
	  if(strncasecmp(dim,".2D",3)==0) fscanf(wpin,"%lf %lf %lf %lf",&X[i],&Y[j],&inWPRe[k][i][j],&inWPIm[k][i][j]);
	  else fscanf(wpin,"%lf %lf %lf",&Y[j],&inWPRe[k][i][j],&inWPIm[k][i][j]);
	}
      }
    }
    fclose(wpin);
    printf("initial wavefunctions read! \n\n");

    //--- Read all wavepackets and generate spline coefficient matrices
    iflth = strlen(file);
    readallspl_(file,&iflth,X,Y,T,bcoefre,bcoefim,xknot,yknot,tknot,np,&nf,&kx,&stfil,dim);

    printf("\n << Starting Correlation section >>  \n");
    val[0] = 0e+0;
    val[1] = 0e+0;
    for(k=0;k<nfunc;k++){
      for(l=0;l<nf;l++){
	//----------------------------------
	for(i=0;i<np[1];i++){
	  for(j=0;j<np[0];j++){
	  
	    if(strncasecmp(dim,".2D",3)==0){
	      val[2] = dbs3vl_ (&Y[j],&X[i],&T[l],&kx,&kx,&kx,yknot,xknot,tknot,&np[0],&np[1],&nf,bcoefre);
	      val[3] = dbs3vl_ (&Y[j],&X[i],&T[l],&kx,&kx,&kx,yknot,xknot,tknot,&np[0],&np[1],&nf,bcoefim);
	    }else if(strncasecmp(dim,".1D",3)==0){
	      val[2] = dbs2vl_ (&Y[j],&T[l],&kx,&kx,yknot,tknot,&np[0],&nf,bcoefre);
	      val[3] = dbs2vl_ (&Y[j],&T[l],&kx,&kx,yknot,tknot,&np[0],&nf,bcoefim);
	    }

	    val[0] = val[0] + inWPRe[k][i][j]*val[2] + inWPIm[k][i][j]*val[3];
	    val[1] = val[1] + inWPRe[k][i][j]*val[3] - inWPIm[k][i][j]*val[2];
	  }
	}
	fcorrelRe[k][l] = val[0];
	fcorrelIm[k][l] = val[1];
	val[0] = 0.0e+0;
	val[1] = 0.0e+0;
	//----------------------------------
      }
    }
  
    printf("All desired correlation functions have been computed! \n");
    printf("------------------\n");
    for(l=0;l<nf;l++){
      printf("%lf ",T[l]/41.3413745758e+0);
      for(k=0;k<nfunc;k++) printf("%lf %lf ",fcorrelRe[k][l],fcorrelIm[k][l]);
      printf("\n");
    }

//----------------END------------------------------
//------------------------------------------------
  } else if(type==1){
    printf("\n << Running 2D+1D model: final cross section >>\n");

    printf("\n << Starting Franck-Condon section >> \n\n");
    printf("number of vc states: %d ,   number of vf states: %d \n",nvc,nvf);
    FCGVc = malloc(nvc*sizeof(double));
    FCVcVf = malloc(nvc*sizeof(double));
    intens = malloc(nvc*sizeof(double));

    for(i=0;i<nvc;i++) FCVcVf[i] = malloc(nvf*sizeof(double));

    printf("reading <0|vc> and <vc|vf> FC matrices from file %s \n",fcnam);

    fcf = fopen(fcnam,"r");
    printf("\n<0|vc>\n   ");
    for(i=0;i<nvc;i++)  printf("      |%d> ",i);
    printf("\n");
    printf("< %d| ",0);
    for(i=0;i<nvc;i++){
      fscanf(fcf,"%lf",&FCGVc[i]); 
      printf("% lf ", FCGVc[i]);
    }

    printf("\n\n<vc|vf>\n   ");
    for(i=0;i<nvf;i++) printf("      |%d> ",i);
    printf("\n");

    for(i=0;i<nvc;i++){
      printf("< %d| ",i);
      for(j=0;j<nvf;j++){
	fscanf(fcf,"%lf",&FCVcVf[i][j]);
	printf("% lf ",FCVcVf[i][j]);
      }
      printf("\n");
    }

    // reading the normalization factors from the wf generated with the fft
    for(i=0;i<nvc;i++){
      fscanf(fcf,"%lf",&intens[i]);
    }

    fclose(fcf);


    printf("\n");
    printf("final state bending energy values:\n");
    for(i=0;i<nvf;i++) printf(" Evf%i = %E , ",i,Evf[i]);
    printf("\n");

    printf("\n> Franck-Condon Section finished \n");


    printf("\n << Starting Auto-Correlation Function section >> \n");

    printf("\n the correlation functions will be read from file %s \n",fcornam);

    T = malloc(nf*sizeof(double));
    fcorrelRe=malloc(nvc*nvc*sizeof(double));
    fcorrelIm=malloc(nvc*nvc*sizeof(double));

    for(i=0;i<nvc*nvc;i++){
      fcorrelRe[i] = malloc(nf*sizeof(double));
      fcorrelIm[i] = malloc(nf*sizeof(double));
    }

    // read all correlation functions
    fcorr = fopen(fcornam,"r");
    for(j=0;j<nvc;j++){
      for(ii=0;ii<nf;ii++){
	fscanf(fcorr,"%lf",&T[ii]);
	T[ii] = T[ii] * 41.3413745758e+0;
	for(l=0;l<nvc;l++){
	  k = l + j*nvc;
	  fscanf(fcorr,"%lf %lf",&fcorrelRe[k][ii],&fcorrelIm[k][ii]);
	}
      }
    }
    fclose(fcorr);

    stept = (T[nf-1] - T[0])/(nf - 1.0);

    // total correlation function
    tfcorrelRe = malloc(nf*sizeof(double));
    tfcorrelIm = malloc(nf*sizeof(double));
    for(ii=0;ii<nf;ii++){
      tfcorrelRe[ii] = 0.0e+0;
      tfcorrelIm[ii] = 0.0e+0;
    }

    
    //debug viktor problem
    deb2=fopen("debug2.dat","w");

    // eq. (51) from file
    for(i=0;i<nvf;i++){
      for(j=0;j<nvc;j++){
	for(l=0;l<nvc;l++){
	  k = l + j*nvc;
	  // franck-condon factors
	  FC = FCGVc[l]*FCVcVf[l][i]*FCVcVf[j][i]*FCGVc[j];
	  // we multiply by the intensity, due to previous normalization of wf in the energy domain
	  fprintf(deb2,"vf = %d, vc = %d vc' = %d (%d), FC = %E \n",i,j,l,k,FC);
	  FC = FC * intens[j]*intens[l];
	  for(ii=0;ii<nf;ii++){
	    tfcorrelRe[ii] = tfcorrelRe[ii] + FC*(fcorrelRe[k][ii]*cos(Evf[i]*T[ii]) + fcorrelIm[k][ii]*sin(Evf[i]*T[ii]) );
	    tfcorrelIm[ii] = tfcorrelIm[ii] + FC*(fcorrelIm[k][ii]*cos(Evf[i]*T[ii]) - fcorrelRe[k][ii]*sin(Evf[i]*T[ii]) );
	    //debug below
	    //tfcorrelRe[ii] = tfcorrelRe[ii] + FC*fcorrelRe[k][ii];
	    //tfcorrelIm[ii] = tfcorrelIm[ii] + FC*fcorrelIm[k][ii];
	  }
	}
      }
    }

    printf("Final Total Auto-Correlation Function \n");
    printf("------------------\n");
    for(ii=0;ii<nf;ii++){
      printf("%lf % lf % lf \n",T[ii],tfcorrelRe[ii],tfcorrelIm[ii]);
    }

    printf("\n> Auto-Correlation Section finished \n");

    printf("\n << Starting Final Cross-Section section >> \n");

    //NTG = pow(2,twopow+1);
    NTG = pow(2,twopow);
    steptspl=(2.0*T[nf-1])/NTG;
    //stepE = 2*M_PI/(NTG*FSAU(steptspl));
    stepE = 2*M_PI/(NTG*steptspl);
    Ei = -NTG*stepE/(2.0E+0);
    nE =  NTG;

    printf("\nAuto-Correlation Function will be adjusted for Fourier Transform using splines \n");
    printf("\nFourier parameters \n");
    printf("npoints: 2^%d = %d\n",twopow,NTG);
    printf("time step: %lf fs \n",steptspl);
    printf("Energy step: %lf a.u. , final energy: %lf a.u. \n",stepE,-Ei);
    //steptspl=(tf - ti)/(NTG);

    //spline arrays
    bcoefre = malloc(nf*sizeof(double));
    bcoefim = malloc(nf*sizeof(double));
    tknot = malloc((nf+kx)*sizeof(double));
    //fft arrays
    workin = fftw_malloc(sizeof(fftw_complex) * NTG);
    workout = fftw_malloc(sizeof(fftw_complex) * NTG);

    dbsnak_ (&nf, T, &kx, tknot);
    dbsint_ (&nf,T,tfcorrelRe,&kx,tknot,bcoefre);
    dbsint_ (&nf,T,tfcorrelIm,&kx,tknot,bcoefim);

    //generate mirrored data
    for(i=0;i<NTG;i++){
      t = -T[nf-1] + i*steptspl;
      fprintf(deb,"%lf ",t);
      if(t < 0) s = -1.0e+0;
      else s = 1.0e+0;

      t = fabs(t);

      if(strncasecmp(windtype,".SGAUSS",7)==0){
	//TAUX = TMAX**2/(LOG(ONE/WP))**2 espec
	//taux = pow(T[nf-1],2)/(log(1.000/width)/log(M_E));
	taux = pow(T[nf-1],2)/pow(log(1.000/width),2);
	window = exp(-pow(t,2)/taux);
      }else if(strncasecmp(windtype,".EXPDEC",7)==0){
	window = exp(-width*t);
      }

      //debug
      //window=1.000;
      workin[i][0] =   window*dbsval_ (&t,&kx,tknot,&nf,bcoefre);
      workin[i][1] = s*window*dbsval_ (&t,&kx,tknot,&nf,bcoefim);
      fprintf(deb,"%lf %lf \n",workin[i][0], workin[i][1]);
    }

    //centers t=0 at the zero frequency position of the fft input vector
    center_fft(workin,NTG);
    //--- do forward fft
    // FFTW_BACKWARD -> sign in the exponent = 1.
    p = fftw_plan_dft_1d(NTG,workin,workout,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p);  
    //shifts the fft result back to the center of the array
    center_fft(workout,NTG);

    printf("\nFinal spectrum computed! \n incoming photon frequency: %lf a.u. \n\n",omega);
    

    printf("\n > RIXS spectrum as function of the emitted photon energy \n\n");
    printf("    E' (eV)          Re           Im \n");
    printf("------------------\n");
    for(i=0;i<nE;i++){
      val[0] = Ei + i*stepE;
      val[0] = omega - val[0] + shift;
      printf("% E % E % E \n",val[0]*27.2114,workout[i][0],workout[i][1]);
    }

    printf("\n\n > RIXS spectrum as function of energy loss \n\n");
    printf("    E - E' (eV)          Re           Im \n");
    printf("------------------\n");
    for(i=0;i<nE;i++){
      val[0] = Ei + i*stepE;
      val[0] = val[0] - shift;
      printf("% E % E % E \n",val[0]*27.2114,workout[i][0],workout[i][1]);
    } 

//------------------------------------------------------------------------
//------------------------------------------------------------------------

  }else if(type==2){
    //------- XAS Spectrum
    printf("\n << Running 2D+1D model: XAS cross-section >>\n");

    printf("\n << Starting Franck-Condon section >> \n\n");
    printf("number of vc states: %d ",nvc);
    FCGVc = malloc(nvc*sizeof(double));
    intens = malloc(nvc*sizeof(double));

    printf("reading <0|vc> FC vector from file %s \n",fcnam);

    fcf = fopen(fcnam,"r"); 
    printf("\n<0|vc>\n   ");
    for(i=0;i<nvc;i++)  printf("      |%d> ",i);
    printf("\n");
    printf("< %d| ",0);
    for(i=0;i<nvc;i++){
      fscanf(fcf,"%lf",&FCGVc[i]); 
      printf("% lf ", FCGVc[i]);
    }

    // reading the normalization factors from the wf generated with the fft 
    //fscanf(fcf,"%lf",&intens[0]);
    
    fclose(fcf);

    printf("\n");
    printf("core-excited state bending energy values:\n");
    for(i=0;i<nvc;i++) printf(" Evc%i = %E , ",i,Evc[i]);
    printf("\n");

    printf("\n> Franck-Condon Section finished \n");

    printf("\n << Starting Auto-Correlation Function section >> \n");

    printf("\n the correlation function will be read from file %s \n",fcornam);

    T = malloc(nf*sizeof(double));

    fcorrelRe = malloc(sizeof(double));
    fcorrelIm = malloc(sizeof(double));

    fcorrelRe[0] = malloc(nf*sizeof(double));
    fcorrelIm[0] = malloc(nf*sizeof(double));

    // read correlation function
    fcorr = fopen(fcornam,"r");
    for(ii=0;ii<nf;ii++){
      fscanf(fcorr,"%lf",&T[ii]);
      T[ii] = T[ii] * 41.3413745758e+0;
      fscanf(fcorr,"%lf %lf",&fcorrelRe[0][ii],&fcorrelIm[0][ii]);
    }
    fclose(fcorr);
    stept = (T[nf-1] - T[0])/(nf - 1.0);

    // total correlation function
    tfcorrelRe = malloc(nf*sizeof(double));
    tfcorrelIm = malloc(nf*sizeof(double));
    for(ii=0;ii<nf;ii++){
      tfcorrelRe[ii] = 0.0e+0;
      tfcorrelIm[ii] = 0.0e+0;
    }

    // eq. (85) from file

    for(j=0;j<nvc;j++){
      //FC = intens[0] * FCGVc[j]* FCGVc[j];
      FC = FCGVc[j]* FCGVc[j];
      for(ii=0;ii<nf;ii++){
	tfcorrelRe[ii] = tfcorrelRe[ii] + FC*(fcorrelRe[0][ii]*cos(Evc[i]*T[ii]) + fcorrelIm[0][ii]*sin(Evc[i]*T[ii]) );
	tfcorrelIm[ii] = tfcorrelIm[ii] + FC*(fcorrelIm[0][ii]*cos(Evc[i]*T[ii]) - fcorrelRe[0][ii]*sin(Evc[i]*T[ii]) );
	//debug below
	//tfcorrelRe[ii] = tfcorrelRe[ii] + FC*fcorrelRe[k][ii];
	//tfcorrelIm[ii] = tfcorrelIm[ii] + FC*fcorrelIm[k][ii];
      }
    }

    printf("Final Total Auto-Correlation Function \n");
    printf("------------------\n");
    for(ii=0;ii<nf;ii++){
      printf("%lf % lf % lf \n",T[ii],tfcorrelRe[ii],tfcorrelIm[ii]);
    }

    printf("\n> Auto-Correlation Section finished \n");

    printf("\n << Starting XAS Cross-Section section >> \n");

    //NTG = pow(2,twopow+1);
    NTG = pow(2,twopow);
    steptspl=(2.0*T[nf-1])/NTG;
    //stepE = 2*M_PI/(NTG*FSAU(steptspl));
    stepE = 2*M_PI/(NTG*steptspl);
    Ei = -NTG*stepE/(2.0E+0);
    nE =  NTG;

    printf("\nAuto-Correlation Function will be adjusted for Fourier Transform using splines \n");
    printf("\nFourier parameters \n");
    printf("npoints: 2^%d = %d\n",twopow,NTG);
    printf("time step: %lf fs \n",steptspl);
    printf("Energy step: %lf a.u. , final energy: %lf a.u. \n",stepE,-Ei);
    //steptspl=(tf - ti)/(NTG);

    //spline arrays
    bcoefre = malloc(nf*sizeof(double));
    bcoefim = malloc(nf*sizeof(double));
    tknot = malloc((nf+kx)*sizeof(double));
    //fft arrays
    workin = fftw_malloc(sizeof(fftw_complex) * NTG);
    workout = fftw_malloc(sizeof(fftw_complex) * NTG);

    dbsnak_ (&nf, T, &kx, tknot);
    dbsint_ (&nf,T,tfcorrelRe,&kx,tknot,bcoefre);
    dbsint_ (&nf,T,tfcorrelIm,&kx,tknot,bcoefim);

    //generate mirrored data
    for(i=0;i<NTG;i++){
      t = -T[nf-1] + i*steptspl;
      fprintf(deb,"%lf ",t);
      if(t < 0) s = -1.0e+0;
      else s = 1.0e+0;

      t = fabs(t);

      if(strncasecmp(windtype,".SGAUSS",7)==0){
	//TAUX = TMAX**2/(LOG(ONE/WP))**2 espec
	//taux = pow(T[nf-1],2)/(log(1.000/width)/log(M_E));
	taux = pow(T[nf-1],2)/pow(log(1.000/width),2);
	window = exp(-pow(t,2)/taux);
      }else if(strncasecmp(windtype,".EXPDEC",7)==0){
	window = exp(-width*t);
      }

      //debug
      //window=1.000;
      workin[i][0] =   window*dbsval_ (&t,&kx,tknot,&nf,bcoefre);
      workin[i][1] = s*window*dbsval_ (&t,&kx,tknot,&nf,bcoefim);
      fprintf(deb,"%lf %lf \n",workin[i][0], workin[i][1]);
    }

    //centers t=0 at the zero frequency position of the fft input vector
    center_fft(workin,NTG);
    //--- do forward fft
    // FFTW_BACKWARD -> sign in the exponent = 1.
    p = fftw_plan_dft_1d(NTG,workin,workout,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p);  
    //shifts the fft result back to the center of the array
    center_fft(workout,NTG);

    printf("\nFinal spectrum computed! \n");

    printf("\n > XAS spectrum as function of detuning from resonance \n\n");
    printf("    Detuning (eV)          Re           Im \n");
    printf("------------------\n");
    for(i=0;i<nE;i++){
      val[0] = Ei + i*stepE - Delta;
      printf("% E % E % E \n",val[0]*27.2114,workout[i][0],workout[i][1]);
    }

    printf("\n > XAS spectrum as function od photon frequency\n\n");
    printf("    E (eV)          Re           Im \n");
    printf("------------------\n");
    for(i=0;i<nE;i++){
      val[0] = Ei + i*stepE - Delta + shift;
      printf("% E % E % E \n",val[0]*27.2114,workout[i][0],workout[i][1]);
    }

//------------------------------------------------------------------------
//------------------------------------------------------------------------
  }else if(type==3){
    printf("\n\n<< REXS cross section with Self-Absorption term calculation >>\n\n");
    printf("original REXS cross section will be read from file %s, npoints = %d. \n",rexfnam,nrexs);
    printf("XAS cross section will be read from file %s, npoints = %d. \n\n",xasfnam,nxas);

    omega = omega * 27.2114;
    printf("\n incident photon energy: %lf a.u., %lf eV \n\n",omega/27.2114,omega);

    rexsfile=fopen(rexfnam,"r");
    xasfile=fopen(xasfnam,"r");

    xas_omega = malloc(nxas * sizeof(double));
    xas_cross = malloc(nxas * sizeof(double));
    xas_bcoef = malloc(nxas * sizeof(double));
    xas_knot  = malloc((nxas + kx) * sizeof(double));

    rexs_omegap = malloc(nrexs * sizeof(double));
    rexs_cross = malloc(nrexs * sizeof(double));

    rexs_cross_sa = malloc(nrexs * sizeof(double));

    //reading XAS cross section from file
    for(i=0;i<nxas;i++) fscanf(xasfile,"%lf %lf",&xas_omega[i],&xas_cross[i]);
    fclose(xasfile);

    //spline interpolation of XAS cross section
    dbsnak_ (&nxas, xas_omega, &kx, xas_knot);
    dbsint_ (&nxas,xas_omega,xas_cross,&kx,xas_knot,xas_bcoef);

    //reading REXS cross section from file
    for(i=0;i<nrexs;i++) fscanf(rexsfile,"%lf %lf",&rexs_omegap[i],&rexs_cross[i]);
    fclose(rexsfile);


    //final REXS cross section including self-absorption
    printf("\n > REXS-SA spectrum as function of emitted photon energy \n\n");
    printf("    E' (eV)       Re           \n");
    printf("------------------\n");
    for(i=0;i<nrexs;i++){
      if(rexs_omegap[i] > xas_omega[0] && rexs_omegap[i] < xas_omega[nxas-1]){
	xas_cross_omp = dbsval_ (&rexs_omegap[i],&kx,xas_knot,&nxas,xas_bcoef);
	xas_cross_om  = dbsval_ (&omega,&kx,xas_knot,&nxas,xas_bcoef);

	//to avoid numerical instabilities due to small negative values that may occur in FT
	xas_cross_omp = fabs(xas_cross_omp);
	xas_cross_om = fabs(xas_cross_om);
	rexs_cross[i] = fabs(rexs_cross[i]);

	rexs_cross_sa[i] = rexs_cross[i] / ( 1.0e+0 +  xas_cross_omp/xas_cross_om  );
	printf("% E % E \n",rexs_omegap[i],rexs_cross_sa[i]);
      }
    }

    printf("\n\n > REXS-SA spectrum as function of energy loss \n\n");
    printf("    E - E' (eV)       Re           \n");
    printf("------------------\n");
    for(i=0;i<nrexs;i++){
      //printf("%lf > %lf , %lf > %lf \n",rexs_omegap[i], xas_omega[0], rexs_omegap[i], xas_omega[nxas-1]);
      if(rexs_omegap[i] > xas_omega[0] && rexs_omegap[i] < xas_omega[nxas-1]){
	printf("% E % E \n",omega - rexs_omegap[i],rexs_cross_sa[i]);
      }
    } 

    free(xas_omega);free(xas_cross);free(xas_bcoef);free(xas_knot);
    free(rexs_cross);free(rexs_cross_sa);free(rexs_omegap);

  }

  printf("\n\n# End of Calculation!\n");
  printf("# Correlation terminated successfully!\n");
  return 0;
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
