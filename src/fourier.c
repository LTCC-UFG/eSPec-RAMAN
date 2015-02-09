#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>
#include<string.h>

#include "splinesurf.h"
#define one 1.0E+0

/*
  gendata_Psic:

  this routine generates the Psic core excited wave packet. According to

  Psic = Psic(t) , t>=0
  Psic = Psic(|t|)*, t <0

 */

int gendata_Psic(fftw_complex *work, double *bcoefre, double *bcoefim, double *xknot, double *yknot, double *tknot, int kx, int *np, int nf, double x, double y, double ti,double tf, double steptspl, int NTG, double width, char *windtype, char *dim){
  int i,j,nt,m;
  double s;
  double t,val,window,taux;
  //FILE *wind=fopen("window.dat","w");


  for(i=0;i<NTG;i++){
    t = -tf + i*steptspl;
    
    if(t < 0) s = -1.0e+0;
    else s = 1.0e+0;

    t = fabs(t);

    if(strncasecmp(windtype,".SGAUSS",7)==0){
      taux = pow(tf,2)/(log(one/width)/log(M_E));
      window = exp(-pow(t,2)/taux);
    }else if(strncasecmp(windtype,".EXPDEC",7)==0) window = exp(-width*t);

    if(strncasecmp(dim,".2D",3)==0){
      work[i][0] = window*dbs3vl_ (&y,&x,&t,&kx,&kx,&kx,yknot,xknot,tknot,&np[0],&np[1],&nf,bcoefre);
      work[i][1] = s*window*dbs3vl_ (&y,&x,&t,&kx,&kx,&kx,yknot,xknot,tknot,&np[0],&np[1],&nf,bcoefim);
    }else if(strncasecmp(dim,".1D",3)==0){
      work[i][0] = window*dbs2vl_ (&y,&t,&kx,&kx,yknot,tknot,&np[0],&nf,bcoefre);
      work[i][1] = s*window*dbs2vl_ (&y,&t,&kx,&kx,yknot,tknot,&np[0],&nf,bcoefim);
    }

  }
}

//----------------------------------------------------------

/*
  gendata_Psic_bar:

  this routine generates the Psic_bar core excited wave packet. According to

  Psic_bar =  Psic(t) , t>=0
  Psic_bar = -Psic(|t|)*, t <0

 */

int gendata_Psic_bar(fftw_complex *work, double *bcoefre, double *bcoefim, double *xknot, double *yknot, double *tknot, int kx, int *np, int nf, double x, double y, double ti,double tf, double steptspl, int NTG, double width, char *windtype, char *dim){
  int i,j,nt,m;
  double s;
  double t,val,window,taux;
  //FILE *wind=fopen("window.dat","w");


  for(i=0;i<NTG;i++){
    t = -tf + i*steptspl;
    
    if(t < 0) s = -1.0e+0;
    else s = 1.0e+0;

    t = fabs(t);

    if(strncasecmp(windtype,".SGAUSS",7)==0){
      taux = pow(tf,2)/(log(one/width)/log(M_E));
      window = exp(-pow(t,2)/taux);
    }else if(strncasecmp(windtype,".EXPDEC",7)==0) window = exp(-width*t);

    if(strncasecmp(dim,".2D",3)==0){
      work[i][0] = s*window*dbs3vl_ (&y,&x,&t,&kx,&kx,&kx,yknot,xknot,tknot,&np[0],&np[1],&nf,bcoefre);
      work[i][1] = window*dbs3vl_ (&y,&x,&t,&kx,&kx,&kx,yknot,xknot,tknot,&np[0],&np[1],&nf,bcoefim);
    }else if(strncasecmp(dim,".1D",3)==0){
      work[i][0] = s*window*dbs2vl_ (&y,&t,&kx,&kx,yknot,tknot,&np[0],&nf,bcoefre);
      work[i][1] = window*dbs2vl_ (&y,&t,&kx,&kx,yknot,tknot,&np[0],&nf,bcoefim);
    }

  }
}





/*
 shifts zero frequency to the center of the array
*/

int center_fft(fftw_complex *out,int N){
  int i;
  double work;

  //centering the fourier transform
  for(i=0;i<N/2;i++){
    work= out[N/2 +i][0];
    out[N/2 +i][0] = out[i][0];
    out[i][0] = work;

    work= out[N/2 +i][1];
    out[N/2 +i][1] = out[i][1];
    out[i][1] = work;  
  }
}



/*
  run_all_fft

  This routine converts a wavepacket propagation from |f(x,y,t)> to |f(x,y,E)> , from a previously calculated cubic spline representation, using the fftw library routines

  main input: bcoefre, bcoefim. real and complex part spline matrices of |f(x,y,t)>
  main output: WPERe, WPEIm. real and complex part vectors of |f(x,y,E)>

 */
int run_all_fft(double *bcoefre,double *bcoefim, double *X, double *Y,double *xknot, double *yknot, double *tknot, int kx,double ti,double tf,double steptspl,int *np,int nf,int NTG, double width,double *WPERe,double *WPEIm, char *windtype, char *dim){
  int i,j,l,ll,k,nE;
  double x,y,E,ke,norm,stepxspl;
  fftw_complex *workin,*workout, *workin_b,*workout_b;
  fftw_plan p,p_b;
  double Ei,stepE;
  FILE *filin=fopen("fft_check_in.dat","w");
  FILE *filout=fopen("fft_check_out.dat","w");

  stepE = 2*M_PI/(NTG*steptspl);
  Ei = -NTG*stepE/(2.0E+0);
  nE =  NTG; //NTG/2
  printf("\nEnergy step: %E \n\n",stepE);

  //--- allocating work vectors
  workin = fftw_malloc(sizeof(fftw_complex) * NTG);
  workout = fftw_malloc(sizeof(fftw_complex) * NTG);
  workin_b = fftw_malloc(sizeof(fftw_complex) * NTG);
  workout_b = fftw_malloc(sizeof(fftw_complex) * NTG);


  for(i=0;i<np[1];i++){
    for(j=0;j<np[0];j++){
      //generates data set from spline coeff.
      gendata_Psic(workin,bcoefre,bcoefim,xknot,yknot,tknot,kx,np,nf,X[i],Y[j],ti,tf,steptspl,NTG,width,windtype,dim);
      gendata_Psic_bar(workin_b,bcoefre,bcoefim,xknot,yknot,tknot,kx,np,nf,X[i],Y[j],ti,tf,steptspl,NTG,width,windtype,dim);
      //centers t=0 at the zero frequency position of the fft input vector(first) due to periodicity requirement of fft
      center_fft(workin,NTG);
      center_fft(workin_b,NTG);

      //--- do forward fft
      // FFTW_BACKWARD -> sign in the exponent = 1.
      p = fftw_plan_dft_1d(NTG,workin,workout,FFTW_BACKWARD,FFTW_ESTIMATE);
      p_b = fftw_plan_dft_1d(NTG,workin_b,workout_b,FFTW_BACKWARD,FFTW_ESTIMATE);

      fftw_execute(p);  
      fftw_execute(p_b);  

      //shifts the fft result back to the center of the array
      center_fft(workout,NTG);
      center_fft(workout_b,NTG);

      //debug ---
      if(j==26){
	//center_fft(workin,NTG);
	fprintf(filin,"# %E %E\n",X[i],Y[j]);
	for(k=0;k<NTG;k++)fprintf(filin,"%E %E %E %E %E\n",-tf + k*steptspl,workin[k][0],workin[k][1],workin_b[k][0],workin_b[k][1]);
	fprintf(filout,"# %E %E\n",X[i],Y[j]);
	//for(k=0;k<NTG;k++)fprintf(filout,"%E %E %E\n",Ei + k*stepE,steptspl*workout[k][0],steptspl*workout[k][1]);
	for(k=0;k<NTG;k++)fprintf(filout,"%E %E %E %E %E\n",Ei + k*stepE,steptspl*workout[k][0],steptspl*workout[k][1],steptspl*workout_b[k][0],steptspl*workout_b[k][1]);
      }
      //----------
      
      //copies fft results to permanent array
      for(l=0;l<NTG;l++){ 
	//result of discrete FT must be multiplied by the time step (see num. recipes)
	ll = j + i*np[0] + l*np[1]*np[0];
	//printf("%d %d %d -> %d \n",j,i,k,ll);
	//fprintf(fil," %E %E %E\n",Ei + l*stepE,steptspl*workout[l][0],steptspl*workout[l][1]);
	//WPERe[ll] = (1.0/sqrt(2*M_PI))*steptspl*workout[l][0];
	//WPEIm[ll] = (1.0/sqrt(2*M_PI))*steptspl*workout[l][1];
	WPERe[ll] = steptspl*workout[l][0] + steptspl*workout_b[l][0];
	WPEIm[ll] = steptspl*workout[l][1] + steptspl*workout_b[l][1];
      }
      //fprintf(fil,"\n");
    }
  }
  //----------


  fftw_destroy_plan(p);
  fftw_destroy_plan(p_b);
  fftw_free(workin); 
  fftw_free(workout);
  fftw_free(workin_b); 
  fftw_free(workout_b);
  return 0;
}


/*
 gausswp

 this routine generates the initial state spatial distribution

 */

double gausswp(double x, double x0, double hwhm, double norm){
  double gauss,X,ln_2;
  X = x - x0;
  ln_2 = log(2);
  gauss = norm*exp(-ln_2*pow(X/hwhm,2));
  return gauss;
}

/*
 gausswp_k

 this routine generates the initial state momentum distribution

 */

double gausswp_k(double k, double k0, double hwhm){
  double gauss,norm,K,ln2;
  K = k - k0;
  // M_LN2 = 6.931471805599453E-01
  ln2 = 0.69314718e+0; //value used in eSPec
  //norm = sqrt(sqrt( hwhm*hwhm/(2.0*M_PI*M_LN2)));
  //gauss = norm*exp(-pow(hwhm*K,2)/(4.0*M_LN2));
  norm = sqrt(sqrt( hwhm*hwhm/(2.0*M_PI*ln2)));
  gauss = norm*exp(-pow(hwhm*K,2)/(4.0*ln2));
  return gauss;
}


/*
  print_bcoef

  this routine prints the spline coefficients into a file for a later calculation.

 */

int print_bcoef(char *jobnam,double *X, double *Y, double *E,double *bcoefre, double *bcoefim, int *np,int nE){
  int i,nt;
  FILE *out=fopen("fft_spline.bcoef","w");
  nt = np[0]*np[1]*nE;

  fprintf(out,"%s\n",jobnam);
  fprintf(out,"%d %d %d\n",np[0],np[1],nE);

  for(i=0;i<np[0];i++)fprintf(out,"%.15E ",Y[i]);
  for(i=0;i<np[1];i++)fprintf(out,"%.15E ",X[i]);
  for(i=0;i<nE;i++)fprintf(out,"%.15E ",E[i]);

  for(i=0;i<nt;i++){
    fprintf(out,"%.15E %.15E ",bcoefre[i],bcoefim[i]);
  }

  fclose(out);

  return 0;

}


/*
  print_bcoef

  this routine reads the spline coefficients from a file from a previous calculation.

 */

int read_bcoef(char *jobnam,double *X, double *Y, double *E, double *bcoefre, double *bcoefim, int *np,int nE){
  int i,nt;
  int npr[3],nEr;
  char jobnamr[50];
  FILE *out=fopen("fft_spline.bcoef","r");
  nt = np[0]*np[1]*nE;

  fscanf(out,"%s\n",jobnamr);
  fscanf(out,"%d %d %d\n",&npr[0],&npr[1],&nEr);

  if(np[0] != npr[0]){
    printf("Error! discrepancies found between the inputed number of points and the bcoef file ny = %d != %d \n you might not want to use this file.",np[0],npr[0]);
    return 66;
  }else if(np[0] != npr[0]){
    printf("Error! discrepancies found between the inputed number of points and the bcoef file nx = %d != %d \n you might not want to use this file.",np[1],npr[1]);
    return 66;
  }else if(np[0] != npr[0]){
    printf("Error! discrepancies found between the inputed number of points and the bcoef file nx = %d != %d \n you might not want to use this file.",nE,nEr);
    return 66;
  }

  for(i=0;i<np[0];i++)fscanf(out,"%lf",&Y[i]);
  for(i=0;i<np[1];i++)fscanf(out,"%lf",&X[i]);
  for(i=0;i<nE;i++)fscanf(out,"%lf",&E[i]);

  for(i=0;i<nt;i++){
    fscanf(out,"%lf %lf",&bcoefre[i],&bcoefim[i]);
  }

  fclose(out);

  return 0;

}


/*
  gendata_function
  routine to test the fourier transform
*/
int gendata_function(fftw_complex *work,double tf,double steptspl,int NTG,double width){
  int i;
  double t,t0,window,taux,w0;
  //FILE *wind=fopen("window.dat","w");
  t0 = 0.0*tf/2;//0.0e+0;//tf/2;
  w0 = 5.0e+0;
  width = 1.0e-2;
  taux = pow(tf/(log(one/width)/log(M_E)),2);
  for(i=0;i<NTG;i++){
    t = -tf + i*steptspl;
    window = exp(-pow(t,2)/taux);
    //fprintf(wind,"%E %E\n",t,window);
    //window = 1.0e+0;
    //work[i][0]=window*(cos(2*t)+cos(t)+5.0*cos(5*t));
    work[i][0]=exp(-M_PI*pow(t-t0,2))*cos(w0*t);
    work[i][1]=exp(-M_PI*pow(t-t0,2))*sin(w0*t); //0.0;//window*sin(t);//0.0E+0;
  }
  return 0;
}


/*

OLD BACKUP


/*
  gendata_complex:

  this routine generates a full complex data set from its [0,tf] spline representation, then mirrors it back to [-tf,0].
  it is also applied a supergaussian window to the data set, according to the width input value.
  the main purpose is to prepare the data for fourier transform.

int gendata_complex(fftw_complex *work, double *bcoefre, double *bcoefim, double *xknot, double *yknot, double *tknot, int nxr, int nyr, int nf, double x, double y, double tf, double steptspl, int NTG, double width){
  int i,j,nt,korder;
  double t,val,window,taux;
  korder = 3.0;
  
  taux = pow(tf,2)/(log(one/width)/log(M_E));

  for(i=0;i<NTG;i++){
    t = -tf + i*steptspl;
    window = exp(-pow(t,2)/taux);
    window = 1.0000;
    if(t>=0.00e+0){
      t = fabs(t);
      work[i][0] = window*dbs3vl_ (&y,&x,&t,&korder,&korder,&korder,yknot,xknot,tknot,&nyr,&nxr,&nf,bcoefre);
      work[i][1] = window*dbs3vl_ (&y,&x,&t,&korder,&korder,&korder,yknot,xknot,tknot,&nyr,&nxr,&nf,bcoefim);
    }else{
      t = fabs(t);
      work[i][0] =  window*dbs3vl_ (&y,&x,&t,&korder,&korder,&korder,yknot,xknot,tknot,&nyr,&nxr,&nf,bcoefre);
      work[i][1] = -window*dbs3vl_ (&y,&x,&t,&korder,&korder,&korder,yknot,xknot,tknot,&nyr,&nxr,&nf,bcoefim);
    }
  }
}
 */

void checkspl1d(double *Y,double *T,double *bcoefre, double *bcoefim, double *xknot,double *yknot,double *tknot,int *np,int nf,int kx){
  int i,j,k;
  double val[2];
  FILE *spl;
  char fnam[20];

  for(k=0;k<nf;k++){
    sprintf(fnam,"spl_%d.dat",k+1);
    printf("opening file %s \n",fnam);
    spl = fopen(fnam,"w");
    for(j=0;j<np[0];j++){
      val[0] = dbs2vl_ (&Y[j],&T[k],&kx,&kx,yknot,tknot,&np[0],&nf,bcoefre);
      val[1] = dbs2vl_ (&Y[j],&T[k],&kx,&kx,yknot,tknot,&np[0],&nf,bcoefim);
      fprintf(spl,"%E %E %E \n",Y[j],val[0],val[1]);
    }
    fclose(spl);
  }
  
}
