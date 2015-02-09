#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>
#include<string.h>

#include "splinesurf.h"
#define one 1.0E+0


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
  }*/

