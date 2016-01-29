#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

int rdinput(FILE *arq,char *dim,int *np,char *file,int *stfil,int *endfil, double *ti, double *tf, double *pti, double *ptf, double *pstept, double *m,char *potfile, int *nf,int *twopow, double *width, int *nEf,double *Ef, int *type, double *crosst,char *windtype,char *jobnam, int *nfunc, char *funam,int *nvc,int *nvf,double *Evf,double *Evc,char *fcnam, char *fcornam,double *Delta,double *shift, char *rexfnam, char *xasfnam, int *nxas, int *nrexs, double *omega){
  int i,j,k,spl;
  int NXG,NYG,jcheck1,jcheck2;
  char workk[50],coment[150],typenam[20];

  k = 0;

  //rdinput
  while(fscanf(arq,"%s", workk)!=EOF){
    //MAIN INPUT GROUP
    //printf("3> %s \n", workk);
  check:
    if(strncasecmp(workk,"*Main",5)==0){   
      //printf("inside #Main \n");
      for(i=0;;i++){
	fscanf(arq,"%s", workk);
	while(workk[0]=='!'){
	  fgets(coment,90,arq);
	  fscanf(arq,"%s", workk);
	  if(feof(arq)) break;
	}
	if(feof(arq) || workk[0]=='*') goto check; //break;
	if(strncasecmp(workk,"runtype",7)==0){
	  fscanf(arq, "%s", typenam);
	  if(strncasecmp(typenam,"correl",6)==0)*type=0;
	  else if(strncasecmp(typenam,"spectrum",8)==0)*type=1;
	  else if(strncasecmp(typenam,"xas-spectrum",8)==0)*type=2;
	  else if(strncasecmp(typenam,"self-abs",8)==0)*type=3;
	  else{
	    printf("invalid run type entered: %s",typenam);
	    return 666;
	  }
	}
	if(strncasecmp(workk,"job",3)==0) fscanf(arq, "%s", jobnam);
	if(strncasecmp(workk,"dimension",6)==0){
	  fscanf(arq,"%s", dim);
	  fscanf(arq,"%d",&np[0]);
	  if(strncasecmp(dim,".2D",3)==0) fscanf(arq,"%d",&np[1]);
	}
	if(strncasecmp(workk,"Mass",4)==0){
	  fscanf(arq, "%lf", &m[0]);
	  if(strncasecmp(dim,".2D",3)==0)fscanf(arq, "%lf", &m[1]);
	}
	if(strncasecmp(workk,"corr_np",7)==0) fscanf(arq, "%d", nf);
	if(strncasecmp(workk,"filename",8)==0) fscanf(arq, "%s", file);
	if(strncasecmp(workk,"nfiles",6)==0){
	  fscanf(arq, "%d %d",stfil, endfil);
	  *nf=(*endfil-*stfil) + 1 ;
	}
	if(strncasecmp(workk,"timeinterval",12)==0) fscanf(arq, "%lf %lf", ti, tf);
      }
    }
    //FLUX GROUP INPUT 
    if(strncasecmp(workk,"*Correlation",12)==0){
      for(i=0;;i++){
	fscanf(arq,"%s", workk);
	//coment structure
	while(workk[0]=='!'){
	  fgets(coment,90,arq);
	  fscanf(arq,"%s", workk);
	  if(feof(arq)) break;
	}
	if(feof(arq) || workk[0]=='*') goto check;// break;

	if(strncasecmp(workk,"time",4)==0)fscanf(arq,"%lf %lf %lf",pti,ptf,pstept);
	if(strncasecmp(workk,"potential",9)==0)fscanf(arq,"%s",potfile);
	if(strncasecmp(workk,"wfunctions",10)==0)fscanf(arq, "%d %s", nfunc,funam);
      }
    }
    //2D + 1D GROUP INPUT 
    if(strncasecmp(workk,"*crosssection",13)==0){
      for(i=0;;i++){
	fscanf(arq,"%s", workk);
	//coment structure
	while(workk[0]=='!'){
	  fgets(coment,90,arq);
	  fscanf(arq,"%s", workk);
	  if(feof(arq)) break;
	}
	if(feof(arq) || workk[0]=='*') goto check;// break;

	if(strncasecmp(workk,"vf",2)==0)fscanf(arq, "%d",nvf);
	if(strncasecmp(workk,"Evf",3)==0){
	  for(i=0;i<*nvf;i++)fscanf(arq, "%lf",&Evf[i]);
	}
	if(strncasecmp(workk,"vc",2)==0)fscanf(arq, "%d",nvc);
	if(strncasecmp(workk,"Evc",3)==0){
	  for(i=0;i<*nvc;i++)fscanf(arq, "%lf",&Evc[i]);
	}
	if(strncasecmp(workk,"franckcondon",12)==0)fscanf(arq, "%s",fcnam);
	if(strncasecmp(workk,"fcorrel",7)==0)fscanf(arq, "%s",fcornam);
	if(strncasecmp(workk,"Delta",5)==0)fscanf(arq,"%lf",Delta);
	if(strncasecmp(workk,"shift",5)==0)fscanf(arq,"%lf",shift);
	if(strncasecmp(workk,"Fourier",7)==0)fscanf(arq,"%d",twopow);
	if(strncasecmp(workk,"Window",6)==0){
	  fscanf(arq,"%s",windtype);
	  fscanf(arq,"%lf",width);
	}
	if(strncasecmp(workk,"rexs-cs",7)==0){
	  fscanf(arq, "%d",nrexs);
	  fscanf(arq, "%s",rexfnam);
	}
	if(strncasecmp(workk,"xas-cs",6)==0){
	  fscanf(arq, "%d",nxas);
	  fscanf(arq, "%s",xasfnam);
	}
	if(strncasecmp(workk,"omega",5)==0)fscanf(arq,"%lf",omega);
	/*if(strncasecmp(workk,"Fourier",7)==0)fscanf(arq,"%d",twopow);
	if(strncasecmp(workk,"Window",6)==0){
	  fscanf(arq,"%s",windtype);
	  fscanf(arq,"%lf",width);
	  }*/	
      }
    }
    //PRINT INPUT GROUP
    /*if(strncasecmp(workk,"*Print",6)==0){
      // printf("inside #Print \n");
      for(i=0;;i++){
	fscanf(arq,"%s", &workk);
	while(workk[0]=='!'){
	  fgets(coment,90,arq);
	  fscanf(arq,"%s", &workk);
	  if(feof(arq)) break;
	}
	if(feof(arq) || workk[0]=='*') goto check;//break;
	if(strncasecmp(workk,"PRTEIGVC",8)==0) prtpot=0;
	if(strncasecmp(workk,"PRTPOT",8)==0) prteig=0;
	if(strncasecmp(workk,"PRTVEFF",11)==0)  prtjacpot=0;
	if(strncasecmp(workk,"Debug",5)==0)  debug=0;
	if(strncasecmp(workk,"DebugF",5)==0)  debugf=0;
      }
    }*/
  }
  //-----END OF INPUT
}
