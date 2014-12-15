#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

int rdinput(FILE *arq,char *dim,int *np,char *file,int *stfil,int *endfil, double *ti, double *tf, double *pti, double *ptf, double *pstept, double *m,char *potfile, int *nf,int *twopow, double *width, double *Ef, int *type, double *crosst){
  int i,j,k,spl;
  int NXG,NYG,jcheck1,jcheck2;
  char workk[50],coment[150],typenam[20];

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
	  if(strncasecmp(typenam,"non-reac",8)==0)*type=0;
	  else if(strncasecmp(typenam,"reac",4)==0)*type=1;
	  else{
	    printf("invalid run type entered: %s",typenam);
	    return 666;
	  }
	}
	if(strncasecmp(workk,"dimension",6)==0){
	  fscanf(arq,"%s", dim);
	  fscanf(arq,"%d",&np[0]);
	  if(strncasecmp(dim,".2D",3)==0) fscanf(arq,"%d",&np[1]);
	}
	if(strncasecmp(workk,"Mass",4)==0){
	  fscanf(arq, "%lf", &m[0]);
	  if(strncasecmp(dim,".2D",3)==0)fscanf(arq, "%lf", &m[1]);
	  if(strncasecmp(dim,".2DCT",5)==0)fscanf(arq, "%lf", &m[2]);
	}
	if(strncasecmp(workk,"filename",8)==0) fscanf(arq, "%s", file);
	if(strncasecmp(workk,"nfiles",6)==0){
	  fscanf(arq, "%d %d",stfil, endfil);
	  *nf=(*endfil-*stfil) + 1 ;
	}
	if(strncasecmp(workk,"timeinterval",12)==0) fscanf(arq, "%lf %lf", ti, tf);
      }
    }
    //FLUX GROUP INPUT 
    if(strncasecmp(workk,"*Propagation",12)==0){
      printf("inside *propag \n");
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
	if(strncasecmp(workk,"crossterm",9)==0)fscanf(arq, "%lf", crosst);
	if(strncasecmp(workk,"finalenergy",11)==0)fscanf(arq, "%lf", Ef);
	if(strncasecmp(workk,"Fourier",7)==0)fscanf(arq,"%d",twopow);
	if(strncasecmp(workk,"Window",6)==0)fscanf(arq,"%lf",width);
	
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
