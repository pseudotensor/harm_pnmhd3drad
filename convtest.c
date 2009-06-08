// compile using:
// gcc -Wall -lm -O3 -o convtest convtest.c
// to run in base dir with c code:
// mkdir ./conv
// ./convtest > ./conv/0_o.out
// 1-D L1-error norm convergence and rate of convergence calculation
// Currently works with sod shock and advection test

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define DATEXT ".dat"
#define PAREXT ".par"

#define NUMSCA 3
#define NUMVEC 2
//#define NUMTOT (NUMSCA+NUMVEC*3)
#define NUMTOT (NUMSCA+(NUMVEC-1)*3)

#define NUMRUNS 6 // number of runs to do
#define STARTNX 5 // what grid size to start at
#define ENDNX  25 // what grid size to end at
#define STEPNX (ENDNX-STARTNX+1)/NUMRUNS // the change in grid size per run
#define SMOOTHSIZE 5 // number of samples to take per rate of convergence point from raw differences, should be an odd number 1,3,5,7,etc.  Only necessary when doing small steps(1-5) in grid size

#define MAXFILENAME 400

#define ZER 0x30 // definition of ascii 0 in hex
#define DIGILEN 10 //  how many digits max allowed by itoa converter at bottom of this code
#define TRUEDIGITS 4 // must agree with above NX(esp ENDNX). #digits(ENDNX)<=this amount
#define TIMEDIGITS 4 // how many digits incoming dumpxxxxxx.dat has where x=this number.
#define ROFILE 0 // 0=do full run over again, computing solutions, etc. 1= already  have solutions, just recompute rate of convergence using new parms


#define N2 64 // needs to be same as in global.h

void itoa(int x,char*p,int lead0);

int main(void)
{

  int i,j,k,l,m,n;
  char change1[MAXFILENAME]="#define EXTERNALRUN 0";
  char part1[MAXFILENAME]="sed \'s/#define N1 ";
  char part2[DIGILEN+1];
  char part3[MAXFILENAME]="/#define N1 ";
  char part4[DIGILEN+1];
  char part5[MAXFILENAME]="/\' global.h > global.h.1";
  char csed[MAXFILENAME];

  char fnames[3][MAXFILENAME];

  char tempc1[MAXFILENAME],tempc2[MAXFILENAME],tempc3[MAXFILENAME];

  double s[NUMTOT][ENDNX*N2];
  double anal[NUMTOT][ENDNX*N2];

  double sum[NUMTOT][NUMRUNS];
  double roc[NUMTOT][NUMRUNS];
  double sums[NUMTOT][NUMRUNS];

  int gsize[NUMRUNS];

  double dum;
  int dumd;
  int dumi;
  double h1,h2;

  FILE * out[3];
  FILE * tempout;

  FILE * in;

  char tch;
  int nn;

  int nx,ny;
  int extra;
  int inan;

  char temps[200];



  system("mkdir ./conv");
  strcpy(fnames[0],"./conv/conv");
  strcat(fnames[0],DATEXT);
  strcpy(fnames[1],"./conv/roconv");
  strcat(fnames[1],DATEXT);
  strcpy(fnames[2],"./conv/sroconv");
  strcat(fnames[2],DATEXT);

#if(ROFILE==0)

  strcpy(temps,"rm ./conv/*");
  strcat(temps,DATEXT);
  system(temps);

  // estimate minutes for anal runs

  dum=0;
  for(i=STARTNX;i<=ENDNX;i+=STEPNX){
    dum=5.817E-6*i*i+dum;
    if( ((ENDNX-i)<STEPNX)&&( (ENDNX-i)!=0)){
      i=ENDNX-STEPNX;
    }
  }
  dum=dum+8.0/60.0*NUMRUNS; // account for compile time

  strcpy(temps,"./conv/time");
  strcat(temps,DATEXT);
  if((tempout=fopen(temps,"wt"))==NULL){
    printf("cannot write %s file\n",temps);
    exit(1);
  }
  fprintf(tempout,"Convergence Test running:\n");
  fprintf(tempout,"%d to %d\n",STARTNX,ENDNX);
  fprintf(tempout,"%d runs taking about %21.15g minutes on dedicated P2-450 1 CPU\n",NUMRUNS,dum);
  fclose(tempout);
  j=0;
  
  extra=0;
  m=0; //run number
  for(i=STARTNX;i<=ENDNX;i+=STEPNX){
    gsize[m]=i;

    strcpy(csed,part1);
    itoa(j,part2,0);    
    strcat(csed,part2);
    strcat(csed,part3);
    itoa(i,part4,0);    
    strcat(csed,part4);
    strcat(csed,part5);


    
    printf("%d\n",i);
    printf("%s\n",csed);
    system(csed);
    //    system("touch global.h");
    system("cp global.h.1 global.h");
    system("sed 's/#define EXTERNALRUN 0/#define EXTERNALRUN 1/' global.h > global.h.1");
    system("cp global.h.1 global.h");
    system("rm *.o");
    system("make -j ");
    system("rm ./bin/*.dat*");
    system("rm ./bin/i/*.dat*");
    
    strcpy(temps,"nice -10 ./bin/twod > ./bin/0_o.out");
    strcat(temps,DATEXT);
    system(temps);

    // determine last file name
    strcpy(temps,"ls -al ./bin/dump*");
    strcat(temps,DATEXT);
    strcat(temps," > ./conv/list");
    strcat(temps,DATEXT);
    system(temps);

    strcpy(temps,"./conv/list");
    strcat(temps,DATEXT);
    if((in=fopen(temps,"rt"))==NULL){
      printf("cannot read %s file\n",temps);
      exit(1);
    }
    nn=0;
    while(!feof(in)){
      tch=fgetc(in);
      if(tch=='\n') nn++;
    }
    fclose(in);
    strcpy(temps,"rm ./conv/list");
    strcat(temps,DATEXT);
    system(temps);

    strcpy(tempc1,"./bin/dump");
    itoa(nn-1,tempc2,TIMEDIGITS);
    strcat(tempc1,tempc2);
    strcat(tempc1,DATEXT);


    
    if((in=fopen(tempc1,"rt"))==NULL){
      printf("cannot read dump file\n");
      exit(1);
    }

    // copy final file for this run into safe place for later analysis
    strcpy(tempc3,"cp ");
    strcat(tempc3,tempc1);

    strcat(tempc3," ./conv/");
    
    itoa(i,tempc1,TRUEDIGITS);
    strcat(tempc3,tempc1);
    strcat(tempc3,"dump");
    strcat(tempc3,tempc2);
    strcat(tempc3,DATEXT);
    system(tempc3);

    while((tch=fgetc(in))!='\n'); // skip to next line(skip first comment line)
    while((tch=fgetc(in))!='\n'); // skip to next line(unneeded data line)
    while((tch=fgetc(in))!='\n'); // skip to next line(skip first comment line)
    fscanf(in," %d",&nx) ;
    while((tch=fgetc(in))!='\n'); // skip to end of line
    while((tch=fgetc(in))!='\n'); // skip to next line(skip comment line)

    if(i!=nx){
      printf("before s: nx does not equal i...nx: %d i: %d\n",nx,i);
      exit(1);
    }
    for(k=0;k<i*N2;k++){
      fscanf(in,"%lf %lf %lf", &s[0][k],&s[1][k],&s[2][k]) ;
      
      fscanf(in,"%lf %lf %lf\n", &s[3][k],&s[4][k],&s[5][k]);
      //      fscanf(in,"%lf %lf %lf\n", &s[6][k],&s[7][k],&s[8][k]);
    }
    fclose(in);

    strcpy(tempc1,"./bin/adump");
    strcat(tempc1,tempc2);
    strcat(tempc1,DATEXT);


    if((in=fopen(tempc1,"rt"))==NULL){
      printf("cannot read analytic dump file\n");
      exit(1);
    }

    // copy final file for this run into safe place for later analysis
    strcpy(tempc3,"cp ");
    strcat(tempc3,tempc1);

    strcat(tempc3," ./conv/");
    itoa(i,tempc1,TRUEDIGITS);
    strcat(tempc3,tempc1);
    strcat(tempc3,"adump");
    strcat(tempc3,tempc2);
    strcat(tempc3,DATEXT);
    system(tempc3);

    while((tch=fgetc(in))!='\n'); // skip to next line(skip first comment line)
    while((tch=fgetc(in))!='\n'); // skip to next line(unneeded data line)
    while((tch=fgetc(in))!='\n'); // skip to next line(skip first comment line)
    fscanf(in," %d",&nx) ;
    while((tch=fgetc(in))!='\n'); // skip to end of line
    while((tch=fgetc(in))!='\n'); // skip to next line(skip comment line)

    if(i!=nx){
      printf("before anal: nx does not equal i...nx: %d i: %d\n",nx,i);
      exit(1);
    }
    for(k=0;k<i*N2;k++){

      fscanf(in,"%lf %lf %lf", &anal[0][k],&anal[1][k],&anal[2][k]) ;
      
      fscanf(in,"%lf %lf %lf\n", &anal[3][k],&anal[4][k],&anal[5][k]);
      //fscanf(in,"%lf %lf %lf\n", &anal[6][k],&anal[7][k],&anal[8][k]);
    }    
    fclose(in);

    if((out[0]=fopen(fnames[0],"at"))==NULL){
      printf("cannot open %s file\n",fnames[0]);
      exit(1);
    }

    fprintf(out[0],"%10d",i);

    
    for(l=0;l<NUMTOT;l++){
      sum[l][m]=0;
      for(k=0;k<i*N2;k++){
	sum[l][m]=sum[l][m]+fabs(s[l][k]-anal[l][k]);
      }
      sum[l][m]=sum[l][m]/(i*(double)(N2));
      //check if nan

      sprintf(tempc1,"%3f",sum[l][m]);
      if(tempc1[0]=='n'){
	printf("sum is nan\n");
	sum[l][m]=0;
	for(k=0;k<i*N2;k++){
	  sum[l][m]=sum[l][m]+fabs(s[l][k]-anal[l][k]);
	  printf("i: %d i*N2: %d sum: %21.15g s[%d][%d]: %21.15g anal[%d][%d]: %21.15g\n",i,i*N2,sum[l][m],l,k,s[l][k],l,k,anal[l][k]);
	}
	sum[l][m]=sum[l][m]/(i*(double)(N2));
	printf("i: %d finalsum: %21.15g\n",i,sum[l][m]);
	//	exit(1);

      }
      //end check if nan
      
      fprintf(out[0]," %21.15g",sum[l][m]);      
    }
    fprintf(out[0],"\n");
    fclose(out[0]);
    m++;
    j=i;

    if( ((ENDNX-i)<STEPNX)&&( (ENDNX-i)!=0)){
      i=ENDNX-STEPNX;
      extra=1;
    }
  } // done all runs

#elif(ROFILE==1)
  
  if((in=fopen(fnames[0],"rt"))==NULL){
    printf("cannot open %s file\n",fnames[0]);
    exit(1);
  }
  for(m=0;m<NUMRUNS;m++){
    fscanf(in,"%d ",&dumi);

    for(l=0;l<NUMTOT;l++){
      fscanf(in,"%lf",&sum[l][m]);
      if(l==NUMTOT-1) fscanf(in,"\n");
      else fscanf(in," ");
    }
  }
  fclose(in);

  j=0;
  extra=0;
  m=0; //immitate run number

  for(i=STARTNX;i<=ENDNX;i+=STEPNX){
    gsize[m]=i;
    //null stuff
    m++;
    if( ((ENDNX-i)<STEPNX)&&( (ENDNX-i)!=0)){
      i=ENDNX-STEPNX;
      extra=1;
    }
  } //done immitate runs

#endif

  // compute rate of convergences

  if((out[1]=fopen(fnames[1],"wt"))==NULL){
    printf("cannot open %s file\n",fnames[1]);
    exit(1);
  }

  for(l=0;l<NUMTOT;l++){
    for(m=1;m<NUMRUNS-1;m++){
      //roc is the slope of the log-log of error vs resolution
      h1    = log10(1.E-20+gsize[m])-log10(1.E-20+gsize[m-1]);
      h2    = log10(1.E-20+gsize[m+1])-log10(1.E-20+gsize[m]);
      roc[l][m]=(log10(1.E-20+sum[l][m+1])*h2*h2
		 -log10(1.E-20+sum[l][m])*(h2*h2-h1*h1)
		 -log10(1.E-20+sum[l][m-1])*h1*h1)/(h1*h2*h2+h2*h1*h1);
    }
  }

  for(m=1;m<NUMRUNS-1;m++){
    //print roc
    inan=0;
    for(l=0;l<NUMTOT;l++){
      sprintf(tempc1,"%3f",roc[l][m]);
      if( (tempc1[0]=='i')||(tempc1[0]=='n') ){
	printf("roc[%d][%d] is inf or nan: %c\n",l,m,tempc1[0]);
	inan=1;
	roc[l][m]=-1000.0; // pick -1000.0 if nan or inf so other good data can be read properly
      }
    }
    //    if(inan==0){
    fprintf(out[1],"%10d ",gsize[m]);
    for(l=0;l<NUMTOT;l++){
      fprintf(out[1],"%21.15g",roc[l][m]);
      if(l==NUMTOT-1) fprintf(out[1],"\n");
      else fprintf(out[1]," ");
    }
  }

  fclose(out[1]);

  if((out[2]=fopen(fnames[2],"wt"))==NULL){
    printf("cannot open %s file\n",fnames[2]);
    exit(1);
  }


  
  // now produce smoothed version: really only valid for h1=h2
  for(l=0;l<NUMTOT;l++){
    for(m=SMOOTHSIZE/2;m<NUMRUNS-SMOOTHSIZE/2;m++){
      sums[l][m]=0;
      for(n=0;n<SMOOTHSIZE;n++){
	sums[l][m]=sums[l][m]+log10(1.E-20+sum[l][m+n-SMOOTHSIZE/2]);
      }
      sums[l][m]=sums[l][m]/((double)SMOOTHSIZE);
    }
  }


  // now compute smooth derivative
  for(l=0;l<NUMTOT;l++){
    for(m=SMOOTHSIZE/2+1;m<NUMRUNS-1-SMOOTHSIZE/2;m++){

      h1    = log10(1.E-20+gsize[m])-log10(1.E-20+gsize[m-1]);
      h2    = log10(1.E-20+gsize[m+1])-log10(1.E-20+gsize[m]);
      roc[l][m]=(sums[l][m+1]*h2*h2
		 -sums[l][m]*(h2*h2-h1*h1)
		 -sums[l][m-1]*h1*h1)/(h1*h2*h2+h2*h1*h1);
    }
  }



  for(m=SMOOTHSIZE/2+1;m<NUMRUNS-1-SMOOTHSIZE/2;m++){
    //print roc smooth version
    inan=0;
    for(l=0;l<NUMTOT;l++){
      sprintf(tempc1,"%3f",roc[l][m]);
      if( (tempc1[0]=='i')||(tempc1[0]=='n') ){
	printf("2:roc[%d][%d] is inf or nan: %c\n",l,m,tempc1[0]);
	inan=1;
      }
    }
    //    if(inan==0){
    fprintf(out[2],"%10d ",gsize[m]);
    for(l=0;l<NUMTOT;l++){
      fprintf(out[2],"%21.15g",roc[l][m]);
      if(l==NUMTOT-1) fprintf(out[2],"\n");
      else fprintf(out[2]," ");
    }
      //    }
  }
  
  fclose(out[2]);

  // copy setup and other files for last run as basis for entire run
  strcpy(temps,"cp ./bin/0_*");
  strcat(temps,DATEXT);
  strcat(temps," ./conv/");
  system(temps);
  strcpy(temps,"cp ./bin/*.m");
  strcat(temps," ./conv/");
  system(temps);

  return(0);
}





void itoa(int x,char*p, int lead0)
{
  // x holds number to be converted
  // p points to string
  // lead0 is number of true digits, to fill with 0s on left.  None if 0.

  int temp1;
  int i;
  int digits;

  temp1=x;
  if(temp1==0) digits=1;
  else{
    for(i=0;i<DIGILEN;i++){
      if(temp1==0){
	digits=i;
	break;
      }
      else temp1=temp1/10;
    }
  }
  
  for(i=0;i<digits;i++){
    temp1=x/10;
    temp1=x-temp1*10;
    p[digits-i-1]=ZER+temp1;
    x=x/10;
  }
  p[digits]='\0';

  if(lead0>digits){
    for(i=0;i<=digits;i++){ // copy string terminator too
      p[lead0-i]=p[digits-i];
    }
    for(i=0;i<lead0-digits;i++){
      p[i]='0';
    }
  }
  if(lead0>0){
    if(lead0<digits){
      printf("number of digits larger than specified number of true digits to fill with 0s when empty\n");
      exit(1);
    }
    if(lead0>DIGILEN){
      printf("number of true digits set larger than number allowed by memory allocation\n");
      exit(1);
    } 
  }
}



