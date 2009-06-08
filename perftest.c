// compile using:
// ccc -Wall -lm -O3 -o perftest perftest.c
// to run in dir with c code:
//  ./perftest > ./perf.out

// make sure only 1 space between  N1/N2 and the number

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "global.h"

// this controls number of cpus perf to be done on, where grid sizes are PER cpu
#define NUMCPUS 4 // >1 assumes USEMPI and USEGM are set appropriately already
// also assumes makefile.ncsa or makefile is set right
// if NUMCPUS>1, use the below machine or general if 0's
#define PHOTON 0
#define RAINMAN 0
#define WISEALPHA 0

// 64,8,512 ->8
// 32,16,512 -> 16
// 10,10,100 -> 10
// 63,4,128 -> 2
// 16,16,512 -> 31

#define NUMTRIALSX 63 // number of runs to do in each direction
#define STARTNX 4 // what grid size to start at
#define ENDNX  128 // what grid size to end at
//#if(NUMTRIALSX<=1)
#define STEPNX (ENDNX-STARTNX)/(NUMTRIALSX-1) // the change in grid size per run  (should really be a perfect integer)
//#else
//#define STEPNX 1 // hack
//#endif

#define NUMTRIALSY NUMTRIALSX // number of runs to do in each direction
#define STARTNY STARTNX // what grid size to start at
#define ENDNY  ENDNX // what grid size to end at
//#if(NUMTRIALSY<=1)
#define STEPNY (ENDNY-STARTNY)/(NUMTRIALSY-1) // the change in grid size per run  (should really be a perfect integer)
//#else
//#define STEPNY (1)  // hack
//#endif
/*
#define NUMTRIALSX 20 // number of runs to do in each direction
#define STARTNX 208 // what grid size to start at
#define ENDNX  512 // what grid size to end at
#define STEPNX (ENDNX-STARTNX)/(NUMTRIALSX-1) // the change in grid size per run  (should really be a perfect integer)

#define NUMTRIALSY 1 // number of runs to do in each direction
#define STARTNY 512 // what grid size to start at
#define ENDNY  512 // what grid size to end at
#define STEPNY 1 // hack
*/

int main(
         int argc,
         char *argv[],
         char *envp[]
         )
{            
  int i,j;
  int ii,jj;
  int lasti,lastj;
  int fromi,fromj;
  char part1NX[MAXFILENAME]="sed \'s/#define N1 ";
  char part3NX[MAXFILENAME]="/#define N1 ";
  char part5NX[MAXFILENAME]="/\' global.h > global.h.1";
  char part1NY[MAXFILENAME]="sed \'s/#define N2 ";
  char part3NY[MAXFILENAME]="/#define N2 ";
  char part5NY[MAXFILENAME]="/\' global.h > global.h.1";
  char csed[1000];
  
  FILE*perfout,*perfin;
  char temps[100];

  int perfgrid[NUMTRIALSX+1][NUMTRIALSY+1];


  // clean up perf dir
  sprintf(temps,"rm %sperf/*",DATADIR);
  system(temps);


  // make output file
  sprintf(temps,"mkdir %sperf",DATADIR);
  system(temps);
  sprintf(temps,"%sperf/perfgrid%s",DATADIR,DATEXT);
  if(!(perfout=fopen(temps,"wt"))){
    fprintf(stderr,"cannot open %s\n",temps);
    exit(1);
  }
  // output some info to the file as a comment
  fprintf(perfout,"#STARTNX: %d ENDNX: %d STEPNX: %d\n",STARTNX,ENDNX,STEPNX);
  fprintf(perfout,"#STARTNY: %d ENDNY: %d STEPNY: %d\n",STARTNY,ENDNY,STEPNY);
  fflush(perfout);

  // backup global.h
  system("cp global.h global.h.orig");

  // do the perf generation sequence
  lasti=STARTNX;
  lastj=STARTNY;
  ii=0;
  jj=0;
  for(j=STARTNY;j<=ENDNY;j+=STEPNY){
    for(i=STARTNX;i<=ENDNX;i+=STEPNX){

      // change N1
      if( (i==STARTNX)&&(j==STARTNY) ){ // first time ever in loop
	fromi=N1;
      }
      else{
	fromi=lasti;
      }
      sprintf(csed,"%s%d %s%d %s",part1NX,fromi,part3NX,i,part5NX);
      printf("%d\n",i);
      printf("%s\n",csed);
      system(csed);
      system("cp global.h.1 global.h");


      // change N2 (doesn't hurt to do each time--need first time)
      if((i==STARTNX)&&(j==STARTNY)){// this means first time here
	fromj=N2;
      }
      else{
	fromj=lastj;
      }
      sprintf(csed,"%s%d %s%d %s",part1NY,fromj,part3NY,j,part5NY);
      printf("%d\n",i);
      printf("%s\n",csed);
      system(csed);
      system("cp global.h.1 global.h");

      // make the binary
      system("rm *.o");
      if(PRODUCTION==0){
	system("make -j 4");
      }
      else{
	system("make -f makefile.ncsa");
      }
      //sprintf(temps,"rm %s*.txt*",DATADIR);
      //system(temps);
      sprintf(temps,"rm %s*.dat*",DATADIR);
      system(temps);
      sprintf(temps,"rm %si/*.dat*",DATADIR);
      system(temps);

      // copy the binary to the datadir
      sprintf(temps,"cp ./bin/twod %s",DATADIR);
      system(temps);

      // run the simulation
      if(PRODUCTION){
	if(NUMCPUS==1){
	  sprintf(temps,"%stwod > %s0_o.out",DATADIR,DATADIR);
	}
	else{
	  sprintf(temps,"/usr/apps/MessPass/mpich/bin/mpirun -np %d %stwod > %s0_o.out",NUMCPUS,DATADIR,DATADIR);
	}
      }
      else{
	if(NUMCPUS==1){
	  sprintf(temps,"%stwod > %s0_o.out",DATADIR,DATADIR);
	}
	else{
	  if(USEGM){
	    sprintf(temps,"mpirun.gm -np %d %stwod > %s0_o.out",NUMCPUS,DATADIR,DATADIR);
	  }
	  else{
	    if(RAINMAN){
	      sprintf(temps,"mpirun -machinefile /usr/local/share/rainman2.1 -np %d %stwod > %s0_o.out",NUMCPUS,DATADIR,DATADIR);
	    }
	    else if(PHOTON){
	      sprintf(temps,"mpirun -machinefile /usr/local/share/photon2.1 -np %d %stwod > %s0_o.out",NUMCPUS,DATADIR,DATADIR);
	    }
	    else if(WISEALPHA){
	      sprintf(temps,"mpirun -machinefile /usr/local/share/wisealpha.1 -np %d %stwod > %s0_o.out",NUMCPUS,DATADIR,DATADIR);
	    }
	    else{ // uses: /usr/local/share/machines.LINUX_ALPHA
	      sprintf(temps,"mpirun -np %d %stwod > %s0_o.out",NUMCPUS,DATADIR,DATADIR);
	    }
	  }
	}
      }
      system(temps);


      // read the perf datapoint
      sprintf(temps,"%sfinalperf%s",DATADIR,".txt");
      if((perfin=fopen(temps,"rt"))==NULL){
	printf("cannot read perf file: %s\n",temps);
	exit(1);
      }
      fscanf(perfin,"%d",&perfgrid[jj][ii]);
      fclose(perfin);
      // output while running
      fprintf(perfout,"%d %d %d\n",j,i,perfgrid[jj][ii]);
      fflush(perfout);
      fprintf(stdout,"%d %d %d\n",j,i,perfgrid[jj][ii]);
      fflush(stdout);
    
      lasti=i;
      ii++;
    } // over i
    lastj=j;
    ii=0;
    jj++;
  }// over j

  /*
  // output data so SM can read it
  for(j=STARTNY;j<=ENDNY;j+=STEPNY){
    for(i=STARTNX;i<=ENDNX;i+=STEPNX){
      fprintf(perfout,"%d %d %d\n",j,i,perfgrid[j][i]);
    }
  }
  */
  fclose(perfout);
  // done
  return(0);
}
