#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif

#include "diag.h"

/* diagnostics subroutine */

#define NUMINPUTS 15 // see how used! (overestimated)
#define NUMMODES 10 // see below how used (number of m modes in disk)

#if(SENSITIVE==1)
#define INPUT3 "%lf"
#define INPUT4 " %lf"
#define INPUT5 " %lf %lf %lf %lf %lf %lf %lf"
#define INPUT6 " %lf %lf %lf %lf %lf %lf %lf"
#define INPUT60 " %lf"
#define INPUT61 " %lf %lf"
#else
#define INPUT3 "%f"
#define INPUT4 " %f"
#define INPUT5 " %f %f %f %f %f %f %f"
#define INPUT6 " %f %f %f %f %f %f %f"
#define INPUT60 " %f"
#define INPUT61 " %f %f"
#endif

// computes averages during a given period
void diagavg(int call_code)
{
	int num1d_31,num1d_32;
	int size;
  int i,j,k,l,m;
  static int firsttime=1;
  FTYPE ftemp0,ftemp1,ftemp2,ftemp3;
  FTYPE ftemp[NUMAVG2_3+1];
  static FILE * avgout1d;
  static char temps[MAXFILENAME];
  static FTYPE d1range[2];
  static FTYPE d2range[2];
  static FTYPE dxtot1[NUMAVG1_31+1][3+1]; // true size determined by what average over [a=1,b=2][dir]
  static FTYPE dxtot2[NUMAVG1_32+1][3+1]; // true size determined by what average over [a=1,b=2][dir]
  static FTYPE dxtot1_full[NUMAVG1_31+1][3+1]; // true size determined by what average over [a=1,b=2][dir]
  static FTYPE dxtot2_full[NUMAVG1_32+1][3+1]; // true size determined by what average over [a=1,b=2][dir]
  static FTYPE (*iq)[INTN2][INTN1];
  static FTYPE *avgtemp1_31[NUMAVG1_31+1];
  static FTYPE *avgtemp1_32[NUMAVG1_32+1];
  int avgjpos;
  static FTYPE (*sigma)[3][N3M][N2M][N1M]; // -2*rho*nu*e_{ij}= sigma_{ij}
#if(USEMPI)
  static MPI_Request request[MAXCPUS]; 
#endif
  int dumpsmdump;


  dumpsmdump=DUMPSM; // always

#if(TVDLF)
  if(firsttimetvdlf==0){
    tvdlf2zeus();
  }
#endif

  if(firsttime==1){
    tavgstart=t;

    fprintf(log_file,"#proc: %s BEGIN AVGERAGE taking.\n",myidtxt);
    fflush(log_file);
    if(myid<=0){
      fprintf(logfull_file,"#BEGIN AVGERAGE taking.\n");
      fflush(logfull_file);
    }

    for(i=1;i<=NUMAVG1_31;i++){
      for(j=1;j<=3;j++){
				dxtot1[i][j]=0;
				dxtot1_full[i][j]=0;
      }
    }
    for(i=1;i<=NUMAVG1_32;i++){
      for(j=1;j<=3;j++){
				dxtot2[i][j]=0;
				dxtot2_full[i][j]=0;
      }
    }

    if(myid<=0){

      for(i=1;i<=NUMAVG1_31;i++){
				avgtemp1_31[i]=(FTYPE*)malloc(sizeof(FTYPE)*(totalsize[2]+2*N2BND));
				avgtemp1_31[i]+=N2BND;
				if(avgtemp1_31[i]==NULL){
					fprintf(fail_file,"Can't allocate memory for temp 1d averages,%d\n",i);
					myexit(1);
				}
      }
      for(i=1;i<=NUMAVG1_32;i++){
				avgtemp1_32[i]=(FTYPE*)malloc(sizeof(FTYPE)*(totalsize[1]+2*N1BND));
				avgtemp1_32[i]+=N1BND;
				if(avgtemp1_32[i]==NULL){
					fprintf(fail_file,"Can't allocate memory for temp 1d averages,%d\n",i);
					myexit(1);
				}
      }
    }
    avgcount=0;
    // initialize averages
    for(l=1;l<=NUMAVG2_3;l++){
      LOOPF{
				avg2_3[l][j][i]=0;
      }
    }
    for(l=1;l<=NUMAVG1_31;l++){
      for(j=INFULL2;j<OUTFULL2;j++){
				avg1_31[l][j]=0;
      }
    }
    for(l=1;l<=NUMAVG1_32;l++){
      for(i=INFULL1;i<OUTFULL1;i++){
				avg1_32[l][i]=0;
      }
    }
    // first guess
    // range over r to average over
    d1range[0]=1.0E-2 * (L[2][1])+L[1][1];
    d1range[1]=3.0E-2 * (L[2][1])+L[1][1];
    // range in theta to average over
    d2range[0]=(85.0/180.0)*M_PI;
    d2range[1]=(95.0/180.0)*M_PI;

    // determine number of elements to average spatially
    num1d_31=0;
    num1d_32=0;
    num1d_31_full=0;
    num1d_32_full=0;


    LOOPDIAGOUTPUT(FULLOUTPUT,FULLOUTPUT,FULLOUTPUT){
      //1-D over 3 and 1
      if( (x[1][1][i]>d1range[0])&& (x[1][1][i]<d1range[1]) ){
				if(j==0){
					num1d_31++; // get # of elements for spatial average
					dxtot1[1][1]+=dx[1][1][i];
					dxtot1[2][1]+=dx[1][1][i];
					dxtot1[3][1]+=dx[1][1][i];
					dxtot1[4][1]+=dx[1][1][i];
					dxtot1[5][1]+=dx[1][1][i];
					if(!dumpsmdump){
						dxtot1[6][1]+=dx[2][1][i];
					}
					else{
						dxtot1[6][1]+=dx[1][1][i];
					}
					dxtot1[7][1]+=dx[1][1][i];
					dxtot1[8][1]+=dx[1][1][i];
	  
					dxtot1[9][1]+=dx[1][1][i];
					if(!dumpsmdump){
						dxtot1[10][1]+=dx[2][1][i];
					}
					else{
						dxtot1[10][1]+=dx[1][1][i];
					}
					if(!dumpsmdump){
						dxtot1[11][1]+=dx[2][1][i];
					}
					else{
						dxtot1[11][1]+=dx[1][1][i];
					}
					dxtot1[12][1]+=dx[1][1][i];
					dxtot1[13][1]+=dx[1][1][i];
					dxtot1[14][1]+=dx[1][1][i];
					dxtot1[15][1]+=dx[1][1][i];
	  
				}
      }
      //1-D over 3 and 2
      if( (x[1][2][j]>d2range[0])&& (x[1][2][j]<d2range[1]) ){
				if(i==0){
					num1d_32++; // get # of elements for spatial average
					dxtot2[1][2]+=dx[1][2][j];
					dxtot2[2][2]+=dx[1][2][j];
					dxtot2[3][2]+=dx[1][2][j];
					dxtot2[4][2]+=dx[1][2][j];
					dxtot2[5][2]+=dx[1][2][j];
					dxtot2[6][2]+=dx[1][2][j];
					if(!dumpsmdump){
						dxtot2[7][2]+=dx[2][2][j];
					}
					else{
						dxtot2[7][2]+=dx[1][2][j];
					}
					dxtot2[8][2]+=dx[1][2][j];
	  
					dxtot2[9][2]+=dx[1][2][j];
					if(!dumpsmdump){
						dxtot2[10][2]+=dx[2][2][j];
					}
					else{
						dxtot2[10][2]+=dx[1][2][j];
					}
					dxtot2[11][2]+=dx[1][2][j];
					dxtot2[12][2]+=dx[1][2][j];
					if(!dumpsmdump){
						dxtot2[13][2]+=dx[2][2][j];
					}
					else{
						dxtot2[13][2]+=dx[1][2][j];
					}
					dxtot2[14][2]+=dx[1][2][j];
					dxtot2[15][2]+=dx[1][2][j];
				}
      }
    }
    if(numprocs>1){
#if(USEMPI)
      // averaged x1-dir ok for now
      num1d_31_full=num1d_31;
      // averaged x2-dir, gathered for myid=0
      MPI_Reduce(&num1d_32,&num1d_32_full, 1, MPI_INT, MPI_SUM,  0,MPI_COMM_WORLD);

      /*      
							for(i=1;i<=NUMAVG1_31;i++){
							for(j=1;j<=3;j++){
							MPI_Reduce(&dxtot1[i][j],&dxtot1_full[i][j], 1, MPI_FTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
							}
							}
      */
      for(i=1;i<=NUMAVG1_31;i++){
				for(j=1;j<=3;j++){
					dxtot1_full[i][j]=dxtot1[i][j];
				}
      }
      for(i=1;i<=NUMAVG1_32;i++){
				for(j=1;j<=3;j++){
					MPI_Reduce(&dxtot2[i][j],&dxtot2_full[i][j], 1, MPI_FTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
				}
      }

#endif
    }
    else{
      num1d_31_full=num1d_31;
      num1d_32_full=num1d_32;
      
      for(i=1;i<=NUMAVG1_31;i++){
				for(j=1;j<=3;j++){
					dxtot1_full[i][j]=dxtot1[i][j];
				}
      }
      for(i=1;i<=NUMAVG1_32;i++){
				for(j=1;j<=3;j++){
					dxtot2_full[i][j]=dxtot2[i][j];
				}
      }
    }// end if 1 cpu
  }// end if firsttime

  if(myid<=0){
    for(i=1;i<=NUMAVG1_31;i++){
      for(j=1;j<=3;j++){
				fprintf(logfull_file,"dxtot1_full[%d][%d]=%15.10g\n",i,j,dxtot1_full[i][j]);
      }
    }
    for(i=1;i<=NUMAVG1_32;i++){
      for(j=1;j<=3;j++){
				fprintf(logfull_file,"dxtot2_full[%d][%d]=%15.10g\n",i,j,dxtot2_full[i][j]);
      }
    }
    fprintf(logfull_file,"num1d_31_full=%d\n",num1d_31_full);
    fprintf(logfull_file,"num1d_32_full=%d\n",num1d_32_full);
    fflush(logfull_file);
  }
  for(i=1;i<=NUMAVG1_31;i++){
    for(j=1;j<=3;j++){
      fprintf(log_file,"dxtot1[%d][%d]=%15.10g\n",i,j,dxtot1[i][j]);
    }
  }
  for(i=1;i<=NUMAVG1_32;i++){
    for(j=1;j<=3;j++){
      fprintf(log_file,"dxtot2[%d][%d]=%15.10g\n",i,j,dxtot2[i][j]);
    }
  }
  fprintf(log_file,"num1d_31=%d\n",num1d_31);
  fprintf(log_file,"num1d_32=%d\n",num1d_32);
  fflush(log_file);

  if( (call_code==1)||( (call_code==2)&&(firsttime==1)) ){
    
    LOOPDIAGOUTPUT(FULLOUTPUT,FULLOUTPUT,FULLOUTPUT){
      // avg2_3 (full 2D time average)
      
      ftemp[1]=s[1][k][j][i]; // rho
      ftemp[2]=s[2][k][j][i]; // en
      
      if(i<OUTFULL1-1) ftemp1=e2z_1(v[1][1],k,j,i);
      else ftemp1=v[1][1][k][j][i];
      if(j<OUTFULL2-1) ftemp2=e2z_2(v[1][2],k,j,i);
      else ftemp2=v[1][2][k][j][i];

      
      ftemp3=v[1][3][k][j][i];
      
      ftemp0=ftemp1*ftemp1+ftemp2*ftemp2+ftemp3*ftemp3; // v^2 on zone center
      ftemp2=gam*s[2][k][j][i]/s[1][k][j][i]; // c_s^2/(gam-1) (zone center)
      
      ftemp[3]=ftemp0*0.5+ftemp2+s[3][k][j][i]; // Be (zone center)
      
      ftemp[4]=(gam-1.)*ftemp2; // cs^2 (zone center)
      
      ftemp[5]=(gam-1.)*s[2][k][j][i]/pow(s[1][k][j][i],gam); // Exp[S/(gam*N)] (zone center) // Exp[entropy per degree of freedom] // Really Exp[S/(gam*K*N)-(gam/(gam-1))]=(m^{gam+1}/(2\pi*k*hbar^2) * P/rho^gam)
      
      if(!dumpsmdump){
				ftemp[6]=v[1][1][k][j][i];
      }
      else ftemp[6]=e2z_1(v[1][1],k,j,i);
      if(!dumpsmdump){
				ftemp[7]=v[1][2][k][j][i];
      }
      else ftemp[7]=e2z_2(v[1][2],k,j,i);

      ftemp[8]=v[1][3][k][j][i];

#if(VISCMEM)
      ftemp[9]=sigma[1][1][k][j][i];
      if(!dumpsmdump){
				ftemp[10]=sigma[1][2][k][j][i];
      }
      else ftemp[10]=c2z_3(sigma[1][2],k,j,i);
      if(!dumpsmdump){
				ftemp[11]=sigma[1][3][k][j][i];
      }
      else ftemp[11]=c2z_2(sigma[1][3],k,j,i);
      ftemp[12]=sigma[2][2][k][j][i];
      if(!dumpsmdump){
				ftemp[13]=sigma[2][3][k][j][i];
      }
      else ftemp[13]=c2z_1(sigma[2][3],k,j,i);
      ftemp[14]=sigma[3][3][k][j][i];
      ftemp[15]=nu_real[k][j][i];
#else
      ftemp[9]=ftemp[10]=ftemp[11]=ftemp[12]=ftemp[13]=ftemp[14]=ftemp[15]=0.0;
#endif
      // 2D always time average, no spatial average
      avg2_3[1][j][i]+=ftemp[1];
      avg2_3[2][j][i]+=ftemp[2];    
      avg2_3[3][j][i]+=ftemp[3];
      avg2_3[4][j][i]+=ftemp[4];
      avg2_3[5][j][i]+=ftemp[5];
      avg2_3[6][j][i]+=ftemp[6];
      avg2_3[7][j][i]+=ftemp[7];
      avg2_3[8][j][i]+=ftemp[8];
      avg2_3[9][j][i]+=ftemp[9];
      avg2_3[10][j][i]+=ftemp[10];
      avg2_3[11][j][i]+=ftemp[11];
      avg2_3[12][j][i]+=ftemp[12];
      avg2_3[13][j][i]+=ftemp[13];
      avg2_3[14][j][i]+=ftemp[14];
      avg2_3[15][j][i]+=ftemp[15];

    }
    avgcount++;
  }
  if( (call_code==2)&&(avgcount>=0)){
    tavgfinal=t;

    fprintf(log_file,"#proc: %d BEGIN AVGERAGE DUMPING.\n",myid);
    fflush(log_file);
    if(myid<=0){
      fprintf(logfull_file,"#BEGIN AVGERAGE DUMPING.\n");
      fflush(logfull_file);
    }
    //if(tavgi>=t) tavgstart=t; else tavgstart=tavgi;
    //if(tavgf>=t) tavgfinal=t; else tavgfinal=tavgf;

    // divide out number of calls
    // assumes times between each call is equal(DTavg)
   

    LOOPDIAGOUTPUT(FULLOUTPUT,FULLOUTPUT,FULLOUTPUT){
      // 2D always time average, no spatial average
      avg2_3[1][j][i]/=((FTYPE)(avgcount));
      avg2_3[2][j][i]/=((FTYPE)(avgcount));
      avg2_3[3][j][i]/=((FTYPE)(avgcount));
      avg2_3[4][j][i]/=((FTYPE)(avgcount));
      avg2_3[5][j][i]/=((FTYPE)(avgcount));
      avg2_3[6][j][i]/=((FTYPE)(avgcount));
      avg2_3[7][j][i]/=((FTYPE)(avgcount));
      avg2_3[8][j][i]/=((FTYPE)(avgcount));
      avg2_3[9][j][i]/=((FTYPE)(avgcount));
      avg2_3[10][j][i]/=((FTYPE)(avgcount));
      avg2_3[11][j][i]/=((FTYPE)(avgcount));
      avg2_3[12][j][i]/=((FTYPE)(avgcount));
      avg2_3[13][j][i]/=((FTYPE)(avgcount));
      avg2_3[14][j][i]/=((FTYPE)(avgcount));
      avg2_3[15][j][i]/=((FTYPE)(avgcount));
    }
    LOOPDIAGOUTPUT(FULLOUTPUT,FULLOUTPUT,FULLOUTPUT){
      //1-D over 3 and 1
      if( (x[1][1][i]>d1range[0])&& (x[1][1][i]<d1range[1]) ){
				avg1_31[1][j]+=avg2_3[1][j][i]*dx[1][1][i];
				avg1_31[2][j]+=avg2_3[2][j][i]*dx[1][1][i];
				avg1_31[3][j]+=avg2_3[3][j][i]*dx[1][1][i];
				avg1_31[4][j]+=avg2_3[4][j][i]*dx[1][1][i];
				avg1_31[5][j]+=avg2_3[5][j][i]*dx[1][1][i];
				if(!dumpsmdump){
					avg1_31[6][j]+=avg2_3[6][j][i]*dx[2][1][i];
				}
				else{
					avg1_31[6][j]+=avg2_3[6][j][i]*dx[1][1][i];
				}
				avg1_31[7][j]+=avg2_3[7][j][i]*dx[1][1][i];
				avg1_31[8][j]+=avg2_3[8][j][i]*dx[1][1][i];
				avg1_31[9][j]+=avg2_3[9][j][i]*dx[1][1][i];
				if(!dumpsmdump){
					avg1_31[10][j]+=avg2_3[10][j][i]*dx[2][1][i];
				}
				else{
					avg1_31[10][j]+=avg2_3[10][j][i]*dx[1][1][i];
				}
				if(!dumpsmdump){
					avg1_31[11][j]+=avg2_3[11][j][i]*dx[2][1][i];
				}
				else{
					avg1_31[11][j]+=avg2_3[11][j][i]*dx[1][1][i];
				}
				avg1_31[12][j]+=avg2_3[12][j][i]*dx[1][1][i];
				avg1_31[13][j]+=avg2_3[13][j][i]*dx[1][1][i];
				avg1_31[14][j]+=avg2_3[14][j][i]*dx[1][1][i];
				avg1_31[15][j]+=avg2_3[15][j][i]*dx[1][1][i];
      }
    }
    LOOPDIAGOUTPUT(FULLOUTPUT,FULLOUTPUT,FULLOUTPUT){
      //1-D over 3 and 2
      if( (x[1][2][j]>d2range[0])&& (x[1][2][j]<d2range[1]) ){
				avg1_32[1][i]+=avg2_3[1][j][i]*dx[1][2][j];
				avg1_32[2][i]+=avg2_3[2][j][i]*dx[1][2][j];
				avg1_32[3][i]+=avg2_3[3][j][i]*dx[1][2][j];
				avg1_32[4][i]+=avg2_3[4][j][i]*dx[1][2][j];
				avg1_32[5][i]+=avg2_3[5][j][i]*dx[1][2][j];
				avg1_32[6][i]+=avg2_3[6][j][i]*dx[1][2][j];
				if(!dumpsmdump){
					avg1_32[7][i]+=avg2_3[7][j][i]*dx[2][2][j];
				}
				else{
					avg1_32[7][i]+=avg2_3[7][j][i]*dx[1][2][j];
				}
				avg1_32[8][i]+=avg2_3[8][j][i]*dx[1][2][j];

				avg1_32[9][i]+=avg2_3[9][j][i]*dx[1][2][j];
				if(!dumpsmdump){
					avg1_32[10][i]+=avg2_3[10][j][i]*dx[2][2][j];
				}
				else{
					avg1_32[10][i]+=avg2_3[10][j][i]*dx[1][2][j];
				}
				avg1_32[11][i]+=avg2_3[11][j][i]*dx[1][2][j];
				avg1_32[12][i]+=avg2_3[12][j][i]*dx[1][2][j];
				if(!dumpsmdump){
					avg1_32[13][i]+=avg2_3[13][j][i]*dx[2][2][j];
				}
				else{
					avg1_32[13][i]+=avg2_3[13][j][i]*dx[1][2][j];
				}
				avg1_32[14][i]+=avg2_3[14][j][i]*dx[1][2][j];
				avg1_32[15][i]+=avg2_3[15][j][i]*dx[1][2][j];
      }
    }
    // now sum up all cpus total, giving total to myid=0
    if(numprocs>1){
#if(USEMPI)
      for(l=1;l<=NUMAVG1_31;l++){
	
				// for myid==0
				if(myid<=0){
					avgjpos=INFULL2;
					for(j=avgjpos;j<sizes[2][0];j++){
						avgtemp1_31[l][j]=avg1_31[l][j];
					}
					avgjpos=sizes[2][0]; // next position to write to
				}
				for(i=1;i<numprocs;i++){
					if(i==numprocs-1) size=sizes[2][i]+N2BND;
					else size=sizes[2][i];

					if(i==myid){
						MPI_Isend(avg1_31[l],size,MPI_FTYPE,0,i,MPI_COMM_WORLD,&request[i]);
						MPI_Wait(&request[i],&mpichstatus);
					}
					if(myid==0){
						MPI_Irecv(&avgtemp1_31[l][avgjpos],size,MPI_FTYPE,i,i,MPI_COMM_WORLD,&request[i]);
						avgjpos+=sizes[2][i];
					}
				}
				if(myid<=0){
					for(i=1;i<numprocs;i++){
						MPI_Wait(&request[i],&mpichstatus);
					}
				}
				if((myid<=0)&&(num1d_31_full>0)){
					// now normalize correctly
					for(j=INFULL2;j<totalsize[2]+N2BND;j++){
						avgtemp1_31[l][j]=avgtemp1_31[l][j]/dxtot1_full[l][1]; // assign total back to myid=0 and get real average result
					}
				}
				/*
					if(myid<=0){
					for(j=0;j<totalsize[2];j++){
					avgtemp1_31[l][j]=0;
					}
					}
				*/
	
      }

      for(l=1;l<=NUMAVG1_32;l++){
	
				for(i=0;i<totalsize[1];i++){
					MPI_Reduce(&(avg1_32[l][i]),   &(ftemp0), 1, MPI_FTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
					if((myid<=0)&&(num1d_32_full>0)){
						// assign and normalize
						avgtemp1_32[l][i]=ftemp0/dxtot2_full[l][2]; // assign total back to myid=0 and get real average result
					}
				}
				/*	
					if(myid<=0){
					for(i=0;i<totalsize[1];i++){
					avgtemp1_32[l][i]=0;
					}
					}
				*/
      }
#endif
    }
    else{
      if((myid<=0)&&(num1d_31_full>0)){
				for(l=1;l<=NUMAVG1_31;l++){
					for(j=INFULL2;j<totalsize[2]+N2BND;j++){
						avgtemp1_31[l][j]=avg1_31[l][j]/dxtot1_full[l][1]; // assign total back to myid=0 and get real average result
					}
				}
      }
      if((myid<=0)&&(num1d_32_full>0)){
				for(l=1;l<=NUMAVG1_32;l++){
					for(i=INFULL1;i<totalsize[1]+N1BND;i++){
						avgtemp1_32[l][i]=avg1_32[l][i]/dxtot2_full[l][2]; // assign total back to myid=0 and get real average result
					}
				}
      }
    }

    // all cpus dump their 2d data
    dump(NULL,-1,AVG2DTYPE,0);// primitive avg2d data

    if(myid<=0){
    
      // open file
      sprintf(temps,"%s0_avg1d%s",DATADIR,DAT2EXT) ;    
      if((avgout1d=fopen(temps,"wt"))==NULL){
				fprintf(fail_file,"avg1d: Cannot open: %s\n",temps);
				myexit(1);
      }
      dump_header(avgout1d,AVG1DTYPE,0,0);
      dump_header2(avgout1d,AVG1DTYPE);
    
      // shouldn't need to interp 1-d since no contour plotting for it      
      LOOPDIAGOUTPUTFULL(FULLOUTPUT,FULLOUTPUT,FULLOUTPUT){
				fprintf(avgout1d," %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",avgtemp1_31[1][j],avgtemp1_31[2][j],avgtemp1_31[3][j],avgtemp1_31[4][j],avgtemp1_31[5][j],avgtemp1_31[6][j],avgtemp1_31[7][j],avgtemp1_31[8][j],avgtemp1_31[9][j],avgtemp1_31[10][j],avgtemp1_31[11][j],avgtemp1_31[12][j],avgtemp1_31[13][j],avgtemp1_31[14][j],avgtemp1_31[15][j],avgtemp1_32[1][i],avgtemp1_32[2][i],avgtemp1_32[3][i],avgtemp1_32[4][i],avgtemp1_32[5][i],avgtemp1_32[6][i],avgtemp1_32[7][i],avgtemp1_32[8][i],avgtemp1_32[9][i],avgtemp1_32[10][i],avgtemp1_32[11][i],avgtemp1_32[12][i],avgtemp1_32[13][i],avgtemp1_32[14][i],avgtemp1_32[15][i]);
      }

      fflush(avgout1d);
      fclose(avgout1d);

      for(i=1;i<=NUMAVG1_31;i++){
				avgtemp1_31[i]-=N2BND;
				free(avgtemp1_31[i]);
      }
      for(i=1;i<=NUMAVG1_32;i++){
				avgtemp1_32[i]-=N1BND;
				free(avgtemp1_32[i]);
      }
    }// end if final call and writting cpu
    fprintf(log_file,"#proc: %d  AVERAGE DUMPING DONE!\n",myid);
    fflush(log_file);
    if(myid<=0){
      fprintf(logfull_file,"#AVERAGE DUMPING DONE!\n");
      fflush(logfull_file);
    }
  }// end if final call
  firsttime=0;
}



void diag(int call_code)
	// call_code
	// -1: do every call to diag(usually each time step as called by main.c)
	// 0: first call setup/dump
	// 1: normal diag dump of variables each DTl (DT for logging)
	// 2: final diag call after done
     
{
  static int firsttime=1;
  static SFTYPE totloss_full[NUMLOSSVAR+1]; // 0: etot
  static FTYPE tener,tloss; // used for multiple functions
  static int enerc,lossc;
  static FTYPE tdump,tfldump,tpdump,timage,tsym,tdivb ;
  static int dumpc,pdumpc,imagec,fldumpc,symc,divbc; // time-based counts for this run
  // for diag_mode()
  static FTYPE tmode ;
  static int modec;
  static FTYPE DTsym;
  static int firstenercall=0; // used to call divb when first calling ener


  //  if(t>32) DTd=1E-5;

  
  if(firsttime){
    // GODMARK -- should give DTsym its own DT
    DTsym=DTmode*5;

    enerc=(int)((t-tstart)/DTener)+ireenter;
    tener = tstart+(FTYPE)(enerc)*DTener ; // next ener time

    lossc=(int)((t-tstart)/DTloss)+ireenter;
    tloss = tstart+(FTYPE)(lossc)*DTloss ; // next loss time
    
    modec=(int)((t-tstart)/DTmode)+ireenter;
    tmode = tstart+(FTYPE)(modec)*DTmode ; // next ener time

    symc=(int)((t-tstart)/DTsym)+ireenter; // GODMARK -- uses DTi for now
    tsym = tstart+(FTYPE)(symc)*DTsym ; // next ener time

    divbc=(int)((t-tstart)/DTdivb)+ireenter;
    tdivb = tstart+(FTYPE)(divbc)*DTdivb ; // next divb time


    if(myid<=0){
      fprintf(logfull_file,"t: %15.10g tstart: %15.10g timereenter: %15.10g enerc: %10d %15.10g lossc: %10d %15.10g modec: %10d %15.10g symc: %10d %15.10g\n",t,tstart,timereenter,enerc,tener,lossc,tloss,modec,tmode,symc,tsym);
      fflush(logfull_file);
    }

    // this is count of when to check next to dump
    // same as in init_reentrance2 except when specifying dump#'s manually
    pdumpc=(int)((t-tstart)/DTpd)+ireenter;
    tpdump = tstart+((FTYPE)(pdumpc)*DTpd);

    dumpc=(int)((t-tstart)/DTd)+ireenter;
    tdump = tstart+((FTYPE)(dumpc)*DTd);

    fldumpc=(int)((t-tstart)/DTfld)+ireenter;
    tfldump = tstart+((FTYPE)(fldumpc)*DTfld);

    imagec=(int)((t-tstart)/DTi)+ireenter;
    timage = tstart+((FTYPE)(imagec)*DTi);


    if(myid<=0){
      fprintf(logfull_file,"t: %15.10g tstart: %15.10g timereenter: %15.10g pdumpc: %10d %15.10g dumpc: %10d %15.10g imagec: %10d %15.10g\n",t,tstart,timereenter,pdumpc,tpdump,dumpc,tdump,imagec,timage);
      fflush(logfull_file);
    }

  }

  if(DOLOSSDIAG){
    // determine flux totals, write flux totals
    if((LOOPTYPE==1)||((t>=tener) ||(t>=tloss)||(call_code==2))){ // do each time if LOOPTYPE==1

      diag_loss(call_code,tener,tloss,totloss_full);
      
      if((t>=tloss)||(call_code==2)){ // only increment time if t>tloss (different than rest due to summing each timestep up inside diag_loss for rect version
				lossc++;
				tloss = tstart+(FTYPE)(lossc)*DTloss ;
      }
    }
  }


  if(call_code>=0){
    
#if(TVDLF==1)
    if((
				 ( (call_code>=0)&& ((t>=tener)||(call_code==2) ))||
				 ( ((dt<DTLOWDUMP)&&(CHECKDTLOW==1))||(t >= tpdump) || (call_code == 2) )||
				 ( ((dt<DTLOWDUMP)&&(CHECKDTLOW==1))||(t >= tdump) || (call_code == 2) ) ||
				 ( (t >= timage)|| (call_code == 2) )
				 )){
      if(firsttimetvdlf==0){
				tvdlf2zeus();
      }
    }
#endif

    // determine integrals, write integrals

    if((firstenercall==0)&&((t>=tener)||(call_code==2))){ // needed since restart needs to call divb before ener
      firstenercall=1;
    }

    if(TVDLF==0){
      if(DODIVBDIAG==1){
				//	if((t>=tdivb)||(t>=tener)||(call_code==2)){
				if((t>=tdivb)||(firstenercall==1)||(call_code==2)){
	  
					divb0check(1); // should come before diag_ener when ener outputs divbmax

					divbc++;
					tdivb = tstart+(FTYPE)(divbc)*DTdivb ; // must before after diag_ener()
				}
      }    
    }
    // determine integrals, write integrals
    if((t>=tener)||(call_code==2)){
      
      diag_ener(call_code,totloss_full);
      firstenercall=2;
      
      enerc++;
      tener = tstart+(FTYPE)(enerc)*DTener ; // must come after diag_loss()
    }

    // DETERMINE MODES, write modes
    if((t>=tmode)||(call_code==2)){
      
      diag_mode(call_code);
      
      modec++;
      tmode = tstart+(FTYPE)(modec)*DTmode ;
    }


    // DETERMINE symmetry
    if((t>=tsym)||(call_code==2)){
      
      symmetry_check(0); // 0=assumes normal scalar position
      
      symc++;
      tsym = tstart+(FTYPE)(symc)*DTsym ; // GODMARK -- uses DTi for now
    }

    
    // DATA DUMPS  
    diag_dumps(call_code,tdump,tfldump,tpdump,timage); // GODMARK -- should seperate dump calls
    

    if(PDUMPFLAG&&( ((dt<DTLOWDUMP)&&(CHECKDTLOW==1))||(t >= tpdump) || (call_code == 2) ) ) {
      pdumpc++;
      tpdump = tstart+(FTYPE)(pdumpc)*DTpd ;
    }
    if(((DUMPFLAG)||(NPDUMPFLAG)||(FLOORDUMPFLAG)||(ADUMPFLAG&&analoutput>0))&&(((dt<DTLOWDUMP)&&(CHECKDTLOW==1))||(t >= tdump) || (call_code == 2))  ) {
      dumpc++;
      tdump = tstart+(FTYPE)(dumpc)*DTd ;
    }
    if(FLDUMPFLAG&&((t >= tfldump) || (call_code == 2) ) ) {
      fldumpc++;
      tfldump = tstart+(FTYPE)(fldumpc)*DTfld ;
    }
    if(IMAGEFLAG&&( (t >= timage)|| (call_code == 2) ) ) {
      imagec++;
      timage = tstart+(FTYPE)(imagec)*DTi ;
    }
  }


  
  firsttime=0;
}



void diag_loss_rect(int call_code,FTYPE tener,FTYPE tloss, SFTYPE* totloss_full)
	// call_code
	// 0: first call setup/dump
	// 1: normal diag dump of variables each DTl (DT for logging)
	// 2: final diag call after done
     
{
  static SFTYPE told=-1000.00;
  fpos_t fposgod;
  SFTYPE tempf,tempf2;
  SFTYPE tempfns;
  SFTYPE masstempc,masstemp1,masstemp2,masstemp3;
  char dfnam[MAXFILENAME];
  char dfheader[MAXFILENAME];
  char dfnamtemp[MAXFILENAME];
  char dfnamback[MAXFILENAME];


  SFTYPE vxa,vya,vza ;
  SFTYPE Bxa,Bya,Bza ;
  SFTYPE dxdyc,dxdy1,dxdy2,dxdy3;
  int i,j,k,l,m,ii,ll ;
  int loopdido;
  int dumi[10];
  SFTYPE tcheck; // used to check time for loss append

  int which, comp;
  SFTYPE realloss; // need for enthalpy and sm
  SFTYPE ftemp,ftemp1,ftemp2,ftemp3;
  SFTYPE floattemp[NUMLOSSVAR+1][NUMINPUTS+1]; // NUMINPUTS types of output for 0_ener.dat per variable
  
  static int dump_cnt,pdump_cnt,im_cnt ; // used to enumerate files
  SFTYPE tempsumloss[NUMLOSSVAR+1][6];
  static SFTYPE sumloss[NUMLOSSVAR+1][6];
  static SFTYPE sumloss_full[NUMLOSSVAR+1][6];
  static SFTYPE totloss[NUMLOSSVAR+1]; // 0: etot
  static FILE *loss_file;
  static FILE *loss_file_temp;
  static FILE *loss_file2;
  static FILE *final_output;
  static int timeless=0;
  static int firsttime=1;
#if(USEMPI)
  static MPI_Request request[6]; 
#endif
  int typenum;
  char temps[MAXFILENAME];
  char filename[MAXFILENAME];
  char filenametemp[MAXFILENAME];
  char filenameback[MAXFILENAME];
  long fpos0;
  int gotit;
  // for easy coord conv

 
  // tvdlf has no flux counters yet

  // if every dt or initial or final, do special diags not otherwise done
  if(firsttime==1){ // the check for if first time here
    for(i=1;i<=NUMLOSSVAR;i++){
      totloss[i]=0.0;
    }
    if(myid<=0){
      
      if(DETAILMLOSS>=0){
	
				sprintf(temps,DATADIR);
	
				sprintf(dfnam,"%s0_loss%s",temps,DAT2EXT);
				if((loss_file = fopen(dfnam,WRITETYPE))==NULL) {
					fprintf(fail_file,"error opening loss output file %s\n",dfnam) ;
					myexit(1) ;
				}
				if((restartloss==0)||(appendold==0)){
					// version header
					fprintf(loss_file,"#%10s\n%10d %10d\n","LOSSVER",LOSSVER,LOSSTYPE);
					if(DETAILMLOSS==0){
						fprintf(loss_file,"#%21s %21s %21s %21s %21s\n","time","Total Mass Loss","Total Energy Loss","Total AngMom Loss","Total ViscE Loss");
					}
					else if(DETAILMLOSS>=1){
						if(REALNUMVEC==NUMVEC){
							fprintf(loss_file,"#%21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s\n"
											,"time"
											,"Total Mass Loss"   ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total IEnergy Loss","Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total PEnergy Loss","Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total KEnergy Loss","Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x1-Mom Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x2-Mom Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x3-Mom Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total MEnergy Loss","Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x1-Bflux Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x2-Bflux Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x3-Bflux Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total visc-e Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x1-AMom Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x2-AMom Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x3-AMom Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
								);
						}
						else{
							fprintf(loss_file,"#%21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s"
											" %21s %21s %21s %21s %21s %21s %21s\n"
											,"time"
											,"Total Mass Loss"   ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total IEnergy Loss","Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total PEnergy Loss","Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total KEnergy Loss","Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x1-Mom Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x2-Mom Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x3-Mom Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total visc-e Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x1-AMom Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x2-AMom Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
											,"Total x3-AMom Loss" ,"Inner x1","Outer x1","Inner x2","Outer x2","Inner x3","Outer x3"
								);
						}
					}
				}//if appendold==0
				else{//if appendold==1
					if(DETAILMLOSS<1){
						if(myid<=0){
							fprintf(logfull_file,"appendold==1 with detailmassloss<1 won't allow correct append reentrance to loss data file\n");
							fflush(logfull_file);
						}
					}
					else{//if detailmassloss>=1
						if(myid<=0){
							fprintf(logfull_file,"Start setup of loss file append\n");
							fflush(logfull_file);
						}
						// need to read in totloss and sumloss and stick on right cpu
						// need to make sure getting right time.  If DTloss<DTd||DTi, then should be able to get easily.  You should generally have this true anyways
	    
						rewind(loss_file); // go to start
						// check version info
						while(fgetc(loss_file)!='\n'); // skip comment line
						fscanf(loss_file,"%d %d",&dumi[0],&dumi[1]);
						if( (dumi[0]!=LOSSVER)||(dumi[1]!=LOSSTYPE) ){
							fprintf(fail_file,"Expected lossver/losstype: %d %d got %d %d\n",LOSSVER,LOSSTYPE,dumi[0],dumi[1]);
							myexit(6);
						}
						while(fgetc(loss_file)!='\n'); // skip to next line
						while(fgetc(loss_file)!='\n'); // skip comment line	      
						gotit=0;
						while( (!feof(loss_file))&&(gotit==0) ){
	      
							fscanf(loss_file,INPUT3,&tcheck);
							//fprintf(stderr,"%15.10g %15.10g %d\n",t,tcheck,gotit);
							if(fabs(tcheck-t)<0.5*DTloss){
								gotit=1;
								for(l=1;l<=NUMLOSSVAR;l++){
									if(REALNUMVEC!=NUMVEC){
										if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
											l=NUMSCA+NUMVEC*(3+1)+1;
										}
									}
									fscanf(loss_file,INPUT6,&totloss_full[l],&sumloss_full[l][0],&sumloss_full[l][1],&sumloss_full[l][2],&sumloss_full[l][3],&sumloss_full[l][4],&sumloss_full[l][5]);
								}
							}
							else{
								while((fgetc(loss_file)!='\n')&&(!feof(loss_file))); // skip this bad line
							}
							// continue after successful get since successful get is good data and should keep since corresponds to dump one is keeping
							fpos0=ftell(loss_file); // position to continue writting at if successful get
						}
						if(gotit==0){
							fprintf(fail_file,"Never found right time in loss file when appending: looking for t=%21.15g lastt=%21.15g\n",timereenter,tcheck);
							myexit(1);
						}
						else{
							fprintf(logfull_file,"found goodtime t=%21.15g (wanted %21.15g) to restart loss file\n",tcheck,timereenter);
							sprintf(temps,DATADIR);
							sprintf(dfnamtemp,"%s0_loss%s.temp",temps,DAT2EXT);
							sprintf(dfnam,"%s0_loss%s",temps,DAT2EXT);
							sprintf(dfnamback,"%s0_loss%s.back",temps,DAT2EXT);
	      
							// now that done, fix up file
							if( (loss_file_temp=fopen(dfnamtemp,"wt"))==NULL){
								fprintf(fail_file,"Cannot open temp loss file for appending: %s\n",dfnamtemp);
								myexit(1);
							}
							else{
								rewind(loss_file);
								while(ftell(loss_file)<fpos0+1){// +1 is for '\n' at end of line
									fputc(fgetc(loss_file),loss_file_temp);
								}
								fclose(loss_file_temp);
								fclose(loss_file);
								rename(dfnam,dfnamback); // move old to backup location
								rename(dfnamtemp,dfnam); // move new to old name(normal name)
								// reopen loss_file (now normal name)
								if((loss_file = fopen(dfnam,"at"))==NULL) {
									fprintf(fail_file,"2: error opening loss output file %s\n",dfnam) ;
									myexit(1) ;
								}
								if(myid<=0){
									fprintf(logfull_file,"End setup of loss file append\n");
									fflush(logfull_file);
								}
							}
						}
					}//end else if detailmassloss>=1
				}//end else if appendold==1
      }// end if doing normal mloss output
      
      if(DETAILMLOSS==2){
	
	
				sprintf(temps,DATADIR);
	
				sprintf(dfnam,"%s0_lossd%s",temps,DAT2EXT);
				if((loss_file2 = fopen(dfnam,WRITETYPE))==NULL) {
					fprintf(fail_file,"error opening loss detail output file %s\n",dfnam) ;
					myexit(1) ;
				}
	
				if((restartloss==0)||(appendold==0)){
					if(REALNUMVEC==NUMVEC){
						fprintf(loss_file2,"#%21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s\n"
										,"time"
										,"i-coord","j-coord","k-coord","Mass Loss"
										,"i-coord","j-coord","k-coord","IEnergy Loss"
										,"i-coord","j-coord","k-coord","PEnergy Loss"
										,"i-coord","j-coord","k-coord","KEnergy Loss"
										,"i-coord","j-coord","k-coord","s1 loss"
										,"i-coord","j-coord","k-coord","s2 loss"
										,"i-coord","j-coord","k-coord","s3 loss"
										,"i-coord","j-coord","k-coord","BEnergy Loss"
										,"i-coord","j-coord","k-coord","Bf1 loss"
										,"i-coord","j-coord","k-coord","Bf2 loss"
										,"i-coord","j-coord","k-coord","Bf3 loss"
										,"i-coord","j-coord","k-coord","visc-e loss"
										,"i-coord","j-coord","k-coord","a1 loss"
										,"i-coord","j-coord","k-coord","a2 loss"
										,"i-coord","j-coord","k-coord","a3 loss"
							);
					}
					else{
						fprintf(loss_file2,"#%21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s"
										" %7s %7s %7s %21s\n"
										,"time"
										,"i-coord","j-coord","k-coord","Mass Loss"
										,"i-coord","j-coord","k-coord","IEnergy Loss"
										,"i-coord","j-coord","k-coord","PEnergy Loss"
										,"i-coord","j-coord","k-coord","KEnergy Loss"
										,"i-coord","j-coord","k-coord","s1 loss"
										,"i-coord","j-coord","k-coord","s2 loss"
										,"i-coord","j-coord","k-coord","s3 loss"
										,"i-coord","j-coord","k-coord","visc-e loss"
										,"i-coord","j-coord","k-coord","a1 loss"
										,"i-coord","j-coord","k-coord","a2 loss"
										,"i-coord","j-coord","k-coord","a3 loss"
							);
					}
				}
      }
    }// end if write cpu
    if((appendold==1)&&(restartloss==1)){
      // now must distribute loss data to appropriate cpu from root==0
      // given only x2-dir splitting
      if(myid<=0){
				fprintf(logfull_file,"Begin transfer setup of loss file append\n");
				fflush(logfull_file);
      }
      for(l=1;l<=NUMLOSSVAR;l++){
				if(REALNUMVEC!=NUMVEC){
					if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
						l=NUMSCA+NUMVEC*(3+1)+1;
					}
	  
				}
				if(myid<=0){
					// cpu=0 can keep total losses
					totloss[l]=totloss_full[l];
					// cpu=0 can keep inner radial fluxes(no way to recapture that seperate data anyways), and cpu=0 needs to keep theta=0 surface flux
					sumloss[l][0]=sumloss_full[l][0];
					//	    sumloss[l][1]=sumloss_full[l][1];
					sumloss[l][2]=sumloss_full[l][2];
	  
					sumloss[l][4]=sumloss_full[l][4];
				}
				if(numprocs>1){
#if(USEMPI)
					// send outer radial flux to an outer cpu
					if(myid==0){
						MPI_Isend(&sumloss_full[l][1],1,MPI_SFTYPE,numprocs-1,99,MPI_COMM_WORLD,&request[0]);
						MPI_Wait(&request[0],&mpichstatus);
					}
					else if(myid==numprocs-1){
						MPI_Irecv(&sumloss[l][1],1,MPI_SFTYPE,0,99,MPI_COMM_WORLD,&request[1]);
						MPI_Wait(&request[1],&mpichstatus);
					}
#endif
				}
				else{
					sumloss[l][1]=sumloss_full[l][1];
				}
				if(numprocs>1){
#if(USEMPI)
					// cpu=numprocs-1 needs theta=M_PI surface flux
					if(myid==0){
						MPI_Isend(&sumloss_full[l][3],1,MPI_SFTYPE,numprocs-1,99,MPI_COMM_WORLD,&request[0]);
						MPI_Wait(&request[0],&mpichstatus);
					}
					else if(myid==numprocs-1){
						MPI_Irecv(&sumloss[l][3],1,MPI_SFTYPE,0,99,MPI_COMM_WORLD,&request[1]);
						MPI_Wait(&request[1],&mpichstatus);
					}
#endif
				}
				else{
					sumloss[l][3]=sumloss_full[l][3];
				}
				if(numprocs>1){
#if(USEMPI)
					// cpu=numprocs-1 needs phi=2*M_PI surface flux
					if(myid==0){
						MPI_Isend(&sumloss_full[l][5],1,MPI_SFTYPE,numprocs-1,99,MPI_COMM_WORLD,&request[0]);
						MPI_Wait(&request[0],&mpichstatus);
					}
					else if(myid==numprocs-1){
						MPI_Irecv(&sumloss[l][5],1,MPI_SFTYPE,0,99,MPI_COMM_WORLD,&request[1]);
						MPI_Wait(&request[1],&mpichstatus);
					}
#endif
				}
				else{
					sumloss[l][5]=sumloss_full[l][5];
				}
      }// over loss loop
      if(myid<=0){
				fprintf(logfull_file,"End transfer setup of loss file append\n");
				fflush(logfull_file);
      }
    }// if appendold==1
  }// end if first time in here      
  

  ///////////////////////////////////////////////////////////////
  //
  // DONE EVERY TIMESTEP
  // compute total mass lost in this timestep
  //
  //
  for(l=1;l<=NUMLOSSVAR;l++){
    
    if(l<=NUMSCA){ m=-1; typenum=0;}
    else if(l==NUMSCA+1){ m=0; typenum=1; }
    else if(l==NUMSCA+1+1){ m=1; typenum=1; }
    else if(l==NUMSCA+1+2){ m=2; typenum=1; }
    else if(l==NUMSCA+1+3){ m=3; typenum=1; }
    else if(l==NUMSCA+1*(3+1)+1){ m=0; typenum=2; }
    else if(l==NUMSCA+1*(3+1)+1+1){ m=1; typenum=2; }
    else if(l==NUMSCA+1*(3+1)+2+1){ m=2; typenum=2; }
    else if(l==NUMSCA+1*(3+1)+3+1){ m=3; typenum=2; }
    else if(l==NUMSCA+2*(3+1)+1){ m=-2; }
    if(REALNUMVEC!=NUMVEC){
      if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
				l=NUMSCA+NUMVEC*(3+1)+1;
				m=-2;  
      }
    }
    
    j=1; // x1-direction
    for(k=0;k<2;k++){//i/o
      // add up stuff that just went through boundary
      tempsumloss[l][k]=0;
      if(N1>1){
				// only add those cpus that are on inner or outer boundary for right loop
				if( ((k==0)&&(mycpupos[1]==0))||((k==1)&&(mycpupos[1]==ncpux1-1)) ) {
					for(ii=0;ii<N3;ii++) for(i=0;i<N2;i++){//list
						if(m==-1) tempf=losss[l][j][k][ii][i];
						else if(m==-2) tempf=lossvisc[1][j][k][ii][i];
						else tempf=lossv[typenum][m][j][k][ii][i];
						tempsumloss[l][k]+=tempf;
					}
					// add new quantity passing through boundary
					sumloss[l][k]+=tempsumloss[l][k];
				}
      }
      else sumloss[l][k]=0;
    }
    j=2; // x2-direction
    for(k=0;k<2;k++){//i/o
      tempsumloss[l][k+2]=0;
      if(N2>1){
				// only outer cpus are on boundary
				if( ((k==0)&&(mycpupos[2]==0))||((k==1)&&(mycpupos[2]==ncpux2-1)) ) {
					for(ii=0;ii<N3;ii++) for(i=0;i<N1;i++){//list
						if(m==-1) tempf=losss[l][j][k][ii][i];
						else if(m==-2) tempf=lossvisc[1][j][k][ii][i];
						else tempf=lossv[typenum][m][j][k][ii][i];
						tempsumloss[l][k+2]+=tempf;
					}
					sumloss[l][k+2]+=tempsumloss[l][k+2];
				}
      }
      else sumloss[l][k+2]=0;
    }
    j=3; // x3-direction
    for(k=0;k<2;k++){//i/o
      tempsumloss[l][k+4]=0;
      if(N3>1){
				// only outer cpus are on boundary
				if( ((k==0)&&(mycpupos[3]==0))||((k==1)&&(mycpupos[3]==ncpux3-1)) ) {
					for(ii=0;ii<N2;ii++) for(i=0;i<N1;i++){//list
						if(m==-1) tempf=losss[l][j][k][ii][i];
						else if(m==-2) tempf=lossvisc[1][j][k][ii][i];
						else tempf=lossv[typenum][m][j][k][ii][i];
						tempsumloss[l][k+4]+=tempf;
					}
					sumloss[l][k+4]+=tempsumloss[l][k+4];
				}
      }
      else sumloss[l][k+4]=0;
    }
    // only add new stuff to total sum ( 6 surfaces in 3D)
    for(k=0;k<6;k++){
      totloss[l]+=tempsumloss[l][k];
    }
  }


  ////////////////////////
  //  
  // must compute sum over cpus everytime since needed for diagnostics on DTl period, or for loss diagnostics
  //
  //
  if((t>=tener) ||(t>=tloss)||(call_code==2)){
    
    for(l=1;l<=NUMLOSSVAR;l++){
      
      if(REALNUMVEC!=NUMVEC){
				if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
					l=NUMSCA+NUMVEC*(3+1)+1;
					m=-2;  
				}
      }
      
      if(numprocs>1){
#if(USEMPI)
				MPI_Reduce(&totloss[l], &totloss_full[l], 1, MPI_SFTYPE, MPI_SUM, 0, MPI_COMM_WORLD);
	
				// this sum only works, because of check above when setting these values(irrelevant cpus are 0, and added anyways currently.)
				// SUPERMARK -- make add only needed cpus
				for(k=0;k<6;k++){
					MPI_Reduce(&sumloss[l][k], &sumloss_full[l][k], 1, MPI_SFTYPE, MPI_SUM, 0, MPI_COMM_WORLD);
				}
#endif
      }
      else{
				totloss_full[l]=totloss[l];
				for(k=0;k<6;k++){
					sumloss_full[l][k]=sumloss[l][k];
				}
      }
    }
  }

  /////////////////////////////////
  //
  // only dump loss data when wanted
  //
  //
  if((t!=told)&&( (t>=tloss)||(call_code==2)) ){
    told=t;
    if(myid<=0){
      // now output to file
      if(DETAILMLOSS==0){
				fprintf(loss_file," %21.15g",t);
				for(l=1;l<=NUMLOSSVAR;l++){
					if(REALNUMVEC!=NUMVEC){
						if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
							l=NUMSCA+NUMVEC*(3+1)+1;
						}
					}
	  
					// totloss_full[l] holds total loss so far of variable l on entire grid through boundaries.
					fprintf(loss_file," %21.15g",totloss_full[l]);
				}
				fprintf(loss_file,"\n");
				fflush(loss_file);
      }
      else if(DETAILMLOSS>=1){
				fprintf(loss_file," %21.15g",t);
				for(l=1;l<=NUMLOSSVAR;l++){
					if(REALNUMVEC!=NUMVEC){
						if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
							l=NUMSCA+NUMVEC*(3+1)+1;
						}
					}
					// sumloss_full[l][k] holds total loss so far of variable l on boundary k where k=0 is inner r-edge, k=1 outer r-edge, k=2 inner theta edge, k=3 outer theta edge
					fprintf(loss_file," %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g",totloss_full[l],sumloss_full[l][0],sumloss_full[l][1],sumloss_full[l][2],sumloss_full[l][3],sumloss_full[l][4],sumloss_full[l][5]);
				}
				fprintf(loss_file,"\n");
				fflush(loss_file);
      }
      
      if(DETAILMLOSS==2){
				j=1;
				for(k=0;k<2;k++){//i/o
					for(ii=0;ii<N3;ii++) for(i=0;i<N2;i++){//list
						fprintf(loss_file2," %21.15g",t);
						for(l=1;l<=NUMLOSSVAR;l++){

							if(l<=NUMSCA){ m=-1; typenum=0;}
							else if(l==NUMSCA+1){ m=0; typenum=1; }
							else if(l==NUMSCA+1+1){ m=1; typenum=1; }
							else if(l==NUMSCA+1+2){ m=2; typenum=1; }
							else if(l==NUMSCA+1+3){ m=3; typenum=1; }
							else if(l==NUMSCA+1*(3+1)+1){ m=0; typenum=2; }
							else if(l==NUMSCA+1*(3+1)+1+1){ m=1; typenum=2; }
							else if(l==NUMSCA+1*(3+1)+2+1){ m=2; typenum=2; }
							else if(l==NUMSCA+1*(3+1)+3+1){ m=3; typenum=2; }
							else if(l==NUMSCA+2*(3+1)+1){ m=-2; }
							if(REALNUMVEC!=NUMVEC){
								if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
									l=NUMSCA+NUMVEC*(3+1)+1;
									m=-2;  
								}
							}
	      
							// loss holds this dt transported variable l across boundary j=direction(1=x1, 2=x2, 3=x3) k(=0 inner =1, outer)
							if(m==-1) tempf=losss[l][j][k][ii][i];
							else if(m==-2) tempf=lossvisc[1][j][k][ii][i];
							else tempf=lossv[typenum][m][j][k][ii][i];
							fprintf(loss_file2," %7d %7d %7d %21.15g",k*N1,i,ii,tempf);
						}
						fprintf(loss_file2,"\n");
					}
				}
	
				j=2;
				for(k=0;k<2;k++){//i/o
					for(ii=0;ii<N3;ii++) for(i=0;i<N1;i++){//list
						fprintf(loss_file2," %21.15g",t);
						for(l=1;l<=NUMLOSSVAR;l++){

							if(l<=NUMSCA){ m=-1; typenum=0;}
							else if(l==NUMSCA+1){ m=0; typenum=1; }
							else if(l==NUMSCA+1+1){ m=1; typenum=1; }
							else if(l==NUMSCA+1+2){ m=2; typenum=1; }
							else if(l==NUMSCA+1+3){ m=3; typenum=1; }
							else if(l==NUMSCA+1*(3+1)+1){ m=0; typenum=2; }
							else if(l==NUMSCA+1*(3+1)+1+1){ m=1; typenum=2; }
							else if(l==NUMSCA+1*(3+1)+2+1){ m=2; typenum=2; }
							else if(l==NUMSCA+1*(3+1)+3+1){ m=3; typenum=2; }
							else if(l==NUMSCA+2*(3+1)+1){ m=-2; }
							if(REALNUMVEC!=NUMVEC){
								if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
									l=NUMSCA+NUMVEC*(3+1)+1;
									m=-2;  
								}
							}
	      
							if(m==-1) tempf=losss[l][j][k][ii][i];
							else if(m==-2) tempf=lossvisc[1][j][k][ii][i];
							else tempf=lossv[typenum][m][j][k][ii][i];
							fprintf(loss_file2," %7d %7d %7d %21.15g",i,k*N2,ii,tempf);
						}
						fprintf(loss_file2,"\n");
					}
				}
				j=3;
				for(k=0;k<2;k++){//i/o
					for(ii=0;ii<N2;ii++) for(i=0;i<N1;i++){//list
						fprintf(loss_file2," %21.15g",t);
						for(l=1;l<=NUMLOSSVAR;l++){

							if(l<=NUMSCA){ m=-1; typenum=0;}
							else if(l==NUMSCA+1){ m=0; typenum=1; }
							else if(l==NUMSCA+1+1){ m=1; typenum=1; }
							else if(l==NUMSCA+1+2){ m=2; typenum=1; }
							else if(l==NUMSCA+1+3){ m=3; typenum=1; }
							else if(l==NUMSCA+1*(3+1)+1){ m=0; typenum=2; }
							else if(l==NUMSCA+1*(3+1)+1+1){ m=1; typenum=2; }
							else if(l==NUMSCA+1*(3+1)+2+1){ m=2; typenum=2; }
							else if(l==NUMSCA+1*(3+1)+3+1){ m=3; typenum=2; }
							else if(l==NUMSCA+2*(3+1)+1){ m=-2; }
							if(REALNUMVEC!=NUMVEC){
								if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
									l=NUMSCA+NUMVEC*(3+1)+1;
									m=-2;  
								}
							}
	      
							if(m==-1) tempf=losss[l][j][k][ii][i];
							else if(m==-2) tempf=lossvisc[1][j][k][ii][i];
							else tempf=lossv[typenum][m][j][k][ii][i];
							fprintf(loss_file2," %7d %7d %7d %21.15g",i,ii,k*N3,tempf);
						}
						fprintf(loss_file2,"\n");
					}
				}
				fflush(loss_file2);
      }
    }// end if cpu that writes
  }// end if time to output or final output
  // initialize loss data, so don't have to worry about how added in routines!
  // if any size is 1 in a certain direction, I assume one doesn't care about loss then(see init_loss())
  init_loss();
  firsttime=0;


  // cleanup if final call
  if(call_code==2){
    if(myid<=0){
      
      if(DOLOSSDIAG){
				fclose(loss_file);
				if(DETAILMLOSS==2){
					fclose(loss_file2);
				}
      }
    }// end over write cpu
  }// end if final call


  // we now have totloss_full for other functions
}


void diag_loss_gen(int call_code,FTYPE tener,FTYPE tloss, SFTYPE* totloss_full)
	// call_code
	// -1: do every call to diag(usually each time step as called by main.c)
	// 0: first call setup/dump
	// 1: normal diag dump of variables each DTl (DT for logging)
	// 2: final diag call after done
     
{
  static SFTYPE told=-1000.00;
  fpos_t fposgod;
  SFTYPE tempf,tempf2;
  SFTYPE tempfns;
  SFTYPE masstempc,masstemp1,masstemp2,masstemp3;
  char dfnam[MAXFILENAME];
  char dfheader[MAXFILENAME];
  char dfnamtemp[MAXFILENAME];
  char dfnamback[MAXFILENAME];

  int i,j,k,l,m,ii,ll ;
  int dumi[10];
  SFTYPE tcheck; // used to check time for loss append

  int which, comp;
  SFTYPE ftemp,ftemp1,ftemp2,ftemp3;
  
  static int dump_cnt,pdump_cnt,im_cnt ; // used to enumerate files
  static SFTYPE totloss[NUMLOSSVAR+1]; // 0: etot
  static SFTYPE lossflux_full[NUMLOSSVAR+1][2]; // 0: etot
  static FILE *loss_file;
  static FILE *loss_file_temp;
  static int firsttime=1;
#if(USEMPI)
  static MPI_Request request[6]; 
#endif
  int typenum;
  char temps[MAXFILENAME];
  char filename[MAXFILENAME];
  char filenametemp[MAXFILENAME];
  char filenameback[MAXFILENAME];
  long fpos0;
  int gotit;
  // for easy coord conv

  if(told==t) return; // to avoid duplicate data output

  // tvdlf has no flux counters yet

  if(firsttime==1){ // the check for if first time here
    for(i=1;i<=NUMLOSSVAR;i++){
      totloss[i]=0.0;
    }
    
    if(myid<=0){
      
      if(DETAILMLOSS>=0){
	
				sprintf(temps,DATADIR);
	
				sprintf(dfnam,"%s0_loss%s",temps,DAT2EXT);
				if((loss_file = fopen(dfnam,WRITETYPE))==NULL) {
					fprintf(fail_file,"error opening loss output file %s\n",dfnam) ;
					myexit(1) ;
				}
	
				if((restartloss==0)||(appendold==0)){
					// version header
					fprintf(loss_file,"#%10s\n%10d %10d\n","LOSSVER",LOSSVER,LOSSTYPE);
					if(DETAILMLOSS==0){
						if(REALNUMVEC==NUMVEC){	  
							fprintf(loss_file,"#%21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s\n","time","Total Mass Loss","Total IntEnergy Loss","Total GravEnergy Loss","Total ke Loss","Total mx1 Loss","Total mx2 Loss","Total mx3 Loss","Total be Loss","Total bx1 Loss","Total bx2 Loss","Total bx3 Loss","Total ViscE Loss","Total ax1 Loss","Total ax2 Loss","Total ax3 Loss");
						}
						else{
							fprintf(loss_file,"#%21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s\n","time","Total Mass Loss","Total IntEnergy Loss","Total GravEnergy Loss","Total ke Loss","Total mx1 Loss","Total mx2 Loss","Total mx3 Loss","Total ViscE Loss","Total ax1 Loss","Total ax2 Loss","Total ax3 Loss");
						}
					}
					else if(DETAILMLOSS>=1){
						if(REALNUMVEC==NUMVEC){
							fprintf(loss_file,"#%21s"
											" %21s %21s" // 1
											" %21s %21s"
											" %21s %21s" // 3
											" %21s %21s"
											" %21s %21s"
											" %21s %21s" // 6
											" %21s %21s"
											" %21s %21s"
											" %21s %21s" // 9
											" %21s %21s"
											" %21s %21s"
											" %21s %21s" // 12
											" %21s %21s"
											" %21s %21s"
											" %21s %21s" // 15
											,"time"
											,"Mass Loss Inner"    ,"Outer"
											,"IEnergy Loss Inner","Outer"
											,"PEnergy Loss Inner","Outer"
											,"KEnergy Loss Inner","Outer"
											,"x1-Mom Loss Inner" ,"Outer"
											,"x2-Mom Loss Inner" ,"Outer"
											,"x3-Mom Loss Inner" ,"Outer"
											,"MEnergy Loss Inner","Outer"
											,"x1-Bflux Loss Inner" ,"Outer"
											,"x2-Bflux Loss Inner" ,"Outer"
											,"x3-Bflux Loss Inner","Outer"
											,"visc-e Loss Inner" ,"Outer"
											,"x1-AMom Loss Inner" ,"Outer"
											,"x2-AMom Loss Inner" ,"Outer"
											,"x3-AMom Loss Inner" ,"Outer"

								);
						}
						else{
							fprintf(loss_file,"#%21s"
											" %21s %21s"
											" %21s %21s"
											" %21s %21s"
											" %21s %21s"
											" %21s %21s"
											" %21s %21s"
											" %21s %21s"
											" %21s %21s"
											" %21s %21s"
											" %21s %21s"
											" %21s %21s"
											,"time"
											,"Mass Loss Inner","Outer"
											,"IEnergy Loss Inner","Outer"
											,"PEnergy Loss Inner","Outer"
											,"KEnergy Loss Inner","Outer"
											,"x1-Mom Loss Inner" ,"Outer"
											,"x2-Mom Loss Inner" ,"Outer"
											,"x3-Mom Loss Inner" ,"Outer"
											,"visc-e Loss Inner" ,"Outer"
											,"x1-AMom Loss Inner" ,"Outer"
											,"x2-AMom Loss Inner" ,"Outer"
											,"x3-AMom Loss Inner" ,"Outer"
								);
						}
						fprintf(loss_file,"\n");
					}
					fflush(loss_file);
				}//if appendold==0
				else{//if appendold==1
					if(myid<=0){
						fprintf(logfull_file,"Start setup of loss file append\n");
						fflush(logfull_file);
					}
					// need to read in totloss and sumloss and stick on right cpu
					// need to make sure getting right time.  If DTloss<DTd||DTi, then should be able to get easily.  You should generally have this true anyways
	  
					rewind(loss_file); // go to start
					// check version info
					while(fgetc(loss_file)!='\n'); // skip comment line
					fscanf(loss_file,"%d %d",&dumi[0],&dumi[1]);
					if( (dumi[0]!=LOSSVER)||(dumi[1]!=LOSSTYPE) ){
						fprintf(fail_file,"Expected lossver/losstype: %d %d got %d %d\n",LOSSVER,LOSSTYPE,dumi[0],dumi[1]);
						myexit(6);
					}
					while(fgetc(loss_file)!='\n'); // skip to next line
					while(fgetc(loss_file)!='\n'); // skip comment line	      
					gotit=0;
					while( (!feof(loss_file))&&(gotit==0) ){
	    
						fscanf(loss_file,INPUT3,&tcheck);
						//fprintf(stderr,"%15.10g %15.10g %d\n",t,tcheck,gotit);
						if(fabs(tcheck-t)<0.5*DTloss){
							gotit=1;
							for(l=1;l<=NUMLOSSVAR;l++){
								if(REALNUMVEC!=NUMVEC){
									if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
										l=NUMSCA+NUMVEC*(3+1)+1;
									}
								}
								if(DETAILMLOSS==0){
									fscanf(loss_file,INPUT60,&totloss_full[l]);
								}
								else if(DETAILMLOSS==1){
									fscanf(loss_file,INPUT61,&lossflux_full[l][0],&lossflux_full[l][1]);
								}
							}
						}
						else{
							while((fgetc(loss_file)!='\n')&&(!feof(loss_file))); // skip this bad line
						}
						// continue after successful get since successful get is good data and should keep since corresponds to dump one is keeping
						fpos0=ftell(loss_file); // position to continue writting at if successful get
					}
					if(gotit==0){
						fprintf(fail_file,"Never found right time in loss file when appending: looking for t=%21.15g lastt=%21.15g\n",timereenter,tcheck);
						myexit(1);
					}
					else{
						fprintf(logfull_file,"found goodtime t=%21.15g (wanted %21.15g) to restart loss file\n",tcheck,timereenter);
						sprintf(temps,DATADIR);
						sprintf(dfnamtemp,"%s0_loss%s.temp",temps,DAT2EXT);
						sprintf(dfnam,"%s0_loss%s",temps,DAT2EXT);
						sprintf(dfnamback,"%s0_loss%s.back",temps,DAT2EXT);
	    
						// now that done, fix up file
						if( (loss_file_temp=fopen(dfnamtemp,"wt"))==NULL){
							fprintf(fail_file,"Cannot open temp loss file for appending: %s\n",dfnamtemp);
							myexit(1);
						}
						else{
							rewind(loss_file);
							while(ftell(loss_file)<fpos0+1){
								fputc(fgetc(loss_file),loss_file_temp);
							}
							fclose(loss_file_temp);
							fclose(loss_file);
							rename(dfnam,dfnamback); // move old to backup location
							rename(dfnamtemp,dfnam); // move new to old name(normal name)
							// reopen loss_file (now normal name)
							if((loss_file = fopen(dfnam,"at"))==NULL) {
								fprintf(fail_file,"2: error opening loss output file %s\n",dfnam) ;
								myexit(1) ;
							}
							if(myid<=0){
								fprintf(logfull_file,"End setup of loss file append\n");
								fflush(logfull_file);
							}
						}
					}
				}//end else if appendold==1
      }// end if doing normal mloss output
    }// end if write cpu
    // no distribution to other cpus necessary
  }// end if first time in here      
  
  
  // must compute sum over cpus everytime since needed for diagnostics on DTl period, or for loss diagnostics
  if((t>=tener) ||(t>=tloss)||(call_code==2)){

    for(l=1;l<=NUMLOSSVAR;l++){

      if(REALNUMVEC!=NUMVEC){
				if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
					l=NUMSCA+NUMVEC*(3+1)+1;
					m=-2;  
				}
      }
      // sum up for this cpu
      totloss[l]=lossflux[l][0]+lossflux[l][1];

      // sum up over cpus
      if(numprocs>1){
#if(USEMPI)
				MPI_Reduce(&totloss[l], &totloss_full[l], 1, MPI_SFTYPE, MPI_SUM, 0, MPI_COMM_WORLD);

				MPI_Reduce(&lossflux[l][0], &lossflux_full[l][0], 1, MPI_SFTYPE, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(&lossflux[l][1], &lossflux_full[l][1], 1, MPI_SFTYPE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
      }
      else{
				totloss_full[l]=totloss[l];
				lossflux_full[l][0]=lossflux[l][0];
				lossflux_full[l][1]=lossflux[l][1];
      }
    }
  }

  // only dump loss data when wanted
  if((t>=tloss)||(call_code==2)){
    if(myid<=0){
      // now output to file
      if(DETAILMLOSS==0){
				fprintf(loss_file," %21.15g",t);
				for(l=1;l<=NUMLOSSVAR;l++){
					if(REALNUMVEC!=NUMVEC){
						if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
							l=NUMSCA+NUMVEC*(3+1)+1;
						}
					}
	  
					// totloss_full[l] holds total loss so far of variable l on entire grid through boundaries.
					fprintf(loss_file," %21.15g",totloss_full[l]);
				}
				fprintf(loss_file,"\n");
				fflush(loss_file);
      }
      else if(DETAILMLOSS>=1){
				fprintf(loss_file," %21.15g",t);
				for(l=1;l<=NUMLOSSVAR;l++){
					if(REALNUMVEC!=NUMVEC){
						if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
							l=NUMSCA+NUMVEC*(3+1)+1;
						}
					}
					fprintf(loss_file," %21.15g %21.15g",lossflux_full[l][0],lossflux_full[l][1]);
				}
				fprintf(loss_file,"\n");
				fflush(loss_file);
      }
    }// end if cpu that writes
  }// end if time to output or final output


  // cleanup if final call
  if(call_code==2){
    if(myid<=0){
      if(DOLOSSDIAG){
				fclose(loss_file);
      }
    }// end over write cpu
  }// end if final call


  // we now have totloss_full for other functions
  firsttime=0;
  told=t;
}




void diag_ener(int call_code,SFTYPE*totloss_full)
	// call_code
	// -1: do every call to diag(usually each time step as called by main.c)
	// 0: first call setup/dump
	// 1: normal diag dump of variables each DTl (DT for logging)
	// 2: final diag call after done
     
{
  static SFTYPE told=-1000.00;
  fpos_t fposgod;
  SFTYPE tempf,tempf2;
  SFTYPE tempfns;
  SFTYPE masstempc,masstemp1,masstemp2,masstemp3;
  char dfnam[MAXFILENAME];
  char dfheader[MAXFILENAME];
  char dfnamtemp[MAXFILENAME];
  char dfnamback[MAXFILENAME];

  SFTYPE mass_full,eth_full,ek_full,eg_full,eb_full,se_full,smom_full[3+1],angmom_full[3+1],magflux_full[3+1],cmode_amp_full,smode_amp_full,floors_full[NUMLOSSVAR+1],inflows_full[NUMLOSSVAR+1],radiations_full[NUMLOSSVAR+1]; // 0: etot
  SFTYPE ek,eth,eg,eb,se,cmode_amp,smode_amp,mass,smom[3+1],angmom[3+1],magflux[3+1]; // 0: etot

  SFTYPE b2pole,b2pole_full;
  SFTYPE polevolume,polevolume_full;

  SFTYPE vx2i,vy2i,vz2i,bx2i,by2i,bz2i;
  SFTYPE bx20_full,by20_full,bz20_full;// mean initial field
  SFTYPE bx2,bx2_full,by2,by2_full,bz2,bz2_full;
  SFTYPE vx2,vx2_full,vy2,vy2_full,vz2,vz2_full;

  SFTYPE eki,ebi,ethi,egi ;
  SFTYPE vxa,vya,vza ;
  SFTYPE Bxa,Bya,Bza ;
  SFTYPE dxdyc,dxdy1,dxdy2,dxdy3;
  FILE *efboth;
  int i,j,k,l,m,ii,ll ;
  int loopdido;
  int dumi[10];
  SFTYPE tcheck; // used to check time for loss append
  
  int which, comp;
  
  static SFTYPE varinit[NUMLOSSVAR+1]; // 0: etot
  static SFTYPE varfinal[NUMLOSSVAR+1]; // 0: etot
  
  static char losstext[NUMLOSSVAR+1][50]; // 0: etot
  SFTYPE realloss; // need for enthalpy and sm
  SFTYPE ftemp,ftemp1,ftemp2,ftemp3;
  SFTYPE floattemp[NUMLOSSVAR+1][NUMINPUTS+1]; // NUMINPUTS types of output for 0_ener.dat per variable


  static FILE *ener_file;
  static FILE *ener_file_temp;
  static FILE *final_output;
  static int timeless=0;
  static int firsttime=1;
#if(USEMPI)
  static MPI_Request request[6]; 
#endif
  int typenum;
  char temps[MAXFILENAME];
  char filename[MAXFILENAME];
  char filenametemp[MAXFILENAME];
  char filenameback[MAXFILENAME];
  long fpos0;
  int gotit;
  // for easy coord conv
  SFTYPE posx,posy,posz;
  SFTYPE radius,theta,phi;
  // for collaboration
  SFTYPE IPLUS,IMINUS,MASSPLUS,MASSMINUS,VOLPLUS,VOLMINUS,IPLUS_full,IMINUS_full,MASSPLUS_full,MASSMINUS_full,VOLPLUS_full,VOLMINUS_full;

  if(t==told) return; // since need to avoid duplicate data output


  // SETUP ENERGY DIAGNOSTICS
  
  if(firsttime==1){
    
    
    if(myid<=0){
      // setup labels for each varinit[] and varfinal[]
      strcpy(losstext[0],"ETotal");
      strcpy(losstext[1],"Mass");
      strcpy(losstext[2],"IntE");
      strcpy(losstext[3],"PotE");
      strcpy(losstext[4],"KinE");
      strcpy(losstext[5],"SMom1");
      strcpy(losstext[6],"SMom2");
      strcpy(losstext[7],"SMom3");
      strcpy(losstext[8],"MagE");
      strcpy(losstext[9],"Magflux1");
      strcpy(losstext[10],"Magflux2");
      strcpy(losstext[11],"Magflux3");
      strcpy(losstext[12],"ViscEn");
      strcpy(losstext[13],"AngMom1");
      strcpy(losstext[14],"AngMom2");
      strcpy(losstext[15],"AngMom3");

      
      sprintf(temps,DATADIR);
      sprintf(filename,"%s0_ener%s",temps,DAT2EXT);
      if((ener_file = fopen(filename,WRITETYPE))==NULL) {
				fprintf(fail_file,"error opening energy output file %s\n",filename) ;
				myexit(1) ;
      }

      sprintf(temps,DATADIR);	
      sprintf(filename,"%s0_final%s",temps,DAT2EXT);
      if((final_output = fopen(filename,"w"))==NULL) {
				fprintf(fail_file,"error opening final output file %s\n",filename) ;
				myexit(1) ;
      }
      
      
      if((restartener==0)||(appendold==0)){
				// volume integral file:
				// version header
				fprintf(ener_file,"#%10s\n%10d %10d\n","ENERVER",ENERVER,ENERTYPE);
				fprintf(ener_file,"#%21s","t") ;
				for(l=0;l<=NUMLOSSVAR;l++){
					if(REALNUMVEC!=NUMVEC){
						if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
							l=NUMSCA+NUMVEC*(3+1)+1;
						}
					}
					fprintf(ener_file," %20s%1d %20s%1d %20s%1d %20s%1d %20s%1d %20s%1d %20s%1d",losstext[l],l,"Totalchange",l,"b_lost",l,"fl_add",l,"inject_add",l,"radiate",l,"totdel+bloss-fladd",l);
				}
				fprintf(ener_file,"\n");
      }
      else{//if appendold==1


				// volume integral file append
				if(myid<=0){
					fprintf(logfull_file,"Start setup of energy file append\n");
					fflush(logfull_file);
				}
				// need to read in first(or any) data line of 0_ener.dat file for varinit[l] values
	
				rewind(ener_file); // go to start
				// check on version info
				while(fgetc(ener_file)!='\n'); // skip comment line	      
				fscanf(ener_file,"%d %d\n",&dumi[0],&dumi[1]);
				if((dumi[0]!=ENERVER)||(dumi[1]!=ENERTYPE) ){
					fprintf(fail_file,"Expected enerver/enertype: %d %d got %d %d\n",ENERVER,ENERTYPE,dumi[0],dumi[1]);
					myexit(6);
				}

				while(fgetc(ener_file)!='\n'); // skip comment line	
				gotit=0;
				while( (!feof(ener_file))&&(gotit==0) ){
	  
					fscanf(ener_file,INPUT4,&tcheck);
					//fprintf(stderr,"%21.15g %21.15g\n",tcheck,t);
					if(fabs(tcheck-timereenter)<0.5*DTener){
						gotit=1;
						// read in init values
						for(l=0;l<=NUMLOSSVAR;l++){
							if(REALNUMVEC!=NUMVEC){
								if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
									l=NUMSCA+NUMVEC*(3+1)+1;
								}
							}
							fscanf(ener_file,INPUT5,&floattemp[l][0],&floattemp[l][1],&floattemp[l][2],&floattemp[l][3],&floattemp[l][4],&floattemp[l][5],&floattemp[l][6]);
							// now compute init values from last value and difference
							// varfinal computed already above, and correct
							varfinal[l]=varinit[l]=floattemp[l][0]-floattemp[l][1];
							// floors/inflows isn't really just on this CPU, but the final counting doesn't care
							// and _full version is correct(and overrides collective above), and needed for later calcs
							floors[l]=floors_full[l]=floattemp[l][3];
							inflows[l]=inflows_full[l]=floattemp[l][4];
							radiations[l]=radiations_full[l]=floattemp[l][5];
						}
						while((fgetc(ener_file)!='\n')&&(!feof(ener_file))); // skip rest of line
						fpos0=ftell(ener_file); // position to continue writting at if successful get
					}// if good time
					else{
						while((fgetc(ener_file)!='\n')&&(!feof(ener_file))); // skip this bad line
					}
				}// while finding good time or hitting end of file
				if(gotit==0){
					fprintf(fail_file,"Never found right time in energy file when appending: looking for t=%21.15g lastt=%21.15g\n",timereenter,tcheck);
					myexit(1);
				}
				else{
					fprintf(logfull_file,"found goodtime t=%21.15g (wanted %21.15g) to restart energy file (current time=%15.10g)\n",tcheck,timereenter,t);
					sprintf(filename,"%s0_ener%s",DATADIR,DAT2EXT);
					sprintf(filenametemp,"%s0_ener%s.temp",DATADIR,DAT2EXT);
					sprintf(filenameback,"%s0_ener%s.back",DATADIR,DAT2EXT);
	  
					if((ener_file_temp = fopen(filenametemp,"wt"))==NULL) {
						fprintf(fail_file,"error opening temp energy output file %s\n",filenametemp) ;
						myexit(1) ;
					}
					else{
						rewind(ener_file);
						while(ftell(ener_file)<fpos0){ // 0 through fpos0-1
							fputc(fgetc(ener_file),ener_file_temp);
						}
						fclose(ener_file_temp);
						fclose(ener_file);
						rename(filename,filenameback); // move old to backup location
						rename(filenametemp,filename); // move new to old(normal)
						// reopen ener_file
						if((ener_file = fopen(filename,"at"))==NULL) {
							fprintf(fail_file,"error opening energy output file %s\n",filename) ;
							myexit(1) ;
						}
						if(myid<=0){
							fprintf(logfull_file,"End setup of energy file append\n");
							fflush(logfull_file);
						}
					}
				}// end else if gotit==1 on ener file

      }// end else if appendold==1
    }// end if myid<=0
  }// end if firsttime

  // generally do:
  
  // initialize variables for this dt step
  mass=0;
  smom[1]=0;
  smom[2]=0;
  smom[3]=0;
  angmom[1]=0;
  angmom[2]=0;
  angmom[3]=0;
  magflux[1]=0;
  magflux[2]=0;
  magflux[3]=0;
  
  /* calculate energies */
  se = 0. ;
  eb = 0. ;
  ek = 0. ;
  eth = 0. ;
  eg = 0. ;
  b2pole=0.;
  polevolume=0.;
  bx2=by2=bz2=0;
  vx2=vy2=vz2=0;

  // collaboration
  IPLUS=0;
  IMINUS=0;
  MASSPLUS=0;
  MASSMINUS=0;
  VOLPLUS=0;
  VOLMINUS=0;
  

#define ALLZONECENTERED 1
  // ==0 not really right since vectors offset from scalars
  // ==1 offset right, but not precise for machine conserved quantities such as ang mom flux in spc
  
  LOOPINT{    // begin loop over real domain to get diagnostics
    
    vxa = e2z_1(v[1][1],k,j,i);      
    vya = e2z_2(v[1][2],k,j,i);
    vza = e2z_3(v[1][3],k,j,i) ;
    
    eki = 0.5*s[1][k][j][i]*(vxa*vxa + vya*vya + vza*vza) ;
    
    Bxa = e2z_1(v[2][1],k,j,i);
    Bya = e2z_2(v[2][2],k,j,i);
    Bza = e2z_3(v[2][3],k,j,i);
    
    ebi = 0.5*(Bxa*Bxa + Bya*Bya + Bza*Bza) ;

    if( (analoutput==13)||(analoutput==12)){
      vx2i=0.5*s[1][k][j][i]*vxa*vxa;
      vy2i=0.5*s[1][k][j][i]*vya*vya;
      vz2i=0.5*s[1][k][j][i]*vza*vza;

      // this creeps?
      //      bx2i=0.5*Bxa*Bxa;
      //by2i=0.5*Bya*Bya;
      //bz2i=0.5*Bza*Bza;

      // differential energy : field component of wave energy (1.0 here is mean Bx field
      bx2i=0.5*(v[2][1][k][j][i]-1.0)*(v[2][1][k][j][i]-1.0);
      by2i=0.5*(v[2][2][k][j][i])*(v[2][2][k][j][i]);
      bz2i=0.5*(v[2][3][k][j][i])*(v[2][3][k][j][i]);
    }
    
    //egi = 0.5*s[1][k][j][i]*s[3][k][j][i] ; // for self-gravity
    egi = s[1][k][j][i]*s[3][k][j][i] ;
    
    /* Equation of state */
    if(press==1){
      if(wgam) ethi = s[2][k][j][i] ;
      else ethi = cs*cs*s[1][k][j][i]*log(s[1][k][j][i]) ;
    }
    else{
      ethi=0;
    }
    
    dxdyc=DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k) ;    
    masstempc = s[1][k][j][i]*dxdyc; // mass in center of a zone

    // average polar field b^2
    // uses same as sm average
    // for BZ-effect estimation (only directly around black hole)
#define B2THETARANGE (0.25*M_PI)
#define B2RRANGE (1.1*L[1][1])
 
    if((analoutput==5)||(analoutput==14)){
      if(COORD==3){
				if((fabs(x[1][2][j]-0.5*M_PI)>=B2THETARANGE)&&(x[1][1][i]>=L[1][1])&&((x[1][1][i]<=B2RRANGE)||(L[1][1]+dx[2][1][0]>=B2RRANGE)) ){
					// account if either within L[1][1] to 1.1*L[1][1] or 
					b2pole+=(Bxa*Bxa+Bya*Bya+Bza*Bza)*dxdyc;
					polevolume+=dxdyc;
				}
      }
      else if(COORD==1){
				radius=cart2spc(1,x[1][1][i],x[1][2][j],x[1][3][k]);
				theta=cart2spc(2,x[1][1][i],x[1][2][j],x[1][3][k]);
				//phi=cart2spc(3,x[1][1][i],x[1][2][j],x[1][3][k]);
				if((fabs(theta-0.5*M_PI)>=0.25*M_PI)&&(radius>Rinner)&&(radius<1.1*Rinner) ){
					b2pole+=(Bxa*Bxa+Bya*Bya+Bza*Bza)*dxdyc;
					polevolume+=dxdyc;
				}
      }
    }

    // sum up the zone values to totals
    ek  += eki*dxdyc ;
    eb  += ebi*dxdyc ;
    eth += ethi*dxdyc ;
    eg  += egi*dxdyc ;
    se+=s[2][k][j][i]/s[1][k][j][i];
    mass += masstempc;

    if( (analoutput==13)||(analoutput==12)){
      vx2+=vx2i*dxdyc;
      vy2+=vy2i*dxdyc;
      vz2+=vz2i*dxdyc;

      bx2+=bx2i*dxdyc;
      by2+=by2i*dxdyc;
      bz2+=bz2i*dxdyc;
    }

#if(ALLZONECENTERED==0)
    dxdy1=DVL(2,1,i)*DVL(1,2,j)*DVL(1,3,k) ;
    dxdy2=DVL(1,1,i)*DVL(2,2,j)*DVL(1,3,k) ;
    dxdy3=DVL(1,1,i)*DVL(1,2,j)*DVL(2,3,k) ;
    masstemp1 = z2e_1(s[1],k,j,i)*dxdy1; // mass in zone around edge 1
    masstemp2 = z2e_2(s[1],k,j,i)*dxdy2; // mass in zone around edge 2
    masstemp3 = z2e_3(s[1],k,j,i)*dxdy3; // mass in zone around edge 3

    smom[1] += masstemp1*v[1][1][k][j][i];
    smom[2] += masstemp2*v[1][2][k][j][i]*G2(2,i);
    smom[3] += masstemp3*v[1][3][k][j][i]*G3(2,i)*G4(2,j);
    
    magflux[1] += v[2][1][k][j][i];
    magflux[2] += v[2][2][k][j][i]/G2(2,i);
    magflux[3] += v[2][3][k][j][i]/(G3(2,i)*G4(2,j));
#else
    smom[1] += masstempc*vxa;
    smom[2] += masstempc*vya*G2(2,i);
    smom[3] += masstempc*vza*G3(2,i)*G4(2,j);
    
    magflux[1] += Bxa;
    magflux[2] += Bya/G2(2,i);
    magflux[3] += Bza/(G3(2,i)*G4(2,j));
#endif

#if(COORD==3)
    angmom[1] += 0; // not yet
    angmom[2] += 0; // not yet
    angmom[3] += masstempc*vza*G3(2,i)*G4(2,j); // about z
#elif(COORD==1)
    posx=x[2][1][i];
    posy=x[2][2][j];
    posz=x[2][3][k];
    angmom[1] += masstempc*sqrt(vza*vza*posy*posy+vya*vya*posz*posz); 
    angmom[2] += masstempc*sqrt(vxa*vxa*posz*posz+vza*vza*posx*posx); 
    angmom[3] += masstempc*sqrt(vxa*vxa*posy*posy+vya*vya*posx*posx); 
#endif

    if((RMAX!=0.0)&&((analoutput==5)||(analoutput==14))){ // then doing mode analysis for anti-symmetric epicyclic mode

      if(COORD==3){
				//	fprintf(stderr,"%15.10g %15.10g %15.10g\n",RMAX-DELTAR,x[1][1][i],RMAX+DELTAR); fflush(stderr);
				if((x[1][1][i]>=RMAX-DELTAR)&&(x[1][1][i]<=RMAX+DELTAR)){
	  
					if((x[1][2][j]>M_PI*0.5)&&(x[1][2][j]<=M_PI*0.5+2.0*HOR)){ // IPLUS
	    
						IPLUS+=s[1][k][j][i]*vxa*dxdyc;
						VOLPLUS+=dxdyc;
						MASSPLUS+=s[1][k][j][i]*dxdyc;
					}
					else if((x[1][2][j]<M_PI*0.5)&&(x[1][2][j]>=M_PI*0.5-2.0*HOR)){ // IMINUS
						IMINUS+=s[1][k][j][i]*vxa*dxdyc;
						VOLMINUS+=dxdyc;
						MASSMINUS+=s[1][k][j][i]*dxdyc;
					}
				}
      }
      else{
				fprintf(fail_file,"no such method yet\n");
				myexit(1);
      }
    }// end if doing collaboration computation


  }// end of volume loop

  
  // add up volume integrals over all cpus
  if(numprocs>1){
#if(USEMPI)
    MPI_Reduce(&mass, &mass_full, 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&eth,  &eth_full, 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&eg,   &eg_full, 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&ek, &ek_full, 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&eb,   &eb_full, 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&se,   &se_full, 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    
    MPI_Reduce(&(angmom[1]), &(angmom_full[1]), 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&(angmom[2]), &(angmom_full[2]), 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&(angmom[3]), &(angmom_full[3]), 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&(smom[1]), &(smom_full[1]), 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&(smom[2]), &(smom_full[2]), 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&(smom[3]), &(smom_full[3]), 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&(magflux[1]), &(magflux_full[1]), 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&(magflux[2]), &(magflux_full[2]), 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&(magflux[3]), &(magflux_full[3]), 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    for(l=1;l<=NUMLOSSVAR;l++){
      MPI_Reduce(&(floors[l]),   &(floors_full[l]), 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
      MPI_Reduce(&(inflows[l]),   &(inflows_full[l]), 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
      MPI_Reduce(&(radiations[l]),   &(radiations_full[l]), 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    }

    MPI_Reduce(&b2pole,   &b2pole_full, 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&polevolume,   &polevolume_full, 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);

    MPI_Reduce(&vx2,   &vx2_full, 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&vy2,   &vy2_full, 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&vz2,   &vz2_full, 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);

    MPI_Reduce(&bx2,   &bx2_full, 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&by2,   &by2_full, 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&bz2,   &bz2_full, 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);

    MPI_Reduce(&IPLUS,    &IPLUS_full,     1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&IMINUS,   &IMINUS_full,    1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&VOLPLUS,  &VOLPLUS_full,   1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&VOLMINUS, &VOLMINUS_full,  1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&MASSPLUS,  &MASSPLUS_full,   1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    MPI_Reduce(&MASSMINUS, &MASSMINUS_full,  1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
    
    
#endif
  }
  else{
    mass_full=mass;
    eth_full=eth;
    eg_full=eg;
    ek_full=ek;
    eb_full=eb;
    se_full=se;
    smom_full[1]=smom[1];
    smom_full[2]=smom[2];
    smom_full[3]=smom[3];
    angmom_full[1]=angmom[1];
    angmom_full[2]=angmom[2];
    angmom_full[3]=angmom[3];
    magflux_full[1]=magflux[1];
    magflux_full[2]=magflux[2];
    magflux_full[3]=magflux[3];
    for(l=1;l<=NUMLOSSVAR;l++){
      floors_full[l]=floors[l];
      inflows_full[l]=inflows[l];
      radiations_full[l]=radiations[l];
    }
    b2pole_full=b2pole;
    polevolume_full=polevolume;

    vx2_full=vx2;
    vy2_full=vy2;
    vz2_full=vz2;

    bx2_full=bx2;
    by2_full=by2;
    bz2_full=bz2;

    IPLUS_full=IPLUS;
    IMINUS_full=IMINUS;
    VOLPLUS_full=VOLPLUS;
    VOLMINUS_full=VOLMINUS;
    MASSPLUS_full=MASSPLUS;
    MASSMINUS_full=MASSMINUS;
  }


  if( (analoutput==13)||(analoutput==12)){
    if(call_code==0){
      bx20_full=bx2_full;
      by20_full=by2_full;
      bz20_full=bz2_full;
    }
  }


  if(myid<=0){
		if((analoutput==5)||(analoutput==14)){
			// divide out polar volume
			if(polevolume_full>0.0)    b2pole_full/=polevolume_full;
			else{
	      b2pole_full=0.0;
	      fprintf(logfull_file,"warning: 0 volume for b2omdot calculation\n");
			}
		}
  }
  
  // now write results using cpu=0 (which received the totals)
  if(myid<=0){
    if(call_code == 0) {
      if((restartener==0)||(appendold==0)){
				// for call_code==0, varfinal final so far!
				varinit[1] = mass_full+SSMALL ;
				varinit[2] = eth_full+SSMALL ;
				varinit[3] = eg_full+SSMALL ;
				varinit[4] = ek_full+SSMALL ;
				varinit[5] = smom_full[1]+SSMALL ;
				varinit[6] = smom_full[2]+SSMALL ;
				varinit[7] = smom_full[3]+SSMALL ;
				varinit[8] = eb_full+SSMALL ;
				varinit[9] = magflux_full[1]+SSMALL ;
				varinit[10] = magflux_full[2]+SSMALL ;
				varinit[11] = magflux_full[3]+SSMALL ;
				varinit[13] = angmom_full[3]+SSMALL ;
				varinit[14] = angmom_full[3]+SSMALL ;
				varinit[15] = angmom_full[3]+SSMALL ;
				// sum up all energy terms from totals
				varinit[0]=varinit[2]+varinit[3]+varinit[4]+varinit[8]; // mass energy?
      }// end if callcode==0
    }
    
    // final so far, in any call_code
    varfinal[1] = mass_full+SSMALL ;
    varfinal[2] = eth_full+SSMALL ;
    varfinal[3] = eg_full+SSMALL ;
    varfinal[4] = ek_full+SSMALL ;
    varfinal[5] = smom_full[1]+SSMALL ;
    varfinal[6] = smom_full[2]+SSMALL ;
    varfinal[7] = smom_full[3]+SSMALL ;
    varfinal[8] = eb_full+SSMALL ;
    varfinal[9] = magflux_full[1]+SSMALL ;
    varfinal[10] = magflux_full[2]+SSMALL ;
    varfinal[11] = magflux_full[3]+SSMALL ;
    varfinal[13] = angmom_full[1]+SSMALL ;
    varfinal[14] = angmom_full[2]+SSMALL ;
    varfinal[15] = angmom_full[3]+SSMALL ;

    // sum up all energy terms from totals
    varfinal[0]=varfinal[2]+varfinal[3]+varfinal[4]+varfinal[8]; // mass energy?
    
    // need to do both final out and ener out when call_code==2
    if(call_code==2){
      loopdido=2;
    }
    else loopdido=1;
    
    for(i=1;i<=loopdido;i++){
      if(i==1) efboth=ener_file;
      else if(i==2) efboth=final_output;
      
      if(i==2) fprintf(efboth,"#%21s\n","t") ;
      fprintf(efboth," %21.15g",t) ;
      if(i==2) fprintf(efboth,"\n");
      if(i==2) fprintf(efboth,"#%21s %21s %21s %21s %21s %21s %21s %21s\n","what","value","tot_d","b_lost","fl_add","minject","radiate","totd+bl-fladd");
      
      // compute total energy variables from all data
      totloss_full[0]=totloss_full[2]+totloss_full[3]+totloss_full[4]+totloss_full[8]+totloss_full[12]; // [2] is really enthalpy as it should be
      floors_full[0]=floors_full[2]+floors_full[3]+floors_full[4];
      inflows_full[0]=inflows_full[2]+inflows_full[3]+inflows_full[4];
      radiations_full[0]=radiations_full[2]+radiations_full[3]+radiations_full[4];
      
      for(l=0;l<=NUMLOSSVAR;l++){
				if(REALNUMVEC!=NUMVEC){
					if(l==NUMSCA+1*(3+1)+1){ // skip B field for now, go directly to visc energy
						l=NUMSCA+NUMVEC*(3+1)+1;
					}
				}
				if(i==2) fprintf(efboth," %21s",losstext[l]);
				// below line only true for ideal EOS
				if(l==2) realloss=totloss_full[l]/gam; // totloss_full holds enthalpy across boundary, so change to ie
				else  realloss=totloss_full[l];
				if(fabs(realloss)<1.0E-99) realloss=SSMALL;// needed for sm to not die when using float version of sm(float sm can't handle >3 exponential digits)
				fprintf(efboth," %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g"
								,varfinal[l]
								,(varfinal[l]-varinit[l])
								//,totloss_full[l]
								,realloss
								,floors_full[l]
								,inflows_full[l]
								,radiations_full[l]
								//,(varfinal[l]+totloss_full[l]+radiations_full[l]-floors_full[l]-inflows_full[l]-varinit[l]));
								,(varfinal[l]+realloss+radiations_full[l]-floors_full[l]-inflows_full[l]-varinit[l])
					);
				if(i==2) fprintf(efboth,"\n");
      }
      /////////////
      //
      // now add in any extra stuff (don't treat as new version, just read in additional columns seperately
      // starts with 7*NUMLOSSVAR+1 th column (106 currently with version 7)
      // assumes need no restart of this either
      fprintf(efboth,
							" %ld %21.15g %21.15g"
							,nstep
							,divbmax_full // just last saved value from divb0check() (don't worry about forcing sync)
							,divbavg_full // just last saved value from divb0check() (don't worry about forcing sync)
	      );
      if((analoutput==5)||(analoutput==14)){
				fprintf(efboth,
								" %21.15g %21.15g"
								,b2pole_full
								,polevolume_full
					);
      }
      if((RMAX!=0.0)&&((analoutput==5)||(analoutput==14))){
				fprintf(efboth,
								" %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g"
								,IPLUS_full
								,IMINUS_full
								,VOLPLUS_full
								,VOLMINUS_full
								,MASSPLUS_full
								,MASSMINUS_full
					);
      }
      if( (analoutput==13)||(analoutput==12)){
				fprintf(efboth,
								" %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g"
								,vx2_full
								,vy2_full
								,vz2_full
								,bx2_full
								,by2_full
								,bz2_full
					);
				// just some random point on the grid
				k=N3/2;
				j=N2/2;
				i=N1/2;
				fprintf(efboth,
								" %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g"
								,s[1][k][j][i]
								,s[2][k][j][i]
								,v[1][1][k][j][i]
								,v[1][2][k][j][i]
								,v[1][3][k][j][i]
								,v[2][1][k][j][i]
								,v[2][2][k][j][i]
								,v[2][3][k][j][i]
					);
      }

      //
      //
      ///////////////
      fprintf(efboth,"\n");
      fflush(efboth) ;
      
    }// end over ener or final output of integrated variables
    
    
  }// end if cpu to write to file

  // cleanup if final call
  if(call_code==2){
    if(myid<=0){
      fclose(ener_file);
      fclose(final_output);      
    }// end over write cpu
  }// end if final call
  
  firsttime=0;
  told=t;
}



#define MODEMETHOD 2
// see below

void diag_mode(int call_code)
	// call_code
	// -1: do every call to diag(usually each time step as called by main.c)
	// 0: first call setup/dump
	// 1: normal diag dump of variables each DTl (DT for logging)
	// 2: final diag call after done
     
{
  static SFTYPE told=-1000.00;
  fpos_t fposgod;
  SFTYPE tempf,tempf2;
  SFTYPE tempfns;
  SFTYPE masstempc,masstemp1,masstemp2,masstemp3;
  char dfnam[MAXFILENAME];
  char dfheader[MAXFILENAME];
  char dfnamtemp[MAXFILENAME];
  char dfnamback[MAXFILENAME];

  SFTYPE mode_full[NUMSCA+NUMVEC*3+1][NUMMODES]; // m=0...NUMMODES-1
  SFTYPE mode[NUMSCA+NUMVEC*3+1][NUMMODES];
  int modeloop;

  SFTYPE modei1[NUMSCA+NUMVEC*3+1][NUMMODES];
  SFTYPE modei2[NUMSCA+NUMVEC*3+1][NUMMODES];

#define MAXSHELLS 50

  SFTYPE modei1cart[MAXSHELLS][MAXSHELLS][NUMSCA+NUMVEC*3+1][NUMMODES];
  SFTYPE modei2cart[MAXSHELLS][MAXSHELLS][NUMSCA+NUMVEC*3+1][NUMMODES];
  SFTYPE vxa,vya,vza ;
  SFTYPE Bxa,Bya,Bza ;
  SFTYPE dxdyc,dxdy1,dxdy2,dxdy3;
  int i,j,k,l,m,ii,ll ;
  int loopdido;
  int dumi[10];
  SFTYPE tcheck; // used to check time for loss append

  int which, comp;
  SFTYPE realloss; // need for enthalpy and sm
  SFTYPE ftemp,ftemp1,ftemp2,ftemp3;
  
  static FILE *mode_file;
  static FILE *mode_file_temp ;
  static int timeless=0;
  static int firsttime=1;
#if(USEMPI)
  static MPI_Request request[6]; 
#endif
  int typenum;
  char temps[MAXFILENAME];
  char filename[MAXFILENAME];
  char filenametemp[MAXFILENAME];
  char filenameback[MAXFILENAME];
  long fpos0;
  int gotit;
  SFTYPE Fun,Jac,dV,radius,theta,phi,TOTALVOLUME;
  int ri,thetai;
  // for easy coord conv



  if((COORD==3)&&(N3==1)) return; // since not relevant
  if(told==t) return; // to avoid duplicate data output


  if(firsttime==1){
    
    
    if(myid<=0){

      sprintf(temps,DATADIR);
      sprintf(filename,"%s0_mode%s",temps,DAT2EXT);
      if((mode_file = fopen(filename,WRITETYPE))==NULL) {
				fprintf(fail_file,"error opening mode output file %s\n",filename) ;
				myexit(1) ;
      }

      if((restartmode==0)||(appendold==0)){
				// mode file:
				fprintf(mode_file,"#%10s\n%10d %10d\n","MODEVER",MODEVER,MODETYPE);
				fprintf(mode_file,"#%21s","t") ;
				for(ll=1;ll<=NUMSCA+NUMVEC*3;ll++){
					for(modeloop=0;modeloop<NUMMODES;modeloop++){
						fprintf(mode_file," %17s%02d%02d","lm=",ll,modeloop);
					}
				}
				fprintf(mode_file,"\n");
      }
      else{//if appendold==1

				// mode file append
				if(myid<=0){
					fprintf(logfull_file,"Start setup of mode file append\n");
					fflush(logfull_file);
				}
	
				rewind(mode_file); // go to start
				// check on version info
				while(fgetc(mode_file)!='\n'); // skip comment line	      
				fscanf(mode_file,"%d %d\n",&dumi[0],&dumi[1]);
				if((dumi[0]!=MODEVER)||(dumi[1]!=MODETYPE) ){
					fprintf(fail_file,"Expected modever/modetype: %d %d got %d %d\n",MODEVER,MODETYPE,dumi[0],dumi[1]);
					myexit(6);
				}

				while(fgetc(mode_file)!='\n'); // skip comment line	
				gotit=0;
				while( (!feof(mode_file))&&(gotit==0) ){
	  
					fscanf(mode_file,INPUT4,&tcheck);
					//fprintf(stderr,"%21.15g %21.15g\n",tcheck,t);
					if(fabs(tcheck-timereenter)<0.5*DTmode){
						gotit=1;
					}// if good time (no need to continue data, just wanted time to start)
					else{
						while((fgetc(mode_file)!='\n')&&(!feof(mode_file))); // skip this bad line
					}
					fpos0=ftell(mode_file); // position to continue writting at if successful get
				}// while finding good time or hitting end of file
				if(gotit==0){
					fprintf(fail_file,"Never found right time in mode file when appending\n");
					myexit(1);
				}
				else{		
					sprintf(filename,"%s0_mode%s",DATADIR,DAT2EXT);
					sprintf(filenametemp,"%s0_mode%s.temp",DATADIR,DAT2EXT);
					sprintf(filenameback,"%s0_mode%s.back",DATADIR,DAT2EXT);
	  
					if((mode_file_temp = fopen(filenametemp,"wt"))==NULL) {
						fprintf(fail_file,"error opening temp mode output file %s\n",filenametemp) ;
						myexit(1) ;
					}
					else{
						rewind(mode_file);
						while(ftell(mode_file)<fpos0+1){
							fputc(fgetc(mode_file),mode_file_temp);
						}
						fclose(mode_file_temp);
						fclose(mode_file);
						rename(filename,filenameback); // move old to backup location
						rename(filenametemp,filename); // move new to old(normal)
						// reopen ener_file
						if((mode_file = fopen(filename,"at"))==NULL) {
							fprintf(fail_file,"error opening mode output file %s\n",filename) ;
							myexit(1) ;
						}
						if(myid<=0){
							fprintf(logfull_file,"End setup of mode file append\n");
							fflush(logfull_file);
						}
					}
				}// end else if gotit==1 on modefile
	
      }// end else if appendold==1
    }// end if myid<=0
  }// end if firsttime here


  // do always when called
  for(ll=1;ll<=NUMSCA+NUMVEC*3;ll++){
    for(modeloop=0;modeloop<NUMMODES;modeloop++){
      mode[ll][modeloop]=0;
      modei1[ll][modeloop]=0;
      modei2[ll][modeloop]=0;
    }
  }
  
  
#if(((COORD==3)&&(N3>1))) // otherwise not relevant and just slows things down
  // only applies to LOOP/BOUNDTYPE>1 in current format of looping, could redo, but not needed right now
  for(ll=1;ll<=NUMSCA+NUMVEC*3;ll++){
    which=(int)((ll-NUMSCA-1)/3)+1;
    comp=(int)(ll-NUMSCA-1)%3+1;
    // do volume mode loop
    LOOPINT2 LOOPINT1{
      for(modeloop=0;modeloop<NUMMODES;modeloop++){
				modei1[ll][modeloop]=0;
				modei2[ll][modeloop]=0;
      }
      LOOPINT3{ // phi integral
				if(ll<=NUMSCA){
					ftemp=s[ll][k][j][i];
				}
				else if(ll<=NUMSCA+3){
					if(comp==1){
						ftemp1=z2e_1(s[1],k,j,i);
					}
					else if(comp==2){
#if(COORD==3)
						ftemp1=z2e_2(s[1],k,j,i)*G2(2,i);
#elif(COORD==1)
						ftemp1=z2e_2(s[1],k,j,i)*sqrt(x[2][1][i]*x[2][1][i]+x[2][2][j]*x[2][2][j]+x[2][3][k]*x[2][3][k]); //(r)
#endif
					}
					else if(comp==3){
#if(COORD==3)
						ftemp1=z2e_3(s[1],k,j,i)*G3(2,i)*G4(2,j);
#elif(COORD==1)
						ftemp1=z2e_3(s[1],k,j,i)*sqrt(x[2][1][i]*x[2][1][i]+x[2][2][j]*x[2][2][j]); // r*sin(theta)
#endif
					}
					ftemp=ftemp1*v[which][comp][k][j][i]; // ang momentum (so density weighted)
				}
				else if(ll<=NUMSCA+2*3){
					ftemp=v[which][comp][k][j][i]; // weighting for B?
				}
#if(COORD==3)
				ftemp1=ftemp*dx[1][3][k]; // no need for metric here since not phi-dep, just do below
				ftemp2=x[2][3][k];
#elif(COORD==1)
				ftemp1=ftemp*dx[1][3][k]; // dz really now
				// cart -> phi in spcoord
				if(x[2][1][i]>0.0){
					if(x[1][2][j]>0.0) ftemp2=atan(x[1][2][j]/x[2][1][i]);
					else if(x[1][2][j]<0.0) ftemp2=2.0*M_PI+atan(x[1][2][j]/x[2][1][i]);
					else ftemp2=0.0;
				}
				else if(x[2][1][i]<0.0){
					ftemp2=M_PI+atan(x[1][2][j]/x[2][1][i]);
				}
				else{
					if(x[1][2][j]>0) ftemp2=M_PI*0.5;
					else if(x[1][2][j]<0) ftemp2=M_PI*1.5;
					else ftemp2=0.0;
				}
				// done getting phi coord
#endif
				for(modeloop=0;modeloop<NUMMODES;modeloop++){
					ftemp3=ftemp2*(SFTYPE)(modeloop);
					modei1[ll][modeloop]+=cos(ftemp3)*ftemp1;
					modei2[ll][modeloop]+=sin(ftemp3)*ftemp1;
					//	    fprintf(stdout,"%02d %02d %15.10g %15.10g %15.10g %15.10g\n",ll,modeloop,modei1[ll][modeloop],modei2[ll][modeloop],ftemp3,ftemp1); fflush(stdout);
				}
      }
      /*
				for(modeloop=0;modeloop<NUMMODES;modeloop++){
				fprintf(stdout,"%02d %02d %15.10g %15.10g\n",ll,modeloop,modei1[ll][modeloop],modei2[ll][modeloop]); fflush(stdout);
				}
      */
      // assign true mode so far ( can split this since r/theta/phi are independent coordinates and no phi dependent metric components, and all metrics>0)
      for(modeloop=0;modeloop<NUMMODES;modeloop++){
				mode[ll][modeloop]+=sqrt(modei1[ll][modeloop]*modei1[ll][modeloop]+modei2[ll][modeloop]*modei2[ll][modeloop])
#if(COORD==3)
					*dx[1][1][i]*dx[1][2][j]*G2(2,i)*G3(2,i)*G4(2,j)
#elif(COORD==1)
					*dx[1][1][i]*dx[1][2][j]
#endif
					;
      }
    }// end of volume mode loop
    
    // finish off with a normalization
    for(modeloop=0;modeloop<NUMMODES;modeloop++){
      mode[ll][modeloop]/=
#if(COORD==3)
				(((pow(x1out,3.0)-pow(x1in,3.0))*(-cos(x2out)+cos(x2in))*(x3out-x3in))) // mode=1.0 if fun=1 and m=0
#elif(COORD==1)
				((x1out-x1in)*(x2out-x2in)*(x3out-x3in))
#endif
				;
    }
  }
#endif
  
  
  
  
#if(COORD==1)
  for(ll=1;ll<=NUMSCA+NUMVEC*3;ll++){
    for(modeloop=0;modeloop<NUMMODES;modeloop++){
      for(i=0;i<MAXSHELLS;i++){
				for(j=0;j<MAXSHELLS;j++){
					modei1cart[i][j][ll][modeloop]=0;
					modei2cart[i][j][ll][modeloop]=0;
				}
      }
    }
  }
  // generate map for phi
  // this tells mode if in same phi-loop or not
  for(ll=1;ll<=NUMSCA+NUMVEC*3;ll++){
    which=(int)((ll-NUMSCA-1)/3)+1;
    comp=(int)(ll-NUMSCA-1)%3+1;
    // do volume mode loop
    TOTALVOLUME=0;
    LOOPINT{
      if(ll<=NUMSCA){
				radius=cart2spc(1,x[2][1][i],x[2][2][j],x[2][3][k]);
				theta=cart2spc(2,x[2][1][i],x[2][2][j],x[2][3][k]);
				phi=cart2spc(3,x[2][1][i],x[2][2][j],x[2][3][k]);
	
				dV=dx[1][1][i]*dx[1][2][j]*dx[1][3][k];
				Jac=fabs(1.0/(radius*radius*sin(theta)));
				if(Jac>1E10) Jac=0.0;
				Jac=1.0;
				Fun=s[ll][k][j][i]*dV*Jac;
				TOTALVOLUME+=dV;
      }
      else{
				if(comp==1){
					radius=cart2spc(1,x[1][1][i],x[2][2][j],x[2][3][k]);
					theta=cart2spc(2,x[1][1][i],x[2][2][j],x[2][3][k]);
					phi=cart2spc(3,x[1][1][i],x[2][2][j],x[2][3][k]);

					dV=dx[2][1][i]*dx[1][2][j]*dx[1][3][k];
					Jac=fabs(1.0/(radius*radius*sin(theta)));
					if(Jac>1E10) Jac=0.0;
					Jac=1.0;
					TOTALVOLUME+=dV;
				}
				else if(comp==2){
					radius=cart2spc(1,x[2][1][i],x[1][2][j],x[2][3][k]);
					theta=cart2spc(2,x[2][1][i],x[1][2][j],x[2][3][k]);
					phi=cart2spc(3,x[2][1][i],x[1][2][j],x[2][3][k]);

					dV=dx[1][1][i]*dx[2][2][j]*dx[1][3][k];
					Jac=fabs(1.0/(radius*radius*sin(theta)));
					if(Jac>1E10) Jac=0.0;
					Jac=1.0;
					TOTALVOLUME+=dV;
				}
				else if(comp==3){
					radius=cart2spc(1,x[2][1][i],x[2][2][j],x[1][3][k]);
					theta=cart2spc(2,x[2][1][i],x[2][2][j],x[1][3][k]);
					phi=cart2spc(3,x[2][1][i],x[2][2][j],x[1][3][k]);
	  
					dV=dx[1][1][i]*dx[1][2][j]*dx[2][3][k];
					Jac=fabs(1.0/(radius*radius*sin(theta)));
					if(Jac>1E10) Jac=0.0;
					Jac=1.0;
					TOTALVOLUME+=dV;
				}
				if(ll<=NUMSCA+3){ // density scale
					if(comp==1) Fun=z2e_1(s[1],k,j,i)*Jac*dV*v[which][comp][k][j][i];
					else if(comp==2) Fun=z2e_2(s[1],k,j,i)*Jac*dV*v[which][comp][k][j][i];
					else if(comp==3) Fun=z2e_3(s[1],k,j,i)*Jac*dV*v[which][comp][k][j][i];
				}
				else if(ll<=NUMSCA+2*3){
					Fun=Jac*dV*v[which][comp][k][j][i]; // weight B?
				}
      }
      // do each shell, like coord==3 case

#if(MODEMETHOD==1)
      for(modeloop=0;modeloop<NUMMODES;modeloop++){
        ftemp3=phi*(SFTYPE)(modeloop);
        modei1[ll][modeloop]+=cos(ftemp3)*Fun;
        modei2[ll][modeloop]+=sin(ftemp3)*Fun;
      }
#elif(MODEMETHOD==2)
      // Rinner and Router must be defined
      for(modeloop=0;modeloop<NUMMODES;modeloop++){
				ftemp3=phi*(SFTYPE)(modeloop);
				ri=(int)(MAXSHELLS*.999/(Router-Rinner)*(radius-Rinner));
				thetai=(int)(MAXSHELLS*0.999/M_PI*theta);
				if(ri<0) ri=0;
				if(ri>MAXSHELLS-1) ri=MAXSHELLS-1;
				if(thetai<0) thetai=0;
				if(thetai>MAXSHELLS-1) thetai=MAXSHELLS-1;
				modei1cart[ri][thetai][ll][modeloop]+=cos(ftemp3)*Fun;
				modei2cart[ri][thetai][ll][modeloop]+=sin(ftemp3)*Fun;
      }
#endif
    }// end volume loop
    
#if(MODEMETHOD==1)
    for(modeloop=0;modeloop<NUMMODES;modeloop++){
      // mode[ll][modeloop]=sqrt(modei1[ll][modeloop]*modei1[ll][modeloop]+modei2[ll][modeloop]*modei2[ll][modeloop])/((x1out-x1in)*(x2out-x2in)*(x3out-x3in));
      mode[ll][modeloop]=sqrt(modei1[ll][modeloop]*modei1[ll][modeloop]+modei2[ll][modeloop]*modei2[ll][modeloop]);
    }
#elif(MODEMETHOD==2)
    for(modeloop=0;modeloop<NUMMODES;modeloop++){
      for(i=0;i<MAXSHELLS;i++){
				for(j=0;j<MAXSHELLS;j++){
					// sum up shells
					mode[ll][modeloop]+=sqrt(modei1cart[i][j][ll][modeloop]*modei1cart[i][j][ll][modeloop]+modei2cart[i][j][ll][modeloop]*modei2cart[i][j][ll][modeloop]);
				}
      }
      // divide out spacial domain
      //      mode[ll][modeloop]=M_PI*(pow(Router,3.0)-pow(Rinner,3.0))/((x1out-x1in)*(x2out-x2in)*(x3out-x3in));
      //      mode[ll][modeloop]/=((x1out-x1in)*(x2out-x2in)*(x3out-x3in));
      mode[ll][modeloop]/=1.0;
    }
#endif
  } // end over variables

#endif
  
  // now have mode for other functions
  
  // add up volume integrals over all cpus
  if(numprocs>1){
#if(USEMPI)
    for(ll=1;ll<=NUMSCA+NUMVEC*3;ll++){
      for(modeloop=0;modeloop<NUMMODES;modeloop++){
				MPI_Reduce(&(mode[ll][modeloop]),&(mode_full[ll][modeloop]), 1, MPI_SFTYPE, MPI_SUM,  0,MPI_COMM_WORLD);
      }
    }    
#endif
  }
  else{
    for(ll=1;ll<=NUMSCA+NUMVEC*3;ll++){
      for(modeloop=0;modeloop<NUMMODES;modeloop++){
				mode_full[ll][modeloop]=mode[ll][modeloop];
      }
    }
  }
  
  
  // now write results using cpu=0 (which received the totals)
  if(myid<=0){
    // now output mode data into seperate file
    fprintf(mode_file," %21.15g",t) ;
    for(ll=1;ll<=NUMSCA+NUMVEC*3;ll++){
      for(modeloop=0;modeloop<NUMMODES;modeloop++){
				fprintf(mode_file," %21.15g"
								,mode_full[ll][modeloop]
					);
      }
    }
    fprintf(mode_file,"\n");
    fflush(mode_file) ;
    
  }// end if cpu to write to file
  
  // cleanup if final call
  if(call_code==2){
    if(myid<=0){
      fclose(mode_file);
    }// end over write cpu
  }// end if final call
  
  firsttime=0;
  told=t;
}




void diag_dumps(int call_code,FTYPE tdump,FTYPE tfldump,FTYPE tpdump,FTYPE timage)
	// call_code
	// -1: do every call to diag(usually each time step as called by main.c)
	// 0: first call setup/dump
	// 1: normal diag dump of variables each DTl (DT for logging)
	// 2: final diag call after done
     
{
  static SFTYPE told[4]={-1000.00111,-1000.00111,-1000.00111,-1000.00111};
  fpos_t fposgod;
  SFTYPE tempf,tempf2;
  SFTYPE tempfns;
  SFTYPE masstempc,masstemp1,masstemp2,masstemp3;
  char dfnam[MAXFILENAME];
  char dfheader[MAXFILENAME];
  char dfnamtemp[MAXFILENAME];
  char dfnamback[MAXFILENAME];

  FILE *dump_file ;
  FILE *image_file ;
  int i,j,k,l,m,ii,ll ;
  int loopdido;
  int dumi[10];
  SFTYPE tcheck; // used to check time for loss append

  int which, comp;
  SFTYPE realloss; // need for enthalpy and sm
  SFTYPE ftemp,ftemp1,ftemp2,ftemp3;
  
  static int dump_cnt,pdump_cnt,im_cnt ; // used to enumerate files
  static int adump_cnt,npdump_cnt,floordump_cnt,fldump_cnt;
  static int timeless=0;
#if(USEMPI)
  static MPI_Request request[6]; 
#endif
  int typenum;
  char temps[MAXFILENAME];
  char filename[MAXFILENAME];
  char filenametemp[MAXFILENAME];
  char filenameback[MAXFILENAME];
  long fpos0;
  int gotit;

  //  fprintf(stdout,"callcode: %d timage: %15.10g\n",call_code,timage); fflush(stdout);

  //  if(told==t) return; // to avoid duplicate data output


  // setup dump stuff
  if(call_code==0){
    // these are dump#'s used next
    pdump_cnt  = pdump_start ;
    dump_cnt   = dump_start ;
    fldump_cnt   = fldump_start ;
    npdump_cnt = npdump_start ;
    adump_cnt  = adump_start;
    floordump_cnt  = floor_start;
    im_cnt     = image_start ;

  }

  // dump
  if(call_code>=0){

    if((told[0]!=t)&&(PDUMPFLAG&&( ((dt<DTLOWDUMP)&&(CHECKDTLOW==1))||(t >= tpdump) || (call_code == 2) ) ) ) {
      dump(NULL,pdump_cnt,PDTYPE,0);// primitive variable dumps
      pdump_cnt++ ;
      
      
      fprintf(log_file,"proc: %2d pdump: %5d, cc: %5d ...\n",myid,pdump_cnt,call_code);
      fflush(log_file);
      
      if(myid<=0){
				sprintf(dfnam,"%s0_numpdumps%s",DATADIR,DAT2EXT);
				if((dump_file = fopen(dfnam,"w"))==NULL) {
					fprintf(fail_file,"error opening dump output file %s\n",dfnam) ;
					myexit(1) ;
				}
				fprintf(dump_file,"#%30s\n","Number of pdumps");
				fprintf(dump_file,"%30d\n",pdump_cnt);
				fclose(dump_file);
      }

      fprintf(log_file,"proc: %2d pdump done: %d t=%15.10g.\n",myid,pdump_cnt,tpdump);
      fflush(log_file);
      told[0]=t;

    }
    // GODMARK
    if((told[1]!=t)&&(((DUMPFLAG)||(NPDUMPFLAG)||(FLOORDUMPFLAG)||(ADUMPFLAG&&analoutput>0))&&( ((dt<DTLOWDUMP)&&(CHECKDTLOW==1))||(t >= tdump) || (call_code == 2) ) )  ){
      
      if((DUMPFLAG)||(NPDUMPFLAG)||(FLOORDUMPFLAG)||(ADUMPFLAG&&analoutput>0)){
				fprintf(log_file,"proc: %2d dump: %5d, cc: %5d ...\n",myid,dump_cnt,call_code);
				fflush(log_file);
      }
      if(DUMPFLAG==1){
				dump(NULL,dump_cnt,DTYPE,0);// normal variable dumps
				dump_cnt++ ;
      }
      if(NPDUMPFLAG==1){
				dump(NULL,npdump_cnt,NPTYPE,0);// non-primitive variable dumps(those things complicated to compute)
				npdump_cnt++ ;
      }
      if(FLOORDUMPFLAG==1){
				dump(NULL,floordump_cnt,FLTYPE,0);// floor dump
				floordump_cnt++ ;
      }
      if(analoutput>0){
				if( (ADUMPFLAG==1)||( (call_code==0)&&(ADUMPFLAG==-1))){	  
					dump(NULL,adump_cnt,ADTYPE,0) ;
					adump_cnt++;
				}
      }
      if((DUMPFLAG)||(NPDUMPFLAG)||(FLOORDUMPFLAG)||(ADUMPFLAG&&analoutput>0)){
				fprintf(log_file,"proc: %2d dump done: %d t=%15.10g.\n",myid,dump_cnt,tdump);
				fflush(log_file);
      }
      // now write number of dumps out to file for pp control
      if(myid<=0){
				// normal dump
				sprintf(dfnam,"%s0_numdumps%s",DATADIR,DAT2EXT);
				if((dump_file = fopen(dfnam,"w"))==NULL) {
					fprintf(fail_file,"error opening dump output file %s\n",dfnam) ;
					myexit(1) ;
				}
				fprintf(dump_file,"#%30s\n","Number of dumps");
				fprintf(dump_file,"%30d\n",dump_cnt);
				fclose(dump_file);

				// np dumps
				sprintf(dfnam,"%s0_numnpdumps%s",DATADIR,DAT2EXT);
				if((dump_file = fopen(dfnam,"w"))==NULL) {
					fprintf(fail_file,"error opening dump output file %s\n",dfnam) ;
					myexit(1) ;
				}
				fprintf(dump_file,"#%30s\n","Number of npdumps");
				fprintf(dump_file,"%30d\n",npdump_cnt);
				fclose(dump_file);

				// floor dumps
				sprintf(dfnam,"%s0_numfloordumps%s",DATADIR,DAT2EXT);
				if((dump_file = fopen(dfnam,"w"))==NULL) {
					fprintf(fail_file,"error opening dump output file %s\n",dfnam) ;
					myexit(1) ;
				}
				fprintf(dump_file,"#%30s\n","Number of floordumps");
				fprintf(dump_file,"%30d\n",floordump_cnt);
				fclose(dump_file);

				// adump
				sprintf(dfnam,"%s0_numadumps%s",DATADIR,DAT2EXT);
				if((dump_file = fopen(dfnam,"w"))==NULL) {
					fprintf(fail_file,"error opening dump output file %s\n",dfnam) ;
					myexit(1) ;
				}
				fprintf(dump_file,"#%30s\n","Number of adumps");
				fprintf(dump_file,"%30d\n",adump_cnt);
				fclose(dump_file);

      }
      told[1]=t;
    }

    if((told[2]!=t)&&(FLDUMPFLAG&&((t >= tfldump) || (call_code == 2) ) ) ) {

      fprintf(log_file,"proc: %2d fldump: %5d, cc: %5d ...\n",myid,fldump_cnt,call_code);
      fflush(log_file);
      
      dump(NULL,fldump_cnt,FLINETYPE,0);
      fldump_cnt++;

      fprintf(log_file,"proc: %2d fldump done: %d t=%15.10g.\n",myid,fldump_cnt,tfldump);
      fflush(log_file);

      // now write number of dumps out to file for pp control
      if(myid<=0){
				// normal dump
				sprintf(dfnam,"%s0_numfldumps%s",DATADIR,DAT2EXT);
				if((dump_file = fopen(dfnam,"w"))==NULL) {
					fprintf(fail_file,"error opening fldump output file %s\n",dfnam) ;
					myexit(1) ;
				}
				fprintf(dump_file,"#%30s\n","Number of fldumps");
				fprintf(dump_file,"%30d\n",fldump_cnt);
				fclose(dump_file);
      }
      told[2]=t;
    }

    /* make image of variables at regular intervals */
    if((told[3]!=t)&&(IMAGEFLAG&&( (t >= timage)|| (call_code == 2) ) ) ) {
      fprintf(log_file,"proc: %2d image: %5d, cc: %5d ...\n",myid,im_cnt,call_code);      fflush(log_file);
      image(im_cnt,-1,-1,call_code,0); // set which images to produce here
      fprintf(log_file,"proc: %2d image done: %d t=%15.10g.\n",myid,im_cnt,timage);      fflush(log_file);
      im_cnt++ ;
      
      if(myid<=0){
				sprintf(temps,"%s%s",DATADIR,IMAGEDIR);	
	
				sprintf(dfnam,"%s0_numimages%s",temps,DAT2EXT);
				if((image_file = fopen(dfnam,"w"))==NULL) {
					fprintf(fail_file,"error opening dump output file %s\n",dfnam) ;
					myexit(1) ;
				}
				fprintf(image_file,"#%30s\n","Number of images");
				fprintf(image_file,"%30d\n",im_cnt);
				fclose(image_file);
      }
      told[3]=t;
    }

  }//end if call_code>=0
  
}





// if change any headers, change init.c's input stuff!!

// which sign used to determine sample type, and so size/sample/zonec output
// currently used by all except rwhich==11 (avg1d)
void dump_header(FILE *fp,int which, int realsampled,int nogridchoice)
{
  int realzone;
  int realn1,realn2,realn3;
  if(nogridchoice==1){
    realzone=0;
  }
  else{
    realzone=DUMPSM;
  }
  if(mpicombine==1){
    realn1=totalsize[1];
    realn2=totalsize[2];
    realn3=totalsize[3];
  }
  else if(which==AVG1DTYPE){
    if( (realsampled==0)||(nogridchoice==1) ){ // no interp
      if(FULLOUTPUT==0){
				realn1=totalsize[1];
				realn2=totalsize[2];
				realn3=totalsize[3];
      }
      else if(FULLOUTPUT==1){
				realn1=totalsize[1]+N1BND*N1NOT1;
				realn2=totalsize[2]+N2BND*N2NOT1;
				realn3=totalsize[3]+N3BND*N3NOT1;
      }
      else if(FULLOUTPUT==2){
				realn1=totalsize[1]+2*N1BND*N1NOT1;
				realn2=totalsize[2]+2*N2BND*N2NOT1;
				realn3=totalsize[3]+2*N3BND*N3NOT1;
      }
    }
    else{
      realn1=itotalsize[1];
      realn2=itotalsize[2];
      realn3=itotalsize[3];
    }
  }
  else{
    if( (realsampled==0)||(nogridchoice==1) ){ // no interp
      if(FULLOUTPUT==0){
				realn1=N1;
				realn2=N2;
				realn3=N3;
      }
      else if(FULLOUTPUT==1){
				realn1=N1+N1BND*N1NOT1;
				realn2=N2+N2BND*N2NOT1;
				realn3=N3+N3BND*N3NOT1;
      }
      else if(FULLOUTPUT==2){
				realn1=N1+2*N1BND*N1NOT1;
				realn2=N2+2*N2BND*N2NOT1;
				realn3=N3+2*N3BND*N3NOT1;
      }
    }
    else{
      realn1=DUMN1;
      realn2=DUMN2;
      realn3=DUMN3;
    }
  }

  // version header
  if(which==DTYPE){
    fprintf(fp,"#%10s %10s\n%10d %10d\n","DVER","DTYPE",DVER,DTYPE);
  }
  else if(which==ADTYPE){
    fprintf(fp,"#%10s %10s\n%10d %10d\n","ADVER","ADTYPE",ADVER,ADTYPE);
  }
  else if(which==PDTYPE){
    fprintf(fp,"#%10s %10s\n%10d %10d\n","PDVER","PDTYPE",PDVER,PDTYPE);
  }
  else if(which==FLTYPE){
    fprintf(fp,"#%10s %10s\n%10d %10d\n","FLVER","FLTYPE",FLVER,FLTYPE);
  }
  else if(which==NPTYPE){
    fprintf(fp,"#%10s %10s\n%10d %10d\n","NPVER","NPTYPE",NPVER,NPTYPE);
  }
  else if(which==AVG2DTYPE){
    fprintf(fp,"#%10s %10s\n%10d %10d\n","AVG2DVER","AVG2DTYPE",AVG2DVER,AVG2DTYPE);
  }
  else if(which==AVG1DTYPE){
    fprintf(fp,"#%10s %10s\n%10d %10d\n","AVG1DVER","AVG1DTYPE",AVG1DVER,AVG1DTYPE);
  }
  else if(which==CALCTYPE){
    fprintf(fp,"#%10s %10s\n%10d %10d\n","CALCVER","CALCTYPE",CALCVER,CALCTYPE);
  }
  else if(which==FLINETYPE){
    fprintf(fp,"#%10s %10s\n%10d %10d\n","FLINEVER","FLINETYPE",FLINEVER,FLINETYPE);
  }
  fprintf(fp,"#%21s %6s %6s\n","t","SAMPLE","ZONEC") ;
  fprintf(fp," %21.15g %6d %6d\n",t,realsampled,realzone) ;
  fprintf(fp,"#%6s %6s %6s\n","N1","N2","N3") ;
  fprintf(fp," %6d %6d %6d\n",realn1,realn2,realn3);

  if(which==AVG2DTYPE){ // extra needed info
    fprintf(fp,"# Averagefromto %21.15g %21.15g\n",tavgstart,tavgfinal);
    fprintf(fp,"# SAMPLENUM %d\n",avgcount);
  }
  if(which==AVG1DTYPE){ // extra needed info
    fprintf(fp,"# Averagefromto %21.15g %21.15g\n",tavgstart,tavgfinal);
    fprintf(fp,"# SAMPLETCNT %d SAMPLE1DNUM %d %d\n",avgcount,num1d_31_full,num1d_32_full);
  }
}



void dump_header2(FILE *fp, int which)
{

  if( (which==DTYPE)||(which==PDTYPE)||(which==ADTYPE)||(which==FLTYPE) ){
    fprintf(fp,"#%21s %21s %21s ", "rho","u","pot");
  }
  if(which==FLTYPE){ // add in ke for floor
    fprintf(fp,"%21s ", "ke");
  }
  if( (which==DTYPE)||(which==PDTYPE)||(which==ADTYPE)||(which==FLTYPE) ){
    fprintf(fp,"%21s %21s %21s %21s %21s %21s\n", "vx1","vx2","vx3","bx1","bx2","bx3") ;
  }
  if(which==NPTYPE){
    fprintf(fp,"#%21s %21s %21s %21s %21s %21s %21s\n", "sigma11","sigma12","sigma13","sigma22","sigma23","sigma33","nuvisc") ;
  }
  if(which==AVG2DTYPE){
    fprintf(fp,"#%21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s\n","2D_rho","2D_en","2D_Be","2D_cs2","2D_Exp[S]","2D_vx1","2D_vx2","2D_vx3","2D_sig11","2D_sig12","2D_sig13","2D_sig22","2D_sig23","2D_sig33","2D_nuvisc");
  }
  if(which==AVG1DTYPE){
    fprintf(fp,"#%21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s\n","1D1_rho","1D1_en","1D1_Be","1D1_cs2","1D1_S","1D1_vx1","1D1_vx2","1D1_vx3","1D1_sig11","1D1_sig12","1D1_sig13","1D1_sig22","1D1_sig23","1D1_sig33","1D1_nuvisc","1D2_rho","1D2_en","1D2_Be","1D2_cs2","1D2_S","1D2_vx1","1D2_vx2","1D2_vx3","1D2_sig11","1D2_sig12","1D2_sig13","1D2_sig22","1D2_sig23","1D2_sig33","1D2_nuvisc");
  }      
  if(which==CALCTYPE){
    fprintf(fp,"#%21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s\n", "LRx","LRy","LRz","LVx","LVy","LVz","ERx","ERy","ERz","EVx","EVy","EVz");
  }
  if(which==FLINETYPE){
    fprintf(fp,"#%21s\n","Ab-x3");
  }
}


//wrapper to control region real data is outputted
void dump1data(int which, int i, int j, int k, FILE* fileptr, FTYPE var)
{

  if(!inside_accretor(which,i,j,k)){
    fprintf(fileptr,"%21.15g ",var);
  }
  else{
    fprintf(fileptr,"%21.15g ",0.0);
  }
}

//wrapper to control region real data is outputted (same order as fwrite)
void dump1databin(int which, int i, int j, int k, FTYPE *var, int datasize, int numdata, FILE* fileptr)
{
  FTYPE zerovar=0.0;
  
  if(!inside_accretor(which, i,j,k)){
    fwrite(var,sizeof(FTYPE),1,fileptr);
  }
  else{
    fwrite(&zerovar,sizeof(FTYPE),1,fileptr);
  }
}



//wrapper to control region real data is outputted (for writebuf)
void writebuf1data(int which, int i, int j, int k, FTYPE var, FTYPE *writebuf)
{
  if(!inside_accretor(which,i,j,k)){
    *writebuf=var;
  }
  else{
    *writebuf=0.0;
  }
}

void dumpcheckfprintf(int which, int i, int j, int k, FILE* fileptr, char *format, ...)
{
  va_list arglist;

  va_start (arglist, format);

  if(fileptr==NULL){
    fprintf(stderr,"tried to print to null file pointer: %s\n",format);
    fflush(stderr);
  }
  else{
    if(!inside_accretor(which, i,j,k)){
      vfprintf (fileptr, format, arglist);
      fflush(fileptr);
    }
    else{
      // for now not sure how to restrict values in arglist
      vfprintf (fileptr, format, arglist);
      fflush(fileptr);      
    }
  }
  va_end (arglist);
}



// fp: file pointer to write to if applicable
// dump_cnt: file # if applicable
// which: see types in global.h
// outtype: zoom or not for interp
void dump(FILE * fp, int dump_cnt, int which,int outtype)
{
  static FTYPE (*iq)[INTN2][INTN1];
  static FTYPE (*viq)[3][INTN2][INTN1];
  static FTYPE (*iq1)[INTN1];
  static FTYPE (*iq2)[INTN1];
  static FTYPE (*iq3)[INTN1];

  static FTYPE (*iqT00)[INTN1];
  static FTYPE (*iqT01)[INTN1];
  static FTYPE (*iqT02)[INTN1];
  static FTYPE (*iqT11)[INTN1];
  static FTYPE (*iqT12)[INTN1];
  static FTYPE (*iqT22)[INTN1];

  int realn1,realn2;
  int computesigma;
  int i,j,k,l,m,p,q ;
  int kk,jj,ii;
  FTYPE Bxa,Bya,Bza ;
  FTYPE tempf;

  FTYPE pos3[3+1];
  FTYPE ftemp[NUMVEC+1][3+1];
  FTYPE ftempT[6+1];

  FTYPE ftempC[12+1];

  FTYPE (*ax3)[N2M][N1M];
  FTYPE (*ax3a)[N2M][N1M];
  FTYPE (*iax3)[N2M][N1M];

  FTYPE (*scain)[N3M][N2M][N1M];
  FTYPE (*vecin)[3][N3M][N2M][N1M];
  char dheader[MAXFILENAME];
  char dfheader[MAXFILENAME];
  char temps[MAXFILENAME];
  char dfnam[MAXFILENAME];
  FTYPE odef[NUMAVG2_3+1];

  static FTYPE (*isigma)[3][N3M][N2M][N1M]; // -2*rho*nu*e_{ij}= sigma_{ij}
  static FTYPE (*irost)[3][N3M][N2M][N1M]; // e_{ij}
  static FTYPE (*irostnu)[3][N3M][N2M][N1M]; // e_{ij}*nu
  static FTYPE (*nurho_real)[N2M][N1M];
  static FTYPE (*ippcalcs)[N3M][N2M][N1M];
  static FTYPE (*differva)[N3M][N2M][N1M];
  static FTYPE (*differvb)[N3M][N2M][N1M];
  static FTYPE (*inu_real)[N2M][N1M]; // nu_real alias
  static FTYPE (*iavg2_3)[N2M][N1M] ;
 
  static FTYPE (*delv)[N2M][N1M]; // deldotv
  int realsampled;
  FTYPE vxa,vya,vza,sig11a,sig12a,sig21a,sig13a,sig31a,sig22a,sig23a,sig32a,sig33a,vtot2a,enthalpy,pot;
  static int firsttime=1;
  int fpnull;
  int nogridchoice;
  int it;
  int itemp;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
#if(USEMPI&&USEROMIO)
  MPI_Datatype newtype;
  MPI_File fh;
  MPI_Status status;
  MPI_Request request;
#endif
#if(USEMPI)
  FTYPE *jonio;
#endif
  FTYPE totalnorm,recvnorm,sendnorm;
  int ndims, array_of_gsizes[4], array_of_distribs[4];
  int order, len;
  int array_of_dargs[4], array_of_psizes[4];
  int bufcount, array_size;
  FTYPE *writebuf,*tempbuf;
  int numcolumns;
  int dumi;
  FTYPE dumf;
  SFTYPE dumlf;
  char truemyidtxt[100];
  int nextbuf;
  FILE * headerptr;
  int dumpsmdump;
  

  // grab work space for interp(note workspaced used below too)
  iq=workiq;
  viq=workviq;

  // defaults
  computesigma=0;
  if(outtype>=0){
    realsampled=SAMPLED;
  }
  else{
    realsampled=0;
    outtype=0;
  }
  nogridchoice=0;
  realn1=N1;
  realn2=N2;

  if(which==NPCOMPUTETYPE){ // otherwise same as NPTYPE
    computesigma=1;
    realsampled=0;
    which=NPTYPE;
  }
  if(which==PDTYPE){
		nogridchoice=1;
		dumpsmdump=0;
  }
  else dumpsmdump=DUMPSM;

  if(fp==NULL){
    fpnull=1;

    if(which==DTYPE){
      if(realsampled>0) strcpy(dfheader,"idump");
      else strcpy(dfheader,"dump");
      if(analoutput==13) numcolumns=7;
      else numcolumns=9;
    }
    else if(which==PDTYPE){
      if(realsampled>0) strcpy(dfheader,"ipdump");
      else strcpy(dfheader,"pdump");
      if(analoutput==13) numcolumns=7;
      else numcolumns=9;
    }
    else if(which==ADTYPE){
      if(realsampled>0) strcpy(dfheader,"aidump");
      else strcpy(dfheader,"adump");
      if(analoutput==13) numcolumns=7;
      else numcolumns=9;
    }
    else if(which==FLTYPE){
      if(realsampled>0) strcpy(dfheader,"ifloor");
      else strcpy(dfheader,"floor");
      if(analoutput==13) numcolumns=7;
      else numcolumns=11;
    }
    else if(which==NPTYPE){
      if(realsampled>0) strcpy(dfheader,"npidump"); // not inp due to sm read method
      else strcpy(dfheader,"npdump");
			if(VISCMEM){
				numcolumns=7;
			}
			else numcolumns=0;
			if(1) numcolumns+=72;
    }
    else if(which==AVG2DTYPE){
      if(realsampled>0) strcpy(dfheader,"i0_avg2d");
      else strcpy(dfheader,"0_avg2d");
      numcolumns=15;
    }
    else if(which==CALCTYPE){
      if(realsampled>0) strcpy(dfheader,"icdump"); // not inp due to sm read method
      else strcpy(dfheader,"cdump");
      numcolumns=12;
    }
    else if(which==FLINETYPE){
      if(realsampled>0) strcpy(dfheader,"ifldump"); // not inp due to sm read method
      else strcpy(dfheader,"fldump");
      numcolumns=1;
    }

    if(mpicombine==0){
      strcpy(truemyidtxt,myidtxt);
    }
    else strcpy(truemyidtxt,"");

    if(dump_cnt==-1){ // e.g. when which==AVG2DTYPE
      sprintf(dfnam,"%s%s%s%s",DATADIR,dfheader,DATEXT,truemyidtxt) ;        
    }
    else{
      sprintf(dfnam,"%s%s%04d%s%s",DATADIR,dfheader,dump_cnt,DATEXT,truemyidtxt);
    }


    if(mpicombine==1){
#if(USEMPI)
			if(USEROMIO){
#if(USEROMIO)
				//create the distributed array filetype
				ndims = 4;
				order = MPI_ORDER_C;
      
				array_of_gsizes[3] = numcolumns;
				array_of_gsizes[2] = totalsize[1];
				array_of_gsizes[1] = totalsize[2];
				array_of_gsizes[0] = totalsize[3];
      
				array_of_distribs[3] = MPI_DISTRIBUTE_BLOCK;
				array_of_distribs[2] = MPI_DISTRIBUTE_BLOCK;
				array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
				array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
      
				array_of_dargs[3] = MPI_DISTRIBUTE_DFLT_DARG;
				array_of_dargs[2] = MPI_DISTRIBUTE_DFLT_DARG;
				array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
				array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;

				array_of_psizes[3]=1;
				array_of_psizes[2]=ncpux1;
				array_of_psizes[1]=ncpux2;
				array_of_psizes[0]=ncpux3;
      
				MPI_Type_create_darray(numprocs, myid, ndims, array_of_gsizes, array_of_distribs, array_of_dargs, array_of_psizes, order, MPI_FTYPE, &newtype);
				MPI_Type_commit(&newtype);
				MPI_Type_size(newtype, &bufcount);
				bufcount = bufcount/sizeof(FTYPE);
				writebuf = (FTYPE *) malloc(bufcount * sizeof(FTYPE));
				// setup file handler

				MPI_File_open(MPI_COMM_WORLD, dfnam, MPI_MODE_CREATE | MPI_MODE_RDWR, 
											MPI_INFO_NULL, &fh);
				MPI_File_set_view(fh, 0, MPI_FTYPE, newtype, "native", MPI_INFO_NULL);
				// all that needs to be done now is fill writebuf with the data
#endif
			}
			else if(USEJONIO){
				jonio_init_fp(&fp,1,dfnam);
				if(sizeof(FTYPE)==sizeof(double)) jonio_init_mem(numcolumns,sizeof(FTYPE),0,0,&jonio,0,0,&writebuf,0,0,&tempbuf);
				else if(sizeof(FTYPE)==sizeof(float)) jonio_init_mem(numcolumns,sizeof(FTYPE),0,&jonio,0,0,&writebuf,0,0,&tempbuf,0);
				else if(sizeof(FTYPE)==sizeof(unsigned char)) jonio_init_mem(numcolumns,sizeof(FTYPE),&jonio,0,0,&writebuf,0,0,&tempbuf,0,0);
			}
#endif
    }
    else{
      if(!binaryoutput){
				if((fp=fopen(dfnam,"wt"))==NULL){
					fprintf(fail_file,"Cannot open: %s\n",dfnam);
					myexit(1);
				}
				headerptr=fp;
      }
      else{ // use binary mode
				if((fp=fopen(dfnam,"w"))==NULL){
					fprintf(fail_file,"Cannot open: %s\n",dfnam);
					myexit(1);
				}
      }
    }
  }
  else{
    fpnull=0;
    mpicombine=0; // can't send in file pointer in MPI mode when doing parallel i/o yet, so force this
  }

  if( ((myid==0)&&(mpicombine==1))||((mpicombine==0)&&(binaryoutput==1)) ){
    // header
    strcat(dfnam,".head");
    if((headerptr=fopen(dfnam,"wt"))==NULL){
      fprintf(fail_file,"Can't open %s for writting\n",dfnam);
      myexit(1);
    }
  }
  
  if( ((myid==0)&&(mpicombine==1))||(mpicombine==0) ){
    dump_header(headerptr,which,realsampled,nogridchoice);
    dump_header2(headerptr,which);
  }

  if( (realsampled==0)||(nogridchoice==1) ){ // no interp
    realn1=N1;
    realn2=N2;
  }
  else{// sizes for input data
    if(FULLOUTPUT==0){
      realn1=N1;
      realn2=N2;
    }
    else if(FULLOUTPUT==1){
      // assumes 2 boundary zones
      realn1=OUTFULL1;
      realn2=OUTFULL2;
    }
    else if(FULLOUTPUT==2){
      realn1=N1M;
      realn2=N2M;
    }
  }

  if( (which==DTYPE)||(which==PDTYPE) ){
    if( (realsampled==0)||(nogridchoice==1) ){ // no interp
      scain=s;
      vecin=v;
    }
    else{ // must reference off 0 for interp function
      // for interp of dumps
      if(FULLOUTPUT==0){
				scain=s;
				vecin=v;
      }
      else if(FULLOUTPUT==1){
				// assumes 2 boundary zones
				scain=(FTYPE (*) [N3M][N2M][N1M])(&(sa[-1][N3BND/2][N2BND/2][N1BND/2]));
				vecin=(FTYPE (*) [3][N3M][N2M][N1M])(&(va[-1][-1][N3BND/2][N2BND/2][N1BND/2]));
      }
      else if(FULLOUTPUT==2){
				scain=(FTYPE (*) [N3M][N2M][N1M])(&(sa[-1][0][0][0]));
				vecin=(FTYPE (*) [3][N3M][N2M][N1M])(&(va[-1][-1][0][0][0]));
      }
    }
  }
  else if(which==ADTYPE){
    if( (realsampled==0)||(nogridchoice==1) ){ // no interp
      scain=sanal;
      vecin=vanal;
    }
    else{ // must reference off 0 for interp function
      // for interp of dumps
      if(FULLOUTPUT==0){
				scain=sanal;
				vecin=vanal;
      }
      else if(FULLOUTPUT==1){
				// assumes 2 boundary zones
				scain=(FTYPE (*) [N3M][N2M][N1M])(&(sanala[-1][N3BND/2][N2BND/2][N1BND/2]));
				vecin=(FTYPE (*) [3][N3M][N2M][N1M])(&(vanala[-1][-1][N3BND/2][N2BND/2][N1BND/2]));
      }
      else if(FULLOUTPUT==2){
				scain=(FTYPE (*) [N3M][N2M][N1M])(&(sanala[-1][0][0][0]));
				vecin=(FTYPE (*) [3][N3M][N2M][N1M])(&(vanala[-1][-1][0][0][0]));
      }
    }
  }
  else if(which==FLTYPE){
#if(FLOORDUMPFLAGMEMORY)
    if( (realsampled==0)||(nogridchoice==1) ){ // no interp
      scain=floorvars;
      vecin=floorvarv;
    }
    else{ // must reference off 0 for interp function
      // for interp of dumps
      if(FULLOUTPUT==0){
				scain=floorvars;
				vecin=floorvarv;
      }
      else if(FULLOUTPUT==1){
				// assumes 2 boundary zones
				scain=(FTYPE (*) [N3M][N2M][N1M])(&(floorvarsa[-1][N3BND/2][N2BND/2][N1BND/2]));
				vecin=(FTYPE (*) [3][N3M][N2M][N1M])(&(floorvarva[-1][-1][N3BND/2][N2BND/2][N1BND/2]));
      }
      else if(FULLOUTPUT==2){
				scain=(FTYPE (*) [N3M][N2M][N1M])(&(floorvarsa[-1][0][0][0]));
				vecin=(FTYPE (*) [3][N3M][N2M][N1M])(&(floorvarva[-1][-1][0][0][0]));
      }
    }
#endif
    // floorvar0 directly
  }
  else if(which==NPTYPE){
#if(VISCMEM)
    // compute sigma
    delv = work1;
    nurho_real = work2;
    if((computesigma==1)||(npoldschool==0)){ // must compute if not read in, must have parameters same for run
      tdep_compute();
      nu_compute();
    }
    if(computesigma==1){
      compute_sigma_gen(sigma,rost,rostnu,nurho_real,delv);
    }
    if( (realsampled==0)||(nogridchoice==1) ){ // no interp
			inu_real = nu_real;
    }
    else{ // must reference off 0 for interp function
      // for interp of dumps
      if(FULLOUTPUT==0){
				inu_real = nu_real;
      }
      else if(FULLOUTPUT==1){
				// assumes 2 boundary zones
				isigma=(FTYPE (*) [3][N3M][N2M][N1M])(&(sigmaa[-1][-1][N3BND/2][N2BND/2][N1BND/2]));
				inu_real = (FTYPE (*) [N2M][N1M])(&(nu_reala[N3BND/2][N2BND/2][N1BND/2])); 
      }
      else if(FULLOUTPUT==2){
				isigma=(FTYPE (*) [3][N3M][N2M][N1M])(&(sigmaa[-1][-1][0][0][0]));
				inu_real = (FTYPE (*) [N2M][N1M])(&(nu_reala[0][0][0])); 
      }
    } 
#endif
  }
  else if(which==AVG2DTYPE){
#if(DOAVGDIAGMEMORY)
    // use avg2_3[NUMAVG2_3][..][..], 2D data with NUMAVG2_3 types
    delv = work1;
    nurho_real = work2;
    if(avg2doldschool==1){ // hack to compute from older avg2d version data, uses avg data which was written to normal scalars/vectors on input
      // compute sigma
      //alpha_real=alpha_real0; // assumes average is while viscosity is on!
      tdep_compute();
      nu_compute();
      compute_sigma_gen(sigma,rost,rostnu,nurho_real,delv);
      // now assign back to avg data
      LOOPF{
				avg2_3[9][j][i]=sigma[1][1][k][j][i];
				avg2_3[10][j][i]=sigma[1][2][k][j][i];
				avg2_3[11][j][i]=sigma[1][3][k][j][i];
				avg2_3[12][j][i]=sigma[2][2][k][j][i];
				avg2_3[13][j][i]=sigma[2][3][k][j][i];
				avg2_3[14][j][i]=sigma[3][3][k][j][i];
				avg2_3[15][j][i]=nu_real[k][j][i];
      }
    }
    if( (realsampled==0)||(nogridchoice==1) ){ // no interp
      iavg2_3=avg2_3;
    }
    else{ // must reference off 0 for interp function
      // for interp of dumps
      if(FULLOUTPUT==0){
				iavg2_3=avg2_3;
      }
      else if(FULLOUTPUT==1){
				// assumes 2 boundary zones
				iavg2_3=(FTYPE (*) [N2M][N1M])(&(avg2_3a[-1][N2BND/2][N1BND/2])); 
      }
      else if(FULLOUTPUT==2){
				iavg2_3=(FTYPE (*) [N2M][N1M])(&(avg2_3a[-1][0][0])); 
      }
    }
#endif
  }
  else if(which==CALCTYPE){
    delv = work1;
    nurho_real = work2;
    // compute sigma
#if(!PPCLEAN)
    tdep_compute();
    nu_compute();
    compute_sigma_gen(sigma,rost,rostnu,nurho_real,delv);
#endif
    ///////////////////
    // 
    //  Compute various interesting things pp only, no read-in possible except using sm
    //
    // assume all values wanted as zone centered.  Since mixing occurs, no native grid location
    LOOPF{
      if(i<=N1){
				vxa = e2z_1(v[1][1],k,j,i);
      }
      else vxa=v[1][1][k][j][i];
      if(j<=N2){
				vya = e2z_2(v[1][2],k,j,i);
      }
      else vya = v[1][2][k][j][i];
      vza = v[1][3][k][j][i] ;  // When fake-3d
      vtot2a=vxa*vxa+vya*vya+vza*vza;

      enthalpy=gam*s[2][k][j][i]/s[1][k][j][i]; // really specific enthalpy: enthalpy per unit mass, needed as used for this calculation
      pot=s[3][k][j][i];

#if(VISCMEM)
      sig11a=sigma[1][1][k][j][i];
      sig12a=sig21a=c2z_3(sigma[1][2],k,j,i);
      sig13a=sig31a=c2z_2(sigma[1][3],k,j,i);
      sig22a=sigma[2][2][k][j][i];
      sig23a=sig32a=c2z_1(sigma[2][3],k,j,i);
      sig33a=sigma[3][3][k][j][i];
#endif
      //fprintf(stdout,"%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",sig11a,sig12a,sig13a,sig22a,sig23a,sig33a);
      //fprintf(stdout,"%15.10g %15.10g %15.10g\n",nu_real[k][j][i],alpha_real,alpha_real0);
      // argument of div for ang mom flux, no visc
      ppcalcs[1][k][j][i]=G3(2,i)*G4(2,j)*s[1][k][j][i]*vxa*vza; // r-dir
      ppcalcs[2][k][j][i]=G3(2,i)*G4(2,j)*s[1][k][j][i]*vya*vza; // theta dir
      ppcalcs[3][k][j][i]=G3(2,i)*G4(2,j)*s[1][k][j][i]*vza*vza; // phi dir

#if(VISCMEM)
      // argument of div for visc component of ang mom flux
      ppcalcs[4][k][j][i]=G3(2,i)*G4(2,j)*sig13a; // r-dir
      ppcalcs[5][k][j][i]=G3(2,i)*G4(2,j)*sig23a; // theta-dir
      ppcalcs[6][k][j][i]=G3(2,i)*G4(2,j)*sig33a; // phi-dir
#endif

      tempf=0.5*vtot2a+enthalpy+pot;
      // argument of div for energy flux, no visc
      ppcalcs[7][k][j][i]=tempf*s[1][k][j][i]*vxa; // r-dir
      ppcalcs[8][k][j][i]=tempf*s[1][k][j][i]*vya; // theta-dir
      ppcalcs[9][k][j][i]=tempf*s[1][k][j][i]*vza; // phi-dir

#if(VISCMEM)
      // argument of div for energy flux, visc component
      ppcalcs[10][k][j][i]=sig13a*vza+sig12a*vya+sig11a*vxa; // r-dir
      ppcalcs[11][k][j][i]=sig23a*vza+sig22a*vya+sig21a*vxa; // theta-dir
      ppcalcs[12][k][j][i]=sig33a*vza+sig32a*vya+sig31a*vxa; // phi-dir
#endif
    }
#if(POSTPROC)
    if( (realsampled==0)||(nogridchoice==1) ){ // no interp
      ippcalcs=ppcalcs;
    }
    else{ // must reference off 0 for interp function
      // for interp of dumps
      if(FULLOUTPUT==0){
				ippcalcs=ppcalcs;
      }
      else if(FULLOUTPUT==1){
				// assumes 2 boundary zones
				ippcalcs=(FTYPE (*) [N3M][N2M][N1M])(&(ppcalcsa[-1][N3BND/2][N2BND/2][N1BND/2])); 
      }
      else if(FULLOUTPUT==2){
				ippcalcs=(FTYPE (*) [N3M][N2M][N1M])(&(ppcalcsa[-1][0][0][0])); 
      }
    }
#endif
  }
  else if(which==FLINETYPE){
    differva = workv1;
    differvb = workv2;
    ax3 = work1;
    ax3a = work1a;
    // GODMARK -- bad loops for boundtype==2?
    ///////////////////
    // 
    //  Compute phi component of vector potential in coordinate basis, such that contour plot is field lines
    //
    // use ppcalcs arrays here too
    // 
    // only do B-field(not v-field) and only output x3-component(x2 if interpolated in coord==3)
    //
    l=2; // just b-field
    if((POSTPROC)&&(dumpsmdump)){ // all values wanted as zone centered(integral ends up at zone center)
      LOOPF{
				differva[l][k][j][i]=-v[l][2][k][j][i]*G3(2,i)*G4(2,j)*dx[1][1][i]; // a
				differvb[l][k][j][i]=v[l][1][k][j][i]*G2(2,i)*G3(2,i)*G4(2,j)*dx[1][2][j]; // b
				ax3[k][j][i]=0;
      }
    }
    else{ // standard grid location (integral ends up on vector locations, or effectively at corners) (if dumpsmdump, later interpolated to zone center)
      LOOPF{
				differva[l][k][j][i]=-v[l][2][k][j][i]*G3(2,i)*G4(1,j)*dx[1][1][i]; // a
				differvb[l][k][j][i]=v[l][1][k][j][i]*G2(1,i)*G3(1,i)*G4(2,j)*dx[1][2][j]; // b
				ax3[k][j][i]=0;
      }
    }
    if((mycpupos[1]==ncpux1-1)&&(mycpupos[2]==0)){
      for(k=INFULL3;k<OUTFULL3;k++){
				ax3[k][0][N1]=-differva[l][k][0][N1]; // reference for all cpus ( shouldn't just choose 0)
      }
    }
    else{
      for(k=INFULL3;k<OUTFULL3;k++){
				ax3[k][0][N1]=0; // reference, or overridden by other cpus (reference for mycpupos[1]=ncpux1-1 and mycpupos[2]=0 cpu)
      }
    }

    for(k=0;k<N3;k++) for(j=0;j<=N2;j++) for(i=0;i<=N1;i++){
      // do x1 part of integral
      jj=0; kk=k;
      for(ii=N1-1;ii>=i;ii--){ // skips i==N1 since that line is not determined by x1 integral
				ax3[k][j][i]-=differva[l][kk][jj][ii];
      }
      // x2 integral
      ii=i; kk=k;
      for(jj=0;jj<=j-1;jj++){ // skips j==0
				ax3[k][j][i]+=differvb[l][kk][jj][ii];
      }
    }
    if(numprocs>1){
#if(USEMPI)
      // need to ship off ax3 "corners" to other cpus and renormalize the solution
      for(k=0;k<N3;k++){
				// send down j-line only if on outer i-edge
				if((mycpupos[1]==ncpux1-1)&&(mycpupos[2]==0)){
					totalnorm=ax3[k][0][N1];
				}
				else totalnorm=0;

				if(mycpupos[1]==ncpux1-1){
					for(j=0;j<ncpux2-1;j++){
						if(mycpupos[2]==j){
							sendnorm=totalnorm+ax3[k][N2][N1];

							MPI_Isend(&sendnorm,1,MPI_FTYPE,mycpupos[2]+1,0,combound[0],&srequest);
							MPI_Wait(&srequest,&mpichstatus);
						}
						if(mycpupos[2]==j+1){
							MPI_Irecv(&recvnorm,1,MPI_FTYPE,mycpupos[2]-1,0,combound[0],&rrequest);
							MPI_Wait(&rrequest,&mpichstatus);
							totalnorm=totalnorm+recvnorm;
						}
						MPI_Barrier(combound[0]); // must wait so can add up sequentially
					}
				}
				MPI_Barrier(MPI_COMM_WORLD); // hold all cpus till here
				// send down i-line for all j
				for(i=ncpux1-1;i>0;i--){
					if(mycpupos[1]==i){
						sendnorm=totalnorm+ax3[k][0][0];
						MPI_Isend(&sendnorm,1,MPI_FTYPE,myid-1,myid,MPI_COMM_WORLD,&srequest);
						MPI_Wait(&srequest,&mpichstatus);
					}
					if(mycpupos[1]==i-1){
						MPI_Irecv(&recvnorm,1,MPI_FTYPE,myid+1,myid+1,MPI_COMM_WORLD,&rrequest);
						MPI_Wait(&rrequest,&mpichstatus);
						totalnorm=totalnorm+recvnorm;
					}
					MPI_Barrier(MPI_COMM_WORLD); // must wait for all cpus since no per-i (const j/k) communicator
				}
      }
      // now renormalize
      for(k=0;k<N3;k++) for(j=0;j<=N2;j++) for(i=0;i<=N1;i++)
				{
					ax3[k][j][i]+=totalnorm;
				}
      if(dumpsmdump&&(!POSTPROC)){
				for(k=0;k<N3;k++) for(j=0;j<=N2;j++) for(i=0;i<=N1;i++){
					ax3[k][j][i]=c2z_3(ax3,k,j,i); // can do this since no dependencies
				}
      }// otherwise effectively on corners
 
#endif // end if MPI
    }
 

    // the field lines are actually the coord. basis ax3, for which v\cdot\nabla Aphi=0 and B\cdot\nabla Aphi=0
    //LOOP{
    //ax3[k][j][i]/=(G3(2,i)*G4(2,j)); // revert to orthonormal version(assumes zone centered)
    // }
    if( (realsampled==0)||(nogridchoice==1) ){ // no interp
      iax3=ax3;
    }
    else{ // must reference off 0 for interp function
      // for interp of dumps
      if(FULLOUTPUT==0){
				iax3=ax3;
      }
      else if(FULLOUTPUT==1){
				// assumes 2 boundary zones
				iax3=(FTYPE (*) [N2M][N1M])(&(ax3a[N3BND/2][N2BND/2][N1BND/2])); 
      }
      else if(FULLOUTPUT==2){
				iax3=(FTYPE (*) [N2M][N1M])(&(ax3a[0][0][0])); 
      }
    }
  }
  else{
    fprintf(fail_file,"not setup for which=%d in dump\n",which);
    myexit(1);
  }






  //////////////////////////////////////////////////////
  //
  // begin real output of data to a file
  //
  //////////////////////////////////////////////////////




  ///////////////////////////////////
  // OUTPUT NON_INTERP DATA

  if( (realsampled==0)||(nogridchoice==1) ){ // i.e. no interp output
    // loop over choicen output region for input=output data

    LOOPDIAGOUTPUT(FULLOUTPUT,FULLOUTPUT,FULLOUTPUT){
      BUFFERINIT;

      if((mpicombine==0)&&(!binaryoutput)) fprintf(fp," ");
      if((which==DTYPE)||(which==PDTYPE)||(which==ADTYPE)||(which==FLTYPE) ){
				if((mpicombine==0)&&(!binaryoutput)){
					for(l=1;l<=NUMSCA;l++){
					  dump1data(CENT,i,j,k,fp,scain[l][k][j][i]);
						    //fprintf(fp,"%21.15g ",scain[l][k][j][i]);
					}
	  
					for(l=1;l<=REALNUMVEC;l++){
						if(which==FLTYPE){
							// already in middle of zone
							dump1data(CENT,i,j,k,fp,floorvar0[l][k][j][i]);
						}
						for(m=1;m<=3;m++){
							if(dumpsmdump&&(nogridchoice==0) ){
								if(m==1){
									if(i==OUTFULL1-1){
										tempf=vecin[l][m][k][j][i];
									}
									else tempf=e2z_1(vecin[l][m],k,j,i);
								}
								else if(m==2){
									if(j==OUTFULL2-1){
										tempf=vecin[l][m][k][j][i];
									}
									else tempf=e2z_2(vecin[l][m],k,j,i);
								}
								else if(m==3){
									if(k==OUTFULL3-1){
										tempf=vecin[l][m][k][j][i];
									}
									else tempf=vecin[l][m][k][j][i];
								}
							}
							else tempf=vecin[l][m][k][j][i];
							dump1data(m,i,j,k,fp,tempf);
						}
					}
				}
				else if((mpicombine==0)&&(binaryoutput)){
#if(FLOATTYPE==0)
#define WHICHDUM (dumf)
#elif(FLOATTYPE==1)
#define WHICHDUM (dumlf)
#endif
					// binary mode
					if(!dumpsmdump){
						if(analoutput!=13){
							WHICHDUM=scain[1][k][j][i];
							dump1databin(CENT,i,j,k,&WHICHDUM, sizeof(FTYPE), 1, fp);
							//							fwrite(&WHICHDUM,sizeof(FTYPE),1,fp);
						}
						if(analoutput!=13){ // limited variables needed
							WHICHDUM=scain[2][k][j][i];
							dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
							WHICHDUM=scain[3][k][j][i];
							dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
						}
						if(analoutput!=13){ // limited variables needed
							if(which==FLTYPE){
								WHICHDUM=floorvar0[1][k][j][i];
								dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
							}
						}
						WHICHDUM=vecin[1][1][k][j][i];
						dump1databin(VDIR1,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
						WHICHDUM=vecin[1][2][k][j][i];
						dump1databin(VDIR2,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
						WHICHDUM=vecin[1][3][k][j][i];
						dump1databin(VDIR3,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);

						if(analoutput!=13){ // limited variables needed
							if(which==FLTYPE){
								WHICHDUM=floorvar0[2][k][j][i];
								dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
							}
						}
	    
						WHICHDUM=vecin[2][1][k][j][i];
						dump1databin(VDIR1,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
						WHICHDUM=vecin[2][2][k][j][i];
						dump1databin(VDIR2,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
						WHICHDUM=vecin[2][3][k][j][i];
						dump1databin(VDIR3,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);

						if(analoutput==13){
							WHICHDUM=scain[1][k][j][i];
							dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
						}

					}
					else{
						if(analoutput!=13){
							WHICHDUM=scain[1][k][j][i];
							dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
						}
						if(analoutput!=13){
							WHICHDUM=scain[2][k][j][i];
							dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
							WHICHDUM=scain[3][k][j][i];
							dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
						}
						if(analoutput!=13){
							if(which==FLTYPE){
								WHICHDUM=floorvar0[1][k][j][i];
								dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
							}
						}

						WHICHDUM=e2z_1(vecin[1][1],k,j,i);
						dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
						WHICHDUM=e2z_2(vecin[1][2],k,j,i);
						dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
						WHICHDUM=e2z_3(vecin[1][3],k,j,i);
						dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);

						if(analoutput!=13){
							if(which==FLTYPE){
								WHICHDUM=floorvar0[2][k][j][i];
								dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
							}
						}
	    
						WHICHDUM=e2z_1(vecin[2][1],k,j,i);
						dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
						WHICHDUM=e2z_2(vecin[2][2],k,j,i);
						dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
						WHICHDUM=e2z_3(vecin[2][3],k,j,i);
						dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);

						if(analoutput==13){
							WHICHDUM=scain[1][k][j][i];
							dump1databin(CENT,i,j,k,&WHICHDUM,sizeof(FTYPE),1,fp);
						}
					}
				}
				else if(mpicombine==1){// then MPI output in binary mode(forced)
					// fill writebuf with data
					if(analoutput!=13){
						if(!dumpsmdump){
							// rho
    						        writebuf1data(CENT,i,j,k,scain[1][k][j][i],&writebuf[BUFFERMAP]);

							// en
							writebuf1data(CENT,i,j,k,scain[2][k][j][i],&writebuf[BUFFERMAP]);
							// pot
							writebuf1data(CENT,i,j,k,scain[3][k][j][i],&writebuf[BUFFERMAP]);
							if(which==FLTYPE){
								writebuf1data(CENT,i,j,k,floorvar0[1][k][j][i],&writebuf[BUFFERMAP]);
							}
							// v1
							writebuf1data(VDIR1,i,j,k,vecin[1][1][k][j][i],&writebuf[BUFFERMAP]);
							// v2
							writebuf1data(VDIR2,i,j,k,vecin[1][2][k][j][i],&writebuf[BUFFERMAP]);
							// v3
							writebuf1data(VDIR3,i,j,k,vecin[1][3][k][j][i],&writebuf[BUFFERMAP]);
							if(which==FLTYPE){
								writebuf1data(CENT,i,j,k,floorvar0[2][k][j][i],&writebuf[BUFFERMAP]);
							}
							// B1
							writebuf1data(VDIR1,i,j,k,vecin[2][1][k][j][i],&writebuf[BUFFERMAP]);
							// B2
							writebuf1data(VDIR2,i,j,k,vecin[2][2][k][j][i],&writebuf[BUFFERMAP]);
							// B3
							writebuf1data(VDIR3,i,j,k,vecin[2][3][k][j][i],&writebuf[BUFFERMAP]);
						}
						else{
							// rho
							writebuf1data(CENT,i,j,k,scain[1][k][j][i],&writebuf[BUFFERMAP]);
							// en
							writebuf1data(CENT,i,j,k,scain[2][k][j][i],&writebuf[BUFFERMAP]);
							// pot
							writebuf1data(CENT,i,j,k,scain[3][k][j][i],&writebuf[BUFFERMAP]);
							// v1
							writebuf1data(CENT,i,j,k,e2z_1(vecin[1][1],k,j,i),&writebuf[BUFFERMAP]);
							// v2
							writebuf1data(CENT,i,j,k,e2z_2(vecin[1][2],k,j,i),&writebuf[BUFFERMAP]);
							// v3
							writebuf1data(CENT,i,j,k,e2z_3(vecin[1][3],k,j,i),&writebuf[BUFFERMAP]);
							// B1
							writebuf1data(CENT,i,j,k,e2z_1(vecin[2][1],k,j,i),&writebuf[BUFFERMAP]);
							// B2
							writebuf1data(CENT,i,j,k,e2z_2(vecin[2][2],k,j,i),&writebuf[BUFFERMAP]);
							// B3
							writebuf1data(CENT,i,j,k,e2z_3(vecin[2][3],k,j,i),&writebuf[BUFFERMAP]);
						}
					}
					else{
						if(!dumpsmdump){
							// v1
							writebuf1data(VDIR1,i,j,k,vecin[1][1][k][j][i],&writebuf[BUFFERMAP]);
							// v2
							writebuf1data(VDIR2,i,j,k,vecin[1][2][k][j][i],&writebuf[BUFFERMAP]);
							// v3
							writebuf1data(VDIR3,i,j,k,vecin[1][3][k][j][i],&writebuf[BUFFERMAP]);
							// B1
							writebuf1data(VDIR1,i,j,k,vecin[2][1][k][j][i],&writebuf[BUFFERMAP]);
							// B2
							writebuf1data(VDIR2,i,j,k,vecin[2][2][k][j][i],&writebuf[BUFFERMAP]);
							// B3
							writebuf1data(VDIR3,i,j,k,vecin[2][3][k][j][i],&writebuf[BUFFERMAP]);
							// rho
							writebuf1data(CENT,i,j,k,scain[1][k][j][i],&writebuf[BUFFERMAP]);
						}
						else{
							// v1
						        writebuf1data(CENT,i,j,k,e2z_1(vecin[1][1],k,j,i),&writebuf[BUFFERMAP]);
							// v2
							writebuf1data(CENT,i,j,k,e2z_2(vecin[1][2],k,j,i),&writebuf[BUFFERMAP]);
							// v3
							writebuf1data(CENT,i,j,k,e2z_3(vecin[1][3],k,j,i),&writebuf[BUFFERMAP]);
							// B1
							writebuf1data(CENT,i,j,k,e2z_1(vecin[2][1],k,j,i),&writebuf[BUFFERMAP]);
							// B2
							writebuf1data(CENT,i,j,k,e2z_2(vecin[2][2],k,j,i),&writebuf[BUFFERMAP]);
							// B3
							writebuf1data(CENT,i,j,k,e2z_3(vecin[2][3],k,j,i),&writebuf[BUFFERMAP]);
							// rho
							writebuf1data(CENT,i,j,k,scain[1][k][j][i],&writebuf[BUFFERMAP]);
						}
					}// end if analoutput==13
				}// end else mpicombine
      }// end dtypes
      else if(which==NPTYPE){
				if(VISCMEM){
					for(l=1;l<=3;l++){
						for(m=1;m<=3;m++){
							if( ((l==1)&&(m==3))|| ((l==3)&&(m==1)) ){
								if((i==OUTFULL1-1)||(dumpsmdump&&(nogridchoice==0)) ){
									ftempT[3]=sigma[l][m][k][j][i];
								}
								else ftempT[3]=c2z_2(sigma[l][m],k,j,i);
							}
							else if( ((l==2)&&(m==3)) || ((l==3)&&(m==2)) ){
								if((j==OUTFULL2-1)||(dumpsmdump&&(nogridchoice==0)) ){
									ftempT[5]=sigma[l][m][k][j][i];
								}
								else ftempT[5]=c2z_1(sigma[l][m],k,j,i);
							}
							else if( ((l==2)&&(m==1)) || ((l==1)&&(m==2)) ){
								if((i==OUTFULL1-1)||(j==OUTFULL2-1)||(dumpsmdump&&(nogridchoice==0)) ){
									ftempT[2]=sigma[l][m][k][j][i];
								}
								else ftempT[2]=c2z_3(sigma[l][m],k,j,i);
							}
							else if( ((l==1)&&(m==1)) ){
								ftempT[1]=sigma[l][m][k][j][i];
							}
							else if( ((l==2)&&(m==2)) ){
								ftempT[4]=sigma[l][m][k][j][i];
							}
							else if( ((l==3)&&(m==3)) ){
								ftempT[6]=sigma[l][m][k][j][i];
							}
						}// over m
					}// over l
					if(mpicombine==0){
						for(l=1;l<=6;l++){ // number of components in sym tensor
						  dump1data(CENT,i,j,k,fp,ftempT[l]);
						}
						dump1data(CENT,i,j,k,fp,nu_real[k][j][i]);	
					}
					else{
						for(l=1;l<=6;l++){ // number of components in sym tensor
							writebuf1data(CENT,i,j,k,ftempT[l],&writebuf[BUFFERMAP]);
						}
						writebuf1data(CENT,i,j,k,nu_real[k][j][i],&writebuf[BUFFERMAP]);
					}
				}// if sigma dump
				if(1){ //forces,fluxes, and energy
					if(mpicombine==0){
						for(l=1;l<=18;l++){// number of terms
							dump1data(CENT,i,j,k,fp,ffv_calcs(1,l,i,j,k));
						}
						for(l=1;l<=36;l++){// number of terms
							dump1data(CENT,i,j,k,fp,ffv_calcs(2,l,i,j,k));
						}
						for(l=1;l<=18;l++){// number of terms
							dump1data(CENT,i,j,k,fp,ffv_calcs(3,l,i,j,k));
						}
					}
					else{
						for(l=1;l<=18;l++){// number of terms
							writebuf1data(CENT,i,j,k,ffv_calcs(1,l,i,j,k),&writebuf[BUFFERMAP]);
						}
						for(l=1;l<=36;l++){// number of terms
							writebuf1data(CENT,i,j,k,ffv_calcs(2,l,i,j,k),&writebuf[BUFFERMAP]);
						}
						for(l=1;l<=18;l++){// number of terms
							writebuf1data(CENT,i,j,k,ffv_calcs(3,l,i,j,k),&writebuf[BUFFERMAP]);
						}
					}
				}
      }
      else if(which==AVG2DTYPE){// avg2d dump
				// already DUMPSM'ed in average routine if needed
				if(mpicombine==0){
				  dumpcheckfprintf(CENT,i,j,k,fp,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g",avg2_3[1][j][i],avg2_3[2][j][i],avg2_3[3][j][i],avg2_3[4][j][i],avg2_3[5][j][i],avg2_3[6][j][i],avg2_3[7][j][i],avg2_3[8][j][i],avg2_3[9][j][i],avg2_3[10][j][i],avg2_3[11][j][i],avg2_3[12][j][i],avg2_3[13][j][i],avg2_3[14][j][i],avg2_3[15][j][i]);
				}
				else{
					for(l=1;l<=15;l++){
						writebuf1data(CENT,i,j,k,avg2_3[l][j][i],&writebuf[BUFFERMAP]); // only applies in 2D
					}
				}
      }// if avg2d dump
      else if(which==CALCTYPE){// calctype dump
				// already DUMPSM'ed above
				if(mpicombine==0){
					dumpcheckfprintf(CENT,i,j,k,fp,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g",ppcalcs[1][k][j][i],ppcalcs[2][k][j][i],ppcalcs[3][k][j][i],ppcalcs[4][k][j][i],ppcalcs[5][k][j][i],ppcalcs[6][k][j][i],ppcalcs[7][k][j][i],ppcalcs[8][k][j][i],ppcalcs[9][k][j][i],ppcalcs[10][k][j][i],ppcalcs[11][k][j][i],ppcalcs[12][k][j][i]);
				}
				else{
					for(l=1;l<=12;l++){
						writebuf1data(CENT,i,j,k,ppcalcs[l][k][j][i],&writebuf[BUFFERMAP]);
					}
				}
      }// end if calctype dump
      else if(which==FLINETYPE){// fieldlinetype dump
				if(mpicombine==0){
					// SUPERMARK
					//	fprintf(fp,"%21.15g",ppcalcs[1][k][j][i]);
					dumpcheckfprintf(CENT,i,j,k,fp,"%.6g",ax3[k][j][i]); // assume for now don't compute anything from it, just display, so don't need all precision
				}
				else{
					writebuf1data(CENT,i,j,k,ax3[k][j][i],&writebuf[BUFFERMAP]);
				}
      }// end if fieldlinetype dump

      if((mpicombine==0)&&(!binaryoutput)) fprintf(fp,"\n");
    }
  }
  /////////////////////
  //
  //  ELSE IF DOING    INTERPOLATION (2d only right now) (and only in 1 cpu process, so no mpi combine stuff)
  //
  //
  else{
    // here loop on output is always over same chosen region, using input region defined by FULLOUTPUT.  This is done instead of making interpolate() deal with offsetting input data

    if((which==DTYPE)||(which==PDTYPE)||(which==ADTYPE)||(which==FLTYPE) ){
      if(POSTPROC){ printf(".i%d.",which); fflush(stdout); }
      for(i=1;i<=NUMSCA;i++){
				interpolate(0,i,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,scain[i][0],iq[i],douterdefs[i],outtype);
				if(POSTPROC){ printf("s"); fflush(stdout); }
      }
      for(i=1;i<=REALNUMVEC;i++){
				if(which==FLTYPE){
					interpolate(0,i,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,floorvar0[i][0],work0iq[i],douterdef0[i],outtype);
					if(POSTPROC){ printf("s"); fflush(stdout); }
				}
				for(j=1;j<=3;j++){
					if(j==1) itemp=1;
					else if(j==2) itemp=2;
					else if(j==3) itemp=0;
					interpolate(itemp,i,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,vecin[i][j][0],viq[i][j],douterdefv[i][j],outtype);
					if(POSTPROC){ printf("v"); fflush(stdout); }
				}
      }
    }
    else if(which==NPTYPE){
      if(POSTPROC){ printf(".i%d.",which); fflush(stdout); }
      for(p=1;p<=6;p++){ // hack the sigma output, should always have enough in viq for symmetric tensor hack
				if(p==1){
					l=1; m=1; i=1; j=1; itemp=0;
				}
				if(p==2){
					l=1; m=2; i=1; j=2; itemp=3;
				}
				if(p==3){
					l=1; m=3; i=1; j=3; itemp=1;
				}
				if(p==4){
					l=2; m=2; i=2; j=1; itemp=0;
				}
				if(p==5){
					l=2; m=3; i=2; j=2; itemp=2;
				}
				if(p==6){
					l=3; m=3; i=2; j=3; itemp=0;
				}
				interpolate(itemp,p,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,isigma[l][m][0],viq[i][j],0.0,outtype); // force 0.0 outerdef
				if(POSTPROC){ printf("t"); fflush(stdout); }
      }
      itemp=0;
      interpolate(itemp,1,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,inu_real[0],iq[1],SSMALL,outtype); // force SSMALL outerdef(no 0)
      if(POSTPROC){ printf("s"); fflush(stdout); }

    }
    else if(which==AVG2DTYPE){
      if(POSTPROC){ printf(".i%d.",which); fflush(stdout); }
      //      iq=workiqavg; // just use work space natively
      // set new outerdefs for avg2d data
      odef[1]=douterdefs[1];
      odef[2]=douterdefs[2];
      odef[3]=0.0;
      odef[4]=0.0;
      odef[5]=0.0;
      odef[6]=douterdefv[1][1];
      odef[7]=douterdefv[1][2];
      odef[8]=douterdefv[1][3];
      odef[9]=0;
      odef[10]=0;
      odef[11]=0;
      odef[12]=0;
      odef[13]=0;
      odef[14]=0;
      odef[15]=SSMALL;
      // do scalar type stuff
      for(i=1;i<=5;i++){
				interpolate(0,i,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,iavg2_3[i],workiqavg[i],odef[i],outtype);
				if(POSTPROC){ printf("s"); fflush(stdout); }
      }
      // do vector type stuff
      // vx1
      i=6;
      if(dumpsmdump) itemp=0; else itemp=1;
      interpolate(itemp,3,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,iavg2_3[i],workiqavg[i],odef[i],outtype);
      if(POSTPROC){ printf("v"); fflush(stdout); }
      // vx2
      i=7;
      if(dumpsmdump) itemp=0; else itemp=2;
      interpolate(itemp,4,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,iavg2_3[i],workiqavg[i],odef[i],outtype);
      if(POSTPROC){ printf("v"); fflush(stdout); }
      // vx3
      i=8;
      if(dumpsmdump) itemp=0; else itemp=0;
      interpolate(itemp,5,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,iavg2_3[i],workiqavg[i],odef[i],outtype);
      if(POSTPROC){ printf("v"); fflush(stdout); }
      // sig11
      i=9;
      if(dumpsmdump) itemp=0; else itemp=0;
      interpolate(itemp,3,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,iavg2_3[i],workiqavg[i],odef[i],outtype);
      if(POSTPROC){ printf("t"); fflush(stdout); }
      // sig12
      i=10;
      if(dumpsmdump) itemp=0; else itemp=3;
      interpolate(itemp,3,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,iavg2_3[i],workiqavg[i],odef[i],outtype);
      if(POSTPROC){ printf("t"); fflush(stdout); }
      // sig13
      i=11;
      if(dumpsmdump) itemp=0; else itemp=1;
      interpolate(itemp,3,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,iavg2_3[i],workiqavg[i],odef[i],outtype);
      if(POSTPROC){ printf("t"); fflush(stdout); }
      // sig22
      i=12;
      if(dumpsmdump) itemp=0; else itemp=0;
      interpolate(itemp,3,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,iavg2_3[i],workiqavg[i],odef[i],outtype);
      if(POSTPROC){ printf("t"); fflush(stdout); }
      // sig23
      i=13;
      if(dumpsmdump) itemp=0; else itemp=2;
      interpolate(itemp,3,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,iavg2_3[i],workiqavg[i],odef[i],outtype);
      if(POSTPROC){ printf("t"); fflush(stdout); }
      // sig33
      i=14;
      if(dumpsmdump) itemp=0; else itemp=0;
      interpolate(itemp,3,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,iavg2_3[i],workiqavg[i],odef[i],outtype);
      if(POSTPROC){ printf("t"); fflush(stdout); }
      // nu_real
      i=15;
      interpolate(0,1,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,iavg2_3[i],workiqavg[i],odef[i],outtype);
      if(POSTPROC){ printf("s"); fflush(stdout); }
    }// end if doing avg2d
    else if(which==CALCTYPE){
      if(POSTPROC){ printf(".i%d.",which); fflush(stdout); }
      // set new outerdefs for data
      odef[1]=0.0;
      odef[2]=0.0;
      odef[3]=0.0;
      odef[4]=0.0;
      odef[5]=0.0;
      odef[6]=0.0;
      odef[7]=0.0;
      odef[8]=0.0;
      odef[9]=0;
      odef[10]=0;
      odef[11]=0;
      odef[12]=0;
      // do vector type stuff
      for(i=1;i<=10;i+=3){ // do each component in a group
				// vx1
				interpolate(0,3,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,ippcalcs[i][0],workiqavg[i],odef[i],outtype);
				if(POSTPROC){ printf("v"); fflush(stdout); }
				// vx2
				interpolate(0,4,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,ippcalcs[i+1][0],workiqavg[i+1],odef[i+1],outtype);
				if(POSTPROC){ printf("v"); fflush(stdout); }
				// vx3
				interpolate(0,5,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,ippcalcs[i+2][0],workiqavg[i+2],odef[i+2],outtype);
				if(POSTPROC){ printf("v"); fflush(stdout); }
      }
    }// end if doing calctype dump
    else if(which==FLINETYPE){
      if(POSTPROC){ printf(".i%d.",which); fflush(stdout); }
      // set new outerdefs for data
      odef[1]=0.0;
      // do vector type stuff
      i=1;
      interpolate(0,5,realsampled,SAMEGRIDD,realn1,realn2,DUMN1,DUMN2,iax3[0],workiqavg[i],odef[i],outtype);
      if(POSTPROC){ printf("v"); fflush(stdout); }
    }// end if doing fieldlinetype dump

    if(POSTPROC){ printf(".o."); fflush(stdout); }
    // create conversion matrix to speed up things.  Needed to convert vectors when going from 1 coord system to another
#if(POSTPROC==1)
    if((COORD==3)&&(outtype>0)){
      if(firsttime==1){
				LOOPD{	  
				  k=0; // only 2D right now

					// r
					pos3[1]=sqrt(ix[1][1][i]*ix[1][1][i] + ix[1][3][k]*ix[1][3][k] + ix[1][2][j]*ix[1][2][j]);
					// theta
					if(pos3[1]!=0.0) pos3[2]=acos(ix[1][2][j]/pos3[1]);
					else pos3[2]=0.0; // arbitrary;
					// phi
					if(ix[1][1][i]!=0.0) pos3[3]=atan(idx[1][3][k]/ix[1][1][i]);
					else pos3[3]=0.0;
	  
					convmat[0][0][j][i]=sin(pos3[2])*cos(pos3[3]); // dx/dr
					convmat[0][1][j][i]=cos(pos3[2])*cos(pos3[3]); // dx/dtheta
					convmat[0][2][j][i]=(-sin(pos3[3])); // dx/dphi -> 0
	  
					convmat[1][0][j][i]=sin(pos3[2])*sin(pos3[3]); // dy/dr -> 0
					convmat[1][1][j][i]=cos(pos3[2])*sin(pos3[3]); // dy/dtheta -> 0
					convmat[1][2][j][i]=cos(pos3[3]); // dy/dphi
	  
					convmat[2][0][j][i]=cos(pos3[2]); // dz/dr
					convmat[2][1][j][i]=(-sin(pos3[2])); // dz/dtheta
					convmat[2][2][j][i]=0.0; // dz/dphi
				}
				firsttime=0; // do not need to do again unless dynamic image grid(image or real)
      }
    }
#endif
    // OUTPUT INTERPOLATED DATA
    LOOPD{
      k=0; // only 2D right now


      // some vars(vectors and tensors) require transformation when converting from spherical polar coords to cart coords.
      if(SAMEGRIDD>0){
				if( (which==DTYPE)||(which==PDTYPE)||(which==ADTYPE)||(which==FLTYPE)||(which==NPTYPE)||(which==AVG2DTYPE) ){
					for(l=1;l<=REALNUMVEC;l++){

						if(which==AVG2DTYPE){
							iq1=workiqavg[6];
							iq2=workiqavg[7];
							iq3=workiqavg[8];

							iqT00=workiqavg[9];
							iqT01=workiqavg[10];
							iqT02=workiqavg[11];
							iqT11=workiqavg[12];
							iqT12=workiqavg[13];
							iqT22=workiqavg[14];
						}
						else if( (which==DTYPE)||(which==PDTYPE)||(which==ADTYPE)||(which==FLTYPE) ){
							iq1=viq[l][1];
							iq2=viq[l][2];
							iq3=viq[l][3];
						}
						else if(which==NPTYPE){
							iqT00=viq[1][1];
							iqT01=viq[1][2];
							iqT02=viq[1][3];
							iqT11=viq[2][1];
							iqT12=viq[2][2];
							iqT22=viq[2][3];
						}
						if( (which==DTYPE)||(which==PDTYPE)||(which==ADTYPE)||(which==FLTYPE)||(which==NPTYPE)||(which==AVG2DTYPE) ){
							if( (COORD==3)&&(outtype>0) ){
#if(POSTPROC==1)
		
								//vx and Bx
								ftemp[l][1]=iq1[j][i]*convmat[0][0][j][i]+iq2[j][i]*convmat[0][1][j][i]+iq3[j][i]*convmat[0][2][j][i];
								//vy and By
								ftemp[l][2]=iq1[j][i]*convmat[1][0][j][i]+iq2[j][i]*convmat[1][1][j][i]+iq3[j][i]*convmat[1][2][j][i];
								//vz and Bz
								ftemp[l][3]=iq1[j][i]*convmat[2][0][j][i]+iq2[j][i]*convmat[2][1][j][i];	      
#endif
							}
							else{
								//vx
								ftemp[l][1]=iq1[j][i];
								//vy
								ftemp[l][2]=iq2[j][i];
								//vz
								ftemp[l][3]=iq3[j][i];
							}
						}
						if( (which==AVG2DTYPE)||(which==NPTYPE) ){
							if( (COORD==3)&&(outtype>0) ){
#if(POSTPROC==1)
								// now compute tensor stuff
								// 0 0 = 
								ftempT[1]=iqT00[j][i]*convmat[0][0][j][i]*convmat[0][0][j][i]+
									iqT11[j][i]*convmat[0][1][j][i]*convmat[0][1][j][i]+
									iqT22[j][i]*convmat[0][2][j][i]*convmat[0][2][j][i]+
									2.0*iqT01[j][i]*convmat[0][0][j][i]*convmat[0][1][j][i]+
									2.0*iqT02[j][i]*convmat[0][0][j][i]*convmat[0][2][j][i]+
									2.0*iqT12[j][i]*convmat[0][1][j][i]*convmat[0][2][j][i];

								// 0 1 = 
								ftempT[2]=iqT00[j][i]*convmat[0][0][j][i]*convmat[1][0][j][i]+
									iqT11[j][i]*convmat[0][1][j][i]*convmat[1][1][j][i]+
									iqT22[j][i]*convmat[0][2][j][i]*convmat[1][2][j][i]+
									iqT01[j][i]*(convmat[0][0][j][i]*convmat[1][1][j][i]+convmat[1][0][j][i]*convmat[0][1][j][i])+
									iqT02[j][i]*(convmat[1][0][j][i]*convmat[0][2][j][i]+convmat[0][0][j][i]*convmat[1][2][j][i])+
									iqT12[j][i]*(convmat[1][1][j][i]*convmat[0][2][j][i]+convmat[0][1][j][i]*convmat[1][2][j][i]);

								// 0 2 = 
								ftempT[3]=iqT00[j][i]*convmat[0][0][j][i]*convmat[2][0][j][i]+
									iqT11[j][i]*convmat[0][1][j][i]*convmat[2][1][j][i]+
									iqT22[j][i]*convmat[0][2][j][i]*convmat[2][2][j][i]+
									iqT01[j][i]*(convmat[0][0][j][i]*convmat[2][1][j][i]+convmat[2][0][j][i]*convmat[0][1][j][i])+
									iqT02[j][i]*(convmat[2][0][j][i]*convmat[0][2][j][i]+convmat[0][0][j][i]*convmat[2][2][j][i])+
									iqT12[j][i]*(convmat[2][1][j][i]*convmat[0][2][j][i]+convmat[0][1][j][i]*convmat[2][2][j][i]);

								// 1 1 =
								ftempT[4]=iqT00[j][i]*convmat[1][0][j][i]*convmat[1][0][j][i]+
									iqT11[j][i]*convmat[1][1][j][i]*convmat[1][1][j][i]+
									iqT22[j][i]*convmat[1][2][j][i]*convmat[1][2][j][i]+
									2.0*iqT01[j][i]*convmat[1][0][j][i]*convmat[1][1][j][i]+
									2.0*iqT02[j][i]*convmat[1][0][j][i]*convmat[1][2][j][i]+
									2.0*iqT12[j][i]*convmat[1][1][j][i]*convmat[1][2][j][i];
		
								// 1 2 = 
								ftempT[5]=iqT00[j][i]*convmat[1][0][j][i]*convmat[2][0][j][i]+
									iqT11[j][i]*convmat[1][1][j][i]*convmat[2][1][j][i]+
									iqT22[j][i]*convmat[1][2][j][i]*convmat[2][2][j][i]+
									iqT01[j][i]*(convmat[1][0][j][i]*convmat[2][1][j][i]+convmat[2][0][j][i]*convmat[1][1][j][i])+
									iqT02[j][i]*(convmat[2][0][j][i]*convmat[1][2][j][i]+convmat[1][0][j][i]*convmat[2][2][j][i])+
									iqT12[j][i]*(convmat[2][1][j][i]*convmat[1][2][j][i]+convmat[1][1][j][i]*convmat[2][2][j][i]);

								// 2 2 =
								ftempT[6]=iqT00[j][i]*convmat[2][0][j][i]*convmat[2][0][j][i]+
									iqT11[j][i]*convmat[2][1][j][i]*convmat[2][1][j][i]+
									iqT22[j][i]*convmat[2][2][j][i]*convmat[2][2][j][i]+
									2.0*iqT01[j][i]*convmat[2][0][j][i]*convmat[2][1][j][i]+
									2.0*iqT02[j][i]*convmat[2][0][j][i]*convmat[2][2][j][i]+
									2.0*iqT12[j][i]*convmat[2][1][j][i]*convmat[2][2][j][i];
		
#endif
							}
							else{
								ftempT[1]=iqT00[j][i];
								ftempT[2]=iqT01[j][i];
								ftempT[3]=iqT02[j][i];
								ftempT[4]=iqT11[j][i];
								ftempT[5]=iqT12[j][i];
								ftempT[6]=iqT22[j][i];
							}
							l=REALNUMVEC; // only need to be here once
						}// end if avg2d or np
					}// end over vectors
				}// end if normal data or avg2d data
				if(which==CALCTYPE){
					for(m=1;m<=10;m+=3){ // do each component in a group
						iq1=workiqavg[m];
						iq2=workiqavg[m+1];
						iq3=workiqavg[m+2];
						if( (COORD==3)&&(outtype>0) ){
#if(POSTPROC==1)
							//vx and Bx
							ftempC[m]=iq1[j][i]*convmat[0][0][j][i]+iq2[j][i]*convmat[0][1][j][i]+iq3[j][i]*convmat[0][2][j][i];
							//vy and By
							ftempC[m+1]=iq1[j][i]*convmat[1][0][j][i]+iq2[j][i]*convmat[1][1][j][i]+iq3[j][i]*convmat[1][2][j][i];
							//vz and Bz
							ftempC[m+2]=iq1[j][i]*convmat[2][0][j][i]+iq2[j][i]*convmat[2][1][j][i];	      
#endif
						}
						else{
							//vx
							ftempC[m]=iq1[j][i];
							//vy
							ftempC[m+1]=iq2[j][i];
							//vz
							ftempC[m+2]=iq3[j][i];
						}
					}
				}
				if(which==FLINETYPE){
					for(m=1;m<=3;m+=3){ // do each component in a group
						iq3=workiqavg[1];
						if( (COORD==3)&&(outtype>0) ){
#if(POSTPROC==1)
							//vy and By
							ftempC[m+1]=iq3[j][i]*convmat[1][2][j][i];
#endif
						}
						else{
							//vz
							ftempC[m+2]=iq3[j][i];
						}
					}
				}

				// since mapped to image z vs x, map all other variables to z vs x, y, i.e. order as: (r,theta,phi=0)->(x,z,y=0)

				// actually print out the result after any transformation
				fprintf(fp," ");
				if( (which==DTYPE)||(which==PDTYPE)||(which==ADTYPE)||(which==FLTYPE) ){
					for(l=1;l<=NUMSCA;l++){
						dump1data(CENT,i,j,k,fp,iq[l][j][i]);
					}
					for(l=1;l<=REALNUMVEC;l++){
						if(which==FLTYPE){
							dump1data(CENT,i,j,k,fp,work0iq[l][j][i]);
						}
						if( (COORD==3)&&(outtype>0) ){
							// coords go to z vs x in 2D in mapping, y is flat coord
							dump1data(CENT,i,j,k,fp,ftemp[l][1]); // value interpolated and mapped, so need ftemp
							dump1data(CENT,i,j,k,fp,ftemp[l][3]); // value interpolated and mapped, so need ftemp
							dump1data(CENT,i,j,k,fp,ftemp[l][2]); // value interpolated and mapped, so need ftemp
						}
						else{
							for(m=1;m<=3;m++){
								dump1data(m,i,j,k,fp,ftemp[l][m]);
							}
						}	      
					}
				}
				else if(which==NPTYPE){
					if( (COORD==3)&&(outtype>0) ){
						// mapping of r,th,phi to x,z,y leads to this ordering change
						dump1data(CENT,i,j,k,fp,ftempT[1]);
						dump1data(CENT,i,j,k,fp,ftempT[3]);
						dump1data(CENT,i,j,k,fp,ftempT[2]);
						dump1data(CENT,i,j,k,fp,ftempT[6]);
						dump1data(CENT,i,j,k,fp,ftempT[5]);
						dump1data(CENT,i,j,k,fp,ftempT[4]);
					}
					else{
						for(l=1;l<=6;l++){ // number of components in sym tensor
							dump1data(CENT,i,j,k,fp,ftempT[l]);
						}
					}
					dump1data(CENT,i,j,k,fp,iq[1][j][i]);  // interp of nuvisc
				}// if sigma
				else if(which==AVG2DTYPE){
					for(l=1;l<=5;l++){
						dump1data(CENT,i,j,k,fp,workiqavg[l][j][i]);
					}
					if( (COORD==3)&&(outtype>0) ){
						// coords go to z vs x in 2D in mapping, y is flat coord
						dump1data(CENT,i,j,k,fp,ftemp[1][1]); // value interpolated and mapped, so need ftemp
						dump1data(CENT,i,j,k,fp,ftemp[1][3]); // value interpolated and mapped, so need ftemp
						dump1data(CENT,i,j,k,fp,ftemp[1][2]); // value interpolated and mapped, so need ftemp
					}
					else{
						for(m=1;m<=3;m++){
							dump1data(m,i,j,k,fp,ftemp[1][m]);
						}
					}
					if( (COORD==3)&&(outtype>0) ){
						// mapping of r,th,phi to x,z,y leads to this ordering change
						dump1data(CENT,i,j,k,fp,ftempT[1]);
						dump1data(CENT,i,j,k,fp,ftempT[3]);
						dump1data(CENT,i,j,k,fp,ftempT[2]);
						dump1data(CENT,i,j,k,fp,ftempT[6]);
						dump1data(CENT,i,j,k,fp,ftempT[5]);
						dump1data(CENT,i,j,k,fp,ftempT[4]);
					}
					else{
						for(l=1;l<=6;l++){ // number of components in sym tensor
							dump1data(CENT,i,j,k,fp,ftempT[l]);
						}
					}
					dump1data(CENT,i,j,k,fp,workiqavg[15][j][i]);  // interp of nuvisc

				}// end if avg2d data
				else if(which==CALCTYPE){
					for(m=1;m<=10;m+=3){
						if( (COORD==3)&&(outtype>0) ){
							// coords go to z vs x in 2D in mapping, y is flat coord
							dump1data(CENT,i,j,k,fp,ftempC[m]); // value interpolated and mapped, so need ftemp
							dump1data(CENT,i,j,k,fp,ftempC[m+2]); // value interpolated and mapped, so need ftemp
							dump1data(CENT,i,j,k,fp,ftempC[m+1]); // value interpolated and mapped, so need ftemp
						}
						else{
							dump1data(CENT,i,j,k,fp,ftempC[m]);
							dump1data(CENT,i,j,k,fp,ftempC[m+1]);
							dump1data(CENT,i,j,k,fp,ftempC[m+2]);
						}
					}
				}// end if calctype data
				else if(which==FLINETYPE){
					for(m=1;m<=3;m+=3){
						if( (COORD==3)&&(outtype>0) ){
							// SUPERMARK
							// output only needed precision to do contour plot of at most 256 levels
							// coords go to z vs x in 2D in mapping, y is flat coord
							//	      dump1data(m,i,j,k,fp,ftempC[m+1]); // value interpolated and mapped, so need ftemp
							dumpcheckfprintf(CENT,i,j,k,fp,"%.6g ",ftempC[m+1]); // value interpolated and mapped, so need ftemp
						}
						else{
							//	      dump1data(m,i,j,k,fp,ftempC[m+2]);
							dumpcheckfprintf(CENT,i,j,k,fp,"%.6g ",ftempC[m+2]);
						}
					}
				}// end if fieldlinetype data
				fprintf(fp,"\n");
      }// end if want to interp all data to same grid
      else{
				fprintf(fail_file,"Can't have SAMEGRIDD==0 as dump contains both sca and vec\n");
				myexit(1);
      }
    }//end loopd
  }// end else if doing interp
  if( ((myid==0)&&(mpicombine==1))||((mpicombine==0)&&(binaryoutput==1)) ){
    fclose(headerptr);
  }

  if(mpicombine==1){
#if(USEMPI)
		if(USEROMIO){
#if(USEROMIO)
			// now write the buffer:
			MPI_File_write_all(fh, writebuf, bufcount, MPI_FTYPE, &status);
			MPI_File_close(&fh);
			free(writebuf);
			MPI_Type_free(&newtype);
#endif
		}
		else if(USEJONIO){
			jonio_combine(1,MPI_FTYPE,numcolumns,sizeof(FTYPE),fp,jonio,writebuf,tempbuf);
			jonio_combine(2,MPI_FTYPE,numcolumns,sizeof(FTYPE),fp,jonio,writebuf,tempbuf);
			jonio_combine(3,MPI_FTYPE,numcolumns,sizeof(FTYPE),0,jonio,writebuf,tempbuf);
		}
#endif
  }
  else{
    fflush(fp); // always flush
    if(fpnull==1){
      fclose(fp); // otherwise assume caller will close it
      fp=NULL; // return as null since started null(fclose should do this)
    }
  }
}





 
#define SCAHEADER(fp,ll,outtype,iii,sliceloop,cnt,t,nx,ny) \
	  fprintf(fp,"#scalar#: %d ot: %d mt: %d sl#: %d image#: %d t=%15.10g\n",ll,outtype,iii,sliceloop,im_cnts[outtype][ll],t); \
	  fprintf(fp,"%i %i\n", nx, ny); \
	  fprintf(fp,"255\n") 

#define VECHEADER(fp,ll,outtype,iii,sliceloop,q,cnt,t,nx,ny) \
	    fprintf(fp,"#vector#: %d outtype#: %d maptype#: %d sl#: %d comp#: %d image#: %d t=%15.10g\n",ll,outtype,iii,sliceloop,q,cnt,t); \
	    fprintf(fp,"%i %i\n", nx, ny); \
	    fprintf(fp,"255\n")


// mapping from function to image space in value(color)
#if(!GAMMIEIMAGE)
#define FMAPSCA1(x) log10(x)
#define FMAPSCA2(x) log10(x)
#else
#define FMAPSCA1(x) (x)
#define FMAPSCA2(x) (x)
#endif

#define FMAPSCA3(x) (x)
#define FMAPSCAGEN(x) (x)

#define FMAPVEC1(x) (x)
#define FMAPVEC2(x) (x)

// DYNAMICMM description
// 0 : use static values as defined in initial array in image func(see analsol.c/defs.h/init.c)
// 1 : when creating image, use first min/max values as basis for all sequence of images
// 2 : when creating image, use *latest* min/max values as basis for this image
// 3 : when creating image, use latest *peak* min/max values as basis for this image, upto predetermined time: timagescale
// 4 : just read in file and use that for all time
#define ERASETILLT 0
// 1: erase latest min/max until the given time, so don't use crazy-high values that exist in early transient evolution
#define ERASETILLFACTOR 5.0 // timagescale/thisfactor is time at which we stop erasing

#define MAXSLICENUMBER 5 // maximum number of slices to output

#define SHOWBOUNDARYZONES 0
// 0: when BOUNDTYPE==2, whether to show boundary zone data on otherwise normal grid

#if(POSTPROC==0)
#define TOTALMM 1
// 1: find total time-max and time-min of all values and output to file when t=tf
#else
#define TOTALMM 0 // don't care about this in pp, and don't want to overwrite
#endif

void image(int im_cnt,int wsca,int wvec,int call_code,int outtype)
{
  static int distparuse=0;
  // CUT PASTE START 0
  char ifheader[MAXFILENAME*2];
  FTYPE liq,a,b,c,lmax,lmin,rholiq ;
  int iii,dualgo;
  int i,j,l,m,ll,mm ;
  int k;
  int p,q;
  int wshort;
  FTYPE ftempfix;
  int floop; //tells loop when first entered

  static FTYPE (*iq)[INTN2][INTN1];
  static FTYPE (*viq)[3][INTN2][INTN1];
  FILE *im_file[MAXSLICENUMBER];
  static FILE *ipartot ;
  static FILE *iparuse ;
  static char ifnam[MAXFILENAME*2],ifnamback[MAXFILENAME*2], temp[MAXFILENAME*2] ;  

  static int firstfirsttime=1;
  static int lastlasttime=0;
  static int firsttimes[ITYPES][NUMSCA+1];
  static int firsttimev[ITYPES][NUMVEC+1];
  static int nowiparuses[ITYPES][NUMSCA+1];
  static int nowiparusev[ITYPES][NUMVEC+1];
  //static int pal[3][256];
  static unsigned char pal[3][256];
  FILE*pal_file;
  FILE*size_file;
  static int outtypein=0;
  int im_cnts[ITYPES][NUMSCA+1];
  int im_cntv[ITYPES][NUMVEC+1];

  FTYPE ftemp[2];

  char temps[MAXFILENAME*2];
  char command[MAXFILENAME*2];
  static int SLICENUMBER;
  int sliceloop;

  int qstart;
  int itemp;
  static int dynamicmm3outs[ITYPES][CTYPES][NUMSCA+1][2];
  static int dynamicmm3outv[ITYPES][CTYPES][NUMVEC+1][3+1][2];
  static int createuse;
  int sj,ej,si,ei,sk,ek,iterj,iteri,iterk;
  static int PLOTAXIS;
  static int n1,n2,n3,imagen1,imagen2,imagen3;
  static int dir[3+1];
  static FTYPE Aconst,Bconst;
  static int GENERALSLICE,ko[MAXSLICENUMBER],jo[MAXSLICENUMBER],io[MAXSLICENUMBER],color,iimage,jimage;
  static int notimageoutputloop,ie[MAXSLICENUMBER],je[MAXSLICENUMBER],ke[MAXSLICENUMBER],ic[MAXSLICENUMBER],jc[MAXSLICENUMBER],kc[MAXSLICENUMBER];
  int ii,jj,sii,sjj,eii,ejj,iterjj,iterii;
  static int N[3+1],TILEIMAGE,TILEDIR,ntile1,ntile2;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
#if(USEMPI&&USEROMIO)
  MPI_Datatype newtype;
  MPI_File fh[MAXSLICENUMBER];
  MPI_Status status;
  MPI_Request request;
#endif
#if(USEMPI)
  unsigned char *jonio;
#endif
  FTYPE totalnorm,recvnorm,sendnorm;
  static int ndims;
  int array_of_gsizes[4], array_of_distribs[4];
  int order, len;
  int array_of_dargs[4], array_of_psizes[4];
  int bufcount, array_size;
  unsigned char *writebuf,*tempbuf;
  int dumi;
  FTYPE dumf;
  SFTYPE dumlf;
  char truemyidtxt[100];
  static int imagempicombine;
  int nextbuf;
  int numcolumns;


  numcolumns=1;
  //////// BEGIN
  //
  if(call_code==2) lastlasttime++;

  if(firstfirsttime==1){

    /////////////////////
    // setup geometric and MPI dependencies

    // ndims: # of dimensions of image output

    if((N1==1)||(N2==1)||(N3==1)||(POSTPROC==1)){
      TILEIMAGE=0; // no point in tiling, although it works, this confirms mpicombine method
      ndims=2; // no choice
    }
    else{
      if(numprocs==1){
				if(1){ // chooser
					TILEIMAGE=0; ndims=3 ; // choice
				}
				else{
					TILEIMAGE=1; ndims=2; // choice
				}
      }
      else{
				TILEIMAGE=0; // no choice, full 3d image
				ndims=3;
      }
    }

    ////////////////
    // for given setup, determine how to plot

    // 0: use below slicing
    // 1: tile entire data set onto 1 image of sqrt(N1*N2*N3) size
    if(TILEIMAGE==0){
      // below for coordinate line slicings
      // assume this by default until better viz in place
      
      // 3: xy-plane 2: xz-plane 1: yz-plane
      GENERALSLICE=0;
      if(N1==1){
				PLOTAXIS=1;
      }
      else if(N2==1){
				PLOTAXIS=2;
      }
      else if(N3==1){
				PLOTAXIS=3;
      }
      else{ // then 3d
				if(numprocs==1){
					if(ndims==2){
						GENERALSLICE=1; // choice: 1: if 3D case, force choice of below slices
					}
					else GENERALSLICE=0; // choice
				}
				else{ // no choice
					GENERALSLICE=0;
				}
      }

      if(GENERALSLICE==0){
				if(ndims==2){
					if(PLOTAXIS==3){
						n1=N1;
						n2=N2;
						imagen1=IMGN1;
						imagen2=IMGN2;
						imagen3=1;
						// sets up which array values are actually considered
						dir[1]=1;
						dir[2]=1;
						dir[3]=0;
					}
					else if(PLOTAXIS==2){
						n1=N1;
						n2=N3;
						imagen1=IMGN1;
						imagen2=IMGN3;
						imagen3=1;
						// sets up which array values are actually considered
						dir[1]=1;
						dir[2]=0;
						dir[3]=1;
					}
					else if(PLOTAXIS==1){
						n1=N2;
						n2=N3;
						imagen1=IMGN2;
						imagen2=IMGN3;
						imagen3=1;
						// sets up which array values are actually considered
						dir[1]=0;
						dir[2]=1;
						dir[3]=1;
					}
				}
				else{ // ndims==3
					n1=N1;
					n2=N2;
					n3=N3;
					imagen1=IMGN1;
					imagen2=IMGN2;
					imagen3=IMGN3;
					dir[1]=1;
					dir[2]=1;
					dir[3]=1;
				}
				SLICENUMBER=1; // no choice
      }
      else{ // if GENERALSLICE==1
				if(0){ // slanted slice
					SLICENUMBER=1;
					ko[0]=N3-1;
					jo[0]=0;
					io[0]=0;
	  
					ke[0]=N3-1;
					je[0]=N2-1;
					ie[0]=0;
	  
					jc[0]=N2-1;
					ic[0]=N1-1;
					// kc fixed by equations for this image surface
					itemp=copysign(1.,(FTYPE)(ke[0]-ko[0]))*(ie[0]-io[0])*(ic[0]-ie[0])-(je[0]-jo[0])*(jc[0]-je[0])/(fabs(ke[0]-ko[0])+SSMALL);
					if(!((itemp==0)&&(ke[0]-ko[0]==0))){
						kc[0]=ke[0]-copysign(1.,(FTYPE)(ke[0]-ko[0]))*(ie[0]-io[0])*(ic[0]-ie[0])-(je[0]-jo[0])*(jc[0]-je[0])/(fabs(ke[0]-ko[0])+SSMALL);
					}
					else kc[0]=0; // must specify then
				}
				else{ // equatorial and r-z slice
					SLICENUMBER=2;
	  
					/*
					// equatorial
					io[0]=0;
					jo[0]=N2/2;
					ko[0]=0;
	    
					ie[0]=N1-1;
					je[0]=N2/2;
					ke[0]=0;
	    
					ic[0]=N1-1;
					jc[0]=N2/2;
					// kc fixed by equations for this image surface
					itemp=copysign(1.,(FTYPE)(ke[0]-ko[0]))*(ie[0]-io[0])*(ic[0]-ie[0])-(je[0]-jo[0])*(jc[0]-je[0])/(fabs(ke[0]-ko[0])+SSMALL);
					if(!((itemp==0)&&(ke[0]-ko[0]==0))){
					kc[0]=ke[0]-copysign(1.,(FTYPE)(ke[0]-ko[0]))*(ie[0]-io[0])*(ic[0]-ie[0])-(je[0]-jo[0])*(jc[0]-je[0])/(fabs(ke[0]-ko[0])+SSMALL);
					}
					else kc[0]=N3-1; // must specify then
					*/
					// r-theta slice
					ko[0]=N3/8+N3/2;
					jo[0]=0;
					io[0]=0;
	  
					ke[0]=N3/8+N3/2;
					je[0]=0;
					ie[0]=N1-1;
	  
					jc[0]=N2-1;
					ic[0]=N1-1;
					// kc fixed by equations for this image surface
					itemp=copysign(1.,(FTYPE)(ke[0]-ko[0]))*(ie[0]-io[0])*(ic[0]-ie[0])-(je[0]-jo[0])*(jc[0]-je[0])/(fabs(ke[0]-ko[0])+SSMALL);
					if(!((itemp==0)&&(ke[0]-ko[0]==0))){
						kc[0]=ke[0]-copysign(1.,(FTYPE)(ke[0]-ko[0]))*(ie[0]-io[0])*(ic[0]-ie[0])-(je[0]-jo[0])*(jc[0]-je[0])/(fabs(ke[0]-ko[0])+SSMALL);
					}
					else kc[0]=N3/8+N3/2; // must specify then
	  
					// r-theta slice
					ko[1]=N3/8;
					jo[1]=0;
					io[1]=0;
	  
					ke[1]=N3/8;
					je[1]=0;
					ie[1]=N1-1;
	  
					jc[1]=N2-1;
					ic[1]=N1-1;
					// kc fixed by equations for this image surface
					itemp=copysign(1.,(FTYPE)(ke[1]-ko[1]))*(ie[1]-io[1])*(ic[1]-ie[1])-(je[1]-jo[1])*(jc[1]-je[1])/(fabs(ke[1]-ko[1])+SSMALL);
					if(!((itemp==0)&&(ke[1]-ko[1]==0))){
						kc[1]=ke[1]-copysign(1.,(FTYPE)(ke[1]-ko[1]))*(ie[1]-io[1])*(ic[1]-ie[1])-(je[1]-jo[1])*(jc[1]-je[1])/(fabs(ke[1]-ko[1])+SSMALL);
					}
					else kc[1]=N3/8; // must specify then
	  
				}// else if equatorial slice
				// since no general coordinate plane, do:
				dir[3]=1;
				dir[2]=1;
				dir[1]=1;
				n1=imagen1=NBIG;
				n2=imagen2=NBIG;
				n3=imagen3=1;
				// solution is then: k=k0+(int)(Aconst*(FTYPE)(i-i0)+Bconst*(FTYPE)(j-j0));
	
      }// else if generalslice==1
    }// end if TILEIMAGE==0
    else{ // if TILEIMAGE==1
      PLOTAXIS=0;
      GENERALSLICE=0;
      
      // doesn't work with POSTPROC==1
      SLICENUMBER=1;
      // pick the tiled direction
      if((N1>1)&&(N2>1)){
				TILEDIR=3;
				// i->ii
				// j->jj
      }
      else if((N1>1)&&(N3>1)){
				TILEDIR=2;
				// k->ii
				// i->jj
      }
      else if((N2>1)&&(N3>1)){
				TILEDIR=1;
	
      }
      N[1]=N1;
      N[2]=N2;
      N[3]=N3;
      
      dir[3]=1;
      dir[2]=1;
      dir[1]=1;
      ntile1=(int)ceil(sqrt((FTYPE)(N[TILEDIR]))-1E-6);
      ntile2=(int)floor(sqrt((FTYPE)(N[TILEDIR]))+1E-6);
      while(ntile1*ntile2<N[TILEDIR]){
				if(ntile1>ntile2) ntile2++; else ntile1++;
      }
      n1=imagen1=ntile1*N[TILEDIR%3+1];
      n2=imagen2=ntile2*N[3-(4-TILEDIR)%3];
      n3=imagen3=1;

    }


    // determine whether output file of image will be combined or not combined from all CPUs

    if(mpicombine==1){// only ref to mpicombine itself in here
      // only do this if NDIM==2 since too difficult to combine anything but a slice
      if((N3==1)&&(ncpux3==1)&&(TILEIMAGE==0)&&(GENERALSLICE==0)) imagempicombine=1;
      else{
				if(ndims==3){
					imagempicombine=1; // doing 3d image block
				}
				else imagempicombine=0;
      }
    }
    else imagempicombine=0;

    fprintf(log_file,"n1: %d n2: %d n3: %d ntile1: %d ntile2: %d imagen1: %d imagen2: %d imagen3: %d dir[1]: %d dir[2]: %d dir[3]: %d TILEDIR: %d SLICENUMBER: %d PLOTAXIS: %d GENREALSLICE: %d TILEIMAGE: %d\n",n1,n2,n3,ntile1,ntile2,imagen1,imagen2,imagen3,dir[1],dir[2],dir[3],TILEDIR,SLICENUMBER,PLOTAXIS,GENERALSLICE,TILEIMAGE); fflush(log_file);

#if(0)
    fprintf(stdout,"%15.10g %15.10g\n",mms[outtype][0][1][0],mms[outtype][0][1][1]); fflush(stdout);
#endif

  } // end if firstfirsttime
  

  // some MPI stuff

  if(imagempicombine==0){
    strcpy(truemyidtxt,myidtxt);
  }
  else strcpy(truemyidtxt,"");

  if(imagempicombine==1){
#if(USEMPI)
		if(USEROMIO){
#if(USEROMIO)
			if(ndims==3){
				//create the distributed array filetype
				order = MPI_ORDER_C;
      
				array_of_gsizes[2] = totalsize[1];
				array_of_gsizes[1] = totalsize[2];
				array_of_gsizes[0] = totalsize[3];
      
				array_of_distribs[2] = MPI_DISTRIBUTE_BLOCK;
				array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
				array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
      
				array_of_dargs[2] = MPI_DISTRIBUTE_DFLT_DARG;
				array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
				array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
      
				array_of_psizes[2]=ncpux1;
				array_of_psizes[1]=ncpux2;
				array_of_psizes[0]=ncpux3;
      
				MPI_Type_create_darray(numprocs, myid, ndims, array_of_gsizes, 
															 array_of_distribs, array_of_dargs,
															 array_of_psizes, order, MPI_BYTE, &newtype);
				MPI_Type_commit(&newtype);
				MPI_Type_size(newtype, &bufcount);
				bufcount = bufcount/sizeof(unsigned char);
				writebuf = (unsigned char *) malloc(bufcount * sizeof(unsigned char));
			}
			else{ // ndims==2
				//create the distributed array filetype
				order = MPI_ORDER_C;
      
				array_of_gsizes[1] = totalsize[1];
				array_of_gsizes[0] = totalsize[2];
      
				array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
				array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
      
				array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
				array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
      
				array_of_psizes[1]=ncpux1;
				array_of_psizes[0]=ncpux2;
      
				MPI_Type_create_darray(numprocs, myid, ndims, array_of_gsizes, 
															 array_of_distribs, array_of_dargs,
															 array_of_psizes, order, MPI_BYTE, &newtype);
				MPI_Type_commit(&newtype);
				MPI_Type_size(newtype, &bufcount);
				bufcount = bufcount/sizeof(unsigned char);
				writebuf = (unsigned char *) malloc(bufcount * sizeof(unsigned char));
			}
#endif
		}
		else if(USEJONIO){
			jonio_init_mem(numcolumns,sizeof(unsigned char),&jonio,0,0,&writebuf,0,0,&tempbuf,0,0);
		}
#endif
  }


  // short hand
#define INDEX [k*dir[3]][j*dir[2]][i*dir[1]]




  // CUT PASTE END 0

  if(outtype>(ITYPES-1)){
    fprintf(fail_file,"outtype: %d > number of types allocated\n",outtype);
    myexit(1);
  }

  // assign work space
  iq=workiq;
  viq=workviq;		     

  // CUT PASTE START 1

  // for now assume all im_cnt same
  for(i=1;i<=NUMSCA;i++){
    im_cnts[outtype][i]=im_cnt;
  }
  for(i=1;i<=NUMVEC;i++){
    im_cntv[outtype][i]=im_cnt;
  }


  // open image parameter file
  if(firstfirsttime){

    // initialize bitswitch for whether outputted iparuse info or not
    for(i=0;i<ITYPES;i++){
      for(iii=0;iii<CTYPES;iii++){
				for(ll=1;ll<=NUMSCA;ll++){
					for(j=0;j<=1;j++){
						dynamicmm3outs[i][iii][ll][j]=0;
					}
				}
      }
    }
    for(i=0;i<ITYPES;i++){
      for(iii=0;iii<CTYPES;iii++){
				for(ll=1;ll<=NUMVEC;ll++){
					for(q=0;q<=3;q++){
						for(j=0;j<=1;j++){
							dynamicmm3outv[i][iii][ll][q][j]=0;
						}
					}
				}
      }
    }

    // initialize whether first time or not to image a variable
    for(i=0;i<ITYPES;i++){
      for(j=1;j<=NUMSCA;j++){
				firsttimes[i][j]=1;
				nowiparuses[i][j]=0;
      }
      for(j=1;j<=NUMVEC;j++){
				firsttimev[i][j]=1;
				nowiparusev[i][j]=0;
      }
    }

    // create file that says image size

    if(myid<=0){ // only done on myid=0
      createuse=0;
      // setup file output
      sprintf(temps,"%s%s",DATADIR,IMAGEDIR);	
      
      if(SAMPLEI>0) strcpy(ifnam,"i");
      else strcpy(ifnam,"");
      sprintf(temps,"%s%simsize.out",temps,ifnam);
      if( (size_file=fopen(temps,"wt"))==NULL){
				fprintf(fail_file,"cannot open: %s\n",temps);
				myexit(1);
      }
      else{
				fprintf(size_file,"%5d %5d\n",imagen1,imagen2);
				fclose(size_file);
      }

      // create directory structure

      sprintf(temps,"%s%s",DATADIR,IMAGEDIR);	
  
#if(USEGM==0)
      if(!GAMMIEIMAGE){
				for(i=0;i<ITYPES;i++){ // outtype
					for(ll=1;ll<=NUMSCA;ll++){
						for(iii=0;iii<CTYPES;iii++){
							for(k=0;k<SLICENUMBER;k++){
								sprintf(ifnam,"%simx%01d-%01d-%01d-s%01d/",temps,i,iii,k,ll);
								sprintf(command,"mkdir %s",ifnam);
								system(command);
								if(deleteoldimg){
									sprintf(command,"rm %s*%s*",ifnam,IMGEXT);
									system(command);
								}
							}
						}
					}
				}
				for(i=0;i<ITYPES;i++){
					for(ll=1;ll<=REALNUMVEC;ll++){
						for(iii=0;iii<CTYPES;iii++){
							for(q=0;q<=3;q++){
								for(k=0;k<SLICENUMBER;k++){
									sprintf(ifnam,"%simx%01d-%01d-%01d-v%01d-%01d/",temps,i,iii,k,ll,q);
									sprintf(command,"mkdir %s",ifnam);
									system(command);
									if(deleteoldimg){
										sprintf(command,"rm %s*%s*",ifnam,IMGEXT);
										system(command);
									}
								}
							}
						}
					}
				}
      }
      else{
				sprintf(command,"mkdir %s%s",DATADIR,IMAGEDIR);
				system(command);
				sprintf(command,"mkdir %s%s",DATADIR,IMAGEIDIR);
				system(command);
      }
#endif

      /////////////////
      // setup file that holds image max/min data

      sprintf(temps,"%s%s",DATADIR,IMAGEDIR);
      if(TOTALMM==1){
				sprintf(ifnam,"%s0_ipartot%s",temps,PAREXT);
				if((ipartot = fopen(ifnam,"w"))==NULL) {
					fprintf(fail_file,"error opening ipartot output file %s\n",ifnam) ;
					myexit(1) ;
				}
      }

      // check if already have paruse file
      if((DYNAMICMM==3)||(DYNAMICMM==4)){
	
				sprintf(temps,"%s%s",DATADIR,IMAGEDIR);
				sprintf(ifnam,"%s0_iparuse%s",temps,PAREXT);
				sprintf(ifnamback,"%s0_iparuse%s.back",temps,PAREXT);
				if((POSTPROC==1)||((POSTPROC==0)&&(runtype>0)&&(directinput>0))){ // then should have no completed use file, no matter what time, and want to use it for seemless continuation of image if not postproc.
					// if pp, then only have file if time is later than timagescale.  rest is unreliable for computation no matter what since dynamic anyways.  assumes timagescale is set correctly........ Just assume file exists if pp for now
					distparuse=1;
					fprintf(logfull_file,"#reading iparuse file\n");
					fflush(logfull_file);
					if((iparuse = fopen(ifnam,"r"))==NULL) {
						fprintf(fail_file,"error reading iparuse output file %s\n",ifnam) ;
						myexit(1) ;
					}
	  
					qstart=0;
					for(outtypein=0;outtypein<ITYPES;outtypein++){
						for(ll=1;ll<=NUMSCA;ll++){
							for(i=0;i<=1;i++){ // min and max
								for(m=0;m<CTYPES;m++){
		  
									if(outtypein==0){ // only outputted this, so only read in this
										while(fgetc(iparuse)!='=');
										fscanf(iparuse,INPUTPAR,&mms[outtypein][m][ll][i]);
									}
									else mms[outtypein][m][ll][i]=mms[0][m][ll][i];
									if( (t>=timagescale)||(POSTPROC==1)) dynamicmm3outs[outtypein][m][ll][i]=1; // tells already outputted
								}
							}
						}
						for(ll=1;ll<=REALNUMVEC;ll++){
							for(q=qstart;q<=3;q++){
								for(i=0;i<=1;i++){ // min and max
									for(m=0;m<CTYPES;m++){
										if((outtypein==0)&&(q>=1)){ // only outputted this, so only read in this
											while(fgetc(iparuse)!='=');
											fscanf(iparuse,INPUTPAR,&mmv[outtypein][m][ll][q][i]);
										}
										else{
											if(q!=0) mmv[outtypein][m][ll][q][i]=mmv[0][m][ll][q][i];
										}
										if((t>=timagescale)||(POSTPROC==1)) dynamicmm3outv[outtypein][m][ll][q][i]=1;
									}
								}
							}
							// q==0
							for(i=0;i<=1;i++){ // min and max
								for(m=0;m<CTYPES;m++){
									if(i==1){ // certainly max, but should really set or read in q=0
										mmv[outtypein][m][ll][0][i]=mmv[0][m][ll][1][i]*mmv[0][m][ll][1][i]+mmv[0][m][ll][2][i]*mmv[0][m][ll][2][i]+mmv[0][m][ll][3][i]*mmv[0][m][ll][3][i];
									}
									if(i==0){// certainly min
										mmv[outtypein][m][ll][0][i]=0;
									}
								}
							}
						}
					}
					fclose(iparuse);
					fprintf(logfull_file,"#done reading iparuse file\n");
					fflush(logfull_file);
				} // end if read paruse file
				else distparuse=0;
        // then no usefile or done reading in data, create it if need to
				if( ((t<timagescale)||((distparuse==0)&&(t>timagescale)))  &&(POSTPROC==0)){
					createuse=1;
#if(USEGM==0)
					sprintf(temps,"cp %s %s",ifnam,ifnamback);
					system(temps);
#else
					rename(ifnam,ifnamback);
#endif
				}
				else{
					createuse=0;
#if(USEGM==0)
					sprintf(temps,"cp %s %s",ifnam,ifnamback);
					system(temps);
#endif
					iparuse=NULL;
				}
				// should write immediately
				if( (((distparuse==0)&&(t>timagescale)))  &&(POSTPROC==0) ){
					for(i=0;i<ITYPES;i++){
						for(j=1;j<=NUMSCA;j++){
							nowiparuses[i][j]=1;
						}
						for(j=1;j<=NUMVEC;j++){
							nowiparusev[i][j]=1;
						}
					}
				}
      }// end if dynamicmm==3


      /////////////////////
      //  open palette

      if(IMAGEFORMAT>0){
				// load pallette
				sprintf(temps,"%s%s",DATADIR,IMAGEDIR);
				sprintf(ifnam,"%sgenimages%s",temps,".pal");
				if( (pal_file = fopen(ifnam,"rb"))==NULL){
					fprintf(fail_file,"error opening %s\n",ifnam) ;
					myexit(2) ;
				}
				for(j=0;j<3;j++){
					for(i=0;i<256;i++){
						fread(&pal[j][i],sizeof(unsigned char),1,pal_file); 
						//fprintf(stderr,"%d %d %d\n",j,i,pal[j][i]);
					}
				}
				fclose(pal_file);
      }
    }// end if myid<=0


    ////////////
    // distribute min/max image data to other cpus
    if(numprocs>1){
#if(USEMPI)
      MPI_Bcast(&distparuse,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
    }
  
    if((distparuse==1)&&(numprocs>1)){
      outtypein=0;
      qstart=1;
#if(USEMPI)
      for(ll=1;ll<=3;ll++){
				for(i=0;i<=1;i++){ // min and max
					for(m=0;m<CTYPES;m++){
						MPI_Bcast(&mms[outtypein][m][ll][i],1,MPI_FTYPE,0,MPI_COMM_WORLD);
					}
				}
      }
      for(ll=1;ll<=REALNUMVEC;ll++){
				for(q=qstart;q<=3;q++){
					for(i=0;i<=1;i++){ // min and max
						for(m=0;m<CTYPES;m++){
							MPI_Bcast(&mmv[outtypein][m][ll][q][i],1,MPI_FTYPE,0,MPI_COMM_WORLD);
						}
					}
				}
      }
#endif
    }

    ///////////////////
    // distribute pallette to other cpus
    if(numprocs>1){
			// need to transfer pallette to other cpus(better than having each read file)
#if(USEMPI)
      MPI_Bcast(&pal[0][0],3*256,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
#endif
    }

    // for images (sets up outer values of interpolated data
    for(i=0;i<ITYPES;i++){ // both views
      for(j=0;j<CTYPES;j++){ // both comps
				outerdefs[i][j][1]=mms[i][j][1][0]; // rho
				outerdefs[i][j][2]=mms[i][j][2][0]; // en
				outerdefs[i][j][3]=mms[i][j][3][0]; // pot
	
				outerdefv[i][j][1][0]=mmv[i][j][1][0][0]; // magnitude of v
				outerdefv[i][j][1][1]=mmv[i][j][1][1][0]; // v1
				outerdefv[i][j][1][2]=mmv[i][j][1][2][0]; // v2
				outerdefv[i][j][1][3]=mmv[i][j][1][3][0]; // v3
	
				outerdefv[i][j][2][0]=mmv[i][j][2][0][0];
				outerdefv[i][j][2][1]=mmv[i][j][2][1][0];
				outerdefv[i][j][2][2]=mmv[i][j][2][2][0];
				outerdefv[i][j][2][3]=mmv[i][j][2][3][0];
      }
    }                 

  }// end if firstfirsttime
  // CUT PASTE END 1
  



  //////////////////////////////////////////
  //////////////////////////////////////////
  ///////////  SCALARS
  if(wsca!=0){
    for(l=1;l<=NUMSCA;l++){
      /* if not to do all, pick */
      if(wsca!=-1){ ll=wsca; }
      else ll=l;

      if(SAMPLEI>0){
				interpolate(0,ll,SAMPLEI,SAMEGRIDI,n1,n2,imagen1,imagen2,s[ll][0],iq[ll],outerdefs[0][0][ll],outtype);
      }

      //////////////
      // write current min/max to temp values
      if( ((DYNAMICMM==1)&&(firsttimes[outtype][ll]==1) )||(TOTALMM==1)||(DYNAMICMM==2)||( (DYNAMICMM==3)&&(t<timagescale) ) ){
				if(SAMPLEI>0){
					if( (firsttimes[outtype][ll]==1)||(DYNAMICMM==2)||(ERASETILLT&&(t<timagescale/ERASETILLFACTOR) ) ){
						for(m=0;m<CTYPES;m++){
							mmst[outtype][m][ll][0]=iq[ll][IMGN2-1][0];
							mmst[outtype][m][ll][1]=iq[ll][IMGN2-1][0];
						}
					}
#include "imageloophead.h"
					notimageoutputloop=1;
					while(1){
#include "imageloopinside.h"
						for(m=0;m<CTYPES;m++){
							if(mmst[outtype][m][ll][0]>iq[ll][j][i]) mmst[outtype][m][ll][0]=iq[ll][j][i];
							if(mmst[outtype][m][ll][1]<iq[ll][j][i]) mmst[outtype][m][ll][1]=iq[ll][j][i];
						}
					}
        }
        else{
					if( (firsttimes[outtype][ll]==1)||(DYNAMICMM==2)||(ERASETILLT&&(t<timagescale/ERASETILLFACTOR) ) ){
						for(m=0;m<CTYPES;m++){
							mmst[outtype][m][ll][0]=s[ll][0][IMGN2-1][0]; // not right for slices
							mmst[outtype][m][ll][1]=s[ll][0][IMGN2-1][0]; // not right for slices
						}
					}
#include "imageloophead.h"
					notimageoutputloop=1;
					while(1){
#include "imageloopinside.h"
						for(m=0;m<CTYPES;m++){
							if(mmst[outtype][m][ll][0]>s[ll]INDEX) mmst[outtype][m][ll][0]=s[ll]INDEX;
							if(mmst[outtype][m][ll][1]<s[ll]INDEX) mmst[outtype][m][ll][1]=s[ll]INDEX;
						}
					}
				}
				/////////////////
				// set used values from current/temp values
				if(DYNAMICMM>0){
					if( (DYNAMICMM==2)||( (firsttimes[outtype][ll]==1)&&(DYNAMICMM==1) )||( (DYNAMICMM==3)&&(t<timagescale) ) ){
						for(m=0;m<CTYPES;m++){
							mms[outtype][m][ll][0]=mmst[outtype][m][ll][0];
							mms[outtype][m][ll][1]=mmst[outtype][m][ll][1];
						}
						if(numprocs>1){
#if(USEMPI)
							for(i=0;i<=1;i++){ // min and max
								for(m=0;m<CTYPES;m++){
									if(i==0){
										MPI_Allreduce(&mmst[outtype][m][ll][i],&mms[outtype][m][ll][i],1,MPI_FTYPE,MPI_MIN,MPI_COMM_WORLD);
									}
									else if(i==1){
										MPI_Allreduce(&mmst[outtype][m][ll][i],&mms[outtype][m][ll][i],1,MPI_FTYPE,MPI_MAX,MPI_COMM_WORLD);			
									}
								}
							}
#endif
						}
					}
				}
				////////////////
				// output min/max data to file
				if(myid<=0){
					if(TOTALMM==1){
						if((call_code==2)&&(lastlasttime==1)){
							for(i=0;i<=1;i++){ // min-max
								for(m=0;m<CTYPES;m++){// compute types
									fprintf(ipartot,"mms[%d][%d][%d][%d]= %21.15g ;\n",outtype,m,ll,i,mmst[outtype][m][ll][i]);
								}
							}
						}
					}
					if(DYNAMICMM==3){
						if(((call_code==2)&&(lastlasttime==1)) || (t>=timagescale)||(nowiparuses[outtype][ll]==1) ){
							if(createuse){
								if((iparuse = fopen(ifnam,"w"))==NULL) {
									fprintf(fail_file,"error opening iparuse output file %s\n",ifnam) ;
									myexit(1) ;
								}
								createuse=0;
							}
							for(i=0;i<=1;i++){ // min and max
								for(m=0;m<CTYPES;m++){ // compute types
									if(dynamicmm3outs[outtype][m][ll][i]==0){
										fprintf(iparuse,"mms[%d][%d][%d][%d]= %21.15g ;\n",outtype,m,ll,i,mms[outtype][m][ll][i]);
										dynamicmm3outs[outtype][m][ll][i]=1;
									}
								}
							}
							nowiparuses[outtype][ll]=0;
						}
					}
				}
				firsttimes[outtype][ll]=0; // done first time with this scalar
      }
      /*
				if(POSTPROC==1){
        if(ll<=2) dualgo=1;
        else dualgo=0;
				}
				else dualgo=0; // no need for dualgo==1 if not pp (not true if want good info for linear AND log data)
      */
      dualgo=1; // force so can use linear for linear and log for log, and linear+log to recapture data stream at both levels

      for(iii=0;iii<=dualgo;iii++){

				// create all slices file pointers
				for(sliceloop=0;sliceloop<SLICENUMBER;sliceloop++){

					///////////////////////////
					// setup and open image file
					if(!GAMMIEIMAGE){
						sprintf(temps,"%s%s",DATADIR,IMAGEDIR);
						if(SAMPLEI>0) strcpy(ifheader,"imx"); // was iimx, no need for difference since interp is pp only and only then is outtype>0
						else strcpy(ifheader,"imx");
						sprintf(temps,"%simx%01d-%01d-%01d-s%01d/",temps,outtype,iii,sliceloop,ll);
						sprintf(ifnam,"%s%s%01d-%01d-%01d-s%01d-%04d%s%s",temps,ifheader,outtype,iii,sliceloop,ll,im_cnts[outtype][ll],IMGEXT,truemyidtxt);
					}
					else{
						sprintf(ifnam,"%s%sim%01dp%04d",DATADIR,IMAGEIDIR,ll-1,im_cnts[outtype][ll]);
						sprintf(temp,"gzip > %s.gz",ifnam);
					}
					if(imagempicombine==0){
						if(IMAGEFORMAT==0){
							strcat(ifnam,".r8");
							// MARK
							//fprintf(stderr,"got here: sliceloop: %d\n",sliceloop); fflush(stderr);
							if(GZIPIMAGE!=3) im_file[sliceloop] = fopen(ifnam,"wb");
							else{
								sprintf(temp,"gzip > %s.gz",ifnam);
								strcpy(ifnam,temp); // for below fprintf
								im_file[sliceloop] = popen(ifnam,"w");
							}
							if(im_file[sliceloop]==NULL){
								fprintf(fail_file,"error opening image file: %s\n",ifnam) ;
								myexit(2) ;
							}
							if(ndims==2){
								fprintf(im_file[sliceloop],"RAW\n");	     
								SCAHEADER(im_file[sliceloop],ll,outtype,iii,sliceloop,im_cnts[outtype][ll],t,imagen1,imagen2);
							}// otherwise pure binary
						}
	    
						if(IMAGEFORMAT==1){
							strcat(ifnam,".ppm");
							if(GZIPIMAGE!=3) im_file[sliceloop] = fopen(ifnam,"wt");
							else{
								sprintf(temp,"gzip > %s.gz",ifnam);
								strcpy(ifnam,temp); // for below fprintf
								im_file[sliceloop] = popen(ifnam,"w");
							}
							if(im_file[sliceloop]==NULL){
								fprintf(fail_file,"error opening image file: %s\n",ifnam) ;
								myexit(2) ;
							}
							if(ndims==2){
								fprintf(im_file[sliceloop],"P6\n");
								SCAHEADER(im_file[sliceloop],ll,outtype,iii,sliceloop,im_cnts[outtype][ll],t,imagen1,imagen2);
							}
							if(GZIPIMAGE!=3){
								fclose(im_file[sliceloop]);
								//reopen in binary append mode
								if( (im_file[sliceloop] = fopen(ifnam, "ab" ))==NULL){
									fprintf(fail_file,"error opening image file binary append: %s\n",ifnam) ;
									myexit(2) ;
								}
							}
						}
					}// end if imagempicombine==0
					else{
						// setup file handler
						strcat(ifnam,".r8");
#if(USEMPI)
						if(USEROMIO){
#if(USEROMIO)
							MPI_File_open(MPI_COMM_WORLD, ifnam, MPI_MODE_CREATE | MPI_MODE_RDWR,MPI_INFO_NULL, &fh[sliceloop]);
							MPI_File_set_view(fh[sliceloop], 0, MPI_BYTE, newtype, "native", MPI_INFO_NULL);
#endif
						}
						else if(USEJONIO){
							jonio_init_fp(&im_file[sliceloop],1,ifnam);
						}
#endif
						// all that needs to be done now is fill writebuf with the data
	    
					}// end if imagempicombine==1
				}// end over create files for all slices

				////////////////////////
				// now set image map using min/max
				if(ll==1){
					if(iii==0){
						//fprintf(stdout,"dodo: %d %d %d %15.10g %15.10g\n",outtype,iii,ll,mms[outtype][iii][ll][0],mms[outtype][iii][ll][1]);
						b=FMAPSCA1(mms[outtype][iii][ll][0]);
						if(fabs(b-FMAPSCA1(mms[outtype][iii][ll][1]))>1E-10){
							a=255./(FMAPSCA1(mms[outtype][iii][ll][1])-b);
						}
						else{
							a=0; // since then undefined
						}
					}
					if(iii==1){
						b=FMAPSCAGEN(mms[outtype][iii][ll][0]);
						if(fabs(b-FMAPSCAGEN(mms[outtype][iii][ll][1]))>1E-10){
							a=255./(FMAPSCAGEN(mms[outtype][iii][ll][1])-b);
						}
						else{
							a=0; // since then undefined
						}
					}
				}
				else if(ll==2){
					if(iii==0){
						b=FMAPSCA2(mms[outtype][iii][ll][0]);
						if(fabs(b-FMAPSCA2(mms[outtype][iii][ll][1]))>1E-10){
							a=255./(FMAPSCA2(mms[outtype][iii][ll][1])-b);
						}
						else{
							a=0; // since then undefined
						}
					}
					if(iii==1){
						b=FMAPSCAGEN(mms[outtype][iii][ll][0]);
						if(fabs(b-FMAPSCAGEN(mms[outtype][iii][ll][1]))>1E-10){
							a=255./(FMAPSCAGEN(mms[outtype][iii][ll][1])-b);
						}
						else{
							a=0; // since then undefined
						}
					}
				}
				else if(ll==3){
					b=FMAPSCA3(mms[outtype][iii][ll][0]);
					if(fabs(b-FMAPSCA3(mms[outtype][iii][ll][1]))>1E-10){
						a=255./(FMAPSCA3(mms[outtype][iii][ll][1])-b);
					}
					else{
						a=0; // since then undefined
					}
				}
				c=-a*b;

				//////////////////////////
				// output image to file


      	
#include "imageloophead.h"
				notimageoutputloop=0;
				if(imagempicombine==1) nextbuf=0;
				while(1){
#include "imageloopinside.h"
	  
					//fprintf(stderr,"%d %d %d %d\n",ll,k,j,i); fflush(stderr);

					// choose sample
					if(SAMPLEI>0) liq=iq[ll][j][i];
					else liq=s[ll][k][j][i];

					// scale the sample
					if(ll==1){
						if(iii==0){
							ftempfix=FMAPSCA1(liq);
							//	      fprintf(stdout,"got here: %15.10g %15.10g\n",liq,FMAPSCA1(liq),ftempfix); fflush(stdout);
						}
						if(iii==1){
							ftempfix=FMAPSCAGEN(liq);
						}
					}
					else if(ll==2){
						if(iii==0){
							ftempfix=FMAPSCA2(liq);
						}
						if(iii==1){
							ftempfix=FMAPSCAGEN(liq);
						}
					}
					else if(ll==3){
						ftempfix=FMAPSCA3(liq);
					}

					//fprintf(stdout,"%d %d %15.10g %15.10g %15.10g %15.10g %15.10g\n",j,i,a,b,c,ftempfix,liq);

					// fix the sample
					if(ftempfix>=b){
						liq = a*ftempfix+c ;
					}
					else liq=0.0;
	  
					if(liq > 255.) liq = 255. ;
					if(liq < 0.) liq = 0. ;

					if(!color) liq=0;


#if(0)
					fprintf(stdout,"%d %d %15.10g %15.10g %15.10g %15.10g %15.10g %d %15.10g %15.10g %d\n",i,j,mms[outtype][0][1][0],mms[outtype][0][1][1],a,b,c,(int)liq,ftempfix,iq[ll][j][i],SAMPLEI); fflush(stdout);
#endif


					// output sample
					if(imagempicombine==0){
						if(IMAGEFORMAT==0){
							//fprintf(stderr,"%d : %d %d %d %d : %d %d %d\n",(int)liq,sliceloop,ii,jj,(int)im_file[sliceloop],k,j,i); fflush(stderr);
							//	      fprintf(stderr,"a %d : %d %d : %d %d %d\n",(int)liq,sliceloop,(int)im_file[sliceloop],k,j,i); fflush(stderr);
							fputc( (int)liq , im_file[sliceloop]);   /* write value */
							//	      fprintf(stderr,"b %d : %d %d : %d %d %d\n",(int)liq,sliceloop,(int)im_file[sliceloop],k,j,i); fflush(stderr);
						}
						if(IMAGEFORMAT==1){
							fputc((int)pal[0][(int)liq] , im_file[sliceloop]);   /* write red */
							fputc((int)pal[1][(int)liq] , im_file[sliceloop]);   /* write green */
							fputc((int)pal[2][(int)liq] , im_file[sliceloop]);   /* write blue */
						}
					}
					else{
						writebuf[nextbuf++]=(unsigned char)liq;
					}	  
				}// end loop

				// close file
				if(imagempicombine==0){
					for(sliceloop=0;sliceloop<SLICENUMBER;sliceloop++){
						if(GZIPIMAGE==0){
							fclose(im_file[sliceloop]) ;
						}
						if(GZIPIMAGE==1){
							fclose(im_file[sliceloop]) ;
							strcpy(temp,"gzip ");
							strcat(temp,ifnam);
							system(temp);
						}
						else if(GZIPIMAGE==2){
							fclose(im_file[sliceloop]) ;
							mysys("gzip",ifnam);
						}
						else if(GZIPIMAGE==3){
							pclose(im_file[sliceloop]);
						}
						if(POSTPROC==1){
							printf("s");
							fflush(stdout);
						}
					}
				}
				else{
#if(USEMPI)
					for(sliceloop=0;sliceloop<SLICENUMBER;sliceloop++){
						if(USEROMIO){
#if(USEROMIO)
							// now write the buffer:
							MPI_File_write_all(fh[sliceloop], writebuf, bufcount, MPI_BYTE, &status);
							MPI_File_close(&fh[sliceloop]);
#endif
						}
						else if(USEJONIO){
							jonio_combine(1,MPI_BYTE,numcolumns,sizeof(unsigned char),im_file[sliceloop],jonio,writebuf,tempbuf);
							jonio_combine(2,MPI_BYTE,numcolumns,sizeof(unsigned char),im_file[sliceloop],jonio,writebuf,tempbuf);
						}
					}
#endif
				}
      }// over dualgo
      /* cut short loop if only to do one */
      if(wsca!=-1) l=NUMSCA;
      if(GAMMIEIMAGE&&(ll==2)) break;
    }
  }
  // fflush(stdout);
  //   if(t>1.0)  myexit(0);

  

  ///////////////////////////
  //////////// VECTORS
  if(wvec!=0){
    if(POSTPROC==0){ // no need for 0 components when not pp
      qstart=1;
    }
    else{
      if(!GAMMIEIMAGE) qstart=0;
      else qstart=1;
    }

    for(l=1;l<=REALNUMVEC;l++){
      /* if not to do all, pick */
      if(wvec!=-1){ ll=wvec; }
      else ll=l;
      
      if(SAMPLEI>0){
				for(q=1;q<=3;q++){ // should really interp q=0 if want right looking
					// determine coordinate position
					if(q==1) itemp=1;
					else if(q==2) itemp=2;
					else if(q==3) itemp=0;
					interpolate(itemp,ll,SAMPLEI,SAMEGRIDI,n1,n2,imagen1,imagen2,v[ll][q][0],viq[ll][q],outerdefv[0][0][ll][q],outtype);
					//	  interpolate(itemp,ll,SAMPLEI,SAMEGRIDI,n1,n2,imagen1,imagen2,v[ll][q][0],viq[ll][q],0.0,outtype);
				}
      }

      if( ((DYNAMICMM==1)&&(firsttimev[outtype][ll]==1) )||(TOTALMM==1)||(DYNAMICMM==2)||( (DYNAMICMM==3)&&(t<timagescale) ) ){      
				for(q=qstart;q<=3;q++){ // 0 is mag
					if(SAMPLEI>0){
						if(q>0){
							if((firsttimev[outtype][ll]==1)||(DYNAMICMM==2)||(ERASETILLT&&(t<timagescale/ERASETILLFACTOR)) )  floop=1;
#include "imageloophead.h"
							notimageoutputloop=1;
							while(1){
#include "imageloopinside.h"
								for(m=0;m<CTYPES;m++){
									if( (mmvt[outtype][m][ll][q][0]>viq[ll][q][j][i])||(floop==1)) mmvt[outtype][m][ll][q][0]=viq[ll][q][j][i];
									if( (mmvt[outtype][m][ll][q][1]<viq[ll][q][j][i])||(floop==1)) mmvt[outtype][m][ll][q][1]=viq[ll][q][j][i];
								}
								floop=0;
							}
						}
						else{
							if((firsttimev[outtype][ll]==1)||(DYNAMICMM==2)||(ERASETILLT&&(t<timagescale/ERASETILLFACTOR)) )  floop=1;
#include "imageloophead.h"
							notimageoutputloop=1;
							while(1){
#include "imageloopinside.h"
								ftemp[0]=0.0;
								for(m=1;m<=3;m++){
									ftemp[0]+=viq[ll][m][j][i]*viq[ll][m][j][i];
								}
								//		ftemp[0]=sqrt(ftemp[0]);
								for(m=0;m<CTYPES;m++){
									if( (mmvt[outtype][m][ll][q][0]>ftemp[0])||(floop==1)) mmvt[outtype][m][ll][q][0]=ftemp[0];
									if( (mmvt[outtype][m][ll][q][1]<ftemp[0])||(floop==1)) mmvt[outtype][m][ll][q][1]=ftemp[0];
								}
								floop=0;
							}
						}
					} // end if interpolating
					else{
						if(q>0){
							if((firsttimev[outtype][ll]==1)||(DYNAMICMM==2)||(ERASETILLT&&(t<timagescale/ERASETILLFACTOR)) )  floop=1;
#include "imageloophead.h"
							notimageoutputloop=1;
							while(1){
#include "imageloopinside.h"
								for(m=0;m<CTYPES;m++){
									if( (mmvt[outtype][m][ll][q][0]>v[ll][q]INDEX)||(floop==1)) mmvt[outtype][m][ll][q][0]=v[ll][q]INDEX;
									if( (mmvt[outtype][m][ll][q][1]<v[ll][q]INDEX)||(floop==1)) mmvt[outtype][m][ll][q][1]=v[ll][q]INDEX;
								}
								floop=0;
							}
						}
						else{
							if((firsttimev[outtype][ll]==1)||(DYNAMICMM==2)||(ERASETILLT&&(t<timagescale/ERASETILLFACTOR)) )  floop=1;
#include "imageloophead.h"
							notimageoutputloop=1;
							while(1){
#include "imageloopinside.h"
								ftemp[0]=0.0;
								for(m=1;m<=3;m++){
									ftemp[0]+=v[ll][m]INDEX*v[ll][m]INDEX;
								}
								//		ftemp[0]=sqrt(ftemp[0]);
								for(m=0;m<CTYPES;m++){
									if( (mmvt[outtype][m][ll][q][0]>ftemp[0])||(floop==1)) mmvt[outtype][m][ll][q][0]=ftemp[0];
									if( (mmvt[outtype][m][ll][q][1]<ftemp[0])||(floop==1)) mmvt[outtype][m][ll][q][1]=ftemp[0];
								}
								floop=0;
							}
						}
					}// end if non-interpolating
				}// over mag+dir

	
				if(DYNAMICMM>0){
					// mmv is what's actually used.
					if( (DYNAMICMM==2)||( (firsttimev[outtype][ll]==1)&&(DYNAMICMM==1) )||( (DYNAMICMM==3)&&(t<timagescale) ) ){
						for(q=qstart;q<=3;q++){
							for(m=0;m<CTYPES;m++){
								mmv[outtype][m][ll][q][0]=mmvt[outtype][m][ll][q][0];
								mmv[outtype][m][ll][q][1]=mmvt[outtype][m][ll][q][1];
							}
						}
						if(numprocs>1){
#if(USEMPI)
							for(q=qstart;q<=3;q++){
								for(i=0;i<=1;i++){ // min and max
									for(m=0;m<CTYPES;m++){
										if(i==0){
											MPI_Allreduce(&mmvt[outtype][m][ll][q][i],&mmv[outtype][m][ll][q][i],1,MPI_FTYPE,MPI_MIN,MPI_COMM_WORLD);
										}
										else if(i==1){
											MPI_Allreduce(&mmvt[outtype][m][ll][q][i],&mmv[outtype][m][ll][q][i],1,MPI_FTYPE,MPI_MAX,MPI_COMM_WORLD);
										}
									}
								}
							}
#endif
						}
					}
				}

				if(myid<=0){
					// need mmvt diff from mmv so can find all-time min/max here
					if(TOTALMM==1){
						if((call_code==2)&&(lastlasttime==1)){
							for(q=qstart;q<=3;q++){ // component
								for(i=0;i<=1;i++){ // min and max
									for(m=0;m<CTYPES;m++){ // compute types
										fprintf(ipartot,"mmv[%d][%d][%d][%d][%d]= %21.15g ;\n",outtype,m,ll,q,i,mmvt[outtype][m][ll][q][i]);
									}
								}
							}
						}
					}
					if(DYNAMICMM==3){
						if(((call_code==2)&&(lastlasttime==1)) || (t>=timagescale)||(nowiparusev[outtype][ll]==1) ){
							if(createuse){
								if((iparuse = fopen(ifnam,"w"))==NULL) {
									fprintf(fail_file,"error opening iparuse output file %s\n",ifnam) ;
									myexit(1) ;
								}
								createuse=0;
							}
							for(q=qstart;q<=3;q++){ // component
								for(i=0;i<=1;i++){ // min and max
									for(m=0;m<CTYPES;m++){ // compute types
										if(dynamicmm3outv[outtype][m][ll][q][i]==0){
											fprintf(iparuse,"mmv[%d][%d][%d][%d][%d]= %21.15g ;\n",outtype,m,ll,q,i,mmv[outtype][m][ll][q][i]);
											dynamicmm3outv[outtype][m][ll][q][i]=1;
										}
									}
								}
							}
							nowiparusev[outtype][ll]=0;
						}
					}
				}
				firsttimev[outtype][ll]=0;
      }// if firsttime[outtype][ll]==1 or dynamicmm==2 or totalmm==1
      
      
      for(q=qstart;q<=3;q++){ // over components





				if(POSTPROC==1){
					dualgo=0; // for now till really need and can fix problems with ==1 part
					//dualgo=1;
				}
				else dualgo=0;	 // no need for dualgo==1 if not pp.  Definitely, just use linear scalar and linear vel to compute this in pp.

				for(iii=0;iii<=dualgo;iii++){

					for(sliceloop=0;sliceloop<SLICENUMBER;sliceloop++){
						// setup file output
						if(!GAMMIEIMAGE){
							sprintf(temps,"%s%s",DATADIR,IMAGEDIR);	    
	      
							if(SAMPLEI>0) strcpy(ifheader,"imx");
							else strcpy(ifheader,"imx");
							sprintf(temps,"%simx%01d-%01d-%01d-v%01d-%01d/",temps,outtype,iii,sliceloop,ll,q);
							sprintf(ifnam,"%s%s%01d-%01d-%01d-v%01d-%01d-%04d%s%s",temps,ifheader,outtype,iii,sliceloop,ll,q,im_cnts[outtype][ll],IMGEXT,truemyidtxt);
						}
						else{
							sprintf(ifnam,"%s%sim%01dp%04d",DATADIR,IMAGEIDIR,(ll-1)*3+(q-1)+2,im_cnts[outtype][ll]);// gets VX,VY,VZ,BX,BY,BZ
						}
	    
						if(imagempicombine==0){
							if(IMAGEFORMAT==0){
								strcat(ifnam,".r8");
								if(GZIPIMAGE!=3) im_file[sliceloop] = fopen(ifnam,"wb");
								else{
									sprintf(temp,"gzip > %s.gz",ifnam);
									strcpy(ifnam,temp); // for below fprintf
									im_file[sliceloop] = popen(ifnam,"w");
								}
								if(im_file[sliceloop]==NULL){
									fprintf(fail_file,"error opening image file: %s\n",ifnam) ;
									myexit(2) ;
								}
								if(ndims==2){
									fprintf(im_file[sliceloop],"RAW\n");
									VECHEADER(im_file[sliceloop],ll,outtype,iii,sliceloop,q,im_cnts[outtype][ll],t,imagen1,imagen2);
								}
							}
							if(IMAGEFORMAT==1){
								strcat(ifnam,".ppm");
								if(GZIPIMAGE!=3) im_file[sliceloop] = fopen(ifnam,"wt");
								else{
									sprintf(temp,"gzip > %s.gz",ifnam);
									strcpy(ifnam,temp); // for below fprintf
									im_file[sliceloop] = popen(ifnam,"w");
								}
								if(im_file[sliceloop]==NULL){
									fprintf(fail_file,"error opening image file: %s\n",ifnam) ;
									myexit(2) ;
								}
								if(ndims==2){
									fprintf(im_file[sliceloop],"P6\n");
									VECHEADER(im_file[sliceloop],ll,outtype,iii,sliceloop,q,im_cnts[outtype][ll],t,imagen1,imagen2);
								}
								if(GZIPIMAGE!=3){
									fclose(im_file[sliceloop]);
									//reopen in binary append mode
									if( (im_file[sliceloop] = fopen(ifnam, "ab" ))==NULL){
										fprintf(fail_file,"error opening image file binary append: %s\n",ifnam) ;
										myexit(2) ;
									}
								}
							}
						}
						else{
							// setup file handler
							strcat(ifnam,".r8");
#if(USEMPI)
							if(USEROMIO){
#if(USEROMIO)
								MPI_File_open(MPI_COMM_WORLD, ifnam, MPI_MODE_CREATE | MPI_MODE_RDWR, 
															MPI_INFO_NULL, &fh[sliceloop]);
								MPI_File_set_view(fh[sliceloop], 0, MPI_BYTE, newtype, "native", MPI_INFO_NULL);
								// all that needs to be done now is fill writebuf with the data
#endif
							}
							else if(USEJONIO){
								jonio_init_fp(&im_file[sliceloop],1,ifnam); 
							}
#endif
						}
					}
					// setup map
					if(ll==1){
						if(iii==0){
							// warning:  silent fail if 0,1 are both 0 for q=0 ?
							b=FMAPVEC1(mmv[outtype][iii][ll][q][0]);
							if(fabs(b-FMAPVEC1(mmv[outtype][iii][ll][q][1]))>1E-10){
								a=255./(FMAPVEC1(mmv[outtype][iii][ll][q][1])-b);
							}
							else{
								a=0;
							}
						}
						if(iii==1){ // rho*v // only good if v range is - to + values
							b=FMAPVEC1(mms[outtype][iii][1][1]*mmv[outtype][iii][ll][q][0]);
							if(fabs(b-FMAPVEC1(mms[outtype][iii][1][1]*mmv[outtype][iii][ll][q][1]))>1E-10){
								a=255./(FMAPVEC1(mms[outtype][iii][1][1]*mmv[outtype][iii][ll][q][1])-b);
							}
							else{
								a=0;
							}
						}
					}
					if(ll==2){
						if(iii==0){// B
							b=FMAPVEC2(mmv[outtype][iii][ll][q][0]);
							if(fabs(b-FMAPVEC2(mmv[outtype][iii][ll][q][1]))>1E-10){
								a=255./(FMAPVEC2(mmv[outtype][iii][ll][q][1])-b);
							}
							else{
								a=0;
							}
						}
						if(iii==1){ // rho*B // only good if v range is - to + values
							b=FMAPVEC2(mms[outtype][iii][1][1]*mmv[outtype][iii][ll][q][0]);
							if(fabs(b-FMAPVEC2(mms[outtype][iii][1][1]*mmv[outtype][iii][ll][q][1]))>1E-10){
								a=255./(FMAPVEC2(mms[outtype][iii][1][1]*mmv[outtype][iii][ll][q][1])-b);
							}
							else{
								a=0;
							}
						}
					}
					c=-a*b;


#include "imageloophead.h"
					notimageoutputloop=0;
					if(imagempicombine==1) nextbuf=0;
					while(1){
#include "imageloopinside.h"
						liq=0;
						if(q==0){
							for(m=1;m<=3;m++){
								if(SAMPLEI>0){
									liq =viq[ll][m][j][i]*viq[ll][m][j][i]+liq ;
								}
								else{
									liq =v[ll][m]INDEX*v[ll][m]INDEX+liq ;
								}
							}
							if(SAMPLEI>0){
								rholiq=iq[1][j][i];
							}
							else{
								rholiq=s[1]INDEX;
							}
							//	      liq=sqrt(liq);
						}
						else{
							if(SAMPLEI>0){
								liq =viq[ll][q][j][i] ;
								rholiq=iq[1][j][i];
							}
							else{
								liq =v[ll][q]INDEX ;
								rholiq=s[1]INDEX;
							}
						}
						if(iii==1){
							liq=rholiq*liq;
						}


						if(ll==1){
							ftempfix=FMAPVEC1(liq);
						}
						else if(ll==2){
							ftempfix=FMAPVEC2(liq);
						}
						if(ftempfix>=b){
							liq = a*ftempfix+c ;
						}
						else liq=0.0;
	    
						if(liq > 255.) liq = 255. ;
						if(liq < 0.) liq = 0. ;


						if(!color) liq=0;

						if(imagempicombine==0){
	      
							if(IMAGEFORMAT==0){
								fputc( (int)liq , im_file[sliceloop]);   /* write value */
							}
							if(IMAGEFORMAT==1){
								fputc((int)pal[0][(int)liq] , im_file[sliceloop]);   /* write red */
								fputc((int)pal[1][(int)liq] , im_file[sliceloop]);   /* write green */
								fputc((int)pal[2][(int)liq] , im_file[sliceloop]);   /* write blue */
							}
						}
						else{
							//	      writebuf[j*N1+i]=(int)liq;
							writebuf[nextbuf++]=(unsigned char)liq;
						}
					} // end loop

					if(imagempicombine==0){
						for(sliceloop=0;sliceloop<SLICENUMBER;sliceloop++){
							if(GZIPIMAGE==0){
								fclose(im_file[sliceloop]) ;
							}
							if(GZIPIMAGE==1){
								fclose(im_file[sliceloop]) ;
								strcpy(temp,"gzip ");
								strcat(temp,ifnam);
								system(temp);
							}
							else if(GZIPIMAGE==2){
								fclose(im_file[sliceloop]) ;
								mysys("gzip",ifnam);
							}
							else if(GZIPIMAGE==3){
								pclose(im_file[sliceloop]);
							}
							if(POSTPROC==1){
								printf("v");
								fflush(stdout);
							}
						}
					}
					else{
#if(USEMPI)
						for(sliceloop=0;sliceloop<SLICENUMBER;sliceloop++){
							if(USEROMIO){
#if(USEROMIO)
								// now write the buffer:
								MPI_File_write_all(fh[sliceloop], writebuf, bufcount, MPI_BYTE, &status);
								MPI_File_close(&fh[sliceloop]);
#endif
							}
							else if(USEJONIO){
								jonio_combine(1,MPI_BYTE,numcolumns,sizeof(unsigned char),im_file[sliceloop],jonio,writebuf,tempbuf);
								jonio_combine(2,MPI_BYTE,numcolumns,sizeof(unsigned char),im_file[sliceloop],jonio,writebuf,tempbuf);
	  
							}
						}
#endif
					}
				}//dualgo
      } // over components of the vector
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMVEC;
    }// over vectors
  }// if any vectors
  firstfirsttime=0;

  if(imagempicombine==1){
#if(USEMPI)
		if(USEROMIO){
#if(USEROMIO)
			free(writebuf);
			MPI_Type_free(&newtype);
#endif
		}
		else if(USEJONIO){ // only 1 call to stage==3
			jonio_combine(3,MPI_BYTE,numcolumns,sizeof(unsigned char),im_file[sliceloop],jonio,writebuf,tempbuf);
		}
#endif
  }
  if(myid<=0){
    if((call_code==2)&&(lastlasttime==1)){
      // should be done now
      if(TOTALMM==1){ if(ipartot!=NULL) fclose(ipartot); } // assumes 1 outtype per run
      if(DYNAMICMM==3){ if(iparuse!=NULL) fclose(iparuse); } // assumes 1 outtype per run
    }
  }
  
  if(call_code==2) lastlasttime++;
  if(POSTPROC==1){
		printf("\n");
    fflush(stdout);
  }
}


FTYPE ffv_calcs(int type, int term, int i, int j, int k) // forces,fluxes, and volume terms
{
  SFTYPE tempxx,tempxy,tempxz;
  SFTYPE tempyx,tempyy,tempyz;
  SFTYPE tempzx,tempzy,tempzz;
  SFTYPE temp1,temp2,plus,minus;
  SFTYPE vxa,vya,vza,bxa,bya,bza;
  SFTYPE ethi,eki,egi,ebi;
  SFTYPE posx,posy,posz;
  SFTYPE rho;
  

#if(DOINVSOL2FUN)
  invsol2=invsol2fun[k][j][i];
#endif

  
  // now begins the forces (really accelerations)
  if(type==1){
    if(term==1){
      // x1 component of -(v.grad)v in cart or spherical coords
      tempxx=(v[1][1][k][j][ip1]-v[1][1][k][j][i])*OARC21(k,j,i) ; // x-component of gradient of x-component of v
      tempxy=(v[1][1][k][jp1][i]-v[1][1][k][j][i])*OARC22(k,j,i) ; // y-component of gradient of x-component of v
      tempxz=(v[1][1][kp1][j][i]-v[1][1][k][j][i])*OARC23(k,j,i)*ODX(2,3,k) ; // z-component of gradient of x-component of v
      temp1=e2z_1(v[1][1],k,j,i)*tempxx+e2z_2(v[1][2],k,j,i)*tempxy+e2z_3(v[1][3],k,j,i)*tempxz;
      if(COORD==3){
				temp2=-e2z_2(v[1][2],k,j,i)*e2z_2(v[1][2],k,j,i)/G2(2,i)-e2z_3(v[1][3],k,j,i)*e2z_3(v[1][3],k,j,i)/G3(2,i);
      }
      else temp2=0;
      return(-(temp1+temp2));
    }
    else if(term==2){
      // x2 component of -(v.grad)v in cart or spherical coords
      tempyx=(v[1][2][k][j][ip1]-v[1][2][k][j][i])*OARC31(k,j,i) ; // x-component of gradient of y-component of v
      tempyy=(v[1][2][k][jp1][i]-v[1][2][k][j][i])*OARC32(k,j,i) ; // y-component of gradient of y-component of v
      tempyz=(v[1][2][kp1][j][i]-v[1][2][k][j][i])*OARC33(k,j,i)*ODX(2,3,k) ; // z-component of gradient of y-component of v
      temp1=e2z_1(v[1][1],k,j,i)*tempyx+e2z_2(v[1][2],k,j,i)*tempyy+e2z_3(v[1][3],k,j,i)*tempyz;
      if(COORD==3){
				temp2=e2z_2(v[1][2],k,j,i)*e2z_1(v[1][1],k,j,i)/G2(2,i)-e2z_3(v[1][3],k,j,i)*e2z_3(v[1][3],k,j,i)*DG4(2,j)/(G2(2,i)*G4(2,j));
      }
      else temp2=0;
      return(-(temp1+temp2));
    }
    else if(term==3){
      // x3 component of -(v.grad)v in cart or spherical coords
      tempzx=(v[1][3][k][j][ip1]-v[1][3][k][j][i])*OARC11(k,j,i) ; // x-component of gradient of z-component of v
      tempzy=(v[1][3][k][jp1][i]-v[1][3][k][j][i])*OARC12(k,j,i) ; // y-component of gradient of z-component of v
      tempzz=(v[1][3][kp1][j][i]-v[1][3][k][j][i])*OARC13(k,j,i)*ODX(1,3,k) ; // z-component of gradient of z-component of v
      temp1=e2z_1(v[1][1],k,j,i)*tempzx+e2z_2(v[1][2],k,j,i)*tempzy+e2z_3(v[1][3],k,j,i)*tempzz;
      if(COORD==3){
				temp2=e2z_3(v[1][3],k,j,i)*e2z_1(v[1][1],k,j,i)/G3(2,i)+e2z_3(v[1][3],k,j,i)*e2z_2(v[1][2],k,j,i)*DG4(2,j)/(G2(2,i)*G4(2,j));
      }
      else temp2=0;
      return(-(temp1+temp2));      
    }
    else if(term==4){
			// x1 component of -grad(P)/rho
			if(wgam){
				plus=(gam-1.0)*z2e_1(s[2],k,j,ip1mac(i));
				minus=(gam-1.0)*z2e_1(s[2],k,j,i);
      }
      else{
				plus=cs*cs*z2e_1(s[1],k,j,ip1mac(i));
				minus=cs*cs*z2e_1(s[1],k,j,i);
      }
      temp1=(plus-minus)*OARC21(k,j,i);
			if(CSLIMIT==0) rho=s[1][k][j][i];
			else rho=s[1][k][j][i]+gam*(gam-1.0)*s[2][k][j][i]*invsol2;
      return(-temp1/rho);
    }
    else if(term==5){
      // x2 component of -grad(P)/rho
      if(wgam){
				plus=(gam-1.0)*z2e_2(s[2],k,jp1mac(j),i);
				minus=(gam-1.0)*z2e_2(s[2],k,j,i);
      }
      else{
				plus=cs*cs*z2e_2(s[1],k,jp1mac(j),i);
				minus=cs*cs*z2e_2(s[1],k,j,i);
      }
      temp1=(plus-minus)*OARC32(k,j,i);
			if(CSLIMIT==0) rho=s[1][k][j][i];
			else rho=s[1][k][j][i]+gam*(gam-1.0)*s[2][k][j][i]*invsol2;
      return(-temp1/rho);
    }
    else if(term==6){
      // x3 component of -grad(P)/rho
      if(wgam){
				plus=(gam-1.0)*z2e_3(s[2],kp1mac(k),j,i);
				minus=(gam-1.0)*z2e_3(s[2],k,j,i);
      }
      else{
				plus=cs*cs*z2e_3(s[1],kp1mac(k),j,i);
				minus=cs*cs*z2e_3(s[1],k,j,i);
      }
      temp1=(plus-minus)*OARC13(k,j,i)*ODX(1,3,k);
			if(CSLIMIT==0) rho=s[1][k][j][i];
			else rho=s[1][k][j][i]+gam*(gam-1.0)*s[2][k][j][i]*invsol2;
      return(-temp1/rho);      
    }
    else if(term==7){
      // x1 component of -grad(phi)
      plus=z2e_1(s[3],k,j,ip1mac(i));
      minus=z2e_1(s[3],k,j,i);
      temp1=(plus-minus)*OARC21(k,j,i);
      return(-temp1);
    }
    else if(term==8){
      // x2 component of -grad(phi)
      plus=z2e_2(s[3],k,jp1mac(j),i);
      minus=z2e_2(s[3],k,j,i);
      temp1=(plus-minus)*OARC32(k,j,i);
      return(-temp1);
    }
    else if(term==9){
			// x3 component of -grad(phi)
			plus=z2e_3(s[3],kp1mac(k),j,i);
      minus=z2e_3(s[3],k,j,i);
      temp1=(plus-minus)*OARC13(k,j,i)*ODX(1,3,k);
      return(-temp1);
    }
    else if(term==10){
			// x1 component of (curl(b)xb)/rho
			temp1=(curlcv21(v[2],k,j,i)+curlcv22(v[2],k,j,i))*e2z_3(v[2][3],k,j,i);
			temp2=-(curlcv31(v[2],k,j,i)+curlcv32(v[2],k,j,i))*e2z_2(v[2][2],k,j,i);
			if(ALFVENLIMIT==0) rho=s[1][k][j][i];
			else{
				bxa = e2z_1(v[2][1],k,j,i);
				bya = e2z_2(v[2][2],k,j,i);
				bza = e2z_3(v[2][3],k,j,i);    				
				rho=s[1][k][j][i]+(bxa*bxa + bya*bya + bza*bza)*invsol2;
			}
			return((temp1+temp2)/rho);
    }
    else if(term==11){
      // x2 component of (curl(b)xb)/rho
			temp1=-(curlcv11(v[2],k,j,i)+curlcv12(v[2],k,j,i))*e2z_3(v[2][3],k,j,i);
			temp2=(curlcv31(v[2],k,j,i)+curlcv32(v[2],k,j,i))*e2z_1(v[2][1],k,j,i);
			if(ALFVENLIMIT==0) rho=s[1][k][j][i];
			else{
				bxa = e2z_1(v[2][1],k,j,i);
				bya = e2z_2(v[2][2],k,j,i);
				bza = e2z_3(v[2][3],k,j,i);    				
				rho=s[1][k][j][i]+(bxa*bxa + bya*bya + bza*bza)*invsol2;
			}
			return((temp1+temp2)/rho);
    }
    else if(term==12){
      // x3 component of (curl(b)xb)/rho
			temp1=(curlcv11(v[2],k,j,i)+curlcv12(v[2],k,j,i))*e2z_2(v[2][2],k,j,i);
			temp2=-(curlcv21(v[2],k,j,i)+curlcv22(v[2],k,j,i))*e2z_1(v[2][1],k,j,i);
			if(ALFVENLIMIT==0) rho=s[1][k][j][i];
			else{
				bxa = e2z_1(v[2][1],k,j,i);
				bya = e2z_2(v[2][2],k,j,i);
				bza = e2z_3(v[2][3],k,j,i);    				
				rho=s[1][k][j][i]+(bxa*bxa + bya*bya + bza*bza)*invsol2;
			}
			return((temp1+temp2)/rho);
    }
    else if(term==13){
      // x1 component of (-grad(b^2/2))/rho
			bxa = v[2][1][k][j][ip1];
			bya = v2tov1(v[2][2],k,j,ip1mac(i));
			bza = v3tov1(v[2][3],k,j,ip1mac(i));			
			plus = 0.5*(bxa*bxa + bya*bya + bza*bza) ;
			bxa = v[2][1][k][j][i];
			bya = v2tov1(v[2][2],k,j,i);
			bza = v3tov1(v[2][3],k,j,i);			
			minus = 0.5*(bxa*bxa + bya*bya + bza*bza) ;
      temp1=(plus-minus)*OARC21(k,j,i);
			if(ALFVENLIMIT==0) rho=s[1][k][j][i];
			else{
				bxa = e2z_1(v[2][1],k,j,i);
				bya = e2z_2(v[2][2],k,j,i);
				bza = e2z_3(v[2][3],k,j,i);    				
				rho=s[1][k][j][i]+(bxa*bxa + bya*bya + bza*bza)*invsol2;
			}
      return(-temp1/rho);
    }
    else if(term==14){
      // x2 component of (-grad(b^2/2))/rho
			bxa = v1tov2(v[2][1],k,jp1mac(j),i);
			bya = v[2][2][k][jp1][i];
			bza = v3tov2(v[2][3],k,jp1mac(j),i);			
			plus = 0.5*(bxa*bxa + bya*bya + bza*bza) ;
			bxa = v1tov2(v[2][1],k,j,i);
			bya = v[2][2][k][j][i];
			bza = v3tov2(v[2][3],k,j,i);			
			minus = 0.5*(bxa*bxa + bya*bya + bza*bza) ;
      temp1=(plus-minus)*OARC32(k,j,i);
			if(ALFVENLIMIT==0) rho=s[1][k][j][i];
			else{
				bxa = e2z_1(v[2][1],k,j,i);
				bya = e2z_2(v[2][2],k,j,i);
				bza = e2z_3(v[2][3],k,j,i);    				
				rho=s[1][k][j][i]+(bxa*bxa + bya*bya + bza*bza)*invsol2;
			}
      return(-temp1/rho);
    }
    else if(term==15){
      // x3 component of (-grad(b^2/2))/rho
			bxa = v1tov3(v[2][1],kp1mac(k),j,i);
			bya = v2tov3(v[2][2],kp1mac(k),j,i);
			bza = v[2][3][kp1][j][i];
			plus = 0.5*(bxa*bxa + bya*bya + bza*bza) ;
			bxa = v1tov3(v[2][1],k,j,i);
			bya = v2tov3(v[2][2],k,j,i);
			bza = v[2][3][k][j][i];
			minus = 0.5*(bxa*bxa + bya*bya + bza*bza) ;
      temp1=(plus-minus)*OARC13(k,j,i)*ODX(1,3,k);
			if(ALFVENLIMIT==0) rho=s[1][k][j][i];
			else{
				bxa = e2z_1(v[2][1],k,j,i);
				bya = e2z_2(v[2][2],k,j,i);
				bza = e2z_3(v[2][3],k,j,i);    				
				rho=s[1][k][j][i]+(bxa*bxa + bya*bya + bza*bza)*invsol2;
			}
      return(-temp1/rho);
    }
    else if(term==16){
      // x1 component of (b.grad)b/rho in cart or spherical coords
			// (just copy of v version above with v[1]->v[2]
      tempxx=(v[2][1][k][j][ip1]-v[2][1][k][j][i])*OARC21(k,j,i) ; // x-component of gradient of x-component of v
      tempxy=(v[2][1][k][jp1][i]-v[2][1][k][j][i])*OARC22(k,j,i) ; // y-component of gradient of x-component of v
      tempxz=(v[2][1][kp1][j][i]-v[2][1][k][j][i])*OARC23(k,j,i)*ODX(2,3,k) ; // z-component of gradient of x-component of v
      temp1=e2z_1(v[2][1],k,j,i)*tempxx+e2z_2(v[2][2],k,j,i)*tempxy+e2z_3(v[2][3],k,j,i)*tempxz;
      if(COORD==3){
				temp2=-e2z_2(v[2][2],k,j,i)*e2z_2(v[2][2],k,j,i)/G2(2,i)-e2z_3(v[2][3],k,j,i)*e2z_3(v[2][3],k,j,i)/G3(2,i);
      }
      else temp2=0;
			if(ALFVENLIMIT==0) rho=s[1][k][j][i];
			else{
				bxa = e2z_1(v[2][1],k,j,i);
				bya = e2z_2(v[2][2],k,j,i);
				bza = e2z_3(v[2][3],k,j,i);    				
				rho=s[1][k][j][i]+(bxa*bxa + bya*bya + bza*bza)*invsol2;
			}
      return((temp1+temp2)/rho);      
    }
    else if(term==17){
      // x2 component of (b.grad)b/rho in cart or spherical coords
      tempyx=(v[2][2][k][j][ip1]-v[2][2][k][j][i])*OARC31(k,j,i) ; // x-component of gradient of y-component of v
      tempyy=(v[2][2][k][jp1][i]-v[2][2][k][j][i])*OARC32(k,j,i) ; // y-component of gradient of y-component of v
      tempyz=(v[2][2][kp1][j][i]-v[2][2][k][j][i])*OARC33(k,j,i)*ODX(2,3,k) ; // z-component of gradient of y-component of v
      temp1=e2z_1(v[2][1],k,j,i)*tempyx+e2z_2(v[2][2],k,j,i)*tempyy+e2z_3(v[2][3],k,j,i)*tempyz;
      if(COORD==3){
				temp2=e2z_2(v[2][2],k,j,i)*e2z_1(v[2][1],k,j,i)/G2(2,i)-e2z_3(v[2][3],k,j,i)*e2z_3(v[2][3],k,j,i)*DG4(2,j)/(G2(2,i)*G4(2,j));
      }
      else temp2=0;
			if(ALFVENLIMIT==0) rho=s[1][k][j][i];
			else{
				bxa = e2z_1(v[2][1],k,j,i);
				bya = e2z_2(v[2][2],k,j,i);
				bza = e2z_3(v[2][3],k,j,i);    				
				rho=s[1][k][j][i]+(bxa*bxa + bya*bya + bza*bza)*invsol2;
			}
      return((temp1+temp2)/rho);      
    }
    else if(term==18){
      // x3 component of (b.grad)b/rho in cart or spherical coords
      tempzx=(v[2][3][k][j][ip1]-v[2][3][k][j][i])*OARC11(k,j,i) ; // x-component of gradient of y-component of v
      tempzy=(v[2][3][k][jp1][i]-v[2][3][k][j][i])*OARC12(k,j,i) ; // y-component of gradient of y-component of v
      tempzz=(v[2][3][kp1][j][i]-v[2][3][k][j][i])*OARC13(k,j,i)*ODX(1,3,k) ; // z-component of gradient of y-component of v
      temp1=e2z_1(v[2][1],k,j,i)*tempzx+e2z_2(v[2][2],k,j,i)*tempzy+e2z_3(v[2][3],k,j,i)*tempzz;
      if(COORD==3){
				temp2=e2z_3(v[2][3],k,j,i)*e2z_1(v[2][1],k,j,i)/G3(2,i)+e2z_3(v[2][3],k,j,i)*e2z_2(v[2][2],k,j,i)*DG4(2,j)/(G2(2,i)*G4(2,j));
      }
      else temp2=0;
			if(ALFVENLIMIT==0) rho=s[1][k][j][i];
			else{
				bxa = e2z_1(v[2][1],k,j,i);
				bya = e2z_2(v[2][2],k,j,i);
				bza = e2z_3(v[2][3],k,j,i);    				
				rho=s[1][k][j][i]+(bxa*bxa + bya*bya + bza*bza)*invsol2;
			}
      return((temp1+temp2)/rho);      
    }
  }
  // now begins the fluxes
  else if(type==2){
    if(term==1){
      // energy flux of kinetic energy (x1 component)
			vxa = e2z_1(v[1][1],k,j,i);      
			vya = e2z_2(v[1][2],k,j,i);
			vza = e2z_3(v[1][3],k,j,i) ;
			eki = 0.5*s[1][k][j][i]*(vxa*vxa + vya*vya + vza*vza) ;
			temp1=eki*e2z_1(v[1][1],k,j,i);
			return(temp1*G2(2,i)*G3(2,i)*DVL(1,2,j)*dx[1][3][k]);
    }
    else if(term==2){
      // energy flux of kinetic energy (x2 component)
			vxa = e2z_1(v[1][1],k,j,i);      
			vya = e2z_2(v[1][2],k,j,i);
			vza = e2z_3(v[1][3],k,j,i) ;
			eki = 0.5*s[1][k][j][i]*(vxa*vxa + vya*vya + vza*vza) ;
			temp1=eki*e2z_2(v[1][2],k,j,i);
			return(temp1*G3(2,i)*G4(2,j)*dx[1][3][k]*dx[1][1][i]);
    }
    else if(term==3){
      // energy flux of kinetic energy (x3 component)
			vxa = e2z_1(v[1][1],k,j,i);      
			vya = e2z_2(v[1][2],k,j,i);
			vza = e2z_3(v[1][3],k,j,i) ;
			eki = 0.5*s[1][k][j][i]*(vxa*vxa + vya*vya + vza*vza) ;
			temp1=eki*e2z_3(v[1][3],k,j,i);
			return(temp1*G2(2,i)*dx[1][2][j]*dx[1][1][i]);
    }
    else if(term==4){
			// note this is not quite right for CSLIMIT==1, all other terms types correct.

      // energy flux of enthalpy energy (x1 component)
			/* Equation of state : enthalpy per unit mass times density */
			if(press==1){
				if(wgam) ethi = gam*s[2][k][j][i] ;
				else ethi = gam*cs*cs*s[1][k][j][i]*log(s[1][k][j][i]) ; // ? GODMARK
			}
			else ethi=0;
			temp1=ethi*e2z_1(v[1][1],k,j,i);
			return(temp1*G2(2,i)*G3(2,i)*DVL(1,2,j)*dx[1][3][k]);
    }
    else if(term==5){
      // energy flux of enthalpy energy (x2 component)
			/* Equation of state : enthalpy per unit mass times density */
			if(press==1){
				if(wgam) ethi = gam*s[2][k][j][i] ;
				else ethi = gam*cs*cs*s[1][k][j][i]*log(s[1][k][j][i]) ; // ? GODMARK
			}
			else ethi=0;
			temp1=ethi*e2z_2(v[1][2],k,j,i);
			return(temp1*G3(2,i)*G4(2,j)*dx[1][3][k]*dx[1][1][i]);
    }
    else if(term==6){
      // energy flux of enthalpy energy (x3 component)
			/* Equation of state : enthalpy per unit mass times density */
			if(press==1){
				if(wgam) ethi = gam*s[2][k][j][i] ;
				else ethi = gam*cs*cs*s[1][k][j][i]*log(s[1][k][j][i]) ; // ? GODMARK
			}
			else ethi=0;

			temp1=ethi*e2z_3(v[1][3],k,j,i);
			return(temp1*G2(2,i)*dx[1][2][j]*dx[1][1][i]);
    }
    else if(term==7){
      // energy flux of potential energy (x1 component)
			//egi = 0.5*s[1][k][j][i]*s[3][k][j][i] ; // for self-gravity
			egi = s[1][k][j][i]*s[3][k][j][i] ;			
			temp1=egi*e2z_1(v[1][1],k,j,i);
			return(temp1*G2(2,i)*G3(2,i)*DVL(1,2,j)*dx[1][3][k]);
    }
    else if(term==8){
      // energy flux of potential energy (x2 component)
			//egi = 0.5*s[1][k][j][i]*s[3][k][j][i] ; // for self-gravity
			egi = s[1][k][j][i]*s[3][k][j][i] ;			
			temp1=egi*e2z_2(v[1][2],k,j,i);
			return(temp1*G3(2,i)*G4(2,j)*dx[1][3][k]*dx[1][1][i]);
    }
    else if(term==9){
      // energy flux of potential energy (x3 component)
			//egi = 0.5*s[1][k][j][i]*s[3][k][j][i] ; // for self-gravity
			egi = s[1][k][j][i]*s[3][k][j][i] ;			
			temp1=egi*e2z_3(v[1][3],k,j,i);
			return(temp1*G2(2,i)*dx[1][2][j]*dx[1][1][i]);
    }
    else if(term==10){
      // energy flux of magnetic tension energy (x1 component)
			ebi = (e2z_1(v[1][1],k,j,i)*e2z_1(v[2][1],k,j,i)+e2z_2(v[1][2],k,j,i)*e2z_2(v[2][2],k,j,i)+e2z_3(v[1][3],k,j,i)*e2z_3(v[2][3],k,j,i));
			temp1=-ebi*e2z_1(v[2][1],k,j,i);
			return(temp1*G2(2,i)*G3(2,i)*DVL(1,2,j)*dx[1][3][k]);
    }
    else if(term==11){
      // energy flux of magnetic tension energy (x2 component)
			ebi = (e2z_1(v[1][1],k,j,i)*e2z_1(v[2][1],k,j,i)+e2z_2(v[1][2],k,j,i)*e2z_2(v[2][2],k,j,i)+e2z_3(v[1][3],k,j,i)*e2z_3(v[2][3],k,j,i));

			temp1=-ebi*e2z_2(v[2][2],k,j,i);
			return(temp1*G3(2,i)*G4(2,j)*dx[1][3][k]*dx[1][1][i]);
    }
    else if(term==12){
      // energy flux of magnetic tension energy (x3 component)
			ebi = (e2z_1(v[1][1],k,j,i)*e2z_1(v[2][1],k,j,i)+e2z_2(v[1][2],k,j,i)*e2z_2(v[2][2],k,j,i)+e2z_3(v[1][3],k,j,i)*e2z_3(v[2][3],k,j,i));

			temp1=-ebi*e2z_3(v[2][3],k,j,i);
			return(temp1*G2(2,i)*dx[1][2][j]*dx[1][1][i]);
    }
    else if(term==13){
      // energy flux of magnetic pressure energy (x1 component)
			bxa = e2z_1(v[2][1],k,j,i);
			bya = e2z_2(v[2][2],k,j,i);
			bza = e2z_3(v[2][3],k,j,i);    
			ebi = (bxa*bxa + bya*bya + bza*bza) ;

			temp1=ebi*e2z_1(v[1][1],k,j,i);			
			return(temp1*G2(2,i)*G3(2,i)*DVL(1,2,j)*dx[1][3][k]);
    }
    else if(term==14){
      // energy flux of magnetic pressure energy (x2 component)
			bxa = e2z_1(v[2][1],k,j,i);
			bya = e2z_2(v[2][2],k,j,i);
			bza = e2z_3(v[2][3],k,j,i);    
			ebi = (bxa*bxa + bya*bya + bza*bza) ;

			temp1=ebi*e2z_2(v[1][2],k,j,i);
			return(temp1*G3(2,i)*G4(2,j)*dx[1][3][k]*dx[1][1][i]);
    }
    else if(term==15){
      // energy flux of magnetic pressure energy (x3 component)
			bxa = e2z_1(v[2][1],k,j,i);
			bya = e2z_2(v[2][2],k,j,i);
			bza = e2z_3(v[2][3],k,j,i);    
			ebi = (bxa*bxa + bya*bya + bza*bza) ;

			temp1=ebi*e2z_3(v[1][3],k,j,i);
			return(temp1*G2(2,i)*dx[1][2][j]*dx[1][1][i]);
    }
    else if(term==16){
      // flux of mass (x1 component)
			return(s[1][k][j][i]*e2z_1(v[1][1],k,j,i)*G2(2,i)*G3(2,i)*DVL(1,2,j)*dx[1][3][k]);
    }
    else if(term==17){
      // flux of mass (x2 component)
			return(s[1][k][j][i]*e2z_1(v[1][1],k,j,i)*G3(2,i)*G4(2,j)*dx[1][3][k]*dx[1][1][i]);
    }
    else if(term==18){
      // flux of mass (x3 component)
			return(s[1][k][j][i]*e2z_1(v[1][1],k,j,i)*G2(2,i)*dx[1][2][j]*dx[1][1][i]);
    }
    else if(term==19){ // totals for now, not broken up into terms
      // flux of x1 momentum (x1 component)
			bxa = e2z_1(v[2][1],k,j,i);
			bya = e2z_2(v[2][2],k,j,i);
			bza = e2z_3(v[2][3],k,j,i);    
			ebi = (bxa*bxa + bya*bya + bza*bza) ;
			if(press==1){
				if(wgam) temp1 = (gam-1.0)*s[2][k][j][i] ;
				else temp1 = cs*cs*s[1][k][j][i];
			}
			else temp1=0;

			return(s[1][k][j][i]*e2z_1(v[1][1],k,j,i)*e2z_1(v[1][1],k,j,i)+(temp1+s[1][k][j][i]*s[3][k][j][i]+0.5*ebi)-e2z_1(v[2][1],k,j,i)*e2z_1(v[2][1],k,j,i)*G2(2,i)*G3(2,i)*DVL(1,2,j)*dx[1][3][k]);
    }
    else if(term==20){
      // flux of x1 momentum (x2 component)
			return(s[1][k][j][i]*e2z_1(v[1][1],k,j,i)*e2z_2(v[1][2],k,j,i)-e2z_1(v[2][1],k,j,i)*e2z_2(v[2][2],k,j,i)*G3(2,i)*G4(2,j)*dx[1][3][k]*dx[1][1][i]);
    }
    else if(term==21){
      // flux of x1 momentum (x3 component)
			return(s[1][k][j][i]*e2z_1(v[1][1],k,j,i)*e2z_3(v[1][3],k,j,i)-e2z_1(v[2][1],k,j,i)*e2z_3(v[2][3],k,j,i)*G2(2,i)*dx[1][2][j]*dx[1][1][i]);
    }
    else if(term==22){
      // flux of x2 momentum (x1 component)
			return(s[1][k][j][i]*e2z_2(v[1][2],k,j,i)*e2z_1(v[1][1],k,j,i)-e2z_2(v[2][2],k,j,i)*e2z_1(v[2][1],k,j,i)*G2(2,i)*G3(2,i)*DVL(1,2,j)*dx[1][3][k]);
    }
    else if(term==23){
      // flux of x2 momentum (x2 component)
			bxa = e2z_1(v[2][1],k,j,i);
			bya = e2z_2(v[2][2],k,j,i);
			bza = e2z_3(v[2][3],k,j,i);    
			ebi = (bxa*bxa + bya*bya + bza*bza) ;
			if(press==1){
				if(wgam) temp1 = (gam-1.0)*s[2][k][j][i] ;
				else temp1 = cs*cs*s[1][k][j][i];
			}
			else temp1=0;

			return(s[1][k][j][i]*e2z_2(v[1][2],k,j,i)*e2z_2(v[1][2],k,j,i)+(temp1+s[1][k][j][i]*s[3][k][j][i]+0.5*ebi)-e2z_2(v[2][2],k,j,i)*e2z_2(v[2][2],k,j,i)*G3(2,i)*G4(2,j)*dx[1][3][k]*dx[1][1][i]);
    }
    else if(term==24){
      // flux of x2 momentum (x3 component)
			return(s[1][k][j][i]*e2z_2(v[1][2],k,j,i)*e2z_3(v[1][3],k,j,i)-e2z_2(v[2][2],k,j,i)*e2z_3(v[2][3],k,j,i)*G2(2,i)*dx[1][2][j]*dx[1][1][i]);
    }
    else if(term==25){
      // flux of x3 momentum (x1 component)
			return(s[1][k][j][i]*e2z_3(v[1][3],k,j,i)*e2z_1(v[1][1],k,j,i)-e2z_3(v[2][3],k,j,i)*e2z_1(v[2][1],k,j,i)*G2(2,i)*G3(2,i)*DVL(1,2,j)*dx[1][3][k]);
    }
    else if(term==26){
      // flux of x3 momentum (x2 component)
			return(s[1][k][j][i]*e2z_3(v[1][3],k,j,i)*e2z_2(v[1][2],k,j,i)-e2z_3(v[2][3],k,j,i)*e2z_2(v[2][2],k,j,i)*G3(2,i)*G4(2,j)*dx[1][3][k]*dx[1][1][i]);
    }
    else if(term==27){
      // flux of x3 momentum (x3 component)
			bxa = e2z_1(v[2][1],k,j,i);
			bya = e2z_2(v[2][2],k,j,i);
			bza = e2z_3(v[2][3],k,j,i);    
			ebi = (bxa*bxa + bya*bya + bza*bza) ;
			if(press==1){
				if(wgam) temp1 = (gam-1.0)*s[2][k][j][i] ;
				else temp1 = cs*cs*s[1][k][j][i];
			}
			else temp1=0;

			return(s[1][k][j][i]*e2z_3(v[1][3],k,j,i)*e2z_3(v[1][3],k,j,i)+(temp1+s[1][k][j][i]*s[3][k][j][i]+0.5*ebi)-e2z_3(v[2][3],k,j,i)*e2z_3(v[2][3],k,j,i)*G3(2,i)*G4(2,j)*dx[1][3][k]*dx[1][1][i]);
    }
    else if(term==28){
      // flux of x1 magnetic flux (x1 component)
			return(0);
    }
    else if(term==29){
      // flux of x1 magnetic flux (x2 component)
			return(0);
    }
    else if(term==30){
      // flux of x1 magnetic flux (x3 component)
			return(0);
    }
    else if(term==31){
      // flux of x2 magnetic flux (x1 component)
			return(0);
    }
    else if(term==32){
      // flux of x2 magnetic flux (x2 component)
			return(0);
    }
    else if(term==33){
      // flux of x2 magnetic flux (x3 component)
			return(0);
    }
    else if(term==34){
      // flux of x3 magnetic flux (x1 component)
			return(0);
    }
    else if(term==35){
      // flux of x3 magnetic flux (x2 component)
			return(0);
    }
    else if(term==36){
      // flux of x3 magnetic flux (x3 component)
			return(0);
    }
  }
  // now begins the volume quantities
  else if(type==3){
    if(term==1){
      // mass
			return(s[1][k][j][i]*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k)) ;
    }
    else if(term==2){
      // kinetic energy (x1 term)
			temp1 = e2z_1(v[1][1],k,j,i);			
			return(0.5*s[1][k][j][i]*temp1*temp1*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k)) ;
    }
    else if(term==3){
      // kinetic energy (x2 term)
			temp1 = e2z_2(v[1][2],k,j,i);			
			return(0.5*s[1][k][j][i]*temp1*temp1*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k)) ;
    }
    else if(term==4){
      // kinetic energy (x3 term)
			temp1 = e2z_3(v[1][3],k,j,i);			
			return(0.5*s[1][k][j][i]*temp1*temp1*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k)) ;
    }
    else if(term==5){
      // internal energy
			/* Equation of state */
			if(press==1){
				if(wgam) ethi = s[2][k][j][i] ;
				else ethi = cs*cs*s[1][k][j][i]*log(s[1][k][j][i]) ;
			}
			else{
				ethi=0;
			}
			return(ethi*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k)) ;
    }
    else if(term==6){
      // potential energy
			//egi = 0.5*s[1][k][j][i]*s[3][k][j][i] ; // for self-gravity
			egi = s[1][k][j][i]*s[3][k][j][i] ;
			return(egi*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k)) ;
    }
    else if(term==7){
      // magnetic energy (x1 term)
			temp1 = e2z_1(v[2][1],k,j,i);    
		  return(0.5*temp1*temp1*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k));
    }
    else if(term==8){
      // magnetic energy (x2 term)
			temp1 = e2z_2(v[2][2],k,j,i);    
		  return(0.5*temp1*temp1*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k));
    }
    else if(term==9){
      // magnetic energy (x3 term)
			temp1 = e2z_3(v[2][3],k,j,i);    
		  return(0.5*temp1*temp1*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k));
    }
    else if(term==10){
      // momentum (x1 term)
			return(s[1][k][j][i]*e2z_1(v[2][1],k,j,i)*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k));
    }
    else if(term==11){
      // momentum (x2 term)
			return(s[1][k][j][i]*e2z_2(v[2][2],k,j,i)*G2(2,i)*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k));
    }
    else if(term==12){
      // momentum (x3 term)
			return(s[1][k][j][i]*e2z_3(v[2][3],k,j,i)*G3(2,i)*G4(2,j)*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k));
    }
    else if(term==13){
      // angular momentum (x1)
			vxa = e2z_1(v[1][1],k,j,i);      
			vya = e2z_2(v[1][2],k,j,i);
			vza = e2z_3(v[1][3],k,j,i) ;
#if(COORD==3)
			return(0); // not yet GODMARK
#elif(COORD==1)
			posx=x[2][1][i];
			posy=x[2][2][j];
			posz=x[2][3][k];
			return(s[1][k][j][i]*sqrt(vza*vza*posy*posy+vya*vya*posz*posz)*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k));
#endif
    }
    else if(term==14){
      // angular momentum (x2)
			vxa = e2z_1(v[1][1],k,j,i);      
			vya = e2z_2(v[1][2],k,j,i);
			vza = e2z_3(v[1][3],k,j,i) ;
#if(COORD==3)
			return(0); // not yet
#elif(COORD==1)
			posx=x[2][1][i];
			posy=x[2][2][j];
			posz=x[2][3][k];
			return(s[1][k][j][i]*sqrt(vxa*vxa*posz*posz+vza*vza*posx*posx)*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k)); 
#endif
    }
    else if(term==15){
      // angular momentum (x3)
			vxa = e2z_1(v[1][1],k,j,i);      
			vya = e2z_2(v[1][2],k,j,i);
			vza = e2z_3(v[1][3],k,j,i) ;
#if(COORD==3)
			return(s[1][k][j][i]*vza*G3(2,i)*G4(2,j)*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k)); // about z
#elif(COORD==1)
    posx=x[2][1][i];
    posy=x[2][2][j];
    posz=x[2][3][k];
    return(s[1][k][j][i]*sqrt(vxa*vxa*posy*posy+vya*vya*posx*posx)*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k)); 
#endif
    }
    else if(term==16){
      // magnetic flux (x1 term)
			return(e2z_1(v[2][1],k,j,i));
    }
    else if(term==17){
      // magnetic flux (x2 term)
			return(e2z_2(v[2][2],k,j,i)/G2(2,i));
    }
    else if(term==18){
      // magnetic flux (x3 term)
			return(e2z_3(v[2][3],k,j,i)/(G3(2,i)*G4(2,j)));
    }
  }
    
}

void divb0check(int which) // should make like rest of diags w.r.t. time and whether it's done
{
  int i,j,k;
  static int switchbc=0;
  int error=0;
  static int firsttime=1;
  int divbmaxpos[3+1];

  static FTYPE DTdivb=0;
  static FTYPE tdivb=0;

  SFTYPE divb1avg,divb2avg;
  SFTYPE divb1max,divb2max;
  SFTYPE divb1avg_full,divb2avg_full;
  SFTYPE divb1max_full,divb2max_full;
  SFTYPE ftemp;



  /////////////////////
  // check on divB=0
  //////
  // now determine divergence of B and output to each logfile
  divb1avg=0;
  divb1max=0;
  divbmaxpos[1]=-3;
  divbmaxpos[2]=-3;
  divbmaxpos[3]=-3;

  LOOPDIVB{
    ftemp=fabs(deldotv(v,2,k,j,i));
    if(ftemp>divb1max){
      divb1max=ftemp;
      divbmaxpos[1]=i;
      divbmaxpos[2]=j;
      divbmaxpos[3]=k;
    }
    divb1avg+=fabs(ftemp);
  }
  /*
    divb2avg=0;
    divb2max=0;
    LOOPDIVB{
    ftemp=fabs(deldotv2(v,2,k,j,i));
    if(ftemp>divb2max) divb2max=ftemp;
    divb2avg+=ftemp;
    }
  */
  if(firsttime){
    //      fprintf(log_file,"IC: t=%21.15g divb1avg: %21.15g divb2avg: %21.15g divb1max: %21.15g divb2max: %21.15g\n",t,divb1avg/((SFTYPE)(realtotalzones)),divb2avg/((SFTYPE)(realtotalzones)),divb1max,divb2max);
    fprintf(log_file,"IC: t=%21.15g divbavg: %21.15g  divbmax: %21.15g k=%d j=%d i=%d\n",t,divb1avg/((SFTYPE)(realtotalzones)),divb1max,divbmaxpos[3],divbmaxpos[2],divbmaxpos[1]);
  }
  else{
    //      fprintf(log_file,"ST: t=%21.15g divb1avg: %21.15g divb2avg: %21.15g divb1max: %21.15g divb2max: %21.15g\n",t,divb1avg/((SFTYPE)(realtotalzones)),divb2avg/((SFTYPE)(realtotalzones)),divb1max,divb2max);
    fprintf(log_file,"ST: t=%21.15g divb1avg: %21.15g divb1max: %21.15g   k=%d j=%d i=%d\n",t,divb1avg/((SFTYPE)(realtotalzones)),divb1max,divbmaxpos[3],divbmaxpos[2],divbmaxpos[1]);
  }
  fflush(log_file);
  // now sum it up and put total in full log file
  if(numprocs>1){
#if(USEMPI)
    MPI_Reduce(&(divb1avg), &(divb1avg_full), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
    //      MPI_Reduce(&(divb2avg), &(divb2avg_full), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
    MPI_Reduce(&(divb1max), &(divb1max_full), 1, MPI_SFTYPE, MPI_MAX, 0,MPI_COMM_WORLD);
    //MPI_Reduce(&(divb2max), &(divb2max_full), 1, MPI_SFTYPE, MPI_MAX, 0,MPI_COMM_WORLD);
    
#endif
  }
  else{
    divb1avg_full=divb1avg;
    //divb2avg_full=divb2avg;
    divb1max_full=divb1max;
    // divb2max_full=divb2max;
  }
  if(myid<=0){
    divb1avg_full/=(SFTYPE)(totalzones);
    //divb2avg_full/=(SFTYPE)(totalzones);
    
    // for diag_ener
    divbavg_full=divb1avg_full;
    divbmax_full=divb1max_full;
    
    if(firsttime){
      //fprintf(logfull_file,"0IC: t=%21.15g divb1avg: %21.15g divb2avg: %21.15g divb1max: %21.15g divb2max: %21.15g\n",t,divb1avg_full,divb2avg_full,divb1max_full,divb2max_full);
      fprintf(logfull_file,"0IC: t=%21.15g divbavg: %21.15g divbmax: %21.15g\n",t,divb1avg_full,divb1max_full);
    }
    else{
      //	fprintf(logfull_file,"0ST: t=%21.15g divb1avg: %21.15g divb2avg: %21.15g divb1max: %21.15g divb2max: %21.15g\n",t,divb1avg_full,divb2avg_full,divb1max_full,divb2max_full);
      fprintf(logfull_file,"0ST: t=%21.15g divbavg: %21.15g divbmax: %21.15g\n",t,divb1avg_full,divb1max_full);
    }
    fflush(logfull_file);
    
  }

  //
  firsttime=0;
}


void sp_compute(void)
{
  int i,j,k;
  FTYPE ftemp,bxa,bya,bza,velfastm2;
  static int firsttime=1;
  int error=0;
  FTYPE spx1,spx2,spx3,sptot;
  static FTYPE Mach2x1l,Mach2x2l,Mach2x3l,Mach2totl;
  static FTYPE Mach2x1l_full,Mach2x2l_full,Mach2x3l_full,Mach2totl_full;
  static FTYPE Mach2x1h,Mach2x2h,Mach2x3h,Mach2toth;
  static FTYPE Mach2x1h_full,Mach2x2h_full,Mach2x3h_full,Mach2toth_full;
  FTYPE Mach2x1,Mach2x2,Mach2x3,Mach2tot;
  char temps[MAXFILENAME];
  FTYPE cs2,vx1,vx2,vx3,vx12,vx22,vx32,vtot2;
  FTYPE signx1,signx2,signx3,signtot;
  static FTYPE tdumpsp;
  int dumpflag;
  static int dumpspc=0;


  // not setup for seemless reentrance
  if(firsttime==1){
    //dumpspc=(int)((t-tstart)/DTsp)+ireenter;
    dumpspc=0;
    //    tdumpsp = tstart+((FTYPE)(dumpspc)*DTsp)-1.E-6;
    //tdumpsp = tstart-1.0E-6;
    tdumpsp = t-1.0E-6;
    //    fprintf(stderr,"%15.10g %d %15.10g\n",t,dumpspc,tdumpsp);

    Mach2x1l=Mach2x2l=Mach2x3l=Mach2totl=100000.0;
    Mach2x1h=Mach2x2h=Mach2x3h=Mach2toth=-100000.0;

    if(myid<=0){
      Mach2x1l_full=Mach2x2l_full=Mach2x3l_full=Mach2totl_full=100000.0;
      Mach2x1h_full=Mach2x2h_full=Mach2x3h_full=Mach2toth_full=-100000.0;

      sprintf(temps,"%s0_logsp%s",DATADIR,DAT2EXT) ;
      if((logsp_file = fopen(temps,WRITETYPE))==NULL) { // just naively append if appendold==1
				fprintf(fail_file,"error opening sp log output file %s\n",temps) ;
				myexit(1) ;
      }
      if((restartsonicpoint==0)||(appendold==0)){
				fprintf(logsp_file,"#%10s\n%10d %10d\n","SPVER",SPVER,SPTYPE);
				fprintf(logsp_file,"#%16s %16s %16s %16s %16s %16s %16s %16s %16s\n","time","Mach2x1l","Mach2x1h","Mach2x2l","Mach2x2h","Mach2x3l","Mach2x3h","Mach2totl","Mach2toth"); 
				fflush(logsp_file);
      }
    }

  }

  if(t>=tdumpsp){
    dumpflag=1;
  }
  else dumpflag=0;


  //////// Check on sonic point on inner edge
  //
  k=0;
  i=0;  // don't choose intix1 since really want what inner grid is doing, which is really what matters
  for(j=0;j<N2;j++){
    if(mag==1){
      /* alfven velocity */
      bxa = e2z_1(v[2][1],k,j,i);
      bya = e2z_2(v[2][2],k,j,i);
      bza = e2z_3(v[2][2],k,j,i);
      ftemp=bxa*bxa + bya*bya + bza*bza; // b^2
      if(wgam) velfastm2 = (gam*(gam-1.)*s[2][k][j][i]+ftemp)/s[1][k][j][i]; // fast magneto-sonic speed
      else velfastm2=cs*cs;
    }
    else{
      if(wgam) velfastm2 = gam*(gam-1.)*s[2][k][j][i]/s[1][k][j][i];
      else velfastm2=cs*cs;
    }
    vx1=v[1][1][k][j][i];
    if(vx1<0) signx1=-1.0; else signx1=1.0;
    vx2=v[1][2][k][j][i];
    if(vx2<0) signx2=-1.0; else signx2=1.0;
    vx3=v[1][3][k][j][i];
    if(vx3<0) signx3=-1.0; else signx3=1.0;

    signtot=1;
    vx12=vx1*vx1;
    vx22=vx2*vx2;
    vx32=vx3*vx3;
    vtot2=vx12+vx22+vx32;
    Mach2x1=signx1*vx12/velfastm2; // sign tells sign of sqrt of this value actually
    Mach2x2=signx2*vx22/velfastm2; // sign tells sign of sqrt of this value actually
    Mach2x3=signx3*vx32/velfastm2; // sign tells sign of sqrt of this value actually
    Mach2tot=signtot*vtot2/velfastm2;


    if(Mach2x1<Mach2x1l){ Mach2x1l=Mach2x1;}
    if(Mach2x2<Mach2x2l){ Mach2x2l=Mach2x2;}
    if(Mach2x3<Mach2x3l){ Mach2x3l=Mach2x3;}
    if(Mach2tot<Mach2totl){ Mach2totl=Mach2tot;}
    
    if(Mach2x1>Mach2x1h){ Mach2x1h=Mach2x1;}
    if(Mach2x2>Mach2x2h){ Mach2x2h=Mach2x2;}
    if(Mach2x3>Mach2x3h){ Mach2x3h=Mach2x3;}
    if(Mach2tot>Mach2toth){ Mach2toth=Mach2tot;}
    
  }

  if(dumpflag==1){
    if(numprocs>1){
#if(USEMPI)
      // send max/min to cpu=0 (since cpu=0 will always be on inner x1-edge, this is ok
      MPI_Reduce(&Mach2x1l,  &Mach2x1l_full , 1, MPI_FTYPE, MPI_MIN, 0,combound[2]);
      MPI_Reduce(&Mach2x2l,  &Mach2x2l_full , 1, MPI_FTYPE, MPI_MIN, 0,combound[2]);
      MPI_Reduce(&Mach2x3l,  &Mach2x3l_full , 1, MPI_FTYPE, MPI_MIN, 0,combound[2]);
      MPI_Reduce(&Mach2totl, &Mach2totl_full, 1, MPI_FTYPE, MPI_MIN, 0,combound[2]);

      MPI_Reduce(&Mach2x1h,  &Mach2x1h_full , 1, MPI_FTYPE, MPI_MAX, 0,combound[2]);
      MPI_Reduce(&Mach2x2h,  &Mach2x2h_full , 1, MPI_FTYPE, MPI_MAX, 0,combound[2]);
      MPI_Reduce(&Mach2x3h,  &Mach2x3h_full , 1, MPI_FTYPE, MPI_MAX, 0,combound[2]);
      MPI_Reduce(&Mach2toth, &Mach2toth_full, 1, MPI_FTYPE, MPI_MAX, 0,combound[2]);
#endif
    }
    else{
      Mach2x1l_full=Mach2x1l;
      Mach2x2l_full=Mach2x2l;
      Mach2x3l_full=Mach2x3l;
      Mach2totl_full=Mach2totl;

      Mach2x1h_full=Mach2x1h;
      Mach2x2h_full=Mach2x2h;
      Mach2x3h_full=Mach2x3h;
      Mach2toth_full=Mach2toth;
    }
    
    if(myid<=0){
      fprintf(logsp_file,"%16.10g %16.10g %16.10g %16.10g %16.10g %16.10g %16.10g %16.10g %16.10g\n",t,Mach2x1l_full,Mach2x1h_full,Mach2x2l_full,Mach2x2h_full,Mach2x3l_full,Mach2x3h_full,Mach2totl_full,Mach2toth_full); 
      fflush(logsp_file);
    }

    
    // clean up for next round of check
    if(myid<=0){
      Mach2x1l_full=Mach2x2l_full=Mach2x3l_full=Mach2totl_full=100000.0;
      Mach2x1h_full=Mach2x2h_full=Mach2x3h_full=Mach2toth_full=-100000.0;
    }
    Mach2x1l=Mach2x2l=Mach2x3l=Mach2totl=100000.0;
    Mach2x1h=Mach2x2h=Mach2x3h=Mach2toth=-100000.0;

    dumpspc++;
    tdumpsp = tstart+(FTYPE)(dumpspc)*DTsp ; 

  }
  //
  ////////

  //
  firsttime=0;
}



// need to make more 1D/2D friendly (N3==1?)

void symmetry_check(int pos)
	// pos:
	// 0: zone center
	// 1: vx pos
	// 2: vy pos
	// 3: vz pos
	// above only setup for scalars
	// -1: analyze but dump out badness into final arrays (used for visual output)
{
  int signchange;
  int countedoffenses;
  int globalsym; // -1: none, else that #
  int symmetrytypes[NUMSCA+1];
  int symmetrytypev[NUMVEC+1][3+1];
  SFTYPE symbreaks[1+NUMSCA];
  SFTYPE symbreakv[NUMVEC+1][1+3];
  int i,j,k,ll,m;
  int counter;
  SFTYPE offense,offenseold;
  int whooffends[3+1];
  FTYPE (*stouse)[N3M][N2M][N1M] ;
  FTYPE (*vtouse)[3][N3M][N2M][N1M] ;
  int badnessoutput;


  if(COORD!=1) return; // only needed for coord==1



  if(pos>=0){
    badnessoutput=0;
    stouse=s;
    vtouse=v;
  }
  else{
    pos=0; // assume normal variables
    badnessoutput=1;
    stouse=sanal;
    vtouse=vanal;
    for(ll=1;ll<=NUMSCA;ll++){
      LOOPF{
				sanal[ll][k][j][i]=s[ll][k][j][i];
      }
    }
    for(ll=1;ll<=NUMVEC;ll++){
      for(m=1;m<=3;m++){
				LOOPF{
					vanal[ll][m][k][j][i]=v[ll][m][k][j][i];
				}
      }
    }
    // now 0 out symmetry badness holders(use these for easy outputting)
    for(ll=1;ll<=NUMSCA;ll++){
      LOOPF{
				s[ll][k][j][i]=0;
      }
    }
    for(ll=1;ll<=NUMVEC;ll++){
      for(m=1;m<=3;m++){
				LOOPF{
					v[ll][m][k][j][i]=0;
				}
      }
    }

  }
  // user defined things dependent on IC/BC/GRID
  globalsym=0;

  for(ll=1;ll<=NUMSCA;ll++){
    symmetrytypes[ll]=globalsym;
  }
  for(ll=1;ll<=NUMVEC;ll++){
    for(m=1;m<=3;m++){
      symmetrytypev[ll][m]=globalsym;
    }
  }  

  // 0: full symmetry
  // 1: only 180rot w/ parity, not 90deg

  fprintf(log_file,"t=%10.5g ",t); fflush(log_file);


  if(COORD==1){

    /////////////
    // check for z + 90degree rotation symmetry assuing even grids
    /////////////
    for(ll=1;ll<=NUMSCA;ll++){
      symbreaks[ll]=0.0;
    }
    for(ll=1;ll<=NUMVEC;ll++){
      // v_1
      for(m=1;m<=3;m++){
				symbreakv[ll][m]=0.0;
      }
    }


    // scalars
    for(ll=1;ll<=NUMSCA;ll++){
      if(pos==0){
				counter=0;
				countedoffenses=0;
				offenseold=0; whooffends[1]=-3;whooffends[2]=-3;whooffends[3]=-3;
				// LOOPF{ // for full grid test ( assumes nothing bad happens in nowhere zones!
				LOOPSUPERGEN(5){ // good for compgrid test
					if((symmetrytypes[ll]==0)||(symmetrytypes[ll]==1)){
						offense=(
							// for real scalars
							+fabs(stouse[ll][k][j][i]-stouse[ll][k][N2-1-j][N1-1-i]) // 180rot
							+fabs(stouse[ll][k][j][i]-stouse[ll][N3-1-k][j][i]) // z-parity
							+fabs(stouse[ll][k][j][i]-stouse[ll][N3-1-k][N2-1-j][N1-1-i]) // 180rot + z-parity
							);
						counter+=3;
					}
					if(symmetrytypes[ll]==0){
						if(N1==N2){
							offense+=(
								// for real scalars
								+fabs(stouse[ll][k][j][i]-stouse[ll][k][N2-1-i][j]) // 270rot
								+fabs(stouse[ll][k][j][i]-stouse[ll][k][i][N1-1-j]) // 90rot
								+fabs(stouse[ll][k][j][i]-stouse[ll][N3-1-k][N2-1-i][j]) // 270rot + z-parity
								+fabs(stouse[ll][k][j][i]-stouse[ll][N3-1-k][i][N1-1-j]) // 90rot + z-parity
								);
							counter+=4;
						}
					}
					symbreaks[ll]+=offense;
					if(offense>offenseold){
						offenseold=offense;
						whooffends[1]=i;
						whooffends[2]=j;
						whooffends[3]=k;
					}
					if(offense>100.0*NUMEPSILON){
						countedoffenses++;
					}
					if(badnessoutput) s[ll][k][j][i]=offense;
				}
				symbreaks[ll]/=((SFTYPE)(counter));
      }
      else if(pos==1){ // scalars at v_x position
				counter=0;
				countedoffenses=0;
				offenseold=0; whooffends[1]=-3;whooffends[2]=-3;whooffends[3]=-3;
				//	LOOPF3 LOOPF2 LOOPHMFP{  // assumes nothing bad happens in nowhere zones!
				LOOPSUPERGEN(5){ // doesn't directly loop over outer edge, but gets it anyways indirectly
					offense=0;
					if((symmetrytypes[ll]==0)||(symmetrytypes[ll]==1)){
						offense=(
							+fabs(stouse[ll][k][j][i]-stouse[ll][k][N2-1-j][N1-i]) // 180rot
							+fabs(stouse[ll][k][j][i]-stouse[ll][N3-1-k][j][i]) // z-parity
							+fabs(stouse[ll][k][j][i]-stouse[ll][N3-1-k][N2-1-j][N1-i]) // 180rot + z-parity
							);
					}
					if(symmetrytypes[ll]==0){
						if(N1==N2){
							offense+=(
								+fabs(stouse[ll][k][j][i]-stouse[ll][k][N2-1-i][j]) // 270rot
								+fabs(stouse[ll][k][j][i]-stouse[ll][k][i][N1-j]) // 90rot
								+fabs(stouse[ll][k][j][i]-stouse[ll][N3-1-k][N2-1-i][j]) // 270rot + z-parity
								+fabs(stouse[ll][k][j][i]-stouse[ll][N3-1-k][i][N1-j]) // 90rot + z-parity
								);
						}
					}
					symbreaks[ll]+=offense;
					if(offense>offenseold){
						offenseold=offense;
						whooffends[1]=i;
						whooffends[2]=j;
						whooffends[3]=k;
					}
					if(offense>100.0*NUMEPSILON){
						countedoffenses++;
					}
					if(badnessoutput) s[ll][k][j][i]=offense;
				}
				symbreaks[ll]/=((SFTYPE)(counter));
      }
      else if(pos==2){ // scalars at v_y position
				counter=0;
				countedoffenses=0;
				offenseold=0; whooffends[1]=-3;whooffends[2]=-3;whooffends[3]=-3;
				//	LOOPF3 LOOPHMFP2 LOOPF1{ // as above
				LOOPSUPERGEN(5){
					offense=0;
					if((symmetrytypes[ll]==0)||(symmetrytypes[ll]==1)){
						offense=(
							// for real scalars
							+fabs(stouse[ll][k][j][i]-stouse[ll][k][N2-j][N1-1-i]) // 180rot
							+fabs(stouse[ll][k][j][i]-stouse[ll][N3-1-k][j][i]) // z-parity
							+fabs(stouse[ll][k][j][i]-stouse[ll][N3-1-k][N2-j][N1-1-i]) // 180rot + z-parity
							);
						counter+=3;
					}
					if(symmetrytypes[ll]==0){
						if(N1==N2){
							offense+=(
								// for real scalars
								+fabs(stouse[ll][k][j][i]-stouse[ll][k][N2-i][j]) // 270rot
								+fabs(stouse[ll][k][j][i]-stouse[ll][k][i][N1-1-j]) // 90rot
								+fabs(stouse[ll][k][j][i]-stouse[ll][N3-1-k][N2-i][j]) // 270rot + z-parity
								+fabs(stouse[ll][k][j][i]-stouse[ll][N3-1-k][i][N1-1-j]) // 90rot + z-parity
								);
							counter+=4;
						}
					}
					symbreaks[ll]+=offense;
					if(offense>offenseold){
						offenseold=offense;
						whooffends[1]=i;
						whooffends[2]=j;
						whooffends[3]=k;
					}
					if(offense>100.0*NUMEPSILON){
						countedoffenses++;
					}
					if(badnessoutput) s[ll][k][j][i]=offense;
				}
				symbreaks[ll]/=((SFTYPE)(counter));
      }
      else if(pos==3){ // scalars at v_z position
				counter=0;
				countedoffenses=0;
				offenseold=0; whooffends[1]=-3;whooffends[2]=-3;whooffends[3]=-3;
				//	LOOPHMFP3 LOOPF2 LOOPF1{ // as above
				LOOPSUPERGEN(5){
					offense=0;
					if((symmetrytypes[ll]==0)||(symmetrytypes[ll]==1)){
						offense=(
							// for real scalars
							+fabs(stouse[ll][k][j][i]-stouse[ll][k][N2-1-j][N1-1-i]) // 180rot
							+fabs(stouse[ll][k][j][i]-stouse[ll][N3-k][j][i]) // z-parity
							+fabs(stouse[ll][k][j][i]-stouse[ll][N3-k][N2-1-j][N1-1-i]) // 180rot + z-parity
							);
						counter+=3;
					}
					if(symmetrytypes[ll]==0){
						if(N1==N2){
							offense+=(
								// for real scalars
								+fabs(stouse[ll][k][j][i]-stouse[ll][k][N2-1-i][j]) // 270rot
								+fabs(stouse[ll][k][j][i]-stouse[ll][k][i][N1-1-j]) // 90rot
								+fabs(stouse[ll][k][j][i]-stouse[ll][N3-k][N2-1-i][j]) // 270rot + z-parity
								+fabs(stouse[ll][k][j][i]-stouse[ll][N3-k][i][N1-1-j]) // 90rot + z-parity
								);
							counter+=4;
						}
					}
					symbreaks[ll]+=offense;
					if(offense>offenseold){
						offenseold=offense;
						whooffends[1]=i;
						whooffends[2]=j;
						whooffends[3]=k;
					}
					if(offense>100.0*NUMEPSILON){
						countedoffenses++;
					}
					if(badnessoutput) s[ll][k][j][i]=offense;
				}
				symbreaks[ll]/=((SFTYPE)(counter));
      }
      fprintf(log_file,"ll=%1d sb: %15.10g  (%7d) (%3d %3d %3d) ",ll,symbreaks[ll],countedoffenses,whooffends[3],whooffends[2],whooffends[1]);
    }

    // vectors
    for(ll=1;ll<=NUMVEC;ll++){
      // v_1
      if(ll==1) signchange=1;
      if(ll==2) signchange=-1; // for parity of field vx vy
      counter=0;
      countedoffenses=0;
      offenseold=0; whooffends[1]=-3;whooffends[2]=-3;whooffends[3]=-3;
      //      LOOPF3 LOOPF2 LOOPHMFP{ // as above
      LOOPSUPERGEN(20){ // good for b or v really
				offense=0;
				if((symmetrytypev[ll][1]==0)||(symmetrytypev[ll][1]==1)){
					offense=(
						+fabs(vtouse[ll][1][k][j][i]+vtouse[ll][1][k][N2-1-j][N1-i]) // 180deg rotation
						+fabs(vtouse[ll][1][k][j][i]+(signchange)*vtouse[ll][1][N3-1-k][N2-1-j][N1-i]) // 180deg rotation + parity in z
						+fabs(vtouse[ll][1][k][j][i]-(signchange)*vtouse[ll][1][N3-1-k][j][i]) // parity in z only (no sign change)
						);
					counter+=3;
				}
				if(symmetrytypev[ll][1]==0){
					// these terms only true if N1==N2
					if(N1==N2){
						offense+=(
							+fabs(vtouse[ll][1][k][j][i]-vtouse[ll][2][k][i][N1-1-j]) // 90deg rotation (no sign change)
							+fabs(vtouse[ll][1][k][j][i]+vtouse[ll][2][k][N2-i][j]) // 270deg rotation (sign change)
							// + parity in z (no sign change)
							+fabs(vtouse[ll][1][k][j][i]-(signchange)*vtouse[ll][2][N3-1-k][i][N1-1-j]) // 90deg rotation
							+fabs(vtouse[ll][1][k][j][i]+(signchange)*vtouse[ll][2][N3-1-k][N2-i][j]) // 270deg rotation
							);
						counter+=4;
					}
				}
				symbreakv[ll][1]+=offense;
				if(offense>offenseold){
					offenseold=offense;
					whooffends[1]=i;
					whooffends[2]=j;
					whooffends[3]=k;
				}
				if(offense>100.0*NUMEPSILON){
					countedoffenses++;
				}
				if(badnessoutput) v[ll][1][k][j][i]=offense;
      }
      symbreakv[ll][1]/=((SFTYPE)(counter));
      fprintf(log_file,"ll=%1d c=%d sb: %15.10g (%7d) (%3d %3d %3d) ",ll,1,symbreakv[ll][1],countedoffenses,whooffends[3],whooffends[2],whooffends[1]);


      // v_2
      if(ll==1) signchange=1;
      if(ll==2) signchange=-1; // for parity of field vx vy
      countedoffenses=0;
      counter=0;
      offenseold=0; whooffends[1]=-3;whooffends[2]=-3;whooffends[3]=-3;
      //LOOPF3 LOOPHMFP2 LOOPF1{ // as above
      LOOPSUPERGEN(21){ // good for b or v really
				offense=0;
				if((symmetrytypev[ll][2]==0)||(symmetrytypev[ll][2]==1)){
					offense=( // the 180's
						+fabs(vtouse[ll][2][k][j][i]+vtouse[ll][2][k][N2-j][N1-1-i]) // 180deg rotation(sign change)
						+fabs(vtouse[ll][2][k][j][i]+(signchange)*vtouse[ll][2][N3-1-k][N2-j][N1-1-i]) // 180deg rotation + parity in z
						+fabs(vtouse[ll][2][k][j][i]-(signchange)*vtouse[ll][2][N3-1-k][j][i]) // parity in z ony (no sign change)
						);
					counter+=3;
				}
				if(symmetrytypev[ll][2]==0){
					// true 90deg rot
					// these terms only true if N1==N2
					if(N1==N2){ // then 90/270's
						offense+=(
							+fabs(vtouse[ll][2][k][j][i]+vtouse[ll][1][k][i][N1-j]) // 90deg rotation (sign change)
							+fabs(vtouse[ll][2][k][j][i]-vtouse[ll][1][k][N2-1-i][j]) // 270deg rotation (no sign change)
							// + parity in z
							+fabs(vtouse[ll][2][k][j][i]+(signchange)*vtouse[ll][1][N3-1-k][i][N1-j]) // 90deg rotation
							+fabs(vtouse[ll][2][k][j][i]-(signchange)*vtouse[ll][1][N3-1-k][N2-1-i][j]) // 270deg rotation
							);
						counter+=4;
					}
				}
				symbreakv[ll][2]+=offense;
				if(offense>offenseold){
					offenseold=offense;
					whooffends[1]=i;
					whooffends[2]=j;
					whooffends[3]=k;
				}
				if(offense>100.0*NUMEPSILON){
					countedoffenses++;
				}
				if(badnessoutput) v[ll][2][k][j][i]=offense;
      }
      symbreakv[ll][2]/=((SFTYPE)(counter));
      fprintf(log_file,"ll=%1d c=%d sb: %15.10g (%7d) (%3d %3d %3d) ",ll,2,symbreakv[ll][2],countedoffenses,whooffends[3],whooffends[2],whooffends[1]);

      // v_3
      countedoffenses=0;
      counter=0;
      if(ll==1) signchange=-1; // sign change for parity
      if(ll==2) signchange=1; // no sign change for parity
      offenseold=0; whooffends[1]=-3;whooffends[2]=-3;whooffends[3]=-3;
      //LOOPHMFP3 LOOPF2 LOOPF1{ // as above
      LOOPSUPERGEN(22){ // good for b or v really
				offense=0;
				if((symmetrytypev[ll][3]==0)||(symmetrytypev[ll][3]==1)){
					offense=(
						+fabs(vtouse[ll][3][k][j][i]-vtouse[ll][3][k][N2-1-j][N1-1-i]) // rot180 (no change in sign)
						+fabs(vtouse[ll][3][k][j][i]-(signchange)*vtouse[ll][3][N3-k][N2-1-j][N1-1-i]) // parity + rot180
						+fabs(vtouse[ll][3][k][j][i]-(signchange)*vtouse[ll][3][N3-k][j][i]) // parity only(change in sign)
						);
					counter+=3;
				}
				if(symmetrytypev[ll][3]==0){
					// true 90deg rot
					if(N1==N2){
						offense+=(
							+fabs(vtouse[ll][3][k][j][i]-vtouse[ll][3][k][i][N1-1-j]) // rot90 (no change in sign)
							+fabs(vtouse[ll][3][k][j][i]-vtouse[ll][3][k][N2-1-i][j]) // rot270 (no change in sign)
							// + parity in z ( change in sign)
							+fabs(vtouse[ll][3][k][j][i]-(signchange)*vtouse[ll][3][N3-k][i][N1-1-j]) // rot90
							+fabs(vtouse[ll][3][k][j][i]-(signchange)*vtouse[ll][3][N3-k][N2-1-i][j]) // rot270
							);
						counter+=4;
					}
				}
				symbreakv[ll][3]+=offense;
				if(offense>offenseold){
					offenseold=offense;
					whooffends[1]=i;
					whooffends[2]=j;
					whooffends[3]=k;
				}
				if(offense>100.0*NUMEPSILON){
					countedoffenses++;
				}
				if(badnessoutput) v[ll][3][k][j][i]=offense;
      }
      symbreakv[ll][3]/=((SFTYPE)(counter));
      fprintf(log_file,"ll=%1d c=%d sb: %15.10g  (%7d) (%3d %3d %3d) ",ll,3,symbreakv[ll][3],countedoffenses,whooffends[3],whooffends[2],whooffends[1]);
    }
    fprintf(log_file,"\n");
  }// end if coord==1

  // nothing if coord==3


}
  

void crazy_check(void)
{
  int i,j,k,l,m;
  int crazy[NUMSCA+NUMVEC*3];

  // a crazy location
  k=31;
  j=32;
  i=26;


  if(s[1][k][j][i]>1.0){
    crazy[0]=1;
  }
  else crazy[0]=0;

  if(s[2][k][j][i]>1.0E-2){
    crazy[1]=1;
  }
  else crazy[1]=0;

  if(v[1][1][k][j][i]>1.0){
    crazy[2]=1;
  }
  else crazy[2]=0;

  if(v[1][2][k][j][i]>1.0){
    crazy[3]=1;
  }
  else crazy[3]=0;

  if(v[1][3][k][j][i]>1.0){
    crazy[4]=1;
  }
  else crazy[4]=0;

  // b-field
  if(v[2][1][k][j][i]>1.0E-2){
    crazy[5]=1;
  }
  else crazy[5]=0;

  if(v[2][2][k][j][i]>1.0E-2){
    crazy[6]=1;
  }
  else crazy[6]=0;

  if(v[2][3][k][j][i]>1.0E-1){
    crazy[7]=1;
  }
  else crazy[7]=0;

  crazy[8]=0;


  // output who went crazy
  for(l=0;l<NUMSCA+NUMVEC*3-1;l++){
    fprintf(stdout,"%d ",crazy[l]);
  }
  fprintf(stdout,"\n");
  fflush(stdout);

}


#if(USEMPI)

void jonio_init_fp(FILE **fp,int which,char *filename)
{
	if(myid==0){ // total on CPU=0	
		if(which==1) *fp = fopen(filename,"w") ;
		else if(which==2) *fp = fopen(filename,"rb") ;
		if(*fp==NULL) {
	    fprintf(fail_file,"error opening file: %s\n",filename) ;
	    myexit(2) ;
		}
	}
}
void jonio_init_mem(int numcolumns,int datatype,unsigned char **jonio1,float **jonio4, double **jonio8,unsigned char **writebuf1,float **writebuf4,double **writebuf8,unsigned char **tempbuf1, float **tempbuf4,double **tempbuf8)
{
	// based on sizeof()
	// 1: unsigned char
	// 4: float
	// 8: double
	void * jonio,*writebuf,*tempbuf;
 
	int sizeofmemory;


//    fprintf(MYOUT,"proc: %d : got here-1\n",myid); fflush(MYOUT);

//    fprintf(MYOUT,"proc: %d : datatype: %d\n",myid,datatype); fflush(MYOUT);

	if(myid==0){ // total on CPU=0
//	fprintf(MYOUT,"proc: %d : got after fp\n"); fflush(MYOUT);
		// initialize memory for jonio
		sizeofmemory=datatype*totalsize[1]*totalsize[2]*totalsize[3]*numcolumns;
		if(datatype==1) *jonio1=jonio=(unsigned char *)malloc(sizeofmemory);
		else if(datatype==4) *jonio4=jonio=(float *)malloc(sizeofmemory);
		else if(datatype==8) *jonio8=jonio=(double *)malloc(sizeofmemory);
		if(jonio==NULL){
	    fprintf(fail_file,"Can't initialize jonio memory\n");
	    myexit(1);
		}
	}
//    if(myid==0){ fprintf(MYOUT,"proc: %d : got here-2: %d\n",myid,jonio); fflush(MYOUT);}
//    else{  fprintf(MYOUT,"proc: %d : got here-2\n",myid); fflush(MYOUT);}

	sizeofmemory=datatype*N1*N2*N3*numcolumns;
	if(datatype==1) *writebuf1=writebuf=(unsigned char *)malloc(sizeofmemory);
	else if(datatype==4) *writebuf4=writebuf=(float *)malloc(sizeofmemory);
	else if(datatype==8) *writebuf8=writebuf=(double *)malloc(sizeofmemory);
	if(writebuf==NULL){
		fprintf(fail_file,"Can't initialize writebuf memory\n");
		myexit(1);
	}
//    fprintf(MYOUT,"proc: %d : got here-3 : %d\n",myid,writebuf); fflush(MYOUT);
	if(myid==0){ // tempbuf only used by CPU=0
		sizeofmemory=datatype*N1*N2*N3*numcolumns;
		if(datatype==1) *tempbuf1=tempbuf=(unsigned char *)malloc(sizeofmemory);
		else if(datatype==4) *tempbuf4=tempbuf=(float *)malloc(sizeofmemory);
		else if(datatype==8) *tempbuf8=tempbuf=(double *)malloc(sizeofmemory);
		if(tempbuf==NULL){
	    fprintf(fail_file,"Can't initialize tempbuf memory\n");
	    myexit(1);
		}
//	fprintf(MYOUT,"proc: %d : got here-4 : %d\n",myid,tempbuf); fflush(MYOUT);
	}
}

void jonio_combine(int stage,MPI_Datatype mpidt,int numcolumns,int datatype,FILE* fp, void * jonio, void * writebuf, void * tempbuf)
{
	// based on sizeof()
	// 1: unsigned char
	// 4: float
	// 8: double
	int i,j,k,l,col,mapvaluejonio,mapvaluetempbuf;
#if(USEMPI)
	MPI_Request rrequest[MAXCPUS];
	MPI_Request srequest;
#endif
	int othercpupos[3+1];
	unsigned char *jonio1;
	float *jonio4;
	double *jonio8;
	unsigned char *writebuf1;
	float *writebuf4;
	double *writebuf8;
	unsigned char *tempbuf1;
	float *tempbuf4;
	double *tempbuf8;
    

	if(datatype==1) jonio1=(unsigned char*)jonio;
	else if(datatype==4) jonio4=(float*)jonio;
	else if(datatype==8) jonio8=(double*)jonio;
	if(datatype==1) writebuf1=(unsigned char*)writebuf;
	else if(datatype==4) writebuf4=(float*)writebuf;
	else if(datatype==8) writebuf8=(double*)writebuf;
	if(datatype==1) tempbuf1=(unsigned char*)tempbuf;
	else if(datatype==4) tempbuf4=(float*)tempbuf;
	else if(datatype==8) tempbuf8=(double*)tempbuf;

	if(stage==1){

//    fprintf(MYOUT,"0proc: %d : %d %d %d : %d %d %d\n",myid,jonio,writebuf,tempbuf,datatype,mpidt==MPI_FTYPE,mpidt==MPI_BYTE); fflush(MYOUT);
#if(USEMPI)
    MPI_Isend(writebuf,N1*N2*N3*numcolumns,mpidt,0,myid,MPI_COMM_WORLD,&srequest);
    if(myid==0){
			for(l=0;l<numprocs;l++){
//	    fprintf(MYOUT,"proc: %d : begin: l=%d\n",myid,l); fflush(MYOUT);
				MPI_Irecv(tempbuf,N1*N2*N3*numcolumns,mpidt,l,l,MPI_COMM_WORLD,&rrequest[l]);
				MPI_Wait(&rrequest[l],&mpichstatus);
	    
				othercpupos[1]=l%ncpux1;
				othercpupos[2]=(int)((l%(ncpux1*ncpux2))/ncpux1);
				othercpupos[3]=(int)(l/(ncpux1*ncpux2));
				// now fill jonio with proper sequence (i.e. tiled mapping)
				for(k=0;k<N3;k++) for(j=0;j<N2;j++) for(i=0;i<N1;i++) for(col=0;col<numcolumns;col++){
					mapvaluejonio=
						+ncpux1*ncpux2*N1*N2*numcolumns*(k+othercpupos[3]*N3)
						+ncpux1       *N1   *numcolumns*(j+othercpupos[2]*N2)
						+                    numcolumns*(i+othercpupos[1]*N1)
						+col;
					mapvaluetempbuf=k*N1*N2*numcolumns+j*N1*numcolumns+i*numcolumns+col;
		
//		fprintf(MYOUT,"proc: %d : inside: %d %d : %d %d %d\n",myid,mapvaluejonio,mapvaluetempbuf,k,j,i); fflush(MYOUT);
 
					if(datatype==1) jonio1[mapvaluejonio]=tempbuf1[mapvaluetempbuf];
					if(datatype==4) jonio4[mapvaluejonio]=tempbuf4[mapvaluetempbuf];		
					if(datatype==8) jonio8[mapvaluejonio]=tempbuf8[mapvaluetempbuf];
				}
//	    fprintf(MYOUT,"proc: %d : end: l=%d\n",myid,l); fflush(MYOUT);
			}
    }
    MPI_Wait(&srequest,&mpichstatus);

//    fprintf(MYOUT,"1proc: %d : %d %d %d\n",myid,jonio,writebuf,tempbuf); fflush(MYOUT);
    if(myid==0){
			// now write out collected data using CPU=0
			fwrite(jonio,datatype,totalsize[1]*totalsize[2]*totalsize[3]*numcolumns,fp);
    }
	}
	else if(stage==2){
		if(myid==0){
	    fclose(fp);
	    fp=NULL;
		}
	}
	else if(stage==3){
		free(writebuf);
		if(myid==0){
	    free(tempbuf);
	    free(jonio);
		}
	}
    
//    fprintf(MYOUT,"2proc: %d : %d %d %d\n",myid,jonio,writebuf,tempbuf); fflush(MYOUT);
#endif

}

void jonio_seperate(int stage,MPI_Datatype mpidt,int numcolumns,int datatype,FILE* fp, void * jonio, void * writebuf, void * tempbuf)
{
	// baesd on sizeof()
	// 1: unsigned char
	// 4: float
	// 8: double
	int i,j,k,l,col,mapvaluejonio,mapvaluetempbuf;
#if(USEMPI)
	MPI_Request rrequest;
	MPI_Request srequest;
#endif
	int othercpupos[3+1];
	unsigned char *jonio1;
	float *jonio4;
	double *jonio8;
	unsigned char *writebuf1;
	float *writebuf4;
	double *writebuf8;
	unsigned char *tempbuf1;
	float *tempbuf4;
	double *tempbuf8;
    

	if(datatype==1) jonio1=(unsigned char*)jonio;
	else if(datatype==4) jonio4=(float*)jonio;
	else if(datatype==8) jonio8=(double*)jonio;
	if(datatype==1) writebuf1=(unsigned char*)writebuf;
	else if(datatype==4) writebuf4=(float*)writebuf;
	else if(datatype==8) writebuf8=(double*)writebuf;
	if(datatype==1) tempbuf1=(unsigned char*)tempbuf;
	else if(datatype==4) tempbuf4=(float*)tempbuf;
	else if(datatype==8) tempbuf8=(double*)tempbuf;

#if(USEMPI)

	if(stage==1){
		if(myid==0){
	    // first let cpu=0 read data
//	    fread(jonio,datatype,totalsize[1]*totalsize[2]*totalsize[3]*numcolumns,fp);
	    for(k=0;k<totalsize[3];k++) for(j=0;j<totalsize[2];j++) for(i=0;i<totalsize[1];i++) for(col=0;col<numcolumns;col++){
				fread(&jonio4[k*numcolumns*totalsize[2]*totalsize[1]+j*numcolumns*totalsize[1]+i*numcolumns+col],sizeof(FTYPE),1,fp);
//	    for(i=0;i<numcolumns*totalsize[3]*totalsize[2]*totalsize[1];i++){
//		fread(&jonio4[i],sizeof(FTYPE),1,fp);
//k*numcolumns*totalsize[2]*totalsize[1]+j*numcolumns*totalsize[1]+i*numcolumns+1E-6;
	    }
		}
		for(l=0;l<numprocs;l++){
	    if(myid==0){
				othercpupos[1]=l%ncpux1;
				othercpupos[2]=(int)((l%(ncpux1*ncpux2))/ncpux1);
				othercpupos[3]=(int)(l/(ncpux1*ncpux2));
		
				// now unfill jonio with proper sequence (i.e. tiled mapping)
				for(k=0;k<N3;k++) for(j=0;j<N2;j++) for(i=0;i<N1;i++) for(col=0;col<numcolumns;col++){
					mapvaluejonio=
						ncpux1*ncpux2*N1*N2*numcolumns*(k+othercpupos[3]*N3)
						+ncpux1       *N1   *numcolumns*(j+othercpupos[2]*N2)
						+                    numcolumns*(i+othercpupos[1]*N1)
						+col;
					//mapvaluejonio=N1*numcolumns*j+numcolumns*i+col;
					mapvaluetempbuf=k*N1*N2*numcolumns+j*N1*numcolumns+i*numcolumns+col;
					if(datatype==1) tempbuf1[mapvaluetempbuf]=jonio1[mapvaluejonio];
					if(datatype==4) tempbuf4[mapvaluetempbuf]=jonio4[mapvaluejonio];		
					if(datatype==8) tempbuf8[mapvaluetempbuf]=jonio8[mapvaluejonio];
		    
					//tempbuf4[mapvaluetempbuf]=mapvaluetempbuf/numcolumns/((SFTYPE)(N1*N2*numcolumns))*256.0/4.0+l*256.0/4.0+1E-6;
//l*256/4+255/4;
//		fprintf(log_file,"%d %d\n",mapvaluejonio,mapvaluetempbuf); fflush(log_file);
				}
		
				MPI_Isend(tempbuf,N1*N2*N3*numcolumns,mpidt,l,l,MPI_COMM_WORLD,&srequest);
	    }
	    if(myid==l) MPI_Irecv(writebuf,N1*N2*N3*numcolumns,mpidt,0,myid,MPI_COMM_WORLD,&rrequest);
	    if(myid==0) MPI_Wait(&srequest,&mpichstatus);
	    if(myid==l) MPI_Wait(&rrequest,&mpichstatus); // writebuf used until stage 2
//	    for(k=0;k<N3;k++) for(j=0;j<N2;j++) for(i=0;i<N1;i++) for(col=0;col<numcolumns;col++){
//		mapvaluetempbuf=k*N1*N2*numcolumns+j*N1*numcolumns+i*numcolumns+col;
			//writebuf4[mapvaluetempbuf]=myid*256/4+255/4;
			//writebuf4[mapvaluetempbuf]=mapvaluetempbuf/numcolumns/((SFTYPE)(N1*N2*numcolumns))*256.0/4.0+myid*256.0/4.0+1E-6;
//	    }
	    if(myid==0) fprintf(log_file,"jonio_sep: from 0 to %d: %d %d %d\n",l,othercpupos[1],othercpupos[2],othercpupos[3]); fflush(log_file);
		}
	}
	else if(stage==2){
		if(myid==0){
	    fclose(fp);
	    fp=NULL;
		}
	}
	else if(stage==3){
		if(myid==0) free(jonio); // put here since above used multiple times for different files but same memory
		free(writebuf);
	}
#endif



}

#endif
