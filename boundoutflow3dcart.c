#include "bound.h"
#include "boundrect.h"


// simplified 3D outflow BC code

void bound_rect1(FTYPE (*vars)[N2M][N1M],
	   FTYPE (*varv)[N3M][N2M][N1M],
	   int wsca,
	   int wvec,
	   int wcom) // 1,2,3, 0: none, 123=all 12=12 13=13, 23=23 (not currently setup to do 13 or 23 since never needed)
{
  static int firsttime=1;
  int i,j,k,l,m,p,q ;
  int ii,jj,kk,ll;
  FTYPE (*works)[N3M][N2M][N1M];
  FTYPE (*workv)[3][N3M][N2M][N1M];
  FTYPE ftemp;
  static int numhit=0;



  numhit++;
  if(wsca<=-2){
    works=(FTYPE (*) [N3M][N2M][N1M])(&vars[0][0][0]);
  }
  else works=s;
  if(wvec<=-2){
    workv=(FTYPE (*) [3][N3M][N2M][N1M])(&varv[0][0][0][0]);
  }
  else workv=v;



  /////////////////////////////////////
  /////////////////////////////////////
  ////////////// SCALARS
  /////////////////////////////////////
  /////////////////////////////////////

  if(wsca!=0){

    for(l=1;l<=NUMSCA;l++){ // bounds potential
    
      if(wsca!=-1){ if(wsca<=-2) ll=0; else ll=wsca; }
      else ll=l;
      
      
      
      // x1 down
      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<0;i++){
	kk=k; jj=j; ii=0;
	works[ll][k][j][i]=works[ll][kk][jj][ii];
      }
      // x1 up
      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=N1;i<N1+N1BND;i++){
	kk=k; jj=j; ii=N1-1;
	works[ll][k][j][i]=works[ll][kk][jj][ii];
      }

      // x2 down
      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<0;j++) for(i=-N1BND;i<N1+N1BND;i++){
	kk=k; jj=0; ii=i;
	works[ll][k][j][i]=works[ll][kk][jj][ii];
      }
      // x2 up
      for(k=-N3BND;k<N3+N3BND;k++) for(j=N2;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	kk=k; jj=N2-1; ii=i;
	works[ll][k][j][i]=works[ll][kk][jj][ii];
      }

      // x3 down
      for(k=-N3BND;k<0;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	kk=0; jj=j; ii=i;
	works[ll][k][j][i]=works[ll][kk][jj][ii];
      }
      // x3 up
      for(k=N3;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	kk=N3-1; jj=j; ii=i;
	works[ll][k][j][i]=works[ll][kk][jj][ii];
      }

    
      // cut short loop if only 1 scalar
      if(wsca!=-1) l=NUMSCA;
    } // end over scalars
  } // end if any scalars
  

  










  ////////////////////////////////////////////////
  /////////////////////////////////////////////////
  ////////////////////////////////////////////////
  //////////////////////////// VECTORS
  ///////////////////////////////////////////////
  /////////////////////////////////////////////////
  ////////////////////////////////////////////////
  // Now do vectors if any
  if(wvec!=0){
    

    for(l=1;l<=REALNUMVEC;l++){
      
      /* if not to do all, pick */
      if(wvec!=-1){ if(wvec<=-2) ll=0; else ll=wvec; }
      else ll=l;
      
	
      // x1 down
      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<=0;i++){
	kk=k; jj=j; ii=1;
	workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	if(ll==1 && workv[ll][1][k][j][i]>0.0) workv[ll][1][k][j][i]=0.0;

	kk=k; jj=j; ii=0;
	workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];

      }
      // x1 up
      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=N1;i<N1+N1BND;i++){
	kk=k; jj=j; ii=N1-1;
	workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];

	if(ll==1 && workv[ll][1][k][j][i]<0.0) workv[ll][1][k][j][i]=0.0;

	kk=k; jj=j; ii=N1-1;
	workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];

      }

      // x2 down
      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<=0;j++) for(i=-N1BND;i<N1+N1BND;i++){
	kk=k; jj=1; ii=i;
	workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];

	if(ll==1 && workv[ll][2][k][j][i]>0.0) workv[ll][2][k][j][i]=0.0;

	kk=k; jj=0; ii=i;
	workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];

      }
      // x2 up
      for(k=-N3BND;k<N3+N3BND;k++) for(j=N2;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	kk=k; jj=N2-1; ii=i;
	workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];

	if(ll==1 && workv[ll][2][k][j][i]<0.0) workv[ll][2][k][j][i]=0.0;

	kk=k; jj=N2-1; ii=i;
	workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];

      }

      // x3 down
      for(k=-N3BND;k<=0;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	kk=1; jj=j; ii=i;
	workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];

	if(ll==1 && workv[ll][3][k][j][i]>0.0) workv[ll][3][k][j][i]=0.0;

	kk=0; jj=j; ii=i;
	workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];

      }
      // x3 up
      for(k=N3;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	kk=N3-1; jj=j; ii=i;
	workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];

	if(ll==1 && workv[ll][3][k][j][i]<0.0) workv[ll][3][k][j][i]=0.0;

	kk=N3-1; jj=j; ii=i;
	workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];


      }


    
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMVEC;
    }// end over vectors
  }// endif vectors to be done
  







#if(DOINTERNALBOUNDARY==1)
  // if doing MPI:
#if(USEMPI)
    bound_mpi(vars,varv,wsca,wvec,wcom);
#endif
#endif

  firsttime=0;
}// end function


