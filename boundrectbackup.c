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



  //  fprintf(stderr,"I got here %g\n",t); fflush(stderr);


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

	// OUTFLOW
	//	works[ll][k][j][i]=works[ll][kk][jj][ii];

	// Avery WIND
	if(ll!=0) works[ll][k][j][i]=sanal[ll][k][j][i];
	else works[ll][k][j][i]=works[ll][kk][jj][ii]; // just outflow passive scalar
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


      //      fprintf(stderr,"1 I got here %g l=%d ll=%d\n",t,l,ll); fflush(stderr);
      
	
      // x1 down
      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<=0;i++){

	// OUTFLOW
	kk=k; jj=j; ii=1;
	//	workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	//	if(ll==1 && workv[ll][1][k][j][i]>0.0) workv[ll][1][k][j][i]=0.0;
	kk=k; jj=j; ii=0;
	//	workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	//	workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];


	// Avery WIND
	if(ll!=0){
	  workv[ll][1][k][j][i]=vanal[ll][1][k][j][i];
	  workv[ll][2][k][j][i]=vanal[ll][2][k][j][i];
	  workv[ll][3][k][j][i]=vanal[ll][3][k][j][i];
	}
	else{
	  // OUTFLOW
	  kk=k; jj=j; ii=1;
	  workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  kk=k; jj=j; ii=0;
	  workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	  workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];
	}

      }
      // x1 up
      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=N1;i<N1+N1BND;i++){
	kk=k; jj=j; ii=N1-1;
	workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];

	if(ll==1 && workv[ll][1][k][j][i]<0.0) workv[ll][1][k][j][i]=0.0;

	kk=k; jj=j; ii=N1-1;
	workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];

	//	fprintf(stderr,"2 I got here %g\n",t); fflush(stderr);
	

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











#define LOOPFBACK3 for(k=N3+N3BND-1;k>=-N3BND;k--)
#define LOOPFBACK2 for(j=N2+N2BND-1;j>=-N2BND;j--)
#define LOOPFBACK1 for(i=N1+N1BND-1;i>=-N1BND;i--)

#define LOOPFBACK LOOPFBACK3 LOOPFBACK2 LOOPFBACK1



#define DOSCALARACCRETOR 1
#define DOVECTORACCRETOR 1




void bound_accretor(int dir,FTYPE (*vars)[N2M][N1M],
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
  int vdir;
  int ri,rj,rk;



  numhit++;
  if(wsca<=-2){
    works=(FTYPE (*) [N3M][N2M][N1M])(&vars[0][0][0]);
  }
  else works=s;
  if(wvec<=-2){
    workv=(FTYPE (*) [3][N3M][N2M][N1M])(&varv[0][0][0][0]);
  }
  else workv=v;


#if(!DO_ACCRETOR)
  // turn on/off accretor
  return;
#endif

  /////////////////////////////////////
  /////////////////////////////////////
  ////////////// SCALARS
  /////////////////////////////////////
  /////////////////////////////////////


  if(wsca!=0){

    for(l=1;l<=NUMSCA;l++){ // bounds potential
    
      if(wsca!=-1){ if(wsca<=-2) ll=0; else ll=wsca; }
      else ll=l;


      // create ad hoc values for inside accretor
      LOOPF{
	if(inside_accretor(CENT,i,j,k)){
	  if(ll==1)  works[ll][k][j][i]=1.0; // mass density
	  else if(ll==2)  works[ll][k][j][i]=0.1; // internal energy density (low sound speed)
	  else if(ll==3)  works[ll][k][j][i]=0.0; // potential (not used)
	}
      }

      // create ad hoc values for inside accretor
      LOOPF{
	if(inside_accretor_general(3,CENT,i,j,k)){ // # zones for passive guy
	  if( (ll==0) && (wsca==-3) )  works[ll][k][j][i]=1.0; // passive scalar
	}
      }



      ////////
      // Treating passive scalar in same way as other scalars causes passive scalar is removed as copying from real domain
      //
      // Try to just set passive scalar and return
      if( (ll==0) && (wsca==-3) ){

      }
      else{ // only do for non-passive scalar

#if(DOSCALARACCRETOR)

	LOOPF{ // loop over all i,j,k

      
	  // bounding in the 1-direction
	  if(dir==1){
	    // if accretor position is RIGHT of the grid location at i, then copy from left INTO grid location
	    ri=i;
	    rj=j;
	    rk=k;
	    if( (inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][ri]<pos_accretor[dir]) ){
	      works[ll][rk][rj][ri] = works[ll][rk][rj][ri-1];
	    }
	    ri=N1-1-i;
	    rj=j;
	    rk=k;
	    if( (inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][ri]>pos_accretor[dir]) ){
	      works[ll][rk][rj][ri] = works[ll][rk][rj][ri+1];
	    }
	    ri=i;
	    rj=j;
	    rk=k;
	    if((inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][ri]==pos_accretor[dir]) ) works[ll][rk][rj][ri] = 0.5*(works[ll][rk][rj][ri+1]+works[ll][rk][rj][ri-1]);

	  }
	  // bounding in the 2-direction
	  else if(dir==2){
	    ri=i;
	    rj=j;
	    rk=k;
	    if( (inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][rj]<pos_accretor[dir]) ){
	      works[ll][rk][rj][ri] = works[ll][rk][rj-1][ri];
	    }
	    ri=i;
	    rj=N2-1-j;
	    rk=k;
	    if( (inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][rj]>pos_accretor[dir]) ){
	      works[ll][rk][rj][ri] = works[ll][rk][rj+1][i];
	    }
	    ri=i;
	    rj=j;
	    rk=k;
	    if((inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][rj]==pos_accretor[dir]) ) works[ll][rk][rj][ri] = 0.5*(works[ll][rk][rj+1][ri]+works[ll][rk][rj-1][ri]);
	  }
	  // bounding in the 3-direction
	  else if(dir==3){
	    ri=i;
	    rj=j;
	    rk=k;
	    if( (inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][rk]<pos_accretor[dir]) ){
	      works[ll][rk][rj][ri] = works[ll][rk-1][rj][ri];
	    }
	    ri=i;
	    rj=j;
	    rk=N3-1-k;
	    if( (inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][rk]>pos_accretor[dir]) ){
	      works[ll][rk][rj][ri] = works[ll][rk+1][rj][ri];
	    }
	    ri=i;
	    rj=j;
	    rk=k;
	    if((inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][rk]==pos_accretor[dir]) ) works[ll][rk][rj][ri] = 0.5*(works[ll][rk+1][rj][ri]+works[ll][rk-1][rj][ri]);
	  }
	}// end over all zones

    
#endif

	// bomb
	//if(ll==2) works[ll][16][16][8] = 1000.0;
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


#define VDIRLOOP for(vdir=1;vdir<=3;vdir++)

      // create ad hoc values for inside accretor
      LOOPF{
	VDIRLOOP{
	  if(inside_accretor(vdir,i,j,k)){
	    if(ll==0) workv[ll][vdir][k][j][i]=0.0; // other vectors (mdot, etc.)
	    else if(ll==1)  workv[ll][vdir][k][j][i]=0.0; // velocity
	    else if(ll==2)  workv[ll][vdir][k][j][i]=0.0; // field
	  }
	}
      }



#if(DOVECTORACCRETOR)

      // GODMARK
      // by doing "2" loops at once, I assume all relevant velocities are set right, which probably means accretor is not near real boundary

      LOOPF{ // loop over all i,j,k

	VDIRLOOP{ // loop over vector directions

	  // bounding in the 1-direction
	  if(dir==1){ // dir corresponds to which direction interpolation is going to be worked on and need to set BCs so that interpolation is like shifting stencil w/ DONOR

	    ri=i;
	    rj=j;
	    rk=k;
	    if((inside_accretor(vdir,ri,rj,rk))&&
	       (
		((x[1][dir][ri]<pos_accretor[dir])&&(vdir==1))||
		((x[2][dir][ri]<pos_accretor[dir])&&(vdir==2))||
		((x[2][dir][ri]<pos_accretor[dir])&&(vdir==3))
		)
	       ){
	      workv[ll][vdir][rk][rj][ri] = workv[ll][vdir][rk][rj][ri-1];
	      if(vdir==dir) if(workv[ll][vdir][rk][rj][ri]<0.0) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]+VHORIZON[vdir]; // choosing VHORIZON!=0 screws up effective shifting!
	      // also not general for all "ll"
	    }
	    ri=N1-1-i;
	    rj=j;
	    rk=k;
	    if((inside_accretor(vdir,ri,rj,rk))&&
	       (
		((x[1][dir][ri]>pos_accretor[dir])&&(vdir==1))||
		((x[2][dir][ri]>pos_accretor[dir])&&(vdir==2))||
		((x[2][dir][ri]>pos_accretor[dir])&&(vdir==3))
		)
	       ){
	      workv[ll][vdir][rk][rj][ri] = workv[ll][vdir][rk][rj][ri+1];
	      if(vdir==dir) if(workv[ll][vdir][rk][rj][ri]>0.0) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]-VHORIZON[vdir]; // choosing VHORIZON!=0 screws up effective shifting!
	      // also not general for all "ll"
	    }
	    ri=i;
	    rj=j;
	    rk=k;
	    if((inside_accretor(vdir,ri,rj,rk))&&
	       (
		((x[1][dir][ri]==pos_accretor[dir])&&(vdir==1))||
		((x[2][dir][ri]==pos_accretor[dir])&&(vdir==2))||
		((x[2][dir][ri]==pos_accretor[dir])&&(vdir==3))
		)
	       ){
	      workv[ll][vdir][rk][rj][ri] =  0.5*(workv[ll][vdir][rk][rj][ri+1]+workv[ll][vdir][rk][rj][ri-1]);
	      if(vdir==dir) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]; // or velocity of accretor itself! // GODMARK
	      // also not general for all "ll"
	    }
	  }
	  // bounding in the 2-direction
	  else if(dir==2){

	    ri=i;
	    rj=j;
	    rk=k;
	    if((inside_accretor(vdir,ri,rj,rk))&&
	       (
		((x[2][dir][rj]<pos_accretor[dir])&&(vdir==1))||
		((x[1][dir][rj]<pos_accretor[dir])&&(vdir==2))||
		((x[2][dir][rj]<pos_accretor[dir])&&(vdir==3))
		)
	       ){
	      workv[ll][vdir][rk][rj][ri] = workv[ll][vdir][rk][rj-1][ri];
	      if(vdir==dir) if(workv[ll][vdir][rk][rj][ri]<0.0) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]+VHORIZON[vdir]; // choosing VHORIZON!=0 screws up effective shifting!
	      // also not general for all "ll"
	    }
	    ri=i;
	    rj=N2-1-j;
	    rk=k;
	    if((inside_accretor(vdir,ri,rj,rk))&&
	       (
		((x[2][dir][rj]>pos_accretor[dir])&&(vdir==1))||
		((x[1][dir][rj]>pos_accretor[dir])&&(vdir==2))||
		((x[2][dir][rj]>pos_accretor[dir])&&(vdir==3))
		)
	       ){
	      workv[ll][vdir][rk][rj][ri] = workv[ll][vdir][rk][rj+1][ri];
	      if(vdir==dir) if(workv[ll][vdir][rk][rj][ri]>0.0) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]-VHORIZON[vdir]; // choosing VHORIZON!=0 screws up effective shifting!
	      // also not general for all "ll"
	    }
	    ri=i;
	    rj=j;
	    rk=k;
	    if((inside_accretor(vdir,ri,rj,rk))&&
	       (
		((x[2][dir][rj]==pos_accretor[dir])&&(vdir==1))||
		((x[1][dir][rj]==pos_accretor[dir])&&(vdir==2))||
		((x[2][dir][rj]==pos_accretor[dir])&&(vdir==3))
		)
	       ){
	      workv[ll][vdir][rk][rj][ri] =  0.5*(workv[ll][vdir][rk][rj+1][ri]+workv[ll][vdir][rk][rj-1][ri]);
	      if(vdir==dir) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]; // or velocity of accretor itself!
	      // also not general for all "ll"
	    }

	  }
	  // bounding in the 3-direction
	  else if(dir==3){

	    ri=i;
	    rj=j;
	    rk=k;
	    if((inside_accretor(vdir,ri,rj,rk))&&
	       (
		((x[2][dir][rk]<pos_accretor[dir])&&(vdir==1))||
		((x[2][dir][rk]<pos_accretor[dir])&&(vdir==2))||
		((x[1][dir][rk]<pos_accretor[dir])&&(vdir==3))
		)
	       ){
	      workv[ll][vdir][rk][rj][ri] = workv[ll][vdir][rk-1][rj][ri];
	      if(vdir==dir) if(workv[ll][vdir][rk][rj][ri]<0.0) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]+VHORIZON[vdir]; // choosing VHORIZON!=0 screws up effective shifting!
	      // also not general for all "ll"
	    }
	    ri=i;
	    rj=j;
	    rk=N3-1-k;
	    if((inside_accretor(vdir,ri,rj,rk))&&
	       (
		((x[2][dir][rk]>pos_accretor[dir])&&(vdir==1))||
		((x[2][dir][rk]>pos_accretor[dir])&&(vdir==2))||
		((x[1][dir][rk]>pos_accretor[dir])&&(vdir==3))
		)
	       ){
	      workv[ll][vdir][rk][rj][ri] = workv[ll][vdir][rk+1][rj][ri];
	      if(vdir==dir) if(workv[ll][vdir][rk][rj][ri]>0.0) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]-VHORIZON[vdir]; // choosing VHORIZON!=0 screws up effective shifting!
	      // also not general for all "ll"
	    }
	    ri=i;
	    rj=j;
	    rk=k;
	    if((inside_accretor(vdir,ri,rj,rk))&&
	       (
		((x[2][dir][rk]==pos_accretor[dir])&&(vdir==1))||
		((x[2][dir][rk]==pos_accretor[dir])&&(vdir==2))||
		((x[1][dir][rk]==pos_accretor[dir])&&(vdir==3))
		)
	       ){
	      workv[ll][vdir][rk][rj][ri] =  0.5*(workv[ll][vdir][rk+1][rj][ri]+workv[ll][vdir][rk-1][rj][ri]);
	      if(vdir==dir) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]; // or velocity of accretor itself!
	      // also not general for all "ll"
	    }


	  }
	} // end over vector directions
      }// end over all zones


#endif
    
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMVEC;
    }// end over vectors
  }// endif vectors to be done
  






}


int inside_accretor(int which, int i, int j, int k)
{
  int inside_accretor_general(int numzones, int which, int i, int j, int k);
  
  return(inside_accretor_general(-1, which, i, j, k));

}


// which = CENT, VDIR1, VDIR2, VDIR3
int inside_accretor_general(int numzones, int which, int i, int j, int k)
{
  FTYPE local_DX_accretor;

#if(DO_ACCRETOR)
  
  if(numzones<0) local_DX_accretor=DX_accretor;
  else local_DX_accretor=((FTYPE)numzones)*DX_single;

  if(which==CENT){
    if(
       (fabs(x[2][1][i]-pos_accretor[1])<local_DX_accretor)&&
       (fabs(x[2][2][j]-pos_accretor[2])<local_DX_accretor)&&
       (fabs(x[2][3][k]-pos_accretor[3])<local_DX_accretor)
       ){
      return(1);
    }
    else return(0);
  }
  else if(which==VDIR1){
    if(
       (fabs(x[1][1][i]-pos_accretor[1])<local_DX_accretor)&&
       (fabs(x[2][2][j]-pos_accretor[2])<local_DX_accretor)&&
       (fabs(x[2][3][k]-pos_accretor[3])<local_DX_accretor)
       ){
      return(1);
    }
    else return(0);
  }
  else if(which==VDIR2){
    if(
       (fabs(x[2][1][i]-pos_accretor[1])<local_DX_accretor)&&
       (fabs(x[1][2][j]-pos_accretor[2])<local_DX_accretor)&&
       (fabs(x[2][3][k]-pos_accretor[3])<local_DX_accretor)
       ){
      return(1);
    }
    else return(0);
  }
  else if(which==VDIR3){
    if(
       (fabs(x[2][1][i]-pos_accretor[1])<local_DX_accretor)&&
       (fabs(x[2][2][j]-pos_accretor[2])<local_DX_accretor)&&
       (fabs(x[1][3][k]-pos_accretor[3])<local_DX_accretor)
       ){
      return(1);
    }
    else return(0);
  }
  else{
    fprintf(stderr,"No such which=%d\n",which);
    exit(1);
  }
  
#else

  return(0); // assume 0 means no change, nothing special to be done

#endif


}

