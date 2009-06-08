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
  int offset1,offset2,offset3;


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


  if(skipix1==1) offset1=0;
  else offset1=-1;

  if(skipix2==1) offset2=0;
  else offset2=-1;

  if(skipix3==1) offset3=0;
  else offset3=-1;




  /////////////////////////////////////
  /////////////////////////////////////
  ////////////// SCALARS
  /////////////////////////////////////
  /////////////////////////////////////

#define LIMITCOPY {if(kk>N3-1) kk=N3-1; else if(kk<0) kk=0; if(jj>N2-1) jj=N2-1; else if(jj<0) jj=0; if(ii>N1-1) ii=N1-1; else if(ii<0) ii=0;}
 

  if(wsca!=0){

    for(l=1;l<=NUMSCA;l++){ // bounds potential
    
      if(wsca!=-1){ if(wsca<=-2) ll=0; else ll=wsca; }
      else ll=l;
 

      
      
      // x1 down
      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<0;i++){
	kk=k; jj=j; ii=0;

	LIMITCOPY;

	// OUTFLOW
	//	works[ll][k][j][i]=works[ll][kk][jj][ii];

	// Avery WIND
	if(ll!=0) works[ll][k][j][i]=sanal[ll][k][j][i];
	else works[ll][k][j][i]=works[ll][kk][jj][ii]; // just outflow passive scalar
      }
      // x1 up
      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=N1;i<N1+N1BND;i++){
	kk=k; jj=j; ii=N1-1;

	LIMITCOPY;
	
	works[ll][k][j][i]=works[ll][kk][jj][ii];

	//	if(ll==1) if(works[ll][k][j][i]<1E-6){
	//	  fprintf(log_file,"Problem at i=%d j=%d k=%d %21.15g\n",i,j,k,works[ll][k][j][i]);
	//	}
      }

      // x2 down
      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<0;j++) for(i=-N1BND;i<N1+N1BND;i++){
	kk=k; jj=0; ii=i;

	LIMITCOPY;

	works[ll][k][j][i]=works[ll][kk][jj][ii];
      }
      // x2 up
      for(k=-N3BND;k<N3+N3BND;k++) for(j=N2;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	kk=k; jj=N2-1; ii=i;

	LIMITCOPY;

	works[ll][k][j][i]=works[ll][kk][jj][ii];
      }

      // x3 down
      for(k=-N3BND;k<0;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	kk=0; jj=j; ii=i;

	LIMITCOPY;

	works[ll][k][j][i]=works[ll][kk][jj][ii];
      }
      // x3 up
      for(k=N3;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	kk=N3-1; jj=j; ii=i;

	LIMITCOPY;

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
  if(wvec!=0 && wvec!=-3){

    // GODMARK: wvec==-3 corresponds to mdot, which is only bounded for periodic or MPI... not normal quantity
    

    for(l=1;l<=REALNUMVEC;l++){
      
      /* if not to do all, pick */
      if(wvec!=-1){ if(wvec<=-2) ll=0; else ll=wvec; }
      else ll=l;


      //      fprintf(stderr,"1 I got here %g l=%d ll=%d\n",t,l,ll); fflush(stderr);
      
	
      // x1 down
      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<offset1;i++){

	// OUTFLOW
	kk=k; jj=j; ii=1;
	LIMITCOPY;
	//	workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	//	if(ll==1 && workv[ll][1][k][j][i]>0.0) workv[ll][1][k][j][i]=0.0;
	kk=k; jj=j; ii=0;
	LIMITCOPY;
	//	workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	//	workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];



	// Avery WIND
	if(ll!=0){
	  workv[ll][1][k][j][i]=vanal[ll][1][k][j][i];
	  workv[ll][2][k][j][i]=vanal[ll][2][k][j][i];
	  workv[ll][3][k][j][i]=vanal[ll][3][k][j][i];
	}
	else{
	  if(wvec!=-4){
	    // OUTFLOW
	    kk=k; jj=j; ii=1;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	    if(ll==1 && workv[ll][1][k][j][i]>0.0) workv[ll][1][k][j][i]=0.0;
	    
	    kk=k; jj=j; ii=0;
	    LIMITCOPY;
	    workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	    workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];
	  }
	  else{
	    // OUTFLOW
	    if(passivepos==VDIR1){
	      kk=k; jj=j; ii=1;
	      LIMITCOPY;
	      workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	    }
	    else if((passivepos==VDIR2)||(passivepos==VDIR3)){
	      kk=k; jj=j; ii=0;
	      LIMITCOPY;
	      workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	    }
	    else if(passivepos==-1){
	      kk=k; jj=j; ii=1;
	      LIMITCOPY;
	      workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	      kk=k; jj=j; ii=0;
	      LIMITCOPY;
	      workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	      workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];

	    }
	  }
	}

      }
      // x1 up
      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=N1;i<N1+N1BND;i++){

	if(wvec!=-4){
	  kk=k; jj=j; ii=N1-1;
	  LIMITCOPY;
	  //	  if(ll==1 || wvec>=-1){
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	    //	  }
	    //	  else{
	    //	    fprintf(log_file,"wvec=%d\n",wvec);
	    //	  }

	  //	  fprintf(log_file,"workv[%d][1][%d][%d][%d]=%21.15g\n",ll,k,j,i,workv[ll][1][k][j][i]);
	  
	  if(ll==1 && workv[ll][1][k][j][i]<0.0) workv[ll][1][k][j][i]=0.0;
	  
	  kk=k; jj=j; ii=N1-1;
	  LIMITCOPY;
	  workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	  workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];
	}
	else{
	  if(passivepos==VDIR1){
	    kk=k; jj=j; ii=N1-1;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  }
	  else if((passivepos==VDIR2)||(passivepos==VDIR3)){
	    kk=k; jj=j; ii=N1-1;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  }
	  else if(passivepos==-1){
	    kk=k; jj=j; ii=N1-1;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	    kk=k; jj=j; ii=N1-1;
	    LIMITCOPY;
	    workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	    workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];
	  }
	}

	//	fprintf(stderr,"2 I got here %g\n",t); fflush(stderr);
	
	//	if(myid==1) fprintf(log_file,"workv[%d][1][%d][%d][%d]=%21.15g\n",ll,k,j,i,workv[ll][1][k][j][i]);


      }

      // x2 down
      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<offset2;j++) for(i=-N1BND;i<N1+N1BND;i++){

	if(wvec!=-4){
	  
	  kk=k; jj=1; ii=i;
	  LIMITCOPY;
	  workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	  
	  if(ll==1 && workv[ll][2][k][j][i]>0.0) workv[ll][2][k][j][i]=0.0;
	  
	  kk=k; jj=0; ii=i;
	  LIMITCOPY;
	  // BAD?
	  //	  workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];
	  
	}
	else{
	  if(passivepos==VDIR2){
	    kk=k; jj=1; ii=i;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  }
	  else if((passivepos==VDIR1)||(passivepos==VDIR3)){
	    kk=k; jj=0; ii=i;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  }
	  else if(passivepos==-1){
	    kk=k; jj=1; ii=i;
	    LIMITCOPY;
	    workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	    kk=k; jj=0; ii=i;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	    workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];
	  }
	}
      }
      // x2 up
      for(k=-N3BND;k<N3+N3BND;k++) for(j=N2;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	
	if(wvec!=-4){
	  kk=k; jj=N2-1; ii=i;
	  LIMITCOPY;
	  workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	  
	  if(ll==1 && workv[ll][2][k][j][i]<0.0) workv[ll][2][k][j][i]=0.0;
	  
	  kk=k; jj=N2-1; ii=i;
	  LIMITCOPY;
	  // BAD?
	  //	  workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];
	}
	else{
	  if(passivepos==VDIR2){
	    kk=k; jj=N2-1; ii=i;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  }
	  else if((passivepos==VDIR1)||(passivepos==VDIR3)){
	    kk=k; jj=N2-1; ii=i;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  }
	  else if(passivepos==-1){
	    kk=k; jj=N2-1; ii=i;
	    LIMITCOPY;
	    workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	    kk=k; jj=N2-1; ii=i;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	    workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];
	  }
	}
      }

      // x3 down
      for(k=-N3BND;k<offset3;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){

	if(wvec!=-4){
	  kk=1; jj=j; ii=i;
	  LIMITCOPY;
	  workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];
	  
	  if(ll==1 && workv[ll][3][k][j][i]>0.0) workv[ll][3][k][j][i]=0.0;
	  
	  kk=0; jj=j; ii=i;
	  LIMITCOPY;
	  // BAD?
	  //	  workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	}
	else{
	  if(passivepos==VDIR3){
	    kk=1; jj=j; ii=i;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  }
	  else if((passivepos==VDIR1)||(passivepos==VDIR2)){
	    kk=0; jj=j; ii=i;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  }
	  else if(passivepos==-1){
	    kk=1; jj=j; ii=i;
	    LIMITCOPY;
	    workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];
	    kk=0; jj=j; ii=i;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	    workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	  }
	}
      }
      // x3 up
      for(k=N3;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	
	if(wvec!=-4){
	  kk=N3-1; jj=j; ii=i;
	  LIMITCOPY;
	  workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];
	  
	  if(ll==1 && workv[ll][3][k][j][i]<0.0) workv[ll][3][k][j][i]=0.0;
	  
	  kk=N3-1; jj=j; ii=i;
	  LIMITCOPY;
	  // BAD?
	  //	  workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	}
	else{
	  if(passivepos==VDIR3){
	    kk=N3-1; jj=j; ii=i;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  }
	  else if((passivepos==VDIR1)||(passivepos==VDIR2)){
	    kk=N3-1; jj=j; ii=i;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	  }
	  else if(passivepos==-1){
	    kk=N3-1; jj=j; ii=i;
	    LIMITCOPY;
	    workv[ll][3][k][j][i]=workv[ll][3][kk][jj][ii];
	    kk=N3-1; jj=j; ii=i;
	    LIMITCOPY;
	    workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
	    workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
	  }
	}
	

      }


    
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMVEC;
    }// end over vectors
  }// endif vectors to be done
  

  /*
  LOOPF{
    if(myid==0){
      if(((i==N1)||(i==N1+1))&&(k==N3/2)){
	fprintf(log_file,"BEFORE   rho=%g u=%g v1=%g v2=%g v3=%g\n",s[1][k][j][i],s[2][k][j][i],v[1][1][k][j][i],v[1][2][k][j][i],v[1][3][k][j][i]);
      }
    }
  }

  LOOPF{
    if(myid==1){
      if(((i==-1)||(i==-2))&&(k==N3/2)){
	fprintf(log_file,"BEFORE   rho=%g u=%g v1=%g v2=%g v3=%g\n",s[1][k][j][i],s[2][k][j][i],v[1][1][k][j][i],v[1][2][k][j][i],v[1][3][k][j][i]);
      }
    }
  }
  */
  
  
  
#if(DOINTERNALBOUNDARY==1)
  // if doing MPI:
#if(USEMPI)
  bound_mpi(vars,varv,wsca,wvec,wcom);
#endif
#endif
  
  
  /*
  LOOPF{
    if(myid==0){
      if(((i==N1)||(i==N1+1))&&(k==N3/2)){
	fprintf(log_file,"AFTER    rho=%g u=%g v1=%g v2=%g v3=%g\n",s[1][k][j][i],s[2][k][j][i],v[1][1][k][j][i],v[1][2][k][j][i],v[1][3][k][j][i]);
      }
    }
  }

  LOOPF{
    if(myid==1){
      if(((i==-1)||(i==-2))&&(k==N3/2)){
	fprintf(log_file,"AFTER   rho=%g u=%g v1=%g v2=%g v3=%g\n",s[1][k][j][i],s[2][k][j][i],v[1][1][k][j][i],v[1][2][k][j][i],v[1][3][k][j][i]);
      }
    }
  }
  */

  firsttime=0;
}// end function











#define LOOPFBACK3 for(k=N3+N3BND-1;k>=-N3BND;k--)
#define LOOPFBACK2 for(j=N2+N2BND-1;j>=-N2BND;j--)
#define LOOPFBACK1 for(i=N1+N1BND-1;i>=-N1BND;i--)

#define LOOPFBACK LOOPFBACK3 LOOPFBACK2 LOOPFBACK1


  // can extend the below into the GHOST cells GODMARK
#define LIMITCOPYACCRETOR {if(kk>N3-1){ kk=N3-1; beyondaccretor=1;} else if(kk<0){ kk=0; beyondaccretor=1;} if(jj>N2-1){ jj=N2-1; beyondaccretor=1;} else if(jj<0){ jj=0; beyondaccretor=1;} if(ii>N1-1){ ii=N1-1;  beyondaccretor=1;} else if(ii<0){ ii=0;  beyondaccretor=1;}}


#define DOSCALARACCRETOR 1
#define DOVECTORACCRETOR 1

//#define NUMPASSIVEZONES 3
#define NUMPASSIVEZONES (NUMZONE_accretor)


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
  FTYPE ftemp2,cs2;
  int beyondaccretor;

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

  // assume not copying beyond accretor
  beyondaccretor=0;



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


      if( (ll==0) && (wsca==-3) ){
	if(passivepos==CENT){
	  // create ad hoc values for inside accretor
	  LOOPF{
	    if(inside_accretor_general(NUMPASSIVEZONES,passivepos,i,j,k)){ // # zones for passive guy
	      works[ll][k][j][i]=1.0; // passive scalar
	      //	      fprintf(stderr,"got here i=%d j=%d k=%d v[2][1]=%g\n",i,j,k,v[2][1][k][j][i]);
	    }
	  }
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
	      kk=rk;
	      jj=rj;
	      ii=ri-1;
	      LIMITCOPYACCRETOR;
	      works[ll][rk][rj][ri] = works[ll][kk][jj][ii];
	    }
	    ri=N1-1-i;
	    rj=j;
	    rk=k;
	    if( (inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][ri]>pos_accretor[dir]) ){
	      kk=rk;
	      jj=rj;
	      ii=ri+1;
	      LIMITCOPYACCRETOR;
	      works[ll][rk][rj][ri] = works[ll][kk][jj][ii];
	    }
	    ri=i;
	    rj=j;
	    rk=k;
	    if((inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][ri]==pos_accretor[dir]) ){
	      kk=rk;
	      jj=rj;
	      ii=ri+1;
	      LIMITCOPYACCRETOR;

	      works[ll][rk][rj][ri] += 0.5*works[ll][kk][jj][ii];

	      kk=rk;
	      jj=rj;
	      ii=ri-1;
	      LIMITCOPYACCRETOR;
	      works[ll][rk][rj][ri] += 0.5*works[ll][kk][jj][ii];
	    }

	  }
	  // bounding in the 2-direction
	  else if(dir==2){
	    ri=i;
	    rj=j;
	    rk=k;
	    if( (inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][rj]<pos_accretor[dir]) ){
	      kk=rk;
	      jj=rj-1;
	      ii=ri;
	      LIMITCOPYACCRETOR;
	      works[ll][rk][rj][ri] = works[ll][kk][jj][ii];
	    }
	    ri=i;
	    rj=N2-1-j;
	    rk=k;
	    if( (inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][rj]>pos_accretor[dir]) ){
	      kk=rk;
	      jj=rj+1;
	      ii=ri;
	      LIMITCOPYACCRETOR;
	      works[ll][rk][rj][ri] = works[ll][kk][jj][ii]; // was i instead of ri before!
	    }
	    ri=i;
	    rj=j;
	    rk=k;
	    if((inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][rj]==pos_accretor[dir]) ){
	      kk=rk;
	      jj=rj+1;
	      ii=ri;
	      LIMITCOPYACCRETOR;

	      works[ll][rk][rj][ri] += 0.5*works[ll][kk][jj][ii];

	      kk=rk;
	      jj=rj-1;
	      ii=ri;
	      LIMITCOPYACCRETOR;

	      works[ll][rk][rj][ri] += 0.5*works[ll][kk][jj][ii];

	    }
	  }
	  // bounding in the 3-direction
	  else if(dir==3){
	    ri=i;
	    rj=j;
	    rk=k;
	    if( (inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][rk]<pos_accretor[dir]) ){
	      kk=rk-1;
	      jj=rj;
	      ii=ri;
	      LIMITCOPYACCRETOR;
	      works[ll][rk][rj][ri] = works[ll][kk][jj][ii];
	    }
	    ri=i;
	    rj=j;
	    rk=N3-1-k;
	    if( (inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][rk]>pos_accretor[dir]) ){
	      kk=rk+1;
	      jj=rj;
	      ii=ri;
	      LIMITCOPYACCRETOR;
	      works[ll][rk][rj][ri] = works[ll][kk][jj][ii];
	    }
	    ri=i;
	    rj=j;
	    rk=k;
	    if((inside_accretor(CENT,ri,rj,rk)) && (x[2][dir][rk]==pos_accretor[dir]) ){
	      kk=rk+1;
	      jj=rj;
	      ii=ri;
	      LIMITCOPYACCRETOR;

	      works[ll][rk][rj][ri] += 0.5*works[ll][kk][jj][ii];

	      kk=rk-1;
	      jj=rj;
	      ii=ri;
	      LIMITCOPYACCRETOR;

	      works[ll][rk][rj][ri] += 0.5*works[ll][kk][jj][ii];

	    }
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
  

  



  if(beyondaccretor){
    fprintf(fail_file,"SCA: copying beyond accretor\n");
    myexit(6262);
  }







  ////////////////////////////////////////////////
  /////////////////////////////////////////////////
  ////////////////////////////////////////////////
  //////////////////////////// VECTORS
  ///////////////////////////////////////////////
  /////////////////////////////////////////////////
  ////////////////////////////////////////////////
  // Now do vectors if any
  if(wvec!=0 && wvec!=-3){
    
    // GODMARK: wvec==-3 corresponds to mdot, which is only bounded for periodic or MPI... not normal quantity


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


      // create ad hoc values for inside accretor
      if(ll==0  && wvec==-4){
	if(passivepos!=CENT && passivepos!=-1){
	  vdir=1; // passive scalar is using x-dir position.  passivepos is position of this quantity
	  LOOPF{
	    //	VDIRLOOP{
	    if(inside_accretor_general(NUMPASSIVEZONES,passivepos,i,j,k)){
	      workv[ll][vdir][k][j][i]=1.0; // passive scalar at vector position
	    }
	    //	}
	  }
	}
	else if(passivepos==-1){
	  LOOPF{
	    VDIRLOOP{
	      if(inside_accretor_general(NUMPASSIVEZONES,vdir,i,j,k)){
		workv[ll][vdir][k][j][i]=1.0; // passive scalar at vector position
	      }
	    }
	  }
	}
      }



#if(DOVECTORACCRETOR)


      if(ll==0  && wvec==-4){

      }
      else{

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
		kk=rk;
		jj=rj;
		ii=ri-1;
		LIMITCOPYACCRETOR;

		workv[ll][vdir][rk][rj][ri] = workv[ll][vdir][kk][jj][ii];
		if(vdir==dir && ll==1){

		  if(MACHHORIZON>=0){
		    ftemp2=gam*(gam-1.);
		    if(wgam) cs2 = ftemp2*s[2][kk][jj][ii]/s[1][kk][jj][ii]; // assumes scalars are set first!
		    else cs2=cs*cs;
		    VHORIZON[vdir]=MACHHORIZON*sqrt(cs2);
		  }

		  if(workv[ll][vdir][rk][rj][ri]<VACCRETOR[vdir]+VHORIZON[vdir]) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]+VHORIZON[vdir]; // choosing VHORIZON!=0 screws up effective shifting!
		}
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
		kk=rk;
		jj=rj;
		ii=ri+1;
		LIMITCOPYACCRETOR;

		workv[ll][vdir][rk][rj][ri] = workv[ll][vdir][kk][jj][ii];
		if(vdir==dir && ll==1){

		  if(MACHHORIZON>=0){
		    ftemp2=gam*(gam-1.);
		    if(wgam) cs2 = ftemp2*s[2][kk][jj][ii]/s[1][kk][jj][ii];
		    else cs2=cs*cs;
		    VHORIZON[vdir]=MACHHORIZON*sqrt(cs2);
		  }
		  
		  if(workv[ll][vdir][rk][rj][ri]>VACCRETOR[vdir]-VHORIZON[vdir]) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]-VHORIZON[vdir]; // choosing VHORIZON!=0 screws up effective shifting!
		}
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
		kk=rk;
		jj=rj;
		ii=ri+1;
		LIMITCOPYACCRETOR;

		workv[ll][vdir][rk][rj][ri] +=  0.5*workv[ll][vdir][kk][jj][ii];

		kk=rk;
		jj=rj;
		ii=ri-1;
		LIMITCOPYACCRETOR;

		workv[ll][vdir][rk][rj][ri] +=  0.5*workv[ll][vdir][kk][jj][ii];

		if(vdir==dir && ll==1) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]; // or velocity of accretor itself! // GODMARK
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
		kk=rk;
		jj=rj-1;
		ii=ri;
		LIMITCOPYACCRETOR;

		workv[ll][vdir][rk][rj][ri] = workv[ll][vdir][kk][jj][ii];
		if(vdir==dir && ll==1){

		  if(MACHHORIZON>=0){
		    ftemp2=gam*(gam-1.);
		    if(wgam) cs2 = ftemp2*s[2][kk][jj][ii]/s[1][kk][jj][ii];
		    else cs2=cs*cs;
		    VHORIZON[vdir]=MACHHORIZON*sqrt(cs2);
		  }
		  
		  if(workv[ll][vdir][rk][rj][ri]<VACCRETOR[vdir]+VHORIZON[vdir]) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]+VHORIZON[vdir]; // choosing VHORIZON!=0 screws up effective shifting!
		}
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
		kk=rk;
		jj=rj+1;
		ii=ri;
		LIMITCOPYACCRETOR;

		workv[ll][vdir][rk][rj][ri] = workv[ll][vdir][kk][jj][ii];
		if(vdir==dir && ll==1){

		  if(MACHHORIZON>=0){
		    ftemp2=gam*(gam-1.);
		    if(wgam) cs2 = ftemp2*s[2][kk][jj][ii]/s[1][kk][jj][ii];
		    else cs2=cs*cs;
		    VHORIZON[vdir]=MACHHORIZON*sqrt(cs2);
		  }
		  
		  if(workv[ll][vdir][rk][rj][ri]>VACCRETOR[vdir]-VHORIZON[vdir]) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]-VHORIZON[vdir]; // choosing VHORIZON!=0 screws up effective shifting!
		}
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
		kk=rk;
		jj=rj+1;
		ii=ri;
		LIMITCOPYACCRETOR;

		workv[ll][vdir][rk][rj][ri] +=  0.5*workv[ll][vdir][kk][jj][ii];

		kk=rk;
		jj=rj-1;
		ii=ri;
		LIMITCOPYACCRETOR;

		workv[ll][vdir][rk][rj][ri] +=  0.5*workv[ll][vdir][kk][jj][ii];

		if(vdir==dir && ll==1) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]; // or velocity of accretor itself!
		// also not general for all "ll"
	      }

	    }
	    // bounding in the 3-direction
	    else if(dir==3){

	      //	      if(vdir==1 || vdir==2){////////////////////////////////////////////////////// still bad
	      //	      if(1){

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
		kk=rk-1;
		jj=rj;
		ii=ri;
		LIMITCOPYACCRETOR;

		workv[ll][vdir][rk][rj][ri] = workv[ll][vdir][kk][jj][ii];
		if(vdir==dir && ll==1){

		  if(MACHHORIZON>=0){
		    ftemp2=gam*(gam-1.);
		    if(wgam) cs2 = ftemp2*s[2][kk][jj][ii]/s[1][kk][jj][ii];
		    else cs2=cs*cs;
		    VHORIZON[vdir]=MACHHORIZON*sqrt(cs2);
		  }
		  
		  if(workv[ll][vdir][rk][rj][ri]<VACCRETOR[vdir]+VHORIZON[vdir]) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]+VHORIZON[vdir]; // choosing VHORIZON!=0 screws up effective shifting!
		}
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
		kk=rk+1;
		jj=rj;
		ii=ri;
		LIMITCOPYACCRETOR;

		workv[ll][vdir][rk][rj][ri] = workv[ll][vdir][kk][jj][ii];
		if(vdir==dir && ll==1){

		  if(MACHHORIZON>=0){
		    ftemp2=gam*(gam-1.);
		    if(wgam) cs2 = ftemp2*s[2][kk][jj][ii]/s[1][kk][jj][ii];
		    else cs2=cs*cs;
		    VHORIZON[vdir]=MACHHORIZON*sqrt(cs2);
		  }
		  
		  if(workv[ll][vdir][rk][rj][ri]>VACCRETOR[vdir]-VHORIZON[vdir]) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]-VHORIZON[vdir]; // choosing VHORIZON!=0 screws up effective shifting!
		}
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
		kk=rk+1;
		jj=rj;
		ii=ri;
		LIMITCOPYACCRETOR;

		workv[ll][vdir][rk][rj][ri] +=  0.5*workv[ll][vdir][kk][jj][ii];

		kk=rk-1;
		jj=rj;
		ii=ri;
		LIMITCOPYACCRETOR;

		workv[ll][vdir][rk][rj][ri] +=  0.5*workv[ll][vdir][kk][jj][ii];

		if(vdir==dir && ll==1) workv[ll][vdir][rk][rj][ri]=VACCRETOR[vdir]; // or velocity of accretor itself!
		// also not general for all "ll"
	      }

	      //	      }

	    }
	  } // end over vector directions
	}// end over all zones
      } // end if normal vector components

#endif
    
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMVEC;
    }// end over vectors
  }// endif vectors to be done
  


  if(beyondaccretor){
    fprintf(fail_file,"VEC: copying beyond accretor\n");
    myexit(6263);
  }





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

