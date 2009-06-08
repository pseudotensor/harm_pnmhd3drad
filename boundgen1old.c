#include "bound.h"
#include "boundgen.h"


// this code is for completely general boundary condition

void bound_gen1(FTYPE (*vars)[N2M][N1M],
	   FTYPE (*varv)[N3M][N2M][N1M],
	   int wsca,
	   int wvec,
	   int wcom) // 1,2,3, 0: none, 123=all 12=12 13=13, 23=23 (not currently setup to do 13 or 23 since never needed)
{
  static int firsttime=1;
  int tagheader,tagsend,tagrecv;
  int othercpu;
  int jumpfactor;
  int i,j,k,l,m,p,q ;
  int i2,j2,k2;
  int ii,jj,kk,ll;
  int jje2,jjje2;
  int iio,jjo,kko,llo;
  int iio2,jjo2,kko2,llo2;
  int iir,jjr,kkr,llr;
  int iii,jjj,kkk,lll;
  int iii2,jjj2,kkk2,lll2;
  int bct,bcdim,bcdir,bcd1,bcd2;
  int bcdir1,bcdir2,bcdir3;
  int bi,bj,bk;
  int bireal,bjreal,bkreal;
  int diri,dirj;
  int startloop,endloop,loopvar;
  int looper,component;
  FTYPE slope;
  FTYPE (*works)[N3M][N2M][N1M];
  FTYPE (*workv)[3][N3M][N2M][N1M];
  static int numhit=0;
  int numhitmin,numhitmax;
  int wbound;
  FTYPE ftemp,ftempv[3+1],ftempv2[3+1];
  char temps[50];
  int tempslen=40;
  int tempi;
  int forces=0, forcev=0, forcecheck;
  int docom[3+1], comlength,comorder[3+1],itemp;

  int dosign;

  int dualloop,divbfix;
  int doingscalars,doingvectors;

  int gval;
  int typeofbc;
  int loopi,loopj,loopk;

  int othercomp;
  int cornerto,corner;
  int edgeto, edge;

  int whichzone;
  int iisum,jjsum,kksum;
  int ininflag,outinflag,innerreflectflag;

  int doinner,doouter;
  int looperstart,looperend;

  // might not want to bound vectors if not advecting them
  //  if(wvec==-1) return;
  //  if(wvec!=0) return;
  //if(wsca!=0) return;
  // GODMARK(commented)
  //return;
  //  wsca=0;
  //  wvec=0;

  numhit++;
  if(wsca<=-2){
    works=(FTYPE (*) [N3M][N2M][N1M])(&vars[0][0][0]);
  }
  else works=s;
  if(wvec<=-2){
    workv=(FTYPE (*) [3][N3M][N2M][N1M])(&varv[0][0][0][0]);
  }
  else workv=v;


  if(firsttime||(!(nstep%NUMOUTERSKIP)) ) doouter=1; // only do firsttime and every NUMOUTERSKIPth time
  else doouter=0;

  doinner=1; // always do it

  /////////////////////////////////////
  /////////////////////////////////////
  ////////////// SCALARS
  /////////////////////////////////////
  /////////////////////////////////////

#if(DOBOUNDSCA)
  if(wsca!=0){
    doingvectors=0;
    doingscalars=1;

    for(l=1;l<=NUMSCA-NOBOUNDPOT;l++){

      /* if not to do all, pick */
      if(wsca!=-1){ if(wsca<=-2) ll=0; else ll=wsca; }
      else ll=l;
      
      // determine how things are bounded for other vars
      if(ll==0){
	if(wsca==-2){//other
	  // bound like rho
	  wbound=1;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0 for scalar case: wsca=%d\n",wsca);
	  myexit(1);
	}
      }
      else wbound=ll;
      
      
      looperstart=6;
      looperend=7;

      for(looper=looperstart;looper<=looperend;looper++){
	if((!doinner)&&(looper==6)) continue;
	if((!doouter)&&(looper==7)) continue;

	LOOPBOUND(looper){

	  bct=bzmask[k][j][i];

	// GODMARK
	//if( (bct>0)&&(bct<10) ){ // skip zone if computational zone or null zone (can get rid of this when satisfied that bound loop is correct)// GODMARK

	  switch(bct){
	    
	    // which zone determined a priori
	  case 4:
	  case 1:
	  case 2:
	  case 5: // assume if periodic done reasonably
	    // find average value for symmetric projection
	    itemp=bzs[looper][temptempi][0];
	    ftemp=0;
	    for(m=0;m<itemp;m++) ftemp+=works[ll][bzs[looper][temptempi][m*3+2+1]][bzs[looper][temptempi][m*3+1+1]][bzs[looper][temptempi][m*3+0+1]];
	    // now assign value
	    works[ll][k][j][i]=ftemp/((FTYPE)(itemp));
	    break;
	  case 3:
	    works[ll][k][j][i]=sanal[ll][k][j][i];
	    break;
	  case 99:
	  case 98:
	    // do nothing
	    break;
	  default:
	    fprintf(fail_file,"bound.c: switch(bct) error, case given is: %d\n\r",bct);
	    myexit(1);
	  } // end switch
	  //	} // end if bct>0
#if(DEBUG)
	if(bct<=0){
	  fprintf(fail_file,"bound.c: never should have bct!>0: %d\n\r",bct);
	  myexit(1);
	}
#endif
      }//end over current scalar boundary zones(minus corners)

      /* cut short loop if only to do one */
      if(wsca!=-1) l=NUMSCA;
    } // end over scalars
  } // end if any scalars

#endif // end if bounding scalars
  








  ////////////////////////////////////////////////
  /////////////////////////////////////////////////
  ////////////////////////////////////////////////
  //////////////////////////// VECTORS
  ///////////////////////////////////////////////
  /////////////////////////////////////////////////
  ////////////////////////////////////////////////
#if(DOBOUNDVEC)  
  // Now do vectors if any
  if(wvec!=0){
    doingvectors=1;
    doingscalars=0;
    
    
    // determine which components to do, if any
    if( (wcom==1)||(wcom==12)||(wcom==13)||(wcom==123) ) docom[1]=1; else docom[1]=0;
    if( (wcom==2)||(wcom==12)||(wcom==23)||(wcom==123) ) docom[2]=1; else docom[2]=0;
    if( (wcom==3)||(wcom==13)||(wcom==23)||(wcom==123) ) docom[3]=1; else docom[3]=0;
    if(wcom==0){ docom[1]=docom[2]=docom[3]=0; }
  }
#endif


  // repeats inflow code at inner edge for no technical reason except ease of code
#if(DOBOUNDVEC)  
  // Now do vectors if any
  if(wvec!=0){
    doingvectors=1;
    doingscalars=0;
    
    for(l=1;l<=REALNUMBOUNDVEC;l++){
      
      /* if not to do all, pick */
      if(wvec!=-1){ if(wvec<=-2) ll=0; else ll=wvec; }
      else ll=l;
      
      if(ll==0){
	if((wvec==-3)||(wvec==-4)){//mdot or flux(2) of v[1][2]
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-5){//mdot or flux(2) of v[1][2] or magemf
	  // bound like emf
	  wbound=3;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0, wvec=%d\n",wvec);
	  myexit(1);
	}
      }
      else wbound=ll;



      if(docom[1]){
	LOOPBOUNDV1{
	  whichzone=temptempi;
	  bct=maskv1[k][j][i];
	  
	  //if( (bct>0)&&(bct<10) ){ // skip zone if computational zone or null zone (can get rid of this when satisfied that bound loop is correct)// GODMARK

	    
	    switch(bct){
	      
	    case 1:
	    case 2:
	      /* reflect and AOS */
	      // order important for if/elseif's
	      
	      // find average value
	      ftemp=0;
	      innerreflectflag=0;
	      for(m=0;m<bzv1[whichzone][0];m++){
		
		ii=bzv1[whichzone][m*3+0+1];
		jj=bzv1[whichzone][m*3+1+1];
		kk=bzv1[whichzone][m*3+2+1];

		// this checks if want to force reflective condition
		if((ii-i==0)&&(jj-j==0)&&(kk-k==0)){ innerreflectflag=1; }

		ftemp+=workv[ll][1][kk][jj][ii];
	      }// end over zones
	      // now assign average value
	      ftemp=ftemp/((FTYPE)(bzv1[whichzone][0]));
	      if(innerreflectflag){
		ftemp=0;
	      }
	      workv[ll][1][k][j][i]=ftemp;

	      break;
	    case 3:
	      /* Fix/Time vary: Dirichlet */
	      if( (wvec>0)||(wvec==-1)){
		workv[ll][1][k][j][i]=vanal[ll][1][k][j][i];
	      }
	      else{
		fprintf(fail_file,"Case 3 bound with wvec==-2 has no definition\n");
		myexit(1);
	      }
	      break;
	    case 4:
	      if(( (wvec==-3)||(wvec==-4))) break;
	      /* outflow */
	      // deal with asymmetry in velocity components on grid w.r.t. inner/outer edges
	      // deals also with inflow checking
	      // order important for if/elseif's
	      
	      //flux can't be outflowed like velocity(i.e. leads to constant density on boundaries).  Thus, only periodic and MPI transfer should assign mass flux on boundaries.  Thus, otherwise should not bound mass flux.  Let computation of mass flux(which is ok out to -1 & N) get new mass, etc.
	      // only overwrite mdot[-1,-2...] and mdot[N+1,N+2...] like terms for consistency
	      
	      // find average value
	      ftemp=0;
	      ininflag=0;
	      outinflag=0;
	      itemp=bzv1[whichzone][0];
	      for(m=0;m<itemp;m++){
		
		ii=bzv1[whichzone][m*3+0+1];
		jj=bzv1[whichzone][m*3+1+1];
		kk=bzv1[whichzone][m*3+2+1];

		// assumes that 1+ boundary zone(along that line) is enough to deal with inflow problem (otherwise no generally valid way to deal with other zones)
		if((wbound==1)&&INFLOWCHECKIX1&&(ii-i>0)&&(jj-j==0)&&(kk-k==0)){ ininflag=1; } // so must force inflow condition
		if((wbound==1)&&INFLOWCHECKOX1&&(ii-i<0)&&(jj-j==0)&&(kk-k==0)){ outinflag=1; } // so must force inflow condition
		      
		ftemp+=workv[ll][1][kk][jj][ii];
	      } // end over zones
	      // inflow check now
	      ftemp=ftemp/((FTYPE)(itemp));// holds average value now
	      if((ininflag)&&(ftemp>0)){
		workv[ll][1][k][j][i]=0.0;
	      }
	      else if((outinflag)&&(ftemp<0)){
		workv[ll][1][k][j][i]=0.0;
	      }
	      else workv[ll][1][k][j][i]=ftemp;
	      
	      break;
	    case 5:
	      // note that loop is like scalar for this then
	      // periodic	    
	      // find average value
	      ftemp=0;
	      for(m=0;m<bzv1[whichzone][0];m++){
		
		ii=bzv1[whichzone][m*3+0+1];
		jj=bzv1[whichzone][m*3+1+1];
		kk=bzv1[whichzone][m*3+2+1];
		
		ftemp+=workv[ll][1][kk][jj][ii];
	      }// end over zones
	      ftemp=ftemp/((FTYPE)(bzv1[whichzone][0]));
	      workv[ll][1][k][j][i]=ftemp;
	      
	      break;
	    case 99:
	    case 98:
	      // no nothing
	      break;
	    default:
	      fprintf(fail_file,"bound.c: switch(bct) error, case given is: %d\n\r",bct);
	      myexit(1);
	    } // end switch
	    //	  } // end if bct>0 for vectors
	}// end main loop
      }// end if vx/bx


      // now do vy: just copy of above with appropriate changes
      if(docom[2]){
	LOOPBOUNDV2{
	  whichzone=temptempi;
	  bct=maskv2[k][j][i];
	  
	  //	if( (bct>0)&&(bct<10) ){ // skip zone if computational zone or null zone (can get rid of this when satisfied that bound loop is correct)// GODMARK
	    
	    switch(bct){
	      
	    case 1:
	    case 2:
	      /* reflect and AOS */
	      // order important for if/elseif's
	      
	      // find average value
	      ftemp=0;
	      innerreflectflag=0;
	      for(m=0;m<bzv2[whichzone][0];m++){
		
		ii=bzv2[whichzone][m*3+0+1];
		jj=bzv2[whichzone][m*3+1+1];
		kk=bzv2[whichzone][m*3+2+1];

		// this checks if want to force reflective condition
		if((ii-i==0)&&(jj-j==0)&&(kk-k==0)){ innerreflectflag=1; }

		ftemp+=workv[ll][2][kk][jj][ii];
	      }// end over zones
	      // now assign average value
	      ftemp=ftemp/((FTYPE)(bzv2[whichzone][0]));
	      if(innerreflectflag){
		ftemp=0;
	      }
	      workv[ll][2][k][j][i]=ftemp;

	      break;
	    case 3:
	      /* Fix/Time vary: Dirichlet */
	      if( (wvec>0)||(wvec==-1)){
		workv[ll][2][k][j][i]=vanal[ll][2][k][j][i];
	      }
	      else{
		fprintf(fail_file,"Case 3 bound with wvec==-2 has no definition\n");
		myexit(1);
	      }
	      break;
	    case 4:
	      if(( (wvec==-3)||(wvec==-4))) break;
	      /* outflow */
	      // deal with asymmetry in velocity components on grid w.r.t. inner/outer edges
	      // deals also with inflow checking
	      // order important for if/elseif's
	      
	      //flux can't be outflowed like velocity(i.e. leads to constant density on boundaries).  Thus, only periodic and MPI transfer should assign mass flux on boundaries.  Thus, otherwise should not bound mass flux.  Let computation of mass flux(which is ok out to -1 & N) get new mass, etc.
	      // only overwrite mdot[-1,-2...] and mdot[N+1,N+2...] like terms for consistency
	      
	      // find average value
	      ftemp=0;
	      ininflag=0;
	      outinflag=0;
	      for(m=0;m<bzv2[whichzone][0];m++){
		
		ii=bzv2[whichzone][m*3+0+1];
		jj=bzv2[whichzone][m*3+1+1];
		kk=bzv2[whichzone][m*3+2+1];

		// assumes that 1+ boundary zone(along that line) is enough to deal with inflow problem (otherwise no generally valid way to deal with other zones)
		if((wbound==1)&&INFLOWCHECKIX2&&(ii-i==0)&&(jj-j>0)&&(kk-k==0)){ ininflag=1; } // so must force inflow condition
		if((wbound==1)&&INFLOWCHECKOX2&&(ii-i==0)&&(jj-j<0)&&(kk-k==0)){ outinflag=1; } // so must force inflow condition
		      
		ftemp+=workv[ll][2][kk][jj][ii];
	      } // end over zones
	      // inflow check now
	      ftemp=ftemp/((FTYPE)(bzv2[whichzone][0]));// holds average value now
	      if((ininflag)&&(ftemp>0)){
		workv[ll][2][k][j][i]=0.0;
	      }
	      else if((outinflag)&&(ftemp<0)){
		workv[ll][2][k][j][i]=0.0;
	      }
	      else workv[ll][2][k][j][i]=ftemp;
	      
	      break;
	    case 5:
	      // note that loop is like scalar for this then
	      // periodic	    
	      // find average value
	      ftemp=0;
	      for(m=0;m<bzv2[whichzone][0];m++){
		
		ii=bzv2[whichzone][m*3+0+1];
		jj=bzv2[whichzone][m*3+1+1];
		kk=bzv2[whichzone][m*3+2+1];
		
		ftemp+=workv[ll][2][kk][jj][ii];
	      }// end over zones
	      ftemp=ftemp/((FTYPE)(bzv2[whichzone][0]));
	      workv[ll][2][k][j][i]=ftemp;
	      
	      break;
	    case 99:
	    case 98:
	      // no nothing
	      break;
	    default:
	      fprintf(fail_file,"bound.c: switch(bct) error, case given is: %d\n\r",bct);
	      myexit(1);
	    } // end switch
	    //	  } // end if bct>0 for vectors
	}// end main loop
      }// end if vy/by
      
      if(docom[3]){
	LOOPBOUNDV3{
	  whichzone=temptempi;
	  bct=maskv3[k][j][i];
	  
	  //	if( (bct>0)&&(bct<10) ){ // skip zone if computational zone or null zone (can get rid of this when satisfied that bound loop is correct)// GODMARK
	    
	    switch(bct){
	      
	    case 1:
	    case 2:
	      /* reflect and AOS */
	      // order important for if/elseif's
	      
	      // find average value
	      ftemp=0;
	      innerreflectflag=0;
	      for(m=0;m<bzv3[whichzone][0];m++){
		
		ii=bzv3[whichzone][m*3+0+1];
		jj=bzv3[whichzone][m*3+1+1];
		kk=bzv3[whichzone][m*3+2+1];

		// this checks if want to force reflective condition
		if((ii-i==0)&&(jj-j==0)&&(kk-k==0)){ innerreflectflag=1; }

		ftemp+=workv[ll][3][kk][jj][ii];
	      }// end over zones
	      // now assign average value
	      ftemp=ftemp/((FTYPE)(bzv3[whichzone][0]));
	      if(innerreflectflag){
		ftemp=0;
	      }
	      workv[ll][3][k][j][i]=ftemp;

	      break;
	    case 3:
	      /* Fix/Time vary: Dirichlet */
	      if( (wvec>0)||(wvec==-1)){
		workv[ll][3][k][j][i]=vanal[ll][3][k][j][i];
	      }
	      else{
		fprintf(fail_file,"Case 3 bound with wvec==-2 has no definition\n");
		myexit(1);
	      }
	      break;
	    case 4:
	      if(( (wvec==-3)||(wvec==-4))) break;
	      /* outflow */
	      // deal with asymmetry in velocity components on grid w.r.t. inner/outer edges
	      // deals also with inflow checking
	      // order important for if/elseif's
	      
	      //flux can't be outflowed like velocity(i.e. leads to constant density on boundaries).  Thus, only periodic and MPI transfer should assign mass flux on boundaries.  Thus, otherwise should not bound mass flux.  Let computation of mass flux(which is ok out to -1 & N) get new mass, etc.
	      // only overwrite mdot[-1,-2...] and mdot[N+1,N+2...] like terms for consistency
	      
	      // find average value
	      ftemp=0;
	      ininflag=0;
	      outinflag=0;
	      for(m=0;m<bzv3[whichzone][0];m++){
		
		ii=bzv3[whichzone][m*3+0+1];
		jj=bzv3[whichzone][m*3+1+1];
		kk=bzv3[whichzone][m*3+2+1];

		// assumes that 1+ boundary zone(along that line) is enough to deal with inflow problem (otherwise no generally valid way to deal with other zones)
		if((wbound==1)&&INFLOWCHECKIX3&&(ii-i==0)&&(jj-j==0)&&(kk-k>0)){ ininflag=1; } // so must force inflow condition
		if((wbound==1)&&INFLOWCHECKOX3&&(ii-i==0)&&(jj-j==0)&&(kk-k<0)){ outinflag=1; } // so must force inflow condition
		      
		ftemp+=workv[ll][3][kk][jj][ii];
	      } // end over zones
	      // inflow check now
	      ftemp=ftemp/((FTYPE)(bzv3[whichzone][0]));// holds average value now
	      if((ininflag)&&(ftemp>0)){
		workv[ll][3][k][j][i]=0.0;
	      }
	      else if((outinflag)&&(ftemp<0)){
		workv[ll][3][k][j][i]=0.0;
	      }
	      else workv[ll][3][k][j][i]=ftemp;
	      
	      break;
	    case 5:
	      // note that loop is like scalar for this then
	      // periodic	    
	      // find average value
	      ftemp=0;
	      for(m=0;m<bzv3[whichzone][0];m++){
		
		ii=bzv3[whichzone][m*3+0+1];
		jj=bzv3[whichzone][m*3+1+1];
		kk=bzv3[whichzone][m*3+2+1];
		
		ftemp+=workv[ll][3][kk][jj][ii];
	      }// end over zones
	      ftemp=ftemp/((FTYPE)(bzv3[whichzone][0]));
	      workv[ll][3][k][j][i]=ftemp;
	      
	      break;
	    case 99:
	    case 98:
	      // no nothing
	      break;
	    default:
	      fprintf(fail_file,"bound.c: switch(bct) error, case given is: %d\n\r",bct);
	      myexit(1);
	    } // end switch
	    // } // end if bct>0 for vectors
	}// end main loop
      }// end if vz/bz

    }// end LOOP over current vector
  

    // now fix divb if necessary
    // should consider how to combine with above
    if(mag==1){
      doingvectors=1;
      doingscalars=0;
      
      
      wbound=l=ll=2;

      if(docom[1]){
	LOOPBOUNDV1{ // now fix for divb
	  whichzone=temptempi;
	  bct=maskv1[k][j][i];
	
	  //	  if( (bct>0)&&(bct<10) ){ // skip zone if computational zone or null zone (can get rid of this when satisfied that bound loop is correct)// GODMARK
	  
	    // for this, just find whether difference between ii and i is generally + or -, sufficient for what's wanted
	    iisum=jjsum=kksum=0;
	    for(m=0;m<bzv1[whichzone][0];m++){
	      
	      ii=bzv1[whichzone][m*3+0+1];
	      jj=bzv1[whichzone][m*3+1+1];
	      kk=bzv1[whichzone][m*3+2+1];
	      
	      iisum+=(ii-i);
	      jjsum+=(jj-j);
	      kksum+=(kk-k);
	    }
	    if((j<=N2)&&(k<=N3)){ // only do if valid memory element exists for b1()
	      if(iisum>0){
		workv[ll][1][k][j][i]=b1(-1,workv[ll],k,j,i);
	      }
	      else workv[ll][1][k][j][i]=b1(1,workv[ll],k,j,i); // true for iisum<=0, 0 case could mean sum'ed to 0 or is 0.  If is 0, then really don't need to do this, but if summed to 0, then should, but doesn't hurt to do it in general (maybe--GODMARK)
	      // otherwise taken care of by other components
	    }
	    //  }
	}
      }
      if(docom[2]){
	LOOPBOUNDV2{ // now fix for divb
	  whichzone=temptempi;
	  bct=maskv2[k][j][i];
	
	  //  if( (bct>0)&&(bct<10) ){ // skip zone if computational zone or null zone (can get rid of this when satisfied that bound loop is correct)// GODMARK

	  
	    // for this, just find whether difference between ii and i is generally + or -, sufficient for what's wanted
	    iisum=jjsum=kksum=0;
	    for(m=0;m<bzv2[whichzone][0];m++){
	      
	      ii=bzv2[whichzone][m*3+0+1];
	      jj=bzv2[whichzone][m*3+1+1];
	      kk=bzv2[whichzone][m*3+2+1];
	      
	      iisum+=(ii-i);
	      jjsum+=(jj-j);
	      kksum+=(kk-k);
	    }
	    if((i<=N1)&&(k<=N3)){ // only do if valid memory element exists for b2()
	      if(jjsum>0){
		workv[ll][2][k][j][i]=b2(-1,workv[ll],k,j,i);
	      }
	      else workv[ll][2][k][j][i]=b2(1,workv[ll],k,j,i);
	      // otherwise taken care of by other components
	    }
	    // }
	}
      }
      if(docom[3]){
	LOOPBOUNDV3{ // now fix for divb
	  whichzone=temptempi;
	  bct=maskv3[k][j][i];
	
	  // if( (bct>0)&&(bct<10) ){ // skip zone if computational zone or null zone (can get rid of this when satisfied that bound loop is correct)// GODMARK

	  
	    // for this, just find whether difference between ii and i is generally + or -, sufficient for what's wanted
	    iisum=jjsum=kksum=0;
	    for(m=0;m<bzv3[whichzone][0];m++){
	      
	      ii=bzv3[whichzone][m*3+0+1];
	      jj=bzv3[whichzone][m*3+1+1];
	      kk=bzv3[whichzone][m*3+2+1];
	      
	      iisum+=(ii-i);
	      jjsum+=(jj-j);
	      kksum+=(kk-k);
	    }
	    if((i<=N1)&&(j<=N2)){ // only do if valid memory element exists for b3()
	      if(kksum>0){
		workv[ll][3][k][j][i]=b3(-1,workv[ll],k,j,i);
	      }
	      else workv[ll][3][k][j][i]=b3(1,workv[ll],k,j,i);
	      // otherwise taken care of by other components
	    }
	  //}
	}
      }
    }// end divb fix
  }// end if doing vectors
#endif // end if do bound vectors



  // if doing MPI:
#if(USEMPI)
  bound_mpi(vars,varv,wsca,wvec,wcom);
#endif

  firsttime=0;
}// end function
