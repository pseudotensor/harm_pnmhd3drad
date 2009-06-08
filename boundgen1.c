#include "bound.h"
#include "boundgen.h"

extern void bound_ppclean(FTYPE (*vars)[N2M][N1M],
		  FTYPE (*varv)[N3M][N2M][N1M],
		  int wsca,
		  int wvec,
		  int wcom)
{
}

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

#if(DOTRUEBOUNDARY==1)

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
	else if(wsca==-11){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-12){//other
	  // bound like en
	  wbound=2;
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
	  
	  bct=bzs[looper][temptempi][-3];
	  
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
	    if(wsca<-1) break;
	    works[ll][k][j][i]=sanal[ll][k][j][i];
	    break;
	  case 90:
	  case 91:
	  case 92:
	  case 93:
	  case 94:
	  case 95:
	  case 99:
	  case 98:
	    // do nothing
	    break;
	  default:
	    fprintf(fail_file,"bound.c: switch(bct) error, case given is: %d\n\r",bct);
	    myexit(1);
	  } // end switch
	}//end over current scalar boundary zones(minus corners)
	/* cut short loop if only to do one */
      } // over various scalars(in/out)
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
	else if((wvec==-5)){//mdot or flux(2) of v[1][2] or magemf
	  // bound like v
	  wbound=3;
	}
	else if(wvec==-11){// v1
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-12){// v2
	  // bound like B
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0, wvec=%d\n",wvec);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(wbound==1){ looperstart=8; looperend=13;}
      if(wbound==2){ looperstart=14; looperend=19;}
      if(wbound==3){ looperstart=8; looperend=13;}
      


      for(looper=looperstart;looper<=looperend;looper++){


	// in? out?
	if( (!doinner)&&((looper==8)||(looper==10)||(looper==12)||(looper==14)||(looper==16)||(looper==18)) ) continue;
	if( (!doouter)&&((looper==9)||(looper==11)||(looper==13)||(looper==15)||(looper==17)||(looper==19)) ) continue;

	// which component?
	if((looper==8)||(looper==9)||(looper==14)||(looper==15)) component=1;
	if((looper==10)||(looper==11)||(looper==16)||(looper==17)) component=2;
	if((looper==12)||(looper==13)||(looper==18)||(looper==19)) component=3;

	if(!docom[component]) continue;

	
	
	LOOPBOUND(looper){
	  
	  bct=bzs[looper][temptempi][-3];
	  
	  switch(bct){
	    
	  case 1:
	  case 2:
	    /* reflect and AOS */
	    // as outflow with no inflow check
	    ftemp=0;
	    itemp=bzs[looper][temptempi][0];
	    for(m=0;m<itemp;m++) ftemp+=workv[ll][component][bzs[looper][temptempi][m*3+2+1]][bzs[looper][temptempi][m*3+1+1]][bzs[looper][temptempi][m*3+0+1]];
	    // inflow check now
	    ftemp=ftemp/((FTYPE)(itemp));// holds average value now
	    // may require reflect flag addition
	    workv[ll][1][k][j][i]=ftemp;
	    break;
	  case 3:
	    if(wvec<-1) break; // since only real primitives can be assigned currently
	    /* Fix/Time vary: Dirichlet */
	    workv[ll][component][k][j][i]=vanal[ll][component][k][j][i];
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
	    itemp=bzs[looper][temptempi][0];
	    for(m=0;m<itemp;m++) ftemp+=workv[ll][component][bzs[looper][temptempi][m*3+2+1]][bzs[looper][temptempi][m*3+1+1]][bzs[looper][temptempi][m*3+0+1]];
	    // inflow check now
	    ftemp=ftemp/((FTYPE)(itemp));// holds average value now
	    // then do an inflow check
	    if(wbound==1){
	      if((bzs[looper][temptempi][-1])&&(ftemp>0)){
		workv[ll][component][k][j][i]=0.0;
	      }
	      else if((bzs[looper][temptempi][-2])&&(ftemp<0)){
		workv[ll][component][k][j][i]=0.0;
	      }
	      else workv[ll][component][k][j][i]=ftemp;
	    }
	    else workv[ll][component][k][j][i]=ftemp;	      
	    break;
	  case 5:
	    // periodic	    
	    // as outflow with no inflow check
	    ftemp=0;
	    itemp=bzs[looper][temptempi][0];
	    for(m=0;m<itemp;m++) ftemp+=workv[ll][component][bzs[looper][temptempi][m*3+2+1]][bzs[looper][temptempi][m*3+1+1]][bzs[looper][temptempi][m*3+0+1]];
	    ftemp=ftemp/((FTYPE)(itemp));// holds average value now
	    workv[ll][1][k][j][i]=ftemp;	      
	    break;
	  case 90:
	  case 91:
	  case 92:
	  case 93:
	  case 94:
	  case 95:
	  case 99:
	  case 98:
	    // no nothing
	    break;
	  default:
	    fprintf(fail_file,"bound.c: switch(bct) error, case given is: %d\n\r",bct);
	    myexit(1);
	  } // end switch	  
	}// end main loop
      }// end over various in/outs and vectors
    }// end over vectors
  }// end if doing vectors
#endif // end if do bound vectors

#endif// end if doing true boundary

#if(DOINTERNALBOUNDARY==1)
  // if doing MPI:
#if(USEMPI)
  bound_mpi(vars,varv,wsca,wvec,wcom);
#endif
#endif
  firsttime=0;
}// end function
