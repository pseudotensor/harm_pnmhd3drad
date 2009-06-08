#include "bound.h"
#include "boundrect.h"

#define NUMCORNERS 8
#define NUMEDGES 12

// this bound code is solid working for any rectangular grid in any dimension, any coordinate system

void bound_rect1(FTYPE (*vars)[N2M][N1M],
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
  FTYPE slope;
  FTYPE (*works)[N3M][N2M][N1M];
  FTYPE (*workv)[3][N3M][N2M][N1M];
  static int numhit=0;
  int numhitmin,numhitmax;
  int wbound;
  FTYPE ftemp;
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
  static int numcorners,numcornerzones;
  static int si[NUMCORNERS],sj[NUMCORNERS],sk[NUMCORNERS],gi[NUMCORNERS],gj[NUMCORNERS],gk[NUMCORNERS];
  static int ix[NUMCORNERS][NUMCORNERS][4][3+1]; // ix[corner][zone in corner][i/p/o/r(a)][dir=1,2,3] (zone order such that inner corner first, then rest of 2D layer, then 3D layer last 4
  static int numedges,numedgezones;
  static int edgesi[NUMEDGES],edgesj[NUMEDGES],edgesk[NUMEDGES],edgegi[NUMEDGES],edgegj[NUMEDGES],edgegk[NUMEDGES],edgeend[NUMEDGES];
  static int edgeix[NUMEDGES][4][4][3+1]; // edgeix[edge#][zone in edge][i/p/o/r(a)][dir=1,2,3] (zone order such that inner corner first, then rest of 2D layer, then 3D layer last 4

  int othercomp;
  int cornerto,corner;
  int edgeto, edge;




  // might not want to bound vectors if not advecting them
  //  if(wvec==-1) return;
  //if(wvec==1) return;

  numhit++;
  if(wsca<=-2){
    works=(FTYPE (*) [N3M][N2M][N1M])(&vars[0][0][0]);
  }
  else works=s;
  if(wvec<=-2){
    workv=(FTYPE (*) [3][N3M][N2M][N1M])(&varv[0][0][0][0]);
  }
  else workv=v;

#if(FORCEOUTFLOW)
  // check to see if user wants to force outflow on some variable
  forces=forcev=0;
  // force outflow for arbitrary incoming variable
  if(wsca==-2) forces++; // force arbitrary to be outflow

  if(forces||forcev) forcecheck=1;
  else forcecheck=0;
#else
  forces=forcev=0;
  forcecheck=0;
#endif


  if(firsttime){
    if(COMPDIM==3){
      // corners always exist, some parts of which may not actually be used, but ok
      numcorners=8;
      numcornerzones=8;
      
      // setup corner zone for simple loops
      si[0]=-N1OFF;
      sj[0]=-N2OFF;
      sk[0]=-N3OFF;
      gi[0]=N1OFF;
      gj[0]=N2OFF;
      gk[0]=N3OFF;
      
      si[1]=N1*N1OFF;
      sj[1]=-N2OFF;
      sk[1]=-N3OFF;
      gi[1]=-N1OFF;
      gj[1]=N2OFF;
      gk[1]=N3OFF;
      
      si[2]=N1*N1OFF;
      sj[2]=N2*N2OFF;
      sk[2]=-N3OFF;
      gi[2]=-N1OFF;
      gj[2]=-N2OFF;
      gk[2]=N3OFF;
      
      si[3]=-N1OFF;
      sj[3]=N2*N2OFF;
      sk[3]=-N3OFF;
      gi[3]=N1OFF;
      gj[3]=-N2OFF;
      gk[3]=N3OFF;
      
      si[4]=-N1OFF;
      sj[4]=-N2OFF;
      sk[4]=N3*N3OFF;
      gi[4]=N1OFF;
      gj[4]=N2OFF;
      gk[4]=-N3OFF;
      
      si[5]=N1*N1OFF;
      sj[5]=-N2OFF;
      sk[5]=N3*N3OFF;
      gi[5]=-N1OFF;
      gj[5]=N2OFF;
      gk[5]=-N3OFF;
      
      si[6]=N1*N1OFF;
      sj[6]=N2*N2OFF;
      sk[6]=N3*N3OFF;
      gi[6]=-N1OFF;
      gj[6]=-N2OFF;
      gk[6]=-N3OFF;
      
      si[7]=-N1OFF;
      sj[7]=N2*N2OFF;
      sk[7]=N3*N3OFF;
      gi[7]=N1OFF;
      gj[7]=-N2OFF;
      gk[7]=-N3OFF;

      if((N1>1)&&(N2>1)&&(N3>1)){ // edges only exist in full 3D, don't do extra work (i.e. corners become edges via above assignments)
	numedges=12;
	numedgezones=4;
      }
      else{
	numedges=0;
	numedgezones=0;
      }
      // always start from "bottom" and loop to higher index

      // no need to do N?OFF's here since only used if N?OFF's are all ==1

      // C0 to C4 : E4
      edgesi[4]=-1;
      edgesj[4]=-1;
      edgesk[4]=0; // start of loop
      edgeend[4]=N3-1; // end of loop
      edgegi[4]=1;
      edgegj[4]=1;
      edgegk[4]=0; // what's looped over

      // C0 to C1 : E0
      edgesi[0]=0;
      edgeend[0]=N1-1;
      edgesj[0]=-1;
      edgesk[0]=-1;
      edgegi[0]=0;
      edgegj[0]=1;
      edgegk[0]=1; 

      // C1 to C5 : E5
      edgesi[5]=N1;
      edgesj[5]=-1;
      edgesk[5]=0;
      edgeend[5]=N3-1;
      edgegi[5]=-1;
      edgegj[5]=1;
      edgegk[5]=0; 

      // C4 to C5 : E8
      edgesi[8]=0;
      edgeend[8]=N1-1;
      edgesj[8]=-1;
      edgesk[8]=N3;
      edgegi[8]=0;
      edgegj[8]=1;
      edgegk[8]=-1;

      // C0 to C3 : E3
      edgesi[3]=-1;
      edgesj[3]=0;
      edgeend[3]=N2-1;
      edgesk[3]=-1;
      edgegi[3]=1;
      edgegj[3]=0;
      edgegk[3]=1; 

      // C1 to C2: E1
      edgesi[1]=N1;
      edgesj[1]=0;
      edgeend[1]=N2-1;
      edgesk[1]=-1;
      edgegi[1]=-1;
      edgegj[1]=0;
      edgegk[1]=1; 

      // C5 to C6 : E9
      edgesi[9]=N1;
      edgesj[9]=0;
      edgeend[9]=N2-1;
      edgesk[9]=N3;
      edgegi[9]=-1;
      edgegj[9]=0;
      edgegk[9]=-1; 

      // C4 to C7 : E11
      edgesi[11]=-1;
      edgesj[11]=0;
      edgeend[11]=N2-1;
      edgesk[11]=N3;
      edgegi[11]=1;
      edgegj[11]=0;
      edgegk[11]=-1; 

      // C3 to C7 : E7
      edgesi[7]=-1;
      edgesj[7]=N2;
      edgesk[7]=0;
      edgeend[7]=N3-1;
      edgegi[7]=1;
      edgegj[7]=-1;
      edgegk[7]=0;

      // C3 to C2 : E2
      edgesi[2]=0;
      edgeend[2]=N1-1;
      edgesj[2]=N2;
      edgesk[2]=-1;
      edgegi[2]=0;
      edgegj[2]=-1;
      edgegk[2]=1;

      // C2 to C6 : E6
      edgesi[6]=N1;
      edgesj[6]=N2;
      edgesk[6]=0;
      edgeend[6]=N3-1;
      edgegi[6]=-1;
      edgegj[6]=-1;
      edgegk[6]=0;

      // C7 to C6 : E10
      edgesi[10]=0;
      edgeend[10]=N1-1;
      edgesj[10]=N2;
      edgesk[10]=N3;
      edgegi[10]=0;
      edgegj[10]=-1;
      edgegk[10]=-1;


    }
    else{
      numcorners=4;
      numcornerzones=4;
      numedges=0;
      numedgezones=0;

      si[0]=-N1OFF;
      sj[0]=-N2OFF;
      sk[0]=0;
      gi[0]=N1OFF;
      gj[0]=N2OFF;
      gk[0]=0;
      
      si[1]=N1*N1OFF;
      sj[1]=-N2OFF;
      sk[1]=0;
      gi[1]=-N1OFF;
      gj[1]=N2OFF;
      gk[1]=0;
      
      si[2]=N1*N1OFF;
      sj[2]=N2*N2OFF;
      sk[2]=0;
      gi[2]=-N1OFF;
      gj[2]=-N2OFF;
      gk[2]=0;
      
      si[3]=-N1OFF;
      sj[3]=N2*N2OFF;
      sk[3]=0;
      gi[3]=N1OFF;
      gj[3]=-N2OFF;
      gk[3]=0;
      
    }

    // now setup zones so know where data comes from for given call type

    // ix[corner][zone#][kk,kko,kkk][i,j,k]
    // just need to make sure that far zone in corner is last, and inner corner zone first(for ordering of divb corrections below to be valid without extra loops)
    for(corner=0;corner<numcorners;corner++){
      i=si[corner];
      j=sj[corner];
      k=sk[corner];
      // inner corner zones
      // for assignment and itself
      ii=i-gi[corner];
      jj=j-gj[corner];
      kk=k-gk[corner];
      
      // for reflective or outflow
      iio=i+gi[corner];
      jjo=j+gj[corner];
      kko=k+gk[corner];

      // below for reflective condition(2nd layer)
      if(N1==1) jumpfactor=1; else jumpfactor=2;
      iio2=i+jumpfactor*gi[corner];
      if(N2==1) jumpfactor=1; else jumpfactor=2;
      jjo2=j+jumpfactor*gj[corner];
      if(N3==1) jumpfactor=1; else jumpfactor=2;
      kko2=k+jumpfactor*gk[corner];
      
      // periodic exchange zones
      iii=i+gi[corner]*N1;
      jjj=j+gj[corner]*N2;
      kkk=k+gk[corner]*N3;

      // 2nd layer of periodic
      iii2=i+gi[corner]*(N1-1);
      jjj2=j+gj[corner]*(N2-1);
      kkk2=k+gk[corner]*(N3-1);
      
      for(p=0;p<numcornerzones;p++){ // zone#
	for(q=1;q<=3;q++){ // component
	  if(q==3){
	    if(p/4==0){ ix[corner][p][0][q]=k; ix[corner][p][1][q]=kkk;  ix[corner][p][2][q]=kko; ix[corner][p][3][q]=kko;  }
	    else{  ix[corner][p][0][q]=kk; ix[corner][p][1][q]=kkk2; ix[corner][p][2][q]=kko; ix[corner][p][3][q]=kko2;  }
	  }
	  if(q==2){
	    if(!((p/2)%2)){ ix[corner][p][0][q]=j;  ix[corner][p][1][q]=jjj; ix[corner][p][2][q]=jjo;  ix[corner][p][3][q]=jjo; }
	    else{  ix[corner][p][0][q]=jj; ix[corner][p][1][q]=jjj2;  ix[corner][p][2][q]=jjo; ix[corner][p][3][q]=jjo2;  }
	  }
	  if(q==1){
	    if(!(p%2)){ ix[corner][p][0][q]=i;  ix[corner][p][1][q]=iii; ix[corner][p][2][q]=iio; ix[corner][p][3][q]=iio;  }
	    else{  ix[corner][p][0][q]=ii;  ix[corner][p][1][q]=iii2;  ix[corner][p][2][q]=iio; ix[corner][p][3][q]=iio2;  }
	  }
	}
      }
    }


    // now setup edges for 3D calc
    // edgeix[edge][zone#][kk,kko,kkk][i,j,k]
    for(edge=0;edge<numedges;edge++){
      i=edgesi[edge];
      j=edgesj[edge];
      k=edgesk[edge];
      // inner corner zones
      // for assignment and itself
      ii=i-edgegi[edge];
      jj=j-edgegj[edge];
      kk=k-edgegk[edge];
      
      // for reflective or outflow
      iio=i+edgegi[edge];
      jjo=j+edgegj[edge];
      kko=k+edgegk[edge];

      // below for reflective condition(2nd layer)
      if(N1==1) jumpfactor=1; else jumpfactor=2;
      iio2=i+jumpfactor*edgegi[edge];
      if(N2==1) jumpfactor=1; else jumpfactor=2;
      jjo2=j+jumpfactor*edgegj[edge];
      if(N3==1) jumpfactor=1; else jumpfactor=2;
      kko2=k+jumpfactor*edgegk[edge];
      
      // periodic exchange zones
      iii=i+edgegi[edge]*N1;
      jjj=j+edgegj[edge]*N2;
      kkk=k+edgegk[edge]*N3;

      // 2nd layer of periodic
      iii2=i+edgegi[edge]*(N1-1);
      jjj2=j+edgegj[edge]*(N2-1);
      kkk2=k+edgegk[edge]*(N3-1);

      // make sure inner zone first, outer zone last for all slices
      
      if(edgegk[edge]==0){// then permute over i,j and just assign either k

	for(p=0;p<numedgezones;p++){ // zone#

	  q=3;
	  edgeix[edge][p][0][q]=k; edgeix[edge][p][1][q]=kkk;  edgeix[edge][p][2][q]=kko; edgeix[edge][p][3][q]=kko;

	  q=2;
	  if(!(p%2)){ edgeix[edge][p][0][q]=j;  edgeix[edge][p][1][q]=jjj; edgeix[edge][p][2][q]=jjo;  edgeix[edge][p][3][q]=jjo; }
	  else{  edgeix[edge][p][0][q]=jj; edgeix[edge][p][1][q]=jjj2;  edgeix[edge][p][2][q]=jjo; edgeix[edge][p][3][q]=jjo2;  }
	  
	  q=1;
	  if(!(p/2)){ edgeix[edge][p][0][q]=i;  edgeix[edge][p][1][q]=iii; edgeix[edge][p][2][q]=iio; edgeix[edge][p][3][q]=iio;  }
	  else{  edgeix[edge][p][0][q]=ii;  edgeix[edge][p][1][q]=iii2;  edgeix[edge][p][2][q]=iio; edgeix[edge][p][3][q]=iio2;  }
	}
      }
      if(edgegj[edge]==0){// then permute over i,k and just assign either j

	for(p=0;p<numedgezones;p++){ // zone#

	  q=2;
	  edgeix[edge][p][0][q]=j;  edgeix[edge][p][1][q]=jjj; edgeix[edge][p][2][q]=jjo;  edgeix[edge][p][3][q]=jjo;

	  q=1;
	  if(!(p%2)){ edgeix[edge][p][0][q]=i;  edgeix[edge][p][1][q]=iii; edgeix[edge][p][2][q]=iio; edgeix[edge][p][3][q]=iio;  }
	  else{  edgeix[edge][p][0][q]=ii;  edgeix[edge][p][1][q]=iii2;  edgeix[edge][p][2][q]=iio; edgeix[edge][p][3][q]=iio2;  }
	  
	  q=3;
	  if(!(p/2)){ edgeix[edge][p][0][q]=k; edgeix[edge][p][1][q]=kkk;  edgeix[edge][p][2][q]=kko; edgeix[edge][p][3][q]=kko; }
	  else{  edgeix[edge][p][0][q]=kk; edgeix[edge][p][1][q]=kkk2; edgeix[edge][p][2][q]=kko; edgeix[edge][p][3][q]=kko2;  }
	}
      }
      if(edgegi[edge]==0){// then permute over j,k and just assign either i

	for(p=0;p<numedgezones;p++){ // zone#

	  q=1;
	  edgeix[edge][p][0][q]=i;  edgeix[edge][p][1][q]=iii; edgeix[edge][p][2][q]=iio; edgeix[edge][p][3][q]=iio;

	  q=2;
	  if(!(p%2)){ edgeix[edge][p][0][q]=j;  edgeix[edge][p][1][q]=jjj; edgeix[edge][p][2][q]=jjo;  edgeix[edge][p][3][q]=jjo; }
	  else{  edgeix[edge][p][0][q]=jj; edgeix[edge][p][1][q]=jjj2;  edgeix[edge][p][2][q]=jjo; edgeix[edge][p][3][q]=jjo2;  }
	  
	  q=3;
	  if(!(p/2)){ edgeix[edge][p][0][q]=k; edgeix[edge][p][1][q]=kkk;  edgeix[edge][p][2][q]=kko; edgeix[edge][p][3][q]=kko; }
	  else{  edgeix[edge][p][0][q]=kk; edgeix[edge][p][1][q]=kkk2; edgeix[edge][p][2][q]=kko; edgeix[edge][p][3][q]=kko2;  }
	}
      }
    }

  }// end if firsttime
















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
      
      
      
      LOOPBOUNDRECT{ // if using 2 boundary zones, still only loop over first layer, assuming 2nd layer is calculatable
#if(!INNERBC)
	SKIPCONDITION; // skip known non-boundary zones-3d version
#endif
       
#if(PUREBC==0)
	bct=bcs[wbound][1][k][j][i];
	bcdim=bcs[wbound][2][k][j][i];
	bcdir=bcs[wbound][3][k][j][i];

#else
	bct=purebccall(k,j,i,&bcdim,&bcdir);	
#endif	
	if(forcecheck){
	  if(bct==3) bct=4; // force dq or mdot to be outflowed(or other) if no assignments put in
	}
	if((bct>0)&&(bct<10)){ /* skip zone if computational zone or null zone(corner zones) */

	  
	  if(bcdim==1){ bcdir1=bcdir; bcdir2=0; bcdir3=0;}
	  else if(bcdim==2){ bcdir1=0; bcdir2=bcdir; bcdir3=0;}
	  else if(bcdim==3){ bcdir1=0; bcdir2=0; bcdir3=bcdir;}
	  
#include "boundpos.h"
	  
	  switch(bct){
	    
	  case 1:
	  case 2:
	    /* reflect/AOS same for scalar */
	    works[ll][k][j][i]=works[ll][kk][jj][ii];
	    // copy 2nd layer of scalar
	    works[ll][k2][j2][i2]=works[ll][kkk][jjj][iii];
	    break;
	  case 4:
	    /* outflow for scalar */
	    works[ll][k][j][i]=works[ll][kk][jj][ii];
	    // copy 2nd layer of scalar
	    works[ll][k2][j2][i2]=works[ll][kk][jj][ii]; // same copy as k,j,i
	    break;
	  case 3:
	    if( (wsca>0)||(wsca==-1)){
	      works[ll][k][j][i]=sanal[ll][k][j][i];
	      works[ll][k2][j2][i2]=sanal[ll][k2][j2][i2];
	    }
	    else if(wsca==-2){
	      fprintf(fail_file,"Case 3 bound with wsca==-2 has no definition\n");
	      fprintf(fail_file,"%d %d %d %d %d\n",forcecheck,forces,forcev,numhit,bct);
	      myexit(1);
	    }
	    break;
	  case 5: // if periodicx2special, then same form as below
	    /* periodic */
	    works[ll][k][j][i]=works[ll][kk][jj][ii];
	    works[ll][k2][j2][i2]=works[ll][kkk][jjj][iii];
	    break;
	  case 90:
	  case 91:
	  case 92:
	  case 93:
	  case 94:
	  case 95:
	  case 96:
	  case 97:
	  case 98:
	  case 99:
	    // do nothing
	    break;
	  default:
	    fprintf(fail_file,"bound.c: switch(bct) error, case given is: %d\n\r",bct);
	    myexit(1);
	  } // end switch
	} // end if bct>0
      }//end over current scalar boundary zones(minus corners)
      
#if(SCABCCORNER)
      // now do corner zones(assumes other boundary zones already done!)
      
      
      // only 3D should hit this (only valid for simplebc edges--not mixed edges--otherwise would need to invert some of the loops which is more costly)
      for(edge=0;edge<numedges;edge++){
	
	i=edgesi[edge];
	j=edgesj[edge];
	k=edgesk[edge];
#if(PUREBC==0)
	// multiply boundary type by 0 if not relevant type(i.e. only 2d ambiguity)
	bireal=bcs[wbound][1][k+edgegk[edge]][j+edgegj[edge]][i]*abs(edgegi[edge]); // i-boundary type
	bjreal=bcs[wbound][1][k+edgegk[edge]][j][i+edgegi[edge]]*abs(edgegj[edge]); // j-boundary type
	bkreal=bcs[wbound][1][k][j+edgegj[edge]][i+edgegi[edge]]*abs(edgegk[edge]); // k-boundary type
#else
	bireal=purebccall(k+edgegk[edge],j+edgegj[edge],i,&bcdim,&bcdir)*abs(edgegi[edge]); // i-boundary type
	bjreal=purebccall(k+edgegk[edge],j,i+edgegi[edge],&bcdim,&bcdir)*abs(edgegj[edge]); // j-boundary type
	bkreal=purebccall(k,j+edgegj[edge],i+edgegi[edge],&bcdim,&bcdir)*abs(edgegk[edge]); // k-boundary type
#endif
	if(bireal==0){
	  bi=bjreal; diri=2;
	  bj=bkreal; dirj=3;
	  startloop=edgesi[edge];
	}
	if(bjreal==0){
	  bi=bireal; diri=1;
	  bj=bkreal; dirj=3;
	  startloop=edgesj[edge];
	}
	if(bkreal==0){
	  bi=bireal; diri=1;
	  bj=bjreal; dirj=2;
	  startloop=edgesk[edge];
	}
	endloop=edgeend[edge];
	
	if(
	   (bi==0)||(bj==0)
	   ){
	  fprintf(fail_file,"s[%d]: k=%d,j=%d,i=%d quad1: undefined boundary bi: %d bj: %d\n",ll,k,j,i,bi,bj);
	  myexit(1);
	}
	else if(
		(bi>=90)||(bj>=90)
		){
	}
	else{ // if good edge
	  
	  if(// i-inflow AND j-inflow or if p-i or i-p 
	     ((bi==3)&&(bj==3))||((bi==5)&&(bj==3))||((bi==3)&&(bj==5))
	     ){
	    
	    typeofbc=1;
	  }
	  // i-periodic AND j-periodic
	  else if(
		  (bi==5)&&(bj==5)
		  ){
	    
	    typeofbc=2;
	  }
	  // i/k-outflow/inflow/periodic AND j-outflow (if o-o, same effect as next case)
	  else if(
		  ((bi==4)||(bi==3)||(bi==5))&&(bj==4) // so "j" is either k or j real dir
		  ){
	    if(dirj==2) typeofbc=3; // so j is j 
	    else typeofbc=4; // so j is k
	  }
	  // j/k-outflow/inflow/periodic AND i-outflow
	  else if(
		  ((bj==4)||(bj==3)||(bj==5))&&(bi==4) // so "i" is either i or j real dir
		  ){
	    if(diri==1) typeofbc=5; // so i is i
	    else typeofbc=6; // so i is j
	  }
	  // i/k-inflow/outflow/periodic/reflect(or AOS) AND j-reflect(or AOS)
	  else if(
		  ((bi==3)||(bi==4)||(bi==5)||(bi==1)||(bi==2))&& ((bj==1)||(bj==2))
		  ){
	    if(dirj==2) typeofbc=7; // so j is j 
	    else typeofbc=8; // so j is k
	  }
	  //  i-reflect(or AOS) AND j/k-inflow/outflow/periodic/reflect(or AOS)
	  else if(
		  ((bj==3)||(bj==4)||(bj==5)||(bj==1)||(bj==2))&& ((bi==1)||(bi==2))
		  ){
	    if(diri==1) typeofbc=9; // so i is i
	    else typeofbc=10; // so i is j
	  }
	  else{
	    fprintf(fail_file,"s[%d]: k=%d j=%d,i=%d quad2: No boundary case setup for bi: %d bj: %d\n",ll,k,j,i,bi,bj);
	    myexit(1);
	  }


	  loopi=!abs(edgegi[edge]);
	  loopj=!abs(edgegj[edge]);
	  loopk=!abs(edgegk[edge]);	  

	  
	  for(p=0;p<numedgezones;p++){
	    // inflow or itself
	    ii=edgeix[edge][p][0][1];
	    jj=edgeix[edge][p][0][2];
	    kk=edgeix[edge][p][0][3];
	    
	    // periodic
	    iii=edgeix[edge][p][1][1];
	    jjj=edgeix[edge][p][1][2];
	    kkk=edgeix[edge][p][1][3];
	    
	    // outflow
	    iio=edgeix[edge][p][2][1];
	    jjo=edgeix[edge][p][2][2];
	    kko=edgeix[edge][p][2][3];
	    
	    // reflective or AOS
	    iir=edgeix[edge][p][3][1];
	    jjr=edgeix[edge][p][3][2];
	    kkr=edgeix[edge][p][3][3];
	    
	    for(loopvar=startloop;loopvar<=endloop;loopvar++){ // perform loop down non-bounded direction

	      //printf("hit: %d %d %d %d %d    %d %d %d\n",edge,p,loopvar,startloop, endloop, kk,jj,ii);

	      // could invert if-loop statement here, but annoying extra code
	      if(typeofbc==1){
		works[ll][kk][jj][ii]=sanal[ll][kk][jj][ii];
	      }
	      else if(typeofbc==2){
		works[ll][kk][jj][ii]=works[ll][kkk][jjj][iii];
	      }
	      
	      else if(typeofbc==3){
		works[ll][kk][jj][ii]=works[ll][kk][jjo][ii];
	      }
	      else if(typeofbc==4){
		works[ll][kk][jj][ii]=works[ll][kko][jj][ii];
	      }
	      else if(typeofbc==5){
		works[ll][kk][jj][ii]=works[ll][kk][jj][iio];
	      }
	      else if(typeofbc==6){
		works[ll][kk][jj][ii]=works[ll][kk][jjo][ii];
	      }
	      
	      else if(typeofbc==7){
		works[ll][kk][jj][ii]=works[ll][kk][jjr][ii];
	      }
	      else if(typeofbc==8){
		works[ll][kk][jj][ii]=works[ll][kkr][jj][ii];
	      }
	      else if(typeofbc==9){
		works[ll][kk][jj][ii]=works[ll][kk][jj][iir];
	      }
	      else if(typeofbc==10){
		works[ll][kk][jj][ii]=works[ll][kk][jjr][ii];
	      }
	      else{
		fprintf(fail_file,"s[%d]: k=%d j=%d,i=%d quad3: No boundary case setup for bi: %d bj: %d\n",ll,k,j,i,bi,bj);
		myexit(1);
	      }
	      // update to next zone over

	      // inflow or itself
	      ii+=loopi;  jj+=loopj;  kk+=loopk;
	      // periodic
	      iii+=loopi; jjj+=loopj; kkk+=loopk;
	      // outflow
	      iio+=loopi; jjo+=loopj; kko+=loopk;
	      // reflective or AOS
	      iir+=loopi; jjr+=loopj; kkr+=loopk;

	    }// over over line
	  }// end over each line
	}// else if good boundary
      }// end over all edges


      // both 2d and 3d hit this
      for(corner=0;corner<numcorners;corner++){
	
	i=si[corner];
	j=sj[corner];
	k=sk[corner];

	// assume if N?OFF==0, then not relevant boundary, -100 signifies irrelevant

	if(N1NOT1){// i-boundary type
#if(PUREBC==0)
	  bi=bcs[wbound][1][k+gk[corner]][j+gj[corner]][i];
#else
	  bi=purebccall(k+gk[corner],j+gj[corner],i,&bcdim,&bcdir);
#endif
	}
	else{
	  bi=5; // periodic like if no direction in this direction
	}

	if(N2NOT1){// j-boundary type
#if(PUREBC==0)
	  bj=bcs[wbound][1][k+gk[corner]][j][i+gi[corner]];
#else
	  bj=purebccall(k+gk[corner],j,i+gi[corner],&bcdim,&bcdir);
#endif
	}
	else{
	  bj=5;
	}

	if((COMPDIM<3)||(N3NOT1==0)){
	  bk=5;
	}
	else{
#if(PUREBC==0)
	  bk=bcs[wbound][1][k][j+gj[corner]][i+gi[corner]];	  
#else
	  bk=purebccall(k,j+gj[corner],i+gi[corner],&bcdim,&bcdir);	  
#endif
	}
	

	if(
#if(COMPDIM<3)
	   (bi==0)||(bj==0)
#else
	   (bi==0)||(bj==0)||(bk==0)
#endif
	   ){
	  fprintf(fail_file,"s[%d]: k=%d,j=%d,i=%d quad4: undefined boundary bi: %d bj: %d bk: %d\n",ll,k,j,i,bi,bj,bk);
	  myexit(1);
	}
	else if(
#if(COMPDIM<3)
		(bi>=90)||(bj>=90)
#else
		(bi>=90)||(bj>=90)||(bk>=90)
#endif
		){
	}
	else{ // if good corner
	  
	  if(// i-inflow AND j-inflow or if p-i or i-p 
#if(COMPDIM<3)
	     ((bi==3)&&(bj==3))||((bi==5)&&(bj==3))||((bi==3)&&(bj==5))
#else
	     ((bi==3)&&((bj==3)||(bj==5))&&((bk==3)||(bk==5)))||
	     ((bj==3)&&((bi==3)||(bi==5))&&((bk==3)||(bk==5)))||
	     ((bk==3)&&((bj==3)||(bj==5))&&((bi==3)||(bi==5)))
#endif
	     ){
	    
	    typeofbc=1;
	  }
	  // i-periodic AND j-periodic
	  else if(
#if(COMPDIM<3)
		  (bi==5)&&(bj==5)
#else
		  (bi==5)&&(bj==5)&&(bk==5)
#endif
		  ){
	    
	    typeofbc=2;
	  }
	  // i/k-outflow/inflow/periodic AND j-outflow (if o-o, same effect as next case)
	  else if(
#if(COMPDIM<3)
		  ((bi==4)||(bi==3)||(bi==5))&&(bj==4)
#else
		  (((bi==4)||(bi==3)||(bi==5))&&(bj==4))&&
		  (((bk==4)||(bk==3)||(bk==5))&&(bj==4))
#endif	  
		  ){
	    typeofbc=3;
	  }
	  // j/k-outflow/inflow/periodic AND i-outflow
	  else if(
#if(COMPDIM<3)
		  ((bj==4)||(bj==3)||(bj==5))&&(bi==4)
#else
		  (((bj==4)||(bj==3)||(bj==5))&&(bi==4))&&
		  (((bk==4)||(bk==3)||(bk==5))&&(bi==4))
#endif
		  ){
	    
	    typeofbc=4;
	  }
	  // j/i-outflow/inflow/periodic AND k-outflow
	  else if(
#if(COMPDIM<3)
		  0
#else
		  (((bj==4)||(bj==3)||(bj==5))&&(bk==4))&&
		  (((bi==4)||(bi==3)||(bi==5))&&(bk==4))
#endif
		  ){
	    
	    typeofbc=5;
	  }
	  // i/k-inflow/outflow/periodic/reflect(or AOS) AND j-reflect(or AOS)
	  else if(
#if(COMPDIM<3)
		  ((bi==3)||(bi==4)||(bi==5)||(bi==1)||(bi==2))&& ((bj==1)||(bj==2))
#else
		  (((bi==3)||(bi==4)||(bi==5)||(bi==1)||(bi==2))&& ((bj==1)||(bj==2)) )&&
		  (((bk==3)||(bk==4)||(bk==5)||(bk==1)||(bk==2))&& ((bj==1)||(bj==2)) )
#endif
		  ){
	    
	    typeofbc=6;
	  }
	  //  i-reflect(or AOS) AND j/k-inflow/outflow/periodic/reflect(or AOS)
	  else if(
#if(COMPDIM<3)
		  ((bj==3)||(bj==4)||(bj==5)||(bj==1)||(bj==2))&& ((bi==1)||(bi==2))
#else
		  (((bj==3)||(bj==4)||(bj==5)||(bj==1)||(bj==2))&& ((bi==1)||(bi==2)) )&&
		  (((bk==3)||(bk==4)||(bk==5)||(bk==1)||(bk==2))&& ((bi==1)||(bi==2)) )
#endif
		  ){
	    
	    typeofbc=7;
	  }
	  //  k-reflect(or AOS) AND j/i-inflow/outflow/periodic/reflect(or AOS)
	  else if(
#if(COMPDIM<3)
		  0
#else
		  (((bj==3)||(bj==4)||(bj==5)||(bj==1)||(bj==2))&& ((bk==1)||(bk==2)) )&&
		  (((bi==3)||(bi==4)||(bi==5)||(bi==1)||(bi==2))&& ((bk==1)||(bk==2)) )
#endif
		  ){
	    
	    typeofbc=8;
	  }
	  else{
	    fprintf(fail_file,"s[%d]: k=%d j=%d,i=%d quad5: No boundary case setup for bi: %d bj: %d bk: %d\n",ll,k,j,i,bi,bj,bk);
	    myexit(1);
	  }	  for(p=0;p<numcornerzones;p++){

	    // inflow or itself
	    ii=ix[corner][p][0][1];
	    jj=ix[corner][p][0][2];
	    kk=ix[corner][p][0][3];

	    // periodic
	    iii=ix[corner][p][1][1];
	    jjj=ix[corner][p][1][2];
	    kkk=ix[corner][p][1][3];

	    // outflow
	    iio=ix[corner][p][2][1];
	    jjo=ix[corner][p][2][2];
	    kko=ix[corner][p][2][3];

	    // reflective or AOS
	    iir=ix[corner][p][3][1];
	    jjr=ix[corner][p][3][2];
	    kkr=ix[corner][p][3][3];

	    if(typeofbc==1){
	      works[ll][kk][jj][ii]=sanal[ll][kk][jj][ii];
	    }
	    else if(typeofbc==2){
	      works[ll][kk][jj][ii]=works[ll][kkk][jjj][iii];
	    }
	    else if(typeofbc==3){
	      works[ll][kk][jj][ii]=works[ll][kk][jjo][ii];
	    }
	    else if(typeofbc==4){
	      works[ll][kk][jj][ii]=works[ll][kk][jj][iio];
	    }
	    else if(typeofbc==5){
	      works[ll][kk][jj][ii]=works[ll][kko][jj][ii];
	    }
	    else if(typeofbc==6){
	      works[ll][kk][jj][ii]=works[ll][kk][jjr][ii];
	    }
	    else if(typeofbc==7){
	      works[ll][kk][jj][ii]=works[ll][kk][jj][iir];
	    }
	    else if(typeofbc==8){
	      works[ll][kk][jj][ii]=works[ll][kkr][jj][ii];
	    }
	    else{
	      fprintf(fail_file,"s[%d]: k=%d j=%d,i=%d quad6: No boundary case setup for bi: %d bj: %d bk: %d\n",ll,k,j,i,bi,bj,bk);
	      myexit(1);
	    }
	  }
	}
      }// end of corners


#endif // end if do corner scalar bound
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
    
    for(l=1;l<=REALNUMVEC;l++){
      
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
	  wbound=1;
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

      if(wbound==2){
	dualloop=1; // for divb=0 fix since dependent on "nonlocal" variables
      }
      else dualloop=0;

      for(divbfix=0;divbfix<=dualloop;divbfix++){
	
      LOOPBOUNDRECT{ // if using 2 boundary zones, still only loop over first layer, assuming 2nd layer is calculatable
#if(!INNERBC )
	SKIPCONDITION;
#endif

#if(PUREBC==0)
	bct=bcv[wbound][1][k][j][i];
	bcdim=bcv[wbound][2][k][j][i];
	bcdir=bcv[wbound][3][k][j][i];

#else
	bct=purebccall(k,j,i,&bcdim,&bcdir);	
#endif	
	
	if(forcecheck){
	  if(bct==3) bct=4; // force dq or mdot to be outflowed(or other) if no assignments put in
	}
	
	if((bct>0)&&(bct<10)){ /* skip zone if computational zone or null zone(corner zones) */
	    
	  
	  
	  if(bcdim==1){ bcdir1=bcdir; bcdir2=0; bcdir3=0;}
	  else if(bcdim==2){ bcdir1=0; bcdir2=bcdir; bcdir3=0;}
	  else if(bcdim==3){ bcdir1=0; bcdir2=0; bcdir3=bcdir;}
	  
#include "boundpos.h"
	  
	  
	  //printf("bct: %d bcdim: %d bcdir: %d\n",bct,bcdim,bcdir);
	  //printf("i: %d j: %d k: %d ii: %d jj: %d kk: %d i2: %d j2: %d k2: %d iii: %d jjj: %d kkk: %d\n",i,j,k,ii,jj,kk,i2,j2,k2,iii,jjj,kkk);
	  
	  switch(bct){
	    
	  case 1:
	  case 2:	      // AOS arbitrary here
	    /* reflect and AOS */
	    if(bct==1) tempi=1;
	    else if(bct==2) tempi=-1;
	    else{
	      fprintf(fail_file,"case2 bct for vec\n");
	      myexit(6);
	    }
	    if(bcdim==1){
	      // correct for singular case
	      if(N1==1) jumpfactor=0;
	      else jumpfactor=1;

	      if(docom[1]){
		workv[ll][1][k][j][i+bcd1]=0; /* reflect inner vector */
		if(i+bcd2 < N1+N1BND) workv[ll][1][k][j][i+bcd2]=-workv[ll][1][kk][jj][ii+bcd1]; /* reflect outer vector */ // gets "2nd" velocity layer on right of grid
		if(i+bcd2-bcdir1 < N1+N1BND) workv[ll][1][k][j][i+bcd2-bcdir1]=-workv[ll][1][kk][jj][ii+bcd1+bcdir1*jumpfactor];
	      }
	      if(docom[2]){
		workv[ll][2][k][j][i]=workv[ll][2][kk][jj][ii];
		workv[ll][2][k][j][i-bcdir1]=workv[ll][2][kk][jj][ii+bcdir1*jumpfactor];
	      }
	      if(docom[3]){
		workv[ll][3][k][j][i]=tempi*workv[ll][3][kk][jj][ii]; /* same for DIM=2 or 3 */
		workv[ll][3][k][j][i-bcdir1]=tempi*workv[ll][3][kk][jj][ii+bcdir1*jumpfactor]; /* same for DIM=2 or 3 */
	      }
	    }
	    if(bcdim==2){
	      // correct for singular case
	      if(N2==1) jumpfactor=0;
	      else jumpfactor=1;

	      if(docom[2]){
		workv[ll][2][k][j+bcd1][i]=0; /* reflect inner vector */
		if(j+bcd2 < N2+N2BND) workv[ll][2][k][j+bcd2][i]=-workv[ll][2][kk][jj+bcd1][ii]; /* reflect outer vector */
		if(j+bcd2-bcdir2 < N2+N2BND) workv[ll][2][k][j+bcd2-bcdir2][i]=-workv[ll][2][kk][jj+bcd1+bcdir2*jumpfactor][ii]; // reflect outer vector
	      }
	      if(docom[1]){
		workv[ll][1][k][j][i]=workv[ll][1][kk][jj][ii];
		workv[ll][1][k][j-bcdir2][i]=workv[ll][1][kk][jj+bcdir2*jumpfactor][ii];
	      }
	      if(docom[3]){
		workv[ll][3][k][j][i]=tempi*workv[ll][3][kk][jj][ii]; /* same for DIM=2 or 3 */
		workv[ll][3][k][j-bcdir2][i]=tempi*workv[ll][3][kk][jj+bcdir2*jumpfactor][ii]; //same for DIM=2 or 3
	      }
	    }
	    if(bcdim==3){ /* never reached unless really in 3-d */
	      // correct for singular case
	      if(N3==1) jumpfactor=0;
	      else jumpfactor=1;
	      
	      if(docom[3]){
		workv[ll][3][k+bcd1][j][i]=0; /* reflect inner vector */
		if(k+bcd2 < N3+N3BND) workv[ll][3][k+bcd2][j][i]=-workv[ll][3][kk+bcd1][jj][ii]; /* reflect outer vector */
		if(k+bcd2-bcdir3 < N3+N3BND) workv[ll][3][k+bcd2-bcdir3][j][i]=-workv[ll][3][kk+bcd1+bcdir3*jumpfactor][jj][ii]; /* reflect outer vector */
	      }
	      if(docom[1]){
		workv[ll][1][k][j][i]=tempi*workv[ll][1][kk][jj][ii];
		workv[ll][1][k-bcdir3][j][i]=tempi*workv[ll][1][kk+bcdir3*jumpfactor][jj][ii];
	      }
	      if(docom[2]){
		workv[ll][2][k][j][i]=tempi*workv[ll][2][kk][jj][ii];
		workv[ll][2][k-bcdir3][j][i]=tempi*workv[ll][2][kk+bcdir3*jumpfactor][jj][ii];
	      }
	    }
	    break;
	  case 3:
	    if(divbfix==0){
	      /* Fix/Time vary: Dirichlet */
	      if( (wvec>0)||(wvec==-1)){
		for(m=1;m<=3;m++){
		  if(docom[m]){
		    workv[ll][m][k][j][i]=vanal[ll][m][k][j][i]; // v[-1] or v[N]
		    workv[ll][m][k2][j2][i2]=vanal[ll][m][k2][j2][i2]; // v[-2] or v[N+1]
		    // if on inner edge, also get element v[0](due to staggered grid)
		    if(m==bcdim){
		      if(bcdir==1){
			workv[ll][m][k+bcdir3][j+bcdir2][i+bcdir1]=vanal[ll][m][k+bcdir3][j+bcdir2][i+bcdir1]; // v[0] 
		      }
		    }
		  }
		}
	      }
	      else if(wvec==-2){
		fprintf(fail_file,"Case 3 bound with wvec==-2 has no definition\n");
		myexit(1);
	      }
	    }
	    break;
	  case 4:
	    if((wvec==-3)||(wvec==-4)) break; // this outflow is causing problems with wvec==-3!!
	    /* outflow */
	    // deal with asymmetry in velocity components on grid w.r.t. inner/outer edges
	    // deals also with inflow checking
	    if(bcdim==1){
	      // stuff independent of bcdir and/or need to do before other things(b2 for b1)
	      if(divbfix==0){
#if(HAWLEYBOUND==0)
		if(docom[3]){
		  workv[ll][3][k][j][i]=workv[ll][3][k][j][ii];
		  workv[ll][3][k][j][i-bcdir1]=workv[ll][3][k][j][i];
		}
		if(docom[2]){
		  // outer boundary zone1
		  workv[ll][2][k][j][i]=workv[ll][2][k][j][ii];
		  // outer boundary zone2
		  workv[ll][2][k][j][i-bcdir1]=workv[ll][2][k][j][i];
		}
#else // hawleybound==1
		if(docom[3]){
		  workv[ll][3][k][j][i]=0;
		  workv[ll][3][k][j][i-bcdir1]=0;
		}
		if(docom[2]){
		  // outer boundary zone1
		  workv[ll][2][k][j][i]=0;
		  // outer boundary zone2
		  workv[ll][2][k][j][i-bcdir1]=0;
		}
#endif // end if HK00 method
	      }
	      if(docom[1]){
		if(bcdir==1){
		  if(divbfix==0){
		    if(wbound==1){
		      if(INFLOWCHECKIX1&&(workv[ll][1][k][j][ii+1]>0.0)){
			workv[ll][1][k][j][ii]=0.0;
			workv[ll][1][k][j][i]=0.0;
			workv[ll][1][k][j][i-bcdir1]=0.0;
		      }
		      else{
			if(N1>1) workv[ll][1][k][j][ii]=workv[ll][1][k][j][ii+1];
			workv[ll][1][k][j][i]=workv[ll][1][k][j][ii];
			workv[ll][1][k][j][i-bcdir1]=workv[ll][1][k][j][i];
		      }// end if not inflow check and velocity
		    }// end if vel type bound
		  }
		  else{
		    if(wbound==2){// assumes 1&2 always together for B field!!!!!!!!!!!!!!  which is currently true.
		      // Br- : forced to have divb=0, so no choice on boundary
		      //workv[ll][1][k][j][ii]=b1m(workv[ll],k,j,ii);
		      workv[ll][1][k][j][ii]=b1(-1,workv[ll],k,j,ii);
		      // Br- : set such that the theta components are unchanged as 2 lines above
		      //workv[ll][1][k][j][i]=b1m(workv[ll],k,j,i);
		      workv[ll][1][k][j][i]=b1(-1,workv[ll],k,j,i);
		      // Br- : set such that the theta components are unchanged as 2 lines above
		      workv[ll][1][k][j][i-bcdir1]=b1(-1,workv[ll],k,j,i-bcdir1);
		      //workv[ll][1][k][j][i-bcdir1]=b1m(workv[ll],k,j,i-bcdir1);
		    }// endif magnetic field type bound
		  }// end else if divbfix==1
		}// endif bcdir==1
		else{
		  if(divbfix==0){
		    if(wbound==1){
		      if(INFLOWCHECKOX1&&(wbound==1)&&(workv[ll][1][k][j][ii]<0.0)){
			workv[ll][1][k][j][i]=0.0;
			workv[ll][1][k][j][i+1]=0.0;
		      }
		      else{
			workv[ll][1][k][j][i]=workv[ll][1][k][j][ii];
			workv[ll][1][k][j][i+1]=workv[ll][1][k][j][i];
		      }// end if not inflowcheck and velocity
		    }// end if vel type bound
		  }
		  else{
		    if(wbound==2){
		      
		      // Br+:Br(i+1,j): forced by divB=0 condition
		      workv[ll][1][k][j][i]=b1(1,workv[ll],k,j,i);
		      // Br+ : as above
		      workv[ll][1][k][j][i-bcdir1]=b1(1,workv[ll],k,j,i-bcdir1);
		    }// end if mag field type bound
		  }// end if divbfix==1
		}// end if bcdir==-1
	      }// end if docom[1]=1
	    }// end if bcdim==1
	    else if(bcdim==2){
	      if(divbfix==0){
		// stuff independent of bcdir and/or need to do before other things(b2 for b1)
#if(HAWLEYBOUND==0)
		if(docom[3]){
		  workv[ll][3][k][j][i]=workv[ll][3][k][jj][i];
		  workv[ll][3][k][j-bcdir2][i]=workv[ll][3][k][j][i];
		}
		if(docom[1]){ // assumes 1&2 always together for B field!!!!!!!!!!!!!!  which is currently true.
		  // outer boundary zone1
		  workv[ll][1][k][j][i]=workv[ll][1][k][jj][i];
		  // outer boundary zone2
		  workv[ll][1][k][j-bcdir2][i]=workv[ll][1][k][j][i];
		}
		
#else // if hawleybound==1
		
		if(docom[3]){
		  workv[ll][3][k][j][i]=0;
		  workv[ll][3][k][j-bcdir2][i]=0;
		}
		if(docom[1]){
		  // outer boundary zone1
		  workv[ll][1][k][j][i]=0;
		  // outer boundary zone2
		  workv[ll][1][k][j-bcdir2][i]=0;
		}
#endif // endif using HK00 method
	      }
	      if(docom[2]){
		if(bcdir==1){
		  if(divbfix==0){
		    if(wbound==1){
		      if(INFLOWCHECKIX2&&(workv[ll][2][k][jj+1][i]>0.0)){
			workv[ll][2][k][jj][i]=0.0;
			workv[ll][2][k][j][i]=0.0;
			workv[ll][2][k][j-bcdir2][i]=0.0;
		      }
		      else{
			if(N2>1) workv[ll][2][k][jj][i]=workv[ll][2][k][jj+1][i];
			workv[ll][2][k][j][i]=workv[ll][2][k][jj][i];
			workv[ll][2][k][j-bcdir2][i]=workv[ll][2][k][j][i];
		      }
		    }// end if vel type bound
		  }
		  else{
		    if(wbound==2){ // assumes 1&2 always together for B field!!!!!!!!!!!!!!  which is currently true.
		      // Btheta- : forced to have divb=0, so no choice on boundary
		      workv[ll][2][k][jj][i] = b2(-1,workv[ll],k,jj,i);
		      
		      // Btheta- : set such that the theta components are unchanged as lines above
		      workv[ll][2][k][j][i] = b2(-1,workv[ll],k,j,i);
		      
		      // Btheta- : set such that the theta components are unchanged as lines above
		      workv[ll][2][k][j-bcdir2][i] = b2(-1,workv[ll],k,j-bcdir2,i);
		    }// endif magnetic field type bound
		  }// end if divbfix==1
		} // end if bcdir==1
		else{
		  if(divbfix==0){
		    if(wbound==1){
		      if(INFLOWCHECKOX2&&(wbound==1)&&(workv[ll][2][k][jj][i]<0.0)){
			workv[ll][2][k][j][i]=0.0;
			workv[ll][2][k][j+1][i]=0.0;
		      }
		      else{
			workv[ll][2][k][j][i]=workv[ll][2][k][jj][i];
			workv[ll][2][k][j+1][i]=workv[ll][2][k][j][i];
		      }
		    }
		  }
		  else{
		    if(wbound==2){
		      // Btheta+:
		      workv[ll][2][k][j][i] = b2(1,workv[ll],k,j,i);
		      // Btheta+:
		      workv[ll][2][k][j+1][i] = b2(1,workv[ll],k,j+1,i);
		    }
		  }// end if divbfix==1
		} // end if bcdir==-1
	      }// end if docom[2]=1
	    }// end if bcdim=2
	    else if(bcdim==3){
	      if(divbfix==0){
#if(HAWLEYBOUND==0)
		if(docom[1]){
		  workv[ll][1][k][j][i]=workv[ll][1][kk][j][i];
		  workv[ll][1][k-bcdir3][j][i]=workv[ll][1][k][j][i];
		}
		if(docom[2]){
		  workv[ll][2][k][j][i]=workv[ll][2][kk][j][i];
		  workv[ll][2][k-bcdir3][j][i]=workv[ll][2][k][j][i];
		}
		
#else // if hawleybound==1
		
		if(docom[1]){
		  workv[ll][1][k][j][i]=0;
		  workv[ll][1][k-bcdir3][j][i]=0;
		}
		if(docom[2]){
		  workv[ll][2][k][j][i]=0;
		  workv[ll][2][k-bcdir3][j][i]=0;
		}
		
#endif // end if HK00 method
	      }
	      if(docom[3]){
		if(bcdir==1){
		  if(divbfix==0){
		    if(wbound==1){
		      if(INFLOWCHECKIX3&&(workv[ll][3][kk+1][j][i]>0.0)){
			workv[ll][3][kk][j][i]=0.0;
			workv[ll][3][k][j][i]=0.0;
			workv[ll][3][k-bcdir3][j][i]=0.0;
		      }
		      else{
			if(N3>1) workv[ll][3][kk][j][i]=workv[ll][3][kk+1][j][i];
			workv[ll][3][k][j][i]=workv[ll][3][kk][j][i];
			workv[ll][3][k-bcdir3][j][i]=workv[ll][3][k][j][i];
		      }
		    }// end if vel type bound
		  }
		  else{
		    if(wbound==2){// assumes 2&3 always together for B field!!!!!!!!!!!!!!  which is currently true in full 3D
		      // Bp- : forced to have divb=0, so no choice on boundary
		      workv[ll][3][kk][j][i]=b3(-1,workv[ll],kk,j,i);
		      // Bp- : set such that the theta components are unchanged as 2 lines above
		      workv[ll][3][kk][j][i]=b3(-1,workv[ll],k,j,i);
		      // Bp- : set such that the theta components are unchanged as 2 lines above
		      workv[ll][3][k-bcdir3][j][i]=b3(-1,workv[ll],k-bcdir3,j,i);
		    }// endif magnetic field type bound
		  }// end if divbfix==1
		}// endif bcdir==1
		else{
		  if(divbfix==0){
		    if(wbound==1){
		      if(INFLOWCHECKOX3&&(wbound==1)&&(workv[ll][3][kk][j][i]<0.0)){
			workv[ll][3][k][j][i]=0.0;
			workv[ll][3][k+1][j][i]=0.0;
		      }
		      else{
			workv[ll][3][k][j][i]=workv[ll][3][kk][j][i];
			workv[ll][3][k+1][j][i]=workv[ll][3][k][j][i];
		      }
		    }// end if vel type bound
		  }
		  else{
		    if(wbound==2){
		      
		      // Bp+:Bp(k+1,i,j): forced by divB=0 condition
		      workv[ll][3][k][j][i]=b3(1,workv[ll],k,j,i);
		      // Bp+ : as above
		      workv[ll][3][k-bcdir3][j][i]=b3(1,workv[ll],k-bcdir3,j,i);
		      
		    }// end if mag field type bound
		  }// end if divbfix==1
		}// end if bcdir==-1
	      }// docom[3]==1
	    }// bcdim==3
	    break;
	  case 5:
	    if(divbfix==0){
	      // periodic
	      for(m=1;m<=3;m++){
		if(docom[m]){
		  if(periodicx2special){
		    if((bcdim==2)&&(m==2)&&((wvec==-3)||(wvec==1)||(wvec==2)||(wvec==-1))){
		      // then on x2-edge
		      // - sign since metric doesn't flip sign
		      if(j2==N2+1){
			workv[ll][m][k2][j2][i2]=-workv[ll][m][kkk][jje2-bcdir2][iii]; // jje2-1=N2-1
		      }
		      else{
			workv[ll][m][k][j][i]=-workv[ll][m][kk][jje2][ii];
			workv[ll][m][k2][j2][i2]=-workv[ll][m][kkk][jjje2][iii];
		      }
		    }
		    else{
		      if((bcdim==2)&&((m==2)||(m==3))) dosign=-1; // applies to wvec==-4 only currently
		      else dosign=1; // otherwise either normal vector or other boundary
		      // then not on x2-edge, so form is same as normal case
		      //if(m==2){
		      //fprintf(stdout,"%d %d : %d %d %d = %d %d %d     %d %d %d = %d %d %d\n",ll,m,k,j,i,kk,jj,ii,k2,j2,i2,kkk,jjj,iii); fflush(stdout);
		      //}
		      workv[ll][m][k][j][i]=dosign*workv[ll][m][kk][jj][ii];
		      workv[ll][m][k2][j2][i2]=dosign*workv[ll][m][kkk][jjj][iii];
		    }
		  }
		  else{
		    workv[ll][m][k][j][i]=workv[ll][m][kk][jj][ii];
		    workv[ll][m][k2][j2][i2]=workv[ll][m][kkk][jjj][iii];

		  }
		}
	      }
	    }
	    break;
	  case 90:
	  case 91:
	  case 92:
	  case 93:
	  case 94:
	  case 95:
	  case 96:
	  case 97:
	  case 98:
	  case 99:
	    // no nothing
	    break;
	  default:
	    fprintf(fail_file,"bound.c: switch(bct) error, case given is: %d\n\r",bct);
	    myexit(1);
	  } // end switch
	} // end if bct>0 for vectors
      }// end over current vector
      }// end over divbfix loop        
    
#if(VECBCCORNER)


      // only 3D should hit this
      for(edge=0;edge<numedges;edge++){
	
	i=edgesi[edge];
	j=edgesj[edge];
	k=edgesk[edge];

#if(PUREBC==0)
	// multiply boundary type by 0 if not relevant type(i.e. only 2d ambiguity)
	bireal=bcv[wbound][1][k+edgegk[edge]][j+edgegj[edge]][i]*abs(edgegi[edge]); // i-boundary type
	bjreal=bcv[wbound][1][k+edgegk[edge]][j][i+edgegi[edge]]*abs(edgegj[edge]); // j-boundary type
	bkreal=bcv[wbound][1][k][j+edgegj[edge]][i+edgegi[edge]]*abs(edgegk[edge]); // k-boundary type
#else
	bireal=purebccall(k+edgegk[edge],j+edgegj[edge],i,&bcdim,&bcdir)*abs(edgegi[edge]); // i-boundary type
	bjreal=purebccall(k+edgegk[edge],j,i+edgegi[edge],&bcdim,&bcdir)*abs(edgegj[edge]); // j-boundary type
	bkreal=purebccall(k,j+edgegj[edge],i+edgegi[edge],&bcdim,&bcdir)*abs(edgegk[edge]); // k-boundary type
#endif

	if(bireal==0){
	  bi=bjreal; diri=2;
	  bj=bkreal; dirj=3;
	  startloop=edgesi[edge];
	}
	if(bjreal==0){
	  bi=bireal; diri=1;
	  bj=bkreal; dirj=3;
	  startloop=edgesj[edge];
	}
	if(bkreal==0){
	  bi=bireal; diri=1;
	  bj=bjreal; dirj=2;
	  startloop=edgesk[edge];
	}
	endloop=edgeend[edge];
	  
	
	if(
	   (bi==0)||(bj==0)
	   ){
	  fprintf(fail_file,"s[%d]: k=%d,j=%d,i=%d quad7: undefined boundary bi: %d bj: %d\n",ll,k,j,i,bi,bj);
	  myexit(1);
	}
	else if(
		(bi>=90)||(bj>=90)
		){
	}
	else{ // if good edge
	  
	  if(// i-inflow AND j-inflow or if p-i or i-p 
	     ((bi==3)&&(bj==3))||((bi==5)&&(bj==3))||((bi==3)&&(bj==5))
	     ){
	    
	    typeofbc=1;
	  }
	  // i-periodic AND j-periodic
	  else if(
		  (bi==5)&&(bj==5)
		  ){
	    
	    typeofbc=2;
	  }
	  // i/k-outflow/inflow/periodic AND j-outflow (if o-o, same effect as next case)
	  else if(
		  ((bi==4)||(bi==3)||(bi==5))&&(bj==4) // so "j" is either k or j real dir
		  ){
	    if(dirj==2) typeofbc=3; // so j is j 
	    else typeofbc=4; // so j is k
	  }
	  // j/k-outflow/inflow/periodic AND i-outflow
	  else if(
		  ((bj==4)||(bj==3)||(bj==5))&&(bi==4) // so "i" is either i or j real dir
		  ){
	    if(diri==1) typeofbc=5; // so i is i
	    else typeofbc=6; // so i is j
	  }
	  // i/k-inflow/outflow/periodic/reflect(or AOS) AND j-reflect(or AOS)
	  else if(
		  ((bi==3)||(bi==4)||(bi==5)||(bi==1)||(bi==2))&& ((bj==1)||(bj==2))
		  ){
	    if(dirj==2) typeofbc=7; // so j is j 
	    else typeofbc=8; // so j is k
	  }
	  //  i-reflect(or AOS) AND j/k-inflow/outflow/periodic/reflect(or AOS)
	  else if(
		  ((bj==3)||(bj==4)||(bj==5)||(bj==1)||(bj==2))&& ((bi==1)||(bi==2))
		  ){
	    if(diri==1) typeofbc=9; // so i is i
	    else typeofbc=10; // so i is j
	  }
	  else{
	    fprintf(fail_file,"s[%d]: k=%d j=%d,i=%d quad8: No boundary case setup for bi: %d bj: %d\n",ll,k,j,i,bi,bj);
	    myexit(1);
	  }


	  loopi=!abs(edgegi[edge]);
	  loopj=!abs(edgegj[edge]);
	  loopk=!abs(edgegk[edge]);

	  for(p=0;p<numedgezones;p++){
	    
	    // inflow or itself
	    ii=edgeix[edge][p][0][1];
	    jj=edgeix[edge][p][0][2];
	    kk=edgeix[edge][p][0][3];
	    
	    // periodic
	    iii=edgeix[edge][p][1][1];
	    jjj=edgeix[edge][p][1][2];
	    kkk=edgeix[edge][p][1][3];
	    
	    // outflow
	    iio=edgeix[edge][p][2][1];
	    jjo=edgeix[edge][p][2][2];
	    kko=edgeix[edge][p][2][3];
	    
	    // reflective or AOS
	    iir=edgeix[edge][p][3][1];
	    jjr=edgeix[edge][p][3][2];
	    kkr=edgeix[edge][p][3][3];
	    
	    for(loopvar=startloop;loopvar<=endloop;loopvar++){ // perform loop down non-bounded direction
	    
	      //printf("hit: %d %d %d %d %d    %d %d %d\n",edge,p,loopvar,startloop, endloop, kk,jj,ii);
	      
	      if(typeofbc==1){
		workv[ll][1][kk][jj][ii]=vanal[ll][1][kk][jj][ii];
		workv[ll][2][kk][jj][ii]=vanal[ll][2][kk][jj][ii];
		workv[ll][3][kk][jj][ii]=vanal[ll][3][kk][jj][ii];
	      }
	      else if(typeofbc==2){
		workv[ll][1][kk][jj][ii]=workv[ll][1][kkk][jjj][iii];
		workv[ll][2][kk][jj][ii]=workv[ll][2][kkk][jjj][iii];
		workv[ll][3][kk][jj][ii]=workv[ll][3][kkk][jjj][iii];
	      }
	      
	      else if((typeofbc==3)||(typeofbc==6)){
		workv[ll][1][kk][jj][ii]=workv[ll][1][kk][jjo][ii];
		workv[ll][2][kk][jj][ii]=workv[ll][2][kk][jjo][ii];
		workv[ll][3][kk][jj][ii]=workv[ll][3][kk][jjo][ii];
	      }
	      else if(typeofbc==4){
		workv[ll][1][kk][jj][ii]=workv[ll][1][kko][jj][ii];
		workv[ll][2][kk][jj][ii]=workv[ll][2][kko][jj][ii];
		workv[ll][3][kk][jj][ii]=workv[ll][3][kko][jj][ii];
	      }
	      else if(typeofbc==5){
		workv[ll][1][kk][jj][ii]=workv[ll][1][kk][jj][iio];
		workv[ll][2][kk][jj][ii]=workv[ll][2][kk][jj][iio];
		workv[ll][3][kk][jj][ii]=workv[ll][3][kk][jj][iio];
	      }
	      else if((typeofbc==7)||(typeofbc==10)){ // j-reflect
		if((edge==0)||(edge==4)||(edge==5)||(edge==8)){ // then fix applies like this:
		  workv[ll][2][kk][jj][ii]=-workv[ll][2][kk][jjr+edgegj[edge]][ii];
		  if(jj==-1){
		    workv[ll][2][kk][jj+1][ii]=0;
		  }
		}
		else if((edge==2)||(edge==6)||(edge==7)||(edge==10)){  // fix applies like this:
		  if(jj==N2){
		    workv[ll][2][kk][jj][ii]=0; //"corner" zone is really on edge 
		  }
		  else{
		    workv[ll][2][kk][jj][ii]=-workv[ll][2][kk][jjr-edgegj[edge]][ii];
		  }
		}
		else{
		  fprintf(fail_file,"Got to undefined spot in edge bc loop: edge: %d typeofbc: %d\n",edge,typeofbc);
		  myexit(1);
		}
		workv[ll][1][kk][jj][ii]=workv[ll][1][kk][jjr][ii];
		workv[ll][3][kk][jj][ii]=workv[ll][3][kk][jjr][ii];
	      }
	      else if(typeofbc==8){ // k-reflect
		if((edge==0)||(edge==1)||(edge==2)||(edge==3)){ // then fix applies like this:
		  workv[ll][3][kk][jj][ii]=-workv[ll][3][kkr+edgegk[edge]][jj][ii];
		  if(kk==-1){
		    workv[ll][3][kk+1][jj][ii]=0;
		  }
		}
		else if((edge==8)||(edge==9)||(edge==10)||(edge==11)){  // fix applies like this:
		  if(kk==N3){ // right on edge
		    workv[ll][3][kk][jj][ii]=0; //"corner" zone is really on edge 
		  }
		  else{
		    workv[ll][3][kk][jj][ii]=-workv[ll][3][kkr-edgegk[edge]][jj][ii];
		  }
		}
		else{
		  fprintf(fail_file,"Got to undefined spot in edge bc loop: edge: %d typeofbc: %d\n",edge,typeofbc);
		  myexit(1);
		}
		workv[ll][1][kk][jj][ii]=workv[ll][1][kkr][jj][ii];
		workv[ll][2][kk][jj][ii]=workv[ll][2][kkr][jj][ii];
	      }
	      else if(typeofbc==9){ // i-reflect
		if((edge==3)||(edge==4)||(edge==7)||(edge==11)){ // then fix applies like this:
		  workv[ll][1][kk][jj][ii]=-workv[ll][1][kk][jj][iir+edgegi[edge]];
		  if(ii==-1){
		    workv[ll][1][kk][jj][ii+1]=0;
		  }
		}
		else if((edge==1)||(edge==5)||(edge==6)||(edge==9)){  // fix applies like this:
		  if(ii==N1){ // right on edge
		    workv[ll][1][kk][jj][ii]=0; //"corner" zone is really on edge 
		  }
		  else{
		    workv[ll][1][kk][jj][ii]=-workv[ll][1][kk][jj][iir-edgegi[edge]];
		  }
		}
		else{
		  fprintf(fail_file,"Got to undefined spot in edge bc loop: edge: %d typeofbc: %d\n",edge,typeofbc);
		  myexit(1);
		}
		workv[ll][2][kk][jj][ii]=workv[ll][2][kk][jj][iir];
		workv[ll][3][kk][jj][ii]=workv[ll][3][kk][jj][iir];
	      }
	      else{
		fprintf(fail_file,"s[%d]: k=%d j=%d,i=%d quad9: No boundary case setup for bi: %d bj: %d\n",ll,k,j,i,bi,bj);
		myexit(1);
	      }

	      // update position to next zone over

	      // inflow or itself
	      ii+=loopi;
	      jj+=loopj;
	      kk+=loopk;
	      
	      // periodic
	      iii+=loopi;
	      jjj+=loopj;
	      kkk+=loopk;
	      
	      // outflow
	      iio+=loopi;
	      jjo+=loopj;
	      kko+=loopk;
	      
	      // reflective or AOS
	      iir+=loopi;
	      jjr+=loopj;
	      kkr+=loopk;

	    } // over zones 
	  }// end loop down the tube of zones
	  
	  if(wbound==2){
	    // now fix for divb constraint(must do after complete upper loop so side zones are defined!)
	    // can do edges indpendently since they are nonlocal to eachother(so no dependence)
	    
	    for(p=0;p<numedgezones;p++){
	      
	      
	      // inflow or itself
	      ii=edgeix[edge][p][0][1];
	      jj=edgeix[edge][p][0][2];
	      kk=edgeix[edge][p][0][3];
	      
	      for(loopvar=startloop;loopvar<=endloop;loopvar++){ // perform loop down non-bounded direction
	      
		if((typeofbc==3)||(typeofbc==6)||(typeofbc==7)||(typeofbc==10)||(typeofbc==5)||(typeofbc==9)||(typeofbc==4)||(typeofbc==8)){
		  if((typeofbc==3)||(typeofbc==6)||(typeofbc==7)||(typeofbc==10)){// j-outflow/reflect: same as scalar
		    if((typeofbc==3)||(typeofbc==6)){ // outflow
		      othercomp=2;
		    }
		    else{ // reflect
		      if(typeofbc==7) othercomp=diri; // k
		      else if(typeofbc==10) othercomp=dirj; // i
		    }
		  }
		  else if((typeofbc==5)||(typeofbc==9)){// i-outflow/reflect: same as scalar
		    // then edge must = 1,3,4,5,6,7,9,or 11
		    if(typeofbc==5){
		      othercomp=1;
		    }
		    else{
		      othercomp=dirj;
		    }
		  }
		  else if((typeofbc==4)||(typeofbc==8)){// k-outflow/reflect: same as scalar
		    if(typeofbc==4){
		      othercomp=3;
		    }
		    else{
		      othercomp=diri;
		    }
		  }
		  // now do the zone(s)
		  if(othercomp==1){
		    gval=edgegi[edge];
		    if(
		       (((edge==8)||(edge==9)||(edge==10)||(edge==11))&&(kk-edgesk[edge]!=0))||
		       (((edge==2)||(edge==6)||(edge==7)||(edge==10))&&(jj-edgesj[edge]!=0))
			) continue;
		  }		      
		  if(othercomp==2){
		    gval=edgegj[edge];
		    if(
		       (((edge==1)||(edge==5)||(edge==6)||(edge==9))&&(ii-edgesi[edge]!=0))||
		       (((edge==8)||(edge==9)||(edge==10)||(edge==11))&&(kk-edgesk[edge]!=0))
			 ) continue;
		  }
		  if(othercomp==3){
		    gval=edgegk[edge];
		    if(
		       (((edge==1)||(edge==5)||(edge==6)||(edge==9))&&(ii-edgesi[edge]!=0))||
		       (((edge==2)||(edge==6)||(edge==7)||(edge==10))&&(jj-edgesj[edge]!=0))
			) continue;
		  }
		  if(ii==N1){
		    workv[ll][othercomp][kk][jj][ii-1]=bdivb(othercomp,-gval,workv[ll],kk,jj,ii-1);
		  }
		  else if(ii==-1){
		    workv[ll][othercomp][kk][jj][ii+1]=bdivb(othercomp,-gval,workv[ll],kk,jj,ii+1);
		  }
		  if(jj==N2){
		    workv[ll][othercomp][kk][jj-1][ii]=bdivb(othercomp,-gval,workv[ll],kk,jj-1,ii);
		  }
		  else if(jj==-1){
		    workv[ll][othercomp][kk][jj+1][ii]=bdivb(othercomp,-gval,workv[ll],kk,jj+1,ii);
		  }
		  if(kk==-1){
		    workv[ll][othercomp][kk+1][jj][ii]=bdivb(othercomp,-gval,workv[ll],kk+1,jj,ii);
		  }
		  else if(kk==N3){
		    workv[ll][othercomp][kk-1][jj][ii]=bdivb(othercomp,-gval,workv[ll],kk-1,jj,ii);
		  }
		  workv[ll][othercomp][kk][jj][ii]=bdivb(othercomp,-gval,workv[ll],kk,jj,ii); // must come last
		  //}
		}
		else if(typeofbc==1){ // inflow: all same as scalar
		  // no fix needed
		}
		else if(typeofbc==2){ // periodic: all same as scalar
		  // no fix needed
		}
		else{
		  fprintf(fail_file,"v[%d]: k=%d j=%d,i=%d quad10: No boundary case setup for bi: %d bj: %d bk: %d\n",ll,k,j,i,bi,bj,bk);
		  myexit(1);
		}
		// update position to next zone over
		ii+=loopi;
		jj+=loopj;
		kk+=loopk;

	      }// over zones in a slice of tube
	    }// over tube
	  }// end if ll==2
	}// end else if good edge(i.e. assignable boundary type)
      }// end over edges

    
      // now do corner zones(assumes other boundary zones already done!)
      // normal boundary code takes care of 0 zone issue(due to staggered grid)
      
      for(corner=0;corner<numcorners;corner++){
	
	i=si[corner];
	j=sj[corner];
	k=sk[corner];

	// assume if N?OFF==0, then not relevant boundary, -100 signifies irrelevant

	if(N1NOT1){// i-boundary type
#if(PUREBC==0)
	  bi=bcv[wbound][1][k+gk[corner]][j+gj[corner]][i];
#else
	  bi=purebccall(k+gk[corner],j+gj[corner],i,&bcdim,&bcdir);
#endif
	}
	else{
	  bi=5; // periodic like if no direction in this direction
	}

	if(N2NOT1){// j-boundary type
#if(PUREBC==0)
	  bj=bcv[wbound][1][k+gk[corner]][j][i+gi[corner]];
#else
	  bj=purebccall(k+gk[corner],j,i+gi[corner],&bcdim,&bcdir);
#endif
	}
	else{
	  bj=5;
	}

	if((COMPDIM<3)||(N3NOT1==0)){
	  bk=5;
	}
	else{
#if(PUREBC==0)
	  bk=bcv[wbound][1][k][j+gj[corner]][i+gi[corner]];	  
#else
	  bk=purebccall(k,j+gj[corner],i+gi[corner],&bcdim,&bcdir);	  
#endif
	}


	if(
#if(COMPDIM<3)
	   (bi==0)||(bj==0)
#else
	   (bi==0)||(bj==0)||(bk==0)
#endif
	   ){
	  fprintf(fail_file,"v[%d]: k=%d,j=%d,i=%d quad11: undefined boundary bi: %d bj: %d bk: %d\n",ll,k,j,i,bi,bj,bk);
	  myexit(1);
	}
	else if(
#if(COMPDIM<3)
		(bi>=90)||(bj>=90)
#else
		(bi>=90)||(bj>=90)||(bk>=90)
#endif
		){
	}
	else{ // if good corner
	  
	  if(// i-inflow AND j-inflow or if p-i or i-p 
#if(COMPDIM<3)
	     ((bi==3)&&(bj==3))||((bi==5)&&(bj==3))||((bi==3)&&(bj==5))
#else
	     ((bi==3)&&((bj==3)||(bj==5))&&((bk==3)||(bk==5)))||
	     ((bj==3)&&((bi==3)||(bi==5))&&((bk==3)||(bk==5)))||
	     ((bk==3)&&((bj==3)||(bj==5))&&((bi==3)||(bi==5)))
#endif
	     ){
	    
	    typeofbc=1;
	  }
	  // i-periodic AND j-periodic
	  else if(
#if(COMPDIM<3)
		  (bi==5)&&(bj==5)
#else
		  (bi==5)&&(bj==5)&&(bk==5)
#endif
		  ){
	    
	    typeofbc=2;
	  }
	  // i/k-outflow/inflow/periodic AND j-outflow (if o-o, same effect as next case)
	  else if(
#if(COMPDIM<3)
		  ((bi==4)||(bi==3)||(bi==5))&&(bj==4)
#else
		  (((bi==4)||(bi==3)||(bi==5))&&(bj==4))&&
		  (((bk==4)||(bk==3)||(bk==5))&&(bj==4))
#endif	  
		  ){
	    typeofbc=3;
	  }
	  // j/k-outflow/inflow/periodic AND i-outflow
	  else if(
#if(COMPDIM<3)
		  ((bj==4)||(bj==3)||(bj==5))&&(bi==4)
#else
		  (((bj==4)||(bj==3)||(bj==5))&&(bi==4))&&
		  (((bk==4)||(bk==3)||(bk==5))&&(bi==4))
#endif
		  ){
	    
	    typeofbc=4;
	  }
	  // j/i-outflow/inflow/periodic AND k-outflow
	  else if(
#if(COMPDIM<3)
		  0
#else
		  (((bj==4)||(bj==3)||(bj==5))&&(bk==4))&&
		  (((bi==4)||(bi==3)||(bi==5))&&(bk==4))
#endif
		  ){
	    
	    typeofbc=5;
	  }
	  // i/k-inflow/outflow/periodic/reflect(or AOS) AND j-reflect(or AOS)
	  else if(
#if(COMPDIM<3)
		  ((bi==3)||(bi==4)||(bi==5)||(bi==1)||(bi==2))&& ((bj==1)||(bj==2))
#else
		  (((bi==3)||(bi==4)||(bi==5)||(bi==1)||(bi==2))&& ((bj==1)||(bj==2)) )&&
		  (((bk==3)||(bk==4)||(bk==5)||(bk==1)||(bk==2))&& ((bj==1)||(bj==2)) )
#endif
		  ){
	    
	    typeofbc=6;
	  }
	  //  i-reflect(or AOS) AND j/k-inflow/outflow/periodic/reflect(or AOS)
	  else if(
#if(COMPDIM<3)
		  ((bj==3)||(bj==4)||(bj==5)||(bj==1)||(bj==2))&& ((bi==1)||(bi==2))
#else
		  (((bj==3)||(bj==4)||(bj==5)||(bj==1)||(bj==2))&& ((bi==1)||(bi==2)) )&&
		  (((bk==3)||(bk==4)||(bk==5)||(bk==1)||(bk==2))&& ((bi==1)||(bi==2)) )
#endif
		  ){
	    
	    typeofbc=7;
	  }
	  //  k-reflect(or AOS) AND j/i-inflow/outflow/periodic/reflect(or AOS)
	  else if(
#if(COMPDIM<3)
		  0
#else
		  (((bj==3)||(bj==4)||(bj==5)||(bj==1)||(bj==2))&& ((bk==1)||(bk==2)) )&&
		  (((bi==3)||(bi==4)||(bi==5)||(bi==1)||(bi==2))&& ((bk==1)||(bk==2)) )
#endif
		  ){
	    
	    typeofbc=8;
	  }
	  else{
	    fprintf(fail_file,"v[%d]: k=%d j=%d,i=%d quad12: No boundary case setup for bi: %d bj: %d bk: %d\n",ll,k,j,i,bi,bj,bk);
	    myexit(1);
	  }
	  for(p=0;p<numcornerzones;p++){

	    // inflow or itself
	    ii=ix[corner][p][0][1];
	    jj=ix[corner][p][0][2];
	    kk=ix[corner][p][0][3];

	    // periodic
	    iii=ix[corner][p][1][1];
	    jjj=ix[corner][p][1][2];
	    kkk=ix[corner][p][1][3];

	    // outflow
	    iio=ix[corner][p][2][1];
	    jjo=ix[corner][p][2][2];
	    kko=ix[corner][p][2][3];

	    // reflective or AOS
	    iir=ix[corner][p][3][1];
	    jjr=ix[corner][p][3][2];
	    kkr=ix[corner][p][3][3];

	    if(typeofbc==1){ // inflow: all same as scalar
	      workv[ll][1][kk][jj][ii]=vanal[ll][1][kk][jj][ii];
	      workv[ll][2][kk][jj][ii]=vanal[ll][2][kk][jj][ii];
	      workv[ll][3][kk][jj][ii]=vanal[ll][3][kk][jj][ii];
	    }
	    else if(typeofbc==2){ // periodic: all same as scalar
	      workv[ll][1][kk][jj][ii]=workv[ll][1][kkk][jjj][iii];
	      workv[ll][2][kk][jj][ii]=workv[ll][2][kkk][jjj][iii];
	      workv[ll][3][kk][jj][ii]=workv[ll][3][kkk][jjj][iii];
	    }
	    else if(typeofbc==3){// j-outflow: same as scalar
	      workv[ll][1][kk][jj][ii]=workv[ll][1][kk][jjo][ii];
	      workv[ll][2][kk][jj][ii]=workv[ll][2][kk][jjo][ii];
	      workv[ll][3][kk][jj][ii]=workv[ll][3][kk][jjo][ii];
	    }
	    else if(typeofbc==4){// i-outflow: same as scalar
	      workv[ll][1][kk][jj][ii]=workv[ll][1][kk][jj][iio];
	      workv[ll][2][kk][jj][ii]=workv[ll][2][kk][jj][iio];
	      workv[ll][3][kk][jj][ii]=workv[ll][3][kk][jj][iio];
	    }
	    else if(typeofbc==5){// k-outflow: same as scalar
	      workv[ll][1][kk][jj][ii]=workv[ll][1][kko][jj][ii];
	      workv[ll][2][kk][jj][ii]=workv[ll][2][kko][jj][ii];
	      workv[ll][3][kk][jj][ii]=workv[ll][3][kko][jj][ii];
	    }
	    else if(typeofbc==6){ // j-reflective: not same as scalar for directional vector
	      if((corner==0)||(corner==1)||(corner==4)||(corner==5)){ // then fix applies like this:
		workv[ll][2][kk][jj][ii]=-workv[ll][2][kk][jjr+gj[corner]][ii];
		if(jj==-1){
		  workv[ll][2][kk][jj+1][ii]=0; //"corner" zone is really on edge and not fixed correctly in normal bc routine since this is dual bc ambiguous value
		}
	      }
	      else{  // fix applies like this:
		if(jj==N2){
		  workv[ll][2][kk][jj][ii]=0; //"corner" zone is really on edge 
		}
		else{
		  workv[ll][2][kk][jj][ii]=-workv[ll][2][kk][jjr-gj[corner]][ii];
		}
	      }
	      workv[ll][1][kk][jj][ii]=workv[ll][1][kk][jjr][ii];
	      workv[ll][3][kk][jj][ii]=workv[ll][3][kk][jjr][ii];
	    }
	    else if(typeofbc==7){ // i-reflective: not same as scalar for directional vector
	      if((corner==0)||(corner==3)||(corner==4)||(corner==7)){ // then fix applies like this:
		workv[ll][1][kk][jj][ii]=-workv[ll][1][kk][jj][iir+gi[corner]];
		if(ii==-1){
		  workv[ll][1][kk][jj][ii+1]=0; // also get this zone
		}
	      }
	      else{
		if(ii==N1){
		  workv[ll][1][kk][jj][ii]=0;
		}
		else{
		  workv[ll][1][kk][jj][ii]=-workv[ll][1][kk][jj][iir-gi[corner]];
		}
	      }
	      workv[ll][2][kk][jj][ii]=workv[ll][2][kk][jj][iir];
	      workv[ll][3][kk][jj][ii]=workv[ll][3][kk][jj][iir];
	    }
	    else if(typeofbc==8){ // k-reflective: not same as scalar for directional vector
	      if((corner==0)||(corner==1)||(corner==2)||(corner==3)){ // then fix applies like this:
		workv[ll][3][kk][jj][ii]=-workv[ll][3][kkr+gk[corner]][jj][ii];
		if(kk==-1){
		  workv[ll][3][kk+1][jj][ii]=0;
		}
	      }
	      else{
		if(kk==N3){
		  workv[ll][3][kk][jj][ii]=0;
		}
		else{
		  workv[ll][3][kk][jj][ii]=-workv[ll][3][kkr-gk[corner]][jj][ii];
		}
	      }
	      workv[ll][1][kk][jj][ii]=workv[ll][1][kkr][jj][ii];
	      workv[ll][2][kk][jj][ii]=workv[ll][2][kkr][jj][ii];
	    }
	    else{
	      fprintf(fail_file,"v[%d]: k=%d j=%d,i=%d quad13: No boundary case setup for bi: %d bj: %d bk: %d\n",ll,k,j,i,bi,bj,bk);
	      myexit(1);
	    }
	  }  // over zones in this corner
	

	  if(wbound==2){
	    // now fix for divb constraint(must do after complete upper loop so side zones are defined!)
	    // can do corners indpendently since they are nonlocal to eachother(so no dependence)
	    for(p=0;p<numcornerzones;p++){
	      
	      // inflow or itself
	      ii=ix[corner][p][0][1];
	      jj=ix[corner][p][0][2];
	      kk=ix[corner][p][0][3];
	      
	      if((typeofbc==3)||(typeofbc==6)||(typeofbc==4)||(typeofbc==7)||(typeofbc==5)||(typeofbc==8)){
		if((typeofbc==3)||(typeofbc==6)){// j-outflow/reflect: same as scalar
		  if(typeofbc==3){
		    othercomp=2;
		  }
		  else{
		    if((bk!=4)||(N3==1)) othercomp=1; // should be fine if 1-dir is reflect also, should be self-consistent
		    else othercomp=3;
		  }
		}
		else if((typeofbc==4)||(typeofbc==7)){// i-outflow/reflect: same as scalar
		  if(typeofbc==4){
		    othercomp=1;
		  }
		  else{
		    if((bk!=4)||(N3==1)) othercomp=2; // should be fine if 2-dir is reflect also, should be self-consistent
		    else othercomp=3;
		  }
		}
		else if((typeofbc==5)||(typeofbc==8)){// k-outflow/reflect: same as scalar
		  if(typeofbc==5){
		    othercomp=3;
		  }
		  else{
		    if((bj!=4)||(N2==1)) othercomp=1; // should be fine if 1-dir is reflect also, should be self-consistent
		    else othercomp=2;
		  }
		}
		// now must avoid different invalid memory regions
		if(othercomp==1){
		  gval=gi[corner];
		  if( (((corner==3)||(corner==2)||(corner==7)||(corner==6))&&(jj-sj[corner]!=0))||
		      (((corner==4)||(corner==5)||(corner==6)||(corner==7))&&(kk-sk[corner]!=0))
		      ) continue;
		}
		if(othercomp==2){
		  gval=gj[corner];
		  if( (((corner==1)||(corner==2)||(corner==5)||(corner==6))&&(ii-si[corner]!=0))||
		      (((corner==4)||(corner==5)||(corner==6)||(corner==7))&&(kk-sk[corner]!=0))
		      ) continue;
		}
		if(othercomp==3){
		  gval=gk[corner];
		  if( (((corner==3)||(corner==2)||(corner==7)||(corner==6))&&(jj-sj[corner]!=0))||
		      (((corner==1)||(corner==5)||(corner==6)||(corner==2))&&(ii-si[corner]!=0))
		      ) continue;
		}
		// extend fix since this other zone used changed corner zone
		if(ii==N1){
		  workv[ll][othercomp][kk][jj][ii-1]=bdivb(othercomp,-gval,workv[ll],kk,jj,ii-1);
		}
		else if(ii==-1){
		  workv[ll][othercomp][kk][jj][ii+1]=bdivb(othercomp,-gval,workv[ll],kk,jj,ii+1);
		}
		if(jj==N2){
		  workv[ll][othercomp][kk][jj-1][ii]=bdivb(othercomp,-gval,workv[ll],kk,jj-1,ii);
		}
		else if(jj==-1){
		  workv[ll][othercomp][kk][jj+1][ii]=bdivb(othercomp,-gval,workv[ll],kk,jj+1,ii);
		}
		if(kk==N3){
		  workv[ll][othercomp][kk-1][jj][ii]=bdivb(othercomp,-gval,workv[ll],kk-1,jj,ii);
		}
		else if(kk==-1){
		  workv[ll][othercomp][kk+1][jj][ii]=bdivb(othercomp,-gval,workv[ll],kk+1,jj,ii);
		}
		workv[ll][othercomp][kk][jj][ii]=bdivb(othercomp,-gval,workv[ll],kk,jj,ii); // must come last
	      }
	      else if(typeofbc==1){ // inflow: all same as scalar
		// no fix needed
	      }
	      else if(typeofbc==2){ // periodic: all same as scalar
		// no fix needed
	      }
	      else{
		fprintf(fail_file,"v[%d]: k=%d j=%d,i=%d quad14: No boundary case setup for bi: %d bj: %d bk: %d\n",ll,k,j,i,bi,bj,bk);
		myexit(1);
	      }
	    }// end over corner zones
	  }// end if ll==2
	}
      }// end over corners

#endif



    
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMVEC;
    }// end over vectors
  }// endif vectors to be done
  
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


int purebccall(int k, int j, int i,int*bcdim,int*bcdir)
{
  int bct;

  if(i>=N1){
    bct=PUREBC;
    *bcdim=1;
    *bcdir=-1;
  }
  else if(i<0){
    bct=PUREBC;
    *bcdim=1;
    *bcdir=1;
  }
  else if(j>=N2){
    bct=PUREBC;
    *bcdim=2;
    *bcdir=-1;
  }
  else if(j<0){
    bct=PUREBC;
    *bcdim=2;
    *bcdir=1;
  }
  else if(k>=N3){
    bct=PUREBC;
    *bcdim=3;
    *bcdir=-1;
  }
  else if(k<0){
    bct=PUREBC;
    *bcdim=3;
    *bcdir=1;
  }
  else{
    bct=0;
    *bcdim=0;
    *bcdir=0;
  }
  return(bct);
}
