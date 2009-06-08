#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif


#define  WIDTH 20

// Sod Shock Tube dump at specific time
// Fill <V/r/P/e>anal[i] with solution to sod shock problem at time t
// set L[1][1]=0 and L[2][1]=1 in init.c
void sodsol(int calltype)
{
  SFTYPE as_m,as_r,mu,rr,Pr,Vr,Pm,rl,Pl,Vl,Vm,rho2,Vs,rm;
  SFTYPE Bxl,Byl,Bzl,Bxr,Byr,Bzr,Vxl,Vyl,Vzl,Vxr,Vyr,Vzr;
  SFTYPE as0;
  SFTYPE regions[4];
  int i,j,k,l,si;
  static FILE*analyticout;
  static int firstsolve=1;
  char temps[100];
  char filename[100];
  SFTYPE ftemp;
  SFTYPE pos;
  int hit;
  SFTYPE SMOOTHER[WIDTH*2+1];
  static int CONVOLVESOD, MAGSOD,MAGDIR,MAGSODTYPE;



  // some global problem parameters
  CONVOLVESOD =0;
  
  MAGSOD= 1; // whether do magnetic version or not
  MAGDIR =1; // test in x1 or x2 or x3 direction ( need to setup x3 dir)
  MAGSODTYPE =6 ;// see below for what the types are



  if(calltype==100){
    // just setup problem

    if(MAGSOD){
      mag=1;
      switch(MAGSODTYPE){
      case 1:
	tf=0.1;// 1=S&N92 t=0.1; // he uses L=800 and tf=80
	gam=2.0;
	break;
      case 2:
	tf=0.08;// 2=RJ95-1a t=.08
	gam=5.0/3.0;
	break;
      case 3:
	tf=0.03;// 3=RJ95-1b t=.03
	gam=5.0/3.0;
	break;
      case 4:
	tf=0.2;// 4=RJ95-2a t=.2
	gam=5.0/3.0;
	break;
      case 5:
	tf=0.035;// 5=RJ95-2b t=.035
	gam=5.0/3.0;
	break;
      case 6:
	tf=0.01;// 6-RJ95-3a t=.01
	gam=5.0/3.0;
	break;
      case 7:
	tf=0.1;// 7-RJ95-3b t=.1
	gam=5.0/3.0;
	break;
      case 8:
	tf=0.15;// 8-RJ95-4a t=.15
	gam=5.0/3.0;
	break;
      case 9:
	tf=0.15;// 9-RJ95-4b t=.15
	gam=5.0/3.0;
	break;
      case 10:
	tf=0.15;// 10-RJ95-4c t=.15
	gam=5.0/3.0;
	break;
      case 11:
	tf=0.16;// 11-RJ95-4d t=.16
	gam=5.0/3.0;
	break;
      case 12:
	tf=0.1; // 12-RJ95-5a t=.1
	gam=5.0/3.0;
	break;
      case 13:
	tf=0.16;// 13-RJ95-5b t=.16
	gam=5.0/3.0;
	break;
      default:
	tf=0.1;
	gam=5.0/3.0;
	break;
      }
    }
    else{
      mag=0;
      tf=0.245; // S&N92
      gam=(7.0/5.0);
    }

    x1in=0;
    x1out=1;
    x2in=0;
    x2out=1;
    x3in=0;
    x3out=1;

    // assuming all physics is on by default:
    visc_real=0;
    res_real=0;
    mdotin=0;
    cool=0;
    nonunigridx1=0;
    nonunigridx2=0;
    nonunigridx3=0;
    simplebc=1;
    if(MAGDIR==1){
      bcix1=4;
      bcox1=4;
      bcix2=5;
      bcox2=5;
      bcix3=5;
      bcox3=5;
    }
    else if(MAGDIR==2){
      bcix1=5;
      bcox1=5;
      bcix2=4;
      bcox2=4;
      bcix3=5;
      bcox3=5;
    }
    else if(MAGDIR==3){
      bcix1=5;
      bcox1=5;
      bcix2=5;
      bcox2=5;
      bcix3=4;
      bcox3=4;
    }

    DTl=tf/100.0;
    DTfloor=DTd=tf/10.0;
    DTpd=DTd*2.0;
    DTi=10000.0;
    DTener=DTloss=DTl;
    DTtimestep=DTl*10;
    DTsp=10000.0;

  }
  else if(calltype==0){

    // want to write some interesting data on solution
    if(firstsolve==1){
    
      strcpy(temps,DATADIR);

      sprintf(filename,"%s0_analdata%s%s",temps,DAT2EXT,myidtxt);
      if((analyticout=fopen(filename,"wt"))==NULL){
	fprintf(fail_file,"Cannot open %s\n",temps);
	myexit(1);
      }
    }

    if((!MAGSOD)||(mag==0)){

      R0=L[1][1]+L[2][1];
      mu=sqrt((gam-1.)/(gam+1.));
      // initial conditions
      // override gam
      Pl=1;
      rl=1;
      Vl=0;
      Pr=.1;
      rr=.125 ;
      Vr=0;

      // Find Vm and Pm


#define V1 (Vr+(P-Pr)*sqrt((1.-mu*mu)/(rr*(P+mu*mu*Pr))))
#define V2 (Vl-sqrt(1.-pow(mu,4.))/(mu*mu)*pow(Pl,1./(2.*gam))/sqrt(rl)*(pow(P,(gam-1.)/(2.*gam))-pow(Pl,(gam-1.)/(2.*gam))))

      // Assume Vr=Vl=0 from now on

      Vm=.928;
      Pm=.303; // for gam=1.4

      rho2=rr*( (gam+1.)*Pm/(gam-1.)/Pr+1.)/( (gam+1.)/(gam-1.)+Pm/Pr);

      Vs=rho2/(rho2-rr)*Vm;

      rm=rl*pow((Pm/Pl),1./gam);


      as0=sqrt(gam*Pl/rl); // not sure
      as_m=sqrt(gam*Pm/rm);
      as_r=sqrt(gam*Pr/rr);

      //printf("rho2: %lf Vs: %lf rm: %lf as0: %lf as_m %lf\n",rho2,Vs,rm,as0,as_m);


#define VR(x) ( (1.-mu*mu)*( (x-0.5)/(t+SSMALL)+as0) )
#define asR(x) ( mu*mu*(x-0.5)/(t+SSMALL)+(1.-mu*mu)*as0 )

#define RHOR(x) (rm*pow(1.-(gam-1.)/2.*(VR(x)-Vm)/as_m,2./(gam-1.)))
#define PR(x) (Pm*pow(1.-(gam-1.)/2.*(VR(x)-Vm)/as_m,2.*gam/(gam-1.)))


      // Find points at which regions are seperated

      // Vl->Vrar

  
      for(i=INFULL1;i<OUTFULL1;i++){
	//printf("t: %lf Vl: %lf VR: %lf\n",t,Vl,VR(x[1][1][i]));
	if(Vl<VR(x[2][1][i])){
	  regions[0]=x[2][1][i];
	  break;
	}
      }
      if(i==OUTFULL1){
	regions[0]=x[2][1][INFULL1]-.01; // fudge it back
      }
  
      //regions[0] = -as0*t + 0.5 ;

      for(i=INFULL1;i<OUTFULL1;i++){
	if(VR(x[2][1][i])>Vm){
	  regions[1]=x[2][1][i];
	  break;
	}
      }
      if(i==OUTFULL1){
	regions[1]=x[2][1][INFULL1]-.01; // fudge it back
      }
      //regions[1]=0.5;
      regions[2]=Vm*t+0.5; // assuming contact discontinuity started at x=0 and moves to right
      regions[3]=Vs*t+0.5;

      //  printf("detalxreal: %lf\n",regions[1]-regions[0]);

      //  if(t>.240) printf("0: %lf 1: %lf 2: %lf 3: %lf\n",regions[0],regions[1],regions[2],regions[3]);
      //printf("0: %lf 1: %lf 2: %lf 3: %lf\n",regions[0],regions[1],regions[2],regions[3]);

      LOOPF{
	pos=x[2][1][i];
	ftemp=VR(pos);
	//    if(pos<=regions[0]){
	if(ftemp<Vl){
	  sanal[1][k][j][i]=rl;
	  //xsanal[1][k][j][i]=Pl;
	  //xsanal[2][k][j][i]=as0; // ?
	  sanal[2][k][j][i]=Pl/(gam-1.);
	}
	//    if( (pos<=regions[1])&& (pos>regions[0]))  {

	if( (ftemp<Vm)&&(ftemp>Vl))   {
	  sanal[1][k][j][i]=RHOR(pos);
	  //xsanal[1][k][j][i]=PR(pos);
	  //xsanal[2][k][j][i]=asR(pos); // ?
	  sanal[2][k][j][i]=PR(pos)/(gam-1.);
	}
	//    if( (pos<=regions[2])&& (pos>regions[1]))  {
	if( (RHOR(pos)<rm)&&(pos<=regions[2]) ){
	  sanal[1][k][j][i]=rm;
	  //xsanal[1][k][j][i]=Pm;
	  //xsanal[2][k][j][i]=as_m; // ?
	  sanal[2][k][j][i]=Pm/(gam-1.);
	}
	if( (pos<=regions[3])&& (pos>regions[2]))  {
	  sanal[1][k][j][i]=rho2;
	  //xsanal[1][k][j][i]=Pm;
	  //xsanal[2][k][j][i]=as_m; // ?
	  sanal[2][k][j][i]=Pm/(gam-1.);
	}
	if(pos>regions[3])  {
	  sanal[1][k][j][i]=rr;
	  //xsanal[1][k][j][i]=Pr;
	  //xsanal[2][k][j][i]=as_r;
	  sanal[2][k][j][i]=Pr/(gam-1.);
	}
      }
      LOOPF{
	hit=0;
	pos=x[1][1][i];
	ftemp=VR(pos);

	//    if( ((pos<=regions[1])||(ftemp<Vm))&&( (pos>regions[0])||(ftemp>Vl)))  {
	//if( (pos<=regions[1])&&( pos>regions[0]))  {// up small spike
	if( (ftemp<Vm)&&(ftemp>Vl)){ // down large spike
	  vanal[1][1][k][j][i]=VR(pos);
	  hit=1;
	}
	// correct if miss above and get a large down spike
	if( (pos<=regions[1])&&( pos>regions[0]))  {// up small spike
	  if(hit==0){
	    vanal[1][1][k][j][i]=Vm;
	  }
	}
	if(pos<=regions[0]){
	  //if(ftemp<Vl){
	  vanal[1][1][k][j][i]=Vl;
	}
	if( (pos<=regions[2])&& (pos>regions[1]))  {
	  vanal[1][1][k][j][i]=Vm;
	}
	if( (pos<=regions[3])&& (pos>regions[2]))  {
	  vanal[1][1][k][j][i]=Vm;
	}
	if(pos>regions[3])  {
	  vanal[1][1][k][j][i]=Vr;
	}
      }

      // mag field
      LOOPF{
	vanal[2][1][k][j][i]=0;
	vanal[2][2][k][j][i]=0;
	vanal[2][3][k][j][i]=0;
      }

      // convolve the data

      if(CONVOLVESOD){
	// get smoother
	ftemp=0;
	for(i=0;i<=2*WIDTH;i++){
	  SMOOTHER[i]=cos(2.0*M_PI/(4.0*(SFTYPE)WIDTH)*((SFTYPE)i-(SFTYPE)WIDTH));
	  ftemp+=SMOOTHER[i]*dx[1][1][0]; // only right for uni grid
	  fprintf(stderr,"%d %15.10g\n",i,SMOOTHER[i]); 
	}
	// normalize
	for(i=0;i<=2*WIDTH;i++){
	  SMOOTHER[i]/=ftemp;
	}
	LOOPF{
	  work1[k][j][i]=0;
	  work2[k][j][i]=0;
	  work3[k][j][i]=0;
	  for(l=i-WIDTH;l<=i+WIDTH;l++){
	    if( (l>=INFULL1)&&(l<OUTFULL1) ){
	      work1[k][j][i]+=vanal[1][1][k][j][l]*SMOOTHER[WIDTH+i-l]*dx[2][1][l];
	      work2[k][j][i]+=sanal[1][k][j][l]*SMOOTHER[WIDTH+i-l]*dx[1][1][l];
	      work3[k][j][i]+=sanal[2][k][j][l]*SMOOTHER[WIDTH+i-l]*dx[1][1][l];
	    }
	  }
	}
	LOOPF{
	  vanal[1][1][k][j][i]=work1[k][j][i];
	  sanal[1][k][j][i]=work2[k][j][i];
	  sanal[2][k][j][i]=work3[k][j][i];
	}
      }

    }
    else{ // Magnetic sod shock problem

      if(MAGSODTYPE==1){
	Pl=1.0;
	rl=1.0;
	Vxl=0.0;
	Vyl=0.0;
	Vzl=0.0;
	Bxl=.75;
	Byl=1.0;
	Bzl=0.0;
      
	Pr=.1;
	rr=.125 ;
	Vxr=0.0;
	Vyr=0.0;
	Vzr=0.0;
	Bxr=0.75;
	Byr=-1.0;
	Bzr=0.0;
	Vr=0;
      }
      else if(MAGSODTYPE==2){
	Pl=20.0;
	rl=1.0;
	Vxl=10.0;
	Vyl=0.0;
	Vzl=0.0;
	Bxl=5.0/sqrt(4.0*M_PI);
	Byl=5.0/sqrt(4.0*M_PI);
	Bzl=0.0;
      
	Pr=0.00352107;
	rr=1.0 ;
	Vxr=-10.0;
	Vyr=0.0;
	Vzr=0.0;
	Bxr=5.0/sqrt(4.0*M_PI);
	Byr=5.0/sqrt(4.0*M_PI);
	Bzr=0.0;
      }
      else if(MAGSODTYPE==3){
	Pl=1.0;
	rl=1.0;
	Vxl=0.0;
	Vyl=0.0;
	Vzl=0.0;
	Bxl=3.0/sqrt(4.0*M_PI);
	Byl=5.0/sqrt(4.0*M_PI);
	Bzl=0.0;
      
	Pr=10.0;
	rr=0.1 ;
	Vxr=0.0;
	Vyr=0.0;
	Vzr=0.0;
	Bxr=3.0/sqrt(4.0*M_PI);
	Byr=2.0/sqrt(4.0*M_PI);
	Bzr=0.0;
      }
      else if(MAGSODTYPE==4){
	Pl=.95;
	rl=1.08;
	Vxl=1.2;
	Vyl=.01;
	Vzl=0.5;
	Bxl=2.0/sqrt(4.0*M_PI);
	Byl=3.6/sqrt(4.0*M_PI);
	Bzl=2.0/sqrt(4.0*M_PI);
      
	Pr=1.0;
	rr=1.0 ;
	Vxr=0.0;
	Vyr=0.0;
	Vzr=0.0;
	Bxr=2.0/sqrt(4.0*M_PI);
	Byr=4.0/sqrt(4.0*M_PI);
	Bzr=2.0/sqrt(4.0*M_PI);
      }
      else if(MAGSODTYPE==5){
	Pl=1.0;
	rl=1.0;
	Vxl=0.0;
	Vyl=0.0;
	Vzl=0.0;
	Bxl=3.0/sqrt(4.0*M_PI);
	Byl=6.0/sqrt(4.0*M_PI);
	Bzl=0.0;
      
	Pr=10.0;
	rr=0.1 ;
	Vxr=0.0;
	Vyr=2.0;
	Vzr=1.0;
	Bxr=3.0/sqrt(4.0*M_PI);
	Byr=1.0/sqrt(4.0*M_PI);
	Bzr=0.0;
      }
      else if(MAGSODTYPE==6){
	Pl=0.4;
	rl=0.1;
	Vxl=50.0;
	Vyl=0.0;
	Vzl=0.0;
	Bxl=0.0;
	Byl=-1.0/sqrt(4.0*M_PI);
	Bzl=-2.0/sqrt(4.0*M_PI);
      
	Pr=0.2;
	rr=0.1 ;
	Vxr=0.0;
	Vyr=0.0;
	Vzr=0.0;
	Bxr=0.0;
	Byr=1.0/sqrt(4.0*M_PI);
	Bzr=2.0/sqrt(4.0*M_PI);
      }
      else if(MAGSODTYPE==7){
	Pl=1.0;
	rl=1.0;
	Vxl=-1.0;
	Vyl=0.0;
	Vzl=0.0;
	Bxl=0.0;
	Byl=1.0;
	Bzl=0.0;
      
	Pr=1.0;
	rr=1.0 ;
	Vxr=1.0;
	Vyr=0.0;
	Vzr=0.0;
	Bxr=0.0;
	Byr=1.0;
	Bzr=0.0;
      }
      else if(MAGSODTYPE==8){
	Pl=1.0;
	rl=1.0;
	Vxl=0.0;
	Vyl=0.0;
	Vzl=0.0;
	Bxl=1.0;
	Byl=1.0;
	Bzl=0.0;
      
	Pr=0.1;
	rr=0.2;
	Vxr=0.0;
	Vyr=0.0;
	Vzr=0.0;
	Bxr=1.0;
	Byr=0.0;
	Bzr=0.0;
      }
      else if(MAGSODTYPE==9){
	Pl=.52467;
	rl=0.4;
	Vxl=-.66991;
	Vyl=.98263;
	Vzl=0.0;
	Bxl=1.3;
	Byl=.0025293;
	Bzl=0.0;
      
	Pr=1.0;
	rr=1.0;
	Vxr=0.0;
	Vyr=0.0;
	Vzr=0.0;
	Bxr=1.3;
	Byr=1.0;
	Bzr=0.0;
      }
      else if(MAGSODTYPE==10){
	Pl=.5;
	rl=0.65;
	Vxl=.667;
	Vyl=-.257;
	Vzl=0.0;
	Bxl=.75;
	Byl=.55;
	Bzl=0.0;
      
	Pr=.75;
	rr=1.0;
	Vxr=0.4;
	Vyr=-.94;
	Vzr=0.0;
	Bxr=.75;
	Byr=0.0;
	Bzr=0.0;
      }
      else if(MAGSODTYPE==11){
	Pl=1.0;
	rl=1.0;
	Vxl=0.0;
	Vyl=0.0;
	Vzl=0.0;
	Bxl=0.7;
	Byl=0.0;
	Bzl=0.0;
      
	Pr=.2;
	rr=0.3;
	Vxr=0.0;
	Vyr=0.0;
	Vzr=1.0;
	Bxr=0.7;
	Byr=1.0;
	Bzr=0.0;
      }
      else if(MAGSODTYPE==12){
	Pl=1.0;
	rl=1.0;
	Vxl=0.0;
	Vyl=0.0;
	Vzl=0.0;
	Bxl=0.75;
	Byl=1.0;
	Bzl=0.0;
      
	Pr=.1;
	rr=0.125;
	Vxr=0.0;
	Vyr=0.0;
	Vzr=0.0;
	Bxr=0.75;
	Byr=-1.0;
	Bzr=0.0;
      }
      else if(MAGSODTYPE==13){
	Pl=1.0;
	rl=1.0;
	Vxl=0.0;
	Vyl=0.0;
	Vzl=0.0;
	Bxl=1.3;
	Byl=1.0;
	Bzl=0.0;
      
	Pr=.4;
	rr=0.4;
	Vxr=0.0;
	Vyr=0.0;
	Vzr=0.0;
	Bxr=1.3;
	Byr=-1.0;
	Bzr=0.0;
      }

      LOOPF{

	if(MAGDIR==1){
	  if(x[2][1][i]<x1in+x1out*0.5){ // first half
	    sanal[1][k][j][i]=rl;
	    sanal[2][k][j][i]=Pl/(gam-1.0);
	    sanal[3][k][j][i]=0.0;
	  
	    vanal[1][1][k][j][i]=Vxl;
	    vanal[1][2][k][j][i]=Vyl;
	    vanal[1][3][k][j][i]=Vzl;
	  
	    vanal[2][1][k][j][i]=Bxl;
	    vanal[2][2][k][j][i]=Byl;
	    vanal[2][3][k][j][i]=Bzl;
	  }
	  else{ // right half
	    sanal[1][k][j][i]=rr;
	    sanal[2][k][j][i]=Pr/(gam-1.0);
	    sanal[3][k][j][i]=0.0;
	  
	    vanal[1][1][k][j][i]=Vxr;
	    vanal[1][2][k][j][i]=Vyr;
	    vanal[1][3][k][j][i]=Vzr;
	  
	    vanal[2][1][k][j][i]=Bxr;
	    vanal[2][2][k][j][i]=Byr;
	    vanal[2][3][k][j][i]=Bzr;
	  }
	}
	else if(MAGDIR==2){
	  if(x[2][2][j]<x2in+x2out*0.5){ // first half
	    sanal[1][k][j][i]=rl;
	    sanal[2][k][j][i]=Pl/(gam-1.0);
	    sanal[3][k][j][i]=0.0;
	  
	    vanal[1][1][k][j][i]=Vzl;
	    vanal[1][2][k][j][i]=Vxl;
	    vanal[1][3][k][j][i]=Vyl;
	  
	    vanal[2][1][k][j][i]=Bzl;
	    vanal[2][2][k][j][i]=Bxl;
	    vanal[2][3][k][j][i]=Byl;
	  }
	  else{ // right half
	    sanal[1][k][j][i]=rr;
	    sanal[2][k][j][i]=Pr/(gam-1.0);
	    sanal[3][k][j][i]=0.0;
	  
	    vanal[1][1][k][j][i]=Vzr;
	    vanal[1][2][k][j][i]=Vxr;
	    vanal[1][3][k][j][i]=Vyr;
	  
	    vanal[2][1][k][j][i]=Bzr;
	    vanal[2][2][k][j][i]=Bxr;
	    vanal[2][3][k][j][i]=Byr;
	  }
	}
	else if(MAGDIR==3){
	  if(x[2][3][k]<x3in+x3out*0.5){ // first half
	    sanal[1][k][j][i]=rl;
	    sanal[2][k][j][i]=Pl/(gam-1.0);
	    sanal[3][k][j][i]=0.0;
	  
	    vanal[1][1][k][j][i]=Vyl;
	    vanal[1][2][k][j][i]=Vzl;
	    vanal[1][3][k][j][i]=Vxl;
	  
	    vanal[2][1][k][j][i]=Byl;
	    vanal[2][2][k][j][i]=Bzl;
	    vanal[2][3][k][j][i]=Bxl;
	  }
	  else{ // right half
	    sanal[1][k][j][i]=rr;
	    sanal[2][k][j][i]=Pr/(gam-1.0);
	    sanal[3][k][j][i]=0.0;
	  
	    vanal[1][1][k][j][i]=Vyr;
	    vanal[1][2][k][j][i]=Vzr;
	    vanal[1][3][k][j][i]=Vxr;
	  
	    vanal[2][1][k][j][i]=Byr;
	    vanal[2][2][k][j][i]=Bzr;
	    vanal[2][3][k][j][i]=Bxr;
	  }
	}
      }

    }
    // output some interesting data
    // and setup initial problem stuff
    if(firstsolve==1){
      firstsolve=0;
  
      // now output some interesting analytic data

      fprintf(analyticout,"#hello!\n");
      LOOPF{
	fprintf(analyticout,"%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g\n",sanal[1][k][j][i],sanal[2][k][j][i],sanal[3][k][j][i],vanal[1][1][k][j][i],vanal[1][2][k][j][i],vanal[1][3][k][j][i],vanal[2][1][k][j][i],vanal[2][2][k][j][i],vanal[2][3][k][j][i]);
      }
      fclose(analyticout);
    }

  }// end else if not just initial problem setup
}

// Advection test dump at specific time
// Fill <V/r/P/e>anal[i] with solution to advection solution for given coord system
// must use dirichlet or periodic BCs in general


void advsol(int calltype)
{ // solution is on b-grid as normal
  // for mass and vx1 advection only

  int i,j,k;
  SFTYPE temp,temp2;
  SFTYPE exponent=0;
  SFTYPE RHOI,RHO0,VI,V0,XI,X0,Const;
  int CASE;
  int SN;

  CASE=0;
  // 0: dv/dt=0
  // 1: dv/dt =/=0

  if(CASE==1){  // CASE1 context

    RHOI=1.0; // rho(t=0,x=0)=RHOI
    V0=0.1;
    VI=0.01; // v(t=0,x=0)=VI
    X0=1.0;  // v(t=0,x=x0)=V0
  }
  else if(CASE==0){
    
    RHOI=0.9;  // rho(t=0,x=0)=RHOI
    RHO0=1.0; // rho(t=0,x=x0)=RHO0
    V0=0.001;
    X0=1.0;  // v(t=0,x=x0)=V0
    
    Const=COORD*1.0; // S&N condition
    
    SN=1;
  }

  R0=L[1][1]+L[2][1];

  LOOPF{
    if(CASE==1){
      if(transv1x1==0){
	fprintf(fail_file,"Case=1 requires velocity transport\n");
	myexit(1);
      }
      if(COORD==1){
	vanal[1][1][k][j][i]=(x[1][1][i]+VI*X0/(V0-VI))/(t+X0/(V0-VI));
	sanal[1][k][j][i]=RHOI*(V0-VI)/(X0*sqrt(VI))*pow( (t+X0/(V0-VI))*(x[2][1][i]+VI*X0/(V0-VI)),-0.5);
      }
      if((COORD==3)||(COORD==2)){
	fprintf(fail_file,"Case=1 has no coord==2or3 case\n");
	myexit(1);	
      }
    }
    else if(CASE==0){
      if(transv1x1==1){
	fprintf(fail_file,"Case=0 requires no velocity transport\n");
	myexit(1);
      }
      if(COORD==1){
	vanal[1][1][k][j][i]=V0;
	sanal[1][k][j][i]=RHOI*exp(V0/X0*log(RHO0/RHOI)*(x[2][1][i]/V0-t));
      }
      if(COORD==2){
	// r is x2 here
	vanal[1][2][k][j][i]=V0*x[1][2][j]; // force this condition to get density
	sanal[1][k][j][i]=RHO0*exp(-Const*t)*pow(x[2][2][j],(Const-2.0*V0)/V0); // if Const=2.0, as in SN
      }
      if(COORD==3){
	vanal[1][1][k][j][i]=V0*x[1][1][i]; // force this condition to get density
	if(SN==1){
	  sanal[1][k][j][i]=RHO0*exp(-Const*V0*t);
	}
	else{
	  sanal[1][k][j][i]=RHO0*exp(-Const*V0*t)*pow(x[2][1][i],(Const-3.0)); // if Const=3.0, as in SN
	}
      }
    }

    sanal[2][k][j][i]=0.0;
    sanal[3][k][j][i]=0.0;
    if(COORD!=2){
      vanal[1][2][k][j][i]=0.0;
    }
    else vanal[1][1][k][j][i]=0.0;
    vanal[1][3][k][j][i]=0.0;
  }
}


// Advection test dump at specific time
// Fill <s/v>anal[k][j][i] with solution to simple gauss adv problem

void gausssol(int calltype)
{ // solution is on b-grid as normal

  int i,j,k;
  SFTYPE temp,temp2;
  SFTYPE exponent=0;
  SFTYPE A,q10;

#define X0(i) (q10+vanal[1][1][k][j][i]*t)
#define OVERLAP(q1,i) floor(q1-(X0(i)-0.5))
#define F(q1,i) (1.0+0.1*exp(-pow( (q1-X0(i)-OVERLAP(q1,i))/A,2.0)))

  A=(1.0/12.0);
  q10= (1.0/6.0);

  /*
    #define X01(i) (q110+vanal[1][1][k][j][i]*t)
    //#define X02(j) (q120+vanal[1][2][k][j][i]*t)
    #define OVERLAP1(q1,i) floor(q1-(X01(i)-0.5))
    //#define OVERLAP2(q2,j) floor(q2-(X02(j)-0.5))
    #define F(q1,q2,i,j) (1.0+0.1*exp( -pow( (q1-X01(i)-OVERLAP1(q1,i)),2.0)/(A*A)))
    //#define F(q1,q2,i,j) (1.0+0.1*exp( (-pow( (q1-X01(i)-OVERLAP1(q1,i)),2.0)-pow( (q2-X02(j)-OVERLAP2(q2,j)),2.0))/(A*A)))
    */
  R0=L[1][1]+L[2][1];
  LOOPF{
    vanal[1][1][k][j][i]=0.1; // on a-grid, but same on b-grid since constant
    vanal[1][2][k][j][i]=0.1; 
    vanal[1][3][k][j][i]=0.0; 

    //    sanal[1][k][j][i]=F(x[2][1][i],x[2][2][j],i,j);
    sanal[1][k][j][i]=F(x[2][1][i],i);
    sanal[2][k][j][i]=0.0;
    sanal[3][k][j][i]=0.0;
  }

  for(i=0;i<ITYPES;i++){
    for(j=0;j<CTYPES;j++){
      mms[i][j][1][0]=0.5;
      mms[i][j][1][1]=1.5;
      mms[i][j][2][0]=0;
      mms[i][j][2][1]=1.0;
      mms[i][j][3][0]=0.0;
      mms[i][j][3][1]=1.0;
	
      mmv[i][j][1][0][0]=0;
      mmv[i][j][1][0][1]=1.0;
      mmv[i][j][1][1][0]=0;
      mmv[i][j][1][1][1]=1.0;
      mmv[i][j][1][2][0]=0;
      mmv[i][j][1][2][1]=1.0;
      mmv[i][j][1][3][0]=0;
      mmv[i][j][1][3][1]=1.0;

      // define outer region when interpolation is used.
      // same order as scalar/vector arrays
      outerdefs[i][j][1]=mms[0][0][1][0]; // rho
      outerdefs[i][j][2]=mms[0][0][2][0]; // en
      outerdefs[i][j][3]=mms[0][0][3][0]; // pot
	
      outerdefv[i][j][1][0]=0.0; // magnitude of v
      outerdefv[i][j][1][1]=0.0; // v1
      outerdefv[i][j][1][2]=0.0; // v2
      outerdefv[i][j][1][3]=0.0; // v3
	
      outerdefv[i][j][2][0]=0.0;
      outerdefv[i][j][2][1]=0.0;
      outerdefv[i][j][2][2]=0.0;
      outerdefv[i][j][2][3]=0.0;

    }
  }
    
  douterdefs[1]=mms[0][0][1][0]; // rho
  douterdefs[2]=mms[0][0][2][0]; // en
  douterdefs[3]=mms[0][0][3][0]; // pot
    
  douterdefv[1][0]=0.0; // magnitude of v
  douterdefv[1][1]=0.0; // v1
  douterdefv[1][2]=0.0; // v2
  douterdefv[1][3]=0.0; // v3
    
  douterdefv[2][0]=0.0;
  douterdefv[2][1]=0.0;
  douterdefv[2][2]=0.0;
  douterdefv[2][3]=0.0;

}

// boundary test, using different pulses at boundaries

// also advection test if source/transv off


void pulsesol(int calltype)
{ // solution is on b-grid as normal

  int i,j,k;
  SFTYPE posx,posy,posz;
  SFTYPE fracx,fracy,fracz;
  SFTYPE ftemp1,ftemp2;
  int PULSENUM, MAGPULSE,VELNUM;
  SFTYPE size1,size2,size3,pos1,pos2,pos3;
  FTYPE (*vectorpot)[N3M][N2M][N1M]; // phi component of vector potential
  SFTYPE magbeta;
  SFTYPE b0,torib0,vertb0,pgtot,pbtot,pgtot_full,pbtot_full;
  SFTYPE rho0,e0;
  int gooddensity;
  int GLOBALFIELD;
  SFTYPE divb1avg,divb2avg;
  SFTYPE divb1max,divb2max;
  SFTYPE divb1avg_full,divb2avg_full;
  SFTYPE divb1max_full,divb2max_full;
  SFTYPE ftemp;
  SFTYPE drho,wavelength;

  


  vectorpot=workv1;

  if(calltype==100){
    
    if(POSTPROC==1){
      DYNAMICMM=0; // already read in file usually
    }
    else{
      DYNAMICMM=2;
    }
    if(COORD==1){
      /*
	x1in=0;
	x1out=1;
	x2in=0;
	x2out=1;
	x3in=0;
	x3out=1;
      */
      x1in=-23;
      x1out=23;
      x2in=-23;
      x2out=23;
      x3in=-23;
      x3out=23;

    }
    
    if(COORD==3){
      x1in=1;
      x1out=20;
      x2in=0;
      x2out=M_PI;
      //x2in=0.00001;
      //x2out=M_PI-.00001;
      x3in=0;
      x3out=2.0*M_PI;
    }

    //    tf=.94868329805;

    //tf=80.0;
    tf=120.0;
    //tf=1000;


    // assuming all physics is on by default:
    visc_real=0;

    // MARK
    //    res_real=1; rreal=2; resheat= 1; // resistivity
    res_real=0;

    mdotin=0;
    cool=0;
    nonunigridx1=0;
    nonunigridx2=0;
    nonunigridx3=0;
    simplebc=1;

    
    bcix1=5;
    bcox1=5;
    bcix2=5;
    bcox2=5;
    bcix3=5;
    bcox3=5;
    
    /*
      bcix1=4;
      bcox1=4;
      bcix2=4;
      bcox2=4;
      //bcix2=5;
      //bcox2=5;
      bcix3=4;
      bcox3=4;
    */
    pressx1=1;
    pressx2=1;
    pressx3=1;
    mag=0;
    ie=1;
    visc_art=1;
    trans=1;
    transx2=1;
    transx1=1;
    transx3=1;

    //transv3x1=0;
    //transmagx1=0;
    //transv3x2=0;
    //transmagx2=0;
    //transbzx=transbzy=0;

    /*    
	  DTl=DTener=DTloss=tf/5000.0;
	  DTtimestep=DTsp=tf/100.0;  DTfloor=DTd=tf/(2.5806452*10.0); DTpd=tf/(2.5806452*5.0); DTi=tf/(2.5806452*300.0);
	  DTfld=DTi*2.0;

	  DTd=1.05409/50.0;
    */
    DTl=DTener=DTloss=tf/1000.0;
    DTtimestep=DTsp=tf/500.0;  DTfloor=DTd=tf/(20.0); DTpd=(tf/20.0); DTi=tf/(300.0);
    DTfld=DTd;



  }
  else if(calltype==0){

    gam=5.0/3.0;
    //  magbeta=100.0; // controls absolute strength of magnetic field
    magbeta=100.0; // controls absolute strength of magnetic field
    torib0=1.0; // controls relative strength to any other fields added
    vertb0=torib0/10.0; // controls relative strength to any other fields added
    rho0=1.0;
    drho=.001;
    wavelength=0.25;
    e0=1.0; // K really
    GLOBALFIELD=0;

    DENSITYFLOOR=1E-5;
    IEFRACFLOOR=DENSITYFRACFLOOR=1E-5;

    MAGPULSE=0;
    // global params
    PULSENUM= 1;
    // 0: point shock
    // 1: square
    // 2: gaussian
    // 3: spherical
    // 4: just sound wave(see waves() for alfven wave)
    VELNUM= 0;
    
    
    R0=L[1][1]+L[2][1];
    
    fracx=.1;
    fracy=.1;
    fracz=.1;
    
    LOOPF{
      posx=x[2][1][i];
      posy=x[2][2][j];
      posz=x[2][3][k];

      if(PULSENUM==0){
	// 3D tests
	if( ((k==N3/2)||(k==N3/2-1))&&((j==N2/2)||(j==N2/2-1))&&((i==N1/2)||(i==N1/2-1)) ){ // symmetric for even grid
	  //	if((k==N3/8)&&(j==N2/8)&&(i==N1/8)){
	  // 1D tests
	  //if((i==N1/2)||(i==N1/2-1)){
	  //	if((j==N2/2)||(j==N2/2-1)){
	  //if((k==N3/2)||(k==N3/2-1)){
	  // 2D tests x-y
	  //	if(((i==N1/2)||(i==N1/2-1))&&((j==N2/2)||(j==N2/2-1))){
	  sanal[1][k][j][i]=1.0;
        }
	else sanal[1][k][j][i]=0.1;
	//else sanal[1][k][j][i]=1.0+0.1*sin(x[1][1][i]);
	sanal[2][k][j][i]=pow(sanal[1][k][j][i],gam);
	sanal[3][k][j][i]=0.0;
	vanal[1][1][k][j][i]=vanal[1][2][k][j][i]=vanal[1][3][k][j][i]=0.0;
	vanal[2][1][k][j][i]=vanal[2][2][k][j][i]=vanal[2][3][k][j][i]=0.0;
      }
      if(PULSENUM==4){
	// initialize background
	sanal[1][k][j][i]=rho0;
	sanal[2][k][j][i]=e0;
	cs=gam*(gam-1.)*sanal[2][k][j][i]/sanal[1][k][j][i];
	sanal[1][k][j][i]+=drho*sin(2.0*M_PI/wavelength*x[2][1][i]);
	sanal[2][k][j][i]+=drho*cs*cs/(gam-1.0)*sin(2.0*M_PI/wavelength*x[2][1][i]); // p=(gam-1)u dP=cs^2 d\rho
	vanal[1][1][k][j][i]=drho*cs*cs/rho0*sin(2.0*M_PI/wavelength*x[2][1][i]); // cs^2/\rho0*drho
	vanal[1][2][k][j][i]=0;
	vanal[1][3][k][j][i]=0;

	vanal[2][1][k][j][i]=0;
	vanal[2][2][k][j][i]=0;
	vanal[2][3][k][j][i]=0;
      }

      if((PULSENUM==1)||(PULSENUM==3)){
	
	pos1=L[1][1]+L[2][1]*.1;
	size1=L[2][1]*fracx*.1; // the 1/2 size
	pos2=L[1][2]+L[2][2]*.1;
	size2=L[2][2]*fracy*.1;
	pos3=L[1][3]+L[2][3]*.1;
	size3=L[2][3]*fracz*.1;

	if(
	   ((PULSENUM==1)&&( 
			    ( (N1==1)||( (posx>pos1-size1)&& (posx<pos1+size1) ))&&
			    ( (N2==1)||( (posy>pos2-size2)&& (posy<pos2+size2) ))&&
			    ( (N3==1)||( (posz>pos3-size3)&& (posz<pos3+size3) ))
			    ))
	   ||
	   ((PULSENUM==3)&&( 
			    ((posx-pos1)*(posx-pos1)+(posy-pos2)*(posy-pos2)+(posz-pos3)*(posz-pos3)<size1*size1)
			    ))
	   ) {
	  sanal[1][k][j][i]=rho0;
	  sanal[2][k][j][i]=e0*pow(sanal[1][k][j][i],gam);
	  if(VELNUM==2) vanal[1][1][k][j][i]=0.1;
	}
	else{
	  //sanal[1][k][j][i]=rho0*0.999;
	  sanal[1][k][j][i]=rho0*0.9;
	  //sanal[1][k][j][i]=1.0;
	  sanal[2][k][j][i]=e0*pow(sanal[1][k][j][i],gam);
	  if(VELNUM==2) vanal[1][1][k][j][i]=0.0;
	}
	
      }
      if(PULSENUM==2){
	ftemp1=x[2][1][i]-(L[1][1]+L[2][1]*0.5);
	ftemp2=x[2][2][j]-(L[1][2]+L[2][2]*0.5);
	
	sanal[1][k][j][i]=1.1*exp(-(ftemp1*ftemp1+ftemp2*ftemp2)/(0.1));
	//sanal[2][k][j][i]=0.001*pow(sanal[1][k][j][i],gam);
	sanal[2][k][j][i]=0.0;
	if(VELNUM==2) vanal[1][1][k][j][i]=0.1;
      }

      if((PULSENUM>0)&&(PULSENUM<4)){
	if(VELNUM==0){
	  vanal[1][1][k][j][i]=0.0; // on a-grid, but same on b-grid since constant
	}
	else if(VELNUM==1){
	  if(x[1][1][i]<L[1][1]+L[2][1]*0.5){
	    vanal[1][1][k][j][i]=1.0/(L[2][1])*(x[1][1][i]-L[1][1]);
	  }
	  else if(x[1][1][i]>=L[1][1]+L[2][1]*0.5){
	    vanal[1][1][k][j][i]=-1.0/(L[2][1])*(x[1][1][i]-(L[1][1]+L[2][1]));
	  }
	}
	vanal[1][2][k][j][i]=0.0; // on a-grid, but same on b-grid since constant
	vanal[1][3][k][j][i]=0.0; // on a-grid, but same on b-grid since constant
	
	sanal[3][k][j][i]=0.0;
      }
    
    } // end loopf


    if((mag==1)&&(MAGPULSE)){// can't do in loop above since use z2c
      LOOPF{// must initialize since only do LOOPH below
	vectorpot[1][k][j][i]=vectorpot[2][k][j][i]=vectorpot[3][k][j][i]=0;
      }
      // setup x3 component of vector potential
      LOOPH{
	// check if density is definitely within torus so alfven velocity isn't crazy at torus edge
	if((sanal[1][k][j][i]>DENSITYFLOOR*10.0)&&
	   (sanal[1][k][j][i+1]>DENSITYFLOOR*10.0)&&
	   (sanal[1][k][j][i-1]>DENSITYFLOOR*10.0)&&
	   (sanal[1][k][j+1][i]>DENSITYFLOOR*10.0)&&
	   (sanal[1][k][j+1][i+1]>DENSITYFLOOR*10.0)&&
	   (sanal[1][k][j+1][i-1]>DENSITYFLOOR*10.0)&&
	   (sanal[1][k][j-1][i]>DENSITYFLOOR*10.0)&&
	   (sanal[1][k][j-1][i+1]>DENSITYFLOOR*10.0)&&
	   (sanal[1][k][j-1][i-1]>DENSITYFLOOR*10.0)){
	  gooddensity=1;
	}
	else gooddensity=0;
	// vector pot should be located as vector, below curl will interp as needed
	if(GLOBALFIELD||gooddensity){
	  // normal rotation
	  vectorpot[3][k][j][i]=torib0*sanal[1][k][j][i]/rho0;
	  // opposite rotation
	  //vectorpot[3][k][j][i]=-torib0*sanal[1][k][j][i]/rho0;
	
	  // SP00 version
	  //vectorpot[3][k][j][i]=torib0*pow(sanal[1][k][j][i]/rho0,2.0);
	  // double circle opposite rotation
	  //	vectorpot[3][k][j][i]*=sin(2.0*x[2][2][j]);
	  // double circle same rotation
	  //vectorpot[3][k][j][i]*=fabs(sin(2.0*x[2][2][j]));
	}
	else vectorpot[3][k][j][i]=0.0;
	// add vertical field (- is down, + is up)
	//      vectorpot[3][k][j][i]+= 0.5*vertb0*x[2][1][i]*G4(2,j);

	vectorpot[1][k][j][i]=vectorpot[2][k][j][i]=0.0;
      }

      // initialze field
      LOOPF{
	vanal[2][1][k][j][i]=vanal[2][2][k][j][i]=vanal[2][3][k][j][i]=0;
      }
      // only use LOOPH so initial divb calc will be right, could use LOOP and calculation will be normal otherwise
      LOOPH{// can't do entire grid since use z2c and differencing
	// assumes bounded before use and not used as fixed bc(shouldn't be used that way anyways)
	// must use normal diff method even if error near r=0, at least divB=0, course then one worries about volume terms used in code!
	vanal[2][1][k][j][i]=curlv1(vectorpot,k,j,i);

	vanal[2][2][k][j][i]=curlv2(vectorpot,k,j,i);

	vanal[2][3][k][j][i]=curlv3(vectorpot,k,j,i);
      }
    }
    else{
      LOOPF{
	vanal[2][1][k][j][i]=0.0;
	vanal[2][2][k][j][i]=0.0;
	vanal[2][3][k][j][i]=0.0;
      }
    }
    // now determine normalization constant from beta=p_gas/B^2
    if(MAGPULSE&&(mag==1)){
      b0=0;
      pgtot=0.0;
      pbtot=0.0;
      LOOP{
	pgtot+=(gam-1.0)*sanal[2][k][j][i]/(OVOL1(k,j,i)*ODX(2,3,k));
	pbtot+=0.5*(e2z_1(vanal[2][1],k,j,i)*e2z_1(vanal[2][1],k,j,i)+e2z_2(vanal[2][2],k,j,i)*e2z_2(vanal[2][2],k,j,i)+vanal[2][3][k][j][i]*vanal[2][3][k][j][i])/(OVOL1(k,j,i)*ODX(2,3,k));
      }

      // need to give all cpus the normalization constant
      if(numprocs>1){
#if(USEMPI)
	MPI_Allreduce(&(pgtot), &(pgtot_full), 1, MPI_SFTYPE, MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&(pbtot), &(pbtot_full), 1, MPI_SFTYPE, MPI_SUM,MPI_COMM_WORLD);
#endif
      }
      else{
	pgtot_full=pgtot;
	pbtot_full=pbtot;
      }

      b0=sqrt(pgtot_full/(pbtot_full*magbeta));
      LOOPF{
	vanal[2][1][k][j][i]*=b0;
	vanal[2][2][k][j][i]*=b0;
	vanal[2][3][k][j][i]*=b0;
      }
  
    }

    /*
      LOOPF{
      // MPI test
      vanal[1][1][k][j][i]=ranc(0);
      vanal[1][2][k][j][i]=ranc(0);
      vanal[1][3][k][j][i]=ranc(0);
      }
    */
  }


  for(i=0;i<ITYPES;i++){
    for(j=0;j<CTYPES;j++){
      mms[i][j][1][0]=0.001;
      mms[i][j][1][1]=1.0;
      mms[i][j][2][0]=0.001;
      mms[i][j][2][1]=.002;
      mms[i][j][3][0]=0.0;
      mms[i][j][3][1]=1.0;
	
      mmv[i][j][1][0][0]=0;
      mmv[i][j][1][0][1]=1.0;
      mmv[i][j][1][1][0]=0;
      mmv[i][j][1][1][1]=1.0;
      mmv[i][j][1][2][0]=0;
      mmv[i][j][1][2][1]=1.0;
      mmv[i][j][1][3][0]=0;
      mmv[i][j][1][3][1]=1.0;

      // define outer region when interpolation is used.
      // same order as scalar/vector arrays
      outerdefs[i][j][1]=mms[0][0][1][0]; // rho
      outerdefs[i][j][2]=mms[0][0][2][0]; // en
      outerdefs[i][j][3]=mms[0][0][3][0]; // pot
	
      outerdefv[i][j][1][0]=0.0; // magnitude of v
      outerdefv[i][j][1][1]=0.0; // v1
      outerdefv[i][j][1][2]=0.0; // v2
      outerdefv[i][j][1][3]=0.0; // v3
	
      outerdefv[i][j][2][0]=0.0;
      outerdefv[i][j][2][1]=0.0;
      outerdefv[i][j][2][2]=0.0;
      outerdefv[i][j][2][3]=0.0;


    }
  }
    
  douterdefs[1]=mms[0][0][1][0]; // rho
  douterdefs[2]=mms[0][0][2][0]; // en
  douterdefs[3]=mms[0][0][3][0]; // pot
    
  douterdefv[1][0]=0.0; // magnitude of v
  douterdefv[1][1]=0.0; // v1
  douterdefv[1][2]=0.0; // v2
  douterdefv[1][3]=0.0; // v3
    
  douterdefv[2][0]=0.0;
  douterdefv[2][1]=0.0;
  douterdefv[2][2]=0.0;
  douterdefv[2][3]=0.0;

}



// assume r goes from 0 to 15
void magbreaksol(int calltype)
{ // assume COORD==2

  static FILE*analyticout;
  static int firstsolve=1;
  char temps[100];
  char filename[100];
  SFTYPE ftemp;
  SFTYPE pos;
  int i,j,k;
  SFTYPE rhod,rhox,Omega0,pd,px,B0,THICK;
  int DIC;

  DIC=1; // 0=CIC 1=DIC


  // for magtest // 3 for S&N MHD Riemman 13 for DIC 15 for CIC

  // want to write some interesting data on solution
  if(firstsolve==1){
    
    strcpy(temps,DATADIR);

    sprintf(filename,"%s0_analdata%s%s",temps,DAT2EXT,myidtxt);
    if((analyticout=fopen(filename,"wt"))==NULL){
      fprintf(fail_file,"Cannot open %s\n",temps);
      myexit(1);
    }
  }


  rhod=10.0;
  pd=10.0;
  rhox=1.0;
  px=1.0;
  Omega0=1.0;
  B0=1.0;
  THICK=1.0;

  gam=2.0;

  LOOPF{
    if(x[2][1][i]<1.0){ // the disk
      sanal[1][k][j][i]=rhod;
      sanal[2][k][j][i]=pd/(gam-1.0);
      sanal[3][k][j][i]=0.0;

      vanal[1][1][k][j][i]=0.0;
      vanal[1][2][k][j][i]=0.0;
      if(DIC){
	vanal[1][3][k][j][i]=x[2][2][j]*Omega0;
      }
      else{
	vanal[1][3][k][j][i]=x[2][2][j]*Omega0*0.5*(1.0+cos(M_PI*x[2][2][j]/THICK));
      }
      vanal[2][1][k][j][i]=B0;
      vanal[2][2][k][j][i]=0.0;
      vanal[2][3][k][j][i]=0.0;
    }
    else{ // ambient medium
      sanal[1][k][j][i]=rhox;
      sanal[2][k][j][i]=px/(gam-1.0);
      sanal[3][k][j][i]=0.0;

      vanal[1][1][k][j][i]=0.0;
      vanal[1][2][k][j][i]=0.0;
      vanal[1][3][k][j][i]=0.0;

      vanal[2][1][k][j][i]=B0;
      vanal[2][2][k][j][i]=0.0;
      vanal[2][3][k][j][i]=0.0;
    }
    
    
  }


  // output some interesting data
  // and setup initial problem stuff
  if(firstsolve==1){
    firstsolve=0;
  
    // now output some interesting analytic data

    fprintf(analyticout,"#hello!\n");
    fprintf(analyticout,"#T=%15.10g alpha=%15.10g\n",THICK*sqrt(rhox)/B0,rhod/rhox);

    LOOPF{
      fprintf(analyticout,"%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g\n",sanal[1][k][j][i],sanal[2][k][j][i],sanal[3][k][j][i],vanal[1][1][k][j][i],vanal[1][2][k][j][i],vanal[1][3][k][j][i],vanal[2][1][k][j][i],vanal[2][2][k][j][i],vanal[2][3][k][j][i]);
    }
    fclose(analyticout);
  }


}


// parameters as in Ryu Jones Frank 1995
void magvortex(int calltype)
{ // assume COORD==1

  static FILE*analyticout;
  static int firstsolve=1;
  char temps[100];
  char filename[100];
  SFTYPE ftemp;
  SFTYPE pos;
  int i,j,k;
  SFTYPE rho0,p0,en0,beta,mach,b0,v0;
  SFTYPE argsin;
  int RYUJONES,VORTEXROTATION,VORTEXPLANE;
  SFTYPE xa,xb,ya,yb;
  int dir[3+1];

  RYUJONES=0; // 1: Use Ryu Jones parameters 0: use Toth parameters
  VORTEXROTATION= 0; // 0: normal 1: rotate 90 CW 2: rotate 180 CW 3: rotate 270 CW
  VORTEXPLANE=3; //3: x-y (3-normal dir) 2: x-z (2-normal dir)  1: y-z (1-normal dir)

  if(calltype==100){

    if(POSTPROC==1){
      DYNAMICMM=0; // already read in file usually
    }
    else{
      DYNAMICMM=2; // to show all structure
    }

    if(RYUJONES){
      x1in=0;
      x1out=1;
      x2in=0;
      x2out=1;
      x3in=0;
      x3out=1;
    }
    else{
      x1in=0.0; x1out=2.0*M_PI;
      x2in=0.0; x2out=2.0*M_PI;
      x3in=0.0; x3out=2.0*M_PI;
    }

    tf=8.0;
    // vortex
    //  tf=3.1;
    //  tf=0.48*100.0/33.0; // to match gammie test of vortex
    //  tf=8.0;
    //  tf = 100.0;
    //  tf = 2.0;
    // assuming all physics is on by default:
    mag=1;
    visc_real=0;

    // MARK
    res_real=1; rreal=2; resheat= 1; // resistivity
    
    resist_real0=.01;
    // .01 : wee bit too small to kill all point shocks(although point shock is small effect), barely changes structure
    // .02 : only 1 extremely small pair of point shocks, but structure changes a bit from otherwise.
    // .1  : too high, kills lots of structure in squashed saddles in vector potential, otherwise ok.
    //res_real=0;

    mdotin=0;
    cool=0;
    nonunigridx1=0;
    nonunigridx2=0;
    simplebc=1;
    bcix1=5;
    bcox1=5;
    bcix2=5;
    bcox2=5;
    bcix3=5;
    bcox3=5;


    //transv3x1=0;
    //transmagx1=0;
    //transv3x2=0;
    //transmagx2=0;
    //transbzx=transbzy=0;
    
    DTl=DTener=DTloss=tf/5000.0;
    DTtimestep=DTsp=tf/100.0;  DTfloor=DTd=tf/(2.5806452*10.0); DTpd=tf/(2.5806452*5.0); DTi=tf/(2.5806452*300.0);
    DTfld=DTi*2.0;
  }
  else if(calltype==0){

    // want to write some interesting data on solution
    if(firstsolve==1){
    
      strcpy(temps,DATADIR);

      //      fprintf(stdout,"%s",temps); // want to see what temps is

      sprintf(filename,"%s0_analdata%s%s",temps,DAT2EXT,myidtxt);
            if((analyticout=fopen(filename,"wt"))==NULL){
	fprintf(fail_file,"Cannot open %s\n",temps);
       	myexit(1);
      }
    }

    if(RYUJONES){
      // use x1,x2=0..1  tf=.48
      v0=1.0;
      b0=1.0/sqrt(4.0*M_PI);
      gam=5.0/3.0;
      mach=1.0;
      beta=10.0/3.0;
      p0=beta*b0*b0*0.5;
    
      rho0=mach*gam*p0/(v0*v0);
      argsin=2.0*M_PI;
    }
    else{
      // use x1,x2=0..M_PI tf=3.1
      v0=1.0;
      b0=1.0;
      gam=5.0/3.0;
      mach=1.0;
      beta=10.0/3.0;
      p0=5.0/3.0;
      rho0=25.0/9.0;
      argsin=1.0;
    }

    en0=p0/(gam-1.0);

    LOOPF{
      sanal[1][k][j][i]=rho0;
      sanal[2][k][j][i]=en0;
      sanal[3][k][j][i]=0.0;

      if(VORTEXPLANE==3){
	dir[1]=1;
	dir[2]=2;
	dir[3]=3;
	xa=x[1][1][i];
	xb=x[2][1][i];
	yb=x[2][2][j];
	ya=x[1][2][j];
      }
      else if(VORTEXPLANE==2){
	dir[1]=1;
	dir[2]=3;
	dir[3]=2;
	xa=x[1][1][i];
	xb=x[2][1][i];
	yb=x[2][3][k];
	ya=x[1][3][k];
      }
      else if(VORTEXPLANE==1){
	dir[1]=2;
	dir[2]=3;
	dir[3]=1;
	xa=x[1][2][j];
	xb=x[2][2][j];
	yb=x[2][3][k];
	ya=x[1][3][k];
      }

      if(VORTEXROTATION==0){
	vanal[1][dir[1]][k][j][i]=v0*(-sin(argsin*yb));
	vanal[1][dir[2]][k][j][i]=v0*(sin(argsin*xb));
	vanal[1][dir[3]][k][j][i]=0.0;
      
	//vanal[2][dir[1]][k][j][i]=0.0;
	//vanal[2][dir[2]][k][j][i]=0.0;
	vanal[2][dir[1]][k][j][i]=b0*(-sin(argsin*yb));
	vanal[2][dir[2]][k][j][i]=b0*(sin(2.0*argsin*xb));
	vanal[2][dir[3]][k][j][i]=0.0;
      }

      if(VORTEXROTATION==1){
	vanal[1][dir[1]][k][j][i]=v0*(-sin(argsin*xa));
	vanal[1][dir[2]][k][j][i]=v0*(-sin(argsin*ya));
	vanal[1][dir[3]][k][j][i]=0.0;
      
	//vanal[2][dir[1]][k][j][i]=0.0;
	//vanal[2][dir[2]][k][j][i]=0.0;
	vanal[2][dir[1]][k][j][i]=b0*(-sin(2.0*argsin*xa));
	vanal[2][dir[2]][k][j][i]=b0*(-sin(argsin*ya));
	vanal[2][dir[3]][k][j][i]=0.0;
      }

      if(VORTEXROTATION==2){
	vanal[1][dir[1]][k][j][i]=v0*(sin(argsin*yb));
	vanal[1][dir[2]][k][j][i]=v0*(-sin(argsin*xb));
	vanal[1][dir[3]][k][j][i]=0.0;
      
	//vanal[2][dir[1]][k][j][i]=0.0;
	//vanal[2][dir[2]][k][j][i]=0.0;
	vanal[2][dir[1]][k][j][i]=b0*(sin(argsin*yb));
	vanal[2][dir[2]][k][j][i]=b0*(-sin(2.0*argsin*xb));
	vanal[2][dir[3]][k][j][i]=0.0;
      }

      if(VORTEXROTATION==3){
	vanal[1][dir[1]][k][j][i]=v0*(sin(argsin*xa));
	vanal[1][dir[2]][k][j][i]=v0*(sin(argsin*ya));
	vanal[1][dir[3]][k][j][i]=0.0;
      
	//vanal[2][dir[1]][k][j][i]=0.0;
	//vanal[2][dir[2]][k][j][i]=0.0;
	vanal[2][dir[1]][k][j][i]=b0*(sin(2.0*argsin*xa));
	vanal[2][dir[2]][k][j][i]=b0*(sin(argsin*ya));
	vanal[2][dir[3]][k][j][i]=0.0;
      }

      if(VORTEXROTATION==5){ // not vortex, but shows rhointerp problem
	vanal[1][dir[1]][k][j][i]=v0*(-sin(argsin*xa));
	vanal[1][dir[2]][k][j][i]=v0*(-sin(argsin*ya));
	vanal[1][dir[3]][k][j][i]=0.0;
      
	//vanal[2][dir[1]][k][j][i]=0.0;
	//vanal[2][dir[2]][k][j][i]=0.0;
	vanal[2][dir[1]][k][j][i]=b0*(-sin(argsin*xa));
	vanal[2][dir[2]][k][j][i]=b0*(-sin(2.0*argsin*ya));
	vanal[2][dir[3]][k][j][i]=0.0;
      }
    
    }


    // output some interesting data
    // and setup initial problem stuff
    if(firstsolve==1){
      firstsolve=0;
      if(RYUJONES){
	mms[0][0][1][0]=     0.103677772 ;
	mms[0][1][1][0]=     0.103677772 ;
	mms[0][0][1][1]=    0.6175035238 ;
	mms[0][1][1][1]=    0.6175035238 ;
	mms[0][0][2][0]=   0.06121816859 ;
	mms[0][1][2][0]=   0.06121816859 ;
	mms[0][0][2][1]=    0.7375584841 ;
	mms[0][1][2][1]=    0.7375584841 ;
	mms[0][0][3][0]=               0 ;
	mms[0][1][3][0]=               0 ;
	mms[0][0][3][1]=               0 ;
	mms[0][1][3][1]=               0 ;
	mmv[0][0][1][1][0]=    -10.211948156 ;
	mmv[0][1][1][1][0]=    -10.211948156 ;
	mmv[0][0][1][1][1]=     1.207605004 ;
	mmv[0][1][1][1][1]=     1.207605004 ;
	mmv[0][0][1][2][0]=   -0.9498059154 ;
	mmv[0][1][1][2][0]=   -0.9498059154 ;
	mmv[0][0][1][2][1]=    0.9503436685 ;
	mmv[0][1][1][2][1]=    0.9503436685 ;
	mmv[0][0][1][3][0]=               0 ;
	mmv[0][1][1][3][0]=               0 ;
	mmv[0][0][1][3][1]=               0 ;
	mmv[0][1][1][3][1]=               0 ;
	mmv[0][0][2][1][0]=   -0.6236243844 ;
	mmv[0][1][2][1][0]=   -0.6236243844 ;
	mmv[0][0][2][1][1]=     0.623854816 ;
	mmv[0][1][2][1][1]=     0.623854816 ;
	mmv[0][0][2][2][0]=    -0.572280705 ;
	mmv[0][1][2][2][0]=    -0.572280705 ;
	mmv[0][0][2][2][1]=    0.5630972981 ;
	mmv[0][1][2][2][1]=    0.5630972981 ;
	mmv[0][0][2][3][0]=               0 ;
	mmv[0][1][2][3][0]=               0 ;
	mmv[0][0][2][3][1]=               0 ;
	mmv[0][1][2][3][1]=               0 ; 
      }
      else{
	mms[0][0][1][0]=    0.8;
	mms[0][1][1][0]=    0.8;
	mms[0][0][1][1]=     8.0;
	mms[0][1][1][1]=     8.0;
	mms[0][0][2][0]=    0.4;
	mms[0][1][2][0]=    0.4;
	mms[0][0][2][1]=     8.0;
	mms[0][1][2][1]=     8.0;
	mms[0][0][3][0]=               0 ;
	mms[0][1][3][0]=               0 ;
	mms[0][0][3][1]=               0 ;
	mms[0][1][3][1]=               0 ;
	mmv[0][0][1][1][0]=    -1.5;
	mmv[0][1][1][1][0]=    -1.5;
	mmv[0][0][1][1][1]=     1.5;
	mmv[0][1][1][1][1]=     1.5;
	mmv[0][0][1][2][0]=    -1.5;
	mmv[0][1][1][2][0]=    -1.5;
	mmv[0][0][1][2][1]=     1.5;
	mmv[0][1][1][2][1]=     1.5;
	mmv[0][0][1][3][0]=               0 ;
	mmv[0][1][1][3][0]=               0 ;
	mmv[0][0][1][3][1]=               0 ;
	mmv[0][1][1][3][1]=               0 ;
	mmv[0][0][2][1][0]=    -3.0;
	mmv[0][1][2][1][0]=    -3.0;
	mmv[0][0][2][1][1]=     3.0;
	mmv[0][1][2][1][1]=     3.0;
	mmv[0][0][2][2][0]=    -3.0;
	mmv[0][1][2][2][0]=    -3.0;
	mmv[0][0][2][2][1]=     3.0;
	mmv[0][1][2][2][1]=     3.0;
	mmv[0][0][2][3][0]=               0 ;
	mmv[0][1][2][3][0]=               0 ;
	mmv[0][0][2][3][1]=               0 ;
	mmv[0][1][2][3][1]=               0 ; 
	/*
	  mms[0][0][1][0]=    0.8747499479 ;
	  mms[0][1][1][0]=    0.8747499479 ;
	  mms[0][0][1][1]=     18.24006259 ;
	  mms[0][1][1][1]=     18.24006259 ;
	  mms[0][0][2][0]=    0.4605765082 ;
	  mms[0][1][2][0]=    0.4605765082 ;
	  mms[0][0][2][1]=     14.68577011 ;
	  mms[0][1][2][1]=     14.68577011 ;
	  mms[0][0][3][0]=               0 ;
	  mms[0][1][3][0]=               0 ;
	  mms[0][0][3][1]=               0 ;
	  mms[0][1][3][1]=               0 ;
	  mmv[0][0][1][1][0]=    -1.475806569 ;
	  mmv[0][1][1][1][0]=    -1.475806569 ;
	  mmv[0][0][1][1][1]=     1.475806569 ;
	  mmv[0][1][1][1][1]=     1.475806569 ;
	  mmv[0][0][1][2][0]=    -1.311625211 ;
	  mmv[0][1][1][2][0]=    -1.311625211 ;
	  mmv[0][0][1][2][1]=     1.311625211 ;
	  mmv[0][1][1][2][1]=     1.311625211 ;
	  mmv[0][0][1][3][0]=               0 ;
	  mmv[0][1][1][3][0]=               0 ;
	  mmv[0][0][1][3][1]=               0 ;
	  mmv[0][1][1][3][1]=               0 ;
	  mmv[0][0][2][1][0]=    -6.009352377 ;
	  mmv[0][1][2][1][0]=    -6.009352377 ;
	  mmv[0][0][2][1][1]=     6.009352472 ;
	  mmv[0][1][2][1][1]=     6.009352472 ;
	  mmv[0][0][2][2][0]=    -6.053225136 ;
	  mmv[0][1][2][2][0]=    -6.053225136 ;
	  mmv[0][0][2][2][1]=     6.053224113 ;
	  mmv[0][1][2][2][1]=     6.053224113 ;
	  mmv[0][0][2][3][0]=               0 ;
	  mmv[0][1][2][3][0]=               0 ;
	  mmv[0][0][2][3][1]=               0 ;
	  mmv[0][1][2][3][1]=               0 ; 
	*/
	/*
	  mms[0][0][1][0]=     1.125885617 ;
	  mms[0][1][1][0]=     1.125885617 ;
	  mms[0][0][1][1]=     13.36235639 ;
	  mms[0][1][1][1]=     13.36235639 ;
	  mms[0][0][2][0]=    0.6858458864 ;
	  mms[0][1][2][0]=    0.6858458864 ;
	  mms[0][0][2][1]=     11.40762333 ;
	  mms[0][1][2][1]=     11.40762333 ;
	  mms[0][0][3][0]=               0 ;
	  mms[0][1][3][0]=               0 ;
	  mms[0][0][3][1]=               0 ;
	  mms[0][1][3][1]=               0 ;
	  mmv[0][0][1][1][0]=    -1.027579457 ;
	  mmv[0][1][1][1][0]=    -1.027579457 ;
	  mmv[0][0][1][1][1]=     1.011192604 ;
	  mmv[0][1][1][1][1]=     1.011192604 ;
	  mmv[0][0][1][2][0]=    -1.053431636 ;
	  mmv[0][1][1][2][0]=    -1.053431636 ;
	  mmv[0][0][1][2][1]=     1.032067915 ;
	  mmv[0][1][1][2][1]=     1.032067915 ;
	  mmv[0][0][1][3][0]=               0 ;
	  mmv[0][1][1][3][0]=               0 ;
	  mmv[0][0][1][3][1]=               0 ;
	  mmv[0][1][1][3][1]=               0 ;
	  mmv[0][0][2][1][0]=    -1.865236614 ;
	  mmv[0][1][2][1][0]=    -1.865236614 ;
	  mmv[0][0][2][1][1]=     1.924464629 ;
	  mmv[0][1][2][1][1]=     1.924464629 ;
	  mmv[0][0][2][2][0]=    -1.933663744 ;
	  mmv[0][1][2][2][0]=    -1.933663744 ;
	  mmv[0][0][2][2][1]=      1.87284829 ;
	  mmv[0][1][2][2][1]=      1.87284829 ;
	  mmv[0][0][2][3][0]=               0 ;
	  mmv[0][1][2][3][0]=               0 ;
	  mmv[0][0][2][3][1]=               0 ;
	  mmv[0][1][2][3][1]=               0 ; 
	*/
      }
      // now output some interesting analytic data

      fprintf(analyticout,"#hello!\n");
      fprintf(analyticout,"%15s %15s %15s %15s %15s %15s %15s %15s\n","v0","b0","gam","mach","beta","p0","en0","rho0");
      fprintf(analyticout,"%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",v0,b0,gam,mach,beta,p0,en0,rho0);
      /*
        LOOPF{
	fprintf(analyticout,"%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g\n",sanal[1][k][j][i],sanal[2][k][j][i],sanal[3][k][j][i],vanal[1][1][k][j][i],vanal[1][2][k][j][i],vanal[1][3][k][j][i],vanal[2][1][k][j][i],vanal[2][2][k][j][i],vanal[2][3][k][j][i]);
	}
      */
      fclose(analyticout);
    }

  }// end else if not calltype==100
}



// parameters as in S&N92
void magcorona(int calltype)
{ // assume COORD==3

  static FILE*analyticout;
  static int firstsolve=1;
  char temps[100];
  char filename[100];
  SFTYPE ftemp;
  SFTYPE pos;
  int i,j,k;
  SFTYPE t0,lambda,A0,P0,escic,eta,nu,M,R0,mp,G;
  SFTYPE phi,ecsi,rc,rs,d0;

  if(calltype==100){
    // corona
    x1in=1.0E11; x1out=5.0E11; // full Pi in theta as above
    // for corona
    tstart=8740;
    tf=2.62E4;
  }
  else if(calltype==0){
    // want to write some interesting data on solution
    if(firstsolve==1){
    
      strcpy(temps,DATADIR);

      sprintf(filename,"%s0_analdata%s%s",temps,DAT2EXT,myidtxt);
      if((analyticout=fopen(filename,"wt"))==NULL){
	fprintf(fail_file,"Cannot open %s\n",temps);
	myexit(1);
      }
    }

    t=t0=8.740E3; // sec
    //tstart=8740; // in init.c
    //tf=2.62E4;
    //x1in=1.0E11; x1out=5.0E11; // in init.c

    lambda=5.54E11; // cm^-1
    A0=1.5E21; // G cm^2
    P0=1.01327;
    escic=1.104E11; // cm
    eta=5.24E-8; // sec^-2
    nu=2.42E20; // cm^2 s^-2 gm^-1/3
    M=2.0E33; // gm
    R0=1.0E11; // cm
    mp=1.673E-24; // gm
    G=6.67E-8; // dynes cm^2 g^-2



    LOOPF{

      phi=sqrt(eta)*t;
      ecsi=R0*pow(phi,1.0/6.0);
      rc=escic*phi;
      rs=R0*pow(phi,7.0/6.0);
      d0=1E8*mp*exp(-2.0*G*M/(3.0*eta*R0*R0*R0));
    

    

      sanal[1][k][j][i]=0.0;
      sanal[2][k][j][i]=0.0;
      sanal[3][k][j][i]=0.0;
    
      vanal[1][1][k][j][i]=0.0;
      vanal[1][2][k][j][i]=0.0;
      vanal[1][3][k][j][i]=0.0;

      vanal[2][1][k][j][i]=0.0;
      vanal[2][2][k][j][i]=0.0;
      vanal[2][3][k][j][i]=0.0;
    
    
    }


    // output some interesting data
    // and setup initial problem stuff
    if(firstsolve==1){
      firstsolve=0;
  
      // now output some interesting analytic data

      fprintf(analyticout,"#hello!\n");
      /*
        LOOPF{
	fprintf(analyticout,"%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g\n",sanal[1][k][j][i],sanal[2][k][j][i],sanal[3][k][j][i],vanal[1][1][k][j][i],vanal[1][2][k][j][i],vanal[1][3][k][j][i],vanal[2][1][k][j][i],vanal[2][2][k][j][i],vanal[2][3][k][j][i]);
	}
      */
      fclose(analyticout);
    }

  }
}




#define  MAGBONDI  0
// 0: normal bondi is full spherical problem (solves for velocity profile)
// 1: MHD bondi(webber-davis) (solves for density profile)

void bondisol(int calltype,FILE*analyticoutreal)
{ // solution is on b-grid as normal

  static int firstsolve=1;
  static FILE*analyticout;
  static FILE*analyticout2;
  SFTYPE time2fall,cscrosstime,upstream,downstream;
  SFTYPE min1a_full,mout1a_full,kein1a_full,keout1a_full,hin1a_full,hout1a_full,pein1a_full,peout1a_full; // analytic values

  SFTYPE min1,mout1,kein1,keout1,hin1,hout1,pein1,peout1;
  SFTYPE mout1_full,keout1_full,hout1_full,peout1_full,min1_full,kein1_full,hin1_full,pein1_full;

  SFTYPE tempf;

  SFTYPE RSTAR,CS_RSTAR,V_RSTAR,RHO_RSTAR,KAPPA,MACH;
  int SELFSIM,TYPEKAPPA;
  SFTYPE Bern,Bec,RATIO,N_x1,N_x2,EPS,XS,VS,SS,NDOT;
  SFTYPE fractioneps;
  int i,j,k;
  int l,m,z;
  SFTYPE ndot;
  SFTYPE ep,epa,epb,dep;
  SFTYPE da;
  SFTYPE which;
  SFTYPE am,amm1,al;
  SFTYPE fold,fnew;
  int gotit,gotgood;
  SFTYPE (*V1grid)[N2M][N1M];
  SFTYPE (*V2grid)[N2M][N1M];
  SFTYPE (*V3grid)[N2M][N1M];
  SFTYPE (*V4grid)[N2M][N1M];
  int gridsize;
  char temps[100];
  char filename[100];
  SFTYPE ftemp,ftemp1,ftemp2;
  int usealf;
  SFTYPE magk,magf,magOmega,magphi,velra,magenergy;
  SFTYPE (*works1a)[N2M][N1M];
  SFTYPE (*works2a)[N2M][N1M];
  SFTYPE (*works3a)[N2M][N1M];
  SFTYPE (*works4a)[N2M][N1M];
  SFTYPE (*works1)[N2M][N1M];
  SFTYPE (*works2)[N2M][N1M];
  SFTYPE (*works3)[N2M][N1M];
  SFTYPE (*works4)[N2M][N1M];
  SFTYPE MAGTHETA,MAGomega,SPFAST,SPALF,SPSLOW,RHORA,MAGRA,MAGK,MAGF,MAGOMEGA,MAGPHI,VELRA,MAGENERGY;
  int GSI,GSM,GSM2,GSM3,MAXZOOMS,MAXSOLNS;

  SFTYPE vr,vtheta,vphi,br,btheta,bphi;
  FTYPE (*sca3)[N3M][N2M][N1M];
  FTYPE (*vx3)[N3M][N2M][N1M];
  FTYPE (*vy3)[N3M][N2M][N1M];
  FTYPE (*vz3)[N3M][N2M][N1M];

  FTYPE (*radiustouse)[N2M][N1M];
  FTYPE (*thetatouse)[N2M][N1M];
  FTYPE (*phitouse)[N2M][N1M];





  
  GSI= 100;
  GSM= 1000;
  GSM2= 10000;
  GSM3= 100000;
  MAXZOOMS= 50;
  // MAXZOOMS= 200;
  MAXSOLNS= 2;

  if(calltype==100){

    if(COORD==3){
      x3in=0;
      x3out=2.0*M_PI;
      if(MAGBONDI){
	mag=1;
	rgp=0.0;
	gam=1.2;
	//x1in=0.4; x1out=3.0;
	//x1in=200.0; x1out=1000.0;
	x1in=290.0; x1out=915.0;
	x2in=M_PI/2-M_PI*.00001;
	x2out=M_PI/2+M_PI*.00001;
      }
      else{
	rgp=rg;
	mag=0;
	gam=5.0/3.0;
	x1in=1.4*rg;
	x1out=100.0*rg;
	x2in=0;
	x2out=M_PI;
      }
      tf=1.0E5;
    }
    else if(COORD==1){
      Rinner=6*rg;
      Router=100*rg;
      gam=4.0/3.0;
      rgp=rg;
      mag=0;
      x1in=-100*rg;
      x1out=100*rg;
      x2in=-100*rg;
      x2out=100*rg;
      x3in=-100*rg;
      x3out=100*rg;
      tf=1.0E4;
    }

    nonunigridx1=0;
    nonunigridx2=0;
    nonunigridx3=0;

    visc_real=0;
    res_real=0;
    mdotin=0;
    cool=0;

    // for BOUNDTYPE==1
    if(BOUNDTYPE==1){
      simplebc=1;
      bcix1=4;
      bcox1=3;
      bcix2=1;
      bcox2=1;
    }
    //    DTl=DTener=DTloss=DTtimestep=DTsp=tf/1000.0;  DTfloor=DTd=tf/(20.0); DTpd=tf/(10.0); DTi=tf;
    DTl=DTener=DTloss=tf/1000.0; DTtimestep=DTsp=tf/100.0; DTfloor=DTd=DTmode=tf/(20.0); DTpd=tf/(10.0); DTi=tf/(100.0);
    // GODMARK
    DTloss=tf/1000000.0;
  }
  else if((calltype==0)||(calltype==2)){


    sca3=workv2;
    vx3=workv3;
    vy3=workv4;
    vz3=workv5;

    grids_cart2spc(sca3,vx3,vy3,vz3);
    
    if((works1a=(SFTYPE (*) [N2M][N1M])malloc(sizeof(SFTYPE)*N3M*N2M*N1M))==NULL){
      fprintf(fail_file,"Can't open bondi works1 memory\n");
      myexit(45);
    }
    if((works2a=(SFTYPE (*) [N2M][N1M])malloc(sizeof(SFTYPE)*N3M*N2M*N1M))==NULL){
      fprintf(fail_file,"Can't open bondi works2 memory\n");
      myexit(45);
    }
    if((works3a=(SFTYPE (*) [N2M][N1M])malloc(sizeof(SFTYPE)*N3M*N2M*N1M))==NULL){
      fprintf(fail_file,"Can't open bondi works3 memory\n");
      myexit(45);
    }
    if((works4a=(SFTYPE (*) [N2M][N1M])malloc(sizeof(SFTYPE)*N3M*N2M*N1M))==NULL){
      fprintf(fail_file,"Can't open bondi works4 memory\n");
      myexit(45);
    }

    works1=(SFTYPE (*) [N2M][N1M])(&(works1a[N3BND][N2BND][N1BND]));
    works2=(SFTYPE (*) [N2M][N1M])(&(works2a[N3BND][N2BND][N1BND]));
    works3=(SFTYPE (*) [N2M][N1M])(&(works3a[N3BND][N2BND][N1BND]));
    works4=(SFTYPE (*) [N2M][N1M])(&(works4a[N3BND][N2BND][N1BND]));

    V1grid=works1;
    V2grid=works2;
    V3grid=works3;
    V4grid=works4;
  

    R0=L[1][1]+L[2][1];

    // want to write some interesting data on solution
    if(firstsolve==1){    
      strcpy(temps,DATADIR);

      if(myid<=0){
      
	sprintf(filename,"%s0_analdata2%s%s",temps,DAT2EXT,myidtxt);
	if((analyticout2=fopen(filename,"wt"))==NULL){
	  fprintf(fail_file,"Cannot open %s\n",temps);
	  myexit(1);
	}
	if(calltype!=2){
	  strcpy(temps,DATADIR);
	
	  sprintf(filename,"%s0_analdata%s",temps,DAT2EXT);
	  if((analyticout=fopen(filename,"wt"))==NULL){
	    fprintf(fail_file,"Cannot open %s\n",temps);
	    myexit(1);
	  }	
	}
	else analyticout=analyticoutreal;
      }
    }
  
  
    if(COORD==2){
      fprintf(fail_file,"Analytic Bondi solution not setup for cylindrical coordinates\n");
      myexit(1);
    }

    ////////////////////////////////
    //
    // normal bondi problem define
    //
    //
    // Fy(y,ep) = 0 is a solution to this system
    // Y is vel
    // S is rho
    // X is cs
    // ep=epsilon=r
  
    // some notes:
    //
    // for  1<gam>=5/3 for Be=1/(gam-1) normal general solution can be used(NDOT/FY,etc)
    //
    // for gam==1 Be=Anything can be converted to Be=0, then ndot xor SS changes functional form, both change numerical value based upon eachother.
    // I choose to change SSs function form, keeping ndots functional form the same

    // given star values may not resolve to a physical solution

    GM=1.0; // general setup

    SELFSIM=0;  // whether to choose gam=5/3 as self-similar solution or not(only relevant when !wpw and !gam53)
  

		// init things
    RSTAR=0.0;
    TYPEKAPPA=0;
    V_RSTAR=0.0;
    CS_RSTAR=0.0;
    KAPPA=0.0;
    Bern=0.0;
    Bec=0.0;
    RATIO=0.0;
    N_x1=0.0;
    N_x2=0.0;
    EPS=0.0;
    XS=VS=SS=ndot=0.0;

    // for gam=5/3 SS,  the CS_RSTAR and V_RSTAR must be nonzero and defined as a ratio for self-similar case
    if( (!wgam53)||(wgam53&&( wpw||( (!wpw)&&(!SELFSIM) ) ) ) ){ // if not chosen self-similar, gam=5/3,non-PW problem

      if(wgam1){// gam==1
      
	RSTAR=(1.E+23); // position which defines Be and KAPPA

	TYPEKAPPA=1; // whether entropy will be 0 of gas or not
	CS_RSTAR=cs; // so set cs like normal for gam=1
  
	if(TYPEKAPPA==1){
	
	  V_RSTAR=(0.0); // velocity at infinity(used to control Be in general), although effect is absorbed into SS's functional form, so as if still no v at infinity // units in number of c_s for gam==1, in general units for other gam
	  KAPPA=(CS_RSTAR*CS_RSTAR); // for gam=1 gas
	  RHO_RSTAR=1.0;
	}
	else if(TYPEKAPPA==0){
	  // simple solution in this case actually, for any gamma, but not done yet
	  KAPPA=(0.0);
	  RHO_RSTAR=(1.0);
	}
      
	Bern=( 0.5*V_RSTAR*V_RSTAR+CS_RSTAR*CS_RSTAR*log(RHO_RSTAR)-GM/(RSTAR-rgp) ); // Bern param for gam=1
	RATIO=(0.25*GM/(cs*cs)); // special reoccuring ratio for gam=1
	EPS=( RATIO*(1.+rgp/RATIO+sqrt(1.+2.*rgp/RATIO)) ); // sonic point for gam=1(EPS would work if Bec==1
	XS=(cs);
	SS=( exp(Bern/(cs*cs)-0.5+GM/(cs*cs*(EPS-rgp))) ); // sonic point density for gam=1
      }
      else{// gam!=1

	RSTAR=20.0; // position which defines Be and KAPPA
	TYPEKAPPA=1; // whether entropy will be 0 of gas or not
      
	if(TYPEKAPPA==1){
	
	  V_RSTAR=(-.5); // velocity at rstar(used to control Be in general), although effect is absorbed into SS's functional form, so as if still no v at infinity // units in number of c_s for gam==1, in general units for other gam
	  CS_RSTAR=(.5);

	  if(calltype!=2){
	    RHO_RSTAR=(1.0);
	  }
	  else RHO_RSTAR=1E-4;
	
	  // set KAPPA which is equivelant to setting CS(r*) and RHO(r*), where r* is some defined location where each is nonzero(for cs_inf=1,rho_inf=1 K=1/gam) (if KAPPA==0, must set CS and RHO instead)
	  //KAPPA=( CS_RSTAR*CS_RSTAR/(gam*pow(RHO_RSTAR,gam-1.)) );
	  //RHO_RSTAR=( pow(CS_RSTAR*CS_RSTAR/(gam*KAPPA),1./(gam-1.)) ); // or set KAPPA
	  //KAPPA=(1.0/gam); // cs_inf^2=gam*K*\rho_inf^{gam-1}, and when inf's are 1 (i.e. in those units which is possible here) then K=1/gam, then v_r_inf=0
	}
	else if(TYPEKAPPA==0){
	  // simple solution in this case actually, for any gamma, but not done yet
	  KAPPA=(0.0);
	  RHO_RSTAR=(1.0);
	  CS_RSTAR=(0.0);
	}
	// one can override this Bern since STAR stuff not used after this
	//      Bern=( 0.5*V_RSTAR*V_RSTAR+CS_RSTAR*CS_RSTAR/(gam-1.)-GM/(RSTAR-rgp) ); // Bernoulli parameter otherwise
	//      KAPPA=100.0;
	//if(calltype==2)      Bern=.036; // for bondi+torus solution
	//else Bern=.00001; // good for gam=5/3 to keep stable for long time
      
	//Bern=.001; // bad for gam=5/3, unstable

	//Bern=2.36295E-4; // rs=10*R0
	//Bern=7.08617E-3; // rs=2*R0
	Bern=3.6523E-2;// rs=R0
	//Bern=2.74348E-1; // rs=R0/2
	//Bern=2.0E4; // rs=R0/10

	//first rs=R0*10 (10^(-5) 10^(-3) 10^(-1) )
	//KAPPA=7.17806;
	//KAPPA=.333176;
	//KAPPA=.0154647;
	//second rs=R0*2
	//KAPPA=43.0521;
	//KAPPA=1.9983;
	//KAPPA=.0927529;
	//third rs=R0
	KAPPA=110.948;
	//KAPPA=5.14974;
	//KAPPA=.23903;
	//fourth rs=R0/2
	//KAPPA=416.701;
	//KAPPA=19.3415;
	//KAPPA=.897755;
	//fifth rs=R0/10
	//KAPPA=1.29912E7;
	//KAPPA=603000;
	//KAPPA=27988.8;


	//other tests
	//      KAPPA=2390.3; // rho(EPS)=1E-7
	//				KAPPA=23.903; // rho(EPS)=1E-4



	Bec=( Bern*(gam-1.) ); // reduced Be parameter

	if(fabs(Bec)>1E-10){
	  if( (Bec<0)&&(wgam53)){
	    fprintf(fail_file,"\nNo physical solution for gam=5/3 and Bern<0, only self-similar -> Bern=0 !\n");
	    myexit(1);
	  }
	  RATIO= (8.*Bec/GM); // special reoccuring ratio
	  // sonic point values
	  N_x1=(5.-3.*gam)/RATIO;
	  N_x2=(N_x1*N_x1+2.*(gam+1.)*rgp/RATIO);
	  EPS=(N_x1+rgp+sqrt(N_x2) ); // sonic point for 1<gam<=5/3
	}
	else{// special case for Bern=0, doesn't work for gam=5/3
	  EPS=4.*rgp*(gam-1.)/(5.-3.*gam);
	  if(wgam53){
	    fprintf(fail_file,"\nNo physical solution for gam=5/3 and Bern=0, only self-similar!\n");
	    myexit(1);
	  }
	}
	XS=(sqrt(0.5*EPS*GM)/(EPS-rgp)); // sonic point sound speed
	SS=( pow(XS,2./(gam-1.))/pow(gam*KAPPA,1./(gam-1.)) ); // sonic point density
      }
      VS=-XS; // velocity at sonic point for any gamma
    
      ndot=(4.*M_PI*EPS*EPS*SS*VS); // for any gamma
      // below form correct, but not necessarily a solution for given star values
      //ndot=4.*M_PI*RSTAR*RSTAR*V_RSTAR*RHO_RSTAR; // simplest form    
    }
    // else if gam=5/3 selfsim solution
    else{// generally different setup, so have different assignment section

      RSTAR=(50.0); // position which defines Be and KAPPA
  
      // Kappa must be nonzero for SS case

      V_RSTAR=(-10.0); // velocity at infinity(used to control Be in general), although effect is absorbed into SS's functional form, so as if still no v at infinity // units in number of c_s for gam==1, in general units for other gam
      CS_RSTAR=(1.0);
      RHO_RSTAR=(1.0);
      KAPPA=( CS_RSTAR*CS_RSTAR/(gam*pow(RHO_RSTAR,gam-1.)) ); // which is equivelant to setting CS(r*) and RHO(r*), where r* is some defined location where each is nonzero(for cs_inf=1,rho_inf=1 K=1/gam)
      //    RHO_RSTAR=( pow(CS_RSTAR*CS_RSTAR/(gam*KAPPA),1./(gam-1.)) ); // or set KAPPA
      MACH=(V_RSTAR/CS_RSTAR); // M=A/B

      XS=VS=SS=1.E+30;

      Bern=0.0; // must be 0 for self-similar solution
      ftemp1 = sqrt(2.0*GM*MACH*MACH/(3.0+MACH*MACH)); // -v_r prefactor
      RATIO=RSTAR*pow(V_RSTAR/ftemp1,2.); // scaling relation for SS solution(x/Q=ep)
      // ratio of 1.0 seems to do very well(last longer before overfilling with mass), may just be time thing
      // higher mach's than 2.0 do very well too indep of ratio.  This probably has to do with material bouncing off the inner grid
      // outflow boundary on inflow case does very poorly

      // ndot
      ndot=4.*M_PI*RSTAR/RATIO*RSTAR/RATIO*V_RSTAR*RHO_RSTAR; // simplest form
    }



    // check if solutions are valid

    if(!wgam1){
      if(N_x2<0){
	fprintf(fail_file,"N_x2 < 0 leading to Im values: N_x2: %15.10g\n",N_x2);
	fprintf(fail_file,"Must choose larger Bern perhaps\n");
	myexit(1);
      }
      if(wpw){
	ftemp=-(N_x1+sqrt(N_x2)/rgp);
	if(RATIO<ftemp){
	  fprintf(fail_file,"RATIO: %15.10g < %15.10g leading to rs<0 : %15.10g\n",RATIO,ftemp,EPS);
	  fprintf(fail_file,"Must choose larger Bern perhaps\n");
	  myexit(1);
	}
      }
    }


#if(MAGBONDI==0)
    //#define NDOT53 (M_PI/4.*pow(2.+2.*sqrt(6.*rgp)+3.*rgp,2.)) // selfsim case

    // used for non-self-sim case
#define FY(ndot,y,ep) ( 0.5*y*y+gam*KAPPA/(gam-1.)*pow(fabs(ndot/(4.*M_PI*ep*ep*y)),(gam-1.))-GM/(ep-rgp)-Bern )
#define FY1(ndot,y,ep) ( 0.5*y*y+cs*cs*log(fabs(ndot/(4.*M_PI*ep*ep*y)))-GM/(ep-rgp)-Bern )
    //#define FY53(ndot,y,ep) (0.5*y*y+1.5*pow(ndot/(4.*M_PI*ep*ep*y),2.0/3.0)-1./(ep-rgp)) // self-sim case
#define S(ndot,y,ep) (ndot/(4.*M_PI*ep*ep*y))
#define X(s) (sqrt(gam*KAPPA)*pow(s,0.5*(gam-1.)))

    // since FY in original form is an ellipse in X and Y, there exist known max and mins for X and Y:
    // We will only investigate the upper X,Y quandrant where the ndot solution exists, and since X must be positive
#define Ymax(ep) ( sqrt(2.*(Bern+GM/(ep-rgp))) )
#define Ymin(ep) (.00000001)
#define Xmax(ep) ( sqrt((gam-1.)*(Bern+GM/(ep-rgp))) )
#define Xmin(ep) (.00000001)
    // override, probably need to increase resolution too
    //#define Ymax(ep) 100
    //#define Xmax(ep) 100

    // for gam=1
    // no known general condition
#define Ymax1(ep) (20.0+( sqrt(2.*(Bern+GM/(ep-rgp))) ) ) // hack
#define Ymin1(ep) (.00000001)
#define Xmax1(ep) (2.0*cs) // really this is cs, no?
#define Xmin1(ep) (.00000001)
#define LOWESTY (.00001)
    // given the above, for each ep there should be 2 y, 2 s, and 2 x, representing the 2 solutions that pass through the sonic point, which is assumed in the derivation of the above.
    // We want the solution that has y(ep) monotonicly increasing.  This the the accretion solution.  The wind solution has y(ep) decreasing monotonically.


    //  End of definition of normal bondi problem
    //
    ////////////////////////
#else
    ///////////////////////////
    // magbondi
    // magnetic bondi problem(WD/Sakuri85/S&N92)
    //
    //

    // some of these should really be determined in general so can have more than 1 case
    MAGTHETA =(0.5); // dimenless
    MAGomega =(.25); // dimenless
    // use mathematica to solve, very fast, easy, and can get *solution* to high precision.
    // magbondi2.nb
    //{xs -> 0.7770330780030363194074308952074505, ys -> 1.940544269119426993709021901238330, xf -> 1.302391299892080014807550863286971, yf -> 0.5139371598659040285268997762736311, En -> 1.738429467866434807018809306133664, beta -> 0.5755932415448330480108551469612932}
    // order from outer to inner radius is FAST->ALF->SLOW in increasing density value
    SPFAST=(1.302391299892080014807550863286971); // y=0.5139371598659040285268997762736311
    SPALF=(1.0); // y=1.0
    SPSLOW=(0.7770330780030363194074308952074505); // y=1.940544269119426993709021901238330

    MAGBETA =(0.5755932415448330480108551469612932); // dimenless
    MAGE =(1.738429467866434807018809306133664); // dimenless
  
    RHORA= (1.0E-13);
    MAGRA =(500.0);
    //	MAGRA =(1.0);
  
    // use magbondi3.nb and use above mathematica solution to solve for these
    // no, use general solution and write to real variables before use
    MAGK =(pow(RHORA,1.0-gam)*MAGTHETA/(gam*MAGRA));
    MAGF =(RHORA*sqrt(MAGBETA*GM*pow(MAGRA,3.0)));
    MAGOMEGA =(sqrt(GM*MAGomega/pow(MAGRA,3.0)));
    //	MAGPHI=(sqrt(4.0*M_PI*MAGBETA*GM*pow(MAGRA,3.0)*RHORA));
    MAGPHI (sqrt(MAGBETA*GM*pow(MAGRA,3.0)*RHORA));
    VELRA =(MAGF/(MAGRA*MAGRA*RHORA));
    MAGENERGY =(MAGE*GM/MAGRA);
  
    magk=MAGK;
    magf=MAGF;
    magOmega=MAGOMEGA;
    magphi=MAGPHI;
    velra=VELRA;
    magenergy=MAGENERGY;
    // for S&N problem
    //#define MAGOMEGA (0.0000447213595499957939281834733746255247088)
    //#define MAGK (0.331755975461247708975210254239793369570)
    //#define MAGF (8.482284786135403683281765211691831E-10)
    //#define MAGPHI (0.009508626344254944156708616148606349)
    //#define VELRA (0.03392913914454161473312706084676733) // from magbondi2.nb
    //#define MAGENERGY (0.003476858935732869614037618612267327)

    // here x is radius and y is density, only depends on dimensionless quantities
#define FYMAG(y,x) ( MAGBETA/(2.0*x*x*x*x*y*y)+MAGTHETA/(gam-1.)*pow(y,gam-1.) -1.0/(x-rgp)+MAGomega*0.5*( pow(x-1.0/x,2.0)/pow(y-1.0,2.0)-x*x)-MAGE )

    // very much like y=1/x^2
    //#define Ymax(x) (100.0/(x*x)) // 100.0 overcompensates
#define LOWESTY (.000001)
#define Ymax(x) ( (x<1.0) ? 1.0/(x*x*x)+0.5 : .9999999) // 1/x^3 is solid upper bound for that region
#define Ymin(x) ( (x<1.0) ? 1.0000001 : LOWESTY)
#define MAGNEAR1ERR (1E-9) // must be smaller than Ymin-1
#define RAFIX 1.000001 // can't hit 1.0 exactly for rho, and even with rho=1 the formulas for vphi/bphi fail(div 0)
    // x range defined by problem statement
    //  End of definition of magnetic bondi problem
    //
    ////////////////////////
#endif

#define DEBUGBONDI 0 // whether to output debug statements in stdout


#if(DEBUGBONDI)
#define TERM1(y,x) (MAGBETA/(2.0*x*x*x*x*y*y))
#define TERM2(y,x) (MAGTHETA/(gam-1.)*pow(y,gam-1.))
#define TERM3(y,x) (-1.0/(x-rgp))
#define TERM4(y,x) (MAGomega*0.5*( pow(x-1.0/x,2.0)/pow(y-1.0,2.0)-x*x))
#define TERM5(y,x) (-MAGE)
#define XTRY 0.826666626
#define YTRY 1.0000001
    printf("%15.10g+%15.10g+%15.10g+%15.10g+%15.10g=%15.10g\n",TERM1(YTRY,XTRY),TERM2(YTRY,XTRY),TERM3(YTRY,XTRY),TERM4(YTRY,XTRY),TERM5(YTRY,XTRY),FYMAG(YTRY,XTRY));
    fflush(stdout);
    //myexit(1);
#endif
    // Now find solutions by finding y which satisifies FY(ep)=0 for each ep for given ndot.


    //  fprintf(stderr,"%15.10g %15.10g\n",ndot,EPS);

    // non self-similar
    if(MAGBONDI||(!wgam53)||(wgam53&&( wpw||( (!wpw)&&(!SELFSIM) ) ) ) ){ // if not chosen self-similar, gam=5/3,non-PW problem
      k=0;
      j=0;
      for(z=1;z<=1+3;z++){ //  1=scalars 2=v_1 3=v_2 4=v_3 , since different positions and want IC to be exactly numerically symmetric
	// going over grid position types
	if(z==1){
	  radiustouse=sca3[1];
	  thetatouse=sca3[2];
	  phitouse=sca3[3];
	}
	else if(z==2){
	  radiustouse=vx3[1];
	  thetatouse=vx3[2];
	  phitouse=vx3[3];
	}
	else if(z==3){
	  radiustouse=vy3[1];
	  thetatouse=vy3[2];
	  phitouse=vy3[3];
	}
	else if(z==4){
	  radiustouse=vz3[1];
	  thetatouse=vz3[2];
	  phitouse=vz3[3];
	}
#if(COORD==3)
	k=j=0;
	LOOPF1
#else
	  LOOPF // may be too expensive
#endif
	  {


#if(MAGBONDI==0)
	    ep=radiustouse[k][j][i];  // assumes ep>rgp!  Make sure!
#else
	    ep=radiustouse[k][j][i]/MAGRA;
#endif
	    if(ep<=rgp){
	      fprintf(fail_file,"ep=%21.15g<=rgp=%21.15g\n",ep,rgp);
	      myexit(1);
	    }
	    // this tells where to start from so first solution hit is right one
	    if(MAGBONDI==0){ // assumes accretion
	      if(ep>EPS) which=1; // starts at lower value
	      else which=-1; // starts at upper value
	    }
	    else{ // assumes wind
	      if(ep>SPFAST) which=1;  // starts at lower value
	      if( (ep<SPFAST)&&(ep>SPALF) ) which=-1; // starts at upper value
	      if( (ep<SPALF)&&(ep>SPSLOW) ) which=1; // starts at lower value
	      if(ep<SPSLOW) which=-1; // starts at upper value
	    }
#if(DEBUGBONDI)
	    printf("starttype: %15.10g\n",which);
#endif
	    // lay down finer mesh for near sonic point so do not step over solutions and cause aliasing problem
	    // this is because there are really 2 solutions, and might accidentally choose other solution if very close to EPS
	    if(MAGBONDI==0){
	      fractioneps=fabs(ep-EPS);
	    }
	    else{
	      // determine closest critical point
	      fractioneps=fabs(ep-SPFAST);
	      ftemp=fabs(ep-SPALF);
	      if(fractioneps>ftemp) fractioneps=ftemp;
	      ftemp=fabs(ep-SPSLOW);
	      if(fractioneps>ftemp) fractioneps=ftemp;
	    }
	    if(fractioneps<.0015) gridsize=GSM3;
	    else if(fractioneps<.015) gridsize=GSM2;
	    else if(fractioneps<.15) gridsize=GSM;
	    else gridsize=GSI;
	    //gridsize=GSM3; //override to see if just resolution problem

#if(MAGBONDI==0)
	    if(!wgam1){ // if gam!=1
	      amm1=Ymin(ep)+(Ymax(ep)-Ymin(ep))/2.0;
	    }
	    else{
	      amm1=Ymin1(ep)+(Ymax1(ep)-Ymin1(ep))/2.0;
	    }
#else
	    amm1=Ymin(ep)+(Ymax(ep)-Ymin(ep))/2.0; // start solution at center of range
#endif	
	    gotgood=0;
	    usealf=0;
#if(MAGBONDI)
	    // need to make sure don't hit al=1 dead on
	    if(fabs(ep-1.0)<MAGNEAR1ERR){
	      ep=RAFIX;
	    }
#endif
	    if(gotgood==0){
	      for(m=1;m<MAXZOOMS;m++){
		am=0;
#if(MAGBONDI==0)
		if(!wgam1){ // if gam!=1
		  da=(pow(2.,m-1)*(Ymax(ep)-Ymin(ep))/pow(gridsize,m));
		}
		else{
		  da=(pow(2.,m-1)*(Ymax1(ep)-Ymin1(ep))/pow(gridsize,m));
		}
#else
		da=pow(2.,m-1)*(Ymax(ep)-Ymin(ep))/pow(gridsize,m);
#endif
		fold=1E30;
		gotit=0;
#if(DEBUGBONDI)
		if(i==112)
		  printf("m: %d gridsize: %d da: %15.10g Ymax(%15.10g): %15.10g\n",m,gridsize,da,ep,Ymax(ep));
#endif
		// use which to control whether to start high or low so get correct solution first
		for(l=0;l<gridsize;l++){
		  al=amm1+which*da*(l-gridsize/2);
		  if(al<LOWESTY) al=LOWESTY; // avoid singularity at y=0 or occasional negative near 0
#if(MAGBONDI==0)
		  if(!wgam1){ // if gam!=1
		    fnew=FY(ndot,al,ep);
		  }
		  else fnew=FY1(ndot,al,ep); // if gam=1
#else
		  if((usealf==0)||( (usealf==1)&&(ep!=1.0) ) ) fnew=FYMAG(al,ep);
		  else  fnew=0.0;
#endif
#if(DEBUGBONDI)
		  if(i==112)
		    printf("i: %d m: %d l: %d ep: %15.10g al: %15.10g newerr: %15.10g olderr: %15.10g\n",i,m,l,ep,al,fabs(fnew),fabs(fold));
#endif

		  if(fabs(fnew)<fabs(fold)){
		    fold=fnew;
		    am=al; // new choice for next zoom center so far
		  }
		  else{ // if error good and travelling back up hill then keep that last point
		    // since once begins to grow again we have past point, assume at point and last am saved is golden
		    if(0&&(MAGBONDI==1)){
		      if(fabs(fold)<1E-2){ // attempt at error check if not monotonic, but always monotonic error so don't use
			amm1=am;
			gotit=1;
			break;
		      }
		    }
		    else{
		      amm1=am;
		      gotit=1;
		      break;
		    }
		  }
		} // end over changing grid size
		// check if we've yet gotten our point.
		if(gotit==0){
		  fprintf(fail_file,"never found solution to bondi problem: m: %d i: %d\n",m,i);
		  myexit(1);
		}
		else{
#if(MAGBONDI==0)
		  if(!wgam1){// if gam!=1
		    ftemp=fabs(FY(ndot,amm1,ep));
		  }
		  else ftemp=fabs(FY1(ndot,amm1,ep)); // if gam=1
#else
		  if((usealf==0)||( (usealf==1)&&(ep!=1.0) ) ) ftemp=fabs(FYMAG(amm1,ep));
		  else  ftemp=0.0;
#endif
		  if(ftemp<(1E-13)){ // limit of accuracy on value of velocity or density if magbondi
		    // have solution so save it(use this even with magbondi for storing a/b grid solutions)
#if(MAGBONDI==0)
		    if(z==1){ V1grid[k][j][i]=-amm1; }
		    else if(z==2){ V2grid[k][j][i]=-amm1; }
		    else if(z==3){ V3grid[k][j][i]=-amm1; }
		    else if(z==4){ V4grid[k][j][i]=-amm1; }
#else
		    if(z==1){ V1grid[k][j][i]=amm1; }
		    else if(z==2){ V2grid[k][j][i]=amm1; }
		    else if(z==3){ V3grid[k][j][i]=amm1; }
		    else if(z==4){ V4grid[k][j][i]=amm1; }
#endif
		    gotgood=1;
		    break; // no need to zoom anymore
		  }
		}
	      }
	    }
	    // should have vanal[1][1][k][j][i] containing velocity solution
	    if(gotgood==0){
	      fprintf(fail_file,"oops, reached end of ZOOM loop! i: %d z: %d ep: %15.10g\n",i,z,ep);
#if(MAGBONDI==0)
	      if(!wgam1){// if gam!=1
		ftemp=fabs(FY(ndot,amm1,ep));
	      }
	      else ftemp=fabs(FY1(ndot,amm1,ep)); // if gam=1
#else
	      if((usealf==0)||( (usealf==1)&&(ep!=1.0) ) ) ftemp=fabs(FYMAG(amm1,ep));
	      else  ftemp=0.0;
#endif
	      fprintf(fail_file,"error down to: %15.10g\n",ftemp);
	      fprintf(fail_file,"Will go ahead and use value\n");
#if(MAGBONDI==0)
	      if(z==1){ V1grid[k][j][i]=-amm1; }
	      else if(z==2){ V2grid[k][j][i]=-amm1; }
	      else if(z==3){ V3grid[k][j][i]=-amm1; }
	      else if(z==4){ V4grid[k][j][i]=-amm1; }
#else
	      if(z==1){ V1grid[k][j][i]=amm1; }
	      else if(z==2){ V2grid[k][j][i]=amm1; }
	      else if(z==3){ V3grid[k][j][i]=amm1; }
	      else if(z==4){ V4grid[k][j][i]=amm1; }
#endif
	    }
	
	  }//end over position
      }// end over each grid

      // copy over if COORD==3 (so not too expensive)
      if(COORD==3){
	LOOPF{
	  V1grid[k][j][i]=V1grid[0][0][i];
	  V2grid[k][j][i]=V2grid[0][0][i];
	  V3grid[k][j][i]=V3grid[0][0][i];
	  V4grid[k][j][i]=V4grid[0][0][i];
	}
      }
	
      
#if(MAGBONDI==0)
      // setup other velocity components and B field
      LOOPF{
	if(COORD==3){
	  vanal[1][1][k][j][i]=V2grid[k][j][i];
	  vanal[1][2][k][j][i]=0.0;
	  vanal[1][3][k][j][i]=0.0;
	}
	else if(COORD==1){// assumes only v_r is nonzero
	  vanal[1][1][k][j][i]=V2grid[k][j][i]*sin(vx3[2][k][j][i])*cos(vx3[3][k][j][i]);
	  vanal[1][2][k][j][i]=V3grid[k][j][i]*sin(vy3[2][k][j][i])*sin(vy3[3][k][j][i]);
	  vanal[1][3][k][j][i]=V4grid[k][j][i]*cos(vz3[2][k][j][i]);
	}
	vanal[2][1][k][j][i]=0.0;
	vanal[2][2][k][j][i]=0.0;
	vanal[2][3][k][j][i]=0.0;
      }
    
      // setup scalars
      LOOPF{
	sanal[1][k][j][i]=S(ndot,V1grid[k][j][i],sca3[1][k][j][i]); // Vb defined from above only on k=j=0
	if(wgam) sanal[2][k][j][i]=KAPPA*pow(sanal[1][k][j][i],gam)/(gam-1.);
	else sanal[2][k][j][i]=alpha*cs*cs*sanal[1][k][j][i];
	sanal[3][k][j][i]=-GM/(sca3[1][k][j][i]-rgp);
      }
#endif
    }
    else{// if gam=5/3 and rgp=0(newtonian pot) and choosing self-similar case
      // for gam=5/3 SS,  the CS_RSTAR and V_RSTAR must be nonzero and defined as a ratio for self-similar case
#if(MAGBONDI==0)

      ftemp1 = -sqrt(2.0*GM*MACH*MACH/(3.0+MACH*MACH)); // (A) coef to v_r set by length RSTAR
      ftemp2 = sqrt(2.0*GM/(3.0+MACH*MACH)); // (B) coef to c_s

      LOOPF{
	if(COORD==3){
	  ftemp=ftemp1*pow((vx3[1][k][j][i]/RATIO),-.5); // vr
	  vanal[1][1][k][j][i]=ftemp;
	}
	else if(COORD==1){
	  ftemp=ftemp1*pow((vx3[1][k][j][i]/RATIO),-.5); // vr
	  vanal[1][1][k][j][i]=ftemp*sin(vx3[2][k][j][i])*cos(vx3[3][k][j][i]);
	  ftemp=ftemp1*pow((vy3[1][k][j][i]/RATIO),-.5); // vr
	  vanal[1][2][k][j][i]=ftemp*sin(vy3[2][k][j][i])*sin(vy3[3][k][j][i]);
	  ftemp=ftemp1*pow((vz3[1][k][j][i]/RATIO),-.5); // vr
	  vanal[1][3][k][j][i]=ftemp*cos(vz3[2][k][j][i]);
	}

	cs = ftemp2*pow((sca3[1][k][j][i]/RATIO),-.5);// or X(sanal[1][k][j][i]) if sanal[1] defined first
	V1grid[k][j][i]=ftemp1*pow((sca3[1][k][j][i]/RATIO),-.5); // define Vb for all zones for below
	V2grid[k][j][i]=ftemp1*pow((vx3[1][k][j][i]/RATIO),-.5); // define Vb for all zones for below
	V3grid[k][j][i]=ftemp1*pow((vy3[1][k][j][i]/RATIO),-.5); // define Vb for all zones for below
	V4grid[k][j][i]=ftemp1*pow((vz3[1][k][j][i]/RATIO),-.5); // define Vb for all zones for below
	// either form should be exactly the same! (good check of code)
#define FORMGEN 1

#if(FORMGEN)
	sanal[1][k][j][i]=S(ndot,V1grid[k][j][i],(sca3[1][k][j][i]/RATIO));
	sanal[2][k][j][i]=KAPPA*pow(sanal[1][k][j][i],gam)/(gam-1.);
#else
	sanal[1][k][j][i]=cs*cs*cs/(pow(gam*KAPPA,1.5));
	sanal[2][k][j][i]=cs*cs*sanal[1][k][j][i]/(gam*(gam-1.));
#endif
	sanal[3][k][j][i]=-GM/((sca3[1][k][j][i]/RATIO)); // rgp should be 0
      }
#endif
    }



#if(MAGBONDI==1)

    ndot=4.0*M_PI*MAGRA*MAGRA*RHORA*VELRA;

    // position is already assumed to be correct scaling in MAGRA

    // copy temp holder of solution into density:
    LOOPF{
      sanal[1][k][j][i]=RHORA*V1grid[k][j][i];
    }
  
#define FIXEDRADIUS(ep) ((fabs(ep-MAGRA)<MAGNEAR1ERR) ? RAFIX*MAGRA : ep)

    LOOPF{
      sanal[2][k][j][i]=magk*pow(sanal[1][k][j][i],gam)/(gam-1.0); // P=K\rho^{\gamma} P=(\gamma-1) u
      sanal[3][k][j][i]=-GM/(sca3[1][k][j][i]-rgp); // \phi = -GM/(r-r_{gp})

      if(COORD==3){
	ep=FIXEDRADIUS(vx3[1][k][j][i]);
	vanal[2][1][k][j][i]=magphi/(ep*ep); // B_{r}=\Phi/r^2
	vanal[1][1][k][j][i]=magf/(RHORA*V2grid[0][0][i]*ep*ep); // vr = f/(rho*r^2)
      
	//    ftemp=epb*(magphi*magphi*rho-4.0*magf*magf*M_PI);
	ep=FIXEDRADIUS(vz3[1][k][j][i]);
	ftemp=ep*(magphi*magphi*(RHORA*V4grid[k][j][i])-magf*magf);
	//    vanal[1][3][k][j][i]=magOmega*(magphi*magphi*epb*epb*sanal[1][k][j][i]-4.0*magf*magf*M_PI*MAGRA*MAGRA)/ftemp;
	vanal[1][3][k][j][i]=magOmega*(magphi*magphi*ep*ep*(RHORA*V4grid[k][j][i])-magf*magf*MAGRA*MAGRA)/ftemp;
	//    vanal[2][3][k][j][i]=4.0*magf*magOmega*magphi*M_PI*(epb-MAGRA)*(epb+MAGRA)*sanal[1][k][j][i]/ftemp;
	vanal[2][3][k][j][i]=magf*magOmega*magphi*(ep-MAGRA)*(ep+MAGRA)*(RHORA*V4grid[k][j][i])/ftemp;
            
	vanal[1][2][k][j][i]=0.0;
	vanal[2][2][k][j][i]=0.0;
      }
      else if(COORD==1){
	// for vx
	ep=FIXEDRADIUS(vx3[1][k][j][i]);
	br=magphi/(ep*ep); // B_{r}=\Phi/r^2
	vr=magf/(RHORA*V2grid[0][0][i]*ep*ep); // vr = f/(rho*r^2)
	//    ftemp=epb*(magphi*magphi*rho-4.0*magf*magf*M_PI);
	ftemp=ep*(magphi*magphi*(RHORA*V2grid[k][j][i])-magf*magf);
	//    vanal[1][3][k][j][i]=magOmega*(magphi*magphi*epb*epb*sanal[1][k][j][i]-4.0*magf*magf*M_PI*MAGRA*MAGRA)/ftemp;
	vphi=magOmega*(magphi*magphi*ep*ep*(RHORA*V2grid[k][j][i])-magf*magf*MAGRA*MAGRA)/ftemp;
	//vanal[2][3][k][j][i]=4.0*magf*magOmega*magphi*M_PI*(epb-MAGRA)*(epb+MAGRA)*sanal[1][k][j][i]/ftemp;
	bphi=magf*magOmega*magphi*(ep-MAGRA)*(ep+MAGRA)*(RHORA*V2grid[k][j][i])/ftemp;
	vtheta=0;
	btheta=0;
      
	vanal[1][1][k][j][i]=vr*sin(vx3[2][k][j][i])*cos(vx3[3][k][j][i])+vtheta*cos(vx3[2][k][j][i])*cos(vx3[3][k][j][i])-vphi*sin(vx3[3][k][j][i]);
	vanal[2][1][k][j][i]=br*sin(vx3[2][k][j][i])*cos(vx3[3][k][j][i])+btheta*cos(vx3[2][k][j][i])*cos(vx3[3][k][j][i])-bphi*sin(vx3[3][k][j][i]);

	// for vy
	ep=FIXEDRADIUS(vy3[1][k][j][i]);
	br=magphi/(ep*ep); // B_{r}=\Phi/r^2
	vr=magf/(RHORA*V3grid[0][0][i]*ep*ep); // vr = f/(rho*r^2)
	//    ftemp=epb*(magphi*magphi*rho-4.0*magf*magf*M_PI);
	ftemp=ep*(magphi*magphi*(RHORA*V3grid[k][j][i])-magf*magf);
	//    vanal[1][3][k][j][i]=magOmega*(magphi*magphi*epb*epb*sanal[1][k][j][i]-4.0*magf*magf*M_PI*MAGRA*MAGRA)/ftemp;
	vphi=magOmega*(magphi*magphi*ep*ep*(RHORA*V3grid[k][j][i])-magf*magf*MAGRA*MAGRA)/ftemp;
	//vanal[2][3][k][j][i]=4.0*magf*magOmega*magphi*M_PI*(epb-MAGRA)*(epb+MAGRA)*sanal[1][k][j][i]/ftemp;
	bphi=magf*magOmega*magphi*(ep-MAGRA)*(ep+MAGRA)*(RHORA*V3grid[k][j][i])/ftemp;
	vtheta=0;
	btheta=0;

	vanal[1][2][k][j][i]=vr*sin(vy3[2][k][j][i])*sin(vy3[3][k][j][i])+vtheta*cos(vy3[2][k][j][i])*sin(vy3[3][k][j][i])+vphi*cos(vy3[3][k][j][i]);
	vanal[2][2][k][j][i]=br*sin(vy3[2][k][j][i])*sin(vy3[3][k][j][i])+btheta*cos(vy3[2][k][j][i])*sin(vy3[3][k][j][i])+bphi*cos(vy3[3][k][j][i]);

	// for vz
	ep=FIXEDRADIUS(vz3[1][k][j][i]);
	br=magphi/(ep*ep); // B_{r}=\Phi/r^2
	vr=magf/(RHORA*V4grid[0][0][i]*ep*ep); // vr = f/(rho*r^2)
	//    ftemp=epb*(magphi*magphi*rho-4.0*magf*magf*M_PI);
	ftemp=ep*(magphi*magphi*(RHORA*V4grid[k][j][i])-magf*magf);
	//    vanal[1][3][k][j][i]=magOmega*(magphi*magphi*epb*epb*sanal[1][k][j][i]-4.0*magf*magf*M_PI*MAGRA*MAGRA)/ftemp;
	vphi=magOmega*(magphi*magphi*ep*ep*(RHORA*V4grid[k][j][i])-magf*magf*MAGRA*MAGRA)/ftemp;
	//vanal[2][3][k][j][i]=4.0*magf*magOmega*magphi*M_PI*(epb-MAGRA)*(epb+MAGRA)*sanal[1][k][j][i]/ftemp;
	bphi=magf*magOmega*magphi*(ep-MAGRA)*(ep+MAGRA)*(RHORA*V4grid[k][j][i])/ftemp;
	vtheta=0;
	btheta=0;

	vanal[1][3][k][j][i]=vr*cos(vz3[2][k][j][i])-vtheta*sin(vz3[2][k][j][i]);
	vanal[2][3][k][j][i]=br*cos(vz3[2][k][j][i])-btheta*sin(vz3[2][k][j][i]);
      }
    }
#endif



    //////////////////
    // now for some IC diagnostics
    if((myid<=0)&&(firstsolve==1)){
      fprintf(analyticout2,"#%21s %21s %21s %21s %21s %21s %21s %21s %21s %21s\n","r","rho","ie","pot","v_r","v_theta","v_phi","b_r","b_theta","b_phi");
#if(MAGBONDI)
      for(i=INFULL1;i<OUTFULL1;i++){
	fprintf(analyticout2," %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",x[2][1][i],sanal[1][0][0][i],sanal[2][0][0][i],sanal[3][0][0][i],V1grid[0][0][i],0.0,vanal[1][3][0][0][i],MAGPHI/(x[2][1][i]*x[2][1][i]),0.0,vanal[2][3][0][0][i]);
      }
#else
#if(COORD==3)
      k=j=0;
      for(i=INFULL1;i<OUTFULL1;i++)
#elif(COORD==1)
	LOOPF
#endif
	{
	  fprintf(analyticout2," %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",sca3[1][k][j][i],sanal[1][k][j][i],sanal[2][k][j][i],sanal[3][k][j][i],V1grid[k][j][i],0.0,0.0,0.0,0.0,0.0);
	}
#endif
      fflush(analyticout2);
      fclose(analyticout2);
    }
    massdot=ndot/(2.0*M_PI);


    if(firstsolve==1){
      firstsolve=0;
    
      // calculate and print some interesting quantities for this analytic solution

      // these calculations only are good for COORD==3
      if(COORD==3){
	k=0;
	j=0;

	// compute time it takes for fluid element to fall from outer edge to inner edge
	time2fall=0;
	for(i=0;i<N1;i++){
	  dep=dx[1][1][i];
      
	  //      time2fall-=2.*dep/(vanal[1][1][k][j][i]+vanal[1][1][k][j][i+1]);
	  time2fall-=dep/V1grid[k][j][i];
	}
	// sound crossing time within moving fluid, but of course fluid is moving, so time to go from outer to inner is
	cscrosstime=0;
	for(i=0;i<N1;i++){
	  dep=dx[1][1][i];
	  // dep/cs
	  cscrosstime+=dep/(pow(sanal[1][k][j][i],0.5*(gam-1.))); // same as -= with -cs
	}
	// times for sound wave travel
	downstream=0;
	for(i=0;i<N1;i++){
	  dep=dx[1][1][i];
      
	  // dep/(v-cs)
	  if(wgam) cs=gam*(gam-1.)*sanal[2][k][j][i]/sanal[1][k][j][i];
	  downstream-=dep/(V1grid[k][j][i]-cs); // -= since going from outer to inner
	}
	upstream=0;
	for(i=0;i<N1;i++){
	  // dep/(v+cs)
	  dep=dx[1][1][i];
	  if(wgam) cs=gam*(gam-1.)*sanal[2][k][j][i]/sanal[1][k][j][i];
	  upstream+=dep/(V1grid[k][j][i]+cs ); // += since going from inner to outer
	}

	min1=0;
	kein1=0;
	hin1=0;
	pein1=0;
	k=0;
	i=0;
	if(mycpupos[1]==0){
	  for(j=0;j<N2;j++){
	    tempf=z2e_1(sanal[1],k,j,i)*vanal[1][1][k][j][i]*G2(1,i)*G3(1,i)*DVL(1,2,j);
	    // inner radial flux of mass
	    min1+=-tempf;
	    kein1+=-0.5*vanal[1][1][k][j][i]*vanal[1][1][k][j][i]*tempf;
	    hin1+=-gam*z2e_1(sanal[2],k,j,i)/z2e_1(sanal[1],k,j,i)*tempf;
	    pein1+=-z2e_1(sanal[3],k,j,i)*tempf;
	  }
	  min1a_full=massdot;
	  kein1a_full=-0.5*min1a_full*vanal[1][1][0][0][0]*vanal[1][1][0][0][0];
	  hin1a_full=-gam*z2e_1(sanal[2],0,0,0)/z2e_1(sanal[1],0,0,0)*min1a_full;
	  pein1a_full=-z2e_1(sanal[3],0,0,0)*min1a_full;
	}
	mout1=0;
	keout1=0;
	hout1=0;
	peout1=0;
	k=0;
	i=N1;
	if(mycpupos[1]==ncpux1-1){
	  mout1a_full=-massdot;
	  keout1a_full=0.5*mout1a_full*vanal[1][1][0][0][N1]*vanal[1][1][0][0][N1];
	  hout1a_full=gam*z2e_1(sanal[2],0,0,N1)/z2e_1(sanal[1],0,0,N1)*mout1a_full;
	  peout1a_full=z2e_1(sanal[3],0,0,N1)*mout1a_full;

	  for(j=0;j<N2;j++){
	    tempf=z2e_1(sanal[1],k,j,i)*vanal[1][1][k][j][i]*G2(1,i)*G3(1,i)*DVL(1,2,j);
	    // inner radial flux of mass
	    mout1+=tempf;
	    keout1+=0.5*vanal[1][1][k][j][i]*vanal[1][1][k][j][i]*tempf;
	    hout1+=gam*z2e_1(sanal[2],k,j,i)/z2e_1(sanal[1],k,j,i)*tempf;
	    peout1+=z2e_1(sanal[3],k,j,i)*tempf;
	  }
	}
      }
      else if(COORD==1){
	// could work on some stuff
      }

#if(MAGBONDI==0)
      if( (!wgam53)||(wgam53&&( wpw||( (!wpw)&&(!SELFSIM) ) ) ) ){ // if not chosen self-similar, gam=5/3,non-PW problem
	if(wgam){
	  ftemp1=Ymax(EPS);
	  ftemp2=Xmax(EPS);
	}
	else{
	  ftemp1=Ymax1(EPS);
	  ftemp2=Xmax1(EPS);
	}
      }
      else{
	ftemp1=vanal[1][1][0][0][-2];
	ftemp2=ftemp1/MACH;
      }
#else
      ftemp1=Ymax(MAGRA);
      ftemp2=0;
#endif
    
      if(numprocs>1){
#if(USEMPI)
	MPI_Reduce(&(min1), &(min1_full), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
	MPI_Reduce(&(mout1), &(mout1_full), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
	MPI_Reduce(&(kein1), &(kein1_full), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
	MPI_Reduce(&(keout1), &(keout1_full), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
	MPI_Reduce(&(hin1), &(hin1_full), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
	MPI_Reduce(&(hout1), &(hout1_full), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
	MPI_Reduce(&(pein1), &(pein1_full), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
	MPI_Reduce(&(peout1), &(peout1_full), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
#endif
      }
      else{
	mout1_full=mout1;
	keout1_full=keout1;
	hout1_full=hout1;
	peout1_full=peout1;
	min1_full=min1;
	kein1_full=kein1;
	hin1_full=hin1;
	pein1_full=pein1;
      }
      if(myid<=0){
#if(MAGBONDI)
	fprintf(analyticout,"#%21s %21s %21s %21s %21s %21s %21s %21s\n","Theta","omega","gamma","beta","En","ra","rhora","GM");
	fprintf(analyticout," %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",MAGTHETA,MAGomega,gam,MAGBETA,MAGE,MAGRA,RHORA,GM);
	fprintf(analyticout,"#%21s %21s %21s %21s %21s %21s\n","magk","magf","magOmega","magphi","magenergy","velra");
	fprintf(analyticout," %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",magk,magf,magOmega,magphi,magenergy,velra);
#endif
	if(COORD==3){ // nothing for otherwise yet
	  // inner radial fluxes (slopes as function of t)
	  // sign is + loss through inner edge
	  fprintf(analyticout,"#%21s %21s %21s %21s\n","min1","kein1","hin1","pein1");
	  fprintf(analyticout," %21.15g %21.15g %21.15g %21.15g\n",min1_full,kein1_full,hin1_full,pein1_full);
	  fprintf(analyticout,"#%21s %21s %21s %21s\n","min1a","kein1a","hin1a","pein1a");
	  fprintf(analyticout," %21.15g %21.15g %21.15g %21.15g\n",min1a_full,kein1a_full,hin1a_full,pein1a_full);
	
	  // outer radial fluxes (slopes as function of t)
	  // sign is + loss through outer edge
	  fprintf(analyticout,"#%21s %21s %21s %21s\n","mout1","keout1","hout1","peout1");
	  fprintf(analyticout," %21.15g %21.15g %21.15g %21.15g\n",mout1_full,keout1_full,hout1_full,peout1_full);
	  fprintf(analyticout,"#%21s %21s %21s %21s\n","mout1a","keout1a","hout1a","peout1a");
	  fprintf(analyticout," %21.15g %21.15g %21.15g %21.15g\n",mout1a_full,keout1a_full,hout1a_full,peout1a_full);
      
	  fprintf(analyticout,"#%21s %21s %21s %21s\n","Time to Fall","cs crosstime","downstream time","upstream time");
	  fprintf(analyticout," %21.15g %21.15g %21.15g %21.15g\n",time2fall,cscrosstime,downstream,upstream);
	}

	fprintf(analyticout,"#%21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s\n","GM","TYPEKAPPA","RSTAR","V_RSTAR","CS_RSTAR","RHO_RSTAR","KAPPA","Bern","RATIO","RS(sonic point)","XS(cs at sp)","VS(v at sp)","SS(rho at sp)","ndot(total mrate)","massdot(asym mrate)");
	fprintf(analyticout," %21.15g %21d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",GM,TYPEKAPPA,RSTAR,V_RSTAR,CS_RSTAR,RHO_RSTAR,KAPPA,Bern,RATIO,EPS,XS,VS,SS,ndot,massdot);
	fprintf(analyticout,"#%21s %21s\n","v_MAX(EPS)","cs_MAX(EPS)");
	fprintf(analyticout," %21.15g %21.15g\n",ftemp1,ftemp2);


	if(calltype!=2) fclose(analyticout);
      }// end if write cpu


      //// VISUALIZATION settings


      for(i=0;i<ITYPES;i++){
	for(j=0;j<CTYPES;j++){
	  mms[i][j][1][0]=1;
	  mms[i][j][1][1]=1.61;
	  mms[i][j][2][0]=.9;
	  mms[i][j][2][1]=2.;
	  mms[i][j][3][0]=-1.96;
	  mms[i][j][3][1]=0;
	
	  mmv[i][j][1][0][0]=0;
	  mmv[i][j][1][0][1]=1.65;
	  mmv[i][j][1][1][0]=0;
	  mmv[i][j][1][1][1]=1.65;
	  mmv[i][j][1][2][0]=0;
	  mmv[i][j][1][2][1]=1.65;
	  mmv[i][j][1][3][0]=0;
	  mmv[i][j][1][3][1]=1.65;

	  // define outer region when interpolation is used.
	  // same order as scalar/vector arrays
	  outerdefs[i][j][1]=mms[0][0][1][0]; // rho
	  outerdefs[i][j][2]=mms[0][0][2][0]; // en
	  outerdefs[i][j][3]=mms[0][0][3][0]; // pot
	
	  outerdefv[i][j][1][0]=0.0; // magnitude of v
	  outerdefv[i][j][1][1]=0.0; // v1
	  outerdefv[i][j][1][2]=0.0; // v2
	  outerdefv[i][j][1][3]=0.0; // v3
	
	  outerdefv[i][j][2][0]=0.0;
	  outerdefv[i][j][2][1]=0.0;
	  outerdefv[i][j][2][2]=0.0;
	  outerdefv[i][j][2][3]=0.0;


	}
      }
    
      douterdefs[1]=mms[0][0][1][0]; // rho
      douterdefs[2]=mms[0][0][2][0]; // en
      douterdefs[3]=mms[0][0][3][0]; // pot
    
      douterdefv[1][0]=0.0; // magnitude of v
      douterdefv[1][1]=0.0; // v1
      douterdefv[1][2]=0.0; // v2
      douterdefv[1][3]=0.0; // v3
    
      douterdefv[2][0]=0.0;
      douterdefv[2][1]=0.0;
      douterdefv[2][2]=0.0;
      douterdefv[2][3]=0.0;

    }// end if first solve

  
    free(works1a);
    free(works2a);
    free(works3a);
    free(works4a);
  }// end else if calltype==100
}// end of used method--more general
















FILE * tori1sol(int calltype)
{
  SFTYPE numvtime,tvisc;
  int outtype;
  SFTYPE rho;
  SFTYPE divb1avg,divb2avg;
  SFTYPE divb1max,divb2max;
  SFTYPE divb1avg_full,divb2avg_full;
  SFTYPE divb1max_full,divb2max_full;
  char check[50];
  int i,j,k,l,m,n;
  SFTYPE temp,temp2;
  SFTYPE radsize2,tempf1,tempf2,xz0;
  SFTYPE lkep0,Cconst,deltatheta;
  SFTYPE C2const,Beta,RBeta;
  SFTYPE ftemp,postemp;
  char temps[100];
  char filename[100];
  SFTYPE magbeta;
  SFTYPE b0,torib0,vertb0,pgtot,pbtot,pgtot_full,pbtot_full;
  SFTYPE Rpot0;
  SFTYPE theta0;
  SFTYPE Rout,l0sq;
  SFTYPE CFRACT;
  int gooddensity;
  SFTYPE phi_eff,phi_tot;
  static FILE*analyticout;
  static int firstsolve=1;
  SFTYPE totalmass[2]; // 0: tori 1: atmosphere
  SFTYPE totalmass_full[2]; // 0: tori 1: atmosphere
  // whether tori is at a point in the grid or not
  int torigrid;
  FTYPE (*vectorpot)[N3M][N2M][N1M]; // phi component of vector potential
  static int ATMONLYT, ADATM,MAGTORI,GLOBALFIELD;
  static int EQATM, UNITS,RANDOMIZE,VRANDOMIZE;
  SFTYPE MODEPURT,MODEAMP_rho,MODEAMP_u,MODEAMP_v1,MODEAMP_v2,MODEAMP_v3,MODEAMP_b1,MODEAMP_b2,MODEAMP_b3;
  SFTYPE QCONST;
  SFTYPE randomization;
  SFTYPE vr,vtheta,vphi,br,btheta,bphi;
  FTYPE (*sca3)[N3M][N2M][N1M];
  FTYPE (*vx3)[N3M][N2M][N1M];
  FTYPE (*vy3)[N3M][N2M][N1M];
  FTYPE (*vz3)[N3M][N2M][N1M];

  FTYPE (*radiustouse)[N2M][N1M];
  FTYPE (*thetatouse)[N2M][N1M];
  FTYPE (*phitouse)[N2M][N1M];
  SFTYPE ftemp1,ftemp2;
  int itemp;
  int TORITYPE;
  FTYPE gravpottemp;
  int component;
  FTYPE RGPFUDGE;
  SFTYPE ptemp,xtemp,ytemp,xtemp1,xtemp2,amp,phase,r0n;
  SFTYPE vrand;
  SFTYPE KAPPA;
  SFTYPE SMOOTHSIZE_rho,SMOOTHSIZE_u,SMOOTHSIZE_v1,SMOOTHSIZE_v2,SMOOTHSIZE_v3;
  SFTYPE ATMFACTOR;
  int OMEGATYPE;
  SFTYPE DISKTYPE,Routt;
  SFTYPE shapeouter;


  TORITYPE=0;
  // 0: my type
  // 1: gammie type
  // 2: sp01 type RUNF (or mod)
  // 3: Dmitry Psaltis type
  // 4: Gammie thin disk

  OMEGATYPE=0; // not used yet, if ever due to phi_eff craziness with rgp=2
  // 0: $\Omega\propto R^{-q}$
  // 1: $\Omega\propto R^{-q} (1-r_g/R)^{-1}$

  ATMONLYT= 0; // 0: normal torus 1: just atmosphere
  ATMFACTOR=1.0; // factor above floor given to atmosphere
  ADATM= 0; // 1: adiabatic atm 0: fraction pot atm
  MAGTORI=0; // 1: include magnetic equadensity rings 0: don't have any magnetic field
  // Fill analytic solution with solution to tori solution, mostly used for initial conditions and comparisons as with bondi
  GLOBALFIELD= 0; // 1: set field all over grid 0: don't

  RANDOMIZE=1; // whether any randomization at all
  VRANDOMIZE=2; // 1==vphi only for velocity, and density, 2=all (see below) 3=m=1 density pert

  // for VRANDOMIZE==1 or 2
  //
  vrand=1.0; // scale for randomizing velocities
  randomization=5.E-2; //  (at N3=32)

  // FOR VRANDOMIZE==3
  //
  MODEPURT=4.0; // mode to purtube for VRANDOMIZE==3

  MODEAMP_rho=3.16E-5; // mode amplitude "" as fraction of density maximum
  MODEAMP_u=MODEAMP_rho*0.4/0.2; // mode amplitude "" as fraction of density maximum

  MODEAMP_v1=MODEAMP_rho*.01/0.2; // mode amplitude "" as fraction of density maximum
  MODEAMP_v2=MODEAMP_rho*.01/0.2; // mode amplitude "" as fraction of density maximum
  MODEAMP_v3=MODEAMP_rho*.01/0.2; // mode amplitude "" as fraction of density maximum

  MODEAMP_b1=0.0; // mode amplitude "" as fraction of density maximum
  MODEAMP_b2=0.0; // mode amplitude "" as fraction of density maximum
  MODEAMP_b3=0.0; // mode amplitude "" as fraction of density maximum

  // SMOOTHING
  //
  SMOOTHSIZE_rho=0; // kernel around zone to combine in smoothing, 0=no smoothing
  SMOOTHSIZE_u=0; // kernel around zone to combine in smoothing, 0=no smoothing
  SMOOTHSIZE_v1=0; // kernel around zone to combine in smoothing, 0=no smoothing
  SMOOTHSIZE_v2=0; // kernel around zone to combine in smoothing, 0=no smoothing
  SMOOTHSIZE_v3=0; // kernel around zone to combine in smoothing, 0=no smoothing

  if(calltype==100){
    if(POSTPROC==1){
      DYNAMICMM=0; 
      // already read in file usually
    }
    else{
      //DYNAMICMM=0; // use data below
      DYNAMICMM=2; // testing
    }
    // PW pot
    rg=rgp=2;
    
    // Newtonian pot
    //    rg=2; rgp=0;

      
    if(MAGTORI){
      mag=1;
      visc_real=0;
      res_real=1; rreal=2; resheat= 1; // resistivity
      resist_real0=.01; // from vortex tests
      //res_real= 1;  rreal=2; resheat= 1; resist_real0=0.1; // resistivity
      //res_real=0;
    }
    else{
      visc_real=1;      vreal=4;
      mag=0;      res_real=0;
    }
    if(TORITYPE==0){
      if(COORD==3){
	x3in=  (0.0);
	x3out= (2.0*M_PI);
	
	//Rinner=.01;
	Rinner=1.2*rg;
	Router=11.5*rg; // normal outer boundary HK-like
       	//Router=10.0*11.5*rg; // 10X HK outer boundary
	//Router=100.0*11.5*rg; // 100X HK outer boundary
	//x1in=3.5*rg; // PP insta-1 (A2/A3)
	//x1in=3.0*rg; // PP insta-1 (C2, and now A2/A3 too)
	//x1in=2.5*rg; // PP insta-1 (C2, and now A2/A3 too Run 25)
	//x1in=1.5*rg;    // H00 MHD axisym runs 1-3
	// avoid superfast issue? --nope
	//x1in=1.02*rg; // avoid superfast point issue?
	//x1in=1.05*rg; // avoid superfast issue? -- nope
	//x1in=2.0*rg; // SP00 RunA
	//x1in=(1.34666*rg);
	//x1in= (3.0*rg);
	//x1in=  (1.4*rg);
	//    x1out=25.0*rg; // PP insta-1 (A2, A3 try2 (run25))
	//x1out=15.0*rg; // PP insta-1 (A2, A3 real)
	//x1out=40.0*rg; // PP insta-1 (C2)
	
	x1in=Rinner;
	//x1out= (4); // SPB99 Model A
	x1out= Router; // H00 first axisym
	//x1out= (160.0*rg); // SP00 RunA
	//x1out= (80.0*rg); // SP00 RunA mod
	//x1out= (40.0*rg); // SP00 RunA mod2
	//x1out= (21.5*rg); // H00 third axisym
	//x2in=  (0.0);
	//x2out= (M_PI);
	
	if(N2==1){
	  x2in=  (M_PI*0.5-1E-6); // PP-insta-1 (A2,C2)
	  x2out= (M_PI*0.5+1E-6); // PP-insta-1 (A2,C2)
	}
	else{ // shouldn't be used if dtheta>>1 since then dvl2 much diff than sin(theta)*dx2
	  x2in=  (0.0); // A3
	  x2out= (M_PI);
	}
	// normal
	nonunigridx1= 5; // 0=unigrid 1+=nonunigrid
	nonunigridx2= 3; // 0=unigrid 1+=nonunigrid
	nonunigridx3= 0;

  
	simplebc=1;
	bcix1=4;
	bcox1=4;
	bcix2=2; // really AOS (see boundrect1.c)
	bcox2=2; // really AOS (see boundrect1.c)
	bcix3=5;
	bcox3=5;

	// for bondi overlay (assumes torus is contained on grid
	if(analoutput==14)	bcox1=3;

      }
      else if(COORD==1){ // cart type 
	Rinner=1.5*rg;

	// tight on torus
	/*
	  Router=11.5*rg;
	  x1in=-11.5*rg;
	  x1out=11.5*rg;
	  x2in=-11.5*rg;
	  x2out=11.5*rg;
	
	  x3in=-11.5*rg*0.25;
	  x3out=11.5*rg*0.25;
	*/
	// loose on torus

	/*	
	  Router=11.5*rg;
	  x1in=-11.5*rg*1.5;
	  x1out=11.5*rg*1.5;
	  x2in=-11.5*rg*1.5;
	  x2out=11.5*rg*1.5;
	
	  x3in=-11.5*rg*0.5;
	  x3out=11.5*rg*0.5;
	*/
	// 2D height integrated
	// n_2D=n_3D+1/2, gam=1+1/n (i.e. gam=5/3 -> gam=3/2)
	//gam=3./2.;
	gam=5.0/3.0; // H90 2D PP-inst
	Rinner=80.0;
	Router=120.0;


	x1in=-Router;
	x1out=Router;
	x2in=-Router;
	x2out=Router;
	
	x3in=-1E-10; // with outflow or periodic
	x3out=1E-10;
	


	nonunigridx1= 0; // 0=unigrid 1+=nonunigrid
	nonunigridx2= 0; // 0=unigrid 1+=nonunigrid
	nonunigridx3= 0;
	
	// only applies if BOUNDTYPE==1
	simplebc=1;
	bcix1=4;
	bcox1=4;
	bcix2=4;
	bcox2=4;
	bcix3=4;
	bcox3=4;
	  
      }
      
    }
    else if(TORITYPE==1){
      // override
      nonunigridx1=nonunigridx2=0;
      Rinner=1.2*rg;
      x1in=Rinner;
      x1out=10.0*rg;
      x2in=  (0.0);
      x2out= (M_PI);
      x3in=  (0.0);
      x3out= (2.0*M_PI);
  
      simplebc=1;
      bcix1=4;
      bcox1=4;
      bcix2=2;
      bcox2=2;
      bcix3=5;
      bcox3=5;
    }
    else if(TORITYPE==2){
      nonunigridx1=5;
      nonunigridx2=3;
      nonunigridx3=0;
      
      //      Rinner=2.0*rg; // normal SP01
      Rinner=2.5;
      x1in=Rinner;
      x1out=400.0*rg;	
      x2in=  (0.0); // A3
      x2out= (M_PI);
      x3in=  (0.0);
      x3out= (2.0*M_PI);
  
      simplebc=1;
      bcix1=4;
      bcox1=4;
      bcix2=2;
      bcox2=2;
      bcix3=5;
      bcox3=5;
    }  
    else if(TORITYPE==3){
      nonunigridx1=0;
      nonunigridx2=4;
      nonunigridx3=0;
      
      // isothermal
      gam=1.0;  transiex1=transiex2=transiex3=ie=0;

      Rinner=4.0;
      Router=15.0;
      x1in=Rinner;
      x1out=Router;
      x2in=  (0.0);
      x2out= (M_PI);
      x3in=  (0.0);
      x3out= (2.0*M_PI);
  
      simplebc=1;
      bcix1=4;
      bcox1=4;
      bcix2=2;
      bcox2=2;
      bcix3=5;
      bcox3=5;
  
    }
    else if(TORITYPE==4){
      nonunigridx1=5;
      nonunigridx2=5;
      nonunigridx3=0;
      
      gam=5.0/3.0;

      cool=1;			coolfact=1;

      Rinner=1.2*rg;
      //      Router=10.0*11.5*rg;
      Router=100.0;
      x1in=Rinner;
      x1out=Router;
      x2in=  (0.0);
      x2out= (M_PI);
      x3in=  (0.0);
      x3out= (2.0*M_PI);
  
      simplebc=1;
      bcix1=4;
      bcox1=4;
      bcix2=2;
      bcox2=2;
      bcix3=5;
      bcox3=5;
  
    }
    
  

    if(TORITYPE==0){
      // viscosity run
      tvisc=2.0*M_PI*pow(12.0,1.5)/alpha_real0; //viscous time scales
      numvtime=2.0;
      tf = numvtime*tvisc; // N viscous time scales
      //tf = 2.0*M_PI*4.325;
      //tf = 10.0;
      //  tf = 4000.0;
      //tf = .245;
      //  tf = 0.240869213846978;
      //tf = .279;
      //tf = 50.0*(2.*M_PI);
      //tf     = 3.0+4.4*(2.*M_PI) ; 
      ////tf     = 3.0+2.3*(2.*M_PI) ; // end time
      //tf     = 100.0*(2.*M_PI) ; // end time
      //  timagescale=tvisc;
      //tavgi  = 3.0+2.0*(2.0*M_PI); // time to start averaging
      //tavgf  = tavgi+0.278*(2.0*M_PI); // time to stop averaging(only approximate)
      //DTl    = (1./1000.)*(2.*M_PI) ;       /* logging period(lower limit on all except DTtimestep, this is also the ener/loss file output period */
      
      // mhd tori up close to rg
      //  DTl    = (1.0/30.0);
      //  DTl    = (1.0/5.0);
      //DTl = DTCONS;
      
      // MHD TORI
      //    tf = 2423.4; // R0=4.7rg 17 orbits 
      //tf=10.0*6000.0; // PP-insta-1 (A2,C2)
      //tf=3500.0; // PP-insta-1 (A2,C2)
      //tf=3195.0; // PP-insta-1 (A3)
      //tf=1500.0; // boring after this, HK00
      //tf=2500.0; // b2mdot-1
      //tf=4000.0; // jet
      //    tf = 3000.0;
      // GODMARK(tf) -- commented
      //    tf = 10.0;
      //tf=5.08*2.0*M_PI*pow(40.0*rg,1.5); // SP00 RunA
      //tf=5.08*2.0*M_PI*pow(20.0*rg,1.5); // SP00 RunA mod
      //tf=5.08*2.0*M_PI*pow(10.0*rg,1.5); // SP00 RunA mod2
      //tf = 1000.0; // R0=4.7rg 5 orbits
      //tf = 5; // test
      // SP00
      //nonunigridx1= 3; // 0=unigrid 1+=nonunigrid
      //nonunigridx2= 3; // 0=unigrid 1+=nonunigrid
      // H90 2D PP-inst
      //      tf=43960.0*2.0; // 14 orbits at R0=100
      // H87 Model 3 2d pp inst
      //      tf=2.*M_PI*pow(100.,1.5)*4.5; // 4.5 orbits at R0
    }
    else if(TORITYPE==1){
      //tf=2000.0;
      tf=2000.0;
    }
    else if(TORITYPE==2){
      tf=4.5*2.0*M_PI*pow(100.0*rg,1.5); // SP01 RunF
    }
    else if(TORITYPE==3){
      tf=2500.0;
    }
    else if(TORITYPE==4){
      tf=2.0*M_PI*pow(Router,1.5)*1.5; // 50% longer than cooling for outer edge
    }


    if(TORITYPE==0){
      //DTl=tf/15000.0; // visc
      //      DTl=tf/50000.0; // H90 2D PP-inst
      //      DTl=tf/50000.0; // H87 Model3
      //    DTl = 1.0;
      //DTl = 1.0;
      DTl = 1.0/30.0; // H00
      //DTl = 1.0/5.0; // 3d cart H00-like run
      //DTl = 1.0/5.0; // SP00 runA
      //DTl = 10.0; // SP00 runA
      // strictly for purposes of reentrant ability, should have dump/floor same DT
      DTd    = DTl*250.0 ;      /* dumping period(data and analytic)*/
      DTfloor=DTd*10.0 ;      /* dumping period(data and analytic)*/
      // strictly for reentract ability, DTi should be integral number of DTds.
      DTpd = 5.0*DTd; // always reentrant data
      DTi    = DTl*30.0 ;        /* image & fieldline(Ax3) period */
      DTfld=DTi*2.0;
      //below not restricted by DTl calls, just each self
      DTener = DTl*1.0 ; // negative means do every time step (must be multiple of DTl)
      // f2's FFT shows could do 20.0 here(was 4.0)
      DTloss = DTl*1.0 ; // negative means do every time step (must be multiple of DTl
      DTmode = DTl*20.0;
      
      DTtimestep = DTl*100.0 ; // how often to output dominate term in timestep to failfile
      DTdivb=DTsp = DTl*100.0 ; // how often to output sonic point info
    }
    else if(TORITYPE==1){
      DTl=1.0/4.0;
      //DTl=10.0;
      DTener=DTloss=DTl;
      DTi=DTfld=4.0*DTl;
      DTd=DTfloor=DTpd=40.0*DTl;
      DTdivb=DTtimestep=DTsp=100.0*DTl;
      DTmode=1000.0*DTl;
    }
    else if(TORITYPE==2){
      DTl=4.0;
      //DTl=10.0;
      DTener=DTloss=DTl;
      DTi=DTfld=4.0*DTl;
      DTd=DTfloor=DTpd=16.0*DTl;
      DTdivb=DTtimestep=DTsp=100.0*DTl;
      DTmode=1000.0*DTl;
    }
    else if(TORITYPE==3){
      DTl=93.8/100.0;
      DTener=DTloss=DTl;
      DTi=DTfld=DTl*4.0;
      DTd=DTfloor=DTpd=DTl*10.0;
      DTdivb=DTtimestep=DTsp=100.0*DTl;
      DTmode=1000.0*DTl;
    }
    else if(TORITYPE==4){
      DTl = 1.0/30.0; // H00
      DTd    = DTl*250.0 ;      /* dumping period(data and analytic)*/
      DTfloor=DTd*10.0 ;      /* dumping period(data and analytic)*/
      DTpd = 5.0*DTd; // always reentrant data
      DTi    = DTl*30.0 ;        /* image & fieldline(Ax3) period */
      DTfld=DTi*2.0;
      DTener = DTl*1.0 ; // negative means do every time step (must be multiple of DTl)
      DTloss = DTl*1.0 ; // negative means do every time step (must be multiple of DTl
      DTmode = DTl*20.0;
      DTtimestep = DTl*100.0 ; // how often to output dominate term in timestep to failfile
      DTdivb=DTsp = DTl*100.0 ; // how often to output sonic point info
    }
    if(TORITYPE==0){
    }
    else if(TORITYPE==3){
			
			
      HOR=1.0; // H/R
      RMAX=2.0*(2.0+sqrt(3.0));
      DELTAR=0.1*RMAX;
    }
    else if(TORITYPE==4){
      HOR=0.1;
    }
  }// everything after this is called only after grid is init'ed, etc.
  else if( (calltype==0)||(calltype==2)){ // normal or injection called it

    vectorpot=workv1;

    sca3=workv2;
    vx3=workv3;
    vy3=workv4;
    vz3=workv5;

    grids_cart2spc(sca3,vx3,vy3,vz3);

    /*
    // GODMARK -- commented
    // symmetry check
    LOOPF{
    // sca3
    //s[1][k][j][i]=sca3[1][k][j][i]; itemp=0;
    //      s[1][k][j][i]=fabs(M_PI/2.0-sca3[2][k][j][i]); itemp=0; // make it z-sym, already rot sym
    //s[1][k][j][i]=fabs(sin(sca3[3][k][j][i])); itemp=0; // make it rot sym, already z-sym

    // vx (still a scalar, but in vx position)
    //s[1][k][j][i]=vx3[1][k][j][i]; itemp=1; // good
    //s[1][k][j][i]=fabs(M_PI/2.0-vx3[2][k][j][i]); itemp=1; // good (10^{-12} variance)
    //s[1][k][j][i]=fabs(sin(vx3[3][k][j][i])); itemp=1; // good (10^{-12} variance)

    // vy (still a scalar, but in vy position)
    //s[1][k][j][i]=vy3[1][k][j][i]; itemp=2; // good
    //s[1][k][j][i]=fabs(M_PI/2.0-vy3[2][k][j][i]); itemp=2; // good
    //s[1][k][j][i]=fabs(sin(vy3[3][k][j][i])); itemp=2; // good

    // vz (still a scalar, but in vz position)
    //s[1][k][j][i]=vz3[1][k][j][i]; itemp=3; // good
    //s[1][k][j][i]=fabs(M_PI/2.0-vz3[2][k][j][i]); itemp=3; // good
    //      s[1][k][j][i]=fabs(sin(vz3[3][k][j][i])); itemp=3; // good

    // check coordinates
    //s[1][k][j][i]=fabs(x[2][1][i]); itemp=0; // good
    //s[1][k][j][i]=fabs(x[1][1][i]); itemp=1; // good
    //s[1][k][j][i]=fabs(x[2][2][j]); itemp=0; // good
    //s[1][k][j][i]=fabs(x[1][2][j]); itemp=2; // good
    //s[1][k][j][i]=fabs(x[2][3][k]); itemp=0; // good
    //s[1][k][j][i]=fabs(x[1][3][k]); itemp=3; // good
    }
    symmetry_check(itemp);
    myexit(0);
    */


    // GODMARK(commented)
    /*
      i=j=32;
      LOOPF3{
      fprintf(stdout,"k=%d vz3[1]=%15.10g vz3[2]=%15.10g vz3[3]=%15.10g x: %15.10g y: %15.10g z: %15.10g\n",k,vz3[1][k][j][i],vz3[2][k][j][i],vz3[3][k][j][i],x[2][1][i],x[2][2][j],x[2][3][k]);
      }
      fflush(stdout);
      myexit(0);
    */
    

    // GODMARK(commented)
    
    /*
      j=32;
      k=0;
      LOOPF1{
      fprintf(stdout,"i=%d vx3[3]=%11.10g vy3[3]=%11.10g vz3[3]=%11.10g xa: %11.10g xb: %11.10g ya: %11.10g yb: %11.10g za: %11.10g zb: %11.10g\n",i,vx3[3][k][j][i],vy3[3][k][j][i],vz3[3][k][j][i],x[1][1][i],x[2][1][i],x[1][2][j],x[2][2][j],x[1][3][k],x[2][3][k]);
      }
      fflush(stdout);
      myexit(0);
    */    

    
    
    // want to write some interesting data on solution
    if(firstsolve==1){
      if(myid<=0){
	fprintf(logfull_file,"torisol\n");
	fflush(logfull_file);
	strcpy(temps,DATADIR);
	
	sprintf(filename,"%s0_analdata%s",temps,DAT2EXT);
	if((analyticout=fopen(filename,"wt"))==NULL){
	  fprintf(fail_file,"Cannot open %s\n",filename);
	  myexit(1);
	}
      }//end if cpu write
    }
    
    
    // setup tori conditions
    // MARK
    if(TORITYPE==0){
      magbeta=100.0; // controls absolute strength of magnetic field
      //magbeta=500.0; // test to see if magbeta changes b^2/mdot
    }
    else if(TORITYPE==1){
      magbeta=400.0;
      // magbeta=2.*20.*4.E2*(3.5365); // includes correction for bug in gammie init.c for test4 case only
      //magbeta=900.0; // includes correction for bug in gammie init.c for test4 case only
      // magbeta=9000.0; // includes correction for bug in gammie init.c for test4 case only
    }
    else if(TORITYPE==2){
      magbeta=200.0;
    }
    else if(TORITYPE==3){
      magbeta=100.0;
    }
    else if(TORITYPE==4){
      magbeta=100.0;
    }


    torib0=1.0; // controls relative strength to any other fields added
    vertb0=torib0/10.0; // controls relative strength to any other fields added
    EQATM=0;

    theta0=M_PI/2.0;
    if(TORITYPE==0){
      //QCONST=1.5;
      QCONST= 2.0; // H00 first axisym
      //    QCONST= 1.88; // H00 first axisym //GODMARK (QCONST)
      //QCONST=1.68; // H00 third axisym
      //      QCONST=1.7; // thinner
    }
    else if(TORITYPE==1){
      QCONST= 2.0; // H00 first axisym
      //    QCONST= 1.88; // H00 first axisym //GODMARK (QCONST)
      //QCONST=1.68; // H00 third axisym
    }
    else if(TORITYPE==2){
      QCONST=2.0;
    }
    else if(TORITYPE==3){
      QCONST=1.5;
    }
    else if(TORITYPE==4){
      QCONST=1.5;
    }

    UNITS= 2;
    // 1: original torus units(R0=1,Omega0=1,M=R0^3*rho0)
    // 2: GM=rg/2=1
    // 3: arbitrary


    RGPFUDGE=1.01; // used INSIDE grav pot when r<=rgp.  Should never be used, but good for plots (should always be inside any real domain value)
    
    if(UNITS==1){
      R0=1.0; // must be 1.0
      Omega0=1.0; // must be 1.0
      rho0=1.0; // must be 1.0
      GM=R0*R0*R0*Omega0*Omega0*pow((1./sin(theta0)-rgp/R0),2.0);
      KAPPA=1.0;
    }
    else if(UNITS==2){
      if(fabs(rg-2.0)>ERR){
	fprintf(fail_file,"These units require rgp to be set to 2.0 in init.c!\n");
	myexit(1);
      }
      //#define MDOTIX1 ((0.8)*(L[1][1]+L[2][1]/4.0))
      //#define MDOTOX1 ((0.85)*(L[1][1]+L[2][1]/4.0))
      rho0=1.0;
      GM=1.0;

      if(TORITYPE==0){
	//    rg=2.0; // make sure rg is 2.0 in init.c!!
	// number of rgp/2's (400 here with units==2 is like rgp=.005 with units==1)
	//R0=(1.0)*(L[1][1]+L[2][1]/4.0);
	//R0=(.825)*(L[1][1]+L[2][1]);
	R0 = 4.7*rg; // H00 first axisym setup
	//R0 = 2.0*rg; // H00 first axisym setup // GODMARK
	//R0 = 7.65135*rg; // PP-insta-1 (A2/A3/C2) (since l0=R0^2*Omega0=4.5) // hawleyconv.nb
	//R0 = 40.0*rg; // SP00 RunA
	//R0 = 20.0*rg; // SP00 RunA mod
	//R0 = 10.0*rg; // SP00 RunA mod2
	//    R0 = 10.0*rg; // H00 third axisym setup
	//    R0=(.6)*(L[1][1]+L[2][1]);
	//R0=166.0*rg;
	//R0=50.*rg; // H87 PP Model3
      }
      else if(TORITYPE==1){
	l0sq=pow(5.0,2.0); // l0=5.0;
	r0n = l0sq/(GM);
	ptemp = rgp/(r0n/3.0);
	ytemp = ptemp*ptemp*ptemp*(4.0/9.0-ptemp);
	//xtemp = 1. - 3.*ptemp + 3./2.*ptemp*ptemp + 3./2.*sqrt(y);
	xtemp1 = 1. - 3.*ptemp + 3./2.*ptemp*ptemp; // real part of xtemp (x=Acos(phase))
	xtemp2 = 3./2.*sqrt(ytemp); // Imaginary part of xtemp (y=Asin(phase))
	//R0 = r0n/3.*(1. + (1. - 2.*ptemp)/pow(x,1./3.) + pow(x,1./3.));
	amp=sqrt(xtemp1*xtemp1+xtemp2*xtemp2);
	phase=atan2(xtemp2,xtemp1);
	R0 = r0n/3.*(1.+(1.-2.*ptemp)*pow(amp,-1./3.)*cos(-phase/3.)+pow(amp,1./3.)*cos(phase/3.));

	// try closer to GR solution of R0
	R0 = 17.0;
	// closer in shape (pointy)
	//R0 = 10.7;
	// trying to get cusp
	//	R0 = 6.5; // generates small thing

      }
      else if(TORITYPE==2){

	R0 = 100.0*rg;
      }
      else if(TORITYPE==3){
	R0 = RMAX; // doesn't have to be
      }
      else if(TORITYPE==4){
	R0=4.7*rg;
	shapeouter=3.0*R0;
	DISKTYPE=1.0; // 1.0 = H/R=const 3.0=H=const (Self-Similar solution)
      }
      

      //Omega0 // number of c^3/GM (from x=R0*Omega/c, x^2 in PHI_EFF, x in V_PHI)
      Omega0=sqrt(GM/(R0*R0*R0*pow(1.0-rgp/R0,2.0))); // forced
      KAPPA=1.0;
      
    }
    else{
      // whatever (make omega0 and GM consistent, or whatever)

      // H90 2D PP-inst
      R0=100.0;
      Omega0=1E-3;
      rho0=1.0/10./sqrt(2.);
      GM=R0*R0*R0*Omega0*Omega0*pow((1./sin(theta0)-rgp/R0),2.0);
      KAPPA=1.E-5;
    }


    
    // rest should be generally true for any units
    
#define PHI_EFF(R) (R0*R0*Omega0*Omega0*pow(fabs((R)/R0),(2.-2.*QCONST))/(2.*QCONST-2.))
#define PHI_G(r) (-GM/((r)-rgp))
#define V_PHIN(R) (R0*Omega0*pow(fabs((R)/R0),1.-QCONST))
#define V_PHIPN(R) (V_PHIN(R)*(1-rgp/R0)/(1-rgp/(R)))  // assuming PN factor is unchanged by QCONST : crazy phieffective anyways, so never use
#define V_PHI(R) (V_PHIN(R)) // for now
#define V_R(R) (0)
#define V_THETA(R) (0)
    // below line on theta=Pi/2 axis only and rgp=0
#define R_P0 (R0*pow(GM*(2.*QCONST-2.)/(Omega0*Omega0*R0*R0*R0),1./(3.-2.*QCONST)))

    if(TORITYPE==0){
      //  Rint=((0.775)*(L[1][1]+L[2][1]/4.0));
      //Rint=.53*R0;
      //Rint= .75*R0; // SPB99 d=1.125 for rgp=0,theta0=Pi/2,q=2
      //Rint=0.633975*R0; // d=1.5 SP00 RunA
      // SPB conversion: d = (2*(R_0/R_int)-1/2*(R_0/R_int)^2)^{-1}
      Rint = 3.0*rg; // H00 first axisym setup
      //      Rint=6.5; // thinner torus

      //    Rint = 0.2*rg; // H00 first axisym setup (test for symmetry problem) // GODMARK
      //Rint = 5.5*rg; // H91 (A2/A3)
      //    Rint = 5.36*rg; // My 26 test of a stable torus against perturbations
      //Rint = 4.75*rg; // H91 (C2)
      //Rint = 4.0*rg; // H00 third axisym setup
      //Rint=.53*R0;
      //Rint=0.6*R0; // standard so far
      //Rint=.85*R0;
      //      Rint=95.0; // H90 2D PP-inst
      //Rint=90.0; // H87 Model3
      Routt=Router;
    }
    else if(TORITYPE==1){
      Rint=5.7;
      //Rint=3.25913;
      Routt=Router;
    }
    else if(TORITYPE==2){
      Rint=0.634*R0;
      Routt=Router;
    }
    else if(TORITYPE==3){
      Rint=6.0;
      Routt=Router;
    }
    else if(TORITYPE==4){
      Rint=6.0;
      Routt=Router;
    }


    if(QCONST>1.5){
    }
    else{
      // set Keplerian disk
      deltatheta=0.05;
    }


    if(TORITYPE!=4){

      if(wgam){
	if(QCONST>1.5){
	  Cconst=PHI_EFF(Rint)+PHI_G(Rint/sin(M_PI/2.0));
	  u0=rho0/gam*(Cconst-(PHI_EFF(R0)+PHI_G(R0/sin(M_PI/2.0))));
	}
	else{
	  // assume want fake Keplerian disk
	  Cconst=0;
	  u0=KAPPA*pow(rho0,gam)/(gam-1);
	}
      }
      else{
	Cconst=PHI_EFF(R0)+PHI_G(R0/sin(M_PI/2.0));
	cs=sqrt(PHI_EFF(RMAX)+PHI_G(RMAX/sin(M_PI*0.5-HOR))-( PHI_EFF(RMAX)+PHI_G(RMAX/sin(M_PI/2.0)) ) );
	u0=cs*cs*rho0; // irrelevant
      }
			
      if(u0<=0){
	fprintf(fail_file,"Cannot have u0<0, choose other parameters so that solution exists\n");
	myexit(1);
      }
      //printf("%15.10g %15.10g %15.10g\n\n",l0sq,Cconst,Cconst+GRAVC*MASS/(Rint/sin(theta0)-rgp)-PHI_EFF(Rint));
    }
    // atmosphere params
    C2const=-1.0; // for equilibrium atmosphere
    if(TORITYPE==0){
      Beta=1.0; // factor on tori KAPPA for any atmosphere
    }
    else if(TORITYPE==1){
      // Beta=1.0E-3; // not really
      Beta=1.E-3;
    }
    else if(TORITYPE==2){
      Beta=1.E-3;
    }
    else if(TORITYPE==3){
      Beta=1.0;
    }
    else if(TORITYPE==4){
      Beta=1.0;
    }
    // used for non adiabatic atm
    IEFRACT=0.2;
    VZFRACT=.95;

    if((TORITYPE==0)||(TORITYPE==3)||(TORITYPE==4)){
      // normal
      DENSITYFRACFLOOR=1.E-4;
      IEFRACFLOOR=DENSITYFRACFLOOR;
      DENSITYFLOOR=DENSITYFRACFLOOR*rho0;
      //      DENSITYFLOOR=1.E-6*rho0;
      //    IEFRACFLOOR=DENSITYFRACFLOOR=1.E-6;
    }
    else if(TORITYPE==1){
      DENSITYFLOOR=1E-4*rho0;
      IEFRACFLOOR=DENSITYFRACFLOOR=1.E-4;
    }
    else if(TORITYPE==2){
      // SP00, SPB99
      DENSITYFLOOR=1.E-4*rho0;
      IEFRACFLOOR=DENSITYFRACFLOOR=1.E-4;
    }

    if((TORITYPE==0)||(TORITYPE==3)||(TORITYPE==4)){
      if(ADATM){
	IEFLOOR=u0*pow(DENSITYFLOOR/rho0,gam);
      }
      else{
	IEFLOOR=IEFRACT*fabs(DENSITYFLOOR*PHI_G(fabs(x1out)));
      }
    }
    else if(TORITYPE==1){
      IEFLOOR=1.0E-6;
    }
    else if(TORITYPE==2){
      IEFLOOR=u0*DENSITYFLOOR/rho0;
    }



    ////////////////////////////////////////
    ////////////////////////////////////////
    /// SET THE SOLUTION
    ///
    ////////////////////////////////////////
    ////////////////////////////////////////
    
    if(myid<=0){
      if(!wpw){
	Rpot0=R_P0; 
	fprintf(analyticout,"GravPotential==0 @ R=%15.10g\n",Rpot0);
	if(Rint<Rpot0) fprintf(analyticout,"Warning, your set inner torus radius is smaller than the 0 of the grav pot, so high dumping will occur at large radius");
      }
      else{
	fprintf(analyticout,"No check to make sure Rint>Rpot0\n");
      }
    }
    if(Rint>R0){
      fprintf(fail_file,"Invalid to have Rint>R0\n");
      myexit(1);
    }
    if(rgp>R0){
      fprintf(fail_file,"Invalid to have rgp>R0\n");
      myexit(1);
    }
    
    totalmass[0]=0.0;
    totalmass[1]=0.0;

    initialize_gravity(sca3,vx3,vy3,vz3); // optimizations
    
    // set the primitive variables, except magnetic field(which is derived from density)
    for(l=1;l<=1+3;l++){ // 1=scalars 2=v_1 3=v_2 4=v_3 , since different positions and want IC to be exactly numerically symmetric

      if(l==1){
	radiustouse=sca3[1];
	thetatouse=sca3[2];
	phitouse=sca3[3];
      }
      else if(l==2){
	radiustouse=vx3[1];
	thetatouse=vx3[2];
	phitouse=vx3[3];
      }
      else if(l==3){
	radiustouse=vy3[1];
	thetatouse=vy3[2];
	phitouse=vy3[3];
      }
      else if(l==4){
	radiustouse=vz3[1];
	thetatouse=vz3[2];
	phitouse=vz3[3];
      }
      fprintf(log_file,"%d here 1\n",l); fflush(log_file);
	
      LOOPF{
	torigrid=1; // initialize
	// setup grav pot and phi_eff for density calculation
	if(radiustouse[k][j][i]>rgp){
	  gravpottemp=PHI_G(radiustouse[k][j][i]);
	  phi_eff=PHI_EFF(radiustouse[k][j][i]*sin(thetatouse[k][j][i]));
	}
	else{ // potential inside never used assuming loops/bc setup right, so just for plotting purposes
	  gravpottemp=PHI_G(rgp*RGPFUDGE);
	  phi_eff=PHI_EFF(rgp*RGPFUDGE*sin(thetatouse[k][j][i]));
	  torigrid=0;
	}

	if(torigrid&&(ATMONLYT==0)&&(fabs(radiustouse[k][j][i]*sin(thetatouse[k][j][i]))>=Rint)&&(fabs(radiustouse[k][j][i]*sin(thetatouse[k][j][i]))<=Routt) ){ // then assume within the torus  //  need since some solutions have solutions beyond tori at inner region

	  if(QCONST>1.5){
	    if(wgam){
	      ftemp=(Cconst-gravpottemp-phi_eff)/(gam*u0/rho0);
	      if(ftemp>=0.0) torigrid=1;
	      else torigrid=0;
	    }
	    else{
	      ftemp=(Cconst-gravpottemp-phi_eff)/(cs*cs);
	      torigrid=1;
	    }
	  }
	  else{
	    ftemp=-pow((thetatouse[k][j][i]-theta0)/deltatheta,2.0)-pow((radiustouse[k][j][i]-R0)/(R0-Rint),2.0);
	  }
	  if(TORITYPE==3){
	    // new mod.
	    ftemp=-pow((radiustouse[k][j][i]-R0)/(R0/3.0),2.0)-pow((thetatouse[k][j][i]-M_PI*0.5)/(HOR),2.0);
	  }
	  if(TORITYPE==4){
	    // new mod also
	    ftemp=pow(radiustouse[k][j][i]*sin(thetatouse[k][j][i])/R0,-QCONST+DISKTYPE)*exp(-0.5/HOR/HOR*pow(radiustouse[k][j][i]*cos(thetatouse[k][j][i])/R0,2.0)*pow(radiustouse[k][j][i]*sin(thetatouse[k][j][i])/R0,-2.0*QCONST+DISKTYPE));
	    if(radiustouse[k][j][i]*sin(thetatouse[k][j][i])<shapeouter){ // already only true if R>Rint
	      ftemp*=(1.0-sqrt(Rint/(radiustouse[k][j][i]*sin(thetatouse[k][j][i]))))/(1.0-sqrt(Rint/(shapeouter))); // shaping function
	    }
	  }
	}
	else torigrid=0;
	
	
	// do scalars
	if(l==1){
	  
	  if(torigrid){// then torus
	    // density
	    if(QCONST>1.5){
	      if(wgam)	    sanal[1][k][j][i] = rho=rho0*pow(ftemp,1./(gam-1.));
	      else 	    sanal[1][k][j][i] = rho=rho0*exp(ftemp);
	    }
	    else sanal[1][k][j][i]=rho=rho0*exp(ftemp);
	    if(TORITYPE==4){
	      sanal[1][k][j][i]=rho=rho0*ftemp;
	    }

	    if(RANDOMIZE){ // uncorrelated(not even adiabatic) random perturbations that scale with resolution such that at N3=32 the perturbation is =randomization of the initial value at R0 for rho, en, and vphi
	      if(VRANDOMIZE==1) sanal[1][k][j][i]+=rho0*randomization*sqrt((SFTYPE)(N3)/32.0)*(ranc(0)-.5);
	      else if(VRANDOMIZE==2) sanal[1][k][j][i]*=(1. + 1.e-2*(ranc(0)-0.5)) ;
	      //	      else if(VRANDOMIZE==3){ sanal[1][k][j][i]+=rho0*MODEAMP*sin(MODEPURT*phitouse[k][j][i]); fprintf(stderr,"%d %d %d %15.10g %15.10g\n",k,j,i,phitouse[k][j][i],rho0*MODEAMP*sin(MODEPURT*phitouse[k][j][i])); fflush(stderr); }
	      else if(VRANDOMIZE==3) sanal[1][k][j][i]+=rho0*MODEAMP_rho*sin(MODEPURT*phitouse[k][j][i]);
	    }
	    // floor check applies to torus too! and make sure randomize doesn't go below floor
	    if(sanal[1][k][j][i]<DENSITYFLOOR) sanal[1][k][j][i]=DENSITYFLOOR;
	    
	    // mass add up
	    if( (i>=0)&&(i<N1)&&(j>=0)&&(j<N2)&&(k>=0)&&(k<N3) ){// only add that which is on real grid
	      totalmass[1]+=sanal[1][k][j][i]*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k);
	    }
	    
	    // internal energy (uses unperturbed rho)
	    if(wgam){
	      sanal[2][k][j][i] = u0*pow(rho/rho0,gam); //en
	    }
	    else sanal[2][k][j][i] = alpha*cs*cs*rho; //en
	    if(TORITYPE==4){
	      sanal[2][k][j][i]=rho*pow(HOR*Omega0*R0,2.0)/(gam*(gam-1.0))*pow(radiustouse[k][j][i]*sin(thetatouse[k][j][i])/R0,-DISKTYPE);
	    }

	    if(RANDOMIZE){
	      if(VRANDOMIZE==1) sanal[2][k][j][i]+=u0*randomization*sqrt((SFTYPE)(N3)/32.0)*(ranc(0)-.5);
	      //	      else if(VANDOMIZE==2) 
	      else if(VRANDOMIZE==3) sanal[2][k][j][i]+=u0*MODEAMP_u*sin(MODEPURT*phitouse[k][j][i]);
	    }
	    // floor check applies to torus too!
	    if(sanal[2][k][j][i]<IEFLOOR) sanal[2][k][j][i]=IEFLOOR;
	    
	  }
	  else{ // then atmosphere
	    if((TORITYPE==0)||(TORITYPE==3)||(TORITYPE==4)){
	      // density
	      
	      if(EQATM==1){
		// implement equilibrium atmosphere
		if(wgam){
		  ftemp=(C2const-gravpottemp)/(gam*u0/rho0);
		}
		else ftemp=-1; //forces other solution

		if(ftemp>=0.0){
		  sanal[1][k][j][i]=pow((C2const-gravpottemp)/(gam*u0/rho0),1./(gam-1.));
		}
		else{
		  sanal[1][k][j][i]=DENSITYFLOOR*ATMFACTOR;
		}
	      }
	      else{
		sanal[1][k][j][i]=DENSITYFLOOR*ATMFACTOR;
	      }
	      
	      // internal energy
	      RBeta=Beta;
	      if(ADATM){
		if(wgam){
		  sanal[2][k][j][i] = RBeta*u0*pow(sanal[1][k][j][i]/rho0,gam); //en
		}
		else sanal[2][k][j][i] = RBeta*u0*sanal[1][k][j][i]/rho0; //en
	      }
	      else{
		sanal[2][k][j][i] = IEFRACT*fabs(sanal[1][k][j][i]*gravpottemp); //en
	      }

	      sanal[1][k][j][i] = DENSITYFLOOR*ATMFACTOR;
	      sanal[2][k][j][i] = IEFLOOR*ATMFACTOR;

	    }
	    else if(TORITYPE==1){
	      sanal[1][k][j][i] = DENSITYFLOOR ;
	      sanal[2][k][j][i] = IEFLOOR ;
	    }
	    else if(TORITYPE==2){
	      sanal[1][k][j][i] = DENSITYFLOOR;
	      sanal[2][k][j][i] = DENSITYFLOOR*GM/radiustouse[k][j][i]/(gam-1.);
	    }

	    if( (i>=0)&&(i<N1)&&(j>=0)&&(j<N2)&&(k>=0)&&(k<N3) ){// only add that which is on real grid
	      totalmass[0]+=sanal[1][k][j][i]*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k);
	    }
	    
	  }// end atmosphere
	  
	  // grav pot
	  sanal[3][k][j][i]=gravpottemp;
	  
	  
	}// end if scalars
	if(l>=2){ // v_1,2,3, using different positions
	  if(l==2) component=1;
	  if(l==3) component=2;
	  if(l==4) component=3;
	  
	  if(torigrid){// torus
	    
	    
	    if(RANDOMIZE){ // scale for randomizations
	      vphi=V_PHI(R0)*randomization*sqrt((SFTYPE)(N3)/32.0)*(ranc(0)-.5);
	    }
	    if(COORD==3){
	      //p[i][j][UR] += 1.e-3*(ranc(0)-0.5) ;
	      if(component==1){
		vanal[1][1][k][j][i]=0.0;
		if(RANDOMIZE){
		  if(VRANDOMIZE==2) vanal[1][1][k][j][i]+=vrand*randomization*(ranc(0)-0.5) ;
		  else if(VRANDOMIZE==3) vanal[1][1][k][j][i]+=vphi*MODEAMP_v1*sin(MODEPURT*phitouse[k][j][i]);
		}
	      }
	      if(component==2){
		vanal[1][2][k][j][i]=0.0;
		if(RANDOMIZE){
		  if(VRANDOMIZE==2) vanal[1][2][k][j][i]+=vrand*randomization*(ranc(0)-0.5) ;
		  else if(VRANDOMIZE==3) vanal[1][2][k][j][i]+=vphi*MODEAMP_v2*sin(MODEPURT*phitouse[k][j][i]);
		}
	      }
	      if(component==3){
		vanal[1][3][k][j][i]=V_PHI(radiustouse[k][j][i]*sin(thetatouse[k][j][i]));
		if(RANDOMIZE){
		  if(VRANDOMIZE==1) vanal[1][3][k][j][i]+=vphi;
		  else if(VRANDOMIZE==2) vanal[1][3][k][j][i]+=vrand*randomization*(ranc(0)-0.5) ;
		  else if(VRANDOMIZE==3) vanal[1][3][k][j][i]+=vphi*MODEAMP_v3*sin(MODEPURT*phitouse[k][j][i]);
		}
	      }
	    }
	    if(COORD==1){
	      if(component==1){
		ftemp=V_PHI(R0*sin(M_PI/2.0))*(-sin(phitouse[k][j][i]));
		vanal[1][1][k][j][i]=V_PHI(radiustouse[k][j][i]*sin(thetatouse[k][j][i]))*(-sin(phitouse[k][j][i]));
		if(RANDOMIZE){
		  if(VRANDOMIZE==1) vanal[1][1][k][j][i]+=vphi*(-sin(phitouse[k][j][i]));
		  else if(VRANDOMIZE==2) vanal[1][1][k][j][i]+=vrand*randomization*(ranc(0)-0.5);
		  else if(VRANDOMIZE==3) vanal[1][1][k][j][i]+=ftemp*MODEAMP_v1*sin(MODEPURT*phitouse[k][j][i]);
		}
	      }
	      if(component==2){
		ftemp=V_PHI(R0*sin(M_PI/2.0))*(cos(phitouse[k][j][i]));
		vanal[1][2][k][j][i]=V_PHI(radiustouse[k][j][i]*sin(thetatouse[k][j][i]))*(cos(phitouse[k][j][i]));
		if(RANDOMIZE){
		  if(VRANDOMIZE==1) vanal[1][2][k][j][i]+=vphi*(cos(phitouse[k][j][i]));
		  else if(VRANDOMIZE==2) vanal[1][2][k][j][i]+=vrand*randomization*(ranc(0)-0.5);
		  else if(VRANDOMIZE==3) vanal[1][2][k][j][i]+=ftemp*MODEAMP_v2*sin(MODEPURT*phitouse[k][j][i]);
		}
	      }
	      // GODMARKGODMARK (commented- symmmetry test)
	      /*
		if(component==3){
		if(x[1][3][k]>0.0) vanal[1][3][k][j][i]=radiustouse[k][j][i];
		else if(x[1][3][k]<0.0) vanal[1][3][k][j][i]=-radiustouse[k][j][i];
		else if(x[1][3][k]==0.0) vanal[1][3][k][j][i]=0.0;
		}
	      */
	      if(component==3){
		vanal[1][3][k][j][i]=0.0;
		// no basis for VRANDOMIZE==1 or 3
		if(RANDOMIZE){
		  if(VRANDOMIZE==2) vanal[1][3][k][j][i]+=vrand*randomization*(ranc(0)-0.5);
		}
	      }
	      // GODMARK(commented)
	      /*
		if((component==2)&&(i==42)&&(k==16)){
		fprintf(stdout,"k=%d j=%d i=%d v2=%15.10g rho: %15.10g torigrid: %d\n",k,j,i,vanal[1][2][k][j][i],rho0*pow(ftemp,1./(gam-1.)),torigrid); fflush(stdout);
		}
	      */
	    }
	  }
	  else{// atmosphere
	    
	    if((TORITYPE==0)||(TORITYPE==3)||(TORITYPE==4)){
	      if(component==1) vanal[1][1][k][j][i]=0.0;
	      if(component==2) vanal[1][2][k][j][i]=0.0;
	      if(component==3) vanal[1][3][k][j][i]=0.0;
	    }
	    else if(TORITYPE==1){
	      if(component==1){
		vanal[1][1][k][j][i]=-rgp/(2.0*radiustouse[k][j][i]);
		if(RANDOMIZE&&(VRANDOMIZE==2)) vanal[1][1][k][j][i]+=vrand*randomization*(ranc(0)-0.5) ;
	      }
	      if(component==1){
		vanal[1][2][k][j][i]=0.0;
		if(RANDOMIZE&&(VRANDOMIZE==2)) vanal[1][2][k][j][i]+=vrand*randomization*(ranc(0)-0.5) ;
	      }
	      if(component==1){
		vanal[1][3][k][j][i]=0.0;
		if(RANDOMIZE&&(VRANDOMIZE==2)) vanal[1][3][k][j][i]+=vrand*randomization*(ranc(0)-0.5) ;
	      }
	    }
	    else if(TORITYPE==2){
	      if(component==1) vanal[1][1][k][j][i]=0.0;
	      if(component==2) vanal[1][2][k][j][i]=0.0;
	      if(component==3) vanal[1][3][k][j][i]=0.0;
	    }

	  }
	}// end if vectors
	
      }// end full loop
      // bound each new variable since successive variables may depend on it (e.g. multiple cpus and random fluctuations in density for field)
      fprintf(log_file,"%d here 2\n",l); fflush(log_file);
      if(l==1){
	smooth(SMOOTHSIZE_rho,DENSITYFLOOR,sanal[1]);
	bound(sanal[1],NULL,-10-1,0,0); // assumes all scalars independently computed
	smooth(SMOOTHSIZE_u,IEFLOOR,sanal[2]);
	bound(sanal[2],NULL,-10-2,0,0); // assumes all scalars independently computed
	// sanal[3] must be independent and so never needs bounding
      }
      if(l==2){
	smooth(SMOOTHSIZE_v1,0.0,vanal[1][1]);
	bound(NULL,vanal[1],0,-10-1,1);
      }
      if(l==3){
	smooth(SMOOTHSIZE_v2,0.0,vanal[1][2]);
	bound(NULL,vanal[1],0,-10-1,2);
      }
      if(l==4){
	smooth(SMOOTHSIZE_v3,0.0,vanal[1][3]);
	bound(NULL,vanal[1],0,-10-1,3);
      }
      fprintf(log_file,"%d here 2.5\n",l); fflush(log_file);    
    }// end loop over variously centered primitive variables
    fprintf(log_file,"here 3\n"); fflush(log_file);    


    // can't just do tori region to make divb=0 correclty
    if(MAGTORI&&(mag==1)){// can't do in loop above since use z2c
      LOOPF{// must initialize since only do LOOPH below
	vectorpot[1][k][j][i]=vectorpot[2][k][j][i]=vectorpot[3][k][j][i]=0;
      }
      // setup vector potential (located at corners like emf) 
      LOOPHP{ // v[2] needed on 0..N, needs vpot at 0..N (non-dir) 0..N (dir), needs rho -1..N
	// check if density is definitely within torus so alfven velocity isn't crazy at torus edge
	// need to fix for 3D (see why rho_min method doesn't work again) -- GODMARK (TODO!)
	// problem is to determine whether should go in that direction, eg. cart case
	if(TORITYPE!=1){
	  if(
	     (sanal[1][k][j][i]>DENSITYFLOOR*10.0)&& // k
	     (sanal[1][k][j][i+1]>DENSITYFLOOR*10.0)&&
	     (sanal[1][k][j][i-1]>DENSITYFLOOR*10.0)&&
	     (sanal[1][k][j+1][i]>DENSITYFLOOR*10.0)&&
	     (sanal[1][k][j+1][i+1]>DENSITYFLOOR*10.0)&&
	     (sanal[1][k][j+1][i-1]>DENSITYFLOOR*10.0)&&
	     (sanal[1][k][j-1][i]>DENSITYFLOOR*10.0)&&
	     (sanal[1][k][j-1][i+1]>DENSITYFLOOR*10.0)&&
	     (sanal[1][k][j-1][i-1]>DENSITYFLOOR*10.0)
	     ){
	    gooddensity=1;
	  }
	  else gooddensity=0;
	  if(TORITYPE==4){
	    ftemp = z2c_3(sanal[1],k,j,i)/rho0 - 0.2 ;
	    if(ftemp > 0.) gooddensity=1;
	    else gooddensity=0;
	  }
	  // vector pot should be located as emf is for symmetry purposes
	  if(GLOBALFIELD||gooddensity){
	    if(COORD==3){
	      if((TORITYPE==0)||(TORITYPE==3)||(TORITYPE==4)){
		// normal rotation
		if(TORITYPE!=4){
		  vectorpot[3][k][j][i]=torib0*z2c_3(sanal[1],k,j,i)/rho0;
		}
		else{
		  vectorpot[3][k][j][i]=torib0*(z2c_3(sanal[1],k,j,i)/rho0-.2);
		}
		vectorpot[1][k][j][i]=vectorpot[2][k][j][i]=0.0; // z2c_1 and z2c_2 should be used for these
		//vectorpot[3][k][j][i]=torib0*sanal[1][k][j][i]/rho0;
		// opposite rotation
		//vectorpot[3][k][j][i]=-torib0*z2c_3(sanal[1],k,j,i)/rho0;
		// double circle opposite rotation
		//	vectorpot[3][k][j][i]*=sin(2.0*x[2][2][j]);
		// double circle same rotation
		//vectorpot[3][k][j][i]*=fabs(sin(2.0*x[2][2][j]));
		// radial loops
		//vectorpot[3][k][j][i]*=sin(2.0*M_PI*4.0/(Routt-Rint)*(x[2][1][i]-Rint));
	      }

	      else if(TORITYPE==2){
		// SP01 version
		vectorpot[3][k][j][i]=torib0*pow(z2c_3(sanal[1],k,j,i)/rho0,2.0);
	      }
	    }
	    if(COORD==1){
	      // normal rotation of field
	      // note that since vpot_1 = (0.5,0,0) vy3[3][k][j][i] works here
	      // note that since vpot_2 = (0,0.5,0) vx3[3][k][j][i] works here
	      vectorpot[1][k][j][i]=torib0*z2c_1(sanal[1],k,j,i)/rho0*(-sin(vy3[3][k][j][i]));
	      vectorpot[2][k][j][i]=torib0*z2c_2(sanal[1],k,j,i)/rho0*(cos(vx3[3][k][j][i]));
	      
	      vectorpot[3][k][j][i]=0.0;
	    }
	  }
	  else{
	    vectorpot[3][k][j][i]=0.0;
	    // add vertical field (- is down, + is up)
	    //      vectorpot[3][k][j][i]+= 0.5*vertb0*x[2][1][i]*G4(2,j);
	    
	    vectorpot[1][k][j][i]=vectorpot[2][k][j][i]=0.0; // z2c_1 and z2c_2 should be used for these
	  }
	  // vectorpot[1] along j=0,j=N2 must be 0 for below to be consistent in avoidance of singularity
	  //if((j==0)||(j==N2)) vectorpot[1][k][j][i]=0; // for curlvnatback22
	}
	else{// if TORITYPE==1
	  ftemp = z2c_3(sanal[1],k,j,i)/rho0 - 0.2 ;
	  if(ftemp > 0.) vectorpot[3][k][j][i] = ftemp*ftemp*x[2][1][i]*x[2][1][i]*G4(2,j) ;
	  else vectorpot[3][k][j][i]=0.0;
	  vectorpot[1][k][j][i]=vectorpot[2][k][j][i]=0.0;
	}
      }
      
      // initialze field
      LOOPF{
	vanal[2][1][k][j][i]=vanal[2][2][k][j][i]=vanal[2][3][k][j][i]=0;
      }

      // only use LOOPH so initial divb calc will be right, could use LOOP and calculation will be normal otherwise
      //LOOPFMHP{ // for curvbacknat version assuming outer region corrected (not needed if bc right)
      // need b-field at 0..N for each dir
      //      LOOPHP{ // for curvbacknat version assuming outer region corrected
      // assumes bounded before use and not used as fixed bc(shouldn't be used that way anyways)
      // must use normal diff method even if error near r=0, at least divB=0, course then one worries about volume terms used in code!
      
      //vanal[2][1][k][j][i]=curlv1(vectorpot,k,j,i);
      //if(G4(1,j)==0) vanal[2][2][k][j][i]=curlvbacknat22(vectorpot,k,j,i);
      //else vanal[2][2][k][j][i]=curlvbacknat2(vectorpot,k,j,i);
      //vanal[2][2][k][j][i]=curlv2(vectorpot,k,j,i);
      //if(G4(1,j)==0) vanal[2][2][k][j][i]=curlv22(vectorpot,k,j,i);
      //else vanal[2][2][k][j][i]=curlv2(vectorpot,k,j,i);
      //	vanal[2][3][k][j][i]=curlv3(vectorpot,k,j,i);
      // to make sure symmetric before bound (doesn't care about divb)
      LOOPN3 LOOPN2 LOOPHP1{
	// general method
	vanal[2][1][k][j][i]=curlvbacknat1(vectorpot,k,j,i);
      }
      LOOPN3 LOOPHP2 LOOPN1{
	vanal[2][2][k][j][i]=curlvbacknat2(vectorpot,k,j,i);
      }
      LOOPHP3 LOOPN2 LOOPN1{
	vanal[2][3][k][j][i]=curlvbacknat3(vectorpot,k,j,i);
      }
      // correction for divb=0 (not needed if bc is right
      /*      
	      i=OUTFULL1-1;
	      LOOPFMHP3 LOOPFMHP2{
	      vanal[2][1][k][j][i]=b1(1,vanal[2],k,j,i);
	      }
	      j=OUTFULL2-1;
	      LOOPFMHP3 LOOPFMHP1{
	      vanal[2][2][k][j][i]=b2(1,vanal[2],k,j,i);
	      }
	      if(COMPDIM==3){ // only fix if on boundary
	      k=OUTFULL3-1;
	      LOOPFMHP2 LOOPFMHP1{
	      vanal[2][3][k][j][i]=b3(1,vanal[2],k,j,i);
	      }
	      }
      */
    }
    else{
      LOOPF{
	vanal[2][1][k][j][i]=0.0;
	vanal[2][2][k][j][i]=0.0;
	vanal[2][3][k][j][i]=0.0;
      }
    }
    // no need to bound since all CPUs the same at this point
    fprintf(log_file,"here 4\n"); fflush(log_file);        
    // now determine normalization constant from beta=p_gas/B^2
    if(MAGTORI&&(mag==1)){
      if(TORITYPE!=1){
	pgtot=0.0;
	pbtot=0.0;
	LOOP{
	  if(wgam)	  pgtot+=(gam-1.0)*sanal[2][k][j][i]/(OVOL1(k,j,i)*ODX(2,3,k));
	  else pgtot+=cs*cs*sanal[1][k][j][i]/(OVOL1(k,j,i)*ODX(2,3,k));
	  pbtot+=0.5*(e2z_1(vanal[2][1],k,j,i)*e2z_1(vanal[2][1],k,j,i)+e2z_2(vanal[2][2],k,j,i)*e2z_2(vanal[2][2],k,j,i)+e2z_3(vanal[2][3],k,j,i)*e2z_3(vanal[2][3],k,j,i))/(OVOL1(k,j,i)*ODX(2,3,k));
	}
	// need to give all cpus the normalization constant
	if(numprocs>1){
#if(USEMPI)
	  MPI_Allreduce(&(pgtot), &(pgtot_full), 1, MPI_SFTYPE, MPI_SUM,MPI_COMM_WORLD);
	  MPI_Allreduce(&(pbtot), &(pbtot_full), 1, MPI_SFTYPE, MPI_SUM,MPI_COMM_WORLD);
#endif
	}
	else{
	  pgtot_full=pgtot;
	  pbtot_full=pbtot;
	}
      }
      else{ // TORITYPE==1
	pgtot=pbtot=0.0;
	LOOP{
	  if(wgam)	  ftemp1=(gam-1.0)*sanal[2][k][j][i];
	  else 	  ftemp1=cs*cs*sanal[1][k][j][i];
	  ftemp2=0.5*(e2z_1(vanal[2][1],k,j,i)*e2z_1(vanal[2][1],k,j,i)+e2z_2(vanal[2][2],k,j,i)*e2z_2(vanal[2][2],k,j,i)+e2z_3(vanal[2][3],k,j,i)*e2z_3(vanal[2][3],k,j,i));
	  if(ftemp1>pgtot) pgtot=ftemp1;
	  if(ftemp2>pbtot) pbtot=ftemp2;
	}
	// need to give all cpus the normalization constant
	if(numprocs>1){
#if(USEMPI)
	  MPI_Allreduce(&(pgtot), &(pgtot_full), 1, MPI_SFTYPE, MPI_MAX,MPI_COMM_WORLD);
	  MPI_Allreduce(&(pbtot), &(pbtot_full), 1, MPI_SFTYPE, MPI_MAX,MPI_COMM_WORLD);
#endif
	}
	else{
	  pgtot_full=pgtot;
	  pbtot_full=pbtot;
	}
      }
      b0=sqrt(pgtot_full/(pbtot_full*magbeta));
    }
    LOOPF{
      vanal[2][1][k][j][i]*=b0;
      vanal[2][2][k][j][i]*=b0;
      vanal[2][3][k][j][i]*=b0;
    }
    //bound(NULL,vanal[2],0,-10-2,123); // assumes all field components computed from some primitive function independent of any field component (no need since last one and gets bounded in a moment)
    
    fprintf(log_file,"here 5\n"); fflush(log_file);        
    
    fprintf(log_file,"DENSITYFLOOR: %21.15g IEFLOOR: %21.15g\n",DENSITYFLOOR,IEFLOOR); fflush(log_file);
    fprintf(log_file,"DENSITYFRACFLOOR: %21.15g IEFRACFLOOR: %21.15g\n",DENSITYFRACFLOOR,IEFRACFLOOR); fflush(log_file);
    
    
    
    
    
    // output some interesting data
    // and setup initial problem stuff
    if(firstsolve==1){
      firstsolve=0;
      
      if(numprocs>1){
#if(USEMPI)
	MPI_Reduce(&(totalmass[0]), &(totalmass_full[0]), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
	MPI_Reduce(&(totalmass[1]), &(totalmass_full[1]), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
#endif
      }
      else{
	totalmass_full[0]=totalmass[0];
	totalmass_full[1]=totalmass[1];
      }
      if(myid<=0){
	// now output some interesting analytic data
	if(MAGTORI&&(mag==1)){
	  fprintf(analyticout,"beta: %15.10g b0: %15.10g\n",magbeta,b0);
	  fprintf(analyticout,"pgtot: %15.10g pbtot: %15.10g\n",pgtot_full,pbtot_full);
	}
	fprintf(analyticout,"rho0: %15.10g R0: %15.10g u0: %15.10g Omega0: %15.10g\n",rho0,R0,u0,Omega0);
	fprintf(analyticout,"Rpot0: %15.10g Rint: %15.10g\n",Rpot0,Rint);
	fprintf(analyticout,"q: %15.10g Cconst: %15.10g\n",QCONST,Cconst);
	if(!wgam) fprintf(analyticout,"cs: %15.10g\n",cs);
	if(RMAX!=0.0) fprintf(analyticout,"RMAX: %15.10g DELTAR: %15.10g HOR: %15.10g\n",RMAX,DELTAR,HOR);
	fprintf(analyticout,"DENSITYFLOOR: %15.10g IEFLOOR: %15.10g\n",DENSITYFLOOR,IEFLOOR);
	fprintf(analyticout,"DENSITYFRACFLOOR: %15.10g IEFRACFLOOR: %15.10g\n",DENSITYFRACFLOOR,IEFRACFLOOR);
	if(1.0*.001<DENSITYFLOOR){
	  fprintf(analyticout,"Warning, density floor is higher than .001 times of density max!\n");
	}
	
	fprintf(analyticout,"total tori mass: %15.10g total atm mass: %15.10g\n",totalmass_full[1],totalmass_full[0]);
	if(calltype!=2){
	  fclose(analyticout);
	}

	if(calltype==2){
	  return(analyticout);
	}
	else{
	  return(NULL);
	}

      }// end if cpu write
    }//end if firstsolve

  }
    
  //////////////////// VISUALIZATION settings

  // use image() to tell you how to set this(TOTALMM==1, and make sure DYNAMICMM==0)
  // for dvr=.1 Rint=.85

  /*
  // for normal test 
  for(i=0;i<ITYPES;i++){ // both views
  for(j=0;j<CTYPES;j++){ // both computations
  mms[i][j][1][0]=1.e-4;
  mms[i][j][1][1]=1.1;
  mms[i][j][2][0]=3e-08;
  mms[i][j][2][1]=14.0;
  mms[i][j][3][0]=-60.0;
  mms[i][j][3][1]=0.0;

  mmv[i][j][1][0][0]=0;
  mmv[i][j][1][0][1]=2.6;
  mmv[i][j][1][1][0]=-2.5;
  mmv[i][j][1][1][1]=.54;
  mmv[i][j][1][2][0]=-1.5;
  mmv[i][j][1][2][1]=1.5;
  mmv[i][j][1][3][0]=-7.e-7;
  mmv[i][j][1][3][1]=2.3;
    
  mmv[i][j][2][0][0]=0;
  mmv[i][j][2][0][1]=0;
  mmv[i][j][2][1][0]=0;
  mmv[i][j][2][1][1]=0;
  mmv[i][j][2][2][0]=0;
  mmv[i][j][2][2][1]=0;
  mmv[i][j][2][3][0]=0;
  mmv[i][j][2][3][1]=0;
  }
  }
  */
  // for art visc with inner radius at .01 and visc=.05
  /*
    for(i=0;i<ITYPES;i++){ // both views
    for(j=0;j<CTYPES;j++){ // both comps
    mms[i][j][1][0]=1.e-4;
    mms[i][j][1][1]=1.1;
    mms[i][j][2][0]=3e-08;
    mms[i][j][2][1]=14.0;
    mms[i][j][3][0]=-60.0;
    mms[i][j][3][1]=0.0;

    mmv[i][j][1][0][0]=0;
    mmv[i][j][1][0][1]=9.0;
    mmv[i][j][1][1][0]=-9.0;
    mmv[i][j][1][1][1]=4.0;
    mmv[i][j][1][2][0]=-.6;
    mmv[i][j][1][2][1]=1.0;
    mmv[i][j][1][3][0]=-.1;
    mmv[i][j][1][3][1]=9.0;
    
    mmv[i][j][2][0][0]=0;
    mmv[i][j][2][0][1]=0;
    mmv[i][j][2][1][0]=0;
    mmv[i][j][2][1][1]=0;
    mmv[i][j][2][2][0]=0;
    mmv[i][j][2][2][1]=0;
    mmv[i][j][2][3][0]=0;
    mmv[i][j][2][3][1]=0;
    }
    }
  */
  //////////////////////////////

  /*
  //for(i=0;i<ITYPES;i++){ // both views
  //  for(j=0;j<CTYPES;j++){ // both comps

  // scalars

  j=0; // normal output
  i=0; // view large

  mms[i][j][1][0]=          (9./10.)*1.e-06;  
  mms[i][j][1][1]=     1.02012908;

  mms[i][j][2][0]=(9./10.)*1.33333333e-11;  
  mms[i][j][2][1]=   0.7428376443;

  mms[i][j][3][0]=   -3.265618839;
  mms[i][j][3][1]=  -0.1881626381;

  j=1;  // second type of comp
  i=0; // view large
  
  mms[i][j][1][0]=          (9./10.)*1.e-06;  
  mms[i][j][1][1]=     1.02012908;

  mms[i][j][2][0]=(9./10.)*1.33333333e-11;  
  mms[i][j][2][1]=   0.7428376443;

  mms[i][j][3][0]=   -3.265618839;
  mms[i][j][3][1]=  -0.1881626381;


  j=0; // normal comp
  i=1; // view zoom
  // for alpha=.01 with both terms
  mms[i][j][1][0]=          .0001;  
  mms[i][j][1][1]=     .2;

  mms[i][j][2][0]=.0001;  
  mms[i][j][2][1]=   0.7428376443;

  mms[i][j][3][0]=   -3.265618839;
  mms[i][j][3][1]=  -0.1881626381;

  j=1; // 2nd comp
  i=1; // view zoom
  // for alpha=.01 with both terms
  mms[i][j][1][0]=          .0001;  
  mms[i][j][1][1]=     .1;

  mms[i][j][2][0]=.0001;  
  mms[i][j][2][1]=   0.7428376443;

  mms[i][j][3][0]=   -3.265618839;
  mms[i][j][3][1]=  -0.1881626381;


  // vectors
  
  i=0; // normal view
  j=0; // 1st comp
  // available from ipar.out

  // v0
  mmv[i][j][1][0][0]=              0.0;
  mmv[i][j][1][0][1]=    10.0;
  // vx1
  mmv[i][j][1][1][0]=   -2.0;
  mmv[i][j][1][1][1]=    .1;
  // vx2
  mmv[i][j][1][2][0]=   -1.5;
  mmv[i][j][1][2][1]=    1.5;
  // vx3
  mmv[i][j][1][3][0]=-.1;
  mmv[i][j][1][3][1]=    10.0;
  
  // B0
  mmv[i][j][2][0][0]=0;
  mmv[i][j][2][0][1]=0;
  // Bx1
  mmv[i][j][2][1][0]=0;
  mmv[i][j][2][1][1]=0;
  // Bx2
  mmv[i][j][2][2][0]=0;
  mmv[i][j][2][2][1]=0;
  // Bx3
  mmv[i][j][2][3][0]=0;
  mmv[i][j][2][3][1]=0;

  i=0;
  j=1;
  // rho times v or B
  // must get from sm plots

  // rho*v0
  mmv[i][j][1][0][0]=              0.0;
  mmv[i][j][1][0][1]=    1.0;
  // rho*vx1
  mmv[i][j][1][1][0]=   -0.01;
  mmv[i][j][1][1][1]=    0.01;
  // rho*vx2
  mmv[i][j][1][2][0]=   -0.01;
  mmv[i][j][1][2][1]=    .01;
  // rho*vx3
  mmv[i][j][1][3][0]=-.1;
  mmv[i][j][1][3][1]=    10.0;

  // rho*B0
  mmv[i][j][2][0][0]=0;
  mmv[i][j][2][0][1]=0;
  // rho*Bx1
  mmv[i][j][2][1][0]=0;
  mmv[i][j][2][1][1]=0;
  // rho*Bx2
  mmv[i][j][2][2][0]=0;
  mmv[i][j][2][2][1]=0;
  // rho*Bx3
  mmv[i][j][2][3][0]=0;
  mmv[i][j][2][3][1]=0;



  //////////////// ZOOM  
  i=1; // zoom view (vectors)
  j=0; // 1st comp
  // available from ipar.out

  // v0
  mmv[i][j][1][0][0]=              0.0;
  mmv[i][j][1][0][1]=    10.0;
  // vx1
  mmv[i][j][1][1][0]=   -3.0;
  mmv[i][j][1][1][1]=    .1;
  // vx2
  mmv[i][j][1][2][0]=   -1.5;
  mmv[i][j][1][2][1]=    1.5;
  // vx3
  mmv[i][j][1][3][0]=-1;
  mmv[i][j][1][3][1]=    10.0;
  
  // B0
  mmv[i][j][2][0][0]=0;
  mmv[i][j][2][0][1]=0;
  // Bx1
  mmv[i][j][2][1][0]=0;
  mmv[i][j][2][1][1]=0;
  // Bx2
  mmv[i][j][2][2][0]=0;
  mmv[i][j][2][2][1]=0;
  // Bx3
  mmv[i][j][2][3][0]=0;
  mmv[i][j][2][3][1]=0;

  i=1; // zoom view (vectors)
  j=1;
  // rho times v or B
  // must get from sm plots

  // rho*v0
  mmv[i][j][1][0][0]=              0.0;
  mmv[i][j][1][0][1]=    1.0;
  // rho*vx1
  mmv[i][j][1][1][0]=   -0.01;
  mmv[i][j][1][1][1]=    0.01;
  // rho*vx2
  mmv[i][j][1][2][0]=   -0.01;
  mmv[i][j][1][2][1]=    .01;
  // rho*vx3
  mmv[i][j][1][3][0]=-.1;
  mmv[i][j][1][3][1]=    1.0;

  // rho*B0
  mmv[i][j][2][0][0]=0;
  mmv[i][j][2][0][1]=0;
  // rho*Bx1
  mmv[i][j][2][1][0]=0;
  mmv[i][j][2][1][1]=0;
  // rho*Bx2
  mmv[i][j][2][2][0]=0;
  mmv[i][j][2][2][1]=0;
  // rho*Bx3
  mmv[i][j][2][3][0]=0;
  mmv[i][j][2][3][1]=0;

  */
  // mhdtorin1
  for(outtype=0;outtype<=1;outtype++){
    mms[outtype][0][1][0]=           9.9e-06 ;
    mms[outtype][1][1][0]=           9.9e-06 ;
    mms[outtype][0][1][1]=    1.0 ;
    mms[outtype][1][1][1]=    1.0 ;
    mms[outtype][0][2][0]=  1E-11 ;
    mms[outtype][1][2][0]=  1E-11 ;
    mms[outtype][0][2][1]=  0.004 ;
    mms[outtype][1][2][1]=  0.004 ;
    mms[outtype][0][3][0]=   -2.56;
    mms[outtype][1][3][0]=   -2.56;
    mms[outtype][0][3][1]=  -0.04818533834 ;
    mms[outtype][1][3][1]=  -0.04818533834 ;  
    mmv[outtype][0][1][1][0]=    -2.5 ;
    mmv[outtype][1][1][1][0]=    -2.5 ;
    mmv[outtype][0][1][1][1]=    1.0 ;
    mmv[outtype][1][1][1][1]=    1.0 ;
    mmv[outtype][0][1][2][0]=   -1.2 ;
    mmv[outtype][1][1][2][0]=   -1.2 ;
    mmv[outtype][0][1][2][1]=   1.2 ;
    mmv[outtype][1][1][2][1]=   1.2 ;
    mmv[outtype][0][1][3][0]=   -1.5 ;
    mmv[outtype][1][1][3][0]=   -1.5 ;
    mmv[outtype][0][1][3][1]=     1.5 ;
    mmv[outtype][1][1][3][1]=     1.5 ;
    mmv[outtype][0][2][1][0]=  -0.05 ;
    mmv[outtype][1][2][1][0]=  -0.05 ;
    mmv[outtype][0][2][1][1]=   0.05 ;
    mmv[outtype][1][2][1][1]=   0.05 ;
    mmv[outtype][0][2][2][0]= -0.02 ;
    mmv[outtype][1][2][2][0]= -0.02 ;
    mmv[outtype][0][2][2][1]= 0.02 ;
    mmv[outtype][1][2][2][1]= 0.02 ;
    mmv[outtype][0][2][3][0]=  -0.06 ;
    mmv[outtype][1][2][3][0]=  -0.06 ;
    mmv[outtype][0][2][3][1]=   0.06 ;
    mmv[outtype][1][2][3][1]=   0.06 ;
  }
  // define outer region when interpolation is used.
  // same order as scalar/vector arrays

  // for images
  for(i=0;i<ITYPES;i++){ // both views
    for(j=0;j<CTYPES;j++){ // both comps
      outerdefs[i][j][1]=mms[i][j][1][0]; // rho
      outerdefs[i][j][2]=mms[i][j][2][0]; // en
      outerdefs[i][j][3]=mms[i][j][3][0]; // pot
      
      outerdefv[i][j][1][0]=0.0; // magnitude of v
      outerdefv[i][j][1][1]=0.0; // v1
      outerdefv[i][j][1][2]=0.0; // v2
      outerdefv[i][j][1][3]=0.0; // v3
      
      outerdefv[i][j][2][0]=0.0;
      outerdefv[i][j][2][1]=0.0;
      outerdefv[i][j][2][2]=0.0;
      outerdefv[i][j][2][3]=0.0;
      
    }
  }
  // for dumps(never use 2nd comp type, and always same for any view)
  douterdefs[1]=mms[0][0][1][0]; // rho
  douterdefs[2]=mms[0][0][2][0]; // en
  douterdefs[3]=mms[0][0][3][0]; // pot

  douterdefv[1][0]=0.0; // magnitude of v
  douterdefv[1][1]=0.0; // v1
  douterdefv[1][2]=0.0; // v2
  douterdefv[1][3]=0.0; // v3

  douterdefv[2][0]=0.0;
  douterdefv[2][1]=0.0;
  douterdefv[2][2]=0.0;
  douterdefv[2][3]=0.0;

  /*
    LOOPF{
    if(sanal[1][k][j][i]<DENSITYFLOOR){
    fprintf(stderr,"problem\n");
    fflush(stderr);
    exit(1);
    }
    }
  */
  return(0);
}



// inject solution(starting with just atmosphere)

void injectsol(int calltype)
{
  char check[50];
  int i,j,k,l,m;
  SFTYPE temp,temp2;
  SFTYPE radsize2,tempf1,tempf2,xz0;
  SFTYPE lkep0,Cconst;
  SFTYPE C2const,Beta,RBeta;
  SFTYPE ftemp,postemp;
  char temps[100];
  char filename[100];
  SFTYPE ftemp1,ftemp2,ftempx1,ftempx2;
  int stemp,rtemp;
  SFTYPE Rpot0;
  SFTYPE R0,theta0;
  SFTYPE Rin,Rout,l0sq;
  SFTYPE Dconst;
  SFTYPE MASSDOT,IEDOT;

  SFTYPE phi_eff,phi_tot;
  static FILE*analyticout;
  static int firstsolve=1;
  SFTYPE totalmass[2]; // 0: tori 1: atmosphere
  SFTYPE totalmass_full[2]; // 0: tori 1: atmosphere
  SFTYPE massdot,massdottot,massdottot_full,massrealdottot,massrealdottot_full,massreal2dottot,massreal2dottot_full;
  SFTYPE iedot,iedottot,iedottot_full,ierealdottot,ierealdottot_full,iereal2dottot,iereal2dottot_full;
  int stoptag;
#if(USEMPI)
  static MPI_Request request[4];
#endif
  SFTYPE DENSITYINJFLOOR,ENINJFLOOR,DIX1,DOX1,DIX2,DOX2,MASSDOTTOTAL,IEDOTTOTAL;
  int ADATM,TORIINIT,PATTERN,POLARINJECT;

  

  DENSITYINJFLOOR= (1.0E-5);
  ENINJFLOOR= (1.0E-8);

  DENSITYFRACFLOOR=IEFRACFLOOR=1E-5;
  
  ADATM= 0; // 1: adiabatic atm 0: fraction pot atm
  
  TORIINIT= 0; // 1: use tori as initial data 0: Don't
  
  PATTERN= 2; // 1: gaussian 2: tori
  
  POLARINJECT= 0; // 0: normal inject 1: do dual polar injection
  
  // relevant parameters for gaussian size, and for both gaussian/tori loop restrictions
  //
  // if this is too far out, most of mass will outflow through outer boundary, making understanding of what fraction of injection the accretion is difficult
  if(PATTERN==1){
    //x1
    DIX1= ((0.8)*(L[1][1]+L[2][1]));
    DOX1= ((0.85)*(L[1][1]+L[2][1]));
    // below to start injection on inside of grid, or at Rin+Rout/4
    //	 DIX1 ((0.8)*(L[1][1]+L[2][1]/4.0));
    //	 DOX1 ((0.85)*(L[1][1]+L[2][1]/4.0));
  }
  else if(PATTERN==2){
    
    DIX1= (3.0); // Rint
    DOX1= ((L[1][1]+L[2][1]));
  }
  if(POLARINJECT==0){
    //x2
    DIX2= ((3.0/8.0)*(L[1][2]+L[2][2]));
    DOX2= ((5.0/8.0)*(L[1][2]+L[2][2]));
    //	 DIX2 ((7.0/16.0)*(L[1][2]+L[2][2]));
    //	 DOX2 ((9.0/16.0)*(L[1][2]+L[2][2]));
  }
  else{
    DIX2= ((0)*(L[1][2]+L[2][2]));
    DOX2= ((1.0/4.0)*(L[1][2]+L[2][2]));
  }
  
  

  if(calltype==100){
    x3in=  (0.0);
    x3out= (2.0*M_PI);

    mag=1;

    // H00 MHD run1
    //x1in=1.5*rg;
    x1in=1.2*rg; // avoid sonic point issue? --nope
    //x1in=(1.34666*rg);
    //x1in= (3.0*rg);
    //x1in=  (1.4*rg);
    x1out= (11.5*rg);
    x2in=  (0.0);
    x2out= (M_PI);

    mdotin=1;
    // viscosity run
    //tvisc=pow(x1out,1.5)/alpha_real0; //viscous time scales
    //numvtime=2.0;
    //tf = numvtime*tvisc; // N viscous time scales
    //tf = 10.0;
    //  tf = 4000.0;
    //tf = .245;
    //  tf = 0.240869213846978;
    //tf = .279;
    //tf = 50.0*(2.*M_PI);
    //tf     = 3.0+4.4*(2.*M_PI) ; 
    ////tf     = 3.0+2.3*(2.*M_PI) ; // end time
    //tf     = 100.0*(2.*M_PI) ; // end time
    //  timagescale=tvisc;
    //tavgi  = 3.0+2.0*(2.0*M_PI); // time to start averaging
    //tavgf  = tavgi+0.278*(2.0*M_PI); // time to stop averaging(only approximate)
    //DTl    = (1./1000.)*(2.*M_PI) ;       /* logging period(lower limit on all except DTtimestep, this is also the ener/loss file output period */

    // mhd tori up close to rg
    //  DTl    = (1.0/30.0);
    //  DTl    = (1.0/5.0);
    //DTl = DTCONS;

    // MHD TORI
    tf = 2423.4; // R0=4.7rg 17 orbits
    //tf = 1000.0; // R0=4.7rg 5 orbits
    
    simplebc=1;
    bcix1=4;
    bcox1=4;
    bcix2=1;
    bcox2=1;
    DTl=DTener=DTloss=DTtimestep=DTsp=tf/100.0;  DTfloor=DTd=tf/(20.0); DTpd=tf/(10.0); DTi=tf/10.0;
  }
  else if(calltype==0){

    if(1||(calltype==0)){// need computed values!

      if(!wgam){
	fprintf(fail_file,"No solution setup for gam=%g yet (i.e. gam==1)\n",gam);
	myexit(1);
      }

      if((PATTERN==2)||(TORIINIT==1)){
	analyticout=tori1sol(2); // submit for solution and file setup
      }
      else{    // want to write some interesting data on solution
	if(firstsolve==1){
	  if(myid<=0){
	    strcpy(temps,DATADIR);
	  
	    sprintf(filename,"%s0_analdata%s",temps,DAT2EXT);
	    if((analyticout=fopen(filename,"wt"))==NULL){
	      fprintf(fail_file,"Cannot open %s\n",temps);
	      myexit(1);
	    }
	  }//end if cpu write
	}
      }
      // setup tori conditions
    
      MASSDOT=1.0;
      IEDOT=1.0; // not really used currently
      MASSDOTTOTAL=MASSDOT; // total amount of mass per unit time injected
      IEDOTTOTAL=IEDOT; // total amount of internal energy per unit time injected
      GM=1;
      //rg=2.0; // make sure set in init.c for these units
      // MARK
      IEFRACT=0.2; // normal IGU thing
      //IEFRACT=0.6; // 200 is too extreme, crashes code eventually
      VZFRACT=.95;

      // in case want to use SPB visc with injection
      rho0=1;
      R0=DIX1;
      //R0=(1.2)*(L[1][1]+L[2][1]/4.0); // make same as tori1sol to compare
      Omega0=sqrt(GM/(R0*R0*R0*pow(1.0-rgp/R0,2.0))); // forced    
    
#define PHI_G(r) (-GM/((r)-rgp))
    
      Beta=1.0; // factor on tori KAPPA for any atmosphere
      RBeta=Beta;
      Dconst=1.0;
      Cconst=1.0;
      // MARK
      DENSITYFLOOR=1.E-10;
      //DENSITYFLOOR=1.E-5;
      //IEFLOOR=Dconst*pow(DENSITYFLOOR,gam);

      DENSITYFRACFLOOR=IEFRACFLOOR=1E-10;

      tagii=0;
      tagfi=N1;
      tagij=0;
      tagfj=N2;
      tagik=0;
      tagfk=N3;
    
      // determine region of loop
      for(i=INFULL1;i<OUTFULL1;i++){
	if(x[2][1][i]>DIX1){
	  tagii=i;
	  break;
	}
      }
      if(i==OUTFULL1){ // then never found it, so not in this domain
	tagii=i;
	tagfi=i;
	stoptag=1;
      }
      else{
	for(i=OUTFULL1-1;i>=INFULL1;i--){
	  if(x[2][1][i]<DOX1){
	    tagfi=i;
	    break;
	  }
	}
	if(i==INFULL1-1){ // then never found it, so not in this domain
	  tagii=i+1;
	  tagfi=i+1;
	  stoptag=1;
	}
	else{
	  tagfi++; // since loop says i<tagfi, not <=
	}
      }

      stoptag=0;
      // now x2-dir
      for(j=INFULL2;j<OUTFULL2;j++){
	if(x[2][2][j]>DIX2){
	  tagij=j;
	  break;
	}
      }
      if(j==OUTFULL2){
	tagij=j;
	tagfj=j;
	stoptag=1;
      }
      else{
	for(j=OUTFULL2-1;j>=INFULL2;j--){
	  if(x[2][2][j]<DOX2){
	    tagfj=j;
	    break;
	  }
	}
	if(j==INFULL2-1){
	  tagij=j+1;
	  tagfj=j+1;
	  stoptag=1;
	}
	else{
	  tagfj++;
	}
      }
      // correct for unresolved regions, moves out 1 or 2 zones(1 if on boundary already, 2 if not)
      if(tagii==tagfi){
	if(tagii>0){
	  tagii--;
	}
	if(tagfi<N1){
	  tagfi++;
	}
      }
      if(tagij==tagfj){
	if(tagij>0){
	  tagij--;
	}
	if(tagfj<N2){
	  tagfj++;
	}
      }
      if(POLARINJECT){ // do entire grid if polar injection since not localized
	tagii=0;
	tagfi=N1;
	tagij=0;
	tagfj=N2;
	tagik=0;
	tagfk=N3;
      }

      t2gii=tagii;
      t2gfi=tagfi;
      t2gij=tagij;
      t2gfj=tagfj;
      t2gik=tagik=0;
      t2gfk=tagfk=N3;

      // now correct normal tag
      if(tagii<0) tagii=0;
      if(tagii>N1) tagii=N1;
      if(tagfi<0) tagfi=0;
      if(tagfi>N1) tagfi=N1;
      if(tagij<0) tagij=0;
      if(tagij>N2) tagij=N2;
      if(tagfj<0) tagfj=0;
      if(tagfj>N2) tagfj=N2;

      // spread for rhoi,rhof in injection velocity part(interp needs rhoi/rhof there)(both since v needs -1,-1 normally, and since v really changed at (+1,+1) need rhoi/rhof there outside it for interp
      if(t2gii!=t2gfi){
	t3gii=t2gii-1;
	t3gfi=t2gfi+1;
      }
      if(t2gij!=t2gfj){
	t3gij=t2gij-1;
	t3gfj=t2gfj+1;
      }
      t3gik=t2gik;
      t3gfk=t2gfk;

      // spread (-1,-1) & (1,1) correction
      if(t3gii<INFULL1) t3gii=INFULL1;
      if(t3gij<INFULL2) t3gij=INFULL2;
      if(t3gfi>OUTFULL1) t3gfi=OUTFULL1;
      if(t3gfj>OUTFULL2) t3gfj=OUTFULL2;

      // spread for v since v really changed by next outer layer of rho interpolated between injection and no injection
      if(t2gii!=t2gfi){
	t4gii=t2gii;
	t4gfi=t2gfi+1;
      }
      if(t2gij!=t2gfj){
	t4gij=t2gij;
	t4gfj=t2gfj+1;
      }
      t4gik=t2gik;
      t4gfk=t2gfk;

      // spread (1,1) correction
      if(t4gfi>OUTFULL1) t4gfi=OUTFULL1;
      if(t4gfj>OUTFULL2) t4gfj=OUTFULL2;

      /*
      // need to setup so don't have to bound when I inject
      //if(numprocs>1){
      
      //      if(myid==0){
      //tagii=0;
      //tagfi=N1;
      //tagij=0;
      //tagfj=OUTFULL2;
      //}
      else if(myid==numprocs-1){
      tagii=0;
      tagfi=N1;
      tagij=INFULL2;
      tagfj=N2;
      }
      else{
      tagii=0;
      tagfi=N1;
      tagij=INFULL2;
      tagfj=OUTFULL2;
      }
      }
      else{
      tagii=0;
      tagfi=N1;
      tagij=0;
      tagfj=N2;
      }
      tagik=0;
      tagfk=N3;
      */
    
      // must be within comp grid, so injection interpolation works
      /*
	tagii=0;
	tagfi=N1;
	tagij=0;
	tagfj=N2;
    
	tagik=0;
	tagfk=N3;
      */    
      // NOW set the injection rate

      LOOPF{ // all zones
	rhoinject[k][j][i]=0;
	eninject[k][j][i]=0;
      }
    
      //    LOOPF{ // need boundary values for interp in step.c and esp. for multiple cpus, injection still only really occurs on grid
      // includes boundary zones so don't have to transfer rhoinject effects to other cpuss
      LOOPFINJ{
	if(PATTERN==1){
	  // half widths
	  ftemp1=0.5*(DOX1-DIX1);
	  ftemp2=0.5*(DOX2-DIX2);
	  // position minus position of central peak
	  ftempx1=x[2][1][i]-0.5*(DOX1+DIX1);
	  ftempx2=x[2][2][j]-0.5*(DOX2+DIX2);
	  // distribution
	  ftemp=-ftempx1*ftempx1/(ftemp1*ftemp1)-ftempx2*ftempx2/(ftemp2*ftemp2);
	  massdot=1.0*exp(ftemp );
	  rhoinject[k][j][i]=massdot*OVOL1(k,j,i)*ODX(2,3,k); // unnormalized
	  eninject[k][j][i]=IEFRACT*fabs(rhoinject[k][j][i]*s[3][k][j][i]);
	  if(POLARINJECT){
	    // half widths
	    ftemp1=0.5*(DOX1-DIX1);
	    ftemp2=0.5*(DOX2-DIX2);
	    // position minus position of central peak
	    ftempx1=x[2][1][i]-0.5*(DOX1+DIX1);
	    ftempx2=x[2][2][j]-(M_PI-0.5*(DOX2+DIX2));
	    // distribution
	    ftemp=-ftempx1*ftempx1/(ftemp1*ftemp1)-ftempx2*ftempx2/(ftemp2*ftemp2);
	    massdot=1.0*exp(ftemp );
	    rhoinject[k][j][i]+=massdot*OVOL1(k,j,i)*ODX(2,3,k); // unnormalized
	    eninject[k][j][i]+=IEFRACT*fabs(rhoinject[k][j][i]*s[3][k][j][i]);
	  }
	}
	else if(PATTERN==2){
	  rhoinject[k][j][i]=sanal[1][k][j][i]; // from tori1sol
	  //eninject[k][j][i]=sanal[2][k][j][i]; // not really well defined since unnormalized
	  eninject[k][j][i]=IEFRACT*fabs(rhoinject[k][j][i]*s[3][k][j][i]); // to be similar
	}
	if(rhoinject[k][j][i]<DENSITYINJFLOOR) rhoinject[k][j][i]=DENSITYINJFLOOR;
	if(eninject[k][j][i]<ENINJFLOOR) eninject[k][j][i]=ENINJFLOOR;
      }
      // ?'s:
      // normalize eninject?  Is en+en really tori en? No.

      // now only over normal grid
      massdottot=0;
      iedottot=0;
      // now only over real physical injection region
      LOOPINJ{ // if this were LOOP-> can't make this LOOPINT because this is really determining the physical amount of mass injected, regardless of accounting region
	massdottot+=rhoinject[k][j][i]*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k); // normalizer for this cpu
	iedottot+=eninject[k][j][i]*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k); // normalizer for this cpu
      }
    
    
      // unlike in most other places, the full value is needed by all CPUs
      if(numprocs>1){
#if(USEMPI)
	MPI_Allreduce(&massdottot, &massdottot_full, 1, MPI_SFTYPE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&iedottot, &iedottot_full, 1, MPI_SFTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif
      }
      else{
	massdottot_full=massdottot;
	iedottot_full=iedottot;
      }
    
      LOOPFINJ{// LOOPF would work here too
	rhoinject[k][j][i]=rhoinject[k][j][i]*(MASSDOTTOTAL/massdottot_full); // renormalize
	eninject[k][j][i]=eninject[k][j][i]*(MASSDOTTOTAL/massdottot_full); // normalized same since eninject\propto rhoinject(same as computing eninject from above rhoinject instead)
	//*(IEDOTTOTAL/iedottot_full); // renormalize
      
	//fprintf(fail_file,"%d %d %d %15.10g %15.10g\n",k,j,i,rhoinject[k][j][i],OVOL1(k,j,i));
      }

      // NOW FOR SOME CAUTIOUS COUNTING
      //
      // now count up real mass per unit time to be injected
      massrealdottot=0;
      ierealdottot=0;
      // now only over real physical injection region
      LOOPINJ{ // if this were LOOP-> can't make this LOOPINT because this is really determining the physical amount of mass injected, regardless of accounting region
	massrealdottot+=rhoinject[k][j][i]*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k); // normalizer for this cpu
	ierealdottot+=eninject[k][j][i]*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k); // normalizer for this cpu
      }
    
    
      // unlike in most other places, the full value is needed by all CPUs
      if(numprocs>1){
#if(USEMPI)
	MPI_Allreduce(&massrealdottot, &massrealdottot_full, 1, MPI_SFTYPE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&ierealdottot, &ierealdottot_full, 1, MPI_SFTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif
      }
      else{
	massrealdottot_full=massrealdottot;
	ierealdottot_full=ierealdottot;
      }

      // now count up real mass per unit time to be injected on full active grid
      massreal2dottot=0;
      iereal2dottot=0;
      // now only over real physical injection region
      LOOP{ // if this were LOOP-> can't make this LOOPINT because this is really determining the physical amount of mass injected, regardless of accounting region
	massreal2dottot+=rhoinject[k][j][i]*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k); // normalizer for this cpu
	iereal2dottot+=eninject[k][j][i]*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k); // normalizer for this cpu
      }
    
    
      // unlike in most other places, the full value is needed by all CPUs
      if(numprocs>1){
#if(USEMPI)
	MPI_Allreduce(&massreal2dottot, &massreal2dottot_full, 1, MPI_SFTYPE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&iereal2dottot, &iereal2dottot_full, 1, MPI_SFTYPE, MPI_SUM, MPI_COMM_WORLD);
#endif
      }
      else{
	massreal2dottot_full=massreal2dottot;
	iereal2dottot_full=iereal2dottot;
      }


      // setup initial atmosphere
      LOOPF{
	sanal[3][k][j][i]=PHI_G(x[2][1][i]);
      }
    
      if(ADATM){
	IEFLOOR=Dconst*pow(DENSITYFLOOR,gam);
      }
      else{
	ftemp=sanal[3][0][0][N1];
	if(numprocs>1){
#if(USEMPI)
	  MPI_Bcast(&ftemp,1,MPI_FTYPE,ncpux1-1,MPI_COMM_WORLD); // ncpux1-1 cpu is an outer radial cpu
#endif
	}
	IEFLOOR=IEFRACT*fabs(DENSITYFLOOR*ftemp);
      }
 
      totalmass[0]=0.0;
      totalmass[1]=0.0;
      LOOPF{
	if(TORIINIT==0){
	  sanal[1][k][j][i]=DENSITYFLOOR;
	  vanal[1][1][k][j][i]=0.0;
	  vanal[1][2][k][j][i]=0.0;      
	  vanal[1][3][k][j][i]=0.0;
	}	
	if( (i>=0)&&(i<N1)&&(j>=0)&&(j<N2)&&(k>=0)&&(k<N3) ){// only add that which is on real grid
	  totalmass[1]+=sanal[1][k][j][i]*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k);
	}
	if(TORIINIT==0){
	  if(ADATM){
	    if(wgam){
	      sanal[2][k][j][i] = RBeta*Dconst*pow(sanal[1][k][j][i],gam); //en
	    }
	    else sanal[2][k][j][i] = RBeta*alpha*cs*cs*sanal[1][k][j][i]; //en
	  }
	  else{
	    sanal[2][k][j][i] = IEFRACT*fabs(sanal[1][k][j][i]*sanal[3][k][j][i]); //en
	  }
	}
      }

      // output some interesting data
      // and setup initial problem stuff
      if(firstsolve==1){
	firstsolve=0;
      
	fprintf(log_file,"%7s %7s %7s %7s %7s %7s\n","tagix1","tagox1","tagix2","tagox2","tagix3","tagox3");
	fprintf(log_file,"%7d %7d %7d %7d %7d %7d\n",tagii,tagfi,tagij,tagfj,tagik,tagfk);
	fprintf(log_file,"%7s %7s %7s %7s %7s %7s\n","t2gix1","t2gox1","t2gix2","t2gox2","t2gix3","t2gox3");
	fprintf(log_file,"%7d %7d %7d %7d %7d %7d\n",t2gii,t2gfi,t2gij,t2gfj,t2gik,t2gfk);
	fprintf(log_file,"%7s %7s %7s %7s %7s %7s\n","t3gix1","t3gox1","t3gix2","t3gox2","t3gix3","t3gox3");
	fprintf(log_file,"%7d %7d %7d %7d %7d %7d\n",t3gii,t3gfi,t3gij,t3gfj,t3gik,t3gfk);
	fprintf(log_file,"%7s %7s %7s %7s %7s %7s\n","t4gix1","t4gox1","t4gix2","t4gox2","t4gix3","t4gox3");
	fprintf(log_file,"%7d %7d %7d %7d %7d %7d\n",t4gii,t4gfi,t4gij,t4gfj,t4gik,t4gfk);
	LOOPINJ{
	  fprintf(log_file,"%d %d %d %21.15g\n",k,j,i,rhoinject[k][j][i]);
	}
	fflush(log_file);

	if(numprocs>1){
#if(USEMPI)
	  MPI_Reduce(&(totalmass[0]), &(totalmass_full[0]), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
	  MPI_Reduce(&(totalmass[1]), &(totalmass_full[1]), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
#endif
	}
	else{
	  totalmass_full[0]=totalmass[0];
	  totalmass_full[1]=totalmass[1];
	}
	if(myid<=0){
	  // now output some interesting global analytic data
	
	  fprintf(analyticout,"DENSITYFLOOR: %15.10g IEFLOOR: %15.10g\n",DENSITYFLOOR,IEFLOOR);
	  fprintf(analyticout,"DENSITYFRACFLOOR: %15.10g IEFRACFLOOR: %15.10g\n",DENSITYFRACFLOOR,IEFRACFLOOR);
	  fprintf(analyticout,"IEFRACT: %15.10g VZFRACT: %15.10g\n",IEFRACT,VZFRACT);
	  fprintf(analyticout,"normalizer: %15.10g massdottotal: %15.10g\n",massdottot_full,MASSDOTTOTAL);
	  fprintf(analyticout,"massdottotal: %21.15g massdottotal2: %21.15g\n",massrealdottot_full,massreal2dottot_full);
	  fprintf(analyticout,"normalizer: %15.10g iedottotal: %15.10g\n",iedottot_full,IEDOTTOTAL);
	  fprintf(analyticout,"iedottotal: %21.15g iedottotal2: %21.15g\n",ierealdottot_full,iereal2dottot_full);
	  if(1.0*.001<DENSITYFLOOR){
	    fprintf(analyticout,"Warning, density floor is higher than .001 times of density max!\n");
	  }
	
	  fprintf(analyticout,"total tori mass: %15.10g total atm mass: %15.10g\n",totalmass_full[0],totalmass_full[1]); 
	  fclose(analyticout);
	}// end if cpu write
      }//end if firstsolve
    
    }// end if calltype==0
  

    //////////////////// VISUALIZATION settings
    //

    // scalars

    j=0; // normal output
    i=0; // outtype

    mms[i][j][1][0]=          DENSITYFLOOR;  
    mms[i][j][1][1]=     .05;

    mms[i][j][2][0]=    IEFLOOR;  
    mms[i][j][2][1]=   .006;

    mms[i][j][3][0]=   sanal[3][0][0][0];
    mms[i][j][3][1]=  0.0;

    j=1;  // second type of comp
    i=0; // view large


    mms[i][j][1][0]=          DENSITYFLOOR;  
    mms[i][j][1][1]=     .05;

    mms[i][j][2][0]=    IEFLOOR;  
    mms[i][j][2][1]=   .006;

    mms[i][j][3][0]=   sanal[3][0][0][0];
    mms[i][j][3][1]=  0.0;


    j=0; // normal comp
    i=1; // view zoom


    mms[i][j][1][0]=          DENSITYFLOOR;  
    mms[i][j][1][1]=     .05;

    mms[i][j][2][0]=    IEFLOOR;  
    mms[i][j][2][1]=   .006;

    mms[i][j][3][0]=   sanal[3][0][0][0];
    mms[i][j][3][1]=  0.0;

    j=1; // 2nd comp
    i=1; // view zoom


    mms[i][j][1][0]=          DENSITYFLOOR;  
    mms[i][j][1][1]=     .05;

    mms[i][j][2][0]=    IEFLOOR;  
    mms[i][j][2][1]=   .006;

    mms[i][j][3][0]=   sanal[3][0][0][0];
    mms[i][j][3][1]=  0.0;

    // vectors
  
    i=0; // normal view
    j=0; // 1st comp
    // available from ipar.out

    // v0
    mmv[i][j][1][0][0]=    0.0;
    mmv[i][j][1][0][1]=    1.0;
    // vx1
    mmv[i][j][1][1][0]=   -1.0;
    mmv[i][j][1][1][1]=    .1;
    // vx2
    mmv[i][j][1][2][0]=   -.3;
    mmv[i][j][1][2][1]=    .3;
    // vx3
    mmv[i][j][1][3][0]=    0.0;
    mmv[i][j][1][3][1]=    1.26;
  
    // B0
    mmv[i][j][2][0][0]=0;
    mmv[i][j][2][0][1]=0;
    // Bx1
    mmv[i][j][2][1][0]=0;
    mmv[i][j][2][1][1]=0;
    // Bx2
    mmv[i][j][2][2][0]=0;
    mmv[i][j][2][2][1]=0;
    // Bx3
    mmv[i][j][2][3][0]=0;
    mmv[i][j][2][3][1]=0;

    i=0;
    j=1;
    // rho times v or B
    // must get from sm plots

    // rho*v0
    mmv[i][j][1][0][0]=   0.0;
    mmv[i][j][1][0][1]=  .1;
    // rho*vx1
    mmv[i][j][1][1][0]=   -0.012;
    mmv[i][j][1][1][1]=    0.001;
    // rho*vx2
    mmv[i][j][1][2][0]=   -0.002;
    mmv[i][j][1][2][1]=    0.002;
    // rho*vx3
    mmv[i][j][1][3][0]=     0.0;
    mmv[i][j][1][3][1]=     0.02;

    // rho*B0
    mmv[i][j][2][0][0]=0;
    mmv[i][j][2][0][1]=0;
    // rho*Bx1
    mmv[i][j][2][1][0]=0;
    mmv[i][j][2][1][1]=0;
    // rho*Bx2
    mmv[i][j][2][2][0]=0;
    mmv[i][j][2][2][1]=0;
    // rho*Bx3
    mmv[i][j][2][3][0]=0;
    mmv[i][j][2][3][1]=0;



  
    i=1; // zoom view (vectors)
    j=0; // 1st comp
    // available from ipar.out


    // v0
    mmv[i][j][1][0][0]=    0.0;
    mmv[i][j][1][0][1]=    1.0;
    // vx1
    mmv[i][j][1][1][0]=   -1.0;
    mmv[i][j][1][1][1]=    .1;
    // vx2
    mmv[i][j][1][2][0]=   -.3;
    mmv[i][j][1][2][1]=    .3;
    // vx3
    mmv[i][j][1][3][0]=    0.0;
    mmv[i][j][1][3][1]=    1.26;
  
    // B0
    mmv[i][j][2][0][0]=0;
    mmv[i][j][2][0][1]=0;
    // Bx1
    mmv[i][j][2][1][0]=0;
    mmv[i][j][2][1][1]=0;
    // Bx2
    mmv[i][j][2][2][0]=0;
    mmv[i][j][2][2][1]=0;
    // Bx3
    mmv[i][j][2][3][0]=0;
    mmv[i][j][2][3][1]=0;

    i=1; // zoom view (vectors)
    j=1;
    // rho times v or B
    // must get from sm plots


    // rho*v0
    mmv[i][j][1][0][0]=   0.0;
    mmv[i][j][1][0][1]=  .1;
    // rho*vx1
    mmv[i][j][1][1][0]=   -0.012;
    mmv[i][j][1][1][1]=    0.001;
    // rho*vx2
    mmv[i][j][1][2][0]=   -0.002;
    mmv[i][j][1][2][1]=    0.002;
    // rho*vx3
    mmv[i][j][1][3][0]=     0.0;
    mmv[i][j][1][3][1]=     0.02;

    // rho*B0
    mmv[i][j][2][0][0]=0;
    mmv[i][j][2][0][1]=0;
    // rho*Bx1
    mmv[i][j][2][1][0]=0;
    mmv[i][j][2][1][1]=0;
    // rho*Bx2
    mmv[i][j][2][2][0]=0;
    mmv[i][j][2][2][1]=0;
    // rho*Bx3
    mmv[i][j][2][3][0]=0;
    mmv[i][j][2][3][1]=0;




    /*
    // used for when inner edge not too close to rgp, like with Igumenshchev comps
    // scalars

    j=0; // normal output
    i=0; // view large

    mms[i][j][1][0]=          1.E-5;  
    mms[i][j][1][1]=     .02;

    mms[i][j][2][0]=    5.E-7;  
    mms[i][j][2][1]=   .003;

    mms[i][j][3][0]=   sanal[3][0][0][0];
    mms[i][j][3][1]=  0.0;

    j=1;  // second type of comp
    i=0; // view large


    mms[i][j][1][0]=          1.E-5;  
    mms[i][j][1][1]=     .02;

    mms[i][j][2][0]=    5.E-7;  
    mms[i][j][2][1]=   .003;

    mms[i][j][3][0]=   sanal[3][0][0][0];
    mms[i][j][3][1]=  0.0;


    j=0; // normal comp
    i=1; // view zoom


    mms[i][j][1][0]=          1.E-5;  
    mms[i][j][1][1]=     .02;

    mms[i][j][2][0]=    5.E-7;  
    mms[i][j][2][1]=   .003;

    mms[i][j][3][0]=   sanal[3][0][0][0];
    mms[i][j][3][1]=  0.0;

    j=1; // 2nd comp
    i=1; // view zoom


    mms[i][j][1][0]=          1.E-5;  
    mms[i][j][1][1]=     .02;

    mms[i][j][2][0]=    5.E-7;  
    mms[i][j][2][1]=   .003;

    mms[i][j][3][0]=   sanal[3][0][0][0];
    mms[i][j][3][1]=  0.0;

    // vectors
  
    i=0; // normal view
    j=0; // 1st comp
    // available from ipar.out

    // v0
    mmv[i][j][1][0][0]=    0.0;
    mmv[i][j][1][0][1]=    0.1;
    // vx1
    mmv[i][j][1][1][0]=   -.04;
    mmv[i][j][1][1][1]=    .06;
    // vx2
    mmv[i][j][1][2][0]=   -.02;
    mmv[i][j][1][2][1]=    .02;
    // vx3
    mmv[i][j][1][3][0]=    0.0;
    if(alpha>0.6){
    mmv[i][j][1][3][1]=    0.1;
    }
    else{
    mmv[i][j][1][3][1]=    0.5;
    }
    // B0
    mmv[i][j][2][0][0]=0;
    mmv[i][j][2][0][1]=0;
    // Bx1
    mmv[i][j][2][1][0]=0;
    mmv[i][j][2][1][1]=0;
    // Bx2
    mmv[i][j][2][2][0]=0;
    mmv[i][j][2][2][1]=0;
    // Bx3
    mmv[i][j][2][3][0]=0;
    mmv[i][j][2][3][1]=0;

    i=0;
    j=1;
    // rho times v or B
    // must get from sm plots

    // rho*v0
    mmv[i][j][1][0][0]=   0.0;
    mmv[i][j][1][0][1]=  -.002;
    // rho*vx1
    mmv[i][j][1][1][0]=   -0.0005;
    mmv[i][j][1][1][1]=    0.0003;
    // rho*vx2
    mmv[i][j][1][2][0]=   -0.0002;
    mmv[i][j][1][2][1]=    .0002;
    // rho*vx3
    mmv[i][j][1][3][0]=     0.0;
    mmv[i][j][1][3][1]=    .02;

    // rho*B0
    mmv[i][j][2][0][0]=0;
    mmv[i][j][2][0][1]=0;
    // rho*Bx1
    mmv[i][j][2][1][0]=0;
    mmv[i][j][2][1][1]=0;
    // rho*Bx2
    mmv[i][j][2][2][0]=0;
    mmv[i][j][2][2][1]=0;
    // rho*Bx3
    mmv[i][j][2][3][0]=0;
    mmv[i][j][2][3][1]=0;



  
    i=1; // zoom view (vectors)
    j=0; // 1st comp
    // available from ipar.out


    // v0
    mmv[i][j][1][0][0]=    0.0;
    mmv[i][j][1][0][1]=    0.1;
    // vx1
    mmv[i][j][1][1][0]=   -.04;
    mmv[i][j][1][1][1]=    .06;
    // vx2
    mmv[i][j][1][2][0]=   -.02;
    mmv[i][j][1][2][1]=    .02;
    // vx3
    mmv[i][j][1][3][0]=    0.0;
    if(alpha>0.6){
    mmv[i][j][1][3][1]=    0.1;
    }
    else{
    mmv[i][j][1][3][1]=    0.5;
    }
  
    // B0
    mmv[i][j][2][0][0]=0;
    mmv[i][j][2][0][1]=0;
    // Bx1
    mmv[i][j][2][1][0]=0;
    mmv[i][j][2][1][1]=0;
    // Bx2
    mmv[i][j][2][2][0]=0;
    mmv[i][j][2][2][1]=0;
    // Bx3
    mmv[i][j][2][3][0]=0;
    mmv[i][j][2][3][1]=0;

    i=1; // zoom view (vectors)
    j=1;
    // rho times v or B
    // must get from sm plots


    // rho*v0
    mmv[i][j][1][0][0]=   0.0;
    mmv[i][j][1][0][1]=  -.002;
    // rho*vx1
    mmv[i][j][1][1][0]=   -0.0005;
    mmv[i][j][1][1][1]=    0.0003;
    // rho*vx2
    mmv[i][j][1][2][0]=   -0.0002;
    mmv[i][j][1][2][1]=    .0002;
    // rho*vx3
    mmv[i][j][1][3][0]=     0.0;
    mmv[i][j][1][3][1]=    .02;

    // rho*B0
    mmv[i][j][2][0][0]=0;
    mmv[i][j][2][0][1]=0;
    // rho*Bx1
    mmv[i][j][2][1][0]=0;
    mmv[i][j][2][1][1]=0;
    // rho*Bx2
    mmv[i][j][2][2][0]=0;
    mmv[i][j][2][2][1]=0;
    // rho*Bx3
    mmv[i][j][2][3][0]=0;
    mmv[i][j][2][3][1]=0;

    */


    // define outer region when interpolation is used.
    // same order as scalar/vector arrays

    // for images
    for(i=0;i<ITYPES;i++){ // both views
      for(j=0;j<CTYPES;j++){ // both comps
	outerdefs[i][j][1]=mms[i][j][1][0]; // rho
	outerdefs[i][j][2]=mms[i][j][2][0]; // en
	outerdefs[i][j][3]=mms[i][j][3][0]; // pot
      
	outerdefv[i][j][1][0]=0.0; // magnitude of v
	outerdefv[i][j][1][1]=0.0; // v1
	outerdefv[i][j][1][2]=0.0; // v2
	outerdefv[i][j][1][3]=0.0; // v3
      
	outerdefv[i][j][2][0]=0.0;
	outerdefv[i][j][2][1]=0.0;
	outerdefv[i][j][2][2]=0.0;
	outerdefv[i][j][2][3]=0.0;
      
      }
    }
    // for dumps(never use 2nd comp type, and always same for any view)
    douterdefs[1]=mms[0][0][1][0]; // rho
    douterdefs[2]=mms[0][0][2][0]; // en
    douterdefs[3]=mms[0][0][3][0]; // pot

    douterdefv[1][0]=0.0; // magnitude of v
    douterdefv[1][1]=0.0; // v1
    douterdefv[1][2]=0.0; // v2
    douterdefv[1][3]=0.0; // v3

    douterdefv[2][0]=0.0;
    douterdefv[2][1]=0.0;
    douterdefv[2][2]=0.0;
    douterdefv[2][3]=0.0;

  }
}

void visctermtest(int calltype)
{
  int i,j,k;

  LOOPF{
    // rho for real
    sanal[1][k][j][i]=1.0/(alpha_real*pow(x[2][1][i],4.0));

    // en for real (note, sum is over r_phi and phi_r)
    sanal[2][k][j][i]=-alpha_real/(2.*x[2][1][i]*x[2][1][i])*exp(-2.*alpha_real*M_PI*M_PI*t)*pow(cos(M_PI*x[2][1][i]),2.0);

    sanal[3][k][j][i]=0.0;

    // sigma_{r,phi} total
    vanal[1][1][k][j][i]=-alpha_real*M_PI/(x[1][1][i]*x[1][1][i]*x[1][1][i])*exp(-alpha_real*M_PI*M_PI*t)*cos(M_PI*x[1][1][i]);
    // sigma on b-grid
    vanal[1][2][k][j][i]=-alpha_real*M_PI/(x[2][1][i]*x[2][1][i]*x[2][1][i])*exp(-alpha_real*M_PI*M_PI*t)*cos(M_PI*x[2][1][i]);
    vanal[1][3][k][j][i]=alpha_real*x[2][1][i]*exp(-alpha_real*M_PI*M_PI*t)*sin(M_PI*x[2][1][i]);
  }
  // test floor
  //  IEFLOOR=1.e-6;
}


void test2sol(int calltype)
{
  int i,j,k,l,m;
  SFTYPE radsize2,tempf1,tempf2,xz0;
  SFTYPE lkep,Cconst;
  SFTYPE ftemp,postemp;
  int blob[2][3];
  int firstblob=1;

  LOOPF{
    // disk shock
    // setup disk shock parms
    //xz0=(L[1][1]+L[2][1]*.2);
    //tempf1=pow(x[2][1][i]*cos(x[2][2][j])-xz0,2.);
    //tempf2=pow(x[2][1][i]*sin(x[2][2][j])-xz0,2.);
    //tempf1=pow(x[2][1][i]-(L[1][1]+L[2][1]*.5),2.);
    //tempf2=pow(x[2][2][j]-(L[1][2]+L[2][2]*.5),2.);
    //radsize2=pow(L[2][1]/4.,2.);

    //if(tempf1+tempf2<radsize2){
    //s[1][k][j][i] += 1.0 ;
    // }
    //else{
    //s[1][k][j][i] = 1.0 ;
    //}
    // point shock
    //if((k==0)&&(j==(int)(N2*0.5))&&(i==(int)(N1*0.5))){
    //s[1][k][j][i] =2 ;
    // s[1][k][j][i-1] =2;
    // s[1][k][j-1][i-1] =2;
    // s[1][k][j-1][i] =2;
    //}
    //else s[1][k][j][i] = 1.0;
    //s[1][k][j][i] = 1.0;
    //if(j==N2/2) s[1][k][j][i]=1.0+.1*sin(2.0*M_PI/.25*x[2][1][i]);
    //else s[1][k][j][i]=1.0;

    //s[1][k][j][i] = gam*(.1-s[3][k][j][i])
  }


  LOOPF{

    //v[1][1][k][j][i]=-cos(x[1][2][j]);
    //v[1][2][k][j][i]=sin(x[1][2][j]);
    
    //xz0=(L[1][1]+L[2][1]*.2);
    //tempf1=pow(x[2][1][i]*cos(x[2][2][j])-xz0,2.);
    //tempf2=pow(x[2][1][i]*sin(x[2][2][j])-xz0,2.);
    //tempf1=pow(x[2][1][i]-(L[1][1]+L[2][1]*.5),2.);
    //tempf2=pow(x[2][2][j]-(L[1][2]+L[2][2]*.5),2.);
    //radsize2=pow(L[2][1]/4.,2.);

    //if(tempf1+tempf2<radsize2){
    //v[1][2][k][j][i]-=1.0;
    //}

    //v[1][1][k][j][i]=-.1*(L[1][1]+L[2][1]-x[1][1][i]);
    //v[1][1][k][j][i]=V0/(L[1][1]+L[2][1])*x[1][1][i];
    //v[1][1][k][j][i]=0.0;
    //    for(m=1;m<=3;m++){ /* 3 regardless of DIM */
    v[2][1][k][j][i]=0.0;
    v[2][2][k][j][i]=0.0;
    v[2][3][k][j][i]=0.0;
    

    //}
  }




}



void waves(int calltype)
{
  static FILE*analyticout;
  int i,j,k,l,m;
  SFTYPE rho0,rho0db,en0,va,dva,b0,k0,v0;
  SFTYPE x1,x2,x3,kx,ky,kz,bx,by,bz,dvx,dvy,dvz;
  SFTYPE thetak,phik,thetab,phib,thetav,phiv,thetadv,phidv,thetadb,phidb;
  SFTYPE nonlinear1,nonlinear2;
  static int firstsolve=1;
  SFTYPE lambda0;
  SFTYPE X,Y,XP,YP;
  int PROBLEMTYPE,GAMMIE,PUREGAMMIE, ROTATION90;
  SFTYPE vx,vy,vz;
  char filename[100];
  int mx,my,mz;
  int mbx,mby,mbz,mvx,mvy,mvz,mdvx,mdvy,mdvz;
  int kparb,genk;
  SFTYPE determinant,determinant2;
  SFTYPE cm2,k2,kdotva,alfperiod,fastperiod,slowperiod;
  int WAVETYPE;
  SFTYPE drho,den,dbx,dby,dbz;
  SFTYPE norm0,norm1;
  int ISOTHERMAL;
  int RANDOMIZE;

  PROBLEMTYPE= (0);
  // 0: alfven waves(rotated or not)
  // 1: current sheet

  GAMMIE= (0); // only relevant with problemtype==0, whether to use gammie gridding method


#define ALFVENWAVE (1)
#define MAGNETICWAVE (2)
#define SOUNDWAVE (3)

  WAVETYPE=2;
  // whether isothermal wave
  ISOTHERMAL=0;



  if(calltype==100){



    if(POSTPROC==1){
      DYNAMICMM=0; // already read in file usually
    }
    else{
      DYNAMICMM=2;
    }

    x1in=0;
    x1out=1;
    x2in=0;
    x2out=1;
    x3in=0;
    x3out=1;

    tf=50.0;
    //    tf=8.0;
    // assuming all physics is on by default:
    if(WAVETYPE==SOUNDWAVE) mag=0;
    else mag=1;

    visc_real=0;
    res_real=0;
    //visc_real=1; vreal=5; vischeat=1; alpha_real0=0.001;
    // MARK
    //res_real=1; rreal=2; resheat= 1; // resistivity
    //res_real=1; rreal=3; resheat= 1; resist_real0=0.001; // resistivity

    mdotin=0;
    cool=0;
    nonunigridx1=0;
    nonunigridx2=0;
    simplebc=1;
    bcix1=5;
    bcox1=5;
    bcix2=5;
    bcox2=5;
    bcix3=5;
    bcox3=5;

    
    if(ISOTHERMAL){
      gam=1.0; ie=0; // isothermal
      transiex1=transiex2=transiex3=vischeat=0;
    }
    else{
      gam=5.0/3.0;
    }

    if(PROBLEMTYPE==1){
      
      DTl=DTener=DTloss=tf/5000.0;
      DTtimestep=DTsp=tf/100.0;  DTfloor=DTd=tf/(2.5806452*10.0); DTpd=tf/(2.5806452*5.0); DTi=tf/(2.5806452*300.0);
    }


    if(PROBLEMTYPE==0){
      DTl=tf/5000.0;
      DTener=DTloss=DTl;
      DTtimestep=DTsp=DTl*100.0;  DTfloor=DTd=DTl*(2.5806452*10.0); DTpd=DTl*(2.5806452*5.0); DTi=DTl*(2.5806452*30.0);
      DTd=.1;
      
      // for good sampling of current wave
      //DTl=DTloss=DTener=.01;
      // for linear waves, good analysis of outputted data may require forcing the dt_{comp}->smaller.
      //cour=cour/2.0;

    }
    // for long cs problem
    //  DTi = tf/(5000.0);
    //DTl=DTener=DTloss=DTtimestep=DTsp=tf/1000.0;  DTfloor=DTd=tf/(2.5806452*10.0); DTpd=tf/(2.5806452*1.0); DTi=tf/(2.5806452*3.0);
    //DTd=0.2;
    //DTdivb=tf/100.0;
    //DTd=tf/10.0;

  }
  else if(calltype==0){ // else normal problem
    
    // want to write some interesting data on solution
    if(firstsolve==1){
      if(myid<=0){
	fprintf(logfull_file,"waves\n");
	fflush(logfull_file);
	
	sprintf(filename,"%s0_analdata%s",DATADIR,DAT2EXT);
	if((analyticout=fopen(filename,"wt"))==NULL){
	  fprintf(fail_file,"Cannot open %s\n",filename);
	  myexit(1);
	}
      }//end if cpu write
    }
    

    if(PROBLEMTYPE==0){
      if(!GAMMIE){

	//nonlinear1=nonlinear2=.01; // N1,N2	
	nonlinear1=1E-4;
	nonlinear2=1E-4; // N1,N2	
	rho0=1.0;	
	b0=1.0;
	v0=0.0;

	//DENSITYFLOOR=1.0; // density floor alfven limiter interaction test


	// when using alfven limiter rho0->rho0+b^2/c^2
	if(ALFVENLIMIT){
	  rho0db=rho0+b0*b0*invsol2_original;
	}
	else{
	  rho0db=rho0;
	}


	// use to scale ANY general wave problem
	va=b0/sqrt(rho0db);	
	dva=va*nonlinear1; // dva/va=N1
	// for non-sound wave
	if(WAVETYPE==SOUNDWAVE){
	  // for sound wave
	  cs=1; // choose
	}
	else{
	  cs=dva/nonlinear2; // dva/cs=N2
	}


	if(wgam){
	  en0=cs*cs*rho0/(gam*(gam-1.0));
	}
	else{
	  en0=cs*cs*rho0;
	}


	// gen, pick 3 numbers, mx,my,mz integers, then kx,ky,kz determined, solve for lambda	
	////////////
	// k
	//
	mx=4;
	my=0;
	mz=0;
	/*
	  if((SFTYPE)mx>(SFTYPE)totalsize[1]*THIRD){
	  fprintf(fail_file,"Unresolved mode mx=%d\n",mx);
	  myexit(1);
	  }
	  if((SFTYPE)my>(SFTYPE)totalsize[2]*THIRD){
	  fprintf(fail_file,"Unresolved mode my=%d\n",my);
	  myexit(1);
	  }
	  if((SFTYPE)mz>(SFTYPE)totalsize[3]*THIRD){
	  fprintf(fail_file,"Unresolved mode mz=%d\n",mz);
	  myexit(1);
	  }
	*/	
	kx=(SFTYPE)(mx)*2.0*M_PI/L[2][1];
	ky=(SFTYPE)(my)*2.0*M_PI/L[2][2];
	kz=(SFTYPE)(mz)*2.0*M_PI/L[2][3];
	k0=sqrt(kx*kx+ky*ky+kz*kz); // can't prefix k0 since require periodic condition generally satisifed
	lambda0=2.0*M_PI/k0;
	thetak=cart2spc(2,kx,ky,kz);	
	phik=cart2spc(3,kx,ky,kz);

	// arbitrary B
	/*
	  thetab=23.0/180.0*M_PI;
	  phib=167.0/180.0*M_PI;
	  bx=b0*(sin(thetab)*cos(phib));
	  by=b0*(sin(thetab)*sin(phib));
	  bz=b0*(cos(thetab));
	*/
	// use numbers like k
	mbx=0;
	mby=1;
	mbz=0;
	
	bx=b0*(SFTYPE)(mbx)/sqrt(mbx*mbx+mby*mby+mbz*mbz);
	by=b0*(SFTYPE)(mby)/sqrt(mbx*mbx+mby*mby+mbz*mbz);
	bz=b0*(SFTYPE)(mbz)/sqrt(mbx*mbx+mby*mby+mbz*mbz);
	
	b0=sqrt(bx*bx+by*by+bz*bz); // a check
	thetab=cart2spc(2,bx,by,bz);	
	phib=cart2spc(3,bx,by,bz);


	
	// indep of everything
	/*
	  thetav=0.0/180.0*M_PI;
	  phiv=0.0/180.0*M_PI;
	  vx=v0*(sin(thetav)*cos(phiv));
	  vy=v0*(sin(thetav)*sin(phiv));
	  vz=v0*(cos(thetav));
	*/
	// use numbers like k
	mvx=0;
	mvy=0;
	mvz=0;
	
	if((mvx!=0)&&(mvy!=0)&&(mvz!=0)){
	  vx=v0*(SFTYPE)(mvx)/sqrt(mvx*mvx+mvy*mvy+mvz*mvz);
	  vy=v0*(SFTYPE)(mvy)/sqrt(mvx*mvx+mvy*mvy+mvz*mvz);
	  vz=v0*(SFTYPE)(mvz)/sqrt(mvx*mvx+mvy*mvy+mvz*mvz);
	}
	else{
	  vx=vy=vz=0;
	}

	v0=sqrt(vx*vx+vy*vy+vz*vz); // a check
	thetav=cart2spc(2,vx,vy,vz);	
	phiv=cart2spc(3,vx,vy,vz);





	// ALFVEN WAVE ONLY

	if(WAVETYPE==ALFVENWAVE){
	  // arbitrary
	  /*
	  //thetadv=90.0/180.0*M_PI;
	    
	  //phidv=90.0/180.0*M_PI;
	    
	  thetadv=90/180.0*M_PI+thetab;
	  phidv=phib;
	    
	  dvx=dva*(sin(thetadv)*cos(phidv));
	  dvy=dva*(sin(thetadv)*sin(phidv));
	  dvz=dva*(cos(thetadv));
	    
	  */
	  // use numbers like k
	  mdvx=0;
	  mdvy=1;
	  mdvz=0;
	  
	  dvx=dva*(SFTYPE)(mdvx)/sqrt(mdvx*mdvx+mdvy*mdvy+mdvz*mdvz);
	  dvy=dva*(SFTYPE)(mdvy)/sqrt(mdvx*mdvx+mdvy*mdvy+mdvz*mdvz);
	  dvz=dva*(SFTYPE)(mdvz)/sqrt(mdvx*mdvx+mdvy*mdvy+mdvz*mdvz);
	  
	  dva=sqrt(dvx*dvx+dvy*dvy+dvz*dvz); // a check

	  dbx=sqrt(rho0db)*dvx;
	  dby=sqrt(rho0db)*dvy;
	  dbz=sqrt(rho0db)*dvz;
	  
	  // no density pertubation
	  drho=0.0;
	  den=0.0;
	  
	}
	else if(WAVETYPE==MAGNETICWAVE){
	  // general wave, use mathematica to setup (/u/jon/math/mhdwaves2.nb)

#if(ALFVENLIMIT==0)	

	  //RightFastMode
	  //omega=35.5431
	  //Period=0.176777
	  //Velocity=Velocity
	  den=-0.00015000000000000007;
	  drho=-0.00010000000000000006;
	  dvx=-0.00014142135623730956;
	  dvy=0.;
	  dvz=0.;
	  dbx=0.;
	  dby=-0.0001;
	  dbz=0.;  
#else


	  //RightFastMode
	  //omega=23.2196
	  //Period=0.270598
	  //Velocity=Velocity
	  den=0.000053033008588991035;
	  drho=0.00007071067811865475;
	  dvx=0.0000653281482438188;
	  dvy=-0.000027059805007309845;
	  dvz=-1.2996866623550045e-21;
	  dbx=0.;
	  dby=0.00007071067811865474;
	  dbz=0.;  
#endif 
	  // could be any pertubation, scaled by any 1 parameter (normalization factor, here nonlinear1.  nonlinear2 has been absorbed into choice for cs in mathematica)

	  // 1D nonlinear1/10000 2=1/10 k=64,0,0
	  /*
	  // 1D N1=128, nonlinear1/10000 2=1/10 k=64,0,0 b=0,0,1 (ie is unstable)
	  // cour=1.2 highly unstable all terms (instability saturates at 0.2=etotdiff)
	  // cour=1.0 etotdiff unstable, initial bz amp 2X high
	  // cour=0.5 strange transition at t=5.9, amp a bit too high near t=0, at trans. the etotdiff average beginds to rise
	  // cour=0.499 strange bumpy amplitude in all, like etotdiff whose average rises alot
	  // cour=0.495 bit better
	  // cour=0.49 smooth transition, bit off still on t=0 value
	  // cour=0.45 bit off on the amplitude t=0 and onward
	  // cour=0.4 good t=0 value
	  // cour=0.3 same as 0.4
	  // cour=0.1 smooth transition, still ie unstable
	  //RightFastMode
	  //omega=402.124
	  //Period=0.015625
	  den=-1.5000000000000015e-10;
	  drho=-0.00009999999999999998;
	  dvx=-0.0001000000499999875;
	  dvy=0.;
	  dvz=0.;
	  dbx=0.;
	  dby=0.;
	  dbz=-0.0001; 
	  */	

	  // alfven limiter test with kx=4 bx=1 dvy=1  
	  //RightAlfvenMode
	  //omega=25.1327
	  //Period=0.25
	  /*
	    den=0.;
	    drho=0.;
	    dvx=0.;
	    dvy=0.0001;
	    dvz=0.0;
	    dbx=0.;
	    dby=-0.0001;
	    dbz=0.;   
	  */
	  
	  // cour=1.0 : totally goofed
	  // cour=.7: wiggly amplitude, etotdiff rises gradually
	  // cour=0.5: wiggles gone, etotdiff rises
	  // cour=0.4: odd, little wiggly amplitude
	  // cour=0.3: odd, wiggles gone but amplitude at t=0 too high
	  // cour=0.2: amp back down and no wiggles
	  // cour=0.1: amp very good and no wiggles
	  // knx=32
	  //RightFastMode
	  //omega=201.062
	  //Period=0.03125
	  /*
	    den=-1.5000000000000015e-10;
	    drho=-0.00009999999999999998;
	    dvx=-0.0001000000499999875;
	    dvy=0.;
	    dvz=0.;
	    dbx=0.;
	    dby=0.;
	    dbz=-0.0001; 
	  */
	  /*
	  // cour=1.0: totally bad (sat~0.1)
	  // cour=0.7: crazy amp oscillations (en still unstable)
	  // cour=0.5: pretty damn good, small t=0 amp jump (no small oscillations!)
	  // cour=0.4: small amp oscillations
	  // cour=0.3: ""
	  // cour=0.1: very good, no small oscilations, no t=0 jump
	  //	  knx=kny=N/4
	  //RightFastMode
	  //omega=568.689
	  //Period=0.0110485
	  den=-1.5000000000000036e-10;
	  drho=-0.00009999999999999999;
	  dvx=-0.00007071071347398495;
	  dvy=-0.00007071071347398496;
	  dvz=0.;
	  dbx=0.;
	  dby=0.;
	  dbz=-0.0001;
	  */

	  /*
	  // cour=1.0: totally goofed
	  // cour=0.7: large amp oscillations
	  // cour=0.5: t=0 jump, but smooth amp
	  // cour=0.4: odd amp oscillations
	  // cour=0.2: long wavelength small oscillations
	  // cour=0.1: very good
	  //	  knx=kny=N/2
	  //RightFastMode
	  //omega=568.689
	  //Period=0.0110485
	  den=-1.5000000000000036e-10;
	  drho=-0.00009999999999999999;
	  dvx=-0.00007071071347398495;
	  dvy=-0.00007071071347398496;
	  dvz=0.;
	  dbx=0.;
	  dby=0.;
	  dbz=-0.0001;
	  */

	  /*
	  //RightFastMode
	  //omega=348.25
	  //Period=0.0180422
	  // knx=kny=knz=N/4	  // unresolved in 3D!
	  // knx=kny=knz=N/8 // decayt=.045 w/ cour=0.4
	  // cour=0.8: bit freaky then goes unstable
	  // cour=0.6: fine for while, then goes unstable
	  den=-1.5000000000000002e-10;
	  drho=-0.00010000000000000002;
	  dvx=-0.000057735055786468844;
	  dvy=-0.00005773505578646884;
	  dvz=-0.00005773505578646884;
	  dbx=-0.0000408248290463863;
	  dby=-0.0000408248290463863;
	  dbz=0.0000816496580927726;
	  */









	  // now normalize by my chosen factor (doesn't assume mathematica gives orthonormal version, but mathematica IS setup for that)
	  //norm0=sqrt(drho*drho+dvx*dvx+dvy*dvy+dvz*dvz+dbx*dbx+dby*dby+dbz*dbz); // not used anyways

	  
	  // normed by dva (allows same eigenvector but different nonlinear1/nonlinear2 factor)
	  // turn off norm for now
	  norm1=sqrt(dbx*dbx+dby*dby+dbz*dbz)/sqrt(rho0db);

	  den*=dva/norm1;
	  drho*=dva/norm1;
	  dvx*=dva/norm1;
	  dvy*=dva/norm1;
	  dvz*=dva/norm1;
	  dbx*=dva/norm1;
	  dby*=dva/norm1;
	  dbz*=dva/norm1;
	  
	}
	else if(WAVETYPE==SOUNDWAVE){
	  // general wave, use mathematica to setup (/u/jon/math/mhdwaves.nb)

	  // could be any pertubation, scaled by any 1 parameter (normalization factor, here nonlinear1.  nonlinear2 has been absorbed into choice for cs in mathematica)


	  den=0.0364302;
	  drho=0.0218581;
	  dvx=0.0218581;
	  dvy=dvz=dbx=dby=dbz=0;

	}
	thetadv=cart2spc(2,dvx,dvy,dvz);	
	phidv=cart2spc(3,dvx,dvy,dvz);

	thetadb=cart2spc(2,dbx,dby,dbz);	
	phidb=cart2spc(3,dbx,dby,dbz);

	// NOW ASSIGN VALUES

	LOOPF{
	  // sin?  re or im on exp(i(kx-wt))
	  x1=x[2][1][i];	  x2=x[2][2][j];	  x3=x[2][3][k];
	  sanal[1][k][j][i]=rho0+drho*sin(kx*x1+ky*x2+kz*x3);
	  sanal[2][k][j][i]=en0+den*sin(kx*x1+ky*x2+kz*x3);
	  sanal[3][k][j][i]=0.0;
	  
	  
	  x1=x[1][1][i];	  x2=x[2][2][j];	  x3=x[2][3][k];
	  vanal[1][1][k][j][i]=vx+dvx*sin(kx*x1+ky*x2+kz*x3);

	  x1=x[2][1][i];	  x2=x[1][2][j];	  x3=x[2][3][k];
	  vanal[1][2][k][j][i]=vy+dvy*sin(kx*x1+ky*x2+kz*x3);

	  x1=x[2][1][i];	  x2=x[2][2][j];	  x3=x[1][3][k];
	  vanal[1][3][k][j][i]=vz+dvz*sin(kx*x1+ky*x2+kz*x3);

	  x1=x[1][1][i];	  x2=x[2][2][j];	  x3=x[2][3][k];	  	  
	  vanal[2][1][k][j][i]=bx+dbx*sin(kx*x1+ky*x2+kz*x3);
	  

	  x1=x[2][1][i];	  x2=x[1][2][j];	  x3=x[2][3][k];
	  vanal[2][2][k][j][i]=by+dby*sin(kx*x1+ky*x2+kz*x3);

	  x1=x[2][1][i];	  x2=x[2][2][j];	  x3=x[1][3][k];
	  vanal[2][3][k][j][i]=bz+dbz*sin(kx*x1+ky*x2+kz*x3);
	}
	//	LOOPF{
	//  fprintf(stderr,"%d %d %d %15.10g\n",k,j,i,vanal[2][1][k][j][i]); fflush(stderr);
	//	}
	// myexit(0);
	
	// OUTPUT INTERESTING DATA
	cm2=cs*cs+va*va;
	k2=kx*kx+ky*ky+kz*kz;
	kdotva=(kx*bx+ky*by+kz*bz)/sqrt(rho0db);

	if(kdotva!=0.0){
	  alfperiod=2.0*M_PI/fabs(kdotva);
	}
	else alfperiod=-1; // tells us instead of inf

	determinant=cm2*cm2*k2*k2-4.*cs*cs*k2*kdotva*kdotva;
	if(determinant>=0.0){
	  determinant2=0.5*(cm2*k2+sqrt(determinant));
	  if(determinant2!=0.0){
	    fastperiod=2.0*M_PI/sqrt(determinant2);
	  }
	  else fastperiod=-1;
	  determinant2=0.5*(cm2*k2-sqrt(determinant));
	  if(determinant2!=0.0){
	    slowperiod=2.0*M_PI/sqrt(determinant2);
	  }
	  else slowperiod=-1;
	}
	else{
	  fastperiod=slowperiod=-1; // to tell us
	}

	fprintf(analyticout," alfperiod=%15.10g fastperiod=%15.10g slowperiod=%15.10g determinant: %15.10g determinant2: %15.10g\n",alfperiod,fastperiod,slowperiod,determinant,determinant2);
	fprintf(analyticout," cs: %15.10g va: %15.10g cm2: %15.10g k2: %15.10g kdotva: %15.10g\n",cs,va,cm2,k2,kdotva);
	fprintf(analyticout," lambda0=%15.10g rho0=%15.10g en0=%15.10g b0=%15.10g k0=%15.10g v0=%15.10g va=%15.10g dva=%15.10g nonlinear1=%15.10g nonlinear2=%15.10g\n",lambda0,rho0db,en0,b0,k0,v0,va,dva,nonlinear1,nonlinear2);
	fprintf(analyticout,"drho: %15.10g den: %15.10g\n",drho,den);
	fprintf(analyticout,"    THETA   PHI       X       Y       Z\n");
	fprintf(analyticout," k: %15.10g %15.10g   %15.10g %15.10g %15.10g\n",thetak,phik,kx,ky,kz);
	fprintf(analyticout," b: %15.10g %15.10g   %15.10g %15.10g %15.10g\n",thetab,phib,bx,by,bz);
	fprintf(analyticout,"db: %15.10g %15.10g   %15.10g %15.10g %15.10g\n",thetadb,phidb,dbx,dby,dbz);
	fprintf(analyticout," v: %15.10g %15.10g   %15.10g %15.10g %15.10g\n",thetav,phiv,vx,vy,vz);
	fprintf(analyticout,"dv: %15.10g %15.10g   %15.10g %15.10g %15.10g\n",thetadv,phidv,dvx,dvy,dvz);
	fprintf(analyticout,"TNL: x: %15.10g y: %15.10g z: %15.10g\n",fabs(1.0/(k0*dvx)),fabs(1.0/(k0*dvy)),fabs(1.0/(k0*dvz)));
	fflush(analyticout);
	fclose(analyticout);

      }
      else if(GAMMIE){
	  
	PUREGAMMIE=(0);
	rho0=1.0;
	b0=1.0;
	nonlinear1=.01; // not really this term
	
	// gammie crazy method
	LOOPF{
	  if(PUREGAMMIE){
	    X = (i + .5)*dx[1][1][i];
	    XP = i*dx[1][1][i];
	    Y = (j + .5)*dx[1][2][j];
	    YP = j*dx[1][2][j];
	  }
	  else{
	    X = x[2][1][i];
	    XP = x[1][1][i];
	    Y = x[2][2][j];
	    YP = x[1][2][j];
	  }

	
	  vanal[2][1][k][j][i] = b0 ;
	  vanal[2][2][k][j][i] = b0 ;
	  vanal[2][3][k][j][i]=0.0;
	  
	  vanal[1][1][k][j][i] = nonlinear1*sin(2.*M_PI*(XP + Y)) ;
	  vanal[1][2][k][j][i] = -nonlinear1*sin(2.*M_PI*(X + YP)) ;
	  vanal[1][3][k][j][i]=0.0;
	  
	  sanal[1][k][j][i]=rho0;    
	  sanal[2][k][j][i] = (1.)/(gam - 1.) ;
	  sanal[3][k][j][i]=0.0;
	}	
      }
    }
    else if(PROBLEMTYPE==1){
      
      rho0=1.0;
      b0=1.0;
      nonlinear1=.1; // not really this term
      lambda0=1.0;
      ROTATION90=(0);
  
      LOOPF{
	
	if(ROTATION90==0){    
	  vanal[2][1][k][j][i] = 0.0 ;
	  vanal[2][3][k][j][i]=0.0;
	  
	  vanal[1][1][k][j][i] = nonlinear1*sin(2.*M_PI*x[2][2][j]/lambda0) ;
	  vanal[1][2][k][j][i] = 0.0;
	  vanal[1][3][k][j][i]=0.0;
	  
	  sanal[1][k][j][i]=rho0;    
	  sanal[2][k][j][i] = (1.)/(gam - 1.) ;
	  sanal[3][k][j][i]=0.0;
	  
	  if(x[2][1][i]>=(L[2][1]*0.5+L[1][1])){
	    vanal[2][2][k][j][i] = -b0 ;
	  }
	  else{
	    vanal[2][2][k][j][i] = b0 ;
	  }
	}
	else{
	  vanal[2][2][k][j][i] = 0.0 ;
	  vanal[2][3][k][j][i]=0.0;
	  
	  vanal[1][1][k][j][i] = 0.0;
	  vanal[1][2][k][j][i] = nonlinear1*sin(2.*M_PI*x[2][1][i]/lambda0) ;
	  vanal[1][3][k][j][i]=0.0;
	  
	  sanal[1][k][j][i]=rho0;    
	  sanal[2][k][j][i] = (1.)/(gam - 1.) ;
	  sanal[3][k][j][i]=0.0;
	  
	  if(x[2][2][j]>=(L[2][2]*0.5+L[1][2])){
	    vanal[2][1][k][j][i] = -b0 ;
	  }
	  else{
	    vanal[2][1][k][j][i] = b0 ;
	  }
	}
      }
    }

    firstsolve=0;
    // GODMARK
    // LOOPF{
    //	  fprintf(stderr,"%d %d %d %15.10g\n",k,j,i,vanal[2][1][k][j][i]); fflush(stderr);
    //	}
    //	 myexit(0);

  }

}



void chandran(int calltype)
{
  int outtype;
  SFTYPE rho0,bx0;
  SFTYPE divb1avg,divb2avg;
  SFTYPE divb1max,divb2max;
  SFTYPE divb1avg_full,divb2avg_full;
  SFTYPE divb1max_full,divb2max_full;
  char check[50];
  int i,j,k,l,m;
  SFTYPE temp,temp2;
  SFTYPE radsize2,tempf1,tempf2,xz0;
  char temps[100];
  char filename[100];
  static FILE*analyticout;
  FILE*velinput;
  static int firstsolve=1;
  SFTYPE ftemp1,ftemp2;
  int dumi[10];
  SFTYPE dumf[10];




  if(calltype==100){
    if(POSTPROC==1){
      DYNAMICMM=0; // already read in file usually
    }
    else{
      //DYNAMICMM=0; // use data below
      DYNAMICMM=2; // testing
    }
    

    gam=1.0;
    cs=1.0;
    transiex1=transiex2=transiex3=ie=vischeat=resheat=0;
    
    mag=1;
    visc_real=0;
    res_real=1; rreal=2; resheat= 1; // resistivity
    resist_real0=.01; // from vortex tests
    //res_real= 1;  rreal=2; resheat= 1; resist_real0=0.1; // resistivity
    //res_real=1;

    if(COORD==1){
      x1in=0;
      x1out=1;
      x2in=0;
      x2out=1;
      x3in=0;
      x3out=1;

      nonunigridx1= 0; // 0=unigrid 1+=nonunigrid
      nonunigridx2= 0; // 0=unigrid 1+=nonunigrid
      nonunigridx3= 0;

      // only applies if BOUNDTYPE==1
      simplebc=1;
      bcix1=5;
      bcox1=5;
      bcix2=5;
      bcox2=5;
      bcix3=5;
      bcox3=5;

    }
    else{
      fprintf(fail_file,"change to coord=1\n");
      myexit(1);
    }
    //    cour=0.25;
    //    tf=10.0;
    tf = 50.0;
    //    tf=2000.0;
    DTl = 1.0E-2;
    // strictly for purposes of reentrant ability, should have dump/floor same DT
    DTd    = tf/500.0 ;      /* dumping period(data and analytic)*/
    DTfloor=DTd*10.0 ;      /* dumping period(data and analytic)*/
    // strictly for reentract ability, DTi should be integral number of DTds.
    DTpd = 1.0*DTd; // always reentrant data
    DTi    = DTl*5.0 ;        /* image & fieldline(Ax3) period */
    DTfld=DTi*2.0;
    //below not restricted by DTl calls, just each self
    DTener = DTl*1.0 ; // negative means do every time step (must be multiple of DTl)
    // f2's FFT shows could do 20.0 here(was 4.0)
    DTloss = DTl*1.0 ; // negative means do every time step (must be multiple of DTl
    DTmode = DTl*500.0;

    DTtimestep = DTl*100.0 ; // how often to output dominate term in timestep to failfile
    DTsp = DTl*500.0 ; // how often to output sonic point info

    DTtimescale=DTtimestep;
    DTdivb=DTl*10.0;

    // debug mode
    //DTd=DTpd=DTtimestep=tf/100.0;

  }
  else if( (calltype==0)||(calltype==2)){ // normal or injection called it

        
    // want to write some interesting data on solution
    if(firstsolve==1){
      fprintf(log_file,"chandran\n"); fflush(log_file);
      if(myid<=0){
	fprintf(logfull_file,"chandran\n"); fflush(logfull_file);
	strcpy(temps,DATADIR);
	
	sprintf(filename,"%s0_analdata%s",temps,DAT2EXT);
	if((analyticout=fopen(filename,"wt"))==NULL){
	  fprintf(fail_file,"Cannot open %s\n",filename);
	  myexit(1);
	}
      }//end if cpu write
    }
    
    
    DENSITYFLOOR=1.0E-10;
    IEFLOOR=1.0E-10;

    DENSITYFRACFLOOR=IEFRACFLOOR=1E-6;

    rho0=1;
    bx0=1;

	
    LOOPF{
      
      sanal[1][k][j][i] =rho0;
      
      // internal energy (uses unperturbed rho)
      if(wgam){
	sanal[2][k][j][i] = pow(sanal[1][k][j][i],gam); //en
      }
      else sanal[2][k][j][i] = cs*cs*rho0; //en
      
      // grav pot
      sanal[3][k][j][i]=0;
      
    }
    bound(sanal[1],NULL,-10-1,0,0); // assumes all scalars independently computed
    bound(sanal[2],NULL,-10-2,0,0); // assumes all scalars independently computed

    // read in velocities from file

#define READINITDATA 1

#if(READINITDATA)

    fprintf(log_file,"proc: %d :  opening velocity file\n",myid); fflush(log_file);
    strcpy(temps,DATADIR);    
    sprintf(filename,"%svelocities%s%s",temps,".in",myidtxt);
    if((velinput=fopen(filename,"r"))==NULL){
      fprintf(fail_file,"Cannot open %s\n",filename);
      myexit(1);
    }


    LOOP{
      //      fscanf(velinput,"%d %d %d %lf %lf %lf\n",&dumi[1],&dumi[2],&dumi[3],&dumf[1],&dumf[2],&dumf[3]);
      if(feof(velinput)){
	fprintf(fail_file,"found end of %s before should have: k=%d j=%d i=%d\n",filename,k,j,i);
	myexit(1);
      }
      fread(&dumi[1],sizeof(int),1,velinput);
      fread(&dumi[2],sizeof(int),1,velinput);
      fread(&dumi[3],sizeof(int),1,velinput);
      fread(&dumf[1],sizeof(double),1,velinput);
      fread(&dumf[2],sizeof(double),1,velinput);
      fread(&dumf[3],sizeof(double),1,velinput);

      if((dumi[3]>N3-1)||(dumi[3]<0)||(dumi[2]>N2-1)||(dumi[2]<0)||(dumi[1]>N1-1)||(dumi[1]<0)){
	fprintf(fail_file,"Problem with reading in velocity file: k=%d j=%d i=%d\n",dumi[3],dumi[2],dumi[1]);
	myexit(1);
      }

      vanal[1][1][dumi[3]][dumi[2]][dumi[1]]=dumf[1];
      vanal[1][2][dumi[3]][dumi[2]][dumi[1]]=dumf[2];
      vanal[1][3][dumi[3]][dumi[2]][dumi[1]]=dumf[3];
    }
    fprintf(log_file,"proc: %d :  done reading velocity file\n",myid); fflush(log_file);
    // velocities will be bounded

#else


#endif    
    LOOPF{
      vanal[2][1][k][j][i]=bx0;
      vanal[2][2][k][j][i]=0.0;
      vanal[2][3][k][j][i]=0.0;
    }
    
    
    // output some interesting data
    // and setup initial problem stuff
    if(firstsolve==1){
      firstsolve=0;
      if(myid<=0){
	fprintf(analyticout,"DENSITYFLOOR: %15.10g IEFLOOR: %15.10g\n",DENSITYFLOOR,IEFLOOR);
	if(calltype!=2){
	  fclose(analyticout);
	}
      }// end if cpu write
    }//end if firstsolve
  }
  
  //////////////////// VISUALIZATION settings
  for(outtype=0;outtype<=1;outtype++){
    mms[outtype][0][1][0]=           9.9e-06 ;
    mms[outtype][1][1][0]=           9.9e-06 ;
    mms[outtype][0][1][1]=    1.0 ;
    mms[outtype][1][1][1]=    1.0 ;
    mms[outtype][0][2][0]=  1E-11 ;
    mms[outtype][1][2][0]=  1E-11 ;
    mms[outtype][0][2][1]=  0.004 ;
    mms[outtype][1][2][1]=  0.004 ;
    mms[outtype][0][3][0]=   -2.56;
    mms[outtype][1][3][0]=   -2.56;
    mms[outtype][0][3][1]=  -0.04818533834 ;
    mms[outtype][1][3][1]=  -0.04818533834 ;  
    mmv[outtype][0][1][1][0]=    -2.5 ;
    mmv[outtype][1][1][1][0]=    -2.5 ;
    mmv[outtype][0][1][1][1]=    1.0 ;
    mmv[outtype][1][1][1][1]=    1.0 ;
    mmv[outtype][0][1][2][0]=   -1.2 ;
    mmv[outtype][1][1][2][0]=   -1.2 ;
    mmv[outtype][0][1][2][1]=   1.2 ;
    mmv[outtype][1][1][2][1]=   1.2 ;
    mmv[outtype][0][1][3][0]=   -1.5 ;
    mmv[outtype][1][1][3][0]=   -1.5 ;
    mmv[outtype][0][1][3][1]=     1.5 ;
    mmv[outtype][1][1][3][1]=     1.5 ;
    mmv[outtype][0][2][1][0]=  -0.05 ;
    mmv[outtype][1][2][1][0]=  -0.05 ;
    mmv[outtype][0][2][1][1]=   0.05 ;
    mmv[outtype][1][2][1][1]=   0.05 ;
    mmv[outtype][0][2][2][0]= -0.02 ;
    mmv[outtype][1][2][2][0]= -0.02 ;
    mmv[outtype][0][2][2][1]= 0.02 ;
    mmv[outtype][1][2][2][1]= 0.02 ;
    mmv[outtype][0][2][3][0]=  -0.06 ;
    mmv[outtype][1][2][3][0]=  -0.06 ;
    mmv[outtype][0][2][3][1]=   0.06 ;
    mmv[outtype][1][2][3][1]=   0.06 ;
  }
  // define outer region when interpolation is used.
  // same order as scalar/vector arrays

  // for images
  for(i=0;i<ITYPES;i++){ // both views
    for(j=0;j<CTYPES;j++){ // both comps
      outerdefs[i][j][1]=mms[i][j][1][0]; // rho
      outerdefs[i][j][2]=mms[i][j][2][0]; // en
      outerdefs[i][j][3]=mms[i][j][3][0]; // pot
      
      outerdefv[i][j][1][0]=0.0; // magnitude of v
      outerdefv[i][j][1][1]=0.0; // v1
      outerdefv[i][j][1][2]=0.0; // v2
      outerdefv[i][j][1][3]=0.0; // v3
      
      outerdefv[i][j][2][0]=0.0;
      outerdefv[i][j][2][1]=0.0;
      outerdefv[i][j][2][2]=0.0;
      outerdefv[i][j][2][3]=0.0;
      
    }
  }
  // for dumps(never use 2nd comp type, and always same for any view)
  douterdefs[1]=mms[0][0][1][0]; // rho
  douterdefs[2]=mms[0][0][2][0]; // en
  douterdefs[3]=mms[0][0][3][0]; // pot

  douterdefv[1][0]=0.0; // magnitude of v
  douterdefv[1][1]=0.0; // v1
  douterdefv[1][2]=0.0; // v2
  douterdefv[1][3]=0.0; // v3

  douterdefv[2][0]=0.0;
  douterdefv[2][1]=0.0;
  douterdefv[2][2]=0.0;
  douterdefv[2][3]=0.0;

}



void averystar(int calltype)
{
  int outtype;
  SFTYPE rho0,bx0;
  SFTYPE divb1avg,divb2avg;
  SFTYPE divb1max,divb2max;
  SFTYPE divb1avg_full,divb2avg_full;
  SFTYPE divb1max_full,divb2max_full;
  char check[50];
  int i,j,k,l,m;
  SFTYPE temp,temp2;
  SFTYPE radsize2,tempf1,tempf2,xz0;
  char temps[100];
  char filename[100];
  static FILE*analyticout;
  FILE*velinput;

  FILE*fp; //TONO
  // static FILE *vels //TONO

  static int firstsolve=1;
  SFTYPE ftemp1,ftemp2;
  int dumi[12];
  FTYPE dumf[13];
  FTYPE (*vectorpot)[N3M][N2M][N1M]; // vector potential


  double myx,myy,myz,myr;

  vectorpot=workv1;


  /////////////////////////////////////////
  //
  // Physics and code setup
  //
  /////////////////////////////////////////



  if(calltype==100){
    if(POSTPROC==1){
      DYNAMICMM=0; // already read in file usually
    }
    else{
      //DYNAMICMM=0; // use data below
      DYNAMICMM=2; // testing
    }
    

    gam=1.4;
    cs=1.0;

    invsol2=invsol2_original=1.0/pow(5.0,2);

    // BAROTROPIC:
    transiex1=transiex2=transiex3=0;
    // ADIABATIC:
    //transiex1=transiex2=transiex3=1;

    ie=1;
    vischeat=resheat=1;
    
    mag=1; // whether magnetic field is evolved
    visc_real=0;
    res_real=1; rreal=4; resheat= 0; // resistivity
    resist_real0=.01; // from vortex tests
    //res_real= 1;  rreal=2; resheat= 1; resist_real0=0.1; // resistivity
    //res_real=1;

    if(COORD==1){
#if(0)
      x1in=-1.1;
      x1out=1.1;
      x2in=-1.1;
      x2out=1.1;
      x3in=-1.1;
      x3out=1.1;
#else
      x1in=-1.5;
      x1out=1.5;
      x2in=-1.5;
      x2out=1.5;
      x3in=-1.5;
      x3out=1.5;
#endif

      nonunigridx1= 0; // 0=unigrid 1+=nonunigrid
      nonunigridx2= 0; // 0=unigrid 1+=nonunigrid
      nonunigridx3= 0;

      // only applies if BOUNDTYPE==1
      simplebc=1;
      bcix1=4;
      bcox1=4;
      bcix2=4;
      bcox2=4;
      bcix3=4;
      bcox3=4;

    }
    else{
      fprintf(fail_file,"change to coord=1\n");
      myexit(1);
    }
    //    cour=0.25;
    //    tf=10.0;
    tf = 100.0;
    //    tf=2000.0;
    DTl = tf/1000.0;
    // strictly for purposes of reentrant ability, should have dump/floor same DT
    DTd    = tf/200.0 ;      /* dumping period(data and analytic)*/
    DTfloor=DTd*10.0 ;      /* dumping period(data and analytic)*/
    // strictly for reentract ability, DTi should be integral number of DTds.
    DTpd = 1.0*DTd; // always reentrant data
    DTi    = DTl*5.0 ;        /* image & fieldline(Ax3) period */
    DTfld=DTi*2.0;
    //below not restricted by DTl calls, just each self
    DTener = DTl*1.0 ; // negative means do every time step (must be multiple of DTl)
    // f2's FFT shows could do 20.0 here(was 4.0)
    DTloss = DTl*1.0 ; // negative means do every time step (must be multiple of DTl
    DTmode = DTl*500.0;

    DTtimestep = DTl*100.0 ; // how often to output dominate term in timestep to failfile
    DTsp = DTl*500.0 ; // how often to output sonic point info

    DTtimescale=DTtimestep;
    DTdivb=DTl*10.0;

    // debug mode
    //DTd=DTpd=DTtimestep=tf/100.0;

  }
  else if( (calltype==0)||(calltype==2)){ // normal or injection called it

        
    // want to write some interesting data on solution
    if(firstsolve==1){
      fprintf(log_file,"averystar\n"); fflush(log_file);
      if(myid<=0){
	fprintf(logfull_file,"averystar\n"); fflush(logfull_file);
	strcpy(temps,DATADIR);
	
	sprintf(filename,"%s0_analdata%s",temps,DAT2EXT);
	if((analyticout=fopen(filename,"wt"))==NULL){
	  fprintf(fail_file,"Cannot open %s\n",filename);
	  myexit(1);
	}
      }//end if cpu write
    }
    
    
    DENSITYFLOOR=1.0E-10;
    IEFLOOR=1.0E-12;
    DENSITYFRACFLOOR=IEFRACFLOOR=1E-10;












    // BEGIN WORRIES


    ////////////////////////////
    //
    // solution
    //
    ////////////////////////////



    // need to define in the code
    // set means you have to read it in

    // set: sanal[1][k][j][i] = mass density
    // set: sanal[2][k][j][i] = internal energy density
    // sanal[3][k][j][i] = 0


    // set: vanal[1][1][k][j][i] = vx
    // set: vanal[1][2][k][j][i] = vy
    // set: vanal[1][3][k][j][i] = vz


    // set: vectorpot[1,2,3][k][j][i]; // A_x A_y A_z

    // vanal[2][1][k][j][i] = Bx
    // vanal[2][2][k][j][i] = By
    // vanal[2][3][k][j][i] = Bz


    // set: gravacc[1,2,3][k][j][i]= gravitational acc

    ////////////////////////////
    //
    // first read in data:
    //
    /////////////////////////////

    fprintf(log_file,"proc: %d :  opening inits file\n",myid); fflush(log_file);

#define READINFILE 1

#if(READINFILE)



    ////////////////////////////
    //
    // READ IN NORMAL FILE


    // char string;
    //FILE *manuel;
    //    FILE *"velocities.in"
    //char manuelfile3[100];
    //    char a = 'velocities.in';
    
    // open inits file
    // fprintf(stdout,"temps = %s",temps);
    // fprintf(stdout,"DATADIR = %s",DATADIR);
    strcpy(temps,DATADIR);    
    sprintf(filename,"%sicnormal%s%s",temps,".in",myidtxt);
    if((velinput=fopen(filename,"r"))==NULL){
      fprintf(fail_file,"Cannot open %s\n",filename);
      myexit(1);
    }

    // initialize vector potential
    LOOPF{
      vectorpot[1][k][j][i]=vectorpot[2][k][j][i]=vectorpot[3][k][j][i]=0;
    }


    



    // only over "active" grid
    LOOP{


      if(feof(velinput)){
        fprintf(fail_file,"found end of %s before should have: k=%d j=%d i=%d\n",filename,k,j,i);
        myexit(1);
      }

#if(FLOATTYPE==0)
      fscanf(velinput,"%d %d %d %f %f %f %f %f %f %f %f %f %f %f %f\n",&dumi[1],&dumi[2],&dumi[3],&dumf[1],&dumf[2],&dumf[3],&dumf[4],&dumf[5],&dumf[6],&dumf[7],&dumf[8],&dumf[9],&dumf[10],&dumf[11],&dumf[12]);  
#elif(FLOATTYPE==1)
      fscanf(velinput,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&dumi[1],&dumi[2],&dumi[3],&dumf[1],&dumf[2],&dumf[3],&dumf[4],&dumf[5],&dumf[6],&dumf[7],&dumf[8],&dumf[9],&dumf[10],&dumf[11],&dumf[12]);
#endif

      /*	     
	        fread(&dumi[1],sizeof(int),1,velinput);
		fread(&dumi[2],sizeof(int),1,velinput);
		fread(&dumi[3],sizeof(int),1,velinput);
      
		fread(&dumf[1],sizeof(double),1,velinput);
		fread(&dumf[2],sizeof(double),1,velinput);
		fread(&dumf[3],sizeof(double),1,velinput);
		fread(&dumf[4],sizeof(double),1,velinput);
		fread(&dumf[5],sizeof(double),1,velinput);
		fread(&dumf[6],sizeof(double),1,velinput);
		fread(&dumf[7],sizeof(double),1,velinput);
		fread(&dumf[8],sizeof(double),1,velinput);
		fread(&dumf[9],sizeof(double),1,velinput);
		fread(&dumf[10],sizeof(double),1,velinput); // don't know why it wasn't defined
		fread(&dumf[11],sizeof(double),1,velinput); 
      */

      //      fprintf(stderr,"%d %d %d : %d %d %d\n",i,j,k,dumi[1],dumi[2],dumi[3]);
      //      fprintf(stdout,"%d %d %d %f %f %f %f %f %f %f %f %f %f %f %f\n",dumi[1],dumi[2],dumi[3],dumf[1],dumf[2],dumf[3],dumf[4],dumf[5],dumf[6],dumf[7],dumf[8],dumf[9],dumf[10],dumf[11],dumf[12]);   

      if((dumi[3]>N3-1)||(dumi[3]<0)||(dumi[2]>N2-1)||(dumi[2]<0)||(dumi[1]>N1-1)||(dumi[1]<0)){
	fprintf(fail_file,"Problem with reading in init file: k=%d j=%d i=%d\n",dumi[3],dumi[2],dumi[1]);
	myexit(1);
      }

      // in code this is ....[k][j][i]
      // dumi[3]=k dumi[2]=j dumi[1]=i
      vanal[1][1][dumi[3]][dumi[2]][dumi[1]]=dumf[1]; // v^x
      vanal[1][2][dumi[3]][dumi[2]][dumi[1]]=dumf[2]; // v^y
      vanal[1][3][dumi[3]][dumi[2]][dumi[1]]=dumf[3]; // v^z

      // test vis5d vector plots
#if(0)
      myx=-1.1 + (dumi[1]+0.5)*2.2/32.0;
      myy=-1.1 + (dumi[2]+0.5)*2.2/32.0;
      myz=-1.1 + (dumi[3]+0.5)*2.2/32.0;
      myr=sqrt(myx*myx+myy*myy+myz*myz);
      vanal[1][1][dumi[3]][dumi[2]][dumi[1]]=myx/myr;
      vanal[1][2][dumi[3]][dumi[2]][dumi[1]]=myy/myr;
      vanal[1][3][dumi[3]][dumi[2]][dumi[1]]=myz/myr;
#endif

      vectorpot[1][dumi[3]][dumi[2]][dumi[1]]=dumf[4]; // A_x
      vectorpot[2][dumi[3]][dumi[2]][dumi[1]]=dumf[5]; // A_y
      vectorpot[3][dumi[3]][dumi[2]][dumi[1]]=dumf[6]; // A_z
      
      sanal[1][dumi[3]][dumi[2]][dumi[1]]=dumf[7]; // rho
      sanal[2][dumi[3]][dumi[2]][dumi[1]]=dumf[8]; // ie den
      sanal[3][dumi[3]][dumi[2]][dumi[1]]=dumf[9]; // grav. pot. (could try using GRAVACC)

      if(GRAVACC&&(!POSTPROC)){
	for(l=1;l<=3;l++){
	  gravacc[l][dumi[3]][dumi[2]][dumi[1]]=dumf[9+l];
	}
      }

      

    }// end LOOP

    // bounds vector pot
    bound(NULL,vectorpot,0,-10-2,123); // assumes all field components computed from some primitive function independent of any field component (no need since last one and gets bounded in a moment)









    ////////////////////////////
    //
    // READ IN VECTORPOT FILE


    // char string;
    //FILE *manuel;
    //    FILE *"velocities.in"
    //char manuelfile3[100];
    //    char a = 'velocities.in';
    
    // open inits file
    // fprintf(stdout,"temps = %s",temps);
    // fprintf(stdout,"DATADIR = %s",DATADIR);
    strcpy(temps,DATADIR);    
    sprintf(filename,"%sicpot%s%s",temps,".in",myidtxt);
    if((velinput=fopen(filename,"r"))==NULL){
      fprintf(fail_file,"Cannot open %s\n",filename);
      myexit(1);
    }

    // initialize vector potential
    LOOPF{
      vectorpot[1][k][j][i]=vectorpot[2][k][j][i]=vectorpot[3][k][j][i]=0;
    }


    



    // loop over active grid PLUS 1
    for(k=0;k<=N3;k++) for(j=0;j<=N2;j++) for(i=0;i<=N1;i++){


      if(feof(velinput)){
        fprintf(fail_file,"found end of %s before should have: k=%d j=%d i=%d\n",filename,k,j,i);
        myexit(1);
      }

#if(FLOATTYPE==0)
      fscanf(velinput,"%d %d %d %f %f %f\n",&dumi[1],&dumi[2],&dumi[3],&dumf[4],&dumf[5],&dumf[6]);  
#elif(FLOATTYPE==1)
      fscanf(velinput,"%d %d %d %lf %lf %lf\n",&dumi[1],&dumi[2],&dumi[3],&dumf[4],&dumf[5],&dumf[6]);
#endif

      if((dumi[3]>N3)||(dumi[3]<0)||(dumi[2]>N2)||(dumi[2]<0)||(dumi[1]>N1)||(dumi[1]<0)){
	fprintf(fail_file,"Problem with reading in init file: k=%d j=%d i=%d\n",dumi[3],dumi[2],dumi[1]);
	myexit(1);
      }

      vectorpot[1][dumi[3]][dumi[2]][dumi[1]]=dumf[4]; // A_x
      vectorpot[2][dumi[3]][dumi[2]][dumi[1]]=dumf[5]; // A_y
      vectorpot[3][dumi[3]][dumi[2]][dumi[1]]=dumf[6]; // A_z

      

    }// end LOOP



#endif







#define TEST3DCODE 0


#if(TEST3DCODE)

    // initialize vector potential
    LOOPF{
      vectorpot[1][k][j][i]=vectorpot[2][k][j][i]=vectorpot[3][k][j][i]=0;
    }


    // only over "active" grid
    LOOP{
      vanal[1][1][k][j][i]=0; // v^x
      vanal[1][2][k][j][i]=0; // v^y
      vanal[1][3][k][j][i]=0; // v^z

      vectorpot[1][k][j][i]=0; // A_x
      vectorpot[2][k][j][i]=0; // A_y
      vectorpot[3][k][j][i]=0; // A_z

      if((i==N1/2)&&(j==N2/2)&&(k==N3/2)){
	sanal[1][k][j][i]=1.0; // rho
	sanal[2][k][j][i]=0.1; // ie den
      }
      else{
	sanal[1][k][j][i]=0.1; // rho
	sanal[2][k][j][i]=0.01; // ie den
      }
      
      sanal[3][k][j][i]=0.0; // grav. pot. (could try using GRAVACC)

      if(GRAVACC&&(!POSTPROC)){
	for(l=1;l<=3;l++){
	  gravacc[l][k][j][i]=0.0;
	}
      }

    }
#endif

    fprintf(log_file,"proc: %d :  done reading velocity, vector potential, density, and grav. pot file\n",myid); fflush(log_file);






    // END WORRIES









    ///////////////////////
    //
    // Use vector potential to define magnetic field
    //
    ////////////////////////





    // set field
    LOOPN3 LOOPN2 LOOPHP1{
      vanal[2][1][k][j][i]=curlvbacknat1(vectorpot,k,j,i);
    }
    LOOPN3 LOOPHP2 LOOPN1{
      vanal[2][2][k][j][i]=curlvbacknat2(vectorpot,k,j,i);
    }
    LOOPHP3 LOOPN2 LOOPN1{
      vanal[2][3][k][j][i]=curlvbacknat3(vectorpot,k,j,i);
    }


    
    


    
    // output some interesting data
    // and setup initial problem stuff
    if(firstsolve==1){
      firstsolve=0;
      if(myid<=0){
	fprintf(analyticout,"DENSITYFLOOR: %15.10g IEFLOOR: %15.10g\n",DENSITYFLOOR,IEFLOOR);
	if(calltype!=2){
	  fclose(analyticout);
	}
      }// end if cpu write
    }//end if firstsolve
  }
  












  //////////////////////////////////////////////
  //
  //////////////////// VISUALIZATION settings
  //
  //////////////////////////////////////////////
  for(outtype=0;outtype<=1;outtype++){
    mms[outtype][0][1][0]=           9.9e-06 ;
    mms[outtype][1][1][0]=           9.9e-06 ;
    mms[outtype][0][1][1]=    1.0 ;
    mms[outtype][1][1][1]=    1.0 ;
    mms[outtype][0][2][0]=  1E-11 ;
    mms[outtype][1][2][0]=  1E-11 ;
    mms[outtype][0][2][1]=  0.004 ;
    mms[outtype][1][2][1]=  0.004 ;
    mms[outtype][0][3][0]=   -2.56;
    mms[outtype][1][3][0]=   -2.56;
    mms[outtype][0][3][1]=  -0.04818533834 ;
    mms[outtype][1][3][1]=  -0.04818533834 ;  
    mmv[outtype][0][1][1][0]=    -2.5 ;
    mmv[outtype][1][1][1][0]=    -2.5 ;
    mmv[outtype][0][1][1][1]=    1.0 ;
    mmv[outtype][1][1][1][1]=    1.0 ;
    mmv[outtype][0][1][2][0]=   -1.2 ;
    mmv[outtype][1][1][2][0]=   -1.2 ;
    mmv[outtype][0][1][2][1]=   1.2 ;
    mmv[outtype][1][1][2][1]=   1.2 ;
    mmv[outtype][0][1][3][0]=   -1.5 ;
    mmv[outtype][1][1][3][0]=   -1.5 ;
    mmv[outtype][0][1][3][1]=     1.5 ;
    mmv[outtype][1][1][3][1]=     1.5 ;
    mmv[outtype][0][2][1][0]=  -0.05 ;
    mmv[outtype][1][2][1][0]=  -0.05 ;
    mmv[outtype][0][2][1][1]=   0.05 ;
    mmv[outtype][1][2][1][1]=   0.05 ;
    mmv[outtype][0][2][2][0]= -0.02 ;
    mmv[outtype][1][2][2][0]= -0.02 ;
    mmv[outtype][0][2][2][1]= 0.02 ;
    mmv[outtype][1][2][2][1]= 0.02 ;
    mmv[outtype][0][2][3][0]=  -0.06 ;
    mmv[outtype][1][2][3][0]=  -0.06 ;
    mmv[outtype][0][2][3][1]=   0.06 ;
    mmv[outtype][1][2][3][1]=   0.06 ;
  }
  // define outer region when interpolation is used.
  // same order as scalar/vector arrays

  // for images
  for(i=0;i<ITYPES;i++){ // both views
    for(j=0;j<CTYPES;j++){ // both comps
      outerdefs[i][j][1]=mms[i][j][1][0]; // rho
      outerdefs[i][j][2]=mms[i][j][2][0]; // en
      outerdefs[i][j][3]=mms[i][j][3][0]; // pot
      
      outerdefv[i][j][1][0]=0.0; // magnitude of v
      outerdefv[i][j][1][1]=0.0; // v1
      outerdefv[i][j][1][2]=0.0; // v2
      outerdefv[i][j][1][3]=0.0; // v3
      
      outerdefv[i][j][2][0]=0.0;
      outerdefv[i][j][2][1]=0.0;
      outerdefv[i][j][2][2]=0.0;
      outerdefv[i][j][2][3]=0.0;
      
    }
  }
  // for dumps(never use 2nd comp type, and always same for any view)
  douterdefs[1]=mms[0][0][1][0]; // rho
  douterdefs[2]=mms[0][0][2][0]; // en
  douterdefs[3]=mms[0][0][3][0]; // pot

  douterdefv[1][0]=0.0; // magnitude of v
  douterdefv[1][1]=0.0; // v1
  douterdefv[1][2]=0.0; // v2
  douterdefv[1][3]=0.0; // v3

  douterdefv[2][0]=0.0;
  douterdefv[2][1]=0.0;
  douterdefv[2][2]=0.0;
  douterdefv[2][3]=0.0;

}








void averywind(int calltype)
{
  int outtype;
  SFTYPE rho0,bx0;
  SFTYPE divb1avg,divb2avg;
  SFTYPE divb1max,divb2max;
  SFTYPE divb1avg_full,divb2avg_full;
  SFTYPE divb1max_full,divb2max_full;
  char check[50];
  int i,j,k,l,m;
  SFTYPE temp,temp2;
  SFTYPE radsize2,tempf1,tempf2,xz0;
  char temps[100];
  char filename[100];
  static FILE*analyticout;
  FILE*velinput;
  FTYPE q,rstar,MASSSTAR,alpha;

  FILE*fp; //TONO
  // static FILE *vels //TONO

  static int firstsolve=1;
  SFTYPE ftemp1,ftemp2;
  int dumi[12];
  FTYPE dumf[13];
  FTYPE (*vectorpot)[N3M][N2M][N1M]; // vector potential
  FTYPE (*sca3)[N3M][N2M][N1M];
  FTYPE (*vx3)[N3M][N2M][N1M];
  FTYPE (*vy3)[N3M][N2M][N1M];
  FTYPE (*vz3)[N3M][N2M][N1M];

  FTYPE (*radiustouse)[N2M][N1M];
  FTYPE (*thetatouse)[N2M][N1M];
  FTYPE (*phitouse)[N2M][N1M];



  double myx,myy,myz,myr;


  double NUMBONDI,G,C,kbT,mpcsq,Msolar,MASSBH,vgas,rhogas,LENGTHUNIT,TIMEUNIT,MASSUNIT,GRAVACCUNIT,cs,MachNum,vnorm,rhonorm,GMnorm,Rbondi,Rbondicode;

  double vsq;
  int find_ijk_from_x(FTYPE *xyz, int *ijk);
  int ijk_accretor[4];
  FTYPE xtilde,ytilde,ztilde,rtilde;
  FTYPE rtildesoft,Rtilde;
  FTYPE send;
  FTYPE DX_accretorpdir[4],DX_accretormdir[4]; // 3 directions in 1,2,3 indicies
  int dir,ioffset;



  /////////////////////////////////////////
  //
  // Units Setup
  //
  /////////////////////////////////////////
  
  dopassive=0;

  //  DTLOWEST=1E-8;
  //  IDTLOWEST=(1.0/DTLOWEST);
  //  SQIDTLOWEST=(1.0/(DTLOWEST*DTLOWEST)); // actually dt over cour then squared
  //  DTLOWDUMP=(DTLOWEST); // lowest dt to initiate dumping


  passivepos=CENT;
  //  passivepos=VDIR1;
  //  passivepos=-1; // all VDIR1,2,3

  NUMBONDI=15.0; // #  of Bondi radii across box

  G = 6.672E-8; // cgs
  C = 2.99792458E10; // cgs
  kbT = 1.0 ; // keV
  mpcsq = 938375.0; // keV
    
  Msolar = 1.989E33; // cgs
  MASSBH = 10.0*Msolar; // cgs
  MASSSTAR=10.0*Msolar; // cgs
  vgas = 1E3*1E5; // cgs
  rhogas = 1E-2; // cgs
  // 1keV = 1.1605E7Kelvin
  gam=1.4;

  // cs^2 = \gamma P/\rho = \gamma(\gamma-1)u/\rho
  cs = sqrt(gam*kbT/mpcsq)*C; // cgs
  MachNum=vgas/cs; // dimensionless    
  Rbondi=G*MASSBH/pow(vgas,2.0); // cgs

  LENGTHUNIT = Rbondi*NUMBONDI; // cgs
  TIMEUNIT = 1.0/sqrt(G*(MASSBH/pow(LENGTHUNIT,3.0)+rhogas)); // cgs
  MASSUNIT = rhogas*pow(LENGTHUNIT,3.0); // cgs

  Rbondicode=Rbondi/LENGTHUNIT;
  GMnorm=G*MASSBH*pow(TIMEUNIT,2.0)/pow(LENGTHUNIT,3.0); // GM/r^2 in normalized units
  vnorm=vgas/(LENGTHUNIT/TIMEUNIT); // velocity in normalized units
  rhonorm = rhogas/(MASSUNIT/pow(LENGTHUNIT,3.0)); // mass-density in normalized units


  // ASSUMPTIONS:
  // 1) pos_accretor[1,2,3]=0
  // 2) R_NS vector minus R_BH vector is along x-direction always

  // setup rotating system
  // set ROTATINGFORCE=1 in global.h and note how used in stepgen.c
  // distance from accretor to star
  rstar = 1E12/LENGTHUNIT; //normalized
  // Mass ratio
  q = MASSSTAR/MASSBH;
  // distance from accretor to center of mass
  rcm = rstar*q/(1.0+q) ; // normalized
  // Omega of system in normalized units
#if(1)
  Omegasystem = sqrt(GMnorm*(1.0+q)/pow(rstar,3.0));
#else
  Omegasystem = 0.0;
#endif
  alpha=2.0;


  //  fprintf(log_file,"Omegasystem=%21.15g\n",Omegasystem); fflush(log_file);



  // setup temporary memory space

  vectorpot=workv1;
  sca3=workv2;
  vx3=workv3;
  vy3=workv4;
  vz3=workv5;


  /////////////////////////////////////////
  //
  // Physics and code setup
  //
  /////////////////////////////////////////



  if(calltype==100){
    if(POSTPROC==1){
      DYNAMICMM=0; // already read in file usually
    }
    else{
      //DYNAMICMM=0; // use data below
      DYNAMICMM=2; // testing
    }
    

    //    gam=1.4;
    //    cs=1.0;
    rgp=0.0;

    invsol2=invsol2_original=0.0; // 1.0/pow(5.0,2);

    // BAROTROPIC:
    //    transiex1=transiex2=transiex3=0;
    // ADIABATIC:
    transiex1=transiex2=transiex3=1;


    ie=1;
    vischeat=resheat=1;
    
    mag=0; // whether magnetic field is evolved
    visc_real=0;
    res_real=1; rreal=4; resheat= 0; // resistivity
    resist_real0=.01; // from vortex tests
    //res_real= 1;  rreal=2; resheat= 1; resist_real0=0.1; // resistivity
    //res_real=1;


    advint=2; // try WOODWARD slope limiter


    if(dopassive){
      if(mag==0){
	transpassivex1=transpassivex2=transpassivex3=1;
	if(passivepos==CENT){
	  // passive scalar (only works without magnetic field mag=0)
	
	}
	else if(passivepos==VDIR1){
	  transpassivev1x1=1;
	  transpassivev1x2=1;
	  transpassivev1x3=1;
	}
	else if(passivepos==VDIR2){
	  transpassivev2x1=1;
	  transpassivev2x2=1;
	  transpassivev2x3=1;
	}
	else if(passivepos==VDIR3){
	  transpassivev3x1=1;
	  transpassivev3x2=1;
	  transpassivev3x3=1;
	}
	else if(passivepos==-1){
	  transpassivev1x1=1;
	  transpassivev1x2=1;
	  transpassivev1x3=1;
	  transpassivev2x1=1;
	  transpassivev2x2=1;
	  transpassivev2x3=1;
	  transpassivev3x1=1;
	  transpassivev3x2=1;
	  transpassivev3x3=1;
	}
      }
    }
    else{
      transpassivex1=0;
      transpassivex2=0;
      transpassivex3=0;

      transpassivev1x1=0;
      transpassivev1x2=0;
      transpassivev1x3=0;
      transpassivev2x1=0;
      transpassivev2x2=0;
      transpassivev2x3=0;
      transpassivev3x1=0;
      transpassivev3x2=0;
      transpassivev3x3=0;
    }

    if(COORD==1){
      x1in=-0.25;
      x1out=0.75;
      x2in=-0.5;
      x2out=0.5;
      x3in=-0.5;
      x3out=0.5;

#if(0)
      // old uniform grid
      nonunigridx1= 0; // 0=unigrid 1+=nonunigrid
      nonunigridx2= 0; // 0=unigrid 1+=nonunigrid
      nonunigridx3= 0;
#else
      // make sure DOPARDIAG==1 in gpar.c
      // new bi-log grid
      nonunigridx1= 7;
      nonunigridx2= 7;
      nonunigridx3= 7;
#endif


      // only applies if BOUNDTYPE==1
      simplebc=1;
      bcix1=4;
      bcox1=4;
      bcix2=4;
      bcox2=4;
      bcix3=4;
      bcox3=4;

    }
    else{
      fprintf(fail_file,"change to coord=1\n");
      myexit(1);
    }
    //    cour=0.25;
    //    tf=10.0;
    tf = 10.0;
    //    tf=2000.0;
    DTl = tf/1000.0;
    // strictly for purposes of reentrant ability, should have dump/floor same DT
    DTd    = tf/200.0 ;      /* dumping period(data and analytic)*/
    DTfloor=DTd*10.0 ;      /* dumping period(data and analytic)*/
    // strictly for reentract ability, DTi should be integral number of DTds.
    DTpd = 1.0*DTd; // always reentrant data
    DTi    = DTl*5.0 ;        /* image & fieldline(Ax3) period */
    DTfld=DTi*2.0;
    //below not restricted by DTl calls, just each self
    DTener = DTl*1.0 ; // negative means do every time step (must be multiple of DTl)
    // f2's FFT shows could do 20.0 here(was 4.0)
    DTloss = DTl*1.0 ; // negative means do every time step (must be multiple of DTl
    DTmode = DTl*500.0;

    DTtimestep = DTl*100.0 ; // how often to output dominate term in timestep to failfile
    DTsp = DTl*500.0 ; // how often to output sonic point info

    DTtimescale=DTtimestep;
    DTdivb=DTl*10.0;

    // debug mode
    //DTd=DTpd=DTtimestep=tf/100.0;


    //////////////////////////////
    //
    // define accretor properties
    // velocity into horizon of accretor
    VHORIZON[1]=0.1;
    VHORIZON[2]=0.1;
    VHORIZON[3]=0.1;
    
    // overrides VHORIZON above.  Mach number at horizon per component (so real velocity even higher).
    MACHHORIZON=1.5;
    //    MACHHORIZON=-1;
    // if MACHHORIZON<0, then uses VHORIZON

    // velocity of accretor itself
    VACCRETOR[1]=0.0;
    VACCRETOR[2]=0.0;
    VACCRETOR[3]=0.0;

    
    // initial position of accretor
    pos_accretor[1] = 0.0;
    pos_accretor[2] = 0.0;
    pos_accretor[3] = 0.0;
    
    radialdist_accretor = sqrt(pow(pos_accretor[1],2.0)+pow(pos_accretor[2],2.0)+pow(pos_accretor[3],2.0));

#if(1)
    NUMZONE_accretor=3; // no asymmetries and no leaking of passive scalars in CENT,VDIR1,2,3
    //NUMZONE_accretor=2; // leads to asymmetries
#else
    NUMZONE_accretor=4; // try more?
#endif



  }
  else if( (calltype==0)||(calltype==2)){ // normal or injection called it


#define MAX(a,b) (((a)>(b)) ? (a) : (b))

        
    // want to write some interesting data on solution
    if(firstsolve==1){
      fprintf(log_file,"averywind\n"); fflush(log_file);
      if(myid<=0){
	fprintf(logfull_file,"averywind\n"); fflush(logfull_file);
	strcpy(temps,DATADIR);
	
	sprintf(filename,"%s0_analdata%s",temps,DAT2EXT);
	if((analyticout=fopen(filename,"wt"))==NULL){
	  fprintf(fail_file,"Cannot open %s\n",filename);
	  myexit(1);
	}
      }//end if cpu write
    }
    
    
    DENSITYFLOOR=1.0E-10;
    IEFLOOR=1.0E-12;
    DENSITYFRACFLOOR=IEFRACFLOOR=1E-10;





    if(1||DO_ACCRETOR){
      // get accretor size
      if(find_ijk_from_x(pos_accretor, ijk_accretor)){
	// if have accretor on grid, then find its size
	DX_accretor=0.0;
	for(dir=1;dir<=3;dir++){
	  DX_accretorpdir[dir]=DX_accretormdir[dir]=0.0;
	  for(ioffset=0;ioffset<NUMZONE_accretor;ioffset++){
	    DX_accretorpdir[dir]+=dx[1][dir][ijk_accretor[dir]+ioffset];
	    DX_accretormdir[dir]+=dx[1][dir][ijk_accretor[dir]-ioffset];
	  }
	  DX_accretor=MAX(MAX(DX_accretor,DX_accretorpdir[dir]),DX_accretormdir[dir]);
	}
	DX_single=DX_accretor/((FTYPE)NUMZONE_accretor);

	// SUPERGODMARK: This only works if accretor on 3D center.  Need 2.0001 for velocity when accretor is at cell corner
	//	DX_single = 1.0001*MAX(MAX(dx[1][1][ijk_accretor[1]+NUMZONE_accretor],dx[1][2][ijk_accretor[2]]),dx[1][3][ijk_accretor[3]]);
	//	DX_accretor=((FTYPE)NUMZONE_accretor)*DX_single;
      }
      else{
	DX_accretor=DX_single=0.0;
	// then have to get accretor from another CPU
      }
      
      for(dir=1;dir<=3;dir++) fprintf(log_file,"proc: %d DX_accretormdir[%d]=%21.15g DX_accretorpdir[%d]=%21.15g\n",myid,dir,DX_accretormdir[dir],dir,DX_accretorpdir[dir]);
      fprintf(log_file,"proc: %d DX_accretor=%21.15g DX_single=%21.15g\n",myid,DX_accretor,DX_single);

#if(USEMPI)
      // assume maximum is ok
      send=DX_accretor;
      MPI_Allreduce(&send, &DX_accretor, 1, MPI_FTYPE, MPI_MAX, MPI_COMM_WORLD);

      send=DX_single;
      MPI_Allreduce(&send, &DX_single, 1, MPI_FTYPE, MPI_MAX, MPI_COMM_WORLD);
#endif
    }
    else{
      // then assume not doing soft accretor, use hard potential
      DX_single=DX_accretor=0.0;
    }
    ////////////////// GODMARK
    //
    // no check to see if accretor is on full grid!
    //
    //////////////////////////////////////////////
    


    if(myid<=0) {
      fprintf(analyticout,"LENGTHUNIT: %15.10g TIMEUNIT: %15.10g MASSUNIT: %15.10g\n",LENGTHUNIT,TIMEUNIT,MASSUNIT);
      fprintf(analyticout,"GMnorm: %15.10g vnorm: %15.10g rhonorm: %15.10g\n",GMnorm,vnorm,rhonorm);
      fprintf(analyticout,"cs(cgs): %15.10g MachNum %15.10g\n",cs,MachNum);
      // cgs per cgs
      fprintf(analyticout,"u/(\\rho\\Phi): %15.10g   \\rho vgas^2/(\\rho\\Phi) %15.10g\n",(cs*cs*rhogas/(gam*(gam-1.0)))/(rhogas*G*MASSBH/LENGTHUNIT),rhogas*pow(vgas,2.0)/(rhogas*G*MASSBH/LENGTHUNIT));
      fprintf(analyticout,"Rbondi/BOXLENGTH = %15.10g   Rbondi/GRIDLENGTH = %15.10g\n",Rbondicode,Rbondicode/(1.0/(double)N1));
      for(i=1;i<=3;i++) fprintf(analyticout,"ijk_accretor[%d]=%d\n",i,ijk_accretor[i]);
      for(i=1;i<=3;i++) fprintf(analyticout,"pos_accretor[%d]=%g : %g [cgs]\n",i,pos_accretor[i],pos_accretor[i]*LENGTHUNIT);   
      fprintf(analyticout,"DX_single=%g : %g [cgs]\n",DX_single,DX_single*LENGTHUNIT);
      fprintf(analyticout,"NUMZONE_accretor=%d\n",NUMZONE_accretor);
      fprintf(analyticout,"DX_single=%g : DX_accretor=%g : %g [cgs]\n",DX_single,DX_accretor,DX_accretor*LENGTHUNIT);
      fprintf(analyticout,"radialdist_accretor=%g : %g [cgs]\n",radialdist_accretor,radialdist_accretor*LENGTHUNIT);
      fprintf(analyticout,"rstar=%g q=%g rcm=%g Omegasystem=%g alpha=%g [all code units]\n",rstar,q,rcm,Omegasystem,alpha);
    }








    /////////////////////////////////////////
    //
    // Solution Setup
    //
    /////////////////////////////////////////


    // need to define in the code
    // set means you have to read it in

    // set: sanal[1][k][j][i] = mass density
    // set: sanal[2][k][j][i] = internal energy density
    // sanal[3][k][j][i] = 0


    // set: vanal[1][1][k][j][i] = vx
    // set: vanal[1][2][k][j][i] = vy
    // set: vanal[1][3][k][j][i] = vz


    // set: vectorpot[1,2,3][k][j][i]; // A_x A_y A_z

    // vanal[2][1][k][j][i] = Bx
    // vanal[2][2][k][j][i] = By
    // vanal[2][3][k][j][i] = Bz


    // set: gravacc[1,2,3][k][j][i]= gravitational acc

    fprintf(log_file,"proc: %d :  opening inits file\n",myid); fflush(log_file);



    // bounds vector pot
    //    bound(NULL,vectorpot,0,-10-2,123); // assumes all field components computed from some primitive function independent of any field component (no need since last one and gets bounded in a moment)















    
    grids_cart2spc(sca3,vx3,vy3,vz3);
    // vx3[1] : radius x-face
    // vy3[1] : radius y-face
    // vz3[1] : radius z-face





    ////////////////////////////////
    //
    // BELOW IS ALL NORMALIZED QUANTITIES
    //
    ////////////////////////////////

    // LOOP over all grid
    LOOPF{




#if(ROTATINGFORCE)
      /////////////////////////
      //
      // VELOCITIES

      
      
      // X-FACE for x-dir
      xtilde=x[1][1][i]-(pos_accretor[1]-rstar);
      ytilde=x[2][2][j]-(pos_accretor[2]);
      ztilde=x[2][3][k]-(pos_accretor[3]);
      // radius vector to gravitational pit
      rtilde=sqrt( pow(xtilde,2.0) + pow(ytilde,2.0) + pow(ztilde,2.0) );

      // gas velocity in black hole frame  = vg(instantaneous inertial stellar frame) - \Omega\times (\vec{r} - \vec{rstar})
      vanal[1][1][k][j][i]=vnorm*pow(rtilde/(rstar-Rbondicode),2.0-alpha)*(xtilde/rtilde) ; // v^x // gas in star's instantaneous inertial frame
      vanal[1][1][k][j][i]+= ytilde*Omegasystem;



      // Y-FACE for y-dir
      xtilde=x[2][1][i]-(pos_accretor[1]-rstar);
      ytilde=x[1][2][j]-(pos_accretor[2]);
      ztilde=x[2][3][k]-(pos_accretor[3]);
      // radius vector to gravitational pit
      rtilde=sqrt( pow(xtilde,2.0) + pow(ytilde,2.0) + pow(ztilde,2.0) );
      vanal[1][2][k][j][i]=vnorm*pow(rtilde/(rstar-Rbondicode),2.0-alpha)*(ytilde/rtilde) ; // v^y // gas in star's instantaneous inertial frame
      vanal[1][2][k][j][i]+= -xtilde*Omegasystem;


      // Z-FACE for z-dir
      xtilde=x[2][1][i]-(pos_accretor[1]-rstar);
      ytilde=x[2][2][j]-(pos_accretor[2]);
      ztilde=x[1][3][k]-(pos_accretor[3]);
      // radius vector to gravitational pit
      rtilde=sqrt( pow(xtilde,2.0) + pow(ytilde,2.0) + pow(ztilde,2.0) );
      vanal[1][3][k][j][i]=vnorm*pow(rtilde/(rstar-Rbondicode),2.0-alpha)*(ztilde/rtilde) ; // v^z // gas in star's instantaneous inertial frame
      vanal[1][3][k][j][i]+= 0.0;


      //      vsq = vanal[1][1][k][j][i]*vanal[1][1][k][j][i]+vanal[1][2][k][j][i]*vanal[1][2][k][j][i]+vanal[1][3][k][j][i]*vanal[1][3][k][j][i];



      ////////////////////////////////
      //
      // DENSITIES

      // CENT for x-dir acceleration
      xtilde=x[2][1][i]-(pos_accretor[1]-rstar);
      ytilde=x[2][2][j]-pos_accretor[2];
      ztilde=x[2][3][k]-pos_accretor[3];
      // radius vector to gravitational pit
      rtilde=sqrt( pow(xtilde,2.0) + pow(ytilde,2.0) + pow(ztilde,2.0) );


      sanal[1][k][j][i]=rhonorm*pow(rtilde/(rstar-Rbondicode),-alpha); // rho
      sanal[2][k][j][i]=sanal[1][k][j][i]*(vnorm*vnorm)/(gam*(gam-1.0)*MachNum*MachNum); // ie den
      sanal[3][k][j][i]=0.0; // grav. pot. (could try using GRAVACC)



#else


      /////////////////////////
      //
      // VELOCITIES

      
      
      vanal[1][1][k][j][i]=vnorm;
      vanal[1][2][k][j][i]=0;
      vanal[1][3][k][j][i]=0;

      vsq = vanal[1][1][k][j][i]*vanal[1][1][k][j][i]+vanal[1][2][k][j][i]*vanal[1][2][k][j][i]+vanal[1][3][k][j][i]*vanal[1][3][k][j][i];



      ////////////////////////////////
      //
      // DENSITIES

      sanal[1][k][j][i]=rhonorm; // rho
      sanal[2][k][j][i]=sanal[1][k][j][i]*(vsq)/(gam*(gam-1.0)*MachNum*MachNum); // ie den
      sanal[3][k][j][i]=0.0; // grav. pot. (could try using GRAVACC)

#endif



      if(DO_ACCRETOR){


	//////////////////////
	//
	// GRAVITY

	// X-FACE for x-dir acceleration
	xtilde=x[1][1][i]-pos_accretor[1];
	ytilde=x[2][2][j]-pos_accretor[2];
	ztilde=x[2][3][k]-pos_accretor[3];
	// radius vector to gravitational pit
	rtilde=sqrt( pow(xtilde,2.0) + pow(ytilde,2.0) + pow(ztilde,2.0) );


	// X-FACE for x-dir acceleration
	Rtilde=sqrt( pow(xtilde+rstar,2.0) + pow(ytilde,2.0) + pow(ztilde,2.0) );

	

	// soften
	rtildesoft = MAX(rtilde,0.2*DX_single);
	// set
	gravacc[1][k][j][i]=-GMnorm/pow(rtildesoft-rgp,2.0)*xtilde/rtildesoft  - GMnorm*q/pow(Rtilde,2.0)*(xtilde+rstar)/Rtilde;


	// Y-FACE for y-dir acceleration
	xtilde=x[2][1][i]-pos_accretor[1];
	ytilde=x[1][2][j]-pos_accretor[2];
	ztilde=x[2][3][k]-pos_accretor[3];
	// radius vector to gravitational pit
	rtilde=sqrt( pow(xtilde,2.0) + pow(ytilde,2.0) + pow(ztilde,2.0) );


	// Y-FACE for y-dir acceleration
	Rtilde=sqrt( pow(xtilde+rstar,2.0) + pow(ytilde,2.0) + pow(ztilde,2.0) );


	// soften
	rtildesoft = MAX(rtilde,0.2*DX_single);
	// set
	gravacc[2][k][j][i]=-GMnorm/pow(rtildesoft-rgp,2.0)*ytilde/rtildesoft   - GMnorm*q/pow(Rtilde,2.0)*(ytilde)/Rtilde;


	// Z-FACE for z-dir acceleration
	xtilde=x[2][1][i]-pos_accretor[1];
	ytilde=x[2][2][j]-pos_accretor[2];
	ztilde=x[1][3][k]-pos_accretor[3];
	// radius vector to gravitational pit
	rtilde=sqrt( pow(xtilde,2.0) + pow(ytilde,2.0) + pow(ztilde,2.0) );


	// Z-FACE for z-dir acceleration
	Rtilde=sqrt( pow(xtilde+rstar,2.0) + pow(ytilde,2.0) + pow(ztilde,2.0) );

	// soften
	rtildesoft = MAX(rtilde,0.2*DX_single);
	// set
	gravacc[3][k][j][i]=-GMnorm/pow(rtildesoft-rgp,2.0)*ztilde/rtildesoft   - GMnorm*q/pow(Rtilde,2.0)*(ztilde)/Rtilde;

      }
      else{
	gravacc[1][k][j][i]=gravacc[2][k][j][i]=gravacc[3][k][j][i]=0;

      }

    }

    fprintf(log_file,"proc: %d :  done reading velocity, vector potential, density, and grav. pot file\n",myid); fflush(log_file);





    ///////////////////////
    //
    // Use vector potential to define magnetic field
    //
    ////////////////////////

    // initialize vector potential
    LOOPF{
      vectorpot[1][k][j][i]=vectorpot[2][k][j][i]=vectorpot[3][k][j][i]=0;
    }





    // set field
    LOOPN3 LOOPN2 LOOPHP1{
      vanal[2][1][k][j][i]=curlvbacknat1(vectorpot,k,j,i);
    }
    LOOPN3 LOOPHP2 LOOPN1{
      vanal[2][2][k][j][i]=curlvbacknat2(vectorpot,k,j,i);
    }
    LOOPHP3 LOOPN2 LOOPN1{
      vanal[2][3][k][j][i]=curlvbacknat3(vectorpot,k,j,i);
    }




    
    // output some interesting data
    // and setup initial problem stuff
    if(firstsolve==1){
      firstsolve=0;
      if(myid<=0){
	fprintf(analyticout,"DENSITYFLOOR: %15.10g IEFLOOR: %15.10g\n",DENSITYFLOOR,IEFLOOR);
	if(calltype!=2){
	  fclose(analyticout);
	}
      }// end if cpu write
    }//end if firstsolve
  }
  












  //////////////////////////////////////////////
  //
  //////////////////// VISUALIZATION settings
  //
  //////////////////////////////////////////////
  for(outtype=0;outtype<=1;outtype++){
    mms[outtype][0][1][0]=           9.9e-06 ;
    mms[outtype][1][1][0]=           9.9e-06 ;
    mms[outtype][0][1][1]=    1.0 ;
    mms[outtype][1][1][1]=    1.0 ;
    mms[outtype][0][2][0]=  1E-11 ;
    mms[outtype][1][2][0]=  1E-11 ;
    mms[outtype][0][2][1]=  0.004 ;
    mms[outtype][1][2][1]=  0.004 ;
    mms[outtype][0][3][0]=   -2.56;
    mms[outtype][1][3][0]=   -2.56;
    mms[outtype][0][3][1]=  -0.04818533834 ;
    mms[outtype][1][3][1]=  -0.04818533834 ;  
    mmv[outtype][0][1][1][0]=    -2.5 ;
    mmv[outtype][1][1][1][0]=    -2.5 ;
    mmv[outtype][0][1][1][1]=    1.0 ;
    mmv[outtype][1][1][1][1]=    1.0 ;
    mmv[outtype][0][1][2][0]=   -1.2 ;
    mmv[outtype][1][1][2][0]=   -1.2 ;
    mmv[outtype][0][1][2][1]=   1.2 ;
    mmv[outtype][1][1][2][1]=   1.2 ;
    mmv[outtype][0][1][3][0]=   -1.5 ;
    mmv[outtype][1][1][3][0]=   -1.5 ;
    mmv[outtype][0][1][3][1]=     1.5 ;
    mmv[outtype][1][1][3][1]=     1.5 ;
    mmv[outtype][0][2][1][0]=  -0.05 ;
    mmv[outtype][1][2][1][0]=  -0.05 ;
    mmv[outtype][0][2][1][1]=   0.05 ;
    mmv[outtype][1][2][1][1]=   0.05 ;
    mmv[outtype][0][2][2][0]= -0.02 ;
    mmv[outtype][1][2][2][0]= -0.02 ;
    mmv[outtype][0][2][2][1]= 0.02 ;
    mmv[outtype][1][2][2][1]= 0.02 ;
    mmv[outtype][0][2][3][0]=  -0.06 ;
    mmv[outtype][1][2][3][0]=  -0.06 ;
    mmv[outtype][0][2][3][1]=   0.06 ;
    mmv[outtype][1][2][3][1]=   0.06 ;
  }
  // define outer region when interpolation is used.
  // same order as scalar/vector arrays

  // for images
  for(i=0;i<ITYPES;i++){ // both views
    for(j=0;j<CTYPES;j++){ // both comps
      outerdefs[i][j][1]=mms[i][j][1][0]; // rho
      outerdefs[i][j][2]=mms[i][j][2][0]; // en
      outerdefs[i][j][3]=mms[i][j][3][0]; // pot
      
      outerdefv[i][j][1][0]=0.0; // magnitude of v
      outerdefv[i][j][1][1]=0.0; // v1
      outerdefv[i][j][1][2]=0.0; // v2
      outerdefv[i][j][1][3]=0.0; // v3
      
      outerdefv[i][j][2][0]=0.0;
      outerdefv[i][j][2][1]=0.0;
      outerdefv[i][j][2][2]=0.0;
      outerdefv[i][j][2][3]=0.0;
      
    }
  }
  // for dumps(never use 2nd comp type, and always same for any view)
  douterdefs[1]=mms[0][0][1][0]; // rho
  douterdefs[2]=mms[0][0][2][0]; // en
  douterdefs[3]=mms[0][0][3][0]; // pot

  douterdefv[1][0]=0.0; // magnitude of v
  douterdefv[1][1]=0.0; // v1
  douterdefv[1][2]=0.0; // v2
  douterdefv[1][3]=0.0; // v3

  douterdefv[2][0]=0.0;
  douterdefv[2][1]=0.0;
  douterdefv[2][2]=0.0;
  douterdefv[2][3]=0.0;

}


// accretor on corner assumed.
// that is, find ijk for face center x,y,z positions
int find_ijk_from_x(FTYPE *xyz, int *ijk)
{
  int i,j,k;

  // see if on grid
  if(x[1][1][-N1BND]>xyz[1]) return(0);
  if(x[1][1][N1-1-N1BND-1]<xyz[1]) return(0);

  if(x[1][2][-N2BND]>xyz[2]) return(0);
  if(x[1][2][N2-1-N2BND-1]<xyz[2]) return(0);

  if(x[1][3][-N3BND]>xyz[3]) return(0);
  if(x[1][3][N3-1-N3BND-1]<xyz[3]) return(0);

  // if on grid, then locate i,j,k

  LOOPF1 if(x[1][1][i]>xyz[1]) break;
  ijk[1]=i-1;

  LOOPF2 if(x[1][2][j]>xyz[2]) break;
  ijk[2]=j-1;

  LOOPF3 if(x[1][3][k]>xyz[3]) break;
  ijk[3]=k-1;

  return(1);

}










void bondiplustorus(int calltype)
{
  int i,j,k;
  static FILE*analyticout;


  if(calltype==100){
    tori1sol(calltype); // uess torus problem setup
  }
  else if( (calltype==0) ){
    fprintf(log_file,"Bondiplustori before tori1sol(2)\n"); fflush(log_file);
    analyticout=tori1sol(2); // submit for solution and file setup (assumes atmosphere is floor value)
    fprintf(log_file,"Bondiplustori after tori1sol(2)\n"); fflush(log_file);

    if((runtype==0)||(runtype==2)){ // needed since using s and gets restart can't be overwritten, don't need torus on outer edge anyways
      fprintf(log_file,"begin transfer\n"); fflush(log_file);
      // transfer to holding cell
      LOOPF{
	s[1][k][j][i]=sanal[1][k][j][i];
	s[2][k][j][i]=sanal[2][k][j][i];
	s[3][k][j][i]=sanal[3][k][j][i];
	v[1][1][k][j][i]=vanal[1][1][k][j][i];
	v[1][2][k][j][i]=vanal[1][2][k][j][i];
	v[1][3][k][j][i]=vanal[1][3][k][j][i];
	v[2][1][k][j][i]=vanal[2][1][k][j][i];
	v[2][2][k][j][i]=vanal[2][2][k][j][i];
	v[2][3][k][j][i]=vanal[2][3][k][j][i];
      }
      
      fprintf(log_file,"end transfer\n"); fflush(log_file);
    }

    fprintf(log_file,"begin bondi\n"); fflush(log_file);
    bondisol(2,analyticout); // submit for solution and file setup
    fprintf(log_file,"end bondi\n"); fflush(log_file);

    if((runtype==0)||(runtype==2)){ // needed since using s and gets restart can't be overwritten, don't need torus on outer edge anyways
      fprintf(log_file,"begin retransfer\n"); fflush(log_file);
      LOOPF{
	if((sanal[1][k][j][i]>s[1][k][j][i])||(sanal[2][k][j][i]>s[2][k][j][i])){
	  // second line gets rid of little spikes in ie, first in rho
	  // bondi is bottom in mass density
	  // the above line isn't necessarily general, requires fine tuning.
	  // fill with bondi solution
	  s[1][k][j][i]=sanal[1][k][j][i];
	  s[2][k][j][i]=sanal[2][k][j][i];
	  s[3][k][j][i]=sanal[3][k][j][i];
	  v[1][1][k][j][i]=vanal[1][1][k][j][i];
	  v[1][2][k][j][i]=vanal[1][2][k][j][i];
	  v[1][3][k][j][i]=vanal[1][3][k][j][i];
	  //v[2][1][k][j][i]=vanal[2][1][k][j][i];
	  //v[2][2][k][j][i]=vanal[2][2][k][j][i];
	  //v[2][3][k][j][i]=vanal[2][3][k][j][i];
	}
	else{
	  // leave alone
	  // leave B alone so divB is right
	}

	// finally assign value
	sanal[1][k][j][i]=s[1][k][j][i];
	sanal[2][k][j][i]=s[2][k][j][i];
	sanal[3][k][j][i]=s[3][k][j][i];
	vanal[1][1][k][j][i]=v[1][1][k][j][i];
	vanal[1][2][k][j][i]=v[1][2][k][j][i];
	vanal[1][3][k][j][i]=v[1][3][k][j][i];
	vanal[2][1][k][j][i]=v[2][1][k][j][i];
	vanal[2][2][k][j][i]=v[2][2][k][j][i];
	vanal[2][3][k][j][i]=v[2][3][k][j][i];
      }
    }
    fprintf(log_file,"end retransfer\n"); fflush(log_file);
    
    if(myid<=0) fclose(analyticout);
  }


}


FILE * diskini(int calltype)
{
  FTYPE rstar,beta,alpha,h0,hei,dm,tcritcgs;
  FTYPE opt,surf,sigtemp,kaptemp,dentemp,ttemp,opatemp;
  SFTYPE numvtime,tvisc;
  int outtype;
  SFTYPE rho;
  char check[50];
  int i,j,k,l,m,n,itest,jtest;
  FILE *testfp;
  SFTYPE temp,temp2;
  SFTYPE radsize2,tempf1,tempf2,xz0;
  char temps[100];
  char filename[100];
  SFTYPE theta0;
  SFTYPE Rout,l0sq;
  SFTYPE CFRACT;
  static FILE*analyticout;
  static int firstsolve=1;
  SFTYPE totalmass[2]; // 0: tori 1: atmosphere
  SFTYPE totalmass_full[2]; // 0: tori 1: atmosphere
  // whether tori is at a point in the grid or not
  int torigrid;
  static int ATMONLYT, ADATM;
  static int RANDOMIZE,VRANDOMIZE;
  SFTYPE MODEPURT,MODEAMP_rho,MODEAMP_u,MODEAMP_v1,MODEAMP_v2,MODEAMP_v3;
  SFTYPE randomization;
  SFTYPE vr,vtheta,vphi,br,btheta,bphi;
  FTYPE (*sca3)[N3M][N2M][N1M];
  FTYPE (*vx3)[N3M][N2M][N1M];
  FTYPE (*vy3)[N3M][N2M][N1M];
  FTYPE (*vz3)[N3M][N2M][N1M];

  FTYPE (*radiustouse)[N2M][N1M];
  FTYPE (*thetatouse)[N2M][N1M];
  FTYPE (*phitouse)[N2M][N1M];
  SFTYPE ftemp1,ftemp2;
  int itemp;
  FTYPE gravpottemp;
  int component;
  SFTYPE vrand;
  SFTYPE KAPPA;
  SFTYPE SMOOTHSIZE_rho,SMOOTHSIZE_u,SMOOTHSIZE_v1,SMOOTHSIZE_v2,SMOOTHSIZE_v3;
  SFTYPE ATMFACTOR;
  SFTYPE DISKTYPE,Routt;
  SFTYPE shapeouter;


  ATMONLYT= 0; // 0: normal disk 1: just atmosphere
  ATMFACTOR=1.0; // factor above floor given to atmosphere
  ADATM= 0; // 1: adiabatic atm 0: fraction pot atm

  RANDOMIZE=0; // whether any randomization at all
  VRANDOMIZE=2; // 1==vphi only for velocity, and density, 2=all (see below) 3=m=1 density pert

  // for VRANDOMIZE==1 or 2
  //
  vrand=1.0; // scale for randomizing velocities
  randomization=5.E-2; //  (at N3=32)

  // FOR VRANDOMIZE==3
  //
  MODEPURT=4.0; // mode to purtube for VRANDOMIZE==3

  MODEAMP_rho=3.16E-5; // mode amplitude "" as fraction of density maximum
  MODEAMP_u=MODEAMP_rho*0.4/0.2; // mode amplitude "" as fraction of density maximum

  MODEAMP_v1=MODEAMP_rho*.01/0.2; // mode amplitude "" as fraction of density maximum
  MODEAMP_v2=MODEAMP_rho*.01/0.2; // mode amplitude "" as fraction of density maximum
  MODEAMP_v3=MODEAMP_rho*.01/0.2; // mode amplitude "" as fraction of density maximum


  // SMOOTHING
  //
  SMOOTHSIZE_rho=0; // kernel around zone to combine in smoothing, 0=no smoothing
  SMOOTHSIZE_u=0; // kernel around zone to combine in smoothing, 0=no smoothing
  SMOOTHSIZE_v1=0; // kernel around zone to combine in smoothing, 0=no smoothing
  SMOOTHSIZE_v2=0; // kernel around zone to combine in smoothing, 0=no smoothing
  SMOOTHSIZE_v3=0; // kernel around zone to combine in smoothing, 0=no smoothing

  if(calltype==100){
    if(POSTPROC==1){
      DYNAMICMM=0; 
      // already read in file usually
    }
    else{
      DYNAMICMM=2; // testing
    }
    
    // Newtonian pot
      rg=2; 
      rgp=0;

      visc_real=1;      vreal=4;
      mag=0;      res_real=0;
  	alpha_real0=0.1; 
	alphafrac=1.e19;
      if(COORD==3){
	x3in=  (0.0);
	x3out= (2.0*M_PI);
	
	Rinner=0.07*rg;
	Router=20.*rg; // normal outer boundary HK-like
	
	x1in=Rinner;
	x1out= Router; // H00 first axisym
	
	
        x2in=  (4.*M_PI/12.); // A3
        x2out= (M_PI/2.);
	
	// normal
	nonunigridx1= 5; // 0=unigrid 1+=nonunigrid
	nonunigridx2= 0; // 0=unigrid 1+=nonunigrid
	nonunigridx3= 0;

  
	simplebc=1;
	bcix1=4;
	bcox1=4;
	bcix2=4; // really AOS (see boundrect1.c)
	bcox2=1; // really AOS (see boundrect1.c)
	bcix3=5;
	bcox3=5;

      }
      else if(COORD==1){
      }
      
    
    
      // viscosity run
      tvisc=2.0*M_PI*pow(12.0,1.5)/alpha_real0; //viscous time scales
      numvtime=1.0;
      tf =1.e4; // N viscous time scale

      DTl = 4./30.0; // H00
      // strictly for purposes of reentrant ability, should have dump/floor same DT
      DTd    = DTl*30. ;      /* dumping period(data and analytic)*/
      DTfloor=DTd*10.0 ;      /* dumping period(data and analytic)*/
      // strictly for reentract ability, DTi should be integral number of DTds.
      DTpd = 10.0*DTd; // always reentrant data
      DTi    = DTl*300.0 ;        /* image & fieldline(Ax3) period */
      DTfld=DTi*20.0;
      //below not restricted by DTl calls, just each self
      DTener = DTl*10.0 ; // negative means do every time step (must be multiple of DTl)
      // f2's FFT shows could do 20.0 here(was 4.0)
      DTloss = DTl*10.0 ; // negative means do every time step (must be multiple of DTl
      DTmode = DTl*20.0;
      
      DTtimestep = DTl*10.0 ; // how often to output dominate term in timestep to failfile
      DTdivb=DTsp = DTl*100.0 ; // how often to output sonic point info
   
  }// everything after this is called only after grid is init'ed, etc.
  else if((calltype==0)||(calltype==2)){ // normal or injection called it


    sca3=workv2;
    vx3=workv3;
    vy3=workv4;
    vz3=workv5;

    grids_cart2spc(sca3,vx3,vy3,vz3);
    
    
    // want to write some interesting data on solution
    if(firstsolve==1){
      if(myid<=0){
	fprintf(logfull_file,"disk init\n");
	fflush(logfull_file);
	strcpy(temps,DATADIR);
	
	sprintf(filename,"%s0_analdata%s",temps,DAT2EXT);
	if((analyticout=fopen(filename,"wt"))==NULL){
	  fprintf(fail_file,"Cannot open %s\n",filename);
	  myexit(1);
	}
      }//end if cpu write
    }
    

    theta0=M_PI/2.0;
  
      R0=1.0; // must be 1.0
      GM=1.0;
      Omega0=sqrt(GM/pow(R0,3.0)); // 
      KAPPA=1.0;
      Rint=Rinner*1.1; 
      Routt=Router*0.9;
      mstarcgs=1.989E33 ; // in cgs
// set unit
      munit=5.9736E27; // in cgs
      lunit=1.49598E13; //in cgs
      rhounit=munit/lunit/lunit/lunit;
      tunit=pow(6.67E-8*mstarcgs/lunit/lunit/lunit,-0.5);
      eneunit=munit*lunit/tunit*lunit/tunit;
      enedenunit=eneunit/lunit/lunit/lunit;
// set constant 
      lsuncgs=3.827E33  ; //in cgs 
      kcontcgs=1.381E-16  ; //in cgs
      sunmasscgs=1.989E33 ; // in cgs
      sigmacontcgs=5.67051E-5 ; // in cgs
      yrcgs=365.24*24.*3600 ; // in cgs
      rsun=6.96E10/lunit;
      mmw=2.3  ;
      mhydrcgs=1.67352E-24  ; //in cgs
      tcritcgs=1500.; // in cgs
// set constant in new unit
      sigmacontnew=sigmacontcgs/eneunit*lunit*lunit*tunit;
      kcontnew=kcontcgs/eneunit;
      mhydrnew=mhydrcgs/munit;
// layer modle
      Siga=50/(munit/lunit/lunit);
      actcontour=1;
// cooling
      coolfact=7.;
      cool=1;
// set constant accretion disk
      beta=1.25;
      alpha=2.25;
      lsunfrac=0.5;
      lstarcgs=1.*lsuncgs; //in cgs
      rstar=1.*rsun;     
      h0=0.01*rstar;
      MASSBH=1.*sunmasscgs/munit; 
      GRAVC=GM/MASSBH;
      cstcrit2=gam*kcontcgs*tcritcgs/mmw/mhydrcgs/(lunit/tunit*lunit/tunit);
      dm=1.0E-5*sunmasscgs/yrcgs/(munit/tunit); 
      rho00=dm*rstar/sqrt(18.0*M_PI*M_PI*M_PI)/0.1/sqrt(GM/rstar)/h0/h0/h0;     	
      rho0=rho00*(1.0-sqrt(rstar/Rinner))*pow(rstar/Rinner,alpha);
      u0=1.*rho0*rhounit*kcontcgs/mmw/mhydrcgs/(gam-1.)/enedenunit;
//
#define PHI_G(r) (-GM/((r)-rgp))
#define V_PHII(R) (R0*Omega0*pow(fabs((R)/R0),-0.5))
#define V_R(R) (0)
#define V_THETA(R) (0)


    // used for non adiabatic atm
    IEFRACT=0.2;
    VZFRACT=.95;

      // normal
//      DENSITYFRACFLOOR=1.E-4;
//      IEFRACFLOOR=DENSITYFRACFLOOR;
//      DENSITYFLOOR=DENSITYFRACFLOOR*rho0;
 	DENSITYFLOOR=Siga/Router*1.E-15;
 
	IEFLOOR=u0*DENSITYFLOOR/rho0;
//	IEFLOOR=IEFRACT*fabs(DENSITYFLOOR*PHI_G(fabs(x1out)));



    ////////////////////////////////////////
    ////////////////////////////////////////
    /// SET THE SOLUTION
    ///
    ////////////////////////////////////////
    ////////////////////////////////////////
    
    
    totalmass[0]=0.0;
    totalmass[1]=0.0;

    if(REENTER){reentrancezz();}

    initialize_gravity(sca3,vx3,vy3,vz3); // optimizations
    
    // set the primitive variables, except magnetic field(which is derived from density)
    for(l=1;l<=1+3;l++){ // 1=scalars 2=v_1 3=v_2 4=v_3 , since different positions and want IC to be exactly numerically symmetric

      if(l==1){
	radiustouse=sca3[1];
	thetatouse=sca3[2];
	phitouse=sca3[3];
      }
      else if(l==2){
	radiustouse=vx3[1];
	thetatouse=vx3[2];
	phitouse=vx3[3];
      }
      else if(l==3){
	radiustouse=vy3[1];
	thetatouse=vy3[2];
	phitouse=vy3[3];
      }
      else if(l==4){
	radiustouse=vz3[1];
	thetatouse=vz3[2];
	phitouse=vz3[3];
      }
      fprintf(log_file,"%d here 1\n",l); fflush(log_file);
	
      LOOPF{
	
	torigrid=1;
	gravpottemp=PHI_G(radiustouse[k][j][i]);
/*	if((ATMONLYT==0)&&(fabs(radiustouse[k][j][i]*sin(thetatouse[k][j][i]))>=Rint)&&(fabs(radiustouse[k][j][i]*sin(thetatouse[k][j][i]))<=Routt) ){ // then assume within the torus  //  need since some solutions have solutions beyond tori at inner region
		torigrid=1;
	}
	else torigrid=0;*/
	torigrid=1;
	
	// do scalars
	if(l==1){
	  
	  if(torigrid){// then torus
             hei=h0*pow(radiustouse[k][j][i]*fabs(sin(thetatouse[k][j][i]))/rstar,beta);
	     sanal[1][k][j][i]=rho=rho00*(1.0-sqrt(rstar/fabs(radiustouse[k][j][i])/fabs(sin(thetatouse[k][j][i]))))*pow(rstar/fabs(radiustouse[k][j][i])/fabs(sin(thetatouse[k][j][i])),alpha)*exp(-0.5*fabs(radiustouse[k][j][i])*cos(thetatouse[k][j][i])/hei*fabs(radiustouse[k][j][i])*cos(thetatouse[k][j][i])/hei);
//            printf("%g %g %g %g %g %g\n",radiustouse[k][j][i],thetatouse[k][j][i],fabs(sin(thetatouse[k][j][i])),1.0-sqrt(rstar/fabs(radiustouse[k][j][i])/fabs(sin(thetatouse[k][j][i]))),pow(rstar/fabs(radiustouse[k][j][i])/fabs(sin(thetatouse[k][j][i])),alpha),exp(-0.5*fabs(radiustouse[k][j][i])*cos(thetatouse[k][j][i])/hei*fabs(radiustouse[k][j][i])*cos(thetatouse[k][j][i])/hei));
	    if(RANDOMIZE){ // uncorrelated(not even adiabatic) random perturbations that scale with resolution such that at N3=32 the perturbation is =randomization of the initial value at R0 for rho, en, and vphi
	      if(VRANDOMIZE==1) sanal[1][k][j][i]+=rho0*randomization*sqrt((SFTYPE)(N3)/32.0)*(ranc(0)-.5);
	      else if(VRANDOMIZE==2) sanal[1][k][j][i]*=(1. + 1.e-2*(ranc(0)-0.5)) ;
	      else if(VRANDOMIZE==3) sanal[1][k][j][i]+=rho0*MODEAMP_rho*sin(MODEPURT*phitouse[k][j][i]);
	    }
	    // floor check applies to torus too! and make sure randomize doesn't go below floor
	    if(sanal[1][k][j][i]<DENSITYFLOOR) sanal[1][k][j][i]=DENSITYFLOOR;
	    
	    // mass add up
	    if( (i>=0)&&(i<N1)&&(j>=0)&&(j<N2)&&(k>=0)&&(k<N3) ){// only add that which is on real grid
	      totalmass[1]+=sanal[1][k][j][i]*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k);
	    }
	    
	    // internal energy (uses unperturbed rho)
// temperature times 0.1 to simulate the viscous heating
	      tirr0[k][j][i]=pow(lsunfrac*lstarcgs/4./M_PI/(radiustouse[k][j][i]*lunit*radiustouse[k][j][i]*lunit)/sigmacontcgs,0.25);// means default temperature unit is kelvin
	      tirr0[k][j][i]=30.;
          }
        }
      }	

     LOOPC3{
	LOOPC1{
	  teff4[i]=3.*GM*5.0E-8*sunmasscgs/yrcgs/(munit/tunit)/8./M_PI/sigmacontnew/pow(radiustouse[k][1][i],3.);
	  LOOPC2{
	      tirr0[k][j][i]=pow(teff4[i],0.25);
	      dentemp=sanal[1][k][j][i]*rhounit;
              ttemp=tirr0[k][j][i];
              opatemp=opac(dentemp, ttemp);
	      kaptemp=opatemp/lunit/lunit*munit*sanal[1][k][j][i];
	      sigtemp=0./lunit/lunit*munit*sanal[1][k][j][i];
	      opaini[k][j][i]=kaptemp+sigtemp;
	  }
	}
     }

     bound_mpi(sanal[1],NULL,-10-1,0,0); 
     tempini(0);
     for(irowmpi=0;irowmpi<ncpux2-1;irowmpi++){
      bound1D_mpi(sigini, optdini, 2, irowmpi, 0);
      tempini(1);
     }
         //no tempin(0) again
     for(irowmpi=0;irowmpi<ncpux2-1;irowmpi++){
      bound1D_mpi(tirr0, tirr0, 1, irowmpi, 0);
      tempini(1);
     }
     LOOPC{
       diagn[3][k][j][i]=sanal[1][k][j][i];
       diagn[1][k][j][i]=sigini[k][j][i];
       diagn[2][k][j][i]=tirr0[k][j][i];
     }

     LOOPF{
	torigrid=1;
        gravpottemp=PHI_G(radiustouse[k][j][i]);
        torigrid=1;

	if(l==1){
          if(torigrid){
	      sanal[2][k][j][i] = tirr0[k][j][i]*sanal[1][k][j][i]*rhounit*kcontcgs/mmw/mhydrcgs/(gam-1.)/enedenunit; 
              cs02[k][j][i]=gam*(gam-1.0)*sanal[2][k][j][i]/sanal[1][k][j][i]; 	   

	    if(RANDOMIZE){
	      if(VRANDOMIZE==1) sanal[2][k][j][i]+=u0*randomization*sqrt((SFTYPE)(N3)/32.0)*(ranc(0)-.5);
	      else if(VRANDOMIZE==3) sanal[2][k][j][i]+=u0*MODEAMP_u*sin(MODEPURT*phitouse[k][j][i]);
	    }
	    // floor check applies to torus too!
	    if(sanal[2][k][j][i]<IEFLOOR) sanal[2][k][j][i]=IEFLOOR;
	    
	  }
	  else{ // then atmospher
	      // density
	     
	      sanal[1][k][j][i] = DENSITYFLOOR*ATMFACTOR;
	      sanal[2][k][j][i] = IEFLOOR*ATMFACTOR;

	    if( (i>=0)&&(i<N1)&&(j>=0)&&(j<N2)&&(k>=0)&&(k<N3) ){// only add that which is on real grid
	      totalmass[0]+=sanal[1][k][j][i]*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k);
	    }
	    
	  }// end atmosphere
	  
	  // grav pot
	  sanal[3][k][j][i]=gravpottemp;
	  
	  
	}// end if scalars
	if(l>=2){ // v_1,2,3, using different positions
	  if(l==2) component=1;
	  if(l==3) component=2;
	  if(l==4) component=3;
	  
	  if(torigrid){// torus
	    
	    
	    if(RANDOMIZE){ // scale for randomizations
	      vphi=V_PHII(R0)*randomization*sqrt((SFTYPE)(N3)/32.0)*(ranc(0)-.5);
	    }
	    if(COORD==3){
	      if(component==1){
		vanal[1][1][k][j][i]=0.0;
		if(RANDOMIZE){
		  if(VRANDOMIZE==2) vanal[1][1][k][j][i]+=vrand*randomization*(ranc(0)-0.5) ;
		  else if(VRANDOMIZE==3) vanal[1][1][k][j][i]+=vphi*MODEAMP_v1*sin(MODEPURT*phitouse[k][j][i]);
		}
	      }
	      if(component==2){
		vanal[1][2][k][j][i]=0.0;
		if(RANDOMIZE){
		  if(VRANDOMIZE==2) vanal[1][2][k][j][i]+=vrand*randomization*(ranc(0)-0.5) ;
		  else if(VRANDOMIZE==3) vanal[1][2][k][j][i]+=vphi*MODEAMP_v2*sin(MODEPURT*phitouse[k][j][i]);
		}
	      }
	      if(component==3){
		vanal[1][3][k][j][i]=V_PHII(radiustouse[k][j][i]*sin(thetatouse[k][j][i]));
		if(RANDOMIZE){
		  if(VRANDOMIZE==1) vanal[1][3][k][j][i]+=vphi;
		  else if(VRANDOMIZE==2) vanal[1][3][k][j][i]+=vrand*randomization*(ranc(0)-0.5) ;
		  else if(VRANDOMIZE==3) vanal[1][3][k][j][i]+=vphi*MODEAMP_v3*sin(MODEPURT*phitouse[k][j][i]);
		}
	      }
	    }
	    if(COORD==1){
	    }
	  }
	  else{// atmosphere   
	      if(component==1) vanal[1][1][k][j][i]=0.0;
	      if(component==2) vanal[1][2][k][j][i]=0.0;
	      if(component==3) vanal[1][3][k][j][i]=0.0;
	  }
	}// end if vectors
	
      }// end full loop
      // bound each new variable since successive variables may depend on it (e.g. multiple cpus and random fluctuations in density for field)
      fprintf(log_file,"%d here 2\n",l); fflush(log_file);

      // reenter one dump file

      if(REENTER) {
        LOOPC{
          if(l==1){
	   if(radiustouse[k][j][i]<0.15*rg){
            sanal[1][k][j][i]=1.*extre(radiustouse[k][j][i],thetatouse[k][j][i],1);
            sanal[2][k][j][i]=1.*extre(radiustouse[k][j][i],thetatouse[k][j][i],2);
/*	    if(extre(radiustouse[k][j][i],thetatouse[k][j][i],8)>5.){
		sanal[1][k][j][i]=1.*extre(radiustouse[k][j][i],thetatouse[k][j][i],1);
		sanal[2][k][j][i]=1.*extre(radiustouse[k][j][i],thetatouse[k][j][i],2);
	    }
*/
	   }else{
	    sanal[1][k][j][i]=1.*extre(radiustouse[k][j][i],thetatouse[k][j][i],1);
            sanal[2][k][j][i]=1.*extre(radiustouse[k][j][i],thetatouse[k][j][i],2);
	   }
          }
          if(l==2){vanal[1][1][k][j][i]=extre(radiustouse[k][j][i],thetatouse[k][j][i],4);}
          if(l==3){
		if(radiustouse[k][j][i]>0.5*rg){
		   vanal[1][2][k][j][i]=extre(radiustouse[k][j][i],thetatouse[k][j][i],5);
		}else{
		   vanal[1][2][k][j][i]=extre(radiustouse[k][j][i],thetatouse[k][j][i],5);
		}
	  }
          if(l==4){
		if(radiustouse[k][j][i]>0.5*rg){
		   vanal[1][3][k][j][i]=extre(radiustouse[k][j][i],thetatouse[k][j][i],6);
		}else{
		   vanal[1][3][k][j][i]=extre(radiustouse[k][j][i],thetatouse[k][j][i],6);
		}
          }
        }
//      fprintf(stderr,"radius %10.3lf %10.3lf\n ",radiustouse[0][0][0],thetatouse[0][0][0]);
      }



      if(l==1){
	smooth(SMOOTHSIZE_rho,DENSITYFLOOR,sanal[1]);
	bound(sanal[1],NULL,-10-1,0,0); // assumes all scalars independently computed
	fprintf(log_file," Sanal %15.10g %15.10g %15.10g %15.10g \n",sanal[1][0][10][N1-2],sanal[1][0][10][N1-1],sanal[1][0][10][N1],sanal[1][0][10][N1+1]);
	fprintf(log_file," X %15.10g %15.10g %15.10g %15.10g \n",x[2][1][N1-2],x[2][1][N1-1],x[2][1][N1],x[2][1][N1+1]);
 	fprintf(log_file," Sanal %15.10g %15.10g %15.10g %15.10g \n",sanal[1][0][10][-2],sanal[1][0][10][-1],sanal[1][0][10][0],sanal[1][0][10][1]);
	fprintf(log_file," X %15.10g %15.10g %15.10g %15.10g \n",x[2][1][-2],x[2][1][-1],x[2][1][0],x[2][1][1]);
	smooth(SMOOTHSIZE_u,IEFLOOR,sanal[2]);
	bound(sanal[2],NULL,-10-2,0,0); // assumes all scalars independently computed
	// sanal[3] must be independent and so never needs bounding
      }
      if(l==2){
	smooth(SMOOTHSIZE_v1,0.0,vanal[1][1]);
	bound(NULL,vanal[1],0,-10-1,1);
      }
      if(l==3){
	smooth(SMOOTHSIZE_v2,0.0,vanal[1][2]);
	bound(NULL,vanal[1],0,-10-1,2);
      }
      if(l==4){
	smooth(SMOOTHSIZE_v3,0.0,vanal[1][3]);
	bound(NULL,vanal[1],0,-10-1,3);
      }
      fprintf(log_file,"%d here 2.5\n",l); fflush(log_file);    
    }// end loop over variously centered primitive variables
    fprintf(log_file,"here 3\n"); fflush(log_file);    


    // can't just do tori region to make divb=0 correclty
      testfp=fopen("initialrho.dat","wt");
	fprintf(testfp,"%g %g %g %g %g %g %g %g %g\n",rhounit,kcontcgs,mmw,mhydrcgs,gam,enedenunit,munit,lunit,tunit);
      for(jtest=0;jtest<64;jtest++){
  	for(itest=0;itest<64;itest++){
	  fprintf(testfp,"%g %g %g %g %g %g %g\n",sanal[1][0][jtest][itest],DENSITYFLOOR,tirr0[0][jtest][itest],sanal[2][0][jtest][itest],IEFLOOR,radiustouse[0][jtest][itest],thetatouse[0][jtest][itest]);    
        }
      }
      LOOPF{
	vanal[2][1][k][j][i]=0.0;
	vanal[2][2][k][j][i]=0.0;
	vanal[2][3][k][j][i]=0.0;
      }
  
    // no need to bound since all CPUs the same at this point
    fprintf(log_file,"here 4\n"); fflush(log_file);        

    
    fprintf(log_file,"DENSITYFLOOR: %21.15g IEFLOOR: %21.15g\n",DENSITYFLOOR,IEFLOOR); fflush(log_file);
    fprintf(log_file,"DENSITYFRACFLOOR: %21.15g IEFRACFLOOR: %21.15g\n",DENSITYFRACFLOOR,IEFRACFLOOR); fflush(log_file);
    
    
    
    
    
    // output some interesting data
    // and setup initial problem stuff
    if(firstsolve==1){
      firstsolve=0;
      
      if(numprocs>1){
#if(USEMPI)
	MPI_Reduce(&(totalmass[0]), &(totalmass_full[0]), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
	MPI_Reduce(&(totalmass[1]), &(totalmass_full[1]), 1, MPI_SFTYPE, MPI_SUM, 0,MPI_COMM_WORLD);
#endif
      }
      else{
	totalmass_full[0]=totalmass[0];
	totalmass_full[1]=totalmass[1];
      }
      if(myid<=0){
	// now output some interesting analytic data
	fprintf(analyticout,"rho0: %15.10g R0: %15.10g u0: %15.10g Omega0: %15.10g\n",rho0,R0,u0,Omega0);
	if(RMAX!=0.0) fprintf(analyticout,"RMAX: %15.10g DELTAR: %15.10g HOR: %15.10g\n",RMAX,DELTAR,HOR);
	fprintf(analyticout,"DENSITYFLOOR: %15.10g IEFLOOR: %15.10g\n",DENSITYFLOOR,IEFLOOR);
	fprintf(analyticout,"DENSITYFRACFLOOR: %15.10g IEFRACFLOOR: %15.10g\n",DENSITYFRACFLOOR,IEFRACFLOOR);
	if(1.0*.001<DENSITYFLOOR){
	  fprintf(analyticout,"Warning, density floor is higher than .001 times of density max!\n");
	}
	
	fprintf(analyticout,"total tori mass: %15.10g total atm mass: %15.10g\n",totalmass_full[1],totalmass_full[0]);
	if(calltype!=2){
	  fclose(analyticout);
	}

	if(calltype==2){
	  return(analyticout);
	}
	else{
	  return(NULL);
	}

      }// end if cpu write
    }//end if firstsolve

  }
    
  //////////////////// VISUALIZATION settings

  // use image() to tell you how to set this(TOTALMM==1, and make sure DYNAMICMM==0)
  // for dvr=.1 Rint=.85

  /*
  // for normal test 
  for(i=0;i<ITYPES;i++){ // both views
  for(j=0;j<CTYPES;j++){ // both computations
  mms[i][j][1][0]=1.e-4;
  mms[i][j][1][1]=1.1;
  mms[i][j][2][0]=3e-08;
  mms[i][j][2][1]=14.0;
  mms[i][j][3][0]=-60.0;
  mms[i][j][3][1]=0.0;

  mmv[i][j][1][0][0]=0;
  mmv[i][j][1][0][1]=2.6;
  mmv[i][j][1][1][0]=-2.5;
  mmv[i][j][1][1][1]=.54;
  mmv[i][j][1][2][0]=-1.5;
  mmv[i][j][1][2][1]=1.5;
  mmv[i][j][1][3][0]=-7.e-7;
  mmv[i][j][1][3][1]=2.3;
    
  mmv[i][j][2][0][0]=0;
  mmv[i][j][2][0][1]=0;
  mmv[i][j][2][1][0]=0;
  mmv[i][j][2][1][1]=0;
  mmv[i][j][2][2][0]=0;
  mmv[i][j][2][2][1]=0;
  mmv[i][j][2][3][0]=0;
  mmv[i][j][2][3][1]=0;
  }
  }
  */
  // for art visc with inner radius at .01 and visc=.05
  /*
    for(i=0;i<ITYPES;i++){ // both views
    for(j=0;j<CTYPES;j++){ // both comps
    mms[i][j][1][0]=1.e-4;
    mms[i][j][1][1]=1.1;
    mms[i][j][2][0]=3e-08;
    mms[i][j][2][1]=14.0;
    mms[i][j][3][0]=-60.0;
    mms[i][j][3][1]=0.0;

    mmv[i][j][1][0][0]=0;
    mmv[i][j][1][0][1]=9.0;
    mmv[i][j][1][1][0]=-9.0;
    mmv[i][j][1][1][1]=4.0;
    mmv[i][j][1][2][0]=-.6;
    mmv[i][j][1][2][1]=1.0;
    mmv[i][j][1][3][0]=-.1;
    mmv[i][j][1][3][1]=9.0;
    
    mmv[i][j][2][0][0]=0;
    mmv[i][j][2][0][1]=0;
    mmv[i][j][2][1][0]=0;
    mmv[i][j][2][1][1]=0;
    mmv[i][j][2][2][0]=0;
    mmv[i][j][2][2][1]=0;
    mmv[i][j][2][3][0]=0;
    mmv[i][j][2][3][1]=0;
    }
    }
  */
  //////////////////////////////

  /*
  //for(i=0;i<ITYPES;i++){ // both views
  //  for(j=0;j<CTYPES;j++){ // both comps

  // scalars

  j=0; // normal output
  i=0; // view large

  mms[i][j][1][0]=          (9./10.)*1.e-06;  
  mms[i][j][1][1]=     1.02012908;

  mms[i][j][2][0]=(9./10.)*1.33333333e-11;  
  mms[i][j][2][1]=   0.7428376443;

  mms[i][j][3][0]=   -3.265618839;
  mms[i][j][3][1]=  -0.1881626381;

  j=1;  // second type of comp
  i=0; // view large
  
  mms[i][j][1][0]=          (9./10.)*1.e-06;  
  mms[i][j][1][1]=     1.02012908;

  mms[i][j][2][0]=(9./10.)*1.33333333e-11;  
  mms[i][j][2][1]=   0.7428376443;

  mms[i][j][3][0]=   -3.265618839;
  mms[i][j][3][1]=  -0.1881626381;


  j=0; // normal comp
  i=1; // view zoom
  // for alpha=.01 with both terms
  mms[i][j][1][0]=          .0001;  
  mms[i][j][1][1]=     .2;

  mms[i][j][2][0]=.0001;  
  mms[i][j][2][1]=   0.7428376443;

  mms[i][j][3][0]=   -3.265618839;
  mms[i][j][3][1]=  -0.1881626381;

  j=1; // 2nd comp
  i=1; // view zoom
  // for alpha=.01 with both terms
  mms[i][j][1][0]=          .0001;  
  mms[i][j][1][1]=     .1;

  mms[i][j][2][0]=.0001;  
  mms[i][j][2][1]=   0.7428376443;

  mms[i][j][3][0]=   -3.265618839;
  mms[i][j][3][1]=  -0.1881626381;


  // vectors
  
  i=0; // normal view
  j=0; // 1st comp
  // available from ipar.out

  // v0
  mmv[i][j][1][0][0]=              0.0;
  mmv[i][j][1][0][1]=    10.0;
  // vx1
  mmv[i][j][1][1][0]=   -2.0;
  mmv[i][j][1][1][1]=    .1;
  // vx2
  mmv[i][j][1][2][0]=   -1.5;
  mmv[i][j][1][2][1]=    1.5;
  // vx3
  mmv[i][j][1][3][0]=-.1;
  mmv[i][j][1][3][1]=    10.0;
  
  // B0
  mmv[i][j][2][0][0]=0;
  mmv[i][j][2][0][1]=0;
  // Bx1
  mmv[i][j][2][1][0]=0;
  mmv[i][j][2][1][1]=0;
  // Bx2
  mmv[i][j][2][2][0]=0;
  mmv[i][j][2][2][1]=0;
  // Bx3
  mmv[i][j][2][3][0]=0;
  mmv[i][j][2][3][1]=0;

  i=0;
  j=1;
  // rho times v or B
  // must get from sm plots

  // rho*v0
  mmv[i][j][1][0][0]=              0.0;
  mmv[i][j][1][0][1]=    1.0;
  // rho*vx1
  mmv[i][j][1][1][0]=   -0.01;
  mmv[i][j][1][1][1]=    0.01;
  // rho*vx2
  mmv[i][j][1][2][0]=   -0.01;
  mmv[i][j][1][2][1]=    .01;
  // rho*vx3
  mmv[i][j][1][3][0]=-.1;
  mmv[i][j][1][3][1]=    10.0;

  // rho*B0
  mmv[i][j][2][0][0]=0;
  mmv[i][j][2][0][1]=0;
  // rho*Bx1
  mmv[i][j][2][1][0]=0;
  mmv[i][j][2][1][1]=0;
  // rho*Bx2
  mmv[i][j][2][2][0]=0;
  mmv[i][j][2][2][1]=0;
  // rho*Bx3
  mmv[i][j][2][3][0]=0;
  mmv[i][j][2][3][1]=0;



  //////////////// ZOOM  
  i=1; // zoom view (vectors)
  j=0; // 1st comp
  // available from ipar.out

  // v0
  mmv[i][j][1][0][0]=              0.0;
  mmv[i][j][1][0][1]=    10.0;
  // vx1
  mmv[i][j][1][1][0]=   -3.0;
  mmv[i][j][1][1][1]=    .1;
  // vx2
  mmv[i][j][1][2][0]=   -1.5;
  mmv[i][j][1][2][1]=    1.5;
  // vx3
  mmv[i][j][1][3][0]=-1;
  mmv[i][j][1][3][1]=    10.0;
  
  // B0
  mmv[i][j][2][0][0]=0;
  mmv[i][j][2][0][1]=0;
  // Bx1
  mmv[i][j][2][1][0]=0;
  mmv[i][j][2][1][1]=0;
  // Bx2
  mmv[i][j][2][2][0]=0;
  mmv[i][j][2][2][1]=0;
  // Bx3
  mmv[i][j][2][3][0]=0;
  mmv[i][j][2][3][1]=0;

  i=1; // zoom view (vectors)
  j=1;
  // rho times v or B
  // must get from sm plots

  // rho*v0
  mmv[i][j][1][0][0]=              0.0;
  mmv[i][j][1][0][1]=    1.0;
  // rho*vx1
  mmv[i][j][1][1][0]=   -0.01;
  mmv[i][j][1][1][1]=    0.01;
  // rho*vx2
  mmv[i][j][1][2][0]=   -0.01;
  mmv[i][j][1][2][1]=    .01;
  // rho*vx3
  mmv[i][j][1][3][0]=-.1;
  mmv[i][j][1][3][1]=    1.0;

  // rho*B0
  mmv[i][j][2][0][0]=0;
  mmv[i][j][2][0][1]=0;
  // rho*Bx1
  mmv[i][j][2][1][0]=0;
  mmv[i][j][2][1][1]=0;
  // rho*Bx2
  mmv[i][j][2][2][0]=0;
  mmv[i][j][2][2][1]=0;
  // rho*Bx3
  mmv[i][j][2][3][0]=0;
  mmv[i][j][2][3][1]=0;

  */
  // mhdtorin1
  for(outtype=0;outtype<=1;outtype++){
    mms[outtype][0][1][0]=           9.9e-06 ;
    mms[outtype][1][1][0]=           9.9e-06 ;
    mms[outtype][0][1][1]=    1.0 ;
    mms[outtype][1][1][1]=    1.0 ;
    mms[outtype][0][2][0]=  1E-11 ;
    mms[outtype][1][2][0]=  1E-11 ;
    mms[outtype][0][2][1]=  0.004 ;
    mms[outtype][1][2][1]=  0.004 ;
    mms[outtype][0][3][0]=   -2.56;
    mms[outtype][1][3][0]=   -2.56;
    mms[outtype][0][3][1]=  -0.04818533834 ;
    mms[outtype][1][3][1]=  -0.04818533834 ;  
    mmv[outtype][0][1][1][0]=    -2.5 ;
    mmv[outtype][1][1][1][0]=    -2.5 ;
    mmv[outtype][0][1][1][1]=    1.0 ;
    mmv[outtype][1][1][1][1]=    1.0 ;
    mmv[outtype][0][1][2][0]=   -1.2 ;
    mmv[outtype][1][1][2][0]=   -1.2 ;
    mmv[outtype][0][1][2][1]=   1.2 ;
    mmv[outtype][1][1][2][1]=   1.2 ;
    mmv[outtype][0][1][3][0]=   -1.5 ;
    mmv[outtype][1][1][3][0]=   -1.5 ;
    mmv[outtype][0][1][3][1]=     1.5 ;
    mmv[outtype][1][1][3][1]=     1.5 ;
    mmv[outtype][0][2][1][0]=  -0.05 ;
    mmv[outtype][1][2][1][0]=  -0.05 ;
    mmv[outtype][0][2][1][1]=   0.05 ;
    mmv[outtype][1][2][1][1]=   0.05 ;
    mmv[outtype][0][2][2][0]= -0.02 ;
    mmv[outtype][1][2][2][0]= -0.02 ;
    mmv[outtype][0][2][2][1]= 0.02 ;
    mmv[outtype][1][2][2][1]= 0.02 ;
    mmv[outtype][0][2][3][0]=  -0.06 ;
    mmv[outtype][1][2][3][0]=  -0.06 ;
    mmv[outtype][0][2][3][1]=   0.06 ;
    mmv[outtype][1][2][3][1]=   0.06 ;
  }
  // define outer region when interpolation is used.
  // same order as scalar/vector arrays

  // for images
  for(i=0;i<ITYPES;i++){ // both views
    for(j=0;j<CTYPES;j++){ // both comps
      outerdefs[i][j][1]=mms[i][j][1][0]; // rho
      outerdefs[i][j][2]=mms[i][j][2][0]; // en
      outerdefs[i][j][3]=mms[i][j][3][0]; // pot
      
      outerdefv[i][j][1][0]=0.0; // magnitude of v
      outerdefv[i][j][1][1]=0.0; // v1
      outerdefv[i][j][1][2]=0.0; // v2
      outerdefv[i][j][1][3]=0.0; // v3
      
      outerdefv[i][j][2][0]=0.0;
      outerdefv[i][j][2][1]=0.0;
      outerdefv[i][j][2][2]=0.0;
      outerdefv[i][j][2][3]=0.0;
      
    }
  }
  // for dumps(never use 2nd comp type, and always same for any view)
  douterdefs[1]=mms[0][0][1][0]; // rho
  douterdefs[2]=mms[0][0][2][0]; // en
  douterdefs[3]=mms[0][0][3][0]; // pot

  douterdefv[1][0]=0.0; // magnitude of v
  douterdefv[1][1]=0.0; // v1
  douterdefv[1][2]=0.0; // v2
  douterdefv[1][3]=0.0; // v3

  douterdefv[2][0]=0.0;
  douterdefv[2][1]=0.0;
  douterdefv[2][2]=0.0;
  douterdefv[2][3]=0.0;

  /*
    LOOPF{
    if(sanal[1][k][j][i]<DENSITYFLOOR){
    fprintf(stderr,"problem\n");
    fflush(stderr);
    exit(1);
    }
    }
  */
  return(0);
}


#define HEADER4_SRE  "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "

void  reentrancezz(void)
{
#if(REENTER)
   FILE *fpre;
   char buffer[1000],*cptr,namere[1000];
   int l,i,j,k,m,tempreint;
   FTYPE tempre;

   xre=(FTYPE (*) [3][NREBIG])(&(xrea[-1][-1][NBIGBND]));
   datre=(FTYPE (*) [N3M][N2RE][N1RE])(&(datrea[-1][N3BND][N2BND][N1BND]));

   kre=1-1;
   jre=N2RE-1;
   ire=N1RE-1;
// read grid
 for(l=1;l<=2;l++){
   if(l==1){
    sprintf(namere,"./allviscburst/0_gridact1.par");
   }else{
    sprintf(namere,"./allviscburst/0_gridact2.par");
   }

   if((fpre=fopen(namere,"rt"))==NULL){
        fprintf(fail_file,"grid1: Cannot open: %s\n",namere);
        myexit(1);
   }

   for(i=1;i<=3;i++){
     cptr=fgets(buffer,1000,fpre);
     if(cptr!=NULL){
//       fprintf(stderr,"read grid %s \n",cptr);
     }
   }
   for(k=0;k<=kre;k++){
   for(j=0;j<=jre;j++){
   for(i=0;i<=ire;i++){
       fscanf(fpre,HEADER4_SRE,&tempre,&tempre,&tempre,&tempre,&tempre,&tempre,&tempre,&xre[l][1][i],&xre[l][2][j],&xre[l][3][k],
              &tempre,&tempre,&tempre,&tempre,&tempre,&tempre,&tempre,&tempre,&tempre,&tempre,&tempre);
//       fprintf(log_file,"xrej %lf %lf %lf %d %d %d \n",xre[l][1][i],xre[l][2][j],xre[l][2][0],i,j,k);
       for(m=1;m<=15;m++){
         fscanf(fpre,"%d",&tempreint);
       }
   }
   }
   }
   fclose(fpre);
//   fprintf(stderr,"xre %lf %d \n",xre[l][1][0],l);
 }
// read data
   sprintf(namere,"./allviscburst/dump0001.dat");
   if((fpre=fopen(namere,"rt"))==NULL){
        fprintf(fail_file,"grid1: Cannot open: %s\n",namere);
        myexit(1);
   }

   for(i=1;i<=7;i++){
     cptr=fgets(buffer,1000,fpre);
     if(cptr!=NULL){
//       fprintf(stderr,"read grid %s",cptr);
     }
   }
   for(k=0;k<=kre;k++){
   for(j=0;j<=jre;j++){
   for(i=0;i<=ire;i++){
     for(l=1;l<=39;l++){
       fscanf(fpre,"%lf",&datre[l][k][j][i]);
     }
   }
   }
   }
   fprintf(log_file," %21.7e xrei xrej %lf %lf",datre[1][0][1][2],xre[1][1][0],xre[1][2][0]);
   fclose(fpre);
#endif
}

FTYPE extre2(FTYPE radius, FTYPE theta, int m)
{
#if(REENTER)
   int l,dir,nre[2],i,j,k;
   FTYPE xretemp[2],extretem1,extretem2,extre;
   int o1[2],o2[2];

   xretemp[0]=radius;
   xretemp[1]=theta;
   nre[0]=ire;
   nre[1]=jre;
   if(m<=3){
     l=2;
   }else{
     l=1;
   }
 for(dir=1;dir<=2;dir++){
   if(xre[l][dir][2]>xre[l][dir][1]){
     o1[dir-1]=0;
     o2[dir-1]=1;
     for(i=0;i<=nre[dir-1]-1;i++){
       if(xretemp[dir-1]>=xre[l][dir][i]){
         o1[dir-1]=i;
         o2[dir-1]=i+1;
       }
     }
//spetial treatment for beyond the boundary
     if(xretemp[dir-1]<=xre[l][dir][0]){xretemp[dir-1]=xre[l][dir][0];}
     if(xretemp[dir-1]>=xre[l][dir][nre[dir-1]]){xretemp[dir-1]=xre[l][dir][nre[dir-1]];}
//end
   }else{
     o1[dir-1]=0;
     o2[dir-1]=1;
     for(i=0;i<=nre[dir-1]-1;i++){
       if(xretemp[dir-1]<=xre[l][dir][i]){
         o1[dir-1]=i;
         o2[dir-1]=i+1;
       }
     }
//spetial treatment for beyond the boundary
     if(xretemp[dir-1]>=xre[l][dir][0]){xretemp[dir-1]=xre[l][dir][0];}
     if(xretemp[dir-1]<=xre[l][dir][nre[dir-1]]){xretemp[dir-1]=xre[l][dir][nre[dir-1]];}
//end
   }
 }
//
//   if(m==1){
//   fprintf(log_file,"%10.3lf %10.3lf %10.3lf %10.3lf %8d %8d\n",xretemp[0],theta,xre[l][1][o1[0]],xre[l][2][o1[1]],o1[0],o1[1]);
//   }
   extretem1=(datre[m][0][o1[1]][o1[0]]-datre[m][0][o1[1]][o2[0]])/(xre[l][1][o1[0]]-xre[l][1][o2[0]])
                *(xretemp[0]-xre[l][1][o1[0]])+datre[m][0][o1[1]][o1[0]];
   extretem2=(datre[m][0][o2[1]][o1[0]]-datre[m][0][o2[1]][o2[0]])/(xre[l][1][o1[0]]-xre[l][1][o2[0]])
                *(xretemp[0]-xre[l][1][o1[0]])+datre[m][0][o2[1]][o1[0]];
   extre=(extretem1-extretem2)/(xre[l][2][o1[1]]-xre[l][2][o2[1]])*(xretemp[1]-xre[l][2][o1[1]])+extretem1;
//   fprintf(stderr,"%lf %lf\n",extre,datre[m][0][o1[1]][o1[0]]);
   return(extre);
#endif
}


FTYPE extre(FTYPE radius, FTYPE theta, int m)
{
#if(REENTER)
   int l[3],dir,nre[2],i,j,k;
   FTYPE xretemp[2],extretem1,extretem2,extre;
   int o1[2],o2[2];

   xretemp[0]=radius;
   xretemp[1]=theta;
   nre[0]=ire;
   nre[1]=jre;
   if(m<=3){
     l[1]=2;
     l[2]=2;
   }
   if(m==4){
     l[1]=1;
     l[2]=2;
   }
   if(m==5){
     l[1]=2;
     l[2]=1;
   }
   if(m==6){
     l[1]=2;
     l[2]=2;
   }

 for(dir=1;dir<=2;dir++){
   if(xre[l[dir]][dir][2]>xre[l[dir]][dir][1]){
     o1[dir-1]=0;
     o2[dir-1]=1;
     for(i=0;i<=nre[dir-1]-1;i++){
       if(xretemp[dir-1]>=xre[l[dir]][dir][i]){
         o1[dir-1]=i;
         o2[dir-1]=i+1;
       }
     }
//spetial treatment for beyond the boundary
     if(xretemp[dir-1]<=xre[l[dir]][dir][0]){xretemp[dir-1]=xre[l[dir]][dir][0];}
     if(xretemp[dir-1]>=xre[l[dir]][dir][nre[dir-1]]){xretemp[dir-1]=xre[l[dir]][dir][nre[dir-1]];}
//end
   }else{
     o1[dir-1]=0;
     o2[dir-1]=1;
     for(i=0;i<=nre[dir-1]-1;i++){
       if(xretemp[dir-1]<=xre[l[dir]][dir][i]){
         o1[dir-1]=i;
         o2[dir-1]=i+1;
       }
     }
//spetial treatment for beyond the boundary
     if(xretemp[dir-1]>=xre[l[dir]][dir][0]){xretemp[dir-1]=xre[l[dir]][dir][0];}
     if(xretemp[dir-1]<=xre[l[dir]][dir][nre[dir-1]]){xretemp[dir-1]=xre[l[dir]][dir][nre[dir-1]];}
//end
   }
 }
//
   if(m==1){
   fprintf(log_file,"%10.3lf %10.3lf %10.3lf %10.3lf %8d %8d\n",xretemp[0],theta,xre[2][1][o1[0]],xre[2][2][o1[1]],o1[0],o2[0]);
   }
   extretem1=(datre[m][0][o1[1]][o1[0]]-datre[m][0][o1[1]][o2[0]])/(xre[l[1]][1][o1[0]]-xre[l[1]][1][o2[0]])
                *(xretemp[0]-xre[l[1]][1][o1[0]])+datre[m][0][o1[1]][o1[0]];
   extretem2=(datre[m][0][o2[1]][o1[0]]-datre[m][0][o2[1]][o2[0]])/(xre[l[1]][1][o1[0]]-xre[l[1]][1][o2[0]])
                *(xretemp[0]-xre[l[1]][1][o1[0]])+datre[m][0][o2[1]][o1[0]];
   extre=(extretem1-extretem2)/(xre[l[2]][2][o1[1]]-xre[l[2]][2][o2[1]])*(xretemp[1]-xre[l[2]][2][o1[1]])+extretem1;
//   fprintf(log_file,"%10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %d\n",extre,xre[2][1][o1[0]],xre[l[1]][1][o2[0]],extretem1,extretem2,datre[1][0][o1[1]][o1[0]],m);
   return(extre);
#endif
}

void tempini(int initmode)
{
  int i,j,k;
  LOOPC3{
     LOOPC1{
       if(initmode==0){
        optdini[k][-1][i]=0.;
        sigini[k][-1][i]=0.;
       }
       LOOPC2{
         optdini[k][j][i]=optdini[k][j-1][i]+dx[2][2][j]*G2(2,i)*opaini[k][j][i];
         sigini[k][j][i]=sigini[k][j-1][i]+dx[2][2][j]*G2(2,i)*sanal[1][k][j][i];
       }
     }
  } 
    LOOPC3{
    LOOPC1{
    LOOPC2{
      if(sigini[k][j][i]<=Siga){
        tirr0[k][j][i]=pow(3./4.*teff4[i]*(optdini[k][j][i]+2./3.),0.25);
      }else{
        tirr0[k][j][i]=tirr0[k][j-1][i];
      }
    }   
    }
    }
  fprintf(log_file,"sigma1ini %15.10g %15.10g %15.10g %15.10g \n",sigini[0][-N2BND][10],sigini[0][-N2BND+1][10],sigini[0][0][10],sigini[0][1][10]);
  fprintf(log_file,"sigma2ini %15.10g %15.10g %15.10g %15.10g \n",sigini[0][N2-N2BND][10],sigini[0][N2-N2BND+1][10],sigini[0][N2][10],sigini[0][N2+1][10]);

}



// General analytic dump for given analytic solution
void analsolve_normal(int calltype)
{
  // if calltype>=1 or <0, then always do the call.  ==0, only do if timeless==0
  static int timeless=0; // =0 not computed stationary solution yet 1= have, so skip

  if(calltype==100){
    fprintf(log_file,"begin: analsolve_normal(100) ... "); fflush(log_file);
  }

  if(analoutput==1)      sodsol(calltype);
  else if(analoutput==2) advsol(calltype);
  else if(analoutput==3) gausssol(calltype);
  else if(analoutput==4){
    if((timeless==0)||(calltype)){
      if(calltype!=100) timeless=1;
      bondisol(calltype,NULL);
    }
  }
  else if(analoutput==5){
    if((timeless==0)||(calltype)){
      if(calltype!=100) timeless=1;
      tori1sol(calltype);
    }
  }
  else if(analoutput==6){
    visctermtest(calltype);
  }
  else if(analoutput==7){
    pulsesol(calltype);
  }
  else if(analoutput==8){
    if((timeless==0)||(calltype)){
      if(calltype!=100) timeless=1;
#if(MDOTMEM)
      if(myid<=0){
	fprintf(logfull_file,"injectsol\n");
	fflush(logfull_file);
      }
      injectsol(calltype);
#endif
    }
  }
  else if(analoutput==9){
    magbreaksol(calltype);
  }
  else if(analoutput==10){
    if((timeless==0)||(calltype)){
      if(calltype!=100) timeless=1;
      magvortex(calltype);
    }
  }
  else if(analoutput==11){
    magcorona(calltype);
  }
  else if(analoutput==12){
    if((timeless==0)||(calltype)){
      if(calltype!=100) timeless=1;
      waves(calltype);
    }
  }
  else if(analoutput==13){
    if((timeless==0)||(calltype)){
      if(calltype!=100) timeless=1;
      chandran(calltype);
    }
  }
  else if(analoutput==14){
    if((timeless==0)||(calltype)){
      if(calltype!=100) timeless=1;
      bondiplustorus(calltype);
    }
  }
  else if(analoutput==15){
    if((timeless==0)||(calltype)){
      if(calltype!=100) timeless=1;
      averystar(calltype);
    }
  }
  else if(analoutput==16){
    if((timeless==0)||(calltype)){
      if(calltype!=100) timeless=1;
      averywind(calltype);
    }
  }
  else if(analoutput==17){
    if((timeless==0)||(calltype)){
      if(calltype!=100) timeless=1;
      diskini(calltype);
    }
  }

      

  if(calltype==100){
    fprintf(log_file," end: analsolve_normal(100)\n"); fflush(log_file);
  }


}











