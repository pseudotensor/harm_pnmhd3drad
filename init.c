#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif

#include "radiation.h"

//MARK need to fix since not all are sensitive!
#if(SENSITIVE==1)
#define INPUT1 "%lf %d %d"
#define INPUT1OLD "%lf %d"
#define INPUTAVGH1 "%lf %lf"
#define HEADER3_S2 "%lf %lf %lf %lf %lf %lf\n"
#define HEADER3_S3 "%lf %lf %lf %lf %lf %lf %lf\n"
#define HEADER3_S4 "%lf %lf %lf %lf %lf\n"
#define HEADER3_S5 "%lf %lf %lf %lf\n"
#define HEADER3_S6 "%lf %lf %lf %lf %d %lf\n"
#define HEADER3_S7 "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n"
#define HEADER3_S8 "%lf %lf\n"
#define HEADER3_S9 "%lf %lf %lf\n"
#else
#define INPUT1 "%f %d %d"
#define INPUT1OLD "%f %d"
#define INPUTAVGH1 "%f %f"
#define HEADER3_S2 "%f %f %f %f %f %f\n"
#define HEADER3_S3 "%f %f %f %f %f %f %f\n"
#define HEADER3_S4 "%f %f %f %f %f\n"
#define HEADER3_S5 "%f %f %f %f\n"
#define HEADER3_S6 "%f %f %f %f %d %f\n"
#define HEADER3_S7 "%f %f %f %f %f %f %f %f %f %f\n"
#define HEADER3_S8 "%f %f\n"
#define HEADER3_S9 "%f %f %f\n"
#endif

#define PARAMETERPARLIST int seed, FTYPE beta, FTYPE nj

// 1/2/3-d valid
int init_paramstot(PARAMETERPARLIST) // anything here can be overridden in analsol.c if desired, and the code will output the parameter you choose in the parm file.
{
  SFTYPE ftemp,tvisc;
  SFTYPE numvtime;
  /* BEGIN: Assign global parameters */
  fprintf(log_file,"begin: init_paramstot ... "); fflush(log_file);

  // stuff only needed if not reading in parameter data
  if((runtype==1)||(runtype==0)){ // stuff here stored in par file

    init_paramsphysics(seed,beta,nj);
    init_paramsnumerical(seed,beta,nj);
    init_paramstime(seed,beta,nj);
    init_paramsgeom(seed,beta,nj);
    analsolve(100); // 100 signifies to analytic solution to just setup problem as required, adjusting parameters that are just above set to defaults
    // just comment out above to force to use parameters set here in init.c instead of ever using analsol.c sets
  }
  fprintf(log_file," end: init_paramstot\n"); fflush(log_file);
  return(0);
}


int init_paramsphysics(PARAMETERPARLIST) // physics stuff
{
  fprintf(log_file,"begin: init_paramsphysics ... "); fflush(log_file);
  //////////////////////////////
  //
  // physics

  // MARK
  trans=      1; // 0->no trans 1->do trans step at all

  transx1=    1; // 0->no trans 1->do trans step in x1-dir
  transrhox1= 1; // 0->do not transport mass density in x1-dir 1->do
  transiex1=  1; // 0->do not transport internal energy density in x1-dir 1->do
  transv1x1=  1; // 0->do not transport velocity1 in x1-dir 1->do
  transv2x1=  1; // 0->do not transport velocity2 in x1-dir 1->do
  transv3x1=  1; // 0->do not transport velocity3 in x1-dir 1->do

  transx2=    1; // 0->no trans 1->do trans step in x2-dir
  transrhox2= 1; // 0->do not transport mass density in x2-dir 1->do
  transiex2=  1; // 0->do not transport internal energy density in x2-dir 1->do
  transv1x2=  1; // 0->do not transport velocity1 in x2-dir 1->do
  transv2x2=  1; // 0->do not transport velocity2 in x2-dir 1->do
  transv3x2=  1; // 0->do not transport velocity3 in x2-dir 1->do

  if(COMPDIM<3){
    // modified in stepvar_2d for other algorithms
    transmagx1= 1; // 0->do not transport magnetic 3-flux in x1-dir 1->do
    transmagx2= 1; // 0->do not transport magnetic 3-flux in x2-dir 1->do
  }
  else{
    transmagx1= 0;
    transmagx2= 0;
  }

  transx3=    1; // 0->no trans 1->do trans step in x2-dir
  transrhox3= 1; // 0->do not transport mass density in x2-dir 1->do
  transiex3=  1; // 0->do not transport internal energy density in x2-dir 1->do
  transv1x3=  1; // 0->do not transport velocity1 in x2-dir 1->do
  transv2x3=  1; // 0->do not transport velocity2 in x2-dir 1->do
  transv3x3=  1; // 0->do not transport velocity3 in x2-dir 1->do
  

  // use Bx as passive scalar to advect (see sweep.c)
  // has to be turned on by user in analsol.c
  // assume not restarting with this on (so not present in init.c for restart file)
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


  press=      1; // 0->no pressure at all 1->pressure at all (includes curvature! -- doesn't go away when N?==1) (check pressure terms in global.h)
  pressx1=1;
  pressx2=1;
  pressx3=1;
  
  mag= 1;        // 0->no mag steps 1->do mag steps(check pressure in global.h)
  if(COMPDIM<3){
    transbzx=1;    // really: global.h:TORV/B X1/X2 for v3/b3-comp
    transbzy=1;
  }
  else{
    transbzx=0;
    transbzy=0;
  }
  stepmocct=1;
  mocctvx1=1;
  mocctvx2=1;
  mocctvx3=1; // only applies if COMPDIM==3 
  mocctbx1=1;
  mocctbx2=1;
  mocctbx3=1; // only applies if COMPDIM==3 

  res_real= 1; // 0->no resistivity step 1->do res step
  rreal=2;
  // 1: gammie
  // 2: SP00
  // 3: constant
  resheat= 1; // only relevant if res_real>0

  // DEBUG REMARKS
  //
  // can turn off velocity source terms and just transport everything
  //
  // press=0; mocctvx1=0; mocctvx2=0; // and TORX1=TORX2=0 in global.h
  // and turn off transv1x1/x2/x3 v2x1/x2/x3
  // then turn trans of v back on to see if that's the problem
  //
  // can turn off each component of v
  //
  // turn on/off ie and transiex1/x2
  //
  // with multiple terms of above, can turn off each one to see which dominates error, if any
  //
  // can check symmetry by turning off any code component!

  // MARK
  ie=         1; // 0-> no ie step 1->do ie step (controls art viscosity steps too)
  visc_art=    1; // 0->no artificial visc step 1->do visc step
  // MARK
  visc_real=   0; // 0-> do not do real visc 1-> do
  vreal=       3;
  // 1: gammie1(visce13 on)
  // 2: igumenshchev1(all visce terms on)
  // 3: stone et al. 1(visce13,visce23 on)
  // 4: mg (igu with sin)
  // 5: constant
  vischeat= 1; // only relevant if visc_real>0
  mdotin=0; // 0->no mass influx on grid, 1: yes (see injection()) (make sure MDOTMEM is set right in global.h)

  cool=0; // 0->no cooling 1-> cooling

  // algorithm/implementation details
  advint=1;
  // 0=DONOR
  // 1=VANLEER
  // 2 = WOODWARD

  kever=1;
  // 0: use native fluxes when possible for ke (not good for periodic)
  // 1: use mass flux as basis for flux of ke (best)
  


  // MARK
  gam    = (5.0/3.0) ; // adiabatic index
  alpha_real0 = 0.01 ; // alpha(see step.c:tdep_compute())
  n_real = 2.0 ;       // power of sin(theta) for Length scale for Gammie visc
  coolfact=0; // factor by which cooling function is multiplied(see step.c:cooling())
  alpha  = 5.0/2.0; // equipartition parameter(for gam=1 case)
  cs     = (1.0) ; // sound speed, only used for gam=1.0
  RMAX=0.0; // tells whether doing this or not
  // SP00 value for RunA
  //  resist_real0 = 0.1 ; // dimensionless resistivity constant (pretty high for Orzag-Tang vortex--kills current sheet within 10s of zones and no bubbles of field)  but does kill point shocks even with lorentz before emfs
  resist_real0 = .01; // preserves most structure with minimal point shocks in vortex problem, so use this
  //MARK for boundary problem testing(whether resistivity causes problem at boundary)
  //resist_real0 = 0.1*100.0 ; // dimensionless resistivity constant
  blackholejz=1.0; // number of J/M's, where maximum is 1
  // inverse speed of light squared for alfven trick
  invsol2=invsol2_original=1.0;

  fprintf(log_file,"end: init_paramsphysics\n"); fflush(log_file);
  return(0);
}

int init_paramsgeom(PARAMETERPARLIST) // geometry stuff
{
  fprintf(log_file,"begin: init_paramsgeom ... "); fflush(log_file);
  // zone to start/stop on, in which to integrate over for integrated values, and flux terms
  // used to avoid problems with boundary conditions
  // NOTE: Outer boundary is loop ending value where <INTOX, not <=INTOX !
  // restrict boundary to avoid error due to outflow boundary condition, factor of 20 for injection case in error rate decrease
  // global, need to correct per cpu
  // only need 2 here if coord==3
#if((COORD==3)&&(N1>3))
  intix1= 2; // must be less than size of first cpu in radius
#else
  intix1=0;
#endif
  intox1= N1; // should keep close to outer edge since otherwise too different for different res, and injection won't add up
  intix2= 0;
  intox2= N2;
  intix3= 0;
  intox3= N3;


  // MARK
  nonunigridx1= 5; // 0=unigrid 1+=nonunigrid
  nonunigridx2= 3; // 0=unigrid 1+=nonunigrid
  nonunigridx3= 0; // 0=unigrid 1+=nonunigrid
//x1:  1 : Cos in r(puts highest res in middle)
//     2 : Log in r (puts highest res in inner r-edge)
//     3 : Power law per decade (unconstrained Rin)
//     4 : Power law per decade (dx/x=constant - constrains Rin)
//     5 : Power law per decade equal, spacing in log(r-rgp)
//x2:  1 : Cos in theta(puts highest res in middle of theta=pi/2
//     2 : Log split at middle
//     3 : Power law per decade

  // must be correct for whatever parameter setup or runtype
  simplebc=1; // choice: whether to use simple boundary conditions or not
  // those simple preferences
  // MARK
  bcix1=4;
  bcox1=4;
  bcix2=1; // 5=special periodic (0+dphi -> Pi+dphi)
  bcox2=1;
  bcix3=5;
  bcox3=5;

  // remember to set REFLECTIX1=1 if you set below to 0 for spc.
  // and set SKIPIX1 to 1
  rg    = 2.0;  // rg=2GM/c^2
  rgp   = 0.;  // PW-potential gravitational radius(set by units==2 if chosen) (rg==0 if Newtonian pot)
  //rgp=0.0;
  Rinner=1.2*rg;
  Router=11.5*rg;

  x1in=Rinner;
  x1out= Router;
  x2in=  (0.0);
  x2out= (M_PI);
  x3in=  (0.0);
  x3out= (2.0*M_PI);

  fprintf(log_file,"end: init_paramsgeom\n"); fflush(log_file);
  return(0);
}

int init_paramstime(PARAMETERPARLIST) // time stuff
{
  fprintf(log_file,"begin: init_paramstime ... "); fflush(log_file);
  tstart = 0.0;  //start unless reading in start time

  tavgi  = tf*.5; // time to start averaging
  tavgf  = tf*.8; // time to stop averaging(only approximate)
  numavg = 139;

  timagescale=tf*0.5;
  //#define DTDUMP (1.0)
  //#define DTCONS (0.1)

  //  DTl = 1.0;
  DTl = 1.0/30.0;
  // strictly for purposes of reentrant ability, should have dump/floor same DT
  DTfld=DTfloor=DTd    = DTl*300.0 ;      /* dumping period(data and analytic)*/
  // MARK
  //DTfloor=DTd    = 100.0 ;      /* dumping period(data and analytic)*/
  //DTfloor=DTd = DTDUMP;
  // strictly for reentract ability, DTi should be integral number of DTds.
  DTpd = 10.0*DTd; // always reentrant data
  //DTpd = DTDUMP;
  DTi    = DTl*30.0 ;        /* image period */
  //DTi=DTDUMP;
  //DTi=30.0;

  //below not restricted by DTl calls, just each self
  DTener = DTl*1.0 ; // negative means do every time step (must be multiple of DTl)
  // f2's FFT shows could do 20.0 here(was 4.0)
  DTloss = DTl*1.0 ; // negative means do every time step (must be multiple of DTl
  DTmode = DTl*100.0 ; // negative means do every time step (must be multiple of DTl)
  //DTloss=DTCONS;
  DTtimestep = DTl*100.0 ; // how often to output dominate term in timestep to failfile
  DTsp = DTl*100.0 ; // how often to output sonic point info
  //DTtimestep=DTCONS;
  DTtimescale = DTl*100.0;
  DTdivb = DTl*100.0;

  // MARK for mag test
  //  DTl=DTener=DTloss=DTtimestep=DTsp=tf/100.0;  DTfloor=DTd=tf/(2.5806452*10.0); DTpd=tf/(2.5806452*5.0); DTi=tf/(2.5806452*300.0);
  //DTl=DTener=DTloss=DTtimestep=DTsp=tf/100.0;  DTfloor=DTd=tf/(2.5806452*10.0); DTpd=tf/(2.5806452*5.0); DTi=tf/(2.5806452*30.0);
  // for long cs problem
  //  DTi = tf/(5000.0);
  //  DTl=DTener=DTloss=DTtimestep=DTsp=tf/1000.0;  DTfloor=DTd=tf/(2.5806452*10.0); DTpd=tf/(2.5806452*1.0); DTi=tf/(2.5806452*3.0);
  //DTd=0.2;

  fprintf(log_file,"end: init_paramstime\n"); fflush(log_file);
  return(0);
}  


int init_paramsnumerical(PARAMETERPARLIST) // usually never need to modify this
{

  fprintf(log_file,"begin: init_paramsnumerical ... "); fflush(log_file);
  // nu_l==0.05(min) to .1 to .2(max) works best for oscillations near shocks
  // the bigger you make this, the lower the oscillations, but the longer time it takes and the worse errors will become around shocks due to non-oscilatting parts.  Can let some oscillations appear for sake of those non-oscillating errors.
  // nu_vnr==3.0 works best for shocks
  //nu_vnr = 3.0 ; // Von Neuman Rictmeyer artificial viscosity
  nu_vnr = 2.0 ; // Von Neuman Rictmeyer artificial viscosity
  nu_l   = 0.1 ; // linear viscosity(not needed unless want pretty shocks)(see global.h)
  nu_ten = 3.0 ; // tensor viscosity
  GRAVC  = 1.0; // Gravitational Constant
  MASSBH   = 1.0; // Mass of central object
  GM     = 1.0; // as units, should be set in analsol.c, not here
  dt     = 1.e10 ; // start time step, used as reference, keep large!

  // courant factor
  // <1.0 works for linear waves without shocks
  // <0.5 works well for all shocks
  // <

  // the lower the courant number the less diffusive the evolution

  // the nonsplit moc_ct routine requires a lower courant per dimension

#if(TVDLF==0)
  // MARK
  if(SPLITMETHOD==0){
    cour = 0.4999;
    if(N1OFF+N2OFF+N3OFF==1){
      cour/=1.0;
    }
  }
  else{
    cour = 0.4999;
  }
  // overrides
  //cour   = 0.4999 ; // courant condition
  //  cour=1.0; // for 1D linear waves
  //cour=0.45;
    // GODMARK
  //  cour=0.1;
  //  cour=0.25;
#else
  cour   = 0.25;
#endif
  cour2   = 0.24999 ; // diffuse courant condition

  nu_sh  = 0.0 ; // not used
  vg[1]     = 0.0 ; // velocity of grid in x1-dir, not used
  vg[2]     = 0.0 ; // velocity of grid in x2-dir, not used
  vg[3]     = 0.0 ; // velocity of grid in x3-dir, not used


  if(POSTPROC==1){
    DYNAMICMM=0; // already read in file usually
  }
  else{
    DYNAMICMM=2; //MARK // problem with 3 and vortex?
  }


  fprintf(log_file,"end: init_paramsnumerical\n"); fflush(log_file);
  return(0);
}


// below values are never stored, so always called so get these presets
int init_paramspresets(PARAMETERPARLIST) // probably never need to modify this
{
  SFTYPE ftemp,tvisc;
  SFTYPE numvtime;
  /* BEGIN: Assign global parameters */
  fprintf(log_file,"begin: init_paramspresets ... "); fflush(log_file);

  ////////////// no need to modify/////

  // must be correct for whatever parameter setup or runtype
  if((runtype==2)||(runtype==3) ){
    simplebc=0; // no choice: since boundary conditions will be read in
  }

  if(myid<=0){
    if(TVDLF){
      fprintf(logfull_file,"Using TVDLF Scheme\n");
    }
    else{
      fprintf(logfull_file,"Using ZEUS Scheme\n");
    }

    if(tf<tavgf){
      fprintf(logfull_file,"Final average time is beyond run time\n");
    }
    if(tstart>tavgi){
      fprintf(logfull_file,"Initial average time is before start run time\n");
    }
  }

  if(COMPDIM<=2){ // so all volumes in 3D revert to 2D form
    x3in=0.0;
    x3out=1.0;
  }

  // should never modify the below
  // x1
  L[1][1]=  x1in;
  //L[1][1]=  -1.0;
  L[2][1]= (x1out-x1in) ; // width of x1-direction
  //L[2][1]= 2.0;

  // x2
  // remember to set reflectix2=1 if you set below to 0 for spc.
  L[1][2]= x2in ; // x2 inner edge position
  //L[1][2]= -1.0 ;
  L[2][2]= (x2out-x2in); // width of x2-direction
  //L[2][2]= 2.0;

  // x3
  L[1][3]= x3in ; // x3 inner edge position
  L[2][3]= (x3out-x3in) ; // width of x3-direction

  if(BOUNDTYPE==1){
    ///
    // both inner and outer must be 5 above!
    periodicx1=((bcix1==5)||(bcox1==5) ? 1 : 0); // 0: not periodic on x1-dir 1: is
    // if 1, implies global reflection on that boundary
    reflectix1=((bcix1==1)||(bcix1==2) ? 1 : 0); // if r=0 is inner edge, so reflecting in spc
    reflectox1=((bcox1==1)||(bcox1==2) ? 1 : 0); // outer radial reflection
    // skips used to avoid computing inner edge values for that direction.  Needed to ensure compute all relevant stuff and only relevant stuff.  Generally 0 for periodic and 1 otherwise.
    skipix1=((periodicx1==1) ? 0 : 1); // what zone to start at on inner edge
    // allows not to compute inner r-edge zones.  Set to 0 if not reflecting at r=0 OR 1 if reflecting at r=0 or if don't want to calculate inner r-edge zone because it's a boundary zone(e.g. outflow)
    
    periodicx2=((bcix2==5)||(bcox2==5) ? 1 : 0);
    // if 1, implies global reflection on that boundary
    reflectix2=((bcix2==1)||(bcix2==2) ? 1 : 0); // as above but with x2 dir
    reflectox2=((bcox2==1)||(bcox2==2) ? 1 : 0); // as above but with x2 dir
    skipix2=((periodicx2==1) ? 0 : 1); // what zone to start at on inner edge of x2-grid
    // skip inner x2 edge when is boundary zone.
    
    periodicx3=((bcix3==5)||(bcox3==5) ? 1 : 0);
    reflectix3=((bcix3==1)||(bcix3==2) ? 1 : 0);
    reflectox3=((bcox3==1)||(bcox3==2) ? 1 : 0);
    skipix3=((periodicx3==1) ? 0 : 1);


    // overrides for limited dimensions (assumes valid choice for grid edge positions(e.g. N2==1 and reflecti/ox2 makes no sense with location of v2)
    if(N1==1){ skipix1=0;}
    if(N2==1){ skipix2=0;}
    if(N3==1){ skipix3=0;}

  }


  // compute surfaces if BOUNDTYPE==2
  if((BOUNDTYPE==2)||(BOUNDTYPE==3)){
    // still care about skips, but used differently by init_bc_gen now when BOUNDTYPE>1
    // start out at 1
    skipix1=1;
    skipix2=1;
    skipix3=1;
    // assume don't care about these anymore
    intix1=0;
    intox1=N1;
    intix2=0;
    intox2=N2;
    intix3=0;
    intox3=N3;
    periodicx2special=0;
  }

  // no point in doing transport if 1-D in that direction
  if(N1==1){ transx1=0;}
  if(N2==1){ transx2=0;}
  if(N3==1){ transx3=0;}

  ///
  ////////////// no need to modify/////


  // probably don't need to modify below
  if(POSTPROC==0){
    if(runtype==0){
      deleteolddat=1; // 0: keep old data 1: delete them
      if(IMAGEFLAG==1){
	deleteoldimg=1;  // 0: keep old images 1: delete them
      }
      else{
	deleteoldimg=0;  // assume when don't want to make, want to keep
      }
      if(DOPARDIAG==1){
	deleteoldpar=1; // 0: keep old par files 1: delete them
      }
      else deleteoldpar=0; // assume when don't want to make that want to keep
    }
    else{
      deleteolddat=0; // assume want to keep when reentrant
      deleteoldimg=0;
      deleteoldpar=0;
    }
  }
  else{
    deleteolddat=0; // keep the data when pp!!
    deleteoldimg=0;
    deleteoldpar=0;
  }
  // assume that directinput>0 means normally want to continue ener/loss data sets too, and that directinput==0 means starting fresh from seperate data set as initial data, so no appending wanted
  if(runtype==0){
    appendold=0; // 0: no append time series data file 1: do append
  }
  else{
    if(directinput>0){
      // assume want to append if doing reentrance of any kind, unless specified  not to
    }
    else{
      appendold=0; // assume no append wanted when inputting files using directinput==0
    }
  }

  /* END: global params */
  fprintf(log_file,"end: init_paramspresets\n"); fflush(log_file);
  return(0);
}

// 1/2/3-d valid
int init_reentrance2(SFTYPE time,int which)
{
  int ifresult1,ifresult2;
  fprintf(log_file,"begin: init_reentrance2 ... "); fflush(log_file);

  // check if close to dump
  if(DUMPSM==1){
    ifresult1=(fabs(fmod(time-tstart,DTpd))<1E-2);
    ifresult2=(fabs(fmod(time-tstart,DTpd))>0.99);
    //ifresult3=(fabs(fmod(time-tstart,DTpd))==DTpd);
  }
  else{
    ifresult1=(fabs(fmod(time-tstart,DTd))<1E-2);
    ifresult2=(fabs(fmod(time-tstart,DTd))>0.99);
    //ifresult3=(fabs(fmod(time-tstart,DTd))==DTd);
  }
  if(ifresult1||ifresult2){
    if(ifresult1){ // same if ifresult3
      ireenter=1;
      // so close that should keep current and skip to next
    }
    if(ifresult2){
      ireenter=2;
      // so close that should keep that file and skip to next(which after (int) is 2 away)
    }
  }
  else{
    ireenter=1; // overwrite next file
  }

  // must have time to do this, so call after time set
  // this is dump # to start at
  pdump_start=(int)((time-tstart)/DTpd)+ireenter;
  dump_start=(int)((time-tstart)/DTd)+ireenter;
  fldump_start=(int)((time-tstart)/DTfld)+ireenter;
  npdump_start=(int)((time-tstart)/DTd)+ireenter;
  adump_start=(int)((time-tstart)/DTd)+ireenter;
  floor_start=(int)((time-tstart)/DTfloor)+ireenter;
  image_start=(int)((time-tstart)/DTi)+ireenter;
  
  if(ADUMPFLAG==-1){
    adump_start=0; // force since only 1
  }
  if(PDUMPFLAG==0) pdump_start=0;
  if(DUMPFLAG==0) dump_start=0;
  if(FLDUMPFLAG==0) fldump_start=0;
  if(NPDUMPFLAG==0) npdump_start=0;
  if(FLOORDUMPFLAG==0) floor_start=0;
  if(ADUMPFLAG==0) adump_start=0;
  if(IMAGEFLAG==0) image_start=0;

  if(which==1){ // correct if changed number of files from figured
    // make sure not over file # by ireenter
    if(pdump_start>numpdumps-1+ireenter){ pdump_start=numpdumps-1+ireenter; }
    if(dump_start>numdumps-1+ireenter){ dump_start=numdumps-1+ireenter; }
    if(fldump_start>numfldumps-1+ireenter){ fldump_start=numfldumps-1+ireenter; }
    if(npdump_start>numnpdumps-1+ireenter){ npdump_start=numnpdumps-1+ireenter; }
    if(adump_start>numadumps-1+ireenter){ adump_start=numadumps-1+ireenter; }
    if(floor_start>numfloordumps-1+ireenter){ floor_start=numfloordumps-1+ireenter; }
    if(image_start>numimages-1+ireenter){ image_start=numimages-1+ireenter; }

    // force a start sequence
    if((runtype==1)&&(directinput==1)){
	pdump_start=PDUMP_START;
	numpdumps=pdump_start+1;
	dump_start=DUMP_START;
	numdumps=dump_start+1;
	fldump_start=FLDUMP_START;
	numfldumps=fldump_start+1;
	npdump_start=NPDUMP_START;
	numnpdumps=npdump_start+1;
	adump_start=ADUMP_START;
	numadumps=adump_start+1;
	floor_start=FLOOR_START;
	numfloordumps=floor_start+1;
        image_start = IMAGE_START;
	numimages=image_start+1;
    }

    if(pdump_start<0){ pdump_start=0; }
    if(dump_start<0){ dump_start=0; }
    if(fldump_start<0){ fldump_start=0; }
    if(npdump_start<0){ npdump_start=0; }
    if(adump_start<0){ adump_start=0; }
    if(floor_start<0){ floor_start=0; }
    if(image_start<0){ image_start=0; }
  }

  if(myid<=0){
    fprintf(logfull_file,"starts: %d %d %d %d %d %d\n",pdump_start,dump_start,npdump_start,adump_start,floor_start,image_start);
    fflush(logfull_file);
  }

  fprintf(log_file,"end: init_reentrance2\n"); fflush(log_file);
  return(0);
}

// 1/2/3-d valid
int init_reentrance(void)
{
  FILE *in;
  char temps[MAXFILENAME];

  fprintf(log_file,"begin: init_reentrance ... "); fflush(log_file);

  if(runtype>0){

      // determine number of primitive dumps
      if(PDUMPFLAG){
	  if(myid==0){
	      sprintf(temps,"%s0_numpdumps%s",DATADIR,DAT2EXT);
	      if((in = fopen(temps,"r"))==NULL) {
		  fprintf(fail_file,"error opening dump output file %s\n",temps) ;
		  myexit(1) ;
	      }
	      while(fgetc(in)!='\n');    
	      fscanf(in,"%d",&numpdumps);
	      fclose(in);
	      fprintf(logfull_file,"number of pdumps: %d\n",numpdumps);
	  }
#if(USEMPI)
	  MPI_Bcast(&numpdumps,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
      }
	  
      // determine number of normal dump files
      if(DUMPFLAG){
	  if(myid==0){
	      sprintf(temps,"%s0_numdumps%s",DATADIR,DAT2EXT);
	      if((in = fopen(temps,"r"))==NULL) {
		  fprintf(fail_file,"error opening dump output file %s\n",temps) ;
		  myexit(1) ;
	      }
	      while(fgetc(in)!='\n');    
	      fscanf(in,"%d",&numdumps);
	      fclose(in);
	      fprintf(logfull_file,"number of dumps: %d\n",numdumps);
	  }
#if(USEMPI)
	  MPI_Bcast(&numdumps,1,MPI_INT,0,MPI_COMM_WORLD);
#endif

      }
	  
      // determine number of fieldline dump files
      if(FLDUMPFLAG){
	  if(myid==0){
	      sprintf(temps,"%s0_numfldumps%s",DATADIR,DAT2EXT);
	      if((in = fopen(temps,"r"))==NULL) {
		  fprintf(fail_file,"error opening fldump output file %s\n",temps) ;
		  myexit(1) ;
	      }
	      while(fgetc(in)!='\n');    
	      fscanf(in,"%d",&numfldumps);
	      fclose(in);
	      fprintf(logfull_file,"number of fldumps: %d\n",numfldumps);
	  }
#if(USEMPI)
	  MPI_Bcast(&numfldumps,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
      }
	  
      if(NPDUMPFLAG){
	  if(myid==0){
	      sprintf(temps,"%s0_numnpdumps%s",DATADIR,DAT2EXT);
	      if((in = fopen(temps,"r"))==NULL) {
		  fprintf(fail_file,"error opening dump output file %s\n",temps) ;
		  myexit(1) ;
	      }
	      while(fgetc(in)!='\n');    
	      fscanf(in,"%d",&numnpdumps);
	      fclose(in);
	      
	      fprintf(logfull_file,"number of np dumps: %d\n",numnpdumps);
	  }
#if(USEMPI)
	  MPI_Bcast(&numnpdumps,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
      }
	  
      if(ADUMPFLAG){
	  if(myid==0){
	      sprintf(temps,"%s0_numadumps%s",DATADIR,DAT2EXT);
	      if((in = fopen(temps,"r"))==NULL) {
		  fprintf(fail_file,"error opening dump output file %s\n",temps) ;
		  myexit(1) ;
	      }
	      while(fgetc(in)!='\n');    
	      fscanf(in,"%d",&numadumps);
	      fclose(in);
	      
	      fprintf(logfull_file,"number of adumps: %d\n",numadumps);
	  }
#if(USEMPI)
	  MPI_Bcast(&numadumps,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
      }
	  
      if(FLOORDUMPFLAG){
	  if(myid==0){
	      sprintf(temps,"%s0_numfloordumps%s",DATADIR,DAT2EXT);
	      if((in = fopen(temps,"r"))==NULL) {
		  fprintf(fail_file,"error opening dump output file %s\n",temps) ;
		  myexit(1) ;
	      }    
	      while(fgetc(in)!='\n');    
	      fscanf(in,"%d",&numfloordumps);
	      fclose(in);
	      
	      fprintf(logfull_file,"number of floor dumps: %d\n",numfloordumps);
	  }
#if(USEMPI)
	  MPI_Bcast(&numfloordumps,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
      }

    if(IMAGEFLAG){
	if(myid==0){
	    // get # of images
	    sprintf(temps,"%s%s0_numimages%s",DATADIR,IMAGEDIR,DAT2EXT);
	    if((in = fopen(temps,"r"))==NULL) {
		fprintf(fail_file,"error opening dump output file %s\n",temps) ;
		myexit(1) ;
	    }
	    while(fgetc(in)!='\n'); // skip comment line
	    fscanf(in,"%d",&numimages);
	    fclose(in);
	    
	    fprintf(logfull_file,"number of image dumps: %d\n",numimages);
	}
#if(USEMPI)
	MPI_Bcast(&numimages,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
    }

    // first guess
    pdump_start=numpdumps-1;
    dump_start=numdumps-1;
    fldump_start=numfldumps-1;
    npdump_start=numnpdumps-1;
    adump_start=numadumps-1;
    floor_start=numfloordumps-1;
    image_start=numimages-1;
    ireenter=1;

  }
  else{
  // only exactly valid for good DT if match t on dump and images, and only can do that if dumps occur at image dump interval.
  //So dump DT needs to be integral of image DT.

    pdump_start=0;
    dump_start=npdump_start=0;
    fldump_start=0;
    adump_start=0;
    floor_start=0;
    image_start=0;
    ireenter=0;
  }


  /* END: reentrance params */
  fflush(log_file);
  if(myid<=0){
    fflush(logfull_file);
  }

  fprintf(log_file,"end: init_reentrance\n"); fflush(log_file);
  return(0);
}



int init_otherparams(void) // numerical stuff not stored to gparm file
{

  fprintf(log_file,"begin: init_otherparams ... \n"); fflush(log_file);


  if( (runtype==3)||(runtype==1)){
    // t set by input of data
  }
  else{
    t      = tstart ; // start time
  }
  // compute things that always need but formula doesn't change
  invcour= 1.0/cour; // inverse of courant number
  invcour2= 1.0/cour2; // inverse of diffuse courant number
  DTavg  = (tavgf-tavgi)/((SFTYPE)(numavg));
  wgam53 = (int)(fabs(gam - 5./3.) < ERR); // switch to see if gam==5/3
  wgam1   = (int)(fabs(gam - 1.) < ERR); // switch to see if gam==1
  wgam   = (int)(fabs(gam - 1.) > ERR); // switch to see if need to compute pressure
  wpw = (int)(fabs(rgp - 0.) > ERR); // switch to see if newtonian or PW pot
  
  // here since don't want to dump and still want reentrance
  IOBound[0]=0.1; // first transition from out to inflow
  IOBound[1]=0.9; // second transition from inflow to outflow
  IOBound[2]=M_PI-IOBound[1]; // third transition from inflow to outflow
  IOBound[3]=M_PI-IOBound[0]; // fourth transition from inflow to outflow

  // defaults, keep low, must be lower than problem demands
  IEFLOOR=1E-30;
  DENSITYFLOOR=1E-30;
  // something very demanding in general
  DENSITYFRACFLOOR=IEFRACFLOOR=10.0*NUMEPSILON;

  fprintf(log_file,"end: init_otherparams\n"); fflush(log_file);
  return(0);
}


// periodicx?, skip's, and int's only used by old loop method(LOOPTYPE==1)
int init_mainbc(int px1,int six1,int rix1,int rox1,int px2,int six2,int rix2,int rox2,int px3,int six3,int rix3,int rox3)
{
  // intix2 stuff not really right if account region goes inside inner cpus, but this should never be done anyways!
  int i,j,k;
  int m;

  fprintf(log_file,"begin: init_mainbc ... "); fflush(log_file);
  //x1
  periodicx1=px1; // global condition, figured out in bound.c
  // sets inner and outer i(x1) boundary
  if(mycpupos[1]==0){ // then inner i boundary
    intix1=intix1;
    skipix1=six1;
    reflectix1=rix1;
    if(ncpux1>1){ // must check if cpu is also outer i boundary
      reflectox1=0;
      intox1=N1;
      mpiperiodicx1=px1;
    }
    else{
      reflectox1=rox1;
      intox1=intox1;
      mpiperiodicx1=0;
    }
    if(skipix1==1){
      if(intix1<=1){
	skipintix1=1;
      }
      else skipintix1=intix1;
    }
    else skipintix1=intix1;
  }
  else if(mycpupos[1]==ncpux1-1){// then outer i boundary and NOT inner boundary(due to else if and not just if)
    // only reaches here if ncpux1>1, so definitely outer i cpu and not on inner i
    intix1=0;
    intox1=intox1;
    skipix1=0;
    reflectix1=0;
    reflectox1=rox1;
    skipintix1=0;
    mpiperiodicx1=px1; // ncpux1>1 to be here
  }
  else{ // then must be core zone, not on boundary
    intix1=0;
    intox1=N1;
    skipix1=0;
    reflectix1=0;
    reflectox1=0;
    skipintix1=0;
    periodicx1=0;// doesn't need to worry about this
    mpiperiodicx1=0; // ncpux1>1 to be here
  }


  // x2
  periodicx2=px2; // global condition figured out in bound.c

  if(mycpupos[2]==0){ // then inner j boundary
    intix2=intix2;
    skipix2=six2;
    reflectix2=rix2;
    if(ncpux2>1){ // i.e. this cpu is not also on outer j boundary
      reflectox2=0;
      intox2=N2;
      mpiperiodicx2=px2;
    }
    else{
      reflectox2=rox2;
      intox2=intox2;
      mpiperiodicx2=0;
    }
    if(skipix2==1){
      if(intix2<=1){
	skipintix2=1;
      }
      else skipintix2=intix2;
    }
    else skipintix2=intix2;
  }
  else if(mycpupos[2]==ncpux2-1){ // then outer j boundary and NOT j inner boundary
    // only reaches here if ncpux2>1, so definitely not on inner j boundary
    intix2=0; // go all the way inside
    intox2=intox2; // normal
    skipix2=0;
    reflectix2=0;
    reflectox2=rox2;
    skipintix2=0;
    mpiperiodicx2=px2; // must have ncpux2>1 to be here
  }
  else{ // then must be core zone, not on boundary
    intix2=0;
    intox2=N2;
    skipix2=0;
    reflectix2=0;
    reflectox2=0;
    skipintix2=0;
    periodicx2=0;// doesn't need to worry about this
    mpiperiodicx2=0;
  }


  // x3
  periodicx3=px3; // global condition figured out in bound.c

  if(mycpupos[3]==0){ // then inner k boundary
    intix3=intix3;
    skipix3=six3;
    reflectix3=rix3;
    if(ncpux3>1){ // i.e. this cpu is not also on outer k boundary
      reflectox3=0;
      intox3=N3;
      mpiperiodicx3=px3;
    }
    else{
      reflectox3=rox3;
      intox3=intox3;
      mpiperiodicx3=0;
    }
    if(skipix3==1){
      if(intix3<=1){
	skipintix3=1;
      }
      else skipintix3=intix3;
    }
    else skipintix3=intix3;
  }
  else if(mycpupos[3]==ncpux3-1){ // then outer k boundary and NOT k inner boundary
    // only reaches here if ncpux2>1, so definitely not on inner k boundary
    intix3=0; // go all the way inside
    intox3=intox3; // normal
    skipix3=0;
    reflectix3=0;
    reflectox3=rox3;
    skipintix3=0;
    mpiperiodicx3=px3; // must have ncpux3>1 to be here
  }
  else{ // then must be core zone, not on boundary
    intix3=0;
    intox3=N3;
    skipix3=0;
    reflectix3=0;
    reflectox3=0;
    skipintix3=0;
    periodicx3=0;// doesn't need to worry about this
    mpiperiodicx3=0;
  }

  // output setup
  fprintf(log_file,"x1: %d %d %d %d %d %d\n",intix1,intox1,skipix1,reflectix1,reflectox1,skipintix1);
  fprintf(log_file,"x2: %d %d %d %d %d %d\n",intix2,intox2,skipix2,reflectix2,reflectox2,skipintix2);
  fprintf(log_file,"x3: %d %d %d %d %d %d\n",intix3,intox3,skipix3,reflectix3,reflectox3,skipintix3);
  fprintf(log_file,"per: %d %d %d\n",periodicx1,periodicx2,periodicx3);
  fprintf(log_file,"mpiper: %d %d %d\n",mpiperiodicx1,mpiperiodicx2,mpiperiodicx3);
  for(m=1;m<=3;m++){
    fprintf(log_file,"mycpupos[%d]: %d\n",m,mycpupos[m]);
    fprintf(log_file,"startpos[%d]: %d\n",m,startpos[m]);
    fprintf(log_file,"endpos[%d]: %d\n",m,endpos[m]);
    fprintf(log_file,"totalsize[%d]: %d\n",m,totalsize[m]);
  }
  for(m=0;m<=3;m++){
    fprintf(log_file,"srdir[%d]: %d\n",m,srdir[m]);
  }

  fprintf(log_file,"end: init_mainbc\n"); fflush(log_file);
  return(0);
}



// for BOUNDTYPE==1 case (old general rectangular format)
// 1/2/3-d valid
int init_bc_rect1(int simple,int ix1,int ox1,int ix2,int ox2,int ix3,int ox3)
{
  int i,j,k,l,m,n;
  int ii;
  int N[3+1];
  int numbc[3+1];
  int skipi[3+1];
  int typei,typeo; // inner/outer temps
  int ranki,ranko; // rank of cpu to get info from for this boundary
  int itemp,itemp2i,itemp2o;
  int looper;


  if(PUREBC>0) return(0); // no need/can't initialize bcs/bcv

  fprintf(log_file,"begin: init_bc_rect1 ... "); fflush(log_file);

  N[1]=N1;
  N[2]=N2;
  N[3]=N3;
  // use optimized idea of boundary
  numbc[1]=N1BND*N1OFF;
  numbc[2]=N2BND*N2OFF;
  numbc[3]=N3BND*N3OFF;
  skipi[1]=skipix1;
  skipi[2]=skipix2;
  skipi[3]=skipix3;

  for(l=1;l<=NUMSCA;l++){
    LOOPF{
      bcs[l][1][k][j][i]=5; // default is periodic for 3D (i.e. 1-zone would be just like periodic)
    }
    LOOP{
      /* Stick here the assignment of boundary conditions for active zones */
      bcs[l][1][k][j][i]=0; // other entries do not matter
    }
  }
  
  for(l=1;l<=REALNUMVEC;l++){
    LOOPF{
      bcv[l][1][k][j][i]=5; // default is periodic for 3D (i.e. 1-zone would be just like periodic)
    }
    LOOP{
      /* Stick here the assignment of boundary conditions for active zones */
      bcv[l][1][k][j][i]=0; // other entries do not matter
    }
  }

  /* Now assign bc for boundary zones--good for normal boundary zones in any situation - just change type from here with conditions */

  // go over all except potential
  for(l=1;l<=NUMSCA-NOBOUNDPOT;l++){
    for(m=1;m<=3;m++){
      for(j=-numbc[3-(4-m)%3];j<N[3-(4-m)%3]+numbc[3-(4-m)%3];j++){
	for(i=-numbc[m%3+1];i<N[m%3+1]+numbc[m%3+1];i++){ 
	  /* m%3+1 gives next 1->2,2->3,3->1
	     3-(4-m)%3 gives previous 1->3,2->1,3->2
	  */
	  if(numbc[m]==0) looper=1; else looper=numbc[m]; // account for degenerate case

	  // but don't want to overwrite comp domain!
	  if(m==3){
	    if((numbc[m]==0)&&(i>=0)&&(i<N1)&&(j>=0)&&(j<N2)) continue;
	  }
	  if(m==2){
	    if((numbc[m]==0)&&(i>=0)&&(i<N1)&&(k>=0)&&(k<N3)) continue;
	  }
	  if(m==1){
	    if((numbc[m]==0)&&(k>=0)&&(k<N3)&&(j>=0)&&(j<N2)) continue;
	  }

	  for(n=0;n<looper;n++){
	    
	    if(numbc[m]==0){ itemp2i=0; itemp2o=0;} else{ itemp2i=-numbc[m]+n; itemp2o=N[m]+n; }

	    if((m==1)&&(N1OFF==1)){ /* Assign over x=const boundaries */
	      if((mycpupos[1]==0)&&(!mpiperiodicx1)){ // then inner i boundary
		if(simple==1) itemp=ix1;
		else itemp=4;

		if(bcs[l][1][j][i][itemp2i]<90){ // otherwise don't change
		  bcs[l][1][j][i][itemp2i]=itemp;
		  bcs[l][2][j][i][itemp2i]=m;
		  bcs[l][3][j][i][itemp2i]=1;
		}
	      }
	      else bcs[l][1][j][i][itemp2i]=99;
	      
	      if((mycpupos[1]==ncpux1-1)&&(!mpiperiodicx1)){ // then outer boundary
		//outer boundary
		if(IOBOUNDARY==0){
		  if(simple==1) itemp=ox1;
		  else itemp=4;
		}
		else{
		  if(x[2][2][j]<IOBound[0]){
		    itemp=3;
		  }
		  else if( (x[2][2][j]>=IOBound[0])&& (x[2][2][j]<IOBound[1]) ){
		    itemp=4;
		  }
		  else if( (x[2][2][j]>=IOBound[1])&& (x[2][2][j]<IOBound[2]) ){
		    itemp=3;
		  }
		  else if( (x[2][2][j]>=IOBound[2])&& (x[2][2][j]<IOBound[3]) ){
		    itemp=4;
		  }
		  else if(x[2][2][j]>=IOBound[3]){
		    itemp=3;
		  }
		}
		if(bcs[l][1][j][i][itemp2o]<90){ // otherwise don't change		
		  bcs[l][1][j][i][itemp2o]=itemp;
		  bcs[l][2][j][i][itemp2o]=m;
		  bcs[l][3][j][i][itemp2o]=-1;
		}
	      }
	      else bcs[l][1][j][i][itemp2o]=99;
	    }
	    if((m==2)&&(N2OFF==1)){  /* Assign over y=const boundaries */
	      
	      // setup for global grid, only need to change this
	      
	      //inner boundary
	      if((mycpupos[2]==0)&&(!mpiperiodicx2)){ // then inner j boundary
		if(simple==1) itemp=ix2;
		else itemp=1;
		if(bcs[l][1][i][itemp2i][j]<90){ // otherwise don't change
		  bcs[l][1][i][itemp2i][j]=itemp;
		  bcs[l][2][i][itemp2i][j]=m;
		  bcs[l][3][i][itemp2i][j]=1;
		}
	      }
	      else bcs[l][1][i][itemp2i][j]=99;

	      if((mycpupos[2]==ncpux2-1)&&(!mpiperiodicx2)){ // then outer j boundary
		//outer boundary
		if(simple==1) itemp=ox2;
		else itemp=1;
		if(bcs[l][1][i][itemp2o][j]<90){ // otherwise don't change
		  bcs[l][1][i][itemp2o][j]=itemp;
		  bcs[l][2][i][itemp2o][j]=m;
		  bcs[l][3][i][itemp2o][j]=-1;
		}
	      }
	      else bcs[l][1][i][itemp2o][j]=99;
	    }
	    if((m==3)&&(N3OFF==1)){ /* Assign over z=const boundaries */

	      //inner boundary
	      if((mycpupos[3]==0)&&(!mpiperiodicx3)){ // then inner k boundary
		if(simple==1) itemp=ix3;
		else itemp=1;
		if(bcs[l][1][itemp2i][j][i]<90){ // otherwise don't change
		  bcs[l][1][itemp2i][j][i]=itemp;
		  bcs[l][2][itemp2i][j][i]=m;
		  bcs[l][3][itemp2i][j][i]=1;
		}
	      }
	      else bcs[l][1][itemp2i][j][i]=99;

	      if((mycpupos[3]==ncpux3-1)&&(!mpiperiodicx3)){ // then outer k boundary
		//outer boundary
		if(simple==1) itemp=ox3;
		else itemp=1;
		if(bcs[l][1][itemp2o][j][i]<90){ // otherwise don't change
		  bcs[l][1][itemp2o][j][i]=itemp;
		  bcs[l][2][itemp2o][j][i]=m;
		  bcs[l][3][itemp2o][j][i]=-1;
		}
	      }
	      else bcs[l][1][itemp2o][j][i]=99;
	    }
	  }
	}
      }
    }
  }

  if(NOBOUNDPOT==1){
    // potential is static everywhere for all time
    l=3;
    LOOPF{
      /* Stick here the assignment of boundary conditions for active zones */
      bcs[l][1][k][j][i]=0; // other entries do not matter
    }
  }


  // now do vectors
  for(l=1;l<=NUMVEC;l++){
    for(m=1;m<=3;m++){
      for(j=-numbc[3-(4-m)%3];j<N[3-(4-m)%3]+numbc[3-(4-m)%3];j++){
	for(i=-numbc[m%3+1];i<N[m%3+1]+numbc[m%3+1];i++){ 
	  /* m%3+1 gives next 1->2,2->3,3->1
	     3-(4-m)%3 gives previous 1->3,2->1,3->2
	  */
	  if(numbc[m]==0) looper=1; else looper=numbc[m]; // account for degenerate case

	  // but don't want to overwrite comp domain!
	  if(m==3){
	    if((numbc[m]==0)&&(i>=0)&&(i<N1)&&(j>=0)&&(j<N2)) continue;
	  }
	  if(m==2){
	    if((numbc[m]==0)&&(i>=0)&&(i<N1)&&(k>=0)&&(k<N3)) continue;
	  }
	  if(m==1){
	    if((numbc[m]==0)&&(k>=0)&&(k<N3)&&(j>=0)&&(j<N2)) continue;
	  }

	  for(n=0;n<looper;n++){

	    if(numbc[m]==0){ itemp2i=0; itemp2o=0;} else{ itemp2i=-numbc[m]+n; itemp2o=N[m]+n; }

	    if((m==1)&&(N1OFF==1)){ /* Assign over x=const boundaries */

	      if((mycpupos[1]==0)&&(!mpiperiodicx1)){ // then inner i boundary  
		//inner boundary
		if(simple==1) itemp=ix1;
		else itemp=4;	 
		if(bcv[l][1][j][i][itemp2i]<90){ // otherwise don't change
		  bcv[l][1][j][i][itemp2i]=itemp;
		  bcv[l][2][j][i][itemp2i]=m;
		  bcv[l][3][j][i][itemp2i]=1;
		}
	      }
	      else bcv[l][1][j][i][itemp2i]=99;

	      if((mycpupos[1]==ncpux1-1)&&(!mpiperiodicx1)){
		//outer boundary
		if(IOBOUNDARY==0){	      
		  if(simple==1) itemp=ox1;
		  else itemp=4;
		}
		else{
		  if(x[2][2][j]<IOBound[0]){
		    itemp=3;
		  }
		  else if( (x[2][2][j]>=IOBound[0])&& (x[2][2][j]<IOBound[1]) ){
		    itemp=4;
		  }
		  else if( (x[2][2][j]>=IOBound[1])&& (x[2][2][j]<IOBound[2]) ){
		    itemp=3;
		  }
		  else if( (x[2][2][j]>=IOBound[2])&& (x[2][2][j]<IOBound[3]) ){
		    itemp=4;
		  }
		  else if(x[2][2][j]>=IOBound[3]){
		    itemp=3;
		  }
		}
		if(bcv[l][1][j][i][itemp2o]<90){ // otherwise don't change
		  bcv[l][1][j][i][itemp2o]=itemp;
		  bcv[l][2][j][i][itemp2o]=m;
		  bcv[l][3][j][i][itemp2o]=-1;
		}
	      }
	      else bcv[l][1][j][i][itemp2o]=99;
	    }
	    if((m==2)&&(N2OFF==1)){  /* Assign over y=const boundaries */

	      // setup for global grid, only need to change this

	      if((mycpupos[2]==0)&&(!mpiperiodicx2)){ // then inner j boundary 
		if(simple==1) itemp=ix2;
		else itemp=1;	      
		if(bcv[l][1][i][itemp2i][j]<90){ // otherwise don't change
		  bcv[l][1][i][itemp2i][j]=itemp;
		  bcv[l][2][i][itemp2i][j]=m;
		  bcv[l][3][i][itemp2i][j]=1;
		}
	      }
	      else bcv[l][1][i][itemp2i][j]=99;

	      if((mycpupos[2]==ncpux2-1)&&(!mpiperiodicx2)){
		//outer boundary
		if(simple==1) itemp=ox2;
		else itemp=1;	
		if(bcv[l][1][i][itemp2o][j]<90){ // otherwise don't change
		  bcv[l][1][i][itemp2o][j]=itemp;
		  bcv[l][2][i][itemp2o][j]=m;
		  bcv[l][3][i][itemp2o][j]=-1;
		}
	      }
	      else bcv[l][1][i][itemp2o][j]=99;
	    }
	    if((m==3)&&(N3OFF==1)){ /* Assign over z=const boundaries */

	      if((mycpupos[3]==0)&&(!mpiperiodicx3)){ // then inner k boundary 
		if(simple==1) itemp=ix3;
		else itemp=1;	
		if(bcv[l][1][itemp2i][j][i]<90){ // otherwise don't change
		  bcv[l][1][itemp2i][j][i]=itemp;
		  bcv[l][2][itemp2i][j][i]=m;
		  bcv[l][3][itemp2i][j][i]=1;
		}
	      }
	      else bcv[l][1][itemp2i][j][i]=99;

	      if((mycpupos[3]==ncpux3-1)&&(!mpiperiodicx3)){
		//outer boundary
		if(simple==1) itemp=ox3;
		else itemp=1;
		if(bcv[l][1][itemp2o][j][i]<90){ // otherwise don't change
		  bcv[l][1][itemp2o][j][i]=itemp;
		  bcv[l][2][itemp2o][j][i]=m;
		  bcv[l][3][itemp2o][j][i]=-1;
		}
	      }
	      else bcv[l][1][itemp2o][j][i]=99;

	    }
	  }
	}
      }
    }
  }

  fprintf(log_file,"end: init_bc_rect1\n"); fflush(log_file);
  return(0);
}



#define MAXCOPY 100
// for BOUNDTYPE==2 || 3 case
// 1/2/3-d valid
int init_bc_gen1(int simple,int ix1,int ox1,int ix2,int ox2,int ix3,int ox3)
{
  int skipfactor[2];
  int gotitinner,gotitouter;
  int firsttime;
  int i,j,k,l,m,n;
  int ii,jj,kk;
  int N[3+1];
  int not1[3+1];
  int numbc[3+1];
  int skipi[3+1];
  int typei,typeo; // inner/outer temps
  int ranki,ranko; // rank of cpu to get info from for this boundary
  int itemp,itemp_full;
  SFTYPE radius;
  SFTYPE Rin,Rout;
  SFTYPE radiusa,radiusb,thetaa,thetab,phia,phib;
  int nowhere,gotcompzone,inneri,innerj,innerk,outeri,outerj,outerk,SIZE;
  int nowherelimits[3+1][2]; // 0: in 1:out
  FTYPE (*diffmask)[N2M][N1M];
  FTYPE (*debugmask)[N2M][N1M];
  int limits[NUMINDEX][3+1][2]; // NUMINDEX cases per direction, 3 directions, inner and outer
  int lasti,lastj,lastk,goodi,goodj,goodk;
  int BCTYPEIN,BCTYPEOUT;
  int DONOWHERE,RANGE;
  int otherii,otherjj,otherkk;
  int zonetocopy[MAXCOPY][3+1]; // [0][0] holds # of zones to copy, rest hold zone#'s position to copy, upto 26 zones to copy normally, or upto some larger number for k,j,i=-2 in velocity
  int whichzone;
  int GRIDLINEFIRST,GRIDLINEONLY;
  int ininflag,outinflag;
  short (**mask)[N2M][N1M];
  int iisum,jjsum,kksum;
  int countfaces[3];
  int maskfix[2];

  fprintf(log_file,"begin: init_bc_gen1 ... "); fflush(log_file);

  diffmask=work1;
  debugmask=work2;

  // initialize mask memory
  if( (mask=(short (**)[N2M][N1M])malloc(sizeof(short*)*NUMINDEX))==NULL){
    fprintf(stdout,"cannot initialze mask memory\n");
    myexit(1);
  }
  for(i=0;i<NUMINDEX;i++){
    if( (mask[i]=(short (*)[N2M][N1M])malloc(sizeof(short)*N3M*N2M*N1M))==NULL){
      fprintf(stdout,"cannot initialze mask memory\n");
      myexit(1);
    }
    else{
      // shift pointer
      mask[i]=(short (*) [N2M][N1M])(&(mask[i][N3BND][N2BND][N1BND]));
    }
  }

  bzs=(short***)(malloc(sizeof(short**)*NUMINDEX));
  if(bzs==NULL){
    fprintf(fail_file,"couldn't allocate bzs array\n");
    myexit(1);
  }
  




  N[1]=N1;
  N[2]=N2;
  N[3]=N3;
  not1[1]=N1NOT1;
  not1[2]=N2NOT1;
  not1[3]=N3NOT1;

  numbc[1]=N1BND;
  numbc[2]=N2BND;
  numbc[3]=N3BND;
  // below used according to boundary type
  skipi[1]=skipix1;
  skipi[2]=skipix2;
  skipi[3]=skipix3;

  DONOWHERE=1; // whether to completely avoid those zones that should never be use in calculation
  GRIDLINEFIRST=1; // whether to check along grid lines first for samples to take for boundary zone (little less samples)
  GRIDLINEONLY=1; // whether to only check nowhereness along grid lines(for flux)

  if(analoutput==4){
    Rin=Rinner;  // with PW pot must have Rin>rgp=rg
    Rout=Router;
    BCTYPEIN=3; // could be 4
    BCTYPEOUT=3;
  }
  if(analoutput==5){
    Rin=Rinner;  // with PW pot must have Rin>rgp=rg
    Rout=Router;
    BCTYPEIN=4;
    BCTYPEOUT=4;
  }
  if(analoutput==7){
    Rin=0; // not really used
    // tests
    //    Rout=x1out*3;
    //    Rout=x1out*0.9;
    Rout=Router;
    BCTYPEIN=4;
    BCTYPEOUT=4;
  }
  // below tells routines below whether to compute inner velocity edges or not based upon skipi[x] where x is the direction of interest(1,2,3)
  if(BCTYPEIN==5) skipfactor[0]=0; else skipfactor[0]=1;
  if(BCTYPEOUT==5) skipfactor[1]=0; else skipfactor[1]=1;
  

  // form the mask for the normal boundary zones using coordinates as basis for creation
  LOOPF{ 
    if(COORD==3){
      radiusa=x[1][1][i];
      radiusb=x[2][1][i];
    }
    if(COORD==1){
      // r=0 is on the corner if N1,N2,N3 is even
      radiusa=sqrt(x[1][1][i]*x[1][1][i]+x[1][2][j]*x[1][2][j]+x[1][3][k]*x[1][3][k]);
      radiusb=sqrt(x[2][1][i]*x[2][1][i]+x[2][2][j]*x[2][2][j]+x[2][3][k]*x[2][3][k]);
    }

    // radiusb used for scalar location reference
    if(analoutput==4){
      if(radiusb<Rin) mask[0][k][j][i]=BCTYPEIN;
      else if(radiusb>Rout) mask[0][k][j][i]=BCTYPEOUT;
      else mask[0][k][j][i]=0;

      // inside only
      if(radiusb<Rin) mask[6][k][j][i]=BCTYPEIN;
      else mask[6][k][j][i]=0; // computational zone

      // outside only
      if(radiusb>Rout) mask[7][k][j][i]=BCTYPEOUT;
      else mask[7][k][j][i]=0; // computational zone
    }
    // radiusb used for scalar location reference
    if(analoutput==5){ // in-out same
      // 0 mask
      if((radiusb<Rin)||(radiusb>Rout)) mask[0][k][j][i]=BCTYPEIN; // analoutput==5
      //if(radiusb<Rin) mask[0][k][j][i]=BCTYPE; // analoutput==5
      else mask[0][k][j][i]=0; // computational zone

      // inside only
      if(radiusb<Rin) mask[6][k][j][i]=BCTYPEIN; // analoutput==5
      else mask[6][k][j][i]=0; // computational zone

      // outside only
      if(radiusb>Rout) mask[7][k][j][i]=BCTYPEIN; // analoutput==5
      else mask[7][k][j][i]=0; // computational zone
    }
    if(analoutput==7){
      // 0 and outside are same
      if((radiusb>Rout)) mask[7][k][j][i]=mask[0][k][j][i]=BCTYPEOUT; // analoutput==7 or something else
      else mask[7][k][j][i]=mask[0][k][j][i]=0; // computational zone

      // inside only (i.e. no inside)
      mask[6][k][j][i]=0;
    }
  }

  if((analoutput==4)||(analoutput==5)){ // bondi or tori
    itemp=0;
    LOOPF{
      if(mask[6][k][j][i]!=0) itemp++;
    }
    if(numprocs>1){
#if(USEMPI)
      MPI_Allreduce(&itemp,&itemp_full,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
    }
    else itemp_full=itemp;
    // make sure at least 2^DIM zones inside Rin
    if(itemp_full<pow(2,N1NOT1+N2NOT1+N3NOT1)){
      fprintf(fail_file,"Not enough zones to cover central region.  Got %d, need at least %d\n",itemp_full,(int)pow(2,COMPDIM));
      myexit(1);
    }
  }

  // should adjust factors in below 2 processes(-1, 2 ,etc) for other masks for other loops
  // and reduce [1] mask

  // now make sure 2 outer boundary zones, and assign according to if real or MPI boundary zone
  // assumes nowhere NOT set yet
  
  for(l=0;l<=7;l++){ // only for 0,6,7
    if(l==1) l=6;
    LOOPF{
      // INTERIORS
      if((i<=-1)&&(k>=0)&&(k<=N3-1)&&(j>=0)&&(j<=N2-1)){ // INT#1 INT(i1,j,k)
	if(mycpupos[1]==0){
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i>=N1)&&(k>=0)&&(k<=N3-1)&&(j>=0)&&(j<=N2-1)){ // INT#2 INT(o1,j,k)
	if(mycpupos[1]==ncpux1-1){
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i>=0)&&(i<=N1-1)&&(k>=0)&&(k<=N3-1)&&(j<=-1)){ // INT#3 INT(i2,i,k)
	if(mycpupos[2]==0){
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i>=0)&&(i<=N1-1)&&(k>=0)&&(k<=N3-1)&&(j>=N2)){ // INT#4 INT(o2,i,k)
	if(mycpupos[2]==ncpux2-1){
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i>=0)&&(i<=N1-1)&&(k<=-1)&&(j>=0)&&(j<=N2-1)){ // INT#5 INT(i3,i,j)
	if(mycpupos[3]==0){
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i>=0)&&(i<=N1-1)&&(k>=N3)&&(j>=0)&&(j<=N2-1)){ // INT#6 INT(o3,i,j)
	if(mycpupos[3]==ncpux3-1){
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      // NOW START BARS
      else if((i>=0)&&(i<=N1-1)&&(k<=-1)&&(j<=-1)){ // BAR(1,i2,i3) - 1
	if((mycpupos[2]==0)&&(mycpupos[3]==0)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i>=0)&&(i<=N1-1)&&(k>=N3)&&(j<=-1)){ // BAR(1,i2,o3) - 2
	if((mycpupos[2]==0)&&(mycpupos[3]==ncpux3-1)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i>=0)&&(i<=N1-1)&&(k<=-1)&&(j>=N2)){ // BAR(1,o2,i3) - 3
	if((mycpupos[2]==ncpux2-1)&&(mycpupos[3]==0)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i>=0)&&(i<=N1-1)&&(k>=N3)&&(j>=N2)){ // BAR(1,o2,o3) - 4
	if((mycpupos[2]==ncpux2-1)&&(mycpupos[3]==ncpux3-1)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((j>=0)&&(j<=N2-1)&&(k<=-1)&&(i<=-1)){ // BAR(2,i1,i3) - 5
	if((mycpupos[1]==0)&&(mycpupos[3]==0)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((j>=0)&&(j<=N2-1)&&(k>=N3)&&(i<=-1)){ // BAR(2,i1,o3) - 6
	if((mycpupos[1]==0)&&(mycpupos[3]==ncpux3-1)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((j>=0)&&(j<=N2-1)&&(k<=-1)&&(i>=N1)){ // BAR(2,o1,i3) - 7
	if((mycpupos[1]==ncpux1-1)&&(mycpupos[3]==0)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((j>=0)&&(j<=N2-1)&&(k>=N3)&&(i>=N1)){ // BAR(2,o1,o3) - 8
	if((mycpupos[1]==ncpux1-1)&&(mycpupos[3]==ncpux3-1)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((k>=0)&&(k<=N3-1)&&(i<=-1)&&(j<=-1)){ // BAR(3,i1,i2) - 9
	if((mycpupos[1]==0)&&(mycpupos[2]==0)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((k>=0)&&(k<=N3-1)&&(i<=-1)&&(j>=N2)){ // BAR(3,i1,o2) - 10
	if((mycpupos[1]==0)&&(mycpupos[2]==ncpux2-1)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((k>=0)&&(k<=N3-1)&&(i>=N1)&&(j<=-1)){ // BAR(3,o1,i2) - 11
	if((mycpupos[1]==ncpux1-1)&&(mycpupos[2]==0)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((k>=0)&&(k<=N3-1)&&(i>=N1)&&(j>=N2)){ // BAR(3,o1,o2) - 12
	if((mycpupos[1]==ncpux1-1)&&(mycpupos[2]==ncpux2-1)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      // NOW START CORNERS
      else if((i<=-1)&&(j<=-1)&&(k<=-1)){ // CORNER(i1,i2,i3) - 1
	if((mycpupos[1]==0)&&(mycpupos[2]==0)&&(mycpupos[3]==0)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i<=-1)&&(j<=-1)&&(k>=N3)){ // CORNER(i1,i2,o3) - 2
	if((mycpupos[1]==0)&&(mycpupos[2]==0)&&(mycpupos[3]==ncpux3-1)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i<=-1)&&(j>=N2)&&(k<=-1)){ // CORNER(i1,o2,i3) - 3
	if((mycpupos[1]==0)&&(mycpupos[2]==ncpux2-1)&&(mycpupos[3]==0)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i<=-1)&&(j>=N2)&&(k>=N3)){ // CORNER(i1,o2,o3) - 4
	if((mycpupos[1]==0)&&(mycpupos[2]==ncpux2-1)&&(mycpupos[3]==ncpux3-1)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i>=N1)&&(j<=-1)&&(k<=-1)){ // CORNER(o1,i2,i3) - 5
	if((mycpupos[1]==ncpux1-1)&&(mycpupos[2]==0)&&(mycpupos[3]==0)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i>=N1)&&(j<=-1)&&(k>=N3)){ // CORNER(o1,i2,o3) - 6
	if((mycpupos[1]==ncpux1-1)&&(mycpupos[2]==0)&&(mycpupos[3]==ncpux3-1)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i>=N1)&&(j>=N2)&&(k<=-1)){ // CORNER(o1,o2,i3) - 7
	if((mycpupos[1]==ncpux1-1)&&(mycpupos[2]==ncpux2-1)&&(mycpupos[3]==0)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
      else if((i>=N1)&&(j>=N2)&&(k>=N3)){ // CORNER(o1,o2,o3) - 8
	if((mycpupos[1]==ncpux1-1)&&(mycpupos[2]==ncpux2-1)&&(mycpupos[3]==ncpux3-1)){ // only real if both real, otherwise transfers will assign
	  if(!( (l==6) ) ){ // don't add outer zones to inner boundary loop!
	    mask[l][k][j][i]=BCTYPEOUT;
	  }
	}
	else mask[l][k][j][i]+=90; // always assign to outer boundary the MPI type if not real outer boundary
      }
    }
  }
  

  // GODMARK(commented)
  /*
  l=7;
  LOOPF{
    fprintf(stdout,"k=%d j=%d i=%d mask=%d\n",k,j,i,mask[l][k][j][i]);
  }
  fflush(stdout);
  */
  // now use this mask to setup all the possible loop types(involves contraction/expansion of mask(not unless skip's are on, which we assume for now they aren't), then assignment of loop iterates)
  // first do boundary zone loop
  // contraction must use original mask, then final assignment can be made using subtraction of differential mask.  This decouples the directions
  // differential mask (diffmask[][][])
  // 1: change to 0
  // 0: don't change

  // the limits specify: [0]: how far in to look [1]: how far out to look.
  // If comp-loop: finds a bzone and looks for comp zone, if found then change bzone into a compzone
  // or if bzone-loop: finds a comp-zone and looks for a bzone, if found then change comp-zone into bzone.

  // all possible cases of loops used (value specifies how many more zones in that direction)
  l=0; // LOOPBOUNDGEN (no limits)
  l=6; // scalar in
  l=7; // scalar out

  l=1; // case 1, LOOPFC
  limits[l][1][0]=limits[l][1][1]=limits[l][2][0]=limits[l][2][1]=limits[l][3][0]=limits[l][3][1]=2;
  l=2; // LOOPHC
  limits[l][1][0]=limits[l][1][1]=limits[l][2][0]=limits[l][2][1]=limits[l][3][0]=limits[l][3][1]=1;
  l=3; // LOOPFMHPC
  limits[l][1][0]=limits[l][2][0]=limits[l][3][0]=1;
  limits[l][1][1]=limits[l][2][1]=limits[l][3][1]=2;
  l=4; // LOOPHMFPC
  limits[l][1][0]=limits[l][2][0]=limits[l][3][0]=2;
  limits[l][1][1]=limits[l][2][1]=limits[l][3][1]=1;
  l=5; // LOOPC
  limits[l][1][0]=limits[l][1][1]=limits[l][2][0]=limits[l][2][1]=limits[l][3][0]=limits[l][3][1]=0; // uses original mask

  // this boundary type is also used for fluxes(albeit 1 bzone)
  // LOOPBOUNDV1
  for(l=8;l<=24;l++){
    if(l==10) l=23; // skip to it

    if(l==8) limits[l][1][0]=skipfactor[0]; // global statement about all real boundaries
    else if(l==23) limits[l][1][0]=1; // must get surface flux
    else if(l==9) limits[l][1][0]=skipfactor[1];
    else if(l==24) limits[l][1][0]=1; // must get surface flux
    // if per or not, outer ok
    limits[l][1][1]=0;
    limits[l][2][0]=limits[l][2][1]=0;
    limits[l][3][0]=limits[l][3][1]=0;
  }
  // LOOPBOUNDV2
  for(l=10;l<=26;l++){
    if(l==12) l=25; // skip to it

    limits[l][1][0]=limits[l][1][1]=0;
    if(l==10) limits[l][2][0]=skipfactor[0];
    else if(l==25) limits[l][2][0]=1;
    else if(l==11) limits[l][2][0]=skipfactor[1];
    else if(l==26) limits[l][2][0]=1;
    limits[l][2][1]=0;
    limits[l][3][0]=limits[l][3][1]=0;
  }
  // LOOPBOUNDV3
  for(l=12;l<=28;l++){
    if(l==14) l=27; // skip to it

    limits[l][1][0]=limits[l][1][1]=0;
    limits[l][2][0]=limits[l][2][1]=0;
    if(l==12) limits[l][3][0]=skipfactor[0];
    else if(l==27) limits[l][3][0]=1;
    else if(l==13) limits[l][3][0]=skipfactor[1];
    else if(l==28) limits[l][3][0]=1;
    limits[l][3][1]=0;
  }

  // LOOPBOUNDB1 (always compute edges, so this sets up killing of outer boundary zone since really on edge and should be computed)
  for(l=14;l<=15;l++){
    limits[l][1][0]=1;
    // if per or not, outer ok
    limits[l][1][1]=0;
    limits[l][2][0]=limits[l][2][1]=0;
    limits[l][3][0]=limits[l][3][1]=0;
  }
  // LOOPBOUNDB2
  for(l=16;l<=17;l++){
    limits[l][1][0]=limits[l][1][1]=0;
    limits[l][2][0]=1;
    limits[l][2][1]=0;
    limits[l][3][0]=limits[l][3][1]=0;
  }
  // LOOPBOUNDB3
  for(l=18;l<=19;l++){
    limits[l][1][0]=limits[l][1][1]=0;
    limits[l][2][0]=limits[l][2][1]=0;
    limits[l][3][0]=1;
    limits[l][3][1]=0;
  }

  // comp loops (setups up same condition as LOOPBOUNDB1)
  l=20; // LOOPB1
  limits[l][1][0]=1;
  limits[l][1][1]=0;
  limits[l][2][0]=limits[l][2][1]=0;
  limits[l][3][0]=limits[l][3][1]=0;
  l=21; // LOOPB2
  limits[l][1][0]=limits[l][1][1]=0;
  limits[l][2][0]=1;
  limits[l][2][1]=0;
  limits[l][3][0]=limits[l][3][1]=0;
  l=22; // LOOPB3
  limits[l][1][0]=limits[l][1][1]=0;
  limits[l][2][0]=limits[l][2][1]=0;
  limits[l][3][0]=1;
  limits[l][3][1]=0;



  
  //  l=14; // LOOPINJ (GODMARK) // not done yet

  for(l=0;l<=NUMINDEX-1;l++){ // different cases(except boundary loop itself)
    
    if( (l!=0)&&(l!=6)&&(l!=7) ){
      // initialize differential mask
      LOOPF{
	diffmask[k][j][i]=0;
      }
      
      // expand/contract for this case
      LOOPF{ // from reference of boundary zones
	if(
	   // no <90 for mask check since ok to change to MPI zone to comp zone
	   ((mask[0][k][j][i]>0)&&( ((l>=1)&&(l<=5))||((l>=20)&&(l<=22)) ) )|| // normal loops (not boundary loops)
	   // find comp zone to change to bzone(for v-field which needs 1..N-1 type comp zones
	   // 90 MPI zone is same as comp zone w.r.t. conversions (must then make sure MPI doesn't overwrite <90 zones!)
	   (( (mask[6][k][j][i]==0)||(mask[6][k][j][i]==90))&&((l==8)||(l==10)||(l==12)||(l==23)||(l==25)||(l==27) ))|| // inner boundary loops
	   (( (mask[7][k][j][i]==0)||(mask[7][k][j][i]==90))&&((l==9)||(l==11)||(l==13)||(l==24)||(l==26)||(l==28) ))|| // outer boundary loops
	   // find bzone to change to comp zone (for b-field which needs 0..N type comp zones
	   ((mask[6][k][j][i]>0)&&((l==14)||(l==16)||(l==18)))|| // inner boundary loops
	   ((mask[7][k][j][i]>0)&&((l==15)||(l==17)||(l==19))) // outer boundary loops
	   ){
	  inneri=((i-limits[l][1][0]*N1NOT1>=-N1BND*N1NOT1) ? (i-limits[l][1][0]*N1NOT1) : -N1BND*N1NOT1);
	  innerj=((j-limits[l][2][0]*N2NOT1>=-N2BND*N2NOT1) ? (j-limits[l][2][0]*N2NOT1) : -N2BND*N2NOT1);
	  innerk=((k-limits[l][3][0]*N3NOT1>=-N3BND*N3NOT1) ? (k-limits[l][3][0]*N3NOT1) : -N3BND*N3NOT1);
	  outeri=((i+limits[l][1][1]*N1NOT1<=N1-1+N1BND*N1NOT1) ? (i+limits[l][1][1]*N1NOT1) : N1-1+N1BND*N1NOT1);
	  outerj=((j+limits[l][2][1]*N2NOT1<=N2-1+N2BND*N2NOT1) ? (j+limits[l][2][1]*N2NOT1) : N2-1+N2BND*N2NOT1);
	  outerk=((k+limits[l][3][1]*N3NOT1<=N3-1+N3BND*N3NOT1) ? (k+limits[l][3][1]*N3NOT1) : N3-1+N1BND*N3NOT1);

	  gotcompzone=0;
	  for(kk=innerk;kk<=outerk;kk++){
	    for(jj=innerj;jj<=outerj;jj++){
	      for(ii=inneri;ii<=outeri;ii++){
		if(
		   // comp loops
		   ((mask[0][kk][jj][ii]==0)&&( ((l>=1)&&(l<=5))||((l>=20)&&(l<=22)) ) )||
		   // for v-field and flux
		   // NOT ok for comp zone to change to MPI type boundary zone.  Must compute it so it can be exchanged.
		   ((mask[6][kk][jj][ii]>0)&&(mask[6][kk][jj][ii]<90)&&((l==8)||(l==10)||(l==12)||(l==23)||(l==25)||(l==27) ))|| // inner boundary loops
		   ((mask[7][kk][jj][ii]>0)&&(mask[7][kk][jj][ii]<90)&&((l==9)||(l==11)||(l==13)||(l==24)||(l==26)||(l==28) ))|| // outer boundary loops
		   // for b-field
		   ((mask[6][kk][jj][ii]==0)&&((l==14)||(l==16)||(l==18)))|| // inner boundary loops
		   ((mask[7][kk][jj][ii]==0)&&((l==15)||(l==17)||(l==19))) // outer boundary loops
		   ){
		  /*
		  if(j==N2-1+NBIGBND*N2NOT1){
		    fprintf(stderr,"%d %d %d  %d %d %d  %d %d %d  %d %d %d\n",k,j,i,kk,jj,ii,innerk,innerj,inneri,outerk,outerj,outeri);
		  }
		  */
		  diffmask[k][j][i]=1;
		  gotcompzone=1;
		  // can break now since only needed to know if this *1* boundary zone is close to a comp zone
		  break;
		}
	      }
	      if(gotcompzone) break;
	    }
	    if(gotcompzone) break;
	  }
	}
      }

      // now assign new mask
      LOOPF{
	if(diffmask[k][j][i]==0){
	  if( ((l>=1)&&(l<=5))||((l>=20)&&(l<=22)) ){ // normal loops
	    mask[l][k][j][i]=mask[0][k][j][i];
	  }
	  if((l==8)||(l==10)||(l==12)||(l==14)||(l==16)||(l==18)||(l==23)||(l==25)||(l==27)){ // inner bound
	    mask[l][k][j][i]=mask[6][k][j][i];
	  }
	  if((l==9)||(l==11)||(l==13)||(l==15)||(l==17)||(l==19)||(l==24)||(l==26)||(l==28)){ // outer bound
	    mask[l][k][j][i]=mask[7][k][j][i];
	  }
	}
	else{
	  // expand comp zones
	  if( ((l>=1)&&(l<=5))||((l>=20)&&(l<=22)) ){ mask[l][k][j][i]=0; }
	  // expand bzones
	  else if((l==8)||(l==23))  mask[l][k][j][i]=mask[6][k][j][i-N1NOT1]; 
	  else if((l==9)||(l==24))  mask[l][k][j][i]=mask[7][k][j][i-N1NOT1]; 
	  else if((l==10)||(l==25)) mask[l][k][j][i]=mask[6][k][j-N2NOT1][i]; 
	  else if((l==11)||(l==26)) mask[l][k][j][i]=mask[7][k][j-N2NOT1][i]; 
	  else if((l==12)||(l==27)) mask[l][k][j][i]=mask[6][k-N3NOT1][j][i];
	  else if((l==13)||(l==28)) mask[l][k][j][i]=mask[7][k-N3NOT1][j][i];
	  // contract bzones
	  else if(l==14) mask[l][k][j][i]=0;
	  else if(l==15) mask[l][k][j][i]=0;
	  else if(l==16) mask[l][k][j][i]=0;
	  else if(l==17) mask[l][k][j][i]=0;
	  else if(l==18) mask[l][k][j][i]=0;
	  else if(l==19) mask[l][k][j][i]=0;
	}
      }
    }// if non-scalar(i.e. reference) loop
  }
    
    // do rest separately since want to define masks above and rest is independently done for each mask
  for(l=0;l<=NUMINDEX-1;l++){ // different cases

    if( (l==0)||((l>=6)&&(l<=19))||( (l>=23)&&(l<=28)) ){ // boundary loops only required, meaningless to comp loops
      if(DONOWHERE){ // for symmetry should do after all masks are defined
	// now find nowhere zones
	// nowhere zones are internal boundary zones that are never used.  They have at least 2 true boundary zones till reaching a computational zone for ALL directions
	LOOPF{
	  if(mask[l][k][j][i]>0){ // ok to not have <90 here since ok to change MPI zone into NULL zone
	    
	    // check if ANY computational zone within 2 zones each direction(including corners!), otherwise no-where zone
	    nowhere=1;
	    
	    // was just SIZE=2 here and Newtonian calc was good, only once added field was problem.
	    // Loop does: if I go over to i+-limits from myself(boundary zone) and I find a compzone, I must be somewhere, otherwise nowhere.
	    // so scalars should be 2 in general
	    // velocity should be 2 in general since symmetric
	    // B-field should be 1 both sides in general

	    
	    if((l==0)||(l==6)||(l==7)){ // normal boundaries
	      nowherelimits[1][0]=2;
	      nowherelimits[1][1]=2;
	      nowherelimits[2][0]=2;
	      nowherelimits[2][1]=2;
	      nowherelimits[3][0]=2;
	      nowherelimits[3][1]=2;
	    }
	    if((l==8)||(l==9)){ // v_1 bound type
	      nowherelimits[1][0]=2;
	      nowherelimits[1][1]=2; // was 3, when no problem with MHD.  since when inside boundary, look for 3 zones to right for comp zone before saying nowhere. (shouldn't be necessary if (-2) not needed)
	      nowherelimits[2][0]=2;
	      nowherelimits[2][1]=2;
	      nowherelimits[3][0]=2;
	      nowherelimits[3][1]=2;
	    }
	    if((l==10)||(l==11)){ // v_2 bound type
	      nowherelimits[1][0]=2;
	      nowherelimits[1][1]=2;
	      nowherelimits[2][0]=2;
	      nowherelimits[2][1]=2; // since when inside boundary, look for 3 zones to right for comp zone before saying nowhere. (shouldn't be necessary if (-2) not needed)
	      nowherelimits[3][0]=2;
	      nowherelimits[3][1]=2;
	    }
	    if((l==12)||(l==13)){ // v_3 bound type
	      nowherelimits[1][0]=2;
	      nowherelimits[1][1]=2;
	      nowherelimits[2][0]=2;
	      nowherelimits[2][1]=2;
	      nowherelimits[3][0]=2;
	      nowherelimits[3][1]=2; // since when inside boundary, look for 3 zones to right for comp zone before saying nowhere. (shouldn't be necessary if (-2) not needed)
	    }
	    if((l==14)||(l==15)){ // B_1 bound type
	      nowherelimits[1][0]=1; // since when inside right boundary, only need 1 boundary layer since N is actually computed
	      nowherelimits[1][1]=1; // 2 before: shouldn't need -2 like zone
	      nowherelimits[2][0]=2;
	      nowherelimits[2][1]=2;
	      nowherelimits[3][0]=2;
	      nowherelimits[3][1]=2;
	    }
	    if((l==16)||(l==17)){ // B_2 bound type
	      nowherelimits[1][0]=2;
	      nowherelimits[1][1]=2;
	      nowherelimits[2][0]=1;
	      nowherelimits[2][1]=1;
	      nowherelimits[3][0]=2;
	      nowherelimits[3][1]=2;
	    }
	    if((l==18)||(l==19)){ // B_3 bound type
	      nowherelimits[1][0]=2;
	      nowherelimits[1][1]=2;
	      nowherelimits[2][0]=2;
	      nowherelimits[2][1]=2;
	      nowherelimits[3][0]=1;
	      nowherelimits[3][1]=1;
	    }
	    if((l==23)||(l==24)){ // 1-layer thick boundary for flux
	      nowherelimits[1][0]=1;
	      nowherelimits[1][1]=1;
	      nowherelimits[2][0]=1;
	      nowherelimits[2][1]=1;
	      nowherelimits[3][0]=1;
	      nowherelimits[3][1]=1;
	    }
	    if((l==25)||(l==26)){ // 1-layer thick boundary for flux
	      nowherelimits[1][0]=1;
	      nowherelimits[1][1]=1;
	      nowherelimits[2][0]=1;
	      nowherelimits[2][1]=1;
	      nowherelimits[3][0]=1;
	      nowherelimits[3][1]=1;
	    }
	    if((l==27)||(l==28)){ // 1-layer thick boundary for flux
	      nowherelimits[1][0]=1;
	      nowherelimits[1][1]=1;
	      nowherelimits[2][0]=1;
	      nowherelimits[2][1]=1;
	      nowherelimits[3][0]=1;
	      nowherelimits[3][1]=1;
	    }


	    inneri=((i-nowherelimits[1][0]*N1NOT1>=-N1BND*N1NOT1) ? (i-nowherelimits[1][0]*N1NOT1) : -N1BND*N1NOT1);
	    innerj=((j-nowherelimits[2][0]*N2NOT1>=-N2BND*N2NOT1) ? (j-nowherelimits[2][0]*N2NOT1) : -N2BND*N2NOT1);
	    innerk=((k-nowherelimits[3][0]*N3NOT1>=-N3BND*N3NOT1) ? (k-nowherelimits[3][0]*N3NOT1) : -N3BND*N3NOT1);
	    outeri=((i+nowherelimits[1][1]*N1NOT1<=N1-1+N1BND*N1NOT1) ? (i+nowherelimits[1][1]*N1NOT1) : N1-1+N1BND*N1NOT1);
	    outerj=((j+nowherelimits[2][1]*N2NOT1<=N2-1+N2BND*N2NOT1) ? (j+nowherelimits[2][1]*N2NOT1) : N2-1+N2BND*N2NOT1);
	    outerk=((k+nowherelimits[3][1]*N3NOT1<=N3-1+N3BND*N3NOT1) ? (k+nowherelimits[3][1]*N3NOT1) : N3-1+N3BND*N3NOT1);


	    // just check along grid lines (for flux loops)
	    if((l>=23)&&(l<=28)&&GRIDLINEONLY){
	      kk=k;
	      jj=j;
	      for(ii=inneri;ii<=outeri;ii++){
		//fprintf(stdout,"kk,jj,ii: %d %d %d\n",kk,jj,ii);
		if(mask[l][kk][jj][ii]==0){
		  nowhere=0;
		}
	      }
	      if(nowhere){ // if still not found comp zone
		kk=k;
		ii=i;
		for(jj=innerj;jj<=outerj;jj++){
		  //fprintf(stdout,"kk,jj,ii: %d %d %d\n",kk,jj,ii);
		  if(mask[l][kk][jj][ii]==0){
		    nowhere=0;
		  }
		}
		if(nowhere){ // still?
		  ii=i;
		  jj=j;
		  for(kk=innerk;kk<=outerk;kk++){
		    //fprintf(stdout,"kk,jj,ii: %d %d %d\n",kk,jj,ii);
		    if(mask[l][kk][jj][ii]==0){
		      nowhere=0;
		    }
		  }
		}
	      }
	    }
	    else{
	      for(kk=innerk;kk<=outerk;kk++){
		for(jj=innerj;jj<=outerj;jj++){
		  for(ii=inneri;ii<=outeri;ii++){
		    if(mask[l][kk][jj][ii]==0){
		      nowhere=0;
		      break;
		    }
		  }
		  if(!nowhere) break;
		}
		if(!nowhere) break;
	      }
	    }
	    if(nowhere) mask[l][k][j][i]=-1; // otherwise keep same (ok to set while looping over since ==0 will not change)
	  }
	}
      }
    }
    // find boundaries for indicies(min/max of k,j,i), defines rectangular range of which loops go over at most


    // NOW MASKS are made so can design loops

  
    // now assign loop iterates (only hit indicies if modified "computational" zone)
    firsttime=1;
    numiter[l]=0;
    LOOPF{
      if(
	 ((mask[l][k][j][i]==0)&&( ((l>=1)&&(l<=5))||((l>=20)&&(l<=22))) )||
	 ((mask[l][k][j][i]>0)&&(mask[l][k][j][i]<90)&&( (l==0)||((l>=6)&&(l<=19))||( (l>=23)&&(l<=28) ) ) )
	 // only want to iterate loops over non-boundary and non-mpi zones
	 ){
	numiter[l]++;
      }
    }
    fprintf(log_file,"numiter[%d]=%d\n",l,numiter[l]);

    /*
  if(l==8){
    LOOPF{
      fprintf(log_file,"%d %d %d %d\n",k,j,i,mask[l][k][j][i]);
    }
    
    fflush(log_file);
    // GODMARK
    //  itemp+=init_outgparm(-1); // -1 means out all within define params
    myexit(0);
  }
    */

    if( ((l>=6)&&(l<=19))||( (l>=23)&&(l<=28)) ){ // only needed for boundary loops in bound() code, and flux code
      bzs[l]=(short**)(malloc(sizeof(short*)*numiter[l]*3));
      if(bzs[l]==NULL){
	fprintf(fail_file,"couldn't allocate bzs array: bzs[%d]\n",l);
	myexit(1);
      }
    }
    else{
      // otherwise don't need since not boundary loop
    }

    // allocate array for loop index holder
    // GODMARK force for now to always make indx
    if(1){ // see internal comments
      indx[l]=(short (*))(malloc(sizeof(short)*numiter[l]*3));
      if(indx[l]==NULL){
	fprintf(fail_file,"couldn't allocate # %d index array\n",l);
	myexit(1);
      }
      // now assign actual values
      temptempi=0;
      LOOPF{
	if(
	   ((mask[l][k][j][i]==0)&&( ((l>=1)&&(l<=5))||((l>=20)&&(l<=22))) )|| // comp loops (not needed unless BOUNDTYPE==2)
	   ((mask[l][k][j][i]>0)&&(mask[l][k][j][i]<90)&&(((l>=6)&&(l<=19))||( (l>=23)&&(l<=28)) )) // boundary loops
	   // only want to iterate loops over non-boundary and non-mpi zones
	   ){
	    if(temptempi>numiter[l]){
		fprintf(fail_file,"Beyond defined indx[%d] array: temptempi=%d,numiter[l]=%d : %d %d %d\n",l,temptempi,numiter[l],i,j,k); fflush(fail_file);
	    }
	  indx[l][temptempi*3+0]=i;
	  indx[l][temptempi*3+1]=j;
	  indx[l][temptempi*3+2]=k;
	  temptempi++;
	}
      }
    }
    // GODMARK(commented)
    /*
    for(temptempi=0;temptempi<numiter[l];temptempi++){
      fprintf(stdout,"l=%03d temptempi=%03d i=%03d j=%03d k=%03d\n",l,temptempi,indx[l][temptempi*3+0],indx[l][temptempi*3+1],indx[l][temptempi*3+2]); 
    }
    */    
  
    // now setup BOUNDTYPE==3 stuff
    // currently not setup to use this for boundary/flux things (would need to use mask directly)
    if( (BOUNDTYPE==3)&&(l>0) ){
      if(!DONOWHERE){
	fprintf(fail_file,"must have nowhere on when doing boundtype==3\n");
	myexit(1);
      }
      // use masks to find where, for each j,k, i begins and ends
      LOOPF3 LOOPF2{
	gotitinner=0;
	for(i=-N1BND;i<N1+N1BND;i++){
	  // hunt for negative side of computational boundary
	  if(
	     ( ( ((l>=1)&&(l<=5))||((l>=20)&&(l<=22)) )&&(mask[l][k][j][i]==0) )||
	     ( (((l>=6)&&(l<=19))||( (l>=23)&&(l<=28)))&&(mask[l][k][j][i]>0)&&(mask[l][k][j][i]<90) ) // assumes nowhere is on
	     // only want to iterate loops over non-boundary and non-mpi zones
	     )
	    {
	    iindx[l][0][k][j]=i; // starting point
	    i=N1+N1BND; // force to next k,j
	    gotitinner=1;
	  }
	}
	gotitouter=0;
	for(i=N1+N1BND-1;i>=-N1BND;i--){
	  // hunt for negative side of computational boundary
	  if(
	     ( ( ((l>=1)&&(l<=5))||((l>=20)&&(l<=22)) )&&(mask[l][k][j][i]==0) )||
	     ( (((l>=6)&&(l<=19))||( (l>=23)&&(l<=28)))&&(mask[l][k][j][i]>0)&&(mask[l][k][j][i]<90) ) // assumes nowhere is on
	     // only want to iterate loops over non-boundary and non-mpi zones
	     )
	    { // assumes nowhere is on
	    iindx[l][1][k][j]=i; // ending point
	    i=-N1BND-1; // force to next k,j
	    gotitouter=1;
	  }
	}
	if(!gotitinner){
	  if(gotitouter){
	    fprintf(fail_file,"Bad condition reached where found outer but not inner edge\n");
	    myexit(1);
	  }
	  else{
	    // then a not needed j,k, so set inner/outer interior of loop isn't performed
	    iindx[l][0][k][j]=-1;
	    iindx[l][1][k][j]=-2;
	  }	      
	}
	else{
	  if(!gotitouter){
	    fprintf(fail_file,"Bad condition reached where found inner but not outer edge\n");
	    myexit(1);
	  }
	  else{
	    // then got everything fine
	  }
	}
      }	
    }
  }  // done with all loop iterate issues
  

  // Finally, assign over masks for code use:

  LOOPF{
    bzmask[k][j][i]=mask[0][k][j][i];
    bzmaskin[k][j][i]=mask[6][k][j][i];
  }






  //////////////////////////
  // BOUND SPECIFIC STUFF (l>=6)&&(l<=19)   ||  (l>=23)&(l<=28)  (for flux loops)
  /////////
  // now figure out stuff for bound() using the mask



  // count the number of real computational zones for stats
  realtotalzones=0;
  rtotalzones=0;
  LOOPF{
    if(bzmask[k][j][i]==0) rtotalzones++;
  }
  // now tell rest of cpus the total
#if(USEMPI)
  MPI_Allreduce(&rtotalzones, &realtotalzones, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  realtotalzones=rtotalzones;
#endif

  // From this mask, determine which other zone to copy from for each boundary zone
  // this is purely a bound() code used result


  //l=6; // start (0 just a complete mask)
  if(DONOWHERE) RANGE=NBIGBND;
  else RANGE=NBIGM;
  
  for(l=6;l<=28;l++){
    if(l==20) l=23; // 20-22 aren't bound loops

    // cycle through scalar and vector component boundary copies
    // now use debugmask for some debugging of the boundary exchange info
    LOOPF{
      debugmask[k][j][i]=-2;
    }
    // reset diffmask for flux cases to use to make sure no duplicate comp zones used in flux calc
    LOOPF{
      diffmask[k][j][i]=0;
    }
    
    whichzone=0; // holds which zone among those that are bounded in the order of the loop for boundary zones
    LOOPF{
      //      fprintf(stdout,"shitbag: %d %d %d : mask: %d\n",k,j,i,masktemp[k][j][i]);
      
      if( (mask[l][k][j][i]>0)&&(mask[l][k][j][i]<90) ){ // find boundary zone
	// this part only relevant to real boundary zones, not MPI boundary zones

	// find closet computational zone to copy from(set of them now)
	gotcompzone=0;
	// note the preference by loop order for i copy then j copy then k copy if all else equal.  Assume negligible effect that resolves away
	// don't go beyond grid
	zonetocopy[0][0]=0; // initialize # of zones to copy for this bzone
	if(!(((l>=23)&&(l<=28)))){
	  for(SIZE=1;SIZE<=RANGE;SIZE++){
	    inneri=((i-SIZE*N1NOT1>=-N1BND*N1NOT1) ? (i-SIZE*N1NOT1) : -N1BND*N1NOT1);
	    innerj=((j-SIZE*N2NOT1>=-N2BND*N2NOT1) ? (j-SIZE*N2NOT1) : -N2BND*N2NOT1);
	    innerk=((k-SIZE*N3NOT1>=-N3BND*N3NOT1) ? (k-SIZE*N3NOT1) : -N3BND*N3NOT1);
	    outeri=((i+SIZE*N1NOT1<=N1-1+N1BND*N1NOT1) ? (i+SIZE*N1NOT1) : N1-1+N1BND*N1NOT1);
	    outerj=((j+SIZE*N2NOT1<=N2-1+N2BND*N2NOT1) ? (j+SIZE*N2NOT1) : N2-1+N2BND*N2NOT1);
	    outerk=((k+SIZE*N3NOT1<=N3-1+N3BND*N3NOT1) ? (k+SIZE*N3NOT1) : N3-1+N3BND*N3NOT1);	  
	    //fprintf(stdout,"l=%d SIZE: %d inner/outer: %d %d %d : %d %d , %d %d , %d %d : mask: %d\n",l,SIZE,k,j,i,innerk,outerk,innerj,outerj,inneri,outeri,mask[l][k][j][i]);
	    
	    // first just check along grid lines
	    if(GRIDLINEFIRST){
	      kk=k;
	      jj=j;
	      for(ii=inneri;ii<=outeri;ii++){
		//fprintf(stdout,"kk,jj,ii: %d %d %d\n",kk,jj,ii);
		if((mask[l][kk][jj][ii]==0)&&(diffmask[kk][jj][ii]==0)){
		  zonetocopy[0][0]++; // this holds how many zones to copy, comes first since used below for array value and 1 is first element
		  zonetocopy[zonetocopy[0][0]][1]=ii;
		  zonetocopy[zonetocopy[0][0]][2]=jj;
		  zonetocopy[zonetocopy[0][0]][3]=kk;
		  //fprintf(stdout,"l=%d SIZE: %d : zonetocopy: #=%d k=%d: j=%d i=%d : mask: %d\n",l,SIZE,zonetocopy[0][0],zonetocopy[zonetocopy[0][0]][3],zonetocopy[zonetocopy[0][0]][2],zonetocopy[zonetocopy[0][0]][1],mask[l][k][j][i]);
		  gotcompzone=1;
		  // must get all such zones in this kernel for symmetric average, so don't break
		}
	      }
	      kk=k;
	      ii=i;
	      for(jj=innerj;jj<=outerj;jj++){
		//fprintf(stdout,"kk,jj,ii: %d %d %d\n",kk,jj,ii);
		if((mask[l][kk][jj][ii]==0)&&(diffmask[kk][jj][ii]==0)){
		  zonetocopy[0][0]++; // this holds how many zones to copy, comes first since used below for array value and 1 is first element
		  zonetocopy[zonetocopy[0][0]][1]=ii;
		  zonetocopy[zonetocopy[0][0]][2]=jj;
		  zonetocopy[zonetocopy[0][0]][3]=kk;
		  //fprintf(stdout,"l=%d SIZE: %d : zonetocopy: #=%d k=%d: j=%d i=%d : mask: %d\n",l,SIZE,zonetocopy[0][0],zonetocopy[zonetocopy[0][0]][3],zonetocopy[zonetocopy[0][0]][2],zonetocopy[zonetocopy[0][0]][1],mask[l][k][j][i]);
		  gotcompzone=1;
		  // must get all such zones in this kernel for symmetric average, so don't break
		}
	      }
	      ii=i;
	      jj=j;
	      for(kk=innerk;kk<=outerk;kk++){
		//fprintf(stdout,"kk,jj,ii: %d %d %d\n",kk,jj,ii);
		if((mask[l][kk][jj][ii]==0)&&(diffmask[kk][jj][ii]==0)){
		  zonetocopy[0][0]++; // this holds how many zones to copy, comes first since used below for array value and 1 is first element
		  zonetocopy[zonetocopy[0][0]][1]=ii;
		  zonetocopy[zonetocopy[0][0]][2]=jj;
		  zonetocopy[zonetocopy[0][0]][3]=kk;
		  //fprintf(stdout,"l=%d SIZE: %d : zonetocopy: #=%d k=%d: j=%d i=%d : mask: %d\n",l,SIZE,zonetocopy[0][0],zonetocopy[zonetocopy[0][0]][3],zonetocopy[zonetocopy[0][0]][2],zonetocopy[zonetocopy[0][0]][1],mask[l][k][j][i]);
		  gotcompzone=1;
		  // must get all such zones in this kernel for symmetric average, so don't break
		}
	      }
	    }
	    
	    // now if no comp zones along grid lines, check cross-grid lines
	    if(!gotcompzone){
	      for(kk=innerk;kk<=outerk;kk++){
		for(jj=innerj;jj<=outerj;jj++){
		  for(ii=inneri;ii<=outeri;ii++){
		    //fprintf(stdout,"kk,jj,ii: %d %d %d\n",kk,jj,ii);
		    if((mask[l][kk][jj][ii]==0)&&(diffmask[kk][jj][ii]==0)){		  
		      zonetocopy[0][0]++; // this holds how many zones to copy, comes first since used below for array value and 1 is first element
		      zonetocopy[zonetocopy[0][0]][1]=ii;
		      zonetocopy[zonetocopy[0][0]][2]=jj;
		      zonetocopy[zonetocopy[0][0]][3]=kk;
		      //fprintf(stdout,"l=%d SIZE: %d : zonetocopy: #=%d k=%d: j=%d i=%d : mask: %d\n",l,SIZE,zonetocopy[0][0],zonetocopy[zonetocopy[0][0]][3],zonetocopy[zonetocopy[0][0]][2],zonetocopy[zonetocopy[0][0]][1],mask[l][k][j][i]);
		      gotcompzone=1;
		      // must get all such zones in this kernel for symmetric average, so don't break
		    }
		  }
		}
	      }
	    }
	    if(gotcompzone) break; // can break now since only this layer needed for copy.  More layers may be better or worse, but certainly overly non-local
	  }
	}// done if not flux
	else{ // for flux types
	  gotcompzone=1;
	}
	if(!gotcompzone){
	  fprintf(fail_file,"Never found computational zone to copy from!: l=%d %d %d %d\n",l,k,j,i);
	  myexit(20);
	}
	else{ // good to go
	  // first allocate this boundary zone's copy info and assign over
	  if((l>=23)&&(l<=28)){
	    bzs[l][whichzone]=(short (*))(malloc(sizeof(short)*(1))); // 1: determines flux direction (-1,0,1)
	    if(bzs[l][whichzone]==NULL){
	      fprintf(fail_file,"Couldn't allocate bzs[l][%d]'s memory for 3*%d zones\n",whichzone,zonetocopy[0][0]);
	      myexit(1);
	    }
	  }
	  else{ // need 3 flags for each other zone
	    bzs[l][whichzone]=(short (*))(malloc(sizeof(short)*(3*zonetocopy[0][0]+4))); // 4 is for 0: number of zones, -1: ininflag, -2: outinflag (and 1 more)
	    if(bzs[l][whichzone]==NULL){
	      fprintf(fail_file,"Couldn't allocate bzs[l][%d]'s memory for 3*%d zones\n",whichzone,zonetocopy[0][0]);
	      myexit(1);
	    }
	    bzs[l][whichzone]+=3; //(so -3, -2, -1 hold flags
	    // unlike zonetocopy, bzs[l][] uses only [m][0] (1 element) for # of elements, and no offset for 3-directions
	    bzs[l][whichzone][0]=zonetocopy[0][0]; // go ahead assign before check since NULL check will fail anyways if this assignment is bad
	    //bzs[l][whichzone][0]=1;
	    if(bzs[l][whichzone][0]>MAXCOPY){
	      fprintf(fail_file,"Too many zones?!: l=%d : %d %d %d: %d\n",l,k,j,i,bzs[l][whichzone][0]);
	      myexit(100);
	    }
	  
	    // go over copying info down for this boundary zone
	    for(m=0;m<bzs[l][whichzone][0];m++){
	      ii=zonetocopy[m+1][1]; // zonetocopy[m+1]'s +1 since 0 holds # of zones
	      jj=zonetocopy[m+1][2];
	      kk=zonetocopy[m+1][3];
	      
	      // now adjust this for various boundary condition types	  
	      switch(mask[l][k][j][i]){
	      case 1:
	      case 2:
		// scalar reflector assuming ii'th zone is ii-1's reflective pair
		if((l==6)||(l==7)){
		  if(ii>i){
		    otherii=ii+(ii-i-1);
		  }
		  else if(ii<i){
		    otherii=ii-(i-ii-1);
		  }
		  else if(ii==i){
		    otherii=ii;
		  }
		  if(jj>j){
		    otherjj=jj+(jj-j-1);
		  }
		  else if(jj<j){
		    otherjj=jj-(j-jj-1);
		  }
		  else if(jj==j){
		    otherjj=jj;
		  }
		  if(kk>k){
		    otherkk=kk+(kk-k-1);
		  }
		  else if(kk<k){
		    otherkk=kk-(k-kk-1);
		  }
		  else if(kk==k){
		    otherkk=kk;
		  }
		}
		if((l==8)||(l==9)){
		  if(ii>i){
		    otherii=ii+(ii-i-2);
		  }
		  else if(ii<i){
		    otherii=ii-(i-ii-2);
		  }
		  else if(ii==i){
		    otherii=ii;
		  }
		  if(jj>j){
		    otherjj=jj+(jj-j-1);
		  }
		  else if(jj<j){
		    otherjj=jj-(j-jj-1);
		  }
		  else if(jj==j){
		    otherjj=jj;
		  }
		  if(kk>k){
		    otherkk=kk+(kk-k-1);
		  }
		  else if(kk<k){
		    otherkk=kk-(k-kk-1);
		  }
		  else if(kk==k){
		    otherkk=kk;
		  }
		}
		if((l==10)||(l==11)){
		  if(ii>i){
		    otherii=ii+(ii-i-1);
		  }
		  else if(ii<i){
		    otherii=ii-(i-ii-1);
		  }
		  else if(ii==i){
		    otherii=ii;
		  }
		  if(jj>j){
		    otherjj=jj+(jj-j-2);
		  }
		  else if(jj<j){
		    otherjj=jj-(j-jj-2);
		  }
		  else if(jj==j){
		    otherjj=jj;
		  }
		  if(kk>k){
		    otherkk=kk+(kk-k-1);
		  }
		  else if(kk<k){
		    otherkk=kk-(k-kk-1);
		  }
		  else if(kk==k){
		    otherkk=kk;
		  }
		}
		if((l==12)||(l==13)){
		  if(ii>i){
		    otherii=ii+(ii-i-1);
		  }
		  else if(ii<i){
		    otherii=ii-(i-ii-1);
		  }
		  else if(ii==i){
		    otherii=ii;
		  }
		  if(jj>j){
		    otherjj=jj+(jj-j-1);
		  }
		  else if(jj<j){
		    otherjj=jj-(j-jj-1);
		  }
		  else if(jj==j){
		    otherjj=jj;
		  }
		  if(kk>k){
		    otherkk=kk+(kk-k-2);
		  }
		  else if(kk<k){
		    otherkk=kk-(k-kk-2);
		  }
		  else if(kk==k){
		  otherkk=kk;
		  }
		}
		if((l>=14)&&(l<=19)){
		  // GODMARK -- not setup yet
		  otherii=ii;
		  otherjj=jj;
		  otherkk=kk;
		}
		if((l>=23)&&(l<=28)){
		  // GODMARK -- not setup yet (no change?)
		  otherii=ii;
		  otherjj=jj;
		  otherkk=kk;
		}
		ii=otherii;
		jj=otherjj;
		kk=otherkk;
		break;
	      case 3:
		// not used
		break;
	      case 4:
		// no modification
		// just copy closest one(modulo other details as in bound() or flux calc
		break;
	      case 5:
		// not yet
		break;
	      default:
		myexit(666);
	      }
	      // now assign after any needed adjustments
	      itemp=3;
	      // ii	      
	      bzs[l][whichzone][m*itemp+0+1]=ii; // last +1 in bzs since [whichzone][0] holds # of elements, 
	      // jj
	      bzs[l][whichzone][m*itemp+1+1]=jj;
	      // kk
	      bzs[l][whichzone][m*itemp+2+1]=kk;
	      
	      // debug check (newmask will hold all computational zones used for boundary zones)
	      debugmask[kk][jj][ii]=mask[l][k][j][i];
	      if((l==0)&&(mask[l][kk][jj][ii]!=0)){
		fprintf(stdout,"ERROR: you are assigning a boundary zone a non-computational zone value!: kk/jj/ii: %d %d %d k/j/i: %d %d %d\n",kk,jj,ii,k,j,i);
		myexit(777);
	      }
	      // GODMARK(commented)
	      //fprintf(stdout,"l=%d whichzone: %d m=%d bz k: %d->%d j: %d->%d i: %d->%d\n",l,whichzone,m,kk,k,jj,j,ii,i); fflush(stdout);
	    }// done with real assignment of what to copy for each boundary zone
	  }// done with: else if normal boundary stuff
	}// done with else if got good computational zone	
	whichzone++; // next zone
      } // end if boundary zone
    }// end over all zones
  }// end over index types




  // now figure out flags
  for(l=6;l<=7;l++){
    LOOPBOUND(l){
      itemp=bzs[l][temptempi][0];
      for(m=0;m<itemp;m++){
	
	ii=bzs[l][temptempi][m*3+0+1];
	jj=bzs[l][temptempi][m*3+1+1];
	kk=bzs[l][temptempi][m*3+2+1];

	
      } // end over zones
      // now set
      bzs[l][temptempi][-1]=0;
      bzs[l][temptempi][-2]=0;
      bzs[l][temptempi][-3]=mask[l][k][j][i]; // holds boundary type
    }// end main loop
  }
  // now figure out flags
  for(l=8;l<=9;l++){
    LOOPBOUND(l){
      itemp=bzs[l][temptempi][0];
      ininflag=0;
      outinflag=0;
      for(m=0;m<itemp;m++){
	
	ii=bzs[l][temptempi][m*3+0+1];
	jj=bzs[l][temptempi][m*3+1+1];
	kk=bzs[l][temptempi][m*3+2+1];

	// assumes that 1+ boundary zone(along that line) is enough to deal with inflow problem (otherwise no generally valid way to deal with other zones)
	// need to worry about 4 and 94?
	if(INFLOWCHECKIX1&&(ii-i>0)&&(jj-j==0)&&(kk-k==0)){ ininflag=1; } // so must force inflow condition
	if(INFLOWCHECKOX1&&(ii-i<0)&&(jj-j==0)&&(kk-k==0)){ outinflag=1; } // so must force inflow condition      
      } // end over zones
      // now set
      if(mask[l][k][j][i]==4){
	bzs[l][temptempi][-1]=ininflag;
	bzs[l][temptempi][-2]=outinflag;
      }
      else{
	bzs[l][temptempi][-1]=0;
	bzs[l][temptempi][-2]=0;
      }
      bzs[l][temptempi][-3]=mask[l][k][j][i]; // holds boundary type
    }// end main loop
  }
  // GODMARK(commented)
  //  fflush(stdout);
  
  for(l=10;l<=11;l++){
    LOOPBOUND(l){
      ininflag=0;
      outinflag=0;
      itemp=bzs[l][temptempi][0];
      for(m=0;m<itemp;m++){
	
	ii=bzs[l][temptempi][m*3+0+1];
	jj=bzs[l][temptempi][m*3+1+1];
	kk=bzs[l][temptempi][m*3+2+1];

	// assumes that 1+ boundary zone(along that line) is enough to deal with inflow problem (otherwise no generally valid way to deal with other zones)
	if(INFLOWCHECKIX2&&(ii-i==0)&&(jj-j>0)&&(kk-k==0)){ ininflag=1; } // so must force inflow condition
	if(INFLOWCHECKOX2&&(ii-i==0)&&(jj-j<0)&&(kk-k==0)){ outinflag=1; } // so must force inflow condition
	
      } // end over zones
      // now set
      if(mask[l][k][j][i]==4){
	bzs[l][temptempi][-1]=ininflag;
	bzs[l][temptempi][-2]=outinflag;
      }
      else{
	bzs[l][temptempi][-1]=0;
	bzs[l][temptempi][-2]=0;
      }
      bzs[l][temptempi][-3]=mask[l][k][j][i]; // holds boundary type
    }// end main loop
  }
  for(l=12;l<=13;l++){
    LOOPBOUND(l){
      ininflag=0;
      outinflag=0;
      itemp=bzs[l][temptempi][0];
      for(m=0;m<itemp;m++){
	
	ii=bzs[l][temptempi][m*3+0+1];
	jj=bzs[l][temptempi][m*3+1+1];
	kk=bzs[l][temptempi][m*3+2+1];
	
	// assumes that 1+ boundary zone(along that line) is enough to deal with inflow problem (otherwise no generally valid way to deal with other zones)
	if(INFLOWCHECKIX3&&(ii-i==0)&&(jj-j==0)&&(kk-k>0)){ ininflag=1; } // so must force inflow condition
	if(INFLOWCHECKOX3&&(ii-i==0)&&(jj-j==0)&&(kk-k<0)){ outinflag=1; } // so must force inflow condition
      } // end over zones
      // now set
      if(mask[l][k][j][i]==4){
	bzs[l][temptempi][-1]=ininflag;
	bzs[l][temptempi][-2]=outinflag;
      }
      else{
	bzs[l][temptempi][-1]=0;
	bzs[l][temptempi][-2]=0;
      }
      bzs[l][temptempi][-3]=mask[l][k][j][i]; // holds boundary type
    }// end main loop
  }

  // b-field bctype set
  for(l=14;l<=19;l++){
    LOOPBOUND(l){
      itemp=bzs[l][temptempi][0];
      for(m=0;m<itemp;m++){
	
	ii=bzs[l][temptempi][m*3+0+1];
	jj=bzs[l][temptempi][m*3+1+1];
	kk=bzs[l][temptempi][m*3+2+1];
	
      } // end over zones
      // now set
      bzs[l][temptempi][-1]=0;
      bzs[l][temptempi][-2]=0;
      bzs[l][temptempi][-3]=mask[l][k][j][i]; // holds boundary type
    }// end main loop
  }

  // GODMARK(commented)
  // LOOPF{
  //  fprintf(stdout,"%d %d %d %d\n",k,j,i,newmask[k][j][i]);
  // }

  // the below checks whether 20-22 is rim-to-rim as referenced off LOOPC, which is zone-centered comp zones.
  // this should only show problems with inner bzones when using LOOPTYPE==3 if end up not using indx but iindx, but normally always uses indx for correct diagnostics
  LOOPDIVB{
    if(
       (mask[20][k][j][i]!=0)||
       (mask[21][k][j][i]!=0)||
       (mask[22][k][j][i]!=0)||
       (mask[20][k][j][i+1]!=0)||
       (mask[21][k][j+1][i]!=0)||
       (mask[22][k+1][j][i]!=0)
       ){
      fprintf(stdout,"problem at: k=%d j=%d i=%d\n",k,j,i);
      fprintf(stdout," %d %d %d %d %d %d\n",mask[20][k][j][i],mask[21][k][j][i],mask[22][k][j][i],mask[20][k][j][i+1],mask[21][k][j+1][i],mask[22][k+1][j][i]);

      fflush(stdout);
    }
  }
  // now figure out flags for flux (currently only uses 1 flag)
  // flux need not know boundary type or other inflow/reflect flags
  // need to determine if this zone is on the + or - face of each direction, this is probably not the most general method
  // most general method should use ii-i, etc. for each positioned variable

  // no need for global flags currently
  for(l=23;l<=28;l++){
    countfaces[0]=0;
    countfaces[1]=0;
    countfaces[2]=0;
    
    LOOPF{
      debugmask[k][j][i]=0;
    }

    LOOPBOUND(l){
      // see if left or right, etc.
      // inners
      if(l==23){
	// this flux surface needs to know whether an MPI zone is really a boundary/null zone or comp zone so surface is smooth/correct across CPUs
	if(mask[l][k][j][i+1]>=90) maskfix[0]=mask[l][k][j][i+1]-90;
	else maskfix[0]=mask[l][k][j][i+1];

	if(mask[l][k][j][i-1]>=90) maskfix[1]=mask[l][k][j][i-1]-90;
	else maskfix[1]=mask[l][k][j][i-1];

	// check if flux goes is on surface or not
	if(!( ((maskfix[0]==0)&&(maskfix[1]!=0))||((maskfix[0]!=0)&&(maskfix[1]==0)) ) ) bzs[l][temptempi][0]=0;
	else{ // then really a flux on surface
	  if(x[1][1][i]<0.0) bzs[l][temptempi][0]=1;
	  else bzs[l][temptempi][0]=-1;
	}
      }
      if(l==25){
	// this flux surface needs to know whether an MPI zone is really a boundary/null zone or comp zone so surface is smooth/correct across CPUs
	if(mask[l][k][j+1][i]>=90) maskfix[0]=mask[l][k][j+1][i]-90;
	else maskfix[0]=mask[l][k][j+1][i];

	if(mask[l][k][j-1][i]>=90) maskfix[1]=mask[l][k][j-1][i]-90;
	else maskfix[1]=mask[l][k][j-1][i];

	// check if flux goes is on surface or not
	if(!( ((maskfix[0]==0)&&(maskfix[1]!=0))||((maskfix[0]!=0)&&(maskfix[1]==0)) ) ) bzs[l][temptempi][0]=0;
	else{ // then really a flux on surface
	  if(x[1][2][j]<0.0) bzs[l][temptempi][0]=1;
	  else bzs[l][temptempi][0]=-1;
	}
      }
      if(l==27){
	// this flux surface needs to know whether an MPI zone is really a boundary/null zone or comp zone so surface is smooth/correct across CPUs
	if(mask[l][k+1][j][i]>=90) maskfix[0]=mask[l][k+1][j][i]-90;
	else maskfix[0]=mask[l][k+1][j][i];

	if(mask[l][k-1][j][i]>=90) maskfix[1]=mask[l][k-1][j][i]-90;
	else maskfix[1]=mask[l][k-1][j][i];

	// check if flux goes is on surface or not
	if(!( ((maskfix[0]==0)&&(maskfix[1]!=0))||((maskfix[0]!=0)&&(maskfix[1]==0)) ) ) bzs[l][temptempi][0]=0;
	else{ // then really a flux on surface
	  if(x[1][3][k]<0.0) bzs[l][temptempi][0]=1;
	  else bzs[l][temptempi][0]=-1;
	}
      }
      // outers
      if(l==24){
	// this flux surface needs to know whether an MPI zone is really a boundary/null zone or comp zone so surface is smooth/correct across CPUs
	if(mask[l][k][j][i+1]>=90) maskfix[0]=mask[l][k][j][i+1]-90;
	else maskfix[0]=mask[l][k][j][i+1];

	if(mask[l][k][j][i-1]>=90) maskfix[1]=mask[l][k][j][i-1]-90;
	else maskfix[1]=mask[l][k][j][i-1];

	// check if flux goes is on surface or not
	if(!( ((maskfix[0]==0)&&(maskfix[1]!=0))||((maskfix[0]!=0)&&(maskfix[1]==0)) ) ) bzs[l][temptempi][0]=0;
	else{ // then really a flux on surface
	  if(x[1][1][i]<0.0) bzs[l][temptempi][0]=-1;
	  else bzs[l][temptempi][0]=1;
	}
      }
      if(l==26){
	// this flux surface needs to know whether an MPI zone is really a boundary/null zone or comp zone so surface is smooth/correct across CPUs
	if(mask[l][k][j+1][i]>=90) maskfix[0]=mask[l][k][j+1][i]-90;
	else maskfix[0]=mask[l][k][j+1][i];

	if(mask[l][k][j-1][i]>=90) maskfix[1]=mask[l][k][j-1][i]-90;
	else maskfix[1]=mask[l][k][j-1][i];

	// check if flux goes is on surface or not
	if(!( ((maskfix[0]==0)&&(maskfix[1]!=0))||((maskfix[0]!=0)&&(maskfix[1]==0)) ) ) bzs[l][temptempi][0]=0;
	else{ // then really a flux on surface
	  if(x[1][2][j]<0.0) bzs[l][temptempi][0]=-1;
	  else bzs[l][temptempi][0]=1;
	}
      }
      if(l==28){
	// this flux surface needs to know whether an MPI zone is really a boundary/null zone or comp zone so surface is smooth/correct across CPUs
	if(mask[l][k+1][j][i]>=90) maskfix[0]=mask[l][k+1][j][i]-90;
	else maskfix[0]=mask[l][k+1][j][i];

	if(mask[l][k-1][j][i]>=90) maskfix[1]=mask[l][k-1][j][i]-90;
	else maskfix[1]=mask[l][k-1][j][i];

	// check if flux goes is on surface or not
	if(!( ((maskfix[0]==0)&&(maskfix[1]!=0))||((maskfix[0]!=0)&&(maskfix[1]==0)) ) ) bzs[l][temptempi][0]=0;
	else{ // then really a flux on surface
	  if(x[1][3][k]<0.0) bzs[l][temptempi][0]=-1;
	  else bzs[l][temptempi][0]=1;
	}
      }
      //      fprintf(stdout,"%2d %2d %2d %2d %2d\n",l,k,j,i,bzs[l][temptempi][0]);
      debugmask[k][j][i]=bzs[l][temptempi][0];
      if(bzs[l][temptempi][0]==1) countfaces[0]++;
      else if(bzs[l][temptempi][0]==-1) countfaces[1]++;
      else if(bzs[l][temptempi][0]==0) countfaces[2]++;
    }// end main loop
    fprintf(log_file,"numfaces l=%d 1: %d -1: %d 0: %d total!=0: %d\n",l,countfaces[0],countfaces[1],countfaces[2],countfaces[0]+countfaces[1]); fflush(log_file);

    // do symmetry checks on -1,1
    if(numprocs==1){
      LOOPF{
	if((l==23)||(l==24)){
	  if((int)debugmask[k][j][i]!=-(int)debugmask[k][j][N1-i]){
	    fprintf(stdout,"bad l=%d: %d %d %d: %d  %d %d %d: %d :: -1m: %d 0m: %d 1m: %d\n",l,k,j,i,(int)debugmask[k][j][i],k,j,N1-i,(int)debugmask[k][j][N1-i],mask[l][k][j][i-1],mask[l][k][j][i],mask[l][k][j][i+1]);
	  }
	}
	if((l==25)||(l==26)){
	  if((int)debugmask[k][j][i]!=-(int)debugmask[k][N2-j][i]){
	    fprintf(stdout,"bad l=%d: %d %d %d: %d  %d %d %d: %d :: -1m: %d 0m: %d 1m: %d\n",l,k,j,i,(int)debugmask[k][j][i],k,N2-j,i,(int)debugmask[k][N2-j][i],mask[l][k][j-1][i],mask[l][k][j][i],mask[l][k][j+1][i]);
	  }
	}
	if((l==27)||(l==28)){
	  if((int)debugmask[k][j][i]!=-(int)debugmask[N3-k][j][i]){
	    fprintf(stdout,"bad l=%d: %d %d %d: %d  %d %d %d: %d :: -1m: %d 0m: %d 1m: %d\n",l,k,j,i,(int)debugmask[k][j][i],N3-k,j,i,(int)debugmask[N3-k][j][i],mask[l][k-1][j][i],mask[l][k][j][i],mask[l][k+1][j][i]);
	  }
	}
      }    
    }
    /*
    LOOPF{
      fprintf(stdout,"%2d %2d %2d %2d %2d\n",l,k,j,i,(int)debugmask[k][j][i]);
    }
    */
  }
  fflush(stdout);
  //  fflush(stdout); myexit(0);
  


  


  //  myexit(0);


  // GODMARK GODMARK(commented)
  //  numiter[0]=1; // test
  //numiter[6]=1; // test
  //numiter[7]=1; // test
  //numiter[8]=1; // test
  //  fflush(stdout);
  // GODMARK(commented)
  //myexit(0);


  // do symmetry check on masks
  
  // GODMARK -- commented
  /*
  // symmetry check
  LOOPF{
    // mask[0]
    //s[1][k][j][i]=mask[0][k][j][i]; itemp=0; // good
    // mask[6]
    //s[1][k][j][i]=mask[6][k][j][i]; itemp=0; // good
    // mask[7]
    //s[1][k][j][i]=mask[7][k][j][i]; itemp=0; // good
    // mask[1]
    //    s[1][k][j][i]=mask[1][k][j][i]; itemp=0; // good
    // mask[2]
    //    s[1][k][j][i]=mask[2][k][j][i]; itemp=0; // good
    // mask[5]
    //s[1][k][j][i]=mask[5][k][j][i]; itemp=0; // good

    // mask[8]
    //s[1][k][j][i]=mask[8][k][j][i]; itemp=1;
    // mask[9]
    //s[1][k][j][i]=mask[9][k][j][i]; itemp=1;


    // mask[12]
    //s[1][k][j][i]=mask[12][k][j][i]; itemp=3;
    // mask[13]
    s[1][k][j][i]=mask[13][k][j][i]; itemp=3;
  }
  
  //LOOPF{
  //  fprintf(stdout,"%2d ",mask[9][k][j][i]);
  //  if(i==N1+N1BND-1) fprintf(stdout,"\n");
  // }
  //fflush(stdout);
  
  //  dump(NULL,666,DTYPE,0);
  fprintf(log_file,"maskcheck\n"); fflush(stdout);
  symmetry_check(-itemp-1);
  //  myexit(0);
  */
  
  // GODMARK -- commented
  /*
  LOOPF{
    for(l=0;l<NUMINDEX;l++){
      fprintf(log_file,"%d ",mask[l][k][j][i]);
    }
    fprintf(log_file,"\n");
  }
  fflush(log_file);
  itemp+=init_outgparm(-1); // -1 means out all within define params
  myexit(0);
  */

  /*
  // see if bound works
  LOOPF{
    for(l=1;l<=NUMSCA;l++){
      s[l][k][j][i]=1;
    }
    for(l=1;l<=NUMVEC;l++){
      for(m=1;m<=3;m++){
	v[l][m][k][j][i]=1;
      }
    }
  }

  // now dump bad stuff into bzones and see if bound() corrects them all
  l=6;
  LOOPF{
    if(mask[l][k][j][i]==BCTYPE){
      s[1][k][j][i]=10;
      s[2][k][j][i]=10;
      s[3][k][j][i]=10;
    }
    else if(mask[l][k][j][i]==-1){
      s[1][k][j][i]=5;
      s[2][k][j][i]=5;
      s[3][k][j][i]=5;
    }
  }
  l=7;
  LOOPF{
    if(mask[l][k][j][i]==BCTYPE){
      s[1][k][j][i]=10;
      s[2][k][j][i]=10;
      s[3][k][j][i]=10;
    }
    else if(mask[l][k][j][i]==-1){
      s[1][k][j][i]=5;
      s[2][k][j][i]=5;
      s[3][k][j][i]=5;
    }
  }

  l=8;
  LOOPF{
    if(mask[l][k][j][i]==BCTYPE){
      v[1][1][k][j][i]=10;
    }
    else if(mask[l][k][j][i]==-1){
      v[1][1][k][j][i]=5;
    }
  }
  l=9;
  LOOPF{
    if(mask[l][k][j][i]==BCTYPE){
      v[1][1][k][j][i]=10;
    }
    else if(mask[l][k][j][i]==-1){
      v[1][1][k][j][i]=5;
    }
  }

  l=10;
  LOOPF{
    if(mask[l][k][j][i]==BCTYPE){
      v[1][2][k][j][i]=10;
    }
    else if(mask[l][k][j][i]==-1){
      v[1][2][k][j][i]=5;
    }
  }
  l=11;
  LOOPF{
    if(mask[l][k][j][i]==BCTYPE){
      v[1][2][k][j][i]=10;
    }
    else if(mask[l][k][j][i]==-1){
      v[1][2][k][j][i]=5;
    }
  }

  l=12;
  LOOPF{
    if(mask[l][k][j][i]==BCTYPE){
      v[1][3][k][j][i]=10;
    }
    else if(mask[l][k][j][i]==-1){
      v[1][3][k][j][i]=5;
    }
  }
  l=13;
  LOOPF{
    if(mask[l][k][j][i]==BCTYPE){
      v[1][3][k][j][i]=10;
    }
    else if(mask[l][k][j][i]==-1){
      v[1][3][k][j][i]=5;
    }
  }

  l=14;
  LOOPF{
    if(mask[l][k][j][i]==BCTYPE){
      v[2][1][k][j][i]=10;
    }
    else if(mask[l][k][j][i]==-1){
      v[2][1][k][j][i]=5;
    }
  }
  l=15;
  LOOPF{
    if(mask[l][k][j][i]==BCTYPE){
      v[2][1][k][j][i]=10;
    }
    else if(mask[l][k][j][i]==-1){
      v[2][1][k][j][i]=5;
    }
  }

  l=16;
  LOOPF{
    if(mask[l][k][j][i]==BCTYPE){
      v[2][2][k][j][i]=10;
    }
    else if(mask[l][k][j][i]==-1){
      v[2][2][k][j][i]=5;
    }
  }
  l=17;
  LOOPF{
    if(mask[l][k][j][i]==BCTYPE){
      v[2][2][k][j][i]=10;
    }
    else if(mask[l][k][j][i]==-1){
      v[2][2][k][j][i]=5;
    }
  }

  l=18;
  LOOPF{
    if(mask[l][k][j][i]==BCTYPE){
      v[2][3][k][j][i]=10;
    }
    else if(mask[l][k][j][i]==-1){
      v[2][3][k][j][i]=5;
    }
  }
  l=19;
  LOOPF{
    if(mask[l][k][j][i]==BCTYPE){
      v[2][3][k][j][i]=10;
    }
    else if(mask[l][k][j][i]==-1){
      v[2][3][k][j][i]=5;
    }
  }
  DYNAMICMM=1;
  image(998,-1,-1,0,0);
  dump(NULL,998,DTYPE,0);
#if(!PPCLEAN)
  bound(NULL,NULL,-1,-1,123);
#endif
  image(999,-1,-1,0,0);
  dump(NULL,999,DTYPE,0);
  myexit(0);
  // end to see if bound works
  */




  // free mask memory
  for(i=0;i<NUMINDEX;i++){
    // reset position
    mask[i]=(short (*) [N2M][N1M])(&(mask[i][-N3BND][-N2BND][-N1BND]));    
    free(mask[i]);
  }
  free(mask);



  fprintf(log_file,"end: init_bc_gen1\n"); fflush(log_file);
  return(0);
}



int init(int argc,
	 char *argv[],
	 char*envp[])
{
  int k,j,i,l,m;
  int seed;
  FTYPE beta,nj;
  int error=0;
  char temps[100];
  int itemp;
  char extension[MAXFILENAME];


  fprintf(log_file,"begin: init ... "); fflush(log_file);

  // left over parameters
  seed=13;
  ranc(seed) ;   /* fire up random number generator */
  beta=0.1;
  nj=1.0;
  npoldschool=NPVER; // start at current version
  avg2doldschool=AVG2DVER; // start at current version
  pdump_start=0;
  dump_start=npdump_start=0;
  adump_start=0;
  floor_start=0;
  image_start=0;
  ireenter=0;


  init_checks();

  error+=init_mem();
  error+=init_pointers();

  // initialize general stuff
  init_general(); // should be before init_paramstot()

  error+=init_paramstot(seed,beta,nj); // must come before most everything else


  init_placeongrid();


  // get or set parameters needed for data(get/set)
  if( (runtype==3)||(runtype==2)||(runtype==22)){
    error+=init_runpar(); //function for par files
  }

  error+=init_paramspresets(seed,beta,nj); // must come before grid assignment/etc, must come after parameter set/read
  error+=init_otherparams(); // other parameter init



  // Setup directory structure in case do not exist
  // (must come after parampresets())
  if(myid<=0){
    if(USEGM==0){
      sprintf(temps,"mkdir %s%s",DATADIR,IMAGEDIR);
      system(temps);
      if(deleteolddat){
	sprintf(temps,"rm %s*",DATADIR);
	strcat(temps,DATEXT);
	strcat(temps,"*");
	system(temps);
      }
      if(deleteoldpar){
	sprintf(temps,"rm %s*",DATADIR);
	strcat(temps,PAREXT);
	strcat(temps,"*");
	system(temps);
      }
    }
  }
#if(USEMPI)
  // pause so can create directory structure before going on
  MPI_Barrier(MPI_COMM_WORLD);
#endif



#if(TVDLF==0)
  error+=init_mainbc(periodicx1,skipix1,reflectix1,reflectox1,periodicx2,skipix2,reflectix2,reflectox2,periodicx3,skipix3,reflectix3,reflectox3);
#endif

  
  if((runtype==1)||(runtype==0)||(runtype==22)){
    error+=init_dx(N1,N2,N3,N1BND,N2BND,N3BND,startx,0,0);
    error+=init_x(N1,N2,N3,N1BND,N2BND,N3BND,startx,0,0);
  }

  // don't store these to file, so must always compute
  if(LOWMEMMODE==0) error+=init_diffs();
  error+=init_reduceddimension();


#if((TVDLF==0)&&(POSTPROC==0))
  if(LOWMEMMODE==0) init_compsave();
#endif

#if((TVDLF==0)&&(POSTPROC==0)) // no need if POSTPROC==1, don't want if PPCLEAN=1
 // requires grid to be setup
  error+=init_bc(simplebc,bcix1,bcox1,bcix2,bcox2,bcix3,bcox3);
#endif

#if(TVDLF==1)
  init_tvdlfgrid();
#endif

  // get data
  if( (runtype==3)||(runtype==1)){
    if( (directinput>0)&&(directinput!=3) ){
      error+=init_reentrance(); // gets # of dumps and sets dump_start or pdump_start
    }
    else if(directinput==3){ // then don't care about number of dumps
      init_reentrance2(timereenter,1); // base start dump on time given and DTs set
    }
    if(DUMPSM==0){ // only valid to read in dump files if not in SM interpolated format(in general)
      if(directinput>0) itemp=dump_start; else itemp=-1;
      error+=init_rundat(itemp,DTYPE);//function for dump dat files
    }
    else{
      if(directinput>0) itemp=pdump_start; else itemp=-1;
      error+=init_rundat(itemp,PDTYPE);//function for dump dat files
    }
    if(directinput>0){
      init_reentrance2(t,1); // make sure actual time in file consistent with guess or directinput==3 timereenter
    }
    timereenter=t; // this is the real time that the run was restarted
  }
  if(myid<=0){
    fprintf(stderr,"real starts: %d %d %d %d %d %d ir: %d\n",pdump_start,dump_start,npdump_start,adump_start,floor_start,image_start,ireenter);
    fflush(stderr);
  }
    // set data
  if((runtype==0)||(runtype==2)||(runtype==22)){
    error+=init_data();
    // GODMARK
    // only init below 2 if not reentrant on data
    init_inflows();
    init_radiations();
  }
  else{
    // need to solve analytic solution for many cases: need bc data, etc...just as with nonloading case.
    analsolve(POSTPROC);
#if(RAD)
    initrad();
#endif

    tdep_compute(); // compute time dep stuff
    // place any overrides here.  e.g. file loads potential, but potential is not bounded since static and needs proper static dependence on boundary zones, so force s[3] to be sanal[3], ignoring file data, for wherever you assign the override.  Can use init_data as template.
    
    LOOPF{
      s[3][k][j][i] = sanal[3][k][j][i] ;
    }
  }

    
  /* enforce boundary conditions on all scalars and vectors */
  if((TVDLF==0)&&(POSTPROC==0)){
    //dump(NULL,998,DTYPE,0); // GODMARK 
    //image(998,-1,-1,0,0); // GODMARK (to show variables before bounding)
    fprintf(log_file,"before bound symmetry check\n"); fflush(log_file);
    symmetry_check(0); // before bound symmetry check// GODMARK (use -1 for image output of badness)
    fprintf(log_file,"before bound divb=0 check\n"); fflush(log_file);
    divb0check(0); // check before bound
    //image(666,-1,-1,0,0); // GODMARK  (to show asymmetry)
    //boundtest(0);
    //boundtest(1);
    //boundtest(2);
    fprintf(log_file,"at bound\n"); fflush(log_file);
#if(!PPCLEAN)
    bound(NULL,NULL,-1,-1,123);
#endif
    fprintf(log_file,"after bound divb=0 check\n"); fflush(log_file);
    divb0check(0); // check after bound
    //dump(NULL,999,DTYPE,0); // GODMARK
    //image(999,-1,-1,0,0); // GODMARK (to show variables after bounding)
    //myexit(0);
  }

  if(POSTPROC==0){ // otherwise do manually in postproc.c
    if(!((PRODUCTION==1)&&(N1>1)&&(N2>1)&&(N3>1))){
      error+=init_outgparm(-1); // -1 means out all within define params
    }
  }


  // init final things
  accountstoreset(); // must be before any accounting
  if((POSTPROC==0)&&(TVDLF==0)){
    init_floor();
    if(VISCMEM&&(visc_real==1)){
      init_visc();
    }
    if(RESMEM&&(res_real==1)){
      init_res();
    }  
    init_loss();
  }

  // more checks(doesn't check during runtime, so memory use may be invalid if variables turn on and mem off)
  if((MDOTMEM==0)&&(mdotin==1)){
    fprintf(fail_file,"Must turn on MDOTMEM==1 to use mdotin=1\n");
    myexit(1);
  }
  if((VISCMEM==0)&&(visc_real==1)){
    fprintf(fail_file,"Must turn on VISCMEM==1 to use visc_real=1\n");
    myexit(1);
  }
  if((RESMEM==0)&&(res_real==1)){
    fprintf(fail_file,"Must turn on RESMEM==1 to use res_real=1\n");
    myexit(1);
  }

  
  if(appendold==1){
    sprintf(WRITETYPE,"at+");
  }
  else{
    sprintf(WRITETYPE,"wt");
  }
  strcpy(extension,OUTEXT);
  
  if(POSTPROC==0){
    // setup computational log files
    if(DODTDIAG){
      sprintf(temps,"%s0_logdt%s%s",DATADIR,extension,myidtxt) ;
      
      if((logdt_file=fopen(temps,WRITETYPE))==NULL){ // naive append if appendold==1
	fprintf(stderr,"dtdiag: Cannot open: %s\n",temps);
	myexit(1);
      }
    }
    if(DOLOGSTEP){
      if(myid<=0){
	sprintf(temps,"%s0_logstep%s",DATADIR,extension) ;
	
	if((logstep_file=fopen(temps,WRITETYPE))==NULL){ // naive append if appendold==1
	  fprintf(stderr,"logstep: Cannot open: %s\n",temps);
	  myexit(1);
	}
      }
    }
    if(DOLOGPERF){
      if(myid<=0){
	sprintf(temps,"%s0_logperf%s",DATADIR,extension) ;
	
	if((logperf_file=fopen(temps,WRITETYPE))==NULL){ // naive append if appendold==1
	  fprintf(stderr,"logperf: Cannot open: %s\n",temps);
	  myexit(1);
	}
      }
    }
    
    if(DOFLOORDIAG>0){
      
      sprintf(temps,"%s0_logfl%s%s",DATADIR,extension,myidtxt) ;
      
      if((logfl_file=fopen(temps,WRITETYPE))==NULL){ // naive append if appendold==1
	fprintf(stderr,"floordiag: Cannot open: %s\n",temps);
	myexit(1);
      }
    }
  }
  if(myid<=0){
    fprintf(stderr,"#gocont: %d runtype: %d directinput: %d\n",gocont,runtype,directinput);
    fprintf(logfull_file,"#gocont: %d runtype: %d directinput: %d\n",gocont,runtype,directinput);
  }
  fprintf(log_file,"end: init\n"); fflush(log_file);
  //LOOPF{
  //	  fprintf(stderr,"%d %d %d %15.10g\n",k,j,i,v[2][1][k][j][i]); fflush(stderr);
  //	}
    //	myexit(0);



  return(error);
}// end init()


// very expensive for 1-D problems to do each time step
void init_loss_rect(void)
{
  int i,j,k,l,m,ll;
  static int firsttime=1;
  static int NS1[3+1],NS2[3+1];


  if(firsttime==1){
    fprintf(log_file,"begin: init_loss_rect ... "); fflush(log_file);
    NS1[1]=N2; NS2[1]=N3;
    NS1[2]=N1; NS2[2]=N3;
    NS1[3]=N1; NS2[3]=N2;
    if(N1==1){ NS1[1]=1; NS2[1]=1; } // blow off since likely don't care
    if(N2==1){ NS1[2]=1; NS2[2]=1; } // blow off since likely don't care
    if(N3==1){ NS1[3]=1; NS2[3]=1; } // blow off since likely don't care
  }

  if(DOLOSSDIAG){
    // initialize loss data
    for(i=1;i<=NUMSCA;i++){
      for(j=1;j<=3;j++){
	for(k=0;k<=1;k++){
	  for(l=0;l<NS2[j];l++){ // only need to go over needed zones
	    for(ll=0;ll<NS1[j];ll++){
	      losss[i][j][k][l][ll]=0.0;
	    }
	  }
	}
      }
    }
    for(i=1;i<=NUMVEC;i++){
      for(m=0;m<=3;m++){  // 0 for KE
	for(j=1;j<=3;j++){
	  for(k=0;k<=1;k++){
	    for(l=0;l<NS2[j];l++){
	      for(ll=0;ll<NS1[j];ll++){
		lossv[i][m][j][k][l][ll]=0.0;
	      }
	    }
	  }
	}
      }
    }
    for(i=1;i<=1;i++){
      for(j=1;j<=3;j++){
	for(k=0;k<=1;k++){
	  for(l=0;l<NS2[j];l++){
	    for(ll=0;ll<NS1[j];ll++){
	      lossvisc[i][j][k][l][ll]=0.0;
	    }
	  }
	}
      }
    }
  }
  if(firsttime==1){
    fprintf(log_file,"end: init_loss_rect\n"); fflush(log_file);
  }

  firsttime=0;
}


// very expensive for 1-D problems to do each time step
void init_loss_gen(void)
{
  int i,j,k,l,m,ll;
  static int firsttime=1;
  static int NS1[3+1],NS2[3+1];


  if(firsttime==1){
    fprintf(log_file,"begin: init_loss_gen ... "); fflush(log_file);
    NS1[1]=N2; NS2[1]=N3;
    NS1[2]=N1; NS2[2]=N3;
    NS1[3]=N1; NS2[3]=N2;
    if(N1==1){ NS1[1]=1; NS2[1]=1; } // blow off since likely don't care
    if(N2==1){ NS1[2]=1; NS2[2]=1; } // blow off since likely don't care
    if(N3==1){ NS1[3]=1; NS2[3]=1; } // blow off since likely don't care
  }
  
  if(DOLOSSDIAG){
    // initialize loss data
    for(i=1;i<=NUMLOSSVAR;i++){
      for(k=0;k<=1;k++){
	lossflux[i][k]=0.0;
      }
    }
  }
  if(firsttime==1){
    fprintf(log_file,"end: init_loss_gen\n"); fflush(log_file);
  }

  firsttime=0;
}

void init_general(void)// as is setup for a tori case.  Needed to avoid general floating point issues in image routines if not assigned
{
  int i,j,k,l;

  fprintf(log_file,"begin: init_general ... "); fflush(log_file);

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
  mmv[i][j][2][0][0]=-1.0;
  mmv[i][j][2][0][1]=1.0;
  // Bx1
  mmv[i][j][2][1][0]=-1.0;
  mmv[i][j][2][1][1]=1.0;
  // Bx2
  mmv[i][j][2][2][0]=-1.0;
  mmv[i][j][2][2][1]=1.0;
  // Bx3
  mmv[i][j][2][3][0]=-1.0;
  mmv[i][j][2][3][1]=1.0;

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
  mmv[i][j][2][0][0]=-1.0;
  mmv[i][j][2][0][1]=1.0;
  // rho*Bx1
  mmv[i][j][2][1][0]=-1.0;
  mmv[i][j][2][1][1]=1.0;
  // rho*Bx2
  mmv[i][j][2][2][0]=-1.0;
  mmv[i][j][2][2][1]=1.0;
  // rho*Bx3
  mmv[i][j][2][3][0]=-1.0;
  mmv[i][j][2][3][1]=1.0;



  
  i=1; // zoom view (vectors)
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
  mmv[i][j][2][0][0]=-1.0;
  mmv[i][j][2][0][1]=1.0;
  // Bx1
  mmv[i][j][2][1][0]=-1.0;
  mmv[i][j][2][1][1]=1.0;
  // Bx2
  mmv[i][j][2][2][0]=-1.0;
  mmv[i][j][2][2][1]=1.0;
  // Bx3
  mmv[i][j][2][3][0]=-1.0;
  mmv[i][j][2][3][1]=1.0;

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
  mmv[i][j][2][0][0]=-1.0;
  mmv[i][j][2][0][1]=1.0;
  // rho*Bx1
  mmv[i][j][2][1][0]=-1.0;
  mmv[i][j][2][1][1]=1.0;
  // rho*Bx2
  mmv[i][j][2][2][0]=-1.0;
  mmv[i][j][2][2][1]=1.0;
  // rho*Bx3
  mmv[i][j][2][3][0]=-1.0;
  mmv[i][j][2][3][1]=1.0;




  // define outer region when interpolation is used.
  // same order as scalar/vector arrays

  // for images
  for(i=0;i<ITYPES;i++){ // both views
    for(j=0;j<CTYPES;j++){ // both comps
      outerdefs[i][j][1]=mms[i][j][1][0]; // rho
      outerdefs[i][j][2]=mms[i][j][2][0]; // en
      outerdefs[i][j][3]=mms[i][j][3][0]; // pot
      
      outerdefv[i][j][1][0]=0; // magnitude of v
      outerdefv[i][j][1][1]=0; // v1
      outerdefv[i][j][1][2]=0; // v2
      outerdefv[i][j][1][3]=0; // v3
      
      outerdefv[i][j][2][0]=0;
      outerdefv[i][j][2][1]=0;
      outerdefv[i][j][2][2]=0;
      outerdefv[i][j][2][3]=0;
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

  if(GRAVACC&&(!POSTPROC)){
    for(l=1;l<=3;l++){
      LOOPF{
	gravacc[l][k][j][i]=0;
      }
    }
  }

  if(GRAVITOMAGNETIC&&(!POSTPROC)){
    LOOPF{
      gravitom[1][k][j][i]=0;
      gravitom[2][k][j][i]=0;
    }
  }


  fprintf(log_file,"end: init_general\n"); fflush(log_file);
}


#if(USEMPI)

int init_MPI(int argc, char *argv[])
{
  SFTYPE numcpux1low,numcpux1high,holdrealncpux1;
  SFTYPE tempx1,tempx2,nz1,nz2;
  SFTYPE n1,n2,n3,nc,A,B,C;

  fprintf(stderr,"begin: init_MPI\n"); fflush(stderr);

  MPI_Init(&argc,&argv);
  myargs(argc,argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  sprintf(myidtxt,CPUTXT,myid);
  MPI_Get_processor_name(processor_name,&procnamelen);
  MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);

  if(MAXCPUS<numprocs){
    fprintf(stderr,"Must increase MAXCPUS in global.h, %d is too many\n",numprocs);
    myexit(1);
  }

  // determine how grid should be partitioned -- optimal setting do or report
  init_optimalmpi(cpugeompick);

  if(myid<=0){
    fprintf(stderr,"numprocs=%d ncpux1=%d ncpux2=%d ncpux3=%d\n",numprocs,ncpux1,ncpux2,ncpux3);
  }

  init_MPIgroup();

  fprintf(stderr,"proc: %s on %s\n", myidtxt, processor_name);

  fprintf(stderr,"end: init_MPI\n"); fflush(stderr);
  return(0);
}

void init_optimalmpi(int which)
{           
  int i,lowesti;
  SFTYPE numcpux1low,numcpux1high,fncpux1,fncpux2,fncpux3,testncpux2,holdrealncpux1,lowestnz;
  SFTYPE tempx1[4+1],tempx2[4+1],nz[4+1];
  SFTYPE n1,n2,n3,nc,A,B,C;
  SFTYPE shouldn1,shouldn2,shouldn3,shouldncpux1,shouldncpux2,shouldncpux3;
  // determine how grid should be partitioned -- optimal setting

  fprintf(stderr,"begin: proc: %s init_optimalmpi\n",myidtxt); fflush(stderr);
  // if not setting geom, then say at least what's optimal

  // assumes all cpus are same size
  n1=(SFTYPE)N1;
  n2=(SFTYPE)N2;
  n3=(SFTYPE)N3;
  nc=(SFTYPE)numprocs;

  // assumes all cpus are same size
  n1=(SFTYPE)N1;
  n2=(SFTYPE)N2;
  n3=(SFTYPE)N3;
  nc=(SFTYPE)numprocs;
  
  if((COMPDIM<=2)&&(which==0) ){
    B=n2*nc*(16.*nc-1.);
    A=pow(2.0,2.0/3.0)*pow( (54.*nc*nc*n2*n1*(n2-n1)+sqrt(pow(54.*n1*(n1-n2)*n2*nc*nc,2.0)+4.*pow(3.*n1*B,3.0))),THIRD);
    C=sqrt((3.+4.*A/n1)/48.-B/A);
    
    fncpux1=1./8.-0.5*C+0.5*sqrt(1./8.+B/A-A/(12.*n1)+(-1./8.+4.*n2*nc/n1)/(4.*C));
    fncpux2=nc/fncpux1;
    
#define NZONES(x1,x2) (2.0*(n1*x1*(x1-1.)+n2*x2*(x2-1.)))
    // now force to be optimal integer, assumes NZONES is a smooth function in ncpux1 & ncpux2
    // finds optimal cost per total problem size for fixed N1,N2 per cpu
    tempx1[1]=floor(fncpux1);
    tempx2[1]=floor(nc/tempx1[1]);
    nz[1]=NZONES(tempx1[1],tempx2[1])/(tempx1[1]*tempx2[1]*n1*n2);
    
    tempx1[2]=ceil(fncpux1);
    tempx2[2]=ceil(nc/tempx1[2]);
    nz[2]=NZONES(tempx1[2],tempx2[2])/(tempx1[2]*tempx2[2]*n1*n2);
    
    tempx1[3]=floor(fncpux1);
    tempx2[3]=ceil(nc/tempx1[3]);
    nz[3]=NZONES(tempx1[3],tempx2[3])/(tempx1[3]*tempx2[3]*n1*n2);
    
    tempx1[4]=ceil(fncpux1);
    tempx2[4]=floor(nc/tempx1[4]);
    nz[4]=NZONES(tempx1[4],tempx2[4])/(tempx1[4]*tempx2[4]*n1*n2);
    
    lowestnz=nz[1]; lowesti=1;
    for(i=2;i<=4;i++){
      if(lowestnz>nz[i]){ lowestnz=nz[i]; lowesti=i;}
    }
    
    for(i=1;i<=4;i++){
      fprintf(stderr,"%d: %15.10g\n",i,nz[i]);
    }
    // choose optimal
    holdrealncpux1=fncpux1;
    fncpux1=tempx1[lowesti];
    fncpux2=nc/fncpux1;
    testncpux2=tempx2[lowesti];
    
    fprintf(stderr,"ncpux1: %15.10g ncpux2: %15.10g testncpux2: %15.10g\n",fncpux1,fncpux2,testncpux2);
    
    if(fabs(testncpux2-fncpux2)>1E-3){
      fprintf(stderr,"ncpux1: %15.10g ncpux2: %15.10g\n",fncpux1,fncpux2);
      fprintf(stderr,"give me %15.10g cpus not %15.10g cpus to be optimal mpi\n",fncpux1*testncpux2,nc);
      fprintf(stderr,"I'm just forcing you to choose a CPU# which is optimal and factors properly so grid has no holes");
      myexit(1);
    }
    else fncpux2=testncpux2;
    
    fprintf(stderr,"%d\n",lowesti);
    fprintf(stderr,"%15.10g %15.10g\n",fncpux1,fncpux2);
    
    // now set
    ncpux1=(int)fncpux1;
    ncpux2=(int)fncpux2;
    fprintf(stderr,"%d %d\n",ncpux1,ncpux2);
  }
  else if(which==1){
    // set and report whether choice was optimal

    if(COMPDIM<=2){

      fncpux1=(SFTYPE)ncpux1;
      fncpux2=(SFTYPE)ncpux2;
      fncpux3=(SFTYPE)ncpux3;
      
      if(fabs(nc/(fncpux1*fncpux2*fncpux3)-1.0)>1E-4){
	fprintf(stderr,"Number of cpus doesn't correspond correctly to cpu tile numbers\n");
	fprintf(stderr,"CPUS: %d  ncpux1: %d  ncpux2: %d ncpux3: %d\n",numprocs,ncpux1,ncpux2,ncpux3);
	myexit(1);
      }
      
      // report on size of grid and choice of N1, N2 for wanted global size
      shouldncpux1=sqrt(n2*fncpux2*nc/(n1*fncpux1));
      shouldncpux2=sqrt(n1*fncpux1*nc/(n2*fncpux2));
      
      shouldn1=(n1*fncpux1)/shouldncpux1;
      shouldn2=(n2*fncpux2)/shouldncpux2;
      
      if((shouldncpux1!=fncpux1)||(shouldncpux2!=ncpux2)||(shouldn1!=n1)||(shouldn2!=n2)){
	
	fprintf(stderr,"should use N1=%d N2=%d ncpux1=%d ncpux2=%d\n",(int)shouldn1,(int)shouldn2,(int)shouldncpux1,(int)shouldncpux2);
      }
      else{
	fprintf(stderr,"Good choice of N1,N2,ncpux1,ncpux2 for optimal MPI\n");
      }
    }
  }

  if((ncpux3==0)||(ncpux2==0)||(ncpux1==0)){
    fprintf(stderr,"You meant by 0 cpus, 1 cpu in that direction, irregardless of COMPDIM.  Fix it.\n");
    myexit(1);
  }
  
  fprintf(stderr,"Remember to consider CPU perf too!\n");
  fprintf(stderr,"4MB L2 cache min: min of N1~40 N2~40, best at N1~100+ N2~50+\n");

  fprintf(stderr,"end: proc: %s init_optimalmpi\n",myidtxt); fflush(stderr);  
}


void init_MPIgroup(void)
{
  int ranks[MAXCPUS] = {0}; 
  int i,j,k,numranks;

  MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);

  fprintf(stderr,"begin: proc: %s init_MPIgroup\n",myidtxt); fflush(stderr);

  // x1-inner
  j=0;
  for(i=0;i<numprocs;i++){
    if(i%ncpux1==0){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux2*ncpux3){
    fprintf(stderr,"problem with inner x1-group: numranks: %d ncpux2: %d ncpux3: %d ncpux2*ncpux3: %d\n",numranks,ncpux2,ncpux3,ncpux2*ncpux3);
    myexit(1);
  }
  // now ranks holds inner x1 boundary of cpus, and numranks holds number of such ranks

  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[2]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[2], &combound[2]); 

  // x1 - outer
  j=0;
  for(i=0;i<numprocs;i++){
    if(i%ncpux1==ncpux1-1){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux2*ncpux3){
    fprintf(stderr,"problem with outer x1-group: numranks: %d ncpux2*ncpux3: %d\n",numranks,ncpux2*ncpux3);
    myexit(1);
  }
  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[0]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[0], &combound[0]); 

  // x2 - inner
  j=0;
  for(i=0;i<numprocs;i++){
    if((i%(ncpux1*ncpux2))/ncpux1==0){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux1*ncpux3){
    fprintf(stderr,"problem with inner x2-group: numranks: %d ncpux1*ncpux3: %d\n",numranks,ncpux1*ncpux3);
    myexit(1);
  }
  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[1]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[1], &combound[1]); 

  // x2 - outer
  j=0;
  for(i=0;i<numprocs;i++){
    if((i%(ncpux1*ncpux2))/ncpux1==ncpux2-1){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux1*ncpux3){
    fprintf(stderr,"problem with outer x2-group: numranks: %d ncpux1*ncpux3: %d\n",numranks,ncpux1*ncpux3);
    myexit(1);
  }

  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[3]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[3], &combound[3]); 



  // x3 - inner
  j=0;
  for(i=0;i<numprocs;i++){
    if(i/(ncpux1*ncpux2)==0){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux1*ncpux2){
    fprintf(stderr,"problem with inner x3-group: numranks: %d ncpux1*ncpux2: %d\n",numranks,ncpux1*ncpux2);
    myexit(1);
  }
  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[5]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[5], &combound[5]); 

  // x3 - outer
  j=0;
  for(i=0;i<numprocs;i++){
    if(i/(ncpux1*ncpux2)==ncpux3-1){
      ranks[j]=i;
      j++;
    }
  }
  numranks=j;
  if(numranks!=ncpux1*ncpux2){
    fprintf(stderr,"problem with outer x3-group: numranks: %d ncpux1*ncpux2: %d\n",numranks,ncpux1*ncpux2);
    myexit(1);
  }

  // now create group and communicator
  MPI_Group_incl(MPI_GROUP_WORLD, numranks, ranks, &grprem[4]);
  MPI_Comm_create(MPI_COMM_WORLD, grprem[4], &combound[4]); 

  // 0: right 1: up 2: left 3: down 4: out 5: in(as in bound.c)

  // when using these communicators, must make sure the call to a communication using it isn't done by a non-member cpu!
  // (above: stupid, I know, should just skip if non-member cpu tries a function)
  fprintf(stderr,"end: proc: %s init_MPIgroup\n",myidtxt); fflush(stderr);
}



#endif

void init_placeongrid(void)
{
  int i,m;
  int N[3+1];
  int numbercpu[3+1];

  fprintf(log_file,"begin: init_placeongrid ... "); fflush(log_file);

  N[1]=N1;
  N[2]=N2;
  N[3]=N3;

  numbercpu[1]=ncpux1;
  numbercpu[2]=ncpux2;
  numbercpu[3]=ncpux3;

  mycpupos[1]=myid%ncpux1;
  mycpupos[2]=(int)((myid%(ncpux1*ncpux2))/ncpux1);
  mycpupos[3]=(int)(myid/(ncpux1*ncpux2));
    
  for(m=1;m<=3;m++){
    startpos[m]=mycpupos[m]*N[m];
    endpos[m]=(mycpupos[m]+1)*N[m]-1;

    // add up sizes for total size of grid
    totalsize[m]=0;
    itotalsize[m]=0;
    for(i=0;i<numbercpu[m];i++){
      totalsize[m]+=N[m];
      itotalsize[m]+=N[m];
    }
  }

  realtotalzones=totalzones=totalsize[1]*totalsize[2]*totalsize[3];
  itotalzones=itotalsize[1]*itotalsize[2]*itotalsize[3];

  // send/receive directions for this cpu(interior only)
 /////////////// standard interior MPI data transfer setup
  srdir[0]=srdir[1]=srdir[2]=srdir[3]=srdir[4]=srdir[5]=0;
  // see where this cpu needs to send/recv

  // figure out left/right send/recv
  if(mycpupos[1]>0){
    srdir[2]=1; // do left
  }
  if(mycpupos[1]<ncpux1-1){
    srdir[0]=1; // do right
  }

  // figure out up/down send/recv
  if(mycpupos[2]>0){
    srdir[1]=1; // do up
  }
  if(mycpupos[2]<ncpux2-1){
    srdir[3]=1; // do down
  }
  // figure out in/out send/recv
  if(mycpupos[3]>0){
    srdir[5]=1; // do in
  }
  if(mycpupos[3]<ncpux3-1){
    srdir[4]=1; // do out
  }

  fprintf(log_file,"totalzones: %d\n",totalzones);

  fprintf(log_file,"end: init_placeongrid\n"); fflush(log_file);
}

void init_visc(void)
{
  int i,j,k,l,m;

  fprintf(log_file,"begin: init_visc ... "); fflush(log_file);

#if(VISCMEM)
  // init visc param(part the never changes)(must be consistent with step.c)
  LOOPF{
    if(vreal==1){
      // alpha-parameterization
      // $\nu = \alpha r sin(\theta)^(n_real)
      nu_fact[k][j][i]=x[2][1][i]*pow(G4(2,j),n_real);
    }

    if(vreal==2){
      //Igumenshchev alpha param: $\nu=\alpha*cs^2/Omega_{k}$
      // here by Keplerian we mean spherical rotation, even if fluid is meant to rotate on cylinders.
      //nu_fact[k][j][i]=pow(x[2][1][i],1.5)/sqrt(GM);
      // below was being used for all IGU runs, when IGU really uses above!
      //nu_fact[k][j][i]=pow(x[2][1][i]*G4(2,j),1.5)/sqrt(GM);
      // below for PW pot, use IGU like below
      nu_fact[k][j][i]=pow(x[2][1][i],1.5)*(1.0-rgp/x[2][1][i])/sqrt(GM);
    }

    if(vreal==3){ 
      // Stone et al.
      // Run A through J with different alpha (J diff nu_real)
      // assumes using units=1 in tori1sol
      nu_fact[k][j][i]=R0*R0*Omega0/rho0 ;
      
      // Run K with different alpha
      //nu_fact[k][j][i]=alpha_real*sqrt(x[2][1][i]*G4(2,j)) ;
    }
    if(vreal==4){
      // below was being used for all IGU runs, when IGU really uses above!
      nu_fact[k][j][i]=pow(x[2][1][i]*G4(2,j),1.5)/sqrt(GM);
      // below for PW pot, use IGU like below
      //nu_fact[k][j][i]=pow(x[2][1][i],1.5)*(1.0-rgp/x[2][1][i])/sqrt(GM);
    }
    if(vreal==5){ // constant
      nu_fact[k][j][i]=1.0;
    }

    if(analoutput==6){
      // for test of the visc real code
      nu_fact[k][j][i]=1.0;
    }

  }
#endif
  fprintf(log_file,"end: init_visc\n"); fflush(log_file);
}



void init_res(void)
{
  int i,j,k,l,m;

  fprintf(log_file,"begin: init_res ... "); fflush(log_file);

#if(RESMEM)
  // init visc param(part the never changes)(must be consistent with step.c)
  LOOPF{
    if(rreal==1){
      // gammie factor
      //nu_res_fact[k][j][i]=dx[2][1][i];
      if(OARC11(k,j,i)<OARC12(k,j,i)){ // i.e. dx1>dx2
	nu_res_fact[k][j][i]=1.0/pow(OARC12(k,j,i),1.0);
      }
      else{
	nu_res_fact[k][j][i]=1.0/pow(OARC11(k,j,i),1.0);
      }
    }

    if(rreal==2){
      // SP00 factor
      // radial control only
      //nu_res_fact[k][j][i]=(dx[1][2][i]*dx[1][2][i]);
      // theta control only
      //      nu_res_fact[k][j][i]=(x[2][1][i]*dx[2][2][j]*x[2][1][i]*dx[2][2][j]);
      // do a zone by zone measure(whatever smaller)
      if(OARC11(k,j,i)<OARC12(k,j,i)){ // i.e. dx1>dx2
	nu_res_fact[k][j][i]=1.0/pow(OARC12(k,j,i),2.0);
      }
      else{
	if(N3==1){
	  nu_res_fact[k][j][i]=1.0/pow(OARC11(k,j,i),2.0);
	}
	else{
	  if(OARC11(k,j,i)*ODX(2,3,k)<OARC13(k,j,i)*ODX(2,3,k)){ // i.e. dx1>dx3
	    nu_res_fact[k][j][i]=1.0/pow(OARC13(k,j,i)*ODX(2,3,k),2.0);
	  }
	}
      }
      // may want to try average(nah)
    }
    if(rreal==3){ // constant
      nu_res_fact[k][j][i]=1.0;
    }

  }
#endif
  fprintf(log_file,"end: init_res\n"); fflush(log_file);
}


void init_floor(void)
{
  int i,j,k,l,m;

  fprintf(log_file,"begin: init_floor ... "); fflush(log_file);

  if(FLOORDUMPFLAG){
    // init floor variables
    LOOP{
      for(l=1;l<=NUMSCA;l++){
	floorvars[l][k][j][i]=0.0;
      }
      for(l=1;l<=NUMVEC;l++){
	floorvar0[l][k][j][i]=0.0;
	for(m=1;m<=3;m++){
	  floorvarv[l][m][k][j][i]=0.0;
	}
      }
    }
  }
  // initialize floor checks
  for(i=0;i<=NUMFLOORVAR;i++){
    floorlowest[i]=1.0;
    wherelowest[i]=-1;
  }
  for(i=0;i<=NUMLOSSVAR;i++){
    floors[i]=0.0;
  }

  fprintf(log_file,"end: init_floor\n"); fflush(log_file);
}

void init_inflows(void)
{
  int i,j,k,l,m;

  fprintf(log_file,"begin: init_inflows ... "); fflush(log_file);

  // init inflows variables
  for(i=0;i<=NUMLOSSVAR;i++){
    inflows[i]=0.0;
  }
  fprintf(log_file,"end: init_inflows\n"); fflush(log_file);
}
void init_radiations(void)
{
  int i,j,k,l,m;

  fprintf(log_file,"begin: init_radiations ... "); fflush(log_file);

  // init inflows variables
  for(i=0;i<=NUMLOSSVAR;i++){
    radiations[i]=0.0;
  }

  fprintf(log_file,"end: init_radiations\n"); fflush(log_file);
}


// 1/2/3-d valid
int init_mem(void)
{
  fprintf(log_file,"begin: init_mem ... "); fflush(log_file);
  /* END: Allocate memory for parameter structure */
  fprintf(log_file,"end: init_mem\n"); fflush(log_file);
  return(0);
}

void init_genfiles(int gopp)
{
  char temps[MAXFILENAME];
  char extension[MAXFILENAME];

  fprintf(stderr,"begin: init_genfiles ... "); fflush(stderr);

  if(gopp==1){
    strcpy(extension,PPEXT);
  }
  else if(gopp==0){
    strcpy(extension,OUTEXT);
  }

  // always have fail and general log open

  sprintf(temps,"%s0_fail%s%s",DATADIR,extension,myidtxt) ;
    
  
  if((fail_file=fopen(temps,"wt"))==NULL){
    fprintf(stderr,"fail: Cannot open: %s\n",temps);
    exit(1);
  }
 
  sprintf(temps, "%s0_gridsigma%s", DATADIR, extension);
  if ((grid_file = fopen(temps, "wt")) == NULL) {
    fprintf(stderr, "fail: Cannot open: %s\n", temps);
    exit(1);
  }

#if(RADTEST)
  sprintf(temps, "%s0_radtest%s", DATADIR, extension);
  if ((rad_file = fopen(temps, "wt")) == NULL) {
    fprintf(stderr, "fail: Cannot open: %s\n", temps);
    exit(1);
  }
#endif
 
  sprintf(temps,"%s0_log%s%s",DATADIR,extension,myidtxt) ;
    
  if((log_file=fopen(temps,"wt"))==NULL){
    fprintf(stderr,"log: Cannot open: %s\n",temps);
    exit(1);
  }
  if(myid<=0){
    sprintf(temps,"%s0_logfull%s",DATADIR,extension) ;
    
    if((logfull_file=fopen(temps,"wt"))==NULL){
      fprintf(stderr,"logfull: Cannot open: %s\n",temps);
      exit(1);
    }
  }

  fprintf(stderr,"end: init_genfiles\n"); fflush(stderr);
}

// assume periodicix1,skipix1,etc. are consistent with read boundary conditions on grid as needed by correct code
int init_runpar(void)
{
  FILE* data;
  FILE* in[2];
  char fname[2][200];
  int dumi[50];
  SFTYPE dumf[50];
#if(FLOATTYPE==0)
  FTYPE dumfread;
#else
  SFTYPE dumfread;
#endif
  char ch;
  int i,j,k,l,m;
  char crap[50],testcrap[50];
  char extension[MAXFILENAME];
  char extension2[MAXFILENAME];
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
#if(USEMPI&&USEROMIO)
  MPI_Datatype newtype;
  MPI_File fh[2];
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
  char truemyidtxt[100];
  int pdiddy,pdiddystart,pdiddyto;
  int nextbuf;
  int range;

  // assume: grid1, grid2, and gparam for parameter data

  // remember to make same changes to this and output in diag.c/init.c!!
  fprintf(log_file,"begin: init_runpar ... "); fflush(log_file);

  if(directinput>0){
    strcpy(extension,PAREXT);
    strcpy(extension2,PAR2EXT);
  }
  else{
    strcpy(extension,INEXT);
  }

  if(mpicombine==0){
    strcpy(truemyidtxt,myidtxt);
  }
  else strcpy(truemyidtxt,"");


  if((runtype==22)||(runtype==2)){

  for(i=0;i<numprocs;i++){ // READ of parameter file, one cpu at a time
    if(i==myid){
      sprintf(fname[0],"%s0_gparam%s",DATADIR,extension2) ;
      
      if((in[0]=fopen(fname[0],"rt"))==NULL){
	fprintf(fail_file,"gparam: Cannot open: %s\n",fname[0]);
	myexit(1);
      }
      fprintf(log_file,"reading par file\n"); fflush(log_file);
      
      
      l=0;
      ch=fgetc(in[l]);
      if(ch=='#'){
	while(fgetc(in[l])!='\n');
      }
      else{
	fprintf(fail_file,"1: gparm: global parm file doesn't match expected format\n");
	myexit(1);
      }
      
      fscanf(in[l],HEADER3_S0,&dumi[0],&dumi[1]);
      if((dumi[0]!=PVER)||(dumi[1]!=PTYPE)){
	fprintf(fail_file,"Expected pver/ptype: %d %d got %d %d\n",PVER,PTYPE,dumi[0],dumi[1]);
	myexit(6);
      }
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S1,&dumi[0],&dumi[1],&dumi[2]);
      if( (dumi[0]!=totalsize[1])|| (dumi[1]!=totalsize[2])|| (dumi[2]!=totalsize[3])){
	fprintf(fail_file,"2: gparm: Gridsize in gparm file needs to match global.h file\n");
	fprintf(fail_file,"got: N1: %d N2: %d N3: %d\n",dumi[0],dumi[1],dumi[2]);
	fprintf(fail_file,"expected: N1: %d N2: %d N3: %d\n",totalsize[1],totalsize[2],totalsize[3]);
	myexit(1);
      }
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S2,&x1in,&x2in,&x3in,&x1out,&x2out,&x3out);
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S3,&rg,&rgp,&cs,&coolfact,&gam,&alpha_real0,&n_real);
      //fprintf(stdout,"%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",rg,rgp,cs,coolfact,gam,alpha_real0,n_real);
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S4,&nu_vnr,&nu_l,&nu_ten,&cour,&cour2);
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S5,&GRAVC,&MASSBH,&invsol2,&blackholejz);
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S6,&tstart,&tf,&tavgi,&tavgf,&numavg,&timagescale);
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S7,&DTl,&DTd,&DTi,&DTloss,&DTfloor,&DTtimestep,&DTpd,&DTener,&DTfld,&DTmode);
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S8,&resist_real0,&nu_sh);
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S9,&vg[1],&vg[2],&vg[3]);
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S10,&dumi[0],&dumi[1],&analoutput,&DYNAMICMM);
      if( (dumi[0]!=COORD)|| (dumi[1]!=(REALNUMVEC==NUMVEC)) ){
	fprintf(fail_file,"2.1: gparm: coord needs to match that in global.h\n");
	fprintf(fail_file,"got: coord: %d RNV==NV: %d\n",dumi[0],dumi[1]);
	fprintf(fail_file,"expected: coord: %d RNV=NV: %d \n",COORD,(REALNUMVEC==NUMVEC));
	myexit(1);
      }
      
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S11,&trans,&transx1,&transrhox1,&transiex1,&transv1x1,&transv2x1,&transv3x1,&transmagx1,&transx2,&transrhox2,&transiex2,&transv1x2,&transv2x2,&transv3x2,&transmagx2,&transx3,&transrhox3,&transiex3,&transv1x3,&transv2x3,&transv3x3);
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S12,&press,&pressx1,&pressx2,&pressx3);
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S13,&mag,&transbzx,&transbzy,&stepmocct,&mocctvx1,&mocctvx2,&mocctvx3,&mocctbx1,&mocctbx2,&mocctbx3);
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S14,&ie,&visc_art,&visc_real,&vreal,&vischeat,&mdotin,&cool,&res_real,&rreal,&resheat,&advint,&kever);
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S15,&intix1,&intox1,&intix2,&intox2,&intix3,&intox3);
      while((ch=fgetc(in[l]))!='\n');
      fscanf(in[l],HEADER3_S16,&nonunigridx1,&nonunigridx2,&nonunigridx3,&simplebc,&bcix1,&bcox1,&bcix2,&bcox2,&bcix3,&bcox3);
      
      fclose(in[0]);
    }
#if(USEMPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }



  }
  if(runtype==2){

    sprintf(fname[0],"%s0_grid1%s%s",DATADIR,extension,myidtxt) ;
    sprintf(fname[1],"%s0_grid2%s%s",DATADIR,extension,myidtxt) ;
    
    if(mpicombine==0){
      if((in[0]=fopen(fname[0],"rt"))==NULL){
	fprintf(fail_file,"grid1: Cannot open: %s\n",fname[0]);
	myexit(1);
      }
      
      if((in[1]=fopen(fname[1],"rt"))==NULL){
	fprintf(fail_file,"grid2: Cannot open: %s\n",fname[1]);
	myexit(1);
      }
    }
    else{

      numcolumns=21+NUMSCA*3+NUMVEC*3;
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
      
      MPI_Type_create_darray(numprocs, myid, ndims, array_of_gsizes, 
			     array_of_distribs, array_of_dargs,
			     array_of_psizes, order, MPI_FTYPE, &newtype);
      MPI_Type_commit(&newtype);
      MPI_Type_size(newtype, &bufcount);
      bufcount = bufcount/sizeof(FTYPE);
      writebuf = (FTYPE *) malloc(bufcount * sizeof(FTYPE));
      // setup file handler
      
      MPI_File_open(MPI_COMM_WORLD, fname[0], MPI_MODE_CREATE | MPI_MODE_RDWR, 
		    MPI_INFO_NULL, &fh[0]);
      MPI_File_set_view(fh[0], 0, MPI_FTYPE, newtype, "native", MPI_INFO_NULL);
      
      MPI_File_open(MPI_COMM_WORLD, fname[1], MPI_MODE_CREATE | MPI_MODE_RDWR, 
		    MPI_INFO_NULL, &fh[1]);
      MPI_File_set_view(fh[1], 0, MPI_FTYPE, newtype, "native", MPI_INFO_NULL);
      // all that needs to be done now is fill writebuf with the data and write back to normal data arrays
#endif
      }
      else if(USEJONIO){
	  jonio_init_fp(&in[0],2,fname[0]);
	  jonio_init_fp(&in[1],2,fname[1]);
	  if(sizeof(FTYPE)==sizeof(double)){
	      jonio_init_mem(numcolumns,sizeof(double),0,0,&jonio,0,0,&writebuf,0,0,&tempbuf);
	  }
	  else if(sizeof(FTYPE)==sizeof(float)){
	      jonio_init_mem(numcolumns,sizeof(float),0,&jonio,0,0,&writebuf,0,0,&tempbuf,0);
	  }
	  else if(sizeof(FTYPE)==sizeof(unsigned char)){
	      jonio_init_mem(numcolumns,sizeof(unsigned char),&jonio,0,0,&writebuf,0,0,&tempbuf,0,0);
	  }
      }
#endif
    }
    
    
    for(l=1;l<=NUMGRID;l++){ 
      fprintf(log_file,"reading grid file %d\n",l); fflush(log_file);
      
      if(mpicombine==0){
	ch=fgetc(in[l-1]);
	if(ch=='#'){
	  while(fgetc(in[l-1])!='\n');
	}
	else{
	  fprintf(fail_file,"1: grid data file doesn't match expected format: %d\n",l);
	  myexit(1);
	}
	fscanf(in[l-1],"%d %d\n",&dumi[0],&dumi[1]);
	if((dumi[0]!=GRIDVER)||(dumi[1]!=GRIDTYPE)){
	  fprintf(fail_file,"expected gridver/gridtype: %d %d got %d %d\n",GRIDVER,GRIDTYPE,dumi[0],dumi[1]);
	  myexit(6);
	}
	ch=fgetc(in[l-1]);
	if(ch=='#'){
	  while(fgetc(in[l-1])!='\n');
	}
	else{
	  fprintf(fail_file,"1.1: grid data file doesn't match expected format: %d\n",l);
	  myexit(1);
	}
      }
      if(mpicombine==0){
	LOOPF{
	  fscanf(in[l-1],HEADER4_S,
		 &dumi[0],&dumi[1],&dumi[2],&dumi[3],
		 &dx[l][1][i],&dx[l][2][j],&dx[l][3][k],
		 &x[l][1][i],&x[l][2][j],&x[l][3][k],
#if(LOWMEMMODE==0)
		 &g[l][1][i],&dg[l][1][i],&g[l][2][i],&dg[l][2][i],&g[l][3][i],&dg[l][3][i],&g[l][4][j],&dg[l][4][j],
		 &dvl[l][1][i],&dvl[l][2][j],&dvl[l][3][k]
#else
		 &dumfread,&dumfread,&dumfread,&dumfread,&dumfread,&dumfread,&dumfread,&dumfread,
		 &dumfread,&dumfread,&dumfread
#endif
		 );
	  if((PUREBC==0)&&(BOUNDTYPE==1)){
	    for(m=1;m<=NUMSCA;m++){
	      fscanf(in[l-1],"%d ",&dumi[4+(m-1)*3]);
	      fscanf(in[l-1],"%d ",&dumi[5+(m-1)*3]);
	      fscanf(in[l-1],"%d ",&dumi[6+(m-1)*3]);
	    }
	    for(m=1;m<=NUMVEC;m++){
	      fscanf(in[l-1],"%d ",&dumi[4+NUMSCA*3+(m-1)*3]);
	      fscanf(in[l-1],"%d ",&dumi[5+NUMSCA*3+(m-1)*3]);
	      fscanf(in[l-1],"%d",&dumi[6+NUMSCA*3+(m-1)*3]);
	      if(m==NUMVEC) fscanf(in[l-1],"\n"); else fscanf(in[l-1]," ");
	    }
	    // must do it this way to read in properly and get memory writting correct too
	    for(m=1;m<=NUMSCA;m++){
	      bcs[m][1][k][j][i]=(short)dumi[4+(m-1)*3];
	      bcs[m][2][k][j][i]=(short)dumi[5+(m-1)*3];
	      bcs[m][3][k][j][i]=(short)dumi[6+(m-1)*3];
	    }
	    for(m=1;m<=NUMVEC;m++){
	      bcv[m][1][k][j][i]=(short)dumi[4+NUMSCA*3+(m-1)*3];
	      bcv[m][2][k][j][i]=(short)dumi[5+NUMSCA*3+(m-1)*3];
	      bcv[m][3][k][j][i]=(short)dumi[6+NUMSCA*3+(m-1)*3];
	    }
	  }
	  else if((BOUNDTYPE==2)||(BOUNDTYPE==3)){
	    for(m=1;m<=15;m++){
	      fscanf(in[l-1],"%d ",&dumi[4+m-1]);
	    }
	    bzmask[k][j][i]=(short)dumi[4];
	    bzmaskin[k][j][i]=(short)dumi[5];
	    //	maskv2[k][j][i]=(short)dumi[6];
	    //maskv3[k][j][i]=(short)dumi[7];
	  }
	  else{
	    for(m=1;m<=15;m++){
	      fscanf(in[l-1],"%d ",&dumi[4+m-1]);
	    }
	  }
	}
      }
      else{	
#if(USEMPI)
	  if(USEROMIO){
#if(USEROMIO)
	      MPI_File_read_all(fh[l-1], writebuf, bufcount, MPI_FTYPE, &status);
	      MPI_File_close(&fh[l-1]);
#endif
	  }
	  else if(USEJONIO){
	      jonio_seperate(1,MPI_FTYPE,numcolumns,sizeof(FTYPE),in[l-1],jonio,writebuf,tempbuf);
	  }
#endif
	LOOP{// must figure out how to read in properly the full grid
	  BUFFERINIT;
	  dumi[0]=writebuf[BUFFERMAP];
	  dumi[1]=writebuf[BUFFERMAP];
	  dumi[2]=writebuf[BUFFERMAP];
	  dx[l][1][i]=writebuf[BUFFERMAP];
	  dx[l][2][i]=writebuf[BUFFERMAP];
	  dx[l][3][i]=writebuf[BUFFERMAP];
	  x[l][1][i]=writebuf[BUFFERMAP];
	  x[l][2][i]=writebuf[BUFFERMAP];
	  x[l][3][i]=writebuf[BUFFERMAP];
#if(LOWMEMMODE==0)
	  g[l][1][i]=writebuf[BUFFERMAP];
	  dg[l][2][i]=writebuf[BUFFERMAP];
	  g[l][2][i]=writebuf[BUFFERMAP];
	  dg[l][2][i]=writebuf[BUFFERMAP];
	  g[l][3][i]=writebuf[BUFFERMAP];
	  dg[l][3][i]=writebuf[BUFFERMAP];
	  g[l][4][i]=writebuf[BUFFERMAP];
	  dg[l][4][i]=writebuf[BUFFERMAP];
	  dvl[l][1][i]=writebuf[BUFFERMAP];
	  dvl[l][2][i]=writebuf[BUFFERMAP];
	  dvl[l][3][i]=writebuf[BUFFERMAP];
#else
	  dumf[1]=writebuf[BUFFERMAP];
	  dumf[1]=writebuf[BUFFERMAP];
	  dumf[1]=writebuf[BUFFERMAP];
	  dumf[1]=writebuf[BUFFERMAP];
	  dumf[1]=writebuf[BUFFERMAP];
	  dumf[1]=writebuf[BUFFERMAP];
	  dumf[1]=writebuf[BUFFERMAP];
	  dumf[1]=writebuf[BUFFERMAP];
	  dumf[1]=writebuf[BUFFERMAP];
	  dumf[1]=writebuf[BUFFERMAP];
	  dumf[1]=writebuf[BUFFERMAP];
#endif
	  if((PUREBC==0)&&(BOUNDTYPE==1)){
	    for(m=1;m<=NUMSCA;m++){
	      bcs[m][1][k][j][i]=writebuf[BUFFERMAP];
	      bcs[m][2][k][j][i]=writebuf[BUFFERMAP];
	      bcs[m][3][k][j][i]=writebuf[BUFFERMAP];
	    }
	    for(m=1;m<=NUMVEC;m++){
	      bcv[m][1][k][j][i]=writebuf[BUFFERMAP];
	      bcv[m][2][k][j][i]=writebuf[BUFFERMAP];
	      bcv[m][3][k][j][i]=writebuf[BUFFERMAP];
	    }
	  }
	  else if((BOUNDTYPE==2)||(BOUNDTYPE==3)){
	    for(m=1;m<=15;m++){
	      dumi[4+m-1]=writebuf[BUFFERMAP];
	    }
	    bzmask[k][j][i]=(short)dumi[4];
	    bzmaskin[k][j][i]=(short)dumi[5];
	    //	maskv2[k][j][i]=(short)dumi[6];
	  //maskv3[k][j][i]=(short)dumi[7];
	  }
	  else{
	    for(m=1;m<=15;m++){
	      dumi[4+m-1]=writebuf[BUFFERMAP];
	    }
	  }
	}
      }
      if( (dumi[0]!=l)|| (dumi[1]!=k)|| (dumi[2]!=j+startpos[2])|| (dumi[3]!=i)){
	fprintf(fail_file,"2: grid data file doesn't match expected format: %d %d %d %d\n",l,k,j+startpos[2],i);
	fprintf(fail_file,"%d %d %d %d\n",dumi[0],dumi[1],dumi[2],dumi[3]);
	myexit(1);
	}
    }// over grids
  }// end if runtype==2
  if(mpicombine==0){
    for(l=1;l<=NUMGRID;l++){
      fclose(in[l-1]);
    }
  }
  else{
#if(USEMPI)
      if(USEROMIO){
#if(USEROMIO)
	  free(writebuf);
	  MPI_Type_free(&newtype);
#endif
      }
      else if(USEJONIO){
	  for(l=1;l<=NUMGRID;l++){
	      jonio_seperate(2,MPI_FTYPE,numcolumns,sizeof(FTYPE),in[l-1],jonio,writebuf,tempbuf); 
	  }
	  jonio_seperate(3,MPI_FTYPE,numcolumns,sizeof(FTYPE),0,jonio,writebuf,tempbuf); 
      }
#endif
  }

  fprintf(log_file,"end: init_runpar\n"); fflush(log_file);  
  return(0);
}



int init_rundat(int dump_cnt, int which)
{
  int N[4];
  FILE *data,*header;
  char fname[MAXFILENAME];
  int dumi[50];
  FTYPE dumf[50];
  SFTYPE dumfs[10];
  char ch;
  int i,j,k,l,m,p;
  char crap[50],testcrap[50];
  FTYPE (*sin)[N3M][N2M][N1M] ;
  FTYPE (*vin)[3][N3M][N2M][N1M] ;
  char dfheader[MAXFILENAME];
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
  char truemyidtxt[100];
  int range;
  int nextbuf;

  fprintf(log_file,"begin: init_rundat ... "); fflush(log_file);
  // for now, assume filename is: dump for data

  // remember to make same changes to this and output in diag.c/init.c!!

  // can only read in non-interp'ed data
  if(which==DTYPE){
    sin=s;
    vin=v;
    strcpy(dfheader,"dump");
    if(analoutput==13) numcolumns=7;
    else numcolumns=9;
  }
  else if(which==EXPANDTYPE){ // for post proc feature(doexpand)
    sin=s;
    vin=v;
    strcpy(dfheader,"toedump");
    which=DTYPE; // otherwise treat as dump
  }
  else if(which==PDTYPE){
    sin=s;
    vin=v;
    strcpy(dfheader,"pdump");
    if(analoutput==13) numcolumns=7;
    else numcolumns=9;
  }  
  else if(which==ADTYPE){
    sin=sanal;
    vin=vanal;
    strcpy(dfheader,"adump");
    if(analoutput==13) numcolumns=7;
    else numcolumns=9;
  }
  else if(which==FLTYPE){
    sin=floorvars;
    vin=floorvarv;
    //floorvar0 used natively
    strcpy(dfheader,"floor");
    if(analoutput==13) numcolumns=7;
    else numcolumns=11;
  }
  else if(which==NPTYPE){
    strcpy(dfheader,"npdump");
    numcolumns=7;
  }
  else if(which==AVG2DTYPE){
    // use avg2_3[NUMAVG2_3][..][..], 2D data with NUMAVG2_3 types
    strcpy(dfheader,"0_avg2d");
    numcolumns=15;
  }
  else{
    fprintf(fail_file,"not setup for which=%d in init_rundat\n",which);
    myexit(1);
  }


  if(mpicombine==0){
    strcpy(truemyidtxt,myidtxt);
  }
  else strcpy(truemyidtxt,"");


  if(directinput==0){
    sprintf(fname,"%s0_%s%s%s",DATADIR,dfheader,INEXT,truemyidtxt) ; 
  }
  else{
    if(dump_cnt>=0){
      sprintf(fname,"%s%s%04d%s%s",DATADIR,dfheader,dump_cnt,DATEXT,truemyidtxt) ; 
    }
    else{
      sprintf(fname,"%s%s%s%s",DATADIR,dfheader,DATEXT,truemyidtxt) ; 
    }
  }
  fprintf(log_file,"Start Reading in %s\n",fname); fflush(log_file);

  if(mpicombine==0){
    if(binaryoutput==0){
      if((data=fopen(fname,"rt"))==NULL){
	fprintf(fail_file,"data: Cannot open: %s\n",fname);
	myexit(1);
      }
      if(myid==0) header=data;
    }
    else{
      if((data=fopen(fname,"rb"))==NULL){
	fprintf(fail_file,"data: Cannot open: %s\n",fname);
	myexit(1);
      }
    }
  }
  else{
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
    
    MPI_Type_create_darray(numprocs, myid, ndims, array_of_gsizes, 
			   array_of_distribs, array_of_dargs,
			   array_of_psizes, order, MPI_FTYPE, &newtype);
    MPI_Type_commit(&newtype);
    MPI_Type_size(newtype, &bufcount);
    bufcount = bufcount/sizeof(FTYPE);
    writebuf = (FTYPE *) malloc(bufcount * sizeof(FTYPE));
    // setup file handler
    
    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_RDWR, 
		  MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, 0, MPI_FTYPE, newtype, "native", MPI_INFO_NULL);
      // all that needs to be done now is fill writebuf with the data and write back to normal data arrays
#endif
      }
      else if(USEJONIO){
	  jonio_init_fp(&data,2,fname); // 2 for read
	  if(sizeof(FTYPE)==sizeof(double)){
	      jonio_init_mem(numcolumns,sizeof(double),0,0,&jonio,0,0,&writebuf,0,0,&tempbuf);
	  }
	  else if(sizeof(FTYPE)==sizeof(float)){
	      jonio_init_mem(numcolumns,sizeof(float),0,&jonio,0,0,&writebuf,0,0,&tempbuf,0);
	  }
	  else if(sizeof(FTYPE)==sizeof(unsigned char)){
	      jonio_init_mem(numcolumns,sizeof(unsigned char),&jonio,0,0,&writebuf,0,0,&tempbuf,0,0);
	  }
      }
#endif
  }


  fprintf(log_file,"Start Reading in header\n"); fflush(log_file);
  if(myid==0){
      ///////////////
      // BEGIN PROCESS HEADER
      ///////////////
      
      if((mpicombine==1)||((mpicombine==0)&&(binaryoutput==1))){
	  // open header file
	  strcat(fname,".head");
	  if((header=fopen(fname,"rt"))==NULL){
	      fprintf(fail_file,"header: Cannot open: %s\n",fname);
	      myexit(1);
	  }
      }

      // check version info
      ch=fgetc(header);
      if(ch=='#'){
	  while(fgetc(header)!='\n');
      }
      else{
	  fprintf(fail_file,"1:read: data file doesn't match expected format\n");
	  myexit(1);
      }
      // version header
      fscanf(header,"%d %d\n",&dumi[0],&dumi[1]);
      
      if(which==DTYPE){
	  if((dumi[0]!=DVER)||(dumi[1]!=DTYPE)){
	      fprintf(fail_file,"read: Expected DVER/DTYPE %d %d got %d %d\n",DVER,DTYPE,dumi[0],dumi[1]);
	      myexit(6);
	  }
      }
      else if(which==ADTYPE){
	  if((dumi[0]!=ADVER)||(dumi[1]!=ADTYPE)){
	      fprintf(fail_file,"read: Expected ADVER/ADTYPE %d %d got %d %d\n",ADVER,ADTYPE,dumi[0],dumi[1]);
	      myexit(6);
	  }
      }
      else if(which==PDTYPE){
	  if((dumi[0]!=PDVER)||(dumi[1]!=PDTYPE)){
	      fprintf(fail_file,"read: Expected PDVER/PDTYPE %d %d got %d %d\n",PDVER,PDTYPE,dumi[0],dumi[1]);
	      myexit(6);
	  }
      }
      else if(which==FLTYPE){
	  if((dumi[0]!=FLVER)||(dumi[1]!=FLTYPE)){
	      fprintf(fail_file,"read: Expected FLVER/FLTYPE %d %d got %d %d\n",FLVER,FLTYPE,dumi[0],dumi[1]);
	      myexit(6);
	  }
      }
      else if(which==NPTYPE){
	  if((dumi[0]!=NPVER)||(dumi[1]!=NPTYPE)){
	      fprintf(fail_file,"read: Expected NPVER/NPTYPE %d %d got %d %d\n",NPVER,NPTYPE,dumi[0],dumi[1]);
	      if( (dumi[1]==NPTYPE)&&(dumi[0]==1)){
		  npoldschool=1;
		  fprintf(fail_file,"read: Determined got %d %d\n",dumi[0],dumi[1]);
	      }
	      else if( (dumi[1]==NPTYPE)&&(dumi[0]==0)){
		  npoldschool=0;
		  fprintf(fail_file,"read: Determined got %d %d\n",dumi[0],dumi[1]);
	      }
	      else{
		  myexit(6);
	      }
	  }
      }
      else if(which==AVG2DTYPE){
	  if((dumi[0]!=AVG2DVER)||(dumi[1]!=AVG2DTYPE)){
	      fprintf(fail_file,"read: Expected AVG2DVER/AVG2DTYPE %d %d got %d %d\n",AVG2DVER,AVG2DTYPE,dumi[0],dumi[1]);
	      if( (dumi[1]==AVG2DTYPE)&&(dumi[0]==1)){
		  avg2doldschool=1;
		  fprintf(fail_file,"read: Determined got %d %d\n",dumi[0],dumi[1]);
	      }
	      else{
		  myexit(6);
	      }
	  }
      }
      else{
	  fprintf(fail_file,"not setup for which=%d in init_rundat\n",which);
	  myexit(1);
      }
      
      
      // continue with data file and normal header
      ch=fgetc(header);
      if(ch=='#'){
	  while(fgetc(header)!='\n');
      }
      else{
	  fprintf(fail_file,"1.1:read: data file doesn't match expected format\n");
	  myexit(1);
      }
      if(POSTPROC==0){
	  if(OLDSCHOOL==1){
	      if( (which==DTYPE)||(which==PDTYPE) ) fscanf(header,INPUT1OLD,&t,&dumi[0]);
	      else fscanf(header,INPUT1OLD,&dumfs[0],&dumi[0]);
	  }
	  else if(OLDSCHOOL==0){
	      if( (which==DTYPE)||(which==PDTYPE) ) fscanf(header,INPUT1,&t,&dumi[0],&dumi[1]);
	      else fscanf(header,INPUT1,&dumfs[0],&dumi[0],&dumi[1]);
	  }
      }
      else{ // always read in time if postproc==1
	  if(OLDSCHOOL==1){
	      fscanf(header,INPUT1OLD,&t,&dumi[0]);
	  }
	  else if(OLDSCHOOL==0){
	      fscanf(header,INPUT1,&t,&dumi[0],&dumi[1]);
	  }
	  
      }
      if(POSTPROC==0){
	  fprintf(log_file,"restarting at t=%15.10g\n",t);
      }
      if(dumi[0]!=0){
	  fprintf(log_file,"1.5: Warning: cannot take in interpolated data: expected: 0 got: %d\n",dumi[0]);
      }
      if(OLDSCHOOL==0){
	  if(dumi[1]!=0){
	      fprintf(log_file,"1.6: Warning, data is not gridded correctly: expected: 0 got: %d\n",dumi[1]);
	      fflush(log_file);
	  }
      }
      while(fgetc(header)!='\n'); // skip to next line
      
      ch=fgetc(header);
      if(ch=='#'){
	  while(fgetc(header)!='\n'); // skip comment line
      }
      else{
	  fprintf(fail_file,"2: data file doesn't match expected format\n");
	  myexit(1);
      }
      
      if(mpicombine==1){
	  N[1]=totalsize[1];
	  N[2]=totalsize[2];
	  N[3]=totalsize[3];
      }
      else{
	  if(FULLINPUT==0){
	      N[1]=N1;
	      N[2]=N2;
	      N[3]=N3;
	  }
	  else if(FULLINPUT==1){
	      N[1]=OUTFULL1;
	      N[2]=OUTFULL2;
	      N[3]=OUTFULL3;
	  }
	  else if(FULLINPUT==2){
	      N[1]=N1M;
	      N[2]=N2M;
	      N[3]=N3M;
	  }
      }
      
      fscanf(header,"%d %d %d",&dumi[0],&dumi[1],&dumi[2]);
      if( (dumi[0]!=N[1])||(dumi[1]!=N[2])||(dumi[2]!=N[3])){
	  fprintf(fail_file,"2.4: wrong data size: got: %d %d %d, expected %d %d %d\n",dumi[0],dumi[1],dumi[2],N[1],N[2],N[3]);
	  myexit(1);
      }
      
      
      while(fgetc(header)!='\n'); // skip to next line
      
      
      if(which==AVG2DTYPE){
	  // get avg type header data
	  for(i=1;i<=16;i++) fgetc(header); // skip text
	  fscanf(header,INPUTAVGH1,&tavgstart,&tavgfinal);
	  while(fgetc(header)!='\n'); // skip to next line
	  for(i=1;i<=12;i++) fgetc(header); // skip text
	  fscanf(header,"%d",&avgcount);
	  while(fgetc(header)!='\n'); // skip to next line
      }
      
      ch=fgetc(header);
      if(ch=='#'){
	  while(fgetc(header)!='\n'); // skip comment line
      }
      else{
	  fprintf(fail_file,"2.5: data file doesn't match expected format\n");
	  myexit(1);
      }
      if((mpicombine==1)||((mpicombine==0)&&(binaryoutput==1))){
        fclose(header);
      }

  }
#if(USEMPI)
  MPI_Bcast(&t,1,MPI_SFTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(&npoldschool,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&avg2doldschool,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&tavgstart,1,MPI_SFTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(&tavgfinal,1,MPI_SFTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(&avgcount,1,MPI_INT,0,MPI_COMM_WORLD);
  
#endif
  fprintf(log_file,"End Reading in header\n"); fflush(log_file);
  ///////////////
  // END PROCESS HEADER
  ///////////////
  

  ///////////////
  // BEGIN PROCESS DATA
  ///////////////

  fprintf(log_file,"reading data file %d\n",which); fflush(log_file);


  fprintf(log_file,"1proc=%d : god: sin=%d\n",myid,sin); fflush(log_file);
  fprintf(stderr,"1proc=%d : god: sin=%d\n",myid,sin); fflush(stderr);
#if(USEMPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if(mpicombine==1){
#if(USEMPI)
      if(USEROMIO){
#if(USEROMIO)
	  MPI_File_read_all(fh, writebuf, bufcount, MPI_FTYPE, &status);
	  MPI_File_close(&fh);
#endif
      }
      else if(USEJONIO){
	  jonio_seperate(1,MPI_FTYPE,numcolumns,sizeof(FTYPE),data,jonio,writebuf,tempbuf);
      }
#endif
    range=0;
  }
  else{
    range=FULLINPUT;
  }

  fprintf(log_file,"2proc=%d : god: sin=%d\n",myid,sin); fflush(log_file);
  fprintf(stderr,"2proc=%d : god: sin=%d\n",myid,sin); fflush(stderr);
#if(USEMPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  LOOPDIAGOUTPUT(range,range,range){
    BUFFERINIT;

    if((which==DTYPE)||(which==PDTYPE)||(which==ADTYPE)||(which==FLTYPE) ){
      if(analoutput!=13){
	for(l=1;l<=NUMSCA;l++){
	  if(mpicombine==0) fscanf(data,INPUT2,&sin[l][k][j][i]);
	  else if(mpicombine==1){ sin[l][k][j][i]=writebuf[BUFFERMAP];}
	}
      }
      for(l=1;l<=REALNUMVEC;l++){
	if(which==FLTYPE){
	  if(mpicombine==0) fscanf(data,INPUT2,&floorvar0[l][k][j][i]);
	  else if(mpicombine==1){ floorvar0[l][k][j][i]=writebuf[BUFFERMAP];}
	}
	for(m=1;m<=3;m++){
	  if(mpicombine==0){
	    fscanf(data,INPUT3I,&vin[l][m][k][j][i]);
	    if( (l==REALNUMVEC)&&(m==3)) fscanf(data,"\n");
	    else fscanf(data," ");
	  }
	  else{
	    vin[l][m][k][j][i]=writebuf[BUFFERMAP];
	  }
	}
      }
      if((analoutput==13)){
	for(l=1;l<=1;l++){
	  if(mpicombine==0) fscanf(data,INPUT2,&sin[l][k][j][i]);
	  else if(mpicombine==1){ sin[l][k][j][i]=writebuf[BUFFERMAP];}
	}
	// dummy, not used, but avoid <0 condition below
	sin[2][k][j][i]=pow(sin[1][k][j][i],gam);
	sin[3][k][j][i]=0.0;	
      }
    }
    else if(which==NPTYPE){ // sigma read(only used during post process) (don't worry about nuvisc)
      if(npoldschool<=1){	
	for(l=1;l<=3;l++){
	  for(m=1;m<=3;m++){
	    if(mpicombine==0) fscanf(data,INPUT2,&sigma[l][m][k][j][i]);
	    else if(mpicombine==1){ sigma[l][m][k][j][i]=writebuf[BUFFERMAP];}
	  }
	}
      }
      else{
	if(mpicombine==0) fscanf(data,INPUT2,&sigma[1][1][k][j][i]);
	else if(mpicombine==1){ sigma[1][1][k][j][i]=writebuf[BUFFERMAP];}
	if(mpicombine==0) fscanf(data,INPUT2,&sigma[1][2][k][j][i]);
	else if(mpicombine==1){ sigma[1][2][k][j][i]=writebuf[BUFFERMAP];}
	if(mpicombine==0) fscanf(data,INPUT2,&sigma[1][3][k][j][i]);
	else if(mpicombine==1){ sigma[1][3][k][j][i]=writebuf[BUFFERMAP];}
	if(mpicombine==0) fscanf(data,INPUT2,&sigma[2][2][k][j][i]);
	else if(mpicombine==1){ sigma[2][2][k][j][i]=writebuf[BUFFERMAP];}
	if(mpicombine==0) fscanf(data,INPUT2,&sigma[2][3][k][j][i]);
	else if(mpicombine==1){ sigma[2][3][k][j][i]=writebuf[BUFFERMAP];}
	if(mpicombine==0) fscanf(data,INPUT2,&sigma[3][3][k][j][i]);
	else if(mpicombine==1){ sigma[3][3][k][j][i]=writebuf[BUFFERMAP];}
	// for completeness
	sigma[2][1][k][j][i]=sigma[1][2][k][j][i];
	sigma[3][1][k][j][i]=sigma[1][3][k][j][i];
	sigma[3][2][k][j][i]=sigma[2][3][k][j][i];
      }
      if(npoldschool>0){
	if(mpicombine==0) fscanf(data,INPUT2B,&nu_real[k][j][i]);
	else if(mpicombine==1){ nu_real[k][j][i]=writebuf[BUFFERMAP];}
      }
    }// end if reading in sigma
    else if(which==AVG2DTYPE){ // avg2d read(only used during post process)
      if(avg2doldschool>1){
	for(l=1;l<=NUMAVG2_3;l++){
	  if(mpicombine==0) fscanf(data,INPUT7,&avg2_3[l][j][i]);
	  else if(mpicombine==1){ avg2_3[l][j][i]=writebuf[BUFFERMAP];}
	}
      }
      else if(avg2doldschool==1){
	for(l=1;l<=8;l++){
	  if(mpicombine==0) fscanf(data,INPUT7,&avg2_3[l][j][i]);
	  else if(mpicombine==1){ avg2_3[l][j][i]=writebuf[BUFFERMAP];}
	}
	for(l=9;l<=NUMAVG2_3;l++){
	  avg2_3[l][j][i]=0; // we'll compute the sigmas just before output
	}
      }
      // hack for sigma computation in diag.c from average data if avg2doldschool==1, or for other calcs
      s[1][k][j][i]=avg2_3[1][j][i];
      s[2][k][j][i]=avg2_3[2][j][i];
      v[1][1][k][j][i]=avg2_3[6][j][i];
      v[1][2][k][j][i]=avg2_3[7][j][i];
      v[1][3][k][j][i]=avg2_3[8][j][i];
    }// end if reading in avg2d

  }// end over loop of data


  fprintf(log_file,"3proc=%d : god: sin=%d\n",myid,sin); fflush(log_file);
  fprintf(stderr,"3proc=%d : god: sin=%d\n",myid,sin); fflush(stderr);
#if(USEMPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif


  if(mpicombine==0) fclose(data);
  else{
#if(USEMPI)
      if(USEROMIO){
#if(USEROMIO)
	  free(writebuf);
	  MPI_Type_free(&newtype);
#endif
      }
      else if(USEJONIO){
	  jonio_seperate(2,MPI_FTYPE,numcolumns,sizeof(FTYPE),data,jonio,writebuf,tempbuf); 
	  jonio_seperate(3,MPI_FTYPE,numcolumns,sizeof(FTYPE),0,jonio,writebuf,tempbuf); 
      }
#endif
  }

  fprintf(log_file,"4proc=%d : god: sin=%d\n",myid,sin); fflush(log_file);
  fprintf(stderr,"4proc=%d : god: sin=%d\n",myid,sin); fflush(stderr);
#if(USEMPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif


  // check inputted data for correctness
  if((which==DTYPE)||(which==PDTYPE)||(which==ADTYPE)||(which==FLTYPE) ){  
    LOOPDIAGOUTPUT(range,range,range){
      if(s[1][k][j][i]<1E-20){
	fprintf(fail_file,"input data has mass density <1E-20: %d %d %d %15.10g\n",k,j,i,s[1][k][j][i]);
	myexit(1);
      }
      if(s[2][k][j][i]<1E-20 ){
	fprintf(fail_file,"input data has ie density <1E-20: %d %d %d %15.10g\n",k,j,i,s[2][k][j][i]);
	myexit(1);
      }
    }
  }

  fprintf(log_file,"Done Reading in %s file\n",dfheader);
  fflush(log_file);

  fprintf(log_file,"end: init_rundat\n"); fflush(log_file);
  return(0);
}

// mapping from image space to function space
#if(!GAMMIEIMAGE)
#define CMAPSCA1(x) pow(10.0,x)
#define CMAPSCA2(x) pow(10.0,x)
#else
#define CMAPSCA1(x) (x)
#define CMAPSCA2(x) (x)
#endif

#define CMAPSCA3(x) (x)
#define CMAPSCAGEN(x) (x)

#define CMAPVEC1(x) (x)
#define CMAPVEC2(x) (x)

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

#if(GAMMIEIMAGE)
#define LINLOGFIX 0 // whether to use both linear and log10 input data to recreate old data
#else
#define LINLOGFIX 1
#endif
// 0: use log only
// 1: use linear and log
// 2: use linear only

#if(GAMMIEIMAGE)
#define READMM 0 // whether to read file for min/max
#else
#define READMM 1 // whether to read file for min/max
#endif
// used to input image data(much extracted from image()
// assume one doesn't care about data scaling, just used to interpolate during post process run.  Only need to invert color mapping for ppm
int init_runimage(int im_cnt,int wsca,int wvec,int call_code,int outtype)
{
  // CUT PASTE START 0 
  FILE* data;
  FILE* in[3];
  char fname[200];
  int dumi[50];
  SFTYPE dumf[50];
  char ch;
  int i,j,k,l,m,p;
  char crap[50],testcrap[50];
  FTYPE (*sin)[N3M][N2M][N1M] ;
  FTYPE (*vin)[3][N3M][N2M][N1M] ;
  char dfheader[MAXFILENAME];
  char ifheader[MAXFILENAME];
  SFTYPE liq,a,b,c,lmax,lmin,rholiq ;
  unsigned char tempuch; // int doesn't work // GODMARK
  
  int iii,dualgostart,dualgoend;
  int ll,mm ;
  int q;
  int wshort;
  SFTYPE ftempfix;
  int floop; //tells loop when first entered

  FILE *im_file;
  static FILE *iparout ;
  static char ifnam[MAXFILENAME], temp[MAXFILENAME] ;  

  static int firstfirsttime=1;
  static int lastlasttime=0;
  static int firsttimes[ITYPES][NUMSCA+1];
  static int firsttimev[ITYPES][NUMVEC+1];
  FILE*size_file;

  int im_cnts[ITYPES][NUMSCA+1];
  int im_cntv[ITYPES][NUMVEC+1];

  int noiparusefile;
  SFTYPE ftemp[2];
  FILE * iparuse;
  static int outtypein;
  char temps[MAXFILENAME];
  char command[MAXFILENAME];
  // CUT PASTE END 0
  char dumc[50];
  int realouttype;
  static int qstart;


  fprintf(log_file,"begin: init_runimage ... "); fflush(log_file);

  // CUT PASTE START 1
  // for now assume all im_cnt same
  for(i=1;i<=NUMSCA;i++){
    im_cnts[outtype][i]=im_cnt;
  }
  for(i=1;i<=NUMVEC;i++){
    im_cntv[outtype][i]=im_cnt;
  }
  // CUT PASTE END 1


  noiparusefile=0; // GODMARK

  if(myid<=0){
    if(READMM){
      sprintf(temps,"%s%s",DATADIR,IMAGEDIR);
      sprintf(ifnam,"%s0_iparuse%s",temps,PAR2EXT);
      if((POSTPROC==1)||((POSTPROC==0)&&(runtype>0)&&(directinput>0))){ // then should have no completed use file, no matter what time, and want to use it for seemless continuation of image if not postproc.
	// if pp, then only have file if time is later than timagescale.  rest is unreliable for computation no matter what since dynamic anyways.  assumes timagescale is set correctly........ Just assume file exists if pp for now
	fprintf(logfull_file,"#reading iparuse file(init.c)\n");
	fflush(logfull_file);
	if((iparuse = fopen(ifnam,"r"))==NULL) {
	  fprintf(fail_file,"error reading iparuse output file %s\n",ifnam) ;
	  noiparusefile=1;
	}
	if(!noiparusefile){
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
		  }
		}
	      }// end over components(q)
	      // q==0
	      for(i=0;i<=1;i++){ // min and max
		for(m=0;m<CTYPES;m++){
		  // really just a hack, should read in or set  q=0
		  if(i==1){ // certainly the maximum possible
		    mmv[outtypein][m][ll][0][i]=mmv[outtypein][m][ll][1][i]*mmv[outtypein][m][ll][1][i]+mmv[outtypein][m][ll][2][i]*mmv[outtypein][m][ll][2][i]+mmv[outtypein][m][ll][3][i]*mmv[outtypein][m][ll][3][i];
		  }
		  if(i==0){
		    mmv[outtypein][m][ll][0][i]=0; // certainly minimum possible
		  }
		}
	      }
	    }// over ll
	  } // over outtype
	  fclose(iparuse);
	  fprintf(logfull_file,"#done reading iparuse file(init.c)\n");
	  fflush(logfull_file);
	} // end if read paruse file
      }
    }
  }
#if(USEMPI)
  // tell other cpus if file existed so can stop if necessary
  MPI_Bcast(&noiparusefile,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
  if(noiparusefile){
    fprintf(fail_file,"No iparuse file as needed in init_runimage\n");
    myexit(7);
  }

  ///////////  SCALARS
  if(wsca!=0){
    for(l=1;l<=NUMSCA;l++){
      /* if not to do all, pick */
      if(wsca!=-1){ ll=wsca; }
      else ll=l;

      fprintf(log_file,"reading image file %d\n",ll); fflush(log_file);
      
      // assume samplei originally 0, no interp from interp'ed data

      if(LINLOGFIX==0){ dualgostart=0; dualgoend=0;}
      if(LINLOGFIX==1){ dualgostart=0; dualgoend=1;}
      if(LINLOGFIX==2){ dualgostart=1; dualgoend=1;}

      for(iii=dualgostart;iii<=dualgoend;iii++){ // assume only the 0th tile(sliceloop)
	// setup file input
	sprintf(temps,"%s%s",DATADIR,IMAGEDIR);
	if(!GAMMIEIMAGE){
	  strcpy(ifheader,"imx");
	  if(!OLDIMAGEFORMAT){
	    sprintf(temps,"%simx%01d-%01d-0-s%01d/",temps,outtype,iii,ll);
	    sprintf(ifnam,"%s%s%01d-%01d-0-s%01d-%04d%s%s",temps,ifheader,outtype,iii,ll,im_cnts[outtype][ll],IMGEXT,myidtxt);
	  }
	  else{
	    sprintf(temps,"%simx%01d-%01d-s%01d/",temps,outtype,iii,ll);
	    sprintf(ifnam,"%s%s%01d-%01d-s%01d-%04d%s%s",temps,ifheader,outtype,iii,ll,im_cnts[outtype][ll],IMGEXT,myidtxt);
	  }
	}
	else{
	  sprintf(ifnam,"%s%sim%01dp%04d",DATADIR,IMAGEDIR,ll-1,im_cnts[outtype][ll]); // gets RHO and UU
	}
	if(IMAGEFORMATINPUT==0){
	  strcat(ifnam,".r8");
	  if(GZIPIMAGEINPUT==0){
	    im_file = fopen(ifnam,"rb");
	  }
	  else if(GZIPIMAGEINPUT>0){
	    strcat(ifnam,".gz");
	    
	    sprintf(temp,"gzip -d < %s",ifnam);
	    im_file = popen(temp,"r");
	    //fprintf(stdout,"%s %s\n",ifnam,temp);
	  }
	  if(im_file==NULL){
	    fprintf(fail_file,"error opening image file: %s\n",ifnam) ;
	    myexit(2) ;
	  }
	  // skip 4 comment lines
	  if(OLDSCHOOL2==0){
	    while(fgetc(im_file)!='\n'); // skip to next line
	    while(fgetc(im_file)!='\n'); // skip to next line
	    while(fgetc(im_file)!='\n'); // skip to next line
	    while(fgetc(im_file)!='\n'); // skip to next line
	  }
	}
	if(IMAGEFORMATINPUT==1){
	  fprintf(fail_file,"Can't input this format of image, since not simple to reverse lookup function value to interpolate\n");
	  myexit(1);
	}
	
	// now set function map using min/max
	// same code as in image()
	if(ll==1){
	  if(iii==0){
	    b=FMAPSCA1(mms[outtype][iii][ll][0]);
	    a=255./(FMAPSCA1(mms[outtype][iii][ll][1])-b);
	  }
	  if(iii==1){
	    b=FMAPSCAGEN(mms[outtype][iii][ll][0]);
	    a=255./(FMAPSCAGEN(mms[outtype][iii][ll][1])-b);
	  }
	}
	else if(ll==2){
	  if(iii==0){
	    b=FMAPSCA2(mms[outtype][iii][ll][0]);
	    a=255./(FMAPSCA2(mms[outtype][iii][ll][1])-b);
	  }
	  if(iii==1){
	    b=FMAPSCAGEN(mms[outtype][iii][ll][0]);
	    a=255./(FMAPSCAGEN(mms[outtype][iii][ll][1])-b);
	  }
	}
	else if(ll==3){
	  b=FMAPSCA3(mms[outtype][iii][ll][0]);
	  a=255./(FMAPSCA3(mms[outtype][iii][ll][1])-b);
	}
	c=-a*b;

	//fprintf(stdout,"%d %d %d %15.10g %15.10g\n",outtype,iii,ll,mms[outtype][iii][ll][0],mms[outtype][iii][ll][1]);
	
	LOOPINI{
	  /* read value */
	  //	  fread(&tempuch,sizeof(unsigned char),1,im_file);
	  fread(&tempuch,sizeof(unsigned char),1,im_file);
	  //tempuch=(unsigned char)fgetc(im_file);
	  liq=((SFTYPE)(tempuch))/a+b; // inverse value
	  //liq=(FTYPE)tempuch;
	  //if((ll==1)&&(iii==0)){
	  //fprintf(stdout,"%d %d %15.10g %u %15.10g %15.10g %15.10g\n",j,i,liq,tempuch,a,b,c);
	  //} 
	  // now invert from color to function space
	  if(ll==1){
	    if(iii==0){
	      ftempfix=CMAPSCA1(liq);
	    }
	    if(iii==1){
	      ftempfix=CMAPSCAGEN(liq);
	    }
	  }
	  else if(ll==2){
	    if(iii==0){
	      ftempfix=CMAPSCA2(liq);
	    }
	    if(iii==1){
	      ftempfix=CMAPSCAGEN(liq);
	    }
	  }
	  else if(ll==3){
	    ftempfix=CMAPSCA3(liq);
	  }
	  liq=ftempfix;	

	  if(LINLOGFIX==1){
	    // only modify in a way which uses the linear and log10 data best
	    if((iii==0)&&(liq<mms[outtype][iii][ll][1]/(256.0))){
	      s[ll][0][j][i]=liq; // assume interp only uses k=0
	    }
	    if((iii==1)&&(liq>=mms[outtype][iii][ll][1]/(256.0))){
	      s[ll][0][j][i]=liq; // assume interp only uses k=0
	    }
	  }
	  else if(LINLOGFIX==0){
	    if(iii==0){
	      s[ll][0][j][i]=liq;
	    }
	  }
	  else if(LINLOGFIX==2){
	    if(iii==1){
	      s[ll][0][j][i]=liq;
	    }
	  }
	  // GODMARK
	  //fprintf(stdout,"%d %d %15.10g %15.10g %15.10g %15.10g\n",j,i,(SFTYPE)(tempuch),liq,ftempfix,s[ll][0][j][i]);
	}


	if(GZIPIMAGEINPUT==0){
	  fclose(im_file) ;
	}
	else if(GZIPIMAGEINPUT>0){
	  pclose(im_file);
	}
      }// over dualgo
      /* cut short loop if only to do one */
      if(wsca!=-1) l=NUMSCA;
      if(GAMMIEIMAGE&&(ll==2)) break; // stop at UU
    }// over scalars
  }// if any scalars


  //fflush(stdout); myexit(0);

  //////////// VECTORS
  if(wvec!=0){
    for(l=1;l<=REALNUMVEC;l++){
      /* if not to do all, pick */
      if(wvec!=-1){ ll=wvec; }
      else ll=l;
      

      for(q=1;q<=3;q++){ // over components (no need to do q=0)
	fprintf(log_file,"reading image file %d %d\n",ll,q); fflush(log_file);

	dualgoend=0; // force for input, no need for other data since can just deduce this from linear rho/v in pp
	for(iii=0;iii<=dualgoend;iii++){
	  // setup file output
	  if(!GAMMIEIMAGE){
	    sprintf(temps,"%s%s",DATADIR,IMAGEDIR);	    
	    
	    strcpy(ifheader,"imx");

	    if(!OLDIMAGEFORMAT){
	      sprintf(temps,"%simx%01d-%01d-0-v%01d-%01d/",temps,outtype,iii,ll,q);
	      sprintf(ifnam,"%s%s%01d-%01d-0-v%01d-%01d-%04d%s%s",temps,ifheader,outtype,iii,ll,q,im_cnts[outtype][ll],IMGEXT,myidtxt);
	    }
	    else{
	      sprintf(temps,"%simx%01d-%01d-v%01d-%01d/",temps,outtype,iii,ll,q);
	      sprintf(ifnam,"%s%s%01d-%01d-v%01d-%01d-%04d%s%s",temps,ifheader,outtype,iii,ll,q,im_cnts[outtype][ll],IMGEXT,myidtxt);
	    }
	  }
	  else{
	    sprintf(ifnam,"%s%sim%01dp%04d",DATADIR,IMAGEDIR,(ll-1)*3+(q-1)+2,im_cnts[outtype][ll]); // gets VX,VY,VZ,BX,BY,BZ
	  }

	  if(IMAGEFORMATINPUT==0){
	    strcat(ifnam,".r8");
	    if(GZIPIMAGEINPUT==0){
	      im_file = fopen(ifnam,"rb");
	    }
	    else if(GZIPIMAGEINPUT>0){
	      strcat(ifnam,".gz");
	      sprintf(temp,"gzip -d < %s",ifnam);
	      strcpy(ifnam,temp); // for below fprintf
	      im_file = popen(ifnam,"r");
	    }
	    if(im_file==NULL){
	      fprintf(fail_file,"error opening image file: %s\n",ifnam) ;
	      myexit(2) ;
	    }
	    // skip 4 comment lines
	    if(OLDSCHOOL2==0){
	      while(fgetc(im_file)!='\n'); // skip to next line
	      while(fgetc(im_file)!='\n'); // skip to next line
	      while(fgetc(im_file)!='\n'); // skip to next line
	      while(fgetc(im_file)!='\n'); // skip to next line
	    }
	  }
	  if(IMAGEFORMATINPUT==1){
	    fprintf(fail_file,"Can't input this format of image, since not simple to reverse lookup function value to interpolate\n");
	    myexit(1);
	  }

	  // setup map
	  // same code as in image()
	  if(ll==1){
	    if(iii==0){
	      // warning:  silent fail if 0,1 are both 0 for q=0 ?
	      b=FMAPVEC1(mmv[outtype][iii][ll][q][0]);
	      a=255./(FMAPVEC1(mmv[outtype][iii][ll][q][1])-b);
	    }
	    if(iii==1){ // rho*v // only good if v range is - to + values
	      b=FMAPVEC1(mmv[outtype][iii][ll][q][0]);
	      a=255./(FMAPVEC1(mmv[outtype][iii][ll][q][1])-b);
	    }
	  }
	  if(ll==2){
	    if(iii==0){// B
	      b=FMAPVEC2(mmv[outtype][iii][ll][q][0]);
	      a=255./(FMAPVEC2(mmv[outtype][iii][ll][q][1])-b);
	    }
	    if(iii==1){ // B*v // only good if v range is - to + values
	      b=FMAPVEC2(mmv[outtype][iii][ll][q][0]);
	      a=255./(FMAPVEC2(mmv[outtype][iii][ll][q][1])-b);
	    }
	  }
	  c=-a*b;

	  // input image without map
	  LOOPINI{
	    
	    fread(&tempuch,sizeof(unsigned char),1,im_file);
	    liq=(SFTYPE)((int)(tempuch))/a+b; // inverse value

	    // invert
	    if(ll==1){
	      ftempfix=CMAPVEC1(liq);
	    }
	    else if(ll==2){
	      ftempfix=CMAPVEC2(liq);
	    }
	    liq=ftempfix;

	    v[ll][q][0][j][i]=liq; // assume interp only uses k=0
	  }
      
	      
	  if(GZIPIMAGEINPUT==0){
	    fclose(im_file) ;
	  }
	  else if(GZIPIMAGEINPUT>0){
	    pclose(im_file);
	  }
	}//dualgo
      } // over components of the vector
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMVEC;
    }// over vectors
  }// if any vectors

#if(0)
  LOOPF{
    fprintf(stdout,"%d %d %d %15.10g \n",k,j,i,s[1][k][j][i]);
  }
  fflush(stdout);
#endif

  fprintf(log_file,"end: init_runimage\n"); fflush(log_file);
  return(0);
}



// which: -1: all 0: params 1: grid no interp 2: grid interp 3: interp params
int init_outgparm(int which)
{
  int error;
  int i,j,k,l,m;
  FILE* in[4];
  FILE* out;
  char filename[MAXFILENAME];
  char temps[MAXFILENAME];
  int realdumn1,realdumn2;
  // for interp grid output(like in numerics.c's interpolate)
  SFTYPE n1,n2,N[3+1];
  int outtype;
  SFTYPE width[2];
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
  char truemyidtxt[100];
  int pdiddy,pdiddystart,pdiddyto;
  int nextbuf;
  int range;
  FILE * headerptr;
  int dumi;
  float dumf;
  double dumlf;


  fprintf(log_file,"begin: init_outgparm ... "); fflush(log_file);


  if(DOGPARDIAG){
  //output global parameters
  if(myid<=0){
    fprintf(logfull_file,"outputting global parameters\n");
    fflush(logfull_file);
    if( (which==-1)||(which==0)){
      sprintf(filename,"%s0_gparam%s",DATADIR,PAR2EXT);
      
      if((out=fopen(filename,"wt"))==NULL){
	fprintf(fail_file,"Can't open %s for writting\n",filename);
	myexit(1);
      }
      else{
	fprintf(logfull_file,"outputting parameter file\n"); fflush(log_file);
	fprintf(out,HEADER3_P
		,"PVER","TYPE"
                ,PVER,PTYPE
		,"N1","N2","N3"
		,totalsize[1],totalsize[2],totalsize[3]
		,"x1-start","x2-start","x3-start","x1-end","x2-end","x3-end"
		,x1in,x2in,x3in,x1out,x2out,x3out
		,"rg","rgp","cs","coolfact","gam","alpha_real0","n_real"
		,rg,rgp,cs,coolfact,gam,alpha_real0,n_real
		,"nu_vnr","nu_l","nu_ten","cour","cour2"
		,nu_vnr,nu_l,nu_ten,cour,cour2
		,"GRAVC","MASSBH","invsol2","blackholejz"
		,GRAVC,MASSBH,invsol2,blackholejz
		,"tstart","tf","tavgi","tavgf","numavg","timagescale"
		,tstart,tf,tavgi,tavgf,numavg,timagescale
		,"DTl","DTd","DTi","DTloss","DTfloor","DTtimestep","DTpd","DTener","DTfld","DTmode"
		,DTl,DTd,DTi,DTloss,DTfloor,DTtimestep,DTpd,DTener,DTfld,DTmode
		,"resist_real0","nu_sh"
		,resist_real0,nu_sh
		,"vg-x1","vg-x2","vg-x3"
		,vg[1],vg[2],vg[3]
		,"coord","fullvec","analoutput","DYNAMICMM"
		,COORD,REALNUMVEC==NUMVEC,analoutput,DYNAMICMM
		,"trans","transx1","transrhox1","transiex1","transv1x1","transv2x1","transv3x1","transmagx1","transx2","transrhox2","transiex2","transv1x2","transv2x2","transv3x2","transmagx2","transx3","transrhox3","transiex3","transv1x3","transv2x3","transv3x3"
		,trans,transx1,transrhox1,transiex1,transv1x1,transv2x1,transv3x1,transmagx1,transx2,transrhox2,transiex2,transv1x2,transv2x2,transv3x2,transmagx2,transx3,transrhox3,transiex3,transv1x3,transv2x3,transv3x3
		,"press","pressx1","pressx2","pressx3"
		,press,pressx1,pressx2,pressx3
		,"mag","transbzx","transbzy","stepmocct","mocctvx1","mocctvx2","mocctvx3","mocctbx1","mocctbx2","mocctbx3"
		,mag,transbzx,transbzy,stepmocct,mocctvx1,mocctvx2,mocctvx3,mocctbx1,mocctbx2,mocctbx3
		,"ie","visc_art","visc_real","vreal","vischeat","mdotin","cool","res_resl","rreal","resheat","advint","kever"
		,ie,visc_art,visc_real,vreal,vischeat,mdotin,cool,res_real,rreal,resheat,advint,kever
		,"intix1","intox1","intix2","intox2","intix3","intox3"
		,intix1,intox1,intix2,intox2,intix3,intox3
		,"nonunigridx1","nonunigridx2","nonunigridx3","simplebc","bcix1","bcox1","bcix2","bcox2","bcix3","bcox3"
		,nonunigridx1,nonunigridx2,nonunigridx3,simplebc,bcix1,bcox1,bcix2,bcox2,bcix3,bcox3
		);
      }
      
      fclose(out);
    }
  }
  }
  if(mpicombine==0){
    strcpy(truemyidtxt,myidtxt);
  }
  else strcpy(truemyidtxt,"");
  

  if(DOPARDIAG){
    fprintf(log_file,"output grid data\n"); fflush(log_file);
    // always need general grid data
    //output dx x diffs to file--active and ghost zones
    if( (which==-1)||(which==1) ){

      numcolumns=21+NUMSCA*3+NUMVEC*3;
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
	
	MPI_Type_create_darray(numprocs, myid, ndims, array_of_gsizes, 
			       array_of_distribs, array_of_dargs,
			       array_of_psizes, order, MPI_FTYPE, &newtype);
	MPI_Type_commit(&newtype);
	MPI_Type_size(newtype, &bufcount);
	bufcount = bufcount/sizeof(FTYPE);
	writebuf = (FTYPE *) malloc(bufcount * sizeof(FTYPE));
#endif
	  }
	  else if(USEJONIO){
	      if(sizeof(FTYPE)==sizeof(double)) jonio_init_mem(numcolumns,sizeof(FTYPE),0,0,&jonio,0,0,&writebuf,0,0,&tempbuf);
	      else if(sizeof(FTYPE)==sizeof(float)) jonio_init_mem(numcolumns,sizeof(FTYPE),0,&jonio,0,0,&writebuf,0,0,&tempbuf,0);
	      else if(sizeof(FTYPE)==sizeof(unsigned char)) jonio_init_mem(numcolumns,sizeof(FTYPE),&jonio,0,0,&writebuf,0,0,&tempbuf,0,0);
	  }
#endif
      }
      


      if(mpicombine==0){
	if(FULLOUTPUT==0){ // only need if not doing fulloutput
	  pdiddyto=2; // do LOOP output
	}
	else pdiddyto=1; // only do LOOPF output
	pdiddystart=1;
      }
      else{
	// only do LOOP output
	pdiddystart=2;
	pdiddyto=2;
      }

      for(pdiddy=pdiddystart;pdiddy<=pdiddyto;pdiddy++){ // loop over type of output (loop or loopf)
	if(mpicombine==0){
	  if(pdiddy==1) range=2;
	  else range=0;
	}
	else if(mpicombine==1){ range=0;}


	for(l=1;l<=NUMGRID;l++){ // loop over grids (a/b)
	  fprintf(log_file,"outputting grid file: act/full: %d grid a/b: %d\n",pdiddy,l); fflush(log_file);

	  if(pdiddy==1) sprintf(temps,"%s0_grid",DATADIR);
	  else sprintf(temps,"%s0_gridact",DATADIR);
	  sprintf(filename,"%s%01d%s%s",temps,l,PAREXT,truemyidtxt);
	  
	  
	  if(mpicombine==0){
	      if(binaryoutput==0){
		  if((out=fopen(filename,"wt"))==NULL){
		      fprintf(fail_file,"Can't open %s for writting\n",filename);
		      myexit(1);
		  }
		  headerptr=out;
	      }
	      else if(binaryoutput==1){
		  if((out=fopen(filename,"wb"))==NULL){
		      fprintf(fail_file,"Can't open %s for writting\n",filename);
		      myexit(1);
		  }
	      }
	  }
	  else if(mpicombine==1){
#if(USEMPI)
	      if(USEROMIO){
#if(USEROMIO)
		  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, 
				MPI_INFO_NULL, &fh);
		  MPI_File_set_view(fh, 0, MPI_FTYPE, newtype, "native", MPI_INFO_NULL);
		  // all that needs to be done now is fill writebuf with the data
#endif
	      }
	      else if(USEJONIO){
		  jonio_init_fp(&out,1,filename);
	      }
#endif
	  }
	  fprintf(log_file,"."); fflush(log_file);

	  if( ((myid==0)&&(mpicombine==1))||((binaryoutput==1)&&(mpicombine==0)) ){
	    // write seperate header file header
	    strcat(filename,".head");
	    if((headerptr=fopen(filename,"wt"))==NULL){
	      fprintf(fail_file,"Can't open %s for writting\n",filename);
	      myexit(1);
	    }
	  }
	  fprintf(log_file,"."); fflush(log_file);
	  if( ((myid==0)&&(mpicombine==1))||(mpicombine==0)){
	    fprintf(headerptr,"#%10s\n%10d %10d\n","GRIDVER",GRIDVER,GRIDTYPE);
	    fprintf(headerptr,"#%4s %4s %4s %4s "
		    "%21s %21s %21s "
		    "%21s %21s %21s "
		    "%21s %21s %21s %21s %21s %21s %21s %21s "
		    "%21s %21s %21s "
		    ,"grid","k","j","i","dx1","dx2","dx3","x1","x2","x3","g1","dg1","g2","dg2","g3","dg3","g4","dg4","dvl1","dvl2","dvl3");
	    for(m=1;m<=NUMSCA;m++){
	      fprintf(headerptr,"%6s%02d %6s%02d %6s%02d ","bsty",m,"bsdm",m,"bsdr",m);
	    }
	    for(m=1;m<=NUMVEC;m++){
	      fprintf(headerptr,"%6s%02d %6s%02d %6s%02d","bvty",m,"bvdm",m,"bvdr",m);
	      if(m==NUMVEC) fprintf(headerptr,"\n"); else fprintf(headerptr," ");
	    }
	  }
	  fprintf(log_file,"|"); fflush(log_file);
	  if(mpicombine==0){
	    LOOPDIAGOUTPUT(range,range,range){
	      if(binaryoutput==0){
		fprintf(out,HEADER4_P,
			l,startpos[3]+k,startpos[2]+j,startpos[1]+i,
			dx[l][1][i],dx[l][2][j],dx[l][3][k],
			x[l][1][i],x[l][2][j],x[l][3][k],
#if(LOWMEMMODE==0)
			g[l][1][i],dg[l][1][i],g[l][2][i],dg[l][2][i],g[l][3][i],dg[l][3][i],g[l][4][j],dg[l][4][j],
			dvl[l][1][i],dvl[l][2][j],dvl[l][3][k]
#else
			0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			0.0,0.0,0.0
#endif
			);
		if((PUREBC==0)&&(BOUNDTYPE==1)){
		  for(m=1;m<=NUMSCA;m++){
		    fprintf(out,"%6d ",bcs[m][1][k][j][i]);
		    fprintf(out,"%6d ",bcs[m][2][k][j][i]);
		    fprintf(out,"%6d ",bcs[m][3][k][j][i]);
		  }
		  for(m=1;m<=NUMVEC;m++){
		    fprintf(out,"%6d ",bcv[m][1][k][j][i]);
		    fprintf(out,"%6d ",bcv[m][2][k][j][i]);
		    fprintf(out,"%6d",bcv[m][3][k][j][i]);
		  }
		}
		else if((BOUNDTYPE==2)||(BOUNDTYPE==3)){
		  fprintf(out,"%6d ",bzmask[k][j][i]);
		  fprintf(out,"%6d ",bzmaskin[k][j][i]);
		  //	      fprintf(out,"%6d ",maskv2[k][j][i]);
		  //fprintf(out,"%6d ",maskv3[k][j][i]);
		  //	      fprintf(out,"%6d ",bzs[1][k][j][i]);
		  //fprintf(out,"%6d ",bzs[2][k][j][i]);
		  //fprintf(out,"%6d ",bzs[3][k][j][i]);
		  for(m=1;m<=14-1;m++){ // GODMARK for now
		    fprintf(out,"%6d ",0);
		  }
		}
		else{
		  for(m=1;m<=14-1+2;m++){ // GODMARK for now
		    fprintf(out,"%6d ",0);
		  }
		}
		fprintf(out,"\n");
	      }// end text output
	      else{
#if(FLOATTYPE==0)
#define WHICHDUM (dumf)
#elif(FLOATTYPE==1)
#define WHICHDUM (dumlf)
#endif
		dumi=l;
		fwrite(&dumi,sizeof(int),1,out);
		dumi=startpos[3]+k;
		fwrite(&dumi,sizeof(int),1,out);
		dumi=startpos[2]+j;
		fwrite(&dumi,sizeof(int),1,out);
		dumi=startpos[1]+i;
		fwrite(&dumi,sizeof(int),1,out);
		WHICHDUM=dx[l][1][i];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=dx[l][2][j];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=dx[l][3][k];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=x[l][1][i];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=x[l][2][j];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=x[l][3][k];
#if(LOWMEMMODE==0)
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=g[l][1][i];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=dg[l][1][i];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=g[l][2][i];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=dg[l][2][i];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=g[l][3][i];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=dg[l][3][i];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=g[l][4][j];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=dg[l][4][j];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=dvl[l][1][i];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=dvl[l][2][j];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		WHICHDUM=dvl[l][3][k];
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
#else
		WHICHDUM=0;
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
		fwrite(&WHICHDUM,sizeof(FTYPE),1,out);
#endif
		if((PUREBC==0)&&(BOUNDTYPE==1)){
		  for(m=1;m<=NUMSCA;m++){
		    dumi=bcs[m][1][k][j][i]; fwrite(&dumi,sizeof(int),1,out);
		    dumi=bcs[m][2][k][j][i];fwrite(&dumi,sizeof(int),1,out);
		    dumi=bcs[m][3][k][j][i];fwrite(&dumi,sizeof(int),1,out);
		  }
		  for(m=1;m<=NUMVEC;m++){
		    dumi=bcv[m][1][k][j][i]; fwrite(&dumi,sizeof(int),1,out);
		    dumi=bcv[m][2][k][j][i]; fwrite(&dumi,sizeof(int),1,out);
		    dumi=bcv[m][3][k][j][i]; fwrite(&dumi,sizeof(int),1,out);
		  }
		}
		else if((BOUNDTYPE==2)||(BOUNDTYPE==3)){
		  dumi=bzmask[k][j][i]; fwrite(&dumi,sizeof(int),1,out);
		  dumi=bzmaskin[k][j][i]; fwrite(&dumi,sizeof(int),1,out);
		  //	      fprintf(out,"%6d ",maskv2[k][j][i]);
		  //fprintf(out,"%6d ",maskv3[k][j][i]);
		  //	      fprintf(out,"%6d ",bzs[1][k][j][i]);
		  //fprintf(out,"%6d ",bzs[2][k][j][i]);
		  //fprintf(out,"%6d ",bzs[3][k][j][i]);
		  for(m=1;m<=14-1;m++){ // GODMARK for now
		    dumi=0; fwrite(&dumi,sizeof(int),1,out);
		  }
		}
		else{
		  for(m=1;m<=14-1+2;m++){ // GODMARK for now
		    dumi=0; fwrite(&dumi,sizeof(int),1,out);
		  }
		}
	      }// end binary output
	    }	    
	    fclose(out);
	  }// end if mpicombine==0
	  else{
	    LOOP{
	      BUFFERINIT;
	      writebuf[BUFFERMAP]=l;
	      writebuf[BUFFERMAP]=startpos[3]+k;
	      writebuf[BUFFERMAP]=startpos[2]+j;
	      writebuf[BUFFERMAP]=startpos[1]+i;
	      writebuf[BUFFERMAP]=dx[l][1][i];
	      writebuf[BUFFERMAP]=dx[l][2][j];
	      writebuf[BUFFERMAP]=dx[l][3][k];
	      writebuf[BUFFERMAP]=x[l][1][i];
	      writebuf[BUFFERMAP]=x[l][2][j];
	      writebuf[BUFFERMAP]=x[l][3][k];
#if(LOWMEMMODE==0)
	      writebuf[BUFFERMAP]=g[l][1][i];
	      writebuf[BUFFERMAP]=dg[l][1][i];
	      writebuf[BUFFERMAP]=g[l][2][i];
	      writebuf[BUFFERMAP]=dg[l][2][i];
	      writebuf[BUFFERMAP]=g[l][3][i];
	      writebuf[BUFFERMAP]=dg[l][3][i];
	      writebuf[BUFFERMAP]=g[l][4][j];
	      writebuf[BUFFERMAP]=dg[l][4][j];
	      writebuf[BUFFERMAP]=dvl[l][1][i];
	      writebuf[BUFFERMAP]=dvl[l][2][j];
	      writebuf[BUFFERMAP]=dvl[l][3][k];
#else
	      writebuf[BUFFERMAP]=0;
	      writebuf[BUFFERMAP]=0;
	      writebuf[BUFFERMAP]=0;
	      writebuf[BUFFERMAP]=0;
	      writebuf[BUFFERMAP]=0;
	      writebuf[BUFFERMAP]=0;
	      writebuf[BUFFERMAP]=0;
	      writebuf[BUFFERMAP]=0;
	      writebuf[BUFFERMAP]=0;
	      writebuf[BUFFERMAP]=0;
	      writebuf[BUFFERMAP]=0;
#endif
	      // 21
	      if((PUREBC==0)&&(BOUNDTYPE==1)){
		for(m=1;m<=NUMSCA;m++){
		  writebuf[BUFFERMAP]=bcs[m][1][k][j][i];
		  writebuf[BUFFERMAP]=bcs[m][2][k][j][i];
		  writebuf[BUFFERMAP]=bcs[m][3][k][j][i];
		}
		// 30
		for(m=1;m<=NUMVEC;m++){
		  writebuf[BUFFERMAP]=bcv[m][1][k][j][i];
		  writebuf[BUFFERMAP]=bcv[m][2][k][j][i];
		  writebuf[BUFFERMAP]=bcv[m][3][k][j][i];
		}
		// 36
	      }
	      else if((BOUNDTYPE==2)||(BOUNDTYPE==3)){
		writebuf[BUFFERMAP]=bzmask[k][j][i];
		writebuf[BUFFERMAP]=bzmaskin[k][j][i];
		for(m=1;m<=14-1;m++){ // GODMARK for now
		  writebuf[BUFFERMAP]=0;
		}
	      }	      
	      else{
		for(m=1;m<=14-1+2;m++){ // GODMARK for now
		  writebuf[BUFFERMAP]=0;
		}
	      }	      
	    }// end over loop
	    fprintf(log_file,"end over loop\n"); fflush(log_file);
#if(USEMPI)
	    if(USEROMIO){
#if(USEROMIO)
		// now write the buffer:
		MPI_File_write_all(fh, writebuf, bufcount, MPI_FTYPE, &status);
		MPI_File_close(&fh);
#endif
	    }
	    else if(USEJONIO){
		jonio_combine(1,MPI_FTYPE,numcolumns,sizeof(FTYPE),out,jonio,writebuf,tempbuf);
		jonio_combine(2,MPI_FTYPE,numcolumns,sizeof(FTYPE),out,jonio,writebuf,tempbuf);
	    }
#endif
	  }// end if mpicombine==1
	  if( ((myid==0)&&(mpicombine==1))||((mpicombine==0)&&(binaryoutput==1)) ) fclose(headerptr);
	}// end over grids (l)
      }// end pdiddy
      fprintf(log_file,"X"); fflush(log_file);
      if(mpicombine==1){
#if(USEMPI)
	  if(USEROMIO){
#if(USEROMIO)
	      free(writebuf);
	      MPI_Type_free(&newtype);
#endif
	  }
	  else if(USEJONIO){
	      jonio_combine(3,MPI_FTYPE,numcolumns,sizeof(FTYPE),0,jonio,writebuf,tempbuf);
	  }
#endif
      }
    }// end if which=-1,1 to do this grid output
  }// end if wanting grid output




  if(DOIPARDIAG||(POSTPROC==1)){
    // no MPI with interp

    // not much diff from above, just:
    
    if( (which==-1)||(which==2)){
      
      for(l=1;l<=2;l++){
	fprintf(log_file,"outputting interp grid file %d\n",l); fflush(log_file);
	sprintf(temps,"%s0_igridact",DATADIR);
	
	sprintf(filename,"%s%01d%s%s",temps,l,PAREXT,myidtxt);
	
	if((out=fopen(filename,"wt"))==NULL){
	  fprintf(fail_file,"Can't open %s for writting\n",filename);
	  myexit(1);
	}
	else{
	  fprintf(out,"#%10s\n%10d %10d\n","GRIDVER",GRIDVER,GRIDTYPE);
	  fprintf(out,"#%4s %4s %4s %4s "
		  "%21s %21s %21s "
		  "%21s %21s %21s "
		  "%21s %21s %21s %21s %21s %21s %21s %21s "
		  "%21s %21s %21s "
		  ,"grid","k","j","i","dx1","dx2","dx3","x1","x2","x3","g1","dg1","g2","dg2","g3","dg3","g4","dg4","dvl1","dvl2","dvl3");
	  for(m=1;m<=NUMSCA;m++){
	    fprintf(out,"%6s%02d %6s%02d %6s%02d ","bsty",m,"bsdm",m,"bsdr",m);
	  }
	  for(m=1;m<=NUMVEC;m++){
	    fprintf(out,"%6s%02d %6s%02d %6s%02d","bvty",m,"bvdm",m,"bvdr",m);
	    if(m==NUMVEC) fprintf(out,"\n"); else fprintf(out," ");
	  }
	  
	  k=0;
	  LOOPD{
	    fprintf(out,HEADER4_P,
		    l,k,j,i,
		    idx[l][1][i],idx[l][2][j],idx[l][3][k],
		    ix[l][1][i],ix[l][2][j],ix[l][3][k],
		    1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,
		    1.0,1.0,1.0
		    );
	    for(m=1;m<=NUMSCA;m++){
	      fprintf(out,"%6d ",0);
	      fprintf(out,"%6d ",0);
	      fprintf(out,"%6d ",0);
	    }
	    for(m=1;m<=NUMVEC;m++){
	      fprintf(out,"%6d ",0);
	      fprintf(out,"%6d ",0);
	      fprintf(out,"%6d",0);
	      if(m==NUMVEC) fprintf(out,"\n"); else fprintf(out," ");
	    }
	  }
	}
	fclose(out);
      }
    }

    

    if( (which==-1)||(which==3)){
      //output global interp parameters(using data from above)
      if(myid<=0){
	sprintf(filename,"%s0_igparam%s",DATADIR,PAR2EXT);
	
	if((out=fopen(filename,"wt"))==NULL){
	  fprintf(fail_file,"Can't open %s for writting\n",filename);
	  myexit(1);
	}
	else{
	  fprintf(log_file,"outputting interp par file\n"); fflush(log_file);
	  fprintf(out,HEADER3_P
		  ,"PVER","TYPE"
		  ,PVER,PTYPE
		  ,"N1","N2","N3"
		  ,(int)(DUMN1),(int)(DUMN2),1
		  ,"x1-start","x2-start","x3-start","x1-end","x2-end","x3-end"
		  ,ix1in,ix2in,ix3in,ix1out,ix2out,ix3out
		  ,"rg","rgp","cs","coolfact","gam","alpha_real0","n_real"
		  ,rg,rgp,cs,coolfact,gam,alpha_real0,n_real
		  ,"nu_vnr","nu_l","nu_ten","cour","cour2"
		  ,nu_vnr,nu_l,nu_ten,cour,cour2
		  ,"GRAVC","MASSBH","invsol2","blackholejz"
		  ,GRAVC,MASSBH,invsol2,blackholejz
		  ,"tstart","tf","tavgi","tavgf","numavg","timagescale"
		  ,tstart,tf,tavgi,tavgf,numavg,timagescale
		  ,"DTl","DTd","DTi","DTloss","DTfloor","DTtimestep","DTpd","DTener","DTfld","DTmode"
		  ,DTl,DTd,DTi,DTloss,DTfloor,DTtimestep,DTpd,DTener,DTfld,DTmode
		  ,"resist_real0","nu_sh"
		  ,resist_real0,nu_sh
		  ,"vg-x1","vg-x2","vg-x3"
		  ,vg[1],vg[2],vg[3]
		  ,"coord","fullvec","analoutput","DYNAMICMM"
		  ,COORD,REALNUMVEC==NUMVEC,analoutput,DYNAMICMM
		,"trans","transx1","transrhox1","transiex1","transv1x1","transv2x1","transv3x1","transmagx1","transx2","transrhox2","transiex2","transv1x2","transv2x2","transv3x2","transmagx2","transx3","transrhox3","transiex3","transv1x3","transv2x3","transv3x3"
		,trans,transx1,transrhox1,transiex1,transv1x1,transv2x1,transv3x1,transmagx1,transx2,transrhox2,transiex2,transv1x2,transv2x2,transv3x2,transmagx2,transx3,transrhox3,transiex3,transv1x3,transv2x3,transv3x3
		  ,"press","pressx1","pressx2","pressx3"
		  ,press,pressx1,pressx2,pressx3
		  ,"mag","transbzx","transbzy","stepmocct","mocctvx1","mocctvx2","mocctvx3","mocctbx1","mocctbx2","mocctbx3"
		  ,mag,transbzx,transbzy,stepmocct,mocctvx1,mocctvx2,mocctvx3,mocctbx1,mocctbx2,mocctbx3
		  ,"ie","visc_art","visc_real","vreal","vischeat","mdotin","cool","res_resl","rreal","resheat","advint","kever"
		  ,ie,visc_art,visc_real,vreal,vischeat,mdotin,cool,res_real,rreal,resheat,advint,kever
		  ,"intix1","intox1","intix2","intox2","intix3","intox3"
		  ,intix1,intox1,intix2,intox2,intix3,intox3
		  ,"nonunigridx1","nonunigridx2","nonunigridx3","simplebc","bcix1","bcox1","bcix2","bcox2","bcix3","bcox3"
		  ,nonunigridx1,nonunigridx2,nonunigridx3,simplebc,bcix1,bcox1,bcix2,bcox2,bcix3,bcox3
		  );
	}
	
	fclose(out);
      }
    }
  }

  fprintf(log_file,"end: init_outgparm\n"); fflush(log_file);
  return(0);
}



// 1/2/3-d valid
int init_pointers(void)
{

  fprintf(log_file,"begin: init_pointers ... "); fflush(log_file);

  dx=(FTYPE (*) [3][NBIGM]) (&(dxa[-1][-1][NBIGBND])); 
  x=(FTYPE (*) [3][NBIGM])(&(xa[-1][-1][NBIGBND]));

#if(POSTPROC==1)
  idx=(FTYPE (*) [3][INTNBIG]) (&(idxa[-1][-1][0])); 
  ix=(FTYPE (*) [3][INTNBIG])(&(ixa[-1][-1][0]));
  ppcalcs=(FTYPE (*) [N3M][N2M][N1M])(&(ppcalcsa[-1][N3BND][N2BND][N1BND]));
#else
  // fake assignment
  idx=(FTYPE (*) [3][INTNBIG])(0);
  ix=(FTYPE (*) [3][INTNBIG])(0);
  ppcalcs=(FTYPE (*) [N3M][N2M][N1M])(0);
#endif

  s=(FTYPE (*) [N3M][N2M][N1M])(&(sa[-1][N3BND][N2BND][N1BND]));
  v=(FTYPE (*) [3][N3M][N2M][N1M])(&(va[-1][-1][N3BND][N2BND][N1BND]));

  sanal=(FTYPE (*) [N3M][N2M][N1M])(&(sanala[-1][N3BND][N2BND][N1BND]));
  vanal=(FTYPE (*) [3][N3M][N2M][N1M])(&(vanala[-1][-1][N3BND][N2BND][N1BND]));

  tirr0=(FTYPE (*) [N2M][N1M])(&(tirr0a[N3BND][N2BND][N1BND]));
  opaini=(FTYPE (*) [N2M][N1M])(&(opainia[N3BND][N2BND][N1BND]));
  sigini=(FTYPE (*) [N2M][N1M])(&(siginia[N3BND][N2BND][N1BND]));
  optdini=(FTYPE (*) [N2M][N1M])(&(optdinia[N3BND][N2BND][N1BND]));
  cs02=(FTYPE (*) [N2M][N1M])(&(cs02a[N3BND][N2BND][N1BND]));
#if(DOINVSOL2FUN)
  invsol2fun=(FTYPE (*) [N2M][N1M])(&(invsol2funa[N3BND][N2BND][N1BND]));
#else
  invsol2fun=(FTYPE (*) [N2M][N1M])(0);
#endif

#if(GRAVACC&&(POSTPROC==0))
  gravacc=(FTYPE (*) [N3M][N2M][N1M])(&(gravacca[-1][N3BND][N2BND][N1BND]));
#else
  gravacc=(FTYPE (*) [N3M][N2M][N1M])(0);
#endif

#if(GRAVITOMAGNETIC&&(POSTPROC==0))
  gravitom=(FTYPE (*) [N3M][N2M][N1M])(&(gravitoma[-1][N3BND][N2BND][N1BND]));
#else
  gravitom=(FTYPE (*) [N3M][N2M][N1M])(0);
#endif

#if(DOAVGDIAGMEMORY)
  avg2_3=(FTYPE (*) [N2M][N1M])(&(avg2_3a[-1][N2BND][N1BND]));
  avg1_31=(FTYPE (*) [N2M])(&(avg1_31a[-1][N2BND]));
  avg1_32=(FTYPE (*) [N1M])(&(avg1_32a[-1][N1BND]));
#else
  avg2_3=(FTYPE (*) [N2M][N1M])(0);
  avg1_31=(FTYPE (*) [N2M])(0);
  avg1_32=(FTYPE (*) [N1M])(0);
#endif

#if(FLOORDUMPFLAGMEMORY)
  floorvars=(FTYPE (*) [N3M][N2M][N1M])(&(floorvarsa[-1][N3BND][N2BND][N1BND]));
  floorvarv=(FTYPE (*) [3][N3M][N2M][N1M])(&(floorvarva[-1][-1][N3BND][N2BND][N1BND]));
  floorvar0=(FTYPE (*) [N3M][N2M][N1M])(&(floorvar0a[-1][N3BND][N2BND][N1BND]));
#else
  floorvars=(FTYPE (*) [N3M][N2M][N1M])(0);
  floorvarv=(FTYPE (*) [3][N3M][N2M][N1M])(0);
  floorvar0=(FTYPE (*) [N3M][N2M][N1M])(0);
#endif

#if((LOWMEMMODE==0)&&(TVDLF==0))
  ds=(FTYPE (*) [3][N3M][N2M][N1M]) (&(dsa[-1][-1][N3BND][N2BND][N1BND])); 
  oarcl=(FTYPE (*) [3][N3M][N2M][N1M]) (&(oarcla[-1][-1][N3BND][N2BND][N1BND]));
  odx=(FTYPE (*) [3][NBIGM]) (&(odxa[-1][-1][NBIGBND]));
  ovol=(FTYPE (*) [N3M][N2M][N1M]) (&(ovola[-1][N3BND][N2BND][N1BND])); 
#else
  ds=(FTYPE (*) [3][N3M][N2M][N1M]) (0);
  oarcl=(FTYPE (*) [3][N3M][N2M][N1M]) (0);
  odx=(FTYPE (*) [3][NBIGM]) (0);
  ovol=(FTYPE (*) [N3M][N2M][N1M]) (0);
#endif

#if(LOWMEMMODE==0)
  g=(FTYPE (*) [NUMMETRIC][NBIGM])(&(ga[-1][-1][NBIGBND]));
  og=(FTYPE (*) [NUMMETRIC][NBIGM])(&(oga[-1][-1][NBIGBND]));
  dg=(FTYPE (*) [NUMMETRIC][NBIGM])(&(dga[-1][-1][NBIGBND]));
  dvl=(FTYPE (*) [3][NBIGM])(&(dvla[-1][-1][NBIGBND]));
  odvl=(FTYPE (*) [3][NBIGM])(&(odvla[-1][-1][NBIGBND]));
#else
  g=(FTYPE (*) [NUMMETRIC][NBIGM])(0);
  og=(FTYPE (*) [NUMMETRIC][NBIGM])(0);
  dg=(FTYPE (*) [NUMMETRIC][NBIGM])(0);
  dvl=(FTYPE (*) [3][NBIGM])(0);
  odvl=(FTYPE (*) [3][NBIGM])(0);
#endif

  
#if((PUREBC==0)&&(BOUNDTYPE==1))
  bcs=(short (*) [3][N3M][N2M][N1M])(&(bcsa[-1][-1][N3BND][N2BND][N1BND]));
  bcv=(short (*) [3][N3M][N2M][N1M])(&(bcva[-1][-1][N3BND][N2BND][N1BND]));

  // dummy
  bzmask=(short (*) [N2M][N1M])(0);
  bzmaskin=(short (*) [N2M][N1M])(0);

#elif((BOUNDTYPE==2)||(BOUNDTYPE==3))
  bzmask=(short (*) [N2M][N1M])(&(bzmaska[N3BND][N2BND][N1BND]));
  bzmaskin=(short (*) [N2M][N1M])(&(bzmaskina[N3BND][N2BND][N1BND]));

  // dummy
  bcs=(short (*) [3][N3M][N2M][N1M])(0);
  bcv=(short (*) [3][N3M][N2M][N1M])(0);

#else
  // dummy
  bzmask=(short (*) [N2M][N1M])(0);
  bzmaskin=(short (*) [N2M][N1M])(0);
  // dummy
  bcs=(short (*) [3][N3M][N2M][N1M])(0);
  bcv=(short (*) [3][N3M][N2M][N1M])(0);
#endif

#if((BOUNDTYPE==3)&&(POSTPROC==0))
  iindx=(short (*) [2][N3M][N2M])(&(iindxa[0][0][N3BND][N2BND]));
#else
  iindx=(short (*) [2][N3M][N2M])(0);
#endif

#if((LOOPTYPE==1)&&((DOLOSSDIAGMEMORY)&&(POSTPROC==0)) )
  losss=(SFTYPE (*) [3][2][NBIG][NBIG])(&(losssa[-1][-1][0][0]));
  lossv=(SFTYPE (*) [3+1][3][2][NBIG][NBIG])(&(lossva[-1][0][-1][0][0]));
  lossvisc=(SFTYPE (*) [3][2][NBIG][NBIG])(&(lossvisca[-1][-1][0][0]));
#else
  losss=(SFTYPE (*) [3][2][NBIG][NBIG])(0);
  lossv=(SFTYPE (*) [3+1][3][2][NBIG][NBIG])(0);
  lossvisc=(SFTYPE (*) [3][2][NBIG][NBIG])(0);
#endif

#if(MDOTMEM&&(POSTPROC==0))
  rhoinject=(FTYPE (*) [N2M][N1M])(&(rhoinjecta[N3BND][N2BND][N1BND]));
  eninject=(FTYPE (*) [N2M][N1M])(&(eninjecta[N3BND][N2BND][N1BND]));
#else
  rhoinject=(FTYPE (*) [N2M][N1M])(0);
  eninject=(FTYPE (*) [N2M][N1M])(0);
#endif

#if((TVDLF==0)&&(POSTPROC==0))
  worksbc=(FTYPE (*) [2][COMPDIM*2][NBIGBND*NBIGSM])(&(worksbca[-1][-1][0][0]));
  workvbc=(FTYPE (*) [2][COMPDIM*2][3*NBIGBND*NBIGSM])(&(workvbca[-1][-1][0][0]));
#else
  worksbc=(FTYPE (*) [2][COMPDIM*2][NBIGBND*NBIGSM])(0);
  workvbc=(FTYPE (*) [2][COMPDIM*2][3*NBIGBND*NBIGSM])(0);
#endif



  accountstore=(short (*) [N2M][N1M])(&(accountstorea[N3BND][N2BND][N1BND]));

#if(VISCMEM)
  nu_fact=(FTYPE (*) [N2M][N1M])(&(nu_facta[N3BND][N2BND][N1BND]));
  nu_real=(FTYPE (*) [N2M][N1M])(&(nu_reala[N3BND][N2BND][N1BND]));
  alpha_reala = (FTYPE(*)[N2M][N1M]) (&(alpha_realar[N3BND][N2BND][N1BND]));
#else
  nu_fact=(FTYPE (*) [N2M][N1M])(0);
  nu_real=(FTYPE (*) [N2M][N1M])(0);
#endif
  

#if(RESMEM&&(POSTPROC==0))
  nu_res_fact=(FTYPE (*) [N2M][N1M])(&(nu_res_facta[N3BND][N2BND][N1BND]));
  nu_res_real=(FTYPE (*) [N2M][N1M])(&(nu_res_reala[N3BND][N2BND][N1BND]));
  jcurrent=(FTYPE (*) [N3M][N2M][N1M])(&(jcurrenta[-1][N3BND][N2BND][N1BND]));
#else
  nu_res_fact=(FTYPE (*) [N2M][N1M])(0);
  nu_res_real=(FTYPE (*) [N2M][N1M])(0);
  jcurrent=(FTYPE (*) [N3M][N2M][N1M])(0);
#endif

  work1=(FTYPE (*) [N2M][N1M])(&(work1a[N3BND][N2BND][N1BND]));
  work2=(FTYPE (*) [N2M][N1M])(&(work2a[N3BND][N2BND][N1BND]));
  work3=(FTYPE (*) [N2M][N1M])(&(work3a[N3BND][N2BND][N1BND]));
#if(TVDLF==0)
  work4=(FTYPE (*) [N2M][N1M])(&(work4a[N3BND][N2BND][N1BND]));
  work5=(FTYPE (*) [N2M][N1M])(&(work5a[N3BND][N2BND][N1BND]));
  work6=(FTYPE (*) [N2M][N1M])(&(work6a[N3BND][N2BND][N1BND]));
  work7=(FTYPE (*) [N2M][N1M])(&(work7a[N3BND][N2BND][N1BND]));
  work8=(FTYPE (*) [N2M][N1M])(&(work8a[N3BND][N2BND][N1BND]));
  work9=(FTYPE (*) [N2M][N1M])(&(work9a[N3BND][N2BND][N1BND]));
  work10=(FTYPE (*) [N2M][N1M])(&(work10a[N3BND][N2BND][N1BND]));
#else
  work4=(FTYPE (*) [N2M][N1M])(0);
  work5=(FTYPE (*) [N2M][N1M])(0);
  work6=(FTYPE (*) [N2M][N1M])(0);
  work7=(FTYPE (*) [N2M][N1M])(0);
  work8=(FTYPE (*) [N2M][N1M])(0);
  work9=(FTYPE (*) [N2M][N1M])(0);
  work10=(FTYPE (*) [N2M][N1M])(0);
#endif

#if(MDOTMEMANAL&&(POSTPROC==0))
  mdotanal=(FTYPE (*) [N2M][N1M])(&(mdotanala[N3BND][N2BND][N1BND]));
#else
  mdotanal=(FTYPE (*) [N2M][N1M])(0);
#endif

  osqrtrho=(FTYPE (*) [N2M][N1M])(&(osqrtrhoa[N3BND][N2BND][N1BND]));

#if(ALFVENLIMIT)
  rholimited=(FTYPE (*) [N3M][N2M][N1M])(&(rholimiteda[0][N3BND][N2BND][N1BND]));
#else
  rholimited=(FTYPE (*) [N3M][N2M][N1M])(0);
#endif

#if(TVDLF==0)
  workv1=(FTYPE (*) [N3M][N2M][N1M])(&(workv1a[-1][N3BND][N2BND][N1BND]));
  workv2=(FTYPE (*) [N3M][N2M][N1M])(&(workv2a[-1][N3BND][N2BND][N1BND]));
  workv3=(FTYPE (*) [N3M][N2M][N1M])(&(workv3a[-1][N3BND][N2BND][N1BND]));
  workv4=(FTYPE (*) [N3M][N2M][N1M])(&(workv4a[-1][N3BND][N2BND][N1BND]));
  workv5=(FTYPE (*) [N3M][N2M][N1M])(&(workv5a[-1][N3BND][N2BND][N1BND]));
#else
  workv1=(FTYPE (*) [N3M][N2M][N1M])(0);
  workv2=(FTYPE (*) [N3M][N2M][N1M])(0);
  workv3=(FTYPE (*) [N3M][N2M][N1M])(0);
  workv4=(FTYPE (*) [N3M][N2M][N1M])(0);
  workv5=(FTYPE (*) [N3M][N2M][N1M])(0);
#endif

#if(VISCMEM)
  sigma=(FTYPE (*) [3][N3M][N2M][N1M])(&(sigmaa[-1][-1][N3BND][N2BND][N1BND]));
  rost=(FTYPE (*) [3][N3M][N2M][N1M])(&(rosta[-1][-1][N3BND][N2BND][N1BND]));
  rostnu=(FTYPE (*) [3][N3M][N2M][N1M])(&(rostnua[-1][-1][N3BND][N2BND][N1BND]));
#else
  sigma=(FTYPE (*) [3][N3M][N2M][N1M])(0);
  rost=(FTYPE (*) [3][N3M][N2M][N1M])(0);
  rostnu=(FTYPE (*) [3][N3M][N2M][N1M])(0);
#endif

  // only interp/image averages/ke floor if pp
#if(POSTPROC==1)
  workiq=(FTYPE (*) [INTN2][INTN1])(&(workiqa[-1][0][0]));
  workviq=(FTYPE (*) [3][INTN2][INTN1])(&(workviqa[-1][-1][0][0]));
  workiqavg=(FTYPE (*) [INTN2][INTN1])(&(workiqavga[-1][0][0]));
  work0iq=(FTYPE (*) [INTN2][INTN1])(&(work0iqa[-1][0][0]));
#else
  // fake assignment
  workiq=(FTYPE (*) [INTN2][INTN1])(0);
  workviq=(FTYPE (*) [3][INTN2][INTN1])(0);
  workiqavg=(FTYPE (*) [INTN2][INTN1])(0);
  work0iq=(FTYPE (*) [INTN2][INTN1])(0);
#endif

#if(POSTPROC==1)
  convmat=(FTYPE (*) [3][INTN2][INTN1])(&(convmata[0][0][0][0]));
#else
  // fake assignment
  convmat=(FTYPE (*) [3][INTN2][INTN1])(0);
#endif

  //ptraddr(-1);
#if(RAD)
// init radiation array
  liib= (int (*))(&(liiba[N2BND]));
  loib= (int (*))(&(loiba[N2BND]));
  lijb= (int (*))(&(lijba[N1BND]));
  lojb= (int (*))(&(lojba[N1BND]));
  niib= (int (*))(&(niiba[N2BND]));
  noib= (int (*))(&(noiba[N2BND]));
  nijb= (int (*))(&(nijba[N1BND]));
  nojb= (int (*))(&(nojba[N1BND]));
  er= (FTYPE (*)[N2M][N1M])(&(era[N3BND][N2BND][N1BND]));
  ern= (FTYPE (*)[N2M][N1M])(&(erna[N3BND][N2BND][N1BND]));
  e= (FTYPE (*)[N2M][N1M])(&(ea[N3BND][N2BND][N1BND]));
  en= (FTYPE (*)[N2M][N1M])(&(ena[N3BND][N2BND][N1BND]));
  der = (FTYPE (*)[N2M][N1M])(&(dera[N3BND][N2BND][N1BND]));
  de = (FTYPE (*)[N2M][N1M])(&(dea[N3BND][N2BND][N1BND]));
  fn1 = (FTYPE (*)[N2M][N1M])(&(fn1a[N3BND][N2BND][N1BND]));
  fn2 = (FTYPE (*)[N2M][N1M])(&(fn2a[N3BND][N2BND][N1BND]));
  dfn1de = (FTYPE (*)[N2M][N1M])(&(dfn1dea[N3BND][N2BND][N1BND]));
  dfn2de = (FTYPE (*)[N2M][N1M])(&(dfn2dea[N3BND][N2BND][N1BND]));
  dfn1der = (FTYPE (*)[N2M][N1M])(&(dfn1dera[N3BND][N2BND][N1BND]));
  dfn2der = (FTYPE (*)[N2M][N1M])(&(dfn2dera[N3BND][N2BND][N1BND]));
  dv = (FTYPE (*)[3][N3M][N2M][N1M])(&(dva[-1][-1][N3BND][N2BND][N1BND]));
  divv= (FTYPE (*)[N2M][N1M])(&(divva[N3BND][N2BND][N1BND]));
  f= (FTYPE (*)[2][N3M][N2M][N1M])(&(fa[-1][-1][N3BND][N2BND][N1BND]));
  dr = (FTYPE (*)[N3M][N2M][N1M])(&(dra[-1][N3BND][N2BND][N1BND]));
  fr = (FTYPE (*)[N3M][N2M][N1M])(&(fra[-1][N3BND][N2BND][N1BND]));
  kap = (FTYPE (*)[N2M][N1M])(&(kapa[N3BND][N2BND][N1BND]));
  kapn = (FTYPE (*)[N2M][N1M])(&(kapna[N3BND][N2BND][N1BND]));
  sig = (FTYPE (*)[N2M][N1M])(&(siga[N3BND][N2BND][N1BND]));
  dkapde = (FTYPE (*)[N2M][N1M])(&(dkapdea[N3BND][N2BND][N1BND]));
  bb = (FTYPE (*)[N2M][N1M])(&(bba[N3BND][N2BND][N1BND]));
  bbn = (FTYPE (*)[N2M][N1M])(&(bbna[N3BND][N2BND][N1BND]));
  dbbde = (FTYPE (*)[N2M][N1M])(&(dbbdea[N3BND][N2BND][N1BND]));
  pre = (FTYPE (*)[N2M][N1M])(&(prea[N3BND][N2BND][N1BND]));
  pren =(FTYPE (*)[N2M][N1M])(&(prena[N3BND][N2BND][N1BND]));
  dpde = (FTYPE (*)[N2M][N1M])(&(dpdea[N3BND][N2BND][N1BND]));
  work1y= (FTYPE (*))(&(work1ya[N2BND]));
  work2y= (FTYPE (*))(&(work2ya[N2BND]));
  work3y= (FTYPE (*))(&(work3ya[N2BND]));
#endif
  work1x= (FTYPE (*))(&(work1xa[N1BND]));
  work2x= (FTYPE (*))(&(work2xa[N1BND]));
  work3x= (FTYPE (*))(&(work3xa[N1BND]));
  work4x= (FTYPE (*))(&(work4xa[N1BND]));
  work5x= (FTYPE (*))(&(work5xa[N1BND]));
  work6x= (FTYPE (*))(&(work6xa[N1BND]));
#if(KAICOOL)
   thingrid = (int (*)[N2M][N1M])(&(thingrida[N3BND][N2BND][N1BND]));
   sthin = (FTYPE (*)[N2M][N1M])(&(sthina[N3BND][N2BND][N1BND]));
   tot = (FTYPE (*)[N3M][N2M][N1M])(&(tota[0][N3BND][N2BND][N1BND]));
   optd = (FTYPE (*)[N2M][N1M])(&(optda[N3BND][N2BND][N1BND]));
   tatm4 = (FTYPE (*)[N2M][N1M])(&(tatm4a[N3BND][N2BND][N1BND]));
   kapk = (FTYPE (*)[N2M][N1M])(&(kapka[N3BND][N2BND][N1BND]));
   sigk = (FTYPE (*)[N2M][N1M])(&(sigka[N3BND][N2BND][N1BND]));
   tkai = (FTYPE (*)[N2M][N1M])(&(tkaia[N3BND][N2BND][N1BND]));
   fbdry = (FTYPE (*)[N2M][N1M])(&(fbdrya[N3BND][N2BND][N1BND]));
   taoup = (FTYPE (*)[N2M][N1M])(&(taoupa[N3BND][N2BND][N1BND]));
   tenv = (FTYPE (*))(&(tenva[N1BND]));
   coolerKai =(FTYPE (*)[N2M][N1M])(&(coolerKaia[N3BND][N2BND][N1BND]));
#endif
   diagn = (FTYPE (*)[N3M][N2M][N1M])(&(diagna[-1][N3BND][N2BND][N1BND]));

#if((TVDLF==1)&&(POSTPROC==0))
  p =   (FTYPE (*) [N2+2][NP])(& (  a_p[1][1][0])) ;
  ph =  (FTYPE (*) [N2+2][NP])(& ( a_ph[1][1][0])) ;
  Fx =  (FTYPE (*) [N2+2][NP])(& ( a_Fx[1][1][0])) ;
  Fy =  (FTYPE (*) [N2+2][NP])(& ( a_Fy[1][1][0])) ;
  dqx = (FTYPE (*) [N2+2][NP])(& (a_dqx[1][1][0])) ;
  dqy = (FTYPE (*) [N2+2][NP])(& (a_dqy[1][1][0])) ;
#else
  p =   (FTYPE (*) [N2+2][NP])(0);
  ph =  (FTYPE (*) [N2+2][NP])(0);
  Fx =  (FTYPE (*) [N2+2][NP])(0);
  Fy =  (FTYPE (*) [N2+2][NP])(0);
  dqx = (FTYPE (*) [N2+2][NP])(0);
  dqy = (FTYPE (*) [N2+2][NP])(0);
#endif
//surface density,sigma
  Sigar = (FTYPE (*)[N2M][N1M])(&(Sigara[N3BND][N2BND][N1BND]));
  cs2a = (FTYPE (*)[N2M][N1M])(&(cs2aa[N3BND][N2BND][N1BND]));

  // END: Setup arrays giving NBND boundary zones--makes 0,0 first active zone 
  fprintf(log_file,"end: init_pointers\n"); fflush(log_file);
  return(0);
}


// Notes on dx(dvl):
// -----
// Without any interpolation, must setup dx to be consistent at boundary zones with spatial structure.
// So I use BC interpolation to avoid worrying about this :)  Kinda ingore below
// ------
// Periodic-> spatial condition.  dx must be periodic-> dx(ghost)=dx(what copied), dx as scalar
// Reflect-> spatial condition. dx(-1)=dx(0), etc., dx as scalar
// AOS-> as with Reflect, dx as scalar, do not invert dx->-dx for z-axis
//
// Rest do not give spatial condition, all spatial info indeterminate.
// Outflow: Soft Edge: ~Neuman, Like bigger domain existed or true outflow to edge.
// Inflow: Hard Edge: Dirichlet: or not really hard edge, like bigger domain but with influx of matter.

// Outflow using copy needs Reflect-like spatial condition(i.e. local copy of dx to boundary)
// Outflow using interpolation using position does not care what dx is on bc.  Probably reasonable to make sure dx does not put middle of bc zone at x=0 for non-cart coords.

// Inflow is just Dirichlet, does not care what dx is, as with outflow interpolated.  Make sure to use correct position in condition(center/edge of zone, etc.).

// Below I have chosen out/in-flow to be treated same as with reflect for dx/dvl unless I decide otherwise
// Might be interesting to set dx large for outflow, as if domain existed, but only as average value.  Probably need v to be linearly interpolated within zone to be really meaninglful.

// Necessary dxs:
// besides active grid, must at least have dx on:
// for advection terms(for fluxes(sweepy only)) and dq for sweepx:
// advection(sweepy):
// dx1a(i) <> i=0..N1-1, j=-1..N2
// dx1b(i) <> i=0..N1-1, j=0..N2
// dq(sweepx):
// dx1b(i) <> i=0..N1, j=0..N2-1
// dx1a(i) <> i=-1..N1-1, j=0..N2-1
// dq(sweepy):
// dx2b(i) <> i=0..N1-1, j=0..N2
// dx2a(i) <> i=0..N1-1, j=-1..N2-1


// 1/2/3-d valid
// which: 0 = real grid 1=interp grid
// outtype: for interp grid only: what kind of grid:
// 0=native type
// 1=linear cart
int init_dx(int n1,int n2,int n3,int n1b,int n2b,int n3b, SFTYPE*sx,int which, int outtype)
{
  SFTYPE frac,Fracin,Rin,Width,FullLength;
  int i,j,k,l,m;
  int N[3+1],NB[3+1];
  int NC;
  SFTYPE dellogx;
  SFTYPE cc,cw,dx1,dx2,dx3,Ne1,Ne2,Ne3,dxe,dx0;
  SFTYPE bcoef,bcoefm,bcoefp,acoef,ccoef,dcoef,truea;
  SFTYPE acoefi;
  SFTYPE acoeff;
  SFTYPE ftemp;
  SFTYPE xasofar;
  SFTYPE ftemp1,ftemp2;
  SFTYPE Ref,deltheta,delj,j0,EE,Sum,dtf,dt0;
  SFTYPE N2U;
  int itemp;
  int realsize[3+1];
  
  // for nonunigrid==7
  int dir;
  FTYPE AA,BB,N0,NN,R0,Rminus,Rplus,dRupper0,dRlower0;


  fprintf(log_file,"begin: init_dx ... "); fflush(log_file);

  fprintf(log_file,"nonunigridx1: %d totalsize[1]: %d L[2][1]: %21.15g\n",nonunigridx2,totalsize[1],L[2][1]);  fflush(log_file);


  if((which==1)&&(POSTPROC==0)){
    fprintf(fail_file,"Not setup to do postproc during a run: init_dx: which=%d\n",which);
    myexit(1);
  }

  N[1]=n1;
  N[2]=n2;
  N[3]=n3;
  NB[1]=n1b;
  NB[2]=n2b;
  NB[3]=n3b;

  /* BEGIN: Apply global parameters to cell structure */

  /* Assign dxs--uniform grid */
  /* Note this data would be dynamic if dynamic grid */

  // a-grid


  // X1
  // add up sizes for total size of grid

  if((which==0)||(outtype==0)){

    if(which==1){
      iL[1][1]=L[1][1];
      iL[2][1]=L[2][1];
    }

    if(nonunigridx1==0){

      xasofar=L[1][1];
      
      for(i=-NB[1];i<totalsize[1]+NB[1];i++){
	
        ftemp=L[2][1]/(SFTYPE)totalsize[1];
 
#if(USEMPI==0)
        if(i==0) sx[1]=xasofar;
        if(which==0) dx[1][1][i]=ftemp;
        else if(which==1) idx[1][1][i]=ftemp;
#else
        if( (i>=startpos[1]-NB[1])&&(i<=endpos[1]+NB[1]) ){
 
          if(i==startpos[1]) sx[1]=xasofar;
          if(which==0) dx[1][1][i-startpos[1]]=ftemp;
          else if(which==1) idx[1][1][i-startpos[1]]=ftemp;
 
        }
#endif
        if( (i>=0)&&(i<totalsize[1]) ){
          xasofar+=ftemp;
        }
      }
      // MARK
      //dx[1][1][-1]=dx[1][1][0]/1000.0;
      //dx[1][1][-2]=dx[1][1][0]/1000.0;
    }
    else if(nonunigridx1==1){
      // use mathematica page to get these values or pick from domain determined below:
      NC=(int)(totalsize[1]/3); // specify how many zones for cosine region
      if(NC!=totalsize[1]){
	cw=1.0;
	cc=1.0;
	dx1=(cc-cw*.5)-L[1][1];
	dx2=L[2][1]+L[1][1]-(cs+cw*.5);
	dxe=(dx1+dx2)/((SFTYPE)(totalsize[1]-NC));
	Ne1=rint(dx1/dxe);
	Ne2=rint(dx2/dxe);
	if(cw>L[2][1]){
	  fprintf(fail_file,"Can't have coswidth greater than real width\n");
	  myexit(1);
	}
	if( (cc-cw/2<L[1][1])||(cc+cw/2>L[1][1]+L[2][1])){
	  fprintf(fail_file,"Can't have edges of cos outside domain!\n");
	  myexit(1);
	}
      }

      xasofar=L[1][1];
      
      for(i=-NB[1];i<totalsize[1]+NB[1];i++){
	
	if(NC!=totalsize[1]){
	  if((i>=-NB[1])&&(i<(int)(Ne1))){
	    ftemp=dxe;
	  }
	  else if( (i<=totalsize[1]-1+NB[1])&&(i>=totalsize[1]-(int)(Ne2)) ){
	    ftemp=dxe;
	  }
	  else if( (i>=Ne1)&&(i<totalsize[1]-Ne2) ){
	    // no acoef here(use l21-dx2-dx1 instead of cw for fix on aliasing accuracy
	    bcoef=((SFTYPE)(NC)*dxe-(L[2][1]-dx1-dx2))/((SFTYPE)(NC-1.));
	    //fprintf(stderr,"%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",dx1,dx2,dxe,Ne1,Ne2,bcoef);
	    if(bcoef<0){
	      fprintf(fail_file,"bcoef<0 does not produce intended results!\n");
	      myexit(1);
	    }
	    if(bcoef>dxe*0.5){
	      fprintf(fail_file,"bcoef>dxe/2 does not produce intended results!\n");
	      myexit(1);
	    }
	    ftemp=(dxe-bcoef)+bcoef*cos(2.*M_PI*(SFTYPE)(i-Ne1)/((SFTYPE)(NC-1)));
	  }
	}
	else{
	  acoeff= (L[2][1]/((SFTYPE)totalsize[1]));
	  acoefi= (L[2][1]/((SFTYPE)(totalsize[1]+1)));
	  // must have a>b for nonnegative dx1
	  acoef=(acoeff-acoefi)*.1+acoefi; // add some% off inner acoef to get smallest dx in center
	  bcoef=L[2][1]-(SFTYPE)(totalsize[1])*acoef;
	  
	  ftemp=acoef+bcoef*cos(2.*M_PI*i/((SFTYPE)(totalsize[1]-1)));
	}
	
	
#if(USEMPI==0)
	if(i==0) sx[1]=xasofar;
	if(which==0) dx[1][1][i]=ftemp;
	else if(which==1) idx[1][1][i]=ftemp;
#else
	if( (i>=startpos[1]-NB[1])&&(i<=endpos[1]+NB[1]) ){
	  
	  if(i==startpos[1]) sx[1]=xasofar;
	  if(which==0) dx[1][1][i-startpos[1]]=ftemp;
	  else if(which==1) idx[1][1][i-startpos[1]]=ftemp;
	  
	}
#endif
	if( (i>=0)&&(i<totalsize[1]) ){
	  xasofar+=ftemp;
	}
      }
    }
    else if(nonunigridx1==2){
      ccoef=0.0; // leave at 0.0
      // bcoef starts at >3.000000 and goes to infinity(100.0~uniform) 
      bcoef=3.1; // free
      acoef=(L[2][1]-ccoef*(SFTYPE)totalsize[1])/(alnfact(bcoef+(SFTYPE)totalsize[1]-2)-alnfact(bcoef-2));
      
      xasofar=L[1][1];
 
      for(i=-NB[1];i<totalsize[1]+NB[1];i++){
 
        ftemp=ccoef+acoef*log(i+bcoef);
 
#if(USEMPI==0)
        if(i==0) sx[1]=xasofar;
        if(which==0) dx[1][1][i]=ftemp;
        else if(which==1) idx[1][1][i]=ftemp;
#else
        if( (i>=startpos[1]-NB[1])&&(i<=endpos[1]+NB[1]) ){
 
          if(i==startpos[1]) sx[1]=xasofar;
          if(which==0) dx[1][1][i-startpos[1]]=ftemp;
          else if(which==1) idx[1][1][i-startpos[1]]=ftemp;
 
        }
#endif
        if( (i>=0)&&(i<totalsize[1]) ){
          xasofar+=ftemp;
        }
      }                  

      /*
      for(i=-NB[1];i<N[1]+NB[1];i++){    
	ftemp=ccoef+acoef*log(i+bcoef);
	if(which==0) dx[1][1][i]=ftemp;
	else if(which==1) idx[1][1][i]=ftemp;
      }
      */
    }
    else if((nonunigridx1==3)||(nonunigridx1==4)){
      acoef=10.0; // Factor by which each larger decade dx is to last decade.
      ftemp1=L[1][1]+L[2][1];
      ftemp2=L[1][1];
      bcoef=(SFTYPE)(totalsize[1])/(log10(ftemp1)-log10(ftemp2)); // Number of zones per decade      
      ccoef=(pow(acoef,(SFTYPE)(totalsize[1])/(bcoef-1.))-1.)/(pow(acoef,1./(bcoef-1.))-1.);
      dx0=L[2][1]/ccoef;

      if(nonunigridx1==3) xasofar=L[1][1];
      else xasofar=dx0*(bcoef-1.0)/log(acoef); // forces inner value
      for(i=-NB[1];i<totalsize[1]+NB[1];i++){
 
        ftemp=pow(acoef,(SFTYPE)(i)/(bcoef-1.))*dx0;
 
#if(USEMPI==0)
        if(i==0) sx[1]=xasofar;
        if(which==0) dx[1][1][i]=ftemp;
        else if(which==1) idx[1][1][i]=ftemp;
#else
        if( (i>=startpos[1]-NB[1])&&(i<=endpos[1]+NB[1]) ){
 
          if(i==startpos[1]) sx[1]=xasofar;
          if(which==0) dx[1][1][i-startpos[1]]=ftemp;
          else if(which==1) idx[1][1][i-startpos[1]]=ftemp;
 
        }
#endif
        if( (i>=0)&&(i<totalsize[1]) ){
          xasofar+=ftemp;
        }
      }

      if(myid<=0){
	fprintf(logfull_file,"Factor: %15.10g Nr: %15.10g ccoef: %15.10g dx0: %15.10g\n",acoef,bcoef,ccoef,dx0);
      }
      /*      
      for(i=-NB[1];i<N[1]+NB[1];i++){ // compute as factor of dx0 so don't accumulate error
	ftemp=pow(acoef,(SFTYPE)(i)/(bcoef-1.))*dx0;
	if(which==0) dx[1][1][i]=ftemp;
	else if(which==1) idx[1][1][i]=ftemp;
      }
      */
    }
    else if(nonunigridx1==5){
      // global stuff
      //Rin=3.0*rgp; // location of inner log edge
      Rin=L[1][1]; // normal old method
      if(Rin>L[1][1]){
	Fracin=ceil(totalsize[1]/20.0); // number of zones to resolve inside Rin if Rin!=L[2][1]
      }
      else Fracin=0;
      FullLength=L[1][1]+L[2][1];
      Width=L[1][1]+L[2][1]-Rin;

      acoef=10.0; // Factor by which each larger decade dx is to last decade.
      bcoef=10.0; // definition of decade in this base in log(useless)
      ftemp1=log(Rin-rgp)/log(acoef); // xin
      ftemp2=log(FullLength-rgp)/log(acoef); // xout
      ccoef=-log( (Rin-rgp)/(FullLength-rgp))/log(bcoef); // Number of decades
      dellogx=(ftemp2-ftemp1)/((SFTYPE)(totalsize[1]-Fracin)); // dx in log
      if(myid<=0){
	fprintf(logfull_file,"xin: %15.10g xout: %15.10g F: %15.10g D: %15.10g Nd: %15.10g dellogx: %15.10g\n",ftemp1,ftemp2,acoef,bcoef,ccoef,dellogx);
	fflush(logfull_file);
      }

      xasofar=L[1][1];
      for(i=-NB[1];i<totalsize[1]+NB[1];i++){

	if((i<=Fracin)&&(Rin>L[1][1])){
	  if(i>=0){
	    ftemp=(Rin-L[1][1])/((SFTYPE)Fracin); // uniform grid up to Rin
	  }
	  else ftemp=(Rin-L[1][1])/((SFTYPE)Fracin)/10.0; // want ghost zones to not push too deep to rgp
	}
	else ftemp=pow(acoef,ftemp1+(i-Fracin)*dellogx)*(pow(acoef,dellogx)-1.0);

#if(USEMPI==0)
        if(i==0) sx[1]=xasofar;
        if(which==0) dx[1][1][i]=ftemp;
        else if(which==1) idx[1][1][i]=ftemp;
#else
        if( (i>=startpos[1]-NB[1])&&(i<=endpos[1]+NB[1]) ){
	  
          if(i==startpos[1]) sx[1]=xasofar;
          if(which==0) dx[1][1][i-startpos[1]]=ftemp;
          else if(which==1) idx[1][1][i-startpos[1]]=ftemp;
	  
        }
#endif
        if( (i>=0)&&(i<totalsize[1]) ){
          xasofar+=ftemp;
        }
      }
    }
    else if(nonunigridx1==7){

      dir=1;
      xasofar=L[1][dir]; // inner total grid position

      // see mathematica file averywind_gridsetup.nb to choose parameters
      // Constants below are true for N=128, x1in=-1/4, x1out=3/4, posaccretor=0
      NN = (FTYPE)totalsize[dir];
      R0 = pos_accretor[dir];
      if(totalsize[dir]==128){
	AA = 41.136481433951424;
	BB = 122.08219827806326;
	N0 = 56.0;
      }
      else if(totalsize[dir]==64){
	AA = 94.7668863043728;
	BB = 268.9467546261696;
	N0 = 28.0;
      }
      else if(totalsize[dir]==32){
#if(0)
	// C = 2500
	AA = 197.54656881541703;
	BB = 538.5208544314766;
	N0 = 14.0;
#else
	// C = 320
	AA = 18.17473826949555;
	BB = 38.84343841422464;
	N0 = 12.0;
#endif
      }
	
	
      Rminus = L[1][dir]; // inner edge position (like x1in)
      Rplus = L[1][dir]+L[2][dir]; // outer edge position (like x1out)
	
      dRupper0=(Rplus-R0)*(pow(BB,1.0/(NN-N0))-1.0)/(BB-1.0);
      dRlower0=(R0-Rminus)*(pow(AA,1.0/N0)-1.0)/(AA-1.0);
      
      for(i=-NB[dir];i<totalsize[dir]+NB[dir];i++){
	
	if(i>=N0){
	  // then use dRupper
	  ftemp = dRupper0*pow(BB,(i-N0)/(NN-N0));
	}
	else{
	  // then use dRlower
	  ftemp = dRlower0*pow(AA,(N0-i-1)/N0);
	}
	  

 
#if(USEMPI==0)
	if(i==0) sx[1]=xasofar;
	if(which==0) dx[1][1][i]=ftemp;
	else if(which==1) idx[1][1][i]=ftemp;
#else
	if( (i>=startpos[1]-NB[1])&&(i<=endpos[1]+NB[1]) ){
 
	  if(i==startpos[1]) sx[1]=xasofar;
	  if(which==0) dx[1][1][i-startpos[1]]=ftemp;
	  else if(which==1) idx[1][1][i-startpos[1]]=ftemp;
 
	}
#endif
	if( (i>=0)&&(i<totalsize[1]) ){
	  xasofar+=ftemp;
	}
      }
    }
    else{
      fprintf(fail_file,"No known X1 grid type specified\n");
      myexit(1);
    }

    // override to make reflective copy work right in nonunigrid if reflective
    if(reflectix1){ // fix the boundary zones to be "reflective", so boundary code is less complicated for field(divB=0 vs just copy)  (really needed for couple corner zones)
      // this will only catch for correct cpus
      if(which==0){
	dx[1][1][-1]=dx[1][1][0];
	dx[1][1][-2]=dx[1][1][1];
      }
      if(which==1){
	idx[1][1][-1]=idx[1][1][0];
	idx[1][1][-2]=idx[1][1][1];
      }
    }
    if(reflectox1){
      // this will only catch for correct cpus
      if(which==0){
	dx[1][1][N1]=dx[1][1][N1-1];
	dx[1][1][N1+1]=dx[1][1][N1-2];
      }
      if(which==1){
	idx[1][1][N1]=idx[1][1][N1-1];
	idx[1][1][N1+1]=idx[1][1][N1-2];
      }
    }


      
    
    // X2

    fprintf(log_file,"nonunigridx2: %d totalsize[2]: %d L[2][2]: %21.15g\n",nonunigridx2,totalsize[2],L[2][2]);  fflush(log_file);

    if(which==1){
      iL[1][2]=L[1][2];
      iL[2][2]=L[2][2];
    }
          
    if(nonunigridx2==0){  
      xasofar=L[1][2];
      for(j=-NB[2];j<totalsize[2]+NB[2];j++){
	
	ftemp=L[2][2]/((SFTYPE)(totalsize[2]));
	
#if(USEMPI==0)
	if(j==0) sx[2]=xasofar;
	if(which==0) dx[1][2][j]=ftemp;
	else if(which==1) idx[1][2][j]=ftemp;
#else
	if( (j>=startpos[2]-NB[2])&&(j<=endpos[2]+NB[2]) ){
	  
	  if(j==startpos[2]) sx[2]=xasofar;
	  if(which==0) dx[1][2][j-startpos[2]]=ftemp;
	  else if(which==1) idx[1][2][j-startpos[2]]=ftemp;
	  
	}
#endif
	if( (j>=0)&&(j<totalsize[2]) ){
	  xasofar+=ftemp;
	}
      }
    }
    else if(nonunigridx2==1){
      acoeff=(L[2][2]/((SFTYPE)totalsize[2]));
      acoefi=(L[2][2]/((SFTYPE)(totalsize[2]+1)));
      // must have a>b for nonnegative dx2
      acoef=(acoeff-acoefi)*.9+acoefi; // add some% off inner acoef to get smallest dx in center
      bcoef=L[2][2]-(SFTYPE)(totalsize[2])*acoef;
      
      xasofar=L[1][2];
      for(j=-NB[2];j<totalsize[2]+NB[2];j++){    
	ftemp=acoef+bcoef*cos(2.*M_PI*j/((SFTYPE)(totalsize[2]-1)));
	
#if(USEMPI==0)
	if(j==0) sx[2]=xasofar;
	if(which==0) dx[1][2][j]=ftemp;
	else if(which==1) idx[1][2][j]=ftemp;
#else 
	if( (j>=startpos[2]-NB[2])&&(j<=endpos[2]+NB[2]) ){
	  
	  if(j==startpos[2]) sx[2]=xasofar;
	  if(which==0) dx[1][2][j-startpos[2]]=ftemp;
	  else if(which==1) idx[1][2][j-startpos[2]]=ftemp;
	}
#endif
	if( (j>=0)&&(j<totalsize[2]) ){
	  xasofar+=ftemp;
	} 
      }
    }
    else if(nonunigridx2==2){
      ccoef=0.0; // leave at 0.0
      // bcoef starts at >2.000000 and goes to infinity(100.0~uniform) 
      bcoef=2.1; // free
      acoef=(0.52*L[2][2]-0.5*ccoef*totalsize[2])/(alnfact(bcoef+0.5*totalsize[2]-1.)-alnfact(bcoef-1.));
      
      //printf("acoef: %21.10g\n",acoef);
      xasofar=L[1][2];
      for(j=-NB[2];j<totalsize[2]+NB[2];j++){
	if(j>=totalsize[2]/2){
	  ftemp=ccoef+acoef*log((SFTYPE)j-(SFTYPE)totalsize[2]*(SFTYPE)numprocs*0.5-1.0+bcoef);
	}
	else if(j<=totalsize[2]/2-1){
	  ftemp=ccoef+acoef*log((SFTYPE)totalsize[2]*(SFTYPE)numprocs*0.5-1.0-(SFTYPE)j-1.0+bcoef);
	}
#if(USEMPI==0)
	if(j==0) sx[2]=xasofar;
	if(which==0) dx[1][2][j]=ftemp;
	else if(which==1) idx[1][2][j]=ftemp;
#else 
	if( (j>=startpos[2]-NB[2])&&(j<=endpos[2]+NB[2]) ){
	  
	  if(j==startpos[2]) sx[2]=xasofar;
	  if(which==0) dx[1][2][j-startpos[2]]=ftemp;
	  else if(which==1) idx[1][2][j-startpos[2]]=ftemp;
	}
#endif
	if( (j>=0)&&(j<totalsize[2]) ){
	  xasofar+=ftemp;
	}  
      }
    }
    else if(nonunigridx2==3){
      // N[2] needs to be factor of 2 right now
      
      // F // GODMARK(usually 10)
      acoef=4; // Factor by which each larger decade dx is to last decade. 
      dcoef=(SFTYPE)(totalsize[2])*0.5; // number of zones per decade for half the region (usually N/2 or something less)
      if(dcoef<=1.0){
	fprintf(fail_file,"Can't have this gridtype(nonunigridx2==%d) and totalN2=%d\n",nonunigridx2,totalsize[2]);
	myexit(1);
      }
      truea=pow(acoef,1.0/(dcoef-1.)); // true factor of increase per entire range
      // N
      bcoef=(SFTYPE)(totalsize[2]); // number of zones in Pi
      
      // dx0
      //ccoef=L[2][2]/( ( (pow(acoef,bcoef/(bcoef-2.0))-1.0)/(pow(acoef,2.0/(bcoef-2.0))-1.0) )+ ( (acoef-1.0)/(pow(acoef,2.0/bcoef)-1.0) ) );
      ccoef=L[2][2]*0.5*(truea-1.)/(pow(truea,bcoef*0.5)-1.0);
      xasofar=L[1][2];
      for(j=-NB[2];j<totalsize[2]+NB[2];j++){
	
	// real dxa(j) grid value as if whole grid here      
	if(j<=totalsize[2]/2-1){ // if in -j half
	  //ftemp=ccoef*pow(acoef,(bcoef*0.5-1.0-(SFTYPE)j)/(bcoef*0.5-1.0));
	  ftemp=ccoef*pow(truea,bcoef*0.5-1.0-(SFTYPE)j);
	}
	else if(j>=totalsize[2]/2){
	  //ftemp=ccoef*pow(acoef,(j-bcoef*0.5)/(bcoef*0.5));
	  ftemp=ccoef*pow(truea,(SFTYPE)(j)-bcoef*0.5);
	}
	
	// can use this below part in general for any grid!  Just make sure above calcs are based upon full grid, not each cpu
#if(USEMPI==0)
	if(j==0) sx[2]=xasofar;
	if(which==0) dx[1][2][j]=ftemp;
	else if(which==1) idx[1][2][j]=ftemp;
#else 
	if( (j>=startpos[2]-NB[2])&&(j<=endpos[2]+NB[2]) ){
	  
	  if(j==startpos[2]) sx[2]=xasofar;
	  if(which==0) dx[1][2][j-startpos[2]]=ftemp;
	  else if(which==1) idx[1][2][j-startpos[2]]=ftemp;
	}
#endif
	if( (j>=0)&&(j<totalsize[2]) ){
	  xasofar+=ftemp;
	}    
      }
    }
    else if(nonunigridx2==4){
      // Inverse Gaussian refinement (centers refinement around equator in Gaussian shape)

      // N[2] needs to be factor of 2 right now

      // user parameters
      Ref=10.0; // refinement factor (factor by which central region is refined over uniform grid)
      deltheta=.1*M_PI; // 1/e of Gaussian in theta
      
      // auto-defined parameters
      j0=((SFTYPE)totalsize[2]-1.0)/2.0; // centroid of Gaussian in j      
      delj=2.0*j0*deltheta/(M_PI*0.5); // approximate width of Gaussian in j, just can't solve for delj directly from deltheta without blowing a brain vessel)
      EE=exp(-pow(j0/delj,2.0));
      Sum=0.0; // Useful Sum
      for(j=0;j<=totalsize[2]-1;j++){
	Sum+=exp(-pow((j-j0)/delj,2.0));
      }
      dtf=(M_PI*EE*(1.0/Ref-1.0)+M_PI*(1.0-Sum/((SFTYPE)totalsize[2]*Ref)))/((SFTYPE)totalsize[2]-Sum);
      dt0=M_PI/((SFTYPE)totalsize[2]*Ref);
      bcoef=(dtf-dt0)/(EE-1.0);
      acoef=dt0-bcoef;

      fprintf(log_file,"%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",Ref,deltheta,delj,j0,EE,Sum,dtf,dt0,bcoef,acoef); fflush(log_file);
      // here begins general algorithm for assignment, ftemp below is dy(j)
      xasofar=L[1][2];
      for(j=-NB[2];j<totalsize[2]+NB[2];j++){

	ftemp=acoef+bcoef*exp(-pow((j-j0)/delj,2.0)); // specific value of dy(j)
	
	// can use this below part in general for any grid!  Just make sure above calcs are based upon full grid, not each cpu
#if(USEMPI==0)
	if(j==0) sx[2]=xasofar;
	if(which==0) dx[1][2][j]=ftemp;
	else if(which==1) idx[1][2][j]=ftemp;
#else 
	if( (j>=startpos[2]-NB[2])&&(j<=endpos[2]+NB[2]) ){
	  
	  if(j==startpos[2]) sx[2]=xasofar;
	  if(which==0) dx[1][2][j-startpos[2]]=ftemp;
	  else if(which==1) idx[1][2][j-startpos[2]]=ftemp;
	}
#endif
	if( (j>=0)&&(j<totalsize[2]) ){
	  xasofar+=ftemp;
	}    
      }
    }
    else if(nonunigridx2==5){
      // Gammie sinusoidal nonunigrid

      // N[2] needs to be factor of 2 right now

      // user parameters
      Ref=0.2; // refinement factor

      // here begins general algorithm for assignment, ftemp below is dy(j)
      xasofar=L[1][2];
      for(j=-NB[2];j<totalsize[2]+NB[2];j++){

	ftemp=(1.0+(1.0-Ref)*cos(2.0*M_PI/totalsize[2]*j))*M_PI/totalsize[2]; // specific value of dy(j)
	
	// can use this below part in general for any grid!  Just make sure above calcs are based upon full grid, not each cpu
#if(USEMPI==0)
	if(j==0) sx[2]=xasofar;
	if(which==0) dx[1][2][j]=ftemp;
	else if(which==1) idx[1][2][j]=ftemp;
#else 
	if( (j>=startpos[2]-NB[2])&&(j<=endpos[2]+NB[2]) ){
	  
	  if(j==startpos[2]) sx[2]=xasofar;
	  if(which==0) dx[1][2][j-startpos[2]]=ftemp;
	  else if(which==1) idx[1][2][j-startpos[2]]=ftemp;
	}
#endif
	if( (j>=0)&&(j<totalsize[2]) ){
	  xasofar+=ftemp;
	}    
      }
    }
    else if(nonunigridx2==6){ // split log with uniform central grid
      // N[2] needs to be factor of 2 right now
      
      // F // GODMARK(usually 10)
// Factor by which each larger decade dx is to last decade.
//      acoef=707.359;
      acoef=707.359; // HOR=.1
      //acoef=707.359; // HOR=.05
      //acoef=707.359; // HOR=.025
			//acoef=707.359; // HOR=.0125
      //acoef=1551.66; // HOR=.00625
			//acoef=3355.35; // HOR=.003125
// find above using mathematica (loghr.nb)
			N2U=N2/4; // specific fraction of zones taken by linear uniform central region
			// the above fraction must be exact! (i.e. 64/4 is ok)
      dcoef=(SFTYPE)(totalsize[2]*0.5-N2U*1.0); // number of zones per decade for half the region (usually N/2 or something less)
      if(dcoef<=1.0){
	fprintf(fail_file,"Can't have this gridtype(nonunigridx2==%d) and totalN2=%d\n",nonunigridx2,totalsize[2]);
	myexit(1);
      }
      truea=pow(acoef,1.0/(dcoef-1.)); // true factor of increase per entire range
      // N
      bcoef=(SFTYPE)(totalsize[2]*0.5-N2U*1.0); // # zones
      bcoefm=(SFTYPE)(totalsize[2]*0.5-N2U*1.0); // position of maximimal res (minus)
      bcoefp=(SFTYPE)(totalsize[2]*0.5+N2U*1.0); // position of maximimal res (plus)
      
      ccoef=(L[2][2]*0.5-HOR)*(truea-1.)/(pow(truea,bcoef)-1.0);

      xasofar=L[1][2];
      for(j=-NB[2];j<totalsize[2]+NB[2];j++){
	
	// real dxa(j) grid value as if whole grid here      
	if(j<=totalsize[2]/2-1){ // if in -j half
		if(j>totalsize[2]/2-1-N2U) ftemp=HOR/(SFTYPE)N2U;
	  else ftemp=ccoef*pow(truea,bcoefm-1.0-(SFTYPE)j);
	}
	else if(j>=totalsize[2]/2){
		if(j<totalsize[2]/2+N2U) ftemp=HOR/(SFTYPE)N2U;
	  else ftemp=ccoef*pow(truea,(SFTYPE)(j)-bcoefp);
	}
	
	// can use this below part in general for any grid!  Just make sure above calcs are based upon full grid, not each cpu
#if(USEMPI==0)
	if(j==0) sx[2]=xasofar;
	if(which==0) dx[1][2][j]=ftemp;
	else if(which==1) idx[1][2][j]=ftemp;
#else 
	if( (j>=startpos[2]-NB[2])&&(j<=endpos[2]+NB[2]) ){
	  
	  if(j==startpos[2]) sx[2]=xasofar;
	  if(which==0) dx[1][2][j-startpos[2]]=ftemp;
	  else if(which==1) idx[1][2][j-startpos[2]]=ftemp;
	}
#endif
	if( (j>=0)&&(j<totalsize[2]) ){
	  xasofar+=ftemp;
	}    
      }
    }
    else if(nonunigridx1==7){

      dir=2;
      xasofar=L[1][dir]; // inner total grid position

      // see mathematica file averywind_gridsetup.nb to choose parameters
      // Constants below are true for N=128, x1in=-1/2, x1out=1/2, posaccretor=0
      NN = (FTYPE)totalsize[dir];
      R0 = pos_accretor[dir];
      if(totalsize[dir]==128){
	AA = BB = 84.78019362357591;
	N0 = 64.0;
      }
      else if(totalsize[dir]==64){
	AA = BB = 190.048594960792;
	N0 = 32.0;
      }
      else if(totalsize[dir]==32){
#if(0)
	// C = 2500
	AA = BB = 390.0797453302086;
	N0 = 16.0;
#else
	// C = 320
	AA = BB = 32.21457825571438;
	N0 = 16.0;
#endif
      }	
	
      Rminus = L[1][dir]; // inner edge position (like x1in)
      Rplus = L[1][dir]+L[2][dir]; // outer edge position (like x1out)
	
      dRupper0=(Rplus-R0)*(pow(BB,1.0/(NN-N0))-1.0)/(BB-1.0);
      dRlower0=(R0-Rminus)*(pow(AA,1.0/N0)-1.0)/(AA-1.0);
      
      for(j=-NB[dir];j<totalsize[dir]+NB[dir];j++){
	
	if(j>=N0){
	  // then use dRupper
	  ftemp = dRupper0*pow(BB,(j-N0)/(NN-N0));
	}
	else{
	  // then use dRlower
	  ftemp = dRlower0*pow(AA,(N0-j-1)/N0);
	}
	  
 
#if(USEMPI==0)
	if(j==0) sx[2]=xasofar;
	if(which==0) dx[1][2][j]=ftemp;
	else if(which==1) idx[1][2][j]=ftemp;
#else 
	if( (j>=startpos[2]-NB[2])&&(j<=endpos[2]+NB[2]) ){
	  
	  if(j==startpos[2]) sx[2]=xasofar;
	  if(which==0) dx[1][2][j-startpos[2]]=ftemp;
	  else if(which==1) idx[1][2][j-startpos[2]]=ftemp;
	}
#endif
	if( (j>=0)&&(j<totalsize[2]) ){
	  xasofar+=ftemp;
	}    
      }
    }
    else{
      fprintf(fail_file,"No known X2 grid type specified\n");
      myexit(1);
    }
    // override to make reflective copy work right in nonunigrid if reflective
    if(reflectix2){ // fix the boundary zones to be "reflective", so boundary code is less complicated for field(divB=0 vs just copy)  (really needed for couple corner zones)
      // this will only catch for correct cpus
      if(which==0){
	dx[1][2][-1]=dx[1][2][0];
	dx[1][2][-2]=dx[1][2][1];
      }
      if(which==1){
	idx[1][2][-1]=idx[1][2][0];
	idx[1][2][-2]=idx[1][2][1];
      }
    }
    if(reflectox2){
      // this will only catch for correct cpus
      if(which==0){
	dx[1][2][N2]=dx[1][2][N2-1];
	dx[1][2][N2+1]=dx[1][2][N2-2];
      }
      if(which==1){
	idx[1][2][N2]=idx[1][2][N2-1];
	idx[1][2][N2+1]=idx[1][2][N2-2];
      }
    }






    
    // X3

    fprintf(log_file,"nonunigridx3: %d totalsize[3]: %d L[2][3]: %21.15g\n",nonunigridx3,totalsize[3],L[2][3]);  fflush(log_file);



    if(which==1){
      iL[1][3]=L[1][3];
      iL[2][3]=L[2][3];
    }
          
    if(nonunigridx3==0){  
      xasofar=L[1][3];
      for(k=-NB[3];k<totalsize[3]+NB[3];k++){
	
	ftemp=L[2][3]/((SFTYPE)(totalsize[3]));
	
#if(USEMPI==0)
	if(k==0) sx[3]=xasofar;
	if(which==0) dx[1][3][k]=ftemp;
	else if(which==1) idx[1][3][k]=ftemp;
#else
	if( (k>=startpos[3]-NB[3])&&(k<=endpos[3]+NB[3]) ){
	  
	  if(k==startpos[3]) sx[3]=xasofar;
	  if(which==0) dx[1][3][k-startpos[3]]=ftemp;
	  else if(which==1) idx[1][3][k-startpos[3]]=ftemp;
	  
	}
#endif
	if( (k>=0)&&(k<totalsize[3]) ){
	  xasofar+=ftemp;
	}
      }
    }
    else if(nonunigridx1==7){


      dir=3;
      xasofar=L[1][dir]; // inner total grid position

      // see mathematica file averywind_gridsetup.nb to choose parameters
      // Constants below are true for N=128, x1in=-1/2, x1out=1/2, posaccretor=0
      NN = (FTYPE)totalsize[dir];
      R0 = pos_accretor[dir];
      if(totalsize[dir]==128){
	AA = BB = 84.78019362357591;
	N0 = 64.0;
      }
      else if(totalsize[dir]==64){
	AA = BB = 190.048594960792;
	N0 = 32.0;
      }
      else if(totalsize[dir]==32){
#if(0)
	// C = 2500
	AA = BB = 390.0797453302086;
	N0 = 16.0;
#else
	// C = 320
	AA = BB = 32.21457825571438;
	N0 = 16.0;
#endif
      }	
	
	
      Rminus = L[1][dir]; // inner edge position (like x1in)
      Rplus = L[1][dir]+L[2][dir]; // outer edge position (like x1out)
	
      dRupper0=(Rplus-R0)*(pow(BB,1.0/(NN-N0))-1.0)/(BB-1.0);
      dRlower0=(R0-Rminus)*(pow(AA,1.0/N0)-1.0)/(AA-1.0);
      
      for(k=-NB[dir];k<totalsize[dir]+NB[dir];k++){
	
	if(k>=N0){
	  // then use dRupper
	  ftemp = dRupper0*pow(BB,(k-N0)/(NN-N0));
	}
	else{
	  // then use dRlower
	  ftemp = dRlower0*pow(AA,(N0-k-1)/N0);
	}
	  
 
#if(USEMPI==0)
	if(k==0) sx[3]=xasofar;
	if(which==0) dx[1][3][k]=ftemp;
	else if(which==1) idx[1][3][k]=ftemp;
#else
	if( (k>=startpos[3]-NB[3])&&(k<=endpos[3]+NB[3]) ){
	  
	  if(k==startpos[3]) sx[3]=xasofar;
	  if(which==0) dx[1][3][k-startpos[3]]=ftemp;
	  else if(which==1) idx[1][3][k-startpos[3]]=ftemp;
	  
	}
#endif
	if( (k>=0)&&(k<totalsize[3]) ){
	  xasofar+=ftemp;
	}
      }
    }
    else{
      fprintf(fail_file,"No known X3 grid type specified\n");
      myexit(1);
    }

      
    // override to make reflective copy work right in nonunigrid if reflective
    if(reflectix3){ // fix the boundary zones to be "reflective", so boundary code is less complicated for field(divB=0 vs just copy)  (really needed for couple corner zones)
      // this will only catch for correct cpus
      if(which==0){
	dx[1][3][-1]=dx[1][3][0];
	dx[1][3][-2]=dx[1][3][1];
      }
      if(which==1){
	idx[1][3][-1]=idx[1][3][0];
	idx[1][3][-2]=idx[1][3][1];
      }
    }
    if(reflectox3){
      // this will only catch for correct cpus
      if(which==0){
	dx[1][3][N3]=dx[1][3][N3-1];
	dx[1][2][N3+1]=dx[1][2][N3-2];
      }
      if(which==1){
	idx[1][3][N3]=idx[1][3][N3-1];
	idx[1][3][N3+1]=idx[1][3][N3-2];
      }
    }
    
  }// end if normal grid or interp wants native grid type
  else{ // use non-native grid type for interpolated output NONNINTERP
    if(outtype==1){  // CHANGE OUTPUT GEOMETRY HERE when converting from cart to spc
      // L[1][1] is input geom: Rin, L[2][1] is input geom: Rout-Rin.
      frac=1.0;
      if(COORD==3){ // z vs x, y=0 plane
	startix[1]=0.0;
	startix[2]=-frac*(L[1][1]+L[2][1]);
	startix[3]=L[1][3];
	iL[1][1]=startix[1];
	//iL[2][1]=2.0*frac*(L[1][1]+L[2][1]); // choose for square output for square image size
	iL[2][1]=frac*(L[1][1]+L[2][1]); // choose for rect/square aspect ratio with rect image size
	iL[1][2]=startix[2];
	iL[2][2]=2.0*frac*(L[1][1]+L[2][1]);
	iL[1][3]=startix[3];
	iL[2][3]=2.0*M_PI;
	dx1=iL[2][1]/N[1];
	dx2=iL[2][2]/N[2];
	dx3=iL[2][3]/N[3];
      }
      else{
	startix[1]=L[1][1];
	startix[2]=L[1][2];
	startix[3]=L[1][3];
	iL[1][1]=startix[1];
	iL[2][1]=L[2][1];
	iL[1][2]=startix[2];
	iL[2][2]=L[2][2];
	iL[1][3]=startix[3];
	iL[2][3]=L[2][3];
	dx1=iL[2][1]/N[1];
	dx2=iL[2][2]/N[2];
	dx3=iL[2][3]/N[3];
      }
    }
    else{
      fprintf(fail_file,"no outtype: %d here\n",outtype);
      myexit(1);
    }// if outtype==1
    for(i=-NB[1];i<N[1]+NB[1];i++){
      idx[1][1][i]=dx1;
    }
    for(j=-NB[2];j<N[2]+NB[2];j++){
      idx[1][2][j]=dx2;
    }
    for(k=-NB[3];k<N[3]+NB[3];k++){
      idx[1][3][k]=dx3;
    }
  }// if other outtypes for interp grid




  // now assign b-grid, which merely depends on a-grid except for dxb(-1) which is never used, but set anyways
  if(which==0){
    for(m=1;m<=3;m++){
      for(i=-NB[m]+1;i<N[m]+NB[m];i++){
	dx[2][m][i]=0.5*(dx[1][m][i]+dx[1][m][i-1]);
      }
      if(COMPDIM<=2){
	if(m!=3){
	  dx[2][m][-NB[m]]=dx[2][m][-NB[m]+1]; // for definiteness, really not defined from a-grid
	}
	else{
	  dx[2][m][0]=dx[1][m][0]; // since x3-dir has no boundary zones and has only 1 grid zone for now
	}
      }
      else{
	dx[2][m][-NB[m]]=dx[2][m][-NB[m]+1]; // for definiteness, really not defined from a-grid
      }
    }
  }
  else if(which==1){
    for(m=1;m<=3;m++){
      for(i=-NB[m]+1;i<N[m]+NB[m];i++){
	idx[2][m][i]=0.5*(idx[1][m][i]+idx[1][m][i-1]);
      }
      if(COMPDIM<=2){
	if(m!=3){
	  idx[2][m][-NB[m]]=idx[2][m][-NB[m]+1]; // for definiteness, really not defined from a-grid
	}
	else{
	  idx[2][m][0]=idx[1][m][0]; // since x3-dir has no boundary zones and has only 1 grid zone for now
	}
      }
      else{
	idx[2][m][-NB[m]]=idx[2][m][-NB[m]+1]; // for definiteness, really not defined from a-grid
      }
    }
  }

  // for outgparm
  ix1in=iL[1][1];
  ix2in=iL[1][2];
  ix3in=iL[1][3];
  ix1out=iL[1][1]+iL[2][1];
  ix2out=iL[1][2]+iL[2][2];
  ix3out=iL[1][3]+iL[2][3];

  fprintf(log_file,"end: init_dx\n"); fflush(log_file);
  return(0);
}


#define RGMOVE (.001) // how far to move beyond rgp if start off below rgp
#define STARTF (1.0/4.0) // what factor or less 2nd condition below should be

/* Once dxs fixed, below applies to any grid */
/* BEGIN: Setup x1,x2,x3 positions from dx data */
/* I break this loop from the previous because x depends on diffs in dx */
/* That is to say I do not have to worry about order problems */
/* Note this data would be dynamic if dynamic grid */
// below accumulates small error, could set up grid more symmetrically or define x exactly like dx as function of grid position
// 1/2/3-d valid
int init_x(int n1,int n2,int n3,int n1b,int n2b,int n3b, SFTYPE*sx, int which, int outtype)
{
  int i,j,k,l,m;
  SFTYPE startbc[3+1]={0};
  SFTYPE startdiff=0;
  int N[3+1],NB[3+1];
  SFTYPE ftemp;

  fprintf(log_file,"begin: init_x ... "); fflush(log_file);

  N[1]=n1;
  N[2]=n2;
  N[3]=n3;
  NB[1]=n1b;
  NB[2]=n2b;
  NB[3]=n3b;

  // determine first BC's starting a-grid position
  startbc[1]=sx[1];
  if((which==0)||(outtype==0)){

    for(i=-1;i>=-NB[1];i--){
      if(which==0) startbc[1]-=dx[1][1][i];
      else if(which==1) startbc[1]-=idx[1][1][i];
    }
    startdiff=startbc[1]; // starting value


    //////
    //
    // see if need to correct starting position
    //
    if((COORD==3)&&( (reflectix1==0)||((reflectix1==1)&&(L[1][1]>1.0E-3)) ) ){
      if(startbc[1]<=rgp){
	if(myid<=0){
	  fprintf(logfull_file,"r coordinate assigned at boundary value to less than rg!\n");
	}
	if(nonunigridx1==4){
	  fprintf(fail_file,"This grid requires sx[1] to be set and not changed\n");
	  myexit(1);
	}
	else{
	  if(myid<=0){
	    fprintf(logfull_file,"... moving from: %21.10g to ",startbc[1]);
	  }
	  // so force to be good
	  if(which==0) startbc[1]=dx[2][1][-NB[1]]/(STARTF*.9999)+rgp;
	  else if(which==1) startbc[1]=idx[2][1][-NB[1]]/(STARTF*.9999)+rgp;
	  //sx[1]+=(rgp-startbc[1]+RGMOVE);
	  //startbc[1]=sx[1];
	  //for(i=-1;i>=-NB[1];i--){
	  //if(which==0) startbc[1]-=dx[1][1][i];
	  //else if(which==1) startbc[1]-=idx[1][1][i];
	  //}
	  if(startbc[1]<=rgp){
	    fprintf(fail_file,"Failed to correct starting position(1)\n");
	    myexit(1);
	  }
	  else{
	    if(myid<=0){
	      fprintf(logfull_file," %21.10g\n",startbc[1]);
	    }
	  }
	  if(myid<=0){
	    fprintf(logfull_file,"Can only have negative (r-rgp) with reflecting BC and rgp=0\n");
	  }
	}
      } 
      // now check more conditions (zone i=0 used since bzones may be artificially low and)
      if(which==0) ftemp=dx[2][1][0]/(startbc[1]-rgp);
      else if(which==1) ftemp=idx[2][1][0]/(startbc[1]-rgp);
      else{
        fprintf(fail_file,"no which: %d here\n",which);
	myexit(1);
      }
      if(myid<=0){
	if(ftemp>STARTF){
	  fprintf(logfull_file,"Bad R=dr/(r-rgp) inner boundary ratio: %15.10g\n",ftemp);
	  fprintf(logfull_file,"This bad R has known to produce instabilities with gam=5/3\n"); 
	}
	else{
	  fprintf(logfull_file,"Good R=dr[0]/(r[0]-rgp) inner boundary ratio: %15.10g\n",ftemp);
	}
      }
    }

    // now transmit new startbc[1] to all cpus, assuming cpu=0 holds it if any do
#if(USEMPI)
    ftemp=startbc[1]-startdiff;
    MPI_Bcast(&ftemp,1,MPI_SFTYPE,0,MPI_COMM_WORLD);
    if(mycpupos[1]!=0) startbc[1]+=ftemp; // otherwise no need to correct under normal circumstances
#endif
    // now get comp grid start value
    sx[1]=startbc[1];
    for(i=-1;i>=-NB[1];i--){
      if(which==0) sx[1]+=dx[1][1][i];
      else if(which==1) sx[1]+=idx[1][1][i];
    }
    // only if no CPU chop at different fixed radii
    //    L[1][1]=sx[1];
  }
    
    
    
  startbc[2]=sx[2];
  if((which==0)||(outtype==0)){
    for(j=-1;j>=-NB[2];j--){
      if(which==0) startbc[2]-=dx[1][2][j];
      else if(which==1) startbc[2]-=idx[1][2][j];
    }
    startdiff=startbc[2]; // starting value
    /*  
	if(COORD==2){
	if(!reflectix2){
	if(startbc[2]<0.0){
	fprintf(log_file,"r coordinate assigned at boundary value to less than zero! ... moving from: %21.10g to ",startbc[2]);
	startbc[2]=0.0;
	fprintf(log_file,"%21.10g\n",startbc[2]);
	}
	}
	}
	if(COORD==3){
	if(startbc[2]<0.0){
	fprintf(log_file,"theta coordinate assigned at boundary value to less than zero! ... moving from: %21.10g to ",startbc[2]);
	startbc[2]=0.0;
	fprintf(log_file,"%21.10g\n",startbc[2]);
	}
	}
    */
#if(USEMPI)
    ftemp=startbc[2]-startdiff;
    MPI_Bcast(&ftemp,1,MPI_SFTYPE,0,MPI_COMM_WORLD);
    if(mycpupos[2]!=0) startbc[2]+=ftemp; // otherwise no need to correct under normal circumstances
#endif
    // now get comp grid start value
    sx[2]=startbc[2];
    for(j=-1;j>=-NB[2];j--){
      if(which==0) sx[2]+=dx[1][2][j];
      else if(which==1) sx[2]+=idx[1][2][j];
    }

  }  
  startbc[3]=sx[3];

  if((which==0)||(outtype==0)){
    for(k=-1;k>=-NB[3];k--){
      if(which==0) startbc[3]-=dx[1][3][k];
      else if(which==1) startbc[3]-=idx[1][3][k];
    }
    startdiff=startbc[3]; // starting value
    /*  
	if(COORD>=2){
	if(startbc[3]<0.0){
	fprintf(log_file,"phi coordinate assigned at boundary value out of normal 0-2*Pi domain: %21.10g: no stop condition\n",startbc[3]);
	}
	
    */
    for(l=1;l<=3;l++){
      fprintf(log_file,"sx[%d]: %21.10g startbc[%d]: %21.10g\n",l,sx[l],l,startbc[l]);
    }

#if(USEMPI)
    ftemp=startbc[3]-startdiff;
    MPI_Bcast(&ftemp,1,MPI_SFTYPE,0,MPI_COMM_WORLD);
    if(mycpupos[3]!=0) startbc[3]+=ftemp; // otherwise no need to correct under normal circumstances
#endif
    // now get comp grid start value
    sx[3]=startbc[3];
    for(k=-1;k>=-NB[3];k--){
      if(which==0) sx[3]+=dx[1][3][k];
      else if(which==1) sx[3]+=idx[1][3][k];
    }

  }
  // myexit(1);
  
  if(which==0){
    // allow r<0 so volume factor is correct?
    for(m=1;m<=3;m++){
      for(i=-NB[m];i<N[m]+NB[m];i++){    
	if(i==-NB[m]){
	  x[1][m][i]=startbc[m]; /* a-grid */
	  x[2][m][i]=x[1][m][i]+dx[1][m][i]*0.5; /* b-grid */
	}
	/* Applying rule for a-grid: x(i)=x(i-1)+dx(i-1) */
	/* Applying rule for b-grid: x(i)=x(i-1)+dx(i) */
	else{
	  x[1][m][i]=x[1][m][i-1]+dx[1][m][i-1];
	  x[2][m][i]=x[2][m][i-1]+dx[2][m][i];
	}
      }
    } /* END: Setup x1,x2,x3 positions from dx data */
  }
  else if(which==1){
    // allow r<0 so volume factor is correct?
    for(m=1;m<=3;m++){
      for(i=-NB[m];i<N[m]+NB[m];i++){    
	if(i==-NB[m]){
	  ix[1][m][i]=startbc[m]; /* a-grid */
	  ix[2][m][i]=ix[1][m][i]+idx[1][m][i]*0.5; /* b-grid */
	}
	/* Applying rule for a-grid: x(i)=x(i-1)+dx(i-1) */
	/* Applying rule for b-grid: x(i)=x(i-1)+dx(i) */
	else{
	  ix[1][m][i]=ix[1][m][i-1]+idx[1][m][i-1];
	  ix[2][m][i]=ix[2][m][i-1]+idx[2][m][i];
	}
      }
    } /* END: Setup x1,x2,x3 positions from dx data */

  }

  fprintf(log_file,"end: init_x\n"); fflush(log_file);  
  return(0);
}



/* Currently assume g1=1 in rest of code, h=high, l=low */
#if(COORD==1)
/* x1=x, x2=y, x3=z */
#define g1(x1) 1
#define dg1(x1) 0  // always zero currently, not used currently
#define g2(x1) 1
#define dg2(x1) 0  // par(g2i,dx1)
#define g31(x1) 1
#define dg31(x1) 0  // par(g31,dx1)
#define g32(x2) 1
#define dg32(x2) 0  // par(g32,dx2)

#define dvl1(x1h,x1l) (x1h-x1l)
#define dvl2(x2h,x2l) (x2h-x2l) 
#define dvl3(x2h,x2l) (x2h-x2l)

#endif
#if(COORD==2)
/* x1=z, x2=r, x3=phi */
#define g1(x1) 1
#define dg1(x1) 0
#define g2(x1) 1
#define dg2(x1) 0
#define g31(x1) 1
#define dg31(x1) 0
#define g32(x2) fabs(x2)
#define dg32(x2) 1

#define dvl1(x1h,x1l) (x1h-x1l)
#define dvl2(x2h,x2l) (0.5*(x2h*x2h-x2l*x2l)) 
#define dvl3(x3h,x3l) (x3h-x3l)

#endif
#if(COORD==3)
/* x1=r, x2=theta, x3=phi */
// by making metric >0 always one then deals with primitive variables properly for given geometry.  This method allows access to primitive variables properly and any metric with primitive variable properly.
#define g1(x1) 1
#define dg1(x1) 0
#define g2(x1) fabs(x1)
#define dg2(x1) 1
#define g31(x1) fabs(x1)
#define dg31(x1) 1
#define g32(x2) fabs(sin(x2))
#define dg32(x2) ( (x2<0) ? -cos(x2) : ( (x2>M_PI) ? -cos(x2) : ( ((x2==0.0)||(x2==M_PI)) ? 0.0 : cos(x2)) ) )
// then dg32/g32 is exactly 1/tan(x2), as needed for step_pgc's curvature term

// for volume terms, only require are positive and correct size for across boundaries (like polar axis)
#define dvl1(x1h,x1l) (THIRD*(x1h*x1h*x1h-x1l*x1l*x1l))
#define dvl2(x2h,x2l) (-(cos(x2h)-cos(x2l)))
#define dvl3(x3h,x3l) (x3h-x3l)

#endif

// don't really need dvl3 at all, even in 3D


// 1/2/3-d valid


int init_diffs(void)
{

  int i,j,k,l,m;
  SFTYPE xm,xml,xmh;
  char temps[100];

  int N[3+1],NB[3+1];

  fprintf(log_file,"begin: init_diffs ... "); fflush(log_file);

  N[1]=N1;
  N[2]=N2;
  N[3]=N3;
  NB[1]=N1BND;
  NB[2]=N2BND;
  NB[3]=N3BND;

  /* BEGIN: assign values for matric scales and volume factors */
  /* I break this loop from the previous because metric may need diffs in x */
  /* That is to say I do not have to worry about order problems */
  /* Note this data would be dynamic if dynamic grid */


  // Not general!! --must modify below for different coords than given above
  for(l=1;l<=NUMGRID;l++){
    for(m=1;m<=2;m++){ // no need for 3
      for(i=-NB[m];i<N[m]+NB[m];i++){    
	/* Setup metric */
	xm=x[l][m][i];

	if(m==1){
	  g[l][1][i]=g1(xm);
	  dg[l][1][i]=dg1(xm);
	  g[l][2][i]=g2(xm);
	  dg[l][2][i]=dg2(xm);
	  g[l][3][i]=g31(xm);
	  dg[l][3][i]=dg31(xm);
	}
	else if(m==2){
	  g[l][4][i]=g32(xm);	  
	  dg[l][4][i]=dg32(xm);
	}
      }
    }
  }
  /* END: assign values for matric scales and volume factors */

  /* Assign volume differences */
  for(l=1;l<=NUMGRID;l++){
    for(m=1;m<=3;m++){
      for(i=-NB[m];i<N[m]+NB[m];i++){
	if(l==1){
	  if(i==N[m]+NB[m]-1){
	    xml=x[l][m][i];
	    xmh=xml+dx[l][m][i];
	  }
	  else{
	    xml=x[l][m][i];
	    xmh=x[l][m][i+1];
	  }
	}
	else if(l==2){
	  if(i==-NB[m]){
	    xmh=x[l][m][i];
	    xml=xmh-dx[l][m][i];
	  }
	  else{
	    xml=x[l][m][i-1];
	    xmh=x[l][m][i];
	  }
	}
	if(m==1)      dvl[l][m][i]=fabs(dvl1(xmh,xml)); //fabs for boundary zones.  Shouldn't have negative volumes
	else if(m==2) dvl[l][m][i]=fabs(dvl2(xmh,xml));
	else if(m==3) dvl[l][m][i]=fabs(dvl3(xmh,xml));
      }
      
    }
  }
  //for(i=0;i<=N2/2-1;i++){
  //  printf("%15.10g %15.10g %15.10g\n",dvl[1][2][i],dvl[1][2][N2-i-1],dvl[1][2][i]-dvl[1][2][N2-i-1]);
  //}
  //fflush(stdout);

  // correct for volume diff across theta=0,Pi boundary on b-grid when reflecting...r is already fine
  // below needed for x1-dir vy-adv at theta=0,theta=Pi
  // but not really used for v, just a place holder for advecting vy in x-dir.  Should really take care of there instead.

  // volume at theta edges are nonzero!
  // if used, general makes code easier since don't have to worry about 0/0 for reflecting bc
  // general coordinate singularity fix
  periodicx2special=0; // initialize
  if((COORD==3)&&(L[1][2]<1.0E-3) ){
    if(mycpupos[2]==0){
      dvl[2][2][0]=fabs(dvl2(x[2][2][0],0.0))+fabs(dvl2(0,x[2][2][-1]));
      g[1][4][0]=0; // forced so singuarlity checks will pick up
      if(periodicx2) periodicx2special=1; // tells bound() to use special routine
    }
  }
  if((COORD==3)&&(L[1][2]+L[2][2]>M_PI-1.0E-3) ){
    if(mycpupos[2]==ncpux2-1){
      dvl[2][2][N2]=fabs(dvl2(M_PI,x[2][2][N2-1]))+fabs(dvl2(x[2][2][N2],M_PI));
      g[1][4][N2]=0; // forced so singuarlity checks will pick up
      if(periodicx2) periodicx2special=1; // tells bound() to use special routine
    }
  }

  if(numprocs>1){
    // check to see whether MPI is setup right
    if(periodicx2special){
      if(! ( (ncpux3==1)||(ncpux3%2==0) ) ){ // ncpux3 must be 1 or even
	// then no general setup for this case, so fail
	fprintf(fail_file,"No setup for non-1 and odd cpu number in X3 direction when periodicx2special=%d, ncpux3=%d\n",periodicx2special,ncpux3);
	myexit(10);
      }
    }
  }
  if(periodicx2special){
    if(! ( (N3==1)||(N3%2==0))){
      fprintf(fail_file,"With periodicx2special must have N3==1 or even\n");
      myexit(10);
    }
  }

  fprintf(log_file,"end: init_diffs\n"); fflush(log_file);
  return(0);
}



// this forces code to situate the 1-zone dimensions to have all components zone centered in that direction
int init_reduceddimension(void)
{
  int i,j,k,l,m;

  int N[3+1],NB[3+1];

  fprintf(log_file,"begin: init_reduceddimension ... "); fflush(log_file);

  N[1]=N1;
  N[2]=N2;
  N[3]=N3;
  NB[1]=N1BND;
  NB[2]=N2BND;
  NB[3]=N3BND;
  
  // now correct for reduced dimension cases(compsaves will be correct from these)
  for(m=1;m<=3;m++){
    for(i=-NB[m];i<N[m]+NB[m];i++){    
      if(N[m]==1){
	// generic terms
	dx[2][m][i]=dx[1][m][i];
	if(LOWMEMMODE==0) dvl[2][m][i]=dvl[1][m][i];

	x[1][m][i]=x[2][m][i];
	// special directional terms
	if(m==1){
	  if(LOWMEMMODE==0) g[1][1][i]=g[2][1][i];
	  if(LOWMEMMODE==0) g[1][2][i]=g[2][2][i];
	  if(LOWMEMMODE==0) g[1][3][i]=g[2][3][i];
	  if(LOWMEMMODE==0) dg[1][1][i]=dg[2][1][i];
	  if(LOWMEMMODE==0) dg[1][2][i]=dg[2][2][i];
	  if(LOWMEMMODE==0) dg[1][3][i]=dg[2][3][i];
	}
	if(m==2){
	  if(LOWMEMMODE==0) g[1][4][i]=g[2][4][i];
	  if(LOWMEMMODE==0) dg[1][4][i]=dg[2][4][i];
	}
	if(m==3){
	  // no metric terms
	}
      }
    }
  } /* END: Setup x1,x2,x3 positions from dx data */
  
  // do interpolated part if POSTPROC==1
  // interpolation only requires x and dx
  if(POSTPROC){
    for(m=1;m<=3;m++){
      for(i=-NB[m];i<N[m]+NB[m];i++){    
	if(N[m]==1){
	  idx[2][m][i]=idx[1][m][i];
	  ix[1][m][i]=ix[2][m][i];	  
	}
      }
    } /* END: Setup x1,x2,x3 positions from dx data */
  }

  // force central position if doing reduced COORD=3 N2=1 case:
  if((COORD==3)&&(N2==1)&&(fabs(x[2][2][0]-M_PI*0.5)<1E-3)){ // then assume at equator and force central value to be exactly such that cos(theta)=0
    x[1][2][0]=x[2][2][0]=M_PI*0.5;
    if(LOWMEMMODE==0) g[2][4][0]=g[1][4][0]=1.0;
    if(LOWMEMMODE==0) dg[2][4][0]=dg[1][4][0]=0.0;
  }


  fprintf(log_file,"end: init_reduceddimension\n"); fflush(log_file);
  return(0);
}


void init_compsave(void)
{
  int i,j,k;

  fprintf(log_file,"begin: init_compsave ... "); fflush(log_file);

  // ds[type][direction]...
  LOOPF{
    
    
    // really * dx[1/2][3][k], but don't want to eat that much memory just for 1 multiply
    
    // .5,.5,.5 (used with centered variables)
    ds[1][1][k][j][i]=g[1][2][i]*g[1][3][i]*dvl[1][2][j]; // *dx[2][3][k];
    ds[1][2][k][j][i]=g[2][3][i]*dx[1][1][i]*g[1][4][j]; // *dx[2][3][k];
    ds[1][3][k][j][i]=g[2][2][i]*dx[1][1][i]*dx[1][2][j];
    
    // 0,.5,.5(eg. vx)
    ds[2][1][k][j][i]=g[2][2][i]*g[2][3][i]*dvl[1][2][j]; // *dx[2][3][k];
    ds[2][2][k][j][i]=g[1][3][i]*dx[2][1][i]*g[1][4][j]; // *dx[2][3][k];
    ds[2][3][k][j][i]=g[1][2][i]*dx[2][1][i]*dx[1][2][j];
    
    // .5,0,.5(eg. vy)
    ds[3][1][k][j][i]=g[1][2][i]*g[1][3][i]*dvl[2][2][j]; // *dx[2][3][k];
    ds[3][2][k][j][i]=g[2][3][i]*dx[1][1][i]*g[2][4][j]; // *dx[2][3][k];
    ds[3][3][k][j][i]=g[2][2][i]*dx[1][1][i]*dx[2][2][j];
    
    // 0,0,.5
    ds[4][1][k][j][i]=g[2][2][i]*g[2][3][i]*dvl[2][2][j]; // *dx[2][3][k];
    ds[4][2][k][j][i]=g[1][3][i]*dx[2][1][i]*g[2][4][j]; // *dx[2][3][k];
    ds[4][3][k][j][i]=g[1][2][i]*dx[2][1][i]*dx[2][2][j];
    
    //oarcl[type][dir]...
    // assumes divide of dx3 in oarcl[][3] term(so don't have to store a/b version in x3 dir
    
    // .5,.5,.5 (used with centered variables)
    // length along dir direction, used for the given variable differences
    // (i.e. centered BETWEEN such a variable)
    oarcl[1][1][k][j][i]=1.0/dx[2][1][i] ;
    oarcl[1][2][k][j][i]=1.0/(g[2][2][i]*dx[2][2][j]) ;
    oarcl[1][3][k][j][i]=1.0/(g[2][3][i]*g[2][4][j]) ;// 1/dx[2][3][k]
    
    // 0,.5,.5(eg. vx, i.e. used between 2 vx's)
    oarcl[2][1][k][j][i]=1.0/dx[1][1][i] ;
    if( (COORD==3)&&(reflectix1==1)&&(L[1][1]<1.0E-3)&&(i==0) ){
      oarcl[2][2][k][j][i]=1.0; // set to this to avoid singularity ( this one is no big deal if COMPDIM<=2)
    }
    else{
      oarcl[2][2][k][j][i]=1.0/(g[1][2][i]*dx[2][2][j]) ; // real coord sing if COMPDIM==3 and r=0
    }
    oarcl[2][3][k][j][i]=1.0/(g[1][3][i]*g[2][4][j]) ; // 1/dx[2][3][k]
    
    // .5,0,.5(eg. vy)
    oarcl[3][1][k][j][i]=1.0/dx[2][1][i] ;
    oarcl[3][2][k][j][i]=1.0/(g[2][2][i]*dx[1][2][j]) ;
    //    if((COORD==3)&&( ((reflectix2==1)&&(j==0))||((reflectox2==1)&&(j==N2)) )){
    if(g[1][4][j]==0.0){
      oarcl[3][3][k][j][i]=0.0;
    }
    else{
      oarcl[3][3][k][j][i]=1.0/(g[2][3][i]*g[1][4][j]) ; // /dx[2][3][k] // real coord sing if COMPDIM==3  singularity if theta=0,Pi
    }
    // 0,0,.5(eg. sigma_{r,theta})
    oarcl[4][1][k][j][i]=1.0/dx[1][1][i] ;
    if( (COORD==3)&&(reflectix1==1)&&(L[1][1]<1.0E-3)&&(i==0) ){
      oarcl[4][2][k][j][i]=1.0 ; // to avoid sing
    }
    else{
      oarcl[4][2][k][j][i]=1.0/(g[1][2][i]*dx[1][2][j]) ;  // real coord sing if COMPDIM==3 and r=0
    }
    //    if((COORD==3)&&( ((reflectix2==1)&&(j==0))||((reflectox2==1)&&(j==N2)) )){
    if(g[1][4][j]==0.0){
      oarcl[4][3][k][j][i]=0.0;
    }
    else{
      oarcl[4][3][k][j][i]=1.0/(g[1][3][i]*g[1][4][j]) ;// /dx[1][3][k] // real coord sing if COMPDIM==3  singularity if theta=0,Pi
    }      
    
    //ovol[type]...
    
    // assumes divided by dx3 term in code(FYI dx3=dvl3)
    
    // .5,.5,.5(eg. rho)
    ovol[1][k][j][i]=1.0/(dvl[1][1][i]*dvl[1][2][j]); // /dvl[1][3][k]
    // .5,0,.5(eg. vy)
    ovol[3][k][j][i]=1.0/(dvl[1][1][i]*dvl[2][2][j]) ; // /dvl[1][3][k]
    // 0,.5,.5(eg. vx)
    ovol[2][k][j][i]=1.0/(dvl[2][1][i]*dvl[1][2][j]); // /dvl[1][3][k]
    
    // 0,0,.5
    ovol[4][k][j][i]=1.0/(dvl[2][1][i]*dvl[2][2][j]) ; // /dvl[1][3][k]
    
    // redundant loop, but who cares.
    odx[1][1][i]=1.0/dx[1][1][i];
    odx[2][1][i]=1.0/dx[2][1][i];
    
    odx[1][2][j]=1.0/dx[1][2][j];
    odx[2][2][j]=1.0/dx[2][2][j];

    // below needed for DS and oarcl and OVOL in 3D
    odx[1][3][k]=1.0/dx[1][3][k];
    odx[2][3][k]=1.0/dx[2][3][k];

    odvl[1][1][i]=1.0/dvl[1][1][i];
    odvl[2][1][i]=1.0/dvl[2][1][i];

    odvl[1][2][j]=1.0/dvl[1][2][j];
    odvl[2][2][j]=1.0/dvl[2][2][j];

    odvl[1][3][k]=1.0/dvl[1][3][k];
    odvl[2][3][k]=1.0/dvl[2][3][k];

    og[1][1][i]=1.0/g[1][1][i];
    og[2][1][i]=1.0/g[2][1][i];

    og[1][2][i]=1.0/g[1][2][i];
    og[2][2][i]=1.0/g[2][2][i];

    og[1][3][i]=1.0/g[1][3][i];
    og[2][3][i]=1.0/g[2][3][i];

    if(g[1][4][j]==0.0){
      og[1][4][j]=0.0; // forced since when used it corrects singularity in correct way(since other terms should give 0)
    }
    else og[1][4][j]=1.0/g[1][4][j];

    og[2][4][j]=1.0/g[2][4][j];
    
  }
  fprintf(log_file,"end: init_compsave\n"); fflush(log_file);

}

// 1/2/3-d valid
int init_data(void)
{

  int i,j,k,l,m;
  SFTYPE radsize2,tempf1,tempf2,xz0;
  SFTYPE lkep,Cconst;
  SFTYPE ftemp,postemp;
  int blob[2][3];
  int firstblob=1;


  fprintf(log_file,"begin: init_data ... "); fflush(log_file);
  /* BEGIN: zero out(assign) data */

  /* I break this loop from the previous because sca/vev/metric may need diffs in x */
  /* That is to say I do not have to worry about order problems */
  /* Note this data would be dynamic if dynamic grid */


  if(analoutput!=0){
  
    // this is normal call to analsolve
    analsolve(POSTPROC);
#if(RAD)
    initrad();
#endif
    for(l=1;l<=NUMSCA;l++){
      LOOPF{
        s[l][k][j][i] = sanal[l][k][j][i] ;
      }
    }


    if(POSTPROC==0){
      tdep_compute(); // compute time dep stuff
    }

    for(l=1;l<=NUMSCA;l++){
      LOOPF{
	s[l][k][j][i] = sanal[l][k][j][i] ;
      }
    }
    
    for(m=1;m<=NUMVEC;m++){
      for(l=1;l<=3;l++){
	LOOPF{
	  v[m][l][k][j][i] = vanal[m][l][k][j][i] ;
	}
      }
    }
#if(POSTPROC==0)
    floor_correct(1,-1);
    floor_correct(2,-1);
#endif

  }
  else{
    //POTENTIAL
    LOOPF{
      s[3][k][j][i] = -GRAVC*MASSBH/(x[2][1][i]-rgp) ;
      //s[3][k][j][i] = 0.0 ;
    }
    
    
    //MASS DENSITY    
    LOOPF{
      s[1][k][j][i] = 1.0;
    }
    
    //INTERNAL ENERGY DENSITY
    LOOPF{
      if(press==1){
	if(wgam){
	  cs=pow(s[1][k][j][i],0.5*(gam-1.)); // usually not needed here
	  s[2][k][j][i] = pow(s[1][k][j][i],gam)/gam/(gam-1.);
	}
	else{
	  s[2][k][j][i] = cs*cs*s[1][k][j][i] ;
	}
      }
      else{
	s[2][k][j][i]=0;
      }
    }
    
    //VECTORS(v and B)
    for(m=1;m<=3;m++){  
      LOOPF{
	v[1][m][k][j][i]=0.0;
      }
    }
    
    for(m=1;m<=3;m++){
      LOOPF{
	v[2][m][k][j][i]=0.0;
      }
    }
  }
  /* END: zero out(assign) data */
  fprintf(log_file,"end: init_data\n"); fflush(log_file);
  return(0);
}



void init_tvdlfgrid(void)
{
  // now done with a/b grid stuff.
  // now force all positions to be on b-grid for TVD....don't worry about rest since not used(metric, oarcl,OVOL, etc.)
  int k,j,i;

  fprintf(log_file,"begin: init_tvdlfgrid ... "); fflush(log_file);

  LOOPF{
    // positions as if on b-grid
    x[1][1][i]=x[2][1][i];
    x[1][2][j]=x[2][2][j];
    x[1][3][k]=x[2][3][k];
    //dx as if on b-grid
    dx[1][1][i]=dx[2][1][i];
    dx[1][2][j]=dx[2][2][j];
    dx[1][3][k]=dx[2][3][k];
    // volume should be as on a-grid
    if(LOWMEMMODE==0) dvl[2][1][i]=dvl[1][1][i];
    if(LOWMEMMODE==0) dvl[2][2][j]=dvl[1][2][j];
    if(LOWMEMMODE==0) dvl[2][3][k]=dvl[1][3][k];
  }

  fprintf(log_file,"end: init_tvdlfgrid\n"); fflush(log_file);
}



void init_checks(void)
{
  // general checks
  if( (GZIPIMAGE>0)&&(USEGM==1)){
    fprintf(fail_file,"Cannot use gzip and use gm\n");
    fprintf(fail_file,"Use external gzipper watcher program\n");
    myexit(1);
  }

  if( (POSTPROC>1)||(POSTPROC<0)){
    fprintf(fail_file,"invalid POSTPROC: %d\n",POSTPROC);
    myexit(1);
  }

  if( (PRODUCTION>1)||(PRODUCTION<0)){
    fprintf(fail_file,"invalid PRODUCTION: %d\n",PRODUCTION);
    myexit(1);
  }

  if( (LINUXCLUSTER>1)||(LINUXCLUSTER<0)){
    fprintf(fail_file,"invalid LINUXCLUSTER: %d\n",LINUXCLUSTER);
    myexit(1);
  }

  if( (PERFTEST>1)||(PERFTEST<0)){
    fprintf(fail_file,"invalid PERFTEST: %d\n",PERFTEST);
    myexit(1);
  }

  if( (USEMPI>1)||(USEMPI<0)){
    fprintf(fail_file,"invalid USEMPI: %d\n",USEMPI);
    myexit(1);
  }

  if( (USEGM>1)||(USEGM<0)){
    fprintf(fail_file,"invalid USEGM: %d\n",USEGM);
    myexit(1);
  }

  if( (cpugeompick>1)||(cpugeompick<0)){
    fprintf(fail_file,"invalid cpugeompick: %d\n",cpugeompick);
    myexit(1);
  }

  if( (N3>1)&&(COMPDIM==2)){
    fprintf(fail_file,"invalid COMPDIM: %d for N3=%d \n",COMPDIM,N3);
    myexit(1);
  }

  if( (N3==1)&&(COMPDIM==3)){
    fprintf(fail_file,"You probably didn't mean to set N3==1 and COMPDIM==3.  It's inefficient in memory use: COMPDIM: %d N3=%d\n",COMPDIM,N3);
    //myexit(1);
  }

  if((DOLOSSDIAGMEMORY==0)&&(DOLOSSDIAG==1)&&(LOOPTYPE==1)){
    fprintf(fail_file,"You must turn on memory for old loss diagnostics\n");
    myexit(1);
  }

  if((FLOORDUMPFLAGMEMORY==0)&&(FLOORDUMPFLAG==1)){
    fprintf(fail_file,"You must turn on memory for floor dump diagnostic\n");
    myexit(1);
  }


}


void accountstoreset(void)
{
  int i,j,k;
  
#if(LOOPTYPE==1)
  LOOPF{
    if((i>=intix1)&&(i<intox1)&&(j>=intix2)&&(j<intox2)&&(k>=intix3)&&(k<intox3)) accountstore[k][j][i]=1; else accountstore[k][j][i]=0;
  }
#else
  LOOPF{
    accountstore[k][j][i]=0;
  }
  LOOPINT{
    accountstore[k][j][i]=1;
  }
#endif
}





// tests bound call
void boundtest(int which)
{

  int i,j,k;

  if(which==0){
    // setup variables for testing
    LOOPF{
      s[1][k][j][i]=-1-1*100;
      s[2][k][j][i]=-1-2*100;
      s[3][k][j][i]=-1-3*100;

      v[1][1][k][j][i]=-1-(3+0*3+1)*100;
      v[1][2][k][j][i]=-1-(3+0*3+2)*100;
      v[1][3][k][j][i]=-1-(3+0*3+3)*100;

      v[2][1][k][j][i]=-1-(3+1*3+1)*100;
      v[2][2][k][j][i]=-1-(3+1*3+2)*100;
      v[2][3][k][j][i]=-1-(3+1*3+3)*100;
    }
    LOOP{
      s[1][k][j][i]=myid-1*100;
      s[2][k][j][i]=myid-2*100;
      s[3][k][j][i]=myid-3*100;

      v[1][1][k][j][i]=myid-(3+0*3+1)*100;
      v[1][2][k][j][i]=myid-(3+0*3+2)*100;
      v[1][3][k][j][i]=myid-(3+0*3+3)*100;

      v[2][1][k][j][i]=myid-(3+1*3+1)*100;
      v[2][2][k][j][i]=myid-(3+1*3+2)*100;
      v[2][3][k][j][i]=myid-(3+1*3+3)*100;

    }

    dump(NULL,666,DTYPE,0);
    
  }
  if(which==1){
    //    bound(NULL,NULL,-1,-1,123);
#if(USEMPI)
  bound_mpi(NULL,NULL,-1,-1,123);
#endif  
  }

  if(which==2){
    dump(NULL,667,DTYPE,0);
  }

}



void tdep_compute_ppc(void)
{
  // just minimal assignments for POSTPROC setting


  int i,j,k;
  static int firsttime=1;
  static int switchbc=0;
  int error=0;

  SFTYPE ftemp;

  
  if(runtype>0){
    if(t>=0.0) alpha_real=alpha_real0; // go ahead and turn on when using reentrant data as initial data
  }

  

  // res
  if(t>=0.0) resist_real=resist_real0;
  
  cool=0;

  firsttime=0;
}


void analsolve_ppc(int calltype)
{
  int outtype;
  int i,j,k,l,m;
  SFTYPE arot;

  fprintf(log_file,"begin: analsolve_ppc(100) ... "); fflush(log_file);
  
  DYNAMICMM=0; // already read in file usually
  
  if(COORD==3){
    x3in=  (0.0);
    x3out= (2.0*M_PI);
    
    Rinner=1.2*rg;
    Router=10.0*11.5*rg; // H00 first axisym

    x1in=Rinner;
    x1out= Router;
    
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
    
    if(GAMMIEIMAGE){
      // Rinner is innermost edge, NOT including bounary zones
      //Rinner=1.686836; / 256^2 a=0.5 Rout=20
      rgp=0.0; // since sometimes inside event horizon
      Rinner=1.741087; // 512^2 a=0.5 Rout=20
      arot=.75;
      Rinner=0.98*(1. + sqrt(1. - arot*arot));
      nonunigridx1= 5; // 0=unigrid 1+=nonunigrid
      nonunigridx2= 5; // 0=unigrid 1+=nonunigrid
      nonunigridx3= 0;
      x1in=Rinner;
      Router=20;
      Router=40;
      x1out=Router;
    }
    
    
  }
  if(COORD==1){
    Rinner=1.5*rg;
    Router=11.5*rg;
    
    x1in=-Router;
    x1out=Router;
    x2in=-Router;
    x2out=Router;
    
    x3in=-Router;
    x3out=Router;
    
    nonunigridx1= 0; // 0=unigrid 1+=nonunigrid
    nonunigridx2= 0; // 0=unigrid 1+=nonunigrid
    nonunigridx3= 0;
    
  }

  // below is used as basis for converting images to real data and back
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


