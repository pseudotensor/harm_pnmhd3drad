#include "step.h"
#include "radiation.h"
// MANUELTODO: choose EOS
#define ADIABATICEOS 0
#define BAROTROPICEOS 1
//#define STAREOS BAROTROPICEOS
#define STAREOS ADIABATICEOS

// 0: normal ZEUS
// 1: enthalpy evolution for barotropic EOS
#define EVOLVEENTHALPY 0

// set KAPPA
#define KAPPASTAR 0.35044


// add source terms due to grid expansion and contraction 

// MANUELTODO
#define OUTPUTFORCES 0


void stepvar_3d(void)
{
  // source steps
  int i,j,k;
  // kill the nowhere zones


  fracfloor_correct(); // adjust floor to give a set dynamic range in density

  transmagx1=transmagx2=0; // force since only applies in 2D

  if(subcyclen<=1){

#if(MDOTMEM)
    if(mdotin){
      injection_3d();
    }
#endif
    if(cool){
      cooling();
    }

#if(DIAGF)
   LOOPF{
     diagn[4][k][j][i]=s[1][k][j][i];
     diagn[5][k][j][i]=s[2][k][j][i];
     diagn[6][k][j][i]=v[1][1][k][j][i];
     diagn[7][k][j][i]=v[1][2][k][j][i];
     diagn[8][k][j][i]=v[1][3][k][j][i];
   }
#endif

#if(RAD)
      radstart();
      fld();
#endif

    // compute 1/c^2
    compute_invsol2fun();

    // substep 1
    if(ALFVENLIMIT==1){ // must do if mag==1 or not since step_pgc requires results of magprepare either way.  If no mag, just turn off ALFVENLIMIT.
      magprepare(); // need new rho (rho*) for gas pressure and magnetic pressure terms
      // rho doesn't change between here and transport
      // B doesn't change between here and moc
    }
    if(press){
      step_pgc() ; //includes curvature terms from div(rho*v*v)
    }

#if(DIAGF)
   LOOPF{
     diagn[9][k][j][i]=s[1][k][j][i];
     diagn[10][k][j][i]=s[2][k][j][i];
     diagn[11][k][j][i]=v[1][1][k][j][i];
     diagn[12][k][j][i]=v[1][2][k][j][i];
     diagn[13][k][j][i]=v[1][3][k][j][i];
   }
#endif

  }


  if(subcyclen>=1){

    if(visc_real==1){
#if(VISCMEM)
      nu_compute();
      step_visc_real() ;
      if(DOLOSSDIAG){
	viscous_flux();
      }
#endif
#if(DIAGF)
   LOOPF{
     diagn[14][k][j][i]=s[1][k][j][i];
     diagn[15][k][j][i]=s[2][k][j][i];
     diagn[16][k][j][i]=v[1][1][k][j][i];
     diagn[17][k][j][i]=v[1][2][k][j][i];
     diagn[18][k][j][i]=v[1][3][k][j][i];
   }
#endif
    }
  }


  if(subcyclen<=1){
    // substep 2
    if(visc_art==1){
      step_visc() ;
    }
#if(DIAGF)
   LOOPF{
     diagn[19][k][j][i]=s[1][k][j][i];
     diagn[20][k][j][i]=s[2][k][j][i];
     diagn[21][k][j][i]=v[1][1][k][j][i];
     diagn[22][k][j][i]=v[1][2][k][j][i];
     diagn[23][k][j][i]=v[1][3][k][j][i];
   }
#endif 

#if(RAD)
    gradv();
    moment();
#else
#if(STAREOS!=BAROTROPICEOS)    
    // otherwise don't need since done after density modified
    // assume t=0 data set correctly
    // substep 3
    if(ie){
      if(RELIE){
	step_relie();
      }
      else{
	if(wgam) step_ie() ;
      }
    }
#if(DIAGF)
   LOOPF{
     diagn[24][k][j][i]=s[1][k][j][i];
     diagn[25][k][j][i]=s[2][k][j][i];
     diagn[26][k][j][i]=v[1][1][k][j][i];
     diagn[27][k][j][i]=v[1][2][k][j][i];
     diagn[28][k][j][i]=v[1][3][k][j][i];
   }
#endif
 

#endif
#endif
    // CT and Alfven wave source terms
    if(mag==1){
      if(ALFVENLIMIT==0){ // need for normal rho
	magprepare();
      }
      if(MOC3DVER==1){
	// needs rholimited if alfvenlimit==1, otherwise need osqrtrho
	lorentz_3d(); // compute lorentz force first as suggested by Hawley(see ZEUS-3D)
      }

      if(RESMEM&&(res_real==1)){
	current_compute(123); // need full current for abs(J) in nu_res_compute(rreal==2) || step_res
	nu_res_compute();
	step_res() ; // needs all 3 components
	// need to redo magprepare if alfvenlimit==1 and COMPDIM==2 and MOC2DVER==0
	if((ALFVENLIMIT==1)&&(COMPDIM<=2)&&(MOC2DVER==0)){
	  magprepare(); // since bz changed
	}
      }
      if(stepmocct){
	// above step_res didn't modify field, only ie, so no need to compute current again
	if(MOC3DVER==0) moc_ct_3d_v1(); // older non-split version(includes lorentz force, using old velocities for field)
	if(MOC3DVER==1) moc_ct_3d_v2() ; // uses all 3 components(no lorentz force here, using ALL new velocities for field)
	// rholimited or sqrtrho no longer needed until timestep() (but rho changes in transport, so wait till before timestep to do magprepare() )
      }
      if(DOLOSSDIAG){
	magnetic_flux();
      }
    }
  }
  /*
  // GODMARK -- commented
  fprintf(log_file,"%ld: before step_trans_3d symmetry check\n",nstep); fflush(stdout);
  symmetry_check(0); // before bound symmetry check// GODMARK -- commented
  */  
  if(subcyclen<=1){    
    // transport steps
    if(trans==1){
      step_trans_3d() ;
    }
#if(DIAGF)
   LOOPF{
     diagn[29][k][j][i]=s[1][k][j][i];
     diagn[30][k][j][i]=s[2][k][j][i];
     diagn[31][k][j][i]=v[1][1][k][j][i];
     diagn[32][k][j][i]=v[1][2][k][j][i];
     diagn[33][k][j][i]=v[1][3][k][j][i];
   }
#endif
 

#if(LOOPTYPE!=1)
    if((COMPUTELOSSDIAG==1)&&(DOLOSSDIAG)){
      hydro_flux();
    }
#endif
  }

#if(STAREOS==BAROTROPICEOS)
  // fix internal energy to be consistent with density
  if(ie){
    if(RELIE){
      step_relie();
      }
    else{
      if(wgam) step_ie() ;
    }
  }
#endif


  /*
  // GODMARK -- commented
  fprintf(log_file,"%ld: after step_trans_3d symmetry check\n",nstep); fflush(stdout);
  symmetry_check(-1); // before bound symmetry check// GODMARK -- commented
  image(888,-1,-1,0,0);
  myexit(0);
  */

  t += dt ;

  //  if(OUTPUTFORCES) exit(0);
}



//compute things that are functions of time besides the analytic solution

#define OTAU 100.0 // 1 over tau, the time that alpha becomes more like alpha_real0

// 1d/2d/3d valid
void tdep_compute_normal(void)
{
  int i,j,k;
  static int firsttime=1;
  static int switchbc=0;
  int error=0;

  SFTYPE ftemp;


  //////// MODIFY floor
  //
  /*
    if(t>.7){
    DTd=DTl=DTfloor=.001;
    }

    if(t>.8){
    DTd=DTl=DTfloor=.0005;
    }
  */
  //
  ///////////////////////


  //////MODIFY alpha_real
  //
  //#if(COORD==3)
  //alpha_real=alpha_real0*(1.-exp(-t*OTAU))+1.E-7;
  //if(t>=10000.0) alpha_real=0.0;
  //  if(t>=100.0) alpha_real=alpha_real0;
  //if(t>=20000.0) alpha_real=alpha_real0;
  //else alpha_real=0.0;
  //  if(t>1.0) alpha_real=alpha_real0;
  //  else alpha_real=0.0;
  //#elif(COORD==1)
  if(t>=0.0) {
    alpha_real=alpha_real0;


    sigactmpi(0);
#if(USEMPI)
    for(irowmpi=0;irowmpi<ncpux2-1;irowmpi++){
      bound1D_mpi(Sigar, Sigar, 1, irowmpi, 0);
      sigactmpi(1);
    }
#endif
    bound(alpha_reala,NULL,-10-1,0,0);

    qparampi(0);
#if(USEMPI)
    for(irowmpi=ncpux2-1;irowmpi>0;irowmpi--){
      bound1D_mpi(alpha_reala, alpha_reala , 1, irowmpi, 1);
      qparampi(1);
    }
#endif

   LOOPC{
     if(optd[k][j][i]<0.003){
//        alpha_reala[k][j][i]=alpha_real0*pow(1./fabs(log10(optd[k][j][i])),0.2/x[2][1][i]);
	 alpha_reala[k][j][i]=1.e-20;
     }
   } 

/*
#if(THINCOOL)
   LOOPC{
     if(coolerKai[k][j][i]==SSMALL*0.1){
        alpha_reala[k][j][i]=1.e-10;
     }
   }
#endif
*/
    bound(alpha_reala,NULL,-10-1,0,0);
  }
  //#endif

  if(runtype>0){
    if(t>=0.0) alpha_real=alpha_real0; // go ahead and turn on when using reentrant data as initial data
    sigactmpi(0);
#if(USEMPI)
    for(irowmpi=0;irowmpi<ncpux2-1;irowmpi++){
      bound1D_mpi(Sigar, Sigar, 1, irowmpi, 0);
      sigactmpi(1);
    }
#endif
    bound(alpha_reala,NULL,-10-1,0,0);
    qparampi(0);
#if(USEMPI)
    for(irowmpi=ncpux2-1;irowmpi>0;irowmpi--){
      bound1D_mpi(alpha_reala, alpha_reala , 1, irowmpi, 1);
      qparampi(1);
    }
#endif

   LOOPC{
    if(optd[k][j][i]<0.003){
//        alpha_reala[k][j][i]=alpha_real0*pow(1./fabs(log10(optd[k][j][i])),0.2/x[2][1][i]);
	alpha_reala[k][j][i]=1.e-20;
     }
   }
 
/*

#if(THINCOOL)
   LOOPC{
     if(coolerKai[k][j][i]==SSMALL*0.1){
        alpha_reala[k][j][i]=1.e-10;
     }
   }
#endif
*/
   bound(alpha_reala,NULL,-10-1,0,0);
  }


  // res
  if(t>=0.0) resist_real=resist_real0;
  

  // cooling
  //  if(t<=50) cool=0; else cool=1;
  //if(t>=30000.0) cool=1;
  //else cool=0;
  //  fprintf(stderr,"\n%15.10g\n",alpha_real);
  //
  /////////////////////////////




  /////MODIFY BOUNDARY CONDITIONS! (assumes grid doesn't change, so don't need to run through anything but the below 2 functions)
  //
  // start with initial bc

  //
  // just switch inner boundary to reflective
  /*
    if((t>=4000.0)&&(switchbc==0)){

    error+=init_mainbc(periodicx1,skipix1,!reflectix1,reflectox1,periodicx2,skipix2,reflectix2,reflectix2,periodicx3,skipix3,reflectix3,reflectox3);

    error+=init_bc(1,1,4,1,1,1,1);

    switchbc++;
    fprintf(log_file,"swichbc: %d\n",switchbc);
    }
    // switch back
    else if((t>=6000.0)&&(switchbc==1)){

    error+=init_mainbc(periodicx1,skipix1,!reflectix1,reflectox1,periodicx2,skipix2,reflectix2,reflectix2,periodicx3,skipix3,reflectix3,reflectox3);

    error+=init_bc(1,4,4,1,1,1,1);

    switchbc++;
    fprintf(log_file,"swichbc: %d\n",switchbc);
    }
    // switch outer boundary to reflective now
    else if((t>=8000.0)&&(switchbc==2)){

    error+=init_mainbc(periodicx1,skipix1,reflectix1,!reflectox1,periodicx2,skipix2,reflectix2,reflectix2,periodicx3,skipix3,reflectix3,reflectox3);

    error+=init_bc(1,4,1,1,1,1,1);

    switchbc++;
    fprintf(log_file,"swichbc: %d\n",switchbc);
    }
    // switch back
    else if((t>=10000.0)&&(switchbc==3)){

    error+=init_mainbc(periodicx1,skipix1,reflectix1,!reflectox1,periodicx2,skipix2,reflectix2,reflectix2,periodicx3,skipix3,reflectix3,reflectox3);

    error+=init_bc(1,4,4,1,1,1,1);

    switchbc++;
    fprintf(log_file,"swichbc: %d\n",switchbc);
    }
  */
  //
  ///////////////////////////////


  /////////////////////////////
  //
  // Change kinetic energy computation
  /*
    if((t>=20000)&&(t<23000.0)){
    kever=0;
    }
    else if(t>=23000.0){
    kever=1;
    }
  */  



  //
  firsttime=0;
}


void injection_3d(void)
{

}
void cooling_nocooling(void)
{

}

void cooling_simple1(void)
{
  int i,j,k,l;
  FTYPE HORmod;
  static int firsttime=1;
  FTYPE zcoord,Radius,omega,cs2,H2,H20,omega02;
  FTYPE volume;
  FTYPE eninit;
  FTYPE coolerden,coolerie;
  static FTYPE max=-1E30,min=1E30;
  static FTYPE maxallow=0, minallow=0;
  static FTYPE maxtime,starttime;
  // compute

  HORmod=1.0; // 1.0 would keep H/R constant with cooling

  if(firsttime==1){
    maxtime=100.0;
    starttime=t;
    if(starttime>maxtime){
      starttime=maxtime; // so initial data is used as reference
      maxallow=minallow=0;
    }
  }

  LOOPC{
    Radius=x[2][1][i]*sin(x[2][2][j]);
    zcoord=x[2][1][i]*cos(x[2][2][j]);
    omega02=1.0/(Radius*(Radius-2.0)*(Radius-2.0));
    omega=v[1][3][k][j][i]/Radius;
    cs2=gam*(gam-1.0)*s[2][k][j][i]/s[1][k][j][i]; // cs^2
    H2=cs2/(omega*omega);
    H20=(HOR*HORmod)*(HOR*HORmod)*Radius*Radius; // H/R=constant=HOR (from analsol.c tori1sol)
    //		coolerden=coolfact*omega*s[2][k][j][i]*(H2-H20)/H20*dt; // energy density lost
    coolerden=coolfact*omega*s[2][k][j][i]*(cs2/(H20*omega02)-1)*dt; // energy density lost

    coolerden*=exp(-zcoord*zcoord/(2.0*H20));

    // limit heating or cooling
    if(t<=maxtime){
      if(coolerden>maxallow) maxallow=coolerden;
      if(coolerden<minallow) minallow=coolerden;
    }
    else{
      //if(coolerden>maxallow) coolerden=maxallow; // limits cooling
      if(coolerden<minallow) coolerden=minallow; // limits heating
    }

    // cool it!
    eninit=s[2][k][j][i]; // internal energy density before loss
    s[2][k][j][i]-=coolerden; // internal energy density lost
		
    volume=DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k);

    // floor correction here means don't radiate as much...As opposed to still saying I radiate but fill back in with floor
    if(FORCEIEINTERNAL){
      if(s[2][k][j][i]<IEFLOOR){	
	if(DOFLOORD2==1){
	  fprintf(logfl_file,"corrected en in cooling: t: %15.10g %d %d %d bad: %15.10g\n",t,k,j,i,s[2][k][j][i]);
	}
	s[2][k][j][i]=IEFLOOR;
	coolerie=(eninit-IEFLOOR)*volume;
      }
      else{
	// normal case
	coolerie=coolerden*volume; // full is ok then
      }
    }
    else{// have to check if no floor, or might get floating point exceptions!      
      if(s[2][k][j][i]<=0.0){
	fprintf(fail_file,"catastrophic cooling, new ie density at %d %d %d is %15.10g due to a change of %15.10g\n",k,j,i,s[2][k][j][i],-coolerden);
	myexit(1);
      }
      else{
	// normal case
	coolerie=coolerden*volume; // full is ok then
      }
    }
		
    if(accountstore[k][j][i]) radiations[2]+=coolerie; // internal energy lost
    if(coolerie>max){
      max=coolerie;
      fprintf(log_file,"t=%15.10g: new maximum cool amount: %15.10g @ %d %d %d\n",t,max,k,j,i);
      fprintf(log_file,"%15.10g %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g\n",s[1][k][j][i],s[2][k][j][i],s[3][k][j][i],v[1][1][k][j][i],v[1][2][k][j][i],v[1][3][k][j][i],v[2][1][k][j][i],v[2][2][k][j][i],v[2][3][k][j][i]);
      fflush(log_file);
    }
    if(coolerie<min){
      min=coolerie;
      fprintf(log_file,"t=%15.10g: new minimum cool amount: %15.10g @ %d %d %d\n",t,min,k,j,i);
      fprintf(log_file,"%15.10g %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g\n",s[1][k][j][i],s[2][k][j][i],s[3][k][j][i],v[1][1][k][j][i],v[1][2][k][j][i],v[1][3][k][j][i],v[2][1][k][j][i],v[2][2][k][j][i],v[2][3][k][j][i]);
      fflush(log_file);
    }
  }



  bound(NULL,NULL,2,0,0); // bound scalar just changed
  

  firsttime=0;


}


void cooling_twolayer(void)
{
  int i,j,k,l;
  static int firsttime=1;
  FTYPE zcoord,Radius,omega,cs2,H2,omega02;
  FTYPE volume;
  FTYPE eninit;
  FTYPE coolerden,coolerie;
  static FTYPE max=-1E30,min=1E30;
  static FTYPE maxallow=0, minallow=0;
  static FTYPE maxtime,starttime;
  // compute

  LOOPC{
    Radius=x[2][1][i]*sin(x[2][2][j]);
    zcoord=x[2][1][i]*cos(x[2][2][j]);
    omega02=1.0/(Radius*(Radius-rgp)*(Radius-rgp));
    omega=v[1][3][k][j][i]/Radius;
    cs2=gam*(gam-1.0)*s[2][k][j][i]/s[1][k][j][i]; // cs^2
    H2=cs2/(omega*omega);
    //		coolerden=coolfact*omega*s[2][k][j][i]*(H2-H20)/H20*dt; // energy density lost
    coolerden=coolfact*omega*s[2][k][j][i]*((cs2-cs02[k][j][i])/(cs2+cs02[k][j][i]))*alpha_reala[k][j][i]*dt; // energy density lost 

     // cool it!
    eninit=s[2][k][j][i]; // internal energy density before loss
    s[2][k][j][i]-=coolerden; // internal energy density lost
		
    volume=DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k);

    // floor correction here means don't radiate as much...As opposed to still saying I radiate but fill back in with floor
    if(FORCEIEINTERNAL){
      if(s[2][k][j][i]<IEFLOOR){	
	if(DOFLOORD2==1){
	  fprintf(logfl_file,"corrected en in cooling: t: %15.10g %d %d %d bad: %15.10g\n",t,k,j,i,s[2][k][j][i]);
	}
	s[2][k][j][i]=IEFLOOR;
	coolerie=(eninit-IEFLOOR)*volume;
      }
      else{
	// normal case
	coolerie=coolerden*volume; // full is ok then
      }
    }
    else{// have to check if no floor, or might get floating point exceptions!      
      if(s[2][k][j][i]<=0.0){
	fprintf(fail_file,"catastrophic cooling, new ie density at %d %d %d is %15.10g due to a change of %15.10g\n",k,j,i,s[2][k][j][i],-coolerden);
	myexit(1);
      }
      else{
	// normal case
	coolerie=coolerden*volume; // full is ok then
      }
    }
		
    if(accountstore[k][j][i]) radiations[2]+=coolerie; // internal energy lost
    if(coolerie>max){
      max=coolerie;
      fprintf(log_file,"t=%15.10g: new maximum cool amount: %15.10g @ %d %d %d\n",t,max,k,j,i);
      fprintf(log_file,"%15.10g %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g\n",s[1][k][j][i],s[2][k][j][i],s[3][k][j][i],v[1][1][k][j][i],v[1][2][k][j][i],v[1][3][k][j][i],v[2][1][k][j][i],v[2][2][k][j][i],v[2][3][k][j][i]);
      fflush(log_file);
    }
    if(coolerie<min){
      min=coolerie;
      fprintf(log_file,"t=%15.10g: new minimum cool amount: %15.10g @ %d %d %d\n",t,min,k,j,i);
      fprintf(log_file,"%15.10g %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g\n",s[1][k][j][i],s[2][k][j][i],s[3][k][j][i],v[1][1][k][j][i],v[1][2][k][j][i],v[1][3][k][j][i],v[2][1][k][j][i],v[2][2][k][j][i],v[2][3][k][j][i]);
      fflush(log_file);
    }
  }



  bound(NULL,NULL,2,0,0); // bound scalar just changed
  

  firsttime=0;


}


void optdepth1(int initmode)
{
  int i,j,k;

  LOOPC3{
     LOOPC1{
      if(initmode==0){
       optd[k][-1][i]=0.;
       tatm4[k][-1][i]=0.;
      }
       LOOPC2{
         optd[k][j][i]=optd[k][j-1][i]+dx[2][2][j]*G2(2,i)*0.5*(kapk[k][j][i]+sigk[k][j][i]+kapk[k][j-1][i]+sigk[k][j-1][i]);
	 tatm4[k][j][i]=tatm4[k][j-1][i]+4.*kapk[k][j][i]*tkai[k][j][i]*tkai[k][j][i]*tkai[k][j][i]*tkai[k][j][i]*dx[1][2][j]*G2(2,i);
         if(optd[k][j][i]<optthick){
           thingrid[k][j][i]=1;
         }else{
           thingrid[k][j][i]=0;
         }
       }
     }
  }
//  fprintf(log_file,"optd %15.10g %15.10g %15.10g %15.10g \n",optd[0][-N2BND+1][10],optd[0][0][10],optd[0][N2-N2BND+1][10],optd[0][N2][10]);
//  fprintf(log_file,"tatm4 %15.10g %15.10g %15.10g %15.10g \n",tatm4[0][-N2BND+1][10],tatm4[0][0][10],tatm4[0][N2-N2BND+1][10],tatm4[0][N2][10]);
}

void optdepth2(int initmode)
{
   int i,j,k;
   FTYPE tbdry,taobdry;
   LOOPC3{
     LOOPC1{
       if(initmode==0){
         taoup[k][N2][i]=0.;
	 fbdry[k][N2][i]=0.;
//	 diagn[3][k][N2][i]=0.;
       }
       for(j=N2-1;j>=-1;j--){
         if(optd[k][j][i]<optthick){
	   if(optd[k][j+1][i]>optthick){
             taoup[k][j][i]=(kapk[k][j][i]+sigk[k][j][i])*dx[1][2][j]*G2(2,i)/2.;
             tbdry=tkai[k][j+1][i];
             taobdry=optd[k][j+1][i];
             fbdry[k][j][i]=4.*sigmacontnew*(tbdry*tbdry*tbdry*tbdry-tenv[i]*tenv[i]*tenv[i]*tenv[i])/3./(taobdry+2./3.); // here assume optical depth at different wavelength are the same,gray case, otherwise after tenv[i] there should be exp(-tao)
	     fbdry[k][j+1][i]=fbdry[k][j][i];
//	     diagn[3][k][j][i]=fbdry[k][j][i];
           }else{
             taoup[k][j][i]=taoup[k][j+1][i]+0.5*(kapk[k][j][i]+sigk[k][j][i]+kapk[k][j+1][i]+sigk[k][j+1][i])*dx[2][2][j+1]*G2(2,i);
	     fbdry[k][j][i]=fbdry[k][j+1][i];
//	     diagn[3][k][j][i]=diagn[3][k][j+1][i];
           }
         }
       }
     }
   }
//   fprintf(log_file,"diag3 %g",diag[0][0][0]);
//   fprintf(log_file,"taoup %15.10g %15.10g %15.10g %15.10g \n",taoup[0][-N2BND+1][10],taoup[0][0][10],taoup[0][N2-N2BND+1][10],taoup[0][N2][10]);
//  fprintf(log_file,"fbdry %15.10g %15.10g %15.10g %15.10g \n",fbdry[0][-N2BND+1][10],fbdry[0][0][10],fbdry[0][N2-N2BND+1][10],fbdry[0][N2][10]);

}
/*     LOOPC1{
        if(fbdry[i]==0.){
          fbdry[i]=fbdry[i-1];
        }
     }
     for(i=N1-1;i>=0;i--){
        if(fbdry[i]==0.){
          fbdry[i]=fbdry[i+1];
        }
     }
*/
   


void coolingKaistep(int m)
{
#if(KAICOOL)
  int i,j,k,l,istart,iend,firstthin,mini,mfac,mfac2;
  FTYPE der1,der2,dtdx,yp,coolerden,dttemp,cs2,t4,dtsmallSP,qa,qb,qc,qd,qe;
  static FTYPE drta[2][N3M][N2M][N1M],tatm4a[N2M][N1M],betaa[N3M][N2M][N1M];
  static FTYPE (*tt),(*dtde),(*dkapdt),(*tdr)[N2M][N1M],(*drt)[N3M][N2M][N1M],(*beta)[N2M][N1M];
  drt = (FTYPE (*)[N3M][N2M][N1M])(&(drta[-1][N3BND][N2BND][N1BND]));
  beta = (FTYPE (*)[N2M][N1M])(&(betaa[N3BND][N2BND][N1BND]));
  tdr=work5;
  dtde=work5x;
  dkapdt=work6x;
// y, beta and dr term  
  istart=N1BND;
  iend=N1BND+N1+1;
  optthick=10.;
  mini=5.;
  mfac=10.;
  mfac2=10.;

  if(m<=mini||m%mfac==0){
    LOOPC1{
      tenv[i]=1.;
    }
    LOOPF3 {
      LOOPF2 {
        temp_(&s[2][k][j][-N1BND],&s[1][k][j][-N1BND],&gam,&istart,&iend,&tkai[k][j][-N1BND],&dtde[-N1BND]);
        absorp_(&tkai[k][j][-N1BND],&s[1][k][j][-N1BND],&istart,&iend,&kapk[k][j][-N1BND],&dkapdt[-N1BND]);
        scatt_(&tkai[k][j][-N1BND],&s[1][k][j][-N1BND],&istart,&iend,&sigk[k][j][-N1BND]);
      }
    }
  }else{
    LOOPF3 {
      LOOPF2 {
        temp_(&s[2][k][j][-N1BND],&s[1][k][j][-N1BND],&gam,&istart,&iend,&tkai[k][j][-N1BND],&dtde[-N1BND]);
      }
    }
  }

#if(LIMITER)
  if(m<=mini||m%mfac2==0){
    LOOPRC3 {
      LOOPRC2 {
        LOOPRC1 {
         der1=(tkai[k][j][i+1]-tkai[k][j][i-1])/(dx[2][1][i]+dx[2][1][i+1]);
         der2=(tkai[k][j+1][i]-tkai[k][j-1][i])/((dx[2][2][j]+dx[2][2][j+1])*G2(2,i));
         dtdx=sqrt(der1*der1+der2*der2);
         yp=4./(kapk[k][j][i]+sigk[k][j][i])/tkai[k][j][i]*dtdx;
         beta[k][j][i]=(2.+yp)/(6.+3.*yp+yp*yp);
        }
      }
    }
  }
#else
   LOOPRC3 {
      LOOPRC2 {
        LOOPRC1 {
	  beta[k][j][i]=1./3.;
	}
      }
   }
#endif

   LOOPRC3 {
      LOOPRC2 {
        LOOPRC1 {
         tdr[k][j][i]=-16.*sigmacontnew*beta[k][j][i]*tkai[k][j][i]*tkai[k][j][i]*tkai[k][j][i]/(kapk[k][j][i]+sigk[k][j][i]);
        }
      }
    }


    LOOPC3{
      LOOPC2{
        LOOPFLD1C{
      drt[1][k][j][i] = 0.5*(tdr[k][j][i-1] + tdr[k][j][i]);
        }
      }
    }
    LOOPC3{
      LOOPC1{
        LOOPFLD2C{
      drt[2][k][j][i] = 0.5*(tdr[k][j-1][i] + tdr[k][j][i]);
        }
      }
    }

  if(m<=mini||m%mfac==0){
  optdepth1(0);
#if(USEMPI)
  for(irowmpi=0;irowmpi<ncpux2-1;irowmpi++){
      bound1D_mpi(optd, tatm4, 2, irowmpi, 0);
      optdepth1(1);
  }
  bound_mpi(optd,NULL,-10-1,0,0);
#endif

  optdepth2(0);
#if(USEMPI)
  for(irowmpi=ncpux2-1;irowmpi>0;irowmpi--){
      bound1D_mpi(taoup, fbdry, 2, irowmpi, 1);
      optdepth2(1);
  }
#endif
  } 
/*     LOOPC1{
	if(fbdry[i]==0.){
	  fbdry[i]=fbdry[i-1];
	}
     }
     for(i=N1-1;i>=0;i--){
	if(fbdry[i]==0.){
	  fbdry[i]=fbdry[i+1];
	}
     }
*/

// DIV(F) term
  dtsmallSP=1.;
  LOOPC3 {
  LOOPC1 {
  LOOPC2 {
   if(thingrid[k][j][i]==0){
          qa = G2(1,i+1)*G3(1,i+1)/dx[2][1][i+1]*drt[1][k][j][i+1];
          qb = G2(1,i  )*G3(1,i  )/dx[2][1][i]*drt[1][k][j][i];
          qc = G2(2,i)*G2(2,i) ;
          qd = G4(1,j+1)/dx[2][2][j+1]*drt[2][k][j+1][i];
          qe = G4(1,j  )/dx[2][2][j  ]*drt[2][k][j][i] ;
	  if(x[2][1][i-1]<Rinner||thingrid[k][j][i-1]==1){qb=0.;}
          if(x[2][1][i+1]>Router||thingrid[k][j][i+1]==1){qa=0.;}
	  if(optd[k][j][i]>optthick&&optd[k][j-1][i]<optthick){
	    coolerden =(qa*(tkai[k][j][i+1]-tkai[k][j][i])-qb*(tkai[k][j][i]-tkai[k][j][i-1]))/dvl[1][1][i]+
     (G4(1,j+1)/dx[2][2][j+1]/G2(2,i)*drt[2][k][j+1][i]*(tkai[k][j+1][i]-tkai[k][j][i])+G4(1,j)*fbdry[k][j][i])/G2(2,i)/dvl[1][2][j];
	  diagn[4][k][j][i]=G4(1,j)*fbdry[k][j][i]/G2(2,i)/dvl[1][2][j];
	  diagn[5][k][j][i]=qb*(tkai[k][j][i]-tkai[k][j][i-1])/dvl[1][1][i];
	  diagn[6][k][j][i]=((qa*(tkai[k][j][i+1]-tkai[k][j][i])-qb*(tkai[k][j][i]-tkai[k][j][i-1]))/dvl[1][1][i])/
     ((G4(1,j+1)/dx[2][2][j+1]/G2(2,i)*drt[2][k][j+1][i]*(tkai[k][j+1][i]-tkai[k][j][i])+G4(1,j)*fbdry[k][j][i])/G2(2,i)/dvl[1][2][j]);
	  }else{
            coolerden =
       (qa*(tkai[k][j][i+1]-tkai[k][j][i])-qb*(tkai[k][j][i]-tkai[k][j][i-1]))/dvl[1][1][i]+
       (qd*(tkai[k][j+1][i]-tkai[k][j][i])-qe*(tkai[k][j][i]-tkai[k][j-1][i]))/qc/dvl[1][2][j];
	   diagn[4][k][j][i]=qe*(tkai[k][j][i]-tkai[k][j-1][i])/qc/dvl[1][2][j];
	   diagn[5][k][j][i]=qb*(tkai[k][j][i]-tkai[k][j][i-1])/dvl[1][1][i];
	   diagn[6][k][j][i]=(qa*(tkai[k][j][i+1]-tkai[k][j][i])-qb*(tkai[k][j][i]-tkai[k][j][i-1]))/dvl[1][1][i]/
       ((qd*(tkai[k][j+1][i]-tkai[k][j][i])-qe*(tkai[k][j][i]-tkai[k][j-1][i]))/qc/dvl[1][2][j]);
	  }
   }else{
    	  if(fbdry[k][j][i]<0.){fbdry[k][j][i]=0.;}
#if(DURISEN)
 	  coolerden = 4.*kapk[k][j][i]*sigmacontnew*tkai[k][j][i]*tkai[k][j][i]*tkai[k][j][i]*tkai[k][j][i]
   -4.*kapk[k][j][i]*sigmacontnew*(3./4.*fbdry[k][j][i]/sigmacontnew*(optd[k][j][i]+2./3.)+tenv[i]*tenv[i]*tenv[i]*tenv[i]);
#else
  	  coolerden = 4.*kapk[k][j][i]*sigmacontnew*tkai[k][j][i]*tkai[k][j][i]*tkai[k][j][i]*tkai[k][j][i]
   -kapk[k][j][i]*sigmacontnew*tenv[i]*tenv[i]*tenv[i]*tenv[i]*exp(-optd[k][j][i])
   -(kapk[k][j][i]+sigk[k][j][i])*2.*fbdry[k][j][i]-kapk[k][j][i]*tatm4[k][j][i]*sigmacontnew;
#endif
#if(THINCOOL)
	  if(optd[k][j][i]<1.e-3){
            coolerden=SSMALL*0.1;
	   if(m<=mini||m%mfac==0){
            t4=fbdry[k][j][i]/4./sigmacontnew*3.*(optd[k][j][i]+2./3.)+tenv[i]*tenv[i]*tenv[i]*tenv[i];
            cs2=gam*8.314472e7/eneunit/2.4*sqrt(sqrt(t4));
            sthin[k][j][i]=s[1][k][j][i]*cs2/gam/(gam-1);
            sthin[k][j][i]=s[1][k][j][i]/164782.*sqrt(sqrt(t4));
	   }
	  }
#endif
   }

   diagn[1][k][j][i]=alpha_reala[k][j][i];
   diagn[2][k][j][i]=optd[k][j][i];
/*
   if(thingrid[k][j][i]==0){
      diagn[3][k][j][i]=qd*(tkai[k][j+1][i]-tkai[k][j][i])/G2(2,i)/G4(1,j+1);
   }else{
      diagn[3][k][j][i]=-fbdry[k][j][i];
   }
*/
   dttemp=0.1*s[2][k][j][i]/fabs(coolerden);
   diagn[3][k][j][i]=1./dttemp;
   coolerKai[k][j][i]=coolerden;
   if(dttemp<dtsmallSP){
	dtsmallSP=dttemp;
	jsmall=j;
	ismall=i;
   }
  }
  }
  }
#if(SUBCOOL)
  if(numprocs>1){
#if(USEMPI)
    MPI_Allreduce(&dtsmallSP, &dtsmall, 1, MPI_FTYPE, MPI_MIN, MPI_COMM_WORLD);
#endif
  }else{
    dtsmall=dtsmallSP;
  }
//  if(t<10.){
//	dtsmall=dtsmall/10.;
//  }
#endif
#endif
}

void coolingthin(void)
{
#if(KAICOOL)
  int k,j,i;
#if(THINCOOL)
   LOOPC{
     if(coolerKai[k][j][i]==SSMALL*0.1){
        s[2][k][j][i]=sthin[k][j][i];
	if(FORCEIEINTERNAL){
      	  if(s[2][k][j][i]<IEFLOOR){
            if(DOFLOORD2==1){
              fprintf(logfl_file,"corrected en in cooling: t: %15.10g %d %d %d bad: %15.10g\n",t,k,j,i,s[2][k][j][i]);
            }
            s[2][k][j][i]=IEFLOOR;
          }
        }
     }
   }
#endif


#endif
}


void coolingKai(void)
{
#if(KAICOOL)
  int k,j,i,m,iend,msub;
  FTYPE eninit,coolerie,volume,dtsub;
  static FTYPE max=-1E30,min=1E30;

dtsub=0.;
iend=0;
LOOPC{
 tot[0][k][j][i]=s[2][k][j][i];
// diagn[2][k][j][i]=s[2][k][j][i];
}
#if(SUBCOOL)
msub=1000000;
for(m=1;m<=msub;m++){

 if(m!=1){
	coolingKaistep(m);
 }
 dtsub=dtsub+dtsmall;
 if(dtsub>=dt){
  dtsmall=dt-(dtsub-dtsmall);
  iend=1;
 }
#else
  dtsmall=dt;
#endif
    
  LOOPC{
   eninit=s[2][k][j][i]; // internal energy density before loss
#if(THINCOOL)
   if(coolerKai[k][j][i]!=SSMALL*0.1){
        s[2][k][j][i]-=coolerKai[k][j][i]*dtsmall;
   }else{
	coolerKai[k][j][i]=(s[2][k][j][i]-sthin[k][j][i])/dtsmall;
	s[2][k][j][i]=sthin[k][j][i]; 	
   }
#else
   s[2][k][j][i]-=coolerKai[k][j][i]*dtsmall;
#endif
   volume=DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k);
    // floor correction here means don't radiate as much...As opposed to still saying I radiate but fill back in with floor
    if(FORCEIEINTERNAL){
      if(s[2][k][j][i]<IEFLOOR){
        if(DOFLOORD2==1){
          fprintf(logfl_file,"corrected en in cooling: t: %15.10g %d %d %d bad: %15.10g\n",t,k,j,i,s[2][k][j][i]);
        }
        s[2][k][j][i]=IEFLOOR;
        coolerie=(eninit-IEFLOOR)*volume*dtsmall;
      }
      else{
        // normal case
        coolerie=coolerKai[k][j][i]*volume*dtsmall; // full is ok then
      }
    }
    else{// have to check if no floor, or might get floating point exceptions!   
      if(s[2][k][j][i]<=0.0){
        fprintf(fail_file,"catastrophic cooling, new ie density at %d %d %d is %15.10g due to a change of %15.10g\n",k,j,i,s[2][k][j][i],-coolerKai[k][j][i]);
        myexit(1);
      }
      else{
        // normal case
        coolerie=coolerKai[k][j][i]*volume*dtsmall; // full is ok then
      }
    }

    if(accountstore[k][j][i]) radiations[2]+=coolerie; // internal energy lost
    if(coolerie>max){
      max=coolerie;
//      fprintf(log_file,"t=%15.10g: new maximum cool amount: %15.10g @ %d %d %d\n",t,max,k,j,i);
//      fprintf(log_file,"%15.10g %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g\n",s[1][k][j][i],s[2][k][j][i],s[3][k][j][i],v[1][1][k][j][i],v[1][2][k][j][i],v[1][3][k][j][i],v[2][1][k][j][i],v[2][2][k][j][i],v[2][3][k][j][i]);
      fflush(log_file);
    }
    if(coolerie<min){
      min=coolerie;
//      fprintf(log_file,"t=%15.10g: new minimum cool amount: %15.10g @ %d %d %d\n",t,min,k,j,i);
//      fprintf(log_file,"%15.10g %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g  %15.10g\n",s[1][k][j][i],s[2][k][j][i],s[3][k][j][i],v[1][1][k][j][i],v[1][2][k][j][i],v[1][3][k][j][i],v[2][1][k][j][i],v[2][2][k][j][i],v[2][3][k][j][i]);
      fflush(log_file);
    }
  }


  bound(NULL,NULL,2,0,0); // bound scalar just changed
#if(SUBCOOL)
  if(iend==1)goto last;
//  fprintf(stderr,"dt %g dtsmall %g\n",dt,dtsmall);
//   fprintf(log_file,"cool %d %g",m,dtsub);
 }
   fprintf(log_file,"failed cooling dt %d",m);
   myexit(1);
   last:
	substepm=m;
#endif
//   fprintf(log_file,"cooling %d %g\n",m, dt);
//  fprintf(log_file,"t dt %g %g \n",t,dt);
   LOOPC{
    tot[1][k][j][i]=s[2][k][j][i];
//    diagn[1][k][j][i]=coolerKai[k][j][i];
//    diagn[3][k][j][i]=s[2][k][j][i];
   }
#endif
}


// cooling routine (compute every timestep??)
// this fixes 

// h*nu in log10 (in eV)
#define NUI (1.0)
#define NUO (6.0)

// initial attempt of theta range
#define THETAS (0.0)
#define THETAF (300.0)

#define COOLBOTTOM (1.0) // in theta

// linear interp
#define FUNCTIONCOOL(i,ratio) (funcool[i]*(1.-ratio)+funcool[ip1]*ratio)


// 1d/2d/3d valid
void cooling_brem(void)
{
  int i,j,k,l;
  static int firsttime=1;
  //  static FTYPE dlognu;
  //static FTYPE JCONST,JCONSTINT;
  //static FTYPE SQRT3;
  //FTYPE gaunt;
  //FTYPE temperature,ratioe,freq;
  FTYPE volume;
  FTYPE ftemp,ftemp1,ftemp2;
  //static FTYPE IECONVERT,IECONVERTDEN;
  static FTYPE theta,mpome;
  static FTYPE dtheta,odtheta,thetas,thetaf;
  int cooli;
  FTYPE ratio;
  FTYPE eninit;
  FTYPE coolerden,coolerie;
  // compute 

  /* cgs units
     consts = {Z -> 1.06718, e -> 4.803*10^(-10), me -> 9.110*10^(-28), 
     mp -> 1.673*10^(-24), k -> 1.381*10^(-16), G -> 6.670*10^(-8), 
     c -> 2.99792458*10^(10), h -> 6.626*10^(-27), Mstar -> 1.989*10^33, 
     MBH -> 10^8}
  */

#define ZNUM (1.06718)
#define ECHARGE (4.803E-10)
#define MELECTRON (9.110E-28)
#define MPROTON (1.673E-24)
#define KBOLTZ (1.381E-16)
#define GRAVCONST (6.670E-8)
#define SPEEDLIGHT (2.99792458E10)
#define PLANK (6.626E-27)
#define MASSSOLAR (1.989E33)
#define MH1 (1.000545*MPROTON) //Mass of hydrogen atom
#define GAMMA 1.781

  if(firsttime==1){

    // no 2*Pi if only working with axisym assumption(fluxes/quantities are really 2*Pi bigger)
    //    JCONST=16./3.*sqrt(M_PI/6)*pow(ECHARGE,6.)*pow(ZNUM,2.0)/pow(MELECTRON*SPEEDLIGHT*SPEEDLIGHT,1.5)/pow(MPROTON*0.5,2.);
    // 16/3*(Pi/6)^(1/2)*(e^3*Z)^2 / (me*c^2)^(3/2) * 1/(mi^2*h)
    //    JCONSTINT=JCONST/PLANK;

    //dlognu=(NUO-NUI)/((FTYPE)NUMSPECTRUM); // equal spacing in log(h*nu)

    /*
      for(i=0;i<NUMSPECTRUM;i++){
      spectrum[i]=0.0;
      }
    */
    // multiply by this to convert internal energy density in cgs to computational units
    //IECONVERTDEN=pow(GRAVCONST*MASSBH*MASSSOLAR,2.0)/MASSDOT/pow(SPEEDLIGHT,3.0);

    //    SQRT3=sqrt(3.0);

    coolbottom=COOLBOTTOM;
    
    mpome=MPROTON/MELECTRON;

    compute_funcool(funcool,thetai,NUMFUNCOOL,THETAS,THETAF,&dtheta);
    odtheta=1.0/dtheta;

    thetas=THETAS;
    thetaf=THETAF;

  }


  // from Allen's Astronomical Quantities p.115-116
  // Z -> 1.06718 
  // degeneracy not important: g~Sqrt(3)/Pi*LN(4/Gamma*u) Gamma=1.781 for u<<1 u=h*nu/KT

  // Ne*Ni ~ (rho/mi)^2

  // compute cooling for each zone, subtracting off ie and putting in counter


  // compute spectrum, adding up spectrum from each zone (expensive!)

  ////////////////////////
  //
  //
  LOOPC{
    //temperature=0.5*MPROTON*(gam-1.)*s[2][k][j][i]/s[1][k][j][i]; // k*T actually
    //    ratioe=temperature/(1.0); // WTF, 1.0->h*nu, but supposed to be integrated form!(i.e. no nu)
    //    gaunt=SQRT3/M_PI*log(4.0/GAMMA/ratioe); // gaunt correction factor
    // "integrated gaunt factor" whatever that means.
    // http://www-hpcc.astro.washington.edu/papers/neal/CSTreeSPH/node9.html (missing Z^2 factor)
    // page 116 Allen's Astronomical Quantities (integrated incorrectly, or stated incorrectly since g(nu)
    // very expensive, should probably make table and linearly interp the gaunt factor
    //ftemp=(5.5-log10(temperature/KBOLTZ));
    //gaunt = 1.1+.34*exp(-ftemp*ftemp*THIRD);

    // new method

    ftemp1=s[1][k][j][i];
    // Lambda=JCONSTINT*gaunt*rho^2*sqrt(k*T)*Volume
    // ftemp=Lambda/Volume
    //coolerden=JCONSTINT*gaunt*ftemp1*ftemp1*sqrt(temperature)*dt; // energy density lost
    
    theta=(gam-1.)*s[2][k][j][i]/s[1][k][j][i]*mpome;
    if(theta>coolbottom){ // only cool if above bottom theta
      // make sure lookup table can deal with this theta, if not, recompute table
      if((theta<thetas)||(theta>thetaf)){
	
	if(theta<thetas){
	  thetas=theta*0.1;
	}
	if(theta>thetaf){
	  thetaf=theta*10.0;
	}
	compute_funcool(funcool,thetai,NUMFUNCOOL,thetas,thetaf,&dtheta);
	odtheta=1.0/dtheta;
	
	/*    
	      fprintf(fail_file,"problem with range, consider dynamic correction\n");
	      fprintf(fail_file,"theta: %15.10g thetas: %15.10g thetaf: %15.10g\n",theta,thetas,thetaf);
	      myexit(1);
	*/
      }
      
      // always do, since now corrected lookup, so always good
      cooli=(int)((theta-thetas)*odtheta);
      ratio=(theta-thetai[cooli])*odtheta;
      coolerden=coolfact*ftemp1*ftemp1*FUNCTIONCOOL(cooli,ratio)*dt; // energy density lost
      
      //radiations[2]+=IECONVERT*coolerie; // internal energy lost (in comp units)
      //s[2][k][j][i]-=IECONVERTDEN*coolerden; // internal energy density lost (in comp units)
      
      eninit=s[2][k][j][i]; // internal energy density before loss
      s[2][k][j][i]-=coolerden; // internal energy density lost
      
      volume=DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k);
      
      // floor correction here means don't radiate as much...As opposed to still saying I radiate but fill back in with floor
      if(FORCEIEINTERNAL){
	if(s[2][k][j][i]<IEFLOOR){	
	  if(DOFLOORD2==1){
	    fprintf(logfl_file,"corrected en in cooling: t: %15.10g %d %d %d bad: %15.10g\n",t,k,j,i,s[2][k][j][i]);
	  }
	  s[2][k][j][i]=IEFLOOR;
	  coolerie=(eninit-IEFLOOR)*volume;
	}
	else{
	  // normal case
	  coolerie=coolerden*volume; // full is ok then
	}
      }
      else{// have to check if no floor, or might get floating point exceptions!      
	if(s[2][k][j][i]<=0.0){
	  fprintf(fail_file,"catastrophic cooling, new ie density at %d %d %d is %15.10g due to a change of %15.10g\n",k,j,i,s[2][k][j][i],-coolerden);
	  myexit(1);
	}
	else{
	  // normal case
	  coolerie=coolerden*volume; // full is ok then
	}
      }
      
      if(accountstore[k][j][i]) radiations[2]+=coolerie; // internal energy lost
      
      
      
      /*
      // compute spectrum
      for(l=0;l<NUMSPECTRUM;l++){
      freq=pow(10.0,NUI+l*dlognu); // h*nu actually
      ratioe=temperature/freq;
      // problem with gaunt, what to choose exactly?  assume non-degenerate form
      gaunt=SQRT3/M_PI*log(4.0/GAMMA/ratioe); // gaunt correction factor(p.116 Allen's Astronomical Quantities)
      spectrum[l]=coolerie*exp(-ratioe)*PLANK/temperature; // above is integral of this
      }
      */
    }// end if not below bottom
  }
  bound(NULL,NULL,2,0,0); // bound scalar just changed
  

  //
  ///////////////////////




  firsttime=0;
}

// 1d/2d/3d valid
void compute_funcool(FTYPE *fun,FTYPE *thetap,int num,FTYPE thetas,FTYPE thetaf,FTYPE *dthetap)
{
  int i;
  FTYPE theta,dtheta;
  
  dtheta=*dthetap=(thetaf-thetas)/(FTYPE)(num-1);
  
  for(i=0;i<num+1;i++){ // extra 1 for interp to work perfectly
    
    theta=thetap[i]=thetas+i*dtheta; // theta(i)
    
    if(theta<1.0){
      
      fun[i]=4.0*sqrt(2.*theta/pow(M_PI,3.0))*(1.+1.781*pow(theta,1.34))+1.73*pow(theta,1.5)*(1.+1.1*theta+theta*theta-1.25*pow(theta,5./2.));
      
    }
    else if(theta>=1.0){
      
      fun[i]=9.0*theta/(2.*M_PI)*(log(1.123*theta+0.48)+1.5)+2.30*theta*(log(1.123*theta)+1.28);
      
    }
    
  }

  fprintf(log_file,"funcool: thetas: %15.10g thetaf: %15.10g dtheta: %15.10g\n",thetas,thetaf,dtheta);

}




// 1d/2d/3d valid
// pressure and gravity step (includes longitudinal lortenz force terms)
void step_pgc(void)
{
#if(OUTPUTFORCES)
  FTYPE grav,pressure;
  FILE *out;
  static int forcecount=0;
  static FTYPE tforcedump=0;
  char filename[200];
  int dodump;
#endif
  FTYPE rhoa,orhoa,bxa,bxp,bxm,bya,byp,bym,bza,bzp,bzm ;
  int i,j,k ;
  FTYPE tempfx1=0,tempfx2=0,tempfx3=0,tempfs=0;
  FTYPE ftemp1,ftemp2;
  FTYPE oarclen;
  FTYPE (*rho)[N2M][N1M];
  FTYPE xtilde,ytilde,ztilde,rtilde;


#if(OUTPUTFORCES)
    if(t>= tforcedump){
      dodump=1;
      tforcedump+=DTd;
    }
    else dodump=0;
#endif


  // first pressure gradient
  // x1
  if(pressx1){

#if(OUTPUTFORCES)
    if(dodump){
      sprintf(filename,"forcex%04d.dat",forcecount);
      out=fopen(filename,"wb");
      if(out==NULL){
	fprintf(stderr,"Can't open force file %d\n",forcecount);
	exit(0);
      }
    }
#endif

#if(OUTPUTFORCES)
    LOOPC
#else
    LOOPV1
#endif
      {

#if(N1>1)

#if(ALFVENLIMIT==0)
      rho=s[1];
#else
      rho=rholimited[1];
#endif
      rhoa=z2e_1(rho,k,j,i);
      orhoa=1.0/rhoa;
      oarclen=OARC11(k,j,i);

#if(PDEN)
#if(RELIE)
      // might be problem with division of linear interp
      tempfx1=-(rho[k][j][i]*s[2][k][j][i]-rho[k][j][im1]*s[2][k][j][im1])*oarclen*orhoa;
#else
#if(CSLIMIT==0)
      if(wgam) {
#if((STAREOS==BAROTROPICEOS)&&(EVOLVEENTHALPY))
	tempfx1=-gam*(s[2][k][j][i]/s[1][k][j][i]-s[2][k][j][im1]/s[1][k][j][im1])*oarclen;
#else
	tempfx1=-(gam-1.)*(s[2][k][j][i]-s[2][k][j][im1])*oarclen*orhoa;
#endif
      }
      else{
	tempfx1=-cs*cs*(rho[k][j][i]-rho[k][j][im1])*oarclen*orhoa;
      }
#else

#if(DOINVSOL2FUN)
      invsol2=invsol2fun[k][j][i];
#endif
      if(wgam) {
	tempfx1=-(gam-1.)*(s[2][k][j][i]-s[2][k][j][im1])*oarclen/(z2e_1(s[1],k,j,i)+gam*(gam-1.0)*z2e_1(s[2],k,j,i)*invsol2);
      }
      else{
	tempfx1=-cs*cs*(s[1][k][j][i]+gam*(gam-1.0)*s[2][k][j][i]*invsol2-(s[1][k][j][im1]+gam*(gam-1.0)*s[2][k][j][im1]*invsol2))*oarclen/(z2e_1(s[1],k,j,i)+gam*(gam-1.0)*z2e_1(s[2],k,j,i)*invsol2) ;
      }
#endif // CSLIMIT

#endif
#else
      tempfx1=0;
#endif


#if(RAD)
      tempfx1=tempfx1-oarclen*orhoa*fr[1][k][j][i]*(er[k][j][i]-er[k][j][im1]);
#endif

#if(OUTPUTFORCES)
      if(i>0) pressure=tempfx1;
      else pressure=0;
#endif




#if(PGRAV)

#if((GRAVACC==0)||(COORD==2))

#if(OUTPUTFORCES)
      if(i>0) grav=-(s[3][k][j][i]-s[3][k][j][im1])*oarclen;
      else grav=0;
#endif

      tempfx1+=-(s[3][k][j][i]-s[3][k][j][im1])*oarclen;  // grav from pot
#else
      // use analytic derivative(assumes PW form)


#if(OUTPUTFORCES)
      if(i>0) grav=gravacc[1][k][j][i]; // includes - sign
      else grav=0;
#endif
      tempfx1+=gravacc[1][k][j][i]; // includes - sign

#endif

#endif


#if(ROTATINGFORCE)

      // X-FACE for x-dir
      xtilde=x[1][1][i]-(pos_accretor[1]-rcm);
      ytilde=x[2][2][j]-(pos_accretor[2]);
      ztilde=x[2][3][k]-(pos_accretor[3]);
      // radius distance from this position to center of mass
      //      rtilde=sqrt( pow(xtilde,2.0) + pow(ytilde,2.0) + pow(ztilde,2.0) );
      
      tempfx1+= 2.0*v2tov1(v[1][2],k,j,i)*Omegasystem +xtilde*Omegasystem*Omegasystem;

#endif






#if(PMAG==1)
      //MARK
      if(mag==1){
	bya = v2tov1(v[2][2],k,j,i);
	bza = v3tov1(v[2][3],k,j,i);

	byp = e2z_2(v[2][2],k,j,i);
	bym = e2z_2(v[2][2],k,j,im1);
      
	bzp = e2z_3(v[2][3],k,j,i);
	bzm = e2z_3(v[2][3],k,j,im1);
      
	tempfx1 +=
	  //	-(oarclen*orhoa)*(
	  //				bya/G2(1,i)*(G2(2,i)*byp - G2(2,im1)*bym) +
	  //				bza/G3(1,i)*(G3(2,i)*bzp - G3(2,im1)*bzm)
	  // use volume differencing
	  -(orhoa*ODVL(2,1,i))*(
				bya*G3(1,i)*(G2(2,i)*byp - G2(2,im1)*bym) +
				bza*G2(1,i)*(G3(2,i)*bzp - G3(2,im1)*bzm)
				) ;
      }
#endif


#endif // endif not 1d

#if((COORD==3)&&(GRAVITOMAGNETIC==1)) // only allowed if in spc, like curvature term

      tempfx1+=blackholejz*gravitom[1][k][j][i]*G4(2,j)*z2e_1(v[1][3],k,j,i);

#endif // end if gravitomagnetic==1



#if(CURVE==1) // done even if 1-d    
#if(COORD==3)
      tempfs=v2tov1(v[1][2],k,j,i);
      tempfs*=tempfs;
      tempfx1+=tempfs/G2(1,i); // * DG2(1,i)=1
      tempfs=z2e_1(v[1][3],k,j,i);
      tempfs*=tempfs;
      tempfx1+=tempfs/G3(1,i); // * DG3(1,i)=1
#endif
#endif


#if(OUTPUTFORCES)
    if(dodump){
      if(fabs(grav)>1E-30) fprintf(out,"%3d %3d %3d %21.15g %21.15g %21.15g %21.15g %21.15g\n",i,j,k,fabs(pressure+grav)/(fabs(grav)+1E-15),pressure,grav,log(s[1][k][j][i]),log(s[2][k][j][i]));
      else fprintf(out,"%3d %3d %3d %21.15g %21.15g %21.15g %21.15g %21.15g\n",i,j,k,fabs(pressure+grav),pressure,grav,log(s[1][k][j][i]),log(s[2][k][j][i]));
    }
#endif


    
      v[1][1][k][j][i] += dt*tempfx1;
    
    }

#if(OUTPUTFORCES)
    if(dodump){
      fclose(out);
    }
#endif

  }



  if(pressx2){

#if(OUTPUTFORCES)
    if(dodump){
      sprintf(filename,"forcey%04d.dat",forcecount);
      out=fopen(filename,"wb");
      if(out==NULL){
	fprintf(stderr,"Can't open force file %d\n",forcecount);
	exit(0);
      }
    }
#endif



#if(OUTPUTFORCES)
    LOOPC
#else
    LOOPV2
#endif
      {
      // x2
#if(N2>1)

#if(ALFVENLIMIT==0)
      rho=s[1];
#else
      rho=rholimited[2];
#endif
      rhoa=z2e_2(rho,k,j,i);
      orhoa=1.0/rhoa;
      oarclen=OARC12(k,j,i);

#if(PDEN==1)
#if(RELIE)
      // might be problem with division of linear interp
      tempfx2=-(rho[k][j][i]*s[2][k][j][i]-rho[k][jm1][i]*s[2][k][jm1][i])*oarclen*orhoa;
#else

#if(CSLIMIT==0)
      if(wgam) {
#if((STAREOS==BAROTROPICEOS)&&(EVOLVEENTHALPY))
	tempfx2=-gam*(s[2][k][j][i]/s[1][k][j][i]-s[2][k][jm1][i]/s[1][k][jm1][i])*oarclen;
#else
	tempfx2=-(gam-1.)*(s[2][k][j][i]-s[2][k][jm1][i])*oarclen*orhoa;
#endif
      }
      else{
	tempfx2=-cs*cs*(rho[k][j][i]-rho[k][jm1][i])*oarclen*orhoa;
      }
#else

#if(DOINVSOL2FUN)
      invsol2=invsol2fun[k][j][i];
#endif

      if(wgam) {
	tempfx2=-(gam-1.)*(s[2][k][j][i]-s[2][k][jm1][i])*oarclen/(z2e_2(s[1],k,j,i)+gam*(gam-1.0)*z2e_2(s[2],k,j,i)*invsol2);
      }
      else{
	tempfx2=-cs*cs*(s[1][k][j][i]+gam*(gam-1.0)*s[2][k][j][i]*invsol2-(s[1][k][jm1][i]+gam*(gam-1.0)*s[2][k][jm1][i]*invsol2))*oarclen/(z2e_2(s[1],k,j,i)+gam*(gam-1.0)*z2e_2(s[2],k,j,i)*invsol2);
      }
#endif // CSLIMIT


#endif
    
#else
      tempfx2=0;
#endif
    
#if(RAD)
      tempfx2=tempfx2-oarclen*orhoa*fr[2][k][j][i]*(er[k][j][i]-er[k][jm1][i]);
#endif

#if(OUTPUTFORCES)
      if(j>0) pressure=tempfx2;
      else pressure=0;
#endif


#if(PGRAV)

#if((GRAVACC==0)||(COORD==2))

#if(OUTPUTFORCES)
      if(j>0) grav=-(s[3][k][j][i]-s[3][k][jm1][i])*oarclen;
      else grav=0;
#endif

      tempfx2 += -(s[3][k][j][i]-s[3][k][jm1][i])*oarclen; // should probably use analytical derivative
#else

#if(OUTPUTFORCES)
      if(j>0) grav=gravacc[2][k][j][i];
      else grav=0;
#endif

      tempfx2+=gravacc[2][k][j][i]; // includes - sign
#endif

#endif


#if(ROTATINGFORCE)

      // Y-FACE for y-dir
      xtilde=x[2][1][i]-(pos_accretor[1]-rcm);
      ytilde=x[1][2][j]-(pos_accretor[2]);
      ztilde=x[2][3][k]-(pos_accretor[3]);
      // radius vector from here to center of mass
      //      rtilde=sqrt( pow(xtilde,2.0) + pow(ytilde,2.0) + pow(ztilde,2.0) );

      tempfx2+= -2.0*v1tov2(v[1][1],k,j,i)*Omegasystem +ytilde*Omegasystem*Omegasystem;

#endif


        
#if(PMAG==1)
      // MARK
      if(mag==1){	
	bxa = v1tov2(v[2][1],k,j,i);
	bza = v3tov2(v[2][3],k,j,i);

	bxp = e2z_1(v[2][1],k,j,i);
	bxm = e2z_1(v[2][1],k,jm1,i);
      
	bzp = e2z_3(v[2][3],k,j,i);
	bzm = e2z_3(v[2][3],k,jm1,i);
      
	tempfx2 +=
	  //		-(oarclen*orhoa)*(
	  //		bxa*(bxp - bxm) +
	  //			bza/G4(1,j)*(G4(2,j)*bzp - G4(2,jm1)*bzm)
	  -(orhoa*ODVL(2,2,j)*OG2(2,i))*(
					 G4(1,j)*bxa*(bxp - bxm) +
					 bza*(G4(2,j)*bzp - G4(2,jm1)*bzm)
					 ) ;
      }
#endif


#endif // end if not 1d


#if((COORD==3)&&(GRAVITOMAGNETIC==1)) // only allowed if in spc, like curvature term

      tempfx2+=-2.0*blackholejz*gravitom[2][k][j][i]*DG4(1,j)*z2e_2(v[1][3],k,j,i);

#endif // end if gravitomagnetic==1



#if(CURVE==1) // done even if 1-d.
#if( (COORD==2)||(COORD==3))
      tempfs=z2e_2(v[1][3],k,j,i);
      tempfs*=tempfs;
      tempfx2+=tempfs*OG2(2,i)*OG4(1,j)*DG4(1,j);
      // true coordinate singularity problem, fixed with inverse of G4(1,j)->0, since v[1][3]->0 here.
#endif
#endif


#if(OUTPUTFORCES)
    if(dodump){
      if(fabs(grav)>1E-30) fprintf(out,"%3d %3d %3d %21.15g %21.15g %21.15g %21.15g %21.15g\n",i,j,k,fabs(pressure+grav)/(fabs(grav)+1E-15),pressure,grav,log(s[1][k][j][i]),log(s[2][k][j][i]));
      else fprintf(out,"%3d %3d %3d %21.15g %21.15g %21.15g %21.15g %21.15g\n",i,j,k,fabs(pressure+grav),pressure,grav,log(s[1][k][j][i]),log(s[2][k][j][i]));
    }
#endif

    
      v[1][2][k][j][i] += dt*tempfx2;
    
    }

#if(OUTPUTFORCES)
    if(dodump){
      fclose(out);
    }
#endif

  }

  if(pressx3){

#if(OUTPUTFORCES)
    if(dodump){
      sprintf(filename,"forcez%04d.dat",forcecount);
      out=fopen(filename,"wb");
      if(out==NULL){
	fprintf(stderr,"Can't open force file %d\n",forcecount);
	exit(0);
      }
    }
#endif

#if(OUTPUTFORCES)
    LOOPC
#else
    LOOPV3
#endif
      {
      
      // x3
#if(N3>1)


#if(ALFVENLIMIT==0)
      rho=s[1];
#else
      rho=rholimited[3];
#endif
      rhoa=z2e_3(rho,k,j,i);
      orhoa=1.0/rhoa;
      oarclen=OARC13(k,j,i)*ODX(2,3,k);

#if(PDEN==1)
#if(RELIE)
      // might be problem with division of linear interp
      tempfx3=-(rho[k][j][i]*s[2][k][j][i]-rho[km1][j][i]*s[2][km1][j][i])*oarclen*orhoa;
#else

#if(CSLIMIT==0)
      if(wgam) {
#if((STAREOS==BAROTROPICEOS)&&(EVOLVEENTHALPY))
	tempfx3=-gam*(s[2][k][j][i]/s[1][k][j][i]-s[2][km1][j][i]/s[1][km1][j][i])*oarclen;
#else
	tempfx3=-(gam-1.)*(s[2][k][j][i]-s[2][km1][j][i])*oarclen*orhoa;
#endif
      }
      else{
	tempfx3=-cs*cs*(rho[k][j][i]-rho[km1][j][i])*oarclen*orhoa;
      }
#else

#if(DOINVSOL2FUN)
      invsol2=invsol2fun[k][j][i];
#endif

      if(wgam) {
	tempfx3=-(gam-1.)*(s[2][k][j][i]-s[2][km1][j][i])*oarclen/(z2e_3(s[1],k,j,i)+gam*(gam-1.0)*z2e_3(s[2],k,j,i)*invsol2);
      }
      else{
	tempfx3=-cs*cs*(s[1][k][j][i]+gam*(gam-1.0)*s[2][k][j][i]*invsol2-(s[1][km1][j][i]+gam*(gam-1.0)*s[2][km1][j][i]*invsol2))*oarclen/(z2e_3(s[1],k,j,i)+gam*(gam-1.0)*z2e_3(s[2],k,j,i)*invsol2);
      }
#endif // CSLIMIT



#endif
    
#else
      tempfx3=0;
#endif

#if(OUTPUTFORCES)
      if(k>0) pressure=tempfx3;
      else pressure=0;
#endif
    
#if(PGRAV)


#if((GRAVACC==0)||(COORD==2))

#if(OUTPUTFORCES)
      if(k>0) grav=-(s[3][k][j][i]-s[3][km1][j][i])*oarclen;
      else grav=0;
#endif
      tempfx3 += -(s[3][k][j][i]-s[3][km1][j][i])*oarclen; // should probably use analytical derivative

#else

#if(OUTPUTFORCES)
      if(k>0) grav=gravacc[3][k][j][i];
      else grav=0;
#endif

      tempfx3+=gravacc[3][k][j][i]; // includes - sign    
#endif

#endif // end if PGRAV


     

#if(PMAG==1)
      if(mag==1){	
	bxa = v1tov3(v[2][1],k,j,i);
	bya = v2tov3(v[2][2],k,j,i);

	bxp = e2z_1(v[2][1],k,j,i);
	bxm = e2z_1(v[2][1],km1,j,i);
      
	byp = e2z_2(v[2][2],k,j,i);
	bym = e2z_2(v[2][2],km1,j,i);
      
	tempfx3 += -(oarclen*orhoa)*(
				     bxa*(bxp - bxm) +
				     bya*(byp - bym)
				     ) ;
      }
#endif

#endif // on N3>1      


      // no hydro curvature terms      

      // below can apply even when N3M==1      
#if((COORD==3)&&(GRAVITOMAGNETIC==1)) // only allowed if in spc, like curvature term
      tempfx3+= blackholejz*gravitom[2][k][j][i]*(2.0*e2z_2(v[1][2],k,j,i)*DG4(2,j)-e2z_1(v[1][1],k,j,i)*G4(2,j));
#endif // end if gravitomagnetic==1

#if(OUTPUTFORCES)
    if(dodump){
      if(fabs(grav)>1E-30) fprintf(out,"%3d %3d %3d %21.15g %21.15g %21.15g %21.15g %21.15g\n",i,j,k,fabs(pressure+grav)/(fabs(grav)+1E-15),pressure,grav,log(s[1][k][j][i]),log(s[2][k][j][i]));
      else fprintf(out,"%3d %3d %3d %21.15g %21.15g %21.15g %21.15g %21.15g\n",i,j,k,fabs(pressure+grav),pressure,grav,log(s[1][k][j][i]),log(s[2][k][j][i]));
    }
#endif


      v[1][3][k][j][i]+=dt*tempfx3;
    }
#if(OUTPUTFORCES)
    if(dodump){
      fclose(out);
    }
#endif
  }


  // bound only if changed
  
  if( ((CURVE==1)&&(COORD==3))||(N1>1) ){
    if( ((CURVE==1)&&(COORD==3))||(N2>1) ){
      if((N3>1)||(((COORD==3)&&(GRAVITOMAGNETIC==1))||(COMPDIM==3))){
	bound(NULL,NULL,0,1,123) ;
      }
      else{
	bound(NULL,NULL,0,1,12) ;
      }
    }
    else{
      if((N3>1)||(((COORD==3)&&(GRAVITOMAGNETIC==1))||(COMPDIM==3))){
	bound(NULL,NULL,0,1,13) ;
      }
      else{
	bound(NULL,NULL,0,1,1) ;
      }
    }
  }
  else{
    if( ((CURVE==1)&&(COORD==3))||(N2>1) ){
      if((N3>1)||(((COORD==3)&&(GRAVITOMAGNETIC==1))||(COMPDIM==3))){
	bound(NULL,NULL,0,1,23) ;
      }
      else{
	bound(NULL,NULL,0,1,2) ;
      }
    }
    else{
      if((N3>1)||(((COORD==3)&&(GRAVITOMAGNETIC==1))||(COMPDIM==3))){
	bound(NULL,NULL,0,1,3) ;
      }
      else{
	// nothing!
      }
    }
  }

#if(OUTPUTFORCES)
  if(dodump){
    forcecount++;
  }
#endif

}

#define CRAPOLA 1.0
#define CRAPOLA2 1.0

// 1d/2d/3d
void step_visc(void)
{
  // GODMARK(commented)
  /*
    int counter;
    #if(TIMEMETHOD==0)
    time_t timestart,timestop;
    time_t gtimestart,gtimestop;
    #elif(TIMEMETHOD==1)
    struct timeval timestart,timestop, gtimestart,gtimestop; struct timezone tz;
    #elif(TIMEMETHOD==2)
    clock_t timestart,timestop, gtimestart,gtimestop;
    #endif
    SFTYPE walltime=0,walltimelocal=0,walltot=0;
  */
  static FTYPE (*visc)[N3M][N2M][N1M];
  static FTYPE (*delv)[N2M][N1M];
  static FTYPE (*gradv)[N3M][N2M][N1M];
  static FTYPE (*l2_ten)[N2M][N1M];
  FTYPE dvx,dvy,dvz,qlx,qvnrx,qly,qvnry,qlz,qvnrz;
  FTYPE rhoa,rhob ;
  int i,j,k,l ;
  FTYPE ftemp,ftemp1,ftemp2,ftemp3;
  FTYPE ftempdelv;
  static FTYPE (*dv)[N3M][N2M][N1M];
  FTYPE fltemp11,fltemp22,fltemp33;
  FTYPE odx1,odx2,odx3,ods,odl;
  
  visc = workv1 ;
  dv = workv2;
  
  gradv = workv3 ;
  delv = work1 ;
  l2_ten = work2 ;



  // find x1-dir viscous stresses 
  LOOPVISC{ // no need to bound visc, all good
#if(VISC_TENSOR==1) // not right in 3D, but doesn't seem to work anyways in 1D or 2D

    odx1=OARC21(k,j,i);
    odx2=OARC32(k,j,i);
    if( odx1 < odx2 ){
      ods = odx2;
      odl = odx1;
    }
    else{
      ods = odx1 ;
      odl = odx2 ;
    }
#if(COMPDIM==3)
    odx3=OARC13(k,j,i)*ODX(1,3,k);
    if( ods < odx3 ){
      ods = odx3;
    }
    else if( odl > odx3 ){
      odl = odx3 ;
    }
#endif
#endif
    
    // del v, at zone center
    dvx = dv[1][k][j][i] = v[1][1][k][j][ip1] - v[1][1][k][j][i] ;
    // del v, at zone center 
    dvy = dv[2][k][j][i] = v[1][2][k][jp1][i] - v[1][2][k][j][i] ;

#if(COMPDIM==3)
    dvz = dv[3][k][j][i] = v[1][3][kp1][j][i] - v[1][3][k][j][i] ;
#else
    dvz=0;
#endif
    // linear viscosity
#if(VISC_LINEAR)
    if(wgam){
      cs = sqrt(gam*(gam-1.)*s[2][k][j][i]/s[1][k][j][i]) ;
    }
    qlx = -nu_l*s[1][k][j][i]*cs*dvx ;
    qly = -nu_l*s[1][k][j][i]*cs*dvy ;
#if(COMPDIM==3)
    qlz = -nu_l*s[1][k][j][i]*cs*dvz ;
#else
    qlz=0;
#endif


#else
    qlx=qly=qlz=0;
#endif

#if(VISC_TENSOR==0)
    // von neumann,richtmyer viscosity 
    if(dvx  < 0) {
      qvnrx = nu_vnr*s[1][k][j][i]*dvx*dvx ;
    }
    else qvnrx = 0. ;
    
    // von neumann,richtmyer viscosity 
    if(dvy  < 0) {
      qvnry = nu_vnr*s[1][k][j][i]*dvy*dvy ;
    }
    else qvnry = 0. ;

#if(COMPDIM==3)
    if(dvz  < 0) {
      qvnrz = nu_vnr*s[1][k][j][i]*dvz*dvz ;
    }
    else qvnrz = 0. ;
#else
    qvnrz=0;
#endif

#else
    // don't use this
    ftempdelv=delv[k][j][i]=deldotv(v,1,k,j,i);
    gradv[1][k][j][i]=gradv11(1,k,j,i);
    gradv[2][k][j][i]=gradv22(1,k,j,i);
    gradv[3][k][j][i]=gradv33(1,k,j,i);
    
    l2_ten[k][j][i]=nu_ten/(ods*ods);
    if(gradv[1][k][j][i]<0){
      qvnrx=CRAPOLA*l2_ten[k][j][i]*s[1][k][j][i]*ftempdelv*(gradv[1][k][j][i]-THIRD*ftempdelv);
    }
    else qvnrx=0.0;
    if(gradv[2][k][j][i]<0){
      qvnry=CRAPOLA2*l2_ten[k][j][i]*s[1][k][j][i]*ftempdelv*(gradv[2][k][j][i]-THIRD*ftempdelv);
    }
    else qvnry=0.0;
    //printf("%15.10g %15.10g %15.10g %15.10g\n",nu_vnr*s[1][k][j][i]*dvx*dvx,qvnrx, nu_vnr*s[1][k][j][i]*dvy*dvy,qvnry);
    /*
      if(dvx<0) ftemp1=nu_vnr*s[1][k][j][i]*dvx*dvx;
      else ftemp1=0;
      if(dvy<0) ftemp2=nu_vnr*s[1][k][j][i]*dvy*dvy;
      else ftemp2=0;
      
      if(fabs(ftemp1-qvnrx)>1.E-6){
      printf("s-diff1\n");
      printf("%15.10g %15.10g\n",ftemp1,qvnrx);
      exit(1);
      }
      if(fabs(ftemp2-qvnry)>1.E-6){
      printf("s-diff2\n");
      printf("%15.10g %15.10g\n",ftemp2,qvnry);
      exit(1);
      }
    */
#endif
    
    visc[1][k][j][i] = qlx + qvnrx ;
    visc[2][k][j][i] = qly + qvnry ;
    visc[3][k][j][i] = qlz + qvnrz ;
    
    //printf("%d %d %d %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",k,j,i,ftempdelv,gradv[1][k][j][i],gradv[2][k][j][i],gradv[3][k][j][i],l2_ten[k][j][i],qvnrx,qvnry);
  }
  
#define NUMGODS 100

#if(0&&(LOOPTYPE==3))
  // If LOOPTYPE==3, problem with boundary zones since outflow sets ->0 alot, never smooth, so just ignore those viscosities

  // i.e. ignore viscosity correction inside and on the boundary
  LOOPBOUNDV1in{
    visc[1][k][j][i]=0.0;
  }
  LOOPBOUNDV2in{
    visc[2][k][j][i]=0.0;
  }
  LOOPBOUNDV3in{
    visc[3][k][j][i]=0.0;
  }
  LOOPBOUNDV1out{
    visc[1][k][j][i]=0.0;
  }
  LOOPBOUNDV2out{
    visc[2][k][j][i]=0.0;
  }
  LOOPBOUNDV3out{
    visc[3][k][j][i]=0.0;
  }
#endif


  //counter=0;
  //GETTIME(&timestart);
  // update velocity, internal energy 
  if(ie&&(STAREOS==ADIABATICEOS)){
    if(wgam) {
      //for(l=1;l<=NUMGODS;l++){
      //LOOPSUPERGEN(5){
      LOOPC{
	//i=indx[5][0];j=indx[5][1];k=indx[5][2];
	//for(temptempi=0;temptempi<numiter[5];temptempi++){
	// i=indx[5][temptempi*3];j=indx[5][temptempi*3+1];k=indx[5][temptempi*3+2];
	//LOOP{
	//if(bzmask[k][j][i]!=0) continue;
	//if(bzmask[k][j][i]!=-600) continue;
	//if(bzmask[k][j][i]==-600) continue;

	//	counter++;

#if(VISC_TENSOR==0)
	s[2][k][j][i] += -dt*(
			      visc[1][k][j][i]*dv[1][k][j][i]*OARC21(k,j,i)
			      + visc[2][k][j][i]*dv[2][k][j][i]*OARC32(k,j,i)
#if(COMPDIM==3)
			      + visc[3][k][j][i]*dv[3][k][j][i]*OARC13(k,j,i)*ODX(1,3,k)
#endif
	
			      // GODMARK(commented)
			      /*
				+ visc[3][k][j-1][i]*dv[3][k][j][i]*OARC13(k,j,i)*ODX(1,3,k)
				+ visc[3][k][j][i-1]*dv[3][k][j][i]*OARC13(k,j,i)*ODX(1,3,k)
				+ visc[3][k-1][j][i]*dv[3][k][j][i]*OARC13(k,j,i)*ODX(1,3,k)
				+ visc[3][k][j+1][i]*dv[3][k][j][i]*OARC13(k,j,i)*ODX(1,3,k)
				+ visc[3][k][j][i+1]*dv[3][k][j][i]*OARC13(k,j,i)*ODX(1,3,k)
				+ visc[3][k+1][j][i]*dv[3][k][j][i]*OARC13(k,j,i)*ODX(1,3,k)
				+ visc[3][k+1][j+1][i]*dv[3][k][j][i]*OARC13(k,j,i)*ODX(1,3,k)
				+ visc[3][k][j+1][i+1]*dv[3][k][j][i]*OARC13(k,j,i)*ODX(1,3,k)
				+ visc[3][k+1][j][i+1]*dv[3][k][j][i]*OARC13(k,j,i)*ODX(1,3,k)
			      */
			      );
      
#else// don't use
	ftemp1= -dt*(visc[1][k][j][i]*dv[1][k][j][i]*OARC21(k,j,i) + visc[2][k][j][i]*dv[2][k][j][i]*OARC32(k,j,i) ) ;
	
	ftempdelv=delv[k][j][i];
	fltemp11=gradv[1][k][j][i];
	fltemp22=gradv[2][k][j][i];
	fltemp33=gradv[3][k][j][i];
	ftemp2=- CRAPOLA*dt*(l2_ten[k][j][i]*s[1][k][j][i]*ftempdelv*THIRD*(
									    ((fltemp11-fltemp22)*(fltemp11-fltemp22))+
									    ((fltemp11-fltemp33)*(fltemp11-fltemp33))+
									    ((fltemp33-fltemp22)*(fltemp33-fltemp22)) ) );
	s[2][k][j][i] += ftemp2;
	/*
	  if(fabs(ftemp1-ftemp2)>1.E-6){
	  printf("s-en-diff1\n");
	  printf("%15.10g %15.10g\n",ftemp1,ftemp2);
	  printf("%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",ftempdelv,dv[1][k][j][i],visc[1][k][j][i],fltemp11,fltemp22,fltemp33);
	  exit(1);
	  }
	*/
	
#endif
	
	//}
      }
      if(FORCEIEINTERNAL){
	floor_correct(1,0);
      }
    }
  }
  /*
    GETTIME(&timestop); walltime=(SFTYPE) DELTATIME(timestop,timestart);
    fprintf(stdout,"step: %ld counted: %d walltime: %g\n",nstep,counter/NUMGODS,walltime);  fflush(stdout);
    fprintf(stdout,"floorcount: %d\n",floorcnt[0][2]);  fflush(stdout);
    //myexit(0);
    */  
  if(N1>1){
    LOOPV1{
#if(VISC_TENSOR==0)
      
      v[1][1][k][j][i] += -dt*(visc[1][k][j][i]-visc[1][k][j][im1])*OARC11(k,j,i)/z2e_1(s[1],k,j,i) ;
      
#else
      ftemp1= -dt*(visc[1][k][j][i]-visc[1][k][j][im1])*OARC11(k,j,i)/z2e_1(s[1],k,j,i) ;
      ftemp2= -dt*(G2(2,i)*G2(2,i)*G3(2,i)*visc[1][k][j][i]-G2(2,im1)*G2(2,im1)*G3(2,im1)*visc[1][k][j][im1])/(G2(1,i)*G2(1,i)*G3(1,i)*dx[2][1][i]*z2e_1(s[1],k,j,i)) ;
      v[1][1][k][j][i] += ftemp2; 
      /*
	if(fabs(ftemp1-ftemp2)>1.E-6){
	printf("s-v1-diff1\n");
	printf("%15.10g %15.10g\n",ftemp1,ftemp2);
	exit(1);
	}
      */
#endif
      
    }
  }

  if(N2>1){
    LOOPV2{
      
#if(VISC_TENSOR==0)
      
      v[1][2][k][j][i] += -dt*(visc[2][k][j][i]-visc[2][k][jm1][i])*OARC12(k,j,i)/z2e_2(s[1],k,j,i) ;
      
#else
      ftemp1=-dt*(visc[2][k][j][i]-visc[2][k][jm1][i])*OARC12(k,j,i)/z2e_2(s[1],k,j,i) ;
      
      ftemp2=-(dt/z2e_2(s[1],k,j,i))*(  (G4(2,j)*G4(2,j)*visc[2][k][j][i]-G4(2,jm1)*G4(2,jm1)*visc[2][k][jm1][i])/(G2(2,i)*dx[2][2][j]*G4(1,j)*G4(1,j)) + 0.5*(visc[1][k][j][i]-visc[1][k][jm1][i])/(G2(2,i)*G4(1,j))*DG4(1,j) ) ;
      v[1][2][k][j][i] += ftemp2;
      /*
	if(fabs(ftemp1-ftemp2)>1.E-6){
	printf("s-v2-diff1\n");
	printf("%15.10g %15.10g\n",ftemp1,ftemp2);
	exit(1);
	}
      */
#endif
      
    }
  }
  if(N3>1){
    LOOPV3{
      
      v[1][3][k][j][i] += -dt*(visc[3][k][j][i]-visc[3][km1][j][i])*OARC13(k,j,i)*ODX(2,3,k)/z2e_3(s[1],k,j,i) ;
      
    }
  }
  
  // bound everything
  if( (ie)&&(wgam)){
    if(N1>1){
      if(N2>1){
	if(N3>1){
	  bound(NULL,NULL,2,1,123) ; // both velocities and ie (step_visc)
	}
	else{
	  bound(NULL,NULL,2,1,12) ; // both velocities and ie (step_visc)
	}
      }
      else{
	if(N3>1){
	  bound(NULL,NULL,2,1,13) ; // both velocities and ie (step_visc)
	}
	else{
	  bound(NULL,NULL,2,1,1) ; // both velocities and ie (step_visc)
	}
      }
    }
    else{
      if(N2>1){
	if(N3>1){
	  bound(NULL,NULL,2,1,23) ; // both velocities and ie (step_visc)
	}
	else{
	  bound(NULL,NULL,2,1,2) ; // both velocities and ie (step_visc)
	}
      }
      else{
	if(N3>1){
	  bound(NULL,NULL,2,1,3) ; // both velocities and ie (step_visc)
	}
	else{
	  // nothing!
	}
      }
    }
  }
  else{
    if(N1>1){
      if(N2>1){
	if(N3>1){
	  bound(NULL,NULL,0,1,123) ; // both velocities and ie (step_visc)
	}
	else{
	  bound(NULL,NULL,0,1,12) ; // both velocities and ie (step_visc)
	}
      }
      else{
	if(N3>1){
	  bound(NULL,NULL,0,1,13) ; // both velocities and ie (step_visc)
	}
	else{
	  bound(NULL,NULL,0,1,1) ; // both velocities and ie (step_visc)
	}
      }
    }
    else{
      if(N2>1){
	if(N3>1){
	  bound(NULL,NULL,0,1,23) ; // both velocities and ie (step_visc)
	}
	else{
	  bound(NULL,NULL,0,1,2) ; // both velocities and ie (step_visc)
	}
      }
      else{
	if(N3>1){
	  bound(NULL,NULL,0,1,3) ; // both velocities and ie (step_visc)
	}
	else{
	  // nothing!
	}
      }
    }
  }
  
}



// 1d/2d/3d valid
// must be consistent with nu_fact calc in init.c
void nu_compute(void)
{
  int i,j,k;
  FTYPE ftemp;
  FTYPE ftemp2;
  FTYPE cs2;

  ftemp=sqrt(gam*(gam-1.));
  ftemp2=gam*(gam-1.);

  // no need to bound since can do LOOPFC since nu_real only depends on zone centered quantities

  if(vreal==1){    
    LOOPFC{
      // alpha-param
      if(wgam) cs = ftemp*sqrt(s[2][k][j][i]/s[1][k][j][i]);
      nu_real[k][j][i]=alpha_reala[k][j][i]*nu_fact[k][j][i]*cs;
    }
  }
  else if(vreal==2){
    LOOPFC{
      if(wgam) cs2 = ftemp2*s[2][k][j][i]/s[1][k][j][i];
      else cs2=cs*cs;
      //Igumenshchev alpha param: $\nu=\alpha*cs^2/Omega_{k}$
      nu_real[k][j][i]=alpha_reala[k][j][i]*nu_fact[k][j][i]*cs2;
    }
  }
  else if(vreal==3){
    LOOPFC{
      // Stone et al.
      // make consistent with nu_fact in init.c
      // no need for nu_fact since just constant for now
      
      // run A through I'
      nu_real[k][j][i]=alpha_reala[k][j][i]*s[1][k][j][i]*nu_fact[k][j][i];
      
      // run J
      //nu_real[k][j][i]=alpha_real*s[1][k][j][i]*nu_fact[k][j][i];
      
      // run K
      //nu_real[k][j][i]=alpha_real*nu_fact[k][j][i];
    }
  }
  else if(vreal==4){
    LOOPFC{
      if(wgam) cs2 = ftemp2*s[2][k][j][i]/s[1][k][j][i];
      else cs2=cs*cs;
      //Igumenshchev alpha param: $\nu=\alpha*cs^2/Omega_{k}$
      nu_real[k][j][i]=alpha_reala[k][j][i]*nu_fact[k][j][i]*cs2;
    }
  }
  else if(vreal==5){
    LOOPFC{
      nu_real[k][j][i]=alpha_reala[k][j][i]*nu_fact[k][j][i];
    }
  }
  if(analoutput==6){
    LOOPFC{
      // for test of this code
      nu_real[k][j][i]=nu_fact[k][j][i];
    }
  }

}


/*
  I don't compute momentum boundary flux due to v_r and v_theta components since they involve volume terms for spc
*/

// 1d/2d/3d valid
void step_visc_real(void)
{
  static FTYPE (*nurho_real)[N2M][N1M];

  static FTYPE (*delv)[N2M][N1M]; // deldotv
  int i,j,k,l,m ;
  FTYPE ftemp,ftemp1,ftemp2;
  FTYPE subftemp;
  FTYPE Length;
  FTYPE flux;

  FTYPE crapf;
  
  delv = work1;
  nurho_real = work2;

  compute_sigma_gen(sigma,rost,rostnu,nurho_real,delv);
  
  // update velocity, internal energy 
  //crapf=0;
  if(vischeat&&(STAREOS==ADIABATICEOS)){
    if((ie)||(analoutput==6)){
      if(wgam) {
	LOOPC{
	  ftemp=0;
	  // this is some-what expensive, but I optimized for velocity terms and sigma, not this internal energy term!
	  
#if(VISCE11)
	  ftemp+=rost[1][1][k][j][i]*rost[1][1][k][j][i];
#endif
	  
#if(VISCE22)
	  ftemp+=rost[2][2][k][j][i]*rost[2][2][k][j][i];
#endif
	  
#if(VISCE33)
	  ftemp+=rost[3][3][k][j][i]*rost[3][3][k][j][i];
#endif
	  
#if(VISCE13)
	  subftemp=c2z_2(rost[1][3],k,j,i);
	  //#if(ANALOUTPUT==6)
	  //v[1][2][k][j][i]=sigma[1][3][k][j][i];
	  //#endif      
	  ftemp+=2.0*(subftemp*subftemp); // 2 accounts for rost[3][1]
#endif
	  
#if(VISCE23)
	  subftemp=c2z_1(rost[2][3],k,j,i);
	  ftemp+=2.0*(subftemp*subftemp); // 2 accounts for rost[3][2]
#endif
	  
#if(VISCE12)
	  subftemp=c2z_3(rost[1][2],k,j,i);
	  ftemp+=2.0*(subftemp*subftemp); // 2 accounts for rost[2][1]
#endif
	  
	  s[2][k][j][i] += dt*nurho_real[k][j][i]*ftemp;
	  
	}
	if(FORCEIEINTERNAL){
	  floor_correct(2,1);
	}
      }
    } //endif doing ie
  } //endif vischeat==1

  
  if(VISCE11||VISCE12||VISCE22||VISCE33){
    LOOPV1{
#if(COORD>=2)
	
      v[1][1][k][j][i] +=
	-dt/(z2e_1(s[1],k,j,i))*(
				 (G2(2,i)*G3(2,i)*sigma[1][1][k][j][i]-G2(2,im1)*G3(2,im1)*sigma[1][1][k][j][im1])/DVL(2,1,i)
				 +(G4(1,jp1)*sigma[1][2][k][jp1][i]-G4(1,j)*sigma[1][2][k][j][i])/(G2(1,i)*DVL(1,2,j))
#if(COMPDIM==3)
				 +(sigma[1][3][kp1][j][i]-sigma[1][3][k][j][i])*OARC23(k,j,i)*ODX(1,3,k)
#endif
				 -(z2e_1(sigma[2][2],k,j,i)+z2e_1(sigma[3][3],k,j,i))/G2(1,i)
				 );
	
#elif(COORD==1)
	
      v[1][1][k][j][i] +=
	-dt/(z2e_1(s[1],k,j,i))*(
				 (sigma[1][1][k][j][i]-sigma[1][1][k][j][im1])*ODX(2,1,i)
				 +(sigma[1][2][k][jp1][i]-sigma[1][2][k][j][i])*ODX(1,2,j)
#if(COMPDIM==3)
				 +(sigma[1][3][kp1][j][i]-sigma[1][3][k][j][i])*ODX(1,3,k)
#endif
				 );
	
#endif
    }
  }

  if(VISCE12||VISCE22||VISCE33){
    LOOPV2{
#if(COORD>=2)
	
      v[1][2][k][j][i] +=
	-dt/(z2e_2(s[1],k,j,i))*(
				 (G2(1,ip1)*G2(1,ip1)*G3(1,ip1)*sigma[1][2][k][j][ip1]-G2(1,i)*G2(1,i)*G3(1,i)*sigma[1][2][k][j][i])/(G2(2,i)*DVL(1,1,i))
				 +(G4(2,j)*sigma[2][2][k][j][i]-G4(2,jm1)*sigma[2][2][k][jm1][i])/(G2(2,i)*DVL(2,2,j))
#if(COMPDIM==3)
				 +(sigma[3][2][kp1][j][i]-sigma[3][2][k][j][i])*OARC33(k,j,i)*ODX(1,3,k)
#endif
				 -z2e_2(sigma[3][3],k,j,i)*DG4(1,j)*OG4(1,j)*OG2(2,i)
				 // true coordinate singularity problem, fixed with OG4(1,j)->0 at poles, since sigma[3][3]->0 here really
				 );      
#elif(COORD==1)
	
      v[1][2][k][j][i] +=
	-dt/(z2e_2(s[1],k,j,i))*(
				 +(sigma[1][2][k][j][ip1]-sigma[1][2][k][j][i])*ODX(1,1,i)
				 +(sigma[2][2][k][j][i]-sigma[2][2][k][jm1][i])*ODX(2,2,j)
#if(COMPDIM==3)
				 +(sigma[3][2][kp1][j][i]-sigma[3][2][k][j][i])*ODX(1,3,k)
#endif
				 );      
#endif
    }
  }

  if(VISCE13||VISCE23){  
    LOOPV3{
#if(COORD>=2)
      v[1][3][k][j][i] +=
	-dt/(z2e_3(s[1],k,j,i))*(
				 (G2(1,ip1)*G3(1,ip1)*G3(1,ip1)*sigma[1][3][k][j][ip1]-G2(1,i)*G3(1,i)*G3(1,i)*sigma[1][3][k][j][i])/(G3(2,i)*DVL(1,1,i))
				 +(G4(1,jp1)*G4(1,jp1)*sigma[2][3][k][jp1][i]-G4(1,j)*G4(1,j)*sigma[2][3][k][j][i])/(G2(2,i)*G4(2,j)*DVL(1,2,j))
#if(COMPDIM==3)
				 +(sigma[3][3][k][j][i]-sigma[3][3][km1][j][i])*OARC13(k,j,i)*ODX(2,3,k)
#endif
				 );
	
#elif(COORD==1)
      v[1][3][k][j][i] +=
	-dt/(z2e_3(s[1],k,j,i))*(
				 (sigma[1][3][k][j][ip1]-sigma[1][3][k][j][i])*ODX(1,1,i)
				 +(sigma[2][3][k][jp1][i]-sigma[2][3][k][j][i])*ODX(1,2,j)
#if(COMPDIM==3)
				 +(sigma[3][3][k][j][i]-sigma[3][3][km1][j][i])*ODX(2,3,k)
#endif
				 );
#endif
    }
  }
  
  // bound everything
  if( ((ie)&&(wgam))||(analoutput==6)){
    bound(NULL,NULL,2,1,123) ; // both velocities and ie (step_visc)
  }
  else bound(NULL,NULL,0,1,123); // otherwise just velocities (step_visc)

}


// 1d/2d/3d valid
// assume nu_real, s, v, and geom terms have already been computed as per timestep or initially
void compute_sigma_3(FTYPE (*sigma)[3][N3M][N2M][N1M],FTYPE (*rost)[3][N3M][N2M][N1M],FTYPE (*rostnu)[3][N3M][N2M][N1M],FTYPE (*nurho_real)[N2M][N1M],FTYPE (*delv)[N2M][N1M])
{
  int i,j,k,l,m ;
  FTYPE ftemp,ftemp1,ftemp2;
  FTYPE subftemp;
  FTYPE Length;
  FTYPE flux;

  FTYPE crapf;
  
  // get important coefficient(nu*rho)
  LOOPFC{ // needs full loop for interp for eij on half-full loop    
    nurho_real[k][j][i]= 2.0*s[1][k][j][i]*nu_real[k][j][i] ; // 2*rho*nu really
  }
  
  // find sigma=-2*rho*nu*e_{ij}
  LOOPHC{// ok even if periodicx2special==1
    delv[k][j][i]=deldotv(v,1,k,j,i); // no need to interp since this centered and used to compute only centered sigma's


#if(VISCE11)
    ftemp=( (v[1][1][k][j][ip1]-v[1][1][k][j][i])*OARC21(k,j,i)-THIRD*delv[k][j][i] );
    rost[1][1][k][j][i]=ftemp;
    rostnu[1][1][k][j][i]=ftemp*nu_real[k][j][i];
    sigma[1][1][k][j][i]=-nurho_real[k][j][i]*ftemp;
#endif
    
#if(VISCE22)
    ftemp=( (v[1][2][k][jp1][i]-v[1][2][k][j][i])*OARC32(k,j,i)+e2z_1(v[1][1],k,j,i)/G2(2,i) - THIRD*delv[k][j][i] );
    rost[2][2][k][j][i]=ftemp;
    rostnu[2][2][k][j][i]=ftemp*nu_real[k][j][i];
    sigma[2][2][k][j][i]=-nurho_real[k][j][i]*ftemp;
#endif
    
#if(VISCE33)
    ftemp=( ( e2z_1(v[1][1],k,j,i)+e2z_2(v[1][2],k,j,i)*DG4(2,j)/G4(2,j))/G3(2,i) - THIRD*delv[k][j][i] );
#if(COMPDIM==3)
    ftemp+=OARC13(k,j,i)*ODX(1,3,k)*(v[1][3][kp1][j][i]-v[1][3][k][j][i]);
#endif
    rost[3][3][k][j][i]=ftemp;
    rostnu[3][3][k][j][i]=ftemp*nu_real[k][j][i];
    sigma[3][3][k][j][i]=-nurho_real[k][j][i]*ftemp;
#endif
    
#if(VISCE12)
    ftemp=0.5*( G2(1,i)*(v[1][2][k][j][i]/G2(2,i)-v[1][2][k][j][im1]/G2(2,im1))*OARC31(k,j,i)+ (v[1][1][k][j][i]-v[1][1][k][jm1][i])*OARC22(k,j,i) ) ;
    rost[1][2][k][j][i]=rost[2][1][k][j][i]=ftemp;
    rostnu[1][2][k][j][i]=rostnu[2][1][k][j][i]=ftemp*z2c_3(nu_real,k,j,i);
    sigma[1][2][k][j][i]=sigma[2][1][k][j][i]=-z2c_3(nurho_real,k,j,i)*ftemp;
#endif
    
#if(VISCE13)
    ftemp=0.5*(G3(1,i)*(v[1][3][k][j][i]/G3(2,i)-v[1][3][k][j][im1]/G3(2,im1))*OARC11(k,j,i) );
#if(COMPDIM==3)
    ftemp+=OARC23(k,j,i)*ODX(2,3,k)*(v[1][1][k][j][i]-v[1][1][km1][j][i]);
#endif
    rost[1][3][k][j][i]=rost[3][1][k][j][i]=ftemp;
    rostnu[1][3][k][j][i]=rostnu[3][1][k][j][i]=ftemp*z2e_1(nu_real,k,j,i);
    sigma[1][3][k][j][i]=sigma[3][1][k][j][i]=-z2e_1(nurho_real,k,j,i)*ftemp;
#endif
    
#if(VISCE23)
    ftemp=(0.5*G4(1,j)*OARC12(k,j,i))*( v[1][3][k][j][i]/G4(2,j)-v[1][3][k][jm1][i]/G4(2,jm1)) ;
#if(COMPDIM==3)
    ftemp+=OARC33(k,j,i)*ODX(2,3,k)*(v[1][2][k][j][i]-v[1][2][km1][j][i]);
#endif
    rost[2][3][k][j][i]=rost[3][2][k][j][i]=ftemp;
    rostnu[2][3][k][j][i]=rostnu[3][2][k][j][i]=ftemp*z2e_2(nu_real,k,j,i);
    sigma[2][3][k][j][i]=sigma[3][2][k][j][i]=-z2e_2(nurho_real,k,j,i)*ftemp;
#endif

    //#if(ANALOUTPUT==6)
    //v[1][1][k][j][i]=sigma[1][3][k][j][i]; // for testing real visc
    //#endif
    

  }// end e_{ij} generation

}



// 1d/2d/3d valid
// assume nu_real, s, v, and geom terms have already been computed as per timestep or initially
void compute_sigma_1(FTYPE (*sigma)[3][N3M][N2M][N1M],FTYPE (*rost)[3][N3M][N2M][N1M],FTYPE (*rostnu)[3][N3M][N2M][N1M],FTYPE (*nurho_real)[N2M][N1M],FTYPE (*delv)[N2M][N1M])
{
  int i,j,k,l,m ;
  FTYPE ftemp,ftemp1,ftemp2;
  FTYPE subftemp;
  FTYPE Length;
  FTYPE flux;

  FTYPE crapf;
  
  // get important coefficient(nu*rho)
  LOOPFC{ // needs full loop for interp for eij on half-full loop    
    nurho_real[k][j][i]= 2.0*s[1][k][j][i]*nu_real[k][j][i] ; // 2*rho*nu really
  }
  
  // find sigma=-2*rho*nu*e_{ij}
  LOOPHC{
    delv[k][j][i]=deldotv(v,1,k,j,i); // no need to interp since this centered and used to compute only centered sigma's


#if(VISCE11)
    ftemp=( (v[1][1][k][j][ip1]-v[1][1][k][j][i])*OARC21(k,j,i)-THIRD*delv[k][j][i] );
    rost[1][1][k][j][i]=ftemp;
    rostnu[1][1][k][j][i]=ftemp*nu_real[k][j][i];
    sigma[1][1][k][j][i]=-nurho_real[k][j][i]*ftemp;
#endif
    
#if(VISCE22)
    ftemp=( (v[1][2][k][jp1][i]-v[1][2][k][j][i])*OARC32(k,j,i) - THIRD*delv[k][j][i] );
    rost[2][2][k][j][i]=ftemp;
    rostnu[2][2][k][j][i]=ftemp*nu_real[k][j][i];
    sigma[2][2][k][j][i]=-nurho_real[k][j][i]*ftemp;
#endif
    
#if(VISCE33)
    ftemp=( - THIRD*delv[k][j][i]) ;
#if(COMPDIM==3)
    ftemp+=OARC13(k,j,i)*ODX(1,3,k)*(v[1][3][kp1][j][i]-v[1][3][k][j][i]);
#endif
    rost[3][3][k][j][i]=ftemp;
    rostnu[3][3][k][j][i]=ftemp*nu_real[k][j][i];
    sigma[3][3][k][j][i]=-nurho_real[k][j][i]*ftemp;
#endif
    
#if(VISCE12)
    ftemp=0.5*( (v[1][2][k][j][i]-v[1][2][k][j][im1])*OARC31(k,j,i)+ (v[1][1][k][j][i]-v[1][1][k][jm1][i])*OARC22(k,j,i) ) ;
    rost[1][2][k][j][i]=rost[2][1][k][j][i]=ftemp;
    rostnu[1][2][k][j][i]=rostnu[2][1][k][j][i]=ftemp*z2c_3(nu_real,k,j,i);
    sigma[1][2][k][j][i]=sigma[2][1][k][j][i]=-z2c_3(nurho_real,k,j,i)*ftemp;
#endif
    
#if(VISCE13)
    ftemp=0.5*((v[1][3][k][j][i]-v[1][3][k][j][im1])*OARC11(k,j,i) );
#if(COMPDIM==3)
    ftemp+=OARC23(k,j,i)*ODX(2,3,k)*(v[1][1][k][j][i]-v[1][1][km1][j][i]);
#endif
    rost[1][3][k][j][i]=rost[3][1][k][j][i]=ftemp;
    rostnu[1][3][k][j][i]=rostnu[3][1][k][j][i]=ftemp*z2e_1(nu_real,k,j,i);
    sigma[1][3][k][j][i]=sigma[3][1][k][j][i]=-z2e_1(nurho_real,k,j,i)*ftemp;
#endif
    
#if(VISCE23)
    ftemp=(0.5*OARC12(k,j,i))*( v[1][3][k][j][i]-v[1][3][k][jm1][i]) ;
#if(COMPDIM==3)
    ftemp+=OARC33(k,j,i)*ODX(2,3,k)*(v[1][2][k][j][i]-v[1][2][km1][j][i]);
#endif
    rost[2][3][k][j][i]=rost[3][2][k][j][i]=ftemp;
    rostnu[2][3][k][j][i]=rostnu[3][2][k][j][i]=ftemp*z2e_2(nu_real,k,j,i);
    sigma[2][3][k][j][i]=sigma[3][2][k][j][i]=-z2e_2(nurho_real,k,j,i)*ftemp;
#endif

    //#if(ANALOUTPUT==6)
    //v[1][1][k][j][i]=sigma[1][3][k][j][i]; // for testing real visc
    //#endif
    

  }// end e_{ij} generation

}




// 1d/2d/3d valid
void step_ie(void)
{
  FTYPE dv,dvg ;
  int i,j,k ;
  FTYPE ftemp;


  // update internal energy 
  LOOPC{
#if(STAREOS==ADIABATICEOS)
    dv =deldotv(v,1,k,j,i);
    dvg=0.5*dt*(gam-1.0)*dv;
    s[2][k][j][i] *= (1. - dvg)/(1. + dvg) ;

    //s[2][k][j][i]+=dt*(gam-1.)*s[2][k][j][i]*dv;
#elif(STAREOS==BAROTROPICEOS)
    s[2][k][j][i]=KAPPASTAR*pow(s[1][k][j][i],gam)/(gam-1.0);
#endif

  }
  if(FORCEIEINTERNAL){
    floor_correct(2,2);
  }

  bound(NULL,NULL,2,0,0) ;

}

// initial attempt of theta range
#define THETARS (0.0)
#define THETARF (300.0)

// linear interp
#define FUNCTIONR1(i,ratio) (funrelie[i]*(1.-ratio)+funrelie[ip1]*ratio)


// 1d/2d/3d valid
void step_relie (void)
{
  FTYPE ftemp;
  FTYPE dv;
  int i,j,k ;
  static int firsttime=1;
  static FTYPE theta,thetaguess,dtheta,odtheta,thetas,thetaf;
  int reli;
  FTYPE ratio;
  FTYPE fun1,fun2;

  if(firsttime==1){

    compute_funrelie(funrelie,thetareli,NUMFUNRELIE,THETARS,THETARF,&dtheta);

    odtheta=1.0/dtheta;
    thetas=THETAS;
    thetaf=THETAF;
  }

  // update internal energy 
  LOOPC{
    dv =deldotv(v,1,k,j,i);
    theta=s[2][k][j][i];

    // make sure lookup table can deal with this theta, if not, recompute table
    if((theta<thetas)||(theta>thetaf)){
      
      if(theta<thetas){
	thetas=theta*0.1;
      }
      if(theta>thetaf){
	thetaf=theta*10.0;
      }
      compute_funrelie(funrelie,thetareli,NUMFUNRELIE,thetas,thetaf,&dtheta);      
      odtheta=1.0/dtheta;
    }

    // always do, since now corrected lookup, so always good
    reli=(int)((theta-thetas)*odtheta);
    ratio=(theta-thetareli[reli])*odtheta;

    fun1=FUNCTIONR1(reli,ratio);

    // first guess at new theta
    thetaguess=theta-dt*dv*fun1;

    // make sure lookup table can deal with this theta, if not, recompute table
    // just cut paste of above
    if((thetaguess<thetas)||(thetaguess>thetaf)){
      
      if(thetaguess<thetas){
	thetas=thetaguess*0.1;
      }
      if(thetaguess>thetaf){
	thetaf=thetaguess*10.0;
      }
      compute_funrelie(funrelie,thetareli,NUMFUNRELIE,thetas,thetaf,&dtheta);      
      odtheta=1.0/dtheta;
    }

    // compute new functionr1 based on guess
    reli=(int)((thetaguess-thetas)*odtheta);
    ratio=(thetaguess-thetareli[reli])*odtheta;
    fun2=FUNCTIONR1(reli,ratio);

    // now correct
    s[2][k][j][i]=theta-dt*dv*0.5*(fun1+fun2);

  }

  if(FORCEIEINTERNAL){
    floor_correct(2,2);
  }

  bound(NULL,NULL,2,0,0) ;
  firsttime=0;
}


// 1d/2d/3d valid
// need to make tables for step_relie and advection
void compute_funrelie(FTYPE *fun,FTYPE *thetap,int num,FTYPE thetas,FTYPE thetaf,FTYPE *dthetap)
{
  int i;
  FTYPE theta,dtheta;
  //FTYPE gtheta;
  
  dtheta=*dthetap=(thetaf-thetas)/(FTYPE)(num-1);
  
  for(i=0;i<num+1;i++){ // extra 1 for interp to work perfectly
    
    theta=thetap[i]=thetas+i*dtheta; // theta(i)
    
    //gtheta=(12.+45.*theta+45.*theta*theta)/(8.+20.*theta+15.*theta*theta);

    // theta/D[c^2*t*gtheta,theta]
    fun[i]=theta*pow(8.+20.*theta+15.*theta*theta,2.0)/(3.*(32.+240.*theta+600.*theta*theta+600.*theta*theta*theta+225.*theta*theta*theta*theta));
  }

  fprintf(log_file,"funrelie: thetas: %15.10g thetaf: %15.10g dtheta: %15.10g\n",thetas,thetaf,dtheta);

}

// 1d/2d/3d valid
// compute 1/c^2 functional dependence
void compute_invsol2fun(void)
{
  int i,j,k;

#if(DOINVSOL2FUN==0)
  // invsol2 never changes from invsol2_original
  //  LOOPF{ // no dependence on anything
  //    invsol2fun[k][j][i]=invsol2_original;
  //  }
#elif(DOINVSOL2FUN==1)
  LOOPF{ // no dependence on anything
    if(s[1][k][j][i]>1E4*DENSITYFLOOR) invsol2fun[k][j][i]=1.1E-7; // true speed of light for averystar problem
    else invsol2fun[k][j][i]=invsol2_original;
  }
#endif
  

}



void fracfloor_correct(void)
{

  int i,j,k,l;
  FTYPE ftemp,mins[2],maxs[2];
  FTYPE minssend[2],maxssend[2];
  FTYPE DENSITYFLOORSEND,IEFLOORSEND;
  
  mins[0]=mins[1]=1E30;
  maxs[0]=maxs[1]=0.0;

  LOOPC{
    for(l=0;l<=1;l++){
      ftemp=s[l+1][k][j][i];
      if(ftemp<mins[l]) mins[l]=ftemp;
      if(ftemp>maxs[l]) maxs[l]=ftemp;
    }
  }
  if(numprocs>1){
    minssend[0]=mins[0];
    minssend[1]=mins[1];
    maxssend[0]=maxs[0];
    maxssend[1]=maxs[1];
#if(USEMPI)
    MPI_Allreduce(&minssend[0], &mins[0], 1, MPI_FTYPE, MPI_MIN, MPI_COMM_WORLD);            
    MPI_Allreduce(&minssend[1], &mins[1], 1, MPI_FTYPE, MPI_MIN, MPI_COMM_WORLD);            
    MPI_Allreduce(&maxssend[0], &maxs[0], 1, MPI_FTYPE, MPI_MAX, MPI_COMM_WORLD);            
    MPI_Allreduce(&maxssend[1], &maxs[1], 1, MPI_FTYPE, MPI_MAX, MPI_COMM_WORLD);            
#endif
  }
  // now have min and max densities so can set floor value

  // GODMARK  override this for now
  if(0){
    DENSITYFLOOR=maxs[0]*DENSITYFRACFLOOR;
    IEFLOOR=maxs[1]*IEFRACFLOOR;
  }

}



// set the minimum value of density
// wloc=0 for t=0
void floor_correct(int wsca,int wloc)
{
  int i,j,k;
  FTYPE ftemp1,ftemp2;
  FTYPE v1xa,v1ya,v1za;
  FTYPE ftemp0;
  FTYPE ftemp;


  // only account for accountable region, but checks entire computational grid (bzones included)

  if(wsca==1){

    LOOPFLOOR{ // must avoid non-comp. zones

      if(s[1][k][j][i]<DENSITYFLOOR){

	if(wloc>=0){
	
	  ftemp1=(DENSITYFLOOR-s[1][k][j][i]);
	  ftemp2=ftemp1*s[3][k][j][i];
	  
	  
	  v1xa=e2z_1(v[1][1],k,j,i);
	  v1ya=e2z_2(v[1][2],k,j,i);
	  v1za=e2z_3(v[1][3],k,j,i);
	  // floor ke located at zone center
	  ftemp0=0.5*(DENSITYFLOOR-s[1][k][j][i])*(v1xa*v1xa+v1ya*v1ya+v1za*v1za);

	  if(accountstore[k][j][i]) floors[1]+=ftemp1*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k) ;
	  if(accountstore[k][j][i]) floors[3]+=ftemp2*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k) ;
	  if(accountstore[k][j][i]) floors[NUMSCA+1]+=ftemp0*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k) ;
	  
	  if(FLOORDUMPFLAG==1){
	    floorvars[1][k][j][i]+=ftemp1;
	    floorvars[3][k][j][i]+=ftemp2;
	    floorvar0[1][k][j][i]+=ftemp0;
	  }
	  
	  if(DOFLOORDIAG==1){
	    floorcnt[wloc][1]++;
	    if(s[1][k][j][i]<floorlowest[1]){
	      floorlowest[1]=s[1][k][j][i];
	      wherelowest[1]=wloc;
	    }
	  }
	}
	if(DOFLOORD2==1){
	  fprintf(logfl_file,"corrected rho: t: %15.10g %d %d %d %15.10g\n",t,k,j,i,(DENSITYFLOOR-s[1][k][j][i]));
	}
	s[1][k][j][i] = DENSITYFLOOR;
      }
    }
  }



  if(wsca==2){

    LOOPFLOOR{ // must avoid non-comp. zones
      
      if(s[2][k][j][i]<IEFLOOR){
	
	if(wloc>=0){
	  ftemp=(IEFLOOR-s[2][k][j][i]);
	  if(accountstore[k][j][i]) floors[2]+=ftemp*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k) ;
	  if(FLOORDUMPFLAG==1){
	    floorvars[2][k][j][i]+=ftemp;
	  }
	  if(DOFLOORDIAG==1){
	    floorcnt[wloc][2]++;
	    if(s[2][k][j][i]<floorlowest[2]){
	      floorlowest[2]=s[2][k][j][i];
	      wherelowest[2]=wloc;
	    }
	  }
	}
	if(DOFLOORD2==1){
	  fprintf(logfl_file,"corrected en: t: %15.10g %d %d %d %15.10g\n",t,k,j,i,(IEFLOOR-s[2][k][j][i]));
	}
	s[2][k][j][i]=IEFLOOR;
      }
      
    }
  }




}


// purpose is for optimital gravity terms in speed and quality
void initialize_gravity(  FTYPE (*sca3)[N3M][N2M][N1M], FTYPE (*vx3)[N3M][N2M][N1M],  FTYPE (*vy3)[N3M][N2M][N1M],  FTYPE (*vz3)[N3M][N2M][N1M])
{
  int i,j,k;
  SFTYPE min[3+1],max[3+1];
  FTYPE (*maskdiff)[N2M][N1M];


  maskdiff=work1;
  


  if(!POSTPROC){
    min[1]=min[2]=min[3]=1E30;
    max[1]=max[2]=max[3]=-1E30;
    
    if(GRAVACC){
      LOOPF{
	maskdiff[k][j][i]=0;
      }
      LOOPFC{
	maskdiff[k][j][i]=1;
      }
      // now know what is null space or not (could have used mask if set them to memory in init_bc_gen()
      
      LOOPFC{
	// gravacc is only used inside LOOPV1+LOOPV2+LOOPV3, so for LOOPTYPE==3 where we scan over possible bad regions, let's project out the rest, which is at leat all scalar-like boundary zones.
	if((COORD==3)){
	  gravacc[1][k][j][i]=-GRAVC*MASSBH/(pow(vx3[1][k][j][i]-rgp,2.0));
	  gravacc[2][k][j][i]=0;
	  gravacc[3][k][j][i]=0;
	}
	else if(COORD==1){
	  gravacc[1][k][j][i]=-GRAVC*MASSBH/(pow(vx3[1][k][j][i]-rgp,2.0))*x[1][1][i]/vx3[1][k][j][i];
	  gravacc[2][k][j][i]=-GRAVC*MASSBH/(pow(vy3[1][k][j][i]-rgp,2.0))*x[1][2][j]/vy3[1][k][j][i];
	  gravacc[3][k][j][i]=-GRAVC*MASSBH/(pow(vz3[1][k][j][i]-rgp,2.0))*x[1][3][k]/vz3[1][k][j][i];
	}
	if(gravacc[1][k][j][i]>max[1]) max[1]=gravacc[1][k][j][i];
	if(gravacc[2][k][j][i]>max[2]) max[2]=gravacc[2][k][j][i];
	if(gravacc[3][k][j][i]>max[3]) max[3]=gravacc[3][k][j][i];
	
	if(gravacc[1][k][j][i]<min[1]) min[1]=gravacc[1][k][j][i];
	if(gravacc[2][k][j][i]<min[2]) min[2]=gravacc[2][k][j][i];
	if(gravacc[3][k][j][i]<min[3]) min[3]=gravacc[3][k][j][i];
      }
      if(LOOPTYPE==3){
	LOOPF{
	  // below not a computational issue, just a numerical issue to avoid nan in non-used regions with LOOPTYPE==3
	  if(maskdiff[k][j][i]==0){
	    gravacc[1][k][j][i]=min[1];
	    gravacc[2][k][j][i]=min[2];
	    gravacc[3][k][j][i]=min[3];
	  }
	}
      }
    }


    if(GRAVITOMAGNETIC){
      LOOPF{
	if(COORD==3){
	  // holds 1/r^3 for gravitomagnetic terms
	  gravitom[1][k][j][i]=1.0/(pow(vx3[1][k][j][i],3.0)); // a-grid
	  gravitom[2][k][j][i]=1.0/(pow(vy3[1][k][j][i],3.0)); // b-grid
	}
	else if(COORD==1){
	  // not yet
	}
      }
    }


  }// end if not postproc
}
