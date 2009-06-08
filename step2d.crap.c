#include "step.h"


// add source terms due to grid expansion and contraction 

void stepvar_2d(void)
{


  fracfloor_correct(); // adjust floor to give a set dynamic range in density

  // source steps

  if(MOC2DVER==0){
    transmagx1=1;
    transmagx2=1;
    transx1=transx2=1; // must do for alfven wave advection
  }
  else{
    transmagx1=0;
    transmagx2=0;
  }


  if(subcyclen<=1){

#if(MDOTMEM)
    if(mdotin){
      injection_2d();
    }
#endif
    if(cool){
      cooling();
    }

    // substep 1
    if(ALFVENLIMIT==1){
      magprepare(); // need new rho (rho*) for gas pressure and magnetic pressure terms
      // rho doesn't change between here and transport
      // B doesn't change between here and moc
    }
    if(press){
      step_pgc() ; //includes curvature terms from div(rho*v*v)
    }
    if(mag==1){
      if(MOC2DVER==0){
	if(ALFVENLIMIT==0){
	  magprepare(); // rho doesn't change between here and transport
	}
	step_bz_2d() ; // toroidal field evolution
	if(ALFVENLIMIT==1) magprepare(); // since bz changed
      }
      else{ // must make sure this is off, do this instead of in sweep.c having a condition
	transmagx1=0;
	transmagx2=0;
      }
    }
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
    }
  }

  if(subcyclen<=1){
    // substep 2
    if(visc_art==1){
      step_visc() ;
    }
    
    
    // substep 3
    if(ie){
      if(RELIE){
	step_relie();
      }
      else{
	if(wgam) step_ie() ;
      }
    }

    // CT and Alfven wave source terms
    if(mag==1){
      // magprepare used as above
      if((MOC2DVER==1)&&(SPLITMETHOD==1)){
	lorentz_3d();
      }
      if(stepmocct){
	if(RESMEM&&(res_real==1)){
	  current_compute(123); // need full current for abs(J) in nu_res_compute(rreal==2) || step_res
	  nu_res_compute();
	  step_res() ; // needs all 3 currents, and don't recompute current, just use old current that was used here if 3D.  If 2D then emf3 uses old current
	  if(ALFVENLIMIT==1) magprepare(); // since bz changed
	}
	if(MOC2DVER==0) moc_ct_2d() ; // uses only 3 component of current
	else if((MOC2DVER==1)&&(SPLITMETHOD==0)) moc_ct_3d_v1();
	else if((MOC2DVER==1)&&(SPLITMETHOD==1)) moc_ct_3d_v2();
	// b changed, but wait till timestep to magprepare if alfvenlimit==1
      }
      if(DOLOSSDIAG){
	magnetic_flux();
      }
    }
  }
  
  if(subcyclen<=1){    
    // transport steps
    if(trans==1){
      step_trans_2d() ;
    }
  }
  t += dt ;
}


// only 2d

// mass/etc injection routine

// inflows[] only takes what injected on accountable grid, not outside.  So this doesn't match physical injection rate, which is the true reference
void injection_2d(void)
{
  int i,j,k;
  static int firsttime=1;
  static FTYPE (*rhoi)[N2M][N1M];
  static FTYPE (*rhof)[N2M][N1M];
  FTYPE vxa,vya,vza;
  FTYPE rhoiv[3+1],rhofv[3+1],massi,massf,die,den,dpot;
  FTYPE drhov[3+1],kei[3+1],kef[3+1],dpd[3+1];
  FTYPE dxdyc,dxdy1,dxdy2;
  FTYPE drho, dmass, volume;
  short storeit;
  FTYPE realtemp1,realtemp2;

  rhoi=work1;
  rhof=work2;
  ////////////////////////
  // compute mass inflow as function of time over grid
  LOOPRINJ{ // so no need to bound (LOOPFC would work too)
    // can't use LOOPFCINJ because need rhoi,rhof -1,-1 from original LOOPFCINJ modulo the BC
    storeit=accountstore[k][j][i];

    // below 2 only defined once, needed over whole grid due to interpolation for velocities
    rhoi[k][j][i]=s[1][k][j][i];
    drho=rhoinject[k][j][i]*dt;
    rhof[k][j][i]=rhoi[k][j][i]+drho; // PRECISION

    volume=DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k);
    dmass=drho*volume;

    // below 2 only defined once
    //massi=rhoi[k][j][i]*volume;
    //massf=massi+dmass;
    
    // mass
    if(storeit) inflows[1]+=dmass;
    s[1][k][j][i]+=drho; // PRECISION
    //      fprintf(fail_file,"%d %d %d %15.10g\n",k,j,i,drho);
    /*
    // only check where there should exist some mass added
    if( (i>=tagii)&&(i<tagfi)&&(j>=tagij)&&(j<tagfj)&&(k>=tagik)&&(k<tagfk) ){
      if(rhoi[k][j][i]/drho>1.1E7){
	fprintf(stderr,"problem: rhoi: %21.15g drho: %21.15g dmass: %21.15g inflows[1]: %21.15g\n",rhoi[k][j][i],drho,dmass,inflows[1]);
      }
      fprintf(stdout,"dmass: %21.15g inflows[1]: %21.15g\n",dmass,inflows[1]);
    }  
    */    
  
    //ie
    //iei=s[2][k][j][i]*DVL(1,1,i)*DVL(1,2,j)*DVL(1,3,k);
    //die=IEFRACT*fabs(dmass*s[3][k][j][i]);
    den=eninject[k][j][i]*dt;
    die=den*volume;
    if(storeit) inflows[2]+=die;
    s[2][k][j][i]+=den;
    
    // pot energy
    //poti=massi*s[3][k][j][i];
    dpot=dmass*s[3][k][j][i];
    if(storeit) inflows[3]+=dpot;
    
    //fprintf(stdout,"cpu%02d %d %d %d %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",myid,k,j,i,drho,dmass,massi,massf,inflows[1],s[1][k][j][i],die,inflows[2],s[2][k][j][i],dpot,inflows[3]);
    //fflush(stdout);
    
  }
  // must bound now if multiple cpus since data across boundaries has changed(only when loop had unconstrained inflows[] could I not do LOOPFC, with LOOPFC, all boundary terms are modified appropriately, and no interp needed so no bad accesses)
  //if(numprocs>1){
  //bound(NULL,NULL,-1,0,0); // bound scalars just changed
  //}
  //    for(k=tagik;k<tagfk;k++) for(j=tagij;j<tagfj;j++) for(i=tagii;i<tagfi;i++){
  //LOOPC{ // can't do LOOPFC here due to interp, so need to bound v's later
  LOOPVINJ{ // LOOP would work here too if use LOOPFC/LOOPH above

    storeit=accountstore[k][j][i];
    
    dxdyc=DVL(1,1,i)*DVL(1,2,j) ;
    dxdy1=DVL(2,1,i)*DVL(1,2,j) ;
    dxdy2=DVL(1,1,i)*DVL(2,2,j) ;
    
    
    
    // for v's
    drhov[3]=rhoinject[k][j][i]*dt;
    drhov[2]=z2e_2(rhoinject,k,j,i)*dt;
    drhov[1]=z2e_1(rhoinject,k,j,i)*dt;
    
    // zone centered for ke only
    massi=rhoi[k][j][i]*dxdyc;
    massf=massi+drhov[3]*dxdyc; // since v3 zone centered
    
    // for v's
    rhoiv[3]=rhoi[k][j][i];
    rhofv[3]=rhof[k][j][i];
    rhoiv[2]=z2e_2(rhoi,k,j,i);
    rhofv[2]=z2e_2(rhof,k,j,i);
    rhoiv[1]=z2e_1(rhoi,k,j,i);
    rhofv[1]=z2e_1(rhof,k,j,i);
    
    // initial ke
    vxa = e2z_1(v[1][1],k,j,i);      
    vya = e2z_2(v[1][2],k,j,i);
    vza = v[1][3][k][j][i] ;  // When fake-3d
    
    kei[1]=0.5*massi*vxa*vxa;
    kei[2]=0.5*massi*vya*vya;
    kei[3]=0.5*massi*vza*vza;
    
    // vz
    dpd[3]=VZFRACT*drhov[3]*pow(G3(2,i),-0.5);// s3=r*v
    if(storeit) inflows[7]+=G3(2,i)*G4(2,j)*dpd[3]*dxdyc; // r*sin(theta)*m*v(angmom3)
    v[1][3][k][j][i]=(rhoiv[3]*v[1][3][k][j][i]+dpd[3])/rhofv[3]; // cons of L, but r*sint same
    
    // vy
    dpd[2]=0.0;// s2=r*v (interpolated to vx2s edge)
    if(storeit) inflows[8]+=G3(2,i)*dpd[2]*dxdy2; // r*m*v (angmom2)
    v[1][2][k][j][i]=(rhoiv[2]*v[1][2][k][j][i]+dpd[2])/rhofv[2]; // cons of L, but r same
    
    
    // vx
    dpd[1]=0.0;// s1=r*v (interpolated to vx1s edge)
    if(storeit) inflows[9]+=dpd[1]*dxdy1; // m*v(angmom1)
    v[1][1][k][j][i]=(rhoiv[1]*v[1][1][k][j][i]+dpd[1])/rhofv[1]; // cons of L, but r same
    
    //      fprintf(fail_file,"%d %d %d %15.10g %15.10g\n",k,j,i,rhoiv[2]/rhofv[2],rhoiv[1]/rhofv[1]);
    
    // ke
    // compute final values
    vxa = e2z_1(v[1][1],k,j,i);      
    vya = e2z_2(v[1][2],k,j,i);
    vza = v[1][3][k][j][i] ;  // When fake-3d
    
    kef[1]=0.5*massf*vxa*vxa;
    kef[2]=0.5*massf*vya*vya;
    kef[3]=0.5*massf*vza*vza;
    
    if(storeit) inflows[NUMSCA+1]+=(kef[3]-kei[3])+(kef[2]-kei[2])+(kef[1]-kei[1]);
    
  }
  // have to bound vectors if multiple cpus, since data changed across boundaries
  if(numprocs>1){ // could even isolate bound to necessary transfers only
    bound(NULL,NULL,0,1,123); // just velocity for now
  }
  //    fflush(fail_file);

  //
  ///////////////////////


  firsttime=0;
}

// only for 2d
// evolve toroidal field comp x1
void bzsweepx_2d(void)
{
  FTYPE fraczone,dqminus,dqplus;
  static FTYPE (*bstar)[N2M][N1M],(*vstar)[N2M][N1M],(*dqv)[N2M][N1M],(*dqb)[N2M][N1M] ;
  static FTYPE (*omegav)[N2M][N1M],(*magflux)[N2M][N1M],(*osrhoh3)[N2M][N1M];
  static FTYPE bm,bp,vm,vp ;
  static FTYPE sgn_va ;
  static FTYPE dqvm,dqvp,dqbm,dqbp,pr,va,vap,vam ;
  register int i,j,k ;
  FTYPE dist;
  FTYPE ftemp0,ftemp1;
  FTYPE ftemp;
#if(RHOINTERP==0)
  FTYPE D1a,srhoh3a,srhoa;
#elif(RHOINTERP==1)
  FTYPE D1m,D1p,srhoh3m,srhoh3p;
#endif



  vstar = work1 ;	
  bstar = work2 ;
  dqv = work3 ;
  dqb = work4 ;


#if(SNCODE)
  omegav=v[1][3];
  magflux=v[2][3];
  osrhoh3=osqrtrho;
#else
  omegav=work5;
  magflux=work6;
  osrhoh3=work7;
  // compute rotational velocity and magnetic flux
  LOOPC3 LOOPC2 LOOPH1{ // only need these from i=-1..N1,j=0..N2-1
    ftemp1=1.0/(G3(2,i)*G4(2,j));
    omegav[k][j][i]=v[1][3][k][j][i]*ftemp1;
    magflux[k][j][i]=v[2][3][k][j][i]/ftemp1;
    osrhoh3[k][j][i]=osqrtrho[k][j][i]*ftemp1*ftemp1; // 1/(sqrt(rho)*h3*h3)
  }
#endif

  // find vanleer slopes for omegav,magflux
  
  dqx_calc(omegav,dqv); // same as dqvx wcom=3
  dqx_calc(magflux,dqb); // same as dqvx wcom=3

  // vstar,bstar, located at zone boundary
  // bz,vz: first calculate vstar, bstar
  LOOPT1i{ // needs to go j=0..N2-1 and i=0..N1 for below, which is ok given LOOPH for dqv above since no dqv[ip1] is needed here(dqv[N1+1=N1M] doesn't exist from above)

#if(RHOINTERP==0)
    srhoa=z2e_1(osqrtrho,k,j,i); // 1/sqrt(rho) for center
    srhoh3a=z2e_1(osrhoh3,k,j,i); // 1/(sqrt(rho)*h3^2) for center
#elif(RHOINTERP==1)
    srhoh3p=osrhoh3[k][j][im1];
    srhoh3m=osrhoh3[k][j][i]; // 1/(sqrt(rho)*h3^2)
#endif
    ftemp=v[2][1][k][j][i];
    sgn_va = copysign(1.,ftemp) ;
    ftemp=fabs(ftemp)*dt;
#if(RHOINTERP==0)
    D1a = ftemp*srhoa*OARC11(k,j,i) ; // fraction of zone off center
#else
    D1p= ftemp*osqrtrho[k][j][im1]*OARC11(k,j,i);
    D1m= ftemp*osqrtrho[k][j][i]*OARC11(k,j,i);
#endif

    dist=dx[2][1][i];
    
    // values at the foot of the plus characteristic
    //  are, by convention, in im1 zone
#if(RHOINTERP==0)
    // D1p=D1m=D1a
#if(SYMFORCEMAG==0)
    vp = omegav[k][j][im1]  + (1.0-D1a)*dqv[k][j][im1]*dist ;
    bp = magflux[k][j][im1] + (1.0-D1a)*dqb[k][j][im1]*dist ;
    vm = omegav[k][j][i]    + (-1.0+D1a)*dqv[k][j][i]*dist ;
    bm = magflux[k][j][i]   + (-1.0+D1a)*dqb[k][j][i]*dist ;
#else
    if(D1a>DTOTLIMIT){
      vp = omegav[k][j][im1]  + (1.0-D1a)*dqv[k][j][im1]*dist ;
      bp = magflux[k][j][im1] + (1.0-D1a)*dqb[k][j][im1]*dist ;
      vm = omegav[k][j][i]    + (-1.0+D1a)*dqv[k][j][i]*dist ;
      bm = magflux[k][j][i]   + (-1.0+D1a)*dqb[k][j][i]*dist ;
    }
    else{
#if(SYMFUNKED==0)
      vp = vm = 0.5*(omegav[k][j][im1]  + (1.0-D1a)*dqv[k][j][im1]*dist+
                     omegav[k][j][i]    + (-1.0+D1a)*dqv[k][j][i]*dist);
      bp = bm = 0.5*(magflux[k][j][im1] + (1.0-D1a)*dqb[k][j][im1]*dist+
                     magflux[k][j][i]   + (-1.0+D1a)*dqb[k][j][i]*dist);  
#else
      // this is a little different than sweep.c
      fraczone=D1a*INVDTOTLIMIT;

      dqminus=omegav[k][j][im1]  + (1.0-DTOTLIMIT)*dqv[k][j][im1]*dist ;
      dqplus=omegav[k][j][i]    + (-1.0+DTOTLIMIT)*dqv[k][j][i]*dist ;
      vm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      vp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));

      dqminus=magflux[k][j][im1] + (1.0-DTOTLIMIT)*dqb[k][j][im1]*dist ;
      dqplus=magflux[k][j][i]   + (-1.0+DTOTLIMIT)*dqb[k][j][i]*dist ;
      bm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      bp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
    // no strong sym forced
#if(SYMFORCEMAG==0)
      vp = omegav[k][j][im1]  + (1.0-D1p)*dqv[k][j][im1]*dist ;
      bp = magflux[k][j][im1] + (1.0-D1p)*dqb[k][j][im1]*dist ;
      vm = omegav[k][j][i]    + (-1.0+D1m)*dqv[k][j][i]*dist ;
      bm = magflux[k][j][i]   + (-1.0+D1m)*dqb[k][j][i]*dist ;
#else
    if(D1p>DTOTLIMIT){
      vp = omegav[k][j][im1]  + (1.0-D1p)*dqv[k][j][im1]*dist ;
      bp = magflux[k][j][im1] + (1.0-D1p)*dqb[k][j][im1]*dist ;
    }
    else{
      vp = 0.5*(omegav[k][j][im1]  + (1.0-D1p)*dqv[k][j][im1]*dist+
		omegav[k][j][i]    + (-1.0+D1m)*dqv[k][j][i]*dist);
      bp = 0.5*(magflux[k][j][im1] + (1.0-D1p)*dqb[k][j][im1]*dist+
		magflux[k][j][i]   + (-1.0+D1m)*dqb[k][j][i]*dx);
    }
    if(D1m>DTOTLIMIT){
      vm = omegav[k][j][i]    + (-1.0+D1m)*dqv[k][j][i]*dist ;
      bm = magflux[k][j][i]   + (-1.0+D1m)*dqb[k][j][i]*dist ;
    }
    else{
      vm = 0.5*(omegav[k][j][im1]  + (1.0-D1p)*dqv[k][j][im1]*dist+
		omegav[k][j][i]    + (-1.0+D1m)*dqv[k][j][i]*dist);
      bm = 0.5*(magflux[k][j][im1] + (1.0-D1p)*dqb[k][j][im1]*dist+
		magflux[k][j][i]   + (-1.0+D1m)*dqb[k][j][i]*dx);
    }
#endif // end if force sym

#endif // end if rhointerp=1
    // solution to Stone & Norman eqtn 43,44-- use constant rho
#if(RHOINTERP==0)
    bstar[k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoh3a) ; // bstar is SN's Phistar
    vstar[k][j][i] = 0.5*(vm + vp + sgn_va*(bm - bp)*srhoh3a) ; // vstar is SN's Omegastar
#elif(RHOINTERP==1)
    bstar[k][j][i] = (sgn_va*(vm-vp)+bp*srhoh3p+bm*srhoh3m)/(srhoh3m+srhoh3p);
    vstar[k][j][i] = 0.5*(vm+vp+sgn_va*bstar[k][j][i]*(srhoh3p-srhoh3m) +sgn_va*(-bp*srhoh3p+bm*srhoh3m) );
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif
    
  }

  // (coord fixed)
  // step bz,vz (eq 45,46)
  // ZEUS-2D code does this different, but ZEUS-2D paper is same as this
  // need bstar/vstar on i=0..N1,j=0..N2-1
  LOOPC{
    ftemp0=dt*e2z_1(v[2][1],k,j,i); // S&N's bx in code
    ftemp1=ftemp0*OARC11(k,j,i); // S&N's bx/dista(i) in code
#if(ALFVENLIMIT==0)
    ftemp=s[1][k][j][i];
#else
    ftemp=rholimited[0][k][j][i]; // rho+b^2/c^2
#endif

#if(TORVX1)
#if(SNCODE)
    // as S&N code
    v[1][3][k][j][i] += ftemp1/(ftemp*G3(2,i))*(G3(1,ip1mac(i))*bstar[k][j][ip1] - G3(1,i)*bstar[k][j][i]) ;
    //v3(i,j) = v3(i,j) + bx/d(i,j)*(g31a(ip1mac(i))*b3star(ip1)-g31a(i)*b3star(i))/(g31b(i)*dista(i)) 
#else
    // as S&N paper
    v[1][3][k][j][i] += ftemp1/(ftemp*G3(2,i))*(bstar[k][j][ip1] - bstar[k][j][i]) ;
#endif
#endif
#if(TORBX1)
#if(SNCODE)
    // as S&N code
    v[2][3][k][j][i] += ftemp1*(vstar[k][j][ip1] - vstar[k][j][i])-ftemp0*e2z_1(vstar,k,j,i)*DG3(2,i)/G3(2,i) ;
#else
    // as S&N paper
    v[2][3][k][j][i] += ftemp1*G3(2,i)*(vstar[k][j][ip1] - vstar[k][j][i]) ;
#endif
#endif
  }
  bound(NULL,NULL,0,-1,3) ; // really only need to bound v[1][3] and v[2][3].
}

// only for 2d
// evolve toroidal field comp x2
void bzsweepy_2d(void)
{
  FTYPE fraczone,dqminus,dqplus;
  static FTYPE (*bstar)[N2M][N1M],(*vstar)[N2M][N1M],(*dqv)[N2M][N1M],(*dqb)[N2M][N1M];
  static FTYPE (*omegav)[N2M][N1M],(*magflux)[N2M][N1M],(*osrhoh3)[N2M][N1M];
  static FTYPE bm,bp,vm,vp ;
  static FTYPE sgn_va;
  static FTYPE dqvm,dqvp,dqbm,dqbp,pr,va ;
  register int i,j,k ;
  FTYPE dist;
  FTYPE ftemp0,ftemp1;
  FTYPE ftemp;
#if(RHOINTERP==0)
  FTYPE D1a,srhoa,srhoh3a;
#elif(RHOINTERP==1)
  FTYPE D1m,D1p,srhoh3m,srhoh3p;
#endif



  vstar = work1 ;	
  bstar = work2 ;
  dqv = work3 ;
  dqb = work4 ;



#if(SNCODE)
    omegav=v[1][3];
    magflux=v[2][3];
    osrhoh3=osqrtrho;
#else
  // compute rotational velocity and magnetic flux
  omegav=work5;
  magflux=work6;
  osrhoh3=work7;
  LOOPC3 LOOPH2 LOOPC1{ // only need these from i=0..N1-1,j=-1..N2
    ftemp1=1.0/(G3(2,i)*G4(2,j));
    omegav[k][j][i]=v[1][3][k][j][i]*ftemp1;
    magflux[k][j][i]=v[2][3][k][j][i]/ftemp1;
    osrhoh3[k][j][i]=osqrtrho[k][j][i]*ftemp1*ftemp1; // 1/(sqrt(rho)*h3*h3)
  }
#endif


  // find vanleer slopes for omegav,magflux
  
  dqy_calc(omegav,dqv); // same as dqvy wcom=3
  dqy_calc(magflux,dqb); // same as dqvy wcom=3

  // vstar,bstar, located at zone boundary
  // bz,vz: first calculate vstar, bstar
  LOOPT1j{ // see x-version, but applies to jp1 instead here
    // need rho, etc: i=0..N1-1,j=-1..N2
#if(RHOINTERP==0)
    srhoa=z2e_2(osqrtrho,k,j,i); // 1/sqrt(rho) for center
    srhoh3a=z2e_2(osrhoh3,k,j,i); // 1/(sqrt(rho)*h3^2) for center
#elif(RHOINTERP==1)
    srhoh3m=osrhoh3[k][j][i]; // 1/(sqrt(rho)*h3^2)
    srhoh3p=osrhoh3[k][jm1][i];
#endif

    ftemp=v[2][2][k][j][i];
    sgn_va = copysign(1.,ftemp) ;
    ftemp=fabs(ftemp)*dt;

#if(RHOINTERP==0)
    D1a = ftemp*srhoa*OARC12(k,j,i) ;
#else
    D1p= ftemp*osqrtrho[k][jm1][i]*OARC12(k,j,i);
    D1m= ftemp*osqrtrho[k][j][i]*OARC12(k,j,i);
#endif

    dist=G2(2,i)*dx[2][2][j];    



    // values at the foot of the plus characteristic
    //  are, by convention, in jm1 zone
#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    vp = omegav[k][jm1][i]  + (1.0-D1a)*dqv[k][jm1][i]*dist ;
    bp = magflux[k][jm1][i] + (1.0-D1a)*dqb[k][jm1][i]*dist ;
    vm = omegav[k][j][i]    + (-1.0+D1a)*dqv[k][j][i]*dist ;
    bm = magflux[k][j][i]   + (-1.0+D1a)*dqb[k][j][i]*dist ;
#else
    if(D1a>DTOTLIMIT){
      vp = omegav[k][jm1][i]  + (1.0-D1a)*dqv[k][jm1][i]*dist ;
      bp = magflux[k][jm1][i] + (1.0-D1a)*dqb[k][jm1][i]*dist ;
      vm = omegav[k][j][i]    + (-1.0+D1a)*dqv[k][j][i]*dist ;
      bm = magflux[k][j][i]   + (-1.0+D1a)*dqb[k][j][i]*dist ;
    }
    else{
#if(SYMFUNKED==0)
      vp = vm = 0.5*(omegav[k][jm1][i]  + (1.0-D1a)*dqv[k][jm1][i]*dist+
                     omegav[k][j][i]    + (-1.0+D1a)*dqv[k][j][i]*dist);
      bp = bm = 0.5*(magflux[k][jm1][i] + (1.0-D1a)*dqb[k][jm1][i]*dist+
                     magflux[k][j][i]   + (-1.0+D1a)*dqb[k][j][i]*dist); 
#else
      // this is a little different than sweep.c
      fraczone=D1a*INVDTOTLIMIT;

      dqminus=omegav[k][jm1][i]  + (1.0-DTOTLIMIT)*dqv[k][jm1][i]*dist ;
      dqplus=omegav[k][j][i]    + (-1.0+DTOTLIMIT)*dqv[k][j][i]*dist ;
      vm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      vp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));

      dqminus=magflux[k][jm1][i] + (1.0-DTOTLIMIT)*dqb[k][jm1][i]*dist ;
      dqplus=magflux[k][j][i]   + (-1.0+DTOTLIMIT)*dqb[k][j][i]*dist ;
      bm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      bp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
    // no strong sym forced
#if(SYMFORCEMAG==0)
    vp = omegav[k][jm1][i]  + (1.0-D1p)*dqv[k][jm1][i]*dist ;
    bp = magflux[k][jm1][i] + (1.0-D1p)*dqb[k][jm1][i]*dist ;
    vm = omegav[k][j][i]    + (-1.0+D1m)*dqv[k][j][i]*dist ;
    bm = magflux[k][j][i]   + (-1.0+D1m)*dqb[k][j][i]*dist ;
#else
    if(D1p>DTOTLIMIT){
      vp = omegav[k][jm1][i]  + (1.0-D1p)*dqv[k][jm1][i]*dist ;
      bp = magflux[k][jm1][i] + (1.0-D1p)*dqb[k][jm1][i]*dist ;
    }
    else{
      vp = 0.5*(omegav[k][jm1][i]  + (1.0-D1p)*dqv[k][jm1][i]*dist+
		omegav[k][j][i]    + (-1.0+D1m)*dqv[k][j][i]*dist);
      bp = 0.5*(magflux[k][jm1][i] + (1.0-D1p)*dqb[k][jm1][i]*dist+
		magflux[k][j][i]   + (-1.0+D1m)*dqb[k][j][i]*dist);
    }
    if(D1m>DTOTLIMIT){
      vm = omegav[k][j][i]    + (-1.0+D1m)*dqv[k][j][i]*dist ;
      bm = magflux[k][j][i]   + (-1.0+D1m)*dqb[k][j][i]*dist ;
    }
    else{
      vm = 0.5*(omegav[k][jm1][i]  + (1.0-D1p)*dqv[k][jm1][i]*dist+
		omegav[k][j][i]    + (-1.0+D1m)*dqv[k][j][i]*dist);
      bm = 0.5*(magflux[k][jm1][i] + (1.0-D1p)*dqb[k][jm1][i]*dist+
		magflux[k][j][i]   + (-1.0+D1m)*dqb[k][j][i]*dist);
    }

#endif

#endif

    
#if(RHOINTERP==0)
    // solution to Stone & Norman eqtn 43,44-- use constant rho 
    bstar[k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoh3a) ; // bstart is SN's Phistar
    vstar[k][j][i] = 0.5*(vm + vp + sgn_va*(bm - bp)*srhoh3a) ; // vstar is SN's Omegastar
#elif(RHOINTERP==1)
    bstar[k][j][i] = (sgn_va*(vm-vp)+bp*srhoh3p+bm*srhoh3m)/(srhoh3m+srhoh3p);
    vstar[k][j][i] = 0.5*(vm+vp+sgn_va*bstar[k][j][i]*(srhoh3p-srhoh3m) +sgn_va*(-bp*srhoh3p+bm*srhoh3m) );
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif

  }

  // coord fixed
  // step bz,vz(eq 47, 48)
  // ZEUS-2D code does this different, but ZEUS-2D paper is same as this
  LOOPC{    // need b*,v*: i=0..N1-1, j=0..N2
    ftemp0=dt*e2z_2(v[2][2],k,j,i); // S&N code's bx
    ftemp1=ftemp0*OARC32(k,j,i);
#if(ALFVENLIMIT==0)
    ftemp=s[1][k][j][i];
#else
    ftemp=rholimited[0][k][j][i]; // rho+b^2/c^2
#endif

#if(TORVX2)
#if(SNCODE)
    // S&N code
    v[1][3][k][j][i] += ftemp0/(ftemp*G2(2,i)*DVL(1,2,j))*(G4(1,jp1mac(j))*bstar[k][jp1][i] - G4(1,j)*bstar[k][j][i]);
#else
    // S&N paper
    v[1][3][k][j][i] += ftemp1/(ftemp*G4(2,j))*(bstar[k][jp1][i] - bstar[k][j][i]);
#endif
#endif

#if(TORBX2)
#if(SNCODE)
    // S&N code
    v[2][3][k][j][i] += ftemp1*(vstar[k][jp1][i] - vstar[k][j][i])-ftemp0*e2z_2(vstar,k,j,i)*DG4(2,j)/(G2(2,i)*G4(2,j)) ;
#else
    // S&N paper
    v[2][3][k][j][i] += ftemp1*G4(2,j)*(vstar[k][jp1][i] - vstar[k][j][i]) ;
#endif
#endif
  }
  bound(NULL,NULL,0,-1,3) ; // really only need to bound v[1][3] and v[2][3].
}

// step z magnetic field, velocity using MOC 

// only for 2d
void step_bz_2d(void)
{
  static int nsteps = 0 ;
  
#if((FLOATTYPE==0)&&(SYMFORCEMAG==1))
  fprintf(fail_file,"step_bz: floats and symmetry are known to fail--no idea why yet\n");
  myexit(1);
#endif

  if(nsteps%2 == 0) {
    if(transbzx){
      bzsweepx_2d() ;
      if(ALFVENLIMIT==1) magprepare(); // since bz changed
    }
    if(transbzy){
      bzsweepy_2d() ;
    }
  }
  else {
    if(transbzy){
      bzsweepy_2d() ;
      if(ALFVENLIMIT==1) magprepare(); // since bz changed
    }
    if(transbzx){
      bzsweepx_2d() ;
    }
  }
  
  nsteps++ ;

}


// step x/y magnetic field and velocity using MOC, CT 
//
// CT and Alfven wave source terms


// prim, star, vaa, etc. on corner

void moc_ct_2d(void)
{
  FTYPE fraczone,dqminus,dqplus;
  static FTYPE (*bxstar)[N2M][N1M],(*bystar)[N2M][N1M],(*bzstar)[N2M][N1M],(*vxstar)[N2M][N1M],(*vystar)[N2M][N1M] ,(*vzstar)[N2M][N1M];
  static FTYPE (*dqv)[N3M][N2M][N1M],(*dqb)[N3M][N2M][N1M],(*bxprim)[N2M][N1M],(*byprim)[N2M][N1M],(*bzprim)[N2M][N1M];
  static FTYPE (*magemf)[N2M][N1M],bxa,bya ;
  static FTYPE (*rhoa)[N2M][N1M],(*srhoa)[N2M][N1M];
  static FTYPE rhoa2,srhom,srhop;
  static FTYPE D2,vaa,vxa,vya,vp,vm,bp,bm;
  static FTYPE sgn_va, vaareal ;
  FTYPE dist;
  FTYPE ftemp;
  register int i,j,k ;
  FTYPE ftemp1,ftemp2,ftemp3;
#if(RHOINTERP==0)
  FTYPE Dtota,D1a;
#elif(RHOINTERP==1)
  FTYPE Dtotp,Dtotm,D1m,D1p;
#endif



#if((FLOATTYPE==0)&&(SYMFORCEMAG==1))
  fprintf(fail_file,"moc_ct: floats and symmetry are known to fail--no idea why yet\n");
  myexit(1);
#endif

  dqv = workv1 ;
  dqb = workv2 ;
  magemf = work1 ;
  rhoa=work3;
  srhoa=work4;



  // used twice
  LOOPHC{ // needed since for rhointerp==1, need rho at i=-1..N1, j=-1..N2, so find osqrtrho on LOOPFMHPC!
    srhoa[k][j][i] = z2c_3(osqrtrho,k,j,i) ; //1/sqrt(rho)
  }

  // sweep in x-direction 
  // first get slopes 
  dqvx_calc(2,v[1],dqv) ;
  dqvx_calc(2,v[2],dqb) ;



  vxstar = work5 ;
  vystar = work6 ;
  bxstar = work7 ;
  bystar = work8 ;
  bxprim = work9 ;
  byprim = work10 ;

  LOOPHPC{  // needed since emf needed here, and ok since dq[N1+1] not accessed which isn't created above
    ftemp=z2e_2(v[2][1],k,j,i);
    sgn_va=copysign(1.0,ftemp);
    ftemp=fabs(ftemp)*dt*OARC31(k,j,i);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i] ;
#else
    D1p=ftemp*z2e_2(osqrtrho,k,j,im1mac(i));
    D1m=ftemp*z2e_2(osqrtrho,k,j,i);
#endif

    // determine which formula is for C+ and which is C-
    dist=dx[2][1][i];
    D2 = dt*z2e_2(v[1][1],k,j,i)*OARC31(k,j,i); // fraction of grid

    // Dtot violates symmetry since using fabs, causing mixed vector/scalar, otherwise symmetric like a vector    
    // vp and vm are "located" on corners on full comp grid+border
    // v=vx+|va| // trial Dtot, located on corner
#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    Dtota=D2+D1a;
    if(Dtota > 0) { // C+ from -i region
      vp = v[1][2][k][j][im1] + (1.0 - Dtota )*dqv[2][k][j][im1]*dist ;
      bp = v[2][2][k][j][im1] + (1.0 - Dtota )*dqb[2][k][j][im1]*dist ;
    }
    else{ // C+ is from + region
      vp = v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist ;
      bp = v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist ;
    }

    // v=vx-|va| // trial Dtot
    Dtota=D2-D1a;
    if(Dtota>0) { // C- is from -i region
      vm = v[1][2][k][j][im1] + (1.0 - Dtota )*dqv[2][k][j][im1]*dist ;
      bm = v[2][2][k][j][im1] + (1.0 - Dtota )*dqb[2][k][j][im1]*dist ;
    }
    else{ // C- is from +i region
      vm = v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist ;
    }
#else
    Dtota=D2+D1a;
    if(Dtota > DTOTLIMIT) { // C+ from -i region
      vp = v[1][2][k][j][im1] + (1.0 - Dtota )*dqv[2][k][j][im1]*dist ;
      bp = v[2][2][k][j][im1] + (1.0 - Dtota )*dqb[2][k][j][im1]*dist ;
    }
    else if( Dtota <-DTOTLIMIT) { // C+ is from + region
      vp = v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist ;
      bp = v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist ;
    }
    else{
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][2][k][j][im1] + (1.0 - DTOTLIMIT )*dqv[2][k][j][im1]*dist ;
      dqplus=v[1][2][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[2][k][j][i]*dist ;
      vp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][2][k][j][im1] + (1.0 - DTOTLIMIT )*dqb[2][k][j][im1]*dist ;
      dqplus=v[2][2][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[2][k][j][i]*dist ;
      bp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
    }

    // v=vx-|va| // trial Dtot
    Dtota=D2-D1a;
    if(Dtota > DTOTLIMIT) { // C- is from -i region
      vm = v[1][2][k][j][im1] + (1.0 - Dtota )*dqv[2][k][j][im1]*dist ;
      bm = v[2][2][k][j][im1] + (1.0 - Dtota )*dqb[2][k][j][im1]*dist ;
    }
    else if(Dtota <-DTOTLIMIT) { // C- is from +i region
      vm = v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // Dtotp and Dtotm must have been already evaluated to get here
      vp = 0.5*(v[1][2][k][j][im1] + (1.0 - Dtota )*dqv[2][k][j][im1]*dist+
                v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist) ;
      bp = 0.5*(v[2][2][k][j][im1] + (1.0 - Dtota )*dqb[2][k][j][im1]*dist+
                v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist) ;

#else
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][2][k][j][im1] + (1.0 - DTOTLIMIT )*dqv[2][k][j][im1]*dist ;
      dqplus=v[1][2][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[2][k][j][i]*dist ;
      vm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][2][k][j][im1] + (1.0 - DTOTLIMIT )*dqb[2][k][j][im1]*dist ;
      dqplus=v[2][2][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[2][k][j][i]*dist ;
      bm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }

#endif

#elif(RHOINTERP==1)
    // no strong sym
#if(SYMFORCEMAG==0)
    if( (Dtotp=D2+D1p) > 0) { // C+ from -i region
      vp = v[1][2][k][j][im1] + (1.0 - Dtotp )*dqv[2][k][j][im1]*dist ;
      bp = v[2][2][k][j][im1] + (1.0 - Dtotp )*dqb[2][k][j][im1]*dist ;
    }
    else{ // C+ is from + region
      Dtotm=D2+D1m;
      vp = v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist ;
      bp = v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist ;
    }

    // v=vx-|va| // trial Dtot
    if( (Dtotp=D2-D1p) > 0) { // C- is from -i region
      vm = v[1][2][k][j][im1] + (1.0 - Dtotp )*dqv[2][k][j][im1]*dist ;
      bm = v[2][2][k][j][im1] + (1.0 - Dtotp )*dqb[2][k][j][im1]*dist ;
    }
    else{ // C- is from +i region
      Dtotm=D2-D1m;
      vm = v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist ;
    }

#else

    if( (Dtotp=D2+D1p) > DTOTLIMIT) { // C+ from -i region
      vp = v[1][2][k][j][im1] + (1.0 - Dtotp )*dqv[2][k][j][im1]*dist ;
      bp = v[2][2][k][j][im1] + (1.0 - Dtotp )*dqb[2][k][j][im1]*dist ;
    }
    else if( (Dtotm=D2+D1m) <-DTOTLIMIT) { // C+ is from + region
      vp = v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist ;
      bp = v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist ;
    }
    else{
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // Dtotp and Dtotm must have been already evaluated to get here
      vp = 0.5*(v[1][2][k][j][im1] + (1.0 - Dtotp )*dqv[2][k][j][im1]*dist+
		v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist) ;
      bp = 0.5*(v[2][2][k][j][im1] + (1.0 - Dtotp )*dqb[2][k][j][im1]*dist+
		v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist) ;
    }

    // v=vx-|va| // trial Dtot
    if( (Dtotp=D2-D1p) > DTOTLIMIT) { // C- is from -i region
      vm = v[1][2][k][j][im1] + (1.0 - Dtotp )*dqv[2][k][j][im1]*dist ;
      bm = v[2][2][k][j][im1] + (1.0 - Dtotp )*dqb[2][k][j][im1]*dist ;
    }
    else if( (Dtotm=D2-D1m) <-DTOTLIMIT) { // C- is from +i region
      vm = v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist ;
    }
    else{
      // Dtotp and Dtotm must have been already evaluated to get here
      vm = 0.5*(v[1][2][k][j][im1] + (1.0 - Dtotp )*dqv[2][k][j][im1]*dist+
		v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist) ;
      bm = 0.5*(v[2][2][k][j][im1] + (1.0 - Dtotp )*dqb[2][k][j][im1]*dist+
		v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist);
    }

#endif


#endif // end if rhointerp==1

		
    // should use midpoint of C+/C- for srho
    // solution from characteristics given the footprints from above
#if(RHOINTERP==0)
    // solution to Stone & Norman eqtn 43,44-- use constant rho at C+ C- crossing
    // eq 43,44 should really read: sign(va\dot k)dv=dB
    bystar[k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
    vystar[k][j][i] = 0.5*(vm + vp + sgn_va*(bm - bp)*srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    srhom=z2e_2(osqrtrho,k,j,i); // C- char rho value
    srhop=z2e_2(osqrtrho,k,j,im1mac(i)); // C+ char rho value
    bystar[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
    vystar[k][j][i] = 0.5*(vm+vp+sgn_va*bystar[k][j][i]*(srhop-srhom) +sgn_va*(-bp*srhop+bm*srhom) );
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif

    // find bprim, for updating velocities, no need for vyprime
    // lagrangian frame since op-split and other direction already moving
    // not in source code for more accurate evolution of alphven waves using upwinding
    
#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    // in this case D1p and D1m are equal to D1a
      vp = v[1][2][k][j][im1] + (1.0-D1a)*dqv[2][k][j][im1]*dist ; // C+
      bp = v[2][2][k][j][im1] + (1.0-D1a)*dqb[2][k][j][im1]*dist ; // C+
      vm = v[1][2][k][j][i] + (-1.0+D1a)*dqv[2][k][j][i]*dist ; // C-
      bm = v[2][2][k][j][i] + (-1.0+D1a)*dqb[2][k][j][i]*dist ; // C- 
#else
    // in this case D1p and D1m are equal to D1a
    if(D1a>DTOTLIMIT){ // keep p m names so like above
      vp = v[1][2][k][j][im1] + (1.0-D1a)*dqv[2][k][j][im1]*dist ; // C+
      bp = v[2][2][k][j][im1] + (1.0-D1a)*dqb[2][k][j][im1]*dist ; // C+
      vm = v[1][2][k][j][i] + (-1.0+D1a)*dqv[2][k][j][i]*dist ; // C-
      bm = v[2][2][k][j][i] + (-1.0+D1a)*dqb[2][k][j][i]*dist ; // C- 
    }
    else{
#if(SYMFUNKEDM1==0)
      // true, but only need difference in vp and vm for byprim
      //      vp = vm = 0.5*(v[1][2][k][j][im1] + (1.0-D1a)*dqv[2][k][j][im1]*dist+
      //		     v[1][2][k][j][i] + (-1.0+D1a)*dqv[2][k][j][i]*dist);
      vm = vp = 0.0; // since only need difference to be 0
      bp = bm = 0.5*(v[2][2][k][j][im1] + (1.0-D1a)*dqb[2][k][j][im1]*dist+
                     v[2][2][k][j][i] + (-1.0+D1a)*dqb[2][k][j][i]*dist);

#else
      // this is a little different than sweep.c
      fraczone=D1a*INVDTOTLIMIT;

      dqminus=v[1][2][k][j][im1] + (1.0-DTOTLIMIT)*dqv[2][k][j][im1]*dist ; 
      dqplus=v[1][2][k][j][i] + (-1.0+DTOTLIMIT)*dqv[2][k][j][i]*dist ;
      vm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      vp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));

      //      vm = vp = 0.0; // since only need difference to be 0

      dqminus=v[2][2][k][j][im1] + (1.0-DTOTLIMIT)*dqb[2][k][j][im1]*dist ;
      dqplus=v[2][2][k][j][i] + (-1.0+DTOTLIMIT)*dqb[2][k][j][i]*dist ;
      bm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      bp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
#endif
      //vp = v[1][2][k][j][im1] + (1.0-DTOTLIMIT)*dqv[2][k][j][im1]*dist ; // C+
      //vm = v[1][2][k][j][i] + (-1.0+DTOTLIMIT)*dqv[2][k][j][i]*dist ; // C-
    }
#endif

#elif(RHOINTERP==1)
    // no strong sym
#if(SYMFORCEMAG==0)
      vp = v[1][2][k][j][im1] + (1.0-D1p)*dqv[2][k][j][im1]*dist ; // C+
      bp = v[2][2][k][j][im1] + (1.0-D1p)*dqb[2][k][j][im1]*dist ; // C+
      vm = v[1][2][k][j][i] + (-1.0+D1m)*dqv[2][k][j][i]*dist ; // C-
      bm = v[2][2][k][j][i] + (-1.0+D1m)*dqb[2][k][j][i]*dist ; // C- 
#else
    if(D1p>DTOTLIMIT){
      vp = v[1][2][k][j][im1] + (1.0-D1p)*dqv[2][k][j][im1]*dist ; // C+
      bp = v[2][2][k][j][im1] + (1.0-D1p)*dqb[2][k][j][im1]*dist ; // C+
    }
    else{
      vp = 0.5*(v[1][2][k][j][im1] + (1.0-D1p)*dqv[2][k][j][im1]*dist+
		v[1][2][k][j][i] + (-1.0+D1m)*dqv[2][k][j][i]*dist);
      bp = 0.5*(v[2][2][k][j][im1] + (1.0-D1p)*dqb[2][k][j][im1]*dist+
		v[2][2][k][j][i] + (-1.0+D1m)*dqb[2][k][j][i]*dist);
    }
    if(D1m>DTOTLIMIT){
      vm = v[1][2][k][j][i] + (-1.0+D1m)*dqv[2][k][j][i]*dist ; // C-
      bm = v[2][2][k][j][i] + (-1.0+D1m)*dqb[2][k][j][i]*dist ; // C- 
    }
    else{
      vm = 0.5*(v[1][2][k][j][i] + (-1.0+D1m)*dqv[2][k][j][i]*dist+
		v[1][2][k][j][im1] + (1.0-D1p)*dqv[2][k][j][im1]*dist );
      bm = 0.5*(v[2][2][k][j][i] + (-1.0+D1m)*dqb[2][k][j][i]*dist+
		v[2][2][k][j][im1] + (1.0-D1p)*dqb[2][k][j][im1]*dist);
    }
#endif

#endif // end if rhointerp==1


#if(RHOINTERP==0)
    byprim[k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    // same srhom/srhop from above
    //srhom=z2e_2(osqrtrho,k,j,i); // C- char rho value
    //srhop=z2e_2(osqrtrho,k,j,im1mac(i)); // C+ char rho value
    byprim[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif

  }
  // now vy/bystar live on corners and fill entire grid and border of comp grid


  // sweep in y-direction 
  // first get slopes 
  // dq's act like scalars on rotation
  dqvy_calc(1,v[1],dqv) ;
  dqvy_calc(1,v[2],dqb) ;


  LOOPHPC{ // needed since emf needed here, and ok since dq[N2+1] not accessed which isn't created above
    ftemp=z2e_1(v[2][2],k,j,i);
    sgn_va = copysign(1.,ftemp) ;
    ftemp=fabs(ftemp)*dt*OARC22(k,j,i);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i];
#else
    D1p=ftemp*z2e_1(osqrtrho,k,jm1mac(j),i);
    D1m=ftemp*z2e_1(osqrtrho,k,j,i);
#endif

    // values at the foot of the plus characteristic
    //  are, by convention, in im1 zone
    D2 = dt*z2e_1(v[1][2],k,j,i)*OARC22(k,j,i) ;
    dist=dx[2][2][j]*G2(1,i);


#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    // C+
    Dtota=D2+D1a;
    if(Dtota>0.) { // C+ from -j region
      vp = v[1][1][k][jm1][i] + (1.0 - Dtota )*dqv[1][k][jm1][i]*dist ;
      bp = v[2][1][k][jm1][i] + (1.0 - Dtota )*dqb[1][k][jm1][i]*dist ;
    }
    else{ // C+ from +j region
      vp = v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist ;
      bp = v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist ;
    }
    // C-
    Dtota=D2-D1a;
    if(Dtota>0.) { // C- from -j region
      vm = v[1][1][k][jm1][i] + (1.0 - Dtota )*dqv[1][k][jm1][i]*dist ;
      bm = v[2][1][k][jm1][i] + (1.0 - Dtota )*dqb[1][k][jm1][i]*dist ;
    }
    else{ // C- from +j region
      vm = v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist ;
    }
#else
    // C+
    Dtota=D2+D1a;
    if( Dtota > DTOTLIMIT) { // C+ from -j region
      vp = v[1][1][k][jm1][i] + (1.0 - Dtota )*dqv[1][k][jm1][i]*dist ;
      bp = v[2][1][k][jm1][i] + (1.0 - Dtota )*dqb[1][k][jm1][i]*dist ;
    }
    else if( Dtota < -DTOTLIMIT) { // C+ from +j region
      vp = v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist ;
      bp = v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      vp = 0.5*(v[1][1][k][jm1][i] + (1.0 - Dtota )*dqv[1][k][jm1][i]*dist+
                v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist);
      bp = 0.5*(v[2][1][k][jm1][i] + (1.0 - Dtota )*dqb[1][k][jm1][i]*dist+
                v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist);
#else
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][1][k][jm1][i] + (1.0 - DTOTLIMIT )*dqv[1][k][jm1][i]*dist ;
      dqplus=v[1][1][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[1][k][j][i]*dist ;
      vp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][1][k][jm1][i] + (1.0 - DTOTLIMIT )*dqb[1][k][jm1][i]*dist ;
      dqplus=v[2][1][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[1][k][j][i]*dist ;
      bp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }
    // C-
    Dtota=D2-D1a;
    if(Dtota > DTOTLIMIT) { // C- from -j region
      vm = v[1][1][k][jm1][i] + (1.0 - Dtota )*dqv[1][k][jm1][i]*dist ;
      bm = v[2][1][k][jm1][i] + (1.0 - Dtota )*dqb[1][k][jm1][i]*dist ;
    }
    else if(Dtota < -DTOTLIMIT) { // C- from +j region
      vm = v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      vm = 0.5*(v[1][1][k][jm1][i] + (1.0 - Dtota )*dqv[1][k][jm1][i]*dist+
                v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist);
      bm = 0.5*(v[2][1][k][jm1][i] + (1.0 - Dtota )*dqb[1][k][jm1][i]*dist+
                v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist);
#else
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][1][k][jm1][i] + (1.0 - DTOTLIMIT )*dqv[1][k][jm1][i]*dist ;
      dqplus=v[1][1][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[1][k][j][i]*dist ;
      vm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][1][k][jm1][i] + (1.0 - DTOTLIMIT )*dqb[1][k][jm1][i]*dist ;
      dqplus=v[2][1][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[1][k][j][i]*dist ;
      bm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
#if(SYMFORCEMAG==0)
    // C+
    if( (Dtotp=D2+D1p) > 0.) { // C+ from -j region
      vp = v[1][1][k][jm1][i] + (1.0 - Dtotp )*dqv[1][k][jm1][i]*dist ;
      bp = v[2][1][k][jm1][i] + (1.0 - Dtotp )*dqb[1][k][jm1][i]*dist ;
    }
    else{ // C+ from +j region
      Dtotm=D2+D1m;
      vp = v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist ;
      bp = v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist ;
    }
    // C-
    if( (Dtotp=D2-D1p) > 0.) { // C- from -j region
      vm = v[1][1][k][jm1][i] + (1.0 - Dtotp )*dqv[1][k][jm1][i]*dist ;
      bm = v[2][1][k][jm1][i] + (1.0 - Dtotp )*dqb[1][k][jm1][i]*dist ;
    }
    else{ // C- from +j region
      Dtotm=D2-D1m;
      vm = v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist ;
    }
#else
    // C+
    if( (Dtotp=D2+D1p) > DTOTLIMIT) { // C+ from -j region
      vp = v[1][1][k][jm1][i] + (1.0 - Dtotp )*dqv[1][k][jm1][i]*dist ;
      bp = v[2][1][k][jm1][i] + (1.0 - Dtotp )*dqb[1][k][jm1][i]*dist ;
    }
    else if( (Dtotm=D2+D1m) < -DTOTLIMIT) { // C+ from +j region
      vp = v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist ;
      bp = v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist ;
    }
    else{
      vp = 0.5*(v[1][1][k][jm1][i] + (1.0 - Dtotp )*dqv[1][k][jm1][i]*dist+
		v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist);
      bp = 0.5*(v[2][1][k][jm1][i] + (1.0 - Dtotp )*dqb[1][k][jm1][i]*dist+
		v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist);
    }
    // C-
    if( (Dtotp=D2-D1p) > DTOTLIMIT) { // C- from -j region
      vm = v[1][1][k][jm1][i] + (1.0 - Dtotp )*dqv[1][k][jm1][i]*dist ;
      bm = v[2][1][k][jm1][i] + (1.0 - Dtotp )*dqb[1][k][jm1][i]*dist ;
    }
    else if( (Dtotm=D2-D1m) < -DTOTLIMIT) { // C- from +j region
      vm = v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist ;
    }
    else{
      vm = 0.5*(v[1][1][k][jm1][i] + (1.0 - Dtotp )*dqv[1][k][jm1][i]*dist+
		v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist);
      bm = 0.5*(v[2][1][k][jm1][i] + (1.0 - Dtotp )*dqb[1][k][jm1][i]*dist+
		v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist);
    }
#endif

#endif    // end if rhointerp==1
    // solution to Stone & Norman eqtn 43,44
#if(RHOINTERP==0)
    //-- use constant rho
    bxstar[k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
    vxstar[k][j][i] = 0.5*(vm + vp + sgn_va*(bm - bp)*srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    srhom=z2e_1(osqrtrho,k,j,i); // C- char rho value
    srhop=z2e_1(osqrtrho,k,jm1mac(j),i); // C+ char rho value
    bxstar[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
    vxstar[k][j][i] = 0.5*(vm+vp+sgn_va*bxstar[k][j][i]*(srhop-srhom) +sgn_va*(-bp*srhop+bm*srhom) );
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif

    // find bprim, to update velocities, no need for vxprim
#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    vp = v[1][1][k][jm1][i] + (1.0-D1a)*dqv[1][k][jm1][i]*dist ;
    bp = v[2][1][k][jm1][i] + (1.0-D1a)*dqb[1][k][jm1][i]*dist ;
    vm = v[1][1][k][j][i] + (-1.0+D1a)*dqv[1][k][j][i]*dist ;
    bm = v[2][1][k][j][i] + (-1.0+D1a)*dqb[1][k][j][i]*dist ;
#else
    // D1p=D1m=D1a
    if(D1a>DTOTLIMIT){
      vp = v[1][1][k][jm1][i] + (1.0-D1a)*dqv[1][k][jm1][i]*dist ;
      bp = v[2][1][k][jm1][i] + (1.0-D1a)*dqb[1][k][jm1][i]*dist ;
      vm = v[1][1][k][j][i] + (-1.0+D1a)*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0+D1a)*dqb[1][k][j][i]*dist ;
    }
    else{
#if(SYMFUNKEDM2==0)
      // true, but only need diff (vm-vp) for bxprim
      //      vp = vm = 0.5*(v[1][1][k][jm1][i] + (1.0-D1a)*dqv[1][k][jm1][i]*dist+
      //                v[1][1][k][j][i] + (-1.0+D1a)*dqv[1][k][j][i]*dist);
      vp=vm=0.0;
      bp = bm = 0.5*(v[2][1][k][jm1][i] + (1.0-D1a)*dqb[1][k][jm1][i]*dist+
                     v[2][1][k][j][i] + (-1.0+D1a)*dqb[1][k][j][i]*dist);
#else
      // this is a little different than sweep.c
      fraczone=D1a*INVDTOTLIMIT;

      dqminus=v[1][1][k][jm1][i] + (1.0-DTOTLIMIT)*dqv[1][k][jm1][i]*dist ;
      dqplus=v[1][1][k][j][i] + (-1.0+DTOTLIMIT)*dqv[1][k][j][i]*dist ;
      vm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      vp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));

      dqminus=v[2][1][k][jm1][i] + (1.0-DTOTLIMIT)*dqb[1][k][jm1][i]*dist ;
      dqplus=v[2][1][k][j][i] + (-1.0+DTOTLIMIT)*dqb[1][k][j][i]*dist ;
      bm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      bp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
    // no strong sym
#if(SYMFORCEMAG==0)
      vp = v[1][1][k][jm1][i] + (1.0-D1p)*dqv[1][k][jm1][i]*dist ;
      bp = v[2][1][k][jm1][i] + (1.0-D1p)*dqb[1][k][jm1][i]*dist ;
      vm = v[1][1][k][j][i] + (-1.0+D1m)*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0+D1m)*dqb[1][k][j][i]*dist ;
#else
    if(D1p>DTOTLIMIT){
      vp = v[1][1][k][jm1][i] + (1.0-D1p)*dqv[1][k][jm1][i]*dist ;
      bp = v[2][1][k][jm1][i] + (1.0-D1p)*dqb[1][k][jm1][i]*dist ;
    }
    else{
      vp = 0.5*(v[1][1][k][jm1][i] + (1.0-D1p)*dqv[1][k][jm1][i]*dist+
		v[1][1][k][j][i] + (-1.0+D1m)*dqv[1][k][j][i]*dist);
      bp = 0.5*(v[2][1][k][jm1][i] + (1.0-D1p)*dqb[1][k][jm1][i]*dist+
		v[2][1][k][j][i] + (-1.0+D1m)*dqb[1][k][j][i]*dist);
    }
    if(D1m>DTOTLIMIT){
      vm = v[1][1][k][j][i] + (-1.0+D1m)*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0+D1m)*dqb[1][k][j][i]*dist ;
    }
    else{
      vm = 0.5*(v[1][1][k][jm1][i] + (1.0-D1p)*dqv[1][k][jm1][i]*dist+
		v[1][1][k][j][i] + (-1.0+D1m)*dqv[1][k][j][i]*dist);
      bm = 0.5*(v[2][1][k][jm1][i] + (1.0-D1p)*dqb[1][k][jm1][i]*dist+
		v[2][1][k][j][i] + (-1.0+D1m)*dqb[1][k][j][i]*dist);      
    }
#endif

#endif // end if rhointerp==1

#if(RHOINTERP==0)
    bxprim[k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    //srhom=z2e_1(osqrtrho,k,j,i); // C- char rho value
    //  srhop=z2e_1(osqrtrho,k,jm1mac(j),i); // C+ char rho value
    bxprim[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif

  }
  // now vx/bxstar live on corners on full grid, act like vectors on rotation


  // calculate emf=v cross B, e3 component is only that survives


  // emf lives on corners on full comp grid including border
  // emf is 3-vector, so acts like scalar on rotations
  if(RESMEM&&(res_real==1)){
    LOOPHPC{ // needed for below 2 loops.
      magemf[k][j][i] = vxstar[k][j][i]*bystar[k][j][i] - vystar[k][j][i]*bxstar[k][j][i] - jcurrent[3][k][j][i]*z2c_3(nu_res_real,k,j,i) ;
    }    
  }
  else{
    LOOPHPC{ // needed for below 2 loops.
      magemf[k][j][i] = vxstar[k][j][i]*bystar[k][j][i] - vystar[k][j][i]*bxstar[k][j][i] ;
    }
  }
  
  // STONE AND MILLER
  // update field, CT with emf as conserved flux, hence divB=0.  eq 11,12 for B, eq 30,31 for v
  if(mocctvx1){
    LOOPV1{ // needs bxprim on i=0..N1-1 and j=0..N2
#if(ALFVENLIMIT==0)
      rhoa2 = z2e_1(s[1],k,j,i);
#else
      rhoa2 = z2e_1(rholimited[0],k,j,i);
#endif
      bya = v2tov1(v[2][2],k,j,i);
      v[1][1][k][j][i] += dt*bya*OARC42(k,j,i)/rhoa2*(bxprim[k][jp1][i] - bxprim[k][j][i]);
    }
  }
  // update velocity(2nd partial update of Lorentz force terms resulting in Alfven wave motion) (comoving)
  if(mocctvx2){
    LOOPV2{ // needs byprim at i=0..N1 and j=0..N2-1+periodicx2special (LOOPH above is fine)
#if(ALFVENLIMIT==0)
      rhoa2 =z2e_2(s[1],k,j,i);
#else
      rhoa2 =z2e_2(rholimited[0],k,j,i);
#endif
      bxa = v1tov2(v[2][1],k,j,i);
      v[1][2][k][j][i] += dt*bxa*OARC41(k,j,i)/(rhoa2*G2(2,i))*(G2(1,ip1mac(i))*byprim[k][j][ip1] - G2(1,i)*byprim[k][j][i]) ;
    }
  }
  // must do mag field second since above v's use old B (like new Zeus says.  Old zeus says do field first)
  if(mocctbx1){
    LOOPV1{ // needs emf on i=0..N1-1 and j=0..N2
      v[2][1][k][j][i] += dt*OARC42(k,j,i)/G4(2,j)*(G4(1,jp1mac(j))*magemf[k][jp1][i] - G4(1,j)*magemf[k][j][i]) ;
    }
  }
  if(mocctbx2){
    LOOPV2{ // needs emf at i=0..N1 and j=0..N2-1+periodicx2special (LOOPH above is fine)
      v[2][2][k][j][i] += dt*OARC41(k,j,i)/G3(2,i)*(G3(1,i)*magemf[k][j][i] - G3(1,ip1mac(i))*magemf[k][j][ip1]) ;
      
    }
  }

  bound(NULL,NULL,0,-1,12); // only need to bound 1 and 2 components

}



