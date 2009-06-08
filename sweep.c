#define SYMFORCEHD 0 // whether to force symmetry on dq in HD(sweepx/y.c) stuff


//Checked by Jon

#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif


/*
	Notes:

	Rate of convergence for cart is at best == 2 for nearly constant velocity

	Rate of convergence for spc is at best == 1.5 for nearly constant velocity.  Why not 2?  Both have ROC of 1 for general velocity
	Perhaps could do double velocity loop in advection to see if can get better ROC when vel not constant

	The floor may be implemented better somehow

*/

/*
	Notes on what needs what:
	x1-dir:
	mdot= [s2-1,N2-1],[s1-1,N1]
	en/rho/vz/Bz= fl=[0,N2-1],[0,N1] , mdot=fl, dq=[fl][-1,N1]
	vx = fl=[0,N2-1],[s1-1,N1-1], mdot=[0,N2-1],[s1-1,N1], dq=[fl][-1,N1]
	vy = fl=[s2,N2-1],[0,N1], mdot=[s2-1,N2-1],[0,N1], dq=[fl][-1,N1]

	x2-dir:
	mdot= [s2-1,N2],[s1-1,N1-1]
	en/rho/vz/Bz= fl=[0,N2],[0,N1-1], mdot=fl, dq=[-1,N2][fl]
	vy = fl=[s2-1,N2-1],[0,N1-1], mdot=[s2-1,N2],[0,N1-1], dq=[-1,N2][fl]
	vx = fl=[0,N2-1][s1,N1-1], mdot=[0,N2],[s1-1,N1-1], dq=[-1,N2][fl]

	So when s1=0, need to at least bound mdot1, when s2=0, need to at least bound mdot2 due to dq[fl][-2,N1] needed for mdot1 and dq[-2,N2][fl] needed for mdot2

	This means when periodic in a direction, bound that direction, and when MPI in that direction, must bound in that direction.  Must use dual if(  if((skipix2==0)||(numprocs>1))) because, say, 0 processor needs to transfer +j info to 1, even though skipix2=0 for 1.  May desire to tell 0 not to do anything but transfer(i.e. no regular bound on mdot)....same for numprocs-1.

*/


//   Do transport step.  Order of sweep direction is
//   varied from step to step. 
void step_trans_2d(void)
{
  int i,j,k;
  static int nstep = 0 ;
 
#if((FLOATTYPE==0)&&(SYMFORCEHD==1))
  fprintf(fail_file,"floats and symmetry are known to fail--no idea why yet\n");
  myexit(1);
#endif

 
  if(nstep%2 == 0) {
    if(transx1) sweepx() ;
    if(transx2) sweepy() ;
  }
  else {
    if(transx2) sweepy() ;
    if(transx1) sweepx() ;
  }
  nstep++ ;
}



//   Do transport step.  Order of sweep direction is
//   varied from step to step. 
void step_trans_3d(void)
{// transx? is set to 0 in init.c if Nx=1
  int i,j,k;
  static int nstep = 0 ;
 
#if((FLOATTYPE==0)&&(SYMFORCEHD==1))
  fprintf(fail_file,"floats and symmetry are known to fail--no idea why yet\n");
  myexit(1);
#endif
 
  if(nstep%6 == 0) {
    if(transx1) sweepx() ;
    if(transx2) sweepy() ;
    if(transx3) sweepz() ;
  }
  else if(nstep%6 == 1) {
    if(transx1) sweepx() ;
    if(transx3) sweepz() ;
    if(transx2) sweepy() ;
  }
  else if(nstep%6 == 2) {
    if(transx2) sweepy() ;
    if(transx1) sweepx() ;
    if(transx3) sweepz() ;
  }
  else if(nstep%6 == 3) {
    if(transx2) sweepy() ;
    if(transx3) sweepz() ;
    if(transx1) sweepx() ;
  }
  else if(nstep%6 == 4) {
    if(transx3) sweepz() ;
    if(transx1) sweepx() ;
    if(transx2) sweepy() ;
  }
  else if(nstep%6 == 5) {
    if(transx3) sweepz() ;
    if(transx2) sweepy() ;
    if(transx1) sweepx() ;
  }

  
  //  LOOPFC{
  // fprintf(stderr,"bob %d %d %d val=%g\n",i,j,k,v[2][1][k][j][i]); fflush(stderr);
  //}
  


  nstep++ ;
}


void sweepx(void)
{
  // BEGIN variables
  FTYPE vpp,mdotp ;
  FTYPE dqminus,dqplus,fraczone;
  static FTYPE (*dq)[N2M][N1M],(*dqv)[N3M][N2M][N1M],(*vp)[N3M][N2M][N1M];
  static FTYPE (*mdot)[N3M][N2M][N1M] ;
  static FTYPE (*p)[N3M][N2M][N1M];
  static FTYPE (*u)[N2M][N1M],(*Bzr)[N2M][N1M];
  static FTYPE (*uer)[N2M][N1M];
  static FTYPE (*workv)[N3M][N2M][N1M];
  static FTYPE (*fl)[N3M][N2M][N1M];
  //static SFTYPE mdot[3+1][N3M][N2M][N1M];
  //static SFTYPE fl[3+1][N3M][N2M][N1M];
  //static SFTYPE scatemp[3+1][N3M][N2M][N1M];
  
  FTYPE Dtot;
  register int i;
  int j,k ;
  FTYPE ftemp;
  FTYPE ftemp0,ftemp1,ftemp2;
  FTYPE ftempv;
  FTYPE v1xa,v1ya,v1za,momold,momnew;
  FILE * out;

  static FTYPE sumix1=0,sumox1=0;
  int wcom;
  // END variables


#if(DO_ACCRETOR)
  bound_accretor(1,NULL,NULL,-1,-1,123);
#endif



  // pointer assignments
  p    = workv1 ;
  vp   = workv2 ;
  fl   = workv3 ;
  mdot = workv4 ;
  dqv  = workv5 ;

  dq  = work1 ;
  u   = work2 ;
  uer = work4 ;
  Bzr = work3 ;

  /* transform to momenta */
  LOOPC{ // no need for LOOPV here since multiplying and that's ok. 
    p[1][k][j][i] = z2e_1(s[1],k,j,i)*v[1][1][k][j][i] ;
    p[2][k][j][i] = z2e_2(s[1],k,j,i)*v[1][2][k][j][i]*G2(2,i) ;
    p[3][k][j][i] = z2e_3(s[1],k,j,i)*v[1][3][k][j][i]*G3(2,i)*G4(2,j) ;
  } // no need to bound

  LOOPHC{ // only need half-full loop for sweep, no need to bound
    /* create spacial tranport variable */
    vp[1][k][j][i] = (v[1][1][k][j][i]-vg[1])*dt;
  }

  
  /* first do mass: rho */
  dqx_calc(s[1],dq) ;

  // mdot is really mdot=Mdot*dt
  LOOPT0i{
#if(SYMFORCEHD==0)
    if(vp[1][k][j][i] > 0.){
      mdot[1][k][j][i] = (s[1][k][j][im1] + (dx[2][1][i] - vp[1][k][j][i])*dq[k][j][im1])*vp[1][k][j][i];
    }
    else{
      mdot[1][k][j][i] = (s[1][k][j][i] + (-dx[2][1][i] - vp[1][k][j][i])*dq[k][j][i])*vp[1][k][j][i];
    }
#else
    Dtot=vp[1][k][j][i]*OARC11(k,j,i);
    if(Dtot > DTOTLIMIT){
      mdot[1][k][j][i] = (s[1][k][j][im1] + (1.0 - Dtot)*dq[k][j][im1]*dx[2][1][i])*vp[1][k][j][i];
    }
    else if(Dtot < -DTOTLIMIT) {
      mdot[1][k][j][i] = (s[1][k][j][i] + (-1.0 - Dtot)*dq[k][j][i]*dx[2][1][i])*vp[1][k][j][i];
    }
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      dqminus=(s[1][k][j][im1] + (1.0 - DTOTLIMIT)*dq[k][j][im1]*dx[2][1][i])*vp[1][k][j][i];
      dqplus=(s[1][k][j][i] + (-1.0 + DTOTLIMIT)*dq[k][j][i]*dx[2][1][i])*vp[1][k][j][i];
      mdot[1][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
#endif
  }
  // only need to bound if doing periodic boundary conditions.  This is identically when you don't skip the first zone

  if((skipix1==0)||(ncpux1>1)){
    bound(NULL,mdot,0,-3,1);
  }

#include "sweeppassive1.h"

  /*
  LOOPFC{
    fprintf(stderr,"bob %d %d %d val=%g\n",i,j,k,v[2][1][k][j][i]); fflush(stderr);
  }
  */


  if(transiex1){
    if(wgam){
      /* then specific internal energy */
      LOOPFC{ // full due to dqx_calc requirements
				u[k][j][i] = s[2][k][j][i]/s[1][k][j][i] ;
      }
      
      dqx_calc(u,dq) ;

      // fl is really fl=Flux*dt/dx (because of mdot)
      LOOPT1i{ // same as density w.r.t. loop/bound
#if(SYMFORCEHD==0)
				if(vp[1][k][j][i] > 0.) 
					fl[1][k][j][i] = (u[k][j][im1] + (dx[2][1][i] - vp[1][k][j][i])*dq[k][j][im1])*mdot[1][k][j][i]*DS11(k,j,i); // dx[1][3][k] cancels with below volume
				else	
					fl[1][k][j][i] = (u[k][j][i] + (-dx[2][1][i] - vp[1][k][j][i])*dq[k][j][i])*mdot[1][k][j][i]*DS11(k,j,i);
#else
				Dtot=vp[1][k][j][i]*OARC11(k,j,i);
				if(Dtot > DTOTLIMIT )
					fl[1][k][j][i] = (u[k][j][im1] + (1.0 - Dtot)*dq[k][j][im1]*dx[2][1][i])*mdot[1][k][j][i]*DS11(k,j,i);
				else if(Dtot < -DTOTLIMIT)
					fl[1][k][j][i] = (u[k][j][i] + (-1.0 - Dtot)*dq[k][j][i]*dx[2][1][i])*mdot[1][k][j][i]*DS11(k,j,i);
				else{
					fraczone=Dtot*INVDTOTLIMIT;
					ftemp = mdot[1][k][j][i]*DS11(k,j,i);
					dqminus=(u[k][j][im1] + (1.0 - DTOTLIMIT)*dq[k][j][im1]*dx[2][1][i])*ftemp;
					dqplus=(u[k][j][i] + (-1.0 + DTOTLIMIT)*dq[k][j][i]*dx[2][1][i])*ftemp;
					fl[1][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
				}
#endif
      }      
      // assign internal energy density from surface fluxes
      LOOPC{
				s[2][k][j][i] += (fl[1][k][j][i]-fl[1][k][j][ip1])*OVOL1(k,j,i); // *odx[1][3][k]
      }
      if(FORCEIEINTERNAL){
				floor_correct(2,4);
      }
      bound(NULL,NULL,2,0,0) ;
    }
  }
  else{
    if(DOLOSSDIAG){
      LOOPT1i{
				fl[1][k][j][i]=z2e_1(s[2],k,j,i)/z2e_1(s[1],k,j,i)*mdot[1][k][j][i]*DS11(k,j,i); // division by interp may be problem
      }
    }
  }
  
  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(1,2,fl); // account for "missing" dx[?][3][k] here
  }
  /*
  // GODMARK -- commented
  fprintf(log_file,"%ld: after sweepx(en) : symmetry check\n",nstep); fflush(stdout);
  symmetry_check(0); // before bound symmetry check// GODMARK
  //image(888,-1,-1,0,0);
  //myexit(0);
  */


#if(RAD)
      /* then specific internal energy */
      LOOPFC{ // full due to dqx_calc requirements
                                uer[k][j][i] = er[k][j][i]/s[1][k][j][i] ;
      }

      dqx_calc(uer,dq) ;

      // fl is really fl=Flux*dt/dx (because of mdot)
      LOOPT1i{ // same as density w.r.t. loop/bound
#if(SYMFORCEHD==0)
                                if(vp[1][k][j][i] > 0.)
                                        fl[1][k][j][i] = (uer[k][j][im1] + (dx[2][1][i] - vp[1][k][j][i])*dq[k][j][im1])*mdot[1][k][j][i]*DS11(k,j,i); // dx[1][3][k] cancels with below volume
                                else
                                        fl[1][k][j][i] = (uer[k][j][i] + (-dx[2][1][i] - vp[1][k][j][i])*dq[k][j][i])*mdot[1][k][j][i]*DS11(k,j,i);
#else
                                Dtot=vp[1][k][j][i]*OARC11(k,j,i);
                                if(Dtot > DTOTLIMIT )
                                        fl[1][k][j][i] = (uer[k][j][im1] + (1.0 - Dtot)*dq[k][j][im1]*dx[2][1][i])*mdot[1][k][j][i]*DS11(k,j,i);
                                else if(Dtot < -DTOTLIMIT)
                                        fl[1][k][j][i] = (uer[k][j][i] + (-1.0 - Dtot)*dq[k][j][i]*dx[2][1][i])*mdot[1][k][j][i]*DS11(k,j,i);
                                else{
                                        fraczone=Dtot*INVDTOTLIMIT;
                                        ftemp = mdot[1][k][j][i]*DS11(k,j,i);
                                        dqminus=(uer[k][j][im1] + (1.0 - DTOTLIMIT)*dq[k][j][im1]*dx[2][1][i])*ftemp;
                                        dqplus=(uer[k][j][i] + (-1.0 + DTOTLIMIT)*dq[k][j][i]*dx[2][1][i])*ftemp;
                                        fl[1][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
                                }
#endif
      }
      // assign internal energy density from surface fluxes
      LOOPC{
                                er[k][j][i] += (fl[1][k][j][i]-fl[1][k][j][ip1])*OVOL1(k,j,i); // *odx[1][3][k]
      }
      if(FORCEIEINTERNAL){
                                floor_correct(2,4);
      }
      bvaler();
#endif








  if(COMPDIM==2){
    if(transmagx1&&mag){ // no need to change in 3D since not used
      // magnetic field: advect flux, not field, eq 49, 50.
      LOOPFC{// need full loop due to dq calc
				Bzr[k][j][i] = v[2][3][k][j][i]/s[1][k][j][i] ;
      }
      
      dqx_calc(Bzr,dq); // same as dqvx wcom==3
      
      LOOPT1i{ // same as density w.r.t. loop/bound
#if(SYMFORCEHD==0)
				if(vp[1][k][j][i] > 0) 
					fl[1][k][j][i] = (Bzr[k][j][im1] + (dx[2][1][i] - vp[1][k][j][i])*dq[k][j][im1])*mdot[1][k][j][i]*G2(1,i);
				else	
					fl[1][k][j][i] = (Bzr[k][j][i] + (-dx[2][1][i] - vp[1][k][j][i])*dq[k][j][i])*mdot[1][k][j][i]*G2(1,i);
#else
				Dtot=vp[1][k][j][i]*OARC11(k,j,i);
				if(Dtot > DTOTLIMIT)
					fl[1][k][j][i] = (Bzr[k][j][im1] + (1.0 - Dtot)*dq[k][j][im1]*dx[2][1][i])*mdot[1][k][j][i]*G2(1,i);
				else if(Dtot< -DTOTLIMIT)
					fl[1][k][j][i] = (Bzr[k][j][i] + (-1.0 - Dtot)*dq[k][j][i]*dx[2][1][i])*mdot[1][k][j][i]*G2(1,i);
				else{
					fraczone=Dtot*INVDTOTLIMIT;
					ftemp = mdot[1][k][j][i]*G2(1,i);
					dqminus=(Bzr[k][j][im1] + (1.0 - DTOTLIMIT)*dq[k][j][im1]*dx[2][1][i])*ftemp;
					dqplus=(Bzr[k][j][i] + (-1.0 + DTOTLIMIT)*dq[k][j][i]*dx[2][1][i])*ftemp;
					fl[1][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
				}
#endif
      }
      
      LOOPC{
				v[2][3][k][j][i] += (fl[1][k][j][i]-fl[1][k][j][ip1])*OARC21(k,j,i)/G2(2,i);
      }
      bound(NULL,NULL,0,2,3) ; // bound as vector, only needed 3-comp bounded 
      

      // put here until simulate flux like rest
      if(transmagx1&&mag&&DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
				hydro_flux_adv(1,3,fl);
      }
    }
  }


  if(transv1x1){

    /* vx */
    dqvx_calc(1,v[1],dqv) ;

    
    LOOPT2i{ // only need to loop over -1 to N-1(1 bnd zone) and should not bound fl!
      vpp = e2z_1(vp[1],k,j,i);
      mdotp = e2z_1(mdot[1],k,j,i);
#if(SYMFORCEHD==0)
      if(vpp > 0)
				fl[1][k][j][i] = (v[1][1][k][j][i]   + (dx[1][1][i] - vpp)*dqv[1][k][j][i])  *mdotp*DS21(k,j,i);
      else
				fl[1][k][j][i] = (v[1][1][k][j][ip1] + (-dx[1][1][i] - vpp)*dqv[1][k][j][ip1])*mdotp*DS21(k,j,i);
#else
      Dtot=vpp*OARC21(k,j,i);
      if(Dtot > DTOTLIMIT) 
				fl[1][k][j][i] = (v[1][1][k][j][i]   + (1.0 - Dtot)*dqv[1][k][j][i]*dx[1][1][i])*mdotp*DS21(k,j,i);
      else if(Dtot< -DTOTLIMIT)
				fl[1][k][j][i] = (v[1][1][k][j][ip1] + (-1.0 - Dtot)*dqv[1][k][j][ip1]*dx[1][1][i])*mdotp*DS21(k,j,i);
      else{
				fraczone=Dtot*INVDTOTLIMIT;
				ftemp = mdotp*DS21(k,j,i);
				dqminus=(v[1][1][k][j][i]   + (1.0 - DTOTLIMIT)*dqv[1][k][j][i]*dx[1][1][i])*ftemp;
				dqplus=(v[1][1][k][j][ip1] + (-1.0 + DTOTLIMIT)*dqv[1][k][j][ip1]*dx[1][1][i])*ftemp;
				fl[1][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
      }
#endif    
    }
    
    LOOPV1{
      p[1][k][j][i] += (fl[1][k][j][im1]-fl[1][k][j][i])*OVOL2(k,j,i);
    }
  }
  else{
    if(DOLOSSDIAG){
      LOOPT2i{
				fl[1][k][j][i] = v[1][1][k][j][i]*mdot[1][k][j][i]*DS11(k,j,i);
      }
    }
  }

  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(1,5,fl);
  }


  if(transv2x1){
    /* vy */
    dqvx_calc(2,v[1],dqv) ;
    
    LOOPT3i{ // flux here lives on lower left corner
      mdotp = z2e_2(mdot[1],k,j,i);
      vpp = z2e_2(vp[1],k,j,i);
#if(SYMFORCEHD==0)     
      if(vpp > 0.)
				fl[1][k][j][i] = (v[1][2][k][j][im1] + (dx[2][1][i] - vpp)*dqv[2][k][j][im1])*G2(1,i)*mdotp*DS31(k,j,i);
      else	
				fl[1][k][j][i] = (v[1][2][k][j][i]   + (-dx[2][1][i] - vpp)*dqv[2][k][j][i])  *G2(1,i)*mdotp*DS31(k,j,i);
#else
      Dtot=vpp*OARC11(k,j,i);
      if(Dtot > DTOTLIMIT)
				fl[1][k][j][i] = (v[1][2][k][j][im1] + (1.0 - Dtot)*dqv[2][k][j][im1]*dx[2][1][i])*G2(1,i)*mdotp*DS31(k,j,i);
      else if(Dtot< -DTOTLIMIT)
				fl[1][k][j][i] = (v[1][2][k][j][i]   + (-1.0 - Dtot)*dqv[2][k][j][i]*dx[2][1][i])*G2(1,i)*mdotp*DS31(k,j,i);
      else{
				fraczone=Dtot*INVDTOTLIMIT;
				ftemp = G2(1,i)*mdotp*DS31(k,j,i);
				dqminus= (v[1][2][k][j][im1] + (1.0 - DTOTLIMIT)*dqv[2][k][j][im1]*dx[2][1][i])*ftemp;
				dqplus=(v[1][2][k][j][i]   + (-1.0 + DTOTLIMIT)*dqv[2][k][j][i]*dx[2][1][i])*ftemp;
				fl[1][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
      }
#endif
    }
    
    LOOPV2{
      p[2][k][j][i] += (fl[1][k][j][i]-fl[1][k][j][ip1])*OVOL3(k,j,i);
    }
  }
  else{
    if(DOLOSSDIAG){
      LOOPT3i{
				fl[1][k][j][i] = v2tov1(v[1][2],k,j,i)*G2(1,i)*mdot[1][k][j][i]*DS11(k,j,i);
      }
    }
  }


  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(1,6,fl);
  }

  if(transv3x1){
    
    /* vz */
    dqvx_calc(3,v[1],dqv) ;
    
    LOOPT4i{ //  dx[?][3][k] again cancels
      mdotp = z2e_3(mdot[1],k,j,i);
      vpp = z2e_3(vp[1],k,j,i);
#if(SYMFORCEHD==0)
      if(vpp > 0) 
				fl[1][k][j][i] = (v[1][3][k][j][im1] + (dx[2][1][i] - vpp)*dqv[3][k][j][im1])*G3(1,i)*G4(2,j)*mdotp*DS11(k,j,i);
      else
				fl[1][k][j][i] = (v[1][3][k][j][i] + (-dx[2][1][i] - vpp)*dqv[3][k][j][i])*G3(1,i)*G4(2,j)*mdotp*DS11(k,j,i);
#else
      Dtot = vpp*OARC11(k,j,i);
      if(Dtot > DTOTLIMIT) 
				fl[1][k][j][i] = (v[1][3][k][j][im1] + (1.0 - Dtot)*dqv[3][k][j][im1]*dx[2][1][i])*G3(1,i)*G4(2,j)*mdotp*DS11(k,j,i);
      else if(Dtot< -DTOTLIMIT)
				fl[1][k][j][i] = (v[1][3][k][j][i] + (-1.0 - Dtot)*dqv[3][k][j][i]*dx[2][1][i])*G3(1,i)*G4(2,j)*mdotp*DS11(k,j,i);
      else{
				fraczone=Dtot*INVDTOTLIMIT;
				ftemp =G3(1,i)*G4(2,j)*mdotp*DS11(k,j,i);
				dqminus=(v[1][3][k][j][im1] + (1.0 - DTOTLIMIT)*dqv[3][k][j][im1]*dx[2][1][i])*ftemp;
				dqplus=(v[1][3][k][j][i] + (-1.0 + DTOTLIMIT)*dqv[3][k][j][i]*dx[2][1][i])*ftemp;
				fl[1][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
      }
			      
#endif
    }
    
    LOOPV3{
      p[3][k][j][i] += (fl[1][k][j][i]-fl[1][k][j][ip1])*OVOL1(k,j,i);
    }// do not bound!
  }
  else{
    if(DOLOSSDIAG){
      LOOPT1i{
				fl[1][k][j][i] =v3tov1(v[1][3],k,j,i)*G3(1,i)*G4(2,j)*mdot[1][k][j][i]*DS11(k,j,i);
      }
    }
  }


  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(1,7,fl);
  }


  if(transrhox1){
    // fl is really fl=Flux*dt
    LOOPT1i{ // loops so gets only needed flux values for density calculation below, no BC application should be applied
      fl[1][k][j][i] =  mdot[1][k][j][i]*DS11(k,j,i);
    }
    
    // correction and calculation of rho
    
    // compute local density from surface fluxes
    LOOPC{
      s[1][k][j][i] += (fl[1][k][j][i] - fl[1][k][j][ip1])*OVOL1(k,j,i);
    }     
    if(FORCERHOINTERNAL){
      floor_correct(1,4);
    }
    bound(NULL,NULL,1,0,0) ;
  }
  else{
    if(DOLOSSDIAG){
      LOOPT1i{
				fl[1][k][j][i] =  mdot[1][k][j][i]*DS11(k,j,i);
      }
    }
  }


  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(1,1,fl);
  }


  /* return to velocities */
  
  if(transv1x1){
    LOOPV1{
      v[1][1][k][j][i] = p[1][k][j][i]/z2e_1(s[1],k,j,i); // division by interp may be a problem

    }
  }
  if(transv2x1){
    LOOPV2{
      v[1][2][k][j][i] = p[2][k][j][i]/(z2e_2(s[1],k,j,i)*G2(2,i)) ; // division by interp may be a problem
    }
  }
  if(transv3x1){
    LOOPV3{
      v[1][3][k][j][i] = p[3][k][j][i]/(z2e_3(s[1],k,j,i)*G3(2,i)*G4(2,j));
    }
  }
  if( (transv1x1)||(transv2x1)||(transv3x1)){
    bound(NULL,NULL,0,1,123);  // bound all grid velocity-vector components
  }

}


void sweepy(void)
{
  FTYPE vpp,mdotp ;
  static FTYPE (*dq)[N2M][N1M],(*dqv)[N3M][N2M][N1M],(*vp)[N3M][N2M][N1M] ;
  static FTYPE (*mdot)[N3M][N2M][N1M] ;
  static FTYPE (*p)[N3M][N2M][N1M];
  static FTYPE (*u)[N2M][N1M],(*Bzr)[N2M][N1M];
  static FTYPE (*uer)[N2M][N1M];
  static FTYPE (*workv)[N3M][N2M][N1M];
  static FTYPE (*fl)[N3M][N2M][N1M];
  register int i,j,k ;
  FTYPE Dtot;
  FTYPE dqminus,dqplus,fraczone;

  FTYPE v1xa,v1ya,v1za,momold,momnew,ftemp;
  FTYPE ftempv;
  FTYPE ftemp0,ftemp1,ftemp2;
  FILE * out;
  int wcom;
  /* transform to momenta */


#if(DO_ACCRETOR)
  bound_accretor(2,NULL,NULL,-1,-1,123);
#endif


  // pointer assignments
  p    = workv1 ;
  vp   = workv2 ;
  fl   = workv3 ;
  mdot = workv4 ;
  dqv  = workv5 ;

  dq  = work1 ;
  u   = work2 ;
  uer = work4 ;
  Bzr = work3 ;

  LOOPC{ // no need for LOOPV here since multiplying and that's ok. 
    p[1][k][j][i] = z2e_1(s[1],k,j,i)*v[1][1][k][j][i] ;
    p[2][k][j][i] = z2e_2(s[1],k,j,i)*v[1][2][k][j][i]*G2(2,i) ;
    p[3][k][j][i] = z2e_3(s[1],k,j,i)*v[1][3][k][j][i]*G3(2,i)*G4(2,j) ;
  } // no need to bound


  LOOPHC{ // only need half full loop, no need to bound
    /* create spacial tranport variable */
    vp[2][k][j][i] = (v[1][2][k][j][i]-vg[2])*dt;
  }



  /* first do mass: rho */
  dqy_calc(s[1],dq) ;


  // mdot/fl is really fl=Flux*dt and mdot=Mdot*dt
  LOOPT0j{
#if(SYMFORCEHD==0)
    if(vp[2][k][j][i] > 0.){
      mdot[2][k][j][i] = (s[1][k][jm1][i] + (dx[2][2][j]*G2(2,i) - vp[2][k][j][i])*dq[k][jm1][i])*vp[2][k][j][i]; // if periodic, uses rho(-3), hence bound below
    }
    else{
      mdot[2][k][j][i] = (s[1][k][j][i] + (-dx[2][2][j]*G2(2,i) - vp[2][k][j][i])*dq[k][j][i])*vp[2][k][j][i]; // dq here needs s[1][k][j=N2+2][i] if periodicx2special, but bound takes care of that since was needed anyways for inner edge if periodic (see global.h on LOOPT0j for more notes)
    }
#else
    Dtot=vp[2][k][j][i]*OARC12(k,j,i);
    if(Dtot > DTOTLIMIT){
      mdot[2][k][j][i] = (s[1][k][jm1][i] + (1.0 - Dtot)*dq[k][jm1][i]*dx[2][2][j]*G2(2,i))*vp[2][k][j][i];
    }
    else if(Dtot< -DTOTLIMIT) {
      mdot[2][k][j][i] = (s[1][k][j][i] + (-1.0 - Dtot)*dq[k][j][i]*dx[2][2][j]*G2(2,i))*vp[2][k][j][i];
    }
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      ftemp=vp[2][k][j][i];
      ftemp2=dx[2][2][j]*G2(2,i);
      dqminus=(s[1][k][jm1][i] + (1.0 - DTOTLIMIT)*dq[k][jm1][i]*ftemp2)*ftemp;
      dqplus=(s[1][k][j][i] + (-1.0 + DTOTLIMIT)*dq[k][j][i]*ftemp2)*ftemp;
      mdot[2][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
#endif
  }


  if((skipix2==0)||(ncpux2>1)){
    bound(NULL,mdot,0,-3,2); // since needed corner zone info on rho beyond 2 zones, and for periodic need rho(-3) and for periodicx2special need rho(N2+3)
  }



#include "sweeppassive2.h"



  if(transiex2){
    if(wgam){
      
      /* then specific internal energy */
      LOOPFC{
				u[k][j][i] = s[2][k][j][i]/s[1][k][j][i] ;
      }
      
      dqy_calc(u,dq) ;
      
      // fl is really fl=Flux*dt/dx (because of mdot)
      LOOPT1j{
#if(SYMFORCEHD==0)
				if(vp[2][k][j][i] > 0.) 
					fl[2][k][j][i] = (u[k][jm1][i] + (dx[2][2][j]*G2(2,i) - vp[2][k][j][i])*dq[k][jm1][i])*mdot[2][k][j][i]*DS12(k,j,i);
				else
					fl[2][k][j][i] = (u[k][j][i] + (-dx[2][2][j]*G2(2,i) - vp[2][k][j][i])*dq[k][j][i])*mdot[2][k][j][i]*DS12(k,j,i);
#else
				Dtot=vp[2][k][j][i]*OARC12(k,j,i);
				if(Dtot > DTOTLIMIT) 
					fl[2][k][j][i] = (u[k][jm1][i] + (1.0 - Dtot)*dq[k][jm1][i]*dx[2][2][j]*G2(2,i))*mdot[2][k][j][i]*DS12(k,j,i);
				else if(Dtot< -DTOTLIMIT)
					fl[2][k][j][i] = (u[k][j][i] + (-1.0 - Dtot)*dq[k][j][i]*dx[2][2][j]*G2(2,i))*mdot[2][k][j][i]*DS12(k,j,i);
				else{
					fraczone=Dtot*INVDTOTLIMIT;
					ftemp=mdot[2][k][j][i]*DS12(k,j,i);
					ftemp2=dx[2][2][j]*G2(2,i);
					dqminus=(u[k][jm1][i] + (1.0 - DTOTLIMIT)*dq[k][jm1][i]*ftemp2)*ftemp;
					dqplus=(u[k][j][i] + (-1.0 + DTOTLIMIT)*dq[k][j][i]*ftemp2)*ftemp;
					fl[2][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
				}
#endif
      }
      
      LOOPC{
				s[2][k][j][i] += (fl[2][k][j][i]-fl[2][k][jp1][i])*OVOL1(k,j,i);
      }
      if(FORCEIEINTERNAL){
				floor_correct(2,5);
      }
      bound(NULL,NULL,2,0,0) ;
    }
}
  else{
    if(DOLOSSDIAG){
      LOOPT1j{
				fl[2][k][j][i] = z2e_2(s[2],k,j,i)/z2e_2(s[1],k,j,i)*mdot[2][k][j][i]*DS12(k,j,i); // division by interp may be a problem
      }
    }
  }

  
  
  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(2,2,fl);
  }




#if(RAD)
      /* then specific internal energy */
      LOOPFC{
                                uer[k][j][i] = er[k][j][i]/s[1][k][j][i] ;
      }

      dqy_calc(uer,dq) ;

      // fl is really fl=Flux*dt/dx (because of mdot)
      LOOPT1j{
#if(SYMFORCEHD==0)
                                if(vp[2][k][j][i] > 0.)
                                        fl[2][k][j][i] = (uer[k][jm1][i] + (dx[2][2][j]*G2(2,i) - vp[2][k][j][i])*dq[k][jm1][i])*mdot[2][k][j][i]*DS12(k,j,i);
                                else
                                        fl[2][k][j][i] = (uer[k][j][i] + (-dx[2][2][j]*G2(2,i) - vp[2][k][j][i])*dq[k][j][i])*mdot[2][k][j][i]*DS12(k,j,i);
#else
                                Dtot=vp[2][k][j][i]*OARC12(k,j,i);
                                if(Dtot > DTOTLIMIT)
                                        fl[2][k][j][i] = (uer[k][jm1][i] + (1.0 - Dtot)*dq[k][jm1][i]*dx[2][2][j]*G2(2,i))*mdot[2][k][j][i]*DS12(k,j,i);
                                else if(Dtot< -DTOTLIMIT)
                                        fl[2][k][j][i] = (uer[k][j][i] + (-1.0 - Dtot)*dq[k][j][i]*dx[2][2][j]*G2(2,i))*mdot[2][k][j][i]*DS12(k,j,i);
                                else{
                                        fraczone=Dtot*INVDTOTLIMIT;
                                        ftemp=mdot[2][k][j][i]*DS12(k,j,i);
                                        ftemp2=dx[2][2][j]*G2(2,i);
                                        dqminus=(uer[k][jm1][i] + (1.0 - DTOTLIMIT)*dq[k][jm1][i]*ftemp2)*ftemp;
                                        dqplus=(uer[k][j][i] + (-1.0 + DTOTLIMIT)*dq[k][j][i]*ftemp2)*ftemp;
                                        fl[2][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
                                }
#endif
      }

      LOOPC{
                                er[k][j][i] += (fl[2][k][j][i]-fl[2][k][jp1][i])*OVOL1(k,j,i);
      }
      if(FORCEIEINTERNAL){
                                floor_correct(2,5);
      }
     bvaler();
#endif    














  if(COMPDIM==2){
    if(transmagx2&&mag){
      
      // magnetic field: bz, update FLUX of field, not field, eq. 51, 52
      LOOPFC{
				Bzr[k][j][i] = v[2][3][k][j][i]/s[1][k][j][i] ;
      }
      
      dqy_calc(Bzr,dq) ; //same as dqvy with wcom=3
      
      LOOPT1j{
#if(SYMFORCEHD==0)
				if(vp[2][k][j][i] > 0)
					fl[2][k][j][i] = (Bzr[k][jm1][i] + (dx[2][2][j]*G2(2,i) - vp[2][k][j][i])*dq[k][jm1][i])*mdot[2][k][j][i];
				else
					fl[2][k][j][i] = (Bzr[k][j][i] + (-dx[2][2][j]*G2(2,i) - vp[2][k][j][i])*dq[k][j][i])*mdot[2][k][j][i];
#else
				Dtot=vp[2][k][j][i]*OARC12(k,j,i);
				if(Dtot > DTOTLIMIT) 
					fl[2][k][j][i] = (Bzr[k][jm1][i] + (1.0 - Dtot)*dq[k][jm1][i]*dx[2][2][j]*G2(2,i))*mdot[2][k][j][i];
				else if(Dtot < -DTOTLIMIT)
					fl[2][k][j][i] = (Bzr[k][j][i] + (-1.0 - Dtot)*dq[k][j][i]*dx[2][2][j]*G2(2,i))*mdot[2][k][j][i];
				else{
					fraczone=Dtot*INVDTOTLIMIT;
					ftemp=mdot[2][k][j][i];
					ftemp2=dx[2][2][j]*G2(2,i);
					dqminus=(Bzr[k][jm1][i] + (1.0 - DTOTLIMIT)*dq[k][jm1][i]*ftemp2)*ftemp;
					dqplus=(Bzr[k][j][i] + (-1.0 + DTOTLIMIT)*dq[k][j][i]*ftemp2)*ftemp;
					fl[2][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
				}
#endif
      }
      LOOPC{
				v[2][3][k][j][i] += (fl[2][k][j][i]-fl[2][k][jp1][i])*OARC32(k,j,i);
      }
      bound(NULL,NULL,0,2,3); // only 3-comp needed
      
      
      if(transmagx2&&mag&&DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
				hydro_flux_adv(2,3,fl);
      }
    }

  }


  if(transv2x2){
    
    /* vy */
    dqvy_calc(2,v[1],dqv) ;

    
    LOOPT2j{
      vpp = e2z_2(vp[2],k,j,i);
      mdotp = e2z_2(mdot[2],k,j,i);
#if(SYMFORCEHD==0)
      if(vpp > 0)
				fl[2][k][j][i] = (v[1][2][k][j][i] + (dx[1][2][j]*G2(2,i) - vpp)*dqv[2][k][j][i])*G2(2,i)*mdotp*DS32(k,j,i);
      else
				fl[2][k][j][i] = (v[1][2][k][jp1][i] + (-dx[1][2][j]*G2(2,i) - vpp)*dqv[2][k][jp1][i])*G2(2,i)*mdotp*DS32(k,j,i); // dqv[2][jp1] here needs v[1][2][k][j=N2+2][i] if periodicx2special.
#else
      Dtot=vpp*OARC32(k,j,i);
      if(Dtot > DTOTLIMIT)
				fl[2][k][j][i] = (v[1][2][k][j][i] + (1.0 - Dtot)*dqv[2][k][j][i]*dx[1][2][j]*G2(2,i))*G2(2,i)*mdotp*DS32(k,j,i);
      else if(Dtot< -DTOTLIMIT)
				fl[2][k][j][i] = (v[1][2][k][jp1][i] + (-1.0 - Dtot)*dqv[2][k][jp1][i]*dx[1][2][j]*G2(2,i))*G2(2,i)*mdotp*DS32(k,j,i);
      else{
				fraczone=Dtot*INVDTOTLIMIT;
				ftemp=G2(2,i)*mdotp*DS32(k,j,i);
				ftemp2=dx[1][2][j]*G2(2,i);
				dqminus=(v[1][2][k][j][i] + (1.0 - DTOTLIMIT)*dqv[2][k][j][i]*ftemp2)*ftemp;
				dqplus=(v[1][2][k][jp1][i] + (-1.0 + DTOTLIMIT)*dqv[2][k][jp1][i]*ftemp2)*ftemp;
				fl[2][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
      }
#endif
    }

    if(periodicx2special){ // instead of special dq code and vy variable on j=N2+2, or otherwise
      bound(NULL,fl,0,-4,2);
    }

    LOOPV2{
      p[2][k][j][i] += (fl[2][k][jm1][i]-fl[2][k][j][i])*OVOL3(k,j,i);
    }
  }
  else{
    if(DOLOSSDIAG){
      LOOPT2j{
				fl[2][k][j][i] =v[1][2][k][j][i]*G2(2,i)*mdot[2][k][j][i]*DS12(k,j,i);
      }
    }
  }
  

  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(2,6,fl);
  }


  if(transv1x2){

    /* vx */
    dqvy_calc(1,v[1],dqv) ;
    
    LOOPT3j{
      // below avgs reason for loopt1j over other bzones
      mdotp = z2e_1(mdot[2],k,j,i);
      vpp = z2e_1(vp[2],k,j,i);
#if(SYMFORCEHD==0)
      if(vpp > 0.) 
				fl[2][k][j][i] = (v[1][1][k][jm1][i] + (dx[2][2][j]*G2(1,i) - vpp)*dqv[1][k][jm1][i])*mdotp*DS22(k,j,i);
      else	
				fl[2][k][j][i] = (v[1][1][k][j][i] + (-dx[2][2][j]*G2(1,i) - vpp)*dqv[1][k][j][i])*mdotp*DS22(k,j,i);
#else
      Dtot=vpp*OARC12(k,j,i);
      if(Dtot > DTOTLIMIT) 
				fl[2][k][j][i] = (v[1][1][k][jm1][i] + (1.0 - Dtot)*dqv[1][k][jm1][i]*dx[2][2][j]*G2(1,i))*mdotp*DS22(k,j,i);
      else if(Dtot < -DTOTLIMIT)
				fl[2][k][j][i] = (v[1][1][k][j][i] + (-1.0 - Dtot)*dqv[1][k][j][i]*dx[2][2][j]*G2(1,i))*mdotp*DS22(k,j,i);
      else{
				fraczone=Dtot*INVDTOTLIMIT;
				ftemp=mdotp*DS22(k,j,i);
				ftemp2=dx[2][2][j]*G2(1,i);
				dqminus=(v[1][1][k][jm1][i] + (1.0 - DTOTLIMIT)*dqv[1][k][jm1][i]*ftemp2)*ftemp;
				dqplus=(v[1][1][k][j][i] + (-1.0 + DTOTLIMIT)*dqv[1][k][j][i]*ftemp2)*ftemp;
				fl[2][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
      }
#endif
    }
    
    LOOPV1{
      p[1][k][j][i] += (fl[2][k][j][i]-fl[2][k][jp1][i])*OVOL2(k,j,i);
    }
  }
  else{
    if(DOLOSSDIAG){
      LOOPT3j{
				fl[2][k][j][i] =v1tov2(v[1][1],k,j,i)*mdot[2][k][j][i]*DS12(k,j,i);
      }
    }
  }

  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(2,5,fl);
  }


  if(transv3x2){
    
    /* vz */
    dqvy_calc(3,v[1],dqv) ;
    
    LOOPT4j{
      vpp = z2e_3(vp[2],k,j,i);
      mdotp = z2e_3(mdot[2],k,j,i);
#if(SYMFORCEHD==0)
      if(vpp > 0) 
				fl[2][k][j][i] = (v[1][3][k][jm1][i] + (dx[2][2][j]*G2(2,i) - vpp)*dqv[3][k][jm1][i])*G3(2,i)*G4(1,j)*mdotp*DS12(k,j,i);
      else
				fl[2][k][j][i] = (v[1][3][k][j][i] + (-dx[2][2][j]*G2(2,i) - vpp)*dqv[3][k][j][i])*G3(2,i)*G4(1,j)*mdotp*DS12(k,j,i);
#else
      Dtot = vpp*OARC12(k,j,i);
      if(Dtot > DTOTLIMIT) 
				fl[2][k][j][i] = (v[1][3][k][jm1][i] + (1.0 - Dtot)*dqv[3][k][jm1][i]*dx[2][2][j]*G2(2,i))*G3(2,i)*G4(1,j)*mdotp*DS12(k,j,i);
      else if(Dtot< -DTOTLIMIT)
				fl[2][k][j][i] = (v[1][3][k][j][i] + (-1.0 - Dtot)*dqv[3][k][j][i]*dx[2][2][j]*G2(2,i))*G3(2,i)*G4(1,j)*mdotp*DS12(k,j,i);
      else{
				fraczone=Dtot*INVDTOTLIMIT;
				ftemp=G3(2,i)*G4(1,j)*mdotp*DS12(k,j,i);
				ftemp2=dx[2][2][j]*G2(2,i);
				dqminus=(v[1][3][k][jm1][i] + (1.0 - DTOTLIMIT)*dqv[3][k][jm1][i]*ftemp2)*ftemp;
				dqplus=(v[1][3][k][j][i] + (-1.0 + DTOTLIMIT)*dqv[3][k][j][i]*ftemp2)*ftemp;
				fl[2][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
      }
#endif
    }
    
    LOOPV3{
      p[3][k][j][i] += (fl[2][k][j][i]-fl[2][k][jp1][i])*OVOL1(k,j,i);
    }//do not bound!
  }
  else{
    if(DOLOSSDIAG){
      LOOPT1j{
				fl[2][k][j][i] = v3tov2(v[1][3],k,j,i)*G3(2,i)*G4(1,j)*mdot[2][k][j][i]*DS12(k,j,i);
      }
    }
  }


  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(2,7,fl);
  }


  if(transrhox2){
    
    // fl is really fl=Flux*dt
    LOOPT1j{
      fl[2][k][j][i] =  mdot[2][k][j][i]*DS12(k,j,i);
    }
    
    
    LOOPC{
      s[1][k][j][i] += (fl[2][k][j][i] - fl[2][k][jp1][i])*OVOL1(k,j,i);
    }
    if(FORCERHOINTERNAL){
      floor_correct(1,5);
    }
    bound(NULL,NULL,1,0,0) ;
  }
  else{
    if(DOLOSSDIAG){
      LOOPT1j{
				fl[2][k][j][i] =  mdot[2][k][j][i]*DS12(k,j,i);
      }
    }
  }

  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(2,1,fl);
  }


  /* return to velocities */
  
  if(transv1x2){
    LOOPV1{
      v[1][1][k][j][i] = p[1][k][j][i]/z2e_1(s[1],k,j,i); // division by interp may be a problem
    }
  }

  if(transv2x2){
    LOOPV2{
      v[1][2][k][j][i] = p[2][k][j][i]/(z2e_2(s[1],k,j,i)*G2(2,i)) ; // division by interp may be a problem
    }
  }


  if(transv3x2){
    LOOPV3{
      v[1][3][k][j][i] = p[3][k][j][i]/(z2e_3(s[1],k,j,i)*G3(2,i)*G4(2,j));
    }
  }

  if( (transv1x2)||(transv2x2)||(transv3x2)){
    bound(NULL,NULL,0,1,123);  // bound all grid velocity-vector components
  }

}







void sweepz(void)
{
  FTYPE vpp,mdotp ;
  static FTYPE (*dq)[N2M][N1M],(*dqv)[N3M][N2M][N1M],(*vp)[N3M][N2M][N1M] ;
  static FTYPE (*mdot)[N3M][N2M][N1M] ;
  static FTYPE (*p)[N3M][N2M][N1M];
  static FTYPE (*u)[N2M][N1M],(*Bzr)[N2M][N1M];
  static FTYPE (*workv)[N3M][N2M][N1M];
  static FTYPE (*fl)[N3M][N2M][N1M];
  register int i,j,k ;
  FTYPE Dtot;
  FTYPE dqminus,dqplus,fraczone;

  FTYPE v1xa,v1ya,v1za,momold,momnew,ftemp;
  FTYPE ftempv;
  FTYPE ftemp0,ftemp1,ftemp2;
  FILE * out;
  int wcom;
  /* transform to momenta */


#if(DO_ACCRETOR)
  bound_accretor(3,NULL,NULL,-1,-1,123);
#endif


  // pointer assignments
  p    = workv1 ;
  vp   = workv2 ;
  fl   = workv3 ;
  mdot = workv4 ;
  dqv  = workv5 ;

  dq  = work1 ;
  u   = work2 ;
  Bzr = work3 ;

  LOOPC{ // no need for LOOPV here since multiplying and that's ok. 
    p[1][k][j][i] = z2e_1(s[1],k,j,i)*v[1][1][k][j][i] ;
    p[2][k][j][i] = z2e_2(s[1],k,j,i)*v[1][2][k][j][i]*G2(2,i) ;
    p[3][k][j][i] = z2e_3(s[1],k,j,i)*v[1][3][k][j][i]*G3(2,i)*G4(2,j) ;
  } // no need to bound


  LOOPHC{ // only need half full loop, no need to bound
    /* create spacial tranport variable */
    vp[3][k][j][i] = (v[1][3][k][j][i]-vg[3])*dt;
  }



  /* first do mass: rho */
  dqz_calc(s[1],dq) ;


  // mdot/fl is really fl=Flux*dt and mdot=Mdot*dt
  LOOPT0k{
#if(SYMFORCEHD==0)
    if(vp[3][k][j][i] > 0.){
      mdot[3][k][j][i] = (s[1][km1][j][i] + (dx[2][3][k]*G3(2,i)*G4(2,j) - vp[3][k][j][i])*dq[km1][j][i])*vp[3][k][j][i];
    }
    else{
      mdot[3][k][j][i] = (s[1][k][j][i] + (-dx[2][3][k]*G3(2,i)*G4(2,j) - vp[3][k][j][i])*dq[k][j][i])*vp[3][k][j][i];
    }
#else
    Dtot=vp[3][k][j][i]*OARC13(k,j,i)*ODX(2,3,k);
    if(Dtot > DTOTLIMIT){
      mdot[3][k][j][i] = (s[1][km1][j][i] + (1.0 - Dtot)*dq[km1][j][i]*dx[2][3][k]*G3(2,i)*G4(2,j))*vp[3][k][j][i];
    }
    else if(Dtot< -DTOTLIMIT) {
      mdot[3][k][j][i] = (s[1][k][j][i] + (-1.0 - Dtot)*dq[k][j][i]*dx[2][3][k]*G3(2,i)*G4(2,j))*vp[3][k][j][i];
    }
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      ftemp=vp[3][k][j][i];
      ftemp2=dx[2][3][k]*G3(2,i)*G4(2,j);
      dqminus=(s[1][km1][j][i] + (1.0 - DTOTLIMIT)*dq[km1][j][i]*ftemp2)*ftemp;
      dqplus=(s[1][k][j][i] + (-1.0 + DTOTLIMIT)*dq[k][j][i]*ftemp2)*ftemp;
      mdot[3][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
#endif
  }


  if((skipix3==0)||(ncpux3>1)){
    bound(NULL,mdot,0,-3,3);
  }


#include "sweeppassive3.h"



  if(transiex3){
    if(wgam){
      
      /* then specific internal energy */
      LOOPFC{
				u[k][j][i] = s[2][k][j][i]/s[1][k][j][i] ;
      }
      
      dqz_calc(u,dq) ;
      
      // fl is really fl=Flux*dt/dx (because of mdot)
      LOOPT1k{
#if(SYMFORCEHD==0)
				if(vp[3][k][j][i] > 0.) 
					fl[3][k][j][i] = (u[km1][j][i] + (dx[2][3][k]*G3(2,i)*G4(2,j) - vp[3][k][j][i])*dq[km1][j][i])*mdot[3][k][j][i]*DS13(k,j,i);
				else
					fl[3][k][j][i] = (u[k][j][i] + (-dx[2][3][k]*G3(2,i)*G4(2,j) - vp[3][k][j][i])*dq[k][j][i])*mdot[3][k][j][i]*DS13(k,j,i);
#else
				Dtot=vp[3][k][j][i]*OARC13(k,j,i)*ODX(2,3,k);
				if(Dtot > DTOTLIMIT) 
					fl[3][k][j][i] = (u[km1][j][i] + (1.0 - Dtot)*dq[km1][j][i]*dx[2][3][k]*G3(2,i)*G4(2,j))*mdot[3][k][j][i]*DS13(k,j,i);
				else if(Dtot< -DTOTLIMIT)
					fl[3][k][j][i] = (u[k][j][i] + (-1.0 - Dtot)*dq[k][j][i]*dx[2][3][k]*G3(2,i)*G4(2,j))*mdot[3][k][j][i]*DS13(k,j,i);
				else{
					fraczone=Dtot*INVDTOTLIMIT;
					ftemp=mdot[3][k][j][i]*DS13(k,j,i);
					ftemp2=dx[2][3][k]*G3(2,i)*G4(2,j);
					dqminus=(u[km1][j][i] + (1.0 - DTOTLIMIT)*dq[km1][j][i]*ftemp2)*ftemp;
					dqplus=(u[k][j][i] + (-1.0 + DTOTLIMIT)*dq[k][j][i]*ftemp2)*ftemp;
					fl[3][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
				}
#endif
      }
      
      LOOPC{
				s[2][k][j][i] += (fl[3][k][j][i]-fl[3][kp1][j][i])*OVOL1(k,j,i)*ODX(1,3,k); // no cancelling here unlike in sweepx and sweepy	
      }
      if(FORCEIEINTERNAL){
				floor_correct(2,6);
      }
      bound(NULL,NULL,2,0,0) ;
    }
  }
  else{
    if(DOLOSSDIAG){
      LOOPT1k{
				fl[3][k][j][i] = z2e_3(s[2],k,j,i)/z2e_3(s[1],k,j,i)*mdot[3][k][j][i]*DS13(k,j,i); // division by interp may be a problem
      }
    }
  }

  
  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(3,2,fl);
  }







  if(transv3x3){
    
    /* vz */
    dqvz_calc(3,v[1],dqv) ;
    
    LOOPT2k{
      vpp = e2z_3(vp[3],k,j,i);
      mdotp = e2z_3(mdot[3],k,j,i);
#if(SYMFORCEHD==0)
      if(vpp > 0) 
				fl[3][k][j][i] = (v[1][3][k][j][i] + (dx[1][3][k]*G3(2,i)*G4(2,j) - vpp)*dqv[3][k][j][i])*G3(2,i)*G4(2,j)*mdotp*DS13(k,j,i);
      else
				fl[3][k][j][i] = (v[1][3][kp1][j][i] + (-dx[1][3][k]*G3(2,i)*G4(2,j) - vpp)*dqv[3][kp1][j][i])*G3(2,i)*G4(2,j)*mdotp*DS13(k,j,i);
#else
      Dtot = vpp*OARC13(k,j,i)*ODX(1,3,k);
      if(Dtot > DTOTLIMIT) 
				fl[3][k][j][i] = (v[1][3][k][j][i] + (1.0 - Dtot)*dqv[3][k][j][i]*dx[2][3][k]*G3(2,i)*G4(2,j))*G3(2,i)*G4(2,j)*mdotp*DS13(k,j,i);
      else if(Dtot< -DTOTLIMIT)
				fl[3][k][j][i] = (v[1][3][kp1][j][i] + (-1.0 - Dtot)*dqv[3][kp1][j][i]*dx[2][3][k]*G3(2,i)*G4(2,j))*G3(2,i)*G4(2,j)*mdotp*DS13(k,j,i);
      else{
				fraczone=Dtot*INVDTOTLIMIT;
				ftemp=G3(2,i)*G4(2,j)*mdotp*DS13(k,j,i);
				ftemp2=dx[1][3][k]*G3(2,i)*G4(2,j);
				dqminus=(v[1][3][k][j][i] + (1.0 - DTOTLIMIT)*dqv[3][k][j][i]*ftemp2)*ftemp;
				dqplus=(v[1][3][kp1][j][i] + (-1.0 + DTOTLIMIT)*dqv[3][kp1][j][i]*ftemp2)*ftemp;
				fl[3][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
      }
#endif
    }

    LOOPV3{
      p[3][k][j][i] += (fl[3][km1][j][i]-fl[3][k][j][i])*OVOL1(k,j,i)*ODX(2,3,k);
    }//do not bound!
  }
  else{
    if(DOLOSSDIAG){
      LOOPT2k{
				fl[3][k][j][i] = v[1][3][k][j][i]*G3(2,i)*G4(2,j)*mdot[3][k][j][i]*DS13(k,j,i);
      }
    }
  }


  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(3,7,fl);
  }





  if(transv1x3){

    /* vx */
    dqvz_calc(1,v[1],dqv) ;
    
    LOOPT3k{
      // below avgs reason for loopt1j over other bzones
      mdotp = z2e_1(mdot[3],k,j,i);
      vpp = z2e_1(vp[3],k,j,i);
#if(SYMFORCEHD==0)
      if(vpp > 0.) 
				fl[3][k][j][i] = (v[1][1][km1][j][i] + (dx[2][3][k]*G3(1,i)*G4(2,j) - vpp)*dqv[1][km1][j][i])*mdotp*DS23(k,j,i);
      else	
				fl[3][k][j][i] = (v[1][1][k][j][i] + (-dx[2][3][k]*G3(1,i)*G4(2,j) - vpp)*dqv[1][k][j][i])*mdotp*DS23(k,j,i);
#else
      Dtot=vpp*OARC13(k,j,i)*ODX(2,3,k);
      if(Dtot > DTOTLIMIT) 
				fl[3][k][j][i] = (v[1][1][km1][j][i] + (1.0 - Dtot)*dqv[1][km1][j][i]*dx[2][3][k]*G3(1,i)*G4(2,j))*mdotp*DS23(k,j,i);
      else if(Dtot < -DTOTLIMIT)
				fl[3][k][j][i] = (v[1][1][k][j][i] + (-1.0 - Dtot)*dqv[1][k][j][i]*dx[2][3][k]*G3(1,i)*G4(2,j))*mdotp*DS23(k,j,i);
      else{
				fraczone=Dtot*INVDTOTLIMIT;
				ftemp=mdotp*DS23(k,j,i);
				ftemp2=dx[2][3][k]*G3(1,i)*G4(2,j);
				dqminus=(v[1][1][km1][j][i] + (1.0 - DTOTLIMIT)*dqv[1][km1][j][i]*ftemp2)*ftemp;
				dqplus=(v[1][1][k][j][i] + (-1.0 + DTOTLIMIT)*dqv[1][k][j][i]*ftemp2)*ftemp;
				fl[3][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
      }
#endif
    }
    
    LOOPV1{
      p[1][k][j][i] += (fl[3][k][j][i]-fl[3][kp1][j][i])*OVOL2(k,j,i)*ODX(1,3,k);
    }
  }
  else{
    if(DOLOSSDIAG){
      LOOPT3k{
				fl[3][k][j][i] =v1tov3(v[1][1],k,j,i)*mdot[3][k][j][i]*DS13(k,j,i);
      }
    }
  }

  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(3,5,fl);
  }

  if(transv2x3){
    
    /* vy */
    dqvz_calc(2,v[1],dqv) ;

    
    LOOPT4k{
      vpp = z2e_2(vp[3],k,j,i);
      mdotp = z2e_2(mdot[3],k,j,i);
#if(SYMFORCEHD==0)
      if(vpp > 0) 
				fl[3][k][j][i] = (v[1][2][km1][j][i] + (dx[2][3][k]*G3(2,i)*G4(1,j) - vpp)*dqv[2][km1][j][i])*G2(2,i)*mdotp*DS33(k,j,i);
      else
				fl[3][k][j][i] = (v[1][2][k][j][i] + (-dx[2][3][k]*G3(2,i)*G4(1,j) - vpp)*dqv[2][k][j][i])*G2(2,i)*mdotp*DS33(k,j,i);
#else
      Dtot=vpp*OARC33(k,j,i);
      if(Dtot > DTOTLIMIT)
				fl[3][k][j][i] = (v[1][2][km1][j][i] + (1.0 - Dtot)*dqv[2][km1][j][i]*dx[2][3][k]*G3(2,i)*G4(1,j))*G2(2,i)*mdotp*DS33(k,j,i);
      else if(Dtot< -DTOTLIMIT)
				fl[3][k][j][i] = (v[1][2][k][j][i] + (-1.0 - Dtot)*dqv[2][k][j][i]*dx[2][3][k]*G3(2,i)*G4(1,j))*G2(2,i)*mdotp*DS33(k,j,i);
      else{
				fraczone=Dtot*INVDTOTLIMIT;
				ftemp=G2(2,i)*mdotp*DS33(k,j,i);
				ftemp2=dx[2][3][k]*G3(2,i)*G4(1,j);
				dqminus=(v[1][2][km1][j][i] + (1.0 - DTOTLIMIT)*dqv[2][km1][j][i]*ftemp2)*ftemp;
				dqplus=(v[1][2][k][j][i] + (-1.0 + DTOTLIMIT)*dqv[2][k][j][i]*ftemp2)*ftemp;
				fl[3][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
      }
#endif
    }


    LOOPV2{
      p[2][k][j][i] += (fl[3][k][j][i]-fl[3][kp1][j][i])*OVOL3(k,j,i)*ODX(1,3,k);
    }
  }
  else{
    if(DOLOSSDIAG){
      LOOPT4k{
				fl[3][k][j][i] =v2tov3(v[1][2],k,j,i)*G2(2,i)*mdot[3][k][j][i]*DS13(k,j,i);
      }
    }
  }

  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(3,6,fl);
  }



  if(transrhox3){
    
    // fl is really fl=Flux*dt
    LOOPT1k{
      fl[3][k][j][i] =  mdot[3][k][j][i]*DS13(k,j,i);
    }
    
    
    LOOPC{
      s[1][k][j][i] += (fl[3][k][j][i] - fl[3][kp1][j][i])*OVOL1(k,j,i)*ODX(1,3,k);
    }
    if(FORCERHOINTERNAL){
      floor_correct(1,6);
    }
    bound(NULL,NULL,1,0,0) ;
  }
  else{
    if(DOLOSSDIAG){
      LOOPT1k{
				fl[3][k][j][i] =  mdot[3][k][j][i]*DS13(k,j,i);
      }
    }
  }

  if(DOLOSSDIAG&&(COMPUTELOSSDIAG==0)){
    hydro_flux_adv(3,1,fl);
  }


  /* return to velocities */
  
  if(transv1x3){
    LOOPV1{
      v[1][1][k][j][i] = p[1][k][j][i]/z2e_1(s[1],k,j,i); // division by interp may be a problem
    }
  }

  if(transv2x3){
    LOOPV2{
      v[1][2][k][j][i] = p[2][k][j][i]/(z2e_2(s[1],k,j,i)*G2(2,i)) ; // division by interp may be a problem
    }
  }


  if(transv3x3){
    LOOPV3{
      v[1][3][k][j][i] = p[3][k][j][i]/(z2e_3(s[1],k,j,i)*G3(2,i)*G4(2,j));
    }
  }

  if( (transv1x2)||(transv2x2)||(transv3x2)){
    bound(NULL,NULL,0,1,123);  // bound all grid velocity-vector components
  }

}
