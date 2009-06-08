#include "step.h"


// used so only have to compute sqrt of rho once
// 1d/2d/3d valid
void magprepare(void)
{
  register int i,j,k;
  FTYPE ftemp0,ftemp1,ftemp2,ftemp3,bxa,bya,bza;
  // for 4 times used, just do sqrt part once since costly compared to interp cost
  if(ALFVENLIMIT==0){
    LOOPHC{ // osqrtrho needed on -1..N, see below
      osqrtrho[k][j][i] = 1.0/sqrt(s[1][k][j][i]) ; // only used by magnetic terms
    }
  }
  else{ // osqrtrho only used by magnetic terms and pressure terms only use rholimited
    LOOPHC{ // for D1a's need osqrtrho at (-1,-1)...(N,N) (other dir(say 3 of 1,2) not checked yet)
      bxa = e2z_1(v[2][1],k,j,i);
      bya = e2z_2(v[2][2],k,j,i);
      bza = e2z_3(v[2][3],k,j,i);
      ftemp1=bxa*bxa; // bx^2
      ftemp2=bya*bya; // by^2
      ftemp3=bza*bza; // bz^2
      ftemp0=ftemp1+ftemp2+ftemp3; // b^2
      rholimited[0][k][j][i] = s[1][k][j][i]+ftemp0*invsol2; // rho+b^2/c_lim^2
      rholimited[1][k][j][i] = s[1][k][j][i]+(ftemp2+ftemp3)*invsol2; // rho+(b2^2+b3^2)/c^2
      rholimited[2][k][j][i] = s[1][k][j][i]+(ftemp1+ftemp3)*invsol2; // rho+(b1^2+b3^2)/c^2
      rholimited[3][k][j][i] = s[1][k][j][i]+(ftemp1+ftemp2)*invsol2; // rho+(b1^2+b2^2)/c^2
      osqrtrho[k][j][i] = 1.0/sqrt(rholimited[0][k][j][i]) ; // 1/sqrt(rho+b^2/c^2)
    }
  }
}



// 1d/2d/3d valid
void moc_ct_3d_v1(void)
{
  FTYPE fraczone,dqminus,dqplus;
  static FTYPE (*bxstar)[N2M][N1M],(*bystar)[N2M][N1M],(*bzstar)[N2M][N1M],(*vxstar)[N2M][N1M],(*vystar)[N2M][N1M] ,(*vzstar)[N2M][N1M];
  static FTYPE (*dqv)[N3M][N2M][N1M],(*dqb)[N3M][N2M][N1M];
  static FTYPE (*b1prim)[N3M][N2M][N1M],(*b2prim)[N3M][N2M][N1M];
  static FTYPE (*magemf)[N3M][N2M][N1M];
  static FTYPE bxa,bya,bza ;
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
  magemf = workv3 ;
  b1prim=workv4;
  b2prim=workv5;

  rhoa=work3;
  srhoa=work4;
  // work5,6,7,8,9,10 used below


  //////////////////////////////
  //////////
  ///////////    emf3
  ////////////
  ///////////  
  ////////////////////////////

  if(! ( (N1==1)&&(N2==1) ) ){ // otherwise no need for emf3 or b1prim[1] or b1prim[2]

    if(RHOINTERP==0){
      // used twice
      LOOPHPC{ // since srhoa needed (0..N), needs osqrtrho at -1..N
	srhoa[k][j][i] = z2c_3(osqrtrho,k,j,i) ; //1/sqrt(rho)
      }
    }
    
  // sweep in x-direction 
  // first get slopes 
  dqvx_calc(2,v[1],dqv) ;
  dqvx_calc(2,v[2],dqb) ;
    
    
  vystar = work6 ;
  bystar = work8 ;
    
  LOOPHPC{  // needed since emf needed here, and ok since dq[N1+1] not accessed which isn't created above
    ftemp=z2e_2(v[2][1],k,j,i);
    sgn_va=copysign(1.0,ftemp);
    ftemp=fabs(ftemp)*dt*OARC31(k,j,i);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i] ;
#else
    D1p=ftemp*z2e_2(osqrtrho,k,j,im1);
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
    srhop=z2e_2(osqrtrho,k,j,im1); // C+ char rho value
    bystar[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
    vystar[k][j][i] = 0.5*(vm+vp+sgn_va*bystar[k][j][i]*(srhop-srhom) +sgn_va*(-bp*srhop+bm*srhom) );
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif

    // find bprim, for updating velocities, no need for vyprime
    // lagrangian frame since op-split and other direction already moving
    // not in source code for more accurate evolution of alphven waves using upwinding


#if(N1>1)    // only need b1prim[2] if N1>1
    
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
    b1prim[2][k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    // same srhom/srhop from above
    //srhom=z2e_2(osqrtrho,k,j,i); // C- char rho value
    //srhop=z2e_2(osqrtrho,k,j,im1); // C+ char rho value
    b1prim[2][k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif


#endif // endif N1>1
  }
  // now vy/bystar live on corners and fill entire grid and border of comp grid




  // sweep in y-direction 
  // first get slopes 
  // dq's act like scalars on rotation
  dqvy_calc(1,v[1],dqv) ;
  dqvy_calc(1,v[2],dqb) ;


  bxstar = work7 ;
  vxstar = work5 ;

  LOOPHPC{ // needed since emf needed here, and ok since dq[N2+1] not accessed which isn't created above
    ftemp=z2e_1(v[2][2],k,j,i);
    sgn_va = copysign(1.,ftemp) ;
    ftemp=fabs(ftemp)*dt*OARC22(k,j,i);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i];
#else
    D1p=ftemp*z2e_1(osqrtrho,k,jm1,i);
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
    srhop=z2e_1(osqrtrho,k,jm1,i); // C+ char rho value
    bxstar[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
    vxstar[k][j][i] = 0.5*(vm+vp+sgn_va*bxstar[k][j][i]*(srhop-srhom) +sgn_va*(-bp*srhop+bm*srhom) );
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif



    // find bprim, to update velocities, no need for vxprim

#if(N2>1) // only need b1prim[1] if N2>1

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
    b1prim[1][k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    //srhom=z2e_1(osqrtrho,k,j,i); // C- char rho value
    //  srhop=z2e_1(osqrtrho,k,jm1,i); // C+ char rho value
    b1prim[1][k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif

#endif // endif N2>1
  }
  // now vx/bxstar live on corners on full grid, act like vectors on rotation


  // calculate emf=v cross B, e3 component is only that survives


  // emf lives on corners on full comp grid including border
  // emf is 3-vector, so acts like scalar on rotations
  if(RESMEM&&(res_real==1)){
    LOOPHPC{ // needed for below 2 loops.
      magemf[3][k][j][i] = vxstar[k][j][i]*bystar[k][j][i] - vystar[k][j][i]*bxstar[k][j][i] - jcurrent[3][k][j][i]*z2c_3(nu_res_real,k,j,i) ;
    }    
  }
  else{
    LOOPHPC{ // needed for below 2 loops.
      magemf[3][k][j][i] = vxstar[k][j][i]*bystar[k][j][i] - vystar[k][j][i]*bxstar[k][j][i] ;
    }
  }
  }// end if not true that both N2==1 and N1==1



  //////////////////////////////////
  //////////////////
  /////////////////   emf2
  ////////////////
  /////////////// 
  //////////////////////////////


  if(!( (N1==1)&&(N3==1) ) ){

#if(RHOINTERP==0)
    // used twice
    //  LOOPHC{ // needed since for rhointerp==1, need rho at i=-1..N1, j=-1..N2, so find osqrtrho on LOOPFMHPC!
  LOOPHPC{ // srhoa needed at 0..N, osqrtrho at -1..N
    srhoa[k][j][i] = z2c_2(osqrtrho,k,j,i) ; //1/sqrt(rho)
  }
#endif

  // sweep in x-direction (on 3-components)
  // first get slopes 
  dqvx_calc(3,v[1],dqv) ;
  dqvx_calc(3,v[2],dqb) ;



  vzstar = work6 ;
  bzstar = work8 ;

  LOOPHPC{  // needed since emf needed here, and ok since dq[N1+1] not accessed which isn't created above
    ftemp=z2e_3(v[2][1],k,j,i);
    sgn_va=copysign(1.0,ftemp);
    ftemp=fabs(ftemp)*dt*OARC11(k,j,i);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i] ;
#else
    D1p=ftemp*z2e_3(osqrtrho,k,j,im1);
    D1m=ftemp*z2e_3(osqrtrho,k,j,i);
#endif

    // determine which formula is for C+ and which is C-
    dist=dx[2][1][i];
    D2 = dt*z2e_3(v[1][1],k,j,i)*OARC11(k,j,i); // fraction of grid

    // Dtot violates symmetry since using fabs, causing mixed vector/scalar, otherwise symmetric like a vector    
    // vp and vm are "located" on corners on full comp grid+border
    // v=vx+|va| // trial Dtot, located on corner
#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    Dtota=D2+D1a;
    if(Dtota > 0) { // C+ from -i region
      vp = v[1][3][k][j][im1] + (1.0 - Dtota )*dqv[3][k][j][im1]*dist ;
      bp = v[2][3][k][j][im1] + (1.0 - Dtota )*dqb[3][k][j][im1]*dist ;
    }
    else{ // C+ is from + region
      vp = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }

    // v=vx-|va| // trial Dtot
    Dtota=D2-D1a;
    if(Dtota>0) { // C- is from -i region
      vm = v[1][3][k][j][im1] + (1.0 - Dtota )*dqv[3][k][j][im1]*dist ;
      bm = v[2][3][k][j][im1] + (1.0 - Dtota )*dqb[3][k][j][im1]*dist ;
    }
    else{ // C- is from +i region
      vm = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }
#else
    Dtota=D2+D1a;
    if(Dtota > DTOTLIMIT) { // C+ from -i region
      vp = v[1][3][k][j][im1] + (1.0 - Dtota )*dqv[3][k][j][im1]*dist ;
      bp = v[2][3][k][j][im1] + (1.0 - Dtota )*dqb[3][k][j][im1]*dist ;
    }
    else if( Dtota <-DTOTLIMIT) { // C+ is from + region
      vp = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }
    else{
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][3][k][j][im1] + (1.0 - DTOTLIMIT )*dqv[3][k][j][im1]*dist ;
      dqplus=v[1][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[3][k][j][i]*dist ;
      vp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][3][k][j][im1] + (1.0 - DTOTLIMIT )*dqb[3][k][j][im1]*dist ;
      dqplus=v[2][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[3][k][j][i]*dist ;
      bp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
    }

    // v=vx-|va| // trial Dtot
    Dtota=D2-D1a;
    if(Dtota > DTOTLIMIT) { // C- is from -i region
      vm = v[1][3][k][j][im1] + (1.0 - Dtota )*dqv[3][k][j][im1]*dist ;
      bm = v[2][3][k][j][im1] + (1.0 - Dtota )*dqb[3][k][j][im1]*dist ;
    }
    else if(Dtota <-DTOTLIMIT) { // C- is from +i region
      vm = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // Dtotp and Dtotm must have been already evaluated to get here
      vp = 0.5*(v[1][3][k][j][im1] + (1.0 - Dtota )*dqv[3][k][j][im1]*dist+
                v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist) ;
      bp = 0.5*(v[2][3][k][j][im1] + (1.0 - Dtota )*dqb[3][k][j][im1]*dist+
                v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist) ;

#else
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][3][k][j][im1] + (1.0 - DTOTLIMIT )*dqv[3][k][j][im1]*dist ;
      dqplus=v[1][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[3][k][j][i]*dist ;
      vm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][3][k][j][im1] + (1.0 - DTOTLIMIT )*dqb[3][k][j][im1]*dist ;
      dqplus=v[2][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[3][k][j][i]*dist ;
      bm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }

#endif

#elif(RHOINTERP==1)
    // no strong sym
#if(SYMFORCEMAG==0)
    if( (Dtotp=D2+D1p) > 0) { // C+ from -i region
      vp = v[1][3][k][j][im1] + (1.0 - Dtotp )*dqv[3][k][j][im1]*dist ;
      bp = v[2][3][k][j][im1] + (1.0 - Dtotp )*dqb[3][k][j][im1]*dist ;
    }
    else{ // C+ is from + region
      Dtotm=D2+D1m;
      vp = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }

    // v=vx-|va| // trial Dtot
    if( (Dtotp=D2-D1p) > 0) { // C- is from -i region
      vm = v[1][3][k][j][im1] + (1.0 - Dtotp )*dqv[3][k][j][im1]*dist ;
      bm = v[2][3][k][j][im1] + (1.0 - Dtotp )*dqb[3][k][j][im1]*dist ;
    }
    else{ // C- is from +i region
      Dtotm=D2-D1m;
      vm = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }

#else

    if( (Dtotp=D2+D1p) > DTOTLIMIT) { // C+ from -i region
      vp = v[1][3][k][j][im1] + (1.0 - Dtotp )*dqv[3][k][j][im1]*dist ;
      bp = v[2][3][k][j][im1] + (1.0 - Dtotp )*dqb[3][k][j][im1]*dist ;
    }
    else if( (Dtotm=D2+D1m) <-DTOTLIMIT) { // C+ is from + region
      vp = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }
    else{
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // Dtotp and Dtotm must have been already evaluated to get here
      vp = 0.5*(v[1][3][k][j][im1] + (1.0 - Dtotp )*dqv[3][k][j][im1]*dist+
		v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist) ;
      bp = 0.5*(v[2][3][k][j][im1] + (1.0 - Dtotp )*dqb[3][k][j][im1]*dist+
		v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist) ;
    }

    // v=vx-|va| // trial Dtot
    if( (Dtotp=D2-D1p) > DTOTLIMIT) { // C- is from -i region
      vm = v[1][3][k][j][im1] + (1.0 - Dtotp )*dqv[3][k][j][im1]*dist ;
      bm = v[2][3][k][j][im1] + (1.0 - Dtotp )*dqb[3][k][j][im1]*dist ;
    }
    else if( (Dtotm=D2-D1m) <-DTOTLIMIT) { // C- is from +i region
      vm = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }
    else{
      // Dtotp and Dtotm must have been already evaluated to get here
      vm = 0.5*(v[1][3][k][j][im1] + (1.0 - Dtotp )*dqv[3][k][j][im1]*dist+
		v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist) ;
      bm = 0.5*(v[2][3][k][j][im1] + (1.0 - Dtotp )*dqb[3][k][j][im1]*dist+
		v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist);
    }

#endif


#endif // end if rhointerp==1

		
    // should use midpoint of C+/C- for srho
    // solution from characteristics given the footprints from above
#if(RHOINTERP==0)
    // solution to Stone & Norman eqtn 43,44-- use constant rho at C+ C- crossing
    // eq 43,44 should really read: sign(va\dot k)dv=dB
    bzstar[k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
    vzstar[k][j][i] = 0.5*(vm + vp + sgn_va*(bm - bp)*srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    srhom=z2e_3(osqrtrho,k,j,i); // C- char rho value
    srhop=z2e_3(osqrtrho,k,j,im1); // C+ char rho value
    bzstar[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
    vzstar[k][j][i] = 0.5*(vm+vp+sgn_va*bzstar[k][j][i]*(srhop-srhom) +sgn_va*(-bp*srhop+bm*srhom) );
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif



    // find bprim, for updating velocities, no need for vyprime
    // lagrangian frame since op-split and other direction already moving
    // not in source code for more accurate evolution of alphven waves using upwinding

#if(N1>1) // otherwise no need for b1prim[3]
    
#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    // in this case D1p and D1m are equal to D1a
      vp = v[1][3][k][j][im1] + (1.0-D1a)*dqv[3][k][j][im1]*dist ; // C+
      bp = v[2][3][k][j][im1] + (1.0-D1a)*dqb[3][k][j][im1]*dist ; // C+
      vm = v[1][3][k][j][i] + (-1.0+D1a)*dqv[3][k][j][i]*dist ; // C-
      bm = v[2][3][k][j][i] + (-1.0+D1a)*dqb[3][k][j][i]*dist ; // C- 
#else
    // in this case D1p and D1m are equal to D1a
    if(D1a>DTOTLIMIT){ // keep p m names so like above
      vp = v[1][3][k][j][im1] + (1.0-D1a)*dqv[3][k][j][im1]*dist ; // C+
      bp = v[2][3][k][j][im1] + (1.0-D1a)*dqb[3][k][j][im1]*dist ; // C+
      vm = v[1][3][k][j][i] + (-1.0+D1a)*dqv[3][k][j][i]*dist ; // C-
      bm = v[2][3][k][j][i] + (-1.0+D1a)*dqb[3][k][j][i]*dist ; // C- 
    }
    else{
#if(SYMFUNKEDM1==0)
      // true, but only need difference in vp and vm for byprim
      //      vp = vm = 0.5*(v[1][3][k][j][im1] + (1.0-D1a)*dqv[3][k][j][im1]*dist+
      //		     v[1][3][k][j][i] + (-1.0+D1a)*dqv[3][k][j][i]*dist);
      vm = vp = 0.0; // since only need difference to be 0
      bp = bm = 0.5*(v[2][3][k][j][im1] + (1.0-D1a)*dqb[3][k][j][im1]*dist+
                     v[2][3][k][j][i] + (-1.0+D1a)*dqb[3][k][j][i]*dist);

#else
      // this is a little different than sweep.c
      fraczone=D1a*INVDTOTLIMIT;

      dqminus=v[1][3][k][j][im1] + (1.0-DTOTLIMIT)*dqv[3][k][j][im1]*dist ; 
      dqplus=v[1][3][k][j][i] + (-1.0+DTOTLIMIT)*dqv[3][k][j][i]*dist ;
      vm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      vp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));

      //      vm = vp = 0.0; // since only need difference to be 0

      dqminus=v[2][3][k][j][im1] + (1.0-DTOTLIMIT)*dqb[3][k][j][im1]*dist ;
      dqplus=v[2][3][k][j][i] + (-1.0+DTOTLIMIT)*dqb[3][k][j][i]*dist ;
      bm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      bp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
#endif
      //vp = v[1][3][k][j][im1] + (1.0-DTOTLIMIT)*dqv[3][k][j][im1]*dist ; // C+
      //vm = v[1][3][k][j][i] + (-1.0+DTOTLIMIT)*dqv[3][k][j][i]*dist ; // C-
    }
#endif

#elif(RHOINTERP==1)
    // no strong sym
#if(SYMFORCEMAG==0)
      vp = v[1][3][k][j][im1] + (1.0-D1p)*dqv[3][k][j][im1]*dist ; // C+
      bp = v[2][3][k][j][im1] + (1.0-D1p)*dqb[3][k][j][im1]*dist ; // C+
      vm = v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist ; // C-
      bm = v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist ; // C- 
#else
    if(D1p>DTOTLIMIT){
      vp = v[1][3][k][j][im1] + (1.0-D1p)*dqv[3][k][j][im1]*dist ; // C+
      bp = v[2][3][k][j][im1] + (1.0-D1p)*dqb[3][k][j][im1]*dist ; // C+
    }
    else{
      vp = 0.5*(v[1][3][k][j][im1] + (1.0-D1p)*dqv[3][k][j][im1]*dist+
		v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist);
      bp = 0.5*(v[2][3][k][j][im1] + (1.0-D1p)*dqb[3][k][j][im1]*dist+
		v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist);
    }
    if(D1m>DTOTLIMIT){
      vm = v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist ; // C-
      bm = v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist ; // C- 
    }
    else{
      vm = 0.5*(v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist+
		v[1][3][k][j][im1] + (1.0-D1p)*dqv[3][k][j][im1]*dist );
      bm = 0.5*(v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist+
		v[2][3][k][j][im1] + (1.0-D1p)*dqb[3][k][j][im1]*dist);
    }
#endif

#endif // end if rhointerp==1

#if(RHOINTERP==0)
    b1prim[3][k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    // same srhom/srhop from above
    //srhom=z2e_3(osqrtrho,k,j,i); // C- char rho value
    //srhop=z2e_3(osqrtrho,k,j,im1); // C+ char rho value
    b1prim[3][k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif


#endif // end if N1>1

  }
  // now vy/bystar live on corners and fill entire grid and border of comp grid




  // sweep in z-direction (on x-components)
  // first get slopes 
  // dq's act like scalars on rotation
  dqvz_calc(1,v[1],dqv) ;
  dqvz_calc(1,v[2],dqb) ;


  bxstar = work7 ;
  vxstar = work5 ;

  LOOPHPC{ // needed since emf needed here, and ok since dq[N2+1] not accessed which isn't created above
    ftemp=z2e_1(v[2][3],k,j,i);
    sgn_va = copysign(1.,ftemp) ;
    ftemp=fabs(ftemp)*dt*OARC23(k,j,i)*ODX(2,3,k);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i];
#else
    D1p=ftemp*z2e_1(osqrtrho,km1,j,i);
    D1m=ftemp*z2e_1(osqrtrho,k,j,i);
#endif

    // values at the foot of the plus characteristic
    //  are, by convention, in im1 zone
    D2 = dt*z2e_1(v[1][3],k,j,i)*OARC23(k,j,i)*ODX(2,3,k) ;
    dist=dx[2][3][k]*G3(1,i)*G4(2,j);


#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    // C+
    Dtota=D2+D1a;
    if(Dtota>0.) { // C+ from -k region
      vp = v[1][1][km1][j][i] + (1.0 - Dtota )*dqv[1][km1][j][i]*dist ;
      bp = v[2][1][km1][j][i] + (1.0 - Dtota )*dqb[1][km1][j][i]*dist ;
    }
    else{ // C+ from +k region
      vp = v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist ;
      bp = v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist ;
    }
    // C-
    Dtota=D2-D1a;
    if(Dtota>0.) { // C- from -k region
      vm = v[1][1][km1][j][i] + (1.0 - Dtota )*dqv[1][km1][j][i]*dist ;
      bm = v[2][1][km1][j][i] + (1.0 - Dtota )*dqb[1][km1][j][i]*dist ;
    }
    else{ // C- from +k region
      vm = v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist ;
    }
#else
    // C+
    Dtota=D2+D1a;
    if( Dtota > DTOTLIMIT) { // C+ from -k region
      vp = v[1][1][km1][j][i] + (1.0 - Dtota )*dqv[1][km1][j][i]*dist ;
      bp = v[2][1][km1][j][i] + (1.0 - Dtota )*dqb[1][km1][j][i]*dist ;
    }
    else if( Dtota < -DTOTLIMIT) { // C+ from +k region
      vp = v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist ;
      bp = v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      vp = 0.5*(v[1][1][km1][j][i] + (1.0 - Dtota )*dqv[1][km1][j][i]*dist+
                v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist);
      bp = 0.5*(v[2][1][km1][j][i] + (1.0 - Dtota )*dqb[1][km1][j][i]*dist+
                v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist);
#else
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][1][km1][j][i] + (1.0 - DTOTLIMIT )*dqv[1][km1][j][i]*dist ;
      dqplus=v[1][1][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[1][k][j][i]*dist ;
      vp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][1][km1][j][i] + (1.0 - DTOTLIMIT )*dqb[1][km1][j][i]*dist ;
      dqplus=v[2][1][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[1][k][j][i]*dist ;
      bp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }
    // C-
    Dtota=D2-D1a;
    if(Dtota > DTOTLIMIT) { // C- from -k region
      vm = v[1][1][km1][j][i] + (1.0 - Dtota )*dqv[1][km1][j][i]*dist ;
      bm = v[2][1][km1][j][i] + (1.0 - Dtota )*dqb[1][km1][j][i]*dist ;
    }
    else if(Dtota < -DTOTLIMIT) { // C- from +k region
      vm = v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      vm = 0.5*(v[1][1][km1][j][i] + (1.0 - Dtota )*dqv[1][km1][j][i]*dist+
                v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist);
      bm = 0.5*(v[2][1][km1][j][i] + (1.0 - Dtota )*dqb[1][km1][j][i]*dist+
                v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist);
#else
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][1][km1][j][i] + (1.0 - DTOTLIMIT )*dqv[1][km1][j][i]*dist ;
      dqplus=v[1][1][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[1][k][j][i]*dist ;
      vm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][1][km1][j][i] + (1.0 - DTOTLIMIT )*dqb[1][km1][j][i]*dist ;
      dqplus=v[2][1][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[1][k][j][i]*dist ;
      bm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
#if(SYMFORCEMAG==0)
    // C+
    if( (Dtotp=D2+D1p) > 0.) { // C+ from -k region
      vp = v[1][1][km1][j][i] + (1.0 - Dtotp )*dqv[1][km1][j][i]*dist ;
      bp = v[2][1][km1][j][i] + (1.0 - Dtotp )*dqb[1][km1][j][i]*dist ;
    }
    else{ // C+ from +k region
      Dtotm=D2+D1m;
      vp = v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist ;
      bp = v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist ;
    }
    // C-
    if( (Dtotp=D2-D1p) > 0.) { // C- from -k region
      vm = v[1][1][km1][j][i] + (1.0 - Dtotp )*dqv[1][km1][j][i]*dist ;
      bm = v[2][1][km1][j][i] + (1.0 - Dtotp )*dqb[1][km1][j][i]*dist ;
    }
    else{ // C- from +k region
      Dtotm=D2-D1m;
      vm = v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist ;
    }
#else
    // C+
    if( (Dtotp=D2+D1p) > DTOTLIMIT) { // C+ from -k region
      vp = v[1][1][km1][j][i] + (1.0 - Dtotp )*dqv[1][km1][j][i]*dist ;
      bp = v[2][1][km1][j][i] + (1.0 - Dtotp )*dqb[1][km1][j][i]*dist ;
    }
    else if( (Dtotm=D2+D1m) < -DTOTLIMIT) { // C+ from +k region
      vp = v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist ;
      bp = v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist ;
    }
    else{
      vp = 0.5*(v[1][1][km1][j][i] + (1.0 - Dtotp )*dqv[1][km1][j][i]*dist+
		v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist);
      bp = 0.5*(v[2][1][km1][j][i] + (1.0 - Dtotp )*dqb[1][km1][j][i]*dist+
		v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist);
    }
    // C-
    if( (Dtotp=D2-D1p) > DTOTLIMIT) { // C- from -k region
      vm = v[1][1][km1][j][i] + (1.0 - Dtotp )*dqv[1][km1][j][i]*dist ;
      bm = v[2][1][km1][j][i] + (1.0 - Dtotp )*dqb[1][km1][j][i]*dist ;
    }
    else if( (Dtotm=D2-D1m) < -DTOTLIMIT) { // C- from +k region
      vm = v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist ;
    }
    else{
      vm = 0.5*(v[1][1][km1][j][i] + (1.0 - Dtotp )*dqv[1][km1][j][i]*dist+
		v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist);
      bm = 0.5*(v[2][1][km1][j][i] + (1.0 - Dtotp )*dqb[1][km1][j][i]*dist+
		v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist);
    }
#endif

#endif    // end if rhointerp==1

    // bx is really bz and vx is really vz here
    // solution to Stone & Norman eqtn 43,44
#if(RHOINTERP==0)
    //-- use constant rho
    bxstar[k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
    vxstar[k][j][i] = 0.5*(vm + vp + sgn_va*(bm - bp)*srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    srhom=z2e_1(osqrtrho,k,j,i); // C- char rho value
    srhop=z2e_1(osqrtrho,km1,j,i); // C+ char rho value
    bxstar[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
    vxstar[k][j][i] = 0.5*(vm+vp+sgn_va*bxstar[k][j][i]*(srhop-srhom) +sgn_va*(-bp*srhop+bm*srhom) );
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif



    // find bprim, to update velocities, no need for vxprim

#if(N3>1) // otherwise b2prim[1] not needed

#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    vp = v[1][1][km1][j][i] + (1.0-D1a)*dqv[1][km1][j][i]*dist ;
    bp = v[2][1][km1][j][i] + (1.0-D1a)*dqb[1][km1][j][i]*dist ;
    vm = v[1][1][k][j][i] + (-1.0+D1a)*dqv[1][k][j][i]*dist ;
    bm = v[2][1][k][j][i] + (-1.0+D1a)*dqb[1][k][j][i]*dist ;
#else
    // D1p=D1m=D1a
    if(D1a>DTOTLIMIT){
      vp = v[1][1][km1][j][i] + (1.0-D1a)*dqv[1][km1][j][i]*dist ;
      bp = v[2][1][km1][j][i] + (1.0-D1a)*dqb[1][km1][j][i]*dist ;
      vm = v[1][1][k][j][i] + (-1.0+D1a)*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0+D1a)*dqb[1][k][j][i]*dist ;
    }
    else{
#if(SYMFUNKEDM2==0)
      // true, but only need diff (vm-vp) for bxprim
      //      vp = vm = 0.5*(v[1][1][km1][j][i] + (1.0-D1a)*dqv[1][km1][j][i]*dist+
      //                v[1][1][k][j][i] + (-1.0+D1a)*dqv[1][k][j][i]*dist);
      vp=vm=0.0;
      bp = bm = 0.5*(v[2][1][km1][j][i] + (1.0-D1a)*dqb[1][km1][j][i]*dist+
                     v[2][1][k][j][i] + (-1.0+D1a)*dqb[1][k][j][i]*dist);
#else
      // this is a little different than sweep.c
      fraczone=D1a*INVDTOTLIMIT;

      dqminus=v[1][1][km1][j][i] + (1.0-DTOTLIMIT)*dqv[1][km1][j][i]*dist ;
      dqplus=v[1][1][k][j][i] + (-1.0+DTOTLIMIT)*dqv[1][k][j][i]*dist ;
      vm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      vp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));

      dqminus=v[2][1][km1][j][i] + (1.0-DTOTLIMIT)*dqb[1][km1][j][i]*dist ;
      dqplus=v[2][1][k][j][i] + (-1.0+DTOTLIMIT)*dqb[1][k][j][i]*dist ;
      bm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      bp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
    // no strong sym
#if(SYMFORCEMAG==0)
      vp = v[1][1][km1][j][i] + (1.0-D1p)*dqv[1][km1][j][i]*dist ;
      bp = v[2][1][km1][j][i] + (1.0-D1p)*dqb[1][km1][j][i]*dist ;
      vm = v[1][1][k][j][i] + (-1.0+D1m)*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0+D1m)*dqb[1][k][j][i]*dist ;
#else
    if(D1p>DTOTLIMIT){
      vp = v[1][1][km1][j][i] + (1.0-D1p)*dqv[1][km1][j][i]*dist ;
      bp = v[2][1][km1][j][i] + (1.0-D1p)*dqb[1][km1][j][i]*dist ;
    }
    else{
      vp = 0.5*(v[1][1][km1][j][i] + (1.0-D1p)*dqv[1][km1][j][i]*dist+
		v[1][1][k][j][i] + (-1.0+D1m)*dqv[1][k][j][i]*dist);
      bp = 0.5*(v[2][1][km1][j][i] + (1.0-D1p)*dqb[1][km1][j][i]*dist+
		v[2][1][k][j][i] + (-1.0+D1m)*dqb[1][k][j][i]*dist);
    }
    if(D1m>DTOTLIMIT){
      vm = v[1][1][k][j][i] + (-1.0+D1m)*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0+D1m)*dqb[1][k][j][i]*dist ;
    }
    else{
      vm = 0.5*(v[1][1][km1][j][i] + (1.0-D1p)*dqv[1][km1][j][i]*dist+
		v[1][1][k][j][i] + (-1.0+D1m)*dqv[1][k][j][i]*dist);
      bm = 0.5*(v[2][1][km1][j][i] + (1.0-D1p)*dqb[1][km1][j][i]*dist+
		v[2][1][k][j][i] + (-1.0+D1m)*dqb[1][k][j][i]*dist);      
    }
#endif

#endif // end if rhointerp==1

#if(RHOINTERP==0)
    b2prim[1][k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    //srhom=z2e_1(osqrtrho,k,j,i); // C- char rho value
    //  srhop=z2e_1(osqrtrho,km1,j,i); // C+ char rho value
    b2prim[1][k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif

#endif // endif N3>1
  }
  // now vx/bxstar live on corners on full grid, act like vectors on rotation



  // calculate emf=v cross B, e3 component is only that survives


  // emf lives on corners on full comp grid including border
  // emf is 3-vector, so acts like scalar on rotations
  if(RESMEM&&(res_real==1)){
    LOOPHPC{ // needed for below 2 loops.
      magemf[2][k][j][i] = vzstar[k][j][i]*bxstar[k][j][i] - vxstar[k][j][i]*bzstar[k][j][i] - jcurrent[2][k][j][i]*z2c_2(nu_res_real,k,j,i) ;
    }    
  }
  else{
    LOOPHPC{ // needed for below 2 loops.
      magemf[2][k][j][i] = vzstar[k][j][i]*bxstar[k][j][i] - vxstar[k][j][i]*bzstar[k][j][i] ;
    }
  }
  }// end if not true that both N1==1 and N3==1  



  /////////////////////
  /////////////////
  /////////////////  emf1
  //////////////
  //////////////  
  //////////////////////////////


  if(!( (N2==1)&&(N3==1) ) ){ // otherwise no emf1 or b2prim[3] or b2prim[2] needed

    if(RHOINTERP==0){
      // used twice
      //  LOOPHC{ // needed since for rhointerp==1, need rho at i=-1..N1, j=-1..N2, so find osqrtrho on LOOPFMHPC!
      LOOPHPC{ // need srhoa: 0..N  need osqrtrho: -1..N
	srhoa[k][j][i] = z2c_1(osqrtrho,k,j,i) ; //1/sqrt(rho)
      }
    }
  // sweep in y-direction (on 3-components)
  // first get slopes 
  // dq's act like scalars on rotation
  dqvy_calc(3,v[1],dqv) ;
  dqvy_calc(3,v[2],dqb) ;


  vzstar = work6 ;
  bzstar = work8 ;

  LOOPHPC{ // needed since emf needed here, and ok since dq[N2+1] not accessed which isn't created above
    ftemp=z2e_3(v[2][2],k,j,i);
    sgn_va = copysign(1.,ftemp) ;
    ftemp=fabs(ftemp)*dt*OARC12(k,j,i);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i];
#else
    D1p=ftemp*z2e_3(osqrtrho,k,jm1,i);
    D1m=ftemp*z2e_3(osqrtrho,k,j,i);
#endif

    // values at the foot of the plus characteristic
    //  are, by convention, in im1 zone
    D2 = dt*z2e_3(v[1][2],k,j,i)*OARC12(k,j,i);
    dist=dx[2][2][j]*G2(2,i);


#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    // C+
    Dtota=D2+D1a;
    if(Dtota>0.) { // C+ from -j region
      vp = v[1][3][k][jm1][i] + (1.0 - Dtota )*dqv[3][k][jm1][i]*dist ;
      bp = v[2][3][k][jm1][i] + (1.0 - Dtota )*dqb[3][k][jm1][i]*dist ;
    }
    else{ // C+ from +j region
      vp = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }
    // C-
    Dtota=D2-D1a;
    if(Dtota>0.) { // C- from -j region
      vm = v[1][3][k][jm1][i] + (1.0 - Dtota )*dqv[3][k][jm1][i]*dist ;
      bm = v[2][3][k][jm1][i] + (1.0 - Dtota )*dqb[3][k][jm1][i]*dist ;
    }
    else{ // C- from +j region
      vm = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }
#else
    // C+
    Dtota=D2+D1a;
    if( Dtota > DTOTLIMIT) { // C+ from -j region
      vp = v[1][3][k][jm1][i] + (1.0 - Dtota )*dqv[3][k][jm1][i]*dist ;
      bp = v[2][3][k][jm1][i] + (1.0 - Dtota )*dqb[3][k][jm1][i]*dist ;
    }
    else if( Dtota < -DTOTLIMIT) { // C+ from +j region
      vp = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      vp = 0.5*(v[1][3][k][jm1][i] + (1.0 - Dtota )*dqv[3][k][jm1][i]*dist+
                v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist);
      bp = 0.5*(v[2][3][k][jm1][i] + (1.0 - Dtota )*dqb[3][k][jm1][i]*dist+
                v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist);
#else
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][3][k][jm1][i] + (1.0 - DTOTLIMIT )*dqv[3][k][jm1][i]*dist ;
      dqplus=v[1][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[3][k][j][i]*dist ;
      vp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][3][k][jm1][i] + (1.0 - DTOTLIMIT )*dqb[3][k][jm1][i]*dist ;
      dqplus=v[2][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[3][k][j][i]*dist ;
      bp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }
    // C-
    Dtota=D2-D1a;
    if(Dtota > DTOTLIMIT) { // C- from -j region
      vm = v[1][3][k][jm1][i] + (1.0 - Dtota )*dqv[3][k][jm1][i]*dist ;
      bm = v[2][3][k][jm1][i] + (1.0 - Dtota )*dqb[3][k][jm1][i]*dist ;
    }
    else if(Dtota < -DTOTLIMIT) { // C- from +j region
      vm = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      vm = 0.5*(v[1][3][k][jm1][i] + (1.0 - Dtota )*dqv[3][k][jm1][i]*dist+
                v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist);
      bm = 0.5*(v[2][3][k][jm1][i] + (1.0 - Dtota )*dqb[3][k][jm1][i]*dist+
                v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist);
#else
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][3][k][jm1][i] + (1.0 - DTOTLIMIT )*dqv[3][k][jm1][i]*dist ;
      dqplus=v[1][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[3][k][j][i]*dist ;
      vm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][3][k][jm1][i] + (1.0 - DTOTLIMIT )*dqb[3][k][jm1][i]*dist ;
      dqplus=v[2][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[3][k][j][i]*dist ;
      bm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
#if(SYMFORCEMAG==0)
    // C+
    if( (Dtotp=D2+D1p) > 0.) { // C+ from -j region
      vp = v[1][3][k][jm1][i] + (1.0 - Dtotp )*dqv[3][k][jm1][i]*dist ;
      bp = v[2][3][k][jm1][i] + (1.0 - Dtotp )*dqb[3][k][jm1][i]*dist ;
    }
    else{ // C+ from +j region
      Dtotm=D2+D1m;
      vp = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }
    // C-
    if( (Dtotp=D2-D1p) > 0.) { // C- from -j region
      vm = v[1][3][k][jm1][i] + (1.0 - Dtotp )*dqv[3][k][jm1][i]*dist ;
      bm = v[2][3][k][jm1][i] + (1.0 - Dtotp )*dqb[3][k][jm1][i]*dist ;
    }
    else{ // C- from +j region
      Dtotm=D2-D1m;
      vm = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }
#else
    // C+
    if( (Dtotp=D2+D1p) > DTOTLIMIT) { // C+ from -j region
      vp = v[1][3][k][jm1][i] + (1.0 - Dtotp )*dqv[3][k][jm1][i]*dist ;
      bp = v[2][3][k][jm1][i] + (1.0 - Dtotp )*dqb[3][k][jm1][i]*dist ;
    }
    else if( (Dtotm=D2+D1m) < -DTOTLIMIT) { // C+ from +j region
      vp = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }
    else{
      vp = 0.5*(v[1][3][k][jm1][i] + (1.0 - Dtotp )*dqv[3][k][jm1][i]*dist+
		v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist);
      bp = 0.5*(v[2][3][k][jm1][i] + (1.0 - Dtotp )*dqb[3][k][jm1][i]*dist+
		v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist);
    }
    // C-
    if( (Dtotp=D2-D1p) > DTOTLIMIT) { // C- from -j region
      vm = v[1][3][k][jm1][i] + (1.0 - Dtotp )*dqv[3][k][jm1][i]*dist ;
      bm = v[2][3][k][jm1][i] + (1.0 - Dtotp )*dqb[3][k][jm1][i]*dist ;
    }
    else if( (Dtotm=D2-D1m) < -DTOTLIMIT) { // C- from +j region
      vm = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }
    else{
      vm = 0.5*(v[1][3][k][jm1][i] + (1.0 - Dtotp )*dqv[3][k][jm1][i]*dist+
		v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist);
      bm = 0.5*(v[2][3][k][jm1][i] + (1.0 - Dtotp )*dqb[3][k][jm1][i]*dist+
		v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist);
    }
#endif

#endif    // end if rhointerp==1
    // solution to Stone & Norman eqtn 43,44
#if(RHOINTERP==0)
    //-- use constant rho
    bzstar[k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
    vzstar[k][j][i] = 0.5*(vm + vp + sgn_va*(bm - bp)*srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    srhom=z2e_3(osqrtrho,k,j,i); // C- char rho value
    srhop=z2e_3(osqrtrho,k,jm1,i); // C+ char rho value
    bzstar[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
    vzstar[k][j][i] = 0.5*(vm+vp+sgn_va*bzstar[k][j][i]*(srhop-srhom) +sgn_va*(-bp*srhop+bm*srhom) );
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif



    // find bprim, to update velocities, no need for vxprim

#if(N2>1) // otherwise no b2prim[3] needed

#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    vp = v[1][3][k][jm1][i] + (1.0-D1a)*dqv[3][k][jm1][i]*dist ;
    bp = v[2][3][k][jm1][i] + (1.0-D1a)*dqb[3][k][jm1][i]*dist ;
    vm = v[1][3][k][j][i] + (-1.0+D1a)*dqv[3][k][j][i]*dist ;
    bm = v[2][3][k][j][i] + (-1.0+D1a)*dqb[3][k][j][i]*dist ;
#else
    // D1p=D1m=D1a
    if(D1a>DTOTLIMIT){
      vp = v[1][3][k][jm1][i] + (1.0-D1a)*dqv[3][k][jm1][i]*dist ;
      bp = v[2][3][k][jm1][i] + (1.0-D1a)*dqb[3][k][jm1][i]*dist ;
      vm = v[1][3][k][j][i] + (-1.0+D1a)*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0+D1a)*dqb[3][k][j][i]*dist ;
    }
    else{
#if(SYMFUNKEDM2==0)
      // true, but only need diff (vm-vp) for bxprim
      //      vp = vm = 0.5*(v[1][3][k][jm1][i] + (1.0-D1a)*dqv[3][k][jm1][i]*dist+
      //                v[1][3][k][j][i] + (-1.0+D1a)*dqv[3][k][j][i]*dist);
      vp=vm=0.0;
      bp = bm = 0.5*(v[2][3][k][jm1][i] + (1.0-D1a)*dqb[3][k][jm1][i]*dist+
                     v[2][3][k][j][i] + (-1.0+D1a)*dqb[3][k][j][i]*dist);
#else
      // this is a little different than sweep.c
      fraczone=D1a*INVDTOTLIMIT;

      dqminus=v[1][3][k][jm1][i] + (1.0-DTOTLIMIT)*dqv[3][k][jm1][i]*dist ;
      dqplus=v[1][3][k][j][i] + (-1.0+DTOTLIMIT)*dqv[3][k][j][i]*dist ;
      vm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      vp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));

      dqminus=v[2][3][k][jm1][i] + (1.0-DTOTLIMIT)*dqb[3][k][jm1][i]*dist ;
      dqplus=v[2][3][k][j][i] + (-1.0+DTOTLIMIT)*dqb[3][k][j][i]*dist ;
      bm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      bp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
    // no strong sym
#if(SYMFORCEMAG==0)
      vp = v[1][3][k][jm1][i] + (1.0-D1p)*dqv[3][k][jm1][i]*dist ;
      bp = v[2][3][k][jm1][i] + (1.0-D1p)*dqb[3][k][jm1][i]*dist ;
      vm = v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist ;
#else
    if(D1p>DTOTLIMIT){
      vp = v[1][3][k][jm1][i] + (1.0-D1p)*dqv[3][k][jm1][i]*dist ;
      bp = v[2][3][k][jm1][i] + (1.0-D1p)*dqb[3][k][jm1][i]*dist ;
    }
    else{
      vp = 0.5*(v[1][3][k][jm1][i] + (1.0-D1p)*dqv[3][k][jm1][i]*dist+
		v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist);
      bp = 0.5*(v[2][3][k][jm1][i] + (1.0-D1p)*dqb[3][k][jm1][i]*dist+
		v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist);
    }
    if(D1m>DTOTLIMIT){
      vm = v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist ;
    }
    else{
      vm = 0.5*(v[1][3][k][jm1][i] + (1.0-D1p)*dqv[3][k][jm1][i]*dist+
		v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist);
      bm = 0.5*(v[2][3][k][jm1][i] + (1.0-D1p)*dqb[3][k][jm1][i]*dist+
		v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist);      
    }
#endif

#endif // end if rhointerp==1

#if(RHOINTERP==0)
    b2prim[3][k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    //srhom=z2e_3(osqrtrho,k,j,i); // C- char rho value
    //  srhop=z2e_3(osqrtrho,k,jm1,i); // C+ char rho value
    b2prim[3][k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif

#endif // end if N2>1

  }
  // now vx/bxstar live on corners on full grid, act like vectors on rotation


  // sweep in z-direction (on y-components)
  // first get slopes 
  // dq's act like scalars on rotation
  dqvz_calc(2,v[1],dqv) ;
  dqvz_calc(2,v[2],dqb) ;



  bystar = work7 ;
  vystar = work5 ;

  LOOPHPC{ // needed since emf needed here, and ok since dq[N2+1] not accessed which isn't created above
    ftemp=z2e_2(v[2][3],k,j,i);
    sgn_va = copysign(1.,ftemp) ;
    ftemp=fabs(ftemp)*dt*OARC33(k,j,i)*ODX(2,3,k);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i];
#else
    D1p=ftemp*z2e_2(osqrtrho,km1,j,i);
    D1m=ftemp*z2e_2(osqrtrho,k,j,i);
#endif

    // values at the foot of the plus characteristic
    //  are, by convention, in im1 zone
    D2 = dt*z2e_2(v[1][3],k,j,i)*OARC33(k,j,i)*ODX(2,3,k) ;
    dist=dx[2][3][k]*G3(2,i)*G4(1,j);


#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    // C+
    Dtota=D2+D1a;
    if(Dtota>0.) { // C+ from -k region
      vp = v[1][2][km1][j][i] + (1.0 - Dtota )*dqv[2][km1][j][i]*dist ;
      bp = v[2][2][km1][j][i] + (1.0 - Dtota )*dqb[2][km1][j][i]*dist ;
    }
    else{ // C+ from +k region
      vp = v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist ;
      bp = v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist ;
    }
    // C-
    Dtota=D2-D1a;
    if(Dtota>0.) { // C- from -k region
      vm = v[1][2][km1][j][i] + (1.0 - Dtota )*dqv[2][km1][j][i]*dist ;
      bm = v[2][2][km1][j][i] + (1.0 - Dtota )*dqb[2][km1][j][i]*dist ;
    }
    else{ // C- from +k region
      vm = v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist ;
    }
#else
    // C+
    Dtota=D2+D1a;
    if( Dtota > DTOTLIMIT) { // C+ from -k region
      vp = v[1][2][km1][j][i] + (1.0 - Dtota )*dqv[2][km1][j][i]*dist ;
      bp = v[2][2][km1][j][i] + (1.0 - Dtota )*dqb[2][km1][j][i]*dist ;
    }
    else if( Dtota < -DTOTLIMIT) { // C+ from +k region
      vp = v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist ;
      bp = v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      vp = 0.5*(v[1][2][km1][j][i] + (1.0 - Dtota )*dqv[2][km1][j][i]*dist+
                v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist);
      bp = 0.5*(v[2][2][km1][j][i] + (1.0 - Dtota )*dqb[2][km1][j][i]*dist+
                v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist);
#else
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][2][km1][j][i] + (1.0 - DTOTLIMIT )*dqv[2][km1][j][i]*dist ;
      dqplus=v[1][2][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[2][k][j][i]*dist ;
      vp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][2][km1][j][i] + (1.0 - DTOTLIMIT )*dqb[2][km1][j][i]*dist ;
      dqplus=v[2][2][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[2][k][j][i]*dist ;
      bp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }
    // C-
    Dtota=D2-D1a;
    if(Dtota > DTOTLIMIT) { // C- from -k region
      vm = v[1][2][km1][j][i] + (1.0 - Dtota )*dqv[2][km1][j][i]*dist ;
      bm = v[2][2][km1][j][i] + (1.0 - Dtota )*dqb[2][km1][j][i]*dist ;
    }
    else if(Dtota < -DTOTLIMIT) { // C- from +k region
      vm = v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      vm = 0.5*(v[1][2][km1][j][i] + (1.0 - Dtota )*dqv[2][km1][j][i]*dist+
                v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist);
      bm = 0.5*(v[2][2][km1][j][i] + (1.0 - Dtota )*dqb[2][km1][j][i]*dist+
                v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist);
#else
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][2][km1][j][i] + (1.0 - DTOTLIMIT )*dqv[2][km1][j][i]*dist ;
      dqplus=v[1][2][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[2][k][j][i]*dist ;
      vm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][2][km1][j][i] + (1.0 - DTOTLIMIT )*dqb[2][km1][j][i]*dist ;
      dqplus=v[2][2][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[2][k][j][i]*dist ;
      bm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
#if(SYMFORCEMAG==0)
    // C+
    if( (Dtotp=D2+D1p) > 0.) { // C+ from -k region
      vp = v[1][2][km1][j][i] + (1.0 - Dtotp )*dqv[2][km1][j][i]*dist ;
      bp = v[2][2][km1][j][i] + (1.0 - Dtotp )*dqb[2][km1][j][i]*dist ;
    }
    else{ // C+ from +k region
      Dtotm=D2+D1m;
      vp = v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist ;
      bp = v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist ;
    }
    // C-
    if( (Dtotp=D2-D1p) > 0.) { // C- from -k region
      vm = v[1][2][km1][j][i] + (1.0 - Dtotp )*dqv[2][km1][j][i]*dist ;
      bm = v[2][2][km1][j][i] + (1.0 - Dtotp )*dqb[2][km1][j][i]*dist ;
    }
    else{ // C- from +k region
      Dtotm=D2-D1m;
      vm = v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist ;
    }
#else
    // C+
    if( (Dtotp=D2+D1p) > DTOTLIMIT) { // C+ from -k region
      vp = v[1][2][km1][j][i] + (1.0 - Dtotp )*dqv[2][km1][j][i]*dist ;
      bp = v[2][2][km1][j][i] + (1.0 - Dtotp )*dqb[2][km1][j][i]*dist ;
    }
    else if( (Dtotm=D2+D1m) < -DTOTLIMIT) { // C+ from +k region
      vp = v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist ;
      bp = v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist ;
    }
    else{
      vp = 0.5*(v[1][2][km1][j][i] + (1.0 - Dtotp )*dqv[2][km1][j][i]*dist+
		v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist);
      bp = 0.5*(v[2][2][km1][j][i] + (1.0 - Dtotp )*dqb[2][km1][j][i]*dist+
		v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist);
    }
    // C-
    if( (Dtotp=D2-D1p) > DTOTLIMIT) { // C- from -k region
      vm = v[1][2][km1][j][i] + (1.0 - Dtotp )*dqv[2][km1][j][i]*dist ;
      bm = v[2][2][km1][j][i] + (1.0 - Dtotp )*dqb[2][km1][j][i]*dist ;
    }
    else if( (Dtotm=D2-D1m) < -DTOTLIMIT) { // C- from +k region
      vm = v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist ;
    }
    else{
      vm = 0.5*(v[1][2][km1][j][i] + (1.0 - Dtotp )*dqv[2][km1][j][i]*dist+
		v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist);
      bm = 0.5*(v[2][2][km1][j][i] + (1.0 - Dtotp )*dqb[2][km1][j][i]*dist+
		v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist);
    }
#endif

#endif    // end if rhointerp==1
    // solution to Stone & Norman eqtn 43,44
#if(RHOINTERP==0)
    //-- use constant rho
    bystar[k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
    vystar[k][j][i] = 0.5*(vm + vp + sgn_va*(bm - bp)*srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    srhom=z2e_2(osqrtrho,k,j,i); // C- char rho value
    srhop=z2e_2(osqrtrho,km1,j,i); // C+ char rho value
    bystar[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
    vystar[k][j][i] = 0.5*(vm+vp+sgn_va*bystar[k][j][i]*(srhop-srhom) +sgn_va*(-bp*srhop+bm*srhom) );
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif



    // find bprim, to update velocities, no need for vxprim

#if(N3>1) // otherwise no b2prim[2] needed

#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    vp = v[1][2][km1][j][i] + (1.0-D1a)*dqv[2][km1][j][i]*dist ;
    bp = v[2][2][km1][j][i] + (1.0-D1a)*dqb[2][km1][j][i]*dist ;
    vm = v[1][2][k][j][i] + (-1.0+D1a)*dqv[2][k][j][i]*dist ;
    bm = v[2][2][k][j][i] + (-1.0+D1a)*dqb[2][k][j][i]*dist ;
#else
    // D1p=D1m=D1a
    if(D1a>DTOTLIMIT){
      vp = v[1][2][km1][j][i] + (1.0-D1a)*dqv[2][km1][j][i]*dist ;
      bp = v[2][2][km1][j][i] + (1.0-D1a)*dqb[2][km1][j][i]*dist ;
      vm = v[1][2][k][j][i] + (-1.0+D1a)*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0+D1a)*dqb[2][k][j][i]*dist ;
    }
    else{
#if(SYMFUNKEDM2==0)
      // true, but only need diff (vm-vp) for bxprim
      //      vp = vm = 0.5*(v[1][2][km1][j][i] + (1.0-D1a)*dqv[2][km1][j][i]*dist+
      //                v[1][2][k][j][i] + (-1.0+D1a)*dqv[2][k][j][i]*dist);
      vp=vm=0.0;
      bp = bm = 0.5*(v[2][2][km1][j][i] + (1.0-D1a)*dqb[2][km1][j][i]*dist+
                     v[2][2][k][j][i] + (-1.0+D1a)*dqb[2][k][j][i]*dist);
#else
      // this is a little different than sweep.c
      fraczone=D1a*INVDTOTLIMIT;

      dqminus=v[1][2][km1][j][i] + (1.0-DTOTLIMIT)*dqv[2][km1][j][i]*dist ;
      dqplus=v[1][2][k][j][i] + (-1.0+DTOTLIMIT)*dqv[2][k][j][i]*dist ;
      vm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      vp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));

      dqminus=v[2][2][km1][j][i] + (1.0-DTOTLIMIT)*dqb[2][km1][j][i]*dist ;
      dqplus=v[2][2][k][j][i] + (-1.0+DTOTLIMIT)*dqb[2][k][j][i]*dist ;
      bm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      bp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
    // no strong sym
#if(SYMFORCEMAG==0)
      vp = v[1][2][km1][j][i] + (1.0-D1p)*dqv[2][km1][j][i]*dist ;
      bp = v[2][2][km1][j][i] + (1.0-D1p)*dqb[2][km1][j][i]*dist ;
      vm = v[1][2][k][j][i] + (-1.0+D1m)*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0+D1m)*dqb[2][k][j][i]*dist ;
#else
    if(D1p>DTOTLIMIT){
      vp = v[1][2][km1][j][i] + (1.0-D1p)*dqv[2][km1][j][i]*dist ;
      bp = v[2][2][km1][j][i] + (1.0-D1p)*dqb[2][km1][j][i]*dist ;
    }
    else{
      vp = 0.5*(v[1][2][km1][j][i] + (1.0-D1p)*dqv[2][km1][j][i]*dist+
		v[1][2][k][j][i] + (-1.0+D1m)*dqv[2][k][j][i]*dist);
      bp = 0.5*(v[2][2][km1][j][i] + (1.0-D1p)*dqb[2][km1][j][i]*dist+
		v[2][2][k][j][i] + (-1.0+D1m)*dqb[2][k][j][i]*dist);
    }
    if(D1m>DTOTLIMIT){
      vm = v[1][2][k][j][i] + (-1.0+D1m)*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0+D1m)*dqb[2][k][j][i]*dist ;
    }
    else{
      vm = 0.5*(v[1][2][km1][j][i] + (1.0-D1p)*dqv[2][km1][j][i]*dist+
		v[1][2][k][j][i] + (-1.0+D1m)*dqv[2][k][j][i]*dist);
      bm = 0.5*(v[2][2][km1][j][i] + (1.0-D1p)*dqb[2][km1][j][i]*dist+
		v[2][2][k][j][i] + (-1.0+D1m)*dqb[2][k][j][i]*dist);      
    }
#endif

#endif // end if rhointerp==1

#if(RHOINTERP==0)
    b2prim[2][k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    //srhom=z2e_2(osqrtrho,k,j,i); // C- char rho value
    //  srhop=z2e_2(osqrtrho,km1,j,i); // C+ char rho value
    b2prim[2][k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif

#endif // end if N3>1

  }
  // now yx/bystar live on corners on full grid, act like vectors on rotation

  // calculate emf=v cross B


  // emf lives on corners on full comp grid including border
  if(RESMEM&&(res_real==1)){
    LOOPHPC{ // needed for below 2 loops.
      magemf[1][k][j][i] =
	vystar[k][j][i]*bzstar[k][j][i] 
	- vzstar[k][j][i]*bystar[k][j][i] 
	- jcurrent[1][k][j][i]*z2c_1(nu_res_real,k,j,i) ;
    }    
  }
  else{
    LOOPHPC{ // needed for below 2 loops.
      magemf[1][k][j][i] = vystar[k][j][i]*bzstar[k][j][i] - vzstar[k][j][i]*bystar[k][j][i] ;
    }
  }
  }// end if not true that both N2==1 and N3==1



  // now all needed emfs are computed, can compute new field


  // now all needed prim quantities are computed, can compute transverse lorentz velocities

  // update velocity(2nd partial update of Lorentz force terms resulting in Alfven wave motion) (comoving)
  if(!( (N2==1)&&(N3==1))){
    if(mocctvx1){
      LOOPV1{ // needs b1prim[1]/b2prim[1] on i=0..N1-1 and j=0..N2, k=0..N3 (overdone, but ok)
#if(ALFVENLIMIT==0)
	rhoa2 = z2e_1(s[1],k,j,i);
#else
	rhoa2 = z2e_1(rholimited[0],k,j,i);
#endif
#if(N2!=1)
	bya = v2tov1(v[2][2],k,j,i);
#endif
#if(N3!=1)
	bza = v3tov1(v[2][3],k,j,i);
#endif
	v[1][1][k][j][i] += dt/rhoa2*(
#if(N2!=1)
				      (bya*OARC42(k,j,i)*(b1prim[1][k][jp1][i] - b1prim[1][k][j][i]))
#else
				      0
#endif
#if(N3!=1)
				      + (bza*OARC23(k,j,i)*ODX(1,3,k)*(b2prim[1][kp1][j][i]-b2prim[1][k][j][i]))
#endif
				      );
      }
    }
  }
  if(!( (N1==1)&&(N3==1) )){
    if(mocctvx2){
      LOOPV2{ // needs b1prim[2]/b2prim[2] at i=0..N1 and j=0..N2-1+periodicx2special (LOOPH above is fine) k=0..N3 (overdone, but ok)
#if(ALFVENLIMIT==0)
	rhoa2 =z2e_2(s[1],k,j,i);
#else
	rhoa2 =z2e_2(rholimited[0],k,j,i);
#endif
#if(N1!=1)
	bxa = v1tov2(v[2][1],k,j,i);
#endif
#if(N3!=1)
	bza = v3tov2(v[2][3],k,j,i);
#endif
	v[1][2][k][j][i] += dt/rhoa2*(
#if(N1!=1)
				      (bxa*OARC41(k,j,i)/(G2(2,i))*(G2(1,ip1)*b1prim[2][k][j][ip1] - G2(1,i)*b1prim[2][k][j][i]))
#else
				      0
#endif
#if(N3!=1)
				      +(bza*OARC33(k,j,i)*ODX(1,3,k)*(b2prim[2][kp1][j][i]-b2prim[2][k][j][i]))
#endif
				      );
      }
    }
  }
  if(!( (N1==1)&&(N2==1) )){
    if(mocctvx3){
      LOOPV3{ // needs b1prim[3]/b2prim[3] at k=0..N3-1 and i/j=0..N1/2 (overdone, but ok)
#if(ALFVENLIMIT==0)
	rhoa2 =z2e_3(s[1],k,j,i);
#else
	rhoa2 =z2e_3(rholimited[0],k,j,i);
#endif
#if(N1!=1)
	bxa = v1tov3(v[2][1],k,j,i);
#endif
#if(N2!=1)
	bya = v2tov3(v[2][2],k,j,i);
#endif
	v[1][3][k][j][i] += dt/rhoa2*(
#if(N1!=1)
				      (bxa*OARC41(k,j,i)/(G3(2,i))*(G3(1,ip1)*b1prim[3][k][j][ip1] - G3(1,i)*b1prim[3][k][j][i]))
#else
				      0
#endif
#if(N2!=1)
				      +(bya*OARC32(k,j,i)/(G4(2,j))*(G4(1,jp1)*b2prim[3][k][jp1][i]-G4(1,j)*b2prim[3][k][j][i]))
#endif
				      );
      }
    }
  }

  if(BOUNDFIELD==2){
    // BOUND EMFS SO DIVB=0 CONSERVED EXACTLY
    bound(NULL,magemf,0,-5,123);
  }


  // STONE AND MILLER
  // update field, CT with emf as conserved flux, hence divB=0.  eq 11,12 for B, eq 30,31 for v



  // must do mag field second since above v's use old B (like new Zeus says.  Old zeus says do field first)
  // newest zeus says split lorentz forces all before emf calcs, unlike in this _v1 function

  if(VOLUMEDIFF==0){ // old differential way

  if(!( (N2==1)&&(N3==1) ) ){
    if(mocctbx1){
      LOOPB1{ // needs emf on i=0..N1-1 and j=0..N2
	v[2][1][k][j][i] += dt*(
#if(N2!=1)
				OARC42(k,j,i)/G4(2,j)*(G4(1,jp1)*magemf[3][k][jp1][i] - G4(1,j)*magemf[3][k][j][i])
#else
				0
#endif
#if(N3!=1)
				+OARC23(k,j,i)*ODX(1,3,k)*(magemf[2][k][j][i]-magemf[2][kp1][j][i])
#endif
				) ;
      }
    }
  }
  if(!( (N1==1)&&(N3==1) ) ){
    if(mocctbx2){
      LOOPB2{ // needs emf at i=0..N1 and j=0..N2-1+periodicx2special (LOOPH above is fine)
	v[2][2][k][j][i] += dt*(
#if(N1!=1)
				OARC41(k,j,i)/G3(2,i)*(G3(1,i)*magemf[3][k][j][i] - G3(1,ip1)*magemf[3][k][j][ip1])
#else
				0
#endif
#if(N3!=1)
				+OARC33(k,j,i)*ODX(1,3,k)*(magemf[1][kp1][j][i]-magemf[1][k][j][i])
#endif
				) ;
	
      }
    }
  }
  if(!( (N1==1)&&(N2==1) )){
    if(mocctbx3){
      LOOPB3{ // needs emf at i=0..N1 and j=0..N2-1
	v[2][3][k][j][i] += dt*(
#if(N1!=1)
				OARC41(k,j,i)/G2(2,i)*(G2(1,ip1)*magemf[2][k][j][ip1] - G2(1,i)*magemf[2][k][j][i])
#else
				0
#endif
#if(N2!=1)
				+OARC32(k,j,i)*(magemf[1][k][j][i]-magemf[1][k][jp1][i])
#endif
				) ;
	
      }
    }
  }

  }
  else{ // volume way, which still preserves divB=0

  if(!( (N2==1)&&(N3==1) ) ){
    if(mocctbx1){
      LOOPB1{ // needs emf on i=0..N1-1 and j=0..N2
	v[2][1][k][j][i] += dt*(
#if(N2!=1)
				ODVL(1,2,j)/G2(1,i)*(G4(1,jp1)*magemf[3][k][jp1][i] - G4(1,j)*magemf[3][k][j][i])
#else
				0
#endif
#if(N3!=1) // no change for volume way
				+OARC23(k,j,i)*ODX(1,3,k)*(magemf[2][k][j][i]-magemf[2][kp1][j][i])
#endif
				) ;
      }
    }
  }
  if(!( (N1==1)&&(N3==1) ) ){
    if(mocctbx2){
      LOOPB2{ // needs emf at i=0..N1 and j=0..N2-1+periodicx2special (LOOPH above is fine)
	v[2][2][k][j][i] += dt*(
#if(N1!=1)
				G2(2,i)*ODVL(1,1,i)*(G3(1,i)*magemf[3][k][j][i] - G3(1,ip1)*magemf[3][k][j][ip1])
#else
				0
#endif
#if(N3!=1) // no change for volume way
				+OARC33(k,j,i)*ODX(1,3,k)*(magemf[1][kp1][j][i]-magemf[1][k][j][i])
#endif
				) ;
	
      }
    }
  }
  if(!( (N1==1)&&(N2==1) )){
    if(mocctbx3){
      LOOPB3{ // needs emf at i=0..N1 and j=0..N2-1
	v[2][3][k][j][i] += dt*(
#if(N1!=1)
				G3(2,i)*ODVL(1,1,i)*(G2(1,ip1)*magemf[2][k][j][ip1] - G2(1,i)*magemf[2][k][j][i])
#else
				0
#endif
#if(N2!=1)
				+G4(2,j)*ODVL(1,2,j)/G2(2,i)*(magemf[1][k][j][i]-magemf[1][k][jp1][i])
#endif
				) ;
	
      }
    }
  }


  } // end if volume way

  if((BOUNDTYPE>1)&&(BOUNDFIELD==1)){
  // limit what's bounded to only necessary changing items for faster reduced dimensional speeds(only relevant in 1D)
  // note that divB=0 correction in bound() is fine since other components don't change anyways(so didn't need to be updated, so are ready for divB=0 correction already)
  if(!( (N2==1)&&(N3==1) ) ){ // then do 1-comp
    if(!( (N1==1)&&(N3==1) ) ){ // then do 2-comp
      if(!( (N1==1)&&(N2==1) )){ // then do 3-comp
	bound(NULL,NULL,0,-1,123);
      }
      else{
	bound(NULL,NULL,0,-1,12); // no 3
      }
    }
    else{ // no 2
      if(!( (N1==1)&&(N2==1) )){
	bound(NULL,NULL,0,-1,13);
      }
      else{
	bound(NULL,NULL,0,-1,1); // no 3
      }
    }
  }
  else{// no 1
    if(!( (N1==1)&&(N3==1) ) ){ // then do 2-comp
      if(!( (N1==1)&&(N2==1) )){ // then do 3-comp
	bound(NULL,NULL,0,-1,23);
      }
      else{
	bound(NULL,NULL,0,-1,2); // no 3
      }
    }
    else{ // no 2
      if(!( (N1==1)&&(N2==1) )){
	bound(NULL,NULL,0,-1,3);
      }
      else{
	// never get here // no 3
      }
    }
  }
  }
  else if((BOUNDTYPE>1)&&(BOUNDFIELD==2)){ // must still bound velocities
  // limit what's bounded to only necessary changing items for faster reduced dimensional speeds(only relevant in 1D)
  // note that divB=0 correction in bound() is fine since other components don't change anyways(so didn't need to be updated, so are ready for divB=0 correction already)
  if(!( (N2==1)&&(N3==1) ) ){ // then do 1-comp
    if(!( (N1==1)&&(N3==1) ) ){ // then do 2-comp
      if(!( (N1==1)&&(N2==1) )){ // then do 3-comp
	bound(NULL,NULL,0,1,123);
      }
      else{
	bound(NULL,NULL,0,1,12); // no 3
      }
    }
    else{ // no 2
      if(!( (N1==1)&&(N2==1) )){
	bound(NULL,NULL,0,1,13);
      }
      else{
	bound(NULL,NULL,0,1,1); // no 3
      }
    }
  }
  else{// no 1
    if(!( (N1==1)&&(N3==1) ) ){ // then do 2-comp
      if(!( (N1==1)&&(N2==1) )){ // then do 3-comp
	bound(NULL,NULL,0,1,23);
      }
      else{
	bound(NULL,NULL,0,1,2); // no 3
      }
    }
    else{ // no 2
      if(!( (N1==1)&&(N2==1) )){
	bound(NULL,NULL,0,1,3);
      }
      else{
	// never get here // no 3
      }
    }
  }
  }
  else{// then just normal
  // limit what's bounded to only necessary changing items for faster reduced dimensional speeds(only relevant in 1D)
  // note that divB=0 correction in bound() is fine since other components don't change anyways(so didn't need to be updated, so are ready for divB=0 correction already)
  if(!( (N2==1)&&(N3==1) ) ){ // then do 1-comp
    if(!( (N1==1)&&(N3==1) ) ){ // then do 2-comp
      if(!( (N1==1)&&(N2==1) )){ // then do 3-comp
	bound(NULL,NULL,0,-1,123);
      }
      else{
	bound(NULL,NULL,0,-1,12); // no 3
      }
    }
    else{ // no 2
      if(!( (N1==1)&&(N2==1) )){
	bound(NULL,NULL,0,-1,13);
      }
      else{
	bound(NULL,NULL,0,-1,1); // no 3
      }
    }
  }
  else{// no 1
    if(!( (N1==1)&&(N3==1) ) ){ // then do 2-comp
      if(!( (N1==1)&&(N2==1) )){ // then do 3-comp
	bound(NULL,NULL,0,-1,23);
      }
      else{
	bound(NULL,NULL,0,-1,2); // no 3
      }
    }
    else{ // no 2
      if(!( (N1==1)&&(N2==1) )){
	bound(NULL,NULL,0,-1,3);
      }
      else{
	// never get here // no 3
      }
    }
  }
  }

  
}


// 1d/2d/3d valid
// this is just cut-pasted out of moc_ct_3d_v1()
// this version has no lorentz transverse calculations, as suggested by Hawley as in ZEUS-3D
// cut/paste from _v1 then:
// 1) comment out/delete prims allocated and used, and velocity calc
void moc_ct_3d_v2(void)
{
  FTYPE fraczone,dqminus,dqplus;
  static FTYPE (*bxstar)[N2M][N1M],(*bystar)[N2M][N1M],(*bzstar)[N2M][N1M],(*vxstar)[N2M][N1M],(*vystar)[N2M][N1M] ,(*vzstar)[N2M][N1M];
  static FTYPE (*dqv)[N3M][N2M][N1M],(*dqb)[N3M][N2M][N1M];
  static FTYPE (*magemf)[N3M][N2M][N1M];
  static FTYPE bxa,bya,bza ;
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



  if((FLOATTYPE==0)&&(SYMFORCEMAG==1)){
    fprintf(fail_file,"moc_ct: floats and symmetry are known to fail--no idea why yet\n");
    myexit(1);
  }

  dqv = workv1 ;
  dqb = workv2 ;
  magemf = workv3 ;

  rhoa=work3;
  srhoa=work4;
  // work5,6,7,8,9,10 used below



  //////////////////////////////
  //////////
  ///////////    emf3
  ////////////
  ///////////  
  ////////////////////////////

  if(! ( (N1==1)&&(N2==1) ) ){ // otherwise no need for emf3 or b1prim[1] or b1prim[2]

    if(RHOINTERP==0){
      // used twice
      LOOPHPC{ // see above
	// needed since for rhointerp==1, need rho at i=-1..N1, j=-1..N2, so find osqrtrho on LOOPFMHPC!
	srhoa[k][j][i] = z2c_3(osqrtrho,k,j,i) ; //1/sqrt(rho)
      }
    }
  // sweep in x-direction 
  // first get slopes 
  dqvx_calc(2,v[1],dqv) ;
  dqvx_calc(2,v[2],dqb) ;
    
    
  vystar = work6 ;
  bystar = work8 ;
    
  LOOPHPC{  // needed since emf needed here, and ok since dq[N1+1] not accessed which isn't created above
    ftemp=z2e_2(v[2][1],k,j,i);
    sgn_va=copysign(1.0,ftemp);
    ftemp=fabs(ftemp)*dt*OARC31(k,j,i);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i] ;
#else
    D1p=ftemp*z2e_2(osqrtrho,k,j,im1);
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
    srhop=z2e_2(osqrtrho,k,j,im1); // C+ char rho value
    bystar[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
    vystar[k][j][i] = 0.5*(vm+vp+sgn_va*bystar[k][j][i]*(srhop-srhom) +sgn_va*(-bp*srhop+bm*srhom) );
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


  bxstar = work7 ;
  vxstar = work5 ;

  LOOPHPC{ // needed since emf needed here, and ok since dq[N2+1] not accessed which isn't created above
    ftemp=z2e_1(v[2][2],k,j,i);
    sgn_va = copysign(1.,ftemp) ;
    ftemp=fabs(ftemp)*dt*OARC22(k,j,i);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i];
#else
    D1p=ftemp*z2e_1(osqrtrho,k,jm1,i);
    D1m=ftemp*z2e_1(osqrtrho,k,j,i);
#endif

    // values at the foot of the plus characteristic
    //  are, by convention, in im1 zone
    D2 = dt*z2e_1(v[1][2],k,j,i)*OARC22(k,j,i);
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
    srhop=z2e_1(osqrtrho,k,jm1,i); // C+ char rho value
    bxstar[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
    vxstar[k][j][i] = 0.5*(vm+vp+sgn_va*bxstar[k][j][i]*(srhop-srhom) +sgn_va*(-bp*srhop+bm*srhom) );
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
      magemf[3][k][j][i] = vxstar[k][j][i]*bystar[k][j][i] - vystar[k][j][i]*bxstar[k][j][i] - jcurrent[3][k][j][i]*z2c_3(nu_res_real,k,j,i) ;
    }    
  }
  else{
    LOOPHPC{ // needed for below 2 loops.
      magemf[3][k][j][i] = vxstar[k][j][i]*bystar[k][j][i] - vystar[k][j][i]*bxstar[k][j][i] ;
    }
  }
  }// end if not true that both N2==1 and N1==1





  //////////////////////////////////
  //////////////////
  /////////////////   emf2
  ////////////////
  /////////////// 
  //////////////////////////////


  if(!( (N1==1)&&(N3==1) ) ){

    if(RHOINTERP==0){
      // used twice
      LOOPHPC{ // see above
	// needed since for rhointerp==1, need rho at i=-1..N1, j=-1..N2, so find osqrtrho on LOOPFMHPC!
	srhoa[k][j][i] = z2c_2(osqrtrho,k,j,i) ; //1/sqrt(rho)
      }
    }


  // sweep in x-direction (on 3-components)
  // first get slopes 
  dqvx_calc(3,v[1],dqv) ;
  dqvx_calc(3,v[2],dqb) ;



  vzstar = work6 ;
  bzstar = work8 ;

  LOOPHPC{  // needed since emf needed here, and ok since dq[N1+1] not accessed which isn't created above
    ftemp=z2e_3(v[2][1],k,j,i);
    sgn_va=copysign(1.0,ftemp);
    ftemp=fabs(ftemp)*dt*OARC11(k,j,i);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i] ;
#else
    D1p=ftemp*z2e_3(osqrtrho,k,j,im1);
    D1m=ftemp*z2e_3(osqrtrho,k,j,i);
#endif

    // determine which formula is for C+ and which is C-
    dist=dx[2][1][i];
    D2 = dt*z2e_3(v[1][1],k,j,i)*OARC11(k,j,i); // fraction of grid

    // Dtot violates symmetry since using fabs, causing mixed vector/scalar, otherwise symmetric like a vector    
    // vp and vm are "located" on corners on full comp grid+border
    // v=vx+|va| // trial Dtot, located on corner
#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    Dtota=D2+D1a;
    if(Dtota > 0) { // C+ from -i region
      vp = v[1][3][k][j][im1] + (1.0 - Dtota )*dqv[3][k][j][im1]*dist ;
      bp = v[2][3][k][j][im1] + (1.0 - Dtota )*dqb[3][k][j][im1]*dist ;
    }
    else{ // C+ is from + region
      vp = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }

    // v=vx-|va| // trial Dtot
    Dtota=D2-D1a;
    if(Dtota>0) { // C- is from -i region
      vm = v[1][3][k][j][im1] + (1.0 - Dtota )*dqv[3][k][j][im1]*dist ;
      bm = v[2][3][k][j][im1] + (1.0 - Dtota )*dqb[3][k][j][im1]*dist ;
    }
    else{ // C- is from +i region
      vm = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }
#else
    Dtota=D2+D1a;
    if(Dtota > DTOTLIMIT) { // C+ from -i region
      vp = v[1][3][k][j][im1] + (1.0 - Dtota )*dqv[3][k][j][im1]*dist ;
      bp = v[2][3][k][j][im1] + (1.0 - Dtota )*dqb[3][k][j][im1]*dist ;
    }
    else if( Dtota <-DTOTLIMIT) { // C+ is from + region
      vp = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }
    else{
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][3][k][j][im1] + (1.0 - DTOTLIMIT )*dqv[3][k][j][im1]*dist ;
      dqplus=v[1][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[3][k][j][i]*dist ;
      vp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][3][k][j][im1] + (1.0 - DTOTLIMIT )*dqb[3][k][j][im1]*dist ;
      dqplus=v[2][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[3][k][j][i]*dist ;
      bp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
    }

    // v=vx-|va| // trial Dtot
    Dtota=D2-D1a;
    if(Dtota > DTOTLIMIT) { // C- is from -i region
      vm = v[1][3][k][j][im1] + (1.0 - Dtota )*dqv[3][k][j][im1]*dist ;
      bm = v[2][3][k][j][im1] + (1.0 - Dtota )*dqb[3][k][j][im1]*dist ;
    }
    else if(Dtota <-DTOTLIMIT) { // C- is from +i region
      vm = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // Dtotp and Dtotm must have been already evaluated to get here
      vp = 0.5*(v[1][3][k][j][im1] + (1.0 - Dtota )*dqv[3][k][j][im1]*dist+
                v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist) ;
      bp = 0.5*(v[2][3][k][j][im1] + (1.0 - Dtota )*dqb[3][k][j][im1]*dist+
                v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist) ;

#else
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][3][k][j][im1] + (1.0 - DTOTLIMIT )*dqv[3][k][j][im1]*dist ;
      dqplus=v[1][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[3][k][j][i]*dist ;
      vm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][3][k][j][im1] + (1.0 - DTOTLIMIT )*dqb[3][k][j][im1]*dist ;
      dqplus=v[2][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[3][k][j][i]*dist ;
      bm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }

#endif

#elif(RHOINTERP==1)
    // no strong sym
#if(SYMFORCEMAG==0)
    if( (Dtotp=D2+D1p) > 0) { // C+ from -i region
      vp = v[1][3][k][j][im1] + (1.0 - Dtotp )*dqv[3][k][j][im1]*dist ;
      bp = v[2][3][k][j][im1] + (1.0 - Dtotp )*dqb[3][k][j][im1]*dist ;
    }
    else{ // C+ is from + region
      Dtotm=D2+D1m;
      vp = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }

    // v=vx-|va| // trial Dtot
    if( (Dtotp=D2-D1p) > 0) { // C- is from -i region
      vm = v[1][3][k][j][im1] + (1.0 - Dtotp )*dqv[3][k][j][im1]*dist ;
      bm = v[2][3][k][j][im1] + (1.0 - Dtotp )*dqb[3][k][j][im1]*dist ;
    }
    else{ // C- is from +i region
      Dtotm=D2-D1m;
      vm = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }

#else

    if( (Dtotp=D2+D1p) > DTOTLIMIT) { // C+ from -i region
      vp = v[1][3][k][j][im1] + (1.0 - Dtotp )*dqv[3][k][j][im1]*dist ;
      bp = v[2][3][k][j][im1] + (1.0 - Dtotp )*dqb[3][k][j][im1]*dist ;
    }
    else if( (Dtotm=D2+D1m) <-DTOTLIMIT) { // C+ is from + region
      vp = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }
    else{
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // Dtotp and Dtotm must have been already evaluated to get here
      vp = 0.5*(v[1][3][k][j][im1] + (1.0 - Dtotp )*dqv[3][k][j][im1]*dist+
		v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist) ;
      bp = 0.5*(v[2][3][k][j][im1] + (1.0 - Dtotp )*dqb[3][k][j][im1]*dist+
		v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist) ;
    }

    // v=vx-|va| // trial Dtot
    if( (Dtotp=D2-D1p) > DTOTLIMIT) { // C- is from -i region
      vm = v[1][3][k][j][im1] + (1.0 - Dtotp )*dqv[3][k][j][im1]*dist ;
      bm = v[2][3][k][j][im1] + (1.0 - Dtotp )*dqb[3][k][j][im1]*dist ;
    }
    else if( (Dtotm=D2-D1m) <-DTOTLIMIT) { // C- is from +i region
      vm = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }
    else{
      // Dtotp and Dtotm must have been already evaluated to get here
      vm = 0.5*(v[1][3][k][j][im1] + (1.0 - Dtotp )*dqv[3][k][j][im1]*dist+
		v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist) ;
      bm = 0.5*(v[2][3][k][j][im1] + (1.0 - Dtotp )*dqb[3][k][j][im1]*dist+
		v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist);
    }

#endif


#endif // end if rhointerp==1

		
    // should use midpoint of C+/C- for srho
    // solution from characteristics given the footprints from above
#if(RHOINTERP==0)
    // solution to Stone & Norman eqtn 43,44-- use constant rho at C+ C- crossing
    // eq 43,44 should really read: sign(va\dot k)dv=dB
    bzstar[k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
    vzstar[k][j][i] = 0.5*(vm + vp + sgn_va*(bm - bp)*srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    srhom=z2e_3(osqrtrho,k,j,i); // C- char rho value
    srhop=z2e_3(osqrtrho,k,j,im1); // C+ char rho value
    bzstar[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
    vzstar[k][j][i] = 0.5*(vm+vp+sgn_va*bzstar[k][j][i]*(srhop-srhom) +sgn_va*(-bp*srhop+bm*srhom) );
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif


  }
  // now vy/bystar live on corners and fill entire grid and border of comp grid




  // sweep in z-direction (on x-components)
  // first get slopes 
  // dq's act like scalars on rotation
  dqvz_calc(1,v[1],dqv) ;
  dqvz_calc(1,v[2],dqb) ;


  bxstar = work7 ;
  vxstar = work5 ;

  LOOPHPC{ // needed since emf needed here, and ok since dq[N2+1] not accessed which isn't created above
    ftemp=z2e_1(v[2][3],k,j,i);
    sgn_va = copysign(1.,ftemp) ;
    ftemp=fabs(ftemp)*dt*OARC23(k,j,i)*ODX(2,3,k);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i];
#else
    D1p=ftemp*z2e_1(osqrtrho,km1,j,i);
    D1m=ftemp*z2e_1(osqrtrho,k,j,i);
#endif

    // values at the foot of the plus characteristic
    //  are, by convention, in im1 zone
    D2 = dt*z2e_1(v[1][3],k,j,i)*OARC23(k,j,i)*ODX(2,3,k) ;
    dist=dx[2][3][k]*G3(1,i)*G4(2,j);


#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    // C+
    Dtota=D2+D1a;
    if(Dtota>0.) { // C+ from -k region
      vp = v[1][1][km1][j][i] + (1.0 - Dtota )*dqv[1][km1][j][i]*dist ;
      bp = v[2][1][km1][j][i] + (1.0 - Dtota )*dqb[1][km1][j][i]*dist ;
    }
    else{ // C+ from +k region
      vp = v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist ;
      bp = v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist ;
    }
    // C-
    Dtota=D2-D1a;
    if(Dtota>0.) { // C- from -k region
      vm = v[1][1][km1][j][i] + (1.0 - Dtota )*dqv[1][km1][j][i]*dist ;
      bm = v[2][1][km1][j][i] + (1.0 - Dtota )*dqb[1][km1][j][i]*dist ;
    }
    else{ // C- from +k region
      vm = v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist ;
    }
#else
    // C+
    Dtota=D2+D1a;
    if( Dtota > DTOTLIMIT) { // C+ from -k region
      vp = v[1][1][km1][j][i] + (1.0 - Dtota )*dqv[1][km1][j][i]*dist ;
      bp = v[2][1][km1][j][i] + (1.0 - Dtota )*dqb[1][km1][j][i]*dist ;
    }
    else if( Dtota < -DTOTLIMIT) { // C+ from +k region
      vp = v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist ;
      bp = v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      vp = 0.5*(v[1][1][km1][j][i] + (1.0 - Dtota )*dqv[1][km1][j][i]*dist+
                v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist);
      bp = 0.5*(v[2][1][km1][j][i] + (1.0 - Dtota )*dqb[1][km1][j][i]*dist+
                v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist);
#else
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][1][km1][j][i] + (1.0 - DTOTLIMIT )*dqv[1][km1][j][i]*dist ;
      dqplus=v[1][1][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[1][k][j][i]*dist ;
      vp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][1][km1][j][i] + (1.0 - DTOTLIMIT )*dqb[1][km1][j][i]*dist ;
      dqplus=v[2][1][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[1][k][j][i]*dist ;
      bp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }
    // C-
    Dtota=D2-D1a;
    if(Dtota > DTOTLIMIT) { // C- from -k region
      vm = v[1][1][km1][j][i] + (1.0 - Dtota )*dqv[1][km1][j][i]*dist ;
      bm = v[2][1][km1][j][i] + (1.0 - Dtota )*dqb[1][km1][j][i]*dist ;
    }
    else if(Dtota < -DTOTLIMIT) { // C- from +k region
      vm = v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      vm = 0.5*(v[1][1][km1][j][i] + (1.0 - Dtota )*dqv[1][km1][j][i]*dist+
                v[1][1][k][j][i] + (-1.0 - Dtota )*dqv[1][k][j][i]*dist);
      bm = 0.5*(v[2][1][km1][j][i] + (1.0 - Dtota )*dqb[1][km1][j][i]*dist+
                v[2][1][k][j][i] + (-1.0 - Dtota )*dqb[1][k][j][i]*dist);
#else
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][1][km1][j][i] + (1.0 - DTOTLIMIT )*dqv[1][km1][j][i]*dist ;
      dqplus=v[1][1][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[1][k][j][i]*dist ;
      vm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][1][km1][j][i] + (1.0 - DTOTLIMIT )*dqb[1][km1][j][i]*dist ;
      dqplus=v[2][1][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[1][k][j][i]*dist ;
      bm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
#if(SYMFORCEMAG==0)
    // C+
    if( (Dtotp=D2+D1p) > 0.) { // C+ from -k region
      vp = v[1][1][km1][j][i] + (1.0 - Dtotp )*dqv[1][km1][j][i]*dist ;
      bp = v[2][1][km1][j][i] + (1.0 - Dtotp )*dqb[1][km1][j][i]*dist ;
    }
    else{ // C+ from +k region
      Dtotm=D2+D1m;
      vp = v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist ;
      bp = v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist ;
    }
    // C-
    if( (Dtotp=D2-D1p) > 0.) { // C- from -k region
      vm = v[1][1][km1][j][i] + (1.0 - Dtotp )*dqv[1][km1][j][i]*dist ;
      bm = v[2][1][km1][j][i] + (1.0 - Dtotp )*dqb[1][km1][j][i]*dist ;
    }
    else{ // C- from +k region
      Dtotm=D2-D1m;
      vm = v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist ;
    }
#else
    // C+
    if( (Dtotp=D2+D1p) > DTOTLIMIT) { // C+ from -k region
      vp = v[1][1][km1][j][i] + (1.0 - Dtotp )*dqv[1][km1][j][i]*dist ;
      bp = v[2][1][km1][j][i] + (1.0 - Dtotp )*dqb[1][km1][j][i]*dist ;
    }
    else if( (Dtotm=D2+D1m) < -DTOTLIMIT) { // C+ from +k region
      vp = v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist ;
      bp = v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist ;
    }
    else{
      vp = 0.5*(v[1][1][km1][j][i] + (1.0 - Dtotp )*dqv[1][km1][j][i]*dist+
		v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist);
      bp = 0.5*(v[2][1][km1][j][i] + (1.0 - Dtotp )*dqb[1][km1][j][i]*dist+
		v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist);
    }
    // C-
    if( (Dtotp=D2-D1p) > DTOTLIMIT) { // C- from -k region
      vm = v[1][1][km1][j][i] + (1.0 - Dtotp )*dqv[1][km1][j][i]*dist ;
      bm = v[2][1][km1][j][i] + (1.0 - Dtotp )*dqb[1][km1][j][i]*dist ;
    }
    else if( (Dtotm=D2-D1m) < -DTOTLIMIT) { // C- from +k region
      vm = v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist ;
    }
    else{
      vm = 0.5*(v[1][1][km1][j][i] + (1.0 - Dtotp )*dqv[1][km1][j][i]*dist+
		v[1][1][k][j][i] + (-1.0 - Dtotm )*dqv[1][k][j][i]*dist);
      bm = 0.5*(v[2][1][km1][j][i] + (1.0 - Dtotp )*dqb[1][km1][j][i]*dist+
		v[2][1][k][j][i] + (-1.0 - Dtotm )*dqb[1][k][j][i]*dist);
    }
#endif

#endif    // end if rhointerp==1

    // bx is really bz and vx is really vz here
    // solution to Stone & Norman eqtn 43,44
#if(RHOINTERP==0)
    //-- use constant rho
    bxstar[k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
    vxstar[k][j][i] = 0.5*(vm + vp + sgn_va*(bm - bp)*srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    srhom=z2e_1(osqrtrho,k,j,i); // C- char rho value
    srhop=z2e_1(osqrtrho,km1,j,i); // C+ char rho value
    bxstar[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
    vxstar[k][j][i] = 0.5*(vm+vp+sgn_va*bxstar[k][j][i]*(srhop-srhom) +sgn_va*(-bp*srhop+bm*srhom) );
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
      magemf[2][k][j][i] = vzstar[k][j][i]*bxstar[k][j][i] - vxstar[k][j][i]*bzstar[k][j][i] - jcurrent[2][k][j][i]*z2c_2(nu_res_real,k,j,i) ;
    }    
  }
  else{
    LOOPHPC{ // needed for below 2 loops.
      magemf[2][k][j][i] = vzstar[k][j][i]*bxstar[k][j][i] - vxstar[k][j][i]*bzstar[k][j][i] ;
    }
  }
  }// end if not true that both N1==1 and N3==1  



  /////////////////////
  /////////////////
  /////////////////  emf1
  //////////////
  //////////////  
  //////////////////////////////


  if(!( (N2==1)&&(N3==1) ) ){ // otherwise no emf1 or b2prim[3] or b2prim[2] needed

    if(RHOINTERP==0){
      // used twice
      LOOPHPC{ // see above
	// needed since for rhointerp==1, need rho at i=-1..N1, j=-1..N2, so find osqrtrho on LOOPFMHPC!
	srhoa[k][j][i] = z2c_1(osqrtrho,k,j,i) ; //1/sqrt(rho)
      }
    }
  
  // sweep in y-direction (on 3-components)
  // first get slopes 
  // dq's act like scalars on rotation
  dqvy_calc(3,v[1],dqv) ;
  dqvy_calc(3,v[2],dqb) ;


  vzstar = work6 ;
  bzstar = work8 ;

  LOOPHPC{ // needed since emf needed here, and ok since dq[N2+1] not accessed which isn't created above
    ftemp=z2e_3(v[2][2],k,j,i);
    sgn_va = copysign(1.,ftemp) ;
    ftemp=fabs(ftemp)*dt*OARC12(k,j,i);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i];
#else
    D1p=ftemp*z2e_3(osqrtrho,k,jm1,i);
    D1m=ftemp*z2e_3(osqrtrho,k,j,i);
#endif

    // values at the foot of the plus characteristic
    //  are, by convention, in im1 zone
    D2 = dt*z2e_3(v[1][2],k,j,i)*OARC12(k,j,i) ;
    dist=dx[2][2][j]*G2(2,i);


#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    // C+
    Dtota=D2+D1a;
    if(Dtota>0.) { // C+ from -j region
      vp = v[1][3][k][jm1][i] + (1.0 - Dtota )*dqv[3][k][jm1][i]*dist ;
      bp = v[2][3][k][jm1][i] + (1.0 - Dtota )*dqb[3][k][jm1][i]*dist ;
    }
    else{ // C+ from +j region
      vp = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }
    // C-
    Dtota=D2-D1a;
    if(Dtota>0.) { // C- from -j region
      vm = v[1][3][k][jm1][i] + (1.0 - Dtota )*dqv[3][k][jm1][i]*dist ;
      bm = v[2][3][k][jm1][i] + (1.0 - Dtota )*dqb[3][k][jm1][i]*dist ;
    }
    else{ // C- from +j region
      vm = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }
#else
    // C+
    Dtota=D2+D1a;
    if( Dtota > DTOTLIMIT) { // C+ from -j region
      vp = v[1][3][k][jm1][i] + (1.0 - Dtota )*dqv[3][k][jm1][i]*dist ;
      bp = v[2][3][k][jm1][i] + (1.0 - Dtota )*dqb[3][k][jm1][i]*dist ;
    }
    else if( Dtota < -DTOTLIMIT) { // C+ from +j region
      vp = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      vp = 0.5*(v[1][3][k][jm1][i] + (1.0 - Dtota )*dqv[3][k][jm1][i]*dist+
                v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist);
      bp = 0.5*(v[2][3][k][jm1][i] + (1.0 - Dtota )*dqb[3][k][jm1][i]*dist+
                v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist);
#else
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][3][k][jm1][i] + (1.0 - DTOTLIMIT )*dqv[3][k][jm1][i]*dist ;
      dqplus=v[1][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[3][k][j][i]*dist ;
      vp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][3][k][jm1][i] + (1.0 - DTOTLIMIT )*dqb[3][k][jm1][i]*dist ;
      dqplus=v[2][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[3][k][j][i]*dist ;
      bp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }
    // C-
    Dtota=D2-D1a;
    if(Dtota > DTOTLIMIT) { // C- from -j region
      vm = v[1][3][k][jm1][i] + (1.0 - Dtota )*dqv[3][k][jm1][i]*dist ;
      bm = v[2][3][k][jm1][i] + (1.0 - Dtota )*dqb[3][k][jm1][i]*dist ;
    }
    else if(Dtota < -DTOTLIMIT) { // C- from +j region
      vm = v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      vm = 0.5*(v[1][3][k][jm1][i] + (1.0 - Dtota )*dqv[3][k][jm1][i]*dist+
                v[1][3][k][j][i] + (-1.0 - Dtota )*dqv[3][k][j][i]*dist);
      bm = 0.5*(v[2][3][k][jm1][i] + (1.0 - Dtota )*dqb[3][k][jm1][i]*dist+
                v[2][3][k][j][i] + (-1.0 - Dtota )*dqb[3][k][j][i]*dist);
#else
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][3][k][jm1][i] + (1.0 - DTOTLIMIT )*dqv[3][k][jm1][i]*dist ;
      dqplus=v[1][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[3][k][j][i]*dist ;
      vm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][3][k][jm1][i] + (1.0 - DTOTLIMIT )*dqb[3][k][jm1][i]*dist ;
      dqplus=v[2][3][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[3][k][j][i]*dist ;
      bm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
#if(SYMFORCEMAG==0)
    // C+
    if( (Dtotp=D2+D1p) > 0.) { // C+ from -j region
      vp = v[1][3][k][jm1][i] + (1.0 - Dtotp )*dqv[3][k][jm1][i]*dist ;
      bp = v[2][3][k][jm1][i] + (1.0 - Dtotp )*dqb[3][k][jm1][i]*dist ;
    }
    else{ // C+ from +j region
      Dtotm=D2+D1m;
      vp = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }
    // C-
    if( (Dtotp=D2-D1p) > 0.) { // C- from -j region
      vm = v[1][3][k][jm1][i] + (1.0 - Dtotp )*dqv[3][k][jm1][i]*dist ;
      bm = v[2][3][k][jm1][i] + (1.0 - Dtotp )*dqb[3][k][jm1][i]*dist ;
    }
    else{ // C- from +j region
      Dtotm=D2-D1m;
      vm = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }
#else
    // C+
    if( (Dtotp=D2+D1p) > DTOTLIMIT) { // C+ from -j region
      vp = v[1][3][k][jm1][i] + (1.0 - Dtotp )*dqv[3][k][jm1][i]*dist ;
      bp = v[2][3][k][jm1][i] + (1.0 - Dtotp )*dqb[3][k][jm1][i]*dist ;
    }
    else if( (Dtotm=D2+D1m) < -DTOTLIMIT) { // C+ from +j region
      vp = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bp = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }
    else{
      vp = 0.5*(v[1][3][k][jm1][i] + (1.0 - Dtotp )*dqv[3][k][jm1][i]*dist+
		v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist);
      bp = 0.5*(v[2][3][k][jm1][i] + (1.0 - Dtotp )*dqb[3][k][jm1][i]*dist+
		v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist);
    }
    // C-
    if( (Dtotp=D2-D1p) > DTOTLIMIT) { // C- from -j region
      vm = v[1][3][k][jm1][i] + (1.0 - Dtotp )*dqv[3][k][jm1][i]*dist ;
      bm = v[2][3][k][jm1][i] + (1.0 - Dtotp )*dqb[3][k][jm1][i]*dist ;
    }
    else if( (Dtotm=D2-D1m) < -DTOTLIMIT) { // C- from +j region
      vm = v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist ;
    }
    else{
      vm = 0.5*(v[1][3][k][jm1][i] + (1.0 - Dtotp )*dqv[3][k][jm1][i]*dist+
		v[1][3][k][j][i] + (-1.0 - Dtotm )*dqv[3][k][j][i]*dist);
      bm = 0.5*(v[2][3][k][jm1][i] + (1.0 - Dtotp )*dqb[3][k][jm1][i]*dist+
		v[2][3][k][j][i] + (-1.0 - Dtotm )*dqb[3][k][j][i]*dist);
    }
#endif

#endif    // end if rhointerp==1
    // solution to Stone & Norman eqtn 43,44
#if(RHOINTERP==0)
    //-- use constant rho
    bzstar[k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
    vzstar[k][j][i] = 0.5*(vm + vp + sgn_va*(bm - bp)*srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    srhom=z2e_3(osqrtrho,k,j,i); // C- char rho value
    srhop=z2e_3(osqrtrho,k,jm1,i); // C+ char rho value
    bzstar[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
    vzstar[k][j][i] = 0.5*(vm+vp+sgn_va*bzstar[k][j][i]*(srhop-srhom) +sgn_va*(-bp*srhop+bm*srhom) );
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif




  }
  // now vx/bxstar live on corners on full grid, act like vectors on rotation


  // sweep in z-direction (on y-components)
  // first get slopes 
  // dq's act like scalars on rotation
  dqvz_calc(2,v[1],dqv) ;
  dqvz_calc(2,v[2],dqb) ;



  bystar = work7 ;
  vystar = work5 ;

  LOOPHPC{ // needed since emf needed here, and ok since dq[N2+1] not accessed which isn't created above
    ftemp=z2e_2(v[2][3],k,j,i);
    sgn_va = copysign(1.,ftemp) ;
    ftemp=fabs(ftemp)*dt*OARC33(k,j,i)*ODX(2,3,k);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i];
#else
    D1p=ftemp*z2e_2(osqrtrho,km1,j,i);
    D1m=ftemp*z2e_2(osqrtrho,k,j,i);
#endif

    // values at the foot of the plus characteristic
    //  are, by convention, in im1 zone
    D2 = dt*z2e_2(v[1][3],k,j,i)*OARC33(k,j,i)*ODX(2,3,k) ;
    dist=dx[2][3][k]*G3(2,i)*G4(1,j);


#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    // C+
    Dtota=D2+D1a;
    if(Dtota>0.) { // C+ from -k region
      vp = v[1][2][km1][j][i] + (1.0 - Dtota )*dqv[2][km1][j][i]*dist ;
      bp = v[2][2][km1][j][i] + (1.0 - Dtota )*dqb[2][km1][j][i]*dist ;
    }
    else{ // C+ from +k region
      vp = v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist ;
      bp = v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist ;
    }
    // C-
    Dtota=D2-D1a;
    if(Dtota>0.) { // C- from -k region
      vm = v[1][2][km1][j][i] + (1.0 - Dtota )*dqv[2][km1][j][i]*dist ;
      bm = v[2][2][km1][j][i] + (1.0 - Dtota )*dqb[2][km1][j][i]*dist ;
    }
    else{ // C- from +k region
      vm = v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist ;
    }
#else
    // C+
    Dtota=D2+D1a;
    if( Dtota > DTOTLIMIT) { // C+ from -k region
      vp = v[1][2][km1][j][i] + (1.0 - Dtota )*dqv[2][km1][j][i]*dist ;
      bp = v[2][2][km1][j][i] + (1.0 - Dtota )*dqb[2][km1][j][i]*dist ;
    }
    else if( Dtota < -DTOTLIMIT) { // C+ from +k region
      vp = v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist ;
      bp = v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      vp = 0.5*(v[1][2][km1][j][i] + (1.0 - Dtota )*dqv[2][km1][j][i]*dist+
                v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist);
      bp = 0.5*(v[2][2][km1][j][i] + (1.0 - Dtota )*dqb[2][km1][j][i]*dist+
                v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist);
#else
      // use symmetrized Dtot
      //      Dtot=D2+D1a;
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][2][km1][j][i] + (1.0 - DTOTLIMIT )*dqv[2][km1][j][i]*dist ;
      dqplus=v[1][2][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[2][k][j][i]*dist ;
      vp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][2][km1][j][i] + (1.0 - DTOTLIMIT )*dqb[2][km1][j][i]*dist ;
      dqplus=v[2][2][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[2][k][j][i]*dist ;
      bp=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }
    // C-
    Dtota=D2-D1a;
    if(Dtota > DTOTLIMIT) { // C- from -k region
      vm = v[1][2][km1][j][i] + (1.0 - Dtota )*dqv[2][km1][j][i]*dist ;
      bm = v[2][2][km1][j][i] + (1.0 - Dtota )*dqb[2][km1][j][i]*dist ;
    }
    else if(Dtota < -DTOTLIMIT) { // C- from +k region
      vm = v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist ;
    }
    else{
#if(SYMGOOD==0)
      vm = 0.5*(v[1][2][km1][j][i] + (1.0 - Dtota )*dqv[2][km1][j][i]*dist+
                v[1][2][k][j][i] + (-1.0 - Dtota )*dqv[2][k][j][i]*dist);
      bm = 0.5*(v[2][2][km1][j][i] + (1.0 - Dtota )*dqb[2][km1][j][i]*dist+
                v[2][2][k][j][i] + (-1.0 - Dtota )*dqb[2][k][j][i]*dist);
#else
      // like sweep.c strong sym:
      fraczone=Dtota*INVDTOTLIMIT;
      dqminus=v[1][2][km1][j][i] + (1.0 - DTOTLIMIT )*dqv[2][km1][j][i]*dist ;
      dqplus=v[1][2][k][j][i] + (-1.0 + DTOTLIMIT )*dqv[2][k][j][i]*dist ;
      vm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));

      dqminus=v[2][2][km1][j][i] + (1.0 - DTOTLIMIT )*dqb[2][km1][j][i]*dist ;
      dqplus=v[2][2][k][j][i] + (-1.0 + DTOTLIMIT )*dqb[2][k][j][i]*dist ;
      bm=0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
#if(SYMFORCEMAG==0)
    // C+
    if( (Dtotp=D2+D1p) > 0.) { // C+ from -k region
      vp = v[1][2][km1][j][i] + (1.0 - Dtotp )*dqv[2][km1][j][i]*dist ;
      bp = v[2][2][km1][j][i] + (1.0 - Dtotp )*dqb[2][km1][j][i]*dist ;
    }
    else{ // C+ from +k region
      Dtotm=D2+D1m;
      vp = v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist ;
      bp = v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist ;
    }
    // C-
    if( (Dtotp=D2-D1p) > 0.) { // C- from -k region
      vm = v[1][2][km1][j][i] + (1.0 - Dtotp )*dqv[2][km1][j][i]*dist ;
      bm = v[2][2][km1][j][i] + (1.0 - Dtotp )*dqb[2][km1][j][i]*dist ;
    }
    else{ // C- from +k region
      Dtotm=D2-D1m;
      vm = v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist ;
    }
#else
    // C+
    if( (Dtotp=D2+D1p) > DTOTLIMIT) { // C+ from -k region
      vp = v[1][2][km1][j][i] + (1.0 - Dtotp )*dqv[2][km1][j][i]*dist ;
      bp = v[2][2][km1][j][i] + (1.0 - Dtotp )*dqb[2][km1][j][i]*dist ;
    }
    else if( (Dtotm=D2+D1m) < -DTOTLIMIT) { // C+ from +k region
      vp = v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist ;
      bp = v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist ;
    }
    else{
      vp = 0.5*(v[1][2][km1][j][i] + (1.0 - Dtotp )*dqv[2][km1][j][i]*dist+
		v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist);
      bp = 0.5*(v[2][2][km1][j][i] + (1.0 - Dtotp )*dqb[2][km1][j][i]*dist+
		v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist);
    }
    // C-
    if( (Dtotp=D2-D1p) > DTOTLIMIT) { // C- from -k region
      vm = v[1][2][km1][j][i] + (1.0 - Dtotp )*dqv[2][km1][j][i]*dist ;
      bm = v[2][2][km1][j][i] + (1.0 - Dtotp )*dqb[2][km1][j][i]*dist ;
    }
    else if( (Dtotm=D2-D1m) < -DTOTLIMIT) { // C- from +k region
      vm = v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist ;
    }
    else{
      vm = 0.5*(v[1][2][km1][j][i] + (1.0 - Dtotp )*dqv[2][km1][j][i]*dist+
		v[1][2][k][j][i] + (-1.0 - Dtotm )*dqv[2][k][j][i]*dist);
      bm = 0.5*(v[2][2][km1][j][i] + (1.0 - Dtotp )*dqb[2][km1][j][i]*dist+
		v[2][2][k][j][i] + (-1.0 - Dtotm )*dqb[2][k][j][i]*dist);
    }
#endif

#endif    // end if rhointerp==1
    // solution to Stone & Norman eqtn 43,44
#if(RHOINTERP==0)
    //-- use constant rho
    bystar[k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
    vystar[k][j][i] = 0.5*(vm + vp + sgn_va*(bm - bp)*srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    srhom=z2e_2(osqrtrho,k,j,i); // C- char rho value
    srhop=z2e_2(osqrtrho,km1,j,i); // C+ char rho value
    bystar[k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
    vystar[k][j][i] = 0.5*(vm+vp+sgn_va*bystar[k][j][i]*(srhop-srhom) +sgn_va*(-bp*srhop+bm*srhom) );
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif




  }
  // now yx/bystar live on corners on full grid, act like vectors on rotation

  // calculate emf=v cross B

  // emf lives on corners on full comp grid including border
  if(RESMEM&&(res_real==1)){
    LOOPHPC{ // needed for below 2 loops.
      magemf[1][k][j][i] =
	vystar[k][j][i]*bzstar[k][j][i] 
	- vzstar[k][j][i]*bystar[k][j][i] 
	- jcurrent[1][k][j][i]*z2c_1(nu_res_real,k,j,i) ;
    }    
  }
  else{
    LOOPHPC{ // needed for below 2 loops.
      magemf[1][k][j][i] = vystar[k][j][i]*bzstar[k][j][i] - vzstar[k][j][i]*bystar[k][j][i] ;
    }
  }
  }// end if not true that both N2==1 and N3==1



  if(BOUNDFIELD==2){
    // BOUND EMFS SO DIVB=0 CONSERVED EXACTLY
    bound(NULL,magemf,0,-5,123);
  }



  // STONE AND MILLER
  // update field, CT with emf as conserved flux, hence divB=0.  eq 11,12 for B, eq 30,31 for v



  // must do mag field second since above v's use old B (like new Zeus says.  Old zeus says do field first)
  // newest zeus says split lorentz forces all before emf calcs, unlike in this _v1 function

  if(VOLUMEDIFF==0){ // old differential way

  if(!( (N2==1)&&(N3==1) ) ){
    if(mocctbx1){
      LOOPB1{ // needs emf on i=0..N1-1 and j=0..N2
	v[2][1][k][j][i] += dt*(
#if(N2!=1)
				OARC42(k,j,i)/G4(2,j)*(G4(1,jp1)*magemf[3][k][jp1][i] - G4(1,j)*magemf[3][k][j][i])
#else
				0
#endif
#if(N3!=1)
				+OARC23(k,j,i)*ODX(1,3,k)*(magemf[2][k][j][i]-magemf[2][kp1][j][i])
#endif
				) ;
      }
    }
  }
  if(!( (N1==1)&&(N3==1) ) ){
    if(mocctbx2){
      LOOPB2{ // needs emf at i=0..N1 and j=0..N2-1+periodicx2special (LOOPH above is fine)
	// 
	v[2][2][k][j][i] += dt*(
#if(N1!=1)
				OARC41(k,j,i)/G3(2,i)*(G3(1,i)*magemf[3][k][j][i] - G3(1,ip1)*magemf[3][k][j][ip1])
#else
				0
#endif
#if(N3!=1)
				+OARC33(k,j,i)*ODX(1,3,k)*(magemf[1][kp1][j][i]-magemf[1][k][j][i])
#endif
				) ;
	
      }
    }
  }
  if(!( (N1==1)&&(N2==1) )){
    if(mocctbx3){
      LOOPB3{ // needs emf at i=0..N1 and j=0..N2-1
	v[2][3][k][j][i] += dt*(
#if(N1!=1)
				OARC41(k,j,i)/G2(2,i)*(G2(1,ip1)*magemf[2][k][j][ip1] - G2(1,i)*magemf[2][k][j][i])
#else
				0
#endif
#if(N2!=1)
				+OARC32(k,j,i)*(magemf[1][k][j][i]-magemf[1][k][jp1][i])
#endif
				) ;
	
      }
    }
  }

  }
  else{ // volume way, which still preserves divB=0

  if(!( (N2==1)&&(N3==1) ) ){
    if(mocctbx1){
      LOOPB1{ // needs emf on i=0..N1-1 and j=0..N2
	v[2][1][k][j][i] += dt*(
#if(N2!=1)
				ODVL(1,2,j)/G2(1,i)*(G4(1,jp1)*magemf[3][k][jp1][i] - G4(1,j)*magemf[3][k][j][i])
#else
				0
#endif
#if(N3!=1) // no change for volume way
				+OARC23(k,j,i)*ODX(1,3,k)*(magemf[2][k][j][i]-magemf[2][kp1][j][i])
#endif
				) ;
      }
    }
  }
  if(!( (N1==1)&&(N3==1) ) ){
    if(mocctbx2){
      LOOPB2{ // needs emf at i=0..N1 and j=0..N2-1+periodicx2special (LOOPH above is fine)
	v[2][2][k][j][i] += dt*(
#if(N1!=1)
				G2(2,i)*ODVL(1,1,i)*(G3(1,i)*magemf[3][k][j][i] - G3(1,ip1)*magemf[3][k][j][ip1])
#else
				0
#endif
#if(N3!=1) // no change for volume way
				+OARC33(k,j,i)*ODX(1,3,k)*(magemf[1][kp1][j][i]-magemf[1][k][j][i])
#endif
				) ;
	
      }
    }
  }
  if(!( (N1==1)&&(N2==1) )){
    if(mocctbx3){
      LOOPB3{ // needs emf at i=0..N1 and j=0..N2-1
	v[2][3][k][j][i] += dt*(
#if(N1!=1)
				G3(2,i)*ODVL(1,1,i)*(G2(1,ip1)*magemf[2][k][j][ip1] - G2(1,i)*magemf[2][k][j][i])
#else
				0
#endif
#if(N2!=1)
				+G4(2,j)*ODVL(1,2,j)/G2(2,i)*(magemf[1][k][j][i]-magemf[1][k][jp1][i])
#endif
				) ;
	
      }
    }
  }


  } // end if volume way


  if((BOUNDTYPE>1)&&(BOUNDFIELD==1)){
  // limit what's bounded to only necessary changing items for faster reduced dimensional speeds(only relevant in 1D)
  // note that divB=0 correction in bound() is fine since other components don't change anyways(so didn't need to be updated, so are ready for divB=0 correction already)
  if(!( (N2==1)&&(N3==1) ) ){ // then do 1-comp
    if(!( (N1==1)&&(N3==1) ) ){ // then do 2-comp
      if(!( (N1==1)&&(N2==1) )){ // then do 3-comp
	bound(NULL,NULL,0,-1,123);
      }
      else{
	bound(NULL,NULL,0,-1,12); // no 3
      }
    }
    else{ // no 2
      if(!( (N1==1)&&(N2==1) )){
	bound(NULL,NULL,0,-1,13);
      }
      else{
	bound(NULL,NULL,0,-1,1); // no 3
      }
    }
  }
  else{// no 1
    if(!( (N1==1)&&(N3==1) ) ){ // then do 2-comp
      if(!( (N1==1)&&(N2==1) )){ // then do 3-comp
	bound(NULL,NULL,0,-1,23);
      }
      else{
	bound(NULL,NULL,0,-1,2); // no 3
      }
    }
    else{ // no 2
      if(!( (N1==1)&&(N2==1) )){
	bound(NULL,NULL,0,-1,3);
      }
      else{
	// never get here // no 3
      }
    }
  }
  }
  else{ // then just normal
  // limit what's bounded to only necessary changing items for faster reduced dimensional speeds(only relevant in 1D)
  // note that divB=0 correction in bound() is fine since other components don't change anyways(so didn't need to be updated, so are ready for divB=0 correction already)
  if(!( (N2==1)&&(N3==1) ) ){ // then do 1-comp
    if(!( (N1==1)&&(N3==1) ) ){ // then do 2-comp
      if(!( (N1==1)&&(N2==1) )){ // then do 3-comp
	bound(NULL,NULL,0,-1,123);
      }
      else{
	bound(NULL,NULL,0,-1,12); // no 3
      }
    }
    else{ // no 2
      if(!( (N1==1)&&(N2==1) )){
	bound(NULL,NULL,0,-1,13);
      }
      else{
	bound(NULL,NULL,0,-1,1); // no 3
      }
    }
  }
  else{// no 1
    if(!( (N1==1)&&(N3==1) ) ){ // then do 2-comp
      if(!( (N1==1)&&(N2==1) )){ // then do 3-comp
	bound(NULL,NULL,0,-1,23);
      }
      else{
	bound(NULL,NULL,0,-1,2); // no 3
      }
    }
    else{ // no 2
      if(!( (N1==1)&&(N2==1) )){
	bound(NULL,NULL,0,-1,3);
      }
      else{
	// never get here // no 3
      }
    }
  }
  }
}



// 1d/2d/3d valid
// this is just cut-pasted out of moc_ct_3d_v1()
// cut-paste then delete all star&emf*"capital" D references(e.g. D2) references, field calc, and use more strict loop control since no coupled dependencies
void lorentz_3d(void)
{
  FTYPE fraczone,dqminus,dqplus;
  static FTYPE (*dqv)[N3M][N2M][N1M],(*dqb)[N3M][N2M][N1M];
  static FTYPE (*b1prim)[N3M][N2M][N1M],(*b2prim)[N3M][N2M][N1M];
  static FTYPE bxa,bya,bza ;
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



  if((FLOATTYPE==0)&&(SYMFORCEMAG==1)){
    fprintf(fail_file,"moc_ct: floats and symmetry are known to fail--no idea why yet\n");
    myexit(1);
  }

  dqv = workv1 ;
  dqb = workv2 ;
  b1prim=workv4;
  b2prim=workv5;

  rhoa=work3;
  srhoa=work4;
  // work5,6,7,8,9,10 used below


  //////////////////////////////
  //////////
  ///////////    b1prim[2] and b1prim[1] : like emf3
  ////////////
  ///////////  
  ////////////////////////////

  if(! ( (N1==1)&&(N2==1) ) ){ // otherwise no need for emf3 or b1prim[1] or b1prim[2]

    if(RHOINTERP==0){
      // used twice
      LOOPHPC{ // see above
	// needed since for rhointerp==1, need rho at i=-1..N1, j=-1..N2, so find osqrtrho on LOOPFMHPC!
	srhoa[k][j][i] = z2c_3(osqrtrho,k,j,i) ; //1/sqrt(rho)
      }
    }

  if(N1>1){ // otherwise no need for b1prim[2]

  // sweep in x-direction 
  // first get slopes 
  dqvx_calc(2,v[1],dqv) ;
  dqvx_calc(2,v[2],dqb) ;
    
  LOOPHPC{  // needed since emf needed here, and ok since dq[N1+1] not accessed which isn't created above
    ftemp=z2e_2(v[2][1],k,j,i);
    sgn_va=copysign(1.0,ftemp);
    ftemp=fabs(ftemp)*dt*OARC31(k,j,i);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i] ;
#else
    D1p=ftemp*z2e_2(osqrtrho,k,j,im1);
    D1m=ftemp*z2e_2(osqrtrho,k,j,i);
#endif

    // determine which formula is for C+ and which is C-
    dist=dx[2][1][i];

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
    b1prim[2][k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    // same srhom/srhop from above
    //srhom=z2e_2(osqrtrho,k,j,i); // C- char rho value
    //srhop=z2e_2(osqrtrho,k,j,im1); // C+ char rho value
    b1prim[2][k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif


  }
  }// end if N1>1
  // now vy/bystar live on corners and fill entire grid and border of comp grid



  if(N2>1){ // otherwise no b1prim[1] needed
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
    D1p=ftemp*z2e_1(osqrtrho,k,jm1,i);
    D1m=ftemp*z2e_1(osqrtrho,k,j,i);
#endif

    // values at the foot of the plus characteristic
    //  are, by convention, in im1 zone
    dist=dx[2][2][j]*G2(1,i);


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
    b1prim[1][k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    //srhom=z2e_1(osqrtrho,k,j,i); // C- char rho value
    //  srhop=z2e_1(osqrtrho,k,jm1,i); // C+ char rho value
    b1prim[1][k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif

  }
  // now vx/bxstar live on corners on full grid, act like vectors on rotation
  }// end if N2>1
  }// end if not true that both N2==1 and N1==1





  //////////////////////////////////
  //////////////////
  /////////////////   find b1prim[3] and b2prim[1]: like emf2
  ////////////////
  ///////////////
  //////////////////////////////


  if(!( (N1==1)&&(N3==1) ) ){

    if(RHOINTERP==0){
      // used twice
      LOOPHPC{ // see above
	// needed since for rhointerp==1, need rho at i=-1..N1, j=-1..N2, so find osqrtrho on LOOPFMHPC!
	srhoa[k][j][i] = z2c_2(osqrtrho,k,j,i) ; //1/sqrt(rho)
      }
    }

  if(N1>1){
  // sweep in x-direction (on 3-components)
  // first get slopes 
  dqvx_calc(3,v[1],dqv) ;
  dqvx_calc(3,v[2],dqb) ;



  LOOPHPC{  // needed since emf needed here, and ok since dq[N1+1] not accessed which isn't created above
    ftemp=z2e_3(v[2][1],k,j,i);
    sgn_va=copysign(1.0,ftemp);
    ftemp=fabs(ftemp)*dt*OARC11(k,j,i);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i] ;
#else
    D1p=ftemp*z2e_3(osqrtrho,k,j,im1);
    D1m=ftemp*z2e_3(osqrtrho,k,j,i);
#endif

    // determine which formula is for C+ and which is C-
    dist=dx[2][1][i];

    // find bprim, for updating velocities, no need for vyprime
    // lagrangian frame since op-split and other direction already moving
    // not in source code for more accurate evolution of alphven waves using upwinding

    
#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    // in this case D1p and D1m are equal to D1a
      vp = v[1][3][k][j][im1] + (1.0-D1a)*dqv[3][k][j][im1]*dist ; // C+
      bp = v[2][3][k][j][im1] + (1.0-D1a)*dqb[3][k][j][im1]*dist ; // C+
      vm = v[1][3][k][j][i] + (-1.0+D1a)*dqv[3][k][j][i]*dist ; // C-
      bm = v[2][3][k][j][i] + (-1.0+D1a)*dqb[3][k][j][i]*dist ; // C- 
#else
    // in this case D1p and D1m are equal to D1a
    if(D1a>DTOTLIMIT){ // keep p m names so like above
      vp = v[1][3][k][j][im1] + (1.0-D1a)*dqv[3][k][j][im1]*dist ; // C+
      bp = v[2][3][k][j][im1] + (1.0-D1a)*dqb[3][k][j][im1]*dist ; // C+
      vm = v[1][3][k][j][i] + (-1.0+D1a)*dqv[3][k][j][i]*dist ; // C-
      bm = v[2][3][k][j][i] + (-1.0+D1a)*dqb[3][k][j][i]*dist ; // C- 
    }
    else{
#if(SYMFUNKEDM1==0)
      // true, but only need difference in vp and vm for byprim
      //      vp = vm = 0.5*(v[1][3][k][j][im1] + (1.0-D1a)*dqv[3][k][j][im1]*dist+
      //		     v[1][3][k][j][i] + (-1.0+D1a)*dqv[3][k][j][i]*dist);
      vm = vp = 0.0; // since only need difference to be 0
      bp = bm = 0.5*(v[2][3][k][j][im1] + (1.0-D1a)*dqb[3][k][j][im1]*dist+
                     v[2][3][k][j][i] + (-1.0+D1a)*dqb[3][k][j][i]*dist);

#else
      // this is a little different than sweep.c
      fraczone=D1a*INVDTOTLIMIT;

      dqminus=v[1][3][k][j][im1] + (1.0-DTOTLIMIT)*dqv[3][k][j][im1]*dist ; 
      dqplus=v[1][3][k][j][i] + (-1.0+DTOTLIMIT)*dqv[3][k][j][i]*dist ;
      vm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      vp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));

      //      vm = vp = 0.0; // since only need difference to be 0

      dqminus=v[2][3][k][j][im1] + (1.0-DTOTLIMIT)*dqb[3][k][j][im1]*dist ;
      dqplus=v[2][3][k][j][i] + (-1.0+DTOTLIMIT)*dqb[3][k][j][i]*dist ;
      bm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      bp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
#endif
      //vp = v[1][3][k][j][im1] + (1.0-DTOTLIMIT)*dqv[3][k][j][im1]*dist ; // C+
      //vm = v[1][3][k][j][i] + (-1.0+DTOTLIMIT)*dqv[3][k][j][i]*dist ; // C-
    }
#endif

#elif(RHOINTERP==1)
    // no strong sym
#if(SYMFORCEMAG==0)
      vp = v[1][3][k][j][im1] + (1.0-D1p)*dqv[3][k][j][im1]*dist ; // C+
      bp = v[2][3][k][j][im1] + (1.0-D1p)*dqb[3][k][j][im1]*dist ; // C+
      vm = v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist ; // C-
      bm = v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist ; // C- 
#else
    if(D1p>DTOTLIMIT){
      vp = v[1][3][k][j][im1] + (1.0-D1p)*dqv[3][k][j][im1]*dist ; // C+
      bp = v[2][3][k][j][im1] + (1.0-D1p)*dqb[3][k][j][im1]*dist ; // C+
    }
    else{
      vp = 0.5*(v[1][3][k][j][im1] + (1.0-D1p)*dqv[3][k][j][im1]*dist+
		v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist);
      bp = 0.5*(v[2][3][k][j][im1] + (1.0-D1p)*dqb[3][k][j][im1]*dist+
		v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist);
    }
    if(D1m>DTOTLIMIT){
      vm = v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist ; // C-
      bm = v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist ; // C- 
    }
    else{
      vm = 0.5*(v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist+
		v[1][3][k][j][im1] + (1.0-D1p)*dqv[3][k][j][im1]*dist );
      bm = 0.5*(v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist+
		v[2][3][k][j][im1] + (1.0-D1p)*dqb[3][k][j][im1]*dist);
    }
#endif

#endif // end if rhointerp==1

#if(RHOINTERP==0)
    b1prim[3][k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    // same srhom/srhop from above
    //srhom=z2e_3(osqrtrho,k,j,i); // C- char rho value
    //srhop=z2e_3(osqrtrho,k,j,im1); // C+ char rho value
    b1prim[3][k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif



  }
  }// end if N1>1
  // now vy/bystar live on corners and fill entire grid and border of comp grid



  if(N3>1){
  // sweep in z-direction (on x-components)
  // first get slopes 
  // dq's act like scalars on rotation
  dqvz_calc(1,v[1],dqv) ;
  dqvz_calc(1,v[2],dqb) ;


  LOOPHPC{ // needed since emf needed here, and ok since dq[N2+1] not accessed which isn't created above
    ftemp=z2e_1(v[2][3],k,j,i);
    sgn_va = copysign(1.,ftemp) ;
    ftemp=fabs(ftemp)*dt*OARC23(k,j,i)*ODX(2,3,k);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i];
#else
    D1p=ftemp*z2e_1(osqrtrho,km1,j,i);
    D1m=ftemp*z2e_1(osqrtrho,k,j,i);
#endif

    // values at the foot of the plus characteristic
    //  are, by convention, in im1 zone
    dist=dx[2][3][k]*G3(1,i)*G4(2,j);

    // find bprim, to update velocities, no need for vxprim

#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    vp = v[1][1][km1][j][i] + (1.0-D1a)*dqv[1][km1][j][i]*dist ;
    bp = v[2][1][km1][j][i] + (1.0-D1a)*dqb[1][km1][j][i]*dist ;
    vm = v[1][1][k][j][i] + (-1.0+D1a)*dqv[1][k][j][i]*dist ;
    bm = v[2][1][k][j][i] + (-1.0+D1a)*dqb[1][k][j][i]*dist ;
#else
    // D1p=D1m=D1a
    if(D1a>DTOTLIMIT){
      vp = v[1][1][km1][j][i] + (1.0-D1a)*dqv[1][km1][j][i]*dist ;
      bp = v[2][1][km1][j][i] + (1.0-D1a)*dqb[1][km1][j][i]*dist ;
      vm = v[1][1][k][j][i] + (-1.0+D1a)*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0+D1a)*dqb[1][k][j][i]*dist ;
    }
    else{
#if(SYMFUNKEDM2==0)
      // true, but only need diff (vm-vp) for bxprim
      //      vp = vm = 0.5*(v[1][1][km1][j][i] + (1.0-D1a)*dqv[1][km1][j][i]*dist+
      //                v[1][1][k][j][i] + (-1.0+D1a)*dqv[1][k][j][i]*dist);
      vp=vm=0.0;
      bp = bm = 0.5*(v[2][1][km1][j][i] + (1.0-D1a)*dqb[1][km1][j][i]*dist+
                     v[2][1][k][j][i] + (-1.0+D1a)*dqb[1][k][j][i]*dist);
#else
      // this is a little different than sweep.c
      fraczone=D1a*INVDTOTLIMIT;

      dqminus=v[1][1][km1][j][i] + (1.0-DTOTLIMIT)*dqv[1][km1][j][i]*dist ;
      dqplus=v[1][1][k][j][i] + (-1.0+DTOTLIMIT)*dqv[1][k][j][i]*dist ;
      vm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      vp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));

      dqminus=v[2][1][km1][j][i] + (1.0-DTOTLIMIT)*dqb[1][km1][j][i]*dist ;
      dqplus=v[2][1][k][j][i] + (-1.0+DTOTLIMIT)*dqb[1][k][j][i]*dist ;
      bm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      bp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
    // no strong sym
#if(SYMFORCEMAG==0)
      vp = v[1][1][km1][j][i] + (1.0-D1p)*dqv[1][km1][j][i]*dist ;
      bp = v[2][1][km1][j][i] + (1.0-D1p)*dqb[1][km1][j][i]*dist ;
      vm = v[1][1][k][j][i] + (-1.0+D1m)*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0+D1m)*dqb[1][k][j][i]*dist ;
#else
    if(D1p>DTOTLIMIT){
      vp = v[1][1][km1][j][i] + (1.0-D1p)*dqv[1][km1][j][i]*dist ;
      bp = v[2][1][km1][j][i] + (1.0-D1p)*dqb[1][km1][j][i]*dist ;
    }
    else{
      vp = 0.5*(v[1][1][km1][j][i] + (1.0-D1p)*dqv[1][km1][j][i]*dist+
		v[1][1][k][j][i] + (-1.0+D1m)*dqv[1][k][j][i]*dist);
      bp = 0.5*(v[2][1][km1][j][i] + (1.0-D1p)*dqb[1][km1][j][i]*dist+
		v[2][1][k][j][i] + (-1.0+D1m)*dqb[1][k][j][i]*dist);
    }
    if(D1m>DTOTLIMIT){
      vm = v[1][1][k][j][i] + (-1.0+D1m)*dqv[1][k][j][i]*dist ;
      bm = v[2][1][k][j][i] + (-1.0+D1m)*dqb[1][k][j][i]*dist ;
    }
    else{
      vm = 0.5*(v[1][1][km1][j][i] + (1.0-D1p)*dqv[1][km1][j][i]*dist+
		v[1][1][k][j][i] + (-1.0+D1m)*dqv[1][k][j][i]*dist);
      bm = 0.5*(v[2][1][km1][j][i] + (1.0-D1p)*dqb[1][km1][j][i]*dist+
		v[2][1][k][j][i] + (-1.0+D1m)*dqb[1][k][j][i]*dist);      
    }
#endif

#endif // end if rhointerp==1

#if(RHOINTERP==0)
    b2prim[1][k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    //srhom=z2e_1(osqrtrho,k,j,i); // C- char rho value
    //  srhop=z2e_1(osqrtrho,km1,j,i); // C+ char rho value
    b2prim[1][k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif

  }
  }// end if N3>1
  // now vx/bxstar live on corners on full grid, act like vectors on rotation
  }// end if not true that both N1==1 and N3==1  



  /////////////////////
  /////////////////
  /////////////////  find b2prim[3] and b2prim[2]: like emf1
  //////////////
  //////////////  
  //////////////////////////////


  if(!( (N2==1)&&(N3==1) ) ){ // otherwise no emf1 or b2prim[3] or b2prim[2] needed

    if(RHOINTERP==0){
      // used twice
      LOOPHPC{ // see above
	// needed since for rhointerp==1, need rho at i=-1..N1, j=-1..N2, so find osqrtrho on LOOPFMHPC!
	srhoa[k][j][i] = z2c_1(osqrtrho,k,j,i) ; //1/sqrt(rho)
      }
    }

  if(N2>1){

  // sweep in y-direction (on 3-components)
  // first get slopes 
  // dq's act like scalars on rotation
  dqvy_calc(3,v[1],dqv) ;
  dqvy_calc(3,v[2],dqb) ;



  LOOPHPC{ // needed since emf needed here, and ok since dq[N2+1] not accessed which isn't created above
    ftemp=z2e_3(v[2][2],k,j,i);
    sgn_va = copysign(1.,ftemp) ;
    ftemp=fabs(ftemp)*dt*OARC12(k,j,i);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i];
#else
    D1p=ftemp*z2e_3(osqrtrho,k,jm1,i);
    D1m=ftemp*z2e_3(osqrtrho,k,j,i);
#endif

    // values at the foot of the plus characteristic
    //  are, by convention, in im1 zone
    dist=dx[2][2][j]*G2(2,i);


    // find bprim, to update velocities, no need for vxprim


#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    vp = v[1][3][k][jm1][i] + (1.0-D1a)*dqv[3][k][jm1][i]*dist ;
    bp = v[2][3][k][jm1][i] + (1.0-D1a)*dqb[3][k][jm1][i]*dist ;
    vm = v[1][3][k][j][i] + (-1.0+D1a)*dqv[3][k][j][i]*dist ;
    bm = v[2][3][k][j][i] + (-1.0+D1a)*dqb[3][k][j][i]*dist ;
#else
    // D1p=D1m=D1a
    if(D1a>DTOTLIMIT){
      vp = v[1][3][k][jm1][i] + (1.0-D1a)*dqv[3][k][jm1][i]*dist ;
      bp = v[2][3][k][jm1][i] + (1.0-D1a)*dqb[3][k][jm1][i]*dist ;
      vm = v[1][3][k][j][i] + (-1.0+D1a)*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0+D1a)*dqb[3][k][j][i]*dist ;
    }
    else{
#if(SYMFUNKEDM2==0)
      // true, but only need diff (vm-vp) for bxprim
      //      vp = vm = 0.5*(v[1][3][k][jm1][i] + (1.0-D1a)*dqv[3][k][jm1][i]*dist+
      //                v[1][3][k][j][i] + (-1.0+D1a)*dqv[3][k][j][i]*dist);
      vp=vm=0.0;
      bp = bm = 0.5*(v[2][3][k][jm1][i] + (1.0-D1a)*dqb[3][k][jm1][i]*dist+
                     v[2][3][k][j][i] + (-1.0+D1a)*dqb[3][k][j][i]*dist);
#else
      // this is a little different than sweep.c
      fraczone=D1a*INVDTOTLIMIT;

      dqminus=v[1][3][k][jm1][i] + (1.0-DTOTLIMIT)*dqv[3][k][jm1][i]*dist ;
      dqplus=v[1][3][k][j][i] + (-1.0+DTOTLIMIT)*dqv[3][k][j][i]*dist ;
      vm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      vp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));

      dqminus=v[2][3][k][jm1][i] + (1.0-DTOTLIMIT)*dqb[3][k][jm1][i]*dist ;
      dqplus=v[2][3][k][j][i] + (-1.0+DTOTLIMIT)*dqb[3][k][j][i]*dist ;
      bm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      bp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
    // no strong sym
#if(SYMFORCEMAG==0)
      vp = v[1][3][k][jm1][i] + (1.0-D1p)*dqv[3][k][jm1][i]*dist ;
      bp = v[2][3][k][jm1][i] + (1.0-D1p)*dqb[3][k][jm1][i]*dist ;
      vm = v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist ;
#else
    if(D1p>DTOTLIMIT){
      vp = v[1][3][k][jm1][i] + (1.0-D1p)*dqv[3][k][jm1][i]*dist ;
      bp = v[2][3][k][jm1][i] + (1.0-D1p)*dqb[3][k][jm1][i]*dist ;
    }
    else{
      vp = 0.5*(v[1][3][k][jm1][i] + (1.0-D1p)*dqv[3][k][jm1][i]*dist+
		v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist);
      bp = 0.5*(v[2][3][k][jm1][i] + (1.0-D1p)*dqb[3][k][jm1][i]*dist+
		v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist);
    }
    if(D1m>DTOTLIMIT){
      vm = v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist ;
      bm = v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist ;
    }
    else{
      vm = 0.5*(v[1][3][k][jm1][i] + (1.0-D1p)*dqv[3][k][jm1][i]*dist+
		v[1][3][k][j][i] + (-1.0+D1m)*dqv[3][k][j][i]*dist);
      bm = 0.5*(v[2][3][k][jm1][i] + (1.0-D1p)*dqb[3][k][jm1][i]*dist+
		v[2][3][k][j][i] + (-1.0+D1m)*dqb[3][k][j][i]*dist);      
    }
#endif

#endif // end if rhointerp==1

#if(RHOINTERP==0)
    b2prim[3][k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    //srhom=z2e_3(osqrtrho,k,j,i); // C- char rho value
    //  srhop=z2e_3(osqrtrho,k,jm1,i); // C+ char rho value
    b2prim[3][k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif


  }
  }// end if N2>1
  // now vx/bxstar live on corners on full grid, act like vectors on rotation


  if(N3>1){
  // sweep in z-direction (on y-components)
  // first get slopes 
  // dq's act like scalars on rotation
  dqvz_calc(2,v[1],dqv) ;
  dqvz_calc(2,v[2],dqb) ;



  LOOPHPC{ // needed since emf needed here, and ok since dq[N2+1] not accessed which isn't created above
    ftemp=z2e_2(v[2][3],k,j,i);
    sgn_va = copysign(1.,ftemp) ;
    ftemp=fabs(ftemp)*dt*OARC33(k,j,i)*ODX(2,3,k);
#if(RHOINTERP==0)
    D1a=ftemp*srhoa[k][j][i];
#else
    D1p=ftemp*z2e_2(osqrtrho,km1,j,i);
    D1m=ftemp*z2e_2(osqrtrho,k,j,i);
#endif

    // values at the foot of the plus characteristic
    //  are, by convention, in im1 zone
    dist=dx[2][3][k]*G3(2,i)*G4(1,j);



    // find bprim, to update velocities, no need for vxprim

#if(RHOINTERP==0)
#if(SYMFORCEMAG==0)
    vp = v[1][2][km1][j][i] + (1.0-D1a)*dqv[2][km1][j][i]*dist ;
    bp = v[2][2][km1][j][i] + (1.0-D1a)*dqb[2][km1][j][i]*dist ;
    vm = v[1][2][k][j][i] + (-1.0+D1a)*dqv[2][k][j][i]*dist ;
    bm = v[2][2][k][j][i] + (-1.0+D1a)*dqb[2][k][j][i]*dist ;
#else
    // D1p=D1m=D1a
    if(D1a>DTOTLIMIT){
      vp = v[1][2][km1][j][i] + (1.0-D1a)*dqv[2][km1][j][i]*dist ;
      bp = v[2][2][km1][j][i] + (1.0-D1a)*dqb[2][km1][j][i]*dist ;
      vm = v[1][2][k][j][i] + (-1.0+D1a)*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0+D1a)*dqb[2][k][j][i]*dist ;
    }
    else{
#if(SYMFUNKEDM2==0)
      // true, but only need diff (vm-vp) for bxprim
      //      vp = vm = 0.5*(v[1][2][km1][j][i] + (1.0-D1a)*dqv[2][km1][j][i]*dist+
      //                v[1][2][k][j][i] + (-1.0+D1a)*dqv[2][k][j][i]*dist);
      vp=vm=0.0;
      bp = bm = 0.5*(v[2][2][km1][j][i] + (1.0-D1a)*dqb[2][km1][j][i]*dist+
                     v[2][2][k][j][i] + (-1.0+D1a)*dqb[2][k][j][i]*dist);
#else
      // this is a little different than sweep.c
      fraczone=D1a*INVDTOTLIMIT;

      dqminus=v[1][2][km1][j][i] + (1.0-DTOTLIMIT)*dqv[2][km1][j][i]*dist ;
      dqplus=v[1][2][k][j][i] + (-1.0+DTOTLIMIT)*dqv[2][k][j][i]*dist ;
      vm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      vp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));

      dqminus=v[2][2][km1][j][i] + (1.0-DTOTLIMIT)*dqb[2][km1][j][i]*dist ;
      dqplus=v[2][2][k][j][i] + (-1.0+DTOTLIMIT)*dqb[2][k][j][i]*dist ;
      bm = 0.5*(dqminus*(1.0-fraczone)+dqplus*(1.0+fraczone));
      bp = 0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
#endif
    }
#endif

#elif(RHOINTERP==1)
    // no strong sym
#if(SYMFORCEMAG==0)
      vp = v[1][2][km1][j][i] + (1.0-D1p)*dqv[2][km1][j][i]*dist ;
      bp = v[2][2][km1][j][i] + (1.0-D1p)*dqb[2][km1][j][i]*dist ;
      vm = v[1][2][k][j][i] + (-1.0+D1m)*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0+D1m)*dqb[2][k][j][i]*dist ;
#else
    if(D1p>DTOTLIMIT){
      vp = v[1][2][km1][j][i] + (1.0-D1p)*dqv[2][km1][j][i]*dist ;
      bp = v[2][2][km1][j][i] + (1.0-D1p)*dqb[2][km1][j][i]*dist ;
    }
    else{
      vp = 0.5*(v[1][2][km1][j][i] + (1.0-D1p)*dqv[2][km1][j][i]*dist+
		v[1][2][k][j][i] + (-1.0+D1m)*dqv[2][k][j][i]*dist);
      bp = 0.5*(v[2][2][km1][j][i] + (1.0-D1p)*dqb[2][km1][j][i]*dist+
		v[2][2][k][j][i] + (-1.0+D1m)*dqb[2][k][j][i]*dist);
    }
    if(D1m>DTOTLIMIT){
      vm = v[1][2][k][j][i] + (-1.0+D1m)*dqv[2][k][j][i]*dist ;
      bm = v[2][2][k][j][i] + (-1.0+D1m)*dqb[2][k][j][i]*dist ;
    }
    else{
      vm = 0.5*(v[1][2][km1][j][i] + (1.0-D1p)*dqv[2][km1][j][i]*dist+
		v[1][2][k][j][i] + (-1.0+D1m)*dqv[2][k][j][i]*dist);
      bm = 0.5*(v[2][2][km1][j][i] + (1.0-D1p)*dqb[2][km1][j][i]*dist+
		v[2][2][k][j][i] + (-1.0+D1m)*dqb[2][k][j][i]*dist);      
    }
#endif

#endif // end if rhointerp==1

#if(RHOINTERP==0)
    b2prim[2][k][j][i] = 0.5*(bm + bp + sgn_va*(vm - vp)/srhoa[k][j][i]) ;
#elif(RHOINTERP==1)
    //srhom=z2e_2(osqrtrho,k,j,i); // C- char rho value
    //  srhop=z2e_2(osqrtrho,km1,j,i); // C+ char rho value
    b2prim[2][k][j][i] = (sgn_va*(vm-vp)+bp*srhop+bm*srhom)/(srhom+srhop);
#elif(RHOINTERP==2)
    fprintf(fail_file,"no such rho interp ready: %d",RHOINTERP);
    myexit(1);
#endif


  }
  }// endif N3>1
  // now yx/bystar live on corners on full grid, act like vectors on rotation

  }// end if not true that both N2==1 and N3==1



  // now all needed prim quantities are computed, can compute transverse lorentz velocities

  // update velocity(2nd partial update of Lorentz force terms resulting in Alfven wave motion) (comoving)
  if(!( (N2==1)&&(N3==1))){
    if(mocctvx1){
      LOOPV1{ // needs b1prim[1]/b2prim[1] on i=0..N1-1 and j=0..N2, k=0..N3 (overdone, but ok)
#if(ALFVENLIMIT==0)
	rhoa2 = z2e_1(s[1],k,j,i);
#else
	rhoa2 = z2e_1(rholimited[0],k,j,i);
#endif
#if(N2!=1)
	bya = v2tov1(v[2][2],k,j,i);
#endif
#if(N3!=1)
	bza = v3tov1(v[2][3],k,j,i);
#endif
	v[1][1][k][j][i] += dt/rhoa2*(
#if(N2!=1)
				      (bya*OARC42(k,j,i)*(b1prim[1][k][jp1][i] - b1prim[1][k][j][i]))
#else
				      0
#endif
#if(N3!=1)
				      + (bza*OARC23(k,j,i)*ODX(1,3,k)*(b2prim[1][kp1][j][i]-b2prim[1][k][j][i]))
#endif
				      );
      }
    }
  }
  if(!( (N1==1)&&(N3==1) )){
    if(mocctvx2){
      LOOPV2{ // needs b1prim[2]/b2prim[2] at i=0..N1 and j=0..N2-1+periodicx2special (LOOPH above is fine) k=0..N3 (overdone, but ok)
#if(ALFVENLIMIT==0)
	rhoa2 =z2e_2(s[1],k,j,i);
#else
	rhoa2 =z2e_2(rholimited[0],k,j,i);
#endif
#if(N1!=1)
	bxa = v1tov2(v[2][1],k,j,i);
#endif
#if(N3!=1)
	bza = v3tov2(v[2][3],k,j,i);
#endif
	v[1][2][k][j][i] += dt/rhoa2*(
#if(N1!=1)
				      (bxa*OARC41(k,j,i)/(G2(2,i))*(G2(1,ip1)*b1prim[2][k][j][ip1] - G2(1,i)*b1prim[2][k][j][i]))
#else
				      0
#endif
#if(N3!=1)
				      +(bza*OARC33(k,j,i)*ODX(1,3,k)*(b2prim[2][kp1][j][i]-b2prim[2][k][j][i]))
#endif
				      );
      }
    }
  }
  if(!( (N1==1)&&(N2==1) )){
    if(mocctvx3){
      LOOPV3{ // needs b1prim[3]/b2prim[3] at k=0..N3-1 and i/j=0..N1/2 (overdone, but ok)
#if(ALFVENLIMIT==0)
	rhoa2 =z2e_3(s[1],k,j,i);
#else
	rhoa2 =z2e_3(rholimited[0],k,j,i);
#endif
#if(N1!=1)
	bxa = v1tov3(v[2][1],k,j,i);
#endif
#if(N2!=1)
	bya = v2tov3(v[2][2],k,j,i);
#endif
	v[1][3][k][j][i] += dt/rhoa2*(
#if(N1!=1)
				      (bxa*OARC41(k,j,i)/(G3(2,i))*(G3(1,ip1)*b1prim[3][k][j][ip1] - G3(1,i)*b1prim[3][k][j][i]))
#else
				      0
#endif
#if(N2!=1)
				      +(bya*OARC32(k,j,i)/(G4(2,j))*(G4(1,jp1)*b2prim[3][k][jp1][i]-G4(1,j)*b2prim[3][k][j][i]))
#endif
				      );
      }
    }
  }



  // limit what's bounded to only necessary changing items for faster reduced dimensional speeds
  if(!( (N2==1)&&(N3==1) ) ){ // then do 1-comp
    if(!( (N1==1)&&(N3==1) ) ){ // then do 2-comp
      if(!( (N1==1)&&(N2==1) )){ // then do 3-comp
	bound(NULL,NULL,0,1,123);
      }
      else{
	bound(NULL,NULL,0,1,12); // no 3
      }
    }
    else{ // no 2
      if(!( (N1==1)&&(N2==1) )){
	bound(NULL,NULL,0,1,13);
      }
      else{
	bound(NULL,NULL,0,1,1); // no 3
      }
    }
  }
  else{// no 1
    if(!( (N1==1)&&(N3==1) ) ){ // then do 2-comp
      if(!( (N1==1)&&(N2==1) )){ // then do 3-comp
	bound(NULL,NULL,0,1,23);
      }
      else{
	bound(NULL,NULL,0,1,2); // no 3
      }
    }
    else{ // no 2
      if(!( (N1==1)&&(N2==1) )){
	bound(NULL,NULL,0,1,3);
      }
      else{
	// never get here // no 3
      }
    }
  }


}



// 1d/2d/3d valid
// currents are located as each component of emf, such that curl of curl is a natural vector location
void current_compute(int wcom)
{
  int i,j,k;
  int docom[3+1];

  // determine which components to do, if any
  if( (wcom==1)||(wcom==12)||(wcom==13)||(wcom==123) ) docom[1]=1; else docom[1]=0;
  if( (wcom==2)||(wcom==12)||(wcom==23)||(wcom==123) ) docom[2]=1; else docom[2]=0;
  if( (wcom==3)||(wcom==13)||(wcom==23)||(wcom==123) ) docom[3]=1; else docom[3]=0;
  if(wcom==0){ docom[1]=docom[2]=docom[3]=0; }

  // doesn't use Bx(i=-2) for example, just offsets like By(i=-2)
  if(docom[1]){
    LOOPHMFPC{// needed for nu_res_compute below( bit more than needed)
      jcurrent[1][k][j][i]=curlvfornat1(v[2],k,j,i);
    }
  }
  if(docom[2]){
    LOOPHMFPC{ // bit more than needed
      jcurrent[2][k][j][i]=curlvfornat2(v[2],k,j,i);
    }
  }
  if(docom[3]){
    LOOPHMFPC{ // needed for nu_res_compute below(exactly what needed)
      jcurrent[3][k][j][i]=curlvfornat3(v[2],k,j,i);
    }
  }
}



// 1d/2d/3d valid
void nu_res_compute(void)
{
  int i,j,k;
  FTYPE ftemp;
  FTYPE ftemp2;
  FTYPE cs2;
  FTYPE jx,jy,jz;



  ftemp=sqrt(gam*(gam-1.));
  ftemp2=gam*(gam-1.);

  // no need to bound since can do LOOPFC since nu_real only depends on zone centered quantities

  if(rreal==1){    
    LOOPHC{
      // gammie version
      if(wgam) cs = ftemp*sqrt(s[2][k][j][i]/s[1][k][j][i]);
      nu_res_real[k][j][i]=resist_real*nu_res_fact[k][j][i]*cs;
    }
  }
  else if(rreal==2){
    LOOPHC{ // most needed(by moc_ct)
      // assumes have all currents updated
      // current at emf location
      jx=c2z_1(jcurrent[1],k,j,i);
      jy=c2z_2(jcurrent[2],k,j,i);
      jz=c2z_3(jcurrent[3],k,j,i);
      // SP00 version
      nu_res_real[k][j][i]=resist_real*nu_res_fact[k][j][i]*sqrt((jx*jx+jy*jy+jz*jz)/s[1][k][j][i]);
    }
  }
  else if(rreal==3){
    LOOPHC{ // most needed(by moc_ct)
      nu_res_real[k][j][i]=resist_real*nu_res_fact[k][j][i];
    }
  }
}


// 1d/2d/3d valid
void step_res(void)
{
  static FTYPE jx,jy,jz ;
  static FTYPE (*resj)[N3M][N2M][N1M];
  register int i,j,k ;
  FTYPE dxdx,dydy;
  FTYPE ftemp;
  

  if(ie&&resheat&&wgam){
    // update internal energy 
    LOOPC{ // jcurrent fine assuming LOOPH below is satisfied
      
      jx=c2z_1(jcurrent[1],k,j,i); // current in natural current location(not vector location)
      jy=c2z_2(jcurrent[2],k,j,i);
      jz=c2z_3(jcurrent[3],k,j,i);
      s[2][k][j][i] += dt*nu_res_real[k][j][i]*(jx*jx + jy*jy + jz*jz) ;
    }

    if(FORCEIEINTERNAL){
      floor_correct(1,0);
    }
    bound(NULL,NULL,2,0,0) ;
  }

  if((COMPDIM<=2)&&(MOC2DVER==0)){ // otherwise done in moc_ct routine completely self-consistently
    resj=workv1; // \nu \times J

    // then do magnetic fields(induction eq.)
    LOOPHPC{
      // needs nu_res_real: i=0..N1, j=0..N2
      // note that using natural current location, thus why z2e_x "odd".
      resj[1][k][j][i]=jcurrent[1][k][j][i]*z2c_1(nu_res_real,k,j,i);
      resj[2][k][j][i]=jcurrent[2][k][j][i]*z2c_2(nu_res_real,k,j,i);
    }
    
    LOOPC{
      // needs resj[2]: i=0..N1, j=0..N2-1 inclusive
      // needs resj[1]: i=0..N1-1, j=0..N2 inclusive (LOOPHPC above bit much, but ok)
      v[2][3][k][j][i]+=-dt*curlvbacknat3(resj,k,j,i);
    }

    bound(NULL,NULL,0,2,3); // bound 3 component
  }
}



