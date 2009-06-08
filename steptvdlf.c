// changes from TVDLF in gammie code
// NX->N1
// NY->N2
// double -> FTYPE
// removed decs.h and put in global.h/defs.h defines
// took defs.h -> global.h for new stuff.  Took decs.h -> defs.h for new stuff
// LOOP -> LOOPTVDLF
// PLOOP -> LOOPP
// set_arrays() is in init_pointers()
// various modifications to my code to translate primitive variables so can diag
// advint instead of #if(1) craziness, and SYMFLUX instead of #if(1)
// dx -> dxTVD dy-> dyTVD dV ->dVTVD
// s->pr

// i-1 -> im1
// i+1 -> ip1
// etc. like in my code
// then bounds() only bounds necessary direction

#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif


#define SYMFLUX 1 // see below
/* Toth/Rusanov/Yee  method */
#define FAC	0.5	

void step_ch()
{
	int i,j,k ;
	FTYPE Dqm,Dqp,Dqc,pr ;
	FTYPE (* Fx_ct)[N2+2][NP] ;
	FTYPE (* Fy_ct)[N2+2][NP] ;
	FTYPE cmax,cs2x,cs2y,va2x,va2y,cmsx,cmsy,vx,vy ;
	static FTYPE cx[N1][N2],cy[N1][N2] ;
	static FTYPE c = 1. ;
	static FTYPE Uh[NP],U[NP] ;
	static FTYPE pix[NP],plx[NP],prx[NP],Flx[NP],Frx[NP],Ulx[NP],Urx[NP] ;
	static FTYPE piy[NP],ply[NP],pry[NP],Fly[NP],Fry[NP],Uly[NP],Ury[NP] ;
	void bounds() ;
	void primtoU(FTYPE *pb, FTYPE *Ub) ;
	void Utoprim(FTYPE *Ua, FTYPE *pa) ;
	void primtoxflux(FTYPE *pa, FTYPE *Fa) ;
	void primtoyflux(FTYPE *pa, FTYPE *Fa) ;


	//	dump(NULL,666+nstep,DTYPE,0);
	/*
	  LOOP
    {
      if(s[1][k][j][i]<1E-20){
	fprintf(fail_file,"s:input data has mass density <1E-20: %d %d %d %15.10g\n",k,j,i,s[1][k][j][i]);
	//myexit(1);
      }
      if(s[2][k][j][i]<1E-20 ){
	fprintf(fail_file,"s:input data has ie density <1E-20: %d %d %d %15.10g\n",k,j,i,s[2][k][j][i]);
	//myexit(1);
      }
    }
    LOOPTVDLF{
      if(p[i][j][0]<1E-20){
	fprintf(fail_file,"p:input data has mass density <1E-20: %d %d %d %15.10g\n",k,j,i,p[i][j][0]);
	//myexit(1);
      }
      if(p[i][j][7]<1E-20 ){
	fprintf(fail_file,"p:input data has ie density <1E-20: %d %d %d %15.10g\n",k,j,i,p[i][j][7]);
	//myexit(1);
      }
    }
//    fflush(fail_file);
 //   fprintf(stdout,"1 n: %ld ie: %15.10g\n",nstep,p[306][1023][7]);
    fflush(stdout);
	*/
	/** set timestep **/
	dt = cour*dxTVD/1. ;

        /* don't step beyond end of run */
        if(t + dt > tf) dt = tf - t ;

	bounds() ;
	if(advint==2){
	/* evaluate Woodward slopes of primitive variables */
	LOOPTVDLF {
		LOOPP {
			Dqm = 2.0*(p[i][j][k] - p[im1][j][k]) ;
			Dqp = 2.0*(p[ip1][j][k] - p[i][j][k]) ;
			Dqc = 0.5*(p[ip1][j][k] - p[im1][j][k]) ;
			pr = Dqm*Dqp ;
			if(pr <= 0.) dqx[i][j][k] = 0. ;
			else {
				if(fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
					dqx[i][j][k] = Dqm ;
				else if(fabs(Dqp) < fabs(Dqc))
					dqx[i][j][k] = Dqp ;
				else
					dqx[i][j][k] = Dqc ;
			}
		}
		LOOPP {
			Dqm = 2.0*(p[i][j][k] - p[i][jm1][k]) ;
			Dqp = 2.0*(p[i][jp1][k] - p[i][j][k]) ;
			Dqc = 0.5*(p[i][jp1][k] - p[i][jm1][k]) ;
			pr = Dqm*Dqp ;
			if(pr <= 0.) dqy[i][j][k] = 0. ;
			else {
				if(fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
					dqy[i][j][k] = Dqm ;
				else if(fabs(Dqp) < fabs(Dqc))
					dqy[i][j][k] = Dqp ;
				else
					dqy[i][j][k] = Dqc ;
			}
		}
	}
	}
	if(advint==1){
	/* evaluate vanleer slopes */
	LOOPTVDLF {
		LOOPP {
			Dqm = p[i][j][k] - p[im1][j][k] ;
			Dqp = p[ip1][j][k] - p[i][j][k] ;
			pr = Dqm*Dqp ;
			if(pr <= 0.) dqx[i][j][k] = 0. ;
			else dqx[i][j][k] = 2.*pr/(Dqm+Dqp) ;
		}
		LOOPP {
			Dqm = p[i][j][k] - p[i][jm1][k] ;
			Dqp = p[i][jp1][k] - p[i][j][k] ;
			pr = Dqm*Dqp ;
			if(pr <= 0.) dqy[i][j][k] = 0. ;
			else dqy[i][j][k] = 2.*pr/(Dqm+Dqp) ;
		}
	}
	}
	bounds() ;

	/* evaluate timestep */
        cmax = 0. ;
        LOOPTVDLF {
                LOOPP {
                        pix[k] = 0.5*(
				p[im1][j][k] + 0.5*dqx[im1][j][k] +
                                p[i][j][k]   - 0.5*dqx[i][j][k]) ;
                        piy[k] = 0.5*(
				p[i][jm1][k] + 0.5*dqy[i][jm1][k] +
                                p[i][j][k]   - 0.5*dqy[i][j][k]) ;
                }
		/*
		if(pix[7]<0.){
		  fprintf(stdout,"%d %d pix[7]: %15.10g\n",i,j,pix[7]);
		}
		if(piy[7]<0.){
		  fprintf(stdout,"%d %d piy[7]: %15.10g\n",i,j,piy[7]);
		}
		if(pix[0]<0.){
		  fprintf(stdout,"%d %d pix[0]: %15.10g\n",i,j,pix[0]);
		}
		if(piy[0]<0.){
		  fprintf(stdout,"%d %d piy[0]: %15.10g\n",i,j,piy[0]);
		}
		fflush(stdout);
		*/
                /* evaluate speeds */
                cs2x = gam*(gam - 1.)*pix[UU]/pix[RHO];
                cs2y = gam*(gam - 1.)*piy[UU]/piy[RHO];
                va2x = (pix[BX]*pix[BX] + pix[BY]*pix[BY] + pix[BZ]*pix[BZ])/
		  pix[RHO] ;
                va2y = (piy[BX]*piy[BX] + piy[BY]*piy[BY] + piy[BZ]*piy[BZ])/
				piy[RHO] ;
                cmsx = sqrt(cs2x + va2x) ;
                cmsy = sqrt(cs2y + va2y) ;
		vx = fabs(pix[UX]) ;
		vy = fabs(piy[UY]) ;

		cx[i][j] = vx + cmsx ;
		cy[i][j] = vy + cmsy ;

                if(cmsx > cmax) cmax = cmsx ;
                if(cmsy > cmax) cmax = cmsy ;
        }
        dt = cour*dxTVD/cmax ;
        /* don't step beyond end of run */
        if(t + dt > tf) dt = tf - t ;


	/* perform Hancock predictor half-step */
	LOOPTVDLF {
		LOOPP {
			plx[k] = p[i][j][k] - 0.5*dqx[i][j][k] ;
			prx[k] = p[i][j][k] + 0.5*dqx[i][j][k] ;
			ply[k] = p[i][j][k] - 0.5*dqy[i][j][k] ;
			pry[k] = p[i][j][k] + 0.5*dqy[i][j][k] ;
		}

		primtoxflux(plx,Flx) ;
		primtoxflux(prx,Frx) ;

		primtoyflux(ply,Fly) ;
		primtoyflux(pry,Fry) ;

		primtoU(p[i][j],Uh) ;

		LOOPP {
			Uh[k] -= 0.5*dt*( (Frx[k] - Flx[k])/dxTVD 
				+ (Fry[k] - Fly[k])/dyTVD ) ;
			ph[i][j][k] = p[i][j][k] ;
		}

		Utoprim(Uh,ph[i][j]) ;

	}
	bounds() ;

	/* evaluate zone boundary fluxes at half-step */
	LOOPTVDLF {
		LOOPP {
			pix[k] = 0.5*(ph[im1][j][k] + 0.5*dqx[im1][j][k] +
			   	      ph[i][j][k]   - 0.5*dqx[i][j][k]) ;
			piy[k] = 0.5*(ph[i][jm1][k] + 0.5*dqy[i][jm1][k] +
			   	      ph[i][j][k]   - 0.5*dqy[i][j][k]) ;
		}
		primtoxflux(pix, Fx[i][j]) ;
		primtoyflux(piy, Fy[i][j]) ;
	}

	/* evaluate diffusive correction to flux */
	LOOPTVDLF {
		LOOPP {
			plx[k] = ph[im1][j][k] + 0.5*dqx[im1][j][k] ;
			prx[k] = ph[i][j][k]   - 0.5*dqx[i][j][k] ;

			ply[k] = ph[i][jm1][k] + 0.5*dqy[i][jm1][k] ;
			pry[k] = ph[i][j][k]   - 0.5*dqy[i][j][k] ;
		}
		primtoU(plx,Ulx) ;
		primtoU(prx,Urx) ;
		primtoU(ply,Uly) ;
		primtoU(pry,Ury) ;
		
		LOOPP {
			Fx[i][j][k] += FAC*cx[i][j]*(Ulx[k] - Urx[k]) ;
			Fy[i][j][k] += FAC*cy[i][j]*(Uly[k] - Ury[k]) ;
		}
	}
	bounds() ;

#if(SYMFLUX)
	/* replace B-field fluxes with flux-ct symmetrized fluxes */
	Fx_ct = dqx ;
	Fy_ct = dqy ;
	LOOPTVDLF {

		Fx_ct[i][j][BX] = 0. ;
		Fy_ct[i][j][BX] = 0.125*( 
			2.*Fy[i][j][BX]
			+ Fy[ip1][j][BX]
			+ Fy[im1][j][BX]
			- Fx[i][j][BY]
			- Fx[i][jm1][BY]
			- Fx[ip1][j][BY]
			- Fx[ip1][jm1][BY] ) ;
		Fx_ct[i][j][BY] = 0.125*( 
			2.*Fx[i][j][BY]
			+ Fx[i][jp1][BY]
			+ Fx[i][jm1][BY]
			- Fy[i][j][BX]
			- Fy[im1][j][BX]
			- Fy[i][jp1][BX]
			- Fy[im1][jp1][BX] ) ;
		Fy_ct[i][j][BY] = 0. ;
	}
	LOOPTVDLF {
		Fx[i][j][BX] = Fx_ct[i][j][BX] ;
		Fy[i][j][BX] = Fy_ct[i][j][BX] ;
		Fx[i][j][BY] = Fx_ct[i][j][BY] ;
		Fy[i][j][BY] = Fy_ct[i][j][BY] ;
	}
	bounds() ;
#endif

	LOOPTVDLF {
		/* calculate conserved quantities */
		primtoU(p[i][j], U) ;

		/* evolve conserved quantities */
		LOOPP   U[k] -= (dt/dxTVD)*(Fx[ip1][j][k] - Fx[i][j][k]) 
			+ (dt/dyTVD)*(Fy[i][jp1][k] - Fy[i][j][k]) ;

		/* recover new primitive variables */
		Utoprim(U, p[i][j]) ;

	}
	/* apply boundary conditions */
	bounds() ;

	/* done! */

	t += dt ;

	//    fprintf(stdout,"9 n: %ld ie: %15.10g\n",nstep,p[306][1023][7]);
	//fflush(stdout);
}

// END OF step_ch()


/** more Physics **/

void Utoprim(FTYPE *U, FTYPE *pr)
{
	FTYPE u1,u2,u3,u4,u5,u6,u7,u8 ;

	u1 = U[0] ;
	u2 = U[1] ;
	u3 = U[2] ;
	u4 = U[3] ;
	u5 = U[4] ;
	u6 = U[5] ;
	u7 = U[6] ;
	u8 = U[7] ;

	pr[RHO] = u1 ;
	pr[UX] = u2/(u1+DENSITYFLOOR) ;
	pr[UY] = u3/(u1+DENSITYFLOOR) ;
	pr[UZ] = u4/(u1+DENSITYFLOOR) ;
	pr[BX] = u5 ;
	pr[BY] = u6 ;
	pr[BZ] = u7 ;
	pr[UU] = u8 - 0.5*(u5*u5 + u6*u6 + u7*u7 + (u2*u2 + u3*u3 + u4*u4)/(u1+DENSITYFLOOR)) ;

}

void primtoxflux(FTYPE *pr, FTYPE *f)
{
	FTYPE r,vx,vy,vz,bx,by,bz,u ;

	r = pr[RHO] ;
	vx = pr[UX] ;
	vy = pr[UY] ;
	vz = pr[UZ] ;
	bx = pr[BX] ;
	by = pr[BY] ;
	bz = pr[BZ] ;
	u = pr[UU] ;

	f[0] = r*vx ;
	f[1] = r*vx*vx + (gam - 1)*u + 0.5*(by*by + bz*bz - bx*bx) ;
	f[2] = -bx*by + r*vx*vy ;
	f[3] = -bx*bz + r*vx*vz ;
	f[4] = 0. ;
	f[5] = by*vx - bx*vy ;
	f[6] = bz*vx - bx*vz ;
	f[7] = vx*(r*0.5*(vx*vx + vy*vy + vz*vz) + by*by + bz*bz + gam*u) 
		- bx*(by*vy + bz*vz) ;
}

void primtoyflux(FTYPE *pr, FTYPE *f)
{
	FTYPE r,vx,vy,vz,bx,by,bz,u ;

	r = pr[RHO] ;
	vx = pr[UX] ;
	vy = pr[UY] ;
	vz = pr[UZ] ;
	bx = pr[BX] ;
	by = pr[BY] ;
	bz = pr[BZ] ;
	u = pr[UU] ;

	f[0] = r*vy ;
	f[1] = -by*bx + r*vy*vx ;
	f[2] = r*vy*vy + (gam - 1)*u + 0.5*(-by*by + bz*bz + bx*bx) ;
	f[3] = -by*bz + r*vy*vz ;
	f[4] = bx*vy - by*vx ;
	f[5] = 0. ;
	f[6] = bz*vy - by*vz ;
	f[7] = vy*(r*0.5*(vx*vx + vy*vy + vz*vz) + bx*bx + bz*bz + gam*u) 
		- by*(bx*vx + bz*vz) ;

}

void primtoU(FTYPE *pr, FTYPE *U)
{
	FTYPE r,vx,vy,vz,bx,by,bz,u ;

	r = pr[RHO] ;
	vx = pr[UX] ;
	vy = pr[UY] ;
	vz = pr[UZ] ;
	bx = pr[BX] ;
	by = pr[BY] ;
	bz = pr[BZ] ;
	u = pr[UU] ;

	U[0] = r ;
	U[1] = r*vx ;
	U[2] = r*vy ;
	U[3] = r*vz ;
	U[4] = bx ;
	U[5] = by ;
	U[6] = bz ;
	U[7] = 0.5*(bx*bx + by*by + bz*bz) + u + 0.5*r*(vx*vx + vy*vy + vz*vz) ;
}


void bounds()
{
	int i,j,k ;

	if(N2>1){
	  if( (bcix2==5)&&(bcox2==5) ){
	    for(i=0;i<N1;i++) {
	      LOOPP {
		/* lower boundary */
		p  [i][-1][k] = p  [i][N2-1][k] ;
		ph [i][-1][k] = ph [i][N2-1][k] ;
		Fx [i][-1][k] = Fx [i][N2-1][k] ;
		Fy [i][-1][k] = Fy [i][N2-1][k] ;
		dqx[i][-1][k] = dqx[i][N2-1][k] ;
		dqy[i][-1][k] = dqy[i][N2-1][k] ;
		
		/* upper boundary */
		p  [i][N2][k] = p  [i][0][k] ;
		ph [i][N2][k] = ph [i][0][k] ;
		Fx [i][N2][k] = Fx [i][0][k] ;
		Fy [i][N2][k] = Fy [i][0][k] ;
		dqx[i][N2][k] = dqx[i][0][k] ;
		dqy[i][N2][k] = dqy[i][0][k] ;
	      }
	    }
	  }
	  else if( (bcix2==4)&&(bcox2==4) ){
	    for(i=0;i<N1;i++) {
	      LOOPP {
		/* lower boundary */
		p  [i][-1][k] = p  [i][0][k] ;
		ph [i][-1][k] = ph [i][0][k] ;
		Fx [i][-1][k] = Fx [i][0][k] ;
		Fy [i][-1][k] = Fy [i][0][k] ;
		dqx[i][-1][k] = dqx[i][0][k] ;
		dqy[i][-1][k] = dqy[i][0][k] ;
		
		/* upper boundary */
		p  [i][N2][k] = p  [i][N2-1][k] ;
		ph [i][N2][k] = ph [i][N2-1][k] ;
		Fx [i][N2][k] = Fx [i][N2-1][k] ;
		Fy [i][N2][k] = Fy [i][N2-1][k] ;
		dqx[i][N2][k] = dqx[i][N2-1][k] ;
		dqy[i][N2][k] = dqy[i][N2-1][k] ;
	      }
	    }
	  }
	  else{
	    fprintf(fail_file,"x2: %d %d %d %d no such bound in tvdlf\n",bcix1,bcox1,bcix2,bcox2);
	    myexit(1);
	  }
	}
	if(N1>1){
	  if( (bcix1==5)&&(bcox1==5) ){
	    for(j=0;j<N2;j++) {
	      LOOPP {
		/* lower boundary */
		p  [-1][j][k] = p  [N1-1][j][k] ;
		ph [-1][j][k] = ph [N1-1][j][k] ;
		Fx [-1][j][k] = Fx [N1-1][j][k] ;
		Fy [-1][j][k] = Fy [N1-1][j][k] ;
		dqx[-1][j][k] = dqx[N1-1][j][k] ;
		dqy[-1][j][k] = dqy[N1-1][j][k] ;
		
		/* upper boundary */
		p  [N1][j][k] = p  [0][j][k] ;
		ph [N1][j][k] = ph [0][j][k] ;
		Fx [N1][j][k] = Fx [0][j][k] ;
		Fy [N1][j][k] = Fy [0][j][k] ;
		dqx[N1][j][k] = dqx[0][j][k] ;
		dqy[N1][j][k] = dqy[0][j][k] ;
	      }
	    }
	  }
	  else if( (bcix1==4)&&(bcox1==4) ){
	    for(j=0;j<N2;j++) {
	      LOOPP {
		/* lower boundary */
		p  [-1][j][k] = p  [0][j][k] ;
		ph [-1][j][k] = ph [0][j][k] ;
		Fx [-1][j][k] = Fx [0][j][k] ;
		Fy [-1][j][k] = Fy [0][j][k] ;
		dqx[-1][j][k] = dqx[0][j][k] ;
		dqy[-1][j][k] = dqy[0][j][k] ;
		
		/* upper boundary */
		p  [N1][j][k] = p  [N1-1][j][k] ;
		ph [N1][j][k] = ph [N1-1][j][k] ;
		Fx [N1][j][k] = Fx [N1-1][j][k] ;
		Fy [N1][j][k] = Fy [N1-1][j][k] ;
		dqx[N1][j][k] = dqx[N1-1][j][k] ;
		dqy[N1][j][k] = dqy[N1-1][j][k] ;
	      }
	    }
	  }
	  else{
	    fprintf(fail_file,"x1: %d %d %d %d no such bound in tvdlf\n",bcix1,bcox1,bcix2,bcox2);
	    myexit(1);
	  }
	}
	if( (bcix1==5)&&(bcox1==5)&&(bcix2==5)&&(bcox2==5) ){
	  /* corners */
	  LOOPP {
	    p  [-1][-1][k] = p  [N1-1][N2-1][k] ;
	    ph [-1][-1][k] = ph [N1-1][N2-1][k] ;
	    Fx [-1][-1][k] = Fx [N1-1][N2-1][k] ;
	    Fy [-1][-1][k] = Fy [N1-1][N2-1][k] ;
	    dqx[-1][-1][k] = dqx[N1-1][N2-1][k] ;
	    dqy[-1][-1][k] = dqy[N1-1][N2-1][k] ;
	    
	    p  [-1][N2][k] = p  [N1-1][0][k] ;
	    ph [-1][N2][k] = ph [N1-1][0][k] ;
	    Fx [-1][N2][k] = Fx [N1-1][0][k] ;
	    Fy [-1][N2][k] = Fy [N1-1][0][k] ;
	    dqx[-1][N2][k] = dqx[N1-1][0][k] ;
	    dqy[-1][N2][k] = dqy[N1-1][0][k] ;
	    
	    p  [N1][N2][k] = p  [0][0][k] ;
	    ph [N1][N2][k] = ph [0][0][k] ;
	    Fx [N1][N2][k] = Fx [0][0][k] ;
	    Fy [N1][N2][k] = Fy [0][0][k] ;
	    dqx[N1][N2][k] = dqx[0][0][k] ;
	    dqy[N1][N2][k] = dqy[0][0][k] ;
	    
	    p  [N1][-1][k] = p  [0][N2-1][k] ;
	    ph [N1][-1][k] = ph [0][N2-1][k] ;
	    Fx [N1][-1][k] = Fx [0][N2-1][k] ;
	    Fy [N1][-1][k] = Fy [0][N2-1][k] ;
	    dqx[N1][-1][k] = dqx[0][N2-1][k] ;
	    dqy[N1][-1][k] = dqy[0][N2-1][k] ;
	  }
	}
	else if( (bcix1==4)&&(bcox1==4)&&(bcix2==4)&&(bcox2==4) ){
	  /* corners */
	  LOOPP {
	    p  [-1][-1][k] = p  [0][0][k] ;
	    ph [-1][-1][k] = ph [0][0][k] ;
	    Fx [-1][-1][k] = Fx [0][0][k] ;
	    Fy [-1][-1][k] = Fy [0][0][k] ;
	    dqx[-1][-1][k] = dqx[0][0][k] ;
	    dqy[-1][-1][k] = dqy[0][0][k] ;
	    
	    p  [-1][N2][k] = p  [0][N2-1][k] ;
	    ph [-1][N2][k] = ph [0][N2-1][k] ;
	    Fx [-1][N2][k] = Fx [0][N2-1][k] ;
	    Fy [-1][N2][k] = Fy [0][N2-1][k] ;
	    dqx[-1][N2][k] = dqx[0][N2-1][k] ;
	    dqy[-1][N2][k] = dqy[0][N2-1][k] ;
	    
	    p  [N1][N2][k] = p  [N1-1][N2-1][k] ;
	    ph [N1][N2][k] = ph [N1-1][N2-1][k] ;
	    Fx [N1][N2][k] = Fx [N1-1][N2-1][k] ;
	    Fy [N1][N2][k] = Fy [N1-1][N2-1][k] ;
	    dqx[N1][N2][k] = dqx[N1-1][N2-1][k] ;
	    dqy[N1][N2][k] = dqy[N1-1][N2-1][k] ;
	    
	    p  [N1][-1][k] = p  [N1-1][0][k] ;
	    ph [N1][-1][k] = ph [N1-1][0][k] ;
	    Fx [N1][-1][k] = Fx [N1-1][0][k] ;
	    Fy [N1][-1][k] = Fy [N1-1][0][k] ;
	    dqx[N1][-1][k] = dqx[N1-1][0][k] ;
	    dqy[N1][-1][k] = dqy[N1-1][0][k] ;
	  }

	}
	else{
	  fprintf(fail_file,"corner: %d %d %d %d no such bound in tvdlf\n",bcix1,bcox1,bcix2,bcox2);
	  myexit(1);
	}
}



void steptvdlf(void)
{
  step_ch();
}

void zeus2tvdlf(void)
{
  int i,j,k;
  // transfer jon to gammie variable format
  LOOPH{
    p[i][j][0]=s[1][k][j][i];
    p[i][j][1]=v[1][1][k][j][i];
    p[i][j][2]=v[1][2][k][j][i];
    p[i][j][3]=v[1][3][k][j][i];
    p[i][j][4]=v[2][1][k][j][i];
    p[i][j][5]=v[2][2][k][j][i];
    p[i][j][6]=v[2][3][k][j][i];
    p[i][j][7]=s[2][k][j][i];
    // no potential
  }
  // TVD code is uniform and cartesian
  dxTVD=dx[1][1][0];
  dyTVD=dx[1][2][0];
  dVTVD=dxTVD*dyTVD;
}

void tvdlf2zeus(void)
{
  int i,j,k;
  // transfer jon to gammie variable format
  LOOPH{
    s[1][k][j][i]=p[i][j][0];
    v[1][1][k][j][i]=p[i][j][1];
    v[1][2][k][j][i]=p[i][j][2];
    v[1][3][k][j][i]=p[i][j][3];
    v[2][1][k][j][i]=p[i][j][4];
    v[2][2][k][j][i]=p[i][j][5];
    v[2][3][k][j][i]=p[i][j][6];
    s[2][k][j][i]=p[i][j][7];
    // no potential, wouldn't change anyways.
  }
}









