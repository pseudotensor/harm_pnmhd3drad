#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif


//=======================================================================
/////////////////////////  SUBROUTINE RHS  \\\\\\\\\\\\\\\\\\\\\\\\\\\\

void rhs(void)
{
#if(RAD)

/*  PURPOSE: Computes the finite differenced radiation and gas energy
  equations as fn1(e,en,er,ern) and fn2(e,en,er,ern). These functions
  are used as the RHS of the Newton-Raphson matrix eqn.  Note the -ve
  sign appearing in the definition of the Jacobian is added when
  the moment eqn matrix is set up (RICCG), ie Ax = -rhs.
  The N-R scheme will continue to iterate en and ern until fn1=fn2=0.

  INPUT ARGUMENTS:[none]

  OUTPUT ARGUMENTS:
   fn1 = radiation energy equation 
   fn2 = gas       energy equation

  EXTERNALS:[none]

  LOCALS:
-----------------------------------------------------------------------
*/

  static FTYPE (*erb)[N2M][N1M] ;
  int i,j,k ;
  FTYPE qa,qb,qc,qd,qe;
/*=======================================================================
  INLINE FUNCTIONS
   TH BAR: time centers a variable according to the value of radth
=======================================================================
*/

#define THEBAR(xn,xnp1) (radth *xnp1 + (1.0-radth )*xn)

  erb=work1 ;  

  LOOPRC {

	  erb[k][j][i] = THEBAR(er[k][j][i],ern[k][j][i]);
  }

//  Compute the GRAD(v):P term

  LOOPC {
          fn1[k][j][i] = (dv[1][1][k][j][i]*f[1][1][k][j][i]+
	 dv[2][2][k][j][i]*f[2][2][k][j][i] + 
	 dv[3][3][k][j][i]*(1.0-f[1][1][k][j][i]-f[2][2][k][j][i])+ 
	(dv[1][2][k][j][i]+ dv[2][1][k][j][i])*f[1][2][k][j][i])*erb[k][j][i];
  }

//  Now sum DIV(F) terms into moment eqns. 

  LOOPC {
          qa = G2(1,i+1)*G3(1,i+1)/dx[2][1][i+1]*dr[1][k][j][i+1];
          qb = G2(1,i  )*G3(1,i  )/dx[2][1][i]*dr[1][k][j][i];
          qc = G2(2,i)*G2(2,i) ;
          qd = G4(1,j+1)/dx[2][2][j+1]*dr[2][k][j+1][i];
          qe = G4(1,j  )/dx[2][2][j  ]*dr[2][k][j][i] ;
	  fn1[k][j][i] = fn1[k][j][i]
       -(qa*(erb[k][j][i+1]-erb[k][j][i])-qb*(erb[k][j][i]-erb[k][j][i-1]))/dvl[1][1][i]
       -(qd*(erb[k][j+1][i]-erb[k][j][i])-qe*(erb[k][j][i]-erb[k][j-1][i]))/qc/dvl[1][2][j];
  }


//  Now sum remaining terms into moment eqns.

  LOOPC {
	  qa = 4.0*M_PI*THEBAR(kap[k][j][i],kapn[k][j][i])*THEBAR(bb[k][j][i],bbn[k][j][i])
                - CCON*THEBAR(kap[k][j][i],kapn[k][j][i])*erb[k][j][i];
          fn1[k][j][i]= ern[k][j][i] - er[k][j][i] - dt*(qa - fn1[k][j][i]);
	  fn2[k][j][i] = en[k][j][i] - e[k][j][i]
                 + dt*(qa + THEBAR(pre[k][j][i],pren[k][j][i])*divv[k][j][i]);
  }
#endif
}
