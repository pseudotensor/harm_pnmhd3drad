#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif


//=====================================================================
///////////////////////////  SUBROUTINE DERIVS  \\\\\\\\\\\\\\\\\\\\\\\

      void derivs(void)
{
#if(RAD)
/*
  PURPOSE:  Computes the derivatives of the moment eqns with respect
  to (wrt) the material energy density (e) and the radiation energy
  density (er) at the advanced time level.

  INPUT ARGUMENTS:[none]

  OUTPUT ARGUMENTS:
   dfn1de,dfn1der = derivatives of the "fn1" moment equation (radiation
    energy eqn) wrt material energy and radiation energy density (e,er) 
   dfn2de,dfn2der = derivatives of the "fn2" moment equation (gas
    energy eqn) wrt material energy and radiation energy density (e,er) 

  EXTERNALS: [none]

  LOCALS:
//-----------------------------------------------------------------------
*/
  int i,j,k;
  FTYPE qa,tempera,tempera3,dtemperade;

#define THEBAR(xn,xnp1) (radth *xnp1 + (1.0-radth )*xn)

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
/*=======================================================================
  Compute of derivatives of moment eqns (fn1,fn2) wrt
  material internal energy density (e).
*/
  LOOPC {
	  qa = 4.0*M_PI*radth*dkapde[k][j][i]*THEBAR(bb[k][j][i],bbn[k][j][i])
        +  4.0*M_PI*radth*THEBAR(kap[k][j][i],kapn[k][j][i])*dbbde[k][j][i]
        - CCON*radth*dkapde[k][j][i]*THEBAR(er[k][j][i],ern[k][j][i]);
        dfn1de[k][j][i] = -dt*qa;
        dfn2de[k][j][i] = 1.0 + dt*(qa+radth*dpde[k][j][i]*divv[k][j][i]);
  }
/*
  Compute derivatives of moment eqns(fn1,f2) wrt radiation
  internal energy density (er).
  Note the derivative of GRAD(v) wrt er term is included in dfn1der,
  but the derivative of DIV(F) term wrt er is added in RICCG
*/
  LOOPC {
        qa = CCON*THEBAR(kap[k][j][i],kapn[k][j][i]);
        dfn1der[k][j][i] = 1.0 + dt*(qa
          + radth*(dv[1][1][k][j][i]*f[1][1][k][j][i]+dv[2][2][k][j][i]*f[2][2][k][j][i]
          + dv[3][3][k][j][i]*(1.0-f[1][1][k][j][i]-f[2][2][k][j][i])
          + (dv[1][2][k][j][i] + dv[2][1][k][j][i])*f[1][2][k][j][i]));
        dfn2der[k][j][i] = -dt*qa;
  }
#if(RADTEST)
 fprintf(rad_file,"dfn1de,dfn2de,dfn1der,dfn2der");
 for(i=1;i<=1;i++){
 fprintf(rad_file," %g %g %g %g",dfn1de[0][0][i],dfn2de[0][0][i],dfn1der[0][0][i],dfn2der[0][0][i]);
 tempera=(gam - 1.f) * en[0][0][i]/s[1][0][0][i]*mmw*mhydrnew/kcontnew;
 tempera3=pow(tempera,3);
 dtemperade=tempera/en[0][0][i];
 fprintf(rad_file,"\n T T3 dtde dbbde %g %g %g %g",tempera,tempera3,dtemperade,dbbde[0][0][i]);
 }
 fprintf(rad_file,"\n");
#endif
#endif
}
