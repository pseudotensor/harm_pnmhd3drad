#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif

////////////////////////  SUBROUTINE GRADV  \\\\\\\\\\\\\\\\\\\\\\\\\\

    void gradv(void)
{
#if(RAD)
/*
  PURPOSE: This routine computes the components of the velocity
  gradient tensor.  The dv3i and dvi3 [i=1,2] componants are not
  needed since the corresponding tensor
  variable Eddington factor components are zero, and the velocity
  gradient always appears as the product GRAD(v):fE.  The velocity
  divergence DIV(v) is also computed here.  Note that although
  DIV(v) = dv11 + dv22 + dv33, we use the full coordinate independent
  formula to compute DIV(v) separately for accuracy.  The dv11,dv22,
  dv33 and divv are zone centered.  The dv12 and dv21 are averaged to
  zone centers.

  EXTERNALS: [none]

  LOCALS:
   dv** = ** component of the velocity gradient tensor
   divv = DIV(v)
----------------------------------------------------------------------
*/
  int i,j,k,ilower,iupper ;
  FTYPE  qa,qb,v1zc,(*t12)[N2M][N1M],(*t21)[N2M][N1M] ;

  t12=work1;
  t21=work2;

  LOOPC {
          v1zc    = 0.5*(v[1][1][k][j][i] + v[1][1][k][j][i+1]);
          qa = G2(1,i+1)*G3(1,i+1)*v[1][1][k][j][i+1];
          qb = G2(1,i  )*G3(1,i  )*v[1][1][k][j][i];

          dv[1][1][k][j][i] = (v[1][1][k][j][i+1]-v[1][1][k][j][i])/dx[1][1][i];
          dv[2][2][k][j][i] =  v1zc/G2(2,i)*DG2(2,i)
          + (v[1][2][k][j+1][i]-v[1][2][k][j][i])/(G2(2,i)*dx[1][2][j]);
          dv[3][3][k][j][i] =  v1zc/G3(2,i)*DG3(2,i)
          + 0.5*(v[1][2][k][j][i] + v[1][2][k][j+1][i] )/(G2(2,i)*G4(2,j))*DG4(2,j);
          divv[k][j][i] = (qa - qb)/dvl[1][1][i]
          + (G4(1,j+1)*v[1][2][k][j+1][i] -G4(1,j)*v[1][2][k][j][i] )/(G2(2,i)*dvl[1][2][j]);
  }
/*
  Note for the off-diagonal components of GRAD(v) special formula are
  needed in RT geometry to prevent divide by zero at i=ii(j).
*/
  LOOPC3{
  LOOPDIVJ {   
//change this when looptype==2 or 3 
        ilower = 0;
        iupper = N1;
#if(RT)
        if (G2(1,0)==0.0) {ilower = 1;}
#endif
        for(i=ilower;i<=iupper;i++){
          t12[k][j][i] = (v[1][2][k][j][i] - v[1][2][k][j][i-1])/dx[2][1][i];
          t21[k][j][i] = (v[1][1][k][j][i] - v[1][1][k][j-1][i])/(G2(1,i)*dx[2][2][j])
           -0.5*(v[1][2][k][j][i] + v[1][2][k][j][i-1])/G2(1,i)*DG2(1,i);
        }
#if(RT)
        if (ilower!=0){
          t12[k][j][0] = (v[1][2][k][j][0] - v[1][2][k][j][-1])/dx[2][1][0];
          t21[k][j][0] = 0.0 ;
        }
#endif
  }
  }
/*
  average dv12 and dv21 to zone centers
*/
  LOOPC {
	  dv[1][2][k][j][i] =  c2z_3(t12,k,i,j);
	  dv[2][1][k][j][i] =  c2z_3(t21,k,i,j);
  }
#endif
}
