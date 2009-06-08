#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif

//=======================================================================
///////////////////////////  SUBROUTINE FLD  \\\\\\\\\\\\\\\\\\\\\\\\\\

void fld(void)
{
#if(RAD)

/*  PURPOSE:  Computes the flux limited diffusion coeficient in Fick's
  law for the radiation flux.  Also computes the components of the
  tensor variable Eddington factors from the flux limiters.  In this
  implementation, two forms for the flux limiter are implemented:
   (1) Chapman-Enskog theory    [selected via ifld=1]
   (2) piecewise linear Minerbo [selected via ifld=2]
  See Levermore & Pomraning, Ap.J., 248, 321 (1981) and Levermore,
  JQSRT, 31, 149 (1984) for more details of the forms of the
  flux limiters.  Note the diffusion coefficient are face centered,
  while all components of the tensor Eddington factors are zone
  centered.

  EXTERNALS: [none]

  LOCALS:
-----------------------------------------------------------------------
*/
  static FTYPE *tt,*dtde,*dkapdt,(*n1)[N2M][N1M],(*n2)[N2M][N1M],
		(*r)[N2M][N1M],(*t12)[N2M][N1M],(*tdr)[N2M][N1M],(*tfr)[N2M][N1M];
  int i,j,k,istart,iend;
  FTYPE   der1,der2,dernorm,chi,lmda,qa,qb;
     
  tt=work1x;
  dtde=work2x;
  dkapdt=work3x;
  n1=work1;
  n2=work2;
  r=work3;
  t12=work4;
  tdr=work5;
  tfr=work6;
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
/*=======================================================================
  Compute opacity, radiation energy gradient vector, and R
*/
  istart=N1BND;
  iend=N1BND+N1+1;
  LOOPRC3 {
    LOOPRC2 {
      temp_(&e[k][j][-N1BND],&s[1][k][j][-N1BND],&gam,&istart,&iend,&tt[-N1BND],&dtde[-N1BND]);
      absorp_(&tt[-N1BND],&s[1][k][j][-N1BND],&istart,&iend,&kap[k][j][-N1BND],&dkapdt[-N1BND]);
      scatt_(&tt[-N1BND],&s[1][k][j][-N1BND],&istart,&iend,&sig[k][j][-N1BND]);
    }
  }
  LOOPRC3 {
    LOOPRC2 {
      LOOPRC1 {
        qa = dx[2][1][i]/(er[k][j+1][i+1]+er[k][j+1][i-1]+er[k][j-1][i+1]+er[k][j-1][i-1]);
        der1 = (qa*(er[k][j][i+1]-er[k][j][i-1]))/ (dx[2][1][i]+dx[2][1][i+1]);
        der2 = (qa*(er[k][j+1][i]-er[k][j-1][i]))/((dx[2][2][j]+dx[2][2][j+1])*G2(2,i));
	dernorm = sqrt(der1*der1 + der2*der2);
	n1[k][j][i]= der1/(dernorm+SMALL);
	n2[k][j][i]= der2/(dernorm+SMALL);
	r[k][j][i] = dernorm/((kap[k][j][i]+sig[k][j][i])*0.25*dx[2][1][i]) + SMALL;
      }
    }
  }
/*  Compute FLD constant and components of Eddington tensor for 
  Chapman-Enskog theory (ifld=1)
*/
  if (ifld==1) {
    LOOPRC {    
	lmda = (2.0 + r[k][j][i])/(6.0 + 3.0*r[k][j][i]+ r[k][j][i]*r[k][j][i]);
	chi  = lmda + lmda*lmda*r[k][j][i]*r[k][j][i];
	f[1][1][k][j][i] = 0.5*((1.-chi) + (3.*chi-1.)*n1[k][j][i]*n1[k][j][i]);
	f[2][2][k][j][i] = 0.5*((1.-chi) + (3.*chi-1.)*n2[k][j][i]*n2[k][j][i]);
	f[1][2][k][j][i] = 0.5*(           (3.*chi-1.)*n1[k][j][i]*n2[k][j][i]);
	tdr[k][j][i] = CCON*lmda/(kap[k][j][i]+sig[k][j][i]); 
        tfr[k][j][i] = lmda ;
    }
  }

/*
  Compute FLD constant and components of Eddington tensor for
  piecewise-linear Minerbo theory (ifld=2)
*/
  if (ifld==2) {
    LOOPRC{
	  qa   = 2.0/(3.0 + sqrt(9.0+12.0*r[k][j][i]*r[k][j][i]));
	  qb   = 1.0/(1.0 + r[k][j][i] + sqrt(1.0+2.0*r[k][j][i]));
          if (r[k][j][i]<=1.5){
            lmda = qa;
          }else {
            lmda = qb;
          }
	  chi  = lmda + lmda*lmda*r[k][j][i]*r[k][j][i];
        f[1][1][k][j][i] = 0.5*((1.-chi) + (3.*chi-1.)*n1[k][j][i]*n1[k][j][i]);
        f[2][2][k][j][i] = 0.5*((1.-chi) + (3.*chi-1.)*n2[k][j][i]*n2[k][j][i]);
        f[1][2][k][j][i] = 0.5*(           (3.*chi-1.)*n1[k][j][i]*n2[k][j][i]);
        tdr[k][j][i] = CCON*lmda/(kap[k][j][i]+sig[k][j][i]);
        tfr[k][j][i] = lmda ;
    }	
  }
/*
  Average diffusion constant to face centers.
*/
  LOOPC3{
    LOOPC2{
      LOOPFLD1C{
    dr[1][k][j][i] = 0.5*(tdr[k][j][i-1] + tdr[k][j][i]);
    fr[1][k][j][i] = 0.5*(tfr[k][j][i-1] + tfr[k][j][i]);
      }
    }
  }
  LOOPC3{
    LOOPC1{
      LOOPFLD2C{
    dr[2][k][j][i] = 0.5*(tdr[k][j-1][i] + tdr[k][j][i]);
    fr[2][k][j][i] = 0.5*(tfr[k][j-1][i] + tfr[k][j][i]);
      }
    }
  }
#endif
}
