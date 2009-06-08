#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif

void initrad(void)
{
  int i,j,k;
  erfloor=SMALL;
  ifld=1;
// maximum number of NR iteration allowed
  nmeiter=200; 
// error limit for ICCGAF
  epsrad=1.0e-10;
// level of cyclic reduction in ICCGAF
  ks0rad=99 ;
//maximum iteration count for ICCGAF
  maxrad=2000 ;
// iorad=1 to dump convergence data for NR iteration
  iorad=0;
// error criteria for NR convergence
  epsme=1.0e-5;
// max relative change allowed in e, er in a single NR iter
  demax=0.2;
  dermax=0.2;
// max total relative change allowed over all NR iterations
  dtotmax=0.2;
//constant
  radth=1.;
  CCON=2.998e10/lunit*tunit;
  nred=0;

  LOOPFC{
    er[k][j][i]=pow(sanal[2][k][j][i]/sanal[1][k][j][i]/rhounit/kcontcgs*mmw*mhydrcgs*(gam-1.)*enedenunit,4.)*sigmacontnew*4./CCON;
    der[k][j][i]=1.0e-6*er[k][j][i];
    de[k][j][i]=1.0e-6*sanal[2][k][j][i];
  }
// initialize boundary condition
  LOOPF2{
    liib[j]=1;
    loib[j]=1;
    niib[j]=1;
    noib[j]=1;
  }
  LOOPF1{
    lijb[i]=1;
    lojb[i]=1;
    nijb[i]=1;
    nojb[i]=1;
  }

}
