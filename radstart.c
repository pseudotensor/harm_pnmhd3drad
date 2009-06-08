#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif


void radstart(void)
{
#if(RAD)
  int i,j,k;
  ifld=1;
  LOOPFC {
    e[k][j][i]=s[2][k][j][i];
//change v ,g,dx too, probably just change c and t
    pre[k][j][i]=(gam-1.)*e[k][j][i];
  }
#endif
}


void radend(void)
{
#if(RAD)
  int i,j,k;
  LOOPFC{
    s[2][k][j][i]=e[k][j][i];
  }
#endif
}
