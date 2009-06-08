#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif

/*
=======================================================================
//////////////////////////  SUBROUTINES BVAL*  \\\\\\\\\\\\\\\\\\\\\\\\

  PURPOSE:  These routines set the boundary values for the radiation
  energy density at the old and new timestep needed during the NR
  iterations of the moment equations.  The same boundary tyes as
  implemented for the hydrodynamical variables are implemented here.

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
=======================================================================

---------------  radiation energy density boundary values  ------------
*/
void bvaler(void)
{
#if(RAD)
  int i,j,k;
  erfloor=SMALL;
//  inner i boundary 
  LOOPC2{
        if (liib[j]==1)  er[0][j][-1] = er[0][j][0] ;
        if (liib[j]==2)  er[0][j][-1] = er[0][j][0] ;
        if (liib[j]==3)  er[0][j][-1] = erfloor ;
        if (liib[j]==4)  er[0][j][-1] = er[0][j][N1-1] ;
                         er[0][j][-2] = er[0][j][-1];
        if (liib[j]==1)  er[0][j][-2] = er[0][j][1];
        if (liib[j]==3)  er[0][j][-2] = erfloor ;
        if (liib[j]==4)  er[0][j][-2] = er[0][j][N1-2];
  }

//  outer i boundary

  LOOPC2{
        if (loib[j]==1)  er[0][j][N1] = er[0][j][N1-1];
        if (loib[j]==2)  er[0][j][N1] = er[0][j][N1-1];
        if (loib[j]==3)  er[0][j][N1] = erfloor ;
        if (loib[j]==4)  er[0][j][N1] = er[0][j][0];
                         er[0][j][N1+1] = er[0][j][N1];
        if (loib[j]==1)  er[0][j][N1+1] = er[0][j][N1-2];
        if (loib[j]==3)  er[0][j][N1+1] = erfloor ;
        if (loib[j]==4)  er[0][j][N1+1] = er[0][j][1];
  }
//  inner j boundary 

   LOOPF1{
        if (lijb[i]==1)  er[0][-1][i] = er[0][0][i];
        if (lijb[i]==2)  er[0][-1][i] = er[0][0][i];
        if (lijb[i]==3)  er[0][-1][i] = erfloor ;
        if (lijb[i]==4)  er[0][-1][i] = er[0][N2-1][i];
                         er[0][-2][i] = er[0][-1][i];
        if (lijb[i]==1)  er[0][-2][i] = er[0][1][i];
        if (lijb[i]==3)  er[0][-2][i] = erfloor ;
        if (lijb[i]==4)  er[0][-2][i] = er[0][N2-2][i];
   }

//  outer j boundary

   LOOPF1{
        if (lojb[i]==1)  er[0][N2][i] = er[0][N2-1][i];
        if (lojb[i]==2)  er[0][N2][i] = er[0][N2-1][i];
        if (lojb[i]==3)  er[0][N2][i] = erfloor;
        if (lojb[i]==4)  er[0][N2][i] = er[0][0][i];
                         er[0][N2+1][i] = er[0][N2][i];
        if (lojb[i]==1)  er[0][N2+1][i] = er[0][N2-2][i];
        if (lojb[i]==3)  er[0][N2+1][i] = erfloor;
        if (lojb[i]==4)  er[0][N2+1][i] = er[0][1][i];
   }
#endif
} 
/*
------------  new radiation energy density boundary values  -----------
*/
    void bvalern()
{
#if(RAD)
  int i,j,k;
  erfloor=SMALL;
//  inner i boundary 
  LOOPC2{
        if (liib[j]==1)  ern[0][j][-1] = ern[0][j][0] ;
        if (liib[j]==2)  ern[0][j][-1] = ern[0][j][0] ;
        if (liib[j]==3)  ern[0][j][-1] = erfloor ;
        if (liib[j]==4)  ern[0][j][-1] = ern[0][j][N1-1] ;
                         ern[0][j][-2] = ern[0][j][-1];
        if (liib[j]==1)  ern[0][j][-2] = ern[0][j][1];
        if (liib[j]==3)  ern[0][j][-2] = erfloor ;
        if (liib[j]==4)  ern[0][j][-2] = ern[0][j][N1-2];
  }

//  outer i boundary

  LOOPC2{
        if (loib[j]==1)  ern[0][j][N1] = ern[0][j][N1-1];
        if (loib[j]==2)  ern[0][j][N1] = ern[0][j][N1-1];
        if (loib[j]==3)  ern[0][j][N1] = erfloor ;
        if (loib[j]==4)  ern[0][j][N1] = ern[0][j][0];
                         ern[0][j][N1+1] = ern[0][j][N1];
        if (loib[j]==1)  ern[0][j][N1+1] = ern[0][j][N1-2];
        if (loib[j]==3)  ern[0][j][N1+1] = erfloor ;
        if (loib[j]==4)  ern[0][j][N1+1] = ern[0][j][1];
  }
//  inner j boundary 

   LOOPF1{
        if (lijb[i]==1)  ern[0][-1][i] = ern[0][0][i];
        if (lijb[i]==2)  ern[0][-1][i] = ern[0][0][i];
        if (lijb[i]==3)  ern[0][-1][i] = erfloor ;
        if (lijb[i]==4)  ern[0][-1][i] = ern[0][N2-1][i];
                         ern[0][-2][i] = ern[0][-1][i];
        if (lijb[i]==1)  ern[0][-2][i] = ern[0][1][i];
        if (lijb[i]==3)  ern[0][-2][i] = erfloor ;
        if (lijb[i]==4)  ern[0][-2][i] = ern[0][N2-2][i];
   }

//  outer j boundary

   LOOPF1{
        if (lojb[i]==1)  ern[0][N2][i] = ern[0][N2-1][i];
        if (lojb[i]==2)  ern[0][N2][i] = ern[0][N2-1][i];
        if (lojb[i]==3)  ern[0][N2][i] = erfloor;
        if (lojb[i]==4)  ern[0][N2][i] = ern[0][0][i];
                         ern[0][N2+1][i] = ern[0][N2][i];
        if (lojb[i]==1)  ern[0][N2+1][i] = ern[0][N2-2][i];
        if (lojb[i]==3)  ern[0][N2+1][i] = erfloor;
        if (lojb[i]==4)  ern[0][N2+1][i] = ern[0][1][i];
   }
#endif
} 
////////////////////////// bvale /////////////
/*
nflo =  1  for reflecting
           =  2  for flow out
           =  3  for flow in
           =  4  for periodic
*/
void bvale(void)
{
#if(RAD)
  int i,j,k;
  erfloor=SMALL;
//  inner i boundary 
  LOOPC2{
        if (niib[j]==1)  e[0][j][-1] = e[0][j][0] ;
        if (niib[j]==2)  e[0][j][-1] = e[0][j][0] ;
        if (niib[j]==3)  e[0][j][-1] = erfloor ;
        if (niib[j]==4)  e[0][j][-1] = e[0][j][N1-1] ;
        if (niib[j]==1)  e[0][j][-2] = e[0][j][1];
	  if (niib[j]==2)  e[0][j][-2] = e[0][j][1] ;
        if (niib[j]==3)  e[0][j][-2] = erfloor ;
        if (niib[j]==4)  e[0][j][-2] = e[0][j][N1-2];
  }

//  outer i boundary

  LOOPC2{
        if (noib[j]==1)  e[0][j][N1] = e[0][j][N1-1];
        if (noib[j]==2)  e[0][j][N1] = e[0][j][N1-1];
        if (noib[j]==3)  e[0][j][N1] = erfloor ;
        if (noib[j]==4)  e[0][j][N1] = e[0][j][0];
        if (noib[j]==1)  e[0][j][N1+1] = e[0][j][N1-2];
	  if (noib[j]==2)  e[0][j][N1+1] = e[0][j][N1-2];
        if (noib[j]==3)  e[0][j][N1+1] = erfloor ;
        if (noib[j]==4)  e[0][j][N1+1] = e[0][j][1];
  }
//  inner j boundary 

   LOOPF1{
        if (nijb[i]==1)  e[0][-1][i] = e[0][0][i];
        if (nijb[i]==2)  e[0][-1][i] = e[0][0][i];
        if (nijb[i]==3)  e[0][-1][i] = erfloor ;
        if (nijb[i]==4)  e[0][-1][i] = e[0][N2-1][i];
        if (nijb[i]==1)  e[0][-2][i] = e[0][1][i];
	  if (nijb[i]==2)  e[0][-2][i] = e[0][1][i];
        if (nijb[i]==3)  e[0][-2][i] = erfloor ;
        if (nijb[i]==4)  e[0][-2][i] = e[0][N2-2][i];
   }

//  outer j boundary

   LOOPF1{
        if (nojb[i]==1)  e[0][N2][i] = e[0][N2-1][i];
        if (nojb[i]==2)  e[0][N2][i] = e[0][N2-1][i];
        if (nojb[i]==3)  e[0][N2][i] = erfloor;
        if (nojb[i]==4)  e[0][N2][i] = e[0][0][i];
        if (nojb[i]==1)  e[0][N2+1][i] = e[0][N2-2][i];
        if (nojb[i]==2)  e[0][N2+1][i] = e[0][N2-2][i];
        if (nojb[i]==3)  e[0][N2+1][i] = erfloor;
        if (nojb[i]==4)  e[0][N2+1][i] = e[0][1][i];
   }
#endif
} 
