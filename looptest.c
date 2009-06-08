#include <math.h>
#include <stdio.h>

#define N 32


#define LOOPSUPERGEN for(temptempi=0,i=indx[0],j=indx[1],k=indx[2];temptempi<numiter;temptempi++,i=indx[temptempi*3],j=indx[temptempi*3+1],k=indx[temptempi*3+2])

int main(void)
{
  int temptempi,i,j,k;
  int indx[N*3+1];

  numiter=50;
  for(k=0;k<N;k++)  for(j=0;j<N;j++)  for(i=0;i<N;i++){

    indx[k*3
  }

}
