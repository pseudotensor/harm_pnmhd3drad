#include <stdio.h>



#define im1(i) i-1
#define ip1 i+1

#define fun(i) (im1(i))

int main(void)
{

  int i,ii;
  int result;

  i=5;

  result=fun(i+1);
  printf("%d\n",result);

  i=6;
  result=fun(i);
  printf("%d\n",result);

  return(0);
}
