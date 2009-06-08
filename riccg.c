#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif

 FTYPE av0b[2*N1M*N2M+1],av1b[2*N1M*N2M+1],bv0b[2*N1M*N2M+1],xvb[2*N1M*N2M+1],yvb[2*N1M*N2M+1],erry[2*N1M*N2M+1] ;
 FTYPE bv1b[2*N1M*N2M+1],bm1b[2*N1M*N2M+1],workb[4*N1M*N2M+1],epst,ks0t,err;
//=======================================================================
////////////////////////  SUBROUTINE RICCG  \\\\\\\\\\\\\\\\\\\\\\\\\\\

  void riccg(FTYPE *epss,int *ks00,int *maxitt)
{
#if(RAD)
/*
  PURPOSE:

  INPUT ARGUMENTS:
      eps   = minimum L2 error for ICCGAF
      ks0   = max number of cyclic reductions in ICCGAF
      maxit = max number of cg iterations     in ICCGAF

  OUTPUT ARGUMENTS:
      eps   = actual L2 error achieved         in ICCGAF
      maxit = actual number of cg iterations used "   "

  EXTERNALS: ICCGAF

  LOCALS:
  av0,av1,bv0   = vector arrays of matrix diagonal bands. See ICCGAF
    bv1,bm1       documentation. In our application, bv1=bm1=0
  xv            = vector of potential values (solution vector)
  yv            = RHS of matrix eqn
  work          = work space used by ICCGAF
  qa...qe       = dummy variables
-----------------------------------------------------------------------
*/
  int nz,i,j,k,m,n11,n22,ntest;
  FTYPE ttp;

  static FTYPE qa,qb,qc,qd,qe
        ,*av0,*av1,*bv0, *bv1
        ,*bm1,*xv,*yv,*work;

  av0=(FTYPE (*))(&(wa[0][-1]));
  av1=(FTYPE (*))(&(wc[0][-1]));
  bv0=(FTYPE (*))(&(wcg1[-1]));
  bv1=(FTYPE (*))(&(wcg1[2*N1M*N2M-1]));
  bm1=(FTYPE (*))(&(wcg1[4*N1M*N2M-1]));
  xv=(FTYPE (*))(&(wcg1[6*N1M*N2M-1]));
  yv=(FTYPE (*))(&(wcg1[8*N1M*N2M-1]));
  work=(FTYPE (*))(&(wcg1[10*N1M*N2M-1]));
 

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////////////
/*=======================================================================
  Set up matrix elements
*/
  derivs();
  rhs();
  LOOPC {
    m  = (j-0)*N1 + i - 0 + 1 ;
    qa = G2(1,i+1)*G3(1,i+1)/dx[2][1][i+1]*dr[1][k][j][i+1];
    qb = G2(1,i  )*G3(1,i  )/dx[2][1][i]*dr[1][k][j][i];
    qc = G2(2,i)*G2(2,i) ;
    qd = G4(1,j+1)/dx[2][2][j+1]*dr[2][k][j+1][i];
    qe = G4(1,j  )/dx[2][2][j  ]*dr[2][k][j][i] ;


    av0[m] = -dt*(dvl[1][1][i]/qc*(qd+qe) + dvl[1][2][j]*(qa+qb));
    av0[m] = av0[m] - dvl[1][1][i]*dvl[1][2][j]*(dfn1der[k][j][i]
           - dfn1de[k][j][i]*dfn2der[k][j][i]/dfn2de[k][j][i]);
    bv0[m] = dt*dvl[1][1][i]*qd/qc;
    av1[m] = dt*dvl[1][2][j]*qa;
    bv1[m] = 0.0;
    bm1[m] = 0.0;

    xv [m] = der[k][j][i];
    yv [m] = -dvl[1][1][i]*dvl[1][2][j]
        *(dfn1de[k][j][i]*fn2[k][j][i]/dfn2de[k][j][i]-fn1[k][j][i]);
    
  }
/*  Add boundary conditions to matrix elements or yv, as appropriate.
  Order is iib, oib, ijb, ojb. 
*/
  qa = dt*G2(1,0)*G3(1,0)/dx[2][1][0];
  LOOPC3 {
    LOOPC2 {
        m = (j-0)*N1 + 1;
        qb = qa*dr[1][k][j][0];
        if (liib[j]==1||liib[j]==2) {av0[m]=av0[m]+qb*dvl[1][2][j];}
    }
  }

  qa = dt*G2(1,N1)*G3(1,N1)/dx[2][1][N1];
  LOOPC3 {
    LOOPC2 {
        m = (j-0+1)*N1;
        av1[m] = 0.0;
        qb = qa*dr[1][k][j][N1];
        if (loib[j]==1||loib[j]==2) {av0[m]=av0[m]+qb*dvl[1][2][j];}
    }
  }

  qa = dt*G4(1,0)/dx[2][2][0];
  LOOPC3{
    LOOPC1{
        m = i - 0 + 1;
        qb = dr[2][k][0][i]*dvl[1][1][i]/(G2(2,i)*G2(2,i));
        if (lijb[i]==1||lijb[i]==2) av0[m] = av0[m] + qa*qb;
    }
  }

  qa = dt*G4(1,N2)/dx[2][2][N2];
  LOOPC3 {
    LOOPC1 {
        m = (N2-1)*N1 + i - 0 + 1;
        bv0[m] = 0.0;
        qb = dr[2][k][N2][i]*dvl[1][1][i]/(G2(2,i)*G2(2,i));
        if (lojb[i]==1||lojb[i]==2) {av0[m] = av0[m] + qa*qb;}
    }
  }

/*
  Matrix elements are finished, solve system.

      *epss = *epss/sasum(nx1z*nx2z,yv,1)
*/
  nz = N1 * N2;
  n11=N1;
  n22=N2;
#if(RADTEST)
  ntest=10;
  fprintf(rad_file,"\n av0b");
  for(i=1;i<=ntest;i++){
    fprintf(rad_file," %g",av0b[i]);
  }
  fprintf(rad_file,"\n av1b");
  for(i=1;i<=ntest;i++){
    fprintf(rad_file," %g",av1b[i]);
  }
  fprintf(rad_file,"\n bv0b");
  for(i=1;i<=ntest;i++){
    fprintf(rad_file," %g",bv0b[i]);
  }
  fprintf(rad_file,"\n yvb");
  for(i=3170;i<=ntest+3170;i++){
    fprintf(rad_file," %g",yvb[i]);
  }
/*  fprintf(rad_file,"\n fn2");
  for(i=1;i<=ntest;i++){
    ttp=fn2[0][0][i];
    fprintf(rad_file," %g",ttp);
  }
  fprintf(rad_file,"\n fn1");
  for(i=1;i<=ntest;i++){
    fprintf(rad_file," %g",fn1[0][0][i]);
  }*/
#endif
  iccgaf_(&n11,&n22,epss,ks00,maxitt,&xv[nz+1],&yv[nz+1],&av0[1],
                 &work[1],&work[1],&work[1],&work[nz+1],&work[nz+1],
                 &work[2*nz+1],&work[3*nz+1],
                       &av0[1],&av1[1],&bv0[1],&bv1[1],&bm1[1],&xv[1],&yv[1],&work[1]);

#if(RADTEST)
/*  fprintf(rad_file,"\n xv");
  for(i=nz+1;i<=nz+ntest;i++){
    fprintf(rad_file," %g",xv[i]);
  }
*/
   *maxitt=200;
   *epss=1.e-20;
   *ks00=99;
   iccgaf_(&n11,&n22,epss,ks00,maxitt,&xvb[nz+1],&yvb[nz+1],&av0[1],
                 &work[1],&work[1],&work[1],&work[nz+1],&work[nz+1],
                 &work[2*nz+1],&work[3*nz+1],
                       &av0[1],&av1[1],&bv0[1],&bv1[1],&bm1[1],&xvb[1],&yvb[1],&work[1]);
   fprintf(rad_file," \n xvb");
   for(i=1;i<=ntest;i++){
     fprintf(rad_file," %g",xvb[i]);
   }
    fprintf(rad_file,",maxit, %d",*maxitt);
   *maxitt=200;
   *epss=1.e-20;
   *ks00=99;
   iccgaf_(&n11,&n22,epss,ks00,maxitt,&xv[nz+1],&yv[nz+1],&av0[1],
                 &work[1],&work[1],&work[1],&work[nz+1],&work[nz+1],
                 &work[2*nz+1],&work[3*nz+1],
                       &av0[1],&av1[1],&bv0[1],&bv1[1],&bm1[1],&xv[1],&yv[1],&work[1]);
  err=0.;
  LOOPC {
    m  = (j-0)*N1 + i - 0 + 1 ;
    av0b[m]=av0[m];
    bv0b[m]=bv0[m];
    av1b[m]=av1[m];
    xvb[m]=xv[m];
    erry[m]=(yv[m]-yvb[m])/yv[m];
    if(fabs(erry[m])>err)err=fabs(erry[m]);
    yvb[m]=yv[m]/2.2864;
  }

    fprintf(rad_file,"\n xv2");
  for(i=1;i<=ntest;i++){
    fprintf(rad_file," %g",xv[i]);
  }

  fprintf(rad_file,",maxit2, %d\n",*maxitt);
  fprintf(rad_file," av0");
  for(i=1;i<=ntest;i++){
    fprintf(rad_file," %g",av0[i]);
  }
  fprintf(rad_file,"\n av1");
  for(i=1;i<=ntest;i++){
    fprintf(rad_file," %g",av1[i]);
  }
  fprintf(rad_file,"\n bv0");
  for(i=1;i<=ntest;i++){
    fprintf(rad_file," %g",bv0[i]);
  }
  fprintf(rad_file,"\n yv");
  for(i=3170;i<=ntest+3170;i++){
    fprintf(rad_file," %g",yv[i]);
  }
   fprintf(rad_file,"\n err");
    fprintf(rad_file," %g\n",err);
#endif
/*
  deconvolve xv values to get der(i,j), compute de(i,j)
*/
  LOOPC {
    der[k][j][i] = xv[(j-0)*N1 + i - 0 + 1];
    de [k][j][i] = (-fn2[k][j][i] -dfn2der[k][j][i]*der[k][j][i] )/dfn2de[k][j][i] ;
  }
/* if(*maxitt>=200){
  LOOPC {
    der[k][j][i]=0.;
    de[k][j][i]=0.;
  }
 }
*/
#endif
}
