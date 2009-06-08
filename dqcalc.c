#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif

// not sure if dq is right in spc, need grad?->curvature terms?
#define DQCURVATURE 0
// 0: no curvature terms
// 1: curvature terms from gradient of vector


/* Van Leer formula on staggered grid(scalars on zone centers, vectors on zone faces):

scalar and non-directional vectors:

S[i]* = S[im1] + 0.5*dxb[i]*dq[im1]-0.5*v[i]*dt*dq[im1] if v[i]>0
S[i]* = S[i] - 0.5*dxb[i]*dq[i]-0.5*v[i]*dt*dq[i] if v[i]<0

directional vectors:

V[i]* = V[i] + 0.5*dxa[i]*dq[i]-0.5*((v[ip1]+v[i])/2.0)*dt*dq[i] if ((v[ip1]+v[i])/2.0) > 0
V[i]* = V[ip1] - 0.5*dxa[i]*dq[ip1]+0.5*((v[ip1]+v[i])/2.0)*dt*dq[ip1] if ((v[ip1]+v[i])/2.0) < 0

Note: dq has 2* in it which cancels with one .5 above



*/
#define MINMODB(a,b) ( (fabs(a)<fabs(b)) ? (a) : (b) )
#define MINMOD(a,b) ( ((a)*(b)<=0) ? 0.0 : MINMODB(a,b) )



// really returns 0.5*dq
// slope limiter
 // slope is limited for central region to help connect the slopes up for smoother solution
#define WOODDQOLD(DQ) \
	pr = Dqp*Dqm ;\
	if(pr<=0){\
	  DQ[k][j][i]=0.0;\
	}\
	else{\
	  if( (fabs(Dqm)<fabs(Dqp))&&(fabs(Dqm)<fabs(Dqc)) ){\
	    DQ[k][j][i]= Dqm;\
	  }\
	  else if(fabs(Dqp)<fabs(Dqc)){\
	    DQ[k][j][i] = Dqp;\
	  }\
	  else{\
	    DQ[k][j][i] = 0.5*Dqc;\
	  }\
	}


// really returns 0.5*dq
// Note that if choses Dqm/Dqp, that chooses twice slope -- i.e. steepens slope by factor of 2
#define WOODDQ(DQ) \
     DQ[k][j][i]=(0.5*MINMOD(Dqc,MINMOD(2.0*Dqm,2.0*Dqp)));


// really returns 0.5*dq
// reduced slope for monotonicity
// should really average Dqm and Dqp correctly for nonunigrid and adjust factor of 2 issue
#define VLDQOLD(DQ) \
      pr = Dqp*Dqm ;\
      if(pr>0.0){\
	DQ[k][j][i]= pr/(Dqm + Dqp);\
      }\
      else{\
	DQ[k][j][i]=0.0;\
      }

#define VLDQ(DQ) \
      pr = Dqp*Dqm ;\
      DQ[k][j][i]=((pr>0.0) ? pr/(Dqm + Dqp) : 0.0);




void dqx_calc(FTYPE (*var)[N2M][N1M],FTYPE (*dq)[N2M][N1M])
{
  register int i,j,k ;
  static FTYPE Dqp,Dqm,Dqc,pr;

  if(advint==2){
    LOOPHC{
      Dqm = (var[k][j][i] - var[k][j][im1])*OARC11(k,j,i); // real derivative
      Dqc = (var[k][j][ip1] - var[k][j][im1])/(dx[2][1][i]+dx[2][1][ip1]); // real derivative with all factors of 2 included
      Dqp = (var[k][j][ip1] - var[k][j][i])*OARC11(k,j,ip1); // real derivative
      WOODDQ(dq)
    }// end looph
  } // end if woodward slope
  else if(advint==1){
    LOOPHC{
      Dqp = (var[k][j][ip1] - var[k][j][i])*OARC11(k,j,ip1);
      Dqm = (var[k][j][i] - var[k][j][im1])*OARC11(k,j,i);
      VLDQ(dq)
    }
  }
  else if(advint==0){
    LOOPHC{
      dq[k][j][i]=0;
    }
  }
}

// x/r has no curvature terms
void dqvx_calc(int wcom, FTYPE (*var)[N3M][N2M][N1M],FTYPE (*dqv)[N3M][N2M][N1M])
{
  register int i;
  int j,k;
  static FTYPE Dqp,Dqm,Dqc,pr ;
  int vecwcom;
  int memwcom;

  // set real component of array to address correct location
  if(wcom<0){
    memwcom=0;
    vecwcom=-wcom;
  }
  else{
    memwcom=wcom;
    vecwcom=wcom;
  }

  

  if(advint==2){
    if(vecwcom==1){
      LOOPHC{
	Dqm = (var[memwcom][k][j][i] - var[memwcom][k][j][im1])*OARC21(k,j,im1);
	Dqc = (var[memwcom][k][j][ip1]-var[memwcom][k][j][im1])/(dx[1][1][i]+dx[1][1][im1]);
	Dqp = (var[memwcom][k][j][ip1] - var[memwcom][k][j][i])*OARC21(k,j,i);
	WOODDQ(dqv[vecwcom])
      }// end looph
    } // end if wcom==1
    else if(vecwcom==2){
      LOOPHC{
	Dqm = (var[memwcom][k][j][i] - var[memwcom][k][j][im1])*OARC31(k,j,i);
	Dqc = (var[memwcom][k][j][ip1]-var[memwcom][k][j][im1])/(dx[2][1][ip1]+dx[2][1][i]);
	Dqp = (var[memwcom][k][j][ip1] - var[memwcom][k][j][i])*OARC31(k,j,ip1);
	WOODDQ(dqv[vecwcom])
      }// end looph
    } // end if wcom==2
    else if(vecwcom==3){
      LOOPHC{
	Dqm = (var[memwcom][k][j][i] - var[memwcom][k][j][im1])*OARC11(k,j,i);
	Dqc = (var[memwcom][k][j][ip1] - var[memwcom][k][j][im1])/(dx[2][1][ip1]+dx[2][1][i]);
	Dqp = (var[memwcom][k][j][ip1] - var[memwcom][k][j][i])*OARC11(k,j,ip1);
	WOODDQ(dqv[vecwcom])
      }// end looph
    } // end if wcom==3
  } // end if woodward slope
  else if(advint==1){
    if(vecwcom==1){
      LOOPHC{
	Dqp = (var[memwcom][k][j][ip1] - var[memwcom][k][j][i])*OARC21(k,j,i);
	Dqm = (var[memwcom][k][j][i] - var[memwcom][k][j][im1])*OARC21(k,j,im1);
	VLDQ(dqv[vecwcom])
      }
    }
    else if(vecwcom==2){
      // OARC31=OARC11, so ok for scalar, v2 or v3
      LOOPHC{
	Dqp = (var[memwcom][k][j][ip1] - var[memwcom][k][j][i])*OARC31(k,j,ip1);
	Dqm = (var[memwcom][k][j][i] - var[memwcom][k][j][im1])*OARC31(k,j,i);
	VLDQ(dqv[vecwcom])
      }
    }
    else if(vecwcom==3){
      LOOPHC{
	Dqp = (var[memwcom][k][j][ip1] - var[memwcom][k][j][i])*OARC11(k,j,ip1);
	Dqm = (var[memwcom][k][j][i] - var[memwcom][k][j][im1])*OARC11(k,j,i);
	VLDQ(dqv[vecwcom])
      }
    }
  }// end if van leer
  else if(advint==0){
    LOOPHC{
      dqv[vecwcom][k][j][i]=0;
    }
  } // end if donor cell
}


void dqy_calc(FTYPE (*var)[N2M][N1M],FTYPE (*dq)[N2M][N1M])
{
  register int i,j,k ;
  static FTYPE Dqp,Dqm,Dqc,pr ;

  if(advint==2){
    LOOPHC{   // special for periodicx2special in sweepy(vy) (hence bound)
      Dqm = (var[k][j][i] - var[k][jm1][i])*OARC12(k,j,i);
      Dqc = (var[k][jp1][i]-var[k][jm1][i])/(G2(2,i)*(dx[2][2][jp1]+dx[2][2][j]));
      Dqp = (var[k][jp1][i] - var[k][j][i])*OARC12(k,jp1,i);
      WOODDQ(dq)
    }
  }
  else if(advint==1){   // special for periodicx2special in sweepy(vy) (hence bound)
    LOOPHC{
      Dqp = (var[k][jp1][i] - var[k][j][i])*OARC12(k,jp1,i);
      Dqm = (var[k][j][i] - var[k][jm1][i])*OARC12(k,j,i);
      VLDQ(dq)
    }
  }
  else if(advint==0){
    LOOPHC{
      dq[k][j][i]=0;
    }
  }
}

void dqvy_calc(int wcom,FTYPE (*var)[N3M][N2M][N1M],FTYPE (*dqv)[N3M][N2M][N1M])
{
  register int i,j,k ;
  static FTYPE Dqp,Dqm,Dqc,pr ;
  int vecwcom;
  int memwcom;

  // set real component of array to address correct location
  if(wcom<0){
    memwcom=0;
    vecwcom=-wcom;
  }
  else{
    memwcom=wcom;
    vecwcom=wcom;
  }



  if(advint==2){
    if(vecwcom==1){
      LOOPHC{
	Dqm = (var[memwcom][k][j][i] - var[memwcom][k][jm1][i])*OARC22(k,j,i);
	Dqc = (var[memwcom][k][jp1][i] - var[memwcom][k][jm1][i])/(G2(1,i)*(dx[2][2][jp1]+dx[2][2][j]));
	Dqp = (var[memwcom][k][jp1][i] - var[memwcom][k][j][i])*OARC22(k,jp1,i);
	
	WOODDQ(dqv[vecwcom])
      }
    }
    else if(vecwcom==2){
      LOOPHC{ // special for periodicx2special, hence bound of flux(2) in sweepy(vy) (see notes there and in global.h)
	Dqm = (var[memwcom][k][j][i] - var[memwcom][k][jm1][i])*OARC32(k,jm1,i);
	Dqc = (var[memwcom][k][jp1][i] - var[memwcom][k][jm1][i])/(G2(2,i)*(dx[1][2][j]+dx[1][2][jm1]));
	Dqp = (var[memwcom][k][jp1][i] - var[memwcom][k][j][i])*OARC32(k,j,i);
	
	WOODDQ(dqv[vecwcom])
      }
    }
    else if(vecwcom==3){
      LOOPHC{
	Dqm = (var[memwcom][k][j][i] - var[memwcom][k][jm1][i])*OARC12(k,j,i);
	Dqc = (var[memwcom][k][jp1][i] - var[memwcom][k][jm1][i])/(G2(2,i)*(dx[2][2][jp1]+dx[2][2][j]));
	Dqp = (var[memwcom][k][jp1][i] - var[memwcom][k][j][i])*OARC12(k,jp1,i);
	WOODDQ(dqv[vecwcom])
      }
    }
  }// end if woodward
  else if(advint==1){
    if(vecwcom==1){
      LOOPHC{
	Dqp = (var[memwcom][k][jp1][i] - var[memwcom][k][j][i])*OARC22(k,jp1,i);
	Dqm = (var[memwcom][k][j][i] - var[memwcom][k][jm1][i])*OARC22(k,j,i);
	VLDQ(dqv[vecwcom])
#if(DQCURVATURE&&(COORD>1))
	dqv[vecwcom][k][j][i]+=-v2tov1(var[2],k,j,i)/G2(1,i)*DG2(1,i);
#endif
      }
    }
    else if(vecwcom==2){  // special for periodicx2special in sweepy(vy)
      LOOPHC{
	Dqp = (var[memwcom][k][jp1][i] - var[memwcom][k][j][i])*OARC32(k,j,i);
	Dqm = (var[memwcom][k][j][i] - var[memwcom][k][jm1][i])*OARC32(k,jm1,i);
	VLDQ(dqv[vecwcom])
#if(DQCURVATURE&&(COORD>1))
	dqv[vecwcom][k][j][i]+=v1tov2(var[1],k,j,i)/G2(2,i)*DG2(2,i);
#endif
      }
    }
    else if(vecwcom==3){
      LOOPHC{
	Dqp = (var[memwcom][k][jp1][i] - var[memwcom][k][j][i])*OARC12(k,jp1,i);
	Dqm = (var[memwcom][k][j][i] - var[memwcom][k][jm1][i])*OARC12(k,j,i);
	VLDQ(dqv[vecwcom])
      }
    }
  }// end if vanleer
  else if(advint==0){
    LOOPHC{
      dqv[vecwcom][k][j][i]=0;
    }
  }// end if donor cell
}



void dqz_calc(FTYPE (*var)[N2M][N1M],FTYPE (*dq)[N2M][N1M])
{
  register int i,j,k ;
  static FTYPE Dqp,Dqm,Dqc,pr ;

  if(advint==2){
    LOOPHC{
      Dqm = (var[k][j][i] - var[km1][j][i])*OARC13(k,j,i)*ODX(2,3,k);
      Dqc = (var[kp1][j][i]-var[km1][j][i])*OARC13(k,j,i)*(ODX(2,3,kp1)+ODX(2,3,k));
      Dqp = (var[kp1][j][i] - var[k][j][i])*OARC13(kp1,j,i)*ODX(2,3,kp1);
      WOODDQ(dq)
    }
  }
  else if(advint==1){
    LOOPHC{
      Dqp = (var[kp1][j][i] - var[k][j][i])*OARC13(kp1,j,i)*ODX(2,3,kp1);
      Dqm = (var[k][j][i] - var[km1][j][i])*OARC13(k,j,i)*ODX(2,3,k);
      VLDQ(dq)
    }
  }
  else if(advint==0){
    LOOPHC{
      dq[k][j][i]=0;
    }
  }
}

void dqvz_calc(int wcom,FTYPE (*var)[N3M][N2M][N1M],FTYPE (*dqv)[N3M][N2M][N1M])
{
  register int i,j,k ;
  static FTYPE Dqp,Dqm,Dqc,pr ;
  int vecwcom;
  int memwcom;

  // set real component of array to address correct location
  if(wcom<0){
    memwcom=0;
    vecwcom=-wcom;
  }
  else{
    memwcom=wcom;
    vecwcom=wcom;
  }



  if(advint==2){
    if(vecwcom==1){
      LOOPHC{
	Dqm = (var[memwcom][k][j][i] - var[memwcom][km1][j][i])*OARC23(k,j,i)*ODX(2,3,k);
	Dqc = (var[memwcom][kp1][j][i] - var[memwcom][km1][j][i])*OARC23(k,j,i)*(ODX(2,3,kp1)+ODX(2,3,k));
	Dqp = (var[memwcom][kp1][j][i] - var[memwcom][k][j][i])*OARC23(kp1,j,i)*ODX(2,3,kp1);

	WOODDQ(dqv[vecwcom])
      }
    }
    else if(vecwcom==2){
      LOOPHC{
	Dqm = (var[memwcom][k][j][i] - var[memwcom][km1][j][i])*OARC33(k,j,i)*ODX(2,3,k);
	Dqc = (var[memwcom][kp1][j][i] - var[memwcom][km1][j][i])*OARC33(k,j,i)*(ODX(2,3,kp1)+ODX(2,3,k));
	Dqp = (var[memwcom][kp1][j][i] - var[memwcom][k][j][i])*OARC33(kp1,j,i)*ODX(2,3,kp1);

	WOODDQ(dqv[vecwcom])
      }
    }
    else if(vecwcom==3){
      LOOPHC{
	Dqm = (var[memwcom][k][j][i] - var[memwcom][km1][j][i])*OARC13(km1,j,i)*ODX(1,3,km1);
	Dqc = (var[memwcom][kp1][j][i] - var[memwcom][km1][j][i])*OARC13(k,j,i)*(ODX(1,3,k)+ODX(1,3,km1));
	Dqp = (var[memwcom][kp1][j][i] - var[memwcom][k][j][i])*OARC13(k,j,i)*ODX(1,3,k);

	WOODDQ(dqv[vecwcom])
      }
    }
  }// end if woodward
  else if(advint==1){
    if(vecwcom==1){
      LOOPHC{
	Dqp = (var[memwcom][kp1][j][i] - var[memwcom][k][j][i])*OARC23(kp1,j,i)*ODX(2,3,kp1);
	Dqm = (var[memwcom][k][j][i] - var[memwcom][km1][j][i])*OARC23(k,j,i)*ODX(2,3,k);
	VLDQ(dqv[vecwcom])
#if(DQCURVATURE&&(COORD>1))
	dqv[vecwcom][k][j][i]+=-v3tov1(var[3],k,j,i)/G3(1,i)*DG3(1,i);
#endif
      }
    }
    else if(vecwcom==2){
      LOOPHC{
	Dqp = (var[memwcom][kp1][j][i] - var[memwcom][k][j][i])*OARC33(kp1,j,i)*ODX(2,3,kp1);
	Dqm = (var[memwcom][k][j][i] - var[memwcom][km1][j][i])*OARC33(k,j,i)*ODX(2,3,k);
	VLDQ(dqv[vecwcom])
#if(DQCURVATURE&&(COORD>1))
	dqv[vecwcom][k][j][i]+=-v3tov2(var[3],k,j,i)*OG4(1,j)*DG4(1,j)/(G2(2,i));
#endif

      }
    }
    else if(vecwcom==3){
      LOOPHC{
	Dqp = (var[memwcom][kp1][j][i] - var[memwcom][k][j][i])*OARC33(k,j,i)*ODX(1,3,k);
	Dqm = (var[memwcom][k][j][i] - var[memwcom][km1][j][i])*OARC33(km1,j,i)*ODX(1,3,km1);
	VLDQ(dqv[vecwcom])
#if(DQCURVATURE&&(COORD>1))
	dqv[vecwcom][k][j][i]+=v1tov3(var[1],k,j,i)/G3(2,i)*DG3(2,i)+v2tov3(var[2],k,j,i)*OG4(2,j)*DG4(2,j)/(G2(2,i));
#endif
      }
    }
  }// end if vanleer
  else if(advint==0){
    LOOPHC{
      dqv[vecwcom][k][j][i]=0;
    }
  }// end if donor cell
}
