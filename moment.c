#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif


//=======================================================================
//////////////////////////  SUBROUTINE MOMENT  \\\\\\\\\\\\\\\\\\\\\\\\

void moment(void)
{
#if(RAD)
/*
  PURPOSE:  Controls the update of all the dynamical variables (e,er,
  fr1,fr2) due to source terms in the moment equations. This requires
  solving the dynamical (moment) equations including (only) the
  source/sink terms. We use Newton-Raphson iteration of the coupled,
  implicit equations to update e and er, and then use the algebraic
  automatic flux limiting (AFL) relation to update fr1 and fr2.

  Note that the full solution of the radiation dynamical equations also
  requires including the transport terms.  This is done in the routines
  TRANX* for er and MOMX* for fr1 and fr2, which are executed in the
  transport step.  Since fr1 and fr2 are vectors, there are additional
  "source" terms due to derivatives of the unit vectors added in the
  AFL routine (FLUX).

  EXTERNALS:

  LOCALS:
-----------------------------------------------------------------------
*/
  int boost;
  int iter,miticcg,i,j,k,di,dj,NTOT,istart,iend,n1mm;
  FTYPE one,max;
  FTYPE dthydro,trad,eiccg,epsemx,epsermx,q1
         ,demx,dermx,enorm,ernorm;
  static FTYPE *tt,*dtde,*dbbdt,*dkapdt;
  static FTYPE (*detot)[N2M][N1M],(*dertot)[N2M][N1M];
  static FTYPE *derel,*derrel,*derelj,*derrelj;

/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
=======================================================================
  implicit solution of the coupled moment equations
  Initial guesses for radiation and material energy densities are 
  current values of these variables.
    line 1 = target if timestep reduced
    line 2 = target if timestep not yet finished, timestep not reduced
*/
  tt=work1x;
  dtde=work2x;
  dbbdt=work3x;
  dkapdt=work4x;
  derel=work5x;
  derrel=work6x;
  derelj=work1y;
  derrelj=work2y;
  detot=work1;
  dertot=work2;
  di=1;
  dj=1;
  n1mm=N1M;
  one = 1.0;
  nred = 0;
  dthydro = dt;
L1:    
      trad = 0.0;
L2:     
  boost = 0;
  LOOPC{
          detot[k][j][i] = 0.0;
          dertot[k][j][i] = 0.0;
  }
  LOOPFC {
          ern[k][j][i] = er[k][j][i] + der[k][j][i] ;
          en[k][j][i]  = e [k][j][i]  + de[k][j][i] ;
  }
  printf("e de %g %g\n",e[0][1][1],de[0][1][1]);
/*
  Compute gas pressure and absorption coefficients at old and new
   material energy density (which are initially identical), and
   derivatives (evaluated at new material energy density)
*/
  istart=N1BND+1;
  iend=N1BND+N1;
  LOOPC3 {
    LOOPC2 {
      temp_(&e[k][j][-N1BND],&s[1][k][j][-N1BND],&gam,&istart,&iend,&tt[-N1BND],&dtde[-N1BND]);
      eos_(&e[k][j][-N1BND],&s[1][k][j][-N1BND],&gam,&istart,&iend,&pre[k][j][-N1BND],&dpde[k][j][-N1BND]);  
      planck_(&tt[-N1BND],&istart,&iend,&bb[k][j][-N1BND],&dbbdt[-N1BND]);
      absorp_(&tt[-N1BND],&s[1][k][j][-N1BND],&istart,&iend,&kap[k][j][-N1BND],&dkapdt[-N1BND]);      
      scopy_(&n1mm,&pre[k][j][-N1BND],&di,&pren[k][j][-N1BND],&dj);
      scopy_(&n1mm,&bb[k][j][-N1BND],&di,&bbn[k][j][-N1BND],&dj);
      scopy_(&n1mm,&kap[k][j][-N1BND],&di,&kapn[k][j][-N1BND],&dj);
      LOOPC1 {
          dbbde[k][j][i]  = dbbdt[i]*dtde[i];
          dkapde[k][j][i] = dkapdt[i]*dtde[i];
      }
    }
  }

//---------------  Start Newton-Raphson iteration loop  -----------------
  for(iter=1;iter<=nmeiter;iter++){
      if(0>=(3-iter)){max=0;}
	if((3-iter)>0){max=3-iter;}
	eiccg   = epsrad*pow(10.0,max);
	miticcg = maxrad;
	riccg(&eiccg,&ks0rad,&miticcg);
      if (miticcg>=maxrad){
          printf("ICCGAF did not converge with dt= %g  epsrad  eiccg = %g  %g \n",dt,epsrad,eiccg);
          goto L1001;
      }  
      else{
          if (iorad==1){ printf(" nhy=, %d, MOMENT: ICCG  converged in ,  %d, iterations, eps=, %g \n",
         nhy,miticcg,eiccg);}
      }
//  restrict maximum change in a variable by q1
//  sum changes in e and er into running totals
      epsemx = 0.0;
	epsermx = 0.0;
      enorm = 0.0;
      ernorm = 0.0;
      LOOPC{
            enorm = enorm + fabs(en[k][j][i]);
            ernorm = ernorm + fabs(ern[k][j][i]);
      }
      LOOPC{
            detot[k][j][i] = detot[k][j][i] + de[k][j][i];
            dertot[k][j][i] = dertot[k][j][i] + der[k][j][i] ;
		if(fabs(de[k][j][i]/enorm)>=epsemx){
		  epsemx=fabs(de[k][j][i]/enorm);}
		if(fabs(der[k][j][i]/ernorm)>=epsermx){
		  epsermx=fabs(der[k][j][i]/ernorm);}
	}

      if(demax/epsemx<=dermax/epsermx){
	  q1= demax/epsemx;
	}else{
	  q1= dermax/epsermx;
      }
	if(q1>one){q1=one;}
//  apply iterates to update variables // the energy change amount can not reachdemax or dermax of the total energy
      LOOPC{
            en[k][j][i] = en[k][j][i] + q1*de[k][j][i];
            ern[k][j][i] = ern[k][j][i] + q1*der[k][j][i];
      }
      bvalern();
#if(RADTEST)
      fprintf(rad_file,"dt= %g, boost=%d, epsemx epsermx<epsme=  %g %g %g \n",boost,dt,epsemx,epsermx,epsme);
#endif
//  Apply boost iteration, exit loop if iterative solution has converged

      if (epsemx<=epsme&&epsermx<epsme){
        if (boost==1) {goto L1002;}
          boost =1;
      }

//  Else, prepare for another iteration

      LOOPC3{
        LOOPC2{
	  temp_(&en[k][j][-N1BND],&s[1][k][j][-N1BND],&gam,&istart,&iend,&tt[-N1BND],&dtde[-N1BND]);
          eos_(&en[k][j][-N1BND],&s[1][k][j][-N1BND],&gam,&istart,&iend,&pren[k][j][-N1BND],&dpde[k][j][-N1BND]);  
          planck_(&tt[-N1BND],&istart,&iend,&bbn[k][j][-N1BND],&dbbdt[-N1BND]);
          absorp_(&tt[-N1BND],&s[1][k][j][-N1BND],&istart,&iend,&kapn[k][j][-N1BND],&dkapdt[-N1BND]);      

          LOOPC1{
            dbbde[k][j][i] = dbbdt[i]*dtde[i];
            dkapde[k][j][i]  = dkapdt[i]*dtde[i];
	    }
        }
      }
//      printf("OK %d %g \n",iter,en[0][1][1]);

  }
//------------------  End Newton-Raphson iteration loop  ----------------
/*
  If this point is reached, then Newton-Raphson did not converge within
  the maximum allowed number of iterations.  Try reducing timestep.
    1001 = target if total correction exceeds limit or ICCG failed
           to converge => reduce timestep 
*/
      printf("**********  NR failed to converge with dt=,%g , **********,(iter,de,der) = %d %g %g",dt,iter
       ,epsemx,epsermx);
L1001:
      nred = nred+1;
      if (nred > 6){
        printf("**********  TIMESTEP REDUCTION FAILED:  ABORTING  **********");
	myexit(1);
        return;
      }
      dt = dt/4.0;
      goto L1;
/*
  NR apparently converged.  Check that total changes in variables are
  less than max allowed
*/
L1002:  
      if (iorad==1 ){ printf("NR iteration, %d , finished: (q1,de,der)=, %g , %g, %g", iter,q1,epsemx,epsermx);}
      demx = 0.0;
      dermx = 0.0;
      enorm = 0.0;
      ernorm = 0.0;
      LOOPC{
           enorm = enorm + fabs(e[k][j][i] );
           ernorm = ernorm + fabs(er[k][j][i] );
      }
      LOOPC{
	  if(demx<=abs(detot[k][j][i]/enorm)){demx=abs(detot[k][j][i]/enorm);}
	  if(dermx<=abs(dertot[k][j][i]/ernorm)){dermx=abs(dertot[k][j][i]/ernorm);}
      }
//  total changes too large, reduce timestep
      if (demx > dtotmax || dermx > dtotmax) {
        printf("**********  Total changes too large with dt=, %g ********** (demx,dermx) = %g, %g", dt,demx,dermx);
        goto L1001;
      }

//  final energy density update

	NTOT=N1M*N2M*1;
      scopy_(&NTOT,&ern[0][-N2BND][-N1BND],&di,&er[0][-N2BND][-N1BND],&dj);
      scopy_(&NTOT,&en[0][-N2BND][-N1BND],&di,&e[0][-N2BND][-N1BND],&dj);
      bvaler();
      bvale();

//  Check that entire timestep is finished
      printf("dt %g trad in radiation %g, t %g thydro %g e  %g  de %g\n",dt,trad,t,dthydro,e[0][1][1],de[0][1][1]);
      trad = trad + dt;
      if (trad < dthydro) {
        if(dt>(dthydro-trad)){dt=dthydro-trad;}
        goto L2;
      }
      dt = dthydro;
      printf("nhy= %d NR converged with %d dt reductions (de,der)=,%e ,%e ",nhy,nred,demx,dermx);
      radend();
#endif
}
