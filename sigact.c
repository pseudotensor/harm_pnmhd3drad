#include "global.h"
#include "defs.h"

FTYPE r2,the2;
int pflag,twolayer,N2tem;
void sigact(void);
FTYPE siginte(int iinit, int jinit);
FTYPE gridl(FTYPE xini,int i, int j, FTYPE r1, FTYPE the1);
FTYPE gridld(FTYPE xini,int i, int j, FTYPE r1, FTYPE the1);

void sigact(void)
{
  int i,j,k;
  int intsimp;
  FTYPE Sig,Sigar[N2][N1],ijudge,cs2a[N2][N1];
  FTYPE cs2temp,qpara2;
  static FTYPE tgriddiag;
  static int igriddiag;
  static int first=0;
  intsimp=1;
  twolayer=1;
  if (first==0){
    igriddiag=(int)((t-tstart)/DTi)+ireenter;
    tgriddiag=tstart+((FTYPE)(igriddiag)*DTi)-1.E-12;
    first=1;
  } 
  LOOPF{
    alpha_reala[0][j][i]=alpha_real0;
//    fprintf(grid_file,"x:%10.4g y:%10.4g rho:%12.6e siga:%12.6e ",x[2][1][i]*fabs(cos(x[2][2][j])),x[2][1][i]*sin(x[2][2][j]),s[1][0][j][i],Siga);
  }
//  fprintf(grid_file,"\n");
 for(j=0;j<N2;j++){
	if(x[2][2][j]<=M_PI/2.)N2tem=j;
 }
 if(twolayer==1){
  if(intsimp==1){
   for(i=0;i<N1;i++){
     Sig=0.;	
     for(j=0;j<=N2tem;j++){
         Sig=Sigar[j][i]=Sig+0.5*(s[1][0][j][i]+s[1][0][j-1][i])*x[2][1][i]*(x[2][2][j]-x[2][2][j-1]);
	 cs2temp=cs2a[j][i]=gam*(gam-1.0)*s[2][0][j][i]/s[1][0][j][i];
	 if (Sig<Siga||cs2temp>cstcrit2){ 
           alpha_reala[0][j][i]=alpha_real0;
         }
         else {
           alpha_reala[0][j][i]=0.;
         }
     }
     Sig=0.;
     for(j=N2-1;j>N2tem;j--){
         Sig=Sigar[j][i]=Sig+0.5*(s[1][0][j][i]+s[1][0][j+1][i])*x[2][1][i]*(x[2][2][j+1]-x[2][2][j]); 
         cs2temp=cs2a[j][i]=gam*(gam-1.0)*s[2][0][j][i]/s[1][0][j][i];
         if (Sig<Siga||cs2temp>cstcrit2){
           alpha_reala[0][j][i]=alpha_real0;
         }
         else {
           alpha_reala[0][j][i]=0.;
         }
     }
     qpara2=gam*(gam-1.)*s[2][0][N2tem][i]/s[1][0][N2tem][i]*(1./x[2][1][i]/x[2][1][i]/x[2][1][i])/M_PI/M_PI/2./Sigar[N2tem][i]/2./Sigar[N2tem][i]/1.*MASSBH*MASSBH;
     for(j=0;j<=N2-1;j++){
	if(alpha_reala[0][j][i]==0.){
	  alpha_reala[0][j][i]=exp(-qpara2);
	}
     }
   }
  } else if (intsimp==0){
   for(i=0;i<N1-1;i++) {
     ijudge=0.;
     for(j=N2/2;j<N2;j++){
       cs2temp=cs2a[j][i]=gam*(gam-1.0)*s[2][0][j][i]/s[1][0][j][i];
       if (ijudge<=N2/20.){
         Sig=Sigar[j][i]=siginte(i,j);
       }
       else {
         Sig=Sigar[j][i]=0.;
       }
       if (Sig<Siga||cs2temp>cstcrit2){
	ijudge=ijudge+1.0;
    	alpha_reala[0][j][i]=alpha_real0;
       }	
       else {
 	alpha_reala[0][j][i]=alpha_real0/alphafrac;
      }
                   //	fprintf(grid_file,"x:%10.4g y:%10.4g rho:%10.4g sig:%12.6g alpha:%12.6g ",x[2][1][i]*fabs(cos(x[2][2][j])),x[2][1][i]*sin(x[2][2][j]),s[1][0][j][i],Sig,alpha_reala[0][j][i]);
                   //       printf("%d y %f sig %f \n",iab[0][j][i],sqrt(pow(pow(10.,(N1-1)*2./N1),2)-x[2][1][i]*sin(x[2][2][j])*x[2][1][i]*sin(x[2][2][j]))-x[2][1][i]*fabs(cos(x[2][2][j])),Sig);
     }
                   //    fprintf(grid_file,"\n i: %d\n",i);
     ijudge=0.;
     for(j=N2/2-1;j>=0;j--){
       cs2temp=cs2a[j][i]=gam*(gam-1.0)*s[2][0][j][i]/s[1][0][j][i];
       if (ijudge<=N2/20.){
         Sig=Sigar[j][i]=siginte(i,j);
       }
       else {
         Sig=Sigar[j][i]=0.;
       }
       if (Sig<Siga||cs2temp>cstcrit2){
	 ijudge=ijudge+1.0;
         alpha_reala[0][j][i]=alpha_real0;
       } 
       else {
         alpha_reala[0][j][i]=alpha_real0/alphafrac;
      }
    }
   } 
  }
  bound(alpha_reala,NULL,-10-1,0,0);	
 }
  if(t>=tgriddiag&&igriddiag<=10){
    for(i=0;i<N1-1;i++){
      for(j=0;j<N2;j++){
       fprintf(grid_file,"x:%10.4g y:%10.4g rho:%15.6e sig:%15.6e t:%15.6e alpha:%12.6g ",x[2][1][i]*fabs(cos(x[2][2][j])),x[2][1][i]*sin(x[2][2][j]),s[1][0][j][i],Sigar[j][i],tirr0[0][j][i],alpha_reala[0][j][i]);
      }   
    fprintf(grid_file,"\n i: %d t %f igrid %d Siga: %f Cs2crit %f\n",i,t,igriddiag,Siga,cstcrit2);
    }
    igriddiag++;
    tgriddiag=tstart+(FTYPE)(igriddiag)*DTi;
  }
}





FTYPE siginte(int iinit, int jinit)
{
  FTYPE Sig=0.;
  int i,j;
  FTYPE xini,len=0.;
  FTYPE r1=0.,the1=0.;
  xini=x[2][1][iinit]*sin(x[2][2][jinit]);
  r1=x[2][1][iinit];
  the1=x[2][2][jinit];
  if (x[2][2][jinit]< M_PI/2.) {
    len=gridl(xini,iinit,jinit,r1,the1);
    Sig=Sig+s[1][0][jinit][iinit]*len ; 
    i=iinit;
    j=jinit;
    r1=r2;
    the1=the2;
    while (! (i>=(N1-2)&&pflag==2))  {
      if(pflag==1){
	j--;
      }
      if(pflag==2){
        i++;
      }
      len=gridl(xini,i,j,r1,the1);
      Sig=Sig+s[1][0][j][i]*len ;
      r1=r2;
      the1=the2;
    }    
  } else {
    len=gridld(xini,iinit,jinit,r1,the1);
    Sig=Sig+s[1][0][jinit][iinit]*len ;
    i=iinit;
    j=jinit;
    r1=r2;
    the1=the2;
    while(! (i>=(N1-2)&&pflag==2) )  {
      if(pflag==1){
        j++;
      }
      if(pflag==2){
        i++;
      }
      len=gridld(xini,i,j,r1,the1);
      Sig=Sig+s[1][0][j][i]*len ;
      r1=r2;
      the1=the2;
    }
  }
  
  return Sig;
}

FTYPE gridl(FTYPE xini,int i, int j,FTYPE r1, FTYPE the1)
{
  FTYPE len;
  if(xini>x[1][1][i]*sin(x[1][2][j])&&xini<x[1][1][i+1]*sin(x[1][2][j])){
    pflag=1;
    r2=xini/sin(x[1][2][j]);
    the2=x[1][2][j];
  }else if (xini>x[1][1][i+1]*sin(x[1][2][j])&&xini<x[1][1][i+1]*sin(x[1][2][j+1])){
    pflag=2;
    r2=x[1][1][i+1];
    the2=asin(xini/r2);
  }else {
    fprintf(fail_file,"wrong grid calculation in active layer at xini: %f, i: %d, j: %d",xini,i,j);
    exit(0);
  }
    len=(r2-r1)*cos(0.5*(the1+the2))-0.5*(r1+r2)*sin(0.5*(the1+the2))*(the2-the1);
    return len;
}

FTYPE gridld(FTYPE xini,int i, int j,FTYPE r1, FTYPE the1)
{
  FTYPE len;
  if(xini>x[1][1][i]*sin(x[1][2][j+1])&&xini<x[1][1][i+1]*sin(x[1][2][j+1])){
    pflag=1;
    r2=xini/sin(x[1][2][j+1]);
    the2=x[1][2][j+1];
  }else if (xini>x[1][1][i+1]*sin(x[1][2][j+1])&&xini<x[1][1][i+1]*sin(x[1][2][j])){
    pflag=2;
    r2=x[1][1][i+1];
    the2=M_PI-asin(xini/r2);
  }else {
    fprintf(fail_file,"wrong grid calculation in active layer at xini: %f, i: %d, j: %d",xini,i,j);
    exit(0);
  }
    len=-((r2-r1)*cos(0.5*(the1+the2))-0.5*(r1+r2)*sin(0.5*(the1+the2))*(the2-the1));
    return len; 
}

