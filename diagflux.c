#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif



void hydro_flux_rect(int dir,int which,FTYPE (*fl)[N3M][N2M][N1M]) // hydro_flux() for LOOPTYPE==1
{
  // 1: rho
  // 2: en
  // 3: bz
  // 4: ke
  // 5: vx
  // 6: vy
  // 7: vz

  int i,j,k;
  FTYPE ftemp,ftempv;
  FTYPE ftemp2;
  FTYPE sign;
  int io;
  int doit;

  // rho-x1
  if((dir==1)&&(which==1)){
    // capture radial flux of mass
    for(io=0;io<2;io++){

      if((io==0)&&(reflectix1==0)){
	i=intix1;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox1==0)){
	if(N1>1) i=intox1; else i=intix1;
	sign=1;
	doit=1;
      }
      else doit=0;

      if(doit){
	for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){
	  ftemp=fl[1][k][j][i]*dx[1][3][k]; // corrected since didn't do dx3 in sweep flux since it cancelled with volume division there
	  // mass
	  losss[1][dir][io][k][j]+=sign*ftemp;
	  // grav pot energy
	  losss[3][dir][io][k][j]+=sign*ftemp*z2e_1(s[3],k,j,i);
	
	  // enthalpy (computed above instead)
	  //losss[2][dir][io][k][j]+=-gam*z2e_1(s[2],k,j,i)/z2e_1(s[1],k,j,i)*ftemp;
	  
	  // adding up loss like this intead of as commented out in velocities gives better results in bondi runs
	  
	  if(kever==1){
	    // vx1-inner-ke
	    ftempv=v[1][1][k][j][i];
	    lossv[1][0][dir][io][k][j]+=sign*0.5*ftemp*ftempv*ftempv;
	    // probably should do j=0 if skipix2=0
	    ftempv=v2tov1(v[1][2],k,j,i);
	    lossv[1][0][dir][io][k][j]+=sign*0.5*ftemp*ftempv*ftempv; // vx2-inner-ke
	    ftempv=v3tov1(v[1][3],k,j,i);
	    lossv[1][0][dir][io][k][j]+=sign*0.5*ftemp*ftempv*ftempv; // vx3-inner-ke
	  }
	}
      }
    }
  }


  // en-x1
  if((dir==1)&&(which==2)){

    for(io=0;io<2;io++){
    
      if((io==0)&&(reflectix1==0)){
	i=intix1;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox1==0)){
	if(N1>1) i=intox1; else i=intix1;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){
	  
	  // capture flux of enthalpy(only true for ideal gas EOS)
	  losss[2][dir][io][k][j]+=sign*gam*fl[1][k][j][i]*dx[1][3][k];
	}
      }
    }
  }

  if((COMPDIM==2)&&(mag)&&(HYDROBZFLUX==1)){
    // bz-x1
    if((dir==1)&&(which==3)){
      for(io=0;io<2;io++){
	
	if((io==0)&&(reflectix1==0)){
	  i=intix1;
	  sign=-1;
	  doit=1;
	}
	else if((io==1)&&(reflectox1==0)){
	  if(N1>1) i=intox1; else i=intix1;
	  sign=1;
	  doit=1;
	}
	else doit=0;
	
	if(doit){
	  for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){
	    
	    // capture magflux3 -- actually bz/(rsin(theta))
	    lossv[2][3][dir][io][k][j]+=sign*fl[1][k][j][i]*dx[1][2][j]/(DS11(k,j,i)*dx[1][1][i]);
	  }
	}
      }
    }
  }
  // ke is taken care of in other terms to get better calculation

  // note that for both x1 and x2 that there are volume terms in spc for x1/x2 momentum fluxes that aren't actually computed here!  Thus, only really valid for cartesian coords.

  // due to fact that fl[N1] doesn't exist, we use fl[N1-1] when periodicx1==0 and otherwise fl[N1]=fl[0]!
  // when using MPI, this works too, and don't care about flux on "interior" surfaces so no special conditions

  // vx-x1
  if((dir==1)&&(which==5)){

    for(io=0;io<2;io++){
      
      if((io==0)&&(reflectix1==0)){
	if(periodicx1) i=-N1OFF; else i=intix1;

	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox1==0)){
	if(periodicx1) i=N1-1; else i=intox1; // intox1==N1 should be true when periodicx1
	if(N1==1) i=intix1; // overrules

	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){
	  
	  ftemp=fl[1][k][j][i]*dx[1][3][k]; // rho*vr*vr term // don't move to boundary since fl[1] may not be defined other places
	  ftemp2=((gam-1.0)*s[2][k][j][i]+s[1][k][j][i]*s[3][k][j][i])*G2(2,i)*G3(2,i)*G4(2,j)*dx[1][2][j]*dx[1][3][k]; // p+rho*phi flux
	  lossv[1][1][dir][io][k][j]+=sign*(ftemp + ftemp2); // total flux
	  if(kever==0){
	    lossv[1][0][dir][io][k][j]+=sign*ftemp*e2z_1(v[1][1],k,j,i)*0.5;
	  }
	}
      }
    }
  }

  // vy-x1
  if((dir==1)&&(which==6)){

    for(io=0;io<2;io++){
      
      if((io==0)&&(reflectix1==0)){
	i=intix1;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox1==0)){
	i=intox1;
	if(N1==1) i=intix1;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){

	  ftemp=fl[1][k][j][i]*dx[1][3][k];
	  lossv[1][2][dir][io][k][j]+=sign*ftemp;
	  // add in kinetic energy loss
	  if(kever==0){
	    lossv[1][0][dir][io][k][j]+=sign*ftemp*(z2e_1(v[1][2],k,j,i)*0.5/G2(1,i));
	  }
	}
      }
    }
  }


  // vz-x1
  if((dir==1)&&(which==7)){

    for(io=0;io<2;io++){
      
      if((io==0)&&(reflectix1==0)){
	i=intix1;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox1==0)){
	i=intox1;
	if(N1==1) i=intix1;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){
	  // capture ang mom(s3) boundary loss

	  ftemp=fl[1][k][j][i]*dx[2][3][k];
	  lossv[1][3][dir][io][k][j]+=sign*ftemp;
	  // add in kinetic energy
	  if(kever==0){
	    lossv[1][0][dir][io][k][j]+=sign*ftemp*(z2e_1(v[1][3],k,j,i)*0.5/(G3(1,i)*G4(2,j)));
	  }
	}
      }
    }
  }

  ////////////////////////////
  ///////////////////////////
  // x2
  //////////////////////////
  //////////////////////////


  // rho-x2
  if((dir==2)&&(which==1)){

    for(io=0;io<2;io++){
      
      if((io==0)&&(reflectix2==0)){
	j=intix2;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox2==0)){
	if(N2>1) j=intox2; else j=intix2;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){

	for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
	  ftemp=fl[2][k][j][i]*dx[1][3][k];
	  //mass
	  losss[1][dir][io][k][i]+=sign*ftemp;
	  //losss[2][dir][io][k][i]+=0; // don't use
	  
	  // grav pot energy
	  losss[3][dir][io][k][i]+=sign*ftemp*z2e_2(s[3],k,j,i);

	  if(kever==1){	
	    // see sweepx for comment on why this done here instead of above commented ke's
	    ftempv=v[1][2][k][j][i];
	    lossv[1][0][dir][io][k][i]+=sign*0.5*ftemp*ftempv*ftempv; // vy
	    ftempv=v1tov2(v[1][1],k,j,i);
	    lossv[1][0][dir][io][k][i]+=sign*0.5*ftemp*ftempv*ftempv; // vx
	    ftempv=v3tov2(v[1][3],k,j,i);
	    lossv[1][0][dir][io][k][i]+=sign*0.5*ftemp*ftempv*ftempv; // vz
	  }
	}
      }
    }
  }

  // en-x2
  if((dir==2)&&(which==2)){

    for(io=0;io<2;io++){
      
      if((io==0)&&(reflectix2==0)){
	j=intix2;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox2==0)){
	if(N2>1) j=intox2; else j=intix2;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){

	for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
	  // capture flux of enthalpy
	  ftemp=fl[2][k][j][i]*dx[1][3][k];
	  losss[2][dir][io][k][i]+=sign*gam*ftemp;
	}
      }
    }
  }
  
  if((COMPDIM==2)&&(mag)&&(HYDROBZFLUX==1)){
    // bz-x2
    if((dir==2)&&(which==3)){
      
      for(io=0;io<2;io++){
	
	if((io==0)&&(reflectix2==0)){
	  j=intix2;
	  sign=-1;
	  doit=1;
	}
	else if((io==1)&&(reflectox2==0)){
	  if(N2>1) j=intox2; else j=intix2;
	  sign=1;
	  doit=1;
	}
	else doit=0;
	
	if(doit){
	  
	  for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
	    // capture magnetic flux 3 component (bz/(r*sin(theta)) actually)
	    
	    lossv[2][3][dir][io][k][i]+=sign*fl[2][k][j][i]*dx[1][1][i]/(DS12(k,j,i)*x[2][1][i]*dx[1][2][j]);
	  }
	}
      }
    }
  }

  // ke taken care of in other terms

  // vx-x2
  if((dir==2)&&(which==5)){

    for(io=0;io<2;io++){
      
      if((io==0)&&(reflectix2==0)){
	j=intix2;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox2==0)){
	if(N2>1) j=intox2; else j=intix2;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	
	for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){

	  ftemp=fl[2][k][j][i]*dx[1][3][k];
	  lossv[1][1][dir][io][k][i]+=sign*ftemp;
	  // ke loss
	  if(kever==0){
	    lossv[1][0][dir][io][k][i]+=sign*ftemp*(z2e_2(v[1][1],k,j,i)*0.5);
	  }
	}
      }
    }
  }



  // vy-x2
  if((dir==2)&&(which==6)){

    for(io=0;io<2;io++){
      
      if((io==0)&&(reflectix2==0)){
	if(periodicx2) j=-N2OFF; else j=intix2;

	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox2==0)){
	if(periodicx2) j=N2-1; else j=intox2; // intox2==N2 should be true when periodicx1
	if(N2==1) j=intix2; // overrules

	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	
	for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){

	  ftemp=fl[2][k][j][i]*dx[1][3][k]; // not on boundary!
	  ftemp2=((gam-1.0)*s[2][k][j][i]+s[1][k][j][i]*s[3][k][j][i])*G2(2,i)*G3(2,i)*G4(2,j)*dx[1][3][k]*dx[1][1][i]; // p+rho*phi flux
	  lossv[1][2][dir][io][k][i]+=sign*(ftemp+ftemp2);
	  if(kever==0){
	    lossv[1][0][dir][io][k][i]+=sign*ftemp*e2z_2(v[1][2],k,j,i)*0.5/G2(2,i);
	  }
	}
      }
    }
  }


  // vz-x2
  if((dir==2)&&(which==7)){

    for(io=0;io<2;io++){
      
      if((io==0)&&(reflectix2==0)){
	j=intix2;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox2==0)){
	j=intox2;
	if(N2==1) j=intix2;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	
	for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
	  
	  ftemp=fl[2][k][j][i]*dx[2][3][k];
	  lossv[1][3][dir][io][k][i]+=sign*ftemp;
	  // ke loss
	  if(kever==0){
	    lossv[1][0][dir][io][k][i]+=sign*ftemp*z2e_2(v[1][3],k,j,i)*0.5*OG4(1,j)/G3(2,i);
	  }
	}
      }
    }
  }




  ////////////////////////////
  ///////////////////////////
  // x3
  //////////////////////////
  //////////////////////////


  // rho-x3
  if((dir==3)&&(which==1)){

    for(io=0;io<2;io++){
      
      if((io==0)&&(reflectix3==0)){
	k=intix3;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox3==0)){

	k=intox3;
	if(N3==1) k=intix3;

	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){

	for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
	  ftemp=fl[3][k][j][i];
	  //mass
	  losss[1][dir][io][j][i]+=sign*ftemp;
	  //losss[2][dir][io][j][i]+=0; // don't use
	  
	  // grav pot energy
	  losss[3][dir][io][j][i]+=sign*ftemp*z2e_3(s[3],k,j,i);

	  if(kever==1){	
	    // see sweepx for comment on why this done here instead of above commented ke's
	    ftempv=v2tov3(v[1][2],k,j,i);
	    lossv[1][0][dir][io][j][i]+=sign*0.5*ftemp*ftempv*ftempv; // vy
	    ftempv=v1tov3(v[1][1],k,j,i);
	    lossv[1][0][dir][io][j][i]+=sign*0.5*ftemp*ftempv*ftempv; // vx
	    ftempv=v[1][3][k][j][i];
	    lossv[1][0][dir][io][j][i]+=sign*0.5*ftemp*ftempv*ftempv; // vz
	  }
	}
      }
    }
  }

  // en-x3
  if((dir==3)&&(which==2)){

    for(io=0;io<2;io++){
      
      if((io==0)&&(reflectix3==0)){
	k=intix3;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox3==0)){

	k=intox3;
	if(N3==1) k=intix3;

	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){

	for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
	  // capture flux of enthalpy
	  
	  losss[2][dir][io][j][i]+=sign*gam*fl[3][k][j][i];
	}
      }
    }
  }
  
  // no bz-x3

  // ke taken care of in other terms

  // vx-x3
  if((dir==3)&&(which==5)){


    for(io=0;io<2;io++){
      
      if((io==0)&&(reflectix3==0)){
	k=intix3;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox3==0)){
	k=intox3;
	if(N3==1) k=intix3;

	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	
	for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
	  
	  lossv[1][1][dir][io][j][i]+=sign*fl[3][k][j][i];
	  // ke loss
	  if(kever==0){
	    lossv[1][0][dir][io][j][i]+=sign*fl[3][k][j][i]*(z2e_3(v[1][1],k,j,i)*0.5);
	  }
	}
      }
    }
  }



  // vy-x3
  if((dir==3)&&(which==6)){


   for(io=0;io<2;io++){
      
      if((io==0)&&(reflectix3==0)){
	k=intix3;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox3==0)){
	k=intox3;
	if(N3==1) k=intix3;

	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	
	for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){

	  ftemp=fl[3][k][j][i]; // not on boundary!
	  lossv[1][2][dir][io][j][i]+=sign*ftemp;
	  if(kever==0){
	    lossv[1][0][dir][io][j][i]+=sign*ftemp*(z2e_3(v[1][2],k,j,i)*0.5/G2(2,i));
	  }
	}
      }
    }
  }


  // vz-x3
  if((dir==3)&&(which==7)){

   for(io=0;io<2;io++){
      
      if((io==0)&&(reflectix3==0)){
	if(periodicx3) k=-N3OFF; else k=intix3;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox3==0)){
	if(periodicx3) k=N3-1; else k=intox3; // intox3==N3 should be true when periodicx3
	if(N3==1) k=intix3; // overrules

	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	
	for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
	  
	  ftemp=fl[3][k][j][i]; // not really on boundary
	  ftemp2=((gam-1.0)*s[2][k][j][i]+s[1][k][j][i]*s[3][k][j][i])*G2(2,i)*G3(2,i)*G4(2,j)*dx[1][2][j]*dx[1][1][i]; // p+rho*phi flux
	  lossv[1][3][dir][io][j][i]+=sign*(ftemp+ftemp2);
	  // ke loss
	  if(kever==0){
	    lossv[1][0][dir][io][j][i]+=sign*ftemp*e2z_3(v[1][3],k,j,i)*0.5/(G3(2,i)*G4(2,j));
	  }
	}
      }
   }
  }

} 



// reconsider reflecting conditions here

// in 3D the locations are:
// (x,y,z)
// sig11,22,33: (.5,.5,.5)
// sig13: (0,.5,0)
// sig12: (0,0,.5)
// sig23: (.5,0,0)
void viscous_flux_rect(void){ // viscous_flux() for LOOPTYPE==1
  // right before going to advection, compute fluxes of energy and momentum due to viscosity
  
  int i,j,k,l,m ;
  int dir,io;
  FTYPE sign;
  int doit;
  FTYPE ftemp,ftemp1,ftemp2;
  FTYPE subftemp;
  FTYPE Length;
  FTYPE flux;
  FTYPE crapf;
  
  // the viscous calculation is flux conservative(for momentum).
  if(vischeat){
    // capture flux of internal energy
    
    if(VISCE31){
      dir=1;
      for(io=0;io<2;io++){
	if(N1==1) break;
	if(io==0){
	  i=intix1;
	  sign=-1;
	  doit=1;
	}
	else if(io==1){
	  if(N1>1) i=intox1; else i=intix1;
	  sign=1;
	  doit=1;
	}
	else doit=0;
	
	if(doit){
	  for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){
	    flux=dt*sigma[3][1][k][j][i]*z2e_1(v[1][3],k,j,i)*DS11(k,j,i)*dx[2][3][k];
	    lossvisc[1][dir][io][k][j]+=sign*flux;
	  }
	}
      }
    }
    if(VISCE13){
      dir=3;
      for(io=0;io<2;io++){
	if(N3==1) break;
	
	if(io==0){
	  k=intix3;
	  sign=-1;
	  doit=1;
	}
	else if(io==1){
	  if(N3>1) k=intox3; else k=intix3;
	  sign=1;
	  doit=1;
	}
	else doit=0;
	
	if(doit){
	  for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
	    flux=dt*sigma[1][3][k][j][i]*z2e_3(v[1][1],k,j,i)*DS23(k,j,i);
	    lossvisc[1][dir][io][j][i]+=sign*flux;
	  }
	}
      }
    }
    
    if(VISCE23){
      dir=2;
      for(io=0;io<2;io++){
	if(N2==1) break;
	
	if(io==0){
	  j=intix2;
	  sign=-1;
	  doit=1;
	}
	else if(io==1){
	  if(N2>1) j=intox2; else j=intix2;
	  sign=1;
	  doit=1;
	}
	else doit=0;
	
	if(doit){
	  
	  for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
	    flux=dt*sigma[2][3][k][j][i]*z2e_2(v[1][3],k,j,i)*DS12(k,j,i)*dx[2][3][k];
	    lossvisc[1][dir][io][k][i]+=sign*flux;
	  }
	}
      }
    }

    if(VISCE32){
      dir=3;
      for(io=0;io<2;io++){
	if(N3==1) break;
	
	if(io==0){
	  k=intix3;
	  sign=-1;
	  doit=1;
	}
	else if(io==1){
	  if(N3>1) k=intox3; else k=intix3;
	  sign=1;
	  doit=1;
	}
	else doit=0;
	
	if(doit){
	  
	  for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
	    flux=dt*sigma[3][2][k][j][i]*z2e_3(v[1][2],k,j,i)*DS33(k,j,i);
	    lossvisc[1][dir][io][j][i]+=sign*flux;
	  }
	}
      }
    }

    
    if(VISCE11){
      dir=1;
      for(io=0;io<2;io++){
	if(N1==1) break;
	
	if(io==0){
	  i=intix1;
	  sign=-1;
	  doit=1;
	}
	else if(io==1){
	  if(N1>1) i=intox1; else i=intix1;
	  sign=1;
	  doit=1;
	}
	else doit=0;
	
	if(doit){
	  for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){
	    flux=dt*z2e_1(sigma[1][1],k,j,i)*v[1][1][k][j][i]*DS11(k,j,i)*dx[1][3][k];
	    lossvisc[1][dir][io][k][j]+=sign*flux;
	  }
	}
      }
    }
    
    if(VISCE22){
      dir=2;
      for(io=0;io<2;io++){
	if(N2==1) break;
	
	if(io==0){
	  j=intix2;
	  sign=-1;
	  doit=1;
	}
	else if(io==1){
	  if(N2>1) j=intox2; else j=intix2;
	  sign=1;
	  doit=1;
	}
	else doit=0;
	
	if(doit){
	  
	  for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
	    flux=dt*z2e_2(sigma[2][2],k,j,i)*v[1][2][k][j][i]*DS12(k,j,i)*dx[1][3][k];
	    lossvisc[1][dir][io][k][i]+=sign*flux;
	  }
	}
      }
    }

    if(VISCE33){
      dir=3;
      for(io=0;io<2;io++){
	if(N3==1) break;
	
	if(io==0){
	  k=intix3;
	  sign=-1;
	  doit=1;
	}
	else if(io==1){
	  if(N3>1) k=intox3; else k=intix3;
	  sign=1;
	  doit=1;
	}
	else doit=0;
	
	if(doit){
	  
	  for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
	    flux=dt*z2e_3(sigma[3][3],k,j,i)*v[1][3][k][j][i]*DS13(k,j,i);
	    lossvisc[1][dir][io][j][i]+=sign*flux;
	  }
	}
      }
    }
    
    if(VISCE12){
      // only flux is radial!

      dir=1;
      for(io=0;io<2;io++){
	if(N1==1) break;
	
	if(io==0){
	  i=intix1;
	  sign=-1;
	  doit=1;
	}
	else if(io==1){
	  if(N1>1) i=intox1; else i=intix1;
	  sign=1;
	  doit=1;
	}
	else doit=0;
	
	if(doit){
	  for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){
	    flux=dt*sigma[1][2][k][j][i]*z2e_1(v[1][2],k,j,i)*DS31(k,j,i)*dx[1][3][k];
	    lossvisc[1][dir][io][k][j]+=sign*flux;
	  }
	}
      }
    }
    
    if(VISCE21){
      dir=2;
      for(io=0;io<2;io++){
	if(N2==1) break;
	
	if(io==0){
	  j=intix2;
	  sign=-1;
	  doit=1;
	}
	else if(io==1){
	  if(N2>1) j=intox2; else j=intix2;
	  sign=1;
	  doit=1;
	}
	else doit=0;
	
	if(doit){
	  
	  for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
	    flux=dt*sigma[2][1][k][j][i]*z2e_2(v[1][1],k,j,i)*DS22(k,j,i)*dx[1][3][k];
	    lossvisc[1][dir][io][k][i]+=sign*flux;
	  }
	}
      }
    }
  }// end capture energy flux  


    
  ///////////
  // momentum fluxes due to viscosity
  // in vector positions
  // just copy of magnetic version and stick in sigma in right vector location
  // since interp on sigma sometimes, assume sigma exists at least on LOOPH

  //vx-x1 flux  without volume terms
  if(VISCE11){
    dir=1;
    for(io=0;io<2;io++){
      if(N1==1) break;
      
      if((io==0)&&(reflectix1==0)){
	i=intix1;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox1==0)){
	if(N1>1) i=intox1; else i=intix1;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){        
	  lossv[1][1][dir][io][k][j]+=sign*(z2e_1(sigma[1][1],k,j,i))*DS11(k,j,i)*dx[1][3][k]*dt;
	}
      }
    }
  }

  // vx-x2 flux   without volume terms
  if(VISCE21){
    dir=2;
    for(io=0;io<2;io++){
      if(N2==1) break;
      
      if((io==0)&&(reflectix2==0)){
	j=intix2;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox2==0)){
	if(N2>1) j=intox2; else j=intix2;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	
	for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
	  
	  lossv[1][1][dir][io][k][i]+=sign*(e2z_1(sigma[2][1],k,j,i))*DS12(k,j,i)*dx[1][3][k]*dt;
	}
      }
    }
  }

  // vx-x3 flux   without volume terms
  if(VISCE31){
    dir=3;
    for(io=0;io<2;io++){
      if(N3==1) break;
      
      if((io==0)&&(reflectix3==0)){
	k=intix3;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox3==0)){
	if(N3>1) k=intox3; else k=intix3;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	
	for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
	  
	  lossv[1][1][dir][io][j][i]+=sign*(e2z_1(sigma[3][1],k,j,i))*DS13(k,j,i)*dt;
	}
      }
    }
  }


  //vy-x1 flux  without volume terms
  if(VISCE12){
    dir=1;
    for(io=0;io<2;io++){
      if(N1==1) break;
      
      if((io==0)&&(reflectix1==0)){
	i=intix1;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox1==0)){
	if(N1>1) i=intox1; else i=intix1;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){
	  
	  lossv[1][2][dir][io][k][j]+=sign*(e2z_2(sigma[1][2],k,j,i))*G2(1,i)*DS11(k,j,i)*dx[1][3][k]*dt; // g12i is from absorption of one volume term
	}
      }
    }
  }

  // vy-x2 flux   without volume terms
  if(VISCE22){
    dir=2;
    for(io=0;io<2;io++){
      if(N2==1) break;
      
      if((io==0)&&(reflectix2==0)){
	j=intix2;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox2==0)){
	if(N2>1) j=intox2; else j=intix2;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	
	for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
	  
	  lossv[1][2][dir][io][k][i]+=sign*(z2e_2(sigma[2][2],k,j,i))*G2(2,i)*DS12(k,j,i)*dx[1][3][k]*dt; // G2(2,i) from definition of s2=h2*rho*v2
	}
      }
    }
  }

  //vy-x3 flux  without volume terms
  if(VISCE32){
    dir=3;
    for(io=0;io<2;io++){
      if(N3==1) break;
      
      if((io==0)&&(reflectix3==0)){
	k=intix3;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox3==0)){
	if(N3>1) k=intox3; else k=intix3;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
	  
	  lossv[1][2][dir][io][j][i]+=sign*(e2z_2(sigma[3][2],k,j,i))*G2(2,i)*DS13(k,j,i)*dt; // G2(2,i) from definition of s2=h2*rho*v2
	}
      }
    }
  }


  //vz-x1 flux
  if(VISCE13){
    dir=1;
    for(io=0;io<2;io++){
      if(N1==1) break;
      
      if((io==0)&&(reflectix1==0)){
	i=intix1;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox1==0)){
	if(N1>1) i=intox1; else i=intix1;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){
	  
	  lossv[1][3][dir][io][k][j]+=sign*(e2z_3(sigma[1][3],k,j,i))*G3(1,i)*G4(2,j)*DS11(k,j,i)*dx[1][3][k]*dt; // g13i*g24j is from absorbed volume term
	}
      }
    }
  }
  
  // vz-x2 flux
  if(VISCE23){
    dir=2;
    for(io=0;io<2;io++){
      if(N2==1) break;
      
      if((io==0)&&(reflectix2==0)){
	j=intix2;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox2==0)){
	if(N2>1) j=intox2; else j=intix2;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	
	for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
	  
	  lossv[1][3][dir][io][k][i]+=sign*(e2z_3(sigma[2][3],k,j,i))*G3(2,i)*G4(1,j)*DS12(k,j,i)*dx[1][3][k]*dt;// g23i*g14j is from absorbed volume term
	}
      }
    }
  }


  // vz-x3 flux   without volume terms
  if(VISCE33){
    dir=3;
    for(io=0;io<2;io++){
      if(N3==1) break;

      if((io==0)&&(reflectix3==0)){
	k=intix3;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox3==0)){
	if(N3>1) k=intox3; else k=intix3;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
	  
	  lossv[1][3][dir][io][j][i]+=sign*(z2e_3(sigma[3][3],k,j,i))*G3(2,i)*G4(2,j)*DS13(k,j,i)*dt; // g23i*g24j is from definition of s3=g3*g4*v3 
	}
      }
    }
  }
  
}


// all magnetic fluxes are on vector positions
// very similar to viscous flux in form for momentum terms

// energy, momentum, and magnetic
void magnetic_flux_rect(void) // magnetic_flux() for LOOPTYPE==1
{
  FTYPE b2,vdotb,bbtensor,emf,ftemp;
  int i,j,k;
  int dir,io;
  FTYPE sign;
  int doit;

  // done all together in readable form since not in flux conservative form, unlike viscous or hydro flux
  
  // Energy flux due to magnetic field
  
  //be-x1
  dir=1;
  for(io=0;io<2;io++){
    if(N1==1) break;
    
    if((io==0)&&(reflectix1==0)){
      i=intix1;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox1==0)){
      if(N1>1) i=intox1; else i=intix1;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){      
	vdotb=v[1][1][k][j][i]*v[2][1][k][j][i]+v2tov1(v[1][2],k,j,i)*v2tov1(v[2][2],k,j,i)+v3tov1(v[1][3],k,j,i)*v3tov1(v[2][3],k,j,i);
	b2=v[2][1][k][j][i]*v[2][1][k][j][i]+v2tov1(v[2][2],k,j,i)*v2tov1(v[2][2],k,j,i)+v3tov1(v[2][3],k,j,i)*v3tov1(v[2][3],k,j,i);
	
	lossv[2][0][dir][io][k][j]+=sign*(b2*v[1][1][k][j][i]-vdotb*v[2][1][k][j][i])*DS11(k,j,i)*dx[1][3][k]*dt;
      }
    }
  }
  
  // be-x2
  dir=2;
  for(io=0;io<2;io++){
    if(N2==1) break;
      
    if((io==0)&&(reflectix2==0)){
      j=intix2;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox2==0)){
      if(N2>1) j=intox2; else j=intix2;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      
      for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
	
	vdotb=v1tov2(v[1][1],k,j,i)*v1tov2(v[2][1],k,j,i)+v[1][2][k][j][i]*v[2][2][k][j][i]+v3tov2(v[1][3],k,j,i)*v3tov2(v[2][3],k,j,i);
	b2=v1tov2(v[2][1],k,j,i)*v1tov2(v[2][1],k,j,i)+v[2][2][k][j][i]*v[2][2][k][j][i]+v3tov2(v[2][3],k,j,i)*v3tov2(v[2][3],k,j,i);
	
	lossv[2][0][dir][io][k][i]+=sign*(b2*v[1][2][k][j][i]-vdotb*v[2][2][k][j][i])*DS12(k,j,i)*dx[1][3][k]*dt;
      }
    }
  }

  // be-x3
  dir=3;
  for(io=0;io<2;io++){
    if(N3==1) break;
      
    if((io==0)&&(reflectix3==0)){
      k=intix3;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox3==0)){
      if(N3>1) k=intox3; else k=intix3;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      
      for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
	
	vdotb=v1tov3(v[1][1],k,j,i)*v1tov3(v[2][1],k,j,i)+v2tov3(v[1][2],k,j,i)*v2tov3(v[2][2],k,j,i)+v[1][3][k][j][i]*v[2][3][k][j][i];
	b2=v1tov3(v[2][1],k,j,i)*v1tov3(v[2][1],k,j,i)+v2tov3(v[2][2],k,j,i)*v2tov3(v[2][2],k,j,i)+v[2][3][k][j][i]*v[2][3][k][j][i];
	
	lossv[2][0][dir][io][j][i]+=sign*(b2*v[1][3][k][j][i]-vdotb*v[2][3][k][j][i])*DS13(k,j,i)*dt;
      }
    }
  }
  
  ///////////
  // momentum fluxes due to magnetic field

  //vx-x1 flux  without volume terms
  dir=1;
  for(io=0;io<2;io++){
    if(N1==1) break;
    
    if((io==0)&&(reflectix1==0)){
      i=intix1;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox1==0)){
      if(N1>1) i=intox1; else i=intix1;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){        
	bbtensor=-v[2][1][k][j][i]*v[2][1][k][j][i];
	b2=v[2][1][k][j][i]*v[2][1][k][j][i]+v2tov1(v[2][2],k,j,i)*v2tov1(v[2][2],k,j,i)+v3tov1(v[2][3],k,j,i)*v3tov1(v[2][3],k,j,i);
	
	lossv[1][1][dir][io][k][j]+=sign*(bbtensor+0.5*b2)*DS11(k,j,i)*dx[1][3][k]*dt;
      }
    }
  }
  
  // vx-x2 flux   without volume terms
  dir=2;
  for(io=0;io<2;io++){
    if(N2==1) break;
      
    if((io==0)&&(reflectix2==0)){
      j=intix2;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox2==0)){
      if(N2>1) j=intox2; else j=intix2;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      
      for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
      
	bbtensor=-v1tov2(v[2][1],k,j,i)*v[2][2][k][j][i];
	//      b2=0;
      
	lossv[1][1][dir][io][k][i]+=sign*(bbtensor)*DS12(k,j,i)*dx[1][3][k]*dt;
      }
    }
  }

  // vx-x3 flux   without volume terms
  dir=3;
  for(io=0;io<2;io++){
    if(N3==1) break;
      
    if((io==0)&&(reflectix3==0)){
      k=intix3;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox3==0)){
      if(N3>1) k=intox3; else k=intix3;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      
      for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
      
	bbtensor=-v1tov3(v[2][1],k,j,i)*v[2][3][k][j][i];
	//      b2=0;
      
	lossv[1][1][dir][io][j][i]+=sign*(bbtensor)*DS13(k,j,i)*dt;
      }
    }
  }


  //vy-x1 flux  without volume terms
  dir=1;
  for(io=0;io<2;io++){
    if(N1==1) break;
    
    if((io==0)&&(reflectix1==0)){
      i=intix1;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox1==0)){
      if(N1>1) i=intox1; else i=intix1;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){
      
	bbtensor=-v[2][1][k][j][i]*v2tov1(v[2][2],k,j,i);
	b2=0;
	
	lossv[1][2][dir][io][k][j]+=sign*(bbtensor)*G2(1,i)*DS11(k,j,i)*dx[1][3][k]*dt; // g12i is from absorption of one volume term
      }
    }
  }
  
  // vy-x2 flux   without volume terms
  dir=2;
  for(io=0;io<2;io++){
    if(N2==1) break;
     
    if((io==0)&&(reflectix2==0)){
      j=intix2;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox2==0)){
      if(N2>1) j=intox2; else j=intix2;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      
      for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
      
	bbtensor=-v[2][2][k][j][i]*v[2][2][k][j][i];
	b2=v1tov2(v[2][1],k,j,i)*v1tov2(v[2][1],k,j,i)+v[2][2][k][j][i]*v[2][2][k][j][i]+v3tov2(v[2][3],k,j,i)*v3tov2(v[2][3],k,j,i);
	
	lossv[1][2][dir][io][k][i]+=sign*(bbtensor+0.5*b2)*G2(2,i)*DS12(k,j,i)*dx[1][3][k]*dt; // G2(2,i) from definition of s2=h2*rho*v2
      }
    }
  }

  //vy-x3 flux  without volume terms
  dir=3;
  for(io=0;io<2;io++){
    if(N3==1) break;
    
    if((io==0)&&(reflectix3==0)){
      k=intix3;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox3==0)){
      if(N3>1) k=intox3; else k=intix3;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
      
	bbtensor=-v[2][3][k][j][i]*v2tov3(v[2][2],k,j,i);
	//b2=0;
	
	lossv[1][2][dir][io][j][i]+=sign*(bbtensor)*G2(2,i)*DS13(k,j,i)*dt; // G2(2,i) from definition of s2=h2*rho*v2
      }
    }
  }


  //vz-x1 flux
  dir=1;
  for(io=0;io<2;io++){
    if(N1==1) break;

    if((io==0)&&(reflectix1==0)){
      i=intix1;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox1==0)){
      if(N1>1) i=intox1; else i=intix1;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){
      
	bbtensor=-v[2][1][k][j][i]*v3tov1(v[2][3],k,j,i);
	//b2=0;
	
	lossv[1][3][dir][io][k][j]+=sign*(bbtensor)*G3(1,i)*G4(2,j)*DS11(k,j,i)*dx[1][3][k]*dt; // g13i*g24j is from absorbed volume term
      }
    }
  }
  
  // vz-x2 flux
  dir=2;
  for(io=0;io<2;io++){
    if(N2==1) break;
      
    if((io==0)&&(reflectix2==0)){
      j=intix2;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox2==0)){
      if(N2>1) j=intox2; else j=intix2;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      
      for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
      
	bbtensor=-v[2][2][k][j][i]*v3tov2(v[2][3],k,j,i);
	//b2=0;
	
	lossv[1][3][dir][io][k][i]+=sign*(bbtensor)*G3(2,i)*G4(1,j)*DS12(k,j,i)*dx[1][3][k]*dt;// g23i*g14j is from absorbed volume term
      }
    }
  }


  // vz-x3 flux   without volume terms
  dir=3;
  for(io=0;io<2;io++){
    if(N3==1) break;
      
    if((io==0)&&(reflectix3==0)){
      k=intix3;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox3==0)){
      if(N3>1) k=intox3; else k=intix3;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
      
	bbtensor=-v[2][3][k][j][i]*v[2][3][k][j][i];
	b2=v1tov3(v[2][1],k,j,i)*v1tov3(v[2][1],k,j,i)+v2tov3(v[2][2],k,j,i)*v2tov3(v[2][2],k,j,i)+v[2][3][k][j][i]*v[2][3][k][j][i];
	
	lossv[1][3][dir][io][j][i]+=sign*(bbtensor+0.5*b2)*G3(2,i)*G4(2,j)*DS13(k,j,i)*dt; // g23i*g24j is from definition of s3=g3*g4*v3 
      }
    }
  }
  

  // flux of magnetic flux through boundaries
  
  // bx-x1 flux is 0.

  // bx-x2 flux
  dir=2;
  for(io=0;io<2;io++){
    if(N2==1) break;
      
    if((io==0)&&(reflectix2==0)){
      j=intix2;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox2==0)){
      if(N2>1) j=intox2; else j=intix2;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      
      for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
	// emf is really -emf, so that sign is same as usual
	emf=-(v1tov2(v[1][1],k,j,i)*v[2][2][k][j][i]-v[1][2][k][j][i]*v1tov2(v[2][1],k,j,i)); // -emf_z
	lossv[2][1][dir][io][k][i]+=sign*(emf)*DS12(k,j,i)*dx[1][3][k]*dt;
      }
    }
  }

  // bx-x3 flux
  dir=3;
  for(io=0;io<2;io++){
    if(N3==1) break;
      
    if((io==0)&&(reflectix3==0)){
      k=intix3;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox3==0)){
      if(N3>1) k=intox3; else k=intix3;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      
      for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
	emf=v[1][3][k][j][i]*v1tov3(v[2][1],k,j,i)-v1tov3(v[1][1],k,j,i)*v[2][3][k][j][i]; // emf_y
	lossv[2][1][dir][io][j][i]+=sign*(emf)*DS13(k,j,i)*dt;
      }
    }
  }


  //by/g12 -x1 flux
  dir=1;
  for(io=0;io<2;io++){
    if(N1==1) break;
    
    if((io==0)&&(reflectix1==0)){
      i=intix1;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox1==0)){
      if(N1>1) i=intox1; else i=intix1;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){
      
	emf=(v[1][1][k][j][i]*v2tov1(v[2][2],k,j,i)-v2tov1(v[1][2],k,j,i)*v[2][1][k][j][i]); // emf_z
	lossv[2][2][dir][io][k][j]+=sign*(emf)*DS11(k,j,i)*dx[1][3][k]*dt/G2(1,i); // g12i since different for the field
      }
    }
  }
  
  // by/g22 -x2 flux is 0.

  //by/g22 -x3 flux
  dir=3;
  for(io=0;io<2;io++){
    if(N3==1) break;
    
    if((io==0)&&(reflectix3==0)){
      k=intix3;
      sign=-1;
      doit=1;
    }
    else if((io==1)&&(reflectox3==0)){
      if(N3>1) k=intox3; else k=intix3;
      sign=1;
      doit=1;
    }
    else doit=0;
    
    if(doit){
      for(j=intix2;j<intox2;j++) for(i=intix1;i<intox1;i++){
      
	emf=-(v2tov3(v[1][2],k,j,i)*v[2][3][k][j][i]-v[1][3][k][j][i]*v2tov3(v[2][2],k,j,i)); // -emf_x
	lossv[2][2][dir][io][j][i]+=sign*(emf)*DS13(k,j,i)*dt/G2(2,i); // g12i since different for the field
      }
    }
  }


  // let hydro_flux() compute this term since flux conservative there.
  if((HYDROBZFLUX==0)||(COMPDIM>2)){
  //bz/g13/g42 -x1 flux
    dir=1;
    for(io=0;io<2;io++){
      if(N1==1) break;
      
      if((io==0)&&(reflectix1==0)){
	i=intix1;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox1==0)){
	if(N1>1) i=intox1; else i=intix1;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	for(k=intix3;k<intox3;k++) for(j=intix2;j<intox2;j++){
      
	  emf=-(-v[1][1][k][j][i]*v3tov1(v[2][3],k,j,i)+v3tov1(v[1][3],k,j,i)*v[2][1][k][j][i]); // -emf_y
	  
	  lossv[2][3][dir][io][k][j]+=sign*(emf)*DS11(k,j,i)*dx[1][3][k]*dt/(G3(1,i)*G4(2,j)); // g13i*g24j is since field different
	}
      }
    }
  
  
    // bz/g13/g42 -x2 flux
    dir=2;
    for(io=0;io<2;io++){
      if(N2==1) break;
      
      if((io==0)&&(reflectix2==0)){
	j=intix2;
	sign=-1;
	doit=1;
      }
      else if((io==1)&&(reflectox2==0)){
	if(N2>1) j=intox2; else j=intix2;
	sign=1;
	doit=1;
      }
      else doit=0;
      
      if(doit){
	
	for(k=intix3;k<intox3;k++) for(i=intix1;i<intox1;i++){
      
	  emf=(v[1][2][k][j][i]*v3tov2(v[2][3],k,j,i)-v3tov2(v[1][3],k,j,i)*v[2][2][k][j][i]);
	  
	  lossv[2][3][dir][io][k][i]+=sign*(emf)*DS12(k,j,i)*dx[1][3][k]*dt*OG4(1,j)/G3(2,i);// g23i*g14j is since field different
	}
      }
    }
  }
}



#define USEADVFLUX 1
// 0: compute flux
// 1: use advection flux

// general loop version
// use when COMPUTELOSSDIAG==0 and LOOPTYPE>1
void hydro_flux_gen_adv(int dir,int which,FTYPE (*fl)[N3M][N2M][N1M]) // hydro_flux()
{
  // 1: rho
  // 2: en
  // 3: bz
  // 4: ke
  // 5: vx
  // 6: vy
  // 7: vz

  int i,j,k,m;
  int itemp;
  int looper;
  SFTYPE ftemp1,ftemp2,ftemp3,ftemp4,ftemp5;
  SFTYPE inftemp1,inftemp2;
  int io;
  SFTYPE theflux;
  SFTYPE fluxdir;
  int looperstart,looperend;
  //  return;

  if(dir==1){
    looperstart=23;
    looperend=24;
  }
  else if(dir==2){
    looperstart=25;
    looperend=26;
  }
  else if(dir==3){
    looperstart=27;
    looperend=28;
  }
  
  for(looper=looperstart;looper<=looperend;looper++){
    // enumerate surfaces
    if((looper==23)||(looper==25)||(looper==27)) io=0; else io=1;
    // loop over surfaces zones
    if( ((looper==23)||(looper==24))&&(N1==1) ) break;
    if( ((looper==25)||(looper==26))&&(N2==1) ) break;
    if( ((looper==27)||(looper==28))&&(N3==1) ) break;
    
    LOOPSUPERGEN(looper){
      ftemp1=0;
      ftemp2=0;
      ftemp3=0;
      ftemp4=0;
      ftemp5=0;
      
      fluxdir=bzs[looper][temptempi][0];

      if(which==1){
	if(dir==1){
	  // mass flux
	  if(USEADVFLUX==0){
	    theflux=dt*z2e_1(s[1],k,j,i)*v[1][dir][k][j][i]*G2(1,i)*G3(1,i)*G4(2,j)*dx[1][2][j]*dx[1][3][k];
	  }
	  else{
	    theflux=fl[dir][k][j][i]*dx[1][3][k];
	  }
	  theflux*=fluxdir;
	  
	  ftemp1+=theflux;
	  // corrected since didn't do dx3 in sweep flux since it cancelled with volume division there
	  // grav. pot. energy flux
	  ftemp2+=theflux*z2e_1(s[3],k,j,i);
	  // ke flux
	  if(kever==1){
	    ftemp3+=theflux*0.5*v[1][1][k][j][i]*v[1][1][k][j][i];
	    ftemp4+=theflux*0.5*v2tov1(v[1][2],k,j,i)*v2tov1(v[1][2],k,j,i);
	    ftemp5+=theflux*0.5*v3tov1(v[1][3],k,j,i)*v3tov1(v[1][3],k,j,i);
	  }
	}
	else if(dir==2){
	  // mass flux
	  if(USEADVFLUX==0){
	    theflux=dt*z2e_2(s[1],k,j,i)*v[1][dir][k][j][i]*G3(2,i)*G4(1,j)*dx[1][1][i]*dx[1][3][k];
	  }
	  else{
	    theflux=fl[dir][k][j][i]*dx[1][3][k];
	  }
	  theflux*=fluxdir;
	  
	  ftemp1+=theflux;
	  // grav pot energy
	  ftemp2+=theflux*z2e_2(s[3],k,j,i);
	  
	  if(kever==1){	
	    // see sweepx for comment on why this done here instead of above commented ke's
	    ftemp3+=theflux*0.5*v1tov2(v[1][1],k,j,i)*v1tov2(v[1][1],k,j,i);
	    ftemp4+=theflux*0.5*v[1][2][k][j][i]*v[1][2][k][j][i];
	    ftemp5+=theflux*0.5*v3tov2(v[1][3],k,j,i)*v3tov2(v[1][3],k,j,i);
	  }
	}
	else if(dir==3){
	  // mass flux
	  if(USEADVFLUX==0){
	    theflux=dt*z2e_3(s[1],k,j,i)*v[1][dir][k][j][i]*G2(2,i)*dx[1][2][j]*dx[1][1][i];
	  }
	  else{
	    theflux=fl[dir][k][j][i];
	  }
	  theflux*=fluxdir;
	  
	  ftemp1+=theflux;
	  // grav. pot. energy flux
	  ftemp2+=theflux*z2e_3(s[3],k,j,i);
	  // ke flux
	  if(kever==1){
	    ftemp3+=theflux*0.5*v1tov3(v[1][1],k,j,i)*v1tov3(v[1][1],k,j,i);
	    ftemp4+=theflux*0.5*v2tov3(v[1][2],k,j,i)*v2tov3(v[1][2],k,j,i);
	    ftemp5+=theflux*0.5*v[1][3][k][j][i]*v[1][3][k][j][i];
	  }
	}
      }
      else if(which==2){
	if(dir==1){
	  // capture flux of enthalpy
	  if(USEADVFLUX==0){
	    theflux=dt*gam*z2e_1(s[2],k,j,i)*v[1][dir][k][j][i]*G2(1,i)*G3(1,i)*G4(2,j)*dx[1][2][j]*dx[1][3][k];
	  }
	  else{
	    theflux=gam*fl[dir][k][j][i]*dx[1][3][k];
	  }
	  theflux*=fluxdir;
	  
	  ftemp1+=theflux;
	}
	else if(dir==2){
	  // capture flux of enthalpy
	  if(USEADVFLUX==0){
	    theflux=dt*gam*z2e_2(s[2],k,j,i)*v[1][dir][k][j][i]*G3(2,i)*G4(1,j)*dx[1][1][i]*dx[1][3][k];
	  }
	  else{
	    theflux=gam*fl[dir][k][j][i]*dx[1][3][k];
	  }
	  theflux*=fluxdir;
	  
	  ftemp1+=theflux;
	}
	else if(dir==3){
	  // capture flux of enthalpy
	  if(USEADVFLUX==0){
	    theflux=dt*gam*z2e_3(s[2],k,j,i)*v[1][dir][k][j][i]*G2(2,i)*dx[1][2][j]*dx[1][1][i];
	  }
	  else{
	    theflux=gam*fl[dir][k][j][i];
	  }
	  theflux*=fluxdir;
	  
	  ftemp1+=theflux;
	}
      }
      else if(which==3){
	// capture magflux3 -- actually bz/(rsin(theta))
	if(dir==1){
	  theflux=fl[dir][k][j][i]*dx[1][2][j]/(DS11(k,j,i)*dx[1][1][i]);
	  theflux*=fluxdir;
	  
	  ftemp1+=theflux;
	}
	else if(dir==2){
	  theflux=fl[dir][k][j][i]*dx[1][1][i]/(DS12(k,j,i)*x[2][1][i]*dx[1][2][j]);
	  theflux*=fluxdir;
	  
	  ftemp1+=theflux;
	}
	// no dir==3
      }
      else if(which==5){ // flux of vx momentum
	if(dir==1){
	  // fl[outer N1-like] doesn't exist, so don't ever use ADV flux here
	  if(1||(USEADVFLUX==0)){ // on edge
	    theflux=dt*(z2e_1(s[1],k,j,i)*v[1][1][k][j][i]*v[1][1][k][j][i]+(gam-1.0)*z2e_1(s[2],k,j,i)+z2e_1(s[1],k,j,i)*z2e_1(s[3],k,j,i))*G2(1,i)*dx[1][2][j]*G3(1,i)*G4(2,j)*dx[1][3][k];
	    theflux*=fluxdir;
	    ftemp1+=theflux; // total flux
	    if(kever==0){
	      ftemp2+=theflux*v[1][1][k][j][i]*0.5;
	    }
	  }
	  else{
	    inftemp1=fl[dir][k][j][i]*dx[1][3][k]; // rho*vr*vr term // don't move to boundary since fl[1] may not be defined other places	    
	    inftemp1+=((gam-1.0)*s[2][k][j][i]+s[1][k][j][i]*s[3][k][j][i])*G2(2,i)*G3(2,i)*G4(2,j)*dx[1][2][j]*dx[1][3][k]; // p+rho*phi flux
	    theflux=inftemp1;
	    theflux*=fluxdir;
	    ftemp1+=theflux; // total flux
	    if(kever==0){
	      ftemp2+=theflux*e2z_1(v[1][1],k,j,i)*0.5;
	    }
	  }
	}
	else if(dir==2){
	  if(USEADVFLUX==0){ // on edge
	    theflux=dt*z2e_2(s[1],k,j,i)*v1tov2(v[1][1],k,j,i)*v[1][2][k][j][i]*dx[1][1][i]*G3(2,i)*G4(1,j)*dx[1][3][k];
	    theflux*=fluxdir;
	    ftemp1+=theflux;
	    if(kever==0){
	      // ke loss
	      ftemp2+=theflux*v1tov2(v[1][1],k,j,i)*0.5;
	    }
	  }
	  else{
	    theflux=fl[dir][k][j][i]*dx[1][3][k];// rho*vr*vr term // don't move to boundary since fl[1] may not be defined other places
	    theflux*=fluxdir;	  
	    ftemp1+=theflux;
	    if(kever==0){
	      // ke loss
	      ftemp2+=theflux*z2e_2(v[1][1],k,j,i)*0.5;
	    }
	  }
	}
	else if(dir==3){
	  if(USEADVFLUX==0){ // on edge
	    theflux=dt*z2e_3(s[1],k,j,i)*v1tov3(v[1][1],k,j,i)*v[1][3][k][j][i]*dx[1][1][i]*G2(2,i)*dx[1][2][j];
	    theflux*=fluxdir;	  
	    ftemp1+=theflux;
	    if(kever==0){
	      ftemp2+=theflux*v1tov3(v[1][1],k,j,i)*0.5;
	    }
	  }
	  else{
	    theflux=fl[dir][k][j][i];
	    theflux*=fluxdir;	  
	    ftemp1+=theflux;
	    if(kever==0){
	      ftemp2+=theflux*z2e_3(v[1][1],k,j,i)*0.5;
	    }
	  }
	}
      }
      else if(which==6){ // flux of vy momentum
	if(dir==1){
	  if(USEADVFLUX==0){ // on edge
	    theflux=dt*G2(1,i)*z2e_1(s[1],k,j,i)*v[1][1][k][j][i]*v2tov1(v[1][2],k,j,i)*G2(1,i)*dx[1][2][j]*G3(1,i)*G4(2,j)*dx[1][3][k];
	    theflux*=fluxdir;
	    ftemp1+=theflux;
	    if(kever==0){
	      ftemp2+=theflux*(v2tov1(v[1][2],k,j,i)*0.5/G2(1,i));
	    }
	  }
	  else{
	    theflux=fl[dir][k][j][i]*dx[1][3][k];
	    theflux*=fluxdir;
	    ftemp1+=theflux;
	    if(kever==0){
	      ftemp2+=theflux*(z2e_1(v[1][2],k,j,i)*0.5/G2(1,i));
	    }
	  }
	}
	else if(dir==2){	    
	  if(1||(USEADVFLUX==0)){ // on edge
	    theflux=dt*(z2e_2(s[1],k,j,i)*v[1][2][k][j][i]*v[1][2][k][j][i]+(gam-1.0)*z2e_2(s[2],k,j,i)+z2e_2(s[1],k,j,i)*z2e_2(s[3],k,j,i))*G2(2,i)*dx[1][1][i]*G3(2,i)*G4(1,j)*dx[1][3][k];
	    theflux*=fluxdir;
	    ftemp1+=theflux; // total flux
	    if(kever==0){
	      ftemp2+=theflux*v[1][2][k][j][i]*0.5/G2(2,i);
	    }
	  }
	  else{
	    inftemp1=fl[dir][k][j][i]*dx[1][3][k]; // not on boundary!
	    inftemp1+=((gam-1.0)*s[2][k][j][i]+s[1][k][j][i]*s[3][k][j][i])*G2(2,i)*G3(2,i)*G4(2,j)*dx[1][3][k]*dx[1][1][i]; // p+rho*phi flux
	    theflux=inftemp1;
	    theflux*=fluxdir;
	    ftemp1+=theflux; // total flux
	    if(kever==0){
	      ftemp2+=theflux*e2z_2(v[1][2],k,j,i)*0.5/G2(2,i);
	    }
	  }
	}
	else if(dir==3){
	  if(USEADVFLUX==0){ // on edge
	    theflux=dt*G2(2,i)*z2e_3(s[1],k,j,i)*v[1][3][k][j][i]*v2tov3(v[1][2],k,j,i)*dx[1][1][i]*G2(2,i)*dx[1][2][j];
	    theflux*=fluxdir;	  
	    ftemp1+=theflux;
	    if(kever==0){
	      ftemp2+=theflux*v2tov3(v[1][2],k,j,i)*0.5/G2(2,i);
	    }
	  }
	  else{
	    theflux=fl[dir][k][j][i];
	    theflux*=fluxdir;	  
	    ftemp1+=theflux;
	    if(kever==0){
	      ftemp2+=theflux*(z2e_3(v[1][2],k,j,i)*0.5/G2(2,i));
	    }
	  }
	}
      }
      else if(which==7){ // flux of vz momentum
	if(dir==1){
	  if(USEADVFLUX==0){ // on edge
	    theflux=dt*G3(1,i)*G4(2,j)*z2e_1(s[1],k,j,i)*v[1][1][k][j][i]*v3tov1(v[1][3],k,j,i)*G2(1,i)*dx[1][2][j]*G3(1,i)*G4(2,j)*dx[1][3][k];
	    theflux*=fluxdir;	  
	    ftemp1+=theflux;
	    // add in kinetic energy
	    if(kever==0){
	      ftemp2+=theflux*v3tov1(v[1][3],k,j,i)*0.5/(G3(1,i)*G4(2,j));
	    }
	  }
	  else{
	    theflux=fl[dir][k][j][i]*dx[2][3][k];
	    theflux*=fluxdir;	  
	    ftemp1+=theflux;
	    // add in kinetic energy
	    if(kever==0){
	      ftemp2+=theflux*(z2e_1(v[1][3],k,j,i)*0.5/(G3(1,i)*G4(2,j)));
	    }
	  }
	}
	else if(dir==2){
	  if(USEADVFLUX==0){ // on edge
	    theflux=dt*G3(2,i)*G4(1,j)*z2e_2(s[1],k,j,i)*v[1][2][k][j][i]*v3tov2(v[1][3],k,j,i)*dx[1][1][i]*G3(2,i)*G4(1,j)*dx[1][3][k];
	    theflux*=fluxdir;
	    ftemp1+=theflux;
	    // add in kinetic energy
	    if(kever==0){
	      ftemp2+=theflux*v3tov2(v[1][3],k,j,i)*0.5*OG4(1,j)/G3(2,i);
	    }
	  }
	  else{
	    theflux=fl[dir][k][j][i]*dx[2][3][k];
	    theflux*=fluxdir;
	    ftemp1+=theflux;
	    // add in kinetic energy
	    if(kever==0){
	      ftemp2+=theflux*(z2e_2(v[1][3],k,j,i)*0.5*OG4(1,j)/G3(2,i));
	    }
	  }
	}
	else if(dir==3){
	  if(1||(USEADVFLUX==0)){ // on edge
	    theflux=dt*(z2e_3(s[1],k,j,i)*v[1][3][k][j][i]*v[1][3][k][j][i]+(gam-1.0)*z2e_3(s[2],k,j,i)+z2e_3(s[1],k,j,i)*z2e_3(s[3],k,j,i))*G3(2,i)*G4(2,j)*dx[1][1][i]*G2(2,i)*dx[1][2][j];
	    theflux*=fluxdir;	  
	    ftemp1+=theflux; // total flux
	    if(kever==0){
	      ftemp2+=theflux*v[1][3][k][j][i]*0.5/(G3(2,i)*G4(2,j));
	    }
	  }
	  else{
	    inftemp1=fl[dir][k][j][i]; // not on boundary!
	    inftemp1+=((gam-1.0)*s[2][k][j][i]+s[1][k][j][i]*s[3][k][j][i])*G2(2,i)*G3(2,i)*G4(2,j)*dx[1][2][j]*dx[1][1][i]; // p+rho*phi flux
	    theflux=inftemp1;
	    theflux*=fluxdir;	  
	    ftemp1+=theflux; // total flux
	    if(kever==0){
	      ftemp2+=theflux*e2z_3(v[1][3],k,j,i)*0.5/(G3(2,i)*G4(2,j));
	    }
	  }
	}
      }
      // now assign value
      if(which==1){
	// mass loss
	lossflux[1][io]+=ftemp1;
	
	// grav. pot. energy loss
	lossflux[3][io]+=ftemp2;
	
	// ke flux
	if(kever==1){
	  lossflux[4][io]+=ftemp3;
	  lossflux[4][io]+=ftemp4;
	  lossflux[4][io]+=ftemp5;
	}
      }
      else if(which==2){
	// capture flux of enthalpy(only true for ideal gas EOS)
	lossflux[2][io]+=ftemp1;
      }
      else if(which==3){
	// capture magflux3 -- actually bz/(rsin(theta))
	lossflux[11][io]+=ftemp1;
      }
      else if(which==5){
	lossflux[5][io]+=ftemp1;
	if(kever==0){
	  lossflux[4][io]+=ftemp2;
	}
      }
      else if(which==6){
	lossflux[6][io]+=ftemp1;
	if(kever==0){
	  lossflux[4][io]+=ftemp2;
	}
      }
      else if(which==7){
	lossflux[7][io]+=ftemp1;
	if(kever==0){
	  lossflux[4][io]+=ftemp2;
	}
      }
    }
  }

} 



// general loop version, but compute only and all at once
// used when COMPUTELOSSDIAG==1 and LOOPTYPE>1

// note this isn't goint to preserve mass across all boundaries exactly, but still uses mdot (an unconserved form) as base for various quantities.
void hydro_flux_gen_gen(void) // hydro_flux()
{
  int i,j,k,m;
  int itemp;
  int looper;
  SFTYPE ftemp1,ftemp2,ftemp3,ftemp4,ftemp5;
  SFTYPE inftemp1,inftemp2;
  int io;
  SFTYPE theflux;
  SFTYPE fluxdir;
  int looperstart,looperend;
  SFTYPE ftempv[3+1];
  SFTYPE smom[3+1];
  

  looperstart=23;
  looperend=28;
  
  for(looper=looperstart;looper<=looperend;looper++){
    // enumerate surfaces
    if((looper==23)||(looper==25)||(looper==27)) io=0; else io=1;

    if( ((looper==23)||(looper==24))&&(N1==1) ) break;
    if( ((looper==25)||(looper==26))&&(N2==1) ) break;
    if( ((looper==27)||(looper==28))&&(N3==1) ) break;

    // loop over surfaces zones
    
    LOOPSUPERGEN(looper){
      ftemp1=0;
      ftemp2=0;
      ftemp3=0;
      ftemp4=0;
      ftempv[1]=ftempv[2]=ftempv[3]=0;
      
      fluxdir=bzs[looper][temptempi][0];

      // mass flux

      if((looper==23)||(looper==24)){
	// dir==1
	theflux=dt*z2e_1(s[1],k,j,i)*v[1][1][k][j][i]*G2(1,i)*G3(1,i)*G4(2,j)*dx[1][2][j]*dx[1][3][k];
	theflux*=fluxdir;
	
	ftemp1+=theflux;
	
	// grav. pot. energy flux
	ftemp3+=theflux*z2e_1(s[3],k,j,i);
	// ke flux
	if(kever==1){
	  ftemp4+=theflux*0.5*v[1][1][k][j][i]*v[1][1][k][j][i];
	  ftemp4+=theflux*0.5*v2tov1(v[1][2],k,j,i)*v2tov1(v[1][2],k,j,i);
	  ftemp4+=theflux*0.5*v3tov1(v[1][3],k,j,i)*v3tov1(v[1][3],k,j,i);
	}
      }
      if((looper==25)||(looper==26)){
	// dir==2
	theflux=dt*z2e_2(s[1],k,j,i)*v[1][2][k][j][i]*G3(2,i)*G4(1,j)*dx[1][1][i]*dx[1][3][k];
	theflux*=fluxdir;
	
	ftemp1+=theflux;
	// grav pot energy
	ftemp3+=theflux*z2e_2(s[3],k,j,i);
	
	if(kever==1){	
	  // see sweepx for comment on why this done here instead of above commented ke's
	  ftemp4+=theflux*0.5*v1tov2(v[1][1],k,j,i)*v1tov2(v[1][1],k,j,i);
	  ftemp4+=theflux*0.5*v[1][2][k][j][i]*v[1][2][k][j][i];
	  ftemp4+=theflux*0.5*v3tov2(v[1][3],k,j,i)*v3tov2(v[1][3],k,j,i);
	}
      }
      
      if((looper==27)||(looper==28)){
	// dir==3
	theflux=dt*z2e_3(s[1],k,j,i)*v[1][3][k][j][i]*G2(2,i)*dx[1][2][j]*dx[1][1][i];
	theflux*=fluxdir;
	
	ftemp1+=theflux;
	// grav. pot. energy flux
	ftemp3+=theflux*z2e_3(s[3],k,j,i);
	
	// ke flux
	if(kever==1){
	  ftemp4+=theflux*0.5*v1tov3(v[1][1],k,j,i)*v1tov3(v[1][1],k,j,i);
	  ftemp4+=theflux*0.5*v2tov3(v[1][2],k,j,i)*v2tov3(v[1][2],k,j,i);
	  ftemp4+=theflux*0.5*v[1][3][k][j][i]*v[1][3][k][j][i];
	}
      }
      // capture flux of enthalpy
      if((looper==23)||(looper==24)){
	// dir==1
	theflux=dt*gam*z2e_1(s[2],k,j,i)*v[1][1][k][j][i]*G2(1,i)*G3(1,i)*G4(2,j)*dx[1][2][j]*dx[1][3][k];
	theflux*=fluxdir;
	
	ftemp2+=theflux;
      }
      if((looper==25)||(looper==26)){
	// dir==2
	theflux=dt*gam*z2e_2(s[2],k,j,i)*v[1][2][k][j][i]*G3(2,i)*G4(1,j)*dx[1][1][i]*dx[1][3][k];
	theflux*=fluxdir;
	
	ftemp2+=theflux;
      }
      if((looper==27)||(looper==28)){
	// dir==3
	theflux=dt*gam*z2e_3(s[2],k,j,i)*v[1][3][k][j][i]*G2(2,i)*dx[1][2][j]*dx[1][1][i];
	theflux*=fluxdir;
	
	ftemp2+=theflux;
      }
      if((looper==23)||(looper==24)){
	// flux of vx momentum
	// dir==1
	theflux=dt*(z2e_1(s[1],k,j,i)*v[1][1][k][j][i]*v[1][1][k][j][i]+(gam-1.0)*z2e_1(s[2],k,j,i)+z2e_1(s[1],k,j,i)*z2e_1(s[3],k,j,i))*G2(1,i)*dx[1][2][j]*G3(1,i)*G4(2,j)*dx[1][3][k];
	theflux*=fluxdir;
	ftempv[1]+=theflux; // total flux
	if(kever==0){
	  ftemp4+=theflux*v[1][1][k][j][i]*0.5;
	}
      }
      if((looper==25)||(looper==26)){
	// dir==2
	theflux=dt*z2e_2(s[1],k,j,i)*v1tov2(v[1][1],k,j,i)*v[1][2][k][j][i]*dx[1][1][i]*G3(2,i)*G4(1,j)*dx[1][3][k];
	theflux*=fluxdir;
	ftempv[1]+=theflux;
	if(kever==0){
	  // ke loss
	  ftemp4+=theflux*v1tov2(v[1][1],k,j,i)*0.5;
	}
      }
      if((looper==27)||(looper==28)){
	// dir==3
	theflux=dt*z2e_3(s[1],k,j,i)*v1tov3(v[1][1],k,j,i)*v[1][3][k][j][i]*dx[1][1][i]*G2(2,i)*dx[1][2][j];
	theflux*=fluxdir;	  
	ftempv[1]+=theflux;
	if(kever==0){
	  ftemp4+=theflux*v1tov3(v[1][1],k,j,i)*0.5;
	}
      }
      // flux of vy momentum
      if((looper==23)||(looper==24)){
	// dir==1
	theflux=dt*G2(1,i)*z2e_1(s[1],k,j,i)*v[1][1][k][j][i]*v2tov1(v[1][2],k,j,i)*G2(1,i)*dx[1][2][j]*G3(1,i)*G4(2,j)*dx[1][3][k];
	theflux*=fluxdir;
	ftempv[2]+=theflux;
	if(kever==0){
	  ftemp4+=theflux*(v2tov1(v[1][2],k,j,i)*0.5/G2(1,i));
	}
      }
      if((looper==25)||(looper==26)){
	// dir==2
	theflux=dt*(z2e_2(s[1],k,j,i)*v[1][2][k][j][i]*v[1][2][k][j][i]+(gam-1.0)*z2e_2(s[2],k,j,i)+z2e_2(s[1],k,j,i)*z2e_2(s[3],k,j,i))*G2(2,i)*dx[1][1][i]*G3(2,i)*G4(1,j)*dx[1][3][k];
	theflux*=fluxdir;
	ftempv[2]+=theflux; // total flux
	if(kever==0){
	  ftemp4+=theflux*v[1][2][k][j][i]*0.5/G2(2,i);
	}
      }
      if((looper==27)||(looper==28)){
	// dir==3
	theflux=dt*G2(2,i)*z2e_3(s[1],k,j,i)*v[1][3][k][j][i]*v2tov3(v[1][2],k,j,i)*dx[1][1][i]*G2(2,i)*dx[1][2][j];
	theflux*=fluxdir;	  
	ftempv[2]+=theflux;
	if(kever==0){
	  ftemp4+=theflux*v2tov3(v[1][2],k,j,i)*0.5/G2(2,i);
	}
      }
      // flux of vz momentum
      if((looper==23)||(looper==24)){
	// dir==1
	theflux=dt*G3(1,i)*G4(2,j)*z2e_1(s[1],k,j,i)*v[1][1][k][j][i]*v3tov1(v[1][3],k,j,i)*G2(1,i)*dx[1][2][j]*G3(1,i)*G4(2,j)*dx[1][3][k];
	theflux*=fluxdir;	  
	ftempv[3]+=theflux;
	// add in kinetic energy
	if(kever==0){
	  ftemp4+=theflux*v3tov1(v[1][3],k,j,i)*0.5/(G3(1,i)*G4(2,j));
	}
      }
      if((looper==25)||(looper==26)){
	// dir==2
	theflux=dt*G3(2,i)*G4(1,j)*z2e_2(s[1],k,j,i)*v[1][2][k][j][i]*v3tov2(v[1][3],k,j,i)*dx[1][1][i]*G3(2,i)*G4(1,j)*dx[1][3][k];
	theflux*=fluxdir;
	ftempv[3]+=theflux;
	// add in kinetic energy
	if(kever==0){
	  ftemp4+=theflux*v3tov2(v[1][3],k,j,i)*0.5*OG4(1,j)/G3(2,i);
	}
      }
      if((looper==27)||(looper==28)){
	// dir==3
	theflux=dt*(z2e_3(s[1],k,j,i)*v[1][3][k][j][i]*v[1][3][k][j][i]+(gam-1.0)*z2e_3(s[2],k,j,i)+z2e_3(s[1],k,j,i)*z2e_3(s[3],k,j,i))*G3(2,i)*G4(2,j)*dx[1][1][i]*G2(2,i)*dx[1][2][j];
	theflux*=fluxdir;	  
	ftempv[3]+=theflux; // total flux
	if(kever==0){
	  ftemp4+=theflux*v[1][3][k][j][i]*0.5/(G3(2,i)*G4(2,j));
	}
      }
    
      // now assign value
      // mass loss
      lossflux[1][io]+=ftemp1;

      // capture flux of enthalpy(only true for ideal gas EOS)
      lossflux[2][io]+=ftemp2;
      
      // grav. pot. energy loss
      lossflux[3][io]+=ftemp3;
      
      // ke flux
      lossflux[4][io]+=ftemp4;

      // s mom flux (mixes directions and so positions)
      lossflux[5][io]+=ftempv[1];
      lossflux[6][io]+=ftempv[2];
      lossflux[7][io]+=ftempv[3]; // ang mom about z-axis if COORD==3

      // general ang mom flux (mixes directions even more by assuming located at center when not really)
#if(COORD==1)
      smom[1]=ftempv[1];
      smom[2]=ftempv[2];
      smom[3]=ftempv[3];

      lossflux[13][io]+=sqrt(smom[2]*smom[2]*x[2][3][k]*x[2][3][k]+smom[3]*smom[3]*x[2][2][j]*x[2][2][j]);
      lossflux[14][io]+=sqrt(smom[1]*smom[1]*x[2][3][k]*x[2][3][k]+smom[3]*smom[3]*x[2][1][i]*x[2][1][i]);
      lossflux[15][io]+=sqrt(smom[2]*smom[2]*x[2][1][i]*x[2][1][i]+smom[1]*smom[1]*x[2][2][j]*x[2][2][j]);
#endif
    }
  }

} 



// reconsider reflecting conditions here

// in 3D the locations are:
// (x,y,z)
// sig11,22,33: (.5,.5,.5)
// sig13: (0,.5,0)
// sig12: (0,0,.5)
// sig23: (.5,0,0)
void viscous_flux_gen(void){
  // right before going to advection, compute fluxes of energy and momentum due to viscosity
  
  int i,j,k,m;
  int itemp;
  SFTYPE ftemp1,ftemp2,ftemp3;
  SFTYPE ftempv[3+1];
  int io;
  int dir;
  int looper;
  SFTYPE theflux;
  SFTYPE fluxdir;
  int looperstart,looperend;
  FTYPE subftemp;
  FTYPE Length;
  FTYPE flux;
  FTYPE crapf;
  SFTYPE smom[3+1];


  // do all in 1 call
  looperstart=23;
  looperend=28;

  for(looper=looperstart;looper<=looperend;looper++){
    // enumerate surfaces
    if((looper==23)||(looper==25)||(looper==27)) io=0; else io=1;

    if( ((looper==23)||(looper==24))&&(N1==1) ) break;
    if( ((looper==25)||(looper==26))&&(N2==1) ) break;
    if( ((looper==27)||(looper==28))&&(N3==1) ) break;

    // loop over surfaces zones

    LOOPSUPERGEN(looper){
      ftemp1=ftemp2=ftemp3=0;
      ftempv[1]=ftempv[2]=ftempv[3]=0;

      fluxdir=bzs[looper][temptempi][0];

      if(vischeat){
	// capture flux of internal energy
	if((looper==23)||(looper==24)){
	  if(VISCE31){	    // dir==1
	    ftemp1+=fluxdir*dt*sigma[3][1][k][j][i]*z2e_1(v[1][3],k,j,i)*DS11(k,j,i)*dx[2][3][k];
	  }
	}
	if((looper==27)||(looper==28)){
	  if(VISCE13){      // dir==3
	    ftemp3+=fluxdir*dt*sigma[1][3][k][j][i]*z2e_3(v[1][1],k,j,i)*DS23(k,j,i);
	  }
	}
	if((looper==25)||(looper==26)){
	  if(VISCE23){      // dir==2
	    ftemp2+=fluxdir*dt*sigma[2][3][k][j][i]*z2e_2(v[1][3],k,j,i)*DS12(k,j,i)*dx[2][3][k];
	  }
	}
	if((looper==27)||(looper==28)){
	  if(VISCE32){      // dir==3
	    ftemp3+=fluxdir*dt*sigma[3][2][k][j][i]*z2e_3(v[1][2],k,j,i)*DS33(k,j,i);
	  }
	}
	if((looper==23)||(looper==24)){
	  if(VISCE11){      // dir==1
	    ftemp1+=fluxdir*dt*z2e_1(sigma[1][1],k,j,i)*v[1][1][k][j][i]*DS11(k,j,i)*dx[1][3][k];
	  }
	}
	if((looper==25)||(looper==26)){
	  if(VISCE22){      // dir==2
	    ftemp2+=fluxdir*dt*z2e_2(sigma[2][2],k,j,i)*v[1][2][k][j][i]*DS12(k,j,i)*dx[1][3][k];
	  }
	}
	if((looper==27)||(looper==28)){
	  if(VISCE33){      // dir==3
	    ftemp3+=fluxdir*dt*z2e_3(sigma[3][3],k,j,i)*v[1][3][k][j][i]*DS13(k,j,i);
	  }
	}
	if((looper==23)||(looper==24)){
	  if(VISCE12){      // dir==1
	    ftemp1+=fluxdir*dt*sigma[1][2][k][j][i]*z2e_1(v[1][2],k,j,i)*DS31(k,j,i)*dx[1][3][k];
	  }
	}
	if((looper==25)||(looper==26)){
	  if(VISCE21){ // dir==2
	    ftemp2+=fluxdir*dt*sigma[2][1][k][j][i]*z2e_2(v[1][1],k,j,i)*DS22(k,j,i)*dx[1][3][k];
	  }
	}
      }
      // momentum fluxes (without volume terms)
      ///////////
      // momentum fluxes due to viscosity
      // in vector positions
      // just copy of magnetic version and stick in sigma in right vector location
      // since interp on sigma sometimes, assume sigma exists at least on LOOPH
      if((looper==23)||(looper==24)){
	if(VISCE11){ // dir==1 (vx)
	  ftempv[1]+=fluxdir*(z2e_1(sigma[1][1],k,j,i))*DS11(k,j,i)*dx[1][3][k]*dt;
	}
      }
      if((looper==25)||(looper==26)){
	if(VISCE21){ // dir==2 (vx)
	  ftempv[1]+=fluxdir*(e2z_1(sigma[2][1],k,j,i))*DS12(k,j,i)*dx[1][3][k]*dt;
	}
      }
      if((looper==27)||(looper==28)){
	if(VISCE31){ // dir==3 (vx)
	  ftempv[1]+=fluxdir*(e2z_1(sigma[3][1],k,j,i))*DS13(k,j,i)*dt;
	}
      }
      if((looper==23)||(looper==24)){
	if(VISCE12){ // dir==1 (vy)
	  ftempv[2]+=fluxdir*(e2z_2(sigma[1][2],k,j,i))*G2(1,i)*DS11(k,j,i)*dx[1][3][k]*dt; // g12i is from absorption of one volume term
	}
      }
      if((looper==25)||(looper==26)){
	if(VISCE22){ // dir==2 (vy)
	  ftempv[2]+=fluxdir*(z2e_2(sigma[2][2],k,j,i))*G2(2,i)*DS12(k,j,i)*dx[1][3][k]*dt; // G2(2,i) from definition of s2=h2*rho*v2
	}
      }
      if((looper==27)||(looper==28)){
	if(VISCE32){ // dir==3 (vy)
	  ftempv[2]+=fluxdir*(e2z_2(sigma[3][2],k,j,i))*G2(2,i)*DS13(k,j,i)*dt; // G2(2,i) from definition of s2=h2*rho*v2
	}
      }
      if((looper==23)||(looper==24)){
	if(VISCE13){ // dir==1 (vz)
	  ftempv[3]+=fluxdir*(e2z_3(sigma[1][3],k,j,i))*G3(1,i)*G4(2,j)*DS11(k,j,i)*dx[1][3][k]*dt; // g13i*g24j is from absorbed volume term
	}
      }
      if((looper==25)||(looper==26)){
	if(VISCE23){ // dir==2 (vz)
	  ftempv[3]+=fluxdir*(e2z_3(sigma[2][3],k,j,i))*G3(2,i)*G4(1,j)*DS12(k,j,i)*dx[1][3][k]*dt;// g23i*g14j is from absorbed volume term
	}
      }
      if((looper==27)||(looper==28)){
	if(VISCE33){ // dir==3 (vz)
	  ftempv[3]+=fluxdir*(z2e_3(sigma[3][3],k,j,i))*G3(2,i)*G4(2,j)*DS13(k,j,i)*dt; // g23i*g24j is from definition of s3=g3*g4*v3 
	}
      }
      // visc heating
      lossflux[12][io]+=ftemp1;
      lossflux[12][io]+=ftemp2;
      lossflux[12][io]+=ftemp3;
      // momentum fluxes
      lossflux[5][io]+=ftempv[1];
      lossflux[6][io]+=ftempv[2];
      lossflux[7][io]+=ftempv[3];

      // general ang mom flux (mixes directions even more by assuming located at center when not really)
#if(COORD==1)
      smom[1]=ftempv[1];
      smom[2]=ftempv[2];
      smom[3]=ftempv[3];

      lossflux[13][io]+=sqrt(smom[2]*smom[2]*x[2][3][k]*x[2][3][k]+smom[3]*smom[3]*x[2][2][j]*x[2][2][j]);
      lossflux[14][io]+=sqrt(smom[1]*smom[1]*x[2][3][k]*x[2][3][k]+smom[3]*smom[3]*x[2][1][i]*x[2][1][i]);
      lossflux[15][io]+=sqrt(smom[2]*smom[2]*x[2][1][i]*x[2][1][i]+smom[1]*smom[1]*x[2][2][j]*x[2][2][j]);
#endif

    }
  }
  
}


// all magnetic fluxes are on vector positions
// very similar to viscous flux in form for momentum terms

// energy, momentum, and magnetic
void magnetic_flux_gen(void)
{
  FTYPE b2,vdotb,bbtensor,emf,ftemp;
  int i,j,k,m;
  int dir,io;
  int looper;
  int itemp;
  int looperstart,looperend;
  SFTYPE ftemp1,ftemp2,ftemp3;
  SFTYPE ftempv[3+1],ftempb[3+1];
  SFTYPE theflux;
  SFTYPE fluxdir;
  SFTYPE smom[3+1];

  // done all together in readable form since not in flux conservative form, unlike viscous or hydro flux

  // do all in 1 call
  looperstart=23;
  looperend=28;


  for(looper=looperstart;looper<=looperend;looper++){
    // enumerate surfaces
    if((looper==23)||(looper==25)||(looper==27)) io=0; else io=1;

    if( ((looper==23)||(looper==24))&&(N1==1) ) break;
    if( ((looper==25)||(looper==26))&&(N2==1) ) break;
    if( ((looper==27)||(looper==28))&&(N3==1) ) break;

    // loop over surfaces zones

    LOOPSUPERGEN(looper){
      ftemp1=ftemp2=ftemp3=0;
      ftempv[1]=ftempv[2]=ftempv[3]=0;
      ftempb[1]=ftempb[2]=ftempb[3]=0;

      fluxdir=bzs[looper][temptempi][0];

      ////////////////////////////
      // magnetic flux of energy
      
      // dir==1 (be-x1)
      if((looper==23)||(looper==24)){
	vdotb=v[1][1][k][j][i]*v[2][1][k][j][i]+v2tov1(v[1][2],k,j,i)*v2tov1(v[2][2],k,j,i)+v3tov1(v[1][3],k,j,i)*v3tov1(v[2][3],k,j,i);
	b2=v[2][1][k][j][i]*v[2][1][k][j][i]+v2tov1(v[2][2],k,j,i)*v2tov1(v[2][2],k,j,i)+v3tov1(v[2][3],k,j,i)*v3tov1(v[2][3],k,j,i);
	
	ftemp1+=fluxdir*(b2*v[1][1][k][j][i]-vdotb*v[2][1][k][j][i])*DS11(k,j,i)*dx[1][3][k]*dt;
      }
      // dir==2 (be-x2)
      if((looper==25)||(looper==26)){
	vdotb=v1tov2(v[1][1],k,j,i)*v1tov2(v[2][1],k,j,i)+v[1][2][k][j][i]*v[2][2][k][j][i]+v3tov2(v[1][3],k,j,i)*v3tov2(v[2][3],k,j,i);
	b2=v1tov2(v[2][1],k,j,i)*v1tov2(v[2][1],k,j,i)+v[2][2][k][j][i]*v[2][2][k][j][i]+v3tov2(v[2][3],k,j,i)*v3tov2(v[2][3],k,j,i);
	
	ftemp2+=fluxdir*(b2*v[1][2][k][j][i]-vdotb*v[2][2][k][j][i])*DS12(k,j,i)*dx[1][3][k]*dt;
      }
      // dir==3 (be-x3)
      if((looper==27)||(looper==28)){
	vdotb=v1tov3(v[1][1],k,j,i)*v1tov3(v[2][1],k,j,i)+v2tov3(v[1][2],k,j,i)*v2tov3(v[2][2],k,j,i)+v[1][3][k][j][i]*v[2][3][k][j][i];
	b2=v1tov3(v[2][1],k,j,i)*v1tov3(v[2][1],k,j,i)+v2tov3(v[2][2],k,j,i)*v2tov3(v[2][2],k,j,i)+v[2][3][k][j][i]*v[2][3][k][j][i];
	
	ftemp3+=fluxdir*(b2*v[1][3][k][j][i]-vdotb*v[2][3][k][j][i])*DS13(k,j,i)*dt;
      }
      /////////////////////
      // momentum fluxes due to magnetic field (without volume terms)
      
      // vx-x1 flux (dir==1)
      if((looper==23)||(looper==24)){
	bbtensor=-v[2][1][k][j][i]*v[2][1][k][j][i];
	b2=v[2][1][k][j][i]*v[2][1][k][j][i]+v2tov1(v[2][2],k,j,i)*v2tov1(v[2][2],k,j,i)+v3tov1(v[2][3],k,j,i)*v3tov1(v[2][3],k,j,i);
	
	ftempv[1]+=fluxdir*(bbtensor+0.5*b2)*DS11(k,j,i)*dx[1][3][k]*dt;
      }
      // vx-x2 flux (dir==2)
      if((looper==25)||(looper==26)){
	bbtensor=-v1tov2(v[2][1],k,j,i)*v[2][2][k][j][i];
	//      b2=0;
	
	ftempv[1]+=fluxdir*(bbtensor)*DS12(k,j,i)*dx[1][3][k]*dt;
      }
      // vx-x3 flux (dir==3)
      if((looper==27)||(looper==28)){
	bbtensor=-v1tov3(v[2][1],k,j,i)*v[2][3][k][j][i];
	//      b2=0;
	
	ftempv[1]+=fluxdir*(bbtensor)*DS13(k,j,i)*dt;
      }
      // vy-x1 (dir==1)
      if((looper==23)||(looper==24)){
	bbtensor=-v[2][1][k][j][i]*v2tov1(v[2][2],k,j,i);
	b2=0;
	
	ftempv[2]+=fluxdir*(bbtensor)*G2(1,i)*DS11(k,j,i)*dx[1][3][k]*dt; // g12i is from absorption of one volume term
      }
      // vy-x2 (dir==2)
      if((looper==25)||(looper==26)){
	bbtensor=-v[2][2][k][j][i]*v[2][2][k][j][i];
	b2=v1tov2(v[2][1],k,j,i)*v1tov2(v[2][1],k,j,i)+v[2][2][k][j][i]*v[2][2][k][j][i]+v3tov2(v[2][3],k,j,i)*v3tov2(v[2][3],k,j,i);
	
	ftempv[2]+=fluxdir*(bbtensor+0.5*b2)*G2(2,i)*DS12(k,j,i)*dx[1][3][k]*dt; // G2(2,i) from definition of s2=h2*rho*v2
      }
      // vy-x3 (dir==3)
      if((looper==27)||(looper==28)){
	bbtensor=-v[2][3][k][j][i]*v2tov3(v[2][2],k,j,i);
	//b2=0;
	
	ftempv[2]+=fluxdir*(bbtensor)*G2(2,i)*DS13(k,j,i)*dt; // G2(2,i) from definition of s2=h2*rho*v2
      }
      // vz-x1 (dir==1)
      if((looper==23)||(looper==24)){
	bbtensor=-v[2][1][k][j][i]*v3tov1(v[2][3],k,j,i);
	//b2=0;
	
	ftempv[3]+=fluxdir*(bbtensor)*G3(1,i)*G4(2,j)*DS11(k,j,i)*dx[1][3][k]*dt; // g13i*g24j is from absorbed volume term
      }
      // vz-x2 (dir==2)
      if((looper==25)||(looper==26)){
	bbtensor=-v[2][2][k][j][i]*v3tov2(v[2][3],k,j,i);
	//b2=0;
	
	ftempv[3]+=fluxdir*(bbtensor)*G3(2,i)*G4(1,j)*DS12(k,j,i)*dx[1][3][k]*dt;// g23i*g14j is from absorbed volume term
      }
      // vz-x3 (dir==3)
      if((looper==27)||(looper==28)){
	bbtensor=-v[2][3][k][j][i]*v[2][3][k][j][i];
	b2=v1tov3(v[2][1],k,j,i)*v1tov3(v[2][1],k,j,i)+v2tov3(v[2][2],k,j,i)*v2tov3(v[2][2],k,j,i)+v[2][3][k][j][i]*v[2][3][k][j][i];
	
	ftempv[3]+=fluxdir*(bbtensor+0.5*b2)*G3(2,i)*G4(2,j)*DS13(k,j,i)*dt; // g23i*g24j is from definition of s3=g3*g4*v3 
      }
      /////////////////////
      // magnetic fluxes
      
      // bx-x1 flux (dir==1) // is 0 by definition
      
      // bx-x2 flux (dir==2)
      if((looper==25)||(looper==26)){
	// emf is really -emf, so that sign is same as usual
	emf=-(v1tov2(v[1][1],k,j,i)*v[2][2][k][j][i]-v[1][2][k][j][i]*v1tov2(v[2][1],k,j,i)); // -emf_z
	ftempb[1]+=fluxdir*(emf)*DS12(k,j,i)*dx[1][3][k]*dt;
      }
	// bx-x3 flux (dir==3)
      if((looper==27)||(looper==28)){
	emf=v[1][3][k][j][i]*v1tov3(v[2][1],k,j,i)-v1tov3(v[1][1],k,j,i)*v[2][3][k][j][i]; // emf_y
	ftempb[1]+=fluxdir*(emf)*DS13(k,j,i)*dt;
      }
      // by/g12 -x1 flux (dir==1)
      if((looper==23)||(looper==24)){
	emf=(v[1][1][k][j][i]*v2tov1(v[2][2],k,j,i)-v2tov1(v[1][2],k,j,i)*v[2][1][k][j][i]); // emf_z
	ftempb[2]+=fluxdir*(emf)*DS11(k,j,i)*dx[1][3][k]*dt/G2(1,i); // g12i since different for the field
      }
      // by/g22 - x2 flux is 0 by definition
      
      // by/g22 - x3 flux
      if((looper==27)||(looper==28)){
	emf=-(v2tov3(v[1][2],k,j,i)*v[2][3][k][j][i]-v[1][3][k][j][i]*v2tov3(v[2][2],k,j,i)); // -emf_x
	ftempb[2]+=fluxdir*(emf)*DS13(k,j,i)*dt/G2(2,i); // g12i since different for the field
      }
      if((HYDROBZFLUX==0)||(COMPDIM>2)){
	// bz/g13/g42 - x1 flux
	if((looper==23)||(looper==24)){
	  emf=-(-v[1][1][k][j][i]*v3tov1(v[2][3],k,j,i)+v3tov1(v[1][3],k,j,i)*v[2][1][k][j][i]); // -emf_y
	  ftempb[3]+=fluxdir*(emf)*DS11(k,j,i)*dx[1][3][k]*dt/(G3(1,i)*G4(2,j)); // g13i*g24j is since field different
	}
	// bz/g13/g42 - x2 flux
	if((looper==25)||(looper==26)){
	  emf=(v[1][2][k][j][i]*v3tov2(v[2][3],k,j,i)-v3tov2(v[1][3],k,j,i)*v[2][2][k][j][i]);
	  ftempb[3]+=fluxdir*(emf)*DS12(k,j,i)*dx[1][3][k]*dt*OG4(1,j)/G3(2,i);// g23i*g14j is since field different
	}
	// bz - x3 flux is 0 by definition
      }
      
      // momentum fluxes
      lossflux[5][io]+=ftempv[1];
      lossflux[6][io]+=ftempv[2];
      lossflux[7][io]+=ftempv[3];
      // magnetic energy flux
      lossflux[8][io]+=ftemp1;
      lossflux[8][io]+=ftemp2;
      lossflux[8][io]+=ftemp3;

      // general ang mom flux (mixes directions even more by assuming located at center when not really)
#if(COORD==1)
      smom[1]=ftempv[1];
      smom[2]=ftempv[2];
      smom[3]=ftempv[3];

      lossflux[13][io]+=sqrt(smom[2]*smom[2]*x[2][3][k]*x[2][3][k]+smom[3]*smom[3]*x[2][2][j]*x[2][2][j]);
      lossflux[14][io]+=sqrt(smom[1]*smom[1]*x[2][3][k]*x[2][3][k]+smom[3]*smom[3]*x[2][1][i]*x[2][1][i]);
      lossflux[15][io]+=sqrt(smom[2]*smom[2]*x[2][1][i]*x[2][1][i]+smom[1]*smom[1]*x[2][2][j]*x[2][2][j]);
#endif

      // magnetic fluxes
      lossflux[9][io]+=ftempb[1];
      lossflux[10][io]+=ftempb[2];
      lossflux[11][io]+=ftempb[3];
    }
  }
  
}

