
if(transpassivex2){
  // add passive scalar to inside accretor
  // must come after "field" set
#if(DO_ACCRETOR)
  if(passivepos==CENT){
    bound_accretor(2,v[2][1],NULL,-3,0,0); // direction (first argument) doesn't matter
  }
  else{
    if(passivepos!=-1) bound_accretor(2,NULL,v[2],0,-4,1); // direction (first argument) doesn't matter
    else bound_accretor(2,NULL,v[2],0,-4,123);
  }
#endif
}



workv=(FTYPE (*) [N3M][N2M][N1M])(&u[0][0][0]);



if(transpassivex2&&(passivepos==CENT)){
  LOOPFC{
    u[k][j][i] = v[2][1][k][j][i]/s[1][k][j][i] ;
  }
      
  dqy_calc(u,dq) ;
      
  // fl is really fl=Flux*dt/dx (because of mdot)
  LOOPT1j{
#if(SYMFORCEHD==0)
    if(vp[2][k][j][i] > 0.) 
      fl[2][k][j][i] = (u[k][jm1][i] + (dx[2][2][j]*G2(2,i) - vp[2][k][j][i])*dq[k][jm1][i])*mdot[2][k][j][i]*DS12(k,j,i);
    else
      fl[2][k][j][i] = (u[k][j][i] + (-dx[2][2][j]*G2(2,i) - vp[2][k][j][i])*dq[k][j][i])*mdot[2][k][j][i]*DS12(k,j,i);
#else
    Dtot=vp[2][k][j][i]*OARC12(k,j,i);
    if(Dtot > DTOTLIMIT) 
      fl[2][k][j][i] = (u[k][jm1][i] + (1.0 - Dtot)*dq[k][jm1][i]*dx[2][2][j]*G2(2,i))*mdot[2][k][j][i]*DS12(k,j,i);
    else if(Dtot< -DTOTLIMIT)
      fl[2][k][j][i] = (u[k][j][i] + (-1.0 - Dtot)*dq[k][j][i]*dx[2][2][j]*G2(2,i))*mdot[2][k][j][i]*DS12(k,j,i);
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      ftemp=mdot[2][k][j][i]*DS12(k,j,i);
      ftemp2=dx[2][2][j]*G2(2,i);
      dqminus=(u[k][jm1][i] + (1.0 - DTOTLIMIT)*dq[k][jm1][i]*ftemp2)*ftemp;
      dqplus=(u[k][j][i] + (-1.0 + DTOTLIMIT)*dq[k][j][i]*ftemp2)*ftemp;
      fl[2][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
#endif
  }
      
  LOOPC{
    v[2][1][k][j][i] += (fl[2][k][j][i]-fl[2][k][jp1][i])*OVOL1(k,j,i);
  }
  bound(v[2][1],NULL,-3,0,0) ;
}










if(transpassivev2x2&&( (passivepos==VDIR2)||(passivepos==-1))){

  if(passivepos==VDIR2)  wcom=1;
  else wcom=2;


  LOOPF3 LOOPFP12 LOOPF1{
    u[k][j][i] = v[2][wcom][k][j][i]/z2e_2(s[1],k,j,i);
  }

    
  /* vy */
  dqvy_calc(-2,workv,dqv) ;

    
  LOOPT2j{
    vpp = e2z_2(vp[2],k,j,i);
    mdotp = e2z_2(mdot[2],k,j,i);
#if(SYMFORCEHD==0)
    if(vpp > 0)
      fl[2][k][j][i] = (u[k][j][i] + (dx[1][2][j]*G2(2,i) - vpp)*dqv[2][k][j][i])*G2(2,i)*mdotp*DS32(k,j,i);
    else
      fl[2][k][j][i] = (u[k][jp1][i] + (-dx[1][2][j]*G2(2,i) - vpp)*dqv[2][k][jp1][i])*G2(2,i)*mdotp*DS32(k,j,i); // dqv[2][jp1] here needs u[k][j=N2+2][i] if periodicx2special.
#else
    Dtot=vpp*OARC32(k,j,i);
    if(Dtot > DTOTLIMIT)
      fl[2][k][j][i] = (u[k][j][i] + (1.0 - Dtot)*dqv[2][k][j][i]*dx[1][2][j]*G2(2,i))*G2(2,i)*mdotp*DS32(k,j,i);
    else if(Dtot< -DTOTLIMIT)
      fl[2][k][j][i] = (u[k][jp1][i] + (-1.0 - Dtot)*dqv[2][k][jp1][i]*dx[1][2][j]*G2(2,i))*G2(2,i)*mdotp*DS32(k,j,i);
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      ftemp=G2(2,i)*mdotp*DS32(k,j,i);
      ftemp2=dx[1][2][j]*G2(2,i);
      dqminus=(u[k][j][i] + (1.0 - DTOTLIMIT)*dqv[2][k][j][i]*ftemp2)*ftemp;
      dqplus=(u[k][jp1][i] + (-1.0 + DTOTLIMIT)*dqv[2][k][jp1][i]*ftemp2)*ftemp;
      fl[2][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
#endif
  }

  if(periodicx2special){ // instead of special dq code and vy variable on j=N2+2, or otherwise
    bound(NULL,fl,0,-4,2);
  }

  LOOPV2{
    v[2][wcom][k][j][i] += (fl[2][k][jm1][i]-fl[2][k][j][i])*OVOL3(k,j,i);
  }
  bound(NULL,v[2],0,-4,wcom) ; // bound passive scalar
}




if(transpassivev1x2&&( (passivepos==VDIR1)||(passivepos==-1))){

  if(passivepos==VDIR1)  wcom=1;
  else wcom=1;


  LOOPF3 LOOPF2 LOOPFP11{
    u[k][j][i] = v[2][wcom][k][j][i]/z2e_1(s[1],k,j,i);
  }


  /* vx */
  dqvy_calc(-1,workv,dqv) ;
    
  LOOPT3j{
    // below avgs reason for loopt1j over other bzones
    mdotp = z2e_1(mdot[2],k,j,i);
    vpp = z2e_1(vp[2],k,j,i);
#if(SYMFORCEHD==0)
    if(vpp > 0.) 
      fl[2][k][j][i] = (u[k][jm1][i] + (dx[2][2][j]*G2(1,i) - vpp)*dqv[1][k][jm1][i])*mdotp*DS22(k,j,i);
    else	
      fl[2][k][j][i] = (u[k][j][i] + (-dx[2][2][j]*G2(1,i) - vpp)*dqv[1][k][j][i])*mdotp*DS22(k,j,i);
#else
    Dtot=vpp*OARC12(k,j,i);
    if(Dtot > DTOTLIMIT) 
      fl[2][k][j][i] = (u[k][jm1][i] + (1.0 - Dtot)*dqv[1][k][jm1][i]*dx[2][2][j]*G2(1,i))*mdotp*DS22(k,j,i);
    else if(Dtot < -DTOTLIMIT)
      fl[2][k][j][i] = (u[k][j][i] + (-1.0 - Dtot)*dqv[1][k][j][i]*dx[2][2][j]*G2(1,i))*mdotp*DS22(k,j,i);
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      ftemp=mdotp*DS22(k,j,i);
      ftemp2=dx[2][2][j]*G2(1,i);
      dqminus=(u[k][jm1][i] + (1.0 - DTOTLIMIT)*dqv[1][k][jm1][i]*ftemp2)*ftemp;
      dqplus=(u[k][j][i] + (-1.0 + DTOTLIMIT)*dqv[1][k][j][i]*ftemp2)*ftemp;
      fl[2][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
#endif
  }
    
  LOOPV1{
    v[2][wcom][k][j][i] += (fl[2][k][j][i]-fl[2][k][jp1][i])*OVOL2(k,j,i);
  }
  bound(NULL,v[2],0,-4,wcom) ; // bound passive scalar
}




if(transpassivev3x2&&( (passivepos==VDIR3)||(passivepos==-1))){

  if(passivepos==VDIR3)  wcom=1;
  else wcom=3;

  LOOPFP13 LOOPF2 LOOPF1{
    u[k][j][i] = v[2][wcom][k][j][i]/z2e_3(s[1],k,j,i);
  }
    
  /* vz */
  dqvy_calc(-3,workv,dqv) ;
    
  LOOPT4j{
    vpp = z2e_3(vp[2],k,j,i);
    mdotp = z2e_3(mdot[2],k,j,i);
#if(SYMFORCEHD==0)
    if(vpp > 0) 
      fl[2][k][j][i] = (u[k][jm1][i] + (dx[2][2][j]*G2(2,i) - vpp)*dqv[3][k][jm1][i])*G3(2,i)*G4(1,j)*mdotp*DS12(k,j,i);
    else
      fl[2][k][j][i] = (u[k][j][i] + (-dx[2][2][j]*G2(2,i) - vpp)*dqv[3][k][j][i])*G3(2,i)*G4(1,j)*mdotp*DS12(k,j,i);
#else
    Dtot = vpp*OARC12(k,j,i);
    if(Dtot > DTOTLIMIT) 
      fl[2][k][j][i] = (u[k][jm1][i] + (1.0 - Dtot)*dqv[3][k][jm1][i]*dx[2][2][j]*G2(2,i))*G3(2,i)*G4(1,j)*mdotp*DS12(k,j,i);
    else if(Dtot< -DTOTLIMIT)
      fl[2][k][j][i] = (u[k][j][i] + (-1.0 - Dtot)*dqv[3][k][j][i]*dx[2][2][j]*G2(2,i))*G3(2,i)*G4(1,j)*mdotp*DS12(k,j,i);
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      ftemp=G3(2,i)*G4(1,j)*mdotp*DS12(k,j,i);
      ftemp2=dx[2][2][j]*G2(2,i);
      dqminus=(u[k][jm1][i] + (1.0 - DTOTLIMIT)*dqv[3][k][jm1][i]*ftemp2)*ftemp;
      dqplus=(u[k][j][i] + (-1.0 + DTOTLIMIT)*dqv[3][k][j][i]*ftemp2)*ftemp;
      fl[2][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
#endif
  }
    
  LOOPV3{
    v[2][wcom][k][j][i] += (fl[2][k][j][i]-fl[2][k][jp1][i])*OVOL1(k,j,i);
  }
  bound(NULL,v[2],0,-4,wcom) ; // bound passive scalar
}


