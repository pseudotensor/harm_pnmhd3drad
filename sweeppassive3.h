
if(transpassivex3){
  // add passive scalar to inside accretor
  // must come after "field" set
#if(DO_ACCRETOR)
  if(passivepos==CENT){
    bound_accretor(3,v[2][1],NULL,-3,0,0); // direction (first argument) doesn't matter
  }
  else{
    if(passivepos!=-1) bound_accretor(3,NULL,v[2],0,-4,1); // direction (first argument) doesn't matter
    else bound_accretor(3,NULL,v[2],0,-4,123);
  }
#endif
}


workv=(FTYPE (*) [N3M][N2M][N1M])(&u[0][0][0]);



if(transpassivex3&&(passivepos==CENT)){
  LOOPFC{
    u[k][j][i] = v[2][1][k][j][i]/s[1][k][j][i] ;
  }
      
  dqz_calc(u,dq) ;
      
  // fl is really fl=Flux*dt/dx (because of mdot)
  LOOPT1k{
#if(SYMFORCEHD==0)
    if(vp[3][k][j][i] > 0.) 
      fl[3][k][j][i] = (u[km1][j][i] + (dx[2][3][k]*G3(2,i)*G4(2,j) - vp[3][k][j][i])*dq[km1][j][i])*mdot[3][k][j][i]*DS13(k,j,i);
    else
      fl[3][k][j][i] = (u[k][j][i] + (-dx[2][3][k]*G3(2,i)*G4(2,j) - vp[3][k][j][i])*dq[k][j][i])*mdot[3][k][j][i]*DS13(k,j,i);
#else
    Dtot=vp[3][k][j][i]*OARC13(k,j,i)*ODX(2,3,k);
    if(Dtot > DTOTLIMIT) 
      fl[3][k][j][i] = (u[km1][j][i] + (1.0 - Dtot)*dq[km1][j][i]*dx[2][3][k]*G3(2,i)*G4(2,j))*mdot[3][k][j][i]*DS13(k,j,i);
    else if(Dtot< -DTOTLIMIT)
      fl[3][k][j][i] = (u[k][j][i] + (-1.0 - Dtot)*dq[k][j][i]*dx[2][3][k]*G3(2,i)*G4(2,j))*mdot[3][k][j][i]*DS13(k,j,i);
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      ftemp=mdot[3][k][j][i]*DS13(k,j,i);
      ftemp2=dx[2][3][k]*G3(2,i)*G4(2,j);
      dqminus=(u[km1][j][i] + (1.0 - DTOTLIMIT)*dq[km1][j][i]*ftemp2)*ftemp;
      dqplus=(u[k][j][i] + (-1.0 + DTOTLIMIT)*dq[k][j][i]*ftemp2)*ftemp;
      fl[3][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
#endif
  }
      
  LOOPC{
    v[2][1][k][j][i] += (fl[3][k][j][i]-fl[3][kp1][j][i])*OVOL1(k,j,i)*ODX(1,3,k); // no cancelling here unlike in sweepx and sweepy	
  }
  bound(v[2][1],NULL,-3,0,0) ;
}




if(transpassivev3x3&&( (passivepos==VDIR3)||(passivepos==-1))){

  if(passivepos==VDIR3)  wcom=1;
  else wcom=3;



  LOOPFP13 LOOPF2 LOOPF1{
    u[k][j][i] = v[2][wcom][k][j][i]/z2e_3(s[1],k,j,i);
  }


    
  /* vz */
  dqvz_calc(-3,workv,dqv) ;
    
  LOOPT2k{
    vpp = e2z_3(vp[3],k,j,i);
    mdotp = e2z_3(mdot[3],k,j,i);
#if(SYMFORCEHD==0)
    if(vpp > 0) 
      fl[3][k][j][i] = (u[k][j][i] + (dx[1][3][k]*G3(2,i)*G4(2,j) - vpp)*dqv[3][k][j][i])*G3(2,i)*G4(2,j)*mdotp*DS13(k,j,i);
    else
      fl[3][k][j][i] = (u[kp1][j][i] + (-dx[1][3][k]*G3(2,i)*G4(2,j) - vpp)*dqv[3][kp1][j][i])*G3(2,i)*G4(2,j)*mdotp*DS13(k,j,i);
#else
    Dtot = vpp*OARC13(k,j,i)*ODX(1,3,k);
    if(Dtot > DTOTLIMIT) 
      fl[3][k][j][i] = (u[k][j][i] + (1.0 - Dtot)*dqv[3][k][j][i]*dx[2][3][k]*G3(2,i)*G4(2,j))*G3(2,i)*G4(2,j)*mdotp*DS13(k,j,i);
    else if(Dtot< -DTOTLIMIT)
      fl[3][k][j][i] = (u[kp1][j][i] + (-1.0 - Dtot)*dqv[3][kp1][j][i]*dx[2][3][k]*G3(2,i)*G4(2,j))*G3(2,i)*G4(2,j)*mdotp*DS13(k,j,i);
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      ftemp=G3(2,i)*G4(2,j)*mdotp*DS13(k,j,i);
      ftemp2=dx[1][3][k]*G3(2,i)*G4(2,j);
      dqminus=(u[k][j][i] + (1.0 - DTOTLIMIT)*dqv[3][k][j][i]*ftemp2)*ftemp;
      dqplus=(u[kp1][j][i] + (-1.0 + DTOTLIMIT)*dqv[3][kp1][j][i]*ftemp2)*ftemp;
      fl[3][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
#endif
  }

  LOOPV3{
    v[2][wcom][k][j][i] += (fl[3][km1][j][i]-fl[3][k][j][i])*OVOL1(k,j,i)*ODX(2,3,k);
  }
  bound(NULL,v[2],0,-4,wcom) ; // bound passive scalar
}





if(transpassivev1x3&&( (passivepos==VDIR1)||(passivepos==-1))){

  if(passivepos==VDIR1)  wcom=1;
  else wcom=1;


  LOOPF3 LOOPF2 LOOPFP11{
    u[k][j][i] = v[2][wcom][k][j][i]/z2e_1(s[1],k,j,i);
  }


  /* vx */
  dqvz_calc(-1,workv,dqv) ;
    
  LOOPT3k{
    // below avgs reason for loopt1j over other bzones
    mdotp = z2e_1(mdot[3],k,j,i);
    vpp = z2e_1(vp[3],k,j,i);
#if(SYMFORCEHD==0)
    if(vpp > 0.) 
      fl[3][k][j][i] = (u[km1][j][i] + (dx[2][3][k]*G3(1,i)*G4(2,j) - vpp)*dqv[1][km1][j][i])*mdotp*DS23(k,j,i);
    else	
      fl[3][k][j][i] = (u[k][j][i] + (-dx[2][3][k]*G3(1,i)*G4(2,j) - vpp)*dqv[1][k][j][i])*mdotp*DS23(k,j,i);
#else
    Dtot=vpp*OARC13(k,j,i)*ODX(2,3,k);
    if(Dtot > DTOTLIMIT) 
      fl[3][k][j][i] = (u[km1][j][i] + (1.0 - Dtot)*dqv[1][km1][j][i]*dx[2][3][k]*G3(1,i)*G4(2,j))*mdotp*DS23(k,j,i);
    else if(Dtot < -DTOTLIMIT)
      fl[3][k][j][i] = (u[k][j][i] + (-1.0 - Dtot)*dqv[1][k][j][i]*dx[2][3][k]*G3(1,i)*G4(2,j))*mdotp*DS23(k,j,i);
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      ftemp=mdotp*DS23(k,j,i);
      ftemp2=dx[2][3][k]*G3(1,i)*G4(2,j);
      dqminus=(u[km1][j][i] + (1.0 - DTOTLIMIT)*dqv[1][km1][j][i]*ftemp2)*ftemp;
      dqplus=(u[k][j][i] + (-1.0 + DTOTLIMIT)*dqv[1][k][j][i]*ftemp2)*ftemp;
      fl[3][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
#endif
  }
    
  LOOPV1{
    v[2][wcom][k][j][i] += (fl[3][k][j][i]-fl[3][kp1][j][i])*OVOL2(k,j,i)*ODX(1,3,k);
  }
  bound(NULL,v[2],0,-4,wcom) ; // bound passive scalar
}



if(transpassivev2x3&&( (passivepos==VDIR2)||(passivepos==-1))){

  if(passivepos==VDIR2)  wcom=1;
  else wcom=2;


  LOOPF3 LOOPFP12 LOOPF1{
    u[k][j][i] = v[2][wcom][k][j][i]/z2e_2(s[1],k,j,i);
  }
    
  /* vy */
  dqvz_calc(-2,workv,dqv) ;

    
  LOOPT4k{
    vpp = z2e_2(vp[3],k,j,i);
    mdotp = z2e_2(mdot[3],k,j,i);
#if(SYMFORCEHD==0)
    if(vpp > 0) 
      fl[3][k][j][i] = (u[km1][j][i] + (dx[2][3][k]*G3(2,i)*G4(1,j) - vpp)*dqv[2][km1][j][i])*G2(2,i)*mdotp*DS33(k,j,i);
    else
      fl[3][k][j][i] = (u[k][j][i] + (-dx[2][3][k]*G3(2,i)*G4(1,j) - vpp)*dqv[2][k][j][i])*G2(2,i)*mdotp*DS33(k,j,i);
#else
    Dtot=vpp*OARC33(k,j,i);
    if(Dtot > DTOTLIMIT)
      fl[3][k][j][i] = (u[km1][j][i] + (1.0 - Dtot)*dqv[2][km1][j][i]*dx[2][3][k]*G3(2,i)*G4(1,j))*G2(2,i)*mdotp*DS33(k,j,i);
    else if(Dtot< -DTOTLIMIT)
      fl[3][k][j][i] = (u[k][j][i] + (-1.0 - Dtot)*dqv[2][k][j][i]*dx[2][3][k]*G3(2,i)*G4(1,j))*G2(2,i)*mdotp*DS33(k,j,i);
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      ftemp=G2(2,i)*mdotp*DS33(k,j,i);
      ftemp2=dx[2][3][k]*G3(2,i)*G4(1,j);
      dqminus=(u[km1][j][i] + (1.0 - DTOTLIMIT)*dqv[2][km1][j][i]*ftemp2)*ftemp;
      dqplus=(u[k][j][i] + (-1.0 + DTOTLIMIT)*dqv[2][k][j][i]*ftemp2)*ftemp;
      fl[3][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
#endif
  }


  LOOPV2{
    v[2][wcom][k][j][i] += (fl[3][k][j][i]-fl[3][kp1][j][i])*OVOL3(k,j,i)*ODX(1,3,k);
  }
  bound(NULL,v[2],0,-4,wcom) ; // bound passive scalar
}


