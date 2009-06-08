
if(transpassivex1){
  // add passive scalar to inside accretor
  // must come after "field" set
#if(DO_ACCRETOR)
  if(passivepos==CENT){
    bound_accretor(1,v[2][1],NULL,-3,0,0); // direction (first argument) doesn't matter
  }
  else{
    if(passivepos!=-1) bound_accretor(1,NULL,v[2],0,-4,1); // direction (first argument) doesn't matter
    else bound_accretor(1,NULL,v[2],0,-4,123);
  }
#endif
}

workv=(FTYPE (*) [N3M][N2M][N1M])(&u[0][0][0]);



////////////////////////////////////////////////////////
//
////////////////// PASSIVE SCALAR WHEN DOING HYDRO (use Bx as scalar)
//
///////////////////////////////////////////////////////
if(transpassivex1&&(passivepos==CENT)){
  LOOPFC{ // full due to dqx_calc requirements
    u[k][j][i] = v[2][1][k][j][i]/s[1][k][j][i] ;
  }
      
  dqx_calc(u,dq) ;

  // fl is really fl=Flux*dt/dx (because of mdot)
  LOOPT1i{ // same as density w.r.t. loop/bound
#if(SYMFORCEHD==0)
    if(vp[1][k][j][i] > 0.) 
      fl[1][k][j][i] = (u[k][j][im1] + (dx[2][1][i] - vp[1][k][j][i])*dq[k][j][im1])*mdot[1][k][j][i]*DS11(k,j,i); // dx[1][3][k] cancels with below volume
    else	
      fl[1][k][j][i] = (u[k][j][i] + (-dx[2][1][i] - vp[1][k][j][i])*dq[k][j][i])*mdot[1][k][j][i]*DS11(k,j,i);
#else
    Dtot=vp[1][k][j][i]*OARC11(k,j,i);
    if(Dtot > DTOTLIMIT )
      fl[1][k][j][i] = (u[k][j][im1] + (1.0 - Dtot)*dq[k][j][im1]*dx[2][1][i])*mdot[1][k][j][i]*DS11(k,j,i);
    else if(Dtot < -DTOTLIMIT)
      fl[1][k][j][i] = (u[k][j][i] + (-1.0 - Dtot)*dq[k][j][i]*dx[2][1][i])*mdot[1][k][j][i]*DS11(k,j,i);
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      ftemp = mdot[1][k][j][i]*DS11(k,j,i);
      dqminus=(u[k][j][im1] + (1.0 - DTOTLIMIT)*dq[k][j][im1]*dx[2][1][i])*ftemp;
      dqplus=(u[k][j][i] + (-1.0 + DTOTLIMIT)*dq[k][j][i]*dx[2][1][i])*ftemp;
      fl[1][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
#endif
  }      
  // assign internal energy density from surface fluxes
  LOOPC{
    v[2][1][k][j][i] += (fl[1][k][j][i]-fl[1][k][j][ip1])*OVOL1(k,j,i); // *odx[1][3][k]
  }
  bound(v[2][1],NULL,-3,0,0) ; // bound passive scalar
}






if(transpassivev1x1&& ((passivepos==VDIR1)||(passivepos==-1))){

  if(passivepos==VDIR1)  wcom=1;
  else wcom=1;

  LOOPF3 LOOPF2 LOOPFP11{ // full due to dqx_calc requirements
    u[k][j][i] = v[2][wcom][k][j][i]/z2e_1(s[1],k,j,i);
  }




  /* vx */
  dqvx_calc(-1,workv,dqv) ;

    
  LOOPT2i{ // only need to loop over -1 to N-1(1 bnd zone) and should not bound fl!
    vpp = e2z_1(vp[1],k,j,i);
    mdotp = e2z_1(mdot[1],k,j,i);
#if(SYMFORCEHD==0)
    if(vpp > 0)
      fl[1][k][j][i] = (u[k][j][i]   + (dx[1][1][i] - vpp)*dqv[1][k][j][i])  *mdotp*DS21(k,j,i);
    else
      fl[1][k][j][i] = (u[k][j][ip1] + (-dx[1][1][i] - vpp)*dqv[1][k][j][ip1])*mdotp*DS21(k,j,i);
#else
    Dtot=vpp*OARC21(k,j,i);
    if(Dtot > DTOTLIMIT) 
      fl[1][k][j][i] = (u[k][j][i]   + (1.0 - Dtot)*dqv[1][k][j][i]*dx[1][1][i])*mdotp*DS21(k,j,i);
    else if(Dtot< -DTOTLIMIT)
      fl[1][k][j][i] = (u[k][j][ip1] + (-1.0 - Dtot)*dqv[1][k][j][ip1]*dx[1][1][i])*mdotp*DS21(k,j,i);
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      ftemp = mdotp*DS21(k,j,i);
      dqminus=(u[k][j][i]   + (1.0 - DTOTLIMIT)*dqv[1][k][j][i]*dx[1][1][i])*ftemp;
      dqplus=(u[k][j][ip1] + (-1.0 + DTOTLIMIT)*dqv[1][k][j][ip1]*dx[1][1][i])*ftemp;
      fl[1][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
#endif    
  }
    
  LOOPV1{
    v[2][wcom][k][j][i] += (fl[1][k][j][im1]-fl[1][k][j][i])*OVOL2(k,j,i);
  }
  bound(NULL,v[2],0,-4,wcom) ; // bound passive scalar
}



if(transpassivev2x1&&( (passivepos==VDIR2)||(passivepos==-1))){

  if(passivepos==VDIR2)  wcom=1;
  else wcom=2;


  LOOPF3 LOOPFP12 LOOPF1{ // full due to dqx_calc requirements
    u[k][j][i] = v[2][wcom][k][j][i]/z2e_2(s[1],k,j,i);
  }

  /* vy */
  dqvx_calc(-2,workv,dqv) ;
    
  LOOPT3i{ // flux here lives on lower left corner
    mdotp = z2e_2(mdot[1],k,j,i);
    vpp = z2e_2(vp[1],k,j,i);
#if(SYMFORCEHD==0)     
    if(vpp > 0.)
      fl[1][k][j][i] = (u[k][j][im1] + (dx[2][1][i] - vpp)*dqv[2][k][j][im1])*G2(1,i)*mdotp*DS31(k,j,i);
    else	
      fl[1][k][j][i] = (u[k][j][i]   + (-dx[2][1][i] - vpp)*dqv[2][k][j][i])  *G2(1,i)*mdotp*DS31(k,j,i);
#else
    Dtot=vpp*OARC11(k,j,i);
    if(Dtot > DTOTLIMIT)
      fl[1][k][j][i] = (u[k][j][im1] + (1.0 - Dtot)*dqv[2][k][j][im1]*dx[2][1][i])*G2(1,i)*mdotp*DS31(k,j,i);
    else if(Dtot< -DTOTLIMIT)
      fl[1][k][j][i] = (u[k][j][i]   + (-1.0 - Dtot)*dqv[2][k][j][i]*dx[2][1][i])*G2(1,i)*mdotp*DS31(k,j,i);
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      ftemp = G2(1,i)*mdotp*DS31(k,j,i);
      dqminus= (u[k][j][im1] + (1.0 - DTOTLIMIT)*dqv[2][k][j][im1]*dx[2][1][i])*ftemp;
      dqplus=(u[k][j][i]   + (-1.0 + DTOTLIMIT)*dqv[2][k][j][i]*dx[2][1][i])*ftemp;
      fl[1][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
#endif
  }
    
  LOOPV2{
    v[2][wcom][k][j][i] += (fl[1][k][j][i]-fl[1][k][j][ip1])*OVOL3(k,j,i);
  }
  bound(NULL,v[2],0,-4,wcom) ; // bound passive scalar
}




if(transpassivev3x1&&((passivepos==VDIR3)||(passivepos==-1))){

  if(passivepos==VDIR3)  wcom=1;
  else wcom=3;


  LOOPFP13 LOOPF2 LOOPF1{ // full due to dqx_calc requirements
    u[k][j][i] = v[2][wcom][k][j][i]/z2e_3(s[1],k,j,i);
  }

    
  /* vz */
  dqvx_calc(-3,workv,dqv) ;
    
  LOOPT4i{ //  dx[?][3][k] again cancels
    mdotp = z2e_3(mdot[1],k,j,i);
    vpp = z2e_3(vp[1],k,j,i);
#if(SYMFORCEHD==0)
    if(vpp > 0) 
      fl[1][k][j][i] = (u[k][j][im1] + (dx[2][1][i] - vpp)*dqv[3][k][j][im1])*G3(1,i)*G4(2,j)*mdotp*DS11(k,j,i);
    else
      fl[1][k][j][i] = (u[k][j][i] + (-dx[2][1][i] - vpp)*dqv[3][k][j][i])*G3(1,i)*G4(2,j)*mdotp*DS11(k,j,i);
#else
    Dtot = vpp*OARC11(k,j,i);
    if(Dtot > DTOTLIMIT) 
      fl[1][k][j][i] = (u[k][j][im1] + (1.0 - Dtot)*dqv[3][k][j][im1]*dx[2][1][i])*G3(1,i)*G4(2,j)*mdotp*DS11(k,j,i);
    else if(Dtot< -DTOTLIMIT)
      fl[1][k][j][i] = (u[k][j][i] + (-1.0 - Dtot)*dqv[3][k][j][i]*dx[2][1][i])*G3(1,i)*G4(2,j)*mdotp*DS11(k,j,i);
    else{
      fraczone=Dtot*INVDTOTLIMIT;
      ftemp =G3(1,i)*G4(2,j)*mdotp*DS11(k,j,i);
      dqminus=(u[k][j][im1] + (1.0 - DTOTLIMIT)*dqv[3][k][j][im1]*dx[2][1][i])*ftemp;
      dqplus=(u[k][j][i] + (-1.0 + DTOTLIMIT)*dqv[3][k][j][i]*dx[2][1][i])*ftemp;
      fl[1][k][j][i]=0.5*(dqminus*(1.0+fraczone)+dqplus*(1.0-fraczone));
    }
			      
#endif
  }
    
  LOOPV3{
    v[2][wcom][k][j][i] += (fl[1][k][j][i]-fl[1][k][j][ip1])*OVOL1(k,j,i);
  }
  bound(NULL,v[2],0,-4,wcom) ; // bound passive scalar
}

