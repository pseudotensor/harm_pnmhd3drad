// choose smallest and largest

#define DIRSPLIT (1)

// ARC LENGTHS used should really change for zone centered and velocities, but close enough
#if(N1>1)
odx1=ods=odl=OARC21(k,j,i);
#else
odx1=ods=odl=0;
#endif
// ods=1/dsmallest
// odl=1/dlargest

#if(N2>1)
odx2=OARC32(k,j,i);
if( ods < odx2 ){
  ods = odx2;
}
else if( odl > odx2 ){
  odl = odx2 ;
}
#else
odx2=0;
#endif

#if(N3>1)
odx3=OARC13(k,j,i)*ODX(1,3,k);
if( ods < odx3 ){
  ods = odx3;
}
else if( odl > odx3 ){
  odl = odx3 ;
}
#endif

u = s[2][k][j][i];
rho=s[1][k][j][i];

// sound speed
idt2[1] = SSMALL; // don't use anymore
ftemp=gam*(gam-1.)*s[2][k][j][i];

#if(CSLIMIT==0)
if(wgam) cs2 = ftemp/rho;
else cs2=cs*cs;
#else
if(wgam) cs2 = ftemp/(rho+ftemp*invsol2);
else cs2=cs*cs; // assumed less than c=1
#endif


if(mag==1){
  /* alfven velocity */
  bxa = e2z_1(v[2][1],k,j,i);
  bya = e2z_2(v[2][2],k,j,i);
  bza = e2z_3(v[2][3],k,j,i);

  ftemp=(bxa*bxa + bya*bya + bza*bza); // b^2
#if(ALFVENLIMIT==0)
  valphen2=ftemp/rho; // va^2=B^2/rho
#else
  valphen2=ftemp/(rho+ftemp*invsol2) ; // va=b^2/(rho+b^2/c^2)
#endif


}
else{
  valphen2=SSMALL;  
}

velfastm=sqrt(valphen2+cs2);

#if(DIRSPLIT==1)
// x-velocity
#if(N1>1)
vel1 =e2z_1(v[1][1],k,j,i)-vg[1] ;
ftemp = invcour*(fabs(vel1)+velfastm)*odx1 ;
idt2[2]=ftemp*ftemp;
#else
idt2[2]=SSMALL;
#endif

#if(N2>1)
// y-velocity
vel2 =e2z_2(v[1][2],k,j,i)-vg[2] ;
ftemp = invcour*(fabs(vel2)+velfastm)*odx2 ;
idt2[3]=ftemp*ftemp;
#else
idt2[3]=SSMALL;
#endif

#if((COMPDIM==3)&&(N3>1))
// z-velocity
vel3 =e2z_3(v[1][3],k,j,i)-vg[3] ;
ftemp = invcour*(fabs(vel3)+velfastm)*odx3 ;
idt2[4]=ftemp*ftemp;
#else
idt2[4]=SSMALL;
#endif

#else // else if DIRSPLIT==0

#if(N1>1)
vel1 =e2z_1(v[1][1],k,j,i)-vg[1] ;
#else
vel1=0; // since then independent of that direction
#endif
#if(N2>1)
vel2 =e2z_2(v[1][2],k,j,i)-vg[2] ;
#else
vel2=0;
#endif
#if(N3>1)
vel3 =e2z_3(v[1][3],k,j,i)-vg[3] ;
#else
vel3=0;
#endif

ftemp = invcour*(sqrt(vel1*vel1*odx1*odx1+vel2*vel2*odx2*odx2+vel3*vel3*odx3*odx3)+velfastm*ods) ;
idt2[2]=ftemp*ftemp;
idt2[3]=SSMALL;
idt2[4]=SSMALL;

#endif// DIRSPLIT==0


if(visc_art==1){
  /* linear viscosity */
#if(VISC_LINEAR)
  ftemp = invcour2*nu_l*ods ;
  idt2[5] = ftemp*ftemp*cs2;
#else
  idt2[5]=SSMALL;
#endif

#if(VISC_TENSOR==0)
  /* VNR viscosity: x-dir */
#if(N1>1)
  dvx = v[1][1][k][j][ip1]-v[1][1][k][j][i] ;
  dvdx=dvx*odx1;
#else
  dvdx=0;
#endif

#if(N2>1)
  dvy = v[1][2][k][jp1][i]-v[1][2][k][j][i] ;
  ftemp=dvy*odx2;
  if(dvdx<ftemp) dvdx=ftemp;
#endif

#if(N3>1)
  dvz = v[1][3][kp1][j][i]-v[1][3][k][j][i] ;
  ftemp=dvz*odx3;
  if(dvdx<ftemp) dvdx=ftemp;
#endif

  ftemp=invcour2*nu_vnr*dvdx ;
  idt2[6] = ftemp*ftemp;
  idt2[7]=SSMALL;

#else  // else if tensor art visc

  delv = deldotv(1,k,j,i);
  l2_ten = nu_ten*ds*ds;
  ftemp= invcour2*l2_ten*delv*ods*ods ;
  idt2[6] = ftemp*ftemp;

  idt2[7]=SSMALL;
  
  /*
    printf("%15.10g %15.10g %15.10g %15.10g\n"
    ,(invcour2*nu_vnr*(v[1][1][k][j][ip1]-v[1][1][k][j][i]))*odx1
    ,idt[6]
    ,(invcour2*nu_vnr*(v[1][1][k][jp1][i]-v[1][1][k][j][i]))*odx2
    ,idt[7]);
  */
#endif // end if tens art visc
}
else{ // else if no art visc
  idt2[5]=SSMALL;
  idt2[6]=SSMALL;
  idt2[7]=SSMALL;
}

if(visc_real==1){
  
  ftemp=invcour2*nu_real[k][j][i]*ods*ods ; // just constrain by smaller zone
  idt2[8] = ftemp*ftemp;
}
else{ // else if no real visc
  idt2[8]=SSMALL;
}
idt2[9] = SSMALL; // no longer used

if(RESMEM&&(res_real==1)){
  /* resistivity */
  ftemp= invcour2*nu_res_real[k][j][i]*ods*ods ; // constrain by smaller zone dimension
  idt2[10] = ftemp*ftemp;
}
else{
  idt2[10] = SSMALL;
}
//printf("k: %3d j: %3d i: %3d [1]: %5.5g [2]: %5.5g [3]: %5.5g [4]: %5.5g [5]: %5.5g [6]: %5.5g [7]: %5.5g [8]: %5.5g [9]: %5.5g [10]: %5.5g\n",k,j,i,idt[1],idt[2],idt[3],idt[4],idt[5],idt[6],idt[7],idt[8],idt[9],idt[10]);


#if(!SUBCOOL)
  if(dx[2][2][j]*G2(2,i)*(kapk[k][j][i]+sigk[k][j][i])<1.&&thingrid[k][j][i]==0){
  idt2[11]=1000./(s[2][k][j][i]*dx[2][2][j]*G2(2,i)*(kapk[k][j][i]+sigk[k][j][i])/coolerKai[k][j][i])/(s[2][k][j][i]*dx[2][2][j]*G2(2,i)*(kapk[k][j][i]+sigk[k][j][i])/coolerKai[k][j][i]);
  }else{
  idt2[11]=1000./(s[2][k][j][i]/coolerKai[k][j][i])/(s[2][k][j][i]/coolerKai[k][j][i]);
  }
#endif
