// emulate functions as macros for speed, they aren't inlined for some stupid reason
// see numerics.c for details
#if(TVDLF==0)

#if( (POSTPROC==1)||(LINEARINTERP==2) ) // use best for pp
#define z2e_1(f,k,j,i) (0.5*(dx[1][1][i]*f[k][j][im1mac(i)]+dx[1][1][im1mac(i)]*f[k][j][i])/dx[2][1][i])
#define z2e_2(f,k,j,i) (0.5*(dx[1][2][j]*f[k][jm1mac(j)][i]+dx[1][2][jm1mac(j)]*f[k][j][i])/dx[2][2][j])
#define e2z_1(f,k,j,i) (0.5*(f[k][j][i] + f[k][j][ip1mac(i)]))
#define e2z_2(f,k,j,i) (0.5*(f[k][j][i] + f[k][jp1mac(j)][i]))
#define v2tov1(f,k,j,i) (0.25* ((f[k][j][i] + f[k][jp1mac(j)][i])*dx[1][1][im1mac(i)] + (f[k][j][im1mac(i)] + f[k][jp1mac(j)][im1mac(i)])* dx[1][1][i])/(dx[2][1][i]))
#define v1tov2(f,k,j,i) (0.25* ((f[k][j][i] + f[k][j][ip1mac(i)])* dx[1][2][jm1mac(j)] + (f[k][jm1mac(j)][i] + f[k][jm1mac(j)][ip1mac(i)])* dx[1][2][j])/dx[2][2][j])
#define z2c_3(f,k,j,i) (0.25*((f[k][j][i]*dx[1][1][im1mac(i)] + f[k][j][im1mac(i)]*dx[1][1][i])*dx[1][2][jm1mac(j)] + (f[k][jm1mac(j)][i]*dx[1][1][im1mac(i)] + f[k][jm1mac(j)][im1mac(i)]*dx[1][1][i])*dx[1][2][j])/(dx[2][1][i]*dx[2][2][j]))
#define c2z_3(f,k,j,i) (0.25*(f[k][j][i] + f[k][j][ip1mac(i)] + f[k][jp1mac(j)][i] + f[k][jp1mac(j)][ip1mac(i)]))

#elif(LINEARINTERP==1)

#define z2e_1(f,k,j,i) (0.5*(f[k][j][im1mac(i)]+f[k][j][i]))
#define z2e_2(f,k,j,i) (0.5*(f[k][jm1mac(j)][i]+f[k][j][i]))
#define e2z_1(f,k,j,i) (0.5*(f[k][j][i] + f[k][j][ip1mac(i)]))
#define e2z_2(f,k,j,i) (0.5*(f[k][j][i] + f[k][jp1mac(j)][i]))
#define v2tov1(f,k,j,i) (0.25* (f[k][j][i] + f[k][jp1mac(j)][i] + f[k][j][im1mac(i)] + f[k][jp1mac(j)][im1mac(i)] ))
#define v1tov2(f,k,j,i) (0.25* (f[k][j][i] + f[k][j][ip1mac(i)] + f[k][jm1mac(j)][i] + f[k][jm1mac(j)][ip1mac(i)] ))
#define z2c_3(f,k,j,i) (0.25*(f[k][j][i] + f[k][j][im1mac(i)] + f[k][jm1mac(j)][i]+ f[k][jm1mac(j)][im1mac(i)] ))
#define c2z_3(f,k,j,i) (0.25*(f[k][j][i] + f[k][j][ip1mac(i)] + f[k][jp1mac(j)][i] + f[k][jp1mac(j)][ip1mac(i)]))

#elif(LINEARINTERP==0)

#define z2e_1(f,k,j,i) (f[k][j][i])
#define z2e_2(f,k,j,i) (f[k][j][i])
#define e2z_1(f,k,j,i) (f[k][j][i])
#define e2z_2(f,k,j,i) (f[k][j][i])
#define v2tov1(f,k,j,i) (f[k][j][i])
#define v1tov2(f,k,j,i) (f[k][j][i])
#define z2c_3(f,k,j,i) (f[k][j][i])
#define c2z_3(f,k,j,i) (f[k][j][i])

#endif

#elif(TVDLF==1)

#define z2e_1(f,k,j,i) (f[k][j][i])
#define z2e_2(f,k,j,i) (f[k][j][i]) 
#define e2z_1(f,k,j,i)  (f[k][j][i]) 
#define e2z_2(f,k,j,i)  (f[k][j][i]) 
#define v2tov1(f,k,j,i)  (f[k][j][i]) 
#define v1tov2(f,k,j,i)  (f[k][j][i]) 
#define z2c_3(f,k,j,i)  (f[k][j][i]) 
#define c2z_3(f,k,j,i)  (f[k][j][i]) 

#endif


// 3D -> 2D conversions
#define z2e_3(f,k,j,i) (f[k][j][i])
#define e2z_3(f,k,j,i) (f[k][j][i])
#define v2to050(f,k,j,i) v2tov1(f,k,j,i)
#define v3to005(f,k,j,i) z2c_3(f,k,j,i)
#define v1to500(f,k,j,i) v1tov2(f,k,j,i)
#define v3tov1(f,k,j,i) z2e_1(f,k,j,i)
#define v3tov2(f,k,j,i) z2e_2(f,k,j,i)
#define v1tov3(f,k,j,i) e2z_1(f,k,j,i)
#define v2tov3(f,k,j,i) e2z_2(f,k,j,i)
#define c2z_3d(f,k,j,i) c2z_3(f,k,j,i)
#define z2c_3d(f,k,j,i) z2c_3(f,k,j,i)
#define z2c_2(f,k,j,i) z2e_1(f,k,j,i)
#define c2z_2(f,k,j,i) e2z_1(f,k,j,i)
#define z2c_1(f,k,j,i) z2e_2(f,k,j,i)
#define c2z_1(f,k,j,i) e2z_2(f,k,j,i)

#if((COORD==3)||(COORD==2))
// 2d only
#define deldotv(name,wvec,k,j,i) \
( (G2(1,ip1mac(i))*G3(1,ip1mac(i))*name[wvec][1][k][j][ip1mac(i)]-G2(1,i)*G3(1,i)*name[wvec][1][k][j][i])/DVL(1,1,i)+\
(G4(1,jp1mac(j))*name[wvec][2][k][jp1mac(j)][i]-G4(1,j)*name[wvec][2][k][j][i])/(G2(2,i)*DVL(1,2,j)) )


// 2d only: general coord
#define deldotv2(name,wvec,k,j,i) \
( (G2(1,ip1mac(i))*G3(1,ip1mac(i))*name[wvec][1][k][j][ip1mac(i)]-G2(1,i)*G3(1,i)*name[wvec][1][k][j][i])/(dx[1][1][i]*G2(2,i)*G3(2,i))+\
(G4(1,jp1mac(j))*name[wvec][2][k][jp1mac(j)][i]-G4(1,j)*name[wvec][2][k][j][i])/(G2(2,i)*G4(2,j)*dx[1][2][j]) )

// 2d spc only
#define deldotv2spc(name,wvec,k,j,i) \
( (G2(1,ip1mac(i))*G3(1,ip1mac(i))*name[wvec][1][k][j][ip1mac(i)]-G2(1,i)*G3(1,i)*name[wvec][1][k][j][i])/(dx[1][1][i]*x[2][1][i]*x[2][1][i])+\
(G4(1,jp1mac(j))*name[wvec][2][k][jp1mac(j)][i]-G4(1,j)*name[wvec][2][k][j][i])/(G2(2,i)*DVL(1,2,j)) )



#define gradv11(name,wvec,k,j,i) \
( (name[wvec][1][k][j][ip1mac(i)]-name[wvec][1][k][j][i])/dx[1][1][i] )

#define gradv22(name,wvec,k,j,i) \
( (name[wvec][2][k][jp1mac(j)][i]-name[wvec][2][k][j][i])/(G2(2,i)*dx[1][2][j])+\
0.5*(name[wvec][1][k][j][i]+name[wvec][1][k][j][ip1mac(i)])/(G2(2,i))*DG2(2,i) )

#define gradv33(name,wvec,k,j,i) \
(0.5*( (name[wvec][1][k][j][i]+name[wvec][1][k][j][ip1mac(i)])/(G3(2,i))*DG3(2,i)+\
(name[wvec][2][k][j][i]+name[wvec][2][k][jp1mac(j)][i])/(G2(2,i)*G4(2,j))*DG4(2,j) ) )

// below assumes d/dphi->0 and normal vector positions for input and output
// assumes vector exists in local and 3 lower-corner&left&right zones in general
#define curlv1(name,k,j,i) \
(1.0/(G2(1,i)*G4(2,j))  *(  ( G4(1,jp1mac(j))*v3to005(name[3],k,jp1mac(j),i)-G4(1,j)*v3to005(name[3],k,j,i) )/dx[1][2][j] -0))

#define curlv2(name,k,j,i) \
(1.0/(G3(2,i)) *(0- ( G3(1,ip1mac(i))*v3to005(name[3],k,j,ip1mac(i))-G3(1,i)*v3to005(name[3],k,j,i) )/dx[1][1][i] ))

#define curlv3(name,k,j,i) \
(1.0/(G2(2,i)) * ( (G2(1,ip1mac(i))*v2to050(name[2],k,j,ip1mac(i))-G2(1,i)*v2to050(name[2],k,j,i))/dx[1][1][i] - (v1to500(name[1],k,jp1mac(j),i)-v1to500(name[1],k,j,i))/dx[1][2][j] ) )


// curl forced to result in centered quantity

#define curlcv11(name,k,j,i) (1.0/(G4(2,j)*G2(2,i))*( G4(1,jp1mac(j))*v3tov2(name[3],k,jp1mac(j),i)-G4(1,j)*v3tov2(name[3],k,j,i) )*ODX(1,2,j))

#define curlcv12(name,k,j,i) (-1.0/(G4(2,j)*G3(2,i))*(v2tov3(name[2],kp1mac(k),j,i)-v2tov3(name[2],k,j,i))*ODX(1,3,k))

#define curlcv21(name,k,j,i) (1.0/(G3(2,i)*G4(2,j))*(v1tov3(name[1],kp1mac(k),j,i)-v1tov3(name[1],k,j,i))*ODX(1,3,k))

#define curlcv22(name,k,j,i) (-1.0/(G3(2,i))*(G3(1,ip1mac(i))*v3tov1(name[3],k,j,ip1mac(i))-G3(1,i)*v3tov1(name[3],k,j,i) )*ODX(1,1,i))

#define curlcv31(name,k,j,i) (1.0/(G2(2,i))*(G2(1,ip1mac(i))*v2tov1(name[2],k,j,ip1mac(i))-G2(1,i)*v2tov1(name[2],k,j,i))*ODX(1,1,i))

#define curlcv32(name,k,j,i) (-1.0/(G2(2,i))*(v1tov2(name[1],k,jp1mac(j),i)-v1tov2(name[1],k,j,i))*ODX(1,2,j))

// natural curls, such that a first curl of vector will have curl(v) in curlv1: (0.5,0) curlv2: (0.5,0) curlv3: (0,0)
// this is such that a curl of a curl gives back natural vector positions.
// d/dphi->0

// problem if reflect?x?=1 for that ? boundary where term should go to 0 but gets 0/0, so skip that term and just set to 0
// problem with the below commented out stuff since uses volume terms, won't preserve divB=0

//(1.0/(G2(2,i)*G4(1,j))  *(  ( G4(2,j)*name[3][k][j][i]-G4(2,jm1mac(j))*name[3][k][jm1mac(j)][i] )/dx[2][2][j] ))

//#define curlvfornat1(name,k,j,i) \
//(1.0/(G2(2,i))  *(  ( G4(2,j)*name[3][k][j][i]-G4(2,jm1mac(j))*name[3][k][jm1mac(j)][i] )/DVL(2,2,j) ))

//#define curlvfornat2(name,k,j,i) \
//(1.0/(G3(1,i)) *( 0-( G3(2,i)*name[3][k][j][i]-G3(2,im1mac(i))*name[3][k][j][im1mac(i)] )/dx[2][1][i] ))

//#define curlvfornat3(name,k,j,i) \
//(1.0/(G2(1,i)) * ( (G2(2,i)*name[2][k][j][i]-G2(2,im1mac(i))*name[2][k][j][im1mac(i)])/dx[2][1][i] - (name[1][k][j][i]-name[1][k][jm1mac(j)][i])/dx[2][2][j] ) )

// the curl to be used on a previously curled thing(like emf), ends up at normal vector positions
//(1.0/(G2(1,i)*G4(2,j))  *(  ( G4(1,jp1mac(j))*name[3][k][jp1mac(j)][i]-G4(1,j)*name[3][k][j][i] )/dx[1][2][j] ))
//#define curlvbacknat1(name,k,j,i) \
//(1.0/(G2(1,i))  *(  ( G4(1,jp1mac(j))*name[3][k][jp1mac(j)][i]-G4(1,j)*name[3][k][j][i] )/DVL(1,2,j) ))

//#define curlvbacknat2(name,k,j,i) \
//(1.0/(G3(2,i)) *( 0-( G3(1,ip1mac(i))*name[3][k][j][ip1mac(i)]-G3(1,i)*name[3][k][j][i] )/dx[1][1][i] ))

//#define curlvbacknat3(name,k,j,i) \
//(1.0/(G2(2,i)) * ( (G2(1,ip1mac(i))*name[2][k][j][ip1mac(i)]-G2(1,i)*name[2][k][j][i])/dx[1][1][i] - (name[1][k][jp1mac(j)][i]-name[1][k][j][i])/dx[1][2][j] ) )

// corrected for reflective singularity at polar axis
#define curlvfornat1(name,k,j,i) \
( (G4(1,j)==0.0) ? 0.0  : (1.0/(G2(2,i)*G4(1,j))  *(  ( G4(2,j)*name[3][k][j][i]-G4(2,jm1mac(j))*name[3][k][jm1mac(j)][i] )/dx[2][2][j] )) )
 
#define curlvfornat2(name,k,j,i) \
(1.0/(G3(1,i)) *( 0-( G3(2,i)*name[3][k][j][i]-G3(2,im1mac(i))*name[3][k][j][im1mac(i)])/dx[2][1][i] ))
 
#define curlvfornat3(name,k,j,i) \
(1.0/(G2(1,i)) * ( (G2(2,i)*name[2][k][j][i]-G2(2,im1mac(i))*name[2][k][j][im1mac(i)])/dx[2][1][i] - (name[1][k][j][i]-name[1][k][jm1mac(j)][i])/dx[2][2][j] ) )
 
// the curl to be used on a previously curled thing(like emf), ends up at normal vector positions
#define curlvbacknat1(name,k,j,i) \
(1.0/(G2(1,i)*G4(2,j))  *(  ( G4(1,jp1mac(j))*name[3][k][jp1mac(j)][i]-G4(1,j)*name[3][k][j][i] )/dx[1][2][j] ))
 
#define curlvbacknat2(name,k,j,i) \
(1.0/(G3(2,i)) *( 0-( G3(1,ip1mac(i))*name[3][k][j][ip1mac(i)]-G3(1,i)*name[3][k][j][i])/dx[1][1][i] ))
 
#define curlvbacknat3(name,k,j,i) \
(1.0/(G2(2,i)) * ( (G2(1,ip1mac(i))*name[2][k][j][ip1mac(i)]-G2(1,i)*name[2][k][j][i])/dx[1][1][i] - (name[1][k][jp1mac(j)][i]-name[1][k][j][i])/dx[1][2][j] ) )
              
#elif(COORD==1)


// 2d only
#define deldotv(name,wvec,k,j,i) \
( (name[wvec][1][k][j][ip1mac(i)]-name[wvec][1][k][j][i])/DVL(1,1,i)+\
(name[wvec][2][k][jp1mac(j)][i]-name[wvec][2][k][j][i])/(DVL(1,2,j)) )


// 2d only: general coord
#define deldotv2(name,wvec,k,j,i) \
( (name[wvec][1][k][j][ip1mac(i)]-name[wvec][1][k][j][i])/(dx[1][1][i])+\
(name[wvec][2][k][jp1mac(j)][i]-name[wvec][2][k][j][i])/(dx[1][2][j]) )

// 2d spc only
#define deldotv2spc(name,wvec,k,j,i) \
( (name[wvec][1][k][j][ip1mac(i)]-name[wvec][1][k][j][i])/(dx[1][1][i])+\
(name[wvec][2][k][jp1mac(j)][i]-name[wvec][2][k][j][i])/(DVL(1,2,j)) )



#define gradv11(name,wvec,k,j,i) \
( (name[wvec][1][k][j][ip1mac(i)]-name[wvec][1][k][j][i])/dx[1][1][i] )

#define gradv22(name,wvec,k,j,i) \
( (name[wvec][2][k][jp1mac(j)][i]-name[wvec][2][k][j][i])/(dx[1][2][j]))

#define gradv33(name,wvec,k,j,i) (0)

// below assumes d/dphi->0 and normal vector positions for input and output
// assumes vector exists in local and 3 lower-corner&left&right zones in general
#define curlv1(name,k,j,i) \
((  ( z2c_3(name[3],k,jp1mac(j),i)-z2c_3(name[3],k,j,i) )/dx[1][2][j] -0))

#define curlv2(name,k,j,i) \
((0- ( z2c_3(name[3],k,j,ip1mac(i))-z2c_3(name[3],k,j,i) )/dx[1][1][i] ))

#define curlv3(name,k,j,i) \
( ( (z2c_2(name[2],k,j,ip1mac(i))-z2c_2(name[2],k,j,i))/dx[1][1][i] - (z2c_1(name[1],k,jp1mac(j),i)-z2c_1(name[1],k,j,i))/dx[1][2][j] ) )


// natural curls, such that a first curl of vector will have curl(v) in curlv1: (0.5,0) curlv2: (0.5,0) curlv3: (0,0)
// this is such that a curl of a curl gives back natural vector positions.
// d/dphi->0

// problem if reflect?x?=1 for that ? boundary where term should go to 0 but gets 0/0, so skip that term and just set to 0

//(1.0/(G2(2,i)*G4(1,j))  *(  ( G4(2,j)*name[3][k][j][i]-G4(2,jm1mac(j))*name[3][k][jm1mac(j)][i] )/dx[2][2][j] ))

#define curlvfornat1(name,k,j,i) \
((  ( name[3][k][j][i]-name[3][k][jm1mac(j)][i] )/DVL(2,2,j) ))

#define curlvfornat2(name,k,j,i) \
(( 0-( name[3][k][j][i]-name[3][k][j][im1mac(i)] )/dx[2][1][i] ))

#define curlvfornat3(name,k,j,i) \
( ( (name[2][k][j][i]-name[2][k][j][im1mac(i)])/dx[2][1][i] - (name[1][k][j][i]-name[1][k][jm1mac(j)][i])/dx[2][2][j] ) )

// the curl to be used on a previously curled thing(like emf), ends up at normal vector positions
//(1.0/(G2(1,i)*G4(2,j))  *(  ( G4(1,jp1mac(j))*name[3][k][jp1mac(j)][i]-G4(1,j)*name[3][k][j][i] )/dx[1][2][j] ))
#define curlvbacknat1(name,k,j,i) \
((  ( name[3][k][jp1mac(j)][i]-name[3][k][j][i] )/DVL(1,2,j) ))

#define curlvbacknat2(name,k,j,i) \
(( 0-( name[3][k][j][ip1mac(i)]-name[3][k][j][i] )/dx[1][1][i] ))

#define curlvbacknat3(name,k,j,i) \
(( (name[2][k][j][ip1mac(i)]-name[2][k][j][i])/dx[1][1][i] - (name[1][k][jp1mac(j)][i]-name[1][k][j][i])/dx[1][2][j] ) )


#endif


// divB=0 constraint solutions
// btheta for bcdir=1 bcdim=2
// corrected for reflective singularity at polar axis
#if(COORD>1)
#define b2m(name,k,j,i) \
( (G4(1,j)==0.0) ? ( 0.0 ) : ( (name[2][k][jp1mac(j)][i]*dx[1][1][i]*G4(1,jp1mac(j))*G3(2,i) + dx[1][2][j]*(-name[1][k][j][i]*G2(1,i)*G3(1,i) + name[1][k][j][ip1mac(i)]*G2(1,ip1mac(i))*G3(1,ip1mac(i)))*G4(2,j) ) / ( dx[1][1][i]*G4(1,j)*G3(2,i)) ) )

// Btm = (Btp dx11[i] g14[jp1mac(j)] g23[i] + dx12[j] (-Brm g12[i] g13[i] + Brp g12[ip1mac(i)] g13[ip1mac(i)]) g24[j])/(dx11[i] g14[j] g23[i]] g23[i])

// btheta for bcdir=-1 bcdim=2 (corrected so really at k,j,i)
// corrected for reflective singularity at polar axis
#define b2p(name,k,j,i) \
( (G4(1,j)==0.0) ? (0.0) : ( ( name[2][k][jm1mac(j)][i]*dx[1][1][i]*G4(1,jm1mac(j))*G3(2,i) + dx[1][2][jm1mac(j)]*(name[1][k][jm1mac(j)][i]*G2(1,i)*G3(1,i) - name[1][k][jm1mac(j)][ip1mac(i)]*G2(1,ip1mac(i))*G3(1,ip1mac(i)))*G4(2,jm1mac(j)) ) / ( dx[1][1][i]*G4(1,j)*G3(2,i) ) ) )

//Btp = (Btm dx11[i] g14[j] g23[i] + dx12[j] (Brm g12[i] g13[i] - Brp g12[ip1mac(i)] g13[ip1mac(i)]) g24[j])/(dx11[i] g14[jp1mac(j)] g23[i]])

// br for bcdir=1 bcdim=1
#define b1m(name,k,j,i) ( (dx[1][1][i]*(-name[2][k][j][i]*G4(1,j) + name[2][k][jp1mac(j)][i]*G4(1,jp1mac(j)))*G3(2,i) + name[1][k][j][ip1mac(i)]*dx[1][2][j]*G2(1,ip1mac(i))*G3(1,ip1mac(i))*G4(2,j))/( dx[1][2][j]*G2(1,i)*G3(1,i)*G4(2,j)) )

// Brm = (dx11[i] (-Btm g14[j] + Btp g14[jp1mac(j)]) g23[i] + Brp dx12[j] g12[ip1mac(i)] g13[ip1mac(i)] g24[j])/(dx12[j] g12[i] g13[i] g24[j])

// br for bcdir=-1 bcdim=1 (fixed so b1p is really at k,j,i)
#define b1p(name,k,j,i) ( (dx[1][1][im1mac(i)]*(name[2][k][j][im1mac(i)]*G4(1,j) - name[2][k][jp1mac(j)][im1mac(i)]*G4(1,jp1mac(j)))*G3(2,im1mac(i)) + name[1][k][j][im1mac(i)]*dx[1][2][j]*G2(1,im1mac(i))*G3(1,im1mac(i))*G4(2,j))/(dx[1][2][j]*G2(1,i)*G3(1,i)*G4(2,j)) )

// Brp=(dx11[i] (Btm g14[j] - Btp g14[jp1mac(j)]) g23[i] + Brm dx12[j] g12[i] g13[i] g24[j])/(dx12[j] g12[ip1mac(i)] g13[ip1mac(i)] g24[j])

#elif(COORD==1) // just a performance optimization

#define b2m(name,k,j,i) ( (name[2][k][jp1mac(j)][i]*dx[1][1][i]+ dx[1][2][j]*(-name[1][k][j][i] + name[1][k][j][ip1mac(i)]) ) / ( dx[1][1][i]) )

// btheta for bcdir=-1 bcdim=2 (corrected so really at k,j,i)
#define b2p(name,k,j,i) ( ( name[2][k][jm1mac(j)][i]*dx[1][1][i] + dx[1][2][jm1mac(j)]*(name[1][k][jm1mac(j)][i] - name[1][k][jm1mac(j)][ip1mac(i)]) ) / ( dx[1][1][i] ) )

// br for bcdir=1 bcdim=1
#define b1m(name,k,j,i) ( (dx[1][1][i]*(-name[2][k][j][i] + name[2][k][jp1mac(j)][i]) + name[1][k][j][ip1mac(i)]*dx[1][2][j])/( dx[1][2][j]) )

// br for bcdir=-1 bcdim=1 (fixed so b1p is really at k,j,i)
#define b1p(name,k,j,i) ( (dx[1][1][im1mac(i)]*(name[2][k][j][im1mac(i)] - name[2][k][jp1mac(j)][im1mac(i)]) + name[1][k][j][im1mac(i)]*dx[1][2][j])/(dx[1][2][j]) )

#endif

// dummy assignments, never used
#define b3m(name,k,j,i) (name[3][k][j][i])
#define b3p(name,k,j,i) (name[3][k][j][i])
