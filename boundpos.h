// i is boundary zone
// ii is other zone#1 -- matches with i
// i2 is 2nd boundary zone
// iii is other zone#2 -- matches with i2

//void  legacyblock1(int bct,int bcdim,int bcdir,int k,int j,int i,int*bcd1p,int*bcd2p,int*kkp,int*jjp,int*iip,int*kkkp,int*jjjp,int*iiip){
//legacyblock1(bct,bcdim,bcdir,k,j,i,&bcd1,&bcd2,&kk,&jj,&ii,&kkk,&jjj,&iii);
 
// bcd1/2 account for fact that scalar bzone "index" is offset for directional vectors(and inner/outer dir-vec diff too)
if(bcdir==1){ (bcd1)=1; (bcd2)=0;}
else if(bcdir==(-1)){ (bcd1)=0; (bcd2)=1;}
else{fprintf(fail_file,"error: ll: %d k,j,i: %d %d %d bound.c: bcdir out of bounds: %d\n\r",ll,k,j,i,bcdir); myexit(1);}

if(bct!=5){
  ii=i+bcdir1;
  jj=j+bcdir2;
  kk=k+bcdir3;
#if(NBIGBND==2)
  i2=i-bcdir1;
  j2=j-bcdir2;
  k2=k-bcdir3;
  iii=i+2*bcdir1;
  jjj=j+2*bcdir2;
  kkk=k+2*bcdir3;
#endif
}
else{
    
  /* Determine which zone to copy from  in period conditions */
  
  /* These non-local conditions are not general to any problem */
  /* The below assumes edges of rect-domain are periodic--standard */

    ii=i+bcdir1*N1;
    i2=i-bcdir1;
    iii=i+bcdir1*(N1-1);
  if((bcdim!=2)||(periodicx2special==0)){
    jj=j+bcdir2*N2;
    j2=j-bcdir2;
    jjj=j+bcdir2*(N2-1);
    kk=k+bcdir3*N3;
    k2=k-bcdir3;
    kkk=k+bcdir3*(N3-1);
  }
  else{ // then dealing with x2 boundary when periodicx2special==1
    j2=j-bcdir2; // this is first layer if outer N2 and periodicx2special and edge quantity
    jj=j+bcdir2;
    jjj=j+bcdir2*2;

    jje2=j+2*bcdir2;
    jjje2=j+bcdir2*3; // only used if j2!=N2+2 (i.e. j2==-2)

    kk=(k+N3/2)%(N3);
    k2=k;
    kkk=(k+N3/2)%(N3);

  }

}


