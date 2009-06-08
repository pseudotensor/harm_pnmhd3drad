#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif

// A cell(zone)
//////////////////////////////////////////
//(0,1)    (0.5,1)     (1,1)   //       //
//                             //       //
//                             //       //
//(0,.5)   (0.5,0.5)   (1,0.5) //       //
//                             //       //
//                             //       //
//(0,0)    (0.5,0)     (1,0)   //       //
//////////////////////////////////////////
//---> x1-dir
//^
//|
//|
//x2-dir
// (x1,x2) --> (x1,x2)
/*
// (0.5,0.5) --> (0,.5)
// (.5,0) --> (0,0)
FTYPE z2e_1(FTYPE (*var)[N1M],int j,int i)
{
  return(0.5*(dx[1][1][i]*var[j][i-1]+dx[1][1][i-1]*var[j][i])/dx[2][1][i]);
}

// (0.5,0.5) --> (0.5,0)
// (0,0.5) --> (0,0)
FTYPE z2e_2(FTYPE (*var)[N1M],int j,int i)
{
  // metric drops out, so used for both above
  return(0.5*(dx[1][2][j]*var[j-1][i]+dx[1][2][j-1]*var[j][i])/dx[2][2][j]);
}




// (0,.5) --> (0.5,0.5)
FTYPE e2z_1(FTYPE (*var)[N1M],int j,int i)
{
  return(0.5*(var[j][i] + var[j][i + 1]));
}

// (.5,0) --> (0.5,0.5)
FTYPE e2z_2(FTYPE (*var)[N1M],int j,int i)
{
  return(0.5*(var[j][i] + var[j + 1][i]));
}

// (0.5,0) --> (0,0.5)
FTYPE e2e_v2(FTYPE (*var)[N1M],int j, int i) // v2 label means usually operated on v2
{
  // analytically depends on (z2e_1 and e2z_2) or (e2z_2 and z2e_1), either order gives same result
  return(0.25* ((var[j][i] + var[1 + j][i])*dx[1][1][-1 + i] + (var[j][-1 + i] + var[1 + j][-1 + i])* dx[1][1][i])/(dx[2][1][i]));
}

// (0,0.5) --> (.5,0)
FTYPE e2e_v1(FTYPE (*var)[N1M],int j, int i)  // v1 label means usually operated on v1
{
  // analytically depends on (z2e_2 and e2z_1) or (e2z_1 and z2e_2), either order gives same result
  return(0.25* ((var[j][i] + var[j][1 + i])* dx[1][2][-1 + j] + (var[-1 + j][i] + var[-1 + j][1 + i])* dx[1][2][j])/dx[2][2][j]);
}

// (.5,.5) -> (0,0)
FTYPE z2c(FTYPE (*var)[N1M],int j, int i)
{
  // analytically depends on (z2e_1 and z2e_2) or (z2e_2 and z2e_1), either order
  return(0.25*((var[j][i]*dx[1][1][-1 + i] + var[j][-1 + i]*dx[1][1][i])*dx[1][2][-1 + j] + (var[-1 + j][i]*dx[1][1][-1 + i] + var[-1 + j][-1 + i]*dx[1][1][i])*dx[1][2][j])/(dx[2][1][i]*dx[2][2][j]));
}
// (0,0) -> (.5,.5)
FTYPE c2z(FTYPE (*var)[N1M],int j, int i)
{
  return(0.25 (var[j][i] + var[j][1 + i] + var[1 + j][i] + var[1 + j][1 + i]))
  }
*/

// extrapolate left middle edge value(vector)
// dir: which direction to extrapolate to
// assumes other 2 points exist in normal array
// which:
// 0: linear extrap (just (0,.5) -> (0,-.5) for dim=1 and (.5,0)->(-.5,0) for dim=2)
// 1: reflect extrap(same, just make edge 0)(i,j here is center point when dir=-1 to account for staggered mesh

void ex_v(int dim,int dir,FTYPE (*var)[N1M], int k,int j, int i,int which) // not correct in nonuni grid and curv coords
{
  static FTYPE tempf;
  FTYPE temp1;
  
  temp1=var[j][i];

  if(dim==1){
    if(dir==1){
      tempf=dx[1][1][i]/dx[1][1][i + 1];
      if(which==1){
	var[j][i+1]=0.0;
	var[j][i] = - var[j][i + 2]*tempf;
      }
      else var[j][i]=var[j][i + 1]*(1.0 + tempf) - var[j][i + 2]*tempf;
    }
    else{
      if(which==1) i=i+1; // as an artifact of the staggered mesh and what boundary zone numbers define where the boundary is
      tempf=dx[1][1][i-1]/dx[1][1][i - 2];
      if(which==1){
	var[j][i-1]=0.0;
	if(i<=N1) var[j][i] = - var[j][i - 2]*tempf;
	// need if, since no outer vector memory exists and i=N1+1 is proper boundary zone outer edge
      }
      else var[j][i]=var[j][i - 1]*(1.0 + tempf) - var[j][i - 2]*tempf;
    }      
  }
  else{
    if(dir==1){
      tempf=dx[1][2][j]/dx[1][2][j + 1];
      if(which==1){
	var[j+1][i]=0.0;
	var[j][i] = - var[j + 2][i]*tempf;
      }
      else var[j][i] = var[j + 1][i]*(1.0 + tempf) - var[j + 2][i]*tempf;
    }
    else{
      if(which==1) j=j+1; // as an artifact of the staggered mesh and what boundary zone numbers define where the boundary is
      tempf=dx[1][2][j-1]/dx[1][2][j - 2];
      if(which==1){
	var[j-1][i]=0.0;
	if(j<=N2) var[j][i] = - var[j - 2][i]*tempf;
      }
      else var[j][i] = var[j - 1][i]*(1.0 + tempf) - var[j - 2][i]*tempf;
    }      
  }
  /*
  if(fabs(var[j][i])>5.0){
    printf("problem in ex_v: %15.10g %15.10g %15.10g %15.10g %d %d\n",var[j][i],var[j][i+1],var[j][i+2],temp1,j,i);
    //    var[j][i+1]=var[j][i+2];
  }
  */
}

// extrapolate zone center value(scalar)
// dir: which direction to extrapolate to
// assumes other 2 points exist in normal array

// which(make sure: 0 or 1):
// 0: linear extrap
// 1: reflecting extrap
void ex_s(int dim,int dir,FTYPE (*var)[N1M], int k,int j, int i,int which)
{
  FTYPE tempf;

  if(dim==1){
    if(dir==1){
      if(which==0) tempf=dx[2][1][i+1]/dx[2][1][i + 2];
      else if(which==1) tempf=(dx[1][1][i+1]-dx[2][1][i+1])/dx[2][1][i+2];
      var[j][i] = var[j][i + 1]*(1.0 + tempf) - var[j][i + 2]*tempf;
    }
    else{
      if(which==0) tempf=dx[2][1][i]/dx[2][1][i - 1];
      else if(which==1) tempf=(dx[1][1][i]-dx[2][1][i])/dx[2][1][i-1];
      var[j][i] = var[j][i - 1]*(1.0 + tempf) - var[j][i - 2]*tempf;
    }      
  }
  else{
    if(dir==1){
      if(which==0) tempf=dx[2][2][j+1]/dx[2][2][j + 2];
      else if(which==1) tempf=(dx[1][2][j+1]-dx[2][2][j+1])/dx[2][2][j+2];
      var[j][i] = var[j + 1][i]*(1.0 + tempf) - var[j + 2][i]*tempf;
    }
    else{
      if(which==0) tempf=dx[2][2][j]/dx[2][2][j - 1];
      else if(which==1) tempf=(dx[1][2][j]-dx[2][2][j])/dx[2][2][j-1];
      var[j][i] = var[j - 1][i]*(1.0 + tempf) - var[j - 2][i]*tempf;
    }      
  }
}

// periodic extrapolation for scalars
void ex_s_p(int dim,int dir,FTYPE (*var)[N1M], int k,int j, int i)
{
  FTYPE tempf;

  if(dim==1){
    if(dir==1) tempf=(-2.*dx[2][1][0]/(dx[1][1][0]+dx[1][1][N1-1]));
    else tempf=(dx[1][1][N1]-dx[1][1][0])/(dx[1][1][N1-1]+dx[1][1][0]);
    var[j][i] = var[j][0]*(1.0 + tempf) - var[j][N1-1]*tempf;
  }
  else{
    if(dir==1) tempf=(-2.*dx[2][2][0]/(dx[1][2][0]+dx[1][2][N2-1]));
    else tempf=(dx[1][2][N2]-dx[1][2][0])/(dx[1][2][N2-1]+dx[1][2][0]);
    var[j][i] = var[0][i]*(1.0 + tempf) - var[N2-1][i]*tempf;
  }
}


// periodic extrapolation for vectors
// no need for correction due to stagerred mesh since explicit about indicies
void ex_v_p(int dim,int dir,FTYPE (*var)[N1M], int k,int j, int i)
{
  FTYPE tempf;

  if(dim==1){
    if(dir==1){
      tempf=-dx[1][1][-1]/dx[1][1][N1-1];
      var[j][i] = var[j][0]*(1.0 + tempf) - var[j][N1-1]*tempf;
    }
    else var[j][i] = var[j][0];    
  }
  else{
    if(dir==1){
      tempf=-dx[1][2][-1]/dx[1][2][N2-1];
      var[j][i] = var[0][i]*(1.0 + tempf) - var[N2-1][i]*tempf;
    }
    else var[j][i] = var[0][i];    
  }
}
   






// assume never interp when postproc==0
#if(POSTPROC==1)

// uses wg1 for only 2 grid stores, should have 4 if using wg1/wg2 sources


// could just use 3d card or opengl for square interpolation, but point of routine is for spherical coords--do 3d cards support that texture type?  Boy that would be nice as the below is error prone with so much stuff

// wtype: 0: scalar input 1: vx-like input 2: vy-like input.  3) sig12 type input : Determines where source data is located
// stype: 1: point sampling 2: planar interpolation
// whichg: 1: put pixel centered onto a-grid centers 2: b-grid
// inn1: input size in x1-dir
// inn2: intput size in x2-dir
// outn1: output size in x1-dir
// outn2: output size in x2-dir
// in: input data
// out: output data
// outerdef: what to set value of function to if in non-mapped domain
// outtype:  type of interp(see init_dx/x)

#define DYNAMICCAM 0
#define EXTRAPOLATE 0 // whether to extrapolate or to set to outerdef value when pixel is outside grid
// extrapolation can be bad if done far outside range and linear interpolation is on

void interpolate(int wtype,int wsv,int stype, int whichg, int inn1, int inn2, int outn1, int outn2, FTYPE in[N2M][N1M], FTYPE out[INTN2][INTN1],FTYPE outerdef, int outtype)
{

  int i,j,k,l,m,ll,mm ;
  int p,q;
  int wshort;
  int wg1, wg2;
  int which,whichtemp;
  int gotx1;
  int gi,gj;
  int skip;
  FTYPE xp1,xp2;
  FTYPE x1,x2,x3,y1,y2,y3,z1,z2,z3;
  FTYPE dist[4];
  int N[3+1];
  int NB[3+1];
  FTYPE a3,b3,c3;
  FTYPE boc,aoc;
  FTYPE Cval;
  static int firstfirsttime=1;
  static int firsttime[2][2];
  int error;

  //if(!( (SAMPLED==1)&&(SAMPLEI==1)&&( (IMGN1!=DUMN1)||(IMGN2!=DUMN2))))

  // stored values to speed up calc since mostly doing nonmoving grids
  static FTYPE Cvala[2][2][INTN2][INTN1];
  static int whicha[2][2][INTN2][INTN1];
  static FTYPE Aco[2][2][INTN2][INTN1];
  static FTYPE Bco[2][2][INTN2][INTN1];
  static FTYPE a[2][2][2][INTN2][INTN1];
  static FTYPE b[2][2][2][INTN2][INTN1];

  static int la[2][INTN2][INTN1];
  static int ma[2][INTN2][INTN1];
  static int lla[2][INTN2][INTN1];
  static int mma[2][INTN2][INTN1];


  // can request refresh of interp, whatever predetermined factors accounted for here.
  if(globalinterpmod==1){
    globalinterpmod=0;
    firstfirsttime=1;
  }



  if(firstfirsttime==1){    
    firsttime[0][0]=1;
    firsttime[1][0]=1;
    firsttime[1][1]=1;
    firsttime[0][1]=1;
  }
  firsttime[0][0]=1;
  firsttime[1][0]=1;
  firsttime[1][1]=1;
  firsttime[0][1]=1;


  N[1]=outn1;
  N[2]=outn2;
  N[3]=1;
  NB[1]=0;
  NB[2]=0;
  NB[3]=0;

  // where source data located
  if(wtype==0){ wg1=2; wg2=2;} // pure scalar type grid
  else if(wtype==1){ wg1=1; wg2=2;} // vx-type grid
  else if(wtype==2){ wg1=2; wg2=1;} // vy-type grid
  else if(wtype==3){ wg1=1; wg2=1;} // 0,0 grid
  else{
    fprintf(fail_file,"error, no defined wtype: %d\n",wtype);
    myexit(1);
  }
  
  // overrides:
  // revert to non-memory method if no consistency in interpolated grid size
  // can be fixed by using more memory
  if( (SAMPLED==1)&&(SAMPLEI==1)&&( (IMGN1!=DUMN1)||(IMGN2!=DUMN2))) { firstfirsttime=1; firsttime[wg2-1][wg1-1]=1; }
  // also don't do memory method if grid changing each time or doing 2 kinds of sampling for image and data each
  if( (SAMPLEI!=SAMPLED)&&(SAMPLED>0)&&(SAMPLEI>0))  { firstfirsttime=1; firsttime[wg2-1][wg1-1]=1; }
  // don't use memory if doing moving camera
  if(DYNAMICCAM==1)  { firstfirsttime=1; firsttime[wg2-1][wg1-1]=1; }
  // end overrides
  
  //  printf("firsttime[wg2-1][wg1-1]: %d\n",firsttime[wg2-1][wg1-1]);
  
  //  firsttime[wg2-1][wg1-1]=1; // override for now to do dual imaging // no longer doing dual imaging

  
  if(0&&(firstfirsttime==1)){ // XXXX----  assume taken care of in postproc.c

    // this gets both a and b grids
    error=0;
    error+=init_dx(N[1],N[2],N[3],NB[1],NB[2],NB[3],startix,1,outtype);
    error+=init_x(N[1],N[2],N[3],NB[1],NB[2],NB[3],startix,1,outtype);
    if(error>0){
      fprintf(fail_file,"error>0 in numerics.c\n");
      myexit(1);
    }

  }//end firstfirsttime  //----XXXX
  // only really good if dr and dtheta are on the same order(i.e. grid is close to square(256x256))
  // Interpolate physical values on physics grid to physical value on pixel grid
  for(j=0;j<N[2];j++){
    for(i=0;i<N[1];i++){
      if(firsttime[wg2-1][wg1-1]==1){ // ****----
	// determine where pixel center is w.r.t. physical grid
	if((COORD==3)&&(outtype>0)){// outtype==0 means going to same grid type
	  if( (outtype==1)||(outtype==2)){
	    xp1=sqrt(ix[wg1][1][i]*ix[wg1][1][i]+ix[wg2][2][j]*ix[wg2][2][j]); // xp1 is r-coord
	    xp2=acos(ix[wg2][2][j]/(xp1+1.E-6)); // y-axis is along theta=0, so xp2=theta
	  }
	}
	else{
	  xp1=ix[wg1][1][i];
	  xp2=ix[wg2][2][j];
	}
	// loop over input data
	skip=0; // assume will find a good point

	for(l=0;l<inn1;l++){
	  if(x[wg1][1][l]>=xp1) break;
	}
	for(m=0;m<inn2;m++){
	  if(x[wg2][2][m]>=xp2) break;
	}
	//fprintf(log_file,"%d %d %d %d\n",j,i,m,l);
	
	//now have upper value of a quadrant bounding current pixel center: x[wg1/2][2/1][m/l]
        if(l==0){
          skip=1;
        }
        if(m==0){
          skip=1;
        }
	if(l==inn1){
	  skip=1;
	}
	if(m==inn2){
	  skip=1;
	}
	if(EXTRAPOLATE){
	  skip=0; // force extrapolation instead of outdef setting
	}
	//if((i==0)&&(j==N[2]/2)){
	//fprintf(stderr,"%d %d %d %d %d\n",j,i,m,l,skip);
	  //myexit(1);
	//}
	  // reach here if need pixel value and "inside" active domain or forcing extrapolation
	if(!skip){
	  //determine closest value and assign function value for point sampling
	  if(stype==1){ // clip points sampled
	    if(l==0) ll=0; else ll=l-1;
	    if(m==0) mm=0; else mm=m-1;
	  }
	  else if(stype==2){ // need those 4 points to be different! does not hurt if pixel center is not within quad
	    if(l==0){ // need to deal with 1-D problems!!  below assumes 2-D
	      l=1;
	      ll=0;
	    }
	    else ll=l-1;
	    if(m==0){
	      m=1;
	      mm=0;
	    }
	    else mm=m-1;
	  }

	  // determine distance to pixel physical position from quadrant boundings physical grid positions 
	  if((COORD==3)&&(outtype>0)){
	    dist[0]=xp1*xp1+x[wg1][1][l]  *x[wg1][1][l]  -2.0*xp1*x[wg1][1][l]  *cos(xp2-x[wg2][2][m]  ); // (l,m) -> (xp1,xp2)
	    dist[0]=sqrt((double)dist[0]);
	    dist[1]=xp1*xp1+x[wg1][1][ll] *x[wg1][1][ll] -2.0*xp1*x[wg1][1][ll] *cos(xp2-x[wg2][2][m] ); // (ll,m) -> (xp1,xp2)
	    dist[1]=sqrt((double)dist[1]);
	    dist[2]=xp1*xp1+x[wg1][1][l] *x[wg1][1][l] -2.0*xp1*x[wg1][1][l] *cos(xp2-x[wg2][2][mm] ); // (l,mm) -> (xp1,xp2)
	    dist[2]=sqrt((double)dist[2]);
	    dist[3]=xp1*xp1+x[wg1][1][ll]*x[wg1][1][ll]-2.0*xp1*x[wg1][1][ll]*cos(xp2-x[wg2][2][mm]); // (ll,mm) -> (xp1,xp2)
	    dist[3]=sqrt((double)dist[3]);
	  }
	  else{
	    //      dist[0]=sqrt(pow(xp1-x[wg1][1][l],2.)  +pow(xp2-x[wg2][2][m],2.)  );
	    dist[0]=sqrt((xp1-x[wg1][1][l])*(xp1-x[wg1][1][l])
			 +(xp2-x[wg2][2][m])*(xp2-x[wg2][2][m])  );
	    //      dist[1]=sqrt(pow(xp1-x[wg1][1][ll],2) +pow(xp2-x[wg2][2][m],2) );
	    dist[1]=sqrt((xp1-x[wg1][1][ll])*(xp1-x[wg1][1][ll])
			 +(xp2-x[wg2][2][m])*(xp2-x[wg2][2][m]) );
	    //      dist[2]=sqrt(pow(xp1-x[wg1][1][l],2) +pow(xp2-x[wg2][2][mm],2) );
	    dist[2]=sqrt((xp1-x[wg1][1][l])*(xp1-x[wg1][1][l])
			 +(xp2-x[wg2][2][mm])*(xp2-x[wg2][2][mm]) );
	    
	    //      dist[3]=sqrt(pow(xp1-x[wg1][1][ll],2)+pow(xp2-x[wg2][2][mm],2));
	    dist[3]=sqrt((xp1-x[wg1][1][ll])*(xp1-x[wg1][1][ll])
			 +(xp2-x[wg2][2][mm])*(xp2-x[wg2][2][mm]));
	    //fprintf(stderr,"%15.10g %15.10g %15.10g %15.10g\n",dist[0],dist[1],dist[2],dist[3]);
	  }    

	}// if !skip
	else{ // if skip==1
	  l=ll=m=mm=0;
	}
	la[wg1-1][j][i]=l;
	ma[wg2-1][j][i]=m;
	lla[wg1-1][j][i]=ll;
	mma[wg2-1][j][i]=mm;
      }//end 2nd group of firsttime  //----****

      // if found to be outside physical grid, assign to default value
      if( (la[wg1-1][j][i]==0)||(ma[wg2-1][j][i]==0)||(la[wg1-1][j][i]==inn1)||(ma[wg2-1][j][i]==inn2)){
      //      if( (la[wg1-1][j][i]==0)||(ma[wg2-1][j][i]==0)||(la[wg1-1][j][i]==inn1)||(ma[wg2-1][j][i]==inn2)){
	out[j][i]=outerdef;
	//if((i==0)&&(j==N[2]/2)){
	//fprintf(stderr,"%d %d %d %d %d %15.10g %15.10g\n",j,i,m,l,skip,outerdef,out[j][i]);
	  //myexit(1);
	//}
	continue; // skip since not relevant point
      }
      
      // with point sampling just determine closest neighbor and assign that function value to pixel function value
      if(stype==1){
	if(firsttime[wg2-1][wg1-1]==1){  // ****----
	  which=0;
	  for(p=1;p<4;p++){
	    if(dist[p]<dist[which]) which=p;
	  }
	  //now which contains the closest point
	  whicha[wg2-1][wg1-1][j][i]=which;
	}  //----****
	
	// now assign correct function value to pixel function value
	if(whicha[wg2-1][wg1-1][j][i]==0){
	  gi=la[wg1-1][j][i];
	  gj=ma[wg2-1][j][i];
	}
	else if(whicha[wg2-1][wg1-1][j][i]==1){
	  gi=lla[wg1-1][j][i];
	  gj=ma[wg2-1][j][i];
	}
	else if(whicha[wg2-1][wg1-1][j][i]==2){
	  gi=la[wg1-1][j][i];
	  gj=mma[wg2-1][j][i];
	}
	else if(whicha[wg2-1][wg1-1][j][i]==3){
	  gi=lla[wg1-1][j][i];
	  gj=mma[wg2-1][j][i];
	}
	out[j][i]=in[gj][gi];
	
      } //end if stype==1
      else if(stype==2){
	if(firsttime[wg2-1][wg1-1]==1){  // ****----
	  //with plane interpolation we first determine the 3 closest points so we can form a plane
	  which=0;
	  for(p=1;p<4;p++){
	    if(dist[p]>dist[which]) which=p;
	  }
	  whicha[wg2-1][wg1-1][j][i]=which; // the point farthest away
	}  //----****

	// now which contains the ignorable point
	whichtemp=whicha[wg2-1][wg1-1][j][i];
	if(whichtemp==0){ // choose points 1->1 2->2 3->3   that is: (point#->dist#s index)
	  a3=in[mma[wg2-1][j][i]][la[wg1-1][j][i]]-in[ma[wg2-1][j][i]][lla[wg1-1][j][i]]; //2-1
	  b3=in[mma[wg2-1][j][i]][lla[wg1-1][j][i]]-in[ma[wg2-1][j][i]][lla[wg1-1][j][i]]; //3-1

	  x1=x[wg1][1][ll]*sin(x[wg2][2][m]); //1
	  x2=x[wg1][1][l]*sin(x[wg2][2][mm]); //2
          x3=x[wg1][1][ll]*sin(x[wg2][2][mm]); //3
	  y1=x[wg1][1][ll]*cos(x[wg2][2][m]);
	  y2=x[wg1][1][l]*cos(x[wg2][2][mm]);
	  y3=x[wg1][1][ll]*cos(x[wg2][2][mm]);
	  z1=in[ma[wg2-1][j][i]][lla[wg1-1][j][i]];
	  z2=in[mma[wg2-1][j][i]][la[wg1-1][j][i]];
	  z3=in[mma[wg2-1][j][i]][lla[wg1-1][j][i]];
	}
	else if(whichtemp==1){ // choose points 1->0 2->2 3->3
	  a3=in[mma[wg2-1][j][i]][la[wg1-1][j][i]]-in[ma[wg2-1][j][i]][la[wg1-1][j][i]]; //2-0
	  b3=in[mma[wg2-1][j][i]][lla[wg1-1][j][i]]-in[ma[wg2-1][j][i]][la[wg1-1][j][i]]; //3-0

	  x1=x[wg1][1][l]*sin(x[wg2][2][m]); //0  l-m
	  x2=x[wg1][1][l]*sin(x[wg2][2][mm]); //2 l-mm
          x3=x[wg1][1][ll]*sin(x[wg2][2][mm]); //3 ll-mm
	  y1=x[wg1][1][l]*cos(x[wg2][2][m]);
	  y2=x[wg1][1][l]*cos(x[wg2][2][mm]);
	  y3=x[wg1][1][ll]*cos(x[wg2][2][mm]);
	  z1=in[ma[wg2-1][j][i]][la[wg1-1][j][i]];
	  z2=in[mma[wg2-1][j][i]][la[wg1-1][j][i]];
	  z3=in[mma[wg2-1][j][i]][lla[wg1-1][j][i]];

	}
	else if(whichtemp==2){ // choose points 1->1 2->0 3->3
	  a3=in[ma[wg2-1][j][i]][la[wg1-1][j][i]]-in[ma[wg2-1][j][i]][lla[wg1-1][j][i]]; //0-1
	  b3=in[mma[wg2-1][j][i]][lla[wg1-1][j][i]]-in[ma[wg2-1][j][i]][lla[wg1-1][j][i]]; //3-1


	  x1=x[wg1][1][l]*sin(x[wg2][2][m]); //0  l-m
	  x2=x[wg1][1][ll]*sin(x[wg2][2][m]); //1 ll-m
          x3=x[wg1][1][ll]*sin(x[wg2][2][mm]); //3 ll-mm
	  y1=x[wg1][1][l]*cos(x[wg2][2][m]);
	  y2=x[wg1][1][ll]*cos(x[wg2][2][m]);
	  y3=x[wg1][1][ll]*cos(x[wg2][2][mm]);
	  z1=in[ma[wg2-1][j][i]][la[wg1-1][j][i]];
	  z2=in[ma[wg2-1][j][i]][lla[wg1-1][j][i]];
	  z3=in[mma[wg2-1][j][i]][lla[wg1-1][j][i]];

	}
	else if(whichtemp==3){ // choose points 1->1 2->2 3->0
	  a3=in[mma[wg2-1][j][i]][la[wg1-1][j][i]]-in[ma[wg2-1][j][i]][lla[wg1-1][j][i]]; //2-1
	  b3=in[ma[wg2-1][j][i]][la[wg1-1][j][i]]-in[ma[wg2-1][j][i]][lla[wg1-1][j][i]]; //0-1


	  x1=x[wg1][1][l]*sin(x[wg2][2][m]); //0  l-m
	  x2=x[wg1][1][ll]*sin(x[wg2][2][m]); //1 ll-m
          x3=x[wg1][1][l]*sin(x[wg2][2][mm]); //2 l-mm
	  y1=x[wg1][1][l]*cos(x[wg2][2][m]);
	  y2=x[wg1][1][ll]*cos(x[wg2][2][m]);
	  y3=x[wg1][1][l]*cos(x[wg2][2][mm]);
	  z1=in[ma[wg2-1][j][i]][la[wg1-1][j][i]];
	  z2=in[ma[wg2-1][j][i]][lla[wg1-1][j][i]];
	  z3=in[mma[wg2-1][j][i]][la[wg1-1][j][i]];

	}

	a3=y2* z1 - y3* z1 - y1* z2 + y3* z2 + y1* z3 - y2* z3;
	b3=-(x2* z1 - x3* z1 - x1* z2 + x3* z2 + x1* z3 - x2* z3);
	c3=x2* y1 - x3* y1 - x1* y2 + x3* y2 + x1* y3 - x2* y3;
	aoc=a3/c3;
	boc=b3/c3;

	out[j][i]=z1-aoc*(xp1*sin(xp2)-x1)-boc*(xp1*cos(xp2)-y1);
	//out[j][i]=aoc*Aco[wg2-1][wg1-1][j][i]+boc*Bco[wg2-1][wg1-1][j][i];
	//out[j][i]+=in[ma[wg2-1][j][i]][lla[wg1-1][j][i]];

	// if scalar, force to never be negative for pretty output
	if( (wtype==0)&&((wsv==1)||(wsv==2)) ) if(out[j][i]<outerdef) out[j][i]=outerdef;
	  
      } //end if stype==2
    }//end over x1-dir of interpolated data
  }//end over x2-dir


  firsttime[wg2-1][wg1-1]=0;

  firstfirsttime=0;
}//end interp function

#else
// so compiler won't complain
void interpolate(int wtype,int wsv,int stype, int whichg, int inn1, int inn2, int outn1, int outn2, FTYPE in[N2M][N1M], FTYPE out[INTN2][INTN1],FTYPE outerdef, int outtype)
{
}
#endif

FTYPE alnfact(FTYPE N)
{
  return(N*log(N/exp(1.0))+0.5*log(2.*M_PI*N));
}


// which=1, give r
// which=2, give theta
// which=3, give phi
SFTYPE cart2spc(int which, SFTYPE xx, SFTYPE yy, SFTYPE zz)
{
  
  if(which==1){
    return(sqrt(xx*xx+yy*yy+zz*zz));
  }
  else if(which==2){
    if(zz>0.0){
      return(atan(sqrt(xx*xx+yy*yy)/zz));
    }
    else if(zz<0.0){
      return(M_PI+atan(sqrt(xx*xx+yy*yy)/zz));
    }
    else{
      return(M_PI*0.5);
    }
  }
  else if(which==3){ // could also use atan2(y,x) (still requires changes)
    if(xx>0.0){
      if(yy>0.0) return(atan(yy/xx));
      else if(yy<0.0) return(2.0*M_PI+atan(yy/xx));
      else return(0.0);
    }
    else if(xx<0.0){
      return(M_PI+atan(yy/xx));
    }
    else{
      if(yy>0.0) return(M_PI*0.5);
      if(yy<0.0) return(M_PI*1.5);
      else return(0.0);
    }
  }
  else{
    fprintf(fail_file,"no such which=%d\n",which);
    myexit(1);
  }
  return(0);
}


void grids_cart2spc(FTYPE (*sca3) [N3M][N2M][N1M],FTYPE (*vx3) [N3M][N2M][N1M],FTYPE (*vy3) [N3M][N2M][N1M],FTYPE (*vz3) [N3M][N2M][N1M])
{
  int i,j,k;

  LOOPF{
    if(COORD==3){// a=agrid b=bgrid, z=bgrid
      
      sca3[1][k][j][i]=x[2][1][i];
      sca3[2][k][j][i]=x[2][2][j];
      sca3[3][k][j][i]=x[2][3][k];
      
      vx3[1][k][j][i]=x[1][1][i];
      vx3[2][k][j][i]=x[2][2][j];
      vx3[3][k][j][i]=x[2][3][k];
      
      vy3[1][k][j][i]=x[2][1][i];
      vy3[2][k][j][i]=x[1][2][j];
      vy3[3][k][j][i]=x[2][3][k];
      
      vz3[1][k][j][i]=x[2][1][i];
      vz3[2][k][j][i]=x[2][2][j];
      vz3[3][k][j][i]=x[1][3][k];
      
    }
    if(COORD==1){ // a means on vx position, b on vy position, z zone center, c vz position(not needed yet)
      // radius
      sca3[1][k][j][i]=cart2spc(1,x[2][1][i],x[2][2][j],x[2][3][k]);
      vx3[1][k][j][i] =cart2spc(1,x[1][1][i],x[2][2][j],x[2][3][k]);
      vy3[1][k][j][i] =cart2spc(1,x[2][1][i],x[1][2][j],x[2][3][k]);
      vz3[1][k][j][i] =cart2spc(1,x[2][1][i],x[2][2][j],x[1][3][k]);
      
      // theta
      sca3[2][k][j][i]=cart2spc(2,x[2][1][i],x[2][2][j],x[2][3][k]);
      vx3[2][k][j][i] =cart2spc(2,x[1][1][i],x[2][2][j],x[2][3][k]);
      vy3[2][k][j][i] =cart2spc(2,x[2][1][i],x[1][2][j],x[2][3][k]);
      vz3[2][k][j][i] =cart2spc(2,x[2][1][i],x[2][2][j],x[1][3][k]);
      
      // phi
      sca3[3][k][j][i]=cart2spc(3,x[2][1][i],x[2][2][j],0);
      vx3[3][k][j][i] =cart2spc(3,x[1][1][i],x[2][2][j],0);
      vy3[3][k][j][i] =cart2spc(3,x[2][1][i],x[1][2][j],0);
      vz3[3][k][j][i] =cart2spc(3,x[2][1][i],x[2][2][j],0);
      
    }
  }// end coordinate conversions
  
}


// GODMARK
// not correct for multiple cpus

#define SMOOTHLOOP(SMOOTHSIZE) LOOPF for(n=-SMOOTHSIZE*N3NOT1;n<=SMOOTHSIZE*N3NOT1;n++) for(m=-SMOOTHSIZE*N2NOT1;m<=SMOOTHSIZE*N2NOT1;m++) for(l=-SMOOTHSIZE*N1NOT1;l<=SMOOTHSIZE*N1NOT1;l++)

#define SMOOTHCONDITION(SMOOTHSIZE,k,j,i,n,m,l,var,BAD,goodvar) if(\
	     (n+k>=INFULL3)&& \
	     (m+j>=INFULL2)&& \
	     (l+i>=INFULL1)&& \
	     (n+k<OUTFULL3)&& \
	     (m+j<OUTFULL2)&& \
	     (l+i<OUTFULL1) \
	     ){ \
  var[k][j][i]+=goodvar[n+k][m+j][l+i]; \
	  } \
	  else{ \
	    var[k][j][i]+=BAD; \
	  }


void smooth(int SMOOTHSIZE,FTYPE BAD,FTYPE (*var)[N2M][N1M])
{
  int k,j,i,l,m,n;
  FTYPE ftemp;
  FTYPE (*working)[N2M][N1M];

  if(USEMPI) return; // not setup for MPI yet

  working=work1;

  LOOPF{
    working[k][j][i]=0;
  }

  SMOOTHLOOP(SMOOTHSIZE){
    SMOOTHCONDITION(SMOOTHSIZE,k,j,i,n,m,l,working,BAD,var);
  }
  LOOPF{
    // number of points sampled
    ftemp=pow(1.+SMOOTHSIZE*2.0,N1NOT1+N2NOT1+N3NOT1);
    var[k][j][i]=working[k][j][i]/ftemp;
  }

}
