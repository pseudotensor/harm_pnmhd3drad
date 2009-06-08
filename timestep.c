#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif



void timestep(void) 
{
  FTYPE dtother;
  int i,j,k,l ;
  FTYPE idt2[NUMDTCHECKS+1];
  int ks[NUMDTCHECKS+1], js[NUMDTCHECKS+1], is[NUMDTCHECKS+1];
  FTYPE dt2inv_max[NUMDTCHECKS+1] ;
  int didfail,didfail_full;
  int bigger;
  FTYPE finaln;
  static int firsttime=1;
  static FTYPE ttimestep=0,ttimescale=0;
  static FTYPE dtlast;
  char tempc1[50];
  // for slow idtcreate change
  FTYPE bxa,bya,bza,dv,dvdx,delv;
  FTYPE l2_ten;
  FTYPE rho,u ;
  FTYPE odx1,odx2,odx3,ods,odl;
  FTYPE valphen,velfastm,valphen2,cs2 ;
  int reall,viscl,nonvl;
  FTYPE vel1,vel2,vel3;
  FTYPE ftemp;
  FTYPE dvx,dvy,dvz,ftemp1,ftemp2,ftemp3;
  static FTYPE dtrecv;
  int nstepmin,nstepmax;
  int gosub,gosup;
  FTYPE dt2invl[3];
  static FTYPE dtotherlowest;
  static int laststep;



  if(visc_real==1){
    nu_compute();
  }
  if(RESMEM&&(res_real==1)){
    if(rreal==2) current_compute(123);
    nu_res_compute();
  }

  for(l=2;l<=NUMDTCHECKS;l++){
    dt2inv_max[l]=0.;
    ks[l]=js[l]=is[l]=0;
  }
  dtlast = dt ;

  didfail=0;
  didfail_full=0;
  if(firsttime==1){
    ttimestep=t-1.E-12;
    ttimescale=t-1.E-12;
    laststep=0;
  }

  LOOPTIMESTEP{

#if(BOUNDTYPE==3)
    if(bzmask[k][j][i]!=0) continue;
#endif

#if(TS0CHECK)
    if(s[2][k][j][i] < 0) { // actually detects nan too
      sprintf(tempc1,"%3f",s[2][k][j][i]);
      if(tempc1[0]=='n'){
	fprintf(fail_file,"nan internal energy density error: k: %3d j: %3d i: %3d u: %15.10g \n",k,j,i,s[2][k][j][i]) ;
      }
      else if(s[2][k][j][i]<0) fprintf(fail_file,"negative internal energy density error: k: %3d j: %3d i: %3d en: %15.10g \n",k,j,i,s[2][k][j][i]) ;
      didfail=1;
    }
    if(s[1][k][j][i] < 0) {
      sprintf(tempc1,"%3f",s[1][k][j][i]);
      if(tempc1[0]=='n'){
	fprintf(fail_file,"nan mass density error: k: %3d j: %3d i: %3d rho: %15.10g \n",k,j,i,s[1][k][j][i]) ;
      }
      else if(s[1][k][j][i]<0) fprintf(fail_file,"negative mass density error: k: %3d j: %3d i: %3d rho: %15.10g \n",k,j,i,s[1][k][j][i]) ;
      didfail=1;
    }
#endif



    // not inlining this function for some reason, so slow in loop, so make .h file
    // idtcreate(idt2,k,j,i);
#include "timestep.h"

    for(l=2;l<=NUMDTCHECKS;l++){
      ftemp=idt2[l];
      if(ftemp > dt2inv_max[l]){
	dt2inv_max[l] = ftemp ;
	ks[l]=k;
	js[l]=j;
	is[l]=i;
      }
      if(CHECKDTLOW==1){
	if(ftemp>SQIDTLOWEST){
	  timecheck(-l,idt2,k,j,i,0);
	  didfail=1;
	  fflush(fail_file);
	}
      }
    }



  }// end loop over domain


  if(CHECKDTLOW==1){
    // check if any cpu has failure
    if(numprocs>1){
#if(USEMPI)
      MPI_Allreduce(&didfail, &didfail_full, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
    }
    else{
      didfail_full=didfail;
    }
    if(didfail_full){
      if(myid<=0){
	fprintf(log_file,"timestep failure\n");
      }
      if(DOGENDIAG){
	diag(2);
      }
      if(DOAVGDIAG){
	diagavg(2);
      }
#if(USEMPI)
      MPI_Barrier(MPI_COMM_WORLD); // allow to finish diags
#endif
      myexit(5);
    }
  }



  // find lowest constrainer on dt
  reall=2;
  for(l=3;l<=NUMDTCHECKS;l++){
    if(dt2inv_max[l]>dt2inv_max[reall]){
      reall=l;
    }
  }

  if(DODTDIAG){
    // do check up on dominates of timestep for each type
    if((t>ttimestep)||(dt<DTLOWEST)){ // per cpu pure dt data
      for(l=2;l<=NUMDTCHECKS;l++){
	timecheck(l,idt2,ks[l],js[l],is[l],reall);
      }
      fflush(logdt_file);
      
      ttimestep=t+DTtimestep;
    }
  }

  if(DOTSTEPDIAG){
    if(t>ttimescale){
      if(numprocs==1){
	timescale();
	// SUPERMARK -- need to fix timescale to be correct in new cpu setup
	ttimescale=t+DTtimescale;
      }
    }    
  }



  // find lowest constrainer on dt due to visc
  if(dt2inv_max[8]>dt2inv_max[9]){
    viscl=8;
  }
  else viscl=9;

  // find lowest constrainer on dt of non-viscosity type (next highest dt^2)
  nonvl=2;
  for(l=3;l<=NUMDTCHECKS;l++){
    if(l==8) l=10; // skip viscosity
    if(dt2inv_max[l]>dt2inv_max[nonvl]){
      nonvl=l;
    }
  }
  // communicate the lowest dt values to all cpus
  if(numprocs>1){
#if(USEMPI)
    MPI_Allreduce(&(dt2inv_max[reall]), &dt2invl[0], 1, MPI_FTYPE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&(dt2inv_max[viscl]), &dt2invl[1], 1, MPI_FTYPE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&(dt2inv_max[nonvl]), &dt2invl[2], 1, MPI_FTYPE, MPI_MAX, MPI_COMM_WORLD);
#endif
  }
  else{
    dt2invl[0]=dt2inv_max[reall];
    dt2invl[1]=dt2inv_max[viscl];
    dt2invl[2]=dt2inv_max[nonvl];
  }

  if(TRYSUBCYCLE){

    finaln=sqrt(dt2invl[1]/dt2invl[2]);  // fraction of other dt to viscosity dt
    
    if(subcyclen<=1){
      dtotherlowest=1.E+9; // used to check on subcycling below
      gosub=0; // assume by default won't be able to subcycle
      // if visc limits entire comp grid by factor of 2 in dt or more on other dts, then do subcycle
      if((!laststep)&&(finaln>=2.0) ){
	
	// check if stable enough to subcycle
	//      if(fabs( (dt2invl[1]-dtlast)/dtlast)<5.0){
	// so stable, now find next lowest dt
	
	dt=1.0/sqrt(dt2invl[1]); // dt to be used on viscosity
	// number of subcycles over viscosity allowed given current estimates of dt
	subcyclen=(int)(floor(finaln)); // floor to be conservative on fraction
	
	// check if subcycle possible
	if(subcyclen>=2){// should be true!
	  // setup subcycle
	  tscycleto=t+dt*(FTYPE)(subcyclen); // time to cycle to if stable cycle
	  tscyclefrom=t; // time starting subcycle
	  dtlastscycle=1.0/sqrt(dt2invl[2]); // need to make sure not cycling past newest other dts
	  nthsubcycle=1;  // first subcycle is now
	  gosub=1;
	}
	else{
	  fprintf(fail_file,"Unexpected failure in subcycle code: finaln: %15.10g subcyclen: %d\n",finaln,subcyclen);
	  myexit(1);
	}
	//}// end if stable to subcycle
      }// end if viscosity limit and want to try to subcycle on it
      if(gosub==0){ // if no subcycling possible
	dt=1.0/sqrt(dt2invl[0]); // normal case of no subcycling
	subcyclen=1;
	nthsubcycle=0;
      }
    }// endif not subcycling
    else{ // if currently subcycling
      gosup=0; // assume no need to supercycle yet
      nthsubcycle++;
      // check to see if done with subcycling or need to quit
      
      // check to see if visc no longer limit and so supercycle(do all but visc up to viscs time)
      if(finaln<=2.0){
	gosup=1;
      }
      else{// visc still limit
	dt=1.0/sqrt(dt2invl[0]); // trial dt assuming still going to subcycle(should be same as [1])
	dtother=1.0/sqrt(dt2invl[2]);
	if(dtother<dtotherlowest){
	  dtotherlowest=dtother;
	}
	// check if prospective timestep for viscosity still keeps other terms t0+dt further down t to avoid overstepping the other limits based on current data
	if( (t+dt)>(tscyclefrom+dtotherlowest) ){
	  gosup=1;
	}
      }
      if(gosup){ // general setup for supercycle
	dt=(t-tscyclefrom);
	subcyclen=-1;
	tscycleto=t;
	t=tscyclefrom;
      nthsubcycle=0;
      }
    }// endif was/still are subcycling
  }
  else{
    dt=1.0/sqrt(dt2invl[0]); // normal case of no subcycling
  }
  if(analoutput==6){
    // for checking visc code
    dt=pow(invcour2*alpha_real/(dx[2][1][0]*dx[2][1][0]),-1.0);
    ftemp=pow(invcour2*alpha_real/(x[2][1][0]*x[2][1][0]*dx[2][2][N2/2]*dx[2][2][N2/2]),-1.0);
    if(ftemp<dt) ftemp=dt;
  }
 
#if(SUBCOOL)
   coolingKaistep(1);
   if(1.>pow((dtsmall/dt),0.3))dt=pow((dtsmall/dt),0.3)*dt; 
#endif


#if(USEMPI)
#if(DEBUGMPI>0)
      // first check if correct place in code(debug)
      MPI_Allreduce(&nstep, &nstepmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&nstep, &nstepmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (nstep!=nstepmin)||(nstep!=nstepmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d dt nstep: %d\n",myid,nstep);
	fflush(fail_file);
      }      
#endif

      // don't need anymore
      //  MPI_Allreduce(&dt, &dtrecv, 1, MPI_FTYPE, MPI_MIN, MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
      //dt=dtrecv;
#endif

  // because first time step is bad if e.g. viscosity on and no v at first
  if(firsttime==1){ // tweak for given problem so starts out good
    if(dt>1.E-6) dt=1.E-6;
    firsttime=0;
  }


  // don't increase timestep by too much on subcycle or normal cycle.
  // don't check if supercycle since need to force non-visc back to visc time.

  if(subcyclen>=0){
    if(dt > 1.3*dtlast) dt = 1.3*dtlast ;
  }

  /* don't step beyond end of run */
  if(t + dt >= tf){
    // last timestep
    laststep=1;
    if(subcyclen==1){
      dt = tf - t ;
      reallaststep=1;
    }
    if(subcyclen==-1){
      fprintf(fail_file,"shouldn't be here at end of run at super cycle\n");
      myexit(1);
    }
    if(subcyclen>=2){
      // just end subcycle and let next timestep() figure final dt, never subcycling again
      dt=(t-tscyclefrom);    
      subcyclen=-1;
      tscycleto=t;
      t=tscyclefrom;
      nthsubcycle=0;
      reallaststep=0;
    }
    // make sure don't get pathological case of dt=0 on last step
    if(dt<SSMALL){
      reallaststep=1;
      laststep=1;
      dt=SSMALL;
    }
  }
}



void idtcreate(FTYPE*idt2,int k, int j, int i)
{
  FTYPE bxa,bya,bza,dv,dvdx,delv;
  FTYPE rho,u ;
  FTYPE l2_ten;
  FTYPE Length;
  FTYPE odx1,odx2,odx3,ods,odl;
  FTYPE valphen,velfastm,valphen2,cs2 ;
  FTYPE ftemp;
  FTYPE vel1,vel2,vel3;
  FTYPE dvx,dvy,dvz,ftemp1,ftemp2,ftemp3;

#include "timestep.h"
}


#define NUMPATHS 5

void timescale(void)
{
  int i,j,k,l,m;
  FTYPE bxa,bya,bza,dv,dvdx,delv;
  FTYPE rho,u ;
  FTYPE l2_ten;
  FTYPE Length;
  FTYPE odx1,odx2,odx3,ods,odl;
  FTYPE valphen,velfastm,valphen2,cs2 ;
  FTYPE ftemp;
  FTYPE vel1,vel2,vel3;
  FTYPE dvx,dvy,dvz,ftemp1,ftemp2,ftemp3;
  FTYPE idt2[NUMDTCHECKS+1];
  static int firsttime=1;
  static FTYPE paths[NUMPATHS];
  static FTYPE limits[2]; // inner and outer x1 limits
  FTYPE timescales[NUMDTCHECKS+1][NUMPATHS]; // 5 time scale paths
  static char filename[MAXFILENAME];
  char temps[50];
  static FILE * timescale_file;
  static int startm,lastm;
  int gotfirst,gotlast,testm;
#if(USEMPI)
  static MPI_Request request[2];
#endif

  if(firsttime==1){

    // limits on x1-range
    limits[0]=L[1][1];
    limits[1]=R0; // torus center or injection center

    // locations on x2
    paths[0]=M_PI*0.5-M_PI*0.5*2.0/3.0;
    paths[1]=M_PI*0.5-M_PI*0.5/3.0;
    paths[2]=M_PI*0.5;
    paths[3]=M_PI*0.5+M_PI*0.5/3.0;
    paths[4]=M_PI*0.5+M_PI*0.5*2.0/3.0;

    // determine what cpu gets what trajectories
    gotfirst=-1;
    gotlast=-1;
    startm=0;
    lastm=4;
    k=0;
    for(m=0;m<=NUMPATHS-1;m++){
      for(j=0;j<N2;j++){

	// 0.6 might get 2 trajectories, but so close that's ok anyways since overwritten by 2nd one
	if(fabs(x[2][2][j]-paths[m])<0.6*dx[1][2][j]){ // if within zone

	  if(gotfirst<0) gotfirst=m;// gets which m is first m for this cpu
	  gotlast=m; // gets last gotten
	
	} // end over this trajectory if found
	
      } // end over x2-dir seeking trajectories
    }// end over all trajectories
      
    // now determine which cpus should do what trajectories
    if(numprocs>1){
      lastm=gotlast;
      if(myid==0){
	startm=gotfirst;
      }
#if(USEMPI)
      if(myid<numprocs-1){
	MPI_Isend(&lastm,1,MPI_INT,myid+1,myid+1,MPI_COMM_WORLD,&request[0]);
      }
      if(myid>0){
	MPI_Irecv(&testm,1,MPI_INT,myid-1,myid,MPI_COMM_WORLD,&request[1]);
	MPI_Wait(&request[1],&mpichstatus);
	fprintf(log_file,"testm: %d\n",testm);
	fflush(log_file);
	startm=testm+1;
	// if startm>lastm it won't go anywhere below anyways
      }
      MPI_Wait(&request[0],&mpichstatus);
#endif
    }
    else{
      startm=gotfirst;
      lastm=gotlast;
    }
    fprintf(log_file,"startm: %d lastm: %d\n",startm,lastm);
    fflush(log_file);
    sprintf(filename,"%s0_timescales%s",DATADIR,DATEXT) ;


    // write file header
    if(myid<=0){
      if((timescale_file=fopen(filename,WRITETYPE))==NULL){ // naively append if appendold=1
	fprintf(fail_file,"timescales: Cannot open: %s\n",filename);
	myexit(1);
      }
      if(appendold==0){
	fprintf(timescale_file,"#%10s\n%10d %10d\n","TSVER",TSVER,TSTYPE);
	fprintf(timescale_file,"#%15s","time");
	for(m=startm;m<=lastm;m++){
	  for(l=2;l<=NUMDTCHECKS;l++){
	    sprintf(temps,"p%1d-c%2d",m,l);
	    fprintf(timescale_file," %15s",temps);
	  }
	}
	fprintf(timescale_file,"\n");
      }
      fclose(timescale_file);
    }
  }// endif firsttime==1


  // now find different timescales along different trajectories
  // reset timescales
  for(l=2;l<=NUMDTCHECKS;l++){
    for(m=0;m<=NUMPATHS-1;m++){
      timescales[l][m]=0;
    }
  }
  k=0;
  for(m=startm;m<=lastm;m++){
    for(j=0;j<N2;j++){

      if(fabs(x[2][2][j]-paths[m])<0.6*dx[1][2][j]){ // if within zone
	
	for(i=0;i<N1;i++){
	  
	  if( (x[2][1][i]>limits[0])&&(x[2][1][i]<limits[1]) ){ // if within bounds of trajectory

#include "timestep.h"
	    
	    for(l=2;l<=NUMDTCHECKS;l++){
	      timescales[l][m]+=1.0/sqrt(SSMALL+idt2[l]); // add up dt for this zone
	    }
	    
	  } // end if within bounds of trajectory
	  
	}// end loop over this trajectory
	
      } // end over this trajectory if found

    } // end over x2-dir seeking trajectories

  }// end over all trajectories

  // now write timescales to file
  for(i=0;i<numprocs;i++){
    if(myid==i){
      if((timescale_file=fopen(filename,"at"))==NULL){
	fprintf(fail_file,"timescales: Cannot open: %s\n",filename);
	myexit(1);
      }
      if(i==0){
	fprintf(timescale_file," %15.10g",t);
      }
      for(m=startm;m<=lastm;m++){
	for(l=2;l<=NUMDTCHECKS;l++){
	  fprintf(timescale_file," %15.10g",timescales[l][m]);
	}
      }
      if(i==numprocs-1){
	fprintf(timescale_file,"\n");
      }
      fclose(timescale_file);
    }
#if(USEMPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  firsttime=0;

}

// failmode==0 fail
// failmode=# which to check
void timecheck(int failmode,FTYPE*idt2,int k, int j, int i,int reall)
{
  FTYPE ftemp;
  int wfail,l;
  FTYPE dtinv;
  FTYPE valphen,velfastm,valphen2,cs2 ;
  FILE* out;
  static int firsttime=1;
  FTYPE vel1,vel2,vel3;
  FTYPE ftemp1,ftemp2,ftemp3;
  FTYPE u,rho;
  FTYPE bxa,bya,bza,dv,dvdx,delv;  


  if(failmode>0){
    idtcreate(idt2,k,j,i);
    dtinv = sqrt(idt2[failmode]);
  }
  else{
    // no need to create idt's if failmode<0 since just created it in loop
    dtinv=sqrt(idt2[-failmode]);
    failmode=0;//reactivate fail mode
  }

  u = s[2][k][j][i];
  rho=s[1][k][j][i];
  if(wgam) cs2 = gam*(gam-1.)*s[2][k][j][i]/rho;
  else cs2=cs*cs;

  if(mag==1){
    /* alfven velocity */
    bxa = e2z_1(v[2][1],k,j,i);
    bya = e2z_2(v[2][2],k,j,i);
    bza = e2z_3(v[2][3],k,j,i);
    
#if(ALFVENLIMIT==0)
    valphen2 = (bxa*bxa + bya*bya + bza*bza)/rho ; // max alfven speed, squared, in the zone along length of field line
#else

#if(DOINVSOL2FUN)
    invsol2=invsol2fun[k][j][i];
#endif

    ftemp=bxa*bxa + bya*bya + bza*bza; // b^2
    valphen2 = ftemp/(rho+ftemp*invsol2) ; // max alfven speed, squared
#endif
  }
  else{
    valphen2=SSMALL;
  }
  velfastm=sqrt(valphen2+cs2);

  vel1 =e2z_1(v[1][1],k,j,i)-vg[1];
  vel2 =e2z_2(v[1][2],k,j,i)-vg[2] ;
  vel3 =e2z_3(v[1][3],k,j,i)-vg[3] ;
  // effectives
  vel1=fabs(vel1)+velfastm;
  vel2=fabs(vel2)+velfastm;
  vel3=fabs(vel3)+velfastm;
  
  // see who failed or see what dominates for checkup
  wfail=2;
  ftemp=fabs(idt2[wfail]);
  for(l=3;l<=NUMDTCHECKS;l++){
    if(fabs(idt2[l])>ftemp){
      ftemp=fabs(idt2[l]);
      wfail=l;
    }
  }
  ftemp=1.0/(sqrt(ftemp+SSMALL)); // real dt of failure 

  if(wfail>NUMDTCHECKS){
    fprintf(fail_file,"unexpected failure in timestep.c when checking dt\n");
    myexit(1);
  }
  
  
  if((failmode==0)||(DODTDIAG==0)){
    out=fail_file;
    fprintf(fail_file,"dt has dropped below set lowest threshold.  dt=%15.10g\n", 1.0/dtinv);    
  }
  else out=logdt_file;
  

  if(failmode==0){
    fprintf(out,"Time check at t=%15.10g #%d(%d) js%d is%d has dt=%15.10g\n",t,wfail,reall,jsmall,ismall,ftemp);
    if( (COORD==3)&&(SAMPLED>0) ){
      fprintf(out,"Velocities at k=%d j=%d i=%d x=%15.10g z=%15.10g\n",k,j,i,x[1][1][i]*sin(x[1][2][j]),x[1][1][i]*cos(x[1][2][j]));
    }
    else{
      fprintf(out,"Velocities at  k=%d j=%d i=%d x1=%15.10g x2=%15.10g\n",k,j,i,x[1][1][i],x[1][2][j]);
    }
    
    fprintf(out,
	    "#  Vel Checked   %15s : %15s\n"
	    "-------------------------------------\n"
	    "1: sound speed : %15.10g : %15.10g\n"
	    "2: x1-vel+-cs  : %15.10g : %15.10g\n"
	    "3: x2-vel+-cs  : %15.10g : %15.10g\n"
	    "4: v3-vel+-cs  : %15.10g : %15.10g\n"
	    "5: Linear visc : %15.10g : %15.10g\n"
	    "6: visc(x1)    : %15.10g : %15.10g\n"
	    "7:             : %15.10g : %15.10g\n"
	    "8: rvisc(x1):  : %15.10g : %15.10g\n"
	    "9: rvisc(x2):  : %15.10g : %15.10g\n"
	    "10: resist:    : %15.10g : %15.10g\n\n"
						,"dt","v||dv",
						1.0/SSMALL,cs,
						1./(fabs(sqrt(idt2[2]))+SSMALL),vel1,
						1./fabs(sqrt(idt2[3])+SSMALL),vel2,
						1./fabs(sqrt(idt2[4])+SSMALL),vel3,
						1./fabs(sqrt(idt2[5])+SSMALL),invcour2*nu_l*cs,
						1./fabs(sqrt(idt2[6])+SSMALL),v[1][1][k][j][ip1]-v[1][1][k][j][i],
						1./fabs(sqrt(idt2[7])+SSMALL), v[1][2][k][jp1][i]-v[1][2][k][j][i],
						1./fabs(sqrt(idt2[8])+SSMALL),s[1][k][j][i],s[2][k][j][i],nu_real[k][j][i],s[1][k][j][i+1],s[2][k][j][i+1],nu_real[k][j][i+1],v[1][1][k][j][ip1]-v[1][1][k][j][i],
						1./fabs(sqrt(idt2[9])+SSMALL), v[1][2][k][jp1][i]-v[1][2][k][j][i],
						1./fabs(sqrt(idt2[10])+SSMALL),invcour2*resist_real0*cs);
    
	fprintf(out, "# s[1][k][j][i]  s[1][k][j+1][i]  s[1][k][j][i+1]  s[1][k][j-1][i]  s[1][k][j][i-1] \n"
		     "%16.10g:%16.10g%16.10g:%16.10g:%16.10g\n"
		     "# s[2][k][j][i]  s[2][k][j+1][i]  s[2][k][j][i+1]  s[2][k][j-1][i]  s[2][k][j][i-1] \n"
                     "%16.10g:%16.10g%16.10g:%16.10g:%16.10g\n",
			s[1][k][j][i],s[1][k][j+1][i],s[1][k][j][i+1],s[1][k][j-1][i],s[1][k][j][i-1],
			s[2][k][j][i],s[2][k][j+1][i],s[2][k][j][i+1],s[2][k][j-1][i],s[2][k][j][i-1]);
  }
  else{
    if( (firsttime==1)&&(appendold==0)){
      fprintf(out,"#%10s\n%10d %10d\n","LOGDTVER",LOGDTVER,LOGDTTYPE);
      fprintf(out,"#%15s %15s %5s %5s %5s %5s %5s %5s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n","t","dt","r","w","l","k","j","i","x1","x2","cs_dt","cs_v","x1v_dt","x1v_(v+-cs)","x2v_dt","x2v_(v+-cs)","x3v__dt","x3v__(v+-cs)","lv_dt","lv_dv","vx1_dt","vx1_dv","NONE","NONE","rvx1_dt","rvx1_nu","rvx2_dt","rvx2_nu","resist_dt","resist_v");
    }

    fprintf(out," %15.10g %15.10g %5d %5d %5d",t,1./(fabs(sqrt(idt2[failmode]))+SSMALL),reall,wfail,failmode);
    if( (COORD==3)&&(SAMPLED>0) ){
      fprintf(out," %5d %5d %5d %15.10g %15.10g",k,j,i,x[1][1][i]*sin(x[1][2][j]),x[1][1][i]*cos(x[1][2][j]));
    }
    else{
      fprintf(out," %5d %5d %5d %15.10g %15.10g",k,j,i,x[1][1][i],x[1][2][j]);
    }
    fprintf(out," %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n"

#if(VISCMEM)
	    ,1.0/SSMALL,cs
	    ,1./(fabs(sqrt(idt2[2]))+SSMALL),vel1
	    ,1./fabs(sqrt(idt2[3])+SSMALL),vel2
	    ,1./fabs(sqrt(idt2[4])+SSMALL),vel3
	    ,1./fabs(sqrt(idt2[5])+SSMALL),invcour2*nu_l*cs
	    ,1./fabs(sqrt(idt2[6])+SSMALL),v[1][1][k][j][ip1]-v[1][1][k][j][i]
	    ,1./fabs(sqrt(idt2[7])+SSMALL), v[1][2][k][jp1][i]-v[1][2][k][j][i]
	    ,1./fabs(sqrt(idt2[8])+SSMALL),nu_real[k][j][i]
	    ,1./fabs(sqrt(idt2[9])+SSMALL), nu_real[k][j][i]
	    ,1./fabs(sqrt(idt2[10])+SSMALL),invcour2*resist_real0*cs,diagn[2][k][j][i],tot[0][k][j][i],tot[1][k][j][i]
#else
	    ,1.0/SSMALL,cs
	    ,1./(fabs(sqrt(idt2[2]))+SSMALL),vel1
	    ,1./fabs(sqrt(idt2[3])+SSMALL),vel2
	    ,1./fabs(sqrt(idt2[4])+SSMALL),vel3
	    ,1./fabs(sqrt(idt2[5])+SSMALL),invcour2*nu_l*cs
	    ,1./fabs(sqrt(idt2[6])+SSMALL),v[1][1][k][j][ip1]-v[1][1][k][j][i]
	    ,1./fabs(sqrt(idt2[7])+SSMALL),v[1][2][k][jp1][i]-v[1][2][k][j][i]
	    ,1./fabs(sqrt(idt2[6])+SSMALL),v[1][3][kp1][j][i]-v[1][3][k][j][i]
	    ,1./fabs(sqrt(idt2[9])+SSMALL),0.0
	    ,1./fabs(sqrt(idt2[10])+SSMALL),invcour2*resist_real0*cs,diag2[k][j][i],tot[0][k][j][i],tot[1][k][j][i]
#endif
  );
  }


  firsttime=0;
}

