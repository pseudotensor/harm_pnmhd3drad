#include "global.h"

#if(LINUXCLUSTER==1)
#include "decs.h"
#endif

#include "defs.h"


int main(
	int argc,
	char *argv[],
	char *envp[]
	)
{
  FTYPE ftemp;
  int itemp;
  int dostepout;
  int firstptr=1;
  int k,j,i;
  int ii;
  SFTYPE tlasttime;
  SFTYPE tnext,tnextavg ;
  static int nsteptot=0 ;
  static int gnstep = 0,gnsteptot=0 ;
  static long subnstep = 0,supnstep=0,nornstep=0 ;
  static int gsubnstep = 0 ;
  int error=0,errorcode=0;
  FILE *floor_file,*gogo_file,*gocont_file;
  FILE *fileit;
  char floorfname[50],stemp[50];
  int flooroutc;
  SFTYPE tfloor;
  char goch;
  int goend;
  int avgdone=0;
#if(TIMEMETHOD==0)
  time_t timestart,timestop;
  time_t gtimestart,gtimestop,checktime;
#elif(TIMEMETHOD==1)
  struct timeval timestart,timestop, gtimestart,gtimestop,checktime; struct timezone tz;
#elif((TIMEMETHOD==2)||(TIMEMETHOD==3))
  clock_t timestart,timestop, gtimestart,gtimestop,checktime;
#endif
  SFTYPE walltime=0,walltimelocal=0,walltot=0;
	// general time reports
	clock_t usertmstimestart,usertmstimestop,systmstimestart,systmstimestop;
	struct timeval wttimestart,wttimestop;
#if(TIMEMETHOD!=1)
	struct timezone tz;
#endif
	
  static int logcount=0;
  static int avgcount=0;
  static int floorcount=0;
  int numzones;
  FILE* cpuout;
  FILE* perfout;
  SFTYPE comptstart;
  int imagescaleflag=0;
  char temps[MAXFILENAME];
  int badlog[NUMLOGCHECKS];

  ///////////////////  INITIALIZATION

  // uncomment if get SIGFPE 8 error at run with floats and no result problem
  signal(8,SIG_IGN);


  fprintf(stderr,"START\n");


  fflush(stderr);
  if(POSTPROC){
    fprintf(stderr,"Not designed to postprocess during a run: forgot to turn off POSTPROC in global.h?\n");
    myexit(1);
  }



  // global parameter variables
  gpar_init();




#if(USEMPI)
  // initialize MPI
  init_MPI(argc,argv);
#endif

  init_genfiles(0); // setup general log files

  // output num cpus and prep finalperf file
  if(myid<=0){
    sprintf(temps,"%snumcpus%s",DATADIR,".txt") ;
    if((cpuout=fopen(temps,"wt"))==NULL){
      fprintf(stderr,"fail: Cannot open: %s\n",temps);
      myexit(1);
    }
    else{
      fprintf(cpuout,"%d\n",numprocs);
      fclose(cpuout);
    }
    if(PERFTEST){
      sprintf(temps,"%sfinalperf%s",DATADIR,".txt") ;
      if(!(perfout=fopen(temps,"wt"))){
				fprintf(fail_file,"Can't open %s\n",temps);
				myexit(1);
      }
    }

  }


  if(myid<=0){
    if(CHECKCONT){
      
      sprintf(stemp,"%sgo.cont",DATADIR);
      
      if((gocont_file=fopen(stemp,"rt"))==NULL){
				fprintf(fail_file,"Could not open go.cont file: %s\n",stemp);
				fprintf(fail_file,"Did you copy over base files?\n");
				//	myexit(1); // can't exit yet if want clean MPI exit
				goch='z';
      }
      else goch='a';
    }
  }
#if(USEMPI)
  MPI_Bcast(&goch,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
  if(goch=='z') myexit(1);
  if(myid<=0){
    if(CHECKCONT){
      goch=fgetc(gocont_file);
      if( (goch=='y')||(goch=='Y')){
				gocont=1;
				fprintf(stderr,"go.cont called\n");
				fprintf(logfull_file,"#go.cont called\n");
				fflush(logfull_file);

				fscanf(gocont_file,"%d %d",&runtype,&directinput);
#if(SENSITIVE==1)
#define INPUTCONT1 "%lf"
#else
#define INPUTCONT1 "%f"
#endif
				if(directinput==3){
					fscanf(gocont_file,INPUTCONT1,&timereenter);
				}
      }
      fclose(gocont_file);
    }
  }
  if(numprocs>1){
#if(USEMPI)
    MPI_Bcast(&gocont,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&runtype,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&directinput,1,MPI_INT,0,MPI_COMM_WORLD);
    if(directinput==3){
      MPI_Bcast(&timereenter,1,MPI_SFTYPE,0,MPI_COMM_WORLD);
    }
#endif
  }

  // start interps fresh no matter what.
  globalinterpmod=1;
  /* Initialize ALL */
  if(myid<=0){
    fprintf(logfull_file,"#Start Initialization.\n") ;
    fflush(logfull_file);
  }
  fprintf(stderr,"proc: %s Start Initialization.\n",myidtxt) ;
  if((error=init(argc,argv,envp))>0){
    fprintf(fail_file,"\nError: initialize_par: %d\n",error);
    myexit(1);
  }
  if(TVDLF){
    firsttimetvdlf=1;
  }
  floorcount=0;
  logcount=0;
  nstep=0;
  totstep=0.;
  flooroutc=floor_start;
  avgdone=0;
  avgcount=0;
  logcount=0;
  diagstep=0.;
  tnext=tstart;
  tnext = tstart+(FTYPE)(logcount)*DTl ;
  tnextavg = tavgi ;
  tlasttime=0;
  comptstart=t;

  fprintf(stderr,"proc: %s   End Initialization.\n",myidtxt) ;
  if(myid<=0){
    fprintf(logfull_file,"#End Initialization.\n") ;
    fflush(logfull_file);
  }

  // write version header info for perf and step logs
  if((myid<=0)&&(appendold==0)){
    if(DOLOGSTEP){
      fprintf(logstep_file,"#%10s\n%10d %10d\n","STEPVER",STEPVER,STEPTYPE); fflush(logstep_file);
    }
    if(DOLOGPERF){
      fprintf(logperf_file,"#%10s\n%10d %10d\n","PERFVER",PERFVER,PERFTYPE); fflush(logperf_file);
    }
  }
  /* do initial diagnostics */
  if(DOGENDIAG){
    diag(0) ;
  }
  logcount++;
  tnext = tstart+(FTYPE)(logcount)*DTl ;

  if(DOFLOORDIAG>=1){
    tfloor=t-1.E-12;
    for(i=0;i<NUMFLOOROUT;i++){
      for(j=1;j<=NUMFLOORVAR;j++){
				floorcnt[i][j]=0;
      }
    }
  }


  dtsmall=0.;
  //////////////////   BEGIN COMPUTATION

  fprintf(stderr,"Starting Computation on proc: %s . . .\n",myidtxt);
  if(myid<=0){
    fprintf(logfull_file,"#Starting Computation\n");
    fflush(logfull_file);
  }
  GETTIME(&timestart);
  GETTIME(&gtimestart);
	microtime(&wttimestart);
	myustimes2(&usertmstimestart,&systmstimestart);
	
  goend=0;
  subcyclen=1; // start fresh
  reallaststep=0;
  if(DOLOGSTEP){
    if(myid<=0){ fprintf(logstep_file,"#");}
  }
  //  while(t < tf) {

  tdep_compute(); // find new time dep stuff

  while(reallaststep==0){
    
    /*
			if( (t>.4926)&&(firstptr==1)){
      firstptr=0;
			}
    */
    //    ptraddr(nstep);
/////////////// for diag /////////
//    if(t>0.8){DTd=1.e-10;}
///////////////////


    if(TVDLF==0){
      /* find timestep */
#if(!SUBCOOL)
      coolingKaistep(1);
#endif


      timestep() ;
      //dt=3.78E-5;
      /* step variables forward in time */
      if(COMPDIM<=2){
				stepvar_2d() ; // sets vars and t
      }
      else{
				stepvar_3d(); // sets vars and t
      }
    }
    else{
      if(firsttimetvdlf==1){ // only need to do this once at first
				zeus2tvdlf();
				firsttimetvdlf=0; // assumes no inflow condition where analytic's might be needed to loop here.
      }
      //    fprintf(stderr,"%21.15g %21.15g\n",s[1][0][0][0],p[0][0][0]);
      //    fflush(stderr);
      // sets dt and vars and t
      steptvdlf();
      if(t>=tf) reallaststep=1;
    }
    
    tdep_compute(); // find new time dep stuff
    
    if(DOSPDIAG){
      if(TVDLF){
				//tvdlf2zeus(); // don't care about sp stuff right now
      }
      if(mycpupos[1]==0){ // only inner x1 cpus relevant
				sp_compute(); // check on sonic point
      }
    }
//    analsolve(0); // find new analytic values
    if(TVDLF){
      //zeus2tvdlf(); // don't care about analytic time dependence right now
    }

  
    //    if(subcyclen<=1) nstep++ ;
    nstep++ ; // any type of step
    totstep+=substepm;
    if(subcyclen>1) subnstep++;
    if(subcyclen==-1) supnstep++;
    if(subcyclen==1) nornstep++;

    /* Perform required diagnostics that need to be done per timestep */
    // only do diags if just supercycled or normal cycle
    if(subcyclen<=1){

      //perform general diagnostics      
      if(DOGENDIAG){
				//      if( (t >= tnext)||(dt<DTLOWDUMP)||(DTl<0.0)) {
//			if(nstep*1.>diagstep){
				diag(1) ;
//				diagstep+=5.e4;
//			}
				//logcount++;
				//tnext = tstart+(FTYPE)(logcount)*DTl ;
				// }
      }


      //perform average diagnostics
      if(0&&DOAVGDIAG){ //SUPERMARK -- fix avgdiag for new cpu setup
				if( (t>=tnextavg)&&(t<tavgf)){
					diagavg(1) ;
					avgcount++;
					tnextavg = tavgi+(FTYPE)(avgcount)*DTavg ;
				}
				if( (t>=tavgf)&&(avgdone==0)){
					diagavg(2);
					avgdone=1;
				}
      }

      // floor diagnostics
      if(DOFLOORDIAG>=1){
				if(t>tfloor){
					flooroutc++;
					floorcount++;
					tfloor = tstart+(FTYPE)(floorcount)*DTfloor ;
	  
					fprintf(logfl_file,"floorcnt: t=%15.10g\n",t);
					for(i=0;i<NUMFLOOROUT;i++){
						for(j=1;j<=NUMFLOORVAR;j++){
							fprintf(logfl_file,"%12d ",floorcnt[i][j]);
						}
						fprintf(logfl_file,"\n");
					}
					for(j=1;j<=NUMFLOORVAR;j++){
						fprintf(logfl_file,"lowest %d %15.10g\n",wherelowest[j],floorlowest[j]);
					}
	  
	  
				}
				if(FLUSHFAILDT){
					fflush(logfl_file);
				}
      }
    }// end if normal or super cycle

    // check up on how many timesteps per second so can calibrate period of step, perf, and gocheck
    if(myid<=0){
      if((nstep==1)||( (!(nstep%NTIMECHECK))) ){
				// no need to synch cpus since should be close and no MPI calls used so don't require exact synch
				//just use average time, instead of local time, to do a timestep
				GETTIME(&checktime);
				walltime=(SFTYPE) DELTATIME(checktime,timestart);
				if(walltime<1E-5) walltime=1E-5;
				// now calibrate everything
				NTIMECHECK=(int)((SFTYPE)nstep*DTtimecheck/walltime); // took walltime to go nsteps, so check every so steps corresponding to desired time.
				if(DTstep!=0.0) NDTCCHECK=(int)((SFTYPE)nstep*DTstep/walltime);
				if(DTstepdot!=0.0) NDTDOTCCHECK=(int)((SFTYPE)nstep*DTstepdot/walltime);
				if(DTperf!=0.0) NZCCHECK=(int)((SFTYPE)nstep*DTperf/walltime);
				if(DTgocheck!=0.0) NGOCHECK=(int)((SFTYPE)nstep*DTgocheck/walltime);
	
				if(NTIMECHECK<1) NTIMECHECK=1;
				if(NDTCCHECK<1) NDTCCHECK=1;
				if(NDTDOTCCHECK<1) NDTDOTCCHECK=1;
				if(NZCCHECK<1) NZCCHECK=1;
				if(NGOCHECK<1) NGOCHECK=1;
	
				//fprintf(stderr,"%d %d %d %d %d\n",NTIMECHECK,NDTCCHECK,NDTDOTCCHECK,NZCCHECK,NGOCHECK);
      }

      // check if time to output step/time/dt info
      // setup so can plot in sm
      
      if(DOLOGSTEP&&( (!(nstep%NDTDOTCCHECK)) )){
				fprintf(logstep_file,"."); fflush(logstep_file);
      }
      if(DOLOGSTEP&&( (!(nstep%NDTCCHECK))||(t>=tf-1.0E-7)||(t<=tstart+1.0E-7) ) ){
				fprintf(logstep_file,"\n");
				for(i=1;i<=1;i++){ // ==0 not done since really not needed
					if(i==0){
						fileit=log_file;
						dostepout=1;
					}
					else if(i==1){
						fileit=logstep_file;
						dostepout=1;
					}
					if(dostepout){
						fprintf(fileit,"#step,norm,sup,sub,t,dt,upto,i,N\n"
										"%10ld %10ld %10ld %10ld %15.10g %15.10g ",nstep,nornstep,supnstep,subnstep,t,dt) ;
						if(subcyclen>1) fprintf(fileit,"%15.10g %5d %5d\n",tscycleto,nthsubcycle,subcyclen) ;
						else if(subcyclen==1) fprintf(fileit,"%15.10g %15.10g %10ld %10.2g 0\n",t+dt,dtsmall,substepm,totstep);
						else if(subcyclen==-1) fprintf(fileit,"%15.10g 0 0\n",tscycleto) ;
						fflush(fileit);
					}
				}
				fprintf(logstep_file,"#");
	
				// check image scale too
				if(imagescaleflag==0){
					if(t>=timagescale){
						fprintf(logfull_file,"#t>timagescale hit: t: %21.15g timagescale: %21.15g\n",t,timagescale);
						imagescaleflag=1;
					}
				}
      }
    
  
      // check if user wants to stop or not(go.go)
      if(!(nstep%NGOCHECK)){
				sprintf(stemp,"%sgo.go",DATADIR);
	
				if((gogo_file=fopen(stemp,"rt"))==NULL){
					fprintf(fail_file,"Could not open go file: %s\n",stemp);
					myexit(1);
				}
				goch=fgetc(gogo_file);
				if( (goch=='n')||(goch=='N')){
					goend=1;
					fprintf(stderr,"go.go called\n");
					fprintf(logfull_file,"#go.go called\n");
					fflush(logfull_file);
				}
				fclose(gogo_file);      
				// if myid!=0 and numprocs>1 then could deal with this, but messy due to timers
	
				if(goend){
					reallaststep=1;
					fprintf(log_file,"#proc: %s go.go called\n",myidtxt);
					fflush(log_file);
				}
      }
  
      // speed check
      // setup so can plot in sm
      if(DOLOGPERF&&(!(nstep%NZCCHECK))){
				GETTIME(&gtimestop);
				// running average
				// running average zonecycle rate
				walltime=(SFTYPE) DELTATIME(gtimestop,timestart);
				if(walltime<1E-5) walltime=1E-5;
				// local zonecycle rate
				walltimelocal=(SFTYPE) DELTATIME(gtimestop,gtimestart);
				if(walltimelocal<1E-5) walltimelocal=1E-5;
				for(i=1;i<=1;i++){ // don't really want perf for i==0
					if(i==0){
						fileit=log_file;
						dostepout=1;
						numzones=N1*N2*N3; // GODMARK not really right, but not used
						strcpy(stemp,"");
					}
					else if(i==1){
						fileit=logperf_file;
						dostepout=1;
						numzones=realtotalzones;
						strcpy(stemp,"all");
					}
					if(dostepout){
						fprintf(fileit,"#t, ete, n, wt, zc, tu/hr,  lete, ln, lwt, lzc, ltu/hr\n");
						fprintf(fileit,"%15.10g %15.10g %10ld %15.10g %10d %10.5g  %15.10g %5d %15.10g %10d %10.5g\n"
										,t
										,((tf-t+1.0E-6)/(t-comptstart+1.0E-6)*walltime*SEC2HOUR)
										,nstep
										,walltime*SEC2HOUR
										,(int)((FTYPE)(numzones)*(FTYPE)(nstep)/walltime)
										,((t-comptstart)/(walltime*SEC2HOUR))
		    
										,((tf-t+1.0E-6)/(t-tlasttime+1.0E-6)*walltimelocal*SEC2HOUR)		  
										,NZCCHECK
										,walltimelocal*SEC2HOUR
										,(int)((FTYPE)(numzones)*(FTYPE)(NZCCHECK)/walltimelocal)
										,((t-tlasttime)/(walltimelocal*SEC2HOUR))
							);
						ftemp=walltimelocal/(t-tlasttime); // seconds per time unit
						if(ftemp<LOGLWTCHECK) badlog[0]=1; else badlog[0]=0; 
						if(ftemp<LOGDWTCHECK) badlog[1]=1;else badlog[1]=0;
						if(ftemp<LOGIWTCHECK) badlog[2]=1;else badlog[2]=0;
						if(ftemp<LOGENERWTCHECK) badlog[3]=1;else badlog[3]=0;
						if(ftemp<LOGLOSSWTCHECK) badlog[4]=1;else badlog[4]=0;
						if(ftemp<LOGTIMESTEPWTCHECK) badlog[5]=1;else badlog[5]=0;
						if(ftemp<LOGSPWTCHECK) badlog[6]=1;else badlog[6]=0;
						for(ii=0;ii<NUMLOGCHECKS;ii++){
							// output to perf file if that log is impacting performance
							if(badlog[ii]){
								fprintf(fileit,"#LOGPERF!: %d\n",ii);
							} 
						}
						fflush(fileit);
					}
				}
				GETTIME(&gtimestart);
				tlasttime=t;
      }// end if output speed
    }

    if(PERFTEST){
			GETTIME(&checktime);
			walltime=(SFTYPE) DELTATIME(checktime,timestart);
			if(walltime<1E-5) walltime=1E-5;
			// setup so each turn is about the same WALLtime
      //itemp=(int)((SFTYPE)((SFTYPE)PERFWALLTIME/(walltime/(SFTYPE)nstep)));
			// setup so each turn is same as estimated time and speed, but fixed timesteps for all runs: best for benchmark if you know ahead of time ZCPSESTIMATE, and use same value of this across all tests
			itemp=(int)((SFTYPE)(PERFWALLTIME*ZCPSESTIMATE) /( (SFTYPE)realtotalzones) );
      //    fprintf(stdout,"itemp: %d PWT: %d ZCPS: %d realtotalzones: %d\n",itemp,PERFWALLTIME,ZCPSESTIMATE,realtotalzones);
      if(itemp<1) itemp=1;
      // set final time so steps to desired # of steps
      tf=1.0*(SFTYPE)(itemp)/(SFTYPE)(nstep)*t;
      if(myid==0){ fprintf(stderr,"PERFTEST: nstep=%d/%d t=%15.10g/%15.10g wt=%15.10g/%15.10g\n",nstep,itemp,t,tf,walltime,(SFTYPE)itemp*(SFTYPE)walltime/(SFTYPE)nstep); fflush(stderr);}
      
      //    if(itemp>1000) itemp=1000;
      if(nstep==itemp) reallaststep=1;
			if(walltime>PERFWALLTIME) reallaststep=1;
      // GODMARK(commented) 
      //if(nstep==100) reallaststep=1;
    }
    // GODMARK(commented)
    //    if(nstep==1) reallaststep=1;
    
    if(FLUSHFAILDT){
      fflush(fail_file);
      //    fflush(log_file);
    }

    //    if(nstep==1) exit(0);

  }// end over all time


  /////////////////  END COMPUTATION




  GETTIME(&timestop);
	microtime(&wttimestop);
	myustimes2(&usertmstimestop,&systmstimestop);

  if(DOGENDIAG){
    /* do final diagnostics */
    diag(2);
  }

  if(0&&DOAVGDIAG){// SUPERMARK
    if(avgdone==0){// finish average if didn't actually get past final average time
      diagavg(2);
      avgdone=1;
    }
  }

  if(goend==1){
    fprintf(log_file,"proc: %s Go end called(go.go)\n",myidtxt) ;
    fflush(log_file);
    if(myid<=0){
      fprintf(logfull_file,"#Go end called(go.go)\n") ;
      fflush(logfull_file);
    }
  }
  fprintf(log_file,"nstep: %ld\n",nstep) ;
  fprintf(log_file,"subnstep: %ld\n",subnstep) ;
  fflush(log_file);

  // running average zonecycle rate
  walltime=(SFTYPE) DELTATIME(timestop,timestart);
  if(walltime<1E-5) walltime=1E-5;
  
  if(myid<=0){
    fprintf(logfull_file,"#allproc: steps: %10ld wtime: %10.2g zcycles: %10d t: %10.2g\n",nstep,walltime*SEC2HOUR,(int)((FTYPE)(realtotalzones)*(FTYPE)nstep/walltime),(t-comptstart)) ;
    fprintf(logfull_file,"#(sec) walltime: %21.15g usertime: %21.15g systime: %21.15g\n",diffmicrotime(wttimestop,wttimestart),diffmyustimes(usertmstimestop,usertmstimestart),diffmyustimes(systmstimestop,systmstimestart));
		fprintf(stderr,"#(sec) walltime: %21.15g usertime: %21.15g systime: %21.15g\n",diffmicrotime(wttimestop,wttimestart),diffmyustimes(usertmstimestop,usertmstimestart),diffmyustimes(systmstimestop,systmstimestart));

    if(DOLOGPERF){
      fprintf(logperf_file,"#done: steps: %10ld wtime: %10.2g zcycles: %10d t: %10.2g tu/hour: %10.5g\n",nstep,walltime*SEC2HOUR,(int)((FTYPE)(realtotalzones)*(FTYPE)nstep/walltime),(t-comptstart),(t-comptstart)/(walltime*SEC2HOUR)) ;
      fflush(logperf_file);
    }
    if(PERFTEST){
      fprintf(perfout,"%10d\n",(int)((FTYPE)(realtotalzones)*(FTYPE)nstep/walltime)) ;
      fprintf(stderr,"perf: N3: %d N2: %d N1: %d RTZ: %d ZCPS: %d steps: %ld walltime: %15.10g\n",N3,N2,N1,realtotalzones,(int)((FTYPE)(realtotalzones)*(FTYPE)nstep/walltime),nstep,walltime) ;
      fflush(stderr);
      fclose(perfout);
    }
    if(DOLOGSTEP){
      fprintf(logstep_file,"#done: steps: %10ld wtime: %10.2g zcycles: %10d t: %10.2g\n",nstep,walltime*SEC2HOUR,(int)((FTYPE)(realtotalzones)*(FTYPE)nstep/walltime),(t-comptstart)) ;
			fflush(logstep_file);
    }
    fflush(logfull_file);
  }
  
  // end cut pasted from in loop

  if(DOFLOORDIAG>=1){
    fprintf(logfl_file,"floorcnt: t=%15.10g\n",t);
    for(i=0;i<NUMFLOOROUT;i++){
      for(j=1;j<=NUMFLOORVAR;j++){
				fprintf(logfl_file,"%12d ",floorcnt[i][j]);
      }
      fprintf(logfl_file,"\n");
    }
    for(j=1;j<=NUMFLOORVAR;j++){
      fprintf(logfl_file,"lowest %d %15.10g\n",wherelowest[j],floorlowest[j]);
    }
    fflush(logfl_file);
  }
    

  return(myexit(0));
}




void mycpuclock(clock_t *time)
{
  *time=clock();
}

void myustimes(clock_t *time) // returns number of microseconds
{
	struct tms mytimes;
	clock_t ret;
	long clockspersecond;

	clockspersecond=sysconf(_SC_CLK_TCK);
	ret=times(&mytimes);
	*time=(clock_t) (1000000.0*(SFTYPE)(mytimes.tms_utime+mytimes.tms_stime+mytimes.tms_cutime+mytimes.tms_cstime)/(SFTYPE)clockspersecond );
}
void myustimes2(clock_t *usertime,clock_t *systime) // returns number of microseconds
{
	struct tms mytimes;
	clock_t ret;
	long clockspersecond;

	clockspersecond=sysconf(_SC_CLK_TCK);
	ret=times(&mytimes);
	*usertime=(clock_t) (1000000.0*(SFTYPE)(mytimes.tms_utime+mytimes.tms_stime)/(SFTYPE)clockspersecond );
	*systime=(clock_t) (1000000.0*(SFTYPE)(mytimes.tms_cutime+mytimes.tms_cstime)/(SFTYPE)clockspersecond );
}

void myargs(int argc, char *argv[])
{
	if(GPARREAD==2){

//    fprintf(stderr,"argc: %d argv[0]: %s argv[1]: %s\n",argc,argv[0],argv[1]);
    if(argc!=4){
      fprintf(stderr,"incorrect command line: argc: %d\n",argc);
      exit(1);
    }
    ncpux1=atoi(argv[1]);
    ncpux2=atoi(argv[2]);
    ncpux3=atoi(argv[3]);
  }


}

