#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif


// stuff included in init.c that can be defined as a variable instead of constant in global.h
void gpar_init(void)
{
	FILE *in,*out;
	char temps[MAXFILENAME];


	if(POSTPROC==0){ // otherwise assume set by postproc.c
		if(!gocont){
			// MARK
			runtype=0; // choice
			directinput=0; // choice: only applies if runtype>0
			timereenter=2.0E4;
			// time to reenter calculation
			// time only used when directinput==3
		}
	}
	//
	// runtype(determines from what source data comes, or 0 implies no reentrance at all):
	// change numbers at botton of function too!
	// 0: nonreentrant
	// 1: data only reentrant
	// 2: par only reentrant
	// 22: gpar file but not grid file reentrant
	// 3: data and par reentrant
	// check for runtime parameter file, and if exists, use it.
	// Think binary, 0=scratch, 1=do from file
	// if 00=0 or no number(argc==1) -> par->scratch data->scratch
	// for rest, argc==2
	// if 01=1 -> par->scratch data->file
	// if 10=2 -> par->file    data->scratch
	// if 11=3 -> par->file    data->file
	// runtype is the decimal equiv of that binary code
	// 
	//
	// directinput(controls format of file that's inputted--only used if runtype>0 when pp=0):
	// 0: user uses .in files for all inputs(good for distinct start from other data set as initial data set) (currently, data: 0_pdump.in.xx when DUMPSM==1, 0_dump.in.xx when DUMPSM==0)
	// 1: use given numbers for dumps in init.c and normal .par ext for inputs, directly starting where left off when doing runtype>0
	// 2: Total reentrance based upon assumption that want last data outputted to be used (most useful for startup after crash/go.go stop, or continuing after run done but want to continue(extend time)
	// 3: like #1, but specify time instead of dump # (most useful for random startup)


	// only used if runtype==1 or runtype==0
	//
	// -1=test solution 1
	// 0=no analytic solution output
	// 1=sod solution(MHD Riemann problem too)
	// 2=advection test solution(includes coord=3 test)
	// 3=gaussian advection(setup for periodic boundary conditions)
	// 4=bondi solution(including magnetic WD solution)
	// 5=tori1
	// 6=visc check(must manually turn off all components except viscreal, and make all boundaries dirichlet, also turn off floor checks!). (NOT SETUP, MUST UNCOMMENT ANALOUTPUT==6 lines in step.c)
	// 7=pulsesol(general pulse type tests)
	// 8=inject: initial setup of just atmosphere
	// 9 = magnetic breaking solution
	// 10 = Orszag-Tang Vortex
	// 11 = Low1984 Multidimensional Solar Coronal Transient Model(not yet)
	// 12 = wave tests(including current sheet and rotated alfven waves)
	// 13 = Chandran turbulence
	// 14 = bondi plus torus
	// 15 = avery star
	// 16 = avery wind
        // 17=  layered disk 
	// INITIALDATASET
	analoutput=17;
  
	if(USEMPI==1){
		// set MPI output method
		if(USEMPI&&USEROMIO){
			mpicombine=1; // nochoice
		}
		else if(USEMPI&&USEJONIO){
			mpicombine=1; // no choice
		}
		else mpicombine=0; // nochoice
	}
	else{
		mpicombine=0; // no choice
	}  
   


	// 0: write as seperate files with CPU# as suffix
	// 1: get MPI to automagically write to 1 file (forced binary format)
	// (use bin2txt to convert to/from binary/text)
	if(mpicombine==1){
		binaryoutput=1; // no choice currently
	}
	else{
		// set output method
		binaryoutput=0; // choice
	}
  
	if(COMPDIM==3){
		cpugeompick=0; // no choice
	}
	else cpugeompick=0;  // choice
	// 0: computer pick geometry: see init_MPI, global problem size artibtrary
	// 1: user pick geometry: make sure NCPUX1*NCPUX2=numprocs(i.e. mpirun -np $numprocs), assumes fixed global problem size wanted
	// can run ./testmcpu.c which is same as init_optimalmpi() to check before running
	if(USEMPI){
		if(GPARREAD==1){
			// then read in from file during runtime
			sprintf(temps,"%s0_gpar%s",DATADIR,".in");
			if((in=fopen(temps,"r"))==NULL) {
				fprintf(fail_file,"error opening gpar input file %s\n",temps);
				myexit(1) ;
			}
			fscanf(in,"%d %d %d",&ncpux1,&ncpux2,&ncpux3);
			fclose(in);
		}
		if(GPARREAD==0){
			ncpux1=1;
			ncpux2=2;
			ncpux3=1;
		}

			// rest determined by function calls
	}
	else{
		// no choice...
		ncpux1=1;
		ncpux2=1;
		ncpux3=1;
		myid=0; // defines single process run
		sprintf(myidtxt,"");
		numprocs=1;
	}

  
	// make sure has trailing /
	if((POSTPROC==1)||(PRODUCTION==1) ) sprintf(DATADIR,"./"); // no choice
	else sprintf(DATADIR,"./"); // choice
  
	CHECKCONT=1;
	// 1: check if user wants to continue run 0: don't
  
	if(FLOATTYPE==1){
		// should be no smaller(and a bit greater) than machine precision!(i.e. more precisely, the round off of 1.0+DTOTLIMIT==1.0)
		// lower limit upon which symmetry point is taken so no symmetry violation from machine precision interference
		// really DTOTLIMIT/2 is size of region where sym-interp-fix occurs
		DTOTLIMIT=(1.0E-6);
		INVDTOTLIMIT=1./DTOTLIMIT;
		//#define DTOTLIMIT (1.0E-10) // 2.2204460492503131E-016 is round off epsilon(DBL_EPSILON)
	}
	else{
		// 1.192092896E-7F is round off epsilon(FLT_EPSILON)
		DTOTLIMIT=(1.0E-3);
		INVDTOTLIMIT=(1.0E+3);
	}
  
	// force flush of fail_file data each dt
	FLUSHFAILDT=1;
  
  
	// force floor on varibles
	// MARK
	FORCERHO=1; // 1: force DENSITYFLOOR on rho
	FORCEIE=1; // 1: force IEFLOOR on ie(s2 really)
	FORCERHOINTERNAL=1; // whether to force floor issue on each function each timestep (strictly better than above)
	FORCEIEINTERNAL=1; // whether to force floor issue on each function each timestep (strictly better than above)
  
	TRYSUBCYCLE=0; // whether to try to subcycle on viscosity or not
  
  
	CHECKDTLOW=1; // whether to check if dt went below critical value
	DTLOWEST=(1.E-12);
	IDTLOWEST=(1.0/DTLOWEST);
	SQIDTLOWEST=(1.0/(DTLOWEST*DTLOWEST)); // actually dt over cour then squared
	DTLOWDUMP=(DTLOWEST); // lowest dt to initiate dumping

  
	///////////////////
	// init diagnstics
	//
	// generally, should use sm to grid data correctly, but if don't care about reentering at all, can do here.
	if(mpicombine==1){
		DUMPSM=1; // choice
		FULLOUTPUT=FULLINPUT=0; // forced, no choice
	}
	else if(!GAMMIEIMAGE){
		DUMPSM=1; // 0: dump data where gridded 1: dump in center of zone for all
		FULLOUTPUT=0; // 0: output data for active grid, 1: half bzones out 2: full grid output
		// can only take fullinput on fulloutput data
		FULLINPUT=0; // 0: input data for active grid, 1: half bzones 2: full grid input
	}
	else{
		DUMPSM=1;
		FULLOUTPUT=0;
		FULLINPUT=0;
	}



	if(PERFTEST==0){ // can choose
		// 0: don't generate parameter output 1: do
		if(POSTPROC==0){
			DOGPARDIAG=1; // the flag-setup file
			DOPARDIAG=1; // generally do want to overwrite when pp=0 (the grid file)
		}
		else{
			DOGPARDIAG=1; // GODMARK // GODMARK
			DOPARDIAG=1; // generally don't want to overwrite diag when pp=1, but can.
		}
		// 0: don't generate interpolated parameter output 1: do
		if(POSTPROC==0){
			DOIPARDIAG=0; // should remain 0
		}
		else{
			DOIPARDIAG=1; // need 1 if interp data in pp
		}
    
		//////////
		// runtime diagnostics:

		appendold=0; // if option, whether to append loss/ener/etc files

		// if runtype>0, directinput==1, what file #'s to start at
		PDUMP_START=11;
		DUMP_START=55;
		FLDUMP_START=229;
		NPDUMP_START=55;
		ADUMP_START=55;
		FLOOR_START=5;
		IMAGE_START=458;

    
		// 0: don't output general diag 1: do output
		DOGENDIAG=1; // includes images, dumps, energy, etc.
		// 0: don't output loss diag 1: do output/counting(required in general)
    
    
		/////////////
		// LOSS STUFF
		DOLOSSDIAG=1;
		// 0: don't output avg diag 1: do output
		// if any size(N1/N2/N3) is 1 in a certain direction, I assume one doesn't care about loss then(see init_loss())
		//
		if(LOOPTYPE>1){
			COMPUTELOSSDIAG=1; // choice
			// 0: for adv use advection fluxes or whatever specified by _gen_adv function
			// 1: use _gen_gen function for adv fluxes
		}
		else{
			COMPUTELOSSDIAG=0; // no choice
		}
		//
		HYDROBZFLUX=0;
		// 0: let magnetic_flux do bz flux calc
		// 1: let hydro_flux do bz flux calc(not "working" right now)
		//
		// keep below at 1 or 2, not 0 for all version stuff to work
		DETAILMLOSS=1; // 0=only show total massloss each time
		// 1= show mass loss total and through each boundary(4 in 2D, 8 in 3D)
		// 2= Do 1 but in a seperate file output each grid points output each time (only for LOOPTYPE==1)
		//
		/////////////
    
    
		DOAVGDIAG=0;
		// 0: don't output floor diag 1: do output level 1 2: do full floor diag output (floor counters always active, this just controls diagnostics of counters)
		DOFLOORDIAG=1; // output details on what algorithm caused most lows and how low values got
		DOFLOORD2=0; // VERY DETAILED:  output where and what corrected in logfl_file
		//0: don't do dt diag 1: do
		DODTDIAG=1; // output of dt constraints
		DOTSTEPDIAG=0; // output of timescales(true global data of just above)
		if(COORD==3){
			DOSPDIAG=1; // sonic point location check log output
		}
		else{
			DOSPDIAG=0; // sonic point location check log output
		}
    
		DOLOGSTEP=1; // log the time step as per NDTCCHECK
		DOLOGPERF=1; // log the performance as per NZCCHECK
    
		DODIVBDIAG=1; // whether to do divb diagnostic
    
		// dumps
		if(DUMPSM==1){ // assume wanted, since need non-interp for restart
			PDUMPFLAG=1; // 0: don't dump primitive on own grid for reentrance 1: do
		}
		else PDUMPFLAG=0;
		DUMPFLAG=1; // 0: don't create dump files 1: do create
		if(VISCMEM==1){
			NPDUMPFLAG=0; // 0: don't create np dump files 1: do create (just pp now!)
		}
		else{
			NPDUMPFLAG=0; // no choice
		}
		FLOORDUMPFLAG=0; // 0: don't create floor dump files 1: do create
		if((N3==1)&&(N1>1)&&(N2>1)){ // then assume wanted
			FLDUMPFLAG=1; // 0: don't create fieldline dump files 1: do create
		}
		else{
			FLDUMPFLAG=0; // no choice
		}
		ADUMPFLAG=0; // 0: don't create analytic dump files 1: do create -1: only create first one
    
		if(!((N1==N1*N2*N3)||(N2==N1*N2*N3)||(N3==N1*N2*N3))){
			IMAGEFLAG=1; // choice: 0: don't create images 1: do create
		}
		else{
			IMAGEFLAG=0; // assume don't want image if 1D
		}
    
    
	}
	else if(PERFTEST==1){
    
		DOPARDIAG=0;
		DOGPARDIAG=0;
		DOIPARDIAG=0;
		DOGENDIAG=0;
		DOLOSSDIAG=0;
		DOAVGDIAG=0;
		DOFLOORDIAG=0; 
		DOFLOORD2=0;
		DODTDIAG=0;
		DOTSTEPDIAG=0;
		DOSPDIAG=0;
		DOLOGSTEP=0;
		DOLOGPERF=0;
    
		PDUMPFLAG=0; // 0: don't dump primitive on own grid for reentrance 1: do
		DUMPFLAG=0; // 0: don't create dump files 1: do create
		NPDUMPFLAG=0; // 0: don't create np dump files 1: do create (just pp now!)
		FLOORDUMPFLAG=0; // 0: don't create floor dump files 1: do create
		FLDUMPFLAG=0; // 0: don't create fieldline dump files 1: do create
		ADUMPFLAG=0; // 0: don't create analytic dump files 1: do create -1: only create first one
		IMAGEFLAG=0; // choice: 0: don't create images 1: do create
    
	}
  
	// whether to restart these diagnostics if restarting run
	restartener=1;
	restartloss=1;
	restartmode=0;
	restartsonicpoint=0;
  
  
	// how often in REAL *seconds* to dump 0_logstep.out file (unless 0, which then uses below)
	DTstep=10.0;
	DTstepdot=1.0;
	DTperf=DTstep;
	DTgocheck=30.0;
	DTtimecheck=60.0;
	NTIMECHECK=1000; // initial guess for how often to check time per timestep

	// how often in steps to output step/dt/t data (controlled by above if above are nonzero, else use below numbers)
	// MARK 100 100 20 500
	// MARK 10 10 1 100 for 1024x1024 vortex
	// MARK 1D bondi: 10000 10000 1000 20000
	// MARK 2D MHD Tori128128: 500 500 10 1000
	NDTCCHECK=100;
	// how often in steps to check speed in zonecycles/sec
	NZCCHECK=100;
	NDTDOTCCHECK=10;
	NGOCHECK=1000; // how often in steps to check the go.go file to see if to continue running
  
  
  
	OLDIMAGEFORMAT=0; // 0: newest 0-0-0-s1.r8 (tileable) 1: where 0-0-s1.r8 (no tile)
  
  
	OLDSCHOOL=0; // 0: assume as current data sets 1: use older format for data sets
  
	// whether header(0) or not(1) 
	// assume no header since consistent with MPI ROMIO output
	OLDSCHOOL2=1; // 0: assume as current r8's 1: as old r8's
  
	if(!GAMMIEIMAGE){
		sprintf(DUMPDIR,"");
		sprintf(IMAGEDIR,"i/"); // relative to DATADIR
		sprintf(IMAGEIDIR,"i/"); // interp relative dir
	}
	else{
		sprintf(DUMPDIR,"dumps/");
		sprintf(IMAGEDIR,"images/"); // relative to DATADIR
		sprintf(IMAGEIDIR,"iimages/"); // interp relative dir
	}
  
  
  
  
	// right now below 2 not used as all forced to be centered on a-grid
	// only applies to interpolation
	SAMEGRIDI=2; // 0=a-grid for sca b-grid for vec 1=all on a-grid 2=all on b-grid: Basis for where to interpolate to for image data
	SAMEGRIDD=2; // 0=a-grid for sca b-grid for vec 1=all on a-grid 2=all on b-grid: Basis for where to interpolate to for dump data
  
  
	if(POSTPROC==0){
		IMAGEFORMAT=0; // should always be 0 for post processing ability
		IMAGEFORMATINPUT=0; // dummy
	}
	else{
		if(GAMMIEIMAGE){
			IMAGEFORMAT=0;
		}
		else{
			IMAGEFORMAT=0; // choose final format for post proc(only interp)
		}
		IMAGEFORMATINPUT=0; // must be 0 
	}
	// 0: r8(best in general except for immediate viewing)
	// 1: ppm (best, esp. when used with gzip below)(can't be used during post process since can't reverse lookup easily, so don't use for now)

	if(mpicombine==1) IMAGEFORMAT=0; // forced for now
  
  
	if(POSTPROC==0){
    
		if(USEGM==1){
			GZIPIMAGE=0; // can't gzip with GM
		}
		else{
			GZIPIMAGE=3; // might as well always gzip otherwise
		}
		GZIPIMAGEINPUT=0; // dummy value since never input image normally
	}
	else{
		GZIPIMAGE=3; // best to always use this when post proc
		GZIPIMAGEINPUT=3; // choose was input type for post processing
	}
	// 0: don't gzip(best if need not zipped)(necessary for gm over mpich)
	// 1: do gzip using system call(best for compatibility)
	// 2: do gzip using fork call (best for small files)
	// 3: do gzip using popen (best for large files)
  
  
	// whether to bound potential(generally should be kept as analytic solution)
	NOBOUNDPOT=1; // 0: bound potential like other scalars, 1: no bound of pot
  
  
  
	// Make sure boundary zones do not allow inflow on outflow boundary(see bound.c)
	// only applies if outflow condition
	INFLOWCHECKIX1=1; // inner x1-edge
	INFLOWCHECKOX1=1; // outer x1-edge
	INFLOWCHECKIX2=1; // inner x2-edge // only use in 2D+
	INFLOWCHECKOX2=1; // outer x2-edge
	INFLOWCHECKIX3=1; // inner x3-edge // only use in 3D+
	INFLOWCHECKOX3=1; // outer x3-edge
  
  
	// file versions numbers(use sm for backwards compat)
	PVER= 11;
	GRIDVER= 3; // 3 is without cotangent
	DVER= 1 ;   // dumps same as for pdumps, adumps
	FLVER= 2;
	NPVER= 2;
	AVG1DVER= 2;
	AVG2DVER= 2;
	ENERVER= 7; // 6 is without c/s mode amp in ener, 7 has new ang mom
	MODEVER= 2; // 2 is all vars for 9 modes
	LOSSVER= 7; // 6 has x3 losses, 5 doesn't, 7 has new ang mom losses
	SPVER=   1;
	TSVER=   1;
	LOGDTVER= 1;
	STEPVER= 1;
	PERFVER= 3;
	ADVER= DVER;
	PDVER= DVER;
	CALCVER= 1;
	FLINEVER= 1;
	// type designations for sm automagical read in correct format for similar things
	PTYPE=     1; // global par file
	GRIDTYPE=  2;
	DTYPE=     3 ;// dump
	FLTYPE=    4; // floor
	NPTYPE=    5; // np
	AVG2DTYPE= 6;
	AVG1DTYPE= 7;
	ENERTYPE=  8;
	LOSSTYPE=  9;
	SPTYPE=    10;
	TSTYPE=    11;
	LOGDTTYPE= 12 ;
	STEPTYPE=  13;
	PERFTYPE=  14;
	ADTYPE=    15 ;// same as dump except filename
	PDTYPE=    16; // same as dump except filename
	CALCTYPE=  17; // arbitrary calcs during pp
	FLINETYPE=  18; // field line during pp
	MODETYPE=  19;
	EXPANDTYPE= 50 ;// used to signify doing pp expansion
	NPCOMPUTETYPE= 33; // used to signify want to compute np before output
  


// extention for data files
	strcpy(DATEXT,".dat");
	strcpy(PAREXT,".par");
	strcpy(INEXT,".in");
  
	if(binaryoutput){ // binaryoutput forced if mpicombine==1
		strcat(DATEXT,".bin");
		strcat(PAREXT,".bin");
		strcat(INEXT,".bin");
	}
	else{
		// just leave normal
		/*
			strcat(DATEXT,".txt");
			strcat(PAREXT,".txt");
			strcat(INEXT,".txt");
		*/
	}
	strcpy(PPEXT,".pp"); // always text
	strcpy(OUTEXT,".out"); // always text
	strcpy(IMGEXT,".dat"); // always binary
	strcpy(DAT2EXT,".dat"); // always text
	strcpy(PAR2EXT,".par"); // always text

// number of wallseconds per perf run(see: main.c)
	PERFWALLTIME=30.0;
// 1 linux cluster cpu
// ZCPSESTIMATE (50000)
// 25 linux cluster cpu(550Mhz Xeon's connected via Myrinet)
// ZCPSESTIMATE (1250000)
// 36 linux cluster cpu(550Mhz Xeon's connected via Myrinet)
// ZCPSESTIMATE (1800000)
// 64 linux cluster cpu(550Mhz Xeon's connected via Myrinet)
// ZCPSESTIMATE (3200000)
// 121 linux cluster cpu(550Mhz Xeon's connected via Myrinet)
// ZCPSESTIMATE (6050000)
// photon MHD for one cpu alone
// ZCPSESTIMATE (265000)
// photon HD 1 cpu
// ZCPSESTIMATE (400000)
// rainman MHD for one cpu alone
// ZCPSESTIMATE (220000)
// ZCPSESTIMATE (100000)
// sgi r10000 for one cpu alone(195Mhz)
// ZCPSESTIMATE (80000)
// 4cpu mpigm
// ZCPSESTIMATE (800000)
// 4cpu r10000
// ZCPSESTIMATE (343000)
// 9cpu r10000
// ZCPSESTIMATE (745000)
// 16cpu r10000
// ZCPSESTIMATE (1325000)
// 25cpu r10000
// ZCPSESTIMATE (1200000)
// 36cpu r10000
// ZCPSESTIMATE (1700000)
// 49cpu r10000
// ZCPSESTIMATE (4021000)
// 64 r10000's 64x64 tile
// ZCPSESTIMATE (4309000)
// 121 r10000's
// ZCPSESTIMATE (10943000)
// 256 r10000's
// ZCPSESTIMATE (20000000)
// kerr 2DMHD
// ZCPSESTIMATE (50000)
// rainman 2DMHD
	ZCPSESTIMATE=(200000);



	// then read in from file during runtime
	//sprintf(temps,"%s0_gpar%s",DATADIR,".out");
	//if((in=fopen(temps,"w"))==NULL) {
//		fprintf(fail_file,"error opening gpar output file %s\n",temps);
//		myexit(1) ;
//	}

}
