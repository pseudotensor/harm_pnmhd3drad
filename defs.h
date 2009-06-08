//unit conversion
SFTYPE munit,lunit,rhounit,tunit,eneunit,enedenunit;
//Constant
SFTYPE lsuncgs,kcontcgs,mhydrcgs,sigmacontcgs,sunmasscgs,yrcgs,mstarcgs;
SFTYPE sigmacontnew,kcontnew,mhydrnew;
SFTYPE rsun,mmw;
//Central star
SFTYPE lsunfrac,lstarcgs,rstar;
// layer model
SFTYPE Siga,cstcrit2,alphafrac;
FTYPE iab[N3M][N2M][N1M];
int actcontour;
//Temperature,opaini
FTYPE tirr0a[N3M][N2M][N1M],opainia[N3M][N2M][N1M],siginia[N3M][N2M][N1M],optdinia[N3M][N2M][N1M];
FTYPE cs02a[N3M][N2M][N1M];
FTYPE (*tirr0)[N2M][N1M],(*opaini)[N2M][N1M],(*sigini)[N2M][N1M],(*optdini)[N2M][N1M];
FTYPE (*cs02)[N2M][N1M];
FTYPE teff4[N1];
// make sure all this matches init.c: init_pointers()


/* 1: a: + dx in 1/2/3-direction from cell face: x(i)=x(i-1)+dx(i-1) */
/* 2: b: + dx in 1/2/3-direction from cell center: x(i)=x(i-1)+dx(i) */
FTYPE dxa[NUMGRID][3][NBIGM]; // [agrid=1,bgrid=2][dir=1,2,3][i=-NBIGBND..NBIG+NBIGBND]
FTYPE (*dx)[3][NBIGM] ;

/* Positions and sizes of cell */
/* cell face position 1: a */
/* cell center position 2: b*/
FTYPE xa[NUMGRID][3][NBIGM]; // [agrid=1,bgrid=2][dir=1,2,3][i=-NBIGBND..NBIG+NBIGBND]
FTYPE (*x)[3][NBIGM] ;

// layered 
FILE *grid_file;
FTYPE alpha_realar[N3M][N2M][N1M];
FTYPE(*alpha_reala)[N2M][N1M];

// only interp if pp
#if(POSTPROC==1)
FTYPE idxa[NUMGRID][3][INTNBIG];
FTYPE ixa[NUMGRID][3][INTNBIG];
FTYPE ppcalcsa[NUMPPCALCS][N3M][N2M][N1M] ;
#endif
// keep below so code is intact
FTYPE (*idx)[3][INTNBIG] ;
FTYPE (*ix)[3][INTNBIG] ;
FTYPE (*ppcalcs)[N3M][N2M][N1M];

/* at cell center 1: d 2: e 3: phi */
FTYPE sa[NUMSCA][N3M][N2M][N1M] ; // [rho=1,en=2,pot=3][k=-N3BND..N3+B3BND][j=-N2BND..N2+B2BND][i=-N1BND..N1+B1BND]
FTYPE (*s)[N3M][N2M][N1M] ; /* 3 scalars */

/* at cell face, where position of cell-face var is less than cell center var at same index */
/* 1: v 2: B*/
FTYPE va[NUMVEC][3][N3M][N2M][N1M] ; //[velocity=1,B-field=2][1=x1dir,2=x2dir,3=x3dir][k=-N3BND..N3+B3BND][j=-N2BND..N2+B2BND][i=-N1BND..N1+B1BND]
FTYPE (*v)[3][N3M][N2M][N1M] ;

////// analytic storage
/* at cell center 1: d 2: e 3: phi */
FTYPE sanala[NUMSCA][N3M][N2M][N1M] ; // same as s
FTYPE (*sanal)[N3M][N2M][N1M] ;

#if(GRAVACC&&(!POSTPROC))
FTYPE gravacca[3][N3M][N2M][N1M] ; // analytic gravitational accel for each direction: assumes only 1 direction currently
#endif
FTYPE (*gravacc)[N3M][N2M][N1M] ;

#if(GRAVITOMAGNETIC&&(!POSTPROC))
FTYPE gravitoma[2][N3M][N2M][N1M] ; // geometry terms in gravitomagnetic source term
#endif
FTYPE (*gravitom)[N3M][N2M][N1M] ;

/* at cell face, where position of cell-face var is less than cell center var at same index */
/* 1: v 2: B*/
FTYPE vanala[NUMVEC][3][N3M][N2M][N1M] ; // same as v
FTYPE (*vanal)[3][N3M][N2M][N1M] ;

#if(DOINVSOL2FUN)
FTYPE invsol2funa[N3M][N2M][N1M] ; // 1/c^2 as a function of something
#endif
FTYPE (*invsol2fun)[N2M][N1M] ;

#if(DOAVGDIAGMEMORY)
// 2-D, averaged over 3
/* at original cell location */
// [rho=1,en=2,Be=3,cs^2=4,S=5,|vx1|=6,|vx2|=7,|vx3|=8]
//[j=-N2BND..N2+B2BND][i=-N1BND..N1+B1BND]
FTYPE avg2_3a[NUMAVG2_3][N2M][N1M] ;

// 1-D, averaged over 3 and 1 direction
/* at original cell location */
// [rho=1,en=2,Be=3,cs^2=4,S=5,|vx1|=6,|vx2|=7,|vx3|=8]
//[j=-N2BND..N2+B2BND]
FTYPE avg1_31a[NUMAVG1_31][N2M] ;

// 1-D, averaged over 3 and 2 direction
/* at original cell location */
// [rho=1,en=2,Be=3,cs^2=4,S=5,|vx1|=6,|vx2|=7,|vx3|=8]
//[i=-N1BND..N1+B1BND]
FTYPE avg1_32a[NUMAVG1_32][N1M] ;
#endif
FTYPE (*avg2_3)[N2M][N1M] ;
FTYPE (*avg1_31)[N2M] ;
FTYPE (*avg1_32)[N1M] ;

#if(FLOORDUMPFLAGMEMORY)
// used for floor on grid for scalar values.
// 1: rho->rho
// 2: u->u
// 3: phi->rho*phi (gravitational potential energy density)
FTYPE floorvarsa[NUMSCA][N3M][N2M][N1M] ;

// used for floor on grid for vectors(note: need 0 element here)
// all due to rho changes
// 1: momentum density in radial direction
// 2: theta-angular momentum density
// 3: phi-angular momentum density
FTYPE floorvarva[NUMVEC][3][N3M][N2M][N1M] ;

// 0: KE
// 1: BE
FTYPE floorvar0a[NUMVEC][N3M][N2M][N1M] ;

#endif
FTYPE (*floorvars)[N3M][N2M][N1M] ;
FTYPE (*floorvarv)[3][N3M][N2M][N1M] ;
FTYPE (*floorvar0)[N3M][N2M][N1M] ;

#if((LOWMEMMODE==0)&&(TVDLF==0))
// comp savers (first[] is location of variable used on)
// direction specifies direction of surface normal vector
FTYPE dsa[4][3][N3M][N2M][N1M]; // [1=(.5,5),2=(0,.5),3=(.5,0),4=(0,0)][dir=1,2,3][k=-N3BND..N3+N3BND-1][j=-N2BND..N2+N2BND-1][i=-N1BND..N1+N1BND-1]

FTYPE oarcla[4][3][N3M][N2M][N1M]; // [1=(.5,5),2=(0,.5),3=(.5,0),4=(0,0),5=(.5,.5,0)][dir=1,2,3][k=-N3BND..N3+N3BND-1][j=-N2BND..N2+N2BND-1][i=-N1BND..N1+N1BND-1]

// just dx term, not arclength
FTYPE odxa[NUMGRID][3][NBIGM]; // [agrid=1,bgrid=2][dir=1,2,3][i=-NBIGBND..NBIG+NBIGBND]

FTYPE ovola[4][N3M][N2M][N1M]; // [1=(.5,5),2=(0,.5),3=(.5,0),4=(0,0)][k=-N3BND..N3+N3BND-1][j=-N2BND..N2+N2BND-1][i=-N1BND..N1+N1BND-1]
#endif


  /* metric scales and their dervs--only applies to limited coord systems */
  /* 1: a-mesh */
  /* 2: b-mesh */
  /* 1: g1 2: g2 3: g31 4: g32 */
#if(LOWMEMMODE==0)
FTYPE ga[NUMGRID][NUMMETRIC][NBIGM]; // [1,2][1,2,3,4][...]
FTYPE dga[NUMGRID][NUMMETRIC][NBIGM];  // [1,2][1,2,3,4][...]
FTYPE oga[NUMGRID][NUMMETRIC][NBIGM]; // inverse of g

  /* Volume factors: dvl[a/b][dir][location] */
FTYPE dvla[NUMGRID][3][NBIGM]; //[1,2][1,2,3][...]
FTYPE odvla[NUMGRID][3][NBIGM]; //[1,2][1,2,3][...]
#endif

#if((PUREBC==0)&&(BOUNDTYPE==1))
short bcsa[NUMSCA][3][N3M][N2M][N1M]; // [which scalar][1=type:assign=(1,2,3,4,5),2=dim:assign=(1,2,3),3=direction:assign=(-1,+1)][...][...][...]
//bcv assumes all components bounded with same *type* of boundary condition
short bcva[NUMVEC][3][N3M][N2M][N1M]; // [which vector][1=type:assign=(1,2,3,4,5),2=dim:assign=(1,2,3),3=direction:assign=(-1,+1)][...][...][...]

#elif((BOUNDTYPE==2)||(BOUNDTYPE==3))

short bzmaska[N3M][N2M][N1M];
short bzmaskina[N3M][N2M][N1M];
// 0: comp zone
// 1: normal bzone
// -1: no-where zone(never used)
// 99: MPI zone
// 3 directions, 8 maximal values to average
// specifies location of zone to copy from

#endif

short (*bcv)[3][N3M][N2M][N1M];
short (*bcs)[3][N3M][N2M][N1M];

short (*bzmask)[N2M][N1M];
short (*bzmaskin)[N2M][N1M];
short (***bzs); // worst case scenario is 26 positions per zone to average
// pointer to pointer for an array of pointers

FTYPE (*ds)[3][N3M][N2M][N1M] ;
FTYPE (*oarcl)[3][N3M][N2M][N1M] ;
FTYPE (*odx)[3][NBIGM] ;
FTYPE (*ovol)[N3M][N2M][N1M] ;
FTYPE (*g)[NUMMETRIC][NBIGM];
FTYPE (*og)[NUMMETRIC][NBIGM];
FTYPE (*dg)[NUMMETRIC][NBIGM];
FTYPE (*dvl)[3][NBIGM];
FTYPE (*odvl)[3][NBIGM];


  /* type=1-4 are local, type=5 is nonlocal copy */
  /* see bound.c for types etc. */
// -1: corner
// 0: comp
// 1: reflect
// 2: AOS
// 3: Dirichlet
// 4: outflow
// 5: periodic
// 99: MPI(dirichlet from other cpus)

// 0: bzone loop
// 1: LOOPFC
// 2: LOOPHC
// 3: LOOPFMHP
// 4: LOOPHMFP
// 5: LOOPC/etc
// 6: LOOPV1
// 7: LOOPV2
// 8: LOOPV3

#if((BOUNDTYPE==3)&&(POSTPROC==0))
short iindxa[NUMINDEX][2][N3M][N2M];
#endif
short (*iindx)[2][N3M][N2M];



#if(VISCMEM)
FTYPE nu_facta[N3M][N2M][N1M] ;
FTYPE nu_reala[N3M][N2M][N1M] ;
#endif
FTYPE (*nu_fact)[N2M][N1M] ;
FTYPE (*nu_real)[N2M][N1M] ;


#if(RESMEM&&(POSTPROC==0))
FTYPE nu_res_facta[N3M][N2M][N1M] ;
FTYPE nu_res_reala[N3M][N2M][N1M] ;
FTYPE jcurrenta[3+1][N3M][N2M][N1M] ;
#endif
FTYPE (*nu_res_fact)[N2M][N1M] ;
FTYPE (*nu_res_real)[N2M][N1M] ;
FTYPE (*jcurrent)[N3M][N2M][N1M] ;


#if((LOOPTYPE==1)&&(DOLOSSDIAGMEMORY&&(POSTPROC==0)))
// only over active grid
SFTYPE losssa[NUMSCA][3][2][NBIG][NBIG]; // [which=1...][dir=1,2,3][i/o=0,1][list=0...]
SFTYPE lossva[NUMVEC][3+1][3][2][NBIG][NBIG]; // [which][v-dir][dir][i/o][list]
SFTYPE lossvisca[1][3][2][NBIG][NBIG]; // [which][dir][i/o][list]
#endif
SFTYPE (*losss)[3][2][NBIG][NBIG];
SFTYPE (*lossv)[3+1][3][2][NBIG][NBIG];
SFTYPE (*lossvisc)[3][2][NBIG][NBIG];
// v-dir: vector direction // 0=energy, 1,2,3
// dir: 1=x1 2=x2 3=x3
// i/o: 0=inner boundary 1=outer boundary
// list: var along that side: on active grid: 0..N-1(i.e. don't shift pointer!)


short accountstorea[N3M][N2M][N1M];
short (*accountstore)[N2M][N1M];

#if(MDOTMEM&&(POSTPROC==0))
FTYPE rhoinjecta[N3M][N2M][N1M];
FTYPE eninjecta[N3M][N2M][N1M];
#endif
FTYPE (*rhoinject)[N2M][N1M];
FTYPE (*eninject)[N2M][N1M];

// used to hold boundary data(packed)
#if((TVDLF==0)&&(POSTPROC==0))
FTYPE worksbca[NUMSCA][2][COMPDIM*2][NBIGBND*NBIGSM]; //scalar: [scalar#][1=out/2=in][0=right,1=up,2=left,3=down,4=out,5=in][datawidth]
FTYPE workvbca[NUMVEC][2][COMPDIM*2][3*NBIGBND*NBIGSM]; //vector: [vector#][as above][as above][datawidth]
#endif
FTYPE (*worksbc)[2][COMPDIM*2][NBIGBND*NBIGSM];
FTYPE (*workvbc)[2][COMPDIM*2][3*NBIGBND*NBIGSM];




// workXa is just like s[i] (i.e. one particular scalar)
FTYPE work1a[N3M][N2M][N1M];
FTYPE work2a[N3M][N2M][N1M];
FTYPE work3a[N3M][N2M][N1M];
#if(TVDLF==0)
FTYPE work4a[N3M][N2M][N1M];
FTYPE work5a[N3M][N2M][N1M];
FTYPE work6a[N3M][N2M][N1M];
FTYPE work7a[N3M][N2M][N1M];
FTYPE work8a[N3M][N2M][N1M];
FTYPE work9a[N3M][N2M][N1M];
FTYPE work10a[N3M][N2M][N1M];
#endif
FTYPE (*work1)[N2M][N1M];
FTYPE (*work2)[N2M][N1M];
FTYPE (*work3)[N2M][N1M];
FTYPE (*work4)[N2M][N1M];
FTYPE (*work5)[N2M][N1M];
FTYPE (*work6)[N2M][N1M];
FTYPE (*work7)[N2M][N1M];
FTYPE (*work8)[N2M][N1M];
FTYPE (*work9)[N2M][N1M];
FTYPE (*work10)[N2M][N1M];

#if(MDOTMEMANAL&&(POSTPROC==0))
FTYPE mdotanala[N3M][N2M][N1M];
#endif
FTYPE (*mdotanal)[N2M][N1M];

FTYPE osqrtrhoa[N3M][N2M][N1M];
FTYPE (*osqrtrho)[N2M][N1M];

#if(ALFVENLIMIT)
FTYPE rholimiteda[3+1][N3M][N2M][N1M];
#endif
FTYPE (*rholimited)[N3M][N2M][N1M];

#if(TVDLF==0)
// workvXa is just like v[i] (i.e. like one particular scalar)
FTYPE workv1a[3][N3M][N2M][N1M];
FTYPE workv2a[3][N3M][N2M][N1M];
FTYPE workv3a[3][N3M][N2M][N1M];
FTYPE workv4a[3][N3M][N2M][N1M];
FTYPE workv5a[3][N3M][N2M][N1M];
#endif
FTYPE (*workv1)[N3M][N2M][N1M];
FTYPE (*workv2)[N3M][N2M][N1M];
FTYPE (*workv3)[N3M][N2M][N1M];
FTYPE (*workv4)[N3M][N2M][N1M];
FTYPE (*workv5)[N3M][N2M][N1M];

// [1,2,3][1,2,3][...][...][...]
// below 3 tensors only used if real viscosity is present
#if(VISCMEM)
FTYPE sigmaa[3][3][N3M][N2M][N1M];
FTYPE rosta[3][3][N3M][N2M][N1M];
FTYPE rostnua[3][3][N3M][N2M][N1M];
#endif
FTYPE (*sigma)[3][N3M][N2M][N1M]; // -2*rho*nu*e_{ij}= sigma_{ij}
FTYPE (*rost)[3][N3M][N2M][N1M];// e_{ij}
FTYPE (*rostnu)[3][N3M][N2M][N1M]; // e_{ij}*nu

// workiqa is like s but in 2D and only on active grid [0..n2-1][0..n1-1]
// only postproc interpolates or images averages or kinetic energy floor
#if(POSTPROC==1)
FTYPE workiqa[NUMSCA][INTN2][INTN1];
// workviqa is like v but in 2D and only on active grid
FTYPE workviqa[NUMVEC][3][INTN2][INTN1];
FTYPE workiqavga[NUMWORKIQ][INTN2][INTN1];
FTYPE work0iqa[NUMVEC][INTN2][INTN1];
#endif
// still allocate pointers so code is intact except init.c pointer stuff
FTYPE (*workiq)[INTN2][INTN1];
FTYPE (*workviq)[3][INTN2][INTN1];
FTYPE (*workiqavg)[INTN2][INTN1];
FTYPE (*work0iq)[INTN2][INTN1];


// to speed up things use matrix for const grids
// assumes never interpolate when postproc=0
#if(POSTPROC==1)
FTYPE convmata[3][3][INTN2][INTN1];  
#endif
FTYPE (*convmat)[3][INTN2][INTN1];



#if((TVDLF==1)&&(POSTPROC==0))
/* primitive variables */
FTYPE a_p[N1+2][N2+2][NP] ;
FTYPE a_dqx[N1+2][N2+2][NP] ;
FTYPE a_dqy[N1+2][N2+2][NP] ;
FTYPE a_Fx[N1+2][N2+2][NP] ;
FTYPE a_Fy[N1+2][N2+2][NP] ;
FTYPE a_ph[N1+2][N2+2][NP] ;
#endif
FTYPE (*  p)[N2+2][NP] ;
FTYPE (*dqx)[N2+2][NP] ;
FTYPE (*dqy)[N2+2][NP] ;
FTYPE (* Fx)[N2+2][NP] ;
FTYPE (* Fy)[N2+2][NP] ;
FTYPE (* ph)[N2+2][NP] ;


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// EVERYTHING BELOW HERE NEEDS NO POINTER SHIFTING OR MEMORY ALLOCATION OR DEPENDS ON SIZE OF SYSTEM
//
//

short *indx[NUMINDEX]; // array of pointers
int numiter[NUMINDEX]; // 0-NUMINDEX-1


// counts at different 5 different places in code(see step.c and sweep.c)
// these 2 vars are only sources for floor changes in above arrays
int floorcnt[NUMFLOOROUT+1][NUMFLOORVAR+1]; //[0,1,2,3,4,5,6][1=rho,2=en]
// [0][] : step_visc
// [1][] : step_real_visc
// [2][] : step_ie
// [3][] : step_res
// [4][] : sweepx
// [5][] : sweepy
// [6][] : sweepz


// stores lowest value of rho and en before floor fixed, over entire run
FTYPE floorlowest[NUMFLOORVAR+1]; //[1=rho,2=en]

// what routine has the above lowest value generated
// same routine numbers as above
int wherelowest[NUMFLOORVAR+1]; //[1=rho,2=en]


//int N[3+1] ; /* Number of cells in each direction */
SFTYPE L[2+1][3+1] ; /* L[1] holds position of 1,1,1 cell outer edges, L[2] holds size of grid in units */
SFTYPE iL[2+1][3+1];

// outerdefs for images
SFTYPE outerdefs[ITYPES][CTYPES][NUMSCA+1];
SFTYPE outerdefv[ITYPES][CTYPES][NUMVEC+1][3+1]; // [][0] here is magnitude

// outerdef for dumps for interpolation
SFTYPE douterdefs[NUMSCA+1];
SFTYPE douterdefv[NUMVEC+1][3+1]; // [][0] here is magnitude
SFTYPE douterdef0[NUMVEC+1]; // KE for floor(not set yet)


SFTYPE IOBound[NUMBTRANS]; // used to have inflow/outflow transitions on outer theta edge

// average stuff needed elsewhere too
SFTYPE tavgstart,tavgfinal; // real start and final
int avgcount; //number of averages in time
int num1d_31_full,num1d_32_full;

SFTYPE floors[NUMLOSSVAR+1]; // to keep track of floor input of stuff
SFTYPE inflows[NUMLOSSVAR+1]; // to keep track of injection of stuff
SFTYPE radiations[NUMLOSSVAR+1]; // to keep track of radiation of stuff
// 0: etot
// 1: mass
// 2: enthalpy
// 3: grav pot energy
 
// 4: ke 
// 5: s1 (m*v_r in COORD==3)
// 6: s2 (r*m*v_theta in COORD==3)
// 7: s3 (r*sin(theta)*m*v_phi in COORD==3)

// 8: mag energy
// 9: B1
// 10: B2
// 11: B3
 
// 12: etot(visc)

SFTYPE invsol2,invsol2_original; // inverse speed of light, squared, for alfven limit trick

// next block of changable global parameters(global.h just sets initial value)
int kever;
int mdotin,cool,press,pressx1,pressx2,pressx3,visc_art,ie,mag,pmag,res_real,visc_real,trans;
int stepmocct,mocctvx1,mocctvx2,mocctvx3,mocctbx1,mocctbx2,mocctbx3,transbzy,transbzx;
int vischeat;
int resheat;
int vreal;
int rreal;
int transx1,transx2,transx3;
int transiex1,transmagx1,transv1x1,transv2x1,transv3x1,transrhox1;
int transiex2,transmagx2,transv1x2,transv2x2,transv3x2,transrhox2;
int transiex3,transv1x3,transv2x3,transv3x3,transrhox3;
int dopassive;
int transpassivex1,transpassivex2,transpassivex3;
int transpassivev1x1,transpassivev1x2,transpassivev1x3,transpassivev2x1,transpassivev2x2,transpassivev2x3,transpassivev3x1,transpassivev3x2,transpassivev3x3;
int advint;

int ireenter;

SFTYPE coolbottom; // bottom on cooling in theta
SFTYPE rho0,rho00,u0,R0,Omega0; // for tori problem and for SPB rho visc
SFTYPE coolfact;
char WRITETYPE[10];
int reallaststep;
int tagik,tagfk,tagij,tagfj,tagii,tagfi;
int t2gik,t2gfk,t2gij,t2gfj,t2gii,t2gfi;
int t3gik,t3gfk,t3gij,t3gfj,t3gii,t3gfi;
int t4gik,t4gfk,t4gij,t4gfj,t4gii,t4gfi;
char myidtxt[MAXFILENAME];
int totalzones,realtotalzones;
int rtotalzones;
int itotalzones;
int sizes[3+1][MAXCPUS];
int isizes[3+1][MAXCPUS];
int totalsize[3+1];
int itotalsize[3+1];
int mycpupos[3+1]; // my position amongst the cpus
int srdir[6]; // which direction this cpu sends/receives normal interior data 
int startpos[3+1];
int endpos[3+1]; // startj and endj are where this CPU located on full grid
long nstep;
int runtype,directinput,appendold,deleteolddat,deleteoldpar,deleteoldimg;
int PDUMP_START,DUMP_START,NPDUMP_START,ADUMP_START,FLOOR_START,FLDUMP_START,IMAGE_START;
int simplebc,bcix1,bcox1,bcix2,bcox2,bcix3,bcox3;
int nonunigridx1,nonunigridx2,nonunigridx3;
int analoutput;
SFTYPE x1in,x1out,x2in,x2out,x3in,x3out;
SFTYPE ix1in,ix1out,ix2in,ix2out,ix3in,ix3out;
SFTYPE timereenter;
int pdump_start,dump_start,fldump_start,npdump_start,adump_start,floor_start,image_start;
FILE* fail_file;
FILE* logfull_file;
FILE* log_file;
FILE* logstep_file;
FILE* logperf_file;
FILE* logdt_file;
FILE* logfl_file;
FILE* logsp_file;
// Equation of state coefficients
SFTYPE tscycleto; // time to subcycle to
SFTYPE tscyclefrom; // time to subcycle from
SFTYPE dtlastscycle; // dt used as basis for subcycle
int subcyclen; // number of subycycles for viscosity
int nthsubcycle; // number of subcycles so far
SFTYPE rg,rgp;
SFTYPE cour,invcour ; // invcour used in timestep1.h
SFTYPE cour2,invcour2 ; // invcour2 used in timestep1.h
SFTYPE cs ;
SFTYPE alpha;
SFTYPE gam ;
int wgam;
int wgam1;
int wgam53;
int wpw;
int periodicx1,periodicx2,periodicx3;
int mpiperiodicx1,mpiperiodicx2,mpiperiodicx3;
int periodicx2special;
int skipix1,reflectix1,reflectox1;
int skipix2,reflectix2,reflectox2;
int skipix3,reflectix3,reflectox3;
int intix1,intox1,intix2,intox2,intix3,intox3;
int skipintix1,skipintix2,skipintix3;
SFTYPE dt ;
SFTYPE tstart,TSTART;
SFTYPE t,tf,timagescale ;
SFTYPE tavgi,tavgf;
int numavg;
SFTYPE nu_vnr ;
SFTYPE nu_l ;
SFTYPE nu_ten ;
SFTYPE alpha_real ;
SFTYPE alpha_real0 ;
SFTYPE n_real ;
SFTYPE GRAVC;
SFTYPE MASSBH ;
SFTYPE GM;
SFTYPE DTd ;
SFTYPE DTfld ;
SFTYPE DTpd ;
SFTYPE DTavg ;
SFTYPE DTi ;
SFTYPE DTl ;
SFTYPE DTener;
SFTYPE DTloss;
SFTYPE DTmode;
SFTYPE DTfloor;
SFTYPE DTtimestep;
SFTYPE DTdivb;
SFTYPE DTtimescale;
SFTYPE DTsp;
SFTYPE DTstep,DTstepdot,DTperf,DTgocheck,DTtimecheck;
SFTYPE resist_real0 ;
SFTYPE resist_real ;
SFTYPE nu_sh ;
//int numbc[3+1];
SFTYPE startx[3+1]; // starting location for each cpu for each direction
SFTYPE startix[3+1]; // starting location for each cpu for each direction
/* grid velocity */
SFTYPE vg[3+1];
SFTYPE RHO0,V0; // for advection
SFTYPE DENSITYFLOOR,IEFLOOR,DENSITYFRACFLOOR,IEFRACFLOOR;
SFTYPE IEFRACT,VZFRACT;
SFTYPE massdot;

SFTYPE divbmax_full,divbavg_full;

int DYNAMICMM;

SFTYPE blackholejz;

int mpicombine;
int binaryoutput,cpugeompick,CHECKCONT,FLUSHFAILDT,FORCERHO,FORCEIE,FORCERHOINTERNAL,FORCEIEINTERNAL,TRYSUBCYCLE,CHECKDTLOW;
SFTYPE DTOTLIMIT,INVDTOTLIMIT,DTLOWEST,IDTLOWEST,SQIDTLOWEST,DTLOWDUMP;
int DOGPARDIAG,DOPARDIAG,DOIPARDIAG,DOGENDIAG,DOLOSSDIAG,COMPUTELOSSDIAG,HYDROBZFLUX,DETAILMLOSS,DOAVGDIAG,DOFLOORDIAG,DOFLOORD2,DODTDIAG,DOTSTEPDIAG,DOSPDIAG,DOLOGSTEP,DOLOGPERF;

int DODIVBDIAG;

int NDTCCHECK,NZCCHECK,NDTDOTCCHECK,NGOCHECK,NTIMECHECK;


int PDUMPFLAG,DUMPFLAG,NPDUMPFLAG,FLOORDUMPFLAG,FLDUMPFLAG,ADUMPFLAG,IMAGEFLAG;
int OLDIMAGEFORMAT,OLDSCHOOL,OLDSCHOOL2;

char DATADIR[MAXFILENAME],DUMPDIR[MAXFILENAME],IMAGEDIR[MAXFILENAME],IMAGEIDIR[MAXFILENAME];

int DUMPSM,FULLOUTPUT,FULLINPUT;

int SAMEGRIDI,SAMEGRIDD;
int IMAGEFORMAT,IMAGEFORMATINPUT,GZIPIMAGE,GZIPIMAGEINPUT;

int NOBOUNDPOT;

int INFLOWCHECKIX1,INFLOWCHECKOX1,INFLOWCHECKIX2,INFLOWCHECKOX2,INFLOWCHECKIX3,INFLOWCHECKOX3;

int PVER,GRIDVER,DVER,FLVER,NPVER,AVG1DVER,AVG2DVER,ENERVER,MODEVER,LOSSVER,SPVER,TSVER,LOGDTVER,STEPVER,PERFVER,ADVER,PDVER,CALCVER,FLINEVER;

int PTYPE,GRIDTYPE,DTYPE,FLTYPE,NPTYPE,AVG2DTYPE,AVG1DTYPE,ENERTYPE,LOSSTYPE,SPTYPE,TSTYPE,LOGDTTYPE,STEPTYPE,PERFTYPE,ADTYPE,PDTYPE,CALCTYPE,FLINETYPE,MODETYPE,EXPANDTYPE,NPCOMPUTETYPE;

int restartener,restartloss,restartmode,restartsonicpoint;

SFTYPE PERFWALLTIME,ZCPSESTIMATE;

// MPI stuff
int numprocs, myid,procnamelen;
int ncpux1,ncpux2,ncpux3;
#if(USEMPI)
char processor_name[MPI_MAX_PROCESSOR_NAME];
MPI_Status mpichstatus;
MPI_Group MPI_GROUP_WORLD;
MPI_Group grprem[6];
MPI_Comm combound[6]; // communicator for physical boundary cpus
#endif


FTYPE mms[ITYPES][CTYPES][NUMSCA+1][2];
FTYPE mmv[ITYPES][CTYPES][NUMVEC+1][3+1][2]; // 0 is magnitude
FTYPE mmst[ITYPES][CTYPES][NUMSCA+1][2];
FTYPE mmvt[ITYPES][CTYPES][NUMVEC+1][3+1][2];


// much simpler for LOOPTYPE>1:
// only over active grid
SFTYPE lossflux[NUMLOSSVAR+1][2];

FTYPE funcool[NUMFUNCOOL+1]; // 0..N-1 cooling function holder, like tot_loss holder
FTYPE thetai[NUMFUNCOOL+1]; // 0..N-1 cooling parameter holder, like tot_loss holder
FTYPE spectrum[NUMSPECTRUM]; // 0..N-1 spectrum holder, like tot_loss holder

FTYPE funrelie[NUMFUNRELIE+1]; // 0..N-1 
FTYPE thetareli[NUMFUNRELIE+1]; // 0..N-1 

int globalinterpmod;
int gocont; // used to continue running(runtype, directinput, timereenter)
int npoldschool; // identify npdump file as older version
int avg2doldschool; 
int avg1doldschool; 
SFTYPE Rint; // torus inner edge
SFTYPE Rinner,Router; // blackhole inner edge

int numpdumps,numdumps,numfldumps,numadumps,numfloordumps,numimages,numnpdumps,numcalcdumps;


FTYPE dxTVD,dyTVD,dVTVD ;
int firsttimetvdlf;

int temptempi;
char DATEXT[20];
char DAT2EXT[20];
char PAREXT[20];
char PAR2EXT[20];
char INEXT[20];
char OUTEXT[20];
char PPEXT[20];
char IMGEXT[20];


// for collaboration TORITYPE==3
SFTYPE DELTAR,HOR,RMAX;

FTYPE pos_accretor[4];
FTYPE radialdist_accretor;
FTYPE DX_single,DX_accretor;
int NUMZONE_accretor;
FTYPE VHORIZON[4],VACCRETOR[4],MACHHORIZON;

int passivepos;

FTYPE Omegasystem,rcm;

#if(KAICOOL)
int (*thingrid)[N2M][N1M],thingrida[N3M][N2M][N1M];
FTYPE dtsmall;
FTYPE (*coolerKai)[N2M][N1M],coolerKaia[N3M][N2M][N1M];
FTYPE (*sthin)[N2M][N1M],sthina[N3M][N2M][N1M];
FTYPE (*tot)[N3M][N2M][N1M],tota[2][N3M][N2M][N1M];
FTYPE optda[N3M][N2M][N1M],tatm4a[N3M][N2M][N1M],tkaia[N3M][N2M][N1M],kapka[N3M][N2M][N1M],sigka[N3M][N2M][N1M],fbdrya[N3M][N2M][N1M],taoupa[N3M][N2M][N1M],tenva[N1M];
FTYPE (*optd)[N2M][N1M],(*tatm4)[N2M][N1M],(*tkai)[N2M][N1M],(*kapk)[N2M][N1M],(*sigk)[N2M][N1M],(*fbdry)[N2M][N1M],(*taoup)[N2M][N1M],(*tenv);

int jsmall,ismall;
#endif
#if(DIAGF)
FTYPE (*diagn)[N3M][N2M][N1M],diagna[33][N3M][N2M][N1M];
#else
FTYPE (*diagn)[N3M][N2M][N1M],diagna[6][N3M][N2M][N1M];
#endif
//reenter zz, need to be changed when using different dump to reenter
FTYPE xrea[NUMGRID][3][NREBIG];
FTYPE (*xre)[3][NREBIG];
FTYPE datrea[39][N3M][N2RE][N1RE];
FTYPE (*datre)[N3M][N2RE][N1RE];
int ire,jre,kre;

//surface density, sigma
FTYPE Sigara[N3M][N2M][N1M],cs2aa[N3M][N2M][N1M];
FTYPE (*Sigar)[N2M][N1M],(*cs2a)[N2M][N1M];
int irowmpi;
//change step
FTYPE totstepchange,stepchange,totstep,diagstep;
// substep
int substepm;
FTYPE optthick;
#include "radiation2.h"
