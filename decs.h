//unit conversion
extern SFTYPE munit,lunit,rhounit,tunit,eneunit,enedenunit;
//Constant
extern SFTYPE lsuncgs,kcontcgs,mhydrcgs,sigmacontcgs,sunmasscgs,yrcgs,mstarcgs;
extern SFTYPE sigmacontnew,kcontnew,mhydrnew;
extern SFTYPE rsun,mmw;
//Central star
extern SFTYPE lsunfrac,lstarcgs,rstar;
// layer model
extern SFTYPE Siga,cstcrit2,alphafrac;
extern FTYPE iab[N3M][N2M][N1M];
extern int actcontour;
//Temperature,opaini
extern FTYPE tirr0a[N3M][N2M][N1M],opainia[N3M][N2M][N1M],siginia[N3M][N2M][N1M],optdinia[N3M][N2M][N1M];
extern FTYPE cs02a[N3M][N2M][N1M];
extern FTYPE (*tirr0)[N2M][N1M],(*opaini)[N2M][N1M],(*sigini)[N2M][N1M],(*optdini)[N2M][N1M];
extern FTYPE (*cs02)[N2M][N1M];
extern FTYPE teff4[N1];
// make sure all this matches init.c: init_pointers()


/* 1: a: + dx in 1/2/3-direction from cell face: x(i)=x(i-1)+dx(i-1) */
/* 2: b: + dx in 1/2/3-direction from cell center: x(i)=x(i-1)+dx(i) */
extern FTYPE dxa[NUMGRID][3][NBIGM]; // [agrid=1,bgrid=2][dir=1,2,3][i=-NBIGBND..NBIG+NBIGBND]
extern FTYPE (*dx)[3][NBIGM] ;

/* Positions and sizes of cell */
/* cell face position 1: a */
/* cell center position 2: b*/
extern FTYPE xa[NUMGRID][3][NBIGM]; // [agrid=1,bgrid=2][dir=1,2,3][i=-NBIGBND..NBIG+NBIGBND]
extern FTYPE (*x)[3][NBIGM] ;

// layered 
FILE *grid_file;
extern FTYPE alpha_realar[N3M][N2M][N1M];
FTYPE(*alpha_reala)[N2M][N1M];

// only interp if pp
#if(POSTPROC==1)
extern FTYPE idxa[NUMGRID][3][INTNBIG];
extern FTYPE ixa[NUMGRID][3][INTNBIG];
extern FTYPE ppcalcsa[NUMPPCALCS][N3M][N2M][N1M] ;
#endif
// keep below so code is intact
extern FTYPE (*idx)[3][INTNBIG] ;
extern FTYPE (*ix)[3][INTNBIG] ;
extern FTYPE (*ppcalcs)[N3M][N2M][N1M];

/* at cell center 1: d 2: e 3: phi */
extern FTYPE sa[NUMSCA][N3M][N2M][N1M] ; // [rho=1,en=2,pot=3][k=-N3BND..N3+B3BND][j=-N2BND..N2+B2BND][i=-N1BND..N1+B1BND]
extern FTYPE (*s)[N3M][N2M][N1M] ; /* 3 scalars */

/* at cell face, where position of cell-face var is less than cell center var at same index */
/* 1: v 2: B*/
extern FTYPE va[NUMVEC][3][N3M][N2M][N1M] ; //[velocity=1,B-field=2][1=x1dir,2=x2dir,3=x3dir][k=-N3BND..N3+B3BND][j=-N2BND..N2+B2BND][i=-N1BND..N1+B1BND]
extern FTYPE (*v)[3][N3M][N2M][N1M] ;

////// analytic storage
/* at cell center 1: d 2: e 3: phi */
extern FTYPE sanala[NUMSCA][N3M][N2M][N1M] ; // same as s
extern FTYPE (*sanal)[N3M][N2M][N1M] ;

#if(GRAVACC&&(!POSTPROC))
extern FTYPE gravacca[3][N3M][N2M][N1M] ; // analytic gravitational accel for each direction: assumes only 1 direction currently
#endif
extern FTYPE (*gravacc)[N3M][N2M][N1M] ;

#if(GRAVITOMAGNETIC&&(!POSTPROC))
extern FTYPE gravitoma[2][N3M][N2M][N1M] ; // geometry terms in gravitomagnetic source term
#endif
extern FTYPE (*gravitom)[N3M][N2M][N1M] ;

/* at cell face, where position of cell-face var is less than cell center var at same index */
/* 1: v 2: B*/
extern FTYPE vanala[NUMVEC][3][N3M][N2M][N1M] ; // same as v
extern FTYPE (*vanal)[3][N3M][N2M][N1M] ;

#if(DOINVSOL2FUN)
extern FTYPE invsol2funa[N3M][N2M][N1M] ; // 1/c^2 as a function of something
#endif
extern FTYPE (*invsol2fun)[N2M][N1M] ;

#if(DOAVGDIAGMEMORY)
// 2-D, averaged over 3
/* at original cell location */
// [rho=1,en=2,Be=3,cs^2=4,S=5,|vx1|=6,|vx2|=7,|vx3|=8]
//[j=-N2BND..N2+B2BND][i=-N1BND..N1+B1BND]
extern FTYPE avg2_3a[NUMAVG2_3][N2M][N1M] ;

// 1-D, averaged over 3 and 1 direction
/* at original cell location */
// [rho=1,en=2,Be=3,cs^2=4,S=5,|vx1|=6,|vx2|=7,|vx3|=8]
//[j=-N2BND..N2+B2BND]
extern FTYPE avg1_31a[NUMAVG1_31][N2M] ;

// 1-D, averaged over 3 and 2 direction
/* at original cell location */
// [rho=1,en=2,Be=3,cs^2=4,S=5,|vx1|=6,|vx2|=7,|vx3|=8]
//[i=-N1BND..N1+B1BND]
extern FTYPE avg1_32a[NUMAVG1_32][N1M] ;
#endif
extern FTYPE (*avg2_3)[N2M][N1M] ;
extern FTYPE (*avg1_31)[N2M] ;
extern FTYPE (*avg1_32)[N1M] ;

#if(FLOORDUMPFLAGMEMORY)
// used for floor on grid for scalar values.
// 1: rho->rho
// 2: u->u
// 3: phi->rho*phi (gravitational potential energy density)
extern FTYPE floorvarsa[NUMSCA][N3M][N2M][N1M] ;

// used for floor on grid for vectors(note: need 0 element here)
// all due to rho changes
// 1: momentum density in radial direction
// 2: theta-angular momentum density
// 3: phi-angular momentum density
extern FTYPE floorvarva[NUMVEC][3][N3M][N2M][N1M] ;

// 0: KE
// 1: BE
extern FTYPE floorvar0a[NUMVEC][N3M][N2M][N1M] ;

#endif
extern FTYPE (*floorvars)[N3M][N2M][N1M] ;
extern FTYPE (*floorvarv)[3][N3M][N2M][N1M] ;
extern FTYPE (*floorvar0)[N3M][N2M][N1M] ;

#if((LOWMEMMODE==0)&&(TVDLF==0))
// comp savers (first[] is location of variable used on)
// direction specifies direction of surface normal vector
extern FTYPE dsa[4][3][N3M][N2M][N1M]; // [1=(.5,5),2=(0,.5),3=(.5,0),4=(0,0)][dir=1,2,3][k=-N3BND..N3+N3BND-1][j=-N2BND..N2+N2BND-1][i=-N1BND..N1+N1BND-1]

extern FTYPE oarcla[4][3][N3M][N2M][N1M]; // [1=(.5,5),2=(0,.5),3=(.5,0),4=(0,0),5=(.5,.5,0)][dir=1,2,3][k=-N3BND..N3+N3BND-1][j=-N2BND..N2+N2BND-1][i=-N1BND..N1+N1BND-1]

// just dx term, not arclength
extern FTYPE odxa[NUMGRID][3][NBIGM]; // [agrid=1,bgrid=2][dir=1,2,3][i=-NBIGBND..NBIG+NBIGBND]

extern FTYPE ovola[4][N3M][N2M][N1M]; // [1=(.5,5),2=(0,.5),3=(.5,0),4=(0,0)][k=-N3BND..N3+N3BND-1][j=-N2BND..N2+N2BND-1][i=-N1BND..N1+N1BND-1]
#endif


  /* metric scales and their dervs--only applies to limited coord systems */
  /* 1: a-mesh */
  /* 2: b-mesh */
  /* 1: g1 2: g2 3: g31 4: g32 */
#if(LOWMEMMODE==0)
extern FTYPE ga[NUMGRID][NUMMETRIC][NBIGM]; // [1,2][1,2,3,4][...]
extern FTYPE dga[NUMGRID][NUMMETRIC][NBIGM];  // [1,2][1,2,3,4][...]
extern FTYPE oga[NUMGRID][NUMMETRIC][NBIGM]; // inverse of g

  /* Volume factors: dvl[a/b][dir][location] */
extern FTYPE dvla[NUMGRID][3][NBIGM]; //[1,2][1,2,3][...]
extern FTYPE odvla[NUMGRID][3][NBIGM]; //[1,2][1,2,3][...]
#endif

#if((PUREBC==0)&&(BOUNDTYPE==1))
extern short bcsa[NUMSCA][3][N3M][N2M][N1M]; // [which scalar][1=type:assign=(1,2,3,4,5),2=dim:assign=(1,2,3),3=direction:assign=(-1,+1)][...][...][...]
//bcv assumes all components bounded with same *type* of boundary condition
extern short bcva[NUMVEC][3][N3M][N2M][N1M]; // [which vector][1=type:assign=(1,2,3,4,5),2=dim:assign=(1,2,3),3=direction:assign=(-1,+1)][...][...][...]

#elif((BOUNDTYPE==2)||(BOUNDTYPE==3))

extern short bzmaska[N3M][N2M][N1M];
extern short bzmaskina[N3M][N2M][N1M];
// 0: comp zone
// 1: normal bzone
// -1: no-where zone(never used)
// 99: MPI zone
// 3 directions, 8 maximal values to average
// specifies location of zone to copy from

#endif

extern short (*bcv)[3][N3M][N2M][N1M];
extern short (*bcs)[3][N3M][N2M][N1M];

extern short (*bzmask)[N2M][N1M];
extern short (*bzmaskin)[N2M][N1M];
extern short (***bzs); // worst case scenario is 26 positions per zone to average
// pointer to pointer for an array of pointers

extern FTYPE (*ds)[3][N3M][N2M][N1M] ;
extern FTYPE (*oarcl)[3][N3M][N2M][N1M] ;
extern FTYPE (*odx)[3][NBIGM] ;
extern FTYPE (*ovol)[N3M][N2M][N1M] ;
extern FTYPE (*g)[NUMMETRIC][NBIGM];
extern FTYPE (*og)[NUMMETRIC][NBIGM];
extern FTYPE (*dg)[NUMMETRIC][NBIGM];
extern FTYPE (*dvl)[3][NBIGM];
extern FTYPE (*odvl)[3][NBIGM];


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
extern short iindxa[NUMINDEX][2][N3M][N2M];
#endif
extern short (*iindx)[2][N3M][N2M];



#if(VISCMEM)
extern FTYPE nu_facta[N3M][N2M][N1M] ;
extern FTYPE nu_reala[N3M][N2M][N1M] ;
#endif
extern FTYPE (*nu_fact)[N2M][N1M] ;
extern FTYPE (*nu_real)[N2M][N1M] ;


#if(RESMEM&&(POSTPROC==0))
extern FTYPE nu_res_facta[N3M][N2M][N1M] ;
extern FTYPE nu_res_reala[N3M][N2M][N1M] ;
extern FTYPE jcurrenta[3+1][N3M][N2M][N1M] ;
#endif
extern FTYPE (*nu_res_fact)[N2M][N1M] ;
extern FTYPE (*nu_res_real)[N2M][N1M] ;
extern FTYPE (*jcurrent)[N3M][N2M][N1M] ;


#if((LOOPTYPE==1)&&(DOLOSSDIAGMEMORY&&(POSTPROC==0)))
// only over active grid
extern SFTYPE losssa[NUMSCA][3][2][NBIG][NBIG]; // [which=1...][dir=1,2,3][i/o=0,1][list=0...]
extern SFTYPE lossva[NUMVEC][3+1][3][2][NBIG][NBIG]; // [which][v-dir][dir][i/o][list]
extern SFTYPE lossvisca[1][3][2][NBIG][NBIG]; // [which][dir][i/o][list]
#endif
extern SFTYPE (*losss)[3][2][NBIG][NBIG];
extern SFTYPE (*lossv)[3+1][3][2][NBIG][NBIG];
extern SFTYPE (*lossvisc)[3][2][NBIG][NBIG];
// v-dir: vector direction // 0=energy, 1,2,3
// dir: 1=x1 2=x2 3=x3
// i/o: 0=inner boundary 1=outer boundary
// list: var aextern long that side: on active grid: 0..N-1(i.e. don't shift pointer!)


extern short accountstorea[N3M][N2M][N1M];
extern short (*accountstore)[N2M][N1M];

#if(MDOTMEM&&(POSTPROC==0))
extern FTYPE rhoinjecta[N3M][N2M][N1M];
extern FTYPE eninjecta[N3M][N2M][N1M];
#endif
extern FTYPE (*rhoinject)[N2M][N1M];
extern FTYPE (*eninject)[N2M][N1M];

// used to hold boundary data(packed)
#if((TVDLF==0)&&(POSTPROC==0))
extern FTYPE worksbca[NUMSCA][2][COMPDIM*2][NBIGBND*NBIGSM]; //scalar: [scalar#][1=out/2=in][0=right,1=up,2=left,3=down,4=out,5=in][datawidth]
extern FTYPE workvbca[NUMVEC][2][COMPDIM*2][3*NBIGBND*NBIGSM]; //vector: [vector#][as above][as above][datawidth]
#endif
extern FTYPE (*worksbc)[2][COMPDIM*2][NBIGBND*NBIGSM];
extern FTYPE (*workvbc)[2][COMPDIM*2][3*NBIGBND*NBIGSM];




// workXa is just like s[i] (i.e. one particular scalar)
extern FTYPE work1a[N3M][N2M][N1M];
extern FTYPE work2a[N3M][N2M][N1M];
extern FTYPE work3a[N3M][N2M][N1M];
#if(TVDLF==0)
extern FTYPE work4a[N3M][N2M][N1M];
extern FTYPE work5a[N3M][N2M][N1M];
extern FTYPE work6a[N3M][N2M][N1M];
extern FTYPE work7a[N3M][N2M][N1M];
extern FTYPE work8a[N3M][N2M][N1M];
extern FTYPE work9a[N3M][N2M][N1M];
extern FTYPE work10a[N3M][N2M][N1M];
#endif
extern FTYPE (*work1)[N2M][N1M];
extern FTYPE (*work2)[N2M][N1M];
extern FTYPE (*work3)[N2M][N1M];
extern FTYPE (*work4)[N2M][N1M];
extern FTYPE (*work5)[N2M][N1M];
extern FTYPE (*work6)[N2M][N1M];
extern FTYPE (*work7)[N2M][N1M];
extern FTYPE (*work8)[N2M][N1M];
extern FTYPE (*work9)[N2M][N1M];
extern FTYPE (*work10)[N2M][N1M];

#if(MDOTMEMANAL&&(POSTPROC==0))
extern FTYPE mdotanala[N3M][N2M][N1M];
#endif
extern FTYPE (*mdotanal)[N2M][N1M];

extern FTYPE osqrtrhoa[N3M][N2M][N1M];
extern FTYPE (*osqrtrho)[N2M][N1M];

#if(ALFVENLIMIT)
extern FTYPE rholimiteda[3+1][N3M][N2M][N1M];
#endif
extern FTYPE (*rholimited)[N3M][N2M][N1M];

#if(TVDLF==0)
// workvXa is just like v[i] (i.e. like one particular scalar)
extern FTYPE workv1a[3][N3M][N2M][N1M];
extern FTYPE workv2a[3][N3M][N2M][N1M];
extern FTYPE workv3a[3][N3M][N2M][N1M];
extern FTYPE workv4a[3][N3M][N2M][N1M];
extern FTYPE workv5a[3][N3M][N2M][N1M];
#endif
extern FTYPE (*workv1)[N3M][N2M][N1M];
extern FTYPE (*workv2)[N3M][N2M][N1M];
extern FTYPE (*workv3)[N3M][N2M][N1M];
extern FTYPE (*workv4)[N3M][N2M][N1M];
extern FTYPE (*workv5)[N3M][N2M][N1M];

// [1,2,3][1,2,3][...][...][...]
// below 3 tensors only used if real viscosity is present
#if(VISCMEM)
extern FTYPE sigmaa[3][3][N3M][N2M][N1M];
extern FTYPE rosta[3][3][N3M][N2M][N1M];
extern FTYPE rostnua[3][3][N3M][N2M][N1M];
#endif
extern FTYPE (*sigma)[3][N3M][N2M][N1M]; // -2*rho*nu*e_{ij}= sigma_{ij}
extern FTYPE (*rost)[3][N3M][N2M][N1M];// e_{ij}
extern FTYPE (*rostnu)[3][N3M][N2M][N1M]; // e_{ij}*nu

// workiqa is like s but in 2D and only on active grid [0..n2-1][0..n1-1]
// only postproc interpolates or images averages or kinetic energy floor
#if(POSTPROC==1)
extern FTYPE workiqa[NUMSCA][INTN2][INTN1];
// workviqa is like v but in 2D and only on active grid
extern FTYPE workviqa[NUMVEC][3][INTN2][INTN1];
extern FTYPE workiqavga[NUMWORKIQ][INTN2][INTN1];
extern FTYPE work0iqa[NUMVEC][INTN2][INTN1];
#endif
// still allocate pointers so code is intact except init.c pointer stuff
extern FTYPE (*workiq)[INTN2][INTN1];
extern FTYPE (*workviq)[3][INTN2][INTN1];
extern FTYPE (*workiqavg)[INTN2][INTN1];
extern FTYPE (*work0iq)[INTN2][INTN1];


// to speed up things use matrix for const grids
// assumes never interpolate when postproc=0
#if(POSTPROC==1)
extern FTYPE convmata[3][3][INTN2][INTN1];  
#endif
extern FTYPE (*convmat)[3][INTN2][INTN1];



#if((TVDLF==1)&&(POSTPROC==0))
/* primitive variables */
extern FTYPE a_p[N1+2][N2+2][NP] ;
extern FTYPE a_dqx[N1+2][N2+2][NP] ;
extern FTYPE a_dqy[N1+2][N2+2][NP] ;
extern FTYPE a_Fx[N1+2][N2+2][NP] ;
extern FTYPE a_Fy[N1+2][N2+2][NP] ;
extern FTYPE a_ph[N1+2][N2+2][NP] ;
#endif
extern FTYPE (*  p)[N2+2][NP] ;
extern FTYPE (*dqx)[N2+2][NP] ;
extern FTYPE (*dqy)[N2+2][NP] ;
extern FTYPE (* Fx)[N2+2][NP] ;
extern FTYPE (* Fy)[N2+2][NP] ;
extern FTYPE (* ph)[N2+2][NP] ;


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// EVERYTHING BELOW HERE NEEDS NO POINTER SHIFTING OR MEMORY ALLOCATION OR DEPENDS ON SIZE OF SYSTEM
//
//

extern short *indx[NUMINDEX]; // array of pointers
extern int numiter[NUMINDEX]; // 0-NUMINDEX-1


// counts at different 5 different places in code(see step.c and sweep.c)
// these 2 vars are only sources for floor changes in above arrays
extern int floorcnt[NUMFLOOROUT+1][NUMFLOORVAR+1]; //[0,1,2,3,4,5,6][1=rho,2=en]
// [0][] : step_visc
// [1][] : step_real_visc
// [2][] : step_ie
// [3][] : step_res
// [4][] : sweepx
// [5][] : sweepy
// [6][] : sweepz


// stores lowest value of rho and en before floor fixed, over entire run
extern FTYPE floorlowest[NUMFLOORVAR+1]; //[1=rho,2=en]

// what routine has the above lowest value generated
// same routine numbers as above
extern int wherelowest[NUMFLOORVAR+1]; //[1=rho,2=en]


//extern int N[3+1] ; /* Number of cells in each direction */
extern SFTYPE L[2+1][3+1] ; /* L[1] holds position of 1,1,1 cell outer edges, L[2] holds size of grid in units */
extern SFTYPE iL[2+1][3+1];

// outerdefs for images
extern SFTYPE outerdefs[ITYPES][CTYPES][NUMSCA+1];
extern SFTYPE outerdefv[ITYPES][CTYPES][NUMVEC+1][3+1]; // [][0] here is magnitude

// outerdef for dumps for interpolation
extern SFTYPE douterdefs[NUMSCA+1];
extern SFTYPE douterdefv[NUMVEC+1][3+1]; // [][0] here is magnitude
extern SFTYPE douterdef0[NUMVEC+1]; // KE for floor(not set yet)


extern SFTYPE IOBound[NUMBTRANS]; // used to have inflow/outflow transitions on outer theta edge

// average stuff needed elsewhere too
extern SFTYPE tavgstart,tavgfinal; // real start and final
extern int avgcount; //number of averages in time
extern int num1d_31_full,num1d_32_full;

extern SFTYPE floors[NUMLOSSVAR+1]; // to keep track of floor input of stuff
extern SFTYPE inflows[NUMLOSSVAR+1]; // to keep track of injection of stuff
extern SFTYPE radiations[NUMLOSSVAR+1]; // to keep track of radiation of stuff
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

extern SFTYPE invsol2,invsol2_original; // inverse speed of light, squared, for alfven limit trick

// next block of changable global parameters(global.h just sets initial value)
extern int kever;
extern int mdotin,cool,press,pressx1,pressx2,pressx3,visc_art,ie,mag,pmag,res_real,visc_real,trans;
extern int stepmocct,mocctvx1,mocctvx2,mocctvx3,mocctbx1,mocctbx2,mocctbx3,transbzy,transbzx;
extern int vischeat;
extern int resheat;
extern int vreal;
extern int rreal;
extern int transx1,transx2,transx3;
extern int transiex1,transmagx1,transv1x1,transv2x1,transv3x1,transrhox1;
extern int transiex2,transmagx2,transv1x2,transv2x2,transv3x2,transrhox2;
extern int transiex3,transv1x3,transv2x3,transv3x3,transrhox3;
extern int dopassive;
extern int transpassivex1,transpassivex2,transpassivex3;
extern int transpassivev1x1,transpassivev1x2,transpassivev1x3,transpassivev2x1,transpassivev2x2,transpassivev2x3,transpassivev3x1,transpassivev3x2,transpassivev3x3;
extern int advint;

extern int ireenter;

extern SFTYPE coolbottom; // bottom on cooling in theta
extern SFTYPE rho0,rho00,u0,R0,Omega0; // for tori problem and for SPB rho visc
extern SFTYPE coolfact;
extern char WRITETYPE[10];
extern int reallaststep;
extern int tagik,tagfk,tagij,tagfj,tagii,tagfi;
extern int t2gik,t2gfk,t2gij,t2gfj,t2gii,t2gfi;
extern int t3gik,t3gfk,t3gij,t3gfj,t3gii,t3gfi;
extern int t4gik,t4gfk,t4gij,t4gfj,t4gii,t4gfi;
extern char myidtxt[MAXFILENAME];
extern int totalzones,realtotalzones;
extern int rtotalzones;
extern int itotalzones;
extern int sizes[3+1][MAXCPUS];
extern int isizes[3+1][MAXCPUS];
extern int totalsize[3+1];
extern int itotalsize[3+1];
extern int mycpupos[3+1]; // my position amongst the cpus
extern int srdir[6]; // which direction this cpu sends/receives normal interior data 
extern int startpos[3+1];
extern int endpos[3+1]; // startj and endj are where this CPU located on full grid
extern long nstep;
extern int runtype,directinput,appendold,deleteolddat,deleteoldpar,deleteoldimg;
extern int PDUMP_START,DUMP_START,NPDUMP_START,ADUMP_START,FLOOR_START,FLDUMP_START,IMAGE_START;
extern int simplebc,bcix1,bcox1,bcix2,bcox2,bcix3,bcox3;
extern int nonunigridx1,nonunigridx2,nonunigridx3;
extern int analoutput;
extern SFTYPE x1in,x1out,x2in,x2out,x3in,x3out;
extern SFTYPE ix1in,ix1out,ix2in,ix2out,ix3in,ix3out;
extern SFTYPE timereenter;
extern int pdump_start,dump_start,fldump_start,npdump_start,adump_start,floor_start,image_start;
extern FILE*  fail_file;
extern FILE*  logfull_file;
extern FILE*  log_file;
extern FILE*  logstep_file;
extern FILE*  logperf_file;
extern FILE*  logdt_file;
extern FILE*  logfl_file;
extern FILE*  logsp_file;
// Equation of state coefficients
extern SFTYPE tscycleto; // time to subcycle to
extern SFTYPE tscyclefrom; // time to subcycle from
extern SFTYPE dtlastscycle; // dt used as basis for subcycle
extern int subcyclen; // number of subycycles for viscosity
extern int nthsubcycle; // number of subcycles so far
extern SFTYPE rg,rgp;
extern SFTYPE cour,invcour ; // invcour used in timestep1.h
extern SFTYPE cour2,invcour2 ; // invcour2 used in timestep1.h
extern SFTYPE cs ;
extern SFTYPE alpha;
extern SFTYPE gam ;
extern int wgam;
extern int wgam1;
extern int wgam53;
extern int wpw;
extern int periodicx1,periodicx2,periodicx3;
extern int mpiperiodicx1,mpiperiodicx2,mpiperiodicx3;
extern int periodicx2special;
extern int skipix1,reflectix1,reflectox1;
extern int skipix2,reflectix2,reflectox2;
extern int skipix3,reflectix3,reflectox3;
extern int intix1,intox1,intix2,intox2,intix3,intox3;
extern int skipintix1,skipintix2,skipintix3;
extern SFTYPE dt ;
extern SFTYPE tstart,TSTART;
extern SFTYPE t,tf,timagescale ;
extern SFTYPE tavgi,tavgf;
extern int numavg;
extern SFTYPE nu_vnr ;
extern SFTYPE nu_l ;
extern SFTYPE nu_ten ;
extern SFTYPE alpha_real ;
extern SFTYPE alpha_real0 ;
extern SFTYPE n_real ;
extern SFTYPE GRAVC;
extern SFTYPE MASSBH ;
extern SFTYPE GM;
extern SFTYPE DTd ;
extern SFTYPE DTfld ;
extern SFTYPE DTpd ;
extern SFTYPE DTavg ;
extern SFTYPE DTi ;
extern SFTYPE DTl ;
extern SFTYPE DTener;
extern SFTYPE DTloss;
extern SFTYPE DTmode;
extern SFTYPE DTfloor;
extern SFTYPE DTtimestep;
extern SFTYPE DTdivb;
extern SFTYPE DTtimescale;
extern SFTYPE DTsp;
extern SFTYPE DTstep,DTstepdot,DTperf,DTgocheck,DTtimecheck;
extern SFTYPE resist_real0 ;
extern SFTYPE resist_real ;
extern SFTYPE nu_sh ;
//extern int numbc[3+1];
extern SFTYPE startx[3+1]; // starting location for each cpu for each direction
extern SFTYPE startix[3+1]; // starting location for each cpu for each direction
/* grid velocity */
extern SFTYPE vg[3+1];
extern SFTYPE RHO0,V0; // for advection
extern SFTYPE DENSITYFLOOR,IEFLOOR,DENSITYFRACFLOOR,IEFRACFLOOR;
extern SFTYPE IEFRACT,VZFRACT;
extern SFTYPE massdot;

extern SFTYPE divbmax_full,divbavg_full;

extern int DYNAMICMM;

extern SFTYPE blackholejz;

extern int mpicombine;
extern int binaryoutput,cpugeompick,CHECKCONT,FLUSHFAILDT,FORCERHO,FORCEIE,FORCERHOINTERNAL,FORCEIEINTERNAL,TRYSUBCYCLE,CHECKDTLOW;
extern SFTYPE DTOTLIMIT,INVDTOTLIMIT,DTLOWEST,IDTLOWEST,SQIDTLOWEST,DTLOWDUMP;
extern int DOGPARDIAG,DOPARDIAG,DOIPARDIAG,DOGENDIAG,DOLOSSDIAG,COMPUTELOSSDIAG,HYDROBZFLUX,DETAILMLOSS,DOAVGDIAG,DOFLOORDIAG,DOFLOORD2,DODTDIAG,DOTSTEPDIAG,DOSPDIAG,DOLOGSTEP,DOLOGPERF;

extern int DODIVBDIAG;

extern int NDTCCHECK,NZCCHECK,NDTDOTCCHECK,NGOCHECK,NTIMECHECK;


extern int PDUMPFLAG,DUMPFLAG,NPDUMPFLAG,FLOORDUMPFLAG,FLDUMPFLAG,ADUMPFLAG,IMAGEFLAG;
extern int OLDIMAGEFORMAT,OLDSCHOOL,OLDSCHOOL2;

extern char DATADIR[MAXFILENAME],DUMPDIR[MAXFILENAME],IMAGEDIR[MAXFILENAME],IMAGEIDIR[MAXFILENAME];

extern int DUMPSM,FULLOUTPUT,FULLINPUT;

extern int SAMEGRIDI,SAMEGRIDD;
extern int IMAGEFORMAT,IMAGEFORMATINPUT,GZIPIMAGE,GZIPIMAGEINPUT;

extern int NOBOUNDPOT;

extern int INFLOWCHECKIX1,INFLOWCHECKOX1,INFLOWCHECKIX2,INFLOWCHECKOX2,INFLOWCHECKIX3,INFLOWCHECKOX3;

extern int PVER,GRIDVER,DVER,FLVER,NPVER,AVG1DVER,AVG2DVER,ENERVER,MODEVER,LOSSVER,SPVER,TSVER,LOGDTVER,STEPVER,PERFVER,ADVER,PDVER,CALCVER,FLINEVER;

extern int PTYPE,GRIDTYPE,DTYPE,FLTYPE,NPTYPE,AVG2DTYPE,AVG1DTYPE,ENERTYPE,LOSSTYPE,SPTYPE,TSTYPE,LOGDTTYPE,STEPTYPE,PERFTYPE,ADTYPE,PDTYPE,CALCTYPE,FLINETYPE,MODETYPE,EXPANDTYPE,NPCOMPUTETYPE;

extern int restartener,restartloss,restartmode,restartsonicpoint;

extern SFTYPE PERFWALLTIME,ZCPSESTIMATE;

// MPI stuff
extern int numprocs, myid,procnamelen;
extern int ncpux1,ncpux2,ncpux3;
#if(USEMPI)
extern char processor_name[MPI_MAX_PROCESSOR_NAME];
extern MPI_Status mpichstatus;
extern MPI_Group MPI_GROUP_WORLD;
extern MPI_Group grprem[6];
extern MPI_Comm combound[6]; // communicator for physical boundary cpus
#endif


extern FTYPE mms[ITYPES][CTYPES][NUMSCA+1][2];
extern FTYPE mmv[ITYPES][CTYPES][NUMVEC+1][3+1][2]; // 0 is magnitude
extern FTYPE mmst[ITYPES][CTYPES][NUMSCA+1][2];
extern FTYPE mmvt[ITYPES][CTYPES][NUMVEC+1][3+1][2];


// much simpler for LOOPTYPE>1:
// only over active grid
extern SFTYPE lossflux[NUMLOSSVAR+1][2];

extern FTYPE funcool[NUMFUNCOOL+1]; // 0..N-1 cooling function holder, like tot_loss holder
extern FTYPE thetai[NUMFUNCOOL+1]; // 0..N-1 cooling parameter holder, like tot_loss holder
extern FTYPE spectrum[NUMSPECTRUM]; // 0..N-1 spectrum holder, like tot_loss holder

extern FTYPE funrelie[NUMFUNRELIE+1]; // 0..N-1 
extern FTYPE thetareli[NUMFUNRELIE+1]; // 0..N-1 

extern int globalinterpmod;
extern int gocont; // used to continue running(runtype, directinput, timereenter)
extern int npoldschool; // identify npdump file as older version
extern int avg2doldschool; 
extern int avg1doldschool; 
extern SFTYPE Rint; // torus inner edge
extern SFTYPE Rinner,Router; // blackhole inner edge

extern int numpdumps,numdumps,numfldumps,numadumps,numfloordumps,numimages,numnpdumps,numcalcdumps;


extern FTYPE dxTVD,dyTVD,dVTVD ;
extern int firsttimetvdlf;

extern int temptempi;
extern char DATEXT[20];
extern char DAT2EXT[20];
extern char PAREXT[20];
extern char PAR2EXT[20];
extern char INEXT[20];
extern char OUTEXT[20];
extern char PPEXT[20];
extern char IMGEXT[20];


// for collaboration TORITYPE==3
extern SFTYPE DELTAR,HOR,RMAX;

extern FTYPE pos_accretor[4];
extern FTYPE radialdist_accretor;
extern FTYPE DX_single,DX_accretor;
extern int NUMZONE_accretor;
extern FTYPE VHORIZON[4],VACCRETOR[4],MACHHORIZON;

extern int passivepos;

extern FTYPE Omegasystem,rcm;

#if(KAICOOL)
extern int (*thingrid)[N2M][N1M],thingrida[N3M][N2M][N1M];
extern FTYPE dtsmall;
extern FTYPE (*coolerKai)[N2M][N1M],coolerKaia[N3M][N2M][N1M];
extern FTYPE (*sthin)[N2M][N1M],sthina[N3M][N2M][N1M];
extern FTYPE (*tot)[N3M][N2M][N1M],tota[2][N3M][N2M][N1M];
extern FTYPE optda[N3M][N2M][N1M],tatm4a[N3M][N2M][N1M],tkaia[N3M][N2M][N1M],kapka[N3M][N2M][N1M],sigka[N3M][N2M][N1M],fbdrya[N3M][N2M][N1M],taoupa[N3M][N2M][N1M],tenva[N1M];
extern FTYPE (*optd)[N2M][N1M],(*tatm4)[N2M][N1M],(*tkai)[N2M][N1M],(*kapk)[N2M][N1M],(*sigk)[N2M][N1M],(*fbdry)[N2M][N1M],(*taoup)[N2M][N1M],(*tenv);

extern int jsmall,ismall;
#endif
#if(DIAGF)
extern FTYPE (*diagn)[N3M][N2M][N1M],diagna[33][N3M][N2M][N1M];
#else
extern FTYPE (*diagn)[N3M][N2M][N1M],diagna[6][N3M][N2M][N1M];
#endif
//reenter zz, need to be changed when using different dump to reenter
extern FTYPE xrea[NUMGRID][3][NREBIG];
extern FTYPE (*xre)[3][NREBIG];
extern FTYPE datrea[39][N3M][N2RE][N1RE];
extern FTYPE (*datre)[N3M][N2RE][N1RE];
extern int ire,jre,kre;

//surface density, sigma
extern FTYPE Sigara[N3M][N2M][N1M],cs2aa[N3M][N2M][N1M];
extern FTYPE (*Sigar)[N2M][N1M],(*cs2a)[N2M][N1M];
extern int irowmpi;
//change step
extern FTYPE totstepchange,stepchange,totstep,diagstep;
// substep
extern int substepm;
extern FTYPE optthick;
#include "radiation2.h"
