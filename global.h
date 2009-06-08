// compile using make -f makefile where makefile can be makefile or makefile.pp or makefile.ncsa, the first being the normal program, the second being the postprocess program, the 3rd the ncsa version of the normal program

/* Notes:
When memory is allocated for variables of size N I add the pointer address by NBND so the memory the pointer points to can be adDressed with the index list: -NBND..-1 , N+1..N+NBND for Boundary zones and 0..N for active zones.

BE CAREFUL with macros, always add parethesis around non-singular defines in case used later with multiplations, etc.

To start mpich run, first do:

ssh-agent /bin/bash
ssh-add 

then do:

mpirun -np 4 ./mhd3d

4 cpus using mpich:

rm nohup.out ; nohup sh -c 'mpirun -np 4 /usr/bin/time -v ./mhd3d > 0_o.out' &

2 cpus using mpich:

rm nohup.out ; nohup sh -c 'mpirun -np 2 /usr/bin/time -v ./mhd3d > 0_o.out' &


2 rainman cpus:
rm nohup.out ; nohup sh -c 'mpirun -machinefile /usr/local/share/rainman2.1 -np 2 ./mhd3d > 0_o.out' &

2 photon cpus:
rm nohup.out ; nohup sh -c 'mpirun -machinefile /usr/local/share/photon2.1 -np 2 ./mhd3d > 0_o.out' &

alphadog + wiseguy
rm nohup.out ; nohup sh -c 'mpirun -machinefile /usr/local/share/wisealpha.1 -np 2 ./mhd3d > 0_o.out' &

1 processor:

rm nohup.out ; nohup sh -c  '/usr/bin/time -v ./mhd3d > 0_o.out' &


Use program "nm" to list objects and symbols, or use objdump -axfhp 

A // MARK in the code means currently being tested as a modification
A // PRECISION in the code means there could be a precision problem(addition machine precision error)

*/

/*
  Generally, switches that are performance related are here, while others are in init with variable complements for flexibility during runtime.

  Not everything is in here so can compile quicker on simple changes.

*/

#include <errno.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
//#include "nrutil.h"

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
// whether in post processing mode or not
#define POSTPROC 0

#if(POSTPROC)
// GAMMIE SETTINGS
////////////////////
#define PPCLEAN 1 // 0: uses all elements of code 1: least needed elements of code
#define GAMMIEIMAGE 0 //whether setup for gammie config
#define INPUTN1 80 // TOTAL SIZE for MPI or not
#define INPUTN2 60
// output not relevant for combine, only interp
#define OUTPUTN1 80
#define OUTPUTN2 60
/////////////////////
#else
#define PPCLEAN 0
#define GAMMIEIMAGE 0
#endif

// whether doing a performance test(see: perftest.c, must come before mytime.h)
#define PERFTEST 0 // set PERFWALLTIME and ZCPSESTIMATE in gpar.c

#include "mytime.h"

#define GPARREAD 2
// 1 from file
// 2 from command line

// whether running in production mode
// turns off grid output in 3D since large and reproducible if necessary
#define PRODUCTION 0

#define LINUXCLUSTER 1 // fix for stupid pgCC compiler(uses extern-OMITTED method instead of COMMON-name method)

// whether to use MPI
// must change in makefile too!
#if(POSTPROC==0)
#define USEMPI 1 // choice
#else
#define USEMPI 0 // always 0
#endif


// whether to use MPI over GM
// This forces no forks to be called.  Can be used for non-gm too
#if(USEMPI==1)
#include "mpi.h"
#define USEGM 1 // choice
#define USEROMIO 1 // choice

#if(!USEROMIO)
#define USEJONIO 0 // choice
#else
#define USEJONIO 0 // no choice
#endif

#else
#define USEGM 0 // always 0, can't have GM without MPI
#define USEROMIO 0 // no choice
#define USEJONIO 0 // no choice
#endif

#define NCSA 0 // only for PERFTEST=1


#define COMPDIM 3 // how many computational, bounded dimensions
// 2: uses some specialized 2D-only code, like bzsweeps
// 3: general code.  Quite optimized even for 1D/2D cases, so should be generally used

#define SPLITMETHOD 0
// whether to use moc_ct split (1) or nonsplit (0)
// changes courant factor scaling with dimensions and diffusivity(split is much more diffusive, etc.)

#if(!GAMMIEIMAGE)

#define COORD 3 // normal choice

#else
#define COORD 1 // gammie
#endif
// 1 for cartesian coords
// 2 for cyl
// 3 for sph
// see init.c init_diffs() for usage

#define SUBCOOL 1
// Under MPI use ./math/mpichmcpusgood.nb to get per cpu size and distribution that's optimal.  Too messy of an equation to do PER CPU to minimize transfered elements on NCPUX1
#if(POSTPROC==0)
//MARK
// grid size PER CPU
#define N1 20  // must be even for totally general consistency
#define N2 32 // must be even for totally general consistency
#define N3 1 // should be 1 if COMPDIM==2, should be even if periodicx2special
// reenter grid size total
#define N1RE 80
#define N2RE 128
#define NREBIG ((N1RE>N2RE) ? N1RE : N2RE)
// consider IC/VIZ issues
#else
#if(!PPCLEAN)
#define N1 1 // should be total size, or make small for postproc.h hack(combine only)
#define N2 1 // should be total size
#define N3 1
#else

#define N1 INPUTN1 // should be total size, or make small for postproc.h hack(combine only)
#define N2 INPUTN2 // should be total size
#define N3 1

//#define N1 1
//#define N2 1
//#define N3 1

#endif

#endif











/* NBIG is bigger of N1 and N2 and N3 */
#define NBIG1 ((N1>N2) ? N1 : N2)
#define NBIG  ((NBIG1>N3) ? NBIG1 : N3)


#define CPUTXT ".%04d" // format for cpu data files

// whether to use TVDLF method
#define TVDLF 0
// only currently works for cartesian grid and periodic/outflow conditions

// size of data type used for all floats
// you must change all % f's to % lf's in global.h when using doubles
// you must change all % lf's to % f's in global.h when using floats!
#define FLOATTYPE 1 // 0: float 1: double (normal non-sensitive or performance critical datatypes)
#define SENSITIVE 1 // 0: float 1: double (non-perf critical or sensitive data types)

#ifndef FLT_EPSILON
#define FLT_EPSILON (1.E-7)
#endif
#ifndef DBL_EPSILON
#define DBL_EPSILON (2.E-16)
#endif

// need not change below datatype stuff
#if(FLOATTYPE==0)
#define NUMEPSILON FLT_EPSILON
#define FTYPE float
#define MPI_FTYPE MPI_FLOAT
#else
#define NUMEPSILON DBL_EPSILON
#define FTYPE double
#define MPI_FTYPE MPI_DOUBLE
#endif
#if(SENSITIVE==0) // for sensitive counters
#define SFTYPE float
#define MPI_SFTYPE MPI_FLOAT
#else
#define SFTYPE double
#define MPI_SFTYPE MPI_DOUBLE
#endif


#define DEBUG 0 // (temp stuff, none right now, moved to real seperate defines)
// 0: no debug statements: real run
// 1: normal run with some checks/statements
// 2: as 1 but add output on changes to data
#define DEBUGMPI 0
// 0: no debug statements: real run
// 1: check if runtime mpich calls are synched across all cpus


// with these memory flags, make sure you don't use mem if 0!!
#define MDOTMEM 0 // whether to allocate memory for injection
#define VISCMEM 1 // whether to allocate memory for viscosity tensors
// MARK
#define RESMEM 1 // whether to allocate memory for resistivity
#define MDOTMEMANAL 0 // whether to allocate memory for analytic mdot(see analsol.c)

// these pressure switches here for performance reasons
// check press==0/1 in init.c
#define PDEN                 1 // 0->no pressure 1->pressure from density/energy
#define PGRAV                1 // 0->no ext pot 1->ext pot
#define CURVE                1 // 0->no curvature terms 1-> do
#define PMAG                 1 // 0->no mag fields 1-> do mag fields
#define ROTATINGFORCE 0 // rotating force (see averywind())

// whether to have an accretor on the grid (like for averywind())
#define DO_ACCRETOR 0

#define CENT 0
#define VDIR1 1
#define VDIR2 2
#define VDIR3 3




// special opts or new computes for pressure step
// MARK
#if(COORD==3)
#define GRAVACC              0 // choice : 1->use analytic acceleration(see step.c and init.c) (must setup in analsol.c)
#else
#define GRAVACC              0 // choice -- user setup
#endif

#define GRAVITOMAGNETIC      0 // 1->gravitomagnetic field 0->no (must setup in analsol.c) // may be missing a term(sigma?)

// in bzsweepx/y:
#if(COMPDIM<3)
#define TORVX1    1 // bzsweepx v3-comp
#define TORVX2    1 // bzsweepy v3-comp
#define TORBX1    1 // bzsweepx b3-comp
#define TORBX2    1 // bzsweepy b3-comp
#else
#define TORVX1    0
#define TORVX2    0
#define TORBX1    0
#define TORBX2    0
#endif

#define RELIE      0 // 0-> Newtonian EOS 1-> Relativistic EOS

#define COOL      4
// 0: no cooling
// 1: brem
// 2: gammie thin disk simple
// 3: gammie two layer cooling
// 4: Kai cooling

#define NUMSPECTRUM 20  // number of frequency bands
#define NUMFUNCOOL 5120  // number of samples of cooling function
#define NUMFUNRELIE (NUMFUNCOOL) // number of samples of relativistic step_ie function

#define VISC_LINEAR 0 // 0-> linear viscosity off/1=on
#define VISC_TENSOR 0 // 0-> vnr visc 1-> tensor visc(not working)

// turn on/off different terms (1 or 0)
#define VISCE11 1
#define VISCE22 1
#define VISCE33 1
#define VISCE12 1
#define VISCE21 VISCE12
#define VISCE13 1
#define VISCE31 VISCE13
#define VISCE23 1
#define VISCE32 VISCE23

// MARK
#define ALFVENLIMIT 0
// 0: normal newtonian EOM
// 1: rho->rho+b^2*(sol^{-2})=rho+(b*invsol)^2

#define CSLIMIT 0
// 0: normal newtonian EOM
// 1: rho->rho+en*(sol^{-2})=rho+(en*invsol)^2

// 0: invsol2 set as constant
// 1: invsol2 set in stepgen.c by compute_invsol2fun()
#define DOINVSOL2FUN 0


#if(COMPDIM==3)
#define VOLUMEDIFF 1  // choice
// 0: use normal diffs or volume diffs, but normal diffs for divB=0 stuff
// 1: use volume diffs for everything such that divB=0 is conserved properly and no singularity issues in magnetic evolution
#else
#define VOLUMEDIFF 0 // no choice
#endif


#define TS0CHECK 0 // 1: check if rho or ie<0 in timestep.c  not needed if forcing floor
// number of dt checks for timestep.c
#if(SUBCOOL)
#define NUMDTCHECKS 10
#else
#define NUMDTCHECKS 11
#endif

// max number of CPUS for various allocations
// if want beyond this need to change CPUTXT on myid outputs in *.c
#define MAXCPUS 1000

// number of ghost zones
// should either use 1 for all non-zero or 2 for all non-zero, don't mix 1 and 2.
#define N1BND 2
#define N2BND 2
#if(COMPDIM==3)
#define N3BND 2
#else
#define N3BND 0
#endif

#define NBIGBND1 ((N1BND>N2BND) ? N1BND : N2BND)
#define NBIGBND  ((NBIGBND1>N3BND) ? NBIGBND1 : N3BND)

// N?OFF and N?NOT1 are a bit redundant
#define N1OFF (((N1BND>0)&&(N1>1)) ? 1 : 0)
#define N2OFF (((N2BND>0)&&(N2>1)) ? 1 : 0)
#define N3OFF (((N3BND>0)&&(N3>1)) ? 1 : 0)

#define N1NOT1 ((N1>1) ? 1 : 0)
#define N2NOT1 ((N2>1) ? 1 : 0)
#define N3NOT1 ((N3>1) ? 1 : 0)

/* allocated memory uses this for active zones 0-N1-1 and bc beyond that */
#define N1M (N1+N1BND*2)
#define N2M (N2+N2BND*2)
#define N3M (N3+N3BND*2)

/* NBIGM is bigger of N1M and N2M and N3M */
#define NBIG1M ((N1M>N2M) ? N1M : N2M)
#define NBIGM  ((NBIG1M>N3M) ? NBIG1M : N3M)

// maximal surface of boundary exchange
#if(COMPDIM==3)
#define NBIGS1M ((N1M*N2M>N1M*N3M) ? N1M*N2M : N1M*N3M)
#define NBIGSM ((NBIGS1M>N2M*N3M) ? NBIGS1M : N2M*N3M)
#else
#define NBIGSM (NBIGM)
#endif

#if(COORD==1)
#define LOWMEMMODE 1 // choice (good bit faster and alot less memory)
#else
#define LOWMEMMODE 0 // no choice (at least right now)
#endif


// GODMARK -- LOOPTYPE
#if(!GAMMIEIMAGE)
#define LOOPTYPE 1
// 1: old regular loop style
// 2: general loop over randomly organized grid (requires boundtype>1)
// 3: loop with outer/inner shell (requires boundtype>1)

// GODMARK -- BOUNDTYPE
#define BOUNDTYPE 1
// 1: general purpose rectangular grid(no holes or special features), although uses mask bits in semi-general way
// 2: general purpose cartesian-based bound code(e.g. sphere within cube, all outflow boundaries(with inflow check choosable)).  Can have complete swiss cheese and will work
// 3: general purpose cart-based, goes over interior boundary zones, but not exterior in comp loops (a compromise--timestep specially avoids interior bzones).  Can't have swiss-cheese, just inner/outer boundary zones for effectiveness.

#define PUREBC (0) // applies 
// 0: use general boundary condition
// 1-5: assume all boundaries are same and have type PUREBC

#else
// GAMMIE IMAGE SETTINGS
#define LOOPTYPE 1
#define BOUNDTYPE 1
#define PUREBC 0
#endif




#define BOUNDFIELD 0 // only applies if BOUNDTYPE>1
// 1: outflow field and no cleaning(bound() actually bounds field)
// 2: outflow emf and compute field from emf evolution(bound() does nothing to field, only emf before field calc from emfs)

// whether older loss diagnostic memory is initialized
#if(LOOPTYPE==1)
#define DOLOSSDIAGMEMORY (1) // choose whether want diagnostics (memory)
#else
#define DOLOSSDIAGMEMORY (0) // normally 0 unless you're stupid
#endif

// whether memoroy is initialized for averaging
#define DOAVGDIAGMEMORY (0)

// whether floor dump memory is initialized
#define FLOORDUMPFLAGMEMORY (0)


// how many seconds, below which if logging is taking place could impact performance and so is reported to user in the perf file
// comment shows what the time check is based upon
// the number inside is the # of seconds, the divisor is a scale factor(see main.c)
#define NUMLOGCHECKS 7 // number of these log checks

// really want walltime(to do step)/dt << walltime(to do dump)/Dt
#define LOGLWTCHECK (1.0/DTl) // 0
#define LOGDWTCHECK (20.0/DTd) // 1
#define LOGIWTCHECK (20.0/DTi) // 2
#define LOGENERWTCHECK (2.0/DTener) // 3
#define LOGLOSSWTCHECK (2.0/DTloss) // 4
#define LOGTIMESTEPWTCHECK (60.0/DTtimestep) // 5
#define LOGSPWTCHECK (60.0/DTsp) // 6


// see diag.c and init.c for more image and dump options.
#if(POSTPROC==0)
// 0=image raw rectangular grid data 1=point sample 2=plane interpolation
#define SAMPLEI 0 // should be 0 in general
#define SAMPLED 0 // should be 0 in general

#else

#define SAMPLEI 1
#define SAMPLED 1  // should choose 2 for pretty dumps, choose 1 for use as input data so no odd extrapolation at edges
#endif




#if(POSTPROC==0)

#if(SAMPLEI==0)
#define IMGN1 N1
#define IMGN2 N2
#define IMGN3 N3
#elif(SAMPLEI>0)
#define IMGN1 256 // size of image outputted in x-direction in image
#define IMGN2 512 // size of image outputted in y-direction in image
#define IMGN3 N3
#endif

#if(SAMPLED==0)
#define DUMN1 N1
#define DUMN2 N2
#define DUMN3 N3
#elif(SAMPLED>0)
#define DUMN1 256 // size of dump outputted in x-direction in image
#define DUMN2 512 // size of dump outputted in y-direction in image
#define DUMN3 N3
#endif

#else // POSTPROC==1

#if(SAMPLEI==0)
#define IMGN1 N1 // no choice
#define IMGN2 N2 // no choice
#define IMGN3 N3 // no choice
#else
#if(!PPCLEAN)
#define IMGN1 256 // choice
#define IMGN2 512 // choice
#define IMGN3 1 // no choice
#else
#define IMGN1 OUTPUTN1 // choice
#define IMGN2 OUTPUTN2 // choice
#define IMGN3 1 // no choice
#endif

#endif
#if(SAMPLED==0)
#define DUMN1 N1 // no choice
#define DUMN2 N2 // no choice
#define DUMN3 N3 // no choice
#else
#define DUMN1 256 // choice
#define DUMN2 512 // choice
#define DUMN3 N3 // no choice
#endif

#endif

#define ITYPES 2 // number of types of image ranges, 0, 1
#define CTYPES 2 // number of types of computed image ranges, 0, 1

// determine largest interpolated grid sizes for memory allocation of working space
#define INTN1 ((DUMN1>IMGN1) ? DUMN1 : IMGN1)
#define INTN2 ((DUMN2>IMGN2) ? DUMN2 : IMGN2)
#define INTNBIG ((INTN1>INTN2) ? INTN1 : INTN2)




#define NUMINDEX 29 // number of index situations per direction in super-reduced/memory optimized but not necessarily performance optimized situation
// 0: scalar total
// 1: LOOPFC
// 2: LOOPHC
// 3: LOOPFMHPC
// 4: LOOPHMFPC
// 5: LOOPC
// 6: scalar inner
// 7: scalar outer
// 8: v_1 inner
// 9: v_1 outer
// 10: v_2 inner
// 11: v_2 outer
// 12: v_3 inner
// 13: v_3 outer
// 14: B_1 inner
// 15: B_1 outer
// 16: B_2 inner
// 17: B_2 outer
// 18: B_3 inner
// 19: B_3 outer
// 20: LOOPB1
// 21: LOOPB2
// 22: LOOPB3
// 23: LOOPFLUXINNERX1
// 24: LOOPFLUXOUTERX1
// 25: LOOPFLUXINNERX2
// 26: LOOPFLUXOUTERX2
// 27: LOOPFLUXINNERX3
// 28: LOOPFLUXOUTERX3

#define IOBOUNDARY 0 // 0: normal boundary, 1: Inflow/outflow split outer x2-boundary
#define NUMBTRANS 4 // number of boundary transitions(see defs.h/init.c)

#define LINEARINTERP 1
// 0: just choose local value (for major speed, but large diffusion)
// 1: simple average or truely correct for uniform grid (for speed or uniform grid) (sorta corresponds to how differencing is done for 1st derivative in sweep.c)
// 2: correct linear average for generaly nonuniform grid (slower)

#define LINEXT 0
// 0: just copy
// 1: do linear extrap based with nonuniform ability, but not setup for 2 bzones

// major constants below
#define NUMSCA 3
#define NUMVEC 2
#define REALNUMVEC 2 // can use if not using MHD(then set to 1)
#define REALNUMBOUNDVEC ((BOUNDFIELD==2) ? (REALNUMVEC-1) : REALNUMVEC)
#define NUMSV (NUMSCA+NUMVEC)
#define NUMGRID 2
#define NUMMETRIC 4
#define DIRVEC 1
#define NUMFLOORVAR 2 // 1: mass density 2: ie density
#define NUMFLOOROUT 8 // number of routines checked for floor
// 0 through 7 with 1 extra space

// number of time(space) averaged quantities
#define NUMAVG2_3 15 // only averaged over time
#define NUMAVG1_32 15 // time and x2
#define NUMAVG1_31 15 // time and x1
// 1: rho
// 2: en
// 3: Be
// 4: cs^2
// 5: entropy (Exp[S])
// 6: vx1
// 7: vx2
// 8: vx3
// 9: sig11
// 10: sig12
// 11: sig13
// 12: sig22
// 13: sig23
// 14: sig33
// 15: nu_real(nuvisc)

#define NUMPPCALCS 12
// see diag.c(dump())

#define NUMWORKIQ ((NUMPPCALCS>NUMAVG2_3) ? NUMPPCALCS : NUMAVG2_3)

#define NUMLOSSVAR (NUMSCA+NUMVEC*(3+1)+1+3)
// 1: mass
// 2: enthalpy
// 3: grav pot energy

// 4: kinE

// 5: s1
// 6: s2
// 7: s3

// 8: magE

// 9: B1
// 10: B2
// 11: B3

// 12: etot(visc)

// below used for general coords
// 13: ang mom about x axis
// 14: ang mom about y axis
// 15: ang mom about z axis


#define INPUTPAR "%f"
#define INPUT2 "%f "
#define INPUT2B "%f "
#define INPUT3I "%f"
#define INPUTIMGT "%f"
// INPUT1/1old/4/5/6/7/avgh1 defines in i/i/diag.c/diag.c/diag.c/init.c/init.c respectivly because sensitive
#define INPUT7 "%f" // this not sensitive for now
#define INPUTRAD "%f"

#define ARGS 5
#define PARMTYPEARGS "%d %d %d %f %f"
#define MAXFILENAME 400

#define ZER 0x30
#define DIGILEN 10  

// if you change these, make sure you know where they are used!!
#define ERR     1.e-6
#define SMALL	1.e-10
#define GMIN 	1.e-10
#define MIN	1.e-20	/* minimum density */
#define SSMALL  1.E-30

#define THIRD (0.3333333333333333333333333333333333333333333333333333333333333333333)

#define HEADER3_S0 "%d %d\n"
#define HEADER3_S1 "%d %d %d\n"
#define HEADER3_S10 "%d %d %d %d\n"
#define HEADER3_S11 "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n"
#define HEADER3_S12 "%d %d %d %d\n"
#define HEADER3_S13 "%d %d %d %d %d %d %d %d %d %d\n"
#define HEADER3_S14 "%d %d %d %d %d %d %d %d %d %d %d %d\n"
#define HEADER3_S15 "%d %d %d %d %d %d\n"
#define HEADER3_S16 "%d %d %d %d %d %d %d %d %d %d\n"

// version of file format
// N1, N2, N3
 // L11, L12, L13, L21, L22, L23
 // rg rgp cs coolfact gam alpha_real0 n_real
 // nu_vnr nu_l nu_ten cour cour2
 // gravc massbh invsol2 blackholejz
 // tstart tf tavgi tavgf numavg
 // dtl dtd dti dtloss dtfloor dttimestep dtpd dtener dtfld
 // resist nu_sh
 // vg1 vg2 vg3
// coord fullvec analoutput DYNAMICMM
// trans,transx1,transrhox1,transiex1,transv1x1,transv2x1,transv3x1,transmagx1,transx2,transrhox2,transiex2,transv1x2,transv2x2,transv3x2,transmagx2,transx3,transrhox3,transiex3,transv1x3,transv2x3,transv3x3
// press,pressx1,pressx2,pressx3
// mag,transbzx,transbzy,stepmocct,mocctvx1,mocctvx2,mocctvx3,mocctbx1,mocctbx2,mocctbx3
// ie,visc_art,visc_real,vreal,vischeat,mdotin,cool,res_real,rreal,resheat,advint,kever
// intix1,intox1,intix2,intox2,intix3,intox3
// nonunigridx1,nonunigridx2,nonunigridx3,simplebc,bcix1,bcox1,bcix2,bcox2,bcix3,bxoc3
#define HEADER3_P  "# %4s %4s\n"\
                   "  %4d %4d\n"\
                   "# %4s %4s %4s\n"\
                   "  %4d %4d %d\n"\
                   "# %21s %21s %21s %21s %21s %21s\n"\
                   "  %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n"\
                   "# %21s %21s %21s %21s %21s %21s %21s\n"\
                   "  %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n"\
		   "# %21s %21s %21s %21s %21s\n"\
		   "  %21.15g %21.15g %21.15g %21.15g %21.15g\n"\
		   "# %21s %21s %21s %21s\n"\
		   "  %21.15g %21.15g %21.15g %21.15g\n"\
		   "# %21s %21s %21s %21s %21s %21s\n"\
		   "  %21.15g %21.15g %21.15g %21.15g %21d %21.15g\n"\
		   "# %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s\n"\
		   "  %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n"\
		   "# %21s %21s\n"\
		   "  %21.15g %21.15g\n"\
		   "# %21s %21s %21s\n"\
		   "  %21.15g %21.15g %21.15g\n"\
                   "# %21s %21s %21s %21s\n"\
                   "  %21d %21d %21d %21d\n"\
                   "# %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s\n"\
                   "  %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d\n"\
                   "# %21s %21s %21s %21s\n"\
                   "  %21d %21d %21d %21d\n"\
                   "# %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s\n"\
                   "  %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d\n"\
                   "# %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s\n"\
                   "  %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d\n"\
                   "# %21s %21s %21s %21s %21s %21s\n"\
                   "  %21d %21d %21d %21d %21d %21d\n"\
                   "# %21s %21s %21s %21s %21s %21s %21s %21s %21s %21s\n"\
                   "  %21d %21d %21d %21d %21d %21d %21d %21d %21d %21d\n"

#define HEADER4_P  " %4d %4d %4d %4d "\
                   "%21.15g %21.15g %21.15g "\
                   "%21.15g %21.15g %21.15g "\
                   "%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g "\
                   "%21.15g %21.15g %21.15g "

#define HEADER4_S  "%d %d %d %d "\
                   "%f %f %f "\
                   "%f %f %f "\
                   "%f %f %f %f %f %f %f %f "\
                   "%f %f %f "


/* Function defines */
// main.c
int main(
	 int argc,
	 char *argv[],
	 char *envp[]
	 );

extern void myargs(int argc, char *argv[]);


extern void gpar_init(void);
// utilfun.c
extern int myexit(int call_code);
extern void itoa(int x,char*p);
extern int mysys(char*com1, char*com2);
extern int fork2(void);
extern void ptraddr(int nstep);

// init.c
extern int init(int argc,
		char *argv[],
		char*envp[]);
extern void init_checks(void);
extern void init_genfiles(int gopp);
extern void init_general(void);
extern int  init_MPI(int argc,char *argv[]);
extern void  init_MPIgroup(void);
extern void  init_placeongrid(void);
extern void  init_optimalmpi(int which);
extern void jonio_init_fp(FILE **fp,int which,char *filename);
#if(USEMPI)
extern void jonio_init_mem(int numcolumns,int datatype,unsigned char **jonio1,float **jonio4, double **jonio8,unsigned char **writebuf1,float **writebuf4,double **writebuf8,unsigned char **tempbuf1, float **tempbuf4,double **tempbuf8);
extern void jonio_combine(int stage,MPI_Datatype mpidt,int numcolumns,int datatype,FILE* fp, void * jonio, void * writebuf, void * tempbuf);
extern void jonio_seperate(int stage,MPI_Datatype mpidt,int numcolumns,int datatype,FILE* fp, void * jonio, void * writebuf, void * tempbuf);
#endif
extern void init_loss_rect(void);
extern void init_loss_gen(void);
#if(LOOPTYPE==1)
#define init_loss init_loss_rect // called each timestep
#elif(LOOPTYPE>1)
#define init_loss init_loss_gen // called once at beginning
#endif
extern void init_compsave(void);
extern void init_visc(void);
extern void init_floor(void);
extern void init_inflows(void);
extern void init_radiations(void);
extern void init_res(void);
extern int init_mem(void);
extern int init_pointers(void);
extern int init_runpar(void);
extern int init_rundat(int dump_cnt,int which);
extern int init_runimage(int im_cnt,int wsca,int wvec,int call_code,int outtype);
#define PARAMETERPARLIST int seed, FTYPE beta, FTYPE nj
extern int init_paramstot(PARAMETERPARLIST);
extern int init_paramsphysics(PARAMETERPARLIST);
extern int init_paramsgeom(PARAMETERPARLIST);
extern int init_paramstime(PARAMETERPARLIST);
extern int init_paramsnumerical(PARAMETERPARLIST);
extern int init_paramspresets(PARAMETERPARLIST);
extern int init_reentrance(void);
extern int init_reentrance2(SFTYPE time, int which);
extern int init_otherparams(void);
extern int init_dx(int n1,int n2,int n3,int n1b,int n2b,int n3b, SFTYPE*sx, int which, int outtype);
extern int init_x(int n1,int n2,int n3,int n1b,int n2b,int n3b, SFTYPE*sx, int which, int outtype);
extern int init_reduceddimension(void);
extern int init_diffs(void);
extern int init_data(void);

#if(BOUNDTYPE==1)
#define init_bc init_bc_rect1
#elif((BOUNDTYPE==2)||(BOUNDTYPE==3))
#define init_bc init_bc_gen1
#endif
extern int init_bc_rect1(int simple,int ix1,int ox1,int ix2,int ox2,int ix3,int ox3);
extern int init_bc_gen1(int simple,int ix1,int ox1,int ix2,int ox2,int ix3,int ox3);
extern int init_mainbc(int px1,int six1,int rix1,int rox1,int px2,int six2,int rix2,int rox2,int px3,int six3,int rix3,int rox3);
extern int init_outgparm(int which);
extern FTYPE ranc(int iseed);

extern void accountstoreset(void);

extern void analsolve_ppc(int gopp);
extern void tdep_compute_ppc(void);
#if(!PPCLEAN)
#define analsolve analsolve_normal
#define tdep_compute tdep_compute_normal
#else
#define analsolve analsolve_ppc
#define tdep_compute tdep_compute_ppc
#endif

// diag.c
extern void diag(int call_code);
// the flux calculation is different enough in code to warrant different functions
extern void diag_loss_rect(int call_code,FTYPE tener,FTYPE tloss, SFTYPE* totloss_full);
extern void diag_loss_gen(int call_code,FTYPE tener,FTYPE tloss, SFTYPE* totloss_full);
extern void magnetic_flux_rect(void);
extern void magnetic_flux_gen(void);
extern void hydro_flux_rect(int dir,int which,FTYPE (*fl)[N3M][N2M][N1M]);
extern void hydro_flux_gen_adv(int dir,int which,FTYPE (*fl)[N3M][N2M][N1M]);
extern void hydro_flux_gen_gen(void);
extern void viscous_flux_rect(void);
extern void viscous_flux_gen(void);
#if(LOOPTYPE==1)
#define diag_loss diag_loss_rect
#define magnetic_flux magnetic_flux_rect
#define hydro_flux hydro_flux_rect
#define hydro_flux_adv hydro_flux_rect
#define viscous_flux viscous_flux_rect
#elif(LOOPTYPE>1)
#define diag_loss diag_loss_gen
#define magnetic_flux magnetic_flux_gen
#define hydro_flux_adv hydro_flux_gen_adv
#define hydro_flux hydro_flux_gen_gen
#define viscous_flux viscous_flux_gen
#endif
// ener and mode are fine in rect or gen loop/bound type
extern void diag_ener(int call_code,SFTYPE*totloss_full);
extern void diag_mode(int call_code);
extern void diag_dumps(int call_code,FTYPE tdump,FTYPE tfldump,FTYPE tpdump,FTYPE timage);
extern void diagavg(int call_code);
void divb0check(int which);
void symmetry_check(int pos);
void crazy_check(void);
void dump_header(FILE *fp,int which, int realsampled,int nogridchoice);
extern void dump_header2(FILE* fp, int which);
extern void dump(FILE *fp, int dump_cnt,int which,int outtype);
extern void image(int im_cnt, int wsca,int wvec, int call_code, int outtype);
extern FTYPE ffv_calcs(int type, int term, int i, int j, int k);
// numerics.c
extern void interpolate(int wtype, int wsv, int stype, int whichg, int inn1, int inn2, int outn1, int outn2, FTYPE in[N2M][N1M], FTYPE out[INTN2][INTN1], FTYPE outerdef,int outtype);
SFTYPE cart2spc(int which, SFTYPE xx, SFTYPE yy, SFTYPE zz);
void grids_cart2spc(FTYPE (*sca3) [N3M][N2M][N1M],FTYPE (*vx3) [N3M][N2M][N1M],FTYPE (*vy3) [N3M][N2M][N1M],FTYPE (*vz3) [N3M][N2M][N1M]);
extern FTYPE alnfact(FTYPE N);
extern void smooth(int SMOOTHSIZE,FTYPE BAD,FTYPE (*var)[N2M][N1M]);

extern int purebccall(int k, int j, int i,int*bcdim,int*bcdir);


#if(!PPCLEAN) // things not really used by postproc

#if(BOUNDTYPE==1)
#define bound bound_rect1
#elif((BOUNDTYPE==2)||(BOUNDTYPE==3))
#if(PUREBC==0) // GODMARK was flipped
#define bound bound_gen1
#elif(PUREBC==4)
#define bound bound_gensimple1
#else // no other special routines
#define bound bound_gen1
#endif
#endif


FTYPE e2zdiag(int whichvar, FTYPE (*var)[N2M][N1M], int k, int j, int i);

void dump1data(int which, int whichvar, int i, int j, int k, FILE* fileptr, FTYPE var);

void dump1databin(int which, int whichvar, int i, int j, int k, FTYPE *var, int datasize, int numdata, FILE* fileptr);

void writebuf1data(int which, int whichvar, int i, int j, int k, FTYPE var, FTYPE *writebuf);

void dumpcheckfprintf(int which, int whichvar, int i, int j, int k, FILE* fileptr, char *format, ...);


int inside_accretor(int which, int i, int j, int k);
int inside_accretor_general(int numzones, int which, int i, int j, int k);

extern void bound_accretor(int dir,FTYPE (*vars)[N2M][N1M],
		    FTYPE (*varv)[N3M][N2M][N1M],
		    int wsca,
		    int wvec,
		    int wcom);


extern void bound_rect1(FTYPE (*vars)[N2M][N1M],
		  FTYPE (*varv)[N3M][N2M][N1M],
		  int wsca,// which scalars to update: -1 for all or 1-NUMSCA for individual
		  int wvec,// which vectors to update: -1 for all or 1-NUMVEC for individual, and -2 for any scalar/vector
		  int wcom); // which component of vector, otherwise ignored for scalars
extern void bound_gen1(FTYPE (*vars)[N2M][N1M],
		  FTYPE (*varv)[N3M][N2M][N1M],
		  int wsca,// which scalars to update: -1 for all or 1-NUMSCA for individual
		  int wvec,// which vectors to update: -1 for all or 1-NUMVEC for individual, and -2 for any scalar/vector
		  int wcom); // which component of vector, otherwise ignored for scalars
extern void bound_gensimple1(FTYPE (*vars)[N2M][N1M],
		  FTYPE (*varv)[N3M][N2M][N1M],
		  int wsca,// which scalars to update: -1 for all or 1-NUMSCA for individual
		  int wvec,// which vectors to update: -1 for all or 1-NUMVEC for individual, and -2 for any scalar/vector
		  int wcom); // which component of vector, otherwise ignored for scalars
extern void bound_mpi(FTYPE (*vars)[N2M][N1M],
		  FTYPE (*varv)[N3M][N2M][N1M],
		  int wsca,// which scalars to update: -1 for all or 1-NUMSCA for individual
		  int wvec,// which vectors to update: -1 for all or 1-NUMVEC for individual, and -2 for any scalar/vector
		  int wcom); // which component of vector, otherwise ignored for scalars

extern void boundtest(int which);

// analsol.c
extern void analsolve_normal(int gopp);
extern void sodsol(int gopp);
extern void advsol(int gopp);
extern void gausssol(int gopp);
extern void bondisol(int gopp,FILE*analyticoutreal);
extern FILE * tori1sol(int gopp);
extern void injectsol(int gopp);
extern void visctermtest(int gopp);
extern void pulsesol(int gopp);
extern void test1sol(int gopp);
extern void test2sol(int gopp);
extern void magbreaksol(int gopp);
extern void magvortex(int gopp);
extern void magcorona(int gopp);
extern void chandran(int gopp);
extern void averystar(int gopp);


// stepgen.c stepmag.c step2d.c
void initialize_gravity(  FTYPE (*sca3)[N3M][N2M][N1M], FTYPE (*vx3)[N3M][N2M][N1M],  FTYPE (*vy3)[N3M][N2M][N1M],  FTYPE (*vz3)[N3M][N2M][N1M]);

extern void stepvar_2d(void);
extern void stepvar_3d(void);
#if(COORD==1)
#define compute_sigma_gen compute_sigma_1
#elif((COORD==3)||(COORD==2))
#define compute_sigma_gen compute_sigma_3
#endif
extern void compute_sigma_3(FTYPE (*sigma)[3][N3M][N2M][N1M],FTYPE (*rost)[3][N3M][N2M][N1M],FTYPE (*rostnu)[3][N3M][N2M][N1M],FTYPE (*nurho_real)[N2M][N1M],FTYPE (*delv)[N2M][N1M]);
extern void compute_sigma_1(FTYPE (*sigma)[3][N3M][N2M][N1M],FTYPE (*rost)[3][N3M][N2M][N1M],FTYPE (*rostnu)[3][N3M][N2M][N1M],FTYPE (*nurho_real)[N2M][N1M],FTYPE (*delv)[N2M][N1M]);
extern void timestep(void);
extern void idtcreate(FTYPE*idt,int k, int j, int i);
extern void timescale(void);
extern void timecheck(int failmode,FTYPE*idt,int k, int j, int i,int reall);
extern void magprepare(void) ;
extern void moc_ct_2d(void) ;
extern void moc_ct_3d_v1(void) ;
extern void moc_ct_3d_v2(void) ;
extern void lorentz_3d(void) ;
extern void step_bz_2d(void) ;
extern void bzsweepx_2d(void);
extern void bzsweepy_2d(void);
extern void step_pgc(void) ;
extern void step_visc(void) ;
extern void step_visc_real(void);
extern void tdep_compute_normal(void);
extern void sp_compute(void);
extern void injection_2d(void);
extern void injection_3d(void);
#if(COOL==0)
#define cooling cooling_nocooling
#elif(COOL==1)
#define cooling cooling_brem
#elif(COOL==2)
#define cooling cooling_simple1
#elif(COOL==3)
#define cooling cooling_twolayer
#elif(COOL==4)
#define cooling coolingKai
#endif
extern void cooling_nocooling(void);
extern void cooling_simple1(void);
extern void cooling_brem(void);
extern void cooling_twolayer(void);
extern void coolingKai(void);
extern void coollingthin(void);
extern void optdepth1(int initmode);
extern void optdepth2(int initmode);
extern void tempini(int initmode);
extern void compute_funcool(FTYPE *fun,FTYPE *thetai,int num,FTYPE thetas,FTYPE thetaf,FTYPE *dthetap);
extern void compute_funrelie(FTYPE *fun,FTYPE *thetai,int num,FTYPE thetas,FTYPE thetaf,FTYPE *dthetap);
extern void nu_compute(void);
extern void nu_res_compute(void);
extern void compute_invsol2fun(void);
extern void current_compute(int wcom);
extern void step_ie(void) ;
extern void step_relie (void);
extern void step_res(void) ;

extern void step_trans_2d(void) ;
extern void step_trans_3d(void) ;
extern void sweepx(void);
extern void sweepy(void);
extern void sweepz(void);
extern void dqx_calc(FTYPE (*var)[N2M][N1M], FTYPE (*dq)[N2M][N1M]) ;
extern void dqvx_calc(int wcom, FTYPE (*var)[N3M][N2M][N1M],FTYPE (*dqv)[N3M][N2M][N1M]);
extern void dqy_calc(FTYPE (*var)[N2M][N1M], FTYPE (*dq)[N2M][N1M]) ;
extern void dqvy_calc(int wcom,FTYPE (*var)[N3M][N2M][N1M],FTYPE (*dqv)[N3M][N2M][N1M]);
extern void dqz_calc(FTYPE (*var)[N2M][N1M], FTYPE (*dq)[N2M][N1M]) ;
extern void dqvz_calc(int wcom,FTYPE (*var)[N3M][N2M][N1M],FTYPE (*dqv)[N3M][N2M][N1M]);

void fracfloor_correct(void);
void floor_correct(int wsca,int wloc);

//extern FTYPE z2e_1(FTYPE (* var)[N1M],int j,int i);
//extern FTYPE z2e_2(FTYPE (*var)[N1M],int j,int i);
//extern FTYPE e2z_1(FTYPE (*var)[N1M],int j,int i);
//extern FTYPE e2z_2(FTYPE (*var)[N1M],int j,int i);
//extern FTYPE e2e_v2(FTYPE (*var)[N1M],int j, int i);
//extern FTYPE e2e_v1(FTYPE (*var)[N1M],int j, int i);
extern void ex_v(int dim,int dir,FTYPE (*var)[N1M], int k, int j, int i,int which);
extern void ex_s(int dim,int dir,FTYPE (*var)[N1M],int k, int j, int i, int which);
extern void ex_s_p(int dim,int dir,FTYPE (*var)[N1M],int k, int j, int i);
extern void ex_v_p(int dim,int dir,FTYPE (*var)[N1M],int k, int j, int i);

//cooling diffusion approximation
extern void fld(void);
extern void radstart(void);
extern void radend(void);
extern void initrad(void);
extern void derivs(void);
extern void rhs(void);
extern void riccg(FTYPE *epss,int *ks00,int *maxit);
extern double opac(double rho, double T);
#endif // end if !PPCLEAN

//reenter zz
extern void reentrancezz(void); 
extern FTYPE extre(FTYPE radius, FTYPE theta, int m);

#if(PPCLEAN)
#define bound bound_ppclean
#endif

extern void bound_ppclean(FTYPE (*vars)[N2M][N1M],
		  FTYPE (*varv)[N3M][N2M][N1M],
		  int wsca,
		  int wvec,
		  int wcom);



// Macros

// restrict loops only over relevant domain in reduced dimension case
#if(N1>1)
#define INFULL1 -N1BND
#define OUTFULL1 N1+N1BND
#define INHALF1 -N1BND/2
#define OUTHALF1 N1+N1BND/2
#define SHIFT1 1
#else
#define INFULL1 0
#define OUTFULL1 N1
#define INHALF1 0
#define OUTHALF1 N1
#define SHIFT1 0
#endif

#if(N2>1)
#define INFULL2 -N2BND
#define OUTFULL2 N2+N2BND
#define INHALF2 -N2BND/2
#define OUTHALF2 N2+N2BND/2
#define SHIFT2 1
#else
#define INFULL2 0
#define OUTFULL2 N2
#define INHALF2 0
#define OUTHALF2 N2
#define SHIFT2 0
#endif

#if(N3>1)
#define INFULL3 -N3BND
#define OUTFULL3 N3+N3BND
#define INHALF3 -N3BND/2
#define OUTHALF3 N3+N3BND/2
#define SHIFT3 1
#else
#define INFULL3 0
#define OUTFULL3 N3
#define INHALF3 0
#define OUTHALF3 N3
#define SHIFT3 0
#endif


// these loops used for general purposes
#define LOOPF3 for(k=INFULL3;k<OUTFULL3;k++)
#define LOOPF2 for(j=INFULL2;j<OUTFULL2;j++)
#define LOOPF1 for(i=INFULL1;i<OUTFULL1;i++)

// these loops used for general purposes
#define LOOPFP13 for(k=INFULL3+SHIFT3;k<OUTFULL3;k++)
#define LOOPFP12 for(j=INFULL2+SHIFT2;j<OUTFULL2;j++)
#define LOOPFP11 for(i=INFULL1+SHIFT1;i<OUTFULL1;i++)


#define LOOPH3 for(k=INHALF3;k<OUTHALF3;k++)
#define LOOPH2 for(j=INHALF2;j<OUTHALF2;j++)
#define LOOPH1 for(i=INHALF1;i<OUTHALF1;i++)

#define LOOPN3 for(k=0;k<N3;k++)
#define LOOPN2 for(j=0;j<N2;j++)
#define LOOPN1 for(i=0;i<N1;i++)

#define LOOPFMHP3 for(k=INFULL3;k<OUTHALF3;k++)
#define LOOPFMHP2 for(j=INFULL2;j<OUTHALF2;j++)
#define LOOPFMHP1 for(i=INFULL1;i<OUTHALF1;i++)

#define LOOPHMFP3 for(k=INHALF3;k<OUTFULL3;k++)
#define LOOPHMFP2 for(j=INHALF2;j<OUTFULL2;j++)
#define LOOPHMFP1 for(i=INHALF1;i<OUTFULL1;i++)

#define LOOPHP3 for(k=0;k<OUTHALF3;k++)
#define LOOPHP2 for(j=0;j<OUTHALF2;j++)
#define LOOPHP1 for(i=0;i<OUTHALF1;i++)

// below used for initialization and such, not a computational issue
#define LOOPF LOOPF3 LOOPF2 LOOPF1
#define LOOPH LOOPH3 LOOPH2 LOOPH1
#define LOOP LOOPN3 LOOPN2 LOOPN1
#define LOOPFMHP LOOPFMHP3 LOOPFMHP2 LOOPFMHP1
#define LOOPHMFP LOOPHMFP3 LOOPHMFP2 LOOPHMFP1
#define LOOPHP LOOPHP3 LOOPHP2 LOOPHP1

#define LOOPINT3 for(k=intix3;k<intox3;k++)
#define LOOPINT2 for(j=intix2;j<intox2;j++)
#define LOOPINT1 for(i=intix1;i<intox1;i++)

#define LOOPDIAGOUTPUT(bzones3,bzones2,bzones1)\
for(k=((-bzones3<INFULL3) ? INFULL3 : -bzones3);k<((N3+bzones3>OUTFULL3) ? OUTFULL3 : N3+bzones3);k++)\
for(j=((-bzones2<INFULL2) ? INFULL2 : -bzones2);j<((N2+bzones2>OUTFULL2) ? OUTFULL2 : N2+bzones2);j++)\
for(i=((-bzones1<INFULL1) ? INFULL1 : -bzones1);i<((N1+bzones1>OUTFULL1) ? OUTFULL1 : N1+bzones1);i++)


#define LOOPLIMIT(bzones3,bzones2,bzones1) LOOPDIAGOUTPUT(-bzones3,-bzones2,-bzones1)

#define LOOPDIAGOUTPUTFULL(bzones3,bzones2,bzones1)\
for(k=((-bzones3<INFULL3) ? INFULL3 : -bzones3);k<((totalsize[3]+bzones3>OUTFULL3) ? OUTFULL3 : totalsize[3]+bzones3);k++)\
for(j=((-bzones2<INFULL2) ? INFULL2 : -bzones2);j<((totalsize[2]+bzones2>OUTFULL2) ? OUTFULL2 : totalsize[2]+bzones2);j++)\
for(i=((-bzones1<INFULL1) ? INFULL1 : -bzones1);i<((totalsize[1]+bzones1>OUTFULL1) ? OUTFULL1 : totalsize[1]+bzones1);i++)


//////////////////////////////////
//////////////////////////////////
//
// LOOPTYPE==1

#if(LOOPTYPE==1)

#define LOOPSUPERGEN(which) LOOPF // dummy
#define LOOPSUPERGEN2(which) LOOPF // dummy

#define LOOPFC LOOPF
#define LOOPHC LOOPH
#define LOOPFMHPC LOOPFMHP
#define LOOPHMFPC LOOPHMFP
#define LOOPHPC LOOPHP


#define LOOPC3 LOOPN3
#define LOOPC2 LOOPN2
#define LOOPC1 LOOPN1

#define LOOPC LOOPC3 LOOPC2 LOOPC1


#define LOOPTIMESTEP LOOPC
//#define LOOPDIVB LOOPFMHPC
#define LOOPDIVB LOOPC // GODMARK until can figure out why cleaning not working

#define LOOPINT LOOPINT3 LOOPINT2 LOOPINT1

#define LOOPFLUXINNER LOOPINT // not used
#define LOOPFLUXOUTER LOOPINT // not used

// must avoid non-comp zones (assumes bzones are assigned correctly)
#define LOOPFLOOR LOOPC


#define LOOPBOUND(which) LOOPHC // not used
#define LOOPBOUNDGEN LOOPHC
#define LOOPBOUNDV1 LOOPHC
#define LOOPBOUNDV2 LOOPHC
#define LOOPBOUNDV3 LOOPHC

// works when skipping bound if skip outer and not inner loop with BOUNDTYPE/LOOPTYPE==3
#define LOOPB1 LOOPHPC
#define LOOPB2 LOOPHPC
#define LOOPB3 LOOPHPC


#define LOOPSK3 for(k=skipix3;k<N3;k++)
#define LOOPSK2 for(j=skipix2;j<N2+periodicx2special;j++) // needed and no problem where used
#define LOOPSK1 for(i=skipix1;i<N1;i++)

// used for velocity calculations in step.c since when reflecting r/theta your metric components are 0 and floating point exception will occur, so just let bound do the work
// can either use below 2 loops to seperate x1/x2-directions or use LOOP above then have conditions using defines/ifs in loop
#define LOOPV1  LOOPC3 LOOPC2 LOOPSK1
#define LOOPV2  LOOPC3 LOOPSK2 LOOPC1
#define LOOPV3  LOOPSK3 LOOPC2 LOOPC1

#define LOOPINJ for(k=tagik;k<tagfk;k++) for(j=tagij;j<tagfj;j++) for(i=tagii;i<tagfi;i++)
#define LOOPFINJ for(k=t2gik;k<t2gfk;k++) for(j=t2gij;j<t2gfj;j++) for(i=t2gii;i<t2gfi;i++)
#define LOOPRINJ for(k=t3gik;k<t3gfk;k++) for(j=t3gij;j<t3gfj;j++) for(i=t3gii;i<t3gfi;i++)
#define LOOPVINJ for(k=t4gik;k<t4gfk;k++) for(j=t4gij;j<t4gfj;j++) for(i=t4gii;i<t4gfi;i++)



#define LOOP2  LOOPC2 LOOPC1


// sweepx

// used for zone centered stuff
#define LOOPT1i LOOPC3 LOOPC2 for(i=0;i<N1+N1OFF;i++)

// used for mdot[1]
#define LOOPT0i for(k=-N3OFF+skipix3;k<N3;k++) for(j=-N2OFF+skipix2;j<N2;j++) for(i=-N1OFF+skipix1;i<N1+N1OFF;i++)

// unless SKIPIX1==1, so that vx(0) is boundary zone so don't need fl(-1) for vx advection in x1-dir
// used for sweepx(vx)
#define LOOPT2i LOOPC3 LOOPC2 for(i=skipix1-N1OFF;i<N1;i++) 
// used for fluxes which never have i=-1 accessed, and don't need j=-1 or 0 accessed for non-periodic bcs
// below used on sweepx(vy)
#define LOOPT3i LOOPC3 for(j=skipix2;j<N2;j++) for(i=0;i<N1+N1OFF;i++)
// below used on sweepx(vz)
#define LOOPT4i for(k=skipix3;k<N3;k++) LOOPC2 for(i=0;i<N1+N1OFF;i++)


// sweepy:

// used for fluxes which never have i=-1 accessed
#define LOOPT1j LOOPC3 for(j=0;j<N2+N2OFF;j++) LOOPC1

// mdot(2)(j=-1) determined ok since innermost vector exists vy(-2), unlike vy(N2+2), but bound anyways.  This bounding is ok since mdot(2)(j=N+1) used for Fl(2)(N2), and mdot(2)(j=N+1) requires rho(N2+2), but THIS mdot not on boundary so easily defined by nonlocal copy of mdot(2)(j=N2-1) across in phi
#define LOOPT0j for(k=-N3OFF+skipix3;k<N3;k++) for(j=-N2OFF+skipix2;j<N2+N2OFF;j++) for(i=-N1OFF+skipix1;i<N1;i++)

// sweepy(vy)
// periodicx2special=1 makes set j=N2, which is ok if vy(j=N2+2) is taken care of properly.  Instead of using special code or whatever, I just bound flux, so no need to get fl(N2) which uses dq(N2+1) which uses vy(N2+2)
// no problem with inner edge since vy(j=-2) exists. assymetry due to staggered grid.
#define LOOPT2j LOOPC3 for(j=-N2OFF+skipix2;j<N2;j++) LOOPC1

// used for fluxes which never have i=-1 accessed
// below used on sweepy(vx)
#define LOOPT3j LOOPC3 for(j=0;j<N2+N2OFF;j++) for(i=skipix1;i<N1;i++)
// below used on sweepy(vz)
#define LOOPT4j for(k=skipix3;k<N3;k++) for(j=0;j<N2+N2OFF;j++) LOOPC1



// sweepz

// used for fluxes which never have i=-1 accessed
#define LOOPT1k for(k=0;k<N3+N3OFF;k++) LOOPC2 LOOPC1

#define LOOPT0k for(k=-N3OFF+skipix3;k<N3+N3OFF;k++) for(j=-N2OFF+skipix2;j<N2;j++) for(i=-N1OFF+skipix1;i<N1;i++)

// used for fluxes which never have i=N3-1+1 accessed
// below used on sweepz(vz)
#define LOOPT2k for(k=-N3OFF+skipix3;k<N3;k++) LOOPC2 LOOPC1

// used for fluxes which never have i=-1 accessed
// below used on sweepz(vx)
#define LOOPT3k for(k=0;k<N3+N3OFF;k++) LOOPC2 for(i=skipix1;i<N1;i++)

// below used on sweepz(vy)
#define LOOPT4k for(k=0;k<N3+N3OFF;k++) for(j=skipix2;j<N2;j++) LOOPC1

// loops over more i for j and more j for i than needed, but faster to do all in one lump
#define LOOPVISC for(k=-N3OFF;k<N3;k++) for(j=-N2OFF;j<N2+periodicx2special;j++) for(i=-N1OFF;i<N1;i++) 





////////////////////////////////////
////////////////////////////////////
// LOOPS FOR BOUNDTYPE==2




#elif((LOOPTYPE==2)||(LOOPTYPE==3))

#define LOOPSUPERGEN(which) for(temptempi=0,i=indx[which][0],j=indx[which][1],k=indx[which][2];temptempi<numiter[which];temptempi++,i=indx[which][temptempi*3],j=indx[which][temptempi*3+1],k=indx[which][temptempi*3+2])

#define LOOPSUPERGEN2(which) LOOPF3 LOOPF2 for(i=iindx[which][0][k][j];i<=iindx[which][1][k][j];i++)


// LOOPFC is used for computations, LOOPF for initializations, etc.

// optimized loops
#define LOOPBOUND(which) LOOPSUPERGEN(which)
#define LOOPBOUNDSCAin LOOPSUPERGEN(6) // only loop where loopvar is actually used directly(temptempi), see boundgen1.c
#define LOOPBOUNDSCAout LOOPSUPERGEN(7) // only loop where loopvar is actually used directly(temptempi), see boundgen1.c
#define LOOPBOUNDV1in LOOPSUPERGEN(8) // only loop where loopvar is actually used directly(temptempi), see boundgen1.c
#define LOOPBOUNDV1out LOOPSUPERGEN(9) // only loop where loopvar is actually used directly(temptempi), see boundgen1.c
#define LOOPBOUNDV2in LOOPSUPERGEN(10) // only loop where loopvar is actually used directly(temptempi), see boundgen1.c
#define LOOPBOUNDV2out LOOPSUPERGEN(11) // only loop where loopvar is actually used directly(temptempi), see boundgen1.c
#define LOOPBOUNDV3in LOOPSUPERGEN(12) // only loop where loopvar is actually used directly(temptempi), see boundgen1.c
#define LOOPBOUNDV3out LOOPSUPERGEN(13) // only loop where loopvar is actually used directly(temptempi), see boundgen1.c


#if(LOOPTYPE==2)
#define LOOPFC LOOPSUPERGEN(1)
#define LOOPHC LOOPSUPERGEN(2)
#define LOOPFMHPC LOOPSUPERGEN(3)
#define LOOPHMFPC LOOPSUPERGEN(4)
#define LOOPC LOOPSUPERGEN(5)

#define LOOPMAG1  LOOPSUPERGEN(20)
#define LOOPMAG2  LOOPSUPERGEN(21)
#define LOOPMAG3  LOOPSUPERGEN(22)

#elif(LOOPTYPE==3)
#define LOOPFC LOOPSUPERGEN2(1)
#define LOOPHC LOOPSUPERGEN2(2)
#define LOOPFMHPC LOOPSUPERGEN2(3)
#define LOOPHMFPC LOOPSUPERGEN2(4)
#define LOOPC LOOPSUPERGEN2(5)

#define LOOPMAG1  LOOPSUPERGEN2(20)
#define LOOPMAG2  LOOPSUPERGEN2(21)
#define LOOPMAG3  LOOPSUPERGEN2(22)

// note that while require to not over compute field on outside if using skip boundary method, ok to overcompute inner boundary zones since assign bzones each timestep.  In general, the precision required to only compute certain zones that are valid is true because of this skipping method...otherwise boundary zones change with bad info and don't get assigned due to skipping.

// but for inner zones, which are assigned every timestep always, ok to overdo, and this is the case with the LOOPTYPE==3 method which computes the interior boundary zones but avoids the exterior boundary zones...so this is ok.

#endif



//#define LOOPBOUNDGEN LOOPF
//#define LOOPFC LOOPF
//#define LOOPHC LOOPH
//#define LOOPFMHPC LOOPFMHP
//#define LOOPHMFPC LOOPHMFP
//#define LOOPC LOOP

// may want to try LOOPTIMESTEP with LOOPSUPERGEN(5) since seems to deal with that loop type well
//#define LOOPTIMESTEP LOOPC
#define LOOPTIMESTEP LOOPSUPERGEN(5) // since must avoid crazy regions (or did I handle this another way?)

#if(BOUNDTYPE==1)
#define LOOPDIVB LOOPFMHPC
#else
#if(BOUNDFIELD==1)
#define LOOPDIVB LOOPSUPERGEN(5) // required since need exactly the right loop for diagnostic
#else
#define LOOPDIVB LOOPSUPERGEN(3) // as above
#endif
#endif

#if(BOUNDFIELD==1)
#define LOOPB1  LOOPMAG1
#define LOOPB2  LOOPMAG2
#define LOOPB3  LOOPMAG3

#elif(BOUNDFIELD==2)
#define LOOPB1  LOOPFMHPC
#define LOOPB2  LOOPFMHPC
#define LOOPB3  LOOPFMHPC
#endif

// goes over velocity boundary zone, but ok
#define LOOPV1 LOOPC
#define LOOPV2 LOOPC
#define LOOPV3 LOOPC

// others just expand a bit.  slightly more zones to do, but vastly simpler and less memory used for index array(which itself gives vastly more performance)

#define LOOPHPC LOOPHC
#define LOOPHMC LOOPHC
#define LOOPINT LOOPSUPERGEN(5) // must avoid non-real comp. regions

 //LOOPSUPERGEN(1) // go over entire grid for security
#define LOOPFLOOR LOOPC


// INJ not setup yet (make #6, with largest range for all the below)
#define LOOPINJ for(k=tagik;k<tagfk;k++) for(j=tagij;j<tagfj;j++) for(i=tagii;i<tagfi;i++)
#define LOOPFINJ for(k=t2gik;k<t2gfk;k++) for(j=t2gij;j<t2gfj;j++) for(i=t2gii;i<t2gfi;i++)
#define LOOPRINJ for(k=t3gik;k<t3gfk;k++) for(j=t3gij;j<t3gfj;j++) for(i=t3gii;i<t3gfi;i++)
#define LOOPVINJ for(k=t4gik;k<t4gfk;k++) for(j=t4gij;j<t4gfj;j++) for(i=t4gii;i<t4gfi;i++)

#define LOOP2  LOOPC2 LOOPC1


// sweepx

// used for zone centered stuff
#define LOOPT1i LOOPHC

// used for mdot[1]
#define LOOPT0i LOOPHC

// unless SKIPIX1==1, so that vx(0) is boundary zone so don't need fl(-1) for vx advection in x1-dir
// used for sweepx(vx)
#define LOOPT2i LOOPHC
// used for fluxes which never have i=-1 accessed, and don't need j=-1 or 0 accessed for non-periodic bcs
// below used on sweepx(vy)
#define LOOPT3i LOOPHC
// below used on sweepx(vz)
#define LOOPT4i LOOPHC


// sweepy:

// used for fluxes which never have i=-1 accessed
#define LOOPT1j LOOPHC

// mdot(2)(j=-1) determined ok since innermost vector exists vy(-2), unlike vy(N2+2), but bound anyways.  This bounding is ok since mdot(2)(j=N+1) used for Fl(2)(N2), and mdot(2)(j=N+1) requires rho(N2+2), but THIS mdot not on boundary so easily defined by nonlocal copy of mdot(2)(j=N2-1) across in phi
#define LOOPT0j LOOPHC

// sweepy(vy)
// periodicx2special=1 makes set j=N2, which is ok if vy(j=N2+2) is taken care of properly.  Instead of using special code or whatever, I just bound flux, so no need to get fl(N2) which uses dq(N2+1) which uses vy(N2+2)
// no problem with inner edge since vy(j=-2) exists. assymetry due to staggered grid.
#define LOOPT2j LOOPHC

// used for fluxes which never have i=-1 accessed
// below used on sweepy(vx)
#define LOOPT3j LOOPHC
// below used on sweepy(vz)
#define LOOPT4j LOOPHC



// sweepz

// used for fluxes which never have i=-1 accessed
#define LOOPT1k LOOPHC

#define LOOPT0k LOOPHC

// used for fluxes which never have i=N3-1+1 accessed
// below used on sweepz(vz)
#define LOOPT2k LOOPHC

// used for fluxes which never have i=-1 accessed
// below used on sweepz(vx)
#define LOOPT3k LOOPHC

// below used on sweepz(vy)
#define LOOPT4k LOOPHC

// loops over more i for j and more j for i than needed, but faster to do all in one lump
#define LOOPVISC LOOPHC

#endif // end if looptype==2||3


// use for starting grid in upper left corner
#define LOOPI for(j=0;j<IMGN2;j++) for(i=0;i<IMGN1;i++)
#define LOOPICART #define LOOPI for(j=IMGN2-1;j>=0;j--) for(i=0;i<IMGN1;i++)
#define LOOPINI for(j=0;j<N2;j++) for(i=0;i<N1;i++) // used for post process 
#define LOOPD for(j=0;j<DUMN2;j++) for(i=0;i<DUMN1;i++)



#define DEBUGP1 \
printf("\n\n\n\n"); \
  LOOPF{ \
printf("j: %2d i: %2d x1a: %21.15g x1b: %21.15g x2a: %21.15g x2b: %21.15g rho: %20.15g e: %21.15g pot: %21.15g Bx: %21.15g By: %21.15g Bz: %21.15g vx: %20.15g vy: %21.15g vz: %21.15g\n",j,i,x[1][1][i],x[2][1][i],x[1][2][j],x[2][2][j],s[1][k][j][i],s[2][k][j][i],s[3][k][j][i],v[2][1][k][j][i],v[2][2][k][j][i],v[2][3][k][j][i],v[1][1][k][j][i],v[1][2][k][j][i],v[1][3][k][j][i]); \
  }


#define BOUNDBUG printf("k: %d j: %d i: %d ll: %d kk: %d jj: %d ii: %d\n",k,j,i,ll,kk,jj,ii); printf("bcd1: %d bcd2: %d bcdir: %d\n",bcd1,bcd2,bcdir);



#define NP      8       /* number of primitive variables */

#define RHO     0       /* mnemonics for primitive vars */
#define UX      1
#define UY      2
#define UZ      3
#define BX      4
#define BY      5
#define BZ      6
#define UU      7 
#define POT     8 // aux, doesn't need new space or anything
#define LOOPP for(k=0;k<NP;k++)
#define LOOPTVDLF for(i=0;i<N1;i++)for(j=0;j<N2;j++)
#if(!PPCLEAN)
extern void zeus2tvdlf(void);
extern void tvdlf2zeus(void);
extern void steptvdlf(void);
extern void step_ch(void);
#endif
extern void init_tvdlfgrid(void);

// setup for various boundary situations
// so doesn't produce differences in irrelevant directions, whether boundary zones or not
// mac(?) macros are for use in definitions of other macros since macro with args needs to be directly a function of what's hardcoded, not some "replacement" since nothing is replaced, not to be used inside code for USING a macro.

#if(N1>1)
#define im1 i-1
#define im1mac(i) i-1
#define ip1 i+1
#define ip1mac(i) i+1
#else

#define im1 i
#define im1mac(i) i
#define ip1 i
#define ip1mac(i) i

#endif

#if(N2>1)
#define jm1 j-1
#define jm1mac(j) j-1
#define jp1 j+1
#define jp1mac(j) j+1
#else
#define jm1 j
#define jm1mac(j) j
#define jp1 j
#define jp1mac(j) j
#endif

#if(N3>1)
#define km1 k-1
#define km1mac(k) k-1
#define kp1 k+1
#define kp1mac(k) k+1
#else
#define km1 k
#define km1mac(k) k
#define kp1 k
#define kp1mac(k) k
#endif


// seperate for now (combine loops later, but make new interp for 1d/2d/3d
#if(COMPDIM==1)
#include "global2dsup.h" // no optimized version yet
#elif(COMPDIM==2)
#include "global2dsup.h"
#elif(COMPDIM==3)
#include "global3dsup.h"
#endif


// nifty macros (which is position relative to grid)
#define b1(which,name,k,j,i) ((which==-1) ? b1m(name,k,j,i) : b1p(name,k,j,i))
#define b2(which,name,k,j,i) ((which==-1) ? b2m(name,k,j,i) : b2p(name,k,j,i))
#define b3(which,name,k,j,i) ((which==-1) ? b3m(name,k,j,i) : b3p(name,k,j,i))

#define bdivb(wcom,which,name,k,j,i) ( (wcom==1) ? b1(which,name,k,j,i) : ( (wcom==2) ? b2(which,name,k,j,i) : b3(which,name,k,j,i) ) )

#define BUFFERMAP ((k*N2*N1+j*N1+i)*numcolumns+nextbuf++)
#define BUFFERINIT if(mpicombine==1) nextbuf=0


#if(LOWMEMMODE==1) // only valid for cart coords

#define ODX(grid,dir,index) (1.0/dx[grid][dir][index])

#define DVL(grid,dir,index) (dx[grid][dir][index])

#define ODVL(grid,dir,index) (ODX(grid,dir,index))

// OVOL<gridx3><gridx2><gridx1> (dx3 external)
#define OVOL1(k,j,i) (1.0/(dx[1][1][i]*dx[1][2][j]))
#define OVOL2(k,j,i) (1.0/(dx[2][1][i]*dx[1][2][j]))
#define OVOL3(k,j,i) (1.0/(dx[1][1][i]*dx[2][2][j]))
#define OVOL4(k,j,i) (1.0/(dx[2][1][i]*dx[2][2][j]))

// VOL<gridx3><gridx2><gridx1> (dx3 external)
#define VOL1(k,j,i) (dx[1][1][i]*dx[1][2][j])
#define VOL2(k,j,i) (dx[2][1][i]*dx[1][2][j])
#define VOL3(k,j,i) (dx[1][1][i]*dx[2][2][j])
#define VOL4(k,j,i) (dx[2][1][i]*dx[2][2][j])



// note this doesn't include dx3, computed directly always
// numbers correpond to memory DS order (see init_compsave())
// DSF<center><dir>
//<center>=0.5,0.5 (i.e. rho)
#define DS11(k,j,i) (dx[1][2][j]) // *dx[2][3][k]
#define DS12(k,j,i) (dx[1][1][i]) // *dx[2][3][k]
#define DS13(k,j,i) (dx[1][1][i]*dx[1][2][j])
//<center>=0,0.5 (i.e. vx)
#define DS21(k,j,i) (dx[1][2][j]) // *dx[2][3][k]
#define DS22(k,j,i) (dx[2][1][i]) // *dx[2][3][k]
#define DS23(k,j,i) (dx[2][1][i]*dx[1][2][j])
//<center>=0.5,0 (i.e. vy)
#define DS31(k,j,i) (dx[2][2][j]) // *dx[2][3][k]
#define DS32(k,j,i) (dx[1][1][i]) // *dx[2][3][k]
#define DS33(k,j,i) (dx[1][1][i]*dx[2][2][j])
//<center>=0,0 (i.e. emf)
#define DS41(k,j,i) (dx[2][2][j]) // *dx[2][3][k]
#define DS42(k,j,i) (dx[2][1][i]) // *dx[2][3][k]
#define DS43(k,j,i) (dx[2][1][i]*dx[2][2][j])

// OARCL<center><dir>
//<center>=0.5,0.5 (i.e. rho)
#define OARC11(k,j,i) (ODX(2,1,i))
#define OARC12(k,j,i) (ODX(2,2,j))
#define OARC13(k,j,i) (1)
//<center>=0,0.5 (i.e. vx)
#define OARC21(k,j,i) (ODX(1,1,i))
#define OARC22(k,j,i) (ODX(2,2,j))
#define OARC23(k,j,i) (1)
//<center>=0.5,0 (i.e. vy)
#define OARC31(k,j,i) (ODX(2,1,i))
#define OARC32(k,j,i) (ODX(1,2,j))
#define OARC33(k,j,i) (1)
//<center>=0,0 (i.e. emf)
#define OARC41(k,j,i) (ODX(2,1,i))
#define OARC42(k,j,i) (ODX(2,2,j))
#define OARC43(k,j,i) (1)

#define G2(grid,i) (1)
#define G3(grid,i) (1)
#define G4(grid,j) (1)

#define DG2(grid,i) (0)
#define DG3(grid,i) (0)
#define DG4(grid,j) (0)

#define OG2(grid,i) (1)
#define OG3(grid,i) (1)
#define OG4(grid,j) (1)



#else // use normal memory

//#define ODX(grid,dir,index) (1.0/dx[grid][dir][index])
#define ODX(grid,dir,index) (odx[grid][dir][index])

#define DVL(grid,dir,index) (dvl[grid][dir][index])

#define ODVL(grid,dir,index) (odvl[grid][dir][index])

// OVOL<gridx3><gridx2><gridx1> (dx3 external)
#define OVOL1(k,j,i) (ovol[1][k][j][i])
#define OVOL2(k,j,i) (ovol[2][k][j][i])
#define OVOL3(k,j,i) (ovol[3][k][j][i])
#define OVOL4(k,j,i) (ovol[4][k][j][i])

// VOL<gridx3><gridx2><gridx1> (dx3 external)
// rarely used
#define VOL1(k,j,i) (1.0/ovol[1][k][j][i])
#define VOL2(k,j,i) (1.0/ovol[2][k][j][i])
#define VOL3(k,j,i) (1.0/ovol[3][k][j][i])
#define VOL4(k,j,i) (1.0/ovol[4][k][j][i])



// note this doesn't include dx3, computed directly always
// numbers correpond to memory DS order (see init_compsave())
// DSF<center><dir>
//<center>=0.5,0.5 (i.e. rho)
#define DS11(k,j,i) (ds[1][1][k][j][i])
#define DS12(k,j,i) (ds[1][2][k][j][i])
#define DS13(k,j,i) (ds[1][3][k][j][i])
//<center>=0,0.5 (i.e. vx)
#define DS21(k,j,i) (ds[2][1][k][j][i])
#define DS22(k,j,i) (ds[2][2][k][j][i])
#define DS23(k,j,i) (ds[2][3][k][j][i])
//<center>=0.5,0 (i.e. vy)
#define DS31(k,j,i) (ds[3][1][k][j][i])
#define DS32(k,j,i) (ds[3][2][k][j][i])
#define DS33(k,j,i) (ds[3][3][k][j][i])
//<center>=0,0 (i.e. emf)
#define DS41(k,j,i) (ds[4][1][k][j][i])
#define DS42(k,j,i) (ds[4][2][k][j][i])
#define DS43(k,j,i) (ds[4][3][k][j][i])

// OARCL<center><dir>
//<center>=0.5,0.5 (i.e. rho)
#define OARC11(k,j,i) (oarcl[1][1][k][j][i])
#define OARC12(k,j,i) (oarcl[1][2][k][j][i])
#define OARC13(k,j,i) (oarcl[1][3][k][j][i])
//<center>=0,0.5 (i.e. vx)
#define OARC21(k,j,i) (oarcl[2][1][k][j][i])
#define OARC22(k,j,i) (oarcl[2][2][k][j][i])
#define OARC23(k,j,i) (oarcl[2][3][k][j][i])
//<center>=0.5,0 (i.e. vy)
#define OARC31(k,j,i) (oarcl[3][1][k][j][i])
#define OARC32(k,j,i) (oarcl[3][2][k][j][i])
#define OARC33(k,j,i) (oarcl[3][3][k][j][i])
//<center>=0,0 (i.e. emf)
#define OARC41(k,j,i) (oarcl[4][1][k][j][i])
#define OARC42(k,j,i) (oarcl[4][2][k][j][i])
#define OARC43(k,j,i) (oarcl[4][3][k][j][i])

#define G2(grid,i) (g[grid][2][i])
#define G3(grid,i) (g[grid][3][i])
#define G4(grid,j) (g[grid][4][j])

#define DG2(grid,i) (dg[grid][2][i])
#define DG3(grid,i) (dg[grid][3][i])
#define DG4(grid,j) (dg[grid][4][j])

#define OG2(grid,i) (og[grid][2][i])
#define OG3(grid,i) (og[grid][3][i])
#define OG4(grid,j) (og[grid][4][j])



#endif
// diffusion approximation
// define LOOP used in radiation module :one more grid beyond active zone
#define LOOPRC1 LOOPH1
#define LOOPRC2 LOOPH2
#define LOOPRC3 LOOPH3
#define LOOPRC LOOPH
//define LOOP used in fld.src
#define LOOPFLD1C for(i=0;i<N1+1;i++)
#define LOOPFLD2C for(j=0;j<N2+1;j++)
// define LOOPDIVJ used in gradv.src
#define LOOPDIVJ for(j=0;j<N2+1;j++)
#define RAD 0
#define RADTEST 0
#define RT 0
#define KAICOOL 1
#define REENTER 1
#define THINCOOL 1
#define DURISEN 1
#define LIMITER 0
#define DIAGF 1
