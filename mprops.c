/* mprops.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/
#include "global.h"
#include "f2c.h"
#if(LINUXCLUSTER==0)
SFTYPE munit,lunit,rhounit,tunit,eneunit,enedenunit;
SFTYPE lsuncgs,kcontcgs,mhydrcgs,sigmacontcgs,sunmasscgs,yrcgs,mstarcgs;
SFTYPE sigmacontnew,kcontnew,mhydrnew;
SFTYPE rsun,mmw;

#else
extern SFTYPE munit,lunit,rhounit,tunit,eneunit,enedenunit;
extern SFTYPE lsuncgs,kcontcgs,mhydrcgs,sigmacontcgs,sunmasscgs,yrcgs,mstarcgs;
extern SFTYPE sigmacontnew,kcontnew,mhydrnew;
extern SFTYPE rsun,mmw;

#endif

/* ======================================================================= */
/* ////////////////////////////  FILE MPROPS  \\\\\\\\\\\\\\\\\\\\\\\\\\\\ */

/*  PURPOSE:  This file contains subroutines which compute the material */
/*  properties needed in the radiation hydrodynamics.  Each routine is */
/*  independent, requires no common blocks (all data passed via */
/*  arguments), and works on a vector whose starting and ending index */
/*  are also input via arguments.  The routines include: */
/*  EOS   - equation of state; computes gas pressure and derivative wrt e */
/*  TEMP  - computes material temperature and derivative wrt e */
/*  PLANCK- computes the frequency integrated Planck function and */
/*              derivative wrt temp */
/*  ABSORP - computes absorption coefficient and derivative wrt temp */
/*  SCATT  - computes scattering coefficient */
/* ----------------------------------------------------------------------- */

/* Subroutine */ int eos_(real *e, real *d__, SFTYPE *gam, integer *istrt, 
	integer *iend, real *p, real *dpde)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

    /* Parameter adjustments */
    --dpde;
    --p;
    --d__;
    --e;

    /* Function Body */
    i__1 = *iend;
    for (i__ = *istrt; i__ <= i__1; ++i__) {
	p[i__] = (*gam - 1.f) * e[i__];
	dpde[i__] = *gam - 1.f;
/* L10: */
    }
    return 0;
} /* eos_ */


/* Subroutine */ int temp_(real *e, real *d__, SFTYPE *gam, integer *istrt, 
	integer *iend, real *t, real *dtde)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

    /* Parameter adjustments */
    --dtde;
    --t;
    --d__;
    --e;

    /* Function Body */
    i__1 = *iend;
    for (i__ = *istrt; i__ <= i__1; ++i__) {
	t[i__] = (*gam - 1.f) * e[i__] / (d__[i__] )*mmw*mhydrnew/kcontnew;
	dtde[i__] = (*gam - 1.f) / (d__[i__] )*mmw*mhydrnew/kcontnew;
/* L10: */
    }
    return 0;
} /* temp_ */


/* Subroutine */ int planck_(real *t, integer *istrt, integer *iend, real *b, 
	real *db)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer i__;

    /* Parameter adjustments */
    --db;
    --b;
    --t;

    /* Function Body */
    i__1 = *iend;
    for (i__ = *istrt; i__ <= i__1; ++i__) {
/* Computing 3rd power */
	r__1 = t[i__];
	db[i__] = r__1 * (r__1 * r__1) * sigmacontnew/3.1415926;
	b[i__] = db[i__] * t[i__];
	db[i__] *= 4.f;
/* L10: */
    }
    return 0;
} /* planck_ */


/* Subroutine */ int absorp_(real *t, real *d__, integer *istrt, integer *
	iend, real *xe, real *dxe)
{
    /* System generated locals */
    integer i__1;
    double dentemp;
    double ttemp;
    double opatemp;


    /* Local variables */
    static integer i__;

    /* Parameter adjustments */
    --dxe;
    --xe;
    --d__;
    --t;

    /* Function Body */
    i__1 = *iend;
    for (i__ = *istrt; i__ <= i__1; ++i__) {
	dentemp=d__[i__]*rhounit;
        ttemp=t[i__];
        opatemp=opac(dentemp, ttemp);
	xe[i__] = opatemp/lunit/lunit*munit*d__[i__];
	dxe[i__] = 0.f/lunit/lunit*munit*d__[i__];
/* L10: */
    }
    return 0;
} /* absorp_ */


/* Subroutine */ int scatt_(real *t, real *d__, integer *istrt, integer *iend,
	 real *xe)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

    /* Parameter adjustments */
    --xe;
    --d__;
    --t;

    /* Function Body */
    i__1 = *iend;
    for (i__ = *istrt; i__ <= i__1; ++i__) {
	xe[i__] = 0.f/lunit/lunit*munit*d__[i__];
/* L10: */
    }
    return 0;
} /* scatt_ */

