/* iccgaf.f -- translated by f2c (version 20061008).
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
/* Common Block Declarations */

struct {
    integer kmic, lmic, kmic2x, numelts;
} ciccg_;

#define ciccg_1 ciccg_

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* RAF */
/* RAF Eliminate the use of pointers by aliasing some arrays through */
/* RAF the subroutine call -- yucky, but portable! */
/* RAF */
/* RAF      subroutine iccgaf (km,lm,eps,ks,maxit, */
/* RAF     .                   a0,a1,b0,b1,bm1,x,y,work) */
/* Subroutine */ int iccgaf_(integer *km, integer *lm, real *eps, integer *ks,
	 integer *maxit, real *temp, real *a0save, real *d__, real *c0, real *
	w, real *solr, real *c1, real *cm1, real *r__, real *p, real *a0, 
	real *a1, real *b0, real *b1, real *bm1, real *x, real *y, real *work)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1;


    /* Local variables */
    extern /* Subroutine */ int tsdecomp_(real *, real *, real *, integer *), 
	    altevens_(real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, integer *, integer *);
    static integer inextlev;
    static real b;
    static integer i__;
    static real aa;
    static integer kp, kq, nbb, iter;
    static real dotp, dotr;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real rerr, xerr;
    extern doublereal snrm2_(integer *, real *, integer *);
    static integer noddb;
    extern doublereal sasum_(integer *, real *, integer *);
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), saxpy_(integer *, real *, real *, integer *, real *, 
	    integer *), tsdby1_(real *, real *, real *, integer *);
    static integer noddcb, nevenb, levlen;
    extern /* Subroutine */ int matmul_(real *, real *, real *, real *, real *
	    , real *, real *, integer *, integer *);
    static real yabsum;
    static integer ibegine;
    extern /* Subroutine */ int genbees_(real *, real *, real *, real *, real 
	    *, real *, real *, real *, real *, real *, integer *), gencees_(
	    real *, real *, real *, real *, real *, real *, real *, integer *)
	    ;
    static integer ibegino, nevencb;
    static real dotprev;
    extern /* Subroutine */ int tssolve_(real *, real *, real *, real *, real 
	    *, real *, real *, real *, real *, integer *);

/* ------------------------------------------------------------- */
/* -- incomplete cholesky - conjugate gradient method */
/* -- written by alex friedman, llnl, 415-422-0827 */
/* -- using algorithms of david kershaw, donald wolitzer. */
/* -- input arrays are all twice as long as needed to specify */
/* -- the problem; the second half is workspace.  that is, the */
/* -- input arrays should all have dimension at least */
/* -- 2*kmic*lmic. these arrays are: a0,a1,b0,b1,bm1,x,y. */
/* -- work is an additional array of length 4*kmic*lmic */
/* -- needed for working space. */
/* -- The matrix to be solved is symmetric, and has the structure: */
/* -- */
/* --   a0 */
/* --     1 */
/* -- */
/* --   a1   a0 */
/* --     1    2 */
/* -- */
/* --        a1   a0 */
/* --          2    3 */
/* -- */
/* --       : */
/* --       : */
/* --       : */
/* --       : */
/* -- */
/* --                           a1    a0 */
/* --                             km-2  km-1 */
/* -- */
/* --                                 a1    a0 */
/* --                                   km-1  km */
/* -- */
/* --   b0   bm1                                  a0 */
/* --     1     1                                  km+1 */
/* -- */
/* --   b1    b0   bm1                            a1    a0 */
/* --     1     2     2                             km+1  km+2 */
/* -- */
/* --         b1    b0   bm1                            a1    a0 */
/* --           2     3     3                             km+2  km+3 */
/* -- */
/* --                :                                        : */
/* --                :                                        : */
/* --                :                                        : */
/* -- */
/* -- Note that elements a1(sub)km etc. are zero due to the block structure. */
/* -- the correspondence with lasnex convention is: */
/* --   a0 <--> a */
/* --   a1 <--> b */
/* --   b0 <--> g */
/* --   b1 <--> d */
/* --   bm1<--> e */
/* -- note that the first physical element of bm1 is bm1(1), in */
/* -- contrast with the lasnex convention where it is bm1(2). */
/* -- scalar input variables are: */
/* -- km = kmic - the "short" dimension of the physical system */
/* --             (k is the rapidly varying index of the mesh array). */
/* -- lm = lmic - the "long" dimension of the physical system */
/* --             (l is the slowly varying index of the mesh array). */
/* --      eps  - convergence criterion, l2 normalized to y vector. */
/* --             on return, the actual accuracy achieved. */
/* --      ks   - last complete level of cyclic reduction desired */
/* --             in loops (one more level is done outside loops). */
/* --             the minimum possible value of ks is 0. */
/* --             the maximum possible value of ks is kp-1, where */
/* --             kp is the highest power of 2 with 2**kp .le. lmic. */
/* --             a reasonable choice is ks=4 for "most" problems. */
/* --             ks.le.14 at present due to dimension (15) */
/* --             in routine tssolve, allows lmic of 32k (big enough). */
/* --     maxit - the maximum number of c.g. iterations desired. */
/* --             on return, the number of iterations used. */
/* -- this package uses routines: */
/* --   TSDECOMP    GENCEES    ALTEVENS    GENBEES    MATMUL */
/* --   TSSOLVE     FORWARD    BACKWARD    ALTWS      MOVEXE */
/* --   FORMT1      TSDBY1     FORBY1      BACKBY1 */
/* --------------------------------------------------------------- */
/* RAF      pointer (pa0s,a0save),(pd,d),(pc0,c0),(pw,w),(pc1,c1), */
/* RAF     .        (pcm1,cm1),(psolr,solr),(pp,p),(pr,r),(ptemp,temp) */
/* ---- set common variables */
    /* Parameter adjustments */
    --work;
    --y;
    --x;
    --bm1;
    --b1;
    --b0;
    --a1;
    --a0;
    --p;
    --r__;
    --cm1;
    --c1;
    --solr;
    --w;
    --c0;
    --d__;
    --a0save;
    --temp;

    /* Function Body */
    ciccg_1.kmic = *km;
    ciccg_1.lmic = *lm;
    ciccg_1.kmic2x = ciccg_1.kmic << 1;
    ciccg_1.numelts = ciccg_1.kmic * ciccg_1.lmic;
/* ---- set pointers (note that many arrays share storage). */
/* RAF      ptemp = loc (x(numelts+1)) */
/* RAF      pa0s  = loc (y(numelts+1)) */
/* RAF      pd    = loc (a0) */
/* RAF      pc0   = loc (work) */
/* RAF      pw    = loc (work) */
/* RAF      psolr = loc (work) */
/* RAF      pc1   = loc (work(numelts+1)) */
/* RAF      pcm1  = loc (work(numelts+1)) */
/* RAF      pr    = loc (work(2*numelts+1)) */
/* RAF      pp    = loc (work(3*numelts+1)) */
/* ---- zero out work spaces */
    i__1 = ciccg_1.numelts + ciccg_1.numelts;
    for (i__ = ciccg_1.numelts + 1; i__ <= i__1; ++i__) {
	a0[i__] = 0.f;
	a1[i__] = 0.f;
	b0[i__] = 0.f;
	b1[i__] = 0.f;
	bm1[i__] = 0.f;
	x[i__] = 0.f;
/* L10: */
    }
    i__1 = ciccg_1.numelts << 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[i__] = 0.f;
/* L20: */
    }
/* ---- copy a0 into a0save to preserve it for matrix multiplies later */
    scopy_(&ciccg_1.numelts, &a0[1], &c__1, &a0save[1], &c__1);
/* ---- compute kp (Kershaw's p, the maximum possible level of */
/* ---- cyclic reduction).  this is used as a check on input ks. */
    for (i__ = 1; i__ <= 1000; ++i__) {
	if (pow(c__2, i__) > ciccg_1.lmic) {
	    goto L210;
	}
/* L200: */
    }
L210:
    kp = i__ - 1;
/* Computing MIN */
    i__1 = *ks, i__2 = kp - 1;
    *ks = min(i__1,i__2);
/* ---- set parameters for first pass through decomposition loop */
    nevenb = ciccg_1.lmic;
    inextlev = 1;
/* ----------------------------------------------------------------- */
/* ---- begin incomplete cholesky decomposition */
/* ----------------------------------------------------------------- */
/* ---- note that the c's are overwritten, level by level, */
/* ---- and so the pointers to c0, c1, cm1 do not move. */
    i__1 = *ks;
    for (kq = 0; kq <= i__1; ++kq) {
/* ---- set parameters for this level */
	levlen = ciccg_1.kmic * nevenb;
	noddb = (nevenb + 1) / 2;
	nevenb /= 2;
	ibegino = inextlev;
	ibegine = ibegino + ciccg_1.kmic;
	inextlev = ibegino + levlen;
/* ---- decompose upper left corner (calculate d's) */
	tsdecomp_(&a0[ibegino], &a1[ibegino], &d__[ibegino], &noddb);
/* ---- generate odd c's */
	noddcb = nevenb;
	gencees_(&b0[ibegino], &b1[ibegino], &bm1[ibegino], &a1[ibegino], &
		d__[ibegino], &c0[1], &cm1[1], &noddcb);
/* ---- generate even c's */
	nevencb = noddb - 1;
	gencees_(&b0[ibegine], &bm1[ibegine], &b1[ibegine], &a1[ibegino + 
		ciccg_1.kmic2x], &d__[ibegino + ciccg_1.kmic2x], &c0[
		ciccg_1.kmic + 1], &c1[ciccg_1.kmic + 1], &nevencb);
/* ---- modify even diagonal arrays, lower right corner (calculate */
/* ---- atilde's).  note that c1odd = b1odd, cm1even = bm1even. */
	altevens_(&a0[ibegine], &a1[ibegine], &c0[1], &b1[ibegino], &cm1[1], &
		d__[ibegino], &c0[ciccg_1.kmic + 1], &c1[ciccg_1.kmic + 1], &
		bm1[ibegine], &d__[ibegino + ciccg_1.kmic2x], &a0[inextlev], &
		a1[inextlev], &noddcb, &nevencb);
/* ---- calculate off-diagonal elements, lower right corner (btilde's) */
	nbb = nevenb - 1;
	genbees_(&c0[ciccg_1.kmic2x + 1], &b1[ibegino + ciccg_1.kmic2x], &cm1[
		ciccg_1.kmic2x + 1], &c0[ciccg_1.kmic + 1], &c1[ciccg_1.kmic 
		+ 1], &bm1[ibegine], &d__[ibegino + ciccg_1.kmic2x], &b0[
		inextlev], &b1[inextlev], &bm1[inextlev], &nbb);
/* L1000: */
    }
/* ---- do final level of tridiagonal sym. decomposition */
    tsdby1_(&a0[inextlev], &a1[inextlev], &d__[inextlev], &nevenb);
/* ---------------------------------------------------------------- */
/* ---- end decomposition */
/* ---- begin generalized conjugate gradient */
/* ---------------------------------------------------------------- */
/* ---- form product A*x in work space w */
    matmul_(&a0save[1], &a1[1], &b0[1], &b1[1], &bm1[1], &x[1], &w[1], &
	    ciccg_1.kmic, &ciccg_1.lmic);
/* ---- form residual r = y - A*x (this is r(nought)) */
    i__1 = ciccg_1.numelts;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L1200: */
	r__[i__] = y[i__] - w[i__];
    }
/* ---- compute y norm for all relative error tests */
    yabsum = sasum_(&ciccg_1.numelts, &y[1], &c__1)+SMALL;
/* ---- compute solr = (LLt)-1 on r */
    tssolve_(&solr[1], &r__[1], &d__[1], &a1[1], &b0[1], &b1[1], &bm1[1], &
	    solr[1], &temp[1], ks);
/* ---- ... and set p(nought) to this. */
/*     call zmovewrd (p,solr,numelts) */
    scopy_(&ciccg_1.numelts, &solr[1], &c__1, &p[1], &c__1);
/* ---- set up previous dotr product for iteration */
    dotprev = sdot_(&ciccg_1.numelts, &r__[1], &c__1, &solr[1], &c__1);
/* --------------------------------------------------------------- */
/* ---- begin conjugate gradient iteration loop */
    i__1 = *maxit;
    for (iter = 1; iter <= i__1; ++iter) {
/* ------- compute m*p in temp (this is wolitzer's "evalp").  note - */
/* ------- tssolve uses temp array, but we'll be done with it by then. */
	matmul_(&a0save[1], &a1[1], &b0[1], &b1[1], &bm1[1], &p[1], &temp[1], 
		&ciccg_1.kmic, &ciccg_1.lmic);
/* ------- (p,mp) and exit if done. */
	dotp = sdot_(&ciccg_1.numelts, &p[1], &c__1, &temp[1], &c__1);
	if (dotp == 0.f) {
	    goto L3000;
	}
/* ------- compute aa, the ratio of old dotr product and (p,mp). */
	aa = dotprev / dotp;
/* ------- x = x + aa * p */
	saxpy_(&ciccg_1.numelts, &aa, &p[1], &c__1, &x[1], &c__1);
/* ------- r = r - aa * m*p */
	aa = -aa;
	saxpy_(&ciccg_1.numelts, &aa, &temp[1], &c__1, &r__[1], &c__1);
/* ------- compute error measure */
	rerr = (r__1 = sasum_(&ciccg_1.numelts, &r__[1], &c__1), dabs(r__1)) /
		 yabsum;
	xerr = dabs(aa) * snrm2_(&ciccg_1.numelts, &p[1], &c__1) / snrm2_(&
		ciccg_1.numelts, &x[1], &c__1);
/* ------- test if done */
	if (xerr < *eps && rerr < *eps) {
	    goto L3000;
	}
/* ------- compute new dotr product (r,solr) using solr = (LLt)-1 r */
	tssolve_(&solr[1], &r__[1], &d__[1], &a1[1], &b0[1], &b1[1], &bm1[1], 
		&solr[1], &temp[1], ks);
	dotr = sdot_(&ciccg_1.numelts, &r__[1], &c__1, &solr[1], &c__1);
/* ------- b is the ratio of old to new dotr prods.  reset dotrprev. */
	b = dotr / dotprev;
	dotprev = dotr;
/* ------- p = (llt)-1 r + b*p */
	i__2 = ciccg_1.numelts;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L1600: */
	    p[i__] = solr[i__] + b * p[i__];
	}
/* L2000: */
    }
/* ---- if this point reached, no convergence after maxit passes. */
/* ---------------------------------------------------------------- */
/* ---- end generalized conjugate gradient */
/* ---------------------------------------------------------------- */
L3000:
/* ---- set parameters to the values they must have on return */
    *eps = rerr;
    *maxit = iter;
/* ---- restore a0 array for user convenience */
/*     call zmovewrd (a0,a0save,numelts) */
    scopy_(&ciccg_1.numelts, &a0save[1], &c__1, &a0[1], &c__1);
    return 0;
} /* iccgaf_ */

/* Subroutine */ int tsdecomp_(real *a0, real *a1, real *d__, integer *
	nblocks)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static integer i__, j, lastbl;

/* ---- note that the d's overwrite the odd a0's at each level */
/* ---- negative pivots are eliminated using absolute value. */
/* ---- no check is made for small denominators. */
/* ---- pointer to beginning of last block to be processed */
    /* Parameter adjustments */
    --d__;
    --a1;
    --a0;

    /* Function Body */
    lastbl = (*nblocks - 1) * ciccg_1.kmic2x + 1;
/* ---- do first row of every other block (compute d's) */
    i__1 = lastbl;
    i__2 = ciccg_1.kmic2x;
    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* L100: */
	d__[j] = (r__1 = 1.f / a0[j], dabs(r__1));
    }
/* ---- do ith row of every other block (compute d's) */
    i__2 = ciccg_1.kmic;
    for (i__ = 2; i__ <= i__2; ++i__) {
/* dir$ ivdep */
	i__1 = lastbl - 1;
	i__3 = ciccg_1.kmic2x;
	for (j = 0; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
/* L200: */
	    d__[i__ + j] = (r__1 = 1.f / (a0[i__ + j] - a1[i__ + j - 1] * a1[
		    i__ + j - 1] * d__[i__ + j - 1]), dabs(r__1));
	}
    }
    return 0;
} /* tsdecomp_ */

/* Subroutine */ int gencees_(real *b0, real *b1, real *bm1, real *a1, real *
	d__, real *c0, real *cm1, integer *nblocks)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, lastbl;

/* ---- test <<before>> entering loop */
    /* Parameter adjustments */
    --cm1;
    --c0;
    --d__;
    --a1;
    --bm1;
    --b1;
    --b0;

    /* Function Body */
    if (*nblocks <= 0) {
	return 0;
    }
/* ---- pointer to beginning of last block to be processed */
    lastbl = (*nblocks - 1) * ciccg_1.kmic2x + 1;
/* ---- do first row of every other block (compute c0's, c1's) */
    i__1 = lastbl;
    i__2 = ciccg_1.kmic2x;
    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
	c0[j] = b0[j];
/* L100: */
	cm1[j] = bm1[j] - a1[j] * d__[j] * c0[j];
    }
/* ---- do ith row of every other block (compute c0's, c1's) */
    i__2 = ciccg_1.kmic;
    for (i__ = 2; i__ <= i__2; ++i__) {
	i__1 = lastbl - 1;
	i__3 = ciccg_1.kmic2x;
	for (j = 0; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
	    c0[i__ + j] = b0[i__ + j] - a1[i__ + j - 1] * d__[i__ + j - 1] * 
		    b1[i__ + j - 1];
/* L200: */
	    cm1[i__ + j] = bm1[i__ + j] - a1[i__ + j] * d__[i__ + j] * c0[i__ 
		    + j];
	}
    }
    return 0;
} /* gencees_ */

/* Subroutine */ int altevens_(real *a0k, real *a1k, real *c0km1, real *c1km1,
	 real *cm1km1, real *dkm1, real *c0k, real *c1k, real *cm1k, real *
	dkp1, real *a0tild, real *a1tild, integer *nbkm1, integer *nbkkp1)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer lastaltb, lastnewb, i__, j, jj;

/* ---- test <<before>> entering loop */
    /* Parameter adjustments */
    --a1tild;
    --a0tild;
    --dkp1;
    --cm1k;
    --c1k;
    --c0k;
    --dkm1;
    --cm1km1;
    --c1km1;
    --c0km1;
    --a1k;
    --a0k;

    /* Function Body */
    if (*nbkm1 <= 0) {
	return 0;
    }
/* ---- pointer to beginning of last new block */
    lastnewb = (*nbkm1 - 1) * ciccg_1.kmic + 1;
/* ---------------------------------------------------------------- */
/* ---- first set atilde to:  a - (k-1 term). */
/* ---------------------------------------------------------------- */
/* ---- do first row of every other block */
    j = 1 - ciccg_1.kmic2x;
    i__1 = lastnewb;
    i__2 = ciccg_1.kmic;
    for (jj = 1; i__2 < 0 ? jj >= i__1 : jj <= i__1; jj += i__2) {
	j += ciccg_1.kmic2x;
	a0tild[jj] = a0k[j] - (c0km1[j] * c0km1[j] * dkm1[j] + cm1km1[j] * 
		cm1km1[j] * dkm1[j + 1]);
	a1tild[jj] = a1k[j] - (c0km1[j] * c1km1[j] * dkm1[j] + cm1km1[j] * 
		c0km1[j + 1] * dkm1[j + 1]);
/* L100: */
    }
/* ---- do ith row of every other block */
    i__2 = ciccg_1.kmic;
    for (i__ = 2; i__ <= i__2; ++i__) {
	j = -ciccg_1.kmic2x;
	i__1 = lastnewb - 1;
	i__3 = ciccg_1.kmic;
	for (jj = 0; i__3 < 0 ? jj >= i__1 : jj <= i__1; jj += i__3) {
	    j += ciccg_1.kmic2x;
	    a0tild[i__ + jj] = a0k[i__ + j] - (c0km1[i__ + j] * c0km1[i__ + j]
		     * dkm1[i__ + j] + cm1km1[i__ + j] * cm1km1[i__ + j] * 
		    dkm1[i__ + j + 1] + c1km1[i__ + j - 1] * c1km1[i__ + j - 
		    1] * dkm1[i__ + j - 1]);
	    a1tild[i__ + jj] = a1k[i__ + j] - (c0km1[i__ + j] * c1km1[i__ + j]
		     * dkm1[i__ + j] + cm1km1[i__ + j] * c0km1[i__ + j + 1] * 
		    dkm1[i__ + j + 1]);
/* L200: */
	}
    }
/* ------------------------------------------------------------------ */
/* ---- subtract off k/k+1 term to form true atilde. */
/* ---- this is done in a separate loop to avoid "overreach" */
/* ---- problems.  since the c pointers do not move, we can't */
/* ---- count on the last block of c's being zero. */
/* ------------------------------------------------------------------ */
/* ---- test <<before>> entering loop */
    if (*nbkkp1 <= 0) {
	return 0;
    }
/* ---- pointer to beginning of last new block to alter */
    lastaltb = (*nbkkp1 - 1) * ciccg_1.kmic + 1;
/* ---- do first row of every other block */
    j = 1 - ciccg_1.kmic2x;
    i__3 = lastaltb;
    i__1 = ciccg_1.kmic;
    for (jj = 1; i__1 < 0 ? jj >= i__3 : jj <= i__3; jj += i__1) {
	j += ciccg_1.kmic2x;
	a0tild[jj] -= c0k[j] * c0k[j] * dkp1[j] + c1k[j] * c1k[j] * dkp1[j + 
		1];
	a1tild[jj] -= c0k[j] * cm1k[j] * dkp1[j] + c1k[j] * c0k[j + 1] * dkp1[
		j + 1];
/* L300: */
    }
/* ---- do ith row of every other block */
    i__1 = ciccg_1.kmic;
    for (i__ = 2; i__ <= i__1; ++i__) {
	j = -ciccg_1.kmic2x;
	i__3 = lastaltb - 1;
	i__2 = ciccg_1.kmic;
	for (jj = 0; i__2 < 0 ? jj >= i__3 : jj <= i__3; jj += i__2) {
	    j += ciccg_1.kmic2x;
	    a0tild[i__ + jj] -= c0k[i__ + j] * c0k[i__ + j] * dkp1[i__ + j] + 
		    c1k[i__ + j] * c1k[i__ + j] * dkp1[i__ + j + 1] + cm1k[
		    i__ + j - 1] * cm1k[i__ + j - 1] * dkp1[i__ + j - 1];
	    a1tild[i__ + jj] -= c0k[i__ + j] * cm1k[i__ + j] * dkp1[i__ + j] 
		    + c1k[i__ + j] * c0k[i__ + j + 1] * dkp1[i__ + j + 1];
/* L400: */
	}
    }
    return 0;
} /* altevens_ */

/* Subroutine */ int genbees_(real *c0kp1, real *c1kp1, real *cm1kp1, real *
	c0k, real *c1k, real *cm1k, real *dkp1, real *b0tild, real *b1tild, 
	real *bm1tild, integer *nblocks)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer lastnewb, i__, j, jj;

/* ---- test <<before>> entering loop */
    /* Parameter adjustments */
    --bm1tild;
    --b1tild;
    --b0tild;
    --dkp1;
    --cm1k;
    --c1k;
    --c0k;
    --cm1kp1;
    --c1kp1;
    --c0kp1;

    /* Function Body */
    if (*nblocks <= 0) {
	return 0;
    }
/* ---- pointer to beginning of last new block */
    lastnewb = (*nblocks - 1) * ciccg_1.kmic + 1;
/* ---- do first row of every other block */
    j = 1 - ciccg_1.kmic2x;
    i__1 = lastnewb;
    i__2 = ciccg_1.kmic;
    for (jj = 1; i__2 < 0 ? jj >= i__1 : jj <= i__1; jj += i__2) {
	j += ciccg_1.kmic2x;
	b0tild[jj] = -c0kp1[j] * c0k[j] * dkp1[j] - cm1kp1[j] * c1k[j] * dkp1[
		j + 1];
	b1tild[jj] = -c1kp1[j] * c0k[j] * dkp1[j] - c0kp1[j + 1] * c1k[j] * 
		dkp1[j + 1];
	bm1tild[jj] = -c0kp1[j] * cm1k[j] * dkp1[j] - cm1kp1[j] * c0k[j + 1] *
		 dkp1[j + 1];
/* L100: */
    }
/* ---- do ith row of every block */
    i__2 = ciccg_1.kmic;
    for (i__ = 2; i__ <= i__2; ++i__) {
	j = -ciccg_1.kmic2x;
	i__1 = lastnewb - 1;
	i__3 = ciccg_1.kmic;
	for (jj = 0; i__3 < 0 ? jj >= i__1 : jj <= i__1; jj += i__3) {
	    j += ciccg_1.kmic2x;
	    b0tild[i__ + jj] = -c0kp1[i__ + j] * c0k[i__ + j] * dkp1[i__ + j] 
		    - cm1kp1[i__ + j] * c1k[i__ + j] * dkp1[i__ + j + 1] - 
		    c1kp1[i__ + j - 1] * cm1k[i__ + j - 1] * dkp1[i__ + j - 1]
		    ;
	    b1tild[i__ + jj] = -c1kp1[i__ + j] * c0k[i__ + j] * dkp1[i__ + j] 
		    - c0kp1[i__ + j + 1] * c1k[i__ + j] * dkp1[i__ + j + 1];
	    bm1tild[i__ + jj] = -c0kp1[i__ + j] * cm1k[i__ + j] * dkp1[i__ + 
		    j] - cm1kp1[i__ + j] * c0k[i__ + j + 1] * dkp1[i__ + j + 
		    1];
/* L200: */
	}
    }
    return 0;
} /* genbees_ */

/* Subroutine */ int matmul_(real *a0, real *a1, real *b0, real *b1, real *
	bm1, real *x, real *w, integer *kmic, integer *lmic)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ---- computes w = A*x. */
    /* Parameter adjustments */
    --w;
    --x;
    --bm1;
    --b1;
    --b0;
    --a1;
    --a0;

    /* Function Body */
    w[1] = a0[1] * x[1] + a1[1] * x[2];
    i__1 = *kmic * *lmic;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L100: */
	w[i__] = a1[i__ - 1] * x[i__ - 1] + a0[i__] * x[i__] + a1[i__] * x[
		i__ + 1];
    }
    w[1] = w[1] + b0[1] * x[*kmic + 1] + b1[1] * x[*kmic + 2];
    i__1 = *kmic * (*lmic - 1);
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L200: */
	w[i__] = w[i__] + bm1[i__ - 1] * x[*kmic + i__ - 1] + b0[i__] * x[*
		kmic + i__] + b1[i__] * x[*kmic + i__ + 1];
    }
    w[*kmic + 1] = w[*kmic + 1] + b0[1] * x[1] + bm1[1] * x[2];
    i__1 = *kmic * (*lmic - 1);
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L300: */
	w[*kmic + i__] = w[*kmic + i__] + b1[i__ - 1] * x[i__ - 1] + b0[i__] *
		 x[i__] + bm1[i__] * x[i__ + 1];
    }
    return 0;
} /* matmul_ */

/* Subroutine */ int tssolve_(real *x, real *y, real *d__, real *a1, real *b0,
	 real *b1, real *bm1, real *w, real *temp, integer *ks)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int backward_(real *, real *, real *, real *, 
	    integer *);
    static integer inextlev, i__, j, kq, ibe[15], neb[15], ibo[15], nob[15], 
	    inl[15], noddb;
    extern /* Subroutine */ int altws_(real *, real *, real *, real *, real *,
	     real *, real *, real *, real *, real *, integer *, integer *), 
	    scopy_(integer *, real *, integer *, real *, integer *), forby1_(
	    real *, real *, real *, real *, integer *), formt1_(real *, real *
	    , real *, real *, real *, real *, real *, real *, real *, integer 
	    *);
    static integer nevenb, levlen;
    extern /* Subroutine */ int movexe_(real *, real *, integer *), backby1_(
	    real *, real *, real *, real *, integer *);
    static integer ibegine, ibegino;
    extern /* Subroutine */ int forward_(real *, real *, real *, real *, 
	    integer *);

/* ---- solves LDLt*x=y by (first) L*w=y (then) DLt*x=w */
/* ---- note that w and x use same storage (solr).  they are */
/* ---- kept distinct here for purposes of clarity. */
/* ---------------------------------------------------------------- */
/* ---- begin forward sweep  L*w=y */
/* ---------------------------------------------------------------- */
/* ---- set w at kq=0 level to y */
/*     call zmovewrd (w,y,numelts) */
    /* Parameter adjustments */
    --temp;
    --w;
    --bm1;
    --b1;
    --b0;
    --a1;
    --d__;
    --y;
    --x;

    /* Function Body */
    scopy_(&ciccg_1.numelts, &y[1], &c__1, &w[1], &c__1);
/* ---- set parameters for first pass through forward sweep loop */
    nevenb = ciccg_1.lmic;
    inextlev = 1;
/* ---- begin loop */
    i__1 = *ks;
    for (kq = 0; kq <= i__1; ++kq) {
/* ---- set parameters for this level */
	levlen = ciccg_1.kmic * nevenb;
	noddb = (nevenb + 1) / 2;
	nevenb /= 2;
	ibegino = inextlev;
	ibegine = ibegino + ciccg_1.kmic;
	inextlev = ibegino + levlen;
/* ---- save parameters for this level for future use */
/* ---- in backward sweep (they are tricky to regenerate */
/* ---- when moving backwards). */
/* ---- kershaw, instead, generates the neven's using a circular */
/* ---- shift right and mask, so the process is reversible. */
	nob[kq] = noddb;
	neb[kq] = nevenb;
	ibo[kq] = ibegino;
	ibe[kq] = ibegine;
	inl[kq] = inextlev;
/* ---- perform forward solve on odd blocks, L*w=w */
	forward_(&d__[ibegino], &a1[ibegino], &w[ibegino], &w[ibegino], &
		noddb);
/* ---- generate (L)-t (D)-1 * w terms for odd blocks, */
/* ---- i.e. DLt*temp=w */
	backward_(&d__[ibegino], &a1[ibegino], &w[ibegino], &temp[1], &noddb);
/* ---- create next level of w's */
	i__2 = noddb - 1;
	altws_(&w[ibegine], &b0[ibegino], &b1[ibegino], &bm1[ibegino], &temp[
		1], &b0[ibegine], &b1[ibegine], &bm1[ibegine], &temp[
		ciccg_1.kmic2x + 1], &w[inextlev], &nevenb, &i__2);
/* L1000: */
    }
/* ---- do kq=ks+1 forward solve */
    forby1_(&d__[inextlev], &a1[inextlev], &w[inextlev], &w[inextlev], &
	    nevenb);
/* ------------------------------------------------------------- */
/* ---- end forward sweep */
/* ---- begin backward sweep */
/* ------------------------------------------------------------- */
/* ---- do kq=ks+1 backward solve */
    backby1_(&d__[inextlev], &a1[inextlev], &w[inextlev], &x[inextlev], &
	    nevenb);
/* ---- begin loop */
    for (kq = *ks; kq >= 0; --kq) {
/* ---- set parameters for this level */
	noddb = nob[kq];
	nevenb = neb[kq];
	ibegino = ibo[kq];
	ibegine = ibe[kq];
	inextlev = inl[kq];
/* ---- move even blocks of x to lower level */
	movexe_(&x[ibegine], &x[inextlev], &nevenb);
/* ---- evaluate b(k)t*x(k+1) + b(k-1)*x(k-1) for odd k, store in temp */
	formt1_(&b0[ibegino], &b1[ibegino], &bm1[ibegino], &x[ibegine], &b0[
		ibegine - ciccg_1.kmic2x], &b1[ibegine - ciccg_1.kmic2x], &
		bm1[ibegine - ciccg_1.kmic2x], &x[ibegine - ciccg_1.kmic2x], &
		temp[1], &noddb);
/* ---- perform forward solve on this, keep in temp */
	forward_(&d__[ibegino], &a1[ibegino], &temp[1], &temp[1], &noddb);
/* ---- ... and subtract from w, still in temp */
	i__1 = ciccg_1.kmic - 1;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    i__2 = (noddb - 1) * ciccg_1.kmic2x;
	    i__3 = ciccg_1.kmic2x;
	    for (j = 0; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3) {
/* L1200: */
		temp[i__ + 1 + j] = w[ibegino + i__ + j] - temp[i__ + 1 + j];
	    }
	}
/* ---- finally, perform backward solve into x */
	backward_(&d__[ibegino], &a1[ibegino], &temp[1], &x[ibegino], &noddb);
/* L2000: */
    }
/* ------------------------------------------------------------ */
/* ---- end backward sweep */
/* ------------------------------------------------------------ */
    return 0;
} /* tssolve_ */

/* Subroutine */ int forward_(real *d__, real *a1, real *y, real *w, integer *
	nblocks)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, lastbl;

/* ---- solves L*w=y */
/* ---- pointer to beginning of last block to be processed */
    /* Parameter adjustments */
    --w;
    --y;
    --a1;
    --d__;

    /* Function Body */
    lastbl = (*nblocks - 1) * ciccg_1.kmic2x + 1;
/* ---- do first row of every other block */
    i__1 = lastbl;
    i__2 = ciccg_1.kmic2x;
    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* L100: */
	w[j] = d__[j] * y[j];
    }
/* ---- do ith row of every other block */
    i__2 = ciccg_1.kmic;
    for (i__ = 2; i__ <= i__2; ++i__) {
/* dir$ ivdep */
	i__1 = lastbl - 1;
	i__3 = ciccg_1.kmic2x;
	for (j = 0; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
/* L200: */
	    w[i__ + j] = d__[i__ + j] * (y[i__ + j] - a1[i__ + j - 1] * w[i__ 
		    + j - 1]);
	}
    }
    return 0;
} /* forward_ */

/* Subroutine */ int backward_(real *d__, real *a1, real *w, real *x, integer 
	*nblocks)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, lastbl;

/* ---- solves DLt*x=w */
/* ---- pointer to beginning of last block to be processed */
    /* Parameter adjustments */
    --x;
    --w;
    --a1;
    --d__;

    /* Function Body */
    lastbl = (*nblocks - 1) * ciccg_1.kmic2x + 1;
/* ---- do last row of every other block */
    i__1 = lastbl + ciccg_1.kmic - 1;
    i__2 = ciccg_1.kmic2x;
    for (j = ciccg_1.kmic; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* L100: */
	x[j] = w[j];
    }
/* ---- do ith row of every other block, progressing backwards */
    for (i__ = ciccg_1.kmic - 1; i__ >= 1; --i__) {
/* dir$ ivdep */
	i__2 = lastbl - 1;
	i__1 = ciccg_1.kmic2x;
	for (j = 0; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
/* L200: */
	    x[i__ + j] = w[i__ + j] - a1[i__ + j] * d__[i__ + j] * x[i__ + j 
		    + 1];
	}
    }
    return 0;
} /* backward_ */

/* Subroutine */ int altws_(real *wk, real *b0km1, real *b1km1, real *bm1km1, 
	real *tkm1, real *b0k, real *b1k, real *bm1k, real *tkp1, real *
	wtilde, integer *nblkslo, integer *nblkshi)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer lastnbhi, lastnblo, i__, j, jj;

/* ---- compute next level of w's, the wtildes, for the solve */
/* ---- pointers to beginning of last new block */
    /* Parameter adjustments */
    --wtilde;
    --tkp1;
    --bm1k;
    --b1k;
    --b0k;
    --tkm1;
    --bm1km1;
    --b1km1;
    --b0km1;
    --wk;

    /* Function Body */
    lastnblo = (*nblkslo - 1) * ciccg_1.kmic + 1;
    lastnbhi = (*nblkshi - 1) * ciccg_1.kmic + 1;
/* ---- do ith row of every other block */
    j = -ciccg_1.kmic2x;
    i__1 = lastnblo - 1;
    i__2 = ciccg_1.kmic;
    for (jj = 0; i__2 < 0 ? jj >= i__1 : jj <= i__1; jj += i__2) {
	j += ciccg_1.kmic2x;
	wtilde[jj + 1] = wk[j + 1] - (b0km1[j + 1] * tkm1[j + 1] + bm1km1[j + 
		1] * tkm1[j + 2]);
/* L10: */
    }
    j = -ciccg_1.kmic2x;
    i__2 = lastnbhi - 1;
    i__1 = ciccg_1.kmic;
    for (jj = 0; i__1 < 0 ? jj >= i__2 : jj <= i__2; jj += i__1) {
	j += ciccg_1.kmic2x;
	wtilde[jj + 1] -= b0k[j + 1] * tkp1[j + 1] + b1k[j + 1] * tkp1[j + 2];
/* L20: */
    }
    i__1 = ciccg_1.kmic;
    for (i__ = 2; i__ <= i__1; ++i__) {
	j = -ciccg_1.kmic2x;
	i__2 = lastnblo - 1;
	i__3 = ciccg_1.kmic;
	for (jj = 0; i__3 < 0 ? jj >= i__2 : jj <= i__2; jj += i__3) {
	    j += ciccg_1.kmic2x;
	    wtilde[i__ + jj] = wk[i__ + j] - (b1km1[i__ + j - 1] * tkm1[i__ + 
		    j - 1] + b0km1[i__ + j] * tkm1[i__ + j] + bm1km1[i__ + j] 
		    * tkm1[i__ + j + 1]);
/* L100: */
	}
	j = -ciccg_1.kmic2x;
	i__3 = lastnbhi - 1;
	i__2 = ciccg_1.kmic;
	for (jj = 0; i__2 < 0 ? jj >= i__3 : jj <= i__3; jj += i__2) {
	    j += ciccg_1.kmic2x;
	    wtilde[i__ + jj] -= bm1k[i__ + j - 1] * tkp1[i__ + j - 1] + b0k[
		    i__ + j] * tkp1[i__ + j] + b1k[i__ + j] * tkp1[i__ + j + 
		    1];
/* L200: */
	}
    }
    return 0;
} /* altws_ */

/* Subroutine */ int movexe_(real *xk, real *xtilde, integer *nblocks)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer lastnewb, i__, j, jj;

/* ---- moves block k/2 of xtilde (=x at level kq+1) into block k */
/* ---- of xk (=x at level kq), for k even */
/* ---- pointer to beginning of last new block */
    /* Parameter adjustments */
    --xtilde;
    --xk;

    /* Function Body */
    lastnewb = (*nblocks - 1) * ciccg_1.kmic + 1;
/* ---- do ith row of every other block */
    i__1 = ciccg_1.kmic;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = -ciccg_1.kmic2x;
	i__2 = lastnewb - 1;
	i__3 = ciccg_1.kmic;
	for (jj = 0; i__3 < 0 ? jj >= i__2 : jj <= i__2; jj += i__3) {
	    j += ciccg_1.kmic2x;
/* L100: */
	    xk[i__ + j] = xtilde[i__ + jj];
	}
    }
    return 0;
} /* movexe_ */

/* Subroutine */ int formt1_(real *b0k, real *b1k, real *bm1k, real *xkp1, 
	real *b0km1, real *b1km1, real *bm1km1, real *xkm1, real *temp, 
	integer *nblocks)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, lastbl;

/* ---- compute intermediate term in backward solve */
/* ---- pointer to beginning of last block */
    /* Parameter adjustments */
    --temp;
    --xkm1;
    --bm1km1;
    --b1km1;
    --b0km1;
    --xkp1;
    --bm1k;
    --b1k;
    --b0k;

    /* Function Body */
    lastbl = (*nblocks - 1) * ciccg_1.kmic2x + 1;
/* ---- do part valid for all blocks, ith row of every other block, */
/* ---- then do part not valid for first block */
    i__1 = lastbl - 1;
    i__2 = ciccg_1.kmic2x;
    for (j = 0; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* L10: */
	temp[j + 1] = b0k[j + 1] * xkp1[j + 1] + b1k[j + 1] * xkp1[j + 2];
    }
    i__2 = lastbl - 1;
    i__1 = ciccg_1.kmic2x;
    for (j = ciccg_1.kmic2x; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
/* L20: */
	temp[j + 1] = temp[j + 1] + b1km1[j] * xkm1[j] + b0km1[j + 1] * xkm1[
		j + 1] + bm1km1[j + 1] * xkm1[j + 2];
    }
    i__1 = ciccg_1.kmic;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = lastbl - 1;
	i__3 = ciccg_1.kmic2x;
	for (j = 0; i__3 < 0 ? j >= i__2 : j <= i__2; j += i__3) {
/* L100: */
	    temp[i__ + j] = bm1k[i__ + j - 1] * xkp1[i__ + j - 1] + b0k[i__ + 
		    j] * xkp1[i__ + j] + b1k[i__ + j] * xkp1[i__ + j + 1];
	}
	i__3 = lastbl - 1;
	i__2 = ciccg_1.kmic2x;
	for (j = ciccg_1.kmic2x; i__2 < 0 ? j >= i__3 : j <= i__3; j += i__2) 
		{
/* L200: */
	    temp[i__ + j] = temp[i__ + j] + b1km1[i__ + j - 1] * xkm1[i__ + j 
		    - 1] + b0km1[i__ + j] * xkm1[i__ + j] + bm1km1[i__ + j] * 
		    xkm1[i__ + j + 1];
	}
/* L300: */
    }
    return 0;
} /* formt1_ */

/* Subroutine */ int tsdby1_(real *a0, real *a1, real *d__, integer *nblocks)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static integer i__, j, lastbl;

/* ---- this routine is similar to tsdecomp except that it */
/* ---- does every block, not every other, as needed for */
/* ---- incomplete cyclic reduction. */
/* ---- pointer to beginning of last block to be processed */
    /* Parameter adjustments */
    --d__;
    --a1;
    --a0;

    /* Function Body */
    lastbl = (*nblocks - 1) * ciccg_1.kmic + 1;
/* ---- do first row of every block (compute d's) */
    i__1 = lastbl;
    i__2 = ciccg_1.kmic;
    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* L100: */
	d__[j] = (r__1 = 1.f / a0[j], dabs(r__1));
    }
/* ---- do ith row of every block (compute d's) */
    i__2 = ciccg_1.kmic;
    for (i__ = 2; i__ <= i__2; ++i__) {
/* dir$ ivdep */
	i__1 = lastbl - 1;
	i__3 = ciccg_1.kmic;
	for (j = 0; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
/* L200: */
	    d__[i__ + j] = (r__1 = 1.f / (a0[i__ + j] - a1[i__ + j - 1] * a1[
		    i__ + j - 1] * d__[i__ + j - 1]), dabs(r__1));
	}
    }
    return 0;
} /* tsdby1_ */

/* Subroutine */ int forby1_(real *d__, real *a1, real *y, real *w, integer *
	nblocks)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, lastbl;

/* ---- similar to forward except that it does every block */
/* ---- pointer to beginning of last block to be processed */
    /* Parameter adjustments */
    --w;
    --y;
    --a1;
    --d__;

    /* Function Body */
    lastbl = (*nblocks - 1) * ciccg_1.kmic + 1;
/* ---- do first row of every block */
    i__1 = lastbl;
    i__2 = ciccg_1.kmic;
    for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* L100: */
	w[j] = d__[j] * y[j];
    }
/* ---- do ith row of every block */
    i__2 = ciccg_1.kmic;
    for (i__ = 2; i__ <= i__2; ++i__) {
/* dir$ ivdep */
	i__1 = lastbl - 1;
	i__3 = ciccg_1.kmic;
	for (j = 0; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3) {
/* L200: */
	    w[i__ + j] = d__[i__ + j] * (y[i__ + j] - a1[i__ + j - 1] * w[i__ 
		    + j - 1]);
	}
    }
    return 0;
} /* forby1_ */

/* Subroutine */ int backby1_(real *d__, real *a1, real *w, real *x, integer *
	nblocks)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, lastbl;

/* ---- similar to backward except that it does every block */
/* ---- pointer to beginning of last block to be processed */
    /* Parameter adjustments */
    --x;
    --w;
    --a1;
    --d__;

    /* Function Body */
    lastbl = (*nblocks - 1) * ciccg_1.kmic + 1;
/* ---- do last row of every block */
    i__1 = lastbl + ciccg_1.kmic - 1;
    i__2 = ciccg_1.kmic;
    for (j = ciccg_1.kmic; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* L100: */
	x[j] = w[j];
    }
/* ---- do ith row of every block, progressing backwards */
    for (i__ = ciccg_1.kmic - 1; i__ >= 1; --i__) {
/* dir$ ivdep */
	i__2 = lastbl - 1;
	i__1 = ciccg_1.kmic;
	for (j = 0; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
/* L200: */
	    x[i__ + j] = w[i__ + j] - a1[i__ + j] * d__[i__ + j] * x[i__ + j 
		    + 1];
	}
    }
    return 0;
} /* backby1_ */

