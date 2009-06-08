/* linpck.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

integer isamax_(integer *n, real *sx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;
    real r__1;

    /* Local variables */
    static integer i__, ix;
    static real smax;


/*     finds the index of element having max. absolute value. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --sx;

    /* Function Body */
    ret_val = 0;
    if (*n < 1) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    smax = dabs(sx[1]);
    ix += *incx;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if ((r__1 = sx[ix], dabs(r__1)) <= smax) {
	    goto L5;
	}
	ret_val = i__;
	smax = (r__1 = sx[ix], dabs(r__1));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        code for increment equal to 1 */

L20:
    smax = dabs(sx[1]);
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if ((r__1 = sx[i__], dabs(r__1)) <= smax) {
	    goto L30;
	}
	ret_val = i__;
	smax = (r__1 = sx[i__], dabs(r__1));
L30:
	;
    }
    return ret_val;
} /* isamax_ */


integer ismax_(integer *n, real *sx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer i__, ix;
    static real smax;


/*     finds the index of element having max. value. */
/*     jack dongarra, linpack, 3/11/78. (added to this file by j. stone) */


    /* Parameter adjustments */
    --sx;

    /* Function Body */
    ret_val = 0;
    if (*n < 1) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    smax = sx[1];
    ix += *incx;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (sx[ix] <= smax) {
	    goto L5;
	}
	ret_val = i__;
	smax = sx[ix];
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        code for increment equal to 1 */

L20:
    smax = sx[1];
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (sx[i__] <= smax) {
	    goto L30;
	}
	ret_val = i__;
	smax = sx[i__];
L30:
	;
    }
    return ret_val;
} /* ismax_ */

integer ismin_(integer *n, real *sx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer i__, ix;
    static real smin;


/*     finds the index of element having min. value. */
/*     jack dongarra, linpack, 3/11/78.  (added to this file by j. stone) */


    /* Parameter adjustments */
    --sx;

    /* Function Body */
    ret_val = 0;
    if (*n < 1) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    smin = sx[1];
    ix += *incx;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (sx[ix] >= smin) {
	    goto L5;
	}
	ret_val = i__;
	smin = sx[ix];
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        code for increment equal to 1 */

L20:
    smin = sx[1];
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (sx[i__] >= smin) {
	    goto L30;
	}
	ret_val = i__;
	smin = sx[i__];
L30:
	;
    }
    return ret_val;
} /* ismin_ */

/* Subroutine */ int scopy_(integer *n, real *sx, integer *incx, real *sy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;


/*     copies a vector, x, to a vector, y. */
/*     uses unrolled loops for increments equal to 1. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sy[iy] = sx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sy[i__] = sx[i__];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
	sy[i__] = sx[i__];
	sy[i__ + 1] = sx[i__ + 1];
	sy[i__ + 2] = sx[i__ + 2];
	sy[i__ + 3] = sx[i__ + 3];
	sy[i__ + 4] = sx[i__ + 4];
	sy[i__ + 5] = sx[i__ + 5];
	sy[i__ + 6] = sx[i__ + 6];
/* L50: */
    }
    return 0;
} /* scopy_ */

doublereal sdot_(integer *n, real *sx, integer *incx, real *sy, integer *incy)
{
    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;
    static real stemp;


/*     forms the dot product of two vectors. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    stemp = 0.f;
    ret_val = 0.f;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stemp += sx[ix] * sy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	stemp += sx[i__] * sy[i__];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
	stemp = stemp + sx[i__] * sy[i__] + sx[i__ + 1] * sy[i__ + 1] + sx[
		i__ + 2] * sy[i__ + 2] + sx[i__ + 3] * sy[i__ + 3] + sx[i__ + 
		4] * sy[i__ + 4];
/* L50: */
    }
L60:
    ret_val = stemp;
    return ret_val;
} /* sdot_ */

doublereal snrm2_(integer *n, real *sx, integer *incx)
{
    /* Initialized data */

    static real zero = 0.f;
    static real one = 1.f;
    static real cutlo = 4.441e-16f;
    static real cuthi = 1.304e19f;

    /* Format strings */
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_110[] = "";

    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, nn;
    static real sum, xmax;
    static integer next;
    static real hitest;

    /* Assigned format variables */
    static char *next_fmt;

    /* Parameter adjustments */
    --sx;

    /* Function Body */

/*     euclidean norm of the n-vector stored in sx() with storage */
/*     increment incx . */
/*     if    n .le. 0 return with result = 0. */
/*     if n .ge. 1 then incx must be .ge. 1 */

/*           c.l.lawson, 1978 jan 08 */

/*     four phase method     using two built-in constants that are */
/*     hopefully applicable to all machines. */
/*         cutlo = maximum of  sqrt(u/eps)  over all known machines. */
/*         cuthi = minimum of  sqrt(v)      over all known machines. */
/*     where */
/*         eps = smallest no. such that eps + 1. .gt. 1. */
/*         u   = smallest positive no.   (underflow limit) */
/*         v   = largest  no.            (overflow  limit) */

/*     brief outline of algorithm.. */

/*     phase 1    scans zero components. */
/*     move to phase 2 when a component is nonzero and .le. cutlo */
/*     move to phase 3 when a component is .gt. cutlo */
/*     move to phase 4 when a component is .ge. cuthi/m */
/*     where m = n for x() real and m = 2*n for complex. */

/*     values for cutlo and cuthi.. */
/*     from the environmental parameters listed in the imsl converter */
/*     document the limiting values are as follows.. */
/*     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are */
/*                   univac and dec at 2**(-103) */
/*                   thus cutlo = 2**(-51) = 4.44089e-16 */
/*     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec. */
/*                   thus cuthi = 2**(63.5) = 1.30438e19 */
/*     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec. */
/*                   thus cutlo = 2**(-33.5) = 8.23181d-11 */
/*     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19 */
/*     data cutlo, cuthi / 8.232d-11,  1.304d19 / */
/*     data cutlo, cuthi / 4.441e-16,  1.304e19 / */

    if (*n > 0) {
	goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
    next_fmt = fmt_30;
    sum = zero;
    nn = *n * *incx;
/*                                                 begin main loop */
    i__ = 1;
L20:
    switch (next) {
	case 0: goto L30;
	case 1: goto L50;
	case 2: goto L70;
	case 3: goto L110;
    }
L30:
    if ((r__1 = sx[i__], dabs(r__1)) > cutlo) {
	goto L85;
    }
    next = 1;
    next_fmt = fmt_50;
    xmax = zero;

/*                        phase 1.  sum is zero */

L50:
    if (sx[i__] == zero) {
	goto L200;
    }
    if ((r__1 = sx[i__], dabs(r__1)) > cutlo) {
	goto L85;
    }

/*                                prepare for phase 2. */
    next = 2;
    next_fmt = fmt_70;
    goto L105;

/*                                prepare for phase 4. */

L100:
    i__ = j;
    next = 3;
    next_fmt = fmt_110;
    sum = sum / sx[i__] / sx[i__];
L105:
    xmax = (r__1 = sx[i__], dabs(r__1));
    goto L115;

/*                   phase 2.  sum is small. */
/*                             scale to avoid destructive underflow. */

L70:
    if ((r__1 = sx[i__], dabs(r__1)) > cutlo) {
	goto L75;
    }

/*                     common code for phases 2 and 4. */
/*                     in phase 4 sum is large.  scale to avoid overflow. */

L110:
    if ((r__1 = sx[i__], dabs(r__1)) <= xmax) {
	goto L115;
    }
/* Computing 2nd power */
    r__1 = xmax / sx[i__];
    sum = one + sum * (r__1 * r__1);
    xmax = (r__1 = sx[i__], dabs(r__1));
    goto L200;

L115:
/* Computing 2nd power */
    r__1 = sx[i__] / xmax;
    sum += r__1 * r__1;
    goto L200;


/*                  prepare for phase 3. */

L75:
    sum = sum * xmax * xmax;


/*     for real or d.p. set hitest = cuthi/n */
/*     for complex      set hitest = cuthi/(2*n) */

L85:
    hitest = cuthi / (real) (*n);

/*                   phase 3.  sum is mid-range.  no scaling. */

    i__1 = nn;
    i__2 = *incx;
    for (j = i__; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
	if ((r__1 = sx[j], dabs(r__1)) >= hitest) {
	    goto L100;
	}
/* L95: */
/* Computing 2nd power */
	r__1 = sx[j];
	sum += r__1 * r__1;
    }
    ret_val = sqrt(sum);
    goto L300;

L200:
    i__ += *incx;
    if (i__ <= nn) {
	goto L20;
    }

/*              end of main loop. */

/*              compute square root and adjust for scaling. */

    ret_val = xmax * sqrt(sum);
L300:
    return ret_val;
} /* snrm2_ */

doublereal sasum_(integer *n, real *sx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1, r__2, r__3, r__4, r__5, r__6;

    /* Local variables */
    static integer i__, m, mp1, nincx;
    static real stemp;


/*     takes the sum of the absolute values. */
/*     uses unrolled loops for increment equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --sx;

    /* Function Body */
    ret_val = 0.f;
    stemp = 0.f;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	stemp += (r__1 = sx[i__], dabs(r__1));
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*        code for increment equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 6;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	stemp += (r__1 = sx[i__], dabs(r__1));
/* L30: */
    }
    if (*n < 6) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 6) {
	stemp = stemp + (r__1 = sx[i__], dabs(r__1)) + (r__2 = sx[i__ + 1], 
		dabs(r__2)) + (r__3 = sx[i__ + 2], dabs(r__3)) + (r__4 = sx[
		i__ + 3], dabs(r__4)) + (r__5 = sx[i__ + 4], dabs(r__5)) + (
		r__6 = sx[i__ + 5], dabs(r__6));
/* L50: */
    }
L60:
    ret_val = stemp;
    return ret_val;
} /* sasum_ */

/* Subroutine */ int saxpy_(integer *n, real *sa, real *sx, integer *incx, 
	real *sy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;


/*     constant times a vector plus a vector. */
/*     uses unrolled loop for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*sa == 0.f) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sy[iy] += *sa * sx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sy[i__] += *sa * sx[i__];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
	sy[i__] += *sa * sx[i__];
	sy[i__ + 1] += *sa * sx[i__ + 1];
	sy[i__ + 2] += *sa * sx[i__ + 2];
	sy[i__ + 3] += *sa * sx[i__ + 3];
/* L50: */
    }
    return 0;
} /* saxpy_ */

