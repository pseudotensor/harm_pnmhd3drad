
/* generate opacity table from power-law fits */

/*

$Id: opac.c,v 3.2 2007/02/23 15:35:01 gammie Exp $

*/
#include "f2c.h"
#include <stdio.h>
#include <math.h>
double opac(double rho, double T)
{
	real lop, lrho, lT;
	if (rho < 0. || T < 0.) {
		fprintf(stderr, "error: rho or T negative\n");
		fprintf(stderr, "rho: %g T: %g\n", rho, T);
		myexit(1);
	}
	if (rho == 0 || T == 0)
		return (1.);
	lrho = log(rho);
	lT = log(T);
	if (lT < log(166.810)) {

		/* ice grains */
		lop = -8.517193191 + 2. * lT;
	}

	else if (lT < log(202.67678165)) {

		/* ice sublimation */
		lop = 37.53450867 - 7. * lT;
	}

	else if (lT < log(2286.78) + (2. / 49.) * lrho) {

		/* metal grain */
		lop = -2.302585093 + 0.5 * lT;
	}

	else if (lT < log(2029.75) + (1. / 81.) * lrho) {

		/* metal grain evap */
		lop = 187.2025397 + lrho - 24. * lT;
	}

	else if (lT < log(10000.) + (1. / 21.) * lrho) {

		/* molecules */
		lop = -18.42068074 + (2. / 3.) * lrho + 3. * lT;
	}

	else if (lT < log(31195.) + (4. / 75.) * lrho) {

		/* H-scattering */
		lop = -82.89306335 + (1. / 3.) * lrho + 10. * lT;
	}

	else if (lT < log(1.79395e8) + (2. / 5.) * lrho) {

		/* bound-free, free-free */
		lop = 46.45716697 + lrho - (5. / 2.) * lT;
	}

	else {

		/* electron scattering */
		lop = -1.055552799;
	}
	return (exp(lop));
}

