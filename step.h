// add source terms due to grid expansion and contraction 



#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif


#define MOC2DVER 1
// 0: old version which needs bz_sweepx/y and bz transport
// 1: new version all in moc_ct non-split or split

#define MOC3DVER (SPLITMETHOD) // nonsplit typically
// 0: non-split (requires lowering of courant factor to do rotated alfven waves, but have little point shocks, fairly nondiffusive, requires courant=0.5/NDIM)
// 1: split (does rotated alfven waves without lowering courant factor, but has large point shocks, and is very diffusive in constancy and diffusivity convergence rate is poor, but requires courant=0.5)

// below applies to MHD transport routines
#define SYMFORCEMAG 0 // whether to force symmetry on dq in mag stuff(moc_ct, bzsweepx/y)
// RHOINTERP==1 doesn't have strong SYMFORCEMAG in, waste since won't use RHOINTERP=1 at all probably
#define SYMGOOD 1
#define SYMFUNKED 0
#define SYMFUNKEDM1 0
#define SYMFUNKEDM2 0

#define RHOINTERP 0 // used for bzsweepx,bzsweepy, and moc_ct
// 0: use single rho
// 1: use approximate time centered rho
// 2: iterate to get interpolated rho (rhoh3 or rho should really be at midpoint of characteristic in time/space).

#define SNCODE 1 // 0: paper 1: code version

// some loops touch corner zones and don't use that data later, assumes no floating point error data in there


