#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif

#define DOTRUEBOUNDARY (1)
#define DOINTERNALBOUNDARY (1)

#define DOBOUNDSCA 1
#define DOBOUNDVEC 1


#define REDUCEBOUND (1) // 1: for speed 0: does unnecessary bounding


#if(REDUCEBOUND)
#if(N1==1)
#define BOUNDN1 0 // 0: don't do N1 layer bound, 1: do
#else
#define BOUNDN1 1 // generally should be 1
#endif
// can use below if make N2==1 and don't care about theta structure.  Assumes many things, so could be nasty, but much faster for 1-D problems
#if(N2==1)
#define BOUNDN2 0 // 0: don't do N2 layer bound, 1: do (0 for speed on pure 1D-1D-vector problems)
#else
#define BOUNDN2 1 // generally this should be 1
#endif

#if(N3==1)
#define BOUNDN3 0 // 0: don't do N3 layer bound, 1: do (0 for speed on pure 1D-1D-vector problems)
#else
#define BOUNDN3 1 // generally this should be 1
#endif
#else
#define BOUNDN1 1
#define BOUNDN2 1
#define BOUNDN3 1
#endif


