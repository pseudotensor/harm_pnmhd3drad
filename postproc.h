// ./i/zipperr8.sh ./i/zipperppm.sh (convert for old image combine method)

// can do different parts of either combine/interp as continued
// shouldn't continue with both interp and combine if not in synch already

// combine really only needs sizes and number of cpu

// interp needs all things setup for grid info, or read it in


// globalpp.h is used for a clean postproc util.  Still have to recompile for different size input/output

#include "global.h"
#include "defs.h"
#include "rasterfile.h"

#include <signal.h>
#include <time.h>


#define RUNTYPE 0 // see init.c, same interp as there
#define DIRECTINPUT 2 // see init.c (2 normal)

#define DOCOMBINE 0
#define DOINTERP 1   // must have combined before interp since interp works on combined data only
// interp also can not interp, just make more image types(controlled with global.h)
#define DOEXPAND 0 // whether to change grid size using read in of dump and interp
#define FILETOEXPAND "dump0010.dat" // choice
#define DOSPLIT 0 // whether to split grid into multiple cpu grids
#if(DOEXPAND&&DOSPLIT)
#define FILETOSPLIT "filetosplit.dat" // no choice!
#else
#define FILETOSPLIT "pdump0050.dat" // choice
#endif

// combine picks
#define CDOIMAGE 0
#define CDOPAR 0 // need to do this before for ANY interp stuff when multiple cpus
#define CDOPDUMP 0 // primitive dumps
#define CDODUMP 1
#define CDONPDUMP 0
#define CDOFLOORDUMP 0
#define CDOADUMP 0 // most of the time just 1 file at t=0
#define CDOAVG2D 0
#define CDOFLINEDUMP 0 // whether to combine poloidal field line function(of which contour plot is field lines)
// calcs ( in general should either compute or combine, not both)
#define CDOCALCNPDUMP 0 // (really computes from dumps)
#define CDOCALCFLINEDUMP 0 // whether to compute(from dumps) poloidal field line function(of which contour plot is field lines)
#define CDOCALCDUMP 0 // (really computes from dumps) (not really needed)
#define CDOCALCAVG2D 0 // 1: compute calcs from avg2d data from raw data


// interp picks
#define IDOIMAGE 1
#define IDOPARI 0 // whether to output interp par for images
#define IDOPAR 0 // needed for 2D grid dump plots in sm
#define IDOPDUMP 0 // primitive dumps
#define IDODUMP 0
#define IDONPDUMP 0 // should do npdump if doing dump
#define IDOCALCDUMP 0 // needed for hard interp stuff
#define IDOCALCAVG2D 0 // 1: compute calcs from avg2d data from interp data
#define IDOFLOORDUMP 0
#define IDOADUMP 0 // most of the time just 1 file at t=0
#define IDOAVG2D 0
#define IDOFLINEDUMP 0 // whether to compute poloidal field line function(of which contour plot is field lines)


 // number of dump to do only, instead of all.  -1 for all -2 for none
#define ONEPDUMPONLY -1
#define ONEDUMPONLY 43
#define ONENPDUMPONLY -1
#define ONECALCDUMPONLY -1
#define ONEFLINEDUMPONLY -1 // pretty expensive, so only 1 file normally
#define ONEFLOORDUMPONLY -1
#define ONEADUMPONLY -1

#define STARTIc numimagesppc  // starting image #
#define ENDIc (numimages) // ending image # plus 1
//#define ENDIc (1237-1) // ending image #
//#define STARTIc 264
//#define ENDIc 267

//#define STARTIi numimagesppi
//#define ENDIi (numimages)// ending image # plus 1
#define STARTIi 713
#define ENDIi  714 // ending image # plus 1


#define OLDIMAGEMETHOD 0 // 1: use old method of combining images using cat for raw r8s and convert for ppm(gz) 0: use new method of doing like dumps for all kinds


// be careful that code works in making combined data before removing!
#define REMOVECPUIMAGES 0 // 1: remove old multiple cpu images once done with appending 0: don't remove

#define IDUMPOUTTYPE 1 // for dumps... see init_dx outtype for types and detail modifications
#define IIMGOUTTYPE 1 // for images... see init_dx outtype for types and detail modifications
#if(!GAMMIEIMAGE)
#define TILEDIMAGE 0 // whether output images are tiled in x=1 y=2 or z=3 direction(z is normal) (for combine)
// not really setup for combining anything but z(3) and ncpux3=1
#else
#define TILEDIMAGE 0 // assume gammie is 2D right now
#endif
#define SLICENUMBER 1 // number of slices

// input zip format
#define GZIP 0 // should keep at 0, won't hurt if files already zipped

#define DOCONVERTR8 0 // whether to convert final r8 combines into ppms

// can make different than reality if just combining, then set global.h for smaller size
//#define N1PP 512
//#define N2PP 512
//#define N3PP 512

#define N1PP N1
#define N2PP N2
#define N3PP N3


#define N1PPM (N1PP+N1BND*2)
#define N2PPM (N2PP+N2BND*2)
#define N3PPM (N3PP+N3BND*2)
#define N1PPBND (N1BND)
#define N2PPBND (N2BND)
#define N3PPBND (N3BND)

#define DUMN1PP (DUMN1)
#define DUMN2PP (DUMN2)
#define DUMN3PP (DUMN3)
#define IMGN1PP (IMGN1)
#define IMGN2PP (IMGN2)
#define IMGN3PP (IMGN3)
