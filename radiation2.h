//radiation test
FILE *rad_file;
//constant
FTYPE CCON,radth,erfloor;

//radiation control variables
int ifld;
FTYPE demax,dermax,dtotmax;

//ICCGAF control 
int nmeiter,ks0rad,maxrad,iorad,maxit;
FTYPE epsrad,epsme;

//root
int nhy;

//timestep variable !also used in nudt
int nred;

//radiation boundary 
int  liiba[N2M];
int  loiba[N2M];
int  lijba[N1M];
int  lojba[N1M];
int  *liib;
int  *loib;
int  *lijb;
int  *lojb;
//energy boundary
int  niiba[N2M];
int  noiba[N2M];
int  nijba[N1M];
int  nojba[N1M];
int  *niib;
int  *noib;
int  *nijb;
int  *nojb;


// radiation E and e
FTYPE era[N3M][N2M][N1M] ;
FTYPE erna[N3M][N2M][N1M] ;
FTYPE ea[N3M][N2M][N1M] ;
FTYPE ena[N3M][N2M][N1M] ;
FTYPE (*er)[N2M][N1M] ;
FTYPE (*ern)[N2M][N1M] ;
FTYPE (*e)[N2M][N1M] ;
FTYPE (*en)[N2M][N1M] ;


//radiation dE and de
FTYPE dera[N3M][N2M][N1M] ;
FTYPE dea[N3M][N2M][N1M] ;
FTYPE (*der)[N2M][N1M] ;
FTYPE (*de)[N2M][N1M] ;

//fn1 and fn2
FTYPE fn1a[N3M][N2M][N1M] ;
FTYPE fn2a[N3M][N2M][N1M] ;
FTYPE (*fn1)[N2M][N1M] ;
FTYPE (*fn2)[N2M][N1M] ;

//dfn1de ...
FTYPE dfn1dea[N3M][N2M][N1M] ;
FTYPE dfn2dea[N3M][N2M][N1M] ;
FTYPE dfn1dera[N3M][N2M][N1M] ;
FTYPE dfn2dera[N3M][N2M][N1M] ;
FTYPE (*dfn1de)[N2M][N1M] ;
FTYPE (*dfn2de)[N2M][N1M] ;
FTYPE (*dfn1der)[N2M][N1M] ;
FTYPE (*dfn2der)[N2M][N1M] ;

//dv 
FTYPE dva[3][3][N3M][N2M][N1M] ;
FTYPE (*dv)[3][N3M][N2M][N1M] ;
FTYPE divva[N3M][N2M][N1M] ;
FTYPE (*divv)[N2M][N1M];

//f 
FTYPE fa[2][2][N3M][N2M][N1M] ;
FTYPE (*f)[2][N3M][N2M][N1M] ;

//D
FTYPE dra[2][N3M][N2M][N1M] ;
FTYPE (*dr)[N3M][N2M][N1M] ;
FTYPE fra[2][N3M][N2M][N1M] ;
FTYPE (*fr)[N3M][N2M][N1M] ;

//opacity 
FTYPE kapa[N3M][N2M][N1M] ;
FTYPE kapna[N3M][N2M][N1M] ;
FTYPE siga[N3M][N2M][N1M] ;
FTYPE dkapdea[N3M][N2M][N1M] ;
FTYPE (*kap)[N2M][N1M] ;
FTYPE (*kapn)[N2M][N1M] ;
FTYPE (*sig)[N2M][N1M] ;
FTYPE (*dkapde)[N2M][N1M];

//black body
FTYPE bba[N3M][N2M][N1M] ;
FTYPE bbna[N3M][N2M][N1M] ;
FTYPE dbbdea[N3M][N2M][N1M] ;
FTYPE (*bb)[N2M][N1M] ;
FTYPE (*bbn)[N2M][N1M] ;
FTYPE (*dbbde)[N2M][N1M] ;

//pressure
FTYPE prea[N3M][N2M][N1M];
FTYPE prena[N3M][N2M][N1M];
FTYPE dpdea[N3M][N2M][N1M] ;
FTYPE (*dpde)[N2M][N1M] ;
FTYPE (*pre)[N2M][N1M];
FTYPE (*pren)[N2M][N1M];

//work variable
FTYPE work1xa[N1M];
FTYPE work2xa[N1M];
FTYPE work3xa[N1M];
FTYPE work4xa[N1M];
FTYPE work5xa[N1M];
FTYPE work6xa[N1M];
FTYPE (*work1x);
FTYPE (*work2x);
FTYPE (*work3x);
FTYPE (*work4x);
FTYPE (*work5x);
FTYPE (*work6x);
FTYPE work1ya[N2M];
FTYPE work2ya[N2M];
FTYPE work3ya[N2M];
FTYPE (*work1y);
FTYPE (*work2y);
FTYPE (*work3y);

FTYPE wa[2*N2M][N1M];
FTYPE wc[2*N2M][N1M];
#if(RAD)
FTYPE wcg1[14*(N1M)*(N2M)];
#endif

