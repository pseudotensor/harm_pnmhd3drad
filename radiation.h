//radiation test
extern FILE *rad_file;

//constant
extern FTYPE CCON,radth,erfloor;

//radiation control variables
extern int ifld;
extern FTYPE demax,dermax,dtotmax;

//ICCGAF control 
extern int nmeiter,ks0rad,maxrad,iorad,maxit;
extern FTYPE epsrad,epsme;

//root
extern int nhy;

//timestep variable !also used in nudt
extern int nred;

//radiation boundary 
extern int  liiba[N2M];
extern int  loiba[N2M];
extern int  lijba[N1M];
extern int  lojba[N1M];
extern int  *liib;
extern int  *loib;
extern int  *lijb;
extern int  *lojb;
//energy boundary
extern int  niiba[N2M];
extern int  noiba[N2M];
extern int  nijba[N1M];
extern int  nojba[N1M];
extern int  *niib;
extern int  *noib;
extern int  *nijb;
extern int  *nojb;


// radiation E and e
extern FTYPE era[N3M][N2M][N1M] ;
extern FTYPE erna[N3M][N2M][N1M] ;
extern FTYPE ea[N3M][N2M][N1M] ;
extern FTYPE ena[N3M][N2M][N1M] ;
extern FTYPE (*er)[N2M][N1M] ;
extern FTYPE (*ern)[N2M][N1M] ;
extern FTYPE (*e)[N2M][N1M] ;
extern FTYPE (*en)[N2M][N1M] ;


//radiation dE and de
extern FTYPE dera[N3M][N2M][N1M] ;
extern FTYPE dea[N3M][N2M][N1M] ;
extern FTYPE (*der)[N2M][N1M] ;
extern FTYPE (*de)[N2M][N1M] ;

//fn1 and fn2
extern FTYPE fn1a[N3M][N2M][N1M] ;
extern FTYPE fn2a[N3M][N2M][N1M] ;
extern FTYPE (*fn1)[N2M][N1M] ;
extern FTYPE (*fn2)[N2M][N1M] ;

//dfn1de ...
extern FTYPE dfn1dea[N3M][N2M][N1M] ;
extern FTYPE dfn2dea[N3M][N2M][N1M] ;
extern FTYPE dfn1dera[N3M][N2M][N1M] ;
extern FTYPE dfn2dera[N3M][N2M][N1M] ;
extern FTYPE (*dfn1de)[N2M][N1M] ;
extern FTYPE (*dfn2de)[N2M][N1M] ;
extern FTYPE (*dfn1der)[N2M][N1M] ;
extern FTYPE (*dfn2der)[N2M][N1M] ;

//dv 
extern FTYPE dva[3][3][N3M][N2M][N1M] ;
extern FTYPE (*dv)[3][N3M][N2M][N1M] ;
extern FTYPE divva[N3M][N2M][N1M] ;
extern FTYPE (*divv)[N2M][N1M];

//f 
extern FTYPE fa[2][2][N3M][N2M][N1M] ;
extern FTYPE (*f)[2][N3M][N2M][N1M] ;

//D
extern FTYPE dra[2][N3M][N2M][N1M] ;
extern FTYPE (*dr)[N3M][N2M][N1M] ;
extern FTYPE fra[2][N3M][N2M][N1M] ;
extern FTYPE (*fr)[N3M][N2M][N1M] ;

//opacity 
extern FTYPE kapa[N3M][N2M][N1M] ;
extern FTYPE kapna[N3M][N2M][N1M] ;
extern FTYPE siga[N3M][N2M][N1M] ;
extern FTYPE dkapdea[N3M][N2M][N1M] ;
extern FTYPE (*kap)[N2M][N1M] ;
extern FTYPE (*kapn)[N2M][N1M] ;
extern FTYPE (*sig)[N2M][N1M] ;
extern FTYPE (*dkapde)[N2M][N1M];

//black body
extern FTYPE bba[N3M][N2M][N1M] ;
extern FTYPE bbna[N3M][N2M][N1M] ;
extern FTYPE dbbdea[N3M][N2M][N1M] ;
extern FTYPE (*bb)[N2M][N1M] ;
extern FTYPE (*bbn)[N2M][N1M] ;
extern FTYPE (*dbbde)[N2M][N1M] ;

//pressure
extern FTYPE prea[N3M][N2M][N1M];
extern FTYPE prena[N3M][N2M][N1M];
extern FTYPE dpdea[N3M][N2M][N1M] ;
extern FTYPE (*dpde)[N2M][N1M] ;
extern FTYPE (*pre)[N2M][N1M];
extern FTYPE (*pren)[N2M][N1M];

//work variable
extern FTYPE work1xa[N1M];
extern FTYPE work2xa[N1M];
extern FTYPE work3xa[N1M];
extern FTYPE work4xa[N1M];
extern FTYPE work5xa[N1M];
extern FTYPE work6xa[N1M];
extern FTYPE (*work1x);
extern FTYPE (*work2x);
extern FTYPE (*work3x);
extern FTYPE (*work4x);
extern FTYPE (*work5x);
extern FTYPE (*work6x);
extern FTYPE work1ya[N2M];
extern FTYPE work2ya[N2M];
extern FTYPE work3ya[N2M];
extern FTYPE (*work1y);
extern FTYPE (*work2y);
extern FTYPE (*work3y);

extern FTYPE wa[2*N2M][N1M];
extern FTYPE wc[2*N2M][N1M];
#if(RAD)
extern FTYPE wcg1[14*(N1M)*(N2M)];
#endif
