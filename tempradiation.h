// define LOOP used in radiation module :one more grid beyond active zone
#define LOOPRC1 LOOPH1
#define LOOPRC2 LOOPH2
#define LOOPRC3 LOOPH3
#define LOOPRC LOOPH
//constant
extern FTYPE C_CON,radth,eps,ks0,maxit

//radiation control variables
extern int ifld

//timestep variable !also used in nudt
extern int nred

//radiation boundary 
extern int  liib(N2M),loib(N2M),lijb(N1M),lojb(N1M)

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
extern FTYPE (*der)[N3M][N2M][N1M] ;
extern FTYPE (*de)[N3M][N2M][N1M] ;

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
extern FTYPE (*dpde)[N2M][N1M];
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
extern FTYPE work1ya[N1M];
extern FTYPE work2ya[N1M];
extern FTYPE work3ya[N1M];
extern FTYPE (*work1y);
extern FTYPE (*work2y);
extern FTYPE (*work3y);

