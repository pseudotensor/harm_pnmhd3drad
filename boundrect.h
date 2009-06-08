#define FORCEOUTFLOW (0)
// whether to force outflow on wsca=-2

#define HAWLEYBOUND (0)
// silly =0 of transverse mag components in outflow


// all this is here since usually not changed and want to recompile quickly if changed
#define INNERBC 0 // whether to include inner zones in bound loop


#define SCABCCORNER 1
#define VECBCCORNER 1

// boundary loops over only 1st layer of bc, bound.c takes care of rest
#if( (BOUNDN2==1)&&(BOUNDN1==1)&&(BOUNDN3==1) )
#define LOOPBOUNDRECT LOOPH3 LOOPH2 LOOPH1
#elif( (BOUNDN3==1)&&(BOUNDN2==1)&&(BOUNDN1==0) )
#define LOOPBOUNDRECT i=0; LOOPH3 LOOPH2
#elif( (BOUNDN3==1)&&(BOUNDN2==0)&&(BOUNDN1==1) )
#define LOOPBOUNDRECT j=0; LOOPH3 LOOPH1
#elif( (BOUNDN3==1)&&(BOUNDN2==0)&&(BOUNDN1==0) )
#define LOOPBOUNDRECT j=0; i=0; LOOPH3
#elif( (BOUNDN3==0)&&(BOUNDN2==1)&&(BOUNDN1==0) )
#define LOOPBOUNDRECT k=0; i=0; LOOPH2
#elif( (BOUNDN3==0)&&(BOUNDN2==0)&&(BOUNDN1==1) )
#define LOOPBOUNDRECT k=0; j=0; LOOPH1
#elif( (BOUNDN3==0)&&(BOUNDN2==1)&&(BOUNDN1==1) )
#define LOOPBOUNDRECT k=0; LOOPH2 LOOPH1
#elif( (BOUNDN3==0)&&(BOUNDN2==0)&&(BOUNDN1==0) ) // then assume not really using this function(say postproc call of combine only), so just define dummy
#define LOOPBOUNDRECT LOOPH3 LOOPH2 LOOPH1
#endif

// boundary loops over only 1st layer of bc, bound.c takes care of rest
#if( (N2>1)&&(N1>1)&&(N3>1) )
#define SKIPCONDITION if( ((i>=0)&&(i<N1))&&((j>=0)&&(j<N2))&&((k>=0)&&(k<N3)) ) i=N1
#elif( (N3>1)&&(N2>1)&&(N1==1) )
#define SKIPCONDITION if( ((j>=0)&&(j<N2))&&((k>=0)&&(k<N3)) ) j=N2
#elif( (N3>1)&&(N2==1)&&(N1>1) )
#define SKIPCONDITION if( ((i>=0)&&(i<N1))&&((k>=0)&&(k<N3)) ) i=N1
#elif( (N3>1)&&(N2==1)&&(N1==1) )
#define SKIPCONDITION if( ((k>=0)&&(k<N3)) ) k=N3
#elif( (N3==1)&&(N2>1)&&(N1==1) )
#define SKIPCONDITION if( ((j>=0)&&(j<N2)) ) j=N2
#elif( (N3==1)&&(N2==1)&&(N1>1) )
#define SKIPCONDITION if( ((i>=0)&&(i<N1)) ) i=N1
#elif( (N3==1)&&(N2>1)&&(N1>1) )
#define SKIPCONDITION if( ((i>=0)&&(i<N1))&&((j>=0)&&(j<N2)) ) i=N1
#elif( (N3==1)&&(N2==1)&&(N1==1) ) // then assume not really using this function(say postproc call of combine only), so just define dummy
#define SKIPCONDITION if( ((i>=0)&&(i<N1))&&((j>=0)&&(j<N2))&&((k>=0)&&(k<N3)) ) i=N1
#endif

