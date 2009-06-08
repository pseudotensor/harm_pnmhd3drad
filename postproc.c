// Program to take multiple CPU data and combine into single data set for various kinds of data output including interpolation

// **compile using:
// make -f makefile.pp
//(makefile.pp) (which is just normal makefile with main replaced by postproc and twod by postproc

// **to be run in the data directory**

// **CHECK POSTPROC defines for settable parameters
// for both combining and interp
// 0) number of CPUS for combine
// 1) sizes of each cpu in this code below
// 2) size of combined grid in global.h(normal for combine, interp values for interp) mess up size for combined grid and not good as I don't check that for combining
// for combining:
// for interpolation:
//1)  sampling(forced in global.h for pp)
//2) file output type(forced in global.h for pp)

// needs in the data dir:
#include "postproc.h"

int starti;
int itemp1,itemp2;
FILE* numin;
// determine number of dumps
int im_cnts[ITYPES][NUMSCA+1];
int im_cntv[ITYPES][NUMVEC+1];
int iii,iiito,ii,jj,kk;
int qtip;
int myidc;
int l,m,ll,mm ;
int q;
int outtype;
char dirname[MAXFILENAME*10];
char outname[MAXFILENAME*10];
char fullname[MAXCPUS][MAXFILENAME*20];
char command[MAXFILENAME*100];
char tempcommand[MAXFILENAME*100];
char command2[MAXFILENAME*100];
int firsthit;
int newtype;
int num1d_31,num1d_32;
int avgcount;
float tempfloat;
int dumi[3+1];
SFTYPE dumf[3+1];
int length;
int numpdumpsppc,numdumpsppc,numnpdumpsppc,numfloordumpsppc,numadumpsppc,numimagesppc,numcalcdumpsppc,numfldumpsppc;
int numpdumpsppi,numdumpsppi,numnpdumpsppi,numfloordumpsppi,numadumpsppi,numimagesppi,numcalcdumpsppi,numfldumpsppi;
int i,j,k;
char fileheader[MAXFILENAME*2];
char filename[MAXFILENAME*2];
char temps[MAXFILENAME*2];
char dfnam[MAXFILENAME*2];
FILE * incombine[MAXCPUS];
FILE * in;
FILE * out;
FILE * image_file;
FILE * dump_file;
unsigned char ch,chtemp;
int go,fullgo;
int count;
int dump_cnt;
int im_cnt;
int call_code;
int n1,n2,n3;
int n[3+1][MAXCPUS];

short pal[3][256];
int skip;
int newskip;
FILE *testfp;
int idodumplist[10000];

int llto,prim,qto,pumpsize,primto;
char myidctxt[100];
fpos_t fpos;
FILE * file_expand;
int error;
int tocont;
int doneout;

int imagen1,imagen2,ntile1,ntile2;
int TILEDIR;

int ppncpux1,ppncpux2,ppncpux3,ppnumprocs;

int checktype(FILE * in, FILE * out, int which, int fixout);
void loadpal(void);
void convertr8(int nx,int ny,int gzipout,int gzipin,short (*pal)[256],char* inname,char* outname);
void decint(int *lp, int n);
void r82ras(int nx,int ny,short (*pal)[256],FILE*in,FILE*out);
void r82ppm(int nx,int ny,short (*pal)[256],FILE*in,FILE*out);
void generalcombine(int which);

int main(
         int argc,
         char *argv[],
         char *envp[]
         )
{

  if(!POSTPROC){
    fprintf(stderr,"Forgot to turn on POSTPROC in global.h?\n");
    myexit(1);
  }

  if(argc==2){
    if(argv[1][0]=='c'){
      tocont=1;
    }
    else{
      tocont=0;
    }
  }
  else tocont=0;

  // load global parameter variables
  gpar_init();
  // load cpu geom into pp state
  ppncpux1=ncpux1;
  ppncpux2=ncpux2;
  ppncpux3=ncpux3;
  ppnumprocs=ppncpux1*ppncpux2*ppncpux3;

  // uncomment if get SIGFPE 8 error at run with floats
  //  signal(8,SIG_IGN); 
  
  // this is overwritten if split is called


  // enter as like runtype==3(i.e. full input)
  runtype=RUNTYPE;
  directinput=DIRECTINPUT;
  // must treat like single cpu now (esp for interp)
  numprocs=1;
  ncpux1=1;
  ncpux2=1;
  ncpux3=1;
  myid=0;
  sprintf(myidtxt,"");
  printf("doing init genfiles\n");
  fflush(stdout);
  init_genfiles(1); // init log files
  // assume pp is always on single cpu
  totalsize[1]=N1PP;
  totalsize[2]=N2PP;
  totalsize[3]=N3PP;


  // notify user that sizes don't match in case not expected
  if((N1!=N1PP)||(N2!=N2PP)||(N3!=N3PP)){
    fprintf(stderr,"NOTE: N1=%d != N1PP=%d or N2=%d != N2PP=%d or N3=%d != N3PP=%d\n",N1,N1PP,N2,N2PP,N3,N3PP);
  }

  for(i=0;i<ppnumprocs;i++){
    n[1][i]=N1PP/ppncpux1;
    n[2][i]=N2PP/ppncpux2;
    n[3][i]=N3PP/ppncpux3;
  }



  // cut/pasted from init.c
  // determine number of primitive dumps
  if( (DOCOMBINE&&CDOPDUMP)||(DOINTERP&&IDOPDUMP) ){
    sprintf(temps,"0_numpdumps.pp");
    if( (in=fopen(temps,"rt"))==NULL){
      printf("can't open in pp numpdumps file, assume fresh start\n");
      numdumpsppc=0;
      numdumpsppi=0;
      //myexit(1);
    }
    else{
      while(fgetc(in)!='\n');
      
      fscanf(in,"%d %d",&numpdumpsppc,&numpdumpsppi);
      fclose(in);
    }
    if(tocont==0){
      if(DOCOMBINE&&CDOPDUMP){
	numpdumpsppc=0;
      }
      if(DOINTERP&&IDOPDUMP){
	numpdumpsppi=0;
      }
    }
    printf("number of ppc pdumps: %d ppi pdumps: %d\n",numpdumpsppc,numpdumpsppi);
  }
  // determine number of normal dump files
  if( (DOCOMBINE&&CDODUMP)||(DOINTERP&&IDODUMP) ){
    sprintf(temps,"0_numdumps.pp");
    if( (in=fopen(temps,"rt"))==NULL){
      printf("can't open pp numdumps file, assume fresh start\n");
      numdumpsppc=0;
      numdumpsppi=0;
      //myexit(1);
    }
    else{
      while(fgetc(in)!='\n');
      
      fscanf(in,"%d %d",&numdumpsppc,&numdumpsppi);
      fclose(in);
    }
    if(tocont==0){
      if(DOCOMBINE&&CDODUMP){
	numdumpsppc=0;	
      }
      if(DOINTERP&&IDODUMP){
	numdumpsppi=0;
      }
    }
    printf("number of ppc dumps: %d ppi dumps: %d\n",numdumpsppc,numdumpsppi);
  }
  if( (DOCOMBINE&&CDONPDUMP)||(DOINTERP&&IDONPDUMP) ){
    sprintf(temps,"0_numnpdumps.pp");
    if( (in=fopen(temps,"rt"))==NULL){
      printf("can't open pp numnpdumps file, assuming fresh start\n");
      numnpdumpsppc=0;
      numnpdumpsppi=0;
      //myexit(1);
    }
    else{
      while(fgetc(in)!='\n');
      
      fscanf(in,"%d %d",&numnpdumpsppc,&numnpdumpsppi);
      fclose(in);
    }
    if(tocont==0){
      if(DOCOMBINE&&CDONPDUMP){
	numnpdumpsppc=0;
      }
      if(DOINTERP&&IDONPDUMP){
	numnpdumpsppi=0;
      }
    }
    printf("number of ppc np dumps: %d  ppi np dumps: %d\n",numnpdumpsppc,numnpdumpsppi);
  }

  if( (DOCOMBINE&&CDOCALCDUMP)||(DOINTERP&&IDOCALCDUMP) ){
    sprintf(temps,"0_numcalcdumps.pp");
    if( (in=fopen(temps,"rt"))==NULL){
      printf("can't open pp numcalcdumps file, assuming fresh start\n");
      numcalcdumpsppc=0;
      numcalcdumpsppi=0;
      //myexit(1);
    }
    else{
      while(fgetc(in)!='\n');
      
      fscanf(in,"%d %d",&numcalcdumpsppc,&numcalcdumpsppi);
      fclose(in);
    }
    if(tocont==0){
      if(DOCOMBINE&&CDOCALCDUMP){
	numcalcdumpsppc=0;
      }
      if(DOINTERP&&IDOCALCDUMP){
	numcalcdumpsppi=0;
      }
    }
    printf("number of ppc calc dumps: %d  ppi np dumps: %d\n",numcalcdumpsppc,numcalcdumpsppi);
  }

  if( (DOCOMBINE&&CDOFLINEDUMP)||(DOINTERP&&IDOFLINEDUMP) ){
    sprintf(temps,"0_numfldumps.pp");
    if( (in=fopen(temps,"rt"))==NULL){
      printf("can't open pp numfldumps file, assuming fresh start\n");
      numfldumpsppc=0;
      numfldumpsppi=0;
      //myexit(1);
    }
    else{
      while(fgetc(in)!='\n');
      
      fscanf(in,"%d %d",&numfldumpsppc,&numfldumpsppi);
      fclose(in);
    }
    if(tocont==0){
      if(DOCOMBINE&&CDOFLINEDUMP){
	numfldumpsppc=0;
      }
      if(DOINTERP&&IDOFLINEDUMP){
	numfldumpsppi=0;
      }
    }
    printf("number of ppc fieldline dumps: %d  ppi np dumps: %d\n",numfldumpsppc,numfldumpsppi);
  }

  if( (DOCOMBINE&&CDOADUMP)||(DOINTERP&&IDOADUMP) ){
    sprintf(temps,"0_numadumps.pp");
    if( (in=fopen(temps,"rt"))==NULL){
      printf("can't open pp numadumps file, assuming fresh start\n");
      numadumpsppc=0;
      numadumpsppi=0;
      //myexit(1);
    }
    else{
      while(fgetc(in)!='\n');
      
      fscanf(in,"%d %d",&numadumpsppc,&numadumpsppi);
      fclose(in);
    }
    if(tocont==0){
      if(DOCOMBINE&&CDOADUMP){
	numadumpsppc=0;
      }
      if(DOINTERP&&IDOADUMP){
	numadumpsppi=0;
      }
    }
    printf("number of ppc adumps: %d ppi adumps: %d\n",numadumpsppc,numadumpsppi);
  }

  if( (DOCOMBINE&&CDOFLOORDUMP)||(DOINTERP&&IDOFLOORDUMP) ){
    sprintf(temps,"0_numfloordumps.pp");
    if( (in=fopen(temps,"rt"))==NULL){
      printf("can't open pp numfloordump file, assuming fresh start\n");
      numfloordumpsppc=0;
      numfloordumpsppi=0;
      //myexit(1);
    }
    else{
      while(fgetc(in)!='\n');
      
      fscanf(in,"%d %d",&numfloordumpsppc,&numfloordumpsppi);
      fclose(in);
    }
    if(tocont==0){  
      if(DOCOMBINE&&CDOFLOORDUMP){
	numfloordumpsppc=0;
      }
      if(DOINTERP&&IDOFLOORDUMP){
	numfloordumpsppi=0;
      }
    }
    printf("number of ppc floor dumps: %d ppi floor dumps: %d\n",numfloordumpsppc,numfloordumpsppi);
  }

  if( (DOCOMBINE&&CDOIMAGE)||(DOINTERP&&IDOIMAGE) ){
    // get # of images
    sprintf(temps,"%s%s",DATADIR,IMAGEDIR);	
    sprintf(dfnam,"%s0_numimages%s",temps,PPEXT);
    if((image_file = fopen(dfnam,"r"))==NULL) {
      fprintf(fail_file,"error opening pp image  file %s, assuming fresh start\n",dfnam) ;
      numimagesppc=0;
      numimagesppi=0;
      //myexit(1) ;
    }
    else{
      while(fgetc(image_file)!='\n'); // skip comment line
      fscanf(image_file,"%d %d",&numimagesppc,&numimagesppi);
      fclose(image_file);
    }
    if(tocont==0){
      if(DOCOMBINE&&CDOIMAGE){
	numimagesppc=0;
      }
      if(DOINTERP&&IDOIMAGE){
	numimagesppi=0;
      }
    }
    printf("number of ppc images: %d ppi images: %d\n",numimagesppc,numimagesppi);
  }
  //
  ///////////

  /////////
  // now for current number
  if( (in=fopen("0_numpdumps.dat","rt"))==NULL){
    printf("can't open in numpdumps file\n");
    //    myexit(1);
    // assume # is 0
    numpdumps=0;
  }
  else{
    while(fgetc(in)!='\n');
    fscanf(in,"%d",&numpdumps);
    fclose(in);
  }
  printf("number of pdumps: %d\n",numpdumps);

  // determine number of normal dump files
  sprintf(dfnam,"%s%s0_numdumps.dat",DATADIR,DUMPDIR);

  if( (in=fopen(dfnam,"rt"))==NULL){
    printf("can't open numdumps file\n");
    //myexit(1);
    numdumps=0;
  }
  else{
    while(fgetc(in)!='\n');  
    fscanf(in,"%d",&numdumps);
    fclose(in);
  }
  printf("number of dumps: %d\n",numdumps);
  //if(!GAMMIEIMAGE){
    if( (!(DOCOMBINE&&CDODUMP))&&(DOINTERP&&IDODUMP) ){
      numdumps=numdumpsppc; // don't go over what's been combined
    }
    //}
  if(0){
    if( (in=fopen("0_numnpdumps.dat","rt"))==NULL){
      printf("can't open numnpdumps file\n");
      //      myexit(1);
      numnpdumps=0;
    }
    else{
      while(fgetc(in)!='\n');
      fscanf(in,"%d",&numnpdumps);
      fclose(in);
    }  
    printf("number of np dumps: %d\n",numnpdumps);
    if( (!(DOCOMBINE&&CDONPDUMP))&&(DOINTERP&&IDONPDUMP) ){
      numnpdumps=numnpdumpsppc; // don't go over what's been combined
    }
  }
  else{
    numnpdumps=numdumps; // now compute npdumps from dumps
  }
  if(0){
    if( (in=fopen("0_numcalcdumps.dat","rt"))==NULL){
      printf("can't open numcalcdumps file\n");
      //      myexit(1);
      numcalcdumps=0;
    }
    else{
      while(fgetc(in)!='\n');
      fscanf(in,"%d",&numcalcdumps);
      fclose(in);
    }  
    printf("number of calc dumps: %d\n",numcalcdumps);
    if( (!(DOCOMBINE&&CDOCALCDUMP))&&(DOINTERP&&IDOCALCDUMP) ){
      numcalcdumps=numcalcdumpsppc; // don't go over what's been combined
    }
  }
  else{
    numcalcdumps=numdumps; // now compute calcdumps from dumps
  }

  if(DOCOMBINE&&CDOFLINEDUMP){
    if( (in=fopen("0_numfldumps.dat","rt"))==NULL){
      printf("can't open numfldumps file\n");
      //      myexit(1);
      numfldumps=0;
    }
    else{
      while(fgetc(in)!='\n');
      fscanf(in,"%d",&numfldumps);
      fclose(in);
    }  
    printf("number of fieldline dumps: %d\n",numfldumps);
    if( (!(DOCOMBINE&&CDOFLINEDUMP))&&(DOINTERP&&IDOFLINEDUMP) ){
      numfldumps=numfldumpsppc; // don't go over what's been combined
    }
  }
  else{
    numfldumps=numdumps; // now compute calcdumps from dumps
  }

  if( (in=fopen("0_numadumps.dat","rt"))==NULL){
    printf("can't open numadumps file\n");
    //    myexit(1);
    numadumps=0;
  }
  else{
    while(fgetc(in)!='\n');  
    fscanf(in,"%d",&numadumps);
    fclose(in);
  }
  printf("number of adumps: %d\n",numadumps);
  if( (!(DOCOMBINE&&CDOADUMP))&&(DOINTERP&&IDOADUMP) ){
    numadumps=numadumpsppc; // don't go over what's been combined
  }
  
  if( (in=fopen("0_numfloordumps.dat","rt"))==NULL){
    printf("can't open numfloordump file\n");
    //    myexit(1);
    numfloordumps=0;
  }
  else{
    while(fgetc(in)!='\n');  
    fscanf(in,"%d",&numfloordumps);
    fclose(in);
  }
  printf("number of floor dumps: %d\n",numfloordumps);
  if( (!(DOCOMBINE&&CDOFLOORDUMP))&&(DOINTERP&&IDOFLOORDUMP) ){
    numfloordumps=numfloordumpsppc; // don't go over what's been combined
  }


  // get # of images
  sprintf(temps,"%s%s",DATADIR,IMAGEDIR);	
  sprintf(dfnam,"%s0_numimages%s",temps,DATEXT);
  if((in = fopen(dfnam,"r"))==NULL) {
    fprintf(fail_file,"error opening dump output file %s\n",dfnam) ;
    //    myexit(1) ;
    numimages=0;
  }
  else{
    while(fgetc(in)!='\n'); // skip comment line
    fscanf(in,"%d",&numimages);
    fclose(in);
  }
  printf("number of images: %d\n",numimages);
  //  if(!GAMMIEIMAGE){
    if( (!(DOCOMBINE&&CDOIMAGE))&&(DOINTERP&&IDOIMAGE) ){
      numimages=numimagesppc; // don't go over what's been combined
    }
    // }



  if( (DOCOMBINE&&CDOPDUMP)||(DOINTERP&&IDOPDUMP) ){
    if( (out=fopen("0_numpdumps.pp","wt"))==NULL){
      printf("can't open out numpdumps file\n");
      myexit(1);
    }
    if(DOCOMBINE&&CDOPDUMP) itemp1=numpdumps; else itemp1=numpdumpsppc;
    if(DOINTERP&&IDOPDUMP) itemp2=numpdumps; else itemp2=numpdumpsppi;
    fprintf(out,"#number of pdumps:\n%d %d\n",itemp1,itemp2);
    fclose(out);
  }

  if( (DOCOMBINE&&CDODUMP)||(DOINTERP&&IDODUMP) ){
    if( (out=fopen("0_numdumps.pp","wt"))==NULL){
      printf("can't open pp numdumps file\n");
      myexit(1);
    }
    if(DOCOMBINE&&CDODUMP) itemp1=numdumps; else itemp1=numdumpsppc;
    if(DOINTERP&&IDODUMP) itemp2=numdumps; else itemp2=numdumpsppi;
    fprintf(out,"#number of dumps:\n%d %d\n",itemp1,itemp2);
    fclose(out);
  }

  if( (DOCOMBINE&&CDONPDUMP)||(DOINTERP&&IDONPDUMP) ){
    if( (out=fopen("0_numnpdumps.pp","wt"))==NULL){
      printf("can't open pp numnpdumps file\n");
      myexit(1);
    }
    if(DOCOMBINE&&CDONPDUMP) itemp1=numnpdumps; else itemp1=numnpdumpsppc;
    if(DOINTERP&&IDONPDUMP) itemp2=numnpdumps; else itemp2=numnpdumpsppi;
    fprintf(out,"#number of np dumps:\n%d %d\n",itemp1,itemp2);
    fclose(out);
  }

  if( (DOCOMBINE&&CDOCALCDUMP)||(DOINTERP&&IDOCALCDUMP) ){
    if( (out=fopen("0_numcalcdumps.pp","wt"))==NULL){
      printf("can't open pp numcalcdumps file\n");
      myexit(1);
    }
    if(DOCOMBINE&&CDOCALCDUMP) itemp1=numcalcdumps; else itemp1=numcalcdumpsppc;
    if(DOINTERP&&IDOCALCDUMP) itemp2=numcalcdumps; else itemp2=numcalcdumpsppi;
    fprintf(out,"#number of calc dumps:\n%d %d\n",itemp1,itemp2);
    fclose(out);
  }

  if( (DOCOMBINE&&CDOFLINEDUMP)||(DOINTERP&&IDOFLINEDUMP) ){
    if( (out=fopen("0_numfldumps.pp","wt"))==NULL){
      printf("can't open pp numfldumps file\n");
      myexit(1);
    }
    if(DOCOMBINE&&CDOFLINEDUMP) itemp1=numfldumps; else itemp1=numfldumpsppc;
    if(DOINTERP&&IDOFLINEDUMP) itemp2=numfldumps; else itemp2=numfldumpsppi;
    fprintf(out,"#number of fieldline dumps:\n%d %d\n",itemp1,itemp2);
    fclose(out);
  }

  if( (DOCOMBINE&&CDOADUMP)||(DOINTERP&&IDOADUMP) ){
    if( (out=fopen("0_numadumps.pp","wt"))==NULL){
      printf("can't open pp numadumps file\n");
      myexit(1);
    }
    if(DOCOMBINE&&CDOADUMP) itemp1=numadumps; else itemp1=numadumpsppc;
    if(DOINTERP&&IDOADUMP) itemp2=numadumps; else itemp2=numadumpsppi;
    fprintf(out,"#number of adumps:\n%d %d\n",itemp1,itemp2);
    fclose(out);
  }

  if( (DOCOMBINE&&CDOFLOORDUMP)||(DOINTERP&&IDOFLOORDUMP) ){
    if( (out=fopen("0_numfloordumps.pp","wt"))==NULL){
      printf("can't open numfloordump file\n");
      myexit(1);
    }
    if(DOCOMBINE&&CDOFLOORDUMP) itemp1=numfloordumps; else itemp1=numfloordumpsppc;
    if(DOINTERP&&IDOFLOORDUMP) itemp2=numfloordumps; else itemp2=numfloordumpsppi;
    fprintf(out,"number of floor dumps:\n%d %d\n",itemp1,itemp2);
    fclose(out);
  }

  if( (DOCOMBINE&&CDOIMAGE)||(DOINTERP&&IDOIMAGE) ){
    sprintf(temps,"%s%s",DATADIR,IMAGEDIR);	
    sprintf(dfnam,"%s0_numimages%s",temps,PPEXT);
    if((out = fopen(dfnam,"w"))==NULL) {
      fprintf(fail_file,"error opening pp dump output file %s\n",dfnam) ;
      myexit(1) ;
    }
    if(DOCOMBINE&&CDOIMAGE) itemp1=numimages; else itemp1=numimagesppc;
    if(DOINTERP&&IDOIMAGE) itemp2=numimages; else itemp2=numimagesppi;
    fprintf(out,"#number of images:\n%d %d\n",itemp1,itemp2);
    fclose(out);
  }

            

  
  if(DOCOMBINE){





    if(CDOIMAGE){

      // currently assume all cpus same size/tilesize
      myidc=0;
      if(TILEDIMAGE){
	// copied from image()
	TILEDIR=TILEDIMAGE;
	ntile1=(int)ceil(sqrt((FTYPE)(n[TILEDIR][myidc]))-1E-6);
	ntile2=(int)floor(sqrt((FTYPE)(n[TILEDIR][myidc]))+1E-6);
	n1=imagen1=ntile1*n[TILEDIR%3+1][myidc];
	n2=imagen2=ntile2*n[3-(4-TILEDIR)%3][myidc];
      }
      else{
	// GODMARK (should generalize this to be consistent/connect with image()
	ntile1=ntile2=1;
	n1=imagen1=n[1][myidc];
	n2=imagen2=n[2][myidc];
      }
      fprintf(stdout,"DEBUG: %d %d %d . %d %d . %d %d\n",n[1][0],n[2][0],n[3][0],ntile1,ntile2,n1,n2); fflush(stdout);

      
      printf("Combining images#: %d\n",numimages);
      fflush(stdout);
      
      sprintf(dirname,"");	
      
      firsthit=1;
      for(qtip=STARTIc;qtip<ENDIc;qtip++){
	
	printf("combining image# %d ",qtip);
	fflush(stdout);
	
	if(!GAMMIEIMAGE) primto=2; else primto=1;

	for(prim=1;prim<=primto;prim++){ // go over scalars and vectors
	  if(!GAMMIEIMAGE){
	    if(prim==1){
	      llto=NUMSCA;
	      qto=1;
	      iiito=1;
	    }
	    else if(prim==2){
	      llto=REALNUMVEC;
	      qto=3;
	      iiito=1;
	    }
	  }
	  else{
	    llto=1; 
	    qto=1;
	    iiito=9;// all primitives log and linear of 2
	  }
	  
	  
	  for(outtype=0;outtype<=0;outtype++){ // force only first view format
	    
	    for(ll=1;ll<=llto;ll++){
	      im_cnts[outtype][ll]=qtip;
	      for(iii=0;iii<=iiito;iii++){ // do both linear and log10 scalars
		if((prim==2)&&(iii==1)) break; // don't do rho*v combine 
		for(q=1;q<=qto;q++){ // no default write of v0
		  for(k=0;k<SLICENUMBER;k++){
		    if(!GAMMIEIMAGE){
		      if(!OLDIMAGEFORMAT){
			if(prim==1){
			  sprintf(dirname,"./%simx%01d-%01d-%01d-s%01d/",IMAGEDIR,outtype,iii,k,ll);
			}
			else{
			  sprintf(dirname,"./%simx%01d-%01d-%01d-v%01d-%01d/",IMAGEDIR,outtype,iii,k,ll,q);
			}
		      }
		      else{
			if(prim==1){
			  sprintf(dirname,"./%simx%01d-%01d-s%01d/",IMAGEDIR,outtype,iii,ll);
			}
			else{
			  sprintf(dirname,"./%simx%01d-%01d-v%01d-%01d/",IMAGEDIR,outtype,iii,ll,q);
			}
		      }
		    }
		    else{
		      sprintf(dirname,"./%s",IMAGEDIR);
		    }
		    // gzip all image files in this directory by default.  If already gzipped, won't hurt
		    if( (firsthit==1)&&(GZIP==0)){
		      if(!GAMMIEIMAGE){
			printf("gzipping %d %d %d %d %d %d ... ",prim,outtype,ll,iii,k,q);
			fflush(stdout);
			
			if(IMAGEFORMATINPUT==1){
			  sprintf(command,"./%szipperppm.sh %s",IMAGEDIR,dirname);
			  system(command);
			}
			if(IMAGEFORMATINPUT==0){
			  sprintf(command,"./%szipperr8.sh %s",IMAGEDIR,dirname);
			  system(command);
			}
		      }
		      else{
			if(iii==0){ // all in same dir, so just do once
			  printf("gzipping  ... ");
			  fflush(stdout);
			  
			  sprintf(command,"./%szipperrawr8.sh %s",IMAGEDIR,IMAGEDIR);
			  system(command);
			}
		      }
		    }
		    if(!GAMMIEIMAGE){
		      if(!OLDIMAGEFORMAT){
			if(prim==1){
			  sprintf(outname," %simx%01d-%01d-%01d-s%01d-%04d%s",dirname,outtype,iii,k,ll,im_cnts[outtype][ll],DATEXT);
			}
			else{
			  sprintf(outname," %simx%01d-%01d-%01d-v%01d-%01d-%04d%s",dirname,outtype,iii,k,ll,q,im_cnts[outtype][ll],DATEXT);
			}
		      }
		      else{
			if(prim==1){
			  sprintf(outname," %simx%01d-%01d-s%01d-%04d%s",dirname,outtype,iii,ll,im_cnts[outtype][ll],DATEXT);
			}
			else{
			  sprintf(outname," %simx%01d-%01d-v%01d-%01d-%04d%s",dirname,outtype,iii,ll,q,im_cnts[outtype][ll],DATEXT);
			}
		      }
		      if(IMAGEFORMATINPUT==1){
			strcat(outname,".ppm.gz");
		      }
		      else{
			strcat(outname,".r8.gz");
		      }
		    }
		    else{
		      sprintf(outname," %sim%01dp%04d",dirname,iii,im_cnts[outtype][ll]);
		      strcat(outname,".r8.gz"); // assume gzipped output
		    }
		    sprintf(command2,"rm ");
		    
		    
		    
		    
		    // check to see if all files exist
		    skip=0;
		    for(myidc=0;myidc<ppnumprocs;myidc++){
		      sprintf(myidctxt,CPUTXT,myidc);
		      if(!GAMMIEIMAGE){
			if(!OLDIMAGEFORMAT){
			  if(prim==1){
			    sprintf(fullname[myidc],"%simx%01d-%01d-%01d-s%01d-%04d%s%s",dirname,outtype,iii,k,ll,im_cnts[outtype][ll],DATEXT,myidctxt);
			  }
			  else{
			    sprintf(fullname[myidc],"%simx%01d-%01d-%01d-v%01d-%01d-%04d%s%s",dirname,outtype,iii,k,ll,q,im_cnts[outtype][ll],DATEXT,myidctxt); 
			  }
			}
			else{
			  if(prim==1){
			    sprintf(fullname[myidc],"%simx%01d-%01d-s%01d-%04d%s%s",dirname,outtype,iii,ll,im_cnts[outtype][ll],DATEXT,myidctxt);
			  }
			  else{
			    sprintf(fullname[myidc],"%simx%01d-%01d-v%01d-%01d-%04d%s%s",dirname,outtype,iii,ll,q,im_cnts[outtype][ll],DATEXT,myidctxt); 
			  }
			}
			if(IMAGEFORMATINPUT==1) strcat(fullname[myidc],".ppm.gz");
			else if(IMAGEFORMATINPUT==0) strcat(fullname[myidc],".r8.gz");
		      }
		      else{
			sprintf(fullname[myidc],"%sim%01dp%04d%s",dirname,iii,im_cnts[outtype][ll],myidctxt);
			strcat(fullname[myidc],".r8.gz"); // assume gzipped first
		      }
		      
		      // test for existence of multiple cpu data
		      testfp=fopen(fullname[myidc],"rb");
		      if(testfp==NULL) skip++;
		      else fclose(testfp);
		      
		      strcat(command2,fullname[myidc]);
		      strcat(command2," ");
		    }
		    
		    // now you know whether all files exist, loop over them and append them
		    if(skip==0){
		      // create output file
		      sprintf(command,"gzip > %s",outname);
		      out = popen(command,"w");
		      if(out==NULL){
			fprintf(fail_file,"error opening image file: %s\n",command) ;
			myexit(1) ;
		      }
		      
		      // loop over cpus
		      for(myidc=0;myidc<ppnumprocs;myidc++){
			
			sprintf(command,"gzip -d < %s",fullname[myidc]);
			incombine[myidc] = popen(command,"r");
			if(incombine[myidc]==NULL){
			  fprintf(fail_file,"error opening image file: %s\n",command) ;
			  myexit(1) ;
			}
		      }
		      // now all cpus files are ready for this image

		      if(OLDSCHOOL2==0){
			myidc=0; // pipe header first
			// pipe/skip comment lines
			for(i=0;i<=1;i++){ // pipe first 2 lines
			  while(1){
			    ch=fgetc(incombine[myidc]);
			    fputc(ch,out);
			    if(ch=='\n') break;
			  }
			}
			// now skip image size line
			while(fgetc(incombine[myidc])!='\n'); // skip to next line
			// write out own line
			fprintf(out,"%i %i\n", imagen1*ppncpux1, imagen2*ppncpux2); // not IMGN1PP/etc because no interp yet ever!
			while(1){
			  ch=fgetc(incombine[myidc]);
			  fputc(ch,out);
			  if(ch=='\n') break;
			}
			if(IMAGEFORMATINPUT==1){ // now need to pass palette
			  for(i=1;i<=256*3;i++){
			    ch=fgetc(incombine[myidc]);
			    fputc(ch,out);
			  }
			}
			// done with passing header info
		      }
		      // for other cpu images just skip header
		      // now forward all other cpus
		      for(myidc=1;myidc<ppnumprocs;myidc++){
			if(OLDSCHOOL2==0){
			  while(fgetc(incombine[myidc])!='\n'); // skip to next line
			  while(fgetc(incombine[myidc])!='\n'); // skip to next line
			  while(fgetc(incombine[myidc])!='\n'); // skip to next line
			  while(fgetc(incombine[myidc])!='\n'); // skip to next line
			  if(IMAGEFORMATINPUT==1){
			    fseek(incombine[myidc],SEEK_CUR,256*3);
			  }
			}
		      }
		      // now at image data(just pipe from in to out)
		      if(!GAMMIEIMAGE){
			if(IMAGEFORMATINPUT==1) pumpsize=3;
			else if(IMAGEFORMATINPUT==0) pumpsize=1;
		      }
		      else pumpsize=1;

		      if(1||(!GAMMIEIMAGE)){
			// pipe the columns of this row
			for(k=1;k<=ntile2;k++){ // do each tile up/down (assumes all cpus are same size)
			  for(kk=0;kk<ppncpux2;kk++){ // do each cpu up/down (assumes all cpus are same size)
			    for(jj=1;jj<=n[2][0];jj++){ // each line up/down
			      for(j=1;j<=ntile1;j++){ // do each tile left/right
				for(myidc=kk*ppncpux1;myidc<ppncpux1+kk*ppncpux1;myidc++){ // do each cpu left/right
				  for(i=1;i<=n[1][myidc]*pumpsize;i++){ // do each line left/right
				    ch=fgetc(incombine[myidc]);
				    fputc(ch,out);
				  }
				}
			      }
			    }
			  }
			}
		      }
		      else if(0){ // no longer used, fixed image.c in grmhd to be like my code
			// stupid gammie inverted images
			// pipe the columns of this row
			for(k=1;k<=ntile2;k++){ // do each tile up/down (assumes all cpus are same size)
			  for(kk=ppncpux2-1;kk>=0;kk--){ // do each cpu up/down (assumes all cpus are same size)
			    for(jj=1;jj<=n[2][0];jj++){ // each line up/down
			      for(j=1;j<=ntile1;j++){ // do each tile left/right
				for(myidc=kk*ppncpux1;myidc<ppncpux1+kk*ppncpux1;myidc++){ // do each cpu left/right
				  for(i=1;i<=n[1][myidc]*pumpsize;i++){ // do each line left/right
				    ch=fgetc(incombine[myidc]);
				    fputc(ch,out);
				  }
				}
			      }
			    }
			  }
			}
		      }
		    
		      // now have completed transfer for that row of cpus
		      
		      // close files
		      for(myidc=0;myidc<ppnumprocs;myidc++){
			pclose(incombine[myidc]); // close cpu myidc file
		      }
		      pclose(out); // close output file
		      // now remove old single cpu files
		      if(REMOVECPUIMAGES){
			printf("%s\n",command2);
			system(command2);
		      }	
		    }
		    else{
		      printf("skipped because couldn't find files to make %s\n",outname);
		      fflush(stdout);
		    }
		    if(prim==1) printf("s"); // indicates done an image
		    else if(prim==2) printf("v"); // indicates done an image
		    fflush(stdout);
		  }// loop over slices
		}// loop over components if a vector
	      }// loop over comp types
	    }// loop over scalars/vectors
	  }// loop over view types
	}// loop over scalar and vector
	firsthit=0;
	printf(" done\n");
	fflush(stdout);
      }// loop over images
      printf("done on all combines of images\n");
      fflush(stdout);
    }// end if doing new image combine



    
    if(CDOPAR&&(tocont==0)){
      
      generalcombine(1);
      generalcombine(2);
      
      if(FULLINPUT!=2){
	generalcombine(3);
	generalcombine(4);
      }
    }
    
    // Primitive dumps
    if(CDOPDUMP){
      generalcombine(5); 
    }

  
    if(CDODUMP){
      generalcombine(6);
    }


    if(0&&CDONPDUMP){ // no longer output, just postprocess on already combined files and compute npdumps
      generalcombine(7);
    }



    if(CDOADUMP){
      generalcombine(8);
    }
    


    if(CDOFLOORDUMP){
      generalcombine(9);
    }
    
    
    if(CDOAVG2D){
      generalcombine(10);
    }    

    if(CDOFLINEDUMP){
      generalcombine(11);
    }

  }// end of combine routines

  //////////////////////////
  //////////////////////////
  //////////////////////////
  ///////////////////////
  //////////////////////////
  /////////////////
  //
  //   BEGIN INTERP, EXPAND, SPLIT, and npcompute routines
  //
  //
  //
  ///////////////////////
  //////////////////////////
  /////////////////



  if(DOINTERP||DOEXPAND||(DOCOMBINE&&CDOCALCNPDUMP)||(DOCOMBINE&&CDOCALCDUMP)||(DOCOMBINE&&CDOCALCAVG2D)||(DOCOMBINE&&CDOCALCFLINEDUMP) ){
    printf("doing init interp or expand\n\n");
    fprintf(stderr,"proc: %s Start Initialization.\n",myidtxt) ;
    fflush(stdout);

    if((error=init(argc,argv,envp))>0){
      fprintf(fail_file,"\nError: initialize_par: %d\n",error);
      myexit(1);
    }
    printf("done\n");
    fflush(stdout);
  }


  ////////////////////
  //
  //   NP COMPUTE
  //  
  if(DOCOMBINE&&CDOCALCNPDUMP){
    printf("doing npdump compute\n");
    fflush(stdout);
    
    for(i=numnpdumpsppc;i<numnpdumps;i++){
      if(ONENPDUMPONLY>-1){
	i=ONENPDUMPONLY;
      }
      else if(ONENPDUMPONLY==-2){
	break;
      }
      printf("doing npcompute from dump#: %d  ... ",i);
      fflush(stdout);
      
      // read in data into single array for all variables
      printf("reading ... ");
      fflush(stdout);
      init_rundat(i,DTYPE);
      printf("done ... ");
      fflush(stdout);
      
      printf("dumping npdump ... ");
      fflush(stdout);
      dump(NULL,i,NPCOMPUTETYPE,-1);// not NPTYPE since need to tell to compute sigma
      printf("done\n");
      fflush(stdout);
      
      if(ONENPDUMPONLY!=-1){
	break;
      }
    }
  }

  ////////////////////
  //
  //   CALC COMPUTE
  //  
  if(DOCOMBINE&&CDOCALCDUMP){
    printf("doing calc dump compute\n");
    fflush(stdout);
    
    for(i=numcalcdumpsppc;i<numcalcdumps;i++){
      if(ONECALCDUMPONLY>-1){
	i=ONECALCDUMPONLY;
      }
      else if(ONECALCDUMPONLY==-2){
	break;
      }
      printf("doing calccompute from dump#: %d  ... ",i);
      fflush(stdout);
      
      // read in data into single array for all variables
      printf("reading ... ");
      fflush(stdout);
      init_rundat(i,DTYPE);
      printf("done ... ");
      fflush(stdout);
      
      printf("dumping calcdump ... ");
      fflush(stdout);
      dump(NULL,i,CALCTYPE,-1);
      printf("done\n");
      fflush(stdout);
      
      if(ONECALCDUMPONLY!=-1){
	break;
      }
    }
  }


  ////////////////////
  //
  //   FLINE COMPUTE
  //  
  if(DOCOMBINE&&CDOCALCFLINEDUMP){
    printf("doing fieldline dump compute\n");
    fflush(stdout);
    
    for(i=numfldumpsppc;i<numfldumps;i++){
      if(ONEFLINEDUMPONLY>-1){
	i=ONEFLINEDUMPONLY;
      }
      else if(ONEFLINEDUMPONLY==-2){
	break;
      }
      printf("doing fieldlinecompute from dump#: %d  ... ",i);
      fflush(stdout);
      
      // read in data into single array for all variables
      printf("reading ... ");
      fflush(stdout);
      init_rundat(i,DTYPE);
      printf("done ... ");
      fflush(stdout);
      
      printf("dumping fieldlinedump ... ");
      fflush(stdout);
      dump(NULL,i,FLINETYPE,-1);
      printf("done\n");
      fflush(stdout);
      
      if(ONEFLINEDUMPONLY!=-1){
	break;
      }
    }
  }

  ////////////////////
  //
  //   CALC COMPUTE from avg2d ( could add fieldline calc from avg2d)
  //  
  if(DOCOMBINE&&CDOCALCAVG2D){
    printf("doing calc dump from avg2d compute ... ");
    fflush(stdout);
    
    // read in data into single array for all variables
    printf("reading ... ");
    fflush(stdout);
    init_rundat(-1,AVG2DTYPE);
    printf("done ... ");
    fflush(stdout);
      
    printf("dumping calcdump ... ");
    fflush(stdout);
    
    dump(NULL,-1,CALCTYPE,-1);
    printf("done\n");
    fflush(stdout);
    
  }
  


  if(DOINTERP){
    if(IDOPAR||IDOPDUMP||IDODUMP||IDONPDUMP||IDOCALCDUMP||IDOCALCAVG2D||IDOFLOORDUMP||IDOADUMP||IDOAVG2D||IDOFLINEDUMP){
      printf("doing init interp\n\n");
      fflush(stdout);
      
      printf("doing par interp(creation)  ...  ");
      fflush(stdout);
      // 0 and 1 should already be there since just read it in or just made it
      globalinterpmod=1;
      error+=init_dx(DUMN1PP,DUMN2PP,1,0,0,0,startix,1,IDUMPOUTTYPE);
      error+=init_x(DUMN1PP,DUMN2PP,1,0,0,0,startix,1,IDUMPOUTTYPE);
      printf("done\n ");
      fflush(stdout);
      if(IDOPAR&&(tocont==0)){
	printf("doing par interp(output)  ...  ");
	fflush(stdout);
	
	init_outgparm(2);
	init_outgparm(3);
	
	printf("done\n");
	fflush(stdout);
      }
    }
    call_code=0;    


    // never want to interp the pdumps

    // figure out what loop needs to be done
    for(i=0;i<numdumps;i++){
      idodumplist[i]=0;
    }
    if(ONEDUMPONLY>-1) idodumplist[ONEDUMPONLY]=1;
    if(ONEFLINEDUMPONLY>-1) idodumplist[ONEFLINEDUMPONLY]=1;
    if(ONECALCDUMPONLY>-1) idodumplist[ONECALCDUMPONLY]=1;
    
    if((IDODUMP&&(ONEDUMPONLY==-1))||(IDOFLINEDUMP&&(ONEFLINEDUMPONLY==-1))||(IDOCALCDUMP&&(ONECALCDUMPONLY==-1)) ){
      for(i=numdumpsppi;i<numdumps;i++){
	idodumplist[i]=1;
      }
    }
    
    // dumps
    if(IDODUMP||IDOCALCDUMP||IDOFLINEDUMP){
      printf("doing dump interps\n");
      fflush(stdout);
      
      fprintf(stderr,"%d %d\n",numdumpsppi,numdumps);
      for(i=numdumpsppi;i<numdumps;i++){ 
	if(idodumplist[i]==0) continue; // else do it
	printf("doing interp dump#: %d  ... ",i);
	fflush(stdout);
	
	// read in data into single array for all variables
	printf("reading ... ");
	fflush(stdout);
	init_rundat(i,DTYPE);
	printf("done ... ");
	fflush(stdout);

	if(IDODUMP){
	  printf("dumping ... ");
	  fflush(stdout);
	  dump(NULL,i,DTYPE,IDUMPOUTTYPE);
	  printf("done\n");
	  fflush(stdout);
	}
	if(IDOCALCDUMP){
	  printf("doing calc dump from interp dump compute: %d ... ",i);
	  fflush(stdout);
	  
	  dump(NULL,i,CALCTYPE,IDUMPOUTTYPE);
	  printf("done\n");
	  fflush(stdout);
	}
	if(IDOFLINEDUMP){
	  printf("doing fieldline dump from interp dump compute: %d ... ",i);
	  fflush(stdout);
	  
	  dump(NULL,i,FLINETYPE,IDUMPOUTTYPE);
	  printf("done\n");
	  fflush(stdout);
	}

	if(ONEDUMPONLY!=-1){
	  break;
	}
      }
    }
    if(IDOFLOORDUMP){
      printf("doing floordump interps\n");
      
      for(i=numfloordumpsppi;i<numfloordumps;i++){
	if(ONEFLOORDUMPONLY!=-1){
	  i=ONEFLOORDUMPONLY;
	}
	else if(ONEFLOORDUMPONLY==-2){
	  break;
	}
	printf("doing interp floordump#: %d  ...  ",i);
	fflush(stdout);
	
	init_rundat(i,FLTYPE); // to get floor data
	dump(NULL,i,FLTYPE,IDUMPOUTTYPE);		
	printf("done\n");
	fflush(stdout);
	if(ONEFLOORDUMPONLY!=-1){
	  break;
	}
      }
    }
    if(IDONPDUMP){
      printf("doing npdump interps\n");
      fflush(stdout);

      for(i=numnpdumpsppi;i<numnpdumps;i++){
	if(ONENPDUMPONLY!=-1){
	  i=ONENPDUMPONLY;
	}
	else if(ONENPDUMPONLY==-2){
	  break;
	}
	printf("doing interp npdump#: %d  ...  ",i);
	fflush(stdout);
	
	init_rundat(i,NPTYPE);// not 3 since below, and to keep same
	dump(NULL,i,NPTYPE,IDUMPOUTTYPE);
	printf("done\n");
	fflush(stdout);

	if(ONENPDUMPONLY!=-1){
	  break;
	}
      }
    }
    if(IDOADUMP){
      printf("doing adump interps\n");
      fflush(stdout);
      
      for(i=numadumpsppi;i<numadumps;i++){
	if(ONEADUMPONLY!=-1){
	  i=ONEADUMPONLY;
	}
	else if(ONEADUMPONLY==-2){
	  break;
	}
	printf("doing interp adump#: %d  ...  ",i);
	fflush(stdout);
	

	init_rundat(i,ADTYPE); 	
	dump(NULL,0,ADTYPE,IDUMPOUTTYPE);
	printf("done\n");
	fflush(stdout);

	if(ONEADUMPONLY!=-1){
	  break;
	}
      }
    }
    if(IDOAVG2D){
      printf("doing interp avg2d dump  ...  ");
      fflush(stdout);
      // only 1 avg2d file
      init_rundat(-1,AVG2DTYPE);
      dump(NULL,-1,AVG2DTYPE,IDUMPOUTTYPE);
      printf("done\n");
      fflush(stdout);
      ////////////////////
      //
      //   CALC COMPUTE from avg2d ( could add fieldline calc for iavg2d)
      //  
      if(IDOCALCAVG2D){ // only done if interp avg2d too.  Could also do interp calctype directly from raw avg1d, which is essentially what's going on here, asuming interp stuff inited beforehand.
	printf("doing calc dump from interp avg2d compute ... ");
	fflush(stdout);
	
	dump(NULL,-1,CALCTYPE,IDUMPOUTTYPE);
	printf("done\n");
	fflush(stdout);
      }
    }    

    if(IDOIMAGE){// can be used to automagically go from r8(or r8.gz) to ppm.gz
      printf("doing par interp for images n1=%d n2=%d  ...  ",IMGN1PP,IMGN2PP);
      fflush(stdout);
      // 0 and 1 should already be there since just read it in or just made it
      globalinterpmod=1;
      error+=init_dx(IMGN1PP,IMGN2PP,1,0,0,0,startix,1,IIMGOUTTYPE);
      error+=init_x(IMGN1PP,IMGN2PP,1,0,0,0,startix,1,IIMGOUTTYPE);
      if(IDOPARI&&(tocont==0)){
	printf("doing par interp(output)  ...  ");
	fflush(stdout);
	
	init_outgparm(1);
	init_outgparm(2);
	init_outgparm(3);
	
	printf("done\n");
	fflush(stdout);
      }
      printf("done\n ");
      fflush(stdout);
    
      call_code=0; // as if first time always(not really used for image())
      printf("Begin with interp of %d images: STARTI: %d ENDI: %d\n",numimages,STARTIi,ENDIi);
      // loop over images
      for(i=STARTIi;i<ENDIi;i++){
	// read in data into single array for all variables

	printf("begin read image#: %d ",i);
	fflush(stdout);
	// currently only image the primitive variable or simple computations from them.
	init_runimage(i,-1,-1,call_code,0);
	//t=tstart; // assign time as would be done in init_params
	t = tstart+(FTYPE)(i)*DTi ;

	printf(" ... done read image#: %d t=%15.10g.  Doing interp ... ",i,t);	
	fflush(stdout);
	// interpolate to new data sets(all variables) and output data
	// dump type same since uses same init_dx/x data
	image(i,-1,-1,call_code,IIMGOUTTYPE); // set which images to produce here
	
	printf(" ... done interp image#: %d\n",i);
	fflush(stdout);
      }// done over all images
      printf("Done with all interp of %d images\n",numimages);
      // can also just use interp but don't have interp to get r8->ppm for all files
      if(DOCONVERTR8){
	loadpal();
	
	// some how loop over files in dir trees
	// convertr8(IMGN1PP,IMGN2PP,GZIPIMAGE,GZIPIMAGEINPUT,pal,inname,outname);
      }
    
    }// done if doing image pp


  }// done if doing interp



  // only setup for normal active grid input/output
  // expand first so only have to use combined data when expanding
  // cases:
  // 1) start with multiple cpu data
  //    a) use combine of postprocess
  //    b) expand
  //    c) split if necessary
  // 2) start with single cpu data
  //    a) expand
  //    b) split

  // for expand: N1PP,N2PP,N3PP in global.h are input file grid sizes and IMGN1PP,IMGN2PP,IMGN3PP are output grid sizes

  if(DOEXPAND==1){
    call_code=0;

    printf("doing par interp  ...  ");
    fflush(stdout);
    // 0 and 1 should already be there since just read it in
    outtype=0; // force, unlike with above interp
    globalinterpmod=1;
    if(FULLOUTPUT==2){
      error+=init_dx(DUMN1PP,DUMN2PP,DUMN3PP,N1PPBND,N2PPBND,N3PPBND,startix,1,outtype);
      error+=init_x(DUMN1PP,DUMN2PP,DUMN3PP,N1PPBND,N2PPBND,N3PPBND,startix,1,outtype);
    }
    else     if(FULLOUTPUT==2){
      // assumes 2 boundary zones
      error+=init_dx(DUMN1PP,DUMN2PP,1,1,1,1,startix,1,outtype);
      error+=init_x(DUMN1PP,DUMN2PP,1,1,1,1,startix,1,outtype);
    }
    else     if(FULLOUTPUT==0){
      error+=init_dx(DUMN1PP,DUMN2PP,1,0,0,0,startix,1,outtype);
      error+=init_x(DUMN1PP,DUMN2PP,1,0,0,0,startix,1,outtype);
    }
    // no output needed here
    printf("done\n");

    // below nearly identical to idodump, but forced to only be normal grid, no zoom, etc.
    printf("doing dump expand ... ");
    fflush(stdout);
      
    // read in data into single array for all variables
    sprintf(temps,"cp %s %s\n",FILETOEXPAND,"toedump0000.dat");
    system(temps);
    init_rundat(0,EXPANDTYPE); // special for toedumpxxxx.dat file readin
    //    bound(NULL,NULL,-1,-1,-1);
    if(DOSPLIT==1){
      strcpy(filename,FILETOSPLIT);
    }
    else sprintf(filename,"1_pdump.in");
    if(!(file_expand=fopen(filename,"wt"))){
      fprintf(fail_file,"couldn't open %s for expansion\n",filename);
      myexit(1);
    }
    dump(file_expand,0,PDTYPE,0);
    
    printf("done\n");
    fflush(stdout);
  }


  // only setup for normal active grid input/output
  // for split: n[][] are output grid sizes, input sizes must match ppncpux1*ppncpux2 and input file grid size(which will be related to Nx if no expand or IMGNx if expand)

  if(DOSPLIT==1){

    if(DOEXPAND==1){

      if(ppncpux1*ppncpux2==4){
	n[1][0]=DUMN1PP;
	n[1][1]=DUMN1PP;
	n[1][2]=DUMN1PP;
	n[1][3]=DUMN1PP;

	// don't have to do equal splitting
	n[2][0]=DUMN2PP/ppncpux1*ppncpux2;
	n[2][1]=DUMN2PP/ppncpux1*ppncpux2;
	n[2][2]=DUMN2PP/ppncpux1*ppncpux2;
	n[2][3]=DUMN2PP/ppncpux1*ppncpux2;
	
	n[3][0]=DUMN3PP;
	n[3][1]=DUMN3PP;
	n[3][2]=DUMN3PP;
	n[3][3]=DUMN3PP;
      }
      if(ppncpux1*ppncpux2==2){
	n[1][0]=DUMN1PP;
	n[1][1]=DUMN1PP;
	
	n[2][0]=DUMN2PP/ppncpux1*ppncpux2;
	n[2][1]=DUMN2PP/ppncpux1*ppncpux2;
	
	n[3][0]=DUMN3PP;
	n[3][1]=DUMN3PP;
      }

    }
    else{
      if(ppncpux1*ppncpux2==4){
	n[1][0]=N1PP;
	n[1][1]=N1PP;
	n[1][2]=N1PP;
	n[1][3]=N1PP;
	
	n[2][0]=N2PP/ppncpux1*ppncpux2;
	n[2][1]=N2PP/ppncpux1*ppncpux2;
	n[2][2]=N2PP/ppncpux1*ppncpux2;
	n[2][3]=N2PP/ppncpux1*ppncpux2;
	
	n[3][0]=N3PP;
	n[3][1]=N3PP;
	n[3][2]=N3PP;
	n[3][3]=N3PP;
      }
      if(ppncpux1*ppncpux2==2){
	n[1][0]=N1PP;
	n[1][1]=N1PP;
	
	n[2][0]=N2PP/ppncpux1*ppncpux2;
	n[2][1]=N2PP/ppncpux1*ppncpux2;
	
	n[3][0]=N3PP;
	n[3][1]=N3PP;
      }

    }
    // 1 to many only, use init_rundat and dump of interp to get to different size
    
    printf("splitting ... ");

    strcpy(filename,FILETOSPLIT);
	
    if( (in=fopen(filename,"rt"))==NULL){
      printf("can't open input file: %s\n",filename);
      myexit(1);
    }
    
    for(i=0;i<ppnumprocs;i++){
      sprintf(myidctxt,CPUTXT,i);
      sprintf(filename,"1_pdump.in%s",myidctxt); // 1 so don't overwrite any such files already in dir from old run perhaps, just rename to 0_* when moved to run dir
      
      if( (out=fopen(filename,"wt"))==NULL){
	printf("can't open output file: %s\n",filename);
	myexit(1);
      }
      // assume input and output file will be new type with version and not full grid
      newskip=7;
      count=0;
      length = newskip;
      
      // store position of input file, rewind, dump initial stuff, then restore position

      fgetpos(in,&fpos);
      fseek(in,0,SEEK_SET);
      go=1;      
      while((!feof(in))&&(go==1)){
	ch=fgetc(in);
	if(ch=='\n') count++;
	if(!feof(in)) fputc(ch,out);
	
	/*	      
		      if(count==1){
		      fscanf(in,"%f %d %d",&tempfloat,&dumi[0],&dumi[1]);
		      while(fgetc(in)!='\n');
		      
		      fprintf(out,"%21.10g %6d %6d\n",tempfloat,0,dumi[1]);
		      count++;
		      }
	*/
	if(count==newskip-2){
	  fscanf(in,"%d %d %d",&dumi[1],&dumi[2],&dumi[3]);
	  while(fgetc(in)!='\n');
	  
	  fprintf(out," %6d %6d %6d\n",n[1][i],n[2][i],n[3][i]);
	  count++;
	}
	if(count==length) go=0;
      }      
      if(i>0) fsetpos(in,&fpos); // reset position to where it was

      // now copy info over
      go=1;
      length = newskip+(n[1][i])*(n[2][i]);
      while((!feof(in))&&(go==1)){
	ch=fgetc(in);
	if(ch=='\n') count++;
	if(!feof(in)) fputc(ch,out);
	
	if(i<ppnumprocs-1){
	  if(count==length) go=0;
	}
      }
      fclose(out);
    }
    fclose(in);
    printf("done\n");
    fflush(stdout);
  } // done if split

  myexit(0);
  return(0);
}



void loadpal(void)
{
  FILE * in;
  char ifnam[MAXFILENAME];
  char temps[MAXFILENAME];

  sprintf(temps,"%s%s",DATADIR,IMAGEDIR);	    
  sprintf(ifnam,"%sgenimages%s",temps,".pal");
  if( (in = fopen(ifnam,"rb"))==NULL){
    fprintf(fail_file,"error opening %s\n",ifnam) ;
    myexit(2) ;
  }
  for(j=0;j<3;j++){
    for(i=0;i<256;i++){
      fread(&pal[j][i],sizeof(short),1,in); 
    }
  }
  fclose(in);

}


void convertr8(int nx,int ny,int gzipout,int gzipin,short (*pal)[256],char* inname,char* outname)
{
  FILE * in;
  FILE * out;
  char ifnam[MAXFILENAME];
  char temps[MAXFILENAME];

  if(GZIPIMAGEINPUT==0){
    in = fopen(inname,"rb");
  }
  else if(GZIPIMAGEINPUT>0){
    sprintf(temps,"gzip -d %s >",inname);
    strcpy(ifnam,temps); // for below fprintf
    in = popen(ifnam,"r");
  }
  if(in==NULL){
    fprintf(fail_file,"error opening image file: %s\n",ifnam) ;
    myexit(2) ;
  } 
      // need ifnam, filename for input file
      // need outfile file(out)


  if(GZIPIMAGE==0) out = fopen(outname,"wb");
  else{
    sprintf(temps,"gzip > %s.gz",outname);
    strcpy(ifnam,temps); // for below fprintf
    out = popen(ifnam,"w");
  }
  if(out==NULL){
    fprintf(fail_file,"error opening image file: %s\n",ifnam) ;
    myexit(2) ;
  }
  if(IMAGEFORMAT==0){
    r82ppm(nx,ny,pal,in,out);
  }
  else r82ras(nx,ny,pal,in,out);


}

void r82ppm(int nx,int ny,short (*pal)[256],FILE*in,FILE*out)
{
  short liq;

  fprintf(out,"P6\n");
  while(fgetc(in)!='\n'); // skip to next line
  // just output next 3 lines from input
  i=0;
  while(i<3){
    ch=fgetc(in);
    if(ch=='\n') i++;
    fputc(ch,out);
  }




  LOOPI{
    liq=fgetc(in);
    fputc(pal[0][(short)liq] , out);   /* write red */
    fputc(pal[1][(short)liq] , out);   /* write green */
    fputc(pal[2][(short)liq] , out);   /* write blue */
  }
}







void r82ras(int nx,int ny,short (*pal)[256],FILE*in,FILE*out)
{
  int nmap ;
  char dum ;
  struct rasterfile header ;
  


  // read raw header of input image(just skip it)
  while(fgetc(in)!='\n'); // skip to next line 
  while(fgetc(in)!='\n'); // skip to next line 
  while(fgetc(in)!='\n'); // skip to next line 
  while(fgetc(in)!='\n'); // skip to next line

  /* set up header */
  header.ras_magic = RAS_MAGIC ;
  header.ras_width = nx ;
  header.ras_height = ny ;
  header.ras_depth = 8 ;
  header.ras_length = nx*ny ;
  header.ras_type = RT_STANDARD ;
  header.ras_maptype = RMT_EQUAL_RGB ;
  header.ras_maplength = 3*256 ;
  
  /* swap & write it */
  decint(&(header.ras_magic),8) ;
  fwrite(&(header.ras_magic), sizeof(int), 8, out) ;
  
  /* write the palette file */
  nmap = 0 ;
  for(i=0;i<3;i++){
    for(j=0;j<256;j++){
      fwrite(&pal[i][j], sizeof(char), 1, out) ;
      nmap++ ;
    }
  }
  decint(&(header.ras_maplength),1) ;
  if(nmap != header.ras_maplength) {
    fprintf(stderr,"err: header: %d, nmap: %d\n",
	    header.ras_maplength, nmap) ;
  }
  
  /* now read & write the image */
  while(fread(&dum, sizeof(char), 1, in) == 1) {
    fwrite(&dum, sizeof(char), 1, out) ;
  }
  
}

void decint(int *lp, int n)
{
  unsigned int t;
  static unsigned long lmask = 0x00ff0000, rmask = 0x0000ff00;
  
  for(; n--; lp++) {
    t = *lp;
    *lp = (t >> 24) | (t << 24) | ((t << 8) & lmask) |
      ((t >> 8) & rmask);
  }
}

// returns number of comment lines that exist in a given file type
// check whether file is new or old format
// which:
// type for : adump,pdump,dump,floordump,npdump,(same as far as first check is concerned), etc.
#define NUMPOUNDS 10
int checktype(FILE * in, FILE * out,int which, int fixout)
{
  int shouldfix=0;
  char ch;
  int toreturn=0;
  int pound=0,hitreturn=0;
  // basically go 10 lines and check for any #'s.  Really only need to go to <10, variance among files(this won't work for VERY small grids(<=2x2 perhaps)

  while(1){
    ch=fgetc(in);
    if(ch=='#'){ pound++; }
    if(ch=='\n'){ hitreturn++; }
    if(hitreturn==10) break;
  }
  switch(which){
  case 2: // global or active grid
    if(pound==2){ toreturn=3; }
    else if(pound==1){ toreturn=1; shouldfix=1;}
    break;
  case 3:
  case 4:
  case 5:
    if(pound==4){ toreturn=7; }
    else if(pound==3){ toreturn=5; shouldfix=1;}
    else if(pound==0){ toreturn=1;} // GAMMIE
    break;
  case 6:
    if(pound==6){ toreturn=9; }
    else if(pound==5){ toreturn=7; shouldfix=1;}
    break;
  default:
    fprintf(fail_file,"which: %d not setup\n",which);
    myexit(1);
    break;
  }
  if(shouldfix&&fixout){
    switch(which){
    case 2:
      fprintf(out,"#GRIDVER GRIDTYPE\n%d %d\n",GRIDVER,GRIDTYPE);
      break;
    case 3:
      fprintf(out,"#DVER DTYPE\n%d %d\n",DVER,DTYPE);
      break;
    case 4:
      fprintf(out,"#FLVER FLTYPE\n%d %d\n",FLVER,FLTYPE);
      break;
    case 5:
      fprintf(out,"#NPVER NPTYPE\n%d %d\n",NPVER,NPTYPE);
      break;
    case 6:
      fprintf(out,"#AVG2DVER AVG2DTYPE\n%d %d\n",AVG2DVER,AVG2DTYPE);
      break;
    }
  }
  rewind(in);
  return(toreturn);
}



void generalcombine(int which)
{
  // general combine for everything except images which has more loops

  
  int startloopdata1,endloopdata1,startloopdata2,endloopdata2;
  int qtip,i,j,k,ii,jj,kk;
  char whichname[100],outname[100],inheader[100];
  FILE *incombine[MAXCPUS],*out,*testp;
  char fullname[MAXCPUS][200];
  int skipnum;
  int format,startloop,endloop;
  int skiptype;
  int newskip[MAXCPUS],length[MAXCPUS],count[MAXCPUS];
  int oneonly;
  int myidc;
  char myidctxt[100];
  
  switch(which){
  case 1:
    strcpy(whichname,"0_grid1.par");
    strcpy(outname,whichname);
    strcpy(inheader,outname);
    skipnum=2;
    format=0; // sprintf(filename,"0_gridact2.par.%02d",i);  
    startloop=0;
    endloop=1; // <endloop version of loop
    skiptype=1; // whether fulloutput requries more complicated patchwork
    oneonly=-1;
    break;
  case 2:
    strcpy(whichname,"0_grid2.par");
    strcpy(outname,whichname);
    strcpy(inheader,outname);
    skipnum=2;
    format=0;
    startloop=0;
    endloop=1; // <endloop version of loop
    skiptype=1;
    oneonly=-1;
    break;
  case 3:
    strcpy(whichname,"0_gridact1.par");
    strcpy(outname,whichname);
    strcpy(inheader,outname);
    skipnum=2;
    format=0;
    startloop=0;
    endloop=1; // <endloop version of loop
    skiptype=0;
    oneonly=-1;
    break;
  case 4:
    strcpy(whichname,"0_gridact2.par");
    strcpy(outname,whichname);
    strcpy(inheader,outname);
    skipnum=2;
    format=0; 
    startloop=0;
    endloop=1;
    skiptype=0;
    oneonly=-1;
    break;
  case 5:
    strcpy(whichname,"pdump");
    strcpy(outname,whichname);
    strcpy(inheader,outname);
    skipnum=3;
    format=1; // sprintf(filename,"pdump%04d%s.%02d",k,".dat",i);
    startloop=numpdumpsppc;
    endloop=numpdumps; // <endloop version of loop
    skiptype=1; // need to deal with fulloutput capability
    oneonly=ONEPDUMPONLY;
    break;
  case 6:
    strcpy(whichname,"dump");
    strcpy(outname,whichname);
    strcpy(inheader,outname);
    skipnum=3;
    if(!GAMMIEIMAGE){
     format=1; // sprintf(filename,"pdump%04d%s.%02d",k,".dat",i);
    }
    else{
      format=10; // GAMMIE
    }
    startloop=numdumpsppc;
    endloop=numdumps; // <endloop version of loop
    skiptype=1; // need to deal with fulloutput capability
    oneonly=ONEDUMPONLY;
    break;
  case 7:
    strcpy(whichname,"npdump");
    strcpy(outname,whichname);
    strcpy(inheader,outname);
    skipnum=3;
    format=1; // sprintf(filename,"pdump%04d%s.%02d",k,".dat",i);
    startloop=numnpdumpsppc;
    endloop=numnpdumps; // <endloop version of loop
    skiptype=1; // need to deal with fulloutput capability
    oneonly=ONENPDUMPONLY;
    break;
  case 8:
    strcpy(whichname,"adump");
    strcpy(outname,whichname);
    strcpy(inheader,outname);
    skipnum=3;
    format=1; // sprintf(filename,"pdump%04d%s.%02d",k,".dat",i);
    startloop=numadumpsppc;
    endloop=numadumps; // <endloop version of loop
    skiptype=1; // need to deal with fulloutput capability
    oneonly=ONEADUMPONLY;
    break;
  case 9:
    strcpy(whichname,"floor");
    strcpy(outname,whichname);
    strcpy(inheader,outname);
    skipnum=4;
    format=1; // sprintf(filename,"pdump%04d%s.%02d",k,".dat",i);
    startloop=numfloordumpsppc;
    endloop=numfloordumps; // <endloop version of loop
    skiptype=1; // need to deal with fulloutput capability
    oneonly=ONEFLOORDUMPONLY;
    break;
  case 10:
    strcpy(whichname,"0_avg2d.dat");
    strcpy(outname,whichname);
    strcpy(inheader,outname);
    skipnum=6;
    format=0;
    startloop=0;
    endloop=1; // <endloop version of loop
    skiptype=1; // need to deal with fulloutput capability
    oneonly=-1;
    break;
  case 11:
    strcpy(whichname,"fldump");
    strcpy(outname,whichname);
    strcpy(inheader,outname);
    skipnum=3;
    format=1; // sprintf(filename,"pdump%04d%s.%02d",k,".dat",i);
    startloop=numfldumpsppc;
    endloop=numfldumps; // <endloop version of loop
    skiptype=1; // need to deal with fulloutput capability
    oneonly=ONEFLINEDUMPONLY;
    break;
  default:
    break;
  }



  // must also account for different writes in loop itself with switch

      
  printf("Combining %s\n",whichname);
  fflush(stdout);
  
  sprintf(dirname,"");	
  
  for(qtip=startloop;qtip<endloop;qtip++){
    if(oneonly>-1){
      qtip=oneonly;
    }
    else if(oneonly==-2){
      break;
    }
                      
    printf("combining %s# %d format: %d ",whichname,qtip,format);
    fflush(stdout);
    

    // check to see if all files exist
    skip=0;
    for(myidc=0;myidc<ppnumprocs;myidc++){
      sprintf(myidctxt,CPUTXT,myidc);
      if(format==1){
	sprintf(fullname[myidc],"%s%04d%s%s",inheader,qtip,DATEXT,myidctxt);
      }
      if(format==10){ // GAMMIE
	sprintf(fullname[myidc],"%s%s%03d%s",DUMPDIR,inheader,qtip,myidctxt);
      }
      else if(format==0){
	sprintf(fullname[myidc],"%s%s",inheader,myidctxt);
      }
      // test for existence of multiple cpu data
      testfp=fopen(fullname[myidc],"rb");
      if(testfp==NULL){
	skip++;
	fprintf(stdout,"couldn't find: %s\n",fullname[myidc]); fflush(stdout);
      }
      else fclose(testfp);
    }

    // now you know whether all files exist, loop over them and append them
    if(skip==0){
      // create output file
      //sprintf(command,"gzip > %s",outname);
      //out = popen(command,"w");
      // open output file
      if(format==1){
	sprintf(outname,"%s%s%04d%s",DUMPDIR,inheader,qtip,".dat");
      }
      if(format==10){
	sprintf(outname,"%s%s%03d",DUMPDIR,inheader,qtip);
      }
      else if(format==0){
	sprintf(outname,"%s",inheader);
      }
      
      out=fopen(outname,"wt");
      if(out==NULL){
	fprintf(fail_file,"error opening image file: %s\n",command) ;
	myexit(1) ;
      }
      fprintf(stdout,"output full filename: %s\n",outname);
      // loop over cpus, hopping down rows at a time
      for(myidc=0;myidc<ppnumprocs;myidc+=ppncpux1){
	
	for(jj=0;jj<ppncpux1;jj++){ // open input files for this row of cpus(each column of the row)
	  //sprintf(command,"gzip -d < %s",fullname[myidc+jj]);
	  //incombine[myidc+jj] = popen(command,"r");
	  fprintf(stdout,"doing cpu#: i=%d j=%d\n",myidc,jj);
	  incombine[myidc+jj]=fopen(fullname[myidc+jj],"rt");
	  if(incombine[myidc+jj]==NULL){
	    fprintf(fail_file,"error opening image file: %s\n",command) ;
	    myexit(1) ;
	  }
	  // find number of comment lines(all but raw data)
	  newskip[myidc+jj]=checktype(incombine[myidc+jj],out,skipnum,!(myidc+jj));// only fix if first cpu so writting fresh out
	  count[myidc+jj]=0;
	}

	// now that all cpus on this row are open, do the main loop
	for(jj=0;jj<ppncpux1;jj++){

	  // pipe and fix header from myidc==0
	  if((myidc==0)&&(jj==0)){
	    go=1;      
	    while((!feof(incombine[myidc+jj]))&&(go==1)){
	      
	      switch(which){
	      case 1:
	      case 2:
	      case 3:
	      case 4:
		break;
	      case 5:
	      case 6:
	      case 7:
	      case 8:
	      case 9:
	      case 11:
		if(!GAMMIEIMAGE){
		  if(count[myidc+jj]==newskip[myidc+jj]-2){
		    fscanf(incombine[myidc+jj],"%d %d %d",&dumi[1],&dumi[2],&dumi[3]);
		    while(fgetc(incombine[myidc+jj])!='\n');
		    count[myidc+jj]++;
		    if(FULLINPUT==0){
		      fprintf(out," %6d %6d %6d\n",N1PP,N2PP,N3PP);
		    }
		    else fprintf(out," %6d %6d %6d\n",N1PPM,N2PPM,N3PPM);
		  }
		}
		else{
		  if(count[myidc+jj]==0){
		    fscanf(incombine[myidc+jj],"%d %d %lf %lf %lf",&dumi[1],&dumi[2],&dumf[1],&dumf[2],&dumf[3]);
		    while(fgetc(incombine[myidc+jj])!='\n');
		    count[myidc+jj]++;
		    fprintf(out," %6d %6d %15.10g %15.10g %15.10g\n",N1PP,N2PP,dumf[1],dumf[2],dumf[3]);
		  }
		}
		break;
	      case 10:
		if(count[myidc+jj]==newskip[myidc+jj]-4){
		  fscanf(incombine[myidc+jj],"%d %d %d",&dumi[1],&dumi[2],&dumi[3]);
		  while(fgetc(incombine[myidc+jj])!='\n');
		  if(FULLINPUT==0){
		    fprintf(out," %6d %6d %6d\n",N1PP,N2PP,N3PP);
		  }
		  else fprintf(out," %6d %6d %6d\n",N1PPM,N2PPM,N3PPM);
		  count[myidc+jj]++;
		}
		break;
	      default:
		break;
	      }

	      ch=fgetc(incombine[myidc+jj]);
	      if(ch=='\n') count[myidc+jj]++;
	      if(!feof(incombine[myidc+jj])) fputc(ch,out);

	      if(count[myidc+jj]==newskip[myidc+jj]) go=0;

	    }
	  }
	  else{
	    // otherwise skip the header
	    for(j=1;j<=newskip[myidc+jj];j++){
	      while(fgetc(incombine[myidc+jj])!='\n');
	      count[myidc+jj]++;
	    }
	  }
	}
      
	// now write out(assumes all cpu input files on this row are just passed header)
	// pipe the columns of this row

	// figure out how far down in j should go
	if((which==1)||(which==2)||((FULLOUTPUT==2)&&(skiptype==1))){
	  if(myidc/ppncpux1==0){ // then inner j edge
	    startloopdata2=-N2PPBND;
	    endloopdata2=n[2][myidc]; // must all be same
	  }
	  else if(myidc/ppncpux1==ppncpux2-1){ // then outer j edge
	    startloopdata2=0;
	    endloopdata2=n[2][myidc]+N2PPBND; // must all be same
	  }
	  else{
	    startloopdata2=0;
	    endloopdata2=n[2][myidc]; // must all be same
	  }
	}
	else{
	  startloopdata2=0;
	  endloopdata2 = n[2][myidc]; // must all be same
	}
	
	
	for(kk=-N2PPBND;kk<N2PPBND+n[2][myidc];kk++){ // do each cpu line (assumes all cpus are same size)
	  
	  for(jj=0;jj<ppncpux1;jj++){ // over each column for this row
	    
	    // figure out how far in i should go
	    if((which==1)||(which==2)||((FULLOUTPUT==2)&&(skiptype==1))){
	      if((myidc+jj)%ppncpux1==0){ // then inner i edge
		startloopdata1=-N1PPBND;
		endloopdata1=n[1][(myidc+jj)]; // must all be same
	      }
	      else if((myidc+jj)%ppncpux1==ppncpux1-1){ // then outer i edge
		startloopdata1=0;
		endloopdata1=n[1][(myidc+jj)]+N1PPBND; // must all be same
	      }
	      else{
		startloopdata1=0;
		endloopdata1=n[1][(myidc+jj)]; // must all be same
	      }
	    }
	    else{
	      startloopdata1=0;
	      endloopdata1 = n[1][(myidc+jj)]; // must all be same
	    }
	    
	    for(ii=-N1PPBND;ii<N1PPBND+n[1][myidc+jj];ii++){ // each can have different n1 size
	      
	      // either pump data line, skip it, or ignore index
	      if( (ii>=startloopdata1)&&(ii<endloopdata1)&&(kk>=startloopdata2)&&(kk<endloopdata2) ){
		// then within correct range and can output this line
		while((ch=fgetc(incombine[myidc+jj]))!='\n'){
		  fputc(ch,out);
		}
		fputc('\n',out);
	      }
	      else{
		if((which==1)||(which==2)||((FULLOUTPUT==2)&&(skiptype==1))){// then other data really exists
		  while((ch=fgetc(incombine[myidc+jj]))!='\n'); // skip this line
		}
		///// else ignore since no data
	      }// end else if outside active data range
	    }// over each line per cpu
	  }// over a row of cpus
	}// over columns of data
	
	// close this row of cpus
	for(jj=0;jj<ppncpux1;jj++){
	  fclose(incombine[myidc+jj]);
	}
	//for(i=0;i<ppncpux1;i++){
	// pclose(incombine[myidc+i]); // close cpu myidc file
	//	}
	// flush output
	fflush(out);
      }// end loop over all cpus
      fclose(out); // only close output if output exists
    // pclose(out);
    }// end if all input files exist and can do output file
    else{
      printf("skipped because couldn't find files to make %s\n",outname);
      fflush(stdout);
    }
    if(oneonly>-1) break;

    firsthit=0;
    printf(" done\n");
    fflush(stdout);
  }// loop over dumps 


  printf("done on all combines of %s\n",whichname);
  fflush(stdout);


}
