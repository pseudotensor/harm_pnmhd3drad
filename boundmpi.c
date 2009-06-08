#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif

#include "bound.h"

// for MPI:
// in order to assign corner zones properly, must either order assignments or directly assign corner zones.  I choose to order the assignments so in the end the corner zones are right.  This eliminates any extra transfers and/or latency associated with small transfers of just th e corner zones
// true for interior zone exchanges, and periodic conditions for outer real corner zones
// I choose left-right N1M/N2M first, then up/down N1M/N2M.  Could just do N1/N2 for interior for L/R, but true boundary needs full N1M/N2M exchanged since cpu sets boundary using normal bc code which needs to get transfered to neight(i.e. currently if corner has bctype 99/? then doesn't do corner)

// for vector: ll(wbound)==2, must correct upper and lower layer so that divb=0 fix is correct, I do this instead of doing 2 loops and setting divb=0 cond in 2nd loop

void bound_mpi(FTYPE (*vars)[N2M][N1M],
	   FTYPE (*varv)[N3M][N2M][N1M],
	   int wsca,
	   int wvec,
	   int wcom) // 1,2,3, 0: none, 123=all 12=12 13=13, 23=23 (not currently setup to do 13 or 23 since never needed)
{
  static int firsttime=1;
  int tagheader,tagsend,tagrecv;
  int othercpu;
  int jumpfactor;
  int i,j,k,l,m,p,q ;
  int i2,j2,k2;
  int ii,jj,kk,ll;
  int jje2,jjje2;
  int iio,jjo,kko,llo;
  int iio2,jjo2,kko2,llo2;
  int iir,jjr,kkr,llr;
  int iii,jjj,kkk,lll;
  int iii2,jjj2,kkk2,lll2;
  int bct,bcdim,bcdir,bcd1,bcd2;
  int bcdir1,bcdir2,bcdir3;
  int bi,bj,bk;
  int bireal,bjreal,bkreal;
  int diri,dirj;
  int startloop,endloop,loopvar;
  FTYPE slope;
  FTYPE (*works)[N3M][N2M][N1M];
  FTYPE (*workv)[3][N3M][N2M][N1M];
  static int numhit=0;
  int numhitmin,numhitmax;
  int wbound;
  FTYPE ftemp;
  char temps[50];
  int tempslen=40;
  int tempi;
  int forces=0, forcev=0, forcecheck;
  int docom[3+1], comlength,comorder[3+1],itemp;

  int dosign;

  int dualloop,divbfix;
  int doingscalars,doingvectors;

  int gval;
  int typeofbc;
  int loopi,loopj,loopk;

  int othercomp;
  int cornerto,corner;
  int edgeto, edge;
#if(USEMPI)
  // should change this so just per cpu(modify tagsend/tagrecv to include another variable for these arrays)
  MPI_Request rrequest[COMPDIM*2*(NUMSCA+NUMVEC)*MAXCPUS]; // 2*dim directions, 1 per scalar and vector for each cpu(receive)
  MPI_Request srequest[COMPDIM*2*(NUMSCA+NUMVEC)*MAXCPUS]; // 2*dim directions, 1 per scalar and vector for each cpu(send)
#endif

  // might not want to bound vectors if not advecting them
  //  if(wvec==-1) return;
  //if(wvec==1) return;


  numhit++;
  if(wsca<=-2){
    works=(FTYPE (*) [N3M][N2M][N1M])(&vars[0][0][0]);
  }
  else works=s;
  if(wvec<=-2){
    workv=(FTYPE (*) [3][N3M][N2M][N1M])(&varv[0][0][0][0]);
  }
  else workv=v;


  ///////////////////////////////
  ////////////////////////////////
  // MPI
  ///////////////////////////////
  ///////////////////////////////

#if(USEMPI)

#if(DOBOUNDVEC)   
    // determine which components to do, if any
    if( (wcom==1)||(wcom==12)||(wcom==13)||(wcom==123) ) docom[1]=1; else docom[1]=0;
    if( (wcom==2)||(wcom==12)||(wcom==23)||(wcom==123) ) docom[2]=1; else docom[2]=0;
    if( (wcom==3)||(wcom==13)||(wcom==23)||(wcom==123) ) docom[3]=1; else docom[3]=0;
    if(wcom==0){ docom[1]=docom[2]=docom[3]=0; }

    // determine number of components to send
    if(wcom==0){ comlength=0; comorder[1]=0; comorder[2]=0; comorder[3]=0;} // really never used
    else if(wcom==1){ comlength=1; comorder[1]=1; comorder[2]=0; comorder[3]=0;} // 1
    else if(wcom==2){ comlength=1; comorder[1]=2; comorder[2]=0; comorder[3]=0;}
    else if(wcom==3){ comlength=1; comorder[1]=3; comorder[2]=0; comorder[3]=0;}
    else if(wcom==12){ comlength=2; comorder[1]=1; comorder[2]=2; comorder[3]=0; } // 2
    else if(wcom==13){ comlength=2; comorder[1]=1; comorder[2]=3;  comorder[3]=0; }
    else if(wcom==23){ comlength=2; comorder[1]=2; comorder[2]=3;  comorder[3]=0; }
    else if(wcom==123){ comlength=3; comorder[1]=1; comorder[2]=2; comorder[3]=3; } // 3

#endif // endif doboundvec==1



   
    // now do MPI stuff.  Done seperately so transfers are highly parallel.
    // 1) first did(already) normal bounding
    // 2) then do packing and send/recv commands
    // 3) then do all RECV waits/unpacks per recv transfer(don't do all recv's at once)
    // 4) then in seperate loop do all SEND waits
    
    // need to break this up into LF first then UD next so corner zones are correctly assigned

    // So just 2 of the above, one for LF, one for UD, and in 3D, one for IN/OUT

    // note that receive : [2][x], x is direction as sent by other cpu.
    // sends: [1][x], x is send direction

    /////////////////////////
    // doing left(2) and right(0)

#if(COMPDIM>=1)

#if(DOBOUNDSCA)

  if(wsca!=0){
    for(l=1;l<=NUMSCA-NOBOUNDPOT;l++){
      /* if not to do all, pick */
      if(wsca!=-1){ if(wsca<=-2) ll=0; else ll=wsca; }
      else ll=l;

      // determine how things are bounded for other vars
      if(ll==0){
	if(wsca==-2){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-11){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-12){//other
	  // bound like en
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0 for scalar case: wsca=%d\n",wsca);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*ll;
      // now that normal bcs are assigned, exchange with other cpus
      // tags are defined by sender's ID and direction sent
      // tag=(myid*COMPDIM*2)+{0,1,2,3,4,5}
      // 0=right, 1=up,2=left,3=down,4=out,5=in
      // works/v bc[1=output/2=input][0,1,2,3,4,5]
      // so sends are like: (sendtoid,myid*COMPDIM*2+?) and recv's are like: (fromid,otherid*COMPDIM*2+*) where ? and * are opposites(i.e. 0 and 2, 1 and 3, 4 and 5)

      //////////////// periodicx1 transfer
      if(mpiperiodicx1){      // ncpux1>1 to be here
	if(mycpupos[1]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  // pack data to be sent
	  
	  for(k=-N3BND;k<N3+N3BND;k++)	  for(j=-N2BND;j<N2+N2BND;j++){
	    worksbc[wbound][1][2][j+N2BND+(k+N3BND)*N2M*N1BND]=works[ll][k][j][0];
	    worksbc[wbound][1][2][j+N2M+N2BND+(k+N3BND)*N2M*N1BND]=works[ll][k][j][1];
	  }
	  // now send/receive data
	  othercpu=myid+(ncpux1-1); // the periodic mate for x1
	  tagrecv=tagheader+othercpu*COMPDIM*2+0;
	  tagsend=tagheader+myid*COMPDIM*2+2;

	  MPI_Irecv(worksbc[wbound][2][0],N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);
	  MPI_Isend(worksbc[wbound][1][2],N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);
	}
	else if(mycpupos[1]==ncpux1-1){
	  // pack data to be sent
	  for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
	    worksbc[wbound][1][0][j+N2BND+(k+N3BND)*N2M*N1BND]=works[ll][k][j][N1-2];
	    worksbc[wbound][1][0][j+N2M+N2BND+(k+N3BND)*N2M*N1BND]=works[ll][k][j][N1-1];
	  }
	  othercpu=myid-(ncpux1-1); // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+0;
	  tagrecv=tagheader+othercpu*COMPDIM*2+2;

	  MPI_Isend(worksbc[wbound][1][0],N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);
	  MPI_Irecv(worksbc[wbound][2][2],N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);

	}// end if i+ boundary and not on i- boundary
      }// end if periodicx1
      


      if(srdir[2]){ // do left

	//	fprintf(log_file,"srdir[2]=%d\n",srdir[2]); fflush(log_file);


	// pack data
	 for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
	  worksbc[wbound][1][2][j+N2BND+(k+N3BND)*N2M*N1BND]=works[ll][k][j][0];
	  worksbc[wbound][1][2][j+N2M+N2BND+(k+N3BND)*N2M*N1BND]=works[ll][k][j][1];
	}
	// send/recv data
	othercpu=myid-1;
	tagrecv=tagheader+othercpu*COMPDIM*2+0;
	tagsend=tagheader+myid*COMPDIM*2+2;
	
	MPI_Irecv(worksbc[wbound][2][0],N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);
	MPI_Isend(worksbc[wbound][1][2],N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);

      }

      if(srdir[0]){ // do right

	//	fprintf(log_file,"srdir[0]=%d\n",srdir[0]); fflush(log_file);

	
	for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
	  worksbc[wbound][1][0][j+N2BND+(k+N3BND)*N2M*N1BND]=works[ll][k][j][N1-2];
	  worksbc[wbound][1][0][j+N2M+N2BND+(k+N3BND)*N2M*N1BND]=works[ll][k][j][N1-1];
	}
	// send/recv data
	othercpu=myid+1;
	tagsend=tagheader+myid*COMPDIM*2+0;
	tagrecv=tagheader+othercpu*COMPDIM*2+2;

	MPI_Isend(worksbc[wbound][1][0],N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);
	MPI_Irecv(worksbc[wbound][2][2],N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);

      }
      
#if(DEBUGMPI>0)
      // first check if correct place in code(debug)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d s[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      
#endif

      /* cut short loop if only to do one */
      if(wsca!=-1) l=NUMSCA;
    } // end over scalars
  } // end if any scalars

#endif // end if doboundsca==1


#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after scalar-1\n",numhit,myidtxt);
#endif




  
#if(DOBOUNDVEC)  
  // Now do vectors if any
  if(wvec!=0){
    for(l=1;l<=REALNUMBOUNDVEC;l++){
   
      /* if not to do all, pick */
      if(wvec!=-1){ if(wvec<=-2) ll=0; else ll=wvec; }
      else ll=l;
      
      if(ll==0){
	if((wvec==-3)||(wvec==-4)){//mdot or flux(2) of v[1][2]
	  // bound like v
	  wbound=1;
	}
	else if((wvec==-5)){//mdot or flux(2) of v[1][2] or magemf
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-11){// v1
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-12){// v2
	  // bound like B
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0, wvec=%d\n",wvec);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(NUMSCA+ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*(NUMSCA+ll);

      /////////////////// periodicx1 transfer
      if(mpiperiodicx1){ // ncpux1>1 to be here
	if(mycpupos[1]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  // pack data to be sent
	  
	  for(kk=1;kk<=3;kk++){ // hit 1 to 3 times
	    if(comlength>=kk){
	      itemp=comorder[kk];
	      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
		workvbc[wbound][1][2][j+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))]=workv[ll][itemp][k][j][0];
		workvbc[wbound][1][2][j+N2M+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))]=workv[ll][itemp][k][j][1];
	      }
	    }
	  }
	  // now send/receive data
	  othercpu=myid+(ncpux1-1);
	  tagrecv=tagheader+othercpu*COMPDIM*2+0;
	  tagsend=tagheader+myid*COMPDIM*2+2;

	  MPI_Irecv(workvbc[wbound][2][0],comlength*N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);
	  MPI_Isend(workvbc[wbound][1][2],comlength*N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);

	}
	else if(mycpupos[1]==ncpux1-1){
	  // pack data to be sent
	  
	  for(kk=1;kk<=3;kk++){
	    if(comlength>=kk){
	      itemp=comorder[kk];
	      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
		workvbc[wbound][1][0][j+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))]=workv[ll][itemp][k][j][N1-2];
		workvbc[wbound][1][0][j+N2M+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))]=workv[ll][itemp][k][j][N1-1];
	      }
	    }
	  }
	  // now send/receive data
	  othercpu=myid-(ncpux1-1); // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+0;
	  tagrecv=tagheader+othercpu*COMPDIM*2+2;

	  MPI_Isend(workvbc[wbound][1][0],comlength*N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);
	  MPI_Irecv(workvbc[wbound][2][2],comlength*N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);

	}// end if i+ and not i- boundary
      }// end if transfering periodicx1



      
      ///////////////////// exchange with other cpus, normal interior transfer


      if(srdir[2]){ // do left
	// pack data

	//	fprintf(log_file,"vec srdir[2]=%d\n",srdir[2]); fflush(log_file);

	
	for(kk=1;kk<=3;kk++){
	  if(comlength>=kk){
	    itemp=comorder[kk];
	    for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
	      workvbc[wbound][1][2][j+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))]=workv[ll][itemp][k][j][0];
	      workvbc[wbound][1][2][j+N2M+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))]=workv[ll][itemp][k][j][1];
	    }
	  }
	}
	// send/recv data
        othercpu=myid-1;
	tagrecv=tagheader+othercpu*COMPDIM*2+0;
	tagsend=tagheader+myid*COMPDIM*2+2;

	MPI_Irecv(workvbc[wbound][2][0],comlength*N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);
	MPI_Isend(workvbc[wbound][1][2],comlength*N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);

      }

      if(srdir[0]){ // go right
	// pack data

	//	fprintf(log_file,"vec srdir[2]=%d\n",srdir[2]); fflush(log_file);

	
	for(kk=1;kk<=3;kk++){
	  if(comlength>=kk){
	    itemp=comorder[kk];
	    for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
	      workvbc[wbound][1][0][j+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))]=workv[ll][itemp][k][j][N1-2];
	      workvbc[wbound][1][0][j+N2M+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))]=workv[ll][itemp][k][j][N1-1];
	    }
	  }
	}
	// send/recv data
        othercpu=myid+1;
	tagsend=tagheader+myid*COMPDIM*2+0;
	tagrecv=tagheader+othercpu*COMPDIM*2+2;

	MPI_Isend(workvbc[wbound][1][0],comlength*N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);
	MPI_Irecv(workvbc[wbound][2][2],comlength*N1BND*N2M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);

      }
	

      // first check if correct place in code(debug)

#if(DEBUGMPI>0)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d v[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      

#endif



    
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMBOUNDVEC;
    }// end over vectors
  }// endif vectors to be done
  
#endif // end if do bound vectors


#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after vector-1\n",numhit,myidtxt);
#endif


  ///////////////////////////////////////////
  ///////////////////////////////////////////
  ///////////////////////////////////////////
  ///////////////////////////////////////////
  // now wait/unpack on each set(not wait or unpack all at once!)
  // exact same loops as above
  ///////////////////////////////////////////
  ///////////////////////////////////////////

#if(DOBOUNDSCA)

  if(wsca!=0){
    for(l=1;l<=NUMSCA-NOBOUNDPOT;l++){
      /* if not to do all, pick */
      if(wsca!=-1){ if(wsca<=-2) ll=0; else ll=wsca; }
      else ll=l;
      // determine how things are bounded for other vars
      if(ll==0){
	if(wsca==-2){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-11){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-12){//other
	  // bound like en
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0 for scalar case: wsca=%d\n",wsca);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*ll;
      // now that normal bcs are assigned, exchange with other cpus
      // tags are defined by sender's ID and direction sent
      // tag=(myid*COMPDIM*2)+{0,1,2,3,4,5}
      // 0=right, 1=up,2=left,3=down,4=out,5=in
      // works/v bc[1=output/2=input][0,1,2,3,4,5]
      // so sends are like: (sendtoid,myid*COMPDIM*2+?) and recv's are like: (fromid,otherid*COMPDIM*2+*) where ? and * are opposites(i.e. 0 and 2, 1 and 3)

      //////////////// periodicx1 transfer
      if(mpiperiodicx1){
	if(mycpupos[1]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  othercpu=myid+(ncpux1-1); // the periodic mate for x1
	  tagrecv=tagheader+othercpu*COMPDIM*2+0;
	  tagsend=tagheader+myid*COMPDIM*2+2;

	  MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	  // now assign transfered values back into data
	  
	  for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
	    works[ll][k][j][-2]=worksbc[wbound][2][0][j+N2BND+(k+N3BND)*N2M*N1BND];
	    works[ll][k][j][-1]=worksbc[wbound][2][0][j+N2M+N2BND+(k+N3BND)*N2M*N1BND];
	  }
	}
	else if(mycpupos[1]==ncpux1-1){
	  othercpu=myid-(ncpux1-1); // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+0;
	  tagrecv=tagheader+othercpu*COMPDIM*2+2;

	  MPI_Wait(&rrequest[tagrecv],&mpichstatus);
	  // unpack
	  
	  for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
	    works[ll][k][j][N1]=worksbc[wbound][2][2][j+N2BND+(k+N3BND)*N2M*N1BND];	  
	    works[ll][k][j][N1+1]=worksbc[wbound][2][2][j+N2M+N2BND+(k+N3BND)*N2M*N1BND];
	  }
	}// end if i+ boundary and not on i- boundary
      }// end if periodicx1


      if(srdir[2]){ // do left

	//	fprintf(log_file,"unpack sca srdir[2]=%d\n",srdir[2]); fflush(log_file);


	othercpu=myid-1;
	tagrecv=tagheader+othercpu*COMPDIM*2+0;
	tagsend=tagheader+myid*COMPDIM*2+2;
	
	MPI_Wait(&rrequest[tagrecv],&mpichstatus);
	// unpack
	
	for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
	  works[ll][k][j][-2]=worksbc[wbound][2][0][j+N2BND+(k+N3BND)*N2M*N1BND];
	  works[ll][k][j][-1]=worksbc[wbound][2][0][j+N2M+N2BND+(k+N3BND)*N2M*N1BND];
	}
      }

      if(srdir[0]){ // do right

	//	fprintf(log_file,"unpack sca srdir[0]=%d\n",srdir[0]); fflush(log_file);

	othercpu=myid+1;
	tagsend=tagheader+myid*COMPDIM*2+0;
	tagrecv=tagheader+othercpu*COMPDIM*2+2;

	MPI_Wait(&rrequest[tagrecv],&mpichstatus);
	
	for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
	  works[ll][k][j][N1]=worksbc[wbound][2][2][j+N2BND+(k+N3BND)*N2M*N1BND];
	  works[ll][k][j][N1+1]=worksbc[wbound][2][2][j+N2M+N2BND+(k+N3BND)*N2M*N1BND];
	}
      }
      
#if(DEBUGMPI>0)
      // first check if correct place in code(debug)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d s[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      
#endif

      /* cut short loop if only to do one */
      if(wsca!=-1) l=NUMSCA;
    } // end over scalars
  } // end if any scalars

#endif // endif doboundsca==1

#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after scalar-2\n",numhit,myidtxt);
#endif
  



#if(DOBOUNDVEC)  
  // Now do vectors if any
  if(wvec!=0){
        for(l=1;l<=REALNUMBOUNDVEC;l++){
   
      /* if not to do all, pick */
      if(wvec!=-1){ if(wvec<=-2) ll=0; else ll=wvec; }
      else ll=l;
      
      if(ll==0){
	if((wvec==-3)||(wvec==-4)){//mdot or flux(2) of v[1][2]
	  // bound like v
	  wbound=1;
	}
	else if((wvec==-5)){//mdot or flux(2) of v[1][2] or magemf
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-11){// v1
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-12){// v2
	  // bound like B
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0, wvec=%d\n",wvec);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(NUMSCA+ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*(NUMSCA+ll);

      /////////////////// periodicx1 transfer
      if(mpiperiodicx1){
	if(mycpupos[1]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  othercpu=myid+(ncpux1-1);
	  tagrecv=tagheader+othercpu*COMPDIM*2+0;
	  tagsend=tagheader+myid*COMPDIM*2+2;

	  MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	  // unpack
	  
	  for(kk=1;kk<=3;kk++){
	    if(comlength>=kk){
	      itemp=comorder[kk];
	      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
		workv[ll][itemp][k][j][-2]=workvbc[wbound][2][0][j+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))];
		workv[ll][itemp][k][j][-1]=workvbc[wbound][2][0][j+N2M+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))];
	      }
	    }
	  }
	}
	else if(mycpupos[1]==ncpux1-1){
	  othercpu=myid-(ncpux1-1); // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+0;
	  tagrecv=tagheader+othercpu*COMPDIM*2+2;

	  MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	  // unpack
	  
	  for(kk=1;kk<=3;kk++){
	    if(comlength>=kk){
	      itemp=comorder[kk];
	      for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
		workv[ll][itemp][k][j][N1]=workvbc[wbound][2][2][j+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))];
		workv[ll][itemp][k][j][N1+1]=workvbc[wbound][2][2][j+N2M+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))];
	      }
	    }
	  }
	}// end if i+ and not i- boundary
      }// end if transfering periodicx1


      
      ///////////////////// exchange with other cpus, normal interior transfer

      if(srdir[2]){ // do left

	//	fprintf(log_file,"unpack vec srdir[2]=%d\n",srdir[2]); fflush(log_file);


        othercpu=myid-1;
	tagrecv=tagheader+othercpu*COMPDIM*2+0;
	tagsend=tagheader+myid*COMPDIM*2+2;

	MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	
	for(kk=1;kk<=3;kk++){
	  if(comlength>=kk){
	    itemp=comorder[kk];
	    for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
	      workv[ll][itemp][k][j][-2]=workvbc[wbound][2][0][j+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))];
	      workv[ll][itemp][k][j][-1]=workvbc[wbound][2][0][j+N2M+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))];
	    }
	  }
	}
      }

      if(srdir[0]){ // go right

	//	fprintf(log_file,"unpack vec srdir[0]=%d\n",srdir[0]); fflush(log_file);

        othercpu=myid+1;
	tagsend=tagheader+myid*COMPDIM*2+0;
	tagrecv=tagheader+othercpu*COMPDIM*2+2;

	MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	// unpack data
	
	for(kk=1;kk<=3;kk++){
	  if(comlength>=kk){
	    itemp=comorder[kk];
	    for(k=-N3BND;k<N3+N3BND;k++) for(j=-N2BND;j<N2+N2BND;j++){
	      workv[ll][itemp][k][j][N1]=workvbc[wbound][2][2][j+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))];
	      workv[ll][itemp][k][j][N1+1]=workvbc[wbound][2][2][j+N2M+N2BND+(k+N3BND)*N2M*N1BND+(N1BND*N2M*N3M*(kk-1))];
	    }
	  }
	}
      }
	

      // first check if correct place in code(debug)

#if(DEBUGMPI>0)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d v[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      

#endif



    
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMBOUNDVEC;
    }// end over vectors
  }// endif vectors to be done
  
#endif // end if do bound vectors


#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after vector-2\n",numhit,myidtxt);
#endif


  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  // finally, wait on any sends left over before leaving
  // this is just same loops as above
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////

#if(DOBOUNDSCA)

  if(wsca!=0){
    for(l=1;l<=NUMSCA-NOBOUNDPOT;l++){
      /* if not to do all, pick */
      if(wsca!=-1){ if(wsca<=-2) ll=0; else ll=wsca; }
      else ll=l;
      // determine how things are bounded for other vars
      if(ll==0){
	if(wsca==-2){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-11){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-12){//other
	  // bound like en
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0 for scalar case: wsca=%d\n",wsca);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*ll;
      // now that normal bcs are assigned, exchange with other cpus
      // tags are defined by sender's ID and direction sent
      // tag=(myid*COMPDIM*2)+{0,1,2,3,4,5}
      // 0=right, 1=up,2=left,3=down,4=out,5=in
      // works/v bc[1=output/2=input][0,1,2,3,4,5]
      // so sends are like: (sendtoid,myid*COMPDIM*2+?) and recv's are like: (fromid,otherid*COMPDIM*2+*) where ? and * are opposites(i.e. 0 and 2, 1 and 3)

      //////////////// periodicx1 transfer
      if(mpiperiodicx1){
	if(mycpupos[1]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  othercpu=myid+(ncpux1-1); // the periodic mate for x1
	  tagrecv=tagheader+othercpu*COMPDIM*2+0;
	  tagsend=tagheader+myid*COMPDIM*2+2;

	  MPI_Wait(&srequest[tagsend],&mpichstatus);

	}
	else if(mycpupos[1]==ncpux1-1){
	  othercpu=myid-(ncpux1-1); // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+0;
	  tagrecv=tagheader+othercpu*COMPDIM*2+2;

	  MPI_Wait(&srequest[tagsend],&mpichstatus);

	}// end if i+ boundary and not on i- boundary
      }// end if periodicx1






      if(srdir[2]){ // do left
	othercpu=myid-1;
	tagrecv=tagheader+othercpu*COMPDIM*2+0;
	tagsend=tagheader+myid*COMPDIM*2+2;
	

	MPI_Wait(&srequest[tagsend],&mpichstatus);

      }

      if(srdir[0]){ // do right
	othercpu=myid+1;
	tagsend=tagheader+myid*COMPDIM*2+0;
	tagrecv=tagheader+othercpu*COMPDIM*2+2;

	MPI_Wait(&srequest[tagsend],&mpichstatus);

      }
      
#if(DEBUGMPI>0)
      // first check if correct place in code(debug)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d s[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      
#endif

      /* cut short loop if only to do one */
      if(wsca!=-1) l=NUMSCA;
    } // end over scalars
  } // end if any scalars

#endif // endif do boundsca==1

#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after scalar-3\n",numhit,myidtxt);
#endif

  
#if(DOBOUNDVEC)  
  // Now do vectors if any
  if(wvec!=0){
        for(l=1;l<=REALNUMBOUNDVEC;l++){
   
      /* if not to do all, pick */
      if(wvec!=-1){ if(wvec<=-2) ll=0; else ll=wvec; }
      else ll=l;
      
      if(ll==0){
	if((wvec==-3)||(wvec==-4)){//mdot or flux(2) of v[1][2]
	  // bound like v
	  wbound=1;
	}
	else if((wvec==-5)){//mdot or flux(2) of v[1][2] or magemf
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-11){// v1
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-12){// v2
	  // bound like B
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0, wvec=%d\n",wvec);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(NUMSCA+ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*(NUMSCA+ll);

      /////////////////// periodicx1 transfer
      if(mpiperiodicx1){
	if(mycpupos[1]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  othercpu=myid+(ncpux1-1);
	  tagrecv=tagheader+othercpu*COMPDIM*2+0;
	  tagsend=tagheader+myid*COMPDIM*2+2;

	  MPI_Wait(&srequest[tagsend],&mpichstatus);

	}
	else if(mycpupos[1]==ncpux1-1){
	  othercpu=myid-(ncpux1-1); // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+0;
	  tagrecv=tagheader+othercpu*COMPDIM*2+2;

	  MPI_Wait(&srequest[tagsend],&mpichstatus);

	}// end if i+ and not i- boundary
      }// end if transfering periodicx1





      
      ///////////////////// exchange with other cpus, normal interior transfer

      if(srdir[2]){ // do left
        othercpu=myid-1;
	tagrecv=tagheader+othercpu*COMPDIM*2+0;
	tagsend=tagheader+myid*COMPDIM*2+2;

	MPI_Wait(&srequest[tagsend],&mpichstatus);

      }

      if(srdir[0]){ // go right
        othercpu=myid+1;
	tagsend=tagheader+myid*COMPDIM*2+0;
	tagrecv=tagheader+othercpu*COMPDIM*2+2;

	MPI_Wait(&srequest[tagsend],&mpichstatus);

      }
	

      // first check if correct place in code(debug)

#if(DEBUGMPI>0)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d v[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      

#endif



    
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMBOUNDVEC;
    }// end over vectors
  }// endif vectors to be done
  
#endif // end if do bound vectors

#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after vector-3\n",numhit,myidtxt);
#endif

#endif // endif compdim>=1







#if(COMPDIM>=2)

  ////////////////////////
  ////////////////////////
  ////////////////////////
  // now do up(1) and down(3)!!!!!!!!!!!!!!!!!!!!!! OH YEAH!
  ////////////////////////




#if(DOBOUNDSCA)

  if(wsca!=0){
    for(l=1;l<=NUMSCA-NOBOUNDPOT;l++){
      /* if not to do all, pick */
      if(wsca!=-1){ if(wsca<=-2) ll=0; else ll=wsca; }
      else ll=l;

      // determine how things are bounded for other vars
      if(ll==0){
	if(wsca==-2){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-11){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-12){//other
	  // bound like en
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0 for scalar case: wsca=%d\n",wsca);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*ll;
      // now that normal bcs are assigned, exchange with other cpus
      // tags are defined by sender's ID and direction sent
      // tag=(myid*COMPDIM*2)+{0,1,2,3,4,5}
      // 0=right, 1=up,2=left,3=down,4=out,5=in
      // works/v bc[1=output/2=input][0,1,2,3,4,5]
      // so sends are like: (sendtoid,myid*COMPDIM*2+?) and recv's are like: (fromid,otherid*COMPDIM*2+*) where ? and * are opposites(i.e. 0 and 2, 1 and 3)

      //////////////// periodicx2 transfer
      if(mpiperiodicx2){ // ncpux2>1 to be here
	if(mycpupos[2]==0){ // then on inner j boundary
	  // now we have a cpu which needs to transfer to other boundary
	  // pack data to be sent
	  
	  for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
	    worksbc[wbound][1][1][i+N1BND+(k+N3BND)*N1M*N2BND]=works[ll][k][0][i];
	    worksbc[wbound][1][1][i+N1M+N1BND+(k+N3BND)*N1M*N2BND]=works[ll][k][1][i];
	  }
	  // now send/receive data
	  othercpu=myid+(ncpux2-1)*ncpux1;
	  tagrecv=tagheader+othercpu*COMPDIM*2+3;
	  tagsend=tagheader+myid*COMPDIM*2+1;

	  MPI_Irecv(worksbc[wbound][2][3],N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);
	  MPI_Isend(worksbc[wbound][1][1],N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);

	}
	else if(mycpupos[2]==ncpux2-1){
	  // pack data to be sent
	  
	  for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
	    worksbc[wbound][1][3][i+N1BND+(k+N3BND)*N1M*N2BND]=works[ll][k][N2-2][i];
	    worksbc[wbound][1][3][i+N1M+N1BND+(k+N3BND)*N1M*N2BND]=works[ll][k][N2-1][i];
	  }
	  othercpu=myid-(ncpux2-1)*ncpux1; // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+3;
	  tagrecv=tagheader+othercpu*COMPDIM*2+1;

	  MPI_Isend(worksbc[wbound][1][3],N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);
	  MPI_Irecv(worksbc[wbound][2][1],N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);

	}// end if j+ boundary and not on j- boundary
      }// end if periodicx2



      

      if(srdir[1]){ // do up
	// pack data
	
	for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
	  worksbc[wbound][1][1][i+N1BND+(k+N3BND)*N1M*N2BND]=works[ll][k][0][i];
	  worksbc[wbound][1][1][i+N1M+N1BND+(k+N3BND)*N1M*N2BND]=works[ll][k][1][i];
	}
	// send/recv data
	othercpu=myid-ncpux1;
	tagrecv=tagheader+othercpu*COMPDIM*2+3;
	tagsend=tagheader+myid*COMPDIM*2+1;

	MPI_Irecv(worksbc[wbound][2][3],N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);
	MPI_Isend(worksbc[wbound][1][1],N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);

      }

      if(srdir[3]){ // do down
	
	for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
	  worksbc[wbound][1][3][i+N1BND+(k+N3BND)*N1M*N2BND]=works[ll][k][N2-2][i];
	  worksbc[wbound][1][3][i+N1M+N1BND+(k+N3BND)*N1M*N2BND]=works[ll][k][N2-1][i];
	}
	// send/recv data
	othercpu=myid+ncpux1;
	tagsend=tagheader+myid*COMPDIM*2+3;
	tagrecv=tagheader+othercpu*COMPDIM*2+1;

	MPI_Isend(worksbc[wbound][1][3],N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);
	MPI_Irecv(worksbc[wbound][2][1],N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);

      }
      
#if(DEBUGMPI>0)
      // first check if correct place in code(debug)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d s[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      
#endif

      /* cut short loop if only to do one */
      if(wsca!=-1) l=NUMSCA;
    } // end over scalars
  } // end if any scalars

#endif // end if doboundsca==1


#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after scalar-1\n",numhit,myidtxt);
#endif




  
#if(DOBOUNDVEC)  
  // Now do vectors if any
  if(wvec!=0){
        for(l=1;l<=REALNUMBOUNDVEC;l++){
   
      /* if not to do all, pick */
      if(wvec!=-1){ if(wvec<=-2) ll=0; else ll=wvec; }
      else ll=l;
      
      if(ll==0){
	if((wvec==-3)||(wvec==-4)){//mdot or flux(2) of v[1][2]
	  // bound like v
	  wbound=1;
	}
	else if((wvec==-5)){//mdot or flux(2) of v[1][2] or magemf
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-11){// v1
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-12){// v2
	  // bound like B
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0, wvec=%d\n",wvec);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(NUMSCA+ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*(NUMSCA+ll);

      /////////////////// periodicx2 transfer
      if(mpiperiodicx2){ // has ncpux2>1 to be here
	if(mycpupos[2]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  // pack data to be sent
	  
	  for(kk=1;kk<=3;kk++){ // hit 1 to 3 times
	    if(comlength>=kk){
	      itemp=comorder[kk];
	      for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
		workvbc[wbound][1][1][i+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))]=workv[ll][itemp][k][0][i];
		workvbc[wbound][1][1][i+N1M+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))]=workv[ll][itemp][k][1][i];
	      }
	    }
	  }
	  // now send/receive data
	  othercpu=myid+(ncpux2-1)*ncpux1;
	  tagrecv=tagheader+othercpu*COMPDIM*2+3;
	  tagsend=tagheader+myid*COMPDIM*2+1;

	  MPI_Irecv(workvbc[wbound][2][3],comlength*N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);
	  MPI_Isend(workvbc[wbound][1][1],comlength*N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);

	}
	else if(mycpupos[2]==ncpux2-1){
	  // pack data to be sent
	  
	  for(kk=1;kk<=3;kk++){
	    if(comlength>=kk){
	      itemp=comorder[kk];
	      for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
		workvbc[wbound][1][3][i+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))]=workv[ll][itemp][k][N2-2][i];
		workvbc[wbound][1][3][i+N1M+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))]=workv[ll][itemp][k][N2-1][i];
	      }
	    }
	  }
	  // now send/receive data
	  othercpu=myid-(ncpux2-1)*ncpux1; // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+3;
	  tagrecv=tagheader+othercpu*COMPDIM*2+1;

	  MPI_Isend(workvbc[wbound][1][3],comlength*N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);
	  MPI_Irecv(workvbc[wbound][2][1],comlength*N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);

	}// end if j+ and not j- boundary
      }// end if transfering periodicx2

      
      ///////////////////// exchange with other cpus, normal interior transfer




      if(srdir[1]){ // do up
	// pack data
	
	for(kk=1;kk<=3;kk++){
	  if(comlength>=kk){
	    itemp=comorder[kk];
	    for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
	      workvbc[wbound][1][1][i+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))]=workv[ll][itemp][k][0][i];
	      workvbc[wbound][1][1][i+N1M+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))]=workv[ll][itemp][k][1][i];
	    }
	  }
	}
	// send/recv data
        othercpu=myid-ncpux1;
	tagrecv=tagheader+othercpu*COMPDIM*2+3;
	tagsend=tagheader+myid*COMPDIM*2+1;

	MPI_Irecv(workvbc[wbound][2][3],comlength*N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);
	MPI_Isend(workvbc[wbound][1][1],comlength*N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);

      }

      if(srdir[3]){ // do down

	// pack data
	
	for(kk=1;kk<=3;kk++){
	  if(comlength>=kk){
	    itemp=comorder[kk];
	    for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
	      workvbc[wbound][1][3][i+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))]=workv[ll][itemp][k][N2-2][i];
	      workvbc[wbound][1][3][i+N1M+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))]=workv[ll][itemp][k][N2-1][i];
	    }
	  }
	}
	// send/recv data
        othercpu=myid+ncpux1;
	tagsend=tagheader+myid*COMPDIM*2+3;
	tagrecv=tagheader+othercpu*COMPDIM*2+1;

	MPI_Isend(workvbc[wbound][1][3],comlength*N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);
	MPI_Irecv(workvbc[wbound][2][1],comlength*N2BND*N1M*N3M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);

      }
	

      // first check if correct place in code(debug)

#if(DEBUGMPI>0)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d v[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      

#endif



    
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMBOUNDVEC;
    }// end over vectors
  }// endif vectors to be done
  
#endif // end if do bound vectors


#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after vector-1\n",numhit,myidtxt);
#endif


  ///////////////////////////////////////////
  ///////////////////////////////////////////
  ///////////////////////////////////////////
  ///////////////////////////////////////////
  // now wait/unpack on each set(not wait or unpack all at once!)
  // exact same loops as above
  ///////////////////////////////////////////
  ///////////////////////////////////////////

#if(DOBOUNDSCA)

  if(wsca!=0){
    for(l=1;l<=NUMSCA-NOBOUNDPOT;l++){
      /* if not to do all, pick */
      if(wsca!=-1){ if(wsca<=-2) ll=0; else ll=wsca; }
      else ll=l;
      // determine how things are bounded for other vars
      if(ll==0){
	if(wsca==-2){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-11){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-12){//other
	  // bound like en
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0 for scalar case: wsca=%d\n",wsca);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*ll;
      // now that normal bcs are assigned, exchange with other cpus
      // tags are defined by sender's ID and direction sent
      // tag=(myid*COMPDIM*2)+{0,1,2,3,4,5}
      // 0=right, 1=up,2=left,3=down,4=out,5=in
      // works/v bc[1=output/2=input][0,1,2,3,4,5]
      // so sends are like: (sendtoid,myid*COMPDIM*2+?) and recv's are like: (fromid,otherid*COMPDIM*2+*) where ? and * are opposites(i.e. 0 and 2, 1 and 3)

      //////////////// periodicx2 transfer
      if(mpiperiodicx2){
	if(mycpupos[2]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  othercpu=myid+(ncpux2-1)*ncpux1;
	  tagrecv=tagheader+othercpu*COMPDIM*2+3;
	  tagsend=tagheader+myid*COMPDIM*2+1;

	  MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	  // now assign transfered values back into data
	  
	  for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
	    works[ll][k][-2][i]=worksbc[wbound][2][3][i+N1BND+(k+N3BND)*N1M*N2BND];
	    works[ll][k][-1][i]=worksbc[wbound][2][3][i+N1M+N1BND+(k+N3BND)*N1M*N2BND];
	  }
	}
	else if(mycpupos[2]==ncpux2-1){
	  othercpu=myid-(ncpux2-1)*ncpux1; // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+3;
	  tagrecv=tagheader+othercpu*COMPDIM*2+1;

	  MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	  // unpack
	  
	  for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
	    works[ll][k][N2][i]=worksbc[wbound][2][1][i+N1BND+(k+N3BND)*N1M*N2BND];	  
	    works[ll][k][N2+1][i]=worksbc[wbound][2][1][i+N1M+N1BND+(k+N3BND)*N1M*N2BND];
	  }
	}// end if j+ boundary and not on j- boundary
      }// end if periodicx2



      

      if(srdir[1]){ // do up
	othercpu=myid-ncpux1;
	tagrecv=tagheader+othercpu*COMPDIM*2+3;
	tagsend=tagheader+myid*COMPDIM*2+1;

	MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	// unpack
	
	for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
	  works[ll][k][-2][i]=worksbc[wbound][2][3][i+N1BND+(k+N3BND)*N1M*N2BND];
	  works[ll][k][-1][i]=worksbc[wbound][2][3][i+N1M+N1BND+(k+N3BND)*N1M*N2BND];
	}
      }

      if(srdir[3]){ // do down
	othercpu=myid+ncpux1;
	tagsend=tagheader+myid*COMPDIM*2+3;
	tagrecv=tagheader+othercpu*COMPDIM*2+1;

	MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	// unpack
	
	for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
	  works[ll][k][N2][i]=worksbc[wbound][2][1][i+N1BND+(k+N3BND)*N1M*N2BND];
	  works[ll][k][N2+1][i]=worksbc[wbound][2][1][i+N1M+N1BND+(k+N3BND)*N1M*N2BND];
	}
      }
      
#if(DEBUGMPI>0)
      // first check if correct place in code(debug)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d s[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      
#endif

      /* cut short loop if only to do one */
      if(wsca!=-1) l=NUMSCA;
    } // end over scalars
  } // end if any scalars

#endif // endif doboundsca==1

#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after scalar-2\n",numhit,myidtxt);
#endif
  



#if(DOBOUNDVEC)  
  // Now do vectors if any
  if(wvec!=0){
        for(l=1;l<=REALNUMBOUNDVEC;l++){
   
      /* if not to do all, pick */
      if(wvec!=-1){ if(wvec<=-2) ll=0; else ll=wvec; }
      else ll=l;
      
      if(ll==0){
	if((wvec==-3)||(wvec==-4)){//mdot or flux(2) of v[1][2]
	  // bound like v
	  wbound=1;
	}
	else if((wvec==-5)){//mdot or flux(2) of v[1][2] or magemf
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-11){// v1
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-12){// v2
	  // bound like B
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0, wvec=%d\n",wvec);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(NUMSCA+ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*(NUMSCA+ll);




      /////////////////// periodicx2 transfer
      if(mpiperiodicx2){
	if(mycpupos[2]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  othercpu=myid+(ncpux2-1)*ncpux1;
	  tagrecv=tagheader+othercpu*COMPDIM*2+3;
	  tagsend=tagheader+myid*COMPDIM*2+1;

	  MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	  // unpack
	  
	  for(kk=1;kk<=3;kk++){
	    if(comlength>=kk){
	      itemp=comorder[kk];
	      for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
		workv[ll][itemp][k][-2][i]=workvbc[wbound][2][3][i+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))];
		workv[ll][itemp][k][-1][i]=workvbc[wbound][2][3][i+N1M+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))];
	      }
	    }
	  }
	}
	else if(mycpupos[2]==ncpux2-1){
	  othercpu=myid-(ncpux2-1)*ncpux1; // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+3;
	  tagrecv=tagheader+othercpu*COMPDIM*2+1;

	  MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	  // unpack
	  
	  for(kk=1;kk<=3;kk++){
	    if(comlength>=kk){
	      itemp=comorder[kk];
	      for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
		workv[ll][itemp][k][N2][i]=workvbc[wbound][2][1][i+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))];
		workv[ll][itemp][k][N2+1][i]=workvbc[wbound][2][1][i+N1M+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))];
	      }
	    }
	  }
	}// end if j+ and not j- boundary
      }// end if transfering periodicx2

      
      ///////////////////// exchange with other cpus, normal interior transfer

      if(srdir[1]){ // do up
        othercpu=myid-ncpux1;
	tagrecv=tagheader+othercpu*COMPDIM*2+3;
	tagsend=tagheader+myid*COMPDIM*2+1;

	MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	
	for(kk=1;kk<=3;kk++){
	  if(comlength>=kk){
	    itemp=comorder[kk];
	    for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
	      workv[ll][itemp][k][-2][i]=workvbc[wbound][2][3][i+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))];
	      workv[ll][itemp][k][-1][i]=workvbc[wbound][2][3][i+N1M+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))];
	    }
	  }
	}
      }
      if(srdir[3]){ // do down
        othercpu=myid+ncpux1;
	tagsend=tagheader+myid*COMPDIM*2+3;
	tagrecv=tagheader+othercpu*COMPDIM*2+1;

	MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	// unpack data
	
	for(kk=1;kk<=3;kk++){
	  if(comlength>=kk){
	    itemp=comorder[kk];
	    for(k=-N3BND;k<N3+N3BND;k++) for(i=-N1BND;i<N1+N1BND;i++){
	      workv[ll][itemp][k][N2][i]=workvbc[wbound][2][1][i+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))];
	      workv[ll][itemp][k][N2+1][i]=workvbc[wbound][2][1][i+N1M+N1BND+(k+N3BND)*N1M*N2BND+(N2BND*N1M*N3M*(kk-1))];
	    }
	  }
	}
      }
	

      // first check if correct place in code(debug)

#if(DEBUGMPI>0)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d v[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      

#endif



    
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMBOUNDVEC;
    }// end over vectors
  }// endif vectors to be done
  
#endif // end if do bound vectors


#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after vector-2\n",numhit,myidtxt);
#endif


  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  // finally, wait on any sends left over before leaving
  // this is just same loops as above
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////

#if(DOBOUNDSCA)

  if(wsca!=0){
    for(l=1;l<=NUMSCA-NOBOUNDPOT;l++){
      /* if not to do all, pick */
      if(wsca!=-1){ if(wsca<=-2) ll=0; else ll=wsca; }
      else ll=l;
      // determine how things are bounded for other vars
      if(ll==0){
	if(wsca==-2){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-11){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-12){//other
	  // bound like en
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0 for scalar case: wsca=%d\n",wsca);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*ll;
      // now that normal bcs are assigned, exchange with other cpus
      // tags are defined by sender's ID and direction sent
      // tag=(myid*COMPDIM*2)+{0,1,2,3,4,5}
      // 0=right, 1=up,2=left,3=down,4=out,5=in
      // works/v bc[1=output/2=input][0,1,2,3,4,5]
      // so sends are like: (sendtoid,myid*COMPDIM*2+?) and recv's are like: (fromid,otherid*COMPDIM*2+*) where ? and * are opposites(i.e. 0 and 2, 1 and 3)


      //////////////// periodicx2 transfer
      if(mpiperiodicx2){
	if(mycpupos[2]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  othercpu=myid+(ncpux2-1)*ncpux1;
	  tagrecv=tagheader+othercpu*COMPDIM*2+3;
	  tagsend=tagheader+myid*COMPDIM*2+1;

	  MPI_Wait(&srequest[tagsend],&mpichstatus);

	}
	else if(mycpupos[2]==ncpux2-1){
	  othercpu=myid-(ncpux2-1)*ncpux1; // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+3;
	  tagrecv=tagheader+othercpu*COMPDIM*2+1;

	  MPI_Wait(&srequest[tagsend],&mpichstatus);

	}// end if j+ boundary and not on j- boundary
      }// end if periodicx2



      

      if(srdir[1]){ // do up
	othercpu=myid-ncpux1;
	tagrecv=tagheader+othercpu*COMPDIM*2+3;
	tagsend=tagheader+myid*COMPDIM*2+1;

	MPI_Wait(&srequest[tagsend],&mpichstatus);

      }

      if(srdir[3]){ // do down
	othercpu=myid+ncpux1;
	tagsend=tagheader+myid*COMPDIM*2+3;
	tagrecv=tagheader+othercpu*COMPDIM*2+1;

	MPI_Wait(&srequest[tagsend],&mpichstatus);

      }
      
#if(DEBUGMPI>0)
      // first check if correct place in code(debug)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d s[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      
#endif

      /* cut short loop if only to do one */
      if(wsca!=-1) l=NUMSCA;
    } // end over scalars
  } // end if any scalars

#endif // endif do boundsca==1

#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after scalar-3\n",numhit,myidtxt);
#endif

  
#if(DOBOUNDVEC)  
  // Now do vectors if any
  if(wvec!=0){
        for(l=1;l<=REALNUMBOUNDVEC;l++){
   
      /* if not to do all, pick */
      if(wvec!=-1){ if(wvec<=-2) ll=0; else ll=wvec; }
      else ll=l;
      
      if(ll==0){
	if((wvec==-3)||(wvec==-4)){//mdot or flux(2) of v[1][2]
	  // bound like v
	  wbound=1;
	}
	else if((wvec==-5)){//mdot or flux(2) of v[1][2] or magemf
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-11){// v1
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-12){// v2
	  // bound like B
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0, wvec=%d\n",wvec);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(NUMSCA+ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*(NUMSCA+ll);



      /////////////////// periodicx2 transfer
      if(mpiperiodicx2){
	if(mycpupos[2]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  othercpu=myid+(ncpux2-1)*ncpux1;
	  tagrecv=tagheader+othercpu*COMPDIM*2+3;
	  tagsend=tagheader+myid*COMPDIM*2+1;

	  MPI_Wait(&srequest[tagsend],&mpichstatus);

	}
	else if(mycpupos[2]==ncpux2-1){
	  othercpu=myid-(ncpux2-1)*ncpux1; // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+3;
	  tagrecv=tagheader+othercpu*COMPDIM*2+1;

	  MPI_Wait(&srequest[tagsend],&mpichstatus);

	}// end if j+ and not j- boundary
      }// end if transfering periodicx2

      
      ///////////////////// exchange with other cpus, normal interior transfer

      if(srdir[1]){ // do up
        othercpu=myid-ncpux1;
	tagrecv=tagheader+othercpu*COMPDIM*2+3;
	tagsend=tagheader+myid*COMPDIM*2+1;

	MPI_Wait(&srequest[tagsend],&mpichstatus);

      }
      if(srdir[3]){ // go down
        othercpu=myid+ncpux1;
	tagsend=tagheader+myid*COMPDIM*2+3;
	tagrecv=tagheader+othercpu*COMPDIM*2+1;

	MPI_Wait(&srequest[tagsend],&mpichstatus);

      }

      // first check if correct place in code(debug)

#if(DEBUGMPI>0)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d v[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      

#endif



    
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMBOUNDVEC;
    }// end over vectors
  }// endif vectors to be done
  
#endif // end if do bound vectors

#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after vector-3\n",numhit,myidtxt);
#endif

#endif // endif compdim>=2












#if(COMPDIM>=3)


  ////////////////////////
  ////////////////////////
  ////////////////////////
  // now do in(5) and out(4)!!!!!!!!!!!!!!!!!!!!!! OH YEAH!
  ////////////////////////




#if(DOBOUNDSCA)

  if(wsca!=0){
    for(l=1;l<=NUMSCA-NOBOUNDPOT;l++){
      /* if not to do all, pick */
      if(wsca!=-1){ if(wsca<=-2) ll=0; else ll=wsca; }
      else ll=l;

      // determine how things are bounded for other vars
      if(ll==0){
	if(wsca==-2){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-11){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-12){//other
	  // bound like en
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0 for scalar case: wsca=%d\n",wsca);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*ll;
      // now that normal bcs are assigned, exchange with other cpus
      // tags are defined by sender's ID and direction sent
      // tag=(myid*COMPDIM*2)+{0,1,2,3,4,5}
      // 0=right, 1=up,2=left,3=down,4=out,5=in
      // works/v bc[1=output/2=input][0,1,2,3,4,5]
      // so sends are like: (sendtoid,myid*COMPDIM*2+?) and recv's are like: (fromid,otherid*COMPDIM*2+*) where ? and * are opposites(i.e. 0 and 2, 1 and 3)

      //////////////// periodicx3 transfer
      if(mpiperiodicx3){ // ncpux3>1 to be here
	if(mycpupos[3]==0){ // then on inner k boundary
	  // now we have a cpu which needs to transfer to other boundary
	  // pack data to be sent
	  
	  for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	    worksbc[wbound][1][5][i+N1BND+(j+N2BND)*N1M*N3BND]=works[ll][0][j][i];
	    worksbc[wbound][1][5][i+N1M+N1BND+(j+N2BND)*N1M*N3BND]=works[ll][1][j][i];
	  }
	  // now send/receive data
	  othercpu=myid+(ncpux3-1)*ncpux1*ncpux2;
	  tagrecv=tagheader+othercpu*COMPDIM*2+4;
	  tagsend=tagheader+myid*COMPDIM*2+5;

	  MPI_Irecv(worksbc[wbound][2][4],N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);
	  MPI_Isend(worksbc[wbound][1][5],N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);

	}
	else if(mycpupos[3]==ncpux3-1){
	  // pack data to be sent
	  
	  for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	    worksbc[wbound][1][4][i+N1BND+(j+N2BND)*N1M*N3BND]=works[ll][N3-2][j][i];
	    worksbc[wbound][1][4][i+N1M+N1BND+(j+N2BND)*N1M*N3BND]=works[ll][N3-1][j][i];
	  }
	  othercpu=myid-(ncpux3-1)*ncpux1*ncpux2; // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+4;
	  tagrecv=tagheader+othercpu*COMPDIM*2+5;

	  MPI_Isend(worksbc[wbound][1][4],N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);
	  MPI_Irecv(worksbc[wbound][2][5],N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);

	}// end if k+ boundary and not on k- boundary
      }// end if periodicx3



      

      if(srdir[5]){ // do in
	// pack data
	
	for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	  worksbc[wbound][1][5][i+N1BND+(j+N2BND)*N1M*N3BND]=works[ll][0][j][i];
	  worksbc[wbound][1][5][i+N1M+N1BND+(j+N2BND)*N1M*N3BND]=works[ll][1][j][i];
	}
	// send/recv data
	othercpu=myid-ncpux1*ncpux2;
	tagrecv=tagheader+othercpu*COMPDIM*2+4;
	tagsend=tagheader+myid*COMPDIM*2+5;

	MPI_Irecv(worksbc[wbound][2][4],N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);
	MPI_Isend(worksbc[wbound][1][5],N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);

      }

      if(srdir[4]){ // do out
	
	for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	  worksbc[wbound][1][4][i+N1BND+(j+N2BND)*N1M*N3BND]=works[ll][N3-2][j][i];
	  worksbc[wbound][1][4][i+N1M+N1BND+(j+N2BND)*N1M*N3BND]=works[ll][N3-1][j][i];
	}
	// send/recv data
	othercpu=myid+ncpux1*ncpux2;
	tagsend=tagheader+myid*COMPDIM*2+4;
	tagrecv=tagheader+othercpu*COMPDIM*2+5;

	MPI_Isend(worksbc[wbound][1][4],N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);
	MPI_Irecv(worksbc[wbound][2][5],N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);

      }
      
#if(DEBUGMPI>0)
      // first check if correct place in code(debug)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d s[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      
#endif

      /* cut short loop if only to do one */
      if(wsca!=-1) l=NUMSCA;
    } // end over scalars
  } // end if any scalars

#endif // end if doboundsca==1


#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after scalar-1\n",numhit,myidtxt);
#endif




  
#if(DOBOUNDVEC)  
  // Now do vectors if any
  if(wvec!=0){
        for(l=1;l<=REALNUMBOUNDVEC;l++){
   
      /* if not to do all, pick */
      if(wvec!=-1){ if(wvec<=-2) ll=0; else ll=wvec; }
      else ll=l;
      
      if(ll==0){
	if((wvec==-3)||(wvec==-4)){//mdot or flux(2) of v[1][2]
	  // bound like v
	  wbound=1;
	}
	else if((wvec==-5)){//mdot or flux(2) of v[1][2] or magemf
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-11){// v1
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-12){// v2
	  // bound like B
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0, wvec=%d\n",wvec);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(NUMSCA+ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*(NUMSCA+ll);

      /////////////////// periodicx3 transfer
      if(mpiperiodicx3){ // has ncpux2>1 to be here
	if(mycpupos[3]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  // pack data to be sent
	  
	  for(kk=1;kk<=3;kk++){ // hit 1 to 3 times
	    if(comlength>=kk){
	      itemp=comorder[kk];
	      for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
		workvbc[wbound][1][5][i+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))]=workv[ll][itemp][0][j][i];
		workvbc[wbound][1][5][i+N1M+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))]=workv[ll][itemp][1][j][i];
	      }
	    }
	  }
	  // now send/receive data
	  othercpu=myid+(ncpux3-1)*ncpux1*ncpux2; // the periodic mate for x2
	  tagrecv=tagheader+othercpu*COMPDIM*2+4;
	  tagsend=tagheader+myid*COMPDIM*2+5;

	  MPI_Irecv(workvbc[wbound][2][4],comlength*N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);
	  MPI_Isend(workvbc[wbound][1][5],comlength*N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);

	}
	else if(mycpupos[3]==ncpux3-1){
	  // pack data to be sent
	  
	  for(kk=1;kk<=3;kk++){
	    if(comlength>=kk){
	      itemp=comorder[kk];
	      for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
		workvbc[wbound][1][4][i+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))]=workv[ll][itemp][N3-2][j][i];
		workvbc[wbound][1][4][i+N1M+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))]=workv[ll][itemp][N3-1][j][i];
	      }
	    }
	  }
	  // now send/receive data
	  othercpu=myid-(ncpux3-1)*ncpux1*ncpux2; // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+4;
	  tagrecv=tagheader+othercpu*COMPDIM*2+5;

	  MPI_Isend(workvbc[wbound][1][4],comlength*N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);
	  MPI_Irecv(workvbc[wbound][2][5],comlength*N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);

	}// end if j+ and not k- boundary
      }// end if transfering periodicx3

      
      ///////////////////// exchange with other cpus, normal interior transfer




      if(srdir[5]){ // do in
	// pack data
	
	for(kk=1;kk<=3;kk++){
	  if(comlength>=kk){
	    itemp=comorder[kk];
	    for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	      workvbc[wbound][1][5][i+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))]=workv[ll][itemp][0][j][i];
	      workvbc[wbound][1][5][i+N1M+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))]=workv[ll][itemp][1][j][i];
	    }
	  }
	}
	// send/recv data
        othercpu=myid-ncpux1*ncpux2;
	tagrecv=tagheader+othercpu*COMPDIM*2+4;
	tagsend=tagheader+myid*COMPDIM*2+5;

	MPI_Irecv(workvbc[wbound][2][4],comlength*N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);
	MPI_Isend(workvbc[wbound][1][5],comlength*N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);

      }

      if(srdir[4]){ // do out

	// pack data
	
	for(kk=1;kk<=3;kk++){
	  if(comlength>=kk){
	    itemp=comorder[kk];
	    for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	      workvbc[wbound][1][4][i+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))]=workv[ll][itemp][N3-2][j][i];
	      workvbc[wbound][1][4][i+N1M+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))]=workv[ll][itemp][N3-1][j][i];
	    }
	  }
	}
	// send/recv data
        othercpu=myid+ncpux1*ncpux2;
	tagsend=tagheader+myid*COMPDIM*2+4;
	tagrecv=tagheader+othercpu*COMPDIM*2+5;

	MPI_Isend(workvbc[wbound][1][4],comlength*N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD,&srequest[tagsend]);
	MPI_Irecv(workvbc[wbound][2][5],comlength*N3BND*N1M*N2M,MPI_FTYPE,othercpu,tagrecv,MPI_COMM_WORLD,&rrequest[tagrecv]);

      }
	

      // first check if correct place in code(debug)

#if(DEBUGMPI>0)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d v[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      

#endif



    
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMBOUNDVEC;
    }// end over vectors
  }// endif vectors to be done
  
#endif // end if do bound vectors


#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after vector-1\n",numhit,myidtxt);
#endif


  ///////////////////////////////////////////
  ///////////////////////////////////////////
  ///////////////////////////////////////////
  ///////////////////////////////////////////
  // now wait/unpack on each set(not wait or unpack all at once!)
  // exact same loops as above
  ///////////////////////////////////////////
  ///////////////////////////////////////////

#if(DOBOUNDSCA)

  if(wsca!=0){
    for(l=1;l<=NUMSCA-NOBOUNDPOT;l++){
      /* if not to do all, pick */
      if(wsca!=-1){ if(wsca<=-2) ll=0; else ll=wsca; }
      else ll=l;
      // determine how things are bounded for other vars
      if(ll==0){
	if(wsca==-2){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-11){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-12){//other
	  // bound like en
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0 for scalar case: wsca=%d\n",wsca);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*ll;
      // now that normal bcs are assigned, exchange with other cpus
      // tags are defined by sender's ID and direction sent
      // tag=(myid*COMPDIM*2)+{0,1,2,3,4,5}
      // 0=right, 1=up,2=left,3=down,4=out,5=in
      // works/v bc[1=output/2=input][0,1,2,3,4,5]
      // so sends are like: (sendtoid,myid*COMPDIM*2+?) and recv's are like: (fromid,otherid*COMPDIM*2+*) where ? and * are opposites(i.e. 0 and 2, 1 and 3)

      //////////////// periodicx3 transfer
      if(mpiperiodicx3){
	if(mycpupos[3]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  othercpu=myid+(ncpux3-1)*ncpux1*ncpux2;
	  tagrecv=tagheader+othercpu*COMPDIM*2+4;
	  tagsend=tagheader+myid*COMPDIM*2+5;

	  MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	  // now assign transfered values back into data
	  
	  for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	    works[ll][-2][j][i]=worksbc[wbound][2][4][i+N1BND+(j+N2BND)*N1M*N3BND];
	    works[ll][-1][j][i]=worksbc[wbound][2][4][i+N1M+N1BND+(j+N2BND)*N1M*N3BND];
	  }
	}
	else if(mycpupos[3]==ncpux3-1){
	  othercpu=myid-(ncpux3-1)*ncpux1*ncpux2; // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+4;
	  tagrecv=tagheader+othercpu*COMPDIM*2+5;

	  MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	  // unpack
	  
	  for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	    works[ll][N3][j][i]=worksbc[wbound][2][5][i+N1BND+(j+N2BND)*N1M*N3BND];	  
	    works[ll][N3+1][j][i]=worksbc[wbound][2][5][i+N1M+N1BND+(j+N2BND)*N1M*N3BND];
	  }
	}// end if k+ boundary and not on k- boundary
      }// end if periodicx3



      

      if(srdir[5]){ // do in
	othercpu=myid-ncpux1*ncpux2;
	tagrecv=tagheader+othercpu*COMPDIM*2+4;
	tagsend=tagheader+myid*COMPDIM*2+5;

	MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	// unpack
	
	for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	  works[ll][-2][j][i]=worksbc[wbound][2][4][i+N1BND+(j+N2BND)*N1M*N3BND];
	  works[ll][-1][j][i]=worksbc[wbound][2][4][i+N1M+N1BND+(j+N2BND)*N1M*N3BND];
	}
      }

      if(srdir[4]){ // do out
	othercpu=myid+ncpux1*ncpux2;
	tagsend=tagheader+myid*COMPDIM*2+4;
	tagrecv=tagheader+othercpu*COMPDIM*2+5;

	MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	// unpack
	
	for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	  works[ll][N3][j][i]=worksbc[wbound][2][5][i+N1BND+(j+N2BND)*N1M*N3BND];
	  works[ll][N3+1][j][i]=worksbc[wbound][2][5][i+N1M+N1BND+(j+N2BND)*N1M*N3BND];
	}
      }
      
#if(DEBUGMPI>0)
      // first check if correct place in code(debug)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d s[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      
#endif

      /* cut short loop if only to do one */
      if(wsca!=-1) l=NUMSCA;
    } // end over scalars
  } // end if any scalars

#endif // endif doboundsca==1

#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after scalar-2\n",numhit,myidtxt);
#endif
  



#if(DOBOUNDVEC)  
  // Now do vectors if any
  if(wvec!=0){
        for(l=1;l<=REALNUMBOUNDVEC;l++){
   
      /* if not to do all, pick */
      if(wvec!=-1){ if(wvec<=-2) ll=0; else ll=wvec; }
      else ll=l;
      
      if(ll==0){
	if((wvec==-3)||(wvec==-4)){//mdot or flux(2) of v[1][2]
	  // bound like v
	  wbound=1;
	}
	else if((wvec==-5)){//mdot or flux(2) of v[1][2] or magemf
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-11){// v1
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-12){// v2
	  // bound like B
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0, wvec=%d\n",wvec);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(NUMSCA+ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*(NUMSCA+ll);




      /////////////////// periodicx3 transfer
      if(mpiperiodicx3){
	if(mycpupos[3]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  othercpu=myid+(ncpux3-1)*ncpux1*ncpux2; // the periodic mate for x2
	  tagrecv=tagheader+othercpu*COMPDIM*2+4;
	  tagsend=tagheader+myid*COMPDIM*2+5;

	  MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	  // unpack
	  
	  for(kk=1;kk<=3;kk++){
	    if(comlength>=kk){
	      itemp=comorder[kk];
	      for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
		workv[ll][itemp][-2][j][i]=workvbc[wbound][2][4][i+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))];
		workv[ll][itemp][-1][j][i]=workvbc[wbound][2][4][i+N1M+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))];
	      }
	    }
	  }
	}
	else if(mycpupos[3]==ncpux3-1){
	  othercpu=myid-(ncpux3-1)*ncpux1*ncpux2; // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+4;
	  tagrecv=tagheader+othercpu*COMPDIM*2+5;

	  MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	  // unpack
	  
	  for(kk=1;kk<=3;kk++){
	    if(comlength>=kk){
	      itemp=comorder[kk];
	      for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
		workv[ll][itemp][N3][j][i]=workvbc[wbound][2][5][i+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))];
		workv[ll][itemp][N3+1][j][i]=workvbc[wbound][2][5][i+N1M+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))];
	      }
	    }
	  }
	}// end if k+ and not k- boundary
      }// end if transfering periodicx3

      
      ///////////////////// exchange with other cpus, normal interior transfer

      if(srdir[5]){ // do in
        othercpu=myid-ncpux1*ncpux2;
	tagrecv=tagheader+othercpu*COMPDIM*2+4;
	tagsend=tagheader+myid*COMPDIM*2+5;

	MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	
	for(kk=1;kk<=3;kk++){
	  if(comlength>=kk){
	    itemp=comorder[kk];
	    for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	      workv[ll][itemp][-2][j][i]=workvbc[wbound][2][4][i+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))];
	      workv[ll][itemp][-1][j][i]=workvbc[wbound][2][4][i+N1M+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))];
	    }
	  }
	}
      }
      if(srdir[4]){ // do out
        othercpu=myid+ncpux1*ncpux2;
	tagsend=tagheader+myid*COMPDIM*2+4;
	tagrecv=tagheader+othercpu*COMPDIM*2+5;

	MPI_Wait(&rrequest[tagrecv],&mpichstatus);

	// unpack data
	
	for(kk=1;kk<=3;kk++){
	  if(comlength>=kk){
	    itemp=comorder[kk];
	    for(j=-N2BND;j<N2+N2BND;j++) for(i=-N1BND;i<N1+N1BND;i++){
	      workv[ll][itemp][N3][j][i]=workvbc[wbound][2][5][i+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))];
	      workv[ll][itemp][N3+1][j][i]=workvbc[wbound][2][5][i+N1M+N1BND+(j+N2BND)*N1M*N3BND+(N3BND*N1M*N2M*(kk-1))];
	    }
	  }
	}
      }
	

      // first check if correct place in code(debug)

#if(DEBUGMPI>0)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d v[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      

#endif



    
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMBOUNDVEC;
    }// end over vectors
  }// endif vectors to be done
  
#endif // end if do bound vectors


#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after vector-2\n",numhit,myidtxt);
#endif


  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  // finally, wait on any sends left over before leaving
  // this is just same loops as above
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////

#if(DOBOUNDSCA)

  if(wsca!=0){
    for(l=1;l<=NUMSCA-NOBOUNDPOT;l++){
      /* if not to do all, pick */
      if(wsca!=-1){ if(wsca<=-2) ll=0; else ll=wsca; }
      else ll=l;
      // determine how things are bounded for other vars
      if(ll==0){
	if(wsca==-2){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-11){//other
	  // bound like rho
	  wbound=1;
	}
	else if(wsca==-12){//other
	  // bound like en
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0 for scalar case: wsca=%d\n",wsca);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*ll;
      // now that normal bcs are assigned, exchange with other cpus
      // tags are defined by sender's ID and direction sent
      // tag=(myid*COMPDIM*2)+{0,1,2,3,4,5}
      // 0=right, 1=up,2=left,3=down,4=out,5=in
      // works/v bc[1=output/2=input][0,1,2,3,4,5]
      // so sends are like: (sendtoid,myid*COMPDIM*2+?) and recv's are like: (fromid,otherid*COMPDIM*2+*) where ? and * are opposites(i.e. 0 and 2, 1 and 3)


      //////////////// periodicx3 transfer
      if(mpiperiodicx3){
	if(mycpupos[3]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  othercpu=myid+(ncpux3-1)*ncpux1*ncpux2;
	  tagrecv=tagheader+othercpu*COMPDIM*2+4;
	  tagsend=tagheader+myid*COMPDIM*2+5;

	  MPI_Wait(&srequest[tagsend],&mpichstatus);

	}
	else if(mycpupos[3]==ncpux3-1){
	  othercpu=myid-(ncpux3-1)*ncpux1*ncpux2; // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+4;
	  tagrecv=tagheader+othercpu*COMPDIM*2+5;

	  MPI_Wait(&srequest[tagsend],&mpichstatus);

	}// end if k+ boundary and not on k- boundary
      }// end if periodicx3



      

      if(srdir[5]){ // do in
	othercpu=myid-ncpux1*ncpux2;
	tagrecv=tagheader+othercpu*COMPDIM*2+4;
	tagsend=tagheader+myid*COMPDIM*2+5;

	MPI_Wait(&srequest[tagsend],&mpichstatus);

      }

      if(srdir[4]){ // do out
	othercpu=myid+ncpux1*ncpux2;
	tagsend=tagheader+myid*COMPDIM*2+4;
	tagrecv=tagheader+othercpu*COMPDIM*2+5;

	MPI_Wait(&srequest[tagsend],&mpichstatus);

      }
      
#if(DEBUGMPI>0)
      // first check if correct place in code(debug)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d s[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      
#endif

      /* cut short loop if only to do one */
      if(wsca!=-1) l=NUMSCA;
    } // end over scalars
  } // end if any scalars

#endif // endif do boundsca==1

#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after scalar-3\n",numhit,myidtxt);
#endif

  
#if(DOBOUNDVEC)  
  // Now do vectors if any
  if(wvec!=0){
    for(l=1;l<=REALNUMBOUNDVEC;l++){
   
      /* if not to do all, pick */
      if(wvec!=-1){ if(wvec<=-2) ll=0; else ll=wvec; }
      else ll=l;
      
      if(ll==0){
	if((wvec==-3)||(wvec==-4)){//mdot or flux(2) of v[1][2]
	  // bound like v
	  wbound=1;
	}
	else if((wvec==-5)){//mdot or flux(2) of v[1][2] or magemf
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-11){// v1
	  // bound like v
	  wbound=1;
	}
	else if(wvec==-12){// v2
	  // bound like B
	  wbound=2;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0, wvec=%d\n",wvec);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(ll!=0){
	tagheader=COMPDIM*2*numprocs*(NUMSCA+ll-1);
      }
      else tagheader=COMPDIM*2*numprocs*(NUMSCA+ll);



      /////////////////// periodicx3 transfer
      if(mpiperiodicx3){
	if(mycpupos[3]==0){
	  // now we have a cpu which needs to transfer to other boundary
	  othercpu=myid+(ncpux3-1)*ncpux1*ncpux2; // the periodic mate for x2
	  tagrecv=tagheader+othercpu*COMPDIM*2+4;
	  tagsend=tagheader+myid*COMPDIM*2+5;

	  MPI_Wait(&srequest[tagsend],&mpichstatus);

	}
	else if(mycpupos[3]==ncpux3-1){
	  othercpu=myid-(ncpux3-1)*ncpux1*ncpux2; // the periodic mate
	  tagsend=tagheader+myid*COMPDIM*2+4;
	  tagrecv=tagheader+othercpu*COMPDIM*2+5;

	  MPI_Wait(&srequest[tagsend],&mpichstatus);

	}// end if j+ and not k- boundary
      }// end if transfering periodicx3

      
      ///////////////////// exchange with other cpus, normal interior transfer

      if(srdir[5]){ // do in
        othercpu=myid-ncpux1*ncpux2;
	tagrecv=tagheader+othercpu*COMPDIM*2+4;
	tagsend=tagheader+myid*COMPDIM*2+5;

	MPI_Wait(&srequest[tagsend],&mpichstatus);

      }
      if(srdir[4]){ // do out
        othercpu=myid+ncpux1*ncpux2;
	tagsend=tagheader+myid*COMPDIM*2+4;
	tagrecv=tagheader+othercpu*COMPDIM*2+5;

	MPI_Wait(&srequest[tagsend],&mpichstatus);

      }

      // first check if correct place in code(debug)

#if(DEBUGMPI>0)

      MPI_Allreduce(&numhit, &numhitmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&numhit, &numhitmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if( (numhit!=numhitmin)||(numhit!=numhitmax)){
	fprintf(fail_file,"out of synch!\n");
	fprintf(fail_file,"proc: %d v[%d] numhit: %d\n",myid,ll,numhit);
	fflush(fail_file);
      }      

#endif



    
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMBOUNDVEC;
    }// end over vectors
  }// endif vectors to be done
  
#endif // end if do bound vectors

#if(DEBUGMPI>1)
  fprintf(stderr,"%d proc: %s got after vector-3\n",numhit,myidtxt);
#endif

#endif // endif compdim>=3






#endif// end if mpi



  firsttime=0;
}// end function
