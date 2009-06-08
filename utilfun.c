#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif






int myexit(int call_code)
{
  int i,l;
  int cleanfinish;

  fprintf(stderr,"proc: %s : Exiting cc: %d\n",myidtxt,call_code);
  if(call_code>=0){
    if(fail_file) fclose(fail_file);
    if(log_file) fclose(log_file);
    if(myid<=0){
      if(POSTPROC==0) { if(logfull_file) fclose(logfull_file);}
    }
    if(DOLOGSTEP){
      if(POSTPROC==0) { if(logstep_file) fclose(logstep_file);}
    }
    if(DOLOGPERF){
      if(POSTPROC==0) { if(logperf_file) fclose(logperf_file);}
    }
    if(DODTDIAG){
      if(POSTPROC==0) { if(logdt_file) fclose(logdt_file);}
    }
    if(DOFLOORDIAG>=1){
      if(POSTPROC==0) { if(logfl_file) fclose(logfl_file);}
    }
    if(DOSPDIAG>=1){
      if(POSTPROC==0) { if(logsp_file) fclose(logsp_file);}
    }
  }
  if(call_code>0){
    fprintf(stderr,"proc: %s : Failure.  Please check failure file: cc: %d\n",myidtxt,call_code);

    if(call_code==5){
      cleanfinish=1;
    }
    else{
      cleanfinish=0;
#if(USEMPI)
      // must abort since no clear to communicate to other cpus now
      MPI_Abort(MPI_COMM_WORLD,1);
#endif
    }
  }
  else cleanfinish=1;

  if(cleanfinish){
    fprintf(stderr,"Ending Computation on proc: %s, holding for other cpus\n",myidtxt);

#if(USEMPI)
    // finish up MPI
    MPI_Barrier(MPI_COMM_WORLD); // required!
    for(i=0;i<4;i++){
      //      MPI_Comm_free(&combound[i]); // messy since makes nonmember have NULL comm.  Should make as is, and if not member, then skip
      MPI_Group_free(&grprem[i]);
    } 
    MPI_Group_free(&MPI_GROUP_WORLD); 
    MPI_Finalize();
#endif
    if(POSTPROC==0){
      if((BOUNDTYPE==2)||(BOUNDTYPE==3)){
	// free up memory      
	for(l=0;l<NUMINDEX;l++){
	  if(!( (l==0)||((BOUNDTYPE==3)&&( ((l>=1)&&(l<=5))||((l>=20)&&(l<=22)) ) ) ) ){ 
	    free(indx[l]);
	  }
	  if( (BOUNDTYPE==3)&&(l>0) ){
	    // iindx is global memory currently
	  }
	}
	
	// free boundary zone exchange stuff
	for(l=6;l<=19;l++){
	  for(i=0;i<numiter[l];i++){
	    bzs[l][i]-=3;
	    free(bzs[l][i]);
	  }
	  free(bzs[l]);
	}
	// free boundary zone exchange stuff
	for(l=23;l<=28;l++){
	  for(i=0;i<numiter[l];i++){
	    free(bzs[l][i]);
	  }
	  free(bzs[l]);
	}
      }
    }


    if(myid<=0) fprintf(stderr,"Ended Computation on all processors\n");
  }

  fprintf(stderr,"END\n");
  fflush(stderr);
  exit(0);
  return(0);
}


void itoa(int x,char*p)
{
  int temp1;
  int i;
  int digits=0;
 
  temp1=x;
  if(temp1==0) digits=1;
  else{
    for(i=0;i<DIGILEN;i++){
      if(temp1==0){
        digits=i;
        break;
      }
      else temp1=temp1/10;
    }
  }
  if(digits==0){
    fprintf(fail_file,"problem with itoa function: x: %d p: %s\n",x,p);
    myexit(1);
  }
 
  for(i=0;i<digits;i++){
    temp1=x/10;
    temp1=x-temp1*10;
    p[digits-i-1]=ZER+temp1;
    x=x/10;
  }
  p[digits]='\0';
}                                       


int mysys(char*com1, char*com2){
  
   pid_t pid=-1;
   int trycnt=0;

   while(pid<0){
     trycnt++;
     pid=fork2(); /* forking new process */
     if(pid==0){ /* child do the command */
       execlp(com1,com1,com2,NULL);
       _exit(0);
     }
     else if(pid<0)/* error creating new process */
       {
	 fprintf(fail_file,"can't create new process for: %s %s\n",com1,com2);
	 fprintf(fail_file,"try#: %d\n",trycnt);
	 fprintf(fail_file,"errno: %d\n",errno);
	 //myexit(1);
       }
     //else if(waitpid(pid,&status,0)!=pid) /* wait till child terminates */     status=-1;
   }
   return(pid);
 }

int fork2(void)
{
  pid_t pid;
  int rc;
  int status;
  
  if (!(pid = fork()))
    {
      switch (fork())
	{
	case 0:  return 0;
	case -1: _exit(errno);    /* assumes all errnos are <256 */
	default: _exit(0);
	}
    }
  
  if (pid < 0 || waitpid(pid,&status,0) < 0)
    return -1;
  
  if (WIFEXITED(status))
    if (WEXITSTATUS(status) == 0)
      return 1;
    else
      errno = WEXITSTATUS(status);
  else
    errno = EINTR;  /* well, sort of :-) */
  
  return -1;
}

void ptraddr(int nstep){

  printf("proc: %d nstep: %07d: ",myid,nstep);

#if(TVDLF)
  printf("p: %x dqx: %x dqy: %x Fx: %x Fy: %x ph: %x\n",p,dqx,dqy,Fx,Fy,ph);
#endif

  printf("dx: %x x: %x s: %x v: %x g: %x dg: %x dvl: %x ",
	 &dx[1][1][-NBIGBND],
	 &x[1][1][-NBIGBND],
	 &s[1][-N3BND][-N2BND][-N1BND],
	 &v[1][1][-N3BND][-N2BND][-N1BND],
	 &g[1][1][-NBIGBND],
	 &dg[1][1][-NBIGBND],
	 &dvl[1][1][-NBIGBND]
	 );

  printf("DS: %x oarcl: %x OVOL: %x ",
	 &ds[1][1][-N3BND][-N2BND][-N1BND],
	 &oarcl[1][1][-N3BND][-N2BND][-N1BND],
	 &ovol[1][-N3BND][-N2BND][-N1BND]);
	 

  
  printf("bcs: %x bcv: %x ",&bcs[1][1][-N3BND][-N2BND][-N1BND],&bcv[1][1][-N3BND][-N2BND][-N1BND]);

  printf("sanal: %x vanal: %x ",
	 &sanal[1][-N3BND][-N2BND][-N1BND],
	 &vanal[1][1][-N3BND][-N2BND][-N1BND]);

  printf("floorvars: %x ",&floorvars[1][-N3BND][-N2BND][-N1BND]);
  printf("floorvarv: %x ",&floorvarv[1][1][-N3BND][-N2BND][-N1BND]);
  printf("floorvar0: %x ",&floorvar0[1][-N3BND][-N2BND][-N1BND]);

  printf("losss: %x ",&losss[1][1][0][0]);
  printf("lossv: %x ",&lossv[1][0][1][0][0]);
  printf("lossvisc: %x ",&lossvisc[1][1][0][0]);

  printf("nu_fact: %x nu_real: %x ",
	 &nu_fact[-N3BND][-N2BND][-N1BND],
	 &nu_real[-N3BND][-N2BND][-N1BND]);


  printf("avg2_3: %x avg1_31: %x avg1_32: %x ",
	 &avg2_3[1][-N2BND][-N1BND],
	 &avg1_31[1][-N2BND],
	 &avg1_32[1][-N1BND]);

	 
  printf("work1: %x work2: %x work3: %x work4: %x work5: %x work6: %x work7: %x work8: %x work9: %x work10: %x ",  
	 &work1[-N3BND][-N2BND][-N1BND],
	 &work2[-N3BND][-N2BND][-N1BND],
	 &work3[-N3BND][-N2BND][-N1BND],
	 &work4[-N3BND][-N2BND][-N1BND],
	 &work5[-N3BND][-N2BND][-N1BND],
	 &work6[-N3BND][-N2BND][-N1BND],
	 &work7[-N3BND][-N2BND][-N1BND],
	 &work8[-N3BND][-N2BND][-N1BND],
	 &work9[-N3BND][-N2BND][-N1BND],
	 &work10[-N3BND][-N2BND][-N1BND]);

  printf("workv1: %x workv2: %x workv3: %x workv4: %x workv5: %x ",
	 &workv1[1][-N3BND][-N2BND][-N1BND],
	 &workv2[1][-N3BND][-N2BND][-N1BND],
	 &workv3[1][-N3BND][-N2BND][-N1BND],
	 &workv4[1][-N3BND][-N2BND][-N1BND],
	 &workv5[1][-N3BND][-N2BND][-N1BND]);

  printf("sigma: %x rost: %x ",
	 &sigma[1][1][-N3BND][-N2BND][-N1BND],
	 &rost[1][1][-N3BND][-N2BND][-N1BND]);

  printf("workiq: %x workiqavg: %x workviq: %x ",
	 &workiq[1][0][0],
	 &workiqavg[1][0][0],
	 &workviq[1][1][0][0]);


  printf("\n");

  fflush(stdout);


}


	 

