#include "global.h"

#if(LINUXCLUSTER==0)
#include "defs.h"
#else
#include "decs.h"
#endif

#include "bound.h"


void bound1D_mpi(FTYPE (*vars)[N2M][N1M],FTYPE (*vars2)[N2M][N1M],int n, int row, int dir)
{
  int  i,j,k;
  int tagsend,tagrecv;
  int othercpu;
#if(USEMPI)
  MPI_Status status;

 if(dir==0){  // down
  if(row<ncpux2-1){
    for(i=-N1BND;i<N1+N1BND;i++){
      worksbc[0][1][3][i+N1BND]=vars[0][N2-2][i];
      worksbc[0][1][3][i+N1M+N1BND]=vars[0][N2-1][i];
	if(n==2){
	  worksbc[0][1][3][i+2*N1M+N1BND]=vars2[0][N2-2][i];
          worksbc[0][1][3][i+3*N1M+N1BND]=vars2[0][N2-1][i];
	}
    }
// send/recv data
    if(mycpupos[2]==row){
        othercpu=myid+ncpux1;
        tagsend=myid*COMPDIM*2+3;
        MPI_Send(worksbc[0][1][3],n*N2BND*N1M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD);
    }
    if(mycpupos[2]==row+1){
	othercpu=myid-ncpux1;
	tagrecv=othercpu*COMPDIM*2+1;
	MPI_Recv(worksbc[0][2][1],n*N2BND*N1M,MPI_FTYPE,othercpu,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
        for(i=-N1BND;i<N1+N1BND;i++){
          vars[0][-2][i]=worksbc[0][2][1][i+N1BND];
          vars[0][-1][i]=worksbc[0][2][1][i+N1M+N1BND];
	    if(n==2){
	       vars2[0][-2][i]=worksbc[0][2][1][i+2*N1M+N1BND];
               vars2[0][-1][i]=worksbc[0][2][1][i+3*N1M+N1BND];
	    }
	}
    }
  }
 }else{
  if(row>0){
    for(i=-N1BND;i<N1+N1BND;i++){
      worksbc[0][1][1][i+N1BND]=vars[0][0][i];
      worksbc[0][1][1][i+N1M+N1BND]=vars[0][1][i];
	if(n==2){
	  worksbc[0][1][1][i+2*N1M+N1BND]=vars2[0][0][i];
          worksbc[0][1][1][i+3*N1M+N1BND]=vars2[0][1][i];
	}
    }
// send/recv data
    if(mycpupos[2]==row){
        othercpu=myid-ncpux1;
        tagsend=myid*COMPDIM*2+1;
        MPI_Send(worksbc[0][1][1],n*N2BND*N1M,MPI_FTYPE,othercpu,tagsend,MPI_COMM_WORLD);
    }
    if(mycpupos[2]==row-1){
        othercpu=myid+ncpux1;
        tagrecv=othercpu*COMPDIM*2+3;
        MPI_Recv(worksbc[0][2][3],n*N2BND*N1M,MPI_FTYPE,othercpu,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
        for(i=-N1BND;i<N1+N1BND;i++){
          vars[0][N2][i]=worksbc[0][2][3][i+N1BND];
          vars[0][N2+1][i]=worksbc[0][2][3][i+N1M+N1BND];
	    if(n==2){
	      vars2[0][N2][i]=worksbc[0][2][3][i+2*N1M+N1BND];
              vars2[0][N2+1][i]=worksbc[0][2][3][i+3*N1M+N1BND];
	    }
	}
    }
  }
 }
#endif
}
