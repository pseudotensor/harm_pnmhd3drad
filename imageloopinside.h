if(ndims==2){
if(TILEIMAGE==0){
  if(GENERALSLICE==0){
    if(PLOTAXIS==3){
      i=i+iteri;
      if(i==ei){
	i=si;
	j=j+iterj;
      }
      if(j==ej){
	i=si;
	j=sj;
	k=k+iterk;
      }
      if(k==ek){
	i=si;
	j=sj;
	k=sk;
	sliceloop++;
      }
      if(sliceloop==SLICENUMBER) break;
    }
    else if(PLOTAXIS==2){
      i=i+iteri;
      if(i==ei){
	i=si;
	k=k+iterk;
      }
      if(k==ek){
	i=si;
	k=sk;
	j=j+iterj;
      }
      if(j==ej){
	i=si;
	k=sk;
	j=sj;
	sliceloop++;
      }
      if(sliceloop==SLICENUMBER) break;
    }
    else if(PLOTAXIS==1){
      j=j+iterj;
      if(j==ej){
	j=sj;
	k=k+iterk;
      }
      if(k==ek){
	j=sj;
	k=sk;
	i=i+iteri;
      }
      if(i==ei){
	j=sj;
	k=sk;
	sliceloop++;
      }
      if(sliceloop==SLICENUMBER) break;
    }
    color=1;
    if(POSTPROC==0){
      if(i>N1-1){ i=N1-1; color=0;}
      if(i<0){ i=0; color=0;}
      if(j>N2-1){ j=N2-1; color=0;}
      if(j<0){ j=0; color=0;}
      if(k>N3-1){ k=N3-1; color=0;}
      if(k<0){ k=0; color=0;}
      // find region within accretor
#if(DO_ACCRETOR)
      if(inside_accretor(CENT,i,j,k)) color=0;
#endif

#if((BOUNDTYPE>1)&&(SHOWBOUNDARYZONES==0))
      if(bzmask[k][j][i]!=0){ color=0; } // so ignore boundary or irrelevant zones
#endif
    }
    if(notimageoutputloop&&(color==0)) continue; // otherwise color it nothing
    //  if(notimageoutputloop==0){ printf("%d %d   %d %d %d   %d\n",ii,jj,i,j,k,color); fflush(stdout);}
    //  printf("%d %d   %d %d %d   %d\n",ii,jj,i,j,k,color); fflush(stdout);
  }
  else{
    ii=ii+iterii;
    if(ii==eii){
      ii=sii;
      jj=jj+iterjj;
    }
    if(jj==ejj){
      ii=sii;
      jj=sjj;
      sliceloop++;
    }
    if(sliceloop==SLICENUMBER) break;
    i=(int)((FTYPE)io[sliceloop]+((FTYPE)ie[sliceloop]-(FTYPE)io[sliceloop])/((FTYPE)(imagen1-1))*(FTYPE)ii+((FTYPE)ic[sliceloop]-(FTYPE)ie[sliceloop])/((FTYPE)(imagen2-1))*(FTYPE)jj);
    j=(int)((FTYPE)jo[sliceloop]+((FTYPE)je[sliceloop]-(FTYPE)jo[sliceloop])/((FTYPE)(imagen1-1))*(FTYPE)ii+((FTYPE)jc[sliceloop]-(FTYPE)je[sliceloop])/((FTYPE)(imagen2-1))*(FTYPE)jj);
    k=(int)((FTYPE)ko[sliceloop]+((FTYPE)ke[sliceloop]-(FTYPE)ko[sliceloop])/((FTYPE)(imagen1-1))*(FTYPE)ii+((FTYPE)kc[sliceloop]-(FTYPE)ke[sliceloop])/((FTYPE)(imagen2-1))*(FTYPE)jj);
    color=1;
    if(POSTPROC==0){
      if(i>N1-1){ i=N1-1; color=0;}
      if(i<0){ i=0; color=0;}
      if(j>N2-1){ j=N2-1; color=0;}
      if(j<0){ j=0; color=0;}
      if(k>N3-1){ k=N3-1; color=0;}
      if(k<0){ k=0; color=0;}
#if(DO_ACCRETOR)
      // find region within accretor
      if(inside_accretor(CENT,i,j,k)) color=0;
#endif
    }
    if(notimageoutputloop&&(color==0)) continue; // otherwise color it nothing
    //  if(notimageoutputloop==0){ printf("%d %d   %d %d %d   %d\n",ii,jj,i,j,k,color); fflush(stdout);}
  }
}
else{
  ii=ii+iterii;
  if(ii==eii){
    ii=sii;
    jj=jj+iterjj;
  }
  if(jj==ejj){
    ii=sii;
    jj=sjj;
    sliceloop++;
  }
  if(sliceloop==SLICENUMBER) break;

  if(TILEDIR==3){
    i=ii%N[1];
    j=jj%N[2];
    k=(int)(ii/N[1])+ntile1*(int)(jj/N[2]);
  }
  else if(TILEDIR==2){
    i=jj%N[1];
    k=ii%N[3];
    j=(int)(ii/N[3])+ntile1*(int)(jj/N[1]);
  }
  else if(TILEDIR==1){
    j=ii%N[2];
    k=jj%N[3];
    i=(int)(ii/N[2])+ntile1*(int)(jj/N[3]);
  }
  color=1;
  if(POSTPROC==0){
    if(i>N1-1){ i=N1-1; color=0;}
    if(i<0){ i=0; color=0;}
    if(j>N2-1){ j=N2-1; color=0;}
    if(j<0){ j=0; color=0;}
    if(k>N3-1){ k=N3-1; color=0;}
    if(k<0){ k=0; color=0;}
#if(DO_ACCRETOR)
    // find region within accretor
    if(inside_accretor(CENT,i,j,k)) color=0;
#endif

#if((BOUNDTYPE>1)&&(SHOWBOUNDARYZONES==0))
    if(bzmask[k][j][i]!=0){ color=0; } // so ignore boundary or irrelevant zones
#endif
  }
  if(notimageoutputloop&&(color==0)) continue; // otherwise color it nothing
  //  if(notimageoutputloop==0){ printf("%d %d   %d %d %d   %d\n",ii,jj,i,j,k,color); fflush(stdout);}
  //  printf("%d %d   %d %d %d   %d\n",ii,jj,i,j,k,color); fflush(stdout);
}
}
else{ // basically same as PLOTAXIS=3 above
  color=1;
  i=i+iteri;
  if(i==ei){
    i=si;
    j=j+iterj;
  }
  if(j==ej){
    i=si;
    j=sj;
    k=k+iterk;
  }
  if(k==ek){
    i=si;
    j=sj;
    k=sk;
    sliceloop++;
  }
  if(sliceloop==SLICENUMBER) break;

#if(DO_ACCRETOR)
  // find region within accretor
  if(inside_accretor(CENT,i,j,k)) color=0;
#endif

#if((BOUNDTYPE>1)&&(SHOWBOUNDARYZONES==0))
  if(bzmask[k][j][i]!=0){ color=0; } // so ignore boundary or irrelevant zones
#endif
  if(notimageoutputloop&&(color==0)) continue; // otherwise color it nothing

}
