
    for(l=1;l<=REALNUMVEC;l++){
      
      /* if not to do all, pick */
      if(wvec!=-1){ if(wvec<=-2) ll=0; else ll=wvec; }
      else ll=l;
      
      if(ll==0){
	if((wvec==-3)||(wvec==-4)){//mdot or flux(2) of v[1][2]
	  // bound like v
	  wbound=1;
	}
	else{
	  fprintf(fail_file,"No definition for ll==0, wvec=%d\n",wvec);
	  myexit(1);
	}
      }
      else wbound=ll;

      if(docom[1]){
	LOOPBOUNDV1{
	  whichzone=temptempi;
	  bct=maskv1[k][j][i];
	  
	  if(bct>0){ //skip zone if computational zone or null zone
	    
	    switch(bct){
	      
	    case 1:
	    case 2:
	      /* reflect and AOS */
	      // order important for if/elseif's
	      
	      // find average value
	      ftemp=0;
	      for(m=0;m<bzs[whichzone][0];m++){
		
		ii=bzv1[whichzone][m*3+0+1];
		jj=bzv1[whichzone][m*3+1+1];
		kk=bzv1[whichzone][m*3+2+1];
		
		// some edge-"vx[0]" if reflectix1=1 type zones may not be set to 0 like should be, general way to get them?
		if(ii-1==i){
		  ftempv[1]+=-workv[ll][1][kk][jj][ii+1];
		  workv[ll][1][kk][jj][ii]=0; // should remain 0 anyways
		}
		else if(ii+1==i){
		  ftempv[1]+=0; // then not copy
		}
		else if(ii>i){
		  ftempv[1]+=-workv[ll][1][kk][jj][ii+1];
		}
		else if(ii<i){
		  ftempv[1]+=-workv[ll][1][kk][jj][ii+1];
		}
		else if(ii==i){
		  ftempv[1]+=workv[ll][1][kk][jj][ii];
		}
	      }
	      if(docom[2]){
		if(jj-1==j){
		  ftempv[2]+=-workv[ll][2][kk][jj+1][ii];
		  workv[ll][2][kk][jj][ii]=0; // should remain 0 anyways
		}
		else if(jj+1==j){
		  ftempv[2]+=0; // then not copy
		}
		else if(jj>j){
		  ftempv[2]+=-workv[ll][2][kk][jj+1][ii];
		}
		else if(jj<j){
		  ftempv[2]+=-workv[ll][2][kk][jj+1][ii];
		}
		else if(jj==j){
		  ftempv[2]+=workv[ll][2][kk][jj][ii];
		}
	      }
	      if(docom[3]){
		if(kk-1==k){
		  ftempv[3]+=-workv[ll][3][kk+1][jj][ii];
		  workv[ll][3][kk][jj][ii]=0; // should remain 0 anyways
		}
		else if(kk+1==k){
		  ftempv[3]+=0; // then not copy
		}
		else if(kk>k){
		  ftempv[3]+=-workv[ll][3][kk+1][jj][ii];
		}
		else if(kk<k){
		  ftempv[3]+=-workv[ll][3][kk+1][jj][ii];
		}
		else if(kk==k){
		  ftempv[3]+=workv[ll][3][kk][jj][ii];
		}
	      }
	    }// end over zones
	    ftempv[1]/=(FTYPE)(bzs[whichzone][0]);// holds average value now
	    ftempv[2]/=(FTYPE)(bzs[whichzone][0]);// holds average value now
	    ftempv[3]/=(FTYPE)(bzs[whichzone][0]);// holds average value now
	    // now assign value
	    workv[ll][1][k][j][i]=ftempv[1];
	    workv[ll][2][k][j][i]=ftempv[2];
	    workv[ll][3][k][j][i]=ftempv[3];
	    break;
	  case 3:
	    /* Fix/Time vary: Dirichlet */
	    if( (wvec>0)||(wvec==-1)){
	      for(m=1;m<=3;m++){
		if(docom[m]){
		  workv[ll][m][k][j][i]=vanal[ll][m][k][j][i];
		  // if on inner edge, also get element v[0](due to staggered grid)
		  if((m==1)&&(ii-1==i)){
		    workv[ll][m][k][j][i+1]=vanal[ll][m][k][j][i+1];
		  }
		  else if((m==2)&&(jj-1==j)){
		    workv[ll][m][k][j+1][i]=vanal[ll][m][k][j+1][i];
		  }
		  else if((m==3)&&(kk-1==k)){
		    workv[ll][m][k+1][j][i]=vanal[ll][m][k+1][j][i];
		  }
		}
	      }
	    }
	    else{
	      fprintf(fail_file,"Case 3 bound with wvec==-2 has no definition\n");
	      myexit(1);
	    }
	    break;
	  case 4:
	    if(( (wvec==-3)||(wvec==-4))) break;
	    /* outflow */
	    // deal with asymmetry in velocity components on grid w.r.t. inner/outer edges
	    // deals also with inflow checking
	    // order important for if/elseif's

	    //flux can't be outflowed like velocity(i.e. leads to constant density on boundaries).  Thus, only periodic and MPI transfer should assign mass flux on boundaries.  Thus, otherwise should not bound mass flux.  Let computation of mass flux(which is ok out to -1 & N) get new mass, etc.
	    // only overwrite mdot[-1,-2...] and mdot[N+1,N+2...] like terms for consistency

	    // find average value
	    ftempv[1]=ftempv[2]=ftempv[3]=0;
	    ftempv2[1]=ftempv2[2]=ftempv2[3]=0;
	    for(m=0;m<bzs[whichzone][0];m++){

	      ii=bzs[whichzone][m*3+0+1];
	      jj=bzs[whichzone][m*3+1+1];
	      kk=bzs[whichzone][m*3+2+1];

	      // GODMARK
	      //fprintf(stderr,"wvec: %d whichzone: %d m=%d %d %d %d %d %d %d\n",wvec,whichzone,m,k,kk,j,jj,i,ii);  fflush(stderr);


	      if(docom[1]){
		// must check to see if we are on a zone where velocity is a "vx[0]" like edge position
		// this will be the case if any computatational zone is referenced that has ii>i & jj-j==0 & kk-k==0
		if(ii>i){
		  // now treat that zone in same manner as normal boundary zone
		  // by symmetry, require offset for inner edge vectors, and must treat as normal ith zone although at i+1
		  // deal with both i+1 and i simultaneously since otherwise have to have special super-inner edge routine(redundant, but ok)
		  ftempv2[1]+=workv[ll][1][kk][jj][ii+1]; // "v[0]-like"
		}
		ftempv[1]+=workv[ll][1][kk][jj][ii]; // v[-1,-2,etc.] like (overwritten by above if inside boundary, forcing symmetry)
	      }
	      if(docom[2]){
		if(jj>j){
		  if( (ii-i==0)&&(jj-j==1)&&(kk-k==0) ){
		    if((wbound==1)&&(workv[ll][2][kk][jj+1][ii]>0)&&(INFLOWCHECKIX2)){
		      workv[ll][2][kk][jj][ii]=0.0; // halt the inflow onto the grid by assigning that "boundary" zone
		    }
		    else workv[ll][2][kk][jj][ii]=workv[ll][2][kk][jj+1][ii];
		  }
		  ftempv[2]+=workv[ll][2][kk][jj][ii];
		  // extra stuff above compared to below due to asymmetry in grid w.r.t. vectors on edge in that direction
		}
		else if(jj<j){
		  if((wbound==1)&&(workv[ll][2][kk][jj][ii]<0)&&(INFLOWCHECKOX2)){
		    ftempv[2]+=0.0; // halt the inflow onto the grid by assigning that boundary zone
		  }
		  else ftempv[2]+=workv[ll][2][kk][jj][ii];
		}
		else{
		  ftempv[2]+=workv[ll][2][kk][jj][ii];
		}
	      }
	      if(docom[3]){
		if(kk>k){
		  if( (ii-i==0)&&(jj-j==0)&&(kk-k==1) ){
		    if((wbound==1)&&(workv[ll][3][kk+1][jj][ii]>0)&&(INFLOWCHECKIX3)){
		      workv[ll][3][kk][jj][ii]=0.0; // halt the inflow onto the grid by assigning that "boundary" zone
		    }
		    else workv[ll][3][kk][jj][ii]=workv[ll][3][kk+1][jj][ii];
		  }
		  ftempv[3]+=workv[ll][3][kk][jj][ii];
		  // extra stuff above compared to below due to asymmetry in grid w.r.t. vectors on edge in that direction
		}
		else if(kk<k){
		  if((wbound==1)&&(workv[ll][3][kk][jj][ii]<0)&&(INFLOWCHECKOX3)){
		    ftempv[3]+=0.0; // halt the inflow onto the grid by assigning that boundary zone
		  }
		  else ftempv[3]+=workv[ll][3][kk][jj][ii];
		}
		else{
		  ftempv[3]+=workv[ll][3][kk][jj][ii];
		}
	      }// end 3-comp
	    } // end over zones
	    ftempv[1]/=(FTYPE)(bzs[whichzone][0]);// holds average value now
	    ftempv[2]/=(FTYPE)(bzs[whichzone][0]);// holds average value now
	    ftempv[3]/=(FTYPE)(bzs[whichzone][0]);// holds average value now
	    // now assign value
	    workv[ll][1][k][j][i]=ftempv[1];
	    workv[ll][2][k][j][i]=ftempv[2];
	    workv[ll][3][k][j][i]=ftempv[3];
	    break;
	  case 5:
	    // periodic	    
	    // find average value
	    ftempv[1]=ftempv[2]=ftempv[3]=0;
	    for(m=0;m<bzs[whichzone][0];m++){

	      ii=bzs[whichzone][m*3+0+1];
	      jj=bzs[whichzone][m*3+1+1];
	      kk=bzs[whichzone][m*3+2+1];

	      if(docom[1]){
		ftempv[1]+=workv[ll][1][kk][jj][ii];
	      }
	      if(docom[2]){
		ftempv[2]+=workv[ll][2][kk][jj][ii];
	      }
	      if(docom[3]){
		ftempv[3]+=workv[ll][3][kk][jj][ii];
	      }
	    }// end over zones
	    ftempv[1]/=(FTYPE)(bzs[whichzone][0]);// holds average value now
	    ftempv[2]/=(FTYPE)(bzs[whichzone][0]);// holds average value now
	    ftempv[3]/=(FTYPE)(bzs[whichzone][0]);// holds average value now
	    // now assign value
	    workv[ll][1][k][j][i]=ftempv[1];
	    workv[ll][2][k][j][i]=ftempv[2];
	    workv[ll][3][k][j][i]=ftempv[3];
	    break;
	  case 99:
	  case 98:
	    // no nothing
	    break;
	  default:
	    fprintf(fail_file,"bound.c: switch(bct) error, case given is: %d\n\r",bct);
	    myexit(1);
	  } // end switch
	} // end if bct>0 for vectors
	

      }// end LOOP over current vector
    
      if(wbound==2){
	LOOPBOUNDGEN{ // now fix for divb
	  whichzone=temptempi;
	  bct=bzmask[k][j][i];
	  
	  if(bct>0){ //skip zone if computational zone or null zone ( can get rid....)
	    
	    // for this, just find whether difference between ii and i is generally + or -, sufficient for what's wanted
	    iisum=jjsum=kksum=0;
	    for(m=0;m<bzs[whichzone][0];m++){

	      ii=bzs[whichzone][m*3+0+1];
	      jj=bzs[whichzone][m*3+1+1];
	      kk=bzs[whichzone][m*3+2+1];

	      iisum+=(ii-i);
	      jjsum+=(jj-j);
	      kksum+=(kk-k);
	    }
	    if(docom[1]){
	      if((j<=N2)&&(k<=N3)){ // only do if valid memory element exists for b1()
		if(iisum>0){
		  workv[ll][1][k][j][i]=b1(-1,workv[ll],k,j,i);
		}
		else workv[ll][1][k][j][i]=b1(1,workv[ll],k,j,i); // true for iisum<=0, 0 case could mean sum'ed to 0 or is 0.  If is 0, then really don't need to do this, but if summed to 0, then should, but doesn't hurt to do it in general (maybe--GODMARK)
		// otherwise taken care of by other components
	      }
	    }
	    if(docom[2]){
	      if((i<=N1)&&(k<=N3)){ // only do if valid memory element exists for b2()
		if(jjsum>0){
		  workv[ll][2][k][j][i]=b2(-1,workv[ll],k,j,i);
		}
		else workv[ll][2][k][j][i]=b2(1,workv[ll],k,j,i);
		// otherwise taken care of by other components
	      }
	    }
	    if(docom[3]){
	      if((i<=N1)&&(j<=N2)){ // only do if valid memory element exists for b3()
		if(kksum>0){
		  workv[ll][3][k][j][i]=b3(-1,workv[ll],k,j,i);
		}
		else workv[ll][3][k][j][i]=b3(1,workv[ll],k,j,i);
		// otherwise taken care of by other components
	      }
	    }
	  }
	}
      }// end divb fix
      
      /* cut short loop if only to do one */
      if(wvec!=-1) l=REALNUMVEC;
    }// end over vectors
  }// endif vectors to be done
  
#endif // end if do bound vectors

