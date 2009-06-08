cRAF
cRAF Eliminate the use of pointers by aliasing some arrays through
cRAF the subroutine call -- yucky, but portable!
cRAF
cRAF      subroutine iccgaf (km,lm,eps,ks,maxit,
cRAF     .                   a0,a1,b0,b1,bm1,x,y,work)
      subroutine iccgaf (km,lm,eps,ks,maxit,temp,a0save,d,c0,w,solr,
     .                   c1,cm1,r,p,a0,a1,b0,b1,bm1,x,y,work)
      implicit NONE
c-------------------------------------------------------------
c-- incomplete cholesky - conjugate gradient method
c-- written by alex friedman, llnl, 415-422-0827
c-- using algorithms of david kershaw, donald wolitzer.
 
c-- input arrays are all twice as long as needed to specify
c-- the problem; the second half is workspace.  that is, the
c-- input arrays should all have dimension at least
c-- 2*kmic*lmic. these arrays are: a0,a1,b0,b1,bm1,x,y.
 
c-- work is an additional array of length 4*kmic*lmic
c-- needed for working space.
 
c-- The matrix to be solved is symmetric, and has the structure:
c--
c--   a0
c--     1
c--
c--   a1   a0
c--     1    2
c--
c--        a1   a0
c--          2    3
c--
c--       :
c--       :
c--       :
c--       :
c--
c--                           a1    a0
c--                             km-2  km-1
c--
c--                                 a1    a0
c--                                   km-1  km
c--
c--   b0   bm1                                  a0
c--     1     1                                  km+1
c--
c--   b1    b0   bm1                            a1    a0
c--     1     2     2                             km+1  km+2
c--
c--         b1    b0   bm1                            a1    a0
c--           2     3     3                             km+2  km+3
c--
c--                :                                        :
c--                :                                        :
c--                :                                        :
c--
c-- Note that elements a1(sub)km etc. are zero due to the block structure.
 
c-- the correspondence with lasnex convention is:
c--   a0 <--> a
c--   a1 <--> b
c--   b0 <--> g
c--   b1 <--> d
c--   bm1<--> e
c-- note that the first physical element of bm1 is bm1(1), in
c-- contrast with the lasnex convention where it is bm1(2).
 
c-- scalar input variables are:
c-- km = kmic - the "short" dimension of the physical system
c--             (k is the rapidly varying index of the mesh array).
c-- lm = lmic - the "long" dimension of the physical system
c--             (l is the slowly varying index of the mesh array).
c--      eps  - convergence criterion, l2 normalized to y vector.
c--             on return, the actual accuracy achieved.
c--      ks   - last complete level of cyclic reduction desired
c--             in loops (one more level is done outside loops).
c--             the minimum possible value of ks is 0.
c--             the maximum possible value of ks is kp-1, where
c--             kp is the highest power of 2 with 2**kp .le. lmic.
c--             a reasonable choice is ks=4 for "most" problems.
c--             ks.le.14 at present due to dimension (15)
c--             in routine tssolve, allows lmic of 32k (big enough).
c--     maxit - the maximum number of c.g. iterations desired.
c--             on return, the number of iterations used.
 
c-- this package uses routines:
c--   TSDECOMP    GENCEES    ALTEVENS    GENBEES    MATMUL
c--   TSSOLVE     FORWARD    BACKWARD    ALTWS      MOVEXE
c--   FORMT1      TSDBY1     FORBY1      BACKBY1
c---------------------------------------------------------------
      integer kmic, lmic, kmic2x, numelts
      common /ciccg/ kmic, lmic, kmic2x, numelts
 
      integer km, lm, ks, maxit
      REAL eps
      REAL a0(*),a1(*),b0(*),b1(*),bm1(*),
     .          x(*),y(*),work(*)
      REAL a0save(*),d(*),c0(*),c1(*),cm1(*),w(*),
     .          solr(*),p(*),r(*),temp(*)
 
      integer i, kp, nevenb, inextlev, noddb, ibegino, ibegine, noddcb,
     .        nevencb, nbb, iter, kq, levlen
      REAL rerr, xerr, yabsum, dotprev, dotp, aa, dotr, b
      REAL sasum, sdot, snrm2
      external sasum, sdot, snrm2

cRAF      pointer (pa0s,a0save),(pd,d),(pc0,c0),(pw,w),(pc1,c1),
cRAF     .        (pcm1,cm1),(psolr,solr),(pp,p),(pr,r),(ptemp,temp)
c---- set common variables
      kmic = km
      lmic = lm
      kmic2x = 2 * kmic
      numelts = kmic * lmic
 
c---- set pointers (note that many arrays share storage).
cRAF      ptemp = loc (x(numelts+1))
cRAF      pa0s  = loc (y(numelts+1))
cRAF      pd    = loc (a0)
cRAF      pc0   = loc (work)
cRAF      pw    = loc (work)
cRAF      psolr = loc (work)
cRAF      pc1   = loc (work(numelts+1))
cRAF      pcm1  = loc (work(numelts+1))
cRAF      pr    = loc (work(2*numelts+1))
cRAF      pp    = loc (work(3*numelts+1))
 
c---- zero out work spaces
      do 10 i=numelts+1,numelts+numelts
        a0(i)=0.
        a1(i)=0.
        b0(i)=0.
        b1(i)=0.
        bm1(i)=0.
        x(i)=0.
 10   continue
      do 20 i=1,4*numelts
        work(i)=0.
 20   continue
 
c---- copy a0 into a0save to preserve it for matrix multiplies later
      call scopy(numelts,a0,1,a0save,1)
 
c---- compute kp (Kershaw's p, the maximum possible level of
c---- cyclic reduction).  this is used as a check on input ks.
      do 200 i = 1,1000
         if (2**i .gt. lmic) go to 210
  200 continue
  210 kp = i - 1
      ks = min0 (ks,kp-1)
 
c---- set parameters for first pass through decomposition loop
      nevenb = lmic
      inextlev = 1
 
c-----------------------------------------------------------------
c---- begin incomplete cholesky decomposition
c-----------------------------------------------------------------
 
c---- note that the c's are overwritten, level by level,
c---- and so the pointers to c0, c1, cm1 do not move.
 
      do 1000 kq = 0,ks
 
c---- set parameters for this level
      levlen = kmic * nevenb
      noddb = ( nevenb + 1 ) / 2
      nevenb = nevenb / 2
      ibegino = inextlev
      ibegine = ibegino + kmic
      inextlev = ibegino + levlen
 
c---- decompose upper left corner (calculate d's)
      call tsdecomp (a0(ibegino),a1(ibegino),d(ibegino),noddb)
 
c---- generate odd c's
      noddcb = nevenb
      call gencees (b0(ibegino),b1(ibegino),bm1(ibegino),
     .              a1(ibegino),d(ibegino),c0,cm1,noddcb)
 
c---- generate even c's
      nevencb = noddb - 1
      call gencees (b0(ibegine),bm1(ibegine),b1(ibegine),
     . a1(ibegino+kmic2x),d(ibegino+kmic2x),c0(1+kmic),c1(1+kmic),
     . nevencb)
 
c---- modify even diagonal arrays, lower right corner (calculate
c---- atilde's).  note that c1odd = b1odd, cm1even = bm1even.
      call altevens (a0(ibegine),a1(ibegine),c0,
     .               b1(ibegino),cm1,d(ibegino),
     .               c0(1+kmic),c1(1+kmic),bm1(ibegine),
     .               d(ibegino+kmic2x),
     .               a0(inextlev),a1(inextlev),noddcb,nevencb)
 
c---- calculate off-diagonal elements, lower right corner (btilde's)
      nbb = nevenb - 1
      call genbees (c0(1+kmic2x),b1(ibegino+kmic2x),cm1(1+kmic2x),
     .              c0(1+kmic),c1(1+kmic),bm1(ibegine),
     .              d(ibegino+kmic2x),
     .              b0(inextlev),b1(inextlev),bm1(inextlev),nbb)
 
 1000 continue
 
c---- do final level of tridiagonal sym. decomposition
      call tsdby1 (a0(inextlev),a1(inextlev),
     .               d(inextlev),nevenb)
 
c----------------------------------------------------------------
c---- end decomposition
c---- begin generalized conjugate gradient
c----------------------------------------------------------------
 
c---- form product A*x in work space w
      call matmul (a0save,a1,b0,b1,bm1,x,w,kmic,lmic)
 
c---- form residual r = y - A*x (this is r(nought))
      do 1200 i = 1,numelts
 1200 r(i) = y(i) - w(i)
 
c---- compute y norm for all relative error tests
      yabsum = sasum (numelts,y,1) 
 
c---- compute solr = (LLt)-1 on r
      call tssolve (solr,r,d,a1,b0,b1,bm1,solr,temp,ks)
 
c---- ... and set p(nought) to this.
c     call zmovewrd (p,solr,numelts)
      call scopy(numelts,solr,1,p,1)
 
c---- set up previous dotr product for iteration
      dotprev = sdot (numelts,r,1,solr,1)
 
c---------------------------------------------------------------
 
c---- begin conjugate gradient iteration loop
      do 2000 iter = 1,maxit
 
c------- compute m*p in temp (this is wolitzer's "evalp").  note -
c------- tssolve uses temp array, but we'll be done with it by then.
         call matmul (a0save,a1,b0,b1,bm1,p,temp,kmic,lmic)
 
c------- (p,mp) and exit if done.
         dotp = sdot (numelts,p,1,temp,1)
         if ( dotp.eq.0.0 ) go to 3000
 
c------- compute aa, the ratio of old dotr product and (p,mp).
         aa = dotprev / dotp
 
c------- x = x + aa * p
         call saxpy (numelts,aa,p,1,x,1)
 
c------- r = r - aa * m*p
	 aa = -aa
         call saxpy (numelts,aa,temp,1,r,1)
 
c------- compute error measure
         rerr = abs ( sasum(numelts,r,1) ) / yabsum
         xerr = abs(aa) * snrm2(numelts,p,1) / snrm2(numelts,x,1)
 
c------- test if done
         if (xerr.lt.eps.and.rerr.lt.eps) go to 3000
 
c------- compute new dotr product (r,solr) using solr = (LLt)-1 r
         call tssolve (solr,r,d,a1,b0,b1,bm1,solr,temp,ks)
         dotr = sdot (numelts,r,1,solr,1)
 
c------- b is the ratio of old to new dotr prods.  reset dotrprev.
         b = dotr / dotprev
         dotprev = dotr
 
c------- p = (llt)-1 r + b*p
         do 1600 i = 1,numelts
 1600       p(i) = solr(i) + b * p(i)
 
 2000 continue
 
c---- if this point reached, no convergence after maxit passes.
 
c----------------------------------------------------------------
c---- end generalized conjugate gradient
c----------------------------------------------------------------
 
 3000 continue
 
c---- set parameters to the values they must have on return
      eps = rerr
      maxit = iter
 
c---- restore a0 array for user convenience
c     call zmovewrd (a0,a0save,numelts)
      call scopy(numelts,a0save,1,a0,1)
      return
      end
      subroutine tsdecomp (a0,a1,d,nblocks)
      implicit NONE
      integer kmic, lmic, kmic2x, numelts
      common /ciccg/ kmic, lmic, kmic2x, numelts
 
      REAL a0(*),a1(*),d(*)
      integer nblocks, lastbl, i, j
 
c---- note that the d's overwrite the odd a0's at each level
c---- negative pivots are eliminated using absolute value.
c---- no check is made for small denominators.
 
c---- pointer to beginning of last block to be processed
      lastbl = ( nblocks - 1 ) * kmic2x + 1
 
c---- do first row of every other block (compute d's)
      do 100 j = 1,lastbl,kmic2x
  100 d(j) = abs ( 1. / a0(j) )
 
c---- do ith row of every other block (compute d's)
      do 200 i = 2,kmic
cdir$ ivdep
         do 200 j = 0,lastbl-1,kmic2x
  200    d(i+j) = abs(1./(a0(i+j)-a1(i+j-1)*a1(i+j-1)*d(i+j-1)))
      return
      end

      subroutine gencees (b0,b1,bm1,a1,d,c0,cm1,nblocks)
      implicit NONE
      integer kmic, lmic, kmic2x, numelts
      common /ciccg/ kmic, lmic, kmic2x, numelts
      REAL b0(*),b1(*),bm1(*),a1(*),d(*),c0(*),cm1(*)
      integer nblocks, lastbl, i, j
 
c---- test <<before>> entering loop
      if (nblocks.le.0) return
 
c---- pointer to beginning of last block to be processed
      lastbl = ( nblocks - 1 ) * kmic2x + 1
 
c---- do first row of every other block (compute c0's, c1's)
      do 100 j = 1,lastbl,kmic2x
      c0(j) = b0(j)
  100 cm1(j) = bm1(j) - a1(j)*d(j)*c0(j)
 
c---- do ith row of every other block (compute c0's, c1's)
      do 200 i = 2,kmic
         do 200 j = 0,lastbl-1,kmic2x
         c0(i+j) = b0(i+j) - a1(i+j-1)*d(i+j-1)*b1(i+j-1)
  200    cm1(i+j) = bm1(i+j) - a1(i+j)*d(i+j)*c0(i+j)
 
      return
      end
      subroutine altevens (a0k,a1k,c0km1,c1km1,cm1km1,dkm1,
     .              c0k,c1k,cm1k,dkp1,a0tild,a1tild,nbkm1,nbkkp1)
 
      implicit NONE
      integer kmic, lmic, kmic2x, numelts
      common /ciccg/ kmic, lmic, kmic2x, numelts
 
      REAL a0k(*),a1k(*),c0km1(*),c1km1(*),cm1km1(*),dkm1(*),
     .          c0k(*),c1k(*),cm1k(*),dkp1(*),a0tild(*),a1tild(*)
      integer nbkm1, nbkkp1, lastnewb, i, j, jj, lastaltb
 
c---- test <<before>> entering loop
      if (nbkm1.le.0) return
 
c---- pointer to beginning of last new block
      lastnewb = ( nbkm1-1 ) * kmic + 1
 
c----------------------------------------------------------------
c---- first set atilde to:  a - (k-1 term).
c----------------------------------------------------------------
 
c---- do first row of every other block
      j = 1 - kmic2x
      do 100 jj = 1,lastnewb,kmic
      j = j + kmic2x
      a0tild(jj) = a0k(j) - ( c0km1(j)*c0km1(j)*dkm1(j)
     .                        + cm1km1(j)*cm1km1(j)*dkm1(j+1) )
      a1tild(jj) = a1k(j) - ( c0km1(j)*c1km1(j)*dkm1(j)
     .                        + cm1km1(j)*c0km1(j+1)*dkm1(j+1) )
  100 continue
 
c---- do ith row of every other block
      do 200 i = 2,kmic
         j = -kmic2x
         do 200 jj = 0,lastnewb-1,kmic
         j = j + kmic2x
         a0tild(i+jj) = a0k(i+j) - ( c0km1(i+j)*c0km1(i+j)*dkm1(i+j)
     .                  + cm1km1(i+j)*cm1km1(i+j)*dkm1(i+j+1)
     .                  + c1km1(i+j-1)*c1km1(i+j-1)*dkm1(i+j-1) )
         a1tild(i+jj) = a1k(i+j) - ( c0km1(i+j)*c1km1(i+j)*dkm1(i+j)
     .                  + cm1km1(i+j)*c0km1(i+j+1)*dkm1(i+j+1) )
  200 continue
 
c------------------------------------------------------------------
c---- subtract off k/k+1 term to form true atilde.
c---- this is done in a separate loop to avoid "overreach"
c---- problems.  since the c pointers do not move, we can't
c---- count on the last block of c's being zero.
c------------------------------------------------------------------
 
c---- test <<before>> entering loop
      if (nbkkp1.le.0) return
 
c---- pointer to beginning of last new block to alter
      lastaltb = ( nbkkp1-1 ) * kmic + 1
 
c---- do first row of every other block
      j = 1 - kmic2x
      do 300 jj = 1,lastaltb,kmic
      j = j + kmic2x
      a0tild(jj) = a0tild(jj) - ( c0k(j)*c0k(j)*dkp1(j)
     .                            + c1k(j)*c1k(j)*dkp1(j+1) )
      a1tild(jj) = a1tild(jj) - ( c0k(j)*cm1k(j)*dkp1(j)
     .                            +c1k(j)*c0k(j+1)*dkp1(j+1) )
  300 continue
 
c---- do ith row of every other block
      do 400 i = 2,kmic
         j = -kmic2x
         do 400 jj = 0,lastaltb-1,kmic
         j = j + kmic2x
         a0tild(i+jj) = a0tild(i+jj) - ( c0k(i+j)*c0k(i+j)*dkp1(i+j)
     .                      + c1k(i+j)*c1k(i+j)*dkp1(i+j+1)
     .                      + cm1k(i+j-1)*cm1k(i+j-1)*dkp1(i+j-1) )
         a1tild(i+jj) = a1tild(i+jj) - ( c0k(i+j)*cm1k(i+j)*dkp1(i+j)
     .                      + c1k(i+j)*c0k(i+j+1)*dkp1(i+j+1) )
  400 continue
 
      return
      end

      subroutine genbees (c0kp1,c1kp1,cm1kp1,c0k,c1k,cm1k,dkp1,
     .                    b0tild,b1tild,bm1tild,nblocks)
 
      implicit NONE
      integer kmic, lmic, kmic2x, numelts
      common /ciccg/ kmic, lmic, kmic2x, numelts
 
      REAL c0kp1(*),c1kp1(*),cm1kp1(*),c0k(*),c1k(*),cm1k(*),
     .          dkp1(*),b0tild(*),b1tild(*),bm1tild(*)
      integer nblocks, lastnewb, i, j, jj
 
c---- test <<before>> entering loop
      if (nblocks.le.0) return
 
c---- pointer to beginning of last new block
      lastnewb = ( nblocks-1 ) * kmic + 1
 
c---- do first row of every other block
      j = 1 - kmic2x
      do 100 jj = 1,lastnewb,kmic
      j = j + kmic2x
      b0tild(jj) = - c0kp1(j)*c0k(j)*dkp1(j)
     .             - cm1kp1(j)*c1k(j)*dkp1(j+1)
      b1tild(jj) = - c1kp1(j)*c0k(j)*dkp1(j)
     .             - c0kp1(j+1)*c1k(j)*dkp1(j+1)
      bm1tild(jj) = - c0kp1(j)*cm1k(j)*dkp1(j)
     .              - cm1kp1(j)*c0k(j+1)*dkp1(j+1)
  100 continue
 
c---- do ith row of every block
      do 200 i = 2,kmic
         j = -kmic2x
         do 200 jj = 0,lastnewb-1,kmic
         j = j + kmic2x
         b0tild(i+jj) = - c0kp1(i+j)*c0k(i+j)*dkp1(i+j)
     .                - cm1kp1(i+j)*c1k(i+j)*dkp1(i+j+1)
     .                - c1kp1(i+j-1)*cm1k(i+j-1)*dkp1(i+j-1)
         b1tild(i+jj) = - c1kp1(i+j)*c0k(i+j)*dkp1(i+j)
     .                - c0kp1(i+j+1)*c1k(i+j)*dkp1(i+j+1)
         bm1tild(i+jj) = - c0kp1(i+j)*cm1k(i+j)*dkp1(i+j)
     .                 - cm1kp1(i+j)*c0k(i+j+1)*dkp1(i+j+1)
  200 continue
 
      return
      end

      subroutine matmul (a0,a1,b0,b1,bm1,x,w,kmic,lmic)
 
      implicit NONE
      integer kmic, lmic
      REAL a0(*),a1(*),b0(*),b1(*),bm1(*),x(*),w(*)
      integer i
 
c---- computes w = A*x.
 
      w(1) =                  a0(1)*x(1) + a1(1)*x(1+1)
      do 100 i = 2,kmic*lmic
  100 w(i) = a1(i-1)*x(i-1) + a0(i)*x(i) + a1(i)*x(i+1)
 
      w(1) = w(1)                        + b0(1)*x(kmic+1) 
     .            + b1(1)*x(kmic+1+1)
      do 200 i = 2,kmic*(lmic-1)
  200 w(i) = w(i) + bm1(i-1)*x(kmic+i-1) + b0(i)*x(kmic+i)
     .            + b1(i)*x(kmic+i+1)
 
      w(kmic+1) = w(kmic+1)                  + b0(1)*x(1) 
     .                      + bm1(1)*x(1+1)
      do 300 i = 2,kmic*(lmic-1)
  300 w(kmic+i) = w(kmic+i) + b1(i-1)*x(i-1) + b0(i)*x(i)
     .                      + bm1(i)*x(i+1)
 
      return
      end

      subroutine tssolve (x,y,d,a1,b0,b1,bm1,w,temp,ks)
 
      implicit NONE
      integer kmic, lmic, kmic2x, numelts
      common /ciccg/ kmic, lmic, kmic2x, numelts
 
      REAL x(*),y(*),d(*),a1(*),b0(*),b1(*),bm1(*),w(*),temp(*)
 
      integer ks, nevenb, inextlev, kq, levlen, noddb, ibegino, 
     .        ibegine, i, j
      integer nob(15),neb(15),ibo(15),ibe(15),inl(15)
 
c---- solves LDLt*x=y by (first) L*w=y (then) DLt*x=w
 
c---- note that w and x use same storage (solr).  they are
c---- kept distinct here for purposes of clarity.
 
c----------------------------------------------------------------
c---- begin forward sweep  L*w=y
c----------------------------------------------------------------
 
c---- set w at kq=0 level to y
c     call zmovewrd (w,y,numelts)
      call scopy(numelts,y,1,w,1)
 
c---- set parameters for first pass through forward sweep loop
      nevenb = lmic
      inextlev = 1
 
c---- begin loop
      do 1000 kq = 0,ks
 
c---- set parameters for this level
      levlen = kmic * nevenb
      noddb = ( nevenb + 1 ) / 2
      nevenb = nevenb / 2
      ibegino = inextlev
      ibegine = ibegino + kmic
      inextlev = ibegino + levlen
 
c---- save parameters for this level for future use
c---- in backward sweep (they are tricky to regenerate
c---- when moving backwards).
c---- kershaw, instead, generates the neven's using a circular
c---- shift right and mask, so the process is reversible.
      nob(kq+1) = noddb
      neb(kq+1) = nevenb
      ibo(kq+1) = ibegino
      ibe(kq+1) = ibegine
      inl(kq+1) = inextlev
 
c---- perform forward solve on odd blocks, L*w=w
      call forward (d(ibegino),a1(ibegino),w(ibegino),
     .              w(ibegino),noddb)
 
c---- generate (L)-t (D)-1 * w terms for odd blocks,
c---- i.e. DLt*temp=w
      call backward (d(ibegino),a1(ibegino),w(ibegino),temp,noddb)
 
c---- create next level of w's
      call altws (w(ibegine),b0(ibegino),b1(ibegino),
     .            bm1(ibegino),temp,
     .            b0(ibegine),b1(ibegine),bm1(ibegine),
     .            temp(kmic2x+1),
     .            w(inextlev),nevenb,noddb-1)
 
 1000 continue
 
c---- do kq=ks+1 forward solve
      call forby1 (d(inextlev),a1(inextlev),w(inextlev),
     .              w(inextlev),nevenb)
 
c-------------------------------------------------------------
c---- end forward sweep
c---- begin backward sweep
c-------------------------------------------------------------
 
c---- do kq=ks+1 backward solve
      call backby1 (d(inextlev),a1(inextlev),w(inextlev),
     .               x(inextlev),nevenb)
 
c---- begin loop
      do 2000 kq = ks,0,-1
 
c---- set parameters for this level
      noddb = nob(kq+1)
      nevenb = neb(kq+1)
      ibegino = ibo(kq+1)
      ibegine = ibe(kq+1)
      inextlev = inl(kq+1)
 
c---- move even blocks of x to lower level
      call movexe (x(ibegine),x(inextlev),nevenb)
 
c---- evaluate b(k)t*x(k+1) + b(k-1)*x(k-1) for odd k, store in temp
      call formt1 (b0(ibegino),b1(ibegino),bm1(ibegino),x(ibegine),
     .             b0(ibegine-kmic2x),b1(ibegine-kmic2x),
     .             bm1(ibegine-kmic2x),x(ibegine-kmic2x),
     .             temp,noddb)
 
c---- perform forward solve on this, keep in temp
      call forward (d(ibegino),a1(ibegino),temp,
     .              temp,noddb)
 
c---- ... and subtract from w, still in temp
      do 1200 i = 0,kmic-1
         do 1200 j = 0,(noddb-1)*kmic2x,kmic2x
 1200    temp(1+i+j) = w(ibegino+i+j) - temp(1+i+j)
 
c---- finally, perform backward solve into x
      call backward (d(ibegino),a1(ibegino),temp,
     .               x(ibegino),noddb)
 
 2000 continue
 
c------------------------------------------------------------
c---- end backward sweep
c------------------------------------------------------------
 
      return
      end

      subroutine forward (d,a1,y,w,nblocks)
 
      implicit NONE
      integer kmic, lmic, kmic2x, numelts
      common /ciccg/ kmic, lmic, kmic2x, numelts
 
      REAL d(*),a1(*),y(*),w(*)
      integer nblocks, lastbl, i, j
 
c---- solves L*w=y
 
c---- pointer to beginning of last block to be processed
      lastbl = ( nblocks-1 ) * kmic2x + 1
 
c---- do first row of every other block
      do 100 j = 1,lastbl,kmic2x
  100 w(j) = d(j) * y(j)
 
c---- do ith row of every other block
      do 200 i = 2,kmic
cdir$ ivdep
         do 200 j = 0,lastbl-1,kmic2x
  200    w(i+j) = d(i+j) * ( y(i+j) - a1(i+j-1)*w(i+j-1) )
 
      return
      end

      subroutine backward (d,a1,w,x,nblocks)
 
      implicit NONE
      integer kmic, lmic, kmic2x, numelts
      common /ciccg/ kmic, lmic, kmic2x, numelts
 
      REAL d(*),a1(*),w(*),x(*)
      integer nblocks, lastbl, i, j
 
c---- solves DLt*x=w
 
c---- pointer to beginning of last block to be processed
      lastbl = ( nblocks-1 ) * kmic2x + 1
 
c---- do last row of every other block
      do 100 j = kmic,lastbl+kmic-1,kmic2x
  100 x(j) = w(j)
 
c---- do ith row of every other block, progressing backwards
      do 200 i = kmic-1,1,-1
cdir$ ivdep
         do 200 j = 0,lastbl-1,kmic2x
  200    x(i+j) = w(i+j) - a1(i+j)*d(i+j)*x(i+j+1)
 
      return
      end

      subroutine altws (wk,b0km1,b1km1,bm1km1,tkm1,
     .                  b0k,b1k,bm1k,tkp1,wtilde,nblkslo,nblkshi)
 
      implicit NONE
      integer kmic, lmic, kmic2x, numelts
      common /ciccg/ kmic, lmic, kmic2x, numelts
 
      REAL wk(*),b0km1(*),b1km1(*),bm1km1(*),tkm1(*),
     .          b0k(*),b1k(*),bm1k(*),tkp1(*),wtilde(*)
      integer nblkslo, nblkshi, lastnblo, lastnbhi, i, j, jj
 
c---- compute next level of w's, the wtildes, for the solve
 
c---- pointers to beginning of last new block
      lastnblo = ( nblkslo-1 ) * kmic + 1
      lastnbhi = ( nblkshi-1 ) * kmic + 1
 
c---- do ith row of every other block

         j = -kmic2x
         do 10 jj = 0,lastnblo-1,kmic
         j = j + kmic2x
         wtilde(1+jj) = wk(1+j) - ( 
     .             b0km1(1+j)*tkm1(1+j) + bm1km1(1+j)*tkm1(1+1+j) )
  10     continue
         j = -kmic2x
         do 20 jj = 0,lastnbhi-1,kmic
         j = j + kmic2x
         wtilde(1+jj) = wtilde(1+jj) - ( 
     .             b0k(1+j)*tkp1(1+j) + b1k(1+j)*tkp1(1+1+j) )
  20  continue

      do 200 i = 2,kmic
         j = -kmic2x
         do 100 jj = 0,lastnblo-1,kmic
         j = j + kmic2x
         wtilde(i+jj) = wk(i+j) - ( b1km1(i+j-1)*tkm1(i+j-1)
     .           + b0km1(i+j)*tkm1(i+j) + bm1km1(i+j)*tkm1(i+j+1) )
  100    continue
         j = -kmic2x
         do 200 jj = 0,lastnbhi-1,kmic
         j = j + kmic2x
         wtilde(i+jj) = wtilde(i+jj) - ( bm1k(i+j-1)*tkp1(i+j-1)
     .           + b0k(i+j)*tkp1(i+j) + b1k(i+j)*tkp1(i+j+1) )
  200 continue
 
      return
      end

      subroutine movexe (xk,xtilde,nblocks)
 
      implicit NONE
      integer kmic, lmic, kmic2x, numelts
      common /ciccg/ kmic, lmic, kmic2x, numelts
 
      REAL xk(*),xtilde(*)
      integer nblocks, lastnewb, i, j, jj
 
c---- moves block k/2 of xtilde (=x at level kq+1) into block k
c---- of xk (=x at level kq), for k even
 
c---- pointer to beginning of last new block
      lastnewb = ( nblocks-1 ) * kmic + 1
 
c---- do ith row of every other block
      do 100 i = 1,kmic
         j = -kmic2x
         do 100 jj = 0,lastnewb-1,kmic
         j = j + kmic2x
  100    xk(i+j) = xtilde(i+jj)
 
      return
      end

      subroutine formt1 (b0k,b1k,bm1k,xkp1,
     .                   b0km1,b1km1,bm1km1,xkm1,temp,nblocks)
 
      implicit NONE
      integer kmic, lmic, kmic2x, numelts
      common /ciccg/ kmic, lmic, kmic2x, numelts
 
      REAL b0k(*),b1k(*),bm1k(*),xkp1(*),b0km1(*),b1km1(*),
     .          bm1km1(*),xkm1(*),temp(*)
      integer nblocks, lastbl, i, j
 
c---- compute intermediate term in backward solve
 
c---- pointer to beginning of last block
      lastbl = ( nblocks-1 ) * kmic2x + 1
 
c---- do part valid for all blocks, ith row of every other block,
c---- then do part not valid for first block
 
      do 10 j = 0,lastbl-1,kmic2x
  10     temp(1+j) =                           b0k(1+j)*xkp1(1+j)
     .               + b1k(1+j)*xkp1(1+1+j)
 
      do 20 j = kmic2x,lastbl-1,kmic2x
  20     temp(1+j) = temp(1+j) + b1km1(j)*xkm1(j)
     .         + b0km1(1+j)*xkm1(1+j) + bm1km1(1+j)*xkm1(1+1+j)
 
      do 300 i = 2,kmic
 
         do 100 j = 0,lastbl-1,kmic2x
  100    temp(i+j) = bm1k(i+j-1)*xkp1(i+j-1) + b0k(i+j)*xkp1(i+j)
     .               + b1k(i+j)*xkp1(i+j+1)
 
         do 200 j = kmic2x,lastbl-1,kmic2x
  200    temp(i+j) = temp(i+j) + b1km1(i+j-1)*xkm1(i+j-1)
     .         + b0km1(i+j)*xkm1(i+j) + bm1km1(i+j)*xkm1(i+j+1)
 
  300 continue
 
      return
      end

      subroutine tsdby1 (a0,a1,d,nblocks)
 
      implicit NONE
      integer kmic, lmic, kmic2x, numelts
      common /ciccg/ kmic, lmic, kmic2x, numelts
 
      REAL a0(*),a1(*),d(*)
      integer nblocks, lastbl, i, j
 
c---- this routine is similar to tsdecomp except that it
c---- does every block, not every other, as needed for
c---- incomplete cyclic reduction.
 
c---- pointer to beginning of last block to be processed
      lastbl = ( nblocks - 1 ) * kmic + 1
 
c---- do first row of every block (compute d's)
      do 100 j = 1,lastbl,kmic
  100 d(j) = abs ( 1. / a0(j) )
 
c---- do ith row of every block (compute d's)
      do 200 i = 2,kmic
cdir$ ivdep
         do 200 j = 0,lastbl-1,kmic
  200    d(i+j) = abs(1./(a0(i+j)-a1(i+j-1)*a1(i+j-1)*d(i+j-1)))
 
      return
      end

      subroutine forby1 (d,a1,y,w,nblocks)
 
      implicit NONE
      integer kmic, lmic, kmic2x, numelts
      common /ciccg/ kmic, lmic, kmic2x, numelts
 
      REAL d(*),a1(*),y(*),w(*)
      integer nblocks, lastbl, i, j
 
c---- similar to forward except that it does every block
 
c---- pointer to beginning of last block to be processed
      lastbl = ( nblocks-1 ) * kmic + 1
 
c---- do first row of every block
      do 100 j = 1,lastbl,kmic
  100 w(j) = d(j) * y(j)
 
c---- do ith row of every block
      do 200 i = 2,kmic
cdir$ ivdep
         do 200 j = 0,lastbl-1,kmic
  200    w(i+j) = d(i+j) * ( y(i+j) - a1(i+j-1)*w(i+j-1) )
 
      return
      end

      subroutine backby1 (d,a1,w,x,nblocks)
 
      implicit NONE
      integer kmic, lmic, kmic2x, numelts
      common /ciccg/ kmic, lmic, kmic2x, numelts
 
      REAL d(*),a1(*),w(*),x(*)
      integer nblocks, lastbl, i, j
 
c---- similar to backward except that it does every block
 
c---- pointer to beginning of last block to be processed
      lastbl = ( nblocks-1 ) * kmic + 1
 
c---- do last row of every block
      do 100 j = kmic,lastbl+kmic-1,kmic
  100 x(j) = w(j)
 
c---- do ith row of every block, progressing backwards
      do 200 i = kmic-1,1,-1
cdir$ ivdep
         do 200 j = 0,lastbl-1,kmic
  200    x(i+j) = w(i+j) - a1(i+j)*d(i+j)*x(i+j+1)
 
      return
      end
