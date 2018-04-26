c     ONETWO DIVERGENCE
c
c
c     ONETWO DIVERGENCE
      subroutine frsuppor(a,ni,nj,nk,nl,isupp,ifail)
      implicit integer (i-n), real*8 (a-h,o-z)
c-----------------------------------------------------------------------
c     determines crude 2d support of four-dimensional array
c     a(i,j,k,l)   array searched for support approximation
c     isupp(4,k,l) output.  with k,l fixed, contains (in order),
c     i1,i2,j1,j2, the indices of a corresponding to the
c     vertices of the smallest rectangle supporting a.
c-----------------------------------------------------------------------
c     ONETWO DIVERGENCE
      include 'param.h'
      dimension a(kix2,kjx2,ke,kb),isupp(4,ke,kb)

      zero=0.d0

      do 10 i=1,ni
        do 11 j=1,nj
          if(a(i,j,nk,nl).ne.zero)go to 12
 11     continue
 10   continue
      go to 999
 12   continue
      isupp(1,nk,nl)=i
      do 20 i=ni,1,-1
        do 21 j=1,nj
          if(a(i,j,nk,nl).ne.zero)go to 22
 21     continue
 20   continue
      go to 999
 22   continue
      isupp(2,nk,nl)=i
      do 30 j=1,nj
        do 31 i=1,ni
          if(a(i,j,nk,nl).ne.zero)go to 32
 31     continue
 30   continue
      go to 999
 32   continue
      isupp(3,nk,nl)=j
      do 40 j=nj,1,-1
        do 41 i=1,ni
          if(a(i,j,nk,nl).ne.zero)go to 42
 41     continue
 40   continue
      go to 999
 42   continue
      isupp(4,nk,nl)=j
      ifail=0
      return
 999  continue
      ifail=1
      return
      end
