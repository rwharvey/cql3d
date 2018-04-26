c
c
      subroutine tdtrchk
      implicit integer (i-n), real*8 (a-h,o-z)

c..............................................................
c     This routine checks the accuracy of the radial time advancement
c     split by evaluating the rhs and the lhs of the difference
c     equation and then comparing. Aggregate sums (sumleft and
c     sumright) are also computed. Errors should be near round-off.
c..............................................................

      include 'param.h'
      include 'comm.h'

      character*8 nobind
      common/nob/ nobind
      include 'trans.h'
      if (iactst.eq."disabled") return
      call bcast(fxsp,zero,iyjx2*ngen*lrors)
      call bcast(f,zero,iyjx2*ngen*lrors)
      sumleft=0.
      sumright=0.
      do 10 k=1,ngen
        do 20 j=1,jx
          do 30 i=1,iytr(lrors)
            do 40 l=l_lower(i),lrors-1
C%OS  do 40 l=l_lower(i),lrors
              ilr=lrindx(l)
              ilrm1=lrindx(l-1)

c..............................................................
c     rhs first
c     RECALL!!!!  We have forced cynt2_(id,l)*cosovb(id,l)=
c     cynt2_(ie,l-1)*cosovb(ie,l-1). This is crucial
c     for conservation and relates back to the
c     Jacobian relating e,mu space to v,theta space. If
c     the code is running in the mode meshy="fixed_mu"
c     then cynt2*cosovb stays constant for for fixed
c     i as we move in l just as de*dmu stays constant.
c..............................................................

              id=idx(i,l)
              ie=idx(i,l-1)
              f(id,j,k,l)=(frn(id,j,k,l)-frn_2(id,j,k,l))/dttr*
     1          vpint_(id,ilr)*cint2(j)/zmaxpsi(ilr)*dvol(ilr)
              if (l.ne.lpt(i).or. nobind.eq."enabled") then
c%OS  if (l.ne.lpt(i).or. nobind.eq."enabled".or.l.eq.lrors)then
c%OS  if (l.ne.1) then
                if (l.ne.1 .and. l.ne.lrors) then
                  fxsp(id,j,k,l)=cosovb(id,l)*(h_r(ilr)*bovcos(id,l)
     1              *sfu(i,j,k,l)-h_r(ilrm1)*bovcos(ie,l-1)
     +              *sfu(i,j,k,l-1))
     1              *cynt2_(id,l)*cint2(j)*4.*pi**2*radmaj
                else
c%OS  else if (l .eq. 1) then
                  fxsp(id,j,k,l)=cosovb(id,l)*(h_r(ilr)*bovcos(id,l)
     1              *sfu(i,j,k,l))*cynt2_(id,l)*cint2(j)*4.*pi**2*radmaj
C%OS  
c%OS  else
c%OS  fxsp(id,j,k,l)=-cosovb(id,l)*h_r(ilrm1)*bovcos(ie,l-1)
c%OS  +          *sfu(i,j,k,l-1)*cynt2_(id,l)*cint2(j)*4.*pi**2*radmaj
C%OS  
                endif
              else if (lpt(i).ne.lrors .and. l.eq.lpt(i)) then
                id=idx(i,l)
                ie=idx(i,l-1)
                i_=iytr(lrors)+1-i
                i_d=idx(i_,l)
                i_e=idx(i_,l-1)
                fxsp(id,j,k,l)=cosovb(id,l)*(h_r(l)*bovcos(id,l)
     1            *sfu(i,j,k,l)-h_r(ilrm1)*bovcos(ie,l-1)
     +            *sfu(i,j,k,l-1))*cynt2_(id,l)*cint2(j)
                fxsp(id,j,k,l)=(.5*fxsp(id,j,k,l)+
     1            .5*cosovb(i_d,l)*(h_r(ilr)*bovcos(i_d,l)
     1            *sfu(i_,j,k,l)-h_r(ilrm1)*bovcos(i_e,l-1)
     1            *sfu(i_,j,k,l-1))*cynt2_(i_d,l)*cint2(j))*4.*pi**2
     1            *radmaj
              endif
              sumleft=sumleft+f(id,j,k,l)
              sumright=sumright+fxsp(id,j,k,l)
     +          + velsou(id,j,k,l)*cynt2_(id,l)*cint2(j)/zmaxpsi(ilr)
 40         continue
 30       continue
 20     continue
        error=2*abs(sumleft-sumright)/(sumleft+sumright)
        if (iactst.eq."abort") then
          if (error.gt.1.e-8) call diagwrng(15)
        endif
 10   continue
      return
      end
