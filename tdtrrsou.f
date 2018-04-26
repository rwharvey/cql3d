c
c
      subroutine tdtrrsou
      implicit integer (i-n), real*8 (a-h,o-z)

c.......................................................................
c     This routine computes the source term due to the radial transport
c     operator, R, applied on f_ (=f_n) needed with the ADI (Alternating
c     direction) method. R(f_) is calculated on the transport velocity
c     mesh and then is interpolated on the full velocity mesh, such
c     that int(R(f_)*d3u0) is constant.
c.......................................................................

      include 'param.h'
      include 'comm.h'
      character*8 nobind
      common/nob/ nobind
      include 'trans.h'

c.......................................................................
c     interpolate f_ on the transport velocity mesh, such that:
c     int(vptp*f_*d3u0) is constant
c.......................................................................

      call tdtrvtor2(f(0,0,1,1),frn(0,0,1,1),vpint,vpint_,1)

c.......................................................................
c     compute advective term
c.......................................................................

      call tdtravct(frn,kelecm,kelecg)  !No-op if pinch="disabled"

c.......................................................................
c     compute the source term
c.......................................................................

      do 10 k=1,ngen
        do 20 j=1,jx
          do 30 i=1,iytr(lrors)
            do 40 l=l_lower(i),lrors-1
C%OS  do 40 l=l_lower(i),lrors
              ilr=lrindx(l)
              ilrm1=lrindx(l-1)
              id=idx(i,l)
              ie=idx(i,l-1)
              if (l.ne.lpt(i).or. nobind.eq."enabled") then
c%OS  if (l.ne.lpt(i).or. nobind.eq."enabled".or.l.eq.lrors)then
                if (l.ne.1) then
c%OS  if (l.ne.1 .and. l.ne.lrors) then
                  spasou(id,j,k,l)=ztr(i,l)*zmaxpsi(ilr)*
     *              (h_r(ilr)  *bovcos(id,l)  *sfu(i,j,k,l)-
     -              h_r(ilrm1)*bovcos(ie,l-1)*sfu(i,j,k,l-1))
                else
c%OS  else if (l .eq. 1) then
                  spasou(id,j,k,l)=ztr(i,l)*zmaxpsi(ilr)*
     *              (h_r(ilr)  *bovcos(id,l)  *sfu(i,j,k,l))
C%OS  
C%OS  else
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
                spasou(id,j,k,l)=ztr(i,l)*zmaxpsi(ilr)*
     *            (h_r(ilr)  *bovcos(id,l)  *sfu(i,j,k,l)-
     -            h_r(ilrm1)*bovcos(ie,l-1)*sfu(i,j,k,l-1))
                spasou(id,j,k,l)=.5*spasou(id,j,k,l)+
     1            .5*ztr(i_,l)*zmaxpsi(ilr)*(h_r(ilr)*bovcos(i_d,l)
     1            *sfu(i_,j,k,l)-h_r(ilrm1)*bovcos(i_e,l-1)
     1            *sfu(i_,j,k,l-1))
              endif
 40         continue
 30       continue
 20     continue
 10   continue

c.......................................................................
c     interpolate source on itl-2,...,itl+2 points such that the integral
c     over velocity is conserved
c.......................................................................

      call tdtrrtov2(spasou,spasou,cynt2,cynt2_,3)

      return
      end
