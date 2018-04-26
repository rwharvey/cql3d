c
c
      subroutine tdtrsavf
      implicit integer (i-n), real*8 (a-h,o-z)

c..............................................................
c     This routine saves various versions of the distribution
c     function.
c     After transport the distribution is stored in frn.
c     After velocity advancement the distribution is stored
c     in fvn.
c     One time step back we have frn_1 and fvn_1 respectively.
c     We store in f either (frn_1+fvn_1+fvn*frn)*.25 or
c     (fvn+frn)*.5
c     f is used to compute collisional coefficients and to update
c     macroscopic quantities, such as temperature, current.
c     For purposes of the next velocity step we need an old version
c     of f, frn_2. We use either (frn+frn_1)*.5====>frn_2 or we use
c     frn=====>frn_2.
c     Next we update distributions: fvn===>fvn_1 and
c     frn====>frn_1
c     Then we copy fvn====>fvn_1
c..............................................................

      include 'param.h'
      include 'comm.h'
      

      do 10 k=1,ngen
        do 20 j=1,jx
          do 30 l=1,lrors
            if (relaxtsp.eq."enabled" .and. nonadi.eq.2 .and. n.ge.2) 
     +        then
              do 41 i=1,iy_(l)
                f(i,j,k,l)=(frn(i,j,k,l)+fvn(i,j,k,l)+frn_1(i,j,k,l)
     1            +fvn_1(i,j,k,l))*.25
                frn_2(i,j,k,l)=(frn(i,j,k,l)+frn_1(i,j,k,l))*.5
C%OS  f(i,j,k,l)=frn_2(i,j,k,l)
 41           continue
              go to 42
            endif
C%OS  
            if (relaxtsp.eq."enabled" .and. n.ge.2) then
              do 40 i=1,iy_(l)
C%OS  f(i,j,k,l)=(frn(i,j,k,l)+fvn(i,j,k,l)+frn_1(i,j,k,l)
C%OS  1            +fvn_1(i,j,k,l))*.25
                frn_2(i,j,k,l)=(frn(i,j,k,l)+frn_1(i,j,k,l))*.5
                f(i,j,k,l)=frn_2(i,j,k,l)
 40           continue
            else
              do 50 i=1,iy_(l)
c-out for now...f(i,j,k,l)=(frn(i,j,k,l)+fvn(i,j,k,l))*.5
                f(i,j,k,l)=frn(i,j,k,l)
                frn_2(i,j,k,l)=frn(i,j,k,l)
 50           continue
            endif
C%OS  
 42         continue
C%OS  
            do 60 i=1,iy_(l)
              frn_1(i,j,k,l)=frn(i,j,k,l)
              fvn_1(i,j,k,l)=fvn(i,j,k,l)
 60         continue
 30       continue
 20     continue
 10   continue
      return
      end
