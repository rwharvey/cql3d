c
c
      subroutine sourcee
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     Subroutine sourcee is the controlling routine for the
c     analytic  SOURCE routines, all of which begin with "sou".
c     This is a facile model, it simply computes a source
c     profile (usually Gaussian in nature) with a specified
c     current. A more sophisticated model for ions utilizes NFREYA, the
c     Monte Carlo Beam deposition code. These are the "fr" routines with
c     frmod="enabled"
c     The "sou" routines are independent of the "fr" routines and vice-
c     versa.
c     Also, calls Knock On source modules, if specified.
c..................................................................

      save
      include 'param.h'
      include 'comm.h'


c..................................................................
c     return if source is not to be recomputed this time step
c..................................................................

      if (nsou.eq.1 .or. n.eq.0) then
        continue
      elseif (mod(n,nsou).eq.1 .and. n.ne.1) then
        continue
      else
        return
      endif

c..................................................................
c     Return if the number of sources per species (nso) = 0
c..................................................................

      if (nso.eq.0) then
        return
      endif

c..................................................................
c     Initialization...
c..................................................................

      if (n.eq.0) then

c..................................................................
c     Compute constants used to determine Gaussian source profiles.
c     If soucoord="cart" specify cartesian species parameters,
c     and if soucoord="polar" specify polar parameters.
c     If soucoord="disabled", no gaussian sources.
c..................................................................

      if (soucoord .ne. "disabled") then
        do 20 k=1,ngen
CDIR$ NEXTSCALAR
          do 10 m=1,nso
            if (soucoord.eq."cart") then
c990131              sxllm1(k,m,lr_)=sign(1.,sellm1(k,m))*
              sxllm1(k,m,lr_)=sign(one,sellm1(k,m))*
     1          sqrt(abs(sellm1(k,m))/fions(k))
              sxllm2(k,m,lr_)=sellm2(k,m)/fions(k)
              sxppm1(k,m,lr_)=sqrt(seppm1(k,m)/fions(k))
              sxppm2(k,m,lr_)=seppm2(k,m)/fions(k)
            else
              xem1(k,m,lr_)=sqrt(sem1(k,m)/fions(k))
              xem2(k,m,lr_)=sem2(k,m)/fions(k)
              cosm1(k,m,lr_)=cos(sthm1(k,m)*pi/180.)
              cosm2(k,m,lr_)=scm2(k,m)
            endif
            zm1(k,m,lr_)=zmax(lr_)*szm1(k,m)
            zm2(k,m,lr_)=(zmax(lr_)*szm2(k,m))**2
 10       continue
 20     continue
      endif  ! On soucoord

c..................................................................
c     sounor will contain normalization constants after the call
c     to sounorm
c..................................................................

        call bcast(sounor(1,1,1,lr_),one,lz*ngen*nsoa)
        if (soucoord .ne. "disabled") then
        call sounorm
        do 40 k=1,ngen
          do 30 m=1,nso
            if (nonso(k,m) .eq. 1) nonso(k,m)=0
 30       continue
 40     continue
        endif  ! On soucoord
      endif  ! On n.eq.0

c..................................................................
c     Initialize the source profile to zero.
c..................................................................

      call bcast(source(0,0,1,indxlr_),zero,iyjx2*ngen)

c..................................................................
c     xlncur will contain the source current (/cm**2/sec).
c     In general asor*zmaxpsi(lr_)=xlncur  (asor in units particles/cc)
c..................................................................

      call bcast(xlncur(1,lr_),zero,ngen)

c..................................................................
c     Determine Guassian+knock-on  source(i,j,k,indxlr_) array
c..................................................................

      call sourcef

      call sourceko
     
c..................................................................
c     define source profile uniquely at x=0.
c..................................................................

      call sourc0

c..................................................................
c     Compute the source power.
c..................................................................

      do 70 k=1,ngen
        call sourcpwr(k)
 70   continue
      return
      end
