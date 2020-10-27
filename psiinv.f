
c
      real*8 function psiinv(yval,iupdown)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c---------------------------------------------------------------------
c
cmnt  function psiinv(psi)
cmnt  calculates the arclength s associated with given value of 
cmnt  yval=psi=B(s)/B0,  ie., performs the functional inversion of 
cmnt  psi(s,lr_):   thus psi(psiinv(psi))=psi or psiinv(psi(s))=s.
c
c     yval is a psi(=B(s)/B0) value.
c     Note: Assuming u=u0 (energy conservation), 
c           we have B(s)/B0= [sin(theta)/sin(theta0)]^2
c           For a trapped particle: B(s)/B0= [1/sin(theta0)]^2
c           at the bouncing point.
c
c     iupdown=1 (upper equilibrium) or -1 (lower), for eqsym="none".
c     Arclength is measured counter-clockwise from s=0 at the
c     minimum |B|-pt, in either case.
c
c---------------------------------------------------------------------
      include 'param.h'
      include 'comm.h'

      data toll/1.e-8/,nhalve/10/
      niter=0

      isign=+1
      if (eqsym.eq."none" .and. iupdown.eq.-1) isign=-1

      if (yval.lt.1. .or. yval.gt.psimx(lr_)) then
        if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
        write(*,*)'psiinv: lr_, yval, psimx(lr_)=', lr_, yval,psimx(lr_)
        endif
        call diagwrng(7)
      else if (yval.eq.1.) then
        psiinv=0.
        return
      else if (abs(psimx(lr_)-yval).lt.1.e-13) then
        psiinv=z_bmax(lr_)
        return
      else
        if (isign.eq.1) then
           sold=((yval-1.)/(psimx(lr_)-1.))*z_bmax(lr_)
        else
           sold=((yval-1.)/(psimx(lr_)-1.))*(zmax(lr_)-z_bmax(lr_))
        endif
      endif

      if (isign.eq.1) then
         sl=0.
         sr=zmax(lr_)
      else
         sl=zmax(lr_)
         sr=z_bmax(lr_)
      endif
 15   continue
      do 5 i=1,nhalve
        yo=psif(sold)
        if (yo.lt.yval) then
          sl=sold
        else
          sr=sold
        endif
        sold=(sr+sl)*0.5
    5 continue
      relerr=abs((sl-sr)/sr)
      if(relerr.lt.toll) then
        psiinv=sold
        return
      endif
c     Newton-Raphson interation:
 10   snew=sold-(psif(sold)-yval)/psifp(sold)  !psifp has sign=isign
      if (((snew.gt.sr).or.(snew.lt.sl)).and.isign.eq.1) then
        sold=(sl+sr)*0.5
        go to 15
      endif
      if (((snew.gt.sl).or.(snew.lt.sr)).and.isign.eq.-1) then
        sold=(sl+sr)*0.5
        go to 15
      endif
      relerr=abs((snew-sold)/sold)
      if (relerr.lt.toll) go to 30
      niter=niter+1
      if (niter.gt.100) call diagwrng(2)
 20   continue
      sold=snew
      go to 10
 30   psiinv=snew
      return
      end
