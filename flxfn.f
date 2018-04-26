c
c
      real*8 function flxfn(ra)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

cmnt  This routine calculates the poloidal magnetic field and the
cmnt  flux function at (normalized) radius ra (r/radmin) based on a
cmnt  gaussian current profile of specified width.

      include 'param.h'
      include 'comm.h'

      data idata/0/
cBH151202 Attempting to get code working with machine.eq.'mirror'
cBH151202      if (machine.ne."toroidal")call diagwrng(3)
      if (machine.ne."toroidal".and.machine.ne."mirror")call diagwrng(3)
      if (idata.eq.1) go to 10
      currwth=.5
      currnorm=1.-exp(-(1./currwth)**2)
      psi0=-radmin*radmaj*bth/currnorm*0.5
      wth2=1./currwth/currwth
      idata=1
 10   continue
      xra=ra/currwth
      xrbtemp=xra*xra
      if (xrbtemp.lt.1.e-13) then
        prf=xrbtemp/currnorm
      else
        expfact=exp(-xrbtemp)
        prf=(1.-expfact)/currnorm
      endif
      bthr(lr_)=0.
      atemp=sqrt(radmaj**2-(ra*radmin)**2)
      btemp=atemp/(radmaj**2-radmin**2)**.5
      if (ra.gt.0.) bthr(lr_)=btemp*bth*prf/ra
      robtemp=ra*ra
      wth2n=1.
      robtempn=1.
      denn=1.
      next=1
      last=1
      sum=0.
 20   continue
      wth2n=wth2n*wth2
      robtempn=robtempn*robtemp
      denn=-denn*next*next/last
      term=wth2n*(1.-robtempn)/denn
      sum=sum+term
      last=next
      next=next+1
      if (next.lt.100) go to 20
      flxfn=psi0*sum
      return
      end
