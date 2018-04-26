c
c
      subroutine tdsetsxr(ien)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c......................................................................
c     This subroutine computes the time independent piece of the SXR
c     (BREMSSTRAHLUNG) diagnostic. It involves an integral over the
c     differential scattering cross sections for electron-ion and 
c     electron-electron collisions. 
c     Argument ien of tdsetsxr is the photon energy index, indicating
c     photon energy enk(ien) in units of mc**2.
c     The results of the integration will be put into  
c     sigsxr(1:jx,0:msxr,1:nen,2), where jx refers to the velocity
c     dimension, msxr the Legendre polynomial dimension, nen the
c     dimension of the photon energy spectrum, and the last 
c     argument is for e-i and e-e collisions.
c......................................................................
      include 'param.h'
      include 'comm.h'

c

c     Compute cos and sin arrays. Determine Legendre polynomials out to
c     order msxr.
c
      jval_(ien)=jx
      do 306 j=1,jx
        if (tcsgm1(j)*fions(kelec) .gt. en_(ien)) then
          jval_(ien)=j-1
          go to 307
        endif
 306  continue
 307  continue
      if (ien .gt. 1) go to 70
      temc1(1)=0.
      temc1(iy)=pi
      hy2=temc1(iy)/dfloat(iy-1)
      do 20 i=2,iy-1
        temc1(i)=temc1(i-1)+hy2
 20   continue
      do 30 i=1,iy
        temc2(i)=cos(temc1(i))
        temc3(i)=sin(temc1(i))
 30   continue
      do 40 i=2,iy-1
        temc3(i)=twopi*temc3(i)*(temc1(i+1)-temc1(i-1))*0.5
 40   continue
      temc3(1)=pi*temc1(2)**2*0.25
      temc3(iy)=temc3(1)
      do 60  i=1,iy
        pleg(0,i)=1.
        if (msxr .eq. 0) go to 60
        pleg(1,i)=temc2(i)
        do 50 m=2,msxr
          pleg(m,i)=((2*m-1)*temc2(i)*pleg(m-1,i)-(m-1)*pleg(m-2,i))/m
 50     continue
 60   continue
 70   continue
      call bcast(sigsxr(1,0,ien,1),zero,jx*(msxr+1))
      call bcast(sigsxr(1,0,ien,2),zero,jx*(msxr+1))

c.......................................................................
c     compute the integral.
c     Loop over Legendre coefficients.
c.......................................................................

c.......................................................................
c     For e-ion emission, use Z=1, and adjust in loop
c       using sigsxr for Zeff.ne.1, including Elwert factor
c.......................................................................

      do 11 m=0,msxr

c.......................................................................
c     Loop over theta (angle)
c.......................................................................

        do 12 i=1,iy

c.......................................................................
c     Loop over velocity
c.......................................................................

          do 13 j=jval_(ien),jx
            sigsxr(j,m,ien,1)=sigsxr(j,m,ien,1)+temc3(i)*pleg(m,i)*
     1        tddsig(enk(ien),x(j)*vnorm/clight,temc1(i),one)
     1        *x(j)*cint2(j)*vnorm/gamma(j)*one_
 13       continue
 12     continue
 11   continue

c.......................................................................
c     For e-e emission
c.......................................................................

      if (softxry.ne."e-ion") then
      do 14 m=0,msxr

c.......................................................................
c     Loop over theta (angle)
c.......................................................................

        do 15 i=1,iy

c.......................................................................
c     Loop over velocity
c.......................................................................

          do 16 j=jval_(ien),jx
            sigsxr(j,m,ien,2)=sigsxr(j,m,ien,2)+temc3(i)*pleg(m,i)*
     1      haugyp(gamm1(j)*restmkev,enk(ien)*restmkev,temc1(i))
     1        *x(j)*cint2(j)*vnorm/gamma(j)*one_
            !YuP[01-2011] replaced 512 by restmkev
 16      continue
 15   continue
 14   continue
      endif

      return
      end
