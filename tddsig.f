c
c
      real*8 function tddsig(k,p0,thet0,z)
      implicit integer (i-n), real*8 (a-h,o-z)
c
      real*8 k,l,k2
c
c     Bremsstrahlung differential cross-section in photon energy and angle.
c     Formula 2BN of Koch and Motz, Rev. of Mod Phys., vol. 31, no. 4, pp.
c     920-955(1959).
c     Units are cm**2 per atom per incident electron,  and
c     per photon energy (units of m_0c**2) and per steradian of
c     emitted  photon solid angle.
c
c     k= emitted photon energy in units of m0*c**2
c     p0= initial electron momentum in units of m0*c
c     thet0= radian angle between vector-p0 and momentum vector of photon.
c     z= ionic charge.
c
      tddsig=0.d0

      zero=0.d0
      one=1.d0
c      pi=3.141592653589793d0
      pi=atan2(zero,-one)
      r0=2.82d-13
      const=z*z*r0*r0/(8.*pi*137.)
c
      p02=p0*p0
      e02=p02+1.d0
c     e0 (e) are initial (final) total energies of the electron.
c     t0 is the initial kinetic energy.
      e0=sqrt(e02)
      t0=e0-1.d0
      if((k-t0).ge.0.) then
         tddsig=0.0
         return
      endif
      e=e0-k
      e2=e*e
      p2=e2-1.d0
      if(p2.le.zero) then
         tddsig=0.d0
         return
      endif
      p=sqrt(p2)
      k2=k*k
c990131      l=alog((e*e0-1.0+p*p0)/(e*e0-1.0-p*p0))
      l=log((e*e0-1.0+p*p0)/(e*e0-1.0-p*p0))
      cost=cos(thet0)
      sint=sin(thet0)
      del0=e0-p0*cost
      del02=del0*del0
      del04=del02*del02
c990131      ep_=alog((e+p)/(e-p))
      ep_=log((e+p)/(e-p))
      q2=p02+k2-2.*p0*k*cost
      q=sqrt(q2)
c990131      eq=alog((q+p)/(q-p))
      eq=log((q+p)/(q-p))
c
      tddsig=8.*sint*sint*(2.*e02+1.0)/(p02*del04)
     1  -2.*(5.*e02+2.*e*e0+3.)/(p02*del02)
     2  -2.*(p02-k2)/(q2*del02)+4.*e/(p02*del0)
      tddsig=tddsig
     1  +(l/(p*p0))*(4.*e0*sint*sint*(3.*k-p02*e)/(p02*del04)
     2  +4.*e02*(e02+e2)/(p02*del02)
     3  +(2.-2.*(7.*e02-3.*e*e0+e2))/(p02*del02)
     4  +2.*k*(e02+e*e0-1.)/(p02*del0))
      tddsig=tddsig
     1  -4.*ep_/(p*del0)+eq/(p*q)*(4./del02-6.*k/del0-2.*k*(p02-k2)/
     2  (q2*del0))
c
      tddsig=const/k*p/p0*tddsig
c
c     Elwert factor for lower energy electons (G. Elwert, Ann. Physik
c     34 (1939) 178).
c     Also, given in Koch and Motz, formula II-6.
c     "As a rough guide, the Elwert factor may be expected to give 
c     results that are accurate to about 10% for electron energies
c     below about 0.1 Mev."(Koch and Motz, p. 931.)

      beta0=p0/e0
      beta=p/e
      f_elwert=beta0*(1.-exp(-2.*pi*z/(137.*beta0)))/
     1  (beta*(1.-exp(-2.*pi*z/(137.*beta))))

      tddsig=f_elwert*tddsig

      return

      end
c
c
      real*8 function haugyp(ec,e,t)
      implicit integer (i-n), real*8 (a-h,o-z)

      real*8 mec2



!      function [seHC,seH,fee,ec,e,t] = haugyp(ec,e,t)
!
!	Differential electron-electron bremstrahlung cross-section 
!         (dsigma/dk.domega)
!	Haug model including Coulomb correction factor.
!
!	Input:
!
!	  - ec: kinetic energy of the incoming electron (mec2)
!               [n,p] if m=1 or [p,m] if n=1, or [n,m] if p=1.
!               (keV)
!	  - e: photon energy (mec2)  [n,p] if m=1 or [p,m] 
!              if n=1, or [n,m] if p=1.  
!               (keV)
!	  - t: angle between the direction of displacement of the 
!              incoming electron and the photon emitted by 
!              bremsstrahlung (radian) [n,p] if m=1
!	       or [p,m] if n=1, or [n,m] if p=1.   
!              (radians).
!
!	Output:
!
!	  - haugyp function value = seHC,
!           seHC: Haug + Coulomb factor bremstrahlung 
!                 cross-section (cm**2) [n,m]
!	  - seH: Haug bremstrahlung cross-section (cm**2) [n,m]
!	  - fee: Elwert correction factor [n,m]
!
!	Warning: Cross-section units : cm**2 but energies are 
!                                      in relativistic units. 
!                To get cross-sections in standard cm**2/keV units, 
!                seH or seHC must be divided
!	         by 511 keV.
!
!
!by Y.PEYSSON CEA-DRFC 15/05/1991 <peysson@fedv09.cad.cea.fr>
!revised for MatLab 4.1 (22/08/1994)
!revised for Fortran, Bob Harvey, 11/28/00 ==> scalar i/o.
!
!
      haugyp=0.d0

      zero=0.d0
      one=1.d0
      pi=atan2(zero,-one)
      c0= 2.99792458d10  !YuP[01-2011] was 3.0e+10 speed of light (cm/s)
      r0= 2.82d-13       !classical radius of the electron (cm).
      mec2= 510.998902d0 !YuP[01-2011] was 511. 
      !Rest mass energy of the electron (keV).
!
!
      if (ec .lt. 1.d-10) ec=ec+ 1.d-10
!
      seHC=0.d0      	
      seH=0.d0
      fee=1.d0
!
      if (ec .le. e) then
         haugyp=0.d0
         return
      endif
!
      ec = ec/mec2
      e = e/mec2
      e0 = ec+1.d0
      c = cos(t)
      s = sin(t)
      p0 = sqrt(e0**2-1.d0)
      v0 = p0/e0
!      p = sqrt((e0-e)**2-1)           !Not needed
!      ep = e0-e                       !Not needed
!
!Calculation parameters
!
      x1 = e*(e0-p0*c)
      x2 = e
      w2 = 2*(e0+1.d0)
      x = x1+x2
      rho2 = w2-2*x
      w24 = w2-4.d0
      rho24 = rho2-4.d0

!
!Photon only emitted if rho24.ge.0
      if (rho24 .le. zero) then
         haugyp=0.d0
         return
      endif
!
      w22 = w2-2.d0
      rho22 = rho2-2.d0
!      write(*,*) 'ec,e,e0,w2,x1,x2,x,rho2,rho24',
!     +            ec,e,e0,w2,x1,x2,x,rho2,rho24
!
!Coulomb correction factor
!
      betaw = sqrt(w2)*sqrt(w24)/w22
      betarho = sqrt(rho2)*sqrt(rho24)/rho22
      ap = 1.d0/(137*betaw)
      a2 = 1.d0/(137*betarho)
      fee = (a2/ap)*(exp(2*pi*ap)-1.d0)/(exp(2*pi*a2)-1.d0)
!
      r1 = rho24+4*x1+4*(x1**2)/rho2
      r2 = rho24+2*x1
      gw2 = sqrt(x2*(rho24*x2/4+2*x*x1/rho2))
      gw4 = sqrt(w24*(w24*rho24/4+4*x1*x2/rho2))
      xl = log((r2+sqrt(rho24*r1))*sqrt(rho2)/(4*x1))
      xl1 = log((sqrt(rho2)+sqrt(rho24))/2)
      xl2 = log(1+rho2*(rho24*x2+2*sqrt(rho24)*gw2)/(4*x*x1))
      xl3 = log(((sqrt(w2)*sqrt(rho24)+sqrt(rho2)*sqrt(w24))**2)/(8*x))
      xl4 = log(1+rho2*(w24*rho24+2*sqrt(rho24)*gw4)/(8*x1*x2))
      
!
!Main expression
!
      s0 = sqrt(rho24)
      s1 = (w2+rho2)*(((x1-x2)/x)**2)/(4*x1*x2)
      s2 = -0.25*((1d0/x1-1d0/x2)**2)
      s3 = -rho2/(2*(x**2))
      s4 = 2*(1+1d0/w24)*(rho2/(w24*x1*x2))
      s5 = 4*rho2/(w2*w24*x1*x1*x1*x1)
      s6 = 3*(w22**2)*x1*x2/w24
      s7 = -2*(x**2)*(1+6*1d0/(w2*w24))
      s8 = -4*rho2*x/(w24*x1*x1*x1)
      s9 = rho2*(4*1d0/w2
     +     -1.5-8*1d0/(w2*(w24**2))+w22*x/(w2*w24))/(x1**2)
      s10 = -rho2/(w24*x1)
      s11 = (w2-4*x2)/rho2
      s12 = -w2*w24/(4*x1*x2)
      s13 = -4*1d0/x1
      s14 = (w2-w22*(rho2/2))/(x1**2)
!
      z1 = (s1+s2+s3+s4+s5*(s6+s7)+s8+s9+s10+(s11+s12+s13+s14)/r1)*s0
!
      s0 = xl/sqrt(r1)
      s1 = 4+(10-w2/2)/x2
      s2 = -2*(2*w2-x2+4)/x1
      s3 = 3*w2*w24/(4*x1*x2)
      s4 = 4*1d0/(x1*r2)
      s5 = x2-3*x1+4+2*(rho2-3)/x1
      s6 = (w24-4*x1*x2/rho2)/(x1*r1)
      s7 = w2*r2/(4*x2)-rho22-2*x1
!
      z2 = s0*(s1+s2+s3+s4*s5+s6*s7)
!
      z3 = sqrt(rho2)*xl1*((rho2+2)/(x**2)+8*1d0/(x1**2))
!
      s0 = xl2/gw2
      s1 = 2*w22/x2
      s2 = -(rho22-x2)/x1
      s3 = rho22*(rho2+x1)/(2*x)
      s4 = w24/r2-2+rho22*((w24+rho2)**2)/(8*x*x1)
      s5 = 1d0/(r2*x1)
      s6 = rho22*(x2*(w24+rho2)/2-w22)-2*(x2**2)-4*x2
      s7 = 1d0/(r2*x)
      s8 = 0.5*rho22*(3*rho24-w2*(rho2-5))+2*(x1**2)-6*x1+w22*x2
!
      z4 = s0*(s1+s2+s3+s4+s5*s6+s7*s8)
!
      s0 = 2*sqrt(rho2)*xl3/(sqrt(w2)*sqrt(w24)*x1)
      s1 = 4+8*w22/(w24**2)
      s2 = -3*0.5*w22/x2
      s3 = ((x2**2)-2*w2-(w2-1)*x2+0.5*(w22**2))/x1
      s4 = (0.5*rho2*rho22/x)-(x2*w22*0.5/x)
      s5 = -((w22**2)*(rho2+w24))/(4*x*x2)
      s6 = x/(4*x1*x2)
      s7 = ((w22**2)+rho22*rho24)-8*rho22/w24
      s8 = 4*(x**2)/(w2*w24*x1)
      s9 = -(2*rho22*x2+w24*x-4*w22+w22*(8-rho2)/(2*x1))/r2
      s10 = (rho22*0.5*(3*rho24-w2*(rho2-5))-x1*(w2+4-2*x1))/(x*r2)
      s11 = 2*(w2*w22*rho24-2*rho22+4*rho2/w2)/((w24**2)*x1)
      s12 = 4*w22*x/(w24*(x1**2))
      s13 = (12*x/(w2*w24))-rho22/2
      s14 = 1-x/(w2*x1)
!
      z5 = s0*(s1+s2+s3+s4+s5+s6*s7+s8+s9+s10+s11+s12*s13*s14)
!
      s0 = (-1)*xl4/gw4
      s1 = 1
      s2 = rho22/(8*x1*x2)
      s3 = w22**2+(rho22**2)-6*(w24+rho2)+16*x/w24
      s4 = 2*(1-x1-(x1**2))/(x1*x2)
      s5 = (rho24-8*1d0/w24)/w24
!
      z6 = s0*(s1+s2*s3+s4+s5)
!
      z7a = z1+z2+z3+z4+z5+z6
!
      x0 = x1
      x1 = x2
      x2 = x0
!
      r1 = rho24+4*x1+4*(x1**2)/rho2
      r2 = rho24+2*x1
      gw2 = sqrt(x2*(rho24*x2/4+2*x*x1/rho2))
      gw4 = sqrt(w24*(w24*rho24/4+4*x1*x2/rho2))
      xl = log((r2+sqrt(rho24*r1))*sqrt(rho2)/(4*x1))
      xl1 = log((sqrt(rho2)+sqrt(rho24))/2)
      xl2 = log(1+rho2*(rho24*x2+2*sqrt(rho24)*gw2)/(4*x*x1))
      xl3 = log(((sqrt(w2)*sqrt(rho24)+sqrt(rho2)*sqrt(w24))**2)/(8*x))
      xl4 = log(1+rho2*(w24*rho24+2*sqrt(rho24)*gw4)/(8*x1*x2))

!
!Main expression...
!
      s0 = sqrt(rho24)
      s1 = (w2+rho2)*(((x1-x2)/x)**2)/(4*x1*x2)
      s2 = -0.25*((1d0/x1-1d0/x2)**2)
      s3 = -rho2/(2*(x**2))
      s4 = 2*(1+1d0/w24)*(rho2/(w24*x1*x2))
      s5 = 4*rho2/(w2*w24*x1*x1*x1*x1)
      s6 = 3*(w22**2)*x1*x2/w24
      s7 = -2*(x**2)*(1+6*1d0/(w2*w24))
      s8 = -4*rho2*x/(w24*x1*x1*x1)
      s9 = rho2*(4*1d0/w2
     +     -1.5-8*1d0/(w2*(w24**2))+w22*x/(w2*w24))/(x1**2)
      s10 = -rho2/(w24*x1)
      s11 = (w2-4*x2)/rho2
      s12 = -w2*w24/(4*x1*x2)
      s13 = -4*1d0/x1
      s14 = (w2-w22*(rho2/2))/(x1**2)
!
      z1 = (s1+s2+s3+s4+s5*(s6+s7)+s8+s9+s10+(s11+s12+s13+s14)/r1)*s0
!
      s0 = xl/sqrt(r1)
      s1 = 4+(10-w2/2)/x2
      s2 = -2*(2*w2-x2+4)/x1
      s3 = 3*w2*w24/(4*x1*x2)
      s4 = 4*1d0/(x1*r2)
      s5 = x2-3*x1+4+2*(rho2-3)/x1
      s6 = (w24-4*x1*x2/rho2)/(x1*r1)
      s7 = w2*r2/(4*x2)-rho22-2*x1
!
      z2 = s0*(s1+s2+s3+s4*s5+s6*s7)
!
      z3 = sqrt(rho2)*xl1*((rho2+2)/(x**2)+8*1d0/(x1**2))
!
      s0 = xl2/gw2
      s1 = 2*w22/x2
      s2 = -(rho22-x2)/x1
      s3 = rho22*(rho2+x1)/(2*x)
      s4 = w24/r2-2+rho22*((w24+rho2)**2)/(8*x*x1)
      s5 = 1d0/(r2*x1)
      s6 = rho22*(x2*(w24+rho2)/2-w22)-2*(x2**2)-4*x2
      s7 = 1d0/(r2*x)
      s8 = 0.5*rho22*(3*rho24-w2*(rho2-5))+2*(x1**2)-6*x1+w22*x2
!
      z4 = s0*(s1+s2+s3+s4+s5*s6+s7*s8)
!
      s0 = 2*sqrt(rho2)*xl3/(sqrt(w2)*sqrt(w24)*x1)
      s1 = 4+8*w22/(w24**2)
      s2 = -3*0.5*w22/x2
      s3 = ((x2**2)-2*w2-(w2-1)*x2+0.5*(w22**2))/x1
      s4 = (0.5*rho2*rho22/x)-(x2*w22*0.5/x)
      s5 = -((w22**2)*(rho2+w24))/(4*x*x2)
      s6 = x/(4*x1*x2)
      s7 = ((w22**2)+rho22*rho24)-8*rho22/w24
      s8 = 4*(x**2)/(w2*w24*x1)
      s9 = -(2*rho22*x2+w24*x-4*w22+w22*(8-rho2)/(2*x1))/r2
      s10 = (rho22*0.5*(3*rho24-w2*(rho2-5))-x1*(w2+4-2*x1))/(x*r2)
      s11 = 2*(w2*w22*rho24-2*rho22+4*rho2/w2)/((w24**2)*x1)
      s12 = 4*w22*x/(w24*(x1**2))
      s13 = (12*x/(w2*w24))-rho22/2
      s14 = 1-x/(w2*x1)
!
      z5 = s0*(s1+s2+s3+s4+s5+s6*s7+s8+s9+s10+s11+s12*s13*s14)
!
      s0 = (-1)*xl4/gw4
      s1 = 1
      s2 = rho22/(8*x1*x2)
      s3 = w22**2+(rho22**2)-6*(w24+rho2)+16*x/w24
      s4 = 2*(1-x1-(x1**2))/(x1*x2)
      s5 = (rho24-8*1d0/w24)/w24
!
      z6 = s0*(s1+s2*s3+s4+s5)
!
      z7b = z1+z2+z3+z4+z5+z6
!
      z8 = z7a+z7b
      seH = (r0**2)*e*z8/(pi*137*sqrt(w2)*sqrt(rho2)*sqrt(w24))
!      seHC = (r0**2)*e*z8*fee/(pi*137*sqrt(w2)*sqrt(rho2)*sqrt(w24))
      seHC = seH*fee
!
      haugyp=seHC

      return
      end
