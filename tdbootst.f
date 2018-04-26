c
c
      subroutine tdbootst
      implicit integer (i-n), real*8 (a-h,o-z)
      save
      include 'param.h'
      parameter(itlrza=3*lrza+1)
      include 'comm.h'

      real*8 nedrv,k13,k23
      dimension workk(itlrza),d2ne(lrza),d2te(lrza),d2ti(lrza),
     1  nedrv(lrza),tedrv(lrza),tidrv(lrza),tidst(lrza),tedst(lrza),
     1  a(lrza,3),dndst(lrza),bpo(lrza),avms(lrza),pe(lrza),xi13(lrza)
     1  ,ynue(lrza),ynui(lrza)
      data isel /0/
c..................................................................
c     This routine computes the bootstrap current using Hinton and
c     Hazeltine formulas (Rev Mod Physics 48, 239 (1976)).
c     This version currently is written for eqmod="enabled" -
c     the circular analogue will be added later. The formula for
c     the current is MKS. So all quantities will be put into MKS
c     for the evaluation of the current which will remain in
c     Amps/cm**2
c..................................................................

      if (eqmod.ne. "enabled") return

c..................................................................
c     Will use splines to determine derivatives of densities and
c     temperatures vs. rho (=rz)
c..................................................................

      i1p(1)=4
      i1p(2)=4
      if (bootst.ne."enabled") then
        do 5 ll=1,lrzmax
          bscurm(ll,1,1)=0.
 5      continue
        return
      endif
      do 10 l=1,lrzmax
        tedst(l)=energy(kelec,l)*2./3.*1.e+3
        dndst(l)=reden(kelec,l)*1.e+6
        tr1(l)=rz(l)*.01
 10   continue
      call coeff1(lrzmax,tr1(1),dndst,d2ne,i1p,1,workk)
      call coeff1(lrzmax,tr1(1),tedst,d2te,i1p,1,workk)
      itab(1)=0
      itab(2)=1
      itab(3)=0
      do 11 l=1,lrzmax

c..................................................................
c     derivative of density w.r.t. rho
c..................................................................

        call terp1(lrzmax,tr1(1),dndst,d2ne,tr1(l),1,tab,itab)
        nedrv(l)=tab(2)

c..................................................................
c     derivative of te with respect to rho
c..................................................................

        call terp1(lrzmax,tr1(1),tedst,d2te,tr1(l),1,tab,itab)
        tedrv(l)=tab(2)
 11   continue
      call bcast(tr(1),zero,lrzmax)
      call bcast(tr2(1),zero,lrzmax)
      call bcast(tr3(1),zero,lrzmax)

c..................................................................
c     Determine average ion temperature over all background species
c..................................................................

      kk=0
      do 21 k=ngen+1,ntotal
        if (k.eq.kelecm .or. k.eq.kelecg) go to 21
        kk=kk+1
        do 20 l=1,lrzmax
          tr(l)=tr(l)+reden(k,l)*1.e+6
          tr2(l)=tr2(l)+reden(k,l)*energy(k,l)*2./3.*1.e+3*1.e+6
          tr3(l)=tr3(l)+reden(k,l)*1.e+6*fmass(k)/proton
 20     continue
 21   continue
      do 22 l=1,lrzmax
        tidst(l)=tr2(l)/tr(l)
        avms(l)=tr3(l)/tr(l)
 22   continue
      call coeff1(lrzmax,tr1(1),tidst,d2ti,i1p,1,workk)
      do 23 l=1,lrzmax
        call terp1(lrzmax,tr1(1),tidst,d2ti,tr1(l),1,tab,itab)
        tidrv(l)=tab(2)
        pe(l)=1.e-3*ergtkev*1.e-1*tedst(l)*dndst(l)*1.e-6
        bpo(l)=-1.e-4/rmag*dpsidrho(l)
 23   continue
c..................................................................
c     compute transport coefficients for current in banana-plateau
c     regime with single ion of charge zeff(lr_) with Hinton-Hazeltine
c     formula
c..................................................................
      twosqt=1.41421356
CDIR$ NEXTSCALAR
      do 200 l=1,lrzmax
        epsqrt=sqrt(aspin(l))
        dene=dndst(l)
        te=tedst(l)
        ti=tidst(l)

c..................................................................
c     xi13(l) is i13/1.95*sqrt(aspin(l)) of Hinton & Hazeltine
c     use Kim, Callen & Hamnen fit to fraction of passing electrons in
c     elliptic cross section for time being.
c..................................................................


        if (isel .eq. 1) then
          fc=1.-epsqrt*(1.46-0.46*aspin(l))
          xi13(l)=(1.+0.31*aspin(l))/fc
        else
          xi13(l)=1.0
        endif
        yi13=xi13(l)
c990131        xloge=31.26+alog(te/sqrt(dene))
c990131        xlogi=35.00+0.5*alog(te*ti/dene)
        xloge=31.26+log(te/sqrt(dene))
        xlogi=35.00+0.5*log(te*ti/dene)
        taue=3.44074d11*te**1.5/(zeff(l)*dene*xloge)
        taui=2.08507d13*ti**1.5/(zeff4(l)*1.e+6*xlogi)
        vthe(lr_)=5.931d5*sqrt(te)   !Evidently, sqrt(2Te/me) m/sec
        vthi=1.384d4*sqrt(ti/avms(l))
        xnue=1.e-6*twosqt*rz(l)*btor
     1    /(aspin(l)*epsqrt*bpo(l)*vthe(lr_)*taue)
        xnui=1.e-6*twosqt*rz(l)*btor/(aspin(l)*epsqrt*bpo(l)*vthi*taui)
        ynue(l)=xnue
        ynui(l)=xnui

c..................................................................
c.... In the next formula, H-H list 0.67, but I think this could be 0.57
c.... to agree with their Table III, or is Table III incorrect?
c     We attribute it to the 1,1 compontent of bscurm below:
c     bscurm(l,1,1)
c..................................................................

        k13=yi13*1.46*(1.0+0.670/zeff(l))/(1.0+1.02*sqrt(xnue)+
     1    1.07*xnue)
        k23=yi13*3.65*(1.0+0.148/zeff(l))/(1.0+0.57*sqrt(xnue)+
     1    0.61*xnue)
        bt1g2i=(1.17-0.35*sqrt(xnui))/(1.0+0.7*sqrt(xnui))
        a(l,1)=k13*(1.0+ti/(zeff(l)*te))
        a(l,2)=-k13*(1.0-bt1g2i/(1.0+xnue**2*aspin(l)**3))/zeff(l)
        a(l,3)=k23-1.5*k13
        bscurm(l,1,1)=-(sqrt(aspin(l))*pe(l)/bpo(l))*(a(l,1)*nedrv(l)
     *    /dndst(l)-a(l,2)*tidrv(l)/tedst(l)+a(l,3)*tedrv(l)/tedst(l))
     *    /1.e+4
 200  continue
      return
      end
