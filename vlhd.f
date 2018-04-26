c
c
      real*8 function vlhd(vll,vprp,thp,nmod)
      implicit integer (i-n), real*8 (a-h,o-z)

c     Determine the local phenomenological lower hybrid
c     diffusion coefficient as a function of local (in z)
c     v-parallel, v-perpendicular and poloidal angle.
c
c
c
      include 'param.h'
      include 'comm.h'


      if (vprprop .eq. "enabled") then
        xvpr=1.
      else
        xvpr=0.
      endif
      if (vdalp .lt. 1.e-8) then
         write(*,*)'vlhd: ***WARNING***  vdalp.lt.1.e-8, reset to 0.03'
         vdalp=.03
      endif
c
c     determine the vparallel interval at R=R0
c
      vpmax=vparmax(nmod)*clight
      vpmin=vparmin(nmod)*clight
c
c
c     determine the vparallel range of nonzero D at this
c     orbit point.
c
      rovr0=(radmaj+xvpr*radmin*rovera(lr_)*cos(thp))/radmaj
      vmin=vpmin*rovr0
      vmax=vpmax*rovr0
c
c     determine the vperp interval at R=R0
c
      vprmax=vprpmax(nmod)*clight
      vprmin=vprpmin(nmod)*clight
c
c
c     determine the vperp range of nonzero D at this
c     orbit point.
c
      vprmin=vprmin*rovr0
      vprmax=vprmax*rovr0

c     Determine the local diffusion coefficient.
c
      delv=(vmax-vmin)*vdalp

      if (vmin .gt. vll .or. vmax .lt. vll) then
        vlhd=0.

      elseif (vprmin.gt.vprp .or. vprmax.lt.vprp) then
        vlhd=0.

      elseif (vlhpolmn(nmod)/57.29577.gt.thp 
     +        .or. vlhpolmx(nmod)/57.29577.lt.thp) then
        vlhd=0.

      else
c       the nonzero diffusion coefficient..
c
        dene=reden(kelec,lr_)
        te=temp(kelec,lr_)*1.e3
c990131        xlog=24.-alog(sqrt(dene)/te)
        xlog=24.-log(sqrt(dene)/te)
        if (kelecg.eq.1) then
          taucoll=3.44074e5*te**1.5/(zeff(lr_)*dene*xlog)
        elseif (kiong(1).eq.1) then
          ti=temp(kiong(1),lr_)*1.e3
          taucoll=2.08507e7*ti**1.5
     +        /(zeff(lr_)*dene*bnumb(kiong(1))**2*xlog)
        else
          stop 'in vlhd'
        endif
        difus=dlndau(nmod)/taucoll*(vth(1,lr_)/vnorm)**2

cBH021118:  Karney uses following with vlh_karney=1.,
c           see, Karney and Fisch, Phys. Fl. 28, p. 116 (1985).        
        if (vlh_karney.ne.0.) then
           clight2=clight*clight
           uu=vll**2+vprp**2
           if (relativ.eq.'disabled') then
              uu=sqrt(uu)
           else
              uu=sqrt(uu/(1.-uu/clight2))
           endif
           factor=(1./(1.+uu/vth(1,lr_)))**vlh_karney
           difus=factor*difus
c           write(*,*)'vlhd: factor=',factor
        endif

        if (vmin.le.vll .and. vll.le. vmin+delv) then
          vlhd=difus*(1.-cos((vll-vmin)/delv*pi))*.5
        elseif (vmax-delv .le. vll .and. vll .le. vmax) then
          vlhd=difus*(1.-cos((vmax-vll)/delv*pi))*.5
        else
          vlhd=difus
        endif
      endif

      return
      end
