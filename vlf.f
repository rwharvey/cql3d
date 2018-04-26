c
c
      subroutine vlf(action)
      implicit integer (i-n), real*8 (a-h,o-z)
      character*(*) action
      save

c  *** NOTE: Bunch of write(*,*) statements related to checking
c            out effect of taunew.  For now, use taunew="disabled".
c            See comments in a_change.h.

c.................................................................
c     subroutine vlf is a  cyclotron QL model
c     designed for electron/ion  heating and current
c     drive, according to the relativistic quasilinear operator.
c     (Stix(1993), p. 498, Eq. 41).
c    
c     It is primarily adapted from urf... routines: urfbes,
c     urfpack, and urfb0.
c
cBH060314     Presently, several harmonics at a time are treated,
cBH060314     or several modes (but not both simultaneously), but
c     Presently, several harmonics at a time are treated,
c     for several wave types/modes, but
c     interaction is only with the krfmode=1-species.
c     Can be called for multiple flux surfaces.
c     It is not to be used in conjuction with urfmod="enabled",
c     since some storage is over-lapped.
c
c     Namelist inputs controlling this routine are prefixed by vlf. 
c  
c     subroutine vlf requires (nrf.eq.1 .and. vlfmod.eq."enabled")
c
c     With cqlpmod.ne."enabled", it provides bounce-averaged
c     QL Fokker-Planck coefficients.
c
c     990317:  It has been generalized to work with the
c              cqlpmod="enabled" option, to facilitate examination
c              of finite-length effects on RFCD efficiency.
c.................................................................
c
      include 'param.h'
      include 'comm.h'
c      common/temp_imax_old/ imax_old(lza)   ;for checking taunew.


      complex*16 cwz,cwxyp,cwxym,cei


      if (action.ne."setup") return
      if (ngen.ne.1) stop ' ngen.ne.1 in subroutine vlf'

      call vlfsetup ! determines mrfn

      ! YuP-101220: Now mrfn is known; allocate wcqlb(), cqlb(), etc.
      if (cqlpmod.eq."enabled") then
        allocate(wcqlb(iy,jx,mrfn,lz),STAT=istat) ! ls  or lz ?
        allocate(wcqlc(iy,jx,mrfn,lz),STAT=istat)
        allocate(wcqle(iy,jx,mrfn,lz),STAT=istat)
        allocate(wcqlf(iy,jx,mrfn,lz),STAT=istat)
        call bcast(wcqlb,zero,SIZE(wcqlb))
        call bcast(wcqlc,zero,SIZE(wcqlc))
        call bcast(wcqle,zero,SIZE(wcqle))
        call bcast(wcqlf,zero,SIZE(wcqlf))
      endif

      if(ASSOCIATED(cqlb)) then
        ! cqlb-cqlf are already allocated => do nothing 
      else ! Not allocated yet
        allocate(cqlb(iy,jx,lrz,mrfn),STAT=istat)
        allocate(cqlc(iy,jx,lrz,mrfn),STAT=istat)
        allocate(cqle(iy,jx,lrz,mrfn),STAT=istat)
        allocate(cqlf(iy,jx,lrz,mrfn),STAT=istat)
        cqlb_size=size(cqlb)
        call bcast(cqlb,zero,cqlb_size)
        call bcast(cqlc,zero,cqlb_size)
        call bcast(cqle,zero,cqlb_size)
        call bcast(cqlf,zero,cqlb_size)      
      endif
      


c.................................................................
c     Initialize some variables
c.................................................................

      cei=(0.,1.)

      ks=1
      bcnst=abs(bnumb(ks))*charge/(fmass(ks)*clight)
      wx=0.5

c.................................................................
c     Up-down symmetry factor symm: defined in aingeom
c**bh050820:
c**bh050820:  Trapped-particle bounce-time factor, trapfac=1., 
c**bh050820:  accounting for bounce time twice transitting bounce time.
c.................................................................
      
      trapfac=1.d0

      do i=1,iy
         alfi(i)=1.d0
      enddo
      do i=itl,itu
         alfi(i)=alfi(i)/trapfac
      enddo

c..................................................................
c     Set indicators of electrons or ions
c..................................................................

        if (kionn.eq.1) then
          signi=1.0
        else
          signi=0.0
        endif

        if (kelecg.eq.1) then
          signe=1.0
        else
          signe=0.0
        endif

c        write(*,*)'vlf:lr_',lr_
c        write(*,*)'vlf: (l,imax(l,lr_),l=1,lz)',(l,imax(l,lr_),l=1,lz)

c..................................................................
c     Main loop over harmonics, or modes
c..................................................................

      do 500  krfmode=1,mrfn   ! loop over wave modes. (ends at line~921)

c..................................................................
c     Polarizations
c..................................................................

      cwxyp=vlfeplus(krfmode)
      cwxym=vlfemin(krfmode)
c990131      cwz=sqrt(1.-cabs(cwxyp**2)-cabs(cwxym**2))
      cwz=sqrt(1.-abs(cwxyp**2)-abs(cwxym**2))


c.................................................................
c     The quasilinear diffusion coefficient is normalized by the
c     collisional diffusion  coefficient.
c     Thus the namelist variable
c     vlfdnorm=0.5*pi*bnumb(1)**2*charge**2/fmass(1)**2*abs(E**2)/
c              (vlfdnpar*v_te*D_coll).
c     Function taper below accounts for poloidal dependence of
c     the diffusion coefficient.
c.................................................................

        dene=reden(kelec,lr_)
        te=temp(kelec,lr_)*1.e3
c990131        xlog=24.-alog(sqrt(dene)/te)
        xlog=24.-log(sqrt(dene)/te)
        if (kelecg.eq.ks) then
          taucoll=3.44074e5*te**1.5/(zeff(lr_)*dene*xlog)
        elseif (kiong(1).eq.ks) then
          ti=temp(kiong(1),lr_)*1.e3
          taucoll=2.08507e7*ti**1.5
     +        /(zeff(lr_)*dene*bnumb(kiong(1))**2*xlog)
        else
          stop 'in vlf'
        endif

c..................................................................
c     Main loop over poloidal mesh
c..................................................................

c..................................................................
c     The coding from here down through the do 64-loop closely
c     follows subroutine urfpack.
c..................................................................
      !write(*,*)'vlf: starting l=1:lz at lr,nharm=', lr_,nharm(krfmode)
      do 20 l=1,lz

        if (cqlpmod.eq."enabled")then !should be one surface: indxlr_=1
           call bcast (cqlb(1,1,indxlr_,krfmode),zero,iyjx)
           call bcast (cqlc(1,1,indxlr_,krfmode),zero,iyjx)
           call bcast (cqle(1,1,indxlr_,krfmode),zero,iyjx)
           call bcast (cqlf(1,1,indxlr_,krfmode),zero,iyjx)
        endif

        tap1=taper(pol(l,lr_),vlfpol(krfmode),
     ~             vlfdpol(krfmode),vlfddpol(krfmode))
c        write(*,*) 'vlf:l,lr_,pol,vlfpol,vlfdpol,vlfddpol,taper',
c     +    l,lr_,pol(l,lr_),vlfpol(krfmode),vlfdpol(krfmode),
c     +       vlfddpol(krfmode),tap1 

        diffus=vlfdnorm(krfmode)*vth(ks,lr_)**2/taucoll
     +    *taper(pol(l,lr_),vlfpol(krfmode),
     ~           vlfdpol(krfmode),vlfddpol(krfmode))

        rr=solrz(l,lr_)
        btot=bbpsi(l,lr_)*bmidplne(lr_)
        wce=bcnst*btot

        if (vlfnpvar.eq."constant") then
          vlfnpar=vlfnp(krfmode)
          vlfdnpar=vlfdnp(krfmode)
          vlfddnpa=vlfddnp(krfmode)
          vlfnper=vlfnperp(krfmode)
        elseif (vlfnpvar.eq."1/R") then
          vlfnpar=vlfnp(krfmode)*rpcon(lr_)/rr
          vlfdnpar=vlfdnp(krfmode)*rpcon(lr_)/rr
          vlfddnpa=vlfddnp(krfmode)*rpcon(lr_)/rr
          vlfnper=vlfnperp(krfmode)*rpcon(lr_)/rr
        else
          stop 'error in vlf input'
        endif

        signn=1.
        if (vlfnpar.lt.0) signn=-1.

        omomn=1.0
        signom=1.0

c..................................................................
        if (nharm(krfmode).eq.0 .or. kionn.eq.1)  then
c..................................................................

c..................................................................
c     Lower Hybrid/Fast Wave, or non-relativistic ions
c..................................................................


c..................................................................
c     Determine the local minimum and maximum parallel speed for
c     which we assume a resonance.
c     For nharm(krfmode)=0, if vparl has opposite sign from vparu, 
c     damping will be omitted.  
c     Generally, this occurs when npar goes through
c     zero and hence the wave resonance is beyond the grid.
c     Sign reversal on vl,vu (and cosmz) is to accomodate luf.
c     Velocity integration limits are worked out using the
c     abs(wnpar).  Limits for negative wnpar can be obtained from
c     the abs(wnpar)-results by reflection about theta=pi/2.
c
c     The ion damping case is treated non-relativistically,
c     including the cyclotron damping.  The resonances are similar
c     to the electron case, in that the resonant velocity is
c     independent of perpendicular velocity.
c..................................................................

                if (kionn.eq.1) then
                  wci=bcnst*btot
                  omn=nharm(krfmode)*wci/omega(krfmode)
                  omomn=1.0-omn
c990131                  signom=sign(1.0,omomn)
                  signom=sign(one,omomn)
                endif

                vparu=clight*signom*omomn/
     1            (signn*vlfnpar-(vlfdnpar+vlfddnpa)*wx)
                vparl=clight*signom*omomn/
     1            (signn*vlfnpar+(vlfdnpar+vlfddnpa)*wx)
     
                if(vparl*vparu .gt. 0.)  then
                  vl=-abs(vparl/vnorm)
                  vu=-abs(vparu/vnorm)
                else
                  vl=-1.1
                  vu=-1.2
                endif

c..................................................................
c     Determine the lowest index in the velocity mesh that experiences
c     a resonance. sx is the normalized speed mesh, whereas x in the
c     momentum mesh.
c     jmin will be .gt.1.  If abs(vl) is off the mesh, jmin=jx+1.
c     
c     Max j is at edge of the mesh.
c..................................................................

                jmin=luf(-vl,sx,jx)
                jmax=jx
                if (jmin.gt.jx) go to 20
                if (jmin.eq.1) call urfwrong(4)
                do 25 j=1,jmin-1
                  ilim2(j)=0
                  ilim1(j)=0
 25             continue

c..................................................................
c     Determine the lowest(highest) theta index for which there is a
c     resonance at speed mesh sx(j)
c..................................................................

                ilim2(jmin)=luf(vu/sx(jmin),cosmz(1,l,lr_),iyh)
                ilim1(jmin)=luf(vl/sx(jmin),cosmz(1,l,lr_),iyh)
                call urfedge(ilim1(jmin),ilim2(jmin),vl/sx(jmin),
     +            vu/sx(jmin),l,lr_,jmin)
                ilim2(jx)=luf(vu/sx(jx),cosmz(1,l,lr_),iyh)
                ilim1(jx)=luf(vl/sx(jx),cosmz(1,l,lr_),iyh)
                call urfedge(ilim1(jx),ilim2(jx),vl/sx(jx),vu/sx(jx),
     1            l,lr_,jx)
                iupjx=ilim2(jx)
                ilwjx=ilim1(jx)
                do 30 j=jmin+1,jx-1
                  ilim2(j)=luf(vu/sx(j),cosmz(ilim2(j-1),l,lr_),iupjx-
     1              ilim2(j-1))-1+ilim2(j-1)
                  ilim1(j)=luf(vl/sx(j),cosmz(ilim1(j-1),l,lr_),ilwjx-
     1              ilim1(j-1))-1+ilim1(j-1)
                  call urfedge(ilim1(j),ilim2(j),vl/sx(j),vu/sx(j),
     1              l,lr_,j)
 30             continue

c..................................................................
              elseif (nharm(krfmode).gt.0)  then
c..................................................................

c..................................................................
c     Cyclotron damping case.  rnpar1 is the lower n_parallel and
c     rnpar2 is the upper one.  The associated resonance velocities (not
c     momentum-per-mass) are vpar1 and vpar2.  ilim1(j) and ilim2(j) are
c     the lower and upper limits of the pitch angle index.
c
c     Also, if rnpar1 and rnpar2 have opposite signs, damping is again
c       omitted.  This can give inaccuracies for very close to perpendicular
c       propagation.  It can be mitigated by choosing small 
c       vlfdnpar+vlfddnpa.
c..................................................................

                rnpar1=signn*vlfnpar-wx*(vlfdnpar+vlfddnpa)
                rnpar2=signn*vlfnpar+wx*(vlfdnpar+vlfddnpa)
                r1mn12=1.-rnpar1*rnpar1
                r1mn22=1.-rnpar2*rnpar2
                omn=nharm(krfmode)*bcnst*btot/omega(krfmode)
                omn2=omn*omn
                rad1=omn2-r1mn12
                rad2=omn2-r1mn22
c
                if (rnpar1*rnpar2.le.0.) then
                  jmin=jx+1
                  go to 20
                endif

c..................................................................
c     Ray outside of pinch point. 
c     No resonance (omega(krfmode) will be .gt. omega_ce).
c     Indicate skipping of the ray element by taking jmin=jx+1.
c..................................................................

                if (rad1.le.0.0.and.rad2.le.0.0)  then

c..................................................................
c     set flag to skip this ray element
c..................................................................

                  jmin=jx+1
                  go to 20
                endif

c..................................................................
c     The following case should not be possible:
c..................................................................

                if (rad1.gt.0.0.and.rad2.lt.0.0)  call urfwrong(6)

c..................................................................
c     Note: We let uu01 and uu02 become negative for npar**2 > 1.
c           This helps the logic go through for this case, 
c           with minimum changes from the npar**2<1 case.
c..................................................................

                if (rad2.gt.0.0.and.rad1.le.0.0)  then
                  upar01=0.0
                  uu01=0.0
                else
                  upar01=clight*omn*rnpar1/(r1mn12*vnorm)
                  uu01=clight*sqrt(rad1)/(r1mn12*vnorm)
                endif
                upar02=clight*omn*rnpar2/(r1mn22*vnorm)
                uu02=clight*sqrt(rad2)/(r1mn22*vnorm)

c..................................................................
c     upar11 and upar12 are parallel momentum-per-mass limits of the
c     resonance ellipse associated with rnpar1..., upar21 and upar22,
c     associated with rnpar2.
c..................................................................

                upar11=upar01-uu01
                upar12=upar01+uu01
                upar21=upar02-uu02
                upar22=upar02+uu02

c               Set rnpar1, rnpar2 .gt.1.0, for npar**2>1 case:
                if (rnpar1.gt.1.) upar12=1.1
                if (rnpar2.gt.1.) upar22=1.1

                jmin=luf(abs(upar21),x,jx)
                jmax=luf(upar22,x,jx)-1
c***  if(jmin.gt.jx.or.jmin.eq.jmax)  go to 20
c***  Try a less restricted condition, compatible with urfb0,...

c           write(*,*)'vlf:l,upar11,upar12,upar21,upar22,jmin,jmax',
c     +                   l,upar11,upar12,upar21,upar22,jmin,jmax
                if(jmin.gt.jx)  go to 20

                if(jmin.eq.1)  call urfwrong(4)

                do 41  j=1,jmin-1
                  ilim1(j)=0
                  ilim2(j)=0
 41             continue
                do 42  j=jmax+1,jx
                  ilim1(j)=0
                  ilim2(j)=0
 42             continue

c..................................................................
c     Treat cases of (nharm(krfmode)*omega_ce/omega(krfmode)).le.1, 
c     and .gt.1, separately
c..................................................................

                if (omn.le.1.)  then
c*bh*940227if(x(jmin).le.upar11)then
                  if(x(jmin).le.upar11 .or. rad1.le.0.0)  then
                    vres1=clight/vnorm ! c*bh*940227
                    ilim2(jmin)=1
                  else
                    vres1=clight*(1.-omn/gamma(jmin))/(rnpar1*vnorm)
                    ilim2(jmin)=luf(-vres1/sx(jmin),cosmz(1,l,lr_),iyh)
                  endif
                  vres2=clight*(1.-omn/gamma(jmin))/(rnpar2*vnorm)
                  ilim1(jmin)=luf(-vres2/sx(jmin),cosmz(1,l,lr_),iyh)
                  ii1=ilim1(jmin)
                  call urfedge(ii1,ilim2(jmin),-vres2/sx(jmin),
     +              -vres1/sx(jmin),l,lr_,jmin)
                  do 45  j=jmin+1,jmax
                    if(x(j).le.upar11)  then
                      ilim2(j)=1
                      vres2=clight*(1.-omn/gamma(j))/(rnpar2*vnorm)
                      ilim1(j)=luf(-vres2/sx(j),cosmz(1,l,lr_),iyh)
                    elseif (x(j).lt.upar12)  then
                      vres1=clight*(1.-omn/gamma(j))/(rnpar1*vnorm)
                      vres2=clight*(1.-omn/gamma(j))/(rnpar2*vnorm)
                      ilim2(j)=luf(-vres1/sx(j),cosmz(1,l,lr_),iyh)
                      ilim1(j)=luf(-vres2/sx(j),cosmz(1,l,lr_),iyh)
                    else
                      ilim2(j)=1
                      vres2=clight*(1.-omn/gamma(j))/(rnpar2*vnorm)
                      ilim1(j)=luf(-vres2/sx(j),cosmz(1,l,lr_),iyh)
                    endif
                    ii1=ilim1(j)
                    call urfedge(ii1,ilim2(j),-vres2/sx(j),-vres1/sx(j),
     1                l,lr_,j)
 45               continue


                else

c...............................................................
c     nharm(krfmode)*omega_ce/omega(krfmode) .gt.1 - case:
c...............................................................

                  uperpstr=sqrt(omn*omn-1.)*clight/vnorm
                  if(x(jmin).le.abs(upar11))  then
cbh940227
                    vres1=clight*(1.-omn/gamma(jmin))/(rnpar1*vnorm)
                    ilim1(jmin)=iy
                  else
                    vres1=clight*(1.-omn/gamma(jmin))/(rnpar1*vnorm)
                    ilim1(jmin)=luf(-vres1/sx(jmin),cosmz(1,l,lr_),iy)
                  endif
                  vres2=clight*(1.-omn/gamma(jmin))/(rnpar2*vnorm)
                  ilim2(jmin)=luf(-vres2/sx(jmin),cosmz(1,l,lr_),iy)
                  ii1=ilim1(jmin)
                  call urfedge(ii1,ilim2(jmin),-vres1/sx(jmin),
     +              -vres2/sx(jmin),l,lr_,jmin)
                  do 55  j=jmin+1,jmax
                    vres1=clight*(1.-omn/gamma(j))/(rnpar1*vnorm)
                    vres2=clight*(1.-omn/gamma(j))/(rnpar2*vnorm)
                    if(x(j).le.abs(upar11))  then
                      ilim1(j)=iy
                      ilim2(j)=luf(-vres2/sx(j),cosmz(1,l,lr_),iy)
                      ii1=ilim1(j)
                      call urfedge(ii1,ilim2(j),-vres1/sx(j),
     +                  -vres2/sx(j),l,lr_,j)
                    elseif (x(j).lt.uperpstr)  then
                      ilim2(j)=luf(-vres2/sx(j),cosmz(1,l,lr_),iy)
                      ilim1(j)=luf(-vres1/sx(j),cosmz(1,l,lr_),iy)
                      ii1=ilim1(j)
                      call urfedge(ii1,ilim2(j),-vres1/sx(j),
     +                  -vres2/sx(j),l,lr_,j)
                    elseif(x(j).lt.upar12)  then
                      ilim2(j)=luf(-vres1/sx(j),cosmz(1,l,lr_),iy)
                      ilim1(j)=luf(-vres2/sx(j),cosmz(1,l,lr_),iy)
                      ii1=ilim1(j)
                      call urfedge(ii1,ilim2(j),-vres2/sx(j),
     +                  -vres1/sx(j),l,lr_,j)
                    elseif (x(j).lt.upar22)  then
                      ilim2(j)=1
                      ilim1(j)=luf(-vres2/sx(j),cosmz(1,l,lr_),iy)
                      ii1=ilim1(j)
                      call urfedge(ii1,ilim2(j),-vres2/sx(j),
     +                  -vres1/sx(j),l,lr_,j)
                    endif


c            if (ilim2(j).gt.imax(l,lr_) .and. 
c     +          ilim2(j).lt.(iy+1-imax(l,lr_))) then
c                write(*,*)'vlf:ilim2 problem,ilim2,imax,j,l,',
c     +                          ilim2(j),imax(l,lr_),j,l
c               i=imax(l,lr_)
c                write(*,*)'x(j),upar11,upar12,upar21,upar22',
c     +                      x(j),upar11,upar12,upar21,upar22
c             write(*,*)'-vres1/sx(j),-vres2/sx(j),cosmz(i),cosmz(i+1)',
c     +        -vres1/sx(j),-vres2/sx(j),cosmz(i,l,lr_),cosmz(i+1,l,lr_)
c            endif
c            if (ilim1(j).gt.imax(l,lr_) .and. 
c     +          ilim1(j).lt.(iy+1-imax(l,lr_))) then
c                write(*,*)'vlf:ilim1 problem,ilim1,imax,j,l,',
c     +                          ilim1(j),imax(l,lr_),j,l
c                i=imax(l,lr_)
c                write(*,*)'x(j),upar11,upar12,upar21,upar22',
c     +                      x(j),upar11,upar12,upar21,upar22
c             write(*,*)'-vres1/sx(j),-vres2/sx(j),cosmz(i),cosmz(i+1)',
c     +        -vres1/sx(j),-vres2/sx(j),cosmz(i,l,lr_),cosmz(i+1,l,lr_)
c            endif


 55               continue

                endif
                 
c..................................................................
c     End of nharm(krfmode) if construction:
c..................................................................

              endif


c....................................................................
cBH010316: Fix small sliver of cosmz where vres1(j)/sx or vres2(j)/sx
cBH010316: can get beyond excluded midplane region:                    
c....................................................................

              iimax=imax(l,lr_)
              iimaxp=iy+1-imax(l,lr_)
              do j=jmin,jmax
                 ii1=ilim1(j)
                 ii2=ilim2(j)
c                    write(*,*)'fix ilim1(j),j,ilim1,iimax',
c     +                         j,ilim1(j),iimax
                 if( (ii1.le.iyh) .and. (ii1.gt.iimax)) then
c                    write(*,*)'fix ilim1(j),j,ilim1,iimax',
c     +                         j,ilim1(j),iimax
                    ilim1(j)=iimax
                 elseif( (ii1.ge.iyh+1) .and. (ii1.lt.iimaxp) ) then
c                    write(*,*)'fix ilim1(j),j,ilim1,iimaxp',
c     +                         j,ilim1(j),iimaxp
                    ilim1(j)=iimaxp
                 endif
c                     write(*,*)'fix ilim2(j),j,ilim2,iimax',
c     +                         j,ilim2(j),iimax
                 if( (ii2.le.iyh) .and. (ii2.gt.iimax) ) then
c                     write(*,*)'fix ilim2(j),j,ilim2,iimax',
c     +                         j,ilim2(j),iimax
                   ilim2(j)=iimax
                 elseif( (ii2.ge.iyh+1) .and. (ii2.lt.iimaxp) ) then
c                    write(*,*)'fix ilim2(j),j,ilim2,iimaxp',
c     +                         j,ilim2(j),iimaxp
                    ilim2(j)=iimaxp
                 endif
              enddo
c              write(*,*)'ilim1=',(ilim1(j),j=1,jx)
c              write(*,*)'ilim2=',(ilim2(j),j=1,jx)
c              write(*,*)'ifct1=',(ifct1(j),j=1,jx)
c              write(*,*)'ifct2=',(ifct2(j),j=1,jx)

c..................................................................
c     Check for errors
c..................................................................

              do 60 j=1,jx
                if(ilim1(j).gt.iy) call urfwrong(7)
                if(ilim2(j).gt.iy) call urfwrong(7)
 60           continue

c..................................................................
c     Reflect about pi/2 for negative k_par
c..................................................................

              if (signn*signom.lt.0.) then
                do 61 j=1,jx
                  iota(j)=ilim2(j)
 61             continue
                do 62 j=1,jx
                  ilim2(j)=iy+1-ilim1(j)
                  ilim1(j)=iy+1-iota(j)
 62             continue
                do 63 j=1,jx
                  iota(j)=ifct1(j)
 63             continue
                do 64 j=1,jx
                  ifct1(j)=ifct2(j)
                  ifct2(j)=iota(j)
 64             continue
              endif

c..................................................................
c     The coding from here down closely follows subroutine urfb0.
c..................................................................


c..................................................................
c     Compute three theta arrays which will be used in the evaluation
c     of B_0. These arrays are independent of speed.
c..................................................................

            do 65 i=1,iy
              cosz1(i)=cosz(i,l,lr_)*cwz
              sinz1(i)=sinz(i,l,lr_)*cwxyp
              sinz2(i)=sinz(i,l,lr_)*cwxym
 65         continue

c..................................................................
c     Determine most of the argument of the Bessel functions, also
c     some leading coefficient factors.
c..................................................................
            do 70 j=jmin,jmax
              i2=ilim2(j)
              i1=min0(iy,ilim1(j))
c              write(*,*)'do70:l,j,i2,i1',l,j,i2,i1
              argmnt(j)=x(j)*vnorm*omega(krfmode)/clight*vlfnper/wce

c.................................................................
c     Smooth out the diffusion coefficient in j to avoid
c     dividing by a near zero value in the evaluation of alfag.
c.................................................................

              rrr= (1.-nharm(krfmode)*wce/(omega(krfmode)*gamma(j)))**2
              div=1.
              if (j.gt.1) then
                div=div+1.
                rrr=rrr+
     ~           (1.-nharm(krfmode)*wce/(omega(krfmode)*gamma(j-1)))**2
              endif
              if (j.lt.jx) then
                div=div+1.
                rrr=rrr+
     ~           (1.-nharm(krfmode)*wce/(omega(krfmode)*gamma(j+1)))**2
              endif
              rrr=sqrt(rrr/div)

              alfag(j)=x(j)**2*vnorm2*diffus*(vth(ks,lr_)/clight)/rrr

              tem=nharm(krfmode)*wce/
     ~            (omega(krfmode)*gamma(j)*bbpsi(l,lr_))

c.......................................................................
c     Add in the contribution to urfb (B_0,indxlr_,krfmode) (innermost loop)
c.......................................................................

              do 71 i=i2,i1
                temc2(i)=1.
 71           continue
              temc2(i2)=dfloat(ifct2(j))/65535.d0
              temc2(i1)=dfloat(ifct1(j))/65535.d0

              do 75 i=i2,i1

c     Skip contribution if outside a test window in momentum-space:
                ull0=x(j)*coss(i,l_)
                uperp0=x(j)*sinn(i,l_)
                if (vlfparmn(krfmode).gt.ull0 .or. 
     ~              vlfparmx(krfmode).lt.ull0) then
                  go to 75
                elseif (vlfprpmn(krfmode).gt.uperp0 .or.
     +              vlfprpmx(krfmode).lt.uperp0) then
                  go to 75
                endif
                 

c     Bounce average factor:  (=1. for cqlpmod.eq."enabled"? CHECK IT)
              if (cqlpmod.ne."enabled") then
cBH091031                 if (l.eq.lz.or.lmax(i,lr_).ne.l) then
cBH091031                    ax=dtau(i,l,lr_)/tau(i,lr_)
cBH091031                 else
cBH091031                    ax=(dtau(i,l,lr_)+dtau(i,l+1,lr_))/tau(i,lr_)
cBH091031                 endif
                if (eqsym.ne."none") then       !i.e. up-down symm
                   !if not bounce interval
                   if(l.eq.lz .or. l.ne.lmax(i,lr_)) then
                      ax=dtau(i,l,lr_)/tau(i,lr_)
                   else !bounce interval: additional contribution
                      ax=(dtau(i,l,lr_)+dtau(i,l+1,lr_))/tau(i,lr_)
                   endif
                else  !eqsym="none"
                   if (l.lt.lz_bmax(lr_) .and. l.eq.lmax(i,lr_))
     +                then
                      !trapped, with tips between l and l+1 (above midplane)
                      ax=(dtau(i,l,lr_)+dtau(i,l+1,lr_))/tau(i,lr_)
                      !-YuP  Note: dtau(i,l+1,lr_)=0
                   elseif (l.gt.lz_bmax(lr_) .and. l.eq.lmax(i+iyh,lr_))
     +                then
                      !trapped, with tips between l and l-1 (below midplane)     
                      ax=(dtau(i,l,lr_)+dtau(i,l-1,lr_))/tau(i,lr_) !NB:l-1
                      !-YuP  Note: dtau(i,l-1,lr_)=0
                   else
                      !passing (i<itl), or trapped but with tips at other l;
                      !also, at l=lz_bmax, includes last trapped particle i=itl
                      !(for such particle, lmax(itl)=lz_bmax; see micxinil)
                      ax=dtau(i,l,lr_)/tau(i,lr_)
                   endif
                endif
              elseif (cqlpmod.eq."enabled") then
                 ax=1.
              endif

c     Npar spectrum factor:
              rnpar=clight*
     ~              (1.-nharm(krfmode)*wce/(omega(krfmode)*gamma(j)))
     +              /(x(j)*vnorm*cosz(i,l,lr_)/gamma(j))
c              rnfactor=
c     +        taper(rnpar,signn*vlfnpar,vlfdnpar,vlfddnpar)
              rnfactor=1.0
                indx=sinz(i,l,lr_)*argmnt(j)/bsslstp(krfmode)+1
                tem1(i)=alfag(j)*alfi(i)*abs(rnpar)*vptb(i,lr_)
     +            *ax*rnfactor*temc2(i)
c990131     1            *cabs(cosz1(i)*jb0(indx,krfmode)
     1            *abs(cosz1(i)*jb0(indx,krfmode)
     1            +sinz1(i)*(signe*jbp1(indx,krfmode)
     ~            +signi*jbm1(indx,krfmode))
     1            +sinz2(i)*(signe*jbm1(indx,krfmode)
     1            +signi*jbp1(indx,krfmode)))**2
                cqlb(i,j,indxlr_,krfmode)=cqlb(i,j,indxlr_,krfmode)
     +                                   +tem1(i)
                tem2(i)=tem1(i)/(x(j)*coss(i,l_))
     +                  *(tem-sinn(i,l_)**2)
                cqle(i,j,indxlr_,krfmode)=cqle(i,j,indxlr_,krfmode)
     +                                   +tem2(i)
   
 75           continue  !i=i2,i1 

c..................................................................
c     urfc and urff have special treatment for theta0=0,pi
c..................................................................

              i2=max0(2,i2)
              i1=min0(iy-1,i1)
              do 76  i=i2,i1

c     Skip contribution if outside a test window in momentum-space:
                ull0=x(j)*coss(i,l_)
                uperp0=x(j)*sinn(i,l_)
                if (vlfparmn(krfmode).gt.ull0 .or. 
     ~              vlfparmx(krfmode).lt.ull0) then
                  go to 76
                elseif (vlfprpmn(krfmode).gt.uperp0 .or.
     +              vlfprpmx(krfmode).lt.uperp0) then
                  go to 76
                endif

                cqlc(i,j,indxlr_,krfmode)=cqlc(i,j,indxlr_,krfmode)
     +                                   +tem2(i)/sinn(i,l_)
                cqlf(i,j,indxlr_,krfmode)=cqlf(i,j,indxlr_,krfmode)
     +         +tem2(i)/(x(j)*sinn(i,l_)*coss(i,l_))*(tem-sinn(i,l_)**2)
 76           continue ! i=i2,i1
              cqlc(1,j,indxlr_,krfmode)=0.0
              cqlc(iy,j,indxlr_,krfmode)=0.0
              cqlf(1,j,indxlr_,krfmode)=0.0
              cqlf(iy,j,indxlr_,krfmode)=0.0

c..................................................................
c     Make certain indx is within table bounds
c..................................................................

              if (indx.gt.nbssltbl) call urfwrong(1)
 70         continue ! j loop

c..................................................................
c     Transfer coeffs to wcqlb,..., if cqlpmod.eq."enabled"
c..................................................................

           if (cqlpmod.eq."enabled") then ! should be one surface
           call dcopy(iyjx,cqlb(1,1,1,krfmode),1,wcqlb(1,1,krfmode,l),1)
           call dcopy(iyjx,cqlc(1,1,1,krfmode),1,wcqlc(1,1,krfmode,l),1)
           call dcopy(iyjx,cqle(1,1,1,krfmode),1,wcqle(1,1,krfmode,l),1)
           call dcopy(iyjx,cqlf(1,1,1,krfmode),1,wcqlf(1,1,krfmode,l),1)
           endif

c..................................................................
c     End of main loop over poloidal mesh (l)
c..................................................................
 20   continue

c..................................................................
c     Store coeffs in cqlb... or wcqlb..., depending on cqlpmod.
c..................................................................

      if (cqlpmod.ne."enabled") then


c..................................................................
c     Renormalize for code.
c..................................................................

        do 84  j=1,jx
          do 85  i=1,iy
            cqlb(i,j,indxlr_,krfmode)=cqlb(i,j,indxlr_,krfmode)/vnorm4
            cqlc(i,j,indxlr_,krfmode)=cqlc(i,j,indxlr_,krfmode)/vnorm4
            cqle(i,j,indxlr_,krfmode)=cqle(i,j,indxlr_,krfmode)/vnorm4
            cqlf(i,j,indxlr_,krfmode)=cqlf(i,j,indxlr_,krfmode)/vnorm4
 85       continue
 84     continue
c        write(*,*) 'renormed cqlb:',
c     +             ((i,j,cqlb(i,j,indxlr_,krfmode),i=1,iy),j=21,22)

c.......................................................................
c     Symmetrize about pi/2 in trapped region.
c.......................................................................

        do 86 j=1,jx
          do 87 i=itl,itu
            temp1(i,j)=cqlb(i,j,indxlr_,krfmode)
            temp2(i,j)=cqlc(i,j,indxlr_,krfmode)
            temp3(i,j)=cqle(i,j,indxlr_,krfmode)
            temp4(i,j)=cqlf(i,j,indxlr_,krfmode)
 87       continue
 86     continue
        do 88 j=1,jx
          do 89 i=itl,itu
            cqlb(i,j,indxlr_,krfmode)=(temp1(i,j)+temp1(iy+1-i,j))*.5
            cqlc(i,j,indxlr_,krfmode)=(temp2(i,j)-temp2(iy+1-i,j))*.5
            cqle(i,j,indxlr_,krfmode)=(temp3(i,j)-temp3(iy+1-i,j))*.5
            cqlf(i,j,indxlr_,krfmode)=(temp4(i,j)+temp4(iy+1-i,j))*.5
 89       continue
 88     continue
 
        ! CHECK:  cqlb*cqlf-cqlc*cqle =0 or not?
        BF_CE_0=0.d0
        do j=1,jx
        do i=1,iy
           B0F0= cqlb(i,j,indxlr_,krfmode)*cqlf(i,j,indxlr_,krfmode)
           C0E0= cqlc(i,j,indxlr_,krfmode)*cqle(i,j,indxlr_,krfmode)
           BF_CE= B0F0-C0E0
           if( abs(BF_CE) .gt. 1.e-8*(abs(B0F0)+abs(C0E0)) ) then
              write(*,*)'vlf: abs(BF-CE) > 1e-8*abs(BF) ',j,i,B0F0,C0E0
              BF_CE_0=BF_CE
           endif
        enddo
        enddo
        if(BF_CE_0.ne.0.d0) pause ! Never happens, 
        !so BF-CE=0 with machine accuracy

c.......................................................................
c     Taper diffusion over last 10 point of velocity, 
c     if ineg="trunc_d"
c.......................................................................

      if (ineg.eq."trunc_d") then
        if (jx.le.11) stop 'vlf: Need to set jx>11'
        do 90 j=jx-11,jx
          do 91 i=1,iy
           cqlb(i,j,indxlr_,krfmode)=truncd(j)*cqlb(i,j,indxlr_,krfmode)
           cqlc(i,j,indxlr_,krfmode)=truncd(j)*cqlc(i,j,indxlr_,krfmode)
           cqle(i,j,indxlr_,krfmode)=truncd(j)*cqle(i,j,indxlr_,krfmode)
           cqlf(i,j,indxlr_,krfmode)=truncd(j)*cqlf(i,j,indxlr_,krfmode)
 91       continue
 90     continue
      endif

c..................................................................
c     Store coeffs in cqlb... or wcqlb..., depending on cqlpmod.
c..................................................................

      elseif (cqlpmod.eq."enabled") then


c..................................................................
c     Renormalize for code.
c..................................................................
      do 200 l=1,lz

        do 184  j=1,jx
          do 185  i=1,iy
            wcqlb(i,j,krfmode,l)=wcqlb(i,j,krfmode,l)/vnorm4
            wcqlc(i,j,krfmode,l)=wcqlc(i,j,krfmode,l)/vnorm4
            wcqle(i,j,krfmode,l)=wcqle(i,j,krfmode,l)/vnorm4
            wcqlf(i,j,krfmode,l)=wcqlf(i,j,krfmode,l)/vnorm4
 185       continue
 184     continue

c..................................................................
c     NO Symmetrization about the pass/trapped boundary.
c..................................................................

c..................................................................
c     Taper diffusion over last 10 point of velocity, 
c     if ineg="trunc_d"
c..................................................................

      if (ineg.eq."trunc_d") then
        if (jx.le.11) stop 'vlf: Need to set jx>11'
        do 190 j=jx-11,jx
          do 191 i=1,iy
            wcqlb(i,j,krfmode,l)=truncd(j)*wcqlb(i,j,krfmode,l)
            wcqlc(i,j,krfmode,l)=truncd(j)*wcqlc(i,j,krfmode,l)
            wcqle(i,j,krfmode,l)=truncd(j)*wcqle(i,j,krfmode,l)
            wcqlf(i,j,krfmode,l)=truncd(j)*wcqlf(i,j,krfmode,l)
 191       continue
 190     continue
      endif

 200  continue

c..................................................................
c     End of if on cqlpmod
c..................................................................

      endif
      
c..................................................................
c     End of main loop over modes and harmonics
c..................................................................

 500  continue !  krfmode=1,mrfn ! loop over wave modes. 


c..................................................................
c     Plotting of cqlb diffusion coefficient
c..................................................................

      call vlfbplt

      return
      end


