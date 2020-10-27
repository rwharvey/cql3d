c
c
      subroutine aingeom
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     This routine is called repeatedly from tdinitl, for the
c     range of values of lr_=lrzmax,1,-1.
c     Alternatively, for lrzmax=1, it is called from achief1.
c
c     aingeom controls the flux surface geometry and sets up magnetic
c     fields, bpsi=(B(z))/B(0)); it also calls the "eq"uilibrium module
c     which utilizes an MHD equilibrium code "eqdsk" file, if required.
c     A heuristic mirror scenario (supplied by Gary Smith, LLNL) is also
c     available.
c..................................................................

      include 'param.h'
      include 'comm.h'

            
c..................................................................
c     Toroidal scenario first.
c..................................................................

cBH151202 Attempting to get code working with machine.eq.'mirror'
cBH151202      if (machine.eq."toroidal") then
      if ((machine.eq."toroidal").or.(machine.eq."mirror")) then

cFollowing makes more sense[BH:990816]    if (l_.eq.lrzmax) call eqalloc
        if (lr_.eq.lrzmax) call eqalloc

c..................................................................
c     If the "eq"uilibrium module is utilized...
c..................................................................

        !YuP[2020-01-29] Moved definition of symm here,
        ! before "if (eqmod.eq."disabled") go to 1010"
        ! Need to check how "symm" is used in case of eqmod="disabled".
        !YuP[03/26/2015] Define symm, to be used throughout the code
        if (eqsym.ne."none") then 
           ! up-dn symmetry: everything is evaluated over 1/2 surface
           symm=2.d0
           ! Factor symm=2 is because solrz(l,ll),solzz(l,ll) is
           ! over half-surface, and so is ddarea,ddvol,
           ! and also vtau_fow() or tau() are evaluated over half-orbits.
        else
           symm=1.d0
        endif

        if (eqmod.eq."disabled") go to 1010

c..................................................................
c     Determine orbit information from equilibrium psi
c     eqsource="ellipse" provides a primitive method for forcing
c     elliptical cross sections...it is a debugging tool and not
c     a standard running mode.
c..................................................................

        if (eqsource.eq."ellipse") then
          call eqelpse ! Set er(ir),ez(iz) grids and epsi(ir,iz)
        elseif (eqsource.eq."miller") then
          call eq_miller_set ! Set er(ir),ez(iz) grids and epsi(ir,iz)
        elseif (eqsource.eq."mirror1") then
          STOP 'eqsource=mirror1 is not available in this CQL3D version'
          !call eq_mirror1_set ! Set er(ir),ez(iz) grids and epsi(ir,iz)
        elseif (eqsource.eq."eqdsk"   .or.  !Includes machine="mirror"
     +          eqsource.eq."topeol"  .or.
     +          eqsource.eq."tsc"   ) then
          ! Set Equilibrium by reading corresponding data file:
          dum=0.d0
          call equilib(dum,dum,0,dum,dum,dum,dum) ! setup
          ! Called with index=0 for setup.
        else
          WRITE(*,*)'aingeom: wrong eqsource? eqsource=',eqsource
          stop
        endif

        call eqcoord !-> eqfndpsi -> eqorbit -> trace flux surface

c..................................................................
c     btor is the magnetic field at the nominal magnetic axis - radmaj
c     bpsi is mod B(z) / mod B (midplane)
c.....................................................................

        if (eqsym.ne."none") then
           ! up-dn symmetry: everything is evaluated over 1/2 surface
           psimx(lr_)=bpsi(lorbit(lr_),lr_)
        else
           psimx(lr_)=bpsi_max(lr_)  !=bpsi(lbpsi_max(lr_))
        endif
	
        if (trapmod.eq."enabled") then
c          Reduce B/B_min by the amount trapredc*(B/B_min-1.)
c          This is for purposes of heuristic check of trapping effects.
           psimx(lr_)=psimx(lr_)-trapredc*(psimx(lr_)-1.)
        endif

c..................................................................
c     Define: On the flux surface of interest we have
c     bthr(lr_) - the poloidal magnetic field at theta-poloidal = pi/2.
c     btoru(lr_) - the toroidal magnetic field at the same position. The
c     radial coordinate erhocon(lr_)=rovera(lr_)*rhomax.
c.....................................................................

        if(machine.eq."mirror")then
          etll0(lr_)=1.d0 ! any number; not used.
        else
          etll0(lr_)=1./sqrt(1.+(bthr(lr_)/btoru(lr_))**2)
        endif

c..................................................................
c     btor0(lr_),bthr0(lr_), and bmod0(lr_) are the toroidal, poloidal and
c     total magnetic fields at the outer midplane of the flux
c     surface. They are defined in module "eq".
c..................................................................

        bmod0(lr_)=sqrt(btor0(lr_)**2+bthr0(lr_)**2)

        if(machine.eq."toroidal")then
           ! Define qsafety(lr_) == the safety factor.
C%OS  qsafety(lr_)=rovera(lr_)*(radmin*btoru(lr_))/(radmaj*bthr(lr_))
CBH020914   qsafety(lr_)=rgeom(lr_)*btoru(lr_) / (r0geom(lr_)*bthr(lr_))
CBH020914   Still very inaccurate.  Interpolate from the eqdsk data.
          itab(1)=1
          itab(2)=0
          itab(3)=0
          call terp1(nfp,psiar,qar,d2qar,epsicon(lr_),1,tab,itab)
          !sub terp1(n, x,   f,       w,        y,       int,tab,itab)
          !             x(n),f(n*int),w(n*int), y_scalar,1,  tab(3),itab(3)

          qsafety(lr_)=tab(1)
        else
          qsafety(lr_)=1.d0 ! just any value for a mirror machine
        endif

c....................................................................
c     eps0 and eps(lr_) are the inverse aspect ratios at the plasma
c     edge and for the relevant flux surface respectively.
c     old bug: eps(lr_)=rovera(lr_)*radmin/radmaj (radmin=rhomax)
c     new def: eps(lr_)=(Rmax-Rmin)/(Rmax+Rmin) of each flux surface
c.....................................................................
        if(machine.eq."mirror")then
          eps0= 1.d0      ! Not relevant: any number
          eps(lr_)= 1.d0  ! Not relevant: any number
        else
          eps0=rgeomp / r0geomp
          eps(lr_)=rgeom(lr_) / r0geom(lr_)
        endif
c
        go to 1011

c..................................................................
c     Below:
c     Standard toroidal circular geometry - does not utilize "eq" module
c..................................................................

 1010   continue ! handle to skip lines 51-133, if(eqmod.eq."disabled")
        ! Now - the case of eqmod='disabled', up to 1011_continue
c....................................................................
c     eps0 and eps(lr_) are the inverse aspect ratios at the plasma
c     edge and for the relevant flux surface respectively.
c.....................................................................

        if(lr_.eq.lrzmax)then  !YuP[2020-01-29] Added, for eqmod='disabled'
          !Note: aingeom is called for each lr_, in descending order
          rmag=radmaj
          zmag=0.d0  !YuP[2020-01-29] Added, for eqmod='disabled'
          rmincon=radmaj-0.5*rbox !YuP[2020-01-29] Added, for eqmod='disabled'
          rmaxcon=radmaj+0.5*rbox !YuP[2020-01-29] Added, for eqmod='disabled'
          ! Set er() and ez() grids, for eqmod='disabled' 
          zst=zbox/(nnz-1)
          rst=rbox/(nnr-1)
          er(1)=rmincon
          ez(1)=-0.5*zbox
          do nn=2,nnr
            er(nn)=er(nn-1)+rst
          enddo
          do nn=2,nnz
            ez(nn)=ez(nn-1)+zst
          enddo
          psimag=flxfn(0.d0) !YuP[2020-01-31] Added, for eqmod=disabled
          psilim=flxfn(1.d0) !YuP[2020-01-31] Added, for eqmod=disabled
        endif  !YuP[2020-01-29] Added, for eqmod='disabled'
        
        rhomax=radmin
        rgeomp=radmin
        rgeom(lr_)=rovera(lr_) * radmin
        zgeomp=radmin
        zgeom(lr_)=rovera(lr_) * radmin
        r0geomp=radmaj
        r0geom(lr_)=radmaj
        rpcon(lr_)=r0geom(lr_) + rgeom(lr_)
        rmcon(lr_)=r0geom(lr_) - rgeom(lr_)
        
        zpcon(lr_)=0.
        zmcon(lr_)=0.
        eps0=radmin/radmaj
        eps(lr_)=rovera(lr_)*radmin/radmaj

        erhocon(lr_)=rovera(lr_)*radmin
        areacon(lr_)=pi*erhocon(lr_)**2
        volcon(lr_)=2.*pi*rmag*areacon(lr_)
        
c.....................................................................
c     psimx(lr_) is the maximum value of bpsi (R+a/R-a)
c.....................................................................

        psimx(lr_)=(1.+eps(lr_)) / (1.-eps(lr_))
        if (trapmod.eq."enabled") then
c          Reduce B/B_min by the amount trapredc*(B/B_min-1.)
c          This is for purposes of heuristic check of trapping effects.
           psimx(lr_)=psimx(lr_)-trapredc*(psimx(lr_)-1.)
        endif

c.....................................................................
c     Define - On the flux surface of interest we have:
c     bthr(lr_) is the poloidal magnetic field at theta-poloidal = pi/2.
c     btoru(lr_) is the toroidal magnetic field at the same position. The
c     radial coordinate is rovera(lr_)*radmin.
c     call routine to compute bthr(lr_) btor..
c.....................................................................

        btoru(lr_)=btor
        psir=flxfn(rovera(lr_)) !-> also, get bthr(lr_) [poloidal field]
        etll0(lr_)=1./sqrt(1.+(bthr(lr_)/btoru(lr_))**2)

c.....................................................................
c     btor0(lr_), bthr0(lr_) and bmod0(lr_) are the toroidal, poloidal and
c     magnetic fields at the outer midplane of the flux surface (R+r/a
c.....................................................................

        btor0(lr_)=btor/(1.+eps(lr_))
        bthr0(lr_)=bthr(lr_)/(1.+eps(lr_)) ! bthr is from func.flxfn (in comm.h)
        bmod0(lr_)=sqrt(btor0(lr_)**2+bthr0(lr_)**2)
        bmidplne(lr_)=bmod0(lr_)

c.....................................................................
c     qsafety(lr_) is the safety factor..
c.....................................................................

C%OS  qsafety(lr_)=rovera(lr_)*(radmin*btoru(lr_))/(radmaj*bthr(lr_))
        qsafety(lr_)=rgeom(lr_)*btoru(lr_) / (r0geom(lr_)*bthr(lr_))

c.....................................................................
c     zmax(lr_) is the maximum gyro-orbit field length
c.....................................................................

        zmax(lr_)=qsafety(lr_)*pi*radmaj*bmod0(lr_)/btor0(lr_)

C-----------------------------------------------------------------------
c     psiavg(i,lr_)=<bpsi**i >, onovrp(i,lr_)=<1/R**i>, fpsi(lr_)=R*B_tor
c     psiovr(lr_)=<bpsi/R> (=psifct*onovrp(1,lr_))
c     onovpsir3(lr_)=<1./(bpsi*r**3)>
C-----------------------------------------------------------------------

        psiovr(lr_)=(1.+eps(lr_)) / sqrt(1.-eps(lr_)**2)/r0geom(lr_)
        psiavg(1,lr_)=1. + eps(lr_)
        psiavg(2,lr_)=(1.+eps(lr_))**2 / sqrt(1.-eps(lr_)**2)
        onovrp(1,lr_)=1. / r0geom(lr_)
        onovrp(2,lr_)=1. / r0geom(lr_)**2 / sqrt(1.-eps(lr_)**2)
        fpsi(lr_)=r0geom(lr_) * btoru(lr_)
        flxavgd(lr_)=zmax(lr_) / bmod0(lr_) / psiavg(1,lr_)
        onovpsir3(lr_)=onovrp(2,lr_)/rpcon(lr_)
        epsicon(lr_)=psir  !YuP[2020-01-31] Added, for eqmod=disabled
        equilpsi(lr_)=psir !YuP[2020-01-31] Added, for eqmod=disabled

 1011   continue

        zmaxi(lr_)=1.0/ zmax(lr_)
        pibzmax(lr_)=  pi*zmaxi(lr_)


c.....................................................................
c     If this is a mirror machine use Gary Smith's model for the
c     bpsi (magnetic field) dependence..
c.....................................................................

      else if (machine.eq."mirror") then
      
        eps(lr_)=(rmirror-1.)/(rmirror+1.)
        psimx(lr_)=(1.+eps(lr_))/(1.-eps(lr_))
        if (trapmod.eq."enabled") then
c          Reduce B/B_min by the amount trapredc*(B/B_min-1.)
c          This is for purposes of heuristic check of trapping effects.
           psimx(lr_)=psimx(lr_)-trapredc*(psimx(lr_)-1.)
        endif
        if (psimodel.eq."smith") then
          toll=1.e-12
          niter=0
          gsla2=gsla*gsla
          gslb2=gslb*gslb
          gszb=zmax(lr_)*(1.+gslb2/(zmax(lr_)**2-(rmirror-1.)*gsla2))
          gsa1=zmax(lr_)/gsla2
          gsa2=(rmirror-1.-(zmax(lr_)/gsla)**2)/gslb2
 20       continue
          gszb2=gszb*gszb
          gszm=zmax(lr_)-gszb
          gszm2=gszm*gszm
          gszp=zmax(lr_)+gszb
          gszp2=gszp*gszp
          gsem=exp(-gszm2/gslb2)
          gsep=exp(-gszp2/gslb2)
          gseb=exp(-gszb2/gslb2)
          gsn=gszm*gsem+gszp*gsep
          gsdnzb=(2.*gszm2/gslb2-1.)*gsem-(2.*gszp2/gslb2-1.)*gsep
          gsd=2.*rmirror*gseb-gsem-gsep
          gsddzb=2.*(-rmirror*gszb*gseb-gszm*gsem+gszp*gsep)/gslb2
          gszbcorr=((gsa1*gsd)/(gsa2*gsn)+1.)/(gsdnzb/gsn-gsddzb/gsd)
          gszbn=gszb-gszbcorr
          relerr=abs((gszbn-gszb)/gszb)
          gszb=gszbn
          if (relerr.lt.toll) then
            gszb2=gszb*gszb
            gszm=zmax(lr_)-gszb
            gszm2=gszm*gszm
            gszp=zmax(lr_)+gszb
            gszp2=gszp*gszp
            gsem=exp(-gszm2/gslb2)
            gsep=exp(-gszp2/gslb2)
            gseb=exp(-gszb2/gslb2)
            gsb=(1.+(zmax(lr_)/gsla)**2-rmirror)/(2.*rmirror*gseb
     1        -gsem-gsep)
            gsnm=1.+2.*gsb*gseb
          else
            niter=niter+1
            if (niter.gt.100) call diagwrng(2)
            goto 20
          endif
        endif ! (psimodel.eq."smith")
        
        if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
        write(*,*)'aingeom/psimodel=smith : gsb,gszb,gsnm=',
     1      gsb,gszb,gsnm
        write(*,*)'aingeom/psimodel=smith : ez(1),ez(nnz)=',
     1      ez(1),ez(nnz)
        endif
        
      else  ! machine=...
      
        call diagwrng(3)
        
      endif ! machine=...
      
      return
      end
