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
      if (machine.eq."toroidal".or.machine.eq."mirror") then

cFollowing makes more sense[BH:990816]    if (l_.eq.lrzmax) call eqalloc
        if (lr_.eq.lrzmax) call eqalloc

c..................................................................
c     If the "eq"uilibrium module is utilized...
c..................................................................

        if (eqmod.eq."disabled") go to 1010

c..................................................................
c     Determine orbit information from equilibrium psi
c     eqsource="ellipse" provides a primitive method for forcing
c     elliptical cross sections...it is a debugging tool and not
c     a standard running mode.
c..................................................................

        if (eqsource.ne."ellipse") then
          dum=0.d0
          call equilib(dum,dum,0,dum,dum,dum,dum) ! setup
          ! Call with index=0 for setup.
        else
          call eqelpse
        endif
	
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
        
        call eqcoord !-> eqfndpsi -> eqorbit -> trace flux surface

c..................................................................
c     btor is the magnetic field at the nominal magnetic axis - radmaj
c     bpsi is mod B(z) / mod B (midplane)
c.....................................................................

        if (eqsym.ne."none") then
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

        etll0(lr_)=1./sqrt(1.+(bthr(lr_)/btoru(lr_))**2)

c..................................................................
c     btor0(lr_),bthr0(lr_), and bmod0(lr_) are the toroidal, poloidal and
c     total magnetic fields at the outer midplane of the flux
c     surface. They are defined in module "eq".
c..................................................................

        bmod0(lr_)=sqrt(btor0(lr_)**2+bthr0(lr_)**2)

c..................................................................
c     qsafety(lr_) is the safety factor.
c..................................................................

C%OS  qsafety(lr_)=rovera(lr_)*(radmin*btoru(lr_))/(radmaj*bthr(lr_))
CBH020914   qsafety(lr_)=rgeom(lr_)*btoru(lr_) / (r0geom(lr_)*bthr(lr_))
CBH020914   Still very inaccurate.  Interpolate from the eqdsk data.
        itab(1)=1
        itab(2)=0
        itab(3)=0
        call terp1(nfp,psiar,qar,d2qar,epsicon(lr_),1,tab,itab)
        qsafety(lr_)=tab(1)

c....................................................................
c     eps0 and eps(lr_) are the inverse aspect ratios at the plasma
c     edge and for the relevant flux surface respectively.
c     old bug: eps(lr_)=rovera(lr_)*radmin/radmaj (radmin=rhomax)
c     new def: eps(lr_)=(Rmax-Rmin)/(Rmax+Rmin) of each flux surface
c.....................................................................

        eps0=rgeomp / r0geomp
        eps(lr_)=rgeom(lr_) / r0geom(lr_)
c
        go to 1011

c..................................................................
c     Below:
c     Standard toroidal circular geometry - does not utilize "eq" module
c..................................................................

 1010   continue

c....................................................................
c     eps0 and eps(lr_) are the inverse aspect ratios at the plasma
c     edge and for the relevant flux surface respectively.
c.....................................................................

        rmag=radmaj
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
        psir=flxfn(rovera(lr_))
        etll0(lr_)=1./sqrt(1.+(bthr(lr_)/btoru(lr_))**2)

c.....................................................................
c     btor0(lr_), bthr0(lr_) and bmod0(lr_) are the toroidal, poloidal and
c     magnetic fields at the outer midplane of the flux surface (R+r/a
c.....................................................................

        btor0(lr_)=btor/(1.+eps(lr_))
        bthr0(lr_)=bthr(lr_)/(1.+eps(lr_))
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
        endif
      else
        call diagwrng(3)
      endif
      return
      end
