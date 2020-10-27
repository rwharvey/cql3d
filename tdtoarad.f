c
c
      subroutine tdtoarad
      implicit integer (i-n), real*8 (a-h,o-z)

c     Storing values on whole mesh lrzmax and equilibrium quantities
c     recall:
c     lr_     : radial index
c     l_      : spatial variable index
c     lmdpln_ : midplane at l_=lmdpln_
c
      include 'param.h'
      include 'comm.h'


c..................................................................
c     BH091022:  Looks like a somewhat orphaned idea.  The "-z"
c                variable here a not used much elsewhere, and
c                original (non-z) variables could be used.
c                Seems to create confusion of variable names.
c                Clean it up sometime.....
c                No references to l_ or lmdpln_?
c..................................................................

c..................................................................
c     Copies flux surface dependent quantities into radial arrays
c..................................................................

c..................................................................
c     rfpwrz(k,lr_) is the local rf power deposited in species k
c     at flux surface lr_ in Watts/cc (averaged over the orbit)
c     currz(k,lr_) is the current of species k at lr_ in Amps (averaged
c     currtz(lr_) is the sum of currz(k,lr_) over ionic general species.
c     currmtz(lmdpln_) is the sum of the midplane currents over ionic species
c     curr(m)tpz(lr_) is curr(m)tz plus the electron contribution to
c     the current.
c..................................................................


      do 21 k=1,ngen
        if (nso.gt.0) then
          numsrce=nso
        else
          numsrce=1
        endif
        do 23 m=1,numsrce
          asorz(k,m,lr_)=asor(k,m,lr_)
 23     continue
 21   continue

      if (eqmod.eq."enabled") then
        area(lr_)=areacon(lr_)
        vol(lr_)=volcon(lr_)
        equilpsi(lr_)=epsicon(lr_)
        bmdplne(lr_)=bmidplne(lr_)
        onovrpz(lr_,1)=onovrp(1,lr_)
        onovrpz(lr_,2)=onovrp(2,lr_)
        aspin(lr_)=(rpcon(lr_)-rmcon(lr_))/(rpcon(lr_)+rmcon(lr_))
        rmconz(lr_)=rmcon(lr_)
        rpconz(lr_)=rpcon(lr_)
        fpsiz(lr_)=fpsi(lr_)
        bpolsqaz(lr_)=bpolsqa(lr_)

c..................................................................
c     Compute the plasma pressure (cgs)
c     if (kpress(k).ne. "enabled") this species is not used
c     in calculation of total pressure.
c..................................................................

        phot=0.
        do 10 k=1,ngen
          if (kpress(k).eq."disabled") go to 10
          phot=phot+reden(k,lr_)*energy(k,lr_)*2./3.*ergtkev
 10     continue
        prest(lr_)=phot
        do 20 k=ngen+1,ntotal
          if (kpress(k).eq."disabled") go to 20
          prest(lr_)=prest(lr_)+reden(k,lr_)*energy(k,lr_)*2./3.*
     *      ergtkev
 20     continue
 30     continue
      endif ! (eqmod.eq."enabled")
      
      !YuP[2020-01-29] Make these arrays available in general case
      ! of eqmod :
      rmconz(lr_)=rmcon(lr_)
      rpconz(lr_)=rpcon(lr_)
      
      return
      end
