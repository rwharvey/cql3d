c
c
      subroutine cfpgamma
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c................................................................
c     calculate the Coulomb logarithm 
c     and an energy dependent factor
c................................................................

      include 'param.h'
      include 'comm.h'


c..................................................................
c     If log is set by fiat,  set it and move to gamefac
c..................................................................


      if (gamaset.gt.zero) then
         do 12 i=1,ntotal
            do 13 k=1,ntotal
               gama(i,k)=gamaset
 13         continue
 12      continue
         goto 99
      endif


      energ=ergtkev*energy(kelec,lr_)
      dense=reden(kelec,lr_)
      if (cqlpmod .eq. "enabled") then
        energ=ergtkev*enrgypa(kelec,ls_)
        dense=denpar(kelec,ls_)
      endif

c..................................................................
c     If two electron species exist - do an average for deby calc.
c..................................................................

      if (colmodl.eq.1 .or. (colmodl.eq.3 .and. kelecm.ne.0)) then
        dense=reden(kelecm,lr_)
        energ=reden(kelecm,lr_)*energy(kelecm,lr_)*ergtkev/dense
        if (cqlpmod .eq. "enabled") then
          dense=denpar(kelecm,ls_)
          energ=enrgypa(kelecm,ls_)*ergtkev
        endif

      elseif (kelecg.gt.1 .and. kelecm.gt. 1) then
        dense=reden(kelec,lr_)+reden(kelecm,lr_)
        km=kelecm
        energ=(reden(kelecg,lr_)*energy(kelecg,lr_)+reden(km,lr_)*
     1    energy(km,lr_))*ergtkev/dense
        if (cqlpmod .eq. "enabled") then
          dense=denpar(kelec,ls_)+denpar(kelecm,ls_)
          energ=(denpar(kelecg,ls_)*enrgypa(kelecg,ls_)+denpar(km,ls_)*
     1      enrgypa(km,ls_))*ergtkev/dense
        endif
      endif
      deby=sqrt(energ/(dense*6.*pi*charge**2))
      do 10 i=1,ntotal
        si=energy(i,lr_)/fmass(i)
        if (cqlpmod .eq. "enabled") si=enrgypa(i,ls_)/fmass(i)
        do 11 k=1,ntotal
          sk=energy(k,lr_)/fmass(k)
          if (cqlpmod .eq. "enabled") sk=enrgypa(k,ls_)/fmass(k)
c990131          sf=amax1(si,sk)
          sf=max(si,sk)
          vikf=sqrt(sf*ergtkev*2.)
          gam3=gamt(i,k)*deby*vikf
c990131          gama(i,k)=(alog(gam3)-.5)
          gama(i,k)=(log(gam3)-.5)
 11     continue
 10   continue


 99   continue

c..................................................................
c     If kelecg.eq.1 and gamafac.eq."enabled", set energy
c       dependant factor gamefac(j) for Coulomb logrithm. 
c     The energy dependent Coulomb log will be
c       alog(flamcql + min(sqrt(gamma-1),1)*(flamrp-flamcql)),
c       where flamcql is the argument of cql e-e Coulomb log,
c       and flamrp is the argument of the Rosenbluth-Putvinski
c       (Nucl. Fus. 1997) high energy electron Coulomb log.
c       The sqrt(gamma-1)-factor is chosen in accord with
c       the energy factor in the the CQL Coulomb log
c       (cf., Killeen et al. book).  The CQL Coulomb log
c       reduces to the NRL, Te.gt.10eV expression.  
c     Some further incites into this factor can be obtained
c       from Physics Vade Mecum, H.L. Anderson, Ed.
c       2nd edition, AIP, NY (1989); Sect. 16.07.F, "energy
c       loss due to scattering from atomic electrons ... given
c       by Moller scattering, with I=(electronic charge)**2/(
c       Debye length) [Rosenbluth, personal comm. 1996].
c     (bh, 980505)  
c..................................................................

      do 14 j=1,jx
         gamefac(j)=1.0
 14   continue

      if (kelecg.eq.1 .and. ngen.eq.1 .and. gamafac.eq."enabled") then

         flamcql=exp(gama(1,1))
         omegape=5.64e4*sqrt(reden(1,lr_))
         flamrp1=(2**.25)*fmass(1)*clite2/(1.0546e-27*omegape)

         do  j=1,jx
c990131            flam=flamcql+min(sqrt(gamma(j)-1.),1.)*
            flam=flamcql+min(sqrt(gamma(j)-1.),one)*
     +           (flamrp1*gamma(j)**1.5 - flamcql)
c990131            gamefac(j)=alog(flam)/gama(1,1)
            gamefac(j)=log(flam)/gama(1,1)
         enddo
      endif            


      return
      end
