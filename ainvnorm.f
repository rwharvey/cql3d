c
c
      subroutine ainvnorm
      implicit integer (i-n), real*8 (a-h,o-z)
c.....................................................................
c     Determine mesh normalization constant vnorm.
c     vnorm is the maximum velocity for a non-
c     relativistic mesh. enorm is the maximum energy.
c     For a relativistic mesh vnorm goes over to maximum
c     momentum/unit rest mass.
c     Running electrons and ions in tandem as general species
c     requires special treatment.
c.....................................................................

      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE


      if (tandem.eq."enabled") then
        xfac=-1.
        icase=1
        enorm=enorme
        if (relativ.eq."disabled") then
          vnorm=sqrt(enorm*ergtkev*2./fmass(kelecg))
        else
          vnorm=sqrt((enorm*ergtkev/(fmass(kelecg)*clite2)+1.)**2-1.)*
     &      clight
        endif
        xlwr=sqrt(enormi*fmass(kelecg)/(enorme*fmass(kionn)))
        xpctlwr=.65
        xmdl=1.
        xpctmdl=.35
CMPIINSERT_IF_RANK_EQ_0
        WRITE(*,*)' WARNING/ainvnorm: For tandem=enabled, 
     +   xfac,xlwr,xpctlwr,xmdl,xpctmdl are reset'
        WRITE(*,*)' ainvnorm: xfac,xlwr,xpctlwr,xmdl,xpctmdl=',
     +                        xfac,xlwr,xpctlwr,xmdl,xpctmdl   
CMPIINSERT_ENDIF_RANK
        
      elseif((kenorm.ge.1).and.(kenorm.le.ngen).and.(enorm.gt.0.)) then
        icase=2
        if (relativ.eq."disabled") then
          vnorm=sqrt(enorm*ergtkev*2./fmass(kenorm))
        else
          vnorm=sqrt((enorm*ergtkev/(fmass(kenorm)*clite2)+1.)**2-1.)*
     &      clight
        endif
      else if (vnorm.gt.0.) then
        icase=3
        if (relativ.eq."disabled") then
          enorm=fmass(1)*vnorm*vnorm/ergtkev*0.5
        else
          enorm=fmass(1)*clite2*(sqrt(1.+vnorm**2/clite2)-1.)/ergtkev
        endif
      else
        call diagwrng(1)
      endif
      vnorm2=vnorm*vnorm
      vnorm3=vnorm2*vnorm
      vnorm4=vnorm2*vnorm2
      cnorm=clight/vnorm
      cnorm2=cnorm**2
      cnorm3=cnorm**3
      if (relativ.eq."disabled") then
        cnorm2i=0. !YuP[07-2016] Bad way to control relativ.eq."disabled"
        cnormi=0.  !YuP[07-2016] Bad way to control relativ.eq."disabled"
      else
        cnorm2i=1./cnorm2
        cnormi=1./cnorm
      endif

c.....................................................................
c     energy conversion factor...
c.....................................................................

CDIR$ NEXTSCALAR
      do 12 k=1,ntotal
        fions(k)=.5*fmass(k)*vnorm**2/ergtkev
 12   continue

c..................................................................
c     Determine some mass ratios and normalization constants.
c..................................................................

      alp=7.3e-3
      r0=2.82e-13
      gacon2=2.*alp/(r0*fmass(kelec)*clight)
      do 10 i=1,ntotal
        do 11 k=1,ntotal
          gamt(i,k)=fmass(i)*fmass(k)*gacon2/(fmass(i)+fmass(k))
          !YuP[2019-07-25] Just a comment: This factor, m(i)*m(k)/(m(i)+m(k))
          !is from Killeen book, Eq.(2.1.5). 
          !For e-on-e Coulomb logarithm,
          !it will give me*me/(2*me) = me/2.
          !For e-on-ion, it will give me*mi/mi = me.
          ! The difference between these two types of interaction is
          ! a factor of 2, and ln(2)=0.7, 
          ! so that ln(L_ei) is larger than ln(L_ee) by 0.7.
          ! In NRL Formulary, and in Hesslow-JPP2018, they are almost same.
          ! It seems this factor, m(i)*m(k)/(m(i)+m(k)), 
          ! was only used in the Killeen book. The origin is unclear.
          !We can get a match with NRL if we use me/2 instead of me
          !in equation for ln(L_ei)
          ! In such a case,gamt(i,k)= (fmass(kelec)/2)*gacon2 = 
          ! = alp/(r0*clight)= 0.86348152303068
 11     continue
 10   continue
      do 40 k=1,ngen
        gam1=4.*pi*(charge*bnumb(k))**4/fmass(k)**2
        tnorm(k)=vnorm**3/(gam1*one_)
 40   continue
      return
      end
