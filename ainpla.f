c
c
c
      subroutine ainpla
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..........................................................
c     This routine initialize some plasma parameter profiles
c     on the whole lrzmax radial mesh
c     (some where defined in diaggnde before)
c.............................................................

      include 'param.h'
      include 'comm.h'

c.......................................................................
cl    1. Energy, v-thermal
c.......................................................................

cl    1.1 radial mesh

      do 110 k=1,ntotal
        rstmss=fmass(k)*clite2/ergtkev
        do 111 l=1,lrzmax
          thta=rstmss/temp(k,l)
          if (thta.gt.100. .or. relativ.eq."disabled") then
            energy(k,l)=1.5*temp(k,l)
            !write(*,*)'ainpla.1:  energy(k,l)/1.5=',energy(k,l)/1.5
          else
            call cfpmodbe(thta,bk1,bk2)
            energy(k,l)=rstmss*(bk1/bk2-1.+3./thta)
            !write(*,*)'ainpla.2:  energy(k,l)/1.5=',energy(k,l)/1.5
          endif
c          write(*,*)'ainpla: k,l,energy(k,l):',k,l,energy(k,l)
          vth(k,l)=((temp(k,l)*ergtkev)/fmass(k))**.5
          if (k .eq. kelec) vthe(l)=vth(kelec,l)
 111    continue
 110  continue

c.......................................................................
cl    1.2 parallel mesh
c.......................................................................

      if (cqlpmod .eq. "enabled") then

        do 120 k=1,ntotal
          rstmss=fmass(k)*clite2/ergtkev
          do 121 l=1,lsmax
            thta=rstmss/temppar(k,l)
            if (thta.gt.100. .or. relativ.eq."disabled") then
              enrgypa(k,l)=1.5*temppar(k,l)
            else
              call cfpmodbe(thta,bk1,bk2)
              enrgypa(k,l)=rstmss*(bk1/bk2-1.+3./thta)
            endif
            vthpar(k,l)=((temppar(k,l)*ergtkev)/fmass(k))**.5
 121      continue
          if (sbdry.eq."periodic" .and. transp.eq."enabled") then
            enrgypa(k,0)=enrgypa(k,lsmax)
            enrgypa(k,lsmax+1)=enrgypa(k,1)
            vthpar(k,0)=vthpar(k,lsmax)
            vthpar(k,lsmax+1)=vthpar(k,1)
          endif
 120    continue

      endif
c.......................................................................
c     2. Compute radial Z-effective
c.......................................................................

      if (izeff.eq."ion") then
        k1=ngen+1
      else
        k1=1
      endif
      do 200 l=1,lrzmax
        zeff(l)=0.
        zeff1=0.
        zeff4(l)=0.d0 !Yup[2014-05-27] Initialize to 0.
        xq=0.
        do 210 k=k1,ntotal
          if (k.eq.kelecg .or. k.eq.kelecm) goto 210
cBobH990128          if (k.eq.izeff) goto 210
          xq=xq+1.
          zeff(l)=zeff(l)+bnumb(k)**2*reden(k,l)
          zeff4(l)=bnumb(k)**4*reden(k,l)+zeff4(l)
          zeff1=zeff1+bnumb(k)*reden(k,l)
 210    continue
        zeff4(l)=zeff4(l)/xq
        zeff(l)=zeff(l)/zeff1
 200  continue

      return
      end
