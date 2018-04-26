c
c
      subroutine urfdampa(krf)
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     This routine computes the collisional or additional linear
c     damping, depending on iurfcoll and iurfl, for all
c     ray elements on all rays for UNIT power delpwr. 
c     It only needs to be done once.
c..................................................................

      include 'param.h'
      include 'comm.h'

      character*8 ifirst
      save ifirst,nray0
      data ifirst/"first"/
      data nray0 /1/


c-YuP      if (ifirst.ne."first") go to 110 ! urfpwrl and urfpwrc are defined
c-YuP      ifirst="notfirst"  ! only for the first called krf (usually krf=1)
c      if (nharms.gt.0 .and. krf.gt.1) go to 110

c..................................................................
c     Loop over rays
c..................................................................

      do 10 iray=nray0,nray(krf)

c..................................................................
c     Loop of ray elements
c..................................................................

        do 20 is=1,nrayelt(iray,krf)

c..................................................................
c     Determine the length of the ray element, slngth
c..................................................................

          slngth=0.0
          if (is.eq.1) then
            slngth=ws(1,iray,krf)
          elseif (is.lt.nrayelt(iray,krf)) then
            slngth=0.5*(ws(is+1,iray,krf)-ws(is-1,iray,krf))
          elseif (is.eq.nrayelt(iray,krf)) then
            slngth=0.5*(ws(is,iray,krf)-ws(is-1,iray,krf))
          endif

c.......................................................................
c     Damping due to collisional absorption of the wave:
c.......................................................................

          if (iurfcoll(krfn(krf)).eq."enabled")  then
            urfpwrc(is,iray,krf)=1.d0-exp(-salphac(is,iray,krf)*slngth)
          endif

c.......................................................................
c     damping due to additional linear absorption (e.g., ions)
c     passed from rayop file:
c.......................................................................

          if (iurfl(krfn(krf)).ne."disabled") then !"enabled" or "split"
            urfpwrl(is,iray,krf)=1.d0-exp(-salphal(is,iray,krf)*slngth)
          endif

 20     continue
 10   continue
      iray=1

 110  continue

      return
      end

