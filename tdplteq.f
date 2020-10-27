c
c
      subroutine tdplteq(krf) !YuP[2020-10-19] now it is called for each wave type (was: each mode)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

      character*8 textt
      common /plttextt/ textt(200)
      character*8 pltrays
      character*8 plteq

      REAL*4 rbot,zbot,rtop,ztop
      REAL*4 RTAB1(LFIELDA),RTAB2(LFIELDA)
      REAL*4 PGER1,PGERNNR,PGZERO,PGZMIN,PGZMAX,PGEZNNZ
      REAL*4 wk_ray_r(nrayelts), wk_ray_z(nrayelts)

      REAL*4 :: R40=0.
      REAL*4 :: R4P15=.15,R4P85=.85,R4P9=.9

c..................................................................
c     This routine plots out the contours (flux surfaces) on
c     which the calculations are done.
c..................................................................

CMPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only
 
      pltrays="enabled" ! for plotting rays over cross-sectional view.
                        ! Can be moved to namelist later.

      if (noplots.eq."enabled1") return
      plteq="enabled"
      if (plteq.eq."disabled") return

c     Set up textt to contain number 1 to 200:
      call micfrplt
c     Would be nice to label the flux surfaces, but
c     it isn't done at the moment....

      if (eqsym.ne."none") then ! Up-dn symm
         ztop=.75*ez(nnz)/(er(nnr)-er(1))+.05
      else ! Not up-dn symm
         ztop=0.95
      endif
      
      ! YuP[03-2016][07-2017] Added plotting rays in cross-sectional view
      if (urfmod.ne."disabled" .and. pltrays.eq.'enabled') then 
         ! over-write setting for page size - plot whole cross-section.
         ! For up-dn symmetrical case only half of surfaces are plotted
         ! but rays could be in the other hemisphere.
         ztop=.95
      endif

      
      ! YuP: min and max value of Z-coord over all flux surfaces:
      solz_min=MINVAL(solz) 
      solz_max=MAXVAL(solz)
      PGZMIN=real(solz_min) 
      PGZMAX=real(solz_max) 
      
      CALL PGPAGE
      CALL PGSVP(R4P15,R4P85,R4P15,ztop)

      PGER1=er(1)
      PGERNNR=er(nnr)
      PGZERO=0.
      PGEZNNZ=ez(nnz)
      !In case of eqmod='disabled' solr() and solz() could be empty.
      !write(*,*)'tdplteq:',PGER1,PGERNNR,PGZMIN,PGZMAX 
      !Actually, with psimodel='spline', they are set now !YuP[2020-01-29]
      if( (eqsym.ne."none") .and.
     +    (urfmod.eq."disabled" .or. pltrays.eq.'disabled') ) then 
         !plot half only
c-YuP         CALL PGSWIN(PGER1,PGERNNR,PGZERO,PGEZNNZ)
c-YuP         CALL PGWNAD(PGER1,PGERNNR,PGZERO,PGEZNNZ)
         CALL PGSWIN(PGER1,PGERNNR,PGZMIN,PGZMAX)
         CALL PGWNAD(PGER1,PGERNNR,PGZMIN,PGZMAX)
      else !eqsym=none; and/or urfmod='enabled',
         ! plot upper and lower halves:
         CALL PGSWIN(PGER1,PGERNNR,-PGEZNNZ,PGEZNNZ)
         CALL PGWNAD(PGER1,PGERNNR,-PGEZNNZ,PGEZNNZ)
      endif
     
      CALL PGBOX('BCNST',R40,0,'BCNST',R40,0)
      if ( (urfmod.ne."disabled") .and. (pltrays.eq.'enabled') 
     +                            .and. (krf.gt.0)            ) then
         CALL PGLAB('Major radius (cms)','Vert height (cms)',
     +        'Fokker-Planck Flux Surfaces + Rays')
      else ! urfmod='disabled', or krf=0
         CALL PGLAB('Major radius (cms)','Vert height (cms)',
     +        'Fokker-Planck Flux Surfaces')
      endif

      IF (LRZMAX.GT.200) STOP 'TDPLTEQ: CHECK DIM OF TEXTT'

      do 10 l=1,lrzmax
         IF (LORBIT(L).GT.LFIELDA) STOP'TDPLTEQ: CHECK DIM OF RTAB1/2'
        do 20 j=1,lorbit(l)
           RTAB1(j)=solr(lorbit(l)+1-j,l)
           RTAB2(j)=solz(lorbit(l)+1-j,l)
 20     continue
        text(1)=textt(l)
        CALL PGLINE(LORBIT(L),RTAB1,RTAB2)
        ! YuP[03-2016] Added plotting rays in cross-sectional view
        if ( (urfmod.ne."disabled") .and. (pltrays.eq.'enabled')
     +       .and. (eqsym.ne.'none')) then 
         ! Add surfaces in whole cross-section.
         ! For up-dn symmetrical case only half of surfaces are plotted
         ! but rays could be in the other hemisphere, 
         ! so plot the other half:
         CALL PGLINE(LORBIT(L),RTAB1,-RTAB2)
        endif
 10   continue

      !if(eqmod.eq."enabled")then
      if((eqsym.eq.'none').or. 
     + ((urfmod.ne."disabled") .and. (pltrays.eq.'enabled')) ) then
     
        if(ncontr.gt.1) then
          ! YuP[2015/05/03] Add LCFS, if available
          ncontr_= min(ncontr,LFIELDA)
          do ilim=1,ncontr_
             RTAB1(ilim)=rcontr(ilim)
             RTAB2(ilim)=zcontr(ilim)
          enddo
          CALL PGSCI(3) !green color
          CALL PGLINE(ncontr_,RTAB1,RTAB2)
          CALL PGSCI(1) !black color restore
        endif
      
        if(nlimiter.gt.1) then
          ! YuP[2019-02-22] Add "last surface" (plasma border), if available
          nline= min(nlimiter,LFIELDA)
          do ilim=1,nline
             RTAB1(ilim)=rlimiter(ilim)
             RTAB2(ilim)=zlimiter(ilim)
          enddo
          CALL PGSCI(4) !blue color
          CALL PGSLW(lnwidth*2) ! bold
          CALL PGLINE(nline,RTAB1,RTAB2)
          if(machine.eq."mirror") then
          CALL PGLINE(nline,-RTAB1,RTAB2) !mirror area to the left of Z-axis
          endif
          CALL PGSLW(lnwidth) ! restore
          CALL PGSCI(1) !black color restore
        endif
        
      endif !eqsym.eq.'none' or (urfmod.ne."disabled")&(pltrays.eq.'enabled')
      !endif !(eqmod.eq."enabled")
      
      
      
c..................................................................
c YuP[03-2016] Added plotting rays in cross-sectional view
      if ((urfmod.ne."disabled") .and. (pltrays.eq.'enabled') .and. 
     +    (krf.gt.0)  ) then
        nray_krf= nray(irfn(krf)) !YuP[2020-09-23] YuP[2020-10-19] 
        !Such arrays as freqcy, delpwr, cwexde, etc, -
        !they are not functions of krf (wave type; Example: krf=1:mrf=1:2),
        !but rather they are functions of krfn (wave modes; Example: mrfn=3+3).
        !Initially, they are read for each krf, 
        !but then duplicated into each krfn index.
        !This is done in subr.urfread, at the end.
        !Here, tdplteq(krf) is called for each krf type (after 2020-10-19).
        !So,  krf is wave type (=1 and 2 in our example).
        !Then, we need to use data as
        ! wr(is,iray, irfn(krf)),
        !where irfn(krf=1)=1 and irfn(krf=2)=4 in our example,
        !irfn(krf) is pointing to beginning harmonic of each wave type
        !== index (in 1:mrfn) of the lowest harmonic for each wave type 
        do iray=1,nray_krf  !Loop over rays, for a given wave type krf
           nrayelt00=nrayelt(iray,irfn(krf))
c           write(*,*)'tdplteq: krf, iray, nrayelt00=',
c     +                         krf, iray, nrayelt00
c           write(*,'(a,3i6)')
c     +      'tdplteq: iray,lloc(nrayelt00),llray(nrayelt00)=',
c     +      iray,lloc(nrayelt00,iray,krf),llray(nrayelt00,iray,krf) !local to ray element point.
           wk_ray_r=0.0 ! reset
           wk_ray_z=0.0 ! reset
           do is=1,nrayelt00
             wk_ray_r(is)=wr(is,iray,irfn(krf)) !YuP[2020-10-19] corrected
             wk_ray_z(is)=wz(is,iray,irfn(krf)) !YuP[2020-10-19] corrected
           enddo
c           write(*,*)'tdplteq: krf, minval(wk_ray_r),maxval(wk_ray_r)=',
c     +                         krf, minval(wk_ray_r),maxval(wk_ray_r)
           CALL PGLINE(nrayelt00,wk_ray_r,wk_ray_z)
        enddo  ! iray 
        !Add some info on rays:
        write(t_,'(a,i3, a,2i3, a,e12.5)') 'krf type=',krf, 
     +                      '  nharm1,nharms=', nharm1(krf),nharms(krf),
     +                      '    f[Hz]=', freqcy(irfn(krf))
        !YuP: Careful! nharm1() and nharms() were not duplicated into each mode,
        !so they are still a function of krf (wave type); 
        !But freqcy was duplicated, so use irfn(krf).
        CALL PGMTXT('T',R4P9,R40,R40,t_) ! 0.9=just outside of viewport
        !enddo
      endif  
c..................................................................

      return
      end
