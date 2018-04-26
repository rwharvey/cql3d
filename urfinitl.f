c
c
      subroutine urfinitl
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     This routine does some post-namelist-read initialization
c     and checking of input.
c..................................................................

      include 'param.h'
      include 'comm.h'

CMPIINSERT_INCLUDE

CMPIINSERT_IF_RANK_EQ_0
      ! make plots on mpirank.eq.0 only
      if (noplots.ne."enabled1") then
      write(t_,1000) 
 1000 format("Urf (lower hybrid, fast wave, ech, ebw...) parameters:")
      CALL PGMTXT('T',-7.,0.,0.,t_)

c      write(t_,1001) nrayn,nrayelts
 1001 format("====>NRAYn =",i5,"     ====>NRAYELTs = ",i5)
c      CALL PGMTXT('T',-8.,0.,0.,t_)

      write(t_,1002) nmodsa
 1002 format("====>NMODSA = ", i3)
      CALL PGMTXT('T',-9.,0.,0.,t_)
      endif
CMPIINSERT_ENDIF_RANK

c..................................................................
c     Return if urfmod.eq."disabled"
c..................................................................

      if(urfmod.eq."disabled") return


c..................................................................
c     If eqsource="tsc", then ensure that xbr [Prior to genray]
c     reads extended eqdsk:
c..................................................................

      if(eqsource.eq."tsc")  ieqbrurf=4

c..................................................................
c     Counter for calculation of urf-diffusion coefficients:
c..................................................................
      nurf=0

c.......................................................................
c     Re-run a few flux surfaces => "do not update delpwr" option
c.......................................................................
      if (urfrstrt .eq. "enabled") then
        nrfpwr=0
        nrfitr1=0
        nrfitr2=0
        nrfitr3=0
      endif

c.......................................................................
c     Check some input values
c.......................................................................
      
      if (nurftime.gt.nbctimea) stop "nurftime.gt.nbctimea"

c.......................................................................
c     Check that rftype() or lh/fw/ech wave-type designators are used,
c     but not both.
c     irftype=0, older method, lh/fw/ech wage-type designation
c             1, new method, via rftype(1:nmodsa).ne.'notset'
c.......................................................................

      irftype=0
      do i=1,nmodsa
         if (rftype(i).ne."notset" .and. irftype.eq.0) irftype=1
      enddo
      if ( (lh.eq."enabled" .or. fw.eq."enabled" .or.
     1     ech.eq."enabled") .and. irftype.eq.1) then
         WRITE(*,*)"STOP: Can't set both lh/fw/ech and rftype()"
         STOP
      endif

c.......................................................................
c     If irftype=1, rfread="netcdf" is only choice
c     (If irftype=0, rfread can be either "text" or "netcdf".)
c.......................................................................

      if (irftype.eq.1 .and. rfread.ne."netcdf") then
         write(*,*)
         WRITE(*,*)"urfinitl STOP: Incompatible rftype and rfread"
         STOP
         write(*,*)
      endif

c.......................................................................
c     Count the number of rf types, mrf
c.......................................................................

      mrf=0
      if ( irftype.eq.0 ) then 
         if (lh.eq."enabled") then
            mrf=mrf+1
            rftype(mrf)="lh"
         endif
         if (ech.eq."enabled") then
            mrf=mrf+1
            rftype(mrf)="ech"
         endif
         if (fw.eq."enabled") then
            mrf=mrf+1
            rftype(mrf)="fw"
         endif
      endif
      
      if ( irftype.eq.1 ) then
         do k=1,nmodsa
            if (rftype(k).ne."notset") mrf=mrf+1
            if (rftype(k).eq."notset") go to 5
         enddo
 5       continue
      endif

c     Make sure not attempting to damp/ql diffuse species which
c     is not present.
      do krf=1,mrf
         if (nrfspecies(krf).gt.ngen) STOP 'nrfspecies().gt.ngen'
      enddo

c...................................................................
c     For each wave type:
c     nharms() is the number of cyclotron harmonics calculated,
c     starting at nharm1().
c     (nharms.gt.1 .and. mrf.gt.1) is now permitted by the
c     storage scheme [BH060314].

c     mrfn is the number of wave "modes", that is, the sum over
c     wave types of the number of harmonics for each wave type,
c     i.e., each harmonic for each wave type is counted as a
c     separate "mode" or interation with the plasma distribution(s).
c
c     NOTE: It is possible that mrfn may be reset to a different
c           value in urfread, based on reading old rf input files
c           which had nharms set = 0. This perhaps justifies some
c           repeat coding in urfsetup [BH080919, but I think this
c           is completely obsolete.]
c...................................................................

      mrfn=0
      do k=1,mrf
         if (nharms(k).eq.0) then
            mrfn=mrfn+1   !  I.E., each rf type implies 
                          !  at least 1 harmonic.
         else
            mrfn=mrfn+nharms(k)
         endif
      enddo

      if (mrfn.gt.nmodsa) then
         write(*,*)'urfinitl: mrfn>nmodsa.  mrfn,nmodsa=',mrfn,nmodsa
         STOP 'Increase nmodsa.'
      endif

c.......................................................................
c     Set up URF module file names, for rfread="netcdf".
c     rffile(1:nmodsa) are set by default to "notset",
c     and may be set in the namelist input.
         
c     If rffile(1)="mnemonic", rffile(1:3) is based
c     on the namelist input character variable mnemonic.
c.......................................................................
      
      if (rfread.eq."netcdf") then
         
c     irftype=0 case:
         krf=0
         if (ech.eq."enabled") then
            krf=krf+1
            if (rffile(krf).eq."notset") rffile(krf)="rayech.nc"
         endif
         if (fw.eq."enabled") then
            krf=krf+1
            if (rffile(krf).eq."notset") rffile(krf)="rayfw.nc"
         endif
         if (lh.eq."enabled") then
            krf=krf+1
            if (rffile(krf).eq."notset") rffile(krf)="raylh.nc"
         endif
         
         if (rffile(1).eq."mnemonic") then
            write(t_,110) mnemonic(1:length_char(mnemonic))
 110        format(a,"_rf.nc")
            rffile(1)=t_
            
c     Here we use achar(48+1)="1", etc., to name input files
            do i=1,mrf-1
               if (mrf.gt.10) then
                  WRITE(*,*)'urfinitl:  Expand calc of file name'
                  STOP
               endif
               write(t_,111) mnemonic(1:length_char(mnemonic)),i
 111           format(a,"_rf.",i1,".nc")
               rffile(i)=t_
               
               write(*,*)'urfinitl: i,rffile(i) =',i,rffile(i)
            enddo
         endif
      endif  !On rfread
      
      return
      end
