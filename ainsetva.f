c
      subroutine ainsetva
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c.......................................................................
c     This routine sets auxiliary values which depend on
c     input variables.
c     Should be called after reading the namelist setup 
c     (not after setup0).
c     Also, it checks some namelist variables for consistency, 
c       and resets some if necessary.
c.......................................................................

      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE
c
CMPIINSERT_IF_RANK_EQ_0
      WRITE(*,*)''
      WRITE(*,*)'In ainsetva'
c
c.......................................................................
c     0. Check for consistency of some parameter settings
c        (in param.h).
c......................................................................

      if (lrorsa.ne.max(lrza,lsa)) 
     +     stop 'check consistency of lrorsa in param.h'

cBH081106:Found following condition leads to overwrite.
cBH081106:For cqlpmod.ne.enabled, could check code to see how
cBH081106:essential the condition is.

      if (lsa.lt.lrza) stop 'check that lsa.ge.lrza?'

      if (cqlpmod.eq."enabled" .and. lza.lt.lsa) 
     +     stop 'check consistency of lza,lsa for cqlpmod=enabled'

cBH060314      if (nmodsa.ne.3)
cBH060314     +     stop 'Better check out code for nmodsa.ne.3'

      if (nefitera.lt.2) stop 'Coding assumes nefitera.ge.2'

c.......................................................................
cl    1. check and define some mesh values
c.......................................................................

c     cannot run CQL3D with lrzmax>lrz (only CQL or CQLP)
      if (lrzdiff.eq."enabled" .and. transp.eq."enabled" .and.
     +  cqlpmod.ne."enabled") call diagwrng(17)
      if (lsdiff.eq."enabled" .and.  transp.eq."enabled" .and.
     +  cqlpmod.eq."enabled") call wpwrng(4)
      if (lrzdiff.eq."enabled" .and. partner.eq."selene") WRITE(*,
     +  '("WARNING: partner=selene and lrzdiff=enabled not checked")')

      if (urfmod.ne."disabled" .and. meshy.eq."fixed_mu")
     +  call diagwrng(19)
      if (urfmod.ne."disabled" .and. iy.gt.255)
     +  call diagwrng(22)

c     storage problem for calc of parallel distn function, if
c     lrz.gt.lz
      if ((pltprpp.eq."enabled" .or.  knockon.ne."disabled") .and.
     +     lrz.gt.lz) then
         WRITE(*,*) 'Need lrz.le.lz for pltprpp or knockon enabled'
         STOP
      endif

CMPIINSERT_ENDIF_RANK

      if ((iy/2)*2 .ne. iy) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,220)
CMPIINSERT_ENDIF_RANK
 220     format(//,'WARNING:  Making iy even, increasing by 1',//)
         iy=iy+1
      endif

c%OS  if (iy .lt. 16) iy=16
c%OS  if (jx .lt. 8) jx=8
      iyh=iy/2
      jxm1=jx-1
      jxp1=jx+1
      iyp1=iy+1
      iyjx=iy*jx
      iyjxp1=iy*(jx+1)
      iyp1jx=(iy+1)*jx
      iyjx2=(iy+2)*(jx+2)


      if (ndeltarho.ne."disabled". and. (nt_delta/2)*2.ne.nt_delta) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,223)
CMPIINSERT_ENDIF_RANK
 223     format(//,'WARNING:  Making nt_delta even, increasing by 1',//)
         nt_delta=nt_delta+1
      endif

c     fr_gyrop='disabled' is default, for backwards compatability.
c     Warn on use of fr_gyro for ZOW particles.
      if (frmodp.eq.'enabled') then
         if (fr_gyrop.ne.'disabled') then
CMPIINSERT_IF_RANK_EQ_0
            WRITE(*,*)
            WRITE(*,*)'WARNING: ZOW NB ions, but using fr_gyro=enabled'
            WRITE(*,*)
CMPIINSERT_ENDIF_RANK
         endif
      endif

      if (ndeltarho.ne."disabled") then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,224)
CMPIINSERT_ENDIF_RANK
 224     format(//,'WARNING: ndeltarho.ne."disabled"',/
     +             '         Remember that baviorbt (taunew=enabled)',/
     +             '         is setup with nii=1. Consider nii.gt.1',//)
      endif

      if (ndeltarho.ne."disabled". and. taunew.ne."enabled") then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,225)
CMPIINSERT_ENDIF_RANK
 225     format(//,'ERROR:  Presently ndeltarho calls are only setup',/
     +             '        for taunew=enabled',//)
         STOP
      endif

      if (ndeltarho.ne."disabled". and. lrzdiff.ne."disabled") then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,226)
CMPIINSERT_ENDIF_RANK
 226     format(//,'ERROR:  ndeltarho feature not setup for setup',/
     +             '        for lrzdiff.ne."disabled"',//)
         STOP
      endif

      if (ndeltarho.ne."disabled". and. urfdmp.eq."secondd") then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,227)
 227     format(//,'ERROR:  Presently ndeltarho calls are only setup',/
     +             '        for urfdmp.ne.disabled',//)
         WRITE(*,228)
CMPIINSERT_ENDIF_RANK
 228     format(//,'ERROR:  Resetting urfdmp="firstd"',/ 
     +             '         See cqlinput_help',//)
         urfdmp="firstd"
      endif



      if (ngen.gt.ngena .or. nmax.gt.nmaxa) call diagwrng(12)
      if (lz.gt.lza) then
         lz=lza
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'WARNING: lz reset to lza =', lza
CMPIINSERT_ENDIF_RANK
      endif
cBH110330      if (msxr.gt.mx) then
cBH110330        msxr=mx
cBH110330CMPIINSERT_IF_RANK_EQ_0
cBH110330        WRITE(*,201) mx 
cBH110330CMPIINSERT_ENDIF_RANK
cBH110330 201    format('WARNING: msxr reset to mx =', i5)
cBH110330      endif
cBH110331      if (mmsv.gt.mx) then
cBH110331        mmsv=mx
cBH110331CMPIINSERT_IF_RANK_EQ_0
cBH110331        WRITE(*,202) mx
cBH110331CMPIINSERT_ENDIF_RANK
cBH110331 202    format('WARNING: mmsv reset to mx =', i5)
cBH110331      endif
      mxp1=mx+1
      if (njene.gt.njenea) call diagwrng(16)

      do 100 ll=1,lrors
        iy_(ll)=iy
        iyh_(ll)=iyh
        iyjx_(ll)=iyjx
 100  continue
      iymax=iy

      ipxy=min(51,iy)
      jpxy=min(101,jx+1)
      if (mod(jpxy,2).eq.0) jpxy=jpxy-1

c.......................................................................
c     1.1 Check mx so no overwrite because of
c        tamt1ptr=temp1ptr, tamt2ptr=temp4ptr.
c        (Overwrite has already occured, but better to warn late
c        than not at all).
c        Similarly, will have overwrite if
c        jpxy*ipxy.gt.(iyp1+1)*(jxp1+1).
c.......................................................................

      !-YuP: Why needed ???      COMMENTING it out:
ccc      if(2*jx*(mx+5)*(mx+5).gt.(iy+2)*(jx+2)) then
ccc      write(*,*)'ainsetva: Possible tamt1 and tamt2 overwrite',jx,iy,mx
ccc         stop 'problem with tamt1 and tamt2 overwrite'
ccc      endif

      if (jpxy*ipxy .gt.  iy*(jx+1)  .or.
     +    jpxy*ipxy .gt. (iy+1)*jx   ) then
         stop 'problem with jpxy*ipxy dimensioning'
      endif
c


c.......................................................................
cl    2. define effective time-step according to model
c.......................................................................

      do 50 i=1,ndtr1a
      if (nondtr1(i).ge.0) then
        if (mod(nondtr1(i),nrstrt).ne.0) then
          nondtr1(i)=(mod(nondtr1(i),nrstrt)+1)*nrstrt
CMPIINSERT_IF_RANK_EQ_0
          WRITE(*,101) nondtr1(i),i
CMPIINSERT_ENDIF_RANK
 101      format('WARNING: nondtr1(i) reset to ',i5,'  i=',i2)
        endif
      endif
      if (nondtr1(i).eq.n) dtr=dtr1(i)
 50   continue

      dtreff=dtr
      if (transp.eq."enabled") then
        dttr=dtr*nrstrt
        if (adimeth .eq. "enabled") then
          dtreff=dtr
          nrstrt=1
          dttr=dtreff
        endif
        do k=1,ngen
           if (difus_type(k).eq."neo_smpl" .and. k.eq.kelecg) then
CMPIINSERT_IF_RANK_EQ_0
              WRITE(*,*)"STOP: Not sensible to use ion neocl for e's"
CMPIINSERT_ENDIF_RANK
              STOP
           endif
        enddo
      endif

c     Time-dep radial transport multipliers
      do k=1,ngena
         drrt(k)=one
         drt(k)=one
      enddo

      if (transp .ne. "enabled") nontran=0
      if (transp.eq."enabled" .and. adimeth.eq."enabled") then
        if (nonadi .lt. nontran) nonadi=nontran
      endif

      if (ndifus_io_t.gt.nbctimea) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)"STOP: ndifus_io_t.gt.nbctimea"
CMPIINSERT_ENDIF_RANK
         STOP
      endif


c.......................................................................
cl    3. Check and adjust species if iprozeff.ne."disabled"
c.......................................................................

      if (iprozeff.ne."disabled") then  !endif at line 285
        !  iprozeff = 
        != "parabola","prbola-t", "spline"  or "spline-t", then this 
        ! option gives ion densities so as to achieve a given zeff
        ! radial profile, consistent with charge neutrality and the 
        ! given electron density profile.

c       Check for electrons
        if (kelec.eq.0) stop 'ainsetva: need electron species'
c       Check for at least one maxwellian ion.  For the
c       case of only general species ions, could adjust code slightly.
        if (nionm.eq.0) 
     +     stop 'ainsetva: need at least one maxwellian ion species' 

c
c     Check number of ion Maxwl species with different bnumb.
        if (nionm.lt.1) stop 'ainsetva: ion species problem: nionm.lt.1'
        ndif_bnumb=1
        do k=2,nionm
           if (abs(bnumb(kionm(k))/bnumb(kionm(1))-1.).gt.0.01) 
     +          ndif_bnumb=ndif_bnumb+1
        enddo
           


        if (ndif_bnumb.eq.1) then
c         add a species
          nionm=nionm+1
          if (nionm.gt.nmaxa) stop 'ainsetva:  check nmaxa'
          kionm(nionm)=kionm(nionm-1)+1
          ntotal=ntotal+1
          nmax=nmax+1
          if (kelecm.gt.kionm(1)) then   ! Move up electron data
            kelecm=kelecm+1
            fmass(kelecm)=fmass(kelecm-1)
            bnumb(kelecm)=bnumb(kelecm-1)
            kspeci(1,kelecm)=kspeci(1,kelecm-1)
            kspeci(2,kelecm)=kspeci(2,kelecm-1)
            reden(kelecm,0)=reden(kelecm-1,0)
            reden(kelecm,1)=reden(kelecm-1,1)
            do 60  ll=1,njene
              enein(ll,kelecm)=enein(ll,kelecm-1)
 60         continue
            temp(kelecm,0)=temp(kelecm-1,0)
            temp(kelecm,1)=temp(kelecm-1,1)
            if (nbctime.ne.0) then
               do jtm=1,nbctime
                  redenc(jtm,kelecm)=redenc(jtm,kelecm-1)
                  redenb(jtm,kelecm)=redenb(jtm,kelecm-1)
                  tempc(jtm,kelecm)=tempc(jtm,kelecm-1)
                  tempb(jtm,kelecm)=tempb(jtm,kelecm-1)
               enddo
            endif
          endif
          if (kelec.gt.kionm(1))  kelec=kelec+1
          fmass(kionm(nionm))=2.*50.*1.676e-24
          bnumb(kionm(nionm))=50.
          kspeci(1,kionm(nionm))="impurity"
          kspeci(2,kionm(nionm))="maxwell"
          if (nbctime.ne.0) then
             do jtm=1,nbctime
                redenc(jtm,kionm(nionm))=1.    !Nominal value, to
                redenb(jtm,kionm(nionm))=1.    !be reset with zeff.
                tempc(jtm,kionm(nionm))=tempc(jtm,kionm(nionm)-1)
                tempb(jtm,kionm(nionm))=tempb(jtm,kionm(nionm)-1)
             enddo
          endif
CMPIINSERT_IF_RANK_EQ_0
          WRITE(*,*)'ainsetva: Impurity species added, k=',k
          WRITE(*,*)'ainsetva: bnumb=',bnumb(kionm(nionm))
          WRITE(*,*)'ainsetva: fmass/proton_mass=',
     +                         fmass(kionm(nionm))/proton
CMPIINSERT_ENDIF_RANK
        elseif (ndif_bnumb.eq.2) then
          continue  
        else ! ndif_bnumb ne.1, ne.2
          WRITE(*,*)'ainsetva: iprozeff=', iprozeff
          WRITE(*,*)'ainsetva: ndif_bnumb=',ndif_bnumb
          stop 'ainsetva: ion species problem. ndif_bnumb'
        endif

      endif  !On iprozeff.ne."disabled"

c.......................................................................
c     Reset parabolic profile specifiers for initial profile calc
c     if time-dependent profiles specified
c     BH171124: Actually, should now be using iprone="prbola-t" in
c               nbctime.ne.0 cases where time dependent profiles are
c               desired.  Will ask that cqlinput be checked.
c.......................................................................
      if (nbctime.ne.0) then
         if (iprone.eq."parabola" .or. iprote.eq."parabola"
     1      .or.iproti.eq."parabola" .or. iprozeff.eq."parabola"
     1      .or.iproelec.eq."parabola" .or. iprovphi.eq."parabola"
     1      .or. iprocur.eq."parabola"   ) then
            !YuP[2018-01-02] Added iprocur in the above if()then
CMPIINSERT_IF_RANK_EQ_0
            WRITE(*,*)'------------------------------------------------'
            WRITE(*,*) 'ainsetva: nbctime=',nbctime
            WRITE(*,*) 'If nbctime.ne.0: CHECK ipro*** settings.'
            WRITE(*,*) 'Some of them (or all) are set to "parabola" .'
            WRITE(*,*) 'Resetting those to time-dependent "prbola-t" '
            WRITE(*,*)
CMPIINSERT_ENDIF_RANK
            ! BH,YuP[2018-01-02] Added resetting of ipro*** values 
            ! to "prbola-t", in case when the run 
            ! is done with nbctime>0 (time-dependent profiles),
            ! but ipro*** values are set to "parabola". 
            if(iprone.eq."parabola") then
               iprone="prbola-t"
CMPIINSERT_IF_RANK_EQ_0
               WRITE(*,*)'iprone is reset to "prbola-t" '
CMPIINSERT_ENDIF_RANK
            endif
            if(iprote.eq."parabola") then
               iprote="prbola-t"
CMPIINSERT_IF_RANK_EQ_0
               WRITE(*,*)'iprote is reset to "prbola-t" '
CMPIINSERT_ENDIF_RANK
            endif
            if(iproti.eq."parabola") then
               iproti="prbola-t"
CMPIINSERT_IF_RANK_EQ_0
               WRITE(*,*)'iproti is reset to "prbola-t" '
CMPIINSERT_ENDIF_RANK
            endif
            if(iprozeff.eq."parabola") then
               iprozeff="prbola-t"
CMPIINSERT_IF_RANK_EQ_0
               WRITE(*,*)'iprozeff is reset to "prbola-t" '
CMPIINSERT_ENDIF_RANK
            endif
            if(iproelec.eq."parabola") then
               iproelec="prbola-t"
CMPIINSERT_IF_RANK_EQ_0
               WRITE(*,*)'iproelec is reset to "prbola-t" '
CMPIINSERT_ENDIF_RANK
            endif
            if(iprocur.eq."parabola") then
               iprocur="prbola-t"
CMPIINSERT_IF_RANK_EQ_0
               WRITE(*,*)'iprocur is reset to "prbola-t" '
CMPIINSERT_ENDIF_RANK
            endif
            if(iprovphi.eq."parabola") then
               iprovphi="prbola-t"
CMPIINSERT_IF_RANK_EQ_0
               WRITE(*,*)'iprovphi is reset to "prbola-t" '
CMPIINSERT_ENDIF_RANK
            endif
CMPIINSERT_IF_RANK_EQ_0
            WRITE(*,*)'------------------------------------------------'
CMPIINSERT_ENDIF_RANK
         endif
      endif
      
      if (nbctime.ne.0) then ! YuP[2018-01-02] Since ipro** are reset
        ! to "prbola-t" (above), this part is not needed anymore:
        if (iprone.eq.'parabola') then
           do k=1,ntotal
              reden(k,0)=redenc(1,k)
              reden(k,1)=redenb(1,k)
CMPIINSERT_IF_RANK_EQ_0
              WRITE(*,*)'ainsetva:k,reden(0:1,k)',k,reden(k,0:1)
CMPIINSERT_ENDIF_RANK
           enddo
        endif
        if (iprote.eq.'parabola') then
           do k=1,ntotal
              if (k.eq.kelecg .or. k.eq.kelecm) then
                 temp(k,0)=tempc(1,k)
                 temp(k,1)=tempb(1,k)
              endif
           enddo
        endif
        if (iproti.eq.'parabola') then
           do k=1,ntotal
              if (k.ne.kelecg .and. k.ne.kelecm) then
                 temp(k,0)=tempc(1,k)
                 temp(k,1)=tempb(1,k)
              endif
           enddo
        endif
        if (iprozeff.eq.'parabola') then
           zeffin(0)=zeffc(1)
           zeffin(1)=zeffb(1)
        endif
        if (iprovphi.eq.'parabola') then
           vphiplin(0)=vphic(1)
           vphiplin(1)=vphib(1)
        endif
        if (iproelec.eq.'parabola') then
           elecfld(0)=elecc(1)
           elecfld(1)=elecb(1)
        endif
      endif !  on nbctime.ne.0

c     Check only tmdmeth='method1', in case of "spline-t" profiles
      if (nbctime.gt.0) then
         if (      iprone.eq."spline-t"
     +      .or. iprote.eq."spline-t"
     +      .or. iproti.eq."spline-t"
     +      .or. iprozeff.eq."spline-t"
     +      .or. iproelec.eq."spline-t"
     +      .or. iprocur.eq."spline-t"
     +      .or. iprovphi.eq."spline-t") then
            if (tmdmeth.ne."method1") then
CMPIINSERT_IF_RANK_EQ_0
               WRITE(*,*)'Warning:  Re-setting tmdmeth=method1'
CMPIINSERT_ENDIF_RANK
               tmdmeth="method1"
            endif
         endif
      endif
      
c.......................................................................
cl    4. Numerical parameters
c
cl    4.1 Set non-valid flags and parameters related to CQLP
c.......................................................................

      if (cqlpmod.ne."enabled" .or. transp.ne."enabled") then
        sbdry = "disabled"
      endif
      if (cqlpmod .eq. "enabled") then
        if (lsmax .ge. 5) lz=lsmax
        numclas=nummods/10
        numindx=mod(mod(nummods,10),5)
        if (numclas.eq.1 .and. transp.eq."enabled") then
          if (updown .ne. "symmetry") call wpwrng(7)
          if (sbdry .ne. "periodic") call wpwrng(8)
          if ((ls/2)*2 .ne. ls) call wpwrng(9)
          if ((lz/2)*2 .ne. lz) call wpwrng(10)
          if (jx.lt.iy) stop 'check rhspar/bndmats storage in comm.h'
        endif
c     So far assumes lmidpln=1, otherwise, z, sz,psi, bmidpln, etc 
c     should be redefined
c     
        if (lmidpln .ne. 1) call wpwrng(3)
        do 310 ll=1,lrz
          lmdpln(ll)=lmidpln
 310    continue
c     parallel transport related flags
        if (transp .eq. "enabled") then
          if (nonelpr .lt. nontran) nonelpr=nontran
          if (noffelpr .gt. nofftran) noffelpr=nofftran
          if (numindx.eq.2 .and. numixts.ne.-1 .and. numixts.ne.1) 
     1      call wpwrng(16)
          if (numindx.eq.4 .and. lmidvel.ne.0) call wpwrng(17)
        endif

c     flags for options not yet valid with CQLP
        if (colmodl .eq. 4) call wpwrng(1)
        if (eoved .ne. 0.0) then
CMPIINSERT_IF_RANK_EQ_0
          PRINT *,' eoved changed to disabled'
CMPIINSERT_ENDIF_RANK
          eoved=0.0
        endif
        implct="enabled"
        machine="toroidal"
        npa_diag="disabled"
        softxry="disabled"
        symtrap="disabled"
        urfmod="disabled"
c$$$        vlfmod="disabled"
        vlhmod="disabled"
        lh="disabled"
        ech="disabled"
        fw="disabled"
        syncrad="disabled"
        bremsrad="disabled"
        qsineut="disabled"
        pltlos="disabled"
        nso=0
        pltso="disabled"
        do 315 k=1,ngen
          lossmode(k)="disabled"
          torloss(k)="disabled"
          do 316 l=0,lrz
            tauegy(k,l)=0.0
 316      continue
 315    continue

      endif

c.......................................................................
cl    4.2 Miscellaneous, Checks for Consistency...
c.......................................................................

      if (nstop.gt.99999) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)
         WRITE(*,*) 'If nstop greater than i5 format, then can'
         WRITE(*,*) 'obtain incorrect results from the code.'
         WRITE(*,*) 'Recode output of time step n for greater values.'
         WRITE(*,*)
CMPIINSERT_ENDIF_RANK
         stop
      endif
         

c     This is a zero-orbit-width version of cql3d. Thus, no 
c     finite-orbit-width FOW effects included.  We keep fow
c     related namelist names for backward compatibility.
      if (fow.ne.'disabled') then
         write(*,*)
         write(*,*) "This is a ZOW version of cql3d."
         write(*,*) "fow reset to 'disabled'"
         write(*,*)
      endif
      fow='disabled' ! just in case, reset anyway


      do k=1,ngen
          if(lossmode(k).eq.'simplban' .or. 
     +       lossmode(k).eq.'simplbn1' .or. 
     +       lossmode(k).eq.'simplbn2'     ) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'For lossmode(k) banana losses, with NB it is '// 
     +   'ADVISEABLE to consistently use ndeltarho = enabled or freya'
CMPIINSERT_ENDIF_RANK
          endif
          !YuP[01-2017] Adjusted/disabled simplbn1 and simplbn2 
          !These new options do not work properly yet.
          !One of problems: Failure in ngen=2 multiURF tests
          !because subr. deltar is setup/called for one general species only,
          !and so all deltar* arrays are saved for k=1 gen.species only.
          !As a temporary measure, reset to simplban:
          if(lossmode(k).eq.'simplbn1' .or. 
     +       lossmode(k).eq.'simplbn2' ) then
CMPIINSERT_IF_RANK_EQ_0
             WRITE(*,*)'lossmode(k)=simplbn1 or simplbn2 are not ready'
             WRITE(*,*)'--------- Resetting to simplban --------------'
CMPIINSERT_ENDIF_RANK
             lossmode(k)='simplban'
          endif ! YuP
      enddo ! k

c     Need to set tavg if f4d_out="tavg"
      if (f4d_out.eq."tavg" .and. tavg.eq."disabled") then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'f4d_out=tavg, but tavg=disabled'
CMPIINSERT_ENDIF_RANK
         stop 'Need to set tavg values'
      endif

c     Expand eegy, if lrzmax.gt.1, and values at lr_=2 are zero.
c     (This is for convenience of input).

      do ny=1,negyrg-1
         do ii=1,2
            do kk=1,ngen
               if (eegy(ny,ii,kk,2).eq.zero) then
                  do l=2,lrzmax
                     eegy(ny,ii,kk,l)=eegy(ny,ii,kk,1)
                  enddo
               endif
            enddo
         enddo
      enddo


      ntotal=ngen+nmax
      if (eqmod.ne."disabled") then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,399)
CMPIINSERT_ENDIF_RANK
 399     format(' WARNING: psimodel set =spline, for eqmod.ne.disabled')
         psimodel="spline"
      endif

      if (yreset.eq."enabled") then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,400) 
CMPIINSERT_ENDIF_RANK
 400     format(' WARNING: yreset set disabled, Needs recommissioning')
         yreset="disabled"
      endif
      if (transp.eq."enabled" .and. soln_method.ne."direct" .and.
     +    tfac.ge.0.) then 
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'WARNING:  Resetting tfac=-1. for transp=enabled/'
     +             //'soln_method not =direct'
CMPIINSERT_ENDIF_RANK
         tfac=-1.
      endif
      if (qsineut .eq. "enabled") locquas="disabled"
cBH100517:  Not needed.      if (lrzmax.eq.1) nrstrt=0
      if (lrzmax.eq.1 .and. meshy.ne."free") then
        meshy="free"
CMPIINSERT_IF_RANK_EQ_0
        WRITE(*,'(/"  WARNING: meshy has been changed to ""free"" ",
     +    "as lrzmax=1 (does not work otherwise)"/)')
CMPIINSERT_ENDIF_RANK
      endif
      if (soln_method.eq."it3drv".and. meshy.ne."fixed_y") then
        meshy="fixed_y"
CMPIINSERT_IF_RANK_EQ_0
        WRITE(*,'(/"  WARNING: meshy has been changed to ""fixed_y"" ",
     +    "as soln_method=it3drv (does not work otherwise)"/)')
CMPIINSERT_ENDIF_RANK
      endif
      if (zero.lt.rovera(1) .and. rovera(1).lt.1.e-8) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,398)
CMPIINSERT_ENDIF_RANK
 398     format(' WARNING: Resetting rovera(1) to minumum 1.e-8')
         rovera(1)=1.e-8
      endif
      if (ngauss .le. 0) then
        analegco="enabled"
      else
        analegco="disabled"
      endif

      if (eqmod.eq."enabled") then
CMPIINSERT_IF_RANK_EQ_0
        if (psimodel.eq."axitorus") WRITE(*,403)
CMPIINSERT_ENDIF_RANK
 403    format('WARNING: psimodel.eq.""axitorus"", see cqlinput_help')
        if (psimodel.ne."spline") then
CMPIINSERT_IF_RANK_EQ_0
          WRITE(*,401)
CMPIINSERT_ENDIF_RANK
 401      format(' WARNING: psimodel set =""spline"" '//
     +                      ' for eqmod=""enabled"" ')
          psimodel="spline"
        endif
      endif

      if (niong.eq.0 .and. sigmamod.eq."enabled") 
     +        stop 'set sigmamod=disabled'

      if (sigmamod.eq."enabled" .and. tandem.eq."enabled") then
CMPIINSERT_IF_RANK_EQ_0
          WRITE(*,402)
CMPIINSERT_ENDIF_RANK
 402      format('WARNING:sigmamod set disabled,Commission tandem case')
          sigmamod="disabled"
      endif
 
      if (totcrt(1).ne.zero .and. xjc(1).eq.zero)
     +        stop 'inconsistent totcrt(1),xjc(1)'
 
cBH170630: Changing lbdry="scale" to "consscal", following problem with
cBH170630: recent "scale" modification giving non-SS results, as pointed
cBH170630: out by Syun'ichi Shiraiwa (2017-06-23).
cBH170630: "consscal" conserves density at the v=0 boundary, and then
cBH170630: linearly rescales the distribution to obtain specified density.
         ! Note: it is not always a good idea to rescale the distr.func.
         ! Example: When there is a large-power NBI source, 
         ! comparing to the initial background 
         ! set in cqlinput. The value of xlndn00 
         ! (in ratio(k,lr_)=xlndn00(k,lr_)/runden, 
         ! which is the rescaling factor) is based 
         ! on the initial density, so it does not include particles from NBI.
         ! And the value of runden (or sden) does include all sources.
         ! In such a case, it is better to use lbdry(k)="conserv"
	 ! (Then the density of species k will increase.)
      do k=1,ngen
         if (lbdry(k).eq."scale")  lbdry(k)="consscal"
CMPIINSERT_IF_RANK_EQ_0
          WRITE(*,404)
CMPIINSERT_ENDIF_RANK
      enddo
 404  format('WARNING: reset obsolete (as of 170630) lbdry="scale" to'
     1        //' "consscal"')
      if (kelecg.ne.0) then
         if (redenc(1,kelecg).ne.zero .and. 
     +    lbdry(kelecg).ne."scale") then
CMPIINSERT_IF_RANK_EQ_0
            WRITE(*,203)
CMPIINSERT_ENDIF_RANK
         endif
 203     format('WARNING:lbdry(kelecg).ne.scale,'//
     +        ' redenc(1,kelecg).ne.0.?')
      endif

      if (zeffc(1).ne.zero .and. 
     +   (iprozeff.ne."parabola" .and. iprozeff.ne."prbola-t") )
     +         stop 'inconsistent zeffc(1) and iprozeff'

      if (gamaset.eq.zero .and. (ngen.gt.1.or.kelecg.ne.1)) 
     +         stop 'inconsistent gamaset'

      if (nso.gt.nsoa) 
     +         stop 'Check nso and nsoa'

      if (efiter.eq."enabled" .and. efswtch.eq."method1") then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,204)
CMPIINSERT_ENDIF_RANK
         efiter="disabled"
      endif
 204  format('WARNING: reset efiter=disabled, consistent with efswtch')

      if (noncntrl.ne.0 .and. (efswtch.ne."method2" .and.
     +     efswtch.ne."method3")) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,205)
CMPIINSERT_ENDIF_RANK
         noncntrl=0
      endif
 205  format('WARNING: resetting noncntrl=0, for efswtch as specified')

      if (negyrg.gt.negyrga) stop 'check negyrg'

      if (negyrga.gt.jx) stop 'check negyrga.le.jx, to use tam1,tam2'

      if (scatmod.eq."disabled".and.scatfrac.ne.1.) then
         scatfrac=1.
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,206)
CMPIINSERT_ENDIF_RANK
      endif
 206  format('WARNING: scatfrac set =1., for scatmod.eq.disabled')

      if (jfl.gt.jx) then
         jfl=jx ! if jfl>jx, reset
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,207)
CMPIINSERT_ENDIF_RANK
      endif
 207  format('WARNING: check jfl too big. Had to reset to jx')
      if (mod(jfl,2).eq.0) jfl=jfl-1  
      ! jfl needed to be odd because of jpxyh=(jfl+1)/2 in pltprppr.f

      if (pltra.ne."disabled" .and. knockon.ne."enabled") then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,208)
CMPIINSERT_ENDIF_RANK
      endif
 208  format('WARNING: Could have pltra problem, knockon.ne."enabled"')

      if (tandem.eq."enabled" .and. pltlim.ne.'x') then
         pltlim='x'
         pltlimm=1.0
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,209)
CMPIINSERT_ENDIF_RANK
      endif
 209  format('WARNING: Resetting pltlim=x, for tandem=enabled')

      if (trapmod.eq."enabled" .and. trapredc.gt.0.95) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,210)
CMPIINSERT_ENDIF_RANK
      endif
 210  format('WARNING:  Have had problems for trapredc.lt.0.99')


cBH131104:  ADD BUNCH more qualifications to use of ampfmod,
cBH131104:  ensuring not using conflicting code capabilites,
cBH131104:  for example, using with iterative electric field control.
      if (ampfmod.eq."enabled") then
         !do k=1,ngen
         if(soln_method.eq.'it3drv' .or. 
     +      soln_method.eq.'it3dv'      ) then
            WRITE(*,*)'Warning: Setting soln_method=direct:'
            WRITE(*,*)'      Code not presently set up for it3dv/it3drv'
            soln_method='direct'
         endif
         !enddo ! k=1,ngen
         if (lrz.ne.lrzmax) then
            WRITE(*,*)
            WRITE(*,*)'STOP: Not a good to use lrz.ne.lrzmax'
            WRITE(*,*)
            stop
         endif
cBH131104:  NEED to increase boundary options/include t-dep Vphib.
         if ((iproelec.eq.'parabola' .or. iproelec.eq.'prbola-t').or.
     +        (iproelec.eq.'spline' .or. iproelec.eq.'spline-t')) then
            if (efswtch.eq."method1".and.tmdmeth.eq."method1") then
               continue
            else
               WRITE(*,*)'STOP1: Incorrectly specd boundary elec fld'
               stop
            endif
         else
            WRITE(*,*)'STOP2: Incorrectly specd boundary elec fld'
            stop
         endif

         if (iproelec.eq.'spline' .or.iproelec.eq.'spline-t') then
            if (abs(ryain(1)).gt.em10.or.
     +          abs(ryain(njene)-1.d0).gt.em10) then
               WRITE(*,*)'STOP Incorrectly specd ryain for ampfar case'
               stop
            endif
         endif

         if (efflag.ne."toroidal") then
            write(*,*)'STOP: AmpFar assumes toroidal elec fld calc'
            stop
         endif
     
      endif  !On ampfmod
      
      if (ampfmod.eq."enabled".and.cqlpmod.eq."enabled") then
         write(*,*)
         WRITE(*,*)'STOP: ampfmod.ne.enabled with cqlpmod.eq.enabled'
         write(*,*)
         stop
      endif

      if (eseswtch.eq."enabled") then
         if (.not.(cqlpmod.eq."enabled".and. sbdry.eq."periodic")) then
            stop 'Check use of eseswtch'
         elseif (efswtch.ne."disabled") then
CMPIINSERT_IF_RANK_EQ_0
            WRITE(*,211) 
CMPIINSERT_ENDIF_RANK
 211        format('WARNING:  Re-setting efswtch=disabled, for now')
            efswtch="disabled"
         endif
      endif

      if (bootcalc.ne."disabled".and.lrzdiff.ne."disabled") then
         stop 'bootcalc.ne.disabled .and. lrzdiff.ne.disabled'
      endif

      if (efswtch.eq."method5") then
         if (eqmod.ne."enabled") stop 'Efswtch=method5,eqmod.ne.enabled'
         if (eqsource.ne."eqdsk") then
CMPIINSERT_IF_RANK_EQ_0
            WRITE(*,216)
CMPIINSERT_ENDIF_RANK
 216        format('WARNING: efswtch=method5 with eqsource.ne.eqdsk',/,
     +           'Not checked out')
         endif
         if (efflag.ne."parallel") then
CMPIINSERT_IF_RANK_EQ_0
            WRITE(*,219)
CMPIINSERT_ENDIF_RANK
 219        format('WARNING: efswtch=method5 with efflag.ne.parallel',/,
     +           'Could cause problems with RFP.  See cqlinput_help.')
         endif
      endif

      if (.not.(eqmod.eq."enabled" .and. eqsource.eq."eqdsk"))  then
         if (eqsym.eq."none") then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)
         WRITE(*,*)'WARNING:  Re-setting eqsym="average":'  
         WRITE(*,*)'          It should not ="none" for eqmod/eqsource'
         WRITE(*,*)
CMPIINSERT_ENDIF_RANK
         endif
      endif

      if (nv.gt.nva) then
         stop 'SXR: nv .gt. nva.  Also check length of XR input arrays.'
      endif

      if (transp.eq."enabled" .and. relaxden.ne.1.) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,217)
CMPIINSERT_ENDIF_RANK
 217     format('WARNING: transp.eq.enabled .and. relaxden.ne.1.',/,
     +      '     relaxden defn changed from divisor to multiplier,',/,
     +      '     BH: April 23, 2001')
      endif

c     Checking rya input (rzset.ne."disabled"):
      if (rzset.ne."disabled") then
      if (rya(0).ne.zero) then
         stop 'rya(0).ne.0. ==>  error in namelist input: rya(1)=?'
      endif
      if (lrz.gt.1) then
         do ll=2,lrz
            if ( (rya(ll)-rya(ll-1)).le.0. ) stop 'Check rya() input'
         enddo
      endif
      endif

CMPIINSERT_IF_RANK_EQ_0
      if (transp.eq."enabled" .and. 
     +   (pinch.ne."simple"  .or. pinch.ne."simplen" .or. 
     +   pinch.ne."case1"  .or. pinch.ne."case1n" .or. 
     +   pinch.ne."case2"  .or. pinch.ne."case2n" .or. 
     +   pinch.ne."case3"  .or. pinch.ne."case3n" )) WRITE(*,218)
CMPIINSERT_ENDIF_RANK
 218  format('WARNING: transp=enabled, Check settings of pinch')

      if (transp.eq."enabled" .and. radcoord.ne."sqtorflx") then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'WARNING: transp=enabled not setup for '//
     +                   'radcoord.ne.sqtorflx'
CMPIINSERT_ENDIF_RANK
c         stop
      endif
      

      if (soln_method.eq.'itsol'.or. soln_method.eq.'itsol1' .or.
     +    soln_method.eq.'it3dv'.or. soln_method.eq.'it3drv') then
         if (cqlpmod.eq.'enabled' .or.
     +       symtrap.ne.'enabled' ) then
CMPIINSERT_IF_RANK_EQ_0
                WRITE(*,*)'ainsetva: Check/re-configure itsol storage'//
     +                    ' for cqlpmod/symtrap/nadjoint.ne.0 cases'
CMPIINSERT_ENDIF_RANK
                STOP
         endif
      endif

c     Necessary for some code logic:
      if (transp.eq.'disabled' .and. soln_method.eq.'it3drv') then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'ainsetva: If transp.eq.disabled, then soln_method'
         WRITE(*,*)'ainsetva: should eq direct,itsol,itsol1, or it3dv'
         WRITE(*,*)'ainsetva: Resetting soln_method it3drv ==> it3dv'
CMPIINSERT_ENDIF_RANK
         soln_method='it3dv'
      endif

c     Check on consistency of urfmod
      ikrf=0
      if (urfmod.ne."disabled" .and. lh.eq."disabled" .and.
     +    ech.eq."disabled" .and. fw.eq."disabled") then
         ikrf=0
      else
         ikrf=1
      endif
      if (urfmod.ne."disabled") then
         do krf=1,nmodsa
cBH100829            if (rftype(krf).eq."notset") ikrf=ikrf+1
            if (rftype(krf).ne."notset") ikrf=ikrf+1
         enddo
      endif
      if (ikrf.eq.0) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)
         WRITE(*,*)'Warning: setting urfmod=disabled, since'
         WRITE(*,*)'         lh, ech, and fw, =disabled, '
         WRITE(*,*)' and all rftype(1:nmodsa)= notset'
         WRITE(*,*)
CMPIINSERT_ENDIF_RANK
         urfmod="disabled"
      endif
      do krf=1,nmodsa
         if (nrfspecies(krf).gt.ngen) then
CMPIINSERT_IF_RANK_EQ_0
            WRITE(*,*) 'Check nrfspecies.lt.ngen'
CMPIINSERT_ENDIF_RANK
            stop
         endif
      enddo
         

c     Checking consistency of netcdfvecXX settings:
      if (netcdfvecc.ne."disabled") then
         if (netcdfvece.ne."disabled" .and. netcdfvece.ne.netcdfvecc)
     +                     stop 'STOP:check netcdfvecXX for consistency'
         if (netcdfvecrf.ne."disabled" .and. netcdfvecrf.ne.netcdfvecc)
     +                     stop 'STOP:check netcdfvecXX for consistency'
         if (netcdfvecal.ne."disabled" .and. netcdfvecal.ne.netcdfvecc)
     +                     stop 'STOP:check netcdfvecXX for consistency'
      endif
      if (netcdfvece.ne."disabled") then
         if (netcdfvecrf.ne."disabled" .and. netcdfvecrf.ne.netcdfvece)
     +                     stop 'STOP:check netcdfvecXX for consistency'
         if (netcdfvecal.ne."disabled" .and. netcdfvecal.ne.netcdfvece)
     +                     stop 'STOP:check netcdfvecXX for consistency'
      endif

      if (netcdfvecrf.ne."disabled") then
         if (netcdfvecal.ne."disabled" .and. netcdfvecal.ne.netcdfvecrf)
     +                     stop 'STOP:check netcdfvecXX for consistency'
      endif

         
      if (frmodp.ne."disabled") then
         if (nso.ne.1) stop 'STOP: frmod.ne.disabled, NEED nso=1'
         if (kfrsou.eq.0) 
     +           stop 'STOP: frmod.ne.disabled, NEED kfrsou.ne.0'
      endif

c     In view of checking settings in sub frset, this condition
c     should not (cannot) occur.
      if (beamplsep.ne."disabled") then
         if (frmodp.eq."disabled")
     +        stop 'STOP:frmod.eq.disabled but beamplse.ne.disabled'
      endif
         
c     Checking consistency of rdcmod and lrzdiff
      if (rdcmod.ne."disabled" .and. lrzdiff.ne."disabled") then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)'ainsetva:  STOP, rdcmod requires lrzdiff=disabled'
CMPIINSERT_ENDIF_RANK
         stop
      endif

c     Soft X-ray
      if (softxry.ne."disabled" .and. kelecg.eq.0) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)
         WRITE(*,*)'**** ainsetva: Inappropriate calc of SXR ******'
         WRITE(*,*)
CMPIINSERT_ENDIF_RANK
      endif

c     Check rd(),x_sxr(),z_sxr() for previous input method, in
c     which they were scalars, i.e., only one initial position was
c     specified.  Fill in 2:nv, if necessary.
      if (nv.gt.1) then

         if (x_sxr(1).eq.zero .and. z_sxr(1).eq.zero) then  !rd input
            if (rd(2).eq.zero) then
               do i=2,nv
                  rd(i)=rd(1)
                  thetd(i)=thetd(1)
               enddo
            endif
         else  !x_sxr,z_sxr input
            if (x_sxr(2).eq.zero .and. z_sxr(2).eq.zero) then
               do i=2,nv
                  x_sxr(i)=x_sxr(1)
                  z_sxr(i)=z_sxr(1)
               enddo
            endif
         endif

      endif  ! on nv



c     NPA
      if (npa_diag.ne."disabled" .and. niong.eq.0) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)
         WRITE(*,*)'**** ainsetva: Inappropriate calc of NPA ******'
         WRITE(*,*)
CMPIINSERT_ENDIF_RANK
      endif
      if (npa_diag.ne."disabled" .and. softxry.ne."disabled") then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)
         WRITE(*,*)'ainsetva: STOP'
         WRITE(*,*)'ainsetva: Code not complete set up for simultaneous'
         WRITE(*,*)'ainsetva: calc of SXR and NPA.  Need separate'
         WRITE(*,*)'ainsetva: storage for mx related variables, unless'
         WRITE(*,*)'ainsetva: the related data is recalculed for each'
         WRITE(*,*)'ainsetva: use of the diagnostics.'
CMPIINSERT_ENDIF_RANK
         STOP
      endif

c     Check rd_npa(),x_npa(),z_npa() for previous input method, in
c     which they were scalars, i.e., only one initial position was
c     specified.  Fill in 2:nv_npa, if necessary.
      if (nv_npa.gt.1) then

         if (x_npa(1).eq.zero .and. z_npa(1).eq.zero) then  !rd input
            if (rd_npa(2).eq.zero) then
               do i=2,nv_npa
                  rd_npa(i)=rd_npa(1)
                  thetd_npa(i)=thetd_npa(1)
               enddo
            endif
         else  !x_npa,z_npa input
            if (x_npa(2).eq.zero .and. z_npa(2).eq.zero) then
               do i=2,nv
                  x_npa(i)=x_npa(1)
                  z_npa(i)=z_npa(1)
               enddo
            endif
         endif

      endif  ! on nv_npa



c.......................................................................
c     Check on structure of integer*1, if (urfmod.ne."disabled").
c     The following has been set up for the Absoft fortran
c     compiler.  Other compilers may be different.
c     
c.......................................................................

      if (urfmod.ne."disabled") then
         ip0=-128
         ip255=127
         iu0=ip0+128
         iu255=ip255+128
         if (iu0.ne.0 .or. iu255.ne.255) then
CMPIINSERT_IF_RANK_EQ_0
            WRITE(*,*) 'ip0,ip255,iu0,iu255 ',ip0,ip255,iu0,iu255
CMPIINSERT_ENDIF_RANK
            stop 'check pack unpack in urfpkun.f'
         endif
      endif
	

c.......................................................................
cl    4.3 Warnings about Changed Use of Variables
c         (For example, previously used as numeric and character
c          input.
c.......................................................................

CMPIINSERT_IF_RANK_EQ_0
      do k=1,ngen
         if (fpld(1,k).ne.zero) WRITE(*,212)
      enddo
 212  format(//,'WARNING: Have changed namelist input options',/,
     +           8x,    ' fpld_dsk ==> -1.0',/,
     +           8x,    ' fpld_ds1 ==> -2.0',//)

      if (izeff.ne."backgrnd") WRITE(*,213)
 213  format(//,'WARNING:  Check use of izeff. See cqlinput_help',//)
CMPIINSERT_ENDIF_RANK

      if (nconteq.ne."psigrid".and.nconteqn.eq.0) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,214)
CMPIINSERT_ENDIF_RANK
 214     format(//,'Must set nconteqn, if nconteq.ne.psigrid',//)
         stop 'Must set nconteqn; See cqlinput_help'
      endif

      if (nconteq.eq."psigrid".and.nconteqn.ne.0) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,215)
CMPIINSERT_ENDIF_RANK
 215     format(//,'Must reset nconteq, if nconteqn.ne.0',//)
         stop 'Must reset nconteq: See cqlinput_help'
      endif

      if (nconteq.eq."psigrid".and.nconteqn.ne.0) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,221)
CMPIINSERT_ENDIF_RANK
 221     format(//,'NO LONGER RECOMMENDED: See cqlinput_help (BH171211)'
     1        ,//)
      endif


c.......................................................................
cl    5. Plot flags
c.......................................................................

c     could be problems if nrskip.eq.0 and the irzplt() are not
c     all distinct (or zero), or are .gt. lrz.
      if (nrskip.eq.0) then
         do i=1,lrors
            if (irzplt(i).gt.lrors) then
CMPIINSERT_IF_RANK_EQ_0
               WRITE(*,*)'WARNING: irzplt(i).gt.lrors are set =0, i=',i
CMPIINSERT_ENDIF_RANK
               irzplt(i)=0
            endif
         enddo
         do i=1,lrors-1
            do j=i+1,lrors
               if (irzplt(j).eq.irzplt(i) .and.
     +             irzplt(j).ne.0 ) then
CMPIINSERT_IF_RANK_EQ_0
                  WRITE(*,*)
     +                   'WARNING: Non-distinct irzplt(j) set =0, j=',j
CMPIINSERT_ENDIF_RANK
                  irzplt(j)=0
               endif
            enddo
         enddo
      endif

      if (eqmod.ne."enabled") pltvs="rho"
      if (pltvs.ne."rho".and. pltvs.ne."psi") pltvs="rho"
      nirzplt=0  ! keeps count of number of radial plot positions

c Having an mplot (compiler?) problem on viz machine:
c      write(*,*)'ainsetva before nrskip setting mplot: ll=1:lrors,mplot'
c      do ll=1,lrors
c         write(*,*) ll,mplot(ll)
c      enddo

      do 500 ll=1,lrors
        if (nrskip .ne. 0) then
          if ((ll-1) .eq. (ll-1)/nrskip*nrskip) then
             mplot(ll)="enabled"
             nirzplt=nirzplt+1
          endif
        else
          do 510 lu=1,lrors
            if (ll .eq. irzplt(lu)) then
               mplot(ll)="enabled"
               nirzplt=nirzplt+1
            endif
 510      continue
        endif
 500  continue
      if (lrors .eq. 1) then
         mplot(1)="enabled"
         nirzplt=1
      endif


      do k=1,ngen
        if (   lbdry(k).ne. "conserv".and. lbdry(k).ne."consscal"
     +   .and. lbdry(k).ne. "fixed"  .and. lbdry(k).ne."scale" 
     +   ) then
CMPIINSERT_IF_RANK_EQ_0
        WRITE(*,*)"lbdry(k)=",lbdry(k)
        WRITE(*,*)"lbdry(k) can only be 'conserv', 'consscal'," // 
     +  "'fixed', or 'scale'  "
CMPIINSERT_ENDIF_RANK
           STOP
        endif
      enddo

c Part of printed output re viz problem
c      write(*,*)'ainsetva after nrskip: nirzplt= ',nirzplt
c      write(*,*)'ainsetva after nrskip:   ll=1:lrors,mplot '
c      do ll=1,lrors
c         write(*,*) ll,mplot(ll)
c      enddo

      return
      end
