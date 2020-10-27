c
c
      subroutine frinitz(
     +        nprim,nimp,nion,ibion,namep,namei,atw,fd,smth)
      implicit integer (i-n), real*8 (a-h,o-z)
c
      include 'param.h'
      include 'comm.h'


c..................................................................
c     This routine communicates between common blocks of CQL3D and the
c     segregated CRAY32 (NFREYA ONETWO version) common blocks,
c     by passing subroutine arguments.
c
c     smooth in comm.h and frcomm.h is passed from frcomm to comm,
c       through argument smth.
c     fd specified in frinitl passed here to give dt mixture ratio.
c     nprim,nimp,nion,ibion,namep,namei,atw are in frcomm.h, not comm.h.
c     These variables are determined from comm.h data, passed to 
c     frcomm.h through the arguments of frinitz.
c       
c..................................................................

      character*8 namep,namei
      dimension namep(*),namei(*),atw(*)

      smooth_=smth

c.....................................................................
c     Determine the index of the general species electrons, the
c     background Maxwellian species electrons, and the first ionic
c     species. There must be at least one electron species and at least
c     one ionic species because the code forces quasineutrality.
c   
c     NOTE: Check for overlap with ainspec.f  (BobH, 990704).
c.....................................................................

      kelecg=0
CDIR$ NEXTSCALAR
      do 1000 k=1,ngen
        if (fmass(k) .lt. 1.e-27) then
          kelecg=k
          goto 1001
        endif
 1000 continue
 1001 continue
CDIR$ NEXTSCALAR
      do 1002 k=ngen+1,ntotal
        if (fmass(k) .lt. 1.e-27) then
          kelecm=k
          goto 1003
        endif
 1002 continue
 1003 continue
      if (kelecg.gt.0) then
        kelec=kelecg
      else
        kelec=kelecm
      endif
cBH180908      if (kelec.eq.0) call diagwrng(9)
      if (kelec.eq.0) then
CMPIINSERT_IF_RANK_EQ_0
         WRITE(*,*)
         WRITE(*,*) 'WARNING: Unphysical plasma, only one species.'
         WRITE(*,*) '         Have not checked with NBI turned on.'
         WRITE(*,*)
CMPIINSERT_ENDIF_RANK
      endif


CDIR$ NEXTSCALAR
c     Find first ion species (general or maxwellian).
      do 1005 k=1,ntotal
        if (fmass(k).gt.1.e-27) then
          kionn=k
          goto 1006
        endif
 1005 continue
 1006 continue

c..................................................................
c     Define nprim, nimp and nion as in ONETWO. This involves getting
c     rid of possible electron species which may reside in CQL3D
c     Recall nimp is the number of impurity species input in
c     namelist frstup. ntotal is the total number of species.
c     nprim is the number of primary species, and this number is
c     also input in frstup.
c..................................................................
      nion=nprim+nimp
c..................................................................
c     Now define the namep and namei array as used in ONETWO
c..................................................................

      write(*,*) 'frinitz: nprim,nimp = ',nprim,nimp

      kk=0
      do 500 k=1,ntotal
        ! exclude e for the stopping cross-section:
        if (k.eq.kelecg .or. k.eq.kelecm) go to 500 
        kk=kk+1
        if (kk.gt.nprim) go to 500
        if (k.eq.kfrsou) ibion=kk
        namep(kk)=kspeci(1,k)
 500  continue
      if(ibion.le.0)then
        WRITE(*,*)'frinitz: ibion=0, probably kfrsou=0; Set to proper k'
        STOP
      endif
cBH081022      kk=0
cBH081022      do 510 k=ntotal+1-nimp,ntotal
cBH081022        kk=kk+1
cBH081022        namei(kk)=kspeci(1,k)
cBH081022 510  continue

cBH081022:  Find impurity species name (not previously working)
cBH081022:  Presently set up for only one species. (Else, more coding).
cBH081022:  Take impurities to be distinct species with atw.gt.4.

      if (nimp.gt.1) then
         write(*,*)
         WRITE(*,*)'STOP:  NBI (Freya setup) only for nimp.le.1.'
         STOP
      endif
      do kk=1,nimp
         do k=1,ntotal
            if (fmass(k)/proton .gt. 4.1) then
c BH081022:  Add code here, if need nimp.gt.1.
               namei(kk)=kspeci(1,k)
               ksave=k  !for use below, for atw of impurity
            endif
         enddo
      enddo
               
               

      write(*,*) 'frinitz:namep(1),namep(2),namei(1),ibion   ',
     +            namep(1),namep(2),namei(1),ibion


c----------------------------------------------------------------
c     determine atomic weights of primary ions
c----------------------------------------------------------------
      do 3410 i=1,nprim
        atw(i) = 0.
        if(trim(namep(i)).eq.'h'  .or. 
     +     trim(namep(i)).eq.'H') atw(i)=1.
     
        if(trim(namep(i)).eq.'d'  .or. 
     +     trim(namep(i)).eq.'D') atw(i)=2.
     
        if(trim(namep(i)).eq.'t'  .or. 
     +     trim(namep(i)).eq.'T') atw(i)=3.
     
        if(trim(namep(i)).eq.'dt' .or. 
     +     trim(namep(i)).eq.'DT' .or. 
     +     trim(namep(i)).eq.'Dt' .or. 
     +     trim(namep(i)).eq.'dT') atw(i) = fd*2. + (1.-fd)*3.
     
        if(trim(namep(i)).eq.'he' .or. 
     +     trim(namep(i)).eq.'HE' .or. 
     +     trim(namep(i)).eq.'He') atw(i)=4.
     
        write(*,*) 'frinitz: i, trim(namep(i)), atw(i)',
     +                       i, trim(namep(i)), atw(i)
        if(atw(i).eq.zero) call frwrong(1)
 3410 continue


c----------------------------------------------------------------
c     determine atomic weights of impurity ions
c----------------------------------------------------------------
      if(nimp.eq.0) go to 3430
      do 3420 i=1,nimp
        k = nprim + i
        atw(k) = 0.
        if(trim(namei(i)).eq.'he' .or. 
     +     trim(namei(i)).eq.'HE' .or. 
     +     trim(namei(i)).eq.'He') atw(k)= 4.
     
        if(trim(namei(i)).eq.'b' .or. 
     +     trim(namei(i)).eq.'B' ) atw(k)= 11.  ! YuP added [2015]

        if(trim(namei(i)).eq.'c' .or. 
     +     trim(namei(i)).eq.'C' ) atw(k)= 12.
     
        if(trim(namei(i)).eq.'o' .or. 
     +     trim(namei(i)).eq.'O' ) atw(k)= 16.
        
        if(trim(namei(i)).eq.'si' .or. 
     +     trim(namei(i)).eq.'SI' .or. 
     +     trim(namei(i)).eq.'Si') atw(k)= 28.
     
        if(trim(namei(i)).eq.'ar' .or. 
     +     trim(namei(i)).eq.'AR' .or. 
     +     trim(namei(i)).eq.'Ar') atw(k)= 40.
     
        if(trim(namei(i)).eq.'cr' .or. 
     +     trim(namei(i)).eq.'CR' .or. 
     +     trim(namei(i)).eq.'Cr') atw(k)= 52.
     
        if(trim(namei(i)).eq.'fe' .or. 
     +     trim(namei(i)).eq.'FE' .or. 
     +     trim(namei(i)).eq.'Fe') atw(k)= 56.
        
        if(trim(namei(i)).eq.'ni' .or. 
     +     trim(namei(i)).eq.'NI' .or. 
     +     trim(namei(i)).eq.'Ni') atw(k)= 59.
        
        if(trim(namei(i)).eq.'kr' .or. 
     +     trim(namei(i)).eq.'KR' .or. 
     +     trim(namei(i)).eq.'Kr') atw(k)= 84.
        
        if(trim(namei(i)).eq.'mo' .or. 
     +     trim(namei(i)).eq.'MO' .or. 
     +     trim(namei(i)).eq.'Mo') atw(k)= 96.
        
        if(trim(namei(i)).eq.'w' .or. 
     +     trim(namei(i)).eq.'W' ) atw(k)= 184.
        
        if(trim(namei(i)).eq.'mixt') atw(k)= 20.
        
        if(trim(namei(i)).eq.'a' .or. 
     +     trim(namei(i)).eq.'A' ) atw(k)=int(fmass(ksave)/proton +0.1)
                                            !Compatible with cql3d added
                                            !impurity, for given zeff.
        if(atw(k).eq.zero) call frwrong(2)
 3420 continue
 3430 continue

      if(nimp.gt.0) then
         do i=1,nimp
            write(*,*)'frinitz: trim(namei(i)),atw(nprim+i)',
     +                          trim(namei(i)),atw(nprim+i)
         enddo
      endif

      return
      end
