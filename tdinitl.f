c
c
      subroutine tdinitl ! called only at n=0
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     This routine initializes the arrays required by CQL3D.
c     It also controls the initialization of distribution
c     functions and quasilinear coefficients.
c..................................................................

      include 'param.h'
      include 'comm.h'

      include 'name.h'
CMPIINSERT_INCLUDE     

      character*8 icall,iplotsxr
      character*128 filenm ! template for file name with data 
      real*8 a_new(njene)  ! working array

      if (n.gt.0) return

c..................................................................
c     set defaults for namelisted quantities and other variables
c..................................................................

      call aindflt
      call eqindflt
      call urfindfl
      call aindflt1

c.......................................................................
c     Read namelist from previous run if rerun case. Thus namelist from
c     this run should only reset some diagnostic or numerical variables.
c     This call will position pointer in file distrfunc to next read f.
c.......................................................................

cBH110201      if (nlrestrt.ne."disabled") call tdreadf(1)
      if (nlrestrt.eq."enabled") call tdreadf(1)

c..................................................................
c     read in namelist input for CQL3D
c..................................................................
      open(unit=2,file='cqlinput',status='old')
        read(2,setup)
        read(2,trsetup)
        read(2,sousetup)
        read(2,eqsetup)
        read(2,rfsetup)
      close(2)

      if (partner.eq."selene") then
        open(unit=18,file='kcqlsel',status='old')
        read(18,2) ncount,noplots
        eqdskalt="enabled"
        close(unit=18)
      endif
 2    format(i5,a8)

c..................................................................
c     Call routine which finds electron and ion species indices.
c..................................................................

      call ainspec

c.......................................................................
c     set variables dependent on input variables
c     (May add an impurity ion, increasing species indices, if
c      iprozeff.ne.'disabled').
c.......................................................................

      call ainsetva

c.......................................................................
c     Allocate arrays and radial transport arrays, if required
c.......................................................................

      call ainalloc
      if (ampfmod.eq."enabled" .and.cqlpmod.ne."enabled") call ampfalloc
      
      if (transp.eq."enabled" .and. cqlpmod.ne."enabled") call tdtraloc
c%OS  if (transp.eq."enabled") call tdtraloc
      if (transp.eq."enabled" .and. cqlpmod.eq."enabled") call wpalloc

c.......................................................................
c     print namelists
c.......................................................................

      if (nmlstout.eq."enabled") then
         open(unit=2,file='cqlinput',delim='apostrophe',status='old')
         write(6,*)'  In tdinitl: '
         write(6,setup0)
         write(6,setup)
         write(6,trsetup)
         write(6,sousetup)
         write(6,eqsetup)
         write(6,rfsetup)
         close(2)
      elseif (nmlstout.eq."trnscrib") then
         write(6,*)'  In tdinitl: '
         call ain_transcribe("cqlinput")
      else
         write(6,*)
         write(6,*) 'mnemonic = ',mnemonic
         write(6,*)
      endif

c..................................................................
c     Call the initialization routines for the appended modules..
c..................................................................
      call eqinitl
      call urfinitl ! mrfn= is set here also
      call frinitl
      open(unit=2,file='cqlinput',delim='apostrophe',status='old')
      call frset(lrz,noplots,nmlstout)   ! Uses unit 2
      close(2)

      if (machine .ne. "toroidal") call tdwrng(1)
      if (lrzdiff.eq."enabled" .and. frmodp.eq."enabled") 
     +  call diagwrng(18)

c.....................................................................
c     This routine initializes the normalized theta mesh that is used
c     if the meshes on the various flux surfaces are linked.
c.....................................................................

      call tdtrmuy

c.....................................................................
c     Call routine which initializes certain arrays using input
c     data. Also initialize the radial (rz) mesh and some plasma
c     parameters arrays.
c.....................................................................

      if(lrzmax.gt.1) call tdxinitl

c.....................................................................
c     Call routines to initialize any time-dependent (parabolic)   
c     profiles, and to set them up in place of profiles 
c     otherwise set up in tdxinitl.
c.....................................................................

      if (nbctime.gt.0) then
         call profiles
      endif

c.....................................................................
c     Determine mesh normalization constant vnorm.
c.....................................................................

      call ainvnorm

c.......................................................................
cl    1.1 Loop over all radial mesh points lr=1,..,lrzmax to determine
c     the plasma and equilibrium parameters on the full radial mesh
c     Also determine the mesh parallel to the magnetic field line
c.......................................................................

      if (cqlpmod.eq."enabled" .and. numclas.eq.1 .and. ls.eq.lsmax)then
        lz=lz/2+1
        lsmax=lsmax/2+1
        ls=ls/2+1
      endif

      do 110 ll=lrzmax,1,-1

c.......................................................................
c     Sets up the parameters on each lrzmax flux surface. Thus does not
c     define lr_=lrindx(l_), but lr_=l_, l_=1,lrzmax
c     (Not like in tdnflxs)
c......................................................................

        l_=ll
        lr_=ll
c     these two indices should not be used in loop 110 => set to -1, to detect
c     out of bound errors
        lmdpln_=-1
        indxlr_=-1

c..................................................................
c     Call an initialization routine which determines flux surface
c     geometry and magnetic field structure.
c..................................................................

        call aingeom !-> eqcoord-> eqfndpsi-> eqorbit-> trace flux surf.
             ! solr(l,lr_), solz(l,lr_) are R,Z coords. of flux surface
c.......................................................................
c     Initialize mesh along the magnetic field line, as well as
c     density and temperature profiles if cqlpmod=enabled
c.......................................................................

        call micxiniz

c..................................................................
c     Copy some radial and non-time dependent diagnostic quantities
c..................................................................

        call tdtoarad

 110  continue ! ll=lrzmax,1,-1

c      do ir=1,lrz
c         WRITE(*,'(a,i3,4e13.5)')'ir,rya,rpcon,rmcon,equilpsi=',
c     +                ir,rya(ir),rpcon(ir),rmcon(ir),equilpsi(ir)
c      enddo
             
 

c.......................................................................
c     Determine equilibrium parameters at lower half of cross-section
c     on lrindx(1) flux surface
c.......................................................................

      if (cqlpmod.eq."enabled" .and. numclas.eq.1 .and. ls.eq.lsmax)then
        lz=2*(lz-1)
        lsmax=2*(lsmax-1)
        ls=2*(ls-1)
        call wploweq
      endif

c.....................................................................
c     Redefines mu-mesh at midplane if needed
c.....................................................................

      if (cqlpmod.eq."enabled" .and. meshy.eq."fixed_mu" .and.
     .  tfac.lt.0.0) call wptrmuy

c.......................................................................
cl    1.1.2 Initialize some plasma parameters on entire lrzmax mesh
c.......................................................................

      call ainpla

c.....................................................................
cl    1.1.3 First Loop over spatial variable index for which the equations
c     will be solved to compute y-mesh on all l_ indices
c     (this takes micxinit out of ainitial subroutine)
c.....................................................................

      do 113 ll=lrors,1,-1

c......................................................................
c     determine local variables depending on variable index l_
c......................................................................

        call tdnflxs(ll)

c.......................................................................
c     call a routine to determine meshes y, x and related quantities
c.......................................................................

        call micxinit

 113  continue

      ieq_tot=0
      do ll=1,lrors
         ieq_tot=ieq_tot+inewjx_(ll) ! inewjx_() is defined in micxinit
         if (ll.eq.1) then
            ieq_(ll)=1
         else  ! Eqn no. at beginning of each flux surface:
            ieq_(ll)=ieq_(ll-1)+inewjx_(ll-1)  
         endif
      enddo
      ieq_(lrors+1)=ieq_tot



c.......................................................................
c     adapt dyp5(iyh) and dy(iyh) such that sum(cynt2) over i is constant
c     on all l_ surfaces
c.......................................................................

      if (nchgdy .eq. 1) then
        do 114 ll=1,lrors
          call tdnflxs(ll)
          call wpchgdy
 114    continue
      endif

      ! YuP[Sept-2014] moved this part outside of ainitial: 
      ! needed for proper work of FOW-logic.
      do ll=lrors,1,-1
        call tdnflxs(ll)
        call micxinil ! integration coefficients, mostly for cfpcoefn:
        ! gamma,alphan,cog,cosz,zmaxpsi (uses bbpsi,solrz from micxiniz)
        if (l_ .eq. lmdpln_) then
           ! calculate tau()==v*tauB_ZOW [also dtau(), vptb()]
           if(taunew.eq."enabled") then
              call baviorbt ! must be called after micxiniz
           else
              call baviorbto
           endif
        endif
      enddo

c....................................................................
c     Set some more radial mesh quantities.
      call tdrmshst !must be called AFTER: aingeom,tdtoarad,micxinit,micxinil
c....................................................................

cBH160530
c     Must be called after micxiniz and baviorbt/baviorbto:
cYuP      if (lossmode(1).eq.'simplban' .or. lossmode(1).eq.'simplbn1'
cYuP      +     .or.ndeltarho.ne.'disabled') 
cYuP      +     call deltar
cYuP[07-2017] exclude simplban  - deltar is not needed.
      if (lossmode(1).eq.'simplbn1'
     +     .or.ndeltarho.ne.'disabled') 
     +     call deltar

c.....................................................................
cl    1.2 Loop over spatial variable index for which the equations will be
c     solved. Compute initial distributions and sources etc.
c.....................................................................

      do 120 ll=lrors,1,-1

c......................................................................
c     determine local variables depending on variable index l_
c......................................................................

        call tdnflxs(ll)

c.....................................................................
c     transfer input data from CQL3d to CQL
c.....................................................................

        call tdstin

c.....................................................................
c     Call spatial variable "surface" initialization routine
c..................................................................

        call ainitial ! mrfn= is set here also (in vlf or vlh)

c..................................................................
c     Call a routine to initialize some urf module routines.
c..................................................................

c     Tried following, but gives problem of code stopping?? [BH041111].
c        if (cqlpmod.ne."enabled".and.urfmod.ne."disabled") call urfsetup
        if (cqlpmod.ne."enabled") call urfsetup ! mrfn= is set here also

c..................................................................
c     Copy some diagnostic quantities
c..................................................................

        call tdtoaray

 120  continue
c----------------------------------------------------------------------

c..................................................................
c     If running with the Brambilla code, create and dump and eqdsk
c     file for Brambilla ray tracing code to read in.
c..................................................................

      if (partner.eq."bramb") then
        call tdeqdsk
      endif

c....................................................................
c     Set some more radial mesh quantities.
c....................................................................

      call tdrmshst

c..................................................................
c     Call first-order radial orbit-shift module
c..................................................................

cBH160530      if (ndeltarho.ne."disabled") call deltar
         
c....................................................................
c     Call routines to set up radial integration coefficients, radial
c     diffusion and advection coefficients, certain flags and the radial
c     Chang-Cooper weights.
c....................................................................

      if (transp.eq."enabled") then
        if (cqlpmod .ne. "enabled") then
c     warning: these routines use lrors, not lrzmax
          call tdtrvint

c         Obtain radial diffusion d_rr, and optionally
c         pinch velocity d_r, from file.  In this case, the
c         d_rr and d_r can be multiplied by a time-dependent
c         scale factor at point of use, in trans.h and tdtravct.
c         The read in d_rr/d_r won't be changed.
          if (difus_io(1).eq."drrin") then
             call diffus_io(3) !Reads d_rr
             if (n_d_rr.ne.0) then
                !Getting time-dep scale factor at t=0
                do k=1,n_d_rr   !n_d_rr obtained in diffus_io
                   drrt(k)=difus_io_scale(k,1)
                enddo
             endif
          elseif(difus_io(1).eq."drrdrin") then
             pinch="disabled"   !Since getting d_r from file
             call diffus_io(4) !Reads d_rr and d_r
             if (n_d_rr.ne.0) then
                !Getting time-dep scale factor at t=0
                do k=1,n_d_rr   !n_d_rr obtained in diffus_io
                   drrt(k)=difus_io_scale(k,1)
                   drt(k)=difus_io_scale(k,2)
                enddo
             endif
          else
             call tdtrdfus
          endif

          call tdtrflg
          call tdtrwtl
        else
          call wpinitl
        endif
      endif
      call tddiag

c%OS  
c     check some mesh values
      call wpmshchk
c%OS  


c.....................................................................
c     Set up sigma-v calculation
c.....................................................................

      if (sigmamod.eq."enabled")  then
        icall="first"
        call sigv(icall)
      endif

c.....................................................................
c     Set up soft X-ray and NPA diagnostic.
c.....................................................................

      if (lrzmax .ge. 3) then
c**bh if (softxry.eq."enabled".and.kelecg.gt.0.and.eqmod.ne."enabled")
        if (softxry.ne."disabled".and.kelecg.gt.0)
     1    then
          icall="first"
          if (noplots.ne."enabled1") then
             iplotsxr='yes'
          else
             iplotsxr='no'
          endif
          call tdsxray(icall,iplotsxr)
        endif
        if (npa_diag.ne."disabled".and.niong.ne.0) then
           icall="first"
cBH_to skip first call, per Vincent Tang:     call tdnpadiag(icall)
cBH100816:  No.      Also, actually for npa, icall has no effect,
cBH100816:  as calculation as setup repeated each call to tdnpadiag.
           call tdnpadiag(icall)
        endif
      endif

c..................................................................
c     Call the neutral beam source module
c..................................................................

      write(*,*)'tdinitl:call frnfreya, n=',n,'time=',timet
      call frnfreya(frmodp,fr_gyrop,beamplsep,beamponp,beampoffp,
     .              hibrzp,mfm1p,noplots)
       write(*,*) 'mfm1 (freya) = ',mfm1p
c       write(*,*) 'hibrzp(i,1,1),hibrzp(i,2,1),hibrzp(i,3,1)'
c       do i=1,mfm1p
c         write(*,'(i4,2x,0p9f9.4)') i, hibrzp(i,1,1),
c     >        hibrzp(i,2,1),hibrzp(i,3,1),
c     >        hibrzp(i,1,2),hibrzp(i,2,2),hibrzp(i,3,2),
c     >        hibrzp(i,1,3),hibrzp(i,2,3),hibrzp(i,3,3)
c       enddo
c      stop
c
c     Initialize ibeampon/ibeamponp (present and previous time step)
      if (frmodp.eq."enabled") then
         ibeampon=0   !Reset in tdchief, for pulse beam
         ibeamponp=1  !Indicates beam source calculated
      else
         ibeampon=0
         ibeamponp=0
      endif

c..................................................................
c     Plot out initial radial profiles
c..................................................................
      if (noplots.ne."enabled1") call tdpltmne

c..................................................................
c     Read RF diffusion coefficient data for multiple flux
c     surfaces if rdcmod="format1" or "aorsa" or "format2"
c..................................................................

      if (rdcmod.eq."format1" .or. rdcmod.eq."aorsa" 
     +                        .or. rdcmod.eq."format2") then
        call rdc_multi
      endif

c..................................................................
c     print initial output and plot equilibrium surfaces
c..................................................................
      call tdoutput(1)
      if (noplots.ne."enabled1" .and. eqmod .eq. "enabled") call tdplteq
c$$$      if (eqmod .eq. "enabled") call tdplteq
c
      return
      end
