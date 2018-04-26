c
c

c
c..................................................................
c     ONETWO DIVERGENCE
c     THE CALLING SEQUENCE TO FREYA IS ALTERED. CIVIC WILL ONLY ALLOW
c     61 ELEMENTS IN CALLING SEQUENCE. THEREFORE THE COMMUNICATION WILL
c     BE THROUGH FRCOMM
c..................................................................

      subroutine freya(ipts,mi,mj,codeid,rin,rmax,zax,zmin,zmax)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
      integer, intent(out) :: ipts
c---------------------------------------------------------------------
c
c     this subroutine calculates the particle and energy sources
c     due to neutral beam injection.
c
c---------------------------------------------------------------------
c     BH130329:  Adding capability to read a set of birth point files
c                from NUBEAM, to compare with results from freya 
c                calculation.
c                New input variables in frsetup namelist are:
c                read_birth_pts="enabled", for use of this capability
c                birth_pts_files=input nbirth_pts_files filenames into
c                               the namelist, in the order which they
c                               will be used in cql3d. character*256.
c                               Max list length=24.  default='notset'
c                nbirth_pts_files=number of birth point files.
c                nbirth_pts= should be same as in all the NUBEAM files
c
c                A single beam is assumed, with up to three beam energy 
c                components.
c
c                The shine through, and total power is read from the
c                files. nbirth_points is checked.  Power and shine
c                through numbers may vary from file to file.
c                       
c
c---------------------------------------------------------------------
c     the input quantities to freya are:
c
c     mb              Number of neutral beam injectors (.le.kb)
c     anglev(ib)      Vertical angle (degrees) between optical axis
c     and horizontal plane; a positive value indicates
c     particles move upward
c     angleh(ib)      Horizontal angle (degrees) between optical axis and
c     vertical plane passing through pivot point and
c     toroidal axis; a zero value denotes perpendicular
c     injection, while a positive value indicates par-
c     ticles move in the co-current direction
c     nsourc          Number of sources per beamline.
c     If 1, source is centered on beamline axis.
c     If nsourc=2, distinguish between the beamline
c     axis and the source centerline (optical axis).
c     The two sources are assumed to be mirror images
c     through the beamline axis.
c     In either case, the exit grid plane is perpendicula
c     to the beamline axis, and contains the source
c     exit grid center(s).
c     If nsourc=2, the alignment of the sources w.r.t.
c     the beamline axis is specified through bhofset,
c     bvofset, and bleni (described further below).
c     bvofset(ib)     Vertical offset from beamline axis to center
c     of each source (cm; used only for nsourc=2)
c     bhofset(ib)     Horizontal offset from beamline axis to center
c     of each source (cm; used only for nsourc=2)
c     bleni(ib)       Length along source centerline (source optical axis)
c     source to intersection point with the beamline axis
c     sfrac1(ib)      Fraction of source current per beamline coming
c     from upper source (used only for nsourc=2)
c     bcur(ib)        Total current (a) in ion beam (used only if bptor
c     is zero)
c     bptor(ib)       Total power (w) through aperture into torus; when
c     nonzero, bptor takes precedence over bcur
c     bshape(ib)      Beam shape
c     'circ' : circular
c     'rect' : rectangular
c     bheigh(ib)      Height of source (cm)
c     bwidth(ib)      Width of source (cm); diameter for
c     circular source.
c     bhfoc(ib)       Horizontal focal length of source (cm)
c     bvfoc(ib)       Vertical focal length of source (cm)
c     bhdiv(ib)       Horizontal divergence of source (degrees)
c     bvdiv(ib)       Vertical divergence of source (degrees)
c     ebkev(ib)       Maximum particle energy in source (keV)
c     fbcur(ie,ib)    Fraction of current at energy ebkeV/ie
c     iborb           Flag for modeling orbit effects on beam-
c     generated fast ions
c     1, model orbit effects
c     0, do not model orbit effects
c     npart           Number of particles followed into plasma
c     (suggest 10000)
c     npskip          Ratio of number of particles followed into plasma
c     to number of source particles (suggest 1)
c     naptr           Total number of aperatures encountered by a particle
c     as is moves from the source into the plasma chamber
c     Maximum is specified by parameter nap (=10).
c     First set of aperatures encountered by the particle
c     are assumed centered on the source axis, and subseq
c     aperatures are centered on the beamline axis; the
c     distinction is made through ashape.
c     ashape(iap,ib)  Aperture shape.
c     Prefix 's-' indicates source axis centered.
c     Prefix 'b-' indicates beamline axis centered.
c     's-circ'          'b-circ'
c     's-rect'          'b-rect'
c     's-vert'          'b-vert'
c     's-horiz'         'b-horiz'
c     'b-d3d'
c     (circ=circular aperature, rect=rectagular,
c     vert=limits vertical height of source particles,
c     horiz=limits horizontal height of source particles,
c     d3d= special DIII-D polygonal aperature)
c     'circ' : circular
c     'rect' : rectangular
c     aheigh(iap,ib)  Height of aperture (cm)
c     awidth(iap,ib)  Width of aperture (cm); diameter for circular
c     aperture
c     alen(iap,ib)    Length from source to aperature for 's-type' aperatur
c     and from exit grid plane along beamline axis for
c     'b-type' aperatures.
c     blenp(ib)       Distance along beamline axis from source exit
c     plane to the fiducial "pivot" point.
c     rpivot(ib)      Radial position of pivot (cm)
c     zpivot(ib)      Axial position of pivot (cm)
c     atw(k)          atomic mass of ion species k
c     codeid          flux surface geometry
c     "onedee"      : elliptical
c     anything else : nonelliptical
c     elong           elongation (height/width) of elliptical
c     cross-section plasma
c     ibion           ion index of beam species
c     mf              number of flux zones
c     mi              number of radial mesh points for nonelliptical
c     plasma
c     mj              number of axial mesh points for nonelliptical
c     plasma
c     nion            number of ion species
c     norb            flag controlling output of orbit information
c     >0, write out extensive orbit information
c     to unit norb
c     0, omit output of orbit information
c     <0, write out abbreviated orbit information
c     to unit iabs(norb)
c     pinsid(i)       poloidal magnetic flux (G-cm2) along horizontal
c     chord inside and through the magnetic axis vs.
c     a uniform mesh in major radius; pinsid is needed
c     only if iborb=1
c     potsid(i)       poloidal magnetic flux (G-cm2) along horizontal
c     chord outside and through the magnetic axis vs.
c     a uniform mesh in major radius; potsid(1) and
c     potsid(mf) are needed if codeid.ne.'onedee';
c     the entire potsid array is needed if iborb=1
c     psi(i,j)        poloidal magnetic flux (G-cm2) at mesh point i,j
c     (needed only if codid.ne.'onedee')
c     psivol(i)       volume (cm3) of flux zone i; depending upon
c     whether codeid.eq.'onedee' or not, psivol is
c     chosen such that either r or sqrt(psi-psiax)
c     varies a constant amount from one flux surface
c     to the next
c     rinsid(i)       major radius (cm) along horizontal chord inside
c     and through the magnetic axis vs. a uniform mesh
c     in sqrt(psi-psiax); rinsid is needed only if
c     iborb=1
c     rotsid(i)       major radius (cm) along horizontal chord outside
c     and through the magnetic axis vs. a uniform mesh
c     in sqrt(psi-psiax); rotsid(1) and rotsid(mf) are
c     needed if codeid.eq.'onedee'; the entire rotsid
c     array is needed if iborb=1
c     r(i)            maj radius mesh point i for nonelliptical plasma (cm)
c     rin             major radius of inside of vacuum vessel (cm)
c     rmax            maximum radial position of plasma (cm)
c     rpivot(ib)      radial position of pivot point (cm)
c     sfrac1(ib)      fraction of source current per beamline coming
c     from upper source (used only for nsourc=2)
c     sofset(ib)      vertical offset from optical axis to center
c     of each source (cm; used only for nsourc=2)
c     ONETWO DIVERGENCE: use w instead of z
c     w(j)            axial mesh point j for nonelliptical plasma (cm)
c     zax             axial position of magnetic axis (cm) for
c     elliptical plasma
c     zmin            minimum axial position of plasma (cm)
c     zmax            maximum axial position of plasma (cm)
c     zpivot(ib)      axial position of pivot point (cm)
c     zne(i)          electron density in zone i (cm-3)
c     zni(i,k)        density of ion species k in zone i (cm-3)
c     zte(i)          electron temperature in zone i (kev)
c     zzi(i,k)        charge number of ion species k in zone i
c     ONETWO DIVERGENCE
c     zeffctv(i)      z-effective (used with iexcit=-1)
c     iexcit       switch that controls the following
c     ONETWO DIVERGENCE
c     iexcit=-1 uses Stearn's average cross sections:otw. like =0
c     iexcit = 0   (default mode) normal freya run
c     iexcit = 1   use hexnb instead of crsecs to get cross sections
c     but neglect presence of excited sates in beam
c     iexcit = 2   use hexnb with excited beam states allowed.
c     inubpat     switch controlling 2-d beam deposition calculation
c     0, (default) bypass 2-d deposition
c     1, determine beam deposition on standard (r,w) grid
c     and output deposition of 3rd excited state
c     components to file 'beamdep'
c     2, same as inubpat=1, except bin deposition on
c     modified (r,w) grid whose dimensions are given by
c     'npat' (see below).
c     note:  if inubpat.gt.0, iexcit must be set equal to 2.
c     this will automatically be checked, iexcit rese
c     (if needed), and a message output to the crt.
c     npat        modified (r,w) dimension array for optional use if
c     inubpat.gt.0, containing:
c     npat(1)=number of new 'r' elements
c     npat(2)=number of new 'w' elements
c     default:  npat(1)=mi, npat(2)=mj
c
c     the output quantities are:
c
c     bion(ie,ib)     intensity of ion beam (particles/s)
c                     [BH: evidently out of the ion source, before 
c                      neutralization]
c     bneut(ie,ib)    intensity of neutral beam (particles/s)
c                     [BH: Coming out of the neutralizer, 
c                      a fraction of bion]
c     bpow(ie,ib)     beam power to aperture (w)
c     ebeam(ie,ib)    particle energy (kev)
c     fap(ie,ib)      fraction of beam stopped by aperture
c     fwall(ie,ib)    fraction of beam incident on wall (shinethrough)
c     forb(ie,ib)     fraction of beam lost on orbits
c     ftrapfi(i,ie,ib)fraction of trapped fast ions in each zone
c     fb11(ie,ib)     fraction of ions passing and axis-encircling
c     fb10(ie,ib)     fraction of ions passing and not encircling
c     fb01(ie,ib)     fraction of ions trapped and axis-encircling
c     fb00(ie,ib)     fraction of ions trapped and not encircling
c     fber(ie,ib)     fraction of ions trapped for which error was detec
c     hibrz(i,ie,ib)  normalized hot ion birth rate
c     hdepz(i,ie,ib)  normalized hot ion deposition rate
c     ipts            number of birth points to be plotted
c     xpts(ii)          x coordinate of birth point
c     ypts(ii)          y coordinate of birth point
c     zpts(ii)          w coordinate of birth point
c     vx(ii)          x component of birth velocity
c     vy(ii)          y component of birth velocity
c     vz(ii)          w component of birth velocity
c     wb11(ie,ib)     orbit width of fb11 ions (cm)
c     wb10(ie,ib)     orbit width of fb10 ions (cm)
c     wb01(ie,ib)     orbit width of fb01 ions (cm)
c     wb00(ie,ib)     orbit width of fb00 ions (cm)
c     zetaz(i,ie,ib)  average pitch angle cosine of deposited hot ions
c     angmpz(i,ie,ib) toroidal angular momentum density deposition rate.
c     (i.e. ang. momtm. deposited in shell per second /volume of she
c     angmpz contains only the toroidal component of ang momentum.
c---------------------------------------------------------------------
c     ONETWO DIVERGENCE
      include 'param.h'
      include 'frcomm.h'

c     Automatic arrays for local dynamic mem, for case where NUBEAM
c     particle birth point list is used. nbirth_pts from frcomm.h
      real*8, allocatable :: x_nub(:),y_nub(:),z_nub(:)
      real*8, allocatable :: vx_nub(:),vy_nub(:),vz_nub(:)
      real*8, allocatable :: v_nub(:),en_nub(:)
      integer, allocatable :: ie_nub(:)  !energy cmpt,=1, 2, or 3
      character*128 filenm
      real*8, dimension(3) :: en_avg_cmpts_nub,pabs_cmpts_nub
      integer, dimension(3) :: nbirth_cmpts_nub
      real*8, dimension(3) :: birth_rate_cmpts_nub

      character*8 codeid 
      dimension psi(ki,kj),r(ki),w(kj)
      dimension bpow(ke,kb),eb(ke,kb)
      equivalence(bpow,pbeam)
      equivalence(eb,ebeam)
      equivalence(r,xxx)
      equivalence(w,yyy)
      equivalence(psi,p)
      dimension cangv(kb), cangh(kb), sangv(kb), sangh(kb),
     *  thetp(kb), thetpp(kb),
     *  costp(kb), sintp(kb), costpp(kb), sintpp(kb),
     *  vbeam(ke,kb), iatype(nap,kb)
c      dimension sgvxne(kz),sgxn(kz,ke,kb),sgxnmi(ke,kb),hxfrac(ke,kb)
      real*8 sgvxne(kz),sgxn(kcmp1,kz,kbe,ksge),
     *       sgxnmi(ke,kb),hxfrac(ke,kb)
c      real*8 sgxn1(kcmp1,kz,kbe,ksge)
      dimension wt(kz), zeta_(kz),angmtp(kz),nmbrz(kz)
      dimension rzpat(kix2,kjx2,ke,kb)
      data pio180 /0.017453293d0/
      data ncalls /0/
c
c.. new (JK)
c   kimp-> 2 for namei
c
      integer ne_tk
      real*8 fe_tk, de_tk
      real*8 sgxnloc(kbe)
      real*8 zangrot(kz), znis(kz,kion), zzis(kz,kion)
c      
      if (read_birth_pts.eq."enabled") then
      if (.NOT. ALLOCATED(x_nub)) then     
      allocate(x_nub(nbirth_pts),y_nub(nbirth_pts),z_nub(nbirth_pts))
      allocate(vx_nub(nbirth_pts),vy_nub(nbirth_pts),vz_nub(nbirth_pts))
      allocate(v_nub(nbirth_pts),en_nub(nbirth_pts))
      allocate(ie_nub(nbirth_pts))
      endif
      endif

c----------------------------------------------------------------------
c     general freya initialization
      call bcast(rzpat,zero,kix2*kjx2*ke*kb)   !rzpat not used w cql3d
      mb=nbeams
      elong=0.d0
      norb=0
c
      do ib=1,ke
       do i=1,kb
         bneut(ib,i)=0.d0
         vbeam(ib,i)=0.d0
         sgxnmi(ib,i)=0.d0
         hxfrac(ib,i)=0.d0
       enddo
      enddo
      do i=1,kz
       do j=1,ke
         do k=1,kb
           hibrz(i,j,k)=0.d0
           hdepz(i,j,k)=0.d0
           ftrapfi(i,j,k)=0.d0
         enddo
       enddo
      enddo
      do i=1,kz
       do j=1,ke
         do k=1,kb
          sgxn(1,i,j,k)=0.d0
          sgxn(2,i,j,k)=0.d0
          sgxn(3,i,j,k)=0.d0
          sgxn(4,i,j,k)=0.d0
        enddo
       enddo
      enddo
      pzone=0.D0
      rzone=0.D0
c
c----------------------------------------------------------------------
c     set up some data for optional (r,w) deposition calculation
      inubpat=0  !Not used with cql3d
      if(inubpat.gt.0) then
        call setrz(npat,r(1),r(mi),w(1),w(mj),drpat,dzpat,nrpat,nzpat)
        iexcit=2
      endif

c----------------------------------------------------------------------
c     if nubeam list case, make sure only 1 beam
      if (read_birth_pts.eq."enabled") then
         nbeams=1
         mb=1
         if (codeid.eq.'onedee') then
            WRITE(*,*)'Freya: nubeam list not setup for this codeid'
            WRITE(*,*)'STOP'
            STOP
         endif
      endif

c
c     turn off orbit calculation in case of counter-injection; present
c     orbit model applies only for co-injection
      do 5 ib=1,mb
 5    if(angleh(ib).lt.5.) iborb = 0
c
c     initialize flux surface quantities
c     BH130327:
c          rotsid, rinsid equal 0, potsid(mf),potsid(1)set but other 0.
      mfm1 = mf-1
      rmajor = rotsid(1)
      drot   = (rotsid(mf)-rotsid(1))/mfm1
      drin   = (rinsid(1)-rinsid(mf))/mfm1
      drutp  = sqrt(potsid(mf)-potsid(1))/mfm1
c990131      drot   = amax1(drot,1.e-6)
c990131      drin   = amax1(drin,1.e-6)
c990131      drutp  = amax1(drutp,1.e-6)
c990131      elong  = amax1(elong,1.e-6)
      deps=1.d-6
      drot   = max(drot,deps)
      drin   = max(drin,deps)
      drutp  = max(drutp,deps)
      elong  = max(elong,deps)   !elong not properly set.
      droti  = one/drot
      drini  = one/drin
      drutpi = one/drutp
      elongi = one/elong
      if(codeid.ne.'onedee') then
        mim1  = mi-1
        mjm1  = mj-1
        dr    = (r(mi)-r(1))/mim1
        dz    = (w(mj)-w(1))/mjm1
        dri   = one/dr
        dzi   = one/dz
      endif
c
c     calculate sines and cosines of various angles
      do 15 ib=1,mb
        cangv(ib)=cos(anglev(ib)*pio180)
        cangh(ib)=cos(angleh(ib)*pio180)
        sangv(ib)=sin(anglev(ib)*pio180)
        sangh(ib)=sin(angleh(ib)*pio180)
 15   continue
c
c     calculate beam power; account for neutralizer efficiency
c
      atwb = atw(ibion)
      bntot = zero
      do 20 ib=1,mb
        do 21 ie=1,3
          ebeam(ie,ib) = ebkev(ib)/ie
          ebx = ebeam(ie,ib)/atwb
c     ONETWO DIVERGENCE   - delete reference to ionization efficiency
c     routines.
      call logint(ebx,beff)  ! JK  BH140501: Check if any diff.
c          beff=1.
          ebev = 1.e3*ebx
          vbeam(ie,ib) = 1.384e6*sqrt(ebev)  !adjstd for mass, cm/sec
          bion(ie,ib)  = 0.625e19*fbcur(ie,ib)*bcur(ib) ! molecular src
                                       !rate (as I understand, BH)
                                       !bcur need not be set 
                                       !in nml, if bptor.ne.0.
                                       !In this case, bneut, bpow are
                                       !set for default bcur()=110.,
                                       !and renormalize for bptor below.
          bneut(ie,ib) = ie*beff*bion(ie,ib)  !ion rates at energy/ie.
          bpow(ie,ib)  = ebeam(ie,ib)*bneut(ie,ib)/0.625e16
c          write(*,*)'freya: ibion,atw(ibion),ie,ib,vbeam(ie,ib)',
c     &                      ibion,atw(ibion),ie,ib,vbeam(ie,ib)
 21     bntot = bntot + bneut(ie,ib)
 20   continue
c
c.......................................................................
cBH130329:  If read_birth_pts="enabled", then read in beam birth points 
c           from a NUBEAM generated file, as described above.
c           There may be some unnecessary freya related calculation
c           since the NUBEAM data-read is shoe-horned in on top of freya
c
c           Else:  freya calculation.
c.......................................................................

      if (read_birth_pts.eq."enabled") then
c        read_nubeam_data reads data file and returns it
c        through the argument list.
c        NOTE:  The NUBEAM birth points are actually shifted from
c        the ion birth point, accounting for the distance to the
c        GC starting position of the particle.  That is, these
c        are birth points for particle guiding centers.
         ncalls=ncalls+1
         if (ncalls.gt.nbirth_pts_files) then
            write(*,*)'FREYA WARNING: Insufficient birth_pts_files'
            write(*,*)'FREYA WARNING: Stepping back to last file'
            ncalls=ncalls-1
         else
            write(*,*)'Freya: read_birth_pts case, ncalls=',ncalls   
         endif
         filenm=trim(birth_pts_files(ncalls))
         call read_nubeam_data(filenm,nbirth_pts,atwb,
     +                   nbirth_cmpts_nub,nshine_nub,
     +                   pinj_nub,pabs_nub,x_nub,y_nub,z_nub,
     +                   vx_nub,vy_nub,vz_nub,v_nub,en_nub,ie_nub,
     +                   birth_rate_cmpts_nub,en_avg_cmpts_nub,
     +                   pabs_cmpts_nub)

c        pinj_nub=Injected power (Watts) = freya bptor
         bptor(1)=pinj_nub
c        Set bcur,bion,bneut,bpow from Nubeam list
c        This is setup according to freya bptor, as above.
c        Here, bcur, here, refers to absorbed power, since
c        present Nubeam list doesn't contain a breakdown of the shine
c        through components.  [Fix later].
         bntot=0d0
         ib=1                   !Nubeam list only set up for 1 beam
         do ie=1,3
            bcur(ib)= bcur(ib)+birth_rate_cmpts_nub(ie)*1.6022e-19
            bneut(ie,ib)=birth_rate_cmpts_nub(ie)
            bion(ie,ib)=bneut(ie,ib)/ie
            bpow(ie,ib)=pabs_cmpts_nub(ie)
            bntot=bntot+bpow(ie,ib)  !Prob here, compared to freya
                                     !calc, since for freya should be
                                     !same units as bneut[see below use]
         enddo

      endif    ! On read_birth_pts.eq."enabled"

      if (read_birth_pts.ne."enabled") then  !Skip, if nubeam case
      
c     ONETWO DIVERGENCE
c..................................................................
c     Determine peak electron density, temp and zeff(lr_) for use with
c     iexcit=-1 (Stearn's formula)
c..................................................................

      call aminmx(zne,1,mfm1,1,dnemin,dnemax,kmin,kmax)
      call aminmx(zte,1,mfm1,1,dtemin,dtemax,kmin,kmax)
      call aminmx(zeffctv,1,mfm1,1,dzemin,dzemax,kmin,kmax)
      dtemax=dtemax/10.
c
c...calculate macroscopic cross sections
c     iexcit .le. 0 : Freeman and Jones
c                 5 : ADAS
c                 6 : Boley parameterization
c
c... new coding (JK)
c
      csgn=1.d0
      do i=1,kz
       zangrot(i)=0.d0
      enddo
c
c... In CQL3D, for nprim=1, nimp=1, main ions in zni(j,1), (j,2)
c    and zzi(j,1) and zzi(j,2) with general and maxwellian
c    distributions, respectively. The impurity ion is in 
c    zni(j,3) and zzi(j,3). Also, Zeff has been copied into
c    zzi(j,nion+2=4). So, for this case we slip
c    the arrays such that (j,3) - > (j,2) for the impurity
c
      if(nprim.eq.1 .and. nimp.eq.1) then
       do j=1,kz
         znis(j,1)=zni(j,2)
         znis(j,2)=zni(j,3)
         zzis(j,1)=zzi(j,2)
         zzis(j,2)=zzi(j,3)
         zzis(j,3)=zzi(j,4)  ! Zeff profile
       enddo
      else
       do j=1,kz
         do i=1,kion
           znis(j,i)=zni(j,i)
           zzis(j,i)=zzi(j,i)
         enddo
       enddo
      endif
c
      write(*,*) 'calling nbsgxn ...'
      write(*,*) 'iexcit = ',iexcit
c       write(*,*) 'zte:'
c       write(*,'(5(2x,1pe14.5))') (zte(j),j=1,kz)
c       write(*,*) 'zti:'
c       write(*,'(5(2x,1pe14.5))') (zti(j),j=1,kz)
c       write(*,*) 'zne:'
c       write(*,'(5(2x,1pe14.5))') (zne(j),j=1,kz)
c       write(*,*) 'zni(j,1-3):'
c       write(*,'(5(2x,1pe14.5))') (znis(j,1),j=1,kz)
c       write(*,'(5(2x,1pe14.5))') (znis(j,2),j=1,kz)
c       write(*,'(5(2x,1pe14.5))') (znis(j,3),j=1,kz)
c       write(*,*) 'zzi(j,1-3):'
c       write(*,'(5(2x,1pe14.5))') (zzis(j,1),j=1,kz)
c       write(*,'(5(2x,1pe14.5))') (zzis(j,2),j=1,kz)
c       write(*,'(5(2x,1pe14.5))') (zzis(j,3),j=1,kz)
c       stop
c
        call nbsgxn (iexcit,namep,namei,mb,mfm1,ne_tk,nprim,
     .             nimp,nion,atwb,atw,ebkev,fe_tk,ibion,vbeam,
     .             zne,znis,zte,zti,zzis,de_tk,dtemax,dnemax,dzemax,
     .             hxfrac,sgxn,sgxnmi)
c      do ib=1,mb
c       do j=1,3
c         write(*,'(2i3,2x,1p1e12.4,a12)') j,ib,sgxnmi(j,ib),' sgxnmi'
c       enddo
c      enddo
c      stop
c
c... old coding (JK)
c
c      if(iexcit.le.0)then
c     ONETWO DIVERGENCE
c        call crsecs(atw,ebkev,ibion,ke,kz,mb,mfm1,nion,vbeam,
c     &    zne,zni,zte,zzi,sgvxne,sgxn,sgxnmi,iexcit,
c     &    dtemax,dnemax,dzemax)
c         do ib = 1, mb
c           do j = 1, 3
c             do i=1,mfm1
c        write(*,'(3i3,2x,1p1e10.4,a12)') i,j,ib,sgxn(i,j,ib), ' sgxn'
c             enddo
c           enddo
c         enddo
c         stop
c      else
c     ONETWO DIVERGENCE
c         write(*,*) 'freya: before frhexdrv'
c        call frhexdrv(mb,sgxn,sgxnmi,hxfrac)
c      endif
c... end old coding
c
      endif  !on read_birth_pts.ne."enabled"
c
c     calculate total plasma volume
      volume = 0.
      do i=1,mfm1
        write(*,*)'freya: i,psivol(i)=',i,psivol(i)
        volume = volume + psivol(i)
      enddo
      write(*,*)'freya: volume(sum of psivol)=',volume
      !Note: psivol is based on eqvol() which is calc-ed by subr.eqvolpsi.
      !Can be a little different from alternative definition
      ! of flux surface volume as setup in subr.tdrmshst. 
c----------------------------------------------------------------------
c     begin loop over beams, ib
c----------------------------------------------------------------------

c     Following statment executed in subroutine frset.
c      if (read_birth_pts.eq."enabled") npart=nbirth_pts  !nubeam case

cBH110309      iskip =       1 + (npart-1)/1500
      maxp=1500000 !Formerly parameter giving max # of ions launched, 
                  !according to code comment
      iskip =       1 + (npart-1)/maxp
      ic    = 0
      ipts  = 0
         
      do 200 ib=1,mb  ! 386 lines down to line 989

      if (read_birth_pts.ne."enabled") then  !Skip, if nubeam case
c
c     Determine aperature designators.
c     If nsourc=1, treat source axis centered aperatures equivalent
c     to beam axis centered.
cBH130915: Not sure of effect of this change. Reverting:do 40  i=1,naptr
	 do 40  i=1,nap
            if (nsourc.eq.1) then
               if(ashape(i,ib).eq.'s-circ')  iatype(i,ib)=5
               if(ashape(i,ib).eq.'s-rect')  iatype(i,ib)=6
               if(ashape(i,ib).eq.'s-vert')  iatype(i,ib)=7
               if(ashape(i,ib).eq.'s-horiz') iatype(i,ib)=8
            elseif (nsourc.gt.1) then
               if(ashape(i,ib).eq.'s-circ')  iatype(i,ib)=1
               if(ashape(i,ib).eq.'s-rect')  iatype(i,ib)=2
               if(ashape(i,ib).eq.'s-vert')  iatype(i,ib)=3
               if(ashape(i,ib).eq.'s-horiz') iatype(i,ib)=4
            endif
            if(ashape(i,ib).eq.'b-circ')  iatype(i,ib)=5
            if(ashape(i,ib).eq.'b-rect')  iatype(i,ib)=6
            if(ashape(i,ib).eq.'b-vert')  iatype(i,ib)=7
            if(ashape(i,ib).eq.'b-horiz') iatype(i,ib)=8
            if(ashape(i,ib).eq.'b-d3d')   iatype(i,ib)=9
c        write(*,*) 'ashape(i,ib) = ',i,ib,ashape(i,ib)
c        write(*,*) 'iatype(i,ib) before rotate = ',iatype(i,ib)
 40      continue
c     
c     Some angles for subroutine rotate
        thetp(ib)=atan2(bvofset(ib),sqrt(bleni(ib)**2-bvofset(ib)**2))
        costp(ib)=cos(thetp(ib))
        sintp(ib)=sin(thetp(ib))
        thetpp(ib)=atan2(bhofset(ib),sqrt(bleni(ib)**2-bvofset(ib)**2))
        costpp(ib)=cos(thetpp(ib))
        sintpp(ib)=sin(thetpp(ib))
c
        endif  !on read_birth_pts.ne."enabled"
c
c----------------------------------------------------------------------
c     begin loop over beam energy components, ie
c----------------------------------------------------------------------
        do 201 ie=1,3  !340 lines down
c         beam fractions
          fap(ie,ib)   = 0.
          fwall(ie,ib) = 0.
          forb(ie,ib)  = 0.
          fb11(ie,ib)  = 0.
          fb10(ie,ib)  = 0.
          fb01(ie,ib)  = 0.
          fb00(ie,ib)  = 0.
          fber(ie,ib)  = 0.
          do 110 i=1,mfm1
             ftrapfi(i,ie,ib)=0.
             hibrz(i,ie,ib) = 0.
             hdepz(i,ie,ib) = 0.
             angmpz(i,ie,ib) = 0.
             hicmz(i,ie,ib,1) = 0.
             hicmz(i,ie,ib,2) = 0.
             hicmz(i,ie,ib,3) = 0.
             zetaz(i,ie,ib) = 0.
 110      continue
          do 310 i=1,mfm1
             nmbrz(i)=0
 310      continue
c     --- nmbrz(i) counts ions born in zone i
c----------------------------------------------------------------------
c     begin loop over particles, for each beam(ib),energy(ie)
c----------------------------------------------------------------------
c
          if (read_birth_pts.ne."enabled") then
             npar = (bneut(ie,ib)/bntot)*npart   !Not nubeam list
          else
             npar=nbirth_cmpts_nub(ie)           !nubeam list
             ii=0
          endif
c
cBH130925          npar = (bneut(ie,ib)/bntot) * npart
          if(npar.eq.0) go to 201
          nparx = 0
          newpar = 0  ! JK
c
c.......................................................................
c     Loop over the particles, for each ie,ib
c.......................................................................
         do 180 ipar=1,npar  ! down 218 lines, to line 919
            if(mod(ipar-1,npskip).eq.0) newpar=1
c YuP[171103] Moved this line inside if()            if(newpar.eq.0) go to 120
c
          if (read_birth_pts.ne."enabled") then  !Skip if nubeam list
             ! FREYA normal calculations
             if(newpar.eq.0) go to 120 ! YuP[171103] moved from above
c
c... generate neutral particle at beam source
c
             call sorspt1(bshape,bheigh,bwidth,bhfoc,bvfoc,bhdiv,
     .            bvdiv,ib,ie,isourc,nsourc,sfrac1,vbeam,x0,y0,z0,
     .            vx0,vy0,vz0)
c
c... transform coordinates and advance particle to pivot point
c
            call rotate(naptr,iatype,aheigh,awidth,alen,bhofset,
     &        bvofset,cangv,cangh,ib,isourc,costp,sintp,costpp,
     &        sintpp,blenp,nsourc,sangv,sangh,rpivot,zpivot,mlost,
     &        x0,y0,z0,vx0,vy0,vz0)
c      write(*,*)'freya,rotate:x0,y0,z0,vx0,vy0,vz0',x0,y0,z0,vx0,vy0,vz0
c
c     skip injection if particle is lost at aperture
cCheck inj  write(*,*)'freya: ie,ipar,npar,mlost',ie,ipar,npar,mlost
c 120        if(mlost.ne.0) go to 160
  120       continue
            if (mlost .ne. 0)  go to 160
c
c... inject particle into plasma, i.e., follow particle from pivot
c     point into, through, or around plasma. Use old subroutine
c     inject for crsecs.
c
c       open(42,file='psi.dat',status='unknown')
c       write(42,*) 'psi data for DIII-D #122080 from cql3d'
c       write(42,*) '(nw,nh) = ',65,65
c       do i=1,65
c        do j=1,65
c          write(42,'(2i4,2x,1p1e15.7,a16)') i,j,psi(i,j),' psi'
c        enddo
c       enddo
c       close(42)
c       stop
c
            if (iexcit.le.0) then
             call inject_old(atw,codeid,drutpi,droti,dri,dzi,
     &             elongi,ib,ie,mfm1,mim1,mjm1,newpar,potsid(1),
     &             psi,r,rmajor,rin,rmax,sgxn,sgxnmi,x0,y0,z0,
     &             vx0,vy0,vz0,vbeam,w,zax,zmin,zmax,izone,
     &             pzone,rzone,rpos,xpos,ypos,zpos)
            else
             call inject1(atwb,codeid,de_tk,drutpi,droti,dri,ds_tk,dzi,
     &             elongi,ib,ie,kb,kbe,ksge,ke,kz,ki,mfm1,mim1,
     &             mjm1,ne_tk,newpar,nout,potsid(1),psi,r,rmajor,rin,
     &             rmax,sgxn,sgxnloc,sgxnmi,x0,y0,z0,vx0,vy0,vz0,
     &             vbeam,w,zangrot,zax,zmin,zmax,izone,
     &             pzone,rzone,rpos,xpos,ypos,zpos,myid,tenter,
     &             smax,texit)
            endif

c           Shift from particle birth point to guiding center point


          else    !On read_birth_pts.ne."enabled"  ! NUBEAM list

c           Nubeam list case.
c           Remember, for given ie (energy cmpt),
c           need to skip to next particle for the given ie.
            do
              ii=ii+1
              if (ie_nub(ii).eq.ie) exit
              if (ii.gt.npart) then
                WRITE(*,*)'Freya: Nubeam list inconsistency'
                STOP
              endif
            enddo
c
            xpos=x_nub(ii)
            ypos=y_nub(ii)
            zpos=z_nub(ii)
            vx0=vx_nub(ii)
            vy0=vy_nub(ii)
            vz0=vz_nub(ii)
            rpos=sqrt(xpos**2+zpos**2)
c           get pzone (i.e., psi-value) and radial izone
            call zone(drutpi,ki,mfm1, mim1, mjm1,dri,dzi,potsid(1),psi,
     +           r,w,xpos,ypos,zpos,pzone,izone)
c         
            vbeam(ie,ib)=v_nub(ii)
c
          endif                  !On read_birth_pts.ne."enabled"
c
c... skip birth data if:  particle missed plasma
c
            if(izone.ge.mf) go to 170
c     For removing NBI source at all psi outside of psicutoff: 
          !print*,'Rpos,Zpos, pzone,psicutoff=',rpos,zpos,pzone,psicutoff
          !pause
            if(psicutoff.ne.0.d0  .and. (-pzone).lt.psicutoff) go to 170
            ! pzone is increasing from center to edge.
            ! equilpsi is decreasing, and psicutoff is based on equilpsi
c
c... accumulate hot ion birth rate
c
            hibrz(izone,ie,ib) = hibrz(izone,ie,ib) + one
            hicmz(izone,ie,ib,1) = hicmz(izone,ie,ib,1) + sgxnloc(1)
            hicmz(izone,ie,ib,2) = hicmz(izone,ie,ib,2) + sgxnloc(2)
            hicmz(izone,ie,ib,3) = hicmz(izone,ie,ib,3) + sgxnloc(3)
c           write(*,*) 'hibrz(1,1,2) = ',hibrz(1,1,2)
c          write(*,*) 'hibrz = ',hibrz(izone,ie,ib)
c          write(*,*) 'hicmz-1 = ',hicmz(izone,ie,ib,1)
c
c... Calculate pitch angle cosine at birth point; accumulate average
c    pitch angle cosine. Calculate toroidal angular momentum deposited
c    in shell by each monte carlo ion.
c990131            zetai = amin1(zetai,1.)
c990131            zetai = amax1(zetai,-1.)
            zetai = (-ypos*vx0+xpos*vy0)/(rpos*vbeam(ie,ib))    ! Original
            vplane = SQRT (vx0**2+vy0**2)                       ! Onetwo
cBH130914            zetai  = csgn*(xpos*vy0-ypos*vx0)/(rpos*vplane)     ! Onetwo
            zetai = min(zetai,one)
            zetai = max(zetai,-one)
            vtroid=zetai*vbeam(ie,ib)                           ! Original
cBH130914            vtroid = csgn*zetai*vplane                          ! Onetwo
            vrad   = -SQRT (1.0 - zetai**2) * vplane            ! Onetwo
            angmtm=rpos*vtroid*atwb*1.673e-24
            angmpz(izone,ie,ib) = angmpz(izone,ie,ib) + angmtm  ! Onetwo
            if(iborb.ne.1) then
              zetaz(izone,ie,ib) = zetaz(izone,ie,ib) + zetai
c              angmpz(izone,ie,ib)=angmpz(izone,ie,ib) + angmtm
            endif
c
c... save occasional birth point for subsequent plotting
c
            ic = ic + 1
            if(mod(ic-1,iskip).eq.0) then
              ipts = ipts+1
              xpts(ipts) = xpos
              ypts(ipts) = ypos
              zpts(ipts) = zpos
              rpts(ipts)=sqrt(xpos**2+ypos**2)
              vx(ipts) = vx0
              vy(ipts) = vy0
              vz(ipts) = vz0
            endif
c
c... calculate (r,w) grid location
c
            if(inubpat.gt.0)then
              i=(rpos-r(1))/drpat + one
              j=(zpos-w(1))/dzpat + one
              rzpat(i,j,ie,ib) = rzpat(i,j,ie,ib) + one
            endif
c
c ----------------------------------------------------------------------
c... calculate orbit widths and orbit loss
c    [BH140501:  Might be good to replace this with the G.C. orbit-based
c                calcs.
c
            if(iborb.ne.0)then
              nparx  = nparx + 1
              call freyorb(atwb,codeid,drutpi,drini,droti,ic,iskip,
     +          izone,mfm1,
     &          norb,pinsid,potsid,pzone,rinsid,rotsid,rzone,rpos,
     &          zetai,vbeam(ie,ib),zpos,ipass,iaxis,ier,izp,wid,wt,
     &          zeta_,angmtp)
c
              if(ier.ne.0) then
                fber(ie,ib) = fber(ie,ib) + one
              elseif(izp.gt.mfm1) then
                go to 175
              else
                nmbrz(izp)=nmbrz(izp)+1
                if(ipass.eq.1 .and. iaxis.eq.1) then
                  fb11(ie,ib) = fb11(ie,ib) + one
                  wb11(ie,ib) = wb11(ie,ib) + wid
                elseif(ipass.eq.1 .and. iaxis.eq.0) then
                  fb10(ie,ib) = fb10(ie,ib) + one
                  wb10(ie,ib) = wb10(ie,ib) + wid
                elseif(ipass.eq.0 .and. iaxis.eq.1) then
                  fb01(ie,ib) = fb01(ie,ib) + one
                  wb01(ie,ib) = wb01(ie,ib) + wid
                  ftrapfi(izp,ie,ib)=ftrapfi(izp,ie,ib) + one
                elseif(ipass.eq.0 .and. iaxis.eq.0) then
                  fb00(ie,ib) = fb00(ie,ib) + 1.
                  wb00(ie,ib) = wb00(ie,ib) + wid
                  ftrapfi(izp,ie,ib)=ftrapfi(izp,ie,ib)+1.
                endif
              endif
c
c... accumulate hot ion deposition rate and average pitch angle cosine
c
              do 155 i=1,mfm1
                hdepz(i,ie,ib) = hdepz(i,ie,ib) + wt(i)
                angmpz(i,ie,ib)=angmpz(i,ie,ib) + wt(i)*angmtp(i)
 155          zetaz(i,ie,ib) = zetaz(i,ie,ib) + wt(i)*zeta_(i)
c
            endif    ! iborb=1
c ----------------------------------------------------------------------
c
            go to 180
c
c... accumulate particles that are not deposited in plasma
c
 160        fap(ie,ib) = fap(ie,ib) + one
            go to 180
 170        fwall(ie,ib) = fwall(ie,ib) + one
            go to 180
 175        forb(ie,ib) = forb(ie,ib) + one
c----------------------------------------------------------------------
c     end loop over particles
c----------------------------------------------------------------------
c
            newpar = 0
 180     continue

c
c... normalize average pitch angle cosine. normalize momentum and birth
c    mode in each shell to a single particle.
c
            write(*,*) '** end loop over particles **'
c          do i=1,mfm1
c        write(*,'(3i4,2x,1p1e10.4,a16)')i,ie,ib,hibrz(i,ie,ib),' hibrz'
c          enddo
c          stop
c
          do 182 i=1,mfm1
            xnorm = hibrz(i,ie,ib)
            if (xnorm.ne.0) then
              angmpz(i,ie,ib)  = angmpz(i,ie,ib)/xnorm   ! Onetwo
              hicmz(i,ie,ib,1) = hicmz(i,ie,ib,1)/xnorm  ! Onetwo
              hicmz(i,ie,ib,2) = hicmz(i,ie,ib,2)/xnorm  ! Onetwo
              hicmz(i,ie,ib,3) = hicmz(i,ie,ib,3)/xnorm  ! Onetwo
            endif
            if(iborb.ne.0) xnorm = hdepz(i,ie,ib)
            if(xnorm.ne.zero) xnorm = 1./xnorm
 182      zetaz(i,ie,ib) = xnorm*zetaz(i,ie,ib)
c
c... get the fraction of trapped ions in each zone. if orbit effects are 
c    not turned on or itrapfi=0 then this effect is not included.
c
          do 300 i=1,mfm1
            nmbrz(i) = MAX0(nmbrz(i),1)
            ftrapfi(i,ie,ib) = ftrapfi(i,ie,ib)/nmbrz(i)
 300      continue
c
c... get true toroidal ang. mtm. density deposition rate by
c    multiplying by weight of each m.c. particle
c
          do 183 i=1,mfm1
            angmpz(i,ie,ib)=angmpz(i,ie,ib)*bneut(ie,ib)/
     +        (psivol(i)*npart)
 183      continue
c
c... normalize loss fractions and hot ion birth rate
c
          fap(ie,ib)   = fap(ie,ib) / npar
          fwall(ie,ib) = fwall(ie,ib) / npar
          forb(ie,ib)  = forb(ie,ib) / npar
          xloss1 = fap(ie,ib) + fwall(ie,ib)
          xloss2 = fap(ie,ib) + fwall(ie,ib) + forb(ie,ib)
          if(xloss1.ge.1.) go to 201
c
c         do i = 1, mfm1
c         write(*,'(i4,2x,1p2e12.4,a12)') i, psivol(i), hibrz(i,ie,ib),
c     .        ' hibrz-pre'
c         enddo
          do 190 i=1,mfm1
            hibrz(i,ie,ib) = hibrz(i,ie,ib)*volume
     +        /((one-xloss1)*npar*psivol(i))
            if(iborb.eq.0) then
              hdepz(i,ie,ib) = hibrz(i,ie,ib)
            elseif(xloss2.lt.1.) then
              hdepz(i,ie,ib) = hdepz(i,ie,ib)*volume
     +          /((one-xloss2)*npar*psivol(i))
            endif
 190      continue
c      if(ie.eq.1 .and. ib.eq.2) stop
c
c... normalize orbit widths and fractions
c
          if(iborb.eq.0) go to 201
          if(fb11(ie,ib).ne.zero) wb11(ie,ib) = wb11(ie,ib)/fb11(ie,ib)
          if(fb10(ie,ib).ne.zero) wb10(ie,ib) = wb10(ie,ib)/fb10(ie,ib)
          if(fb01(ie,ib).ne.zero) wb01(ie,ib) = wb01(ie,ib)/fb01(ie,ib)
          if(fb00(ie,ib).ne.zero) wb00(ie,ib) = wb00(ie,ib)/fb00(ie,ib)
          fb11(ie,ib) = fb11(ie,ib)/nparx
          fb10(ie,ib) = fb10(ie,ib)/nparx
          fb01(ie,ib) = fb01(ie,ib)/nparx
          fb00(ie,ib) = fb00(ie,ib)/nparx
          fber(ie,ib) = fber(ie,ib)/nparx
c
 201    continue  ! end loop over energy components
        write(*,*) '** end loop over energy components **'
 200  continue  ! end loop over beams
c
c----------------------------------------------------------------------
c
c     end loop over beams and components
c
c----------------------------------------------------------------------
      write(*,*) 'end loop over beams...'
      write(*,*) 'ie,ib,ipar = ',ie,ib,ipar
c
      if (read_birth_pts.ne."enabled") then
c
       write(*,*) 'hibrz(i,1,1),hibrz(i,2,1),hibrz(i,3,1)'
       do i=1,mfm1
         write(*,'(i4,2x,0p9f9.4)') i, hibrz(i,1,1),
     >        hibrz(i,2,1),hibrz(i,3,1),
     >        hibrz(i,1,2),hibrz(i,2,2),hibrz(i,3,2),
     >        hibrz(i,1,3),hibrz(i,2,3),hibrz(i,3,3)
       enddo
c       stop
c
c... renormalize currents and powers to bptor
c
      do 240 ib=1,mb
        if(bptor(ib).gt.0.) then
          bptorx = 0.
          do 210 ie=1,3
 210      bptorx = bptorx + (1.-fap(ie,ib))*bpow(ie,ib)
          if(bptorx.gt.0.) then
            xnorm = bptor(ib)/bptorx
            bcur(ib) = xnorm*bcur(ib)
            do 220 ie=1,3
              bion(ie,ib) = xnorm*bion(ie,ib)
              bneut(ie,ib) = xnorm*bneut(ie,ib)
 220        bpow(ie,ib) = xnorm*bpow(ie,ib)
            do 230 ie=1,3
              do 231 i=1,mfm1
 231          angmpz(i,ie,ib)=angmpz(i,ie,ib)*xnorm
 230        continue
          endif
        endif
 240  continue
c
      else
c
c     Set bcur,bion,bneut,bpow from Nubeam list [above]
c
      endif  !On read_birth_pts.ne."enabled"
c
c     write hot ion birth rate, and lost franction:
c     BH130407:  hibrz is not used in cql3d, except for
c       following printout.  For that purpose, it will be useful
c       (in the future) to adjust the calculation of izone, perhaps
c       to base it on the rho coordinate and rotsid,  and to
c       fill in the rotsid() array.
      write(*,*)
      write(*,*)'freya, hot ion birth rate vs rotsid, for ib=1, ie=1:3'
      write(*,*) (rotsid(i),i=1,mfm1) 
      write(*,*) (hibrz(i,1,1),i=1,mfm1)
      write(*,*) (hibrz(i,2,1),i=1,mfm1)
      write(*,*) (hibrz(i,3,1),i=1,mfm1)
c
      write(*,*)
      write(*,*) 'fap(ie=1:3)',(fap(i,1),i=1,3)
      write(*,*) 'fwall(ie=1:3)',(fwall(i,1),i=1,3)
      write(*,*) 'forb(ie=1:3)',(forb(i,1),i=1,3)
      write(*,*)
c      write(*,*)
c      write(*,*)'freya:i,xpts(i),ypts(i),zpts(i),rpts(i),vx(i),vy(i),vz(i)'
c      do i=1,5
c         write(*,*) i,xpts(i),ypts(i),zpts(i),rpts(i),vx(i),vy(i),vz(i)
c         write(*,*) i,xpts(i),ypts(i),zpts(i),rpts(i),vx(i),vy(i),vz(i)
c      enddo
c
c     calculate neutral beam deposition density on (r,w) grid
      if (inubpat.eq.1 .and. codeid.ne.'onedee') then
        call frnbdep2(psi,mi,mj,r,w,potsid,mf,rzpat,nrpat,nzpat,
     &  ke,mb,sgxn,vbeam,hxfrac,inubpat)
      endif
c
      return
      end subroutine freya


      subroutine read_nubeam_data(filenm,nbirth_pts,atwb,
     +                   nbirth_cmpts_nub,nshine_nub,
     +                   pinj_nub,pabs_nub,x_nub,y_nub,z_nub,
     +                   vx_nub,vy_nub,vz_nub,v_nub,en_nub,ie_nub,
     +                   birth_rate_cmpts_nub,en_avg_cmpts_nub,
     +                   pabs_cmpts_nub)
c
c     Read NUBEAM data file
c
      implicit integer (i-n), real*8 (a-h,o-z)

c
c     input: filenm, nbirth_pts, atwb
c     output: nbirth_cmpts_nub,nshine_nub,pinj_nub,pabs_nub,
c             x_nub,y_nub,z_nub,vx_nub,vy_nub,vz_nub,
c             en_nub,ie_nub,birth_rate_cmpts_nub,en_avg_cmpts_nub,
c             pabs_cmpts_nub
c
c     Of the nbirth_pts (frset namelist, checked against nubeam list
c     value): nbirth_cmpts_nub(1:3) is breakdown versus energy cmpt
c     In addition, nshine_nub particles shine through (from nubeam list).
c     [Breakdown into energy components, presumably most full energy,
c      is not provided.]
c     pinj_nub is injected power into the torus(Watts)
c     pabs_nub is absorbed (Watts)
c     en_nub/ie_nub(1:nbirth_pts) is energy/energy-component for each
c        particle in the list
c     birth_rate_cmpts_nub,en_avg_cmpts_nub,pabs_cmpts_nub(1:3) are the
c        neutral particle birth rates, average energy (keV) and power 
c        absorbed for each of the three energy components.
c 
  
      character*128 filenm
      real*8 atwb
      real*8, dimension(nbirth_pts) :: x_nub,y_nub,z_nub
      real*8, dimension(nbirth_pts) :: vx_nub,vy_nub,vz_nub,v_nub,en_nub
      integer, dimension(nbirth_pts) :: ie_nub
      real*8, dimension(3) :: en_avg_cmpts_nub,pabs_cmpts_nub
      integer, dimension(3) :: nbirth_cmpts_nub
      real*8, dimension(3) :: birth_rate_cmpts_nub,birth_rate_cmpts_nub1

      iunit=14
      open(unit=iunit,file=filenm,status='old',iostat=kode)
      if (kode.ne.0) then
         WRITE(*,*)
     +        'NUBEAM file', filenm,' cannot be opened'
         STOP
      endif

      read(iunit,100) pinj_nub
      read(iunit,101) pabs_nub
      read(iunit,102) nbirth_pts_nub,total_source_rate
      read(iunit,103) nshine_nub
      read(iunit,104)

 100  format(///////,37x,e13.6)   !Skip 7 lines plus 37 columns.
 101  format(37x,e13.6)
 102  format(/////////,3x,i8,15x,e13.6)
 103  format(19x,i8)
 104  format(///)
 105  format(6e14.7)

c     Check n_birth_pts_nub against namelist value
      if (nbirth_pts_nub .ne. nbirth_pts) then
         WRITE(*,*)'STOP: inconsistency with NUBEAM file n_birth_pts'
         WRITE(*,*)'Namelist n_birth_pts =',nbirth_pts
         WRITE(*,*)'NUBEAM   n_birth_pts =',nbirth_pts_nub
         STOP
      endif

c     Read the NUBEAM birth points
      ergtkev=1.6022d-09
      do i=1,nbirth_pts
         read(iunit,105)x_nub(i),y_nub(i),z_nub(i),
     +                  vx_nub(i),vy_nub(i),vz_nub(i)
c        Calc beam velocity and energy
         v_nub(i)=vx_nub(i)**2+vy_nub(i)**2+vz_nub(i)**2
         en_nub(i)=0.5*atwb*1.6726e-24*v_nub(i)/ergtkev
         v_nub(i)=sqrt(v_nub(i))
      enddo

c     Close the file.
      close(iunit)

c     Find maximum and minimum particle energies, and use them to
c     calculate if full,half, or 1/3 energy cmpt (ie_nub=1,2,3, resp.)
c     Relative beam deposition rate at each component propto total
c     birth points at each energy.
c     [It is assumed that the list contains only particles born
c      within the plasma.  Add coding to check this.]
      enmin=minval(en_nub,nbirth_pts)
      enmax=maxval(en_nub,nbirth_pts)
      end2=0.5*enmax
      nbirth_cmpts_nub=0d0
      en_avg_cmpts_nub=0d0

      do i=1,nbirth_pts
         if (en_nub(i) .gt. 1.5*end2) then
            ie_nub(i)=1
            nbirth_cmpts_nub(1)=nbirth_cmpts_nub(1)+1
            en_avg_cmpts_nub(1)=en_avg_cmpts_nub(1)+en_nub(i)
         elseif (en_nub(i) .lt. 0.833*end2) then
            ie_nub(i)=3
            nbirth_cmpts_nub(3)=nbirth_cmpts_nub(3)+1
            en_avg_cmpts_nub(3)=en_avg_cmpts_nub(3)+en_nub(i)
         else
            ie_nub(i)=2
            nbirth_cmpts_nub(2)=nbirth_cmpts_nub(2)+1
            en_avg_cmpts_nub(2)=en_avg_cmpts_nub(2)+en_nub(i)
         endif
      enddo

c     Avg energy of each of deposited beam components.
      do i=1,3
         en_avg=en_avg_cmpts_nub(i)
         en_avg_cmpts_nub(i)=en_avg/nbirth_cmpts_nub(i)
      enddo

c     To obtain the normalization for birth_rate_cmpts_nub(), invert
c     pabs(Watt)=sum(i=1,3)[birth_rate_cmpts_nub(i)*en_avg_cmpts_nub(i)]
c                *bcnst
c     to obtain bcnst. Adjust keV energy to joules:*ergtkev*1.e-7(J/erg)
      p_avg_tot_nub=0d0
      do i=1,3
         p_avg_tot_nub=p_avg_tot_nub + nbirth_cmpts_nub(i)*
     +                 en_avg_cmpts_nub(i)*ergtkev*1.d-7  !Watts
      enddo
      bcnst=pabs_nub/p_avg_tot_nub

c     Normalize for deposited power and deposition rate at each energy
      do i=1,3
         birth_rate_cmpts_nub(i)=bcnst*nbirth_cmpts_nub(i)
         pabs_cmpts_nub(i)=birth_rate_cmpts_nub(i)*en_avg_cmpts_nub(i)*
     +                     ergtkev*1.e-7   !Watts
      enddo

c     Alternative calc of birth_rate_cmpts_nub(i)
      do i=1,3
         birth_rate_cmpts_nub1(i)=(total_source_rate/nbirth_pts_nub)*
     +                            nbirth_cmpts_nub(i)
      enddo

c     printout birth_rates, as check for consistency
      write(*,*)
      write(*,*)'birth rate(particles/sec) for each cmpt based on'
      write(*,*)' cmpt   total power    total birth rate'
      do i=1,3
        write(*,106) i,birth_rate_cmpts_nub(i),birth_rate_cmpts_nub1(i)
      enddo
 106  format(i3,5x,e13.6,5x,e13.6)

      return
      end subroutine read_nubeam_data


      subroutine zone(drutpi,ki,mfm1, mim1, mjm1,dri,dzi,psiax,psi,r,z,
     +     xpos,ypos,zpos,pzone,izone)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
         
c      dimension psi(ki,*),r(*),z(*)
      dimension psi(ki,*),r(*),z(*)


      rpos=sqrt(xpos**2+ypos**2)
 120  i=(rpos-r(1))*dri+1.
      j=(zpos-z(1))*dzi+1.
      if(i.gt.mim1) i=mim1
      if(j.gt.mjm1) j=mjm1

      psix  = min(psi(i,j), psi(i+1,j), psi(i,j+1), psi(i+1,j+1))
      ptest = (psix-psiax)*(drutpi/mfm1)**2
      if(ptest.lt.0.02) go to 124
      area1=(rpos-r(i))*(zpos-z(j))
      area2=(r(i+1)-rpos)*(zpos-z(j))
      area3=(r(i+1)-rpos)*(z(j+1)-zpos)
      area4=(rpos-r(i))*(z(j+1)-zpos)
      pzone=(area3*psi(i,j)+area4*psi(i+1,j)+area1*psi(i+1,j+1)
     1  +area2*psi(i,j+1))*dri*dzi
      go to 126
 124  call pfit(psi(i-1,j-1), r(i-1), z(j-1), rpos, zpos, ki, pzone,
     *  dum, dum)
 126  pzone = max(pzone,psiax)
      izone = sqrt(pzone-psiax)*drutpi + 1.

      return
      end subroutine zone
