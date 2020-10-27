c
c
      subroutine rdc_multi
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c
      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

c     Pointers for dynamic memory allocation, local usage:
      real*8, dimension(:), pointer :: upar,uprp,rho_a
      real*8, dimension(:,:), pointer :: tmpb,tmpc,tmpe,tmpf
      real*8, dimension(:,:,:), pointer :: rdc_cqlb,rdc_cqlc,
     +                                     rdc_cqle,rdc_cqlf,
     +                                     ddd
      real*8, dimension(:), allocatable:: rin_med

      real*8, dimension(:), allocatable:: dd_in,dd_out  ! for debugger

c-YuP-> to form file name containing i_psi-number, for "format2"
c-YuP-> diffusion coeffs, one file for one radii, e.g., from DC code.
      character*3 i_psi_index

c.................................................................
c.................................................................


      write(*,*) 'Begin rdc_multi'

c.................................................................
c     Allocate rdc diffusion arrays on iy,jx,lrz grid
c     (pointer statements in comm.h so can use variables elsewhere).
c.................................................................
 
      allocate (rdcb(iy,jx,lrz,nrdc),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("rdcb, sub rdc_multi",0,istat)
      call bcast(rdcb,zero,SIZE(rdcb))
      
      allocate (rdcc(iy,jx,lrz,nrdc),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("rdcc, sub rdc_multi",0,istat)
      call bcast(rdcc,zero,SIZE(rdcc))
      
      allocate (rdce(iy,jx,lrz,nrdc),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("rdce, sub rdc_multi",0,istat)
      call bcast(rdce,zero,SIZE(rdce))
      
      allocate (rdcf(iy,jx,lrz,nrdc),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("rdcf, sub rdf_multi",0,istat)
      call bcast(rdcf,zero,SIZE(rdcf))
      
c.................................................................
c     Read rdc diffusion arrays on u_par,u_prp,lrz grid:
c     first read dimensions, then allocate space, then read.
c.................................................................

 3310 format(1p6e18.10)
  309 format(10i10)

c=======================================================================

      do krf=1,nrdc  !Down to line 1020

      IF(rdcmod.eq.'format1'  .or.  rdcmod.eq.'aorsa') THEN


CMPIINSERT_IF_RANK_EQ_0

      iunit=14
      open(unit=iunit,file=rdcfile(krf),status='old',iostat=kode)
      if (kode.ne.0) then
         WRITE(*,*)
     +        'File rdcfile(krf) cannot be opened, rdcmod.ne.disabled'
         STOP
      endif

      read (iunit, 309) n_uprp   !number of uprp points
      read (iunit, 309) n_upar   !number of upar points
      read (iunit, 309) n_psi    !number of radial points
      
      read (iunit, 3310) vc_cgs  !Max vel (or momentum-rest-mass,
                                 ! if relativistic) on the grid
      read (iunit, 3310) upar_min,upar_max  !Generally, -1. +1.
                                            !Max uprp is assumed
                                            ! to = max upar.
                                            ! min uprp=0.
CMPIINSERT_ENDIF_RANK
CMPIINSERT_BARRIER
CMPIINSERT_BCAST_RDC_GRID

      uprp_min=0.d0
      uprp_max=upar_max

      write(*,*)'rdc_multi: n_uprp,n_upar,n_psi',
     &                      n_uprp,n_upar,n_psi
      write(*,*)'rdc_multi: upar_min,upar_max,uprp_min,uprp_max',
     &                      upar_min,upar_max,uprp_min,uprp_max
      write(*,*)'rdc_multi: vc_cgs=',vc_cgs

c     Check enorm for rf full-wave data:
      
      vc_cgs2=vc_cgs*vc_cgs
      gammac=sqrt(1.d0+vc_cgs2/clite2)
      if ((gammac-1.d0).lt.1.e-6) then
         enormc=0.5d0*fmass(1)*(vc_cgs)**2/ergtkev
      else
         enormc=(gammac-1.d0)*fmass(1)*clite2/ergtkev
      endif
      write(*,*)'Enorm in rdc_multi, gammac =',enormc,gammac


c     Allocate local pointered memory:

      allocate (uprp(n_uprp),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("uprp, sub rdc_multi",0,istat)
      call bcast(uprp,zero,SIZE(uprp))
      
      allocate (upar(n_upar),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("upar, sub rdc_multi",0,istat)
      call bcast(upar,zero,SIZE(upar))
      
      allocate (rho_a(n_psi),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("rho_a, sub rdc_multi",0,istat)
      call bcast(rho_a,zero,SIZE(rho_a))
      
      allocate (rdc_cqlb(n_uprp, n_upar, n_psi),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error(",rdc_cqlb sub rdc_multi",0,istat)
      call bcast(rdc_cqlb,zero,SIZE(rdc_cqlb))
      
      allocate (rdc_cqlc(n_uprp, n_upar, n_psi),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error(",rdc_cqlc sub rdc_multi",0,istat)
      call bcast(rdc_cqlc,zero,SIZE(rdc_cqlc))
      
      allocate (rdc_cqle(n_uprp, n_upar, n_psi),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error(",rdc_cqle sub rdc_multi",0,istat)
      call bcast(rdc_cqle,zero,SIZE(rdc_cqle))
      
      allocate (rdc_cqlf(n_uprp, n_upar, n_psi),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error(",rdc_cqlf sub rdc_multi",0,istat)
      call bcast(rdc_cqlf,zero,SIZE(rdc_cqlf))
      
      
      allocate (tmpb(n_uprp, n_upar),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("tmpb, sub rdc_multi",0,istat)
      call bcast(tmpb,zero,SIZE(tmpb))
      
      allocate (tmpc(n_uprp, n_upar),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("tmpc, sub rdc_multi",0,istat)
      call bcast(tmpc,zero,SIZE(tmpc))
      
      allocate (tmpe(n_uprp, n_upar),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("tmpe, sub rdc_multi",0,istat)
      call bcast(tmpe,zero,SIZE(tmpe))
      
      allocate (tmpf(n_uprp, n_upar),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("tmpf, sub rdc_multi",0,istat)
      call bcast(tmpf,zero,SIZE(tmpf))

c     Read in the grids and diffusion coeffs:

CMPIINSERT_IF_RANK_EQ_0
c     Normalized radius rho_a, uprp, upar (normalized to input vc_cgs)
      read (iunit, 3310) (rho_a(i_psi), i_psi = 1, n_psi)
      read (iunit, 3310) (uprp(i_uprp), i_uprp = 1, n_uprp)
      read (iunit, 3310) (upar(i_upar), i_upar = 1,n_upar)
      
      read (iunit, 3310) (((rdc_cqlb(i_uprp, i_upar, i_psi),
     &     i_uprp = 1, n_uprp), i_upar = 1, n_upar),
     &     i_psi = 1, n_psi)
 
      read (iunit, 3310) (((rdc_cqlc(i_uprp, i_upar, i_psi),
     &     i_uprp = 1, n_uprp), i_upar = 1, n_upar),
     &     i_psi = 1, n_psi)
      
      read (iunit, 3310) (((rdc_cqle(i_uprp, i_upar, i_psi),
     &     i_uprp = 1, n_uprp), i_upar = 1, n_upar),
     &     i_psi = 1, n_psi)

      read (iunit, 3310) (((rdc_cqlf(i_uprp, i_upar, i_psi),
     &     i_uprp = 1, n_uprp), i_upar = 1, n_upar),
     &     i_psi = 1, n_psi)

      WRITE(*,*)'rdc_multi w format1: n_psi,rho_a',
     +          n_psi,(rho_a(i_psi),i_psi=1,n_psi)
     
      close(iunit)
CMPIINSERT_ENDIF_RANK
CMPIINSERT_BARRIER
CMPIINSERT_BCAST_RDC

c.......................................................................
      ELSE !(rdcmod.eq.'format2') ! READ diff.coeffs from separate files
                                  ! Coding only set up for nrdc=1.
c.......................................................................      
            
      write(i_psi_index(1:3),'(3I1)') 0,0,0  ! prepare i_psi_index='000'

c-1-> READ general parameters for grid, min/max, etc. ------------------

CMPIINSERT_IF_RANK_EQ_0
      iunit=14
      open(unit=iunit,file='du0u0_grid',status='old',iostat=kode)
      if (kode.ne.0) then
         WRITE(*,*)
     +     'File du0u0_grid cannot be opened, rdcmod.ne.disabled'
         STOP
      endif

      read (iunit, 309) n_uprp   !number of uprp points
      read (iunit, 309) n_upar   !number of upar points
      read (iunit, 309) n_psi    !number of radial points
      
      read (iunit, 3310) vc_cgs  !Max vel on the grid
      read (iunit, 3310) upar_min,upar_max  !Generally, -1. +1.
                                            !Max uprp is assumed
                                            ! to = max upar.
                                            ! min uprp=0.
      close(iunit)
CMPIINSERT_ENDIF_RANK
CMPIINSERT_BARRIER
CMPIINSERT_BCAST_RDC_GRID

      uprp_min=0.d0
      uprp_max=upar_max

c     Check enorm for rf full-wave data:
      
      vc_cgs2=vc_cgs*vc_cgs
      gammac=sqrt(1.d0+vc_cgs2/clite2)
      if ((gammac-1.d0).lt.1.e-6) then
         enormc=0.5d0*fmass(1)*(vc_cgs)**2/ergtkev
      else
         enormc=(gammac-1.d0)*fmass(1)*clite2/ergtkev
      endif
      write(*,*)'Enorm in rdc_multi, gammac =',enormc,gammac

c-2-> ALLOCATE ARRAYS --------------------------------------------------

      allocate (uprp(n_uprp),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("uprp, sub rdc_multi",0,istat)
      call bcast(uprp,zero,SIZE(uprp))
      
      allocate (upar(n_upar),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("upar, sub rdc_multi",0,istat)
      call bcast(upar,zero,SIZE(upar))
      
      allocate (rho_a(n_psi),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("rho_a, sub rdc_multi",0,istat)
      call bcast(rho_a,zero,SIZE(rho_a))
      
      allocate (rdc_cqlb(n_uprp, n_upar, n_psi),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error(",rdc_cqlb sub rdc_multi",0,istat)
      call bcast(rdc_cqlb,zero,SIZE(rdc_cqlb))
      
      allocate (rdc_cqlc(n_uprp, n_upar, n_psi),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error(",rdc_cqlc sub rdc_multi",0,istat)
      call bcast(rdc_cqlc,zero,SIZE(rdc_cqlc))
      
      allocate (rdc_cqle(n_uprp, n_upar, n_psi),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error(",rdc_cqle sub rdc_multi",0,istat)
      call bcast(rdc_cqle,zero,SIZE(rdc_cqle))
      
      allocate (rdc_cqlf(n_uprp, n_upar, n_psi),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error(",rdc_cqlf sub rdc_multi",0,istat)
      call bcast(rdc_cqlf,zero,SIZE(rdc_cqlf))
      
      allocate (tmpb(n_uprp, n_upar),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("tmpb, sub rdc_multi",0,istat)
      call bcast(tmpb,zero,SIZE(tmpb))
      
      allocate (tmpc(n_uprp, n_upar),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("tmpc, sub rdc_multi",0,istat)
      call bcast(tmpc,zero,SIZE(tmpc))
      
      allocate (tmpe(n_uprp, n_upar),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("tmpe, sub rdc_multi",0,istat)
      call bcast(tmpe,zero,SIZE(tmpe))
      
      allocate (tmpf(n_uprp, n_upar),STAT = istat)
      if(istat .ne. 0)
     .     call allocate_error("tmpf, sub rdc_multi",0,istat)
      call bcast(tmpf,zero,SIZE(tmpf))

CMPIINSERT_IF_RANK_EQ_0
      iunit=14
c-3-> READ data for diff.coeffs from separate files --------------------
      DO i_psi = 1,n_psi   ! Loop through i_psi, read files 'du0u0_r###'
      
        !-> To form file name with i_psi number in it:
        IF(i_psi.LE.9) write(i_psi_index(3:3),'(I1)') i_psi
        IF(i_psi.GE.10  .AND. i_psi.LE.99 ) 
     +                 write(i_psi_index(2:3),'(I2)') i_psi
        IF(i_psi.GE.100 .AND. i_psi.LE.999) 
     +                 write(i_psi_index(1:3),'(I3)') i_psi

        open(unit=iunit,file='du0u0_r'//i_psi_index,status='old',
     +       iostat=kode)
        if (kode.ne.0) then
           WRITE(*,*)
     +       'File du0u0_r### for i_psi=',i_psi,' cannot be opened'
           STOP
        endif
      
c       Normalized radius rho_a, uprp, upar (normalized to input vc_cgs)
        read (iunit, 3310) rho_a(i_psi)  !  one number 
        read (iunit, 3310) (uprp(i_uprp), i_uprp = 1,n_uprp)
        read (iunit, 3310) (upar(i_upar), i_upar = 1,n_upar)
      
        read (iunit, 3310) ((rdc_cqlb(i_uprp, i_upar, i_psi),
     &              i_uprp=1,n_uprp), i_upar=1,n_upar)
 
        read (iunit, 3310) ((rdc_cqlc(i_uprp, i_upar, i_psi),
     &              i_uprp=1,n_uprp), i_upar=1,n_upar)
      
        read (iunit, 3310) ((rdc_cqle(i_uprp, i_upar, i_psi),
     &              i_uprp=1,n_uprp), i_upar=1,n_upar)

        read (iunit, 3310) ((rdc_cqlf(i_uprp, i_upar, i_psi),
     &              i_uprp=1,n_uprp), i_upar=1,n_upar)
     
      ENDDO ! loop in i_psi

      WRITE(*,*)'rdc_multi w format2: n_psi,rho_a',
     +          n_psi,(rho_a(i_psi),i_psi=1,n_psi)
      close(iunit)
      
CMPIINSERT_ENDIF_RANK
CMPIINSERT_BARRIER
CMPIINSERT_BCAST_RDC
c-----------------------------------------------------------------------

      ENDIF ! rdcmod
c=======================================================================



      write(*,*)'rdc_multi: after read of rdcfile(krf), krf=',krf



c.................................................................
c     The rho_a radial mesh and associated coeffc will be reduced
c     by a factor of 2, if lrz.eq.n_psi/2, enabling cql3d to run
c     on half the number of flux surfaces used in the full-wave code.
c     Check if n_psi.eq.lrz.   
c     If not, check n_psi/2.eq.lrz.  If so, omit every 2nd radial point
c       of the diffusion coeff grid and coeffs, and reset n_psi.
c       This enables factor of 2 reduction of the cql3d radial mesh.
c       The cql3d and du0u0_input radial meshes are assumed to be
c         the same (or close).
c       Future modification:  Interpolate the du0u0_input radial mesh
c                             to the cql3d radial mesh.
c     Check if n_psi.eq.lrz:  if not, STOP
     
      if (n_psi.ne.lrz .and.
     +    int((n_psi+1.1)/2) .eq. lrz) then

c.................................................................
c     Adjust the radial mesh to use every second one ==> 32 radial
c     points from 64, or 64 radial points from 128.  
c.................................................................

         n_psi=n_psi/2
         write(*,*)'rdc_multi: n_psi,lrz=',n_psi,lrz
         do i_psi=1,n_psi
            rho_a(i_psi)=rho_a(2*i_psi)
            do i_upar=1,n_upar
               do i_uprp=1,n_uprp
                  rdc_cqlb(i_uprp,i_upar,i_psi)=
     &                 rdc_cqlb(i_uprp,i_upar,2*i_psi)
                  rdc_cqlc(i_uprp,i_upar,i_psi)=
     &                 rdc_cqlc(i_uprp,i_upar,2*i_psi)
                  rdc_cqle(i_uprp,i_upar,i_psi)=
     &                 rdc_cqle(i_uprp,i_upar,2*i_psi)
                  rdc_cqlf(i_uprp,i_upar,i_psi)=
     &                 rdc_cqlf(i_uprp,i_upar,2*i_psi)
               enddo
            enddo
         enddo
         
         write(*,*)'rdc_multi: rho_a, no. elements/2:',
     +        (rho_a(i),i=1,n_psi)
         
      endif

      if (n_psi.ne.lrz) then
         WRITE(*,*)'rdc_multi:  n_psi.ne.lrz'
         STOP
      endif

c.................................................................


c.................................................................
c     Printing max diffusion coefficients, as check
c.................................................................

cTMP
cBH090424      vc_cgs2=vc_cgs*vc_cgs
cBH090424      vc_cgs3=vc_cgs2*vc_cgs
cBH090424      vc_cgs4=vc_cgs2*vc_cgs2

      rdc_cqlbmax=0.d0 
      rdc_cqlfmax=0.d0 
      do i_psi=1,n_psi
      do i_upar=1,n_upar
         do i_uprp=1,n_uprp
            rdc_cqlbmax=max(rdc_cqlbmax,rdc_cqlb(i_uprp,i_upar,i_psi))
cTMP     +                  /vnorm4
cBH090424     +                  /vc_cgs4
     +                  /vnorm4
            rdc_cqlfmax=max(rdc_cqlfmax,rdc_cqlf(i_uprp,i_upar,i_psi))
cTMP     +                  /vnorm2
cBH090424     +                  /vc_cgs2
     +                  /vnorm2
         enddo
      enddo
      enddo
      write(*,*)'rdc_multi: rdc_cqlbmax =',rdc_cqlbmax
      write(*,*)'rdc_multi: rdc_cqlfmax =',rdc_cqlfmax



c  Put diffusion coeffs in code units


      do i_psi=1,n_psi
         do i_upar=1,n_upar
            do i_uprp=1,n_uprp
               rdc_cqlb(i_uprp,i_upar,i_psi)=
cTMP     &              rdc_cqlb(i_uprp,i_upar,i_psi)/vnorm4
cBH090424     &              rdc_cqlb(i_uprp,i_upar,i_psi)/vc_cgs4
     &              rdc_cqlb(i_uprp,i_upar,i_psi)/vnorm4
               rdc_cqlc(i_uprp,i_upar,i_psi)=
cTMP     &              rdc_cqlc(i_uprp,i_upar,i_psi)/vnorm3
cBH090424     &              rdc_cqlc(i_uprp,i_upar,i_psi)/vc_cgs3
     &              rdc_cqlc(i_uprp,i_upar,i_psi)/vnorm3
               rdc_cqle(i_uprp,i_upar,i_psi)=
cTMP     &              rdc_cqle(i_uprp,i_upar,i_psi)/vnorm3
cBH090424     &              rdc_cqle(i_uprp,i_upar,i_psi)/vc_cgs3
     &              rdc_cqle(i_uprp,i_upar,i_psi)/vnorm3
               rdc_cqlf(i_uprp,i_upar,i_psi)=
cTMP     &              rdc_cqlf(i_uprp,i_upar,i_psi)/vnorm2
cBH090424     &              rdc_cqlf(i_uprp,i_upar,i_psi)/vc_cgs2
     &              rdc_cqlf(i_uprp,i_upar,i_psi)/vnorm2
            enddo
         enddo
      enddo

c......................................................................
cBH110430:  For rdc_clipping='enabled',
cBH110430:  clip spikes in the diffusion coeffs which sometimes
cBH110430:  occur at the trapped passing boundary.
cBH110430:  Clipping is according to "Running Median Filters and 
cBH110430:  a General Despiker", John R. Evans, Bull. Seismological
cBH110430:  Soc. Amer. (1982).
cBH110430:  If a peak exceeds trim1*median(2*med_size+1 neighbors),
cBH110430:  trim it.
c......................................................................

      if (rdc_clipping .eq. 'enabled') then

         trim1=5.     ! Trimming factor Dij. Recommended: 4-5
         med_size=2  ! number of diff coeff values on either
                     ! side of the median value.
                     ! Recommend 2 for Duu,Dup,Dpu;  4 for Dpp
         med_size=4 ! For Dpp from DC (i.e., cqlf)

         msp1=med_size+1
         med_window=2*med_size+1
         allocate(rin_med(2*med_size4+1),STAT=istat)

         !Allocate/use ddd for coding convenience
         allocate (ddd(n_uprp, n_upar, 4),STAT = istat)
         if(istat .ne. 0)
     .        call allocate_error("ddd, sub rdc_multi",0,istat)
         call bcast(ddd,zero,SIZE(ddd))

!         !Allocate for debugger testing purposes:
!         allocate (dd_in(n_upar),STAT=istat)
!         allocate (dd_out(n_upar),STAT=istat)

         do i_psi=1,n_psi

            call dcopy(n_upar*n_uprp,rdc_cqlb(1,1,i_psi),1,ddd(1,1,1),1)
            call dcopy(n_upar*n_uprp,rdc_cqlc(1,1,i_psi),1,ddd(1,1,2),1)
            call dcopy(n_upar*n_uprp,rdc_cqle(1,1,i_psi),1,ddd(1,1,3),1)
            call dcopy(n_upar*n_uprp,rdc_cqlf(1,1,i_psi),1,ddd(1,1,4),1)

         do kk=1,3
         do j=1,n_uprp
            do i=1,n_upar
               ! min and max index of median window
               iimn= max(1,   i-med_size)       !limit at low i
               iimn =min(iimn,n_upar-med_window+1)  !limit at high i
               iimx= min(n_upar,i+med_size)     !limit at high i
               iimx= max(iimx,med_window)         !limit at low i
!               dd_in(i)=ddd(j,i,kk)   !For debugger testing
               !r8median gives median value out_med of rin_med
               rin_med=ddd(j,iimn:iimx,kk) ! r8median rearranges input
               out_med=r8median(msp1,med_window,rin_med)
               if (abs(ddd(j,i,kk)).gt.abs(trim1*out_med)) 
     +              ddd(j,i,kk)=out_med
!               dd_out(i)=ddd(j,i,kk)   !For debugger testing
            enddo
         enddo
         enddo

!        Added special treatment for Dpp (i.e., cqlf)
         kk=4
         msp1=med_size4+1
         med_window=2*med_size4+1
         do j=1,n_uprp
            do i=1,n_upar
               ! min and max index of median window
               iimn= max(1,   i-med_size4)       !limit at low i
               iimn =min(iimn,n_upar-med_window+1)  !limit at high i
               iimx= min(n_upar,i+med_size4)     !limit at high i
               iimx= max(iimx,med_window)         !limit at low i
!               dd_in(i)=ddd(j,i,kk)   !For debugger testing
               !r8median gives median value out_med of rin_med
               rin_med=ddd(j,iimn:iimx,kk) ! r8median rearranges input
               out_med=r8median(msp1,med_window,rin_med)
               if (abs(ddd(j,i,kk)).gt.abs(trim1*out_med)) 
     +              ddd(j,i,kk)=out_med
!               dd_out(i)=ddd(j,i,kk)   !For debugger testing
            enddo
         enddo


            call dcopy(n_upar*n_uprp,ddd(1,1,1),1,rdc_cqlb(1,1,i_psi),1)
            call dcopy(n_upar*n_uprp,ddd(1,1,2),1,rdc_cqlc(1,1,i_psi),1)
            call dcopy(n_upar*n_uprp,ddd(1,1,3),1,rdc_cqle(1,1,i_psi),1)
            call dcopy(n_upar*n_uprp,ddd(1,1,4),1,rdc_cqlf(1,1,i_psi),1)

         enddo  ! On i_psi

         deallocate(ddd,STAT=istat)
         deallocate(rin_med,STAT=istat)

      endif  ! On rdc_clipping



c  Invert the order of upar coord, if rdc_upar_sign.lt.0.
c  This is required for DC originated coeffs, since upar in DC is 
c  pos in dirn of B-field,  but in cql3d, upar is positive when
c  in toroidal dirn.
c  For AORSA diffusion coeffs, evidently upar pos is in accord with
c  with cql3d, so use rdc_upar_sign=+1., regardless of btor in eqdsk.
c  We assume here that upar grid is symmetric about parallel vel 0.

      if (rdc_upar_sign.lt.0.d0) then
         do i_psi=1,n_psi
            call dcopy(n_upar*n_uprp,rdc_cqlb(1,1,i_psi),1,tmpb,1)
            call dcopy(n_upar*n_uprp,rdc_cqlc(1,1,i_psi),1,tmpc,1)
            call dcopy(n_upar*n_uprp,rdc_cqle(1,1,i_psi),1,tmpe,1)
            call dcopy(n_upar*n_uprp,rdc_cqlf(1,1,i_psi),1,tmpf,1)
            do i_upar=1,n_upar
               i_uparinv=n_upar-(i_upar-1)
               do i_uprp=1,n_uprp
                  rdc_cqlb(i_uprp,i_uparinv,i_psi)=tmpb(i_uprp,i_upar)
                  rdc_cqlc(i_uprp,i_uparinv,i_psi)=-tmpc(i_uprp,i_upar)
                  rdc_cqle(i_uprp,i_uparinv,i_psi)=-tmpe(i_uprp,i_upar)
                  rdc_cqlf(i_uprp,i_uparinv,i_psi)=tmpf(i_uprp,i_upar)
               enddo
            enddo
         enddo
      endif  ! On rdc_upar_sign



      rdc_cqlbmax=0.d0 
      rdc_cqlfmax=0.d0 
      do i_psi=1,n_psi
      do i_upar=1,n_upar
         do i_uprp=1,n_uprp
            rdc_cqlbmax=max(rdc_cqlbmax,rdc_cqlb(i_uprp,i_upar,i_psi))
            rdc_cqlfmax=max(rdc_cqlfmax,rdc_cqlf(i_uprp,i_upar,i_psi))
         enddo
      enddo
      enddo
      write(*,*)'rdc_multi: rdc_cqlbmax.1 =',rdc_cqlbmax
      write(*,*)'rdc_multi: rdc_cqlfmax.1 =',rdc_cqlfmax

         
c  The velocity grids are assumed to be equi-spaced.  Renormalize
c  to account for possible vc_cgs.ne.vnorm

      xpar_max=upar_max*vc_cgs/vnorm
      xpar_min=upar_min*vc_cgs/vnorm
      xprp_max=uprp_max*vc_cgs/vnorm
      xprp_min=uprp_min*vc_cgs/vnorm
      

c  Setup arrays of initial upar,uprp, normalized momentum-per-rest-mass
c   at the minimum B-field point.  Normalization is adjusted to vnorm.
c   These are du0u0_input arrays.

cTMP      dupar=(upar_max-upar_min)/(n_upar-1)*vc_cgs/vnorm
cBH090424      dupar=(upar_max-upar_min)/(n_upar-1)*vc_cgs/vc_cgs
      dxpar=(xpar_max-xpar_min)/(n_upar-1)
      do i=1,n_upar  !-YuP: was 1,n_upar-1
         upar(i)=xpar_min+(i-1)*dxpar
         !write(*,*)'upar(i),i=1,n_upar',i,upar(i)
      enddo
      
cTMP      duprp=(uprp_max-uprp_min)/(n_uprp-1)*vc_cgs/vnorm
cBH090424      duprp=(uprp_max-uprp_min)/(n_uprp-1)*vc_cgs/vc_cgs
      dxprp=(xprp_max-xprp_min)/(n_uprp-1)
      do i=1,n_uprp !-YuP: was 1,n_uprp-1
         uprp(i)=xprp_min+(i-1)*dxprp
         !write(*,*)'uprp(i),i=1,n_uprp',i,uprp(i)
      enddo
      dxpdxp=dxpar*dxprp

c      write(*,*)'dupar,duprp,dupdup=',dupar,duprp,dupdup


c.................................................................
c     The radial mesh of data from the full-wave code is assumed
c     to be the same as used in cql3d [it may be (slightly) different
c     for some technical reasons].
c     The momemtum-per-mass (velocity) dependent data is interpolated
c     on to the cql3d grids.
c     Bilinear Interpolate rdcb,... to bqlm on x,theta mesh:
c     Divide through by factor symm=2., to account for
c     AORSA calc with bounce time from banana tip to banana tip,
c     whereas present up-down symmetric cql3d calc uses
c     bounce time from equatorial plane to banana tip.
c     (Bilinear interpolation may be preferable to other possibilities,
c     as it does not introduce spurious ringing oscillations;
c     continuity of derivatives does now appear to be a strong
c     consideration for diffusion coeffs. BH)
c.................................................................
c     


      !if (eqsym.eq."none") then !YuP[03/26/2015] symm is set in aingeom now
      !   symm=one
      !else
      !   symm=two
      !endif

      do ll=1,lrz
      do j=1,jx
         do i=1,iy
            xupar=x(j)*coss(i,ll)
            xuprp=x(j)*sinn(i,ll)
c            write(*,*)'upar_min,upar_max,uprp_min,uprp_max,xupar,xuprp',
c     +                 upar_min,upar_max,uprp_min,uprp_max,xupar,xuprp
            if (xupar.lt.xpar_max .and. xupar.gt.xpar_min) then
               if (xuprp.lt.xprp_max .and. xuprp.ge.xprp_min) then
                  ipar=(xupar-xpar_min)/dxpar+1
                  iprp=(xuprp-xprp_min)/dxprp+1
                  w00=(upar(ipar+1)-xupar)*(uprp(iprp+1)-xuprp)/dxpdxp
c                  write(*,*)'upar(ipar+1),uprp(iprp+1),dupdup',
c     +                       upar(ipar+1),uprp(iprp+1),dupdup
                  w01=(upar(ipar+1)-xupar)*(xuprp-uprp(iprp))/dxpdxp
                  w10=(xupar-upar(ipar))*(uprp(iprp+1)-xuprp)/dxpdxp
                  w11=(xupar-upar(ipar))*(xuprp-uprp(iprp))/dxpdxp
                  wsum=w00+w01+w10+w11
c                  write(*,*)'ipar,iprp,w00,w01,w10,w11,wsum=',
c    +                       ipar,iprp,w00,w01,w10,w11,wsum
                  rdcb(i,j,ll,krf)=(w00*rdc_cqlb(iprp,ipar,ll)
     +                     +w01*rdc_cqlb(iprp+1,ipar,ll)
     +                     +w10*rdc_cqlb(iprp,ipar+1,ll)
     +                     +w11*rdc_cqlb(iprp+1,ipar+1,ll))/symm
                  
                  rdcc(i,j,ll,krf)=(w00*rdc_cqlc(iprp,ipar,ll)
     +                     +w01*rdc_cqlc(iprp+1,ipar,ll)
     +                     +w10*rdc_cqlc(iprp,ipar+1,ll)
     +                     +w11*rdc_cqlc(iprp+1,ipar+1,ll))/symm

                  rdce(i,j,ll,krf)=(w00*rdc_cqle(iprp,ipar,ll)
     +                     +w01*rdc_cqle(iprp+1,ipar,ll)
     +                     +w10*rdc_cqle(iprp,ipar+1,ll)
     +                     +w11*rdc_cqle(iprp+1,ipar+1,ll))/symm

                  rdcf(i,j,ll,krf)=(w00*rdc_cqlf(iprp,ipar,ll)
     +                     +w01*rdc_cqlf(iprp+1,ipar,ll)
     +                     +w10*rdc_cqlf(iprp,ipar+1,ll)
     +                     +w11*rdc_cqlf(iprp+1,ipar+1,ll))/symm

                  
c                  write(*,*)'diff_coeff_rd: i,j,ll,rdcb(i,j,ll)',
c     +                      i,j,ll,rdcb(i,j,ll)
               endif
            endif
         enddo
      enddo
      enddo

c.................................................................
c     Printing max diffusion coefficient, as check
c.................................................................
      rdcbmax=0.d0 
      do ll=1,lrz
      do j=1,jx
         do i=1,iy
            rdcbmax=max(rdcbmax,rdcb(i,j,ll,krf))
         enddo
      enddo
      enddo
      write(*,*)'rdc_multi: krf=',krf,' rdcbmax.1 =',rdcbmax

c.................................................................
c     Scaling diffusion coeffs
c.................................................................

      do ll=1,lrz
         do j=1,jx
            do i=1,iy
               rdcb(i,j,ll,krf)=rdcscale(krf)*rdcb(i,j,ll,krf)
               rdcc(i,j,ll,krf)=rdcscale(krf)*rdcc(i,j,ll,krf)
               rdcf(i,j,ll,krf)=rdcscale(krf)*rdcf(i,j,ll,krf)
               rdce(i,j,ll,krf)=rdcscale(krf)*rdce(i,j,ll,krf)
            enddo
         enddo
      enddo

c.................................................................
c     Examining degree of non-parabolic diffusion coeffs
c     Summing pos and neg cqlb,cqlf over mesh volumes.
c     [Could consider weighting by f/n, but the distns
c      are initially Maxwln, so might not want that.]
c     
c.................................................................


      sum_pos_duu=0.d0
      sum_neg_duu=0.d0
      sum_pos_dtt=0.d0
      sum_neg_dtt=0.d0
      do ll=1,lrz
      do j=1,jx
         do i=1,iy 
            sum_pos_duu= sum_pos_duu+cynt2(i,ll)*cint2(j)*
     +           max(zero,rdcb(i,j,ll,krf))
            sum_neg_duu= sum_neg_duu+cynt2(i,ll)*cint2(j)*
     +           min(zero,rdcb(i,j,ll,krf))
            sum_pos_dtt= sum_pos_dtt+cynt2(i,ll)*cint2(j)*
     +           max(zero,rdcf(i,j,ll,krf))
            sum_neg_dtt= sum_neg_dtt+cynt2(i,ll)*cint2(j)*
     +           min(zero,rdcf(i,j,ll,krf))
         enddo
      enddo
      enddo

      write(*,*)'rdc_multi: krf=',krf,' sum_pos_duu,sum_neg_duu',
     +                      krf,sum_pos_duu,sum_neg_duu  
      write(*,*)'rdc_multi: sum_pos_dtt,sum_neg_dtt',
     +                      krf,sum_pos_dtt,sum_neg_dtt

c$$$c     Non-symmetry of rdcc,rdce  [Visual exam shows
c$$$c       this non-symmetry is general].
c$$$      do ll=1,lrz
c$$$      do j=1,jx
c$$$         do i=1,iy
c$$$         if (rdcc(i,j,ll).ne.rdce(i,j,ll)) then
c$$$         write(*,*)'Non-sym coeffs: ll,i,j,rdcc(i,j,ll),rdce(i,j,ll)',
c$$$     +                              ll,i,j,rdcc(i,j,ll),rdce(i,j,ll)
c$$$         endif
c$$$         enddo
c$$$      enddo
c$$$      enddo


c     Non-parabolic cross-coeffs
c$$$      do ll=1,lrz
c$$$      do j=1,jx
c$$$         do i=1,iy
c$$$         if (rdcc(i,j,ll)*rdce(i,j,ll) .gt.
c$$$     +             1.0005*rdcb(i,j,ll)*rdcf(i,j,ll)) then
c$$$         write(*,*)'Non-para coeffs: ll,i,j,rdcc(i,j,ll)*rdce(i,j,ll),
c$$$     +1.0005*rdcb(i,j,ll)*rdcf(i,j,ll)',
c$$$     +           ll,i,j,rdcc(i,j,ll)*rdce(i,j,ll),
c$$$     +           1.0005*rdcb(i,j,ll)*rdcf(i,j,ll)
c$$$         endif
c$$$                   
c$$$         enddo
c$$$      enddo
c$$$      enddo


c..................................................................
c     Symmetrize about pi/2 in the pass/trapped boundary.
c..................................................................

      do ll=1,lrz
        call tdnflxs(ll)
        do 56 j=1,jx
          do 57 i=itl,itu
            temp1(i,j)=rdcb(i,j,indxlr_,krf)
            temp2(i,j)=rdcc(i,j,indxlr_,krf)
            temp3(i,j)=rdce(i,j,indxlr_,krf)
            temp4(i,j)=rdcf(i,j,indxlr_,krf)
 57       continue
 56     continue
        do 60 j=1,jx
          do 70 i=itl,itu
            rdcb(i,j,indxlr_,krf)=(temp1(i,j)+temp1(iy+1-i,j))*.5
            rdcc(i,j,indxlr_,krf)=(temp2(i,j)-temp2(iy+1-i,j))*.5
            rdce(i,j,indxlr_,krf)=(temp3(i,j)-temp3(iy+1-i,j))*.5
            rdcf(i,j,indxlr_,krf)=(temp4(i,j)+temp4(iy+1-i,j))*.5
 70       continue
 60     continue
      enddo



c.................................................................
c     Modifying coeffs to ensure parabolic diffusion coeffs
c.................................................................
c      goto 100 !!!!!!!!!!!!!!!!!!!!!!!! skip adjustments ------------------------

      do ll=1,lrz
         call tdnflxs(ll)
         do j=1,jx
            do i=1,iy
               ! rdcb= Duu*|upar|*taub*u^2 /vnorm^4  
               ! rdcf= Dtt*sin(theta)*|upar|*taub*u^2 /vnorm^2 
               ! Diagonal diffusion elements should not be negative:
               rdcb(i,j,ll,krf)=max(rdcb(i,j,ll,krf),zero) 
               rdcf(i,j,ll,krf)=max(rdcf(i,j,ll,krf),zero)
               ! rdcc= Dut*|upar|*taub*u^2 /vnorm^3
               ! rdce= Dtu*sin(theta)|upar|*taub*u^2 /vnorm^3
               ! Where one is zero, another must also be zero:
               if (rdcc(i,j,ll,krf).eq.zero .and. 
     +            rdce(i,j,ll,krf).ne.zero) then
                  rdce(i,j,ll,krf)=zero
c$$$                  write(*,*)'rdcc zero but rdce not,i,j,ll',
c$$$     +                                              i,j,ll
               endif
               if (rdce(i,j,ll,krf).eq.zero .and.
     +            rdcc(i,j,ll,krf).ne.zero) then
                  rdcc(i,j,ll,krf)=zero
c$$$                  write(*,*)'rdce zero but rdce not,i,j,ll',
c$$$     +                                              i,j,ll
               endif
               ! Dut=Dtu, so rdcc and rdce should have same sign: 
               if (rdcc(i,j,ll,krf)*rdce(i,j,ll,krf).lt.zero) then
                  !-YuP: Which to reverse? C or E?
                  rdcc(i,j,ll,krf)=-1.*rdcc(i,j,ll,krf) 
                  !-YuP: Another possibility (~ same result) -
                  !-YuP: set both C and E to zero:
                  rdcc(i,j,ll,krf)=zero 
                  rdce(i,j,ll,krf)=zero 
c$$$                  write(*,*)'Neg rdcc(i,j,ll)*rdce(i,j,ll),i,j,ll',
c$$$     +                                                     i,j,ll
               endif
c     Further ensuring have parabolic equations.
               ! Ensure that  rdcb*rdcf-rdcc*rdce =0 !-YuP: =0 or >0 ??
               B0F0= rdcb(i,j,ll,krf)*rdcf(i,j,ll,krf)
               C0E0= rdcc(i,j,ll,krf)*rdce(i,j,ll,krf)
               if (B0F0 .ne. zero) then
                  !if (C0E0.ne.zero) then !-YuP: to enforce C0E0=B0F0
                  if (C0E0.gt.B0F0) then !-YuP: to enforce C0E0 =< B0F0
                     scaled=sqrt(B0F0/C0E0)
                     rdcc(i,j,ll,krf)=rdcc(i,j,ll,krf)*scaled
                     rdce(i,j,ll,krf)=rdce(i,j,ll,krf)*scaled
                  else ! C0E0=0 => Should have rdcb*rdcf=0.
                     !-YuP: Could be either rdcb=0 or rdcf=0
                     !-YuP: How to choose? Usually rdcf is small, 
                     !-YuP: so set it to zero, but leave Duu as it is:
                     !rdcf(i,j,ll)=zero ! Comment this line when using if(C0E0>B0F0)
                  endif
               else ! B0F0=0 => Should have rdcc*rdce=0.
                  ! Dut=Dtu so that both should be zero:
                  rdcc(i,j,ll,krf)=zero
                  rdce(i,j,ll,krf)=zero
               endif
               !rdce(i,j,ll)=zero ! Un-comment to set E0=0:
               !rdcc(i,j,ll)=-rdcc(i,j,ll) ! Un-comment to reverse C0
               !rdce(i,j,ll)=-rdce(i,j,ll) ! Un-comment to reverse E0
            enddo ! i
         enddo ! j
      enddo ! ll



 100  continue



c..................................................................
c     Taper diffusion over last 10 points of velocity, 
c     if ineg="trunc_d"
c..................................................................

      rdcbmax0=0.d0 
      do ll=1,lrz
      do j=1,jx
         do i=1,iy
            rdcbmax0=max(rdcbmax0,rdcb(i,j,ll,krf))
         enddo
      enddo
      enddo
      write(*,*)'rdc_multi: krf=',krf,' rdcbmax.2 =',rdcbmax0

      if (ineg.eq."trunc_d") then
        if (jx.le.11) stop 'rdc_multi:  Need jx>11'
        do ll=1,lrz
        do 90 j=jx-11,jx
          do 91 i=1,iy
            rdcb(i,j,ll,krf)=truncd(j)*rdcb(i,j,ll,krf)
            rdcc(i,j,ll,krf)=truncd(j)*rdcc(i,j,ll,krf)
            rdce(i,j,ll,krf)=truncd(j)*rdce(i,j,ll,krf)
            rdcf(i,j,ll,krf)=truncd(j)*rdcf(i,j,ll,krf)
 91       continue
 90     continue
        enddo
      endif
c.................................................................
c     Printing max diffusion coefficient, as check
c.................................................................
      rdcbmax=0.d0 
      do ll=1,lrz
      do j=1,jx
         do i=1,iy
            rdcbmax=max(rdcbmax,rdcb(i,j,ll,krf))
         enddo
      enddo
      enddo
      write(*,*)'rdc_multi: krf=',krf,' rdcbmax.3 =',rdcbmax

c.................................................................
c     Plotting diffusion coefficient
c.................................................................
      mrfn=1

      do ll=1,lrz
         call tdnflxs(ll)
         if(pltrdc.ne."disabled") then
            if (pltrdc.eq."one" .or. pltrdc.eq."enabled") then
               call rdc_bplt(1)
            else
               call rdc_bplt(krf)
            endif
         endif
      enddo

c.................................................................
c     Writing diffusion coefficient to netcdf file (comment out if
c     don't need this any longer).
c.................................................................
      
      if(rdc_netcdf.ne."disabled") then
         if (rdc_netcdf.eq."one" .or. rdc_netcdf.eq."enabled") then
            call netcdf_rdcb(1)
         else
            call netcdf_rdcb(krf)
         endif
      endif

      write(*,*)' END of rdc_multi'
c      write(*,*)'STOP: END of rdc_multi' 
c      STOP

c.................................................................
c     Dallocate local storage
c.................................................................

      deallocate(upar,uprp,rho_a,STAT=istat1)
      deallocate(rdc_cqlb,rdc_cqlc,rdc_cqle,rdc_cqlf,STAT=istat2)
      deallocate(tmpb,tmpc,tmpe,tmpf,STAT=istat3)
      if (istat1.ne.0 .or. istat2.ne.0 .or. istat3.ne.0) then
         WRITE(*,*)'rdc_multi: dallocation error'
         STOP
      endif

c.................................................................
c     enddo on krf=1,nrdc, line 66
c.................................................................
      
      enddo

      return
      end


      real*8 function r8median(k,n,arr)
      implicit integer (i-n), real*8 (a-h,o-z)
c
c     Returns the kth smallest value in the array arr(1:n).
c     The input array will be rearranged to have this value
c     in location arr(k), with all smaller elements moved to
c     arr(1:k-1) (in arbitrary order) and all larger elements in
c     arr(k+1:n) (also in arbitrary order).
c     Thus, for odd n and k=n/2+1, arr(k) contains the median value.
c     Value arr(k) is returned in r8median.
c
      dimension arr(n)
      
      l=1
      ir=n
      
 1    if(ir-l .le. 1) then
         if(ir-l .eq. 1) then
            if(arr(ir) .lt. arr(l)) then
               temp=arr(l)
               arr(l)=arr(ir)
               arr(ir)=temp
            endif
         endif
         r8median=arr(k)
         return
      else
         mid=(l+ir)/2
         temp=arr(mid)
         arr(mid)=arr(l+1)
         arr(l+1)=temp
         if (arr(l+1) .gt. arr(ir)) then
            temp=arr(l+1)
            arr(l+1)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l) .gt. arr(ir)) then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l+1).gt.arr(l)) then
            temp=arr(l+1)
            arr(l+1)=arr(l)
            arr(l)=temp
         endif
         i=l+1
         j=ir
         a=arr(l)

 3       continue
         i=i+1
         if(arr(i).lt.a) go to 3

 4       continue
         j=j-1
         if(arr(j) .gt. a) go to 4
         if(j .lt. i) go to 5
         temp=arr(i)
         arr(i)=arr(j)
         arr(j)=temp
         go to 3

 5       arr(l)=arr(j)
         arr(j)=a
         if(j .ge. k) ir=j-1
         if(j .le. k) l=i
      endif
      go to 1

      end
      
