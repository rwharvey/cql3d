      subroutine tdtranspn
      implicit integer (i-n), real*8 (a-h,o-z)
      dimension  zmatcont(12), janelt(12)
      include 'param.h'
      include 'comm.h'
      include 'trans.h'
      save ifirst
      data ifirst /1/
c..............................................................
c     Radial transport matrix coefficients are generated in
c     CSR format, to be added to the velocity space coefficients
c     for full 2DVel-1Drho matrix prepartory to solution by iterative 
c     sparse-matrix methods:    soln_method="it3drv".
c     "l" is the index for the radial mesh bins.
c     Presently, the set up for transport handles only one general 
c     species (k) at a time.
c     tdtranspn provides an alternative (new) soln method to the
c     splitting soln implemented in tdtransp.
c     BH071105
c.......................................................................
      call dcopy(iyjx2*ngen*lrors,f(0,0,1,1),1,fvn(0,0,1,1),1)
c..............................................................
c     The distribution from the previous time step resides in 
c     fvn(i,j,k,l).  This is now (071111) defined on the complete
c     velocity mesh (ipacktp=0,  see tdtrmuy.f).
c     Finite difference equations have been constructed for a radial
c     diffusion + advection equation as given in the cql3d manual.
c     At each time step, the radial advection term is chosen to
c     keep the density of the diffusing particles near to a specified
c     target radial profile [see cqlinput_help and tdtravct.f].
c     For the j=1 (v=0) points, only the theta=pi (i=iy) point
c     is explicity transported, as the other theta points are
c     set equal to f(i=iy,j=1, l) in the velocity-space setting
c     of coefficients (in impanvnc0.f). [lbdry()="conserv"]
c     For each theta point (i) on the mesh at the outermost 
c     radial bin (l=lrz), there is a smallest minor radius
c     l_lower(i) with corresponding i mesh point.
c     [Presently, only use ymesh="fixed_y" (not "fixed_mu) is fully
c     implemented, which gives l_lower(i)=1.]
c     Boundary conditions are zero-radial-flux at the inner 
c     edge of the rya(l=l_lower(i)) radial bin, and a specified 
c     Maxwellian at the outer radius rya(l=lrz).
c     Better BCs might be future mod [BH071107]:
c        Boundary conditions presently are zero-radial-flux at the
c        inner edge of the rya(l=l_lower(i)) radial bin, and a 
c        specified Maxwellian at the outer radius rya(l=lrz+1/2)=1.0
c        [The outer bc is slightly different from the splitting method
c        setup, in which the distribution at rya(l=lrz) is taken
c        to be Maxwellian for the radial split (but general for the
c        velocity split).]
c
c     At the outer radius, l=lrz, the trapping fraction is
c     greater than towards the magnetic axis.  Therefore, some of the
c     i-points will transform from trapped particles (for which
c     the distribution is symetric about theta=pi/2) to passing
c     particles.  At these radii, given by l=lpt(i), the trapped
c     particles at l.ge.lpt,  with given i, are connected (bound) to
c     corresponding transitting particle distributions at l=lpt-1 
c     in the co-  and counter-dirns, and there is a special equation for 
c     f(i,j,k,l=lpt(i)).  At other points in i,j,l space, 
c     a simple 3-point difference equation holds connecting distrns
c     at l-1,l,l+1, for each i,j. [However, the passing particles
c     at lpt(i)-1 with i.gt.iyh are connected to trapped particles
c     at lpt(i) with i.lt.iyh.]
c..............................................................
c.......................................................................
c  Copy fvn to frn as in tdtransp, although they are now on same grid
c.......................................................................

      call tdtrvtor(fvn,frn)
      
c..............................................................
c     Generate advection coefficients - then determine the
c     Chang-Cooper weights (radial direction).
c     Presently only set up for 1 general species,
c       except for pinch.eq."disabled" (adv()=0.).
c..............................................................
      if (ngen.gt.1) then
         stop 'tdtranspn: Presently set up for only 1 gen species: STOP'
      endif
      kprofile=kelecm
      if (niong.ne.0) kprofile=kionm(1)
      if (colmodl.eq.0) kprofile=1
      call tdtravct(frn,kprofile,1)  !Obtain d_r advection coeffs
                                     !No-op if pinch="disabled"
      call tdtrwtl              !  Sets all dl(,,,)=0.5
c..............................................................
c     Keep a copy for accuracy check later... [following tdtransp].
c..............................................................
      call dcopy(iyjx2*ngen*(lrors+1),frn,1,frn_2,1)
c.......................................................................
c.......................................................................
c  Following do loop setup from impavnc0 might be expanded for
c  other cases.  Use same eqn order as impavnc0 in full radius
c  calcs, for compatibility of the resulting CSR matrix.
c.......................................................................
c.......................................................................
c.......................................................................
c  ipofi(i,l)  (Iprime, i.e., indexsww) of i; set up this array.
c  That is, ipofi(i,l) gives the indexsww value corresponding to
c  to each i, i=1,iy.  For i within the neg v-par trapped region,
c  it uses the the indexsww for theta at symmetic points about y=pi/2.
c  The ipofi(i,l) variable facilitates calculation of corresponding 
c  (i.e., given i,j)equation numbers displaced by one radial
c  mesh point.  
c  It also simplifies radial differencing for each i up to the
c  flux surface where l=lpt(i) [if such condition
c  exists].
c  Recall that i counts counter-clockwise in the vpar-vperp plane,
c    i.e., in the direction of increasing pitch angle theta.
c    ipofi(1,l) will be = inew_(l)-1, ipofi(inew_(l),l)=inew_(l).
c  The indexsww values are as in subroutine impavnc0.
c.......................................................................
      if (ifirst.eq.1) then

         do l=1,lrz
            itl=itl_(l)
            itu=itu_(l)
            iyh=iyh_(l)
            iyy=iy_(l) !-YuP-101215: Don't use iy=; it's in common /params/
                       ! Don't let overwrite the cqlinput value!
            if (itl.lt.iyh ) then
               inew=iyh + itl - 1
               ipassbnd=iyh-itl+1
            else                !i.e., itl=iyh
               inew=iyy
               ipassbnd=0       !i.e., no trapped region for iyh=itl. 
                                !      however, this case not checked out.
            endif
            
            do i = 1,iymax
            ipofi(i,l)=0
            enddo
            ihalf=0
            icntsww=1
            if (iyh.eq.itl) icntsww=0
            do  indexsww=1,inew
               if(indexsww.le.ipassbnd) then
                  i=iyh+1-indexsww
                  ipofi(i,l)=indexsww
                  ipofi(iyy-i+1,l)=indexsww !connects counter-i to co-i
                                           !in the trapped region
               elseif(ihalf.eq.0) then
                  i = itl - icntsww
                  ihalf = 1
                  ipofi(i,l)=indexsww
               else
                  i = itu + icntsww
                  ihalf = 0
                  icntsww = icntsww + 1
                  ipofi(i,l)=indexsww
               endif
            enddo               ! on indexsww
         enddo                  ! on l
         
         ifirst=0
      endif                     ! on ifirst
c.......................................................................
c  Main loops setting up radial transport coefficients
c.......................................................................
      k=1  ! only one gen species for time being
      ieq=0     !equation counter
      icoeff=0  !CSR coefficient counter

      do l=1,lrz   !loop over flux surfaces
         
c         write(*,*)'tdtranspn: l,ieq+1,ieq_(l),icoeff',
c     +                         l,ieq+1,ieq_(l),icoeff
         itl=itl_(l)
         itu=itu_(l)
         iyh=iyh_(l)
         iyy=iy_(l)
c.......................................................................
c     define effective number of unknowns in the theta direction and
c     resulting band widths
c     see below for comments on inew, etc.
c.......................................................................
cBH071029:  During work on 3d fully-implicit solve, found a problem
cBH071029:  with iyh=itl cases.  Fixing this case, and moving ipassbnd
cBH071029:  setting to here [See further notes below for ipassbnd]:
cBH071029:  [But subsequently found had problems with iyh=itl elsewhere
cBH071029:  in the code, so suggest bigger rya(1) or iy for time being.]

      if (itl.lt.iyh ) then
         inew=iyh + itl - 1
         ipassbnd=iyh-itl+1
      else  !i.e., itl=iyh
         inew=iyy
         ipassbnd=0  !i.e., no trapped region for iyh=itl. 
                     !      however, this case not checked out.
      endif

      do j=1,jx   
            ihalf=0
            icntsww=1
            if (iyh.eq.itl) icntsww=0
c..................................................................
c     indexsww is the major loop over I (theta).
c..................................................................
            do indexsww=1,inew
              if(indexsww.le.ipassbnd) then
                i=iyh+1-indexsww
              elseif(ihalf.eq.0) then
                i = itl - icntsww
                ihalf = 1
               else
                i = itu + icntsww
                ihalf = 0
                icntsww = icntsww + 1
              endif

              ieq=ieq+1
              iar_csr(ieq)=icoeff+1   !Coeff counter at beginning of row
              
c             index for scaling of coeffs (saved in impavnc0): 
              if (soln_method.eq."itsol" .or. soln_method.eq."itsol1"
     +           .or. soln_method.eq."direct") then
                 kku=ieq+(k-1)*iyjx  !-YuP-> version from impavnc
              else 
                 kku=ieq-(ieq_(l)-1)+(k-1)*iyjx 
                 ! k.gt.1 not enabled for soln_method it3dv or it3drv
              endif

              if (l.eq.1) then

                 if (j.eq.1) then
                    if (ieq.lt.inew .or. lbdry(k).ne."conserv") then  !0. diagonal element for 1st 
                                           !inew-1 equations
                       icoeff=icoeff+1
                       ar_csr(icoeff)=0.  
                       jar_csr(icoeff)=ieq
                    else        !i.e., ieq.eq.inew (l=1,i=iyy),
                                !link to 2nd flux surface, i-iyy
                       icoeff=icoeff+1
                       ar_csr(icoeff)=scal(kku,l)*betrp(i,l)
                       jar_csr(icoeff)=ieq
                       icoeff=icoeff+1
                       ar_csr(icoeff)=-scal(kku,l)*alprp(i,l)
                       jar_csr(icoeff)=(ieq_(2)-1)+ipofi(i,2)
                   endif
                 else  ! i.e., j.gt.1
c                   No need to distinguish passing or trapped for l=1.
                    icoeff=icoeff+1
                    ar_csr(icoeff)=scal(kku,l)*betrp(i,l)
                    jar_csr(icoeff)=ieq
                    icoeff=icoeff+1
                    ar_csr(icoeff)=-scal(kku,l)*alprp(i,l)
                    jar_csr(icoeff)=ieq_(2)-1+(j-1)*inew_(2)+ipofi(i,2)
                 endif


              elseif (l.lt.lrz) then  !i.e., 1<l<lrz
                 if (j.eq.1 .and. ((ieq-(ieq_(l)-1)).lt.inew  
     ~                          .or. lbdry(k).ne."conserv")) then
                       icoeff=icoeff+1
                       ar_csr(icoeff)=0.
                       jar_csr(icoeff)=ieq
                       
                    
                 elseif (l.eq.lpt(i)) then !at inner radius for trapping
                                           !at the given i: couple to passing
                                           !at smaller l.  CHECK OUT MORE!!!  
                                           !-YuP: corrected. Was:  elseif (i.eq.lpt(i))
                    ii=iyy+1-i

                    icoeff=icoeff+1
                    ar_csr(icoeff)=-scal(kku,l)*0.5*gamrp(i,l)
                    jar_csr(icoeff)=ieq_(l-1)-1+
     +                              (j-1)*inew_(l-1)+ipofi(i,l-1)

                    icoeff=icoeff+1
                    ar_csr(icoeff)=-scal(kku,l)*0.5*gamrp(ii,l)
                    jar_csr(icoeff)=ieq_(l-1)-1+
     +                              (j-1)*inew_(l-1)+ipofi(ii,l-1)

                    icoeff=icoeff+1
                    ar_csr(icoeff)=
     +                        scal(kku,l)*0.5*(betrp(i,l)+betrp(ii,l))
                    jar_csr(icoeff)=ieq

                    icoeff=icoeff+1
                    ar_csr(icoeff)=-scal(kku,l)*alprp(i,l)
                    jar_csr(icoeff)=ieq_(l+1)-1+
     +                              (j-1)*inew_(l+1)+ipofi(i,l+1)
                 else  ! i.e., regular points

                    icoeff=icoeff+1
                    ar_csr(icoeff)=-scal(kku,l)*gamrp(i,l)
                    jar_csr(icoeff)=ieq_(l-1)-1+
     +                              (j-1)*inew_(l-1)+ipofi(i,l-1)

                    icoeff=icoeff+1
                    ar_csr(icoeff)=scal(kku,l)*betrp(i,l)
                    jar_csr(icoeff)=ieq

                    icoeff=icoeff+1
                    ar_csr(icoeff)=-scal(kku,l)*alprp(i,l)
                    jar_csr(icoeff)=ieq_(l+1)-1+
     +                              (j-1)*inew_(l+1)+ipofi(i,l+1)

                 endif

              elseif (l.eq.lrz) then  !no transp for now, just vel distn bc.
                                      !Put placeholders in the CSR matrix

                 icoeff=icoeff+1
                 ar_csr(icoeff)=0.0
                 jar_csr(icoeff)=ieq

              endif

c       if(l.eq.10 .and. j.eq.jx/2) then
c       write(*,'(3i9,e13.5)') 
c     ~  i,iar_csr(ieq),jar_csr(icoeff),ar_csr(icoeff)
c       endif

              
           enddo    ! on indexsww
           
        enddo       ! on j
        
      enddo         ! on l

c  Final element for CSR storage
      iar_csr(ieq+1)=icoeff+1
c      pause

      return
      end
