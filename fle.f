      !YuP[2019-09-29] In the three fle_*** subroutines below, 
      ! changed f() to f_(). 
      ! The reason is related to how these subroutines 
      ! (particularly fle()) are called in MPI runs. 
      ! Subr. fle() can be called by srckotot, which is called
      ! through subr. sourcee(). The sourcee is called from achiefn()
      ! inside the ll-loop (radial index), which is parallelized
      ! for finding a new solution f() at each ll (l_) in parallel.
      ! However, the call to sourcee is not in MPI, 
      ! it is made at every mpirank, including mpirank=0.
      ! So, mpirank=0 calls subr.sourcee for every l_
      ! (while subr.impavnc0() that finds the new f() is called 
      ! at mpiranks>0, each mpirank working on certain range of l_).
      ! It would not affect anything if the sources were independent
      ! from distr.func f(). But subr.fle does depend on f(),
      ! so it may happen that mpirank=0 receives a new f(...,l_)
      ! from some other core BEFORE it starts calculation 
      ! of sources at surface l_. So it will use the new f(l_) 
      ! instead of old f(l_) as it does in a serial run.
      ! To prevent this, I replaced f() by f_(), which is 
      ! the distr.func. at previous time step (i.e., the old f()). 
c
      subroutine fle_pol(setup,lp) !called by subr.pltprppr for each lr_
      implicit integer (i-n), real*8 (a-h,o-z)
      character*(*) setup
c
c.......................................................................
c  At each poloidal angle, specified by poloidal index lp:
c  Computes a local reduced distribution function,  
c  fl(0:jfl)=f(v_parallel), by binning line-averaged density in
c  dz(lp,ll) due to each  equatorial plane velocity-space element 
c  into a parallel velocity  space grid xl(1:jfl).
c  fl(1:jfl-1) are the relevant (non-zero) entries, each centered
c    at normalized velocity xlm(1:jfl-1).  The xl(1:jfl) give bin bndries.
c    The xl-grid is independent of poloidal angle
c  Presently set up only for single flux surface runs.
c
c  Normalization of the output distribution fl is such that it equals
c  vnorm*f_cgs(v_parallel).
c.......................................................................
      save
      include 'param.h'
      include 'comm.h'

c     Diagnostic array:
c      dimension den_of_s(lza)  !Now comm.h:den_of_s2

      character*8 ifirst
      data ifirst/"first"/

      if (lrz.gt.1) stop 'in fle_pol: only set up for lrz=1.'

c     At the first call, set up the mesh, binning and weighting factors.
c     The mesh is the same for each poloidal angle, but the factors
c     vary with poloidal angle.
      if (setup.eq."setup") then
      if (ifirst.eq."first") then
         ifirst="notfirst"
 
c     initialize density diagnostic
         do l=1,lza
            den_of_s2(l)=0.0
         enddo
         ll=1 ! YuP[2019-09] Why ll=1 ? 
         !This subr.fle_pol is called by subr.pltprppr for each lr_,
         ! so I would use ll=lr_ here (or l_). [YuP]
         ll=l_ !YuP[2020-01] No effect from this change
         denfl2(ll)=0.0

c     Set up the mesh:
         jflh=(jfl+1)/2
      if (xlfac .ge. 0.) then

        xl(1)=-xmax
        xl(jflh)=0.
        xl(jfl)=xmax
        hx=xmax/(jfl-jflh)
        hxfg=xlfac*hx
        xl(jflh+1)=hxfg
        xl(jflh-1)=-hxfg
        call micgetr(jflh,xmax,hxfg,ram,ksingul)
        if (ksingul .eq. 1) call diagwrng(8)
        do 120 j=jflh+2,jfl-1
          xl(j)=xl(j-1)+ram*(xl(j-1)-xl(j-2))
          jj=jfl+1-j
          xl(jj)=-xl(j)
 120    continue
      else
        jflwr=xlpctlwr*jflh
        if(jflwr .lt. 3) jflwr=3
        hx=xlwr/(jflwr-1)
        xl(jflh)=0.
        do 200 j=jflh+1,jflh-1+jflwr
          xl(j)=hx*(j-jflh)
 200    continue
        jmdl1=xlpctmdl*jflh
        jmdl=jmdl1+jflwr
        hm=(xlmdl-xllwr)/jmdl1
        do 201 j=jflh+jflwr,jflh-1+jmdl
          xl(j)=xl(j-1)+hm
 201    continue
        hu=(xmax-xlmdl)/(jflh-jmdl)
        do 202 j=jflh+jmdl,jfl
          xl(j)=xl(j-1)+hu
 202    continue
        do 203 j=jflh+1,jfl
           jj=jfl+1-j
           xl(jj)=-xl(j)
 203    continue
      endif

      do j=1,jfl-1
         xlm(j)=0.5*(xl(j)+xl(j+1))
         dxl(j)=xl(j+1)-xl(j)
      enddo
      xlm(0)=xl(1)-0.5*dxl(1)
      xlm(jfl)=xl(jfl)+0.5*dxl(jfl-1)
      dxl(0)=dxl(1)
      dxl(jfl)=dxl(jfl-1)


c        Next loop is over the FP'd flux surfaces.
c        A point i,j of the equatorial distribution function will 
c          contribute to parallel velocity bin jflbin(i,j,ll) with
c          weight wtfl0(i,j,ll) and to bin jflbin(i,j,ll)-1 with
c          weight wtflm(i,j,ll). 
c         do 204 ll=1,lrz
         ll=1 ! YuP[2019-09] Why ll=1 ? 
         !This subr.fle_pol is called by subr.pltprppr for each lr_,
         ! so I would use ll=lr_ here (or l_). [YuP]
         ll=l_ !YuP[2020-01] No effect from this change

         if (lrz.ne.1) stop 'lrz.ne.1 in fle_pol'

         do 2041 l=1,lz
         do 205 j=1,jx
         do 206 i=1,min(imax(l,ll)+1,iyh_(ll))
            xll=x(j)*cosz(i,l,ll)
            jbin=luf(xll,xlm(1),jfl)
            wt=(xlm(jbin)-xll)/(xlm(jbin)-xlm(jbin-1))
            wta=dtau(i,l,ll)*coss(i,ll)/dz(l,ll)*cynt2(i,ll)*cint2(j)
            jflbin(i,j,l)=jbin
            wtfl0(i,j,l)=(1.-wt)*wta
            wtflm(i,j,l)=wt*wta
 206     continue
 205     continue

         do 207 j=1,jx
         do 208 i=1,min(imax(l,ll)+1,iyh_(ll))
            ii=iy_(ll)+1-i
            jflbin(ii,j,l)=jfl-jflbin(i,j,l)+1
c            wtfl0(ii,j,l)=wtfl0(i,j,l)
c            wtflm(ii,j,l)=wtflm(i,j,l)
            wtflm(ii,j,l)=wtfl0(i,j,l)
            wtfl0(ii,j,l)=wtflm(i,j,l)
 208     continue
 207     continue

 2041    continue ! l=1,lz
c 204     continue

      endif
         go to 999
      endif

c  Form the parallel distribution function, for a given pol. angle lp
c  and radial bin l_:

      call bcast(fl(0),zero,jfl+1)
      do j=1,jx
         do i=1,min(imax(lp,ll)+1,iyh_(ll))
            itemc1(i)=jflbin(i,j,lp)
            itemc2(i)=itemc1(i)-1
c            temc1(i)=wtfl0(i,j,lp)*f_(i,j,k,l_)
c            temc2(i)=wtflm(i,j,lp)*f_(i,j,k,l_)
            temc1(i)=wtfl0(i,j,lp)*f_(i,j,1,l_)
            temc2(i)=wtflm(i,j,lp)*f_(i,j,1,l_)
         enddo
         do i=1,min(imax(lp,ll)+1,iyh_(ll))
            ii=iy_(ll)+1-i
            itemc1(ii)=jflbin(ii,j,lp)
            itemc2(ii)=itemc1(ii)-1
            temc1(ii)=wtfl0(ii,j,lp)*f_(ii,j,1,l_)
            temc2(ii)=wtflm(ii,j,lp)*f_(ii,j,1,l_)
         enddo

         do i=1,min(imax(lp,ll)+1,iyh_(ll))
            fl(itemc1(i))=fl(itemc1(i))+temc1(i)
            fl(itemc2(i))=fl(itemc2(i))+temc2(i)
         enddo
         do i=1,min(imax(lp,ll)+1,iyh_(ll))
            ii=iy_(ll)+1-i
            fl(itemc1(ii))=fl(itemc1(ii))+temc1(ii)
            fl(itemc2(ii))=fl(itemc2(ii))+temc2(ii)
         enddo
      enddo

      do jf=0,jfl
         fl(jf)=bbpsi(lp,lr_)*fl(jf)/dxl(jf)
      enddo
      fl(1)=fl(1)+fl(0)
      fl(jfl-1)=fl(jfl-1)+fl(jfl)
      fl(0)=0.0
      fl(jfl)=0.0

c  Set minimum fl = 1.e-100* (max. value)
      flmax=0.
      do jf=1,jfl-1
c990131         flmax=amax1(flmax,fl(jf))
         flmax=max(flmax,fl(jf))
      enddo
      flmin=em100*flmax
      do jf=1,jfl-1
c990131         fl(jf)=amax1(flmin,fl(jf))
         fl(jf)=max(flmin,fl(jf))
      enddo

c     Diagnostic check of densities:
      do jf=1,jfl-1
         den_of_s2(lp)=den_of_s2(lp)+dxl(jf)*fl(jf)
      enddo
      denfl2(l_)=denfl2(l_)
     1          +dz(lp,ll)/bbpsi(lp,ll)*den_of_s2(lp)/zmaxpsi(ll)

 999  return
      end subroutine fle_pol




c
c
      subroutine fle_fsa(setup)
      implicit integer (i-n), real*8 (a-h,o-z)
      character*(*) setup
c
c.......................................................................
c  Computes a flux-surface averaged reduced distribution function,  
c  fl(0:jfl)=f(v_parallel), by binning line-averaged density due to each
c  equatorial plane velocity-space element into a parallel velocity  
c  space grid xl(1:jfl).
c  fl(1:jfl-1) are the relevant (non-zero) entries, each centered
c    at normalized velocity xlm(1:jfl-1).  The xl(1:jfl) give bin bndries.
c  This calculation is set up for use with multiple flux surface runs.
c
c  Normalization of the output distribution fl is such that it equals
c  vnorm*f_cgs(v_parallel).
c.......................................................................
      save
      include 'param.h'
      include 'comm.h'

      character*8 ifirst
      data ifirst/"first"/

      if (lrz.gt. lz) stop 'in fle_fsa:  Need to adjust lrz.le.lz'

      if (setup.eq."setup") then
c     At the first call, set up the mesh, binning and weighting factors.
      if (ifirst.eq."first") then
         ifirst="notfirst"

c     Set up the mesh:
         jflh=(jfl+1)/2
      if (xlfac .ge. 0.) then

        xl(1)=-xmax
        xl(jflh)=0.
        xl(jfl)=xmax
        hx=xmax/(jfl-jflh)
        hxfg=xlfac*hx
        xl(jflh+1)=hxfg
        xl(jflh-1)=-hxfg
        call micgetr(jflh,xmax,hxfg,ram,ksingul)
        if (ksingul .eq. 1) call diagwrng(8)
        do 120 j=jflh+2,jfl-1
          xl(j)=xl(j-1)+ram*(xl(j-1)-xl(j-2))
          jj=jfl+1-j
          xl(jj)=-xl(j)
 120    continue
      else
        jflwr=xlpctlwr*jflh
        if(jflwr .lt. 3) jflwr=3
        hx=xlwr/(jflwr-1)
        xl(jflh)=0.
        do 200 j=jflh+1,jflh-1+jflwr
          xl(j)=hx*(j-jflh)
 200    continue
        jmdl1=xlpctmdl*jflh
        jmdl=jmdl1+jflwr
        hm=(xlmdl-xllwr)/jmdl1
        do 201 j=jflh+jflwr,jflh-1+jmdl
          xl(j)=xl(j-1)+hm
 201    continue
        hu=(xmax-xlmdl)/(jflh-jmdl)
        do 202 j=jflh+jmdl,jfl
          xl(j)=xl(j-1)+hu
 202    continue
        do 203 j=jflh+1,jfl
           jj=jfl+1-j
           xl(jj)=-xl(j)
 203    continue
      endif

      do j=1,jfl-1
         xlm(j)=0.5*(xl(j)+xl(j+1))
         dxl(j)=xl(j+1)-xl(j)
      enddo
      xlm(0)=xl(1)-0.5*dxl(1)
      xlm(jfl)=xl(jfl)+0.5*dxl(jfl-1)
      dxl(0)=dxl(1)
      dxl(jfl)=dxl(jfl-1)


c        Next loop is over the FP'd flux surfaces.
c        A point i,j of the equatorial distribution function will 
c          contribute to parallel velocity bin jflbin(i,j,ll) with
c          weight wtfl0(i,j,ll) and to bin jflbin(i,j,ll)-1 with
c          weight wtflm(i,j,ll). 
         do 204 ll=1,lrz

         do 205 j=1,jx
         do 206 i=1,iyh_(ll)
            xll=x(j)*coss(i,ll)
            jbin=luf(xll,xlm(1),jfl)
            wt=(xlm(jbin)-xll)/(xlm(jbin)-xlm(jbin-1))
            wta=vptb(i,ll)*cynt2(i,ll)*cint2(j)
            jflbin(i,j,ll)=jbin
            wtfl0(i,j,ll)=(1.-wt)*wta
            wtflm(i,j,ll)=wt*wta
 206     continue
 205     continue

         do 207 j=1,jx
         do 208 i=1,iyh_(ll)
            ii=iy_(ll)+1-i
            jflbin(ii,j,ll)=jfl-jflbin(i,j,ll)+1
c            wtfl0(ii,j,ll)=wtfl0(i,j,ll)
c            wtflm(ii,j,ll)=wtflm(i,j,ll)
            wtflm(ii,j,ll)=wtfl0(i,j,ll)
            wtfl0(ii,j,ll)=wtflm(i,j,ll)
 208     continue
 207     continue


 204     continue

         go to 999

      endif
      endif

c  Form the parallel distribution function, for a given radius l_:

      call bcast(fl(0),zero,jfl+1)
      do j=1,jx
         do i=1,iy_(l_)
            itemc1(i)=jflbin(i,j,l_)
            itemc2(i)=itemc1(i)-1
c            temc1(i)=wtfl0(i,j,l_)*f_(i,j,k,l_)
c            temc2(i)=wtflm(i,j,l_)*f_(i,j,k,l_)
            temc1(i)=wtfl0(i,j,l_)*f_(i,j,1,l_)
            temc2(i)=wtflm(i,j,l_)*f_(i,j,1,l_)
         enddo
         do i=1,iy_(l_)
            fl(itemc1(i))=fl(itemc1(i))+temc1(i)
            fl(itemc2(i))=fl(itemc2(i))+temc2(i)
         enddo
      enddo

      do jf=0,jfl
         fl(jf)=fl(jf)/(dxl(jf)*zmaxpsi(lr_))
      enddo
      fl(1)=fl(1)+fl(0)
      fl(jfl-1)=fl(jfl-1)+fl(jfl)
      fl(0)=0.0
      fl(jfl)=0.0

c  Set minimum fl = 1.e-100* (max. value)
      flmax=0.
      do jf=1,jfl-1
c990131         flmax=amax1(flmax,fl(jf))
         flmax=max(flmax,fl(jf))
      enddo
      flmin=em100*flmax
      do jf=1,jfl-1
c990131         fl(jf)=amax1(flmin,fl(jf))
         fl(jf)=max(flmin,fl(jf))
      enddo

 999  return
      end subroutine fle_fsa
c
c
      subroutine fle(setup,lp) !called by sourceko, for each lr_
      implicit integer (i-n), real*8 (a-h,o-z)
      character*(*) setup
c
c.......................................................................
c  At each poloidal angle, specified by poloidal index lp:
c  Computes a local reduced distribution function,  
c  fl(0:jfl)=<f(u)>, where <> is a pitch angle average.
c  The fl are computed centered on the velocity grid
c  xl(1:jfl)=x(2:jx-1).  Thus the xl(1:jfl) give bin bndries.
c  fl(1:jfl-1) are the relevant (non-zero) entries, each centered
c    at normalized velocity xlm(1:jfl-1).  
c
c  Normalization of the output distribution fl is such that it equals
c  vnorm*f_cgs(v_parallel).
c.......................................................................
      save
      include 'param.h'
      include 'comm.h'

c     Diagnostic array:
c      dimension den_of_s(lza)  !  !Now comm.h:den_of_s1

      character*8 ifirst
      data ifirst/"first"/

      if (lrz.gt. lz) stop 'in fle:  Need to adjust lrz.le.lz'

      if (setup.eq."setup") then
c     At the first call, set up the mesh.  This is independent of lp.
      if (ifirst.eq."first") then
         ifirst="notfirst"

c     Set up the mesh:

ccc      jfl=jx  ! YuP: Why re-set here? Could result in access viol.
cBH180717: Resetting jfl=jx in ainsetva.  Otherwise, following
cBH180717: dcopy can give an overwrite for jfl<jx (which is not 
cBH180717: detected with gdb/ddd bounds checking).
cBH180717:      jflh=0  Not used.

      call dcopy(jx,x,1,xl(1),1)
      call dcopy(jx,dx,1,dxl(1),1)

      go to 999

      endif ! ifirst.eq."first"
      endif ! setup.eq."setup"


c  Form the reduced distribution function, for a given pol. angle lp
c  and radial bin l_:

c  We assume that lp starts at 1 on each flux surface l_, and initialize
c  a density diagnostic (for use with debugger) to zero.
      if (lp.eq.1) then
         do l=1,lza
            den_of_s1(l)=0.0
         enddo
         denfl1(l_)=0.0
      endif

c     (Temporarily) gave up on the following commented out approach.
c     One thing to watch out for is missing density when preloading
c       at i=1. (BobH 970721)
c        call bcast(fl(0),zero,jfl+1)
c  c      iii=min(imax(lp,lr_)+1,iyh_(l_))
c  c     Here we relie on that yz=pi/2, exactly, beyond the trapping
c  c     boundary in the equatorial plane.  This give d(yz)=0 at these i.
c  c     Also relie here on xl=x mesh, using cint2 and dx
c
c        do jf=1,jfl
c           do i=2,iy_(l_)
c              fl(jf)=fl(jf)+
c       1           twopi*cint2(jf)/dx(jf)*(yz(i,lp,lr_)-yz(i-1,lp,lr_))*
c       2           0.25*(sinz(i,jf,lr_)+sinz(i-1,jf,lr_))*
c       3           (f(i,jf,1,l_)+f(i-1,jf,1,l_))
c           enddo
c        enddo


      call bcast(fl1(0),zero,jfl+1)
      call bcast(fl2(0),zero,jfl+1)
      iii=min(imax(lp,lr_)+1,iyh_(l_))
      do jf=1,jfl
         if(jf.le.jx) then !YuP:Now jfl cannot be larger than jx(ainsetva,L~871)
         do i=1,iii
            fl1(jf)=fl1(jf)+cynt2(i,l_)*dtau(i,lp,lr_)*
     +           abs(coss(i,l_))*f_(i,jf,1,l_)
         enddo
         do i=1,iii
            ii=iy_(l_)+1-i
            fl2(jf)=fl2(jf)+cynt2(ii,l_)*dtau(ii,lp,lr_)*
     +           abs(coss(ii,l_))*f_(ii,jf,1,l_)
         enddo
         endif ! if(jf.le.jx)
         fl1(jf)=fl1(jf)*cint2(jf)*bbpsi(lp,lr_)/(dxl(jf)*dz(lp,lr_))
         fl2(jf)=fl2(jf)*cint2(jf)*bbpsi(lp,lr_)/(dxl(jf)*dz(lp,lr_))
      enddo

      fl1(1)=fl1(1)+fl1(0)
      fl1(jfl-1)=fl1(jfl-1)+fl1(jfl)
      fl1(0)=0.0
      fl1(jfl)=0.0
      fl2(1)=fl2(1)+fl2(0)
      fl2(jfl-1)=fl2(jfl-1)+fl2(jfl)
      fl2(0)=0.0
      fl2(jfl)=0.0

c  Set minimum fl = 1.e-100* (max. value)
      flmax=0.
      do jf=1,jfl-1
c990131         flmax=amax1(flmax,fl1(jf))
c990131         flmax=amax1(flmax,fl2(jf))
         flmax=max(flmax,fl1(jf))
         flmax=max(flmax,fl2(jf))
      enddo
      flmin=em100*flmax
      do jf=1,jfl-1
c990131         fl1(jf)=amax1(flmin,fl1(jf))
c990131         fl2(jf)=amax1(flmin,fl2(jf))
         fl1(jf)=max(flmin,fl1(jf))
         fl2(jf)=max(flmin,fl2(jf))
      enddo

c     Diagnostic check of densities:
      do jf=1,jfl-1
         den_of_s1(lp)=den_of_s1(lp)+dxl(jf)*(fl1(jf)+fl2(jf))
      enddo

      denfl1(l_)=denfl1(l_)
     1          +(dz(lp,lr_)/bbpsi(lp,lr_))*den_of_s1(lp)/zmaxpsi(lr_)

 999  return
      end subroutine fle
