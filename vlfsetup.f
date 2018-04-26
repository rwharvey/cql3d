c
c
      subroutine vlfsetup
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     This routine sets up arrays for vlf urf modules
c..................................................................

      include 'param.h'
      include 'comm.h'

c...................................................................
c     mrf is the number of RF modes  being calculated.
c     nharms(1:mrf)) is the number of cyclotron harmonics calculated,
c     starting at nharm1.
c     Both nharms().gt.1 and mrf.gt.1 are now permitted by  the
c     storage scheme [BH, 060314].
c...................................................................

      pointer besl
      dimension besl(:) ! Now it's a local working array

      mrf=vlfmodes   !i.e., number of wave types

      if (mrf.gt.nmodsa) then
         write(*,*)'STOP:    vlfmodes.gt.nmodsa'
      endif

      deps=1.e-10
      do k=1,mrf
         nharms(k)=vlfharms(k)+deps
         nharm1(k)=vlfharm1(k)+deps
      enddo

cBH060314      one=1.
cBH060314      nharms=max(one,vlfharms)

c     mrfn is the number of wave "modes", that is, the sum over
c     wave types of the number of harmonics for each wave type
      mrfn=0
      do k=1,mrf
         if (nharms(k).eq.0) then
            mrfn=mrfn+1
         else
            mrfn=mrfn+nharms(k)
         endif
      enddo

      if (mrfn.gt.nmodsa) then
         write(*,*)'vlfsetup: mrfn>nmodsa.  mrfn,nmodsa=',mrfn,nmodsa
         STOP 'Increase nmodsa.'
      endif

cBH020913, added following line.
cBH060314      nharm1=vlfharm1(1)

cBH060314      mrfn=max0(mrf,nharms)

cBH060314c     Test for non-implemented options:
cBH060314      if (vlfmodes.gt.1.and.vlfharms.gt.1) then 
cBH060314        print 102
cBH060314 102    format(' vlfmodes.gt.1.and.vlfharms.gt.1') 
cBH060314        stop ' in vlfsetup'
cBH060314      endif
cBH060314      if (nharms.gt.nmodsa)  stop 'nharms.gt.nmodsa'

c.......................................................................
c
c     Set up table krfn(1:mrfn) pointing to wave type index.
c     Set up table irfn(1:mrf) pointing to the wave mode index of
c       the lowest harmonic for each wave type.
c     Set up nharm(1:mrfn).
c     NOTE:  Cannot set these rrays earlier in code because
c            nharm is possibly read above from ray data files,
c            giving value for nharm1.
c
c.......................................................................
      
      k=0
      do krf=1,mrf
         do kk=1,nharms(krf)
            k=k+1
            krfn(k)=krf
            if (kk.eq.1) irfn(krf)=k
            nharm(k)=nharm1(krf)+(kk-1)
         enddo
      enddo

      write(*,*)
      write(*,*)'vlfsetup: mrf = ',mrf
      write(*,*)'vlfsetup: mrfn = ',mrfn
      write(*,*)'vlfsetup: irfn = ',irfn
      write(*,*)'vlfsetup: krfn = ',krfn
      write(*,*)'vlfsetup: nharm1 = ',nharm1
      write(*,*)'vlfsetup: nharms = ',nharms
      write(*,*)'vlfsetup: nharm = ',nharm
 

c     Duplicate ray data in sets nharms(krf) long using the first ray data
c     for each ray type krf, if there is more than one harmonic:
c     Need to do this duplication in reverse order, k=mrfn,1,-1
c     in order not to overwrite data initially stored as krf=1,mrf.

      do k=mrfn,1,-1
         freqcy(k)=vlffreq(krfn(k))
         omega(k)=twopi*vlffreq(krfn(k))
         vlfnperp(k)=vlfnperp(krfn(k))
         vlfnp(k)=vlfnp(krfn(k))
         vlfdnp(k)=vlfdnp(krfn(k))
         vlfddnp(k)=vlfddnp(krfn(k))
         vlfeplus(k)=vlfeplus(krfn(k))
         vlfemin(k)=vlfemin(krfn(k))
         vlfpol(k)=vlfpol(krfn(k))*pi/180.
         vlfdpol(k)=vlfdpol(krfn(k))*pi/180.
         vlfddpol(k)=vlfddpol(krfn(k))*pi/180.
         vlfdnorm(k)=vlfdnorm(krfn(k))
      enddo

      write(*,*)'vlfsetup: vlfdnorm(k),k=1,mrfn ',
     1                    (vlfdnorm(k),k=1,mrfn)
      
      vlfpol_inrange=0.d0
      !YuP[03-2016] Scan all poloidal angles along flux surface:
      ! if vlfpol is outside of range of all pol(), print warning.
      !Note: in a mirror machine the range of pol() is limited, 
      !      usually less than [-90;+90] degrees.
      ! In a tokamak (closed surfaces), it is [0;180] or [0;360] degrees
      ! depending on eqsym value.
      do k=1,mrfn
      do l=1,lz
      if( abs(pol(l,lr_)-vlfpol(k)) .lt. 0.5*vlfdpol(k) )then
        ! ok, this pol() is near vlfpol, within the +/- vlfdpol/2 range
        vlfpol_inrange=1.0 ! value is changed
      endif
      enddo
      enddo
      if(vlfpol_inrange.eq.0.d0)then 
        ! the above condition was not met for any vlfpol(k)
        WRITE(*,*)'vlfsetup: vlfpol is outside of pol.angle range.'
        WRITE(*,*)'vlfsetup: pol(1,lrz)=',pol(1,lrz),
     +                    '  pol(lz,lrz)=',pol(lz,lrz)
        WRITE(*,*)'vlfsetup: vlfpol(k)=',vlfpol ! print for all k
        stop
      endif
      
c..................................................................
c     Allocate arrays
c..................................................................
      jjx=((jx-1)/ibytes)*ibytes+ibytes

      call vlfalloc

      do 15 j=1,jx
        sx(j)=x(j)/gamma(j)
 15   continue

c..................................................................
c     "l" refers to poloidal position z(l,lr_)
c     lr_ is the flux surface label.
c..................................................................

      do 20 l=1,lz
        do 21 i=1,iyh
          cosmz(i,l,lr_)=-.5*(cosz(i+1,l,lr_)+cosz(i,l,lr_))
          cosmz(iy+1-i,l,lr_)=-cosmz(i,l,lr_)
 21     continue
        cosmz(iy,l,lr_)=-cosz(iy,l,lr_)


 20   continue

c..................................................................
c     Fill in bessel function table.
c..................................................................
      nharmx = 1
      do k=1,mrfn
         nharmx=MAX(nharmx,nharm(k))
      enddo
      write(*,*)'vlfsetup:  Before allocate'
      allocate(besl(nharmx+2),STAT=istat) !-YuP->added
      call bcast(besl,0.0,SIZE(besl))     !-YuP->added
      allocate(jbm1(nbssltbl,mrfn),STAT=istat) !-YuP-> modified to real
      call bcast(jbm1,0.0,SIZE(jbm1))          !-YuP-> modified to real
      allocate(jb0(nbssltbl,mrfn),STAT=istat)  !-YuP-> modified to real
      call bcast(jb0,0.0,SIZE(jbm1))           !-YuP-> modified to real
      allocate(jbp1(nbssltbl,mrfn),STAT=istat) !-YuP-> modified to real
      call bcast(jbp1,0.0,SIZE(jbm1))          !-YuP-> modified to real
      write(*,*)'vlfsetup:  After allocate'
c..................................................................
c     Loop over excitation modes
c..................................................................

      do 200 k=1,mrfn

c..................................................................
c     Find maximum k-perp, accounting for aspect ratio
c..................................................................

cBH000416        xkp=vlfnperp(k)*freqcy(k)*twopi/clight*solr(1,lr_)/
cBH000416     +        solr(lorbit(lr_),lr_)
        xkp=vlfnperp(k)*freqcy(k)*twopi/clight*solrz(1,lr_)/
     +        solrz(lz,lr_)

c..................................................................
c     Guess at minimum cyclotron frequency - leave bpol out of estimate
c     this will decrease guess.
c..................................................................

        bvalmin=fpsi(lrindx(lrors))/solr(lorbit(lr_),lr_)
        !YuP[03-2016] But in a mirror machine, Btor=0 (fpsi=0).
        !write(*,*)'vlfsetup: krf,bvalmin=',k,bvalmin
        bvalmin=bmidplne(lr_) !YuP[03-2016]
        
        wcemin=abs(bnumb(1))*charge*bvalmin/(fmass(1)*clight)

c..................................................................
c     Maximum argument in table will be argmax - give it a 40% boost
c..................................................................

        argmax=xkp*vnorm/wcemin*1.4
        write(*,'(a,3i6,2e12.3)')'vlfsetup: n,lr_,krf,bvalmin,argmax=',
     +    n,lr_,k,bvalmin,argmax
        bsslstp(k)=argmax/(nbssltbl-1)
        arg=0.
        do 30 i=1,nbssltbl
          call zzbeslri(arg,nharm(k)+2,0,besl,ncalc)
          if (ncalc.ne.nharm(k)+2) stop ' in vlfsetup'
          if(nharm(k).eq.0)  then
            jbm1(i,k)=-besl(2)
            jb0(i,k)=besl(1)
            jbp1(i,k)=besl(2)
          elseif (nharm(k).gt.0)  then
            jbm1(i,k)=besl(nharm(k))
            jb0(i,k)=besl(nharm(k)+1)
            jbp1(i,k)=besl(nharm(k)+2)
          endif
          arg=arg+bsslstp(k)
 30     continue
        argmax=arg-bsslstp(k)
        
 200  continue ! k=1,mrfn
      
      return
      end
