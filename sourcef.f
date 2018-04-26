c
c
      subroutine sourcef
      implicit integer (i-n), real*8 (a-h,o-z)
      save
      include 'param.h'
      include 'comm.h'


c..................................................................
c     If soucoord.ne."disabled",
c     computes the source m of species k - simple Gaussian profiles.
c     If knockon.ne."disabled",
c     also calc. knock-on source of electrons (Besedin and Pankratov).
c..................................................................

      do 900 k=1,ngen
          if (soucoord .ne. "disabled") then
        do 800 m=1,nso
          if (n .lt. nonso(k,m) .or. n .gt. noffso(k,m)) go to 800
          call bcast(temp2(0,0),zero,iyjx2)
          do 700 l=1,lz

c..................................................................
c     Determine the local source at z(l,lr_) in preparation for bounce-av
c..................................................................

            do 650 i=1,imax(l,lr_)
              ii=iy+1-i
              call soup(cosz(i,l,lr_),l,k,m)
              do 630 j=1,jx
                temp1(i,j)=soupp(j,lr_)
 630          continue
              call soup(cosz(ii,l,lr_),l,k,m)
              do 610 j=1,jx
                temp1(ii,j)=soupp(j,lr_)
 610          continue
 650        continue

c..................................................................
c     Do the bounce average of the source..
c..................................................................

            do 500 i=1,imax(l,lr_)
              ii=iy+1-i
              ax=dtau(i,l,lr_)/tau(i,lr_)
              if ((l.ne.lz .and. lmax(i,lr_).eq.l)) then
                ax=ax+dtau(i,l+1,lr_)/tau(i,lr_)
              endif
              do 450 j=1,jx
                temp2(i,j)=temp2(i,j)+ax*temp1(i,j)
                temp2(ii,j)=temp2(ii,j)+ax*temp1(ii,j)
 450          continue
 500        continue
 700      continue  ! On l

c..................................................................
c     Zero out source where less than small fraction of the
c     peak.  This facilitates use of ineg='enabled' option
c     zeroing out f where f<0 occurs, for all j above max of
c     the source.   BH120317.
c..................................................................

          temp2frac=1.d-8
          temp2max=zero
          do j=1,jx
             do i=1,iy
                temp2max=max(temp2max,temp2(i,j))
             enddo
          enddo
          temp2frac=temp2frac*temp2max
          do j=1,jx
             do i=1,iy
                if (temp2(i,j).lt.temp2frac) temp2(i,j)=zero
             enddo
          enddo
          

c..................................................................
c     Compute the density in preparation for scaling to achieve
c     the desired current.
c..................................................................

          s1=0.
          do 300 i=1,iy
            do 250 j=1,jx
              s1=s1+temp2(i,j)*cynt2(i,l_)*cint2(j)*vptb(i,lr_)
 250        continue
 300      continue

c..................................................................
c     Determine q1 so that the flux surface averaged current
c     will be asor (particles/sec/cc) after scaling.
c..................................................................

          if (s1.ne.zero) q1=asor(k,m,lr_)/s1*zmaxpsi(lr_)/one_

c..................................................................
c     scale the current to be equal to asor
c..................................................................

          call dscal(iyjx2,q1,temp2(0,0),1)
          do 200 j=1,jx
            do 150 i=1,iy
              source(i,j,k,indxlr_)=source(i,j,k,indxlr_)+temp2(i,j)
 150        continue
 200      continue

c..................................................................
c     xlncur is in units of particles/cm**2 (field line density)
c..................................................................

          xlncur(k,lr_)=xlncur(k,lr_)+asor(k,m,lr_)*zmaxpsi(lr_)
          xlncurt(lr_)=xlncurt(lr_)+xlncur(k,lr_)
 800    continue
      endif

 900  continue
      return
      end
