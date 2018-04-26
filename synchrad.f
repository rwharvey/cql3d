c
c
      subroutine synchrad
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     This routine computes  synchrotron radiation loss coefficients 
c     for electrons.
c     syncrad="disabled","gyro","geom", or "gyrogeom". 
c     The relativistic coefficients have been calculated
c     using the hueristic approach of Bernstein and Baxter.
c     (I.B. Bernstein and D.C. Baxter, Phys. Fluids 24, 
c      p. 108 (1981), Appendix A (not equation given in main text).)
c     Plus
c     synchrotron radiation due to toroidal and poloidal rotation
c     (BH + Paul Parks, March 2000).
c..................................................................

      include 'param.h'
      include 'comm.h'


c     Temporary, for writing:
      dimension synca_gyro(jx),syncd_gyro(jx)
      dimension synca_geom(jx),syncd_geom(jx)

      if (syncrad .eq. "disabled" .or. kelecg .eq. 0) return
      call bcast(synca,zero,iyjx*lrz)
      call bcast(syncd,zero,iyjx*lrz)

      if (syncrad.eq."gyro" .or. syncrad.eq."gyrogeom") then

      coefnt=(2.*(charge*bnumb(kelecg))**4*bmod0(lr_)**2)/
     *  (3.*clight**5*fmass(kelecg)**3)
      do  j=1,jx
         do  i=1,iy
            synca(i,j,indxlr_)=coefnt*vptb(i,lr_)*sinn(i,l_)**2*xcu(j)*
     *           psicu(i,lr_)*gamma(j)
            syncd(i,j,indxlr_)=coefnt*vptb(i,lr_)*xsq(j)*gammi(j)*
     *           sinn(i,l_)**2/coss(i,l_)*
     *           (psisq(i,lr_)-sinn(i,l_)**2*psicu(i,lr_))
         enddo
      enddo

      call bcast(synca_gyro,zero,jx)
      call bcast(syncd_gyro,zero,jx)
      do j=1,jx
c         do i=1,iy
          i=11 !   comment from YuP[2018-01-16] Why i=11 only ?
            synca_gyro(j)=synca_gyro(j)+
c     *                    sinn(i,l_)*dy(i,l_)*synca(i,j,indxlr_)
     *                    synca(i,j,indxlr_)
            syncd_gyro(j)=syncd_gyro(j)+
c     *                    sinn(i,l_)*dy(i,l_)*syncd(i,j,indxlr_)
     *                    syncd(i,j,indxlr_)
c         enddo
      enddo

      endif

      if (syncrad.eq."geom" .or. syncrad.eq."gyrogeom") then

      coefnt1=2.*(charge*bnumb(kelecg))**2/(3.*clight**3*fmass(kelecg))
      bsq=bthr(lr_)**2+btoru(lr_)**2
      bphiob4=(btoru(lr_)**2/bsq)**2
      bthob4=(bthr(lr_)**2/bsq)**2
      reffm2=onovrp(2,lr_)*bphiob4 + pi/areacon(lr_)*bthob4
      
      call bcast(synca_geom,zero,jx)
      call bcast(syncd_geom,zero,jx)
      do  j=1,jx
         do  i=1,iy
            
            sync1=  +coefnt1*vptb(i,lr_)*
     *           xcu(j)*xsq(j)*vnorm2*gammi(j)*reffm2*
     *           (1.-2.*sinn(i,l_)**2*psiba(i,lr_)+
     *           sinn(i,l_)**4*psisq(i,lr_))
            synca(i,j,indxlr_)=synca(i,j,indxlr_)+sync1
            sync2=  -coefnt1*vptb(i,lr_)*
     *           xcu(j)*x(j)*vnorm2*gammi(j)*reffm2*
     *           sinn(i,l_)**2/coss(i,l_)*
     *           (1.-2.*sinn(i,l_)**2*psiba(i,lr_)+
     *           sinn(i,l_)**4*psisq(i,lr_))
            syncd(i,j,indxlr_)=syncd(i,j,indxlr_)+sync2

            if (i.eq.1) then
            synca_geom(j)=synca_geom(j)+
c     *                    sinn(i,l_)*dy(i,l_)*sync1
     *                    sync1
            syncd_geom(j)=syncd_geom(j)+
c     *                    sinn(i,l_)*dy(i,l_)*sync1
     *                    sync2
            endif
        enddo
      enddo

      endif

      write(*,*) 'syncrad: synca_gyro: ',synca_gyro
      write(*,*) 'syncrad: syncd_gyro: ',syncd_gyro
      write(*,*) 'syncrad: synca_geom: ',synca_geom
      write(*,*) 'syncrad: syncd_geom: ',syncd_geom


      return
      end
