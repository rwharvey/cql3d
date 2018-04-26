c
c
      subroutine sourcpwr(k)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE


c..................................................................
c     This routine computes the Neutral Beam power in Watts/cc
c..................................................................

      call bcast(tam1,zero,jx)
      do 501 i=1,iy
        do 502 j=1,jx
         if(gone(i,j,k,lr_).eq.zero)then !YuP[2017-11-21] added.
          !YuP: added gone() check: ONLY COUNT not-lost particles.
          !YuP: This is similar to the CQL3D-FOW version (added on 03/23/2015)
          tam1(j)=tam1(j)+source(i,j,k,indxlr_)*cynt2(i,l_)*vptb(i,lr_)
          !Note: tam1 is needed for sorpw_nbi below, which is only used 
          !for plots and NetCDF file. 
          !It has no effect on the source operator or solution.
          !The consequence of adding if(gone..) condition is
          !the drop of NBI power profile at plasma edge 
          !(if lossmode is not disabled) because of banana+Larm.rad losses;
          !see tdpltjop, plot 'FSA SOURCE POWER DEN: ...';
          !Also, 'Power from NBI', and 'Power integrated over radius'
          !become smaller (values are printed in *.ps file). 
         endif
 502    continue
 501  continue
      s=0.
      do 503 j=1,jx
        s=s+tam1(j)*tcsgm1(j)*cint2(j)
 503  continue

c..................................................................
c     Source power of species k in Watts/cc averaged over the flux surfa
c..................................................................

      sorpw_nbi(k,lr_)=s*fions(k)*one_*1.6022e-16*zmaxpsii(lr_) !-YuP->NBI source

CMPIINSERT_IF_RANK_EQ_0
      WRITE(*,'(a,2i5,e12.4)')
     +      'sourcpwr: k,lr, sorpw_nbi',
     +                 k,lr_,sorpw_nbi(k,lr_)
      WRITE(*,*)'----------------------------------------------------'
CMPIINSERT_ENDIF_RANK
      
      return
      end
