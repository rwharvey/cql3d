c
c
      subroutine eqrhs(neq,t,yv,ydot)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'

      dimension yv(2), ydot(2)
c..................................................................
c     This routine is used with O.D.E. solver LSODE. It provides
c     the r.h.s. of the set of coupled O.D.E.'s .
c..................................................................

      dpsidr=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz,
     1  epsirz,nnra,1,0)
      dpsidz=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz,
     1  epsirz,nnra,0,1)
      rbval=sqrt(dpsidr**2+dpsidz**2+(fpsi_)**2)
c
c     Calc Br and Bz for integration along the  poloidal flux surface.
c     But, omit cursign factor, so start in pos-Z dirn.
c
      ydot(1)=dpsidz/(rbval)
      ydot(2)=-dpsidr/(rbval)
      return
      end
