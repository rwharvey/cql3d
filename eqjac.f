c
c
      subroutine eqjac(neq,t,yv,ml,mu,pd,nrowpd)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'

      dimension yv(neq),pd(nrowpd,neq)
c..................................................................
c     This routine is used with the O.D.E. solver LSODE. It
c     provides the Jacobian of the system of equations.
c..................................................................

c..................................................................
c     Provide the various derivatives of epsi (first and second
c     order).
c..................................................................

      dpsidr=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz,
     1  epsirz,nnra,1,0)
      dpsidz=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz,
     1  epsirz,nnra,0,1)
      d2psidrr=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz,
     1  epsirz,nnra,2,0)
      d2psidrz=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz,
     1  epsirz,nnra,1,1)
      d2psidzz=terp2(yv(1),yv(2),nnr,er,nnz,ez,epsi,epsirr,epsizz,
     1  epsirz,nnra,0,2)

c..................................................................
c     Compute B and it's first order derivatives.
c..................................................................

      alpha=(dpsidr**2+dpsidz**2+(fpsi_)**2)
      sqrtalp=sqrt(alpha)
      dalphadr=2.*(dpsidz*d2psidrz+dpsidr*d2psidrr)
      dalphadz=2.*(dpsidz*d2psidzz+dpsidr*d2psidrz)

c..................................................................
c     Compute the 4 elements of the Jacobian.
c..................................................................

      pd(1,1)=(d2psidrz-dpsidz*dalphadr/alpha)/alpha
      pd(1,2)=(d2psidzz-dpsidz/alpha*dalphadz)/alpha
      pd(2,1)=(-d2psidrr+dpsidr/alpha*dalphadr)/alpha
      pd(2,2)=(-d2psidrz+dpsidr/alpha*dalphadz)/alpha
      return
      end
