
c
c
      subroutine eqelpse
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c
      include 'param.h'
      include 'comm.h'
      character*8 ifirst


c..................................................................
c     This routine generates an ad-hoc psi (an "solution" to
c     the equilibrium calculation)  with elliptical contours of psi
c     centered at rmag with ellipticity = ellptcty. The funtional form
c     of psi(E,lr_) where E**2=z**2+(r-rmag)**2/(1.-ellptcty**2) is
c     arbitrary and is handled by function eqfn. The psi function
c     is normalized so that the poloidal field at r=rmag, z=radmin
c     is bth.
c..................................................................

c..................................................................
c     Create the z,r meshes..
c..................................................................

      data ifirst /"first"/
      if (lr_.ne.lrzmax) return
      if (ifirst.eq."first") then
        zst=zbox/(nnz-1)
        rst=rbox/(nnr-1)
        er(1)=rboxdst
        ez(1)=-zbox*.5
        do 10 nn=2,nnr
          er(nn)=er(nn-1)+rst
 10     continue
        do 11 nn=2,nnz
          ez(nn)=ez(nn-1)+zst
 11     continue
        zmag=zero
        ifirst="notfirst"
      endif

c..................................................................
c     Set up the psi array. First determine the normalization
c     constant from the poloidal field constraint.
c..................................................................

      fhp=eqfn(radmin*1.0001,one)
      fhm=eqfn(radmin*.9999,one)
      deriv=(fhp-fhm)/(radmin*2.e-4)
      scalfct=-bth*rmag/deriv

c..................................................................
c     scalfct is the factor, now for the psi array..
c..................................................................

      do 20 ir=1,nnr
        do 21 iz=1,nnz
          epsi(ir,iz)=sqrt(ez(iz)**2+(er(ir)-rmag)**2/(1.-ellptcty**2))
 21     continue
 20   continue
      do 30 ir=1,nnr
        do 31 iz=1,nnz
          epsi(ir,iz)=eqfn(epsi(ir,iz),scalfct)
 31     continue
 30   continue
      return
      end
