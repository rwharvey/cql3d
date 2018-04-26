c
c
      subroutine frsmooth(k,curnorm)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'

      if (smooth_ .lt. .005) return

c..................................................................
c     This routine smooths the source deposition current. The weighting
c     function is a Gaussian with a s.deviation (smooth_) which is
c     normally some small fraction (~.05-.1) of the normalized toroidal
c     coordinate radius radmin. This allows the resulting current
c     profile to be far smoother that what is seen with the model turned
c     off (smooth_<.005). This in turn allows for much less noisy
c     p' and ff' input profiles to the MHD equilibrium code.
c     If utilized, array tr3 and scalar curnorm are altered.
c..................................................................


c..................................................................
c Compute the average source: compute the total source current,curnorm
c..................................................................

      curnorm=0.
      smth=smooth_**2*radmin**2
      do 5 l=1,lrz
        tr5(l)=0.
        if (asorz(k,1,l).gt.em120) tr5(l)=1./asorz(k,1,l)
 5    continue
      do 10 l=1,lrz
        tr(l)=0.
        tr1(l)=0.
        do 20 ll=1,lrz
          tr(l)=exp(-(rz(ll)-rz(l))**2/smth)*dvol(l)+tr(l)
          tr1(l)=exp(-(rz(ll)-rz(l))**2/smth)*asorz(k,1,ll)*dvol(l)
     1      +tr1(l)
 20     continue
        tr4(l)=tr1(l)/tr(l)
        curnorm=curnorm+dvol(l)*tr4(l)
        tr3(l)=tr4(l)*tr5(l)
 10   continue
      return
      end
