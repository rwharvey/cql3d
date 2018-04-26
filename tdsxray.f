c
c
      subroutine tdsxray(icall,iplotsxr)
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     sets up call to soft-x-ray analyzer
c..................................................................

      include 'param.h'
      include 'comm.h'
      character*8 icall,iplotsxr

      do 1 l=1,lrzmax
        tr1(l)=reden(kelec,l)
 1    continue
cBH081106:  In some radcoord cases, rrz contains normalized radial
cBH081106:  coord data, and is not suitable for eqmod.ne."enabled"
cBH081106:  circ plasma model, or the eqmod.eq."enabled" eqdsk 
cBH081106:  equilibria.  
      if (eqmod.ne.'enabled') then
         if (radcoord.ne.'sqtorflx')
     +      write(*,*)'tdsxray: WARNING, check our radial coord further'
         call tdsxr0(rrz,tr1(1),icall,iplotsxr)
      else
         call tdsxr0(rpmconz,tr1(1),icall,iplotsxr)
      endif
      return
      end
