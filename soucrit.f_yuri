c
c      
      subroutine soucrit
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'
c
c   scchiu, 9609..
c  calculate the critical momentum-per-mass.
c  Since jxcrit is used for calculation of relativistic runaways,
c   it is chosen to be no larger than for particles at 3.*clight.
c   This is to reduce a problem of definition which occurs
c   transiently for abs(E) < abs(E_crit_for_runaway).
c
      call cfpgamma
      do 10 k=1,ngen
        fack=abs(elecfld(lr_))/reden(k,lr_)*1.e16/0.0918
     1              *18./gama(k,k)
        eoe0(k,lr_)=fack
        fack1=1./(fack-1.)
        if (fack1.le.0.) then
          ucrit(k,lr_)=1.e20
        else
          ucrit(k,lr_)=clight*sqrt(fack1)/vnorm
        endif
c  Take runaway electrons to have momentum per mass beyond the
c  minimum of 3.*clight or ucrit:
c990131        xcrit=amin1(3.*clight/vnorm,ucrit(k,lr_))
        xcrit=min(3.*clight/vnorm,ucrit(k,lr_))
        jxcrit(k,lr_)=luf(xcrit,x,jx)
10    continue
      return
      end
