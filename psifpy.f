c
      real*8 function psifpy(yval)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c------------------------------------------------------------
c
cmnt  calculates the derivative of bbpsi with respect to arc
cmnt  length along a field line-as a function of bbpsi:
cmnt  psifpy(bbpsi)=dbbpsi/ds (bbpsi)
c
cBH091106:  Only called from psifp, and only for psimodel.eq."axitorus".
c
c------------------------------------------------------------

      include 'param.h'
      include 'comm.h'
      
      iupdown=1

      if (machine.eq."toroidal") then
        if (psimodel.eq."axitorus") then
          psifpy=sqrt((yval-1.)*(1.-(yval/psimx(lr_))))*yval*
     *      pibzmax(lr_)
        else if (psimodel.eq."spline") then !But only called with axitorus
          xt=psiinv(yval,iupdown)
          psifpy=psifp(xt)
        endif
      else if (machine.eq."mirror") then
        if (psimodel.eq."axitorus") then
          psifpy=sqrt((yval-1.)*(1.-(yval/psimx(lr_))))*yval*
     *      pibzmax(lr_)
        else if (psimodel.eq."smith") then
          xz=psiinv(yval,iupdown)
          psifpy=psifp(xz)
        else if (psimodel.eq."spline") then !But only called with axitorus
          xt=psiinv(yval,iupdown) !BUG, missing iupdown. [2019-07-08]
          psifpy=psifp(xt)
        endif
      endif
      return
      end
