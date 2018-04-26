
      real*8 function psifp(xz)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c-------------------------------------------------------------------
c
cmnt  function psifp(arclength)
cmnt  calculates the derivative of psi=B/B0 w.r.t. arclength, as a 
cmnt  function of arclength: psifp(s)=dpsi/ds (s)
cmnt  model fields are as for psif.
c
c     OK with eqsym.eq."none".  xz is measured from es(1) and in 
c       this case covers a whole poloidal turn.  The function
c       will be negative for xz>es_bmax(lr_).
c
c-------------------------------------------------------------------
      include 'param.h'
      include 'comm.h'

      itab(1)=0
      itab(2)=1
      itab(3)=0
      if (machine.eq."toroidal") then
        if (psimodel.eq."axitorus") then
          psifp=psifpy(psif(xz))
        else if (psimodel.eq."spline") then
          call terp1(lorbit(lr_),es(1,lr_),bpsi(1,lr_),d2bpsi(1,lr_),
     1      xz,1,tab,itab)
          psifp=tab(2)
          if (trapmod.eq."enabled") then
c     Reduce B/B_min by the amount trapredc*(B/B_min-1.)
c     This is for purposes of a heuristic check of trapping effects.
             psifp=psifp-trapredc*psifp
          endif
       endif
      else if (machine.eq."mirror") then
         if (psimodel.eq."axitorus") then
            psifp=psifpy(psif(xz))
         else if (psimodel.eq."smith") then
            gszxm=xz-gszb
            gszxp=xz+gszb
            gsexm=exp(-gszxm**2/gslb2)
            gsexp=exp(-gszxp**2/gslb2)
            psifp=2.*(xz/gsla2-gsb*(gszxm*gsexm+gszxp*gsexp)/gslb2)/gsnm
         else if (psimodel.eq."spline") then
            call terp1(lorbit(lr_),es(1,lr_),bpsi(1,lr_),d2bpsi(1,lr_),
     1           xz,1,tab,itab)
            psifp=tab(2)
            if (trapmod.eq."enabled") then
c     Reduce B/B_min by the amount trapredc*(B/B_min-1.)
c     This is for purposes of a heuristic check of trapping effects.
               psifp=psifp-trapredc*psifp
            endif
        endif
      endif
      return
      end
