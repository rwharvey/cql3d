c
c
c
      real*8 function psif(xz)
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c-------------------------------------------------------------------
c
cmnt  function psif(arclength)
cmnt  calculates the ratio of mod B at arclength from the midplane
cmnt  to mod B at the midplane; 1 < psi < psimx(lr_).
cmnt  the value of psimodel controls which model is expressed:
cmnt  options are:
cmnt  psimodel="axitorus", "smith", "spline"
c
c     OK with eqsym.eq."none".  xz is measured from es(1) and in 
c       this case covers a whole poloidal turn.
c
c-------------------------------------------------------------------
      include 'param.h'
      include 'comm.h'

      dimension tabl(3),itabl(3)
      data itabl(1) /1/, itabl(2) /0/, itabl(3) /0/
      if (machine.eq."toroidal") then
        if (psimodel.eq."axitorus") then
          psif=(1.+eps(lr_))/(1.+eps(lr_)*cos(pi*xz*zmaxi(lr_)))
        else if (psimodel.eq."spline") then
          call terp1(lorbit(lr_),es(1,lr_),bpsi(1,lr_),d2bpsi(1,lr_),
     1      xz,1,tabl,itabl)
          psif=tabl(1)
        endif
      else if (machine.eq."mirror") then
        if (psimodel.eq."axitorus") then
          psif=(1.+eps(lr_))/(1.+eps(lr_)*cos(pi*xz*zmaxi(lr_)))
        else if (psimodel.eq."smith") then
          gsexm=exp(-((xz-gszb)/gslb)**2)
          gsexp=exp(-((xz+gszb)/gslb)**2)
          psif=(1.+(xz/gsla)**2+gsb*(gsexm+gsexp))/gsnm
        else if (psimodel.eq."spline") then
          call terp1(lorbit(lr_),es(1,lr_),bpsi(1,lr_),d2bpsi(1,lr_),
     1      xz,1,tabl,itabl)
          psif=tabl(1)
        endif
      endif
      if (trapmod.eq."enabled") then
c        Reduce B/B_min by the amount trapredc*(B/B_min-1.)
c        This is for purposes of a heuristic check of trapping effects.
         psif=psif-trapredc*(psif-1.)
      endif
      return
      end
