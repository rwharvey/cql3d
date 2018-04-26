c
      real*8 function psifppy(yval,iupdown)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c------------------------------------------------------------
c
cmnt  function psifppy(psi)
cmnt  calculates the second derivative of psi w/r/t arc length
cmnt  at the midplane or the throat as a function of psi
cmnt  psifppy(psi)=d(dpsi/ds)/ds (psi)

c     NOTE:  Not presently used.
c            If to be re-commisioned for some purpose,
c            update/check for use with trapmod="enabled".
c            Also, update for iupdown (for eqsym.eq."none").
c
c
c
c------------------------------------------------------------

      include 'param.h'
      include 'comm.h'

      itab(1)=0
      itab(2)=0
      itab(3)=1
      if (machine.eq."toroidal") then
        if (psimodel.eq."axitorus") then
          if (yval.eq.1.) then
            psifppy=.5*(1.-1./psimx(lr_))*pibzmax(lr_)**2
          else if (abs(psimx(lr_)-yval).lt.1.e-13) then
            psifppy=-.5*psimx(lr_)*(psimx(lr_)-1.)*pibzmax(lr_)**2
          else
            call diagwrng(6)
          endif
        else if (psimodel.eq."spline") then
          if (yval.eq.1.) then
            xz=0.
          else if (abs(psimx(lr_)-yval).lt.1.e-13) then
            xz=zmax(lr_)
          else
            call diagwrng(6)
          endif
          call terp1(lorbit(lr_),es(1,lr_),bpsi(1,lr_),d2bpsi(1,lr_),
     1      xz,1,tab,itab)
          psifppy=tab(3)
        endif
      else if (machine.eq."mirror") then
        if (psimodel.eq."axitorus") then
          if (yval.eq.1.) then
            psifppy=.5*(1.-1./psimx(lr_))*pibzmax(lr_)**2
          else if (abs(psimx(lr_)-yval).lt.1.e-13) then
            psifppy=-.5*psimx(lr_)*(psimx(lr_)-1.)*pibzmax(lr_)**2
          else
            call diagwrng(6)
          endif
        else if (psimodel.eq."smith") then
          if (yval.eq.1.) then
            gsa3=2.*(1.+gszb2/gslb2)*gseb
            psifppy=2.*(1./gsla2+2.*gsb*gsa3/gslb2)/gsnm
          else if (abs(psimx(lr_)-yval).lt.1.e-13) then
            gsa3=(1.+gszm2/gslb2)*gsem+(1.+gszp2/gslb2)*gsep
            psifppy=2.*(1./gsla2+2.*gsb*gsa3/gslb2)/gsnm
          else
            call diagwrng(6)
          endif
        else if (psimodel.eq."spline") then
          if (yval.eq.1.) then
            xz=0.
          else if (abs(psimx(lr_)-yval).lt.1.e-13) then
            xz=zmax(lr_)
          else
            call diagwrng(6)
          endif
          call terp1(lorbit(lr_),es(1,lr_),bpsi(1,lr_),d2bpsi(1,lr_),
     1      xz,1,tab,itab)
          psifppy=tab(3)
        endif
      endif
      return
      end
