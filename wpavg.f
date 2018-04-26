c
c
      subroutine wpavg
      implicit integer (i-n), real*8 (a-h,o-z)

c..............................................................
c     Compute flux surface average of various quantities for
c     CQLP case
c..............................................................

      include 'param.h'
      include 'comm.h'

c.......................................................................
cl    1. Surface average of parallel current and parallel electric field
c     <j_par/R>, <E_par/R> with <a>=int(ds a/B) / int(ds/B)
c.......................................................................

      zcuravg=0.0
      zcuravg2=0.0
      zeleavg=0.0
      zeleavg2=0.0
      z1oravg=0.0
      zflxavg=0.0
      zflxavg2=0.0
      ilr=lrindx(1)
      zelcof=elecfld(ilr)/300.*rmag*fpsi(ilr)/bmidplne(ilr)**2/3.e+09
      zelcof2=elecfld(ilr)/300.*rmag*fpsi(ilr)/bmidplne(ilr)/3.e+09
      do 100 l=1,ls
        zcuravg=zcuravg+dsz(l)*currmtpz(l)/psis(l)/bmidplne(ilr)/
     /    solrs(l)
        zcuravg2=zcuravg2+dsz(l)*currmtpz(l)
        zeleavg=zeleavg+dsz(l)*zelcof/(solrs(l)*psis(l))**2/solrs(l)
        zeleavg2=zeleavg2+dsz(l)*zelcof2/solrs(l)**2/psis(l)
        z1oravg=z1oravg+dsz(l)/psis(l)/bmidplne(ilr)/solrs(l)
        zflxavg=zflxavg+dsz(l)/psis(l)/bmidplne(ilr)
        zflxavg2=zflxavg2+dsz(l)
 100  continue
      zcuravg=zcuravg/z1oravg
      zeleavg=zeleavg/z1oravg

      write(6,'(/" surface averages:  <j_par/R>/<1/R>= ",1pe13.4,/
     +  "                    <E_par/R>/<1/R>= ",1pe13.4,/
     +  "                    <j_par*B>      = ",1pe13.4,/
     +  "                    <E_par/B>      = ",1pe13.4,/
     +  "        <E_par/R>/<j_par/R>/sptz(1)= ",1pe13.4,/
     +  "        <E_par*B>/<j_par*B>/sptz(1)= ",1pe13.4,/
     +  "  <1/R>= ",1pe13.4,"   flxavg= ",1pe13.4,
     +  " flxavg2= ",1pe13.4,"  n=",i4)')
     +  zcuravg,zeleavg,zcuravg2,zeleavg2,zeleavg/zcuravg/sptzr(1),
     +  zeleavg2/zcuravg2/sptzr(1),z1oravg,zflxavg,zflxavg2,n

      return
      end
