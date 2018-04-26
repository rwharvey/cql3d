c
c
      subroutine pltcycl(iymn,iymx,ymn,ymx)
      implicit integer (i-n), real*8 (a-h,o-z)
c-------------------------
c  find the minimum and maximum log-cycle for ymn,ymx
c-----------------------------
      iymn=0
      iymx=0
      imn=0
      imx=0
10    imn1=imn+1
      imn2=imn-1
      if (ymn.lt.10.**imn) then
        iymn=imn2
        imn=imn2
        if (ymn.ge.10.**imn2) then
          go to 20
        else
          go to 10
        endif
      endif
      if (ymn.ge.10.**imn) then
        if (ymn.ge.10.**imn1) then
          iymn=imn1
          imn=imn1
          go to 10
        else
          go to 20
        endif
      endif
20    imx1=imx+1
      imx2=imx-1
      if (ymx.gt.10.**imx) then
        iymx=imx1
        imx=imx1
        if (ymx.le.10.**imx1) then
          go to 30
        else
          go to 20
        endif
      endif
      if (ymx.le.10.**imx) then
        if (ymx.le.10.**imx2) then
          iymx=imx2
          imx=imx2
          go to 20
        else
          go to 30
        endif
      endif
30    continue
      return
      end
