c
c
      subroutine wpbdry(kiorj,kspec,kleft,kright,kopt)
      implicit integer (i-n), real*8 (a-h,o-z)

c..............................................................
c     Determine boundary conditions for parallel transport
c
c     kopt = 1: y mesh assumed constant and kiorj is the j index
c     2: y mesh non constant (over l), kiorj is the i index
c     with i.le.iyh
c
c     Note: fnhalf is the value at midpoint l+1/2.
c     For mod(nummods,10).le.4 the value at nodes l is known in f
c     mod(nummods,10)  > 4: linear or quadr. interpolation
c
c     Backward and centered schemes : impose cond. at l=1
c     Forward schemes : impose cond. at l=ls
c
c     Assumes for
c     numindx = 1: backward scheme
c     2: back/forward for cos(theta)>0 and for/backward for i>iyh
c     depending if numixts=-1,+1 resp.
c     3: forward
c     4: centered
c..............................................................

      include 'param.h'
      include 'comm.h'

      include 'wpadvnc.h'
      flin(dx1,dx2,f1,f2)=f1-(f2-f1)*dx1/dx2
c.......................................................................

      iband=kleft+kright+1

      if (kopt .eq. 2) go to 500

c***********************************************************************
c*********************KOPT= 1  (case y-mesh constant) ***************
c***********************************************************************

      if (sbdry .eq. "periodic") go to 200

c.......................................................................
cl    1. Bounded mesh
c.......................................................................

      ivelmid=1
      if (mod(nummods,10) .le. 4) ivelmid=0
      if (sbdry .eq. "zeroslop") go to 120
      if (sbdry .eq. "pseudper") go to 130

c.......................................................................
cl    1.1 Fixed boundary case: f_n+1=f_n+1/2
c.......................................................................

      if (numindx .eq. 3) go to 1120

c     left boundary: l=1
      iystart=1
      iyend=iy_(1)
      if (numindx .eq. 2) then
        if (numixts .eq. +1) iystart=iyh_(1)+1
        if (numixts .eq. -1) iyend=iyh_(1)
      endif
      do 110 i=iystart,iyend
        bndmats(i,1,1,1)=0.0
        bndmats(i,1,kleft+1,1)=1.0
        bndmats(i,1,iband,1)=0.0
        rhspar(1,i,1)=cvmgt(f(i,kiorj,kspec,1),flin(dszp5(1),
     1    dszp5(1)+dszp5(2),fnhalf(i,kiorj,kspec,1),
     1    fnhalf(i,kiorj,kspec,2)),ivelmid.eq.0)
 110  continue

c     right boundary: l=ls

 1120 continue
      if (numindx.eq.2 .or. numindx.eq.3) then
        iystart=1
        iyend=iy_(ls)
        if (numindx .eq. 2) then
          if (numixts .eq. +1) iyend=iyh_(ls)
          if (numixts .eq. -1) iystart=iyh_(ls)+1
        endif
        do 112 i=iystart,iyend
          bndmats(i,ls,1,1)=0.0
          bndmats(i,ls,kleft+1,1)=1.0
          bndmats(i,ls,iband,1)=0.0
          rhspar(ls,i,1)=cvmgt(f(i,kiorj,kspec,ls),flin(dszm5(ls),
     1      dszm5(ls)+dszm5(ls-1),fnhalf(i,kiorj,kspec,ls-1),
     1      fnhalf(i,kiorj,kspec,ls-2)),ivelmid.eq.0)
 112    continue
      endif

      return

c.......................................................................
cl    1.2 Zero slopes: f(1)=f(2) and/or f(ls)=f(ls-1)
c.......................................................................

 120  continue

      if (numindx .eq. 3) go to 1220

c     left boundary: ls=1
      iystart=1
      iyend=iy_(1)
      if (numindx .eq. 2) then
        if (numixts .eq. +1) iystart=iyh_(1)+1
        if (numixts .eq. -1) iyend=iyh_(1)
      endif
      do 121 i=iystart,iyend
        bndmats(i,1,1,1)=0.0
        bndmats(i,1,kleft+1,1)=1.0
        bndmats(i,1,iband,1)=-1.0
        rhspar(1,i,1)=0.0
c%OS  rhspar(1,i,1)=fnhalf(i,kiorj,kspec,1)+
c%OS  +                                   dtreff*velsou(i,kiorj,kspec,1)
 121  continue

c     right boundary: l=ls

 1220 continue
      if (numindx.eq.2 .or. numindx.eq.3) then
        iystart=1
        iyend=iy_(ls)
        if (numindx .eq. 2) then
          if (numixts .eq. +1) iyend=iyh_(ls)
          if (numixts .eq. -1) iystart=iyh_(ls)+1
        endif
        do 122 i=iystart,iyend
c%OS  bndmats(i,ls,1,1)=0.0
          bndmats(i,ls,1,1)=-1.0
          bndmats(i,ls,kleft+1,1)=1.0
          bndmats(i,ls,iband,1)=0.0
          rhspar(ls,i,1)=0.0
c%OS  rhspar(ls,i,1)=fnhalf(i,kiorj,kspec,ls)+
c%OS  *                               dtreff*velsou(i,kiorj,kspec,ls)
 122    continue
      endif

      return

c.......................................................................
cl    1.3 Fixed boundary case but pseudo-periodic boundary condition
c.......................................................................

 130  continue

      return

c.......................................................................
cl    2. Periodic case
c.......................................................................

 200  continue

c.......................................................................
cl    2.1 Fixed value at midplane, l=1
c     Note: l-1=0 -> l_effective=3 ; l+1=2 -> l_effective=2
c.......................................................................

c%OS  do 210 i=1,iy_(1)
c%OS  bndmats(i,1,iband,1)=0.0
c%OS  bndmats(i,1,kleft+1,1)=1.0
c%OS  bndmats(i,1,kleft+2,1)=0.0
c%OS  rhspar(1,i,1)=fnhalf(i,kiorj,kspec,1)
c%OS  210  continue

      return

c.......................................................................
cl    2.3 Other boundary conditions for periodic case
c.......................................................................

 230  continue

      return

c.......................................................................
cl    3. Up/down symmetric case, assumes periodic mesh
c.......................................................................

 300  continue

      return

c***********************************************************************
c*********************KOPT= 2 (case y-mesh non-constant *************
c***********************************************************************

 500  continue

      if (sbdry .eq. "periodic") go to 600

c.......................................................................
cl    5. Bounded mesh
c.......................................................................

      ivelmid=1
      if (mod(nummods,10) .le. 4) ivelmid=0
      if (sbdry .eq. "zeroslop") go to 520

c.......................................................................
cl    5.1 Fixed boundary case: f_n+1=f_n+1/2
c.......................................................................

      if (numindx .eq. 3) go to 510

c.......................................................................
c     left boundary: l=1 (backward or centered schemes)

c     i<iyh+1 : numindx=1, (2,numixts=-1) or 4
      if (numindx.eq.2 .and. numixts.eq.+1) go to 5111
      do 511 j=2,jx
        bndmats(j,1,1,1)=0.0
        bndmats(j,1,kleft+1,1)=1.0
        bndmats(j,1,iband,1)=0.0
        rhspar(1,j,1)=cvmgt(f(kiorj,j,kspec,1),flin(dszp5(1),
     1    dszp5(1)+dszp5(2),fnhalf(kiorj,j,kspec,1),
     1    fnhalf(kiorj,j,kspec,2)),ivelmid.eq.0)
 511  continue
 5111 continue

c     i>iyh : numindx=1, (2,numixts=+1) or 4
      if (numindx.eq.2 .and. numixts.eq.-1) go to 5121
      iieff=iy_(1)+1-kiorj
      do 512 j=2,jx
        bndmats(j,1,1,2)=0.0
        bndmats(j,1,kleft+1,2)=1.0
        bndmats(j,1,iband,2)=0.0
        rhspar(1,j,2)=cvmgt(f(iieff,j,kspec,1),flin(dszp5(1),
     1    dszp5(1)+dszp5(2),fnhalf(iieff,j,kspec,1),
     1    fnhalf(iieff,j,kspec,2)),ivelmid.eq.0)
 512  continue
 5121 continue

      if (numindx.eq.1 .or. numindx.eq.4) return

c.......................................................................
c     right boundary: l=l_upper(kiorj) (forward schemes)
 510  continue
      ll=l_upper(kiorj)

c     i<iyh+1 : numindx=(2,numixts=+1) or 3
      if (numindx.eq.2 .and. numixts.eq.-1) go to 5131
      do 513 j=2,jx
        bndmats(j,ll,1,1)=0.0
        bndmats(j,ll,kleft+1,1)=1.0
        bndmats(j,ll,iband,1)=0.0
        rhspar(ll,j,1)=cvmgt(f(kiorj,j,kspec,ll),flin(dszm5(ll),
     1    dszm5(ll)+dszm5(ll-1),fnhalf(kiorj,j,kspec,ll-1),
     1    fnhalf(kiorj,j,kspec,ll-2)),ivelmid.eq.0)
 513  continue
 5131 continue

c     i>iyh : numindx=(2,numixts=-1) or 3
      if (numindx.eq.2 .and. numixts.eq.+1) go to 5141
      iieff=iy_(ll)+1-kiorj
      do 514 j=2,jx
        bndmats(j,ll,1,2)=0.0
        bndmats(j,ll,kleft+1,2)=1.0
        bndmats(j,ll,iband,2)=0.0
        rhspar(ll,j,2)=cvmgt(f(iieff,j,kspec,ll),flin(dszm5(ll),
     1    dszm5(ll)+dszm5(ll-1),fnhalf(iieff,j,kspec,ll-1),
     1    fnhalf(iieff,j,kspec,ll-2)),ivelmid.eq.0)
 514  continue
 5141 continue

      return

c.......................................................................
cl    5.2 Zero slopes: f(1)=f(2) and/or f(l_upper(i))=f(l_upper(i)-1)
c.......................................................................

 520  continue

      if (numindx .eq. 3) go to 521

c.......................................................................
c     left boundary: l=1 (backward or centered schemes)

c     i<iyh+1 : numindx=1, (2,numixts=-1) or 4
      if (numindx.eq.2 .and. numixts.eq.+1) go to 5221
      do 522 j=2,jx
        bndmats(j,1,1,1)=0.0
        bndmats(j,1,kleft+1,1)=1.0
        bndmats(j,1,iband,1)=-1.0
c%OS  bndmats(j,1,iband,1)=0.0
        rhspar(1,j,1)=0.0
c%OS  rhspar(1,j,1)=cvmgt(f(kiorj,j,kspec,1)+dtreff*
c%OS  1         velsou(kiorj,j,kspec,1),flin(dszp5(1),dszp5(1)+dszp5(2),
c%OS  1          fnhalf(kiorj,j,kspec,1)+dtreff*velsou(kiorj,j,kspec,1),
c%OS  1          fnhalf(kiorj,j,kspec,2)+dtreff*velsou(kiorj,j,kspec,2)),
c%OS  1                                                     ivelmid.eq.0)
 522  continue
 5221 continue

c     i>iyh : numindx=1, (2,numixts=+1) or 4
      if (numindx.eq.2 .and. numixts.eq.-1) go to 5231
      iieff=iy_(1)+1-kiorj
      do 523 j=2,jx
        bndmats(j,1,1,2)=0.0
        bndmats(j,1,kleft+1,2)=1.0
        bndmats(j,1,iband,2)=-1.0
c%OS  bndmats(j,1,iband,1)=0.0
        rhspar(1,j,2)=0.0
c%OS  rhspar(1,j,2)=cvmgt(f(iieff,j,kspec,1)+dtreff*
c%OS  1         velsou(iieff,j,kspec,1),flin(dszp5(1),dszp5(1)+dszp5(2),
c%OS  1          fnhalf(iieff,j,kspec,1)+dtreff*velsou(iieff,j,kspec,1),
c%OS  1          fnhalf(iieff,j,kspec,2)+dtreff*velsou(iieff,j,kspec,2)),
c%OS  1                                                     ivelmid.eq.0)
 523  continue
 5231 continue

      if (numindx.eq.1 .or. numindx.eq.4) return

c.......................................................................
c     right boundary: l=l_upper(i) (forward schemes)
 521  continue
      ll=l_upper(kiorj)

c     i<iyh+1 : numindx=(2,numixts=+1) or 3
      if (numindx.eq.2 .and. numixts.eq.-1) go to 5241
      do 524 j=2,jx
c%OS  bndmats(j,ll,1,1)=0.0
        bndmats(j,ll,1,1)=-1.0
        bndmats(j,ll,kleft+1,1)=1.0
        bndmats(j,ll,iband,1)=0.0
        rhspar(ll,j,1)=0.0
c%OS  rhspar(ll,j,1)=cvmgt(f(kiorj,j,kspec,ll)+dtreff*
c%OS  1  velsou(kiorj,j,kspec,ll),flin(dszm5(ll),dszm5(ll)+dszm5(ll-1),
c%OS  1  fnhalf(kiorj,j,kspec,ll-1)+dtreff*velsou(kiorj,j,kspec,ll-1),
c%OS  1  fnhalf(kiorj,j,kspec,ll-2)+dtreff*velsou(kiorj,j,kspec,ll-2)),
c%OS  1                                                     ivelmid.eq.0)
 524  continue
 5241 continue

c     i>iyh : numindx=(2,numixts=-1) or 3
      if (numindx.eq.2 .and. numixts.eq.+1) go to 5251
      iieff=iy_(ll)+1-kiorj
      do 525 j=2,jx
c%OS  bndmats(j,ll,1,2)=0.0
        bndmats(j,ll,1,2)=-1.0
        bndmats(j,ll,kleft+1,2)=1.0
        bndmats(j,ll,iband,2)=0.0
        rhspar(ll,j,2)=0.0
c%OS  rhspar(ll,j,1)=cvmgt(f(iieff,j,kspec,ll)+dtreff*
c%OS  1    velsou(iieff,j,kspec,ll),flin(dszm5(ll),dszm5(ll)+dszm5(ll-1),
c%OS  1    fnhalf(iieff,j,kspec,ll-1)+dtreff*velsou(iieff,j,kspec,ll-1),
c%OS  1    fnhalf(iieff,j,kspec,ll-2)+dtreff*velsou(iieff,j,kspec,ll-2)),
c%OS  1                                                     ivelmid.eq.0)
 525  continue
 5251 continue

      return

c.......................................................................
cl    6. Periodic case (for non-constant y-mesh)
c.......................................................................

 600  continue

c.......................................................................
cl    6.1 Impose fixed value at l=l_upper(i), if different from ls
c     (Assumes l_lower(i)=1=l_mdpln and 
c     lsbtopr(ls-l_upper(i)+2)=lsbtopr(l_upper(i))+1 for all i)
c     Note: There are 4 boundary points at l=l_upper(i) and
c     l=ls-l_upper(i)+2 for i=kiorj and i=iymax+1-kiorj
c     => by reflection, we impose that f is symmetric at theta=pi/2
c     (that is, impose mid-value at iyh: 0.5*(f(iyh_(l))+f(iyh_(l)+1)))
c
c     New try: impose u*cos(theta)*df/ds same and opposite sign, thus impose:
c     df(iyh)/dt+df(iyh+1)/dt=velsou(iyh)+velsou(iyh+1)
c     and    f(iyh+1)=f(iyh), which gives
c     2*f_n+1(iyh)=f_n+1/2(iyh)+f_n+1/2(iyh+1)+dtreff*(velsou(iyh)+velsou(iyh+1))
c.......................................................................
      inewtry=1
c%OS***************************

c     no boundary conditions for passing particles (apart from periodicity)

      if (kiorj .le. itl_(1)) return

      if (inewtry .eq. 1) go to 620

c     special case if only one point: l_lower(i)=l_upper(i)=1
      if (l_upper(kiorj) .eq. 1) then
        do 610 j=2,jx
          bndmats(j,1,1,1)=0.0
          bndmats(j,1,1,2)=0.0
          bndmats(j,1,kleft+1,1)=1.0
          bndmats(j,1,kleft+1,2)=1.0
          bndmats(j,1,iband,1)=0.0
          bndmats(j,1,iband,2)=0.0
          rhspar(1,j,1)=f(kiorj,j,kspec,1)
          rhspar(1,j,2)=f(iy_(1)+1-kiorj,j,kspec,1)
 610    continue
        return
      endif

c     Dependent on scheme for numindx=2 (thus on numixts)
c     iforbak=1 if forward scheme for theta<pi/2 and backward otherwise
c     0 if backward  "     "     "        "  forward    "
      iforbak=(numixts+1)/2
      if (numindx .ne. 2) iforbak=-99

c     point (i=kiorj,l=l_upper(kiorj))
      if (iforbak.eq.1 .or. numindx.eq.3) then
        ll=l_upper(kiorj)
        iief=iy_(ll)+1-kiorj
        lsrow=lsbtopr(ll)
        do 611 j=2,jx
          bndmats(j,lsrow,1,1)=0.0
          bndmats(j,lsrow,kleft+1,1)=1.0
          bndmats(j,lsrow,iband,1)=0.0
          rhspar(lsrow,j,1)=0.5*(f(kiorj,j,kspec,ll)+f(iief,j,kspec,ll))
 611    continue
      endif

c     point (i=iymax+1-kiorj,l=l_upper(kiorj))
      if (iforbak.eq.0 .or. numindx.eq.3) then
        ll=l_upper(kiorj)
        iief=iy_(ll)+1-kiorj
        lsrow=lsbtopr(ll)
        do 612 j=2,jx
          bndmats(j,lsrow,1,2)=0.0
          bndmats(j,lsrow,kleft+1,2)=1.0
          bndmats(j,lsrow,iband,2)=0.0
          rhspar(lsrow,j,2)=0.5*(f(kiorj,j,kspec,ll)+f(iief,j,kspec,ll))
 612    continue
      endif

c     point (i=kiorj,l=ls+2-l_upper(kiorj))
      if (numindx.eq.1 .or. numindx.eq.4 .or. iforbak.eq.0) then
        ll=ls+2-l_upper(kiorj)
        iief=iy_(ll)+1-kiorj
        lsrow=lsbtopr(ll)
        do 613 j=2,jx
          bndmats(j,lsrow,1,1)=0.0
          bndmats(j,lsrow,kleft+1,1)=1.0
          bndmats(j,lsrow,iband,1)=0.0
          rhspar(lsrow,j,1)=0.5*(f(kiorj,j,kspec,ll)+f(iief,j,kspec,ll))
 613    continue
      endif

c     point (i=iymax+1-kiorj,l=ls+2-l_upper(kiorj))
      if (numindx.eq.1 .or. numindx.eq.4 .or. iforbak.eq.1) then
        ll=ls+2-l_upper(kiorj)
        iief=iy_(ll)+1-kiorj
        lsrow=lsbtopr(ll)
        do 614 j=2,jx
          bndmats(j,lsrow,1,2)=0.0
          bndmats(j,lsrow,kleft+1,2)=1.0
          bndmats(j,lsrow,iband,2)=0.0
          rhspar(lsrow,j,2)=0.5*(f(kiorj,j,kspec,ll)+f(iief,j,kspec,ll))
 614    continue
      endif

      return

c.......................................................................
cl    6.2 New try (See above)
c.......................................................................

 620  continue

c     special case if only one point: l_lower(i)=l_upper(i)=1
      if (l_upper(kiorj) .eq. 1) then
        do 621 j=2,jx
          bndmats(j,1,1,1)=0.0
          bndmats(j,1,1,2)=0.0
          bndmats(j,1,kleft+1,1)=1.0
          bndmats(j,1,kleft+1,2)=1.0
          bndmats(j,1,iband,1)=0.0
          bndmats(j,1,iband,2)=0.0
          rhspar(1,j,1)=f(kiorj,j,kspec,1)+
     !      dtreff*velsou(kiorj,j,kspec,1)
          rhspar(1,j,2)=f(iy_(1)+1-kiorj,j,kspec,1)+
     !      dtreff*velsou(iy_(1)+1-kiorj,j,kspec,1)
 621    continue
        return
      endif

c     Dependent on scheme for numindx=2 (thus on numixts)
c     iforbak=1 if forward scheme for theta<pi/2 and backward otherwise
c     0 if backward  "     "     "        "  forward    "
      iforbak=(numixts+1)/2
      if (numindx .ne. 2) iforbak=-99

c     Replace all "four" equations by corresponding values

c     points (i=kiorj,ll=l_upper(kiorj)) and (i=iy_(ll)+1-kiorj,ll=l_upper(kiorj))
      ll=l_upper(kiorj)
      iief=iy_(ll)+1-kiorj
      lsrow=lsbtopr(ll)
      do 622 ipoint=1,2
        do 623 j=2,jx
          bndmats(j,lsrow,1,ipoint)=0.0
          bndmats(j,lsrow,2,ipoint)=0.0
          bndmats(j,lsrow,kleft+2,ipoint)=0.0
          bndmats(j,lsrow,iband,ipoint)=0.0
          bndmats(j,lsrow,kleft+1,ipoint)=2.0
          rhspar(lsrow,j,ipoint)=f(kiorj,j,kspec,ll)+f(iief,j,kspec,ll)+
     !      dtreff*(velsou(kiorj,j,kspec,ll)+velsou(iief,j,kspec,ll))
 623    continue
 622  continue

c     points (kiorj,ls+2-l_upper(kiorj)) and (iy_(ll)+1-kiorj,ls+2-l_upper(kiorj))
      ll=ls+2-l_upper(kiorj)
      iief=iy_(ll)+1-kiorj
      lsrow=lsbtopr(ll)
      do 624 ipoint=1,2
        do 625 j=2,jx
          bndmats(j,lsrow,1,ipoint)=0.0
          bndmats(j,lsrow,1,ipoint)=0.0
          bndmats(j,lsrow,kleft+2,ipoint)=0.0
          bndmats(j,lsrow,iband,ipoint)=0.0
          bndmats(j,lsrow,kleft+1,ipoint)=2.0
          rhspar(lsrow,j,ipoint)=f(kiorj,j,kspec,ll)+f(iief,j,kspec,ll)+
     !      dtreff*(velsou(kiorj,j,kspec,ll)+velsou(iief,j,kspec,ll))
 625    continue
 624  continue

      return

c.......................................................................
cl    6.3 Other boundary conditions for periodic case
c.......................................................................

 630  continue

      return
      end
