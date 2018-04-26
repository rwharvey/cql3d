c
c
      subroutine wpmshchk
      implicit integer (i-n), real*8 (a-h,o-z)

c.......................................................................
c     Check some mesh related values. e.g. sin(th)dth (s)
c.......................................................................

      include 'param.h'
      include 'comm.h'

      return


c.......................................................................
cl    1. compare sin(th)dth with Bdmu
c
c%OS  do 100 i=2,iymax/2
c%OS  write(6,9100) i,(sinn(i,l+1)*dy(i,l+1)/(sinn(i,l)*dy(i,l)),
c%OS  !         psis(l+1)/psis(l)*coss(i,l)/coss(i,l+1),l=1,l_upper(i)-1)
c%OS  100  continue
c%OS  9100 format(i3,1p10e12.4,/,(3x,1p10e12.4))

c%OS  write(6,9100)
c%OS  do 100 i=1,iymax/2
c%OS  write(6,9101) i,(sinn(i,l)/sqrt(psis(l)),
c%OS  !        coss(i,l)*dy(i,l)/sqrt(psis(l)),l=1,l_upper(i))
c%OS  100  continue
      print *," for each i, coss(i,l)*cynt2(i,l)/psis(l),l=1,l_upper(i)"
      do 101 i=1,iymax/2
        write(6,9101) i,(cynt2(i,l)/psis(l)*coss(i,l)
     !    ,l=1,l_upper(i))
        if (l_upper(i) .lt. lrors) then
          write(6,9101) i,(cynt2(i,l)/psis(l)*coss(i,l)
     !      ,l=lrors-l_upper(i)+2,lrors+1)
        endif
 101  continue
      print *,' i coss(i,l) sin(thp)*dthp*bbpsi*cos(thp)/sin(th)/dth  ',
     !  'sqrt(1-bbpsi*sin(thp)**2)'
      do 1011 i=1,iymax/2
        write(6,9102) i,(coss(i,l),cynt2(i,1)*psis(l)*coss(i,1)
     +    /cynt2(i,l)
     !    ,sqrt(1.-psis(l)*sinn(i,1)**2),l=1,l_upper(i))
 1011 continue
      print *,'  l',' sin*dth','  (dthp)'
      do 102 l=1,ls
        zsum0=0.
        zsums=0.
        do 103 i=1,iyh_(l)
          zsum0=zsum0+cynt2(i,1)*psis(l)*coss(i,1)/coss(i,l)
          zsums=zsums+cynt2(i,l)
          if (l .eq. 2) then
            print *,cynt2(i,l),cynt2(i,1)*psis(l)*coss(i,1)/coss(i,l)
          endif
 103    continue
        write(6,'(i3,1p2e12.3)') l,zsums,zsum0
 102  continue

      if (ls.gt.10 .or. cqlpmod.ne."enabled" .or. transp.ne."enabled") 
     !  go to 999

      print *,' dy adjusted and dy at iyh'
      zsum0=0.0
      do 1031 i=1,iyh_(1)
        zsum0=zsum0+cynt2(i,1)
 1031 continue
      do 1032 l=2,ls/2+1
        zsum1=0.0
        do 1033 i=1,iyh_(l)-1
          zsum1=zsum1+cynt2(i,l)
 1033   continue
        zzdyh=(zsum0-zsum1)/twopi/sinn(i,l)
        print *,zzdyh,dy(iyh_(l),l),(zzdyh-dy(iyh_(l),l))/dy(iyh_(l),l)
     !    *100.
 1032 continue
      l=1
      write(6,'(/," y(i+1)-y(i)(l),dyp5(l); l=",i3,/,(1p10e12.4))') 
     !  l,(y(i+1,l)-y(i,l),dyp5(i,l),i=1,iyh_(l))
      do 104 l=2,ls
        do 105 i=2,iyh_(l)
          ii=iy_(l)-i+1
          dyp5(i,l)=dyp5(i,1)*psis(l)*sinn(i,1)*coss(i,1)/coss(i,l)/
     !      sinn(i,l)
          dym5(i,l)=dym5(i,1)*psis(l)*sinn(i,1)*coss(i,1)/coss(i,l)/
     !      sinn(i,l)
          dy(i,l)=0.5*(dym5(i,l)+dyp5(i,l))
          cynt2(i,l)=dy(i,l)*twopi*sinn(i,l)
          dyp5(ii,l)=dym5(i,l)
          dym5(ii,l)=dyp5(i,l)
          dy(ii,l)=dy(i,l)
          cynt2(ii,l)=cynt2(i,l)
 105    continue
        dyp5(1,l)=dyp5(1,1)*psis(l)
        dym5(1,l)=0.0
        dy(1,l)=0.5*dyp5(1,l)
        cynt2(1,l)=psis(l)*cynt2(1,1)
        dym5(iy,l)=dym5(iy,1)*psis(l)
        dyp5(iy,l)=0.0
        dy(iy,l)=0.5*dyp5(iy,l)
        cynt2(iy,l)=psis(l)*cynt2(iy,1)
        write(6,'(/," y(i+1)-y(i)(l),dyp5(l); l=",i3,/,(1p10e12.4))') 
     !    l,(y(i+1,l)-y(i,l),dyp5(i,l),i=1,iyh_(l))
 104  continue
      print *,' '
      print *,'  l',' sum(2*pi*sin*dth)(l)',' sum (dth0) with new dy(l)'
      do 1021 l=1,ls
        zsum0=0.
        zsums=0.
        do 1022 i=1,iyh_(l)
          zsum0=zsum0+cynt2(i,1)*psis(l)*coss(i,1)/coss(i,l)
          zsums=zsums+cynt2(i,l)
          if (l .eq. 2 .and. i.le.iyh_(2)) then
            print *,cynt2(i,l),cynt2(i,1)*psis(l)*coss(i,1)/coss(i,l)
          endif
 1022   continue
        write(6,'(i3,1p2e12.3)') l,zsums,zsum0
 1021 continue

 9100 format(" i ","[sin/sqrt(bbpsi);cos*dth/sqrt(ps)], l=1,l_upper(i)")
 9101 format(i3,1p10e12.4,/,(/,3x,1p10e12.4))
 9102 format(i3,1p9e12.4,/,(3x,1p9e12.4))

      if (ls .le. 10) stop 'wpmshchk: Need to set ls>10'

 999  return
      end
