c
c
      subroutine coefwti(k)
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     This routine computes w = dy(i,l_)*dd(i,j)/df(i,j) and then sets
c     g(w)=di(i,j,k,l_)=1./w - (1./(exp(w)-1.))
c     Note g(minus infinity)=1.; g(0)=.5; g(plus infinity)=0.
c     This is the Chang-Cooper algorithm (2-d form ala Karney)
c     This routine actually has modified the above by subtracting
c     off the RF contribution to the diffusion (see the coding).
c     This procedure, while not yet
c     justified theoretically, helps to keep the distribution from
c     going negative for strongly driven problems.
c..................................................................

      include 'param.h'
      include 'comm.h'

      include 'advnce.h'
      
c-YuP      call bcast(di(0,1,k,l_),half,(iy+1)*jx) ! could it be error?
      call bcast(di(0,0,k,l_),half,(iy+1)*(jx+2))
      !Note: di(0:iy,0:jx+1,1:ngen,lrors)

      if (chang .ne. "disabled") then
      
        if (chang.eq."noneg") then        
          if (n.eq.1) jchang(k,l_)=jx
          do 100 j=1,jx
            do 1001 i=1,iy
              if (f(i,j,k,l_).le.0.d0) then
                jval=j
                go to 101
              endif
 1001       continue
 100      continue
          jval=jx
 101      continue
          if (jval.ne.jx) then
            if (jchang(k,l_).gt.jval-5) then
              jchang(k,l_)=jval-5
              if (jchang(k,l_).le.0) jchang(k,l_)=1
            endif
          endif
        elseif (chang.eq."enabled") then
          jchang(k,l_)=jx
        endif

       if (jchang(k,l_).le.0) jchang(k,l_)=1  
       !-YuP 101126: added to prevent jchang=0
        
        op=one+em8
        call coefmidt(dff,3)
        do 10 j=1,jchang(k,l_)
          do 2 i=1,iy
            temc1(i)=dy(i,l_)*dd(i,j)*op*df(i,j)/(op*df(i,j)
     1        -dff(i,j))**2
 2        continue

c...............................................................
c     Keep code from blowing up in the 9 loop below.
c..............................................................

          do 8 i=1,iy
            if(abs(temc1(i)).gt.em6) then
              temc3(i)=temc1(i)
            else
              temc3(i)=em6
            endif
 8        continue

          do 4 i=1,iy
            if(temc3(i).lt.sevenhun) then
              temc2(i)=temc3(i)
            else
              temc2(i)=sevenhun
            endif
 4        continue

          do 7 i=1,iy
            if(temc2(i).le.-sevenhun) then
              temc2(i)=-sevenhun
            endif
 7        continue

c..............................................................
c     Evaluate the ratio (modified)
c..............................................................

          do 9 i=1,iy
            di(i,j,k,l_)=1.d0/temc2(i)-(1.d0/(exp(temc2(i))-1.d0))
 9        continue

c.............................................................
c     Correct for errors in 9 above
c.............................................................

          do 21 i=1,iy
            wsub=(three+temc1(i))/(two+temc1(i))/three
            if(temc3(i).eq.em6) then
              di(i,j,k,l_)=wsub
            endif
 21       continue

c..............................................................
c     Limit for large positive or negative temc1 follows
c..............................................................

          do 22 i=1,iy
            if (temc2(i).eq.sevenhun) then
              di(i,j,k,l_)=1.d0/temc1(i)
            elseif (temc2(i).eq.-sevenhun) then
              di(i,j,k,l_)=1.d0+1.d0/temc1(i)
            endif
 22       continue
 10     continue

c...............................................................
c     Now force one-sided differencing if desired.
c...............................................................
       if (jchang(k,l_).le.0) jchang(k,l_)=1  
       !-YuP 101126: added to prevent jchang=0

        if (jchang(k,l_).lt.jx) then
          do 60 j=jchang(k,l_),jx
            do 50 i=1,iy
              if(dd(i,j).ge.0.d0) then
                di(i,j,k,l_)=zero
              else
                di(i,j,k,l_)=one
              endif
 50         continue
 60       continue
        endif
        
      endif

c.......................................................................
c     Ensures correct differentiation at end of intervals
c.......................................................................

      do 3 j=1,jx
        di(0,j,k,l_)=0.d0
        di(iy,j,k,l_)=1.d0
 3    continue
 

      return
      end
