      subroutine tdtransp
      implicit integer (i-n), real*8 (a-h,o-z)
c..............................................................
c     Time advancement splitting step for radial transport.
c..............................................................
      include 'param.h'
      include 'comm.h'
      character*8 nobind
      common/nob/ nobind
      include 'trans.h'
      nobind="disabled"
c.......................................................................
      call dcopy(iyjx2*ngen*lrors,f(0,0,1,1),1,fvn(0,0,1,1),1)
c..............................................................
c     The velocity split for this time step is done and we assume
c     that the results of this split resides in fvn(i,j,k,l) and
c     that this is defined on the velocity or complete mesh. Since
c     the tranport equation is not advanced in the vicinity of the
c     pass/trapped boundary we first interpolate onto a new mesh
c     that that does not include these mesh points, and this is
c     accomplished in a density conserving manner.
c..............................................................
      if (nonadi .eq. 6) then
        call tdtrvtor2(fvn(0,0,1,1),frn(0,0,1,1),vpint,vpint_,1)
      else
        call tdtrvtor(fvn,frn)
      endif
c..............................................................
c     Generate advection coefficients - then determine the
c     Chang-Cooper weights (radial direction).
c     Presently only set up for 1 general species,
c       except for pinch.eq."disabled" (adv()=0.).
c..............................................................
      kprofile=kelecm
      if (niong.ne.0) kprofile=kionm(1)
      if (colmodl.eq.0) kprofile=1
      call tdtravct(frn,kprofile,1)  !No-op if pinch="disabled"
      call tdtrwtl
c..............................................................
c     Keep a copy for accuracy check later...
c..............................................................
      call dcopy(iyjx2*ngen*(lrors+1),frn,1,frn_2,1)
c..............................................................
c     Loop over species index.
c..............................................................
      do 100 k=1,ngen
c..............................................................
c     Loop over momentum
c..............................................................
        do 200 ic=1,iytr(lrors)/2
          if (l_lower(ic).eq.lrors) go to 200
          i_=iytr(lrors)+1-ic
          do 201 ii=1,2
            i=ic
            if (ii.eq.2) i=i_
            lu=l_lower(i)
c..............................................................
c     Initialize the calculation (zero flux) at the innermost
c     radial mesh point seen by a particle with a given adiabatic
c     invariant if meshy="fixed_mu" or at l=1 if meshy="fixed_y".
c     In the event that there is a binding equation ((lpt(i).ne.
c     lrors) meaning a particle with a given "i" is passing at the center
c     (l=1) and at some point (lpt(i_)) becomes trapped) we will
c     also initialize at l=lrors. There will then be two initial
c     sweeps, instead of just one. One sweep will be up in "l" the
c     other down, and both terminate at lpt(i).
c..............................................................
            do 222 j=1,jx
              fg_(j,ii,lu)=delr(i,lu)/betr(i,lu)
              eg_(j,ii,lu)=alpr(i,lu)/betr(i,lu)
              frn(idx(i,lrors),j,k,lrors)=vptb_(idx(i,lrors),
     1          lrindx(lrors))/zmaxpsi(lrindx(lrors))
     +          *frn(idx(i,lrors),j,k,lrors)
              eg_(j,ii,lrors)=0.
              fg_(j,ii,lrors)=frn(idx(i,lrors),j,k,lrors)
 222        continue
c..................................................................
c     Now complete the initial sweep over 1 (or 2) regions.
c..................................................................
            do 240 l=l_lower(i)+1,lpt(i)-1
              do 230 j=1,jx
                eg_(j,ii,l)=alpr(i,l)/(betr(i,l)-gamr(i,l)
     1            *eg_(j,ii,l-1))
                fg_(j,ii,l)=(delr(i,l)+gamr(i,l)*fg_(j,ii,l-1))
     1            /(betr(i,l)-gamr(i,l)*eg_(j,ii,l-1))
 230          continue
 240        continue
            if (lpt(i).eq.lrors) go to 256
c..................................................................
c     second region, if necessary...
c..................................................................
            do 250 l=lrors-1,lpt(i)+1,-1
              do 255 j=1,jx
                eg_(j,ii,l)=gamr(i,l)/(betr(i,l)-alpr(i,l)
     1            *eg_(j,ii,l+1))
                fg_(j,ii,l)=(delr(i,l)+alpr(i,l)
     1            *fg_(j,ii,l+1))
     1            /(betr(i,l)-alpr(i,l)*eg_(j,ii,l+1))
 255          continue
 250        continue
 256        continue
 201      continue
c..................................................................
c     The binding equation (at lpt(i)) is next..
c..................................................................
          i=ic
          if (lpt(i).eq.lrors) go to 260
          lv=lpt(i)
          if (nobind.eq."disabled") then

            do 265 j=1,jx
              tam1(j)=-(gamr(i,lv)*eg_(j,1,lv-1)+gamr(i_,lv)
     1          *eg_(j,2,lv-1))*.5
     1          +(betr(i,lv)+betr(i_,lv))*.5 - alpr(i,lv)*eg_(j,1,lv+1)
              tam2(j)=(gamr(i,lv)*fg_(j,1,lv-1)+gamr(i_,lv)
     1          *fg_(j,2,lv-1))*.5
     1          +alpr(i,lv)*fg_(j,1,lv+1)+delr(i,lv)
              frn(idx(i,lpt(i)),j,k,lpt(i))=tam2(j)/tam1(j)
              frn(idx(i_,lpt(i_)),j,k,lpt(i_))=tam2(j)/tam1(j)
 265        continue
          else
            do 266 j=1,jx
              tam1(j)=-(gamr(i,lv)*eg_(j,1,lv-1))
     1          +(betr(i,lv)) - alpr(i,lv)*eg_(j,1,lv+1)
              tam2(j)=(gamr(i,lv)*fg_(j,1,lv-1))
     1          +alpr(i,lv)*fg_(j,1,lv+1)+delr(i,lv)
              frn(idx(i,lpt(i)),j,k,lpt(i))=tam2(j)/tam1(j)
              tam3(j)=-(gamr(i_,lv)*eg_(j,2,lv-1))
     1          +(betr(i_,lv)) - alpr(i_,lv)*eg_(j,2,lv+1)
              tam4(j)=(gamr(i_,lv)*fg_(j,2,lv-1))
     1          +alpr(i_,lv)*fg_(j,2,lv+1)+delr(i_,lv)
              frn(idx(i_,lpt(i_)),j,k,lpt(i_))=tam4(j)/tam3(j)
 266        continue
          endif
 260      continue
c..................................................................
c     solve for the new distribution...
c..................................................................
          do 285 ij=1,2
            i=ic
            if (ij.eq.2) i=i_
            do 280 l=lpt(i)-1,l_lower(i),-1
              do 270 j=1,jx
                frn(idx(i,l),j,k,l)=eg_(j,ij,l)*frn(idx(i,l+1),j,k,l+1)
     1            +fg_(j,ij,l)
 270          continue
 280        continue
            do 295 l=lpt(i)+1,lrors
              do 290 j=1,jx
                frn(idx(i,l),j,k,l)=eg_(j,ij,l)*frn(idx(i,l-1),j,k,l-1)
     1            +fg_(j,ij,l)
c              write(*,'(a,3i4,2e12.4)') 'transp=== l,j,i,d_r=', 
c     ~         l,j,i,d_r(i,j,k,l),d_rr(i,j,k,l)
              
 290          continue
 295        continue
c         pause
         
            do 310 l=l_lower(i),lrors
              ii=idx(i,l)
              do 312 j=1,jx
                frn(ii,j,k,l)=frn(ii,j,k,l)/vptb_(ii,lrindx(l))*
     *            zmaxpsi(lrindx(l))
 312          continue
 310        continue
 
 285      continue

 200    continue
 100  continue

      call tdtrchk

      if (nobind.eq."enabled") then
        call tdtrsym
      endif
c.................................................................
c     Redefine frn at v=0 so it is unique.
c.................................................................
      if (nonadi .ne. 3) then
        do 400 k=1,ngen
          do 410 l=1,lrors
            zs=0.
            zt=0.
            do 420 i=1,iytr(l)
              zs=zs+vpint_(idx(i,l),lrindx(l))
              zt=zt+vpint_(idx(i,l),lrindx(l))*frn(idx(i,l),1,k,l)
 420        continue
c              write(*,*) 'transp=== l,zs,zt=', l,zs,zt
            do 430 i=1,iytr(l)
              frn(idx(i,l),1,k,l)=zt/zs
 430        continue
 410      continue
 400    continue
      endif
c......................................................................
c     The distribution frn is currently defined on the transport 
c     mesh - interpolate onto the velocity mesh, returning it in frn. 
c......................................................................
      if (nonadi .eq. 6) then
        call tdtrrtov2(frn(0,0,1,1),frn(0,0,1,1),vpint,vpint_,1)
      else
        call tdtrrtov(frn)
      endif
      return
      end
