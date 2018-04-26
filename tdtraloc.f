c
c
      subroutine tdtraloc
      implicit integer (i-n), real*8 (a-h,o-z)

c..............................................................
c     Allocates arrays used by transport routines.
c..............................................................

      include 'param.h'
      include 'comm.h'
c.......................................................................

      lnyxgrz=(iyp1+1)*(jxp1+1)*ngen*(lrz+1)
      lnjxng=jx*ngen*lrz
      lnjxg18=18*lnjxng
      lnyrs=iy*lrz
      lnyrmx=iy*lrzmax
      lnyrmxp=iy*(lrzmax+1)
      lnyprsp=(iy+1)*(lrz+1)
      lnng=ngen*lrzmax
      lnrz=lrz
      lnj2=jx*2*lrz

      lndum=8*lnyxgrz+3*lnjxng+lnjxg18+lnyrs+lnyrmx+lnyrmxp+
     +  2*lnyprsp+lnng+lnrz+2*lnj2

c..................................................................
c     A check on allocations is sucessful entering then exiting
c     the subroutine.
c..................................................................
      write(*,*)'tdtraloc:  Entering tdtraloc'

      allocate(frn_2(0:iy+1,0:jx+1,ngen,0:lrz),STAT=istat)
      call bcast(frn_2,zero,SIZE(frn_2))
      allocate(frn_1(0:iy+1,0:jx+1,ngen,0:lrz),STAT=istat)
      call bcast(frn_1,zero,SIZE(frn_1))
      allocate(frn(0:iy+1,0:jx+1,ngen,0:lrz),STAT=istat)
      call bcast(frn,zero,SIZE(frn))
      allocate(fvn_1(0:iy+1,0:jx+1,ngen,0:lrz),STAT=istat)
      call bcast(fvn_1,zero,SIZE(fvn_1))
      allocate(fvn(0:iy+1,0:jx+1,ngen,0:lrz),STAT=istat)
      call bcast(fvn,zero,SIZE(fvn))
      
      
      allocate(dl(0:iyp1,0:jxp1,ngen,0:lrz),STAT=istat)
      call bcast(dl,zero,SIZE(dl))
      allocate(d_rr(0:iyp1,0:jxp1,ngen,0:lrz),STAT=istat)
      call bcast(d_rr,zero,SIZE(d_rr))
      allocate(d_r(0:iyp1,0:jxp1,ngen,0:lrz),STAT=istat)
      call bcast(d_r,zero,SIZE(d_r))
      
      allocate(f_lm(jx,ngen,lrz),STAT=istat)
      call bcast(f_lm,zero,SIZE(f_lm))
      allocate(f_lp(jx,ngen,lrz),STAT=istat)
      call bcast(f_lp,zero,SIZE(f_lp))
      allocate(f_up(jx,ngen,lrz),STAT=istat)
      call bcast(f_up,zero,SIZE(f_up))
      allocate(f_vtor(jx,ngen,18,lrz),STAT=istat)
      call bcast(f_vtor,zero,SIZE(f_vtor))
      allocate(cynt2_(iy,lrz),STAT=istat)
      call bcast(cynt2_,zero,SIZE(cynt2_))
      allocate(vpint_(iy,lrzmax),STAT=istat)
      call bcast(vpint_,zero,SIZE(vpint_))
      allocate(vptb_(iy,0:lrzmax),STAT=istat)
      call bcast(vptb_,zero,SIZE(vptb_))
      allocate(cosovb(0:iy,0:lrz),STAT=istat)
      call bcast(cosovb,zero,SIZE(cosovb))
      allocate(bovcos(0:iy,0:lrz),STAT=istat)
      call bcast(bovcos,zero,SIZE(bovcos))
      allocate(adv(ngen,lrzmax),STAT=istat)
      call bcast(adv,zero,SIZE(adv))
      allocate(dentarget(lrz),STAT=istat)
      call bcast(dentarget,zero,SIZE(dentarget))
      allocate(eg_(jx,2,lrz),STAT=istat)
      call bcast(eg_,zero,SIZE(eg_))
      allocate(fg_(jx,2,lrz),STAT=istat)
      call bcast(fg_,zero,SIZE(fg_))

      write(*,*)'tdtraloc:  Leaving tdtraloc'

      return
      end
