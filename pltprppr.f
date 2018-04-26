c
c
      subroutine pltprppr
      implicit integer (i-n), real*8 (a-h,o-z)
c...  
cmnt  This routine plots the parallel distribution function.
c...
c
c     Modified from Graflib to pgplot calls by Yuri Petrov, 090727,
c     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
c
c     Need to re-work this a bit, for plot of pitch angle averaged
c     reduced distributions (knockon=.ne."disabled" case).
c  
      include 'param.h'
      include 'comm.h'

      REAL RILIN !-> For PGPLOT (text output positioning)

      character*8 target
      if (noplots.eq."enabled1") return

c     Return, if using a theta average distribution
c     (which may conflict with call fle_fsa or fle_pol, below).
      if (flemodel.ne."pol"  .and. flemodel.ne."fsa") return

      do 20 k=1,ngen
        if (tandem.eq."enabled" .and. k.eq.kionn) then
          xll=-xlwr
          xlu=xlwr
          xpl=0.
          xpu=xlwr
          target="ionmesh"
        else
          xll=-xmax
          xlu=xmax
          xpl=0.
          xpu=xmax
          if (xprpmax.ne.1.) xpu=xprpmax
          target="mainmesh"
        endif
c       Obtain equatorial plane pitch-avg'd distribution for ko model
c       case.  (Else setup interfere's with soureko).
        if (knockon.ne."disabled") then
c           call fle("setup",0)
c           call fle("calc",1)
        elseif (lrz.eq.1) then
           call fle_pol("setup",0)
           call fle_pol("calc",1)
        else
           call fle_fsa("setup")
           call fle_fsa("calc")
        endif

c     plot (fll(+xpar)-fll(-xpar)), and xpar*(fll(+xpar)-fll(-xpar))
        jpxyh=(jfl+1)/2
        jpxyhm=jpxyh-1
        tem1(1)=0.0
        tem2(1)=0.0
        do 30  jp=1,jpxyhm
          tem1(jp)=fl(jpxyh+jp-1)-fl(jpxyh-jp)
          tem2(jp)=xlm(jpxyh+jp-1)*tem1(jp)
 30     continue
        call aminmx(tem1,1,jpxyhm,1,fmin1,fmax1,kmin,kmax)
        call aminmx(tem2,1,jpxyhm,1,fmin2,fmax2,kmin,kmax)
c990131        fmin=amin1(fmin1,fmin2)
c990131        fmax=amax1(fmax1,fmax2)
        fmin=min(fmin1,fmin2)
        fmax=max(fmax1,fmax2)
        
        call GXGLFR(0) ! new page(s)
        call GSVP2D(.2,.8,.25,.95) ! (XLEFT, XRIGHT, YBOT, YTOP)
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlin$",0.d0,xlu,fmin,fmax)
        call GPCV2D(xlm(jpxyh),tem1,jpxyhm)
        call GPCV2D(xlm(jpxyh),tem2,jpxyhm)
        
        write(t_,10011) k
10011 format("asymmetric cmpt of f_par, and xpar*cmpt, species:",1x,i5)
        RILIN=5.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)
        
        write(t_,10012)
10012 format("(f_par normed so int(-1,+1)=equatorial ne)")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)
       
        call GXGLFR(0) ! new page(s)
        call aminmx(fl(1),1,2*jpxyhm,1,fmin,fmax,kmin,kmax)
        fmin=1.d-08*fmax
        call GSVP2D(.2,.8,.25,.95)
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlog$",xll,xlu,fmin,fmax)
        do 10 jj=1,2*jpxyhm
          if (fl(jj) .lt. fmin ) fl(jj)=fmin
 10     continue
      !-YuP:   call GSVTCL("on$")
        call GPCV2D(xlm,fl,2*jpxyhm)
      !-YuP:   call GSVTCL("off$")
        
        write(t_,10013) k
10013 format("parallel distribution function for species:",1x,i5)
        RILIN=3.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)
        
        write(t_,10014) 
10014 format("(normed so int(-1,+1)=equatorial ne)")        
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)
        
        write(t_,10020)
10020 format( "(log plot)")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)
        
        rr=rpcon(lr_) !rovera(lr_)*radmin  ! YuP[03-2016] changed to rpcon
        write(t_,10030) n,timet
10030 format("time step (n) is",i5,5x,"time=",e14.6," secs")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)

        write(t_,10031) rovera(lr_),rr
10031 format("r/a=",e14.6,5x,"radial position (R) =",e14.6," cm")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_)
        
        CALL PGSCH(1.0) ! recover default 1.0 fontsize
        
 20   continue

      return
      end
