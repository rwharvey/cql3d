! YuP[2019-10-02] Made some modifications to avoid log10(0.0) 
!  and to avoid plots with near-zero values.
c
      subroutine pltrun
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'
      dimension yg(nonch),xg(nonch)

c     PGPLOT REAL Variables:
      REAL*4 RPG1,RPG2
      REAL*4 RNONCHA1(nonch),RNONCHA2(nonch)
      REAL*4 RBOUND
      REAL*4 :: R4P2=.2,R4P8=.8,R4P5=.5,R4MP1=-.1,R40=0.
      REAL*4 :: R4P6=.6,R4P9=.9,R46=6.,R47=7.,R4P65=.65

c-----------------
c   This subroutine writes and plots runaway population and current
c------------------

CMPIINSERT_INCLUDE

CMPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      if (pltra.eq."disabled") return

      iounit=35
      open(unit=iounit,file='runaway.out',status='unknown')

cBH070405      if (nch(l_).gt.noncha .or. nch(l_).gt.500)
      if (nch(1).gt.nonch) stop 'check dimensions in pltrun'

      write (iounit,20001)
20001 format(3x,'time',12x,'runaway x'//)
      do ll=1,lrz
        write(iounit,20010) lrindx(ll)
20010   format(//2x,'flux surface ',i3)
        write(iounit,20015)
20015   format(//19x,'denra:'/)
cBH070405        do ntm=1,nch(l_)
        do ntm=1,nch(l_)
          write(iounit,20020) ptime(ntm,ll),pdenra(ntm,ll)
20020     format(2x,1pe12.4,2x,1pe18.10)
        enddo
        write(iounit,20016)
20016   format(//19x,'curra:'/)
cBH070405        do ntm=1,nch(l_)
        do ntm=1,nch(ll)
          write(iounit,20020) ptime(ntm,ll),pcurra(ntm,ll)
        enddo
        write(iounit,20017)
20017   format(//19x,'frac. ra-dens.:'/)
cBH070405        do ntm=1,nch(l_)
        do ntm=1,nch(ll)
          write(iounit,20020) ptime(ntm,ll),pfdenra(ntm,ll)
        enddo
        write(iounit,20018)
20018   format(//19x,'frac. ra-cur.:'/)
cBH070405        do ntm=1,nch(l_)
        do ntm=1,nch(ll)
          write(iounit,20020) ptime(ntm,ll),pfcurra(ntm,ll)
        enddo
        write(iounit,20013)
20013   format(//19x,'ucrit/c.:'/)
cBH070405        do ntm=1,nch(l_)
        do ntm=1,nch(ll)
          write(iounit,20020) ptime(ntm,ll),pucrit(ntm,ll)/cnorm
        enddo
        write(iounit,20014)
20014   format(//19x,'e/ecrit0.:'/)
cBH070405        do ntm=1,nch(l_)
        do ntm=1,nch(ll)
          write(iounit,20020) ptime(ntm,ll),peoe0(ntm,ll)
        enddo
      enddo

      write(iounit,20050)
20050 format(//2x,'Pitch Angle Averaged Distribution At Time Slices:'//)

      nframep=min(iplot,5)
      do ns=1,nframep
        write(iounit,20060) tplot(ns)
20060   format(///2x,'Time=', f9.4,'secs.'/)
        do ll=1,lrz
          write(iounit,20070) ll
20070     format(2x,i3,'th flux surface'/,
     1        5x,'GAMMA',15x,'F(GAMMA)'/)
          do j=1,jx
            write(iounit,20080) gamma(j),f_aveth(j,1,ll,ns)
20080       format(1x,e16.8,2x,e20.10)
          enddo
        enddo
      enddo

      close(unit=iounit)

c-----------------------------------------------------------------------
c     Plotting
c-----------------------------------------------------------------------

      if (noplots.eq."enabled1") return

      dgts=1.e-8

      do ll=1,lrz
      
        call tdnflxs(ll)
        rr=rpcon(lr_) !rovera(lr_)*radmin  ! YuP[03-2016] changed to rpcon
cBH070405        do nt=1,nch(l_)
        do nt=1,nch(ll)
          xg(nt)=ptime(nt,1)
          yg(nt)=abs(pdenra(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          if(RNONCHA2(nt).gt.0.)then
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
          else
          RNONCHA2(nt)=-100. !which means log10(1d-100), a very small value
          endif
        enddo
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
        iskip=0 ! to be changed to iskip=1 if denra=0
        if(YMAX.lt.1d-100) iskip=1 ! case of denra~0
        ymax=10.*ymax !1.03*ymax
        ymin=ymax/1.e4

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)
        
c      Following plot is with pdenra,... with 1st dimension up to nch(l_).
        CALL PGPAGE ! -1-  denra and curra
        
        if(iskip.eq.0)then
        CALL PGSVP(R4P2,R4P8,R4P6,R4P9)
        IF (RPG1.lt.RPG2) THEN
           CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
           CALL PGBOX('BCNST',R40,0,'BCNSTL',R40,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
           CALL PGLAB('time (secs)','RA density (cm\u-3\d)',
     +       'Runaway Density and Current vs. Time')
        ENDIF
        endif ! iskip

cBH070405       do nt=1,nch(l_)
        do nt=1,nch(ll)
          xg(nt)=ptime(nt,1)
          yg(nt)=abs(pcurra(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          if(RNONCHA2(nt).gt.0.)then
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
          else
          RNONCHA2(nt)=-100. !which means log10(1d-100), a very small value
          endif
        enddo
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
        iskip=0 ! to be changed to iskip=1 if curra=0
        if(YMAX.lt.1d-100) iskip=1 ! case of curra~0
        ymax=10.*ymax !1.03*ymax
        ymin=ymax/1.e4

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)

        if(iskip.eq.0)then
        CALL PGSVP(R4P2,R4P8,R4P2,R4P5)
        IF (RPG1.lt.RPG2) THEN
           CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
           CALL PGBOX('BCNST',R40,0,'BCNSTL',R40,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
           CALL PGLAB('time (secs)','RA curr den (A/cm\u2\d)',' ')
        ENDIF
        endif ! iskip
        
        write(t_,10010) n,timet
        CALL PGMTXT('B',R46,R4MP1,R40,t_)
        write(t_,10011) rovera(lr_),ll,rr
        CALL PGMTXT('B',R47,R4MP1,R40,t_)

c-----------------------------------------------------------------------
cBH070405       do nt=1,nch(l_)
        do nt=1,nch(ll)
          xg(nt)=ptime(nt,1)
          yg(nt)=abs(pfdenra(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          if(RNONCHA2(nt).gt.0.)then
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
          else
          RNONCHA2(nt)=-100. !which means log10(1d-100), a very small value
          endif
        enddo
c        write(*,*) 'ptime(ijk,1),pfdenra(ijk,ll),ijk=1,5',
c     +              (ptime(ijk,1),pfdenra(ijk,ll),ijk=1,5)
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
        iskip=0 ! to be changed to iskip=1 if fdenra=0
        if(YMAX.lt.1d-100) iskip=1 ! case of fdenra~0
        ymax=10.*ymax !1.03*ymax
        ymin=ymax/1.e4

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)

        CALL PGPAGE ! -2-  fdenra and fcurra
        
        if(iskip.eq.0)then
        CALL PGSVP(R4P2,R4P8,R4P6,R4P9)
        IF (RPG1.lt.RPG2) THEN
           CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
           CALL PGBOX('BCNST',R40,0,'BCNSTL',R40,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
           CALL PGLAB('time (secs)','Fraction RA density',
     +       'Runaway Fraction of Density and Current vs. Time')
        ENDIF
        endif ! iskip

        do nt=1,nch(l_)
          xg(nt)=ptime(nt,1)
          yg(nt)=abs(pfcurra(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          if(RNONCHA2(nt).gt.0.)then
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
          else
          RNONCHA2(nt)=-100. !which means log10(1d-100), a very small value
          endif
         enddo
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
        iskip=0 ! to be changed to iskip=1 if fcurra=0
        if(YMAX.lt.1d-100) iskip=1 ! case of fcurra~0
        ymax=10.*ymax !1.03*ymax
        ymin=ymax/1.e4

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)

        if(iskip.eq.0)then
        CALL PGSVP(R4P2,R4P8,R4P2,R4P5)
        IF (RPG1.lt.RPG2) THEN
           CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
           CALL PGBOX('BCNST',R40,0,'BCNSTL',R40,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
           CALL PGLAB('time (secs)','Fraction RA curr den',' ')
        ENDIF
        endif ! iskip
        
        write(t_,10010) n,timet
        CALL PGMTXT('B',R46,R4MP1,R40,t_)
        write(t_,10011) rovera(lr_),ll,rr
        CALL PGMTXT('B',R47,R4MP1,R40,t_)

c-----------------------------------------------------------------------
c     Plots of ucrit and e/e_dreicer, versus time:

cBH070405       do nt=1,nch(l_)
        do nt=1,nch(ll)
           xg(nt)=ptime(nt,1)
           yg(nt)=abs(pucrit(nt,ll))
           RNONCHA1(nt)=RBOUND(XG(nt))
           RNONCHA2(nt)=RBOUND(YG(nt))
          if(RNONCHA2(nt).gt.0.)then
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
          else
          RNONCHA2(nt)=-100. !which means log10(1d-100), a very small value
          endif
        enddo
      
        if (pltlim.eq.'u/c') then
         do nt=1,nch(l_)
cYuP            yg(nt)=yg(nt)*cnormi !cnormi=0 when relativ.eq."disabled"
          yg(nt)=yg(nt)/cnorm ! YuP[07-2016]
          RNONCHA2(nt)=RBOUND(YG(nt))
          if(RNONCHA2(nt).gt.0.)then
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
          else
          RNONCHA2(nt)=-100. !which means log10(1d-100), a very small value
          endif
         enddo
        elseif (pltlim.eq.'energy') then
cBH070405         do nt=1,nch(l_)
         do nt=1,nch(ll)
            yg2=yg(nt)*yg(nt)
            if(yg2*cnorm2i.lt.1.e-8 .or. relativ.eq."disabled") then
               g1=.5*yg2/cnorm2
            else
               g1=sqrt(1.+yg2*cnorm2i)-1.
            endif
            yg(nt)=g1*restmkev
          RNONCHA2(nt)=RBOUND(YG(nt))
          if(RNONCHA2(nt).gt.0.)then
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
          else
          RNONCHA2(nt)=-100. !which means log10(1d-100), a very small value
          endif
         enddo
        endif
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
        iskip=0 ! to be changed to iskip=1 if ucrit=0
        if(YMAX.lt.1d-100) iskip=1 ! case of ucrit~0
        ymax=10.*ymax !1.03*ymax
        ymin=ymax/1.e4
        !ymin=0.97*ymin
        !ymax=1.e3*ymin

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)


        if (pltlim.eq.'u/c') then
         write(t_,1020)
 1020    format("Critical runaway vel/c")
        elseif (pltlim.eq.'energy') then
         write(t_,1021)
 1021    format("Critical runaway energy (keV)")
        else
         write(t_,1022)
 1022    format("Critical runaway vel/vnorm")
        endif

        CALL PGPAGE ! -3-  ucrit and eoed
        
        if(iskip.eq.0)then
        CALL PGSVP(R4P2,R4P8,R4P6,R4P9)
        IF (RPG1.lt.RPG2) THEN
           CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
           CALL PGBOX('BCNST',R40,0,'BCNSTL',R40,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
           CALL PGLAB('time (secs)',t_,
     +       'Critical vel(energy) and E/E\dDreicer\u vs. Time')
        ENDIF
        endif ! iskip

cBH070405      do nt=1,nch(l_)
        do nt=1,nch(ll)
          xg(nt)=ptime(nt,1)
          yg(nt)=abs(peoed(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          if(RNONCHA2(nt).gt.0.)then
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
          else
          RNONCHA2(nt)=-100. !which means log10(1d-100), a very small value
          endif
        enddo
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
        iskip=0 ! to be changed to iskip=1 if eoed=0
        if(YMAX.lt.1d-100) iskip=1 ! case of eoed~0
        ymax=10.*ymax !1.03*ymax
        ymin=ymax/1.e3

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)

        write(t_,1023)
 1023   format("E-field/E_Dreicer")
        write(t_,1024)
 1024   format("time(sec)")

        if(iskip.eq.0)then
        CALL PGSVP(R4P2,R4P8,R4P2,R4P5)
        IF (RPG1.lt.RPG2) THEN
           CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
           CALL PGBOX('BCNST',R40,0,'BCNSTL',R40,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
           CALL PGLAB('time (secs)','E-field/E\dDreicer\u',' ')
        ENDIF
        endif ! iskip

        write(t_,10010) n,timet
        CALL PGMTXT('B',R46,R4MP1,R40,t_)
        write(t_,10011) rovera(lr_),ll,rr
        CALL PGMTXT('B',R47,R4MP1,R40,t_)

c-----------------------------------------------------------------------


c     Plots of e/e0 and KO source, versus time:

cBH070405      do nt=1,nch(l_)
        do nt=1,nch(ll)
          xg(nt)=ptime(nt,1)
          yg(nt)=abs(peoe0(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          if(RNONCHA2(nt).gt.0.)then
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
          else
          RNONCHA2(nt)=-100. !which means log10(1d-100), a very small value
          endif
        enddo
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
        iskip=0 ! to be changed to iskip=1 if eoe0=0
        if(YMAX.lt.1d-100) iskip=1 ! case of eoe0~0
        ymax=10.*ymax !1.03*ymax
        ymin=ymax/1.e3

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)

        CALL PGPAGE ! -4-  eoe0 and src (KO)
        
        if(iskip.eq.0)then
        CALL PGSVP(R4P2,R4P8,R4P6,R4P9)
        IF (RPG1.lt.RPG2) THEN
           CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
           CALL PGBOX('BCNST',R40,0,'BCNSTL',R40,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
           CALL PGLAB('time (secs)','E-field/E0',
     +       'E/(Critical E0) and KO Source vs. Time')
        ENDIF
        endif ! iskip

cBH070405      do nt=1,nch(l_)
        do nt=1,nch(ll)
          xg(nt)=ptime(nt,1)
          yg(nt)=abs(psrc(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          if(RNONCHA2(nt).gt.0.)then
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
          else
          RNONCHA2(nt)=-100. !which means log10(1d-100), a very small value
          endif
        enddo
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
        iskip=0 ! to be changed to iskip=1 if src=0
        if(YMAX.lt.1d-100) iskip=1 ! case of src~0
        ymax=10.*ymax !1.03*ymax
cBH180714      ymin=ymax/1.e3
        ymin=ymax/1.e8

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)

        if(iskip.eq.0)then
        CALL PGSVP(R4P2,R4P8,R4P2,R4P5)
        IF (RPG1.lt.RPG2) THEN
           CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
           CALL PGBOX('BCNST',R40,0,'BCNSTL',R40,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
           CALL PGLAB('time (secs)','KO Src(electrons/cm\u3\d/sec)',' ')
        ENDIF
        endif ! iskip

        write(t_,10010) n,timet
        CALL PGMTXT('B',R46,R4MP1,R40,t_)
        write(t_,10011) rovera(lr_),ll,rr
        CALL PGMTXT('B',R47,R4MP1,R40,t_)

c-----------------------------------------------------------------------
      enddo   ! on ll=1,lrz

10010 format("After time step(n)=",i4,5x,"time=",1pe14.6," secs")
10011 format("r/a=",f7.4,",  radial bin=",i3,
     +       ",  radial position(R)=",1pe10.3," cm")


      return
      end subroutine pltrun
