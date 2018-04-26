c
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
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))

        enddo
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
c        if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)
c        call pltcycl(iymin,iymax,ymin,ymax)

        ymax=1.03*ymax
        ymin=ymax/1.e4

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)



c      Following plot is with pdena,... with 1st dimension up to nch(l_).

        CALL PGPAGE
        CALL PGSVP(.2,.8,.6,.9)
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN
           RPG2= RPG1+1.e-16
        ENDIF
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
        IF (RPG1.ne.RPG2) THEN
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
        ENDIF
        CALL PGLAB('time (secs)','RA density (cm\u-3\d)',
     +       'Runaway Density and Current vs. Time')

cBH070405       do nt=1,nch(l_)
        do nt=1,nch(ll)
          xg(nt)=ptime(nt,1)
          yg(nt)=abs(pcurra(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))

        enddo
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)

c        if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)

        ymax=1.03*ymax
        ymin=ymax/1.e4

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)


        CALL PGSVP(.2,.8,.2,.5)
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN
           RPG2= RPG1+1.e-16
        ENDIF
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
        IF (RPG1.ne.RPG2) THEN
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
        ENDIF
        CALL PGLAB('time (secs)','RA curr den (A/cm\u2\d)',' ')


        write(t_,10010) n,timet
        CALL PGMTXT('B',6.,-.1,0.,t_)
        write(t_,10011) rovera(lr_),ll,rr
        CALL PGMTXT('B',7.,-.1,0.,t_)

c-----------------------------------------------------------------------


cBH070405       do nt=1,nch(l_)
        do nt=1,nch(ll)
          xg(nt)=ptime(nt,1)
          yg(nt)=abs(pfdenra(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))

        enddo
c        write(*,*) 'ptime(ijk,1),pfdenra(ijk,ll),ijk=1,5',
c     +              (ptime(ijk,1),pfdenra(ijk,ll),ijk=1,5)
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
c        if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)
c        call pltcycl(iymin,iymax,ymin,ymax)

        ymax=1.03*ymax
        ymin=ymax/1.e4

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)


        CALL PGPAGE
        CALL PGSVP(.2,.8,.6,.9)
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
        IF (RPG1.ne.RPG2) THEN
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
        ENDIF
        CALL PGLAB('time (secs)','Fraction RA density',
     +       'Runaway Fraction of Density and Current vs. Time')

        do nt=1,nch(l_)
          xg(nt)=ptime(nt,1)
          yg(nt)=abs(pfcurra(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))

         enddo
        call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
c        if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)
c        call pltcycl(iymin,iymax,ymin,ymax)

        ymax=1.03*ymax
        ymin=ymax/1.e4

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)



        CALL PGSVP(.2,.8,.2,.5)
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
        IF (RPG1.ne.RPG2) THEN
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
        ENDIF
        CALL PGLAB('time (secs)','Fraction RA curr den',' ')

        write(t_,10010) n,timet
        CALL PGMTXT('B',6.,-.1,0.,t_)
        write(t_,10011) rovera(lr_),ll,rr
        CALL PGMTXT('B',7.,-.1,0.,t_)

c-----------------------------------------------------------------------




c     Plots of ucrit and e/e_dreicer, versus time:

cBH070405       do nt=1,nch(l_)
        do nt=1,nch(ll)
           xg(nt)=ptime(nt,1)
           yg(nt)=abs(pucrit(nt,ll))
           RNONCHA1(nt)=RBOUND(XG(nt))
           RNONCHA2(nt)=RBOUND(YG(nt))
           RNONCHA2(nt)=LOG10(RNONCHA2(nt))
        enddo
      
      
      if (pltlim.eq.'u/c') then
         do nt=1,nch(l_)
cYuP            yg(nt)=yg(nt)*cnormi !cnormi=0 when relativ.eq."disabled"
          yg(nt)=yg(nt)/cnorm ! YuP[07-2016]
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
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
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
         enddo
      endif
             
      call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
c     if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)
c     call pltcycl(iymin,iymax,ymin,ymax)
      
      
      ymin=0.97*ymin
      ymax=1.e3*ymin

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
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

        CALL PGPAGE
        CALL PGSVP(.2,.8,.6,.9)
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
        IF (RPG1.ne.RPG2) THEN
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
        ENDIF
        CALL PGLAB('time (secs)',t_,
     +       'Critical vel(energy) and E/E\dDreicer\u vs. Time')

cBH070405      do nt=1,nch(l_)
      do nt=1,nch(ll)
         xg(nt)=ptime(nt,1)
         yg(nt)=abs(peoed(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
      enddo
      call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
c     if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)
c     call pltcycl(iymin,iymax,ymin,ymax)
      
      ymax=1.03*ymax
      ymin=ymax/1.e3

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)

      write(t_,1023)
 1023 format("E-field/E_Dreicer")
      write(t_,1024)
 1024 format("time(sec)")


        CALL PGSVP(.2,.8,.2,.5)
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
        IF (RPG1.ne.RPG2) THEN
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
        ENDIF
        CALL PGLAB('time (secs)','E-field/E\dDreicer\u',' ')

        write(t_,10010) n,timet
        CALL PGMTXT('B',6.,-.1,0.,t_)
        write(t_,10011) rovera(lr_),ll,rr
        CALL PGMTXT('B',7.,-.1,0.,t_)

c-----------------------------------------------------------------------


c     Plots of e/e0 and KO source, versus time:

cBH070405      do nt=1,nch(l_)
      do nt=1,nch(ll)
         xg(nt)=ptime(nt,1)
         yg(nt)=abs(peoe0(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
      enddo
      
      call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
c     if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)
c     call pltcycl(iymin,iymax,ymin,ymax)
      
c$$$      call gxglfr(0)
      
      ymax=1.03*ymax
      ymin=ymax/1.e3

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)


        CALL PGPAGE
        CALL PGSVP(.2,.8,.6,.9)
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
        IF (RPG1.ne.RPG2) THEN
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
        ENDIF
        CALL PGLAB('time (secs)','E-field/E0',
     +       'E/(Critical E0) and KO Source vs. Time')


cBH070405      do nt=1,nch(l_)
      do nt=1,nch(ll)
         xg(nt)=ptime(nt,1)
         yg(nt)=abs(psrc(nt,ll))
          RNONCHA1(nt)=RBOUND(XG(nt))
          RNONCHA2(nt)=RBOUND(YG(nt))
          RNONCHA2(nt)=LOG10(RNONCHA2(nt))
      enddo
      call aminmx(yg,1,nch(l_),1,ymin,ymax,kmin,kmax)
c     if (abs(ymin-ymax).le.abs(ymax)*dgts) ymax=ymin+0.001*abs(ymin)
c     call pltcycl(iymin,iymax,ymin,ymax)
      
      ymax=1.03*ymax
      ymin=ymax/1.e3

        RPG1=RBOUND(Ymin)
        RPG2=RBOUND(Ymax)
        IF ( YMIN.EQ.YMAX ) THEN
           RPG1=.1
           RPG2=1.
        ENDIF
        RPG1=LOG10(RPG1)
        RPG2=LOG10(RPG2)


        CALL PGSVP(.2,.8,.2,.5)
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
        IF (RPG1.ne.RPG2) THEN
           CALL PGBOX('BCNST',0.0,0,'BCNSTL',0.0,0)
           CALL PGLINE(NCH(L_),RNONCHA1,RNONCHA2)
        ENDIF
        CALL PGLAB('time (secs)','KO Src(electrons/cm\u3\d/sec)',' ')

        write(t_,10010) n,timet
        CALL PGMTXT('B',6.,-.1,0.,t_)
        write(t_,10011) rovera(lr_),ll,rr
        CALL PGMTXT('B',7.,-.1,0.,t_)

c-----------------------------------------------------------------------

      
      enddo   ! on ll:1,lrz

10010 format("After time step(n)=",i4,5x,"time=",1pe14.6," secs")
10011 format("r/a=",f7.4,",  radial bin=",i3,
     +       ",  radial position(R)=",1pe10.3," cm")


      return
      end subroutine pltrun
