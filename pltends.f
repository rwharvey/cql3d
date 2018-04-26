c
c
      subroutine pltends
      implicit integer (i-n), real*8 (a-h,o-z)
c
c     Plot energy, density, parallel current and density conservation
c     constant vs. time., at a given s, distance along the magnetic field.
c     A cqlp.eq.'enabled' routine.
c
c     Modified from Graflib to pgplot calls by Yuri Petrov, 090727,
c     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
c
      include 'param.h'
      include 'comm.h'

      REAL RILIN !-> For PGPLOT (text output positioning)
      dimension wk_nch(nonch)
c
      if (noplots.eq."enabled1") return
      if (pltend.eq."disabled") return
      dgts=1.e-8
      rr=rpcon(lr_) !rovera(lr_)*radmin  ! YuP[03-2016] changed to rpcon
      do 220 k=1,ngen
        if (pltend.eq."notplts") then
          goto 100
        elseif (pltend.eq."last") then
          if (n.lt.nstop) goto 100
        endif
c...  
cmnt  Generate plot "endn"
c...  
        call GXGLFR(0)    
        call aminmx(pdens(1,k,l_),1,nch(l_),1,emin,emax,kmin,kmax)
        if (abs(emin-emax).lt.emax*dgts) emax=emin+.001*abs(emin)
        if(emax.gt.0.) emax=emax*1.05 ! extend the upper range        
        !write(*,*) 'pltends-1:', ptime(1,l_),ptime(nch(l_),l_),emin,emax
        call GSVP2D(.2,.8,.6,.9) !---------------> 1st subplot
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlin$",ptime(1,l_),ptime(nch(l_),l_),emin,
     1       emax*1.05)
        call GPCV2D(ptime(1,l_),pdens(1,k,l_),nch(l_)) ! density(time)
        CALL PGLAB('time (sec)',' ',' ')
        !call GSCPVS(.08,.75) ! set current position for text
        write(t_,10120) k
10120   format("density(s) of species",i2)
        RILIN=0.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('LV',RILIN,0.,0.,t_) ! Left-Vertical

        call aminmx(pengy(1,k,l_),1,nch(l_),1,emin,emax,kmin,kmax)
        if (abs(emin-emax).lt.emax*dgts) emax=emin+.001*abs(emin)
        if(emax.gt.0.) emax=emax*1.05 ! extend the upper range        
        !write(*,*) 'pltends-2:', ptime(1,l_),ptime(nch(l_),l_),emin,emax
        call GSVP2D(.2,.8,.2,.5) !---------------> 2nd subplot
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlin$",ptime(1,l_),ptime(nch(l_),l_),emin,
     1       emax*1.05)
        call GPCV2D(ptime(1,l_),pengy(1,k,l_),nch(l_)) ! energy(time)
        CALL PGLAB('time (sec)',' ',' ')
        !call GSCPVS(.08,.35) ! set current position for text
        write(t_,10110) k
10110   format("energy(s) of species",i2)
        RILIN=0.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('LV',RILIN,0.,0.,t_) ! Left-Vertical

        write(t_,10140) denpar(k,ls_),enrgypa(k,ls_)
10140   format("local density(s) (/cm**2) = ",e15.6,
     +    ";  energy(s) (kev) =",e15.6)
        RILIN=3.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('B',RILIN,0.,0.,t_) ! Bottom
     
        write(t_,10150) n,timet,k
10150   format("time step (n) is",i5,5x,"time=",1pe14.6," secs",
     +         "   Species k=",i2)
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_) ! Bottom

        write(t_,10152) rovera(lr_),rr
10152   format("r/a=",e14.6,5x,"radial position (R)=",e14.6," cm")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_) ! Bottom

        write(t_,10153) sz(l_)
10153   format("parallel position (s) =",e14.6," cm")
        RILIN=RILIN+1.
        CALL PGMTXT('B',RILIN,0.,0.,t_) ! Bottom
        
c...  
cmnt  Generate plot "curr"
c...  
        call GXGLFR(0) ! new page
        call aminmx(pcurr(1,k,l_),1,nch(l_),1,emin,emax,kmin,kmax)
        if (abs(emin-emax).lt.emax*dgts) emax=emin+.001*abs(emin)
        if(emax.gt.0.) emax=emax*1.05 ! extend the upper range
        !write(*,*) 'pltends-3:', ptime(1,l_),ptime(nch(l_),l_),emin,emax
        call GSVP2D(.2,.8,.6,.9) 
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlin$",ptime(1,l_),ptime(nch(l_),l_),emin,emax)
        call GPCV2D(ptime(1,l_),pcurr(1,k,l_),nch(l_)) ! current_dens(time)
        CALL PGLAB('time (sec)',' ',' ')
        curramp=currm(k,l_)/3.e9
        write(t_,10160) k,curramp
10160   format("Local in s current density",";",
     &    "of species ",i2," =",e14.5,";",
     &    "units are Amps/cm**2")
        RILIN=3.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('B',RILIN,0.,0.,t_) ! Bottom
     

 100    continue

        if (pltend.eq."tplts") goto 220
c...  
cmnt  Generate plot "currv"
c...  
        call GXGLFR(0) ! new page
        call aminmx(currv(1,k,l_),1,jx,1,fnmin,fnmax,kmin,kmax)
        if (abs(fnmin-fnmax).lt.fnmax*dgts) fnmax=fnmin+.001*abs(fnmin)
        !write(*,*) 'pltends-4:', x(1),x(jx),fnmin,fnmax
        call GSVP2D(.2,.8,.6,.9) !---------------> 1st subplot
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlin$",x(1),x(jx),fnmin,fnmax)
        call GPCV2D(x,currv(1,k,l_),jx)
        CALL PGLAB('v/vnorm','current: J(v)',' ') ! 
        
        call aminmx(currvs(1,k),1,jx,1,fnmin,fnmax,kmin,kmax)
        if (abs(fnmin-fnmax).lt.fnmax*dgts) fnmax=fnmin+.001*abs(fnmin)
        !write(*,*) 'pltends-5:', x(1),x(jx),fnmin,fnmax
        call GSVP2D(.2,.8,.2,.5) !---------------> 2nd subplot
        CALL PGSCH(1.) ! set character size; default is 1.
        call GSWD2D("linlin$",x(1),x(jx),fnmin,fnmax)
        call GPCV2D(x,currvs(1,k),jx)
        CALL PGLAB('v/vnorm','Integral from 0 to v of J(v)dv',' ') !
         
        write(t_,10183) currvs(jx,k)
10183   format("current =",e10.4,"Amps/Watt")
        RILIN=3.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('B',RILIN,0.,0.,t_) ! Bottom

 220  continue ! k=1,ngen
c
c     Plot the density conservation diagnostic vs time
c
      if (pltend.eq."notplts") then
        if (n.lt.nstop) return
      elseif (pltend.eq."last") then
        if (n.lt.nstop) return
      endif
c...  
cmnt  Generate plot "consn(l_)"
c...  
      call GXGLFR(0) ! new page
      call aminmx(consnp(1,l_),1,nch(l_),1,emin,emax,kmin,kmax)
      !write(*,*) 'pltends-6:', ptime(1,l_),ptime(nch(l_),l_),emin,emax
      call GSVP2D(.2,.8,.6,.9)
      CALL PGSCH(1.) ! set character size; default is 1.
      call GSWD2D("linlin$",ptime(1,l_),ptime(nch(l_),l_),emin,emax)
      call GPCV2D(ptime(1,l_),consnp(1,l_),nch(l_))      
      CALL PGLAB('time (sec)','consn(l_) conservation diag',' ')
      
      write(t_,10250) consn(l_)
10250 format("consn(l_)=",1pe12.4)
      RILIN=5.
      CALL PGSCH(0.8) ! set character size; default is 1.
      CALL PGMTXT('B',RILIN,-.2,0.,t_)
        
      write(t_,10251) 
10251 format("Perfect conservation should yield  machine accuracy,")
      RILIN=RILIN+1.
      CALL PGMTXT('B',RILIN,-.2,0.,t_)
        
      write(t_,10252) 
10252 format("or about 1.e-14:")
      RILIN=RILIN+1.
      CALL PGMTXT('B',RILIN,-.2,0.,t_)

      write(t_,10150) n,timet
      RILIN=RILIN+1.
      CALL PGMTXT('B',RILIN,-.2,0.,t_)
        
      write(t_,10152) rovera(lr_),rr
      RILIN=RILIN+1.
      CALL PGMTXT('B',RILIN,-.2,0.,t_)
        
      write(t_,10153) sz(l_)
      RILIN=RILIN+1.
      CALL PGMTXT('B',RILIN,0.,0.,t_) 
      
      
      do 280 k=1,ngen
        call GXGLFR(0) ! new page(s)
c$$$    Possibly write t_ greater than present dimension character*512:
c$$$        write(t_,10260) k,(sgaint(i,k,l_),i=1,8)
c$$$        call gptx2d(t_)
c$$$10260   format("total gain (+) or loss (-) to date for species",i3,";",
c$$$     &    "(in particles/cm**2 - units of line density)",";",
c$$$     &    "forcing nonnegative f x-sweep (implct=disabled)",e12.5,";",
c$$$     &    "forcing nonnegative f y-sweep (implct=disabled)",e12.5,";",
c$$$     &    "due to particle source=",e12.5,";","due to runaway=",e12.5,
c$$$     +    ";",
c$$$     &    "due to lossmode(k)=",e12.5,";","due to torloss(k)=",e12.5,
c$$$     +    ";",
c$$$     &    "due to fusion losses=",e12.5,";",
c$$$     &    "forcing nonnegative distribution (implct=enabled)=",e12.4,
c$$$     +    ";","$")
 280  continue
        write(t_,10150) n,timet
        RILIN=0.
        CALL PGSCH(0.8) ! set character size; default is 1.
        CALL PGMTXT('T',RILIN,0.,0.,t_) ! Top

        write(t_,10152) rovera(lr_),rr
        RILIN=RILIN-1.
        CALL PGMTXT('T',RILIN,0.,0.,t_) ! Top

        write(t_,10153) sz(l_)
        RILIN=RILIN-1.
        CALL PGMTXT('T',RILIN,0.,0.,t_) ! Top
      
      return
      end
