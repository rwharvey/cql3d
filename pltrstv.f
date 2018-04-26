c
c
      subroutine pltrstv
      implicit integer (i-n), real*8 (a-h,o-z)
c
c     Plot electron resistivity and related quantities.
c

c     Modified from Graflib to pgplot calls by Yuri Petrov, 090727,
c     using PGPLOT + GRAFLIBtoPGPLOT.f routines (put in pltmain.f).
c
      include 'param.h'
      include 'comm.h'

      REAL RILIN !-> For PGPLOT (text output positioning)

c
      if (noplots.eq."enabled1") return

      if (kelecg .eq. 0 .or. abs(elecfld(lr_)) .lt. 1.e-10) go to 190

      call GXGLFR(0)    
      call aminmx(sptzrp(2,lmdpln_),1,nch(l_)-1,1,emin,emax,kmin,kmax)
      call aminmx(restp(2,lr_),1,nch(l_)-1,1,fmin,fmax,kmin,kmax)
      if (fmin .lt. emin) emin=fmin
      if (fmax .gt. emax) emax=fmax
      call GSVP2D(.2,.8,.6,.9) !---------------> 1st subplot
      CALL PGSCH(1.) ! set character size; default is 1.
      call GSWD2D("linlin$",ptime(1,l_),ptime(nch(l_),l_),
     + .95d0*emin,1.05d0*emax)
      text(1)="spitzer"
      !-YuP:   call GSCVLB(1)
      !-YuP:   call GSCVTX(loc(text))
      call GPCV2D(ptime(2,l_),sptzrp(2,lmdpln_),nch(l_)-1)
      text(1)="rstvty"
      call GPCV2D(ptime(2,l_),restp(2,lr_),nch(l_)-1)
      !-YuP:   call GSCVLB(0)
      write(t_,170)
 170  format("upper graph - flux avg and spitzer resistivities")
      RILIN=2.
      CALL PGSCH(0.9) ! set character size; default is 1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
      write(t_,171)
 171  format("lower graph - ratio of resist to spitzer or neo resist")
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      illeff=lr_
      if (cqlpmod .eq. "enabled") illeff=ls_
      call aminmx(rovsp(2,illeff),1,nch(l_)-1,1,emin,emax,kmin,kmax)
      
      call GSVP2D(.2,.8,.2,.5) !---------------> 2nd subplot
      call GSWD2D("linlin$",ptime(1,l_),ptime(nch(l_),l_),
     + .95d0*emin,1.05d0*emax)
      CALL PGLAB('time (sec)',' ',' ')
      call GPCV2D(ptime(2,l_),rovsp(2,illeff),nch(l_)-1)

      if (efswtchn.eq."neo_hh") then
         write(t_,10164) 
10164 format("(efswtchn=neo_hh)")
         RILIN=RILIN-1.
         CALL PGMTXT('T',RILIN,0.,0.,t_)
      endif



      call GXGLFR(0) ! new page
      call GSVP2D(0.2,0.8,0.2,0.7) ! (XLEFT, XRIGHT, YBOT, YTOP)       
      write(t_,70)
 70   format("--calculated resistivity and other related quantities--")
      RILIN=2.
      CALL PGSCH(0.9) ! set character size; default is 1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
     
      write(t_,80) sptzr(l_)
 80   format("spitzer restvty= ",1pe14.5)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
      
      write(t_,81) resist
 81   format("toroidal restvty=",e14.5)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
      
      write(t_,82) resistn
 82   format("neoclass restvty=",e14.5)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
      
      write(t_,83) rovs(lr_)
 83   format("ratio of resistivities=",e14.5)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
      
      write(t_,84) rovsf
 84   format("O(epsilon**.5) expansion for resistivity ratio=",e14.5)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
      
      write(t_,85) elecr(lr_)
 85   format("E-Dreicer=",e14.5,"vlts/cm")
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
      
      write(t_,86) rovsloc(l_)
 86   format("local resistivity over spitzer=",e14.5)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,90) rovscf     
 90   format("Small eps(lr_) fla for resist ratio (connor)=",e16.6)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,91) rovsc(lr_)
 91   format("gen. epsilon fla for resist ratio (connor)=",e14.5)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,92) xconn  
 92   format("^\i(connor)=",e16.6)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
      
      write(t_,93) elecfld(lr_)
 93   format("electric field=",e16.6,"vlts/cm")
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
      
      write(t_,94) eovedd
 94   format("E-parallel/E-Dreicer=",e16.6)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
      
      write(t_,95) tauee(lr_)
 95   format("tauee(lr_)=",e16.6,"secs")
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
c
c     add some other relevant quantities
c
      write(t_,96)(btor0(lr_)/bmod0(lr_))**2
 96   format("b_phi/b at outer midpplane=",1pe14.4)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,97) onovrp(2,lr_)*rpcon(lr_)**2
 97   format("R(z=0)**2 * <1/R**2>      =",  e14.4)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)

      write(t_,98) psiavg(2,lr_)
 98   format("<(B(z)/B(0))**2>          =",  e14.4)
      RILIN=RILIN-1.
      CALL PGMTXT('T',RILIN,0.,0.,t_)
     
     
 190  continue
      return
      end
