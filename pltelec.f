c
c
c
      subroutine pltelec
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c
cmnt  this routine plots electron density as a function of poloidal angl
c
      include 'param.h'
      include 'comm.h'

      REAL RTAM1(LZA),RTAM2(LZA)
      REAL RPGMIN,RPGMAX
      REAL RILIN

      if (noplots.eq."enabled1") return
c$$$      call gxglfr(0)
      call aminmx(densz(1,ngen+1,negyrg,lr_),1,lz,1,fmin,fmax,kmin,kmax)
      if (fmin .eq. fmax) fmin=.9*fmax-1.e-20
c$$$      call gswd2d("linlin$",pol(1,lr_),pol(lz,lr_),fmin,fmax)
c$$$      call gsvp2d(.2,.8,.25,.95)
c$$$      call gpgr80("linlin$")

      CALL PGPAGE
      CALL PGSVP(.2,.8,.45,.95)
        
      DO L=1,LZ
         RTAM1(L)=pol(L,lr_)
      ENDDO

      RPGMIN=fmin
      RPGMAX=fmax

      IF ( RPGMAX-RPGMIN .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPGMAX= RPGMIN+1.e-16
      ENDIF
      CALL PGSWIN(RTAM1(1),RTAM1(LZ),RPGMIN,RPGMAX)

c$$$      do 3001 l=1,lz
c$$$        tz1(l)=densz(l,ngen+1,negyrg,lr_)
c$$$ 3001 continue
c$$$      call gpcv2d(pol(1,lr_),tz1,lz)

      do 3001 l=1,lz
         RTAM2(L)=densz(l,ngen+1,negyrg,lr_)
 3001 continue
      CALL PGLINE(LZ,RTAM1,RTAM2)

 3002 continue
c$$$      call gscvlb(0)
c$$$      call gstxno(100.)
c$$$      call gscpvs(.5,.2)

c$$$      write(t_,610) kelec,n,timet,xlndnz(ngen+1,negyrg)
c$$$      call gptx2d(t_)
      RILIN=1.
      CALL PGMTXT(B,RILIN,0.,0.,T_)
      write(t_,611) kelec,n,timet
      RILIN=RILIN+1.
      CALL PGMTXT(B,RILIN,0.,0.,T_)
      write(t_,612) xlndnz(ngen+1,negyrg)
      RILIN=RILIN+1.
      CALL PGMTXT(B,RILIN,0.,0.,T_)
 610  format("Density as a function of poloidal angle(=pi*z/zmax)")
 611  format("species ",i3, " (electrons)   n= ",i5,"  time= ",1pe14.4)
 612  format("density (line-integration) =",1pe16.5)
      return
      end
