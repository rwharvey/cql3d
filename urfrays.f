c
c
      subroutine urfrays(initrf,nrayptlh,nrayptec,nrayptfw)
      implicit integer (i-n), real*8 (a-h,o-z)

c.................................................................
c     This routine calls ray tracing code(s), if namelist
c     variables, call_lh, call_ech, or call_fw are "enabled".
c.................................................................

      include 'param.h'
      include 'comm.h'

      save icall
      data icall /0/

c..................................................................
c     Determine whether to call the ray tracing code or just to assume
c     that a disk file with r.t. output already exists
c     xbr holds the name of the r.t. controllee.
c..................................................................

      if (call_lh .eq."enabled".and.lh.eq."enabled") then

c..................................................................
c     Set up a disk file to control LH ray tracing code.
c     initrf=0 for first call..
c     nraypts is the maximum number of ray elements allowed per ray.
c..................................................................

        if (nrayptlh.gt.nrayelts) nraypts=nrayelts
        open(unit=21,file='iraylh',status='new')
        write(21,102) initrf,nraypts,ieqbrurf
        close(unit=21)
        write (*,101) n,1
c990131        istat=ishell('xbr')
        write (*,100) n,1
      endif
      if (call_ech.eq."enabled".and.ech.eq."enabled") then

c..................................................................
c     Set up a disk file to control ECH ray tracing code.
c     initrf=0 for first call..
c     nraypts is the maximum number of ray elements allowed per ray.
c..................................................................

        if (nrayptec.gt.nrayelts) nraypts=nrayelts
        open(unit=22,file='irayec',status='new')
        write(22,102) initrf,nraypts,ieqbrurf
        close(unit=22)
        write (*,101) n,2
c990131        istat=ishell('toray')
        write (*,100) n,2
      endif
      if (call_fw.eq."enabled".and.fw.eq."enabled") then

c..................................................................
c     Set up a disk file to control FW ray tracing code.
c     initrf=0 for first call..
c     nraypts is the maximum number of ray elements allowed per ray.
c..................................................................

        if (nrayptfw.gt.nrayelts) nraypts=nrayelts
        open(unit=25,file='irayec',status='new')
        write(25,102) initrf,nraypts,ieqbrurf
        close(unit=25)
        write (*,101) n,2
c990131        istat=ishell('xbr')
        write (*,100) n,2
      endif
 100  format("Call for n=",i5, " mode=",i5,"complete for ray tracing")
 101  format("Call for n=",i5," mode=",i5, "begins for ray tracing")
 102  format(3i5)
      if (icall.eq.0) imprf=1
      icall=1
      return
      end
