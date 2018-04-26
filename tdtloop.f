c
c
      subroutine tdtloop
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     Determines whether or not it is time to halt execution.
c..................................................................

      include 'param.h'
      include 'comm.h'
      include 'name.h'
CMPIINSERT_INCLUDE

      character cptline*80
      integer getpid
      real :: cputime
      
      if (n .ge. nstop) go to 20 ! Finished; Finalize and Stop program. 
      go to 30 ! Not finished yet; Proceed to return/end, next time step
      
 10   itloop=1
 20   continue
      call tdtscout
      call tdeqdsk

CMPIINSERT_IF_RANK_EQ_0
      ! make plots on mpirank.eq.0 only
      if (noplots.ne."enabled1") then
         call pgclos  ! PGCLOS
         PRINT *,'PGPLOT CLOSED at time step n=',n
      endif
CMPIINSERT_ENDIF_RANK
      
c..................................................................
c     Check-point this job, if ichkpnt.ne."disabled" into file ichkpnt.
c     This has been disabled, due to lack of generality of dropfile().
c..................................................................
      if (ichkpnt.ne."disabled") then
cBH080118        i=getpid(0)
        i=9999               !Dummy value
        write(cptline,100) i,ichkpnt
 100    format("chkpnt -p ",i5," -f ",a8)
        print *,cptline
c        i=dropfile(ichkpnt)
      endif

c.......................................................................
c     save distribution fnctn, restore original cqlinput, cputime, STOP
c.......................................................................

      if (nlwritf.ne."disabled") call tdwritef

CMPIINSERT_IF_RANK_EQ_0
         call ainadjnl(1)
         !restore cqlinput if &FSETUP was renamed to &SETUP0 (earlier)
         call ainadjnl_fsetup_setup0(1) 
         !pause
CMPIINSERT_ENDIF_RANK

CMPIINSERT_BARRIER

      call cpu_time(cputime)
CMPIINSERT_IF_RANK_EQ_0
      PRINT *, 'tdtloop: >>> END OF CQL3D PROGRAM <<<'
      PRINT *, 'tdtloop: Execution time (seconds)', cputime
CMPIINSERT_ENDIF_RANK
      

 30   continue
      return
      end
