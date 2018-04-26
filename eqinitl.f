c
c
      subroutine eqinitl
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

      REAL RILIN

c..................................................................
c     This routine does some minor initialization for
c     the "eq" module. Called after the namelist read.
c..................................................................

      if (lfield.gt.lfielda) lfield=lfielda
      nrc=(nnr-1)/2+1
      nzc=(nnz-1)/2+1
      zshift=0.0


CMPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      if (noplots.ne."enabled1") then
      CALL PGPAGE
      RILIN=0.
      CALL PGMTXT('T',-RILIN,0.,0.,"PARAMETER VALUES")
      
      write(t_,1000) 
 1000 format("EQUILIBRIUM model parameters:")
      RILIN=2.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
      
      write(t_,1001)
 1001 format("nnra,nnza give the Maximum size the eqdsk")
      RILIN=3.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
      
      write(t_,1002) nnra,nnza
 1002 format("====>NNRA = ",i5,"        ====>NNZA = ",i5)
      RILIN=4.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)
      
      write(t_,1003) nconteqa
 1003 format("====>NCONTEQA = ",i5)
      RILIN=5.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      endif


      return
      end
