c
c
      subroutine tdsxrplt(en,eflux,nen,nenaa,
     +                    efluxt,nv,inegsxr,softxry,lnwidth)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

      include 'param.h'

CMPIINSERT_INCLUDE

      REAL*4 RTAM1(nena),RTAM2(nena)
      REAL*4 REMAX,REMIN
cBH092022:
      REAL*4 RBOUND

c..................................................................
cmnt  this routine plots SXR/NPA energy/particle flux spectra 
cmnt  versus photon/particle energy.
cmnt  It would be simpler to have a separate tdnpaplt [BH100815].
c..................................................................

      dimension en(nenaa),eflux(nenaa,*),efluxt(nv),inegsxr(nv)
      character*1024 t_
      character*8 softxry

c      write(*,*)'tdsxrplt: en(1:nen)',
c     +     en(1:nen)
c      write(*,*)'tdsxrplt: eflux(1:nen,1:nv)',
c     +    (eflux(1:nen,i),i=1,nv)
c      write(*,*)'tdsxrplt: efluxt(1:nv),inegsxr(1:nv)',
c     +     efluxt(1:nv),inegsxr(1:nv)
c      write(*,*)'tdsxrplt: softxry=',softxry

CMPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      ep100=1.d100
      em100=1.d-100

      fmin=ep100
      fmax=-ep100
      do 100  nn=1,nv
        if (inegsxr(nn).le.1) go to 100
        call aminmx(eflux(1,nn),1,inegsxr(nn),1,fmin1,fmax1,kmin,kmax)
        fmin=min(fmin,fmin1)
        fmax=max(fmax,fmax1)
 100  continue
      if (fmin .eq. fmax) fmin=.9*fmax-1.e-20
      decades=6.1
      if(decades.ge.0.)  fmin=fmax/10**decades
      
      DO J=1,nen
         RTAM1(J)=RBOUND(en(j))
      ENDDO

      REMIN=RBOUND(LOG10(fmin))
      REMAX=RBOUND(LOG10(fmax))
      
      CALL PGPAGE
C     CALL PGENV(Rtam1(1),Rtam1(nen),Remin,Remax,0,20)
      CALL PGSVP(.2,.8,.45,.9)
      CALL PGSWIN(Rtam1(1),Rtam1(nen),Remin,Remax)
      CALL PGBOX('BCNST',0.,0,'BCNSTL',0.,0)
      CALL PGSAVE
      CALL PGSCH(1.44)
      if (softxry.eq.'enabled') then
      CALL PGLAB('Photon Energy k (keV)', 
     +     'd\u2\d\(0555)/dtdk (ergs/cm\u2\d-sec-ster-eV)',
     +     'SXR Energy Flux versus Photon Energy')
      else
      CALL PGLAB('Particle Energy k (keV)', 
     +     'd\u2\dN/dtdk (#/cm\u2\d-sec-ster-eV)',
     +     'NPA Flux versus Energy')
      endif

      CALL PGUNSA

      do 200 nn=1,nv
        if (inegsxr(nn) .le. 1) go to 200
        CALL PGSLW(lnwidth) ! line thickness/width
       DO J=1,inegsxr(nn)
           rtam2(J)=rbound(eflux(j,nn))
c           write(*,*)'tdsxrplt.f: nn,j,eflux(j,nn),rtam2(j):',
c     +                            nn,j,eflux(j,nn),rtam2(j)
           RTAM2(J)=ABS(RTAM2(J))
           if (RTAM2(J).gt.em100) then
              RTAM2(J)=LOG10(RTAM2(J))
           else
              RTAM2(j)=1.0
           endif
        ENDDO
           CALL PGLINE(inegsxr(nn),RTAM1,RTAM2)
 200   continue
       CALL PGSLW(lnwidth) ! restore line thickness/width

      if (softxry.eq.'enabled') then 
         write(t_,610)
      else
         write(t_,613)
      endif
 610  format("total flux, enmin to enmax (ergs/cm**2-sec-ster):")
 613  format("total flux, enmin_npa to enmax_npa (#/cm**2-sec-ster):")
         
      CALL PGMTXT('B',7.,-0.1,0.,t_)

      do 300  nn=1,nv,2
         if (nn.eq.nv .and. ((nv/2)*2 .ne. nv)) then
            write(t_,612) efluxt(nv)
         else
            write(t_,611)  (efluxt(in),in=nn,nn+1)
         endif
         CALL PGMTXT('B',7.5+0.5*nn,0.,0.,t_)
 300  continue
c$$$      do 300  nn=1,nv,4
c$$$         if (nn.eq.nv .and. ((nv/4)*4 .ne. nv)) then
c$$$            write(t_,612) efluxt(nv)
c$$$         else
c$$$            write(t_,611)  (efluxt(in),in=nn,nn+3)
c$$$         endif
c$$$         CALL PGMTXT('B',7.+nn,0.,0.,t_)
c$$$ 300  continue
 611  format(1p4e12.2)
 612  format(1pe12.2,12x)

      return
      end


c
c
      subroutine tdsxrvw(tempp4,tempp5,tempp6)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

      character*8  pltsxrvw

      REAL*4 ZTOP
      REAL*4 RTAB1(LFIELDA),RTAB2(LFIELDA)
      REAL*4 PGER1,PGERNNR,PGEZNNZ

      REAL*4 RRTAB1,RRTAB2
      dimension rrtab1(:), rrtab2(:)
      dimension tempp4(*),tempp5(*),tempp6(*)
      pointer rrtab1, rrtab2 
      allocate(rrtab1(iyjx),STAT=istat) 
      allocate(rrtab2(iyjx),STAT=istat) 


c..................................................................
c     This routine plots out the contours (flux surfaces) in a 
c     poloidal plane, and overplots the SXR view cords.
c..................................................................

CMPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      pltsxrvw="enabled"
      if (pltsxrvw.eq."disabled") return


      CALL PGPAGE

      if (eqsym.ne."none") then
         ztop=2.*.5*ez(nnz)/(er(nnr)-er(1))+.05
      else
         ztop=.95
      endif

      CALL PGSVP(.15,.85,.15,ztop)

      PGER1=er(1)
      PGERNNR=er(nnr)
c      write(*,*)'tdsxrvw: PGER1,PGERNNR=',PGER1,PGERNNR
      PGEZNNZ=ez(nnz)
      CALL PGSWIN(PGER1,PGERNNR,-PGEZNNZ,PGEZNNZ)
      CALL PGWNAD(PGER1,PGERNNR,-PGEZNNZ,PGEZNNZ)
      CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)
      CALL PGSLW(lnwidth) ! line thickness/width

      if (softxry.eq."enabled") then
         CALL PGLAB('Major radius (cms)','Vert height (cms)',
     +        'Flux Surfaces and SXR Chords')
      else  ! NPA
         CALL PGLAB('Major radius (cms)','Vert height (cms)',
     +        'Flux Surfaces and NPA Chords')
      endif

      do 10 l=1,lrzmax
         IF (LORBIT(L).GT.LFIELDA) STOP 'tdsxrvw: CHECK DIM OF RTAB1/2'

        if(eqsym.ne."none") then 
c     Up-Down symmetric flux surfaces are plotted:
           do 20 j=1,lorbit(l)
              RTAB1(j)=solr(lorbit(l)+1-j,l)
              RTAB2(j)=abs(solz(lorbit(l)+1-j,l))
 20        continue
           CALL PGLINE(LORBIT(L),RTAB1,RTAB2)
           do 21 j=1,lorbit(l)
              RTAB1(j)=solr(lorbit(l)+1-j,l)
              RTAB2(j)=-abs(solz(lorbit(l)+1-j,l))
 21        continue
           CALL PGLINE(LORBIT(L),RTAB1,RTAB2)
           
        else
c     Full flux surface contours are plotted:
           write(*,*)
           write(*,*)'tdsxrplt: l,lorbit(l)=',l,lorbit(l)
           do j=1,lorbit(l)
cBH091031             write(*,*)'j,solr(j,l),solz(j,l)=',
cBH091031    +             j,solr(j,l),solz(j,l)
           enddo
           do j=1,lorbit(l)
              RTAB1(j)=solr(lorbit(l)+1-j,l)
              RTAB2(j)=solz(lorbit(l)+1-j,l)
           enddo
           CALL PGLINE(LORBIT(L),RTAB1,RTAB2)
        endif
 10   continue

      if(eqsym.eq.'none' .and. ncontr.gt.1) then
        ! YuP[2015/05/03] Add LCFS, if available
        ncontr_= min(ncontr,LFIELDA)
        do ilim=1,ncontr_
           RTAB1(ilim)=rcontr(ilim)
           RTAB2(ilim)=zcontr(ilim)
        enddo
        CALL PGLINE(ncontr_,RTAB1,RTAB2)
      endif

      iistep=0
      do 30 nn=1,nv

        do 40 j=1,lensxr(nn)
           if(j.le.iyjx)  then
             RRTAB1(j)=sqrt(tempp4(iistep+j)**2+tempp5(iistep+j)**2)
             RRTAB2(j)=tempp6(iistep+j)
           else
             write(*,*)'tdsxrvw: Increase dimension of rrtab1,2'
             STOP
           endif
 40     continue

c        write(*,*) 'tdsxrvw: RRTAB1',(RRTAB1(jj),jj=1,lensxr(nn))
c        write(*,*) 'tdsxrvw: RRTAB2',(RRTAB2(jj),jj=1,lensxr(nn))

        CALL PGLINE(lensxr(nn),RRTAB1,RRTAB2)

        iistep=iistep+lensxr(nn)

 30   continue
 
      deallocate(rrtab1,STAT=istat) 
      deallocate(rrtab2,STAT=istat) 
 
      return
      end


