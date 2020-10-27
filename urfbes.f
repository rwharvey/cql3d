
c
c
      subroutine urfbes
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     Routine computes three Bessel tables, J_nharm-1,J_nharm, and
c     J_nharm+1. Space is allocated. Once tables are created, subsequent
c     calls check to see if the maximum argument has increased beyond
c     table limits. If so, table is augmented.
c..................................................................

      include 'param.h'
      include 'comm.h'
      
      pointer besl
      dimension besl(:) ! Now it's a local working array

      nharmx = 1
      do krf=1,mrfn
         nharmx=MAX(nharmx,nharm(krf))
      enddo
      write(*,*)'urfbes:  Before allocation, nharmx=',nharmx
      allocate(besl(nharmx+2),STAT=istat) !-YuP->added
      call bcast(besl,0.0,SIZE(besl))     !-YuP->added
      allocate(jbm1(nbssltbl,mrfn),STAT=istat) !-YuP->added
      call bcast(jbm1,0.0,SIZE(jbm1))          !-YuP->added
      allocate(jb0(nbssltbl,mrfn),STAT=istat)  !-YuP->added
      call bcast(jb0,0.0,SIZE(jbm1))           !-YuP->added
      allocate(jbp1(nbssltbl,mrfn),STAT=istat) !-YuP->added
      call bcast(jbp1,0.0,SIZE(jbm1))          !-YuP->added
      write(*,*)'urfbes:  After allocation'

c..................................................................
c     Loop over excitation modes
c..................................................................

      do 200 krf=1,mrfn
c         write(*,*)'urfbes: krf(1:mrfn) = ',krf

c..................................................................
c     Find maximum k-perp
c..................................................................

        call aminmx(wnper(1,1,krf),1,nrayelts*nrayn,1,emin,emax,
     +        kmin,kmax)
        em=emax
        if (abs(emin).gt.emax) em=abs(emin)
        xkp=em*freqcy(krf)*twopi/clight

c..................................................................
c     Guess at minimum cyclotron frequency - leave bpol out of estimate
c     this will decrease guess.
c..................................................................

c        bvalmin=fpsiz(lrindx(lrors))/er(nnr)
cc**bh930729wcemin=charge*bvalmin/(fmass(kelecg)*clight)
c        wcemin=abs(bnumb(1))*charge*bvalmin/(fmass(1)*clight)
c**bh970620  Using ray data to get bvalmin:
c**bh970620  (This was necessary, to avoid neg. tor. B problem with 
c**bh970620   the RFP).

        bvalmin=ep100
        do iray=1,nray(krf)
           do is=1,nrayelt(iray,krf)
c990131              bvalmin=amin1(bvalmin,abs(sbtot(is,iray,krf)))
              bvalmin=min(bvalmin,abs(sbtot(is,iray,krf)))
           enddo
c           write(*,*)'urfbes: iray,krf,sbtot(is,iray,krf) = ',
c     1             iray,krf,(sbtot(is,iray,krf),is=1,nrayelt(iray,krf))
        enddo
c       k is species to which this krf mode is to be applied:
        k=nrfspecies(krfn(krf))
        wcemin=abs(bnumb(k))*charge*bvalmin/(fmass(k)*clight)

c..................................................................
c     Maximum argument in table will be argmax - give it a 40% boost
c..................................................................

        argmax=xkp*vnorm/wcemin*1.4
        if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
        write(*,*)'urfbes: argmax,nharm(krf),bvalmin = ',
     1                     argmax,nharm(krf),bvalmin
        endif
        bsslstp(krf)=argmax/(nbssltbl-1)
        arg=0.
        do 10 i=1,nbssltbl
          call zzbeslri(arg,nharm(krf)+2,0,besl,ncalc)
cBH150620: Having problem with bessel calc.
cBH150620: Try letting simple accuracy loss get by ok
cBH150620: 
           if (ncalc.ne.nharm(krf)+2) then
             write(*,*)'zzbeslri accuracy loss: arg,nharm(krf)+2,ncalc=' 
     1            ,arg,nharm(krf)+2,ncalc
             if (ncalc.lt.0) then
                stop ' in urfbes; ncalc.lt.0'
             endif
          endif
          if(nharm(krf).eq.0)  then
            jbm1(i,krf)=-besl(2)
            jb0(i,krf)=besl(1)
            jbp1(i,krf)=besl(2)
          elseif (nharm(krf).gt.0)  then
            jbm1(i,krf)=besl(nharm(krf))
            jb0(i,krf)=besl(nharm(krf)+1)
            jbp1(i,krf)=besl(nharm(krf)+2)
          endif
          arg=arg+bsslstp(krf)
 10     continue
        argmax=arg-bsslstp(krf)
        
 200  continue ! krf=1,mrfn
 
      write(*,*)'End urfbes'
 
      return
      end
