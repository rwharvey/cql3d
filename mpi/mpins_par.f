CMPIINSERT_INCLUDE
      include 'mpilib.h'
CMPIINSERT_START
      call init_mpi
CMPIINSERT_FINISH
      call close_mpi
CMPIINSERT_MPIWORKER
      if(soln_method.eq.'direct' .and. lrzmax.gt.1
     +     ) then
         ! Parallelization for the impavnc0 solver is limited 
         ! for soln_method='direct' (for now)
         if(mpisize.gt.1) then
            mpiworker= MOD(ll-1,mpisize-1)+1  !1...(mpisize-1)
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
      else
         ! In all other cases, perform calculations 
         ! for all flux surfaces on mpirank=0, then broadcast results
         mpiworker=0
      endif
c      if(mpirank.eq.mpiworker) then
c        write(*,*)'n,n_(ll),ll=',n,n_(ll),ll,' mpiworker=',mpiworker
c      endif
CMPIINSERT_MPIWORKER_KRF
         if(mpisize.gt.1) then
            mpiworker= MOD(krf-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
CMPIINSERT_MPIWORKER_IRAYKRF
      iraykrf= iray + nrayn*(krf-1) 
         if(mpisize.gt.1) then
            mpiworker= MOD(iraykrf-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
c      if(mpirank.eq.mpiworker) then
c      PRINT *,'n,iray,krf,mpiworker=',n,iray,krf,mpiworker
c      endif
CMPIINSERT_MPIWORKER_LFCT
         if(mpisize.gt.1) then
            mpiworker= MOD(lfct-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
CMPIINSERT_MPIWORKER_L
         if(mpisize.gt.1) then
            mpiworker= MOD(l-1,mpisize-1)+1 ! poloidal grid
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
CMPIINSERT_MPIWORKER_I
         if(mpisize.gt.1) then
            mpiworker= MOD(i-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
CMPIINSERT_MPIWORKER_J
         if(mpisize.gt.1) then
            mpiworker= MOD(j-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
CMPIINSERT_MPIWORKER_IV
         if(mpisize.gt.1) then
            mpiworker= MOD(iv-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
CMPIINSERT_MPIWORKER_IJ
         !ij= j + jx*(i-1) 
         ij= i + iy*(j-1)
         if(mpisize.gt.1) then
            mpiworker= MOD(ij-1,mpisize-1)+1
         else
            PRINT*, '------- WARNING: mpisize=1 -------'
            mpiworker=0
         endif
CMPIINSERT_BARRIER
      call MPI_BARRIER(MPI_COMM_WORLD,mpiierr)
CMPIINSERT_BCAST_RDC_GRID
      call MPI_BCAST(n_uprp,1,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(n_upar,1,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(n_psi,1,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(vc_cgs,1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(upar_min,1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(upar_max,1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
CMPIINSERT_BCAST_RDC
      call MPI_BCAST(rho_a,n_psi,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(uprp,n_uprp,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(upar,n_upar,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(rdc_cqlb,n_uprp*n_upar*n_psi,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(rdc_cqlc,n_uprp*n_upar*n_psi,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(rdc_cqle,n_uprp*n_upar*n_psi,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(rdc_cqlf,n_uprp*n_upar*n_psi,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
CMPIINSERT_BCAST_NRAYN
      call MPI_BCAST(nrayn,1,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(nray,nmodsa,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(nharm,nmodsa,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(freqcy,nmodsa,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(omega,nmodsa,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
CMPIINSERT_BCAST_NRAYELTS
      call MPI_BCAST(nrayelts,1,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(nrayelt,nrayn*mrfn,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
CMPIINSERT_BCAST_RAYS_DATA
      call MPI_BCAST(nharm1,nmodsa,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(nharms,nmodsa,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(jslofas,nrayn*mrfn,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(nurefls,nrayn*mrfn,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(keiks,nrayn*mrfn,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(jpes,nrayn*mrfn,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(jpis,nrayn*mrfn,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(istarts,nrayn*mrfn,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(iprmt5,nrayn*mrfn,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(jhlfs,nrayn*mrfn,
     +     MPI_INTEGER,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(sxxrt,nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(skpsi,nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(skth,nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(skphi,nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(delpwr,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(fluxn,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(seikon,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(spsi,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(sdpwr,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(sbtot,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(sene,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(salphac,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(salphal,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(ws,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(wr,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(wz,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(wnpar,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(wnper,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(wphi,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(wdnpar,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(cwexde,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(cweyde,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(cwezde,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiierr)
CMPIINSERT_BCAST_DISTRIBUTION
      call MPI_BCAST(f,iyjx2*ngen*lrors,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
CMPIINSERT_BCAST_COLL_COEFFS
      call MPI_BCAST(cal,iyjx*ngen*lrors,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(cbl,iyjx*ngen*lrors,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(ccl,iyjx*ngen*lrors,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(cdl,iyjx*ngen*lrors,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(cel,iyjx*ngen*lrors,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(cfl,iyjx*ngen*lrors,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(eal,iyjx*ngen*2*lrors,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(ebl,iyjx*ngen*2*lrors,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
CMPIINSERT_BCAST_SCAL
      call MPI_BCAST(scal,iyjx*ngen*lrors,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
CMPIINSERT_BCAST_VELSOU
      call MPI_BCAST(velsou,iyjx2*ngen*lrors,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
CMPIINSERT_BCAST_ENTR
      call MPI_BCAST(entr(k,-1,l_),1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,0,l_),1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,1,l_),1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,2,l_),1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,3,l_),1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,4,l_),1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,5,l_),1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,6,l_),1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,7,l_),1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,8,l_),1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,11,l_),1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(entr(k,12,l_),1,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(pwrrf(1,k,l_),jx,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
CMPIINSERT_IF_RANK_NE_0_RETURN
      if(mpirank.ne.0) return 
CMPIINSERT_IF_RANK_EQ_0
      if(mpirank.eq.0) then 
CMPIINSERT_IF_RANK_EQ_MPIWORKER
      if(mpirank.eq.mpiworker) then
CMPIINSERT_ENDIF_RANK
      endif  ! for if(mpirank.eq.***)
CMPIINSERT_SEND_RECV
      if(soln_method.eq.'direct' .and. lrzmax.gt.1) then
      ! Parallelization for the impavnc0 solver is limited 
      ! for soln_method='direct' (for now)
      if(mpirank.eq.0.or.mpirank.eq.mpiworker) then
         call send_data ! send or recv data on f and coll.coeffs.
      endif
      endif
CMPIINSERT_SEND_RECV_ENTR
      if(mpirank.eq.0.or.mpirank.eq.mpiworker) then
         call send_entr(k,lefct) 
         !send/recv entr(k,:,l_),pwrrf,pwrrfs(:,k,l_)
      endif
      
CMPIINSERT_SEND_URFPWR
c         PRINT*,'SEND_urfpwr: mpirank,krf,iray=',mpirank,krf,iray
         ! Count number of elements for this ray
         irayis=0 
         do is=1,nrayelt(iray,krf)! Loop over ray elements
            if(ipwr(is).eq.1) irayis=irayis+1 ! incremented up to nrayis
         enddo
         nrayis=irayis
         mpisz=nrayis ! number of elements 
         irayis=0 
         do is=1,nrayelt(iray,krf)! Loop over ray elements
           if(ipwr(is).eq.1)  then
             irayis=irayis+1 ! incremented up to nrayis
             urftmp(0*mpisz+irayis)= urfpwr(is,iray,krf)
             urftmp(1*mpisz+irayis)= urfpwrc(is,iray,krf)
             urftmp(2*mpisz+irayis)= urfpwrl(is,iray,krf)
             urftmp(3*mpisz+irayis)= scalurf(is,iray,krf)
             urftmp(4*mpisz+irayis)= salphac(is,iray,krf)
           endif
         enddo
         mpitag= iray + nrayn*(krf-1) ! combined: ray-number + wave-mode
         call MPI_SEND(ipwr, nrayelts,
     +        MPI_INTEGER1,
     +        0, mpitag, 
     +        MPI_COMM_WORLD,mpiierr)
         call MPI_SEND(urftmp, 5*mpisz,
     +        MPI_DOUBLE_PRECISION,
     +        0, mpitag, 
     +        MPI_COMM_WORLD,mpiierr)
c         PRINT*,'SEND_urfpwr: krf,iray,mpirank=',krf,iray,mpirank
      !-----------------------------------------------------------
      
CMPIINSERT_RECV_URFPWR
      if(mpirank.eq.0) then !-------------------------------------------
c         PRINT*,'recv_urfpwr: mpirank,krf,iray=',mpirank,krf,iray
         call MPI_RECV(ipwr, nrayelts,
     +        MPI_INTEGER1,
     +        MPI_ANY_SOURCE, MPI_ANY_TAG, 
     +        MPI_COMM_WORLD,mpistatus,mpiierr)
         mpitag=mpistatus(MPI_TAG) 
         call MPI_RECV(urftmp, nrayelts*5,
     +        MPI_DOUBLE_PRECISION,
     +        MPI_ANY_SOURCE, mpitag, 
     +        MPI_COMM_WORLD,mpistatus,mpiierr)
         mpiray=MOD(mpitag-1,nrayn)+1  ! determine which ray
         mpikrf=(mpitag-mpiray)/nrayn +1 ! determine which krf wave mode
         ! Get mpisz5 (Number of elements received)
         call MPI_GET_COUNT(mpistatus, 
     +        MPI_DOUBLE_PRECISION,mpisz5,mpiierr) 
         mpisz= mpisz5/5 ! urftmp contains 5 arrays
         irayis=0 
         do is=1,nrayelt(mpiray,mpikrf)! Loop over ray elements
           if(ipwr(is).eq.1)  then
             irayis=irayis+1 ! incremented up to nrayis==mpisz
             urfpwr(is,mpiray,mpikrf)=  urftmp(0*mpisz+irayis)
             urfpwrc(is,mpiray,mpikrf)= urftmp(1*mpisz+irayis) 
             urfpwrl(is,mpiray,mpikrf)= urftmp(2*mpisz+irayis) 
             scalurf(is,mpiray,mpikrf)= urftmp(3*mpisz+irayis) 
             salphac(is,mpiray,mpikrf)= urftmp(4*mpisz+irayis)
           endif
         enddo
c         PRINT*,'recv_urfpwr: mpikrf,mpiray,mpisz=',
c     +                        mpikrf,mpiray,mpisz
      endif !-----------------------------------------------------------
      
CMPIINSERT_BCAST_URFPWR
      call MPI_BCAST(urfpwr,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(urfpwrc,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(urfpwrl,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(scalurf,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(salphac,nrayelts*nrayn*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
     
CMPIINSERT_SEND_URFB0
      if(mpirank.eq.mpiworker) then !-----------------------------------
        mpisz=iyjx*lrz ! number of elements in urfb(i,j,lr)
        call dcopy(mpisz,urfb(1,1,1,krf),1,  urfbwk(0*mpisz+1),1) !         1 : mpisz
        call dcopy(mpisz,urfc(1,1,1,krf),1,  urfbwk(1*mpisz+1),1) ! 1*mpisz+1 : 2*mpisz
        mpisz3=2*mpisz ! the last elem. in above
        urfbwk(mpisz3+0*lrz+1 : mpisz3+1*lrz) = powrfl(1:lrz,krf) !linear damp.
        urfbwk(mpisz3+1*lrz+1 : mpisz3+2*lrz) = powrfc(1:lrz,krf) !coll.damp.
        mpitag= krf ! wave-mode
        call MPI_SEND(urfbwk, mpisz3+2*lrz,
     +       MPI_DOUBLE_PRECISION,
     +       0, mpitag, 
     +       MPI_COMM_WORLD,mpiierr)
      endif !-----------------------------------------------------------
CMPIINSERT_RECV_URFB0
      if(mpirank.eq.0) then !-------------------------------------------
        mpisz=iyjx*lrz ! number of elements in urfb(i,j,lr)
        mpisz3=2*mpisz ! storage size for urfb,urfc
        call MPI_RECV(urfbwk, mpisz3+2*lrz,
     +       MPI_DOUBLE_PRECISION,
     +       MPI_ANY_SOURCE, MPI_ANY_TAG, 
     +       MPI_COMM_WORLD,mpistatus,mpiierr)
        mpitag=mpistatus(MPI_TAG) 
        mpikrf=mpitag ! determine which krf wave mode
        ij=0 
        do ll=1,lrz
        call tdnflxs(lmdpln(ll))
        do j=1,jx
        do i=1,iy
           ij=ij+1 
           urfb(i,j,indxlr_,mpikrf)=urfb(i,j,indxlr_,mpikrf)
     +                             +urfbwk(0*mpisz+ij)
           urfc(i,j,indxlr_,mpikrf)=urfc(i,j,indxlr_,mpikrf)
     +                             +urfbwk(1*mpisz+ij)
        enddo
        enddo
        enddo ! ll
      endif !-----------------------------------------------------------

CMPIINSERT_BCAST_URFB0
      call MPI_BCAST(urfb,iyjx*lrz*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(urfc,iyjx*lrz*mrfn,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)


CMPIINSERT_SEND_EFLUXWK
      if(mpirank.eq.mpiworker) then !-----------------------------------
        call dcopy(mpisz,efluxwk(1,i),1,tem2(1),1)
        mpitag= i ! i=1,lrzmax
        call MPI_SEND(tem2, mpisz,
     +       MPI_DOUBLE_PRECISION,
     +       0, mpitag, 
     +       MPI_COMM_WORLD,mpiierr)
      endif !-----------------------------------------------------------
CMPIINSERT_RECV_EFLUXWK
      if(mpirank.eq.0) then !-------------------------------------------
        call MPI_RECV(tem2, mpisz,
     +       MPI_DOUBLE_PRECISION,
     +       MPI_ANY_SOURCE, MPI_ANY_TAG, 
     +       MPI_COMM_WORLD,mpistatus,mpiierr)
        mpitag=mpistatus(MPI_TAG) 
        mpil_=mpitag ! determine which radial surface sent the data
        call dcopy(mpisz,tem2(1),1,efluxwk(1,mpil_),1)
      endif !-----------------------------------------------------------
      
CMPIINSERT_SEND_FUS
      if(mpirank.eq.mpiworker) then !-----------------------------------
        mpisz=4
        call dcopy(mpisz, fuspwrv(1,lr_), 1,buff(0*mpisz+1),1)
        call dcopy(mpisz, fuspwrm(1,lr_), 1,buff(1*mpisz+1),1)
        call dcopy(mpisz, sigf(1,lr_),    1,buff(2*mpisz+1),1)
        call dcopy(mpisz, sigm(1,lr_),    1,buff(3*mpisz+1),1)
        mpitag= lr_ ! over flux surfaces =1,lrz
        call MPI_SEND(buff, mpisz*4,
     +       MPI_DOUBLE_PRECISION,
     +       0, mpitag, 
     +       MPI_COMM_WORLD,mpiierr)
      endif !-----------------------------------------------------------
CMPIINSERT_RECV_FUS
      if(mpirank.eq.0) then !-------------------------------------------
        mpisz=4
        call MPI_RECV(buff, mpisz*4,
     +       MPI_DOUBLE_PRECISION,
     +       MPI_ANY_SOURCE, MPI_ANY_TAG, 
     +       MPI_COMM_WORLD,mpistatus,mpiierr)
        mpitag=mpistatus(MPI_TAG) 
        mpil_=mpitag ! determine which radial surface sent the data
        call dcopy(mpisz, buff(0*mpisz+1),1, fuspwrv(1,mpil_),1)
        call dcopy(mpisz, buff(1*mpisz+1),1, fuspwrm(1,mpil_),1)
        call dcopy(mpisz, buff(2*mpisz+1),1, sigf(1,mpil_),   1)
        call dcopy(mpisz, buff(3*mpisz+1),1, sigm(1,mpil_),   1)
      endif !-----------------------------------------------------------

CMPIINSERT_BCAST_EFLUX
      call MPI_BCAST(eflux(1,nn),nena,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)

CMPIINSERT_BCAST_EFLUX_NPA
      call MPI_BCAST(eflux_npa(1,nn),nen_npa,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)

CMPIINSERT_BCAST_FUS
      call MPI_BCAST(fuspwrv(1,1),4*lrorsa,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(fuspwrm(1,1),4*lrorsa,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(sigf(1,1),4*lrorsa,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
      call MPI_BCAST(sigm(1,1),4*lrorsa,
     +     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiierr)
     
CMPIINSERT_STARTTIME
      if(mpirank.eq.0) then
         mpitime = MPI_WTIME()
      endif
CMPIINSERT_ENDTIME
      if(mpirank.eq.0) then
         write(*,*) 'MPI Full time =',MPI_WTIME()-mpitime
      endif
