c-------------------------------------------------------
      subroutine send_data
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'
      include 'mpilib.h'
      
      real*8, allocatable :: 
     +        buff(:) ! (mpifsz+10*mpicsz) Buffer for 11 arrays below
           
      mpifsz= iyjx2*ngen !For send/recv of f(0:iy+1,0:jx+1,1:ngen,lr_), 
                         !and velsou(0:iy+1,0:jx+1,1:ngen,lr_)
      mpicsz= iyjx*ngen  !For send/recv of cal(1:iy,1:jx,1:ngen,lr_)
                         !and other collisional coeffs., 
                         !and scal(1:iyjx*ngen,lr_)
      mpisz= mpifsz+11*mpicsz+mpifsz ! buffer size
      
      if(mpirank.eq.0) then ! receive data from other ranks
         if (.NOT.ALLOCATED(buff)) allocate(buff(mpisz))
         call MPI_RECV(buff, mpisz,
     +        MPI_DOUBLE_PRECISION,
     +        MPI_ANY_SOURCE, MPI_ANY_TAG, 
     +        MPI_COMM_WORLD,mpistatus,mpiierr)
         mpil_=mpistatus(MPI_TAG) ! determine which flux surface
         call dcopy(mpifsz,buff(1),1,f(0,0,1,mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+0*mpicsz+1),1,cal(1,1,1,mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+1*mpicsz+1),1,cbl(1,1,1,mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+2*mpicsz+1),1,ccl(1,1,1,mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+3*mpicsz+1),1,cdl(1,1,1,mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+4*mpicsz+1),1,cel(1,1,1,mpil_),1)
         call dcopy(mpicsz,buff(mpifsz+5*mpicsz+1),1,cfl(1,1,1,mpil_),1)
         call dcopy(mpicsz,
     +                   buff(mpifsz+6*mpicsz+1),1,eal(1,1,1,1,mpil_),1)
         call dcopy(mpicsz,
     +                   buff(mpifsz+7*mpicsz+1),1,eal(1,1,1,2,mpil_),1)
         call dcopy(mpicsz,
     +                   buff(mpifsz+8*mpicsz+1),1,ebl(1,1,1,1,mpil_),1)
         call dcopy(mpicsz,
     +                   buff(mpifsz+9*mpicsz+1),1,ebl(1,1,1,2,mpil_),1)
         call dcopy(mpicsz,
     +                   buff(mpifsz+10*mpicsz+1),1,scal(1,mpil_),1)
         !--- Velocity source for radial transport (Note: size=mpifsz)
         call dcopy(mpifsz,
     +                   buff(mpifsz+11*mpicsz+1),1,
     +                   velsou(0,0,1,mpil_),1)
c         PRINT*,'recv: mpirank,mpil_=',mpirank,mpil_
      else !-> all other ranks send data to rank 0
         if (.NOT.ALLOCATED(buff)) allocate(buff(mpisz))
         call dcopy(mpifsz,f(0,0,1,lr_),  1,  buff(1),1)
         call dcopy(mpicsz,cal(1,1,1,lr_),1,  buff(mpifsz+0*mpicsz+1),1)
         call dcopy(mpicsz,cbl(1,1,1,lr_),1,  buff(mpifsz+1*mpicsz+1),1)
         call dcopy(mpicsz,ccl(1,1,1,lr_),1,  buff(mpifsz+2*mpicsz+1),1)
         call dcopy(mpicsz,cdl(1,1,1,lr_),1,  buff(mpifsz+3*mpicsz+1),1)
         call dcopy(mpicsz,cel(1,1,1,lr_),1,  buff(mpifsz+4*mpicsz+1),1)
         call dcopy(mpicsz,cfl(1,1,1,lr_),1,  buff(mpifsz+5*mpicsz+1),1)
         call dcopy(mpicsz,eal(1,1,1,1,lr_),1,buff(mpifsz+6*mpicsz+1),1)
         call dcopy(mpicsz,eal(1,1,1,2,lr_),1,buff(mpifsz+7*mpicsz+1),1)
         call dcopy(mpicsz,ebl(1,1,1,1,lr_),1,buff(mpifsz+8*mpicsz+1),1)
         call dcopy(mpicsz,ebl(1,1,1,2,lr_),1,buff(mpifsz+9*mpicsz+1),1)
         call dcopy(mpicsz,scal(1,lr_),  1,  buff(mpifsz+10*mpicsz+1),1)
         !--- Velocity source for radial transport (Note: size=mpifsz)
         call dcopy(mpifsz,velsou(0,0,1,lr_),1,
     +                     buff(mpifsz+11*mpicsz+1),1)
         mpitag=lr_ ! tag == flux surface number
         call MPI_SEND(buff, mpisz,
     +        MPI_DOUBLE_PRECISION,
     +        0, mpitag, 
     +        MPI_COMM_WORLD,mpiierr)
c         PRINT*,'SEND: mpirank,lr_=',mpirank,lr_
      endif
      
      return
      end

c-------------------------------------------------------
      subroutine send_entr(k,lefct) 
      !send/recv entr(k,lefct,l_),pwrrf(1:jx,k,l_),pwrrfs(1:jx,k,l_)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'
      include 'mpilib.h'
      dimension buff(1+jx) 
      
      if(mpirank.eq.0) then ! receive data from other ranks
         call MPI_RECV(buff, 1+jx,
     +        MPI_DOUBLE_PRECISION,
     +        MPI_ANY_SOURCE, MPI_ANY_TAG, 
     +        MPI_COMM_WORLD,mpistatus,mpiierr)
         mpitag=mpistatus(MPI_TAG) 
         lefct_=mpitag-2 ! determine which lefct was sent
         entr(k,lefct_,l_)=buff(1) ! for a given lefct
         entr(k,4,l_)=entr(k,4,l_)+buff(1) ! sum
         if (lefct_.eq.3) then
           call dcopy(jx,buff(2),1,pwrrf(1,k,l_),1)
           pwrrfs(1,k,l_)=dx(1)*pwrrf(1,k,l_) 
           do j=2,jx  ! sum over j
             pwrrfs(j,k,l_)=pwrrfs(j-1,k,l_)+dx(j)*pwrrf(j,k,l_)
           enddo
         endif
      else !-> all other ranks send data to rank 0
         buff(1)=entr(k,lefct,l_) ! for a given lefct
         mpisz=1
         if (lefct.eq.3) then
           call dcopy(jx,pwrrf(1,k,l_),1,buff(2),1)
           mpisz=1+jx
         endif
         mpitag=lefct+2 ! 2 added to make mpitag>0 (lefct can be -1)
         call MPI_SEND(buff, mpisz,
     +        MPI_DOUBLE_PRECISION,
     +        0, mpitag, 
     +        MPI_COMM_WORLD,mpiierr)
      endif
      
      return
      end
c-------------------------------------------------------
