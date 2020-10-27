c
c
      subroutine sigv5d
      implicit integer (i-n), real*8 (a-h,o-z)
c
c.......................................................................
c     Called for each flux surface.
c     This routine computes fusion reaction rate (n1*n2*sigma-v-bar),
c     for various combinations of species with density n1 and n2,
c     as specified by input variable isigmas(1:6).
c     written by m.g. mccoy, modified for CQL3D (bobH, 941102)
c
c     isigsgv1=0, reaction rate for general species colliding with self
c                 or other distribution.
c              =1, reaction rate for Max. distn. at same temp and 
c                  density as the general species, colliding with self
c                  and other  distributions.
c.......................................................................

cBH110330:  Might want to do this calculation with the number of
cBH110330:  Legendre poly expansion different from mx, similar
cBH110330:  to SXR case.

      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

      !YuP[06-2016] Moved outside of if[] below.
      !Both arrays should be initialized to 0.
      call bcast(sigf(1,lr_),zero,4)
      call bcast(sigm(1,lr_),zero,4)
        
      if (isigsgv1 .eq. 0) then
        mmx1=mmsv
      else
        mmx1=0
      endif


c     Loop over nuclear reactions
      iq=0

      do 90 knumb=1,4
      if (isigmas(knumb) .eq. 0) go to 90
      iq=iq+1
c     Species 1 (i.e., "a"):
      kk=igenrl(1,knumb)
      k=imaxwln(1,knumb)
      if(k.eq.0 .and. kk.eq.0) go to 90 !-> next reaction type
      !write(*,*)'sigv5d:Reactant#1: knumb,igenrl,imaxwln=',knumb,kk,k

cBH120314: Removed isigsvg2 option, since no sensible use for it.
cBH120314: Also, maxwln distributions added below for temp1 had
cBH120314: different density normalization:
cBH120314: Distributions now normalized int*d^3x{f}=reden.
cBH120314:
cBH120314:c     Put a background Maxwellian in temp1 with energy of k,
cBH120314:c     density of kk,  if isigsvg2.eq.1:
cBH120314:      call bcast(temp3(0,0),zero,iyjx2)
cBH120314:      if (k .ne. 0 .and. isigsgv2 .eq. 1) call sigmaxwl(k,kk)
cBH120314:      call dcopy(iyjx2,temp3(0,0),1,temp1(0,0),1)

      call bcast(temp3(0,0),zero,iyjx2)
      if (kk .eq. 0) then  !No species 'a' general distrn.  Use maxwl.
            !write(*,*)'sigv5d:Reactant#1: Use Maxwellian with k=',k
         bn1=reden(k,lr_)
         call sigmax(k) !gives Maxwl species "a" in temp3
      else  !i.e., kk.ne.0
         bn1=reden(kk,lr_)
         if (isigsgv1 .eq. 0) then
            !write(*,*)'sigv5d:Reactant#1: Use f(..kk..)  for kk=',kk
            do j=0,jx+1
               do i=0,iy+1
                  temp3(i,j)=f(i,j,kk,l_)
               enddo
            enddo
         else
c           Generate Maxwellian in temp3 with energy and density of kk:
            call sigmax(kk)
         endif
      endif
      do j=0,jx+1
         do i=0,iy+1
            temp1(i,j)=temp3(i,j)
         enddo
      enddo

      if (knumb .lt. 3) then  !That is, different species interacting.
                              !Reactant "b" different than "a".
c       Species 2, (i.e., "b"):
         kk=igenrl(2,knumb)
         k=imaxwln(2,knumb)
         if(k.eq.0 .and. kk.eq.0) go to 90 !-> next reaction type
         !write(*,*)'sigv5d:Reactant#2: knumb,igenrl,imaxwln=',knumb,kk,k

         call bcast(temp3(0,0),zero,iyjx2)
         if (kk.eq.0) then
               !write(*,*)'sigv5d:Reactant#2: Use Maxwellian with k=',k
            bn2=reden(k,lr_)
            call sigmax(k)      !Maxwl in temp3
         else                   !i.e., kk.ne.0
            bn2=reden(kk,lr_)
            if (isigsgv1 .eq. 0) then
               !write(*,*)'sigv5d:Reactant#2: Use f(..kk..)  for kk=',kk
               do j=0,jx+1
                  do i=0,iy+1
                     temp3(i,j)=f(i,j,kk,l_)
                  enddo
               enddo
            else
c              Generate Maxwln in temp3 with energy and density of kk:
               call sigmax(kk)
            endif
         endif
         do j=0,jx+1
            do i=0,iy+1
               temp2(i,j)=temp3(i,j)
            enddo
         enddo

      else   !knumb=3 or 4

         bn2=bn1
         do j=0,jx+1
            do i=0,iy+1
               temp2(i,j)=temp1(i,j)
            enddo
         enddo
         
      endif

c     Factor to account for like-like collisions (BobH:950307)
      if (knumb.eq.3 .or. knumb.eq.4) then
        flklk=0.5
      else
        flklk=1.0
      endif


c     Loop along field line
      do 100 l=1,lz
        sum=0.

c     Loop over Legendre polynomials
      do l1=0,mmx1
      
      call bcast(fetb(1,l1),zero,jx)
      call bcast(feta(1,l1),zero,jx)
c     Copy distribution a into temp3, for cfpleg:
c%bh      call dcopy(iyjx2,temp4,1,temp3,1)
      call dcopy(iyjx2,temp1(0,0),1,temp3(0,0),1)
      do 202 j=1,jx
      tam3(j)=temp1(iy+1-imax(l,lr_),j)
202   tam2(j)=temp1(imax(l,lr_),j)
      call cfpleg(l1,l,1)
      call dcopy(jx,tam1,1,feta(1,l1),1)

      if (knumb .ge. 3) go to 210

c%bh      call dcopy(iyjx2,temp5,1,temp3,1)
      call dcopy(iyjx2,temp2(0,0),1,temp3(0,0),1)
      do 201 j=1,jx
      tam3(j)=temp2(iy+1-imax(l,lr_),j)
      tam2(j)=temp2(imax(l,lr_),j)
201   continue
      call cfpleg(l1,l,1)
 210  call dcopy(jx,tam1,1,fetb(1,l1),1)

      call bcast(tam6,zero,jx)
      do 60 jv1=1,jx
      if (jv1 .eq. 1) go to 50
      do 40 jv2=1,jv1-1
      tam6(jv2)=tam6(jv2)+(feta(jv1,l1)*fetb(jv2,l1)+feta(jv2,l1)
     1 *fetb(jv1,l1))*csv(iind(jv1)+jv2,l1,iq)
   40 continue
   50 sum=sum+feta(jv1,l1)*fetb(jv1,l1)*
     1                     csv(iind(jv1)+jv1,l1,iq)
 60   continue

      do 70 jv2=2,jx
 70   sum=sum+tam6(jv2) ! Reaction rate (n_a*n_b*<sigma*v>)
           ! summed over Legendre indices and dual-j indices;
           ! Units: [reactions/sec/cm^3]
           
      enddo ! l1 (Legendre pol. order)


c     sum=sum*bn1*bn2   (Normalization of f in cql3d changed from cql).
c     cql is the 1D in radius, pre-cql3d radial code.

      sum=sum*flklk ! Reaction rate (n_a*n_b*<sigma*v>)
           ! summed over Legendre indices and dual-j indices;
           ! Units: [reactions/sec/cm^3]
           ! This value (sum) is for a given surface (lr_), 
           ! given pol.angle (l),
           ! and given reaction type (knumb).
      
      if (isigsgv1 .eq. 0) then
        !fus_rate_sigf(lr_,l,knumb)= sum !local, at (lr,l) point
        !Flux surface (i.e. volume) averages:
        sigf(knumb,lr_)=sigf(knumb,lr_)
     +    +sum*dz(l,lr_)/bbpsi(l,lr_)/zmaxpsi(lr_)
      else
        !fus_rate_sigm(lr_,l,knumb)= sum !local, at (lr,l) point
        !Flux surface (i.e. volume) averages:
        sigm(knumb,lr_)=sigm(knumb,lr_)
     +    +sum*dz(l,lr_)/bbpsi(l,lr_)/zmaxpsi(lr_)
      endif

 100  continue ! lz

 90   continue ! knumb=1,4 reaction type

c     Watts/cm**3
      fuspwrm(1,lr_)=sigm(1,lr_)*17.6*1.602e-13
      fuspwrv(1,lr_)=sigf(1,lr_)*17.6*1.602e-13
      fuspwrm(2,lr_)=sigm(2,lr_)*18.3*1.602e-13
      fuspwrv(2,lr_)=sigf(2,lr_)*18.3*1.602e-13
      fuspwrm(3,lr_)=sigm(3,lr_)*4.03*1.602e-13
      fuspwrv(3,lr_)=sigf(3,lr_)*4.03*1.602e-13
      fuspwrm(4,lr_)=sigm(4,lr_)*3.27*1.602e-13
      fuspwrv(4,lr_)=sigf(4,lr_)*3.27*1.602e-13
c     17.6, etc : Energies of fusion products [MeV]

      if (ioutput(1).ge.1) then !YuP[2020] Useful diagnostic printout
      write(*,'(a,i5,4e12.3)')
     + 'sigv5d: lr_, fuspwrv(1:4,lr_)',lr_,fuspwrv(1:4,lr_)
      endif
      

c     BH120315:
c     Fusion reaction rates <sigma*v> are defined such that
c     the total volumetric rate (reactions /(sec*cm**3) is 
c     n_i*n_j*/(1+delta_ij)<sigma*v>.  The Dirac delta accounts
c     for double counting of the reactions, in the case of like-like
c     collisions.  See Bosch and Hale, NF, p. 611-631 (1992),
c     Eqns. 10-11.  The Direc delta is the flklk factor, above.
c     The sigf above are related to <sigma*v> by
c     <sigma*v>=(1+delta_ij)*fuspwrv/(n_i*n_j*Rx_en_per_part*1.602e-13)
c     (to be compared with Table VIII, p. 625 of Bosch and Hale.)
c     That is, <sigma*v>=sigf(,)*(1+delta_ij)/(n_i*n_j).

      return
      end
