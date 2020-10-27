c
c
      subroutine tdpltmne
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c
c     This routine plots out data as a function of r/a.
c     (Also a time-dependent plot of fusion rate.)
c
      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE
      dimension trilr(lrza),tr1s(lrorsa),tr2s(lrorsa)
c     real*4 variables (and function rbound) for pgplot:
      REAL*4 RPG1,RPG2, RPGmin, RPGmax
      REAL*4 RLRZAP1(0:LRZA),RLRZAP11(0:LRZA),RLRZAP12(0:LRZA),
     +     RLRZAP13(0:LRZA), RLRZAP14(0:LRZA)
      REAL*4 RLRZAP(0:LRZA)
      REAL*4 RLRORSA(LRORSA),RLRORSA1(LRORSA),RLRORSA2(LRORSA)
      REAL*4 RNONCHA1(nonch),RNONCHA2(nonch)
      REAL*4 RBOUND

      REAL*4 :: R4P2=.2,R4P8=.8,R4P6=.6,R41P3=1.3,R41P44=1.44,R4P5=.5
      REAL*4 :: R40=0.,R41=1.,R42=2.,R43=3.,R44=4.,R45=5.,R46=6.
      REAL*4 :: R47=7.,R48=8.,R45P2=5.2

      character*16 t_horiz
      
      real*8 tr1i(0:lrza),tr2i(0:lrza),tr3i(0:lrza),tr4i(0:lrza) !YuP[2019]
      real*8 tr5i(0:lrza),tr6i(0:lrza),tr7i(0:lrza) !YuP[2019-12-20]
      !local arrays, to keep partially integrated (in rho) values

      data em33/1.d-33/

c.......................................................................

CMPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      if (noplots.eq."enabled1") return
      if (plt3d .ne. "enabled") return
      if (lrzmax .le. 2) return
      lrzevn=(lrzmax/2)*2
      dgts=1.e-8
c..................................................................
c     Determine the x-axis for plots (psi or rho - both normalized).
c..................................................................

      if (pltvs.eq."psi") then
        do 20 l=1,lrzmax
          tr(l)=(equilpsi(0)-equilpsi(l))/equilpsi(0)
 20     continue
        write(t_horiz,'(a3)') 'psi'
      else
        do 30 l=1,lrzmax
          tr(l)=rz(l)/radmin
 30     continue
        write(t_horiz,'(a6,a8,a1)') 'rho (=', radcoord, ')'
      endif

      do 21 l=1,lrz
        trilr(l)=tr(lrindx(l))
 21   continue



      do 300 k=1,ntotal
        if (k.gt.ngen .and. k.ne.kelec .and. n.gt.1) go to 300
        CALL PGPAGE
        CALL PGSVP(R4P2,R4P8,R4P2,R4P6)

        CALL PGSAVE
        CALL PGSCH(R41P44)
        CALL PGMTXT('T',R46,R4P5,R4P5,'DENSITIES (/CC) OF SPECIES')
        CALL PGUNSA

        write(t_,3002) k,kspeci(1,k),kspeci(2,k),n
        CALL PGMTXT('T',R43,R40,R40,t_)
 3002   format("species no. ",i2,4x,a8,2x,a8,2x," time step n=",i4)
        do 19 l=1,lrzmax
          tr1(l)=reden(k,l)
          tr2(l)=0.0
          tr3(l)=0.0
          if (zmaxpsi(l).ne.0.0 .and. k.le.ngen) then
            tr2(l)=xlndnv(k,l)/zmaxpsi(l)
            tr3(l)=xlndnr(k,l)/zmaxpsi(l)
          endif
 19     continue   
        fmin=0.
        fmax=0.
        call aminmx(tr1(1),1,lrzmax,1,fmin,fmax,kmin,kmax)
        if (fmin .ge. .95*fmax) fmin=.95*fmax
        call aminmx(tr2(1),1,lrzmax,1,fmin1,fmax1,kmin,kmax)
        if (fmin1.lt.fmin) fmin=fmin1
        if(fmax1.gt.fmax) fmax=fmax1
        call aminmx(tr3(1),1,lrzmax,1,fmin1,fmax1,kmin,kmax)
        if (fmin1.lt.fmin) fmin=fmin1
        if(fmax1.gt.fmax) fmax=fmax1
        fmax=fmax*1.05 ! extend range, in case the profile is flat
        fmin=0.0 ! Set lower limit to 0.
        if(fmax.le.0.0) fmax=1.0 !YuP[2019-10-28]To prevent Y1=Y2 args in PGSWIN
        !YuP: Sometimes the density of added Z=50 impurity is 0., initially

        DO I=1,LRZMAX
           RLRZAP1(I)=tr(i)
        ENDDO
        RPG1=fmin
        RPG2=fmax
        CALL PGSWIN(RLRZAP1(1),RLRZAP1(lrzmax),RPG1,RPG2)
        CALL PGBOX('BCNST',R40,0,'BCNST',R40,0)
        DO I=1,LRZMAX
           RLRZAP11(I)=tr1(i)
        ENDDO
        CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP11(1))
        CALL PGSLS(2)
        DO I=1,LRZMAX
           RLRZAP12(I)=tr2(i)
        ENDDO
        CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP12(1))
        CALL PGSLS(3)
        DO I=1,LRZMAX
           RLRZAP13(I)=tr3(i)
        ENDDO
        CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP13(1))
        CALL PGSLS(1)
        write(t_,4050) pltvs
 4050   format(a8)
        CALL PGSAVE
        CALL PGSCH(R41P44)
        CALL PGLAB(t_,'density (/cm\u3\d)',' ')
        CALL PGUNSA


 300  continue

      do 400 k=1,ntotal
        if (k.gt.ngen .and. n.ge.1) go to 400

        CALL PGPAGE
        CALL PGSVP(R4P2,R4P8,R4P2,R4P6)

        CALL PGSAVE
        CALL PGSCH(R41P44)
        CALL PGMTXT('T',R46,R4P5,R4P5,'ENERGIES OF SPECIES IN KEV')
        if(k.le.ngen) CALL PGMTXT('T',R45P2,R4P5,R4P5,
     1                '(Solid: <..>_FSA)')
        CALL PGUNSA
        write(t_,3002) k,kspeci(1,k),kspeci(2,k),n
        CALL PGMTXT('T',R43,R40,R40,t_)

        if (n .eq. 0) then
          do 18 l=1,lrzmax
            tr1(l)=energy(k,l) ! FSA, solid line
            tr2(l)=0
c     tr3 set to zero so no need to change following plot statements
            tr3(l)=0
 18       continue
          if (k .le. ngen) then
            do 181 ll=1,lrz
              tr2(lrindx(ll))=energym(k,ll) ! midplane, dashed line
 181        continue
          endif
        else if (k .le. ngen) then ! and n>0
          do 182 l=1,lrzmax
            tr1(l)=energy(k,l) ! FSA, solid line
            tr2(l)=0.0
            tr3(l)=0.0
 182      continue
          do 183 l=1,lrz
            ilr=lrindx(l)
            tr2(ilr)=energyv(k,ilr) !?? FSA, transport related,  dashed line
            tr3(ilr)=energyr(k,ilr) !?? FSA, transport related,  dot-dashed line
 183      continue
        endif

        fmin=0.
        fmax=0.
        call aminmx(tr1(1),1,lrzmax,1,fmin,fmax,kmin,kmax)
        if (fmin .ge. .95*fmax) fmin=.95*fmax
        call aminmx(tr2(1),1,lrzmax,1,fmin1,fmax1,kmin,kmax)
        if (fmin1.lt.fmin) fmin=fmin1
        if(fmax1.gt.fmax) fmax=fmax1
        call aminmx(tr3(1),1,lrzmax,1,fmin1,fmax1,kmin,kmax)
        if (fmin1.lt.fmin) fmin=fmin1
        if(fmax1.gt.fmax) fmax=fmax1
        fmax=fmax*1.05 ! extend range, in case the profile is flat

        DO I=1,LRZMAX
           RLRZAP1(I)=tr(i)
           RLRZAP11(I)=tr1(i)
           RLRZAP12(I)=tr2(i)
           RLRZAP13(I)=tr3(i)
        ENDDO
        RPG1=fmin
        RPG2=fmax
        RPG1=min(fmin,0.)		
        CALL PGSWIN(RLRZAP1(1),RLRZAP1(lrzmax),RPG1,RPG2)
        CALL PGBOX('BCNST',R40,0,'BCNST',R40,0)
        
        CALL PGSLS(1)  !solid
        CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP11(1))
        
        CALL PGSLS(2)  !Set dashed line style
        CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP12(1))
        
        CALL PGSLS(3)  !Set dot-dash line style
        CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP13(1))
        
        CALL PGSLS(1)  !Re-set solid line style for annotation
        
        write(t_,4050) pltvs
        CALL PGSAVE
        CALL PGSCH(R41P44)
        CALL PGLAB(t_,'energy (kev)',' ')
        CALL PGUNSA

 400  continue ! k

      if (cqlpmod .eq. "enabled") then

        do 350 k=1,ntotal
          if (k.gt.ngen .and. k.ne.kelec .and. n.gt.1) go to 350

        CALL PGPAGE
        CALL PGSVP(R4P2,R4P8,R4P2,R4P6)

        CALL PGSAVE
        CALL PGSCH(R41P44)
        CALL PGMTXT('T',R46,R4P5,R4P5,'DENSITIES (/CC) OF SPECIES')
        CALL PGUNSA
          write(t_,3050) k,kspeci(1,k),kspeci(2,k)
        CALL PGMTXT('T',R43,R40,R40,t_)
          write(t_,3051) lrindx(1),rz(lrindx(1))/radmin
        CALL PGMTXT('T',R42,R40,R40,t_)
 3050     format("species no. ",i2,2x,a8,2x,a8)
 3051     format("r(",i2,")/a=",1pe11.4)


          do 352 ll=1,lsmax
            tr1s(ll)=denpar(k,ll)
 352      continue   
          fmin=0.
          fmax=0.
          call aminmx(tr1s(1),1,lsmax,1,fmin,fmax,kmin,kmax)
          if (fmin .ge. .95*fmax) fmin=.95*fmax
          fmin=0.d0 ! YuP: make lower limit =0.0
          fmax=fmax*1.05 ! extend range, in case the profile is flat

        DO I=1,LSMAX
           RLRORSA(I)=z(i,lrindx(1))
           RLRORSA1(I)=tr1s(i)
        ENDDO
        RPG1=fmin
        RPG2=fmax
      ! If the horizontal coord is rho, set the limits to [0.,1.]
      RPGmin=RLRORSA(1)
      RPGmax=RLRORSA(LSMAX)
      if(RPGmin.le.0.2) RPGmin=0. ! Lower limit in plots: extend to 0.
      if(RPGmax.ge.0.8 .and. RPGmax.lt.1.) RPGmax=1. ! Upper limit: extend to 1.
        
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
        CALL PGSWIN(RPGmin,RPGmax,RPG1,RPG2)
        CALL PGBOX('BCNST',R40,0,'BCNST',R40,0)
        CALL PGLINE(lsmax,RLRORSA(1),RLRORSA1(1))

          write(t_,4051)
 4051     format("s (cms)")

          CALL PGLAB(t_,'Density (/cm\u3\d)',' ')

 350    continue ! k

        do 450 k=1,ntotal
          if (k.gt.ngen .and. n.gt.1) go to 450

        CALL PGPAGE
        CALL PGSVP(R4P2,R4P8,R4P2,R4P6)

        CALL PGSAVE
        CALL PGSCH(R41P44)
        CALL PGMTXT('T',R46,R4P5,R4P5,'ENERGIES OF SPECIES IN KEV')
        CALL PGMTXT('T',R45,R40,R40,
     +         'Solid: midplane;  Dashed: <..>_FSA')
        CALL PGUNSA
          write(t_,3050) k,kspeci(1,k),kspeci(2,k)
        CALL PGMTXT('T',R43,R40,R40,t_)
          write(t_,3051) lrindx(1),rz(lrindx(1))/radmin
        CALL PGMTXT('T',R42,R40,R40,t_)

          do 452 ll=1,lsmax
            tr1s(ll)=enrgypa(k,ll) ! Midplane
 452      continue
          fmin=0.
          fmax=0.
          call aminmx(tr1s(1),1,lsmax,1,fmin,fmax,kmin,kmax)
          if (fmin .ge. .95*fmax) fmin=.95*fmax
          fmin=0.d0 ! YuP: make lower limit =0.0
          fmax=fmax*1.05 ! extend range, in case the profile is flat
	  
        DO I=1,LSMAX
           RLRORSA(I)=z(i,lrindx(1))
           RLRORSA1(I)=tr1s(i)
        ENDDO
        RPG1=fmin
        RPG2=fmax
      ! If the horizontal coord is rho, set the limits to [0.,1.]
      RPGmin=RLRORSA(1)
      RPGmax=RLRORSA(LSMAX)
      if(RPGmin.le.0.2) RPGmin=0. ! Lower limit in plots: extend to 0.
      if(RPGmax.ge.0.8 .and. RPGmax.lt.1.) RPGmax=1. ! Upper limit: extend to 1.
        
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
        CALL PGSWIN(RPGmin,RPGmax,RPG1,RPG2)
        CALL PGBOX('BCNST',R40,0,'BCNST',R40,0)
        CALL PGLINE(lsmax,RLRORSA(1),RLRORSA1(1))

          write(t_,4051)

          CALL PGLAB(t_,'Energy (keV)',' ')


 450    continue ! k

      endif

c.......................................................................
c     n > 0
c.......................................................................
      CALL PGSLS(1) ! restore solid line
      CALL PGSLW(lnwidth) ! restore
      CALL PGSCI(1) ! black color restored

      if (n .gt. 0) then

        CALL PGPAGE ! YuP[2019-12-20] Added plot of E-field
        CALL PGSVP(R4P2,R4P8,R4P2,R4P6)
        CALL PGSAVE
        CALL PGSCH(R41P44)
        CALL PGMTXT('T',R48,R4P5,R4P5,'Electric field (V/cm)')
        CALL PGUNSA
        n4=n ! add plot for present time step
        t4=ptime(n4,1)
        n3=max(1,NINT(n*2./3.)) ! add plot for step =n*2/3 
        t3=ptime(n3,1)
        n2=max(1,NINT(n*1./3.)) ! add plot for step =n*1/3 
        t2=ptime(n2,1)
        n1=1 ! add plot for initial time step
        t1=ptime(n1,1)
        fmin=0.
        fmax=0.
        DO I=1,LRZMAX
           RLRZAP1(I)=tr(i) ! rho or psi coord.
           tr1(i)=pefld(n1,I) ! at n=n1 ! V/cm
           tr2(i)=pefld(n2,I) ! at n=n2 ! V/cm
           tr3(i)=pefld(n3,I) ! at n=n3 ! V/cm
           tr4(i)=pefld(n4,I) ! at n=n4 ! V/cm
           RLRZAP11(I)=tr1(i)
           RLRZAP12(I)=tr2(i)
           RLRZAP13(I)=tr3(i)
           RLRZAP14(I)=tr4(i)
           fmin= min(fmin,tr1(i),tr2(i),tr3(i),tr4(i))
           fmax= max(fmax,tr1(i),tr2(i),tr3(i),tr4(i))
        ENDDO
        !Extend vertical range:
        dfrange=fmax-fmin !YuP[2019]
        fmin= fmin-0.02*dfrange -1.e-10 !YuP[2019]
        fmax= fmax+0.02*dfrange +1.e-10 !YuP[2019]
        RPG1=fmin
        RPG2=fmax
        !YuP[2019-12-20] If the horizontal coord is rho, set limits to [0.,1.]
        RPGmin=RLRZAP1(1)
        RPGmax=RLRZAP1(LRZMAX)
        if(RPGmin.le.0.2) RPGmin=0. ! Lower limit in plots: extend to 0.
        if(RPGmax.ge.0.8 .and. RPGmax.lt.1.) RPGmax=1. !Upper limit:extend to 1.
        CALL PGSWIN(RPGmin,RPGmax,RPG1,RPG2)
        CALL PGBOX('BCNST',R40,0,'BCNST',R40,0)
        
        CALL PGSLS(2) !dashed --
        CALL PGSCI(2) !red color
        CALL PGLINE(lrzmax,RLRZAP1(1:lrzmax),RLRZAP11(1:lrzmax)) ! n=n1
        write(t_,5004) n1,t1
        CALL PGMTXT('T',R47,R40,R40,t_) !YuP[2019-12-20]

        CALL PGSLS(3) ! -.-
        CALL PGSCI(3) !green color
        CALL PGLINE(lrzmax,RLRZAP1(1:lrzmax),RLRZAP12(1:lrzmax)) ! n=n2
        write(t_,5004) n2,t2
        CALL PGMTXT('T',R46,R40,R40,t_) !YuP[2019-12-20]
        
        CALL PGSLS(4) ! ...
        CALL PGSCI(4) !blue color
        CALL PGLINE(lrzmax,RLRZAP1(1:lrzmax),RLRZAP13(1:lrzmax)) ! n=n3
        write(t_,5004) n3,t3
        CALL PGMTXT('T',R45,R40,R40,t_) !YuP[2019-12-20]

        CALL PGSLS(1) !solid
        CALL PGSCI(1) !black color
        CALL PGLINE(lrzmax,RLRZAP1(1:lrzmax),RLRZAP14(1:lrzmax)) ! n=n4
        write(t_,5004) n4,t4
        CALL PGMTXT('T',R44,R40,R40,t_) !YuP[2019-12-20]

        CALL PGSLS(1) !solid, restore
        CALL PGSCI(1) !black color, restore
        write(t_,4050) pltvs ! pltvs="rho" or "psi"
        CALL PGLAB(t_,'E (V/cm)',' ')
        ! YuP[2019-12-20] Done E-field
        

        do 500 k=1,ngen

 5001     format(2("curr(",e12.5,")=",1pe15.5," Amps per cm**2",3x),"$")
 5011     format(("curr(",e12.5,")=",1pe15.5," Amps per cm**2",3x),"$")

        CALL PGPAGE ! Current densities
        CALL PGSVP(R4P2,R4P8,R4P2,R4P6)
        CALL PGSAVE
        CALL PGSCH(R41P44)
        CALL PGMTXT('T',R48,R4P5,R4P5,
     +              'FLUX SURF. AV. CURNT. (AMPS/CM\u2\d)')
        CALL PGUNSA
        
        write(t_,5002) k, currza(k)
 5002 format("Species:",i2,
     &        " Current from sum[curr*darea]=",1pe14.6," A")
 5003 format(1pe11.3,"A")
 5004 format("n=",i5, ";  t=",1pe14.6,"sec")
        CALL PGMTXT('T',R48,R40,R40,t_)

        !YuP[2019-12-20] Added j_bs in same plots
        ! Choose which j_bs component to plot: 
        ! thermal(=maxwellian) or non-thermal(=general species)
        if(kelecg.ne.0)then  
         kke=2 ! I_bs for e_general (non-thermal)
        else
         kke=1 ! I_bs for e_maxw
        endif
        if(niong.ne.0)then 
         kki=2 ! I_bs for i_general (non-thermal)
        else
         kki=1 ! I_bs for i_maxw
        endif
        if (jhirsh.eq.0) then
         ! In this case, j_bs is only calculated for maxwellian part.
         ! (In fact, it is done for electrons only)
         kke=1
         kki=1
        endif

        RLRZAP14=0.0 ! YuP[2019-12-20] Added, to sum-up different currents
        tr4=0.0 ! YuP[2019-12-20] Added, to sum-up different currents
        fmin=0.
        fmax=0.
        tr1i=0.d0 ! YuP[2019] to keep partially integrated values
        tr2i=0.d0 ! YuP[2019] to keep partially integrated values
        tr3i=0.d0 ! YuP[2019] to keep partially integrated values
        tr4i=0.d0 ! YuP[2019] to keep partially integrated values
        tr5i=0.d0 ! YuP[2019] to keep partially integrated values
        tr6i=0.d0 ! YuP[2019] to keep partially integrated values
        tr7i=0.d0 ! YuP[2019] to keep partially integrated values
        do l=1,lrzmax
            tr1(l)=currz(k,l)  !=curr(k,l)/3.9   [A/cm^2] 
            tr2(l)=currr(k,l)  !=curr(k,l)/3.e9 in case of transp=enabled
            tr3(l)=currv_(k,l) !=curr(k,l)/3.e9 in case of transp=enabled
            !YuP[2019-12-20] sumup all kind of currents for this k-species:
            tr4(l)= currz(k,l)
     &            +(currpar_starnue_n(l)-currpar_starnue0_n(l))
     &            + bscurm(l,1,kke) + bscurm(l,2,kki)
            !bscurm(*,1,*) is for electrons, bscurm(*,2,*) is for ions
            !currpar_starnue_n is from neoclassical resistivity,
            !                  including collisionality.
            !currpar_starnue0_n is based on starnue-->0 (no coll.limit)
            fmin=min(fmin,tr1(l),tr4(l),bscurm(l,1,kke),bscurm(l,2,kki))
            fmax=max(fmax,tr1(l),tr4(l),bscurm(l,1,kke),bscurm(l,2,kki))
            tr1i(l)= tr1i(l-1) +darea(l)*currz(k,l) ![Amps]
            tr4i(l)= tr4i(l-1) +darea(l)*tr4(l) ![Amps]
            tr5i(l)= tr5i(l-1) +darea(l)*
     &               (currpar_starnue_n(l)-currpar_starnue0_n(l)) ![Amps]
            tr6i(l)= tr6i(l-1) +darea(l)*bscurm(l,1,kke) ![Amps]
            tr7i(l)= tr7i(l-1) +darea(l)*bscurm(l,2,kki) ![Amps]
        enddo
        !call aminmx(tr4(1),1,lrzmax,1,fmin,fmax,kmin,kmax) !YuP[2019] tr4
        !Extend vertical range:
        dfrange=fmax-fmin !YuP[2019]
        fmin= fmin-0.02*dfrange -1.d-6 !YuP[2019]
        fmax= fmax+0.02*dfrange +1.d-6 !YuP[2019]
        RPG1=fmin
        RPG2=fmax
        
        DO I=1,LRZMAX
           RLRZAP1(I)=tr(i)
        ENDDO
        !YuP[2019-12-20] If the horizontal coord is rho, set limits to [0.,1.]
        RPGmin=RLRZAP1(1)
        RPGmax=RLRZAP1(LRZMAX)
        if(RPGmin.le.0.2) RPGmin=0. ! Lower limit in plots: extend to 0.
        if(RPGmax.ge.0.8 .and. RPGmax.lt.1.) RPGmax=1. !Upper limit:extend to 1.
        CALL PGSWIN(RPGmin,RPGmax,RPG1,RPG2)
        CALL PGBOX('BCNST',R40,0,'BCNST',R40,0)
        
        CALL PGSLS(1) ! solid 
        CALL PGSLW(lnwidth)  ! thin line
        CALL PGSCI(1) ! black color
        DO I=1,LRZMAX
           RLRZAP11(I)=tr1(i) !=curr(k,l)/3.e9   [A/cm^2]
           RLRZAP14(I)=RLRZAP14(I)+tr1(i) !YuP[2019-12-20]sum-up currents
        ENDDO
        CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP11(1))
        write(t_,5003) tr1i(lrzmax)
        CALL PGMTXT('T',R44,R40,R40,
     &   'Solid/thin: Integral over f (curr() array)'//t_) !YuP[2019-12-20]
        
        if(transp.eq.'enabled')then
        CALL PGSLS(2) ! dashed
        CALL PGSLW(lnwidth) ! thin line
        CALL PGSCI(1) ! black color
        DO I=1,LRZMAX 
           RLRZAP12(I)=tr2(i) !=curr(k,l)/3.e9  [A/cm^2], transp=enabled
        ENDDO
        CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP12(1))
        CALL PGSLS(3) ! dash-dot-dash
        CALL PGSLW(lnwidth)  ! thin line
        CALL PGSCI(1) ! black color
        DO I=1,LRZMAX
           RLRZAP13(I)=tr3(i) !=curr(k,l)/3.e9  [A/cm^2], transp=enabled
        ENDDO
        CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP13(1))
        endif ! (transp.eq.'enabled')
        
        !YuP[2019-12-20] added plot of currpar_starnue_n()-currpar_starnue0_n()
        ! This current addition is from (sigma_coll-neo - sigma_banana)*E_phi
        ! See sub. tddiag, around YuP[2019-12-19] .
        ! Or should we just add it to curr/3e9 ?
        if(k.eq.kelecg)then ! Show on plots of k=1 species only (electrons
          CALL PGSLS(2) ! dashed
          CALL PGSLW(lnwidth)  ! thin line
          CALL PGSCI(2) ! red color
          DO I=1,LRZMAX
           RLRZAP13(I)=currpar_starnue_n(I)-currpar_starnue0_n(I) ![A/cm^2]
           RLRZAP14(I)=RLRZAP14(I)+RLRZAP13(I) !YuP[2019-12-20]sum-up currents
          ENDDO
          CALL PGLINE(lrzmax,RLRZAP1(1:lrzmax),RLRZAP13(1:lrzmax))
          write(t_,5003) tr5i(lrzmax)
          CALL PGMTXT('T',R46,R40,R40,
     &    'Red-- (sigma_coll_neo-sigma_banana)*Ephi'//t_) !YuP[2019-12-20]
        endif ! k.eq.kelecg
        
        if(k.eq.kelecg)then !this k is for general electrons
          CALL PGSLS(3) ! dash-dot-dash
          CALL PGSLW(lnwidth)  ! thin 
          CALL PGSCI(3) ! green color
          DO I=1,LRZMAX
           RLRZAP13(I)=bscurm(i,1,kke) !bscurm(*,1,*) is for electrons
           RLRZAP14(I)=RLRZAP14(I)+RLRZAP13(I) !YuP[2019-12-20]sum-up currents
          ENDDO
          CALL PGLINE(lrzmax,RLRZAP1(1:lrzmax),RLRZAP13(1:lrzmax))
          write(t_,5003) tr6i(lrzmax)
          CALL PGMTXT('T',R45,R40,R40,
     &    'Green-.- Bootstrap (fit model: bscurm())'//t_) !YuP[2019-12-20]
        endif
        
        if(k.eq.niong)then !this k is for general ions 
          CALL PGSLS(4) ! dotted ....
          CALL PGSLW(lnwidth)  ! thin 
          CALL PGSCI(4) ! blue color
          DO I=1,LRZMAX
           RLRZAP13(I)=bscurm(i,2,kki) !bscurm(*,2,*) is for ions
           RLRZAP14(I)=RLRZAP14(I)+RLRZAP13(I) !YuP[2019-12-20]sum-up currents
          ENDDO
          CALL PGLINE(lrzmax,RLRZAP1(1:lrzmax),RLRZAP13(1:lrzmax))
          write(t_,5003) tr7i(lrzmax)
          CALL PGMTXT('T',R45,R40,R40,
     &    'Blue/dotted: Bootstrap (fit model: bscurm() array)'//t_) !YuP[2019-12-20]
        endif
        
        !ALL TOGETHER, for given species k
        CALL PGSCI(1) ! black color
        CALL PGSLS(1) ! solid 
        CALL PGSLW(lnwidth*2)  ! bold line
        CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP14(1))
        write(t_,5003) tr4i(lrzmax)
        CALL PGMTXT('T',R43,R40,R40,
     &   'Solid/bold: All the above together'//t_) !YuP[2019-12-20]
     
        CALL PGSLS(1) ! restore solid line
        CALL PGSLW(lnwidth) ! restore
        CALL PGSCI(1) ! black color restored
        write(t_,4050) pltvs
        CALL PGLAB(t_,'Curr Den (A/cm\u2\d)',' ')



        CALL PGPAGE !Current (A) INTEGRATED UP TO RHO or PSI !YuP[2019-12-20]
        CALL PGSVP(R4P2,R4P8,R4P2,R4P6)
        CALL PGSAVE
        CALL PGSCH(R41P44)
        CALL PGMTXT('T',R48,R4P5,R4P5,
     +              'Current (A) INTEGRATED UP TO RHO or PSI')
        CALL PGUNSA
        write(t_,5002) k, currza(k) !Species, and total current from curr()
        CALL PGMTXT('T',R48,R40,R40,t_)

        fmin=min(0.,MINVAL(tr1i),MINVAL(tr4i),MINVAL(tr6i),MINVAL(tr7i))
        fmax=max(0.,MAXVAL(tr1i),MAXVAL(tr4i),MAXVAL(tr6i),MAXVAL(tr7i))
        !Extend vertical range:
        dfrange=fmax-fmin !YuP[2019]
        fmin= fmin-0.02*dfrange -1.d-6 !YuP[2019]
        fmax= fmax+0.05*dfrange +1.d-6 !YuP[2019]
        RPG1=fmin
        RPG2=fmax
        CALL PGSWIN(RPGmin,RPGmax,RPG1,RPG2)
        CALL PGBOX('BCNST',R40,0,'BCNST',R40,0)
        
        CALL PGSLS(1) ! solid 
        CALL PGSLW(lnwidth)  ! thin line
        CALL PGSCI(1) ! black color
        DO I=1,LRZMAX
           RLRZAP11(I)=tr1i(I) !=From curr(k,l)/3.e9
        ENDDO
        CALL PGLINE(lrzmax,RLRZAP1(1:lrzmax),RLRZAP11(1:lrzmax))
        write(t_,5003) tr1i(lrzmax)
        CALL PGMTXT('T',R44,R40,R40,
     &   'Solid/thin: From Integral over f (curr())'//t_) !YuP[2019-12-20]
                
        ! Plot of currpar_starnue_n()-currpar_starnue0_n()
        ! This current addition is from (sigma_coll-neo - sigma_banana)*E_phi
        ! See sub. tddiag, around YuP[2019-12-19] .
        ! Or should we just add it to curr/3e9 ?
        if(k.eq.kelecg)then ! Show on plots of k=1 species only (electrons
          CALL PGSLS(2) ! dashed
          CALL PGSLW(lnwidth)  ! thin line
          CALL PGSCI(2) ! red color
          DO I=1,LRZMAX
           RLRZAP13(I)=tr5i(I) !from currpar_starnue_n(I)-currpar_starnue0_n(I)
          ENDDO
          CALL PGLINE(lrzmax,RLRZAP1(1:lrzmax),RLRZAP13(1:lrzmax))
          write(t_,5003) tr5i(lrzmax)
          CALL PGMTXT('T',R46,R40,R40,
     &    'Red-- (sigma_coll_neo-sigma_banana)*Ephi'//t_) !YuP[2019-12-20]
        endif ! k.eq.kelecg
        
        if(k.eq.kelecg)then !this k is for general electrons
          CALL PGSLS(3) ! dash-dot-dash
          CALL PGSLW(lnwidth)  ! thin 
          CALL PGSCI(3) ! green color
          DO I=1,LRZMAX
           RLRZAP13(I)=tr6i(I) !from bscurm(i,1,kke) for electrons
          ENDDO
          CALL PGLINE(lrzmax,RLRZAP1(1:lrzmax),RLRZAP13(1:lrzmax))
          write(t_,5003) tr6i(lrzmax)
          CALL PGMTXT('T',R45,R40,R40,
     &    'Green-.- Bootstrap (fit model: bscurm())'//t_) !YuP[2019-12-20]
        endif
        
        if(k.eq.niong)then !this k is for general ions 
          CALL PGSLS(4) ! dotted ....
          CALL PGSLW(lnwidth)  ! thin 
          CALL PGSCI(4) ! blue color
          DO I=1,LRZMAX
           RLRZAP13(I)=tr7i(I) !from bscurm(i,2,kki) for ions
          ENDDO
          CALL PGLINE(lrzmax,RLRZAP1(1:lrzmax),RLRZAP13(1:lrzmax))
          write(t_,5003) tr7i(lrzmax)
          CALL PGMTXT('T',R45,R40,R40,
     &    'Blue/dotted: Bootstrap (fit model: bscurm() array)'//t_) !YuP[2019-12-20]
        endif
        
        !ALL TOGETHER, for given species k
        CALL PGSCI(1) ! black color
        CALL PGSLS(1) ! solid 
        CALL PGSLW(lnwidth*2)  ! bold line
          DO I=1,LRZMAX
           RLRZAP14(I)=tr4i(I) 
          ENDDO
        CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP14(1))
        write(t_,5003) tr4i(lrzmax)
        CALL PGMTXT('T',R43,R40,R40,
     &   'Solid/bold: From All the above together'//t_) !YuP[2019-12-20]
     
        CALL PGSLS(1) ! restore solid line
        CALL PGSLW(lnwidth) ! restore
        CALL PGSCI(1) ! black color restored
        write(t_,4050) pltvs
        CALL PGLAB(t_,'Current (A)',' ')
        !YuP[2019-12-20] Done 


c
          if (cqlpmod.eq."enabled") then

        CALL PGPAGE
        CALL PGSVP(R4P2,R4P8,R4P2,R4P6)

        CALL PGSAVE
        CALL PGSCH(R41P44)
        CALL PGMTXT('T',R46,R4P5,R4P5,
     +       'LOCAL IN s CURNT. DENS. (AMPS per CM**2)')
        CALL PGUNSA
          write(t_,3050) k,kspeci(1,k),kspeci(2,k)
        CALL PGMTXT('T',R43,R40,R40,t_)
          write(t_,3051) lrindx(1),rz(lrindx(1))/radmin
        CALL PGMTXT('T',R42,R40,R40,t_)
 




            do 551 ll=1,lsmax
              tr1s(ll)=pcurr(nch(ll),k,ll)
              tr2s(ll)=psis(ll)*pcurr(nch(1),k,1)
 551        continue
            ilsmx=lsmax
            if (sbdry .eq. "periodic") then
              ilsmx=lsmax+1
              tr1s(ilsmx)=tr1s(1)
              tr2s(ilsmx)=tr2s(1)
            endif
            fmin=0.
            fmax=0.
            call aminmx(tr1s(1),1,ilsmx,1,fmin,fmax,kmin,kmax)
            if (fmin .ge. fmax) fmin=fmax-.1*abs(fmax)-1.e+1
            call aminmx(tr2s(1),1,ilsmx,1,fmin1,fmax1,kmin,kmax)
            if (fmin1.lt.fmin) fmin=fmin1
            if(fmax1.gt.fmax) fmax=fmax1
            if(fmax.gt.0.) fmax=fmax*1.05 ! extend the upper range


        DO I=1,LSMAX
           RLRORSA(I)=z(i,lrindx(1))
           RLRORSA1(I)=tr1s(i)
           RLRORSA2(I)=tr2s(i)
        ENDDO
        RPG1=fmin
        RPG2=fmax
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
        CALL PGSWIN(RLRORSA(1),RLRORSA(LSMAX),RPG1,RPG2)
        CALL PGBOX('BCNST',R40,0,'BCNST',R40,0)
        CALL PGLINE(lsmax,RLRORSA(1),RLRORSA1(1))
        CALL PGSLS(2)
        CALL PGLINE(lsmax,RLRORSA(1),RLRORSA2(1))
        CALL PGSLS(1)
        write(t_,4051)
        CALL PGLAB(t_,'Curr Den (A/cm\u2\d',' ')



            go to 500
          endif
c
          if (nrf .ne. 0) then

        CALL PGPAGE
        CALL PGSVP(R4P2,R4P8,R4P2,R4P6)

        CALL PGSAVE
        CALL PGSCH(R41P44)
        CALL PGMTXT('T',R46,R4P5,R4P5,
     +       'LOCAL RF POWER (WATTS per CC)')
        CALL PGUNSA
          write(t_,3002) k,kspeci(1,k),kspeci(2,k),n
        CALL PGMTXT('T',R43,R40,R40,t_)

            write(t_,6002) rfpwrt(k)
 6002       format("total rf power =",1pe10.2," Watts")
        CALL PGMTXT('T',R42,R40,R40,t_)

            do 99 ll=1,lrzmax
 99         tr1(ll)=rfpwrz(k,ll)
            fmin=0.
            fmax=0.
            call aminmx(tr1(1),1,lrzmax,1,fmin,fmax,kmin,kmax)
cBH090220            if (fmin .ge. fmax) fmin=fmax-.1*abs(fmax)-1.e+1
        if (abs(fmin-fmax).lt.fmax*dgts) fmax=fmin+.001*abs(fmin)
        if(fmax.gt.0.) fmax=fmax*1.05 ! extend the upper range

        DO I=1,LRZMAX
cBH090220           RLRZAP1(I)=tr(i)
           RLRZAP1(I)=RBOUND(tr(i))
        ENDDO
cBH090220        RPG1=fmin
cBH090220        RPG2=fmax
        RPG1=RBOUND(fmin)
        RPG2=RBOUND(fmax)
        if (abs(RPG1-RPG2).le.(2.1*em33)) RPG2=RPG1+10.*em33
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
           RPG2= RPG1+1.e-16
        ENDIF
        CALL PGSWIN(RLRZAP1(1),RLRZAP1(lrzmax),RPG1,RPG2)
        CALL PGBOX('BCNST',R40,0,'BCNST',R40,0)
        DO I=1,LRZMAX
           RLRZAP11(I)=tr1(i)
        ENDDO
        CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP11(1))
        CALL PGLAB('Radius','RF Power (W/cm\u3\d)',' ')


          endif
c
          if(syncrad.ne."disabled" .and. k.eq.kelecg) then
             
            CALL PGPAGE
            CALL PGSVP(RR4P2,R4P8,R4P2,R4P6)            
            CALL PGSAVE
            CALL PGSCH(R41P44)
            CALL PGMTXT('T',R46,R4P5,R4P5,
     +           'LOCAL RF POWER (WATTS per CC)')
            CALL PGUNSA
            write(t_,3002) k,kspeci(1,k),kspeci(2,k),n
            CALL PGMTXT('T',R43,R40,R40,t_)
            write(t_,6004) psynct
 6004       format(";","synchrotron radiated power =",e16.6," Watts")
            CALL  PGMTXT('T',R42,R40,R40,t_)
            
            do 79 ll=1,lrzmax
 79            tr1(ll)=psyncz(ll)
            fmin=0.
            fmax=0.
            call aminmx(tr1(1),1,lrzmax,1,fmin,fmax,kmin,kmax)
            if (fmin .ge. fmax) fmin=fmax-.1*abs(fmax)-1.e+1
            if(fmax.gt.0.) fmax=fmax*1.05 ! extend the upper range
               
            DO I=1,LRZMAX
               RLRZAP1(I)=tr(i)
            ENDDO
            RPG1=fmin
            RPG2=fmax
            IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
               RPG2= RPG1+1.e-16
            ENDIF
            CALL PGSWIN(RLRZAP1(1),RLRZAP1(lrzmax),RPG1,RPG2)
            CALL PGBOX('BCNST',R40,0,'BCNST',R40,0)
            DO I=1,LRZMAX
               RLRZAP11(I)=tr1(i)
            ENDDO
            CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP11(1))
            CALL PGLAB('Radius','Synch Power (W/cm\u3\d)',' ')
               
          endif ! syncrad.ne."disabled" .and. k.eq.kelecg
c     
          if(sigmamod.eq."enabled" .and. pltsig.eq."enabled" 
     +                             .and. k.eq.kiong(1)) then
            do 800 lsig=1,4
               if(isigmas(lsig).eq.0) go to 800
                  
               CALL PGPAGE
               CALL PGSVP(R4P2,R4P8,R4P2,R4P5)
               CALL PGSAVE
               CALL PGSCH(R41P44)
               CALL PGMTXT('T',R48,R4P5,R4P5,
     +              'Fusion Power (Watts/cm\u3\d)')
               CALL PGUNSA
               write(t_,3002) k,kspeci(1,k),kspeci(2,k),n
               CALL PGMTXT('T',R47,R40,R40,t_)
               write(t_,8001) lsig
 8001          format("Fusion reaction number ",i3)
               CALL PGMTXT('T',R46,R40,R40,t_)
               if (lsig.eq.1) then
                  write(t_,*) "d+t=>alpha(3.5MeV)+n(14.1MeV)"
               elseif (lsig.eq.2) then
                  write(t_,*) "d+he3=>alpha(3.6MeV)+p(14.7MeV)"
               elseif (lsig.eq.3) then
                  write(t_,*) "d+d=>t(1.01Mev)+p(3.02MeV)"
               elseif (lsig.eq.4) then
                  write(t_,*) "d+d=>he3(.82MeV)+n(2.45MeV)"
               endif
               CALL PGMTXT('T',R45,R40,R40,t_)
          
               write(t_,8002) fuspwrvt(lsig)
 8002          format("fusion power =",1pe16.6," Watts")
               CALL PGMTXT('T',R43,R40,R40,t_)

cBH120314:  Have removed the isigsgv2 option and references to it, since
cBH120314:  this functionality appears to have no physical use.
cBH120314               if (isigsgv2.eq.1) then
cBH120314                  write(t_,8003) fuspwrmt(lsig)
cBH120314 8003           format("fusion power(equiv. Maxwln) =",1pe16.6," Watts")
cBH120314                  CALL PGMTXT('T',R42,R40,R40,t_)
cBH120314               endif                          
             
               do 801 ll=1,lrzmax
                  tr1(ll)=fuspwrv(lsig,ll)
                  tr2(ll)=fuspwrm(lsig,ll)
 801           continue
             
               fmin=0.
               fmax=0.
               call aminmx(tr1(1),1,lrzmax,1,fmin,fmax,kmin,kmax)
cBH120314               if (isigsgv2.eq.1) then
cBH120314                 call aminmx(tr2(1),1,lrzmax,1,fmin1,fmax1,kmin,kmax)
cBH120314                 fmin=min(fmin,fmin1)
cBH120314                 fmax=max(fmax,fmax1)
cBH120314               endif
               if (fmin .ge. fmax) fmin=fmax-.1*abs(fmax)-1.e+1
               if(fmax.gt.0.) fmax=fmax*1.05 ! extend the upper range
            
               DO I=1,LRZMAX
                  RLRZAP1(I)=tr(i)
               ENDDO
               RPG1=fmin
               RPG2=fmax
               IF ( RPG2-RPG1 .le. 1.e-16 ) THEN ! YuP [02-23-2016]
                  RPG2= RPG1+1.e-16
               ENDIF
               CALL PGSWIN(RLRZAP1(1),RLRZAP1(lrzmax),RPG1,RPG2)
               CALL PGBOX('BCNST',R40,0,'BCNST',R40,0)
               DO I=1,LRZMAX
                  RLRZAP11(I)=tr1(i)
               ENDDO
               CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP11(1))
               CALL PGSLS(2)
               DO I=1,LRZMAX
                  RLRZAP12(I)=tr2(i)
               ENDDO
               CALL PGLINE(lrzmax,RLRZAP1(1),RLRZAP12(1))
               CALL PGSLS(1)
               CALL PGSAVE
               CALL PGSCH(R41P44)
               CALL PGLAB(t_,'Fusion Power (W/cm/u3/d',' ')
               CALL PGUNSA        
c...           mnt  Generate time-dependent plot of total fusion rates
               CALL PGPAGE
               CALL PGSVP(R4P2,R4P8,R4P2,R4P5)
               CALL PGSAVE
               CALL PGSCH(R41P44)
               CALL PGMTXT('T',R48,R4P5,R4P5,
     +            'Fusion Rx Rate (/sec) Vs Time')
               CALL PGUNSA
               write(t_,3002) k,kspeci(1,k),kspeci(2,k),n
               CALL PGMTXT('T',R47,R40,R40,t_)
               write(t_,8005) lsig
 8005          format('Fusion reaction number ',i3)
               CALL PGMTXT('T',R46,R40,R40,t_)              
               write(t_,8006) sigftt(nch(1),lsig)
 8006          format("Reaction rate =",1pe16.6," /sec")
               CALL PGMTXT('T',R44,R40,R40,t_)
              
cBH120314               if (isigsgv2.eq.1) then
cBH120314                 write(t_,8007) sigmtt(nch(1),lsig)
cBH120314 8007           format("Reaction rate(equiv. Maxwln) =",1pe14.6," /sec")
cBH120314                 CALL PGMTXT('T',R42,R40,R40,t_)
cBH120314               endif
              
               call aminmx(sigftt(1,lsig),1,nch(1),1,
     ~         emin,emax,kmin,kmax)
cBH120314               if (isigsgv2.eq.1) then
cBH120314                  call aminmx(sigmtt(1,lsig),1,nch(1),1,
cBH120314     ~            emin1,emax1,kmin,kmax)
cBH120314                  emin=min(emin,emin1)
cBH120314                  emax=max(emax,emax1)
cBH120314               endif
               dgts=1.e-8
               if (abs(emin-emax).lt.emax*dgts) emax=emin+.001*abs(emin)
               if(emax.gt.0.) emax=emax*1.05 ! extend the upper range
               DO I=1,NCH(1)
                  RNONCHA1(I)=ptime(i,1)
               ENDDO
               RPG1=emin
               RPG2=emax
               IF ( RPG2-RPG1 .le. 1.e-16 ) THEN
                   RPG2= RPG1+1.e-16
               ENDIF
               CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(1)),RPG1,RPG2)
               CALL PGBOX('BCNST',R40,0,'BCNST',R40,0)
               CALL PGLAB('Time (sec)','Fusion Rx Rate (/sec)',' ')
               DO I=1,NCH(1)
                  RNONCHA2(I)=sigftt(I,lsig)
               ENDDO
               CALL PGLINE(nch(1),RNONCHA1,RNONCHA2)
cBH120314               if (isigsgv2.eq.1) then
cBH120314                  DO I=1,NCH(1)
cBH120314                     RNONCHA2(I)=sigmtt(I,lsig)
cBH120314                  ENDDO
cBH120314                  CALL PGLINE(nch(1),RNONCHA1,RNONCHA2)
cBH120314               endif        
        
 800        continue ! lsig=1,4
          endif !if(sigmamod="enabled" .and. pltsig="enabled" .and. k=kiong(1)
          
 500    continue ! k




        call tdpltjop



c     endif n>0, before 500 loop
      endif

      return
      end subroutine tdpltmne
