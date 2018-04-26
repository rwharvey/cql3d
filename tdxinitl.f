c
c
      subroutine tdxinitl
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c..................................................................
c     This routine "fills in" input data parabolically. It also
c     computes the normalized radial mesh.
c..................................................................
      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE     

      dimension titmp(250),detmp(250)

c.......................................................................


c..................................................................
c     generate radial (rho) mesh. rya will be the normalized mesh.
c     rrz will be the intermediate (mid) mesh points. Later rz
c     will be defined to be the non-normalized actual radial
c     toroidal flux mesh.
c..................................................................

      if (rzset .eq. "enabled") then

c..................................................................
c     use prespecified  rya() radial mesh in input deck...
c..................................................................

        rya(0)=0.
        rrz(0)=0.
        do 20 ll=1,lrzmax-1
          if (rya(ll) .lt. rya(ll-1)) then
            call tdwrng(3)
          endif
          rrz(ll)=(rya(ll+1)+rya(ll))*.5
 20     continue
        rrz(lrzmax)=1.
        go to 21

      elseif (roveram.gt.1.e-5) then ! (and rzset='disabled' case)
            ! default is 1.e-6

c..................................................................
c     code determines geometric mesh with smallest inverse aspect
c     ratio mesh point prespecified.
c..................................................................

        roveramp=roveram
        hrz=(1.d0-roveramp)/(lrzmax-1)
        hrz1=hrz*rfacz
        rrz(0)=0.
        rrz(1)=roveramp
        rrz(2)=rrz(1)+hrz1
        rrmax= 1.d0-roveramp
        call micgetr(lrzmax,rrmax,hrz1,ram,ksingul)
        do 2 ll=3,lrzmax
          rrz(ll)=rrz(ll-1)+ram*(rrz(ll-1)-rrz(ll-2))
 2      continue
  
      else  ! (roveram.le.1.e-5 and rzset='disabled' case)

c..................................................................
c     code generates radial (rho) mesh geometrically with no minimum
c     mesh point prespecified.
c..................................................................


        hrz=1./lrzmax
        hrz1=hrz*rfacz
        rrz(0)=0.
        rrz(1)=hrz1
        rrmax=one
        call micgetr(lrzmax+1,rrmax,hrz1,ram,ksingul)
        if (ksingul.eq.1) stop 'tdxinitl: problem with radial mesh'
        do 1 ll=2,lrzmax
          rrz(ll)=rrz(ll-1)+ram*(rrz(ll-1)-rrz(ll-2))
 1      continue

      endif  ! On rzset.eq.'enabled' and two cases of roveram

      rya(0)=0.
      do 6 ll=1,lrzmax
        rya(ll)=(rrz(ll-1)+rrz(ll))*.5
 6    continue
      if (roveram.gt.1.e-5) rya(1)=roveram
      ! Note: this reset of rya(1) may result in non-uniform rya grid
      ! at first three points.
      ! To be sure rya grid is uniform, set roveram to a very small value
      ! like roveram=1.e-8 (any value lt.1.e-5).
      ! Example: For roveram=1.e-8, lrz=50, rfacz=1.0 (rzset='disabled')
      ! rya grid is uniform in range [0.01;0.99], with step = 0.02.

CMPIINSERT_IF_RANK_EQ_0      
      WRITE(*,'(a,2e13.4)')'tdxinitl: rrmax,ram=',rrmax,ram
      do ll=1,lrzmax
      WRITE(*,'(a,i6,2e13.4)')'tdxinitl: ll,rrz,rya=',ll,rrz(ll),rya(ll)
      enddo
CMPIINSERT_ENDIF_RANK

 21   continue

c.......................................................................
c     determine array rovera
c.......................................................................

      do 30 ll=1,lrzmax
        rovera(ll)=rya(ll)
        if (0.lt.rovera(ll) .and. rovera(ll).lt.1.e-8) 
     +    rovera(ll)=1.e-8
 30   continue


c..................................................................
c     fill in input arrays between rya=0. and rya=1.
c     This determines density, temperature, source, etc profiles as
c     functions of the normalized radial coordinate rho. If density
c     and temperatures are to be specified as functions of psi
c     poloidal flux coordinate these arrays will be overwritten
c     in subroutine eqcoord
c     If eqsource="tsc" (running with tsc):
c     iprone=iprote=iproti="disabled" (set in tdtscinp).
c..................................................................

c..................................................................
c     If running with TSC code...
c..................................................................


      if (eqsource.eq."tsc") then
        call tdtscinp
        lrzp=lrzmax+1
        melec=2
        if (colmodl.eq.0 .or. colmodl.eq.2) melec=1
        !YuP[05-2017] Looks like in this case only ngen=1 is allowed
        call tdinterp("free","free",rho_,anecc,npsitm,rya(0),tr(0),lrzp)
        call tdinterp("free","free",rho_,tekev,npsitm,rya(0),tr1(0),
     1    lrzp)
        call tdinterp("free","free",rho_,elecf,npsitm,rya(0),tr2(0),
     1    lrzp)
        ngen=1
        if (melec.eq.2) nmax=nspc+1
        if (melec.eq.1) nmax=nspc
        ntotal=nmax+ngen
        kspeci(1,1)="e"
        kspeci(2,1)="general"
        fmass(1)=9.1095e-28
        bnumb(1)=-1.
        if (melec.eq.2) then
          nmax=nspc+1
          kspeci(1,2)="e"
          kspeci(2,2)="maxwl"
          fmass(2)=9.1095e-28
          bnumb(2)=-1.
        endif
        do 52 m=1,melec
          do 50 l=0,lrzmax
            reden(m,l)=tr(l)
            temp(m,l)=tr1(l)
            elecfld(l)=tr2(l)
 50       continue
 52     continue

c     If nspc.eq.1 then present tscout file does not specify the
c     impurity charge.  For the time being, specify the nspc=1-case
c     throught the cqlinput for species melec+1.

        if (nspc.ne.1)  then
          do 54 k=1,nspc
            kspeci(1,melec+k)="ion"
            kspeci(2,melec+k)="maxwell"
            fmass(melec+k)=amass(k)
            bnumb(melec+k)=achrg(k)
 54       continue
        endif  ! On nspc.ne.1
        do 53 k=melec+1,melec+1+nspc
          ku=k-melec
          do 55 l1=1,npsitm
            titmp(l1)=tikev(ku,l1)
            detmp(l1)=anicc(ku,l1)
 55       continue
          call tdinterp("free","free",rho_,titmp,npsitm,rya(0),tr1(0),
     1      lrzp)
          if (nspc.ne.1)  then
            call tdinterp("free","free",rho_,detmp,npsitm,rya(0)
     1        ,tr(0),lrzp)
          else
            do 56  j=0,lrzmax
 56         reden(melec+1,j)=reden(1,j)/bnumb(melec+1)
          endif
          do 51 l=0,lrzmax
            reden(k,l)=tr(l)
            temp(k,l)=tr1(l)
 51       continue
 53     continue
      endif  ! On eqsource.eq."tsc"

c......................................................................
c     Profiles can be entered with "parabola", "spline", 
c     or "asdex" option.
c     For iprozeff.ne."disabled", there are only two different ion Z's.
c......................................................................

c     density

      do 11  k=1,ntotal

c     If iprozeff.ne."disabled", only set electron profile here.
         if (iprozeff.ne."disabled" .and. 
     +        (k.ne.kelecg .and. k.ne.kelecm)) go to 11
         
         if (iprone.eq."parabola")  then

            call tdxin13d(reden,rya,lrzmax,ntotala,k,npwr(0),mpwr(0))
            
         elseif (iprone.eq."spline")  then

            if (enein(1,k).ne.zero) then

               call tdinterp("zero","free",ryain,enein(1,k),njene,
     +              rya(1),tr(1),lrzmax)
               tr(0)=enein(1,k)
               do 13 ll=0,lrzmax
                  reden(k,ll)=tr(ll)
 13            continue
            else
cBH131029: Following do only correct for tr already set when
cBH131029: k.eq.kelecg??
               write(*,*) 'tdxinitl: Check coding here'
               do 9  ll=0,lrzmax
                  reden(k,ll)=tr(ll)/abs(bnumb(k))
 9             continue
            endif
            
         elseif (iprone.eq."asdex")  then
            
            do 141  ll=0,lrzmax
               reden(k,ll)=1.e14*tdpro(rya(ll),radmin/100.,acoefne)
     +              /abs(bnumb(k))
 141        continue
            
         endif
         
 11   continue

c......................................................................
c     Finish up with zeff and ions if iprozeff.ne."disabled"
c
c     If have 2 ion species with Zi braketing Zeff, with charge 
c     neutrality, we use
c     n1=ne*(Zeff-Z2)/(Z1*(Z1-Z2)), and similarly for n2.
c     Apportion the ion density, if have .ge.2 ions with same Z.
c......................................................................
      if (iprozeff.ne."disabled") then  !endif on line 432

         if (iprozeff.eq."parabola") then
            dratio=zeffin(1)/zeffin(0)
            do 12 ll=1,lrzmax
               call profaxis(rn,npwrzeff,mpwrzeff,dratio,rya(ll))
cBH131029:  Should be zeffscal multiplier, consistent with spline case??
               zeff(ll)=zeffin(0)*rn
 12         continue
         elseif (iprozeff.eq."spline") then
            do ij=1,njene
               zeffin(ij)=zeffscal*zeffin(ij)
            enddo
            call  tdinterp("zero","free",ryain,zeffin(1),njene,rya(1),
     +           zeff(1),lrzmax)
         elseif (iprozeff.eq."prbola-t") then
            do it=1,nbctime
               zeffc(it)=zeffscal*zeffc(it)
               zeffb(it)=zeffscal*zeffb(it)
            enddo
         elseif (iprozeff.eq."prbola-t".or.iprozeff.eq."spline-t") then
            !scaled in profiles
	    continue
         endif
         
CMPIINSERT_IF_RANK_EQ_0      
         WRITE(*,*)'tdxinitl, zeff(1:lrzmax)= ',(zeff(ll),ll=1,lrzmax)
         !pause
CMPIINSERT_ENDIF_RANK
         
c     Check that range of bnumb for Maxl species brackets zeff
c     Skip if iprozeff=prbola-t or spline-t, since profiles called later
         if (iprozeff.ne."prbola-t" .and. iprozeff.ne."spline-t") then
         fmaxx=0.
         fminn=1000.
         do k=1,nionm
CMPIINSERT_IF_RANK_EQ_0      
            WRITE(*,*)'k,kionm(k),bnumb(kionm(k))=',
     +                 k,kionm(k),bnumb(kionm(k))
CMPIINSERT_ENDIF_RANK
            fmaxx=max(fmaxx,bnumb(kionm(k)))
            fminn=min(fminn,bnumb(kionm(k)))
         enddo
         !YuP[03-2017] For iprozeff.eq."spline", zeff() is found from spline
         !             using zeffin
         do 121 ll=1,lrzmax
            if(zeff(ll).gt.fmaxx .or. zeff(ll).lt.fminn) then
CMPIINSERT_IF_RANK_EQ_0      
               WRITE(*,*)
     +          'tdxinitl, max/min of bnumb(kionm(k)), fmaxx,fminn=',
     +           fmaxx,fminn
               WRITE(*,*)'tdxinitl.f: ',
     +              'Adjust bnumb(kion) for compatibility with zeff'
CMPIINSERT_ENDIF_RANK
               stop
            endif
 121      continue         
          endif  !On iprozeff
c         write(*,*)'tdxinitl, fmaxx,fminn= ',fmaxx,fminn
         
c     Check number of ion Maxwl species with different bnumb
c     (Need at least two in order to fit zeff. Check ainsetva.f.)
         if (nionm.lt.2) stop 'tdxinitl: ion species problem'
         ndif_bnumb=1
         do k=2,nionm
            if (abs(bnumb(kionm(k))/bnumb(kionm(1))-1.).gt.0.01) 
     +           ndif_bnumb=ndif_bnumb+1
         enddo

CMPIINSERT_IF_RANK_EQ_0      
         WRITE(*,*)
         WRITE(*,*)'tdxinitl: Number Maxl ion species w diffrnt bnumb'
         WRITE(*,*)'tdxinitl: ndif_bnumb= ',ndif_bnumb
         WRITE(*,*)
CMPIINSERT_ENDIF_RANK

c     Interpolate input ion densities onto rya grid
         if (iprone.eq."parabola")  then
            
            do kk=1,nionm
               k=kionm(kk)
               call tdxin13d(reden,rya,lrzmax,ntotala,k,npwr(0),mpwr(0))
c               write(*,*)'tdxinitl: k,reden(k,1:lrzmax)=',
c     +                              k,(reden(k,i),i=1,lrzmax)
            enddo
            
         elseif (iprone.eq."spline")  then
            
cBH080918            do kk=1,ndif_bnumb
            do kk=1,nionm
               k=kionm(kk)
               if (enein(1,k).ne.zero) then
c$$$                  do ij=1,njene
c$$$                     enein(ij,k)=enescal*enein(ij,k)
c$$$                  enddo
                  call tdinterp("zero","free",ryain,enein(1,k),njene,
     +                 rya(1),tr(1),lrzmax)
                  tr(0)=enein(1,k)
                  do ll=0,lrzmax
                     reden(k,ll)=tr(ll)
                  enddo
               else
                  do ll=0,lrzmax
                     reden(k,ll)=tr(ll)/abs(bnumb(k))
                  enddo
               endif
            enddo
         endif



c      write(*,*)'kionm(1:2)= ',kionm(1),kionm(2)
c      write(*,*)'reden(kionm(1),ll)= ',(reden(kionm(1),ll),ll=0,lrzmax)
c      write(*,*)'reden(kionm(2),ll)= ',(reden(kionm(2),ll),ll=0,lrzmax)
         
c     If all ion species have different bnumb, then set ion densities
c     to 1, compatible with following situation of multiple ion species
c     with the same bnumb().
c     They will be recalculated below from zeff and electron density.
c
c     Save ratios of densities for Maxwl ions with same charge in reden
c     to be used in following ion density calc from zeff.
c     (First nsame_bnumb ions have same charge.  These are ion species
c     at the head of the ion species list.)

      if (nionm.eq.ndif_bnumb) then
         nsame_bnumb=0               !i.e., =0 if all diff bnumb ions
      else
         nsame_bnumb=nionm-ndif_bnumb+1 ! .ge.2, for 2 or more the same
      endif
      if (nsame_bnumb.gt.0) then  !i.e., 2 or more ions have same bnumb.
                                  !Densities will contain density ratios.
                                  !Equal-bnumb() species at beginning of
                                  !Maxl ions.
c        Renormalizing the density ratios as fractions for equal-bnumb:
            do ll=1,lrzmax
               dsum=0.
               do k=1,nsame_bnumb
                  dsum=dsum+reden(kionm(k),ll)
               enddo
               do k=1,nsame_bnumb
                  reden(kionm(k),ll)=reden(kionm(k),ll)/dsum
               enddo
            enddo

c      write(*,*)
c      write(*,*)'After ratios stored for the ndif_bnumb Maxl ions:'
c      write(*,*)'reden(kionm(1),ll)= ',(reden(kionm(1),ll),ll=1,lrzmax)
c      write(*,*)'reden(kionm(2),ll)= ',(reden(kionm(2),ll),ll=1,lrzmax)
         
      endif  ! on nsame_bnumb.gt.0

c     Set rest of ion densities to 1.
      do kk=nsame_bnumb+1,nionm
         k=kionm(kk)
         do l=0,lrzmax
            reden(k,l)=one
         enddo
      enddo
         
        
         
c     NOTE: reden(k, ) on rhs of following reden(,) formulas is
c     density ratio, .le.1. for diff species with same bnumb.
      if (nsame_bnumb.eq.0) then
         nsame_bnumb1=nsame_bnumb
      else         
         nsame_bnumb1=nsame_bnumb-1
      endif
      do 14 k=kionm(1),kionm(nionm)
c     write(*,*)'tdxinitl, do 14 loop, k= ',k
c        For k pointing to equal-bnumb species (2 or more)
c        k1,k2 indicate bnumb values to use
         if (k.le.(kionm(1)+nsame_bnumb1)) then
            k1=k
            k2=kionm(nionm)
c        For k beyond equal-bnumb species
         else
            k1=k
            k2=kionm(1)
         endif
         if (iprone.eq."parabola") then
            reden(k,0)=reden(kelec,0)*reden(k,0)*(zeff(1)
     +           -bnumb(k2))/(bnumb(k1)-bnumb(k2))/bnumb(k1)
            do 142 ll=1,lrzmax
               reden(k,ll)=reden(kelec,ll)*reden(k,ll)*(zeff(ll)
     +              -bnumb(k2))/(bnumb(k1)-bnumb(k2))/bnumb(k1)
 142        continue
            
         elseif(iprone.eq."spline") then
            do ll=1,lrzmax
               reden(k,ll)=reden(kelec,ll)*reden(k,ll)*(zeff(ll)
     +              -bnumb(k2))/(bnumb(k1)-bnumb(k2))/bnumb(k1)
            enddo
         endif
 14   continue
      
      
c      write(*,*)
c      write(*,*)'After calc of densities for Maxl ions:'
c      write(*,*)'reden(kionm(1),ll)= ',(reden(kionm(1),ll),ll=1,lrzmax)
c      write(*,*)'reden(kionm(2),ll)= ',(reden(kionm(2),ll),ll=1,lrzmax)
c      if (nionm.gt.2) then
c      write(*,*)'reden(kionm(3),ll)= ',(reden(kionm(3),ll),ll=1,lrzmax)
c      endif
         
c     Copy Maxwellian ion densities to any ion general species
c     in order of the indexes.

      do k=1,niong
         if (nionm.ge.niong) then
            do l=0,lrzmax
               reden(kiong(k),l)=reden(kionm(k),l)
            enddo
         endif
      enddo
         
         
      endif  ! on iprozeff.eq."parabola" .or. iprozeff.eq."spline"
      
c     Renormalize densities using enescal
      do k=1,ntotal
         do l=0,lrzmax
            reden(k,l)=enescal*reden(k,l)
         enddo
      enddo

c     Diagnostic printout of general species densities:
c      write(*,*)'After setting general species densities:'
c      do k=1,ngen
c         if (k.eq.kelecg) then
c         write(*,*)'reden(kelecg,ll)= ',(reden(kelecg,ll),ll=1,lrzmax)
c         endif
c         if (k.eq.kiongg(k)) then
c         write(*,*)'reden(kiongg(k),ll)= ',
c     +        (reden(kiongg(k),ll),ll=1,lrzmax)
c         endif
c      enddo



c     electron temperature

      if (iprote.eq."parabola")  then
        do 15  k=1,ntotal
          if (bnumb(k).eq.-1.)
     +      call tdxin13d(temp,rya,lrzmax,ntotala,k,npwr(k),mpwr(k))
 15     continue
      elseif (iprote.eq."spline")  then
c100126         do ij=1,njene
c100126            tein(ij)=tescal*tein(ij)
c100126         enddo
         do 16  k=1,ntotal
            if(bnumb(k).eq.-1.)  then
               call tdinterp("zero","free",ryain,tein,njene,rya(1),
     +              tr(1),lrzmax)
               tr(0)=tein(1)
               do 19  ll=0,lrzmax
                  temp(k,ll)=tr(ll)
 19            continue
            endif
 16      continue
      elseif (iprote.eq."asdex")  then
         do 17  k=1,ntotal
            if(bnumb(k).eq.-1.)  then
               do 18  ll=0,lrzmax
                  temp(k,ll)=1.e-3*tdpro(rya(ll),radmin/100.,acoefte)
 18            continue
            endif
 17      continue
      endif


c     ions temperature

      if (iproti.eq."parabola")  then
         do 25  k=1,ntotal
            if (bnumb(k).ne.-1.) call tdxin13d(
     +           temp,rya,lrzmax,ntotala,k,npwr(k),mpwr(k))
 25      continue
      elseif (iproti.eq."spline")  then
c100126         do ij=1,njene
c100126            tiin(ij)=tiscal*tiin(ij)
c100126         enddo
         do 26  k=1,ntotal
            if(bnumb(k).ne.-1.)  then
               call tdinterp("zero","free",ryain,tiin,njene,rya(1),
     +              tr(1),lrzmax)
               tr(0)=tiin(1)
               do 29  ll=0,lrzmax
                  temp(k,ll)=tr(ll)
 29            continue
            endif
 26      continue
      elseif (iproti.eq."asdex")  then
         do 27  k=1,ntotal
            if(bnumb(k).ne.-1.)  then
               do 28  ll=0,lrzmax
                  temp(k,ll)=1.e-3*tdpro(rya(ll),radmin/100.,acoefte)
 28            continue
            endif
 27      continue
      endif

c100126
c     Renormalize temperatures using tescal/tiscal
      do k=1,ntotal
         if (bnumb(k).eq.-1.) then
            do l=0,lrzmax
               temp(k,l)=tescal*temp(k,l)
            enddo
         else
            do l=0,lrzmax
               temp(k,l)=tiscal*temp(k,l)
            enddo
         endif
      enddo
      
c     vphipl, toroidal plasma velocity for use in freya
      call bcast(vphipl,zero,lrzmax)
      if (iprovphi.ne."disabled") then
         if (iprovphi.eq."parabola") then
            dratio=vphiplin(1)/vphiplin(0)
            do 60 ll=1,lrzmax
               call profaxis(rn,npwrvphi,mpwrvphi,dratio,rya(ll))
cBH120701:  Added in vphiscal factor.
               vphipl(ll)=vphiscal*vphiplin(0)*rn
 60         continue
         elseif (iprovphi.eq."spline") then
            do ij=1,njene
               vphiplin(ij)=vphiscal*vphiplin(ij)
            enddo
            call  tdinterp("zero","free",ryain,vphiplin(1),njene,
     +           rya(1),vphipl(1),lrzmax)
         endif
      endif

c     neutral density
      
      call bcast(enn,zero,lrza*npaproca)
      if (ipronn.ne."disabled") then
         do kkk=1,npaproc
         if (npa_process(kkk).ne.'notset'.and.kkk.ne.5) then
         if (ipronn.eq."exp") then
            do ll=1,lrzmax
               enn(ll,kkk)=ennb(kkk)*exp((rya(ll)-1.)*radmin/ennl(kkk))
            enddo
         elseif (ipronn.eq."spline") then
            call  tdinterp("zero","free",ryain,ennin(1,kkk),njene,
     +           rya(1),enn(1,kkk),lrzmax)
         endif  ! On ipronn
         endif  ! On npa_process
         enddo  ! kkk
c        Use electron density for radiation recombination:
         if (npa_process(5).eq.'radrecom') then
            kkk=max(kelecg,kelecm)  !I.e., using bkgrnd distn, if avail.
            do ll=1,lrzmax
               enn(ll,5)=reden(kkk,ll)
            enddo
         endif
         do kkk=1,npaproc 
            do ll=1,lrzmax
               enn(ll,kkk)=ennscal(kkk)*enn(ll,kkk)
            enddo
         enddo
      endif

c     tauegyz, eparc, eperc   ! YuP: tauegyz or tauegy?
      do 3 k=1,ngen
         if(tauegy(k,0).gt.0.) ! YuP-101220: changed tauegyz to tauegy !
     1        call tdxin13d(tauegy,rya,lrzmax,ngena,k,negy(k),megy(k))
         if(torloss(k).ne."disabled")  then
            call tdxin13d(eparc,rya,lrzmax,ngena,k,ntorloss(k),
     +           mtorloss(k))
            call tdxin13d(eperc,rya,lrzmax,ngena,k,ntorloss(k),
     +           mtorloss(k))
         endif
         
c     arrays for source term
         nga=ngena
         do 4 kk=1,nso
            call tdxin23d(sellm1z,rya,lrzmax,nga,nsoa,k,kk,npwrsou(k),
     +           mpwrsou(k))
            call tdxin23d(sellm2z,rya,lrzmax,nga,nsoa,k,kk,npwrsou(k),
     +           mpwrsou(k))
            call tdxin23d(szm1z,rya,lrzmax,nga,nsoa,k,kk,npwrsou(k),
     +           mpwrsou(k))
            call tdxin23d(seppm1z,rya,lrzmax,nga,nsoa,k,kk,npwrsou(k),
     +           mpwrsou(k))
            call tdxin23d(seppm2z,rya,lrzmax,nga,nsoa,k,kk,npwrsou(k),
     +           mpwrsou(k))
            call tdxin23d(szm2z,rya,lrzmax,nga,nsoa,k,kk,npwrsou(k),
     +           mpwrsou(k))
            call tdxin23d(sem1z,rya,lrzmax,nga,nsoa,k,kk,npwrsou(k),
     +           mpwrsou(k))
            call tdxin23d(sem2z,rya,lrzmax,nga,nsoa,k,kk,npwrsou(k),
     +           mpwrsou(k))
            call tdxin23d(sthm1z,rya,lrzmax,nga,nsoa,k,kk,npwrsou(k),
     +           mpwrsou(k))
            call tdxin23d(scm2z,rya,lrzmax,nga,nsoa,k,kk,npwrsou(k),
     +           mpwrsou(k))
            
c     Allow for possiblity to intialize asorz(k,kk,1:lrzmax) directly
            if(asorz(k,kk,0).ne.-1.) then
               call tdxin23d(asorz,rya,lrzmax,nga,nsoa,k,kk,npwrsou(0),
     +              mpwrsou(0))
            endif
            
            do 5  ll=1,lrzmax
               asor(k,kk,ll)=asorz(k,kk,ll)
 5          continue
 4       continue
 3    continue

c     electric field profile
      if (eqsource.ne."tsc") then
         if (iproelec.eq."parabola") then
            elecfldc=elecfld(0)
            elecfldb=elecfld(1)
            call tdxin33d(elecfld,rya,lrzmax,npwrelec,mpwrelec)
c           central electric field elecfld(0) is an input.            
         elseif (iproelec.eq."spline") then
c100126            do ij=1,njene
c100126               elecin(ij)=elecscal*elecin(ij)
c100126            enddo
            call tdinterp("zero","free",ryain,elecin,njene,rya(1),
     +           elecfld(1),lrzmax)
cBH171124: for iproelec=spline-t or prbola-t, elecfld set in profiles
         elseif (iproelec.eq."spline-t".or.iproelec.eq."prbola-t") then
            continue
         else
            stop '   iproelec Problem'
         endif
c     Scale electric field profile
         do ij=0,lrzmax
            elecfld(ij)=elecscal*elecfld(ij)
         enddo
         elecfldc=elecscal*elecfldc
         elecfldb=elecscal*elecfldb

      endif

      return
      end



