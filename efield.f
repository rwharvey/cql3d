c
c
      subroutine efield
      implicit integer (i-n), real*8 (a-h,o-z)
      save
      include 'param.h'
      include 'comm.h'

      character*8 elecset

c..................................................................
c     CQL3D mode:
c     At the n.ge.0 time-step:
c        Calculates sptzr, tauee, taueeh, starnue=nue_eff/nue_bounce
c     At n.ge.0:
c        Calculates elecfld(lr_) according to various methods
c        specified by efswtch and efswtchn.
c        sptzr is recalculated with Z.ne.1 corrections.
c..................................................................

      elecset="eoved"
      if (eoved.eq.zero) elecset="elecfld"
c      if (lrzmax.gt.1)   elecset="elecfld"
      call cfpgamma

c.......................................................................
c     Sets Z=1 sptzr resistivity value according to time step n=0 
c     parameters, to 1/sigma_parallel in Book's NRL plasma formulary.
c.......................................................................

      if (n .eq. 0) then
        if (cqlpmod .ne. "enabled") then
c     
c     formula should be 0.5064*N(Z)*me/taueeh/ne/e**2, with N(Z=1)=1
c     Thus one should use game(e,e) as in taueeh, but as it is not
c     clear which to use we take the average between game(e,e) and game(e,i)
          sptzr(l_)=0.5064*4.*sqrt(2.*pi)/3.*charge**2
     *      *0.5*(gama(kelec,kelec)+gama(kelec,kionn))
     *      *sqrt(fmass(kelec))/(energy(kelec,lr_)/1.5*ergtkev)**1.5
        else
          sptzr(l_)=0.5064*4.*sqrt(2.*pi)/3.*charge**2
     *      *0.5*(gama(kelec,kelec)+gama(kelec,kionn))
     *      *sqrt(fmass(kelec))/(enrgypa(kelec,ls_)/1.5*ergtkev)**1.5
c
c     divide by ne*Zeff = ni*Zeff**2
c     Hinton-Hazeltine p.270 and 297 (same as ONETWO 4.2-20, 4.2-30)
          taueeh(ls_)=vthpar(kelec,ls_)**3*fmass(kelec)**2
     1      /(4.*sqrt(2.*pi)*denpar(kelec,ls_)*charge**4*
     1      gama(kelec,kelec))*3./zeff(lr_)
          starnue(ls_)=rgeom(lr_)*bmod0(lr_)/bthr0(lr_)/
     1      vthpar(kelec,ls_)/taueeh(ls_)/eps(lr_)**1.5
        endif

      endif

c     The remainder of the subroutine computes only radial quantities
c     Thus done only when l_=lmdpln_ (as is the case in cql3d operation).

      if (l_ .ne. lmdpln_) return

c..................................................................
c     Use Connor formula to compute resistivity - should
c     be accurate for most values of invers aspect ratio.
c..................................................................

      call restcon

c..................................................................
c     Reset Z=1 Spitzer resistivity, vth,tauee,taueeh,elecr for n.ne.0
c..................................................................


c..................................................................
c     vth is the thermal velocity = sqrt(T/m) (at t=0 defined in ainpla)
c     But, T==temp(k,lr) can be changed in profiles.f, 
c     in case of iprote (or iproti) equal to "prbola-t" or "spline-t"
c..................................................................

      do 10 k=1,ngen
         vth(k,lr_)=(temp(k,lr_)*ergtkev/fmass(k))**.5
 10   continue

c..................................................................
c     tauee(lr_) is the electron-electron slowing down time at midplane
c     tauee = 0.532/nu^{ee}_s, where nu^{ee}_s is the low velocity
c     expression for ee-slowing in Book's NRL Plasma Formulary p.32,
c     0.532=2**(3/2)/(3*pi**0.5). As nu(ee)_s = 2/taueeh/Zeff, we
c     have tauee = 2**(3/2)/3/sqrt(pi)/2*taueeh*Zeff, thus
c     taueeh = (3*sqrt(pi/2)/Zeff)*tauee. 
c     Note taueeh = taue(NRL p.33)=2.9E-06*...
c..................................................................

      tauee(lr_)=vth(kelec,lr_)**3*fmass(kelec)**2
     1  /(4.*pi*reden(kelec,lr_)*charge**4*
     1  gama(kelec,kelec))
c     Hinton-Hazeltine definition, Eq.(5.4), for e-i collisions.
c     1./taueeh = 4/3 sqrt(2*pi) ne*zeff*charge**4*gamma/(sqrt(me)*Te**1.5)
c     Thus: taueeh=3 sqrt(pi/2) / Zeff * tauee
c     Note: taueeh and nuestar same as ONETWO 4.2-20 and 4.2-30
c     (Eq. 4.2-28 and 4.2-39, in GA-A16178, Pfeiffer, Davidson, Miller, Waltz)
      if (cqlpmod .ne. "enabled") then
        taueeh(lr_)=tauee(lr_)*3.*sqrt(pi/2.)/zeff(lr_)
        starnue(lr_)=rgeom(lr_)*bmod0(lr_)/bthr0(lr_)/
     1    vth(kelec,lr_)/taueeh(lr_)/eps(lr_)**1.5
c     
c     formula should be 0.5064*N(Z)*me/taueeh/ne/e**2, with N(Z=1)=1
c     Thus one should use game(e,e) as in taueeh, but as it is not
c     clear which to use we take the average between game(e,e) and game(e,i)
c     This expression for sptzr is given by Eq. 4.2-77 if ONETWO manual,
c     GA-A16178, the same as Eq. 5.66 of Hinton-Hazeltine, Rev. Mod. Phys.
          sptzr(l_)=zeff(lr_)*(0.29+0.46/(1.08+zeff(lr_)))
     *      *4.*sqrt(2.*pi)/3.*charge**2
     *      *0.5*(gama(kelec,kelec)+gama(kelec,kionn))
     *      *sqrt(fmass(kelec))/(temp(kelec,lr_)*ergtkev)**1.5
      endif

 2    continue

c..................................................................
c     elecr is the Dreicer electric field (as in Kulsrud et al.),
c     converted to volts/cm.
c..................................................................

      elecr(lr_)=300.*fmass(kelec)*vth(kelec,lr_)/(2.*charge*tauee(lr_))

c..................................................................
c     rovsf and rovscf are two small-epsilon approximations to
c     rovs(lr_) (see below)
c..................................................................

      rovsf=1./(1.-1.95*eps(lr_)**.5+.95*eps(lr_))
      rovscf=1./(1.-2.09*eps(lr_)**.5)

c..................................................................
c     Toroidal electric field:
c     In the following,currxj(lr_) is a target current density (A/cm^2).
c                      currtp(lr_) is flux surface area average parallel
c                                   current density, calc'd in diaggnde.
c                                   This quantity is depracated as a
c                                   nonphysical quantity, and should be
c                                   replaced (BH, 010225).
c                      curra(kelec,l_) is flux surface area average 
c                                   parallel current density due to
c                                   runaway electrons (momenturm .gt.
c                                   minimum (ucrit, 3.*clight) 
c                                   (see subs diaggnde and soucrit).
c..................................................................

c     Calculation of plasma current for use in elecfld determination
c     (efswtch.ne."method1").

      if (efswtchn.eq."neo_hh") then
         call restcon
         call resthks(l_,lr_,lmdpln_,
     +        zreshin(lr_),zreskim(lr_),zressau1,zressau2)
         currpar(lr_)=(curra(kelec,l_)+
     +        elecfld(lr_)/300./(zreshin(lr_)*sptzr(l_)))/3.e9
      else
         currpar(lr_)=currtp(lr_)/3.e9
      endif

      if (efswtch.eq."method1") then

c        if input variable eoved (E over E-Dreicer) is .ne. 0.
c        then compute the value of elecfld(lr_)
c        (For lbdry(kelec).eq."fixed", only do it for n=0).

         if (kelecg.ne.0) then !=> kelec=kelecg (see ainspec.f)
            if (lbdry(kelecg).eq."fixed" .and. n.gt.0)  go to 31
         endif
         if (elecset.eq."eoved") elecfld(lr_)=eoved*elecr(lr_)
 31      continue

      elseif (efswtch.eq."method2") then


c        If time step .lt. noncntrl, use specified electric
c          field. Then relaxation towards specified current.
 
         if (n.lt.noncntrl) then

            if (kelecg.ne.0) then !=> kelec=kelecg (see ainspec.f)
               if (lbdry(kelecg).eq."fixed" .and. n.gt.0) go to 32
            endif
            if (elecset.eq."eoved") elecfld(lr_)=eoved*elecr(lr_)
 32         continue
            
         else
            
c           Relaxation towards target parallel currxj:
            
            currxj0(lr_)=currxj(lr_)
            elecfld(lr_)=elecfld(lr_)*(1.-efrelax*
c990131     +         (currpar(lr_) -currxj(lr_))/(sign(1.,currxj(lr_))
     +           (currpar(lr_) -currxj(lr_))/(sign(one,currxj(lr_))
c990131     +           *amax1(abs(currpar(lr_)),abs(currxj(lr_)))))
     +           *max(abs(currpar(lr_)),abs(currxj(lr_)))))
            
         endif
         
      elseif (efswtch.eq."method3") then
         

c        If time step .lt. noncntrl, use specified electric
c          field. Then relaxation towards current obtained
c          up to time noncntrl (which is saved in currxj0).
 
         if (n.lt.noncntrl) then

            currxj0(lr_)=currtp(lr_)/3.e9  
            if (kelecg.ne.0) then !=> kelec=kelecg (see ainspec.f)
               if (lbdry(kelecg).eq."fixed" .and. n.gt.0) go to 33
            endif
            if (elecset.eq."eoved") elecfld(lr_)=eoved*elecr(lr_)
 33         continue
         elseif (n.eq.noncntrl) then

            currxj0(lr_)=currpar(lr_)
            elecfld(lr_)=elecfld(lr_)
            
         else
            
c           Relaxation towards saved target parallel currxj0:

            elecfld(lr_)=elecfld(lr_)*(1.-efrelax*
     +           (currpar(lr_) -currxj0(lr_))/(sign(one,currxj0(lr_))
c990131     +           *amax1(abs(currpar(lr_)),abs(currxj0(lr_)))))
     +           *max(abs(currpar(lr_)),abs(currxj0(lr_)))))

         endif


      elseif (efswtch.eq."method4") then  !-YuP: method4 only.
c        If initial time step, use spitzer+neoclassical to get
c          electric field near what is required for specifed 
c          target currxj (convert to volts/cm from cgs)

         if (n.eq.0) then

            currxj0(lr_)=currxj(lr_) 
            call restcon
            call resthks(l_,lr_,lmdpln_,
     +           zreshin(lr_),zreskim(lr_),zressau1,zressau2)
            elecfld(lr_)=currxj(lr_)*(zreshin(lr_)*sptzr(l_))*300.*3.e9
            
         else ! n>0

c           Relaxation towards target parallel currxj:

c           Presently, the calculated parallel current
c           in diaggnde is averaged over the area cross-section
c           (at cnst toroidal angle) between a flux surface.
c           This quantity is depracated, and should be replaced
c           in the future.  Here, for method5, we adjust currpar 
c           back to parallel current at the minimum B position.
c           See comments in eflditer.f, on psifct.

            psifct1=1.
            if(efswtch.eq."method5") then
               psifct1=psiovr(lr_)/onovrp(1,lr_)
            endif

            write(*,*)'efield:elecfld,currpar/psifct,currxj',
     +            elecfld(lr_),currpar(lr_)/psifct1,currxj(lr_)

            currxj0(lr_)=currxj(lr_)
            elecfld(lr_)=elecfld(lr_)*(1.-efrelax*
     +           (currpar(lr_)/psifct1 -currxj(lr_))/
     +           (sign(one,currxj(lr_))
     +           *max(abs(currpar(lr_)/psifct1),abs(currxj(lr_)))))

            write(*,*)'efield:elecfld',elecfld(lr_)

         endif ! n


      elseif (efswtch.eq."method5"  .and.  n.eq.0) then  
          ! calculate 'seed' electric field before starting iterations
      
      
c-YuP: For method5, n>0, elecfld and currxj are calculated in eflditer.f;
c-YuP: currxj0 is set in tdrmshst (at n=0).
c          If initial time step, use spitzer+neoclassical to get
c          electric field near what is required for specifed 
c          target currxj (convert to volts/cm from cgs)
            call restcon
            call resthks(l_,lr_,lmdpln_,
     +           zreshin(lr_),zreskim(lr_),zressau1,zressau2)
            elecfld(lr_)=currxj(lr_)*(zreshin(lr_)*sptzr(l_))*300.*3.e9
         
         
      endif ! selection of efswtch=
      

      return
      end
