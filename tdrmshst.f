
c
c
      subroutine tdrmshst
      implicit integer (i-n), real*8 (a-h,o-z)
c..................................................................
c     This routine initializes the radial meshes and integration
c     coefficients, based on rya() values.
c     For eqmod="enabled":
c       Bin boundaries are DEFINED by the average outer equatorial
c       plane major radius between those of neighboring Fokker-Planck
c       points.
c     Also, gives  currxj for efswtch="method5" .
c     Sets some radial quantities related to neoclassical transp.
c            
c..................................................................
      save
      include 'param.h'
      integer itlrza
      parameter(itlrza=3*lrza+1)
      include 'comm.h'
CMPIINSERT_INCLUDE

      dimension d2equ(lrza),d2areamd(lrza),d2volmid(lrza),worka(lrza)
      dimension d2volcon(lrza)
      dimension workk(itlrza), wk(lrzmax)

      character*8 ifirst
      data ifirst /"first"/


      if(eqmod.eq."enabled") then
        rmn=rhomax
      else
        rmn=radmin
      endif
cBH050917      if (eqmod.eq."disabled") rmn=radmin

      if (ifirst.eq."first") then 
      !!YuP[07-2017] added: do this only at 1st call
      do 10 ll=1,lrzmax
        rz(ll)=rya(ll)*rmn
        rrz(ll)=rrz(ll)*rmn !YuP[07-2017]: recursive: make sure it is done once!
        ! if this subroutine is called more than once, the above line 
        ! is only performed once.
 10   continue
      rrz(lrzmax)=rmn
      ifirst="notfirst"
      endif !YuP[07-2017] added: do this only at 1st call
      
      if (eqmod.eq."disabled") then
        tpisr0=2.*pi**2*radmaj
        dvol(1)=tpisr0*rrz(1)**2
        darea(1)=dvol(1)/(2.*pi*radmaj)
        do 7 ll=1,lrzmax
          tr(ll)=tpisr0*rrz(ll)**2
 7      continue
        do 700 ll=2,lrzmax
          dvol(ll)=tr(ll)-tr(ll-1)
          darea(ll)=dvol(ll)/(2.*radmaj*pi) !here: for eqmod.eq."disabled"
 700    continue

      elseif (eqmod.eq."enabled") then

        do 25 l=1,lrzmax
          equilpsp(l)=-equilpsi(l)
 25     continue
        equilpsp(0)=-psimag
        areamid(0)=0.
        volmid(0)=0.
        rpmconz(0)=rmag
c       Bin boundaries are DEFINED by the average outer equatorial
c       plane major radius between those of neighboring Fokker-Planck
c       points.
        do 16 l=1,lrzmax
          if (l.lt.lrzmax) then
            rpmconz(l)=(rpconz(l)+rpconz(l+1))*.5
          else
            rpmconz(lrzmax)=rmaxcon
          endif
c     rmincon and rmaxcon are the inner and outer major radii of the 
c     LCFS at ez=0.
          !YuP[03-2016] in the next line, it was terp2(rpmconz(l),zero,...)
          !why zero? changed to zmag
          psivalm(l)=terp2(rpmconz(l),zmag,nnr,er,nnz,ez,epsi,epsirr,
     1      epsizz,epsirz,nnra,0,0)
          epsicon_=psivalm(l)
          call eqorbit(epsicon_)
          call eqvolpsi(epsicon_,volmid(l),areamid(l))
          !!write(*,'(i4,3f15.3)') l,rpmconz(l),rpconz(l),areamid(l)
          ! YuP [Apr/2014]: jump in areamid(l) at last l=lrzmax, 
          ! because of jump in rpmconz(l).
          ! Not really a bug; it depends how you define 
          ! the bin around last FP'ed surface.
 16     continue
        equilpsi(0)=psimag
        do 60 l=1,lrzmax
          worka(l)=-psivalm(l)
 60     continue
        i1p(1)=4
        i1p(2)=4
        itab(1)=1
        itab(2)=0
        itab(3)=0
C     cannot use 4-point interpolation if lrzmax=3
c     (should not run with lrzmax=2, as cubic spline not defined)
        if (lrzmax .eq. 3) then
          i1p(1)=2
          i1p(2)=2
          d2volmid(1)=0.0
          d2volmid(lrzmax)=(volmid(lrzmax) - volmid(lrzmax-1)) /
     /      (worka(lrzmax)   - worka(lrzmax-1))
          d2areamd(1)=0.0
          d2areamd(lrzmax)=(areamid(lrzmax) - areamid(lrzmax-1)) /
     /      (worka(lrzmax)   - worka(lrzmax-1))
        endif
c
        call coeff1(lrzmax,worka(1),volmid(1),d2volmid,i1p,1,workk)
        call terp1(lrzmax,worka(1),volmid(1),d2volmid,-psilim,1,tab,
     +    itab)
        volmid(lrzmax)=tab(1)
        call coeff1(lrzmax,worka(1),areamid(1),d2areamd,i1p,1,workk)
        call terp1(lrzmax,worka(1),areamid(1),d2areamd,-psilim,1,tab,
     +    itab)
        areamid(lrzmax)=tab(1)
        psivalm(lrzmax)=psilim
        do 15 ll=1,lrzmax
          darea(ll)=(areamid(ll)-areamid(ll-1)) !here: for eqmod.eq."enabled"
          dvol(ll)=(volmid(ll)-volmid(ll-1))
 15     continue
       ! write(*,*)'tdrmshst: sum(darea)=',sum(darea)

c..................................................................
c     Set up spline arrays for R(ll) and dpsidr(ll) (both at mag axis)
c..................................................................

c     d2bmdpl is set for use with freya.  The values of bmdplne
c     (from the center of the l-th volume)
c     are here attributed to values of psi at the outside of the
c     l-th volume.  (not sure about accuracy here. BobH, 950221).

        call coeff1(lrzmax,worka(1),bmdplne,d2bmdpl,i1p,1,workk)

        !YuP160304 call coeff1(lrzmax,rz(1),equilpsi(1),d2equ,i1p,1,workk)
        tr(0:lrzmax)= psimag-equilpsi(0:lrzmax)
        !pol.flux, in ascending order needed for coeff1
        call coeff1(lrzmax,rya(1),tr(1),d2equ,i1p,1,workk)
        itab(1)=0
        itab(2)=1
        itab(3)=0
        do 50 l=1,lrzmax
          !YuP160304 call terp1(lrzmax,rz(1),equilpsi(1),d2equ,rz(l),1,tab,itab)
          call terp1(lrzmax,rya(1),tr(1),d2equ,rya(l),1,tab,itab)
          dpsidrho(l)=-tab(2) ! '-' because we used ascending psi function
 50     continue

c.......................................................................
c     Determine the del(psi) array
c     dpsi(l) must equal psivalm(l)-psivalm(l-1) or something is wrong.
c.......................................................................

        do 20 l=1,lrzmax ! YuP[01-2017] was 1,lrz 
          ilr=lrindx(l)
          dpsi(ilr)=dvol(ilr)*bmdplne(ilr)/(twopi*zmaxpsi(ilr))
          if (eqsym.ne."none") dpsi(ilr)=dpsi(ilr)*0.5
 20     continue

      endif  ! On eqmod

CMPIINSERT_IF_RANK_EQ_0
      WRITE(*,*)'tdrmshst: lr, dvol(lr), darea(lr) based on eqvolpsi'
      do ll=1,lrzmax ! YuP[01-2017] was 1,lrz
          WRITE(*,'(i6,2e12.4)') ll,  dvol(ll),  darea(ll)
      enddo
      WRITE(*,'(a,2e12.4)')
     + 'tdrmshst: sum(dvol),sum(darea) based on eqvolpsi',
     +            sum(dvol),sum(darea)
      WRITE(*,*)'----------------------------------------------------'
CMPIINSERT_ENDIF_RANK
c.......................................................................
c     Compute H*rho=dV/drho/(4*pi**2*R_0) (used for transport)
c.......................................................................

      i1p(1)=4
      i1p(2)=4
      itab(1)=0
      itab(2)=1
      itab(3)=0
      call coeff1(lrzmax,rz(1),volcon(1),d2volcon,i1p,1,workk)
      do 70 l=0,lrzmax-1
        drp5(l)=rz(l+1)-rz(l)
        vx=(rz(l)+rz(l+1))*.5
        call terp1(lrzmax,rz(1),volcon(1),d2volcon(1),vx,1,tab,itab)
        h_r(l)=tab(2)/(4.*pi**2*radmaj)
 70   continue
      call terp1(lrzmax,rz(1),volcon(1),d2volcon(1),radmin,1,tab,itab)
C%OS  h_r(lrzmax)=tab(2)
      h_r(lrzmax)=tab(2)/(4.*pi**2*radmaj)

c.......................................................................
c     If efswtch="method5",
c     determine target parallel current (A/cm**2) from eqdsk, 
c     MKS formula is: j_parallel=R*(B_phi/|B|)dp/dpsi+
c                                (|B|/mu_0)*dF/dpsi
c     Also, add wedge of edge current outside r/a=0.8.
c.......................................................................

      do l=1,lrzmax
         curreq(l)=10.*rpcon(l)*(btor0(l)/bmidplne(l))*pppsi(l)
     +             +bmidplne(l)/(4.*pi*1.e-1)*fppsi(l)
      enddo


      if (efswtch.eq."method5") then
      curr_edge_roa=0.8   
      do 30 l=1,lrzmax
         currxj(l)=curreq(l)
         if (rovera(l).gt.0.8) then
            currxj(l)=currxj(l)+
     +           (rovera(l)-curr_edge_roa)/(1.-curr_edge_roa)*curr_edge
         endif
         currxj0(l)=curreq(l) !-YuP: added 
         ! the target parallel current (moved here from efield.f)
 30   continue
      endif

c.......................................................................

c     Set some quantities related to neoclassical radial transport
c     of first ion species:
c     tauii is inverse of nu_perp in NRL tables, evaluated at
c           thermal energy, E=1.5 Ti.
c     drr_gs is thermal Galeev and Sagdeev radial diffusion coeff,
c           per Miyamoto, Plasma Physics for Nuclear Fusion, Eq 8.24.
c     tau_neo=(r-radmin)**2/drr_gs, ~time to diffuse to edge 
c              at local rate.
c     rhol and rhol_pol are the ion Larmor radius and pol Larmor radius.
c     The analogous  quantities are calculated for the maximum neutral
c         beam energy, taubi,drr_gs_b, tau_neo_b, rhol_b, rhol_pol_b,
c         assuming mass of first ions species, 80 keV NB energy.
c     
c.......................................................................

      if (gamaset.ne.zero) then
         gamaset1=gamaset
      else
         gamaset1=17.
      endif

      do l=1,lrzmax

         fmu=fmass(kionn)/proton
         beamengy=80.   !Assumed beam energy (keV)

         tauii(l)=1./(1.4e-7*bnumb(kionn)**2*reden(kelec,l)*zeff(l)*
     1        gamaset1/(fmu**0.5*sqrt(temp(kionn,l)*1.e3)*
     1        1.5*temp(kionn,l)*1.e3))
         rhol(l)= vth(kionn,l)/
     1        (bnumb(kionn)*charge*bmod0(l)/(fmass(kionn)*clight))
         rhol_pol(l)=eps(l)**0.5*vth(kionn,l)/
     1        (bnumb(kionn)*charge*bthr0(l)/(fmass(kionn)*clight))
         drr_gs(l)=eps(l)**0.5*rhol_pol(l)**2/(tauii(l)*eps(l))
         tau_neo(l)=(rz(l)-rmn)**2/drr_gs(l)
         
cBH120519:  Small effect correction, checking with NRL Plasma Formulary
c         taubi(l)=1./(1.4e-7*bnumb(kionn)**2*reden(kelec,l)*zeff(l)*
c         taubi is perp collision time of fast ions on ions
         taubi(l)=1./(1.8e-7*bnumb(kionn)**2*reden(kelec,l)*zeff(l)*
     1        gamaset1/(fmu**0.5*sqrt(beamengy*1.e3)*
     1        beamengy*1.e3))
         rhol_b(l)= sqrt(2.*beamengy*ergtkev/fmass(kionn))/
     1        (bnumb(kionn)*charge*bmod0(l)/(fmass(kionn)*clight))
         rhol_pol_b(l)=eps(l)**0.5*
     1        sqrt(2.*beamengy*ergtkev/fmass(kionn))/
     1        (bnumb(kionn)*charge*bthr0(l)/(fmass(kionn)*clight))
         drr_gs_b(l)=eps(l)**0.5*rhol_pol_b(l)**2/(taubi(l)*eps(l))
         tau_neo_b(l)=(rz(l)-rmn)**2/drr_gs_b(l)

      enddo

c
c.......................................................................
c     For use in iterative soln of Ampere-Faraday eqns.
c     compute drpmconz (distance between bin boundaries), and
c             dlpgpsii (integrated pol dist * (grad psi(l)/grad psi(1)))
c             dlpsii (poloidal distance integrated psi=B/B0)
c.......................................................................

      do ll=1,lrzmax

c$$$         if (ll.eq.1) then
c$$$            do l=1,lorbit(ll)
c$$$               delrho(l,ll)=abs((psivalm(1)-psimag)/
c$$$     +              (solr(l,ll)*eqbpol(l,ll)))
c$$$            enddo
c$$$         else
c$$$            do l=1,lorbit(ll)
c$$$               delrho(l,ll)=abs((psivalm(ll)-psivalm(ll-1))/
c$$$     +              (solr(l,ll)*eqbpol(l,ll)))
c$$$            enddo
c$$$         endif
         
         ! Define the bin width around bin center #ll 
         ![which is at R=rpconz(ll), having upper boundary at R=rpmconz(ll)]
         drpmconz(ll)=rpmconz(ll)-rpmconz(ll-1)
         !Note: rpmconz(0)=rmag, and  rpmconz(lrzmax)=rmaxcon
         !so that   drpmconz(1)= rpmconz(1)-Rmag
         ! and      drpmconz(lrzmax)= rmaxcon-rpmconz(lrzmax-1)

         dlpgpsii(ll)=zero
         do l=2,lorbit(ll)
cBH170304            dlpgpsii(ll)=dlpgpsii(ll)+eqdell(l,ll)*0.25*
cBH170304     +      (solr(l-1,ll)+solr(l,ll))*(eqbpol(l-1,ll)+eqbpol(l,ll))
            dlpgpsii(ll)=dlpgpsii(ll)+eqdell(l,ll)*
     +           0.5*(eqbpol(l-1,ll)+eqbpol(l,ll))
         enddo
cBH170304         dlpgpsii(ll)=dlpgpsii(ll)/(solr(1,ll)*eqbpol(1,ll))
         dlpgpsii(ll)=dlpgpsii(ll)/eqbpol(1,ll)

         dlpsii(ll)=zero 
         ! This is the integral over dl_pol*Btor(pol.dist)/B0 
         ! Constant at each ll surface: 
         btor_r = abs(fpsi(ll)) !This is Btor*R  ! YuP: abs() is ok?
         do l=2,lorbit(ll)
            lm=l-1
cYuP            dlpsii(ll)=dlpsii(ll)+(es(l,ll)-es(l-1,ll))*
cYuP     +           0.5*(bpsi(l-1,ll)+bpsi(l,ll))  !YuP[03-2017] Changed to:
            !YuP: setup the tor field:
            btor_l=  btor_r/solr(l,ll)  !== fpsiar/R at l
            btor_lm= btor_r/solr(lm,ll) !== fpsiar/R at l-1
            dlpsii(ll)= dlpsii(ll) +eqdell(l,ll)*0.5*(btor_lm+btor_l)
            !it uses dl_pol == eqdell(l,ll)
         enddo
         dlpsii(ll)=dlpsii(ll)/bmidplne(ll) ! Divided by B0 (equatorial)
c         write(*,*)'tdrmshst: ll, bmidplne(ll),fpsi(ll)/solr(1,ll)=',
c     +    ll, bmidplne(ll),fpsi(ll)/solr(1,ll)

      enddo ! ll
      
      ! YuP[2018-01-12] Added adjustment of dlpsii(1) and dlpgpsii(1)
      ! which are used in subr.ampfsoln().
      !
      ! Note that drpmconz(ll) is the bin width around surface#ll
      ! at the outboard midplane,
      ! drpmconz(ll)=rpmconz(ll)-rpmconz(ll-1)
      ! Or we can write: drpmconz(ll)= (Rp(ll+1)-Rp(ll-1))/2
      ! for ll= 2 : lrz-1
      ! where Rp(ll) is the R coord at radial bin center #ll.
      ! For ll=1, it is different:
      ! drpmconz(1)= rpmconz(1)-rpmconz(0) = 0.5(Rp(1)+Rp(2)) - Rmag
      ! i.e., it is the distance between magn.axis 
      ! and the midpoint between 1st and 2nd surfaces.
      ! It may happen that the position of the 1st bin center,
      ! Rp(1)==rpconz(1),
      ! is NOT in the CENTER of the 1st bin.
      ! For the proper evaluation of the integral f(1)*dr(1)
      ! where dr(1) is the bin width, dr(1)==0.5(Rp(1)+Rp(2)) - Rmag,
      ! we need to evaluate function f(1) at the CENTER of the 1st bin.
      ! This is not important when f(r) is flat near r~0 
      ! (r is minor radius)
      ! but it is important when f(r) is proportional to r, for example
      ! (such as the dlpsii(r) function).
      !
      ! Adjust dlpsii(1) to the CENTER of the first radial bin.
      ! The 1st bin has the left boundary at R=Rmag  (or r=0),
      ! and the right boundary at R= Rmag + drpmconz(1) = 0.5(Rp(1)+Rp(2))
      dlpsii(1)= ( (dlpsii(1)+dlpsii(2))*0.5 +0.d0 )*0.5 
      ! At small flux surface (~circle), dlpsii is approximately 2pi*r,
      ! where r is the minor radius.
      ! In the above, it is assumed that this function is exacly 0 
      ! at the magnetic axis, and it goes linearly at small r.
      
      ! Another adjustment is for dlpgpsii(1) -
      ! Interpolate values of dlpgpsii to the OUTER (upper) border
      ! of radial bin; For a bin#ll with center at rpconz(ll),
      ! the upper border is at  rpmconz(ll) .
      ! Reason: the cumulative integrals over dlpsii(r),
      ! effectively the integrals over r*dr, 
      ! are done numerically as SUM(r(ll)*dr(ll)), ll=1:lll
      ! This method gives a larger value than exact value
      ! because of the last point -- r(lll)*dr(lll).
      ! It contains an extra 0.5*r(lll)*dr(lll) value,
      ! So the values of such integrals correspond to r(lll)+0.5*dr(lll)
      ! rather to r(lll).
      ! Further in calculations, 
      ! the values of SUM(r(ll)*dr(ll)) are divided by 
      ! dlpgpsii(l), effectively we have terms from
      ! SUM(r(ll)*dr(ll))/r, which should scale as (r^2/2)/r = r/2.
      ! But because the values of SUM(r(ll)*dr(ll)) correspond to
      ! r(lll)+0.5*dr(lll), we should also adjust the denominator
      ! to the same point (r --> r+0.5*dr), which is the upper boundary 
      ! of the radial bin lll.
      do ll=1,lrzmax-1
         rrr= (rpmconz(ll)-rpconz(ll))/(rpconz(ll+1)-rpconz(ll))
         wk(ll)= dlpgpsii(ll) +(dlpgpsii(ll+1)-dlpgpsii(ll))*rrr
      enddo
      do ll=1,lrzmax-1
         dlpgpsii(ll)=wk(ll)
      enddo
      ! Note that the last point (ll=lrzmax) is omitted.
      ! Instead, we adjust the width of the last bin:
      ! Original: drpmconz(lrzmax)= rmaxcon-rpmconz(lrzmax-1)
      !drpmconz(lrzmax)= rpconz(lrzmax)-rpmconz(lrzmax-1) ! Adjusted,
      ! which sets it to half-width of the interior-cell width.
      ! (This adjustment is not very useful: gives a small bend in E(r)
      ! at the edge of r-grid).
      !
      ! The adjustments above give a small improvement in behavior 
      ! of E near R=Rmag  (r~0).
      
      return
      end
