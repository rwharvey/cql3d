c
c
      subroutine cfpcoefn
      implicit integer (i-n), real*8 (a-h,o-z)
      save
ccc      real*8,dimension(iy):: prnt1,prnt2,prnt3
ccc      real*8,dimension(iy):: prnt4,prnt5,prnt6
c..................................................................
c     Subroutine to calculate bounce-averaged Fokker-Planck collision
c     coefficients (relativistic corrections added by MARK FRANZ--
c     U.S.A.F.)
c     If (cqlpmod .eq. "enabled") then compute only the coefficients
c     at the orbit position l=l_ and do not perform the bounce-averages.
c..................................................................
      include 'param.h'
      include 'comm.h'
c
      data naccel/100/
      ialpha=2
      impcoef=0
      if (n.eq.1) then
        continue
      elseif (mod(n,ncoef).eq.1 .or. ncoef.eq.1) then
        continue
      else
        return
      endif
      impcoef=1
      nccoef=nccoef+1
      call bcast(cal(1,1,1,l_),zero,iyjx*ngen)
      call bcast(cbl(1,1,1,l_),zero,iyjx*ngen)
      call bcast(ccl(1,1,1,l_),zero,iyjx*ngen)
      call bcast(cdl(1,1,1,l_),zero,iyjx*ngen)
      call bcast(cel(1,1,1,l_),zero,iyjx*ngen)
      call bcast(cfl(1,1,1,l_),zero,iyjx*ngen)
      call bcast(eal(1,1,1,1,l_),zero,iyjx*ngen*2)
      call bcast(ebl(1,1,1,1,l_),zero,iyjx*ngen*2)

c..................................................................
c     if only gen. species contributions are desired execute jump..
c..................................................................
      if (colmodl.eq.2 ) goto 110
c..................................................................
c     if no background species exist execute jump
c..................................................................
      if (ntotal.eq.ngen) goto 110
      iswwflag=0

      do 100 k=ngen+1,ntotal

c..................................................................
c     If this species is not to be included as a field (background)
c     species for calculation of the collision integral, jump out.
c..................................................................
        if (kfield(k).eq."disabled") go to 100

c..................................................................
c     Determine the Maxwellian distribution associated with
c     background species k. This will be a relativistic Maxwellian
c     for relativistic calculations.
c..................................................................

        rstmss=fmass(k)*clite2/ergtkev
        reltmp=rstmss/temp(k,lr_)
        if (cqlpmod .eq. "enabled") reltmp=rstmss/temppar(k,ls_)
        if (reltmp .gt. 100. .or. relativ .eq. "disabled") then
          ebk2=sqrt(pi/(2.*reltmp))
        else if (reltmp .lt. .01) then
          ebk2=2.*exp(reltmp)/reltmp**2
        else
          call cfpmodbe(reltmp,ebk1,ebk2)
        endif
        rnorm=reltmp/(4.*pi*cnorm3*ebk2)
        call bcast(temp1(0,0),zero,iyjx2)
c..................................................................
c     Need extra mesh points to represent ions on a mesh meant to
c     support electrons. Need more resolution near zero.
c     Split each velocity bin into nintg pieces and integrate.
c..................................................................
c990131        nintg0=41.*amax1(1.,sqrt(fmass(k)/1.67e-24))
        nintg0=41.*max(1.d0,sqrt(fmass(k)/1.67e-24))
c990131        nintg=max0(nintg0,min0(51,int(x(2)/sqrt(.2*cnorm2/reltmp))))
        nintg=max(nintg0,min0(51,int(x(2)/sqrt(.2*cnorm2/reltmp))))
        nintg=2*(nintg/2)+1
        do 10 i=1,nintg
          sfac=4.+2.*(2*(i/2)-i)
          if (i.eq.1 .or. i.eq.nintg) sfac=1.
          sfac=sfac/3.
          do 11 j=2,jx
            xx=x(j)-(i-1)*dxm5(j)/dfloat(nintg-1)
            xx2=xx**2
            gg=sqrt(1.+xx2*cnorm2i)
            gg2=gg**2
            if(cnorm2i*xx2-1.e-5 .gt. 0.d0) then
               expon=gg-1.
            else
               expon=.5*xx2/cnorm2
            endif
            fff=sfac*rnorm*exp(-expon*reltmp)*dxm5(j)/dfloat(nintg-1)
            temp1(1,j)=temp1(1,j)+xx*gg*fff
            temp1(2,j)=temp1(2,j)+xx*fff
            temp1(3,j)=temp1(3,j)+xx2*fff
            temp1(4,j)=temp1(4,j)+xx2*fff/gg
            temp1(5,j)=temp1(5,j)+(xx2/gg)**2*fff
            temp1(6,j)=temp1(6,j)+(xx2/gg2)**2*gg*fff
 11       continue
 10     continue
c..................................................................
        tam2(jx)=0.
        tam3(1)=0.
        tam5(1)=0.
        tam6(jx)=0.
        tam7(1)=0.
        tam9(1)=0.
c..................................................................
c     tam2(jx) and tam6(jx) represent integrals from xmax to infinity
c     The following coding performs this integration. In this case 21*xmax
c     represents infinity. The lack of this piece is most obvious when
c     ions are a general species and electrons are fixed Maxwellians.
c..................................................................
        do 15 ll2=1,21
          sfac=4.+2.*(2*(ll2/2)-ll2)
          if ( ll2.eq.1 .or. ll2.eq.21) sfac=1.
          sfac=sfac*x(jx)/60.
          do 16 ll1=1,20
            xx=(realiota(ll1)+.05*realiota(ll2-1))*x(jx)
            xx2=xx**2
            gg=sqrt(1.+xx2*cnorm2i)
            if(cnorm2i*xx2-1.e-5 .gt. 0.d0) then
               expon=gg-1.
            else
               expon=.5*xx2/cnorm2
            endif
            fff=sfac*rnorm*exp(-expon*reltmp)
            tam2(jx)=tam2(jx)+xx*gg*fff
            tam6(jx)=tam6(jx)+xx*fff
 16       continue
 15     continue
        do 20 j=2,jx
          jj=jx+1-j
          jp=jj+1
          jm=j-1
c..................................................................
c     tam2 - M0; tam3 - N0; tam5 - E0
c     tam6 - M0';  tam7 - N0';   tam9 - E0'
c     see UCRL manual..
c..................................................................
          tam2(jj)=tam2(jp)+temp1(1,jp)
          tam3(j)=tam3(jm)+temp1(3,j)
          tam5(j)=tam5(jm)+temp1(5,j)
          tam6(jj)=tam6(jp)+temp1(2,jp)
          tam7(j)=tam7(jm)+temp1(4,j)
          tam9(j)=tam9(jm)+temp1(6,j)
 20     continue
c..................................................................
        do 30 j=2,jx
          tam10(j)=cog(0,1)*(3.*tam7(j)+cnorm2i*(2.*xm(j,3)*tam6(j)-
     *      tam9(j)))*gamsqr(j)
          tam11(j)=cog(0,1)*(xsq(j)*tam2(j)+
     *      gamsqr(j)*xm(j,-1)*tam5(j))
          tam12(j)=cog(0,1)*(tam2(j)+1.5*xm(j,-1)*tam3(j)-
     *      .5*xm(j,-3)*tam5(j))
 30     continue

c.......................................................................
c     Perform the bounce-average for the background species and introduce
c     the contribution to all general species coeff.
c.......................................................................

        if (cqlpmod .ne. "enabled") then
          call bavdens(k)
        else
          do 59 i=1,iy
            bavdn(i,lr_)=denpar(k,ls_)
            bavpd(i,lr_)=denpar(k,ls_)*sinn(i,l_)
 59       continue
        endif

c     Below, eal and ebl are to save contributions to the collisional
c     coefficients cal() and cbl(i,j,k,l_) for general species k and
c     radial location l_, resulting from electrons (eal and  ebl
c     of (i,j,k,1,l_)), and resulting from the sum of the
c     effects of the background ions, (eal and ebl of (i,j,k,2,l_)). 
c     eal and ebl are used later for calculating powers from the
c     general  to the Max. species.

        do 80 kk=1,ngen
          anr1=gama(kk,k)*satioz2(k,kk)*one_
          anr2=anr1*satiom(kk,k)
          call bcast(tem1,zero,iyjx)
          call bcast(tem2,zero,iyjx)
          
          do 70 j=2,jx
            ttta=anr2*tam10(j)*gamefac(j)
            tttb=anr1*tam11(j)*gamefac(j)
            tttf=anr1*tam12(j)*gamefac(j)
            do 60 i=1,iy
              jj=i+(j-1)*iy 
              tem1(jj)=ttta*vptb(i,lr_)
     1          *bavdn(i,lr_)
              tem2(jj)=tttb*vptb(i,lr_)
     1          *bavdn(i,lr_)
              cal(i,j,kk,l_)=cal(i,j,kk,l_)+tem1(jj)
              cbl(i,j,kk,l_)=cbl(i,j,kk,l_)+tem2(jj)
ccc              if (j.eq.2) prnt1(i)=tttf
ccc              if (j.eq.3) prnt2(i)=tttf
ccc              if (j.eq.4) prnt3(i)=tttf
              cfl(i,j,kk,l_)=cfl(i,j,kk,l_)+tttf*vptb(i,lr_)
     1          *bavpd(i,lr_)
 60         continue
 70       continue

ccc          do i=1,iy
ccc             prnt4(i)=vptb(i,lr_)
ccc             prnt5(i)=bavpd(i,lr_)
ccc          enddo

          if(k.eq.kelecm) then
            call daxpy(iyjx,one,tem1,1,eal(1,1,kk,1,l_),1)
            call daxpy(iyjx,one,tem2,1,ebl(1,1,kk,1,l_),1)
          else
            do 101 i=1,nionm
              if(k.eq.kionm(i)) then
                call daxpy(iyjx,one,tem1,1,eal(1,1,kk,2,l_),1)
                call daxpy(iyjx,one,tem2,1,ebl(1,1,kk,2,l_),1)
              endif
 101        continue
          endif

 80     continue ! kk=1,ngen
 100  continue ! k=ngen+1,ntotal

 110  continue
ccc      kk=1
ccc      do i=1,iy
ccc         prnt1(i)=cal(i,2,kk,l_)
ccc         prnt2(i)=cbl(i,2,kk,l_)
ccc         prnt3(i)=ccl(i,2,kk,l_)
ccc         prnt4(i)=cdl(i,2,kk,l_)
ccc         prnt5(i)=cel(i,2,kk,l_)
ccc         prnt6(i)=cfl(i,2,kk,l_)
ccc      enddo

c..................................................................
c     Calculate coefficients and their bounce-averages for a general species
c..................................................................

      if (colmodl.eq.1) goto 700
      if (colmodl.eq.4 .and. n.ge.naccel) then
        do 250 k=1,ngen
          if (k.ne.ngen) then
            do 220 j=1,jx
              do 210 i=1,iy
                if(2.*f(i,j,k,l_)-f_(i,j,k,l_).gt.0.) then
                  fxsp(i,j,k,l_)= 2.*f(i,j,k,l_)-f_(i,j,k,l_)
                else
                  fxsp(i,j,k,l_)= .5*f(i,j,k,l_)
                endif
 210          continue
 220        continue
          else
            call diagescl(kelecg)
            call dcopy(iyjx2,f(0,0,ngen,l_),1,fxsp(0,0,ngen,l_),1)
          endif
 250    continue
      else
        do 280 k=1,ngen
          if (colmodl.eq.4 .and. k.eq.ngen) call diagescl(kelecg)
 280    continue
      endif
      madd=1
      if (machine.eq."mirror") madd=2

c.......................................................................
c     loop over orbit mesh s(l) (for general species only)
c     CQLP case: compute only s(l_)
c.......................................................................

      if (cqlpmod .ne. "enabled") then
        iorbstr=1
        iorbend=lz  ! Whole flux surface, eqsym="none", else half FS.
      else
        iorbstr=l_
        iorbend=l_
      endif
      do 600 l=iorbstr,iorbend
        ileff=l
        if (cqlpmod .eq. "enabled") ileff=ls_

        do 500 k=1,ngen
c.................................................................
c     Jump out if this species is not to be used as a background species
c..................................................................
          if (kfield(k).eq."disabled") go to 500

c     zero ca, cb,.., cf  : 
CPTR>>>REPLACE PTR-BCASTCACD
          call bcast(ca,zero,iyjx)
          call bcast(cb,zero,iyjx)
          call bcast(cc,zero,iyjx)
          call bcast(cd,zero,iyjx)
          call bcast(ce,zero,iyjx)
          call bcast(cf,zero,iyjx)
CPTR<<<END PTR-BCASTCACD
          mu1=0
          mu2=mx
          mu3=madd
          if (colmodl.eq.3) then
            mu1=1
            mu2=mx
            mu3=1
          endif
          do 400 m=mu1,mu2,mu3
            if (colmodl.eq.4 .and. n.ge.naccel) then
              call dcopy(iyjx2,fxsp(0,0,k,l_),1,temp3(0,0),1)
            else
              call dcopy(iyjx2,f(0,0,k,l_),1,temp3(0,0),1)
            endif

c     compute V_m_b in tam1(j), for given l, m and gen. species b
c     coeff. of Legendre decomposition of f (using temp3)
            call cfpleg(m,ileff,1) !-> tam1

            tam2(jx)=0.
            tam3(1)=0.
            tam4(jx)=0.
            tam5(1)=0.
            tam6(jx)=0.
            tam7(1)=0.
            tam8(jx)=0.
            tam9(1)=0.

c     prepare arrays for computation of M_m, N_m, R_m and E_m
            do 301 j=2,jx
              jm=j-1
              tam20(j)=.5*dxm5(j)*
     *          (gamsqr(j)*xm(j,ix1(m))*tam1(j)+
     *          gamsqr(jm)*xm(jm,ix1(m))*tam1(jm))
              tam13(j)=.5*dxm5(j)*
     *          (gamsqr(j)*xm(j,ix2(m))*tam1(j)+
     *          gamsqr(jm)*xm(jm,ix2(m))*tam1(jm))
              tam14(j)=.5*dxm5(j)*
     *          (gamsqr(j)*xm(j,ix3(m))*tam1(j)+
     *          gamsqr(jm)*xm(jm,ix3(m))*tam1(jm))
              tam15(j)=.5*dxm5(j)*
     *          (gamsqr(j)*xm(j,ix4(m))*tam1(j)+
     *          gamsqr(jm)*xm(jm,ix4(m))*tam1(jm))
              tam16(j)=.5*dxm5(j)*
     *          (gamma(j)*xm(j,ix1(m))*tam1(j)+
     *          gamma(jm)*xm(jm,ix1(m))*tam1(jm))
              tam17(j)=.5*dxm5(j)*
     *          (gamma(j)*xm(j,ix2(m))*tam1(j)+
     *          gamma(jm)*xm(jm,ix2(m))*tam1(jm))
              tam18(j)=.5*dxm5(j)*
     *          (gamma(j)*xm(j,ix3(m))*tam1(j)+
     *          gamma(jm)*xm(jm,ix3(m))*tam1(jm))
              tam19(j)=.5*dxm5(j)*
     *          (gamma(j)*xm(j,ix4(m))*tam1(j)+
     *          gamma(jm)*xm(jm,ix4(m))*tam1(jm))
 301        continue

c     note: ialpha set to 2 at top of subroutine
            if (m.ge.1 .or. ialpha.eq.2) goto 308
c..................................................................
c     The do loops 302 through 307 seek to take advantage of the
c     Maxwellian nature of f for v < vth/2. This is employed in
c     the integrations to obtain the functionals - Thus f is
c     assumed to be Maxwellian between velocity mesh points,
c     not linear. (not used)
c..................................................................
            do 302 j=2,jx
              xs=sqrt(temp(k,lr_)*ergtkev*0.5/fmass(k))
              if (cqlpmod .eq. "enabled")
     +          xs=sqrt(temppar(k,ls_)*ergtkev*0.5/fmass(k))
              if (x(j)*vnorm.gt.xs) goto 303
 302        continue
 303        continue
            jthov2=j-1
c..................................................................
c     Determine the "local" Maxwellian between meshpoints.
c..................................................................
            do 304 j=1,jthov2
c990131              tam21(j)=alog(tam1(j)/tam1(j+1))/(gamm1(j+1)-gamm1(j))
              tam21(j)=log(tam1(j)/tam1(j+1))/(gamm1(j+1)-gamm1(j))
              tam22(j)=tam1(j)/exp(-(gamm1(j)*tam21(j)))
 304        continue
            rstmss=fmass(k)*clite2/ergtkev
            reltmp=rstmss/temp(k,lr_)
            if (cqlpmod .eq. "enabled") reltmp=rstmss/temppar(k,ls_)
            call bcast(tam23,zero,jx*8)
            nintg=max0(21,min0(51,int(x(2)/sqrt(.2*cnorm2/reltmp))))
            nintg=2*(nintg/2)+1
            do 305 i=1,nintg
              sfac=4.+2.*(2*(i/2)-i)
              if (i.eq.1 .or. i.eq.nintg) sfac=1.
              sfac=sfac/3.
              do 306 j=2,jthov2
                xx=x(j)-(i-1)*dxm5(j)/dfloat(nintg-1)
                xx2=xx**2
                gg=sqrt(1.+xx2*cnorm2i)
                gg2=gg**2
                if(cnorm2i*xx2-1.e-5 .gt. 0.d0) then
                   expon=gg-1.
                else
                   expon=.5*xx2/cnorm2
                endif
                fff=sfac*tam22(j-1)*exp(-expon*tam21(j-1))*dxm5(j)
     1            /dfloat(nintg-1)
                tam30(j)=tam30(j)+
     *            gg2*(xx/gg)**ix1(m)*fff
                tam23(j)=tam23(j)+
     *            gg2*(xx/gg)**ix2(m)*fff
                tam24(j)=tam24(j)+
     *            gg2*(xx/gg)**ix3(m)*fff
                tam25(j)=tam25(j)+
     *            gg2*(xx/gg)**ix4(m)*fff
                tam26(j)=tam26(j)+
     *            gg*(xx/gg)**ix1(m)*fff
                tam27(j)=tam27(j)+
     *            gg*(xx/gg)**ix2(m)*fff
                tam28(j)=tam28(j)+
     *            gg*(xx/gg)**ix3(m)*fff
                tam29(j)=tam29(j)+
     *            gg*(xx/gg)**ix4(m)*fff
 306          continue
 305        continue
            do 307 j=1,jthov2
              if (ialpha.eq.0) then
                alpha=(x(j)/x(jthov2))**2
              elseif (ialpha .eq. 1) then
                alpha=0.
              else
                alpha=1.
              endif
              tam20(j)=alpha*tam20(j)+(1.-alpha)*tam30(j)
              tam14(j)=alpha*tam14(j)+(1.-alpha)*tam24(j)
              tam16(j)=alpha*tam16(j)+(1.-alpha)*tam26(j)
              tam18(j)=alpha*tam18(j)+(1.-alpha)*tam28(j)
              tam13(j)=alpha*tam13(j)+(1.-alpha)*tam23(j)
              tam15(j)=alpha*tam15(j)+(1.-alpha)*tam25(j)
              tam17(j)=alpha*tam17(j)+(1.-alpha)*tam27(j)
              tam19(j)=alpha*tam19(j)+(1.-alpha)*tam29(j)
 307        continue

 308        continue

            do 310 j=2,jx
              jj=jx+1-j
              jp=jj+1
              jm=j-1
              tam2(jj)=tam2(jp)+tam20(jp)
              tam3(j) =tam3(jm)+tam13(j)
              tam4(jj)=tam4(jp)+tam14(jp)
              tam5(j) =tam5(jm)+tam15(j)
              tam6(jj)=tam6(jp)+tam16(jp)
              tam7(j) =tam7(jm)+tam17(j)
              tam8(jj)=tam8(jp)+tam18(jp)
              tam9(j) =tam9(jm)+tam19(j)
 310        continue

c..................................................................
c     tam2= v**(2+m) * M_m ; tam3= v**(1-m) * N_m ; tam4= v**m * R_m
c     tam5= v**(1-m) * E_m
c     tam6= v**(2+m) * Mh_m ; tam7= v**(1-m) * Nh_m ; tam8= v**m * Rh_m
c     tam9= v**(1-m) * Eh_m ; where Mh="M / gamma(ksi)", etc.
c.......................................................................
            do 330 j=1,jx
              tam2(j)=xm(j,ix5(m))*tam2(j)
              tam3(j)=xm(j,ix6(m))*tam3(j)
              tam4(j)=xm(j,ix7(m))*tam4(j)
              tam5(j)=xm(j,ix8(m))*tam5(j)
              tam6(j)=xm(j,ix5(m))*tam6(j)
              tam7(j)=xm(j,ix6(m))*tam7(j)
              tam8(j)=xm(j,ix7(m))*tam8(j)
              tam9(j)=xm(j,ix8(m))*tam9(j)
 330        continue

c..................................................................
c     sg(j) = B_m_b, needed for Rosenbluth potential g (g,h with relav corr)
c     sh(j) = A_m_b, needed for Rosenbluth potential h (h "=" g/gamma_prime)
c     sgx(j) = dsg/dv (not du), etc. for sgxx, ...
c.......................................................................
            do 340 j=2,jx
              sg(j)=cog(m,1)*(tam5(j)+tam2(j))-
     *          cog(m,2)*(tam3(j)+tam4(j))
              sgx(j)=cog(m,3)*tam2(j)-cog(m,4)*tam5(j)-
     *          cog(m,5)*tam4(j)+cog(m,6)*tam3(j)
              sgx(j)=gamma(j)*xi(j)*sgx(j)
              sgxx(j)=cog(m,7)*(tam5(j)+tam2(j))-
     *          cog(m,8)*(tam3(j)+tam4(j))
              sgxx(j)=gamsqr(j)*x2i(j)*sgxx(j)
              sh(j)=cog(m,1)*(tam9(j)+tam6(j))-
     *          cog(m,2)*(tam7(j)+tam8(j))
              shx(j)=cog(m,3)*tam6(j)-cog(m,4)*tam9(j)-
     *          cog(m,5)*tam8(j)+cog(m,6)*tam7(j)
              shx(j)=gamma(j)*xi(j)*shx(j)
              shxx(j)=cog(m,7)*(tam9(j)+tam6(j))-
     *          cog(m,8)*(tam7(j)+tam8(j))
              shxx(j)=gamsqr(j)*x2i(j)*shxx(j)
              shxxx(j)=cog(m,9)*tam6(j)-cog(m,10)*tam9(j)-
     *          cog(m,11)*tam8(j)+cog(m,12)*tam7(j)
              shxxx(j)=gamcub(j)*x3i(j)*shxxx(j)
 340        continue

c.......................................................................
c     compute A_a, B_a, ..., F_a as in GA report GA-A20978 p.11, 
c     Nov. 1992  (CQL3D Manual, Harvey and McCoy, IAEA TCM Montreal):
c     Mildly relativistic coll FP coeffs per M. Franz, UCRL-96510, 1987.
c.......................................................................
            fmmp1=m*(m+1)
            do 350 j=1,jx
              tam2(j)=-gamcub(j)*xi(j)*fmmp1*sh(j)
     *          +.5*gamsqr(j)*(2.+fmmp1)*shx(j)
     *          -gammi(j)*x(j)*shxx(j)
     *          -.5*gamm2i(j)*xsq(j)*shxxx(j)
              tam3(j)=.5*xsq(j)*sgxx(j)
              tam4(j)=.5*gamma(j)*(sgx(j)-gamma(j)*xi(j)*sg(j))
              tam5(j)=.5*gamcub(j)*x2i(j)*fmmp1*sh(j)
     *          -gamsqr(j)*xi(j)*shx(j)
     *          -.5*gammi(j)*shxx(j)
              tam6(j)=.5*gamma(j)*xi(j)*sgx(j)
              tam7(j)=.5*gamsqr(j)*x2i(j)*sg(j)
 350        continue

c$$$c     sum over m
c$$$            do 380 iii=1,imax(ileff,lr_)
c$$$              do 370 ii=0,1
c$$$                if (madd.eq.2 .and. ii.eq.0) goto 370
c$$$                i=iii*ii-(iy+1-iii)*(ii-1)
c$$$                do 360 j=2,jx
c$$$                  ca(i,j)=ca(i,j)+ss(i,ileff,m,lr_)*tam2(j)*gamefac(j)
c$$$                  cb(i,j)=cb(i,j)+ss(i,ileff,m,lr_)*tam3(j)*gamefac(j)
c$$$                  cc(i,j)=cc(i,j)+ssy(i,ileff,m,lr_)*tam4(j)*gamefac(j)
c$$$                  cd(i,j)=cd(i,j)+sinz(i,ileff,lr_)*
c$$$     *              ssy(i,ileff,m,lr_)*tam5(j)*gamefac(j)
c$$$                  ce(i,j)=ce(i,j)+sinz(i,ileff,lr_)*
c$$$     *              ssy(i,ileff,m,lr_)*tam4(j)*gamefac(j)
c$$$                  cf(i,j)=cf(i,j)+sinz(i,ileff,lr_)*
c$$$     *              (ss(i,ileff,m,lr_)*tam6(j)+ssyy(i,ileff,m,lr_)
c$$$     +              *tam7(j))*gamefac(j)
c$$$ 360            continue
c$$$ 370          continue
c$$$ 380        continue
c$$$c     end of loop over m
c     sum over m: Add contribution from each m

              do 380 iii=1,imax(ileff,lr_)
                if(iii.ge.itl+1 .and. mod(m,2).eq.1) goto 380 !YuP
                !YuP-111202: This removes a bug in the calculation of
                !YuP-111202: collisional contribution to C,D, and F.
                !YuP-111202: Check YuP Email to BH, 111201
                !YuP-111202: 
                !In trap region: no contribution from m=1,3,5...
                !The contribution from m=0,2,4,... will provide proper
                !parity in theta0-(pi/2): even parity for ca,cb,cf; 
                !odd parity for cc,cd,ce (they are ~ dPm/dtheta). 
                !No further symmetrization needed.
                do 370 ii=0,1
                i=iii*ii-(iy+1-iii)*(ii-1) ! i=iy+1-iii or i=iii
                do 360 j=2,jx
                  ca(i,j)=ca(i,j)+ss(i,ileff,m,lr_)*tam2(j)*gamefac(j)
                  cb(i,j)=cb(i,j)+ss(i,ileff,m,lr_)*tam3(j)*gamefac(j)
                  cc(i,j)=cc(i,j)+ssy(i,ileff,m,lr_)*tam4(j)*gamefac(j)
                  cd(i,j)=cd(i,j)+sinz(i,ileff,lr_)*
     *              ssy(i,ileff,m,lr_)*tam5(j)*gamefac(j)
                  ce(i,j)=ce(i,j)+sinz(i,ileff,lr_)*
     *              ssy(i,ileff,m,lr_)*tam4(j)*gamefac(j)
                  cf(i,j)=cf(i,j)+sinz(i,ileff,lr_)*
     *              (ss(i,ileff,m,lr_)*tam6(j)+ssyy(i,ileff,m,lr_)
     +              *tam7(j))*gamefac(j)
 360            continue ! j
 370          continue ! ii
 380          continue ! iii

 400      continue

          if (madd.eq.2 .or. symtrap.ne."enabled") goto 430

c     symmetrize in trap region
          do 420 i=itl+1,iyh
            iu=iy+1-i
            do 410 j=2,jx
              ca(i,j)=.5*(ca(i,j)+ca(iu,j))
              cb(i,j)=.5*(cb(i,j)+cb(iu,j))
              cf(i,j)=.5*(cf(i,j)+cf(iu,j))
              xq=sign(half,cc(i,j))
              xr=sign(half,cd(i,j))
              xs=sign(half,ce(i,j))
              cd(i,j)=xr*(abs(cd(i,j))+abs(cd(iu,j)))
              cc(i,j)=xq*(abs(cc(i,j))+abs(cc(iu,j)))
              ce(i,j)=xs*(abs(ce(i,j))+abs(ce(iu,j)))
              ca(iu,j)=ca(i,j)
              cb(iu,j)=cb(i,j)
              cc(iu,j)=-cc(i,j)
              cd(iu,j)=-cd(i,j)
              ce(iu,j)=-ce(i,j)
              cf(iu,j)=cf(i,j)
 410        continue
 420      continue
 430      continue

c.......................................................................
c     add contribution to each gen. species A_kk,.., including charge,
c     ln(Lambda) and mass coefficients.
c.......................................................................

          do 490 kk=1,ngen
            if (colmodl.eq.4) then
              if (kk.eq.kelecg .and. k.eq.kelecg) goto 490
              if (kk.ne.kelecg .and. k.eq.ngen) goto 490
            endif
            anr1=gama(kk,k)*satioz2(k,kk)*one_
            if (anr1.lt.em90) goto 490
CPTR>>>REPLACE PTR-DSCALCACD
            call dscal(iyjx,anr1,ca(1,1),1) !-YuP: size of ca..cf: iy*jx
            call dscal(iyjx,anr1,cb(1,1),1)
            call dscal(iyjx,anr1,cc(1,1),1)
            call dscal(iyjx,anr1,cd(1,1),1)
            call dscal(iyjx,anr1,ce(1,1),1)
            call dscal(iyjx,anr1,cf(1,1),1)
            call dscal(iyjx,satiom(kk,k),ca(1,1),1)
            call dscal(iyjx,satiom(kk,k),cd(1,1),1)
CPTR<<<END PTR-DSCALCACD

            
c.......................................................................
c     Note: At this point, ca, ..,cf(i,j) are the coeff. from gen. species k
c     at a given orbit position l.
c.......................................................................

            if (cqlpmod .ne. "enabled") then

c     Perform the bounce averaging
              do 480 i=1,imax(l,lr_)
                ii=iy+1-i
                ax=abs(coss(i,l_))*dtau(i,l,lr_)
                ay=tot(i,l,lr_)/sqrt(bbpsi(l,lr_))
                az=ay*tot(i,l,lr_)

cBH091031                ax1=ax
cBH091031                !i.e, not bounce pt interval:
cBH091031                if (l.eq.lz .or. lmax(i,lr_).ne.l) goto 440
cBH091031                ax1=ax+dtau(i,l+1,lr_)*abs(coss(i,l_))
cBH091031 440            continue
                if (eqsym.ne."none") then       !i.e. up-down symm
                   !if not bounce interval
                   if(l.eq.lz .or. l.ne.lmax(i,lr_)) then
                      ax1=ax
                   else !bounce interval: additional contribution
                      ax1=ax+dtau(i,l+1,lr_)*abs(coss(i,l_))
                   endif
                else  !eqsym="none"
                   if (l.lt.lz_bmax(lr_) .and. l.eq.lmax(i,lr_))then
                      !trapped, with tips between l and l+1 (above midplane)
                      ax1=ax+dtau(i,l+1,lr_)*abs(coss(i,l_)) 
                      !-YuP  Note: dtau(i,l+1,lr_)=0
                   elseif (l.gt.lz_bmax(lr_) .and. l.eq.lmax(i+iyh,lr_))
     +                then 
                      !trapped, with tips between l and l-1 (below midplane)     
                      ax1=ax+dtau(i,l-1,lr_)*abs(coss(i,l_)) !NB:l-1
                      !-YuP  Note: dtau(i,l-1,lr_)=0
                   else 
                      !passing (i<itl), or trapped but with tips at other l;
                      !also, at l=lz_bmax, includes last trapped particle i=itl
                      !(for such particle, lmax(itl)=lz_bmax; see micxinil)
                      ax1=ax
                   endif
                endif

                do 450 j=2,jx
                  cal(i,j,kk,l_)=cal(i,j,kk,l_)+ax1*ca(i,j)
                  cbl(i,j,kk,l_)=cbl(i,j,kk,l_)+ax1*cb(i,j)
                  ccl(i,j,kk,l_)=ccl(i,j,kk,l_)+ax*tot(i,l,lr_)*cc(i,j)
                  cdl(i,j,kk,l_)=cdl(i,j,kk,l_)+ax*ay*cd(i,j)
                  cel(i,j,kk,l_)=cel(i,j,kk,l_)+ax*ay*ce(i,j)
                  cfl(i,j,kk,l_)=cfl(i,j,kk,l_)+ax*az*cf(i,j)
 450            continue
                if (madd.eq.2) goto 470
                do 460 j=2,jx
                  cal(ii,j,kk,l_)=cal(ii,j,kk,l_)+ax1*ca(ii,j)
                  cbl(ii,j,kk,l_)=cbl(ii,j,kk,l_)+ax1*cb(ii,j)
                  ccl(ii,j,kk,l_)=ccl(ii,j,kk,l_)+ax*tot(ii,l,lr_)
     1              *cc(ii,j)
                  cdl(ii,j,kk,l_)=cdl(ii,j,kk,l_)+ax*ay*cd(ii,j)
                  cel(ii,j,kk,l_)=cel(ii,j,kk,l_)+ax*ay*ce(ii,j)
                  cfl(ii,j,kk,l_)=cfl(ii,j,kk,l_)+ax*az*cf(ii,j)
 460            continue
 470            continue
 480          continue ! i=1,imax(l,lr_)

            else ! cqlpmod = "enabled"
              do 485 i=1,iy
                do 486 j=2,jx
                  cal(i,j,kk,l_)=cal(i,j,kk,l_)+ca(i,j)
                  cbl(i,j,kk,l_)=cbl(i,j,kk,l_)+cb(i,j)
                  ccl(i,j,kk,l_)=ccl(i,j,kk,l_)+cc(i,j)
                  cdl(i,j,kk,l_)=cdl(i,j,kk,l_)+cd(i,j)
                  cel(i,j,kk,l_)=cel(i,j,kk,l_)+ce(i,j)
                  cfl(i,j,kk,l_)=cfl(i,j,kk,l_)+cf(i,j)
 486            continue
 485          continue
            endif

CPTR>>>REPLACE PTR-DSCAL2
            fscal= one/anr1
            call dscal(iyjx,fscal,ca(1,1),1)
            call dscal(iyjx,fscal,cb(1,1),1)
            call dscal(iyjx,fscal,cc(1,1),1)
            call dscal(iyjx,fscal,cd(1,1),1)
            call dscal(iyjx,fscal,ce(1,1),1)
            call dscal(iyjx,fscal,cf(1,1),1)
            call dscal(iyjx,one/satiom(kk,k),ca(1,1),1)
            call dscal(iyjx,one/satiom(kk,k),cd(1,1),1)
CPTR<<<END PTR-DSCAL2

 490      continue ! kk=1,ngen

c     end of loop over gen. species k
 500    continue ! k=1,ngen

ccc      kk=1
ccc      do i=1,iy
ccc         prnt1(i)=ca(i,2)
ccc         prnt2(i)=cb(i,2)
ccc         prnt3(i)=cc(i,2)
ccc         prnt4(i)=cd(i,2)
ccc         prnt5(i)=ce(i,2)
ccc         prnt6(i)=cf(i,2)
ccc      enddo
      lll=l_  !to create point to stop


c     end of loop over orbit l
 600  continue

ccc      kk=1
ccc      do i=1,iy
ccc         prnt1(i)=cal(i,2,kk,l_)
ccc         prnt2(i)=cbl(i,2,kk,l_)
ccc         prnt3(i)=ccl(i,2,kk,l_)
ccc         prnt4(i)=cdl(i,2,kk,l_)
ccc         prnt5(i)=cel(i,2,kk,l_)
ccc         prnt6(i)=cfl(i,2,kk,l_)
ccc      enddo

 700  continue
c..................................................................
c     define needed coefficients at pass/trapped boundary
c..................................................................
      if (cqlpmod .ne. "enabled") then
        do 2001 k=1,ngen
          do 2002 j=1,jx
            cal(itl,j,k,l_)=0.25*vptb(itl,lr_) * ( cal(itl-1,j,k,l_)/
     /        vptb(itl-1,lr_)+2.*cal(itl+1,j,k,l_)/vptb(itl+1,lr_) +
     +        cal(itu+1,j,k,l_)/vptb(itu+1,lr_) )
            cbl(itl,j,k,l_)=0.25*vptb(itl,lr_) * ( cbl(itl-1,j,k,l_)/
     /        vptb(itl-1,lr_)+2.*cbl(itl+1,j,k,l_)/vptb(itl+1,lr_) +
     +        cbl(itu+1,j,k,l_)/vptb(itu+1,lr_) )
            cal(itu,j,k,l_)=cal(itl,j,k,l_)
            cbl(itu,j,k,l_)=cbl(itl,j,k,l_)
 2002     continue
 2001   continue
      endif
c..................................................................
      if (madd .eq. 2) call cfpsymt
      

      do 2000 k=1,ngen
        fscal= one/tnorm(k)
        call dscal(iyjx,fscal,cal(1,1,k,l_),1)
        call dscal(iyjx,fscal,cbl(1,1,k,l_),1)
        call dscal(iyjx,fscal,ccl(1,1,k,l_),1)
        call dscal(iyjx,fscal,cdl(1,1,k,l_),1)
        call dscal(iyjx,fscal,cel(1,1,k,l_),1)
        call dscal(iyjx,fscal,cfl(1,1,k,l_),1)
        call dscal(iyjx,fscal,eal(1,1,k,1,l_),1)
        call dscal(iyjx,fscal,ebl(1,1,k,1,l_),1)
        call dscal(iyjx,fscal,eal(1,1,k,2,l_),1)
        call dscal(iyjx,fscal,ebl(1,1,k,2,l_),1)

c..................................................................
c     For the case that colmodl=3, a positive definite operator
c     is not guaranteed. This is a bastardized, hybrid collisional
c     model (see the input information in the input deck or the
c     user manual) and should be used with care. Caveat emptor!
c     In any case if we get negative diffusion from the model,
c     it is set to zero to keep the code from blowing up.
c..................................................................
        if (colmodl.eq.3) then
          do 2004 j=1,jx
            do 2005 i=1,iy
              if(cbl(i,j,k,l_).lt.0.) then
                 cbl(i,j,k,l_)=em100
              endif
              if(cfl(i,j,k,l_).lt.0.) then
                 cfl(i,j,k,l_)=em100
              endif
 2005       continue
 2004     continue
        endif
 2000 continue

c Bug (from Gary Kerbel, Oct 31, 2005) needs checking:
c Take difference of eal (linear part for electron + 2nd part from ions)
c and cal (sum of linear and nonlinear) ==> just NL component, species
c scattering from itself.  If Maxwl, NL is same as L.  In limit of
c small grid, (cal-eal)/eal should equal 1.
c But it doesn't.  Check kerbel email, Oct 31, 2005.

      return
      end
