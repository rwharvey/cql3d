c
c
      subroutine cfpcoefr
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c..................................................................
c     Subroutine to calculate bounce-averaged Fokker-Planck collision
c     coefficients. Used for relativ='fully'. See Mark Franz thesis.
c     If (cqlpmod .eq. "enabled") then compute only the coefficients
c     at the orbit position l=l_ and do not perform the bounce-averages.
c..................................................................
      include 'param.h'
      include 'comm.h'
c
c..................................................................
c     These represent the derivatives of the function gamma^n alpha^m
c170403: Changed index n of gman, etc. statement functions to ng, to avoid
c170403: a recent intel compiler problem which showed up for John Wright.
c..................................................................
      gman(m,ng,j)=gamman(j,m)*alphan(j,ng) ! gamma^m * alpha^-ng
      gmans(m,ng,j)=gman(m,ng,j)*asnha(j)   ! gamma^m * alpha^-ng * ln(alpha+gamma)
      
      ! First derivative: d(gamma^m * alpha^-n) / d(v/vnorm)
!!!      dgawv(m,n,j)=gman(m,n,j)*(m*alphan(j,-1)*gamman(j,-2)*cnormi-
!!!     *  n*alphan(j,1)/cnorm)      
      dgawv(m,ng,j)=gman(m,ng,j)*gamman(j,-2)*alphan(j,1)*
     *  ( (m-ng)*alphan(j,-2) - ng )*cnormi 
      !-YuP-> New version: rearranged to improve accuracy
     
      ! Second derivative: d^2(gamma^m * alpha^-n) / d(v/vnorm)^2
!!!      dgawvv(m,n,j)=gman(m,n,j)*(
!!!     *  m*(m-2)*alphan(j,-2)*gamman(j,-2)*gamman(j,-2)+
!!!     *  m*(1-2*n)*gamman(j,-2)+
!!!     *  n*(n+1)*alphan(j,2) )*cnorm2i
      dgawvv(m,ng,j)=gman(m,ng,j)*alphan(j,2)*gamman(j,-2)*gamman(j,-2)*
     +  (  (ng-m)*(ng-m+1)*alphan(j,-2)*alphan(j,-2) +
     +     (2*ng*(ng-m+1)+m)*alphan(j,-2) + ng*(ng+1)  )*cnorm2i
      !-YuP-> New version: rearranged to improve accuracy
     
      ! First derivative: d(gamma^m * alpha^-n * ln(alpha+gamma)) / d(v/vnorm)     
      dgaswv(m,ng,j)=dgawv(m,ng,j)*asnha(j)+
     *  gman(m,ng,j)*gamman(j,-1)*cnormi
     
      ! Second derivative: d^2(gamma^m * alpha^-n * ln(alpha+gamma)) / d(v/vnorm)^2
!!!      dgaswvv(m,n,j)=dgawvv(m,n,j)*asnha(j)+
!!!     *  2.*dgawv(m,n,j)*gamman(j,-1)*cnormi-
!!!     *  gman(m,n,j)*alphan(j,-1)*gamman(j,-2)*gamman(j,-1)*cnorm2i
      dgaswvv(m,ng,j)=dgawvv(m,ng,j)*asnha(j) +
     +  2.*alphan(j,1)*gamman(j,-2)*gamman(j,-1)*
     * ( (m-ng-0.5)*alphan(j,-2) - ng )*gman(m,ng,j)*cnorm2i
      !-YuP-> New version: rearranged to improve accuracy
      !
      ! Following the message from John Wright[04-2017, 04/03/2017] :
      ! Changed index 'n' in the above statement functions
      ! to 'ng'.  Index 'n' is also in comm.h, which causes 
      ! a compilation error (bug) in new Intel compiler (2017)
      !
c..................................................................
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
cBH080215 Following cfpcoefn example
        nintg0=41.*max(1.d0,sqrt(fmass(k)/1.67e-24))
cBH080215        nintg=max0(21,min0(51,int(x(2)/sqrt(.2*cnorm2/reltmp))))
        nintg=max(nintg0,min0(51,int(x(2)/sqrt(.2*cnorm2/reltmp))))
        nintg=2*(nintg/2)+1
        do 10 i=1,nintg
          sfac=4.+2.*(2*(i/2)-i)
          if(i.eq.1 .or. i.eq.nintg) sfac=1.
          sfac=sfac/3.
          do 11 j=2,jx
            delx=dxm5(j)/dfloat(nintg-1)
            xx=x(j)-(i-1)*delx
            xxc=xx/cnorm
            xxc2=xxc*xxc
                    if(cnorm2i-em8 .gt. 0.d0) then
                       asalp=log(xxc+sqrt(xxc2+1.))
                    else
                       asalp=xxc
                    endif 
            xx2=xx*xx
            gg=sqrt(1.+xx2*cnorm2i)
                    if(cnorm2i*xx2-em8 .gt. 0.d0) then
                      expn=gg-1.
                    else
                      expn=.5*xx2/cnorm2
                    endif
            fff=sfac*rnorm*cnorm2*exp(-expn*reltmp)*delx
            temp1(1,j) =temp1(1,j)+fff*xxc2/gg
            temp1(2,j) =temp1(2,j)+fff*xxc2
            temp1(3,j) =temp1(3,j)+xxc*xxc2*fff/gg
            temp1(4,j) =temp1(4,j)+fff*xxc
            temp1(5,j) =temp1(5,j)+gg*fff*xxc
            temp1(6,j) =temp1(6,j)+fff*xxc/gg
            temp1(7,j) =temp1(7,j)+asalp*xxc2*xxc*fff/gg
            temp1(8,j) =temp1(8,j)+asalp*fff*xxc
            temp1(9,j) =temp1(9,j)+asalp*gg*fff*xxc
            temp1(10,j)=temp1(10,j)+asalp*fff*xxc/gg
            temp1(11,j)=temp1(11,j)+asalp*fff*xxc2
 11       continue
 10     continue
c..................................................................
        tam1(1)=0.
        tam2(1)=0.
        tam3(jx)=0.
        tam4(jx)=0.
        tam5(jx)=0.
        tam6(jx)=0.
        tam7(1)=0.
        tam8(1)=0.
        tam9(1)=0.
        tam10(1)=0.
        tam11(jx)=0.
c......................................................................:
c     tam3(jx)->tam6(jx), tam11(jx) and tam12(jx) represent integrals
c     from xmax of the grid to infinity.
c     The following coding performs this integration. In this case 21*xmax
c     represents infinity. The lack of this piece is most obvious when
c     ions are a general species and electrons are fixed Maxwellians.
c..................................................................
        do 15 ll2=1,21
          sfac=4.+2.*(2*(ll2/2)-ll2)
          if( ll2.eq.1 .or. ll2.eq.21) sfac=1.
          sfac=sfac*x(jx)/60.
          do 16 ll1=1,20
            xx=(realiota(ll1)+.05*realiota(ll2-1))*x(jx)
            xxc=xx/cnorm
            xxc2=xxc*xxc
                    if(cnorm2i-em8 .gt. 0.d0) then
                       asalp=log(xxc+sqrt(xxc2+1.))
                    else
                       asalp=xxc
                    endif 
            xx2=xx*xx
            gg=sqrt(1+xx2*cnorm2i)
                    if(cnorm2i*xx2-em8 .gt. 0.d0) then
                      expon=gg-1.
                    else
                      expon=.5*xx2/cnorm2
                    endif
            fff=sfac*rnorm*cnorm2*exp(-expon*reltmp)
            tam3(jx)=tam3(jx)+xxc*xxc2*fff/gg
            tam4(jx)=tam4(jx)+fff*xxc
            tam5(jx)=tam5(jx)+gg*fff*xxc
            tam6(jx)=tam6(jx)+fff*xxc/gg
            tam11(jx)=tam11(jx)+asalp*fff*xxc2
 16       continue
 15     continue
c......................................................................:
        do 20 j=2,jx
          jj=jx+1-j
          jp=jj+1
          jm=j-1
          tam1(j)  =tam1(jm)+temp1(1,j)
          tam2(j)  =tam2(jm)+temp1(2,j)
          tam3(jj) =tam3(jp)+temp1(3,jp)
          tam4(jj) =tam4(jp)+temp1(4,jp)
          tam5(jj) =tam5(jp)+temp1(5,jp)
          tam6(jj) =tam6(jp)+temp1(6,jp)
          tam7(j)  =tam7(jm)+temp1(7,j)
          tam8(j)  =tam8(jm)+temp1(8,j)
          tam9(j)  =tam9(jm)+temp1(9,j)
          tam10(j) =tam10(jm)+temp1(10,j)
          tam11(jj)=tam11(jp)+temp1(11,jp)
 20     continue
c..................................................................
        do 30 j=2,jx
          temp1(1,j)=(gman(2,1,j)*(2.*tam9(j)-tam2(j))+
     *      gman(1,0,j)*(4.*tam11(j)-tam5(j)-tam3(j))+
     *      gman(0,-1,j)*(2.*tam7(j)-tam2(j))-
     *      gman(0,1,j)*tam10(j)+
     *      gmans(2,1,j)*2.*tam5(j)+
     *      gmans(1,0,j)*4.*tam2(j)+
     *      gmans(0,-1,j)*2.*tam3(j)-
     *      gmans(0,1,j)*tam6(j))/4.
          temp1(2,j)=(dgawv(2,1,j)*(5.*tam2(j)-2.*tam9(j))+
     *      dgawv(1,0,j)*(5.*(tam5(j)+tam3(j))-4.*tam11(j))+
     *      dgawv(0,-1,j)*(5.*tam2(j)-2.*tam7(j))-
     *      dgawv(0,1,j)*3.*tam10(j)-
     *      dgaswv(2,1,j)*2.*tam5(j)-
     *      dgaswv(1,0,j)*4.*tam2(j)-
     *      dgaswv(0,-1,j)*2.*tam3(j)-
     *      dgaswv(0,1,j)*3.*tam6(j))/8.
          temp1(3,j)=(dgawvv(2,1,j)*(5.*tam2(j)-2.*tam9(j))+
     *      dgawvv(1,0,j)*(5.*(tam5(j)+tam3(j))-4.*tam11(j))+
     *      dgawvv(0,-1,j)*(5.*tam2(j)-2.*tam7(j))-
     *      dgawvv(0,1,j)*3.*tam10(j)-
     *      dgaswvv(2,1,j)*2.*tam5(j)-
     *      dgaswvv(1,0,j)*4.*tam2(j)-
     *      dgaswvv(0,-1,j)*2.*tam3(j)-
     *      dgaswvv(0,1,j)*3.*tam6(j))/8.
          temp1(4,j)=dgawv(1,1,j)*(2.*tam1(j)-tam8(j))-
     *      dgaswv(1,1,j)*tam4(j)-
     *      dgaswv(0,0,j)*tam1(j)
 30     continue

c.......................................................................
c     Perform the bounce-average for the background species and introduce
c     the contribution to all general species coeff.
c.......................................................................

        if (cqlpmod .ne. "enabled") then
          call bavdens(k)
        else
          do 59 i=1,iy
            bavdn(i,l_)=denpar(k,ls_)
            bavpd(i,l_)=denpar(k,ls_)*sinn(i,l_)
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
          anr1=4.*pi*cnorm*gama(kk,k)*satioz2(k,kk)*one_
          anr2=anr1*satiom(kk,k)/cnorm2
          call bcast(tem1,zero,iyjx)
          call bcast(tem2,zero,iyjx)
          do 70 j=2,jx
            ttta=-anr2*gamman(j,1)*xsq(j)*temp1(4,j)*gamefac(j,k) !YuP[2019-07-26] k index added
            tttb=anr1*xsq(j)*
     *        (gamman(j,3)*temp1(3,j)+
     *        gamman(j,1)*x(j)*temp1(2,j)*cnorm2i+
     *        .5*gamman(j,1)*temp1(1,j)*cnorm2i)*gamefac(j,k) !YuP[2019-07-26] k index added
            tttf=anr1*(gman(1,1,j)*temp1(2,j)/cnorm+
     *        .5*gamman(j,-1)*temp1(1,j)*cnorm2i)*gamefac(j,k) !YuP[2019-07-26] k index added
            do 60 i=1,iy
              jj=i+(j-1)*iy 
              tem1(jj)=ttta*vptb(i,lr_)
     1          *bavdn(i,lr_)
              tem2(jj)=tttb*vptb(i,lr_)
     1          *bavdn(i,lr_)
              cal(i,j,kk,l_)=cal(i,j,kk,l_)+tem1(jj)
              cbl(i,j,kk,l_)=cbl(i,j,kk,l_)+tem2(jj)
              cfl(i,j,kk,l_)=cfl(i,j,kk,l_)+tttf*vptb(i,lr_)
     1          *bavpd(i,lr_)
 60         continue
 70       continue
 
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

 80     continue
 100  continue

 110  continue

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
        iorbend=lz
      else
        iorbstr=l_
        iorbend=l_
      endif
      do 600 l=iorbstr,iorbend
        ileff=l
        if (cqlpmod .eq. "enabled") ileff=ls_

        do 500 k=1,ngen
c..................................................................
c     Jump out if this species is not to be used as a background species
c..................................................................
          if (kfield(k).eq."disabled") go to 500

c     zero ca, cb,.., cf  : 
          call bcast(ca,zero,iyjx)
          call bcast(cb,zero,iyjx)
          call bcast(cc,zero,iyjx)
          call bcast(cd,zero,iyjx)
          call bcast(ce,zero,iyjx)
          call bcast(cf,zero,iyjx)

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
c     coeff. of Legendre decomposition of f

            call cfpleg(m,ileff,1)

c.................................................................
c     Calculate the integrals prior to using them
c..................................................................
            tam2(1) = 0. !-YuP-> added: 
            tam3(1) = 0. !-YuP-> clean-up tam2&3(1) from previous usage.
            ! Note: Generally, for l2>2, tam2&3(1) can go to inf. (v->0).
            ! But we assume that at v->0 all Legendre coeffs. are zero except m=0.
            ! Then, we are only interested in l2=-2,-1,0,+1. 
            ! For such l2, tam2&3(1)->0
            do 315 l1=-1,m+2   !!! YuP-> Was: 0,m+2
              do 316 l2=-1+l1,m+1
                do 301 j=2,jx
                  jm=j-1
                  tam2(j)=.5*dxm5(j)*
     *              (gman(l1,l2,j)*xsq(j)*tam1(j)*gamman(j,-1)+
     *              gman(l1,l2,jm)*xsq(jm)*tam1(jm)*gamman(jm,-1))
                  tam3(j)=.5*dxm5(j)*
     *              (gmans(l1,l2,j)*xsq(j)*tam1(j)*gamman(j,-1)+
     *              gmans(l1,l2,jm)*xsq(jm)*tam1(jm)*gamman(jm,-1))
 301            continue

c     note: ialpha set to 2 at top of subroutine
                if (m.ge.1 .or. ialpha.eq.2) go to 308
c..................................................................
c     The do loops 302 through 307 seek to take advantage of the
c     Maxwellian nature of f for v < vth/2. This is employed in
c     the integrations to obtain the functionals - Thus f is
c     assumed to be Maxwellian between velocity mesh points,
c     not linear.
c..................................................................
                do 302 j=2,jx
                  xs=sqrt(0.5*temp(k,lr_)*ergtkev/fmass(k))
                  if (cqlpmod .eq. "enabled")
     +              xs=sqrt(temppar(k,ls_)*ergtkev*0.5/fmass(k))
                  if (x(j)*vnorm.gt.xs) go to 303
 302            continue
 303            continue
                jthov2=j-1
c..................................................................
c     Determine the "local" Maxwellian between meshpoints.
c..................................................................
                do 304 j=1,jthov2
c990131                  tam6(j)=alog(tam1(j)/tam1(j+1))/
                  tam6(j)=log(tam1(j)/tam1(j+1))/
     *              (gamm1(j+1)-gamm1(j))
                  tam7(j)=tam1(j)/exp(-(gamm1(j)*tam21(j)))
 304            continue
                rstmss=fmass(k)*clite2/ergtkev
                reltmp=rstmss/temp(k,lr_)
                if (cqlpmod .eq. "enabled") reltmp=rstmss/temppar(k,ls_)
                call bcast(tam8,zero,jx*2)
                nintg=max0(21,min0(51,int(x(2)/
     *            sqrt(.2*cnorm2/reltmp))))
                nintg=2*(nintg/2)+1
                do 305 i=1,nintg
                  sfac=4.+2.*(2*(i/2)-i)
                  if(i.eq.1 .or. i.eq.nintg) sfac=1.
                  sfac=sfac/3.
                  do 306 j=2,jthov2
                    xx=x(j)-(i-1)*dxm5(j)/dfloat(nintg-1)
                    xxc=xx/cnorm
                    xxc2=xxc*xxc
                    if(cnorm2i-em8 .gt. 0.d0) then
                       asalp=log(xxc+sqrt(xxc2+1.))
                    else
                       asalp=xxc
                    endif 
                    xx2=xx*xx
                    gg=sqrt(1.+xx2*cnorm2i)
                    if(cnorm2i*xx2-em8 .gt. 0.d0) then
                      expon=gg-1.
                    else
                      expon=.5*xx2/cnorm2
                    endif
                    fff=cnorm2*sfac*tam7(j-1)*exp(-expon*tam6(j-1))*
     *                dxm5(j)/dfloat(nintg-1)
                    tam8(j)=tam8(j)+gg**(l1-1)*xxc**(-l2+2)*fff
                    tam9(j)=tam9(j)+gg**(l1-1)*xxc**(-l2+2)*asalp*fff
 306              continue
 305            continue
                do 307 j=1,jthov2
                  if (ialpha.eq.0) then
                    alph=(x(j)/x(jthov2))**2
                  elseif (ialpha .eq. 1) then
                    alph=0.
                  else
                    alph=1.
                  endif
                  tam2(j)=alph*tam2(j)+(1.-alph)*tam8(j)
                  tam3(j)=alph*tam3(j)+(1.-alph)*tam9(j)
 307            continue
 308            continue
                tamt1(1,1,l2,l1)=0.
                tamt1(2,jx,l2,l1)=0.
                tamt2(1,1,l2,l1)=0.
                tamt2(2,jx,l2,l1)=0.
                do 310 j=2,jx
                  jj=jx+1-j
                  jp=jj+1
                  jm=j-1
                  ! X-integrals after Eq.(53) Franz Thesis
                  tamt1(1,j,l2,l1) =tamt1(1,jm,l2,l1)+tam2(j) 
                  tamt1(2,jj,l2,l1)=tamt1(2,jp,l2,l1)+tam2(jp)
                  tamt2(1,j,l2,l1) =tamt2(1,jm,l2,l1)+tam3(j)
                  tamt2(2,jj,l2,l1)=tamt2(2,jp,l2,l1)+tam3(jp)
 310            continue
 316          continue ! l2=-1+l1,m+1
 315        continue   ! l1=0,m+2
c..................................................................
c     In this set of loops I need to sum up all the contributions to the
c     integrals and then just use these in the coefficients
c..................................................................
            do 309 j=1,jx  !-YuP->: was j=2,jx
              tam1(j)=0.
              tam2(j)=0.
              tam3(j)=0.
              tam4(j)=0.
              tam5(j)=0.
              tam6(j)=0.
 309        continue
c
            do 330 l1=0,m/2
              do 331 l2=0,m-2*l1
                cons1=((-1.)**(l1+l2))*(2.**(-m))*choose(m,l1)*
cBH090808    *            choose(2*(m-l1),m)*choose(m-2*l1,l2)
cBH090808  Bounds error due to choose(2*(m-l1), ).  
cYP090808  Also, questionable looking df/dt in velocity space.
cYP090808  Notes also should have choose(2*(m-l1),l1), from Franz thesis
cBH090810  BUT, original choose(2*(m-l1),m) gives much more reasonable
cBH090810  looking conductivity, after fixing dims and setting of choose(,,).
cBH090810  So, sticking with original coding here until we finish with 
cBH090810  verification of the physics equations.
     *            choose(2*(m-l1),m)*choose(m-2*l1,l2)
c..................................................................
c     Do the j+1 and j+2 terms first
c..................................................................
                do 317 l3=0,l2+2
                  cons3=cons1*choose(l2+2,l3)/dfloat(l2+2)
                  ig0=m-2*l1-l3+1
                  ig1=ig0+1
                  ia0=m-2*l1-l3+1
                  ia1=ia0-1
                  ivn=mod(l3,2) + 1
                  idd=3-ivn
c..................................................................
c     Psi1:  Sigma1(j+1,m-2k-j,m-2k+1)
c     Psi3:  Sigma2(j+1,m-2k-j,m-2k+1)
c..................................................................
                  if(l3.le.l2+1) then
                    cons2=cons1*choose(l2+1,l3)/dfloat(l2+1)
                    do 318 j=2,jx
                      tam1(j)=tam1(j)+cons2*gman(ig0,ia0,j)*
     *                  ( gamman(j,1)*tamt1(ivn,j,ia1,ig0)+
     *                  alphan(j,-1)*tamt1(idd,j,ia0,ig1))
                      tam3(j)=tam3(j)-cons2*gman(ig0,ia0,j)*
     *                  (asnha(j)*tamt1(idd,j,ia0,ig0)+
     *                  tamt2(ivn,j,ia0,ig0))
                      tam4(j)=tam4(j)+cons2*
     *                  (dgawv(ig1,ia0,j)*tamt1(ivn,j,ia1,ig0)+
     *                  dgawv(ig0,ia1,j)*tamt1(idd,j,ia0,ig1))
                      tam5(j)=tam5(j)+cons2*
     *                  (dgawvv(ig1,ia0,j)*tamt1(ivn,j,ia1,ig0)+
     *                  dgawvv(ig0,ia1,j)*tamt1(idd,j,ia0,ig1))
                      tam6(j)=tam6(j)-cons2*
     *                  (dgaswv(ig0,ia0,j)*tamt1(idd,j,ia0,ig0)+
     *                  dgawv(ig0,ia0,j)*tamt2(ivn,j,ia0,ig0))
 318                continue
                  endif
c..................................................................
c     Psi1:  Sigma2(j+2,m-2k-j,m-2k+1)
c     Psi2:  Sigma2(j+2,m-2k-j,m-2k+1)
c..................................................................
                  do 319 j=2,jx
                    tam1(j)=tam1(j)-.5*cons3*gman(ig1,ia0,j)*
     *                (asnha(j)*tamt1(idd,j,ia0,ig1)+
     *                tamt2(ivn,j,ia0,ig1))
                    tam2(j)=tam2(j)+cons3*gman(ig1,ia0,j)*
     *                (asnha(j)*tamt1(idd,j,ia0,ig1)+
     *                tamt2(ivn,j,ia0,ig1))
                    tam4(j)=tam4(j)-.5*cons3*
     *                (dgaswv(ig1,ia0,j)*tamt1(idd,j,ia0,ig1)+
     *                dgawv(ig1,ia0,j)*tamt2(ivn,j,ia0,ig1))
                    tam5(j)=tam5(j)-.5*cons3*
     *                (dgaswvv(ig1,ia0,j)*tamt1(idd,j,ia0,ig1)+
     *                dgawvv(ig1,ia0,j)*tamt2(ivn,j,ia0,ig1))
 319              continue
 317            continue  ! l3= 0 : l2+2
c.......................................................................:
c     Now need to perform the Sigma 3 terms
c     Need to check for odd or even indicies
c.......................................................................:
c     Case 1:  l2 odd
c.......................................................................:
                if(mod(l2,2).eq.1) then
                  cons21=cons1*(fctrl(l2+1)*(l2+2)*.5**(l2+1))/
     *              (dfloat(l2+1)*fctrl((l2+1)/2)**2)
                  cons22=cons1*fctrl((l2+1)/2)**2*2.**(l2+1)/
     *              (dfloat(l2+2)*fctrl(l2+2))
                  cons23=cons22*(l2+3)/dfloat(2*(l2+1))
c..................................................................
c     Psi3:  Sigma3(j+1,m-2k-j,m-2k+1)=>Sigma2(0,m-2k-j,m-2k+1) term
c..................................................................
                  ig0=m-2*l1-l2
                  ia0=m-2*l1+1
                  do 320 j=2,jx
                    tam3(j)=tam3(j)+cons21*gman(ig0,ia0,j)*
     *                (asnha(j)*tamt1(2,j,ia0,ig0)+
     *                tamt2(1,j,ia0,ig0))
                    tam6(j)=tam6(j)+cons21*
     *                (dgaswv(ig0,ia0,j)*tamt1(2,j,ia0,ig0)+
     *                dgawv(ig0,ia0,j)*tamt2(1,j,ia0,ig0))
 320              continue
                  do 322 l3=0,(l2+1)/2
c..................................................................
c     Psi3:  Sigma3(j+1,m-2k-j,m-2k+1)
c..................................................................
                    if(l3.gt.0) then
                      cons31=cons21*fctrl(l3)*fctrl(l3-1)*
     *                  2.**(2*l3-1)/fctrl(2*l3)
                      do 321 l4=0,2*l3-1
                        cons41=cons31*choose(2*l3-1,l4)
                        ig0=m-2*l1-l2+2*l3-1-l4
                        ig1=ig0+1
                        ia0=m-2*l1-l4+1
                        ia1=ia0-1
                        ivn=mod(l4,2)+1
                        idd=3-ivn
                        do 3211 j=2,jx
                          tam3(j)=tam3(j)+cons41*gman(ig0,ia0,j)*
     *                      ( gamman(j,1)*tamt1(ivn,j,ia1,ig0)+
     *                      alphan(j,-1)*tamt1(idd,j,ia0,ig1))
                          tam6(j)=tam6(j)+cons41*
     *                      (dgawv(ig1,ia0,j)*tamt1(ivn,j,ia1,ig0)+
     *                      dgawv(ig0,ia1,j)*tamt1(idd,j,ia0,ig1))
3211                   continue
 321                  continue ! l4= 0 : 2*l3-1
                    endif ! if(l3.gt.0)
c..................................................................
c     Psi1:  Sigma3(j+2,m-2k-j,m-2k+1)=>Sigma1(j+2,m-2k-j,m-2k+1) term
c     Psi2:  Sigma3(j+2,m-2k-j,m-2k+1)=>Sigma1(j+2,m-2k-j,m-2k+1) term
c..................................................................
                    cons31=cons22*fctrl(2*l3)*.5**(2*l3)/fctrl(l3)**2
                    cons32=cons23*fctrl(2*l3)*.5**(2*l3)/fctrl(l3)**2
                    do 3221 l4=0,2*l3
                      cons41=cons31*choose(2*l3,l4)
                      cons42=cons32*choose(2*l3,l4)
                      ig0=m-2*l1-l2+2*l3-l4
                      ig1=ig0+1
                      ia0=m-2*l1-l4+1
                      ia1=ia0-1
                      ivn=mod(l4,2)+1
                      idd=3-ivn
                      do 3222 j=2,jx
                        tam1(j)=tam1(j)-cons42*gman(ig0,ia0,j)*
     *                    ( gamman(j,1)*tamt1(ivn,j,ia1,ig0)+
     *                    alphan(j,-1)*tamt1(idd,j,ia0,ig1))
                        tam2(j)=tam2(j)-cons41*gman(ig0,ia0,j)*
     *                    ( gamman(j,1)*tamt1(ivn,j,ia1,ig0)+
     *                    alphan(j,-1)*tamt1(idd,j,ia0,ig1))
                        tam4(j)=tam4(j)-cons42*
     *                    (dgawv(ig1,ia0,j)*tamt1(ivn,j,ia1,ig0)+
     *                    dgawv(ig0,ia1,j)*tamt1(idd,j,ia0,ig1))
                        tam5(j)=tam5(j)-cons42*
     *                    (dgawvv(ig1,ia0,j)*tamt1(ivn,j,ia1,ig0)+
     *                    dgawvv(ig0,ia1,j)*tamt1(idd,j,ia0,ig1))
 3222                 continue
 3221               continue ! l4 = 0 : 2*l3
 322              continue ! l3
                else
c.......................................................................:
c     Case 2:  l2 even
c.......................................................................:
                  cons21=cons1*(fctrl(l2+2)*.5**(l2+2))/
     *              (dfloat(l2+2)*fctrl((l2+2)/2)**2)
                  cons22=cons21*(l2+3)/dfloat(2*(l2+1))
                  cons23=cons1*(l2+2)*fctrl(l2/2)**2*2.**l2/
     *              (dfloat(l2+1)*fctrl(l2+1))
c..................................................................
c     Psi1:  Sigma3(j+2,m-2k-j,m-2k+1)=>Sigma2(0,m-2k-j,m-2k+1) term
c     Psi2:  Sigma3(j+2,m-2k-j,m-2k+1)=>Sigma2(0,m-2k-j,m-2k+1) term
c..................................................................
                  ig0=m-2*l1-l2
                  ia0=m-2*l1+1
                  do 325 j=2,jx
                    tam1(j)=tam1(j)-cons22*gman(ig0,ia0,j)*
     *                (asnha(j)*tamt1(2,j,ia0,ig0)+
     *                tamt2(1,j,ia0,ig0))
                    tam2(j)=tam2(j)-cons21*gman(ig0,ia0,j)*
     *                (asnha(j)*tamt1(2,j,ia0,ig0)+
     *                tamt2(1,j,ia0,ig0))
                    tam4(j)=tam4(j)-cons22*
     *                (dgaswv(ig0,ia0,j)*tamt1(2,j,ia0,ig0)+
     *                dgawv(ig0,ia0,j)*tamt2(1,j,ia0,ig0))
                    tam5(j)=tam5(j)-cons22*
     *                (dgaswvv(ig0,ia0,j)*tamt1(2,j,ia0,ig0)+
     *                dgawvv(ig0,ia0,j)*tamt2(1,j,ia0,ig0))
 325              continue
                  do 327 l3=0,l2/2+1
c..................................................................
c     Psi3:  Sigma3(j+1,m-2k-j,m-2k+1)
c..................................................................
                    if (l3.le.l2/2) then
                      cons31=cons23*fctrl(2*l3)*.5**(2*l3)/
     !                  fctrl(l3)**2
                      do 323 l4=0,2*l3
                        cons41=cons31*choose(2*l3,l4)
                        ig0=m-2*l1-l2+2*l3-l4
                        ig1=ig0+1
                        ia0=m-2*l1-l4+1
                        ia1=ia0-1
                        ivn=mod(l4,2) + 1
                        idd=3-ivn
                        do 3231 j=2,jx
                          tam3(j)=tam3(j)+cons41*gman(ig0,ia0,j)*
     *                      ( gamman(j,1)*tamt1(ivn,j,ia1,ig0)+
     *                      alphan(j,-1)*tamt1(idd,j,ia0,ig1))
                          tam6(j)=tam6(j)+cons41*
     *                      (dgawv(ig1,ia0,j)*tamt1(ivn,j,ia1,ig0)+
     *                      dgawv(ig0,ia1,j)*tamt1(idd,j,ia0,ig1))
 3231                   continue
323                  continue ! l4= 0 : 2*l3
                    endif
c..................................................................
c     Psi1:  Sigma3(j+2,m-2k-j,m-2k+1)=>Sigma1(j+2,m-2k-j,m-2k+1) term
c     Psi2:  Sigma3(j+2,m-2k-j,m-2k+1)=>Sigma1(j+2,m-2k-j,m-2k+1) term
c..................................................................
                    if(l3.gt.0) then
                      cons31=cons21*fctrl(l3)*fctrl(l3-1)*2.**(2*l3-1)
     *                  /fctrl(2*l3)
                      cons32=cons22*fctrl(l3)*fctrl(l3-1)*2.**(2*l3-1)
     *                  /fctrl(2*l3)
                      do 326 l4=0,2*l3-1
                        cons41=cons31*choose(2*l3-1,l4)
                        cons42=cons32*choose(2*l3-1,l4)
                        ig0=m-2*l1-l2+2*l3-1-l4
                        ig1=ig0+1
                        ia0=m-2*l1-l4+1
                        ia1=ia0-1
                        ivn=mod(l4,2) + 1
                        idd=3-ivn
                        do 324 j=2,jx
                          tam1(j)=tam1(j)-cons42*gman(ig0,ia0,j)*
     *                      ( gamman(j,1)*tamt1(ivn,j,ia1,ig0)+
     *                      alphan(j,-1)*tamt1(idd,j,ia0,ig1))
                          tam2(j)=tam2(j)-cons41*gman(ig0,ia0,j)*
     *                      ( gamman(j,1)*tamt1(ivn,j,ia1,ig0)+
     *                      alphan(j,-1)*tamt1(idd,j,ia0,ig1))
                          tam4(j)=tam4(j)-cons42*
     *                      (dgawv(ig1,ia0,j)*tamt1(ivn,j,ia1,ig0)+
     *                      dgawv(ig0,ia1,j)*tamt1(idd,j,ia0,ig1))
                          tam5(j)=tam5(j)-cons42*
     *                      (dgawvv(ig1,ia0,j)*tamt1(ivn,j,ia1,ig0)+
     *                      dgawvv(ig0,ia1,j)*tamt1(idd,j,ia0,ig1))
 324                    continue
!!!       write(*,'(a,7i3)') 'm,l1,l2,l3,l4,ia0-1,ig0=', 
!!!     ~ m,l1,l2,l3,l4,ia1,ig0
 326                  continue ! l4= 0 : 2*l3-1
                    endif    ! if(l3>0)
 327              continue ! l3= 0 : l2/2+1
                endif    ! l2 odd/even
 331          continue ! l2= 0 : m-2*l1
!!!         write(*,'(a,2i3,6e10.1)') 'tam1-6, m,l1=', m,l1,
!!!     +tam1(3),tam2(3),tam3(3),tam4(3),tam5(3),tam6(3)
!!!         write(*,'(a,i3,2e10.1)') 'tamt2, m=', m,
!!!     ~   tamt2(1,3,ia1,ig0),tamt2(2,3,ia1,ig0)
 330        continue ! l1= 0 : m/2
!!!        pause

c..................................................................
c     Calculate the Terms in the Local Coefficients
c..................................................................
            anr1=4.*pi*cnorm
            anr2=anr1/cnorm2
            do 350 j=2,jx
              tam7(j)=-anr2*gamman(j,1)*xsq(j)*tam6(j)
              tam8(j)=anr1*xsq(j)*
     *          (gamman(j,3)*tam5(j)+
     *          gamman(j,1)*x(j)*tam4(j)*cnorm2i+
     *          .5*gamman(j,1)*tam2(j)*cnorm2i)
              tam9(j)=anr1*gamman(j,1)*(tam4(j)-xi(j)*tam1(j))
              tam10(j)=-anr2*gamman(j,-1)*tam3(j)
              tam11(j)=anr1*(gman(1,1,j)*tam4(j)/cnorm+
     *          .5*gamman(j,-1)*tam2(j)*cnorm2i)
              tam12(j)=anr1*gman(-1,2,j)*tam1(j)/cnorm2
 350        continue
       
c..................................................................
            do 380 iii=1,imax(ileff,lr_)
              do 370 ii=0,1
                if (madd.eq.2 .and. ii.eq.0) goto 370
                i=iii*ii-(iy+1-iii)*(ii-1)
                do 360 j=2,jx
                 ca(i,j)=ca(i,j)+ss(i,ileff,m,lr_)*tam7(j)*gamefac(j,k) !YuP[2019-07-26] k index added
                 cb(i,j)=cb(i,j)+ss(i,ileff,m,lr_)*tam8(j)*gamefac(j,k) !YuP[2019-07-26] k index added
                 cc(i,j)=cc(i,j)+ssy(i,ileff,m,lr_)*tam9(j)*gamefac(j,k) !YuP[2019-07-26] k index added
                  cd(i,j)=cd(i,j)+sinz(i,ileff,lr_)*
     *              ssy(i,ileff,m,lr_)*tam10(j)*gamefac(j,k) !YuP[2019-07-26] k index added
                  ce(i,j)=ce(i,j)+sinz(i,ileff,lr_)*
     *              ssy(i,ileff,m,lr_)*tam9(j)*gamefac(j,k) !YuP[2019-07-26] k index added
                  cf(i,j)=cf(i,j)+sinz(i,ileff,lr_)*
     *              (ss(i,ileff,m,lr_)*tam11(j)+
     +              ssyy(i,ileff,m,lr_)*tam12(j))*gamefac(j,k) !YuP[2019-07-26] k index added
 360            continue
 370          continue
 380        continue
c...................................................................
       
!!!          pause
 
 400      continue ! m (Legendre)

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
            call dscal(iyjx,anr1,ca(1,1),1) !-YuP: size of ca..cf: iy*jx
            call dscal(iyjx,anr1,cb(1,1),1)
            call dscal(iyjx,anr1,cc(1,1),1)
            call dscal(iyjx,anr1,cd(1,1),1)
            call dscal(iyjx,anr1,ce(1,1),1)
            call dscal(iyjx,anr1,cf(1,1),1)
            call dscal(iyjx,satiom(kk,k),ca(1,1),1)
            call dscal(iyjx,satiom(kk,k),cd(1,1),1)
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
cBH091031                !i.e, not bounce pt interval
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
     &              *cc(ii,j)
                  cdl(ii,j,kk,l_)=cdl(ii,j,kk,l_)+ax*ay*cd(ii,j)
                  cel(ii,j,kk,l_)=cel(ii,j,kk,l_)+ax*ay*ce(ii,j)
                  cfl(ii,j,kk,l_)=cfl(ii,j,kk,l_)+ax*az*cf(ii,j)
 460            continue
 470            continue
 480          continue

            else
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

            call dscal(iyjx,one/anr1,ca(1,1),1) !-YuP: size of ca..cf: iy*jx
            call dscal(iyjx,one/anr1,cb(1,1),1)
            call dscal(iyjx,one/anr1,cc(1,1),1)
            call dscal(iyjx,one/anr1,cd(1,1),1)
            call dscal(iyjx,one/anr1,ce(1,1),1)
            call dscal(iyjx,one/anr1,cf(1,1),1)
            call dscal(iyjx,one/satiom(kk,k),ca(1,1),1)
            call dscal(iyjx,one/satiom(kk,k),cd(1,1),1)
 490      continue

c     end of loop over gen. species k
 500    continue

c     end of loop over orbit l
 600  continue
 700  continue
c..................................................................
c     define needed coefficients at pass/trapped boundary
c..................................................................
      if (cqlpmod .ne. "enabled") then
        do 2001 k=1,ngen
          do 2002 j=1,jx
            cal(itl,j,k,l_)=0.25*vptb(itl,lr_) * ( cal(itl-1,j,k,l_)/
     /        vptb(itl-1,lr_)+2.*cal(itl+1,j,k,l_)/vptb(itl+1,lr_)
     +        +cal(itu+1,j,k,l_)/vptb(itu+1,lr_) )
            cbl(itl,j,k,l_)=0.25*vptb(itl,lr_) * ( cbl(itl-1,j,k,l_)/
     /        vptb(itl-1,lr_)+2.*cbl(itl+1,j,k,l_)/vptb(itl+1,lr_)
     +        +cbl(itu+1,j,k,l_)/vptb(itu+1,lr_) )
            cal(itu,j,k,l_)=cal(itl,j,k,l_)
            cbl(itu,j,k,l_)=cbl(itl,j,k,l_)
 2002     continue
 2001   continue
      endif
c..................................................................
      if (madd .eq. 2) call cfpsymt
      
      do 2000 k=1,ngen
        call dscal(iyjx,one/tnorm(k),cal(1,1,k,l_),1)
        call dscal(iyjx,one/tnorm(k),cbl(1,1,k,l_),1)
        call dscal(iyjx,one/tnorm(k),ccl(1,1,k,l_),1)
        call dscal(iyjx,one/tnorm(k),cdl(1,1,k,l_),1)
        call dscal(iyjx,one/tnorm(k),cel(1,1,k,l_),1)
        call dscal(iyjx,one/tnorm(k),cfl(1,1,k,l_),1)
        call dscal(iyjx,one/tnorm(k),eal(1,1,k,1,l_),1)
        call dscal(iyjx,one/tnorm(k),ebl(1,1,k,1,l_),1)
        call dscal(iyjx,one/tnorm(k),eal(1,1,k,2,l_),1)
        call dscal(iyjx,one/tnorm(k),ebl(1,1,k,2,l_),1)
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
      return
      end
