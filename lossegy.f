c
c
      subroutine lossegy
      implicit integer (i-n), real*8 (a-h,o-z)
      save
c
cmnt  this routine computes fp coefficients corresponding to the
cmnt  bounce average of the phenomonological energy loss term,
cmnt  or a Bremsstrahlung energy loss term:

c     Phenomonological loss:
c     df/dt = u**-2*d/du [1./tau(u,theta,lr_)*u**3*f/2.]
c     where   tau is energy loss time specified by
c     tau=rfact*tauegy*gamma**gamegy, u.le.vte
c     rfact*tauegy*gamma**gamegy/
c     (abs(u_par,lr_)/vte)**paregy*(u_per/vte)**pe
c     *(u/vte)**pegy),  u.gt.vte
c     rfact is a radially dependent modifier, as below.
c
c     Bremsstrahlung loss due to electron collisions with (only)
c     the Maxwellian ion species:
c     df/dt = u**-2*d/du [beta*u**3*f],
c     where beta=reden(kion,lr_)*sigmarad*clight**2/(u/gamma),
c           sigmarad=as given in Physics Vade Mecum, H.L.Anderson, Ed.,
c                    AIP, Sect. 16.07.F.
c     Sigmarad may be enhanced beyond energy brfacgm3*16.43/Z**(1/6),
c     using brfac,brfac1 (see cqlinput_help and coding below).
c     This is for purposes of stopping flow of runaway electrons
c     off the grid.



      include 'param.h'
      include 'comm.h'
      data gam1/2.30097589089281/,   gam2/16.4312912649856/

      do 200 k=1,ngen

         call bcast(egylosa(0,0,k,indxlr_),zero,(iy+2)*(jx+2))

c     Phenomenological:

         if(tauegy(k,lr_).le.0.)  go to 130
         coefnt=1./(2.*tauegy(k,lr_))
c     rfact modifier multiplies tau by (1.-rovera(lr_),lr_):
         if(regy(k).eq."enabled") coefnt=coefnt/(1.-rovera(lr_))
c     
c..................................................................
c     Note: the do loop below uses vth(),
c     vth is the thermal velocity =sqrt(T/m) (at t=0 defined in ainpla).
c     But, T==temp(k,lr) can be changed in profiles.f, 
c     in case of iprote (or iproti) equal to "prbola-t" or "spline-t"
c..................................................................
         do 100  j=1,jx
            vel=x(j)*vnorm
            coefnt1=coefnt*xsq(j)*x(j)/gamma(j)**gamegy(k)
            if(vel.le.vth(k,lr_))  then
               do 110  i=1,iy
                  egylosa(i,j,k,indxlr_)=vptb(i,lr_)*coefnt1
 110           continue
            else
               do 120  i=1,iy
                  egylosa(i,j,k,indxlr_)=
     +                 vptb(i,lr_)*coefnt1*(vel/vth(k,lr_))**
     +                 (paregy(k)+peregy(k)+pegy(k))*sincosba(i,k,lr_)
 120           continue
            endif
 100     continue
 130     continue

c     Bremsstrahlung:
         if (k.eq.kelecg .and. bremsrad.eq."enabled") then
            if(brfacgm3.lt.1.) stop 'Check brfacgm3, see cqlinput_help'
            do 150 kk=1,nionm
               kkk=kionm(kk)
               c1=fmass(kelecg)*clite2
               c2=(charge**2/c1)**2/137.*bnumb(kkk)**2
               c3=1./bnumb(kkk)**(1./3.)
               c5=sqrt(c3)
               c4=183.*c3
               c3=137.*c3
               gam3=c5*gam2
               gam4=brfacgm3*gam3
               sigmarad=0.
               do 160 j=1,jx
                  if (gamm1(j).lt.gam1) then
                     sigmarad=16./3.*c2
                  elseif (gamm1(j).gt.gam1.and.
     1                   gamm1(j).le.gam3) then
c990131                     sigmarad=8.*c2*(alog(gamm1(j))-1./6.)
                     sigmarad=8.*c2*(log(gamm1(j))-1./6.)
c990131                     sigmarad=8.*c2*(alog(gamm1(j))-1./6.)
                     sigmarad=8.*c2*(log(gamm1(j))-1./6.)
                  elseif (gamm1(j).gt.gam3.and.
     1                  gamm1(j).le.gam4) then
c990131                     sigmarad=4.*c2*(alog(c4)+1./18.)
                     sigmarad=4.*c2*(log(c4)+1./18.)
                  elseif (gamm1(j).gt.gam4) then
                     gam5=gamm1(j)-gam4
c990131                     sigmarad=4.*c2*(alog(c4)+1./18.)*(1.+brfac*gam5
                     sigmarad=4.*c2*(log(c4)+1./18.)*(1.+brfac*gam5
     1                    **brfac1)
                  endif
                  betau3=reden(kkk,lr_)*sigmarad*clite2*gamma(j)*
     +                 xsq(j)/vnorm
                  do 170 i=1,iy
                     egylosa(i,j,k,indxlr_)=egylosa(i,j,k,indxlr_)+
     +                    vptb(i,lr_)*betau3
 170              continue
 160           continue
 150        continue
         endif

 200  continue

      return
      end
