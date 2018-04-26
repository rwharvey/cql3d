cBH171014: Commenting, since apparently not called, and gfortran 
cBH171014: is having a problem with it.

c$$$c
c$$$c
c$$$c     ONETWO DIVERGENCE
c$$$      subroutine frhexdrv(mb,sgxn,sgxnmi,hxfrac)
c$$$      implicit integer (i-n), real*8 (a-h,o-z)
c$$$c     write statements commented out
c$$$c
c$$$c---hexdrv is a driver for routine hexnb.
c$$$c---hexnb calculates the total reaction rate,rerate(per second), and
c$$$c---the total mean free path rmfp.
c$$$c---hexdrv gets the cross section,sgxn, and inverse of the
c$$$c---pseudocross-section, sgxnmi. hexnb takes excited states of the
c$$$c---neutral beam into account if iexcit=2. if iexcit=1 hexnb should retu
c$$$c---essentially the same cross sections as crsecs does. note
c$$$c---however that the cross section parameterization is different
c$$$c---inhexnb and crsecs.
c$$$c     ONETWO DIVERGENCE
c$$$      include 'param.h'
c$$$      include 'frcomm.h'
c$$$c
c$$$c
c$$$c     ONETWO DIVERGENCE
c$$$      dimension sgxn(kz,ke,kb),sgxnmi(ke,kb),hxfrac(ke,kb)
c$$$c
c$$$c
c$$$c---setup variables for hexnb
c$$$c
c$$$      iexcm1=iexcit-1
c$$$      istart=1
c$$$      if(nouthx.ne.0)nouthx=37
c$$$c     if(iexcm1.eq.0)write(ncrt,10)
c$$$c     if(iexcm1.eq.1)write(ncrt,15)mstate+1
c$$$c     10 format(2x,'entering hexnb,no beam excitation ')
c$$$c     15 format(2x,'entering hexnb, with beam excitation accounted for',
c$$$c     &/,'max principal assumed quantum number is',i5,/,
c$$$c     &'lorentz ionization limit not implemented at this time')
c$$$c
c$$$c---separate primary and impurity ions for easier reading
c$$$c---note:nion=nimp+nprim
c$$$c
c$$$      do 100 j=1,nion
c$$$        if(j.gt.nprim)go to 120
c$$$        atwpm(j)=atw(j)
c$$$        go to 100
c$$$ 120    i=j-nprim
c$$$        atwim(i)=atw(j)
c$$$ 100  continue
c$$$c
c$$$c---get atomic number of impurities
c$$$c
c$$$      do 180 i=1,nimp
c$$$          if(trim(namei(i)).eq.'he' .or. 
c$$$     +       trim(namei(i)).eq.'HE' .or. 
c$$$     +       trim(namei(i)).eq.'He')  iz(i) =  2
c$$$     
c$$$          if(trim(namei(i)).eq.'b' .or. 
c$$$     +       trim(namei(i)).eq.'B' )  iz(i) =  5  ! YuP added [2015]
c$$$          
c$$$          if(trim(namei(i)).eq.'c' .or. 
c$$$     +       trim(namei(i)).eq.'C' )  iz(i) =  6
c$$$          
c$$$          if(trim(namei(i)).eq.'o' .or. 
c$$$     +       trim(namei(i)).eq.'O' )  iz(i) =  8
c$$$          
c$$$          if(trim(namei(i)).eq.'si' .or. 
c$$$     +       trim(namei(i)).eq.'SI' .or. 
c$$$     +       trim(namei(i)).eq.'Si')  iz(i) = 14
c$$$          
c$$$          if(trim(namei(i)).eq.'ar' .or. 
c$$$     +       trim(namei(i)).eq.'AR' .or. 
c$$$     +       trim(namei(i)).eq.'Ar')  iz(i) = 18
c$$$          
c$$$          if(trim(namei(i)).eq.'cr' .or. 
c$$$     +       trim(namei(i)).eq.'CR' .or. 
c$$$     +       trim(namei(i)).eq.'Cr')  iz(i) = 24
c$$$          
c$$$          if(trim(namei(i)).eq.'fe' .or. 
c$$$     +       trim(namei(i)).eq.'FE' .or. 
c$$$     +       trim(namei(i)).eq.'Fe')  iz(i) = 26
c$$$          
c$$$          if(trim(namei(i)).eq.'ni' .or. 
c$$$     +       trim(namei(i)).eq.'NI' .or. 
c$$$     +       trim(namei(i)).eq.'Ni')  iz(i) = 28
c$$$          
c$$$          if(trim(namei(i)).eq.'kr' .or. 
c$$$     +       trim(namei(i)).eq.'KR' .or. 
c$$$     +       trim(namei(i)).eq.'Kr')  iz(i) = 36
c$$$          
c$$$          if(trim(namei(i)).eq.'mo' .or. 
c$$$     +       trim(namei(i)).eq.'MO' .or. 
c$$$     +       trim(namei(i)).eq.'Mo')  iz(i) = 42
c$$$          
c$$$          if(trim(namei(i)).eq.'w' .or. 
c$$$     +       trim(namei(i)).eq.'W' )  iz(i) = 74
c$$$ 180  continue
c$$$c
c$$$c---create diagnostic output file hexnbout
c$$$c---hexnbout contains detailed cross section info
c$$$c
c$$$      if(nouthx.gt.0)
c$$$     &  open(unit=nouthx,file='hexnbout',status='new')
c$$$c
c$$$c---loop over spatial zones,number of beams and beam components
c$$$c
c$$$      do 200 i=1,mfm1
c$$$        teev=1.e3*zte(i)
c$$$        tiev=1.e3*zti(i)
c$$$        do 220 j=1,nion
c$$$          if(j.gt.nprim)go to 230
c$$$          znipm(j)=zni(i,j)
c$$$          go to 220
c$$$ 230      k=j-nprim
c$$$          zniim(k)=zni(i,j)
c$$$ 220    continue
c$$$        do 240 ib=1,mb
c$$$          do 250 j=1,3
c$$$            er0=1.e3*ebkev(ib)/(j*atw(ibion))
c$$$c---no lorentz ionization limit for now
c$$$            bperp=0.0
c$$$            call hexnb(istart,iexcm1,ilorent,mstate,ncont,
c$$$     &        er0,teev,tiev,nprim,atwpm,znipm,nimp,
c$$$     &        iz,atwim,izstrp,zniim,bperp,kdene,kdeni,
c$$$     &        kdenz,ksvi,ksvz,ksve,krad,ngh,ngl,nouthx,
c$$$     &        ncorin,rerate,rmfp,hexfrac,ihxerr)
c$$$c     if(ihxerr.ne.0)write(ncrt,20)ihxerr
c$$$            if(ihxerr.eq.1)then
c$$$c     write(ncrt,21)
c$$$c     write(nout,21)
c$$$              stop 'frhexdrv: problem'
c$$$            endif
c$$$c     20 format(2x,'hexnb error = ',i5)
c$$$c     21 format(2x,'execution terminated - file coronb absent')
c$$$            sgxn(i,j,ib)=1./rmfp
c$$$            if(i.eq.1)hxfrac(j,ib)=hexfrac
c$$$ 250      continue
c$$$ 240    continue
c$$$ 200  continue
c$$$c---  
c$$$c---get inverse of pseudo cross section
c$$$c---  
c$$$      do 400 ib=1,mb
c$$$        if((ib.gt.1).and.(ebkev(ib).eq.ebkev(1)))go to 450
c$$$        do 420 j=1,3
c$$$          sgxnm=sgxn(1,j,ib)
c$$$          do 430 i=2,mfm1
c$$$c990131 430      sqxnm=amax1(sgxnm,sgxn(i,j,ib))
c$$$ 430      sqxnm=max(sgxnm,sgxn(i,j,ib))
c$$$ 420    sgxnmi(j,ib)=1./sgxnm
c$$$        go to 400
c$$$ 450    do 460 j=1,3
c$$$ 460    sgxnmi(j,ib)=sgxnmi(j,1)
c$$$ 400  continue
c$$$      write(*,*) 'sgxnmi',sgxnmi(1,1),sgxnmi(2,1),sgxnmi(3,1)
c$$$      write(*,*) 'sgxn',sgxn(1,1,1),sgxn(2,1,1)
c$$$      write(ncrt,16)
c$$$ 16   format(2x,'leaving hexnb ')
c$$$      return
c$$$      end
