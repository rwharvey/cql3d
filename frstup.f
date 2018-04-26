c
c
      subroutine frstup(mf,mfm1,mi,mj,nion,potsid,codeid,rin,
     1  rmax,zax,zminn,zmaxx,zni,zne,zte,zti,zzi,frpsi,psivol,
     1  xxx,yyy,nprim,nimp,zeffctv,zshift1)
      implicit integer (i-n), real*8 (a-h,o-z)

      include 'param.h'
      include 'comm.h'

      character*8 codeid 
      dimension psivol(*),potsid(*),frpsi(nnra,*),
     1  zne(kz,*),zni(kz,*),zte(*),zzi(kz,*),xxx(*),yyy(*),
     1  zti(*),zeffctv(*)

c..................................................................
c     This routine should be called just before the call to FREYA.
c     It defines a number of time (or iteration) dependent input
c     variables required in the call to FREYA.
c..................................................................


c..................................................................
c     Set the frpsi (epsi) array as required by FREYA
c     Maximum will be at limiter not the axis.
c..................................................................

      do 10 i=1,nnr
        xxx(i)=er(i)
        do 11 j=1,nnz
          frpsi(i,j)=-epsi(i,j)
 11     continue
 10   continue
      do 15 j=1,nnz
        yyy(j)=ez(j)
 15   continue
c.....................................................................
c     eqpsi(1:nconteqn) is set up in cql3d routine eqrhopsi.f.
c     It is used to give the flux zones for freya returned in the
c     izone argument of subroutine inject (now inject_old and inject1),
c     called in subroutine freya.
c     Change the sign of (positive) eqpsi array to get it in ascending 
c     order for splines.  (Changed back at end of subroutine).
c.....................................................................

      do 2 l=1,nconteqn
        eqpsi(l)=-eqpsi(l)
 2    continue
c..................................................................
c     The number of flux zones (mfm1) used in FREYA
c..................................................................

      mfm1=nconteqn-1
      mf=mfm1+1

      write(*,*)
      write(*,*)'Number of flux zones (mfm1) used in FREYA=',mfm1
      write(*,*)'Adjust mfm1=nconteqn-1 throught NL input of nconteqn'
      write(*,*)

c..................................................................
c     Interpolate densities, temperatures etc over to the 
c     (high resolution?) mesh utilized by FREYA.
c
c..................................................................

      psivol(1)=.5*(eqvol(2)+eqvol(1))
      psivol(mfm1)=eqvol(mf)-.5*(eqvol(mfm1)+eqvol(mfm1-1))
      do 20 ll=2,mfm1-1
        psivol(ll)=.5*(eqvol(ll+1)-eqvol(ll-1))
 20   continue

c..................................................................
c     Electron density, energy
c..................................................................

      do 30 ll=1,lrzmax
        tr1(ll)=reden(kelec,ll)
        tr2(ll)=energy(kelec,ll)*2./3.
 30   continue

c..................................................................
c     Call spline interpolation routine to get values on FREYA 
c     eqspi mesh
c     (which is same as intemediary radial mesh used in CQL3D proper).
c..................................................................

      call frsplft(lrzmax,equilpsp(1),tr1(1),mfm1,eqpsi,zne)
      call frsplft(lrzmax,equilpsp(1),tr2(1),mfm1,eqpsi,zte)

c..................................................................
c     Average ion temperature..
c..................................................................

      kimp_=ntotal-nimp+1
      call bcast(tr1(1),zero,lrzmax)
      do 50 ll=1,lrzmax
        tp=0.
        dn=0.
        kk=0
        do 40 k=1,ntotal
          if (k.eq.kelecg .or. k.eq.kelecm) go to 40
          kk=kk+1
          if (kk.gt.nprim .and. k.lt.kimp_) then
            kk=kk-1
            go to 40
          endif
          tp=tp+2./3.*energy(k,ll)*reden(k,ll)
          dn=dn+reden(k,ll)
 40     continue
        tr1(ll)=tp/dn
 50   continue
      call frsplft(lrzmax,equilpsp(1),tr1(1),mfm1,eqpsi,zti(1))

c..................................................................
c     Ion densities..
c..................................................................
c
      kk=0
      do 70 k=1,ntotal
        call bcast(tr1(1),zero,lrzmax)
        if (k.eq.kelecm .or. k.eq.kelecg) go to 70
        kk=kk+1
c        if (kk.gt.nprim .and. k.lt.kimp_) then
        if (kk.gt.nprim+nimp+1 .and. k.lt.kimp_) then
          kk=kk-1
          go to 70
        endif
        do 60 ll=1,lrzmax
          tr1(ll)=reden(k,ll)
c         write(*,*) k, ll, tr1(ll), ' ion density'
 60     continue
        call frsplft(lrzmax,equilpsp(1),tr1(1),mfm1,eqpsi,zni(1,kk))
 70   continue
c      stop
c..................................................................
c     Average charge state..
c..................................................................
      kk=0
      do 90 k=1,ntotal
        if (k.eq.kelecm .or. k.eq.kelecg) go to 90
        kk=kk+1
c        if (kk.gt.nprim .and. k.lt.kimp_) then
        if (kk.gt.nprim+nimp+1 .and. k.lt.kimp_) then
          kk=kk-1
          go to 90
        endif
        do 80 ll=1,lrzmax
          tr1(ll)=bnumb(k)
c         write(*,*) k, ll, tr1(ll), ' ion charge'
 80     continue
        call frsplft(lrzmax,equilpsp(1),tr1(1),mfm1,eqpsi,zzi(1,kk))
 90   continue
c      stop
c..................................................................
c     Compute zeff(lr_) on Freya mesh.
c..................................................................

      call frsplft(lrzmax,equilpsp(1),zeff(1),mfm1,eqpsi,zeffctv(1))
c
       if(nprim.eq.1 .and. nimp.eq.1) then
         write(*,*) 'copying Zeff into zzi...'
         do k=1,mfm1
           zzi(k,nion+2)=zeff(k)
c           write(*,*) 'zzi(k,nion+2) = ',zzi(k,nion+2)
         enddo
       endif
c
c..................................................................
c     Change eqpsi back to convention at call to this subroutine.
c..................................................................

      do 101 l=1,nconteqn
        eqpsi(l)=-eqpsi(l)
 101  continue

c..................................................................
c     For running NFREYA with CQL3D, we assume iborb=0; so
c     it is necessary to set only the endpoints of potsid
c..................................................................

      potsid(1)=-psimag
      potsid(mf)=-psilim
      elong=0
      codeid="twodee"
      rin=rmincon
      write(*,*)'frsetup: rmincon,rmaxcon=',rmincon,rmaxcon
      rmax=rmaxcon
      zminn=zmincon
      zmaxx=zmaxcon
      zax=0.
      mi=nnr
      mj=nnz
      zshift1=zshift

      return
      end
