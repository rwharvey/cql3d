c
c
      subroutine urfwrite_(krf,i_)
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'param.h'
      include 'comm.h'
c..................................................................
c     Write total number of rays.
c     Write first harmonic number (used in damping).
c     Write in the wave antenna frequency.
c..................................................................

      write(i_,1) nray(krf),nharm(krf),freqcy(krf)

c..................................................................
c     Write in ray data, ray by ray
c     nrayelt(iray,k)=number of ray elements for ray iray,k
c     The next six integers in the write are the final states of the mode,
c     number of refl., integ.steps, number of electron and ion power calcs
c     and the restart status of the ray.
c     sxxrt(iray,k) - final value of integ. variable
c     skpsi(iray,k) - n_psi
c     istarts(iray,k) - n_th
c     ws - poloidal plane projection of distance along ray.
c     seikon - eikon
c     spsi - rho
c     wr - R
c     wphi - toroidal angle phi
c     wz - Z
c     wnpar - n_parallel
c     wnper - n_perp
c     delpwr - power flowing in the ray channel.
c     sdpwr - dpwri (not presently used for calcs. in the code)
c     wdnpar - delta n_parallel (n_parallel width locally)
c     cweyde - -i*cweyde is a real number = (E_y/E)
c     cwezde - the polarization E_z/E
c     fluxn - the normalized flux (electric field=1)
c     sbtot - total B field (gauss)
c     sene - n_e(cc)
c     salphac - perp. coll. absorption coefficient (/cm).
c     salphal - additional linear perp. absorption coeff. (on ions) (/cm).
c     All units are cgs.
c..................................................................

c     Truncating difficult (large neg) exponents, for formattted o/p:
      call urfwr0(sxxrt(1,krf),nray(krf),1,nrayn)
      call urfwr0(skpsi(1,krf),nray(krf),1,nrayn)
      call urfwr0(skth(1,krf),nray(krf),1,nrayn)
      call urfwr0(skphi(1,krf),nray(krf),1,nrayn)

      call urfwr0(ws(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(seikon(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(spsi(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(wr(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(wphi(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(wz(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(wnpar(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(wnper(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(delpwr(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(sdpwr(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(wdnpar(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0c(cwexde(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0c(cweyde(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0c(cwezde(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(fluxn(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(sbtot(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(sene(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(salphac(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)
      call urfwr0(salphal(1,1,krf),nrayelt(1,krf),nray(krf),nrayelts)


c.....................................................................
c     Adjust data to standard external format specified in urfread_.f. 
c     (See urfread_.f. Perform inverse transformation)
c.....................................................................
c
      do iray=1,nray(krf)

      do  is=1,nrayelt(iray,krf)
         wdnpar(is,iray,krf)=abs(wdnpar(is,iray,krf))/wdscale(krf)
      enddo

c     Change sign of sbtot and wnpar, if bsign.lt.0.
c     (see cqlinput_help)
      if (eqsource.eq."eqdsk" .and. bsign.lt.0.) then
         sbsign=sign(one,bsign)
         do is=1,nrayelt(iray,krf)
            sbtot(is,iray,krf)=sbtot(is,iray,krf)/bsign1(krf)
            wnpar(is,iray,krf)=sbsign*wnpar(is,iray,krf)
            cwexde(is,iray,krf)=sbsign*cwexde(is,iray,krf)
            cweyde(is,iray,krf)=sbsign*cweyde(is,iray,krf)
         enddo
      endif
      
c     fluxn is renormalized:
      do is=1,nrayelt(iray,krf)
         fluxn(is,iray,krf)=fluxn(is,iray,krf)*(8.*pi)/clight
      enddo
      
c     shift z-position of ray data, if eqdsk has been
c     vertically shifted in subroutine equilib:
      if (zshift.ne.zero) then
         do is=1,nrayelt(iray,krf)
            wz(is,iray,krf)=wz(is,iray,krf)+zshift
         enddo
      endif

      enddo  ! iray

c..................................................................


      do 10 iray=1,nray(krf)
        write(i_,3) nrayelt(iray,krf),jslofas(iray,krf),
     +    nurefls(iray,krf),keiks(iray,krf),jpes(iray,krf),
     2    jpis(iray,krf),istarts(iray,krf),iprmt5(iray,krf),
     +    jhlfs(iray,krf)

        write(i_,2) sxxrt(iray,krf),skpsi(iray,krf),skth(iray,krf),
     1    skphi(iray,krf)

        write(i_,2) (ws(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (seikon(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (spsi(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (wr(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (wphi(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (wz(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (wnpar(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (wnper(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (delpwr(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (sdpwr(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (wdnpar(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (cwexde(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (cweyde(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (cwezde(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (fluxn(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (sbtot(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (sene(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (salphac(is,iray,krf),is=1,nrayelt(iray,krf))
        write(i_,2) (salphal(is,iray,krf),is=1,nrayelt(iray,krf))
 10   continue
 1    format(2i5,1pe16.7)
 2    format(5(1pe16.7))
 3    format(9i5)

c..................................................................
c     Re-Adjust data to internal format, 
c     in case this is not last step. 
c     (See urfread_.f)
c..................................................................
c
      do iray=1,nray(krf)

      do  is=1,nrayelt(iray,krf)
         wdnpar(is,iray,krf)=wdscale(krf)*abs(wdnpar(is,iray,krf))
      enddo

c     Change sign of sbtot and wnpar, if bsign.lt.0.
c     (see cqlinput_help)
      if (eqsource.eq."eqdsk" .and. bsign.lt.0.) then
         sbsign=sign(one,bsign)
         do is=1,nrayelt(iray,krf)
            sbtot(is,iray,krf)=bsign1(krf)*sbtot(is,iray,krf)
            wnpar(is,iray,krf)=sbsign*wnpar(is,iray,krf)
            cwexde(is,iray,krf)=sbsign*cwexde(is,iray,krf)
            cweyde(is,iray,krf)=sbsign*cweyde(is,iray,krf)
         enddo
      endif
      
c     fluxn is renormalized to be as in Stix or Bekefi:
      do is=1,nrayelt(iray,krf)
         fluxn(is,iray,krf)=fluxn(is,iray,krf)*clight/(8.*pi)
      enddo
      
c     shift z-position of ray data, if eqdsk has been
c     vertically shifted in subroutine equilib:
      if (zshift.ne.zero) then
         do is=1,nrayelt(iray,krf)
            wz(is,iray,krf)=wz(is,iray,krf)-zshift
         enddo
      endif

      enddo  !  iray

c..................................................................


      close(unit=i_)
      return
      end
