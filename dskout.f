c
c
      subroutine dskout(ll)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     At the end of the run, this routine at the option of the user 
c     writes out to disk file 'idskf' the namelist input deck 
c     and various computed quantities, including the distn functions.
c     Also, write out to disk file 'idskrf' rf diffusion coefficients,
c     and related quantities.
c..................................................................

      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

      character*80 line
      
CMPIINSERT_IF_RANK_NE_0_RETURN

      if (lrzmax.le.1) then
        if((idskf.eq. "disabled".and.idskrf.eq."disabled")
     1    . or. n.ne.nstop+1)  return
      else
        if((idskf.eq. "disabled".and.idskrf.eq."disabled")
     1    . or. n.ne.nstop)  return
      endif
      if(ll.gt.1)  go to 4
      if(idskf.ne."disabled")
     1  open(unit=4,file=idskf,delim='apostrophe',status='unknown')
      if(idskrf.ne."disabled")
     1  open(unit=5,file=idskrf,delim='apostrophe',status='unknown')
ccc      close(unit=2) ! YuP: Why here?
      ilen=0

c..................................................................
c     The input namelist file is transcribed onto the beginning
c     of file idskf and/or idskrf
c..................................................................

      open(unit=2,file='cqlinput',delim='apostrophe',status='old')
 1    read(2,1003) line
      if (line(1:3).eq."end") go to 3
      if(idskf.ne."disabled") write(4,1003) line
      if(idskrf.ne."disabled") write(5,1003) line
      go to 1
 3    close(unit=2)
      if(idskf.ne."disabled") 
     +     write(4,1006) '***Begin computed output***'
      if(idskrf.ne."disabled") 
     +     write(5,1006) '***Begin computed output***'
 4    continue 

c..................................................................
c     In the following disk write to file named idskf:
c     This subroutine is called to write data for each FP'd
c          flux surface.
c     ll=  FP flux surface number 
c          (ll=1:lrors, lrors.le.lrzmax, see cqlinput_help))
c          (lrindx(ll) gives flux surface number on the full
c                     radial mesh.(lrindx(ll)=ll if lrors=lrzmax,
c                     and using cql3d mode(cqlpmod="disabled"))).
c     iy,jx= dimensions in theta and u(momentum/mass)
c           (In the case where iy varies with ll, iy(1) will be greatest.)
c     lrors= number of flux surfaces FP'd.
c     lrzmax= number of flux surfaces, including any not FP'd.
c     x = momentum-per-mass(nomalized to maximum 1.)
c         at  each flux surface lrindx(ll)
c     y = theta(radians) mesh at  each flux surface lrindx(ll)
c     rovera= normalized radius (ll)
c             (rho, see Hinton and Haseltine for non-circ).
c             (generally ~sqrt(tor. flux), other coords available.)
c     elecfld = toroidal electric field (volts/cm)
c     bthr(lr_) - the poloidal magnetic field at theta-poloidal = pi/2.
c     btoru(lr_) - the toroidal magnetic field at the same position. 
c     bthr0(lr_), btor0(lr_), are the poloidal and  toroidal
c         magnetic fields at the outer midplane of the flux surface. 
c     reden= electron density at minimum B point on flux surface.
c     temp= initial electron temperature (keV)
c     radmin= plasma minor radius (cms).
c     vnorm= normalization momentum-per-mass (maximum on grid) (cm/sec)
c     vmaxdvt= vnorm/(temp/mass)**0.5
c     eovedd= electric  field, normalized to Driecer field 
c             (calc'd in sub restvty).
c
c     distribution function normalised so that
c         integral( (dx)**3 f) = density at minumum B point.
c
c..................................................................

      if(idskf.ne."disabled") then
        write(4,1004)  ll, iy,jx,lrors,lrzmax,ngen
        write(4,1004)  itl,itu
        write(4,1005)  (x(j),j=1,jx)
        write(4,1005)  (y(i,ll),i=1,iy)
        do 1000 k=1,ngen
          write(4,1005)  bnumb(k),fmass(k)
          write(4,1005)  rovera(lrindx(ll)),elecfld(lrindx(ll)),
     +                   bthr(lrindx(ll)),btoru(lrindx(ll))
          write(4,1005)  bthr0(lrindx(ll)),btor0(lrindx(ll)),
     +                   reden(k,lrindx(ll)),temp(k,lrindx(ll))
          vmaxdvt=vnorm/(4.19e7*sqrt(2.*1000.*temp(k,lrindx(ll))))
          write(4,1005)  radmin,vnorm,vmaxdvt,eovedd
          write(4,1005)  ((f(i,j,k,ll),i=1,iy),j=1,jx)
 1000   continue
      endif


c..................................................................
c     In the following disk write:
c     vptb is lambda=u_parallel_o * tau_bounce.
c     (tau_bounce= portion of bounce time from the outer equatorial
c     plane to the bounce point or the inner equatorial plane).
c     temp1 is the bounce averaged uu-diffusion coefficient (cgs units).
c..................................................................
      if(idskrf.ne."disabled") then
        write(5,1004)  ll,iy,jx,lrors,lrzmax
        write(5,1005)  (y(i,ll),i=1,iy)
        write(5,1005)  (x(j),j=1,jx)
        write(5,1005)  rovera(lrindx(ll)),vnorm
        write(5,1005)  (vptb(i,lrindx(ll)),i=1,iy)
        vn2=vnorm*vnorm
        do 10 k=1,mrfn
          vmaxdvt=vnorm/(4.19e7*sqrt(2.*1000.*temp(k,lrindx(ll))))
          write(5,1005) reden(k,lrindx(ll)),temp(k,lrindx(ll)),vmaxdvt
          do 11 i=1,iy
 11       temp1(i,1)=0.0
          do 12 j=2,jx
            x2=x(j)**2
            do 13 i=1,iy
 13         temp1(i,j)=urfb(i,j,indxlr(lrindx(ll)),k)*vn2/
     /          (vptb(i,lrindx(ll))*x2)
 12       continue
          write(5,1005) ((temp1(i,j),i=1,iy),j=1,jx)
 10     continue
      endif
 1003 format(a80)
 1004 format(16i5)
 1005 format(5e16.8)
 1006 format(a27)
      if(idskf.ne."disabled".and.ll.eq.lrors)  close(unit=4)
      if(idskrf.ne."disabled".and.ll.eq.lrors)  close(unit=5)
      return
      end
