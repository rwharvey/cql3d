c
c
      subroutine tdtscout
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     This routine creates a file 'tscout'; data
c     is read in from the file created by the TSC code. The file
c     created here has the power and the current as a function of
c     rho_, the radial coordinate given by TSC.
c..................................................................

      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

      dimension powtsc(nrada),currtsc(nrada)

CMPIINSERT_IF_RANK_NE_0_RETURN

      if (eqsource.ne."tsc") return

c..................................................................
c     Create output file.
c..................................................................

      open(unit=12,file='tscout',delim='apostrophe',status='unknown')

c..................................................................
c     interpolate to TSC mesh
c..................................................................

      do 10 l=1,lrzmax
        tr(l)=sorpw_rf(kelecg,l)
 10   continue
      call tdinterp("free","free",rya(1),tr(1),lrzmax,rho_,powtsc,
     +  npsitm)
      call tdinterp("free","free",rya(1),currtpz(1),lrzmax,rho_,currtsc
     1  ,npsitm)

c..................................................................
c     Now write to disk.
c     powtsc is the power in Watts/cc
c     currtsc is the current in Amps/cm**2
c..................................................................

      write(12,100) (powtsc(l),l=1,npsitm)
      write(12,100) (currtsc(l),l=1,npsitm)
 100  format(5e16.6)
      close(unit=12)
      return
      end
