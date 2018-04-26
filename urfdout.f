c
c
      subroutine urfdout(ll)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     This subroutine writes a disk file which is read by the Ray Tracin
c     code (Brambilla).
c..................................................................

      include 'param.h'
      include 'comm.h'
CMPIINSERT_INCLUDE

      character*1 line(132)
c jakub urban 110708: commented out for g95 compiler
c the above blanket save statement should do the job
c      save ifirst
      data ifirst /0/

CMPIINSERT_IF_RANK_NE_0_RETURN

c..................................................................
c     Create or open a disk file.
c..................................................................

      isizeg=iy*jx*(lrors+1)*ngen+20000
      if(ll.eq.1.and.ifirst.eq.0)  then
        open(unit=31,file=idskf,delim='apostrophe',status='new')
      else if (ll.eq.1) then
        open(unit=31,file=idskf,delim='apostrophe',status='old')
      endif
      ifirst=1

c..................................................................
c     Write out ll (flux surf label) iy jx lrors (# of fl. surf)
c     Write out y (theta mesh) and x (momentum mesh).
c     Write out rovera(lr_) (r/radmin), vnorm (x mesh normal. factor).
c..................................................................


      write(31,555)  ll, iy,jx,lrors
      write(31,560)  (y(i,l_),i=1,iy)
      write(31,560)  (x(j),j=1,jx)
      do 1000 k=1,ngen
        write(31,560) rovera(lr_),vnorm,elecfld(lr_),eovedd,reden(k,lr_)
        vmaxdvt=vnorm/(4.19e7*sqrt(2.*1000.*temp(k,lr_)))
        write(31,560)  temp(k,lr_),hnis(k,lr_),vmaxdvt
        write(31,560)  ((f(i,j,k,l_),i=1,iy),j=1,jx)
 1000 continue
 555  format(5i16)
 560  format(5e16.8)
      if(ll.eq.lrz)  close(unit=31)

      return
      end
