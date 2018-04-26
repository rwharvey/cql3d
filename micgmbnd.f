c
c
      subroutine micgmbnd(xlbmd2)
      implicit integer (i-n), real*8 (a-h,o-z)

c......................................................................
c     This routine computes the pass/trapped boundary volume element.
c     See the user manual for mathematical details (pp. 101 and 130-131 
c      of Killeen, Kerbel, McCoy, Mirin book).
c     It is called from subroutine micgnbnd.
c......................................................................

      save
      include 'param.h'
      include 'comm.h'

      uupls=((y(itl-1,l_)+y(itl,l_))*.5)
      uumns=((y(itl,l_)+y(itl+1,l_))*.5)
      upls=cos(uupls)
      umns=cos(uumns)
      delmu=(upls-umns)*0.005
      hh=zmax(lr_)-zstar
c      write(*,*)'micgmbnd:uupls,uumns,upls,umns,delmu,hh,zmax(),zstar=',
c     +                    uupls,uumns,upls,umns,delmu,hh,zmax(lr_),zstar  
      psi1=psif(zstar)
      psi2=psif(zmax(lr_))
c      Write(*,*)'micgmbnd: psi1,psi2=',psi1,psi2
      xlbmd2=0.
c      write(*,*)'micgmbnd: l,ul,umx2,ustar2,alph,ab,xlbmd2='
      do 10 l=1,200
        uz=l
        ul=(uz-.5)*delmu+umns
        ustar2=(1.-psif(zstar)*(1.-ul**2))
        ustar=sqrt(ustar2)
        umx2=(1.-psif(zmax(lr_))*(1.-ul**2))
        alph=1./sqrt(1.-umx2/ustar2)
        ab=abs((alph+1)/(alph-1))
c990131        xlbmd2=xlbmd2+twopi*hh*delmu*ul/ustar*alph*alog(ab)
        xlbmd2=xlbmd2+twopi*hh*delmu*ul/ustar*alph*log(ab)
cBH091031        write(*,100) l,ul,umx2,ustar2,alph,ab,xlbmd
 100    format(i4,3(1pe16.8),3(1pe10.2))
 10   continue
      return
      end
