

      program diskin_main
      implicit integer (i-n), real*8 (a-h,o-z)
      include 'dskin.h'
c
c     Program to check out dskin.f read of diskf CQL3D o/p file.
c
c     Compile and load (depending on installed compilers:
c                            g77 -o dskin_main dskin_main.f dskin.f
c                       OR,  make -f dskin_main_makefile
c                       OR,  pgf77 -o dskin_main dskin_main.f dskin.f 
c                       OR,  pgf90 -o dskin_main dskin_main.f dskin.f 
c

      pi=4d0*atan2(1.,1.)

      initial=1               !Initialize data

      call dskin(initial,energy,pitch,rho,fdist)

c  Write some data to check everything got read

      write(*,*) 'rovera=',(rovera(ll),ll=1,lrza)
      write(*,*) 'bthr=',(bthr(ll),ll=1,lrza)
      write(*,*) 'elecfld=',(elecfld(ll),ll=1,lrza)
      write(*,*) 'radmin,vnorm,vmaxdvt,eovedd',
     +     radmin,vnorm,vmaxdvt,eovedd
      write(*,*)'Distn function f(x,ll=1)= ',(f(1,j,1,1),j=1,jxa)
      write(*,*)'Distn function f(x,ll=lrza)= ',(f(1,j,lrza,1),j=1,jxa)

      initial=0               !Interpolate data

      rho=0.2499
      energy=1300.
      pitch=0.45*pi
      call dskin(initial,energy,pitch,rho,fdist1)

      rho=0.2499
      energy=1300.
      pitch=0.55*pi
      call dskin(initial,energy,pitch,rho,fdist2)

      rho=0.2499
      energy=1300.
      pitch=0.1*pi
      call dskin(initial,energy,pitch,rho,fdist3)

      rho=0.2499
      energy=1300.
      pitch=0.9*pi
      call dskin(initial,energy,pitch,rho,fdist4)

      write(*,*) 'diskin_main: fdist1->4:',fdist1,fdist2,fdist3,fdist4 

      call exit
      end
