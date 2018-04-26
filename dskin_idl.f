
      subroutine dskin_idl(argc,argv)  !Called by IDL
      implicit integer (i-n), real*8 (a-h,o-z)
      integer*4 argc,argv(*)
      integer*4 argv_len
      include 'dskin.h'
c
c     IDL setup to read distn function data
c
c     iya=80
c     jxa=80
c     lrza=15
c     ngena=1
c     x=dblarr(jxa)
c     y=dblarr(iya)
c     rovera=dblarr(lrza)
c     elecfld=dblarr(lrza)
c     bthr=dblarr(lrza)
c     btoru=dblarr(lrza)
c     bthr0=dblarr(lrza)
c     btor0=dblarr(lrza)
c     reden=dblarr(lrza,ngena)
c     temp=dblarr(lrza,ngena)
c     radmin=0d0
c     vnorm=0d0
c     vmaxdvt=0d0
c     eovedd=0d0
c     f=dblarr(iya,jxa,lrza,ngena)
c
c     I can't see any way to get the data from the dskin.h
c       arrays to the idl variables:  The idl variables
c       are set up with specific storage.  The dskin.h variables
c       are set up with different storage locations.  Somehow, the
c       dskin.h variables need to be copied to the idl
c       variables.
c
c     If we abandon the dskin.h common variables, and pass the needed
c       storage through subroutine arguments as in the standard
c       IDL/fortran examples, of course the problem will be
c       solved.
c       
c     I'd just like to know how to set up new variables at the
c       the assigned locations.   Maybe with pointers.
c       Need to try pgf77 here.     WORKS!
c
c


      argv_len=loc(argc)         !Obtains number of arguments (argc)
                          !Because argc is passed by value,
                          !and fortran thinks is was a reference.

      open(unit=5,file='dskin_op',status='unknown')

c     The following line gives a segmentation fault
c     write(*,*) 'argc: ',argc
c     This is not too suprising, since from the fortran
c     point of view, argc is an address.  De-Referencing the address
c     (as occurs in the write statement) is probably illegal.
c     All we can do is ask what the address is: i.e., loc(argc).
c


      write(*,*) 'argv_len:', argv_len
      write(*,*) '(argv(j),j=1,argv_len): ',(argv(j),j=1,argv_len)
      write(5,*) 'argv_len:',argv_len
      write(5,*) '(argv(j),j=1,argv_len): ',(argv(j),j=1,argv_len)

c     NONE OF FOLLOWING USE OF VAL() IS LEGAL WITH PGF77. 
C     ALSO CAN'T WRITE(*,*) %VAL WITH PGF77
c     (However %VAL works with g77).
c     Check can get value of initial
c      initial=val(argv(14))
c      energy=val(argv(15))
c      pitch=val(arg(16))
c      rho=val(argv(17))
c
c      write(*,*) initial,energy,pitch,rho
c      write(*,*) %val(argv(14)),%val(argv(15)),
c     +           %val(argv(16)),%val(argv(17))


      call dskin(%VAL(argv(14)),%VAL(argv(15)),
     +           %VAL(argv(16)),%VAL(argv(17)))

      write(*,*) rovera,radmin,vnorm
      write(5,*) rovera,radmin,vnorm

      close(unit=5)


      Setting pointers to addresses for data storage set up in IDL:
      fptr=argv(13)
      roveraptr=argv(3)
      radminptr=argv(9)

c*********************************************************************
c
c     GREAT:  The following works.  Can read f, rovera, radmin in IDL.
c
c*********************************************************************
      do i=1,iya
         do j=1,jxa
            do l=1,lrza
               f_idl(i,j,l,1)=f(i,j,l,1)
            enddo
         enddo
      enddo
      do l=1,lrza
         rovera_idl(l)=rovera(l)
      enddo
      radmin_idl=radmin


      return
      end

