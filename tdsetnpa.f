c
c
      subroutine tdsetnpa(sigmanpa)
      implicit integer (i-n), real*8 (a-h,o-z)

c.......................................................................
c     This routine returns the Charge Exchange cross-section (cm**2) 
c     versus energy en_(1:nen).
c     Uses parameters and routine from IAEA ALLADIN data base:
c     http://www-amdis.iaea.org/.
c     It is much simpler than required for the SXR case, due to the
c     delta-functions in energy and scattering angle in the differential
c     cross-section, assuming zero-velocity background neutrals.
c.......................................................................

      real*8,dimension(nen_npa,npaproc)::sigmanpa
      real*8,dimension(11)::pcf
      real*8,dimension(6)::en_b,cx_b,w
      real*8,dimension(19):: wk  !Dimension 3*6+1
      dimension iop(2)
      
      
      character*30 kermsg

c     
      include 'param.h'
      include 'comm.h'


      call bcast(sigmanpa,zero,nen_npa*npaproc)

c.......................................................................
c     Loop over velocity
c.......................................................................

cBH100719:  Replacing following cross-section with
cBH100719:  formula preferred by Aaron Bader, email to BH, 100719.
c$$$      pcf(1)=-72.6656    
c$$$      pcf(2)=-5.49142
c$$$      pcf(3)=-3.42948
c$$$      pcf(4)=-1.98377
c$$$      pcf(5)=-.878009
c$$$      pcf(6)=-.198932
c$$$      pcf(7)=.0837431
c$$$      pcf(8)=.121252 
c$$$      pcf(9)=.0827182 
c$$$      pcf(10)=1.2e-1 
c$$$      pcf(11)=6.3e5
c$$$       
c$$$
c$$$cvt  do equation 16 on CompX SXR report for NPA...
c$$$cBH100423:  But only need cross-section in enmin_npa,enmax_npa range?
c$$$
c$$$      do 101 ien=1,nen_npa
c$$$
c$$$cBH100423 sigmanpa(j)=twopi*x(j)*cint2(j)*vnorm*temc1(2)**2*cxcs(endum)
c$$$cBH100423 Subroutine from IAEA Alladin web site.
c$$$         call alcheb(en_(ien)*1000d0, pcf, kdummy, cxcross, kermsg)
c$$$         if (kermsg.ne.'') then
c$$$            write(*,*)  trim(kermsg), endum,' keV, set cxcross=0.'
c$$$            cxcross=0d0
c$$$         endif
c$$$
c$$$cBH100425:  Looks like this should simply be sigmanpa(j)=cxcross
c$$$cBH100425:   sigmanpa(j)=twopi*u(j)*cint2(j)*temc1(2)**2*cxcross
c$$$         sigmanpa(ien)=cxcross
c$$$ 101  continue

c     Deuterium charge exchange cross section from Aaron Bader, 100719.
      if (npa_process(1).eq.'cxh') then
      pcf(1)=3.2345d0
      pcf(2)=235.88d0
      pcf(3)=0.038371d0
      pcf(4)=3.8068d-6
      pcf(5)=1.1832d-10
      pcf(6)=2.3713d0

      do ien=1,nen_npa
c        Adjusting energy from deuterium
c        (Agreement between Bader and IAEA Aladding p+H==>H+p is
c        much poorer if remove the factor (proton/fmass(1)) below.)
         e=en_(ien)*(proton/fmass(1))
c        Converting from m**2 to cm**2
         sigmanpa(ien,1)=1.d4*1.d-19*pcf(1)*(log(pcf(2)/e)+pcf(6))/
     +       (1.d0+pcf(3)*e+pcf(4)*e**3.5+pcf(5)*e**5.4)
      enddo
      endif  ! on npa_process(1).eq.'cxh'

      if (npaproc.ge.2 .and. npa_process(2).eq.'cxb4') then
c        From Aaron Bader memo, Bader_cql3d_cxcs_100719.pdf
         en_b(1)=75.d0
         cx_b(1)=8.8d-25
         en_b(2)=100.d0
         cx_b(2)=2.0d-24
         en_b(3)=200.d0
         cx_b(3)=5.5d-24
         en_b(4)=300.d0
         cx_b(4)=6.3e-24
         en_b(5)=400.d0
         cx_b(5)=5.9d-24
         en_b(6)=600.d0
         cx_b(6)=3.6d-24
         ncoefs=6
         iop(1)=4
         iop(2)=4
         istride=1
c        itab() and tab() are in comm.h
         itab(1)=1
         itab(2)=0
         itab(3)=0
c        Cubic spline
         call coeff1(ncoefs,en_b,cx_b,w,iop,istride,wk)
         do ien=1,nen_npa
            e=en_(ien)*(proton/fmass(1))
            if (e.lt.en_b(1)) then
               write(*,*)
               write(*,*)'tdsetnpa:  energy less than Boron table'
               write(*,*)
            endif
            if (e.lt.en_b(6)) then
               call terp1(ncoefs,en_b,cx_b,w,e,istride,tab,itab)
               sigmanpa(ien,2)=1.d4*tab(1)
            else
               sigmanpa(ien,2)=1.d4*3.05d-16*e**(-2.85405)
            endif
         enddo
      endif  ! on npa_process(2).eq.'cxb'

      if (npaproc.ge.5 .and. npa_process(5).eq.'radrecom') then
c        From Aaron Bader memo, Bader_cql3d_cxcs_100719.pdf
c        Adding Bader's n=1 and n=2 recombination CXs [BH, check]
c        Too large!
c        Using Bethe-Salpeter formula, See web, H. Poth and A. Wolf,
c        CERN-EP/82-189 (25 November 1982), 198301145.pdf, Eq. (1).
c        sigma_n=1.96*pi**2*alpha*lambda_e**2/(n*(1+n**2*Te/Ry)*Te/Ry)
c        See formula and NRL table for defn and values of constants.
c        Take Te=electron energy=Electron temperature T_elec.
c           (Perhaps (3/2)*T_elec would be better?)
c        Actually, move this calc to tdnpa, since don't have
c        the radial bin (for temperature) here.
c         do ien=1,nen_npa
c            e=en_(ien)*1.d3  !energy in eV here.
c            sigmanpa(ien,5)=1.d4*(9.95d-14/(1.d0+e/13.6d0)
c     1                           +4.98e-14/(1.d0+e/3.4d0))
c            sigmanpa(ien,5)=2.105d-22
c         enddo
      endif  ! on npa_process(5).eq.'radrecom'

c     Add more CX's here.

      do kk=1,npaproc
         if (npa_process(kk).ne.'notset'.and.kk.ne.5) then
            print*,"npa_process,en_,sigmanpa=",npa_process(kk)
            print*,(en_(ien),sigmanpa(ien,kk),ien=1,nen_npa)
         else
            print*,"npa_process,en_,sigmanpa=set elsewhere",
     +           npa_process(kk)
         endif
      enddo

      return
      end
