      subroutine lam(en_m,ene,stoplamda)
c.............................................
c This subroutine determines the neutral mean free path as a function of 
c its energy.
c It uses beam-stopping cross sections from tdnpabscs.f and assumes that
c the energtic neutral has energy greater than Ti and Te.
c The results are stored in array stoplamda(1:lrzmax,1:n_en_length) in cm
c Inputs: en_ is an array of energetic neutral energies with length nen
c en_ is the left edge of energy bin....
c Version 1.0
c Outputs checked and looks reasonable.
c..............V.Tang 9-25-05...............................................
       implicit integer (i-n), real*8 (a-h,o-z)
       include 'param.h'
       include 'comm.h'
       real*8,dimension(lrzmax,nen_npa)::stoplamda
       real*8,dimension(nen_npa)::en_m 
       real*8::stopcs
       real*8,dimension(lrzmax)::ene
     
c       print*, "lam called"
       print*,lrzmax, nen_npa
       do 200 nn=1,lrzmax !flux surface
       do 100 mm=1,nen_npa !energy
c$$$       print*,'lam: loop enetered'
c$$$       print*,'lam: electron den is:',ene(nn)
c$$$       print*,'lam: temp(kelec,nn)is:',temp(kelec,nn)
c$$$       print*,'lam: zeff(nn)is:',zeff(nn)
c$$$       print*,'lam: fast neut engy:',en_m(mm)/bnumb(1)
       stopcs=0.
       call bscs(en_m(mm)/bnumb(1),ene(nn),
     1  temp(kelec,nn),zeff(nn),stopcs) ! calls bscs to get cross section
c$$$       print*,'cross section:', stopcs
       stoplamda(nn,mm)=1/(ene(nn)*stopcs) ! mean-free path=1/(ne*sigma_stop) in cm
c$$$       print *,"lam: flux surface",nn,"Energy index",mm
c$$$       print *,"lam: Mean free path in cm",stoplamda(nn,mm)
100    continue
200    continue
       return
       END
