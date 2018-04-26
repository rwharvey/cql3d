c
c
      real*8 function diagcfac(epsil,zeff)
      implicit integer (i-n), real*8 (a-h,o-z)
      save

c..................................................................
c     Cordey-Start estimates for current reduction due to electron
c     entrainment. A function of epsil and zeff.
c..................................................................

      parameter(n8=8,n12=12)
      dimension ae(n12),az(n8),af(n8,n12),afx(n8,n12),afy(n8,n12),
     c  afxy(n8,n12),iabd(4),wke(100)
      data itimes / 0 /
      data ae / 0.,.01,.02,.04,.07,.1,.2,.3,.4,.5,.6,.9 /
      data az / 1.,1.1,1.2,1.5,2.,4.,8.,16. /
      data iabd / 2,2,2,2 /
      data af / 0.,.0909,.1667,.3333,.5,.75,.875,.9375,
     c  .2242,.2892,.344,.4673,.5942,.7921,.8947,.947,
     c  .3052,.3615,.4094,.5176,.6303,.8087,.9026,.9509,
     c  .4093,.455,.4942,.5838,.6784,.8314,.9135,.9562,
     c  .5106,.5468,.5779,.65,.7274,.8551,.9251,.9619,
     c  .5824,.6122,.6379,.698,.7633,.8729,.9339,.9663,
     c  .7315,.7492,.7646,.8012,.8421,.9132,.9543,.9765,
     c  .8185,.8298,.8398,.8637,.8907,.939,.9676,.9833,
     c  .8764,.8838,.8903,.9062,.9243,.9573,.9772,.9882,
     c  .9173,.9221,.9263,.9367,.9488,.9709,.9843,.9919,
     c  .947,.95,.9527,.9592,.9669,.981,.9898,.9947,
     c  .9949,.9952,.9954,.996,.9968,.9981,.999,.9995 /
      if (itimes .eq. 1) go to 1
      itimes=1
      call bcast(afx,one,n8*n12)
      call bcast(afy,one,n8*n12)
      call bcast(afxy,zero,n8*n12)
      call coeff2(n8,az,n12,ae,af,afx,afy,afxy,n8,iabd,wke)
 1    continue
      s=terp2(zeff,epsil,n8,az,n12,ae,af,afx,afy,afxy,n8,0,0)
      diagcfac=s
      return
      end
