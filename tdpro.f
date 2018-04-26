c
c
      real*8 function tdpro(psi,rplasm,acoef)
      implicit integer (i-n), real*8 (a-h,o-z)
c..................................................................
c     Calculate ASDEX YAG1 type profiles.
c     acoef(i),i=1,4  must be specified.
c..................................................................

      dimension acoef(4)

      x = rplasm*psi
      x2=x*x
      arg = acoef(4)
      do 10  i=1,3
        arg = arg*x2 + acoef(4-i)
 10   continue
      tdpro = exp(arg)
      return
      end
