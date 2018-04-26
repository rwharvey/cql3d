c
c
      subroutine firstdrv(x,f,fprim,n)
      implicit integer (i-n), real*8 (a-h,o-z)
      dimension x(n),f(n),fprim(n)

c..................................................................
c     Compute the first derivative of the function f(x,l_) at the x(l)
c     mesh points and put the solution into the array fprim(l), for
c     l=1,...,n     Taylor expansion used is:
c     f(x+h(i,l_))=f(x,l_) + h(i)*f'(x) + h(i)**2/2*f"(x) for i=1,2
c     the x mesh need not be evenly spaced.
c..................................................................

      do 10 l=2,n-1
        h1=x(l-1)-x(l)
        h2=x(l+1)-x(l)
        a=h2**2/h1**2
        fprim(l)=(f(l-1)*a - f(l+1) - f(l)*(a-1))/(a*h1-h2)
 10   continue

c..................................................................
c     End points extrapolate, but expansion is same as above.
c..................................................................

      h1=x(2)-x(1)
      h2=x(3)-x(1)
      a=h2**2/h1**2
      fprim(1)=(f(2)*a - f(3) - f(1)*(a-1))/(a*h1-h2)
      h1=x(n-1)-x(n)
      h2=x(n-2)-x(n)
      a=h2**2/h1**2
      fprim(n)=(f(n-1)*a - f(n-2) - f(n)*(a-1))/(a*h1-h2)
      return
      end
