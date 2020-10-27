c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c                   ITERATIVE SOLVERS MODULE                           c
c----------------------------------------------------------------------c
c This Version Dated: August 13, 1996. Warning: meaning of some        c
c ============ arguments have changed w.r.t. earlier versions. Some    c
c              Calling sequences may also have changed                 c
c----------------------------------------------------------------------c 
c Contents:                                                            c
c-------------------------preconditioners------------------------------c 
c                                                                      c
c ILUT    : Incomplete LU factorization with dual truncation strategy  c
c ILUTP   : ILUT with column  pivoting                                 c
c ILUD    : ILU with single dropping + diagonal compensation (~MILUT)  c
c ILUDP   : ILUD with column pivoting                                  c
c ILUK    : level-k ILU                                                c
c ILU0    : simple ILU(0) preconditioning                              c
c MILU0   : MILU(0) preconditioning                                    c
c                                                                      c
c----------sample-accelerator-and-LU-solvers---------------------------c 
c                                                                      c
c PGMRES  : preconditioned GMRES solver                                c
c LUSOL   : forward followed by backward triangular solve (Precond.)   c
c LUTSOL  : solving v = (LU)^{-T} u (used for preconditioning)         c
c                                                                      c
c-------------------------utility-routine------------------------------c
c                                                                      c 
c QSPLIT  : quick split routine used by ilut to sort out the k largest c
c           elements in absolute value                                 c
c                                                                      c
c----------------------------------------------------------------------c
c                                                                      c 
c Note: all preconditioners are preprocessors to pgmres.               c
c usage: call preconditioner then call pgmres                          c
c                                                                      c
c----------------------------------------------------------------------c
      subroutine ilut(n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr)
c-----------------------------------------------------------------------
      implicit none 
      integer n 
cBH070419
      real*8 a(*),alu(*),w(n+1),droptol
cBH070419
      integer ja(*),ia(n+1),jlu(*),ju(n),jw(2*n),lfil,iwk,ierr
cyup      real*8 a(*),alu(*),w(*),droptol
cyup      integer ja(*),ia(*),jlu(*),ju(*),jw(*),lfil,iwk,ierr
c----------------------------------------------------------------------*
c                      *** ILUT preconditioner ***                     *
c      incomplete LU factorization with dual truncation mechanism      *
c----------------------------------------------------------------------*
c     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996  *
c----------------------------------------------------------------------*
c PARAMETERS                                                           
c-----------                                                           
c
c on entry:
c========== 
c n       = integer. The row dimension of the matrix A. The matrix 
c
c a,ja,ia = matrix stored in Compressed Sparse Row format.              
c
c lfil    = integer. The fill-in parameter. Each row of L and each row
c           of U will have a maximum of lfil elements (excluding the 
c           diagonal element). lfil must be .ge. 0.
c           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
c           EARLIER VERSIONS. 
c
c droptol = real*8. Sets the threshold for dropping small terms in the
c           factorization. See below for details on dropping strategy.
c
c  
c iwk     = integer. The lengths of arrays alu and jlu. If the arrays
c           are not big enough to store the ILU factorizations, ilut
c           will stop with an error message. 
c
c On return:
c===========
c
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c
c ju      = integer array of length n containing the pointers to
c           the beginning of each row of U in the matrix alu,jlu.
c
c ierr    = integer. Error message with the following meaning.
c           ierr  = 0    --> successful return.
c           ierr .gt. 0  --> zero pivot encountered at step number ierr.
c           ierr  = -1   --> Error. input matrix may be wrong.
c                            (The elimination process has generated a
c                            row in L or U whose length is .gt.  n.)
c           ierr  = -2   --> The matrix L overflows the array al.
c           ierr  = -3   --> The matrix U overflows the array alu.
c           ierr  = -4   --> Illegal value for lfil.
c           ierr  = -5   --> zero row encountered.
c
c work arrays:
c=============
c jw      = integer work array of length 2*n.
c w       = real work array of length n+1.
c  
c----------------------------------------------------------------------
c w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] 
c jw(n+1:2n)  stores nonzero indicators
c 
c Notes:
c ------
c The diagonal elements of the input matrix must be  nonzero (at least
c 'structurally'). 
c
c----------------------------------------------------------------------* 
c---- Dual drop strategy works as follows.                             *
c                                                                      *
c     1) Theresholding in L and U as set by droptol. Any element whose *
c        magnitude is less than some tolerance (relative to the abs    *
c        value of diagonal element in u) is dropped.                   *
c                                                                      *
c     2) Keeping only the largest lfil elements in the i-th row of L   * 
c        and the largest lfil elements in the i-th row of U (excluding *
c        diagonal elements).                                           *
c                                                                      *
c Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
c keeping  the largest  elements in  each row  of L  and U.   Taking   *
c droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
c (however, fill-in is then mpredictible).                             *
c----------------------------------------------------------------------*
c     locals
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,len 
      real*8 tnorm, t, abs, s, fact 
      if (lfil .lt. 0) goto 998
c-----------------------------------------------------------------------
c     initialize ju0 (points to next element to be added to alu,jlu)
c     and pointer array.
c-----------------------------------------------------------------------
      ju0 = n+2
      jlu(1) = ju0
c
c     initialize nonzero indicator array. 
c
      do 1 j=1,n
         jw(n+j)  = 0
 1    continue
c-----------------------------------------------------------------------
c     beginning of main loop.
c-----------------------------------------------------------------------
      do 500 ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
         tnorm = 0.0d0
         do 501 k=j1,j2
            tnorm = tnorm+abs(a(k))
 501     continue
         if (tnorm .eq. 0.0) then
            write(*,*)'ilut (in SPARSKIT):ii,ia(ii),ia(ii+1),a(j1:j2)=',
     +                                    ii,ia(ii),ia(ii+1),a(j1:j2)
            goto 999
         endif
         tnorm = tnorm/real(j2-j1+1)
c     
c     unpack L-part and U-part of row of A in arrays w 
c     
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = 0.0
         jw(n+ii) = ii
c
         do 170  j = j1, j2
            k = ja(j)
            t = a(j)
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n+k) = lenl
            else if (k .eq. ii) then
               w(ii) = t
            else
               lenu = lenu+1
               jpos = ii+lenu-1 
               jw(jpos) = k
               w(jpos) = t
               jw(n+k) = jpos
            endif
 170     continue
         jj = 0
         len = 0 
c     
c     eliminate previous rows
c     
 150     jj = jj+1
         if (jj .gt. lenl) goto 160
c-----------------------------------------------------------------------
c     in order to do the elimination in the correct order we must select
c     the smallest column index among jw(k), k=jj+1, ..., lenl.
c-----------------------------------------------------------------------
         jrow = jw(jj)
         k = jj
c     
c     determine smallest column index
c     
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            endif
 151     continue
c
         if (k .ne. jj) then
c     exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
c     exchange in jr
            jw(n+jrow) = jj
            jw(n+j) = k
c     exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
c
c     zero out element in row by setting jw(n+jrow) to zero.
c     
         jw(n+jrow) = 0
c
c     get the multiplier for row to be eliminated (jrow).
c     
         fact = w(jj)*alu(jrow)
         if (abs(fact) .le. droptol) goto 150
c     
c     combine current row and row jrow
c
         do 203 k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
            j = jlu(k)
            jpos = jw(n+j)
            if (j .ge. ii) then
c     
c     dealing with upper part.
c     
               if (jpos .eq. 0) then
c
c     this is a fill-in element
c     
                  lenu = lenu+1
                  if (lenu .gt. n) goto 995
                  i = ii+lenu-1
                  jw(i) = j
                  jw(n+j) = i
                  w(i) = - s
               else
c
c     this is not a fill-in element 
c
                  w(jpos) = w(jpos) - s

               endif
            else
c     
c     dealing  with lower part.
c     
               if (jpos .eq. 0) then
c
c     this is a fill-in element
c     
                  lenl = lenl+1
                  if (lenl .gt. n) goto 995
                  jw(lenl) = j
                  jw(n+j) = lenl
                  w(lenl) = - s
               else
c     
c     this is not a fill-in element 
c     
                  w(jpos) = w(jpos) - s
               endif
            endif
 203     continue
c     
c     store this pivot element -- (from left to right -- no danger of
c     overlap with the working elements in L (pivots). 
c     
         len = len+1 
         w(len) = fact
         jw(len)  = jrow
         goto 150
 160     continue
c     
c     reset double-pointer to zero (U-part)
c     
         do 308 k=1, lenu
            jw(n+jw(ii+k-1)) = 0
 308     continue
c     
c     update L-matrix
c     
         lenl = len 
         len = min0(lenl,lfil)
c     
c     sort by quick-split
c
         call qsplit (w,jw,lenl,len)
c
c     store L-part
c 
         do 204 k=1, len 
            if (ju0 .gt. iwk) goto 996
            alu(ju0) =  w(k)
            jlu(ju0) =  jw(k)
            ju0 = ju0+1
 204     continue
c     
c     save pointer to beginning of row ii of U
c     
         ju(ii) = ju0
c
c     update U-matrix -- first apply dropping strategy 
c
         len = 0
         do k=1, lenu-1
            if (abs(w(ii+k)) .gt. droptol*tnorm) then 
               len = len+1
               w(ii+len) = w(ii+k) 
               jw(ii+len) = jw(ii+k) 
            endif
         enddo
         lenu = len+1
         len = min0(lenu,lfil)
c
         call qsplit (w(ii+1), jw(ii+1), lenu-1,len)
c
c     copy
c 
         t = abs(w(ii))
         if (len + ju0 .gt. iwk) goto 997
         do 302 k=ii+1,ii+len-1 
            jlu(ju0) = jw(k)
            alu(ju0) = w(k)
            t = t + abs(w(k) )
            ju0 = ju0+1
 302     continue
c     
c     store inverse of diagonal element of u
c     
         if (w(ii) .eq. 0.0) w(ii) = (0.0001 + droptol)*tnorm
c     
         alu(ii) = 1.0d0/ w(ii) 
c     
c     update pointer to beginning of next row of U.
c     
         jlu(ii+1) = ju0
c-----------------------------------------------------------------------
c     end main loop
c-----------------------------------------------------------------------
 500  continue
      ierr = 0
      return
c
c     incomprehensible error. Matrix must be wrong.
c     
 995  ierr = -1
      return
c     
c     insufficient storage in L.
c     
 996  ierr = -2
      return
c     
c     insufficient storage in U.
c     
 997  ierr = -3
      return
c     
c     illegal lfil entered.
c     
 998  ierr = -4
      return
c     
c     zero row encountered
c     
 999  ierr = -5
      return
c----------------end-of-ilut--------------------------------------------
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
       subroutine pgmres(n, im, rhs, sol, vv, eps, maxits, iout,
     *                    aa, ja, ia, alu, jlu, ju, ierr)
c-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
cBH070419
       integer n, im, maxits, iout, ierr, ja(*), ia(n+1), jlu(*), ju(n)
cBH070419
       real*8 vv(n,*), rhs(n), sol(n), aa(*), alu(*), eps
cyup       integer n, im, maxits, iout, ierr, ja(*), ia(*), jlu(*), ju(*)
cyup       real*8 vv(n,*), rhs(*), sol(*), aa(*), alu(*), eps
c----------------------------------------------------------------------*
c                                                                      *
c                 *** ILUT - Preconditioned GMRES ***                  *
c                                                                      *
c----------------------------------------------------------------------*
c This is a simple version of the ILUT preconditioned GMRES algorithm. *
c The ILUT preconditioner uses a dual strategy for dropping elements   *
c instead  of the usual level of-fill-in approach. See details in ILUT *
c subroutine documentation. PGMRES uses the L and U matrices generated *
c from the subroutine ILUT to precondition the GMRES algorithm.        *
c The preconditioning is applied to the right. The stopping criterion  *
c utilized is based simply on reducing the residual norm by epsilon.   *
c This preconditioning is more reliable than ilu0 but requires more    *
c storage. It seems to be much less prone to difficulties related to   *
c strong nonsymmetries in the matrix. We recommend using a nonzero tol *
c (tol=.005 or .001 usually give good results) in ILUT. Use a large    *
c lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the    *
c more reliable the code is. Efficiency may also be much improved.     *
c Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as *
c Gaussian elimination without pivoting.                               *
c                                                                      *
c ILU(0) and MILU(0) are also provided for comparison purposes         *
c USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and *
c then call pgmres.                                                    *
c----------------------------------------------------------------------*
c Coded by Y. Saad - This version dated May, 7, 1990.                  *
c----------------------------------------------------------------------*
c parameters                                                           *
c-----------                                                           *
c on entry:                                                            *
c==========                                                            *
c                                                                      *
c n     == integer. The dimension of the matrix.                       *
c im    == size of krylov subspace:  should not exceed 50 in this      *
c          version (can be reset by changing parameter command for     *
c          kmax below)                                                 *
c rhs   == real vector of length n containing the right hand side.     *
c          Destroyed on return.                                        *
c sol   == real vector of length n containing an initial guess to the  *
c          solution on input. approximate solution on output           *
c eps   == tolerance for stopping criterion. process is stopped        *
c          as soon as ( ||.|| is the euclidean norm):                  *
c          || current residual||/||initial residual|| <= eps           *
c maxits== maximum number of iterations allowed                        *
c iout  == output unit number number for printing intermediate results *
c          if (iout .le. 0) nothing is printed out.                    *
c                                                                      *
c aa, ja,                                                              *
c ia    == the input matrix in compressed sparse row format:           *
c          aa(1:nnz)  = nonzero elements of A stored row-wise in order *
c          ja(1:nnz) = corresponding column indices.                   *
c          ia(1:n+1) = pointer to beginning of each row in aa and ja.  *
c          here nnz = number of nonzero elements in A = ia(n+1)-ia(1)  *
c                                                                      *
c alu,jlu== A matrix stored in Modified Sparse Row format containing   *
c           the L and U factors, as computed by subroutine ilut.       *
c                                                                      *
c ju     == integer array of length n containing the pointers to       *
c           the beginning of each row of U in alu, jlu as computed     *
c           by subroutine ILUT.                                        *
c                                                                      *
c on return:                                                           *
c==========                                                            *
c sol   == contains an approximate solution (upon successful return).  *
c ierr  == integer. Error message with the following meaning.          *
c          ierr = 0 --> successful return.                             *
c          ierr = 1 --> convergence not achieved in itmax iterations.  *
c          ierr =-1 --> the initial guess seems to be the exact        *
c                       solution (initial residual computed was zero)  *
c                                                                      *
c----------------------------------------------------------------------*
c                                                                      *
c work arrays:                                                         *
c=============                                                         *
c vv    == work array of length  n x (im+1) (used to store the Arnoli  *
c          basis)                                                      *
c----------------------------------------------------------------------*
c subroutines called :                                                 *
c amux   : SPARSKIT routine to do the matrix by vector multiplication  *
c          delivers y=Ax, given x  -- see SPARSKIT/BLASSM/amux         *
c lusol : combined forward and backward solves (Preconditioning ope.) *
c BLAS1  routines.                                                     *
c----------------------------------------------------------------------*
       parameter (kmax=50)
       real*8 hh(kmax+1,kmax), c(kmax), s(kmax), rs(kmax+1),t
c-------------------------------------------------------------
c arnoldi size should not exceed kmax=50 in this version..
c to reset modify paramter kmax accordingly.
c-------------------------------------------------------------
       data epsmac/1.d-16/
       n1 = n + 1
       its = 0

cTEMPORARY PRINTOUT
c$$$       iunit=41
c$$$       open(iunit)
c$$$       nnz=ia(n+1)-ia(n)
c$$$       nnz1=ju(n)-ju(1)
c$$$       write(iunit,*)'pgmres: n,nnz,nnz1,im,maxits,eps',
c$$$     +      n,nnz,nnz1,im,maxits,eps
c$$$       write(iunit,*)'rhs(1:n)'
c$$$       write(iunit,101) (rhs(i),i=1,n)
c$$$       write(iunit,*)'sol(1:n)'
c$$$       write(iunit,101) (sol(i),i=1,n)
c$$$       write(iunit,*)'aa(1:nnz)'
c$$$       write(iunit,101) (aa(i),i=1,nnz)
c$$$       write(iunit,*)'ja(1:nnz)'
c$$$       write(iunit,102) (ja(i),i=1,nnz)
c$$$       write(iunit,*)'ia(1:nnz)'
c$$$       write(iunit,102) (ia(i),i=1,n)
c$$$       write(iunit,*)'alu(1:nnz1)'
c$$$       write(iunit,101) (alu(i),i=1,nnz1)
c$$$       write(iunit,*)'jlu(1:nnz1)'
c$$$       write(iunit,102) (jlu(i),i=1,nnz1)
c$$$       write(iunit,*)'ju(1:n)'
c$$$       write(iunit,102) (ju(i),i=1,n)
c$$$       close(iunit)
c$$$ 101   format(10(1pe9.2))
c$$$ 102   format(10i9)


c-------------------------------------------------------------
c outer loop starts here..
c-------------- compute initial residual vector --------------
       call amux (n, sol, vv, aa, ja, ia)
       do 21 j=1,n
          vv(j,1) = rhs(j) - vv(j,1)
 21    continue
c-------------------------------------------------------------
 20    ro = dnrm2(n, vv, 1)
       if (iout .gt. 0 .and. its .eq. 0)
     *      write(iout, 199) its, ro
       if (ro .eq. 0.0d0) goto 999
       t = 1.0d0/ ro
       do 210 j=1, n
          vv(j,1) = vv(j,1)*t
 210   continue
       if (its .eq. 0) eps1=eps*ro
c     ** initialize 1-st term  of rhs of hessenberg system..
       rs(1) = ro
       i = 0
 4     i=i+1
       its = its + 1
       i1 = i + 1
       call lusol (n, vv(1,i), rhs, alu, jlu, ju)
       call amux (n, rhs, vv(1,i1), aa, ja, ia)
c-----------------------------------------
c     modified gram - schmidt...
c-----------------------------------------
       do 55 j=1, i
          t = ddot(n, vv(1,j),1,vv(1,i1),1)
          hh(j,i) = t
          call daxpy(n, -t, vv(1,j), 1, vv(1,i1), 1)
 55    continue
       t = dnrm2(n, vv(1,i1), 1)
       hh(i1,i) = t
       if ( t .eq. 0.0d0) goto 58
       t = 1.0d0/t
       do 57  k=1,n
          vv(k,i1) = vv(k,i1)*t
 57    continue
c
c     done with modified gram schimd and arnoldi step..
c     now  update factorization of hh
c
 58    if (i .eq. 1) goto 121
c--------perfrom previous transformations  on i-th column of h
       do 66 k=2,i
          k1 = k-1
          t = hh(k1,i)
          hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
          hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
 66    continue
 121   gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)
c
c     if gamma is zero then any small value will do...
c     will affect only residual estimate
c
       if (gam .eq. 0.0d0) gam = epsmac
c
c     get  next plane rotation
c
       c(i) = hh(i,i)/gam
       s(i) = hh(i1,i)/gam
       rs(i1) = -s(i)*rs(i)
       rs(i) =  c(i)*rs(i)
c
c     detrermine residual norm and test for convergence-
c
       hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
       ro = abs(rs(i1))
 131   format(1h ,2e14.4)
       if (iout .gt. 0)
     *      write(iout, 199) its, ro
       if (i .lt. im .and. (ro .gt. eps1))  goto 4
c
c     now compute solution. first solve upper triangular system.
c
       rs(i) = rs(i)/hh(i,i)
       do 30 ii=2,i
          k=i-ii+1
          k1 = k+1
          t=rs(k)
          do 40 j=k1,i
             t = t-hh(k,j)*rs(j)
 40       continue
          rs(k) = t/hh(k,k)
 30    continue
c
c     form linear combination of v(*,i)'s to get solution
c
       t = rs(1)
       do 15 k=1, n
          rhs(k) = vv(k,1)*t
 15    continue
       do 16 j=2, i
          t = rs(j)
          do 161 k=1, n
             rhs(k) = rhs(k)+t*vv(k,j)
 161      continue
 16    continue
c
c     call preconditioner.
c
       call lusol (n, rhs, rhs, alu, jlu, ju)
       do 17 k=1, n
          sol(k) = sol(k) + rhs(k)
 17    continue
c
c     restart outer loop  when necessary
c
       if (ro .le. eps1) goto 990
       if (its .ge. maxits) goto 991
c
c     else compute residual vector and continue..
c
       do 24 j=1,i
          jj = i1-j+1
          rs(jj-1) = -s(jj-1)*rs(jj)
          rs(jj) = c(jj-1)*rs(jj)
 24    continue
       do 25  j=1,i1
          t = rs(j)
          if (j .eq. 1)  t = t-1.0d0
          call daxpy (n, t, vv(1,j), 1,  vv, 1)
 25    continue
 199   format('   its =', i4, ' res. norm =', d20.6)
c     restart outer loop.
       goto 20
 990   ierr = 0
       return
 991   ierr = 1
       return
 999   continue
       ierr = -1
       return
c-----------------end of pgmres ---------------------------------------
c-----------------------------------------------------------------------
       end
c-----------------------------------------------------------------------
	subroutine lusol(n, y, x, alu, jlu, ju)
        real*8 x(n), y(n), alu(*)
	integer n, jlu(*), ju(*)
c-----------------------------------------------------------------------
c
c This routine solves the system (LU) x = y, 
c given an LU decomposition of a matrix stored in (alu, jlu, ju) 
c modified sparse row format 
c
c-----------------------------------------------------------------------
c on entry:
c n   = dimension of system 
c y   = the right-hand-side vector
c alu, jlu, ju 
c     = the LU matrix as provided from the ILU routines. 
c
c on return
c x   = solution of LU x = y.     
c-----------------------------------------------------------------------
c 
c Note: routine is in place: call lusol (n, x, x, alu, jlu, ju) 
c       will solve the system with rhs x and overwrite the result on x . 
c
c-----------------------------------------------------------------------
c local variables
c
        integer i,k
c
c forward solve
c
        do 40 i = 1, n
           x(i) = y(i)
           do 41 k=jlu(i),ju(i)-1
              x(i) = x(i) - alu(k)* x(jlu(k))
 41        continue
 40     continue
c
c     backward solve.
c
	do 90 i = n, 1, -1
	   do 91 k=ju(i),jlu(i+1)-1
              x(i) = x(i) - alu(k)*x(jlu(k))
 91	   continue
           x(i) = alu(i)*x(i)
 90     continue
c
  	return
c----------------end of lusol ------------------------------------------
c-----------------------------------------------------------------------
	end
c-----------------------------------------------------------------------
c----------------------------------------------------------------------- 
        subroutine qsplit(a,ind,n,ncut)
        real*8 a(n)
        integer ind(n), n, ncut
c-----------------------------------------------------------------------
c     does a quick-sort split of a real array.
c     on input a(1:n). is a real array
c     on output a(1:n) is permuted such that its elements satisfy:
c
c     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
c     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
c
c     ind(1:n) is an integer array which permuted in the same way as a(*).
c-----------------------------------------------------------------------
        real*8 tmp, abskey
        integer itmp, first, last
c-----
        first = 1
        last = n
        if (ncut .lt. first .or. ncut .gt. last) return
c
c     outer loop -- while mid .ne. ncut do
c
 1      mid = first
        abskey = abs(a(mid))
        do 2 j=first+1, last
           if (abs(a(j)) .gt. abskey) then
              mid = mid+1
c     interchange
              tmp = a(mid)
              itmp = ind(mid)
              a(mid) = a(j)
              ind(mid) = ind(j)
              a(j)  = tmp
              ind(j) = itmp
           endif
 2      continue
c
c     interchange
c
        tmp = a(mid)
        a(mid) = a(first)
        a(first)  = tmp
c
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
c
c     test for while loop
c
        if (mid .eq. ncut) return
        if (mid .gt. ncut) then
           last = mid-1
        else
           first = mid+1
        endif
        goto 1
c----------------end-of-qsplit------------------------------------------
c-----------------------------------------------------------------------
        end


c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c          BASIC MATRIX-VECTOR OPERATIONS - MATVEC MODULE              c
c         Matrix-vector Mulitiplications and Triang. Solves            c
c----------------------------------------------------------------------c
c contents: (as of Nov 18, 1991)                                       c
c----------                                                            c
c 1) Matrix-vector products:                                           c
c---------------------------                                           c
c amux  : A times a vector. Compressed Sparse Row (CSR) format.        c
c amuxms: A times a vector. Modified Compress Sparse Row format.       c
c atmux : Transp(A) times a vector. CSR format.                        c
c atmuxr: Transp(A) times a vector. CSR format. A rectangular.         c
c amuxe : A times a vector. Ellpack/Itpack (ELL) format.               c
c amuxd : A times a vector. Diagonal (DIA) format.                     c
c amuxj : A times a vector. Jagged Diagonal (JAD) format.              c
c vbrmv : Sparse matrix-full vector product, in VBR format             c
c                                                                      c
c 2) Triangular system solutions:                                      c
c-------------------------------                                       c
c lsol  : Unit Lower Triang. solve. Compressed Sparse Row (CSR) format.c
c ldsol : Lower Triang. solve.  Modified Sparse Row (MSR) format.      c
c lsolc : Unit Lower Triang. solve. Comp. Sparse Column (CSC) format.  c
c ldsolc: Lower Triang. solve. Modified Sparse Column (MSC) format.    c
c ldsoll: Lower Triang. solve with level scheduling. MSR format.       c
c usol  : Unit Upper Triang. solve. Compressed Sparse Row (CSR) format.c
c udsol : Upper Triang. solve.  Modified Sparse Row (MSR) format.      c
c usolc : Unit Upper Triang. solve. Comp. Sparse Column (CSC) format.  c
c udsolc: Upper Triang. solve.  Modified Sparse Column (MSC) format.   c
c----------------------------------------------------------------------c
c 1)     M A T R I X    B Y    V E C T O R     P R O D U C T S         c
c----------------------------------------------------------------------c
      subroutine amux (n, x, y, a,ja,ia) 
      real*8  x(*), y(*), a(*) 
      integer n, ja(*), ia(*)
c-----------------------------------------------------------------------
c         A times a vector
c----------------------------------------------------------------------- 
c multiplies a matrix by a vector using the dot product form
c Matrix A is stored in compressed sparse row storage.
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=Ax
c
c-----------------------------------------------------------------------
c local variables
c
      real*8 t
      integer i, k
c-----------------------------------------------------------------------
      do 100 i = 1,n
c
c     compute the inner product of row i with vector x
c 
         t = 0.0d0
         do 99 k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
 99      continue
c
c     store result in y(i) 
c
         y(i) = t
 100  continue
c
      return
c---------end-of-amux---------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------



c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c        BASIC LINEAR ALGEBRA FOR SPARSE MATRICES. BLASSM MODULE       c
c----------------------------------------------------------------------c
c amub   :   computes     C = A*B                                      c
c aplb   :   computes     C = A+B                                      c
c aplb1  :   computes     C = A+B  [Sorted version: A, B, C sorted]    c
c aplsb  :   computes     C = A + s B                                  c
c aplsb1 :   computes     C = A+sB  [Sorted version: A, B, C sorted]   c
c apmbt  :   Computes     C = A +/- transp(B)                          c
c aplsbt :   Computes     C = A + s * transp(B)                        c
c diamua :   Computes     C = Diag * A                                 c
c amudia :   Computes     C = A* Diag                                  c
c aplsca :   Computes     A:= A + s I    (s = scalar)                  c 
c apldia :   Computes     C = A + Diag.                                c
c----------------------------------------------------------------------c 
c Note: this module still incomplete.                                  c
c----------------------------------------------------------------------c

c-----------------------------------------------------------------------
      subroutine aplb (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1),
     *     iw(ncol)
c-----------------------------------------------------------------------
c performs the matrix sum  C = A+B. 
c-----------------------------------------------------------------------
c on entry:
c ---------
c nrow	= integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c job   = integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c 
c b, 
c jb, 
c ib	=  Matrix B in compressed sparse row format.
c
c nzmax	= integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number 
c         of elements that exceeds exceeds nzmax. See ierr.
c 
c on return:
c----------
c c, 
c jc, 
c ic	= resulting matrix C in compressed sparse row sparse format.
c	    
c ierr	= integer. serving as error message. 
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number 
c         of elements in C exceeds nzmax.
c
c work arrays:
c------------
c iw	= integer work array of length equal to the number of
c         columns in A.
c
c-----------------------------------------------------------------------
      logical values
      values = (job .ne. 0) 
      ierr = 0
      len = 0
      ic(1) = 1 
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
c     
      do 500 ii=1, nrow
c     row i 
         do 200 ka=ia(ii), ia(ii+1)-1 
            len = len+1
            jcol    = ja(ka)
            if (len .gt. nzmax) goto 999
            jc(len) = jcol 
            if (values) c(len)  = a(ka) 
            iw(jcol)= len
 200     continue
c     
         do 300 kb=ib(ii),ib(ii+1)-1
            jcol = jb(kb)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
               len = len+1
               if (len .gt. nzmax) goto 999
               jc(len) = jcol
               if (values) c(len)  = b(kb)
               iw(jcol)= len
            else
               if (values) c(jpos) = c(jpos) + b(kb)
            endif
 300     continue
         do 301 k=ic(ii), len
	    iw(jc(k)) = 0
 301     continue
         ic(ii+1) = len+1
 500  continue
      return
 999  ierr = ii
      return
c------------end of aplb ----------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------




c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c                    FORMAT CONVERSION MODULE                          c
c----------------------------------------------------------------------c
      
c-----------------------------------------------------------------------
      subroutine bndcsr (n,abd,nabd,lowd,ml,mu,a,ja,ia,len,ierr)
      real*8 a(*),abd(nabd,*), t
      integer ia(n+1),ja(*)
c----------------------------------------------------------------------- 
c Banded (Linpack ) format   to    Compressed Sparse Row  format.
c----------------------------------------------------------------------- 
c on entry:
c----------
c n	= integer,the actual row dimension of the matrix.
c
c nabd  = first dimension of array abd.
c
c abd   = real array containing the values of the matrix stored in
c         banded form. The j-th column of abd contains the elements
c         of the j-th column of  the original matrix,comprised in the
c         band ( i in (j-ml,j+mu) ) with the lowest diagonal located
c         in row lowd (see below). 
c    
c lowd  = integer. this should be set to the row number in abd where
c         the lowest diagonal (leftmost) of A is located. 
c         lowd should be s.t.  ( 1  .le.  lowd  .le. nabd).
c         The subroutines dgbco, ... of linpack use lowd=2*ml+mu+1.
c
c ml	= integer. equal to the bandwidth of the strict lower part of A
c mu	= integer. equal to the bandwidth of the strict upper part of A
c         thus the total bandwidth of A is ml+mu+1.
c         if ml+mu+1 is found to be larger than nabd then an error 
c         message is set. see ierr.
c
c len   = integer. length of arrays a and ja. bndcsr will stop if the
c         length of the arrays a and ja is insufficient to store the 
c         matrix. see ierr.
c
c on return:
c-----------
c a,
c ja,
c ia    = input matrix stored in compressed sparse row format.
c
c lowd  = if on entry lowd was zero then lowd is reset to the default
c         value ml+mu+l. 
c
c ierr  = integer. used for error message output. 
c         ierr .eq. 0 :means normal return
c         ierr .eq. -1 : means invalid value for lowd. 
c	  ierr .gt. 0 : means that there was not enough storage in a and ja
c         for storing the ourput matrix. The process ran out of space 
c         (as indicated by len) while trying to fill row number ierr. 
c         This should give an idea of much more storage might be required. 
c         Moreover, the first irow-1 rows are correctly filled. 
c
c notes:  the values in abd found to be equal to zero
c -----   (actual test: if (abd(...) .eq. 0.0d0) are removed.
c         The resulting may not be identical to a csr matrix
c         originally transformed to a bnd format.
c          
c----------------------------------------------------------------------- 
      ierr = 0
c-----------
      if (lowd .gt. nabd .or. lowd .le. 0) then 
         ierr = -1
         return
      endif
c-----------
      ko = 1
      ia(1) = 1
      do 30 irow=1,n
c-----------------------------------------------------------------------
         i = lowd 
          do  20 j=irow-ml,irow+mu
             if (j .le. 0 ) goto 19
             if (j .gt. n) goto 21
             t = abd(i,j) 
             if (t .eq. 0.0d0) goto 19
             if (ko .gt. len) then 
               ierr = irow 
               return
            endif
            a(ko) = t
            ja(ko) = j
            ko = ko+1
 19         i = i-1
 20      continue
c     end for row irow
 21      ia(irow+1) = ko
 30   continue
      return
c------------- end of bndcsr ------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c----------------------------------------------------------------------- 

      real*8 function dnrm2 ( n, dx, incx)
      integer          next
cBH100826  real*8   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one
      real*8   dx(*), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d0, 1.0d0/
c
c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
c         cuthi = minimum of  dsqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
c
      cutlo = dsqrt(EPSILON(one))
      
      if(n .gt. 0) go to 10
         dnrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      nn = n * incx
c                                                 begin main loop
      i = 1
   20 continue   ! next i
      go to next,(30, 50, 70, 110)
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 i = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200
c
  115 continue
      sum = sum + (dx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j =i,nn,incx
      if(dabs(dx(j)) .ge. hitest) go to 100
   95    sum = sum + dx(j)**2
      dnrm2 = dsqrt( sum )
      go to 300
c
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end
c----------------------------------------------------------------------



