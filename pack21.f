c
c
      subroutine pack21(a,ibot,itop,jbot,jtop,b,iy,jx)
      implicit integer (i-n), real*8 (a-h,o-z)

c.......................................................................
c     It sometimes becomes necessary to take a
c     2-D array dimensioned ibot:itop by jbot:jtop
c     and repack it as though it were 
c     dimensioned 1:iy by 1:jx, starting at a(1,1).
c     This routine does this, transfering relevant data 
c     from array a to b.
c.......................................................................

      save
      dimension a(ibot:itop,jbot:jtop)
      dimension b(iy*jx)
      do 1 j=1,jx
        i1=(j-1)*iy+1
        call dcopy(iy,a(1,j),1,b(i1),1)
 1    continue
      return
      end
c
c
      subroutine unpack21(a,ibot,itop,jbot,jtop,b,iy,jx)
      implicit integer (i-n), real*8 (a-h,o-z)

c.......................................................................
c     This is the inverse of subroutine pack21.
c     It sometimes becomes necessary to take a
c     2-D array b dimensioned 1:iy,1:jx starting at b(1,1)
c     and repack it as though it were dimensioned  
c     ibot:itop by jbot:jtop starting at a(ibot,jbot).
c     This routine does this, transfering relevant data 
c     from array b to a.
c.......................................................................

      save
      dimension a(ibot:itop,jbot:jtop)
      dimension b(iy*jx)
      zero=0.
      call bcast(a,zero,(itop-ibot+1)*(jtop-jbot+1))
      do 1 j=1,jx
        i1=(j-1)*iy+1
        call dcopy(iy,b(i1),1,a(1,j),1)
 1    continue
      return
      end
c
c
      subroutine ipack21(ia,ibot,itop,jbot,jtop,ib,iy,jx)
      implicit integer (i-n), real*8 (a-h,o-z)

c.......................................................................
c     It sometimes becomes necessary to take a
c     2-D array dimensioned ibot:itop by jbot:jtop
c     and repack it as though it were 
c     dimensioned 1:iy by 1:jx, starting at a(1,1).
c     This routine does this, transfering relevant data 
c     from array a to b.
c     This routine does this, transfering relevant data from
c     array ia to ib.
c.......................................................................

      save
      dimension ia(ibot:itop,jbot:jtop)
      dimension ib(iy*jx)
      do j=1,jx
        i1=(j-1)*iy+1
        do i=0,iy-1
           ib(i1+i)=ia(1+i,j)
        enddo
      enddo
      return
      end
