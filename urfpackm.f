c
c
	subroutine pack(ip, nbits, iu, nw)
c      
c       PACK - Compresses stored data
c       This version prepared for Absoft compiled data
c       (BobH, 990201).
c       SYNOPSIS
c           CALL PACK(ip, nbits, iu, nw)
c           PACK compresses stored integer data.  The following is 
c           a list of valid arguments for this routine.
c           ip        On exit, vector of packed data in range -128:127,
c                         the range of integer*1 for Absoft fortran.
c           nbits     Number of rightmost bits of data in each 
c                     partial word; must be  8
c           iu        Vector of integers in range 0:255 (beware),
c                         to be compressed.
c           nw        Number of integers to be compressed.
c
c           PACK takes the  8 rightmost bits 
c           of several partial words and concatenates them 
c           into INTEGER*1 words.
c           This implements a Cray builtin subroutine.
c           For alignment, it will be best to take nw to
c           be an integral multiple of 8, for 64-bit words.

	integer iu(nw),iu0,iu255
	integer*1 ip(nw),ip0,ip255
	
	do i=1,nw
	   ip(i)=iu(i)-128
	enddo

	return
	end



c
c
	subroutine unpack(ip, nbits, iu, nw)
c       NAME
c           UNPACK - Compresses stored data
c      SYNOPSIS
c           CALL UNPACK(p, nbits, u, nw)
c           UNPACK compresses stored integer data.  The following is 
c           a list of valid arguments for this routine.
c           ip        On exit, vector of packed data.
c           nbits     Number of rightmost bits of data in each 
c                     partial word; must be  8
c           iu        Vector of partial words to be uncompressed.
c           nw        Number of partial words to be uncompressed.
c
c           This is a modification of a Cray Subroutine.
c           For alignment, it will be best to take nw to
c           be an integral multiple of 8, for 64-bit words.

	integer iu(nw)
	integer*1 ip(nw)
	
	do i=1,nw
	   iu(i)=ip(i)+128
	enddo
	
	return
	end

	
c
c
	subroutine pack16(ip, nbits, iu, nw)
c      
c       PACK - Compresses stored data
c       This version prepared for Absoft compiled data
c       (BobH, 990201).  (lf95, BobH, 050809)
c       SYNOPSIS
c           CALL PACK(ip, nbits, iu, nw)
c           PACK compresses stored integer data.  The following is 
c           a list of valid arguments for this routine.
c           ip        On exit, vector of packed data in range -32768:32767,
c                         the range of integer*2 for lf95 fortran.
c           nbits     Number of rightmost bits of data in each 
c                     partial word; must be  8
c           iu        Vector of integers in range 0:65535 (beware),
c                         to be compressed.
c           nw        Number of integers to be compressed.
c
c           PACK takes the  8 rightmost bits 
c           of several partial words and concatenates them 
c           into INTEGER*2 words.
c           This implements a Cray builtin subroutine.
c           For alignment, it will be best to take nw to
c           be an integral multiple of 8, for 64-bit words.

	integer iu(nw),iu0,iu255
	integer*2 ip(nw),ip0,ip255
	
	do i=1,nw
	   ip(i)=iu(i)-32768
	enddo

	return
	end



c
c
	subroutine unpack16(ip, nbits, iu, nw)
c       NAME
c           UNPACK - Compresses stored data
c      SYNOPSIS
c           CALL UNPACK(p, nbits, u, nw)
c           UNPACK compresses stored integer data.  The following is 
c           a list of valid arguments for this routine.
c           ip        On exit, vector of packed data.
c           nbits     Number of rightmost bits of data in each 
c                     partial word; must be  8
c           iu        Vector of partial words to be uncompressed.
c           nw        Number of partial words to be uncompressed.
c
c           This is a modification of a Cray Subroutine.
c           For alignment, it will be best to take nw to
c           be an integral multiple of 8, for 64-bit words.

	integer iu(nw)
	integer*2 ip(nw)
	
	do i=1,nw
	   iu(i)=ip(i)+32768
	enddo
	
	return
	end

	
