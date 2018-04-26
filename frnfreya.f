c
c
      subroutine frnfreya(frmod_,fr_gyro_,beamplse_,beampon_,beampoff_,
     1     hibrz_,mfm1_,noplots)
      implicit integer (i-n), real*8 (a-h,o-z)

c.................................................................
      include 'param.h'
      include 'frcomm.h'
c     ONETWO DIVERGENCE: SEE COMMENTS AT BEGINNING OF FREYA
      character*8 frmod_,fr_gyro_,beamplse_,noplots,codeid
      real*8, intent(out) :: hibrz_(kz,ke,kb)
      integer, intent(out) :: mfm1_

      character*8 ifirst
      save ifirst
      data ifirst/"first"/

c     To pass these freya namelist in frcomm.h to comm.h
      frmod_=frmod
      fr_gyro_=fr_gyro
      beamplse_=beamplse
      beampon_=beampon
      beampoff_=beampoff
      

c..................................................................
c     This routine controls all of the NFREYA routines.
c     Called from subroutine tdintl and tdchief
c..................................................................

c..................................................................
c     Return if input variable frmod is "disabled"
c..................................................................

      if (frmod.ne."enabled") return
      elong=0

c..................................................................
c     First initialize some (iteration or time dependent) data.
c     [Possible dependence is on equilibrium, average temp, 
c      density, etc.   Check in frstup.]
c..................................................................

      call frstup(mf,mfm1,mi,mj,nion,potsid,codeid,rin,rmax,
     1  zax,zmin,zmax,zni,zne,zte,zti,zzi,p,psivol,xxx,yyy,
     1  nprim,nimp,zeffctv,zshift)

c.......................................................................
c     zshift if frstup above is passed back from comm.h data to 
c     adjust beam pivot height according to any shift in the equilibrim.
c.......................................................................

      if (ifirst.eq."first") then
         do ib=1,nbeams
            zpivot(ib)=zpivot(ib)-zshift
         enddo
         ifirst="notfirst"
      endif

c..................................................................
c     Call NFREYA (cray32 from ONETWO)
c..................................................................

      write(*,*)'frnfreya:mi,mj,codeid',
     +                    mi,mj,codeid
      write(*,*)'frnfreya:rin,rmax',
     +                    rin,rmax
      write(*,*)'frnfreya:zax,zmin,zmax',
     +                    zax,zmin,zmax
      call freya(ipts,mi,mj,codeid,rin,rmax,zax,zmin,zmax)
      write(*,*) 'Done calling freya...'
      if (ipts.eq.0) 
     1       write(*,*)'frnfreya: WARNING, ipts=0, NB missed plasma?'

c..................................................................
c     Compute the total number of particles/sec actually deposited in
c     the plasma. This requires subtracting off from the neutral current
c     the fraction that are lost at the aperture and the fraction
c     lost at the walls. ORBIT EFFECTS ARE NOT CURRENTLY CONSIDERED.
c     [2014-5: cql3d-fow now accounts for gc orbits and gyro-radius 
c      offsets.]
c     Note that in the event that Freya input variable bptor is
c     used, we are specifying the total power injected into the tokamak
c     disregarding any losses at the apertures. The expression for
c     the deposited current below takes this into account. While fap may
c     not be zero, in routine FREYA the total number of neutralized
c     particles has been normalized up to account for losses at
c     the apertures. In the event that bptor is not used, bneut is
c     not renormalized, and the expression below makes sense as it
c     stands.
c..................................................................

      curdep=0.
      do 10 ib=1,nbeams
        do 20 ie=1,3
          curdep=curdep+bneut(ie,ib)*(1.-fap(ie,ib))*(1.-fwall(ie,ib))
 20     continue
 10   continue

c..................................................................
c     Now determine the source arrays used in CQL3D
c..................................................................

      call freyasou(xpts,ypts,zpts,rpts,vx,vy,vz,ipts,curdep,
     1  bmsprd,multiply,multiplyn)


c..................................................................
c     Pass additional data to cql3d through subroutine arguments
c..................................................................

      mfm1_=mfm1
      do i=1,kz
         do ie=1,ke
            do ib=1,kb
               hibrz_(i,ie,ib)=hibrz(i,ie,ib)
            enddo
         enddo
      enddo

c..................................................................
c     Plot the FREYA birth points.
c..................................................................

        call frplteq(xpts,ypts,zpts,rpts,vx,vy,vz,ipts,curdep,
     1    nfrplt,frplt)

      return
      end
