  module cirmsd

  contains

!#######################################################################

    subroutine ci_spawn_rmsd
    
      use sysdef
      use expec, only: cifile

      implicit none

      integer                             :: nspawn,nci
      real*8, dimension(:,:), allocatable :: xspawn,xci

!-----------------------------------------------------------------------
! Read the spawn geometries from file
!-----------------------------------------------------------------------
      call rdspawnfile(xspawn,nspawn)

!-----------------------------------------------------------------------
! Read the conical intersection geometries from file
!-----------------------------------------------------------------------
      call rdcifile(xci,nci)

!-----------------------------------------------------------------------
! Calculate the RMSD between each pair of spawn and conical intersection
! geometries using the Kabsch algorithm
!-----------------------------------------------------------------------
      call calc_rmsd(xspawn,xci,nspawn,nci)

      return

    end subroutine ci_spawn_rmsd

!#######################################################################

    subroutine rdspawnfile(xspawn,nspawn)

      use sysdef

      implicit none

      integer                             :: unit,nspawn,i,j,k
      real*8, dimension(:,:), allocatable :: xspawn

!-----------------------------------------------------------------------
! Determine the number of spawn geometries and allocate the xspawn
! array
!-----------------------------------------------------------------------
      unit=20
      open(unit,file='spawngeom.xyz',form='formatted',status='old')
      
      nspawn=0
10    continue
      do i=1,natm+2
         read(unit,*,end=15)
      enddo
      nspawn=nspawn+1
      goto 10

15    continue
      
      allocate(xspawn(nspawn,natm*3))

!-----------------------------------------------------------------------
! Read the spawn geometries from file
!-----------------------------------------------------------------------
      rewind(unit)

      do i=1,nspawn
         read(unit,*)
         read(unit,*)
         do j=1,natm
            read(unit,'(2x,3(3x,F13.9))') (xspawn(i,k), k=j*3-2,j*3)
         enddo
      enddo

      close(unit)

      return

    end subroutine rdspawnfile

!#######################################################################

    subroutine rdcifile(xci,nci)

      use sysdef
      use expec, only: cifile

      implicit none

      integer                             :: unit,nci,i,j,k
      real*8, dimension(:,:), allocatable :: xci
      character(len=10)                   :: atmp

!-----------------------------------------------------------------------
! Determine the number of CIs geometries and allocate the xci array
!-----------------------------------------------------------------------
      unit=20
      open(unit,file=cifile,form='formatted',status='old')
      
      nci=0
10    continue
      do i=1,natm+2
         read(unit,*,end=15)
      enddo
      nci=nci+1
      goto 10

15    continue
      
      allocate(xci(nci,natm*3))

!-----------------------------------------------------------------------
! Read the CI geometries from file
!-----------------------------------------------------------------------
      rewind(unit)

      do i=1,nci
         read(unit,*)
         read(unit,*)
         do j=1,natm
            read(unit,*) atmp,(xci(i,k), k=j*3-2,j*3)
         enddo
      enddo

      close(unit)

      return

    end subroutine rdcifile

!#######################################################################

    subroutine calc_rmsd(xspawn,xci,nspawn,nci)

      use sysdef
      use kabschmod
      
      implicit none

      integer                          :: nspawn,nci,i,j,imin,nother,unit
      integer, dimension(nci)          :: nmatch
      real*8, dimension(nspawn,natm*3) :: xspawn
      real*8, dimension(nci,natm*3)    :: xci
      real*8, dimension(natm*3)        :: cigeom,spawngeom,rotx
      real*8, dimension(nspawn,nci)    :: rmsd
      real*8                           :: fmin,pc,sum
      real*8, parameter                :: tol=2.5d0
      character(len=120)               :: atmp

      nmatch=0

      nother=0

      ! Loop over spawn geometries
      do i=1,nspawn
         fmin=9999999.9d0
         ! Loop over ci geometries
         do j=1,nci
            ! Put the current spawn geometry into maximum coincidence
            ! with the current ci geometry
            cigeom=xci(j,:)
            spawngeom=xspawn(i,:)
            call kabsch(cigeom,spawngeom,rmsd(i,j),rotx)
            ! If the current ci geometry is the closest yet to the
            ! curret spawn geometry, then set imin=j
            if (rmsd(i,j).lt.fmin.and.rmsd(i,j).lt.tol) then
               imin=j
               fmin=rmsd(i,j)
            endif
         enddo
         if (fmin.lt.tol) then
            nmatch(imin)=nmatch(imin)+1
         else
            nother=nother+1
         endif
      enddo

      ! Output results      
      write(6,*)
      sum=0.0d0
      do i=1,nci
         pc=100.0d0*real(nmatch(i))/real(nspawn)
         sum=sum+pc
         write(6,'(2x,a,x,i2,a1,2x,i3,x,a1,i2,a2)') &
              'Number of spawn geometries closest to CI geometry',&
              i,':',nmatch(i),'(',nint(pc),'%)'
      enddo
      write(6,'(2x,a,38x,i3,x,a1,i2,a2)') &
              'Other geometries:',nother,'(',nint(100.0d0-sum),'%)'
      write(6,*)

      unit=20
      open(unit,file='cirmsd.log',form='formatted',status='unknown')
      atmp='Spawn geometry | RMSD_1 | ... | RMSD_NCI'
      write(unit,'(a)') adjustl(trim(atmp))
      do i=1,nspawn
         write(unit,*) i,(rmsd(i,j),j=1,nci)
      enddo
      close(unit)

      return

    end subroutine calc_rmsd

!#######################################################################

  end module cirmsd
