  module cirmsd

  contains

!#######################################################################

    subroutine ci_spawn_rmsd
    
      use sysdef
      use expec, only: cifile

      implicit none

      integer*8                           :: nspawn,nci
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

      integer*8                           :: unit,nspawn,i,j,k
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

      integer*8                           :: unit,nci,i,j,k
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

      implicit none

      integer*8                        :: nspawn,nci,i,j,imin,nother,unit
      integer*8, dimension(nci)        :: nmatch
      real*8, dimension(nspawn,natm*3) :: xspawn
      real*8, dimension(nci,natm*3)    :: xci
      real*8, dimension(natm*3)        :: cigeom,spawngeom
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
            call kabsch(cigeom,spawngeom,rmsd(i,j))
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
! kabsch: puts the Cartesian coordinates x2 into maximum coincidence
!         with the Cartesian coordinates x1 using the Kabsch algorithm
!#######################################################################

    subroutine kabsch(x1,x2,rmsd)

      use sysdef

      implicit none

      integer*8                        :: i,j,k,l
      real*8, dimension(natm*3)        :: x1,x2,rotx
      real*8, dimension(natm,3)        :: x1mat,x2mat,rmat
      real*8, dimension(3,3)           :: vut,rotmat,tmpmat
      real*8, dimension(3)             :: com1,com2,pm
      real*8                           :: totmass,det,d,rmsd,min

      integer                          :: lwork,info,dim
      double precision, dimension(3)   :: sigma
      double precision, dimension(3,3) :: covar,u,vt,orig
      double precision, dimension(15)  :: work

!-----------------------------------------------------------------------
! Translate the origins of both geometries to the centre of mass
!-----------------------------------------------------------------------
      totmass=0.0d0
      com1=0.0d0
      com2=0.0d0
      do i=1,natm
         totmass=totmass+atmass(i)
         do j=1,3
            com1(j)=com1(j)+x1(i*3-3+j)*atmass(i)
            com2(j)=com2(j)+x2(i*3-3+j)*atmass(i)
         enddo
      enddo

      com1=com1/totmass
      com2=com2/totmass

      do i=1,natm
         do j=1,3
            x1(i*3-3+j)=x1(i*3-3+j)-com1(j)
            x2(i*3-3+j)=x2(i*3-3+j)-com2(j)
         enddo
      enddo

      do i=1,natm
         do j=1,3
            x1mat(i,j)=x1(i*3-3+j)
            x2mat(i,j)=x2(i*3-3+j)
         enddo
      enddo

!-----------------------------------------------------------------------
! Calculate the 3X3 covariance matrix
!-----------------------------------------------------------------------
      covar=0.0d0
      do i=1,3
         do j=1,3
            do k=1,natm
               covar(i,j)=covar(i,j)+x1(k*3-3+i)*x2(k*3-3+j)
            enddo
         enddo
      enddo

!-----------------------------------------------------------------------
! Calculate the SVD of the covariance matrix
! N.B., This DOES work
!-----------------------------------------------------------------------
      orig=covar
      lwork=15
      dim=3

      call dgesvd('A','A',dim,dim,covar,dim,sigma,u,dim,vt,dim,work,lwork,info)

      covar=orig

      if (info.ne.0) then
         write(6,'(/,2x,a,/)') 'SVD of the covariance matrix in &
              subroutine kabsch failed'
         STOP
      endif

!-----------------------------------------------------------------------
! Calculate the rotation matrix
!-----------------------------------------------------------------------
      vut=0.0d0
      do i=1,3
         do j=1,3
            do k=1,3
               vut(i,j)=vut(i,j)+vt(k,i)*u(j,k)
            enddo
         enddo
      enddo

      det=finddet(covar,3)

      if (det.lt.0.0d0) then
         d=-1.0d0
      else
         d=1.0d0
      endif

      pm(1)=1.0d0
      pm(2)=1.0d0
      pm(3)=d

      rotmat=0.0d0
      do i=1,3
         do j=1,3
            do k=1,3
               rotmat(i,j)=rotmat(i,j)+vt(k,i)*pm(k)*u(j,k)
            enddo
         enddo
      enddo

!-----------------------------------------------------------------------
! For debugging purposes, calculate the minimum possible residual
!-----------------------------------------------------------------------
      min=0.0d0
      do k=1,natm
         do j=1,3
            min=min+x1mat(k,j)**2+x2mat(k,j)**2
         enddo
      enddo
      min=(min/natm)-(2.0d0/natm)*(sigma(1)+sigma(2)+sigma(3)*d)

!-----------------------------------------------------------------------
! Rotate the spawn geometry
!-----------------------------------------------------------------------
      do i=1,natm
         do j=1,3
            rotx(i*3-3+j)=0.0
            do k=1,3
               rotx(i*3-3+j)=rotx(i*3-3+j)+&
                    rotmat(k,j)*x2(i*3-3+k)
            enddo
         enddo
      enddo

!      rmsd=0.0d0
!      do i=1,natm*3
!         rmsd=rmsd+(x2(i)-x1(i))**2
!      enddo
!      rmsd=rmsd/natm
!      rmsd=sqrt(rmsd)
!      print*
!      print*,"RMSD before rotation:",rmsd

      rmsd=0.0d0
      do i=1,natm*3
         rmsd=rmsd+(rotx(i)-x1(i))**2
      enddo
      rmsd=rmsd/natm
      rmsd=sqrt(rmsd)
!      print*,"RMSD after rotation: ",rmsd
!      print*,"min. possible RMSD:  ",sqrt(min)

      return

    end subroutine kabsch

!#######################################################################

    function finddet(matrix,n)

      implicit none

      integer*4              :: n,i,j,k,l
      real*8, dimension(n,n) :: matrix
      real*8                 :: m,temp,finddet
      logical(kind=4)        :: detexists = .TRUE.      

      l=1

      !Convert to upper triangular form
      do k=1,n-1
         if (matrix(k,k).eq.0) then
            detexists = .false.
            do i=k+1,n
               if (matrix(i,k).ne.0) then
                  do j=1,n
                     temp=matrix(i,j)
                     matrix(i,j)=matrix(k,j)
                     matrix(k,j)=temp
                  enddo
                  detexists=.true.
                  l=-l
                  exit
               endif
            enddo
            if (.not.detexists) then
               finddet=0
               return
            endif
         endif
         do j=k+1,n
            m=matrix(j,k)/matrix(k,k)
            do i=k+1,n
               matrix(j,i)=matrix(j,i)-m*matrix(k,i)
            enddo
         enddo
      enddo
      
      !Calculate determinant by finding product of diagonal elements
      finddet=l
      do i=1,n
         finddet=finddet*matrix(i,i)
      enddo
   
    end function finddet

!#######################################################################

  end module cirmsd
