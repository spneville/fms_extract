  module projmod

    implicit none

  contains

!#######################################################################

    subroutine rdseamfiles2

      use sysdef
      use expec

      implicit none

      integer                     :: unit,i,j,k
      real*8                      :: norm,dp,b2a
      real*8, dimension(2,natm*3) :: grad
      character(len=2)            :: atmp

!----------------------------------------------------------------------- 
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(cigeom(3*natm))
      allocate(branchvec(2,3*natm))

!-----------------------------------------------------------------------
! CI geometry 
!-----------------------------------------------------------------------
      unit=335
      open(unit,file=cifile,form='formatted',status='old')
      read(unit,*)
      read(unit,*)
      do i=1,natm
         read(unit,*) atmp,(cigeom(j),j=i*3-2,i*3)
      enddo
      close(unit)

      ! Convert to a.u.
      cigeom=cigeom/0.529177249d0

!-----------------------------------------------------------------------
! NACT vector
!-----------------------------------------------------------------------
      open(unit,file=hfile,form='formatted',status='old')
      do i=1,natm
         read(unit,*) (branchvec(1,j),j=i*3-2,i*3)
      enddo
      close(unit)

!-----------------------------------------------------------------------
! Gradient difference vector
!-----------------------------------------------------------------------
      do k=1,2
         open(unit,file=gfile(k),form='formatted',status='old')
         do i=1,natm
            read(unit,*) (grad(k,j),j=i*3-2,i*3)
         enddo
         close(unit)
         norm=sqrt(dot_product(grad(k,:),grad(k,:)))
         grad(k,:)=grad(k,:)/norm
      enddo
      branchvec(2,:)=grad(1,:)-grad(2,:)

!-----------------------------------------------------------------------
! Output the normalised h- and g-vectors
!-----------------------------------------------------------------------
      open(unit,file='branchingvecs.xyz',form='formatted',&
           status='unknown')
      
      b2a=0.529177249d0

      write(unit,'(i2)') natm
      write(unit,'(a)') 'h-vector'
      norm=sqrt(dot_product(branchvec(1,:),branchvec(1,:)))
      do i=1,natm
         write(unit,'(a2,6(2x,F10.7))') atlbl(i),&
              (cigeom(j)*b2a,j=i*3-2,i*3),(branchvec(1,j)/norm,j=i*3-2,i*3)
      enddo

      write(unit,'(i2)') natm
      write(unit,'(a)') 'g-vector'
      norm=sqrt(dot_product(branchvec(2,:),branchvec(2,:)))
      do i=1,natm
         write(unit,'(a2,6(2x,F10.7))') atlbl(i),&
              (cigeom(j)*b2a,j=i*3-2,i*3),(branchvec(2,j)/norm,j=i*3-2,i*3)
      enddo

      close(unit)

!-----------------------------------------------------------------------
! Orthonormalisation
!-----------------------------------------------------------------------
      dp=dot_product(branchvec(2,:),branchvec(1,:))
      branchvec(2,:)=branchvec(2,:)-dp*branchvec(1,:)/dot_product(branchvec(1,:),branchvec(1,:))
      dp=dot_product(branchvec(2,:),branchvec(1,:))
      branchvec(2,:)=branchvec(2,:)-dp*branchvec(1,:)/dot_product(branchvec(1,:),branchvec(1,:))

      do k=1,2
         norm=sqrt(dot_product(branchvec(k,:),branchvec(k,:)))
         branchvec(k,:)=branchvec(k,:)/norm
      enddo

      return

    end subroutine rdseamfiles2

!#######################################################################

    subroutine calc_projector2

      use sysdef
      use expec

      implicit none

      integer :: i,j,k

!-----------------------------------------------------------------------
! Allocate arrays 
!-----------------------------------------------------------------------
      allocate(branchproj(3*natm,3*natm))

!-----------------------------------------------------------------------
! Set up the projector onto the branching space
!-----------------------------------------------------------------------
      branchproj=0.0d0
      do k=1,2
         do i=1,3*natm
            do j=1,3*natm
               branchproj(i,j)=branchproj(i,j)+&
                    branchvec(k,i)*branchvec(k,j)
            enddo
         enddo
      enddo

      return

    end subroutine calc_projector2

!#######################################################################

  end module projmod
