  module mcsplinemod

    implicit none
    
    save
    
    integer                                      :: nfiles,maxdat,&
                                                    maxintvl,npoints
    integer, dimension(:), allocatable           :: ndat,nintvl
    real*8, dimension(:,:), allocatable          :: dat,x,deriv,s,dx
    real*8, dimension(:,:,:), allocatable        :: coeffs
    real*8, dimension(2)                         :: erange
    character(len=80), dimension(:), allocatable :: datfile

  contains

!#######################################################################

    subroutine interpolate_stieltjes(e,f,xsec,negrid,siord,egrid,ip)
      
      implicit none

      integer                      :: negrid,siord
      real*8, dimension(3,siord-1) :: e,f
      real*8, dimension(negrid)    :: xsec
      real*8, dimension(3)         :: egrid
      real*8                       :: ip

!-----------------------------------------------------------------------
! Allocate and intitialise arrays
!-----------------------------------------------------------------------
      maxdat=siord-1
      maxintvl=maxdat-1

      allocate(deriv(3,maxdat))
      deriv=0.0d0

      allocate(s(3,maxintvl))
      s=0.0d0

      allocate(dx(3,maxintvl))
      dx=0.0d0

      allocate(coeffs(3,4,maxintvl))
      coeffs=0.0d0

      allocate(ndat(3))
      ndat=maxdat
      
      allocate(nintvl(3))
      nintvl=maxintvl
      
      dat=f

      x=e

!-----------------------------------------------------------------------
! Perform the monotonicity-constrained cubic Hermite interpolation
!-----------------------------------------------------------------------      
      call interpolate(xsec,negrid,egrid,ip)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(deriv)
      deallocate(s)
      deallocate(dx)
      deallocate(coeffs)
      deallocate(ndat)
      deallocate(nintvl)

      return

    end subroutine interpolate_stieltjes

!#######################################################################

    subroutine interpolate(xsec,negrid,egrid,ip)

      implicit none

      integer                   :: negrid
      real*8, dimension(negrid) :: xsec
      real*8, dimension(3)      :: egrid
      real*8                    :: ip
      
!-----------------------------------------------------------------------
! Compute the approximate derivatives and the slopes (which in the 
! current, crude(ish) incarnation are equal)
!-----------------------------------------------------------------------
      call getderiv

!-----------------------------------------------------------------------
! Constrain the derivatives to the Boor-Swartz monotonicity limit
!-----------------------------------------------------------------------
      call constrain_deriv

!-----------------------------------------------------------------------
! Compute the coefficients
!-----------------------------------------------------------------------
      call getcoeff

!-----------------------------------------------------------------------
! Perform the interpolation
!-----------------------------------------------------------------------
      call calcp(egrid,ip,xsec,negrid)

      return

    end subroutine interpolate

!#######################################################################

    subroutine getderiv

      implicit none
      
      integer :: i,j

      ! Loop over data sets
      do i=1,3

         ! Interval lengths for the current set
         do j=1,nintvl(i)
            dx(i,j)=x(i,j+1)-x(i,j)
         enddo

         ! Derivatives
         deriv(i,:)=0.0d0
         do j=1,nintvl(i)
            deriv(i,j)=(dat(i,j+1)-dat(i,j))/dx(i,j)
         enddo
         deriv(i,ndat(i))=deriv(i,ndat(i-1))
         
         ! Slopes within the intervals
         do j=1,nintvl(i)
            s(i,j)=(dat(i,j+1)-dat(i,j))/dx(i,j)
         enddo

      enddo

      return

    end subroutine getderiv

!#######################################################################

    subroutine constrain_deriv

      implicit none

      integer                   :: i,k
      real*8, dimension(maxdat) :: smin,smax
      real*8                    :: ftmp

      ! Loop over data sets
      do k=1,3
         
         ! Calculate the smin and smax values for the current set
         do i=2,ndat(k)-1
            smin(i)=min(s(k,i-1),s(k,i+1))
            smax(i)=max(s(k,i-1),s(k,i+1))
         enddo

         smin(1)=min(s(k,1),s(k,2))
         smin(ndat)=min(s(k,ndat-1),s(k,ndat))

         smax(1)=max(s(k,1),s(k,2))
         smax(ndat)=max(s(k,ndat-1),s(k,ndat))

         ! Constrain the derivatives for the current set according
         ! to the Boor-Swartz monotonicity limit
         do i=2,ndat(k)-1
            if (smin(i).gt.0) then
               ftmp=max(0.0d0,deriv(k,i))
               deriv(k,i)=min(ftmp,3.0d0*smin(i))
            else if (smax(i).lt.0) then
               ftmp=min(0.0d0,deriv(k,i))
               deriv(k,i)=max(ftmp,3.0d0*smax(i))
            else if (s(k,i-1)*s(k,i).le.0) then
               deriv(k,i)=0.0d0
            endif
         enddo
         
      enddo

      return

    end subroutine constrain_deriv

!#######################################################################

    subroutine getcoeff
      
      implicit none
      
      integer :: k,i

      ! Loop over data sets
      do k=1,3

         ! Loop over the intervals for the current set
         do i=1,nintvl(k)

            ! C1
            coeffs(k,1,i)=dat(k,i)
            
            ! C2
            coeffs(k,2,i)=deriv(k,i)
            
            ! C3
            coeffs(k,3,i)=3.0d0*s(k,i)-deriv(k,i+1)-2.0d0*deriv(k,i)
            coeffs(k,3,i)=coeffs(k,3,i)/dx(k,i)
            
            ! C4
            coeffs(k,4,i)=-(2.0d0*s(k,i)-deriv(k,i+1)-deriv(k,i))
            coeffs(k,4,i)=coeffs(k,4,i)/dx(k,i)**2
            
         enddo

      enddo

      return

    end subroutine getcoeff

!#######################################################################

    subroutine calcp(egrid,ip,xsec,negrid)

      implicit none

      integer                   :: negrid,k,i,iint
      real*8, dimension(3)      :: egrid
      real*8, dimension(negrid) :: xsec
      real*8                    :: ip,diff,p,xcurr
      
      diff=(egrid(2)-egrid(1))/egrid(3)

      ! Loop over the grid points
      npoints=int(egrid(3))
      do i=1,npoints
         
         ! Set energy at the current grid point
         xcurr=egrid(1)+(i-1)*diff

         ! Cycle if the current energy is below the IP
         if (xcurr.lt.ip) cycle

         ! Calculate the interpolant at the current gridpoint
         p=0.0d0
         do k=1,3
            
            ! Cycle if we are out of bounds
            if (xcurr.lt.x(k,1)) cycle
            if (xcurr.gt.x(k,ndat(k))) cycle
            
            ! Determine which interval we are in
            call getiint(iint,xcurr,k)
            
            ! Compute the value of the current interpolant at the 
            ! current value of x
            p=p+coeffs(k,1,iint) &
                 + coeffs(k,2,iint)*(xcurr-x(k,iint)) &
                 + coeffs(k,3,iint)*(xcurr-x(k,iint))**2 &
                 + coeffs(k,4,iint)*(xcurr-x(k,iint))**3

         enddo

         ! Set the cross-section at the current grid point
         xsec(i)=p

      enddo

      return

    end subroutine calcp

!#######################################################################

    subroutine getiint(iint,xcurr,k)

      implicit none

      integer :: iint,k,i
      real*8  :: xcurr

      do i=1,nintvl(k)
         if (xcurr.ge.x(k,i).and.xcurr.lt.x(k,i+1)) then
            iint=i
            exit
         endif
      enddo

      return

    end subroutine getiint

!#######################################################################

  end module mcsplinemod
