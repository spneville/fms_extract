  module iomod
    
    implicit none

    save
    integer :: iin,iout,ilog

  contains

!#######################################################################
    
    subroutine getfreeunit(inextfree)

      implicit none

      integer         :: i,inextfree
      logical(kind=4) :: lopen

! N.B. Save the first 20 io units for standard files

      do i=20,1000
         inquire(unit=i,opened=lopen)
         if (.not.lopen) then
            inextfree=i
            exit
         endif
      enddo

      return

    end subroutine getfreeunit

!#######################################################################

  end module iomod
