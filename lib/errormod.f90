  module errormod

    implicit none

  contains

!#######################################################################
!
! errcntrl: writes the passed string to the screen and, if open, the 
!           log file, then terminates the program
!
!#######################################################################

    subroutine errcntrl(string)

      use iomod

      implicit none
      
      integer         :: ilbl
      character(*)    :: string
      logical(kind=4) :: lopen

      ! Write error message to the screen
      write(6,'(/,2x,a,/)') trim(string)

!      ! If a log file is open, write the error message to the log file
!      inquire(unit=ilog,opened=lopen)
!      if (lopen) write(ilog) string(1:ilbl)

      ! Terminate the program
      STOP

    end subroutine errcntrl

!#######################################################################

  end module errormod
