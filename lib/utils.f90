  module utils
  
    implicit none

  contains

!#######################################################################

    function mass(label)

      implicit none
      
      real*8           :: mass
      character(len=*) :: label

      if (label.eq.'H') then
         mass=1.00782504d0
      else if (label.eq.'C') then
         mass=12.00000000d0
      else
         write(6,'(2(2x,a))') 'Unknown atom type:',trim(label)
         STOP
      endif

      return

    end function mass

!#######################################################################

    function atnum(label)

      implicit none

      real*8           :: atnum
      character(len=*) :: label

      if (label.eq.'H') then
         atnum=1.0d0
      else if (label.eq.'C') then
         atnum=6.0d0
      else
         write(6,'(2(2x,a))') 'Unknown atom type:',trim(label)
         STOP
      endif

      return

    end function atnum

!#######################################################################

  end module utils
