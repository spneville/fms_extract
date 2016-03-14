  module mathlib

  contains

!#######################################################################
    
    function factorial(num) result(fac)

      implicit none

      integer :: num,fac,i

      fac=1
      do i=2,num
         fac=fac*i
      enddo
      
      return
      
    end function factorial

!#######################################################################

    function cross_product(vec1,vec2) result(cp)

      implicit none

      real*8, dimension(3) :: vec1,vec2
      real*8, dimension(3) :: cp
      
      cp(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
      cp(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
      cp(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)

      return

    end function cross_product

!#######################################################################

  end module mathlib
