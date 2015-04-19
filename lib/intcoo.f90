  module intcoo

    implicit none

  contains

!#######################################################################
! x2int: for a given set of Cartesian coordinates xcoo returns the
!        value of a single internal coordinate intcoo
!
! ityp=1 <-> bond length
!      2 <-> angle
!      3 <-> dihedral
!      4 <-> twisting angle
!#######################################################################

    function x2int(xcoo) result(intcoo)
      
      use expec
      use sysdef

      implicit none

      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo

      select case(ityp)

      case(1) ! Bond length
         call calc_length(xcoo,intcoo)

      case(2) ! Bond angle
         call calc_angle(xcoo,intcoo)

      case(3) ! Dihedral angle
         call calc_dihedral(xcoo,intcoo)

      case(4) ! Twisting angle
         call calc_twist(xcoo,intcoo)

      end select

      return

    end function x2int
    
!#######################################################################

    subroutine calc_length(xcoo,intcoo)

      use expec
      use sysdef

      implicit none

      integer*8                 :: i,k1,k2
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo
      
      intcoo=0.0d0

      k1=iatm(1)
      k2=iatm(2)
      
      do i=1,3
         intcoo=intcoo+(xcoo(k1*3-3+i)-xcoo(k2*3-3+i))**2
      enddo

      intcoo=sqrt(intcoo)
      
      return

    end subroutine calc_length

!#######################################################################

    subroutine calc_angle(xcoo,intcoo)

      use expec
      use sysdef

      implicit none

      integer*8                 :: i,k1,k2,k3
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,dp,len1,len2,pi
      real*8, dimension(3)      :: vec1,vec2

      pi=4.0d0*datan(1.0d0)

      k1=iatm(1)
      k2=iatm(2)
      k3=iatm(3)

      len1=0.0d0
      len2=0.0d0

      do i=1,3
         vec1(i)=xcoo(k2*3-3+i)-xcoo(k1*3-3+i)         
         vec2(i)=xcoo(k2*3-3+i)-xcoo(k3*3-3+i)
         len1=len1+vec1(i)**2
         len2=len2+vec2(i)**2
      enddo

      len1=sqrt(len1)
      len2=sqrt(len2)

      dp=dot_product(vec1,vec2)

      intcoo=dp/(len1*len2)
      intcoo=dacos(intcoo)
      intcoo=intcoo*180.0d0/pi

      return

    end subroutine calc_angle

!#######################################################################

    subroutine calc_dihedral(xcoo,intcoo)

      use expec
      use sysdef

      implicit none

      integer*8                 :: i,k1,k2,k3,k4
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,pi,ftmp1,ftmp2
      real*8, dimension(3)      :: vec1,vec2,vec3,tmp1,tmp2,tmp3,tmp4

      pi=4.0d0*datan(1.0d0)

      k1=iatm(1)
      k2=iatm(2)
      k3=iatm(3)
      k4=iatm(4)

      do i=1,3
         vec1(i)=xcoo(k2*3-3+i)-xcoo(k1*3-3+i)
         vec2(i)=xcoo(k3*3-3+i)-xcoo(k2*3-3+i)
         vec3(i)=xcoo(k4*3-3+i)-xcoo(k3*3-3+i)
      enddo

      tmp1=cross_product(vec1,vec2)
      tmp2=cross_product(vec2,vec3)
      
      ftmp1=dot_product(vec2,vec2)
      ftmp1=sqrt(ftmp1)
      tmp3=vec2/ftmp1
      
      tmp4=cross_product(tmp1,tmp2)
      ftmp1=dot_product(tmp4,tmp3)

      ftmp2=dot_product(tmp1,tmp2)

      intcoo=atan2(ftmp1,ftmp2)

      intcoo=intcoo*180.0d0/pi

      return

    end subroutine calc_dihedral

!#######################################################################

    subroutine calc_twist(xcoo,intcoo)

      use expec
      use sysdef
      
      implicit none

      integer*8                 :: i,k1,k2,k3,k4,k5,k6,k7,k8
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,pi,dot,&
                                   len1,len2,len3,len4
      real*8, dimension(3)      :: vec1,vec2,vec3,vec4,cross1,cross2

      pi=4.0d0*datan(1.0d0)

      k1=iatm(1)
      k2=iatm(2)
      k3=iatm(3)
      k4=iatm(4)
      k5=iatm(5)
      k6=iatm(6)
      k7=iatm(7)
      k8=iatm(8)
      
      len1=0.0d0
      len2=0.0d0
      len3=0.0d0
      len4=0.0d0
      
      do i=1,3
         vec1(i)=xcoo(k2*3-3+i)-xcoo(k1*3-3+i)
         vec2(i)=xcoo(k4*3-3+i)-xcoo(k3*3-3+i)
         vec3(i)=xcoo(k6*3-3+i)-xcoo(k5*3-3+i)
         vec4(i)=xcoo(k8*3-3+i)-xcoo(k7*3-3+i)
         len1=len1+vec1(i)**2
         len2=len2+vec2(i)**2
         len3=len3+vec3(i)**2
         len4=len4+vec4(i)**2
      enddo
      len1=sqrt(len1)
      len2=sqrt(len2)
      len3=sqrt(len3)
      len4=sqrt(len4)
      vec1=vec1/len1
      vec2=vec2/len2
      vec3=vec3/len3
      vec4=vec4/len4

      cross1=cross_product(vec1,vec2)
      cross2=cross_product(vec3,vec4)

      dot=dot_product(cross1,cross2)

      intcoo=acos(dot)

      intcoo=intcoo*180.0d0/pi

      return

    end subroutine calc_twist

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

  end module intcoo
