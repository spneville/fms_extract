  module intcoo

    implicit none

    save 
    
    integer, dimension(10) :: atmindx

  contains

!#######################################################################
! x2int: for a given set of Cartesian coordinates xcoo returns the
!        value of a single internal coordinate intcoo
!
! ityp=1 <-> bond length
!      2 <-> angle
!      3 <-> dihedral
!      4 <-> twisting angle
!      5 <-> pyramidalisation angle
!      6 <-> maximum pyramidalisation angle over two groups
!      7 <-> maximum of two angles
!      8 <-> maximum of two angles
!      9 <-> absolute difference between two angles
!     10 <-> absolute difference between to pyramidalisation
!            angles
!     -1 <-> Cartesian vector (momentum representation only)
!     -2 <-> Distance from a CI projected onto the branching space
!#######################################################################

    function x2int(xcoo,intnum) result(intcoo)
      
      use expec
      use sysdef

      implicit none

      integer                   :: intnum,inttyp
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo
      
      if (intnum.eq.1) then
         inttyp=ityp
         atmindx=iatm         
      else if (intnum.eq.2) then
         inttyp=ityp2
         atmindx=iatm2
      endif

      select case(inttyp)

      case(1) ! Bond length
         call calc_length(xcoo,intcoo)

      case(2) ! Bond angle
         call calc_angle(xcoo,intcoo)

      case(3) ! Dihedral angle
         call calc_dihedral(xcoo,intcoo)

      case(4) ! Twisting angle
         call calc_twist(xcoo,intcoo)

      case(5) ! Pyramidalisation angle
         call calc_pyr(xcoo,intcoo)

      case(6) ! Pyramidalisation angle, two groups
         call calc_pyr2(xcoo,intcoo)

      case(7) ! Maxiumum of two angles
         call calc_maxangle(xcoo,intcoo)

      case(8) ! Minimum of two angles
         call calc_minangle(xcoo,intcoo)

      case(9) ! Absolute difference between two angles
         call calc_diffangle(xcoo,intcoo)

      case(10) ! Absolute difference between two pyramidalisation
               ! angles
         call calc_diffpyr(xcoo,intcoo)

      case(-2) ! Distance from a CI projected onto the branching
               ! space
         call calc_dist_seam(xcoo,intcoo)
         
      end select

      return

    end function x2int

!#######################################################################

    function x2int_mom(pcoo,xcoo) result(intmom)

      use expec
      use sysdef

      implicit none
      
      real*8, dimension(natm*3) :: pcoo,xcoo
      real*8                    :: intmom

      select case(ityp)

      case(-1) ! Cartesian vector
         call calc_cart(pcoo,xcoo,intmom)
         return
      end select
         
      write(6,'(/,2x,a,/)') 'The selected nternal coordinate type is &
           not currently supported in the momentum representation'
      STOP

    end function x2int_mom
    
!#######################################################################

    subroutine calc_length(xcoo,intcoo)

      use expec
      use sysdef

      implicit none

      integer                   :: i,k1,k2
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo
      
      intcoo=0.0d0

      k1=atmindx(1)
      k2=atmindx(2)

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

      integer                   :: i,k1,k2,k3
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,dp,len1,len2,pi
      real*8, dimension(3)      :: vec1,vec2

      pi=4.0d0*datan(1.0d0)

      k1=atmindx(1)
      k2=atmindx(2)
      k3=atmindx(3)

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
      use mathlib

      implicit none

      integer                   :: i,k1,k2,k3,k4
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,pi,ftmp1,ftmp2
      real*8, dimension(3)      :: vec1,vec2,vec3,tmp1,tmp2,tmp3,tmp4

      pi=4.0d0*datan(1.0d0)

      k1=atmindx(1)
      k2=atmindx(2)
      k3=atmindx(3)
      k4=atmindx(4)

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
      use mathlib
      
      integer                   :: i,k1,k2,k3,k4,k5,k6,k7,k8
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,pi,dot,len1,len2,len3,len4
      real*8, dimension(3)      :: vec1,vec2,vec3,vec4,cross1,cross2
      
      pi=4.0d0*datan(1.0d0)

      k1=atmindx(1)
      k2=atmindx(2)
      k3=atmindx(3)
      k4=atmindx(4)
      k5=atmindx(5)
      k6=atmindx(6)
      k7=atmindx(7)
      k8=atmindx(8)
      
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

    subroutine calc_pyr(xcoo,intcoo)

      use expec
      use sysdef
      use mathlib

      integer                   :: k1,k2,k3,k4,i
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,pi
      real*8, dimension(3)      :: r14,r23,r13,r12,r14xr23,r13xr12
      
      pi=4.0d0*datan(1.0d0)

      k1=atmindx(1)
      k2=atmindx(2)
      k3=atmindx(3)
      k4=atmindx(4)

      do i=1,3
         r14(i)=xcoo(k4*3-3+i)-xcoo(k1*3-3+i)
         r23(i)=xcoo(k3*3-3+i)-xcoo(k2*3-3+i)
         r13(i)=xcoo(k3*3-3+i)-xcoo(k1*3-3+i)
         r12(i)=xcoo(k2*3-3+i)-xcoo(k1*3-3+i)         
      enddo
      r14=r14/sqrt(dot_product(r14,r14))
      r23=r23/sqrt(dot_product(r23,r23))
      r13=r13/sqrt(dot_product(r13,r13))
      r12=r12/sqrt(dot_product(r12,r12))

      r14xr23=cross_product(r14,r23)
      r13xr12=cross_product(r13,r12)
      r14xr23=r14xr23/sqrt(dot_product(r14xr23,r14xr23))
      r13xr12=r13xr12/sqrt(dot_product(r13xr12,r13xr12))

      intcoo=acos(dot_product(r14xr23,r13xr12))
      intcoo=180.0d0-intcoo*180.0d0/pi
      
      return

    end subroutine calc_pyr

!#######################################################################

    subroutine calc_pyr2(xcoo,intcoo)

      use expec
      use sysdef
      use mathlib

      integer                   :: k1,k2,k3,k4,k5,k6,i
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,val1,val2,pi
      real*8, dimension(3)      :: r14,r23,r13,r12,r41,r56,r46,r45,&
                                   r14xr23,r13xr12,r41xr56,r46xr45
      pi=4.0d0*datan(1.0d0)

      k1=atmindx(1)
      k2=atmindx(2)
      k3=atmindx(3)
      k4=atmindx(4)
      k5=atmindx(5)
      k6=atmindx(6)

      do i=1,3
         r14(i)=xcoo(k4*3-3+i)-xcoo(k1*3-3+i)
         r23(i)=xcoo(k3*3-3+i)-xcoo(k2*3-3+i)
         r13(i)=xcoo(k3*3-3+i)-xcoo(k1*3-3+i)
         r12(i)=xcoo(k2*3-3+i)-xcoo(k1*3-3+i)
         r41(i)=xcoo(k1*3-3+i)-xcoo(k4*3-3+i)
         r56(i)=xcoo(k6*3-3+i)-xcoo(k5*3-3+i)
         r46(i)=xcoo(k6*3-3+i)-xcoo(k4*3-3+i)
         r45(i)=xcoo(k5*3-3+i)-xcoo(k4*3-3+i)
      enddo
      r14=r14/sqrt(dot_product(r14,r14))
      r23=r23/sqrt(dot_product(r23,r23))
      r13=r13/sqrt(dot_product(r13,r13))
      r12=r12/sqrt(dot_product(r12,r12))
      r41=r41/sqrt(dot_product(r41,r41))
      r56=r56/sqrt(dot_product(r56,r56))
      r46=r46/sqrt(dot_product(r46,r46))
      r45=r45/sqrt(dot_product(r45,r45))

      r14xr23=cross_product(r14,r23)
      r13xr12=cross_product(r13,r12)
      r41xr56=cross_product(r41,r56)
      r46xr45=cross_product(r46,r45)
      r14xr23=r14xr23/sqrt(dot_product(r14xr23,r14xr23))
      r13xr12=r13xr12/sqrt(dot_product(r13xr12,r13xr12))
      r41xr56=r41xr56/sqrt(dot_product(r41xr56,r41xr56))
      r46xr45=r46xr45/sqrt(dot_product(r46xr45,r46xr45))

      val1=acos(dot_product(r14xr23,r13xr12))
      val2=acos(dot_product(r41xr56,r46xr45))

      intcoo=max(val1,val2)
      intcoo=180.0d0-intcoo*180.0d0/pi

      return

    end subroutine calc_pyr2

!#######################################################################

    subroutine calc_maxangle(xcoo,intcoo)

      use expec
      use sysdef
      use mathlib

      implicit none

      integer                   :: i,k1,k2,k3
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,dp,len1,len2,pi,ang1,ang2
      real*8, dimension(3)      :: vec1,vec2

!      ! MASSIVE BODGE
!      ! TAKE THE ANGLE OF THE CH2 GROUP IN ETHYLENE FOR WHICH
!      ! THE DEGREE OF PYRAMIDALISATION IS GREATEST
!      integer :: k4,k5,k6
!      real*8, dimension(3)      :: r14,r23,r13,r12,r41,r56,r46,r45,&
!                                   r14xr23,r13xr12,r41xr56,r46xr45
!      real*8 :: val1,val2
!      pi=4.0d0*datan(1.0d0)
!      k1=1
!      k2=3
!      k3=4
!      k4=2
!      k5=5
!      k6=6
!      do i=1,3
!         r14(i)=xcoo(k4*3-3+i)-xcoo(k1*3-3+i)
!         r23(i)=xcoo(k3*3-3+i)-xcoo(k2*3-3+i)
!         r13(i)=xcoo(k3*3-3+i)-xcoo(k1*3-3+i)
!         r12(i)=xcoo(k2*3-3+i)-xcoo(k1*3-3+i)
!         r41(i)=xcoo(k1*3-3+i)-xcoo(k4*3-3+i)
!         r56(i)=xcoo(k6*3-3+i)-xcoo(k5*3-3+i)
!         r46(i)=xcoo(k6*3-3+i)-xcoo(k4*3-3+i)
!         r45(i)=xcoo(k5*3-3+i)-xcoo(k4*3-3+i)
!      enddo
!      r14=r14/sqrt(dot_product(r14,r14))
!      r23=r23/sqrt(dot_product(r23,r23))
!      r13=r13/sqrt(dot_product(r13,r13))
!      r12=r12/sqrt(dot_product(r12,r12))
!      r41=r41/sqrt(dot_product(r41,r41))
!      r56=r56/sqrt(dot_product(r56,r56))
!      r46=r46/sqrt(dot_product(r46,r46))
!      r45=r45/sqrt(dot_product(r45,r45))
!
!      r14xr23=cross_product(r14,r23)
!      r13xr12=cross_product(r13,r12)
!      r41xr56=cross_product(r41,r56)
!      r46xr45=cross_product(r46,r45)
!      r14xr23=r14xr23/sqrt(dot_product(r14xr23,r14xr23))
!      r13xr12=r13xr12/sqrt(dot_product(r13xr12,r13xr12))
!      r41xr56=r41xr56/sqrt(dot_product(r41xr56,r41xr56))
!      r46xr45=r46xr45/sqrt(dot_product(r46xr45,r46xr45))
!
!      val1=acos(dot_product(r14xr23,r13xr12))
!      val2=acos(dot_product(r41xr56,r46xr45))
!      if (val1.gt.val2) then
!         k1=3
!         k2=1
!         k3=4
!      else
!         k1=5
!         k2=2
!         k3=6
!      endif
!      len1=0.0d0
!      len2=0.0d0
!      do i=1,3
!         vec1(i)=xcoo(k2*3-3+i)-xcoo(k1*3-3+i)         
!         vec2(i)=xcoo(k2*3-3+i)-xcoo(k3*3-3+i)
!         len1=len1+vec1(i)**2
!         len2=len2+vec2(i)**2
!      enddo
!      len1=sqrt(len1)
!      len2=sqrt(len2)
!      dp=dot_product(vec1,vec2)
!      ang1=dp/(len1*len2)
!      ang1=dacos(ang1)
!      ang1=ang1*180.0d0/pi
!      intcoo=ang1
!      return
!      ! MASSIVE BODGE
      
      pi=4.0d0*datan(1.0d0)
      
      ! Angle 1
      k1=atmindx(1)
      k2=atmindx(2)
      k3=atmindx(3)

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

      ang1=dp/(len1*len2)
      ang1=dacos(ang1)
      ang1=ang1*180.0d0/pi

      ! Angle 2
      k1=atmindx(4)
      k2=atmindx(5)
      k3=atmindx(6)

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

      ang2=dp/(len1*len2)
      ang2=dacos(ang2)
      ang2=ang2*180.0d0/pi
      
      ! Maximum angle
      intcoo=max(ang1,ang2)

      return

    end subroutine calc_maxangle

!#######################################################################

    subroutine calc_minangle(xcoo,intcoo)

      use expec
      use sysdef

      implicit none

      integer                   :: i,k1,k2,k3
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,dp,len1,len2,pi,ang1,ang2
      real*8, dimension(3)      :: vec1,vec2

      pi=4.0d0*datan(1.0d0)

      ! Angle 1
      k1=atmindx(1)
      k2=atmindx(2)
      k3=atmindx(3)

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

      ang1=dp/(len1*len2)
      ang1=dacos(ang1)
      ang1=ang1*180.0d0/pi

      ! Angle 2
      k1=atmindx(4)
      k2=atmindx(5)
      k3=atmindx(6)

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

      ang2=dp/(len1*len2)
      ang2=dacos(ang2)
      ang2=ang2*180.0d0/pi
      
      ! Maximum angle
      intcoo=min(ang1,ang2)
      
      return

    end subroutine calc_minangle

!#######################################################################

    subroutine calc_diffangle(xcoo,intcoo)

      use expec
      use sysdef

      implicit none

      integer                   :: i,k1,k2,k3
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,dp,len1,len2,pi,ang1,ang2
      real*8, dimension(3)      :: vec1,vec2

      pi=4.0d0*datan(1.0d0)

      ! Angle 1
      k1=atmindx(1)
      k2=atmindx(2)
      k3=atmindx(3)

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

      ang1=dp/(len1*len2)
      ang1=dacos(ang1)
      ang1=ang1*180.0d0/pi

      ! Angle 2
      k1=atmindx(4)
      k2=atmindx(5)
      k3=atmindx(6)

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

      ang2=dp/(len1*len2)
      ang2=dacos(ang2)
      ang2=ang2*180.0d0/pi
      
      ! Absolute angle difference
      intcoo=abs(ang1-ang2)
      
      return

    end subroutine calc_diffangle

!#######################################################################

        subroutine calc_diffpyr(xcoo,intcoo)

      use expec
      use sysdef
      use mathlib

      integer                   :: k1,k2,k3,k4,k5,k6,i
      real*8, dimension(natm*3) :: xcoo
      real*8                    :: intcoo,val1,val2,pi
      real*8, dimension(3)      :: r14,r23,r13,r12,r41,r56,r46,r45,&
                                   r14xr23,r13xr12,r41xr56,r46xr45

      pi=4.0d0*datan(1.0d0)

      k1=atmindx(1)
      k2=atmindx(2)
      k3=atmindx(3)
      k4=atmindx(4)
      k5=atmindx(5)
      k6=atmindx(6)

      do i=1,3
         r14(i)=xcoo(k4*3-3+i)-xcoo(k1*3-3+i)
         r23(i)=xcoo(k3*3-3+i)-xcoo(k2*3-3+i)
         r13(i)=xcoo(k3*3-3+i)-xcoo(k1*3-3+i)
         r12(i)=xcoo(k2*3-3+i)-xcoo(k1*3-3+i)
         r41(i)=xcoo(k1*3-3+i)-xcoo(k4*3-3+i)
         r56(i)=xcoo(k6*3-3+i)-xcoo(k5*3-3+i)
         r46(i)=xcoo(k6*3-3+i)-xcoo(k4*3-3+i)
         r45(i)=xcoo(k5*3-3+i)-xcoo(k4*3-3+i)
      enddo
      r14=r14/sqrt(dot_product(r14,r14))
      r23=r23/sqrt(dot_product(r23,r23))
      r13=r13/sqrt(dot_product(r13,r13))
      r12=r12/sqrt(dot_product(r12,r12))
      r41=r41/sqrt(dot_product(r41,r41))
      r56=r56/sqrt(dot_product(r56,r56))
      r46=r46/sqrt(dot_product(r46,r46))
      r45=r45/sqrt(dot_product(r45,r45))

      r14xr23=cross_product(r14,r23)
      r13xr12=cross_product(r13,r12)
      r41xr56=cross_product(r41,r56)
      r46xr45=cross_product(r46,r45)
      r14xr23=r14xr23/sqrt(dot_product(r14xr23,r14xr23))
      r13xr12=r13xr12/sqrt(dot_product(r13xr12,r13xr12))
      r41xr56=r41xr56/sqrt(dot_product(r41xr56,r41xr56))
      r46xr45=r46xr45/sqrt(dot_product(r46xr45,r46xr45))

      val1=acos(dot_product(r14xr23,r13xr12))
      val2=acos(dot_product(r41xr56,r46xr45))

      ! Absolute difference between the pyramidalisation angles
      intcoo=abs(val1-val2)*180.0d0/pi

      return

    end subroutine calc_diffpyr

!#######################################################################

    subroutine calc_cart(pcoo,xcoo,intmom)

      use expec
      use sysdef

      implicit none

      integer                   :: i,j,k
      real*8, dimension(natm*3) :: pcoo,xcoo,tvec
      real*8                    :: intmom
      real*8, dimension(3,3)    :: rotmat
      
!-----------------------------------------------------------------------
! Exit if we are in the position representation
!-----------------------------------------------------------------------
      if (.not.lmomrep) then
         write(6,'(/,2x,a,/)') 'Reduced densities for Cartesian &
              displacements are currently only available in the &
              momentum representation.'
         STOP
      endif

!-----------------------------------------------------------------------
! Rotate the vector to the same frame as the current geometry
!-----------------------------------------------------------------------
      ! (1) Calculate the transformation matrix rotmat
      call minrmsd(xcoo,refgeom,rotmat)

      ! (2) Rotate the Cartesian vector cartvec
      do i=1,natm
         do j=1,3
            tvec(i*3-3+j)=0.0d0
            do k=1,3
               tvec(i*3-3+j)=tvec(i*3-3+j)+&
                    rotmat(k,j)*cartvec(i*3-3+k)
            enddo
         enddo
      enddo

!-----------------------------------------------------------------------
! Calculate the projection of the momentum onto the vector
!-----------------------------------------------------------------------
      intmom=dot_product(tvec,pcoo)

!      intmom=abs(intmom)

      return

    end subroutine calc_cart

!#######################################################################
! minrmsd: puts the Cartesian coordinates x2 into maximum coincidence
!          with the Cartesian coordinates x1 using the Kabsch algorithm
!#######################################################################

    subroutine minrmsd(x1,x2,rotmat)

      use sysdef

      implicit none

      integer                          :: i,j,k,l
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

      det=finddet2(covar,3)

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

      return

    end subroutine minrmsd

!#######################################################################

    function finddet2(matrix,n)

      implicit none

      integer*4              :: n,i,j,k,l
      real*8, dimension(n,n) :: matrix
      real*8                 :: m,temp,finddet2
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
               finddet2=0
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
      finddet2=l
      do i=1,n
         finddet2=finddet2*matrix(i,i)
      enddo
   
    end function finddet2

!#######################################################################

    subroutine calc_dist_seam(xcoo,intcoo)

      use expec
      use sysdef
      use geomtools
      use mathlib

      real*8, dimension(natm*3) :: xcoo,xcurr,xmin,rotx,dvec
      real*8                    :: intcoo

!-----------------------------------------------------------------------
! Put the input geometry into maximum coincidence with the CI geometry
! subject to any user requested permutations of identical nuclei
!-----------------------------------------------------------------------
      xmin=maxcoinc(cigeom,xcoo)

!-----------------------------------------------------------------------
! Calculate the distance from the minimum RMSD geometry to the CI
! geometry within the branching space
!-----------------------------------------------------------------------
      ! Project the displacement of the minimimum RMSD geometry
      ! from the CI geometry onto the branching space
      dvec=cigeom-xmin
      dvec=matmul(branchproj,dvec)

      ! Set intoo to the length of the projected displacement vector
      intcoo=sqrt(dot_product(dvec,dvec))
      
      return
      
    end subroutine calc_dist_seam

!#######################################################################

  end module intcoo
