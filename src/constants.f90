module constants
implicit none

!-------------precision-----------------!
integer, parameter, public :: dp = kind(1.0d0)


!-------------mathematics---------------!
real(dp), parameter, public :: pi = 3.141592653589793238_dp
real(dp), parameter, public :: tpi = 2.0_dp*pi
real(dp), parameter, public :: fpi = 4.0_dp*pi
real(dp), parameter, public :: sqfpi = 3.544907701811031985_dp
real(dp), parameter, public :: e = 2.71828183_dp



!
real(dp), parameter :: zero = 0.0_dp
real(dp), parameter :: one = 1.0_dp

!------physical constants-si units------!
real(dp), parameter, public :: e_si = 1.6022e-19_dp   !electron charge


!----------random number seed-----------!
real(dp), parameter, public :: seed = 57721566.0_dp


!-------hartree atomic units------------!
real(dp), parameter, public :: kb_hartree = 3.1668157e-6 !HK^-1
real(dp), parameter, public :: c_hartree =  137.03599976536860563373  !roughly?!


!-----------thermodynamic units---------!
real(dp), parameter :: kb = 1.0_dp


end module constants
