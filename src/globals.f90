module globals
use constants


!-------some of the standards used--------!
!3000 + x = lable for all output states   !
!1000 + x = debug lables                  !
!        (x is a positive integer)        !
!-----------------------------------------!

implicit none

public

!Storage for the kohn-sham orbtals
type ks_orbital_type
    !radial eigen vectors and derivatives
   real(dp), dimension(:), pointer :: u
   real(dp), dimension(:), pointer :: up
   real(dp), dimension(:), pointer :: upp
   !other properties
   real(dp) :: e              !kohn sham energy
   integer  :: n              !radial quantum number
   integer  :: m              !angular momentum quantum number
   integer  :: l              !z quantum number
   integer :: lmstate
   integer  :: s              !up = +1 / down = -1
   real(dp) :: f              !occupancy of state (0<=f<=1 for gamma or more for analytic bands)
   real(dp) :: match
end type ks_orbital_type
real(dp), allocatable, dimension(:,:,:), target, public :: vec
integer, allocatable, dimension(:), public :: state_index    !an index of the states in terms of energy

!one type for the orbitals and another for configuration locking
type(ks_orbital_type), public, save, allocatable :: ks(:)


!properties of the radial space
real(dp), allocatable, dimension(:) :: x,x3       !the radial grid
real(dp) :: Xmax                  !upper limit of radial grid
real(dp) :: Xmin                  !lowest point in radial space
integer :: nx                     !number of radial grid points
real(dp) :: alpha                 !x = x_{min} exp(alpha*n)



!for the angular mesh
integer, parameter :: npx=32,nlmx=25
real(dp), dimension(3,npx) :: thephi
real(dp), dimension(npx,nlmx) :: ylm
real(dp), dimension(npx,npx) :: el2
integer, dimension(nlmx) :: llm
real(dp), dimension(3,npx) :: pp




!properties of the atom
real(dp) :: Z                       !nuclear charge
integer :: Ne                       !number of electrons
integer :: Nvirt                    !number of virtuals
integer :: Nstate                   !total number of states
integer :: Nstate_up,Nstate_dn      !total number of states (up and down)
real(dp) :: eft                     !electronic fermi temp
integer :: Nup, Ndn

real(dp) :: sigma
!possible eigenvectors
integer, public :: no_radial_vectors    !count m = 0 states of differing l, n and s
integer, public :: iteration


!control various run parameters
real(dp) :: mix                 !mix parameter for Anderson's of linear
logical :: write_orbs           !true = write orbitals to file
logical :: read_dens_file   
logical :: write_psirays   
integer :: output_level         !levewl of output for efficiency set to zero
integer :: max_iter             !maximum number of self-consistency steps
real(dp) :: tol                 !tolerance for the density deviation


!public kohn-sham potential, Hartree potential and densities
real(dp), allocatable, public :: Vks(:,:)                   !Kohn-Sham potentials without coulomb!!!!
real(dp), allocatable, public :: Vh(:)                       !Hartree
real(dp), allocatable :: n(:,:)                             !(grid point,spin)
real(dp), allocatable :: ntot(:)                             
real(dp), allocatable :: nold(:,:)
real(dp), allocatable :: zeta(:)                            !relative magnetisation



!approximately equal to operator
!for the comparisson of two real numbers
public :: operator(.approx.)

interface operator(.approx.)
   module procedure approx_equal_to
end interface



contains


!define logical approx function for the comparison of two reals
logical function approx_equal_to(x,y)
   real(dp), intent(in) :: x,y
   real(dp), parameter :: real_diff_tol = 1.0e-12_dp

   approx_equal_to = .false.
   if (abs(x-y) <= real_diff_tol) approx_equal_to = .true.

end function approx_equal_to



end module globals
