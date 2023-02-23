!------------------------------------------------!
!    The exchange and correlational functionals  !
!------------------------------------------------!
module xcfuncs

use constants
use globals

implicit none

private


   !public functionals
   public :: mu_xc  !(i,A,s)  grid point, amplitude, spin
   public :: eps_xc !(i,A)    grid point, amplitude
   public :: eps_x
   public :: eps_c,eps_cw


   !non-spin polarized constants according to 
   !Gell-Mann, Brueckner (1957) Au,Bu
   !Ceperley, Alder (1980) Cu,Du,\gamma^u,\beta_1^u,\beta_2^u

   !unpolarised constants
   real(dp), parameter :: Au = 0.0311_dp
   real(dp), parameter :: Bu = -0.048_dp
   real(dp), parameter :: Cu = 0.002_dp
   real(dp), parameter :: Du = -0.0116_dp
   real(dp), parameter :: gamma_u = -0.1423_dp
   real(dp), parameter :: beta1_u = 1.0529_dp
   real(dp), parameter :: beta2_u = 0.3334_dp

   !polarized constants
   real(dp), parameter :: Ap = 0.01555_dp
   real(dp), parameter :: Bp = -0.0269_dp
   real(dp), parameter :: Cp = 0.0007_dp
   real(dp), parameter :: D_p = -0.0048_dp
   real(dp), parameter :: gamma_p = -0.0843_dp
   real(dp), parameter :: beta1_p = 1.3981_dp
   real(dp), parameter :: beta2_p = 0.2611_dp


   !useful numbers
   real(dp), parameter :: thrd = 1.0_dp/3.0_dp
   real(dp), parameter :: fthrd = 4.0_dp/3.0_dp
   real(dp), parameter :: xucoef = -0.75_dp/pi*(9.0/4.0_dp*pi)**thrd
   real(dp), parameter :: xpcoef = -0.75_dp/pi*(9.0/2.0_dp*pi)**thrd


contains


   !-----------------------------------------------------------------------------!
   !                               XC Energy density                             !
   !-----------------------------------------------------------------------------!

   real(dp) function eps_xc(rho,z)
      real(dp), intent(in) :: rho,z

      eps_xc = eps_x(rho,z) + eps_c(rho,z)

   end function eps_xc



   real(dp) function eps_x(rho,z)
      real(dp), intent(in) :: rho,z
      real(dp) :: local_rs

      local_rs = (3.0_dp/4.0_dp/pi/rho)**thrd
      !real(dp), intent(in) :: A

      eps_x = epsxu(rho) + (epsxp(rho)-epsxu(rho))*fvbh(z)

   end function eps_x

   real(dp) function eps_c(rho,z)
      real(dp), intent(in) :: rho,z

      eps_c = epsc_pz(rho,.false.) + (epsc_pz(rho,.true.)-epsc_pz(rho,.false.))*fvbh(z)

   end function eps_c


   real(dp) function epsxp(rho)
      real(dp), intent(in) :: rho
      real(dp) :: local_rs

      local_rs = (3.0_dp/4.0_dp/pi/rho)**thrd

      epsxp = xpcoef/local_rs

   end function epsxp

   real(dp) function epsxu(rho)
      real(dp), intent(in) :: rho
      real(dp) :: local_rs

      local_rs = (3.0_dp/4.0_dp/pi/rho)**thrd

      epsxu = xucoef/local_rs

   end function epsxu


   !the correlation, well at least according to Wigner
   real(dp) function eps_cw(rho)
      real(dp), intent(in) :: rho
      real(dp) :: local_rs
      local_rs = (3.0_dp/4.0_dp/pi/rho)**thrd
      eps_cw = -0.44_dp/(local_rs+7.8_dp)
   end function eps_cw



   !Von Barth & Hedin interpolation function
   !J. Phys. C 5, 1629-1642 (1972) - a local exchange-correlation potential for
   !the spin polarised case
   real(dp) function fvbh(z)
       real(dp), intent(in) :: z
       real(dp) ::  denom

       denom = 2.0_dp**fthrd - 2.0_dp

       fvbh = ((1.0_dp+z)**fthrd+(1.0_dp-z)**fthrd-2.0_dp)/denom

   end function fvbh


   !derivative
   real(dp) function fvbh_prime(z)
       real(dp), intent(in) :: z
       real(dp) :: denom

       denom = 2.0_dp**fthrd - 2.0_dp

       fvbh_prime = fthrd*((1.0_dp+z)**thrd-(1.0_dp-z)**thrd)/denom


   end function fvbh_prime

 
   !------------------------------------------------------------------------------!
   !                                XC Potential                                  !
   !------------------------------------------------------------------------------!

   real(dp) function mu_xc(rho,z,s)
      real(dp), intent(in) :: rho,z
      integer, intent(in) :: s

      mu_xc = mu_x(rho,z,s) + mu_c(rho,z,s)

   end function mu_xc

   real(dp) function mu_x(rho,z,s)
      real(dp), intent(in) :: rho,z
      integer, intent(in) :: s

      real(dp) :: local_rs

      local_rs = (3.0_dp/(4.0_dp*pi*rho))**thrd

      if (s == 1) then
         mu_x = fthrd*eps_x(rho,z) +(1.0_dp-z)*(epsxp(rho)-epsxu(rho))*fvbh_prime(z)
      else if (s == -1) then
         mu_x = fthrd*eps_x(rho,z) -(1.0_dp+z)*(epsxp(rho)-epsxu(rho))*fvbh_prime(z)
      else
         write(*,*) "not a valid spin"
         stop
      end if

   end function mu_x



   real(dp) function mu_c(rho,z,s)
      real(dp), intent(in) :: rho,z
      integer, intent(in) :: s

      real(dp) :: local_rs

      local_rs = (3.0_dp/(4.0_dp*pi*rho))**thrd

      if (s == 1) then
         mu_c = eps_c(rho,z)-local_rs*depsc_drs(rho,z)*thrd+(1.0_dp-z)*depsc_dz(rho,z)
      else if (s == -1) then
         mu_c = eps_c(rho,z)-local_rs*depsc_drs(rho,z)*thrd-(1.0_dp+z)*depsc_dz(rho,z)
      else
         write(*,*) "not a valid spin"
         stop
      end if

   end function mu_c



   real(dp) function depsc_drs(rho,z)
      real(dp), intent(in) :: rho,z

      depsc_drs = epsc_pz_prime(rho,.false.) + &
                        &(epsc_pz_prime(rho,.true.)-epsc_pz_prime(rho,.false.))*fvbh(z)

   end function depsc_drs

   real(dp) function depsc_dz(rho,z)
      real(dp), intent(in) :: rho,z

      depsc_dz = (epsc_pz(rho,.true.)-epsc_pz(rho,.false.))*fvbh_prime(z)

   end function depsc_dz

   !---------------------------------!
   !  construct the correlation      !
   !  energy density according to    !
   !    Perdew and Zunger (1981)     !
   !---------------------------------!
   real(dp) function epsC_PZ(rho,p)
      real(dp), intent(in) :: rho
      logical, intent(in) :: p

      real(dp) :: A,B,C,D,gamma,beta1,beta2
      real(dp) :: local_rs    !wigner radius


      !non-spin polarized constants according to
      !Gell-Mann, Brueckner (1957) Au,Bu
      !Ceperley, Alder (1980) Cu,Du,\gamma^u,\beta_1^u,\beta_2^u

      !polarisation switch
      if (p) then
         A = Ap
         B = Bp
         C = Cp
         D = D_p
         gamma = gamma_p
         beta1 = beta1_p
         beta2 = beta2_p
      else
         A = Au
         B = Bu
         C = Cu
         D = Du
         gamma = gamma_u
         beta1 = beta1_u
         beta2 = beta2_u
      end if

      local_rs = (3.0_dp/(4.0_dp*pi*rho))**thrd
      !setting up the potential
      if (local_rs <= 1.0_dp) then
        epsC_PZ = A*log(local_rs)+ B+C*local_rs*log(local_rs)+D*local_rs
      else 
        epsC_PZ = gamma/(1.0_dp+beta1*sqrt(local_rs)+beta2*local_rs)
      end if

   end function epsC_PZ


   !---------------------------------!
   !     correlation functional      !
   !    Perdew and Zunger (1981)     !
   !---------------------------------!
   real(dp) function epsC_PZ_prime(rho,p)
      real(dp), intent(in) :: rho
      logical, intent(in) :: p

      real(dp) :: local_rs    !mean inter electronic distance
      real(dp) :: A,B,C,D,gamma,beta1,beta2



      if (p) then
         A = Ap
         B = Bp
         C = Cp
         D = D_p
         gamma = gamma_p
         beta1 = beta1_p
         beta2 = beta2_p
      else
         A = Au
         B = Bu
         C = Cu
         D = Du
         gamma = gamma_u
         beta1 = beta1_u
         beta2 = beta2_u
      end if

      local_rs = (3.0_dp/(4.0_dp*pi*rho))**thrd
      !setting up the potential
      if (local_rs <= 1.0_dp) then
        epsC_PZ_prime = A/local_rs+C*(log(local_rs)+1.0_dp)+D
      else
        epsC_PZ_prime = -gamma*(beta1/(2.0_dp*sqrt(local_rs))+beta2)/&
               &(1.0_dp+beta1*sqrt(local_rs)+beta2*local_rs)**2
      end if

   end function epsC_PZ_prime



end module xcfuncs
