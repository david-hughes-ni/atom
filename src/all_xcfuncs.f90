module all_xcfuncs
use constants
implicit none


   private

   public :: eps_xc,eps_x_lda,eps_c_perdew_zunger
   public :: mu_xc

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
   real(dp), parameter :: D_p = -0.0048_dp   !not to be confused with dp = kind(1.0d0)
   real(dp), parameter :: gamma_p = -0.0843_dp
   real(dp), parameter :: beta1_p = 1.3981_dp
   real(dp), parameter :: beta2_p = 0.2611_dp


   !useful numbers - pgi compilers don't like expressions in parameters
   real(dp), parameter :: thrd = 1.0_dp/3.0_dp
   !real(dp), parameter :: thrd = 0.3333333333333333_dp
   real(dp), parameter :: fthrd = 4.0_dp/3.0_dp
   !real(dp), parameter :: fthrd = 1.3333333333333333_dp
   real(dp), parameter :: xucoef = -0.75_dp/pi*(9.0/4.0_dp*pi)**thrd
   !real(dp), parameter :: xucoef = -0.4581652932831429548_dp
   real(dp), parameter :: xpcoef = -0.75_dp/pi*(9.0/2.0_dp*pi)**thrd
   !real(dp), parameter :: xpcoef = -0.5772520973386886409_dp
   real(dp), parameter :: vbh_denom = 2.0_dp**fthrd-2.0_dp
   !real(dp), parameter :: vbh_denom = 0.5198420997897463813_dp
   real(dp), parameter :: ddfzz0 = 8.0_dp/vbh_denom/9.0_dp
   !real(dp), parameter :: ddfzz0 = 1.709920934180517292_dp

contains



!-------------------------------------------------------!
!                   MASTER FUNCTION                     !
!   USER SELECTS THE DESIRED XC FUNCTIONALS FROM WHICH  !
!          TO BUILD THE ENERGYIES AND POTENTIALS        !
!-------------------------------------------------------!
function eps_xc(rho,z)
   real(dp) :: eps_xc
   real(dp), intent(in) :: rho,z

   !eps_xc = eps_x_lda(rho,z) + eps_c_perdew_wang(rho,z)
   eps_xc = eps_x_lda(rho,z) + eps_c_perdew_zunger(rho,z)

end function eps_xc


function mu_xc(rho,z,s)
   real(dp) :: mu_xc
   real(dp), intent(in) :: rho,z
   integer, intent(in) :: s

   mu_xc = mu_x_lda(rho,z,s) + mu_c_perdew_wang(rho,z,s)

end function mu_xc




!-------------------------------------------------------!
!           EXCHANGE ENERGY FUNCTIONALS                 !
!-------------------------------------------------------!
!simple local density approximation for the 
!exchange of the homogenous electron gas
function eps_x_lda(rho,zeta)
   real(dp) :: eps_x_lda
   real(dp), intent(in) :: rho,zeta
   real(dp) :: zzeta,rw
   real(dp) :: fxu,fxp,fzz

   eps_x_lda=0.0_dp
   if(rho.gt.1.0e-8_dp) then
      zzeta = zeta
      if (dabs(zeta).gt.1.d0) zzeta=1.0_dp

      rw=(0.75_dp/pi/rho)**thrd
      fxu=xucoef/rw
      fxp=xpcoef/rw

      !Von Barth Hedin interpolation
      fzz=((1.0_dp+zzeta)**fthrd+ &
     &      (1.0_dp-zzeta)**fthrd-2.0_dp)/vbh_denom


      eps_x_lda=fxu+fzz*(fxp-fxu)

   end if

end function eps_x_lda



!perdew and zungers interpolation for the
!correlation energy
function eps_c_perdew_zunger(rho,zeta)
   real(dp) :: eps_c_perdew_zunger
   real(dp), intent(in) :: rho,zeta
   real(dp) :: zzeta,rw
   real(dp) :: fcu,fcp,fzz

   eps_c_perdew_zunger=0.0_dp
   if(rho.gt.1.0e-8_dp) then
      zzeta = zeta
      if(dabs(zeta).gt.1.d0)zzeta=1.0_dp

      rw=(0.75_dp/pi/rho)**thrd

      !Von Barth Hedin interpolation
      fzz=((1.0_dp+zzeta)**fthrd+ &
     &      (1.0_dp-zzeta)**fthrd-2.d0)/vbh_denom

      if(rw.gt.1.0_dp)then
         fcu=gamma_u/(1.0_dp+beta1_u*sqrt(rw)+beta2_u*rw)
         fcp=gamma_p/(1.0_dp+beta1_p*sqrt(rw)+beta2_p*rw)
      else
         fcu=Au*log(rw)+Bu+Cu*rw*log(rw)+Du*rw
         fcp=Ap*log(rw)+Bp+Cp*rw*log(rw)+D_p*rw
      endif

      eps_c_perdew_zunger=fcu+fzz*(fcp-fcu)

   end if

end function eps_c_perdew_zunger


!perdew and zungers interpolation for the
!correlation energy
function eps_c_perdew_wang(rho,zeta)
   real(dp) :: eps_c_perdew_wang
   real(dp), intent(in) :: rho,zeta
   real(dp) :: zzeta,rw
   real(dp) :: fzz,ecp,ecu,alphac

   eps_c_perdew_wang=0.0_dp
   if(rho.gt.1.0e-8_dp) then
      zzeta = zeta
      if(dabs(zeta).gt.1.d0)zzeta=1.0_dp

      rw=(0.75_dp/pi/rho)**thrd

      !Von Barth Hedin interpolation
      fzz=((1.0_dp+zzeta)**fthrd+ &
     &      (1.0_dp-zzeta)**fthrd-2.d0)/vbh_denom

      ecu = Gpw(rw,1)
      ecp = Gpw(rw,2)
      alphac = -Gpw(rw,3)

      eps_c_perdew_wang=ecp+alphac*fzz*(1.0_dp-zzeta**4)/ddfzz0+(ecp-ecu)*fzz*zzeta**4

   end if

end function eps_c_perdew_wang



!-------------------------------------------------------!
!           EXCHANGE POTENTIAL FUNCTIONALS              !
!-------------------------------------------------------!
function mu_x_lda(rho,zeta,isgn)
   real(dp) :: mu_x_lda
   real(dp), intent(in) :: rho,zeta
   integer, intent(in) :: isgn

   real(dp) :: zzeta,rw
   real(dp) :: fxu,fxp
   real(dp) :: fzz,dfzz
   real(dp) :: rsgn

   mu_x_lda=0.0_dp
   if(rho.gt.1.d-08)then

   rw=(0.75_dp/pi/rho)**thrd
   zzeta = zeta
   if(dabs(zeta).gt.1.d0)zzeta=1.0_dp

   fxu=xucoef/rw
   fxp=xpcoef/rw

   !Von Barth Hedin interpolation
   fzz=((1.0_dp+zzeta)**fthrd+ &
     &      (1.0_dp-zzeta)**fthrd-2.d0)/vbh_denom
   dfzz=(1.0_dp+zzeta)**thrd-(1.0_dp-zzeta)**thrd
   dfzz=fthrd*dfzz/vbh_denom


   mu_x_lda = fxu+fzz*(fxp-fxu)
   mu_x_lda = mu_x_lda*fthrd

   rsgn=1.0_dp
   if(isgn.lt.0)rsgn=-1.0_dp

   mu_x_lda=mu_x_lda+rsgn*(1.0_dp-rsgn*zzeta)*(fxp - fxu)*dfzz

   end if

end function mu_x_lda





!-------------------------------------------------------!
!          CORRELATION POTENTIAL FUNCTIONALS            !
!-------------------------------------------------------!
function mu_c_perdew_zunger(rho,zeta,isgn)
   real(dp) :: mu_c_perdew_zunger
   real(dp), intent(in) :: rho,zeta
   integer, intent(in) :: isgn
   real(dp) :: zzeta,rw
   real(dp) :: fcu,fcp,fzz,dfzz
   real(dp) :: dfcp,dfcu
   real(dp) :: excp,excu
   real(dp) :: rsgn

   mu_c_perdew_zunger=0.0_dp
   if(rho.gt.1.0e-8_dp)then

! rw is rs
   rw=(0.75_dp/pi/rho)**thrd
   zzeta = zeta
   if(dabs(zeta).gt.1.d0)zzeta=1.0_dp


   fzz=((1.0_dp+zzeta)**fthrd+ &
      &      (1.0_dp-zzeta)**fthrd-2.0_dp)/vbh_denom
   dfzz=(1.0_dp+zzeta)**thrd-(1.0_dp-zzeta)**thrd
   dfzz=fthrd*dfzz/vbh_denom


   if(rw.gt.1.0_dp)then


      fcu=gamma_u/(1.0_dp+beta1_u*sqrt(rw)+beta2_u*rw)
      fcp=gamma_p/(1.0_dp+beta1_p*sqrt(rw)+beta2_p*rw)

      dfcu=-gamma_u*(0.5_dp*beta1_u/sqrt(rw)+beta2_u)&
            &/(1.0_dp+beta1_u*sqrt(rw)+beta2_u*rw)**2
      dfcp=-gamma_p*(0.5_dp*beta1_p/sqrt(rw)+beta2_p)&
            &/(1.0_dp+beta1_p*sqrt(rw)+beta2_p*rw)**2

      !potcu=fcu*(1.d0+7.d0/6.d0*b1u*dsqrt(rw)+1.33333333d0*b2u*rw)/&
      !&            (1.d0+b1u*dsqrt(rw)+b2u*rw)
      !potcp=fcp*(1.d0+7.d0/6.d0*b1p*dsqrt(rw)+1.33333333d0*b2p*rw)/&
      !&            (1.d0+b1p*dsqrt(rw)+b2p*rw)

   else
      fcu=Au*log(rw)+Bu+Cu*rw*log(rw)+Du*rw
      fcp=Ap*log(rw)+Bp+Cp*rw*log(rw)+D_p*rw

      dfcu=Au/rw+Cu+Cu*log(rw)+Du
      dfcp=Ap/rw+Cp+Cp*log(rw)+D_p

      !potcu=au*dlog(rw)+(bu-au/3.d0)+&
      !&    (2.d0/3.d0)*cu*rw*dlog(rw)+(du+du-cu)*rw/3.d0
      !potcp=ap*dlog(rw)+(bp-ap/3.d0)+&
      !&    (2.d0/3.d0)*cp*rw*dlog(rw)+(ddp+ddp-cp)*rw/3.d0
   endif


   mu_c_perdew_zunger=fcu+(fcp-fcu)*fzz-thrd/rw*( dfcu+(dfcp-dfcu)*fzz )

   rsgn=1.0_dp
   if(isgn.lt.0)rsgn=-1.0_dp

   mu_c_perdew_zunger=mu_c_perdew_zunger+rsgn*(1.0_dp-rsgn*zzeta)*(fcp-fcu)*dfzz

   end if

end function mu_c_perdew_zunger




!potential due to perdew and wang
function mu_c_perdew_wang(rho,zeta,isgn)
   real(dp) :: mu_c_perdew_wang
   real(dp), intent(in) :: rho,zeta
   integer, intent(in) :: isgn

   real(dp) :: zzeta,rw
   real(dp) :: ders,dez,fzz,dfzz
   real(dp) :: ecu,ecp,alphac

   mu_c_perdew_wang=0.0_dp
   if(rho.gt.1.0e-8_dp)then

   rw=(0.75_dp/pi/rho)**thrd
   zzeta = zeta
   if(dabs(zeta).gt.1.d0)zzeta=1.0_dp

   !Von Barth Hedin interpolation
   fzz=((1.0_dp+zzeta)**fthrd+ &
     &      (1.0_dp-zzeta)**fthrd-2.d0)/vbh_denom
   dfzz=(1.0_dp+zzeta)**thrd-(1.0_dp-zzeta)**thrd
   dfzz=fthrd*dfzz/vbh_denom

   ecu = Gpw(rw,1)
   ecp = Gpw(rw,2)
   alphac = Gpw(rw,3)

   mu_c_perdew_wang = eps_c_perdew_wang(rho,zeta)

   ders = dGpw(rw,1)*(1.0_dp-fzz*zzeta**4)
   ders = ders + dGpw(rw,2)*fzz*zzeta**4
   ders = ders - dGpw(rw,3)*fzz*(1.0_dp-zzeta**4)/ddfzz0

   dez = 4.0_dp*zzeta**3*(ecp-ecu-alphac/ddfzz0)
   dez = dez+dfzz*(zzeta**4*ecp-zzeta**4*ecu+(1.0_dp-zzeta**4)*alphac/ddfzz0)

   mu_c_perdew_wang = mu_c_perdew_wang - thrd*rw*ders -(zzeta-real(isgn))*dez

   end if

end function mu_c_perdew_wang


!---------------------------------------!
!      INTERPOLATION FUNCTIONALS        !
!---------------------------------------!
!the Perdew and Wang interpolation functional
function Gpw(rw,fit)
   real(dp) :: Gpw
   real(dp), intent(in) :: rw
   integer, intent(in) :: fit
   real(dp) :: A_pw,alpha1_pw,beta1_pw
   real(dp) :: beta2_pw,beta3_pw,beta4_pw,p_pw
   real(dp) :: arg

   select case(fit)
   case(1)   !eps_c unpolarised
      p_pw = 1.0_dp
      A_pw = 0.031091_dp
      alpha1_pw = 0.21370_dp
      beta1_pw = 7.5957_dp
      beta2_pw = 3.5876_dp
      beta3_pw = 1.6382_dp
      beta4_pw = 0.49294_dp
   case(2)   !eps_c polarised
      p_pw = 1.0_dp
      A_pw = 0.015545_dp
      alpha1_pw = 0.20548_dp
      beta1_pw = 14.1189_dp
      beta2_pw = 6.1977_dp
      beta3_pw = 3.3662_dp
      beta4_pw = 0.62517_dp
   case(3)   !alpha_c
      p_pw = 1.0_dp
      A_pw = 0.016887_dp
      alpha1_pw = 0.11125_dp
      beta1_pw = 10.357_dp
      beta2_pw = 3.6231_dp
      beta3_pw = 0.88026_dp
      beta4_pw = 0.49671_dp
   case default    !error
      write(*,*) 'invalid Perdew Wang interpolation option'
      stop
   end select

   arg = 1.0_dp+1.0_dp/2.0_dp/A_pw/&
         &(beta1_pw*sqrt(rw)+beta2_pw*rw+beta3_pw*rw**1.5_dp+beta4_pw*rw**(p_pw+1.0_dp))
   Gpw = -2.0_dp*A_pw*(1.0_dp+alpha1_pw*rw)*log(arg)

end function gpw


function dGpw(rw,fit)
   real(dp) :: dGpw
   real(dp), intent(in) :: rw
   integer, intent(in) :: fit
   real(dp) :: A_pw,alpha1_pw,beta1_pw
   real(dp) :: beta2_pw,beta3_pw,beta4_pw,p_pw
   real(dp) :: arg,Q0,Q1,Q1p

   select case(fit)
   case(1)   !eps_c unpolarised
      p_pw = 1.0_dp
      A_pw = 0.031091_dp
      alpha1_pw = 0.21370_dp
      beta1_pw = 7.5957_dp
      beta2_pw = 3.5876_dp
      beta3_pw = 1.6382_dp
      beta4_pw = 0.49294_dp
   case(2)   !eps_c polarised
      p_pw = 1.0_dp
      A_pw = 0.015545_dp
      alpha1_pw = 0.20548_dp
      beta1_pw = 14.1189_dp
      beta2_pw = 6.1977_dp
      beta3_pw = 3.3662_dp
      beta4_pw = 0.62517_dp
   case(3)   !alpha_c
      p_pw = 1.0_dp
      A_pw = 0.016887_dp
      alpha1_pw = 0.11125_dp
      beta1_pw = 10.357_dp
      beta2_pw = 3.6231_dp
      beta3_pw = 0.88026_dp
      beta4_pw = 0.49671_dp
   case default   !error
      write(*,*) 'invalid Perdew Wang interpolation option'
      stop
   end select

   Q0 = -2.0_dp*A_pw*(1.0_dp+alpha1_pw*rw)
   Q1 = 2.0_dp*A_pw*(beta1_pw*sqrt(rw)+beta2_pw*rw+&
                   &beta3_pw*rw**1.5_dp+beta4_pw*rw**(p_pw+1.0_dp))
   Q1p = A_pw*(beta1_pw/sqrt(rw)+2.0_dp*beta2_pw+&
                   &3.0_dp*beta3_pw*sqrt(rw)+2.0_dp*(p_pw+1.0_dp)*beta4_pw*rw**p_pw)

   dGpw = -2.0_dp*A_pw*alpha1_pw*log(1.0_dp+1.0_dp/Q1) - Q0*Q1p/(Q1**2+Q1)

end function dGpw



end module all_xcfuncs
