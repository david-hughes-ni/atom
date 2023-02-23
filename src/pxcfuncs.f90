module pxcfuncs
use constants
implicit none




!   real(dp), parameter :: cfzz=2.d0**(4.d0/3.d0)-2.d0
real(dp), parameter :: cfzz = 0.519842099783923928_dp
!c
   real(dp), parameter :: au=0.0311d0
   real(dp), parameter :: bu=-0.048d0
!c Brueckner Gell-Mann
   real(dp), parameter :: cu=0.0020d0
   real(dp), parameter :: du=-0.0116d0
!c End of Brueckner Gell-Mann
!c
   real(dp), parameter :: gu=-0.1423d0
   real(dp), parameter :: b1u=1.0529d0
   real(dp), parameter :: b2u=0.3334d0
!c
   real(dp), parameter :: ap=0.01555
   real(dp), parameter :: bp=-0.0269d0
!c Brueckner Gell-Mann
   real(dp), parameter :: cp=0.0007d0
   real(dp), parameter :: ddp=-0.0048d0
!c End of Brueckner Gell-Mann
!c
   real(dp), parameter :: gp=-0.0843d0
   real(dp), parameter :: b1p=1.3981d0
   real(dp), parameter :: b2p=0.2611d0


   contains



   function eps_xc(rho,zeta)
      implicit none
      real(dp) :: eps_xc
      real(dp), intent(in) :: rho,zeta
      real(dp) :: zzeta,rw
      real(dp) :: fxu,fxp,fcu,fcp,fzz
      real(dp) :: excp,excu

      eps_xc=0.d0
      if(rho.gt.1.d-08) then
      zzeta = zeta
      if(dabs(zeta).gt.1.d0)zzeta=zeta/(dabs(zeta)+1.d-10)

      rw=(3.d0/fpi/rho)**(1.d0/3.d0)
      fxu=-0.458d0/rw
      fxp=-0.458d0*(2.d0**(1.d0/3.D0))/rw

      fzz=((1.d0+zzeta)**(4.d0/3.d0)+ &
     &      (1.d0-zzeta)**(4.d0/3.d0)-2.d0)/cfzz

      if(rw.gt.1.d0)then
         fcu=gu/(1.d0+b1u*dsqrt(rw)+b2u*rw)
         fcp=gp/(1.d0+b1p*dsqrt(rw)+b2p*rw)
      else
         fcu=au*dlog(rw)+bu+cu*rw*dlog(rw)+du*rw
         fcp=ap*dlog(rw)+bp+cp*rw*dlog(rw)+ddp*rw
      endif

      eps_xc=(fxu+fcu)+fzz*(fxp+fcp-fxu-fcu)

      end if

   end function eps_xc



   function mu_xc(rho,zeta,isgn)
      implicit none

      real(dp) :: mu_xc
      real(dp), intent(in) :: rho,zeta
      integer, intent(in) :: isgn
      real(dp) :: zzeta,rw
      real(dp) :: fxu,fxp,fcu,fcp,fzz
      real(dp) :: dfzz,potxu,potxp
      real(dp) :: potcu,potcp
      real(dp) :: potxcp,potxcu
      real(dp) :: excp,excu
      integer :: isg

      mu_xc=0.d0
      if(rho.lt.1.d-08)return
! rw is rs
      rw=(3.d0/fpi/rho)**(1.d0/3.d0)
      zzeta = zeta
      if(dabs(zeta).gt.1.d0)zzeta=zeta/(dabs(zeta)+1.d-10)

      fxu=-0.458d0/rw
      fxp=fxu*(2.d0**(1.d0/3.D0))

      fzz=((1.d0+zeta)**(4.d0/3.d0)+ &
     &      (1.d0-zeta)**(4.d0/3.d0)-2.d0)/cfzz
      dfzz=(1.d0+zeta)**(1.d0/3.d0)-(1.d0-zeta)**(1.d0/3.d0)
      dfzz=4.d0/3.d0*dfzz/cfzz

      potxu=-0.458d0/rw*(1.3333333333d0)
      potxp= potxu*(2.d0**(1.d0/3.d0))

      if(rw.gt.1.d0)then

         fcu=gu/(1.d0+b1u*dsqrt(rw)+b2u*rw)
         fcp=gp/(1.d0+b1p*dsqrt(rw)+b2p*rw)

         potcu=fcu*(1.d0+7.d0/6.d0*b1u*dsqrt(rw)+1.33333333d0*b2u*rw)/&
     &            (1.d0+b1u*dsqrt(rw)+b2u*rw)
         potcp=fcp*(1.d0+7.d0/6.d0*b1p*dsqrt(rw)+1.33333333d0*b2p*rw)/&
     &            (1.d0+b1p*dsqrt(rw)+b2p*rw)

      else
         fcu=au*dlog(rw)+bu+cu*rw*dlog(rw)+du*rw
         fcp=ap*dlog(rw)+bp+cp*rw*dlog(rw)+ddp*rw

         potcu=au*dlog(rw)+(bu-au/3.d0)+&
     &    (2.d0/3.d0)*cu*rw*dlog(rw)+(du+du-cu)*rw/3.d0
         potcp=ap*dlog(rw)+(bp-ap/3.d0)+&
     &    (2.d0/3.d0)*cp*rw*dlog(rw)+(ddp+ddp-cp)*rw/3.d0
      endif

      excu=fxu+fcu
      excp=fxp+fcp

      potxcu=potxu+potcu
      potxcp=potxp+potcp
      
      isg=1
      if(isgn.lt.0)isg=-1

      mu_xc=potxcu+fzz*(potxcp-potxcu)+(excp-excu)*(real(isg,dp)-zeta)*dfzz

   end function mu_xc



end module pxcfuncs
