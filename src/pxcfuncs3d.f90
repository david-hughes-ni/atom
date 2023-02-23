module pxcfuncs3d
use constants
implicit none






   contains



   function eps_xc(rho,zeta)
      implicit none
      real(dp) :: eps_xc
      real(dp), intent(in) :: rho,zeta
      real(dp) :: rs,sqrs,den,excu,excp,rsl,fz


        if(rho .le. 1.d-20) then
          eps_xc=0.d0
        else

        rs=0.62035049d0*rho**(-0.3333333333333333d0)
        if(rho .lt. 0.23873241d0) then
          sqrs=dsqrt(rs)
          den=1.0d0 + 1.0529d0*sqrs + 0.3334d0*rs
          excu=-0.4582d0/rs - 0.1423d0/den
          excp=-0.577252d0/rs - 0.0843d0/den
        else
          rsl=dlog(rs)
          excu=-0.4582d0/rs - 0.0480d0 + 0.0311d0*rsl - 0.0116d0*rs&
     &      + 0.002d0*rs*rsl
          excp=-0.577252d0/rs - 0.0269d0+0.01555d0*rsl-0.0048d0*rs&
     &      + 0.0007d0*rs*rsl
        end if
        fz=(1.d0+zeta)**(4.d0/3.d0)+(1.d0-zeta)**(4.d0/3.d0)-2.d0
        fz=fz/(2.d0**(4.d0/3.d0)-2.d0)
        eps_xc=excu+fz*(excp-excu)

        end if

   end function eps_xc



   function mu_xc(rho,zeta,isgn)
      implicit none

      real(dp) :: mu_xc
      real(dp), intent(in) :: rho,zeta
      integer, intent(in) :: isgn

      real(dp) :: rs,sqrs,excu,den,fxcu,excp,fxcp,rsl
      real(dp) :: dfdz,fz,sgns


        if(rho .le. 1.d-20) then
          mu_xc=0.d0
        else

        rs=0.62035049d0*rho**(-0.3333333333333333d0)
        if(rho .lt. 0.23873241d0) then
          sqrs=dsqrt(rs)
          den=1.0d0 + 1.0529d0*sqrs + 0.3334d0*rs
          excu=-0.4582d0/rs - 0.1423d0/den
          fxcu=excu - rs*(0.15273333d0/rs**2&
     &      + (0.02497128d0/sqrs + 0.01581427d0)/den**2)
          den=1.0d0 + 1.3981d0*sqrs + 0.2611d0*rs
          excp=-0.577252d0/rs - 0.0843d0/den
          fxcp=excp - rs*(0.19241736d0/rs**2&
     &      + (0.01964331d0/sqrs + 0.00733691d0)/den**2)
        else
          rsl=dlog(rs)
          excu=-0.4582d0/rs - 0.0480d0 + 0.0311d0*rsl - 0.0116d0*rs&
     &      + 0.002d0*rs*rsl
          fxcu=excu - rs*(0.15273333d0/rs**2&
     &      + 0.01036667d0/rs - 0.003866667d0&
     &      + 0.00066667d0*(1.0d0 + rsl))
          excp=-0.577252d0/rs - 0.0269d0+0.01555d0*rsl-0.0048d0*rs&
     &      + 0.0007d0*rs*rsl
          fxcp=excp - rs*(0.19241736d0/rs**2&
     &      + 0.00518333d0/rs - 0.0016d0&
     &      + 0.00023333d0*(1.0d0 + rsl))
        end if
        sgns=1.d0
        if(isgn.eq.-1)sgns=-1.d0
        fz=(1.d0+zeta)**(1.333333333333333d0)+&
     &     (1.d0-zeta)**(1.333333333333333d0)-2.d0
        fz=fz*1.92366105093153617d0
        dfdz=2.56488140124204822d0*&
     &         ((1.d0+zeta)**(0.3333333333333d0)&
     &         -(1.d0-zeta)**(0.3333333333333d0))
        mu_xc=fxcu+fz*(fxcp-fxcu)+&
     &     (excp-excu)*(sgns-zeta)*dfdz

        end if

   end function mu_xc



end module pxcfuncs3d
