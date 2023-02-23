module dfsubs

use constants
use globals
!use all_xcfuncs
use pxcfuncs3d
use output

implicit none

private


   public :: start_density
   public :: update_density
   public :: build_potential
   public :: Energy
   public :: Energy2
   public :: testhartree

contains






   !---------------------------------!
   !    creates initial density      !
   ! For now just a simple gaussian  !
   !  Soon it will be Thomas Fermi   !
   !---------------------------------!
   subroutine start_density
   implicit none

      integer :: i
      real(dp) :: norm
      real(dp) :: mag,maxdens
      real(dp) :: nup_old,ndown_old



      !making the electron density
      if ( read_dens_file ) then

         open(unit=3,file="data/dens.dat",status="old",action="read")
         do i = 1,Nx
            read(unit=3,fmt=*) norm,n(i,1),n(i,2)
         end do
         close(3)

         ntot(:) = n(:,1) + n(:,2)

      else


        sigma = (20.0_dp/xmax)**2
        !make an exponential guess
         do i = 1,Nx
           ntot(i) = exp(-sigma*x(i))
         end do


         !normalisation
         norm = 0.5_dp*(ntot(1)*x3(1) + ntot(Nx)*x3(Nx))
         do i = 2,Nx-1
            norm = norm + ntot(i)*x3(i)
         end do
         norm = norm*alpha*fpi
         norm = 1.0_dp/norm
         do i = 1,Nx
           ntot(i) = ntot(i)*norm*real(Ne,dp)
           n(i,1) = ntot(i)*real(Nup,dp)/real(Ne,dp)
           n(i,2) = ntot(i)*real(Ndn,dp)/real(Ne,dp)
         end do

!        check the analytic expression for the gaussian normalisation
!         norm = 0.0_dp
!         do i = 1,Nx
!            norm = norm + x3(i)*real(Ne,dp)*sqrt(sigma**3/pi)*exp(-sigma*x(i)**2)*4.0_dp
!         end do
!         norm = norm*alpha
!         write(*,*) norm


         !initial density to file
         open(unit=10,file="data/start_density.dat",status="replace")
!         open(unit=11,file="data/analytic_density.dat",status="replace")
         do i = 1,Nx
            write(unit=10,fmt='(4f22.10)') x(i),n(i,1:2),ntot(i)
!            write(unit=11,fmt='(2f22.10)') x(i),real(Ne,dp)*sqrt(sigma**3/pi)*exp(-sigma*x(i)**2)*4.0_dp
         end do
         close(10)
!         close(11)
!stop

      end if


      !total polarisation density
      maxdens = maxval(ntot)
      do i = 1,Nx

         mag = abs(n(i,1) - n(i,2))

         if (ntot(i)/maxdens < 1.0e-8_dp) then
            !very low density => no polarization
            zeta(i) = 0.0_dp
         else if (mag >= 0.0_dp) then
            !never greater than 1!
            zeta(i) = min(mag/ntot(i),1.0_dp)
         end if

      end do


      !electronic charge
      norm = 0.5_dp*( ntot(1)*x3(1) + ntot(Nx)*x3(Nx) )
      do i = 2,Nx-1
         norm = norm + ntot(i)*x3(i)
      end do
      norm = -norm*alpha*fpi
      write(*,'(a,1x,f18.12)') "Initial Electronic Charge:",norm


      !spin up electronic
      norm = 0.5_dp*( n(1,1)*x3(1)+n(Nx,1)*x3(Nx) )
      do i = 2,Nx-1
         norm = norm + n(i,1)*x3(i)
      end do
      norm = -norm*alpha*fpi
      nup_old=-norm
      write(*,'(a,1x,f18.12)') "Initial 'spin up' Charge:",norm


      !spin down electron charge
      norm = 0.5_dp*( n(1,2)*x3(1)+n(Nx,2)*x3(Nx) )
      do i = 2,Nx-1
         norm = norm + n(i,2)*x3(i)
      end do
      norm = -norm*alpha*fpi
      write(*,'(a,1x,f18.12)') "Initial 'spin dn' Charge:",norm
      ndown_old = -norm




      write(*,'(a)') "--------------------------------"
      write(*,*)
      write(*,*)


   end subroutine start_density



  !--------------------------------!
  !       Make the new density     !
  !  from the Kohn-Sham orbitals   !
  !      and mix with the old      !
  !--------------------------------!
  subroutine update_density( equal ) 
    implicit none

    logical, intent(out) :: equal
    real(dp) :: ss,mag,maxdens
    integer :: i,ix,ispin


    !wipe
    nold = n
    n = 0.0_dp
    ntot = 0.0_dp


    !build new density
    do i = 1,Nstate

      if ( ks(i)%s > 0) then

        do ix = 1,Nx
          n(ix,1) = n(ix,1) + ( ks(i)%u(ix)/x(ix) )**2*ks(i)%f
        end do

      else

        do ix = 1,Nx
          n(ix,2) = n(ix,2) + ( ks(i)%u(ix)/x(ix) )**2*ks(i)%f
        end do

      end if

    end do

    !test for equality - add later
    equal = .false.
    ss = 0.0_dp
    do ix = 1,Nx
      ss = ss + x(ix)*( n(ix,1) - nold(ix,1) )
    end do
    ss = ss*alpha*fpi
    write(*,*) "DD:",ss
    if ( abs(ss) <= tol )  equal = .true.



    !linear mixing
    do ispin = 1,2
      do ix = 1,Nx
        n(ix,ispin) = (1.0_dp-mix)*nold(ix,ispin) + mix*n(ix,ispin)
      end do
    end do


    !total number density
    do ix = 1,Nx
      ntot(ix) = n(ix,1) + n(ix,2)
    end do


    !total polarisation density
    maxdens = maxval(ntot)
    do ix = 1,Nx

      mag = abs(n(ix,1) - n(ix,2))

      if (ntot(ix)/maxdens < 1.0e-8_dp) then
        !very low density => no polarization
        zeta(ix) = 0.0_dp
      else if (mag >= 0.0_dp) then
        !never greater than 1!
        zeta(ix) = min(mag/ntot(ix),1.0_dp)
      end if

    end do


!check charge conservation - working
!    ss = 0.0_dp
!    do ix = 1,Nx
!      ss = ss + ntot(ix)*x3(ix)
!      write(300,*) x(ix),ntot(ix)
!    end do
!    write(*,*) "update charge:",ss*fpi*alpha
!    stop


  end subroutine update_density



   !---------------------------------!
   !    creates the updated Kohn-    !
   !        Sham potential           !
   !   No Coulomb term included      !
   !---------------------------------!
   subroutine build_potential

     integer :: i


     !a shooting routine
     call hartree


     !effective potential
     do i = 1,Nx
        VKS(i,1) = Vh(i) + mu_xc(ntot(i),zeta(i),1)
        VKS(i,2) = Vh(i) + mu_xc(ntot(i),zeta(i),-1)
!        if (x(i)>0.8_dp*Xmax) then
!          VKS(i,1) = VKS(i,1) + 20.0_dp*(x(i) - 0.8_dp*Xmax)**2
!          VKS(i,2) = VKS(i,2) + 20.0_dp*(x(i) - 0.8_dp*Xmax)**2
!        end if
     end do


     !for the user
!     if ( output_level > 0 ) then

       open(10,file='data/vks.dat',status="replace")
       open(20,file='data/vxc.dat',status="replace")
       do i = 1,Nx
         write(10,*) x(i),VKS(i,1:2)
         write(20,*) x(i),mu_xc(ntot(i),zeta(i),1),mu_xc(ntot(i),zeta(i),-1)
       end do
       close(10)
       close(20)
!     end if

   end subroutine build_potential




  !make the Hartree potential
  Subroutine hartree
    implicit none

    real(dp) :: rmome
    real(dp), dimension(Nx) :: phi,w,wp
    real(dp) :: const,rho0
    real(dp) :: v0,v2,wm
    real(dp) :: a,b,c,dd
    real(dp) :: rr,php,phppi,xx
    character(3) :: char_iter

    integer :: ir,mp


     do ir = 1,nx
       phi(ir) = 0.d0
       w(ir) = 0.d0
       wp(ir) = 0.d0
     end do


!    short range expansion
     rho0 = ntot(1)
     v2   = -pi*rho0
     do ir = 1,4
       phi(ir)   = x(ir)**2*( v2 )
       w(ir)  = 2.0_dp*alpha*v2*x(ir)**2
       wp(ir) = -4.d0*pi*(alpha*x(ir))**2*ntot(ir) - alpha*w(ir)
     end do

!    long-range expansion
     do ir = nx-4,nx
       phi(ir)   =  real(Ne,dp)/x(ir) 
       w(ir)  =  -alpha*real(Ne,dp)/x(ir)
       wp(ir) = -4.d0*pi*(alpha*x(ir))**2*ntot(ir) - alpha*w(ir)
     end do


!    outwards shoot to mp
     mp = Nx*0.85_dp
!     write(*,*) mp,x(mp)
     do ir = 5,mp-1

!      predict
       a  = -4.d0*pi*(alpha*x(ir-1))**2*ntot(ir-1) - alpha*w(ir-1)
       b  = -4.d0*pi*(alpha*x(ir-2))**2*ntot(ir-2) - alpha*w(ir-2)
       c  = -4.d0*pi*(alpha*x(ir-3))**2*ntot(ir-3) - alpha*w(ir-3)
       dd = -4.d0*pi*(alpha*x(ir-4))**2*ntot(ir-4) - alpha*w(ir-4)
       w(ir) = w(ir-1) + (55.d0*a-59.d0*b+37.d0*c-9.d0*dd)/24.d0

!      correct
       dd = -4.d0*pi*(alpha*x(ir))**2*ntot(ir) - alpha*w(ir)
       w(ir) = w(ir-1) + (9.d0*dd+19.d0*a-5.d0*b+c)/24.d0

       wp(ir) = -4.d0*pi*(alpha*x(ir))**2*ntot(ir) - alpha*w(ir)

    end do

!
!   matching point
    ir = mp

!      predict
       a  = -4.d0*pi*(alpha*x(ir-1))**2*ntot(ir-1) - alpha*w(ir-1)
       b  = -4.d0*pi*(alpha*x(ir-2))**2*ntot(ir-2) - alpha*w(ir-2)
       c  = -4.d0*pi*(alpha*x(ir-3))**2*ntot(ir-3) - alpha*w(ir-3)
       dd = -4.d0*pi*(alpha*x(ir-4))**2*ntot(ir-4) - alpha*w(ir-4)
       wm = w(ir-1) + (55.d0*a-59.d0*b+37.d0*c-9.d0*dd)/24.d0


!      correct
       dd = -4.d0*pi*(alpha*x(ir))**2*ntot(ir) - alpha*wm
       wm = w(ir-1) + (9.d0*dd+19.d0*a-5.d0*b+c)/24.d0


!   inwards shoot
    do ir = nx-5,mp,-1


!      predict
       a  = -4.d0*pi*(alpha*x(ir+1))**2*ntot(ir+1) - alpha*w(ir+1)
       b  = -4.d0*pi*(alpha*x(ir+2))**2*ntot(ir+2) - alpha*w(ir+2)
       c  = -4.d0*pi*(alpha*x(ir+3))**2*ntot(ir+3) - alpha*w(ir+3)
       dd = -4.d0*pi*(alpha*x(ir+4))**2*ntot(ir+4) - alpha*w(ir+4)
       w(ir) = w(ir+1) - (55.d0*a-59.d0*b+37.d0*c-9.d0*dd)/24.d0

!      correct
       dd = -4.d0*pi*(alpha*x(ir))**2*ntot(ir) - alpha*w(ir)
       w(ir) = w(ir+1) - (9.d0*dd+19.d0*a-5.d0*b+c)/24.d0


       wp(ir) = -4.d0*pi*(alpha*x(ir))**2*ntot(ir) - alpha*w(ir)


   end do


!  homogenous solutions
      A = -(wm - w(mp))*dexp(alpha*real(mp,dp))

      do ir = nx,mp,-1
        w(ir)  =  w(ir)  - A*dexp(-alpha*real(ir,dp))
        wp(ir) =  wp(ir) + alpha*A*dexp(-alpha*real(ir,dp))
      end do



! make the Hartree potential
      do ir = Nx-5,1,-1
        phi(ir) = phi(ir+1) - (9.d0*w(ir)+19.d0*w(ir+1)-5.d0*w(ir+2)+w(ir+3))/24.d0
      end do

!shift to match coulomb at large distances 
      A = phi(Nx)-real(Ne,dp)/x(Nx)
      do ir = 1,Nx
        phi(ir) =  phi(ir) - A
        Vh(ir) = phi(ir)
      end do

!      write(*,*) iteration
      write(char_iter,"(i3)") iteration
      do ir = 1,len(char_iter)
        if ( char_iter(ir:ir) == " " ) char_iter(ir:ir) = "0"
      end do
!      open(353,file="data/hartree"// char_iter //".dat",status="replace")
!      open(354,file="data/analyticH.dat",status="replace")
!      do ir = 2,Nx-1
         xx = x(ir)
         !write(350,*) xx,-(wp(ir)+alpha*w(ir))/(alpha)**2/4.0_dp/pi
         !write(351,*) xx,ntot(ir)*xx*xx
         !write(352,*) xx,wp(ir),w(ir)
!         write(353,*) xx,phi(ir),real(Ne,dp)/x(ir)
!         write(354,*) xx,ana(xx)
         !write(355,*) xx,( ana(x(ir+1)) - ana(x(ir-1)) )/2.0_dp,dana(xx)
!      end do
!      close(353)
!      close(354)
!      stop

  end subroutine hartree


  function ana(r)
    implicit none
    real(dp) :: ana
    real(dp), intent(in) :: r

    ana = real(Ne,dp)*erf(sqrt(sigma)*r)/r

  end function ana



  function dana(r)
    implicit none
    real(dp) :: dana
    real(dp), intent(in) :: r

    dana = alpha*real(Ne,dp)*( 2.0_dp*sqrt(sigma/pi)*exp(-sigma*r**2) - erf(sqrt(sigma)*r)/r )

  end function dana



  !--------------------------------------------!
  !   The first method of making the energy    !
  ! Non-interacting Kinetic energy + potential !
  !--------------------------------------------!
  subroutine Energy(Etot,T,Eh,Enuc,Exc)
    implicit none

    real(dp), intent(out) :: Etot,T,Eh,Enuc,Exc
    real(dp), dimension(Nx) :: temp
    integer :: ix,is
    integer :: ll
    real(dp) :: tl,Ti

    Etot = 0.0_dp

    !kinetic term
    T    = 0.0_dp
    do is = 1,Ne

      ll = ks(is)%l
      tl = real(ll*(ll+1),dp)
      temp = 0.0_dp
      do ix = 1,Nx
        temp(ix) = alpha**2*tl*ks(is)%u(ix) + alpha*ks(is)%up(ix) - ks(is)%upp(ix)
      end do

      Ti = 0.0_dp
      do ix = 1,Nx
        Ti = Ti + temp(ix)*ks(is)%u(ix)/x(ix)
      end do
      Ti = Ti*tpi/alpha
      T = T + Ti
write(*,*) ".....",is,Ti
    end do

    !Hartree energy
    Eh   = 0.0_dp
    do ix = 1,Nx
      Eh = Eh + vh(ix)*ntot(ix)*x3(ix)
    end do
    Eh = Eh*tpi*alpha


    !Nuclear energy
    Enuc = 0.0_dp
    do ix = 1,Nx
      Enuc = Enuc - ntot(ix)*x(ix)**2
    end do
    Enuc = Enuc*4.0_dp*pi*alpha*Z


    !Exchange and correlation energy
    Exc  = 0.0_dp
    do ix = 1,Nx
      Exc = Exc + ntot(ix)*eps_xc(ntot(ix),zeta(ix))*x3(ix)
    end do
    Exc = Exc*fpi*alpha

    Etot = T + Eh + Enuc + Exc

  end subroutine Energy




  !--------------------------------------------!
  !  The second method of making the energy    !
  ! Sum over eigenvalues minus double counting !
  !--------------------------------------------!
  subroutine Energy2(Etot,T,Eh,Exc)
    implicit none

    real(dp), intent(out) :: Etot,T,Eh,Exc

    Etot = 0.0_dp
    T    = 0.0_dp
    Eh   = 0.0_dp
    Exc  = 0.0_dp

  end subroutine Energy2



subroutine testhartree
   integer :: i,j
   real(dp) :: xx,temp
   do i = 2,Nx-1
      temp = (Vh(i+1) - 2.0_dp*Vh(i) + Vh(i-1))
      temp = temp + alpha*(Vh(i+1)-Vh(i-1))*0.5_dp
      temp = temp / (alpha *x(i))**2
      temp = -temp/4.0_dp/pi
      write(55,*) x(i),temp,temp + ntot(i)
   end do

end subroutine testhartree





end module dfsubs
