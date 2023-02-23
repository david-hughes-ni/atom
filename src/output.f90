!----------------------------------------------------!
!   All output from the code.  Both file and screen  !
!----------------------------------------------------!
module output

use constants
use globals
use pxcfuncs


implicit none

private


public :: write_run_data
public :: final_data_to_screen
public :: Finalise 
public :: Orbitals_to_file
public :: shell_to_screen
public :: shell_to_screen_latex



contains







!--------------------------------------------!
! write states to file with systematic names !
!  also does an unformatted dump to psi.dat  !
!--------------------------------------------!
subroutine Orbitals_to_file
  implicit none

  character(14) :: file_name
  character(3) :: n_name,l_name,s_name
  integer :: is,j,ix,ilm,ip,ispin
  real(dp) :: orb(nx,npx)
  real(dp) :: dorb(nx,npx)
  real(dp) :: d2orb(nx,npx)

!   if (run%write_orbs) then


   do is = 1,Nstate


      if ( ks(is)%m == 0 ) then

         !wipe the old names
         n_name = "   "
         l_name = "   "
         s_name = "   "
         file_name   = "   "

         !make the name
         write(n_name,"(i3)") ks(is)%n
         write(l_name,"(i3)") ks(is)%l
         do j = 1,3
            if (n_name(j:j).eq." ") n_name(j:j)="_"
            if (l_name(j:j).eq." ") l_name(j:j)="_"
         end do
         if (ks(is)%s == 1) then
            write(s_name,"(a3)") "__u"
         else
            write(s_name,"(a3)") "__d"
         end if
         file_name = n_name// l_name // s_name // ".dat"

         open(unit=(3000),file="states/"//file_name,status="replace")
         do j = 1,Nx
            write(3000,30) x(j),ks(is)%u(j)*sqfpi!/x(j),ks(is)%up(j)/x(j),ks(is)%upp(j)/x(j)
            30 format(4(f15.10,4x))
         end do
         close(3000)


      end if

   end do


 !  end if

   if ( write_psirays ) then


   !make them orthogonal if they are not already
   call gram_radial_vectors


   !read the basis states along with the shell structure
   open(2,file="psi.dat",status="replace",action="write",form="unformatted")
   do is = 1,Nstate

     if (ks(is)%s>0) then
     write(2) ks(is)%n
     write(2) ks(is)%l
     write(2) ks(is)%m
     write(2) ks(is)%s
     write(2) ks(is)%e
     write(2) ks(is)%f
     write(2) ks(is)%u(:)
     write(2) ks(is)%up(:)
     write(2) ks(is)%upp(:)
     end if

   end do
   do is = 1,Nstate

     if (ks(is)%s<0) then
     write(2) ks(is)%n
     write(2) ks(is)%l
     write(2) ks(is)%m
     write(2) ks(is)%s
     write(2) ks(is)%e
     write(2) ks(is)%f
     write(2) ks(is)%u(:)
     write(2) ks(is)%up(:)
     write(2) ks(is)%upp(:)
     end if

   end do
   close(2)




   !wave function along the rays
   open(2,file="psi_rays.dat",status="replace",action="write",form="unformatted")
   do ispin = 1,-1,-2


     do is = 1,Nstate


     if ( ks(is)%s == ispin ) then


     !build  the orbital
     ilm = ks(is)%lmstate
!     write(*,*) is,ilm
       do ix = 1,nx
         do ip = 1,npx
           orb(ix,ip)   = ks(is)%u(ix)*ylm(ip,ilm)*sqfpi
           dorb(ix,ip)  = ks(is)%up(ix)*ylm(ip,ilm)*sqfpi
           d2orb(ix,ip) = ks(is)%upp(ix)*ylm(ip,ilm)*sqfpi
         end do
       end do



         write(2) orb
         write(2) dorb
         write(2) d2orb


     end if


     end do


   end do

   end if


end subroutine Orbitals_to_file



   !------------------------------------!
   !  Orthogonalise the radial vectors  !
   !         from the atom code         !
   !   Do this after the angular mesh   !
   !          has been made!!!          !
   !------------------------------------!
   subroutine gram_radial_vectors
     implicit none

      integer :: ix,is,js
      integer :: i,j,ilm,jlm
      real(dp) :: sum

!
! Spin Up
!
      do is=1,Nstate

         if ( ks(is)%s > 0 ) then
         ilm = ks(is)%lmstate


         do js=1,is-1

           if ( ks(js)%s > 0 ) then

             jlm = ks(js)%lmstate
             if ( ilm == jlm ) then
 
              sum = zero
              do ix = 1,Nx
                sum = sum + ks(is)%u(ix)*ks(js)%u(ix)*x(ix)
              end do
              sum = sum*alpha*fpi

              do ix = 1,Nx
                ks(is)%u(ix)   = ks(is)%u(ix)   - sum*ks(js)%u(ix)
                ks(is)%up(ix)  = ks(is)%up(ix)  - sum*ks(js)%up(ix)
                ks(is)%upp(ix) = ks(is)%upp(ix) - sum*ks(js)%upp(ix)
              end do

             end if

           end if


         end do


         sum = zero
         do ix = 1,Nx
            sum = sum + ( ks(is)%u(ix) )**2*x(ix)
         end do
         sum = sum*alpha*fpi
         sum = one/sqrt(sum)


         do ix = 1,Nx
           ks(is)%u(ix)   = sum*ks(is)%u(ix)
           ks(is)%up(ix)  = sum*ks(is)%up(ix)
           ks(is)%upp(ix) = sum*ks(is)%upp(ix)
         end do

        end if

      end do


!
! Spin Down
!
      do is=1,Nstate
         if ( ks(is)%s < 0 ) then
         ilm = ks(is)%lmstate


         do js=1,is-1

           if ( ks(js)%s < 0 ) then

             jlm = ks(js)%lmstate
             if ( ilm == jlm) then

              sum = zero
              do ix = 1,Nx
                sum = sum + ks(is)%u(ix)*ks(js)%u(ix)*x(ix)
              end do
              sum = sum*alpha*fpi

              do ix = 1,Nx
                ks(is)%u(ix)   = ks(is)%u(ix)   - sum*ks(js)%u(ix)
                ks(is)%up(ix)  = ks(is)%up(ix)  - sum*ks(js)%up(ix)
                ks(is)%upp(ix) = ks(is)%upp(ix) - sum*ks(js)%upp(ix)
              end do

             end if

           end if


         end do


         sum = zero
         do ix = 1,Nx
            sum = sum + ( ks(is)%u(ix) )**2*x(ix)
         end do
         sum = sum*alpha*fpi
         sum = one/sqrt(sum)


         do ix = 1,Nx
           ks(is)%u(ix)   = sum*ks(is)%u(ix)
           ks(is)%up(ix)  = sum*ks(is)%up(ix)
           ks(is)%upp(ix) = sum*ks(is)%upp(ix)
         end do

        end if

      end do



    !  open(100,file="ocoefup.dat",status="replace")
    !  open(101,file="ocoefdn.dat",status="replace")
    !  do i = 1,nup_tot
    !    do j = 1,nup_tot
    !      write(100,*) i,j,(cmat(i,j,1))
    !      write(101,*) i,j,(cmat(i,j,2))
    !    end do
    !    write(100,*)
    !    write(101,*)
    !  end do
    !  close(100)
    !  close(101)



   end subroutine gram_radial_vectors 


   !-------------------------------------------!
   !     write shell structure to screen       !
   !             K = 0 points only             !
   !-------------------------------------------!
   subroutine shell_to_screen(ia,ib)
     implicit none

     integer, intent(in) :: ia,ib
     integer :: i,j,si,ii

      write(*,*) "    ----------------------------------------------------"
      write(*,*) "    |       Atomic Shell Structure                     |"
      write(*,*) "    ----------------------------------------------------"
      write(*,*) "    | n  |  l  |  m  |  s  |    f    |        e        |"
      write(*,*) "    ----------------------------------------------------"

      ii = 1
      do i = ia,ib
         if ( ks(i)%s>0 ) then
            !shell structure at gamma point only
            write(*,'(i5,a,1x,i3,2x,i3,2x,i3,2x,i3,2x,f10.6,2x,f20.12)')&
                  & ii,"|",ks(i)%n,ks(i)%l&
                  &,ks(i)%m,ks(i)%s,ks(i)%f,ks(i)%e
            ii = ii + 1
         end if
      end do
      do i = ia,ib
         if ( ks(i)%s<0 ) then
            !shell structure at gamma point only
            write(*,'(i5,a,1x,i3,2x,i3,2x,i3,2x,i3,2x,f10.6,2x,f20.12)')&
                  & ii,"|",ks(i)%n,ks(i)%l&
                  &,ks(i)%m,ks(i)%s,ks(i)%f,ks(i)%e
            ii = ii + 1
         end if
      end do
      write(*,*) "   -----------------------------------------------------"
      write(*,*)
      write(*,*)


   end subroutine shell_to_screen




   subroutine shell_to_screen_latex(ia,ib)
     implicit none

     integer, intent(in) :: ia,ib
     integer :: i,j,si,ii

      write(*,*)
      write(*,*)


      ii = 1
      do i = ia,ib
         if ( ks(i)%s>0 ) then
            !shell structure at gamma point only
            if (ks(i)%l==0)then
              write(*,'(i5,a,f12.6,a)')&
                  & ks(i)%n, "s$^\uparrow$ & ",ks(i)%e," \\"
            else if(ks(i)%l==1)then
              write(*,'(i5,a,f12.6,a)')&
                  & ks(i)%n, "p$^\uparrow$ & ",ks(i)%e," \\"
            else if(ks(i)%l==2)then
              write(*,'(i5,a,f12.6,a)')&
                  & ks(i)%n, "d$^\uparrow$ & ",ks(i)%e," \\"
            end if
            ii = ii + 1
         end if
      end do
      do i = ia,ib
         if ( ks(i)%s<0 ) then
            !shell structure at gamma point only
            if (ks(i)%l==0)then
              write(*,'(i5,a,f12.6,a)')&
                  & ks(i)%n, "s$^\downarrow$ & ",ks(i)%e," \\"
            else if(ks(i)%l==1)then
              write(*,'(i5,a,f12.6,a)')&
                  & ks(i)%n, "p$^\downarrow$ & ",ks(i)%e," \\"
            else if(ks(i)%l==2)then
              write(*,'(i5,a,f12.6,a)')&
                  & ks(i)%n, "d$^\downarrow$ & ",ks(i)%e," \\"
            end if
            ii = ii + 1
         end if
      end do


   end subroutine shell_to_screen_latex


   !------------------------------------------!
   !     final properties to screen           !
   !      -integrated spin                    !
   !      -integrated spin per electron
   !      -total energy                       !
   !      -energy per electron
   !------------------------------------------!
   subroutine final_data_to_screen
      integer :: i
      real(dp) :: temp

      !integrated spin polarisation
!      temp = 0.5_dp*( x(1)*abs(nup(1,1)-ndn(1,1)) + x(Nx)*abs(nup(1,Nx)-ndn(1,Nx)) )
!      do i = 2,Nx-1
!         temp = temp + x(i)*abs(nup(1,i)-ndn(1,i))
!      end do
!      temp = temp*dx*2.0_dp*pi*Length*Rb*rb
      write(*,'(a)') "---------------------------------------------------"
      write(*,'(a)') "|----- ----- |\    | ----- ----- |    | |---- |-\"
      write(*,'(a)') "|        |   | \   |   |   |     |    | |     |  \"
      write(*,'(a)') "|-----   |   |  \  |   |   ----- |----| |---- |   | "
      write(*,'(a)') "|        |   |   \ |   |       | |    | |     |  / "
      write(*,'(a)') "|      ----- |    \| ----- ----- |    | |---- |-/"
!      write(*,'(a)') "-------OUTPUT PROPERTIES------------"
!      write(*,*) "Fermi energy",Ef
!      write(*,'(a,5x,e18.12)') "Total spin polarisation:",temp
!      write(*,*) "Total s.p. per electron:",temp/real(Ne,dp)
!      write(*,'(a,11x,f10.6)') "Fermi Level:",ef
      write(*,*)
      write(*,*)
      write(*,*)

   end subroutine final_data_to_screen




  !------------------------------------------!
  !         all essential data to file       !
  !           (file unit =  300+i)           !
  !------------------------------------------!
  subroutine Finalise
    implicit none

    integer :: i,j,si
    real(dp) :: rhoback

!   rhoback = 3.0_dp/4.0_dp/pi*(1.0_dp/rs)**3


    !all densities
    open(302,file="data/dens.dat",status="replace")
    do j = 1, nx
      write(302,*) x(j),n(j,1:2)
    end do
    close(302)

    !potentials
    open(302,file="data/hartree.dat",status="replace")
    do j = 1, nx
      write(302,*) x(j),vh(j)
    end do
    close(302)

!    call gram_radial_vectors
    call Orbitals_to_file


   end subroutine Finalise




   !-------------------------------------!
   !   all details of the run to screen  !
   !-------------------------------------!
   subroutine write_run_data

!      write (*,'(a)')  "-----------USER INPUT-----------"
!      write (*,'(a,25x,i5)')     "Ne:",Ne
!      write(*,*) "Full K point calc.",run%full_K
!      write(*,*) "number of K points",Nk
!      write (*,'(a,25x,f10.4)')  "Rb:",Rb
!      write (*,'(a,26x,f10.4)')  "L:",Length
!      write (*,'(a,25x,f10.4)')  "rs:",rs
!      write (*,'(a,23x,f16.10)') "Xmax:",Xmax
!      write (*,'(a,23x,f16.10)') "Rmax:",Xmax*Rb
!      write (*,'(a,18x,i5)')     "Grid points:",Nx
!      write (*,'(a,24x,a10)')    "Mixing:",run%mixing
!      write (*,'(a,8x,f10.4)')   "DFT mix parameter:",run%alpha
!      write (*,*) "Gamma Point only", run%gamma_shell_to_screen
!      write (*,'(a,5x,f10.4)')   "Electronic fermi temp.:", eft
!      write(*,*)
!      write(*,*)
!      write(*,*)

   end subroutine write_run_data






end module output
