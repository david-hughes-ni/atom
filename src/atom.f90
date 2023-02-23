program atom


use constants
use globals
use dfsubs
use solvers
use initmod
use output



implicit none


  logical :: equal = .false.
  integer :: i,j
  real(dp) :: Etot,T,Eh,Enuc,Exc


  !starting off
  call initialise
  call start_density
  call build_potential



  !ground state self-consistency
  do iteration = 1,max_iter

     call get_states(1,Ne)
     call update_density( equal )
     if ( equal ) exit
     call build_potential

  end do


  !data for the user
  write(*,*) "Ground State in",iteration,"iterations."
  write(*,*)
  write(*,*)
  call Shell_to_screen(1,Ne)
  write(*,*)
  write(*,*)
  call Energy(Etot,T,Eh,Enuc,Exc)
  write(*,*) "---------------ENERGY------------"
  write(*,*) "kinetic:",T
  write(*,*) "Hartree:",Eh
  write(*,*) "Nuclear:",Enuc
  write(*,*) "XC:",Exc
  write(*,*) "Total:",Etot
  write(*,*)
  write(*,*)



  !Empty state self-consistency
  call build_potential
  call get_states(Ne+1,Nstate)
     
  call Energy(Etot,T,Eh,Enuc,Exc)
  write(*,*) "---------------ENERGY------------"
  write(*,*) "kinetic:",T
  write(*,*) "Hartree:",Eh
  write(*,*) "Nuclear:",Enuc
  write(*,*) "XC:",Exc
  write(*,*) "Total:",Etot
  write(*,*)
  write(*,*)


  !more data
  call Shell_to_screen(1,Nstate)
  call Shell_to_screen_latex(1,Nstate)

  call Finalise

!  call Energy2(Etot,T,Eh,Exc)




  !energy component files
  close(200)
  close(201)
  close(202)
  close(203)


end program atom
