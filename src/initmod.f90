module initmod
use constants
use globals
use output

implicit none

private

public :: initialise

contains


!---------------------------------------------------------------!
!                  initialize the program                       !
!---------------------------------------------------------------!
subroutine initialise

   integer :: is,i
   integer :: no_states
   real(dp) :: setrmax,gridscale
   real(dp) :: dk,delta,zion,sf


   !read all parameters form file
   open(1,file="inp.dat.atom",status="old",action="read")


   !space properties
   read (unit=1,fmt=*) Nx
   read (unit=1,fmt=*) Xmax
   read (unit=1,fmt=*) Xmin
   read (unit=1,fmt=*) mix
   read (unit=1,fmt=*) max_iter
   read (unit=1,fmt=*) tol

   !wire properties
   read (unit=1,fmt=*) Z
   read (unit=1,fmt=*) Ne
   read (unit=1,fmt=*) Nvirt
   read (unit=1,fmt=*) eft       !electronic fermi temperature


   !space for the orbitals
   Nstate = Ne + Nvirt
   allocate(vec(Nx,3,Nstate))
   allocate(ks(Nstate))
   do i = 1,Nstate
     ks(i)%u   => vec(:,1,i)
     ks(i)%up  => vec(:,2,i)
     ks(i)%upp => vec(:,3,i)
   end do


   !run details
   read (unit=1,fmt=*) read_dens_file
   read (unit=1,fmt=*) output_level
   read (unit=1,fmt=*) write_psirays


   !read the desired ground state shell structure
   read (unit=1,fmt=*)
   read (unit=1,fmt=*)
   do is = 1,ne
     write(*,*) is
     read(unit=1,fmt=*) ks(is)%n,ks(is)%l,ks(is)%m,ks(is)%s
   end do

   !read the desired ground state shell structure
   read (unit=1,fmt=*)
   read (unit=1,fmt=*)
   write(*,*) "here"
   do is = ne+1,ne+nvirt
     read(unit=1,fmt=*) ks(is)%n,ks(is)%l,ks(is)%m,ks(is)%s
   end do
   write(*,*) "here2"



   close(1)


!  occupy the filled
  do is = 1,Nstate
    ks(is)%f = 0.0_dp
  end do
  do is = 1,Ne
    ks(is)%f = 1.0_dp
  end do
  do is = Ne+1,Nstate
    ks(is)%f = 0.0_dp
  end do

!  count the up and down states
  Nup = 0
  do is = 1,Nstate
    if (ks(is)%s>0.and.ks(is)%f>0.99_dp) Nup = Nup + 1
  end do
  Ndn = 0
  do is = 1,Nstate
    if (ks(is)%s<0.and.ks(is)%f>0.99_dp) Ndn = Ndn + 1
  end do
  write(*,*) "Up count:",Nup
  write(*,*) "Dn count:",Ndn
  write(*,*)
  write(*,*)


!  guess the eigenvalues
   sf=0.0d0
   do is=1,Nstate
     sf=sf+ks(is)%f
     zion=z+1.0d0-sf
     ks(is)%e=-0.5d0*((zion)/ks(is)%n)**2
!     if(is.gt.1)then
!       if(n(is).eq.n(is-1).and.l(is).eq.l(is-1))eps(is)=eps(is-1)
!     end if
   enddo
   ks(10)%e = 1.0_dp




   call shell_to_screen(1,Nstate)


   !scaled radial grid
   write(*,*) "Making space grid..."
   allocate(x(Nx),x3(Nx))

   alpha = log(Xmax/Xmin) / real(Nx-1,dp)              !initial spacing


   open(20,file="data/space.dat",status="replace")
   do i = 1,Nx
       x(i) = exp(real(i-1,dp)*alpha)*Xmin
       x3(i) = x(i)**3
       write(20,*) i,x(i),x3(i)    !working
   end do
   close(20)


   write (*,*) "Xmin",Xmin
   write (*,*) "Xmax",Xmax
   write (*,*) "alpha",alpha
   write (*,*) "                   ...done!"
   write (*,*)





   !return input parameters to user for clarity
   call write_run_data


   !make space for arrays
   allocate(Vh(nx))
   allocate(VKS(Nx,2))        !kohn-sham without coulomb



   !energetic ordering of the states with this table
   allocate(state_index(Nstate))


   !component densities
   allocate(ntot(nx))
   allocate(n(nx,2))
   allocate(nold(nx,2))


   !total densities
   allocate(zeta(Nx))


   !open output files
   open(200,file="energy/E.dat"  ,status="replace")
   open(201,file="energy/Eh.dat" ,status="replace")
   open(202,file="energy/Exc.dat",status="replace")
   open(203,file="energy/T.dat"  ,status="replace")


   
   call Angular_mesh


end subroutine initialise





   !------------------------------------!
   !   Make the Sobolov Mesh for the    !
   !        angular integrals           !
   !------------------------------------!
   subroutine Angular_Mesh
      implicit none

      integer :: i,j,k
      integer :: ilm,jlm
      integer :: ip,ipp,ix,il,is

      real(dp) :: r2,rr,s2,ss,u2,uu,v2,vv,t2,tt
      real(dp) :: rmin,the,ph,sum
      real(dp) :: al,alpha,beta
      real(dp) :: sc,r2d,xx,y,zz

!      real(dp), dimension(3,npx) :: pp
      real(dp), dimension(3,3) :: rot
      real(dp), dimension(100) :: P
      real(dp), dimension(3) :: aux
      !dimension pp(3,npx),rot(3,3),aux(3),P(100)

!     
! Mesh in theta and phi
!
       r2=(5.d0+dsqrt(5.d0))/10.d0
       rr=dsqrt(r2)
       s2=(5.d0-dsqrt(5.d0))/10.d0
       ss=dsqrt(s2)
       u2=(3.d0-dsqrt(5.d0))/6.d0
       uu=dsqrt(u2)
       v2=(3.d0+dsqrt(5.d0))/6.d0
       vv=dsqrt(v2)
       t2=1.d0/3.d0
       tt=dsqrt(t2)
!
       pp(1,1)=rr
       pp(2,1)=ss
       pp(3,1)=0.d0
!
       pp(1,2)=-rr
       pp(2,2)=ss
       pp(3,2)=0.d0
!
       pp(1,3)=rr
       pp(2,3)=-ss
       pp(3,3)=0.d0
!
       pp(1,4)=-rr
       pp(2,4)=-ss
       pp(3,4)=0.d0
!
       pp(1,5)=0.d0
       pp(2,5)=rr
       pp(3,5)=ss
!
       pp(1,6)=0.d0
       pp(2,6)=-rr
       pp(3,6)=ss
!
       pp(1,7)=0.d0
       pp(2,7)=rr
       pp(3,7)=-ss
!
       pp(1,8)=0.d0
       pp(2,8)=-rr
       pp(3,8)=-ss
!
       pp(1,9)=ss
       pp(2,9)=0.d0
       pp(3,9)=rr
!
       pp(1,10)=-ss
       pp(2,10)=0.d0
       pp(3,10)=rr
!
       pp(1,11)=ss
       pp(2,11)=0.d0
       pp(3,11)=-rr
!
       pp(1,12)=-ss
       pp(2,12)=0.d0
       pp(3,12)=-rr
!
       pp(1,13)=uu
       pp(2,13)=vv
       pp(3,13)=0.d0
!
       pp(1,14)=-uu
       pp(2,14)=vv
       pp(3,14)=0.d0
!
       pp(1,15)=uu
       pp(2,15)=-vv
       pp(3,15)=0.d0
!
       pp(1,16)=-uu
       pp(2,16)=-vv
       pp(3,16)=0.d0
!
       pp(1,17)=0.d0
       pp(2,17)=uu
       pp(3,17)=vv
!
       pp(1,18)=0.d0
       pp(2,18)=-uu
       pp(3,18)=vv
!
       pp(1,19)=0.d0
       pp(2,19)=uu
       pp(3,19)=-vv
!
       pp(1,20)=0.d0
       pp(2,20)=-uu
       pp(3,20)=-vv
!
       pp(1,21)=vv
       pp(2,21)=0.d0
       pp(3,21)=uu
!
       pp(1,22)=-vv
       pp(2,22)=0.d0
       pp(3,22)=uu
!
       pp(1,23)=vv
       pp(2,23)=0.d0
       pp(3,23)=-uu
!
       pp(1,24)=-vv
       pp(2,24)=0.d0
       pp(3,24)=-uu
!
       pp(1,25)=tt
       pp(2,25)=tt
       pp(3,25)=tt
!
       pp(1,26)=-tt
       pp(2,26)=tt
       pp(3,26)=tt
!
       pp(1,27)=tt
       pp(2,27)=-tt
       pp(3,27)=tt
!
       pp(1,28)=tt
       pp(2,28)=tt
       pp(3,28)=-tt
!
       pp(1,29)=tt
       pp(2,29)=-tt
       pp(3,29)=-tt
!
       pp(1,30)=-tt
       pp(2,30)=tt
       pp(3,30)=-tt
!
       pp(1,31)=-tt
       pp(2,31)=-tt
       pp(3,31)=tt
!
       pp(1,32)=-tt
       pp(2,32)=-tt
       pp(3,32)=-tt
!
! Here apply a small random rotation
!
!        alpha=0.181560d0
!        beta=0.063491d0
        alpha=2.181560d0
        beta=1.063491d0
        rot(1,1)=dcos(alpha)*dcos(beta)
        rot(2,1)=-dsin(alpha)*dcos(beta)
        rot(3,1)=-dsin(beta)
        rot(1,2)=dsin(alpha)
        rot(2,2)=dcos(alpha)
        rot(3,2)=0.d0
        rot(1,3)=dsin(beta)*dcos(alpha)
        rot(2,3)=-dsin(alpha)*dsin(beta)
        rot(3,3)=dcos(beta)
        do  ip=1,npx
          aux(1)=0.d0
          aux(2)=0.d0
          aux(3)=0.d0
          do ix=1,3
            aux(1)=aux(1)+rot(1,ix)*pp(ix,ip)
            aux(2)=aux(2)+rot(2,ix)*pp(ix,ip)
            aux(3)=aux(3)+rot(3,ix)*pp(ix,ip)
          end do
          pp(1,ip)=aux(1)
          pp(2,ip)=aux(2)
          pp(3,ip)=aux(3)
!          all normalised
!          write(*,*) pp(1:3,ip)
!          write(*,*) sqrt(pp(1,ip)**2+pp(2,ip)**2+pp(3,ip)**2)
       end do
!       stop
!
! Set up the L^2 matrix
!
       do ip=1,npx
         do ipp=1,npx
           sc=pp(1,ip)*pp(1,ipp)+pp(2,ip)*pp(2,ipp)+pp(3,ip)*pp(3,ipp)
           P(1)=sc
           P(2)=0.5d0*(3.d0*sc*sc-1.d0)
           P(3)=0.5d0*sc*(5.d0*sc*sc-3.d0)
           P(4)=0.125d0*(35.d0*sc**4-30.d0*sc**2+3.d0)
           P(5)=(63.d0*sc**5-70.d0*sc**3+15.d0*sc)/8.d0
           P(6)=(231.d0*sc**6-315.d0*sc**4+105.d0*sc*sc-5.d0)/16.d0
           do il=6,19
             P(il+1)=(il+il+1.d0)*sc*P(il)-il*P(il-1)
             P(il+1)=P(il+1)/(il+1.d0)
           end do
           el2(ip,ipp)=0.d0
           do il=1,6
             el2(ip,ipp)=el2(ip,ipp)+dfloat((il+il+1)*il*(il+1))*P(il)
           end do
           el2(ip,ipp)=el2(ip,ipp)/fpi
         end do
       end do
!        do 500 ip=1,npx
!        do 501 ipp=1,npx
!        write(76,*)el2(ip,ipp),el2(ipp,ip)
! 501    continue
! 500    continue
!        if(npx.ge.0)stop
!
! Here define the 2D coordinates and the weight
!
     do ip=1,npx
       r2d=dsqrt(pp(1,ip)**2+pp(2,ip)**2)
       the=0.d0
       if (r2d .gt. 1.d-12) the=dacos(pp(3,ip))
       ph=0.d0
       if (dabs(pp(1,ip)) .gt. 1.d-12) ph=datan(pp(2,ip)/pp(1,ip))
       if (dsin(the)*dcos(ph)*pp(1,ip) .le. 0.d0) ph = ph + pi
       thephi(1,ip)=the
       thephi(2,ip)=ph
       if(ip.le.12)then
         thephi(3,ip)=25.d0/840.d0*fpi
       else
         thephi(3,ip)=27.d0/840.d0*fpi
       end if
     end do
!
! Define the REAL spherical harmonics
! In this version Ylm up to l=4 are included.
       do ip=1,npx
       xx=pp(1,ip)
       y=pp(2,ip)
       zz=pp(3,ip)
       ylm(ip,1)=1.d0/sqfpi
       ylm(ip,2)=-dsqrt(3.d0)/sqfpi*pp(1,ip)
       ylm(ip,3)=-dsqrt(3.d0)/sqfpi*pp(2,ip)
       ylm(ip,4)=dsqrt(3.d0)/sqfpi*pp(3,ip)
       ylm(ip,5)=dsqrt(15.d0)/sqfpi/2.d0*(pp(1,ip)**2-pp(2,ip)**2)
       ylm(ip,6)=dsqrt(15.d0)/sqfpi*(pp(1,ip)*pp(2,ip))
       ylm(ip,7)=-dsqrt(15.d0)/sqfpi*(pp(1,ip)*pp(3,ip))
       ylm(ip,8)=-dsqrt(15.d0)/sqfpi*(pp(2,ip)*pp(3,ip))
       ylm(ip,9)=dsqrt(5.d0)/sqfpi/2.d0*(3.d0*pp(3,ip)**2-1.d0)
       ylm(ip,10)=-dsqrt(35.d0/(32.d0*pi))*(xx**3-3.d0*xx*y*y)
       ylm(ip,11)=-dsqrt(35.d0/(32.d0*pi))*(y**3-3.d0*xx*xx*y)
       ylm(ip,12)=dsqrt(105.d0/(16.d0*pi))*zz*(xx**2-y**2)
       ylm(ip,13)=dsqrt(105.d0/(16.d0*pi))*2.d0*xx*y*zz
       ylm(ip,14)=-dsqrt(21.d0/(32.d0*pi))*(5.d0*xx*zz*zz-xx)
       ylm(ip,15)=-dsqrt(21.d0/(32.d0*pi))*(5.d0*y*zz*zz-y)
       ylm(ip,16)=dsqrt(7.d0/(16.d0*pi))*(5.d0*zz**3-3.d0*zz)
       ylm(ip,17)=dsqrt(315.0d0/(256.d0*pi))*(xx**4-6.0d0*xx**2*y**2+y**4)
       ylm(ip,18)=dsqrt(315.0d0/(16.0d0*pi))*(y*xx**3-xx*y**3)
       ylm(ip,19)=-dsqrt(315.0d0/(32.0d0*pi))*zz*(xx**3-3.0d0*xx*y**2)
       ylm(ip,20)=-dsqrt(315.0d0/(32.0d0*pi))*zz*(3.0d0*xx**2*y-y**3)
       ylm(ip,21)=dsqrt(45.0d0/(64.0d0*pi))*(7.0d0*zz**2 - 1.0d0)&
     &      *(xx**2 - y**2)
       ylm(ip,22)=dsqrt(45.0d0/(16.0d0*pi))*xx*y*(7.0d0*zz**2-1.0d0)
       ylm(ip,23)=-dsqrt(45.0d0/(32.0d0*pi))*xx*(7.0d0*zz**3-3.0d0*zz)
       ylm(ip,24)=-dsqrt(45.0d0/(32.0d0*pi))*y*(7.0d0*zz**3-3.0d0*zz)
       ylm(ip,25)=dsqrt(9.0d0/(256.0d0*pi))*(35.0d0*zz**4-30.0d0*&
     &      zz**2+3.0d0)
     end do
!
! Test
!     
!       do 333 ilm=1,25
!          do 334 jlm=1,25
!             sum=0.d0
!             do 335 ip=1,npx
!                sum=sum+ylm(ip,ilm)*ylm(ip,jlm)*thephi(3,ip)
! 335         continue
!             write(6,*)ilm,jlm,sum
! 334      continue
! 333   continue
!       stop
!
! List the angular momentum of each Ylm
!
       llm(1)=0
       llm(2)=1
       llm(3)=1
       llm(4)=1
       llm(5)=2
       llm(6)=2
       llm(7)=2
       llm(8)=2
       llm(9)=2
       llm(10)=3
       llm(11)=3
       llm(12)=3
       llm(13)=3
       llm(14)=3
       llm(15)=3
       llm(16)=3
       llm(17)=4
       llm(18)=4
       llm(19)=4
       llm(20)=4
       llm(21)=4
       llm(22)=4
       llm(23)=4
       llm(24)=4
       llm(25)=4
!
! Identify the lm of each state close to the origin
!
     do is=1,nstate

       ks(is)%lmstate=0
       !l=0
       if(ks(is)%l.eq.0.and.ks(is)%m.eq. 0)ks(is)%lmstate=1
       !l=1
       if(ks(is)%l.eq.1.and.ks(is)%m.eq. 1)ks(is)%lmstate=2
       if(ks(is)%l.eq.1.and.ks(is)%m.eq.-1)ks(is)%lmstate=3
       if(ks(is)%l.eq.1.and.ks(is)%m.eq. 0)ks(is)%lmstate=4
       !l=2
       if(ks(is)%l.eq.2.and.ks(is)%m.eq. 2)ks(is)%lmstate=5
       if(ks(is)%l.eq.2.and.ks(is)%m.eq.-2)ks(is)%lmstate=6
       if(ks(is)%l.eq.2.and.ks(is)%m.eq. 1)ks(is)%lmstate=7
       if(ks(is)%l.eq.2.and.ks(is)%m.eq.-1)ks(is)%lmstate=8
       if(ks(is)%l.eq.2.and.ks(is)%m.eq. 0)ks(is)%lmstate=9
       !l=3
       if(ks(is)%l.eq.3.and.ks(is)%m.eq. 3)ks(is)%lmstate=10
       if(ks(is)%l.eq.3.and.ks(is)%m.eq.-3)ks(is)%lmstate=11
       if(ks(is)%l.eq.3.and.ks(is)%m.eq. 2)ks(is)%lmstate=12
       if(ks(is)%l.eq.3.and.ks(is)%m.eq.-2)ks(is)%lmstate=13
       if(ks(is)%l.eq.3.and.ks(is)%m.eq. 1)ks(is)%lmstate=14
       if(ks(is)%l.eq.3.and.ks(is)%m.eq.-1)ks(is)%lmstate=15
       if(ks(is)%l.eq.3.and.ks(is)%m.eq. 0)ks(is)%lmstate=16
       !l=4
       if(ks(is)%l.eq.4.and.ks(is)%m.eq. 4)ks(is)%lmstate=17
       if(ks(is)%l.eq.4.and.ks(is)%m.eq.-4)ks(is)%lmstate=18
       if(ks(is)%l.eq.4.and.ks(is)%m.eq. 3)ks(is)%lmstate=19
       if(ks(is)%l.eq.4.and.ks(is)%m.eq.-3)ks(is)%lmstate=20
       if(ks(is)%l.eq.4.and.ks(is)%m.eq. 2)ks(is)%lmstate=21
       if(ks(is)%l.eq.4.and.ks(is)%m.eq.-2)ks(is)%lmstate=22
       if(ks(is)%l.eq.4.and.ks(is)%m.eq. 1)ks(is)%lmstate=23
       if(ks(is)%l.eq.4.and.ks(is)%m.eq.-1)ks(is)%lmstate=24
       if(ks(is)%l.eq.4.and.ks(is)%m.eq. 0)ks(is)%lmstate=25
       !otherwise
       if(ks(is)%lmstate.eq.0)then
         write(6,*)' error in input: I cannot identify '
         write(6,*)' the angular character of state ',is
         stop
       end if

     end do 
!
! Set up the small basis for the fit of the charge density
!
!       despo=2.d0*z/nespox
!       do 24 ie=1,nespox
!       espo(ie)=2.d0*z-(ie-0.5d0)*despo
!       write(6,*)ie,espo(ie)
! 24    continue
!
! Set the ``external potential''


!
! Set up the small basis for the fit of the charge density
!
!       despo=2.d0*z/nespox
!       do 24 ie=1,nespox
!       espo(ie)=2.d0*z-(ie-0.5d0)*despo
!       write(6,*)ie,espo(ie)
! 24    continue
!
! Set the ``external potential''
!       call perturbation        !don't worry about the external field for now
!
!      if(nbeg.lt.0)call startrho
!
!       do 555 is=1,nx
!       do 556 ip=1,npx
!       do 557 ir=1,nrx
!          psi(ir,ip,is)=0.d0
! 557   continue
! 556   continue
! 555   continue

   end subroutine Angular_Mesh




end module initmod
