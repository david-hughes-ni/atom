!----------------------------------------------!
!  Routine to get the bound and empty states   !
!----------------------------------------------!
module solvers

   use constants
   use globals
   use dfsubs
   use output

   implicit none



   private
   public :: get_states

   contains




subroutine get_states(ia,ib)
  implicit none

  integer, intent(in) :: ia,ib
  integer :: n,l,m,s
  integer :: is,ix
  real(dp) :: eps
  real(dp), dimension(Nx) :: U,Up,Upp


  do is = ia,ib

    n = ks(is)%n
    l = ks(is)%l
    m = ks(is)%m
    s = ks(is)%s
    eps = ks(is)%e

    call lscheq(n,l,m,s,eps,U,Up,Upp)

    ks(is)%e   = eps
    ks(is)%u(:)   = u(:)
!    do ix = 1,Nx
!      write(100+is,*) x(ix),ks(is)%u(ix),u(ix)
!    end do
    ks(is)%up(:)  = up(:)
    ks(is)%upp(:) = upp(:)

  end do
!  stop


end subroutine get_states






!-----------------------------------------------!
!   make file names to store the eigen value    !
!   sweep information and the output from the   !
!            bisection algorithm                !
!          (this is fantastic!!!!)              !
!-----------------------------------------------!
subroutine open_esweep_file( spin, m )
   implicit none

   integer, intent(in) :: spin, m
   character(len=3) :: iter,m_name,s_name
   character(14) :: sweep_name
   character(20) :: bisect_name
   integer :: j

   iter   = "   "
   m_name = "   "
   s_name = "   "

   write(iter,"(i3)") iteration
   write(m_name,"(i3)") m

   !replace blanks
   do j = 1,3
      if (iter(j:j).eq." ") iter(j:j)="_"
      if (m_name(j:j).eq." ") m_name(j:j)="_"
   end do
   if (spin == 1) then
      write(s_name,"(a3)") "__u"
   else
      write(s_name,"(a3)") "__d"
   end if
   sweep_name = iter// m_name // s_name // ".dat"
   bisect_name = "bisect" // iter// m_name // s_name // ".dat"

   open(unit=1234,file="solvers_debug/"//sweep_name,status="replace")
   open(unit=2345,file="solvers_debug/"//bisect_name,status="replace")

end subroutine open_esweep_file






   subroutine lscheq(nn,ll,mm,ss,epss,u,up,upp)
     implicit none

!
! Adapted from Bachelet
!

     integer, intent(in) :: nn,ll,mm,ss
     real(dp), intent(inout) :: epss
     real(dp), intent(out), dimension(Nx) :: U,Up,Upp

     integer :: id,l1,sls,ir,ii,ix
     integer :: nint,match,node,nodps,nin

     real(dp), dimension(Nx) :: cf
     real(dp), parameter :: accu = 5.0e-8_dp
     real(dp), parameter :: ztol = -1.0e-4_dp
     real(dp) :: emax, emin
     real(dp) :: v0,C2,C3,C4
     real(dp) :: ds,als,de,cap
     real(dp) :: uout,upout
     real(dp) :: Sc,Sn,RO,norm
     real(dp) :: k



     ds = alpha**2 

! accu is the precision for the electron ``eigenvalue''
!
      id=1
      if(ss.lt.0)id=2
      l1=ll+1
      sls=ll*(ll+1)

! Need to do something.
      emax=100.0_dp
      emin=1.0e5_dp
      do 8 ir=1,Nx
         emax=max(emax,Vks(ir,id)+0.5_dp*sls/x(ir)**2-Z/x(ir))
         emin=min(emin,vks(ir,id)+0.5_dp*sls/x(ir)**2-Z/x(ir))
 8    continue
!      write(*,*) nn,ll,mm,ss,emax,emin


!      bound states have to be > emin
       if(epss.lt.0.0d0)then
         if(epss.gt.emax)epss=0.95d0*emax
         if(epss.lt.emin)epss=0.95d0*emin
!         if(epss.gt.emax)epss=0.5d0*(emax+emin)
      else
        epss = 0.5_dp*(emax+emin)
!         if(epss.gt.emax)epss=0.95d0*emax
!         if(epss.lt.emin)epss=1.05d0*emin
!         if(epss.gt.emax)epss=0.5d0*(emax+emin)
      endif
!      if(epss.eq.0.0d0) epss=0.5d0*(emax+emin)
!      write(*,*) nn,ll,ss,epss,emin,emax

      nint=0
      als=alpha**2




!
!  ITERATION STARTS HERE -------------------------
!
      NINT = 0
 10   NINT=NINT+1
!      if (nn==2.and.ll==2.and.mm==0)  then
!      do ir = 1,Nx
!        write(100+Nint,*) x(ix),u(ix)
!      end do
!      end if
      IF(NINT.GT.1000)then
         WRITE(6,999)nn,ll,mm,ss,epss
         write(6,*) 'emin',emin
         write(6,*) 'emax',emax
         Write(6,*) 'Cap',cap
         write(6,*) 'De',de
!         write(6,*) 'ip',ip
         do ir=1,nx
            write(99,*) x(ir),cf(ir),u(ir)
            write(80,*) x(ir),vks(ir,id)
            write(81,*) x(ir),vks(ir,id)+0.5d0*real(sls,dp)/x(ir)**2-Z/x(ir)
            write(82,*) x(ir),cf(ir)
!            write(83,*) r(ir), vxc(id,ir,ip),rho(1,ir,ip),rho(2,ir,ip)
!            write(83,41) r(ir), vxc(id,ir,ip),rho(1,ir,ip),rho(2,ir,ip)
 41   format(4d12.5)
         enddo
 999     FORMAT(4I5,2x,E16.10,' NINT > 1000 ' )
         STOP ' NINT. GT. 1000 '
      end if

 
      do 33 ir=1,nx
        u(ir)=0.d0
        up(ir)=0.d0
        upp(ir)=0.d0
        cf(ir)=0.d0
 33   continue


!
!       DEFINE COEFFICIENT ARRAY FOR U IN DIFF. EQ.
!
        DO ir=1,nx
         CF(Ir)=ALS*(SLS+2.D0*(Vks(ir,id)-z/x(ir)-epss)*x(Ir)**2)
        enddo

!       RADIAL SCH EQUATION READS:
!
!       U'' -AL*U'-CF*U = 0
!
!       FIND CLASSICAL TURNING POINT FOR MATCHING
!
        if ( epss < 0.0_dp ) then

          MATCH=Nx
          do 30 ii=nx,1,-1
            if(cf(ii).lt.0.d0)go to 31
            match=ii
 30       continue
          write(6,*)' no classical matching point ----> stop '
          write(6,*) nn,ll,mm,ss
          write(6,*) 'Epss,emax,emin',epss,emax,emin
          do ir=1,Nx
           write(80,*) x(ir),vks(ir,id)-Z/x(ir)
           write(81,*) x(ir),vks(ir,id)-Z/x(ir)+0.5d0*real(sls,dp)/x(ir)**2
           write(82,*) x(ir),cf(ir)
!           write(83,*) x(ir),vxc(id,ir,ip),rho(1,ir,ip),rho(2,ir,ip)
          enddo

          stop

!       CONTINUUM
        else 

          !this section added - debug
!          write(*,*) nn,ll,mm,ss
!          write(*,*) epss
!          write(*,*) "not a bound state"
          match = log( 0.6_dp*(exp(alpha*real(Nx,dp)) - 1.0_dp )  )/alpha
!          write(*,*) match, x(match)
!           stop


        end if

 31     continue


!        write(6,*) match,r(match),ll
!       START WAVEFUNCTION WITH SERIES
!
!       U=R**(L+1)*(1.D0+C2*R+C3*R**2)
!
! Take care: not exact for all electron computations
! Test.
        v0=vks(1,id)
        C2=-z/(real(ll,dp)+1.d0)
        C3=(v0-epss-z*c2)/real(2*ll+3,dp)
        C4=(c2*(v0-epss) -z*C3)/3.d0/real(ll+4,dp)
!        if(dabs(ylm(ip,ilm)).lt.1.d-10)write(6,*)ilm,ylm(ip,ilm)
        DO 50 Ir=1,4
        U(Ir)=x(Ir)**(ll+1)+C2*x(Ir)**(ll+2)+C3*x(Ir)**(ll+3)&
     &   +C4*x(Ir)**(ll+4)
        UP(Ir)=alpha*((ll+1)*x(Ir)**(ll+1)+C2*(ll+2)*x(Ir)**(ll+2)&
     &  +C3*(ll+3)*x(Ir)**(ll+3) &
     &  +C4*(ll+4)*x(Ir)**(ll+4))
        UPP(Ir)=alpha*UP(Ir)+CF(Ir)*U(Ir)
!        u(ir)=u(ir)*ylm(ip,ilm)
!        up(ir)=up(ir)*ylm(ip,ilm)
!        upp(ir)=upp(ir)*ylm(ip,ilm)
 50     continue
!
!       OUTWARD INTEGRATION USING PREDICTOR ONCE, CORRECTOR TWICE
!
        NODE=0
!
        NODPS=nn-1-ll
        do 70 ir=4,match-1
!
!       PREDICTOR (EXTRAPOLATION)
!
        U(Ir+1)=U(Ir)+(55.d0*UP(Ir)-59.d0*UP(Ir-1)+&
     &   37.d0*UP(Ir-2)-9.d0*UP(Ir-3))/24.d0
        UP(Ir+1)=UP(Ir)+(55.d0*UPP(Ir)-59.d0*UPP(Ir-1)+&
     &   37.d0*UPP(Ir-2)-9.d0*UPP(Ir-3))/24.d0
!
!       CORRECTOR (INTERPOLATION)
!       USED TWICE (IT=1,2) IN ORDER TO GET APPROXIMATE "SELFCONSISTENT"
!       SOLUTION
!
        UPP(Ir+1)=alpha*UP(Ir+1)+CF(Ir+1)*U(Ir+1)
        UP(Ir+1)=UP(Ir)+(9.d0*UPP(Ir+1)+19.d0*UPP(Ir)-&
     &    5.d0*UPP(Ir-1)+UPP(Ir-2))/24.d0
        U(Ir+1)=U(Ir)+(9.d0*UP(Ir+1)+19.d0*UP(Ir)-&
     &    5.d0*UP(Ir-1)+UP(Ir-2))/24.d0

        UPP(Ir+1)=alpha*UP(Ir+1)+CF(Ir+1)*U(Ir+1)
        UP(Ir+1)=UP(Ir)+(9.d0*UPP(Ir+1)+19.d0*UPP(Ir)-&
     &    5.d0*UPP(Ir-1)+UPP(Ir-2))/24.d0
        U(Ir+1)=U(Ir)+(9.d0*UP(Ir+1)+19.d0*UP(Ir)-&
     &    5.d0*UP(Ir-1)+UP(Ir-2))/24.d0
        IF(U(Ir+1)*U(Ir).LE.0)NODE=NODE+1
 70     CONTINUE


        !dealing with bound and continuum
        if ( epss < 0.0_dp ) then

          IF(NODE-NODPS < 0) then
            goto 80
          else if (NODE-NODPS == 0) then
            goto 100
          else
            goto 90
          end if

        else
        
          !for scattering states shoot in without
          !node counting (yet).
          goto 100
    
        end if

!       TOO FEW NODES
!

 80     CONTINUE

!       DEALING WITH BOUND STATE NODES
!
!       IF TOO FEW NODES : INCREASE ENERGY Epss
!
!        if(epss.lt.0.0d0) then
!           EMIN=Epss
!           Epss=0.75D0*Epss
!           IF(Epss.GT.EMAX)Epss=0.5D0*(EMIN+EMAX)
!           if(epss.gt.ztol.and.epss.lt.0.0d0)epss=-epss
!        else
!           EMIN=Epss
!           Epss=1.25D0*Epss
!           IF(Epss.GT.EMAX)Epss=0.5D0*(EMIN+EMAX)
!        endif 
         emin = epss
         epss = 0.5_dp*(emax+emin)  
        goto 10

!       TOO MANY NODES
 90     CONTINUE
!       IF TOO MANY NODES: DECREASE ENERGY E
!        if(epss.lt.0.0d0) then
!           EMAX=Epss
!           Epss=1.25D0*Epss
!           IF(Epss.LT.EMIN)Epss=0.5D0*(EMIN+EMAX)
!        else
!           EMAX=Epss
!           Epss=0.75D0*Epss
!           IF(Epss.LT.EMIN)Epss=0.5D0*(EMIN+EMAX)
!           Epss=0.5D0*(EMIN+EMAX) 
!           if(epss.lt.-ztol.and.epss.gt.0.0d0)epss=-epss
!        endif  
         emax = epss
         epss = 0.5_dp*(emax+emin)  
        goto 10

!       CORRECT NUMBER OF NODES
 100    UOUT=U(Match)
        UPOUT=UP(Match)



!       START INWARD INTEGRATION AT 10*CLASSICAL TURNING POINT
!       WITH SIMPLE EXPONENTIAL   
!                                 
!       ATOMIC BOUNDARY CONDITIONS
!                             
        nin=match+2.3d0/alpha
        do 301 ir=nin+5,nx
          u(ir)=0.d0
          up(ir)=0.d0
          upp(ir)=0.d0
 301    continue

        if ( epss < 0.0_dp) then

          !EXPONENTIAL
          if(nin+4.gt.Nx)nin=Nx-4
!         CAP CHANGE
          CAP=dSQRT(real(SLS,dp)/x(NIN)**2+2.D0*(Vks(NIN,id)-z/x(NIN)-Epss))
!         write(6,*) 'Cap',cap
          DO 110 Ir=NIN,NIN+4

!            original boundary conditions
            U(Ir)=dEXP(-CAP*(x(Ir)-x(NIN)))
            UP(Ir)=-x(Ir)*alpha*CAP*U(Ir)
            UPP(Ir)=alpha*UP(Ir)+CF(Ir)*U(Ir)

!           new boundary conditions   -  cap could go negative in old case
!            U(Ir)=dEXP(-Z*x(Ir)/real(nn,dp))*x(ir)**(nn+1)
!            UP(Ir)=U(Ir)*alpha*(real(nn+1,dp)-Z*x(ir)/real(nn,dp))
!            UPP(Ir)=alpha*UP(Ir)+CF(Ir)*U(Ir)

 110      continue

!          write(*,*)epss,cap,real(SLS,dp)/x(NIN)**2+2.D0*(Vks(NIN,id)-z/x(NIN)-Epss)
        else

          !SINE LIKE AT LARGE DISTANCES
          NIN = Nx-4
          k = sqrt(2.0_dp*epss)
!          write(*,*) k
          DO 111 Ir=NIN,NIN+4
            !U(Ir)=k*( x(ir) - x(Nx) )
            !UP(Ir)=k* x(ir)* alpha
            !UPP(Ir)=alpha*UP(Ir)+CF(Ir)*U(Ir)
            U(Ir)=sin( k*( x(ir) - x(Nx) ) )
            UP(Ir)=k*x(ir)*alpha*cos( k*(x(ir) - x(Nx)) )
            UPP(Ir)=alpha*UP(Ir)+CF(Ir)*U(Ir)
 111      continue
          

        end if
!
!       INTEGRATE INWARD
!       INWARD INTEGRATION USING PREDICTOR ONCE, CORRECTOR TWICE
!
        Ir=NIN+1
 120    Ir=Ir-1
        IF(Ir.LE.Match)GO TO 140
!
!       PREDICTOR
!        AEI=-(55.D0*Y(J)-59.D0*Y(J+1)+
!     $  37.D0*Y(J+2)-9.D0*Y(J+3))/24.0d0
!
        U(Ir-1)=U(Ir)-(55.d0*UP(Ir)-59.d0*UP(Ir+1)+&
     &   37.d0*UP(Ir+2)-9.d0*UP(Ir+3))/24.d0
        UP(Ir-1)=UP(Ir)-(55.d0*UPP(Ir)-59.d0*UPP(Ir+1)+&
     &   37.d0*UPP(Ir+2)-9.d0*UPP(Ir+3))/24.d0
!
!       CORRECTOR
!       USED TWICE (IT=1,2) IN ORDER TO GET APPROXIMATE "SELFCONSISTENT"
!       SOLUTION
!
!        DO IT=1,2
           UPP(Ir-1)=alpha*UP(Ir-1)+CF(Ir-1)*U(Ir-1)
           UP(Ir-1)=UP(Ir)-(9.d0*UPP(Ir-1)+19.d0*UPP(Ir)-&
     &    5.d0*UPP(Ir+1)+UPP(Ir+2))/24.d0
           U(Ir-1)=U(Ir)-(9.d0*UP(Ir-1)+19.d0*UP(Ir)-&
     &    5.d0*UP(Ir+1)+UP(Ir+2))/24.d0
!        enddo
           UPP(Ir-1)=alpha*UP(Ir-1)+CF(Ir-1)*U(Ir-1)
           UP(Ir-1)=UP(Ir)-(9.d0*UPP(Ir-1)+19.d0*UPP(Ir)-&
     &    5.d0*UPP(Ir+1)+UPP(Ir+2))/24.d0
           U(Ir-1)=U(Ir)-(9.d0*UP(Ir-1)+19.d0*UP(Ir)-&
     &    5.d0*UP(Ir+1)+UP(Ir+2))/24.d0
        GO TO 120
 140    CONTINUE
!
!
!       SCALE OUTSIDE WAVEFUNCTION FOR CONTINUITY
!
        SC=UOUT/U(Match)
        DO Ir=Match,NIN+4
           UP(Ir)=SC*UP(Ir)
           UPP(Ir)=SC*UPP(Ir)
           U(Ir)=SC*U(Ir)
        enddo
!
        RO=x(1)/SQRT(exp(alpha))
        SN=RO**(2*ll+3)/real(2*ll+3,dp)+2.D0*C2*RO**(2*ll+4)/real(2*ll+4,dp)&
     &       +(2.D0*C3+C2**2)*RO**(2*ll+5)/real(2*ll+5,dp)
        DO Ir=1,NIN
           SN=SN+alpha*x(Ir)*U(Ir)**2
        enddo



        if ( epss >= 0.0_dp ) then

!
!       NODE COUNTING FOR CONTINUUM STATES TO ADJUST ENERGY
!

          NODE = 0
          do ir = 1,NIN
            IF(U(Ir+1)*U(Ir).LE.0)NODE=NODE+1
          end do
          IF(NODE-NODPS < 0) then
            goto 81
          else if (NODE-NODPS == 0) then
            goto 101
          else
            goto 91
          end if

!
!       IF TOO FEW NODES : INCREASE ENERGY Epss
!
 81     continue
          emin = epss
          epss = 0.5_dp*(emax + emin)
        goto 10
!
!       IF TOO MANY NODES : DECREASE ENERGY Epss
!
 91     continue
          emax = epss
          epss = 0.5_dp*(emax + emin)
        goto 10

        end if

!       CORRECT NODES START FINE TUNING
!     
!       PERTURBATION THEORY PER ENERGY SHIFT
!
!       DE IS THE DIFFERENCE BETWEEN THE LOCAL KINETIC ENERGIES AT R(MCH)
!
!       DE=EKIN(R(MCH),INSIDE)-EKIN(R(MCH),OUTSIDE)
!
!       DE VANISHES FOR THE EXACT SOLUTION
!

 101    continue
        DE=0.5D0*UOUT*(UPOUT-UP(Match))/(SN*alpha*x(Match))
        IF(ABS(DE/Epss).LT.accu) then
           GO TO 170
        else if(ABS(DE).LT.1.d-16) then
           go to 170
        endif
        IF(ABS(DE/Epss).GT.0.25D0)DE=0.25D0*DE*ABS(Epss/DE)
        IF(DE.GT.0.D0)then 
           EMIN=Epss
!           if(emin.gt.ztol.and.emin.lt.0.0d0)emin=-ztol
!           if (nn.eq.1.and.ss.eq.1) write(6,*)'emin',emin 
        endif
        IF(DE.LT.0.D0) then 
           EMAX=Epss
!           if(emax.lt.-ztol.and.emax.gt.0.0d0)emax=-ztol   
!           if (nn.eq.1.and.ss.eq.1) write(6,*) 'emax',emax
        endif
        
!
!       THE TOTAL POTENTIAL IN THE OUTER REGION IS NEARLY CONSTANT
!       THE KINETIC ENERGY OF THE OUTER SOLUTION CHANGES THEREFORE
!       BY NEARLY THE SAME AMOUNT AS THE ENERGY E
!
!       david added
!        if (iter == 5 .and. ip == 9) then
!           write(6,*) epss,nn,ll,mm,ss
!        end if

        Epss=Epss+DE 
        if(epss.gt.ztol.and.epss.lt.0.0d0.and.de.gt.0.0d0) then
           epss=-epss
        endif
        if(epss.lt.-ztol.and.epss.gt.0.0d0.and.de.lt.0.d0) then
           epss=-epss
        endif
        IF(Epss.GT.EMAX.OR.Epss.LT.EMIN) Epss=0.5D0*(EMAX+EMIN)
!        if (nn.eq.1.and.ss.eq.1) write(6,*) 'de',de,epss,emax,emin
        GO TO 10
!
!       END OF ITERATION LOOP
 170   continue


         !normalise the wave functions
         norm = 0.5_dp*(u(1)**2*x(1) + u(Nx)**2*x(Nx))
         do ix = 2,Nx-1
            norm = norm + u(ix)**2*x(ix)
         end do
         norm = norm*alpha*fpi
         norm = 1.0_dp/sqrt(norm)
         do ix = 1,Nx
            u(ix)   = u(ix)*norm
            up(ix)  = up(ix)*norm
            upp(ix) = upp(ix)*norm
         end do


      end subroutine lscheq






!-------------------------------------------!
!      sorting and indexing algorithm       !
!             # Quick Sort #                !
!  Taken from numerical recipes section 8.4 !
!          (tested and working)             !
!-------------------------------------------!
SUBROUTINE indexx(n,arr,indx)
   implicit none


   integer, intent(in) :: n
   integer, intent(inout)  :: indx(n)
   real(dp),intent(in) :: arr(n)
   integer, parameter :: M=7,NSTACK=50
   integer :: i,indxt,ir,itemp,j,jstack,l,istack(NSTACK),kint
   real(dp) ::  a

   do j=1,n
      indx(j)=j
   enddo

   jstack=0
   l=1
   ir=n
   1 if(ir-l.lt.M)then
   do j=l+1,ir
      indxt=indx(j)
      a=arr(indxt)
      do i=j-1,l,-1
      if(arr(indx(i)).le.a)goto 2
         indx(i+1)=indx(i)
      enddo
      i=l-1
      2 indx(i+1)=indxt
   enddo
   if(jstack.eq.0)return
      ir=istack(jstack)
      l=istack(jstack-1)
      jstack=jstack-2
   else
   kint=(l+ir)/2
   itemp=indx(kint)
   indx(kint)=indx(l+1)
   indx(l+1)=itemp
   if(arr(indx(l)).gt.arr(indx(ir)))then
      itemp=indx(l)
      indx(l)=indx(ir)
      indx(ir)=itemp
   endif
   if(arr(indx(l+1)).gt.arr(indx(ir)))then
      itemp=indx(l+1)
      indx(l+1)=indx(ir)
      indx(ir)=itemp
   endif
   if(arr(indx(l)).gt.arr(indx(l+1)))then
      itemp=indx(l)
      indx(l)=indx(l+1)
      indx(l+1)=itemp
   endif
   i=l+1
   j=ir
   indxt=indx(l+1)
   a=arr(indxt)
   3 continue
   i=i+1
   if(arr(indx(i)).lt.a)goto 3
   4 continue
   j=j-1
   if(arr(indx(j)).gt.a)goto 4

   if(j.lt.i)goto 5
      itemp=indx(i)
      indx(i)=indx(j)
      indx(j)=itemp
      goto 3
      5 indx(l+1)=indx(j)
      indx(j)=indxt
      jstack=jstack+2
      if(jstack.gt.NSTACK)then
         write(*,*) "indexing routine, NSTACK too small in indexx"
         stop
      end if
      if(ir-i+1.ge.j-l)then
         istack(jstack)=ir
         istack(jstack-1)=i
         ir=j-1
      else
         istack(jstack)=j-1
         istack(jstack-1)=l
         l=i
      endif
   endif

   goto 1

END SUBROUTINE indexx






end module solvers
