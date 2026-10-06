program mainprogram
!2026/03/21
!symmetric bias version
! global parameters defined in module global_parameters in read_tau.f90    
use hybrid
use global_parameters
use read_tau_module
implicit none
REAL*8 :: TL,TR,TC,Tout(3)
! K_TL2T0, K_TR2T0, K_TLR2TC: within use global_parameters  
!real*8   K_TL2T0    ! thermal conductance between the left electrode and environment at T0
                     ! hould be smaller than Kel(T0) to observe Peltier cooling effect)
!real*8   K_TR2T0    ! thermal conductance between the left electrode and environment at T0
                     ! should be smaller than Kel(T0) to observe Peltier cooling effect)
!real*8   K_TLR2TC   ! take the value of the phonon's thermal conductance of the junction: "Kph" 
                     ! do not include "Kel" because its effect has been included in JL(R)_Peltier 
REAL*8,external :: K0,K1,K2,S,dSdT,SS,tau_Thom,P_Thom,G,Kel,Iel
REAL*8,external :: JR_Peltier, JL_Peltier, J_L2Evn_Fourier, J_R2Evn_Fourier
REAL*8,external :: J_Evn2L_Fourier, J_Evn2R_Fourier
REAL*8,external :: J_R2C_Fourier, J_L2C_Fourier, J_C2R_Fourier, J_C2L_Fourier
!
integer :: I1,I2,I3
REAL*8    Vmin, Vmax,V_max, V_min
REAL*8    alpha_input
REAL*8    V_start, V_end, delta_V, alpha_start, alpha_end, delta_alpha
REAL*8, allocatable :: SSSTL(:,:), SSSTC(:,:), SSSTR(:,:),SSJR(:,:), SSJL(:,:)
REAL*8, allocatable :: alpha(:),VV(:),tempJL(:), tempJR(:), tempTL(:), tempTC(:), tempTR(:)
CHARACTER (LEN =12)::   CH(10)
!!=============================================================================

call read_tau      !Read , calculate and iterpolate tau(E) from tau.inp

!!=============================================================================
!! Set thermal conductance for Fourier's law of heat Conduction:
!! a typical nanojunction with a length of 10 nm and a cross-section of 1 nm^2 (see e.g. Nature 498, 209�V212 (2013))
!! K_TLR2TC=Kph+Kel(T0) is the total thermal conductance of the junction, which includes both phonon and electron contributions.
!! K_TLR2TC=100.d-11+Kel(T0) !(95~105)*10^-11  W/K is the value of the Au phonon's thermal conductance (���s�� LAMPPS)
!!-----------------------------------------------------------------------------
Kph=2.6d-11    ! (2.6~2.8)*10^-11 W/K is the value of the DBDT phonon's thermal conductance (���s�� LAMPPS)
!Kph=100.d-11    !(95~105)*10^-11  W/K is the value of the Au phonon's thermal conductance (���s�� LAMPPS)
K_TLR2TC=Kph+Kel(T0) 

!!-----------------------------------------------------------------------------

! TL,TR,TC: strongly depend on the value of K_TL2T0, which is the rate of heat dissipated to environment.
!! we take K_TL2T0 & K_TR2T0 to be tunable parameters to observe their impact on Peltier cooling effect.


!!=============================================================================
Vb_min=-2.5*dabs(S(300.d0)*300.d0)
Vb_max=2.5*dabs(S(300.d0)*300.d0)
!write(*,"(2(A15,E12.4,3X))") 'Vb_min=', Vb_min, 'Vb_max=', Vb_max
!!
!!=============================================================================
Vb=(S(T0)*T0)/2.d0 !set bias voltage to half of V_th at 300K
TL=T0
TR=T0+0.05d0
TC=T0
write(*,"(5(A5,F12.8,3X))") 'Vb=', Vb, 'TL=', TL, 'TC=', TC, 'TR=', TR, 'T0=', T0
write(*,"(3(A12,E12.4,3X))") 'K0(T=100): ', K0(100.d0), 'K0(T=300): ', K0(300.d0), 'K0(T=500): ', K0(500.d0)
write(*,"(3(A12,E12.4,3X))") 'K1(T=100): ', K1(100.d0), 'K1(T=300): ', K1(300.d0), 'K1(T=500): ', K1(500.d0) 
write(*,"(3(A12,E12.4,3X))") 'K2(T=100): ', K2(100.d0), 'K2(T=300): ', K2(300.d0), 'K2(T=500): ', K2(500.d0) 
write(*,"(3(A12,E12.4,3X))") 'S(T=100): ', S(100.d0), 'S(T=300): ', S(300.d0), 'S(T=500): ', S(500.d0)
write(*,"(3(A12,E12.4,3X))") 'SS(TL,TR): ', SS(TL,TR)
write(*,"(3(A13,E12.4,3X))") 'dSdT(T=100):', dSdT(100.d0), 'dSdT(T=300):', dSdT(300.d0), 'dSdT(T=500):', dSdT(500.d0)
write(*,"(3(A16,E12.4,3X))") 'tau_Thom(T=100):', tau_Thom(100.d0), 'tau_Thom(T=300):', tau_Thom(300.d0), 'tau_Thom(T=500):', tau_Thom(500.d0)
write(*,"(3(A10,E12.4,3X))") 'G(T=100):', G(100.d0), 'G(T=300):', G(300.d0), 'G(T=500):', G(500.d0) 
write(*,"(3(A12,E12.4,3X))") 'Kel(T=100): ', Kel(100.d0), 'Kel(T=300): ', Kel(300.d0), 'Kel(T=500): ', Kel(500.d0)
write(*,"(3(A25,F10.6,3X))") 'V_th=S(T=100K)*T (V)=', S(100.d0)*100.d0, 'V_th=S(T=300K)*T (V)=', S(300.d0)*300.d0,  &
                             'V_th=S(T=500K)*T (V)=', S(500.d0)*500.d0
write(*,"(2(A15,E12.4,3X))") 'Iel(Vb,TL,TR):', Iel(Vb,TL,TR), 'Iel(-Vb,TL,TR):', Iel(-1.d0*Vb,TL,TR)
write(*,"(2(A22,E12.4,3X))") 'JR_Peltier(Vb,TL,TR):', JR_Peltier(Vb,TL,TR),'JR_Peltier(-Vb,TL,TR):', JR_Peltier(-1.d0*Vb,TL,TR)
write(*,"(2(A22,E12.4,3X))") 'JL_Peltier(Vb,TL,TR):', JL_Peltier(Vb,TL,TR), 'JL_Peltier(-Vb,TL,TR):', JL_Peltier(-1.d0*Vb,TL,TR)
write(*,"(2(A22,E12.4,3X))") 'P_Thom(Vb,TL,TC,TR):' , P_Thom(Vb,TL,TC,TR) , 'P_Thom(-Vb,TL,TC,TR):' , P_Thom(-1.d0*Vb,TL,TC,TR) 
write(*,"(2(A30,E12.4,3X))") 'JR_Peltier(Vb=0.3,TL=10,TR=10):', JR_Peltier(.3d0,10.d0,10.d0),'JR_Peltier(-0.3,10,10):', JR_Peltier(-0.3d0,10.d0,10.d0)
write(*,"(2(A30,E12.4,3X))") 'JL_Peltier(Vb=0.3,TL=10,TR=10):', JL_Peltier(.3d0,10.d0,10.d0), 'JL_Peltier(-0.3,10,10):', JL_Peltier(-0.3d0,10.d0,10.d0)
!!------------------------------------------------------------------------------------------------------------
call SeebeckCoefficientOutput           !check ok!20262014! 
call ElectronicThermalConductanceOutput !check ok!20262014! 
call ThomsonCoefficientOutput           !check ok!20262014!
call ElectricConductanceOutput          !check ok!20262014!
call ElectricCurrentOutput              !check ok!20262014!
call JR_Peltier_Output                  !check ok!20262014!
call JL_Peltier_Output                  !check ok!20262014!
call JR_PeltierCoefficient_Output
call JL_PeltierCoefficient_Output
call main_output_Peltier                !check ok!20262021!
call main_output_Thomson                !check ok!20262021!


!!=============================================================================


! !!=============================================================================
! !!  Test region: preparing and testing 
! 777 format(1x, *(g0, ", "))  
!     end program mainprogram
! !!=============================================================================
! !!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
! !!=============================================================================



! subroutine doglegtest
! use myutility; use hybrid
! ! subroutine to test dogleg

! implicit none
! real(kind=db), dimension(3,3) :: jacob, Q, R
! real(kind=db), dimension(3) :: x0, p, Qtf
! external  :: funs2
! real(kind=db), dimension(3) ::  fval
! real(kind=db) :: delta
! integer :: flag
! ! real(kind=db), dimension(2,2) :: GetJacobian

! x0(1) = 0.5_db
! x0(2) = 1.0_db
! x0(3) = 1.5_db
! delta = 0.10_db

! call funs2(x0,fval)
! write(*,*) 'funs2([0.50, 1.00. 1.50 ]):'
! write(*,*) fval

! x0 = x0+1
! call funs2(x0,fval)
! call GetJacobian(jacob,  funs2, x0, 0.0001_db,fval)
! call QRfactorization(jacob,Q,R)

! Qtf = matmul(transpose(Q), fval)
! call dogleg(p,Q,R,delta,Qtf,flag)

! write(*,*) 'p:'
! write(*,*) p
! write(*,*) 'flag:'
! write(*,*) flag

! end subroutine doglegtest



! subroutine finitedifftest()
! use myutility; use hybrid, only : GetJacobian
! ! test finite difference
! implicit none

! real(kind=db), dimension(2,2) :: jacob
! real(kind=db), dimension(2) :: x0
! external  :: funs1
! real(kind=db), dimension(2) ::  fval
! ! real(kind=db), dimension(2,2) :: GetJacobian

! x0(1) = 2.0_db
! x0(2) = 3.0_db

! call funs1(x0, fval)
! call GetJacobian(jacob,  funs1, x0, 0.001_db,fval)

! call MatrixWrite(jacob)

! end subroutine finitedifftest

! subroutine funs1(x, fval0)
! use myutility; 
! implicit none
! real(kind=db), intent(IN), dimension(:) :: x
! real(kind=db), intent(OUT), dimension(:) :: fval0

! fval0(1) = x(2) - x(1)**2
! fval0(2) = 2 - x(1) - x(2)

! end subroutine funs1



! subroutine funs3(x, fval0)
! use myutility; 
! implicit none
! real(kind=db), intent(IN), dimension(:) :: x
! real(kind=db), intent(OUT), dimension(:) :: fval0

! fval0(1) = 10_db*(x(2) - x(1)**2)
! fval0(2) = 2 - x(1) - x(2)

! end subroutine funs3


! subroutine funs2(x, fval0)
! use myutility; 
! ! a little difficult function
! ! solution x = [0.50, 1.00. 1.50 ]  (+ 2pi*n)

! implicit none
! real(kind=db), intent(IN), dimension(:) :: x
! real(kind=db), intent(OUT), dimension(:) :: fval0

! fval0(1) = 1.20_db * sin(x(1)) -1.40_db*cos(x(2))+ 0.70_db*sin(x(3)) &  
! 			- 0.517133908732486_db

! fval0(2) = 0.80_db * cos(x(1)) -0.50_db*sin(x(2))+ 1.00_db*cos(x(3)) &
! 			- 0.352067758776053_db

! fval0(3) = 3.50_db * sin(x(1)) -4.25_db*cos(x(2))+ 2.80_db*cos(x(3))  &
! 			+ 0.4202312501553165_db

! end subroutine funs2


! subroutine qrtest
! use myutility; use hybrid, only : QRfactorization, QRupdate
! ! test QR factorization and update

! implicit none
! real(kind=db), dimension(5,5) ::  Q,  A3
! real(kind=db), dimension(5,5) :: A, R ,A2, A4, A5, A6, B1, B2
! real(kind=db), dimension(5) :: u,v
! INTEGER              :: isize
! INTEGER,ALLOCATABLE  :: iseed(:)

! ! set a seed. 
! CALL RANDOM_SEED(SIZE=isize)
! ALLOCATE( iseed(isize) )
! CALL RANDOM_SEED(GET=iseed)
! iseed = 1
! CALL RANDOM_SEED(PUT = iseed)  

! CALL RANDOM_NUMBER(A)           ! generate  random number
! A = A - 0.5
! write(*,*) 'A:'
! call MatrixWrite( A )

! Q = 1
! R = 1
! call QRfactorization(A,Q,R)

! A2 = matmul(Q , R)
! A3 = matmul(Q , transpose(Q))

! write(*,*) 'Q:'
! call MatrixWrite( Q )
! write(*,*) 'R:'
! call MatrixWrite( R )
! write(*,*) 'Q*R:'
! call MatrixWrite( A2)
! write(*,*) "Q*Q':"
! call MatrixWrite( A3)
! A3 = A2 - A
! write(*,*) "Q*R - A:"
! call MatrixWrite( A3)

! ! update test
! call RANDOM_NUMBER(u)
! call RANDOM_NUMBER(v)
! u = (u-0.5) * 1
! v = (v-0.5) * 1

! write(*,*) ' '
! write(*,*) 'update A '
! call MatrixWrite( A+outer(u,v))
! call QRupdate(Q,R,u,v)


! write(*,*) 'update Q:'
! call MatrixWrite( Q )
! write(*,*) 'update R:'
! call MatrixWrite( R )
! write(*,*) "Q*Q':"
! call MatrixWrite( matmul(Q , transpose(Q)))
! write(*,*) "Q*R - A:"
! call MatrixWrite( matmul(Q , R) - (A + outer(u,v)))


! ! write(*,*) "Q*Q':"
! ! call MatrixWrite( A3)
! ! 

! end subroutine qrtest


! subroutine temp
! ! temp place to put code


! end subroutine temp

