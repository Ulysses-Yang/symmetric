
    subroutine  main_output_Thomson
use hybrid
use global_parameters
use read_tau_module
implicit none    
INTEGER :: I1,I2,I3,I4
REAL*8    TLTCTR(3)
REAL*8   TL, TC, TR, TL2, TC2, TR2
REAL*8    Vmin, Vmax,V_max, V_min
REAL*8    alpha_input
REAL*8    V_start, V_end, delta_V, alpha_start, alpha_end, delta_alpha
REAL*8, allocatable :: SS1JR(:,:), SS1JL(:,:),SS1P_Thom(:,:),SS1IelVb(:,:)
REAL*8, allocatable :: alpha(:),VV(:)
CHARACTER (LEN =12) ::   CH(1:10)
REAL*8,external :: K0,K1,K2,S,dSdT,SS,tau_Thom,P_Thom,G,Kel,Iel
REAL*8,external :: JR_Peltier, JL_Peltier, J_L2Evn_Fourier, J_R2Evn_Fourier
REAL*8,external :: J_Evn2L_Fourier, J_Evn2R_Fourier
REAL*8,external :: J_R2C_Fourier_noKel, J_L2C_Fourier_noKel, J_C2R_Fourier, J_C2L_Fourier
!   
!!! output parameters for Peltier cooling effect    
!    N_Vb=201         !! read from global_parameter, number of bias voltage points for testing Peltier cooling effect
!    N_alpha=41       !!read from global_parameter, number of alpha=K_Fourier^{L(R)2Evronment} 
!!(thermal conductance to environment) points for testing Peltier cooling effect
!!
allocate( VV(N_Vb), alpha(N_alpha) )
allocate(  SS1JL(N_alpha,N_Vb), SS1JR(N_alpha,N_Vb), SS1P_Thom(N_alpha,N_Vb), SS1IelVb(N_alpha,N_Vb))             
!
V_max=2.5d0*dabs(S(300.d0)*300.d0)
V_start=-1.0*dabs(V_max)
V_end=1.0*dabs(V_max)
delta_V=(V_end-V_start)/DFLOAT(N_Vb-1)
do I2=1,N_Vb
    VV(I2)=V_start+(I2-1)*delta_V
enddo
alpha_start=1.d-12
alpha_end=1.d-7
!! take log_10 of alpha to have a better resolution for small alpha, 
!! which is more relevant for observing Peltier cooling effect    
alpha_start=DLOG10(alpha_start)
alpha_end=DLOG10(alpha_end)
delta_alpha=(alpha_end-alpha_start)/DFLOAT(N_alpha-1)
!!
do I1=1,N_alpha
    alpha(I1)=alpha_start+(I1-1)*delta_alpha
enddo
alpha(:)=10.0d0**alpha(:)  ! take log_10 of alpha to have a better resolution for small alpha, 
                           ! which is more relevant for observing Peltier cooling effect
!!-----------------------------------------------------------------------------
!alpha are set to be the thermal conductance which dissipate heat to environment from the left and right electrode,
! which is treated as a tunable parameter to observe its impact on Peltier cooling effect. set alpha as,
do I1=1, N_alpha    !do loop for different alpha (thermal conductance to environment)  
        alpha_input=alpha(I1)
        K_TL2T0=alpha_input   !treated as a parameter
        K_TR2T0=alpha_input   !treated as a parameter    
    do I2=1,N_Vb  !do loop for the applied bias voltage Vb
       Vb=VV(I2)
                      TL_guess=T0
                      TC_guess=T0
                      TR_guess=T0          
!!---------------------------------------------------         
         call FsolveTLTCTR_Thomson(TLTCTR)
!!---------------------------------------------------         
       TL=TLTCTR(1)
       TC=TLTCTR(2)
       TR=TLTCTR(3)
       SS1TL(I1,I2)=TL
       SS1TC(I1,I2)=TC
       SS1TR(I1,I2)=TR 
       SS1JL(I1,I2)=JL_Peltier(Vb,TL,TR)
       SS1JR(I1,I2)=JR_Peltier(Vb,TL,TR) 
       SS1P_Thom(I1,I2)=P_Thom(Vb,TL,TC,TR)
       SS1IelVb(I1,I2)=Iel(Vb,TL,TR)*Vb
    enddo !enddo loop for the applied bias voltage Vb
!!------------------------------------------    
enddo !enddo loop for different alpha (thermal conductance to environment)

!!====================================================================================================
!!     Output data (Thomson effect) for alpha, Vb, TL, TC, TR, JL_Peltier, JR_Peltier
!!====================================================================================================
!! output main date for testing: alpha (thermal conductance to environment), Vb, TL, TC, TR, JL_Peltier, JR_Peltier
!! in unit of nW/K for alpha, mV for Vb, K for TL, TC, and TR, and nW for JL_Peltier and JR_Peltier 
open(20,file='TLTCTR_Thomson_2D.csv')
CH(1)='K_LR2T0'; CH(2)='Vb'; CH(3)='TL'; CH(4)='TC'; CH(5)='TR';CH(6)= 'JL_Peltier'
CH(7)='JR_Peltier'; CH(8)='P_Thomson';CH(9)='JL+JR';CH(10)='P=Iel*Vb' 
!
write(20,777) ( CH(1)  ,  CH(2)  ,   CH(3), CH(4) , CH(5) ,  CH(6) ,  CH(7) ,  CH(8) , CH(9) ,  CH(10) ,I1=1,N_alpha)
write(20,777) ( 'nW/K' ,  'mV'   ,   'K'  ,  'K'  ,  'K'  ,  'nW'  ,  'nW'  ,  'nW'  ,  'nW'  ,  'nW'  ,I1=1,N_alpha)
write(20,777) ('alpha=','',alpha(I1)*1.d9 , 'K_LR2T0=','',alpha(I1)*1.d9,'' ,   ''   ,   ''   ,   ''   ,I1=1,N_alpha)
!
do I2=1,N_Vb
write(20,777) (alpha(I1)*1.d9, VV(I2)*1.d3, SS1TL(I1,I2), SS1TC(I1,I2), SS1TR(I1,I2), &
               SS1JL(I1,I2)*1.d9, SS1JR(I1,I2)*1.d9 , SS1P_Thom(I1,I2)*1.d9,          &
               (SS1JL(I1,I2)+SS1JR(I1,I2))*1.d9 , SS1IelVb(I1,I2)*1.d9, I1=1,N_alpha)
enddo  !do I2=1,N_Vb  
close(20)
!!-------------------------------------------------------------------------------------------
open(21,file='Delta_TLTCTR_Thomson_2D.csv')
CH(1)='K_LR2T0'; CH(2)='Vb'; CH(3)='Delta TL'; CH(4)='Delta TC'; CH(5)='Delta TR'; CH(6)= 'JL_Peltier'
CH(7)='JR_Peltier'; CH(8)='P_Thomson';CH(9)='JL+JR';CH(10)='P=Iel*Vb'
write(21,777) ( CH(1)  ,  CH(2)  ,  CH(3) , CH(4) , CH(5) ,  CH(6) ,  CH(7) , CH(8)  , CH(9) ,  CH(10) ,I1=1,N_alpha)
write(21,777) ( 'nW/K' ,  'mV'   ,   'K'  ,  'K'  ,  'K'  ,  'nW'  ,  'nW'  , 'nW'   ,  'nW' ,  'nW'   ,I1=1,N_alpha)
write(21,777) ('alpha=','',alpha(I1)*1.d9,'','Delta T=T-T0','','','','','', I1=1,N_alpha)
!
do I2=1,N_Vb
write(21,777) (alpha(I1)*1.d9, VV(I2)*1.d3, SS1TL(I1,I2)-T0, SS1TC(I1,I2)-T0, SS1TR(I1,I2)-T0, &
               SS1JL(I1,I2)*1.d9, SS1JR(I1,I2)*1.d9 , SS1P_Thom(I1,I2)*1.d9,  &
               (SS1JL(I1,I2)+SS1JR(I1,I2))*1.d9 , SS1IelVb(I1,I2)*1.d9 , I1=1,N_alpha)
enddo  !do I2=1,N_Vb 
close(21)
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
open(30,file='TLTCTR_Thomson_3D.csv')
CH(1)='K_LR2T0'; CH(2)='Vb'; CH(3)='TL'; CH(4)='TC'; CH(5)='TR';CH(6)= 'JL_Peltier'
CH(7)='JR_Peltier'; CH(8)='P_Thomson';CH(9)='JL+JR';CH(10)='P=Iel*Vb'
!
write(30,777) CH(1)  ,  CH(2)  ,  CH(3) , CH(4) , CH(5) ,  CH(6) ,  CH(7) , CH(8), CH(9) ,  CH(10)
write(30,777) 'nW/K' ,  'mV'   ,   'K'  ,  'K'  ,  'K'  ,  'nW'  ,  'nW'  , 'nW' ,  'nW' ,  'nW'  
write(30,777) 'alpha=K_LR2T0','','','','','','','','',''
!
do I1=1,N_alpha
   do I2=1,N_Vb
      write(30,777) alpha(I1)*1.d9, VV(I2)*1.d3, SS1TL(I1,I2), SS1TC(I1,I2), SS1TR(I1,I2), &
          SS1JL(I1,I2)*1.d9, SS1JR(I1,I2)*1.d9 , SS1P_Thom(I1,I2)*1.d9,  &
          (SS1JL(I1,I2)+SS1JR(I1,I2))*1.d9 , SS1IelVb(I1,I2)*1.d9
   enddo  !do I2=1,N_Vb 
enddo  !do I1=1,N_alpha 
close(30)
!!----------------------------------------------------------------------------------------------
open(31,file='Delta_TLTCTR_Thomson_3D.csv')
CH(1)='K_LR2T0'; CH(2)='Vb'; CH(3)='TL'; CH(4)='TC'; CH(5)='TR'; CH(6)= 'JL_Peltier'
CH(7)='JR_Peltier'; CH(8)='P_Thomson'; CH(9)='JL+JR'; CH(10)='P=Iel*Vb'
!
write(31,777) CH(1)  ,  CH(2)  ,  CH(3) , CH(4) , CH(5) ,  CH(6) ,  CH(7) , CH(8), CH(9) ,  CH(10)
write(31,777) 'nW/K' ,  'mV'   ,   'K'  ,  'K'  ,  'K'  ,  'nW'  ,  'nW'  , 'nW' ,  'nW' ,  'nW'
write(31,777) 'Delta T=T-T0','','','','','','','','',''
!
do I1=1,N_alpha
   do I2=1,N_Vb
      write(31,777) alpha(I1)*1.d9, VV(I2)*1.d3, SS1TL(I1,I2)-T0, SS1TC(I1,I2)-T0, SS1TR(I1,I2)-T0, &
          SS1JL(I1,I2)*1.d9, SS1JR(I1,I2)*1.d9, SS1P_Thom(I1,I2)*1.d9, &
          (SS1JL(I1,I2)+SS1JR(I1,I2))*1.d9 , SS1IelVb(I1,I2)*1.d9
   enddo  !do I2=1,N_Vb 
enddo  !do I1=1,N_alpha 
close(31)

777 format(1x, *(g0, ", "))  
    endsubroutine   main_output_Thomson   
!!==========================================================================================    
!!==========================================================================================    
subroutine FsolveTLTCTR_Thomson(Tout)
	! Solev temperatures TL, TC, and TR in the nanojunction: Thomson effect
    use global_parameters
    use MyUtility
    use hybrid
	implicit none
	real*8, dimension(3,3) :: jacob, Q, R
	real*8, dimension(3) :: xout2, fval,Tout
integer :: fsolveinfo
integer, save :: solve_count = 0

solve_count = solve_count + 1

write(*,'(A,I0)') 'FsolveHybrid ', solve_count

	call FsolveHybrid( &
		fun          = heatmodel, &         ! Function to be solved
        x0           =(/TL_guess, TC_guess, TR_guess/), &       ! Initial value
		xout         =xout2, &              ! output 
		xtol         =1.d-14, &              ! error torelance
		info         =fsolveinfo, &         ! info for the solution
		fvalout      =fval,  &              ! f(xout)
		JacobianOut  =jacob, &              ! Jacobian at x = xout
		JacobianStep =1.d-9,  &             ! Stepsize for the Jacobian
		display      =0,  &                 ! Control for the display
		MaxFunCall   = 1000, &              ! Max number of function call
		factor       =1.0d-9, &             ! Initial value of delta
		NoUpdate     = 0)                   ! control for update of Jacobian

!	write(*,*) ' '
!	write(*,*) 'Solution:'
	! call VectorWrite(xout2)
	! write(*,*) ' '
	! write(*,*) 'Function Value at the solution:'
	! call VectorWrite(fval)
    ! write(*,"(5(A12,E12.4))") 'T0=', T0, 'Vb=', Vb,'K_TLR2TC=', K_TLR2TC, 'K_TL2T0=', K_TL2T0, 'K_TR2T0=', K_TR2T0
    ! Tout(:)=xout2(:)	
!!    
	contains
		subroutine heatmodel(x, fval0)
        use global_parameters
        ! heating model for the nanojunction: solve TL, TC, and TR in the nanojunction
        ! solution x = [x(1), x(2), x(3) ]  =[TL, TC, TR]
		implicit none
		real*8, intent(IN), dimension(:) :: x
		real*8, intent(OUT), dimension(:) :: fval0
        real*8 Vb_in
        REAL*8,external :: JR_Peltier, JL_Peltier, J_L2Evn_Fourier, J_R2Evn_Fourier
        REAL*8,external :: J_Evn2L_Fourier, J_Evn2R_Fourier
        REAL*8,external :: J_R2C_Fourier_noKel, J_L2C_Fourier_noKel, J_C2R_Fourier, J_C2L_Fourier
        REAL*8,external :: P_Thom       
        ! P_Thom(Vb_in,TL_in,TC_in,TR_in)=-1.d0 * tau_Thom(T) * I * \Delta T
        !JR_Peltier(Vb_in,TL_in,TR_in); JL_Peltier(Vb_in,TL_in,TR_in);
        !JL2Evn_Fourier(TLin) ; JR2Evn_Fourier(TRin) 
        ! JR2C_Fourier(TRin,TCin) ; JL2C_Fourier(TLin,TCin)  
        !JR2C_Fourier(TRin,TCin) ; JC2R_Fourier(TCin,TRin) 
        !all the thermal currents in SI unit of W/s; 
        !To enhance the magnitude, we will multiply all the thermal currents in the heat model by 1.d14 to convert them into nW/s
        Vb_in=Vb
        
        fval0(1) = (JL_Peltier(Vb_in,x(1),x(3)) + J_Evn2L_Fourier(x(1)) + J_C2L_Fourier(x(2),x(1)))*1.d14
        
        fval0(2) = (J_L2C_Fourier_noKel(x(1),x(2)) + J_R2C_Fourier_noKel(x(3),x(2)) +P_Thom(Vb_in, x(1), x(2), x(3)))*1.d14
        

        fval0(3) = (JR_Peltier(Vb_in,x(1),x(3)) + J_Evn2R_Fourier(x(3)) + J_C2R_Fourier(x(2),x(3)))*1.d14
        
		end subroutine heatmodel
    end subroutine FsolveTLTCTR_Thomson