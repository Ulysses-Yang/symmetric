
    !!===================================================   
subroutine SeebeckCoefficientOutput
implicit none
REAL*8 :: T_begin, T_end, deltaT ,T ,temp 
integer :: NT,I2
REAL*8, allocatable :: TT(:), SSS(:),SdSdT(:),Stau_Thom(:) 
real*8,external :: S,dSdT,tau_Thom 
open(30,file='S(T).csv')
777 format(1x, *(g0, ", ")) 
T_begin=100.d0
T_end=500.d0
deltaT=1.d0
NT=NINT((T_end-T_begin)/deltaT)+1
!
allocate( TT(NT+1),SSS(NT+1),SdSdT(NT+1),Stau_Thom(NT+1)  )
do I2=1,NT 
    T=T_begin+(I2-1)*deltaT 
    TT(I2)=T
    SSS(I2)=S(T)
    SdSdT(I2)=dSdT(T)
    Stau_Thom(I2)=tau_Thom(T)   
enddo
!
write(30,777) 'T(K)', 'S(T)(V/K)', 'S(T)(microV/K)', 'dS/dT (microV/K^2)', 'tau_Thom(T)(microV/K)'
!
do I2=1,NT
    write(30,777) TT(I2), SSS(I2), SSS(I2)*1.d6, SdSdT(I2)*1.d9, Stau_Thom(I2)*1.d9 
enddo
close(30)
end subroutine SeebeckCoefficientOutput
      
    
!!===================================================   
subroutine ElectronicThermalConductanceOutput
implicit none
REAL*8 :: T_begin, T_end, deltaT ,T ,temp 
integer :: NT,I2
REAL*8, allocatable :: TT(:), SSS(:) 
real*8,external :: Kel
open(30,file='Kel(T).csv')
777 format(1x, *(g0, ", ")) 
T_begin=100.d0
T_end=500.d0
deltaT=1.d0
NT=NINT((T_end-T_begin)/deltaT)+1
!
allocate(TT(NT+1),SSS(NT+1))
do I2=1,NT 
    T=T_begin+(I2-1)*deltaT 
    TT(I2)=T
    SSS(I2)=Kel(T)
enddo
!
write(30,777) 'T(K)', 'Kel(T)(W/K)', 'Kel(T)(nW/K)'
!
do I2=1,NT
    write(30,777) TT(I2), SSS(I2), SSS(I2)*1.d9
enddo
close(30)
    end subroutine ElectronicThermalConductanceOutput 
    
!!===================================================   
subroutine ThomsonCoefficientOutput
implicit none
REAL*8 :: T_begin, T_end, deltaT ,T ,temp 
integer :: NT,I2
REAL*8, allocatable :: TT(:), SSS(:) 
real*8,external :: tau_Thom
open(30,file='tau_Thom(T).csv')
777 format(1x, *(g0, ", ")) 
T_begin=100.d0
T_end=500.d0
deltaT=1.d0
NT=NINT((T_end-T_begin)/deltaT)+1
!
allocate(TT(NT+1),SSS(NT+1))
do I2=1,NT 
    T=T_begin+(I2-1)*deltaT 
    TT(I2)=T
    SSS(I2)=tau_Thom(T)
enddo
!
write(30,777) 'T(K)', 'tau_Thom(T)(V/K)', 'tau_Thom(T)(microV/K)'
!
do I2=1,NT
    write(30,777) TT(I2), SSS(I2), SSS(I2)*1.d9
enddo
close(30)
    end subroutine ThomsonCoefficientOutput     
!!===================================================   
subroutine ElectricConductanceOutput
implicit none
REAL*8 :: T_begin, T_end, deltaT ,T ,temp 
integer :: NT,I2
REAL*8, allocatable :: TT(:), SSS(:) 
real*8,external :: G
open(30,file='G(T).csv')
777 format(1x, *(g0, ", ")) 
T_begin=100.d0
T_end=500.d0
deltaT=1.d0
NT=NINT((T_end-T_begin)/deltaT)+1
!
allocate(TT(NT+1),SSS(NT+1))
do I2=1,NT 
    T=T_begin+(I2-1)*deltaT 
    TT(I2)=T
    SSS(I2)=G(T)
enddo
!
write(30,777) 'T(K)', 'G(T)(S;siemens)', 'G(T)(microS)'
!
do I2=1,NT
    write(30,777) TT(I2), SSS(I2), SSS(I2)*1.d6
enddo
close(30)
    end subroutine ElectricConductanceOutput         
!!===================================================   
subroutine ElectricCurrentOutput
use global_parameters
implicit none
REAL*8 :: T_begin,V_begin, T_end,V_end, deltaT,deltaV,T,V,temp 
integer :: NT,NV,I2,J2
REAL*8, allocatable :: TT(:),VV(:), SSS(:,:) 
real*8,external :: Iel
open(30,file='I(T,Vb).csv')
777 format(1x, *(g0, ", "))     
T_begin=100.d0
T_end=500.d0
deltaT=50.d0
NT=NINT((T_end-T_begin)/deltaT)+1
!   
V_begin=Vb_min
V_end=Vb_max
NV=200
deltaV=(V_end-V_begin)/(NV-1)
!
allocate(TT(NT+1),VV(NV+1),SSS(NT+1,NV+1))
do I2=1,NT 
    T=T_begin+(I2-1)*deltaT 
    TT(I2)=T
enddo
!
do I2=1,NV 
    V=V_begin+(I2-1)*deltaV 
    VV(I2)=V
enddo
!write(*,*)(TT(I2), I2=1,NT)
!write(*,*)(VV(I2), I2=1,NV)
!write(*,*) 'Vb_max= ', Vb_max, 'Vb_min= ', Vb_min
do I2=1,NT 
    T=TT(I2)
    do J2=1,NV
        V=VV(J2)
        temp=Iel(V,T,T)
        SSS(I2,J2)=temp  ! Iel(T) at TL=TR=T
    enddo
enddo
!write(*,*) 'SSS(1,1)=',SSS(1,1), 'SSS(NT,NV)=', SSS(NT,NV)
write(30,777) 'T(K)', 'V (V)', 'Iel(T;V) (A)'
do I2=1,NT
    do J2=1,NV
    write(30,777) TT(I2), VV(J2),  SSS(I2,J2)
    enddo
enddo
close(30)
    end subroutine ElectricCurrentOutput          

!!===================================================   
subroutine JR_Peltier_Output
use global_parameters
implicit none
REAL*8 :: T_begin,V_begin, T_end,V_end, deltaT,deltaV,T,V,temp 
integer :: NT,NV,I2,J2
REAL*8, allocatable :: TT(:),VV(:), SSS(:,:) 
real*8,external :: JR_Peltier
open(30,file='JR_Peltier(T,Vb).csv')
777 format(1x, *(g0, ", "))     
T_begin=100.d0
T_end=500.d0
deltaT=50.d0
NT=NINT((T_end-T_begin)/deltaT)+1
!   
V_begin=Vb_min
V_end=Vb_max
NV=200
deltaV=(V_end-V_begin)/(NV-1)
!
allocate(TT(NT+1),VV(NV+1),SSS(NT+1,NV+1))
do I2=1,NT 
    T=T_begin+(I2-1)*deltaT 
    TT(I2)=T
enddo
!
do I2=1,NV 
    V=V_begin+(I2-1)*deltaV 
    VV(I2)=V
enddo
!write(*,*)(TT(I2), I2=1,NT)
!write(*,*)(VV(I2), I2=1,NV)
!write(*,*) 'Vb_max= ', Vb_max, 'Vb_min= ', Vb_min
do I2=1,NT 
    T=TT(I2)
    do J2=1,NV
        V=VV(J2)
        temp=JR_Peltier(V,T,T)
        SSS(I2,J2)=temp  ! Iel(T) at TL=TR=T
    enddo
enddo
!write(*,*) 'SSS(1,1)=',SSS(1,1), 'SSS(NT,NV)=', SSS(NT,NV)
write(30,777) 'T(K)', 'V (V)', 'JR_Peltier(T;V) (W/s)', 'V (mV)', 'JR_Peltier(T;V) (nW/s)'
do I2=1,NT
    do J2=1,NV
        write(30,777) TT(I2), VV(J2),  SSS(I2,J2), VV(J2)*1.d3, SSS(I2,J2)*1.d9
    enddo
enddo
close(30)
    end subroutine JR_Peltier_Output   
    
!!===================================================   
subroutine JL_Peltier_Output
use global_parameters
implicit none
REAL*8 :: T_begin,V_begin, T_end,V_end, deltaT,deltaV,T,V,temp 
integer :: NT,NV,I2,J2
REAL*8, allocatable :: TT(:),VV(:), SSS(:,:) 
real*8,external :: JL_Peltier
open(30,file='JL_Peltier(T,Vb).csv')
777 format(1x, *(g0, ", "))     
T_begin=100.d0
T_end=500.d0
deltaT=50.d0
NT=NINT((T_end-T_begin)/deltaT)+1
!   
V_begin=Vb_min
V_end=Vb_max
NV=200
deltaV=(V_end-V_begin)/(NV-1)
!
allocate(TT(NT+1),VV(NV+1),SSS(NT+1,NV+1))
do I2=1,NT 
    T=T_begin+(I2-1)*deltaT 
    TT(I2)=T
enddo
!
do I2=1,NV 
    V=V_begin+(I2-1)*deltaV 
    VV(I2)=V
enddo
!write(*,*)(TT(I2), I2=1,NT)
!write(*,*)(VV(I2), I2=1,NV)
!write(*,*) 'Vb_max= ', Vb_max, 'Vb_min= ', Vb_min
do I2=1,NT 
    T=TT(I2)
    do J2=1,NV
        V=VV(J2)
        temp=JL_Peltier(V,T,T)
        SSS(I2,J2)=temp  ! Iel(T) at TL=TR=T
    enddo
enddo
!write(*,*) 'SSS(1,1)=',SSS(1,1), 'SSS(NT,NV)=', SSS(NT,NV)
write(30,777) 'T(K)', 'V (V)', 'JL_Peltier(T;V) (W/s)', 'V (mV)', 'JL_Peltier(T;V) (nW/s)'
do I2=1,NT
    do J2=1,NV
        write(30,777) TT(I2), VV(J2),  SSS(I2,J2), VV(J2)*1.d3, SSS(I2,J2)*1.d9
    enddo
enddo
close(30)
    end subroutine JL_Peltier_Output   
    

!!==============Peltier Coefficient PI_R=========================   
!! PI_R == | JR_Peltier/I_el |    , if JR_Peltier<0    
subroutine JR_PeltierCoefficient_Output
use global_parameters
implicit none
REAL*8 :: T_begin,V_begin, T_end,V_end, deltaT,deltaV,T,V,temp 
integer :: NT,NV,I2,J2
REAL*8, allocatable :: TT(:),VV(:), SSS(:,:) 
real*8,external :: JR_Peltier,Iel
open(30,file='JR_PeltierCoefficient(T,Vb).csv')
777 format(1x, *(g0, ", "))     
T_begin=100.d0
T_end=500.d0
deltaT=50.d0
NT=NINT((T_end-T_begin)/deltaT)+1
!   
V_begin=Vb_min
V_end=Vb_max
NV=200
deltaV=(V_end-V_begin)/(NV-1)
!
allocate(TT(NT+1),VV(NV+1),SSS(NT+1,NV+1))
SSS(:,:)=0.d0 !set initial values equal 0
do I2=1,NT 
    T=T_begin+(I2-1)*deltaT 
    TT(I2)=T
enddo
!
do I2=1,NV 
    V=V_begin+(I2-1)*deltaV 
    VV(I2)=V
enddo
!write(*,*)(TT(I2), I2=1,NT)
!write(*,*)(VV(I2), I2=1,NV)
!write(*,*) 'Vb_max= ', Vb_max, 'Vb_min= ', Vb_min
do I2=1,NT 
    T=TT(I2)
    do J2=1,NV
        V=VV(J2)
        temp=JR_Peltier(V,T,T) ! JR_Peltier(V,T,T) at TL=TR=T
        if (temp<0.d0) then
        SSS(I2,J2)=dabs(temp/Iel(V,T,T))  ! Iel(V,T,T) at TL=TR=T
        endif
    enddo
enddo
!write(*,*) 'SSS(1,1)=',SSS(1,1), 'SSS(NT,NV)=', SSS(NT,NV)
write(30,777) 'T(K)', 'V (V)', 'JR_Peltier_Coefficient=JR_Peltier/I (V)', 'V (mV)', 'JR_Peltier_Coefficient (microV)'
do I2=1,NT
    do J2=1,NV
        write(30,777) TT(I2), VV(J2),  SSS(I2,J2), VV(J2)*1.d3, SSS(I2,J2)*1.d6
    enddo
enddo
close(30)
    end subroutine JR_PeltierCoefficient_Output       
    
    
!!==============Peltier Coefficient PI_L=========================   
!! PI_L == | JL_Peltier/I_el |    , if JL_Peltier<0    
subroutine JL_PeltierCoefficient_Output
use global_parameters
implicit none
REAL*8 :: T_begin,V_begin, T_end,V_end, deltaT,deltaV,T,V,temp 
integer :: NT,NV,I2,J2
REAL*8, allocatable :: TT(:),VV(:), SSS(:,:) 
real*8,external :: JL_Peltier,Iel
open(30,file='JL_PeltierCoefficient(T,Vb).csv')
777 format(1x, *(g0, ", "))     
T_begin=100.d0
T_end=500.d0
deltaT=50.d0
NT=NINT((T_end-T_begin)/deltaT)+1
!   
V_begin=Vb_min
V_end=Vb_max
NV=200
deltaV=(V_end-V_begin)/(NV-1)
!
allocate(TT(NT+1),VV(NV+1),SSS(NT+1,NV+1))
SSS(:,:)=0.d0 !set initial values equal 0
do I2=1,NT 
    T=T_begin+(I2-1)*deltaT 
    TT(I2)=T
enddo
!
do I2=1,NV 
    V=V_begin+(I2-1)*deltaV 
    VV(I2)=V
enddo
!write(*,*)(TT(I2), I2=1,NT)
!write(*,*)(VV(I2), I2=1,NV)
!write(*,*) 'Vb_max= ', Vb_max, 'Vb_min= ', Vb_min
do I2=1,NT 
    T=TT(I2)
    do J2=1,NV
        V=VV(J2)
        temp=JL_Peltier(V,T,T) ! JL_Peltier(V,T,T) at TL=TR=T
        if (temp<0.d0) then
        SSS(I2,J2)=dabs(temp/Iel(V,T,T))  ! Iel(V,T,T) at TL=TR=T
        endif
    enddo
enddo
!write(*,*) 'SSS(1,1)=',SSS(1,1), 'SSS(NT,NV)=', SSS(NT,NV)
write(30,777) 'T(K)', 'V (V)', 'JL_Peltier_Coefficient=JL_Peltier/I (V)', 'V (mV)', 'JL_Peltier_Coefficient (microV)'
do I2=1,NT
    do J2=1,NV
        write(30,777) TT(I2), VV(J2),  SSS(I2,J2), VV(J2)*1.d3, SSS(I2,J2)*1.d6
    enddo
enddo
close(30)
    end subroutine JL_PeltierCoefficient_Output   
    
    
    
    
    
    
    
    
    
    