
    !!===========================================================================
        
function S(T_in) 
use global_parameters
real*8, intent(in) :: T_in
REAL*8,external :: K0,K1,K2 
REAL*8 S,T,A   
T=T_in
A=-1.d0/(ee*T)
S=A*K1(T)/K0(T)
return
end function  S
    
!!===========================================================================
!! First derivative of Seebeck coefficient S(T) at temperature T with the formula: 
!! S'(T)=dS/dT=( S(T + h) - S(T - h) )/(2*h) with h=1.d-2 K
function dSdT(T_in) 
real*8, intent(in) :: T_in
REAL*8,external :: S
REAL*8 dSdT,T,dT  
T=T_in
dT=1.D-2
dSdT=( S(T+dT) - S(T-dT) )/(2.d0*dT)
return
    end function  dSdT
    
!!===========================================================================
!! Thomson coefficient tau_Thom(T)=T*[dS/dT(T)] 
function tau_Thom(T_in) 
real*8, intent(in) :: T_in
REAL*8,external :: S,dSdT
REAL*8 tau_Thom,T,dT  
T=T_in
tau_Thom = T * dSdT(T)
return
end function  tau_Thom 
    
!!===========================================================================
!! generalization of Seebeck coefficient S(T) at temperature T with the formula: S=-1/(ee*T)*K1/K0
!Seebeck coefficient S(TL,TR) of junction at temperature TL and TR
function SS(TL_in,TR_in)
use global_parameters
real*8, intent(in) :: TL_in,TR_in
REAL*8,external :: K0,K1,K2 
REAL*8 SS,A,TL,TR  
TL=TL_in
TR=TR_in
A=-1.d0/ee
SS=A*(   K1(TL)/TL + K1(TR)/TR   )/(   K0(TL) + K0(TR)   )
return
end function  SS    
    
!!===========================================================================        
function G(T_in) 
use global_parameters
real*8, intent(in) :: T_in
REAL*8,external :: K0,K1,K2 
REAL*8 G,T,G0
T=T_in
G0=2.d0*ee*ee/h
G=G0*K0(T)
return
end function  G
!!===========================================================================
        
function Kel(T_in) 
use global_parameters
real*8, intent(in) :: T_in
REAL*8,external :: K0,K1,K2 
REAL*8 Kel,T,A
T=T_in 
A=2.d0/(h*T)
Kel=A*(K2(T)*K0(T)-K1(T)**2)/K0(T)
return
    end function  Kel 

!!===========================================================================
    
function K0(T_in)
use global_parameters
real*8, intent(in) :: T_in
INTEGER:: MSTART=1
integer  NE
REAL*8 dfede,fe,EXPON,DX,K0,T
REAL*8, allocatable:: dFK0(:),X(:),F(:),Z(:),Q(:),R(:),S(:)
T=T_in
NE=NE_0
allocate(dFK0(NE),X(NE),F(NE),Z(NE),Q(NE),R(NE),S(NE))
DO I=1,NE
EXPON=DMin1((EE_0(I)-EF)/(AKB*T),300.d0) !to avoid overflow for exponential at low temperature
fe=1.d0/(DEXP(EXPON)+1.d0) !Fermi-Dirac distribution f
dfede=-1.d0*fe*(1.d0-fe)/(AKB*T) !df(E)/dE
dFK0(I)=-1.d0*dfede*tau_0(I)
!dFK1(I)=-1.d0*((EE_0(I)-EF)**1)*dfede**tau_0(I)
!dFK2(I)=-1.d0*((EE_0(I)-EF)**2)*dfede**tau_0(I)
END DO
!write(*,"(I5,2E16.6)") ((I,EE_0(I),dFK0(I)),I=1,NE)  
!write(*,"(I5,2E16.6)") ((I,EE_0(I),dFK1(I)),I=1,NE)  
!---
!Integrate out Kn: dFKn(1:Ndn), n=0,1,2         
!---
  X=0.d0
  F=0.d0
  DX=EE_0(2)-EE_0(1)    
  X(1:NE)=EE_0(1:NE)
  F(1:NE)=dFK0(1:NE)
  CALL  INX(DX,X,F,Z,NE,MSTART,Q,R,S)
  K0=Z(NE)
  !write(*,"(I5,2E16.6)") ((I,X(I),F(I)),I=1,NE)  
  !write(*,*) 'K0: ', K0
return
    end function K0 
!!===========================================================================     
    
function K1(T_in)
use global_parameters
real*8, intent(in) :: T_in
INTEGER:: MSTART=1
integer  NE
REAL*8 dfede,fe,EXPON,DX,K1,T
REAL*8, allocatable:: dFK1(:),X(:),F(:),Z(:),Q(:),R(:),S(:)
T=T_in
NE=NE_0
allocate(dFK1(NE),X(NE),F(NE),Z(NE),Q(NE),R(NE),S(NE))
DO I=1,NE
EXPON=DMin1((EE_0(I)-EF)/(AKB*T),300.d0) !to avoid overflow for exponential at low temperature
fe=1.d0/(DEXP(EXPON)+1.d0) !Fermi-Dirac distribution f
dfede=-1.d0*fe*(1.d0-fe)/(AKB*T) !df(E)/dE
!dFK0(I)=-1.d0*dfede*tau_0(I)
dFK1(I)=-1.d0*((EE_0(I)-EF)**1)*dfede*tau_0(I)
!dFK2(I)=-1.d0*((EE_0(I)-EF)**2)*dfede*tau_0(I)
END DO
!write(*,"(I5,2E16.6)") ((I,EE_0(I),dFK1(I)),I=1,NE)  
!---
!Integrate out Kn: dFKn(1:Ndn), n=0,1,2         
!---
  X=0.d0
  F=0.d0
  DX=EE_0(2)-EE_0(1)    
  X(1:NE)=EE_0(1:NE)
  F(1:NE)=dFK1(1:NE)
  CALL  INX(DX,X,F,Z,NE,MSTART,Q,R,S)
  K1=Z(NE)
  !write(*,"(I5,2E16.6)") ((I,X(I),F(I)),I=1,NE)  
  !write(*,*) 'K1: ', K1
return
    end function K1 
    
!!===========================================================================
    
function K2(T_in)
use global_parameters
real*8, intent(in) :: T_in
INTEGER:: MSTART=1
integer  NE
REAL*8 dfede,fe,EXPON,DX,K2,T
REAL*8, allocatable:: dFK2(:),X(:),F(:),Z(:),Q(:),R(:),S(:)
T=T_in
NE=NE_0
allocate(dFK2(NE),X(NE),F(NE),Z(NE),Q(NE),R(NE),S(NE))
DO I=1,NE
EXPON=DMin1((EE_0(I)-EF)/(AKB*T),300.d0) !to avoid overflow for exponential at low temperature
fe=1.d0/(DEXP(EXPON)+1.d0) !Fermi-Dirac distribution f
dfede=-1.d0*fe*(1.d0-fe)/(AKB*T) !df(E)/dE
!dFK0(I)=-1.d0*dfede*tau_0(I)
!dFK1(I)=-1.d0*((EE_0(I)-EF)**1)*dfede*tau_0(I)
dFK2(I)=-1.d0*((EE_0(I)-EF)**2)*dfede*tau_0(I)
END DO
!write(*,"(I5,2E16.6)") ((I,EE_0(I),dFK2(I)),I=1,NE)  
!---
!Integrate out Kn: dFKn(1:Ndn), n=0,1,2         
!---
  X=0.d0
  F=0.d0
  DX=EE_0(2)-EE_0(1)    
  X(1:NE)=EE_0(1:NE)
  F(1:NE)=dFK2(1:NE)
  CALL  INX(DX,X,F,Z,NE,MSTART,Q,R,S)
  K2=Z(NE)
  !write(*,"(I5,2E16.6)") ((I,X(I),F(I)),I=1,NE)  
  !write(*,*) 'K2: ', K2
return
end function K2       
    
    
