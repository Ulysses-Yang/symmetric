!2026/03/21
!symmetric bias version for I_el, JR_Peltier, JL_Peltier
!Fourier law of heat conduction: J_Fourier^{A->B}=- K * (T_B-T_A)
!    !
!use global_parameters
!real*8   K_TL2T0    ! thermal conductance between the left electrode and environment at T0
!real*8   K_TR2T0    ! thermal conductance between the left electrode and environment at T0
                     ! should be smaller than Kel(T0) to observe Peltier cooling effect)
!real*8   K_TLR2TC   ! take the value of the phonon's thermal conductance of the junction: "Kph" 
                     ! do not include "Kel" because its effect has been included in JL(R)_Peltier 
function Iel(Vb_in,TL_in,TR_in)   !Electric current I(Vb,TL,TR) at temperature TL and TR
!!output current Iel in unit of A 
use global_parameters
real*8, intent(in) :: Vb_in,TL_in,TR_in
INTEGER:: MSTART=1
integer  NE
REAL*8 dfede,fL,fR,EXPONL,EXPONR,DX,dmuL,dmuR,A,Iel
REAL*8 TL,TR
REAL*8, allocatable:: dIel(:),X(:),F(:),Z(:),Q(:),R(:),S(:)
NE=NE_0
TL=TL_in
TR=TR_in
!old: dmuL=EF
!ols: dmuR=EF+ee*Vb_in
!symetric bisa version:
dmuL=EF-0.5d0*ee*Vb_in
dmuR=EF+0.5d0*ee*Vb_in
!symetric bisa version:
A=2.d0*ee/h
allocate(dIel(NE),X(NE),F(NE),Z(NE),Q(NE),R(NE),S(NE))
DO I=1,NE
EXPONL=DMin1((EE_0(I)-dmuL)/(AKB*TL),300.d0) !to avoid overflow for exponential at low temperature
fL=1.d0/(DEXP(EXPONL)+1.d0) !Fermi-Dirac distribution of lead electrode
EXPONR=DMin1((EE_0(I)-dmuR)/(AKB*TR),300.d0) !to avoid overflow for exponential at low temperature
fR=1.d0/(DEXP(EXPONR)+1.d0) !Fermi-Dirac distribution of lead electrode
!
dIel(I)=A*(fR-fL)*tau_0(I)
!
END DO
!---
!Integrate out dIel(I)         
!---
  X=0.d0
  F=0.d0
  DX=EE_0(2)-EE_0(1)    
  X(1:NE)=EE_0(1:NE)
  F(1:NE)=dIel(1:NE)
  CALL  INX(DX,X,F,Z,NE,MSTART,Q,R,S)
  Iel=Z(NE)
return
    end function Iel
!!===========================================================================

function JR_Peltier(Vb_in,TL_in,TR_in)   !JR_Peltier is the rate of heat flow entering into the right electrode of the nanojunction 
!! JR_Peltier in unit of W/s 
!The sign of JR_Peltier is defined as positive at T=0 K (No Peltier cooling is possible at T=0 K). 
use global_parameters
real*8, intent(in) :: Vb_in,TL_in,TR_in
INTEGER:: MSTART=1
integer  NE
REAL*8 dfede,fL,fR,EXPONL,EXPONR,DX,dmuL,dmuR,A,JR_Peltier
REAL*8 TL,TR
REAL*8, allocatable:: dJR_Peltier(:),X(:),F(:),Z(:),Q(:),R(:),S(:)
NE=NE_0
TL=TL_in
TR=TR_in
!old: dmuL=EF
!ols: dmuR=EF+ee*Vb_in
!symetric bisa version:
dmuL=EF-0.5d0*ee*Vb_in
dmuR=EF+0.5d0*ee*Vb_in
!symetric bisa version:
A=2.d0/h
allocate(dJR_Peltier(NE),X(NE),F(NE),Z(NE),Q(NE),R(NE),S(NE))
DO I=1,NE
EXPONL=DMin1((EE_0(I)-dmuL)/(AKB*TL),300.d0) !to avoid overflow for exponential at low temperature
fL=1.d0/(DEXP(EXPONL)+1.d0) !Fermi-Dirac distribution of lead electrode
EXPONR=DMin1((EE_0(I)-dmuR)/(AKB*TR),300.d0) !to avoid overflow for exponential at low temperature
fR=1.d0/(DEXP(EXPONR)+1.d0) !Fermi-Dirac distribution of lead electrode
!
dJR_Peltier(I)=A*(fR-fL)*tau_0(I) * (dmuR-EE_0(I))
!
END DO
!---
!Integrate out dJR_Peltier(I)         
!---
  X=0.d0
  F=0.d0
  DX=EE_0(2)-EE_0(1)    
  X(1:NE)=EE_0(1:NE)
  F(1:NE)=dJR_Peltier(1:NE)
  CALL  INX(DX,X,F,Z,NE,MSTART,Q,R,S)
  JR_Peltier=Z(NE)
!
return
    end function JR_Peltier
!!===========================================================================

function JL_Peltier(Vb_in,TL_in,TR_in)   !JL_Peltier is the rate of heat flow into the left electrode of the nanojunction 
!! JL_Peltier in unit of W/s
!The sign of JL_Peltier is defined as positive at T=0 K (No Peltier cooling is possible at T=0 K). 
use      global_parameters
real*8, intent(in) :: Vb_in,TL_in,TR_in
INTEGER:: MSTART=1
integer  NE
REAL*8 dfede,fL,fR,EXPONL,EXPONR,DX,dmuL,dmuR,A,JL_Peltier
REAL*8 TL,TR
REAL*8, allocatable:: dJL_Peltier(:),X(:),F(:),Z(:),Q(:),R(:),S(:)
NE=NE_0
TL=TL_in
TR=TR_in
!old: muL=EF
!ols: dmuR=EF+ee*Vb_in
!symetric bisa version:
dmuL=EF-0.5d0*ee*Vb_in
dmuR=EF+0.5d0*ee*Vb_in
!symetric bisa version:
A=2.d0/h
allocate(dJL_Peltier(NE),X(NE),F(NE),Z(NE),Q(NE),R(NE),S(NE))
DO I=1,NE
EXPONL=DMin1((EE_0(I)-dmuL)/(AKB*TL),300.d0) !to avoid overflow for exponential at low temperature
fL=1.d0/(DEXP(EXPONL)+1.d0) !Fermi-Dirac distribution of lead electrode
EXPONR=DMin1((EE_0(I)-dmuR)/(AKB*TR),300.d0) !to avoid overflow for exponential at low temperature
fR=1.d0/(DEXP(EXPONR)+1.d0) !Fermi-Dirac distribution of lead electrode
!
dJL_Peltier(I)=A*(fR-fL)*tau_0(I) * (EE_0(I)-dmuL)
!
END DO
!---
!Integrate out dJL_Peltier(I)         
!---
  X=0.d0
  F=0.d0
  DX=EE_0(2)-EE_0(1)    
  X(1:NE)=EE_0(1:NE)
  F(1:NE)=dJL_Peltier(1:NE)
  CALL  INX(DX,X,F,Z,NE,MSTART,Q,R,S)
  JL_Peltier=Z(NE)
!
return
    end function JL_Peltier

!!===========================================================================
!! The rate of heat absorbed or released in the central channel region due to Thomson effect:
!! definition version 2:    
!! P_Thom = (+) 1.d0 * tau_Thom(T) * I * \Delta T = (+) 1.d0 * tau_Thom(T) * I * (TR- TL)   
!! P_Thom in unit of W/s
! 
function P_Thom(Vb_in, TL_in, TC_in, TR_in) 
use global_parameters
real*8, intent(in) :: Vb_in,TL_in,TC_in,TR_in
REAL*8,external :: Iel,tau_Thom
REAL*8 Vb2,TL,TC,TR,P_Thom 
TL=TL_in
TC=TC_in
TR=TR_in
Vb2=Vb_in    
!P_Thom =(1.d0)*tau_Thom(TC)*Iel(Vb2,TL,TR)*(TR-TL) 
P_Thom =(1.d0)*tau_Thom(T0)*Iel(Vb2,TL,TR)*(TR-TL)
!
return
end function  P_Thom             

!!===============================================================================   
function J_L2Evn_Fourier(TLin)   !Thermal current J^{L->Env}_Fourier(TL,T0,K_L2T0)
!! J_L2Evn_Fourier in unit of W/s 
use global_parameters  !K_TL2T0 ,and T0 from global_parameters
real*8, intent(in) :: TLin
real*8 TL,J_L2Evn_Fourier
TL=TLin
!Fourier's law for thermal current from "Left" electrod to "environment"    
J_L2Evn_Fourier= -1.d0* K_TL2T0 * (T0-TL) 
!
return
    end function J_L2Evn_Fourier
!!===============================================================================   
function J_Evn2L_Fourier(TLin)   !Thermal current J^{L->Env}_Fourier(TL,T0,K_L2T0)
!! J_Evn2L_Fourier in unit of W/s 
use global_parameters  !K_TL2T0 and T0 from global_parameters
real*8, intent(in) :: TLin
real*8 TL,J_Evn2L_Fourier
TL=TLin
!Fourier's law for thermal current from "Left" electrod to "environment"    
J_Evn2L_Fourier= -1.d0* K_TL2T0 * (TL-T0) 
!
return
end function J_Evn2L_Fourier
    
!!===============================================================================    
function J_R2Evn_Fourier(TRin)   !Thermal current J^{R->Env}_Fourier(TR,T0,K_R2T0)
!! J_R2Evn_Fourier in unit of W/s 
use global_parameters  !K_TR2T0 from global_parameters
real*8, intent(in) :: TRin
real*8 TR,J_R2Evn_Fourier
TR=TRin
!Fourier's law for thermal current from "Right" electrod to "environment"    
J_R2Evn_Fourier= -1.d0* K_TR2T0 * (T0-TR) 
!
return
    end function J_R2Evn_Fourier 
!!===============================================================================    
function J_Evn2R_Fourier(TRin)   !Thermal current J^{R->Env}_Fourier(TR,T0,K_R2T0)
!! J_Evn2R_Fourier in unit of W/s 
use global_parameters  !K_TR2T0 from global_parameters
real*8, intent(in) :: TRin
real*8 TR,J_Evn2R_Fourier
TR=TRin
!Fourier's law for thermal current from "Right" electrod to "environment"    
J_Evn2R_Fourier= -1.d0* K_TR2T0 * (TR-T0) 
!
return
end function J_Evn2R_Fourier     
    
!!1===============================================================================    
function J_R2C_Fourier(TRin,TCin)   !Thermal current J^{R->C}_Fourier(TR,TC,K_TLR2TC)
!! J_R2C_Fourier in unit of W/s 
use global_parameters   !K_TLR2TC from global_parameters
real*8, intent(in) :: TRin,TCin
real*8 TR,TC,J_R2C_Fourier
REAL*8, external:: Kel
TR=TRin
TC=TCin 
!Fourier's law for thermal current from the "Right" electrode to "Center Channel"
K_TLR2TC=Kph+Kel(TC)
J_R2C_Fourier= -1.d0* K_TLR2TC * (TC-TR)
!
return
end function J_R2C_Fourier 

function J_R2C_Fourier_noKel(TRin,TCin)   !Thermal current J^{R->C}_Fourier(TR,TC,K_TLR2TC)
!! J_R2C_Fourier in unit of W/s 
use global_parameters   !K_TLR2TC from global_parameters
real*8, intent(in) :: TRin,TCin
real*8 TR,TC,J_R2C_Fourier_noKel
REAL*8, external:: Kel
TR=TRin
TC=TCin 
!Fourier's law for thermal current from the "Right" electrode to "Center Channel"
K_TLR2TC=Kph
J_R2C_Fourier_noKel= -1.d0* K_TLR2TC * (TC-TR)
!
return
end function J_R2C_Fourier_noKel
    
!!2===============================================================================    
function J_C2R_Fourier(TCin,TRin)   !Thermal current J^{C->R}_Fourier(TR,TC,K_TLR2TC)
!! J_C2R_Fourier in unit of W/s 
use global_parameters    !K_TLR2TC from global_parameters
real*8, intent(in) :: TRin,TCin
real*8 TR,TC,J_C2R_Fourier
REAL*8, external:: Kel
TR=TRin
TC=TCin
!Fourier's law for thermal current from "Center Channel" to the "Right" electrode to
K_TLR2TC=Kph !Delete Kel
J_C2R_Fourier= -1.d0* K_TLR2TC * (TR-TC) 
!
return
end function J_C2R_Fourier       
!!3===============================================================================    
function J_L2C_Fourier(TLin,TCin)   !Thermal current J^{L->C}_Fourier(TL,TC,K_TLR2TC)
!! J_L2C_Fourier in unit of W/s 
use global_parameters     !K_TLR2TC from global_parameters
real*8, intent(in) :: TLin,TCin
real*8 TL,TC,J_L2C_Fourier
REAL*8, external:: Kel
TL=TLin
TC=TCin 
!Fourier's law for thermal current from the "Left" electrode to "Center Channel"
K_TLR2TC=Kph+Kel(TC)
J_L2C_Fourier= -1.d0* K_TLR2TC * (TC-TL)
!
return
end function J_L2C_Fourier 

function J_L2C_Fourier_noKel(TLin,TCin)   !Thermal current J^{L->C}_Fourier(TL,TC,K_TLR2TC)
!! J_L2C_Fourier in unit of W/s 
use global_parameters     !K_TLR2TC from global_parameters
real*8, intent(in) :: TLin,TCin
real*8 TL,TC,J_L2C_Fourier_noKel
REAL*8, external:: Kel
TL=TLin
TC=TCin 
!Fourier's law for thermal current from the "Left" electrode to "Center Channel"
K_TLR2TC=Kph
 J_L2C_Fourier_noKel= -1.d0* K_TLR2TC * (TC-TL)
!
return
end function  J_L2C_Fourier_noKel

!!4===============================================================================    
function J_C2L_Fourier(TCin,TLin)   !Thermal current J^{C->L}_Fourier(TL,TC,K_TLR2TC)
!! J_C2L_Fourier in unit of W/s 
use global_parameters     !K_TLR2TC from global_parameters
real*8, intent(in) :: TLin,TCin
real*8 TL,TC,J_C2L_Fourier
REAL*8, external:: Kel
TL=TLin
TC=TCin
!Fourier's law for thermal current from "Center Channel" to the "Left" electrode to 
K_TLR2TC=Kph !Delete Kel
J_C2L_Fourier= -1.d0* K_TLR2TC * (TL-TC) 
!
return
    end function J_C2L_Fourier 
!!===============================================================================     