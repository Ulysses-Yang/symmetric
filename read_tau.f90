module global_parameters
  implicit none
  REAL*8, parameter:: ee=1.602176634D-19,AKB=1.38064853D-23,h=6.626069934D-34
  REAL*8, parameter:: hbar=1.0545718D-34,EF=0.d0,PI=3.1415926535897932384626433832795
  REAL*8, parameter:: m0=9.10938356D-31 
  REAL*8, parameter:: alpha_in=0.83, alpha_out=0.33, delta_E=1.d-4 !interval of energy 1.d-3 eV recommended) 
  REAL*8, parameter:: T0=300.d0 ! Envoronment temperature for in K
  integer,parameter :: N_Vb=51 , N_alpha=51  ! number of bias voltage and alpha for calculation
  REAL*8,save::  SS1TL(N_alpha,N_Vb), SS1TC(N_alpha,N_Vb), SS1TR(N_alpha,N_Vb)
  real*8,save::  K_TL2T0    ! thermal conductance between the left electrode and environment at T0
  real*8,save::  K_TR2T0    ! thermal conductance between the left electrode and environment at T0
                      ! should be smaller than Kel(T0) to observe Peltier cooling effect)
  real*8,save::  K_TLR2TC   ! take the value of the phonon's thermal conductance of the junction: "Kph+Kel" 
  real*8,save::  Kph                         
  integer,save:: NE_0
  real*8,save::  TL_guess, TC_guess, TR_guess 
  integer I, J, K, N
  REAL*8, allocatable, save :: EE_0(:), tau_0(:)
  real*8,save :: Vb, Vb_min, Vb_max
end module global_parameters
    
module read_tau_module
    
contains    
subroutine read_tau
use global_parameters
!!============================== 
! 
character (len=40):: filename="tau.inp"
LOGICAL :: alive_tau_inp !check files if exit
integer :: Ntauin_0, NEmax, NEmin,NE
logical :: alive
REAL*8, allocatable :: x(:),y(:),y2(:)
Real*8    temp1,temp2,y_out,x_in,yp1,ypn
!
!!==========================================================================
inquire(file=filename, exist=alive)
if(alive) then  !check if tau.inp [tau(E)] exist
  open(10,file='tau.inp', status='old', action='read')
  open(20,file='tau_interpolation.csv')
else
  write(*,*) 'tau.inp does not exist. Please check input of tau(E).'
end if
!!!*    tau.inp [tau(E)]exit:
     I=0
  DO WHILE(.TRUE.)     
	 READ(10,*,END=200)temp1,temp2
     I=I+1
  ENDDO
200 CONTINUE
    Ntauin_0=I
!    write(*,*) 'Ntauin_0: ', Ntauin_0
allocate(x(Ntauin_0),y(Ntauin_0),y2(Ntauin_0))   
     rewind(10) 
  DO I=1,Ntauin_0
	 READ(10,*,END=201)x(I),y(I)
  ENDDO
201 CONTINUE 
!!*===============     
    !!*extrapolation with an interval of  delta_E:
   NE_0=NINT((x(Ntauin_0)-x(1))/(delta_E))+1 
   allocate(EE_0(NE_0), tau_0(NE_0))  
   do i=1,NE_0          
     EE_0(i)=x(1)+(i-1)*delta_E
   enddo
!write(*,"(I5,F16.6)") ((I,EE_0(I)),I=1,NE_0)  
!Write(*,*) NE_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interpolate tau.inp (obtained from Nanodcal, spin unploarized case)
! interpolation using cubic spine: "splint" and "spline" from Nuemerical Recipe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           yp1=(y(2)-y(1))/(x(2)-x(1))
           ypn=(y(Ntauin_0)-y(Ntauin_0-1))/(x(Ntauin_0)-x(Ntauin_0-1))
!write(*,*)yp1,ypn         
           call spline(x,y,Ntauin_0,yp1,ypn,y2)
           do i=1, NE_0
               x_in=EE_0(i)
               if(x_in<x(1)) then
                   tau_0(i)=y(1)
               else if(x_in>x(Ntauin_0)) then
                   tau_0(i)=y(Ntauin_0)
               else                       
                   call splint(x,y,y2,Ntauin_0,x_in,y_out)
                   tau_0(i)=y_out
               endif
           enddo          
!!*
777 format(1x, *(g0, ", "))            
                   write(20,777)'E','tau','E','tau_intp'
                   write(20,777)'eV','','eV','no dim.'
                   write(20,777)'Nanodcal','','','cubic spline'
!!*find NEmax=Max(NE_0,Ntauin_0);NEmin=MIN(NE_0,Ntauin_0)  
    NEmax=MAX(NE_0,Ntauin_0)
    NEmin=MIN(NE_0,Ntauin_0)
    do I=1,NEmax
!old           do i=1,NE_0
     if(i.LE.NEmin) then
         write(20,777)x(i),y(i),EE_0(i),tau_0(i) 
     else
         if(NEmin.EQ.Ntauin_0) write(20,777) '','',EE_0(i),tau_0(i)
         if(NEmin.EQ.NE_0) write(20,777) x(i),y(i),'',''         
     endif
    enddo 
    EE_0(:)=ee*EE_0(:)  ! convert energy E in eV to Joule in SI unit
end subroutine read_tau
end module  read_tau_module    
    