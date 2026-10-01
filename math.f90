
    !-------------------------------------------------------
      SUBROUTINE INX(DX,X,F,Z,N,MSTART,Q,R,S)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),F(N),Z(N),Q(N),R(N),S(N)
      CALL    SPLCO(1,N,X,F,Q,R,S)
      B=DX**2/2.D0
      C=DX**3/3.D0
      D=DX**4/4.D0
      NM1=N-1
      IF(MSTART) 20,20,10
10    Z(1)=0.D0
      DO 15  I=1,NM1
15    Z(I+1)=Z(I)+DX*F(I)+B*Q(I)+C*R(I)+D*S(I)
      RETURN
20    Z(N)=0.D0
      DO 25  I=1,NM1
      J=N-I
25    Z(J)=Z(J+1)+DX*F(J)+B*Q(J)+C*R(J)+D*S(J)
      RETURN
      END
!-------------------------------------------------------
      SUBROUTINE SPLCO(N1,N2,X,Y,B,C,D)
!     COMPUTATION OF THE COEFFICIENTS OF A CUBIC SPLINE
!     INTERPOLATING BETWEEN GIVEN DATA POINTS.
!     INPUT.
!     N1,N2  NUMBER OF FIRST AND LAST DATA POINT(0.LT.N1.LT.N2)
!     X(I),Y(I),I=N1,N1+1,...,N2 ARRAYS WITH X(1) AND Y(1) AS ABSCISSA
!          AND ORDINATE OF I-TH DATA POINT. THE COMPONENTS OF THE
!          ARRAY X MUST BE EITHER STRICTLY MONOTONIC INCREQSING
!          OR DECREASING.
!     OUTPUT.
!     B(I),C(I),D(I),I=N1,N1+1,..,N2 ARRAYS COLLECTING THE COEFFICIENTS
!          THE CUBIC SPLINE F(XX). IF XX LIES BETWEEN X(I) AND X(I+1),
!          THEN F(XX)=((D(I)*H+C(I))*H+B(I))*H+Y(I),
!     WHERE H=XX-X(I).
!     FURTHERMORE,C(N2)=0 WHILE B(N2) AND D(N2) ARE LEFT
!          UNDEFINED.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(1),Y(1),B(1),C(1),D(1)
      M1=N1+1
      M2=N2-1
      S=0.D0
      M3=M1+M2
      DO 1 K=N1,M2
      D(K)=X(K+1)-X(K)
      R=(Y(K+1)-Y(K))/D(K)
      C(K)=R-S
1     S=R
      C(N1)=0.D0
      C(N2)=0.D0
      S=0.D0
      R=0.D0
      DO 2 K=M1,M2
      C(K)=C(K)+R*C(K-1)
      B(K)=(X(K-1)-X(K+1))*2.D0-R*S
      S=D(K)
2     R=S/B(K)
      DO 3 K=M1,M2
      L=M3-K
3     C(L)=(D(L)*C(L+1)-C(L))/B(L)
      DO 4 K=N1,M2
      B(K)=(Y(K+1)-Y(K))/D(K)-(C(K)+C(K)+C(K+1))*D(K)
      D(K)=(C(K+1)-C(K))/D(K)
4     C(K)=3.D0*C(K)
      RETURN
    END
!-------------------------------------------------------
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!        Interpolation using Cubic spline  !!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=100000)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+    &
     1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*  &
     u(i-1))/p
       end do
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do
      return
    END
!==================================================================================================     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) write(*,*) 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      return
      END

