SUBROUTINE PBC(r)
    ! Periodic boundary conditions: minimum distance convention
    IMPLICIT NONE
    INTEGER k
    REAL*8 r(1:3), boxl(1:3), boxh(1:3), box
    COMMON/CAPSA/boxl,boxh
    
    do k = 1,3
        box = boxh(k)-boxl(k)
        if (r(k) .lt. -0.5d0*box) r(k) = r(k) + box
        if (r(k) .ge. 0.5d0*box) r(k) = r(k) - box
    enddo
    
    return
END SUBROUTINE PBC

SUBROUTINE PBC_1D(k,r)
    ! Periodic boundary conditions: minimum distance convention
    IMPLICIT NONE
    INTEGER, intent(in) :: k
    REAL*8, intent(inout) :: r(1:3)
    REAL*8 boxl(1:3), boxh(1:3), box
    COMMON/CAPSA/boxl,boxh
    
    box = boxh(k)-boxl(k)
    if (r(k) .lt. -0.5d0*box) r(k) = r(k) + box
    if (r(k) .ge. 0.5d0*box) r(k) = r(k) - box
    
    return
END SUBROUTINE PBC_1D
    
!--------------------------------------------------------------------------------------------------------------------------------------
    
REAL*8 FUNCTION distancia(ri,rj,pbc)
    ! Aquesta funció calcula la distància més propera entre dues posicions ri, rj (PBC).
    IMPLICIT NONE
    REAL*8 ri(3), rj(3), Dr(3), boxl(1:3), boxh(1:3), box
    INTEGER k
    LOGICAL pbc
    COMMON/CAPSA/boxl,boxh
    
    !Dr = ri - rj
    !call pbc(Dr)
    
    if (pbc) then
        do k = 1, 3
            box = boxh(k)-boxl(k)
            ! Dr(k) = min(abs(rj(k)-ri(k)), box-abs(rj(k)-ri(k)))
            Dr(k) = ri(k)-rj(k)
            if (Dr(k) .ge. box*0.5d0) Dr(k) = Dr(k) - box
            if (Dr(k) .le. -box*0.5d0) Dr(k) = Dr(k) + box
        enddo

        distancia = dsqrt(dot_product(Dr,Dr))
        box = ((boxh(1)-boxl(1))*(boxh(2)-boxl(2))*(boxh(3)-boxl(3)))**(1d0/3d0)
        
        if (distancia .gt. dsqrt(3d0)*box*box/2d0) print*, "OJOOOO in distancia" , &
        distancia, dsqrt(dot_product(boxh-boxl,boxh-boxl)), box
    else
        distancia = dsqrt(dot_product(ri(:)-rj(:), ri(:)-rj(:)))
    endif
    return
END FUNCTION distancia

REAL*8 FUNCTION distancia_1D(k, ri, rj, pbc)
    ! Aquesta funció calcula la distància més propera entre dues coordenades xi, xj (PBC).
    REAL*8 ri(3), rj(3), Dx, boxl(1:3), boxh(1:3), box
    INTEGER k
    LOGICAL pbc
    COMMON/CAPSA/boxl,boxh
    
    if (pbc) then
        box = boxh(k) - boxl(k)
        Dx = rj(k) - ri(k)
        if (Dx .ge. box*0.5d0) Dx = Dx - box
        if (Dx .le. -box*0.5d0) Dx = Dx + box

        distancia_1D = Dx ! dsqrt(Dx*Dx)
        
        if (distancia_1D .gt. dsqrt(3d0)*box*box/2d0) print*, "OJOOOO in distancia_1D" , &
        distancia_1D, dsqrt((boxh-boxl)*(boxh-boxl)), box
    else
        distancia_1D = rj(k) - ri(k) !dsqrt(Dx*Dx)
    endif
    return
END FUNCTION distancia_1D

SUBROUTINE REPLICA(k, ri, rj)
    IMPLICIT NONE
    INTEGER, intent(in) :: k
    REAL*8, intent(in), dimension(1:3) :: ri
    REAL*8, intent(inout), dimension(1:3) :: rj
    REAL*8 Dx, boxl(1:3), boxh(1:3), box
    LOGICAL pbc
    COMMON/CAPSA/boxl,boxh

    if (pbc) then
        box = boxh(k) - boxl(k)
        ! Dx = ri(k) - rj(k)
        Dx = rj(k) - ri(k)
        if (Dx .ge. box*0.5d0) rj(k) = rj(k) - box
        if (Dx .le. -box*0.5d0) rj(k) = rj(k) + box

        if (Dx .gt. dsqrt(3d0)*box*box/2d0) print*, "OJOOOO in REPLICA" , &
        Dx, dsqrt((boxh-boxl)*(boxh-boxl)), box
    endif

    return
END SUBROUTINE REPLICA

! REAL*8 FUNCTION SUMA_1D(k, ri, rj, pbc)
!     ! Aquesta funció calcula la distància més propera entre dues coordenades xi, xj (PBC).
!     REAL*8 ri(3), rj(3), boxl(1:3), boxh(1:3), box, Dx
!     INTEGER k
!     LOGICAL pbc
!     COMMON/CAPSA/boxl,boxh
    
!     if (pbc) then
!         box = boxh(k) - boxl(k)
!         SUMA_1D = ri(:) + rj(:)
!         ! Dx = ri(k) - rj(k)
!         Dx = rj(k) - ri(k)
!         if (Dx .ge. box*0.5d0) SUMA_1D = SUMA_1D - box
!         if (Dx .le. -box*0.5d0) SUMA_1D = SUMA_1D + box

!         if (Dx .gt. dsqrt(3d0)*box*box/2d0) print*, "OJOOOO in SUMA_1D" , &
!         Dx, dsqrt((boxh-boxl)*(boxh-boxl)), box
!     else
!         SUMA_1D = ri(k) + rj(k) !dsqrt(Dx*Dx)
!     endif
!     return
! END FUNCTION SUMA_1D

! REAL*8 FUNCTION distancia(ri,rj)
!     IMPLICIT NONE
!     REAL*8, intent(in) :: ri(3), rj(3)
!     REAL*8 Dr(3)

!     distancia = dsqrt(dot_product(ri(:)-rj(:), ri(:)-rj(:)))

! END FUNCTION distancia


REAL*8 FUNCTION COM(x,N,L,m,SETPBC)
! Calcula el centre de masses a la coordenada x (unidimensional - 1D). IMPORTANT: 0<x(i)<L per a tota i=1,...,N !!!!
! x --> Coordenada de cadascuna de les N partícules.
! N --> Nombre de partícules.
! L --> Mida de la caixa.
! m --> Massa de cadascuna de les N partícules.
! SETPBC -->
    IMPLICIT NONE
    integer, intent(in) :: N
    real*8, intent(in) :: x(1:N), L, m(1:N)
    logical, intent(in) :: SETPBC
    integer i
    real*8 COM_PBC, mt

!    if (any(x.lt. 0d0)) then
!        print*, "ERROR! x must be greater than 0 to compute the center of masses..."
!        do i = i, N
!            if (x(i) .lt. 0d0) print*, x(i), " < ", 0
!        enddo
!        stop
!    endif

    COM = 0d0
    mt = 0d0
    if (SETPBC) then
        COM = COM_PBC(x,N,L,m) 
    else
        do i = 1, N
            COM = COM + m(i)*x(i)
            mt = mt + m(i)
        enddo
        COM = COM / mt
    endif
END FUNCTION COM

REAL*8 FUNCTION COM_PBC(x,N,L,m)
! Aquesta funció calcula la component x del centre de masses donats N àtoms i una capsa de costat L.
! Té en compte les PBC.
    IMPLICIT NONE
    integer, intent(in) :: N
    real*8, intent(in) :: x(1:N), m(1:N), L
    real*8 xi, dseta, theta, theta_i, sumM
    integer i
    real*8, parameter :: pi=4d0*datan(1d0)

    xi = 0d0
    dseta = 0d0
    sumM = 0d0
    do i = 1, N
        if ((x(i).lt.0d0).or.(x(i).ge.L)) print*, x(i), "OJO, la coordenada ha d'estar a [0,",L,")"
        theta_i = x(i)/L*2d0*pi
        xi = xi + m(i)*dcos(theta_i)
        dseta = dseta + m(i)*dsin(theta_i)
        sumM = sumM + m(i)
    enddo
    ! xi = xi / dble(N)
    ! dseta = dseta / dble(N)
    ! sumM = sumM / dble(N)
    xi = xi / sumM
    dseta = dseta / sumM

    theta = datan2(-dseta,-xi) + pi
    
    COM_PBC = L*theta/(2d0*pi)
    
    return
END FUNCTION COM_PBC
