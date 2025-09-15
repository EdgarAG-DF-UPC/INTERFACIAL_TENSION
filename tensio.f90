SUBROUTINE CALCUL_TENSIO_SIMETRIA_ESFERICA(INPUT_FILE)
    ! This subroutine computes the interfacial tension given the kinetic and virial (configurational) terms of the stress tensor.
    ! Spherical symmetry.
    ! ARGUMENTS:
    ! ==> 'INPUT_FILE': name of the file containing r, p_K(r) and p_V(r) -- 1st, 5th and 6th column respectively.
    USE CONSTANTS
    IMPLICIT NONE
    CHARACTER*12 INPUT_FILE !="r_rho_pv.dat"
    INTEGER i, j, k
    INTEGER, parameter :: k0=50, k1=300
    REAL*8, dimension(k0:k1) :: r, pk, pv, pn
    REAL*8 a, b, c, INTEGRAL, Rs, GAMMA1, GAMMA2, GAMMA3!, GAMMA4, G4D, G4N, POTENCIA

    open(14, FILE=INPUT_FILE)
        do k = 1, k0-1
            read(14,*)
        enddo
        do k = k0, k1
            read(14,*) r(k), a, b, c, pk(k), pv(k)
            pn(k) = pk(k) + pv(k)
        enddo
        print*, INTEGRAL(r,pn,k0,k1,3d0), "eV"
        print*, INTEGRAL(r,pn,k0,k1,2d0), "eV Å^{-1}"
        print*, INTEGRAL(r,pn,k0,k1,0d0), "eV Å^{-3}"
        Rs = (INTEGRAL(r,pn,k0,k1,3d0) / INTEGRAL(r,pn,k0,k1,0d0))**(1d0/3d0) ! Å
        print*, Rs, "Å <--- Radius of the surface of tension"
        print*, ""
        print*, ""
        GAMMA1 = -(1d0/(2d0*Rs**2d0)) * INTEGRAL(r,pn,k0,k1,3d0) ! eV/Å² 
        print*, GAMMA1, "eV/Å²", &
        GAMMA1 * (1d20 * e), "N/m <--- Baker-Buff"

        GAMMA2 = -(Rs/(2d0)) * INTEGRAL(r,pn,k0,k1,0d0) ! eV/Å² 
        print*, GAMMA2, "eV/Å²", &
        GAMMA2 * (1d20 * e), "N/m <--- Mechanical stability" !equiv to YL

        GAMMA3 = -0.5d0*INTEGRAL(r,pn,k0,k1,0d0)*Rs

        print*, GAMMA3, "eV/Å²", &
        GAMMA3 * (1d20 * e), "N/m", &
        " <--- Young-Laplace"
        print*, (-(1d0/8d0)*INTEGRAL(r,pn,k0,k1,0d0)**2d0 * INTEGRAL(r,pn,k0,k1,3d0))**(1d0/3d0), "eV/Å²", &
        (-(1d0/8d0)*INTEGRAL(r,pn,k0,k1,0d0)**2d0 * INTEGRAL(r,pn,k0,k1,3d0))**(1d0/3d0) * (1d20 * e), "N/m", &
        " <--- Thompson et al"
        
        print*, -(INTEGRAL(r,pn,k0,k1,3d0) /(2d0*Rs**2d0)), "eV/Å²", &
        -(INTEGRAL(r,pn,k0,k1,3d0) /(2d0*Rs**2d0)) * (1d20 * e), "N/m"
        print*, (-INTEGRAL(r,pn,k0,k1,0d0)**2d0 * INTEGRAL(r,pn,k0,k1,3d0))**(1d0/3d0) / 2d0, "eV/Å²", &
        (-INTEGRAL(r,pn,k0,k1,0d0)**2d0 * INTEGRAL(r,pn,k0,k1,3d0))**(1d0/3d0) / 2d0 * (1d20 * e), "N/m"

    close(14)
    RETURN
END SUBROUTINE CALCUL_TENSIO_SIMETRIA_ESFERICA

SUBROUTINE CALCUL_TENSIO_PLANAR()
    ! This subroutine ...
    ! Planar symmetry.
    USE CONSTANTS
    IMPLICIT NONE
    INTEGER i, j, k
    INTEGER, parameter :: k0=50, k1=300
    REAL*8, dimension(k0:k1) :: r, pk, pv, pn
    REAL*8 a, b, c, INTEGRAL, Rs, GAMMA1, GAMMA2, GAMMA3, GAMMA4, G4D, G4N, POTENCIA
    
    open(14, FILE="r_rho_pv.dat")
        do k = 1, k0-1
            read(14,*)
        enddo
        do k = k0, k1
            read(14,*) r(k), a, b, c, pk(k), pv(k)
            pn(k) = pk(k) + pv(k)
            ! write(43,*) r(k), pk(k), pv(k), pn(k)
        enddo

        print*, "INTERFACIAL TENSION..."

    close(14)


END SUBROUTINE CALCUL_TENSIO_PLANAR

!========================================================================================================================
!Other functions:
REAL*8 FUNCTION POTENCIA(A, e1, e2)
    REAL*8, intent(in) :: A, e1, e2
    POTENCIA = dabs(A)**(e1/e2)
    ! IF (ISNAN(POTENCIA)) THEN
    !     POTENCIA = (-A) ** (e1/e2) * (-1d0)**e1 * (-1d0)
    ! ENDIF
    RETURN
END FUNCTION POTENCIA
REAL*8 FUNCTION INTEGRAL(r,pn,nr0,nr1,a)
    IMPLICIT NONE
    INTEGER, intent(in) :: nr0, nr1
    REAL*8, dimension(nr0:nr1), intent(in) :: r, pn
    REAL*8, intent(in) :: a
    INTEGER k

    INTEGRAL = 0d0
    do k = nr0+1, nr1-1
        integral = integral + (r(k)**a*(r(k+1)-r(k))/(r(k)-r(k-1))*(pn(k)-pn(k-1)) + &
	    r(k+1)**a*(pn(k+1)-pn(k))) ! integració per trapecis --- units: [integral] = [pN] [r]^a (typically eV/Å³*Å³ = eV when a=3)
    enddo

    INTEGRAL = INTEGRAL / 2d0

    RETURN
END FUNCTION INTEGRAL

