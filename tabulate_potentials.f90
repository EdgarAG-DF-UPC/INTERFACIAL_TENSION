PROGRAM TAB_POT
    USE omp_lib
    USE CONSTANTS
    USE PARAMETRES_LJ
    USE PARAMETRES_EAM
    USE PARAMETRES_SLADEK
    USE PARAMETRES_ZHOU

    IMPLICIT NONE

    ! INTERFACE
    ! SUBROUTINE PARAMETRES_ZHOU(elem, versio, reO, feO, rhoeO, rhosO, alphaO, betaO, &
    !                             AO, BO, kappaO, lambdaO, FnO, FO, etaO, FFeO)
    !     CHARACTER*2, INTENT(IN) :: elem, versio
    !     REAL*8, INTENT(OUT) :: reO, feO, rhoeO, rhosO, alphaO, betaO
    !     REAL*8, INTENT(OUT) :: AO, BO, kappaO, lambdaO
    !     REAL*8, INTENT(OUT) :: FnO(0:3), FO(0:3)
    !     REAL*8, INTENT(OUT) :: etaO, FFeO
    ! END SUBROUTINE
    ! END INTERFACE


    REAL*8 PSIe(0:6), Ae(1:7), Be(1:7), Ce(1:7), me ! EAM Li-Li
    REAL*8 aaa(1:3,0:8), rp(1:4) ! EAM Pb-Pb
    CHARACTER*3 aim ! ---------     "    "    -------- 
    REAL*8 a1(1:3), a2(1:3), a3(1:3) ! Paràmetres dels models TTS
    REAL*8 alpha(1:3), lambda(1:3), gamma(1:3) ! Paràmetres dels models TTS
    REAL*8 b(1:3), A(1:3), c(1:3,6:16) ! Paràmetres dels models TTS
    REAL*8 Re(1:3), De(1:3), ZA, ZB ! Paràmetres dels models TTS
    REAL*8 U_short, U_long, dU_long, dU_short ! Contribucions al model TTS
    REAL*8 Us, Ul, delr, dUl, dUs, r0 ! Contribucions al model TTS    
    REAL*8 reZ, feZ, rhoeZ, rhosZ, alphaZ, betaZ, AZ, BZ, kappaZ, lambdaZ, FnZ(0:3), F0Z(0:3), etaZ, FFeZ

    REAL*8, dimension(1:1000) :: r, varphi, dvarphidr, rho, drhodr
    REAL*8, dimension(1:200) :: VARRHO, PHI, dPHIdr
    REAL*8 pot_ZHOU_Pb, der_pot_ZHOU_Pb
    REAL*8 Usladek, dUSladek
    INTEGER k
    CHARACTER*100 parent
    !   COMMON BLOCKS:
    COMMON /PbHePARAMS/ aa,bb,cc,dd,ee,ff,gg,hh,ii,jj
    COMMON /PbPbPARAMS/ aaa,rp,aim
    COMMON /LiLiPARAMS/ PSIe,Ae,Be,Ce,me
    COMMON /ZhouPARAMS/ reZ, feZ, rhoeZ, rhosZ, alphaZ, betaZ, AZ, BZ, kappaZ, lambdaZ, FnZ, F0Z, etaZ, FFeZ
    ! REAL*8 re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn(0:3), F(0:3), eta, FFe
    COMMON /CONVERS/ Eh, Rbohr
    COMMON /OTHERS/ parent


    print*, "SLADEK"
    open(42, file="SLADEK.dat")
    do k = 1, 1000
        r(k) = 1d-2 * float(k)
        write(42,*) r(k), USladek(r(k)), dUSladek(r(k))
    enddo
    close(42)

    print*, "ZHOU POTENTIALS"

    parent = "."
    call SET_PARAMETRES_ZHOU("Pb", "01", reZ, feZ, rhoeZ, rhosZ, alphaZ, betaZ, AZ, BZ, kappaZ, lambdaZ, FnZ, F0Z, etaZ, FFeZ)
    open(42, file="POTENTIAL_ZHOU01.dat")
    print*, "2001 version"
    do k = 1, 1000
        r(k) = 1d-2 * float(k)
        varphi(k) = pot_ZHOU_Pb(r(k))
        dvarphidr(k) =  der_pot_ZHOU_Pb(r(k))
        call electro_Zhou(r(k), rho(k), drhodr(k))
        write(42, *) r(k), pot_ZHOU_Pb(r(k)), der_pot_ZHOU_Pb(r(k)), rho(k), drhodr(k)
    enddo
    write(42, *)
    write(42, *)
    do k = 1, 200
        VARRHO(k) = 1d-2 * float(k)
        call embedding_Pb_Zhou(VARRHO(k), PHI(k), dPHIdr(k))
        write(42,*) VARRHO(k), PHI(k), dPHIdr(k)
    enddo
    close(42)

    call SET_PARAMETRES_ZHOU("Pb", "04", reZ, feZ, rhoeZ, rhosZ, alphaZ, betaZ, AZ, BZ, kappaZ, lambdaZ, FnZ, F0Z, etaZ, FFeZ)
    open(42, file="POTENTIAL_ZHOU04.dat")
    print*, "2004 version"
    do k = 1, 1000
        r(k) = 1d-2 * float(k)
        varphi(k) = pot_ZHOU_Pb(r(k))
        dvarphidr(k) =  der_pot_ZHOU_Pb(r(k))
        call electro_Zhou(r(k), rho(k), drhodr(k))
        ! call embedding_Pb_Zhou()
        write(42, *) r(k), pot_ZHOU_Pb(r(k)), der_pot_ZHOU_Pb(r(k)), rho(k), drhodr(k)
    enddo
    write(42, *)
    write(42, *)
    do k = 1, 200
        VARRHO(k) = 1d-2 * float(k)
        call embedding_Pb_Zhou(VARRHO(k), PHI(k), dPHIdr(k))
        write(42,*) VARRHO(k), PHI(k), dPHIdr(k)
    enddo
    close(42)
END PROGRAM TAB_POT