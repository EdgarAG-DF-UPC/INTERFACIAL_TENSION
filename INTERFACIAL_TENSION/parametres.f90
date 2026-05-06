MODULE PARAMETRES_LJ
    IMPLICIT NONE
    ! Lennard - Jones
    REAL*8, parameter :: eps_LiLi = 0.07651545d0
    REAL*8, parameter :: sig_LiLi = 2.7242d0

    REAL*8, parameter :: eps_PbPb = 0.06098297d0
    REAL*8, parameter :: sig_PbPb = 3.689d0

    REAL*8, parameter :: eps_HeHe = 0.00088069d0
    REAL*8, parameter :: sig_HeHe = 2.556d0

    REAL*8 , parameter ::eps_LiHe = 0.00014100d0
    REAL*8 , parameter ::sig_LiHe = 5.3565d0

    REAL*8, parameter :: eps_LiPb = 0.06618842d0
    REAL*8, parameter :: sig_LiPb = 3.244d0

    REAL*8, parameter :: eps_PbHe = 0.01057000d0
    REAL*8, parameter :: sig_PbHe = 3.0667d0

    ! EAM Li-Pb (Lennard-Jones tipus 8-4)
    REAL*8, parameter :: eps_84 = 0.0880d0
    REAL*8, parameter :: sig_84 = 2.553d0
END MODULE PARAMETRES_LJ

MODULE PARAMETRES_SLADEK
    IMPLICIT NONE
    ! Sladek
        !	dU[TQ]Z in Fig.1
        !	HF(QZ) + CORR(~r^{-3}, TQ)
    REAL*8, parameter :: aa = 4.6874d0
    REAL*8, parameter :: bb = 5.4114d0
    REAL*8, parameter :: cc = -20.7496d0
    REAL*8, parameter :: dd = 56.2663d0
    REAL*8, parameter :: ee = -95.6327d0
    REAL*8, parameter :: ff = 112.157d0
    REAL*8, parameter :: gg = -77.8808d0
    REAL*8, parameter :: hh = 39.4116d0
    REAL*8, parameter :: ii = 28.2574d0
    REAL*8, parameter :: jj = -3.7882d0
END MODULE PARAMETRES_SLADEK

MODULE PARAMETRES_EAM
    IMPLICIT NONE
    REAL*8, dimension(3), parameter :: p1=(/ 3.0511d0, 0d0, 5.1531d0 /)
    REAL*8, parameter :: p2=1.2200d0
END MODULE PARAMETRES_EAM

MODULE PARAMETRES_ZHOU
    IMPLICIT NONE
    CONTAINS
    SUBROUTINE SET_PARAMETRES_ZHOU(elem, versio, reO, feO, rhoeO, rhosO, alphaO, betaO, &
                            AO, BO, kappaO, lambdaO, FnO, FO, etaO, FFeO)
        IMPLICIT NONE
        CHARACTER*2, intent(in) :: elem, versio
        CHARACTER*2, dimension(1:16) :: elements
        REAL*8, dimension(1:16) :: re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, eta, FFe
        REAL*8, dimension(1:16,0:3) :: Fn, F
        REAL*8, intent(out) :: reO, feO, rhoeO, rhosO, alphaO, betaO, AO, BO, kappaO, lambdaO,&
                            FnO(0:3), FO(0:3), etaO, FFeO
        CHARACTER*6 cnt    
        INTEGER ielement
        !   PARENT DIRECTORY:
        CHARACTER*100 parent
        COMMON /OTHERS/ parent

        OPEN(14, FILE=trim(parent)//"/TAULES/PARAMS_ZHOU"//versio//".tab")
        READ(14,*) cnt, elements
        ielement = findloc(array=elements, value=elem, dim=1)
        ! ielement = 0
        ! do i = 1, size(elements)
        !     if (trim(elements(i)) == trim(elem)) then
        !         ielement = i
        !         exit
        !     end if
        ! end do

        if (ielement == 0) then
            print *, "ERROR: element not found: '", elem, "'"
            stop
        end if

        READ(14,*) cnt, re
        reO = re(ielement)
        READ(14,*) cnt, fe
        feO = fe(ielement)
        READ(14,*) cnt, rhoe
        rhoeO = rhoe(ielement)
        READ(14,*) cnt, rhos
        rhosO = rhos(ielement)
        READ(14,*) cnt, alpha
        alphaO = alpha(ielement)
        READ(14,*) cnt, beta
        betaO = beta(ielement)
        READ(14,*) cnt, A
        AO = A(ielement)
        READ(14,*) cnt, B
        BO = B(ielement)
        READ(14,*) cnt, kappa
        kappaO = kappa(ielement)
        READ(14,*) cnt, lambda
        lambdaO = lambda(ielement)
        READ(14,*) cnt, Fn(:,0)
        READ(14,*) cnt, Fn(:,1)
        READ(14,*) cnt, Fn(:,2)
        READ(14,*) cnt, Fn(:,3)
        FnO(:) = Fn(ielement, :)
        READ(14,*) cnt, F(:,0)
        READ(14,*) cnt, F(:,1)
        READ(14,*) cnt, F(:,2)
        READ(14,*) cnt, F(:,3)
        FO(:) = F(ielement, :)
        READ(14,*) cnt, eta
        etaO = eta(ielement)
        READ(14,*) cnt, FFe
        FFeO = FFe(ielement)
        CLOSE(14)

        return
    END SUBROUTINE SET_PARAMETRES_ZHOU
END MODULE PARAMETRES_ZHOU