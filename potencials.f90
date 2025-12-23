REAL*8 FUNCTION U_LJ(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    U_LJ = 4d0 * (r**(-12d0) - r**(-6d0))
    return
END FUNCTION U_LJ
REAL*8 FUNCTION dU_LJ(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    dU_LJ = -24d0 * (2d0*r**(-13d0) - r**(-7d0))
    return
END FUNCTION dU_LJ


! TTS model:
REAL*8 FUNCTION UTTS_short(x,a1,a2,a3,alpha)
    IMPLICIT NONE
    REAL*8 x
    REAL*8 a1, a2, a3, alpha
    UTTS_short = (1d0/x) * (1d0 + a1*x + a2*x**2 + a3*x**3)*dexp(-alpha*x)
    return
END FUNCTION UTTS_short
REAL*8 FUNCTION DUTTS_short(x,a1,a2,a3,alpha)
    IMPLICIT NONE
    REAL*8 x
    REAL*8 a1, a2, a3, alpha, UTTS_short

    DUTTS_short = - UTTS_short(x,a1,a2,a3,alpha) * (1d0/x + alpha) + &
    (1d0/x) * ( a1 + 2.d0 * a2 * x + 3.d0* a3 * x**2) *dexp(-alpha*x)

    return
END FUNCTION DUTTS_short
REAL*8 FUNCTION UTTS_long(x,alpha,lambda,gamma,b,A,c)
    IMPLICIT NONE
    REAL*8 x
    REAL*8 a1, a2, a3
    REAL*8 alpha, lambda, gamma, b, A, c(6:16)
    REAL*8 U_vdw
    REAL*8 sumn, sumk
    INTEGER n, k
    REAL*8 factk
    sumn = 0d0
    do n = 3, 8
    sumk = 1d0
    factk = 1d0
    do k = 1, 2*n
        factk = factk * dble(k)
        sumk = sumk + (b*x)**k / factk
    enddo
    sumn = sumn + (1d0 - dexp(-b*x)*sumk)*c(2*n)/(x**(2*n))
    enddo
    U_vdw = A * x**gamma * dexp(-lambda*x) - sumn
    UTTS_long = (1d0 - dexp(-alpha*x)) * U_vdw
    return
END FUNCTION UTTS_long
REAL*8 FUNCTION DUTTS_long(x,alpha,lambda,gamma,b,A,c)
    IMPLICIT NONE
    REAL*8 x
    REAL*8 a1, a2, a3
    REAL*8 alpha, lambda, gamma, b, A, c(6:16)
    REAL*8 U_vdw,D_U_vdw_0,D_U_vdw_1,D_U_vdw
    REAL*8 sumn, sumk, sumn_1,sumn_2,sumn_3,sumk_0
    INTEGER n, k
    INTEGER*8 factk
    sumn = 0d0
    sumn_1 = 0d0
    sumn_2 = 0d0
    sumn_3 = 0d0
    do n = 3, 8
        sumk = 1d0
        sumk_0 = 0d0
        factk = 1
        do k = 1, 2*n
            factk = factk * k
            sumk = sumk + (b*x)**k / dble(factk)
            sumk_0 = sumk_0 + (b*x)**k / dble(factk) / x * dble(k)
        enddo
        sumn = sumn + (1d0 - dexp(-b*x)*sumk)*c(2*n)/(x**(2*n))
        sumn_1 = sumn_1 - (1d0 - dexp(-b*x)*sumk)*c(2*n)/(x**(2*n+1)) * 2 * n
        sumn_2 = sumn_2 +  b * dexp(-b*x) * sumk * c(2*n) / (x**(2*n))
        sumn_3 = sumn_3 - dexp(-b*x) * sumk_0 * c(2*n)/(x**(2*n))
    enddo
    U_vdw = A * x**gamma * dexp(-lambda*x) - sumn
    D_U_vdw_0 = A * dexp(-lambda*x) * x**gamma * ( gamma / x - lambda)
    D_U_vdw_1 = sumn_1 + sumn_2 + sumn_3
    D_U_vdw = D_U_vdw_0 - D_U_vdw_1
    DUTTS_long = (1d0 - dexp(-alpha*x)) * D_U_vdw + alpha * dexp(-alpha*x) * U_vdw
    return
END FUNCTION DUTTS_long

REAL*8 FUNCTION pot_EAM_Li(r)
    ! returns the pairwise potential of lithium-lithium (in eV)
    IMPLICIT NONE
    REAL*8 r, phi
    REAL*8 k(0:5), k1, k2, k3, k4
    INTEGER i
    k(0) = -0.161539351212d1
    k(1) = 0.329193195820d2
    k(2) = -0.245830404172d3
    k(3) = 0.840217873656d3
    k(4) = -0.136938125679d4
    k(5) = 0.905623694715d3
    k1 = 0.252868d0
    k2 = 0.15252d0
    k3 = 0.38d0
    k4 = 1.96d0
    phi = 0d0
    if (r .le. 7.5d0) then
        if (r .gt. 2.45d0) then
            do i=0,5
                phi = phi + k(i)/r**dble(i)
            enddo
        else
            phi = k1+k2*(2.45d0-r)+k3*(dexp(k4*(2.45d0-r))-1d0)
        endif
    endif
    pot_EAM_Li = phi
    return
END FUNCTION pot_EAM_Li
REAL*8 FUNCTION der_pot_EAM_Li(r)
    ! returns the pairwise potential of lithium-lithium (in eV)
    IMPLICIT NONE
    REAL*8 r, dphidr
    REAL*8 k(0:5), k1, k2, k3, k4
    INTEGER i
    k(0) = -0.161539351212d1
    k(1) = 0.329193195820d2
    k(2) = -0.245830404172d3
    k(3) = 0.840217873656d3
    k(4) = -0.136938125679d4
    k(5) = 0.905623694715d3
    k1 = 0.252868d0
    k2 = 0.15252d0
    k3 = 0.38d0
    k4 = 1.96d0
    dphidr = 0d0
    if (r .le. 7.5d0) then
        if (r .gt. 2.45d0) then
            do i=2,6
                dphidr = dphidr - (i-1)*k(i-1)/r**dble(i)
            enddo
        else
            dphidr = -k2 - k3*k4*dexp(k4*(2.45d0-r))
        endif
    endif
    der_pot_EAM_Li = dphidr
    return
END FUNCTION der_pot_EAM_Li

REAL*8 FUNCTION pot_AEAM_Li(r) !eV
    IMPLICIT none
    REAL*8 r, potene
    !Attractive term parameters: 
    REAL*8, parameter :: A1=0.290231d0, A2=1.383589d0
    REAL*8, parameter :: B1=1.692128d0, B2=4.237653d0
    REAL*8, parameter :: C1=1.128630d0, C2=-1.802154d0
    REAL*8, parameter :: D1=0.290181d0, D2=0.290311d0
    REAL*8, parameter :: E1=1.307873d0, E2=-0.000074d0, E3=1.021786d0, E4=0.000066d0
    REAL*8, parameter :: F1=0.740131d0, F2=0.753327d0
    REAL*8, parameter ::  r1=0.312355d0, r2=0.787057d0, r3=1.212659d0, r4=1.641811d0, r5=0.524198d0, r6=1.204905d0
    !Repulsive term parameters:
    REAL*8, parameter ::  A=1.209005d0, B=2.218267d0, C=2.115713d0, D=2.280588d0, E=8.313021d0, F=-0.75781d0, G=-7.187612d0
    !Parametrization coefficients of the NPA and LRO potentials:
    ! REAL*8, parameter :: w=1.00d0, rc=7.50d0, S = 0.00d0, sigma=1.00d0, rs = 0.00d0 !NPA
    REAL*8, parameter :: w=0.825d0, rc=7.50d0, S=-0.00046d0, sigma=1.21d0, rs=0.07d0
    
    potene = 0d0

    if (r .lt. rc) then
        if (r .ge. 2.0d0 + rs) then
            potene = w*A1*(r1/(r-rs))**B1 * dcos(E1*((r-rs)/r2 + C1)) + E2*dexp(D1+F1*(r-rs)/r3) + &
            w*A2*(r4/(r-rs))**B2 * dsin(E3*((r-rs)/r5 + C2)) + E4*dexp(D2+F2*(r-rs)/r6) + &
            S*dexp((rc-r)/sigma)
        else
            potene = A*(r/B)**C + D*(r/E)**F + G
        endif
    endif
    pot_AEAM_Li = potene
    return
END FUNCTION pot_AEAM_Li
REAL*8 FUNCTION der_pot_AEAM_Li(r) !eV
    IMPLICIT none
    REAL*8 r, dphi
    !Attractive term parameters: 
    REAL*8, parameter :: A1=0.290231d0, A2=1.383589d0
    REAL*8, parameter :: B1=1.692128d0, B2=4.237653d0
    REAL*8, parameter :: C1=1.128630d0, C2=-1.802154d0
    REAL*8, parameter :: D1=0.290181d0, D2=0.290311d0
    REAL*8, parameter :: E1=1.307873d0, E2=-0.000074d0, E3=1.021786d0, E4=0.000066d0
    REAL*8, parameter :: F1=0.740131d0, F2=0.753327d0
    REAL*8, parameter ::  r1=0.312355d0, r2=0.787057d0, r3=1.212659d0, r4=1.641811d0, r5=0.524198d0, r6=1.204905d0
    !Repulsive term parameters:
    REAL*8, parameter ::  A=1.209005d0, B=2.218267d0, C=2.115713d0, D=2.280588d0, E=8.313021d0, F=-0.75781d0, G=-7.187612d0
    !Parametrization coefficients of the NPA and LRO potentials:
    ! REAL*8, parameter :: w=1.00d0, rc=7.50d0, S = 0.00d0, sigma=1.00d0, rs = 0.00d0 !NPA
    REAL*8, parameter :: w=0.825d0, rc=7.50d0, S=-0.00046d0, sigma=1.21d0, rs=0.07d0
    
    dphi = 0d0

    if (r .lt. rc) then
        if (r .ge. 2.0d0 + rs) then
            dphi = w*A1*(r1/(r-rs))**B1 * (-B1/(r-rs)*dcos(E1*((r-rs)/r2 + C1)) - E1/r2*dsin(E1*((r-rs)/r2 + C1))) + &
            (E2*F1/r3)*dexp(D1+F1*(r-rs)/r3) + &
            w*A2*(r4/(r-rs))**B2 * (-B2/(r-rs)*dsin(E3*((r-rs)/r5 + C2)) + E3/r5*dcos(E3*((r-rs)/r5 + C2))) + &
            (E4*F2/r6)*dexp(D2+F2*(r-rs)/r6) + &
            (-S/sigma)*dexp((rc-r)/sigma)
        else
            dphi = (C*A/B)*(r/B)**(C-1d0) + (D*F/E)*(r/E)**(F-1d0)
        endif
    endif
    der_pot_AEAM_Li = dphi
    return
END FUNCTION der_pot_AEAM_Li

REAL*8 FUNCTION pot_EAM_LiPb(r)
    IMPLICIT NONE
    REAL*8 r
    pot_EAM_LiPb = 4d0 * (r**(-8d0) - r**(-4d0))
    return
END FUNCTION pot_EAM_LiPb
REAL*8 FUNCTION der_pot_EAM_LiPb(r)
    IMPLICIT NONE
    REAL*8 r
    der_pot_EAM_LiPb = -16d0 * (2d0*r**(-9d0) - r**(-5d0))
    return
END FUNCTION der_pot_EAM_LiPb

REAL*8 FUNCTION pot_AEAM_LiPb(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    REAL*8, parameter :: eps=0.115d0, r0=3.06d0,alpha=4.35d0 ! eV, Å, adim
    ! pot_AEAM_LiPb = eps * (dexp(-2d0*alpha*(r-r0)) - 2d0*dexp(-alpha*(r-r0)))     !!!! Notar error a l'apendix de Awad et al (alpha ha de ser adimensional en comptes de 1/Å!!!!!)
    pot_AEAM_LiPb = eps * (dexp(-2d0*alpha*(r/r0-1d0)) - 2d0*dexp(-alpha*(r/r0-1d0)))
    ! print*, pot_AEAM_LiPb, &
    ! eps * (dexp(-2d0*alpha*(r/r0-1)) - 2d0*dexp(-alpha*(r/r0-1)))
    return
END FUNCTION pot_AEAM_LiPb
REAL*8 FUNCTION der_pot_AEAM_LiPb(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    REAL*8, parameter :: eps=0.115d0, r0=3.06d0,alpha=4.35d0 ! eV, Å, adim
    ! der_pot_AEAM_LiPb = -2d0*alpha*eps * (dexp(-2d0*alpha*(r-r0)) - dexp(-alpha*(r-r0)))
    der_pot_AEAM_LiPb = -2d0*(alpha/r0)*eps * (dexp(-2d0*alpha*(r/r0-1d0)) - dexp(-alpha*(r/r0-1d0)))
    ! print*, der_pot_AEAM_LiPb, &
    ! -2d0*(alpha/r0)*eps * (dexp(-2d0*alpha*(r/r0-1)) - 2d0*dexp(-alpha*(r/r0-1)))
    return
END FUNCTION der_pot_AEAM_LiPb

REAL*8 FUNCTION pot_EAM_Pb(r)
    IMPLICIT NONE
    INTEGER i, n, m, k
    PARAMETER(n=3,k=8)
    REAL*8 r, phi, a(1:n,0:k), rp(1:4)
    CHARACTER*3 aim
    COMMON/PbPbPARAMS/a,rp,aim
    ! print*, rp
    ! do m = 0, k
    !     print*, a(:,m)
    ! enddo
    ! stop
    phi = 0d0
    if ( r .gt. 2.60d0 ) then
       do i = 1, n
          do m = 0, k
          if ((r .gt. rp(i)).and.(r .le. rp(i+1))) then
             phi = phi + a(i,m)*(r-rp(i+1))**m 
          endif
          enddo
       enddo
    else
       phi = 0.438472d0 - 3.99326d0*(2.60d0-r) + 2.8d0 *(dexp(1.96d0*(2.60d0-r)) - 1d0)
    endif
    pot_EAM_Pb = phi
    return
END FUNCTION pot_EAM_Pb
REAL*8 FUNCTION der_pot_EAM_Pb(r)
    IMPLICIT NONE
    INTEGER i, n, m, k
    PARAMETER(n=3,k=8)
    REAL*8 r, dphi, a(1:n,0:k), rp(1:4)
    CHARACTER*3 aim
    COMMON/PbPbPARAMS/a,rp,aim
    dphi = 0d0
    if ( r .gt. 2.60d0 ) then
       do i = 1, n
          do m = 1, k
          if ((r .gt. rp(i)).and.(r .le. rp(i+1))) then
             dphi = dphi + m*a(i,m)*(r-rp(i+1))**(m-1) 
          endif
          enddo
       enddo
    else
       dphi = 3.99326d0 - 5.488d0 * dexp(1.96d0*(2.60d0-r)) 
    endif
    der_pot_EAM_Pb = dphi
    return
END FUNCTION der_pot_EAM_Pb

REAL*8 FUNCTION funcio(x, C, gamma, delta)
    IMPLICIT NONE
    REAL*8, intent(in) :: x, C, gamma, delta

    funcio = C * dexp(-gamma * (x - 1d0)) / (1d0 + (x - delta)**20d0)

    return
END FUNCTION funcio
REAL*8 FUNCTION factor(x, C, gamma, delta)
    IMPLICIT NONE
    REAL*8, intent(in) :: x, C, gamma, delta

    factor = -(gamma + 20d0*(x - delta)**19d0 / (1d0 + (x - delta)**20d0))

    return
END FUNCTION factor
REAL*8 FUNCTION pot_ZHOU_Pb(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    ! REAL*8 re, fe, rhoe, alpha, beta, A, B, kappa, lambda, eta, FFe
    ! REAL*8, dimension(0:3) :: Fn, F
    REAL*8 funcio, x
    REAL*8 re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn(0:3), F(0:3), eta, FFe
    COMMON /ZhouPARAMS/ re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn, F, eta, FFe

    x = r / re
    pot_ZHOU_Pb = (funcio(x,A,alpha,kappa) - funcio(x,B,beta,lambda))

    return
END FUNCTION pot_ZHOU_Pb
REAL*8 FUNCTION der_pot_ZHOU_Pb(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    ! REAL*8 re, fe, rhoe, alpha, beta, A, B, kappa, lambda, eta, FFe
    ! REAL*8, dimension(0:3) :: Fn, F
    REAL*8 funcio, factor, x
    REAL*8 re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn(0:3), F(0:3), eta, FFe
    COMMON /ZhouPARAMS/ re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn, F, eta, FFe

    x = r / re
    der_pot_ZHOU_Pb = (  funcio(x,A,alpha,kappa) * factor(x,A,alpha,kappa) &
                       - funcio(x,B,beta,lambda) * factor(x,B,beta,lambda)   ) / re

    return
END FUNCTION der_pot_ZHOU_Pb

SUBROUTINE embedding_Pb_Zhou(rhored,FZ,dFdrhored)
    IMPLICIT NONE
    REAL*8, intent(in) :: rhored
    REAL*8, intent(out) :: FZ, dFdrhored
    INTEGER i
    REAL*8 rhon, rho0, rho
    REAL*8 re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn(0:3), F(0:3), eta, FFe
    COMMON /ZhouPARAMS/ re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn, F, eta, FFe

    rhon=0.85d0*rhoe
    rho0=1.15d0*rhoe
    rho = rhored*rhoe

    FZ = 0d0
    dFdrhored = 0d0
    if (rho .lt. rhon) then
        do i = 0, 3
            FZ = FZ + Fn(i)*(rho/rhon-1d0)**dble(i)
            dFdrhored = dFdrhored + dble(i)*Fn(i)*(rho/rhon-1d0)**dble(i-1)
        enddo
        dFdrhored = dFdrhored * rhoe / rhon
    elseif (rho .lt. rho0) then
        do i = 0, 3
            FZ = FZ + F(i)*(rho/rhoe-1d0)**dble(i)
            dFdrhored = dFdrhored + dble(i)*F(i)*(rho/rhoe-1d0)**dble(i-1)
        enddo
        ! dFdrhored = dFdrhored * rhoe / rhoe
    elseif (rho .ge. rho0) then
        FZ = FFe * (1d0 - eta*dlog(rho/rhos)) * (rho/rhos)**eta
        dFdrhored = -FFe * eta**2d0 * (rho/rhos)**(eta-1d0) * dlog(rho/rhos) * rhoe / rhos 
    else
        print*, "ERROR: rho out of bounds"
        stop
    endif
    
    return
END SUBROUTINE embedding_Pb_Zhou
SUBROUTINE electro_Zhou(r,psi,dpsidr)
! returns the electronic density of the LM (Zhou model)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    REAL*8, intent(out) ::  psi, dpsidr 
    REAL*8 funcio, factor
    REAL*8 re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn(0:3), F(0:3), eta, FFe
    COMMON /ZhouPARAMS/ re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn, F, eta, FFe
    
    psi = funcio(r/re, Fe, beta, lambda) / rhoe
    dpsidr = funcio(r/re, Fe, beta, lambda) * factor(r/re, Fe, beta, lambda) / rhoe / re
    return
END SUBROUTINE electro_Zhou

REAL*8 FUNCTION USladek(r)
	IMPLICIT NONE
	INTEGER k
	REAL*8 x, r, V
	REAL*8 a(1:8), De, Re
	REAL*8 aa,bb,cc,dd,ee,ff,gg,hh,ii,jj
    REAL*8 Eh, Rbohr
	COMMON/PbHePARAMS/aa,bb,cc,dd,ee,ff,gg,hh,ii,jj
    COMMON/CONVERS/Eh, Rbohr
	
	Re = aa
	De = hh
	a(1) = bb
	a(2) = cc
	a(3) = dd
	a(4) = ee
	a(5) = ff
	a(6) = gg
	a(7) = ii
    a(8) = jj
	
	x = (r-Re)/Re
	
	V = 1d0
	do k=1,8 
	   V = V + a(k)*x**dble(k)
	enddo
	V = -De*exp(-a(1)*x)*V
	
    USladek = V*1d-6*Eh
	return
END FUNCTION USladek
 
	
REAL*8 FUNCTION dUSladek(r)
	IMPLICIT NONE
	INTEGER k
	REAL*8 x, r, dV
	REAL*8 a(1:8), De, Re
	REAL*8 sumk
	REAL*8 aa,bb,cc,dd,ee,ff,gg,hh,ii,jj
	REAL*8 Eh, Rbohr
	COMMON/PbHePARAMS/aa,bb,cc,dd,ee,ff,gg,hh,ii,jj
	COMMON/CONVERS/Eh, Rbohr
	
!	COMMON/PARAM/Re,a(1),a(2),a(3),a(4),a(5),a(6),De,a(7),a(8)
	
	Re = aa
	De = hh
	a(1) = bb
	a(2) = cc
	a(3) = dd
	a(4) = ee
	a(5) = ff
	a(6) = gg
	a(7) = ii
    a(8) = jj

	x = (r-Re)/Re
	
	sumk = 0d0
	do k=1,8
	   sumk = sumk + (a(1)*x - k)*a(k)*x**dble(k-1)
	enddo
	dV = De*exp(-a(1)*x)*(a(1) + sumk)

    dUSladek = dV*1d-6*Eh/Re
	
	return
END FUNCTION dUSladek

SUBROUTINE electro(r,psi,dpsidr,p1,p2)
! returns the electronic density of the LM
    IMPLICIT NONE
    REAL*8, intent(in) :: r, p1, p2
    REAL*8, intent(out) ::  psi, dpsidr 
    
    psi = p1 * dexp(-p2*r)
    dpsidr = -p2 * psi

    return
END SUBROUTINE electro

SUBROUTINE embedding_Li(rho,F,dFdrho)
    IMPLICIT NONE
    REAL*8, intent(in) :: rho
    REAL*8, intent(out) :: F, dFdrho
    REAL*8 param(0:6), a(1:7), b(1:7), c(1:7), m
    INTEGER i
    COMMON/LiLiPARAMS/param,a,b,c,m
    
    F = 0d0
    if (rho .gt. param(6)) then ! hauria de ser .ge., però llavors dona NaN quan rho=param(6) i m=0 (cas Awad). Com que la funcio i la derivada són contínues, és totalment lícit.
        F = a(7) + b(7)*(rho-param(6)) + c(7)*(rho-param(6))**m
        dFdrho = b(7) + m*c(7)*(rho-param(6))**(m-1d0)
    else if (rho .ge. param(1)) then
        F = a(1) + c(1)*(rho-param(0))**2d0
        dFdrho = 2d0*c(1)*(rho-param(0))
    else if (rho .lt. param(5)) then
        F =(a(6)+b(6)*(rho-param(5)) + c(6)*(rho-param(5))**2d0)
        dFdrho = b(6) + 2d0*c(6)*(rho-param(5))
    else
        do i = 2,5
        if ((rho.ge.param(i)).and.(rho.lt.param(i-1))) then
            F = a(i) + b(i)*(rho-param(i-1)) + c(i)*(rho-param(i-1))**2d0
            dFdrho = b(i) + 2d0*c(i)*(rho-param(i-1))
        endif
        enddo
    endif
    return
END SUBROUTINE embedding_Li


SUBROUTINE embedding_Pb(rho,F,dFdrho)
    IMPLICIT NONE
    REAL*8, intent(in) :: rho
    REAL*8, intent(out) :: F, dFdrho
    REAL*8 param(0:6), a(1:6), b(1:6), c(1:6), m
    REAL*8 dif, frac
    INTEGER i
    
    m = 1.60d0 
    param(0) = 1d0 
    param(1) = 0.90d0 
    param(2) = 0.81d0 
    param(3) = 0.77d0 
    param(4) = 0.71d0 
    param(5) = 0.46d0 
    param(6) = 1.40d0

    a(1) = -1.5186d0 
    a(2) = -1.500978d0 
    a(3) = -1.469082d0 
    a(4) = -1.4485844d0 
    a(5) = -1.425158d0 
    a(6) = -1.297423d0 

    b(1) = 0d0 
    b(2) = -0.35244d0 
    b(3) = -0.35244d0 
    b(4) = -0.67244d0 
    b(5) = -0.10844d0 
    b(6) = -0.91344d0 

    c(1) = 1.7622d0 
    c(2) = 0.0000d0 
    c(3) = 4.00d0 
    c(4) = -4.70d0 
    c(5) = 1.61d0
    c(6) = -5.70d0
    
    F = 0d0
                
    if (rho .le. param(5)) then
        dif = rho - param(5)
        frac = rho / param(5)
        F =(a(6)+b(6)*dif+c(6)*dif**2d0)*(2d0*frac-frac**2d0)
        dFdrho = (b(6) + 2d0*c(6)*dif)*(2d0*frac-frac**2d0) &
        + (a(6)+b(6)*dif+c(6)*dif**2d0)*2d0/param(5)*(1d0-frac)
    else if (rho .ge. param(1)) then
        dif = rho - param(0)
        F = a(1)+c(1)*dif**2d0
        dFdrho = 2d0*c(1)*dif
    else
        do i = 2,5
        if ((rho.le.param(i-1))) then
            dif = rho - param(i-1)
            F = a(i) + b(i)*dif + c(i)*dif**2d0
            dFdrho = b(i) + 2d0*c(i)*dif
        endif
        enddo
    endif
    
    return
END SUBROUTINE embedding_Pb

SUBROUTINE embedding(tipus, rho, F, dFdrho)
    IMPLICIT NONE
    CHARACTER*2, intent(in) :: tipus
    REAL*8, intent(in) :: rho
    REAL*8, intent(out) :: F, dFdrho

    if (tipus .eq. "Li") then
        call embedding_Li(rho,F,dFdrho)
    elseif (tipus .eq. "Pb") then
        call embedding_Pb(rho,F,dFdrho)
    else
        F = 0d0
        dFdrho = 0d0
    endif
    return
END SUBROUTINE embedding

REAL*8 FUNCTION CORRECCIO_YUKAWA(r, pot, i, j, X1)
    IMPLICIT NONE
    INTEGER, intent(in) :: i, j
    REAL*8, intent(in) :: r, X1
    CHARACTER*9, intent(in) :: pot
    REAL*8, parameter :: lambda=1.10d0
    REAL*8, dimension(1:2) :: q
    REAL*8 aij

    if (X1 .eq. 0.9d0) then
        q(1) = 0.11d0
        q(2) = -0.99d0
    elseif (X1 .eq. 0.8d0) then
        q(1) = 0.13d0
        q(2) = -0.52d0
    elseif (X1 .eq. 0.7d0) then
        q(1) = 0.141d0
        q(2) = -0.329d0
    elseif (X1 .eq. 0.6d0) then
        q(1) = 0.0926d0
        q(2) = -0.139d0
    elseif (X1 .eq. 0.5d0) then
        q(1) = 0.075d0
        q(2) = -0.075d0
    else
        q(1) = 0d0
        q(2) = 0d0
    endif

    aij = q(i) * q(j) * 14.3996d0 !14.3996 factor conversió (computational units --> eV)

    if (pot .eq. "potencial") then
        CORRECCIO_YUKAWA = (aij/r) * dexp(-lambda*r)
    elseif (pot .eq. "derivada_") then
        CORRECCIO_YUKAWA = - (aij/r) * dexp(-lambda*r) * (1d0/r + lambda)
    else
        print*, "ERROR: la variable 'pot' ha de ser 'potencial' o 'derivada_' instead of "//pot//"..."
        stop
    endif
    return
END FUNCTION CORRECCIO_YUKAWA