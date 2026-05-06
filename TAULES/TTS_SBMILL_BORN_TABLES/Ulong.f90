! REAL*8 FUNCTION U_long(x,alpha,lambda,gamma,b,A,c)
! 	IMPLICIT NONE
! 	REAL*8 x
! 	REAL*8 a1, a2, a3
! 	REAL*8 alpha, lambda, gamma, b, A, c(6:16)
! 	REAL*8 U_vdw
! 	REAL*8 sumn, sumk, factk
! 	INTEGER n, k
! !	COMMON/PARAM1/a1, a2, a3, alpha
! !	COMMON/PARAM2/gamma,lambda,A,b,c
	
	
! 	sumn = 0d0
! 	do n = 3, 8
! 	   sumk = 1d0
! 	   factk = 1d0
! 	   do k = 1, 2*n
! 	      factk = factk * k
! 	      sumk = sumk + (b*x)**k / factk
! 	   enddo
! 	   sumn = sumn + (1d0 - dexp(-b*x)*sumk)*c(2*n)/(x**(2*n))
! 	enddo

! 	U_vdw = A * x**gamma * dexp(-lambda*x) - sumn
	
! 	U_long = (1d0 - dexp(-alpha*x)) * U_vdw

! 	return
! END FUNCTION U_long

! REAL*8 FUNCTION D_U_long(x,alpha,lambda,gamma,b,A,c)
! 	IMPLICIT NONE
! 	REAL*8 x
! 	REAL*8 a1, a2, a3
! 	REAL*8 alpha, lambda, gamma, b, A, c(6:16)
! 	REAL*8 U_vdw,D_U_vdw_0,D_U_vdw_1,D_U_vdw
! 	REAL*8 sumn, sumk, factk, sumn_1,sumn_2,sumn_3,sumk_0
! 	INTEGER n, k
! !	COMMON/PARAM1/a1, a2, a3, alpha
! !	COMMON/PARAM2/gamma,lambda,A,b,c
	
	
! 	sumn = 0d0
! 	sumn_1 = 0d0
! 	sumn_2 = 0d0
! 	sumn_3 = 0d0

! 	do n = 3, 8
! 	   sumk = 1d0
! 	   sumk_0 = 0.0
! 	   factk = 1d0
! 	   do k = 1, 2*n
! 	      factk = factk * k
! 	      sumk = sumk + (b*x)**k / factk
! 	      sumk_0 = sumk_0 + (b*x)**k / factk / x * K
! 	   enddo
!        sumn = sumn + (1d0 - dexp(-b*x)*sumk)*c(2*n)/(x**(2*n))
!        sumn_1 = sumn_1 - (1d0 - dexp(-b*x)*sumk)*c(2*n)/(x**(2*n+1)) * 2 * n

!        sumn_2 = sumn_2 +  b * dexp(-b*x) * sumk * c(2*n) / (x**(2*n))

!        sumn_3 = sumn_3 - dexp(-b*x) * sumk_0 * c(2*n)/(x**(2*n))

! 	enddo

!  	U_vdw = A * x**gamma * dexp(-lambda*x) - sumn

! 	D_U_vdw_0 = A * dexp(-lambda*x) * x**gamma * ( gamma / x - lambda)

! 	D_U_vdw_1 = sumn_1 + sumn_2 + sumn_3

!     D_U_vdw = D_U_vdw_0 - D_U_vdw_1

! 	D_U_long = (1d0 - dexp(-alpha*x)) * D_U_vdw + alpha * dexp(-alpha*x) * U_vdw

! 	return
! END FUNCTION D_U_long


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