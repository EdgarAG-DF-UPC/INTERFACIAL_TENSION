! REAL*8 FUNCTION U_short(x,a1,a2,a3,alpha)
! 	IMPLICIT NONE
! 	REAL*8 x
! 	REAL*8 a1, a2, a3, alpha
! !	COMMON/PARAM1/a1, a2, a3, alpha
	
! 	U_short = (1d0/x) * (1d0 + a1*x + a2*x**2 + a3*x**3)*dexp(-alpha*x)
	
! 	return
! 	END FUNCTION U_short
	
! 	REAL*8 FUNCTION D_U_short(x,a1,a2,a3,alpha)
! 	IMPLICIT NONE
! 	REAL*8 x
! 	REAL*8 a1, a2, a3, alpha, U_short
! !	COMMON/PARAM1/a1, a2, a3, alpha
	
!     D_U_short = - u_short(x,a1,a2,a3,alpha) * (1.D0/x + alpha) &
!     + (1d0/x) * ( a1 + 2.d0 * a2 * x + 3.d0* a3 * x**2) *dexp(-alpha*x)
	
! 	return
! END FUNCTION D_U_short

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