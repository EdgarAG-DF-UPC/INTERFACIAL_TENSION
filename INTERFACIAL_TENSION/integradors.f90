REAL*8 FUNCTION INT_TRAP(f,x,Nx)
    IMPLICIT NONE
    INTEGER, intent(in) :: Nx
    REAL*8, dimension(1:Nx), intent(in) :: f, x
    REAL*8 x0, x1, Dx, sum_k
    INTEGER k

    x0 = x(1)
    x1 = x(Nx)
    Dx = (x1 - x0) / dble(Nx-1) 
    INT_TRAP = 0d0
    sum_k = 0d0

    do k = 2, Nx-1
        sum_k = sum_k + f(k)
    enddo

    INT_TRAP = Dx * (.5d0*f(1) + sum_k + .5d0*f(Nx))

    return
END FUNCTION INT_TRAP

!--------------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE Bisection(xA,xB,eps,fun,niter,xarrel)
    IMPLICIT NONE
    real*8, intent(in) :: xA, xB, eps
    real*8 fun
    integer, intent(out) :: niter
    real*8, intent(out) :: xarrel
    real*8 A,B,C
    integer i, nmax
    real*8 fA, fB, fC
    
    A = xA
    B = xB
    fA = fun(A)
    fB = fun(B)
    nmax = nint(log((B-A)/eps)/log(2.))+1
    
    
    do i = 1, nmax
        if (fA*fB .lt. 0) then
            C = (A+B)/2d0
            fC = fun(C)
            if (fc .eq. 0d0) then
                niter = i
                xarrel = c
                exit
            elseif (fC*fA .lt. 0) then
                B = C
                fB = fC
            elseif (fC*fB .lt. 0) then
                A = C
                fA = fC
            endif
        else
            print*, "PROBLEMA: fa*fb >0"
            exit
        endif
        
        if (dabs(B-A) .lt. eps) then
            niter = i
            xarrel = C
            exit
        endif
    enddo
    
    if (dabs(B-A) .gt. eps) print*, "El mètode no ha convergit"
    
    return
END SUBROUTINE BISECTION