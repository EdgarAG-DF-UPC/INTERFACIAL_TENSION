PROGRAM PbHe_potential
	IMPLICIT NONE
	INTEGER k,m
	REAL*8 r, V, dV, delr
	REAL*8 a,b,c,d,e,f,g,h,i,j
	CHARACTER*8 LABEL1, LABEL2, LABEL3
	REAL FACTOR1, FACTOR2, FACTOR3
	REAL*8 rs
	CHARACTER*8 r0s
	EXTERNAL V, dV
	REAL*8, parameter :: r_ref=2.5d0
	COMMON/PARAM/a,b,c,d,e,f,g,h,i,j
	
		
	delr = 1d-2
	
	! read(*,*) FACTOR1, FACTOR2
	read(*,*) FACTOR1, FACTOR2, FACTOR3
	write(LABEL1, "(F8.6)") FACTOR1
	write(LABEL2, "(F8.6)") FACTOR2
	   
	open(42,file="PbHe_Born.table")

	write(42,'(a)') "# DATE: 2022-09-07 UNITS: metal CONTRIBUTOR: Edgar"
	write(42,'(a)') "# Potencial Pb-He (model Sladek - dU[TQ]Z)"// &
	" scaled by a linear coupling parameter = "//LABEL2// &
	" + Born-Mayer potential with A=12.483 eV, B=1 \AA multiplied by a linear coupling parameter = "//LABEL1
	write(42,'(a)')
	write(42,'(a)') "SLADEK_PbHe"
	write(42,'(a)') "N 750 R 0.02 15.00"
	write(42,'(a)')	
!	dU[TQ]Z in Fig.1
!	HF(QZ) + CORR(~r^{-3}, TQ)
	a = 4.6874d0
	b = 5.4114d0
	c = -20.7496d0
	d = 56.2663d0
	e = -95.6327d0
	f = 112.157d0
	g = -77.8808d0
	h = 39.4116d0
	i = 28.2574d0
	j = -3.7882d0
	
	do k=1,750
		r = float(k)*2d-2
		write(42,*) k, real(r)+FACTOR3*r_ref*0.5291772d0, &
		real(FACTOR2*V(r)*1d-6*27.211386245988 &
		+ FACTOR1 * 12.483d0 * dexp(-r / 1d0)) &
		, real(-FACTOR2*dV(r) * 1d-6 * 27.211386245988/a &
		+ FACTOR1 * 12.483d0 / 1d0 * dexp(-r / 1d0))
	enddo
	
	close(42)

	call BISECTION(V, 2d0, 7.5d0, 1d-4, rs)
	write(r0s, "(F8.6)") rs
	print*, "V("//r0s//" Å) = 0"
	call BISECTION(dV, 2d0, 7.5d0, 1d-4, rs)
	write(r0s, "(F8.6)") rs
	print*, "V("//r0s//" Å) = Vmin"
	
	
	END PROGRAM PbHe_potential
	
	
	REAL*8 FUNCTION V(r)
	IMPLICIT NONE
	INTEGER k
	REAL*8 x, r
	REAL*8 a(1:8), De, Re
	REAL*8 aa,bb,cc,dd,ee,ff,gg,hh,ii,jj
	COMMON/PARAM/aa,bb,cc,dd,ee,ff,gg,hh,ii,jj
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
	
	V = 1d0
	do k=1,8 
	   V = V + a(k)*x**dble(k)
	enddo
	V = -De*exp(-a(1)*x)*V
	
	return
	END FUNCTION 
 
	
	REAL*8 FUNCTION dV(r)
	IMPLICIT NONE
	INTEGER k
	REAL*8 x, r
	REAL*8 a(1:8), De, Re
	REAL*8 sumk
	REAL*8 aa,bb,cc,dd,ee,ff,gg,hh,ii,jj
	COMMON/PARAM/aa,bb,cc,dd,ee,ff,gg,hh,ii,jj
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
	
	return
	END FUNCTION 

	SUBROUTINE BISECTION(funcio,AA,BB,eps,arrel)
		IMPLICIT NONE
		REAL*8 funcio
		REAL*8, intent(in) :: AA, BB, eps
		REAL*8, intent(out) :: arrel
		INTEGER niter, i, nmax
		REAL*8 dif
		REAL*8 fA, fB, fC
		REAL*8 A, B, C
		
		A = AA
		B = BB
		
		fA = funcio(A)
		fB = funcio(B)
		
		nmax = 10000
		
		! print*, nint(dlog((B-A)/eps)/dlog(2d0)) +1
		
		!nmax = nint(dlog((B-A)/eps)/dlog(2d0)) +1
		
		! print*, nmax
		
		do i = 1, nmax
		   dif = dabs(B-A)
		   if (fA*fB .lt. 0d0) then
			  C = 0.5d0*(A+B)
			  fC = funcio(C)
			  if ((fC.eq.0d0).or.(dif.lt.eps)) then
				 niter = i
				 arrel = C
				 exit
			  else if (fC*fA .lt. 0d0) then
				 B = C
			  else if (fC*fB .lt. 0d0) then
				 A = C
			  endif
		   endif
		enddo
		
		print*, "Ha convergit en", niter, "iteracions"
		
		return
	END SUBROUTINE BISECTION
