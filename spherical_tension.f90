! Edgar Alvarez Galera
! Tensió superficial - sistema Li+Pb+He --- EAM (Awad) + TTS
! Data de creació: 19 d'abril del 2024
! Última modificació: 26 / 09 / 2025

PROGRAM ST
    USE omp_lib
    USE CONSTANTS
    USE PARAMETRES_LJ
    USE PARAMETRES_EAM
    USE PARAMETRES_SLADEK
    USE PARAMETRES_ZHOU
    IMPLICIT NONE    
    INTERFACE
        subroutine remove_integer(list, element)
            integer, dimension(:), allocatable, intent(inout) :: list
            integer, intent(in) :: element
        end subroutine remove_integer
        subroutine add_integer(list, element)
            integer, dimension(:), allocatable, intent(inout) :: list
            integer, intent(in) :: element
        end subroutine add_integer
        integer function findloc(list, element)
            integer, dimension(:), allocatable, intent(in) :: list
            integer, intent(in) :: element
        end function findloc
        subroutine ordena(list)
            integer, dimension(:), allocatable, intent(inout) :: list
        end subroutine ordena
    END INTERFACE
    TYPE list
        integer, dimension(:), allocatable :: nbr, group
    END TYPE list
    TYPE peratom
        integer tipus, id
        character*2 nom
        real*8 massa
        real*8, dimension(3) :: posicio, forca
    END TYPE peratom
    TYPE perelem
        real*8 massa
        real*8, dimension(:), allocatable :: edens
        character*2 nom
    END TYPE perelem
    TYPE thermo_data
        integer step
        real*8 temperatura, pressio, volum, enertot, potener, kinener, enertail
    END TYPE thermo_data
    TYPE(thermo_data) RTD
    REAL*8 INT_TRAP, INTERF_ENER, Rs, INTEGRAL
    REAL*8 distancia, dist
    REAL*8 Ekin, Epot, Epi
    REAL*8 COM
    INTEGER N, NHe, NLi, NPb 
    INTEGER timesteps, snapshots, every, sni, snf, tsi, tsf
    CHARACTER*15 snap
    REAL*8 Vol, Temp, Press, Box, xl, xh, yl, yh, zl, zh, boxl(1:3), boxh(1:3)
    CHARACTER*3 VAR ! caràcter inútil
	REAL*8 psiij, dpsiij, Fi, dFi, Fj, dFj!.....
    REAL*8 dGdt, sum_Fij_rij, sum_Fr, Press_VIR ! Teorema del Virial
    INTEGER force_field !0 --> LJ; 1 --> EAM + TTS + Sladek
    LOGICAL FULLLJ, EAMLJ, EAMTTSS
!   INDEXS:
    INTEGER i, j, k, l, m, p, q, s, t, ti, tj, n1, n2, n3, nt(3)
    INTEGER, dimension(:,:), allocatable :: PAR, PARL
!   COMPTADORS (Àtoms totals, àtoms a la droplet, timesteps, etc.):
    INTEGER Nats, Nd, ts
    REAL*8 count
    INTEGER c1, c2, c0, ctot, c1i, c2i, c0i, ctoti
!   COORDENADES, FORCES I VELOCITATS:
    REAL*8, dimension(3) :: posicio, velocitat, force, poscm
    REAL*8, dimension(:), allocatable :: edens, edensp
    REAL*8, dimension(:,:), allocatable :: pos, F, pos_LM, pos_He, Fori
!   QUANTITATS "PER ATOM":
    TYPE(peratom), dimension(:), allocatable :: particle, LM, He
    TYPE(perelem) liquid(2)
    REAL*8, dimension(3), parameter :: massa=(/ massa_Li, massa_He, massa_Pb /)
    REAL*8, dimension(:), allocatable :: masscl
!   PARÀMETRES DELS POTENCIALS NO DEFINITS ALS MODULES
    REAL*8 PSI(0:6), Ae(1:7), Be(1:7), Ce(1:7), me ! EAM Li-Li
    REAL*8 aaa(1:3,0:8), rp(1:4) ! EAM Pb-Pb
    CHARACTER*3 aim ! ---------     "    "    -------- 
    REAL*8 a1(1:3), a2(1:3), a3(1:3) ! Paràmetres dels models TTS
    REAL*8 alpha(1:3), lambda(1:3), gamma(1:3) ! Paràmetres dels models TTS
    REAL*8 b(1:3), A(1:3), c(1:3,6:16) ! Paràmetres dels models TTS
    REAL*8 Re(1:3), De(1:3), ZA, ZB ! Paràmetres dels models TTS
    REAL*8 U_short, U_long, dU_long, dU_short ! Contribucions al model TTS
    REAL*8 Us, Ul, delr, dUl, dUs, r0 ! Contribucions al model TTS    
    REAL*8 reZ, feZ, rhoeZ, rhosZ, alphaZ, betaZ, AZ, BZ, kappaZ, lambdaZ, FnZ(0:3), F0Z(0:3), etaZ, FFeZ
!   MAGNITUDS RADIALS:
    REAL*8, dimension(:), allocatable :: r, gdr, z, r_ann, dens, num, num_he, num_li, num_pb, pn, pk, pv, pvaver, pvp
    INTEGER nr0, nr1, nr
    REAL*8 Dr
!   CLUSTERS:
    INTEGER, dimension(:), allocatable :: ATS_cluster
    REAL*8, dimension(:), allocatable :: coord
    LOGICAL ALL_He_ATS
    PARAMETER(ALL_He_ATS=.False.)
!   RADIS:
    REAL*8 MSR, Rd, Rd2, minDISTTOCM, maxDISTTOCM
!   PERIODIC BOUNDARY CONDITIONS:
    LOGICAL POTS, SETPBC, EAM, Awad, Yukawa
    LOGICAL Zhou_01, Zhou_04!, Zhou
    LOGICAL Belashchenko_Li, Belashchenko_Pb!, Belashchenko
!   OMP:
    INTEGER thread_id
!   CPU TIME and SYSTEM CLOCK:
    REAL temps_ini, temps
    INTEGER*4 :: count_0, count_1, count_rate, count_max
!   PRINT ON EVERY STEP?
    LOGICAL PRINTALL, PRINTSNAP
!   PARENT DIRECTORY:
    CHARACTER*100 parent
!   COMMON BLOCKS:
    COMMON /PbHePARAMS/ aa,bb,cc,dd,ee,ff,gg,hh,ii,jj
    COMMON /PbPbPARAMS/ aaa,rp,aim
    COMMON /LiLiPARAMS/ PSI,Ae,Be,Ce,me
    COMMON /ZhouPARAMS/ reZ, feZ, rhoeZ, rhosZ, alphaZ, betaZ, AZ, BZ, kappaZ, lambdaZ, FnZ, F0Z, etaZ, FFeZ
    ! REAL*8 re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn(0:3), F(0:3), eta, FFe
    COMMON /CAPSA/ boxl,boxh
    COMMON /CONVERS/ Eh, Rbohr
    COMMON /OTHERS/ parent

    liquid(1)%massa = massa_Li
    liquid(2)%massa = massa_Pb
    liquid(1)%nom = "Li"
    liquid(2)%nom = "Pb"
    
    open(14, FILE="st.params")
        read(14,*) N, NHe, NLi, NPb
        read(14,*) Temp, Press
        read(14,*) tsi, tsf, every
        sni = 0
        snf = (tsf-tsi)/every
        snapshots = snf-sni+1
        timesteps = tsf-tsi+1
        read(14,*) Dr, nr0, nr1
        read(14,*) force_field
        read(14,*) POTS, SETPBC
        read(14,*) PRINTALL, PRINTSNAP
    close(14)
    open(14, FILE="parent-directory")
        read(14,*) parent
    close(14)
    !   TRIEM LA FAMILIA DE POTENCIALS:
    call SETFF()
    call PARAMETRES()
    call WRITE_POTS(POTS)
    if (PRINTALL) print*, "ALL POTENTIALS ARE SET!"

    ! FIXEM LES DIMENSIONS DE CADA ARRAY:
    allocate(particle(N))
    allocate(pos(N,3), F(N,3), Fori(N,3))
    ! allocate(pos_LM(1:NLi+NPb,3), pos_He(1:NHe,3))
    allocate(LM(1:NLi+NPb), He(1:NHe))
    if (EAM) then
        allocate(edens(1:NLi+NPb), edensp(1:NLi+NPb))
        allocate(liquid(1)%edens(1:NLi))
        allocate(liquid(2)%edens(1:NPb))
    endif
    allocate(PAR(1:N*(N-1)/2, 2), PARL(1:(NLi+NPb)*(NLi+NPb-1)/2, 2))

    k = 0
    do i = 1, N-1
        do j = i+1, N
            k = k + 1
            PAR(k,1) = i
            PAR(k,2) = j
        enddo
    enddo
    k = 0
    do i = 1, NLi+NPb-1
        do j = i+1, NLi+NPb
            k = k + 1
            PARL(k,1) = i
            PARL(k,2) = j
        enddo
    enddo

    Dr = 0.1d0 !Angstroms
    nr = nint(35d0 / Dr) ! 16.6d0 ~= L/2
    nr0 = 1
    nr1 = nr
    allocate(r(nr0:nr1), gdr(nr0:nr1), z(nr0:nr1), pn(nr0:nr1), pk(nr0:nr1), pv(nr0:nr1), dens(nr0:nr1), num(nr0:nr1))
    allocate(num_he(nr0:nr1),num_li(nr0:nr1),num_pb(nr0:nr1))
    allocate(pvaver(nr0:nr1), pvp(nr0:nr1))
    allocate(r_ann(nr0-1:nr1))
    
    r = 0d0
    r_ann = 0d0
    do k = nr0-1, nr1
        r_ann(k) = Dr*dble(k)
    enddo
    do k = nr0, nr1
        r(k) = r_ann(k) - Dr/2d0
    enddo



!   INICIALITZEM MAGNITUDS A PROMITJAR:
    num = 0
    num_li = 0d0
    num_pb = 0d0
    num_he = 0d0
    Rd = 0d0
    Rd2 = 0d0
    pvaver = 0d0
!   OBRIM ELS ARXIUS DE LECTURA
    open(15, FILE="xyz.out")
    open(16, FILE="simulation_box.out")
    read(16,*)
    open(17, FILE="center_mass.out")
    read(17,*)
    open(18, FILE="thermo_data.out")
    read(18,*)
    temps_ini = 0d0

    do s = sni, snf
        if (PRINTALL) then
            write(*,*)
            write(*,*)
            write(*,*) "============================================================================================"
            write(snap, "(I15)") s
            write(*,*) "================================= SNAPSHOT "//snap//" ================================="
            write(*,*)
            
            CALL SYSTEM_CLOCK(count_1, count_rate, count_max)
            write(*,*)
            write(*,"(A51)") "TEMPS D'EXECUCIÓ (segons)"
            write(*,"(A51)") "---------------------------------------------------"
            write(*,"(A15, A3, A15, A3, A15)") "", " | ", "SYSTEM CLOCK", " | ", "ACCUMULATED"
            write(*,"(A51)") "----------------|-----------------|----------------"
            write(*,"(A15, A3, F15.8, A3, F15.8)") "TIME (s)", " | ", float(count_1)/float(count_rate)/1e6, &
            " | ", temps_ini + float(count_1)/float(count_rate)/1e6
            write(*,"(A51)") "---------------------------------------------------"
            write(*,*)
            temps_ini = temps_ini + float(count_1)/float(count_rate)/1e6
            count_0 = count_1
        endif

        read(16,*) ts, xl, xh, yl, yh, zl, zh
        Vol = (xh-xl)*(yh-yl)*(zh-zl)
        Box = (Vol)**(1d0/3d0)
        boxl = (/ xl, yl, zl /)
        boxh = (/ xh, yh, zh /)
        if (PRINTALL) then
            write(*,"(A51)") "SIMULATION BOX COORDINATES (in Ångström)"
            write(*,"(A51)") "---------------------------------------------------"
            write(*,"(A15, A3, A15, A3, A15)") "Axis", " | ", "LOW", " | ", "HIGH"
            write(*,"(A15, A3, A15, A3, A15)") "------", "-|-", "---------------", "-|-", "---------------"
            write(*,"(A15, A4, F15.8, A3, F15.8)") "x (Å)", " | ", xl, " | ", xh
            write(*,"(A15, A4, F15.8, A3, F15.8)") "y (Å)", " | ", yl, " | ", yh
            write(*,"(A15, A4, F15.8, A3, F15.8)") "z (Å)", " | ", zl, " | ", zh
            write(*,"(A51)") "---------------------------------------------------"
            write(*,"(A34, F15.8, A3)") "BOX WIDTH: ", box, "Å"
            write(*,"(A51)") "---------------------------------------------------"
        endif

        read(15,*) Nats
        if (Nats .ne. N) then
            print*, "El nombre d'àtoms al fitxer xyz no coincideix amb l'indicat al fitxer de configuració..."
            print*, "--> ", Nats
            print*, "--> ", N
            print*, "El programa s'aturarà automàticament."
            open(78, file="LOST_ATOMS")
            write(78, *) "WARNING! Some atoms may have been lost..."
            write(78, *) "--> ", Nats
            write(78, *) "--> ", N
            nt(:) = 0
            read(15,*)
            do i = 1, Nats
                read(15,*) k, pos(i,:)
                nt(k) = nt(k) + 1 
            enddo
            write(78, *) nt(1), " vs ", NLi, " Li atoms"
            write(78, *) nt(3), " vs ", NPb, " Pb atoms"
            write(78, *) nt(2), " vs ", NHe, " He atoms"
            close(78)
            stop
        endif
        read(15,*)
        n1 = 0
        n2 = 0
        n3 = 0
        do i = 1, N
            read(15,*) k, pos(i,:)
            particle(i)%id = i
            if (SETPBC) call PBC(pos(i,:)) !No hauria de ser necessari, però:
            ! NOTE
            !Because periodic boundary conditions are enforced only on timesteps when neighbor lists are rebuilt, the coordinates of an atom written to a dump file may be slightly outside the simulation box. Re-neighbor timesteps will not typically coincide with the timesteps dump snapshots are written. See the dump_modify pbc command if you wish to force coordinates to be strictly inside the simulation box.
            ! (see https://docs.lammps.org/dump.html)
            if (k .eq. 1) then
                particle(i)%nom = "Li"
                n1 = n1 + 1
                ! pos_LM(i-n2,:) = pos(i,:) ! i-n2 = n1+n3
                LM(i-n2)%id = i
                LM(i-n2)%posicio = pos(i,:)
                LM(i-n2)%massa = massa_Li
                LM(i-n2)%nom = "Li"
                LM(i-n2)%tipus = 1
            elseif (k .eq. 2) then
                particle(i)%nom = "He"
                n2 = n2 + 1
                ! pos_He(n2,:) = pos(i,:) ! n2 = i-n1-n3
                He(n2)%id = i
                He(n2)%posicio = pos(i,:)
                He(n2)%massa = massa_He
                He(n2)%nom = "He"
                He(n2)%tipus = 2
            elseif (k.eq. 3) then
                particle(i)%nom = "Pb"
                n3 = n3 + 1
                ! pos_LM(i-n2,:) = pos(i,:) ! i-n2 = n1+n3
                LM(i-n2)%id = i
                LM(i-n2)%posicio = pos(i,:)
                LM(i-n2)%massa = massa_Pb
                LM(i-n2)%nom = "Pb"
                LM(i-n2)%tipus = 3
            endif
            particle(i)%tipus = k
            particle(i)%massa = massa(k)
            particle(i)%posicio = pos(i,:)
        enddo
        if ((n1 .ne. NLi) .or. (n2 .ne. NHe) .or. (n3 .ne. NPb)) then
            print*, "ERROR while reading xyz file"
            print*, n1, "-->", NLi
            print*, n2, "-->", NHe
            print*, n3, "-->", NPb
            stop
        endif

        read(18,*) RTD ! RTD%step, RTD%temperatura, RTD%pressio, RTD%volum, RTD%enertot, &

! Identificació de clústers:
        if (ALL_He_ATS) then
            if (.not. allocated(ATS_cluster)) then
                allocate(ATS_cluster(NHe))
                do i = 1, NHe
                    ATS_cluster(i) = NLi + NPb + i
                enddo
            endif
            Nd = NHe
        else
            call CLUSTERS(4.2d0,Nd,ATS_cluster) ! 1.542*sig_HeHe
            call ORDENA(ATS_cluster)
        endif

        if (PRINTALL) then
            read(17,*) ts, poscm(1), poscm(2), poscm(3)
            write(*,"(A75)") "CENTER OF MASS (in Ångström)"
            write(*,"(A75)") "---------------------------------------------------------------------------"
            write(*,"(A30, A3, A12, A4, A12, A4, A12)") "", " | ", "x (Å)", " | ", "y (Å)", " | ", "z (Å)"
            write(*,"(A75)") "-------------------------------|--------------|--------------|-------------"
            write(*,"(A30, A3, F12.8, A3, F12.8, A3, F12.8)") "all He ats (LAMMPS, whole box)", " | ", &
            poscm(1), " | ", poscm(2), " | ", poscm(3)
        endif
        

        allocate(coord(size(ATS_cluster)), masscl(size(ATS_cluster)))
        do p = 1,3
            do k = 1, size(ATS_cluster)
                i = ATS_cluster(k)
                coord(k) = He(i)%posicio(p) - boxl(p)!coord(k) = pos_He(i,p) - boxl(p)
                masscl(k) = massa_He
            enddo            
            poscm(p) = COM(coord, size(ATS_cluster), Box, masscl, .True.) + boxl(p) !particle(NLi+NPb+1:N)%massa
        enddo
        call PBC(poscm)
        deallocate(coord, masscl)
        if (PRINTALL) then
            write(*,"(A30, A3, F12.8, A3, F12.8, A3, F12.8)") "He ats in main cluster", " | ", &
            poscm(1), " | ", poscm(2), " | ", poscm(3)
            write(*,"(A75)") "---------------------------------------------------------------------------"
        endif

        MSR = 0d0
        minDISTTOCM = 1e8
        maxDISTTOCM = 0d0
        do k = 1, Nd
            i = ATS_cluster(k)
            ! MSR = MSR + distancia(pos_He(i,:),poscm(:), .True.)**2d0
            ! if (distancia(pos_He(i,:),poscm(:), .True.) .gt. maxDISTTOCM) maxDISTTOCM = distancia(pos_He(i,:),poscm(:), .True.)
            ! if (distancia(pos_He(i,:),poscm(:), .True.) .lt. minDISTTOCM) minDISTTOCM = distancia(pos_He(i,:),poscm(:), .True.)
            MSR = MSR + distancia(He(i)%posicio,poscm(:), .True.)**2d0
            if (distancia(He(i)%posicio, poscm(:), .True.) .gt. maxDISTTOCM) maxDISTTOCM = distancia(He(i)%posicio,poscm(:), .True.)
            if (distancia(He(i)%posicio, poscm(:), .True.) .lt. minDISTTOCM) minDISTTOCM = distancia(He(i)%posicio,poscm(:), .True.)
        enddo
        MSR = MSR / dble(Nd)
        Rd = Rd + dsqrt(MSR)
        Rd2 = Rd2 + dsqrt(MSR)*dsqrt(MSR)
        if (PRINTALL) then
            write(*,"(A20, F11.8, A4)") "MEAN SQUARED RADIUS:", dsqrt(MSR), "Å"
            write(*,"(A20, F11.8, A4)") "MIN DIST TO CM:", minDISTTOCM, "Å"
            write(*,"(A20, F11.8, A4)") "MAX DIST TO CM:", maxDISTTOCM, "Å"
            write(*,"(A75)") "---------------------------------------------------------------------------"
            write(*,*)
        endif


        pn = 0d0
        pk = 0d0
        pv = 0d0
        dens = 0d0

! CÀLCUL DEL TERME CINÈTIC: pk = rho(r) k T
        do i = 1, N
                k = floor(distancia(pos(i,:),poscm(:),.True.)/Dr) + 1
                if (k .le. nr) then
                    num(k) = num(k) + 1d0
                    dens(k) = dens(k) + 1d0
                    if (particle(i)%tipus .eq. 1) then
                        num_li(k) = num_li(k) + 1d0
                    elseif (particle(i)%tipus .eq. 2) then
                        num_he(k) = num_he(k) + 1d0
                    elseif (particle(i)%tipus .eq. 3) then
                        num_pb(k) = num_pb(k) + 1d0
                    else
                        print*, "Aixo no hauria de passar mai..."
                    endif
                endif
        enddo
        dens(:) = dens(:) / (4d0 * pi * Dr * (r(:)*r(:) + Dr*Dr/12d0))
        pk(:) = kB * RTD%temperatura * dens(:)
        Ekin = (3d0/2d0) * dble(N) * kB * RTD%temperatura
        if (PRINTALL) then
            write(*,"(A60)") "THERMODYNAMIC PROPERTIES"
            write(*,"(A60)") "------------------------------------------------------------------------------"
            write(*,"(A22, A3, A16, A3, A16)") "", " | ", "", " | ", "LAMMPS"
            write(*,"(A22, A3, A16, A3, A16)") "-------------------", "-|-", "----------------", "-|-", "----------------"
            write(*,"(A22, A3, F16.8, A3, F16.8)") "TEMPERATURE (K)", " | ", temp, " | ", RTD%temperatura
            write(*,"(A22, A5, F16.8, A3, F16.8)") "VOLUME (Å³)", " | ", Vol, " | ", RTD%volum
            write(*,"(A22, A3, F16.8, A3, F15.8)") "KINETIC ENERGY (eV)", " | ", Ekin, " | ", RTD%kinener
        endif
        

! CÀLCUL DEL TERME CONFIGURACIONAL: pv
        ctot = 0
        c0 = 0
        c1 = 0
        c2 = 0
        if (EAM) then
            edens(:) = 0d0
            liquid(1)%edens(:) = 0d0
            liquid(2)%edens(:) = 0d0
          !$OMP PARALLEL PRIVATE(psiij,dpsiij,i,j,ti,tj,edensp) SHARED(edens)
            edensp(:) = 0d0
          !$OMP DO
            do k = 1, (NLi+NPb)*(NLi+NPb-1)/2
                i = PARL(k,1)
                j = PARL(k,2)
                ti = LM(i)%tipus
                tj = LM(j)%tipus
                call electro(distancia(LM(i)%posicio,LM(j)%posicio,.True.),psiij,dpsiij,p1(tj),p2)!call electro(distancia(pos_LM(i,:),pos_LM(j,:),.True.),psiij,dpsiij,p1(tj),p2)
                edensp(i) = edensp(i) + psiij !contribució de l'àtom j a la densitat electrònica a la posició de i
                call electro(distancia(LM(i)%posicio,LM(j)%posicio,.True.),psiij,dpsiij,p1(ti),p2)!call electro(distancia(pos_LM(i,:),pos_LM(j,:),.True.),psiij,dpsiij,p1(ti),p2)
                edensp(j) = edensp(j) + psiij !contribució de l'àtom i a la densitat electrònica a la posició de j
            enddo
           !$OMP END DO
           !$OMP BARRIER
           !$OMP CRITICAL
                edens(:) = edens(:) + edensp(:)            
           !$OMP END CRITICAL            
           !$OMP END PARALLEL
            
        endif
        Epot = 0d0
        sum_Fr = 0d0
        F = 0d0
        
        !$OMP PARALLEL PRIVATE(Epi,pvp,Fi,dFi,sum_Fij_rij,Fori,i,j,ctoti,c1i,c2i,c0i) SHARED(Epot,pv,sum_Fr,F,ctot,c1,c2,c0)
        if (EAM) then
            Epi = 0d0
            sum_Fij_rij = 0d0
            !$OMP DO
            do i = 1, N !i = 1, NLi+NPb
                if (particle(i)%nom .eq. "Li") then
                    call embedding_Li(edens(i),Fi,dFi)
                    Epi = Epi + Fi
                elseif (particle(i)%nom .eq. "Pb") then
                    call embedding_Pb(edens(i),Fi,dFi)
                    Epi = Epi + Fi
                endif
            enddo
            !$OMP END DO
            !$OMP CRITICAL
            Epot = Epot + Epi
            !$OMP END CRITICAL
            !$OMP BARRIER
        endif
        
        Epi = 0d0
        pvp = 0d0
        Fori = 0d0
        ctoti = 0
        c0i = 0
        c1i = 0
        c2i = 0
        !$OMP DO
        do k = 1, N*(N-1)/2
        ! do i = 1, N-1
        !     do j = i+1, N
            i = PAR(k,1)
            j = PAR(k,2)
            if (distancia(pos(i,:),pos(j,:),.True.) .le. 9.01) then
                call P_CONF_PARALLEL(i, j, pvp, c0i, c1i, c2i, ctoti)
                call VIRIAL_PARALLEL(i, j, sum_Fij_rij, Fori)
                Epi = Epi + POTENCIAL(PAIR(i,j), distancia(pos(i,:),pos(j,:),.True.))
            endif
        enddo
        !$OMP END DO
        !$OMP BARRIER
        !$OMP CRITICAL
            Epot = Epot + Epi
            pv = pv + pvp
            sum_FR = sum_FR + sum_Fij_rij
            F = F + Fori
            ctot = ctot + ctoti
            c0 = c0 + c0i
            c1 = c1 + c1i
            c2 = c2 + c2i
        !$OMP END CRITICAL
        !$OMP BARRIER
        !$OMP END PARALLEL
        
        if (PRINTALL) then
            write(*,"(A22, A3, F16.8, A3, F15.8)") "POTENTIAL ENERGY (eV)", " | ", Epot , " | ", RTD%potener
            write(*,"(A60)") "------------------------------------------------------------------------------"
            write(*,*)
            write(*,*)
            write(*,"(A11, A11, A11)") "# talls", "# cops", "%"
            write(*,"(A33)") "---------------------------------"
            write(*,"(I11, I11, F11.3)") 0, c0, float(c0)/float(ctot)*100d0
            write(*,"(I11, I11, F11.3)") 1, c1, float(c1)/float(ctot)*100d0
            write(*,"(I11, I11, F11.3)") 2, c2, float(c2)/float(ctot)*100d0
            write(*,"(A33)") "---------------------------------"
            write(*,"(A11, I11, A11)") "TOTAL", ctot, "talls"
            print*, ""
            print*, ""
            
            print*,"  ============  "
            write(*,"(A, F6.2, A)") "||", dble(s-sni) / dble(snf-sni) * 100d0, "%  DONE ||"
            print*,"  ============  "
            print*, ""
            print*, ""
            print*, "============================================================================================"
        endif
        if (PRINTSNAP) then
            write(*, "(F4.1, A)") dble(s-sni) / dble(snf-sni) * 100d0, "% DONE"
        endif

        if (ctot .ne. c0+c1+c2) stop
        pvaver = pvaver + pv   
        
    enddo
    close(15)

! PROMITJOS:
    !----> Nombre d'àtoms
    num = num / dble(snf-sni)
    num_Li = num_Li / dble(snf-sni)
    num_He = num_He / dble(snf-sni)
    num_Pb = num_Pb / dble(snf-sni)
    !---> Terme configuracional:
    pvaver = pvaver / dble(snf-sni)
    count = 0
    open(97,FILE="r_rho_pv.dat")
    do k = nr0, nr1
        dens(k) = num(k) / (4d0 * pi * Dr * (r(k)*r(k) + Dr*Dr/12d0))
        pvaver(k) = pvaver(k) / (4d0*pi*r(k)*r(k)*r(k))
        write(97,*) r(k), num_li(k), num_pb(k), num_he(k), kB*Temp*dens(k), pvaver(k)
        count = count + num(k)
    enddo
    close(97)

    !----> Radi (arrel del MSR)
    Rd = Rd / dble(snf-sni)
    Rd2 = Rd2 / dble(snf-sni)
    print*, "RADI (PROMIG - MSR):  (", Rd,  "+-", dsqrt((Rd2 - Rd*Rd)/dble(snf-sni)), ")  Å"
    print*, "VARIANCE: <Rd²> - <Rd>² = ",  Rd2, "-", Rd*Rd, " = ", Rd2-Rd*Rd 
    print*, "RECOMPTE:", count

    pn(:) = kB*Temp*dens(:) + pvaver(:)
    print*, ""
    print*, ""
    print*, ""
    print*, ""
    print*, "---------------- CÀLCUL FINAL ----------------"
    print*, "----------------------------------------------"
    call CALCUL_TENSIO_SIMETRIA_ESFERICA("r_rho_pv.dat")
    print*, "----------------------------------------------"
    print*, "----------------------------------------------"
    
CONTAINS
    SUBROUTINE PARAMETRES()
        IMPLICIT NONE
        if (BELASHCHENKO_Li) then
            print*, "Belashchenko EAM"
            me = 1.5d0 !Belashchenko/Fraile
            PSI(0) = 1d0 !Belashchenko/Fraile
            PSI(1) = 0.900d0 !Belashchenko/Fraile
            PSI(2) = 0.840d0 !Belashchenko
            PSI(3) = 0.700d0 !Belashchenko/Fraile
            PSI(4) = 0.550d0 !Belashchenko
            PSI(5) = 0.350d0 !Belashchenko
            PSI(6) = 1.100d0 !Belashchenko/Fraile

            Ae(1) = -0.8948d0 !Belashchenko
            Ae(2) = -0.894474d0 !Belashchenko
            Ae(3) = -0.887963d0 !Belashchenko
            Ae(4) = -0.878482d0 !Belashchenko
            Ae(5) = -0.850369d0 !Belashchenko
            Ae(6) = -0.800385d0 !Belashchenko
            Ae(7) = -0.894474d0 !Belashchenko

            Be(1) = 0d0 !Belashchenko/Fraile
            Be(2) = -0.006520d0 !Belashchenko
            Be(3) = -0.210520d0 !Belashchenko
            Be(4) = 0.07580d0 !Belashchenko
            Be(5) = -0.449920d0 !Belashchenko
            Be(6) = -0.049920d0 !Belashchenko
            Be(7) = 0.006520d0 !Belashchenko/Fraile

            Ce(1) = 0.0326d0 !Belashchenko
            Ce(2) = 1.700d0 !Belashchenko
            Ce(3) =  -1.020d0 !Belashchenko
            Ce(4) =  1.750d0 !Belashchenko
            Ce(5) =  -1.000d0!Belashchenko
            Ce(6) = 11.0d0 !Belashchenko
            Ce(7) =  0.000d0!Belashchenko
        elseif (AWAD) then
            print*, "Awad EAM"
            me = 0d0 !Belashchenko/Fraile
            PSI(0) = 1d0 !Belashchenko/Fraile
            PSI(1) = 0.900d0 !Belashchenko/Fraile
            PSI(2) = 0.840d0 !Belashchenko
            PSI(3) = 0.700d0 !Belashchenko/Fraile
            PSI(4) = 0.550d0 !Belashchenko
            PSI(5) = 0.350d0 !Belashchenko
            PSI(6) = 1.100d0 !Belashchenko/Fraile

            Ae(1) = -1.168d0 !Awad
            Ae(2) = -1.166700d0 !Awad
            Ae(3) = -1.159848d0 !Awad
            Ae(4) = -1.136608d0 !Awad
            Ae(5) = -1.065193d0 !Awad
            Ae(6) = -0.863073d0 !Awad
            Ae(7) = -1.166700d0 !Belashchenko

            Be(1) = 0d0 !Belashchenko/Fraile
            Be(2) = -0.026000d0 !Awad
            Be(3) = -0.202400d0 !Awad
            Be(4) = -0.129600d0 !Awad
            Be(5) = -0.82600d0 !Awad
            Be(6) = -1.198600d0 !Awad
            Be(7) = +0.026000d0 !Awad

            Ce(1) = 0.13d0 !Awad
            Ce(2) = 1.47d0 !Awad
            Ce(3) = -0.26d0 !Awad
            Ce(4) = 2.31d0 !Awad
            Ce(5) = 0.94d0!Awad
            Ce(6) = 2.01d0 !Awad
            Ce(7) = 0.000d0!Awad
        endif

        ! Belashchenko EAM Pb
        if (Belashchenko_Pb) then
            rp(1) = 2.60d0 !Angs
            rp(2) = 4.60d0 !Angs
            rp(3) = 7.60d0 !Angs
            rp(4) = 9.01d0 !Angs
            ! open(97, file = "/users/edgar/BUBBLES_NEW/TAULES/Aim.Pb.table")
            open(97, file=trim(parent)//"/TAULES/Aim.Pb.table")
            do m = 0, 8
                read(97,*) aim, aaa(1,m), aaa(2,m), aaa(3,m)
            enddo
            close(97)
        elseif (Zhou_01) then
            call SET_PARAMETRES_ZHOU("Pb", "01", reZ, feZ, rhoeZ, rhosZ, alphaZ, betaZ, AZ, BZ,&
                                    kappaZ, lambdaZ, FnZ, F0Z, etaZ, FFeZ)
        elseif (Zhou_04) then
            call SET_PARAMETRES_ZHOU("Pb", "04", reZ, feZ, rhoeZ, rhosZ, alphaZ, betaZ,&
                                    AZ, BZ, kappaZ, lambdaZ, FnZ, F0Z, etaZ, FFeZ)
        endif
        ! Toennies-Tang-Sheng
        ! open(17,file="/users/edgar/BUBBLES_NEW/TAULES/TTS.table")
        open(17,file=trim(parent)//"/TAULES/TTS.table")
        read(17,*)
        read(17,*)
        read(17,*) VAR, a1(1), a1(2), a1(3)
        read(17,*) VAR, a2(1), a2(2), a2(3)
        read(17,*) VAR, a3(1), a3(2), a3(3)
        read(17,*) VAR, alpha(1), alpha(2), alpha(3)
        read(17,*) VAR, lambda(1), lambda(2), lambda(3)
        read(17,*) VAR, gamma(1), gamma(2), gamma(3)
        read(17,*) VAR, b(1), b(2), b(3)
        read(17,*) VAR, A(1), A(2), A(3)
        read(17,*) VAR, C(1,6), C(2,6), C(3,6)
        read(17,*) VAR, C(1,8), C(2,8), C(3,8)
        read(17,*) VAR, C(1,10), C(2,10), C(3,10)
        read(17,*) VAR, C(1,12), C(2,12), C(3,12)
        read(17,*) VAR, C(1,14), C(2,14), C(3,14)
        read(17,*) VAR, C(1,16), C(2,16), C(3,16)
        close(17)
        
        Re(2) = 11.47d0
        De(2) = 7.36d-6
        Re(3) = 5.608d0
        De(3) = 3.482d-5
        
        
        
    END SUBROUTINE PARAMETRES

    SUBROUTINE SETFF()
        IMPLICIT NONE

        FULLLJ = .False.
        EAMLJ =  .False.
        EAMTTSS = .False.
        Awad = .False.
        Belashchenko_Li = .False.
        Belashchenko_Pb = .False.
        Zhou_01 = .False.
        Zhou_04 = .False.
        YUKAWA = .False.
        if (force_field .eq. 0) then
            FULLLJ = .True.
        elseif (force_field .eq. 1) then
            EAMLJ = .True.
            EAM = .True.    
        elseif (force_field .eq. 2) then
            EAMTTSS = .True.
            EAM = .True.
            ! BELASHCHENKO = .True.
            Belashchenko_Li = .True.
            Belashchenko_Pb = .True.
        elseif (force_field .eq. 3) then
            EAMTTSS = .True.
            EAM = .True.
            ! BELASHCHENKO = .False.
            Belashchenko_Pb = .True.
            AWAD = .True.
        elseif (force_field .eq. 4) then
            EAMTTSS = .True.
            EAM = .True.
            ! BELASHCHENKO = .True.
            Belashchenko_Li = .True.
            Belashchenko_Pb = .True.
            AWAD = .True.
        elseif ((force_field .eq. 5) .or. (force_field .eq. 6)) then
            EAMTTSS = .True.
            EAM = .True.
            ! BELASHCHENKO = .False.
            AWAD = .True.
            if (force_field .eq. 5) then
                Zhou_01 = .True.
            elseif (force_field .eq. 6) then
                Zhou_04 = .True.
            endif


        else
            print*, "ERROR: Camp de forces no vàlid!"
            stop
        endif

    END SUBROUTINE SETFF

    SUBROUTINE WRITE_POTS(SI)
        IMPLICIT NONE
        LOGICAL SI
        REAL*8 F, dF

        IF (SI) THEN
            open(101, file="potencials.dat")
            do p = 1, 6
                write(101,*) "# ", p, "--->", parella(p)
                do k = 100, 800
                    write(101,*) dble(k)*1d-2, POTENCIAL(parella(p),dble(k)*1d-2), DERIVADA(parella(p),dble(k)*1d-2)
                enddo
                write(101,*)
                write(101,*)
            enddo
            close(101)
            do k = 0, 1000
                call embedding_Li(dble(k)*2d-3, F, dF)
            enddo
            do k = 0, 1000
                call embedding_Pb(dble(k)*2d-3, F, dF)
            enddo
        ELSE
            return
        ENDIF
    END SUBROUTINE WRITE_POTS

    REAL*8 FUNCTION epsilon(i,j)
        IMPLICIT NONE
        INTEGER, intent(in) :: i, j
        INTEGER a, b

        A = min(particle(i)%tipus,particle(j)%tipus)
        B = max(particle(i)%tipus,particle(j)%tipus)

        if ((A.eq.1) .and. (B.eq.1)) then
            epsilon = eps_LiLi
        elseif ((A.eq.1) .and. (B.eq.2)) then
            epsilon = eps_LiHe
        elseif ((A.eq.1) .and. (B.eq.3)) then
            epsilon = eps_LiPb
        elseif ((A.eq.2) .and. (B.eq.2)) then
            epsilon = eps_HeHe
        elseif ((A.eq.2) .and. (B.eq.3)) then
            epsilon = eps_PbHe
        elseif ((A.eq.3) .and. (B.eq.3)) then
            epsilon = eps_PbPb
        else
            print*, A, B, "should be either 1 or 2", i, j
            stop
        endif
        return
    END FUNCTION epsilon
    REAL*8 FUNCTION sigma(i,j)
        IMPLICIT NONE
        INTEGER, intent(in) :: i, j
        INTEGER a, b

        A = min(particle(i)%tipus,particle(j)%tipus)
        B = max(particle(i)%tipus,particle(j)%tipus)

        if ((A.eq.1) .and. (B.eq.1)) then
            sigma = sig_LiLi
        elseif ((A.eq.1) .and. (B.eq.2)) then
            sigma = sig_LiHe
        elseif ((A.eq.1) .and. (B.eq.3)) then
            sigma = sig_LiPb
        elseif ((A.eq.2) .and. (B.eq.2)) then
            sigma = sig_HeHe
        elseif ((A.eq.2) .and. (B.eq.3)) then
            sigma = sig_PbHe
        elseif ((A.eq.3) .and. (B.eq.3)) then
            sigma = sig_PbPb
        endif
        return
    END FUNCTION sigma
    CHARACTER*4 FUNCTION PAIR(i,j)
        IMPLICIT NONE
        INTEGER, intent(in) :: i, j
        INTEGER a, b

        A = min(particle(i)%tipus,particle(j)%tipus)
        B = max(particle(i)%tipus,particle(j)%tipus)

        if ((A.eq.1) .and. (B.eq.1)) then
            PAIR = "LiLi"
        elseif ((A.eq.1) .and. (B.eq.2)) then
            PAIR = "LiHe"
        elseif ((A.eq.1) .and. (B.eq.3)) then
            PAIR = "LiPb"
        elseif ((A.eq.2) .and. (B.eq.2)) then
            PAIR = "HeHe"
        elseif ((A.eq.2) .and. (B.eq.3)) then
            PAIR = "PbHe"
        elseif ((A.eq.3) .and. (B.eq.3)) then
            PAIR = "PbPb"
        endif
        return
    END FUNCTION PAIR


    REAL*8 FUNCTION POTENCIAL(parella, r)
        IMPLICIT NONE
        CHARACTER*4, INTENT(IN) :: parella
        REAL*8, INTENT(IN) :: r
        REAL*8 U_LJ, DUMMY
        REAL*8 pot_eam_Li, pot_eam_Pb, pot_eam_LiPb, pot_aeam_LiPb, pot_zhou_Pb
        REAL*8 pot_aeam_li, der_pot_aeam_li
        REAL*8 UTTS_long, UTTS_short, USladek
        REAL*8 X, CORRECCIO_YUKAWA

        if (Yukawa) X = dble(nint(dble(NLi) / (dble(NLi) + dble(NPb))*10d0)*10d0)

        if (parella .eq. "LiLi") then
            if (FULLLJ) then
                POTENCIAL = eps_LiLi*U_LJ(r/sig_LiLi)
            elseif (EAMLJ .or. EAMTTSS) then
                if (Belashchenko_Li) then
                    POTENCIAL = pot_EAM_Li(r)
                    if (Yukawa) POTENCIAL = POTENCIAL + CORRECCIO_YUKAWA(r, "potencial", 1, 1, X) 
                elseif (Awad) then
                    POTENCIAL = pot_AEAM_Li(r)
                endif
            else
                print*, "PROBLEMA AMB EL POTENCIAL Li-Li..."
                stop
            endif
        elseif (parella .eq. "PbPb") then
            if (FULLLJ) then
                POTENCIAL = eps_PbPb*U_LJ(r/sig_PbPb)
            elseif (EAMLJ .or. EAMTTSS) then
                POTENCIAL = pot_EAM_Pb(r)
                if (Belashchenko_Pb .and. Yukawa) POTENCIAL = POTENCIAL + CORRECCIO_YUKAWA(r, "potencial", 2, 2, X)
            else
                print*, "PROBLEMA AMB EL POTENCIAL Pb-Pb..."
                stop
            endif
        elseif (parella .eq. "HeHe") then
            if (FULLLJ .or. EAMLJ) then
                POTENCIAL = eps_HeHe*U_LJ(r/sig_HeHe)
            elseif (EAMTTSS) then
                POTENCIAL = ((4d0/Re(3))*UTTS_short(r/RBohr/Re(3),a1(3),a2(3),a3(3),alpha(3)) + &
                De(3)*UTTS_long(r/RBohr/Re(3),alpha(3),lambda(3),gamma(3),b(3),A(3),c(3,:)))*Eh
            else
                print*, "PROBLEMA AMB EL POTENCIAL He-He..."
                stop
            endif
        elseif (parella .eq. "LiHe") then
            if (FULLLJ .or. EAMLJ) then
                POTENCIAL = eps_LiHe*U_LJ(r/sig_LiHe)
            elseif (EAMTTSS) then
                POTENCIAL = ((6d0/Re(2))*UTTS_short(r/RBohr/Re(2),a1(2),a2(2),a3(2),alpha(2)) + &
                De(2)*UTTS_long(r/RBohr/Re(2),alpha(2),lambda(2),gamma(2),b(2),A(2),c(2,:)))*Eh
            else
                print*, "PROBLEMA AMB EL POTENCIAL Li-He..."
                stop
            endif
        elseif (parella .eq. "LiPb") then
            if (FULLLJ .or. EAMLJ) then
                POTENCIAL = eps_LiPb*U_LJ(r/sig_LiPb)
            elseif (EAMTTSS) then
                if (Belashchenko_Li .and. Belashchenko_Pb) then
                    POTENCIAL = eps_84*pot_EAM_LiPb(r/sig_84)
                    if (Yukawa) POTENCIAL = POTENCIAL + CORRECCIO_YUKAWA(r, "potencial", 1, 2, X)
                elseif (Awad) then
                    POTENCIAL = pot_AEAM_LiPb(r)
                elseif (Zhou_01 .or. Zhou_04) then
                    POTENCIAL = pot_Zhou_Pb(r)
                endif
            else
                print*, "PROBLEMA AMB EL POTENCIAL Li-Pb..."
                stop
            endif
        elseif (parella .eq. "PbHe") then
            if (FULLLJ .or. EAMLJ) then
                POTENCIAL = eps_PbHe*U_LJ(r/sig_PbHe)
            elseif (EAMTTSS) then
                POTENCIAL = USladek(r)
            else
                print*, "PROBLEMA AMB EL POTENCIAL Pb-He..."
                stop
            endif
        else
            print*, "PASSA ALGUNA COSA DOLENTA"
            print*, parella
            print*, "................"
            stop
        endif

        return
    END FUNCTION POTENCIAL


    REAL*8 FUNCTION DERIVADA(parella, r)
        IMPLICIT NONE
        CHARACTER*4, INTENT(IN) :: parella
        REAL*8, INTENT(IN) :: r
        REAL*8 dU_LJ, DUMMY
        REAL*8 der_pot_eam_Li, der_pot_eam_Pb, der_pot_eam_LiPb, der_pot_aeam_li, der_pot_AEAM_LiPb,&
        der_pot_ZHOU_Pb
        REAL*8 dUTTS_long, dUTTS_short, dUSladek
        REAL*8 X, CORRECCIO_YUKAWA

        if (Yukawa) X = dble(int(dble(NLi) / (dble(NLi) + dble(NPb))*10d0)*10d0)

        if (parella .eq. "LiLi") then
            if (FULLLJ) then
                DERIVADA = (eps_LiLi/sig_LiLi)*dU_LJ(r/sig_LiLi)
            elseif (EAMLJ .or. EAMTTSS) then
                if (Belashchenko_Li) then
                    DERIVADA = der_pot_EAM_Li(r)
                    if (Yukawa) DERIVADA = DERIVADA + CORRECCIO_YUKAWA(r, "derivada_", 1, 1, X)
                elseif (Awad) then
                    DERIVADA = der_pot_AEAM_Li(r)
                endif
            else
                print*, "PROBLEMA AMB EL POTENCIAL Li-Li..."
                stop
            endif
        elseif (parella .eq. "PbPb") then
            if (FULLLJ) then
                DERIVADA = (eps_PbPb/sig_PbPb)*dU_LJ(r/sig_PbPb)
            elseif (EAMLJ .or. EAMTTSS) then
                DERIVADA = der_pot_EAM_Pb(r)
                if (Belashchenko_Pb .and. Yukawa) DERIVADA = DERIVADA + CORRECCIO_YUKAWA(r, "derivada_", 2, 2, X)
            else
                print*, "PROBLEMA AMB EL POTENCIAL Pb-Pb..."
                stop
            endif
        elseif (parella .eq. "HeHe") then
            if (FULLLJ .or. EAMLJ) then
                DERIVADA = (eps_HeHe/sig_HeHe)*dU_LJ(r/sig_HeHe)
            elseif (EAMTTSS) then
                DERIVADA = ((4d0/Re(3))*dUTTS_short(r/RBohr/Re(3),a1(3),a2(3),a3(3),alpha(3)) + &
                De(3)*dUTTS_long(r/RBohr/Re(3),alpha(3),lambda(3),gamma(3),b(3),A(3),c(3,:)))*Eh/RBohr/Re(3) !a.u. (Hartree / radi Bohr / Re) --> eV/Å
            else
                print*, "PROBLEMA AMB EL POTENCIAL He-He..."
                stop
            endif
        elseif (parella .eq. "LiHe") then
            if (FULLLJ .or. EAMLJ) then
                DERIVADA = (eps_LiHe/sig_LiHe)*dU_LJ(r/sig_LiHe)
            elseif (EAMTTSS) then
                DERIVADA = ((6d0/Re(2))*dUTTS_short(r/RBohr/Re(2),a1(2),a2(2),a3(2),alpha(2)) + &
                De(2)*dUTTS_long(r/RBohr/Re(2),alpha(2),lambda(2),gamma(2),b(2),A(2),c(2,:)))*Eh/RBohr/Re(2) !a.u. (Hartree / radi Bohr / Re) --> eV/Å
            else
                print*, "PROBLEMA AMB EL POTENCIAL Li-He..."
                stop
            endif
        elseif (parella .eq. "LiPb") then
            if (FULLLJ .or. EAMLJ) then
                DERIVADA = (eps_LiPb/sig_LiPb)*dU_LJ(r/sig_LiPb)
            elseif (EAMTTSS) then
                if (Belashchenko_Li .and. Belashchenko_Pb) then
                    DERIVADA = (eps_84/sig_84)*der_pot_EAM_LiPb(r/sig_84)
                    if (Yukawa) DERIVADA = DERIVADA + CORRECCIO_YUKAWA(r, "derivada_", 1, 2, X)
                elseif (Awad) then
                    DERIVADA = der_pot_AEAM_LiPb(r)
                elseif (Zhou_01 .or. Zhou_04) then
                    DERIVADA = der_pot_ZHOU_Pb(r)
                endif
            else
                print*, "PROBLEMA AMB EL POTENCIAL Li-Pb..."
                stop
            endif
        elseif (parella .eq. "PbHe") then
            if (FULLLJ .or. EAMLJ) then
                DERIVADA = (eps_PbHe/sig_PbHe)*dU_LJ(r/sig_PbHe)
            elseif (EAMTTSS) then
                DERIVADA = dUSladek(r)
            else
                print*, "PROBLEMA AMB EL POTENCIAL Pb-He..."
                stop
            endif
        else
            print*, "PASSA ALGUNA COSA DOLENTA"
            stop
        endif
    END FUNCTION DERIVADA

	SUBROUTINE CLUSTERS(rcut,lcs,Heincluster)
        IMPLICIT NONE
        !type(list), dimension(NLi+NPb+1:N) :: atom
        type(list), dimension(1:NHe) :: cluster, atom !Com a cas extrem tindrem NHe clusters (1 àtom = 1 cluster).
        integer, dimension(:), allocatable :: csize
        integer nclusters, ic, jc, kk, lcs, lc, aux, ni, nj
        real*8, intent(in) :: rcut
        integer, intent(out), dimension(:), allocatable :: Heincluster
        
            
        nclusters = 0
        
        do i = 1, NHe!i = NLi+NPb+1, NLi+NPb+NHe
           k = 0
           l = 0
           ic = 0
           do l = 1, nclusters
              if (any(cluster(l)%nbr==i)) then !Està ja l'àtom i en un cluster?
                 ic = l
                 exit
              endif
           enddo
           if (ic.eq.0) then
              nclusters = nclusters + 1 !Quan no hem trobat i a cap cluster, en creem un de nou.
              call add_integer(cluster(nclusters)%nbr,i) !Afegim l'àtom i al nou cluster.
              ic = nclusters !Etiquetem el nou cluster que conté l'àtom i.
           endif
           do j = i+1, NHe !j = i+1, NLi+NPb+NHe
              dist = distancia(He(i)%posicio, He(j)%posicio, .True.)
              if (dist .le. rcut) then
                 call add_integer(atom(i)%nbr,j) !Afegim j al veïnat de i.
                 call add_integer(atom(j)%nbr,i) !Afegim i al veïnat de j.
              
                 k = 0
                 l = 0
                 jc = 0
                 do l = 1, nclusters
                    if (any(cluster(l)%nbr==j)) then !Està ja l'àtom j en un cluster?
                       jc = l
                       exit
                    endif 
                 enddo
                 
                 
                 if (jc .eq. 0) then
                    call add_integer(cluster(ic)%nbr,j) !En cas que j no estigui en cap cluster s'hi afegirà al mateix que i.
                    
                 elseif (ic .ne. jc) then !Degut al criteri de cutoff considerem que i & j pertanyen al mateix cluster. Ens quedem amb l'etiqueta més petita.
                    ni = 0
                    if (allocated(cluster(ic)%nbr)) ni = size(cluster(ic)%nbr)
                    nj = 0
                    if (allocated(cluster(jc)%nbr)) nj = size(cluster(jc)%nbr)
                    if (ic .lt. jc) then
                       do kk = 1, nj
                          k = cluster(jc)%nbr(kk)
                          call add_integer(cluster(ic)%nbr,k) ! Afegim tant j com els seus veïns k, al mateix cluster que i.
                       enddo
                    elseif (ic .gt. jc) then
                       do kk = 1, ni
                          k = cluster(ic)%nbr(kk)
                          call add_integer(cluster(jc)%nbr,k) ! Afegim tant i com els seus veïns k, al mateix cluster que j.
                       enddo
                    endif
                    do l = max(ic,jc), nclusters-1
                       deallocate(cluster(l)%nbr)
                       allocate(cluster(l)%nbr(size(cluster(l+1)%nbr)))
                       cluster(l)%nbr = cluster(l+1)%nbr
                    enddo
                    deallocate(cluster(nclusters)%nbr)
                    nclusters = nclusters - 1
                    ic = min(ic,jc)
                 endif
              endif
           enddo
        enddo
        
        if (PRINTALL) then
            write(*,*)
            write(*,"(A50)") "CLUSTER ANALYSIS"
            write(*,"(A50)") "--------------------------------------------------"
            write(*,"(A4, A32, A1, I5)") "   >", "# clusters", ":", nclusters
        endif
        aux = 0
        do l = 1, nhe
           if (allocated(cluster(l)%nbr)) aux = aux + 1 
        enddo 
        allocate(csize(nclusters))
        lcs = 0
        lc = 0
        aux = 0
        do l = 1, nclusters
           csize(l) = size(cluster(l)%nbr)
           aux = aux + csize(l)
           if (csize(l) .gt. lcs) then
              lcs = csize(l) !Mida del cluster més gran
              lc = l !Etiqueta del cluster més gran.
           endif
        enddo
        if (PRINTALL) then
            write(*,"(A4, A32, A1, I5, I5)") "   >", "Largest cluster size", ":", maxval(csize), lcs
            write(*,"(A4, A32, A1, I5, I5)") "   >", "Recompte atoms d'He (total)", ":", aux
            write(*,"(A4, A32, A1)") "   >", "Atoms al cluster principal", ":"
            write(*,"(A24, A12)") "", "------------------------------------------------"
            write(*,"(A1, A24, A1, A2, A2, A6)") "", "", "", "#", "", "Id"
            write(*,"(A24, A12)") "", "------------------------------------------------"
            do l = 1, size(cluster(lc)%nbr)
                write(*,"(A25, I3, A2, I6)") "", l, ") ", cluster(lc)%nbr(l)
            enddo
            write(*,"(A50)") "--------------------------------------------------"
            write(*,*)
        endif
        
        if (allocated(Heincluster)) deallocate(Heincluster)
        allocate(Heincluster(lcs))
        do j = 1, lcs
           i = cluster(lc)%nbr(j)
           Heincluster(j) = i
        enddo
            
        return
    END SUBROUTINE CLUSTERS
    SUBROUTINE VIRIAL_PARALLEL(i, j, sum_Fij_rij, F) !Calcul pressió...................
        IMPLICIT NONE
        integer, intent(in) :: i, j
        real*8, intent(inout) :: sum_Fij_rij, F(N,3)
        real*8, dimension(3) :: ri, rj
        REAL*8 rij, dudr, psi, dpsi, rhoi, rhoj, Fi, Fj, dFi, dFj, dij(3), Fij(3)
        character*4 parella

        parella = PAIR(i,j)
        ri(:) = pos(i,:) - poscm(:)
        rj(:) = pos(j,:) - poscm(:)
        dij(:) = ri(:) - rj(:)
        if (SETPBC) then
            call PBC(dij)
            call PBC(ri)
            call PBC(rj)
        endif
        rij = distancia(ri,rj,.True.)

        if (EAM) then
            if ((parella.eq."LiLi") .or. (parella.eq."PbPb") .or. (parella.eq."LiPb")) then
                call electro(rij, psi, dpsi, p1(tj), p2) ! COMPTE!!!! El terme p1 és creuat!
                rhoi = edens(i)
                call embedding(particle(i)%nom, rhoi, Fi, dFi)
                call electro(rij, psi, dpsi, p1(ti), p2)
                rhoj = edens(j)
                call embedding(particle(j)%nom, rhoj, Fj, dFj)
                dudr = ((dFi + dFj) * dpsi + DERIVADA(parella,rij))
            else
                dudr = DERIVADA(parella,rij)
            endif
        endif

        Fij(:) = -dudr * dij(:)/dsqrt(dot_product(dij,dij))


        sum_Fij_rij = sum_Fij_rij - dudr * rij
        
        F(i,:) = F(i,:) - Fij(:)
        F(j,:) = F(j,:) + Fij(:)

        return
    END SUBROUTINE VIRIAL_PARALLEL

    SUBROUTINE P_CONF_SERIAL(i, j, pv)
        IMPLICIT NONE
        integer, intent(in) :: i, j
        real*8, dimension(nr0:nr1), intent(inout) :: pv
        real*8, dimension(3) :: ri, rj
        real *8 rmin, rmax, lmin
        real*8 ri2, rj2, rij, rij2, dij(3), rrij, dUij
        integer k, kmin, kmax, numtalls
        real*8 lambda_p, lambda_n
        real*8 dU_LJ, dudr
        REAL*8 psi, dpsi, rhoi, rhoj, Fi, Fj, dFi, dFj
        character*4 parella
        REAL*8 dUSladek, Usladek
        REAL*8, dimension(nr0:nr1) :: suma_r

        parella = PAIR(i,j)

        ri = pos(i,:) - poscm(:)
        rj = pos(j,:) - poscm(:)
        if (SETPBC) then
            call PBC(ri)
            call PBC(rj)
        endif
        
        ri2 = dot_product(ri,ri)
        rj2 = dot_product(rj,rj)
        dij(:) = ri(:) - rj(:)
        if (SETPBC) call PBC(dij)
        rij2 = dot_product(dij,dij)
        rij = distancia(ri,rj,.True.)
        dudr = 0d0
        if (force_field .eq. 0)  dudr = (epsilon(i,j)/sigma(i,j))*dU_LJ(rij/sigma(i,j))
        if (EAM) then
            if (parella.eq."LiLi") then
                call electro(rij, psi, dpsi, p1(1), p2)
                rhoi = edens(i)
                call embedding_Li(rhoi,Fi,dFi)
                rhoj = edens(j)
                call embedding_Li(rhoj,Fj,dFj)
                dudr = ((dFi + dFj) * dpsi + DERIVADA(parella,rij) ) 
            elseif (parella.eq."PbPb") then
                call electro(rij, psi, dpsi, p1(3), p2)
                rhoi = edens(i)
                call embedding_Pb(rhoi,Fi,dFi)
                rhoj = edens(j)
                call embedding_Pb(rhoj,Fj,dFj)
                dudr = ((dFi + dFj) * dpsi + DERIVADA(parella,rij) )
            else
                dudr = DERIVADA(parella,rij)
            endif
        else
            dudr = DERIVADA(parella,rij)
        endif
        
        lmin = (ri2-rj2) / rij2
	
        if (dabs(lmin) .le. 1d0) then
            rmin = 0.5d0*dsqrt(2d0*(ri2+rj2)-rij2-(ri2-rj2)**2d0/rij2)
        else
            rmin = min(dsqrt(ri2),dsqrt(rj2))
        endif
        if (isnan(rmin)) then
            print*, "WARNING: rmin is NaN"
            print*, ri2, rj2, 2d0*(ri2+rj2)-rij2-(ri2-rj2)**2d0/rij2
            rmin = min(dsqrt(ri2),dsqrt(rj2))
        endif

        rmax = max(dsqrt(ri2),dsqrt(rj2))

        kmin = floor(rmin / Dr) + 1
        kmax = floor(rmax / Dr) + 1
        kmax = min(kmax, nr) !per a no sobrepassar el maxim de l'array r
        
        suma_r = 0d0

        do k = kmin, kmax
            lambda_p = lmin + dsqrt(lmin*lmin + 1d0 - 2d0*(ri2+rj2)/rij2 + 4d0*r(k)*r(k)/rij2)
            lambda_n = lmin - dsqrt(lmin*lmin + 1d0 - 2d0*(ri2+rj2)/rij2 + 4d0*r(k)*r(k)/rij2)
            
            numtalls = 0
            
            if ((dabs(lambda_p).gt.1d0).or.(dabs(lambda_n).gt.1d0)) then
                numtalls = 1
                c1 = c1 + 1
            elseif ((dabs(lambda_p).lt.1d0).and.(dabs(lambda_n).lt.1d0)) then
                numtalls = 2
                c2 = c2 + 1
            else
                numtalls = 0
                rrij = 0d0
                c0 = c0 + 1
            endif
            ctot = ctot + 1

            rrij = 0.5d0*rij2*dsqrt(lmin*lmin + 1d0 - 2d0*(ri2+rj2)/rij2 + 4d0*r(k)*r(k)/rij2)
            if (isnan(rrij).eqv..false.) then
                if (numtalls.eq.0) print*, lmin, lambda_n, lambda_p, numtalls, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", rrij !Això no hauria de passar mai, ja que numtalls=0 quan rrij es NaN...
                suma_r(k) = suma_r(k) + numtalls*rrij
            endif
        enddo
        
        pv(:) = pv(:) - suma_r(:) * dudr / rij

        return
    END SUBROUTINE P_CONF_SERIAL

    SUBROUTINE P_CONF_PARALLEL(i, j, pv, c0, c1, c2, ctot)
        IMPLICIT NONE
        integer, intent(in) :: i, j
        real*8, dimension(nr0:nr1), intent(inout) :: pv
        integer, intent(inout) :: c0, c1, c2, ctot
        real*8, dimension(3) :: ri, rj
        real *8 rmin, rmax, lmin
        real*8 ri2, rj2, rij, rij2, dij(3), rrij, dUij
        integer k, kmin, kmax, numtalls
        real*8 lambda_p, lambda_n
        real*8 dU_LJ, dudr
        REAL*8 psi, dpsi, rhoi, rhoj, Fi, Fj, dFi, dFj
        character*4 parella
        REAL*8 dUSladek, Usladek
        REAL*8, dimension(nr0:nr1) :: pvpar

        pvpar = 0d0

        parella = PAIR(i,j)

        ri = pos(i,:) - poscm(:)
        rj = pos(j,:) - poscm(:)
        if (SETPBC) then
            call PBC(ri)
            call PBC(rj)
        endif
        
        ri2 = dot_product(ri,ri)
        rj2 = dot_product(rj,rj)
        dij(:) = ri(:) - rj(:)
        if (SETPBC) call PBC(dij)
        rij2 = dot_product(dij,dij)
        rij = distancia(ri,rj,.True.)
        if (force_field .eq. 0)  dudr = (epsilon(i,j)/sigma(i,j))*dU_LJ(rij/sigma(i,j))
        if (EAM) then
            if (parella.eq."LiLi") then
                call electro(rij, psi, dpsi, p1(1), p2)
                rhoi = edens(i)
                call embedding_Li(rhoi,Fi,dFi)
                rhoj = edens(j)
                call embedding_Li(rhoj,Fj,dFj)
                dudr = ((dFi + dFj) * dpsi + DERIVADA(parella,rij) ) 
            elseif (parella.eq."PbPb") then
                call electro(rij, psi, dpsi, p1(3), p2)
                rhoi = edens(i)
                call embedding_Pb(rhoi,Fi,dFi)
                rhoj = edens(j)
                call embedding_Pb(rhoj,Fj,dFj)
                dudr = ((dFi + dFj) * dpsi + DERIVADA(parella,rij) )
            else
                dudr = DERIVADA(parella,rij)
            endif
        else
            dudr = DERIVADA(parella,rij)
        endif
        
        lmin = (ri2-rj2) / rij2
	
        if (dabs(lmin) .le. 1d0) then
            rmin = 0.5d0*dsqrt(2d0*(ri2+rj2)-rij2-(ri2-rj2)**2d0/rij2)
        else
            rmin = min(dsqrt(ri2),dsqrt(rj2))
        endif
        if (isnan(rmin)) then
            print*, "WARNING: rmin is NaN"
            print*, ri2, rj2, 2d0*(ri2+rj2)-rij2-(ri2-rj2)**2d0/rij2
            rmin = min(dsqrt(ri2),dsqrt(rj2))
        endif

        rmax = max(dsqrt(ri2),dsqrt(rj2))

        kmin = floor(rmin / Dr) + 1
        kmax = floor(rmax / Dr) + 1
        kmax = min(kmax, nr) !per a no sobrepassar el maxim de l'array r

        pvpar = 0d0
        do k = kmin, kmax
            lambda_p = lmin + dsqrt(lmin*lmin + 1d0 - 2d0*(ri2+rj2)/rij2 + 4d0*r(k)*r(k)/rij2)
            lambda_n = lmin - dsqrt(lmin*lmin + 1d0 - 2d0*(ri2+rj2)/rij2 + 4d0*r(k)*r(k)/rij2)
            
            numtalls = 0
            
            if ((dabs(lambda_p).gt.1d0).or.(dabs(lambda_n).gt.1d0)) then
                numtalls = 1
                c1 = c1 + 1
            elseif ((dabs(lambda_p).lt.1d0).and.(dabs(lambda_n).lt.1d0)) then
                numtalls = 2
                c2 = c2 + 1
            else
                numtalls = 0
                rrij = 0d0
                c0 = c0 + 1
            endif
            ctot = ctot + 1

            rrij = 0.5d0*rij2*dsqrt(lmin*lmin + 1d0 - 2d0*(ri2+rj2)/rij2 + 4d0*r(k)*r(k)/rij2)
            if (isnan(rrij).eqv..false.) then
                if (numtalls.eq.0) print*, lmin, lambda_n, lambda_p, numtalls, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", rrij !Això no hauria de passar mai, ja que numtalls=0 quan rrij es NaN...
                pvpar(k) = pvpar(k) - numtalls * rrij * dudr / rij
            endif
        enddo
        do k = kmin, kmax
            pv(k) = pv(k) + pvpar(k)
        enddo
        return
    END SUBROUTINE P_CONF_PARALLEL
    SUBROUTINE Parallel_Hello_World
        USE omp_lib
        INTEGER :: thread_id
        INTEGER :: partial_Sum, total_Sum
    
        print*, "Max threads: ", OMP_GET_MAX_THREADS()
    
        !$OMP PARALLEL PRIVATE(thread_id, partial_Sum) SHARED(total_Sum)
        
    
        thread_id = OMP_GET_THREAD_NUM()
    
        PRINT*, "Hello from process: ", thread_id
        !$OMP BARRIER
    
        DO i=0,OMP_GET_MAX_THREADS()
            IF (i == thread_id) THEN
                print*, "---"
                print*, i, OMP_GET_MAX_THREADS()
                print*, "---"
                PRINT *, "Hello again from process: ", thread_id
            END IF
            !$OMP BARRIER
        END DO
    
        partial_Sum = 0
        total_Sum = 0
        !$OMP DO
        DO i=1,1000
            partial_Sum = partial_Sum + i
        END DO
        !$OMP END DO
        print*, "Suma parcial al processador", thread_id," --->", partial_Sum
        !$OMP CRITICAL
        total_Sum = total_Sum + partial_Sum
        !$OMP END CRITICAL
        
        !$OMP END PARALLEL
    
        PRINT *, "Total Sum: ", total_Sum
    
    END SUBROUTINE Parallel_Hello_World
END PROGRAM ST
