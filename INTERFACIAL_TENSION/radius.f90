! Edgar Alvarez Galera
! Desplaçament quadràtic mitjà - sistema Li+Pb+He --- EAM (Awad) + TTS
! Data de creació: 06 de març del 2025
! Última modificació: 30 / 06 / 2025

PROGRAM RADI
    USE CONSTANTS
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
        integer function findmin(list)
            real*8, dimension(:), allocatable, intent(in) :: list
        end function findmin
        subroutine ordena(list)
            integer, dimension(:), allocatable, intent(inout) :: list
        end subroutine ordena
    END INTERFACE
    TYPE peratom
        integer tipus, id
        character*2 nom
        real*8 massa
        real*8, dimension(3) :: posicio, forca
    END TYPE peratom
    TYPE list
        integer, dimension(:), allocatable :: nbr, group
    END TYPE list
    INTEGER N, NHe, NLi, NPb, Nats, Nd, n1, n2, n3, nt(1:3)
    REAL*8 Vol, Temp, Press, Box, xl, xh, yl, yh, zl, zh, boxl(1:3), boxh(1:3)
    INTEGER ts, tsi, tsf, every, sni, snf, snapshots, timesteps, s
    REAL*8 Dr
    INTEGER nr0, nr1, force_field
    LOGICAL POTS, SETPBC, PRINTALL, PRINTSNAP
    INTEGER i, j, k, l, m, p
    REAL*8, dimension(:,:), allocatable :: pos
    REAL*8, dimension(1:3) :: poscm
    LOGICAL ALL_He_ATS
    PARAMETER(ALL_He_ATS=.False.)
    INTEGER, dimension(:), allocatable :: ATS_cluster
    REAL*8, dimension(:), allocatable :: coord
    REAL*8, dimension(3), parameter :: massa=(/ massa_Li, massa_He, massa_Pb /)
    REAL*8, dimension(:), allocatable :: masscl
    TYPE(peratom), dimension(:), allocatable :: particle, LM, He
    REAL*8 MSDR, MCR, mindisttocm, maxdisttocm, distancia, dist, COM
    REAL*8 Rd, Rd2, Rt, Rt2, Rmin, Rmin2, Rmax, Rmax2, Rnew, Rnew2, Reqm
    INTEGER dummy
    REAL*8, dimension(:), allocatable :: dens, r
    INTEGER, dimension(:), allocatable :: number
    REAL*8, dimension(:), allocatable :: num_Li, num_Pb, num_He, dens_Li, dens_Pb, dens_He
    INTEGER, dimension(1) :: ml
    REAL*8 INTEGRAL
    COMMON/CAPSA/boxl,boxh

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

    Dr = 0.1d0 !Angstroms
    nr0 = 1
    nr1 = 350


    allocate(pos(1:N, 1:3))
    allocate(LM(1:NLi+NPb), He(1:NHe), particle(1:N))
    allocate(r(nr0:nr1), dens(nr0:nr1), number(nr0:nr1), &
    num_Li(nr0:nr1), num_He(nr0:nr1), num_Pb(nr0:nr1), &
    dens_Li(nr0:nr1), dens_He(nr0:nr1), dens_Pb(nr0:nr1))

    do k = nr0, nr1
        r(k) = Dr*dble(k) - Dr/2d0
    enddo

    open(15, FILE="xyz.out")
    open(16, FILE="simulation_box.out")
    read(16,*)

    Rd = 0d0
    Rt = 0d0
    Rmin = 0d0
    Rmax = 0d0
    Rnew = 0d0
    Rd2 = 0d0
    Rt2 = 0d0
    Rmin2 = 0d0
    Rmax2 = 0d0
    Rnew2 = 0d0
    dummy = 0
    dens = 0d0
    number = 0
    num_He = 0d0
    num_Li = 0d0
    num_Pb = 0d0

    do s = sni, snf
        read(16,*) ts, xl, xh, yl, yh, zl, zh
        Vol = (xh-xl)*(yh-yl)*(zh-zl)
        Box = (Vol)**(1d0/3d0)
        boxl = (/ xl, yl, zl /)
        boxh = (/ xh, yh, zh /)

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
        nt = 0
        do i = 1, N
            read(15,*) k, pos(i,:)
            nt(k) = nt(k) + 1 
            if (SETPBC) call PBC(pos(i,:))
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
        enddo

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

        MSDR = 0d0
        MCR = 0d0
        minDISTTOCM = 1e8
        maxDISTTOCM = 0d0
        Rnew = 0d0
        do k = 1, Nd
            i = ATS_cluster(k)
            MSDR = MSDR + distancia(He(i)%posicio,poscm(:), .True.)**2d0
            MCR = MCR + distancia(He(i)%posicio,poscm(:), .True.)**3d0
            if (distancia(He(i)%posicio, poscm(:), .True.) .gt. maxDISTTOCM) &
            maxDISTTOCM = distancia(He(i)%posicio, poscm(:), .True.)
        enddo
        do k = 1, NLi+NPb
            i = LM(k)%id
            if (distancia(LM(i)%posicio, poscm(:), .True.) .lt. minDISTTOCM) &
            minDISTTOCM = distancia(LM(i)%posicio, poscm(:), .True.)
        enddo

        do i = 1, N
            k = floor(distancia(pos(i,:),poscm(:),.True.)/Dr) + 1
            if ((k.ge.nr0) .and. (k.le.nr1)) then
                number(k) = number(k) + 1
                if (particle(i)%nom .eq. "Li") then
                    num_Li(k) = num_Li(k) + 1d0
                elseif (particle(i)%nom .eq. "Pb") then
                    num_Pb(k) = num_Pb(k) + 1d0
                elseif (particle(i)%nom .eq. "He") then
                    num_He(k) = num_He(k) + 1d0
                else
                    print*, "ERROR!!!! Particle name is not either 'Li', 'Pb' or 'He'..."
                    stop
                endif
            endif
        enddo
        
        
        MSDR = MSDR / dble(Nd)
        Rd = Rd + dsqrt(MSDR)
        Rd2 = Rd2 + dsqrt(MSDR)*dsqrt(MSDR)
        MCR = MCR / dble(Nd)
        Rt = Rt + (MCR)**(1d0/3d0)
        Rt2 = Rt2 + (MCR)**(2d0/3d0)

        Rmax = Rmax + maxdisttocm
        Rmax2 = Rmax2 + maxdisttocm*maxdisttocm

        Rmin = Rmin + mindisttocm
        Rmin2 = Rmin2 + mindisttocm*mindisttocm

        dummy = dummy + 1        

        if (PRINTALL) then
            write(*,"(A20, F11.8, A4)") "MEAN SQUARED RADIUS:", dsqrt(MSDR), "\\AA"
            write(*,"(A20, F11.8, A4)") "MIN DIST TO CM:", minDISTTOCM, "\\AA"
            write(*,"(A20, F11.8, A4)") "MAX DIST TO CM:", maxDISTTOCM, "\\AA"
            write(*,"(A75)") "---------------------------------------------------------------------------"
            write(*,*)
        endif
    enddo

    Rd = Rd / dble(snf-sni+1)
    Rt = Rt / dble(snf-sni+1)
    Rmin = Rmin / dble(snf-sni+1)
    Rmax = Rmax / dble(snf-sni+1)
    Rd2 = Rd2 / dble(snf-sni+1)
    Rt2 = Rt2 / dble(snf-sni+1)
    Rmin2 = Rmin2 / dble(snf-sni+1)
    Rmax2 = Rmax2 / dble(snf-sni+1)
    Rnew = Rnew / dble(snf-sni+1)
    Rnew2 = Rnew2 / dble(snf-sni+1)

    do k = nr0, nr1
        dens(k) = dble(number(k)) / (4d0 * pi * Dr * (r(k)*r(k) + Dr*Dr/12d0)) / dble(snf-sni+1)
        dens_Li(k) = num_Li(k)  / (4d0 * pi * Dr * (r(k)*r(k) + Dr*Dr/12d0)) / dble(snf-sni+1)
        dens_Pb(k) = num_Pb(k)  / (4d0 * pi * Dr * (r(k)*r(k) + Dr*Dr/12d0)) / dble(snf-sni+1)
        dens_He(k) = num_He(k)  / (4d0 * pi * Dr * (r(k)*r(k) + Dr*Dr/12d0)) / dble(snf-sni+1)
    enddo
    
    m = nr1
    do k = nr0+50, nr1-50
        if (dens(k) .lt. dens(m)) m = k
        write(97,*) r(k), dens(k), dens_Li(k), dens_Pb(k), dens_He(k)
    enddo
    Rnew = r(m)

    Reqm = (abs(INTEGRAL(r,dens,50,300,3d0) / INTEGRAL(r,dens,50,300,0d0)))**(1d0/3d0) ! Å

    write(*, "(A14, F6.3, A5, F6.3, A5)") "AVERAGE MSR: ", Rd, " \pm", dsqrt((Rd2 - Rd*Rd)/dble(snf-sni)), " \AA"
    write(*, "(A14, F6.3, A5, F6.3, A5)") "AVERAGE MCR: ", Rt, " \pm", dsqrt((Rt2 - Rt*Rt)/dble(snf-sni)), " \AA"
    write(*, "(A14, F6.3, A5, F6.3, A5)") "LM RADIUS: ", Rmin, " \pm", dsqrt((Rmin2 - Rmin*Rmin)/dble(snf-sni)), " \AA"
    write(*, "(A14, F6.3, A5, F6.3, A5)") "He RADIUS: ", Rmax, " \pm", dsqrt((Rmax2 - Rmax*Rmax)/dble(snf-sni)), " \AA"
    write(*, "(A14, F6.3, A5, F6.3, A5)") "NEW RADIUS: ", Rnew, " \pm", dsqrt((Rnew2 - Rnew*Rnew)/dble(snf-sni)), " \AA"
    write(*, "(A14, F6.3, A5, F6.3, A5)") "EQUIM RADIUS: ", Reqm, " \pm", dsqrt((Rnew2 - Rnew*Rnew)/dble(snf-sni)), " \AA"

    CONTAINS
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
                !dist = distancia(pos_He(i,:),pos_He(j,:),.True.)
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
END PROGRAM RADI

INTEGER FUNCTION findmin(list)
    implicit none
    real*8, dimension(:), allocatable, intent(in) :: list
    integer :: element
    integer i

    findmin = 0
    element = minval(list)

    if (allocated(list)) then
        do i = 1, size(list)
            if (list(i) .eq. element) then
            findmin = i
            return
            endif
        enddo
    endif

END FUNCTION findmin