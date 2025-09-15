MODULE CONSTANTS
    IMPLICIT NONE
    REAL*8, PARAMETER :: pi=4d0*datan(1d0)
    REAL*8, PARAMETER :: NA=6.02214076e23 !mol^-1 ---> Nombre d'Avogadro
    REAL*8, PARAMETER :: kB=8.6173324d-5 ! [kB] = eV/K ------> cnt. Boltzmann
    REAL*8, PARAMETER :: e=1.602176565d-19 ! [e] = C --------> càrrega electró 
    REAL*8, PARAMETER :: Eh=27.2113845d0 ! [Eh] = eV --------> Hartree
    REAL*8, PARAMETER :: Rbohr=0.5291772108d0  ! [Rbohr] = "Å" --> Radi de Bohr
    REAL*8, parameter :: massa_Li=6.941d0, massa_Pb=207.2d0, massa_He=4.002d0 ! g/mol
    CHARACTER*4, dimension(1:6), parameter :: parella=(/ "LiLi", "LiPb", "LiHe", "PbPb", "PbHe", "HeHe" /)
END MODULE CONSTANTS