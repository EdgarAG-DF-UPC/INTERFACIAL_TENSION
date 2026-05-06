import numpy as np

def Vcav(r,lamb):
    A = 12.483 # eV
    B = 1.
    return A*np.exp(-r/B+lamb)

def U_LJ(r): # Dimensionless units
    return 4*(r**(-12) - r**(-6))

def ULiHe(r):
    epsilon = 0.00014100 #eV
    sigma = 5.3565 # Angstrom
    return epsilon*U_LJ(r/sigma)

def UPbHe(r):
    epsilon = 0.01057000 #eV
    sigma = 3.0667 # Angstrom
    return epsilon*U_LJ(r/sigma)