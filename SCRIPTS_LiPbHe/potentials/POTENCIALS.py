import numpy as np

def Vcav(r,lamb):
    A = 12.483 # eV
    B = 1.
    return A*np.exp(-r/B+lamb)

def U_alk_He(r, alk):
    r = r / 0.5291772
    if alk=="Li":
        R_e = 11.47
        D_e = 7.36e-6
        Z_A = 2.
        Z_B = 3.
    elif alk=="Na":
        R_e = 11.85
        D_e = 6.96e-6
        Z_A = 2.
        Z_B = 11.
    elif alk=="K":
        R_e = 13.50
        D_e = 5e-6
        Z_A = 2.
        Z_B = 19.
    elif alk=="Rb":
        R_e = 13.86
        D_e = 4.76e-6
        Z_A = 2.
        Z_B = 37.
    elif alk=="Cs":
        R_e = 14.89
        D_e = 3.36e-6
        Z_A = 2.
        Z_B = 55.
    
    a1 = 11.49538
    a2 = -25.99090
    a3 = 13.49545
    alpha = 25.06466
    lambd = 14.57866
    gamma = 2.19783
    b = 12.38083
    A = 2.31539e6
    c = np.zeros(20)
    c[6] = 1.43308
    c[8] = 0.47909
    c[10] = 0.22501
    c[12] = 0.14847
    c[14] = 0.13762
    c[16] = 0.17921
    
    def Ushort(x):
        return (1./x) * (1. + a1*x + a2*x**2 + a3*x**3)*np.exp(-alpha*x)
    def Ulong(x):
        sumn = 0.
        for n in range(3,9):
            sumk = 1.
            factk = 1.
            for k in range(1,2*n+1):
                factk = factk * k
                sumk = sumk + (b*x)**k / factk
            sumn = sumn + (1. - np.exp(-b*x)*sumk)*c[2*n]/(x**(2*n))
        U_vdw = A * x**gamma * np.exp(-lambd*x) - sumn
        return (1. - np.exp(-alpha*x)) * U_vdw
    Us = Ushort(r/R_e)
    Ul = Ulong(r/R_e)
    return 27.2114*((Z_A*Z_B/R_e)*Us + D_e*Ul)

def ULiHe(r):
    r = r / 0.5291772
    
    R_e = 11.47
    D_e = 7.36e-6
    Z_A = 2.
    Z_B = 3.
    
    a1 = 11.49538
    a2 = -25.99090
    a3 = 13.49545
    alpha = 25.06466
    lambd = 14.57866
    gamma = 2.19783
    b = 12.38083
    A = 2.31539e6
    c = np.zeros(20)
    c[6] = 1.43308
    c[8] = 0.47909
    c[10] = 0.22501
    c[12] = 0.14847
    c[14] = 0.13762
    c[16] = 0.17921
    
    def Ushort(x):
        return (1./x) * (1. + a1*x + a2*x**2 + a3*x**3)*np.exp(-alpha*x)
    def Ulong(x):
        sumn = 0.
        for n in range(3,9):
            sumk = 1.
            factk = 1.
            for k in range(1,2*n+1):
                factk = factk * k
                sumk = sumk + (b*x)**k / factk
            sumn = sumn + (1. - np.exp(-b*x)*sumk)*c[2*n]/(x**(2*n))
        U_vdw = A * x**gamma * np.exp(-lambd*x) - sumn
        return (1. - np.exp(-alpha*x)) * U_vdw
    Us = Ushort(r/R_e)
    Ul = Ulong(r/R_e)
    return 27.2114*((Z_A*Z_B/R_e)*Us + D_e*Ul)

def UHeHe(r):
    r = r / 0.5291772
    
    R_e = 5.608
    D_e = 3.482e-5
    Z_A = 2
    Z_B = 2
    
    a1 = 10.34329
    a2 = -23.68667
    a3 = 12.34334
    alpha = 22.76733
    lambd = 15.14296
    gamma = 1.59236
    b = 13.55060
    A = 3.5552e6
    c = np.zeros(20)
    c[6] = 1.34992
    c[8] = 0.41469
    c[10] = 0.17155
    c[12] = 0.09557
    c[14] = 0.07170
    c[16] = 0.07244
    
    def Ushort(x):
        return (1./x) * (1. + a1*x + a2*x**2 + a3*x**3)*np.exp(-alpha*x)
    def Ulong(x):
        sumn = 0.
        for n in range(3,9):
            sumk = 1.
            factk = 1.
            for k in range(1,2*n+1):
                factk = factk * k
                sumk = sumk + (b*x)**k / factk
            sumn = sumn + (1. - np.exp(-b*x)*sumk)*c[2*n]/(x**(2*n))
        U_vdw = A * x**gamma * np.exp(-lambd*x) - sumn
        return (1. - np.exp(-alpha*x)) * U_vdw
    Us = Ushort(r/R_e)
    Ul = Ulong(r/R_e)
    return 27.2114*((Z_A*Z_B/R_e)*Us + D_e*Ul)

def UPbHe(r):
        Re = 4.6874
        De = 39.4116
        a = [5.4114, -20.7496, 56.2663, -95.6327, 112.157, -77.8808, 28.2574, -3.7882]
        x = (r-Re)/Re
        V = 1.
        for kk in range(8):
            k = kk + 1
            V += a[kk]*x**k
        V = -De*np.exp(-a[0]*x)*V
        return V*1e-6*27.211386245988

if __name__=="__main__":
    import matplotlib.pyplot as plt
    from ALTRES import CONFIGURA_PLT
    CONFIGURA_PLT()
    r = np.arange(4, 15, 0.001)
    plt.plot(r, UPbHe(r))
    plt.xlabel("$r \; ({\\rm \\AA})$")
    plt.ylabel("$V_{\\rm Pb-He} (r) \\; (eV)$")
    plt.show()

    r = np.arange(5, 15, 0.001)
    plt.plot(r, ULiHe(r))
    plt.xlabel("$r \; ({\\rm \\AA})$")
    plt.ylabel("$V_{\\rm Li-He} (r) \\; (eV)$")
    plt.show()