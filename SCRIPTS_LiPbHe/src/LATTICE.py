import numpy as np

def dens(T,TL,A,B):
    return A - B * (T - TL)

def LATTICE(STYLE, N, x, T, DENSITAT, show):
    x_Li = x / 100.
    NLi = round(x_Li*N)
    NPb = N - NLi
    M = (NLi*6.941 + NPb*207.2) / (NLi+NPb)

    if STYLE=="SC":
        Ncp = 1
    elif STYLE=="BCC":
        Ncp = 2
    elif STYLE=="FCC":
        Ncp = 4
    else:
        print("¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡")
        print(STYLE)
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        quit()

    if DENSITAT=="BELASHCHENKO":
        x = int(float(x_Li))
        #x = min([0, 10, 20, 30, 40, 50], key=lambda x:abs(x-float(x_Pb)))
#                # [Vmol] = cm³/mol
        if x==100:
            Vmol = 14.96
        elif x==90:
            Vmol = 14.64
        elif x==80:
            Vmol = 14.54
        elif x==70:
            Vmol = 14.72
        elif x==60:
            Vmol = 15.15
        elif x==50:
            Vmol = 15.80
        elif x==40:
            Vmol = 16.60
        elif x==30:
            Vmol = 17.50
        elif x==20:
            Vmol = 18.42
        elif x==17:
            Vmol = 18.72
        elif x==10:
            Vmol = 19.39
        elif x==0:
            Vmol = 20.37
        else:
            print("ERROR")
        # Tmin = min(TEMPERS)
        rho = 0.6022 / Vmol # Å^{-3}
        x /= 100.

    elif DENSITAT=="KHAIRULIN":    
        xx_Li = x / 100.
        with open("/users/edgar/BUBBLES_NEW/TAULES/dens.LiPb.param") as params:
            for line in params:
                if "#" not in line:
                    line = line.split()
                    xx_Pb = float(line[0]) / 100.
                    if "{:.4f}".format(xx_Pb) == "{:.4f}".format(1.00-xx_Li):
                        TL = float(line[1])
                        A = float(line[2])
                        B = float(line[3])
                        rho = dens(T,TL,A,B) * 6.022/M * 0.0001 # Angs^-3
                        Tmin = float(line[4])
                        Tmax = float(line[5])
                    elif "{:.2f}".format(xx_Pb) == "100.00":
                        rho = 10.678 - 13.174e-4*(T-600.6) # g / cm³ - Kirshenbaum, Cahill, Grosse - 1961
                        rho *= 6.022/M * 0.1 # Angs^-3
                        Tmin = 600.6 #K
                        Tmax = 2024 #K

        x = xx_Li
    
    elif DENSITAT=="KHAIRULIN_INTERPOLACIO":
        xx_Li = x / 100.
        x_Khairulin = []
        densitat_Khairulin = []
        with open("/users/edgar/BUBBLES_NEW/TAULES/dens.LiPb.param") as params:
            for line in params:
                if "#" not in line:
                    line = line.split()
                    xx_Pb = float(line[0]) / 100.

                    TL = float(line[1])
                    A = float(line[2])
                    B = float(line[3])
                    x_Khairulin.append(xx_Pb)
                    densitat_Khairulin.append(dens(T,TL,A,B))
                    Tmin = float(line[4])
                    Tmax = float(line[5])
        x_Khairulin.append(1.)
        densitat_Khairulin.append((10.678 - 13.174e-4*(T-600.6)) * 1e3)

        rho = np.interp(1.-xx_Li, x_Khairulin, densitat_Khairulin) * 6.022/M * 0.0001 # Angs^-3

        x = xx_Li

    elif DENSITAT=="INIT_Li":
        # x = int(float(x_Li))
        x = x_Li
        # rho = 10.52*(1.-113e-6*T) * 6.022e-1 / 6.941 # Angs^-3
        rho = ( 562. - 0.100*T ) *6.022e-4 / 6.941 # Angs^-3
        Tmin = 453. #K
        # NLi = round(x_Li*N)
        # NPb = N - NLi

    else:
        raise NameError("Variable DENSITAT cannot be "+str(DENSITAT)+"...")
    
    if show: print("Núm. cel·les a cada direcció:   C = L / a = ", round((float(NLi+NPb)/float(Ncp))**(1./3.))) # nombre de cel·les a cada direcció == C = L / a = int{(N/2)^3}
    C = round((float(NLi+NPb)/float(Ncp))**(1./3.)) # C = L/a
    
    a = (float(Ncp)/rho)**(1./3.) # lattice spacing
    if show: print("Lattice spacing = ", a, " \AA")

    # Redefinim el nombre total d'àtoms (de manera que (N/2)^3 sigui un nombre enter): N' = 2 C^3 = 2*(int(((NLi+NPb)/2)**(1./3.)))**3 
    N = Ncp*(round(((NLi+NPb)/Ncp)**(1./3.)))**3 
    # ---> D = int(N'/NLi') és la separació entre les etiquetes dels àtoms de liti (la resta seran de plom)
    #D = int( 2*(int(((NLi+NPb)/2)**(1./3.)))**3/NLi ) # 2 x C^3 / NLi
    # El nombre d'àtoms de liti es redefineix com: NLi' = x N' per tal de preservar la fracció de cada component.
    # Ídem pel plom: NPb' = N' - NLi'
    NLi = round(x*N)
    NPb = N - NLi
    
    # -->
    if NPb == 0: 
        DN = 0 
    else:
        DN = 1./(1.-x)

    if show:
        print("")
        print("")
        print("")

    return C, a, DN, NLi, NPb, Tmin