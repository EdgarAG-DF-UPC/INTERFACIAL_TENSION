#!/usr/bin/env python
# coding: utf-8

import numpy as np
import sys
with open("../../parent-directory", "r") as pdir: PARENTDIR = pdir.read().replace('\n','').replace('"','')
# sys.path.append('../../../../../SCRIPTS_LiPbHe')
sys.path.append(PARENTDIR+"/SCRIPTS_LiPbHe/src")
import LATTICE

def dens(T,TL,A,B):
    return A - B * (T - TL)

# PARÀMETRES GENERALS

#Massa de cada tipus d'àtom (Li, Pb, He):
m_Li=6.941
m_He=4.002
m_Pb=207.20

#Potencials LJ:
#cutoff=9.25
cutoff=8.00
### ORIGINAL PARAMETERS (JM)
###epsilon_LiLi=0.07183770
###epsilon_HeHe=0.00092031
###epsilon_PbPb=0.06098297
###epsilon_LiHe=0.00019626
###epsilon_LiPb=0.06618842
###epsilon_HePb=0.00105286
###sigma_LiLi=2.800	
###sigma_HeHe=2.637	
###sigma_PbPb=3.689	
###sigma_LiHe=5.320	
###sigma_LiPb=3.244	
###sigma_HePb=4.172
###
###AZIZ-FRAILE-POLCAR-JM
#epsilon_LiLi=0.07183770
epsilon_LiLi=0.07651545
epsilon_HeHe=0.00088069
epsilon_PbPb=0.06098297
epsilon_LiHe=0.00014100
epsilon_LiPb=0.06618842
epsilon_HePb=0.01057000
#sigma_LiLi=2.800
sigma_LiLi=2.7242
sigma_HeHe=2.556	
sigma_PbPb=3.689	
sigma_LiHe=5.3565	
sigma_LiPb=3.244
sigma_HePb=3.0667
###

#Pas de temps i nombre de passos d'integració:
print("Pas de temps, i nombre de passos (NVT, NPT):")
delta_t, runs_NVT, runs_NPT, Nevery = input().split()
# delta_t=0.002
# runs_NVT=25000
# runs_NPT=1000000
Tdamp=20*float(delta_t)
Pdamp=20*float(delta_t)


#Col·lectivitat:
print("Quina col·lectivitat vols simular (NVT o NPT) ?")
ensemble=input()
while ensemble not in ["NVT", "NPT"]:
    print("Si us plau, insereix una col·lectivitat vàlida")
    print("NVT o NPT")
    ensemble=input()
print("-->", ensemble)
print("")

#Model aliatge:
model = ""
print("Tria un model per a els aliatges Li-Pb (AWAD o BELASHCHENKO o FRAILE o ...)")
while model not in ["AWAD.LiPb", "BELASHCHENKO.LiPb", "FRAILE.LiPb", "Pb-Bela12_Li-Awad23", "Pb-Bela12_Li-Bela11", "Pb-Zhou01_Li-Awad23", "Pb-Zhou04_Li-Awad23", "Pb-Bela12_Li-Bela11_LJ84_19-Yukawa", "Pb-Zhou01_Li-Awad23_Morse23"]:
    model = input()
print("-->", model)
if model in ["BELASHCHENKO.LiPb", ]:
    print("Se simularà el model de Belashchenko. Vols incloure la correcció de Yukawa? (sel·lecciona SI o NO)")
    incloure=input() # Yukawa correction - if false, then Belashchenko transferable option
    if incloure=="SI":
        print("Has triat l'opció amb correcció de Yukawa.")
        correccio = True
    else:
        if incloure=="NO":
            print("Has triat l'opció transferible.")
        else:
            print("ERROR. Per defecte s'ha assignat l'opció transferible.")
        correccio = False
    print("")
elif model=="AWAD.LiPb":
    print("Se simularà el model Li-Awad-LRO (2023) -- Pb-Belashchenko (2012)")
    dummy = input()
    incloure = "NO"
    correccio = False
elif model=="FRAILE.LiPb":
    print("Se simularà el model Fraile et al. (2011)")
    dummy = input()
    incloure = "NO"
    correccio = False
else:
    print("Se simularà el model", model)
    dummy = input()
    incloure = "NO"
    correccio = False
print("-->", incloure, "s'inclourà la correcció de Yukawa")
print("")

# Densitat (volum) inicialització:
DENSITAT = ""
print("Tria inicialització densitat (BELASHCHENKO o KHAIRULIN)")
while DENSITAT not in ["BELASHCHENKO", "KHAIRULIN", "KHAIRULIN_INTERPOLACIO", "INIT_Li"]:
    DENSITAT = input()
print("-->", DENSITAT)
print("")

# Tipus de xarxa (del sòlid inicial)
LATTICE_STYLE = ""
print("Tria el tipus de xarxa (SC o BCC o FCC)")
while LATTICE_STYLE not in ["SC", "BCC", "FCC"]:
    LATTICE_STYLE = input()
if LATTICE_STYLE=="SC":
    Ncp = 1
elif LATTICE_STYLE=="BCC":
    Ncp = 2
elif LATTICE_STYLE=="FCC":
    Ncp = 4
print("---> S'originarà el sistema amb una xarxa de tipus ", LATTICE_STYLE, "la qual conté", Ncp, "àtoms per cel·la unitat")

#Nombre total d'àtoms a les simulacions:
print("Nombre total d'àtoms de metall líquid a les simulacions:")
try:
    N=int(input())
except ValueError as err: # default -- en cas que hi hagi un error en l'entrada, l'escollim per defecte
    if LATTICE=="SC":
        N = 1000
    elif LATTICE=="BCC":
        N = 1024
    elif LATTICE=="FCC":
        N = 1372
    print("ERROR en l'entrada de N. Es procedirà amb el valor per defecte N=",N)
print(N, "--->", round((N/Ncp)**(1/3)))
print("")


# #Nombre d'àtoms d'heli:
# print("Nombre d'àtoms d'heli:")
# try:
#     NHe=int(input())
# except ValueError as err:
#     print("ERROR en l'entrada de N. Es procedirà amb el valor per defecte N=",N)
#     quit()

MINIMIZATION = ""
print("Vols aplicar la minimització (inicial) de les coordenades - gradient conjugat?")
while MINIMIZATION not in ["SI", "NO"]:
    MINIMIZATION = input()
if MINIMIZATION=="SI":
    MINIMIZATION = True
elif MINIMIZATION=="NO":
    MINIMIZATION = False
else:
    print("Això no hauria de passar mai...")
    quit()
    
#Equilibració amb Lennard-Jones:
print("Vols inicar l'equilibració amb model Lennard-Jones ?")
LJ=input()
if LJ=="SI":
    print("S'inicialitzará amb LJ, i després se substituirà el model a EAM + TTS.")
    equil_LJ = True
    # with open("../SCRIPTS_LiPbHe/EQUILIBRATION_LENNARD_JONES", "r") as eqLJ:
    #     comandes_LJ = eqLJ.read()
else:
    print("S'inicialitzarà directament amb els models EAM + TTS.")
    equil_LJ = False
    comandes_LJ = ""
print("-->", LJ)
print("")

#Pressions parcials:
print("Vols calcular les pressions parcials de cada espècie ?")
PP=input()
print("")
# if PP=="SI":
#     with open("../SCRIPTS_LiPbHe/PRESSIONS_PARCIALS", "r") as ppfile:
#         comandes_PP = ppfile.read()
# else:
#     comandes_PP = ""
print("-->", PP)
print("")

#Energies extra:
print("Vols calcular les energies extra ?")
EE=input()
print("")
if EE=="SI":
    with open("../SCRIPTS_LiPbHe/ENERGIES_EXTRA", "r") as eefile:
        comandes_enersextra = eefile.read()
else:
    comandes_enersextra = ""
print("-->", EE)
print("")

#Difusions (Li, He, Pb):
print("Vols calcular les difusions de cada espècie ?")
DIFS=input()
print("")
# if DIFS=="SI":
#     with open("../SCRIPTS_LiPbHe/DIFUSIONS_Li_He_Pb", "r") as difsfile:
#         comandes_difs = difsfile.read()
# else:
#     comandes_difs = ""
# print("-->", DIFS)

print("Vols sacsejar tèrmicament el sistema?")
SHAKE = input()
while SHAKE not in ["SI","NO"]:
    print("Si us plau, tria SI o NO")
    SHAKE = input()
if SHAKE=="SI":
    with open("../SCRIPTS_LiPbHe/THERMAL_SHAKE", "r") as shakefile:
        comandes_shake = shakefile.read()
    try:
        TS, NSS = input().split()
        float(TS)
        float(NSS)
        
        comandes_shake = comandes_shake.replace("Tshake", TS)
        comandes_shake = comandes_shake.replace("runs_NVT", NSS)
    except ValueError as err:
        print("COMPTE, no s'ha introduit una temperatura vàlida. No es farà l'agitació tèrmica...")
        SHAKE = "NO"
        comandes_shake = ""
else:
    comandes_shake = ""

print("A quina partició vols enviar el càlcul?")
PARTICIO = input()
while PARTICIO not in ["all", "brief", "fast"]:
    print("Introdueix una partició vàlida:")
    print("--> all")
    print("--> brief")
    print("--> fast")
    PARTICIO = input()
print("-->", PARTICIO)
print("")

print("Vols reemplaçar els directoris existents (es perdrà la informació) o fer una còpia de seguretat (moure els arxius) ?")
print("Tria entre 'MOVE' i 'DELETE'")
MOVE=False
DELETE=False
while not MOVE and not DELETE:        
    accio=input()
    if accio=="MOVE":
        MOVE=True
    elif accio=="DELETE":
        DELETE=True
    else:
        print("Si us plau, tria entre 'MOVE' o 'DELETE' ")
print("-->", accio)
print("")


import os
import shutil

FRAC_ATOM = []
with open("FRAC_ATOM", "r") as fracs: # També llegirem totes les fraccions atòmiques de l'arxiu generat a 'input.sh'
    for line in fracs:
        if "#" not in line:
            #FRAC_ATOM.append(float(line))
            FRAC_ATOM.append(line.split()[0])
HeATS = []
with open("HELIUM_ATS", "r") as heats:
    for line in heats:
        if "#" not in line:
            HeATS.append(line.split()[0])

if MOVE: # En cas d'haver triat de fer una còpia dels directoris existents abans de ser reemplaçats:
    try:
        os.mkdir("SIMULACIONS_ANTIGUES") # 1. Si es tracta del primer escombrat de simulacions, intentarem crear un nou directori on guardar les còpies.
        print("S'ha creat el directori 'SIMULACIONS_ANTIGUES'")
    except FileExistsError as err:
        print("Les simulacions antigues es guardaran dins del directori existent 'SIMULACIONS_ANTIGUES'")

    print("Anomena el directori per guardar simulacions antigues:") # 2. Anomenem el subdirectori on fer la còpia. Demanarem el nom tantes vegades com sigui necessari per a evitar sobreescriure un directori ja existent.
    demana_nom = True
    while demana_nom:
        NOM = input()
        try:
            os.mkdir("SIMULACIONS_ANTIGUES/"+NOM)
            demana_nom = False # Romandrà False mentre que no passem per l'excepció.
        except FileExistsError as err:
            print(err)
            print(os.listdir("SIMULACIONS_ANTIGUES"))
            print("Si us plau, torna a introduir un nom:")
            demana_nom = True #Mentre que surti l'error, tornarem a demanar el nom.
    
    try:
        shutil.copyfile("INPUTS", "SIMULACIONS_ANTIGUES/"+NOM+"/INPUTS")
        shutil.copyfile("RUN.INFO", "SIMULACIONS_ANTIGUES/"+NOM+"/RUN.INFO")
        shutil.copyfile("TEMPERATURES", "SIMULACIONS_ANTIGUES/"+NOM+"/TEMPERATURES")
        shutil.copyfile("FRAC_ATOM", "SIMULACIONS_ANTIGUES/"+NOM+"/FRAC_ATOM")
    except FileNotFoundError as err:
        print("COMPTE !!! --- ", err)

try:
    BLOCKS = input()
    SET_BLOCKS = True
except EOFError as err:
    SET_BLOCKS = False
try:
    BORN = input()
    if BORN in ["FORWARD", "+", "POSITIVE", "REPULSIVE"]:
        signe = "+"
    elif BORN in ["BACKWARD", "REVERSE", "-", "NEGATIVE", "ATTRACTIVE"]:
        signe = "-"
    else:
        raise NameError("Invalid variable name...")
    TUNING = True
except (NameError, EOFError) as err:
    TUNING = False

num_simuls_moved = 0 # Comptador del nombre de simulacions prèvies que hem mogut abans de llançar les noves simulacions.

for x_Li in FRAC_ATOM:
    x_Pb = "{:.2f}".format(100. - float(x_Li))
    parent = "Li"+x_Li+"-Pb"+x_Pb
    if DELETE:
        try:
            os.rmdir(parent)
        except FileNotFoundError as o:
            print("No s'ha trobat el directori "+parent+"/")
        except OSError as o:
            print("Directori "+parent+"/: no buit")
            shutil.rmtree(parent)
            os.mkdir(parent)
            shutil.copyfile("TEMPERATURES", parent+"/TEMPERATURES")
    elif MOVE:
        try:
            shutil.move(parent, "SIMULACIONS_ANTIGUES/"+NOM+"/")
            num_simuls_moved += 1 # només se sumarà si l'anterior comanda no dona error
        except FileNotFoundError as err1:
            print(err1)
            print("No s'ha trobat l'arxiu")
        except shutil.Error as err2:
            print(err2)
            print("Ja existia un subdirectori amb aquest nom") #això no hauria de passar mai (ja que abans demanem canviar el nom)!!!
        
    else:
        print("Hauria de ser MOVE o DELETE...") # Això no hauria de passar mai!!!!

    TEMPERS = []
    with open(parent+"/TEMPERATURES", "r") as tempers: # Llegirem totes les temperatures de l'arxiu generat a 'input.sh'
        for line in tempers:
            if "#" not in line:
                TEMPERS.append(float(line))
    
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

    for N_He in HeATS:
        os.mkdir(parent+"/He"+str(N_He)+"/")
        for temp in TEMPERS:
            os.mkdir(parent+"/He"+str(N_He)+"/T_"+"{:.2f}".format(temp)+"K/")
            for pres in [1]:
#            for pres in [500]:
                os.mkdir(parent+"/He"+str(N_He)+"/T_"+"{:.2f}".format(temp)+"K/P_"+str(pres)+"bar/")
                
                directori = parent+"/He"+str(N_He)+"/T_"+"{:.2f}".format(temp)+"K/P_"+str(pres)+"bar/"
                print("------------------")
                print(directori)
                with open(directori+"NxPT","w") as fitxer:
                    fitxer.write(str(N)+"\n")
                    fitxer.write(str(x_Li)+"\n")
                    fitxer.write(str(pres)+"\n")
                    fitxer.write(str(temp)+"\n")
                
                #===========================================================================================================================================

                with open(directori+"/RUN.INPUT", "w") as inputfile:
                    inputfile.write("RUNS_NVT = "+str(runs_NVT)+"\n")
                    inputfile.write("RUNS_NPT = "+str(runs_NPT)+"\n")
                    inputfile.write("TIMESTEP = "+str(delta_t)+"\n")
                    inputfile.write("NEVERY = "+str(Nevery)+"\n")
                    inputfile.write("ENSEMBLE = "+ensemble+"\n")
                    inputfile.write("MODEL = "+model+"\n")
                    if model=="BELASHCHENKO": inputfile.write("YUKAWA_CORRECTION = "+incloure+"\n")
                    inputfile.write("DENSINI = "+DENSITAT+"\n")
                    inputfile.write("LATTICE = "+LATTICE_STYLE+"\n")
                    inputfile.write("LJ = "+LJ+"\n")
                    inputfile.write("PP = "+PP+"\n")
                    inputfile.write("EE = "+EE+"\n")
                    inputfile.write("DD = "+DIFS+"\n")
                    inputfile.write("SHAKE = "+SHAKE+"\n")
                    if SHAKE=="SI": inputfile.write("Temp_Shake = "+str(TS)+"\n"+"Nsteps_Shake = "+str(NSS)+"\n")
                    inputfile.write("MIN = "+str(MINIMIZATION)+"\n")
                    inputfile.write("QUEUE = "+PARTICIO+"\n")
                    if SET_BLOCKS: inputfile.write("NBLOCKS = "+BLOCKS+"\n")
                    if TUNING: inputfile.write("SIGN = "+signe+"\n")
                #===========================================================================================================================================

                # fitxer1 = "../SCRIPTS_LiPbHe/in.LiPb.Hebubbles"
                # fitxer2 = directori+"in.LiHePb"
                # shutil.copyfile(fitxer1,fitxer2)
                
                if correccio:
                    X = min([0, 10, 20, 30, 40, 50], key=lambda x:abs(x-float(x_Pb)))
                    if int(X)==0 or float(x_Pb)>=60:
                        scpm="0"
                        qLi=0.
                        qPb=0.
                    elif int(X) == 10:
                        scpm="-1.56811644" # 14.3996 * 0.11 * 0.99 eV
                        qLi=0.11
                        qPb=-0.99
                    elif int(X) == 20:
                        scpm="-0.97341296" # 14.3996 * 0.13 * 0.52 eV
                        qLi=0.13
                        qPb=-0.52
                    elif int(X) == 30:
                        scpm="-0.6679830444" # 14.3996 * 0.141 * 0.329 eV
                        qLi=0.141
                        qPb=-0.329
                    elif int(X) == 40:
                        scpm="-0.1853430114" # 14.3996 * 0.0926 * 0.139 eV
                        qLi=0.0926
                        qPb=-0.1389
                    elif int(X) == 50:
                        scpm="-0.08099775" # 14.3996 * 0.075 * 0.075 eV
                        qLi=0.075
                        qPb=-0.075
                    else:
                        print("ERROR !!!!!!!")

                    print(qLi, qPb, "--->", scpm, "eV")
                    
                # x = float(x_Li)
                # NLi = round(x*N/100.)
                # NPb = N - NLi
                P = pres
                T = temp
                # M = (NLi*6.941 + NPb*207.2) / (NLi+NPb)

                C, a, DN, NLi, NPb, Tmin = LATTICE.LATTICE(LATTICE_STYLE, N, float(x_Li), temp, DENSITAT, False)
                NHe = int(N_He)

                with open(directori+"NPT","w") as fileout:
                    fileout.write(str(NLi)+"\n"+str(NHe)+"\n"+str(NPb)+"\n"+str(P)+"\n"+str(T)+"\n")
                    print(str(NLi)+"\n"+str(NHe)+"\n"+str(NPb)+"\n"+str(P)+"\n"+str(T)+"\n")

                #RUN = "True"
                if T>=Tmin:
                    RUN = True
                else:
                    RUN = False
                print("La simulació correrà?  --->  ", RUN)
                with open(directori+"RUN", "w") as fitxer:
                    fitxer.write(str(RUN)+"\n")
                
                # with open(fitxer2,"r") as file:
                #     dades = file.read()
                    
                #     dades = dades.replace("EQUILIBRATION_LENNARD_JONES", comandes_LJ)
                #     dades = dades.replace("THERMAL_SHAKE", comandes_shake)
                #     dades = dades.replace("PRESSIONS_PARCIALS", comandes_PP)
                #     dades = dades.replace("ENERGIES_EXTRA", comandes_enersextra)
                #     dades = dades.replace("DIFUSIONS_Li_He_Pb", comandes_difs)

                #     if ensemble=="NVT":
                #         dades = dades.replace("comanda_fix_ensemble", "fix 1 all nvt temp temperatura temperatura Tdamp")
                #     elif ensemble=="NPT":
                #         dades = dades.replace("comanda_fix_ensemble", "fix 1 all npt temp temperatura temperatura Tdamp iso pressio  pressio  Pdamp")
                #     else:
                #         print("ERROR: NVT o NPT !!!!!!!!!") # Això no hauria de passar mai.

                #     if correccio:
                #         with open("SCRIPTS/POTENCIAL_YUKAWA", "r") as yukawa:
                #             potencial = yukawa.read()
                #         dades = dades.replace("AFEGEIX_YUKAWA", potencial)
                #         dades = dades.replace("screening",str(scpm))
                #         dades = dades.replace("qLi",str(qLi))    
                #         dades = dades.replace("qPb",str(qPb))
                #     else:
                #         dades = dades.replace("AFEGEIX_YUKAWA", "")


                #     dades = dades.replace("temperatura","{:.2f}".format(temp))
                #     dades = dades.replace("pressio","{:.2f}".format(pres))
                #     dades = dades.replace("rcutoff","{:.2f}".format(cutoff))
                #     dades = dades.replace("N_Li",str(NLi))
                #     dades = dades.replace("N_He",str(NHe))
                #     dades = dades.replace("N_Pb",str(NPb))
                #     dades = dades.replace("m_Li",str(m_Li))
                #     dades = dades.replace("m_He",str(m_He))
                #     dades = dades.replace("m_Pb",str(m_Pb))
                #     dades = dades.replace("epsilon_LiLi",str(epsilon_LiLi))
                #     dades = dades.replace("epsilon_HeHe",str(epsilon_HeHe))
                #     dades = dades.replace("epsilon_PbPb",str(epsilon_PbPb))
                #     dades = dades.replace("epsilon_LiHe",str(epsilon_LiHe))
                #     dades = dades.replace("epsilon_LiPb",str(epsilon_LiPb))
                #     dades = dades.replace("epsilon_HePb",str(epsilon_HePb))
                #     dades = dades.replace("sigma_LiLi",str(sigma_LiLi))
                #     dades = dades.replace("sigma_HeHe",str(sigma_HeHe))
                #     dades = dades.replace("sigma_PbPb",str(sigma_PbPb))
                #     dades = dades.replace("sigma_LiHe",str(sigma_LiHe))
                #     dades = dades.replace("sigma_LiPb",str(sigma_LiPb))
                #     dades = dades.replace("sigma_HePb",str(sigma_HePb))
                #     dades = dades.replace("delta_t",str(delta_t))
                #     dades = dades.replace("a_lattice",str(a))
                #     dades = dades.replace("x0",str(-C/2.))
                #     dades = dades.replace("xf",str(C/2.))
                #     dades = dades.replace("y0",str(-C/2.))
                #     dades = dades.replace("yf",str(C/2.))
                #     dades = dades.replace("z0",str(-C/2.))
                #     dades = dades.replace("zf",str(C/2.))
                #     dades = dades.replace("DN",str(DN))
                #     dades = dades.replace("runs_NVT",str(runs_NVT))
                #     dades = dades.replace("runs_NPT",str(runs_NPT))
                #     dades = dades.replace("Tdamp",str(Tdamp))
                #     dades = dades.replace("Pdamp",str(Pdamp))
                #     dades = dades.replace("RBOMB", "{:.1f}".format(C*a/9.))

                
                # with open(directori+"in.LiHePb","w") as file:
                #     file.write(dades)
                
                
                fitxer3 = directori+"LiPb.eam.alloy"
                if model=="Pb-Bela12_Li-Bela11_LJ84_19-Yukawa":
                    if float(x_Pb)<=50 and float(x_Pb)!=0.:
                        X = min([10, 20, 30, 40, 50], key=lambda x:abs(x-float(x_Pb)))
                    else:
                        X = 0.
                    shutil.copyfile("../../TAULES/"+model+"{:.1f}".format(1.-X/100.)+".eam.alloy",fitxer3)
                    print("../../TAULES/"+model+"{:.1f}".format(1.-X/100.)+".eam.alloy")
                    print("    ")
                    print("    ")
                    print("    ")
                    print("    ")
                    print("    ")
                # elif correccio and float(x_Pb)<=50 and float(x_Pb)!=0.:
                #     X = min([10, 20, 30, 40, 50], key=lambda x:abs(x-float(x_Pb)))
                #     shutil.copyfile("../TAULES/"+model+".eam.alloy", fitxer3)
                else:
                    shutil.copyfile("../../TAULES/"+model+".eam.alloy",fitxer3)
                with open(fitxer3,"r") as file:
                    dades = file.read()
                    dades = dades.replace("amassLi",str(m_Li))
                    dades = dades.replace("amassPb",str(m_Pb))
                    dades = dades.replace("alat",str(a))
                with open(fitxer3,"w") as file:
                    file.write(dades)
                    
                # shutil.copy("../SCRIPTS_LiPbHe/escombrat.VIRTCCQM.sub", directori)
                # fitxer4 = directori+"escombrat.VIRTCCQM.sub"
                # with open(directori+"escombrat.VIRTCCQM.sub","r") as file:
                #     dades = file.read()
                #     dades = dades.replace("NomCalcul",str(T)+"K_"+x_Li+"_"+x_Pb+"_ALLOY")
                #     dades = dades.replace("NomParticio", PARTICIO)
                # with open(directori+"escombrat.VIRTCCQM.sub","w") as file:
                #     file.write(dades)

                shutil.copy("../../SCRIPTS_LiPbHe/onthefly.sub", directori)
                fitxer4 = directori+"onthefly.sub"
                with open(directori+"onthefly.sub","r") as file:
                    dades = file.read()
                    dades = dades.replace("NomCalcul",str(T)+"K_"+x_Li+"_"+x_Pb+"_ALLOY")
                    dades = dades.replace("NomParticio", PARTICIO)
                with open(directori+"onthefly.sub","w") as file:
                    file.write(dades)

                shutil.copy("../../SCRIPTS_LiPbHe/testarea.sub", directori)
                fitxer4 = directori+"testarea.sub"
                with open(directori+"testarea.sub","r") as file:
                    dades = file.read()
                    dades = dades.replace("NomCalcul",str(T)+"K_"+x_Li+"_"+x_Pb+"_TEST_AREA")
                    dades = dades.replace("NomParticio", PARTICIO)
                with open(directori+"testarea.sub","w") as file:
                    file.write(dades)
                
                print("")
                print("------------------")

if MOVE:
    print("S'han mogut un total de "+str(num_simuls_moved)+" directoris")
    if num_simuls_moved==0: os.rmdir("SIMULACIONS_ANTIGUES/"+NOM)
