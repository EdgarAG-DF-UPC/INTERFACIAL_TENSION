import numpy as np
from scipy.optimize import curve_fit
from mpi4py import MPI
nprocs = MPI.COMM_WORLD.Get_size()
me = MPI.COMM_WORLD.Get_rank()
master = me==0
import os, shutil, sys, subprocess
with open("parent-directory", "r") as pdir: PARENTDIR = pdir.read().replace('\n','').replace('"','')
# sys.path.append('../../../../../SCRIPTS_LiPb')
sys.path.append(PARENTDIR+"/SCRIPTS_LiPbHe/src")
import lmpcmds.COMANDES as COMANDES, lmpcmds.ADDITIONAL as ADDITIONAL, lmpcmds.FUNCIONS as FUNCIONS

from src.CONFIGURACIO import CUA, RESTART, RSTFILE, METHOD as WHICH, BLOCKS, AFEGEIX_YUKAWA, runs_NVT, runs_NPT, Nevery, x0, y0, z0, xf, yf, zf, RBOMB, temperatura, Tdamp, pressio, Pdamp, ensemble, N_Li, N_Pb, N_He, rcutoff, DIFUSIONS_Li_He_Pb, PRESSIONS_PARCIALS, MODEL
RESTART = FUNCIONS.evalua(RESTART)
DELTA = float(runs_NPT) / float(Nevery)

def printonmaster(string):
    if master: print(string)


class METHOD:
    def THOMPSON(lmp,
                #  temperatura: float, pressio: float, Tdamp: float, Pdamp: float,
                #  rcutoff: float,
                #  runs_NVT: int, runs_NPT: int, Nevery: int,
                #  ensemble="NPT", RESTART=False, AFEGEIX_YUKAWA=False, PRESSIONS_PARCIALS=False, DIFUSIONS_Li_He_Pb=False,
                TOT=25*nprocs):
        from src.CONFIGURACIO import BLOCK0, BLOCKF, INTEGRATE, proc0, Nevery, MOVING_AVG
        ##################################################################################################################
        def MVA(data, Nevery):
            Nout = int(len(data) / Nevery)
            assert len(data)%Nevery==0
            return np.array([np.mean(data[k*Nevery:(k+1)*Nevery]) for k in range(Nout)])
        ##################################################################################################################
        if master:
            N_LI = []
            N_PB = []
            N_HE = []
            P_K = []
            P_C = []
            γS = []
            RS = []
            Δp = []
            FUNCIONS.CreateExecutable(PARENTDIR+"/INTERFACIAL_TENSION/",
                                      "spherical_tension.x",
                                      "tensio.x")
            FUNCIONS.CreateExecutable(PARENTDIR+"/INTERFACIAL_TENSION/",
                                      "spherical_integral.x",
                                      "spherical_integral.x")
        ##################################################################################################################
        def ANALITZA(block: int, proc: int, N_LI: list, N_PB: list, N_HE: list, P_K: list, P_C: list, RS: list, Δp: list):
            try:
                with open("TENSIO"+str(block)+"_"+str(proc)+"/r_rho_pv.dat") as filein:
                    r = []
                    N_Li = []
                    N_Pb = []
                    N_He = []
                    p_kine = []
                    p_conf = []
                    for line in filein: 
                        if "#" not in line:
                            line = line.split()
                            r.append(float(line[0]))
                            N_Li.append(float(line[1]))
                            N_Pb.append(float(line[2]))
                            N_He.append(float(line[3]))
                            p_kine.append(float(line[4]))
                            p_conf.append(float(line[5]))
                    N_LI.append(N_Li)
                    N_PB.append(N_Pb)
                    N_HE.append(N_He)
                    P_K.append(p_kine)
                    P_C.append(p_conf)
                    os.chdir("TENSIO"+str(block)+"_"+str(proc)+"/")
                    os.system(PARENTDIR+"/INTERFACIAL_TENSION/spherical_integral.x > ST.out")
                    with open("ST.out", "r") as filein:
                        for line in filein:
                            if "N/m" in line and line.split()[2]!="NaN":
                                #continue
                                γS.append(float(line.split()[2]))
                            elif "Radius" in line and line.split()[0]!="NaN":
                                RS.append(float(line.split()[0]))
                            elif "eV Å^{-3}" in line and line.split()[0]!="NaN":
                                Δp.append(-float(line.split()[0]))
                    os.chdir("../")
            except FileNotFoundError as err:
                print(err)
                print("No s'ha trobat l'arxiu TENSIO"+str(block)+"_"+str(proc)+"/r_rho_pv.dat")
                print("...")
            return r, P_C, P_K, [N_LI, N_PB, N_HE], RS, Δp
        ##################################################################################################################
        CURRENTDIR = os.getcwd()
        
        if RESTART or int(BLOCK0)==1:
            last = int(BLOCK0)
            if master:
                print("Trying to restart from '"+RSTFILE+"'... and running from block "+str(last)+" to "+str(int(BLOCKF)))
                for b in range(1, last-1):#i.e., from restart.1 to restart.{last-2}
                    print("=====>", b)
                    for p in range(nprocs):
                        r, P_C, P_K, [N_LI, N_PB, N_HE], RS, Δp = ANALITZA(b, p, N_LI, N_PB, N_HE, P_K, P_C, RS, Δp)
            MPI.COMM_WORLD.Barrier()
            try:
                # lmp.command("read_restart restart."+str(BLOCK0)+"_"+str(nprocs-1))
                if RESTART:
                    lmp.command("read_restart "+RSTFILE)
                    COMANDES.SET_REGIONS(lmp, x0, xf, y0, yf, z0, zf, RBOMB)
                    COMANDES.SET_NEIGH(lmp)
                    COMANDES.SET_THERMO(lmp, 0, ["step", "temp", "press", "vol", "pe", "ke", "etail"])
                    COMANDES.SET_VARIABLES(lmp)
                    COMANDES.SET_POTENTIAL(lmp, CDEAM=MODEL=="FRAILE")
                # if EQUIL: ADDITIONAL.EQUILIBRACIO_NVT(lmp, temperatura, Tdamp, Nevery, runs_NVT)
                if AFEGEIX_YUKAWA: ADDITIONAL.ADD_YUKAWA(lmp, N_Pb, N_Li, rcutoff, temperatura, Tdamp, runs_NVT)
                COMANDES.FIX_ENSEMBLE(lmp, temperatura, Tdamp, pressio, Pdamp, ensemble)
                COMANDES.GROUPS(lmp, [["Li", [1]], ["He", [2]], ["Pb", [3]], ["liquid", [1, 3]]])
                # printonmaster("PROC", me)
                # quit()
                if FUNCIONS.evalua(INTEGRATE):
                    # for proc in range(proc0, nprocs):
                    for proc in range(nprocs):
                        printonmaster("PROC")
                        print(proc, me)
                        if PRESSIONS_PARCIALS: ADDITIONAL.PRESSIONS_PARCIALS(lmp, Nevery, proc)
                        if DIFUSIONS_Li_He_Pb: ADDITIONAL.DIFUSIONS(lmp,
                                                                    Nevery,
                                                                    str(int(float(runs_NPT)/float(Nevery)/100)),
                                                                    str(float(runs_NPT)/100), proc)
                        COMANDES.SET_OUTPUT(lmp, proc, Nevery, str(int(float(runs_NPT)/float(Nevery))), runs_NPT)

                        lmp.command("run " + str(runs_NPT))
                        COMANDES.UNFIX_UNCOMPUTE(lmp)
                        if DIFUSIONS_Li_He_Pb: ADDITIONAL.UNDIF(lmp)
                        if PRESSIONS_PARCIALS: ADDITIONAL.UNPRES(lmp)
                    lmp.command("write_restart restart."+str(last))
                    try:
                        os.mkdir("TENSIO"+str(last)+"_"+str(me))
                    except FileExistsError as err:
                        pass
                    os.chdir("TENSIO"+str(last)+"_"+str(me))

                    shutil.copyfile("../NPT", "NPT")
                    shutil.move("../"+str(me)+"_npt.xyz", "xyz.out")
                    shutil.move("../"+str(me)+"_simulation_box.out", "simulation_box.out")
                    shutil.move("../"+str(me)+"_center_mass.out", "center_mass.out")
                    shutil.move("../"+str(me)+"_thermo_data.out", "thermo_data.out")
                    shutil.move("../"+str(me)+"_averages.out", "averages.out")
                    try:
                        shutil.move("../"+str(me)+"_partial_pressures.out", "partial_pressures.out")
                        shutil.move("../"+str(me)+"_dif.out", "dif.out")
                    except FileNotFoundError as err:
                        pass
                    shutil.move("../"+str(me)+"_aver.rdf", "aver.rdf")
                else:
                    os.chdir("TENSIO"+str(last)+"_"+str(me))

                FUNCIONS.CreateInputFile("st.params", int((me+nprocs*last)*DELTA),
                                         int((me+nprocs*last+1)*DELTA),
                                         int(Nevery),
                                         0.0050, 1, 500, MODEL, ".False.", ".True.")
                shutil.copyfile("../parent-directory", "parent-directory")
                subprocess.call(["../tensio.x"])
                N_He, LOST = FUNCIONS.lost_atoms(True) #En aquest cas s'haurà actualitzat NHe, pel que crearem de nou el fitxer st.params
                if LOST:
                    os.rename("LOST_ATOMS", "LOST_ATOMS_INI")
                    FUNCIONS.CreateInputFile("st.params",
                                             int((me+nprocs*last)*DELTA),
                                             int((me+nprocs*last+1)*DELTA),
                                             int(Nevery),
                                             0.0050, 1, 500, MODEL, ".False.", ".True.", N_He=N_He)
                    subprocess.call(["../tensio.x"])
                    NHe_NEW = N_He
                os.chdir("../")
            except:
                print("????????????????????????")
                last = int(BLOCK0)
                LOST = False
                os.chdir(CURRENTDIR)
    # for block in range(100):
        else:
            last = int(BLOCK0)
            LOST = False
        ##################################################################################################################
        MPI.COMM_WORLD.Barrier()
        for block in range(last+1, BLOCKF+1): #int(TOT/nprocs)
            MPI.COMM_WORLD.Barrier()
            for proc in range(nprocs):
                if PRESSIONS_PARCIALS: ADDITIONAL.PRESSIONS_PARCIALS(lmp, Nevery, proc)
                    # lmp.file("PRESSIONS_PARCIALS")
                if DIFUSIONS_Li_He_Pb: ADDITIONAL.DIFUSIONS(lmp, Nevery, str(int(float(runs_NPT)/float(Nevery)/100)), str(float(runs_NPT)/100), proc)
                    # lmp.file("DIFUSIONS_Li_He_Pb")
                COMANDES.SET_OUTPUT(lmp, proc, Nevery, str(int(float(runs_NPT)/float(Nevery))), runs_NPT)

                lmp.command("run " + str(runs_NPT))
                COMANDES.UNFIX_UNCOMPUTE(lmp)
                if DIFUSIONS_Li_He_Pb: ADDITIONAL.UNDIF(lmp)
                if PRESSIONS_PARCIALS: ADDITIONAL.UNPRES(lmp)
                lmp.command("write_restart restart."+str(block))
            
            try:
                os.mkdir("TENSIO"+str(block)+"_"+str(me))
            except FileExistsError as err:
                pass
            os.chdir("TENSIO"+str(block)+"_"+str(me))

            
            shutil.copyfile("../NPT", "NPT")
            shutil.move("../"+str(me)+"_npt.xyz", "xyz.out")
            shutil.move("../"+str(me)+"_simulation_box.out", "simulation_box.out")
            shutil.move("../"+str(me)+"_center_mass.out", "center_mass.out")
            shutil.move("../"+str(me)+"_thermo_data.out", "thermo_data.out")
            shutil.move("../"+str(me)+"_averages.out", "averages.out")
            try:
                shutil.move("../"+str(me)+"_partial_pressures.out", "partial_pressures.out")
                shutil.move("../"+str(me)+"_dif.out", "dif.out")
            except FileNotFoundError as err:
                pass
            shutil.move("../"+str(me)+"_aver.rdf", "aver.rdf")
            if LOST:
                FUNCIONS.CreateInputFile("st.params",
                                         int((me+nprocs*block)*DELTA),
                                         int((me+nprocs*block+1)*DELTA),
                                         int(Nevery),
                                         0.0050, 1, 500, MODEL, ".False.", ".True.", N_He=NHe_NEW)
            else:
                FUNCIONS.CreateInputFile("st.params",
                                         int((me+nprocs*block)*DELTA),
                                         int((me+nprocs*block+1)*DELTA),
                                         int(Nevery),
                                         0.0050, 1, 500, MODEL, ".False.", ".True.")
            MPI.COMM_WORLD.Barrier()
            shutil.copyfile("../parent-directory", "parent-directory")
            subprocess.call(["../tensio.x"])
                
            FUNCIONS.lost_atoms(False) #La simulació s'aturarà automàticament (saltarà un warning quan l'argument sigui False i hi hagi un error en el nombre d'àtoms)...
                # Només permetem que s'actualitzi el nombre d'àtoms si es perden al primer block. En cas que es tornin a perdre, llavors la simulació s'aturarà...
            os.remove("xyz.out")
            os.remove("simulation_box.out")
            os.remove("center_mass.out")
            os.remove("thermo_data.out")
            os.remove("averages.out")
            os.chdir("..")
            MPI.COMM_WORLD.Barrier()

            if master:
                for p in range(nprocs):
                    r, P_C, P_K, [N_LI, N_PB, N_HE], RS, Δp = ANALITZA(block, p, N_LI, N_PB, N_HE, P_K, P_C, RS, Δp)
                    print(len(r))
                    print(len(P_C))
                    print(len(P_K))
            MPI.COMM_WORLD.Barrier()
        lmp.commands_list(["unfix 1"])

        
        if master:
            N_Li = np.mean(N_LI, axis=0)
            N_Pb = np.mean(N_PB, axis=0)
            N_He = np.mean(N_HE, axis=0)
            err_P_K = np.std(P_K, axis=0)/np.sqrt(len(P_K))
            P_K = np.mean(P_K, axis=0)
            err_P_C = np.std(P_C, axis=0)/np.sqrt(len(P_C))
            P_C = np.mean(P_C, axis=0)
            print("==========================================================")
            print("[", np.mean(γS), "+-", np.sqrt(np.std(γS)/len(γS)), "] N / m")
            print("[", np.mean(RS), "+-", np.sqrt(np.std(RS)/len(RS)), "] Å")
            print("[", np.mean(Δp), "+-", np.sqrt(np.std(Δp)/len(Δp)), "] eV / Å³")

            def f(x,A,B,C,D):
                return 0.5*(A+B) - 0.5*(A-B) * np.tanh(2.*(x-C)/D)
            
            with open("r_rho_pv.dat", "w") as fileout:
                for k in range(len(r)):
                    fileout.write("{:.5E}".format(r[k])+"\t"
                                +"{:.5E}".format(N_Li[k])+"\t"
                                +"{:.5E}".format(N_Pb[k])+"\t"
                                +"{:.5E}".format(N_He[k])+"\t"
                                +"{:.5E}".format(P_K[k])+"\t"
                                +"{:.5E}".format(P_C[k])+"\n")
            os.system(PARENTDIR+"/INTERFACIAL_TENSION/spherical_integral.x > ST.aver")
            with open("ST.aver", "r") as filein:
                for line in filein:
                    if "Radius" in line:
                        print("Rs = {:.3f}".format(float(line.split()[0]))+" Å")
                        Rs = float(line.split()[0])
                    elif "eV Å^{-3}" in line and line.split()[0]!="NaN":
                        print("Δp = {:.5f}".format(-float(line.split()[0]))+" eV / Å³ =  {:.3f}".format(-float(line.split()[0])*1.602e-19*1e25)+" bar")
            print("==========================================================")

            if MOVING_AVG:
                r = MVA(r, 5)
                P_K = MVA(P_K, 5)
                P_C = MVA(P_C, 5)
                err_P_K = MVA(err_P_K, 5)
                err_P_C = MVA(err_P_C, 5)
            
            po = curve_fit(f, r, np.array(P_K+P_C), p0=[np.array(P_K+P_C)[int(len(r)*(5-r[0])/(r[-1]-r[0]))], np.array(P_K+P_C)[int(len(r)*(30-r[0])/(r[-1]-r[0]))], Rs, 1.])

        return np.mean(γS), np.sqrt(np.std(γS)/len(γS)), "N/m"


    def THERMO_INTEG(lmp,
                    #  temperatura: float, pressio: float, Tdamp: float, Pdamp: float,
                    #  rcutoff: float,
                    #  runs_NVT: int, runs_NPT: int, Nevery: int,
                    #  ensemble="NPT", RESTART=False, AFEGEIX_YUKAWA=False, PRESSIONS_PARCIALS=False, DIFUSIONS_Li_He_Pb=False,
                     TOT=nprocs):
        DELTA = float(runs_NPT) / float(Nevery)
        from src.CONFIGURACIO import SIGNE
        from potentials.POTENCIALS import UHeHe as V_22, ULiHe as V_12, UPbHe as V_23
        for block in range(int(TOT/nprocs)):
            if master:
                FILEOUT = open("TEST_AREA_"+str(block)+".out", "w")
                FILEOUT.write("########################\n")
                FILEOUT.close()
                ΔG = []
                ΔA = []
                Rd = [] #llegir a partir de la simulacio inicial sense potencial repulsiu... (i.e. lambda=-infty)
                Rd_He = []
                Rd_LM = []
                Rt = []
                Rd_rho = []
                Rd_eqm = []
                δ_Rd = []
                δ_Rd_He = []
                δ_Rd_LM = []
                δ_Rt = []
                δ_Rd_rho = []
                δ_Rd_eqm = []
                RMSD = np.zeros(nprocs)
                RMCD = np.zeros(nprocs)
                Rmin = np.zeros(nprocs)
                Rmax = np.zeros(nprocs)
                Rnew = np.zeros(nprocs)
                Reqm = np.zeros(nprocs)
                VBorn = []
                Λ = []

            for k, λ in enumerate(np.concatenate((np.linspace(0, SIGNE*0.3, 25), np.linspace(SIGNE*0.3, 0, 25)))):
                MPI.COMM_WORLD.Barrier()
                if master:
                    UB = np.zeros(nprocs)
                    RMSD = np.zeros(nprocs)
                    RMCD = np.zeros(nprocs)
                    # FUNCIONS.CreateBornTables(λ, 1.)
                    # FUNCIONS.CreateBornTables(0., (1.+λ))
                    FUNCIONS.CreateBornTables(0., 1., kappa=λ)
                MPI.COMM_WORLD.Barrier()
                for proc in range(nprocs):
                    # COMANDES.SET_POTENTIAL(lmp, Born=True, Lambda=λ)
                    COMANDES.SET_POTENTIAL(lmp, BornTable=True, CDEAM=MODEL=="FRAILE")
                    COMANDES.SET_OUTPUT(lmp, "LAM_"+str(k)+"_"+str(proc), Nevery, str(int(0.75*float(runs_NPT)/float(Nevery))), runs_NPT)
                    # ADDITIONAL.BORN(lmp, Nevery, str(int(0.75*float(runs_NPT)/float(Nevery))), runs_NPT, "LAM_"+str(k)+"_"+str(proc))
                    if PRESSIONS_PARCIALS: ADDITIONAL.PRESSIONS_PARCIALS(lmp, Nevery, "LAM_"+str(k)+"_"+str(proc))
                    if DIFUSIONS_Li_He_Pb: ADDITIONAL.DIFUSIONS(lmp, Nevery, str(int(0.75*float(runs_NPT)/float(Nevery))), runs_NPT, "LAM_"+str(k)+"_"+str(proc))
                    lmp.command("run " + str(runs_NPT))
                    COMANDES.UNFIX_UNCOMPUTE(lmp)
                    # ADDITIONAL.UNBORN(lmp)
                    if DIFUSIONS_Li_He_Pb: ADDITIONAL.UNDIF(lmp)
                    if PRESSIONS_PARCIALS: ADDITIONAL.UNPRES(lmp)
                    lmp.command("write_restart restart."+str(block)+"_LAM_"+str(k))

                try:
                    os.mkdir("LAM_"+str(k)+"_"+str(block)+"_"+str(me))
                except FileExistsError as err:
                    pass
                os.chdir("LAM_"+str(k)+"_"+str(block)+"_"+str(me))
                with open("COUPLING_PARAMETER", "w") as cpf:
                    cpf.write("V(r,λ) = λ * (12.483 eV) * exp(-r/(1 AA)) \n")
                    cpf.write("λ = "+str(λ)+"\n")
                    cpf.write("V(r,λ) = ("+str(SIGNE*12.483)+" eV) * exp(-r/(1 AA) + λ) \n")
                    cpf.write("λ = "+str(np.log(λ))+"\n")


                shutil.copyfile("../NPT", "NPT")
                shutil.move("../LAM_"+str(k)+"_"+str(me)+"_npt.xyz", "xyz.out")
                shutil.move("../LAM_"+str(k)+"_"+str(me)+"_simulation_box.out", "simulation_box.out")
                shutil.move("../LAM_"+str(k)+"_"+str(me)+"_center_mass.out", "center_mass.out")
                shutil.move("../LAM_"+str(k)+"_"+str(me)+"_thermo_data.out", "thermo_data.out")
                shutil.move("../LAM_"+str(k)+"_"+str(me)+"_averages.out", "averages.out")
                #shutil.move("../LAM_"+str(k))+"_"+str(me)+"_born_tts_s.out", "born_tts_s.out")
                try:
                    shutil.move("../LAM_"+str(k)+"_"+str(me)+"_partial_pressures.out", "partial_pressures.out")
                    shutil.move("../LAM_"+str(k)+"_"+str(me)+"_dif.out", "dif.out")
                except FileNotFoundError as err:
                    pass
                shutil.move("../LAM_"+str(k)+"_"+str(me)+"_aver.rdf", "aver.rdf")

                if block==0:
                    FUNCIONS.CreateInputFile("st.params", int((me+nprocs*block)*DELTA), int((me+nprocs*block+1)*DELTA), int(Nevery), 0.1, 1, 350, MODEL, ".False.", ".True.")
                    os.system("../radius.x > radius.out")
                    N_He, LOST = FUNCIONS.lost_atoms(True) #En aquest cas s'haurà actualitzat NHe, pel que crearem de nou el fitxer st.params
                    if LOST:
                        os.rename("LOST_ATOMS", "LOST_ATOMS_INI")
                        FUNCIONS.CreateInputFile("st.params", int((me+nprocs*block)*DELTA), int((me+nprocs*block+1)*DELTA), int(Nevery), 0.1, 1, 350, MODEL, ".False.", ".True.", N_He=N_He)
                        os.system("../radius.x > radius.out")
                        NHe_NEW = N_He
                else:
                    if LOST:
                        FUNCIONS.CreateInputFile("st.params", int((me+nprocs*block)*DELTA), int((me+nprocs*block+1)*DELTA), int(Nevery), 0.1, 1, 350, MODEL, ".False.", ".True.", N_He=NHe_NEW)
                    else:
                        FUNCIONS.CreateInputFile("st.params", int((me+nprocs*block)*DELTA), int((me+nprocs*block+1)*DELTA), int(Nevery), 0.1, 1, 350, MODEL, ".False.", ".True.")
                    os.system("../radius.x > radius.out")
                    
                    FUNCIONS.lost_atoms(False) #La simulació s'aturarà automàticament (saltarà un warning quan l'argument sigui False i hi hagi un error en el nombre d'àtoms)...
                    # Només permetem que s'actualitzi el nombre d'àtoms si es perden al primer block. En cas que es tornin a perdre, llavors la simulació s'aturarà...
                    
                MPI.COMM_WORLD.Barrier()
                critical_error = False
                if master:
                    Λ.append(λ)
                    for proc in range(nprocs):
                        try:
                            with open("../LAM_"+str(k)+"_"+str(block)+"_"+str(proc)+"/born_tts_s.out", "r") as bttss:
                                for line in bttss:
                                    if "#" not in line:
                                        UB[proc] = float(line.split()[2])
                        except FileNotFoundError as err:
                            try:
                                with open("../LAM_"+str(k)+"_"+str(block)+"_"+str(proc)+"/averages.out", "r") as avgs:
                                    for line in avgs:
                                        if "#" not in line:
                                            vol = float(line.split()[5])
                                with open("../LAM_"+str(k)+"_"+str(block)+"_"+str(proc)+"/aver.rdf", "r") as rdf:
                                    X = []
                                    Y = []
                                    for line in rdf:
                                        if "#" not in line:
                                            line = line.split()
                                            if len(line) == 14:
                                                r = float(line[1])
                                                g_12 = float(line[4])
                                                g_23 = float(line[10])
                                                g_22 = float(line[8])
                                                X.append(r)
                                                Y.append(r*r*(N_Li*g_12*V_12(r-λ*2.5*0.5291772)) + N_Pb*g_23*V_23(r-λ*2.5*0.5291772) + (N_He-1)*g_22*V_22(r-λ*2.5*0.5291772))
                                                # Y.append(r*r*np.exp(-r/1.)*(N_Li*g_12 + N_Pb*g_23 + (N_He-1)*g_22))
                                                # Y.append(r*r*np.exp(-r/1.)*(N_Li*g_12 + N_Pb*g_23))
                                    X = np.array(X)
                                    Y = np.array(Y)
                                    UB[proc] = np.trapz(Y, X) / vol
                            except FileNotFoundError as inner_error:
                                print("================")
                                print(f"Critical error: {inner_error}")
                                print("================")
                                critical_error = True
                                raise
                        with open("../LAM_"+str(k)+"_"+str(block)+"_"+str(proc)+"/radius.out", "r") as radius:
                            for line in radius:
                                if "AVERAGE MSR:" in line:
                                    RMSD[proc] = float(line.split()[2])
                                elif "AVERAGE MCR:" in line:
                                    RMCD[proc] = float(line.split()[2])
                                elif "LM RADIUS:" in line:
                                    Rmin[proc] = float(line.split()[2])
                                elif "He RADIUS:" in line:
                                    Rmax[proc] = float(line.split()[2])
                                elif "NEW RADIUS:" in line:
                                    Rnew[proc] = float(line.split()[2])
                                elif "EQUIM RADIUS:" in line:
                                    Reqm[proc] = float(line.split()[2])

                        with open("../LAM_"+str(k)+"_"+str(block)+"_"+str(proc)+"/TEST_AREA_NAME.out", "w") as TA_file:
                            TA_file.write("# BLOCK \tPROC \t CP \tRMSD (\AA) \tREQM \tFREE ENER (eV) \n")
                            TA_file.write(str(block)
                                        +"\t"+str(proc)
                                        +"\t{:.5f}".format(λ)
                                        +"\t{:.6f}".format(RMSD[proc])
                                        +"\t{:.6f}".format(Reqm[proc])
                                        +"\t{:.6f}".format(UB[proc])+"\n")
                    try:
                        VBorn.append(np.mean(UB))
                        # ΔG.append(np.trapz(VBorn, Λ)* 4*np.pi * N_He * 12.483)
                        ΔG.append(np.trapz(VBorn, Λ)* 4*np.pi * N_He)
                        Rd.append(np.mean(RMSD))
                        δ_Rd.append(np.std(RMSD) / np.sqrt(len(RMSD)-1))
                        Rt.append(np.mean(RMCD))
                        δ_Rt.append(np.std(RMCD) / np.sqrt(len(RMCD)-1))
                        Rd_He.append(np.mean(Rmax))
                        δ_Rd_He.append(np.std(Rmax) / np.sqrt(len(Rmax)-1))
                        Rd_LM.append(np.mean(Rmin))
                        δ_Rd_LM.append(np.std(Rmin) / np.sqrt(len(Rmin)-1))
                        Rd_rho.append(np.mean(Rnew))
                        δ_Rd_rho.append(np.std(Rnew) / np.sqrt(len(Rnew)-1))
                        Rd_eqm.append(np.mean(Reqm))
                        δ_Rd_eqm.append(np.std(Reqm) / np.sqrt(len(Reqm)-1))

                        ΔA.append(4.*np.pi*(Rd[k]**2. - Rd[0]**2.))
                        # FILEOUT = open("TEST_AREA_"+str(block)+".out", "a")
                        FILEOUT = open("../TEST_AREA_"+str(block)+".out", "a")
                        FILEOUT.write("{:.5f}".format(Λ[k])+
                                    "\t{:.4f}".format(Rd[k])+"\t{:.4f}".format(δ_Rd[k])+
                                    "\t{:.4f}".format(Rd_LM[k])+"\t{:.4f}".format(δ_Rd_LM[k])+
                                    "\t"+"{:.4f}".format(Rd_He[k])+"\t"+"{:.4f}".format(δ_Rd_He[k])+
                                    "\t"+"{:.4f}".format(Rd_rho[k])+"\t"+"{:.4f}".format(δ_Rd_rho[k])+
                                    "\t"+"{:.4f}".format(Rd_eqm[k])+"\t"+"{:.4f}".format(δ_Rd_eqm[k])+
                                    "\t"+"{:.4f}".format(ΔA[k])+"\t"+"{:.4f}".format(ΔG[k])+"\n")
                        FILEOUT.close()
                        print("{:.4f}".format(Λ[k])+"\t"+"{:.4f}".format(Rd[k])+"\t"+"{:.4f}".format(Rd_He[k])+"\t"+"{:.4f}".format(Rd_LM[k])+"\t"+"{:.4f}".format(ΔA[k])+"\t"+"{:.4f}".format(ΔG[k])+"\n")
                    except IndexError as err:
                        critical_error = True
                        print("================")
                        print(err)
                        print("================")

                error_on_master = MPI.COMM_WORLD.bcast(critical_error, root=0)
                if error_on_master:
                    sys.exit("Error on master...Exit")

                os.remove("xyz.out")
                os.remove("simulation_box.out")
                os.remove("center_mass.out")
                os.remove("thermo_data.out")
                os.remove("averages.out")
                os.chdir("..")
                MPI.COMM_WORLD.Barrier()
            if master: FILEOUT.close()

                    
        lmp.commands_list(["unfix 1"])
        return lmp

    def PLANAR_VIRIAL(lmp, TOT=25*nprocs):
        ###############################################################3
        def ANALITZA(block, proc):
            with open("simulation_box.out", "r") as sbox:
                for line in sbox:
                    if "#" not in line:
                        line = line.split()
                        Lx = float(line[2]) - float(line[1])
                        Ly = float(line[4]) - float(line[3])
                        Lz = float(line[6]) - float(line[5])
            with open("ave_stress.out", "r") as avs:
                z = []
                p = []
                Δp = []
                γ = 0.
                for line in avs:
                    if "#" not in line:
                        line = line.split()
                        if len(line)!=2:
                            z.append(Δz*float(line[0]))
                            p.append(float(line[4]))
                            Δp.append(float(line[5]))
                            γ += float(line[5])
                        else:
                            Δz = Lz/float(line[1])
            γ *= Δz / 2.
            print("==============================")
            print("       VIRIAL METHOD ("+str(proc)+")      ")
            print("       ·············       ")
            print(γ*1e-5, "N/m")
            print(np.trapz(Δp,z)/2.*1e-5, "N/m")
            print("==============================")
            return
        ###############################################################

        ###############################################################
        # if master and "tensio.x" not in os.listdir("./"):
        #     FUNCIONS.CreateExecutable(PARENTDIR+"/INTERFACIAL_TENSION/",
        #                               "planar.interftens.x",
        #                               "tensio.x")
        #####
        for block in range(int(TOT/nprocs)):
            for proc in range(nprocs):
                if PRESSIONS_PARCIALS: ADDITIONAL.PRESSIONS_PARCIALS(lmp, Nevery, proc)
                if DIFUSIONS_Li_He_Pb: ADDITIONAL.DIFUSIONS(lmp, Nevery, str(int(float(runs_NPT)/float(Nevery)/100)), str(float(runs_NPT)/100), proc)
                COMANDES.SET_OUTPUT(lmp, proc, Nevery, str(int(float(runs_NPT)/float(Nevery))), runs_NPT, stress=True, velocity=False, xyz_file=False)
                lmp.command("run " + str(runs_NPT))
                COMANDES.UNFIX_UNCOMPUTE(lmp, stress=True, velocity=False, xyz_file=False)
                if DIFUSIONS_Li_He_Pb: ADDITIONAL.UNDIF(lmp)
                if PRESSIONS_PARCIALS: ADDITIONAL.UNPRES(lmp)
                lmp.command("write_restart restart."+str(block))

            try:
                os.mkdir("TENSIO"+str(block)+"_"+str(me))
            except FileExistsError as err:
                pass
            os.chdir("TENSIO"+str(block)+"_"+str(me))
            
            shutil.copyfile("../NPT", "NPT")
            shutil.move("../"+str(me)+"_simulation_box.out", "simulation_box.out")
            shutil.move("../"+str(me)+"_center_mass.out", "center_mass.out")
            shutil.move("../"+str(me)+"_thermo_data.out", "thermo_data.out")
            shutil.move("../"+str(me)+"_averages.out", "averages.out")
            shutil.move("../ave_stress_"+str(me)+".out", "ave_stress.out")
            shutil.move("../pressure_tensor_"+str(me)+".out", "pressure_tensor.out")
            try:
                shutil.move("../"+str(me)+"_partial_pressures.out", "partial_pressures.out")
                shutil.move("../"+str(me)+"_dif.out", "dif.out")
            except FileNotFoundError as err:
                pass
            shutil.move("../"+str(me)+"_aver.rdf", "aver.rdf")
            os.chdir("../")
            if master:
                for proc in range(nprocs):
                    os.chdir("TENSIO"+str(block)+"_"+str(proc))
                    ANALITZA(block, proc)
                    os.chdir("../")
        #####
        return lmp

    def PLANAR_IK(lmp, TOT=25*nprocs):
        if master and "tensio.x" not in os.listdir("./"):
            FUNCIONS.CreateExecutable(PARENTDIR+"/INTERFACIAL_TENSION/",
                                      "planar.interftens.x",
                                      "tensio.x")
        
        #############
        for block in range(int(TOT/nprocs)):
            for proc in range(nprocs):
                if PRESSIONS_PARCIALS: ADDITIONAL.PRESSIONS_PARCIALS(lmp, Nevery, proc)
                if DIFUSIONS_Li_He_Pb: ADDITIONAL.DIFUSIONS(lmp, Nevery, str(int(float(runs_NPT)/float(Nevery)/100)), str(float(runs_NPT)/100), proc)
                COMANDES.SET_OUTPUT(lmp, proc, Nevery, str(int(float(runs_NPT)/float(Nevery))), runs_NPT, stress=True, velocity=True, xyz_file=True)
                lmp.command("run " + str(runs_NPT))
                COMANDES.UNFIX_UNCOMPUTE(lmp, stress=True, velocity=False, xyz_file=False)
                if DIFUSIONS_Li_He_Pb: ADDITIONAL.UNDIF(lmp)
                if PRESSIONS_PARCIALS: ADDITIONAL.UNPRES(lmp)
                lmp.command("write_restart restart."+str(block))

            try:
                os.mkdir("TENSIO"+str(block)+"_"+str(me))
            except FileExistsError as err:
                pass
            os.chdir("TENSIO"+str(block)+"_"+str(me))
            
            shutil.copyfile("../NPT", "NPT")
            shutil.move("../"+str(me)+"_npt.xyz", "xyz.out")
            shutil.move("../"+str(me)+"_simulation_box.out", "simulation_box.out")
            shutil.move("../"+str(me)+"_center_mass.out", "center_mass.out")
            shutil.move("../"+str(me)+"_thermo_data.out", "thermo_data.out")
            shutil.move("../"+str(me)+"_averages.out", "averages.out")
            shutil.move("../ave_stress_"+str(me)+".out", "ave_stress.out")
            shutil.move("../pressure_tensor_"+str(me)+".out", "pressure_tensor.out")
            shutil.move("../"+str(me)+"_velocities.out", "velocities.out")
            # shutil.move("../"+str(me)+"_stress_per_atom.out", "stress_per_atom.out")
            try:
                shutil.move("../"+str(me)+"_partial_pressures.out", "partial_pressures.out")
                shutil.move("../"+str(me)+"_dif.out", "dif.out")
            except FileNotFoundError as err:
                pass
            shutil.move("../"+str(me)+"_aver.rdf", "aver.rdf")

            FUNCIONS.CreateInputFile("st.params", int((me+nprocs*block)*DELTA), int((me+nprocs*block+1)*DELTA), int(Nevery), 0.0050, 1, 500, MODEL, ".False.", ".True.")
            shutil.copyfile("../parent-directory", "parent-directory")
            subprocess.call(["../tensio.x"])
            os.remove("xyz.out")
            os.remove("velocities.out")
            os.chdir("../")
        #############
        return
