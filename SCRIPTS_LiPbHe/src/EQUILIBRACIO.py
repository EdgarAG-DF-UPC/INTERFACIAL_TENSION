# from mpi4py import MPI # Note that without the mpi4py specific lines from test.py running the script with mpirun on P processors would lead to P independent simulations to run parallel, each with a single processor. Therefore, if you use the mpi4py lines and you see multiple LAMMPS single processor outputs, mpi4py is not working correctly.
# from lammps import lammps, PyLammps
# lmp = lammps(name="mpi", cmdargs=["-pk","omp","4","-sf","omp"])
import sys
with open("parent-directory", "r") as pdir: PARENTDIR = pdir.read().replace('\n','').replace('"','')
# sys.path.append('../../../../../SCRIPTS_LiPb')
sys.path.append(PARENTDIR+"/SCRIPTS_LiPbHe/src")
# import LATTICE, FUNCIONS
# from LiPb.Hebubbles.config import LATSTYLE, DENSITAT, delta_t, ENSEMBLE, runs_NVT, runs_NPT, EQUILIBRATION_LENNARD_JONES, MINIMIZATION, THERMAL_SHAKE, Tshake, AFEGEIX_YUKAWA, Nevery, m_Li, m_Pb, m_He
from CONFIGURACIO import runs_NVT, runs_NPT, EQUILIBRATION_LENNARD_JONES, MINIMIZATION, THERMAL_SHAKE, Tshake, AFEGEIX_YUKAWA, temperatura, Tdamp, pressio, Pdamp, ensemble, N_Li, N_Pb, N_He, rcutoff, MODEL
import COMANDES, ADDITIONAL

# lmp = lammps(name="mpi", cmdargs=["-pk","omp","4","-sf","omp"])#,"-screen","none" .... -sf omp indica que es faci servir la versió omp (si existeix): "will automatically append “omp” to styles that support "
# nprocs = MPI.COMM_WORLD.Get_size()
# me = MPI.COMM_WORLD.Get_rank()
# master = me==0
# DELTA = float(runs_NPT) / float(Nevery)#1000

from lammps import lammps, PyLammps
lmp = lammps(name="mpi", cmdargs=["-pk","omp","4","-sf","omp"])

def EQUILIBRATION_BUBBLES():
    #### We start by restarting the original crystalline structure...
    print("Restarting from restart.ini...")
    lmp.command("read_restart restart.ini")

    #### Define the variable names inside LAMMPS:
    lmp.commands_list(["run			0",
                    "variable N equal step",
                    "variable pote equal pe",
                    "variable kine equal ke",
                    "variable Etail equal etail",
                    "variable Etotal equal etotal",
                    "variable E2 equal v_Etotal*v_Etotal",
                    "variable T equal temp",
                    "variable T2 equal v_T*v_T",
                    "variable Press equal press",
                    "variable P2 equal v_Press*v_Press",
                    "variable V equal vol",
                    "variable V2 equal v_V*v_V",
                    "variable PV equal v_Press*v_V",
                    "variable P2V2 equal v_Press*v_Press*v_V*v_V",
                    "variable XXl equal xlo",
                    "variable XXh equal xhi",
                    "variable YYl equal ylo",
                    "variable YYh equal yhi",
                    "variable ZZl equal zlo",
                    "variable ZZh equal zhi",
                    "variable H equal enthalpy",
                    "variable H2 equal v_H*v_H"]) # REQUEST OUTPUT
    lmp.command("reset_timestep 0")

    COMANDES.SET_POTENTIAL(lmp, CDEAM=MODEL=="FRAILE") # SET EAM - TTS PARAMETERS

    #### We properly start the equilibration...    
    lmp.command(COMANDES.LMP_PRINT("EQUILIBRATION"))

    #### 'EQUILIBRATION_LENNARD_JONES':
    #### · Boolean variable, which will determine if the initial steps of the equilibration 
    ####   will be performed using either a Lennard-Jones potential or the EAM+TTS+S force field.
    #### · Setting 'EQUILIBRATION_LENNARD_JONES=True' is useful in order to accelerate the nucleation rate 
    #### (seggregation of helium will occur much faster),
    ####   which may be used when one does not care about the accuracy of the bubble formation process 
    ####   but one needs to generate a configuration in which the helium cluster (bubble) is already formed.
    if EQUILIBRATION_LENNARD_JONES: ADDITIONAL.LENNARD_JONES(lmp, rcutoff)
    #### 'MINIMIZATION':
    #### · Boolean variable.
    #### · If true, it will apply the Conjugated Gradient (CG) algorithm to modify all atomic coordinates
    ####   in order to minimize the energy.
    if MINIMIZATION: ADDITIONAL.MINIMIZATION(lmp)
    
    if EQUILIBRATION_LENNARD_JONES:
        for iter in range(10):
            lmp.commands_list(["dump 1 all xyz "+str(runs_NVT)+" xyz_LJ"+str(iter)+".out",
                            COMANDES.FIX_THERMOSTAT("NVT", temperatura, Tdamp),
                            "run "+str(runs_NVT),
                            "unfix 1",
                            "undump 1"])
        lmp.command("write_restart restart.LJ")
        COMANDES.SET_POTENTIAL(lmp, CDEAM=MODEL=="FRAILE") # RESET TO (A)EAM+TTS
    else:
        lmp.commands_list([COMANDES.FIX_THERMOSTAT("NVT", temperatura, Tdamp),
                        "run "+str(runs_NVT),
                        "unfix 1"])
        
    # if EQUILIBRATION_LENNARD_JONES: COMANDES.SET_POTENTIAL(lmp, CDEAM=MODEL=="FRAILE") # RESET TO (A)EAM+TTS

    #### 'THERMAL_SHAKE':
    #### · Boolean variable.
    #### · If true, it will warm the system up to a temperature 'Tshake', useful to melt the system.
    #### · 'Tshake' should be large enough to mix up all atoms, melting the solid, and avoiding the system gets trapped in an artificial solid.
    #### · After that, the system is cooled down to the temperature of the thermostat. 
    if THERMAL_SHAKE: ADDITIONAL.THERMAL_SHAKE(lmp, temperatura, Tshake, Tdamp, runs_NVT)

    lmp.commands_list(["dump 1 all xyz "+str(runs_NVT)+" xyz_equil_LiHe"+str(iter)+".out",
                    COMANDES.FIX_THERMOSTAT("NVT", temperatura, Tdamp),
                    "run " + str(runs_NVT),
                    "unfix 1",
                    "undump 1"])

    #### 'AFEGEIX_YUKAWA':
    #### · Boolean variable.
    #### · If true, additional Coulomb-Yukawa interactions are included for all pairs.
    if AFEGEIX_YUKAWA:
        x_Pb = 100. * float(N_Pb) / float(N_Li+N_Pb)
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
            raise RuntimeError("X must be 0, 10, 20, 30, 40, 50 or 100...")

        a_LiLi = "{:.8f}".format(qLi*qLi*14.399)
        a_PbPb = "{:.8f}".format(qPb*qPb*14.399)
        a_LiPb = "{:.8f}".format(qLi*qPb*14.399)
        lmp.commands_list([COMANDES.LMP_PRINT("ADDITION OF THE YUKAWA POTENTIAL"),
                        "pair_style	hybrid/overlay  eam/alloy table linear 750 table linear 750 table linear 750 yukawa 1.10 " + str(rcutoff) + " yukawa 1.10 " + str(rcutoff) + " yukawa 1.10 " + str(rcutoff),
                        "pair_coeff * * eam/alloy LiPb.eam.alloy Li NULL Pb",
                        "pair_coeff 1 2 table 1 ../../../../TAULES/LiHe.table SHENGTANG_LiHe",
                        "pair_coeff 2 2 table 2 ../../../../TAULES/HeHe.table SHENGTANG_HeHe",
                        "pair_coeff 2 3 table 3 ../../../../TAULES/PbHe.table SLADEK_PbHe",
                        "pair_coeff	1 1 yukawa 1 " + str(a_LiLi),
                        "pair_coeff	1 3 yukawa 2 " + str(a_LiPb),
                        "pair_coeff	3 3 yukawa 3 " + str(a_PbPb),
                        COMANDES.FIX_THERMOSTAT("NVT", temperatura, Tdamp),
                        "run " + str(runs_NVT),
                        "unfix 1"]) # MÒDUL DE POTENCIAL YUKAWA - CORRECCIÓ BELASHCHENKO 2019
        lmp.commands_list([COMANDES.FIX_THERMOSTAT("NVT", temperatura, Tdamp),
                        "run " + str(runs_NVT),
                        "unfix 1"])

    #### Reminder: we initialized the system with only Li+He atoms.
    #### · Now, we switch the identitity of 'N_Pb' Li atoms, one every 'DN = (N_Li + N_Pb) / N_Pb' Li atoms, to Pb.
    #### · After that, the system is equilibrated.
    #### · We end up with a 'restart.equil' file, which may be used to restart a simulation in which from an equilibrated configuration.
    if N_Pb>0:
        lmp.command(COMANDES.LMP_PRINT("Li --> Pb TRANSFORMATION"))
        DN = (N_Li + N_Pb) / N_Pb 
        # for i in range(1,N_Pb+1):
        #     b = int(DN*i)
        for i in range(N_Pb):
            # b = int(DN*(i+1))
            lmp.command("set atom " + str(int(DN*(i+1))) + " type 3")
        lmp.commands_list([COMANDES.FIX_THERMOSTAT("NVT", temperatura, Tdamp),
                        "run " + str(runs_NVT),
                        "unfix 1"])
        lmp.commands_list(["fix 1 all " + str(ensemble) + " temp " + str(temperatura) + " " + str(temperatura) + " " + str(Tdamp)+ " iso " + str(pressio) + " " + str(pressio) + " " + str(Pdamp),
                        "dump 1 all xyz "+str(runs_NVT)+" xyz_equil_PbLiHe"+str(iter)+".out",
                        "run " + str(runs_NVT),
                        "run " + str(runs_NVT),
                        "run " + str(runs_NVT),
                        "run " + str(runs_NVT),
                        "unfix 1",
                        "undump 1"])
        
    lmp.command("write_restart restart.equil")

    
    #### At this point, we are ready for production runs...    
    lmp.commands_list([COMANDES.LMP_PRINT(""),
                    COMANDES.LMP_PRINT(""),
                    "reset_timestep 0"])

    lmp.command(COMANDES.LMP_PRINT("PRODUCTION RUNS"))
    lmp.command("fix 1 all " + str(ensemble) + " temp " + str(temperatura) + " " + str(temperatura) + " " + str(Tdamp)+ " iso " + str(pressio) + " " + str(pressio) + " " + str(Pdamp))
    lmp.commands_list(["group   Li		type 1",
                    "group   He		type 2",
                    "group   Pb		type 3",
                    "group	liquid	type 1 3"])
    
    return lmp

def EQUILIBRATION_PLANAR():
    print("Restarting from restart.ini...")
    lmp.command("read_restart restart.ini") 
    lmp.commands_list(["run			0",
                    "variable N equal step",
                    "variable pote equal pe",
                    "variable kine equal ke",
                    "variable Etail equal etail",
                    "variable Etotal equal etotal",
                    "variable E2 equal v_Etotal*v_Etotal",
                    "variable T equal temp",
                    "variable T2 equal v_T*v_T",
                    "variable Press equal press",
                    "variable P2 equal v_Press*v_Press",
                    "variable V equal vol",
                    "variable V2 equal v_V*v_V",
                    "variable PV equal v_Press*v_V",
                    "variable P2V2 equal v_Press*v_Press*v_V*v_V",
                    "variable XXl equal xlo",
                    "variable XXh equal xhi",
                    "variable YYl equal ylo",
                    "variable YYh equal yhi",
                    "variable ZZl equal zlo",
                    "variable ZZh equal zhi",
                    "variable H equal enthalpy",
                    "variable H2 equal v_H*v_H"]) # REQUEST OUTPUT
    lmp.command("reset_timestep 0")
    
    COMANDES.SET_POTENTIAL(lmp, CDEAM=MODEL=="FRAILE") # SET EAM - TTS PARAMETERS

    lmp.command(COMANDES.LMP_PRINT("EQUILIBRATION"))

    if EQUILIBRATION_LENNARD_JONES: ADDITIONAL.LENNARD_JONES(lmp, rcutoff)

    if MINIMIZATION: ADDITIONAL.MINIMIZATION(lmp)

    if EQUILIBRATION_LENNARD_JONES:
        for iter in range(10):
            lmp.commands_list(["dump 1 all xyz "+str(runs_NVT)+" xyz_LJ"+str(iter)+".out",
                            COMANDES.FIX_THERMOSTAT("NVT", temperatura, Tdamp),
                            "run "+str(runs_NVT),
                            "unfix 1",
                            "undump 1"])
        lmp.command("write_restart restart.LJ")
    else:
        lmp.commands_list([COMANDES.FIX_THERMOSTAT("NVT", temperatura, Tdamp),
                        "run "+str(runs_NVT),
                        "unfix 1"])
        
    if EQUILIBRATION_LENNARD_JONES: COMANDES.SET_POTENTIAL(lmp, CDEAM=MODEL=="FRAILE") # RESET TO (A)EAM+TTS
        
    if THERMAL_SHAKE: ADDITIONAL.THERMAL_SHAKE(lmp, temperatura, Tshake, Tdamp, runs_NVT)

    lmp.commands_list(["dump 1 all xyz "+str(runs_NVT)+" xyz_equil_LiHe"+str(iter)+".out",
                    COMANDES.FIX_THERMOSTAT("NVT", temperatura, Tdamp),
                    "run " + str(runs_NVT),
                    "unfix 1",
                    "undump 1"])

    if AFEGEIX_YUKAWA:
        x_Pb = 100. * float(N_Pb) / float(N_Li+N_Pb)
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
            raise RuntimeError("X must be 0, 10, 20, 30, 40, 50 or 100...")

        a_LiLi = "{:.8f}".format(qLi*qLi*14.399)
        a_PbPb = "{:.8f}".format(qPb*qPb*14.399)
        a_LiPb = "{:.8f}".format(qLi*qPb*14.399)
        lmp.commands_list([COMANDES.LMP_PRINT("ADDITION OF THE YUKAWA POTENTIAL"),
                        "pair_style	hybrid/overlay  eam/alloy table linear 750 table linear 750 table linear 750 yukawa 1.10 " + str(rcutoff) + " yukawa 1.10 " + str(rcutoff) + " yukawa 1.10 " + str(rcutoff),
                        "pair_coeff * * eam/alloy LiPb.eam.alloy Li NULL Pb",
                        "pair_coeff 1 2 table 1 ../../../../TAULES/LiHe.table SHENGTANG_LiHe",
                        "pair_coeff 2 2 table 2 ../../../../TAULES/HeHe.table SHENGTANG_HeHe",
                        "pair_coeff 2 3 table 3 ../../../../TAULES/PbHe.table SLADEK_PbHe",
                        "pair_coeff	1 1 yukawa 1 " + str(a_LiLi),
                        "pair_coeff	1 3 yukawa 2 " + str(a_LiPb),
                        "pair_coeff	3 3 yukawa 3 " + str(a_PbPb),
                        COMANDES.FIX_THERMOSTAT("NVT", temperatura, Tdamp),
                        "run " + str(runs_NVT),
                        "unfix 1"]) # MÒDUL DE POTENCIAL YUKAWA - CORRECCIÓ BELASHCHENKO 2019
        lmp.commands_list([COMANDES.FIX_THERMOSTAT("NVT", temperatura, Tdamp),
                        "run " + str(runs_NVT),
                        "unfix 1"])


    if N_Pb>0:
        lmp.command(COMANDES.LMP_PRINT("Li --> Pb TRANSFORMATION"))
        DN = (N_Li + N_Pb) / N_Pb 
        for i in range(1,N_Pb+1):
            b = int(DN*i)
            lmp.command("set atom " + str(b) + " type 3")
        lmp.commands_list([COMANDES.FIX_THERMOSTAT("NVT", temperatura, Tdamp),
                        "run " + str(runs_NVT),
                        "unfix 1"])
        lmp.commands_list(["fix 1 all " + str(ensemble) + " temp " + str(temperatura) + " " + str(temperatura) + " " + str(Tdamp)+ " iso " + str(pressio) + " " + str(pressio) + " " + str(Pdamp),
                        "dump 1 all xyz "+str(runs_NVT)+" xyz_equil_PbLiHe"+str(iter)+".out",
                        "run " + str(runs_NVT),
                        "run " + str(runs_NVT),
                        "run " + str(runs_NVT),
                        "run " + str(runs_NVT),
                        "unfix 1",
                        "undump 1"])
        
    lmp.command("write_restart restart.equil")

    lmp.commands_list([COMANDES.LMP_PRINT(""),
                    COMANDES.LMP_PRINT(""),
                    "reset_timestep 0"])


    lmp.command(COMANDES.LMP_PRINT("THE SYSTEM HAS BEEN SUCCESSFULLY EQUILIBRATED :)"))
    lmp.command("fix 1 all " + str(ensemble) + " temp " + str(temperatura) + " " + str(temperatura) + " " + str(Tdamp)+ " iso " + str(pressio) + " " + str(pressio) + " " + str(Pdamp))
    lmp.commands_list(["group   Li		type 1",
                    "group   He		type 2",
                    "group   Pb		type 3",
                    "group	liquid	type 1 3"])
    return lmp

