from mpi4py import MPI
nprocs = MPI.COMM_WORLD.Get_size()
me = MPI.COMM_WORLD.Get_rank()
master = me==0

import os, shutil, sys
with open("parent-directory", "r") as pdir: PARENTDIR = pdir.read().replace('\n','').replace('"','')
# sys.path.append('../../../../../SCRIPTS_LiPb')
sys.path.append(PARENTDIR+"/SCRIPTS_LiPbHe/src")
import lmpcmds.COMANDES as COMANDES, lmpcmds.ADDITIONAL as ADDITIONAL, lmpcmds.FUNCIONS as FUNCIONS

class INITIALIZE:
    # This class initializes the system, prior to production runs.
    nprocs = MPI.COMM_WORLD.Get_size()
    me = MPI.COMM_WORLD.Get_rank()
    master = me==0
    def BUBBLES(N_Li: float, N_Pb: float, N_He: float,
                temperatura: float, pressio: float, Tdamp: float, Pdamp: float,
                rcutoff: float,
                runs_NVT: int, runs_NPT: int, Nevery: int,
                RBOMB: float, boxlims: list,
                ensemble="NPT", RESTART=False, RSTFILE="restart.equil", AFEGEIX_YUKAWA=False, PRESSIONS_PARCIALS=False, DIFUSIONS_Li_He_Pb=False):
        # Initialization to prepare systems in which a helium (spherical) bubble
        # is formed inside the liquid metal (pure or alloy).
        x0, xf, y0, yf, z0, zf = boxlims
        from lammps import lammps
        lmp = lammps(name="mpi", cmdargs=["-pk","omp","4","-sf","omp"])
        if not RESTART:
            # When we choose (in 'CONFIGURACIO') not to 'RESTART' the file from an existing simulation 
            # (which may have been stopped during production runs)
            try: # Try running an existing 'restart.equil' file 
                # (just check if in the current directory there is an equilibrated sample, then production runs will start immeadiately)...
                FITXERS = os.listdir("./")
                if RSTFILE not in FITXERS:
                    print("OOOOOOJOOOOO")
                    raise FileNotFoundError(RSTFILE+" is missing...")
                for p in range(nprocs):
                    if ("EQUILIBRATION_"+str(p)) not in FITXERS:
                        raise FileNotFoundError("EQUILIBRATION_"+str(p)+" is missing...")
                
                lmp.command("read_restart "+RSTFILE)
                COMANDES.SET_REGIONS(lmp, x0, xf, y0, yf, z0, zf, RBOMB)
                COMANDES.SET_NEIGH(lmp)
                COMANDES.SET_THERMO(lmp, 0, ["step", "temp", "press", "vol", "pe", "ke", "etail"])
                COMANDES.SET_VARIABLES(lmp)
                from CONFIGURACIO import MODEL
                COMANDES.SET_POTENTIAL(lmp, CDEAM=MODEL=="FRAILE")
                if AFEGEIX_YUKAWA: ADDITIONAL.ADD_YUKAWA(lmp, N_Pb, N_Li, rcutoff, temperatura, Tdamp, runs_NVT)
                COMANDES.FIX_ENSEMBLE(lmp, temperatura, Tdamp, pressio, Pdamp, ensemble)
                COMANDES.GROUPS(lmp, [["Li", [1]], ["He", [2]], ["Pb", [3]], ["liquid", [1, 3]]])
            except FileNotFoundError as err: # In case the 'restart' file does not exist, then the system will be initialized from the beginning (i.e., melting a crystalline structure, etc.).
                from src.EQUILIBRACIO import EQUILIBRATION_BUBBLES as EQUILIBRATION#, lmp#, master, nprocs, BLOCKS, PRESSIONS_PARCIALS, DIFUSIONS_Li_He_Pb, runs_NPT, DELTA, MODEL, Nevery, SIGNE, N_Li, N_Pb, N_He, rcutoff
                lmp = EQUILIBRATION()
                COMANDES.SET_NEIGH(lmp)
                for proc in range(nprocs):
                    COMANDES.SET_OUTPUT(lmp, proc, Nevery, str(max([1, int(0.75*float(runs_NPT)/float(Nevery))])), runs_NPT)
                    if PRESSIONS_PARCIALS: ADDITIONAL.PRESSIONS_PARCIALS(lmp, Nevery, proc)
                    if DIFUSIONS_Li_He_Pb: ADDITIONAL.DIFUSIONS(lmp,
                                                                Nevery,
                                                                str(int(float(runs_NPT)/float(Nevery)/100)),
                                                                str(float(runs_NPT)/100),
                                                                proc)
                    lmp.command("run " + str(runs_NPT))
                    COMANDES.UNFIX_UNCOMPUTE(lmp)
                    if DIFUSIONS_Li_He_Pb: ADDITIONAL.UNDIF(lmp)
                    if PRESSIONS_PARCIALS: ADDITIONAL.UNPRES(lmp)
                    lmp.command("write_restart restart.0")
                try:
                    os.mkdir("EQUILIBRATION_"+str(me))
                except FileExistsError as err:
                    print("Directory 'EQUILIBRATION_"+str(me)+"' already exists")
                os.chdir("EQUILIBRATION_"+str(me))
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
                os.chdir("../")
        return lmp, RESTART, RSTFILE
    
    def PLANAR(N_Li: float, N_Pb: float, N_He: float,
               temperatura: float, pressio: float, Tdamp: float, Pdamp: float,
               rcutoff: float,
               runs_NVT: int, runs_NPT: int, Nevery: int,
               RBOMB: float, boxlims: list,
               ensemble="NPT", RESTART=False, RSTFILE="restart.0", AFEGEIX_YUKAWA=False, PRESSIONS_PARCIALS=False, DIFUSIONS_Li_He_Pb=False):
        ######
        from lammps import lammps, MPIAbortException
        lmp = lammps(name="mpi", cmdargs=["-pk","omp","4","-sf","omp"])
        if "restart.equil" in os.listdir("./") or "restart.equil.1000K" in os.listdir("./"):
            print(os.listdir("./"))
            from CONFIGURACIO import runs_NVT, runs_NPT, AFEGEIX_YUKAWA, Nevery, x0, y0, z0, xf, yf, zf, temperatura, Tdamp, pressio, Pdamp, ensemble, N_Li, N_Pb, N_He, rcutoff, MODEL

            nprocs = MPI.COMM_WORLD.Get_size()
            me = MPI.COMM_WORLD.Get_rank()
            DELTA = float(runs_NPT) / float(Nevery)#1000
            master = me==0

            if "restart.equil.1000K" in os.listdir("./"):
                lmp.command("read_restart restart.equil.1000K")
            elif "restart.equil" in os.listdir("./"):    
                lmp.command("read_restart restart.equil")
            else:
                raise FileNotFoundError("Falten els fitxers 'restart.equil' i 'restart.equil.1000K'...")

            lmp.commands_list(["region whole block "+str(x0)+" "+str(xf)+" "+str(y0)+" "+str(yf)+" "+str(z0)+" "+str(2.*zf-z0),
                               "region liquid block "+str(x0)+" "+str(xf)+" "+str(y0)+" "+str(yf)+" "+str(z0)+" "+str(zf),
                               "region gas block "+str(x0)+" "+str(xf)+" "+str(y0)+" "+str(yf)+" "+str(zf)+" "+str(2.*zf-z0)])
            COMANDES.SET_NEIGH(lmp)
            COMANDES.SET_THERMO(lmp, 0, ["step", "temp", "press", "vol", "pe", "ke", "etail"])
            COMANDES.SET_VARIABLES(lmp)
            if master: FUNCIONS.CreateBornTables(0., 1.)
            COMANDES.SET_POTENTIAL(lmp, BornTable=True, CDEAM=MODEL=="FRAILE")
            # if EQUIL: ADDITIONAL.EQUILIBRACIO_NVT(lmp, temperatura, Tdamp, Nevery, runs_NVT)
            if AFEGEIX_YUKAWA: ADDITIONAL.ADD_YUKAWA(lmp, N_Pb, N_Li, rcutoff, temperatura, Tdamp, runs_NVT)
            COMANDES.FIX_ENSEMBLE(lmp, temperatura, Tdamp, pressio, Pdamp, ensemble)
            COMANDES.GROUPS(lmp, [["Li", [1]], ["He", [2]], ["Pb", [3]], ["liquid", [1, 3]]])
            EQUIL = False
            if "restart.equil.1000K" in os.listdir("./"):
                lmp.commands_list(["fix 1 all " + str(ensemble) + " temp " + str(1000.00) + " " + str(temperatura) + " " + str(Tdamp)+ " iso " + str(pressio) + " " + str(pressio) + " " + str(Pdamp),
                                "run " + str(runs_NVT),
                                "unfix 1"])
                lmp.commands_list([COMANDES.FIX_THERMOSTAT("NVT", temperatura, Tdamp),
                                "run " + str(runs_NVT),
                                "unfix 1"])
        else:
            print(":)")
            from src.EQUILIBRACIO import EQUILIBRATION_PLANAR as EQUILIBRATION, runs_NPT #lmp,
            COMANDES.SET_NEIGH(lmp)
            lmp = EQUILIBRATION()
        ######
        return lmp, RESTART, RSTFILE