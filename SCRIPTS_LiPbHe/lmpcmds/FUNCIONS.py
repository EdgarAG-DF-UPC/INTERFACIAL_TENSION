import os, subprocess, shutil, sys
from mpi4py import MPI # Note that without the mpi4py specific lines from test.py running the script with mpirun on P processors would lead to P independent simulations to run parallel, each with a single processor. Therefore, if you use the mpi4py lines and you see multiple LAMMPS single processor outputs, mpi4py is not working correctly.
from lammps import lammps, MPIAbortException
# sys.path.append('../../../../../SCRIPTS_LiPb')
with open("parent-directory", "r") as pdir:
    for line in pdir: PARENTDIR = line
sys.path.append(PARENTDIR+"/SCRIPTS_LiPbHe/src")
# import LATTICE, FUNCIONS

def CreateExecutable(directori, order, executable, recompile=False):
    current = os.getcwd()
    print(current)
    if order not in os.listdir(current) or recompile:
        if order not in os.listdir(directori) or recompile:
            # os.mkdir("tmpdir")
            print(os.listdir())
            shutil.copytree(directori, "tmpdir")
            os.chdir("tmpdir")#(directori)        
            subprocess.call(["make", "clean", order]) # order = "TENSIO.x"
            shutil.copy(order, current+"/"+executable)
            os.chdir(current)
            shutil.rmtree("tmpdir")
            if recompile: print("RECOMPILED :)")
        else:
            print("ALREADY COMPILED :)")
            shutil.copy(directori+"/"+order, executable)

def force_field(model):
    if model=="LJ" or model=="1":
        return 1
    elif model in ["2", "BELASHCHENKO", "BELASHCHENKO.LiPb", "Pb-Bela12_Li-Bela11"]:
        return 2
    elif model in ["3", "AWAD", "AWAD.LiPb", "Pb-Bela12_Li-Awad23"]:
        return 3
    elif model in ["4", "Pb-Bela12_Li-Bela11_LJ84_19-Yukawa"]:
        return 4
    elif model in ["5", "Pb-Zhou01_Li-Awad23", "Pb-Zhou01_Li-Awad23_Morse23"]:
        return 5
    elif model in ["6", "Pb-Zhou04_Li-Awad23", "Pb-Zhou04_Li-Awad23_Morse23"]:
        return 6
    elif model in ["7", "FRAILE", "FRAILE.LiPb", "Pb-Zhou01_Li-Bela11"]:
        return 7
    else:
        print("El model "+model+" no és vàlid...")
        quit()

def READ_NPT():
    with open("NPT", "r") as NPT:
        DADES = NPT.read().split()
        NLi, NHe, NPb = [int(DADES[k]) for k in range(3)]
        Ntot = NLi+NHe+NPb
        P, T = [float(DADES[k]) for k in range(3,5)]
    return [NLi, NHe, NPb], P, T

def CreateInputFile(nom, sni, snf, Δtsn, Δr, nr0, nr1, FF, WT, PBC, **kwargs):
    [NLi, NHe, NPb], P, T = READ_NPT()
    NHe = kwargs.get("N_He", NHe) #el valor s'actualitzarà únicament quan indiquem explícitament el valor de N_He
        
    Ntot = NLi+NHe+NPb
    
    with open(nom, "w") as filein:
        filein.write(str(Ntot)+"\t"+str(NHe)+"\t"+str(NLi)+"\t"+str(NPb)+"\n")
        filein.write("{:.2f}".format(T)+"\t"+"{:.2f}".format(P)+"\n")
        filein.write(str(sni*Δtsn)+"\t"+str(snf*Δtsn)+"\t"+str(Δtsn)+"\n")
        filein.write(str(Δr)+"\t"+str(nr0)+"\t"+str(nr1)+"\n")
        filein.write(str(force_field(FF))+"\n")
        filein.write(WT+"\t"+PBC+"\n")
        # filein.write(".False.\t.False.\n")
        filein.write(".False.\t.False.\n")


def CreateBornTables(λ,ξ,recompile=False,**kwargs):
    if "tts.born.x" not in os.listdir("./") or recompile: CreateExecutable("/users/edgar/BUBBLES_NEW/TAULES/TTS_SBMILL_BORN_TABLES", "tts.born.x", "tts.born.x", recompile=recompile)
    # with open("LiHe.Born.INPUT", "w") as filein:
    #     filein.write("He \tLi \t"+str(ξ)+"\n")
    # os.system("./tts.born.x < LiHe.Born.INPUT") #Això crearà la taula LiHe_Born.table
    print(os.getcwd())
    κ = kwargs.get("kappa", 0.)
    LiHe_Born_INPUT = "He \tLi \t"+str(λ)+" \t"+str(ξ)+" \t"+str(κ)
    os.system("echo '"+LiHe_Born_INPUT+"' | ./tts.born.x")
    HeHe_Born_INPUT = "He \tHe \t"+str(λ)+" \t"+str(ξ)+" \t"+str(κ)
    # HeHe_Born_INPUT = "He \tHe \t"+str(0.)+" \t"+str(ξ)
    os.system("echo '"+HeHe_Born_INPUT+"' | ./tts.born.x")
    if "sladek.born.x" not in os.listdir("./") or recompile: CreateExecutable("/users/edgar/BUBBLES_NEW/TAULES/TTS_SBMILL_BORN_TABLES", "sladek.born.x", "sladek.born.x", recompile=recompile)
    # with open("PbHe.Born.INPUT", "w") as filein:
    #     filein.write(str(ξ)+"\n")
    # os.system("./sladek.born.x < PbHe.Born.INPUT") #Això crearà la taula PbHe_Born.table
    PbHe_Born_INPUT = str(λ)+" \t"+str(ξ)+" \t"+str(κ)
    os.system("echo '"+PbHe_Born_INPUT+"' | ./sladek.born.x")
    

def evalua(string):
    if string in [True, "True", "YES", "SI"]:
        return True
    elif string in [False, "False", "NO", None]:
        return False
    else:
        raise RuntimeError(string, "must be either True, YES, SI (True) or False, NO (False)")

def lost_atoms(actualitza):
    MPI.COMM_WORLD.Barrier()
    [NLi, NHe, NPb], P, T = READ_NPT()
    try:
        LOST = True
        with open("LOST_ATOMS", "r") as filein:
            for line in filein:
                if "He" in line:
                    N_He = int(line.split()[0])
                if "Li" in line:
                    N_Li = int(line.split()[0])
                if "Pb" in line:
                    N_Pb = int(line.split()[0])
        if N_He!=NHe and N_Li==NLi and N_Pb==NPb:
            if actualitza:
                print("On process "+str(MPI.COMM_WORLD.Get_rank())+":")
                print("   WARNING! The number of He atoms in the simulation is "+str(N_He)+" instead of "+str(NHe)+".")
                print("   The simulation will continue with:")
                print("      >"+str(N_Li)+" Li atoms")
                print("      >"+str(N_Pb)+" Pb atoms")
                print("      >"+str(N_He)+" He atoms")
            else:
                raise Warning("The number of He atoms in the simulation is "+str(N_He)+" instead of "+str(NHe)+". The simulation will be stopped.")
        else:
            raise RuntimeError("The indicated total number of atoms does not coincide with the current number of atoms... The simulation will be stopped.")
    except FileNotFoundError as err:
        N_He = NHe
        LOST = False
    return N_He, LOST
