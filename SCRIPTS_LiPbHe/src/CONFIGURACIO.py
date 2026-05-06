import os, subprocess, shutil, sys
from mpi4py import MPI # Note that without the mpi4py specific lines from test.py running the script with mpirun on P processors would lead to P independent simulations to run parallel, each with a single processor. Therefore, if you use the mpi4py lines and you see multiple LAMMPS single processor outputs, mpi4py is not working correctly.
with open("parent-directory", "r") as pdir: PARENTDIR = pdir.read().replace('\n','').replace('"','')
# sys.path.append('../../../../../SCRIPTS_LiPb')
sys.path.append(PARENTDIR+"/SCRIPTS_LiPbHe/src")
import LATTICE
sys.path.append(PARENTDIR+"/SCRIPTS_LiPbHe/lmpcmds")
import FUNCIONS

with open("RUN.INPUT", "r") as infile:
    INI_CONF = {}
    for line in infile:
        line = line.split(" = ")
        INI_CONF[line[0]] = line[1].replace("\n","")
    
    # PARENTDIR = INI_CONF.get("PARENT", "/users/edgar/BUBBLES_NEW/")
    try:
        os.listdir(PARENTDIR+"/INTERFACIAL_TENSION")
        os.listdir(PARENTDIR+"/SCRIPTS_LiPbHe/src")
        os.listdir(PARENTDIR+"/SIMULATIONS")
        os.listdir(PARENTDIR+"/TAULES")
    except FileNotFoundError as err:
        print(err)
        raise RuntimeError("The indicated parent directory '"+PARENTDIR+"' may be wrong...")
    
    
    runs_NVT = INI_CONF.get("RUNS_NVT")
    runs_NPT = INI_CONF.get("RUNS_NPT")
    delta_t = INI_CONF.get("TIMESTEP")
    Nevery = INI_CONF.get("NEVERY")
    ENSEMBLE = INI_CONF.get("ENSEMBLE")
    MODEL = INI_CONF.get("MODEL")
    AFEGEIX_YUKAWA = FUNCIONS.evalua(INI_CONF.get("YUKAWA_CORRECTION", "NO"))
    DENSITAT = INI_CONF.get("DENSINI")
#     DENSITAT = infile[6+a]
    LATSTYLE = INI_CONF.get("LATTICE")
    EQUILIBRATION_LENNARD_JONES = FUNCIONS.evalua(INI_CONF.get("LJ", "NO"))
    PRESSIONS_PARCIALS = FUNCIONS.evalua(INI_CONF.get("PP", "NO"))
    ENERGIES_EXTRA = FUNCIONS.evalua(INI_CONF.get("EE", "NO"))
    DIFUSIONS_Li_He_Pb = FUNCIONS.evalua(INI_CONF.get("DD", "NO"))
    THERMAL_SHAKE = FUNCIONS.evalua(INI_CONF.get("SHAKE", "NO"))
    Tshake = INI_CONF.get("Temp_Shake", None)
    NSshake = INI_CONF.get("Nsteps_Shake", None)
    MINIMIZATION = INI_CONF.get("MIN", False)
    CUA = INI_CONF.get("QUEUE", "all")
    BLOCKS = INI_CONF.get("NBLOCKS", 25)
    BLOCK0 = INI_CONF.get("BLOCK0", 1)
    BLOCKF = INI_CONF.get("BLOCKF", BLOCKS)
    INTEGRATE = INI_CONF.get("INTEGRATE", True)
    proc0 = INI_CONF.get("LAST_PROC", 0)
    SIGNE = INI_CONF.get("SIGN", "+") # +1 ==> repulsive ==> forward; -1 ==> attractive ==> backward
    if SIGNE=="+":
        SIGNE = +1
    elif SIGNE=="-":
        SIGNE = -1
    else:
        raise ValueError("Variable SIGN must be either '+' or '-'")
    RESTART = INI_CONF.get("RESTART", False)
    if FUNCIONS.evalua(RESTART):
        last = 0
        for file in os.listdir("./"):
            if ("restart" in file) and ("equil" not in file) and ("ini" not in file) and ("LJ" not in file):
                print(file)
                iter = int(file.replace("restart.", ""))
                if iter>last: last = iter
        RSTFILE = "restart."+str(last-2)
        try:
            FILES = os.listdir(".")
            if RSTFILE not in FILES:
                raise FileNotFoundError("File '"+RSTFILE+"' does not exist...")
            else:
                BLOCK0 = last-1
        except FileNotFoundError as err:
            try:
                if "restart.equil" not in FILES:
                    raise FileNotFoundError("File 'restart.equil' neither exists...")
                else:
                    RSTFILE = "restart.equil"
            except:
                RSTFILE = "restart.ini"
            print("Using '"+RSTFILE+"', instead.")
            BLOCK0 = 1
    else:
        RSTFILE = INI_CONF.get("RSTFILE", "restart.equil")
        BLOCK0 = 1
    METHOD = INI_CONF.get("METHOD", None)
    MOVING_AVG = INI_CONF.get("MOVING_AVG", False)


m_Li = 6.941 # u.m.a. = g/mol
m_He = 4.002 # u.m.a. = g/mol
m_Pb = 207.20 # u.m.a. = g/mol
[N_Li, N_He, N_Pb], pressio, temperatura = FUNCIONS.READ_NPT()

C, ALAT, DN, NLi, NPb, Tmin = LATTICE.LATTICE(LATSTYLE, N_Li+N_Pb, float(N_Li)/float(N_Li+N_Pb), temperatura, DENSITAT, MPI.COMM_WORLD.Get_rank()==0)
x0 = -C/2
xf = C/2
y0 = -C/2
yf = C/2
z0 = -C/2
zf = C/2
# RBOMB = C*ALAT/6.
RBOMB = C*ALAT/10.

Tdamp=str(20*float(delta_t))
Pdamp=str(20*float(delta_t))
ensemble = ENSEMBLE.lower()

rcutoff=8.00


if __name__=="__main__":
    # NOTE: Running this initial section in parallel leads to missing helium atoms, 
    # likely due to a bug in LAMMPS. To avoid this issue, run this part in serial 
    # before starting any equilibration. This will produce the initial 
    # non-equilibrated configuration. Afterwards, you can re-run the 'restart.ini' 
    # file from any script in parallel without problems.

    import random
    from lammps import lammps
    lmp = lammps(name="mpi", cmdargs=["-pk","omp","4","-sf","omp"])
    # import sys
    # sys.path.append('../../../../../SCRIPTS_LiPb')
    # sys.path.append(PARENTDIR+"SCRIPTS_LiPbHe")
    # with open("parent-directory", "r") as pdir: PARENTDIR = pdir.read().replace('\n','').replace('"','')
    import COMANDES, ADDITIONAL

    # Generate two random seeds (6 and 7 digits, respectively)
    # and write them in a file ('SEEDS'):
    seed1 = random.randint(100000, 999999) #987654
    seed2 = random.randint(int(1e6), 9999999) #4928459
    with open("SEEDS", "w") as Seeds:
        Seeds.write("SEED 1:\t"+str(seed1)+"\nSEED 2:\t"+str(seed2)+"\n")
    
    # Number of atoms (lithium, helium and lead), pressure and temperature
    # in the current directory:
    [N_Li, N_He, N_Pb], pressio, temperatura = FUNCIONS.READ_NPT()

    # Create lattice and atoms:
    lmp.commands_list(["units metal",
                    "dimension 3",
                    "boundary p p p",
                    "atom_style atomic",
                    "lattice "+LATSTYLE.lower()+" "+str(ALAT),
                    "region whole block "+str(x0)+" "+str(xf)+" "+str(y0)+" "+str(yf)+" "+str(z0)+" "+str(zf),
                    "create_box 3 whole",
                    "mass 1 "+str(m_Li),
                    "mass 2 "+str(m_He),
                    "mass 3 "+str(m_Pb)])
    
    lmp.command("create_atoms 1 box") # NOTE: The initial crystalline structure is made up exclusively by lithium atoms.
    # For an efficient equilibration, it is faster to form a helium bubble when there are no lead atoms.
    # After having built a single bubbble (all helium atoms have gathered), 'N_Pb' lithium atoms should be reconverted to lead atoms,
    # according to the atomic composition of the current simulation indicated in the configuration file.
    # 'x0', 'xf', 'y0', 'yf', 'z0' and 'zf' should be compatible with 
    # (a) the total number of liquid metal atoms;
    # (b) the indicated type of crystalline structure;
    # (c) the density.

    # A spherical region 'esfera' is defined, inside which all 'N_He' helium atoms will be generated.
    # This region may also contain several lithium atoms, so that both type of atoms will be mixed.
    # One should run a few equilibration runs to segregate the atoms.
    # The segregation will be also speeded up if the first equilibration steps are run using Lennard-Jones potentials instead of 
    # realistic ones 
    # (a) It saves computational time (because the Lennard-Jones potentials are computationally efficient in LAMMPS);
    # (b) It also reduces the time in simulation (because Lennard-Jones are more harshly repulsive than the TTS potentials, increasing segregation rates).
    lmp.commands_list([COMANDES.LMP_PRINT("INSERTION OF HELIUM ATOMS"),
                    "region esfera sphere 0 0 0 "+str(RBOMB),
                    # "region esfera sphere 0 0 0 7.8",
                    # "create_atoms 2 random 512 987654 esfera",
                    "create_atoms 2 random "+str(N_He)+" "+str(seed1)+" esfera overlap 0.5",
                    COMANDES.LMP_PRINT("")])
    
    # As explained above, there was an issue with missing helium atoms in simulations when the initialization was performed directly in the parallel runs.
    # We may check the number of 'existing' helium atoms.
    lmp.command("group existing type 2")
    lmp.command("variable gcount equal count(existing)")
    count = int(lmp.extract_variable("gcount", "existing", 0))
    if count < N_He:
        print(str(count)+" < "+str(N_He))
        print("... attempting to create "+str(N_He-count)+" atoms")
        ADDITIONAL.MINIMIZATION(lmp)
        lmp.command("create_atoms 2 random "+str(N_He-count)+" "+str(seed1)+" esfera overlap 0.5")
        count = int(lmp.extract_variable("gcount", "existing", 0))
    else:
        print(":)")
    print("--->",count)

    # Then, we generate the initial velocities using a Gaussian distribution compatible with the temperature of the simulation -- Maxwell-Boltzmann.
    lmp.command("velocity all create " + str(temperatura) + " " + str(seed2) + " dist gaussian")
    # The 'timestep' must be indicated in picoseconds (as we are working with 'metal' units)...
    lmp.command("timestep "+str(delta_t))
    # Additional setup:
    lmp.commands_list(["thermo 0",
                    "thermo_style custom step temp press vol pe ke etail"])
    lmp.commands_list(["neighbor 2.0 bin",
                    "neigh_modify every 5 delay 0 check  yes",
                    "neigh_modify one 10000",
                    "comm_modify cutoff 30.0"])

    # The initalization is done... Just save the 'restart.ini' file to be re-run somewhere else. 
    print("The system has been initialized...")
    lmp.command("write_restart restart.ini")
    exit()
