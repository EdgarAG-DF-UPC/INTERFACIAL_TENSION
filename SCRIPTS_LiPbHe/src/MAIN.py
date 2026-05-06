# IMPORTING LIBRARIES...
import numpy as np
import os, subprocess, shutil, sys, random
from scipy.optimize import curve_fit
from mpi4py import MPI # Note that without the mpi4py specific lines from test.py running the script with mpirun on P processors would lead to P independent simulations to run parallel, each with a single processor. Therefore, if you use the mpi4py lines and you see multiple LAMMPS single processor outputs, mpi4py is not working correctly.
nprocs = MPI.COMM_WORLD.Get_size()
me = MPI.COMM_WORLD.Get_rank()
master = me==0
from lammps import lammps, MPIAbortException
lmp = lammps(name="mpi", cmdargs=["-pk","omp","4","-sf","omp"])#,"-screen","none" .... -sf omp indica que es faci servir la versió omp (si existeix): "will automatically append “omp” to styles that support "

# IMPORTING MY OWN PYTHON SCRIPTS AND FUNCTIONS, AND SETTING THE PARAMETERS OF THE SIMULATION...
with open("parent-directory", "r") as pdir: PARENTDIR = pdir.read().replace('\n','').replace('"','')
# sys.path.append('../../../../../SCRIPTS_LiPb')
sys.path.append(PARENTDIR+"/SCRIPTS_LiPbHe/")

import lmpcmds.FUNCIONS as FUNCIONS, lmpcmds.COMANDES as COMANDES, lmpcmds.ADDITIONAL as ADDITIONAL
# from POTENCIALS import UHeHe as V_22, ULiHe as V_12, UPbHe as V_23
from src.CONFIGURACIO import CUA, RESTART, METHOD as WHICH, BLOCKS, AFEGEIX_YUKAWA, runs_NVT, runs_NPT, Nevery, x0, y0, z0, xf, yf, zf, RBOMB, temperatura, Tdamp, pressio, Pdamp, ensemble, N_Li, N_Pb, N_He, rcutoff, DIFUSIONS_Li_He_Pb, PRESSIONS_PARCIALS, MODEL, SIGNE, RSTFILE
DELTA = float(runs_NPT) / float(Nevery)
RESTART = FUNCIONS.evalua(RESTART)
try:
    TOT = int(BLOCKS)
except NameError as err:
    if CUA=="brief":
        TOT = nprocs*10
    elif CUA=="all":
        TOT = nprocs*100
    else:
        TOT = nprocs

if WHICH in ["THOMPSON", "THERMO_INTEG"]:
    TYPE = "BUBBLES"
elif WHICH in ["PLANAR_VIRIAL", "PLANAR_IK"]:
    TYPE = "PLANAR"
else:
    print(WHICH)
    raise NameError(f"Variable 'WHICH'={WHICH} should be equal either 'THOMPSON', 'THERMO_INTEG', 'PLANAR_VIRIAL' or 'PLANAR_IK'")
from src.INITIALIZATION import INITIALIZE
INI = getattr(INITIALIZE, TYPE)
lmp, RESTART, RSTFILE = INI(N_Li, N_Pb, N_He,
          temperatura, pressio, Tdamp, Pdamp,
          rcutoff,
          runs_NVT, runs_NPT, Nevery,
          RBOMB, [x0, xf, y0, yf, z0, zf],
          ensemble, RESTART, RSTFILE, AFEGEIX_YUKAWA, PRESSIONS_PARCIALS, DIFUSIONS_Li_He_Pb)
if master: print("====>"+WHICH)
from src.METHODS import METHOD
COSA = getattr(METHOD, WHICH)
COSA(lmp)
if master: print("THE END")

MPI.Finalize()
