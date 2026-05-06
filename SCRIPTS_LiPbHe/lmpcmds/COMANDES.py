from mpi4py import MPI 
nprocs = MPI.COMM_WORLD.Get_size()
me = MPI.COMM_WORLD.Get_rank()
master = me==0

def LMP_PRINT(text):
    return "print '==================================================="+str(text)+"====================================================='"

def SET_REGIONS(lmp, x0, xf, y0, yf, z0, zf, RBOMB):
    lmp.commands_list(["region whole block "+str(x0)+" "+str(xf)+" "+str(y0)+" "+str(yf)+" "+str(z0)+" "+str(zf),
                       "region esfera sphere 0 0 0 "+str(RBOMB)])
    

def SET_NEIGH(lmp):
    lmp.commands_list(["neighbor 2.0 bin",
                       "neigh_modify every 5 delay 0 check  yes",
                       "neigh_modify one 10000",
                       "comm_modify cutoff 30.0"])


def SET_THERMO(lmp, Nout, Vars):
    # lmp.commands_list(["thermo 0",
    #                    "thermo_style custom step temp press vol pe ke etail"])
    VAROUT = ""
    for Var in Vars:
        VAROUT += " "+Var
    lmp.commands_list(["thermo "+str(Nout),
                       "thermo_style custom"+VAROUT])


def SET_VARIABLES(lmp, **kwargs):
    lmp.commands_list(["variable N equal step",
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
                       "variable H2 equal v_H*v_H"])

def SET_POTENTIAL(lmp, **kwargs):
    parentdir = kwargs.get("parent", "/users/edgar/INTERFACIAL_TENSION")
    CMDLST_0 = "pair_style hybrid eam/alloy table linear 750 table linear 750 table linear 750"
    if kwargs.get("Born", False):
        CMDLST_0 = CMDLST_0.replace("hybrid", "hybrid/overlay")
        CMDLST_0 += " born 15.00"
    if kwargs.get("CD-EAM", False):
        CMDLST_0 = CMDLST_0.replace("eam/alloy", "eam/cd")
    CMDLST = [CMDLST_0,
              "pair_coeff * * eam/alloy LiPb.eam.alloy Li NULL Pb",
              "pair_coeff 1 2 table 1 "+parentdir+"/TAULES/LiHe.table SHENGTANG_LiHe",
              "pair_coeff 2 2 table 2 "+parentdir+"/TAULES/HeHe.table SHENGTANG_HeHe",
              "pair_coeff 2 3 table 3 "+parentdir+"/TAULES/PbHe.table SLADEK_PbHe"]
#              "pair_coeff 1 2 table 1 ../../../../../TAULES/LiHe.table SHENGTANG_LiHe",
#              "pair_coeff 2 2 table 2 ../../../../../TAULES/HeHe.table SHENGTANG_HeHe",
#              "pair_coeff 2 3 table 3 ../../../../../TAULES/PbHe.table SLADEK_PbHe"]
    if not kwargs.get("helium", True):
        from FUNCIONS import CreateBornTables
        if master: CreateBornTables(0.00, 0.00)
        CMDLST[2] = "pair_coeff 1 2 table 1 LiHe_Born.table SHENGTANG_LiHe"
        CMDLST[3] = "pair_coeff 2 2 table 2 HeHe_Born.table SHENGTANG_HeHe"
        CMDLST[4] = "pair_coeff 2 3 table 3 PbHe_Born.table SLADEK_PbHe"
    
    if kwargs.get("Born", False):
       A = kwargs.get("A", 12.483)
       B = kwargs.get("B", 1.)
       λ = kwargs.get("Lambda", -2.00)
       CMDLST.append("pair_coeff 1 2 born "+str(A)+" "+str(B)+" "+str(λ)+" 0 0")
       CMDLST.append("pair_coeff 2 3 born "+str(A)+" "+str(B)+" "+str(λ)+" 0 0")

    if kwargs.get("BornTable", False):
        CMDLST = [CMDLST_0,
                  "pair_coeff * * eam/alloy LiPb.eam.alloy Li NULL Pb",
                  "pair_coeff 1 2 table 1 LiHe_Born.table SHENGTANG_LiHe",
                  "pair_coeff 2 2 table 2 HeHe_Born.table SHENGTANG_HeHe",
                  "pair_coeff 2 3 table 3 PbHe_Born.table SLADEK_PbHe"]
    
    lmp.commands_list(CMDLST) # SET EAM - TTS PARAMETERS
       

def FIX_ENSEMBLE(lmp, temperatura, Tdamp, pressio, Pdamp, ensemble):
    lmp.command(LMP_PRINT("PRODUCTION RUNS"))
    lmp.command("fix 1 all " + str(ensemble) + " temp " + str(temperatura) + " " + str(temperatura) + " " + str(Tdamp)+ " iso " + str(pressio) + " " + str(pressio) + " " + str(Pdamp))

def FIX_THERMOSTAT(ensemble, temperatura, Tdamp, **kwargs):
    cmd = "fix 1 all "+ensemble.lower()+" temp " + str(temperatura) + " " + str(temperatura) + " " + str(Tdamp)
    if ensemble=="NPT":
        pressio = kwargs.get("press", None)
        Pdamp = kwargs.get("Pdamp", None)
        cmd += " iso " + str(pressio) + " " + str(pressio) + " " + str(Pdamp)
    elif ensemble.lower()!="nvt":
        raise RuntimeError("The ensemble must be either NPT or NVT...")
    return cmd



def GROUPS(lmp, grp):
    COMS = []
    for name, nums in grp:
        linia = "group "+name+" type"        
        for n in nums:
            linia += " "+str(n)
        COMS.append(linia)
    lmp.commands_list(COMS)
    # lmp.commands_list(["group   Li		type 1",
    #                    "group   He		type 2",
    #                    "group   Pb		type 3",
    #                    "group	liquid	type 1 3"])
    
def SET_OUTPUT(lmp, proc, Nevery, Nrepeat, runs_NPT, xyz_file=True, velocity=False, boxlims=True, rdf=True, centermass=True, thermodata=True, averages=True, stress=False): #, **kwargs
    # if kwargs.get("xyz_file", True):
    if xyz_file: lmp.commands_list(["dump 1 all xyz "+Nevery+" "+str(proc)+"_npt.xyz",
                                    "dump_modify 1 pbc yes"])
    # if kwargs.get("velocity", False):
    if velocity: lmp.commands_list(["dump 2 all custom "+Nevery+" "+str(proc)+"_velocities.out vx vy vz"])
    # if kwargs.get("boxlims", True):
    if boxlims: lmp.command("fix 2 all print "+Nevery+"  '${N} ${XXl} ${XXh} ${YYl} ${YYh} ${ZZl} ${ZZh}' file "+str(proc)+"_simulation_box.out screen no")
    # if kwargs.get("rdf", True):
    if rdf: lmp.commands_list(["compute myRDF all rdf 500 1 1 1 2 1 3 2 2 2 3 3 3 cutoff 25.0",
                           "fix 3 all ave/time "+Nevery+" "+Nrepeat+" "+runs_NPT+" c_myRDF[*] file "+str(proc)+"_aver.rdf mode vector"])
    # if kwargs.get("centermass", True):
    if centermass: lmp.commands_list(["compute CoM He com",
                                 "variable xCM equal c_CoM[1]",
                                 "variable yCM equal c_CoM[2]",
                                 "variable zCM equal c_CoM[3]",
                                 "variable time equal step",
                                 "fix	4 all print "+Nevery+" '${time} ${xCM} ${yCM} ${zCM}' file "+str(proc)+"_center_mass.out screen no"])
    # if kwargs.get("thermodata", True):
    if thermodata: lmp.command("fix	5 all print "+Nevery+" '${N} ${T} ${Press} ${V} ${Etotal} ${pote} ${kine} ${Etail}' file "+str(proc)+"_thermo_data.out screen no")
    # if kwargs.get("averages", True):
    if averages: lmp.command("fix 6 all ave/time "+Nevery+" "+Nrepeat+" "+runs_NPT+"	c_thermo_temp v_T2 c_thermo_press v_P2 v_V v_V2 v_Etotal v_E2 v_PV v_P2V2 v_H v_H2 file "+str(proc)+"_averages.out")
    # if kwargs.get("stress", False):
        # lmp.commands_list(["compute peratom all stress/atom NULL",
        #                    "compute "])
        # lmp.commands_list(["variable nbins index 500",
        #                    "variable fraction equal 1.0/v_nbins",
        #                    "variable zCoM equal xcm(He,z)",
        #                    "compute cchunk all chunk/atom bin/1d z ${zCoM} ${fraction} units reduced",
        #                    "compute stress all stress/atom NULL",
        #                    "variable press atom -(c_stress[1]+c_stress[2]+c_stress[3])/(3.0*vol*${fraction})",
        #                    "variable pxx atom -(c_stress[1])/(vol*${fraction})",
        #                    "variable pyy atom -(c_stress[2])/(vol*${fraction})",
        #                    "variable pzz atom -(c_stress[3])/(vol*${fraction})",
        #                    "variable dpress atom (0.5*(c_stress[1]+c_stress[2])-(c_stress[3]))/(vol*${fraction})",
        #                    "variable dstress atom (0.5*(c_stress[1]+c_stress[2])-(c_stress[3]))",
        #                    "compute binpress all reduce/chunk cchunk sum v_press",
        #                    "compute bindpress all reduce/chunk cchunk sum v_dpress",
        #                    "compute binpxx all reduce/chunk cchunk sum v_pxx",
        #                    "compute binpyy all reduce/chunk cchunk sum v_pyy",
        #                    "compute binpzz all reduce/chunk cchunk sum v_pzz",
        #                    "compute binstr all reduce/chunk cchunk sum v_dstress",
        #                    "fix avg all ave/time "+Nevery+" "+Nrepeat+" "+runs_NPT+" c_binpxx c_binpyy c_binpzz c_binpress c_bindpress mode vector file ave_stress_"+str(proc)+".out",
        #                    "variable pxxtot atom -(c_stress[1])/(vol)",
        #                    "variable pyytot atom -(c_stress[2])/(vol)",
        #                    "variable pzztot atom -(c_stress[3])/(vol)",
        #                    "compute PTx all reduce sum v_pxxtot",
        #                    "compute PTy all reduce sum v_pyytot",
        #                    "compute PTz all reduce sum v_pzztot",
        #                    "fix pt all ave/time "+Nevery+" "+Nrepeat+" "+runs_NPT+" c_PTx c_PTy c_PTz file pressure_tensor_"+str(proc)+".out"])
    if stress:
        lmp.commands_list(["variable nbins index 500",
                           "variable fraction equal 1.0/v_nbins",
                           "variable zCoM equal xcm(He,z)",
                           "compute cchunk all chunk/atom bin/1d z ${zCoM} ${fraction} units reduced",
                           "compute stress all stress/atom NULL",
                           "variable press atom -(c_stress[1]+c_stress[2]+c_stress[3])/(3.0*v_V*${fraction})",
                           "variable pxx atom -(c_stress[1])/(v_V*${fraction})",
                           "variable pyy atom -(c_stress[2])/(v_V*${fraction})",
                           "variable pzz atom -(c_stress[3])/(v_V*${fraction})",
                           "variable dpress atom (0.5*(c_stress[1]+c_stress[2])-(c_stress[3]))/(v_V*${fraction})",
                           "variable dstress atom (0.5*(c_stress[1]+c_stress[2])-(c_stress[3]))",
                           "compute binpress all reduce/chunk cchunk sum v_press",
                           "compute bindpress all reduce/chunk cchunk sum v_dpress",
                           "compute binpxx all reduce/chunk cchunk sum v_pxx",
                           "compute binpyy all reduce/chunk cchunk sum v_pyy",
                           "compute binpzz all reduce/chunk cchunk sum v_pzz",
                           "compute binstr all reduce/chunk cchunk sum v_dstress",
                           "fix avg all ave/time "+Nevery+" "+Nrepeat+" "+runs_NPT+" c_binpxx c_binpyy c_binpzz c_binpress c_bindpress c_binstr mode vector file ave_stress_"+str(proc)+".out",
                           "variable pxxtot atom -(c_stress[1])/(v_V)",
                           "variable pyytot atom -(c_stress[2])/(v_V)",
                           "variable pzztot atom -(c_stress[3])/(v_V)",
                        #    "dump sttot all custom 1 "+str(proc)+"_stress_per_atom.out id type c_stress[1] c_stress[2] c_stress[3]",#"dump sttot all custom "+Nevery+" "+str(proc)+"_stress_per_atom.out id type c_stress[1] c_stress[2] c_stress[3]",
                           "compute PTx all reduce sum v_pxxtot",
                           "compute PTy all reduce sum v_pyytot",
                           "compute PTz all reduce sum v_pzztot",
                           "fix pt all ave/time "+Nevery+" "+Nrepeat+" "+runs_NPT+" c_PTx c_PTy c_PTz file pressure_tensor_"+str(proc)+".out"])


def UNFIX_UNCOMPUTE(lmp,  xyz_file=True, velocity=False, boxlims=True, rdf=True, centermass=True, thermodata=True, averages=True, stress=False): #**kwargs
    # lmp.commands_list(["undump    1",
    #                    "unfix     2",
    #                    "unfix     3",
    #                    "unfix     4",
    #                    "unfix     5",
    #                    "unfix     6",
    #                    "uncompute CoM",
    #                    "uncompute myRDF"])
    # if kwargs.get("xyz_file", True):
    if xyz_file: lmp.command("undump 1")
    # if kwargs.get("velocity", False):
    if velocity: lmp.command("undump 2")
    # if kwargs.get("boxlims", True): 
    if boxlims: lmp.command("unfix 2")
    # if kwargs.get("rdf", True): 
    if rdf: lmp.commands_list(["unfix 3", "uncompute myRDF"])
    # if kwargs.get("centermass", True): 
    if centermass: lmp.commands_list(["unfix 4", "uncompute CoM"])
    # if kwargs.get("thermodata", True): 
    if thermodata: lmp.command("unfix 5")
    # if kwargs.get("averages", True): 
    if averages: lmp.command("unfix 6")
    # if kwargs.get("stress", False): 
    if stress: lmp.commands_list(["uncompute cchunk",
                                  "uncompute stress",
                                  "uncompute binpress",
                                  "uncompute bindpress",
                                  "uncompute binpxx",
                                  "uncompute binpyy",
                                  "uncompute binpzz",
                                  "uncompute binstr",
                                  "uncompute PTx",
                                  "uncompute PTy",
                                  "uncompute PTz",
                                  "unfix pt",
                                  "unfix avg"])#,
                                                    #   "undump sttot"
