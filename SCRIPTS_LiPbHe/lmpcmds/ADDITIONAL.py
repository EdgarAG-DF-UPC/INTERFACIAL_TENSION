def LENNARD_JONES(lmp, rcutoff, **kwargs):
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
    if not kwargs.get("helium", True):
        epsilon_HeHe = 0.
        epsilon_LiHe = 0.
        epsilon_HePb = 0.
        sigma_HeHe = 0.
        sigma_LiHe = 0.
        sigma_HePb = 0.

    lmp.commands_list(["pair_style lj/cut "+ str(rcutoff),
                       "pair_coeff 1 1 "+ str(epsilon_LiLi) + " " + str(sigma_LiLi),
                       "pair_coeff 2 2 "+ str(epsilon_HeHe) + " " + str(sigma_HeHe),
                       "pair_coeff 3 3 "+ str(epsilon_PbPb) + " " + str(sigma_PbPb),
                       "pair_coeff 1 2 "+ str(epsilon_LiHe) + " " + str(sigma_LiHe),
                       "pair_coeff 1 3 "+ str(epsilon_LiPb) + " " + str(sigma_LiPb),
                       "pair_coeff 2 3 "+ str(epsilon_HePb) + " " + str(sigma_HePb)]) # MÒDUL D'EQUILIBRACIÓ LENNARD-JONES PRÈVIA
    
def MINIMIZATION(lmp):
    lmp.commands_list(["min_style cg",
                       "minimize 1.0e-8  1.0e-8  10000     100000",
                       "reset_timestep 0"])
    
def THERMAL_SHAKE(lmp, temperatura, Tshake, Tdamp, runs_shake):
    lmp.commands_list(["fix 1 all nvt temp " + str(temperatura) + " " + str(Tshake) + " " + str(Tdamp),
                       "run "+ str(runs_shake),
                       "unfix 1",
                       "fix	1 all nvt temp " + str(Tshake) + " " + str(Tshake) + " " + str(Tdamp),
                       "run " + str(runs_shake),
                       "unfix 1",
                       "fix	1 all nvt temp " + str(Tshake) + " " + str(temperatura) + " " + str(Tdamp),
                       "run " + str(runs_shake),
                       "unfix 1"])

def EQUILIBRACIO_NVT(lmp, temperatura, Tdamp, Nevery, runs_NVT):
    lmp.command("print '========================================= EQUILIBRACIÓ =========================================='")
    lmp.commands_list(["fix 1 all nvt temp " + str(temperatura) + " " + str(temperatura) + " " + str(Tdamp),
                       "fix 2 all print "+Nevery+"  '${N} ${XXl} ${XXh} ${YYl} ${YYh} ${ZZl} ${ZZh}' file equil_simulation_box.out"
                       "run " + str(runs_NVT),
                       "unfix 1",
                       "unfix 2"])

def ADD_YUKAWA(lmp, N_Pb, N_Li, rcutoff, temperatura, Tdamp, runs_NVT):
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
        print("ERROR !!!!!!!")
        quit()

    a_LiLi = "{:.8f}".format(qLi*qLi*14.399)
    a_PbPb = "{:.8f}".format(qPb*qPb*14.399)
    a_LiPb = "{:.8f}".format(qLi*qPb*14.399)
    lmp.commands_list(["print '--------------------------------- ADDICIÓ POTENCIAL YUKAWA----------------------------------'",
                        "pair_style	hybrid/overlay  eam/alloy table linear 750 table linear 750 table linear 750 yukawa 1.10 " + str(rcutoff) + " yukawa 1.10 " + str(rcutoff) + " yukawa 1.10 " + str(rcutoff),
                        "pair_coeff * * eam/alloy LiPb.eam.alloy Li NULL Pb",
                        "pair_coeff 1 2 table 1 ../../../../TAULES/LiHe.table SHENGTANG_LiHe",
                        "pair_coeff 2 2 table 2 ../../../../TAULES/HeHe.table SHENGTANG_HeHe",
                        "pair_coeff 2 3 table 3 ../../../../TAULES/PbHe.table SLADEK_PbHe",
                        "pair_coeff 1 1 yukawa 1 " + str(a_LiLi),
                        "pair_coeff 1 3 yukawa 2 " + str(a_LiPb),
                        "pair_coeff 3 3 yukawa 3 " + str(a_PbPb),
                        "fix 1 all nvt temp " + str(temperatura) + " " + str(temperatura) + " " + str(Tdamp),
                        "run " + str(runs_NVT),
                        "unfix 1"]) # MÒDUL DE POTENCIAL YUKAWA - CORRECCIÓ BELASHCHENKO 2019
    lmp.commands_list(["fix 1 all nvt temp " + str(temperatura) + " " + str(temperatura) + " " + str(Tdamp),
                        "run " + str(runs_NVT),
                        "unfix 1"])


def PRESSIONS_PARCIALS(lmp, Nevery, proc):
    lmp.commands_list(["compute peratom all stress/atom NULL",
                       "compute ST1 Li reduce sum c_peratom[1] c_peratom[2] c_peratom[3]",
                       "compute ST2 He reduce sum c_peratom[1] c_peratom[2] c_peratom[3]",
                       "compute ST3 Pb reduce sum c_peratom[1] c_peratom[2] c_peratom[3]",
                       "variable ppLi equal -(c_ST1[1]+c_ST1[2]+c_ST1[3])/(3*vol)",
                       "variable ppHe equal -(c_ST2[1]+c_ST2[2]+c_ST2[3])/(3*vol)",
                       "variable ppPb equal -(c_ST3[1]+c_ST3[2]+c_ST3[3])/(3*vol)",
                       "variable STxx equal c_ST1[1]",
                       "variable STyy equal c_ST1[2]",
                       "variable STzz equal c_ST1[3]",
                       "fix extra4 all print "+Nevery+" '${N} ${ppLi} ${ppHe} ${ppPb} ${Press}' file "+str(proc)+"_partial_pressures.out"]) # MÒDUL DE PRESSIONS PARCIALS


def DIFUSIONS(lmp, Nevery, Nrepeat, runs_NPT, proc):
    lmp.commands_list(["compute msd1 Li msd",
                       "variable msdLi equal c_msd1[4]",
                       "compute msd2 He msd",
                       "variable msdHe equal c_msd2[4]",
                       "compute msd3 Pb msd",
                       "variable msdPb equal c_msd3[4]",
                       "compute vacf1 Li vacf",
                       "variable vacfLi equal c_vacf1[4]",
                       "compute vacf2 He vacf",
                       "variable vacfHe equal c_vacf2[4]",
                       "compute vacf3 Pb vacf",
                       "variable vacfPb equal c_vacf3[4]",
                       "fix extra9 all ave/time 1 "+Nevery+" "+Nevery+" v_msdLi v_msdHe v_msdPb v_vacfLi v_vacfHe v_vacfPb file "+str(proc)+"_dif.out"]) # MÒDUL DE DIFUSIONS

def UNDIF(lmp):
    lmp.commands_list(["unfix extra9",
                       "uncompute msd1",
                       "uncompute msd2",
                       "uncompute msd3",
                       "uncompute vacf1",
                       "uncompute vacf2",
                       "uncompute vacf3"])

def UNPRES(lmp):
    lmp.commands_list(["unfix extra4",
                 "uncompute ST1",
                 "uncompute ST2",
                 "uncompute ST3",
                 "uncompute peratom"])



def BORN(lmp, Nevery, Nrepeat, runs_NPT, Lambda):
    lmp.commands_list(["compute born all pair born",
                       "variable born2 equal c_born*c_born",
                       "compute uTTS all pair table 1",
                       "compute	uSla all pair table 3",
                       "variable uTTS2 equal c_uTTS*c_uTTS",
                       "variable uSla2 equal c_uSla*c_uSla",
                       "fix extra10 all ave/time "+Nevery+" "+Nrepeat+" "+runs_NPT+" c_born v_born2 c_uTTS v_uTTS2 c_uSla v_uSla2 file "+str(Lambda)+"_born_tts_s.out"])


def UNBORN(lmp):
    lmp.commands_list(["unfix extra10",
                       "uncompute born",
                       "uncompute uTTS",
                       "uncompute uSla"])