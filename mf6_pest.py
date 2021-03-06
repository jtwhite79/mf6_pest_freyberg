import os
import shutil
import string
import numpy as np
import pandas as pd
import platform
import spnspecs

from matplotlib.patches import Polygon
from matplotlib.backends.backend_pdf import PdfPages as pdf
import matplotlib.pyplot as plt
spnspecs.set_graph_specifications()
spnspecs.set_map_specifications()
import flopy
import pyemu

plt_dir = "plots"
if not os.path.exists(plt_dir):
    os.mkdir(plt_dir)

test_dict = {"master_ies":["freyberg6_run_ies.3.par.csv","freyberg6_run_ies.3.obs.csv","freyberg6_run_ies.phi.group.csv"],
                 "master_glm":["freyberg6_run_glm.ipar","freyberg6_run_glm.iobj"],
                 "master_sen_morris":["freyberg6_run_sen.sen.par.csv","freyberg6_run_sen.mio"],
                 "master_opt_neutral":["freyberg6_run_opt.1.est.rei","freyberg6_run_opt.1.sim.rei"],
                 "master_opt_averse":["freyberg6_run_opt.1.sim+chance.rei","freyberg6_run_opt.1.est+chance.rei"]}
unit_dict = {"head":"sw-gw flux $\\frac{ft^3}{d}$",
                "tail": "sw-gw flux $\\frac{ft^3}{d}$",
                "trgw" : "gw level $ft$",
                "gage" : "sw flux $\\frac{ft^3}{d}$"}
label_dict = {"head": "headwater",
             "tail": "tailwater",
             "trgw_2_2_9": "gw_1",
              "trgw_2_33_7": "gw_2",
              "trgw_0_9_1" : "gw_3",
             "gage": "sw_1"}

#forecasts = ["headwater_20171231","tailwater_20161130","trgw_0_9_1_20171130"]
forecasts = ["headwater_20171130","tailwater_20161130","trgw_0_9_1_20161130"]

abet = string.ascii_uppercase
def prep_mf6_model():
    org_ws = "temp_history"
    if True:#not os.path.exists(os.path.join(org_ws,"freyberg6.nam")):

        #first mod the nam file and the dumbass last reach bottom - sigh
        m = flopy.modflow.Modflow.load("freyberg.nam",model_ws=org_ws,check=False)
        print(m.sfr.reach_data.dtype)
        last_reach_bot = m.sfr.reach_data["strtop"][-1] - m.sfr.reach_data["strthick"][-1]
        cell_bot = m.dis.botm[0].array[m.sfr.reach_data["i"][-1],m.sfr.reach_data["j"][-1]]
        print(cell_bot,last_reach_bot)
        if last_reach_bot <= cell_bot:
            m.sfr.reach_data["strtop"][-1] += 0.001
        m.wel.stress_period_data[7]["flux"] = m.wel.stress_period_data[8]["flux"] * 1.01
        m.external_path = "."

        m.write_input()
        lines = open(os.path.join(org_ws,"freyberg.nam"),'r').readlines()
        with open(os.path.join(org_ws,"freyberg_mod.nam"),'w') as f:
            for line in lines:
                if ".sfr.out" in line and not "REPLACE" in line:
                    line = line.strip() + " REPLACE\n"
                f.write(line)
        pyemu.os_utils.run("mf5to6 freyberg_mod.nam freyberg6",cwd=org_ws)


    new_ws = "test"
    if os.path.exists(new_ws):
        shutil.rmtree(new_ws)
    os.mkdir(new_ws)
    sim = flopy.mf6.MFSimulation.load(sim_ws=org_ws)
    sim.simulation_data.mfpath.set_sim_path("test")
    #sim.set_all_data_external()
    sim.name_file.continue_ = True
    m = sim.get_model("freyberg6")

    redis_fac = m.dis.nrow.data / 40 #in case this is a finely discret version

    obs_df = pd.read_csv("obs_loc.csv")
    obs_df.loc[:,"row"] = (obs_df.row * redis_fac).apply(np.int)
    obs_df.loc[:, "col"] = (obs_df.col * redis_fac).apply(np.int)


    obs_df.loc[:,"layer"] = 3
    obs_df.loc[:,"name"] = obs_df.apply(lambda x: "trgw_{0}_{1}_{2}".\
        format(x.layer-1,x.row-1,x.col-1),axis=1)
    obs_df.loc[:,"obstype"] = "HEAD"
    obs_df.loc[:,"col"] = obs_df.col.apply(np.int)
    obs_df2 = obs_df.copy()
    obs_df2.loc[:,"layer"] = 1
    obs_df2.loc[:,"name"] = obs_df2.apply(lambda x: "trgw_{0}_{1}_{2}".\
        format(x.layer-1,x.row-1,x.col-1),axis=1)

    obs_df = pd.concat([obs_df,obs_df2])

    #head_obs = {"head_obs.csv":[("trgw_{0}_{1}".format(r,c),"NPF",(2,r-1,c-1)) for r,c in zip(obs_df.row,obs_df.col)]}
    with open(os.path.join(new_ws,"head.obs"),'w') as f:
        f.write("BEGIN CONTINUOUS FILEOUT heads.csv\n")
        obs_df.loc[:,["name","obstype","layer","row","col"]].to_csv(f,sep=' ',line_terminator='\n',
            index=False,header=False,mode="a")
        f.write("END CONTINUOUS\n")

    props_3d = [("npf","k"),("npf","k33"),("sto","ss"),("sto","sy")]
    props_trans = [("rch","recharge")]
    #print(dir(m.dis.nlay))
    #print(m.dis.nlay.data)
    print(type(m.rch.recharge))

    for pack,attr in props_3d:
        #print(m.get_package(pack).__getattribute__(attr))
        for k in range(m.dis.nlay.data):
            filename = "{0}_{1}_{2}.dat".format(pack,attr,k)
            m.get_package(pack).__getattribute__(attr).store_as_external_file(filename,layer=k)

    sim.write_simulation()
    lines = open(os.path.join(new_ws,"freyberg6.nam"),'r').readlines()
    new_lines = []
    for line in lines:
        if line.strip() == "END Packages":
            new_lines.append("   obs6 head.obs\n")
        new_lines.append(line)
    with open(os.path.join(new_ws,"freyberg6.nam"),'w') as f:
        [f.write(line) for line in new_lines]
        f.flush()

    # mod sfr
    lines = open(os.path.join(new_ws,"freyberg6.sfr"),'r').readlines()
    with open(os.path.join(new_ws,"freyberg6.sfr"),'w') as f:
        iline = 0
        while True:
            if iline >= len(lines):
                break
            line = lines[iline]
            f.write(line)
            if "begin options" in line.lower():
                f.write("  BOUNDNAMES\n")
                f.write("OBS6 FILEIN sfr.obs\n")
            if "begin packagedata" in line.lower():
                iline += 1
                #line = lines[iline]
                irch = 0
                while "end" not in line.lower():
                    print(iline)
                    line = lines[iline]
                    if "end" in line.lower():
                        break
                    bn = "headwater"
                    if irch > int(m.dis.nrow.data / 2):
                        bn = "tailwater"
                    line = line.strip() + " " + bn + "\n"
                    f.write(line)
                    irch += 1
                    iline += 1
                f.write(line)
            iline += 1

    with open(os.path.join(new_ws,"sfr.obs"),'w') as f:
        f.write("BEGIN CONTINUOUS FILEOUT sfr.csv\nheadwater sfr headwater\n")
        f.write("tailwater sfr tailwater\ngage_1 inflow 40\nEND CONTINUOUS")

    #shutil.copy2(os.path.join(org_ws,"mf6.exe"),os.path.join(new_ws,"mf6.exe"))
    pyemu.os_utils.run("mf6",cwd=new_ws)
    #ext_dict = {""}


def setup_pest_interface():
    org_ws = "test"
    assert os.path.exists(org_ws)

    ws = "template"
    if os.path.exists(ws):
        shutil.rmtree(ws)
    shutil.copytree(org_ws,ws)
    obs_df = _write_instuctions(ws)
    par_df = _write_templates(ws)
    io_files = pyemu.helpers.parse_dir_for_io_files(ws,prepend_path=True)
    pst = pyemu.Pst.from_io_files(*io_files,pst_path='.')
    for col in par_df.columns:
        pst.parameter_data.loc[par_df.parnme,col] = par_df.loc[:,col]
    for col in obs_df.columns:
        pst.observation_data.loc[obs_df.obsnme,col] = obs_df.loc[:,col]
    if "window" in platform.platform().lower():
        pst.model_command = "mf6"
        shutil.copy2(os.path.join("bin","win","mf6.exe"),os.path.join(ws,"mf6.exe"))
    elif "linux" in platform.platform().lower():
        pst.model_command = "./mf6"
        shutil.copy2(os.path.join("bin", "linux", "mf6"), os.path.join(ws, "mf6"))
    else:
        pst.model_command = "./mf6"
        shutil.copy2(os.path.join("bin", "mac", "mf6"), os.path.join(ws, "mf6"))
    pst.control_data.noptmax = 0
    pst.write(os.path.join(ws,"freyberg6.pst"),version=2)
    pyemu.os_utils.run("pestpp-ies freyberg6.pst",cwd=ws)


def build_and_draw_prior():
    t_d = "template"
    sim = flopy.mf6.MFSimulation.load(sim_ws=t_d)
    m = sim.get_model("freyberg6")
    xgrid = m.modelgrid.xcellcenters
    ygrid = m.modelgrid.ycellcenters
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6.pst"))
    par = pst.parameter_data
    static_par = par.loc[par.parnme.apply(lambda x: x[:3] in ["npf","sto"]),:].copy()
    static_par.loc[:, "i"] = static_par.parnme.apply(lambda x: int(x.split('_')[3]))
    static_par.loc[:, "j"] = static_par.parnme.apply(lambda x: int(x.split('_')[4]))
    static_par.loc[:, "x"] = static_par.apply(lambda x: xgrid[x.i,x.j],axis=1)
    static_par.loc[:, "y"] = static_par.apply(lambda x: ygrid[x.i, x.j], axis=1)
    static_par.loc[:,"pargp"] = static_par.parnme.apply(lambda x: "_".join(x.split('_')[:3]))

    wel_par = par.loc[par.parnme.apply(lambda x: x.startswith("wel")),:].copy()
    wel_par.loc[:,"x"] = wel_par.parnme.apply(lambda x: int(x.split('_')[-1]))
    wel_par.loc[:,"y"] = 0.0
    wel_par.loc[:,"pargp"] = wel_par.parnme.apply(lambda x: '_'.join(x.split('_')[:-1]))

    rch_par = par.loc[par.parnme.str.startswith("rch"),:].copy()
    rch_par.loc[:,"x"] = rch_par.parnme.apply(lambda x: int(x.split('_')[-1]))
    rch_par.loc[:,"y"] = 0.0

    spatial_v = pyemu.geostats.ExpVario(contribution=1.0,a=1000.0)
    temporal_v = pyemu.geostats.ExpVario(contribution=1.0,a=3)
    spatial_gs = pyemu.geostats.GeoStruct(variograms=spatial_v)
    temporal_gs = pyemu.geostats.GeoStruct(variograms=temporal_v)

    static_struct_dict = {spatial_gs:[]}
    sgrps = static_par.pargp.unique()
    sgrps.sort()
    for pargp in sgrps:
        static_struct_dict[spatial_gs].append(static_par.loc[static_par.pargp==pargp,["parnme","x","y","i","j"]])
    temporal_struct_dict = {temporal_gs: [rch_par.loc[:, ["parnme", "x", "y"]]]}
    wgrps = wel_par.pargp.unique()
    wgrps.sort()
    for pargp in wgrps:
        temporal_struct_dict[temporal_gs].append(wel_par.loc[wel_par.pargp == pargp, ["parnme", "x", "y"]])

    struct_dict = static_struct_dict
    for k,v in temporal_struct_dict.items():
        struct_dict[k] = v
    print(struct_dict)
    np.random.seed(pyemu.en.SEED)
    pe = pyemu.helpers.geostatistical_draws(pst,struct_dict=struct_dict,num_reals=300)
    pe.to_binary(os.path.join(t_d,"prior.jcb"))


def _write_instuctions(ws):
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
    m = sim.get_model("freyberg6")
    perlen = sim.tdis.perioddata.array["perlen"]
    totim = np.cumsum(perlen).astype(int)
    sp_end = pd.to_datetime("12-31-2015") + pd.to_timedelta(totim,unit='d') - pd.to_timedelta(1,unit='d')
    totim_dict = {t:s.strftime("%Y%m%d") for t,s in zip(totim,sp_end)}

    sfr_obs_df = pd.read_csv(os.path.join(ws, "sfr.csv"), index_col=0)
    names, vals = [], []
    with open(os.path.join(ws, "sfr.csv.ins"), 'w') as f:
        f.write("pif ~\nl1\n")
        for i in sfr_obs_df.index:
            f.write("l1 ")
            for c in sfr_obs_df.columns:
                name = "{0}_{1}".format(c.lower(), totim_dict[int(i)])
                f.write(" ~,~ !{0}! ".format(name))
                names.append(name)
                vals.append(sfr_obs_df.loc[i, c])
            f.write("\n")
    sfr_df = pd.DataFrame({"obsnme": names, "obsval": vals}, index=names)
    sfr_df.loc[:, "obgnme"] = sfr_df.obsnme.apply(lambda x: x.split('_')[0])

    obs_df = pd.read_csv(os.path.join(ws,"heads.csv"), index_col=0)
    names,vals = [],[]
    with open(os.path.join(ws,"heads.csv.ins"),'w') as f:
        f.write("pif ~\nl1\n")
        for i in obs_df.index:
            f.write("l1 ")
            for c in obs_df.columns:
                name = "{0}_{1}".format(c.lower(), totim_dict[int(i)])
                f.write("~,~ !{0}! ".format(name))
                names.append(name)
                vals.append(obs_df.loc[i,c])
            f.write("\n")
    head_df = pd.DataFrame({"obsnme":names,"obsval":vals},index=names)
    head_df.loc[:,"obgnme"] = head_df.obsnme.apply(lambda x: '_'.join(x.split('_')[:4]))

    tag = "VOLUME BUDGET FOR ENTIRE MODEL AT END OF TIME STEP    1, STRESS PERIOD"

    with open(os.path.join(ws,"freyberg6.lst.ins"),'w') as f:
        f.write("pif ~\n")
        for kper in range(25):
            f.write("~{0}{1:4d}~ \n".format(tag,kper+1))
            kper = sp_end[kper].strftime("%Y%m%d")
            f.write("l8 w w w w w w !in_ss_flx_{0}! \n".format(kper))
            f.write("l1 w w w w w w !in_sy_flx_{0}! \n".format(kper))
            f.write("l1 w w w w w w !in_wel_flx_{0}! \n".format(kper))
            f.write("l1 w w w w w w !in_ghb_flx_{0}! \n".format(kper))
            f.write("l1 w w w w w w !in_rch_flx_{0}! \n".format(kper))
            f.write("l1 w w w w w w !in_sfr_flx_{0}! \n".format(kper))
            f.write("l6 w w w w w w !out_ss_flx_{0}! \n".format(kper))
            f.write("l1 w w w w w w !out_sy_flx_{0}! \n".format(kper))
            f.write("l1 w w w w w w !out_wel_flx_{0}! \n".format(kper))
            f.write("l1 w w w w w w !out_ghb_flx_{0}! \n".format(kper))
            f.write("l1 w w w w w w !out_rch_flx_{0}! \n".format(kper))
            f.write("l1 w w w w w w !out_sfr_flx_{0}! \n".format(kper))
    ins = pyemu.pst_utils.InstructionFile(os.path.join(ws,"freyberg6.lst.ins"))
    #ins.read_output(os.path.join(ws,"freyberg6.lst"))
    flx_df = ins.read_output_file(os.path.join(ws,"freyberg6.lst"))
    flx_df.loc[:,"obgnme"] = flx_df.index.map(lambda x: '_'.join(x.split('_')[:2]))
    flx_df.loc[:,"obsnme"] = flx_df.index.values
    #print(flx_df)
    df = pd.concat([head_df,flx_df,sfr_df],sort=False)
    return df

def _write_templates(ws):
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
    m = sim.get_model("freyberg6")
    nrow = m.dis.nrow.data
    ncol = m.dis.ncol.data
    # write a wel file template
    f_in = open(os.path.join(ws,"freyberg6.wel"),'r')
    with open(os.path.join(ws,"freyberg6.wel.tpl"),'w') as f:
        f.write("ptf ~\n")
        names,vals = [],[]
        while True:
            line = f_in.readline()
            if line == '':
                break
            f.write(line)
            if "BEGIN PERIOD" in line:
                kper = int(line.strip().split()[-1]) - 1
                while True:
                    line = f_in.readline()
                    if line == '':
                        raise Exception()
                    if "END PERIOD" in line:
                        f.write(line)
                        break
                    raw = line.strip().split()
                    kij = [int(r)-1 for r in raw[:-1]]
                    val = abs(float(raw[-1]))
                    name = 'welflx_{0}_{1}_{2}_{3}'.format(kij[0],kij[1],kij[2],kper)
                    new_line = ' '.join(raw[:3]) + " ~  {0}   ~\n".format(name)
                    f.write(new_line)
                    names.append(name)
                    vals.append(val)
    df = pd.DataFrame({"parnme":names,"parval1":vals},index=names)
    df.loc[:,"pargp"] = df.parnme.apply(lambda x : "welflux_{0}".format(x.split('_')[-1]))
    df.loc[:,"scale"] = -1
    df.loc[:,"parubnd"] = df.parval1 + (df.parval1.apply(np.abs) * 0.75)
    df.loc[:, "parlbnd"] = df.parval1 - (df.parval1.apply(np.abs) * 0.75)
    dfs = [df]

    # write a recharge template
    f_in = open(os.path.join(ws,"freyberg6.rch"),'r')
    names,vals = [],[]
    with open(os.path.join(ws,"freyberg6.rch.tpl"),'w') as f:
        f.write("ptf ~\n")
        while True:
            line = f_in.readline()
            if line == '':
                break
            f.write(line)
            if "BEGIN PERIOD" in line:
                kper = int(line.strip().split()[-1]) - 1
                f.write(f_in.readline())
                line = f_in.readline()
                val = float(line.strip().split()[-1])
                name = "rch_{0}".format(kper)
                f.write("   constant ~     {0}     ~\n".format(name))
                names.append(name)
                vals.append(val)
    df = pd.DataFrame({"parnme":names,"parval1":vals},index=names)
    df.loc[:,"pargp"] = "rch"
    df.loc[:,"parubnd"] = df.parval1 * 1.5
    df.loc[:,"parlbnd"] = df.parval1 * 0.5
    df.loc[:, "scale"] = 1.0
    dfs.append(df)

    #  write array some templates
    def write_array_template(filename,prefix):
        arr = np.loadtxt(filename).reshape(nrow,ncol)
        f_tpl = open(filename+".tpl",'w')
        f_tpl.write("ptf ~\n")
        names,vals = [],[]
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                name = "{0}_{1:03d}_{2:03d}".format(prefix,i,j)
                if arr[i,j] == 0:
                    f_tpl.write("  0   ")
                else:
                    f_tpl.write(" ~     {0}      ~ ".format(name))
                    names.append(name)
                    vals.append(arr[i,j])
            f_tpl.write("\n")
        return pd.DataFrame({"parnme":names,"parval1":vals,"pargp":prefix},index=names)
    arr_files = [f for f in os.listdir(ws) if ("npf" in f or "sto" in f) and f.endswith(".dat")]

    for arr_file in arr_files:
        df = write_array_template(os.path.join(ws,arr_file),arr_file.split('.')[0])
        if df.shape[0] == 0:
            continue
        df.loc[:, "parubnd"] = df.parval1 * 10.
        df.loc[:, "parlbnd"] = df.parval1 * 0.1
        df.loc[:, "scale"] = 1.0
        dfs.append(df)
    df = pd.concat(dfs)
    return df


def run_prior_sweep():
    t_d = "template"
    assert os.path.exists(t_d)
    pst_file = "freyberg6.pst"
    assert os.path.exists(os.path.join(t_d, pst_file))
    pst = pyemu.Pst(os.path.join(t_d, pst_file))
    pst.control_data.noptmax = -1
    pst.pestpp_options["ies_par_en"] = "prior.jcb"
    pst.pestpp_options.pop("ies_num_reals",None)
    pst_file = "freyberg6_sweep.pst"
    pst.write(os.path.join(t_d,pst_file),version=2)
    m_d = "master_prior"
    pyemu.os_utils.start_workers(t_d,"pestpp-ies",pst_file,num_workers=15,master_dir=m_d)

def set_truth_obs():
    t_d = "template"
    m_d = "master_prior"
    assert os.path.exists(m_d)
    pst = pyemu.Pst(os.path.join(m_d,"freyberg6_sweep.pst"))
    oe = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(m_d,"freyberg6_sweep.0.obs.csv"))
    pe = pyemu.ParameterEnsemble.from_csv(pst=pst,filename=os.path.join(m_d,"freyberg6_sweep.0.par.csv"))

    pv = oe.phi_vector
    oe.sort_values(by=forecasts[1],inplace=True)
    idx = oe.index[-int(oe.shape[0]/20)]
    print(idx)
    truth_name = os.path.join("template","truth.pst")
    if os.path.exists(truth_name):
        os.remove(truth_name)
    try:
        plot_par_vector(pe.loc[str(idx),pst.par_names],"truth.pdf")
    except:
        plot_par_vector(pe.loc[int(idx), pst.par_names], "truth.pdf")

    pst.observation_data.loc[:,"obsval"] = oe.loc[idx,pst.obs_names]
    pst.observation_data.loc[:,"weight"] = 0.0
    obs = pst.observation_data
    obs.loc[obs.obsnme.apply(lambda x: "2016" in x and ("trgw_2_33_7" in x or "trgw_2_2_9" in x)),"weight"] = 15.0

    g_obs = obs.loc[obs.obsnme.apply(lambda x: "gage_1" in x and "2016" in x),"obsnme"]
    obs.loc[g_obs, "weight"] = 1.0 / (obs.loc[g_obs,"obsval"] * 0.15)
    pst.control_data.noptmax = 0

    par = pst.parameter_data
    par.loc[par.parnme.apply(lambda x: "welflx" in x),"parval1"] *= 0.9
    par.loc[par.parnme.apply(lambda x: "rch" in x), "parval1"] *= 1.05
    par.loc[par.parnme.apply(lambda x: "npf_k_" in x), "parval1"] *= 0.85
    #par.loc[par.parnme.apply(lambda x: "sto" in x), "parval1"] *= 1.2

    #par.loc[par.parnme.apply(lambda x: "npf_k33_" in x), "parval1"] *= 1.1

    pst.write(os.path.join(t_d,"freyberg6_run.pst"),version=2)

    pyemu.os_utils.run("pestpp-ies.exe freyberg6_run.pst",cwd=t_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run.pst"))

    print(pst.phi_components)

def run_ies_hypoth_demo():
    t_d = "template"
    assert os.path.exists(t_d)
    pst_file = "freyberg6_run.pst"

    org_m_d = "master_ies"

    assert os.path.exists(org_m_d)
    if not os.path.exists(os.path.join(org_m_d,"freyberg6_run_ies.base.rei")):
        shutil.copy2(os.path.join(org_m_d,"freyberg6_run.base.rei"),os.path.join(org_m_d,"freyberg6_run_ies.base.rei"))
    pst = pyemu.Pst(os.path.join(org_m_d, "freyberg6_run_ies.pst"))

    # use the final prev ensembles as the initial ensembles this time
    final_pe = "freyberg6_run_ies.{0}.par.csv".format(pst.control_data.noptmax)
    final_oe = final_pe.replace(".par",".obs")
    print(final_pe,final_oe)
    assert os.path.exists(os.path.join(org_m_d,final_pe))
    assert os.path.exists(os.path.join(org_m_d,final_oe))
    
    shutil.copy2(os.path.join(org_m_d,final_pe),os.path.join(t_d,"hypoth_par.csv"))
    pst.pestpp_options["ies_par_en"] = "hypoth_par.csv"

    shutil.copy2(os.path.join(org_m_d, final_oe), os.path.join(t_d, "hypoth_obs.csv"))
    pst.pestpp_options["ies_obs_en"] = "hypoth_obs.csv"
    # now mod the obs en and the control file for the obs we want to stress on
    hp_obs = "out_sfr_flx_20170930"
    obs = pst.observation_data
    assert hp_obs in obs.obsnme
    assert obs.loc[hp_obs,"weight"] == 0.0
    oe = pd.read_csv(os.path.join(t_d,"hypoth_obs.csv"),index_col=0)
    oe = pyemu.ObservationEnsemble(pst=pst,df=oe)
    oe.loc[:,hp_obs] = -200.0
    oe.to_csv(os.path.join(t_d,"hypoth_obs.csv"))
    #obs.loc[hp_obs,"weight"] = 1.0
    omn = oe.mean(axis=0).values
    #print(oe.shape,omn.shape,pst.nobs,pst.res.shape)
    pst.res.loc[pst.obs_names,"modelled"] = omn
    print(pst.phi_components)
    pst.adjust_weights(obs_dict={hp_obs:oe.phi_vector.mean() * 10.0})
    print(pst.phi_components)
    pst.pestpp_options["ies_add_base"] = False


    loc = pyemu.Matrix.from_binary(os.path.join(t_d,"temporal_loc.jcb")).to_dataframe()
    loc.loc[hp_obs,:] = 1.0
    pyemu.Matrix.from_dataframe(df=loc).to_coo(os.path.join(t_d,"hypoth_loc.jcb"))

    pst.control_data.noptmax = 5
    #pst.pestpp_options["ies_lambda_mults"] = 10000.0
    #pst.pestpp_options["lambda_scale_fac"] = 1.0
    #pst.pestpp_options["ies_subset_size"] = 50
    pst.pestpp_options["ies_lambda_dec_fac"] = 1.0
    pst.pestpp_options["ies_autoadaloc"] = False
    pst.pestpp_options.pop("ies_localizer",None)
    #pst.pestpp_options["ies_localizer"] = "hypoth_loc.jcb"
    pst.write(os.path.join(t_d, "freyberg6_run_ies_hypoth.pst"), version=2)
    m_d = "master_ies_hypoth_low"
    pyemu.os_utils.start_workers(t_d, "pestpp-ies", "freyberg6_run_ies_hypoth.pst", num_workers=15,
                                 master_dir=m_d,port=4005)
    return
    oe = pd.read_csv(os.path.join(t_d, "hypoth_obs.csv"), index_col=0)
    oe = pyemu.ObservationEnsemble(pst=pst, df=oe)
    oe.loc[:, hp_obs] = 1500
    oe.to_csv(os.path.join(t_d, "hypoth_obs.csv"))
    pst.observation_data.loc[hp_obs,"weight"] = 0.0
    print(pst.phi_components)
    pst.adjust_weights(obs_dict={hp_obs:oe.phi_vector.mean() * 2.5})
    print(pst.phi_components)
    m_d = "master_ies_hypoth_high"
    pyemu.os_utils.start_workers(t_d, "pestpp-ies", "freyberg6_run_ies_hypoth.pst", num_workers=15,
                                 master_dir=m_d, port=4005)

def plot_ies_hypoth_results():
    m_d_h = "master_ies_hypoth_high"
    m_d_l = "master_ies_hypoth_low"
    assert os.path.exists(m_d_h)
    assert os.path.exists(m_d_l)
    pst = pyemu.Pst(os.path.join(m_d_h,"freyberg6_run_ies_hypoth.pst"))
    onames = [n for n in pst.nnz_obs_names if "sfr" not in n]

    oei_l = pd.read_csv(os.path.join(m_d_l, "freyberg6_run_ies_hypoth.0.obs.csv"), index_col=0)
    oef_l = pd.read_csv(os.path.join(m_d_l, "freyberg6_run_ies_hypoth.{0}.obs.csv".format(pst.control_data.noptmax)),
                        index_col=0)
    oei_l = pyemu.ObservationEnsemble(pst=pst, df=oei_l)
    oef_l = pyemu.ObservationEnsemble(pst=pst, df=oef_l)

    oei_h = pd.read_csv(os.path.join(m_d_h,"freyberg6_run_ies_hypoth.0.obs.csv"),index_col=0)
    oef_h = pd.read_csv(os.path.join(m_d_h, "freyberg6_run_ies_hypoth.{0}.obs.csv".format(pst.control_data.noptmax)), index_col=0)
    oei_h = pyemu.ObservationEnsemble(pst=pst,df=oei_h)
    oef_h = pyemu.ObservationEnsemble(pst=pst, df=oef_h)

    # oei = pd.read_csv(os.path.join(m_d, "freyberg6_run_ies_hypoth.0.obs.csv"), index_col=0)
    # oef = pd.read_csv(os.path.join(m_d, "freyberg6_run_ies_hypoth.{0}.obs.csv".format(pst.control_data.noptmax)),
    #                   index_col=0)
    # oei = pyemu.ObservationEnsemble(pst=pst, df=oei)
    # oef = pyemu.ObservationEnsemble(pst=pst, df=oef)

    hp_obs = [n for n in pst.nnz_obs_names if "sfr" in n][0]
    fig,axes = plt.subplots(1,2,figsize=(10,5))
    oei_l.loc[:,hp_obs].hist(ax=axes[0],fc="0.5",ec="none",alpha=0.5,label="posterior")
    oef_l.loc[:,hp_obs].hist(ax=axes[0],fc="b",ec="none",alpha=0.5,label="hypothesis low")
    oef_h.loc[:, hp_obs].hist(ax=axes[0], fc="g", ec="none", alpha=0.5, label="hypothesis high")
    oei_l.loc[:,onames].phi_vector.hist(ax=axes[1],fc="0.5",ec="none",alpha=0.5,label="posterior")
    oef_l.loc[:, onames].phi_vector.hist(ax=axes[1], fc="b", ec="none", alpha=0.5,label="hypothesis low")
    oef_h.loc[:, onames].phi_vector.hist(ax=axes[1], fc="g", ec="none", alpha=0.5, label="hypothesis high")

    axes[0].legend()
    axes[1].legend()
    axes[0].set_title("dry season gw-to-sw discharge")
    axes[1].set_title("objective function")

    plt.savefig("hypoth.pdf")




    # covi = np.dot(oei.values.T,oei.values)
    # print(covi.shape,np.linalg.det(covi))
    # covf = np.dot(oef.values.T,oef.values)
    # print(covf.shape,np.linalg.det(covf))
    # print(covf)
    # diff = oef - oei
    # print(diff.mean())
    # dcov = np.dot(diff.values.T,diff.values)
    # plt.imshow(dcov)
    # plt.show()

def run_ies_demo():
    '''todo: localize
    '''
    t_d = "template"
    assert os.path.exists(t_d)
    pst_file = "freyberg6_run.pst"
    assert os.path.exists(os.path.join(t_d,pst_file))
    pst = pyemu.Pst(os.path.join(t_d,pst_file))
    pst.control_data.noptmax = 4
    pe = pyemu.ParameterEnsemble.from_binary(pst=pst,filename=os.path.join(t_d,"prior.jcb"))
    pe = pe.iloc[:50,:]
    pe.to_binary(os.path.join(t_d,"ies_prior.jcb"))
    pst.pestpp_options = {}
    pst.pestpp_options["ies_par_en"] = "ies_prior.jcb"
    #pst.pestpp_options["ies_no_noise"] = True

    sim = flopy.mf6.MFSimulation.load(sim_ws="template")
    m = sim.get_model("freyberg6")
    par = pst.parameter_data
    tpar = par.loc[par.parnme.apply(lambda x: x.startswith("wel") or x.startswith("rch")),:].copy()
    tpar.loc[:,"kper"] = tpar.parnme.apply(lambda x: int(x.split('_')[-1]))
    perlen = sim.tdis.perioddata.array["perlen"]
    totim = np.cumsum(perlen).astype(int)
    sp_start = pd.to_datetime("12-31-2015") + pd.to_timedelta(totim, unit='d')
    tpar.loc[:, "datetime"] = tpar.kper.apply(lambda x: sp_start[x])
    tgrps = set(tpar.pargp.unique())
    obs = pst.observation_data.loc[pst.nnz_obs_names,:].copy()
    obs.loc[:,"datetime"] = pd.to_datetime(obs.obsnme.apply(lambda x: x.split('_')[-1]))
    ogrps = [g for g in pst.par_groups if g not in tgrps]
    cols = tpar.parnme.tolist()
    cols.extend(ogrps)
    print(cols)
    rows = pst.nnz_obs_names
    loc = pd.DataFrame(index=rows,columns=cols)
    loc.loc[:,:] = 1.0
    tdelt_min = pd.to_timedelta(-90,unit='d')
    tdelt_max = pd.to_timedelta(2,unit='d')
    for oname in pst.nnz_obs_names[:5]:
        odt = obs.loc[oname,"datetime"]
        tdelt = (tpar.datetime - odt)
        too_back = tdelt.loc[tdelt < tdelt_min]
        too_forward = tdelt.loc[tdelt > tdelt_max]
        print(oname,too_back,too_forward)
        loc.loc[oname,too_back.index] = 0.0
        loc.loc[oname,too_forward.index] = 0.0
       #break
    pyemu.Matrix.from_dataframe(df=loc).to_coo(os.path.join(t_d,"temporal_loc.jcb"))

    pst.write(os.path.join(t_d,"freyberg6_run_ies.pst"),version=2)
    m_d = "master_ies_default"
    pyemu.os_utils.start_workers(t_d, "pestpp-ies", "freyberg6_run_ies.pst", num_workers=15, master_dir=m_d)

    pst.pestpp_options["ies_localizer"] = "temporal_loc.jcb"
    pst.pestpp_options["ies_autoadaloc"] = True
    pst.pestpp_options["ies_num_threads"] = 3
    #pst.pestpp_options["ies_no_noise"] = True
    pst.pestpp_options["overdue_giveup_fac"] = 10.0
    pst.write(os.path.join(t_d, "freyberg6_run_ies.pst"), version=2)
    m_d = "master_ies"
    pyemu.os_utils.start_workers(t_d, "pestpp-ies", "freyberg6_run_ies.pst", num_workers=15,
                                master_dir=m_d)


    pst = pyemu.Pst(os.path.join(t_d, pst_file))
    pst.control_data.noptmax = 3
    pst.pestpp_options = {"overdue_giveup_fac":10}
    _block_tie(pst)
    pst.write(os.path.join(t_d,"freyberg6_run_ies.pst"))

    m_d = "master_ies_default"
    pyemu.os_utils.start_workers(t_d, "pestpp-ies", "freyberg6_run_ies.pst", num_workers=15, master_dir=m_d)

def make_ies_figs(m_d="master_ies",plt_case="ies"):
    assert os.path.join(m_d)
    pst_file = "freyberg6_run_ies.pst"
    pst = pyemu.Pst(os.path.join(m_d,pst_file))
    pr_oe = pd.read_csv(os.path.join(m_d,pst_file.replace(".pst",".0.obs.csv")),index_col=0)
    pt_oe = pd.read_csv(os.path.join(m_d, pst_file.replace(".pst", ".{0}.obs.csv".format(pst.control_data.noptmax))),index_col=0)
    pr_oe = pyemu.ObservationEnsemble(pst=pst,df=pr_oe)
    pt_oe = pyemu.ObservationEnsemble(pst=pst,df=pt_oe)
    pr_pv = pr_oe.phi_vector
    pt_pv = pt_oe.phi_vector
    #keep = pt_pv.loc[pt_pv<50].index
    #pt_oe = pt_oe.loc[keep,:]
    pt_pv = pt_oe.phi_vector
    pr_pe = pd.read_csv(os.path.join(m_d, pst_file.replace(".pst", ".0.par.csv")),index_col=0)
    pt_pe = pd.read_csv(os.path.join(m_d, pst_file.replace(".pst", ".{0}.par.csv".format(pst.control_data.noptmax))),index_col=0)
    print(pt_pv.sort_values())

    full_pr_oe = pd.read_csv(os.path.join("master_prior","freyberg6_sweep.0.obs.csv"),index_col=0)

    for real in ["base",pt_pe.index[0]]:
        try:
            phi = pt_pv.loc[str(real)]
        except:
            phi = pt_pv.loc[int(real)]
        #plot_par_vector(pr_pe.loc[real],"{0}_pr_{1}_phi_{2:3.1E}.pdf".format(plt_case,real,phi))
        plot_par_vector(pt_pe.loc[real], "{0}_pt_{1}_phi_{2:3.2E}.pdf".format(plt_case,real,phi))





    obs = pst.observation_data
    print(pst.nnz_obs_groups)
    #print(pst.forecast_names)
    fig = plt.figure(figsize=(8,6))
    ax_count = 0

    for i,nz_grp in enumerate(pst.nnz_obs_groups):
        grp_obs = obs.loc[obs.obgnme==nz_grp,:].copy()
        x = np.arange(1, grp_obs.shape[0] + 1)
        print(grp_obs)
        grp_obs.loc[:,"datetime"] = pd.to_datetime(grp_obs.obsnme.apply(lambda x: x.split('_')[-1]))
        ax = plt.subplot2grid((5,4),(i,0),colspan=4)
        [ax.plot(x,pr_oe.loc[i,grp_obs.obsnme],'0.5',lw=0.1, alpha=0.25) for i in pr_oe.index]
        for i in pt_oe.index:
            vals = pt_oe.loc[i,grp_obs.obsnme].values
            vals[vals<20] = np.NaN
            ax.plot(x, vals, 'b', lw=0.1, alpha=0.35)
        #[ax.plot(x,pt_oe.loc[i,grp_obs.obsnme],'b',lw=0.1,alpha=0.35) for i in pt_oe.index]
        grp_obs = grp_obs.loc[grp_obs.weight > 0,:]
        x = np.arange(2, grp_obs.shape[0]+2)
        ax.plot(x,grp_obs.obsval, 'r',lw=1.5)
        full_upper = full_pr_oe.loc[:,grp_obs.obsnme].max().max()
        full_lower = full_pr_oe.loc[:,grp_obs.obsnme].min().min()
        ax.set_ylim(full_lower,full_upper)
        unit = None
        label = None
        for tag,u in unit_dict.items():
            if tag in nz_grp:
                unit = u
        for tag,l in label_dict.items():
            if tag in nz_grp:
                label = l + " observed vs simulated"
        ax.set_title("{0}) {1}".format(abet[ax_count], label), loc="left")
        ax.set_ylabel(unit)

        ax_count += 1

    ax = plt.subplot2grid((5,4),(3,0),rowspan=2)
    pr_pv.loc[pr_pv <= 0.0] = np.NaN
    pt_pv.loc[pt_pv <= 0.0] = np.NaN

    ax.hist(pr_pv.apply(np.log10),alpha=0.5,facecolor="0.5",edgecolor="none")
    ax.hist(pt_pv.apply(np.log10),alpha=0.5,facecolor="b",edgecolor="none")
    ax.set_title("{0}) objective function".format(abet[ax_count]), loc="left")
    #ax.set_yticks([])
    ax.set_xlabel("$log_{10}$ objective function")
    ax.set_ylabel("number of realizations")
    ax_count += 1

    for i,forecast in enumerate(forecasts):
        ax = plt.subplot2grid((5,4),(3,i+1),rowspan=2)
        pr_oe.loc[:,forecast].hist(ax=ax,density=True,facecolor='0.5',edgecolor="none",alpha=0.5)
        pt_oe.loc[:,forecast].hist(ax=ax,density=True,facecolor='b',edgecolor="none",alpha=0.5)
        ylim = ax.get_ylim()
        oval = obs.loc[forecast,"obsval"]
        ax.plot([oval,oval],ylim,'r', lw=1.5)

        unit,label = None,None
        for tag,u in unit_dict.items():
            if tag in forecast:
                unit = u
        for tag,l in label_dict.items():
            if tag in forecast:
                label = l + " forecast"
        ax.set_title("{0}) {1}".format(abet[ax_count], label), loc="left")
        ax.set_xlabel(unit)
        ax.set_ylabel("increasing probability density")
        ax.set_yticks([])
        ax.grid(False)
        ax_count +=1
        full_upper = full_pr_oe.loc[:, forecast].max()
        full_lower = full_pr_oe.loc[:, forecast].min()
        ax.set_xlim(full_lower,full_upper)

    plt.tight_layout()
    plt.savefig(os.path.join(plt_dir,"{0}.pdf".format(plt_case)))


def make_glm_figs(m_d="master_glm",plt_case = "glm"):
    assert os.path.join(m_d)
    pst_file = "freyberg6_run_glm.pst"
    pst = pyemu.Pst(os.path.join(m_d,pst_file))
    pst.parrep(os.path.join(m_d,pst_file.replace(".pst",".par{0}".format(pst.control_data.noptmax))))

    pt_oe = None
    try:
        pt_oe = pd.read_csv(os.path.join(m_d,pst_file.replace(".pst",".post.obsen.csv")),index_col=0)
        pt_oe = pyemu.ObservationEnsemble(pst=pst,df=pt_oe)
    except:
        pass

    f_df = pd.read_csv(os.path.join(m_d,pst_file.replace(".pst",".pred.usum.csv")),index_col=0)
    f_df.index = f_df.index.map(str.lower)


    if pt_oe is not None:
        pv = pt_oe.phi_vector
        pv.sort_values(inplace=True)
        keep = pv.iloc[:50].index
        #keep = pv.loc[pv < max(60, pst.phi * 2.0)].index
        pt_oe = pt_oe.loc[keep, :]
        pt_pe = pd.read_csv(os.path.join(m_d, pst_file.replace(".pst", ".post.paren.csv")), index_col=0)
        pt_pe = pt_pe.loc[keep,:]
        pv = pt_oe.phi_vector
        for real in [pt_pe.index[0]]:
            phi = pv.loc[real]
            plot_par_vector(pt_pe.loc[real], "{0}_pt_{1}_phi_{2:3.2E}.pdf".format(plt_case,real,phi))
        print(pt_oe.shape,pst.phi)
        print(pv.min(),pv.max())
        return

    plot_par_vector(pst.parameter_data.parval1.copy(), "{0}_pt_base_{1:3.2E}.pdf".format(plt_case,pst.phi))
    obs = pst.observation_data
    #print(pst.nnz_obs_groups)
    #print(pst.forecast_names)

    fig = plt.figure(figsize=(8,6))
    ax_count = 0

    full_pr_oe = pd.read_csv(os.path.join("master_prior", "freyberg6_sweep.0.obs.csv"), index_col=0)

    for i,nz_grp in enumerate(pst.nnz_obs_groups):
        grp_obs = obs.loc[obs.obgnme==nz_grp,:].copy()
        x = np.arange(1, grp_obs.shape[0]+1)
        #print(grp_obs)
        grp_obs.loc[:,"datetime"] = pd.to_datetime(grp_obs.obsnme.apply(lambda x: x.split('_')[-1]))
        grp_obs.sort_values(by="datetime")
        ax = plt.subplot2grid((5,4),(i,0),colspan=4)
        #ax.plot(grp_obs.datetime,grp_obs.obsval, 'r')
        if pt_oe is not None:
            grp_oe = pt_oe.loc[:,grp_obs.obsnme].copy()
            if "trgw" in nz_grp:
                grp_oe.values[grp_oe.values<20] = np.NaN
            [ax.plot(x,grp_oe.loc[i,grp_obs.obsnme].values,'b',lw=0.1,alpha=0.5) for i in grp_oe.index]
        else:
            ax.plot(x,pst.res.loc[grp_obs.obsnme,"modelled"],'b',lw=1.5)
        grp_obs = grp_obs.loc[grp_obs.weight >0,:]
        x = np.arange(2, grp_obs.shape[0] + 2)
        ax.plot(x,grp_obs.obsval.values, 'r',lw=1.5)
        full_upper = full_pr_oe.loc[:, grp_obs.obsnme].max().max()
        full_lower = full_pr_oe.loc[:, grp_obs.obsnme].min().min()
        ax.set_ylim(full_lower, full_upper)
        unit = None
        label = None
        for tag, u in unit_dict.items():
            if tag in nz_grp:
                unit = u
        for tag, l in label_dict.items():
            if tag in nz_grp:
                label = l + " observed vs simulated"
        ax.set_title("{0}) {1}".format(abet[ax_count], label), loc="left")
        ax.set_ylabel(unit)
        ax.grid(False)

        ax_count += 1

    if pt_oe is not None:
        ax = plt.subplot2grid((5,4),(3,0),rowspan=2)
        ax.hist(pv.apply(np.log10),alpha=0.5,facecolor="b",edgecolor="none")
        ax.set_title("{0}) objective function".format(abet[ax_count]), loc="left")
        #ax.set_yticks([])
        ax.set_xlabel("$log_{10}$ objective function")
        ax.set_ylabel("number of realizations")
        ax_count += 1

    for i,forecast in enumerate(forecasts):
        ax = plt.subplot2grid((5,4),(3,i+1),rowspan=2)
        axt = plt.twinx()
        x,y = pyemu.plot_utils.gaussian_distribution(f_df.loc[forecast,"prior_mean"],f_df.loc[forecast,"prior_stdev"])
        axt.fill_between(x,0,y,facecolor="0.5",alpha=0.25)
        x,y = pyemu.plot_utils.gaussian_distribution(f_df.loc[forecast,"post_mean"],f_df.loc[forecast,"post_stdev"])
        axt.fill_between(x,0,y,facecolor="b",alpha=0.25)
        if pt_oe is not None:
            pt_oe.loc[:,forecast].hist(ax=ax,density=True,facecolor='b',edgecolor="none",alpha=0.5)
        ylim = ax.get_ylim()
        oval = obs.loc[forecast,"obsval"]
        ax.plot([oval,oval],ylim,'r',lw=1.5)

        ax.set_ylim(ylim)
        unit, label = None, None
        for tag, u in unit_dict.items():
            if tag in forecast:
                unit = u
        for tag, l in label_dict.items():
            if tag in forecast:
                label = l + " forecast"
        ax.set_title("{0}) {1}".format(abet[ax_count], label), loc="left")
        ax.set_xlabel(unit)
        ax.set_ylabel("increasing probability density")
        ax.set_yticks([])
        full_upper = full_pr_oe.loc[:, forecast].max()
        full_lower = full_pr_oe.loc[:, forecast].min()
        ax.set_xlim(full_lower, full_upper)
        ax.grid(False)
        ylim = axt.get_ylim()
        axt.set_ylim(0,ylim[1])
        axt.set_yticks([])


        ax_count += 1


    plt.tight_layout()
    plt.savefig(os.path.join(plt_dir,"{0}.pdf".format(plt_case)))

def invest():
    pst = pyemu.Pst(os.path.join("template","freyberg6_run.pst"))
    oe = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join("master_prior","freyberg6_sweep.0.obs.csv"))
    pv = oe.phi_vector
    pv.loc[pv<0.01] = np.NaN
    print(pv.min(),pv.max(),pv.mean())
    return
    obs = pst.observation_data.loc[pst.nnz_obs_names,:].copy()
    obs = obs.loc[obs.obgnme=="gage",:]
    #obs.loc[:,"weight"] = 0.005
    obs.loc[:,"cv"] = (1./ obs.weight) / obs.obsval
    print(obs)

def _block_tie(pst,t_d="template"):
    sim = flopy.mf6.MFSimulation.load(sim_ws=t_d)
    m = sim.get_model("freyberg6")
    xgrid = m.modelgrid.xcellcenters
    ygrid = m.modelgrid.ycellcenters

    # fix/tie to reduce adj par numbers
    tie_step = 6
    par = pst.parameter_data
    par.loc[:, "partied"] = np.NaN
    gr_par = par.loc[par.pargp.apply(lambda x: "npf" in x or "sto" in x), :].copy()
    gr_par.loc[:, "k"] = gr_par.parnme.apply(lambda x: int(x.split('_')[2]))
    gr_par.loc[:, "i"] = gr_par.parnme.apply(lambda x: int(x.split('_')[3]))
    gr_par.loc[:, "j"] = gr_par.parnme.apply(lambda x: int(x.split('_')[4]))

    for i in range(tie_step, gr_par.i.max() + tie_step, tie_step):
        for j in range(tie_step, gr_par.j.max() + tie_step, tie_step):
            ir = [ii for ii in range(i - tie_step, i)]
            jr = [jj for jj in range(j - tie_step, j)]
            t_par = gr_par.loc[gr_par.apply(lambda x: x.i in ir and x.j in jr, axis=1), :]
            for grp in t_par.pargp.unique():
                tt_par = t_par.loc[t_par.pargp == grp, :]
                par.loc[tt_par.parnme[1:], "partrans"] = "tied"
                par.loc[tt_par.parnme[1:], "partied"] = tt_par.parnme[0]
    w_par = par.loc[par.parnme.str.startswith("wel"), :]
    for grp in w_par.pargp.unique():
        ww_par = w_par.loc[w_par.pargp == grp, :]
        par.loc[ww_par.parnme[1:], "partrans"] = "tied"
        par.loc[ww_par.parnme[1:], "partied"] = ww_par.parnme[0]
    spars = [n for n in pst.adj_par_names if not "rch" in n and not "wel" in n]
    dist_df = pd.DataFrame({"parnme": spars, "x": np.nan, "y": np.nan}, index=spars)
    dist_df.loc[:, "i"] = gr_par.loc[dist_df.parnme, "i"]
    dist_df.loc[:, "j"] = gr_par.loc[dist_df.parnme, "j"]

    dist_df.loc[:, "x"] = dist_df.apply(lambda x: xgrid[x.i, x.j], axis=1)
    dist_df.loc[:, "y"] = dist_df.apply(lambda x: ygrid[x.i, x.j], axis=1)
    dist_df.loc[:, "layer_group"] = dist_df.parnme.apply(lambda x: "_".join(x.split('_')[:3]))
    print(dist_df.layer_group.unique())
    dfs = []
    for lg in dist_df.layer_group.unique():
        dfs.append(dist_df.loc[dist_df.layer_group == lg, :].copy())
    v = pyemu.geostats.ExpVario(contribution=1.0, a=1000)
    gs = pyemu.geostats.GeoStruct(variograms=v)
    tpars = [n for n in pst.adj_par_names if "rch" in n or "wel" in n]
    tpar = par.loc[tpars, :].copy()
    tpar.loc[:, "x"] = tpar.parnme.apply(lambda x: int(x.split('_')[-1]))
    tpar.loc[:, "y"] = 0.0
    tpar.loc[:, "group"] = tpar.parnme.apply(lambda x: x.split('_')[0])
    print(tpar.group.unique())
    tdfs = []
    for g in tpar.group.unique():
        tdfs.append(tpar.loc[tpar.group == g, :].copy())

    tv = pyemu.geostats.ExpVario(contribution=1.0, a=3)
    tgs = pyemu.geostats.GeoStruct(variograms=v)

    cov = pyemu.helpers.geostatistical_prior_builder(pst=pst, struct_dict={gs: dfs, tgs: tdfs})
    cov.to_ascii(os.path.join(t_d, "glm_prior.cov"))


def run_glm_demo():

    t_d = "template"
    assert os.path.exists(t_d)
    pst_file = "freyberg6_run.pst"
    assert os.path.exists(os.path.join(t_d,pst_file))
    pst = pyemu.Pst(os.path.join(t_d,pst_file))

    _block_tie(pst)
    pst.control_data.noptmax = 4
    pst.pestpp_options["forecasts"] = forecasts
    pst.write(os.path.join(t_d,"freyberg6_run_glm.pst"),version=2)
    m_d = "master_glm_default"
    pyemu.os_utils.start_workers(t_d, "pestpp-glm", "freyberg6_run_glm.pst", num_workers=15, master_dir=m_d)

    pst.pestpp_options["n_iter_super"] = pst.control_data.noptmax
    pst.pestpp_options["n_iter_base"] = -1
    pst.pestpp_options["glm_num_reals"] = 200
    #pst.pestpp_options["parcov"] = "glm_prior.cov"
    pst.pestpp_options["glm_normal_form"] = "prior"
    pst.pestpp_options["max_n_super"] = int(pst.nnz_obs)
    #pst.svd_data.maxsing =  int(pst.nnz_obs)
    pst.write(os.path.join(t_d, "freyberg6_run_glm.pst"), version=2)
    m_d = "master_glm"
    pyemu.os_utils.start_workers(t_d, "pestpp-glm", "freyberg6_run_glm.pst", num_workers=15, master_dir=m_d)


def run_sen_demo():
    t_d = "template"
    assert os.path.exists(t_d)
    pst_file = "freyberg6_run.pst"
    assert os.path.exists(os.path.join(t_d,pst_file))
    pst = pyemu.Pst(os.path.join(t_d,pst_file))
    par = pst.parameter_data
    w_names = par.loc[par.parnme.str.startswith("wel"),"parnme"]
    par.loc[w_names,"pargp"] = "welflux"

    pst.pestpp_options["tie_by_group"] = True
    pst.write(os.path.join(t_d,"freyberg6_run_sen.pst"),version=2)
    m_d = "master_sen_morris"
    pyemu.os_utils.start_workers(t_d, "pestpp-sen", "freyberg6_run_sen.pst", num_workers=15, master_dir=m_d)

def start():
    pyemu.os_utils.start_workers("template", "pestpp-opt", "freyberg6_run_opt.pst", num_workers=25)

def make_sen_figs():
    m_d = "master_sen_morris"
    assert os.path.exists(m_d)
    pst_file = "freyberg6_run_sen.pst"
    pst = pyemu.Pst(os.path.join(m_d,pst_file))
    mio_df = pd.read_csv(os.path.join(m_d,pst_file.replace(".pst",".mio")),index_col=0)
    phi_df = pd.read_csv(os.path.join(m_d,pst_file.replace(".pst",".msn")),skipinitialspace=True).loc[:,["parameter_name","sen_mean_abs", "sen_std_dev"]]

    df = mio_df.loc[mio_df.index.map(lambda x: x in forecasts),["parameter_name","sen_mean_abs", "sen_std_dev"]].copy()
    phi_df.index = ["phi" for _ in range(phi_df.shape[0])]
    df = pd.concat([df,phi_df])
    #print(df)
    #print(phi_df.columns)
    #print(df.columns)
    fig,axes = plt.subplots(2,1,figsize=(5,6))
    x = np.arange(df.parameter_name.unique().shape[0])
    offset = -0.2
    step = 0.1
    name_dict = {"headwater":"headwater forecast","tailwater":"tailwater forecast","trgw":"groundwater forecast","phi":"objective function"}
    par_dict = {"k33_0":"Layer 1 VK","k33_1":"Layer 2 VK","k33_2":"Layer 3 VK",
    "ss_0":"Layer 1 SS","ss_1":"Layer 2 SS","ss_2":"Layer 3 SS",
    "sy":"Layer 1 SY","rch":"Recharge","wel":"Well extraction","k_0":"Layer 1 HK","k_1":"Layer 2 HK","k_2":"Layer 3 HK"}
    df.index = df.index.map(lambda x: name_dict[x.split('_')[0]])
    df.loc[:,"parameter_name"] = df.parameter_name.apply(lambda x: [v for k,v in par_dict.items() if k in x][0])
    ax = axes[0]
    for o in df.index.unique():
        print(o)
        odf = df.loc[o,:].copy()
        odf.loc[:,"scaled"] = odf.sen_mean_abs / odf.sen_mean_abs.sum()
        odf.loc[:,"scaled"] *= 100.0
        odf.sort_values(by="parameter_name",inplace=True)

        ax.bar(x+offset,odf.scaled,width=step,label=o)
        offset += step
    ax.set_xticks(x)
    ax.set_xticklabels(odf.parameter_name.values,rotation=90)
    spnspecs.graph_legend(ax=ax,loc="upper left")
    ax.set_ylabel("percent of total mean\nabsolute sensitivity")
    #ax.set_xlabel("parameter")
    ax.set_title("A) summary of mean absolute sensitivity",loc="left")

    ax = axes[1]
    offset = -0.2
    for o in df.index.unique():
        print(o)
        odf = df.loc[o,:].copy()
        odf.loc[:, "sen_std_dev"] /= odf.sen_mean_abs
        odf.sort_values(by="parameter_name",inplace=True)

        ax.bar(x+offset,odf.sen_std_dev,width=step,label=o)
        offset += step
    ax.set_xticks(x)
    ax.set_xticklabels(odf.parameter_name.values,rotation=90)
    #spnspecs.graph_legend(ax=ax,loc="upper left")#bbox_to_anchor=(0.5,-0.75))
    #ax.set_ylabel("normalized standard deviation")
    ax.set_ylabel("coefficient of variation")
    ax.set_xlabel("parameter")
    ax.set_title("B) summary of sensitivity variability",loc="left")

    # ax = axes[2]
    # ax.set_frame_on(False)
    # ax.set_xticks([])
    # ax.set_yticks([])

    plt.tight_layout()
    plt.savefig(os.path.join(plt_dir,"morris.pdf"))


def run_opt_demo():
    t_d = "template"
    assert os.path.exists(t_d)
    pst_file = "freyberg6_run.pst"
    assert os.path.exists(os.path.join(t_d,pst_file))
    pst = pyemu.Pst(os.path.join(t_d,pst_file))
    #_block_tie(pst)

    # process the wel flux pars into dec  vars
    par = pst.parameter_data
    par.loc[:,"partrans"] = "fixed"
    par.loc[par.parnme.apply(lambda x: "rch" in x),"partrans"] = "log"
    w_par = par.loc[par.parnme.str.startswith("wel"),:].copy()
    w_par.loc[:,"parval1"] = 0.0
    w_par.loc[:,"partrans"] = "none"
    w_par.loc[:,"parubnd"] = 500.0
    w_par.loc[:,"parlbnd"] = 0.0
    w_par.loc[:,"kper"] = w_par.parnme.apply(lambda x: int(x.split('_')[-1]))


    w_par.loc[:,"pargp"] = "welflux" #w_par.parnme.apply(lambda x: "_".join(x.split('_')[:-1]))
    w_par.loc[w_par.kper == 0, "pargp"] = "fixed_welpar"
    w_par.loc[w_par.kper == 0, "partrans"] = "fixed"
    w_par.loc[w_par.kper==0,"parval1"] = par.loc[w_par.loc[w_par.kper==0,"parnme"],"parval1"]
    for col in pst.parameter_data.columns:
        pst.parameter_data.loc[w_par.parnme,col] = w_par.loc[:,col]
    pst.rectify_pgroups()
    pst.parameter_groups.loc["welflux","derinc"] = 100.0
    pst.parameter_groups.loc["welflux","inctyp"] = "absolute"
    print(pst.parameter_groups)

    # mod the headwater and tailwater obs
    min_swgw_flux = -250
    obs = pst.observation_data
    #obs.loc[:,"weight"] = 0.0
    constraints = obs.loc[obs.obsnme.apply(lambda x: "headwater" in x or "tailwater" in x),"obsnme"]
    pst.observation_data.loc[constraints,"obgnme"] = "less_than_flux"
    pst.observation_data.loc[constraints,"obsval"] = min_swgw_flux
    pst.observation_data.loc[constraints,"weight"] = 1.0

    # build up pi constraints
    equations,pilbl = [],[]
    min_well_flux = 750
    for kper in w_par.kper.unique():
        kper_par = w_par.loc[w_par.kper==kper,:]
        eq = ""
        for name in kper_par.parnme:
            eq += " 1.0 * {0} +".format(name)
        eq = eq[:-1] + "= {0}".format(min_well_flux)
        equations.append(eq)
        pilbl.append("pi_flx_{0}".format(kper))

    pi_df = pd.DataFrame({"equation":equations,"weight":1.0,"pilbl":pilbl,"obgnme":"greater_than_flux"},index=pilbl)
    pst.prior_information = pi_df
    pst.pestpp_options["opt_dec_var_groups"] = "welflux"
    pst.pestpp_options["opt_direction"] = "max"
    pst.control_data.noptmax = 1
    pst.write(os.path.join(t_d,"freyberg6_run_opt.pst"),version=2)
    m_d = "master_opt_neutral"
    pyemu.os_utils.start_workers(t_d, "pestpp-opt", "freyberg6_run_opt.pst", num_workers=15, master_dir=m_d)

    pst.pestpp_options["opt_risk"] = 0.95
    pst.write(os.path.join(t_d, "freyberg6_run_opt.pst"))
    m_d = "master_opt_fosm"
    pyemu.os_utils.start_workers(t_d, "pestpp-opt", "freyberg6_run_opt.pst", num_workers=15, master_dir=m_d)

    oe_pt = os.path.join("master_ies","freyberg6_run_ies.3.obs.csv")
    if not os.path.exists(oe_pt):
        raise Exception("couldnt find existing oe_pt:"+oe_pt)
    shutil.copy(oe_pt,os.path.join(t_d,"obs_stack.csv"))
    pe_pt = os.path.join("master_ies","freyberg6_run_ies.3.par.csv")
    if not os.path.exists(pe_pt):
        raise Exception("couldnt find existing pe_pt:"+pe_pt)
    shutil.copy(pe_pt,os.path.join(t_d,"par_stack.csv"))

    #pst.pestpp_options["opt_obs_stack"] = "obs_stack.csv"
    par = pst.parameter_data
    par.loc[par.partrans=="tied","partrans"] = "log"
    pst.pestpp_options["opt_par_stack"] = "par_stack.csv"
    #pst.pestpp_options["opt_recalc_chance_every"] = 100
    pst.pestpp_options["opt_risk"] = 0.95
    pst.write(os.path.join(t_d,"freyberg6_run_opt.pst"))
    m_d = "master_opt_stack"
    pyemu.os_utils.start_workers(t_d, "pestpp-opt", "freyberg6_run_opt.pst", num_workers=15, master_dir=m_d)

def _rebase_results():
    #raise Exception("you better be sure!")
    b_d = "baseline_results"
    if os.path.exists(b_d):
        shutil.rmtree(b_d)
    os.mkdir(b_d)
    for m_d,files in test_dict.items():
        assert os.path.exists(m_d),"master dir {0} missing".format(m_d)
        b_m_d = os.path.join(b_d,m_d)
        os.makedirs(b_m_d)
        for f in files:
            assert os.path.exists(os.path.join(m_d,f)),"file {0} missing from master dir {1}".format(f,m_d)
            shutil.copy2(os.path.join(m_d,f),os.path.join(b_m_d,f))

def compare_to_baseline(should_raise=True,tol=0.1):
    b_d = "baseline_results"
    assert os.path.exists(b_d)
    errors = []
    for m_d,files in test_dict.items():
        if not os.path.exists(m_d):
            errors.append("master dir {0} missing".format(m_d))
            continue
        b_m_d = os.path.join(b_d,m_d)
        if not os.path.exists(b_m_d):
            errors.append("base master dir {0} missing".format(b_m_d))
            continue
        for f in files:

            if not os.path.exists(os.path.join(m_d,f)):
                errors.append("file {0} missing from master dir {1}".format(f,m_d))
                continue
            if not os.path.exists(os.path.join(b_m_d,f)):
                errors.append("file {0} missing from base master dir {1}".format(f,b_m_d))
                continue
            if f.lower().endswith(".par"):
                df1 = pyemu.pst_utils.read_parfile(os.path.join(m_d,f))
                df2 = pyemu.pst_utils.read_parfile(os.path.join(b_m_d,f))
            elif f.lower().endswith(".rei"):
                df1 = pyemu.pst_utils.read_resfile(os.path.join(m_d,f))
                df2 = pyemu.pst_utils.read_resfile(os.path.join(b_m_d,f))
            else:
                df1 = pd.read_csv(os.path.join(m_d,f))
                df2 = pd.read_csv(os.path.join(b_m_d,f))

            if df1.shape != df2.shape:
                errors.append("shape {0} vs {1} mismatch for {2}".format(str(df1.shape),str(df2.shape),f))
                continue
            for col in df1.columns:

                if df1.dtypes[col] == object:
                    continue
                d = ((df1.loc[:,col] - df2.loc[:,col])/df1.loc[:,col]).apply(np.abs).sum()
                #print(f,col,d)
                if d > tol:
                    errors.append("col {0} in file {1} too different: {2}".format(col,f,d))
    if len(errors) > 0:
        if should_raise:
            raise Exception("errors in compare: {0}".format("\n".join(errors)))
        else:
            print("differences found: \n"+'\n'.join(errors))

def test(should_raise=False):
    print(os.listdir("."))
    pyemu.os_utils.run("mf6",cwd="template")
    
    run_ies_demo()
    run_glm_demo()
    run_sen_demo()
    run_opt_demo()
    compare_to_baseline(should_raise=should_raise, tol=1.0e+30)

def plot_par_vector(pval_series=None,plt_name="par.pdf"):

    sim = flopy.mf6.MFSimulation.load(sim_ws="template")
    m = sim.get_model("freyberg6")
    ib = m.dis.idomain.array[0,:,:]

    pst = pyemu.Pst(os.path.join("template","freyberg6.pst"))
    par = pst.parameter_data
    par.loc[:,"i"] = -999
    par.loc[:,"j"] = -999
    par.loc[:,"k"] = -999
    par.loc[:,"kper"] = -999

    names = par.loc[par.parnme.apply(lambda x: "npf" in x or "sto" in x),"parnme"]
    par.loc[names,"i"] = names.apply(lambda x: int(x.split('_')[-2]))
    par.loc[names,"j"] = names.apply(lambda x: int(x.split('_')[-1]))
    par.loc[names,"k"] = names.apply(lambda x: int(x.split('_')[-3]))
    
    names = par.loc[par.parnme.apply(lambda x: "wel" in x or "rch" in x),"parnme"]
    par.loc[names,"kper"] = names.apply(lambda x: int(x.split('_')[-1]))

    if pval_series is None:
        print("using parval1")
        pval_series = par.parval1.copy()

    pst.add_transform_columns()
    ub = pst.parameter_data.parubnd_trans
    lb = pst.parameter_data.parlbnd_trans


    fig = plt.figure(figsize=(6.75,8))
    ax_count = 0
    arr_cmap = "plasma"
    ib_cmap = plt.get_cmap("Greys_r")
    ib_cmap.set_bad(alpha=0.0)
    ib = np.ma.masked_where(ib!=0,ib)
    axes = []
    tags = ["npf_k_","npf_k33_","sto_ss_","sto_sy_"]
    prefixes = ["HK","VK","SS","SY"]
    cb_labels = ["HK $log_{10} \\frac{ft}{d}$","VK $log_{10} \\frac{ft}{d}$","SS $log_{10} \\frac{1}{ft}$","SY $log_{10} \\frac{ft}{ft}$"]
    for irow,(tag,prefix,cb_label) in enumerate(zip(tags,prefixes,cb_labels)):
        ppar = par.loc[par.parnme.str.startswith(tag),:]
        # prefix = "HK"
        # cb_label = "HK $\\frac{ft}{d}$"
        #vmin,vmax = pval_series.loc[ppar.parnme].apply(np.log10).min(),pval_series.loc[ppar.parnme].apply(np.log10).max()
        #vmin = truth_pval_series.loc[ppar.parnme].apply(np.log10).min()
        #vmax = truth_pval_series.loc[ppar.parnme].apply(np.log10).max()
        vmin = lb.loc[ppar.parnme].min()
        vmax = ub.loc[ppar.parnme].max()

        for k in range(ppar.k.max()+1):
            print(k)
            ax = plt.subplot2grid((8,3),(irow*2,k),rowspan=2)
            arr = np.zeros_like(ib,dtype=np.float32)
            kppar = ppar.loc[ppar.k==k,:]

            arr[kppar.i,kppar.j] = pval_series.loc[kppar.parnme].values
            arr = np.log10(arr)
            cb = ax.imshow(arr,cmap=arr_cmap,extent=m.modelgrid.extent,vmin=vmin,vmax=vmax,)
            cb = plt.colorbar(cb)
            cb.set_label(cb_label)
            ax.imshow(ib,cmap=ib_cmap,extent=m.modelgrid.extent)
            ax.set_title("{0}) Layer {1} {2}".format(abet[ax_count],k+1,prefix),loc="left")
            ax.set_xlabel("x $ft$")
            ax.set_ylabel("y $ft$")

            axes.append(ax)
            ax_count += 1
    

    #rch
    ax = plt.subplot2grid((8,3),(irow*2,1),colspan=2)
    ppar = par.loc[par.parnme.str.startswith("rch"),:]
    print(pval_series.loc[ppar.parnme].values)
    ppar.sort_values(by="kper",ascending=True,inplace=True)
    ax.bar(ppar.kper+1,pval_series.loc[ppar.parnme].values)
    ax.set_title("{0}) Recharge".format(abet[ax_count]),loc="left")
    ax.set_ylabel("Recharge $\\frac{ft}{d}$")
    ax.set_xlabel("stress period")
    axes.append(ax)
    ax_count += 1
    #wel flux
    ax = plt.subplot2grid((8,3),(irow*2+1,1),colspan=2)
    ppar = par.loc[par.parnme.str.startswith("wel"),:].copy()
    print(ppar.parval1)
    ppar.loc[:,"parval1"] = pval_series.loc[ppar.parnme]
    print(ppar.parval1)
    cumflux = ppar.groupby("kper").sum().loc[:,"parval1"]
    print(cumflux)
    ax.bar(cumflux.index+1,cumflux.values)
    ax.set_title("{0}) Total well extraction".format(abet[ax_count]),loc="left")
    ax.set_ylabel("extraction $\\frac{ft^3}{d}$")
    ax.set_xlabel("stress period")
    axes.append(ax)
    ax_count += 1

    plt.tight_layout()
    plt.savefig(os.path.join(plt_dir,plt_name))
    plt.close(fig)
    pst.write(os.path.join("template",(plt_name.replace(".pdf",".pst"))),version=2)


def plot_domain():
    sim = flopy.mf6.MFSimulation.load(sim_ws="template")
    m = sim.get_model("freyberg6")
    wel_data = m.wel.stress_period_data.array[0]
    sfr_data = m.sfr.packagedata.array
    ghb_data = m.ghb.stress_period_data.array[0]
    ib = m.dis.idomain.array[0, :, :]
    ib_cmap = plt.get_cmap("Greys_r")
    ib_cmap.set_bad(alpha=0.0)

    pst = pyemu.Pst(os.path.join("template","freyberg6_run.pst"))

    fig,ax = plt.subplots(1,1,figsize=(8,8))

    ib = np.ma.masked_where(ib!=0,ib)

    ax.imshow(ib,cmap=ib_cmap,extent=m.modelgrid.extent)
    for ii,cid in enumerate(wel_data.cellid):
        i,j = cid[1],cid[2]
        verts = m.modelgrid.get_cell_vertices(i, j)
        if ii == wel_data.shape[0] - 1:
            p = Polygon(verts, facecolor='b',label="extraction well")
        else:
            p = Polygon(verts, facecolor='b')
        ax.add_patch(p)
    for ii,cid in enumerate(ghb_data.cellid):
        i,j = cid[1],cid[2]
        verts = m.modelgrid.get_cell_vertices(i, j)
        if ii == ghb_data.shape[0] - 1:
            p = Polygon(verts, facecolor='m',label="GHB cell")
        else:
            p = Polygon(verts, facecolor='m')
        ax.add_patch(p)
    for cid in sfr_data.cellid:
        i, j = cid[1], cid[2]
        verts = m.modelgrid.get_cell_vertices(i, j)
        if i < 20:
            c = "g"
        else:
            c = "c"
        if i == 19:
            p = Polygon(verts, facecolor=c, label="SFR headwater reaches")
        elif i == 39:
            p = Polygon(verts, facecolor=c, label="SFR tailwater reaches")
        else:
            p = Polygon(verts, facecolor=c)
        ax.add_patch(p)
    #ax.text(m.modelgrid.xcellcenters[10,j],m.modelgrid.ycellcenters[10,j],"headwater sw-gw exchange forecast reaches",rotation=90,ha="center",va="center")
    #spnspecs.add_text(ax=ax,x=m.modelgrid.xcellcenters[10,j],y=m.modelgrid.ycellcenters[10,j],text="headwater sw-gw exchange forecast reaches",
    #                  rotation=90,ha="center",va="center",bold=False, italic=False,transform=False,bbox={"facecolor":"none","edgecolor":"none"})
    #ax.text(m.modelgrid.xcellcenters[30, j], m.modelgrid.ycellcenters[30, j], "tailwater sw-gw exchange forecast reaches", rotation=90, ha="center",
    #        va="center")
    #spnspecs.add_text(ax=ax, x=m.modelgrid.xcellcenters[30, j], y=m.modelgrid.ycellcenters[30, j],
    #                  text="tailwater sw-gw exchange forecast reaches", rotation=90, ha="center", va="center")

    x = m.modelgrid.xcellcenters[i, j]
    y = m.modelgrid.ycellcenters[i, j]
    ax.scatter([x],[y],marker="^",c='r',s=100,zorder=10,label="point observation/forecast location")
    ax.text(x + 150, y+150, "sw_1".format(1), zorder=11, bbox=dict(facecolor='w', alpha=1,edgecolor="none",pad=1))
    ylim = ax.get_ylim()

    nz_obs = pst.observation_data.loc[pst.nnz_obs_names,:].copy()
    nz_obs = nz_obs.loc[nz_obs.obsnme.str.startswith("trgw"),:]
    nz_obs.loc[:, "i"] = nz_obs.obsnme.apply(lambda x: int(x.split('_')[2]))
    nz_obs.loc[:, "j"] = nz_obs.obsnme.apply(lambda x: int(x.split('_')[3]))
    for ii,g in enumerate(nz_obs.obgnme.unique()):
        i = nz_obs.loc[nz_obs.obgnme == g, "i"][0]
        j = nz_obs.loc[nz_obs.obgnme == g, "j"][0]
        x = m.modelgrid.xcellcenters[i, j]
        y = m.modelgrid.ycellcenters[i, j]
        ax.scatter([x], [y], marker="^", c='r', s=100, zorder=10)
        ax.text(x + 150, y+150,"gw_{0}".format(ii+1),zorder=11, bbox=dict(facecolor='w', alpha=1,edgecolor="none",pad=1))

    gw_fore = [f for f in forecasts if "trgw" in f][0]
    i = int(gw_fore.split('_')[2])
    j = int(gw_fore.split('_')[3])
    x = m.modelgrid.xcellcenters[i, j]
    y = m.modelgrid.ycellcenters[i, j]
    ax.scatter([x], [y], marker="^", c='r', s=100, zorder=10)
    ax.text(x + 150, y + 150, "gw_3", zorder=11,
            bbox=dict(facecolor='w', alpha=1, edgecolor="none", pad=1))

    top = m.dis.top.array.reshape(ib.shape)
    top[top<0] = np.NaN
    cb = ax.imshow(top,extent=m.modelgrid.extent,cmap="bone")
    cb = plt.colorbar(cb,pad=0.01)
    cb.set_label("top $ft$")
    ax.set_xlabel("x $ft$")
    ax.set_ylabel("y $ft$")

    ax.set_ylim(0,ylim[1])
    spnspecs.graph_legend(ax=ax,bbox_to_anchor=(1.15, 0.9))
    plt.tight_layout()
    plt.savefig(os.path.join(plt_dir,"domain.pdf"))
    plt.close("all")


def make_opt_figs2():

    n_m_d = "master_opt_neutral"
    a_m_d = "master_opt_stack"
    f_m_d = "master_opt_fosm"
    assert os.path.exists(n_m_d)
    assert os.path.exists(a_m_d)
    sim = flopy.mf6.MFSimulation.load(sim_ws=a_m_d)
    m = sim.get_model("freyberg6")
    perlen = sim.tdis.perioddata.array["perlen"]
    totim = np.cumsum(perlen).astype(int)
    sp_start = pd.to_datetime("12-31-2015") + pd.to_timedelta(totim, unit='d')

    pst_file = "freyberg6_run_opt.pst"
    pst = pyemu.Pst(os.path.join(a_m_d, pst_file))
    w_par = pst.parameter_data.loc[pst.parameter_data.pargp == "welflux", :]
    w_par.loc[:, "well"] = w_par.parnme.apply(lambda x: "_".join(x.split("_")[:-1]))
    w_par.loc[:, "kper"] = w_par.parnme.apply(lambda x: int(x.split('_')[-1]))
    w_par.loc[:, "i"] = w_par.parnme.apply(lambda x: int(x.split('_')[-3]))
    w_par.loc[:, "j"] = w_par.parnme.apply(lambda x: int(x.split('_')[-2]))
    w_par.loc[:, "datetime"] = w_par.kper.apply(lambda x: sp_start[x])
    a_par = pyemu.pst_utils.read_parfile(os.path.join(a_m_d, pst_file.replace(".pst", ".1.par")))
    n_par = pyemu.pst_utils.read_parfile(os.path.join(n_m_d, pst_file.replace(".pst", ".1.par")))
    f_par = pyemu.pst_utils.read_parfile(os.path.join(f_m_d, pst_file.replace(".pst", ".1.par")))

    # con = pst.observation_data.loc[pst.nnz_obs_names,:]
    obs = pst.observation_data
    con = obs.loc[obs.obgnme.apply(lambda x: x.startswith("less_"))]
    con.loc[:, "obgnme"] = con.obsnme.apply(lambda x: x.split('_')[0])
    con_groups = con.obgnme.unique()
    con_groups.sort()
    con.loc[:, "datetime"] = pd.to_datetime(con.obsnme.apply(lambda x: x.split('_')[-1]))

    a_res = pyemu.pst_utils.read_resfile(os.path.join(a_m_d, pst_file.replace(".pst", ".1.est.rei")))
    n_res = pyemu.pst_utils.read_resfile(os.path.join(n_m_d, pst_file.replace(".pst", ".1.est.rei")))
    f_res = pyemu.pst_utils.read_resfile(os.path.join(f_m_d, pst_file.replace(".pst", ".1.est.rei")))

    wel_names = w_par.well.unique()
    wel_names.sort()
    # fig,axes = plt.subplots(2,1,figsize=(8,4))
    fig = plt.figure(figsize=(8, 5))
    axes = [plt.subplot2grid((3, 3), (i, 0), colspan=2) for i in range(3)]
    x = np.arange(1,w_par.kper.max()+1)
    cmap = plt.get_cmap('terrain')

    axes_t = []
    cval = float(pst.prior_information.iloc[0].equation.split('=')[-1])

    for ax, pvals, label, res in zip(axes, [n_par, f_par, a_par], ["A) risk neutral", "B) FOSM-based risk averse","C) stack-based risk-averse"], [n_res, f_res, a_res]):
        offset = -0.2
        step = 0.1

        for i, w in enumerate(wel_names):
            ww_par = w_par.loc[w_par.well == w, :].copy()
            ww_par.sort_values(by="kper", inplace=True)
            print(x.shape, pvals.loc[ww_par.parnme, "parval1"].shape)
            ax.bar(x + offset, pvals.loc[ww_par.parnme, "parval1"], width=step, facecolor=cmap(i / len(wel_names)),
                   edgecolor="none")
            offset += step


        ax.set_xticks(x)
        ax.set_xticklabels(np.arange(1, 26))
        # ax.set_xticklabels(sp_start,rotation=90)
        ax.set_xlabel("stress period")
        #ax.set_xticks([])
        ax.set_ylabel("well extraction rate $\\frac{ft^3}{d}$")
        ax.set_title(label, loc="left")
        ylim = ax.get_ylim()
        ax.set_ylim(ylim[0], ylim[1])
        axt = plt.twinx(ax)
        tot = []
        for kper in range(1,w_par.kper.max()+1):
            kw_par = w_par.loc[w_par.kper==kper,"parnme"]

            tot.append(pvals.loc[kw_par,"parval1"].sum())
        axt.plot(x,tot,"k",label="total extraction rate")
        xlim = axt.get_xlim()
        axt.plot(xlim, [cval, cval], "r--", label="required minimum total extraction rate")
        #ylim = axt.get_ylim()
        axt.set_ylim(0,max(tot)*1.2)
        ax.set_xlim(xlim)
        axes_t.append(axt)

    spnspecs.graph_legend(ax=axes_t[0], loc="upper left", facecolor='w', framealpha=1.0, bbox_to_anchor=(1.2, 1.0))
    # axes[0].set_xticklabels([])
    # axes[1].set_xticklabels(sp_start.map(lambda x: x.strftime("%d-%m-%Y")))
    #mx = max([ax.get_ylim()[1] for ax in axes_t])
    #mn = max([ax.get_ylim()[0] for ax in axes_t])
    #for ax in axes_t:
    #    ax.set_ylim(mn * 1.5, mx)
    #    ax.set_ylabel("simulated sw-gw flux $\\frac{ft^3}{d}$")
    for ax in axes_t:
        ax.set_ylabel("total extraction rate $\\frac{ft^3}{d}$")
    # ax = plt.subplot2grid((2,3),(0,2),rowspan=2)
    ax = plt.axes((.75, 0.1, 0.2, 0.6))
    ib = m.dis.idomain.array[0, :, :]
    icmap = plt.get_cmap("Greys_r")
    icmap.set_bad(alpha=0.0)
    ib = np.ma.masked_where(ib != 0, ib)
    ax.imshow(ib, cmap=icmap, extent=m.modelgrid.extent)
    for ii, w in enumerate(wel_names):
        ww_par = w_par.loc[w_par.well == w, :]
        i = ww_par.i[0]
        j = ww_par.j[0]
        verts = m.modelgrid.get_cell_vertices(i, j)
        p = Polygon(verts, facecolor=cmap(ii / len(wel_names)))
        ax.add_patch(p)
    # ax.set_xticks([])
    # ax.set_yticks([])
    ax.set_ylabel("x $ft$")
    ax.set_xlabel("y $ft$")
    ax.set_title("D) location of extraction wells", loc="left")

    plt.tight_layout()
    plt.savefig(os.path.join(plt_dir, "opt1.pdf"))

    fig = plt.figure(figsize=(8, 6))
    axes = [plt.subplot2grid((3, 1), (i, 0)) for i in range(3)]
    x = np.arange(w_par.kper.max() + 1)
    cmap = plt.get_cmap('plasma')

    for axt, pvals, label, res in zip(axes, [n_par, f_par, a_par],
                                     ["A) risk neutral", "B) FOSM-based risk averse", "C) stack-based risk-averse"],
                                     [n_res, f_res, a_res]):

        for g, c in zip(con_groups, ["g", "c"]):
            g_con = con.loc[con.obgnme == g]
            g_con.sort_values(by="datetime", inplace=True)
            cval = -1* g_con.obsval[0]
            axt.plot(-1*res.loc[g_con.obsnme, "modelled"].values, c, label=g + " sw-gw flux")
            xlim = axt.get_xlim()


        axt.plot(xlim, [cval, cval], "r--", label="required minimum sw-gw flux")
        axt.set_xlim(xlim)
        axt.set_xticks(x)
        axt.set_xlabel("stress period")
        axt.set_title(label,loc="left")

    spnspecs.graph_legend(ax=axes[0], loc="upper right", facecolor='w', framealpha=1.0)#, bbox_to_anchor=(1.2, 1.0))
    #axes[0].set_xticklabels([])

    for ax in axes:
       ax.set_ylabel("simulated sw-gw flux $\\frac{ft^3}{d}$")

    plt.tight_layout()
    plt.savefig(os.path.join(plt_dir, "opt2.pdf"))

if __name__ == "__main__":

    prep_mf6_model()
    # setup_pest_interface()
    # build_and_draw_prior()
    # run_prior_sweep()
    # set_truth_obs()
    # #
    # run_ies_demo()
    # make_ies_figs()
    # make_ies_figs(m_d="master_ies_default",plt_case="ies_default")
    # # #
    # # #
    # run_glm_demo()
    #make_glm_figs()
    # make_glm_figs(m_d="master_glm_default",plt_case="glm_default")
    # # # # # #
    # run_sen_demo()
    # make_sen_figs()
    #
    # run_opt_demo()
    # make_opt_figs2()

    # plot_domain()

    #invest()

    # plot_par_vector()
    #start()
    #invest()
    # start()
    # _rebase_results()
    # compare_to_baseline()
    # test(True)
