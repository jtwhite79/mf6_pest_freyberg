import os
import shutil
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages as pdf
import matplotlib.pyplot as plt
import flopy
import pyemu


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
        f.write("BEGIN CONTINUOUS FILEOUT sfr.csv\nheadwater outflow headwater\n")
        f.write("tailwater outflow tailwater\ngage_1 inflow 40\nEND CONTINUOUS")

    shutil.copy2(os.path.join(org_ws,"mf6.exe"),os.path.join(new_ws,"mf6.exe"))
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
    pst.model_command = "mf6"
    pst.control_data.noptmax = 0
    pst.pestpp_options["additional_ins_delimiters"] = ","
    pst.pestpp_options["ies_num_reals"] = 50
    pst.write(os.path.join(ws,"freyberg6.pst"))
    pyemu.os_utils.run("pestpp-ies freyberg6.pst",cwd=ws)

    #pst.control_data.noptmax = -1
    #pst.write(os.path.join(ws, "freyberg6.pst"))
    #pyemu.os_utils.start_workers(ws,"pestpp-ies","freyberg6.pst",num_workers=10,master_dir="master_ies")


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
    temporal_v = pyemu.geostats.ExpVario(contribution=1.0,a=60)
    spatial_gs = pyemu.geostats.GeoStruct(variograms=spatial_v)
    temporal_gs = pyemu.geostats.GeoStruct(variograms=temporal_v)

    static_struct_dict = {spatial_gs:[]}
    for pargp in static_par.pargp.unique():
        static_struct_dict[spatial_gs].append(static_par.loc[static_par.pargp==pargp,["parnme","x","y","i","j"]])
    temporal_struct_dict = {temporal_gs: [rch_par.loc[:, ["parnme", "x", "y"]]]}
    for pargp in wel_par.pargp.unique():
        temporal_struct_dict[temporal_gs].append(wel_par.loc[wel_par.pargp == pargp, ["parnme", "x", "y"]])
    ss = pyemu.geostats.SpecSim2d(m.modelgrid.delr,m.modelgrid.delc,spatial_gs)
    #pe = ss.grid_par_ensemble_helper(pst,static_par,num_reals=300,sigma_range=4.0)
    #temporal_pe = pyemu.helpers.geostatistical_draws(pst,struct_dict=temporal_struct_dict,num_reals=300)
    struct_dict = static_struct_dict
    for k,v in temporal_struct_dict.items():
        struct_dict[k] = v
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
            f.write("l1 w")
            for c in sfr_obs_df.columns:
                name = "{0}_{1}".format(c.lower(), totim_dict[int(i)])
                f.write(" !{0}! ".format(name))
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
            f.write("l1 w")
            for c in obs_df.columns:
                name = "{0}_{1}".format(c.lower(), totim_dict[int(i)])
                f.write(" !{0}! ".format(name))
                names.append(name)
                vals.append(obs_df.loc[i,c])
            f.write("\n")
    head_df = pd.DataFrame({"obsnme":names,"obsval":vals},index=names)
    head_df.loc[:,"obgnme"] = head_df.obsnme.apply(lambda x: '_'.join(x.split('_')[:3]))

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
    df.loc[:,"parubnd"] = df.parval1 + (df.parval1.apply(np.abs) * 0.5)
    df.loc[:, "parlbnd"] = df.parval1 - (df.parval1.apply(np.abs) * 0.5)
    dfs = [df]

    # write a recharge template
    f_in = open(os.path.join(ws,"freyberg6.rch"),'r')
    names,vals = [],[]
    with open(os.path.join(ws,"freybeg6.rch.tpl"),'w') as f:
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
    df.loc[:,"parubnd"] = df.parval1 * 1.2
    df.loc[:,"parlbnd"] = df.parval1 * 0.8
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
    pst.write(os.path.join(t_d,pst_file))
    m_d = "master_prior"
    pyemu.os_utils.start_workers(t_d,"pestpp-ies",pst_file,num_workers=15,master_dir=m_d)

def set_truth_obs():
    t_d = "template"
    m_d = "master_prior"
    assert os.path.exists(m_d)
    pst = pyemu.Pst(os.path.join(m_d,"freyberg6_sweep.pst"))
    oe = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(m_d,"freyberg6_sweep.0.obs.csv"))
    pv = oe.phi_vector
    pv.sort_values(inplace=True)
    #idx = pv.index[int(pv.shape[0]/2)]
    idx = pv.index[-1]
    pst.observation_data.loc[:,"obsval"] = oe.loc[idx,pst.obs_names]
    pst.observation_data.loc[:,"weight"] = 0.0
    obs = pst.observation_data
    obs.loc[obs.obsnme.apply(lambda x: "2016" in x and ("trgw_0_29_15" in x or "trgw_0_2_9" in x)),"weight"] = 5.0
    obs.loc[obs.obsnme.apply(lambda x: "gage_1" in x and "2016" in x), "weight"] = 0.005
    pst.control_data.noptmax = 0
    pst.pestpp_options["forecasts"] = ["headwater_20171231","tailwater_20171231"]
    pst.write(os.path.join(t_d,"freyberg6_run.pst"))

    pyemu.os_utils.run("pestpp-ies.exe freyberg6_run.pst",cwd=t_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run.pst"))
    print(pst.phi_components)

def run_ies_demo():
    t_d = "template"
    assert os.path.exists(t_d)
    pst_file = "freyberg6_run.pst"
    assert os.path.exists(os.path.join(t_d,pst_file))
    pst = pyemu.Pst(os.path.join(t_d,pst_file))
    pst.control_data.noptmax = 3
    pst.pestpp_options = {}
    pst.pestpp_options["ies_par_en"] = "prior.jcb"
    pst.pestpp_options["ies_num_reals"] = 50
    pst.pestpp_options["ies_bad_phi_sigma"] = 1.5
    pst.pestpp_options["additional_ins_delimiters"] = ","
    pst.write(os.path.join(t_d,"freyberg6_run_ies.pst"))
    m_d = "master_ies"
    pyemu.os_utils.start_workers(t_d, "pestpp-ies", "freyberg6_run_ies.pst", num_workers=15, master_dir=m_d)


def make_ies_figs():
    m_d = "master_ies"
    assert os.path.join(m_d)
    pst_file = "freyberg6_run_ies.pst"
    pst = pyemu.Pst(os.path.join(m_d,pst_file))
    pr_oe = pd.read_csv(os.path.join(m_d,pst_file.replace(".pst",".0.obs.csv")))
    pt_oe = pd.read_csv(os.path.join(m_d, pst_file.replace(".pst", ".{0}.obs.csv".format(pst.control_data.noptmax))))
    obs = pst.observation_data
    for nz_grp in pst.nnz_obs_groups:
        grp_obs = obs.loc[obs.obgnme==nz_grp,:].copy()
        grp_obs.loc[:,"dateime"] = pd.to_datetime(grp_obs.obsnme.apply(lambda x: x.split('_')[-1]))


def invest():
    pst = pyemu.Pst(os.path.join("master_ies","freyberg6_run_ies.pst"))
    pyemu.helpers.setup_fake_forward_run(pst,"fake.pst","master_ies",new_cwd="master_ies")

if __name__ == "__main__":
    #prep_mf6_model()
    #setup_pest_interface()
    #build_and_draw_prior()
    #run_prior_sweep()
    #set_truth_obs()
    #run_ies_demo()
    #make_ies_figs()
    #invest()