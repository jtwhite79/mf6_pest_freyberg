import os
import shutil
import string

import numpy as np
import pandas as pd

from matplotlib.backends.backend_pdf import PdfPages as pdf
import matplotlib.pyplot as plt
import flopy
import pyemu

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
        f.write("BEGIN CONTINUOUS FILEOUT sfr.csv\nheadwater sfr headwater\n")
        f.write("tailwater sfr tailwater\ngage_1 inflow 40\nEND CONTINUOUS")

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
    pst.write(os.path.join(ws,"freyberg6.pst"),version=2)
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
    #ss = pyemu.geostats.SpecSim2d(m.modelgrid.delr,m.modelgrid.delc,spatial_gs)
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
    df.loc[:,"parubnd"] = df.parval1 + (df.parval1.apply(np.abs) * 0.5)
    df.loc[:, "parlbnd"] = df.parval1 - (df.parval1.apply(np.abs) * 0.5)
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
    pst.write(os.path.join(t_d,pst_file),version=2)
    m_d = "master_prior"
    pyemu.os_utils.start_workers(t_d,"pestpp-ies",pst_file,num_workers=15,master_dir=m_d)

def set_truth_obs():
    t_d = "template"
    m_d = "master_prior"
    assert os.path.exists(m_d)
    pst = pyemu.Pst(os.path.join(m_d,"freyberg6_sweep.pst"))
    pst.pestpp_options["forecasts"] = ["headwater_20171231","tailwater_20171231","trgw_0_9_1_20171231"]
    oe = pyemu.ObservationEnsemble.from_csv(pst=pst,filename=os.path.join(m_d,"freyberg6_sweep.0.obs.csv"))
    pv = oe.phi_vector
    #pv.sort_values(inplace=True)
    #idx = pv.index[int(pv.shape[0]/2)]
    #idx = pv.index[int(pv.shape[0]/2)]
    oe.sort_values(by=pst.forecast_names[0],inplace=True)
    idx = oe.index[-int(oe.shape[0]/10)]
    pst.observation_data.loc[:,"obsval"] = oe.loc[idx,pst.obs_names]
    pst.observation_data.loc[:,"weight"] = 0.0
    obs = pst.observation_data
    obs.loc[obs.obsnme.apply(lambda x: "2016" in x and ("trgw_0_29_15" in x or "trgw_0_2_9" in x)),"weight"] = 5.0
    obs.loc[obs.obsnme.apply(lambda x: "gage_1" in x and "2016" in x), "weight"] = 0.005
    pst.control_data.noptmax = 0
    
    pst.write(os.path.join(t_d,"freyberg6_run.pst"),version=2)

    pyemu.os_utils.run("pestpp-ies.exe freyberg6_run.pst",cwd=t_d)
    pst = pyemu.Pst(os.path.join(t_d,"freyberg6_run.pst"))
    print(pst.phi_components)

def run_ies_demo():
    '''todo: localize
    '''
    t_d = "template"
    assert os.path.exists(t_d)
    pst_file = "freyberg6_run.pst"
    assert os.path.exists(os.path.join(t_d,pst_file))
    pst = pyemu.Pst(os.path.join(t_d,pst_file))
    pst.control_data.noptmax = 3
    pst.pestpp_options = {"forecasts":pst.pestpp_options["forecasts"]}
    pst.pestpp_options["ies_par_en"] = "prior.jcb"
    pst.pestpp_options["ies_num_reals"] = 50
    pst.pestpp_options["ies_bad_phi_sigma"] = 1.5
    pst.pestpp_options["additional_ins_delimiters"] = ","
    pst.write(os.path.join(t_d,"freyberg6_run_ies.pst"),version=2)
    m_d = "master_ies"
    pyemu.os_utils.start_workers(t_d, "pestpp-ies", "freyberg6_run_ies.pst", num_workers=15, master_dir=m_d)


def make_ies_figs():
    m_d = "master_ies"
    assert os.path.join(m_d)
    pst_file = "freyberg6_run_ies.pst"
    pst = pyemu.Pst(os.path.join(m_d,pst_file))
    pr_oe = pd.read_csv(os.path.join(m_d,pst_file.replace(".pst",".0.obs.csv")))
    pt_oe = pd.read_csv(os.path.join(m_d, pst_file.replace(".pst", ".{0}.obs.csv".format(pst.control_data.noptmax))))
    pr_oe = pyemu.ObservationEnsemble(pst=pst,df=pr_oe)
    pt_oe = pyemu.ObservationEnsemble(pst=pst,df=pt_oe)
    pr_pv = pr_oe.phi_vector
    pt_pv = pt_oe.phi_vector


    obs = pst.observation_data
    print(pst.nnz_obs_groups)
    print(pst.forecast_names)
    fig = plt.figure(figsize=(8,6))
    ax_count = 0
    unit_dict = {"head":"sw-gw flux $\\frac{m}{d}$",
                "tail": "sw-gw flux $\\frac{m}{d}$",
                "trgw" : "gw level $m$",
                "gage" : "sw flux $\\frac{m}{d}$"}
    for i,nz_grp in enumerate(pst.nnz_obs_groups):
        grp_obs = obs.loc[obs.obgnme==nz_grp,:].copy()
        print(grp_obs)
        grp_obs.loc[:,"datetime"] = pd.to_datetime(grp_obs.obsnme.apply(lambda x: x.split('_')[-1]))
        ax = plt.subplot2grid((5,4),(i,0),colspan=4)
        ax.plot(grp_obs.datetime,grp_obs.obsval, 'r')
        [ax.plot(grp_obs.datetime,pr_oe.loc[i,grp_obs.obsnme],'0.5',lw=0.1, alpha=0.25) for i in pr_oe.index]
        [ax.plot(grp_obs.datetime,pt_oe.loc[i,grp_obs.obsnme],'b',lw=0.1,alpha=0.35) for i in pt_oe.index]
        ax.plot(grp_obs.datetime,grp_obs.obsval, 'r')
        ax.set_title("{0}) {1}".format(abet[ax_count],nz_grp),loc="left")
        unit = None
        for tag,u in unit_dict.items():
            if tag in nz_grp:
                unit = u
        ax.set_ylabel(unit)

        ax_count += 1

    ax = plt.subplot2grid((5,4),(3,0),rowspan=2)   
    ax.hist(pr_pv,alpha=0.5,facecolor="0.5",edgecolor="none")
    ax.hist(pt_pv,alpha=0.5,facecolor="b",edgecolor="none")
    ax.set_title("{0}) ensemble $\phi$ distributions".format(abet[ax_count]), loc="left")
    #ax.set_yticks([])
    ax.set_xlabel("$\phi$")
    ax.set_ylabel("number of realizations")
    ax_count += 1

    for i,forecast in enumerate(pst.forecast_names):
        ax = plt.subplot2grid((5,4),(3,i+1),rowspan=2)
        pr_oe.loc[:,forecast].hist(ax=ax,density=True,facecolor='0.5',edgecolor="none",alpha=0.5)
        pt_oe.loc[:,forecast].hist(ax=ax,density=True,facecolor='b',edgecolor="none",alpha=0.5)
        ylim = ax.get_ylim()
        oval = obs.loc[forecast,"obsval"]
        ax.plot([oval,oval],ylim,'r')
        ax.set_title("{0}) {1}".format(abet[ax_count],forecast),loc="left")
        unit = None
        for tag,u in unit_dict.items():
            if tag in forecast:
                unit = u
        ax.set_xlabel(unit)
        ax.set_ylabel("increasing probability density")
        ax.set_yticks([])
        ax.grid(False)
        ax_count +=1 


    plt.tight_layout()
    plt.savefig(os.path.join(m_d,"ies_sum.pdf"))


def make_glm_figs():
    m_d = "master_glm"
    assert os.path.join(m_d)
    pst_file = "freyberg6_run_glm.pst"
    pst = pyemu.Pst(os.path.join(m_d,pst_file))
    pt_oe = pd.read_csv(os.path.join(m_d,pst_file.replace(".pst",".post.obsen.csv")))
    pt_oe = pyemu.ObservationEnsemble(pst=pst,df=pt_oe)
    f_df = pd.read_csv(os.path.join(m_d,pst_file.replace(".pst",".pred.usum.csv")),index_col=0)
    f_df.index = f_df.index.map(str.lower)
    pv = pt_oe.phi_vector
    obs = pst.observation_data
    print(pst.nnz_obs_groups)
    print(pst.forecast_names)
    fig = plt.figure(figsize=(8,6))
    ax_count = 0
    unit_dict = {"head":"sw-gw flux $\\frac{m}{d}$",
                "tail": "sw-gw flux $\\frac{m}{d}$",
                "trgw" : "gw level $m$",
                "gage" : "sw flux $\\frac{m}{d}$"}
    for i,nz_grp in enumerate(pst.nnz_obs_groups):
        grp_obs = obs.loc[obs.obgnme==nz_grp,:].copy()
        print(grp_obs)
        grp_obs.loc[:,"datetime"] = pd.to_datetime(grp_obs.obsnme.apply(lambda x: x.split('_')[-1]))
        ax = plt.subplot2grid((5,4),(i,0),colspan=4)
        ax.plot(grp_obs.datetime,grp_obs.obsval, 'r')
        [ax.plot(grp_obs.datetime,pt_oe.loc[i,grp_obs.obsnme],'b',lw=0.1,alpha=0.5) for i in pt_oe.index]
        ax.plot(grp_obs.datetime,grp_obs.obsval, 'r')
        ax.set_title("{0}) {1}".format(abet[ax_count],nz_grp),loc="left")
        unit = None
        for tag,u in unit_dict.items():
            if tag in nz_grp:
                unit = u
        ax.set_ylabel(unit)
        ax.grid(False)

        ax_count += 1

    ax = plt.subplot2grid((5,4),(3,0),rowspan=2)   
    ax.hist(pv,alpha=0.5,facecolor="b",edgecolor="none")
    ax.set_title("{0}) posterior ensemble $\phi$ distribution".format(abet[ax_count]), loc="left")
    #ax.set_yticks([])
    ax.set_xlabel("$\phi$")
    ax.set_ylabel("number of realizations")
    ax_count += 1
    
    for i,forecast in enumerate(pst.forecast_names):
        ax = plt.subplot2grid((5,4),(3,i+1),rowspan=2)
        axt = plt.twinx()
        x,y = pyemu.plot_utils.gaussian_distribution(f_df.loc[forecast,"prior_mean"],f_df.loc[forecast,"prior_stdev"])
        axt.fill_between(x,0,y,facecolor="0.5",alpha=0.25)
        x,y = pyemu.plot_utils.gaussian_distribution(f_df.loc[forecast,"post_mean"],f_df.loc[forecast,"post_stdev"])
        axt.fill_between(x,0,y,facecolor="b",alpha=0.25)

        pt_oe.loc[:,forecast].hist(ax=ax,density=True,facecolor='b',edgecolor="none",alpha=0.5)
        ylim = ax.get_ylim()
        oval = obs.loc[forecast,"obsval"]
        ax.plot([oval,oval],ylim,'r')
        ax.set_title("{0}) {1}".format(abet[ax_count],forecast),loc="left")
        unit = None
        for tag,u in unit_dict.items():
            if tag in forecast:
                unit = u
        ax.set_xlabel(unit)
        ax.set_ylabel("increasing probability density")
        ax.set_yticks([])
        ax.grid(False)
        ylim = axt.get_ylim()
        axt.set_ylim(0,ylim[1])
        

        ax_count += 1


    plt.tight_layout()
    plt.savefig(os.path.join(m_d,"glm_sum.pdf"))

def invest():
    pst = pyemu.Pst(os.path.join("master_ies","freyberg6_run.pst"))
    par = pst.parameter_data
    r_pars = par.loc[par.parnme.apply(lambda x: 'rch' in x),"parnme"]
    print(r_pars)

    pst.parameter_data.loc[r_pars,"parval1"] = par.loc[r_pars,"parubnd"]
    pst.control_data.noptmax = 0
    pst.write(os.path.join("template","test1.pst"))
    pyemu.os_utils.run("pestpp-ies test1.pst",cwd="template")
    shutil.copy2(os.path.join("template","freyberg6.lst"),os.path.join("template","test1.lst"))
    shutil.copy2(os.path.join("template","freyberg6.rch"),os.path.join("template","test1.rch"))
    pst.parameter_data.loc[r_pars,"parval1"] = par.loc[r_pars,"parlbnd"]
    pst.control_data.noptmax = 0
    pst.write(os.path.join("template","test2.pst"))
    pyemu.os_utils.run("pestpp-ies test2.pst",cwd="template")
    shutil.copy2(os.path.join("template","freyberg6.lst"),os.path.join("template","test2.lst"))
    shutil.copy2(os.path.join("template","freyberg6.rch"),os.path.join("template","test2.rch"))
    pst1 = pyemu.Pst(os.path.join("template","test1.pst"))
    pst2 = pyemu.Pst(os.path.join("template","test2.pst"))
    d = pst1.res.modelled - pst2.res.modelled
    print(d)


    

def run_glm_demo():

    t_d = "template"
    assert os.path.exists(t_d)
    pst_file = "freyberg6_run.pst"
    assert os.path.exists(os.path.join(t_d,pst_file))
    pst = pyemu.Pst(os.path.join(t_d,pst_file))
    # fix/tie to reduce adj par numbers
    tie_step = 6
    par = pst.parameter_data
    par.loc[:,"partied"] = np.NaN
    gr_par = par.loc[par.pargp.apply(lambda x: "npf" in x or "sto" in x),:].copy()
    gr_par.loc[:,"k"] = gr_par.parnme.apply(lambda x: int(x.split('_')[2]))
    gr_par.loc[:,"i"] = gr_par.parnme.apply(lambda x: int(x.split('_')[3]))
    gr_par.loc[:,"j"] = gr_par.parnme.apply(lambda x: int(x.split('_')[4]))
    
    for i in range(tie_step,gr_par.i.max()+tie_step,tie_step):
        for j in range(tie_step,gr_par.j.max()+tie_step,tie_step):
            ir = [ii for ii in range(i-tie_step,i)]
            jr = [jj for jj in range(j-tie_step,j)]
            t_par = gr_par.loc[gr_par.apply(lambda x:x.i in ir and x.j in jr, axis=1),:]
            for grp in t_par.pargp.unique():
                tt_par = t_par.loc[t_par.pargp==grp,:]
                par.loc[tt_par.parnme[1:],"partrans"] = "tied"
                par.loc[tt_par.parnme[1:],"partied"] = tt_par.parnme[0]
    w_par = par.loc[par.parnme.str.startswith("wel"),:]
    for grp in w_par.pargp.unique():
        ww_par = w_par.loc[w_par.pargp==grp,:]
        par.loc[ww_par.parnme[1:],"partrans"] = "tied"
        par.loc[ww_par.parnme[1:],"partied"] = ww_par.parnme[0]


    pst.control_data.noptmax = 3
    pst.pestpp_options = {"forecasts":pst.pestpp_options["forecasts"]}
    pst.pestpp_options["additional_ins_delimiters"] = ","
    pst.pestpp_options["n_iter_super"] = 999
    pst.pestpp_options["n_iter_base"] = -1
    pst.pestpp_options["glm_num_reals"] = 50
    pst.pestpp_options["glm_accept_mc_phi"] = True
    pst.write(os.path.join(t_d,"freyberg6_run_glm.pst"),version=2)
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

    pst.pestpp_options = {"forecasts":pst.pestpp_options["forecasts"]}
    pst.pestpp_options["additional_ins_delimiters"] = ","
    pst.pestpp_options["tie_by_group"] = True
    pst.write(os.path.join(t_d,"freyberg6_run_sen.pst"),version=2)
    m_d = "master_sen_morris"
    pyemu.os_utils.start_workers(t_d, "pestpp-sen", "freyberg6_run_sen.pst", num_workers=15, master_dir=m_d)
    pst.pestpp_options['gsa_method'] = "sobol"
    pst.write(os.path.join(t_d,"freyberg6_run_sen.pst"),version=2)
    m_d = "master_sen_sobol"
    pyemu.os_utils.start_workers(t_d, "pestpp-sen", "freyberg6_run_sen.pst", num_workers=15, master_dir=m_d)

def start():
    pyemu.os_utils.start_workers("template", "pestpp-sen", "freyberg6_run_sen.pst", num_workers=15)

def make_sen_figs():
    m_d = "master_sen_morris"
    assert os.path.exists(m_d)
    pst_file = "freyberg6_run_sen.pst"
    pst = pyemu.Pst(os.path.join(m_d,pst_file))
    mio_df = pd.read_csv(os.path.join(m_d,pst_file.replace(".pst",".mio")),index_col=0)
    phi_df = pd.read_csv(os.path.join(m_d,pst_file.replace(".pst",".msn")),skipinitialspace=True).loc[:,["parameter_name","sen_mean_abs", "sen_std_dev"]]
    
    df = mio_df.loc[mio_df.index.map(lambda x: x in pst.forecast_names),["parameter_name","sen_mean_abs", "sen_std_dev"]].copy()
    phi_df.index = ["phi" for _ in range(phi_df.shape[0])]
    df = pd.concat([df,phi_df])
    #print(df)
    #print(phi_df.columns)
    #print(df.columns)
    fig,ax = plt.subplots(1,1,figsize=(8,4))
    x = np.arange(df.parameter_name.unique().shape[0])
    offset = -0.2
    step = 0.1
    name_dict = {"headwater":"headwater forecast","tailwater":"tailwater forecast","trgw":"groundwater level forecast","phi":"phi"}
    par_dict = {"k33_0":"layer 1 VK","k33_1":"layer 2 VK","k33_2":"layer 3 VK",
    "ss_0":"layer 1 SS","ss_1":"layer 2 SS","ss_2":"layer 3 SS",
    "sy":"layer 1 SY","rch":"recharge","wel":"well extraction rates","k_0":"layer 1 HK","k_1":"layer 2 HK","k_2":"layer 3 HK"}
    df.index = df.index.map(lambda x: name_dict[x.split('_')[0]])
    df.loc[:,"parameter_name"] = df.parameter_name.apply(lambda x: [v for k,v in par_dict.items() if k in x][0])
    for o in df.index.unique():
        print(o)
        odf = df.loc[o,:].copy()
        odf.loc[:,"sen_mean_abs"] /= odf.sen_mean_abs.sum()
        odf.loc[:,"sen_mean_abs"] *= 100.0
        odf.sort_values(by="parameter_name",inplace=True)

        ax.bar(x+offset,odf.sen_mean_abs,width=step,label=o)
        offset += step
    ax.set_xticks(x)
    ax.set_xticklabels(odf.parameter_name.values,rotation=90)
    ax.legend()
    ax.set_ylabel("percent of total mean absolute sensitivity")
    plt.tight_layout()
    plt.savefig(os.path.join(m_d,"morris.pdf"))

if __name__ == "__main__":
    # prep_mf6_model()
    # setup_pest_interface()
    # build_and_draw_prior()
    # run_prior_sweep()
    #set_truth_obs()
    
    #run_ies_demo()
    #run_glm_demo()
    #run_sen_demo()
    #make_ies_figs()
    #make_glm_figs()
    make_sen_figs()
    #invest()
    #start()