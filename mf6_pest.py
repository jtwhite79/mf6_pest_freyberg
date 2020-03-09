import os
import shutil
import numpy as np
import pandas as pd
import flopy
import pyemu
# add "REPLACE" to the sfr out in temp_history/freyberg.nam
# then run mf5to6 in temp_history to make freyberg6.nam
#

def prep_mf6_model():
    org_ws = "temp_history"
    new_ws = "test"
    if os.path.exists(new_ws):
        shutil.rmtree(new_ws)
    os.mkdir(new_ws)
    sim = flopy.mf6.MFSimulation.load(sim_ws=org_ws)
    sim.simulation_data.mfpath.set_sim_path("test")
    m = sim.get_model("freyberg6")
    obs_df = pd.read_csv("obs_loc.csv")
    obs_df.loc[:,"name"] = obs_df.apply(lambda x: "trgw_{0}_{1}".format(x.row,x.col),axis=1)
    obs_df.loc[:,"layer"] = 3
    obs_df.loc[:,"obstype"] = "HEAD"
    obs_df.loc[:,"col"] = obs_df.col.apply(np.int)
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


    # for pack,attr in props_trans:
    #     #print(m.get_package(pack).__getattribute__(attr))
    #     #for kper in range(m.nper.data):
    #     filename = "{0}_{1}.dat".format(pack,attr)
    #     m.get_package(pack).__getattribute__(attr)[3].store_as_external_file(filename)
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
    pst.write(os.path.join(ws,"freyberg6.pst"))
    pyemu.os_utils.run("pestpp-ies freyberg6.pst",cwd=ws)

    pst.control_data.noptmax = -1
    pst.write(os.path.join(ws, "freyberg6.pst"))
    pyemu.os_utils.start_workers(ws,"pestpp-ies","freyberg6.pst",num_workers=10,master_dir="master_ies")


def _write_instuctions(ws):
    obs_df = pd.read_csv(os.path.join(ws,"heads.csv"), index_col=0)
    names,vals = [],[]
    with open(os.path.join(ws,"heads.csv.ins"),'w') as f:
        f.write("pif ~\nl1\n")
        for i in obs_df.index:
            f.write("l1 w")
            for c in obs_df.columns:
                name = "{0}_{1}".format(c.lower(),i)
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
    df = pd.concat([head_df,flx_df],sort=False)
    return df

def _write_templates(ws):

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
        arr = np.loadtxt(filename)
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

if __name__ == "__main__":
    #prep_mf6_model()
    setup_pest_interface()