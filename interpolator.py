import numpy as np
import csv
import os
from astropy.table import Table
from config import *
import subprocess
import glob
from tqdm import tqdm

def make_age_grid(age_intervals,age_steps):
    youngest_interval=np.arange(age_intervals[0],age_intervals[1],age_steps[0])
    intermediate_interval=np.arange(age_intervals[1],age_intervals[2],age_steps[1])
    old_interval=np.arange(age_intervals[2],age_intervals[3],age_steps[2])

    intervals=[]
    intervals.extend(youngest_interval)
    intervals.extend(intermediate_interval)
    intervals.extend(old_interval)

    with open(INTERP_PATH+'/YaPSI.age', 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(intervals))


def make_Fe_Y_grid(FeH,Y,grid):
    if grid == "WL11":
        file_orig="/YaPSI.nml_WL11"
    if grid == "LCB98":
        file_orig="/YaPSI.nml_LCB98"

    with open(INTERP_PATH+file_orig, 'r') as file:
        filedata = file.read()
        filedata = filedata.replace('targtFe=0.065    ', 'targtFe='+str(FeH)+'    ')
        filedata = filedata.replace('targtY=0.278     ', 'targtY='+str(Y)+'    ')

    with open(INTERP_PATH+'/YaPSI.nml', 'w') as file:
        file.write(filedata)


def read_models():
    info={"age":[],"mass":[],"logT":[],"logL":[],"logg":[],"Mv":[]}
    f = open(INTERP_PATH+'/test.out', 'r')
    a=f.readlines()
    lines=np.array(a[3:len(a)])
    check=np.array(map(lambda x: x[35:43].strip(),lines))
    w=np.where((check != 'NaN') & (check != '') & (check != '********'))
    lines=lines[w[0]]
    check=np.array(map(lambda x: x[9:19].strip(),lines))
    w=np.where(check != '**********')
    lines=lines[w[0]]
    check=np.array(map(lambda x: x[43:50].strip(),lines))
    w=np.where(check != '*******')
    lines=lines[w[0]]
    info["age"].append(np.asarray(map(lambda x: float(x[0:9].strip()),lines)))
    info["mass"].append(np.asarray(map(lambda x: float(x[9:19].strip()),lines)))
    info["logT"].append(np.asarray(map(lambda x: float(x[19:27].strip()),lines)))
    info["logL"].append(np.asarray(map(lambda x: float(x[27:35].strip()),lines)))
    info["logg"].append(np.asarray(map(lambda x: float(x[35:43].strip()),lines)))
    info["Mv"].append(np.asarray(map(lambda x: float(x[43:50].strip()),lines)))
    return(info)


def list_models():
    models=glob.glob(GRIDS_PATH+'/*.fits')
    for i in np.arange(len(models)):
        model_name=models[i][models[i].find(GRIDS_PATH)+len(GRIDS_PATH)+1:models[i].find('.fits')]
        print model_name


def remove_model(model):
    os.remove(GRIDS_PATH+'/'+model+'.fits')


def create_model(age_intervals=[0.001,0.200,0.500,15.],age_steps=[0.001,0.01,0.020],Fe_step=0.002,Y=0.280,grid="WL11"):
    model_final={"age":[],"mass":[],"logT":[],"logL":[],"logg":[],"FeH":[],"Mv":[]}
    make_age_grid(age_intervals,age_steps)
    for FeH in tqdm(np.arange(-0.500,0.500,Fe_step), desc="Interpolating models"):
        #    for FeH in np.arange(-0.500,0.500,Fe_step):
        make_Fe_Y_grid(FeH,Y,grid)
        os.system('rm '+INTERP_PATH+'/test.out')
        FNULL = open(os.devnull, 'w')
        p = subprocess.Popen('./interp', cwd=INTERP_PATH,stdout=FNULL, stderr=subprocess.STDOUT)
        p.wait()
        model_FeH=read_models()
        model_final['age'].append(np.concatenate(model_FeH['age']))
        model_final['mass'].append(np.concatenate(model_FeH['mass']))
        model_final['logT'].append(np.concatenate(model_FeH['logT']))
        model_final['logL'].append(np.concatenate(model_FeH['logL']))
        model_final['logg'].append(np.concatenate(model_FeH['logg']))
        model_final['Mv'].append(np.concatenate(model_FeH['Mv']))
        n_lines=len(np.concatenate(model_FeH['age']))
        model_final['FeH'].append(np.zeros(n_lines)+FeH)
    model_final["age"]=np.concatenate(model_final["age"])
    model_final["mass"]=np.concatenate(model_final["mass"])
    model_final["logT"]=np.concatenate(model_final["logT"])
    model_final["logL"]=np.concatenate(model_final["logL"])
    model_final["logg"]=np.concatenate(model_final["logg"])
    model_final["Mv"]=np.concatenate(model_final["Mv"])
    model_final["FeH"]=np.concatenate(model_final["FeH"])
    t=Table(model_final)
    n_models=len(model_final["age"])
    t.write(GRIDS_PATH+'/'+grid+'_Y'+str(Y)+'_models'+str(n_models)+'.fits', format='fits',overwrite=True)

#yapsi.interpolator.create_model(age_steps=[0.05,0.1,1.],Fe_step=0.2,Y=0.280,grid="WL11")
