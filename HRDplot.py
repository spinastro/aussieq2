import numpy as np
import glob
from astropy.table import Table
from config import *
import logging
import matplotlib.pyplot as plt
import datetime
from interpolator import *

logger = logging.getLogger(__name__)


def make_age_grid_HRD(age_steps):
    intervals=[]
    intervals.extend(age_steps)
    
    with open(INTERP_PATH+'/YaPSI.age', 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(intervals))


def create_model_HRD(age_steps=[0.001,0.01,0.020],FeH=0.07,Y=0.280,grid="WL11"):
    model_final={"age":[],"mass":[],"logT":[],"logL":[],"logg":[],"Mv":[]}
    make_age_grid_HRD(age_steps)
    make_Fe_Y_grid(FeH,Y,grid)
    os.system('rm '+INTERP_PATH+'/test.out')
    FNULL = open(os.devnull, 'w')
    p = subprocess.Popen('./interp', cwd=INTERP_PATH,stdout=FNULL, stderr=subprocess.STDOUT)
    p.wait()
    model=read_models()
    model_final['age']=model['age'][0]
    model_final['mass']=model['mass'][0]
    model_final['logT']=model['logT'][0]
    model_final['logL']=model['logL'][0]
    model_final['logg']=model['logg'][0]
    model_final['Mv']=model['Mv'][0]
    return(model_final)

def plot_HRD(age_steps=[0.03,0.035,0.04,0.05,0.1,0.2,0.3,0.5],xlim=[6200,5500],ylim=[4.65,4.4],feh=0.0):
    model=create_model_HRD(age_steps=age_steps,FeH=0.05+feh)
    for i in np.arange(len(age_steps)):
        a = model['age'] == age_steps[i]
        plt.plot(10**model['logT'][a],model['logg'][a], label='Age = {i} Gyr'.format(i=age_steps[i]))
    plt.legend(loc='best')
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel(r'T$_{\rm eff}$')
    plt.ylabel('log g')

def plot_HRD_star(Data,star_name,age_steps=[0.03,0.035,0.04,0.05,0.1,0.2,0.3,0.5],xlim=[6200,5500],ylim=[4.65,4.4]):
    f = Data['id'] == star_name
    star=Data[f]
    plt.errorbar(star['teff'],star['logg'],xerr=star['err_teff'],yerr=star['err_logg'],label=star_name, color='black')
    plot_HRD(age_steps,xlim,ylim,feh=Data['feh'])
    plt.savefig('HRD_star.pdf', bbox_inches='tight')
    plt.close()

