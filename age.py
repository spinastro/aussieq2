import numpy as np
import glob
from astropy.table import Table
from config import *
import logging
from scipy.integrate import simps
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import datetime


logger = logging.getLogger(__name__)

def list_models():
    models=glob.glob(GRIDS_PATH+'/*.fits')
    model_names=[]
    for i in np.arange(len(models)):
        model_names.append(models[i][models[i].find(GRIDS_PATH)+len(GRIDS_PATH)+1:models[i].find('.fits')])
        print model_names[i]
    return(model_names)


def import_models(model_file):
    m=Table.read(GRIDS_PATH+'/'+model_file+'.fits')
    return(m)

def get_isochrone_points(model,param,nsigma=5):
    print param['age_lim'][0]
    if param['type'] == "MS":
        f = (model["logT"] > np.log10(param['teff']-nsigma*param['err_teff'])) & \
            (model["logT"] < np.log10(param['teff']+nsigma*param['err_teff'])) & \
            (model["logg"] > param['logg']-nsigma*param['err_logg']) & \
            (model["logg"] < param['logg']+nsigma*param['err_logg']) & \
            (model["FeH"] > param['feh']-nsigma*param['err_feh']) & \
            (model["FeH"] < param['feh']+nsigma*param['err_feh']) & \
            (model["age"] > param['age_lim'][0])
    if param['type'] == "PMS":
        f = (model["logT"] > np.log10(param['teff']-nsigma*param['err_teff'])) & \
            (model["logT"] < np.log10(param['teff']+nsigma*param['err_teff'])) & \
            (model["logg"] > param['logg']-nsigma*param['err_logg']) & \
            (model["logg"] < param['logg']+nsigma*param['err_logg']) & \
            (model["FeH"] > param['feh']-nsigma*param['err_feh']) & \
            (model["FeH"] < param['feh']+nsigma*param['err_feh']) & \
            (model["age"] < param['age_lim'][0])
    iso_points={"age":[],"mass":[],"logT":[],"logL":[],"logg":[],"FeH":[],"Mv":[]}
    iso_points['age']=model['age'][f]
    iso_points['mass']=model['mass'][f]
    iso_points['logT']=model['logT'][f]
    iso_points['logL']=model['logL'][f]
    iso_points['logg']=model['logg'][f]
    iso_points['Mv']=model['Mv'][f]
    iso_points['FeH']=model['FeH'][f]
    return(iso_points)

def read_stars(file):
    Stars = Table.read(file, format='ascii.csv')
    
    if ('plx' in Stars.colnames) and ('v' in Stars.colnames):
        Stars['Mv'] = Stars['v'] - 5 * np.log10(1000./Stars['plx']) + 5.
        Stars['err_Mv'] = np.sqrt(Stars['err_v']**2 + (np.log10(np.exp(1))**2)*25*(Stars['err_plx']/Stars['plx'])**2)
    
    if "type" not in Stars.keys():
        n=len(Stars['teff'])
        Stars["type"]=['MS']*n

    if "age_lim" not in Stars.keys():
        n=len(Stars['teff'])
        aaa=np.zeros(n)+0.2
        Stars["age_lim"]=aaa
    return(Stars)


def solve_one_star(model,Data,Star,solve_par='logg',nsigma=5, feh_offset=0.05,smooth_age=10.,smooth_mass=0.):
    
    f = Data['id'] == Star
    param=Data[f]
    param['feh']=param['feh']+feh_offset
    if param['err_teff'] == 0.:
        param['err_teff']=5
        param['err_logg']=0.01
        param['err_feh']=0.01
    iso_points=get_isochrone_points(model,param,nsigma)
    
    if len(iso_points["age"]) == 0:
        logging.warning('Could not get any isochrone points.')
        return None
    available_solv_par = ["logg", "plx", "both"]

    if solve_par not in available_solv_par:
        logging.warning('Solve par must be one of %s' %(available_solv_par))
        return None
    print 'Using {0} Y2 isochrone points - {1:3s} models'.format(len(iso_points['age']),param['type'][0])
    values = ['mp', 'll1s', 'ul1s', 'll2s', 'ul2s', 'mean', 'std']
    print "                 {0:2s} |   {1:4s} -   {2:4s} |   {3:4s} -   {4:4s} |   {5:4s} +/-    {6:3s}".format(values[0],values[1],values[2],values[3],values[4],values[5],values[6],)

    iso_points['teff'] = 10**iso_points['logT']
    iso_points['r'] = 10**(0.5*(np.log10(iso_points['mass'])-iso_points['logg']+4.437))

    if solve_par == 'logg':
        prob = np.exp(-1*((iso_points['teff']-param['teff'])/(1.414214*param['err_teff']))**2)* \
            np.exp(-1*((iso_points['logg']-param['logg'])/(1.414214*param['err_logg']))**2)* \
            np.exp(-1*((iso_points['FeH']-param['feh'])/(1.414214*param['err_feh']))**2)

    if solve_par == 'plx':
        prob = np.exp(-1*((iso_points['teff']-param['teff'])/(1.414214*param['err_teff']))**2)* \
            np.exp(-1*((iso_points['Mv']-param['Mv'])/(1.414214*param['err_Mv']))**2)* \
            np.exp(-1*((iso_points['FeH']-param['feh'])/(1.414214*param['err_feh']))**2)

    if solve_par == 'both':
        prob = np.exp(-1*((iso_points['teff']-param['teff'])/(1.414214*param['err_teff']))**2)* \
        np.exp(-1*((iso_points['Mv']-param['Mv'])/(1.414214*param['err_Mv']))**2)* \
        np.exp(-1*((iso_points['logg']-param['logg'])/(1.414214*param['err_logg']))**2)* \
        np.exp(-1*((iso_points['FeH']-param['feh'])/(1.414214*param['err_feh']))**2)

    #age
    if param['type'][0] == 'MS':
        ages = 0.1+np.arange(150)*0.1
    else:
        ages = 0.001+np.arange(200)*0.001
    pdf_age_x = ages[np.logical_and(ages >= min(iso_points['age'])-0.02,ages <= max(iso_points['age'])+0.02)]
    pdf_age_y, pdf_age_y_smooth, Star_yyage = pdf(pdf_age_x, iso_points, prob, 'age', smooth_window_len=smooth_age)
    Star_yypdf_age = {'x': pdf_age_x, 'y': pdf_age_y, 'ys': pdf_age_y_smooth}
    
    #mass
    masses = 0.4+np.arange(211)*0.01
    pdf_mass_x = masses[np.logical_and(masses >= min(iso_points['mass'])-0.02,masses <= max(iso_points['mass'])+0.02)]
    pdf_mass_y, pdf_mass_y_smooth, Star_yymass = pdf(pdf_mass_x, iso_points, prob, 'mass', smooth_window_len=smooth_mass)
    Star_yypdf_mass = {'x': pdf_mass_x, 'y': pdf_mass_y, 'ys': pdf_mass_y_smooth}

    #logg
    if solve_par == 'plx':
        loggs = np.arange(501)*0.01
        pdf_logg_x = loggs[np.logical_and(loggs >= min(iso_points['logg'])-0.05,loggs <= max(iso_points['logg'])+0.05)]
        pdf_logg_y, pdf_logg_y_smooth, Star_yylogg = pdf(pdf_logg_x, iso_points, prob, 'logg', smooth_window_len=0)
        Star_yypdf_logg = {'x': pdf_logg_x, 'y': pdf_logg_y, 'ys': pdf_logg_y_smooth}

    if param['type'][0] == 'MS':
        xlim=[0,15]
    else:
        xlim=[0,0.2]
    make_plot_age(Data['id'][f][0], Star_yyage,Star_yypdf_age['x'],Star_yypdf_age['ys'],xlim)
    make_plot_mass(Data['id'][f][0], Star_yymass,Star_yypdf_mass['x'],Star_yypdf_mass['ys'])
    if solve_par == 'plx':
        make_plot_logg(Data['id'][f][0], Star_yylogg,Star_yypdf_logg['x'],Star_yypdf_logg['ys'])
    return(Star_yyage,Star_yymass)


def pdf(pdf_x, iso_points, prob, par, smooth_window_len):
    '''Calculates a probability distribution function (PDF) for parameter par
        given the x-values for the PDF, the isochrone points ips, and their
        probability. Return PDF and smoothed PDF (using smooth_window_len) if
        possible (otherwise returns two non-smoothed PDFs), as well as a stats
        dictionary with mean, std, most probable value, etc.
        '''
    
    values = ['mp', 'll1s', 'ul1s', 'll2s', 'ul2s', 'mean', 'std']

    dx = 0.5*(pdf_x[1] - pdf_x[0])
    pdf_y = []
    for x in pdf_x:
        pdf_y.append(sum(prob[np.logical_and(iso_points[par] >= x-dx,iso_points[par] <  x+dx)]))
    pdf_y = np.array(pdf_y)
    pdf_y = pdf_y/simps(pdf_y, pdf_x)

    try:
        pdf_y_smooth = smooth(pdf_y, smooth_window_len)
        pdf_y_smooth = pdf_y_smooth/simps(pdf_y_smooth, pdf_x)
    except:
        pdf_y_smooth = pdf_y
        logger.warning('Unable to smooth '+par+' PDF.')

    stats = get_stats(pdf_x, pdf_y_smooth)

    if stats['most_probable']:
        print "{0:10s} = {1:6.3f} | {2:6.3f} - {3:6.3f} | "\
            "{4:6.3f} - {5:6.3f} | {6:6.3f} +/- {7:6.3f}"\
            .format(par,stats['most_probable'],
                    stats['lower_limit_1sigma'],
                    stats['upper_limit_1sigma'],
                    stats['lower_limit_2sigma'],
                    stats['upper_limit_2sigma'],
                    stats['mean'], stats['std'])
    else:
        print "{0:10s} =        |        -        |  "\
            "      -        | {1:6.3f} +/- {2:6.3f}"\
                .format(par, stats['mean'], stats['std'])
        logger.warning("--- Unable to calculate PDF stats for "+par)

    return pdf_y, pdf_y_smooth, stats

def get_stats(pdf_x, pdf_y_smooth):
    stats = {}
    stats['most_probable'] = np.mean(np.array(pdf_x)[pdf_y_smooth == max(pdf_y_smooth)])
    stats['mean'] = simps(pdf_y_smooth*pdf_x, pdf_x)
    stats['std'] = np.sqrt(simps(pdf_y_smooth*(pdf_x-stats['mean'])**2,pdf_x))

    k = pdf_x <= stats['most_probable']
    pdf_y_left = 0.5*pdf_y_smooth[k]/simps(pdf_y_smooth[k], pdf_x[k])
    pdf_x_left = pdf_x[k]
    areas_left = []
    for x in pdf_x_left:
        areas_left.append(simps(pdf_y_left[pdf_x_left <= x],pdf_x_left[pdf_x_left <= x]))
    areas_left = np.array(areas_left)
    if np.mean(areas_left) == 0:
        logger.warning("Left side of distribution is empty")
        stats['most_probable'] = None
        stats['lower_limit_1sigma'] = None
        stats['lower_limit_2sigma'] = None
        stats['upper_limit_1sigma'] = None
        stats['upper_limit_2sigma'] = None
        return stats

    k = pdf_x >= stats['most_probable']
    pdf_y_right = 0.5*pdf_y_smooth[k]/simps(pdf_y_smooth[k], pdf_x[k])
    pdf_x_right = pdf_x[k]
    areas_right = []
    for x in pdf_x_right:
        areas_right.append(simps(pdf_y_right[pdf_x_right <= x],pdf_x_right[pdf_x_right <= x]))
    areas_right = np.array(areas_right)
                             
    try:
        stats['lower_limit_1sigma'] = np.mean(griddata(areas_left, pdf_x_left, 0.158))
        stats['lower_limit_2sigma'] = np.mean(griddata(areas_left, pdf_x_left, 0.022))
        stats['upper_limit_1sigma'] = np.mean(griddata(areas_right, pdf_x_right, 0.341))
        stats['upper_limit_2sigma'] = np.mean(griddata(areas_right, pdf_x_right, 0.477))
    except:
        stats['lower_limit_1sigma'] = -9.999
        stats['lower_limit_2sigma'] = -9.999
        stats['upper_limit_1sigma'] = -9.999
        stats['upper_limit_2sigma'] = -9.999

    return stats

def make_plot_age(starname, stats,pdf_x,pdf_y_smooth,xlim):
    plt.figure(figsize=(7, 4))
    plt.rc("axes", labelsize=15, titlesize=12)
    plt.rc("xtick", labelsize=14)
    plt.rc("ytick", labelsize=14)
    plt.rc("lines", markersize=10, markeredgewidth=2)
    plt.rc("lines", linewidth=2)
    plt.rc("xtick.major", size=6, width=1)
    plt.rc("ytick.major", size=6, width=1)
    plt.xlim(xlim)
    plt.xlabel('Age (Gyr)')
    plt.ylabel('Probability density')
    k2 = np.logical_and(pdf_x >= stats['lower_limit_2sigma'],pdf_x <= stats['upper_limit_2sigma'])
    k1 = np.logical_and(pdf_x >= stats['lower_limit_1sigma'],pdf_x <= stats['upper_limit_1sigma'])
    
    plt.fill_between(pdf_x[k2], 0 , pdf_y_smooth[k2],color='0.8', hatch="/")
    plt.fill_between(pdf_x[k1], 0 , pdf_y_smooth[k1],color='0.6', hatch="X")
    plt.plot([stats['most_probable'], stats['most_probable']],[0, max(pdf_y_smooth)], 'g--')
    plt.plot(pdf_x, pdf_y_smooth, 'g')
    plt.text(0.92*plt.xlim()[1], 0.86*plt.ylim()[1], starname,horizontalalignment='right', size=16)
    plt.savefig(starname+'_age.pdf', bbox_inches='tight')
    plt.close()

def make_plot_mass(starname, stats,pdf_x,pdf_y_smooth):
    plt.figure(figsize=(7, 4))
    plt.rc("axes", labelsize=15, titlesize=12)
    plt.rc("xtick", labelsize=14)
    plt.rc("ytick", labelsize=14)
    plt.rc("lines", markersize=10, markeredgewidth=2)
    plt.rc("lines", linewidth=2)
    plt.rc("xtick.major", size=6, width=1)
    plt.rc("ytick.major", size=6, width=1)
    plt.xlim([0.5,1.5])
    plt.xlabel(r'Mass (M$_{\odot}$)')
    plt.ylabel('Probability density')
    k2 = np.logical_and(pdf_x >= stats['lower_limit_2sigma'],pdf_x <= stats['upper_limit_2sigma'])
    k1 = np.logical_and(pdf_x >= stats['lower_limit_1sigma'],pdf_x <= stats['upper_limit_1sigma'])
    
    plt.fill_between(pdf_x[k2], 0 , pdf_y_smooth[k2],color='0.8', hatch="/")
    plt.fill_between(pdf_x[k1], 0 , pdf_y_smooth[k1],color='0.6', hatch="X")
    plt.plot([stats['most_probable'], stats['most_probable']],[0, max(pdf_y_smooth)], 'g--')
    plt.plot(pdf_x, pdf_y_smooth, 'g')
    plt.text(0.92*plt.xlim()[1], 0.86*plt.ylim()[1], starname,horizontalalignment='right', size=16)
    plt.savefig(starname+'_mass.pdf', bbox_inches='tight')
    plt.close()

def make_plot_logg(starname, stats,pdf_x,pdf_y_smooth):
    plt.figure(figsize=(7, 4))
    plt.rc("axes", labelsize=15, titlesize=12)
    plt.rc("xtick", labelsize=14)
    plt.rc("ytick", labelsize=14)
    plt.rc("lines", markersize=10, markeredgewidth=2)
    plt.rc("lines", linewidth=2)
    plt.rc("xtick.major", size=6, width=1)
    plt.rc("ytick.major", size=6, width=1)
    plt.xlim([3.5,5.0])
    plt.xlabel('log g (dex)')
    plt.ylabel('Probability density')
    k2 = np.logical_and(pdf_x >= stats['lower_limit_2sigma'],pdf_x <= stats['upper_limit_2sigma'])
    k1 = np.logical_and(pdf_x >= stats['lower_limit_1sigma'],pdf_x <= stats['upper_limit_1sigma'])
    
    plt.fill_between(pdf_x[k2], 0 , pdf_y_smooth[k2],color='0.8', hatch="/")
    plt.fill_between(pdf_x[k1], 0 , pdf_y_smooth[k1],color='0.6', hatch="X")
    plt.plot([stats['most_probable'], stats['most_probable']],[0, max(pdf_y_smooth)], 'g--')
    plt.plot(pdf_x, pdf_y_smooth, 'g')
    plt.text(0.92*plt.xlim()[1], 0.86*plt.ylim()[1], starname,horizontalalignment='right', size=16)
    plt.savefig(starname+'_logg.pdf', bbox_inches='tight')
    plt.close()


def solve_all(input_file, output_file, solve_par, model, nsigma=5,feh_offset=0.07, smooth_age=15,smooth_mass=5):
    print '------------------------------------------------------'
    print 'Initializing ...'
    start_time = datetime.datetime.now()
    print '- Date and time: '+start_time.strftime('%d-%b-%Y, %H:%M:%S')
    print '- Star data: '+input_file
    print '------------------------------------------------------'
    
    Data=read_stars(input_file)

    fout = open(output_file, 'wb')
    pars = ['age', 'mass']
    if solve_par == 'plx':
        pars.append('logg')
    values = ['mp', 'll1s', 'ul1s', 'll2s', 'ul2s', 'mean', 'std']

    hd = 'id'
    for par in pars:
        for value in values:
            hd += ','+par+'_'+value
    fout.write(hd+'\n')

    for star_id in Data['id']:
        print ''
        print '*'*len(star_id)
        print star_id
        print '*'*len(star_id)
        try:
            solution=solve_one_star(model,Data,star_id,solve_par,nsigma,feh_offset,smooth_age,smooth_mass)
        except:
            print 'Unable to find isochrone parameters.'
            print 'Input data might be missing or are too far from valid isochrone points.'
            string = "{0}".format(star_id)
            for par in pars:
                string += ",,,,,,,"
            fout.write(string+"\n")
            continue
        string = "{0}".format(star_id)
        for par in [0,1]:
            keys = ['most_probable',
                    'lower_limit_1sigma', 'upper_limit_1sigma',
                    'lower_limit_2sigma', 'upper_limit_2sigma']
            try:
                for key in keys:
                    string += ",{0:.3f}".format(solution[par][key])
            except:
                string += ",,,,,"
            try:
                string += ",{0:.3f},{1:.3f}".format(solution[par]['mean'],solution[par]['std'])
            except:
                string += ",,"
        fout.write(string+"\n")
    fout.close()

    print ''
    print '------------------------------------------------------'
    end_time = datetime.datetime.now()
    print '- Date and time: '+end_time.strftime('%d-%b-%Y, %H:%M:%S')
    delta_t = (end_time - start_time).seconds
    hours, remainder = divmod(delta_t, 3600)
    minutes, seconds = divmod(remainder, 60)
    print '- Time elapsed: %sH %sM %sS' % (hours, minutes, seconds)
    print 'Done!'
    print '------------------------------------------------------'
    print ''


def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.
        
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal
        (with the window size) in both ends so that transient parts are minimized
        in the begining and end part of the output signal.
        
        input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
        flat window will produce a moving average smoothing.
        
        output:
        the smoothed signal
        
        example:
        
        t=linspace(-2,2,0.1)
        x=sin(t)+randn(len(t))*0.1
        y=smooth(x)
        
        see also:
        
        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter
        
        TODO: the window parameter could be the window itself if an array instead of a string
        NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
        """
    
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    
    if window_len<3:
        return x
    
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2):-(window_len/2)]
