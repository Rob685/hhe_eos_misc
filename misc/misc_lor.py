import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import pandas as pd
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline as RBS
from scipy.interpolate import RegularGridInterpolator as RGI

from astropy import units as u
from astropy.constants import k_B, m_e, m_p, m_n
from astropy.constants import G, M_jup, M_earth, R_jup, R_earth
from astropy.constants import u as amu
import pdb
from tqdm import tqdm

"""This file reads the Lorenzen et al (2009, 2011) immiscibility data.
    The original data was scrapped from their papers, and the files are
    different for each isobar. The files contain the demixing temperature
    as a function of helium number fraction (x). 
    
    First, I read the files, then I use a regular grid interpolator to
    obtain T(P, x), or the demixing temperature as a function of pressure
    and helium number fraction. 
    
    Then, I return the miscibility curves in the adequate pressure range 
    with the get_misc_curve function."""

erg_to_kbbar = (u.erg/u.Kelvin/u.gram).to(k_B/amu)
MJ_to_kbbar = (u.MJ/u.Kelvin/u.kg).to(k_B/amu)
dyn_to_bar = (u.dyne/(u.cm)**2).to('bar')
erg_to_MJ = (u.erg/u.Kelvin/u.gram).to(u.MJ/u.Kelvin/u.kg)

mp = amu.to('g') # grams
kb = k_B.to('erg/K') # ergs/KG=G.to('cm^3/g s^2').value # Newton's constant
R_jup = R_jup.to('cm').value
M_jup = M_jup.to('g').value

l1 = pd.read_csv('misc/lorenzen_misc_data/lor_1Mbar_top.csv', names=['x', 'T'])
l1['T'] /= 1000
l2 = pd.read_csv('misc/lorenzen_misc_data/lor_2Mbar_top.csv', names=['x', 'T'])
l2['T'] /= 1000
l4 = pd.read_csv('misc/lorenzen_misc_data/lorenzen_4mbar.csv').sort_values('x')
l4['T'] /= 1000
l10 = pd.read_csv('misc/lorenzen_misc_data/lorenzen_10mbar.csv').sort_values('x')
l10['T'] /= 1000
l24 = pd.read_csv('misc/lorenzen_misc_data/lorenzen_24mbar.csv').sort_values('x')
l24['T'] /= 1000

mh = 1
mhe = 4.0026

he_fracs = [l1, l2, l4, l10, l24]
pressures = np.array([1, 2, 4, 10, 24])

def x_to_Y(x):
    return (mhe*x/(mh*(1-x) + mhe*x))

def Y_to_x(Y):
    return (Y/mhe/(Y/mhe + (1-Y)/mh))

interp_1Mbar = interp1d(l1['x'], l1['T'], kind='linear', fill_value='extrapolate')
interp_2Mbar = interp1d(l2['x'], l2['T'], kind='linear', fill_value='extrapolate')
interp_4Mbar = interp1d(l4['x'], l4['T'], kind='linear', fill_value='extrapolate')
interp_10Mbar = interp1d(l10['x'], l10['T'], kind='linear', fill_value='extrapolate')
interp_24Mbar = interp1d(l24['x'], l24['T'], kind='linear', fill_value='extrapolate')

# The Lorenzen curves are non-monotonic at low pressures (1 Mbar)
# Therefore, I limit the the curves to within 0.4 in number fraction (~0.7 in mass fraction)
xgrid = np.arange(np.min(l1['x']), np.max(l1['x']), 0.001) 

# Each demixing temperature range is different, so I interpolate across all isobars here
interp_list = [interp_1Mbar, interp_2Mbar, interp_4Mbar, interp_10Mbar, interp_24Mbar]

T_res = []
for i, p in enumerate(pressures):

    t = interp_list[i](xgrid)
    T_res.append(t)

get_t_x_rgi_cubic = RGI((pressures, xgrid), T_res, method='cubic', bounds_error=False, fill_value=None)
get_t_x_rgi_linear = RGI((pressures, xgrid), T_res, method='linear', bounds_error=False, fill_value=None)

def get_misc_curve(logp, Y, misc_interp='linear'):
    x = Y_to_x(Y)
    p_prof = 10**(logp-12)

    p0_cut = 1.0
    p1_cut = 4.0 
    p2_cut = 50.0
    if misc_interp=='cubic':
        t_new_low = get_t_x_rgi_cubic(np.array([p_prof[(p_prof > p0_cut) & (p_prof < p1_cut)], x[(p_prof > p0_cut) & (p_prof <= p1_cut)]]).T)
    else:
        t_new_low = get_t_x_rgi_linear(np.array([p_prof[(p_prof > p0_cut) & (p_prof < p1_cut)], x[(p_prof > p0_cut) & (p_prof <= p1_cut)]]).T) # smoothness at lower pressures
    t_new_high = get_t_x_rgi_linear(np.array([p_prof[(p_prof > p1_cut) & (p_prof < p2_cut)], x[(p_prof > p1_cut) & (p_prof < p2_cut)]]).T) # flat at higher pressures

    t_new = np.concatenate([t_new_low, t_new_high])

    p_misc = p_prof[(p_prof > p0_cut) & (p_prof < p2_cut)]

    t_new = np.insert(t_new, 0, 0)
    p_misc = np.insert(p_misc, 0, 0)

    t_interp_smooth = interp1d(p_misc, t_new, kind=misc_interp)

    return p_prof[p_prof < p2_cut], t_interp_smooth(p_prof[p_prof < p2_cut])