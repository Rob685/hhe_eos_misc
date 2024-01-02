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

def get_misc_curve(logp, Y):
    """This function returns the demixing temperatures and provides
    a miscibility curve profile given pressure and helium mass fraction
    profiels"""
    x = Y_to_x(Y)
    p_prof = 10**(logp-12)

    p1_cut = 4.0
    p2_cut = 24.0
    t_new_low = get_t_x_rgi_cubic(np.array([p_prof[(p_prof > 1.0) & (p_prof < p1_cut)], x[(p_prof > 1.0) & (p_prof <= p1_cut)]]).T) # smoothness at lower pressures
    t_new_high = get_t_x_rgi_linear(np.array([p_prof[(p_prof > p1_cut) & (p_prof < p2_cut)], x[(p_prof > p1_cut) & (p_prof < p2_cut)]]).T) # flat at higher pressures

    t_new = np.concatenate([t_new_low, t_new_high])

    p_misc = p_prof[(p_prof > 1.0) & (p_prof < 24.0)]

    t_new = np.insert(t_new, 0, 0)
    p_misc = np.insert(p_misc, 0, 0)

    t_interp_smooth = interp1d(p_misc, t_new, kind='cubic')

    return p_prof[p_prof < 24], t_interp_smooth(p_prof[p_prof < 24])

########## OLD CODE, RETURNS UNMOOTH CURVES #############
# for the 1 and 2 Mbar plots...
# def cut(x, slope, b):
#     return slope*x + b

# def t_mix(p, p1, p2, t1, t2):
#     eta1 = (p2 - p)/(p2 - p1)
#     eta2 = 1-eta1
#     tmix = eta1*t1 + eta2*t2
#     return tmix

# # interpolating for any pressure

# def pressure_interp(p):
#     if p in pressures:
#         if p == 1.0:
#             tbl = he_fracs[int(np.where(pressures == p)[0])]
#             line = cut(tbl['x'], 3.0/0.4, 2.90)
#             tbl_above = tbl[(tbl['T'] >= line) & (tbl['x'] < 0.4)]
#             interp1 = interp1d(tbl_above['x'], tbl_above['T'], kind='quadratic')
#             #xarr1 = np.linspace(np.min(tbl_above['x']), np.max(tbl_above['x']), 200)
#             xgrid = np.arange(0.01, 0.4, 0.01) 
            
#             return xgrid, interp1(xgrid)
        
#         elif p == 2.0:
#             tbl = he_fracs[int(np.where(pressures == p)[0])]
#             line = cut(tbl['x'], 3.0/0.4, 1.7)
#             tbl_above = tbl[(tbl['T'] >= line) & (tbl['x'] < 0.4)]
#             #tbl_below = tbl[tbl['T'] <= line]
#             interp1 = interp1d(tbl_above['x'], tbl_above['T'], kind='quadratic')
#             #interp2 = interp1d(tbl_below['x'], tbl_below['T'], kind='linear')
#             xgrid = np.arange(0.01, 0.4, 0.01) 
            
#             return xgrid, interp1(xgrid)
        
#         else:
#             tbl = he_fracs[int(np.where(pressures == p)[0])]
#             tbl = tbl[tbl['x'] < 0.4]
#             interp1 = interp1d(tbl['x'], tbl['T'], kind='quadratic', fill_value='extrapolate')
#             xarr1 = np.linspace(np.min(tbl['x']), np.max(tbl['x']), 200)
#             tnew1 = interp1(xarr1)
#             return xarr1, tnew1
        
#     elif 1 < p < 2: # dealing with non-monotonic data for these regions...
#             p_idx1 = int(np.argwhere(pressures <= p)[-1])
#             p_idx2 = int(np.argwhere(pressures >= p)[0])
#             tbl1 = he_fracs[p_idx1]
#             tbl2 = he_fracs[p_idx2]
#             #pdb.set_trace()
#             line1 = cut(tbl1['x'],  3.0/0.4, 2.90)
#             line2 = cut(tbl2['x'], 3.0/0.4, 1.7)

#             tbl1_above = tbl1[(tbl1['T'] >= line1)  & (tbl1['x'] < 0.4)].drop_duplicates()
#             #tbl1_below = tbl1[tbl1['T'] <= line1]
#             tbl2_above = tbl2[(tbl2['T'] >= line2)  & (tbl2['x'] < 0.4)].drop_duplicates()
#             #tbl2_below = tbl2[tbl2['T'] <= line2]
#             #pdb.set_trace()

#             # tbl1_above = tbl1_above.drop_duplicates()
#             # tbl2_above = tbl2_above.drop_duplicates()

#             interp1_above = interp1d(tbl1_above['x'], tbl1_above['T'], kind='quadratic', fill_value='extrapolate')
#             interp2_above = interp1d(tbl2_above['x'], tbl2_above['T'], kind='quadratic', fill_value='extrapolate')

#             #xarr11 = np.linspace(np.min(tbl1_above['x']), np.max(tbl1_above['x']), 200)
#             #xarr21 = np.linspace(np.min(tbl2_above['x']), np.max(tbl2_above['x']), 200)
#             xgrid = np.arange(0.01, 0.4, 0.01)

#             tnew1_above = interp1_above(xgrid)
#             tnew2_above = interp2_above(xgrid) # sticking to lower ranges for now

#             tmix_above = t_mix(p, pressures[p_idx1], pressures[p_idx2], tnew1_above, tnew2_above)
            

#             #xave1 = t_mix(p, pressures[p_idx1], pressures[p_idx2], xarr11, xarr21)

#             return xgrid, tmix_above
        
#     elif 2 < p < 4:
#         p_idx1 = int(np.argwhere(pressures <= p)[-1])
#         p_idx2 = int(np.argwhere(pressures >= p)[0])
#         tbl1 = he_fracs[p_idx1]
#         tbl2 = he_fracs[p_idx2]
#         line = cut(tbl1['x'], 3.0/0.4, 1.7)
        
#         tbl1_above = tbl1[(tbl1['T'] >= line)  & (tbl1['x'] < 0.4)]
#         tbl2 = tbl2[tbl2['x'] < 0.4]
#         interp1_above = interp1d(tbl1_above['x'], tbl1_above['T'], kind='quadratic', fill_value='extrapolate')
#         interp2 = interp1d(tbl2['x'], tbl2['T'], kind='quadratic', fill_value='extrapolate')
#         #xarr1 = np.linspace(np.min(tbl1_above['x']), np.max(tbl1_above['x']), 200)
#         #xarr2 = np.linspace(np.min(tbl2['x']), np.max(tbl2['x']), 200)
#         xgrid = np.arange(0.01, 0.4, 0.01)
#         tnew1 = interp1_above(xgrid)
#         tnew2 = interp2(xgrid)
#         tmix = t_mix(p, pressures[p_idx1], pressures[p_idx2], tnew1, tnew2)
#         #xave = t_mix(p, pressures[p_idx1], pressures[p_idx2], xarr1, xarr2)
#         return xgrid, tmix
    
#     elif p < 1.0:
#         xgrid = np.arange(0.01, 0.4, 0.01)
#         return xgrid, np.zeros(len(xgrid))
    
#     else:
#         p_idx1 = int(np.argwhere(pressures <= p)[-1])
#         p_idx2 = int(np.argwhere(pressures >= p)[0])
#         tbl1 = he_fracs[p_idx1]
#         tbl2 = he_fracs[p_idx2]
        
#         tbl1 = tbl1[tbl1['x'] < 0.4]
#         tbl2 = tbl2[tbl2['x'] < 0.4]

#         interp1 = interp1d(tbl1['x'], tbl1['T'], kind='linear', fill_value='extrapolate')
#         interp2 = interp1d(tbl2['x'], tbl2['T'], kind='linear', fill_value='extrapolate')

#         xgrid = np.arange(0.01, 0.4, 0.01)
#         #xarr2 = np.linspace(np.min(tbl2['x']), np.max(tbl2['x']), 200)

#         tnew1 = interp1(xgrid)
#         tnew2 = interp2(xgrid)

#         tmix = t_mix(p, pressures[p_idx1], pressures[p_idx2], tnew1, tnew2)
#         #xave = t_mix(p, pressures[p_idx1], pressures[p_idx2], xarr1, xarr2)
#         return xgrid, tmix

# t_misc = []
# x_misc = []

# # pvals = np.linspace(0.1, 23.5, 150)
# pvals = np.logspace(np.log10(0.1), np.log10(23.5), 150)

# for p in pvals:

#     x_arr, t_arr = pressure_interp(p)
#     x_misc.append(x_arr)
#     t_misc.append(t_arr)

# get_t_x = RBS(pvals, np.array(x_misc)[0], t_misc, kx=1, ky=3)

# def get_misc_curve(logp, Y):
    
#     #Y_vals = Y[()]
#     x = Y_to_x(Y)
#     #p_grid = np.linspace(0.5, 24, 300) # for 1-d interp

#     #p_vals = np.linspace(0.5, 24, len(Y))
#     p_prof = 10**(logp-12)
#     # p_vals = logp[p_prof > 0.5]
#     # x_vals = x[p_prof > 0.5]
    
#     t_new = get_t_x.ev(p_prof[(p_prof > 0.5) & (p_prof < 50)], x[(p_prof > 0.5)  & (p_prof < 50)])
#     #tinterp = interp1d(p_vals, t_new, kind='linear', fill_value='extrapolate')
    
#     #p = np.logspace(np.log10(0.5), np.log10(24), Npts)
#     #p = np.linspace(0.1, 24, Npts) # testing for smoothness
#     #t = tinterp(p_vals) # interpolation with few points to get smoothness

#     return p_prof[(p_prof > 0.5) & (p_prof < 50)], t_new