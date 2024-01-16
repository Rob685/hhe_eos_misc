import numpy as np
from scipy.interpolate import interp1d
from astropy import units as u
from astropy.constants import u as amu
import pandas as pd
from astropy.constants import k_B, m_e, m_p, m_n
from scipy.interpolate import RegularGridInterpolator as RGI
import pdb

"""This file reads the Schottler & Redmer (2018)immiscibility data.
    The original data was scrapped from their papers, and the files are
    different for each isobar. The files contain the demixing temperature
    as a function of helium number fraction (x). 
    
    First, I read the files, then I use a regular grid interpolator to
    obtain T(P, x), or the demixing temperature as a function of pressure
    and helium number fraction. 
    
    Then, I return the miscibility curves in the adequate pressure range 
    with the get_misc_curve function."""

# Reading data

s05 = pd.read_csv('misc/schoettler_nonideal_data/schoettler_0_5Mbar.csv', names = ['x', 'T']).sort_values('x')
s1 = pd.read_csv('misc/schoettler_nonideal_data/schoettler_1Mbar.csv', names = ['x', 'T']).sort_values('x')
s1_2 = pd.read_csv('misc/schoettler_nonideal_data/schoettler_1_2Mbar.csv', names = ['x', 'T']).sort_values('x')
s1_5 = pd.read_csv('misc/schoettler_nonideal_data/schoettler_1_5Mbar.csv', names = ['x', 'T']).sort_values('x')
s2 = pd.read_csv('misc/schoettler_nonideal_data/schoettler_2Mbar.csv', names = ['x', 'T']).sort_values('x')
s4 = pd.read_csv('misc/schoettler_nonideal_data/schoettler_4Mbar.csv', names = ['x', 'T']).sort_values('x')
s10 = pd.read_csv('misc/schoettler_nonideal_data/schoettler_10Mbar.csv', names = ['x', 'T']).sort_values('x')
s24 = pd.read_csv('misc/schoettler_nonideal_data/schoettler_24Mbar.csv', names = ['x', 'T']).sort_values('x')

# list of isobars and corresponding pressures
isobars = [s05, s1, s1_2, s1_5, s2, s4, s10, s24]
pressures = np.array([0.5, 1, 1.2, 1.5, 2, 4, 10, 24])

mh = 1
mhe = 4.0026
def x_to_Y(x):
    return (mhe*x/(mh*(1-x) + mhe*x))

def Y_to_x(Y):
    return (Y/mhe/(Y/mhe + (1-Y)/mh))

interp_list = []
for table in isobars:
    interp = interp1d(table['x'], table['T'], kind='linear', fill_value='extrapolate')
    interp_list.append(interp)
    
xgrid = np.arange(np.min(s05['x']), np.max(s05['x']), 0.001)

T_res = []
for i, p in enumerate(pressures):

    t = interp_list[i](xgrid)
    T_res.append(t)

get_t_x_rgi_cubic = RGI((pressures, xgrid), T_res, method='cubic', bounds_error=False, fill_value=None)
get_t_x_rgi_linear = RGI((pressures, xgrid), T_res, method='linear', bounds_error=False, fill_value=None)

def get_misc_curve(logp, Y, misc_interp='linear'):
    x = Y_to_x(Y)
    p_prof = 10**(logp-12)

    p0_cut = 0.5
    p1_cut = 4.0 # works with 4 but has trouble at 5 Gyr, trying it out with 10
    p2_cut = 24.0
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

########## OLD CODE, RETURNS UNMOOTH CURVES #############
# # list of isobars and corresponding pressures
# he_fracs = [s05, s1, s1_2, s1_5, s2, s4, s10, s24]
# pressures = np.array([0.5, 1, 1.2, 1.5, 2, 4, 10, 24])

# mh = 1
# mhe = 4
# def x_to_Y(x):
#     # x is the helium number fraction
#     # converts number fraction to mass fraction
#     #return (mhe*x/(mh*(1-x) + mhe*x)).value
#     return (mhe*x/(mh*(1-x) + mhe*x))

# def Y_to_x(Y):
#     #return (mh*Y/(mhe*(1-Y) + mh*Y))
#     return (Y/mhe/(Y/mhe + (1-Y)/mh))

# he_massfracs = np.array([0.06, 0.13, 0.20, 0.27, 0.35, 0.42, 0.49])
# he_nfracs = Y_to_x(he_massfracs)

# def t_mix(p, p1, p2, t1, t2):
#     eta1 = (p2 - p)/(p2 - p1)
#     eta2 = 1-eta1
#     tmix = eta1*t1 + eta2*t2
#     return tmix

# def pressure_interp(p):
#     pvals = np.array([0.5, 1.0, 1.2, 1.5, 2.0, 4.0, 10.0, 24.0])
#     if p in pvals: 
#         p_l, y_l, t_l = np.load('misc/schoettler_nonideal_data/sr_lowy_{}.npy'.format(p))
#         p_h, y_h, t_h = np.load('misc/schoettler_nonideal_data/sr_highy_{}.npy'.format(p))
#         y_res = np.append(y_l, y_h)
#         t_res = np.append(t_l, t_h)
#         return y_res, t_res
#     else:
#         p1 = pvals[p >= pvals][-1]
#         p2 = pvals[p <= pvals][0]

#         p_l1, y_l1, t_l1 = np.load('misc/schoettler_nonideal_data/sr_lowy_{}.npy'.format(p1))
#         p_h1, y_h1, t_h1 = np.load('misc/schoettler_nonideal_data/sr_highy_{}.npy'.format(p1))

#         p_l2, y_l2, t_l2 = np.load('misc/schoettler_nonideal_data/sr_lowy_{}.npy'.format(p2))
#         p_h2, y_h2, t_h2 = np.load('misc/schoettler_nonideal_data/sr_highy_{}.npy'.format(p2))

#         interp1_l = interp1d(y_l1, t_l1, 
#                       kind='quadratic', fill_value='extrapolate')
#         interp2_l = interp1d(y_l2, t_l2, 
#                       kind='quadratic', fill_value='extrapolate')

#         interp1_h = interp1d(t_h1, y_h1, 
#                       kind='quadratic', fill_value='extrapolate')
#         interp2_h = interp1d(t_h2, y_h2,
#                       kind='quadratic', fill_value='extrapolate')

#         y_arr1_l = np.linspace(np.min(y_l1), np.max(y_l1), 200)
#         y_arr2_l = np.linspace(np.min(y_l2), np.max(y_l2), 200)

#         t_arr1_h = np.linspace(np.min(t_h1), np.max(t_h1), 200)
#         t_arr2_h = np.linspace(np.min(t_h2), np.max(t_h2), 200)

#         tnew1 = interp1_l(y_arr1_l)
#         tnew2 = interp2_l(y_arr2_l)

#         ynew1 = interp1_h(t_arr1_h)
#         ynew2 = interp2_h(t_arr2_h)

#         tmix_l = t_mix(p, p1, p2, tnew1, tnew2)
#         tmix_h = t_mix(p, p1, p2, t_h1, t_h2)
        
#         yave = np.array([np.mean([y1, y2]) for y1, y2 in zip(ynew1, ynew2)])

#         y_res = np.append(y_l1, yave[::-1])
#         t_res = np.append(tmix_l, tmix_h)

#         return y_res, t_res

# def get_crit_y(p, t):
#     y_arr, tinterp = pressure_interp(float(p))
#     # tc = np.max(tinterp)
#     # yc = y_arr[tinterp == np.max(tinterp)]
    
#     t = np.zeros(len(y_arr)) + t
#     idx = np.argwhere(np.diff(np.sign(tinterp - t))).flatten()
    
#     return y_arr[idx], tinterp[idx] 

# def get_crit_t(p):
#     y_arr, tinterp = pressure_interp(float(p))

#     tmax = np.max(tinterp)
#     ymax = y_arr[tinterp == np.max(tinterp)]

#     return ymax, tmax

# pvals = np.logspace(np.log10(0.5), np.log10(24), 150)

# t_misc = []
# x_misc = []
# for p in pvals:
#     x_arr, t_arr = pressure_interp(p)
#     x_misc.append(x_arr)
#     t_misc.append(t_arr)

# get_t_x = RBS(pvals, np.array(x_misc)[0,:], t_misc, kx=1, ky=3)

# def get_misc_curve(logp, Y):
    
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
    #return p, t