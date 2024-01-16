from misc import misc_sch, misc_lor, misc_brygoo
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator as RGI
from scipy.optimize import root

""" 
    This file calls the individual miscibility curves and 
    contains wrappers to get the miscibility curves, the
    immiscibility regions, and the immiscible helium mass
    fractions. 
"""

mh = 1
mhe = 4.0026

def x_to_Y(x):
    # x is the helium number fraction
    # converts number fraction to mass fraction
    return (mhe*x/(mh*(1-x) + mhe*x))

def Y_to_x(Y):
    return (Y/mhe/(Y/mhe + (1-Y)/mh))

def get_misc_p(logp, Y, misc, delta_T, bry_sigma=0, interp='linear'):
    """ This function returns the P-T miscibility 
        curve. Arrays only"""
    if misc=='s':
        pmisc, tmisc = misc_sch.get_misc_curve(logp, Y, misc_interp=interp)
        return pmisc, tmisc+delta_T
    elif misc=='l':
        pmisc, tmisc = misc_lor.get_misc_curve(logp, Y, misc_interp=interp)
        return pmisc, tmisc+delta_T
    elif misc=='b':
        pmisc, tmisc = misc_brygoo.get_misc_curve(logp, sigma=bry_sigma)
        return pmisc, tmisc+delta_T
    else:
        print('misc. curve must be s, or l')

def get_pgap(logp, logt, y, misc = 's', delta_T = 0, sigma=0, interp='linear'):
    """ This function returns the P-T miscibility gap points. 
    This should return two points, P1 and P2, and the region
    between P1 and P2 is the H-He immiscbility region."""
    # delta_T in 10^3 K
    interp_t = interp1d(10**(logp-12), 10**(logt-3), kind='linear')
    pmisc, tmisc = get_misc_p(logp, y, misc, delta_T, bry_sigma=sigma, interp=interp)
    if pmisc is None:
        return None
    tnew = interp_t(pmisc[pmisc < max(10**(logp-12))])
    tmisc = tmisc[pmisc < max(10**(logp-12))]
    #pnew = pmisc[pmisc < max(10**(logp-12))]

    idx = np.argwhere(np.diff(np.sign(tnew - tmisc))).flatten()
    
    if len(idx) > 2:
        return pmisc[idx[-2:]]*1e12
    elif len(idx) == 2:
        if pmisc[int(idx[0])] < 0.1: # avoids intercepts at low pressures due to curve shape
            return None
        else:
            return pmisc[idx]*1e12
    elif len(idx) == 1:
        return np.array([1e12, pmisc[int(idx)]*1e12])
    elif len(idx) == 0:
        return None

def x_err(xval, pval, tval, misc, interp):
    if misc == 's':
        misc_mod = misc_sch
    elif misc == 'l':
        misc_mod = misc_lor

    if interp=='cubic':
        t_new_low = misc_mod.get_t_x_rgi_cubic(np.array([pval[pval < 4.0], xval[pval < 4.0]]).T)
    else:
        t_new_low = misc_mod.get_t_x_rgi_linear(np.array([pval[pval < 4.0], xval[pval < 4.0]]).T)

    t_new_high = misc_mod.get_t_x_rgi_linear(np.array([pval[pval >= 4.0], xval[pval >= 4.0]]).T)

    t = np.concatenate([t_new_low, t_new_high])

    return (t/tval) - 1

def get_y_misc(logp, logt, misc, interp='linear'):
    """ This function inverts T(P, Y) and returns Y(P, T) for
    the miscibility curves."""
    p, t = 10**(logp-12), 10**(logt-3)
    if np.isscalar(p):
        p, t, = np.array([p]), np.array([t])
        sol = root(x_err, [0.08], tol=1e-8, method='hybr', args=(p, t, misc, interp))
        return float(sol.x)
    sol = root(x_err, np.zeros(len(p))+0.08, tol=1e-10, method='hybr', args=(p, t, misc, interp))
    return x_to_Y(sol.x)

logpgrid = np.linspace(6, 14, 100)
logtgrid = np.linspace(2, 4, 150)

y_misc_res_linear, y_misc_res_cubic = np.load('misc/y_misc_pt_tables.npy')

get_y_misc_rgi_linear = RGI((logpgrid, logtgrid), y_misc_res_linear, method='linear', \
            bounds_error=False, fill_value=None)

get_y_misc_rgi_cubic = RGI((logpgrid, logtgrid), y_misc_res_cubic, method='linear', \
            bounds_error=False, fill_value=None)

def get_y_misc_pt(logp, logt, misc_curve='l', interp='linear', tab=True):
    """
    logp in log10 of dyn/cm^2 and logt in log 10 of K.
    """
    if not tab:
        return misc.get_y_misc(logp, logt, misc=misc_curve, interp=interp)
    else:
        if interp=='linear':
            if np.isscalar(logp):
                return float(get_y_misc_rgi_linear((logp, logt)))
            else:
                return get_y_misc_rgi_linear(np.array([logp, logt]).T)
            
        elif interp=='cubic':
            if np.isscalar(logp):
                return float(get_y_misc_rgi_cubic((logp, logt)))
            else:
                return get_y_misc_rgi_cubic(np.array([logp, logt]).T)


### composition derivatives for implicit scheme ###
def get_dydt_misc(logp, logt, misc_curve, dt=0.1, interp='linear', tab=False):
    T0 = 10**logt
    T1 = T0*(1+dt)
    #if interp == 'linear':
    Y0 = get_y_misc_pt(logp, logt, misc_curve, interp, tab)
    Y1 = get_y_misc_pt(logp, np.log10(T1), misc_curve, interp, tab)

    return (Y1 - Y0)/(T1 - T0)

def get_dyds_misc(logp, logt, misc_curve, dt=0.1):
    T0 = 10**logt
    T1 = T0*(1+dt)
    Y0 = get_y_misc(logp, logt, misc_curve)
    Y1 = get_y_misc(logp, np.log10(T1), misc_curve)

    S0 = cms_eos.get_s_pt(logp, logt, Y0)
    S1 = cms_eos.get_s_pt(logp, logt, Y1)

    return (Y1 - Y0)/(S1 - S0)