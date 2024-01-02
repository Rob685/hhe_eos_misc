from misc import misc_sch, misc_lor, misc_brygoo
from eos import cms_eos
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root

""" This file calls the individual miscibility curves and 
    contains wrappers to get the miscibility curves, the
    immiscibility regions, and the immiscible helium mass
    fractions. 
    
    Authors: Roberto Tejada Arevalo"""

mh = 1
mhe = 4.0026

def x_to_Y(x):
    """Converts from number fraction to mass fraction"""
    return (mhe*x/(mh*(1-x) + mhe*x))

def Y_to_x(Y):
    """Converts from mass fraction to number fraction"""
    return (Y/mhe/(Y/mhe + (1-Y)/mh))

def get_misc_p(logp, Y, misc, delta_T, bry_sigma=0):
    """ This function returns the P-T miscibility 
        curve. Arrays only"""
    if misc=='s':
        pmisc, tmisc = misc_sch.get_misc_curve(logp, Y)
        return pmisc, tmisc+delta_T
    elif misc=='l':
        pmisc, tmisc = misc_lor.get_misc_curve(logp, Y)
        return pmisc, tmisc+delta_T
    elif misc=='b':
        pmisc, tmisc = misc_brygoo.get_misc_curve(logp, sigma=bry_sigma)
        return pmisc, tmisc+delta_T
    else:
        print('misc. curve must be s, or l')

def get_pgap(logp, logt, y, misc = 's', delta_T = 0, sigma=0):
    """ This function returns the P-T miscibility gap points. 
    This should return two points, P1 and P2, and the region
    between P1 and P2 is the H-He immiscbility region."""
    # delta_T in 10^3 K
    interp_t = interp1d(10**(logp-12), 10**(logt-3), kind='linear')
    pmisc, tmisc = get_misc_p(logp, y, misc, delta_T, bry_sigma=sigma)
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
    elif len(idx) == 0:
        return None

def x_err(xval, pval, tval, misc):
    if misc == 's':
        misc_mod = misc_sch
    elif misc == 'l':
        misc_mod = misc_lor
    t = misc_mod.get_t_x_rgi_linear(np.array([pval, xval]).T)
    return (t/tval) - 1

def get_y_misc(logp, logt, misc):
    """ This function returns Y(P, T) for
    the miscibility curves."""
    p, t = 10**(logp-12), 10**(logt-3)
    if np.isscalar(p):
        p, t, = np.array([p]), np.array([t])
        sol = root(x_err, [0.08], tol=1e-8, method='hybr', args=(p, t, misc))
        return float(sol.x)
    sol = root(x_err, np.zeros(len(p))+0.08, tol=1e-10, method='hybr', args=(p, t, misc))
    return x_to_Y(sol.x)


### composition derivatives for implicit scheme ###
def get_dydt_misc(logp, logt, misc_curve, dt=0.1):
    T0 = 10**logt
    T1 = T0*(1+dt)
    Y0 = get_y_misc(logp, logt, misc_curve)
    Y1 = get_y_misc(logp, np.log10(T1), misc_curve)

    return (Y1 - Y0)/(T1 - T0)

def get_dyds_misc(logp, logt, misc_curve, dt=0.1):
    T0 = 10**logt
    T1 = T0*(1+dt)
    Y0 = get_y_misc(logp, logt, misc_curve)
    Y1 = get_y_misc(logp, np.log10(T1), misc_curve)

    S0 = cms_eos.get_s_pt(logp, logt, Y0)
    S1 = cms_eos.get_s_pt(logp, logt, Y1)

    return (Y1 - Y0)/(S1 - S0)

# def get_y_misc(logp, logt, misc_curve='s', delta_T=0, sigma=0):
#     """ This function returns the helium fraction at which any
#         logp, logt would intercept the miscibility curve."""
#     Y_arr = np.linspace(0, 1, len(logp)) # spans all helium fractions to detect the intercept
#     interp_t = interp1d(10**(logp-12), 10**(logt-3), kind='linear')
#     p_prof = 10**(logp-12)
#     p_bound = p_prof[(p_prof > 0.5) & (p_prof < 50)]
#     tnew = interp_t(p_bound)
    
#     y_res = []
#     pmisc_res = []
#     tmisc_res = []
#     for y in Y_arr:
#         y_ = np.full_like(logp, y)
#         pmisc, tmisc = get_misc_p(logp, y_, misc_curve, delta_T=0)

#         tnew = interp_t(pmisc[pmisc < max(10**(logp-12))])
#         tmisc = tmisc[pmisc < max(10**(logp-12))]

#         idx = np.argwhere(np.diff(np.sign(tnew - tmisc))).flatten()
#         if len(idx) == 0:
#             continue
#         elif len(idx) > 0:
#             #plt.plot(pmisc, tmisc)
#             #print(x)
#             y_res.append(y)
#             #pmisc_res.append(pmisc)
#             #tmisc_res.append(tmisc)
#     if len(y_res) == 0:
#         return None
#     else:
#         return y_res[0]

# def get_pgap(logp, logt, y, misc = 's', delta_T = 0, sigma=0, bound='median'):
#     """ This function returns the P-T miscibility gap points. 
#     This should return two points, P1 and P2, and the region
#     between P1 and P2 is the H-He immiscbility region."""
#     # delta_T in 10^3 K
#     interp_t = interp1d(10**(logp-12), 10**(logt-3), kind='cubic')
#     pmisc, tmisc = misc_lor.get_misc_curve(logp, y) # fixing for the lorenzen curve for now...
#     #if delta_T != 0:
#         #pmisc, tmisc = np.insert(pmisc, 0, min(10**(logp-12))), np.insert(tmisc+delta_T, 0, 0)
#     #tmisc += delta_T
#     tnew = interp_t(pmisc[pmisc < max(10**(logp-12))])
#     tmisc = tmisc[pmisc < max(10**(logp-12))]
#     #pnew = pmisc[pmisc < max(10**(logp-12))]

#     idx = np.argwhere(np.diff(np.sign(tnew - tmisc))).flatten()
    
#     return idx