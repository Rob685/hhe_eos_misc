import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline as RBS

from astropy import units as u
from astropy.constants import k_B, m_e, m_p, m_n
from astropy.constants import G, M_jup, M_earth, R_jup, R_earth
from astropy.constants import u as amu
kb = k_B.to('eV / K').value

"""This file provides the miscibility curves from
    Pfaffenzeller et al. (1995)."""

pvals = np.linspace(0.1, 24, 50)

mh = 1
mhe = 4
def x_to_Y(x):
    # x is the helium number fraction
    # converts number fraction to mass fraction
    #return (mhe*x/(mh*(1-x) + mhe*x)).value
    return (mhe*x/(mh*(1-x) + mhe*x))

def Y_to_x(Y):
    #return (mh*Y/(mhe*(1-Y) + mh*Y))
    return (Y/mhe/(Y/mhe + (1-Y)/mh))

# for the 1 and 2 Mbar plots...
def cut(x, slope, b):
    return slope*x + b
#fits from the 2004 paper:
def fortney_hubbard_fit(x, p):
    c1 = 0.234
    c2 = 0.315
    c3 = 4.215
    return -c1*np.log10(p) + c2*np.log10(x) + c3

def A(P):
    return 1*(P/4)**(0.387)
def x(T, P):
    return np.exp(-A(P)/(kb*T))
def T_pf(x, P):
    return -A(P)/(np.log(x)*kb)

xarr = np.linspace(0.001, 0.65, 100)

t_misc_pfa = []
x_misc_pfa = []

for p in pvals: # same pressure range as the lorenzen above to compare
    t_arr = T_pf(xarr, p)/1e3
    x_misc_pfa.append(xarr)
    t_misc_pfa.append(t_arr)

get_t_x_pfa = RBS(pvals, np.array(x_misc_pfa)[0], t_misc_pfa)

def get_misc_curve(Y):
    x = Y_to_x(Y)

    p_vals = np.linspace(0.01, 24, len(Y))
    
    t_new = get_t_x_pfa.ev(p_vals, x)
    tinterp = interp1d(p_vals, t_new, kind='linear', fill_value='extrapolate')
    
    p = np.logspace(-3, np.log10(24), 500)
    t = tinterp(p)
    return p, t