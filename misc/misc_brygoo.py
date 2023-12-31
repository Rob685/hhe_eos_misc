import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

def get_misc_curve(logp, sigma):
    # sigma can now be negative
    T = brygoo['T (kK)'] + sigma*brygoo['t_hi_err']

    bry_interp = interp1d(P[:-2], T[:-2], kind='linear', fill_value='extrapolate')

    p_prof = 10**(logp-12)
    tnew = bry_interp(p_prof[(p_prof > 0.1) & (p_prof < 25)])

    return p_prof[(p_prof > 0.1) & (p_prof < 25)], tnew