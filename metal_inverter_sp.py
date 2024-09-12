import matplotlib.pyplot as plt
import numpy as np
from eos import mixtures_eos
from scipy.optimize import root, root_scalar

erg_to_kbbar = mixtures_eos.erg_to_kbbar

from tqdm import tqdm
from scipy.interpolate import RegularGridInterpolator as RGI
from scipy.optimize import root

cms_eos = mixtures_eos.cms_eos
mls_eos = mixtures_eos.mls_eos
scvh_eos = mixtures_eos.scvh_eos
metals_eos = mixtures_eos.metals_eos
ideal_eos = mixtures_eos.ideal_eos
mh13_eos = mixtures_eos.mh13_eos

def inversion_z(xgrid, ygrid, hegrid, zgrid, basis, xy_eos, z_eos, hg=True):
    sol1 = []
    sol2 = []
    for x in tqdm(xgrid):
        res1_y = []
        res2_y = []
        xarr = np.full_like(zgrid, x)
        for y in ygrid:
            res1_he = []
            res2_he = []
            yarr = np.full_like(zgrid, y)
            for yhe in hegrid:
                hearr = np.full_like(zgrid, yhe)
                if basis == 'sp':
                    res1 = mixtures_eos.get_t_sp(xarr, yarr, hearr, zgrid, \
                                                    hhe_eos=xy_eos, z_eos=z_eos, hg=hg)
                    res2 = mixtures_eos.get_rho_pt(yarr, res1, hearr, zgrid, hhe_eos=xy_eos, z_eos=z_eos, hg=hg)

                elif basis == 'rhot':
                    res1 = mixtures_eos.get_p_rhot(xarr, yarr, hearr, zgrid, \
                                                    hhe_eos=xy_eos, z_eos=z_eos, hg=hg)
                    res2 = mixtures_eos.get_s_pt(res1, yarr, hearr, zgrid, hhe_eos=xy_eos, z_eos=z_eos, hg=hg)

                elif basis == 'srho':
                        try:
                            res1 = mixtures_eos.get_p_srho(xarr, yarr, hearr, zgrid, \
                                                        hhe_eos=xy_eos, z_eos=z_eos, hg=hg)
                            res2 = mixtures_eos.get_t_sp(xarr, res1, hearr, zgrid, hhe_eos=xy_eos, hg=hg)
                                                        #hhe_eos=xy_eos, alg='root', z_eos=z_eos)
                        except ValueError:                        
                            #print(res1)
                            raise Exception('Failed at s = {}, rho = {}, y = {}'.format(xarr[0], yarr[0], hearr[0]))

                elif basis == 'rhop':
                    res1 = mixtures_eos.get_t_rhop(xarr, yarr, hearr, zgrid, \
                                                hhe_eos=xy_eos, alg='root', z_eos=z_eos)
                    res2 = mixtures_eos.get_s_pt(yarr, res1, hearr, zgrid, \
                                                    hhe_eos=xy_eos, alg='root', z_eos=z_eos)

                res1_he.append(res1)
                res2_he.append(res2)

            res1_y.append(res1_he)
            res2_y.append(res2_he)
        
        sol1.append(res1_y)
        sol2.append(res2_y)

    return sol1, sol2

svals_sp_cd = np.arange(0.5, 10.1, 0.1)
logpvals_sp_cd = np.arange(6, 14.1, 0.1)
yvals_sp_cd = np.arange(0.02, 1.0, 0.05)
zvals_sp_cd = np.arange(0.0, 1.05, 0.05) # dense grid of Z

print('Starting S, P process ...')

# I can use the tables produced here to calculate P_srho later since some of the inversions are not converging
# for the P_srho. This is akin to using the P(rho, T) tables for the Tsrho inversions.
# This then by-passes the P(rho,T) inversion.
logt_res_sp, logrho_res_sp = inversion_z(svals_sp_cd, logpvals_sp_cd, yvals_sp_cd, zvals_sp_cd, \
                                        basis='sp', xy_eos='cd', z_eos='fo')

np.save('eos/cd/sp_base_z_forsterite_extended.npy', [logrho_res_sp, logt_res_sp])

print('Finished and saved!')