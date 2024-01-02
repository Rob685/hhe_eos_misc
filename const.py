from astropy import units as u
from astropy.constants import u as amu
from astropy.constants import k_B, m_e, m_p, m_n
"""This file contains frequently used constants"""

R_jup = 7.15e9 #6.99e9
M_jup = 1.89914e30
kb = k_B.to('erg/K') # ergs/K
erg_to_kbbar = (u.erg/u.Kelvin/u.gram).to(k_B/amu)
MJ_to_kbbar = (u.MJ/u.Kelvin/u.kg).to(k_B/amu)
dyn_to_bar = (u.dyne/(u.cm)**2).to('bar')
erg_to_MJ = (u.erg/u.Kelvin/u.gram).to(u.MJ/u.Kelvin/u.kg)
kbbar_to_erg = 1/erg_to_kbbar 


msun = 1.9892e33
rsun = 6.9598e10
lsun = 3.8418e33

mearth = 5.9764e27
rearth = 6.37e8

mjup =  1.89914e30 #1.8986112e30 # (15) ! Guillot+2003 Table 3.1

#radiation constant
a = 7.57e-15

# Seidelmann et al. 2007
rjup = 6.9911e9 # volumetric mean radius
drjup = 6e5
drjup_eq = 4e5
rj_eq = rjup_eq = 7.1492e9
rj_pol = rjup_pol = 6.6854e9
drjup_pol = 10e5

msat = 568.34e27
rsat = 58232e5 # volumetric mean radius
drsat = 6e5
rs_eq = rsat_eq = 60268e5 # equatorial
drsat_eq = 4e5
rs_pol = rsat_pol = 54364e5 # polar
drsat_pol = 10e5

G = 6.67408e-8 # this value is NIST / codata reccommended 2014. mesa has 6.67428.
h = 6.62606896e-27
qe = 4.80320440e-10

c = clight = 2.99792458e10
k = kb = 1.3806504e-16
kb_ev = 8.617385e-5
avogadro = 6.02214179e23
rgas = cgas = k * avogadro

amu = 1.660538782e-24
mn = 1.6749286e-24 # neutron mass (g)
mp = 1.6726231e-24
me = 9.10938291e-28
mh = 1.00794 * amu
mhe = 4.002602 * amu
boltz_sigma = sigma_sb = 5.670400e-5
arad = crad = boltz_sigma * 4 / c
au = 1.495978921e13


s_to_Myr = 3.171e-14
s_to_Gyr = 3.171e-17
Myr_to_s = 3.154e+13
Gyr_to_s = 3.154e+16
