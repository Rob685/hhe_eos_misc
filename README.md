# Python Interface for Hydrogen-Helium Equations for Planetary Evolution Models [![DOI](https://zenodo.org/badge/737449342.svg)](https://zenodo.org/doi/10.5281/zenodo.10659198)

# Start Here

This is a common repository for H-He equations of state and miscibility curves. The H-He miscibility curves come from [Lorenzen et al. (2009, 2011)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.84.235109) and [Schöttler & Redmer (2018)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.115703). The H-He equations of state come from [Saumon et al. (1995; SCvH)](https://ui.adsabs.harvard.edu/abs/1995ApJS...99..713S/abstract), [Militzer & Hubbard (2013; MH13)](https://iopscience.iop.org/article/10.1088/0004-637X/774/2/148/meta), [Chabrier et al. (2019; CMS19)](https://iopscience.iop.org/article/10.3847/1538-4357/aaf99f/meta), [Chabrier & Debras (2021; CD21)](https://iopscience.iop.org/article/10.3847/1538-4357/abfc48/meta), and [Mazevet et al. (2022; MLS22)](https://www.aanda.org/articles/aa/abs/2022/08/aa35764-19/aa35764-19.html). The CMS and MLS EOSes are supplemented by the work of [Howard et al. (2023a)](https://www.aanda.org/articles/aa/pdf/2023/04/aa44851-22.pdf) to account for the non-ideal entropy and volume interactions they calculated from MH13 and CD21. Moreover, we calculate H-He-Z mixutures using a water EOS from [Haldemann et al. (2020)](https://www.aanda.org/articles/aa/full_html/2020/11/aa38367-20/aa38367-20.html), iron, and post-perovskite EOSes from [Jisheng Zhang](https://astro.uchicago.edu/people/jisheng-zhang.php) (private communication). 

## Download 

1. Follow the steps to install git lfs on https://git-lfs.com/
2. ```git lfs install```
3. ```git clone https://github.com/Rob685/hhe_eos_misc.git```
4. ```cd hhe_eos_misc```
5. ```git submodule init && git submodule update```
6. ```cd eos```
7. ```git checkout main```

## Tutorials

Please start with the [eos](https://github.com/Rob685/hhe_eos_misc/blob/main/eos_tutorial.ipynb) and [miscibility](https://github.com/Rob685/hhe_eos_misc/blob/main/misc_tutorial.ipynb) tutorials instructions. 
For a tutorial on the EOS derivatives, visit the  [derivatives](https://github.com/Rob685/hhe_eos_misc/blob/main/eos_derivatives_tutorial.ipynb) tutorial.

### Importing

To access the EOS tables, import the ```eos``` module and select which EOS table to call:

```
from eos import cms_eos
```
to import the CMS19+HG23 EOS, for example, or import the `mixtures_eos` file and set zero metallicity,

```
from eos import mixtures_eos

logpgrid = np.linspace(6, 14, 100) # from 1 bar to 100 Mbar
z_0 = np.zeros(len(logpgrid))

...
```
The `mixtures_eos` file provides H-He-Z mixtures with the AQUA tables. A zero metallicity array input yields the H-He EOS, chosen by the `hh_eos` argument in the functions. 

See the eos tutorial for more information.

Similarly, you can import the miscibility module 

```
from misc import misc
```

which will give access to the miscibility curve of one's choice. See the miscibility tutorial for more details.

## EOS Quantities provided

The following table outlines the quantites provided by the equation of state module:

<img width="603" alt="Screenshot 2023-12-30 at 20 01 17" src="https://github.com/Rob685/eos/assets/48569647/5c18c88b-c64a-425a-ac1b-87cb204fc16c">

In this table, the headers are the independent thermodynamic variables and the quantities are the dependent variables and derivatives. See [Tejada Arevalo et al. (2024; submitted to ApJ)](https://arxiv.org/abs/2401.04172) for a full description of these quantities.

## Table ranges

The H-He tables span a wide section of parameter space relevant to gas giant planets, namely:

### CMS19+HG23

- $(P, T)$ tables (original from CMS19 publication):
  - $\log_{10}{P} \in$ [1.0, 23.0] $\rm{dyne}\ \rm{cm}^{-2}$ by $\Delta \log_{10}{P} = 0.05\ \rm{dyne}\ \rm{cm}^{-2}$
  - $\log_{10}{T} \in$ [2.1, 5.0] $K$ by $\Delta \log_{10}{T} = 0.05\ K$
  - $ Y \in $ [0.0, 1.0] by arbitrary amount (volume addition law)

- $(\rho, T)$ tables:
  - $\log_{10}{\rho} \in$ [-5, 1.5] $\rm{g} \rm{cm}^{-3}$ by $\Delta \log_{10}{\rho} = 0.05$ $\rm{g} \rm{cm}^{-3}$
  - $\log_{10}{T} \in$ [2.1, 5.0] $K$ by $\Delta \log_{10}{T} = 0.05\ K$
  - $ Y \in $ [0.05, 1.0] by $\Delta Y = 0.05$
 
- $(S, P)$ tables:
  - $s \in$ [5.5, 9.0] $k_b/\rm{baryon}$ by $\Delta s = 0.05$ $k_b/\rm{baryon}$
  - $\log_{10}{P} \in$ [5.5, 14.0] $\rm{dyne}\ \rm{cm}^{-2}$ by $\Delta \log_{10}{P} = 0.05\ \rm{dyne}\ \rm{cm}^{-2}$
  - $ Y \in $ [0.05, 1.0] by $\Delta Y = 0.05$
 
- $(S, \rho)$ tables:
  - $s \in$ [5.5, 9.0] $k_b/\rm{baryon}$ by $\Delta s = 0.05$ $k_b/\rm{baryon}$
  - $\log_{10}{\rho} \in$ [-4.5, 1.5] $\rm{g} \rm{cm}^{-3}$ by $\Delta \log_{10}{\rho} = 0.05$ $\rm{g} \rm{cm}^{-3}$
  - $ Y \in $ [0.05, 1.0] by $\Delta Y = 0.05$

### MLS22+HG23

- $(P, T)$ tables (original from CMS19 publication):
  - $\log_{10}{P} \in$ [1.0, 23.0] $\rm{dyne}\ \rm{cm}^{-2}$ by $\Delta \log_{10}{P} = 0.05\ \rm{dyne}\ \rm{cm}^{-2}$
  - $\log_{10}{T} \in$ [2.1, 5.0] $K$ by $\Delta \log_{10}{T} = 0.05\ K$
  - $ Y \in $ [0.0, 1.0] by arbitrary amount (volume addition law)

- $(\rho, T)$ tables:
  - $\log_{10}{\rho} \in$ [-5, 1.5] $\rm{g} \rm{cm}^{-3}$ by $\Delta \log_{10}{\rho} = 0.05$ $\rm{g} \rm{cm}^{-3}$
  - $\log_{10}{T} \in$ [2.1, 5.0] $K$ by $\Delta \log_{10}{T} = 0.05\ K$
  - $ Y \in $ [0.05, 1.0] by $\Delta Y = 0.05$
 
- $(S, P)$ tables:
  - $s \in$ [5.25, 10.0] $k_b/\rm{baryon}$ by $\Delta s = 0.05$ $k_b/\rm{baryon}$
  - $\log_{10}{P} \in$ [5.5, 14.0] $\rm{dyne}\ \rm{cm}^{-2}$ by $\Delta \log_{10}{P} = 0.05\ \rm{dyne}\ \rm{cm}^{-2}$
  - $ Y \in $ [0.05, 1.0] by $\Delta Y = 0.05$
 
- $(S, \rho)$ tables:
  - $s \in$ [5.0, 10.0] $k_b/\rm{baryon}$ by $\Delta s = 0.05$ $k_b/\rm{baryon}$
  - $\log_{10}{\rho} \in$ [-5, 1.5] $\rm{g} \rm{cm}^{-3}$ by $\Delta \log_{10}{\rho} = 0.05$ $\rm{g} \rm{cm}^{-3}$
  - $ Y \in $ [0.05, 1.0] by $\Delta Y = 0.05$
 
### SCvH95

- $(P, T)$ tables:
  - $\log_{10}{P} \in$ [6.0, 14.0] $\rm{dyne}\ \rm{cm}^{-2}$ by $\Delta \log_{10}{P} = 0.01\ \rm{dyne}\ \rm{cm}^{-2}$
  - $\log_{10}{T} \in$ [2.1, 5.0] $K$ by $\Delta \log_{10}{T} = 0.05\ K$
  - $ Y \in $ [0.22, 0.98] by $\Delta Y = 0.02$

- $(\rho, T)$ tables:
  - $\log_{10}{\rho} \in$ [-5, 1.5] $\rm{g} \rm{cm}^{-3}$ by $\Delta \log_{10}{\rho} = 0.03$ $\rm{g} \rm{cm}^{-3}$
  - $\log_{10}{T} \in$ [2.1, 5.0] $K$ by $\Delta \log_{10}{T} = 0.05\ K$
  - $ Y \in $ [0.22, 0.98] by $\Delta Y = 0.02$
 
- $(S, P)$ tables:
  - $s \in$ [5.0, 17.0] $k_b/\rm{baryon}$ by $\Delta s = 0.1$ $k_b/\rm{baryon}$
  - $\log_{10}{P} \in$ [6.0, 14.0] $\rm{dyne}\ \rm{cm}^{-2}$ by $\Delta \log_{10}{P} = 0.01\ \rm{dyne}\ \rm{cm}^{-2}$
  - $ Y \in $ [0.22, 0.98] by $\Delta Y = 0.02$
 
- $(S, \rho)$ tables:
  - $s \in$ [5.0, 17.0] $k_b/\rm{baryon}$ by $\Delta s = 0.1$ $k_b/\rm{baryon}$
  - $\log_{10}{\rho} \in$ [-5, 1.5] $\rm{g} \rm{cm}^{-3}$ by $\Delta \log_{10}{\rho} = 0.03$ $\rm{g} \rm{cm}^{-3}$
  - $ Y \in $ [0.22, 0.98] by $\Delta Y = 0.02$
 
### MH13

- $(P, T)$ tables:
  - $\log_{10}{P} \in$ [5.0, 14.0] $\rm{dyne}\ \rm{cm}^{-2}$ by $\Delta \log_{10}{P} = 0.09\ \rm{dyne}\ \rm{cm}^{-2}$
  - $\log_{10}{T} \in$ [2.1, 5.0] $K$ by $\Delta \log_{10}{T} = 0.03\ K$
  - $ Y \in $ [0.246575, 1.0] by arbitrary amount, AVL. 

- $(\rho, T)$ tables:
  - $\log_{10}{\rho} \in$ [-5, 1.5] $\rm{g} \rm{cm}^{-3}$ by $\Delta \log_{10}{\rho} = 0.03$ $\rm{g} \rm{cm}^{-3}$
  - $\log_{10}{T} \in$ [2.1, 5.0] $K$ by $\Delta \log_{10}{T} = 0.05\ K$
  - $ Y \in $ [0.246575, 1.0] by arbitrary amount, AVL. 
 
- $(S, P)$ tables:
  - $s \in$ [5.0, 10.0] $k_b/\rm{baryon}$ by $\Delta s = 0.1$ $k_b/\rm{baryon}$
  - $\log_{10}{P} \in$ [5.0, 14.0] $\rm{dyne}\ \rm{cm}^{-2}$ by $\Delta \log_{10}{P} = 0.09\ \rm{dyne}\ \rm{cm}^{-2}$
  - $ Y \in $ [0.246575, 1.0] by $\Delta Y = 0.01$
 
- $(S, \rho)$ tables:
  - $s \in$ [5.0, 10.0] $k_b/\rm{baryon}$ by $\Delta s = 0.1$ $k_b/\rm{baryon}$
  - $\log_{10}{\rho} \in$ [-5, 1.5] $\rm{g} \rm{cm}^{-3}$ by $\Delta \log_{10}{\rho} = 0.03$ $\rm{g} \rm{cm}^{-3}$
  - $ Y \in $ [0.246575, 1.0] by $\Delta Y = 0.01$






