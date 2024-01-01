Common repository for H-He equations of state and miscibility curves. The H-He miscibility curves come from [Lorenzen et al. (2009, 2011)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.84.235109) and [Schöttler & Redmer (2018)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.115703). The H-He equations of state come from [Saumon et al. (1995; SCvH)](https://ui.adsabs.harvard.edu/abs/1995ApJS...99..713S/abstract), [Militzer & Hubbard (2013; MH13)](https://iopscience.iop.org/article/10.1088/0004-637X/774/2/148/meta), [Chabrier et al. (2019; CMS19)](https://iopscience.iop.org/article/10.3847/1538-4357/aaf99f/meta), [Chabrier & Debras (2021; CD21)](https://iopscience.iop.org/article/10.3847/1538-4357/abfc48/meta), and [Mazevet et al. (2022; MLS22)](https://www.aanda.org/articles/aa/abs/2022/08/aa35764-19/aa35764-19.html). The CMS and MLS EOSes are supplemented by the work of [Howard et al. (2023a)](https://www.aanda.org/articles/aa/pdf/2023/04/aa44851-22.pdf) to account for the non-ideal entropy and volume interactions they calculated from MH13 and CD21. Moreover, we calculate H-He-Z mixutures using a water EOS from [Haldemann et al. (2020)](https://www.aanda.org/articles/aa/full_html/2020/11/aa38367-20/aa38367-20.html), iron, and post-perovskite EOSes from Jisheng Zhang (private communication). 

To download and use these tables, 

1. ```git lfs install```
2. ```git clone https://github.com/Rob685/hhe_eos_misc.git```
3. ```git submodule init && git submodule update```
4. ```cd eos```
5. ```git checkout main```

Please start with the tutorials instructions to access the equations of state (`eos_tutorial.ipynb`) and their derivatives (`eos_derivatives_turorial.ipynb`), along with the miscibility curves or demixing temperatures (`misc_tutorial.ipynb`).

The following table outlines the quantites provided by the equation of state module:

<img width="603" alt="Screenshot 2023-12-30 at 20 01 17" src="https://github.com/Rob685/eos/assets/48569647/5c18c88b-c64a-425a-ac1b-87cb204fc16c">

In this table, the headers are the independent thermodynamic variables and the quantities are the dependent variables and derivatives. See Tejada Arevalo et al. (2024) for a full description of these quantities.


