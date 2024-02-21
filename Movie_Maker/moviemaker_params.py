import matplotlib.pyplot as plt
import numpy as np
from aart_func import *
from params import *  # The file params.py contains all the relevant parameters for the simulations
from astropy import units as u
import kgeo
import image_tools as tls
import subprocess
import scipy.interpolate
from matplotlib import ticker
import importlib
from functools import partial
import final_paper_params as fpp

# brightparams = {
#     "nu0": 230e9,  # 0 nu0
#     "mass": (MMkg * u.kg).to(u.g).value,  # 1 mass
#     "scale_height": .5,  # 2 scale_height
#     "theta_b": 50.0 * (np.pi / 180),  # 3 theta_b
#     "beta": 1.0,  # 4 beta
#     "r_ie": 10.0,  # 5 rie
#     "rb_0": 2,  # 7 rb_0
#     "n_th0": 7.2e4,  # ,4e5 # 8 n_th0
#     "t_e0": 2e11,  # 9 t_e0 1e12
#     "p_dens": -.7,  # 10 p_dens
#     "p_temp": -.84,  # 11 p_temp
#     "nscale": .4  # Scale of Inoisy
# }
# brightparams = fpp.bp_fiducial230
# brightparams = fpp.bp_steeperT230
brightparams = fpp.bp_shallowT230


funckeys = {
    "emodelkey": 0,  # emodelkey Emission Model choice, 0 = thermal ultrarelativistic, 1 = power law
    "bkey": 2,  # bkey
    "nnoisykey": 0,  # nnoisykey Inoisy density. 0 = no noise, 1 = noise
    "tnoisykey": 0,  # tnoisykey Inoisy temperature
    "bnoisykey": 0  # bnoisykey Inoisy magnetic field
}

# Labeling Image--------------------------------------------------------------------------------------------------------
var_label = {
    'nu0': r'$\nu= $',
    'mass': 'BlkHole Mass= ',
    'scale_height': 'Scale Height= ',
    'theta_b': r'$\theta_b= $',
    'beta': r'$\beta= $',
    'r_ie': r'$R_ie= $',
    'rb_0': r'$R_b= $',
    'n_th0': r'$n_{th,0}= $',
    't_e0': r'$T_{e,0}= $',
    'p_dens': r'$p_{dens}= $',
    'p_temp': r'$p_{temp}= $',
    'nscale': r'$n_{scale}$'
}
scale_label = {
    'nu0': 1e9,  # GHz
    'mass': 1e9 * 1.989e33,  # Billion Solar Masses
    'scale_height': 1,  # Rg
    'theta_b': 1,  # Rads
    'beta': 1,
    'r_ie': 1,
    'rb_0': 1,  # Rg
    'n_th0': 1 / 1e6,  # 1/cm^3
    't_e0': 1e9,  # GK
    'p_dens': 1,
    'p_temp': 1,
    'nscale': 1
}
units_label = {
    'nu0': 'GHz',  # GHz
    'mass': r'Billion $M_{\odot}$',  # Billion Solar Masses
    'scale_height': r'$R_g$',  # Rg
    'theta_b': 'Rads',  # Rads
    'beta': '',
    'r_ie': '',
    'rb_0': r'$R_g$',  # Rg
    'n_th0': r'$1/m^{3}$',  # 1/m^3
    't_e0': 'GK',  # GK
    'p_dens': '',
    'p_temp': '',
    'nscale': ''
}

aart_results = "/home/tej/Desktop/Code_Stuff/Repositories/aart_results/"
high_rez = "high_rez/"
low_rez = "low_rez/"
rez = low_rez
# rez = high_rez
