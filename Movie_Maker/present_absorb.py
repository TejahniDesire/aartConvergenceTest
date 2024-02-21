import os.path
import sys
import subprocess

aartpath = '/home/tej/Desktop/Code_Stuff/Repositories/aart'  # insert path to aart repo
sys.path.append(aartpath)

from matplotlib import ticker
from aart_func import *
from params import *  # The file params.py contains all the relevant parameters for the simulations
from astropy import units as u
import image_tools as tls
import numpy as np

aart_results = "/home/tej/Desktop/Code_Stuff/Repositories/aart_results/"
movie_path = "Movie/"
high_rez = "high_rez/"
low_rez = "low_rez/"
# rez = low_rez
rez = high_rez

vid_path = rez + "video/"
graph_path = rez + "graph/"
temp_path = rez + "temporary/"
speed = 8

parser = argparse.ArgumentParser(description='Movies', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('var', help=
''' 
nu0: (Hertz)
mass: (grams)
scale_height: (rg)
theta_b: (radians)
beta: (dimensionless)
r_ie: (dimensionless)
rb_0: (dimensionless)
n_th0: (1/cm^3)
t_e0: (Kelvin)
p_dens: (dimensionless)
p_temp: (dimensionless) 
nscale: (dimensionless)
''',

                    type=str)

parser.add_argument('start', type=float)
parser.add_argument('stop', help='Inclusive Stop', type=float)
parser.add_argument("step_size", type=float)
parser.add_argument("power_val", type=float,
                    help="Colorbar power norm value. Lower for anticipated higher brightness range")
args = parser.parse_args()
action = {
    "var": args.var,
    "start": args.start,
    "stop": args.stop,
    "step": args.step_size,
    "bright_power": args.power_val
}

'''Reading of the lensing bands----------------------------------'''
fnbands = path + "LensingBands_a_%s_i_%s.h5" % (spin_case, i_case)

print("Reading file: ", fnbands)

h5f = h5py.File(fnbands, 'r')

# Points for the boundary of the BH shadow
alpha_critc = h5f['alpha'][:]
beta_critc = h5f['beta'][:]

# The concave hulls for the lensing bands
hull_0i = h5f['hull_0i'][:]
hull_0e = h5f['hull_0e'][:]
hull_1i = h5f['hull_1i'][:]
hull_1e = h5f['hull_1e'][:]
hull_2i = h5f['hull_2i'][:]
hull_2e = h5f['hull_2e'][:]

# The grid points for each lensing band
supergrid0 = h5f['grid0'][:]
N0 = int(h5f["N0"][0])
mask0 = h5f['mask0'][:]
lim0 = int(h5f["lim0"][0])
supergrid1 = h5f['grid1'][:]
N1 = int(h5f["N1"][0])
mask1 = h5f['mask1'][:]
lim1 = int(h5f["lim1"][0])
supergrid2 = h5f['grid2'][:]
N2 = int(h5f["N2"][0])
mask2 = h5f['mask2'][:]
lim2 = int(h5f["lim2"][0])

h5f.close()

'''Reading Analytical Ray-tracing----------------------------------'''
fnrays = path + "Rays_a_%s_i_%s.h5" % (spin_case, i_case)

print("Reading file: ", fnrays)

h5f = h5py.File(fnrays, 'r')

rs0 = h5f['rs0'][:]
sign0 = h5f['sign0'][:]
t0 = h5f['t0'][:]
phi0 = h5f['phi0'][:]

rs1 = h5f['rs1'][:]
sign1 = h5f['sign1'][:]
t1 = h5f['t1'][:]
phi1 = h5f['phi1'][:]

rs2 = h5f['rs2'][:]
sign2 = h5f['sign2'][:]
t2 = h5f['t2'][:]
phi2 = h5f['phi2'][:]

h5f.close()

'''Computing images----------------------------------'''

args = ' '
cmd1_args = {
    "nu0": '--nu ',
    "mass": '--mass ',
    "scale_height": '--scaleh ',
    "theta_b": '--thetab ',
    "beta": '--beta ',
    "r_ie": '--rie ',
    "rb_0": '--rb0 ',
    "n_th0": '--nth0 ',
    "t_e0": '--te0 ',
    "p_dens": '--pdens ',
    "p_temp": '--ptemp ',
    "nscale": '--nscale ',
}

cmd2_args = {
    "emodelkey": '--emodelkey ',
    "bkey": '--bkey ',
    "nnoisykey": '--nnoisykey ',
    "tnoisykey": '--tnoisykey ',
    "bnoisykey": '--bnoisykey ',
}
# brightparams = [
# 	230e9, # nu0
# 	(MMkg * u.kg).to(u.g).value, # mass
# 	.5, # scale_height
# 	50.0 * (np.pi / 180), # theta_b
# 	1.0, # beta
# 	10.0, # Rie
# 	0, # Bchoice
# 	2, # rb
# 	1.0726e+05, # n_th0
# 	1.2428e+11, # t_e0
# 	-.7, # p_dens
# 	-.84 # p_temp
# ]
brightparams = {
    "nu0": 230e9,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g).value,  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180),  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 2e5,  # ,4e5 # 8 n_th0
    "t_e0": 2e11,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -.84,  # 11 p_temp
    # "p_temp": -.3,  # ptemp higher
    # "p_temp": -1.6, # ptemp lower
    "nscale": .4  # Scale of Inoisy
}

funckeys = {
    "emodelkey": 0,  # emodelkey Emission Model choice, 0 = thermal ultrarelativistic, 1 = power law
    "bkey": 1,  # bkey
    "nnoisykey": 1,  # nnoisykey Inoisy density. 0 = no noise, 1 = noise
    "tnoisykey": 0,  # tnoisykey Inoisy temperature
    "bnoisykey": 1  # bnoisykey Inoisy magnetic field
}

# i = 0
# for arg in cmd1_args:
#     args = args + cmd1_args[arg] + str(brightparams[arg]) + ' '
#     i += 1
#
# i = 0
# for arg in cmd2_args:
#     args = args + cmd2_args[arg] + str(funckeys[arg]) + ' '
#
# # for i in range(len(brightparams)):
# #     args = args + cmd1_args[i] + str(brightparams[i]) + ' '
# aartpath = '/home/tej/Desktop/Code_Stuff/Repositories/aart'

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
# ----------------------------------------------------------------------------------------------------------------------
size = 1000  # array size for radii calcs
num_iterations = int((action["stop"] - action["start"]) / action["step"])
doth5_files = []  # .h5 files
thin_names = []  # thin images to be deleted after movie is made
thick_names = []
thin_I2_names = []
thick_I2_names = []

# Secondary Graphs
x_variable = np.zeros(num_iterations)  # counter for independant variable
janksys_thick = np.ndarray([num_iterations, 4])  # [I0, I1, I2, FullImage]
janksys_thin = np.ndarray([num_iterations, 4])  # [I0, I1, I2, total]

mean_radii_Thin = np.ndarray([num_iterations, 3])  # [I0, I1, I2]
mean_radii_Thick = np.ndarray([num_iterations, 4])  # [I0, I1, I2, FullImage]

radii_I0_Thin = np.zeros(size)
radii_I1_Thin = np.zeros(size)
radii_I2_Thin = np.zeros(size)



radii_FullAbsorption_Thick = np.zeros(size)
radii_I0_Thick = np.zeros(size)
radii_I1_Thick = np.zeros(size)
radii_I2_Thick = np.zeros(size)

for i in range(num_iterations):
    print('Creating Data Set: ' + str(i))
    # Update varing parameter
    brightparams[action["var"]] = action["start"] + i * action["step"]
    x_variable[i] = brightparams[action["var"]]

    for arg in cmd1_args:
        args = args + cmd1_args[arg] + str(brightparams[arg]) + ' '

    for arg in cmd2_args:
        args = args + cmd2_args[arg] + str(funckeys[arg]) + ' '

    subprocess.run(['python3 ' + aartpath + '/radialintensity.py' + args], shell=True)

    # Read created file

    fnrays = aart_results + (
        'path_results/Intensity_a_{}_i_{}_nu_{}_mass_{}_scaleh_{}_thetab_{}_beta_{}_rie_{}_rb_{}_nth0_{}_te0_{}_'
        'pdens_{}_ptemp_{}_nscale_{}_emkey_{}_bkey_{}_nkey_{}_tnkey_{}_bnkey_{}.h5').format(
        spin_case,
        i_case,
        "{:.5e}".format(brightparams["nu0"]),
        "{:.5e}".format(brightparams["mass"]),
        float(brightparams["scale_height"]),
        "{:.3f}".format(brightparams["theta_b"]),
        "{:.2f}".format(float(brightparams["beta"])),
        "{:.1f}".format(float(brightparams["r_ie"])),
        "{:.1f}".format(float(brightparams["rb_0"])),
        "{:.1e}".format(brightparams["n_th0"]),
        "{:.1e}".format(brightparams["t_e0"]),
        float(brightparams["p_dens"]),
        float(brightparams["p_temp"]),
        "{:.1f}".format(brightparams["nscale"]),
        funckeys["emodelkey"],
        funckeys["bkey"],
        funckeys["nnoisykey"],
        funckeys["tnoisykey"],
        funckeys["bnoisykey"]
    )
    doth5_files += [fnrays]

    h5f = h5py.File(doth5_files[i], 'r')
    I0 = h5f['bghts0'][:]
    I1 = h5f['bghts1'][:]
    I2 = h5f['bghts2'][:]

    I2_Absorb = h5f['bghts2_absorbtion'][:]
    I1_Absorb = h5f['bghts1_absorbtion'][:]
    I0_Absorb = h5f['bghts0_absorbtion'][:]
    Absorbtion_Image = h5f['bghts_full_absorbtion'][:]
    # tau2 = h5f['tau2'][:]
    # tau1 = h5f['tau1'][:]
    # tau0 = h5f['tau0'][:]
    # full_profiles0 = h5f['full_profiles0'][:]
    # full_profiles1 = h5f['full_profiles1'][:]
    # full_profiles2 = h5f['full_profiles2'][:]
    # full_profiles_unit = h5f['full_profiles_unit'][:]

    # Thin Radii Calcs--------------------------------------------------------------------------------------------------

    radii_I0_Thin_i, theta = tls.radii_of_theta(I0, size)
    radii_I1_Thin_i, theta = tls.radii_of_theta(I1, size)
    radii_I2_Thin_i, theta = tls.radii_of_theta(I2, size)

    r0_thin = tls.curve_params(theta, radii_I0_Thin_i)
    r1_thin = tls.curve_params(theta, radii_I1_Thin_i)
    r2_thin = tls.curve_params(theta, radii_I2_Thin_i)

    mean_radii_Thin[i, 0] = r0_thin
    mean_radii_Thin[i, 1] = r1_thin
    mean_radii_Thin[i, 2] = r2_thin

    radii_I0_Thin = np.vstack((radii_I0_Thin, radii_I0_Thin_i))
    radii_I1_Thin = np.vstack((radii_I1_Thin, radii_I1_Thin_i))
    radii_I2_Thin = np.vstack((radii_I2_Thin, radii_I2_Thin_i))

    # Thick Radii Calcs-------------------------------------------------------------------------------------------------
    radii_I0_Thick_i, theta = tls.radii_of_theta(I0_Absorb, size)
    radii_I1_Thick_i, theta = tls.radii_of_theta(I1_Absorb, size)
    radii_I2_Thick_i, theta = tls.radii_of_theta(I2_Absorb, size)
    radii_FullAbsorption_Thick_i, theta = tls.radii_of_theta(Absorbtion_Image, size)



    r0_thick = tls.curve_params(theta, radii_I0_Thick_i)
    r1_thick = tls.curve_params(theta, radii_I1_Thick_i)
    r2_thick = tls.curve_params(theta, radii_I2_Thick_i)
    full_thick = tls.curve_params(theta, radii_FullAbsorption_Thick_i)

    mean_radii_Thick[i, 0] = r0_thick
    mean_radii_Thick[i, 1] = r1_thick
    mean_radii_Thick[i, 2] = r2_thick
    mean_radii_Thick[i, 3] = full_thick

    radii_I0_Thick = np.vstack((radii_I0_Thick, radii_I0_Thick_i))
    radii_I1_Thick = np.vstack((radii_I1_Thick, radii_I1_Thick_i))
    radii_I2_Thick = np.vstack((radii_I2_Thick, radii_I2_Thick_i))
    radii_FullAbsorption_Thick = np.vstack((radii_FullAbsorption_Thick, radii_FullAbsorption_Thick_i))

    # Total Flux Calcualtions
    janksys_thin[i,0] = ilp.total_jy(I0, brightparams["nu0"], brightparams["mass"]).value
    janksys_thin[i, 1] = ilp.total_jy(I1, brightparams["nu0"], brightparams["mass"]).value
    janksys_thin[i, 2] = ilp.total_jy(I2, brightparams["nu0"], brightparams["mass"]).value
    janksys_thin[i, 3] = ilp.total_jy(I0 + I1 + I2, brightparams["nu0"], brightparams["mass"]).value

    janksys_thick[i, 0] = ilp.total_jy(I0_Absorb, brightparams["nu0"], brightparams["mass"]).value
    janksys_thick[i, 1] = ilp.total_jy(I1_Absorb, brightparams["nu0"], brightparams["mass"]).value
    janksys_thick[i, 2] = ilp.total_jy(I2_Absorb, brightparams["nu0"], brightparams["mass"]).value
    janksys_thick[i, 3] = ilp.total_jy(Absorbtion_Image, brightparams["nu0"], brightparams["mass"]).value

    h5f.close()

# Black Hole Inner Shadow Calc--------------------------
r_inner = np.load('r_inner_spin_{}_inc_{}.npy'.format(spin_case, i_case))
alphas_inner = np.load('alphas_inner_spin_{}_inc_{}.npy'.format(spin_case,  i_case))
betas_inner = np.load('betas_inner_spin_{}_inc_{}.npy'.format(spin_case,  i_case))

# Black Hole Outer Shadow Calc--------------------------
r_outer = np.load('r_outer_spin_{}_inc_{}.npy'.format(spin_case,  i_case))
alphas_outer = np.load('alphas_outer_spin_{}_inc_{}.npy'.format(spin_case,  i_case))
betas_outer = np.load('betas_outer_spin_{}_inc_{}.npy'.format(spin_case,  i_case))



# ------------------------------------------------------

dim = [16, 8]
xaxis = np.array(x_variable) / scale_label[action['var']]
temp_plots_path = aart_results + movie_path + temp_path
for i in range(num_iterations):
    print('Creating Frame : ' + str(i))

    one_M = ilp.rg_func(brightparams["mass"] * u.g).to(u.m)
    M2uas = np.arctan(one_M.value / dBH) / muas_to_rad

    h5f = h5py.File(doth5_files[i], 'r')
    I0 = h5f['bghts0'][:]
    I1 = h5f['bghts1'][:]
    I2 = h5f['bghts2'][:]

    I2_Absorb = h5f['bghts2_absorbtion'][:]
    I1_Absorb = h5f['bghts1_absorbtion'][:]
    I0_Absorb = h5f['bghts0_absorbtion'][:]
    Absorbtion_Image = h5f['bghts_full_absorbtion'][:]
    # tau2 = h5f['tau2'][:]
    # tau1 = h5f['tau1'][:]
    # tau0 = h5f['tau0'][:]
    # full_profiles0 = h5f['full_profiles0'][:]
    # full_profiles1 = h5f['full_profiles1'][:]
    # full_profiles2 = h5f['full_profiles2'][:]
    # full_profiles_unit = h5f['full_profiles_unit'][:]

    h5f.close()
    subprocess.run(['rm ' + doth5_files[i]], shell=True)

    max_factor = 1.2
    if i == 0:
        # vmax0 = np.max(I0 + I1 + I2) * max_factor
        vmax1 = Absorbtion_Image.max() * max_factor

    vmax0 = np.max(I0 + I1 + I2) * max_factor

    # Thin_FULL---------------------------------------------------------
    fig, ax0 = plt.subplots(figsize=dim, dpi=400)
    # im0 = ax0.imshow(I0 + I1 + I2, origin="lower", cmap="afmhot", extent=[-lim0, lim0, -lim0, lim0],
    #                   norm=matplotlib.colors.PowerNorm(action["bright_power"], vmax=vmax0))
    im0 = ax0.imshow(I0 + I1 + I2, origin="lower", cmap="afmhot", extent=[-lim0, lim0, -lim0, lim0],vmax=vmax0)

    ax0.set_xlim(-10, 10)  # units of M
    ax0.set_ylim(-10, 10)
    ax0.set_xlabel(r"$\alpha$" + " " + r"($\mu as$)")
    ax0.set_ylabel(r"$\beta$" + " " + r"($\mu as$)")

    ax0.set_xticks([-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10], labels=[
        str('{:.3}'.format(-10 * M2uas)),
        str('{:.3}'.format(-7.5 * M2uas)),
        str('{:.3}'.format(-5 * M2uas)),
        str('{:.3}'.format(-2.5 * M2uas)),
        str('{:.3}'.format(0 * M2uas)),
        str('{:.3}'.format(2.5 * M2uas)),
        str('{:.3}'.format(5 * M2uas)),
        str('{:.3}'.format(7.5 * M2uas)),
        str('{:.3}'.format(10 * M2uas))
    ])

    ax0.set_yticks([-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10], labels=[
        str('{:.3}'.format(-10 * M2uas)),
        str('{:.3}'.format(-7.5 * M2uas)),
        str('{:.3}'.format(-5 * M2uas)),
        str('{:.3}'.format(-2.5 * M2uas)),
        str('{:.3}'.format(0 * M2uas)),
        str('{:.3}'.format(2.5 * M2uas)),
        str('{:.3}'.format(5 * M2uas)),
        str('{:.3}'.format(7.5 * M2uas)),
        str('{:.3}'.format(10 * M2uas))
    ])

    colorbar0 = fig.colorbar(im0, fraction=0.046, pad=0.04, format='%.1e', ticks=[
        vmax0 * .8,
        vmax0 * .6,
        vmax0 * .4,
        vmax0 * .2,
        vmax0 * .05
    ],
                             label="Brightness Temperature (K)",
                             ax=ax0
                             )

    ax0.text(-9, 8.5, var_label[action["var"]] + str(round(x_variable[i] / scale_label[action["var"]], 2))
             + ' ' + units_label[action["var"]],fontsize=12, color="w")

    thin_names += [temp_plots_path + 'Thin_Full_Fig_' + str(i)]
    plt.savefig(
        thin_names[i],
        bbox_inches='tight'
    )
    plt.close()

    # Thick_Full Image-----------------------------------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=dim, dpi=400)
    im1 = ax1.imshow(Absorbtion_Image, vmax=vmax1, origin="lower", cmap="afmhot", extent=[-lim0, lim0, -lim0, lim0])

    ax1.set_xlim(-10, 10)  # units of M
    ax1.set_ylim(-10, 10)
    ax1.set_xlabel(r"$\alpha$" + " " + r"($\mu as$)")
    ax1.set_ylabel(r"$\beta$" + " " + r"($\mu as$)")

    ax1.set_xticks([-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10], labels=[
        str('{:.3}'.format(-10 * M2uas)),
        str('{:.3}'.format(-7.5 * M2uas)),
        str('{:.3}'.format(-5 * M2uas)),
        str('{:.3}'.format(-2.5 * M2uas)),
        str('{:.3}'.format(0 * M2uas)),
        str('{:.3}'.format(2.5 * M2uas)),
        str('{:.3}'.format(5 * M2uas)),
        str('{:.3}'.format(7.5 * M2uas)),
        str('{:.3}'.format(10 * M2uas))
    ])

    ax1.set_yticks([-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10], labels=[
        str('{:.3}'.format(-10 * M2uas)),
        str('{:.3}'.format(-7.5 * M2uas)),
        str('{:.3}'.format(-5 * M2uas)),
        str('{:.3}'.format(-2.5 * M2uas)),
        str('{:.3}'.format(0 * M2uas)),
        str('{:.3}'.format(2.5 * M2uas)),
        str('{:.3}'.format(5 * M2uas)),
        str('{:.3}'.format(7.5 * M2uas)),
        str('{:.3}'.format(10 * M2uas))
    ])

    colorbar0 = fig.colorbar(im1, fraction=0.046, pad=0.04, format='%.1e', ticks=[
        vmax1 * .8,
        vmax1 * .6,
        vmax1 * .4,
        vmax1 * .2,
        vmax1 * .05
    ],
                             label="Brightness Temperature (K)",
                             ax=ax1
                             )

    ax1.text(-9, 8.5, var_label[action["var"]] + str(round(x_variable[i] / scale_label[action["var"]], 2))
             + ' ' + units_label[action["var"]],fontsize=12, color="w")

    thick_names += [temp_plots_path + 'Thick_Full_Fig_' + str(i)]
    plt.savefig(
        thick_names[i],
        bbox_inches='tight'
    )
    plt.close()



# FMPEG Stiching--------------------------------------------------------------------------------------------------------
iteration = (
    'Var_{}_Start_{}_Stop_{}_Step_{}_BP_{}_a_{}_i_{}_nu_{}'
    '_mass_{}_scaleh_{}_thetab_{}_beta_{}_rie_{}_rb_{}_nth0_{}_te0_{}'
    '_pdens_{}_ptemp_{}_nscale_{}_keys_{}_{}_{}_{}_{}'
).format(
    action["var"],
    "{:.2e}".format(action["start"]),
    "{:.2e}".format(action["stop"]),
    "{:.2e}".format(action["step"]),
    action["bright_power"],
    spin_case,
    i_case,
    "{:.1e}".format(brightparams["nu0"]),
    "{:.1e}".format(brightparams["mass"]),
    float(brightparams["scale_height"]),
    "{:.3f}".format(brightparams["theta_b"]),
    "{:.2f}".format(float(brightparams["beta"])),
    "{:.1f}".format(float(brightparams["r_ie"])),
    "{:.1f}".format(float(brightparams["rb_0"])),
    "{:.1e}".format(brightparams["n_th0"]),
    "{:.1e}".format(brightparams["t_e0"]),
    float(brightparams["p_dens"]),
    float(brightparams["p_temp"]),
    "{:.1f}".format(brightparams["nscale"]),
    funckeys["emodelkey"],
    funckeys["bkey"],
    funckeys["nnoisykey"],
    funckeys["tnoisykey"],
    funckeys["bnoisykey"]
)

new_directory = iteration + '/'

final_vid_path = aart_results + movie_path + vid_path + new_directory

if os.path.isdir(final_vid_path):
    subprocess.run(['rm -r ' + final_vid_path], shell=True)
    subprocess.run(["mkdir " + final_vid_path], shell=True)
else:
    subprocess.run(["mkdir " + final_vid_path], shell=True)

# Thin Full movie-------------------------------------------------------------------------------------------------------
thin_full_movie_name = (final_vid_path + 'ThinFull_' + iteration + '.mp4')

subprocess.run(["ffmpeg -r " + str(speed) + " -i " + temp_plots_path +
                "Thin_Full_Fig_%d.png -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -vcodec "
                "libx264 -crf 10 -pix_fmt yuv420p " + thin_full_movie_name], shell=True)

# Thick Full movie------------------------------------------------------------------------------------------------------
thick_full_movie_name = final_vid_path + 'ThickFull_' + iteration + '.mp4'

subprocess.run(["ffmpeg -r " + str(
    speed) + " -i " + temp_plots_path + "Thick_Full_Fig_%d.png -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -vcodec "
                                        "libx264 -crf 10 -pix_fmt yuv420p " + thick_full_movie_name], shell=True)

# remove temporary files------------------------------------------------------------------------------------------------
for i in range(num_iterations):
    subprocess.run(['rm ' + thin_names[i] + ".png"], shell=True)
    subprocess.run(['rm ' + thick_names[i] + ".png"], shell=True)

