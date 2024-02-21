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

from moviemaker_params import (
    brightparams,
    funckeys,
    var_label,
    scale_label,
    units_label,
    aart_results,
    rez
                               )

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
parser.add_argument("absorb_power_key", type=int,
                    help="0 for no power normalization. 1 for power normalization")
parser.add_argument("absorb_power_val", type=float,
                    help="Colorbar power norm value in absorption. Lower for anticipated higher brightness range")
args = parser.parse_args()
action = {
    "var": args.var,
    "start": args.start,
    "stop": args.stop,
    "step": args.step_size,
    "bright_power": args.power_val,
    "absorb_power_key": args.absorb_power_key,
    "absorb_power": args.absorb_power_val
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

# ----------------------------------------------------------------------------------------------------------------------

# Data Reading_________________________________________---

iteration = (
    'Var_{}_Start_{}_Stop_{}_Step_{}_a_{}_i_{}_nu_{}'
    '_mass_{}_scaleh_{}_thetab_{}_beta_{}_rie_{}_rb_{}_nth0_{}_te0_{}'
    '_pdens_{}_ptemp_{}_nscale_{}_keys_{}_{}_{}_{}_{}'
).format(
    action["var"],
    "{:.2e}".format(action["start"]),
    "{:.2e}".format(action["stop"]),
    "{:.2e}".format(action["step"]),
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

# doth5_Var_nu0_Start_3.00e+10_Stop_6.00e+10_Step_1.00e+10_BP_0.4_a_0.94_i_17_nu_2.3e+11_mass_1.3e+43_scaleh_0.5_thetab_0.873_beta_1.00_rie_10.0_rb_2.0_nth0_2.0e+05_te0_2.0e+11_pdens_-0.7_ptemp_-0.84_nscale_0.4_keys_0_1_0_0_0
#
# File paths
temporary_path = aart_results + "Movie/" + rez + "temporary/"
final_vid_path = aart_results + "Movie/" + rez + "video/" + iteration + '/'
data_path = final_vid_path + "data_sets/"
final_graph_path = final_vid_path + "graph"

if os.path.isdir(final_vid_path):
    print("File for data of input found")
else:
    raise ValueError("File (" + final_vid_path + ") not found")

print("Loading Data...")
# Loading

hdot5_names = np.load(data_path + "hdot5_names_" + iteration + ".npy")
x_variable = np.load(data_path + "x_variable_" + iteration + ".npy")
janksys_thick = np.load(data_path + "janksys_thick_" + iteration + ".npy")
janksys_thin = np.load(data_path + "janksys_thin_" + iteration + ".npy")
mean_radii_Thin = np.load(data_path + "mean_radii_Thin_" + iteration + ".npy")
mean_radii_Thick = np.load(data_path + "mean_radii_Thick_" + iteration + ".npy")
radii_I0_Thin = np.load(data_path + "radii_I0_Thin_" + iteration + ".npy")
radii_I1_Thin = np.load(data_path + "radii_I1_Thin_" + iteration + ".npy")
radii_I2_Thin = np.load(data_path + "radii_I2_Thin_" + iteration + ".npy")
radii_FullAbsorption_Thick = np.load(data_path + "radii_FullAbsorption_Thick_" + iteration + ".npy")
radii_I0_Thick = np.load(data_path + "radii_I0_Thick_" + iteration + ".npy")
radii_I1_Thick = np.load(data_path + "radii_I1_Thick_" + iteration + ".npy")
radii_I2_Thick = np.load(data_path + "radii_I2_Thick_" + iteration + ".npy")
theta = np.load(data_path + "theta_" + iteration + ".npy")

# --------------------------
size = 1000  # array size for radii calcs
num_iterations = int((action["stop"] - action["start"]) / action["step"])

thin_names = []  # thin images to be deleted after movie is made
thick_names = []
thin_I2_names = []
thick_I2_names = []

# Black Hole Inner Shadow Calc--------------------------
r_inner = np.load('r_inner_spin_{}_inc_{}.npy'.format(spin_case, i_case))
alphas_inner = np.load('alphas_inner_spin_{}_inc_{}.npy'.format(spin_case, i_case))
betas_inner = np.load('betas_inner_spin_{}_inc_{}.npy'.format(spin_case, i_case))

# Black Hole Outer Shadow Calc--------------------------
r_outer = np.load('r_outer_spin_{}_inc_{}.npy'.format(spin_case, i_case))
alphas_outer = np.load('alphas_outer_spin_{}_inc_{}.npy'.format(spin_case, i_case))
betas_outer = np.load('betas_outer_spin_{}_inc_{}.npy'.format(spin_case, i_case))

# ------------------------------------------------------

dim = [16, 8]
xaxis = np.array(x_variable) / scale_label[action['var']]
one_M = ilp.rg_func(brightparams["mass"] * u.g).to(u.m)
M2uas = np.arctan(one_M.value / dBH) / muas_to_rad

for i in range(num_iterations):
    print('Creating Frame : ' + str(i))

    h5f = h5py.File(hdot5_names[i], 'r')

    I0 = h5f['I0'][:]  # This implies I0 is 1 pass
    I1 = h5f['I1'][:]
    I2 = h5f['I2'][:]



    I0_Absorb = h5f['I0_Absorb'][:]
    I1_Absorb = h5f['I1_Absorb'][:]
    I2_Absorb = h5f['I2_Absorb'][:]
    Absorbtion_Image = h5f['Absorbtion_Image'][:]

    h5f.close()
    # tau2 = h5f['tau2'][:]
    # tau1 = h5f['tau1'][:]
    # tau0 = h5f['tau0'][:]
    # full_profiles0 = h5f['full_profiles0'][:]
    # full_profiles1 = h5f['full_profiles1'][:]
    # full_profiles2 = h5f['full_profiles2'][:]
    # full_profiles_unit = h5f['full_profiles_unit'][:]

    max_factor = 1.2
    if i == 0:
        vmax0 = np.max(I0 + I1 + I2) * max_factor
        vmax1 = Absorbtion_Image.max() * max_factor
        vmax2 = np.max(I2) * max_factor * max_factor
        vmax3 = np.max(I2_Absorb) * max_factor

    if i != 0:
        for b in range(len(p_full_thin)):
            p_full_thin[b].remove()
            if b != 6:
                p_full_thick[b].remove()
            # p_I2_thin[i].remove()
            # p_I2_thick[i].remove()
    else:
        p_full_thin = [None, None, None, None, None, None, None]
        p_full_thick = [None, None, None, None, None, None, None, None]
        # p_I2_thin = [None, None, None, None, None, None, None]
        # p_I2_thick = [None, None, None, None, None, None, None]

    # Thin Full Image---------------------------------------------------------------------------------------------------

    # vmax0 = np.max(I0 + I1 + I2) * 1.2
    # fig, ax = plt.subplots(figsize=dim, dpi=400)

    # JANKSKY PLOTS________________________________
    fig = plt.subplots(2, 2, figsize=[16, 8], dpi=400, width_ratios=[1, 1.3])
    ax = [None, None, None]
    ax[0] = plt.subplot(2, 2, 1)
    ax[0].plot(xaxis, janksys_thin[:, 0], '-', label='n=0', color='tab:red', linewidth=3)
    ax[0].plot(xaxis, janksys_thin[:, 1], ':', label='n=1', color='tab:orange', linewidth=3)
    ax[0].plot(xaxis, janksys_thin[:, 2], '--', label='n=2', color='tab:blue', linewidth=3)
    ax[0].plot(xaxis, janksys_thin[:, 3], '-.', label='Total', color='tab:purple', linewidth=3)

    ax[0].axhline(.5, color='k', label=R'.5 $J_y$', linestyle=":")
    ax[0].axvline(230, color='k', linestyle=":")

    # Labels
    ax[0].set_ylabel("Total Flux ({})".format(R'$J_y$'))
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')

    ax[0].xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
    ax[0].xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    ax[0].yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
    ax[0].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))

    ax[0].tick_params('x', which="both", labelbottom=False)
    ax[0].tick_params('y', which="minor", labelleft=False)
    n = 4  # Keeps every 4th label
    [l.set_visible(False) for (i, l) in enumerate(ax[0].xaxis.get_minorticklabels()) if i % n != 0]
    ax[0].tick_params('both', length=10, width=1, which='major', labelrotation=80)
    ax[0].set_xlim(xaxis[0], xaxis[xaxis.size - 1])
    ax[0].legend(loc='lower left')
    # Markers
    p_full_thin[0] = ax[0].scatter(xaxis[i], janksys_thin[i, 3], color='tab:purple')
    p_full_thin[1] = ax[0].scatter(xaxis[i], janksys_thin[i, 0], color='tab:red')
    p_full_thin[2] = ax[0].scatter(xaxis[i], janksys_thin[i, 1], color='tab:orange')
    p_full_thin[3] = ax[0].scatter(xaxis[i], janksys_thin[i, 2], color='tab:blue')

    # RADII PLOTS___________________________________
    ax[1] = plt.subplot(2, 2, 3, sharex=ax[0])
    # ax[1].axhline(r_inner, color='k', linewidth=3, linestyle=":")  # label='Blackhole Inner Shadow'
    ax[1].axhline(r_outer, color='dimgrey', linewidth=5)  # label='Blackhole Outer Shadow'
    ax[1].plot(xaxis, mean_radii_Thin[:, 0], '-', label='n=0', color='tab:red', linewidth=3)
    ax[1].plot(xaxis, mean_radii_Thin[:, 1], ':', label='n=1', color='tab:orange', linewidth=3)
    ax[1].plot(xaxis, mean_radii_Thin[:, 2], '-.', label='n=2', color='tab:blue', linewidth=3)

    ax[1].axvline(230, color='k', linestyle=":")
    # Labels
    ax[1].set_xlabel(var_label[action["var"]].replace('=', '') + ' (' + units_label[action["var"]] + ')')
    ax[1].set_ylabel("Ring Radii ({})".format(R'$R_g$'))

    ax[1].set_xscale('log')
    ax[1].set_yscale('log')

    ax[1].xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.0f'))
    ax[1].xaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
    ax[1].yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
    ax[1].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))

    n = 4  # Keeps every 4th label
    [l.set_visible(False) for (i, l) in enumerate(ax[1].xaxis.get_minorticklabels()) if i % n != 0]

    ax[1].legend(frameon=False)

    ax[1].set_xlim(xaxis[0], xaxis[xaxis.size - 1])

    new_ticks = [xaxis[0], 230, xaxis[xaxis.size - 1]]
    ax[1].set_xticks(new_ticks)
    ax[1].tick_params('x', length=20, width=1, which='major')

    # Markers
    p_full_thin[4] = ax[1].scatter(xaxis[i], mean_radii_Thin[i, 0], color='tab:red')
    p_full_thin[5] = ax[1].scatter(xaxis[i], mean_radii_Thin[i, 1], color='tab:orange')
    p_full_thin[6] = ax[1].scatter(xaxis[i], mean_radii_Thin[i, 2], color='tab:blue')

    # MAIN IMAGE____________________________________
    alphaShadow = r_inner * np.cos(theta)
    betaShadow = r_inner * np.sin(theta)
    alpha0 = radii_I0_Thin[i + 1, :] * np.cos(theta)
    beta0 = radii_I0_Thin[i + 1, :] * np.sin(theta)
    alpha1 = radii_I1_Thin[i + 1, :] * np.cos(theta)
    beta1 = radii_I1_Thin[i + 1, :] * np.sin(theta)
    alpha2 = radii_I2_Thin[i + 1, :] * np.cos(theta)
    beta2 = radii_I2_Thin[i + 1, :] * np.sin(theta)

    ax[2] = plt.subplot(1, 2, 2)

    im = ax[2].imshow(I0 + I1 + I2, origin="lower", cmap="afmhot", extent=[-lim0, lim0, -lim0, lim0],
                      norm=matplotlib.colors.PowerNorm(action["bright_power"], vmax=vmax0))

    ax[2].set_xlim(-10, 10)  # units of M
    ax[2].set_ylim(-10, 10)

    ax[2].set_xlabel(r"$\alpha$" + " " + r"($\mu as$)")
    ax[2].set_ylabel(r"$\beta$" + " " + r"($\mu as$)")
    ax[2].text(-9, 8.5, var_label[action["var"]] + str(round(x_variable[i] / scale_label[action["var"]], 2))
               + ' ' + units_label[action["var"]], fontsize=12, color="w")

    ax[2].set_xlim(-10, 10)  # units of M
    ax[2].set_ylim(-10, 10)

    ax[2].set_xlabel(r"$\alpha$" + " " + r"($\mu as$)")
    ax[2].set_ylabel(r"$\beta$" + " " + r"($\mu as$)")

    ax[2].set_xticks([-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10], labels=[
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

    ax[2].set_yticks([-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10], labels=[
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

    # ax[2].plot(alphaShadow, betaShadow, color='white', linestyle='--')
    ax[2].plot(alpha0, beta0, color='tab:red', linestyle='-')
    ax[2].plot(alpha1, beta1, color='tab:orange', linestyle=':')
    ax[2].plot(alpha2, beta2, color='tab:blue', linestyle='--')
    colorbar0 = plt.colorbar(im, fraction=0.046, pad=0.04, format='%.1e', ticks=[
        vmax0 * .8,
        vmax0 * .6,
        vmax0 * .4,
        vmax0 * .2,
        vmax0 * .05
    ],
                             label="Brightnes Temperature (K)",
                             ax=ax[2]
                             )

    thin_names += [temporary_path + 'Thin_Full_Fig_' + str(i)]
    plt.savefig(
        thin_names[i],
        bbox_inches='tight'
    )
    plt.close()

    # Thick Full Image--------------------------------------------------------------------------------------------------

    conv_1 = action["start"] + action["step"] * ilp.ring_convergance(mean_radii_Thick[:, 2], mean_radii_Thick[:, 3], 5)
    conv_1 = conv_1 / scale_label[action['var']]

    # Janksy Plots------------------------------------
    fig = plt.subplots(2, 2, figsize=[16, 8], dpi=400, width_ratios=[1, 1.3])
    ax1 = [None, None, None]
    ax1[0] = plt.subplot(2, 2, 1)
    ax1[0].plot(xaxis, janksys_thick[:, 0], '-', label='One pass', color='tab:red', linewidth=3)
    ax1[0].plot(xaxis, janksys_thick[:, 1], ':', label='Two passes', color='tab:orange', linewidth=3)
    ax1[0].plot(xaxis, janksys_thick[:, 2], '--', label='Three passes', color='tab:blue', linewidth=3)
    ax1[0].plot(xaxis, janksys_thick[:, 3], '-.', label='Cumulative', color='tab:purple', linewidth=3)

    flux_peak = action["start"] + action["step"] * np.argmax(janksys_thick[:, 3])
    flux_peak = flux_peak / scale_label[action['var']]

    ax1[0].axhline(.5, color='k', label=R'.5 $J_y$', linestyle=":")
    ax1[0].axvline(230, color='k', linestyle=":")
    ax1[0].axvline(conv_1, color='k', linestyle="--", linewidth=3)
    ax1[0].axvline(flux_peak, color='dimgrey', linestyle="-", linewidth=2)

    # Labels
    ax1[0].set_ylabel("Total Flux ({})".format(R'$J_y$'))
    ax1[0].set_xscale('log')
    ax1[0].set_yscale('log')

    ax1[0].xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
    ax1[0].xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    ax1[0].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0e"))

    ax1[0].tick_params('x', which="both", labelbottom=False)
    ax1[0].tick_params('y', which="minor", labelleft=False)
    n = 4  # Keeps every 4th label
    [l.set_visible(False) for (i, l) in enumerate(ax1[0].xaxis.get_minorticklabels()) if i % n != 0]
    ax1[0].tick_params('both', length=10, width=1, which='major')
    ax1[0].set_xlim(xaxis[0], xaxis[xaxis.size - 1])
    ax1[0].legend(loc='lower left')

    # Markers
    p_full_thick[0] = ax1[0].scatter(xaxis[i], janksys_thick[i, 0], color='tab:red')
    p_full_thick[1] = ax1[0].scatter(xaxis[i], janksys_thick[i, 1], color='tab:orange')
    p_full_thick[2] = ax1[0].scatter(xaxis[i], janksys_thick[i, 2], color='tab:blue')
    p_full_thick[3] = ax1[0].scatter(xaxis[i], janksys_thick[i, 3], color='tab:purple')



    # RADII PLOTS___________________________________
    ax1[1] = plt.subplot(2, 2, 3, sharex=ax1[0])
    # ax1[1].axhline(r_inner, color='k', linewidth=3, linestyle=":")  # , label='Blackhole Inner Shadow'
    ax1[1].axhline(r_outer, color='dimgrey', linewidth=5)  # , label='Blackhole Outer Shadow'
    ax1[1].plot(xaxis, mean_radii_Thick[:, 0], '-', label='One pass', color='tab:red', linewidth=3)
    ax1[1].plot(xaxis, mean_radii_Thick[:, 1], ':', label='Two passes', color='tab:orange', linewidth=3)
    ax1[1].plot(xaxis, mean_radii_Thick[:, 2], '--', label='Three passes', color='tab:blue', linewidth=3)
    ax1[1].plot(xaxis, mean_radii_Thick[:, 3], '-.', label='Cumulative', color='tab:purple', linewidth=3)
    ax1[1].axvline(230, color='k', linestyle=":")
    ax1[1].axvline(conv_1, color='k', linestyle="--", linewidth=3)
    ax1[1].axvline(flux_peak, color='dimgrey', linestyle="-", linewidth=2)

    # Labels
    ax1[1].set_xlabel(var_label[action["var"]].replace('=', '') + ' (' + units_label[action["var"]] + ')')
    ax1[1].set_ylabel("Ring Radii ({})".format(R'$R_g$'))
    ax1[1].set_xscale('log')
    ax1[1].set_yscale('log')
    n = 4  # Keeps every 4th label
    [l.set_visible(False) for (i, l) in enumerate(ax1[1].xaxis.get_minorticklabels()) if i % n != 0]

    ax1[1].xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.0f'))
    ax1[1].xaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
    ax1[1].yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
    ax1[1].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))

    ax1[1].legend(frameon=False)
    ax1[1].set_xlim(xaxis[0], xaxis[xaxis.size - 1])
    new_ticks = [xaxis[0], 230, conv_1, flux_peak, xaxis[xaxis.size - 1]]
    ax1[1].set_xticks(new_ticks)
    ax1[1].tick_params('x', length=20, width=1, which='major', labelrotation=80)
    # Markers
    p_full_thick[4] = ax1[1].scatter(xaxis[i], mean_radii_Thick[i, 0], color='tab:red')
    p_full_thick[5] = ax1[1].scatter(xaxis[i], mean_radii_Thick[i, 1], color='tab:orange')
    p_full_thick[6] = ax1[1].scatter(xaxis[i], mean_radii_Thick[i, 2], color='tab:blue')
    p_full_thick[7] = ax1[1].scatter(xaxis[i], mean_radii_Thick[i, 3], color='tab:purple')

    # MAIN IMAGE____________________________________
    # TODO: Should there be other labels for I2 and I1?
    alphaShadow = r_inner * np.cos(theta)
    betaShadow = r_inner * np.sin(theta)
    alpha0 = radii_I0_Thick[i + 1, :] * np.cos(theta)
    beta0 = radii_I0_Thick[i + 1, :] * np.sin(theta)
    alpha1 = radii_I1_Thick[i + 1, :] * np.cos(theta)
    beta1 = radii_I1_Thick[i + 1, :] * np.sin(theta)
    alpha2 = radii_I2_Thick[i + 1, :] * np.cos(theta)
    beta2 = radii_I2_Thick[i + 1, :] * np.sin(theta)
    alpha3 = radii_FullAbsorption_Thick[i + 1, :] * np.cos(theta)
    beta3 = radii_FullAbsorption_Thick[i + 1, :] * np.sin(theta)

    ax1[2] = plt.subplot(1, 2, 2)

    if action["absorb_power_key"] == 1:
        im1 = ax1[2].imshow(Absorbtion_Image, origin="lower", cmap="afmhot", extent=[-lim0, lim0, -lim0, lim0],
                            norm=matplotlib.colors.PowerNorm(action["absorb_power"], vmax=vmax1))
    else:
        im1 = ax1[2].imshow(Absorbtion_Image, vmax=vmax1, origin="lower", cmap="afmhot",extent=[-lim0, lim0, -lim0, lim0])

    ax1[2].set_xlim(-10, 10)  # units of M
    ax1[2].set_ylim(-10, 10)

    ax1[2].set_xlabel(r"$\alpha$" + " " + r"($\mu as$)", fontsize=14)
    ax1[2].set_ylabel(r"$\beta$" + " " + r"($\mu as$)", fontsize=14)
    ax1[2].text(-9, 8.5, var_label[action["var"]] + str(round(x_variable[i] / scale_label[action["var"]], 2))
                + ' ' + units_label[action["var"]], fontsize=12, color="w")

    ax1[2].set_xlim(-10, 10)  # units of M
    ax1[2].set_ylim(-10, 10)

    ax1[2].set_xlabel(r"$\alpha$" + " " + r"($\mu as$)")
    ax1[2].set_ylabel(r"$\beta$" + " " + r"($\mu as$)")

    ax1[2].set_xticks([-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10], labels=[
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

    ax1[2].set_yticks([-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10], labels=[
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

    # ax1[2].plot(alphaShadow, betaShadow, color='white', linestyle='--')
    ax1[2].plot(alpha1, beta1, color='tab:orange', linestyle=':')
    ax1[2].plot(alpha0, beta0, color='tab:red', linestyle='-')
    ax1[2].plot(alpha2, beta2, color='tab:blue', linestyle='--')
    ax1[2].plot(alpha3, beta3, color='tab:purple', linestyle='-.')

    colorbar1 = plt.colorbar(im1, fraction=0.046, pad=0.04, format='%.1e', ticks=[
        vmax1 * .8,
        vmax1 * .6,
        vmax1 * .4,
        vmax1 * .2,
        vmax1 * .05
    ],
                             label="Brightnes Temperature (K)",
                             ax=ax1[2]
                             )

    thick_names += [temporary_path + 'Thick_Full_Fig_' + str(i)]
    plt.savefig(
        thick_names[i],
        bbox_inches='tight'
    )
    plt.close()

    # Thin I2 Image-----------------------------------------------------------------------------------------------------

    fig, ax2 = plt.subplots(figsize=dim, dpi=400)

    im2 = ax2.imshow(I2, origin="lower", cmap="afmhot", extent=[-lim0, lim0, -lim0, lim0],
                     norm=matplotlib.colors.PowerNorm(action["bright_power"], vmax=vmax2))

    ax2.set_xlim(-10, 10)  # units of M
    ax2.set_ylim(-10, 10)
    ax2.set_xlabel(r"$\alpha$" + " " + r"($\mu as$)")
    ax2.set_ylabel(r"$\beta$" + " " + r"($\mu as$)")

    colorbar2 = fig.colorbar(im2, fraction=0.046, pad=0.04, format='%.1e', ticks=[
        vmax2 * .8,
        vmax2 * .6,
        vmax2 * .4,
        vmax2 * .2,
        vmax2 * .05
    ],
                             label="Brightnes Temperature (K)",
                             ax=ax2
                             )

    ax2.text(-9, 8.5, var_label[action["var"]] + str(round(x_variable[i] / scale_label[action["var"]], 2))
             + ' ' + units_label[action["var"]], fontsize=12, color="w")

    thin_I2_names += [temporary_path + 'Thin_I2_Fig_' + str(i)]
    plt.savefig(
        thin_I2_names[i],
        bbox_inches='tight'
    )
    plt.close()

    # Thick I2 Image----------------------------------------------------------------------------------------------------

    fig, ax3 = plt.subplots(figsize=dim, dpi=400)
    im3 = ax3.imshow(I2_Absorb, vmax=vmax3, origin="lower", cmap="afmhot", extent=[-lim0, lim0, -lim0, lim0])

    ax3.set_xlim(-10, 10)  # units of M
    ax3.set_ylim(-10, 10)
    ax3.set_xlabel(r"$\alpha$" + " " + r"($\mu as$)")
    ax3.set_ylabel(r"$\beta$" + " " + r"($\mu as$)")

    colorbar3 = fig.colorbar(im3, fraction=0.046, pad=0.04, format='%.1e', ticks=[
        vmax3 * .8,
        vmax3 * .6,
        vmax3 * .4,
        vmax3 * .2,
        vmax3 * .05
    ],
                             label="Brightnes Temperature (K)",
                             ax=ax3
                             )

    ax3.text(-9, 8.5, var_label[action["var"]] + str(round(x_variable[i] / scale_label[action["var"]], 2))
             + ' ' + units_label[action["var"]], fontsize=12, color="w")

    thick_I2_names += [temporary_path + 'Thick_I2_Fig_' + str(i)]
    plt.savefig(
        thick_I2_names[i],
        bbox_inches='tight'
    )
    plt.close()

    # plt.subplots_adjust(wspace=.6)

    # '''Radii Calc'''

# FMPEG Stiching--------------------------------------------------------------------------------------------------------


# Thin Full movie-------------------------------------------------------------------------------------------------------
thin_full_movie_name = (final_vid_path + 'ThinFull_' + iteration + '.mp4')

if os.path.isfile(thin_full_movie_name):
    subprocess.run(['rm ' + thin_full_movie_name], shell=True)

subprocess.run(["ffmpeg -r " + str(speed) + " -i " + temporary_path +
                "Thin_Full_Fig_%d.png -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -vcodec "
                "libx264 -crf 10 -pix_fmt yuv420p " + thin_full_movie_name], shell=True)

# Thick Full movie------------------------------------------------------------------------------------------------------
thick_full_movie_name = final_vid_path + 'ThickFull_' + iteration + '.mp4'

if os.path.isfile(thick_full_movie_name):
    subprocess.run(['rm ' + thick_full_movie_name], shell=True)

subprocess.run(["ffmpeg -r " + str(
    speed) + " -i " + temporary_path + "Thick_Full_Fig_%d.png -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -vcodec "
                                       "libx264 -crf 10 -pix_fmt yuv420p " + thick_full_movie_name], shell=True)

# Thin I2 movie---------------------------------------------------------------------------------------------------------
thin_I2_movie_name = final_vid_path + 'ThinI2_' + iteration + '.mp4'

if os.path.isfile(thin_I2_movie_name):
    subprocess.run(['rm ' + thin_I2_movie_name], shell=True)

subprocess.run(["ffmpeg -r " + str(
    speed) + " -i " + temporary_path + "Thin_I2_Fig_%d.png -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -vcodec "
                                       "libx264 -crf 10 -pix_fmt yuv420p " + thin_I2_movie_name], shell=True)

# Thick I2 movie--------------------------------------------------------------------------------------------------------
thick_I2_movie_name = final_vid_path + 'ThickI2_' + iteration + '.mp4'

if os.path.isfile(thick_I2_movie_name):
    subprocess.run(['rm ' + thick_I2_movie_name], shell=True)

subprocess.run(["ffmpeg -r " + str(
    speed) + " -i " + temporary_path + "Thick_I2_Fig_%d.png -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' -vcodec "
                                       "libx264 -crf 10 -pix_fmt yuv420p " + thick_I2_movie_name], shell=True)

for i in range(num_iterations):
    subprocess.run(['rm ' + thin_names[i] + ".png"], shell=True)
    subprocess.run(['rm ' + thick_names[i] + ".png"], shell=True)
    subprocess.run(['rm ' + thick_I2_names[i] + ".png"], shell=True)
    subprocess.run(['rm ' + thin_I2_names[i] + ".png"], shell=True)
