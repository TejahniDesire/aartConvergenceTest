import os.path
import sys
import subprocess

import matplotlib.pyplot as plt

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
args = parser.parse_args()
action = {
    "var": args.var,
    "start": args.start,
    "stop": args.stop,
    "step": args.step_size,
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

# File paths
temporary_path = aart_results + "Movie/" + rez + "temporary/"
final_vid_path = aart_results + "Movie/" + rez + "video/" + iteration + '/'
data_path = final_vid_path + "data_sets/"
final_graph_path = final_vid_path + "graph/"

if os.path.isdir(final_vid_path):
    print("File for data of input found")
else:
    raise ValueError("File (" + final_vid_path + ") not found")

if os.path.isdir(final_graph_path):
    subprocess.run(['rm -r ' + final_graph_path], shell=True)
    subprocess.run(["mkdir " + final_graph_path], shell=True)
else:
    subprocess.run(["mkdir " + final_graph_path], shell=True)

print("Loading Data...")
# Loading
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


#--------------------------
size = 1000  # array size for radii calcs
num_iterations = int((action["stop"] - action["start"]) / action["step"])


# Black Hole Inner Shadow Calc--------------------------
r_inner = np.load('r_inner_spin_{}_inc_{}.npy'.format(spin_case, i_case))
alphas_inner = np.load('alphas_inner_spin_{}_inc_{}.npy'.format(spin_case,  i_case))
betas_inner = np.load('betas_inner_spin_{}_inc_{}.npy'.format(spin_case,  i_case))

# Black Hole Outer Shadow Calc--------------------------
r_outer = np.load('r_outer_spin_{}_inc_{}.npy'.format(spin_case,  i_case))
alphas_outer = np.load('alphas_outer_spin_{}_inc_{}.npy'.format(spin_case,  i_case))
betas_outer = np.load('betas_outer_spin_{}_inc_{}.npy'.format(spin_case,  i_case))

# THE STATIC GRAPH------------------------------------------------------------------------------------------------------
dim = [10, 8]
xaxis = np.array(x_variable) / scale_label[action['var']]
one_M = ilp.rg_func(brightparams["mass"] * u.g).to(u.m)
M2uas = np.arctan(one_M.value / dBH) / muas_to_rad


# Thin Full Image---------------------------------------------------------------------------------------------------

# vmax0 = np.max(I0 + I1 + I2) * 1.2
# fig, ax = plt.subplots(figsize=dim, dpi=400)

# JANKSKY PLOTS________________________________
fig = plt.subplots(2, 1, figsize=dim, dpi=400)
ax = [None, None, None]
ax[0] = plt.subplot(2, 1, 1)
ax[0].plot(xaxis, janksys_thin[:, 0], '-', label='n=0', color='tab:red', linewidth=3)
ax[0].plot(xaxis, janksys_thin[:, 1], ':', label='n=1', color='tab:orange', linewidth=3)
ax[0].plot(xaxis, janksys_thin[:, 2], '--', label='n=2', color='tab:blue', linewidth=3)
ax[0].plot(xaxis, janksys_thin[:, 3], '-.', label='Total', color='tab:purple', linewidth=3)


flux_peak = action["start"] + action["step"] * np.argmax(janksys_thin[:, 3])
flux_peak = flux_peak / scale_label[action['var']]

ax[0].axhline(.5, color='k',label=R'.5 $J_y$', linestyle=":")
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
ax[0].tick_params('both', length=10, width=1, which='major')
ax[0].set_xlim(xaxis[0], xaxis[xaxis.size - 1])
ax[0].legend(loc='lower left')

# RADII PLOTS___________________________________
ax[1] = plt.subplot(2, 1, 2, sharex=ax[0])
# ax[1].axhline(r_inner, color='k', linewidth=3, linestyle=":")  # , label='Blackhole Inner Shadow'
ax[1].axhline(r_outer, color='dimgrey', linewidth=5)  # , label='Blackhole Outer Shadow'
ax[1].plot(xaxis, mean_radii_Thin[:, 0], '-', label='n=0', color='tab:red', linewidth=3)
ax[1].plot(xaxis, mean_radii_Thin[:, 1], ':', label='n=1', color='tab:orange', linewidth=3)
ax[1].plot(xaxis, mean_radii_Thin[:, 2], '-.', label='n=2', color='tab:blue', linewidth=3)

ax[1].axvline(flux_peak, color='k', linestyle="-.")
ax[1].axvline(230, color='k', linestyle=":")
# Labels

ax[1].set_xscale('log')
ax[1].set_yscale('log')

ax[1].set_xlabel(var_label[action["var"]].replace('=', '') + ' (' + units_label[action["var"]] + ')')
ax[1].set_ylabel("Ring Radii ({})".format(R'$R_g$'))

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
ax[1].tick_params('x', length=20, width=1, which='major', labelrotation=90)
plt.savefig(final_graph_path + "Thin_" + iteration + ".png")
# Markers
plt.close()

# Thick Full Image--------------------------------------------------------------------------------------------------

conv_1 = action["start"] + action["step"] * ilp.ring_convergance(mean_radii_Thick[:, 2], mean_radii_Thick[:, 3], 5)
conv_1 = conv_1 / scale_label[action['var']]

fig = plt.subplots(2, 1, figsize=dim, dpi=400)
ax1 = [None, None]
ax1[0] = plt.subplot(2, 1, 1)
ax1[0].plot(xaxis, janksys_thick[:, 0], '-', label='One pass', color='tab:red', linewidth=3)
ax1[0].plot(xaxis, janksys_thick[:, 1], ':', label='Two passes', color='tab:orange', linewidth=3)
ax1[0].plot(xaxis, janksys_thick[:, 2], '--', label='Three passes', color='tab:blue', linewidth=3)
ax1[0].plot(xaxis, janksys_thick[:, 3], '-.', label='Cumulative', color='tab:purple', linewidth=3)

flux_peak = action["start"] + action["step"] * np.argmax(janksys_thick[:, 3])
flux_peak = flux_peak / scale_label[action['var']]

ax1[0].axhline(.5, color='k',label=R'.5 $J_y$', linestyle=":")
ax1[0].axvline(230, color='k', linestyle=":")
ax1[0].axvline(conv_1, color='k', linestyle="--", linewidth=3)
ax1[0].axvline(flux_peak, color='dimgrey', linestyle="-", linewidth=2)



# Labels
ax1[0].set_ylabel("Total Flux ({})".format(R'$J_y$'))
ax1[0].set_xscale('log')
ax1[0].set_yscale('log')

ax1[0].xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.1f'))
ax1[0].xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
# ax1[0].yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.4f'))
ax1[0].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0e"))


ax1[0].tick_params('x', which="both", labelbottom=False)
ax1[0].tick_params('y', which="minor", labelleft=False)

n = 4  # Keeps every 4th label
[l.set_visible(False) for (i, l) in enumerate(ax1[0].xaxis.get_minorticklabels()) if i % n != 0]
ax1[0].tick_params('both', length=10, width=1, which='major')
ax1[0].set_xlim(xaxis[0], xaxis[xaxis.size - 1])

ax1[0].legend(loc='lower left')

# RADII PLOTS___________________________________
ax1[1] = plt.subplot(2, 1, 2, sharex=ax1[0])


# ax1[1].axhline(r_inner, color='k', linewidth=2, linestyle=":")  # , label='Blackhole Inner Shadow'
ax1[1].axhline(r_outer, color='dimgrey', linewidth=5)  #, label='Blackhole Outer Shadow'
# ax1[1].scatter(xaxis[0],r_outer, marker="o", linewidth=10)
# ax1[1].scatter(xaxis[0],r_inner, marker="o", linewidth=10)

ax1[1].plot(xaxis, mean_radii_Thick[:, 0], '-', label='One pass', color='tab:red', linewidth=3)
ax1[1].plot(xaxis, mean_radii_Thick[:, 1], ':', label='Two passes', color='tab:orange', linewidth=3)
ax1[1].plot(xaxis, mean_radii_Thick[:, 2], '--', label='Three passes', color='tab:blue', linewidth=3)
ax1[1].plot(xaxis, mean_radii_Thick[:, 3], '-.', label='Cumulative', color='tab:purple', linewidth=3)

ax1[1].axvline(230, color='k', linestyle=":")
ax1[1].axvline(conv_1, color='k', linestyle="--", linewidth=3)
ax1[1].axvline(flux_peak, color='dimgrey', linestyle="-", linewidth=2)

# ax1[1].axvline(conv_1, color='b',label=R'1% diff', linestyle="--")
# ax1[1].axvline(conv_5, color='g',label=R'5% diff', linestyle=":")
# ax1[1].axvline(conv_10, color='r',label=R'10% diff', linestyle="-.")


# Labels
ax1[1].set_xlabel(var_label[action["var"]].replace('=', '') + ' (' + units_label[action["var"]] + ')')
ax1[1].set_ylabel("Ring Radii ({})".format(R'$R_g$'))
ax1[1].set_xscale('log')
ax1[1].set_yscale('log')

ax1[1].xaxis.set_minor_formatter(ticker.FormatStrFormatter('%.0f'))
ax1[1].xaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
ax1[1].yaxis.set_minor_formatter(ticker.FormatStrFormatter('%.0f'))
ax1[1].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))

new_ticks = [xaxis[0], 230, conv_1, flux_peak, xaxis[xaxis.size - 1]]
ax1[1].set_xticks(new_ticks)

# new_ticks = np.append(ax1[1].get_yticks(), r_outer)
# new_ticks = np.append(new_ticks, r_inner)
# print(new_ticks)
# ax1[1].set_yticks(new_ticks)
n = 4  # Keeps every 4th label
[l.set_visible(False) for (i, l) in enumerate(ax1[1].xaxis.get_minorticklabels()) if i % n != 0]
ax1[1].legend(frameon=False)
ax1[1].set_xlim(xaxis[0], xaxis[xaxis.size - 1])


ax1[1].tick_params('x', length=20, width=1, which='major', labelrotation=80)

plt.savefig(final_graph_path + "Thick_" + iteration + ".png")
plt.close()
