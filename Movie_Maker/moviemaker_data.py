import os.path
import sys
import subprocess

import scipy

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
brightparams_label = dict(brightparams)
# brightparam_items = brightparams.items()
# brightparams_init = np.array(list(brightparam_items))[:,1]  # Convert object to a list

# File paths-------------------------
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
    "{:.1e}".format(brightparams_label["nu0"]),
    "{:.1e}".format(brightparams_label["mass"]),
    float(brightparams_label["scale_height"]),
    "{:.3f}".format(brightparams_label["theta_b"]),
    "{:.2f}".format(float(brightparams_label["beta"])),
    "{:.1f}".format(float(brightparams_label["r_ie"])),
    "{:.1f}".format(float(brightparams_label["rb_0"])),
    "{:.1e}".format(brightparams_label["n_th0"]),
    "{:.1e}".format(brightparams_label["t_e0"]),
    float(brightparams_label["p_dens"]),
    float(brightparams_label["p_temp"]),
    "{:.1f}".format(brightparams_label["nscale"]),
    funckeys["emodelkey"],
    funckeys["bkey"],
    funckeys["nnoisykey"],
    funckeys["tnoisykey"],
    funckeys["bnoisykey"]
)

movie_path = "Movie/"
vid_path = rez + "video/"
data_path = aart_results + movie_path + rez + "data_sets/"
new_directory = iteration + '/'

final_vid_path = aart_results + movie_path + vid_path + new_directory
final_data_path = final_vid_path + "data_sets/"

if os.path.isdir(final_vid_path):
    subprocess.run(['rm -r ' + final_vid_path], shell=True)
    subprocess.run(["mkdir " + final_vid_path], shell=True)
else:
    subprocess.run(["mkdir " + final_vid_path], shell=True)

subprocess.run(["mkdir " + final_data_path], shell=True)

# ----------------------------------------------------------------------------------------------------------------------
# size = 1000  # array size for radii calcs
# size = 500  # array size for radii calcs
size = 100  # array size for radii calcs
num_iterations = int((action["stop"] - action["start"]) / action["step"])
doth5_files = []  # .h5 files

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

# Intensity images
# I_thins = np.ndarray([num_iterations, 3])  # [I0, I1, I2]
# I_thicks = np.ndarray([num_iterations, 3])  # [I0, I1, I2]
#

hdot5_names = []

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

    I0_Absorb = h5f['bghts0_absorbtion'][:]
    I1_Absorb = h5f['bghts1_absorbtion'][:]
    I2_Absorb = h5f['bghts2_absorbtion'][:]
    Absorbtion_Image = h5f['bghts_full_absorbtion'][:]

    filename = final_data_path + "doth5" + str(i) + "_" + iteration

    h5f = h5py.File(filename, 'w')

    h5f.create_dataset('I0', data=I0)
    h5f.create_dataset('I1', data=I1)
    h5f.create_dataset('I2', data=I2)

    h5f.create_dataset('I0_Absorb', data=I0_Absorb)
    h5f.create_dataset('I1_Absorb', data=I1_Absorb)
    h5f.create_dataset('I2_Absorb', data=I2_Absorb)
    h5f.create_dataset('Absorbtion_Image', data=Absorbtion_Image)

    h5f.close()

    hdot5_names += [filename]
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
    janksys_thin[i, 0] = ilp.total_jy(I0, brightparams["nu0"], brightparams["mass"]).value
    janksys_thin[i, 1] = ilp.total_jy(I1, brightparams["nu0"], brightparams["mass"]).value
    janksys_thin[i, 2] = ilp.total_jy(I2, brightparams["nu0"], brightparams["mass"]).value
    janksys_thin[i, 3] = ilp.total_jy(I0 + I1 + I2, brightparams["nu0"], brightparams["mass"]).value

    janksys_thick[i, 0] = ilp.total_jy(I0_Absorb, brightparams["nu0"], brightparams["mass"]).value
    janksys_thick[i, 1] = ilp.total_jy(I1_Absorb, brightparams["nu0"], brightparams["mass"]).value
    janksys_thick[i, 2] = ilp.total_jy(I2_Absorb, brightparams["nu0"], brightparams["mass"]).value
    janksys_thick[i, 3] = ilp.total_jy(Absorbtion_Image, brightparams["nu0"], brightparams["mass"]).value

    h5f.close()

    subprocess.run(['rm ' + doth5_files[i]], shell=True)

hdot5_names = np.array(hdot5_names)
# Saving Data--------------------------------------------------------------------------------------------------------
np.save(final_data_path + "hdot5_names_" + iteration, hdot5_names)
np.save(final_data_path + "x_variable_" + iteration, x_variable)
np.save(final_data_path + "janksys_thick_" + iteration, janksys_thick)
np.save(final_data_path + "janksys_thin_" + iteration, janksys_thin)
np.save(final_data_path + "mean_radii_Thin_" + iteration, mean_radii_Thin)
np.save(final_data_path + "mean_radii_Thick_" + iteration, mean_radii_Thick)
np.save(final_data_path + "radii_I0_Thin_" + iteration, radii_I0_Thin)
np.save(final_data_path + "radii_I1_Thin_" + iteration, radii_I1_Thin)
np.save(final_data_path + "radii_I2_Thin_" + iteration, radii_I2_Thin)
np.save(final_data_path + "radii_FullAbsorption_Thick_" + iteration, radii_FullAbsorption_Thick)
np.save(final_data_path + "radii_I0_Thick_" + iteration, radii_I0_Thick)
np.save(final_data_path + "radii_I1_Thick_" + iteration, radii_I1_Thick)
np.save(final_data_path + "radii_I2_Thick_" + iteration, radii_I2_Thick)
np.save(final_data_path + "theta_" + iteration, theta)
