import sys
import subprocess
import os.path
import kgeo

aartpath = '/home/td6241/repositories/aart_convergence/aart_tdwork'  # insert path to aart repo
sys.path.append(aartpath)
import image_tools as tls
from aart_func import *
from params import *  # The file params.py contains all the relevant parameters for the simulations
import params
from astropy import units as u
import numpy as np

def curve_params(varphi, rho):
    """calculate Appendix B parameters for a curve rho(varphi)
       assume varphis are equally spaced!!!"""

    # spacing in varphi
    dvarphi = varphi[-1] - varphi[-2]

    # area
    area = np.trapz(0.5 * rho ** 2, dx=dvarphi)

    # centroid
    mux = np.trapz((rho ** 3 * np.cos(varphi)) / (3 * area), dx=dvarphi)
    muy = np.trapz((rho ** 3 * np.sin(varphi)) / (3 * area), dx=dvarphi)

    # second moment
    Sxx = np.trapz((rho ** 4 * np.cos(varphi) ** 2) / (4 * area), dx=dvarphi) - mux ** 2
    Syy = np.trapz((rho ** 4 * np.sin(varphi) ** 2) / (4 * area), dx=dvarphi) - muy ** 2
    Sxy = np.trapz((rho ** 4 * np.sin(varphi) * np.cos(varphi)) / (4 * area), dx=dvarphi) - mux * muy

    # diagonalize 2nd moment matrix
    D = np.sqrt((Sxx - Syy) ** 2 + 4 * Sxy * Sxy)
    a = np.sqrt(2 * (Sxx + Syy + D))
    b = np.sqrt(2 * (Sxx + Syy - D))

    # radius, eccentricity, position angle
    r = np.sqrt(0.5 * (a ** 2 + b ** 2))
    e = np.sqrt(1 - b ** 2 / a ** 2)
    chi = 0.5 * np.arcsin(2 * Sxy / D)

    return (area, mux, muy, r, e, chi)


parser = argparse.ArgumentParser(description='Convergence test for pixel resolution',
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('start', type=float)
parser.add_argument('stop', help='Inclusive Stop', type=float)
parser.add_argument("step_size", type=float)
args = parser.parse_args()
action = {
    "start": args.start,
    "stop": args.stop,
    "step": args.step_size,
}

trials = int((action["stop"] - action["start"]) / action["step"])

iteration = str(action["start"]) + '_' + str(action["stop"]) + '_' + str(action["step"])

# Create intensity files
iteration_path = '/home/td6241/repositories/aart_convergence/aart_results/convergence_data/' + iteration + '/'
lband_path = iteration_path + 'lbands/'  # lensing bands
rtray_path = iteration_path + 'rbands/'  # raytracing bands

# Create a directory for the results
isExist = os.path.exists(iteration_path)
if not isExist:
    os.makedirs(iteration_path)
    os.makedirs(lband_path)
    os.makedirs(rtray_path)
    print("A directory was created to store intensity h.5 files")

mean_radii_Thin = np.ndarray([trials + 1, 4])  # [I0, I1, I2, cumu]
mean_radii_Thick = np.ndarray([trials + 1, 4])  # [I0, I1, I2, FullImage]


k = action["start"]



brightparams = {
    "nu0": 230e9,  # 0 nu0
    "mass": (MMkg * u.kg).to(u.g).value,  # 1 mass
    "scale_height": .5,  # 2 scale_height
    "theta_b": 50.0 * (np.pi / 180),  # 3 theta_b
    "beta": 1.0,  # 4 beta
    "r_ie": 10.0,  # 5 rie
    "rb_0": 2,  # 7 rb_0
    "n_th0": 1.9e4,  # ,4e5 # 8 n_th0
    "t_e0": 2e11,  # 9 t_e0 1e12
    "p_dens": -.7,  # 10 p_dens
    "p_temp": -.84,  # 11 p_temp
    "nscale": .4  # Scale of Inoisy
}

funckeys = {
    "emodelkey": 0,  # emodelkey Emission Model choice, 0 = thermal ultrarelativistic, 1 = power law
    "bkey": 2,  # bkey
    "nnoisykey": 0,  # nnoisykey Inoisy density. 0 = no noise, 1 = noise
    "tnoisykey": 0,  # tnoisykey Inoisy temperature
    "bnoisykey": 0  # bnoisykey Inoisy magnetic field
}

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

for arg in cmd1_args:
    args = args + cmd1_args[arg] + str(brightparams[arg]) + ' '

for arg in cmd2_args:
    args = args + cmd2_args[arg] + str(funckeys[arg]) + ' '


for i in range(trials + 1):
    print("____________________________________________________________________________________________________")
    print("Iteration number: ", i)
    print("dx value: ", k)

    subprocess.run(['python3 ' + aartpath + '/radialintensity.py' + args], shell=True)

    # TODO: modify name
    # fnrays='./Results/Intensity_a_{}|i_{}|nu_{}|mass_{}|scaleh_{}|thetab_{}|beta_{}|rie_{}|rb_{}|nth0_{}|te0_{}|pdens_{}|ptemp_{}|nscale_{}|emkey_{}|bkey_{}|nkey_{}|tnke
    fnrays = path + 'Intensity_a_{}_i_{}_nu_{}_mass_{}_scaleh_{}_thetab_{}_beta_{}_rie_{}_rb_{}_nth0_{}_te0_{}_pdens_{}_ptemp_{}_nscale_{}_emkey_{}_bkey_{}_nkey_{}_tnkey_{}_bnkey_{}.h5'.format(
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

    # move intensity files

    print("Reading file: ", fnrays)

    h5f = h5py.File(fnrays, 'r')

    I0 = h5f['bghts0'][:]  # This implies I0 is 1 pass
    I1 = h5f['bghts1'][:]
    I2 = h5f['bghts2'][:]

    I2_Absorb = h5f['bghts2_absorbtion'][:]
    I1_Absorb = h5f['bghts1_absorbtion'][:]
    I0_Absorb = h5f['bghts0_absorbtion'][:]
    Absorbtion_Image = h5f['bghts_full_absorbtion'][:]

    h5f.close()

    oldname = fnrays
    newname = iteration_path + 'Intensity_' + iteration + '_iteration_' + str(int(i)) + '.h5'
    subprocess.run(["mv " + oldname + " " + newname], shell=True)
    print('file ' + newname + ' created')

    size = 100

    radii_I0_Thin_i, theta = tls.radii_of_theta(I0, size)
    radii_I1_Thin_i, theta = tls.radii_of_theta(I1, size)
    radii_I2_Thin_i, theta = tls.radii_of_theta(I2, size)
    radii_cumulative_Thin_i, theta = tls.radii_of_theta(I0 + I1 + I2, size)

    r0_thin = tls.curve_params(theta, radii_I0_Thin_i)
    r1_thin = tls.curve_params(theta, radii_I1_Thin_i)
    r2_thin = tls.curve_params(theta, radii_I2_Thin_i)
    cumu_thin = tls.curve_params(theta, radii_cumulative_Thin_i)

    mean_radii_Thin[i, 0] = r0_thin
    mean_radii_Thin[i, 1] = r1_thin
    mean_radii_Thin[i, 2] = r2_thin
    mean_radii_Thin[i, 3] = cumu_thin

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

    k += action['step']