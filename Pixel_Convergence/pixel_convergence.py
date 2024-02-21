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
#
# import importlib
# #
# path_sub = "../aart_results/"
# newpath = "../aart_results/Pixel_Convergence/path_results"
# isExist = os.path.exists(newpath)
# if not isExist:
#     os.makedirs(newpath)
#     print("A directory (Results) was created to store the convergence results")

parser = argparse.ArgumentParser(description='Convergence test for pixel resolution', formatter_class=argparse.RawTextHelpFormatter)
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

mean_radii_Thin = np.ndarray([trials + 1, 4])  # [I0, I1, I2, cumu]
mean_radii_Thick = np.ndarray([trials + 1, 4])  # [I0, I1, I2, FullImage]

# Create intensity files
intensity_path = '/home/td6241/repositories/aart_convergence/aart_results/convergence_data/'+ iteration + '/'

# Create a directory for the results
isExist = os.path.exists(intensity_path)
if not isExist:
    os.makedirs(intensity_path)
    print("A directory was created to store intensity h.5 files")




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
        "emodelkey" : 0, # emodelkey Emission Model choice, 0 = thermal ultrarelativistic, 1 = power law
        "bkey" : 2, # bkey
        "nnoisykey" : 0, # nnoisykey Inoisy density. 0 = no noise, 1 = noise
        "tnoisykey" : 0, # tnoisykey Inoisy temperature
        "bnoisykey" : 0# bnoisykey Inoisy magnetic field
}
    


for i in range(trials+1):
    print("____________________________________________________________________________________________________")
    print("Iteration number: ", i)
    print("dx value: ", k)
    params.dx0 = k
    dx0 = params.dx0

    params.dx1 = dx0
    dx1 = params.dx1

    params.dx2 = dx0
    dx2 = params.dx2
    '''Computation of the lensing bands----------------------------------'''
    subprocess.run(['python3 ' + aartpath + '/lensingbands.py ' + str(dx0) + ' ' + str(dx1) + ' ' + str(dx2) ], shell=True)
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

    '''Analytical Ray-tracing----------------------------------'''
    subprocess.run(['python3 ' + aartpath + '/raytracing.py '], shell=True)
    fnrays1 = path + "Rays_a_%s_i_%s.h5" % (spin_case, i_case)

    print("Reading file: ", fnrays1)

    h5f = h5py.File(fnrays1, 'r')

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

#     a = spin_case
#     inc = i_case * np.pi / 180  # inclination angle
#     rh = 1 + np.sqrt(1 - a ** 2)  # event horizon
#     # angles to sample
#     varphis = np.linspace(-180, 179, 360) * np.pi / 180

#     # generate inner shadow (n=0) curve with kgeo
#     data_inner = kgeo.equatorial_lensing.rho_of_req(a, inc, rh, mbar=0, varphis=varphis)
#     (_, rhos_inner, alphas_inner, betas_inner) = data_inner

#     (area_inner, mux_inner, muy_inner, r_inner, e_inner, chi_inner) = curve_params(varphis, rhos_inner)
#     np.save('r_inner_spin_{}_inc_{}'.format(spin_case, i_case), r_inner)
#     np.save('alphas_inner_spin_{}_inc_{}'.format(spin_case, i_case), alphas_inner)
#     np.save('betas_inner_spin_{}_inc_{}'.format(spin_case, i_case), betas_inner)

#     # generate outer shadow (n=inf) curve with kgeo
#     data_outer = kgeo.equatorial_lensing.rho_of_req(a, inc, rh, mbar=5, varphis=varphis)
#     (_, rhos_outer, alphas_outer, betas_outer) = data_outer

#     (area_outer, mux_outer, muy_outer, r_outer, e_outer, chi_outer) = curve_params(varphis, rhos_outer)
#     np.save('r_outer_spin_{}_inc_{}'.format(spin_case, i_case), r_outer)
#     np.save('alphas_outer_spin_{}_inc_{}'.format(spin_case, i_case), alphas_outer)
#     np.save('betas_outer_spin_{}_inc_{}'.format(spin_case, i_case), betas_outer)
    
    
#     oldname = path+"LensingBands_a_%s_i_%s.h5"%(spin_case,i_case)
#     newname = oldname[:-3] + '_' + str(k) + '.h5'
#     subprocess.run(["mv " + oldname + " " + newname], shell=True)
    
#     oldname = path+"Rays_a_%s_i_%s.h5"%(spin_case,i_case)
#     newname = oldname[:-3] + '_' + str(k) + '.h5'
#     subprocess.run(["mv " + oldname + " " + newname], shell=True)

    args = ' '
    cmd1_args = {
        "nu0" : '--nu ',
        "mass" : '--mass ',
        "scale_height" : '--scaleh ',
        "theta_b" : '--thetab ',
        "beta" : '--beta ',
        "r_ie" : '--rie ',
        "rb_0" : '--rb0 ',
        "n_th0" : '--nth0 ',
        "t_e0" : '--te0 ',
        "p_dens" : '--pdens ',
        "p_temp" : '--ptemp ',
        "nscale" : '--nscale ',
    }

    cmd2_args = {
        "emodelkey" : '--emodelkey ',
        "bkey" : '--bkey ',
        "nnoisykey" : '--nnoisykey ',
        "tnoisykey" : '--tnoisykey ',
        "bnoisykey" : '--bnoisykey ',
    }


    for arg in cmd1_args:
        args = args + cmd1_args[arg] + str(brightparams[arg]) + ' ' 

    for arg in cmd2_args:
        args = args + cmd2_args[arg] + str(funckeys[arg]) + ' ' 
        
    subprocess.run(['python3 ' + aartpath + '/radialintensity.py' + args], shell=True)

        # TODO: modify name
    # fnrays='./Results/Intensity_a_{}|i_{}|nu_{}|mass_{}|scaleh_{}|thetab_{}|beta_{}|rie_{}|rb_{}|nth0_{}|te0_{}|pdens_{}|ptemp_{}|nscale_{}|emkey_{}|bkey_{}|nkey_{}|tnke
    fnrays= path + 'Intensity_a_{}_i_{}_nu_{}_mass_{}_scaleh_{}_thetab_{}_beta_{}_rie_{}_rb_{}_nth0_{}_te0_{}_pdens_{}_ptemp_{}_nscale_{}_emkey_{}_bkey_{}_nkey_{}_tnkey_{}_bnkey_{}.h5'.format(
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
    

    print("Reading file: ",fnrays)

    h5f = h5py.File(fnrays,'r')

    I0=h5f['bghts0'][:] # This implies I0 is 1 pass
    I1=h5f['bghts1'][:]
    I2=h5f['bghts2'][:]

    I2_Absorb = h5f['bghts2_absorbtion'][:]
    I1_Absorb = h5f['bghts1_absorbtion'][:]
    I0_Absorb = h5f['bghts0_absorbtion'][:]
    Absorbtion_Image =h5f['bghts_full_absorbtion'][:]

    h5f.close()
    
    oldname = fnrays
    newname = intensity_path + 'Intensity_' + iteration + '_iteration_' + str(int(i)) + '.h5'
    subprocess.run(["mv " + oldname + " " + newname], shell=True)
    print('file ' + newname + ' created')
    
    
    size = 100
    
    
    
    radii_I0_Thin_i, theta = tls.radii_of_theta(I0, size)
    radii_I1_Thin_i, theta = tls.radii_of_theta(I1, size)
    radii_I2_Thin_i, theta = tls.radii_of_theta(I2, size)
    radii_cumulative_Thin_i, theta = tls.radii_of_theta(I0 +I1 +I2, size)
    
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
    
    # remove uneeded files
    
    
    subprocess.run(["rm " + fnbands], shell=True)
    subprocess.run(["rm " + fnrays1], shell=True)

    k += action['step']

# k = action["start"]
# for i in range(trials+1):
    
#     k += action['step']
final_data_path = '/home/td6241/repositories/aart_convergence/aart_results/convergence_data/'


np.save(final_data_path + "mean_radii_Thin_" + iteration, mean_radii_Thin)
np.save(final_data_path + "mean_radii_Thick_" + iteration, mean_radii_Thick)

