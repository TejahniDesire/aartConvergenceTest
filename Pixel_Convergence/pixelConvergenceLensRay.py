import sys
import subprocess

import os.path
import kgeo

import PCPaths

aartpath = PCPaths.aartPath  # insert path to aart repo
sys.path.append(aartpath)
import image_tools as tls
from aart_func import *
from params import *  # The file params.py contains all the relevant parameters for the simulations
import params
from astropy import units as u
import numpy as np


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

# iteration_path = '/home/td6241/repositories/aart_convergence/aart_results/convergence_data/' + iteration + '/'
iteration_path = PCPaths.aartPathResults + iteration + '/'
lband_path = iteration_path + 'lbands/'  # lensing bands
rtray_path = iteration_path + 'rbands/'  # raytracing bands


# Create a directory for the results
isExist = os.path.exists(iteration_path)
if not isExist:
    os.makedirs(iteration_path)
    print("A directory was created to store intensity h.5 files")

isExist = os.path.exists(lband_path)
if not isExist:
    os.makedirs(lband_path)
    print("A directory was created to store lensing band h.5 files")

isExist = os.path.exists(rtray_path)
if not isExist:
    os.makedirs(rtray_path)
    print("A directory was created to store raytracing h.5 files")


k = action["start"]


for i in range(trials + 1):
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
    subprocess.run(['python3 ' + aartpath + '/lensingbands.py ' + str(dx0) + ' ' + str(dx1) + ' ' + str(dx2)],
                   shell=True)
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

    # Move lensing bands and ratracing bands
    subprocess.run(["mv " + fnbands + ' ' + lband_path + 'lensingband_' + str(i) + '.h5'], shell=True)
    subprocess.run(["mv " + fnrays1 + ' ' + rtray_path + 'raytracing_' + str(i) + '.h5'], shell=True)

    k += action['step']





