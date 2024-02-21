import sys
import subprocess
aartpath = '/home/tej/Desktop/Code_Stuff/Repositories/aart' #insert path to aart repo
sys.path.append(aartpath)

import kgeo
from aart_func import *
from params import * # The file params.py contains all the relevant parameters for the simulations
from astropy import units as u

import importlib 
path_sub = "../aart_results/"

def curve_params(varphi, rho):
    """calculate Appendix B parameters for a curve rho(varphi)
       assume varphis are equally spaced!!!"""
          
    # spacing in varphi  
    dvarphi = varphi[-1]-varphi[-2]
    
    # area
    area = np.trapz(0.5*rho**2,dx=dvarphi)
    
    # centroid
    mux = np.trapz((rho**3*np.cos(varphi)) / (3*area), dx=dvarphi)
    muy = np.trapz((rho**3*np.sin(varphi)) / (3*area), dx=dvarphi)  

    # second moment
    Sxx = np.trapz((rho**4*np.cos(varphi)**2) / (4*area), dx=dvarphi) - mux**2
    Syy = np.trapz((rho**4*np.sin(varphi)**2) / (4*area), dx=dvarphi) - muy**2
    Sxy = np.trapz((rho**4*np.sin(varphi)*np.cos(varphi)) / (4*area), dx=dvarphi) - mux*muy
    
    # diagonalize 2nd moment matrix
    D = np.sqrt((Sxx-Syy)**2 + 4*Sxy*Sxy)
    a = np.sqrt(2*(Sxx + Syy + D))
    b = np.sqrt(2*(Sxx + Syy - D))
    
    #radius, eccentricity, position angle
    r = np.sqrt(0.5*(a**2 + b**2))
    e = np.sqrt(1-b**2/a**2)
    chi = 0.5*np.arcsin(2*Sxy/D)
    
    return (area, mux, muy, r, e, chi)    


'''Computation of the lensing bands----------------------------------'''
subprocess.run(['python3 ' + aartpath + '/lensingbands.py'], shell=True)
fnbands= path + "LensingBands_a_%s_i_%s.h5"%(spin_case,i_case)

print("Reading file: ",fnbands)

h5f = h5py.File(fnbands,'r')

#Points for the boundary of the BH shadow
alpha_critc=h5f['alpha'][:]
beta_critc=h5f['beta'][:]

#The concave hulls for the lensing bands
hull_0i=h5f['hull_0i'][:]
hull_0e=h5f['hull_0e'][:]
hull_1i=h5f['hull_1i'][:]
hull_1e=h5f['hull_1e'][:]
hull_2i=h5f['hull_2i'][:]
hull_2e=h5f['hull_2e'][:]

#The grid points for each lensing band
supergrid0=h5f['grid0'][:]
N0=int(h5f["N0"][0])
mask0=h5f['mask0'][:]
lim0=int(h5f["lim0"][0])
supergrid1=h5f['grid1'][:]
N1=int(h5f["N1"][0])
mask1=h5f['mask1'][:]
lim1=int(h5f["lim1"][0])
supergrid2=h5f['grid2'][:]
N2=int(h5f["N2"][0])
mask2=h5f['mask2'][:]
lim2=int(h5f["lim2"][0])

h5f.close()


'''Analytical Ray-tracing----------------------------------'''
subprocess.run(['python3 ' + aartpath + '/raytracing.py'], shell=True)
fnrays= path + "Rays_a_%s_i_%s.h5"%(spin_case,i_case)

print("Reading file: ",fnrays)

h5f = h5py.File(fnrays,'r')

rs0=h5f['rs0'][:]
sign0=h5f['sign0'][:]
t0=h5f['t0'][:]
phi0=h5f['phi0'][:]

rs1=h5f['rs1'][:]
sign1=h5f['sign1'][:]
t1=h5f['t1'][:]
phi1=h5f['phi1'][:]

rs2=h5f['rs2'][:]
sign2=h5f['sign2'][:]
t2=h5f['t2'][:]
phi2=h5f['phi2'][:]

h5f.close()




a = spin_case        
inc = i_case*np.pi/180 # inclination angle
rh = 1 + np.sqrt(1-a**2) # event horizon
# angles to sample
varphis = np.linspace(-180,179,360)*np.pi/180


# generate inner shadow (n=0) curve with kgeo
data_inner = kgeo.equatorial_lensing.rho_of_req(a,inc,rh,mbar=0,varphis=varphis)
(_, rhos_inner, alphas_inner, betas_inner) = data_inner

(area_inner, mux_inner, muy_inner, r_inner, e_inner, chi_inner) = curve_params(varphis,rhos_inner)
np.save('r_inner_spin_{}_inc_{}'.format(spin_case,  i_case), r_inner)
np.save('alphas_inner_spin_{}_inc_{}'.format(spin_case,  i_case), alphas_inner)
np.save('betas_inner_spin_{}_inc_{}'.format(spin_case,  i_case), betas_inner)

# generate outer shadow (n=inf) curve with kgeo
data_outer = kgeo.equatorial_lensing.rho_of_req(a,inc,rh,mbar=5,varphis=varphis)
(_, rhos_outer, alphas_outer, betas_outer) = data_outer

(area_outer, mux_outer, muy_outer, r_outer, e_outer, chi_outer) = curve_params(varphis,rhos_outer)
np.save('r_outer_spin_{}_inc_{}'.format(spin_case,  i_case), r_outer)
np.save('alphas_outer_spin_{}_inc_{}'.format(spin_case,  i_case), alphas_outer)
np.save('betas_outer_spin_{}_inc_{}'.format(spin_case,  i_case), betas_outer)

   