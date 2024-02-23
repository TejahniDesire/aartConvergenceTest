import sys
import subprocess
import os.path
import kgeo

aartPath = '/home/td6241/repositories/aart_convergence/aart_tdwork'  # insert path to aart repo
sys.path.append(aartPath)
import PCPaths
import image_tools as tls
from aart_func import *
from params import *  # The file params.py contains all the relevant parameters for the simulations
import params
from astropy import units as u

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

iteration = str(action["start"]) + '_' + str(action["stop"]) + '_' + str(action["step"])

iteration_path = PCPaths.aartPathResults + iteration + '/'
intent_path = iteration_path + 'intent/'
radial_data_path = iteration_path + 'radii/'
images_path = iteration_path + 'images/'

# Create a directory for the results
isExist = os.path.exists(images_path)
if not isExist:
    os.makedirs(images_path)
    print("A directory was created to store images")

mean_radii_thick = np.load(radial_data_path + "mean_radii_Thick_" + iteration + ".npy")
mean_radii_thin = np.load(radial_data_path + "mean_radii_Thin_" + iteration + ".npy")

trials = int((action["stop"] - action["start"]) / action["step"])

x_var = []
k = action["start"]
for i in range(trials+1):
    x_var += [k]
    k += action["step"]
    
x_var = 2 * limits / np.array(x_var)


fig, (ax1,ax2) = plt.subplots(1,2,figsize=[10,5])
fig.tight_layout(pad=5.0)

ax1.plot(x_var, mean_radii_thin[:,0],label="n_0", linewidth=3)
ax1.plot(x_var, mean_radii_thin[:,1],label="n_1", linewidth=3)
ax1.plot(x_var, mean_radii_thin[:,2],label="n_2", linewidth=3)
ax1.plot(x_var, mean_radii_thin[:,3],label="cumulative", linewidth=2)
ax1.axvline(555.555, color="purple")

ax1.set_xlabel('Pixels')
ax1.set_ylabel('Mean radii for Optically Thick Model')
ax1.legend()
ax1.title.set_text('Optically Thin Assumption')


# ax2.plot(x_var,mean_radii_thick[:,0],label="One Pass")
# ax2.plot(x_var,mean_radii_thick[:,1],label="Two Pass")
# ax2.plot(x_var,mean_radii_thick[:,2],label="Three Pass")
ax2.plot(x_var,mean_radii_thick[:,0],label="n_0", linewidth=3)
ax2.plot(x_var,mean_radii_thick[:,1],label="n_1", linewidth=3)
ax2.plot(x_var,mean_radii_thick[:,2],label="n_2", linewidth=3)
ax2.plot(x_var,mean_radii_thick[:,3],label="cumulative", linewidth=2)
ax2.axvline(555.555, color="purple")


ax2.set_xlabel('Pixels')
ax2.set_ylabel('Mean radii for Optically Thick Model')
ax2.legend()
ax2.title.set_text('Full Solution')


plt.savefig(images_path + 'conv_' + iteration + ".jpeg",bbox_inches='tight')

plt.close()



k = action["start"]


for i in range(trials+1):

    fnrays = intent_path + 'Intensity_innerIteration_' + str(int(i)) + '.h5'
    lim0 = 25
        
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
    2 * limits

    vmax0 = np.nanmax(I0+I1+I2)*1.2
    fig, (ax0, ax1) = plt.subplots(1,2,figsize=[15,7],dpi=400)
    
    # Optically Thin
    
    

    im0 = ax0.imshow(I0 + I1 + I2,vmax=vmax0, origin="lower",cmap="afmhot",extent=[-lim0,lim0,-lim0,lim0])

    ax0.set_xlim(-10,10) # units of M
    ax0.set_ylim(-10,10) 


    ax0.set_xlabel(r"$\alpha$"+" "+r"($\mu as$)")
    ax0.set_ylabel(r"$\beta$"+" "+r"($\mu as$)")
    ax0.title.set_text('Optically Thin Assumption')

    # Optically thick

    im1 = ax1.imshow(Absorbtion_Image, origin="lower",cmap="afmhot",extent=[-lim0,lim0,-lim0,lim0])

    # 
    ax1.set_xlim(-10,10) # units of M
    ax1.set_ylim(-10,10) 


    ax1.set_xlabel(r"$\alpha$"+" "+r"($\mu as$)")
    ax1.set_ylabel(r"$\beta$"+" "+r"($\mu as$)")
    
    ax1.text(-9, 8.5, 'dx= ' + str(k), fontsize=12, color="w")
    ax1.text(-9, 7.5, 'Pixels= ' + str(2 * limits/k), fontsize=12, color="w")
    ax1.title.set_text('Full Solution')
    
    colorbar0=fig.colorbar(im1, fraction=0.046, pad=0.04, format='%.1e', ticks=[
    vmax0*.8,
    vmax0*.6,
    vmax0*.4,
    vmax0*.2,
    vmax0*.05
    ],
    label="Brightnes Temperature (K)",
    ax=ax1
    )

    
    size = 100
    '''Radii Calc______________________'''
    # Thin
    radius, theta = tls.radii_of_theta(I0,size)
    # theta = np.arange(2*np.pi + 1, step=(2*np.pi + 1) / 100)
    radius1, theta1 = tls.radii_of_theta(I1,size)
    radius2, theta2 = tls.radii_of_theta(I2,size)
    radius3, theta = tls.radii_of_theta(I2 + I1 + I0,size)

    alpha0 =  radius * np.cos(theta)
    beta0 =  radius * np.sin(theta)
    alpha1 =  radius1 * np.cos(theta)
    beta1 =  radius1 * np.sin(theta)
    alpha2 =  radius2 * np.cos(theta)
    beta2 =  radius2 * np.sin(theta)
    alpha3 =  radius3 * np.cos(theta)
    beta3 =  radius3 * np.sin(theta)
    
    ax0.plot(alpha0, beta0, color='tab:blue', linestyle='-')
    ax0.plot(alpha1, beta1, color='tab:orange', linestyle=':')
    ax0.plot(alpha2, beta2, color='tab:green', linestyle='--')
    ax0.plot(alpha3, beta3, color='tab:red', linestyle='--')
    

    # Thick
    radius0, theta = tls.radii_of_theta(I0_Absorb, size)
    radius1, theta = tls.radii_of_theta(I1_Absorb, size)
    radius2, theta = tls.radii_of_theta(I2_Absorb, size)
    radius3, theta = tls.radii_of_theta(Absorbtion_Image,size)
    # theta = np.arange(2*np.pi + 1, step=(2*np.pi + 1) / 100)
    alpha0 =  radius0 * np.cos(theta)
    beta0 =  radius0 * np.sin(theta)
    alpha1 =  radius1 * np.cos(theta)
    beta1 =  radius1 * np.sin(theta)
    alpha2 =  radius2 * np.cos(theta)
    beta2 =  radius2 * np.sin(theta)
    alpha3 =  radius3 * np.cos(theta)
    beta3 =  radius3 * np.sin(theta)

    ax1.plot(alpha0, beta0, color='tab:blue', linestyle='-')
    ax1.plot(alpha1, beta1, color='tab:orange', linestyle=':')
    ax1.plot(alpha2, beta2, color='tab:green', linestyle='--')
    ax1.plot(alpha3, beta3, color='tab:red', linestyle='--')
    plt.subplots_adjust(wspace=.3)
    
    k += action['step']
    imagename = images_path + 'FullImage_' + str(i) +  ".jpeg"
    plt.savefig(imagename,bbox_inches='tight')
    print("Jpeg Created:  " + imagename)


    
