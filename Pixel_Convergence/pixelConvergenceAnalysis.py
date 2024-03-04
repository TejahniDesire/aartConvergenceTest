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
radii_images_path = images_path + 'radiiCalc/'
full_images_path = images_path + 'fullImage/'

# Create a directory for the results
isExist = os.path.exists(images_path)
if not isExist:
    os.makedirs(images_path)
    os.makedirs(radii_images_path)
    os.makedirs(full_images_path)
    print("A directory was created to all store images")


isExist = os.path.exists(full_images_path)
if not isExist:
    os.makedirs(full_images_path)
    print("A directory was created to store full images")

isExist = os.path.exists(radii_images_path)
if not isExist:
    os.makedirs(radii_images_path)
    print("A directory was created to store radii calc images")


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
ax1.axvline(1500, color="purple")

ax1.set_xlabel('Pixels')
ax1.set_ylabel('Mean radii for Optically Thick Model')
ax1.legend()
ax1.title.set_text('Optically Thin Assumption')

ax2.plot(x_var,mean_radii_thick[:,0],label="n_0", linewidth=3)
ax2.plot(x_var,mean_radii_thick[:,1],label="n_1", linewidth=3)
ax2.plot(x_var,mean_radii_thick[:,2],label="n_2", linewidth=3)
ax2.plot(x_var,mean_radii_thick[:,3],label="cumulative", linewidth=2)
ax2.axvline(1500, color="purple")

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

    

    '''Radii Calc______________________'''
    # Thin
    radius, theta = tls.radii_of_theta(I0)
    # theta = np.arange(2*np.pi + 1, step=(2*np.pi + 1) / 100)
    radius1, theta1 = tls.radii_of_theta(I1)
    radius2, theta2 = tls.radii_of_theta(I2)
    radius3, theta = tls.radii_of_theta(I2 + I1 + I0)

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
    radius0, theta = tls.radii_of_theta(I0_Absorb)
    radius1, theta = tls.radii_of_theta(I1_Absorb)
    radius2, theta = tls.radii_of_theta(I2_Absorb)
    radius3, theta = tls.radii_of_theta(Absorbtion_Image)
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
    imagename = full_images_path + 'FullImage_' + str(i) +  ".jpeg"
    plt.savefig(imagename ,bbox_inches='tight')
    print("Jpeg Created:  " + imagename)

    plt.close()
    
    # Radii calc graphs__________________________________________
    nu0 = 230e9
    mass = (MMkg * u.kg).to(u.g).value
    rsize = tls.rsize
    rmax = I0.shape[0] * .4

    peak012, interp012 = tls.radii_of_theta_data(I0 + I1 + I2)
    # peak0, interp0  = tls.radii_of_theta_data(I0)
    # peak1, interp1  = tls.radii_of_theta_data(I1)
    # peak2, interp2  = tls.radii_of_theta_data(I2)
    peakAbsorb, interpAbsorb = tls.radii_of_theta_data(Absorbtion_Image)

    peaks = [peak012, peakAbsorb]
    interps = [interp012, interpAbsorb]

    for L in range(len(interps)):
        interps[L] = ilp.to_bright_temp(interps[L], nu0, mass)

    fig, dum = plt.subplots(2, 2, figsize=[15, 7], dpi=400)
    ax0 = plt.subplot(2, 2, 1)
    ax1 = plt.subplot(2, 2, 2)
    ax2 = plt.subplot(2, 2, 3)
    ax3 = plt.subplot(2, 2, 4)

    axes_0 = [ax0, ax2]
    axes_1 = [ax1, ax3]

    images = [I0 + I1 + I2, Absorbtion_Image]
    model = ["for Thin Assumption", "for Full Solution"]
    for J in range(2):
        x = np.linspace(0, rmax - 1, rsize) * k
        # x = np.linspace(0,rsize-1, rsize) * dx_0
        ptheta = [0, np.pi / 2, np.pi]
        colors = ['tab:blue', 'tab:green', 'tab:red']
        parg = []
        for L in range(len(ptheta)):
            parg += [tls.rad_to_arg(ptheta[L])]
            axes_0[J].plot(x, interps[J][parg[L]], linewidth=2, color=colors[L],
                           label=R"$\theta= $" + f"{ptheta[L]:.2f}")
            axes_0[J].axvline(peaks[J][parg[L]], color=colors[L])

        axes_0[J].set_xlim([0, 10])
        axes_0[J].legend()
        axes_0[J].set_xlabel(R"$R_g$")
        axes_0[J].set_ylabel(R"Flux Value " + model[J])

        im1 = axes_1[J].imshow(images[J], origin="lower", cmap="afmhot", extent=[-lim0, lim0, -lim0, lim0])

        axes_1[J].set_xlim(-10, 10)  # units of M
        axes_1[J].set_ylim(-10, 10)

        axes_1[J].set_xlabel(r"$\alpha$" + " " + r"($\mu as$)")
        axes_1[J].set_ylabel(r"$\beta$" + " " + r"($\mu as$)")

        # Plot lines
        rline1 = np.array([0, 10])
        theta1 = np.array([0, 0])

        alpha1 = rline1 * np.cos(theta1)
        beta1 = rline1 * np.sin(theta1)

        for L in range(len(ptheta)):
            rline = np.array([0, 10])
            theta = np.array([ptheta[L], ptheta[L]])

            alpha = rline * np.cos(theta)
            beta = rline * np.sin(theta)

            axes_1[J].plot(alpha, beta, color=colors[L], linestyle='--')

        colorbar0 = fig.colorbar(im1, fraction=0.046, pad=0.04, format='%.1e', ticks=[
            vmax0 * .8,
            vmax0 * .6,
            vmax0 * .4,
            vmax0 * .2,
            vmax0 * .05
        ],
                                 label="Brightnes Temperature (K)",
                                 ax=axes_1[J]
                                 )

    imagename = radii_images_path + 'RadiiImage_' + str(i) + ".jpeg"
    plt.savefig(imagename,bbox_inches='tight')
    print("Jpeg Created:  " + imagename)

    plt.close()
        
    


    
