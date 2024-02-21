import matplotlib.pyplot as plt
from astropy import units as u
from aart_func import *
from params import * 
from astropy import constants as const
from astropy import units as u
from lmfit import Parameters,minimize, fit_report

'''Physical'''
tls_nu0 = 230e9 * ilp.Hz
tls_mass = (MMkg * u.kg).to(u.g)  # m87 mass as put in AART
tls_scale_height = .5
tls_theta_b = 60 * (np.pi / 180) * ilp.rads
tls_beta = 1
tls_R_ie = 10
tls_redshift = 1

'''Math Coeff'''
tls_rb_0 = 2
tls_n_th0 = 1.0e5 * ilp.cmcubed
tls_t_e0 = 1.2e11 * ilp.kelv

#tls_n_th0 = 1.23e6 * cmcubed


'''Exponents'''
tls_p_temp = -.84
tls_p_dens = -.7

tls_Bchoice = 0


def plot_with_units(x_values, y_values, label=None, xlabel='' ,ylabel=''):
    plt.plot(x_values, y_values,  label=label)
    if str(x_values[0].unit) != u.dimensionless_unscaled:
        plt.xlabel(str(xlabel)+ " (" + str(x_values[0].unit) + ")")
    else:
        plt.xlabel(str(xlabel))
        
    if str(y_values[0].unit) != u.dimensionless_unscaled:
        plt.ylabel(str(ylabel)+ " (" + str(y_values[0].unit) + ")")
    else:
        plt.ylabel(str(ylabel))

def scatter_with_units(x_values, y_values, label=None, xlabel='' ,ylabel=''):
    plt.scatter(x_values, y_values,  label=label)
    if str(x_values[0].unit) != u.dimensionless_unscaled:
        plt.xlabel(str(xlabel)+ " (" + str(x_values[0].unit) + ")")
    else:
        plt.xlabel(str(xlabel))
        
    if str(y_values[0].unit) != u.dimensionless_unscaled:
        plt.ylabel(str(ylabel)+ " (" + str(y_values[0].unit) + ")")
    else:
        plt.ylabel(str(ylabel))


def center_finder(I0):
 
    horizontalline = np.zeros(2)
    verticalline = np.zeros(2)

    length = I0.shape[0]
    midpoint = int(length / 2)
    
    horizontalline[0] = I0[midpoint,0:midpoint].argmax()
    horizontalline[1] = midpoint + I0[midpoint,midpoint:length].argmax()

    verticalline[0] = I0[0:midpoint,midpoint].argmax()
    verticalline[1] = midpoint + I0[midpoint:length,midpoint].argmax()

    y = (verticalline[1] - verticalline[0]) / 2 + verticalline[0]
    x = (horizontalline[1] - horizontalline[0]) / 2 + horizontalline[0]
    return x, y

# line is 100 x 2000
def radii_finder(line):
    peaks = np.zeros(2)

    length = line.shape[0]
    midpoint = int(length / 2)
    
    peaks[0] = line[0:midpoint].argmax()
    peaks[1] = midpoint + line[midpoint:length].argmax()
    return (peaks[1] - peaks[0]) / 2


def radii_of_theta(I0, thetapointsamount):
    x = np.arange(I0.shape[0])
    y = x
    interp = RegularGridInterpolator((x,y), I0.T)
    rmax = 1000

    theta = np.matrix(np.arange(0, 2*np.pi + 1, (2*np.pi + 1) / thetapointsamount)) # 1 x thetapointsamount
    r = np.matrix(np.arange(rmax + 1)).T # (rmax + 1) x 1


    onest = np.matrix(np.ones(r.shape[0])).T # (rmax+1) x 1
    onesr = np.matrix(np.ones(thetapointsamount)) # 1 x thetapointsamount



    thetarray= onest @ theta # (2rmax+1) x thetapointsamount
    rarray = r @ onesr # (2rmax+1) x thetapointsamount

    # Convert to pixel coords from aart plot coords
    xaart = np.multiply(rarray , np.cos(thetarray))
    yaart = np.multiply(rarray , np.sin(thetarray))
    xprime = xaart + I0.shape[0] / 2
    yprime = yaart + I0.shape[0] / 2


    coords = np.array([xprime,yprime]).T
    peak = np.argmax(interp(coords), 1)
    return peak * ((limits*2) / I0.shape[0]) # units of Rg
    


def profile_simple(params, r, y):
    brightparams = [
        params['redshift'],
        params['nu0'],
        params['mass'],
        params['scaleh'],
        params['thetab'],
        params['beta'], 
        params['Rie'],
        params['Bchoice'],
        params['rb0'],
        params['nth0'],
        params['te0'],
        params['pden'],
        params['ptemp'],
    ]
    yfit = ilp.profile(r, brightparams[0],brightparams[1] * ilp.Hz,brightparams[2] * ilp.grams,brightparams[3],brightparams[4] * ilp.rads,brightparams[5],
        brightparams[6],brightparams[7],brightparams[8],brightparams[9] * ilp.cmcubed,brightparams[10] * ilp.kelv,brightparams[11], brightparams[12]).value
    return yfit - y


def best_fit(redshift, nu0=tls_nu0.value,mass=tls_mass.value, scale_height=tls_scale_height, theta_b=tls_theta_b.value, 
            beta=tls_beta,R_ie=tls_R_ie, Bchoice=tls_Bchoice, rb_0=tls_rb_0,p_dens=tls_p_dens,p_temp=tls_p_temp):
    
    params = Parameters()
    params.add('redshift', value=redshift, vary=False)
    params.add('nu0', value=nu0, vary=False)
    params.add('mass', value=mass, vary=False)
    params.add('scaleh', value=scale_height, vary=False)
    params.add('thetab', value=theta_b, vary=False)
    params.add('beta', value=beta, vary=False)
    params.add('Rie', value=R_ie, vary=False)
    params.add('Bchoice', value=Bchoice, vary=False)
    params.add('rb0', value=rb_0, vary=False)
    params.add('nth0', value=tls_n_th0.value)
    params.add('te0', value=tls_t_e0.value)
    params.add('pden', value=p_dens, vary=False)
    params.add('ptemp', value=p_temp, vary=False)

    brightparams = [
        tls_redshift,
        tls_nu0, # keep this at 230?
        tls_mass,
        tls_scale_height,
        tls_theta_b,
        tls_beta, 
        tls_R_ie,
        tls_Bchoice,
        tls_rb_0,
        tls_n_th0,
        tls_t_e0,
        tls_p_dens,
        tls_p_temp,
    ]
    r = np.array([2, 15, 50])
    y = ilp.profile(r, brightparams[0],brightparams[1],brightparams[2],brightparams[3],brightparams[4],brightparams[5],
        brightparams[6],brightparams[7],brightparams[8],brightparams[9],brightparams[10],brightparams[11], brightparams[12]).value
    fitted_params = minimize(profile_simple, params, args=(r,y), method='least_squares')
    return fitted_params.params['nth0'].value, fitted_params.params['te0'].value
