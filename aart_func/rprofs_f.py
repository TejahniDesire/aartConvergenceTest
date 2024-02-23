import numpy as np

from aart_func import *
from params import *
from functools import partial
from astropy import constants as const
from astropy import units as u
from lmfit import Parameters, minimize, fit_report
from scipy.special import gamma

# All r in units of Rg
# Constants

# TODO Units of the images pixel value
# TODO find values to change j coeff peak to occur between 3-20rg

inoisy_path = "/home/tej/Desktop/Code_Stuff/Repositories/aart/"
# Constants
G = const.G.cgs
c = const.c.cgs
me = const.m_e.cgs
mp = const.m_p.cgs
e = const.e.esu
kB = const.k_B.cgs
h = const.h.cgs

# Units
cm = u.cm
cmcubed = 1 / (cm ** 3)
kelv = 1 * u.K
grams = 1 * u.g
gauss = 1 * u.cm ** (-1 / 2) * u.g ** (1 / 2) * 1 / (1 * u.s)  # Gauss units in cgs
rads = 1 * u.rad
Hz = 1 * u.Hz
specific_int_units = u.Fr ** 2 * u.Hz * u.s / (u.cm ** 3)

'''KWARGS------------------------------'''
'''Physical'''
kw_nu0 = 230e9 * Hz
# kw_nu0 = 86e9 * Hz
# kw_nu0 = 345e9 * Hz
kw_mass = (MMkg * u.kg).to(u.g)  # m87 mass as put in AART
kw_scale_height = .5  # Angle
kw_theta_b = 60 * (np.pi / 180) * rads
kw_beta = 1
kw_r_ie = 10

'''Math Coeff'''
kw_rb_0 = 2
kw_n_th0 = 1.0726e+05 * cmcubed
kw_t_e0 = 1.2428e+11 * kelv
# TODO: define below
kw_bv_0 = 8.131273135591028  # set from fitting b_func power law

# kw_n_th0 = 1.23e6 * cmcubed


'''Exponents'''
kw_p_temp = -.84
kw_p_dens = -.7
# kw_p_temp = -.05
# kw_p_dens = -.05
# kw_p_bv = -0.8499999999965503  # set from fitting b_func power law
kw_p_bv = -.3 # set from fitting b_func power law

'''B_field '''
b_0 = None
p_b = None
kw_Bchoice = 0

'''Power Model'''
p_power_model = 3  # Others: 3.5, 7.0
gamma_1 = 1
gamma_2 = 10 ** 6

'''Noisy'''
kw_nscale = .4  # 0.4
kw_nnoisykey = 1  # 1 for on, 0 for off
kw_tnoisykey = 0
kw_bnoisykey = 0

'''Parameters'''
kw_emodelkey = 0
kw_bkey = 0

kw_brightparams = {
    "nu0": kw_nu0,  # 1
    "mass": kw_mass,  # 2
    "scale_height": kw_scale_height,  # 3
    "theta_b": kw_theta_b,  # 4
    "beta": kw_beta,  # 5
    "r_ie": kw_r_ie,  # 6
    "rb_0": kw_rb_0,  # 7
    "n_th0": kw_n_th0,  # 8
    "t_e0": kw_t_e0,  # 9
    "p_dens": kw_p_dens,  # 10
    "p_temp": kw_p_temp,  # 11
    "nscale": kw_nscale,  # 12
}

kw_funckeys = {
    "emodelkey": kw_emodelkey,
    "bkey": kw_bkey,
    "nnoisykey": kw_nnoisykey,
    "tnoisykey": kw_tnoisykey,
    "bnoisykey": kw_bnoisykey
}


# B Function through best fit --------------------------------
def full_b_func(r, mass=kw_mass, beta=kw_beta, rb_0=kw_rb_0, n_th0=kw_n_th0, p_dens=kw_p_dens):
    """Full magnetic field strength equation for electons 

    Args:
        r (_Float_): radial distance from black hole
        mass (_Float_, optional): Mass of black hole. Defaults to kw_mass.
        beta (_Float_, optional): Ratio of gass pressure to Ion pressure. Defaults to kw_beta.
        rb_0 (_Float_, optional): Radius at which n = nth0. Defaults to kw_rb_0.
        n_th0 (_Float_, optional): Proportionality constant for density power law. Defaults to kw_n_th0.
        p_dens (_Float_, optional): Exponent for density power law. Defaults to kw_p_dens.

    Returns:
        _type_: Magnetic field strength at that radius
    """
    nth = nth_func(r, mass, rb_0, n_th0, p_dens)
    if n_th0.unit != u.cm ** -3:
        raise ValueError('n_th0 needs units of number density')
    return np.sqrt(nth * 8 * np.pi * mp * c ** 2 / (6 * beta * r))


def b_simple_func(params, r, y):
    rg = params['rg']
    rb = params['rb']
    b_0 = params['b_0']
    p_b = params['p_b']
    y_fit = b_0 * (r * rg / rb) ** p_b
    return y_fit - y


# TODO Check if least_squares assumes squares
def b_best_fit(mass=kw_mass, beta=kw_beta, rb_0=kw_rb_0, n_th0=kw_n_th0, p_dens=kw_p_dens):
    # Create r and y data
    r = np.array([2, 15, 50])
    y = full_b_func(r, mass, beta, rb_0, n_th0, p_dens).value
    params = Parameters()
    params.add('b_0', value=1)
    params.add('p_b', value=1)
    params.add('rg', value=rg_func(mass).value, vary=False)
    params.add('rb', value=rb_func(mass).value, vary=False)
    fitted_params = minimize(b_simple_func, params, args=(r, y), method='least_squares')
    return fitted_params.params['b_0'].value, fitted_params.params['p_b'].value


def set_b_params(mass=kw_mass, beta=kw_beta, rb_0=kw_rb_0, n_th0=kw_n_th0, p_dens=kw_p_dens):
    global b_0
    global p_b
    b_0, p_b = b_best_fit(mass, beta, rb_0, n_th0, p_dens)
    b_0 = b_0 * gauss


def b_func_power(r, mass=kw_mass, rb_0=kw_rb_0):
    rg = rg_func(mass)
    rb = rb_func(mass, rb_0)
    return b_0 * (r * rg / rb) ** p_b
    # return np.sqrt(nth * 8 * np.pi * mp * c ** 2 / (6 * beta * r))


# TODO: pick better name, fix parameters
def b_func_power_variable(r, mass=kw_mass, rb_0=kw_rb_0, bv_0=kw_bv_0, p_bv=kw_p_bv):
    rg = rg_func(mass)
    rb = rb_func(mass, rb_0)
    bv_0 = bv_0 * gauss
    return bv_0 * (r * rg / rb) ** p_bv


# -----------------------------------------------------------------------------------------------------


def rg_func(mass=kw_mass):
    return G * mass / c ** 2


def rb_func(mass=kw_mass, rb_0=kw_rb_0):
    """Value at which the power laws take on the value of their constants of proportionality 

    Args:
        mass (_type_, optional): _description_. Defaults to kw_mass.
        rb_0 (_type_, optional): _description_. Defaults to kw_rb_0.

    Returns:
        _type_: _description_
    """
    return rb_0 * G * mass / (c ** 2)


def te_func(r, mass=kw_mass, rb_0=kw_rb_0, t_e0=kw_t_e0, p_temp=kw_p_temp):
    """Temperature as a function of distance from the black hole

    Args:
        r (Float): radial distance from black hole
        mass (Float, optional): Mass of black hole. Defaults to kw_mass.
        rb_0 (Float, optional): Radius at which n = nth0. Defaults to kw_rb_0.
        t_e0 (Float, optional): Proportionality constant for temperature power law. Defaults to kw_t_e0.
        p_temp (Float, optional): Exponent for temperature power law. Defaults to kw_p_temp.

    Returns:
        _type_: Electron temperature at that distance
    """
    if t_e0.unit != u.K:
        raise ValueError('n_th0 needs units of Kelvin')
    rg = rg_func(mass)
    rb = rb_func(mass, rb_0)
    return t_e0 * (r * rg / rb) ** p_temp


def theta_e_func(temp):
    """Dimensionless temperature value

    Args:
        r (Float): radial distance from black hole
        mass (Float, optional): Mass of black hole. Defaults to kw_mass.
        rb_0 (Float, optional): Radius at which n = nth0. Defaults to kw_rb_0.
        t_e0 (Float, optional): Proportionality constant for temperature power law. Defaults to kw_t_e0.
        p_temp (Float, optional): Exponent for temperature power law. Defaults to kw_p_temp.

    Returns:
        _type_: Electron Dimensionless temperature at that distance 
    """
    return (kB * temp / (me * c ** 2)).to(u.dimensionless_unscaled)


def nth_func(r, mass=kw_mass, rb_0=kw_rb_0, n_th0=kw_n_th0, p_dens=kw_p_dens):
    """Density at as a function of distance from the black hole

    Args:
        r (Float): radial distance from black hole
        mass (Float, optional): Mass of black hole. Defaults to kw_mass.
        rb_0 (Float, optional): Radius at which n = nth0. Defaults to kw_rb_0.
        n_th0 (Float, optional): Proportionality constant for density power law. Defaults to kw_n_th0.
        p_dens (Float, optional): Exponent for density power law. Defaults to kw_p_dens.

    Returns:
        _type_: Density at that radial distance
    """
    if n_th0.unit != u.cm ** -3:
        raise ValueError('n_th0 needs units of number density')
    rg = rg_func(mass)
    rb = rb_func(mass, rb_0)
    return n_th0 * (r * rg / rb) ** p_dens


def b_func_true(beta, r_ie, theta_e, nth):
    """Full magnetic field equation

    Args:
        r (Float): radial distance from black hole
        mass (Float, optional): Mass of black hole. Defaults to kw_mass.
        beta (Float, optional): Ratio of gass pressure to Ion pressure. Defaults to kw_beta.
        r_ie (Float, optional): _description_. Defaults to kw_r_ie.
        rb_0 (Float, optional): Radius at which n = nth0. Defaults to kw_rb_0.
        n_th0 (Float, optional): Proportionality constant for density power law. Defaults to kw_n_th0.
        t_e0 (Float, optional): Proportionality constant for temperature power law. Defaults to kw_t_e0.
        p_dens (Float, optional): Exponent for density power law. Defaults to kw_p_dens.
        p_temp (Float, optional): Exponent for temperature power law. Defaults to kw_p_temp.


    Returns:
        _type_: Magnetic field value in gauss at that radial distance
    """
    return np.sqrt(8 * np.pi * me * c ** 2 * (1 + r_ie) * beta ** -1 * theta_e * nth)


def nu_c_func(b_field, theta_e, theta_b):
    """ Frequnecy scaler

    Args:
        r (Float): radial distance from black hole
        mass (Float, optional): Mass of black hole. Defaults to kw_mass.
        theta_b (Float, optional): Angle between magnetic field and wave vector. Defaults to kw_theta_b.
        beta (Float, optional): Ratio of gass pressure to Ion pressure. Defaults to kw_beta.
        r_ie (Float, optional): _description_. Defaults to kw_r_ie.
        Bchoice (Int, optional): True magnetic field equaiton 0 or Power law magnetic field 1. Defaults to kw_Bchoice.
        rb_0 (Float, optional): Radius at which n = nth0. Defaults to kw_rb_0.
        n_th0 (Float, optional): Proportionality constant for density power law. Defaults to kw_n_th0.
        t_e0 (Float, optional): Proportionality constant for temperature power law. Defaults to kw_t_e0.
        p_dens (Float, optional): Exponent for density power law. Defaults to kw_p_dens.
        p_temp (Float, optional): Exponent for temperature power law. Defaults to kw_p_temp.


    Returns:
        _type_: Frequnecy scaler
    """

    nu_b = (e * b_field / (2 * np.pi * me * c)).to(u.Hz)
    return 3 / 2 * nu_b * np.sin(theta_b) * theta_e ** 2


def synchrotron_func_I(x):
    """Synchrotron emission function
    Args:
        x (_type_): dimensionless frequnecy counter

    Returns:
        _type_: synchrotron emisison at that frequency
    """
    return 2.5651 * (1 + 1.92 * (x ** (-1 / 3)) + (0.9977 * x ** (-2 / 3))) * np.exp(-1.8899 * x ** (1 / 3))


def synchrotron_func_Q(x):
    """Synchrotron emission function
    Args:
        x (_type_): dimensionless frequnecy counter

    Returns:
        _type_: synchrotron emisison at that frequency
    """
    return 2.5651 * (1 + 0.932 * (x ** (-1 / 3)) + (0.4998 * x ** (-2 / 3))) * np.exp(-1.8899 * x ** (1 / 3))


def synchrotron_func_V(x):
    """Synchrotron emission function
    Args:
        x (_type_): dimensionless frequnecy counter

    Returns:
        _type_: synchrotron emisison at that frequency
    """
    return (1.8138 * (x ** -1) + 3.423 * (x ** (-2 / 3)) + 0.02955 * x ** (-1 / 2) + 2.0377 * x ** (-1 / 3)) * np.exp(
        -1.8899 * x ** (1 / 3))


# beginning of Inoisy----------------------------------------------------------------------------------


def inoisy_radius():
    lowerbound = -30
    upperbound = 30
    gridsize = 512
    Xs = np.arange(lowerbound, upperbound, (upperbound - lowerbound) / gridsize)
    Ys = np.arange(lowerbound, upperbound, (upperbound - lowerbound) / gridsize)

    xx, yy = np.meshgrid(Xs, Ys)
    xx, yy = np.meshgrid(Xs, Ys)
    return np.sqrt(xx ** 2 + yy ** 2), Xs, Ys


def inoisy_interp(envelope, scale):  # given an envelope, return noisy version to be evaluated at x and y grid
    GRF = np.load(inoisy_path + "Inoisysnapshot.npy")
    radius, Xs, Ys = inoisy_radius()
    density = envelope * np.exp(scale * GRF - scale ** 1 / 2)
    return RegularGridInterpolator((Xs, Ys), density, fill_value=1, bounds_error=False, method='linear')


def inoisy_value(x, y, envelope, scale, units):  # return value of noisy parameter
    """

    :param x:
    :param y:
    :param envelope:
    :param scale:
    :param units:
    :return:
    """
    interpolation = inoisy_interp(envelope, scale)
    return interpolation(np.vstack([x, y]).T) * units


# Ultrarelativistic
def thermal_profile(coords, redshift, cosAng, bp=kw_brightparams, fk=kw_funckeys):
    """

    Calculate the radial profile emission according to a thermal distribution
    :param coords: dictionary containg "r" "x" and "y" values for each pixel
    :param redshift: radius's corrsponding redshift
    :param bp: dictionary with the following key-value pairs
        nu0: Obeservation Frequency. Defaults to kw_nu0.
        mass: Mass of black hole. Defaults to kw_mass.
        scale_height (Float, optional): Slope of acretion disk height vs radius. Defaults to kw_scale_height.
        theta_b (Float, optional): Angle between magnetic field and wave vector. Defaults to kw_theta_b.
        beta (Float, optional): Ratio of gass pressure to Ion pressure. Defaults to kw_beta.
        r_ie (Float, optional): _description_. Defaults to kw_r_ie.
        Bchoice (Int, optional): True magnetic field equaiton 0 or Power law magnetic field 1. Defaults to kw_Bchoice.
        rb_0 (Float, optional): Radius at which n = nth0. Defaults to kw_rb_0.
        n_th0 (Float, optional): Proportionality constant for density power law. Defaults to kw_n_th0.
        t_e0 (Float, optional): Proportionality constant for temperature power law. Defaults to kw_t_e0.
        p_dens (Float, optional): Exponent for density power law. Defaults to kw_p_dens.
        p_temp (Float, optional): Exponent for temperature power law. Defaults to kw_p_temp.
        nscale: scale factor of inoisy
    :param fk: dictionary with the following key-value pairs
        "absorbkey": absorption or not
        "emodelkey": thermal profile or power model
        "bkey":      true magnetic field model or power model
        "nnoisykey": unperturbed density or inoisy density
        "tnoisykey": unperturbed temperature or inoisy temperature
        "bnoisykey": unperturbed magnetic field or inoisy magnetic field
    :return: brightness temperature and specific intensity

    """

    rnoisy, Xs, Yx = inoisy_radius()
    # Temperature and Theta_e-------------------------------------------------------------------------------------------
    tempnoisy = te_func(rnoisy, bp["mass"], bp["rb_0"], bp["t_e0"], bp["p_temp"])
    te_noisy_funcs = {
        0: partial(te_func, coords["r"], bp["mass"], bp["rb_0"], bp["t_e0"], bp["p_temp"]),
        1: partial(inoisy_value, coords["x"], coords["y"], tempnoisy, bp["nscale"], kelv)
    }

    temp = te_noisy_funcs[fk["tnoisykey"]]()
    # temp[np.all(temp.value)] = 10 ** -3 * kelv
    theta_e = theta_e_func(temp)

    # Density-----------------------------------------------------------------------------------------------------------
    nthnoisy = nth_func(rnoisy, bp["mass"], bp["rb_0"], bp["n_th0"], bp["p_dens"])
    n_noisy__funcs = {
        0: partial(nth_func, coords["r"], bp["mass"], bp["rb_0"], bp["n_th0"], bp["p_dens"]),
        1: partial(inoisy_value, coords["x"], coords["y"], nthnoisy, bp["nscale"], cmcubed)
    }
    n = n_noisy__funcs[fk["nnoisykey"]]()

    # Magnetic Field----------------------------------------------------------------------------------------------------
    b_fieldnoisyfuncs = {
        0: partial(b_func_true, bp["beta"], bp["r_ie"], theta_e_func(tempnoisy), nthnoisy),
        1: partial(b_func_power, rnoisy, bp["mass"], bp["rb_0"]),
        2: partial(b_func_power_variable, rnoisy, bp["mass"], bp["rb_0"], kw_bv_0, kw_p_bv)
    }
    bfieldnoisy = b_fieldnoisyfuncs[fk["bkey"]]()

    b_field_funcs = {
        0: partial(b_func_true, bp["beta"], bp["r_ie"], theta_e, n),
        1: partial(b_func_power, coords["r"], bp["mass"], bp["rb_0"]),
        2: partial(b_func_power_variable, coords["r"], bp["mass"], bp["rb_0"], kw_bv_0, kw_p_bv)
    }

    b_field_noisy_funcs = {
        0: b_field_funcs[fk["bkey"]],
        1: partial(inoisy_value, coords["x"], coords["y"], bfieldnoisy, bp["nscale"], gauss)
    }
    b_field = b_field_noisy_funcs[fk["bnoisykey"]]()

    # ------------------------------------------------------------------------------------------------------------------
    nu_c = nu_c_func(b_field, theta_e, bp["theta_b"])

    nu = bp["nu0"] / redshift
    x = nu / nu_c

    runits = coords["r"] * rg_func(bp["mass"])
    # l' (l in the fluid frame)
    path_length_fluid = runits * bp["scale_height"] / cosAng

    # J Coeff Calculations Returns units of [u.erg / (u.cm ** 3 * u.s * u.Hz)]------------------------------------------
    jcoeff_I_fluid = n * e ** 2 * nu * synchrotron_func_I(x) / (2 * np.sqrt(3) * c * theta_e ** 2)

    # jcoeff_Q_fluid = n * e ** 2 * nu * synchrotron_func_Q(x) / (2 * np.sqrt(3) * c * theta_e ** 2)
    # jcoeff_V_fluid = 2 * n * e ** 2 * nu * synchrotron_func_Q(x) / (3 * np.sqrt(3) * np.tan(theta_b) * c * theta_e ** 3)

    # Absorption--------------------------------------------------------------------------------------------------------

    b_nu_fluid = ((2 * h * nu ** 3) / c ** 2) / (np.exp(h * nu / (theta_e * me * c ** 2)) - 1)

    acoeff_I_fluid = jcoeff_I_fluid / b_nu_fluid
    tau = acoeff_I_fluid * path_length_fluid

    specific_intensity_thin = path_length_fluid * jcoeff_I_fluid * redshift ** 3
    # brightness = ((c ** 2 / (2 * nu ** 2 * kB)) * specific_intensity_thin).to(u.K)

    # specific_intensity_thick = inplus * np.exp(-tau) + redshift ** 3 * b_nu_fluid * (
    #             1 - np.exp(- tau))
    specific_intensity_thick = redshift ** 3 * b_nu_fluid * (1 - np.exp(- tau))

    # Convert planck to brightness radial--------------------
    b_nu_fluid = brightness_temp(b_nu_fluid, bp["nu0"])
    # Packing radial curves--------------------
    jcoeff_I_fluid_packed = jcoeff_I_fluid.reshape(1, jcoeff_I_fluid.shape[0])
    theta_e = theta_e.reshape(1, theta_e.shape[0])
    n = n.reshape(1, n.shape[0])
    b_field = b_field.reshape(1, b_field.shape[0])
    b_nu_fluid = b_nu_fluid.reshape(1, b_nu_fluid.shape[0])
    acoeff_I_fluid = acoeff_I_fluid.reshape(1, acoeff_I_fluid.shape[0])
    tau_curve = tau.reshape(1, tau.shape[0])
    r = coords["r"].reshape(1, coords["r"].shape[0])

    specific_intensity_thin_packed = specific_intensity_thin.reshape(1, specific_intensity_thin.shape[0])
    specific_intensity_thick_packed = specific_intensity_thick.reshape(1, specific_intensity_thick.shape[0])

    # full_profiles = np.concatenate([r, theta_e.value, n.value, b_field.value, b_nu_fluid.value,
    #                                 acoeff_I_fluid.value, tau_curve.value], axis=0)

    full_profiles = np.concatenate([r, theta_e.value, n.value, b_field.value, b_nu_fluid.value,
                                    acoeff_I_fluid.value, tau_curve.value, jcoeff_I_fluid_packed.value], axis=0)

    full_profiles_units = [str(theta_e.unit), str(n.unit), str(b_field.unit), str(b_nu_fluid.unit),
                           str(acoeff_I_fluid.unit), str(tau_curve.unit)]
    # print("b_0: ", b_0)
    # print("p: ", p_b)
    return specific_intensity_thin, specific_intensity_thick, tau, full_profiles, full_profiles_units


# TODO: remove redshift
def brightness_temp(specific_intensity, nu0):
    return ((c ** 2 / (2 * nu0 ** 2 * kB)) * specific_intensity).to(u.K)


def tau_tester(coords, bp=kw_brightparams, fk=kw_funckeys):
    redshift = sqrt(1 - (2 * G * bp['mass'] / (c ** 2 * coords['r'] * rg_func(bp["mass"]))))
    rnoisy, Xs, Yx = inoisy_radius()
    # Temperature and Theta_e---------------------
    tempnoisy = te_func(rnoisy, bp["mass"], bp["rb_0"], bp["t_e0"], bp["p_temp"])
    te_noisy_funcs = {
        0: partial(te_func, coords["r"], bp["mass"], bp["rb_0"], bp["t_e0"], bp["p_temp"]),
        1: partial(inoisy_value, coords["x"], coords["y"], tempnoisy, bp["nscale"], kelv)
    }

    temp = te_noisy_funcs[fk["tnoisykey"]]()
    theta_e = theta_e_func(temp)

    # Density---------------------
    nthnoisy = nth_func(rnoisy, bp["mass"], bp["rb_0"], bp["n_th0"], bp["p_dens"])
    n_noisy__funcs = {
        0: partial(nth_func, coords["r"], bp["mass"], bp["rb_0"], bp["n_th0"], bp["p_dens"]),
        1: partial(inoisy_value, coords["x"], coords["y"], nthnoisy, bp["nscale"], cmcubed)
    }
    n = n_noisy__funcs[fk["nnoisykey"]]()

    # Magnetic Field---------------------

    b_fieldnoisyfuncs = {
        0: partial(b_func_true, bp["beta"], bp["r_ie"], theta_e_func(tempnoisy), nthnoisy),
        1: partial(b_func_power, rnoisy, bp["mass"], bp["rb_0"])
    }
    bfieldnoisy = b_fieldnoisyfuncs[fk["bkey"]]()

    b_field_funcs = {
        0: partial(b_func_true, bp["beta"], bp["r_ie"], theta_e, n),
        1: partial(b_func_power, coords["r"], bp["mass"], bp["rb_0"])
    }

    b_field_noisy_funcs = {
        0: b_field_funcs[fk["bkey"]],
        1: partial(inoisy_value, coords["x"], coords["y"], bfieldnoisy, bp["nscale"], gauss)
    }
    b_field = b_field_noisy_funcs[fk["bnoisykey"]]()

    # -------------------------------------

    nu_c = nu_c_func(b_field, theta_e, bp["theta_b"])
    nu = bp["nu0"] / redshift
    x = nu / nu_c

    # J Coeff Calculations Returns units of [u.erg / (u.cm ** 3 * u.s * u.Hz)]

    jcoeff_I = n * e ** 2 * nu * synchrotron_func_I(x) / (2 * np.sqrt(3) * c * theta_e ** 2)
    # jcoeff_Q = n * e ** 2 * nu * synchrotron_func_Q(x) / (2 * np.sqrt(3) * c * theta_e ** 2)
    # jcoeff_V = 2 * n * e ** 2 * nu * synchrotron_func_Q(x) / (3 * np.sqrt(3) * np.tan(theta_b) * c * theta_e ** 3)

    # Absorption---------------

    b_nu = (2 * h * nu ** 3) / c ** 2 * (np.exp(h * nu / (theta_e * me * c ** 2))) ** -1
    absorptionCoeff = jcoeff_I / b_nu
    runits = coords["r"] * rg_func(bp["mass"])
    tau = absorptionCoeff * runits * bp["scale_height"]

    # WARNING,REDSHIFT IN NU

    # specific_intensity_thin = runits * bp["scale_height"] * jcoeff_I * redshift ** 3
    # specific_intensity_thick = inplus * np.exp(-tau) + redshift ** 3 * b_nu * (
    #             1 - np.exp(- tau))
    # brightness = ((c ** 2 / (2 * nu ** 2 * kB)) * specific_intensity_thin).to(u.K)
    # brightness = ((c ** 2 / (2 * nu ** 2 * kB)) * specific_intensity_absorp).to(u.K)
    return tau


# Powerlaw, produce table mathematica
def power_profile():
    pass


# Assume I is ndarray without astropy units, but in kelvin
def total_jy(I, nu, mass):
    """Calculate totalt jasnky flux of a photon ring

    Args:
        I (ndarray): Ndarray of brightness temperature values
        nu (Float): Observation frequency
        mass (Float): Black hole mass

    Returns:
        _type_: _description_
    """
    I = I * u.K
    nu = nu * u.Hz
    mass = mass * u.g
    one_M = rg_func(mass).to(u.m)
    M2rads = np.arctan(one_M.value / dBH)
    rads2pxls = (M2rads * 2 * limits) / I.shape[0]  # total amount M length units rads / total pixels
    return (I * nu ** 2 * (2 * kB / c ** 2)).to(u.Jy).sum() * rads2pxls ** 2  # was this orifianlly in per radians?


def ring_radius(I0):
    """Calculates radius of each photon ring

    Args:
        I0 (ndarray): Ndarray of brightness temperature values

    Returns:
        _type_: Photon radius
    """
    plx2rg = (limits * 2) / I0.shape[0]  # pixels over M

    horizontalline = np.zeros(2)
    verticalline = np.zeros(2)

    length = I0.shape[0]
    midpoint = int(length / 2)

    horizontalline[0] = I0[midpoint, 0:midpoint].argmax()
    horizontalline[1] = midpoint + I0[midpoint, midpoint:length].argmax()

    verticalline[0] = I0[0:midpoint, midpoint].argmax()
    verticalline[1] = midpoint + I0[midpoint:length, midpoint].argmax()

    return np.mean([verticalline[1] - verticalline[0], horizontalline[1] - horizontalline[0]]) * plx2rg / 2


# Returns index of first occurance of below percent diff
def ring_convergance(ring1, ring2, percent_diff):
    diff = (np.abs(ring1 - ring2) / ring1) * 100
    return np.argmax(diff <= percent_diff)
