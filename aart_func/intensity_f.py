from aart_func import *
from params import *

def Delta(r,a):
    """
    Calculates the Kerr metric function \Delta(t)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return r**2-2*r+a**2

def PIF(r,a):
    """
    Calculates PI(r) (Eq. B6 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return (r**2+a**2)**2-a**2*Delta(r,a)

def urbar(r,a):
    """
    Calculates the r (contravariant) component of the four velocity for radial infall
    (Eq. B34b P1)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return -np.sqrt(2*r*(r**2+a**2))/(r**2)

def Omegabar(r,a):
    """
    Calculates the angular velocity of the radial infall
    (Eq. B32a P1)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return (2*a*r)/PIF(r,a)

def Omegahat(r,a,laux):
    """
    Calculates the angular velocity of the sub-Keplerian orbit
    (Eq. B39 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return (a+(1-2/r)*(laux-a))/(PIF(r,a)/(r**2)-(2*a*laux)/r)

def uttilde(r, a,urT,OT):
    """
    Calculates the t (contravariant) component of the general four velocity
    (Eq. B52 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param urT: r (contravariant) component of the general four velocity
    :param OT: Angular velocity of the general four velocity
    """
    return np.sqrt((1 + urT**2*r**2/Delta(r,a))/(1-(r**2+a**2)*OT**2-(2/r)*(1-a*OT)**2))

def Ehat(r,a,laux):
    """
    Calculates the orbital energy of the sub-Keplerian flow
    (Eq. B44a P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param laux: sub-Keplerian specific angular momentum
    """
    return np.sqrt(Delta(r,a)/(PIF(r,a)/(r**2)-(4*a*laux)/r-(1-2/r)*laux**2))

def nuhat(r,a,laux,Ehataux):
    """
    Calculates the radial velocity of the sub-Keplerian flow
    (Eq. B45 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param laux: sub-Keplerian specific angular momentum
    :param Ehataux: sub-Keplerian orbital energy
    """
    return r/Delta(r,a)*np.sqrt(np.abs(PIF(r,a)/(r**2)-(4*a*laux)/r-(1-2/r)*laux**2-Delta(r,a)/(Ehataux**2)))

def lhat(r,a):
    """
    Calculates the rspecific angular momentum of the sub-Keplerian flow
    (Eq. B44b P1)
    :param r: radius of the source
    :param a: spin of the black hole
    """
    return sub_kep*(r**2+a**2-2*a*np.sqrt(r))/(np.sqrt(r)*(r-2)+a)

def Rint(r,a,lamb,eta):
    """
    Evaluates the "radial potential", for calculating the redshift factor for infalling material
    :param r: radius of the source
    :param a: spin of the black hole
    :param lamb: angular momentum
    :param eta: carter constant

    :return: radial potential evaluated at the source
    """
    #Eqns (P2 5)
    return (r**2 + a**2 - a*lamb)**2 - (r**2 - 2*r + a**2)*(eta + (lamb - a)**2)

def gDisk(r,a,b,lamb,eta):
    """
    Calculates the redshift factor for a photon outside the inner-most stable circular orbit(isco) (assume circular orbit)
    (Eq. B13 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param lamb: angular momentum
    :param eta: Carter constant

    :return: the redshift factor associated with the ray
    """

    OH=Omegahat(r,a,lhat(r,a))
    OT=OH+(1-betaphi)*(Omegabar(r,a)-OH)
    ur=(1-betar)*urbar(r,a)
    ut=uttilde(r,a,ur,OT)
    uphi=ut*OT

    return 1/(ut*(1-b*np.sign(ur)*sqrt(np.abs(Rint(r,a,lamb,eta)*ur**2))/Delta(r,a)/ut-lamb*uphi/ut))

def gGas(r,a,b,lamb,eta):
    """
    Calculates the redshift factor for a photon inside the isco (assume infalling orbit)
    (Eq. B13 P1)
    :param r: radius of the source
    :param a: spin of the black hole
    :param b: sign for the redshift
    :param lamb: angular momentum
    :param eta: carter constant

    :return: the redshift factor associated with the ray
    """
    #Calculate radius of the inner-most stable circular orbit
    isco=rms(a)

    lms=lhat(isco,a)
    OH=Omegahat(r,a,lms)
    OT=OH+(1-betaphi)*(Omegabar(r,a)-OH)

    Ems=Ehat(isco,a,lms)
    urhat=-Delta(r,a)/(r**2)*nuhat(r, a, lms ,Ems)*Ems
    ur=urhat+(1-betar)*(urbar(r,a)-urhat)
    ut=uttilde(r,a,ur,OT)
    uphi=OT*ut

    return 1/(ut*(1-b*np.sign(ur)*sqrt(np.abs(Rint(r,a,lamb,eta)*ur**2))/Delta(r,a)/ut-lamb*uphi/ut))


#calculate the observed brightness for a purely radial profile
def bright_radial(grid,mask,redshift_sign,a,rs,isco,thetao,brightparams,funckeys,phi):
    """
    Calculate the brightness of a rotationally symmetric disk
    (Eq. 50 P1)
    :param grid: alpha and beta grid on the observer plane on which we evaluate the observables
    :param mask: mask out the lensing band, see lb_f.py for detail
    :param redshift_sign: sign of the redshift
    :param a: black hole spin
    :param rs: source radius
    :param isco: radius of the inner-most stable circular orbit
    :param thetao: observer inclination
    :param brightparams: physical parameters for the plasma
    :param funckeys: keys on deciding between equation models
    :param phi: source angle
    :param inplus: specific intensity from previous iterations

    :return: image of a lensed equitorial source with only radial dependence.
    """
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]

    lamb,eta = rt.conserved_quantities(alpha,beta,thetao,a)

    si_thin = np.zeros(rs.shape[0])
    si_thick = np.zeros(rs.shape[0])
    # full_profiles = np.zeros(shape=(7,rs.shape[0]))
    full_profiles = np.zeros(shape=(8, rs.shape[0]))
    tau = np.zeros(rs.shape[0])
    cosAngReturn = np.zeros(rs.shape[0])
    redshift_sign = redshift_sign[mask]

    x_aux=rs*np.cos(phi)
    y_aux=rs*np.sin(phi)

    # rs_inner = rs[rs>=isco]
    # rs_outer = rs[rs<isco]
    redshift_inner = gDisk(rs[rs>=isco],a,redshift_sign[rs>=isco],lamb[rs>=isco],eta[rs>=isco])
    redshift_outter = gGas(rs[rs<isco],a,redshift_sign[rs<isco],lamb[rs<isco],eta[rs<isco])

    CosAng_inner = CosAng(rs[rs>=isco],a,redshift_inner,lamb[rs>=isco],eta[rs>=isco])
    CosAng_outter = CosAng(rs[rs<isco],a,redshift_outter,lamb[rs<isco],eta[rs<isco])

    # r_p = 1 + np.sqrt(1 - a ** 2)
    # CosAng_inner[rs_inner<=r_p] = 1
    # CosAng_outter[rs_outer<=r_p] = 1

    # CosAng_inner[np.isnan(CosAng_inner)] = 1
    # CosAng_outter[np.isnan(CosAng_outter)] = 1

    ilp.set_b_params(brightparams["mass"],brightparams["beta"],brightparams["rb_0"],brightparams["n_th0"],brightparams["p_dens"])

    coords_inner = {
        "r": rs[rs>=isco],
        "x": x_aux[rs>=isco],
        "y": y_aux[rs>=isco]
    }
    coords_outter = {
        "r": rs[rs<isco],
        "x": x_aux[rs<isco],
        "y": y_aux[rs<isco]
    }

    emissionmodel = {
        0: ilp.thermal_profile,
        1: ilp.power_profile
    }
    si_thin[rs>=isco], si_thick[rs>=isco], tau[rs>=isco], full_profiles[:,rs>=isco], full_profiles_unit = emissionmodel[
        funckeys["emodelkey"]](coords_inner,redshift_inner,CosAng_inner,brightparams,funckeys)
    si_thin[rs<isco], si_thick[rs<isco], tau[rs<isco], full_profiles[:,rs<isco], full_profiles_unit = emissionmodel[
        funckeys["emodelkey"]](coords_outter,redshift_outter,CosAng_outter,brightparams,funckeys)

    r_p = 1+np.sqrt(1-a**2)
    # si_thin[rs<=r_p] = 0
    si_thin[rs <= r_p] = -np.inf

    cosAngReturn[rs>=isco] = CosAng_inner
    cosAngReturn[rs<isco] = CosAng_outter

    testa = np.sqrt(eta)
    testb = np.sqrt(eta+a**2*np.cos(thetao)**2-lamb**2/(np.tan(thetao)**2))

    return si_thin, si_thick, tau, full_profiles, full_profiles_unit, cosAngReturn, testa, testb, eta


#calculate the observed brightness for an arbitrary profile, passed in as the interpolation object
#but ignoring the time delay due to lensing
def fast_light(grid,mask,redshift_sign,a,isco,rs,th,interpolation,thetao):
    """
    Calculate the black hole image ignoring the time delay due to lensing or geometric effect
    (Eq. 116 P1)
    :param grid: alpha and beta grid on the observer plane on which we evaluate the observables
    :param mask: mask out the lensing band, see lb_f.py for detail
    :param redshift_sign: sign of the redshift
    :param a: black hole spin
    :param isco: radius of the inner-most stable circular orbit
    :param rs: source radius
    :param th: source angle, polar coordinate
    :param interpolation: 2 dimensional brightness function of the source, interpolation object
    :param thetao: observer inclination

    :return: image of a lensed equitorial source with only radial dependence. 
    """
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]
    th = th[mask]
    lamb,eta = rt.conserved_quantities(alpha,beta,thetao,a)
    brightness = np.zeros(rs.shape[0])
    redshift_sign = redshift_sign[mask]

    x_aux=rs*np.cos(th)
    y_aux=rs*np.sin(th)

    brightness[rs>=isco]= gDisk(rs[rs>=isco],a,redshift_sign[rs>=isco],lamb[rs>=isco],eta[rs>=isco])**gfactor*interpolation(np.vstack([x_aux[rs>=isco],y_aux[rs>=isco]]).T)
    brightness[rs<isco]= gGas(rs[rs<isco],a,redshift_sign[rs<isco],lamb[rs<isco],eta[rs<isco])**gfactor*interpolation(np.vstack([x_aux[rs<isco],y_aux[rs<isco]]).T)

    r_p = 1+np.sqrt(1-a**2)
    brightness[rs<=r_p] = 0

    I = np.zeros(mask.shape)
    I[mask] = brightness
    return(I)


#calculate the observed brightness for an arbitrary, evolving profile, passed in as the interpolation object
def slow_light(grid,mask,redshift_sign,a,isco,rs,th,ts,interpolation,thetao):
    """
    Calculate the black hole image including the time delay due to lensing and geometric effect
    (Eq. 50 P1)

    :param grid: alpha and beta grid on the observer plane on which we evaluate the observables
    :param mask: mask out the lensing band, see lb_f.py for detail
    :param redshift_sign: sign of the redshift
    :param a: black hole spin
    :param isco: radius of the inner-most stable circular orbit
    :param rs: source radius
    :param th: source angle, polar coordinate
    :param ts: time of emission at the source
    :param interpolation: a time series of 2 dimensional brightness function of the source, 3d interpolation object
    :param thetao: observer inclination

    :return: image of a lensed equitorial source with only radial dependence. 
    """
    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]
    th = th[mask]
    ts = ts[mask]
    lamb,eta = rt.conserved_quantities(alpha,beta,thetao,a)
    brightness = np.zeros(rs.shape[0])
    redshift_sign = redshift_sign[mask]

    x_aux=rs*np.cos(th)
    y_aux=rs*np.sin(th)

    brightness[rs>=isco]= gDisk(rs[rs>=isco],a,redshift_sign[rs>=isco],lamb[rs>=isco],eta[rs>=isco])**gfactor*interpolation(np.vstack([ts[rs>=isco],x_aux[rs>=isco],y_aux[rs>=isco]]).T)
    brightness[rs<isco]= gGas(rs[rs<isco],a,redshift_sign[rs<isco],lamb[rs<isco],eta[rs<isco])**gfactor*interpolation(np.vstack([ts[rs<isco],x_aux[rs<isco],y_aux[rs<isco]]).T)

    r_p = 1+np.sqrt(1-a**2)
    brightness[rs<=r_p] = 0

    I = np.zeros(mask.shape)
    I[mask] = brightness
    return(I)


def br(supergrid0,mask0,N0,rs0,sign0,supergrid1,mask1,N1,rs1,sign1,supergrid2,mask2,N2,rs2,sign2,brightparams,funckeys,phi012):
    """
    Calculate and save the radial brightness profile
    """

    # Warning: si_thin is a decomposition while si_thick is cumulative

    full_intensity = np.zeros(rs2.shape)
    # I2-------------
    rs2 = rs2[mask2]
    phi2 = phi012[2][mask2]

    si_thin2, si_thick2, tau2mask2, full_profiles2, full_profiles_unit, cosAngReturn2, testa, testb, eta= bright_radial(
        supergrid2,mask2,sign2,spin_case,rs2,isco,thetao,brightparams,funckeys,phi2)

    cosAngReturn2Full = np.zeros(mask2.shape)
    cosAngReturn2Full[mask2] = cosAngReturn2

    I2_thin = np.zeros(mask2.shape)
    I2_thin[mask2] = si_thin2
    I2_thick = np.zeros(mask2.shape)
    I2_thick[mask2] = si_thick2
    tau2 = np.zeros(mask2.shape)
    tau2[mask2] = tau2mask2
    full_intensity = full_intensity * np.exp(-tau2) + I2_thick

    I2_temp_thin = ilp.brightness_temp(I2_thin*ilp.specific_int_units, brightparams["nu0"])
    I2_temp_thick = ilp.brightness_temp(I2_thick*ilp.specific_int_units, brightparams["nu0"])

    # I1-------------
    rs1 = rs1[mask1]
    phi1 = phi012[1][mask1]

    si_thin1, si_thick1, tau1mask1, full_profiles1, full_profiles_unit, cosAngReturn1, testa, testb, eta = bright_radial(
        supergrid1, mask1, sign1, spin_case,rs1, isco, thetao, brightparams, funckeys, phi1)

    cosAngReturn1Full = np.zeros(mask1.shape)
    cosAngReturn1Full[mask1] = cosAngReturn1

    I1_thin = np.zeros(mask1.shape)
    I1_thin[mask1] = si_thin1
    I1_thick = np.zeros(mask1.shape)
    I1_thick[mask1] = si_thick1
    tau1 = np.zeros(mask2.shape)
    tau1[mask1] = tau1mask1
    full_intensity = full_intensity * np.exp(-tau1) + I1_thick

    I1_temp_thin = ilp.brightness_temp(I1_thin * ilp.specific_int_units, brightparams["nu0"])
    I1_temp_thick = ilp.brightness_temp(I1_thick*ilp.specific_int_units, brightparams["nu0"])

    # I0-------------
    rs0 = rs0[mask0]
    phi0 = phi012[0][mask0]

    si_thin0, si_thick0, tau0mask0, full_profiles0, full_profiles_unit, cosAngReturn0, testa, testb, eta = bright_radial(
        supergrid0, mask0, sign0, spin_case,rs0, isco, thetao, brightparams, funckeys, phi0)


    cosAngReturn0Full = np.zeros(mask0.shape)
    cosAngReturn0Full[mask0] = cosAngReturn0

    I0_thin = np.zeros(mask0.shape)
    I0_thin[mask0] = si_thin0
    I0_thick = np.zeros(mask0.shape)
    I0_thick[mask0] = si_thick0
    tau0 = np.zeros(mask0.shape)
    tau0[mask0] = tau0mask0
    full_intensity = full_intensity * np.exp(-tau0) + I0_thick

    I0_temp_thin = ilp.brightness_temp(I0_thin * ilp.specific_int_units, brightparams["nu0"])
    I0_temp_thick = ilp.brightness_temp(I0_thick*ilp.specific_int_units, brightparams["nu0"])
    full_temp = ilp.brightness_temp(full_intensity*ilp.specific_int_units, brightparams["nu0"])

    cosAngReturn2Full = cosAngReturn2Full.reshape(N0, N0).T
    cosAngReturn1Full = cosAngReturn1Full.reshape(N0, N0).T
    cosAngReturn0Full = cosAngReturn0Full.reshape(N0, N0).T
    I1_temp_thick = I1_temp_thick.reshape(N0, N0).T
    I2_temp_thick = I2_temp_thick.reshape(N0, N0).T
    I0_temp_thick = I0_temp_thick.reshape(N0, N0).T 
    I0_temp_thin = I0_temp_thin.reshape(N0, N0).T
    I1_temp_thin = I1_temp_thin.reshape(N1, N1).T
    I2_temp_thin = I2_temp_thin.reshape(N2, N2).T
    tau0 = tau0.reshape(N0, N0).T
    tau1 = tau1.reshape(N0, N0).T
    tau2 = tau2.reshape(N0, N0).T
    full_temp = full_temp.reshape(N0, N0).T

    filename = path + ('Intensity_a_{}_i_{}_nu_{}_mass_{}_scaleh_{}_thetab_{}_beta_{}_rie_{}_rb_{}_nth0_{}_te0_{}_'
                       'pdens_{}_ptemp_{}_nscale_{}_emkey_{}_bkey_{}_nkey_{}_tnkey_{}_bnkey_{}.h5').format(
        spin_case,
        i_case,
        "{:.5e}".format(brightparams["nu0"].value),
        "{:.5e}".format(brightparams["mass"].value),
        float(brightparams["scale_height"]),
        "{:.3f}".format(brightparams["theta_b"].value),
        "{:.2f}".format(float(brightparams["beta"])),
        "{:.1f}".format(float(brightparams["r_ie"])),
        "{:.1f}".format(float(brightparams["rb_0"])),
        "{:.1e}".format(brightparams["n_th0"].value),
        "{:.1e}".format(brightparams["t_e0"].value),
        float(brightparams["p_dens"]),
        float(brightparams["p_temp"]),
        "{:.1f}".format(brightparams["nscale"]),
        funckeys["emodelkey"],
        funckeys["bkey"],
        funckeys["nnoisykey"],
        funckeys["tnoisykey"],
        funckeys["bnoisykey"]
    )

    h5f = h5py.File(filename, 'w')

    h5f.create_dataset('bghts0', data=I0_temp_thin)
    h5f.create_dataset('bghts1', data=I1_temp_thin)
    h5f.create_dataset('bghts2', data=I2_temp_thin)
    h5f.create_dataset('bghts2_absorbtion', data=I2_temp_thick)
    h5f.create_dataset('bghts1_absorbtion', data=I1_temp_thick)
    h5f.create_dataset('bghts0_absorbtion', data=I0_temp_thick)
    h5f.create_dataset('bghts_full_absorbtion', data=full_temp)
    h5f.create_dataset('tau2', data=tau2)
    h5f.create_dataset('tau1', data=tau1)
    h5f.create_dataset('tau0', data=tau0)
    h5f.create_dataset('full_profiles2', data=full_profiles2)
    h5f.create_dataset('full_profiles1', data=full_profiles1)
    h5f.create_dataset('full_profiles0', data=full_profiles0)
    h5f.create_dataset('full_profiles_unit', data=full_profiles_unit)

    h5f.create_dataset('cosAngReturn2', data=cosAngReturn2Full)
    h5f.create_dataset('cosAngReturn1', data=cosAngReturn1Full)
    h5f.create_dataset('cosAngReturn0', data=cosAngReturn0Full)

    h5f.create_dataset('testa0', data=testa)
    h5f.create_dataset('testb0', data=testb)
    h5f.create_dataset('eta0', data=eta)

    h5f.close()

    print("File ", filename, " created.")

def br_bv(supergrid0,mask0,N0,rs0,sign0):
    """
    Calculate and save the radial brightness profile
    """
    bghts0 = bright_radial(supergrid0,mask0,sign0,spin_case,rs0,isco,thetao)

    I0 = bghts0.reshape(N0,N0).T

    filename=path+"Intensity_bv_a_%s_i_%s.h5"%(spin_case,i_case)
    h5f = h5py.File(filename, 'w')

    h5f.create_dataset('bghts0', data=I0)

    h5f.close()

    print("File ",filename," created.")

def gfactorf(grid,mask,redshift_sign,a,isco,rs,th,thetao):
    """
    Calculate the redshift factor
    :param grid: alpha and beta grid on the observer plane on which we evaluate the observables
    :param mask: mask out the lensing band, see lb_f.py for detail
    :param redshift_sign: sign of the redshift
    :param a: black hole spin
    :param isco: radius of the inner-most stable circular orbit
    :param rs: source radius
    :param th: source angle, polar coordinate
    :param thetao: observer inclination

    :return: redshift factor at each point.

    """

    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]
    th = th[mask]
    lamb,eta = rt.conserved_quantities(alpha,beta,thetao,a)
    gfact = np.zeros(rs.shape[0])
    redshift_sign = redshift_sign[mask]

    x_aux=rs*np.cos(th)
    y_aux=rs*np.sin(th)

    gfact[rs>=isco]= gDisk(rs[rs>=isco],a,redshift_sign[rs>=isco],lamb[rs>=isco],eta[rs>=isco])
    gfact[rs<isco]= gGas(rs[rs<isco],a,redshift_sign[rs<isco],lamb[rs<isco],eta[rs<isco])

    r_p = 1+np.sqrt(1-a**2)
    gfact[rs<=r_p] = 0

    gs = np.zeros(mask.shape)
    gs[mask] = gfact
    return(gs)

# orbit for the centroid with r=8 and phidot = 0.01
# one may put arbitrary orbit here.
def x0(t):
    return(radhs*np.sin(t*velhs))

def y0(t):
    return(radhs*np.cos(t*velhs))

def flare_model(grid,mask,redshift_sign,a,rs,th,ts,thetao,rwidth,delta_t):

    """
    Calculate the black hole image including the time delay due to lensing and geometric effect
    :param grid: alpha and beta grid on the observer plane on which we evaluate the observables
    :param mask: mask out the lensing band, see lb_f.py for detail
    :param redshift_sign: sign of the redshift
    :param mbar: lensing band index 0,1,2,...
    :param a: black hole spin
    :param isco: radius of the inner-most stable circular orbit
    :param rs: source radius
    :param th: source angle, polar coordinate
    :param ts: time of emission at the source
    :param interpolation: a time series of 2 dimensional brightness function of the source, 3d interpolation object
    :param thetao: observer inclination
    
    :return: image of a lensed equitorial source with only radial dependence. 
    """

    alpha = grid[:,0][mask]
    beta = grid[:,1][mask]
    rs = rs[mask]
    th = th[mask]
    ts = ts[mask]
    lamb,eta = rt.conserved_quantities(alpha,beta,thetao,a)
    brightness = np.zeros(rs.shape[0])
    redshift_sign = redshift_sign[mask]

    x_aux = rs*np.cos(th)
    y_aux = rs*np.sin(th)
    # x0 and y0 is now a function of t, where one can specify an arbitrary equitorial orbit
    brightness = np.exp(-(x_aux-x0(ts+delta_t))**2/rwidth**2-(y_aux-y0(ts+delta_t))**2/rwidth**2)

    brightness[rs>=isco]*= gDisk(rs[rs>=isco],a,redshift_sign[rs>=isco],lamb[rs>=isco],eta[rs>=isco])**gfactor
    brightness[rs<isco]*= gGas(rs[rs<isco],a,redshift_sign[rs<isco],lamb[rs<isco],eta[rs<isco])**gfactor

    r_p = 1+np.sqrt(1-a**2)
    brightness[rs<=r_p] = 0

    I = np.zeros(mask.shape)
    I[mask] = brightness
    return(np.nan_to_num(I))


# def CosAng(r,a,b,lamb,eta):
#     """
#     Calculates the cosine of the emission angle
#     :param r: radius of the source
#     :param a: spin of the black hole
#     :param lamb: angular momentum
#
#     :return: the  cosine of the emission angle
#     """
#     #From eta, solve for Sqrt(p_\theta/p_t)
#     kthkt=np.sqrt(eta+a**2-lamb**2/(np.tan(thetao)**2))
#     #Sqrt(g^{\theta\theta}) Evaluated at the equatorial plane
#     thth=1/r
#     return thth*gDisk(r,a,b,lamb,eta)*kthkt

# def CosAng(r,a,b,lamb,eta):
#     """
#     Calculates the cosine of the emission angle
#     :param r: radius of the source
#     :param a: spin of the black hole
#     :param lamb: angular momentum
#     :return: the  cosine of the emission angle
#     """
#     #From eta, solve for Sqrt(p_\theta/p_t)
#     kthkt=np.sqrt(eta)
#     #Sqrt(g^{\theta\theta}) Evaluated at the equatorial plane
#     thth=1/r
#     return thth*gDisk(r,a,b,lamb,eta)*kthkt

def CosAng(r,a,redshift,lamb,eta):
    """
    Calculates the cosine of the emission angle
    :param r: radius of the source
    :param a: spin of the black hole
    :param lamb: angular momentum

    :return: the  cosine of the emission angle
    """
    #From eta, solve for Sqrt(p_\theta/p_t)

    kthkt=np.sqrt(eta)
    #kthkt=np.sqrt(eta+a**2*np.cos(thetao)**2-lamb**2/(np.tan(thetao)**2))

    #Sqrt(g^{\theta\theta}) Evaluated at the equatorial plane
    thth=1/r

    return thth*redshift*kthkt

