import pylab as plt
import numpy as np
import sys
from scipy.integrate import quad
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
from scipy import interpolate
Hubble_h = 0.673
Omega_M = 0.315
Omega_Lambda = 0.683
cosmo = FlatLambdaCDM(H0=Hubble_h*100, Om0=Omega_M)


def integrand(y):
    return pow(y,5./3.) / (np.exp(y)-1)

def X_m(m_small):
    return pow(10, 3.6 + (1./4.)*np.log10(m_small))

def Retrive_theta(Temperature):
    T_Tabulated = [1e9,1.5e9,2.0e9,2.5e9,3.0e9,3.5e9,4.0e9,4.5e9,5.0e9,5.5e9,6.0e9,6.5e9,7.0e9,7.5e9,8.0e9,8.5e9,9.0e9,9.5e9,10.0e9]
    Theta_Tabl  = [0.1686, 0.2530,0.3373,0.4216,0.5059,0.5902,0.6746,0.7589,0.8432,0.9275,1.0118,1.0961,1.1805,1.2648,1.3491,1.4334,1.5177,1.6021,1.6864]
    funcT = interp1d(T_Tabulated,Theta_Tabl,fill_value='extrapolate')
    #print(len(T_Tabulated), len(Theta_Tabl))
    #print(Temperature, func(10e9))
    #sys.exit()
    return funcT(Temperature)#np.interp(Temperature,T_Tabulated,Theta_Tabl)

def Retrive_g_theta(Temperature):
    T_Tabulated   = [1e9,1.5e9,2.0e9,2.5e9,3.0e9,3.5e9,4.0e9,4.5e9,5.0e9,5.5e9,6.0e9,6.5e9,7.0e9,7.5e9,8.0e9,8.5e9,9.0e9,9.5e9,10.0e9]
    G_Theta_Tabl  = [12.003, 6.7292, 4.5134, 3.3386, 2.6261, 2.1540, 1.8209, 1.5746, 1.3859, 1.2369, 1.1166, 1.0175, 0.9345, 0.8640, 0.8035, 0.7509, 0.7048, 0.6641, 0.6278]
    funcG = interp1d(T_Tabulated,G_Theta_Tabl,fill_value='extrapolate')
    return funcG(Temperature)

def Retrive_theta_3_K(Temperature):
    T_Tabulated   = [1e9,1.5e9,2.0e9,2.5e9,3.0e9,3.5e9,4.0e9,4.5e9,5.0e9,5.5e9,6.0e9,6.5e9,7.0e9,7.5e9,8.0e9,8.5e9,9.0e9,9.5e9,10.0e9]
    K_Theta_Tabl  = [8.783e-6, 2.982e-4,2.472e-3,1.092e-2,3.408e-2,8.550e-2,1.849e-1,3.593e-1,6.438e-1,1.083,1.731,2.654,3.930,5.650,7.922,1.086e1,1.463e1,1.933e1,2.519e1]
    funcKT = interp1d(T_Tabulated,K_Theta_Tabl,fill_value='extrapolate')
    return funcKT(Temperature)

def Retrive_F_theta(theta):
    if(theta <= 1.0): F_theta = 4.0*pow(2 * theta / np.pi*np.pi*np.pi ,1./2.) * (1 + 1.781*pow(theta,1.34)) + 1.73*pow(theta,3./2.)*(1 + ( 1.1*theta) +pow(theta,2.0) - (1.25*pow(theta,5./2.)) )
    else: F_theta = ( (9*theta/(2.*np.pi)) * (np.log((1.123*theta) + 0.48) + 1.5 )) + 2.30*theta*(np.log(1.123*theta) + 1.28)
    return F_theta


def Brems_lum(Temperature, frequency, rM, rm, MB, m_d, alpha, c1, h_Planck, k_Boltzman):
    theta = Retrive_theta(Temperature)
    F_theta = Retrive_F_theta(theta)
    return 2.29*1e24 * pow(alpha,-2.0) * pow(c1,-2) * np.log(rM/rm)* F_theta * pow(Temperature,-1) * MB * pow(m_d,2.0) * np.exp(-1.0* h_Planck * frequency / (k_Boltzman*Temperature) )

def Sync_lum(Temperature, frequency, m_d, MB, s1, s2, s3):
    return s3 * pow(s1*s2,8./5.) * pow(MB,6./5.) * pow(m_d,4./5.) * pow(Temperature,21./5.) * pow(frequency,2./5.)


def Q_e_plus(alpha, c1, c3, beta, mbh, mdot, delta, ff, rmin, Temp):
    #print("alpha", alpha, "c1", c1, "c3", c3, "beta", beta, "mbh", mbh, "mdot", mdot, "delta", delta, "ff", ff, "rmin", rmin, "Temp", Temp)
    Term1 = 1.2e38  * Retrive_g_theta(Temp)*pow(alpha,-2)*pow(c1,-2) * c3 * beta * mbh * pow(mdot,2) * pow(rmin,-1)
    Term2 = 9.39e38 * delta * ( (1.-beta) / ff ) * c3 * mbh * mdot * pow(rmin,-1)
    return (Term1 + Term2)/1e38

def Power_Syn(alpha, beta, c1, c3, rmin, Temp, mbh, mdot):
    Term = 5.3e35*pow(X_m(mdot)/1000,3)*pow(alpha/0.3,-3./2.)*pow((1.-beta)/0.5,3./2.)*pow(c1/0.5,-3./2.)*pow(c3/0.3,3./2.)*pow(rmin/3,-7./4.)*pow(Temp/1e9,7)*pow(mbh,1./2.)*pow(mdot,3./2.)
    return Term / 1e38

def Power_Brems(alpha, c1, rmax, rmin, Temp, mbh, mdot):
    theta = Retrive_theta(Temp)
    Term  = 4.78e34*pow(alpha, -2)*pow(c1,-2)*np.log(rmax/rmin)*Retrive_F_theta(theta)*mbh*pow(mdot,2)
    return Term / 1e38

def Power_Compt(s1, s2, s3, Temp, mdot, mbh, alpha, c1, rmin):
    nu_peak  = s1*s2*pow(mbh,-1./2.) * pow(mdot,1./2.) * pow(Temp,2.0) * pow(rmin,-5./4.)
    Lnu_peak = Sync_lum(Temp, nu_peak, mdot, mbh, s1, s2, s3)
    theta = Retrive_theta(Temp)
    A = 1 + (4.*theta) + (16*theta*theta)
    tau_es = (23.87*mdot) * pow(alpha/0.3,-1) * pow(c1/0.5,-1) * pow(rmin/3,-1./2.)
    alpha_c = -np.log(tau_es)/np.log(A) 
    #if(alpha_c > 1.): return nu_peak*Lnu_peak/(alpha_c-1)
    #else: return (nu_peak*Lnu_peak/(1-alpha_c)) * pow(6.2e7*Temp*1e-9/(nu_peak*1e-12),1-alpha_c) 
    Term = (nu_peak*Lnu_peak/(1-alpha_c)) * ( pow(6.2e7*Temp*1e-9/(nu_peak*1e-12),1-alpha_c) -1)
    return Term / 1e38


def Temperature_Equilibrium(TT, s1, s2, s3, mdot, mbh, alpha, beta, c1, c3, rmin, rmax, ff, delta):
    Q_e_p       = Q_e_plus(alpha, c1, c3, beta, mbh, mdot, delta, ff, rmin, TT)
    Psyn        = Power_Syn(alpha, beta, c1, c3, rmin, TT, mbh, mdot)
    Pbrem       = Power_Brems(alpha, c1, rmax, rmin, TT, mbh, mdot)
    Pcomp       = Power_Compt(s1, s2, s3, TT, mdot, mbh, alpha, c1, rmin)
    Diff        = (Psyn + Pbrem + Pcomp - Q_e_p)
    #print("Q_e_p", Q_e_p, "Psyn", Psyn, "Pbrem", Pbrem, "Pcomp", Pcomp, "Temp", TT, "Diff", Diff)
    return Diff

def compute_L_nu(MBH, M_dot, Lbol, fmin, fmax):
    G = 6.67e-11 # m3 kg-1 s-2
    c = 3e8 # m/s
    Sigma_Boltzman = 5.67e-8  # kg s-2 K-4
    k_Boltzman     = 1.38e-23 # J K-1  --> kg m2 s-1 K-1
    h_Planck       = 6.63e-34 # J Hz-1 --> kg m2 s-1 Hz-1
    M_Edd   =        2.2 * (MBH/1e8)   # Msun/yr
    m_dot   = M_dot/M_Edd
    if(m_dot>0.01): frec, L_nu = compute_ThinDisc_L_nu(MBH, M_dot, Lbol,  G , c, Sigma_Boltzman, k_Boltzman, h_Planck, fmin, fmax)
    else: frec, L_nu = compute_ADAF_L_nu(m_dot, MBH , k_Boltzman, h_Planck, fmin, fmax)
    return frec, L_nu 
    

def compute_ThinDisc_L_nu(MBH, M_dot, Lbol, G , c, Sigma_Boltzman, k_Boltzman, h_Planck, fmin, fmax):
    #print ("Thin disc regimen")
    From_Msun_To_Kg = 5.02e-31
    From_yr_To_sec  = 3.154e7
    eV_to_Hz        = 2.41e14 
    MBH   = MBH /  From_Msun_To_Kg
    M_dot = M_dot / (From_Msun_To_Kg * From_yr_To_sec)
    Rs      = 2.0 * ( G * MBH )  / pow(c,2.0) # m
    A = pow(3 * M_dot * pow(pow(c,3)/(G*MBH),2.0) / (64.0 * np.pi * Sigma_Boltzman),1./4.)
    Constant = (32./3.) * pow(np.pi*Rs,2.0) * (h_Planck/pow(c,2.0)) * pow(k_Boltzman*A/h_Planck,8./3)

    frec = np.logspace(np.log10(fmin),np.log10(fmax),100)
    L_nu = np.empty(len(frec))

    xmin = 3
    xmax = 3000

    # To compute the (Hard) X-ray corona (Cocchiararo et al. 2024 & Regan 2019)
    #Lbol = epsilon * M_dot * pow(c,2.0) # kg m2 s-3
    Lsun = 3.82e26 # kg m2 s-3
    # Marconi transformation:
    LL = np.log10(Lbol/Lsun) - 12
    L_X_rays = pow(10,np.log10(Lbol) - 1.54 - (0.24*LL) - (0.012*pow(LL, 2.0)) + (0.0015*pow(LL, 3.0)))
    Estart   = 200 #eV
    Eend     = 10000 # eV
    Ec       = 300000 # ev # Cut-off Shen et al. 2020,  Dadina 2008; Ueda et al. 2014; Aird et al. 2015a
    Estart   = Estart * eV_to_Hz 
    Eend     = Eend   * eV_to_Hz 
    Ec       = Ec     * eV_to_Hz 
    Slope_Xrays = -1.7
    BHnorm = (1+Slope_Xrays) * L_X_rays/ (pow(Eend,Slope_Xrays+1) - pow(Estart,Slope_Xrays+1) )

    rmax = 49./36.
    Tmax = 3.*M_dot*pow(pow(c,3)/(G*MBH),2.0)/(64.0 * np.pi * Sigma_Boltzman * pow(rmax,3.0))
    Tmax = 3.0 * pow(Tmax,1./4.)
    Emin_Xrays = k_Boltzman * Tmax * 6.242e18 # ev
    Emin_Xrays = Emin_Xrays * eV_to_Hz # Hz


    for i in np.arange(0,len(frec),1):
        ymin = ( (h_Planck/k_Boltzman) *  frec[i]  / A ) * pow(xmin,3./4.)
        ymax = ( (h_Planck/k_Boltzman) *  frec[i]  / A ) * pow(xmax,3./4.)
        I, err = quad(integrand, ymin, ymax)
        L_nu[i] = Constant * I * pow(frec[i],1./3.) * 1e7 # erg s-1 Hz-1
        #if(frec[i]>=0.1*Estart):# and frec[i]<=Eend):
        if(frec[i]>=Emin_Xrays):# and frec[i]<=Eend):
            L_nu[i] +=  BHnorm * 1e7 * pow(frec[i],Slope_Xrays) * np.exp(-frec[i] / Ec )# erg s-1 Hz-1
    return frec, L_nu



def compute_ADAF_L_nu(m_dot, MBH ,k_Boltzman,h_Planck, fmin, fmax):
    alpha   = 0.3
    delta   = 1./2000
    beta    = 0.5
    rmin    = 3.0
    rmax    = 1e3
    ff      = 0.05

    c1    = 0.5
    c3    = 0.3
    s1    = 1.42e9    * pow(alpha,-1./2.) * pow(1-beta,1./2.) * pow(c1,-1./2.) * pow(c3,1./2)
    s2    = 1.19e-13 * X_m(m_dot)
    s3    = 1.05e-24
    data = (s1, s2, s3, m_dot, MBH, alpha, beta, c1, c3, rmin, rmax, ff, delta)
    if(m_dot<1e-8): Tguess = 5e11
    else: Tguess = 5e10
    root =fsolve(Temperature_Equilibrium, Tguess, args=data, xtol=1e-4, maxfev=1000)
    Te = root[0]
    #print (m_dot, Te)
    #print ("Temperature [1e9 K]", Te/1e9, "Individual cooling process", Temperature_Equilibrium(Te,s1, s2, s3, m_dot, MBH, alpha, beta, c1, c3, rmin, rmax, ff, delta))
    #print ("Temperature", Te)

    th = Retrive_theta(Te)
    A = 1 + (4.*th) + (16*th*th)
    tau_es = (23.87*m_dot) * pow(alpha/0.3,-1) * pow(c1/0.5,-1) * pow(rmin/3,-1./2.)
    alpha_c = -np.log(tau_es)/np.log(A)
    #print (m_dot, 1-alpha_c)
    nu_p = s1*s2*pow(MBH,-1./2.) * pow(m_dot,1./2.) * pow(Te,2.0) * pow(rmin,-5./4.) # Hz

    frec = np.logspace(np.log10(fmin),np.log10(fmax),100)
    L_nu = np.empty(len(frec))
    L_nu.fill(0.0)
    #print ("nu_p [1e11 Hz]",nu_p/1e11)


    for i in np.arange(0,len(frec),1):
        L_nu[i] += Brems_lum(Te, frec[i], rmax, rmin, MBH, m_dot, alpha, c1, h_Planck, k_Boltzman)
        if(frec[i]<1*nu_p):
            L_nu[i] += Sync_lum(Te, frec[i],m_dot,MBH,s1,s2,s3)
        else:
            if(frec[i] < 3*k_Boltzman*Te/h_Planck): L_nu[i] += Sync_lum(Te, nu_p, m_dot,MBH,s1,s2,s3) * pow(frec[i]/nu_p,-alpha_c)


    return frec, L_nu


def lambdaPivot(longi,trans):
    LS = (trans*longi)
    S_L= (trans/longi)
    lPivot=(np.trapz(LS)/np.trapz(S_L))**0.5
    return lPivot

def Convolve_L_nu_with_Filters(frec,L_nu,redshift,FilterWavelength,FilterTransmision):


    dl = cosmo.luminosity_distance(redshift).value # Mpc
    dl = dl * 3.086e24 # cm
    c = 3e8 # m/s
    c = 3e8 * 1e10 # A/s
    f_nu = (L_nu/pow(dl,2.0)) / (4.0*np.pi)
    
    frec = frec * ( 1 + redshift )
    lamb = c / frec # A
    f_lambda = (f_nu / pow(lamb,2.0)) * c
    func_interp = interpolate.interp1d(np.log10(lamb),np.log10(f_lambda))
    f_lambda_interpolated = func_interp(np.log10(FilterWavelength))
    f_lambda_interpolated = pow(10,f_lambda_interpolated)

    Numerator   = FilterTransmision*FilterWavelength*f_lambda_interpolated
    Denominator = c * FilterTransmision/FilterWavelength
    f_nu_avg    = (np.trapz(Numerator)/np.trapz(Denominator))
    magnitude = -2.5*np.log10(f_nu_avg) - 48.6
    lambda_pivot = lambdaPivot(FilterWavelength,FilterTransmision)
    f_lamb_avg = f_nu_avg * c / lambda_pivot #-2.5*np.log10(f_nu_avg) - 48.6

    return lambda_pivot, magnitude, f_lamb_avg


if __name__ == '__main__':

    c_light = 3e8 # m/s
    Msun_To_Kg = 5.02e-31
    yr_To_sec  = 3.154e7


    #fedd    = pow(10,np.arange(-10,-2.5,0.25))
    #fedd    = [ 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1.0]
    fedd    = [3e-4, 6e-4, 12e-4, 24e-4 ] #0.001
    redshift = 0.1
    MBH     = 1e9 # Msun
    epsilon = 0.01
    scolor = plt.cm.rainbow_r(np.linspace(0, 1.0, len(fedd)))

    for i in np.arange(0,len(fedd),1):
        M_dot   = fedd[i] * 2.2 * (MBH/1e8)   # Msun/yr
        M_dot_cgs = M_dot / (Msun_To_Kg * yr_To_sec)
        Lbol    = epsilon * M_dot_cgs * pow(c_light,2.0) # kg m2 s-3
        #M_dot   = fedd * 1.26e38 * 1e-7 * MBH  / ( epsilon *  pow(c_light,2.0) )  # kg  s-1
        #M_dot = M_dot * 5.02785e-31 / 3.17098e-8 # Msun yr-1
        fmin = 1e7
        fmax = 1e22

        # Computation of the spectrum
        f, L_v = compute_L_nu(MBH, M_dot, Lbol, fmin, fmax)

        #plt.axhline(47,color = "blue")
        plt.plot(np.log10(f),np.log10(L_v*f), color = scolor[i], marker = "o", markersize = 1.0, linestyle = "-")

    yticks = np.arange(30,48,1)
    plt.yticks(yticks)
    xticks = np.arange(5,22,1)
    plt.xticks(xticks)
    plt.xlim(7.5,21.5)
    plt.ylim(30,48)
    plt.xlabel(r"$ \rm log_{10}(\nu \, /Hz)$", fontsize= 26)
    plt.ylabel(r"$ \rm log_{10} ( \nu L_{\nu} \, /erg \, s^{-1})$", fontsize= 20)
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    sys.exit()

    filterNames = ["LSST_LSST.u.dat", "LSST_LSST.g.dat", "LSST_LSST.r.dat", "LSST_LSST.i.dat", "LSST_LSST.z.dat", "LSST_LSST.y.dat"]
    color = ["blue", "green", "red", "orange", "red", "black"]
    for i in np.arange(0,len(filterNames),1):
        f_filter, S_filter = np.loadtxt("/home/david/Desktop/SED_AGNs/LSST_filters/" + filterNames[i], unpack = True)
        lam, mag, f_lam_avg = Convolve_L_nu_with_Filters(f, L_v, redshift, f_filter, S_filter) 
        #plt.plot(lam, f_lam_avg, color = color[i], marker = "o", linestyle = " ")
        plt.plot(lam, mag, color = color[i], marker = "o", linestyle = " ")
    plt.show()
