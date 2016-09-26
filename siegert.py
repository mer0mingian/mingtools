"""siegert.py: Function calculating the firing rates of leaky
integrate-and-fire neurons given their parameter and mean and variance
of the input. Rates rates for delta shaped PSCs after Brunel & Hakim 1999. 
Rate of neuron with synaptic filtering with time constant tau_s after
Fourcoud & Brunel 2002.

Authors: Moritz Helias, Jannis Schuecker, Hannah Bos

Adapted for a generalized integrate-and-fire neuron with two variables
as in Richardson et al. (2003) and Brunel et al. (2003) by Daniel Mingers.
Essentially, just redefined the integration borders and the mu.
"""

from scipy.integrate import quad
from scipy.special import erf
from scipy.special import zetac
import numpy as np
import math

"""
Variables used in this module:
tau_m: membrane time constant
tau_r: refractory time constant
tau_s: synaptic time constant
tau_N: time constant of the GWN
V_th: threshold
V_r: reset potential
mu: mean input / average membrane potential
sigma: std of equivalent GWN input
"""

def mu_ (I_0, g_1, I_N, V_th, V_r, C):
    """
    calculates average voltage for the siegert for a gif2 neuron.
    see eq. (52) in Brunel et al. (2003)
    """
    dis = V_th + V_r
    dis2 = dis / 2. * g_1

    h1 = (I_0 - dis2) * (I_0 - dis2)
    h2 = 2. * g_1 * I_N / C
    h3 = I_0 + dis2
    result = (h3 - np.sqrt(h1 +h2))/2. / g_1
    return result


def nu0_gif2(I_0, I_N, tau_m, tau_s, tau_r, tau_N=1., g=25., g_1=25., V_th=20., V_r=14.):
    """Calculates stationary firing rates for exponential PSCs using
    expression with taylor expansion in k = sqrt(tau_s/tau_m) (Eq. 433
    in Fourcoud & Brunel 2002) for the generalized integrate-and-fire
    model in Richardson et al. (2003), see eq. (A10).
    This only holds for tau_1 >> tau, effectively with a factor > 5.
    """
    # tau_1 = 100.

    C = tau_m / g
    mu = mu_(I_0, g_1, I_N, V_th, V_r, C)
    gamma = g_1 / g
    sigma = I_N / g * np.sqrt(tau_N / tau_m)  # (A9)
    # instead use equation 10 from Richardson:
    # sigma = I_N * np.sqrt(tau_N*(C+tau_1*g +g_1*tau_1)/2./C/(g+g_1)/(g*tau_1+C))

    alpha = np.sqrt(2.) * abs(zetac(0.5) + 1)
    x_th = np.sqrt(2.) * (V_th - gamma * mu - I_0/g) / sigma
    x_r = np.sqrt(2.) * (V_r - gamma * mu - I_0/g) / sigma

    # preventing overflow in np.exponent in Phi(s)
    if x_th > 20.0 / np.sqrt(2.):
        result = nu_0(tau_m, tau_r, V_th, V_r, mu, sigma)
    else:
        r = nu_0(tau_m, tau_r, V_th, V_r, mu, sigma)
        dPhi = Phi(x_th) - Phi(x_r)
        result = r - np.sqrt(tau_s / tau_m) * alpha / \
                     (tau_m * np.sqrt(2)) * dPhi * (r * tau_m) ** 2
    if math.isnan(result):
        print mu, sigma, x_th, x_r
    return result


def nu_0(tau_m, tau_r, V_th, V_r, mu, sigma):
    """ Calculates stationary firing rates for delta shaped PSCs."""

    if mu <= V_th + (0.95 * abs(V_th) - abs(V_th)):
        return siegert1(tau_m, tau_r, V_th, V_r, mu, sigma)
    else:
        return siegert2(tau_m, tau_r, V_th, V_r, mu, sigma)


def nu0_fb(tau_m, tau_s, tau_r, V_th, V_r, mu, sigma):
    
    alpha = np.sqrt(2)*abs(zetac(0.5)+1)
    # effective threshold    
    V_th1 = V_th + sigma*alpha/2.*np.sqrt(tau_s/tau_m)
    # effective reset    
    V_r1 = V_r + sigma*alpha/2.*np.sqrt(tau_s/tau_m)
    # use standard Siegert with modified threshold and reset
    return nu_0(tau_m, tau_r, V_th1, V_r1, mu, sigma)


# stationary firing rate of neuron with synaptic low-pass filter
# of time constant tau_s driven by Gaussian noise with mean mu and 
# standard deviation sigma, from Fourcaud & Brunel 2002
def nu0_fb433(tau_m, tau_s, tau_r, V_th, V_r, mu, sigma):
    """Calculates stationary firing rates for exponential PSCs using
    expression with taylor expansion in k = sqrt(tau_s/tau_m) (Eq. 433
    in Fourcoud & Brunel 2002)
    """
    
    alpha = np.sqrt(2.) * abs(zetac(0.5) + 1)
    x_th = np.sqrt(2.) * (V_th - mu) / sigma
    x_r = np.sqrt(2.) * (V_r - mu) / sigma

    # preventing overflow in np.exponent in Phi(s)
    if x_th > 20.0 / np.sqrt(2.):
        result = nu_0(tau_m, tau_r, V_th, V_r, mu, sigma)
    else:
        r = nu_0(tau_m, tau_r, V_th, V_r, mu, sigma)
        dPhi = Phi(x_th) - Phi(x_r)
        result = r - np.sqrt(tau_s / tau_m) * alpha / \
            (tau_m * np.sqrt(2)) * dPhi * (r * tau_m)**2
    if math.isnan(result):
        print mu, sigma, x_th, x_r
    return result


def Phi(s):
    return np.sqrt(np.pi / 2.) * (np.exp(s**2 / 2.) * (1 + erf(s / np.sqrt(2))))


def Phi_prime_mu(s, sigma):
    return -np.sqrt(np.pi) / sigma * (s * np.exp(s**2 / 2.) * (1 + erf(s / np.sqrt(2)))
                                      + np.sqrt(2) / np.sqrt(np.pi))


def siegert1(tau_m, tau_r, V_th, V_r, mu, sigma):
    # for mu < V_th
    y_th = (V_th - mu) / sigma
    y_r = (V_r - mu) / sigma

    def integrand(u):
        if u == 0:
            return np.exp(-y_th**2) * 2 * (y_th - y_r)
        else:
            return np.exp(-(u - y_th)**2) * (1.0 - np.exp(2 * (y_r - y_th) * u)) / u

    lower_bound = y_th
    err_dn = 1.0
    while err_dn > 1e-12 and lower_bound > 1e-16:
        err_dn = integrand(lower_bound)
        if err_dn > 1e-12:
            lower_bound /= 2

    upper_bound = y_th
    err_up = 1.0
    while err_up > 1e-12:
        err_up = integrand(upper_bound)
        if err_up > 1e-12:
            upper_bound *= 2

    # check preventing overflow
    if y_th >= 20:
        out = 0.
    if y_th < 20:
        out = 1.0 / (tau_r + np.exp(y_th**2)
                     * quad(integrand, lower_bound, upper_bound)[0] * tau_m)

    return out

def siegert2(tau_m, tau_r, V_th, V_r, mu, sigma):
    # for mu > V_th
    y_th = (V_th - mu) / sigma
    y_r = (V_r - mu) / sigma

    def integrand(u):
        if u == 0:
            return 2 * (y_th - y_r)
        else:
            return (np.exp(2 * y_th * u - u**2) - np.exp(2 * y_r * u - u**2)) / u

    upper_bound = 1.0
    err = 1.0
    while err > 1e-12:
        err = integrand(upper_bound)
        upper_bound *= 2

    return 1.0 / (tau_r + quad(integrand, 0.0, upper_bound)[0] * tau_m)

def d_nu_d_mu_fb433(tau_m, tau_s, tau_r, V_th, V_r, mu, sigma):
    alpha = np.sqrt(2) * abs(zetac(0.5) + 1)
    x_th = np.sqrt(2) * (V_th - mu) / sigma
    x_r = np.sqrt(2) * (V_r - mu) / sigma
    integral = 1. / (nu_0(tau_m, tau_r, V_th, V_r, mu, sigma) * tau_m)
    prefactor = np.sqrt(tau_s / tau_m) * alpha / (tau_m * np.sqrt(2))
    dnudmu = d_nu_d_mu(tau_m, tau_r, V_th, V_r, mu, sigma)
    dPhi_prime = Phi_prime_mu(x_th, sigma) - Phi_prime_mu(x_r, sigma)
    dPhi = Phi(x_th) - Phi(x_r)
    phi = dPhi_prime * integral + (2 * np.sqrt(2) / sigma) * dPhi**2
    return dnudmu - prefactor * phi / integral**3

def d_nu_d_mu(tau_m, tau_r, V_th, V_r, mu, sigma):
    y_th = (V_th - mu)/sigma
    y_r = (V_r - mu)/sigma
    nu0 = nu_0(tau_m, tau_r, V_th, V_r, mu, sigma)
    return np.sqrt(np.pi) * tau_m * nu0**2 / sigma * (np.exp(y_th**2) * (1 + erf(y_th)) - np.exp(y_r**2) * (1 + erf(y_r)))
