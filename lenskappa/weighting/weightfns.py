import numpy as np

def gal(cat, **kwargs):
    """
    Simply counts the number of galaxies
    """
    return np.ones(len(cat), dtype=np.float64)

def zweight(cat, z_s, *args, **kwargs):
    z_gal = cat['z_gal']
    return z_s*z_gal - z_gal*z_gal

def mass(cat, **kwargs):
    return cat['m_gal']

def mass2(cat, **kwargs):
    masses = mass(cat)
    return masses**2

def mass3(cat, **kwargs):
    masses = mass(cat)
    return masses**3

def oneoverr(cat, **kwargs):
    rs = cat['r'].value
    return 1.0/rs

def zoverr(cat, **kwargs):
    zweights = zweight(cat)
    rs = cat['r'].value
    return zweights/rs

def massoverr(cat, **kwargs):
    masses = cat['m_gal']
    rs = cat['r'].value
    return masses/rs

def mass2overr(cat, **kwargs):
    masses = cat['m_gal']
    rs = cat['r'].value
    return masses**2/rs

def mass3overr(cat, **kwargs):
    masses = cat['m_gal']
    rs = cat['r'].value
    return masses**3/rs
def mass2rms(cat, **kwargs):
    # This weight can only return a single value. There are no relative weights
    weights = mass2(cat)
    return weights
def mass3rms(cat, **kwargs):
    # This weight can only return a single value. There are no relative weights
    weights = mass3(cat)
    return weights

def mass2overrms(cat, **kwargs):
    # This weight can only return a single value. There are no relative weights
    weights = mass2overr(cat)
    return weights

def mass3overrms(cat, **kwargs):
    # This weight can only return a single value. There are no relative weights
    weights = mass3overr(cat)
    return weights

def flexion(cat, **kwargs):
    m_gals = cat['m_gal']
    rs = cat['r'].value
    return m_gals/(rs**3)

def tidal(cat, **kwargs):
    m_gals = cat['m_gal']

    rs = cat['r'].value
    return m_gals/(rs**2)

def convergence(cat, **kwargs):
    m_gals = cat['m_gal']
    rs = cat['r'].value
    return np.sqrt(m_gals)/rs

def convergencehalo(cat, **kwargs):
    halo_masses = compute_halomass(cat)
    rs = cat['r'].value
    return np.sqrt(halo_masses)/rs

def compute_halomass(cat, **kwargs):
    # z <= 1
    M10_ = 12.35
    M1a_ = 0.28
    Ms00_ = 10.72
    Ms0a_ = 0.55
    b0_ = 0.44
    ba_ = 0.18
    d0_ = 0.57
    da_ = 0.17
    g0_ = 1.56
    ga_ = 2.51
    # z >= 1:
    M10 = 12.27
    M1a = -0.84
    Ms00 = 11.09
    Ms0a = 0.56
    b0 = 0.65
    ba = 0.31
    d0 = 0.55
    da = -0.12
    g0 = 1.12
    ga = -0.53

    z_gal = cat['z_gal']
    m_gal = cat['m_gal']
    halo_masses = np.zeros(len(m_gal))
    i = 0
    for index, z in enumerate(z_gal, **kwargs):
        if z <= 1:
            a = 1/(1+z)

            logM1a = M10_ + M1a_ * (a - 1)
            logMs0a = Ms00_ + Ms0a_ * (a-1)
            notlogMs0a = 10 ** logMs0a
            b = b0_ + ba_ * (a-1)
            d = d0_ + da_ * (a-1)
            g = g0_ + ga_ * (a-1)
            halo_mass = logM1a + b * m_gal[index] + ((10 ** m_gal[index]/notlogMs0a)**d)/(1+(10 ** m_gal[index]/notlogMs0a)**(-g)) - 1/2
        else:
            a = 1/(1+z)
            logM1a = M10 + M1a * (a-1)
            logMs0a = Ms00 + Ms0a * (a-1)
            notlogMs0a = 10 ** logMs0a
            b = b0 + ba * (a-1)
            d = d0 + da * (a-1)
            g = g0 + ga * (a-1)
            halo_mass = logM1a + b * (m_gal[index] - logMs0a) + ((10 ** m_gal[index]/notlogMs0a)**d)/(1+(10 ** m_gal[index]/notlogMs0a)**(-g)) - 1/2
        halo_masses[i] = halo_mass
        i+=1
    return halo_masses
