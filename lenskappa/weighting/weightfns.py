import numpy as np

def gal(cat, data):
    """
    Simply counts the number of galaxies
    """
    return len(cat)

def zweight(cat, data):
    z_gal = cat[data['z_gal']]
    z_s = data['z_s']
    return z_s*z_gal - z_gal*z_gal

def mass(cat, data, is_log = True):
    if is_log:
        m_gal_logs = cat[data['m_gal']]
        return 10**(m_gal_logs)
    else:
        return cat[data['m_gal']]

def mass2(cat, data, is_log = True):
    masses = mass(cat, data, is_log)
    return masses**2

def mass3(cat, data, is_log = True):
    masses = mass(cat, data, is_log)
    return masses**3

def oneoverr(cat, data):
    rs = cat[data['r']]
    return 1.0/rs

def zoverr(cat, data):
    zweights = zweight(cat, data)
    rs = cat[data['r']]
    return zweights/rs

def massoverr(cat, data, is_log = True):
    if is_log:
        masses = 10**(cat[data['m_gal']])
    else:
        masses = cat[data['m_gal']]
    rs = cat[data['r']]
    return masses/rs

def mass2overr(cat, data, is_log = True):
    if is_log:
        masses = 10**(cat[data['m_gal']])
    else:
        masses = cat[data['m_gal']]
    rs = cat[data['r']]
    return masses**2/rs

def mass3overr(cat, data, is_log = True):
    if is_log:
        masses = 10**(cat[data['m_gal']])
    else:
        masses = cat[data['m_gal']]
    rs = cat[data['r']]
    return masses**3/rs
def mass2rms(cat, data, is_log = True):
    # This weight can only return a single value. There are no relative weights
    weights = mass2(cat, data, is_log)
    return np.sqrt(sum(weights))
def mass3rms(cat, data, is_log = True):
    # This weight can only return a single value. There are no relative weights
    weights = mass3(cat, data, is_log)
    return np.cbrt(sum(weights))

def mass2overrms(cat, data, is_log = True):
    # This weight can only return a single value. There are no relative weights
    weights = mass2overr(cat, data, is_log = True)
    return np.sqrt(sum(weights))

def mass3overrms(cat, data, is_log = True):
    # This weight can only return a single value. There are no relative weights
    weights = mass3overr(cat, data, is_log = True)
    return np.cbrt(sum(weights))

def flexion(cat, data, is_log = True):
    if is_log:
        m_gals = 10**(cat[data['m_gal']])
    else:
        m_gals = cat[data['m_gal']]
    
    rs = cat[data['r']]
    return m_gals/(rs**3)

def tidal(cat, data, is_log = True):
    if is_log:
        m_gals = 10**(cat[data['m_gal']])
    else:
        m_gals = cat[data['m_gal']]
    
    rs = cat[data['r']]
    return m_gals/(rs**2)

def convergence(cat, data, is_log = True):
    if is_log:
        m_gals = 10**(cat[data['m_gal']])
    else:
        m_gals = cat[data['m_gal']]
    rs = cat[data['r']]
    return np.sqrt(m_gals)/rs

def convergencehalo(cat, data, is_log = True):
    halo_masses = compute_halomass(cat, data, islog=True)
    if is_log:
        m_halos = 10**halo_masses
    rs = cat[data['r']]
    return np.sqrt(m_halos)/rs

def compute_halomass(cat, data, islog=True):
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

    z_gal = cat[data['z_gal']]
    m_gal = cat[data['m_gal']]
    halo_masses = np.zeros(len(m_gal))
    print(z_gal)
    print(m_gal)
    i = 0
    for index, z in z_gal.iteritems():
        print("HI")
        print(z)
        print(m_gal[index])
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
