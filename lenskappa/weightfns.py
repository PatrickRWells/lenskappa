import numpy as np

def gal(cat, data):
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
    if is_log:
        m_halos = 10**(cat[data['m_halo']])
    else:
        m_halos = cat[data['m_halo']]
    rs = cat[data['r']]
    return np.sqrt(m_halos)/rs