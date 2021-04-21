import numpy as np

def zweight(cat, zs, z_gal_key):
    z_gal = cat[z_gal_key]
    return zs*z_gal - z_gal*z_gal

def mass(cat, m_gal_key, is_log = True):
    if is_log:
        m_gal_logs = cat[m_gal_key]
        return 10**(m_gal_logs)
    else:
        return cat[m_gal_key]

def mass2(cat, m_gal_key, is_log = True):
    masses = mass(cat, m_gal_key, is_log)
    return masses**2

def mass3(cat, m_gal_key, is_log = True):
    masses = mass(cat, m_gal_key, is_log)
    return masses**3

def oneoverr(cat, r_key):
    rs = cat[r_key]
    return 1.0/rs

def zoverr(cat, zs, z_gal_key, r_key):
    zweights = zweight(cat, zs, z_gal_key)
    rs = cat[r_key]
    return zweights/rs

def massoverr(cat, m_gal_key, r_key, is_log = True):
    if is_log:
        masses = 10**(cat[m_gal_key])
    else:
        masses = cat[m_gal_key]
    rs = cat[r_key]
    return masses/rs

def mass2overr(cat, m_gal_key, r_key, is_log = True):
    if is_log:
        masses = 10**(cat[m_gal_key])
    else:
        masses = cat[m_gal_key]
    rs = cat[r_key]
    return masses**2/rs

def mass3overr(cat, m_gal_key, r_key, is_log = True):
    if is_log:
        masses = 10**(cat[m_gal_key])
    else:
        masses = cat[m_gal_key]
    rs = cat[r_key]
    return masses**3/rs
def mass2rms(cat, m_gal_key, is_log = True):
    # This weight can only return a single value. There are no relative weights
    weights = mass2(cat, m_gal_key, is_log)
    return np.sqrt(sum(weights))
def mass3rms(cat, m_gal_key, is_log = True):
    # This weight can only return a single value. There are no relative weights
    weights = mass3(cat, m_gal_key, is_log)
    return np.cbrt(sum(weights))

def mass2overrms(cat, m_gal_key, r_key, is_log = True):
    # This weight can only return a single value. There are no relative weights
    weights = mass2overr(cat, m_gal_key, r_key, is_log = True)
    return np.sqrt(sum(weights))

def mass3overrms(cat, m_gal_key, r_key, is_log = True):
    # This weight can only return a single value. There are no relative weights
    weights = mass3overr(cat, m_gal_key, r_key, is_log = True)
    return np.cbrt(sum(weights))

def flexion(cat, m_gal_key, r_key, is_log = True):
    if is_log:
        m_gals = 10**(cat[m_gal_key])
    else:
        m_gals = cat[m_gal_key]
    
    rs = cat[r_key]
    return m_gals/(rs**3)

def tidal(cat, m_gal_key, r_key, is_log = True):
    if is_log:
        m_gals = 10**(cat[m_gal_key])
    else:
        m_gals = cat[m_gal_key]
    
    rs = cat[r_key]
    return m_gals/(rs**2)

def convergence(cat, m_gal_key, r_key, is_log = True):
    if is_log:
        m_gals = 10**(cat[m_gal_key])
    else:
        m_gals = cat[m_gal_key]
    rs = cat[r_key]
    return np.sqrt(m_gals)/rs

def convergencehalo(cat, m_halo_key, r_key, is_log = True):
    if is_log:
        m_halos = 10**(cat[m_halo_key])
    else:
        m_halos = cat[m_halo_key]
    rs = cat[r_key]
    return np.sqrt(m_halos)/rs