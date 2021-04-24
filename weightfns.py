import numpy as np

def zweight(cat, cat_pars, other):
    z_gal = cat[cat_pars['z_gal']]
    z_s = other['z_s']
    return z_s*z_gal - z_gal*z_gal

def mass(cat, cat_pars, other, is_log = True):
    if is_log:
        m_gal_logs = cat[cat_pars['m_gal']]
        return 10**(m_gal_logs)
    else:
        return cat[cat_pars['m_gal']]

def mass2(cat, cat_pars, other, is_log = True):
    masses = mass(cat, cat_pars, other, is_log)
    return masses**2

def mass3(cat, cat_pars, other, is_log = True):
    masses = mass(cat, cat_pars, other, is_log)
    return masses**3

def oneoverr(cat, cat_pars, other):
    rs = cat[cat_pars['r']]
    return 1.0/rs

def zoverr(cat, cat_pars, other):
    zweights = zweight(cat, cat_pars, other)
    rs = cat[cat_pars['r']]
    return zweights/rs

def massoverr(cat, cat_pars, other, is_log = True):
    if is_log:
        masses = 10**(cat[cat_pars['m_gal']])
    else:
        masses = cat[cat_pars['m_gal']]
    rs = cat[cat_pars['r']]
    return masses/rs

def mass2overr(cat, cat_pars, other, is_log = True):
    if is_log:
        masses = 10**(cat[cat_pars['m_gal']])
    else:
        masses = cat[cat_pars['m_gal']]
    rs = cat[cat_pars['r']]
    return masses**2/rs

def mass3overr(cat, cat_pars, other, is_log = True):
    if is_log:
        masses = 10**(cat[cat_pars['m_gal']])
    else:
        masses = cat[cat_pars['m_gal']]
    rs = cat[cat_pars['r']]
    return masses**3/rs
def mass2rms(cat, cat_pars, other, is_log = True):
    # This weight can only return a single value. There are no relative weights
    weights = mass2(cat, cat_pars, other, is_log)
    return np.sqrt(sum(weights))
def mass3rms(cat, cat_pars, other, is_log = True):
    # This weight can only return a single value. There are no relative weights
    weights = mass3(cat, cat_pars, other, is_log)
    return np.cbrt(sum(weights))

def mass2overrms(cat, cat_pars, other, is_log = True):
    # This weight can only return a single value. There are no relative weights
    weights = mass2overr(cat, cat_pars, other, is_log = True)
    return np.sqrt(sum(weights))

def mass3overrms(cat, cat_pars, other, is_log = True):
    # This weight can only return a single value. There are no relative weights
    weights = mass3overr(cat, cat_pars, other, is_log = True)
    return np.cbrt(sum(weights))

def flexion(cat, cat_pars, other, is_log = True):
    if is_log:
        m_gals = 10**(cat[cat_pars['m_gal']])
    else:
        m_gals = cat[cat_pars['m_gal']]
    
    rs = cat[cat_pars['r']]
    return m_gals/(rs**3)

def tidal(cat, cat_pars, other, is_log = True):
    if is_log:
        m_gals = 10**(cat[cat_pars['m_gal']])
    else:
        m_gals = cat[cat_pars['m_gal']]
    
    rs = cat[cat_pars['r']]
    return m_gals/(rs**2)

def convergence(cat, cat_pars, other, is_log = True):
    if is_log:
        m_gals = 10**(cat[cat_pars['m_gal']])
    else:
        m_gals = cat[cat_pars['m_gal']]
    rs = cat[cat_pars['r']]
    return np.sqrt(m_gals)/rs

def convergencehalo(cat, cat_pars, other, is_log = True):
    if is_log:
        m_halos = 10**(cat[cat_pars['m_halo']])
    else:
        m_halos = cat[cat_pars['m_halo']]
    rs = cat[cat_pars['r']]
    return np.sqrt(m_halos)/rs