# Functions for use with atmospheric pwv noise estimates
import numpy as np

def logistic_bandmodel(nuvec, nu0, dnu,a,n):
    '''Returns a logistic-function band model,
       peaking at unity, given input central frequency
       and bandwidth in GHz.

       Note that we are going to replace all values < 1e-5 with 1e-5,
       to keep some divisions from blowing up.

       nuvec:  vector of frequencies at which to return the function
       nu0:  center frequency
       dnu:  bandwidth
       a: prefactor used in logistic model
       n: exponent used in logistic model

       Note:  a and n control the steepness of the edges.

       band = f1*f2, where f1 is a highpass, f2 is a lowpass
       k1 = a*(20/lowedge)**n
       k2 = a*(20/highedge)**m
       f1 =     1/(1+np.exp(-k1*(nu-nu_lowedge )))
       f2 = 1 - 1/(1+np.exp(-k2*(nu-nu_highedge)))
    '''
    nu_lowedge  = nu0-dnu/2
    nu_highedge = nu0+dnu/2
    # I tuned aa and nn by hand to "best fit by eye" (globally)
    # the bands Sara provided in the first plot below.
    #aa = 2     # smaller aa gives broader tails at all frequency bands
    #nn = 0.7   # larger nn gives broader tails at higher frequency bands
    k1 = a*(20/nu_lowedge)**n
    f1 = 1/(1+np.exp(-k1*(nuvec-nu_lowedge )))
    k2 = a*(20/nu_highedge)**n
    f2 = 1-1/(1+np.exp(-k2*(nuvec-nu_highedge)))
    f = f1*f2
    f = f/np.max(f) # normalize to unity
    f = np.where(f<1e-5,1e-5,f)
    return f

def read_bandmodel(fname):
    # Band file must be a tab-delimited file with
    #   column1 = frequency in GHz, 
    #   column2 = transmission, number between 0 and 1.
    # Note: this function can be used to read in either a detector band, or an optics band.
    print('hey')

def alpha_bandmodel(nuvec, nu0, alpha):
    # Intended for use to describe optics transmission.
    # alpha = 2 for a highly overilluminated Lyot stop
    # alpha = 0 for a very underilluminated Lyot stop
    #   transmission = (nuvec/nu0)**alpha
    #   nuvec:  a numpy vector with the frequencies to calculate
    #   nu0: the frequency at which transmission = 1
    #   alpha:  the exponent of the power law model
    transmission = (nuvec/nu0)**alpha
    return transmission

def read_atmospheres(atmos):
    # This reads in the atmospheric models for South Pole and Atacama.
    # These models were created by "am", (insert link to credit to Scott Paine's DOI thing here.)
    # (include these models as text files in this repository)
    # atmos:  a dictionary in which all four atmospheric models will be loaded.
    # Files must be tab delimited with these columns:
    #  Column 1:  frequency in GHz
    #  Column 2:  Tb (Planck brightness)
    #  Column 3:  transmission
    atmos['Pole'] = {}
    atmos['Atacama'] = {}
    atmos['Pole'][300] = np.loadtxt('Pole_300u_50deg.txt',unpack=True)
    atmos['Pole'][400] = np.loadtxt('Pole_400u_50deg.txt',unpack=True)
    atmos['Atacama'][900]  = np.loadtxt('Atacama_900u_50deg.txt',unpack=True)
    atmos['Atacama'][1000] = np.loadtxt('Atacama_1000u_50deg.txt',unpack=True)

def calc_dPdTcmb():
    print('put stuff here')
    
def calc_dPdpwv():
    print('put stuff here')
    
