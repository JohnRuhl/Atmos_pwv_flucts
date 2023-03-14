# Functions for use with atmospheric pwv noise estimates
import numpy as np

def dBdT(T, freq):
    x= 2*(h**2)*(freq**4)
    y= np.exp((h*freq)/(k*T)) 
    z= 1/(k*((T)**2)*(c**2)) 
    w= (1/((-1)+np.exp((h*freq)/(k*T)))**2)
    dB_dT= x*y*z*w
    return dB_dT
    
def bnu_aomega(freq, t):
    a = (2*h*freq**3)/(c**2)
    b = 1/(np.exp(((h*freq)/(k*t)))-1)
    B_v= a*b
    a_omega= (c**2/freq**2)
    b_v= B_v*a_omega
    return b_v

#constants
c= 3e+8 #m/s
h= 6.62607015e-34 #Plank constant in J/Hz
k= 1.380649e-23 #Boltzmann constant in J/K

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

def flat_band(nuvec, band_int, nu0,dnu):
    '''Returns a flat band model given inputs
    
      nuvec:  vector of frequencies at which to return the function
      band_int: the desired integral of the band.
      nu0:  center frequency
      dnu:  bandwidth
      
    '''
    low = nu0-dnu/2
    high = nu0+dnu/2

    f = np.where((nuvec>low)&(nuvec<high), 1, 0)
    band_integral = np.trapz(f,x=nuvec)
    f = (band_int/band_integral)*f
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
    atmos['Pole'][300] = np.loadtxt('SPole_300um_wojacobian.txt',unpack=True)
    atmos['Pole'][400] = np.loadtxt('SPole_400um_wojacobian.txt',unpack=True)
    atmos['Atacama'][900]  = np.loadtxt('Atacama_900um_wojacobian.txt',unpack=True)
    atmos['Atacama'][1000] = np.loadtxt('Atacama_1000um_wojacobian.txt',unpack=True)

def calc_dPdTcmb(nu_ghz, bandmodel, atmos_trans):
    # Calculates dP_optical/dT_cmb for a single-moded, two-polarization detector, 
    # given the input instrument band (defined by nuvec and bandmodel) and atmospheric transmission.
    #
    # Inputs: (all numpy arrays)
    #   nuvec:  the frequencies to be integrated over, in GHz
    #   bandmodel:  the instrument band, as a function of frequencies in nuvec
    #   atmos_trans:  the atmospheric transmission, as a function of frequencies in nuvec.
    #
    # Output: 
    #   dPdT_cmb, in Watts/K_cmb, at the detector.
    nu= nu_ghz*1e9
    dB_dT = dBdT(2.7, nu)
    AOmega= (c**2/nu**2)
    dPdTcmb=  np.trapz(bandmodel*atmos_trans*dB_dT*AOmega, nu) 
    return dPdTcmb
    
    #need to adjust by integrating Tb_atm over f_inst to get dPopt_atm/dpwv
    #w/ f_tot=f_inst * f_atm and f_inst= f_detect * f_lyot
    
def calc_dPdpwv(nu_ghz, tb1, tb2, dpwv, inst_band):
    # Calculates dP_optical/d_pwv for a single-moded, two-polarization detector, 
    # given the input instrument band (defined by nuvec and bandmodel) and atmospheric Tb(nu) for two neighboring pwvs.
    #
    # Inputs: (all numpy arrays)
    #   nughz:  the frequencies to be integrated over, in GHz
    #   tb1, tb2: atmospheric Tb's as a function of frequencies in nuvec
    #   dpwv: difference in pwv (in units of mm) at which those tb's were calculated.
    #   bandmodel:  the instrument band, as a function of frequencies in nuvec
    #
    # Output:
    #   dPdpwv, in Watts/mm_pwv.
    #
    nu = nu_ghz*1e9
    #
    P_atm0 = np.trapz(inst_band*bnu_aomega(nu, tb1), nu) 
    P_atm1 = np.trapz(inst_band*bnu_aomega(nu, tb2), nu)   
    dPdpwv= (P_atm1-P_atm0)/dpwv
    #
    return dPdpwv
    
    
def calc_gpwv(dPdpwv, dPdTcmb):
    gpwv= dPdpwv/dPdTcmb
    return gpwv


    ''' (COPIED OVER FROM PREVIOUS DICT)
    nu_ghz= np.array(atm_array[0])
                    nu_ghz2= np.array(atm_array2[0])
                    nu= nu_ghz*1e9
                    nu2=nu_ghz2*1e9
                    outdict[site][pwv][t_type][band]['freqvec']= nu_ghz
                    model1= logistic_easy(nu_ghz, amp, nulow, nuhigh)
                    
                    outdict[site][pwv][t_type][band]['spectrum'] = model1
                    dB_dT= dBdT(2.7, nu)
                    dPdTcmb= np.trapz(model1*dB_dT*(c**2/nu**2), nu)
                       
                    P_atm0 = np.trapz(model1*bnu_aomega(nu, atm_array[2]), nu) #outdict[site][pwv]['tb_array']), nu)
                    P_atm1 = np.trapz(model1*bnu_aomega(nu, atm_array2[2]), nu2)
                        
                    dPdpwv= (P_atm1-P_atm0)/0.1
                '''
    
'''(functions that are now used inside of the dictionary itself)
def calc_highedgevalues(nu_low, nu_high, nuvec, tb, tb2, bandmodel, atmos_trans): #dnu, a, n, alpha):
    npts=50
    freqvec_high= np.linspace(nu_high-(nu_high*0.1), nu_high+(nu_high*0.1), npts)
    dPdTcmb_h= np.array([])
    dPdpwv_h= np.array([])
    #gpwv_h= np.array([])
   
    for freqhigh in freqvec_high: 
        nu_ghz= np.array(nuvec)
        nu_ghz2= np.array(nuvec)
        nu= nu_ghz*1e9
        nu2=nu_ghz2*1e9
        tb = np.array(tb)
        tb2 = np.array(tb2)
        #previous_model= logistic_easy(nu_ghz, amp, nulow, freqhigh)
        model1 = bandmodel*atmos_trans
        dB_dT= dBdT(2.7, nu)

        dPdTcmb_h= np.append(dPdTcmb_h, np.trapz(model1*dB_dT*(c**2/nu**2), nu))
        #print(dPdTcmb_h)
        P_atm0 = np.trapz(model1*bnu_aomega(nu, tb), nu) 
        P_atm1 = np.trapz(model1*bnu_aomega(nu, tb2), nu)                 
        dPdpwv_h= np.append(dPdpwv_h, (P_atm1-P_atm0)/0.1)
        print(dPdpwv_h)
        #gpwv_h= np.append(gpwv_h, dPdpwv_h/dPdTcmb_h)
        
    return dPdTcmb_h, dPdpwv_h #, gpwv_h

def calc_lowedgevalues(nu_low, nu_high, nuvec, tb, tb2, bandmodel, atmos_trans):
    npts=50
    freqvec_low= np.linspace(nu_low-(nu_low*0.1), nu_low+(nu_low*0.1), npts)
    dPdTcmb_l=np.array([])
    dPdpwv_l=np.array([])
    #gpwv_l=np.array([])
    
    for freqlow in freqvec_low: 
        nu_ghz= np.array(nuvec)
        nu_ghz2= np.array(nuvec2)
        nu= nu_ghz*1e9
        nu2=nu_ghz2*1e9
        tb = np.array(tb)
        tb2 = np.array(tb2)
        model1 = bandmodel*atmos_trans
        dB_dT= dBdT(2.7, nu)

        dPdTcmb_l= np.append(dPdTcmb_l, np.trapz(model1*dB_dT*(c**2/nu**2), nu))
        P_atm0 = np.trapz(model1*bnu_aomega(nu, tb), nu) 
        P_atm1 = np.trapz(model1*bnu_aomega(nu, tb2), nu)                 
        dPdpwv_l= np.append(dPdpwv_l, (P_atm1-P_atm0)/0.1)
        #gpwv_l= np.append(gpwv_l, dPdpwv_l/dPdTcmb_l)
        #print(gpwv_l)
        
    return dPdTcmb_l, dPdpwv_l #, gpwv_l

'''