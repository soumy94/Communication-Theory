
# coding: utf-8

# In[3]:

# File: filtspecs.py
# Module for filter specfifications 
# FIR: Determine filter taps 
# IIR: Determine numerator (b) and denominator (a) 
# polynomial coefficients 
from pylab import * 
from scipy.signal import butter,lfilter
def trapfilt_taps(N, phiL, alfa):
    """ Returns taps for order N FIR LPF with trapezoidal frequency
    response, normalized cutoff frequency phiL = fL/Fs, and rolloff
    parameter alfa. 
    >>>>> hLn = trapfilt_taps(N, phiL, alfa) <<<<< 
    where N: filter order 
    phiL: normalized cutoff frequency (-6 dB) 
    alfa: frequency rolloff parameter, linear rolloff
    over range (1-alfa)phiL <= |f| <= (1+alfa)phiL 
    """
   
    tt = arange(-(N/2),(N/2)+1)                         # Time axis for h(t) 
    ht = zeros(len(tt))
    ix = where(logical_and(tt>=-N/2,tt<(N/2)+1))[0]
    ht[int(len(ix)/2)] = 1
    ixn = ix[0:N/2]
    ixp = ix[(N/2)+1:N+1]
    ix = hstack((ixn,ixp))
    ht[ix] = (sin(2*pi*phiL*tt[ix])/(pi*2*phiL*tt[ix])) * (sin(2*pi*alfa*phiL*tt[ix])/(2*pi*alfa*phiL*tt[ix]))
    
    if alfa == 0 :
        ix = where(logical_and(tt>=-N/2,tt<(N/2)+1))[0]
        ixn = ix[0:N/2]
        ixp = ix[(N/2)+1:N+1]
        ix = hstack((ixn,ixp))
        ht[int(len(ix)/2)] = 1     # At exception t=0, assign value of sinc directly at t =0 point
        ht[ix] = sin(pi*tt[ix]*2*phiL)/(pi*tt[ix]*2*phiL)
   
    
    return ht                                              # Return h(t)                                             






