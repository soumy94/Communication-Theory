from pylab import * 
from scipy.signal import butter, lfilter
def trapfilt(xt, Fs, fL, k, alfa):
    """ Delay compensated FIR LPF with trapezoidal frequency response.
    >>>>> yt, n = trapfilt(xt, Fs, fL, k, alfa) <<<<< 
    where yt: filter output y(t), sampling rate Fs 
    n: filter order 
    xt: filter input x(t), sampling rate Fs 
    Fs: sampling rate of x(t), y(t) 
    fL: cutoff frequency (-6 dB) in Hz 
    k: h(t) is truncated to |t| <= k/(2fL) 
    alfa: frequency rolloff parameter, linear rolloff 
    over range (1-alfa)fL <= |f| <= (1+alfa)fL 
    """ 
    ixk = round(Fs*k/float(2*fL))                             # Tail cutoff index 
    tt = arange(-ixk,ixk+1)/float(Fs)                         # Time axis for h(t) 
    n = len(tt)-1 # Filter order 
    ht = zeros(len(tt))
    ix = where(logical_and(tt>=-ixk,tt<ixk+1))[0]
    ht[int(len(ix)/2)] = 2*fL
    ixn = ix[0:n/2]
    ixp = ix[(n/2)+1:n+1]
    ix = hstack((ixn,ixp))
    ht[ix] = (sin(2*pi*fL*tt[ix])/(pi*tt[ix])) * (sin(2*pi*alfa*fL*tt[ix])/(2*pi*alfa*fL*tt[ix]))
    #ht[int(len(ix)/2)] = 2*fL
    if alfa == 0 :
        ixk = round(Fs*k/float(2*fL))
        ix = where(logical_and(tt>=-ixk,tt<ixk+1))[0]
        ixn = ix[0:160]
        ixp = ix[161:321]
        ix = hstack((ixn,ixp))
        TL = 1/float(2*fL)
        ht[int(len(ix)/2)] = 1     # At exception t=0, assign value of sinc directly at t =0 point
        ht[ix] = sin(pi*tt[ix]/TL)/(pi*tt[ix]/TL)
        
    yt = lfilter(ht, 1, hstack((xt, zeros(ixk)))) 
                                                              # Compute filter output y(t) 
    yt = yt[ixk:]                                             # Filter delay compensation 
    return yt, n                                              # Return y(t) and filter order





