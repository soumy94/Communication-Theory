from pylab import *

def showft(tt,xt,ff_lim):
    
    k1 = ff_lim[0]
    k2 = ff_lim[1]
    llim = ff_lim[2]

    N = len(tt)
    Fs = int((N-1)/float(tt[-1]-tt[0]))
    ixp = where(tt>=0)[0]
    ixn = where(tt<0)[0]
    xt = hstack((xt[ixp],xt[ixn]))
#compute freq component of the signal
    xf = fft(xt)/float(Fs)
#print('the fourier transform is', xf[1:10])
    ff = Fs*arange(-round(N/2),round(N/2)+1)/float(N)
#***** Compute |X(f)|, arg[X(f)] *****

    absXf = abs(xf)
    argXf = angle(xf)
#********swap the frequency and magnitude spectrum******

    a = absXf[0:len(absXf)/2]
    b = absXf[len(absXf)/2:len(absXf)]
    a,b = b,a
    absXf = append(a,b)

    a = argXf[0:len(argXf)/2]
    b = argXf[len(argXf)/2:len(argXf)]
    a,b = b,a
    argXf = append(a,b)
#if llim is less than 0, calculate log(magnitude spectrum)
    if llim<0:
        absXf = (20*log10(absXf/max(absXf)))#calculate log if llim<0
        for i in range(0,len(absXf)):
            if absXf[i]<=llim:
                absXf[i]=llim
                argXf[i]=0 
    else:
        for i in range(0,len(absXf)):#if llim >0
            if absXf[i]<=llim:
                absXf[i]=llim
                argXf[i]=0
#calculate the frequency bins for plotting in a given range

    con1 = round((k1*N)/Fs)
    con2 = round((k2*N)/Fs)
#***** plot magnitude/phase *****
    f1 = figure()
    af11 = f1.add_subplot(211)
    af11.plot(ff[(len(ff)/2)+(con1):(len(ff)/2)+(con2)],absXf[(len(ff)/2)+(con1):(len(ff)/2)+(con2)])         # Plot magnitude
    af11.grid()
    af11.set_ylabel('|X(f)|')
    strgt = 'FT Approximation, $F_s=$' + str(Fs) + ' Hz'
    strgt = strgt + ', N=' + str(N)
    strgt = strgt + ', $\Delta_f$={0:3.2f}'.format(Fs/float(N)) + ' Hz'
    af11.set_title(strgt)
    af12 = f1.add_subplot(212)
    af12.plot(ff[(len(ff)/2)+(con1):(len(ff)/2)+(con2)],180/pi*argXf[(len(ff)/2)+(con1):(len(ff)/2)+(con2)])  # Plot phase in degrees
    af12.grid()
    af12.set_ylabel('arg[X(f)] [deg]')
    af12.set_xlabel('f [Hz]')
    show()

def showeye(rt, Fs, FB, NTd=50, dispparms=[]): 
    """ Display eye diagram of digital PAM signal r(t) 
    >>>>> showeye(rt, Fs, FB, NTd, dispparms) <<<<<
    where rt: received PAM signal r(t)=sum_n a_n*q(t-nTB) 
    Fs: sampling rate for r(t) 
    FB: Baud rate of DT sequence a_n, TB = 1/FB
    NTd: Number of traces displayed
    dispparms = [delay, width, ylim1, ylim2]
    delay: trigger delay (in TB units, e.g., 0.5) 
    width: display width (in TB units, e.g., 3)
    ylim1: lower display limit, vertical axis
    ylim2: upper display limit, vertical axis
    """
    t0 = dispparms[0]/float(FB)                                  # Delay in sec 
    tw = dispparms[1]/float(FB)                                  # Display width in sec 
    dws = floor(Fs*tw)                                           # Display width in samples 
    tteye = arange(dws)/float(Fs)                                # Time axis for eye 
    trix = around(Fs*(t0+arange(NTd)/float(FB)))
    ix = where(logical_and(trix>=0, trix<=len(rt)-dws))[0] 
    trix = trix[ix]                                              # Trigger indexes within r(t)
    
    TM = rt[trix[0]:trix[0]+dws]                                # First trace
    for i in arange(NTd-1):
        TM = vstack((TM, rt[trix[i+1]:trix[i+1]+dws]))
                                                                 # Second trace 
    plot(FB*tteye, TM.T, '-b')
    ylim(dispparms[2],dispparms[3])
    grid() 
    show()
    
def showpsd0(xt, Fs, ff_lim, N):
    """ Plot (DFT/FFT approximation to) power spectral density (PSD) of x(t). 
    Displays S_x(f) either linear and absolute or normalized in dB. 
    >>>>> showpsd(xt, Fs, ff_lim, N) <<<<< 
    where xt: sampled CT signal x(t) 
    Fs: sampling rate of x(t) 
    ff_lim = [f1,f2,llim] 
    f1: lower frequency limit for display 
    f2: upper frequency limit for display 
    llim = 0: display S_x(f) linear and absolute 
    llim < 0: display 10*log_{10}(S_x(f))/max(S_x(f)) 
    in dB with lower display limit llim dB 
    N: blocklength 
    """
    
    # ***** Determine number of blocks, prepare x(t) ***** 
    N = int(min(N, len(xt)))                                               # N <= length(xt) needed 
    NN = int(floor(len(xt)/float(N))) 
                                                                           # Number of blocks of length N 
    xt = xt[0:N*NN]                                                        # Truncate x(t) to NN blocks 
    xNN = reshape(xt,(NN,N))                                               # NN row vectors of length N 
    # ***** Compute DFTs/FFTs, average over NN blocks ***** 
    Sxf = np.power(abs(fft(xNN)),2.0)                                      # NN FFTs, mag squared 
    if NN > 1:
        Sxf = sum(Sxf, axis=0)/float(NN) 
    Sxf = Sxf/float(N*Fs)                                                  # Correction factor DFT -> PSD 
    Sxf = reshape(Sxf,size(Sxf)) 
    ff = Fs*array(arange(N),int64)/float(N)                                # Frequency axis 
    if ff_lim[0] < 0:                                                      # Negative f1 case 
        ixp = where(ff<0.5*Fs)[0]                                          # Indexes of pos frequencies 
        ixn = where(ff>=0.5*Fs)[0]                                         # Indexes of neg frequencies 
        ff = hstack((ff[ixn]-Fs,ff[ixp]))                                  # New freq axis 
        Sxf = hstack((Sxf[ixn],Sxf[ixp]))                                  # Corresponding S_x(f) 
    # ***** Determine maximum, trim to ff_lim ***** 
    Px = cumsum(Sxf)            # calculating total Px
    #print(Px)
    maxSxf = max(Sxf)                                                      # Maximum of S_x(f) 
    ixf = where(logical_and(ff>=ff_lim[0], ff<ff_lim[1]))[0] 
    ff = ff[ixf]                                                           # Trim to ff_lim specs 
    Sxcf = Sxf[ixf]
    Pxf1f2 = cumsum(Sxcf)            # calculating Px between f1 and f2
    #print(Pxf1f2)
    
    if ff_lim[2]<0:
        Sxcf = (10*log10(Sxcf/max(Sxcf)))                   #calculate log if llim<0
        '''for i in range(len(Sxcf)):
            if Sxcf[i] < ff_lim[2]:
                Sxcf[i] = 0
            else:
                Sxcf[i] = Sxcf[i]'''
      
   
    
    # ***** Plot PSD ***** 
    strgt = 'PSD Approximation, $F_s=${:d} Hz'.format(Fs) 
    strgt = strgt + ', $\\Delta_f=${:.3g} Hz'.format(Fs/float(N)) 
    strgt = strgt + ', $NN=${:d}, $N=${:d}'.format(NN, N) 
    strgt = strgt + ', $P_x=${:f}'.format(Px[-1])
    strgt = strgt + ', $P_x(f1,f2)=${:f}'.format(Pxf1f2[-1])
    f1 = figure() 
    af1 = f1.add_subplot(111) 
    af1.plot(ff, Sxcf, '-b') 
    af1.grid() 
    af1.set_xlabel('f [Hz]') 
    af1.set_ylabel('$S_x$(f)') 
    if ff_lim[2]<0:
        af1.set_ylabel('$S_x$(f)[dB]')
        af1.set_ylim((ff_lim[2]))
        strgt = strgt + ', $P_x(f1,f2)=${:f} %'.format(Pxf1f2[-1]*100/float(Px[-1]))
    af1.set_title(strgt) 
    show()
    
    
    