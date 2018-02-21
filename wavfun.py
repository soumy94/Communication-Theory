
#file: wavfun.py
#functions for writing wav files in python

from pylab import * 
import struct 
import wave
def wavread(fname): 
    fh = wave.open(fname,'rb')
    (nchannels, sampwidth, framerate, nframes, comptype, compname) = fh.getparams()
    if sampwidth == 2: 
        frames = fh.readframes(nframes * nchannels) 
        dn = struct.unpack_from('%dh' % nframes*nchannels, frames)
        if nchannels > 1:
            out = array([dn[i::nchannels] for i in range(nchannels)])/float(2**15) 
        else:
            out = array(dn)/float(2**15) 
    else:
        print('not a 16 bit wav-file') 
        out = [0] 
    fh.close() 
    return out, framerate

def wavwrite(data, framerate, fname): 
    fh = wave.open(fname, 'wb') 
    if len(data.shape) == 1: 
        m = data.size 
        n = 1; 
    else: 
        m, n = data.shape 
        if m < n:                             # make column vectors 
            data = data.transpose() 
            m, n = data.shape 
    dn = reshape(data,data.size)        # interleave channels 
    dn = around(dn*2**15).astype(dtype='int16') 
    fh.setparams((n, 2, framerate, data.size, 'NONE', 'not compressed')) 
    frames = struct.pack('h'*data.size, *dn) 
    fh.writeframesraw(frames) 
    fh.close()



