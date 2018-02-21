#File: pcmfun.py
# Functions for conversion between m(t) and PCM representation
from pylab import *
def mt2pcm(mt, bits=8):
    """
    Message signal m(t) to binary PCM conversion
    >>>>> dn = mt2pcm(mt, bits) <<<<<
    where mt normalized (A=1) "analog" message signal
    bits number of bits used per sample
    dn binary output sequence in sign-magnitude
    form, MSB (sign) first
    """
    step_size = 2/2**(bits)  # total dynamic range/number of levels
    qt = (mt*2**(bits-1))   #Scaling the signal by dividing step size to make step size 1
    corr = array([int(sign(c)) for c in mt]) # Since we cannot represent 2**(bits-1) value, we apply different correction factor
                                             # for positive values and negative values 
    dt = qt -corr/2 
    dt = array([int(floor(c)) for c in dt]) 
    #plot(arange(0,int(Fs/f0)),4*st[0:int(Fs/f0)],arange(0,int(Fs/f0)),qt[0:int(Fs/f0)],'--k',arange(0,int(Fs/f0)),dt[0:int(Fs/f0)],'*-b')
    #ylim([-4,4])
    #grid ()
    bit_seq = []
    # Decimal to binay conversion - first MSB Rule
    for index in arange(size(qt)):
        if sign(dt[index]) <0:
            sign_bit = 1
            dt[index] = abs(dt[index]+1);
        else:
            sign_bit = 0
        p2 = np.power(2.0,arange(-bits+1,0,1)+1)
        B = array(mod(floor(outer(dt[index],p2)),2),int)
        B= reshape(B,size(B))
        C = insert(B,0,sign_bit)
        bit_seq = insert(bit_seq,len(bit_seq),C)
    return bit_seq

def pcm2mt(dn, bits=8):
    """
    Binary PCM to message signal m(t) conversion
    >>>>> mt = pcm2mt(dn, bits) <<<<<
    where dn binary output sequence in sign-magnitude
    form, MSB (sign) first
    bits number of bits used per sample
    mt normalized (A=1) "analog" message signal
    """
    # n-bit Binary to  Decimal conversion - first MSB Rule
    no_of_dec = int(floor(size(dn)/abs(bits)))
    dec =[0]*(no_of_dec) #Initializing array for decimal values with zeros 
    # Pos powers of 2, decreasing exp
    p2 = np.power(2.0,-1+arange(bits-1,0,-1))
    for index in arange(no_of_dec):# loop for converting nbit binary value to decimal value sequenctially from the input bit sequence
        dec_num = int(inner(dn[bits*index+arange(1,bits)],p2))
        sign_bit = dn[bits*index]
        if sign_bit>0:
            dec[index] = -dec_num -1
        else:
            dec[index] = dec_num
    dec = reshape(dec,size(dec))
    corr_factor = array([int(sign(c)) for c in dec]) # Since we cannot represent 2**(bits-1) value, we apply different correction factor
                                             # for positive values and negative values 
    mthat = dec
    mthat = array([int(ceil(c)) for c in mthat]) 
    mthat = mthat/2**(bits-1)
    return mthat
