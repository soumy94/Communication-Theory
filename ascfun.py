# File: ascfun.py
# Functions for conversion between ASCII and bits
from pylab import *
def asc2bin(txt, bits=8):
    """
ASCII text to serial binary conversion
>>>>> dn = asc2bin(txt, bits) <<<<<
where txt input text string
abs(bits) bits per char, default=8
bits > 0 LSB first parallel to serial
bits < 0 MSB first parallel to serial
dn binary output sequence
"""
    txtnum = array([ord(c) for c in txt]) # int array
    if bits > 0: # Neg powers of 2, increasing exp
        p2 = np.power(2.0,arange(0,-bits,-1))
    else: # Neg powers of 2, decreasing exp
        p2 = np.power(2.0,1+arange(bits,0))
    B = array(mod(array(floor(outer(txtnum,p2)),int),2),int8)
    # Rows of B are bits of chars
    dn = reshape(B,B.size)
    return dn # Serial binary output

def bin2asc(dn, bits=8, flg=1):
    """
Serial binary to ASCII text conversion
>>>>> txt = bin2asc(dn, bits, flg) <<<<<
where dn binary input sequence
abs(bits) bits per char, default=8
bits > 0 LSB first parallel to serial
bits < 0 MSB first parallel to serial
flg != 0 limit range to [0...127]
txt output text string
"""
    no_of_char = int(floor(size(dn)/abs(bits)))
    #if (mod(len(dn),8)!=0): # If length of input data bits not a multiple of 8, add zeros to make it a multiple of 8  
    #   no_of_zeros = len(dn) - 8*no_of_char
    #  if bits>0: #Append the zeros
    #     np.lib.pad(dn,(0,no_of_zeros),'constant',constant_values = (0))
    dec =[0]*(no_of_char) #Initializing array for decimal values with zeros 
    if bits > 0: # Pos powers of 2, increasing exp
        p2 = np.power(2.0,arange(0,bits,1))
    else: # Pos powers of 2, decreasing exp
        p2 = np.power(2.0,-1+arange(-bits,0,-1))
    for index in arange(no_of_char):# loop for converting 8 binary bits to decimal value sequenctially from the input bit sequence
        dec[index] = int(inner(dn[abs(bits)*index+arange(0,abs(bits))],p2))
    data_string = ''.join(chr(c) for c in dec) # Converting decimal value to ASCII for each element in the array
    return data_string # string returned
	
