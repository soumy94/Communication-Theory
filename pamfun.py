# File: pamfun.py
# Functions for pulse amplitude modulation (PAM)
from pylab import *
def pam10(an, FB, Fs, ptype, pparms=[]):
	"""
	Pulse amplitude modulation: a_n -> s(t), -TB/2<=t<(N-1/2)*TB,
	V1.0 for ’rect’, ’sinc’, and ’tri’ pulse types.
	>>>>> tt, st = pam10(an, FB, Fs, ptype, pparms) <<<<<
	where an:
	N-symbol DT input sequence a_n, 0 <= n < N
	FB:
	Baud rate of a_n, TB=1/FB
	Fs:
	sampling rate of s(t)
	ptype: pulse type (’rect’,’sinc’,’tri’)
	pparms not used for ’rect’,’tri’
	pparms = [k, beta] for ’sinc’
	k:
	"tail" truncation parameter for ’sinc’
	(truncates p(t) to -k*TB <= t < k*TB)
	beta: Kaiser window parameter for ’sinc’
	tt:
	time axis for s(t), starts at -TB/2
	st:
	CT output signal s(t), -TB/2<=t<(N-1/2)*TB,
	with sampling rate Fs
	"""
	N = len(an)			# Number of data symbols
	TB = 1/float(FB)		# Time per symbol
	ixL = ceil(-Fs*0.5*TB)		# Left index for time axis
	ixR = ceil(Fs*(N-0.5)*TB) 	# Right index for time axis
	tt = arange(ixL,ixR)/float(Fs) 	# Time axis for s(t)
	
	# ***** Conversion from DT a_n to CT a_s(t) *****
	ast = zeros(len(tt))		# Initialize a_s(t)
	ix = array(around(Fs*arange(0,N)*TB),int)	# Symbol center indexes
	ast[ix-int(ixL)] = Fs*an	# delta_n -> delta(t) conversion
	
	# ***** Set up PAM pulse p(t) *****
	ptype = ptype.lower()	# Convert ptype to lowercase
	
	# Set left/right limits for p(t)
	if (ptype=='rect'):
		kL = -0.5; kR = -kL
	else:
		kL = -1.0; kR = -kL
	
	# Default left/right limits
	ixpL = ceil(Fs*kL*TB)		# Left index for p(t) time axis
	ixpR = ceil(Fs*kR*TB)		# Right index for p(t) time axis
	ttp = arange(ixpL,ixpR)/float(Fs)	 # Time axis for p(t)
	pt = zeros(len(ttp))		# Initialize pulse p(t)
	if (ptype=='rect'):		# Rectangular p(t)
		ix = where(logical_and(ttp>=kL*TB, ttp<kR*TB))[0]
		pt[ix] = ones(len(ix))
	elif (ptype == 'tri'):
		pt = array([(1+i*1/TB) if i*1/TB+1 < 1.0 else (1-i*1/TB) for i in list(ttp)])
	elif (ptype=='sinc'):
		k=pparms[0]
		kL = -1.0*k; kR = -kL*k
		ixpL = ceil(Fs*kL*TB)		# Left index for p(t) time axis
		ixpR = ceil(Fs*kR*TB)		# Right index for p(t) time axis
		ttp = arange(ixpL,ixpR)/float(Fs)	 # Time axis for p(t)
		pt = zeros(len(ttp))
		beta=pparms[1]
		pt = array([sin(pi*t/TB)/(pi*t/TB) if t!= 0  else  1.0 for t in list(ttp)])
		pt=pt*kaiser(len(pt), beta)
	else:
		print("ptype ’%s’ is not recognized" % ptype)
	
	# ***** Filter with h(t) = p(t) *****
	st = convolve(ast,pt)/float(Fs) 	# s(t) = a_s(t)*p(t)
	st = st[-ixpL:ixR-ixL-ixpL] 		# Trim after convolution
	return tt, st
