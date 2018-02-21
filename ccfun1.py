
from pylab import *



def ccencod10(di, GD, trim=True): 
    """ Binary rate 1/n convolutional encoder with transfer
    function matrix G(D), V 1.0 
    >>>>> ci = ccencod10(di, GD, trim) <<<<< 
    where ci: multiplexed binary code sequence 
    di: binary (unipolar) data sequence 
    GD: array of n encoder polynomials 
    trim: if True, trim ci to n*len(di)
    
    Examples: GD = [[1,0],[1,1]] for G(D) = [1 1+D] 
    GD = [[1,0,1],[1,1,1]] for G(D) = [1+D^2 1+D+D^2] 
    GD = [[1,0,1,1],[1,1,0,1],[1,1,1,1]] for G(D) = [1+D^2+D^3 1+D+D^3 1+D+D^2+D^3]
    """
    n = len(GD)                                                  # Number of output bits for single bit in input
    k = len(di)                                                  # No. of bits in input sequence
    m = len(GD[0])                                               # memory constraint of encoder
    ci = []                                                      # initializing output sequence as empty array
    
    for i in range(k) :
        for l in range(n):
            xor = 0
            for z in range(0,m):
                if i-z < 0 :                                     # putting condition for initial state
                    state = 0
                else:
                    state = di[i-z]*GD[l][z] 
                xor = xor ^ state                               # Computing output
            ci.append(xor)                                      # Adding output bits to ci 
    
    return ci                                                      # returning output 
                


def ccdecod10(bi):
    """ Viterbi decoder for binary rate 1/2 convolutional code received
    as polar (binary 0 -> -1, binary 1 -> +1) PAM signal from AWGN
    channel. 
    Version 1.0 for encoder transfer function matrix G(D) = [1 1+D]
    >>>>> dihat, DM = ccdecod10(bi) <<<<<
    where dihat: ML estimate of binary data sequence 
    DM: array of final distance metrics 
    bi: received noisy polar binary (+A/-A) sequence 
    """ 
    n, K, m = 2, 2, 1                                       # Rate 1/n, constraint len K, memory m 
    N = int(floor(len(bi)/float(n)))                        # Number of codeword frames 
    CBM = [[-1,-1],[+1,+1],[-1,+1],[+1,-1]]                 # Code-bit (-1/+1) matrix 
    DM = array([0,1000])                                    # (Initial) distance metrics 
    dA = array(zeros((2,1)),int)                            # Competing data sequence array 
    ix2 = array([0,0,1,1],int)                              # State indexes (doubled) 
    for i in range(N):
        bDM = np.power(outer(ones(4),bi[n*i:n*(i+1)]) - CBM,2.0) 
        bDM = dot(bDM,ones((n,1)))                          # Branch distance metrics 
        bDM = reshape(bDM,size(bDM))                        # Convert to 1d array 
        tDM = DM[ix2] + bDM                                 # Tentative distance metrics 
        tDM = reshape(tDM,(2,2))                            # Reshape for path elimination 
        DM = amin(tDM,axis=0)                               # Select paths with smaller metric 
        ix = argmin(tDM,axis=0)                             # Indexes of smaller metric paths 
        dA = hstack((dA[ix,:],array([[0],[1]]))) 
                                                            # Competing data sequence update 
    dA = dA[:,1:]                                           # Discard first (dummy) column 
    ix = argmin(DM)                                         # Index of smallest metric 
    dihat = dA[ix,:]                                        # ML-decoded data sequence 
    dihat = reshape(array(dihat),size(dihat))               # Convert to 1d array 
    return dihat, DM


def ccdecod20(bi, GD):
    """
    Viterbi decoder for binary rate 1/n convolutional code with
    transfer function matrix G(D), received as polar (-1,+1)
    PAM signal from AWGN channel. Version 2.0
    >>>>> dihat, DM = ccdecod20(bi, GD) <<<<<
    where dihat: ML estimate of binary data sequence
    DM: array of final distance metrics
    bi: received noisy polar binary (+A/-A) sequence
    GD: array of n encoder polynomials
    Examples: GD = [[1,0],[1,1]] for G(D) = [1 1+D]
    GD = [[1,0,1],[1,1,1]] for G(D) = [1+D^2 1+D+D^2]
    GD = [[1,0,1,1],[1,1,0,1],[1,1,1,1]]
    for G(D) = [1+D^2+D^3 1+D+D^3 1+D+D^2+D^3]
    """
    k = len(GD[1])
    m = k-1
    n = len(GD)
    N = int(len(bi)/n)
    no_of_states = pow(2,m)
    prev_bit = zeros(no_of_states,int)
    CBM = zeros((2*no_of_states,n),int)
    prev_states = zeros((no_of_states,n),int)
    for index in arange(no_of_states):
        bits = [int(x) for x in binary_repr(index,m)]
        state_bits = bits[::-1]
        
        prev_bit[index] = state_bits[0]
        for transition in arange(2):
            prev_states[index,transition] = sum(pow(2,arange(m))*array(r_[state_bits[1:m],transition]))        
            CBM[n*index+transition,0] = 2*mod(sum(r_[state_bits,transition]*GD[0]),2)-1
            CBM[n*index+transition,1] = 2*mod(sum(r_[state_bits,transition]*GD[1]),2)-1
    DM = zeros(no_of_states) # (Initial) distance metrics
    tDM = zeros(no_of_states) # (Initial) distance metrics
    DM[0] = 1000
    dA = zeros((no_of_states,N),int) # Competing data sequence array
    tdA = zeros((no_of_states,N)) # Competing data sequence array
    for i in range(N):
        bDM = CBM*bi[n*i:n*(i+1)]
        bDM = dot(bDM,ones((n,1))) # Branch distance metrics
        for state in arange(no_of_states):
            temp1 = DM[prev_states[state,0]]+bDM[2*state]
            temp2 = DM[prev_states[state,1]]+bDM[2*state+1]
            tDM[state] = max(temp1,temp2)
            ix = argmax((temp1,temp2))
            #print('tda =',tdA)
            dA[state,:] = r_[prev_bit[state],tdA[prev_states[state,ix],0:N-1]] 
        tdA = 1*dA
        DM = 1*tDM
    ix = argmax(tDM) 
    dihat = dA[ix,:] 
    dihat= dihat[::-1]
    return dihat


