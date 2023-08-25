
"""Defines the PESCy functions."""
__all__ = [
    "calcS",
    "calcH",
    "calcCofH",
    "Cmaxmin",
    "generateCurves",
    "constructPatternCount",
    "calcS_fromPatternCount",
    "calcPESCcurves" 
]

import numpy as np
from math import factorial
from collections import Counter
import matplotlib.pylab as plt

def calcS(data,n=5,delay=1):		
    '''
    function Cjs - Returns the Shannon Permutation Energy
    Input:
        data  - 1-D array
        n     - permutation number (default=5)
        delay - integeter delay (default=1)
    Output:
        Sp - Shannon permuation entropy
        Se - Shannon + Uniform permutation entropy
    ''' 
    N=factorial(n)
    T=np.array(data)
    if len(T.shape)>1:
        raise TypeError( 'Data must be a 1-D array')
    t = len(T)
    Ptot = t - delay*(n - 1)    #Total number of order n permutations in T
    print('Number of permutations = ', Ptot)
    invPtot=1./Ptot     #Inverse for later calcuations
    A = []			 #Array to store each permutation
    S = 0.
    Se = 0.
    
    for i in range(Ptot):	#Will run through all possible n segments of T
        A.append(''.join(T[i:i+(n-1)*delay+1:delay].argsort().astype(str)))
    #Count occurance of patterns
    count=Counter(A)
    #print len(count)
    #Calculate S from the count
    for q in iter(count.values()):
        q*=invPtot #convert to probability
        S += -q * np.log2(q)
        q+=1./N
        q/=2
        Se += -q * np.log2(q)
    for i in range(len(count),N):
        q=1./2/N
        Se += -q * np.log2(q)
    return S,Se

def calcH(data,n=5,delay=1):		
    '''
    function Cjs - Returns the Normalized Shannon Permutation Energy, H
    Input:
        data  - 1-D array
        n     - embedding dimension (default=5)
        delay - embedding delay (default=1)
    Output:
        Sp - Shannon permuation entropy
        Se - Shannon + Uniform permutation entropy
    ''' 
    N=factorial(n)
    T=np.array(data)
    if len(T.shape)>1:
        raise TypeError( 'Data must be a 1-D array')
    t = len(T)
    Ptot = t - delay*(n - 1)    #Total number of order n permutations in T
    print('Number of permutations = ', Ptot)
    invPtot=1./Ptot     #Inverse for later calcuations
    A = []			 #Array to store each permutation
    S = 0.
    Se = 0.
    
    for i in range(Ptot):	#Will run through all possible n segments of T
        A.append(''.join(T[i:i+(n-1)*delay+1:delay].argsort().astype(str)))
    #Count occurance of patterns
    count=Counter(A)
    #print len(count)
    #Calculate S from the count
    for q in iter(count.values()):
        q*=invPtot #convert to probability
        S += -q * np.log2(q)
        q+=1./N
        q/=2
        Se += -q * np.log2(q)
    for i in range(len(count),N):
        q=1./2/N
        Se += -q * np.log2(q)
    return S/np.log2(N),Se/np.log2(N)

def calcCofH(data,n,delay=1):
    '''
    function Cjs - Returns the normalized Jensen-Shannon statistical complexity
    Input:
        data  - array
        n     - embedding dimension (default=5)
        delay - embedding delay (default=1)
    Output:
        C - Normalized Jensen-Shannon complexity
        H - Normalized Shannon Perumation Entropy
    '''		
    N  = factorial(n)
    S, Se  = calcS(data,n,delay)   
    C = -2.*((Se - 0.5*S - 0.5*np.log2(N))
            /((1 + 1./N)*np.log2(N+1) - 2*np.log2(2*N) 
            + np.log2(N))*(S/np.log2(N)))

    return S/np.log2(N), C

def Cmaxmin(nsteps,n):	
	
    N = factorial(n)
    Cmaxx = np.zeros((N-1)*nsteps)
    Cmaxy = np.zeros((N-1)*nsteps)
    Cminx = np.zeros(nsteps)
    Cminy = np.zeros(nsteps)
	
    for i in np.arange(nsteps):
        pk = 1./N + i*(1.-(1./N))/nsteps
        pj = (1. - pk)/(N - 1.)
        S = -pk * np.log2(pk) - (N - 1.) * pj * np.log2(pj)
        qk = pk/2. + 1./(2.*N)
        qj = pj/2. + 1./(2.*N)
        Scom = -qk * np.log2(qk) - (N - 1.) * qj * np.log2(qj)
        Cminx[i] = S / np.log2(N)
        Cminy[i] = -2. * (S/np.log2(N)) * (Scom - 0.5*S - 0.5*np.log2(N)) \
        /((1 + 1./N)*np.log2(N+1) - 2*np.log2(2*N) + np.log2(N))	
		
    for i in np.arange(1,N):
        for l in np.arange(nsteps):
            pk = l*(1./(N-i+1.))/nsteps
            pj = (1. - pk)/(N - i)
            if pk ==0.:
                S = -(N - i) * pj * np.log2(pj)
            else:
                S = -pk * np.log2(pk) - (N - i) * pj * np.log2(pj)
            qk = pk/2. + 1./(2.*N)
            qj = pj/2. + 1./(2.*N)
            Scom = -qk * np.log2(qk) - (N - i) * qj * np.log2(qj) - \
            (i-1)*(1./(2.*N))*np.log2(1./(2.*N))
            #print (i-1.)*nsteps+l
            Cmaxx[(i-1)*nsteps+l] = S / np.log2(N)
            Cmaxy[(i-1)*nsteps+l] = -2.*(S/np.log2(N))*(Scom - 0.5*S - 0.5*np.log2(N)) \
            /((1. + 1./N)*np.log2(N+1.) - 2.*np.log2(2.*N) + np.log2(N))
			
    return Cminx, Cminy, Cmaxx, Cmaxy

def generateCurves(n=5):				# Creates a blank CH plane with maximum and minimum curves for the given embedding dimension, with n=5 as the default
	
	Cminx, Cminy, Cmaxx, Cmaxy = Cmaxmin(1000,n)

	plt.figure(1)
	plt.plot(Cminx,Cminy,'k-',Cmaxx,Cmaxy,'k-')
	plt.xlabel(r"Normalized Permutation Entropy, $H$", fontsize=12)
	plt.ylabel(r"Statistical Complexity, $C$", fontsize=12)
	#plt.axis([0,1.0,0,0.45])
	#plt.xticks(np.arange(0,1.1,0.1))
	#plt.yticks(np.arange(0,0.45,0.05))
	#savefile='CHn5_blank'
	#plt.savefig(str(savefile)+'.png')
    
    
def constructPatternCount(data,n=5,delay=1):#originally called PE_dist(), but added n output for new version
    '''
    Parameters
    ----------
    data : 1-D array
        time-series data.
    n : integer,optional
        embedding dimension. The default is 5.
    delay : integer, optional
        embedday delay. The default is 1.

    Returns
    -------
    count - A Count occurance of patterns
    Ptot - total number of permutations
    n - embedding delay used

    '''
    
    T=np.array(data)
    if len(T.shape)>1:
        raise TypeError('Data must be a 1-D array')
    t = len(T)
    Ptot = t - delay*(n - 1)    #Total number of order n permutations in T
    #print 'Number of permutations = ', Ptot
    A = []			 #Array to store each permutation
    
    for i in range(Ptot):	#Will run through all possible n segments of T
        A.append(''.join(T[i:i+(n-1)*delay+1:delay].argsort().astype(str)))
    #Count occurance of patterns
    count=Counter(A)
    return count,Ptot,n

def calcS_fromPatternCount(count,tot_perms,n):#originally called PE_calc_only(), though n input was added for this new version
    '''
    

    Parameters
    ----------
    count : A Counter object
        Count occurance result from constructPatternCount()
    tot_perms : integer
        total number of permutations from constructPatternCount()
    n : integer
        embedding dimension from constructPatternCount()

    Returns
    -------
    S - shannon permutation entropy
    Se - Shannon + Uniform permuation Entropy

    '''
    Ptot=tot_perms
    N=factorial(n)
    invPtot=1./Ptot     #Inverse for later calcuations
    S = 0.
    Se = 0.
    for q in iter(count.values()):
        q*=invPtot #convert to probability
        S += -q * np.log2(q)
        q+=1./N
        q/=2
        Se += -q * np.log2(q)
    for i in range(len(count),N):
        q=1./2/N
        Se += -q * np.log2(q)
    return S,Se

def calcPESCcurves(data, n=5, min_delay=1, max_delay=100):
    """
    function calcPESCcurves - Returns PE(tau) and SC(tau)

    Parameters
    ----------
    data : array
        Timeseries data for PESCy
    n : int, optional
        Embedding dimension, by default 5
    min_delay : int, optional
        Smallest embedding delay to loop through, by default 1
    max_delay : int, optional
        Largest embedding delay to loop through, by default 100

    Returns
    -------
    C(tau) : array
        Normalized Jensen-Shannon complexity as function of embedding delay tau
    H : array
        Normalized Shannon Perumation Entropy as function of embedding delay tau
    """
    nfac = factorial(n)
    delay_array = np.arange(min_delay,max_delay)
    num_delays=len(delay_array)
    PEs=np.zeros([num_delays])
    SCs=np.zeros([num_delays])
    for loop_delay in np.arange(len(delay_array)):		
        if (loop_delay%100)==0: print( 'On Delay ',delay_array[loop_delay])
        permstore_counter = []
        permstore_counter = Counter(permstore_counter)
        tot_perms = 0
        arr,nperms,n0 = constructPatternCount(data,n=5,delay=delay_array[loop_delay])
        permstore_counter = permstore_counter+arr
        tot_perms = tot_perms+nperms
        PE_tot,PE_tot_Se = calcS_fromPatternCount(permstore_counter,tot_perms,n0)
        C =  -2.*((PE_tot_Se - 0.5*PE_tot - 0.5*np.log2(nfac))
                    /((1 + 1./nfac)*np.log2(nfac+1) - 2*np.log2(2*nfac) 
                    + np.log2(nfac))*(PE_tot/np.log2(nfac)))
        PEs[loop_delay]=PE_tot/np.log2(nfac)
        SCs[loop_delay]=C
    return PEs,SCs