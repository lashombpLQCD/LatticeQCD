# -*- coding: utf-8 -*-
"""
Lattice_QCD

This is a Metropolis Monte Carlo algorithm for Lattice QCD calculations
"""

import numpy as np
import matplotlib.pyplot as plt 
from math import exp


'''
Declaring initial parameters and constants
'''

a = 1/2                         # lattice spacing size
N = 20                          # Number of lattice sites
N_COR = 10                      # Number of correlations 
N_CF = 10000                      # Number of configurations/paths
EPS = 1.4                       # Metropolis sweep step size
BIN_S = 4                       # Bin size
LEN_G = int(N_CF/BIN_S)         # Length of binned G 
BTS_N = 100                     # Number of bootstrapping copies 


'''
Randomly generates a new configuration x and returns it if it meets the 
correct criteria. 
'''

def update(x): 
    for j in range(0, N): 
        old_x = x[j]
        old_Sj = S(j,x)
        x[j] += np.random.uniform(-EPS,EPS)
        dS = S(j,x) - old_Sj
        if dS > 0 and exp(-dS) < np.random.uniform(0,1): 
            x[j] = old_x
    return x


''' 
Defining the action S for harmonic oscillator 
'''             
def S(j,x): 
    jp = (j+1)%N                                
    jm = (j-1)%N
    return a*x[j]**2/2 + x[j]*(x[j]-x[jp]-x[jm])/a  


'''
Computes the correlation function for a time slice for one configuration
'''

#def compute_G(x,n):
#    g = 0; 
#    for j in range(0, N):
#        g = g + x[j]*x[(j+n)%N]
#    return g/N 

'''
Computes the correlation function using cubic sources and sinks. Comment out 
one compute_G in order to test the other correlation funciton so that only 
one is running. 
'''

def compute_G(x,n):
    g = 0; 
    for j in range(0, N):
        g = g + (x[j]**3)*(x[(j+n)%N]**3)
    return g/N 


'''
Paths are randomly generated using the update(x) and the correlation functions
are calculated for each path and time slice and stored in a 2D array G
'''

def MCpaths(x,G): 
    
    for j in range(0,N):                        #Empties array 
        x[j]=0; 
    for j in range(0, 5*N_COR):                 #Thermalizes array
        update(x) 
    for alpha in range(0, N_CF):                #Calculates G[alpha][n]
        for j in range(0, N_COR): 
            update(x)
        for n in range(0, N): 
            G[alpha, n] = compute_G(x,n) 





'''
Creates binned copies of G of a specified binsize. Notice that if binsize=N, 
then since the dimensions of G is N_cf*N and for the binned version it is 
(1/binsize)*N_cf*N, then if binsize = N_cf, then one simply obtains a matrix 
of size N which just corresponds to the average of each time slice over all 
alpha configurations. 
'''

def binning(G):
    G_binned = []
    
    for i in range(0, len(G), BIN_S):
        sumb_G = 0 
        
        for j in range(BIN_S):
            sumb_G += G[i+j]
        G_binned.append(sumb_G/BIN_S)
    
    return G_binned


'''
Creates bootstrapped copies of the binned G. 
'''

def bootstrap(G):
    G_bootstrap = []
    for i in range(LEN_G):
        G_bootstrap.append(G[np.random.randint(0, LEN_G)])
    return G_bootstrap        



'''
Calculates, for each time slice, G and averages over all configurations for 
each time slice and stores in an array cf_aveG
'''

def ave_t_G(G):
    t_aveG = []

    for n in range(0, N):
        aveG = 0
        for alpha in range(LEN_G):
            aveG += G[alpha][n]
        aveG /= LEN_G 
        t_aveG.append(aveG)
    return t_aveG            



'''
Calculates the difference excitation energy delta E(t) for each time slice t
and the uncertainty in G and delta E(t) and stores each as an array
'''


def stdCalc_dE(ave_tot_G, btsp):
    unc_dE = []
    unc_G = []
    dE = []  
    
    for n in range(N):
        err_1 = []
        err_2 = []
        dE.append((1/a)*np.log(np.abs(ave_tot_G[n]/ave_tot_G[(n+1)%N])))
        
        for beta in range(BTS_N):
            err_1.append((1/a)*(np.log(np.abs(btsp[beta][n] / btsp[beta][(n + 1)%N]))))
            err_2.append(btsp[beta][n])
            
        unc_dE.append(np.std(err_1))
        unc_G.append(np.std(err_2))
    
    return unc_dE, unc_G, dE 

    

'''
Declaration of Main function 
'''


def main(): 
    
    '''
    Defining path x as an array and the correlation function as a 2D array
    '''
    
    x = np.zeros(N);
    G_cor = np.zeros((N_CF, N)); 
    btsp = []
    ave_tot_G = []
    

    ##Creates the correlation functions for each time slice and configuration
    
    MCpaths(x, G_cor)
    
    
    binnedG = binning(G_cor)


    for _ in range(BTS_N):    
        btsp.append(ave_t_G(bootstrap(binnedG)))
        
    
    ave_tot_G = ave_t_G(binnedG)  

    unc_dE, unc_G, dE = stdCalc_dE(ave_tot_G, btsp)
    
    t = np.arange(0.0, 20.0, 0.01)
    exp_t = (1/2)*np.exp(-t)
    
    plt.figure()
    plt.plot(range(0, N), ave_tot_G, 'ro', markersize = 4)
    plt.plot(t, exp_t, lw=1, color = 'Black')
    plt.errorbar(range(0, N), ave_tot_G, yerr=unc_G, fmt ='o', markersize = 4, capsize = 3)
    plt.xlim(-.5, 10)
    plt.xlabel('t')
    plt.ylabel('G(t)')
    plt.show()
            
    #
    #
    #
    plt.plot(range(0, N), dE, 'bo', markersize = 4)
    plt.errorbar(range(0, N), dE, yerr = unc_dE, fmt = 'o', markersize = 4, capsize = 3)
    plt.ylabel('$\Delta$E(t)')
    plt.xlabel('t')
    plt.ylim(0.5, 1.5)
    plt.xlim(-.5, 5.5)
    plt.plot((-.5, 6), (1, 1), lw = 1, color = 'Black') 
    plt.show()        

   
if __name__ == "__main__": 
    main()
        
        
        
        
