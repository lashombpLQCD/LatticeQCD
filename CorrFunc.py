"""
Python Calculation of Correlation Functions

Author: Paul Lashomb

"""

"""

Code calculates the correlation function for a simple harmonic oscillator from Quantum Mechanics
and computes the difference in energy levels from the ground state and first excited state
from the correlation functions

"""


from scipy import integrate 
from scipy.integrate import nquad
from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt 
from math import exp, log


'''
result_array = np.array([])

for x in range(1, 4):     
    #print(quad(lambda x: x**2 + i, 0, 2))
    result, err = quad(lambda x1: np.exp(-((1/2)*(x1-x)**2 +((x**2)/2 + (x1**2)/2))) , -np.inf, np.inf)
    result_array = np.append(result_array, result)

print(result_array)


plt.plot([1, 2, 3], result_array, 'ro')
plt.axis([0, 4, 0, 1])
plt.ylabel('some numbers')
plt.show()


#variable_array = np.array([x0, x1, x2, x3, x4, x5, x6])

#print(variable_array)


x = np.array([2, 2, 2, 2, 2, 2, 3, 5, 3, 5, 2, 6, 1, 2, 3, 4, 5, 2, 6, 3])

'''


a = 1/2;
N = 20
N_cor = 10;
N_cf = 1000; 
eps = 1.4;


z = np.array([])
x = np.zeros(N);
G = np.zeros((N_cf, N)); 


def update(x): 
    for j in range(0, N): 
        old_x = x[j]
        old_Sj = S(j,x)
        x[j]=x[j]+np.random.uniform(-eps,eps)
        #print(x[j])
        dS = S(j,x) - old_Sj
        if dS > 0 and exp(-dS) < np.random.uniform(0,1): 
            x[j] = old_x
    #print(x)
    return x


''' 
Defining action S for harmonic oscillator 
'''             
def S(j,x): 
    jp = (j+1)%N
    jm = (j-1)%N 
    return a*x[j]**2/2 + x[j]*(x[j]-x[jp]-x[jm])/a 


''' 
Defining improved action Simp for harmonic oscillator 
'''             
def S_imp(j,x): 
    jp = (j+1)%N
    jm = (j-1)%N 
    return a*x[j]**2/2 + x[j]*(x[j]-x[jp]-x[jm])/a 


def compute_G(x,n):
    g = 0; 
    for j in range(0, N):
        g = g + x[j]*x[(j+n)%N]
    return g/N 

def MCpaths(x,G): 
    
    for j in range(0,N):
        x[j]=0; 
    for j in range(0, 5*N_cor):
        update(x) 
    for alpha in range(0, N_cf): 
        for j in range(0, N_cor): 
            update(x)
        for n in range(0, N): 
            G[alpha, n] = compute_G(x,n) 
#    for n in range(0, N): 
#        avg_G = 0
#        for alpha in range(0, N_cf): 
#            avg_G = avg_G + G[alpha][n] 
#        avg_G = avg_G/N_cf
#        z = np.append(z, avg_G)
#      # print(z)
#      # print("G(%d) = %g" % (n, avg_G))


'''
Calculates, for each time slice, G and averages over all configurations for 
each time slice and stores in an array cf_aveG
'''

def ave_t_G(G):
    t_aveG = np.zeros(N)
    for n in range(0, N):
        t_aveG[n] = (G.sum(axis = 0))[n]
        t_aveG[n] /= G.shape[0]
    return t_aveG            


def bootstrap(G):
    N_cf = len(G) 
    G_bootstrap = []
    for i in range(0, N_cf):
        alpha = int(np.random.uniform(0, N_cf))
        G_bootstrap.append(G[alpha])
    return G_bootstrap        
        
def bin(G, binsize):
    G_binned = []
    for i in range(0, len(G), binsize):
        G_avg = 0 
        for j in range(0, binsize):
            G_avg = G_avg + G[i+j]
        G_binned.append(G_avg/binsize)
    return G_binned

        
        
MCpaths(x, G)


ave_t_G(G)

#print(G[:-1])

#print(G.sum(axis=1))
       

def avg(G): 
    return np.sum(G)/len(G)

def sdev(G): 
    g = np.asarray(G)
    return np.absolute(avg(g**2)-avg(g)**2)**0.5


print("avg G\n", avg(G))
print('avg G (binned)\n', avg(bin(G,4)))
print("avg G (bootstrap) \n", avg(bootstrap(G)))





#def deltaE(G): 
#    avgG = avg(G)
#    adE = np.log(np.absolute(avgG[:-1]/avgG[1:]))
#    return adE/a
#
#
#
#print("Delta D/n", deltaE(G))
#
#print(excitation_E(z))
#
#     
#k= MCpaths(x, G)
#

'''
Calculates the difference excitation energy delta E(t) for each time slice t
and stores as an array
'''
def dE(): 
    E = np.array([])
    for i in range(0, N):
        ip = (i+1)%N
        E = np.append(E, (1/a)*np.log(np.abs(((ave_t_G(G))[i])/((ave_t_G(G))[ip]))))
    return E
#    
#h = np.array([])
#
#p = (1/a)*test(k)
#
#for j in range(0, 6):
#    h = np.append(h, p[j])
#    print(h)
#
#
#
#    
plt.plot(range(0, N), ave_t_G(G), 'ro')
plt.xlabel('t')
plt.ylabel('G(t)')
plt.show()
#        
#
#
#
plt.plot(range(0, N), dE(), 'ro')
plt.ylabel('some numbers')
plt.ylim(0.5, 1.5)
plt.xlim(-.5, 5.5)
plt.plot((-.5, 6), (1, 1), color = 'b') 
plt.show()        
#       
#        
#        
        
        
        
        
        
