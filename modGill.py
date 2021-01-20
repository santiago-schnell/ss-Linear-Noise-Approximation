#! /usr/bin/python3
######################################################################################################################################################################

import numpy as np
from numpy import linalg as LA
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random


######################################################################################################################################################################


def mEvolve(Sol,r1,k1,k2,k3,e0,s0):   #Sol = [t,S,P]

    
        Km = (k2+k3)/k1
        KS =  k2/k1
        #a0 = k3*Sol[1]*e0/(KS + Sol[1])
        #a0 = k1*e0*Sol[1]
        #a0  = k3*Sol[1]*e0*(Sol[1]+KS)/(e0*KS + (Sol[1]+KS)*(Sol[1]+KS))
       

        alpha = e0+s0+KS-Sol[2]
        beta  = alpha**2 - 4.0*e0*(s0-Sol[2])
        rate     = alpha-math.sqrt(beta)

        a0 = 0.5*k3*rate

        tau = (math.log(1.0/r1))/a0
        nSol = [Sol[0]+tau,Sol[1]-1,Sol[2]+1]

        return nSol   #return updated row vector




s0 = 100
e0 = 1
Sol = [0,s0,0]
k1  = 0.01
k2  = 1.0 # this is km1
k3  = 1.0  # this is the calalytic rate constant



def sQSSGillespie(Sol,k1,k2,k3,e0,s0):

    tSeries = Sol

    while Sol[2]<s0:
        r1 = random.uniform(0, 1)
        newSol = mEvolve(Sol,r1,k1,k2,k3,e0,s0)  
        tSeries = np.vstack((tSeries, newSol))  # Store results                                   
        Sol = newSol
    

        t = tSeries[:,0]
        S = tSeries[:,1]
        P = tSeries[:,2]

    return[t,P,S]
     

#######################################################################################################################################################################


#[t,P] = sQSSGillespie(Sol,k1,k2,k3,e0,s0)

#fig=plt.figure(1)
#ax = plt.subplot(111, xlabel='$t$', ylabel='$N_M$')
#for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
  #            ax.get_xticklabels() + ax.get_yticklabels()):
#              item.set_fontsize(15),
#plt.step(t, S, 'b', data=None, where='post')
#plt.step(t, P, 'r', data=None, where='post')
#plt.tight_layout()
#plt.show()

