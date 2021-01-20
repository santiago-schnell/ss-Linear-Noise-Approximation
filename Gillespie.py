#! /usr/bin/python3
######################################################################################################################################################################

import numpy as np
from numpy import linalg as LA
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random
from Sorter import sort
from Sorter import sortDRIVER


######################################################################################################################################################################


def Evolve(Sol,r1,r2,k1,k2,k3):   #Sol = [t,S,E,C,P]

        a1 = Sol[1]*Sol[2]*k1
        a2 = Sol[3]*k2
        a3 = Sol[3]*k3
        a0 = a1 + a2 + a3
        p1 = a1/a0
        p2 = (a1+a2)/a0
        p3 = (a2+a3)/a0

        

        tau = (math.log(1.0/r1))/a0


        if 0 <= r2 < p1:
           nSol = [Sol[0]+tau,Sol[1]-1,Sol[2]-1,Sol[3]+1,Sol[4]]
           flag = 1

        if p1 <= r2 < p2:
            nSol = [Sol[0]+tau,Sol[1]+1,Sol[2]+1,Sol[3]-1,Sol[4]]
            flag = 1

        if p2 <= r2 < 1:
            nSol = [Sol[0]+tau,Sol[1],Sol[2]+1,Sol[3]-1,Sol[4]+1]
            flag = 0

        return [nSol,flag]   #return updated row vector


s0 = 100
e0 = 10
Sol = [0,s0,e0,0,0]
k1  = 0.01
k2  = 2.0    # this is k-1
k3  = 1.0    # this is the calalytic rate constant, k2





def Gillespie(Sol,k1,k2,k3,s0,stop):

   

    tSeries = Sol

    pTime=np.ones(stop+1)
    Prod=np.ones(stop+1)

   
    pTime[0]   = 0.0
    Prod[0]    = 0.0
    i = 0

    while Sol[4]<stop:
       r1 = random.uniform(0, 1)
       r2 = random.uniform(0, 1)
       [newSol,flag] = Evolve(Sol,r1,r2,k1,k2,k3)
       tSeries = np.vstack((tSeries, newSol))  # Store results                                  
       Sol = newSol
    
       
    
       if flag ==0:
            i = i+1
            pTime[i] = newSol[0]
            Prod[i]  = i

            
   
    t = tSeries[:,0]
    S = tSeries[:,1]
    E = tSeries[:,2]
    C = tSeries[:,3]
    P = tSeries[:,4]
    sTime = sortDRIVER(S,t,stop)

    


    return [pTime,Prod,sTime]
    
#[t,Prod] = Gillespie(Sol,k1,k2,k3,s0,s0)



#######################################################################################################################################################################



#fig=plt.figure(1)
##ax = plt.subplot(111, xlabel='$t$', ylabel='$N_M$')
#for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
 #            ax.get_xticklabels() + ax.get_yticklabels()):
 #            item.set_fontsize(15),
#plt.step(t, S, 'b', data=None, where='pre',linewidth=0.25)
#plt.step(t, Prod, '-k', data=None, where='post',linewidth=1.0)
#plt.step(t, C, 'g', data=None, where='pre',linewidth=0.25)
#plt.tight_layout()
#plt.show()
