#! /usr/bin/python3
######################################################################################################################################################################

import numpy as np
from numpy import linalg as LA
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random
import warnings



###############################################################################

def sort(S,t,num):
     
     tindex = np.where(S == num)[0]
     times = np.ones(len(tindex))

     for i in range(len(tindex)):
         times[i] = t[tindex[i]]
         
       
     return times


def sortDRIVER(S,t,s0):

    j=s0
   
    
    Atimes = np.ones(s0+1)

    for i in range(len(Atimes)):
        Stimes=sort(S,t,j)
        Atimes[i] = np.mean(Stimes)
        j = j-1
       
      

    return Atimes
    



#S = np.array([2, 1, 2, 3, 2, 5, 6, 2, 8, 2])

#print(len(S))
#t = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 0])
#num =2

#y = sort(S,t,num)
#print(y)
