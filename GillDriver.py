#! /usr/bin/python3

import numpy as np
from numpy import linalg as LA
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random
from Gillespie import Gillespie
from Gillespie import Evolve
from modGill import sQSSGillespie
from modGill import mEvolve
from Sorter import sort
from Sorter import sortDRIVER


######################################################################################
#################################### Initial Data ####################################
######################################################################################

Max = 10000  #number of runs

s0 = 40
e0 = 40
k1  = 0.05
k2  = 10.0   # this is k-1
k3  = 0.05   # this is the calalytic rate constant, k2

Km = (k2+k3)/k1
eps = e0/(Km+s0)

lam = int(0.5*(e0+s0+Km)-0.5*math.sqrt((e0+s0+Km)**2 - 4*e0*s0))

u=0.0

sInitial = s0-u*lam
eInitial = e0-u*lam
cInitial = u*int(lam)
pInitial = 0


Sol = [0,sInitial,eInitial,cInitial,pInitial]



######################################################################################
############### First Loop; computes mean and STD directly from Gillespie ############
######################################################################################

[times,pSeries,Stimes] = Gillespie(Sol,k1,k2,k3,sInitial,s0)


j = 0

while j < Max:

  Sol = [0,sInitial,eInitial,cInitial,pInitial]

  [t,P,SubsT] = Gillespie(Sol,k1,k2,k3,sInitial,s0)

  times = np.vstack((times, t))
  Stimes = np.vstack((Stimes, SubsT))


  j = j + 1



Tave = times.mean(0)
Tdev = times.std(0)

Stave = Stimes.mean(0)
Stead = Stimes.std(0)


######################################################################################
################ Second Loop; compute mean and STD from reduced CME ##################
######################################################################################


si = s0-0*lam
print(si)

qSol = [0,si,0]

[QSSt,QSSp,QSSs] = sQSSGillespie(qSol,k1,k2,k3,e0,si)



k = 0


while k < Max:

  

  [Qt,QP,QS] = sQSSGillespie(qSol,k1,k2,k3,e0,si)

  QSSt = np.vstack((QSSt, Qt))
  

  k = k + 1


QTave = QSSt.mean(0)
QTdev = QSSt.std(0)

sampleT = QSSt[0,:]

Sseries1=np.ones(len(Tave))
Sseries2=np.ones(len(QTave))

j = s0
for i in range(len(Sseries1)):
               Sseries1[i]=j
               j = j-1
j = si
for i in range(len(Sseries2)):
               Sseries2[i]=j
               j = j-1

n1=len(Sseries1)
n2=len(Stave)
print(n1,n2)
######################################################################################
############################ Plot Averages and Means #################################
######################################################################################

fig=plt.figure(1)
ax = plt.subplot(111, xlabel='$t$', ylabel='$n_P$')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
             item.set_fontsize(15),
#plt.step(t, S, 'b', data=None, where='pre',linewidth=0.25)
#plt.step(Tave, pSeries, '-k', data=None, where='post',linewidth=1.0)
#plt.step(Tave-Tdev, Pave, '-r', data=None, where='post',linewidth=1.0)
#plt.step(Tave+Tdev, Pave, '-b', data=None, where='post',linewidth=1.0)

plt.plot(Tave, pSeries, '-k', linewidth=1.0)
plt.plot(Tave+Tdev, pSeries, '-.k', linewidth=1.0)
plt.plot(Tave-Tdev, pSeries, '-.k', linewidth=1.0)

plt.plot(QTave, QP, '-r', linewidth=1.0)
plt.plot(QTave+QTdev, QP, '-.r', linewidth=1.0)
plt.plot(QTave-QTdev, QP, '-.r', linewidth=1.0)
             

#plt.plot(Tave, Plw, '-.r', linewidth=1.0) 
plt.step(sampleT, pSeries, 'b', data=None, where='post',linewidth=1.0)
#ax.set_xlim((0, 500.0))
#ax.set_ylim((0, 100.0))
plt.tight_layout()
plt.show()






