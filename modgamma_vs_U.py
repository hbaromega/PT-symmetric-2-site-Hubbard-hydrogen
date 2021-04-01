import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
#from numpy import linalg
import os
import subprocess
from pathlib import Path
import cmath

# Fontsize
bsize=20
msize=18

#t = float(input('Enter t\n'))
t=1.0
eps=0.5

plt.xlim([0, 10.1])
plt.ylim(0, 1.05)

style = ['ro-','b*--','g-']
col = ['grey', 'lightgrey']
m = [8,12]
plt.figure(1)
plt.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
plt.tick_params(axis='both',which='minor',width=2,length=5)  
plt.xlabel('$U$', size=bsize)
plt.ylabel('$|\gamma|$', size=bsize)

i=0
for l in [0, 0.5]:
     print('l=',l)
     data=np.loadtxt('L{}/gamma_e_vs_U.data'.format(l))
     x=data[:,0]; y=data[:,1]
     plt.plot(x, y, style[i], ms=m[i], lw=4.0)
     plt.fill_between(x, y, where=y<=y, interpolate=True, color=col[i])
     i += 1
plt.title('$t=${}, $\epsilon=${}'.format(t,eps), size=msize)
#plt.legend()
# Coloring bg
background_color = 'k'
plt.annotate('$PT$ broken', xy=(6,0.8), fontsize=19) 
plt.annotate('$PT$ unbroken', xy=(1,0.075), fontsize=19)
plt.annotate('$|\gamma_e|(\lambda=0)$', xy=(0.49, 0.68), xytext=(1.53, 0.84),
            arrowprops=dict(arrowstyle='->', lw=1.5, facecolor='black'),  fontsize=18) 
plt.annotate('$|\gamma_e|(\lambda=0.5)$', xy=(2, 0.29), xytext=(2.89, 0.37),
            #arrowprops=dict(facecolor='black', shrink=0.05),  fontsize=18) 
            arrowprops=dict(arrowstyle='->', lw=1.5, facecolor='black'),  fontsize=18) 
plt.tight_layout()
plt.savefig('gamma_e_vs_U.png')
plt.savefig('gamma_e_vs_U.eps')
plt.show()
