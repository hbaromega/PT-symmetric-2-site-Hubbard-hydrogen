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
f1 = open('lambda_e_vs_U.dat', 'w')
min=-10; max=10; step=0.01
for U in np.arange(min,max,step):
    lambe=0.25*np.sqrt(16*t*t+U*U)
    print ('{:4.2f} \t {:4.2f} '.format(U,lambe), file=f1)
f1.close()

plt.xlim([-10, 10])
plt.ylim(0, 3)

plt.figure(1)
plt.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
plt.tick_params(axis='both',which='minor',width=2,length=5)  
plt.xlabel('$U$', size=bsize)
plt.ylabel('$|\lambda|$', size=bsize)
data=np.loadtxt('lambda_e_vs_U.dat')
x=data[:,0]; y=data[:,1]
plt.plot(x, y, 'b-', lw=4.0)
plt.title('$t=${}, $\epsilon=${}'.format(t,eps), size=msize)
#plt.legend()
# Coloring bg
background_color = 'k'
plt.fill_between(x, y, where=y<=y, interpolate=True, color='grey')
plt.annotate('$PT$ broken', xy=(-1.5,2), fontsize=20) 
plt.annotate('$PT$ unbroken', xy=(-1.5,0.5), fontsize=20)
plt.annotate('$|\lambda_e|$', xy=(-8.62, 2.37), xytext=(-7.2, 2.5),
            #arrowprops=dict(facecolor='black', shrink=0.05),  fontsize=18) 
            arrowprops=dict(arrowstyle='->', lw=1.5, facecolor='black'),  fontsize=18) 
plt.tight_layout()
plt.savefig('lambda_e_vs_U.png')
plt.savefig('lambda_e_vs_U.eps')
plt.show()
