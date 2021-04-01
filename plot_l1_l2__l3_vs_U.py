import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
#from numpy import linalg
import os
import subprocess
from pathlib import Path
import cmath


#mydir = '/Users/hbar/Documents/'
#mydir = subprocess.check_output("echo $PWD")
# mydir = os.path.dirname(os.path.abspath(__file__)) # works
mydir = os.getcwd()
print('mydir =',mydir)

# Fontsize
bsize=20
msize=18

t=1.0
eps=0.5



style=['ro','b*','gx']

# Plot EP1 vs U 
fig, ax = plt.subplots()
plt.figure(1)
i=0
ax.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
ax.tick_params(axis='both',which='minor',width=2,length=5) 
plt.xlabel('$U$', size=bsize)
plt.ylabel('$|\lambda_{e1}|$', size=bsize)
plt.title('$t=${}, $\epsilon=${}'.format(t,eps), size=msize)
ax.set_xlim(0,6)
for gamma in [0.1, 0.2]:
    print('gamma=',gamma)
    data=np.loadtxt('G{}/l1_vs_U.data'.format(gamma))
    ax.plot(data[:,0], data[:,1], style[i], lw=2.0, label='$\gamma=${}'.format(gamma))
    i += 1
ax.legend(loc='upper right', prop=dict(size=18), fancybox=True, framealpha=0.8)
#xannote=0.0; yannote=1.0
#bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.8)
#ax.text(xannote, yannote, "$\lambda=${}".format(l), ha="center", va="center", size=msize, bbox=bbox_props)
plt.tight_layout()
ax.set_rasterized(True)
plt.savefig('lamdae1_vs_U.png')
plt.savefig('lamdae1_vs_U.eps')
#
#

# Plot EP2 and EP'2 vs U 
fig, ax = plt.subplots()
plt.figure(2)
ax.set_xlim(0,6)
i=0
ax.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
ax.tick_params(axis='both',which='minor',width=2,length=5) 
plt.xlabel('$U$', size=bsize)
plt.ylabel('$|\lambda_{e2}|$', size=bsize)
plt.title('$t=${}, $\epsilon=${}'.format(t,eps), size=msize)
#for gamma in [0.1, 0.2, 0.3]:
for gamma in [0.1, 0.2]:
    data=np.loadtxt('G{}/l2_vs_U.data'.format(gamma))
    ax.plot(data[:,0], data[:,1], style[i], lw=2.0, label='$\gamma=${}'.format(gamma))
    i += 1
ax.legend(loc='upper left', prop=dict(size=18), fancybox=True, framealpha=0.8)
#xannote=0.0; yannote=1.0
#bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.8)
#ax.text(xannote, yannote, "$\lambda=${}".format(l), ha="center", va="center", size=msize, bbox=bbox_props)
plt.tight_layout()
ax.set_rasterized(True)
plt.savefig('lamdae2_vs_U.png')
plt.savefig('lamdae2_vs_U.eps')

fig, ax = plt.subplots()
plt.figure(3)
ax.set_xlim(0,6)
i=0
ax.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
ax.tick_params(axis='both',which='minor',width=2,length=5) 
plt.xlabel('$U$', size=bsize)
plt.ylabel('$|\lambda_{e3}|$', size=bsize)
plt.title('$t=${}, $\epsilon=${}'.format(t,eps), size=msize)
#for gamma in [0.1, 0.2, 0.3]:
for gamma in [0.1, 0.2]:
    data=np.loadtxt('G{}/l3_vs_U.data'.format(gamma))
    ax.plot(data[:,0], data[:,1], style[i], lw=2.0, label='$\gamma=${}'.format(gamma))
    i += 1
ax.legend(loc='upper left', prop=dict(size=18), fancybox=True, framealpha=0.8)
#xannote=0.0; yannote=1.0
#bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.8)
#ax.text(xannote, yannote, "$\lambda=${}".format(l), ha="center", va="center", size=msize, bbox=bbox_props)
plt.tight_layout()
ax.set_rasterized(True)
plt.savefig('lamdae3_vs_U.png')
plt.savefig('lamdae3_vs_U.eps')


plt.show()

    



