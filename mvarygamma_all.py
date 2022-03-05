# Description: Solves cubic eigenenergy eq. for double PT symmetric Hubbard hydrogen molecule
# Varying paramter: gamma

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
#fig, ax = ax.subplots(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
#fig1, ax1 = plt.subplots(figsize=(8, 6))
#fig2, ax2 = plt.subplots(figsize=(8, 6))
fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

#dim=6

t=1.0
E=[0,0,0]
tol = 0.0001
#min=-10; max=10; step=0.5
min=-2; max=2; step=0.01
#min=1.0; max=1.1; step=0.1
#t = np.linspace(min,max,1000)
    
# Inputs
U = float(input('Enter U\n'))
eps = float(input('Enter eps\n'))
l = float(input('Enter lambda\n'))

flag = 1

if flag == 0:
   U = 0; l = 0.6; npoints = 501; tag=0
else:
   U = 2.0; l = 0; npoints = 301; tag=2

#gamma=0.0 # Reproduces single dissipative results
#gamma = complex(0,0.01)
# gamma=0.1
#U = 2.0
#eps = 0.0
#AbsE=np.zeros(6,float)

f2 = open('reroots_vs_gamma.dat', 'w')
f3 = open('imroots_vs_gamma.dat', 'w')

print('eps, lambda, U =',eps, l, U)

tm=t-l; tp=t+l
S=2.0*eps
M=4*tp*tm

iter = 0

def cubicroot(U,K,L):
         coeff=[1,-U,-K,L]
         return np.roots(coeff) 

for gamma in np.linspace(-1.5, 1.5, npoints): 


    Emid=2.0*(eps-gamma)+U/2.0
    igamma=complex(0,gamma)
    epsp=eps+igamma
    epsm=eps-igamma
    Dsq=(epsp-epsm)**2
    K=Dsq+M
    L=Dsq*U
    
    roots=cubicroot(U,K,L)
    roots=S+U-roots 
    print('type=',type(roots), 'cubic roots=',roots)
    print('roots=',roots)

    # Check if all the roots are real of not
    for x in roots:
        #if x.imag != 0:
        if abs(x.imag) >= 0.0001:
             flag_imag = 1 # Complex roots exist


    #if (flag_imag == 1):  
    if abs(roots[0].real - roots[1].real) < tol:
       E[0]=roots[0]
       E[1]=roots[1]
       E[2]=roots[2]
    elif abs(roots[0].real - roots[2].real) < tol:
       E[0]=roots[0]
       E[1]=roots[2]
       E[2]=roots[1]
    else:
       E[0]=roots[1]
       E[1]=roots[2]
       E[2]=roots[0]
 
    for i in range(2):
        # AbsE[i] = abs(E[i]) # Wrong syntax
        if (abs(E[i])<1.e-6):
            E[i]=0.0


    print ('{:4.2f} \t {:4.2f} \t{:4.2f} \t{:4.2f}'.format(gamma,E[0].real,E[1].real,E[2].real), file=f2)
    print ('{:4.2f} \t {:4.2f} \t{:4.2f} \t{:4.2f}'.format(gamma,E[0].imag,E[1].imag,E[2].imag), file=f3)


    if flag == 0:
      ReE1 = E[2].real # b
      ReE2 = E[0].real # r
      ReE3 = E[1].real # c

      ImE1 = E[0].imag
      ImE2 = E[1].imag
      ImE3 = E[2].imag

    else: 
      ReE1 = E[2].real # b
      ReE2 = E[0].real # r
      ReE3 = E[1].real # c

      ImE1 = E[2].imag
      ImE2 = E[1].imag
      ImE3 = E[0].imag


    print("lambda=",l,"Eigenvalues=",E)

  
    # Storing output data into files
    #print ('{:4.2f} \t {:4.2f} \t{:4.2f}'.format(gamma,E[0].real,E[0].imag), file=f0)
    #print ('{:4.2f} \t {:4.2f} \t{:4.2f}'.format(gamma,E[1].real,E[1].imag), file=f1)
    ax1.plot(gamma, ReE1, 'bx', lw=2.0, label='$E=E^-$' if iter == 0 else '')
    ax1.plot(gamma, ReE2, 'r+', lw=2.0, label='$E=E^+$' if iter == 0 else '')
    ax1.plot(gamma, ReE3, 'c.', lw=2.0, label='$E=E^0$' if iter == 0 else '')

    ax2.plot(gamma, ImE1, 'gx', label='$E=E^-$' if iter == 0 else '')
    ax2.plot(gamma, ImE2, 'm+', label='$E=E^+$' if iter == 0 else '')
    ax2.plot(gamma, ImE3, 'orange', marker='.', ls='None', label='$E=E^0$' if iter == 0 else '')


    #iter += 1
    iter = 1 # We make iter=1 since only iter=0 recorded for legends


f2.close()
f3.close()

                                          # Roots should be equally distant from Emid
'''
    print ('{:4.2f} \t {:4.2f} \t{:4.2f} \t{:4.2f}'.format(gamma,roots[0].real,roots[1].real,roots[2].real), file=f2)
    print ('{:4.2f} \t {:4.2f} \t{:4.2f} \t{:4.2f}'.format(gamma,roots[0].imag,roots[1].imag,roots[2].imag), file=f3)
'''


# Plot formatting

ax1.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
ax1.tick_params(axis='both',which='minor',width=2,length=5) 
ax1.set_xlabel('$\gamma$', size=bsize)
ax1.set_ylabel('Re $E$', size=bsize)
ax1.set_title('$t=${}, $U$={}, $\epsilon=${}'.format(t,U,eps), size=msize)
ax1.legend(loc='upper right', prop=dict(size=18), fancybox=True, framealpha=0.8)
xannote=0.0; yannote=2.0
bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.8)
ax1.text(xannote, yannote, "$\lambda=${}".format(l), ha="center", va="center", size=msize, bbox=bbox_props)
fig1.gca().xaxis.set_major_locator(plt.MultipleLocator(0.8))
#fig1.tight_layout()
ax1.set_rasterized(True)
# Adjust layout
fig1.subplots_adjust(left=0.15, right=0.95, bottom=0.2,  top=0.9)
fig1.savefig('ReE_vs_gamma_fixed_lambda_Ueq{}.png'.format(tag))
fig1.savefig('ReE_vs_gamma_fixed_lambda_Ueq{}.pdf'.format(tag))
fig1.savefig('ReE_vs_gamma_fixed_lambda_Ueq{}.eps'.format(tag))

#
ax2.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
ax2.tick_params(axis='both',which='minor',width=2,length=5) 
ax2.set_xlabel('$\gamma$', size=bsize)
ax2.set_ylabel('Im $E$', size=bsize)
ax2.set_title('$t=${}, $U$={}, $\epsilon=${}'.format(t,U,eps), size=msize)
ax2.legend(loc='upper right', prop=dict(size=18), fancybox=True, framealpha=0.8)
xannote=0.0; yannote=2.0
bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.8)
ax2.text(xannote, yannote, "$\lambda=${}".format(l), ha="center", va="center", size=msize, bbox=bbox_props)
fig1.gca().xaxis.set_major_locator(plt.MultipleLocator(0.8))
#fig1.tight_layout()
ax2.set_rasterized(True)
# Adjust layout
fig2.subplots_adjust(left=0.15, right=0.95, bottom=0.2,  top=0.9)
fig2.savefig('ImE_vs_gamma_fixed_lambda_Ueq{}.png'.format(tag))
fig2.savefig('ImE_vs_gamma_fixed_lambda_Ueq{}.pdf'.format(tag))
fig2.savefig('ImE_vs_gamma_fixed_lambda_Ueq{}.eps'.format(tag))
# Check: https://stackoverflow.com/questions/18619880/matplotlib-adjust-figure-margin

plt.show()

    



