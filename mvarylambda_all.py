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


#mydir = '/Users/hbar/Documents/'
#mydir = subprocess.check_output("echo $PWD")
# mydir = os.path.dirname(os.path.abspath(__file__)) # works
mydir = os.getcwd()
print('mydir =',mydir)

# Fontsize
bsize=20
msize=18
#fig, ax = ax.subplots(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()


t=1.0
E=[0,0,0]

#min=-10; max=10; step=0.5
min=-2; max=2; step=0.01
#min=1.0; max=1.1; step=0.1
#t = np.linspace(min,max,1000)
    
# Inputs
U = float(input('Enter U\n'))
eps = float(input('Enter eps\n'))
gamma = float(input('Enter gamma\n'))

f2 = open('reroots_vs_lambda.dat', 'w')
f3 = open('imroots_vs_lambda.dat', 'w')

print('eps, gamma, U =',eps, gamma, U)

tol = 0.0001
gamma_c = 0.19

iter = 0


def cubicroot(U,K,L):
         coeff=[1,-U,-K,L]
         return np.roots(coeff) 

def swap(var1,var2):
          tmp = var1
          var1 = var2 
          var2 = tmp
          return var1, var2


igamma=complex(0,gamma)
epsp=eps+igamma
epsm=eps-igamma
Emid=2.0*(eps-gamma)+U/2.0
S=2.0*eps
Dsq=(epsp-epsm)**2
L=Dsq*U
Emid =2.0	

flag_imag = 0



for l in np.arange(min,max,step):


    tm=t-l; tp=t+l
    M=4*tp*tm
    K=Dsq+M

    roots=cubicroot(U,K,L)
    roots=S+U-roots
    print('type=',type(roots), 'cubic roots=',roots)
    print('roots=',roots)



    # Check if all the roots are real of not
    for x in roots:
        #if x.imag != 0:
        if abs(x.imag) >= 0.0001:
             flag_imag = 1 # Complex roots exist

     
   

    print ('{:4.2f} \t {:4.2f} \t{:4.2f} \t{:4.2f}'.format(l,roots[0].real,roots[1].real,roots[2].real), file=f2)
    print ('{:4.2f} \t {:4.2f} \t{:4.2f} \t{:4.2f}'.format(l,roots[0].imag,roots[1].imag,roots[2].imag), file=f3)



    print('flag_imag =',flag_imag) 
   
    for i in range(1):
        # AbsE[i] = abs(E[i]) # Wrong syntax
        if (abs(E[i])<1.e-6):
            E[i]=0.0


    if gamma < gamma_c:


       if flag_imag == 0: # No complex roots, before hitting EP 
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
       else:  # complex roots arise
          E[0]=roots[0] # b
          E[1]=roots[1] # r
          E[2]=roots[2] # c
          flag_imag = 1

       if E[0].real > E[1].real: # If blue line comes over red line
          # Swap
        ctmp  = E[0] 
        E[0] = E[1]
        E[1] = ctmp

       
       ReE1 = E[0].real # b
       ReE2 = E[1].real # r
       ReE3 = E[2].real # c



    else: #  gam >= gamma_c:
      

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

       ReE1 = E[2].real # b
       ReE2 = E[0].real # r
       ReE3 = E[1].real # c

    
    ImE1 = E[0].imag
    ImE2 = E[1].imag
    ImE3 = E[2].imag

    # Color convention
    # Im E1 --> green (E^-)
    # Im E2 ---> magenta (E^+) [always positive]
    # Im E3 ----> orange (E^0)

    # Make sure ImE2 always positive definite
    if ImE1 > 0 and ImE2 < 0 and abs(ImE3) < tol:
          # Swap ImE1 & ImE2
          ImE1, ImE2 = swap(ImE1,ImE2)
      
    if ImE3 > 0 and ImE2 < 0 and abs(ImE1) < tol:
          # Swap ImE3 & ImE2
          ImE3, ImE2 = swap(ImE3,ImE2)

    if gamma >= gamma_c: 
          # locate conjugate pairs # swap ImE1 & ImE3
          if ImE2 < tol and ImE3 > 0:
             # Swap ImE3 & ImE2
             ImE3, ImE2 = swap(ImE3,ImE2)
          if abs(ImE3) < tol and ImE1 < 0: # Make ImE3 always negative
             # Swap ImE3 & ImE1
             ImE3, ImE1 = swap(ImE3,ImE1)


    print("lambda=",l,"Eigenvalues=",E)

  
    # Storing output data into files
    #print ('{:4.2f} \t {:4.2f} \t{:4.2f}'.format(l,E[0].real,E[0].imag), file=f0)
    #print ('{:4.2f} \t {:4.2f} \t{:4.2f}'.format(l,E[1].real,E[1].imag), file=f1)
    ax1.plot(l, ReE1, 'bx', lw=2.0, label='$E=E^-$' if iter == 0 else '')
    ax1.plot(l, ReE2, 'r+', lw=2.0, label='$E=E^+$' if iter == 0 else '')
    ax1.plot(l, ReE3, 'c*', lw=2.0, label='$E=E^0$' if iter == 0 else '')

    ax2.plot(l, ImE1, 'gx', label='$E=E^-$' if iter == 0 else '')
    ax2.plot(l, ImE2, 'm+', label='$E=E^+$' if iter == 0 else '')
    ax2.plot(l, ImE3, 'orange', marker='*', ls='None', label='$E=E^0$' if iter == 0 else '')


    #iter += 1
    iter = 1 # We make iter=1 since only iter=0 recorded for legends


f2.close()
f3.close()


# Plot formatting

ax1.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
ax1.tick_params(axis='both',which='minor',width=2,length=5) 
ax1.set_xlabel('$\lambda$', size=bsize)
ax1.set_ylabel('Re $E$', size=bsize)
ax1.set_title('$t=${}, $U$={}, $\epsilon=${}'.format(t,U,eps), size=msize)
ax1.legend(loc='upper right', prop=dict(size=18), fancybox=True, framealpha=0.8)
xannote=0.0; yannote=1.0
bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.8)
ax1.text(xannote, yannote, "$\gamma=${}".format(gamma), ha="center", va="center", size=msize, bbox=bbox_props)
fig1.gca().xaxis.set_major_locator(plt.MultipleLocator(0.8))
fig1.tight_layout()
ax1.set_rasterized(True)
fig1.savefig('ReE_vs_lambda_fixed_gamma.png')
fig1.savefig('ReE_vs_lambda_fixed_gamma.eps')

#
ax2.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
ax2.tick_params(axis='both',which='minor',width=2,length=5) 
ax2.set_xlabel('$\lambda$', size=bsize)
ax2.set_ylabel('Im $E$', size=bsize)
ax2.set_title('$t=${}, $U$={}, $\epsilon=${}'.format(t,U,eps), size=msize)
ax2.legend(loc='upper right', prop=dict(size=18), fancybox=True, framealpha=0.8)
xannote=0.0; yannote=1.0
bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.8)
ax2.text(xannote, yannote, "$\gamma=${}".format(gamma), ha="center", va="center", size=msize, bbox=bbox_props)
fig1.gca().xaxis.set_major_locator(plt.MultipleLocator(0.8))
fig1.tight_layout()
ax2.set_rasterized(True)
fig2.savefig('ImE_vs_lambda_fixed_gamma.png')
fig2.savefig('ImE_vs_lambda_fixed_gamma.eps')

# Adjust layout
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15,  top=0.95)
# Check: https://stackoverflow.com/questions/18619880/matplotlib-adjust-figure-margin

plt.show()

    



