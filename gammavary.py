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

#dim=6

t=1.0
E=[0,0]

#min=-10; max=10; step=0.5
min=-2; max=2; step=0.01
#min=1.0; max=1.1; step=0.1
#t = np.linspace(min,max,1000)
    

# Output  files
f0 = open('E0_vs_gamma.dat', 'w')
f1 = open('E1_vs_gamma.dat', 'w')
U = float(input('Enter U\n'))
eps = float(input('Enter eps\n'))
l = float(input('Enter lambda\n'))

#gamma=0.0 # Reproduces single dissipative results
#gamma = complex(0,0.01)
# gamma=0.1
#U = 2.0
#eps = 0.0
#AbsE=np.zeros(6,float)

print('eps, lambda, U =',eps, l, U)

tm=t-l; tp=t+l

for gamma in np.arange(min,max,step):

    Emid=2.0*(eps-gamma)+U/2.0
    igamma=complex(0,gamma)
    #epsp=eps+gamma
    #epsm=eps-gamma
    epsp=eps+igamma
    epsm=eps-igamma
    #H=np.zeros((dim,dim),float)

    #H[0][0]=H[2][2]=H[3][3]=H[5][5]=2.0*eps
    #H[1][1]=H[4][4]=2.0*eps+U
    #H[1][2]=H[2][4]=tm
    #H[2][1]=H[4][2]=tp
    #H[1][3]=H[3][4]=-tm
    #H[3][1]=H[4][3]=-tp


    #E[0]=E[1]=E[2]=2*eps
    #E[0]=E[1]=E[2]=2*eps
    #E[3]=2*eps+U

    #S=epsm+epsm
    S=2.0*eps
    M=4*tp*tm
    Dsq=(epsp-epsm)**2
    K=Dsq+M
    L=Dsq*U

    def cubicroot(U,K,L):
         coeff=[1,-U,-K,L]
         return np.roots(coeff) 

    

    roots=cubicroot(U,K,L)
    print('type=',type(roots), 'cubic roots=',roots)
   
    #roots =  roots.tolist() 
    #print('type=',type(roots), 'cubic roots=',roots)
    
    print('roots=',roots)

    flag_imag = 0
    Emid =2.0	


    # Check if all the roots are real of not
    for x in roots:
        #if x.imag != 0:
        if abs(x.imag) >= 0.0001:
             flag_imag = 1 # Complex roots exist

 

    # Now treating some-complex and all-real roots differently
    if (flag_imag == 1): 
         z = [x for x in roots for y in roots if x != y and abs(abs(x) - abs(y)) < 0.001 ]
                                          # Selecting 
         z = np.array(z)
         z=S+U-z
    else:
         print('All real roots:',roots)
         print('S+U=',S+U)
         roots=S+U-roots
         print(roots[0],roots[1],roots[2])
         #z=S+U-z
         if abs(roots[0].real-Emid-(Emid-roots[1].real)) < 0.001:
             print('roots[0] and roots[1] picked up ...')
             z = np.array([roots[0],roots[1]])
         elif abs(roots[0].real-Emid-(Emid-roots[2].real)) < 0.001:
             print('roots[0] and roots[2] picked up ...')
             z = np.array([roots[0],roots[1]])
         elif abs(roots[1].real-Emid-(Emid-roots[2].real)) < 0.001:
             print('roots[1] and roots[2] picked up ...')
             z = np.array([roots[1],roots[2]])
         #z = complex(z.real,0)
         #z = [x for x in roots for y in roots if abs( (Emid - x.real) -(Emid - y.real) ) < 0.001 ] 
         #z = np.array(z)
         print('z=',z)
         z[0] = roots[0]; z[1] = roots[1] # Not guranteed 
                                          # Roots should be equally distant from Emid


    print('flag_imag =',flag_imag) 
   
    ##z=np.array(z) 
    ##z=z-U
    ##z=S-z

 
    print('z=',z)
    #exit()

    E[0]=z[0]
    E[1]=z[1]
    #E[4]=0.5*(4*eps+U-cmath.sqrt(16*tp*tm+U*U))
    #E[5]=0.5*(4*eps+U+cmath.sqrt(16*tp*tm+U*U))
    #E[4]=0.5*(4*eps-np.sqrt(16*tp*tm+U*U)+U)
    #E[5]=0.5*(4*eps+np.sqrt(16*tp*tm+U*U)+U) #Seems np.sqrt fails to get sqrt of -ve number
    #E,psi=np.linalg.eig(H)
   

 
    for i in range(1):
        # AbsE[i] = abs(E[i]) # Wrong syntax
        if (abs(E[i])<1.e-6):
            E[i]=0.0



    print("lambda=",l,"Eigenvalues=",E)

    #print("Expected eigen values:")
    #print(2.0*eps,"degeneracy=4")
    #print(2.0*(eps+np.sqrt(tp*tm)),"degeneracy=1")
    #print(2.0*(eps-np.sqrt(tp*tm)),"degeneracy=1\n")
    #print("Eigenvector matrix=\n",psi)
    #AbsE=E
    #AbsE=np.absolute(E)
    #AbsE=[abs(x) for x in E]

    #print("E=\n",E)
    #print("AbsE=\n",AbsE)


    # Sorting to new arrays 
    #index = E.argsort()
    #index = AbsE.argsort()
    #index = np.argsort(AbsE) # Alternative to above 
    #index = np.argsort(E) # Alternative to above 
    #print("Index=\n",index)
    #E = E[index]
    #psi = psi[:, index]
    #print("E after indexing =",E)
    #print("GS eigenvalue=",E[0],"\n")
    #print("E=\n",E)

    # Storing output data into files
    print ('{:4.2f} \t {:4.2f} \t{:4.2f}'.format(gamma,E[0].real,E[0].imag), file=f0)
    print ('{:4.2f} \t {:4.2f} \t{:4.2f}'.format(gamma,E[1].real,E[1].imag), file=f1)


f0.close()
f1.close()


# Plotting

# Set plot features 
fig, ax = plt.subplots()
#fig, ax = plt.subplots(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')



plt.figure(1)
ax.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
ax.tick_params(axis='both',which='minor',width=2,length=5) 
plt.xlabel('$\gamma$', size=bsize)
plt.ylabel('Re $E$', size=bsize)
data=np.loadtxt('E0_vs_gamma.dat')
ax.plot(data[:,0], data[:,1], 'bx', lw=2.0, label='$E=E^-$')
data=np.loadtxt('E1_vs_gamma.dat')
ax.plot(data[:,0], data[:,1], 'r+', lw=2.0, label='$E=E^+$')
plt.title('$t=${}, $U$={}, $\epsilon=${}'.format(t,U,eps), size=msize)
ax.legend(loc='upper right', prop=dict(size=18), fancybox=True, framealpha=0.8)
xannote=0.0; yannote=1.0
bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.8)
ax.text(xannote, yannote, "$\lambda=${}".format(l), ha="center", va="center", size=msize, bbox=bbox_props)
plt.gca().xaxis.set_major_locator(plt.MultipleLocator(0.8))
plt.tight_layout()
ax.set_rasterized(True)
plt.savefig('ReE_vs_gamma_fixed_lambda.png')
plt.savefig('ReE_vs_gamma_fixed_lambda.eps')
#
#
#
#
# --------------------------------------------------------------------------- 
#
#
#
#fig, ax = plt.subplots(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
fig, ax = plt.subplots()
#
#
plt.figure(2)
ax.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
ax.tick_params(axis='both',which='minor',width=2,length=5) 
plt.xlabel('$\gamma$', size=bsize)
plt.ylabel('Im $E$', size=bsize)
data=np.loadtxt('E0_vs_gamma.dat')
ax.plot(data[:,0],data[:,2], 'gx', label='$E=E^-$')
data=np.loadtxt('E1_vs_gamma.dat')
ax.plot(data[:,0],data[:,2], 'm+', label='$E=E^+$')
ax.legend(loc='upper right', prop=dict(size=18), fancybox=True, framealpha=0.8)
plt.title('$t=${}, $U$={}, $\epsilon=${}'.format(t,U,eps), size=msize)
bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.8)
xannote=0.0; yannote=3.0 
ax.text(xannote, yannote, "$\lambda=${}".format(l), ha="center", va="center", size=msize, bbox=bbox_props)
#bbox_props = dict(boxstyle="rarrow,pad=0.3", fc="cyan", ec="b", lw=2)
plt.gca().xaxis.set_major_locator(plt.MultipleLocator(0.8))
plt.tight_layout()
ax.set_rasterized(True)
plt.savefig('ImE_vs_gamma_fixed_lambda.png')
plt.savefig('ImE_vs_gamma_fixed_lambda.eps')



plt.show()

    



