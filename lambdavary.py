# Description: Solves cubic eigenenergy eq. for double PT symmetric Hubbard hydrogen molecule
# Varying paramter: lambda

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


#fig, ax = plt.subplots(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
#fig, ax = plt.subplots(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')


#dim=6

t=1.0
E=[0,0]

#min=-10; max=10; step=0.5
min=-1.3; max=1.3; step=0.005
#min=1.0; max=1.1; step=0.1
#t = np.linspace(min,max,1000)
    

# Output  files
f0 = open('E0_vs_lambda.dat', 'w')
f1 = open('E1_vs_lambda.dat', 'w')
U = float(input('Enter U\n'))
eps = float(input('Enter eps\n'))
gamma = float(input('Enter gamma\n'))

#gamma=0.0 # Reproduces single dissipative results
#gamma = complex(0,0.01)
#gamma=0.1
#U = 2.0
#eps = 0.0
#AbsE=np.zeros(6,float)
igamma=complex(0,gamma)
#epsp=eps+gamma
#epsm=eps-gamma
epsp=eps+igamma
epsm=eps-igamma

Emid=2.0*(eps-gamma)+U/2.0
print('eps, gamma, U, Emid =',eps, gamma, U, Emid)


iter=0; zprev=0
for l in np.arange(min,max,step):
    tm=t-l; tp=t+l
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
         print('Shifted roots:',roots[0],roots[1],roots[2])
         #z=S+U-z
                 #z = complex(z.real,0)
         #z = [x for x in roots for y in roots if abs( (Emid - x.real) -(Emid - y.real) ) < 0.001 ]       
         #z = np.array(z)
         roots[0] = complex(roots[0].real,0)
         roots[1] = complex(roots[1].real,0)
         z[0] = roots[0]; z[1] = roots[1] # Not guranteed 
         print('iter=',iter, 'zprev=',zprev)
         #if iter>0 and abs(abs(zprev.real)-abs(roots[0].real))>0.5:
          #    print('If condn satisfied, swapping ...') 
           #   z[0] = roots[1]; z[1] = roots[0]
                                          # Roots should be equally distant from Emid
         zprev=z[0]
         iter = iter + 1
    print('flag_imag =',flag_imag) 
    print('z=',z)
   
    ##z=np.array(z) 
    ##z=z-U
    ##z=S-z

 
    #exit()
     
    E[0]=z[0]
    E[1]=z[1]
    
    if E[0].imag > 0: # E[0].imag always has to be -ve, else swap
        E[0]=z[1]; E[1]=z[0]

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
    print('================ X ==============')
    print('')

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
    print ('{:4.3f} \t {:4.3f} \t{:4.3f}'.format(l,E[0].real,E[0].imag), file=f0)
    print ('{:4.3f} \t {:4.3f} \t{:4.3f}'.format(l,E[1].real,E[1].imag), file=f1)


f0.close()
f1.close()


# Plotting
#
#
# REAL PARTS
#
fig, ax = plt.subplots()
plt.figure(1) 
plt.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
plt.tick_params(axis='both',which='minor',width=2,length=5)
plt.xlabel('$\lambda$', size=bsize)
plt.ylabel('Re $E$', size=bsize)
data=np.loadtxt('E0_vs_lambda.dat')
plt.plot(data[:,0], data[:,1], 'bx', lw=2.0, label='$E=E^-$')
data=np.loadtxt('E1_vs_lambda.dat')
plt.plot(data[:,0], data[:,1], 'r+', lw=2.0, label='$E=E^+$')
plt.title('$t=${}, $U$={}, $\epsilon=${}'.format(t,U,eps), size=msize)
#ax.legend(loc='upper right', prop=dict(size=18), fancybox=True, framealpha=0.7, shadow=True)
ax.legend(loc='upper right', prop=dict(size=18), fancybox=True, framealpha=0.8)
# Annotate text
xannote=0.0; yannote=1.9
bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.8)
plt.text(xannote, yannote, "$\gamma=${}".format(gamma), ha="center", va="center", size=msize, bbox=bbox_props)
ax.annotate("$\lambda_{e3}$",
            xy=(-1.13, 1.98), xycoords='data',
            xytext=(0, 1), textcoords='data',
            size=20, va="center", ha="center",
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3,rad=-0.2"),
            )
ax.annotate("",
            xy=(1.13, 1.98), xycoords='data',
            xytext=(0.13, 1), textcoords='data',
            size=20, va="center", ha="center",
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3,rad=0.2"),
            )

ax.annotate("$\lambda_{e2}$",
            xy=(-1.09, 2.72), xycoords='data',
            xytext=(0, 2.5), textcoords='data',
            size=20, va="center", ha="center",
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3,rad=-0.05"),
            )
ax.annotate("",
            xy=(1.09, 2.72), xycoords='data',
            xytext=(0.13, 2.5), textcoords='data',
            size=20, va="center", ha="center",
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3,rad=0.05"),
            )
ax.annotate("$\lambda_{e1}$",
            xy=(-0.88, 3.15), xycoords='data',
            xytext=(0, 3.5), textcoords='data',
            size=20, va="center", ha="center",
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3,rad=0.1"),
            )
ax.annotate("",
            xy=(0.88, 3.15), xycoords='data',
            xytext=(0.13, 3.5), textcoords='data',
            size=20, va="center", ha="center",
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3,rad=-0.1"),
            )

plt.tight_layout()
ax.set_rasterized(True)
plt.savefig('ReE_vs_lambda_fixed_gamma.png')
plt.savefig('ReE_vs_lambda_fixed_gamma.eps')
#
#
#
# IMAGINARY PARTS
#
fig, ax = plt.subplots()
plt.figure(2) 
plt.tick_params(axis='both',which='major',width=2,length=10,labelsize=18)
plt.tick_params(axis='both',which='minor',width=2,length=5)
plt.xlabel('$\lambda$', size=bsize)
plt.ylabel('Im $E$', size=bsize)
data=np.loadtxt('E0_vs_lambda.dat')
plt.plot(data[:,0],data[:,2], 'gx', lw=2.0, label='$E=E^-$')
data=np.loadtxt('E1_vs_lambda.dat')
plt.plot(data[:,0],data[:,2], 'm+', lw=1.5, label='$E=E^+$')
plt.title('$t=${}, $U$={}, $\epsilon=${}'.format(t,U,eps), size=msize)
ax.legend(loc='upper right', prop=dict(size=18), fancybox=True, framealpha=0.8)
# Annotate text
xannote=0.0; yannote=0.25
bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.8) 
plt.text(xannote, yannote, "$\gamma=${}".format(gamma), ha="center", va="center", size=msize, bbox=bbox_props)
#ax.annotate( '$\lambda_{e1}$', xy=(-1, 0), xytext=(0,-2), ha='center', arrowprops={'shrink':0.05})
ax.annotate("$\lambda_{e3}$",
            xy=(-1.13, 0), xycoords='data',
            xytext=(0, -1), textcoords='data',
            size=20, va="center", ha="center",
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3,rad=-0.2"),
            )
ax.annotate("",
            xy=(1.13, 0), xycoords='data',
            xytext=(0.13, -1), textcoords='data',
            size=20, va="center", ha="center",
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3,rad=0.2"),
            )

ax.annotate("$\lambda_{e2}$",
            xy=(-1.09, 0), xycoords='data',
            xytext=(0, .75), textcoords='data',
            size=20, va="center", ha="center",
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3,rad=0.2"),
            )
ax.annotate("",
            xy=(1.09, 0), xycoords='data',
            xytext=(0.13, .75), textcoords='data',
            size=20, va="center", ha="center",
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3,rad=-0.2"),
            )

ax.annotate("$\lambda_{e1}$",
            xy=(-0.88, 0), xycoords='data',
            xytext=(0, 1), textcoords='data',
            size=20, va="center", ha="center",
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3,rad=0.4"),
            )
ax.annotate("",
            xy=(0.88, 0), xycoords='data',
            xytext=(0.13, 1), textcoords='data',
            size=20, va="center", ha="center",
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3,rad=-0.4"),
            )


plt.tight_layout()
ax.set_rasterized(True)
plt.savefig('ImE_vs_lambda_fixed_gamma.png')
plt.savefig('ImE_vs_lambda_fixed_gamma.eps')

plt.show()

    



