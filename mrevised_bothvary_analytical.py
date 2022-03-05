# Description: Solves cubic eigenenergy eq. for double PT symmetric Hubbard hydrogen molecule
#   and surface-plots (3D) against non-Hermiticity parameters lambda and gamma
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable 
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
my_cmap=cm.get_cmap('autumn')

# Input parameters 
#U = float(input('Enter U\n'))
#eps = float(input('Enter eps\n'))

t = 1.0
eps = 0.5

flag_surf = 0
flagU = 0
if flagU ==0: #flagU = 0 # Set 0 when U = 0
   U = 0
else:
   U = 2.0

print('eps, U =',eps, U)

# Plot setup
bsize=26
msize=22
#fig, ax = ax.subplots(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
#fig1, ax1 = plt.subplots()
#fig2, ax2 = plt.subplots()
fig1 = plt.figure(0, figsize=(10, 8))
ax1 = fig1.add_subplot(111, projection='3d')
fig2 = plt.figure(1, figsize=(10, 8))
ax2 = fig2.add_subplot(111, projection='3d')
#fig = plt.figure(figsize=(6,6))


#elev = 30; azim = -60
#elev1 = 6; azim1 = 160
if flagU == 0: 
   azim1 = 142; elev1 = 9
   azim2 = 120; elev2 = 10 
else:
   azim1 = 172; elev1 = 33
   azim1 = 126; elev1 = 19
   #azim1 = 101; elev1 = 18 #**
   #azim1 = 165; elev1 = 28
   #azim2 = 139; elev2 = 29
   azim2 = 161; elev2 = 9
   azim2 = 107; elev2 = 14
   azim2 = 102; elev2 = 7
   azim2 = 104; elev2 = 2


#E=[0,0]
#E=np.zeros(2, dtype='complex_')
#E=np.zeros(2, dtype=np.complex_)
E = [ 0.+0.j,  0.+0.j, 0.+0.j]
#min=-10; max=10; step=0.5
lmin=-1.5; lmax=1.5; lstep=.2
gmin=-1.5; gmax=1.5; gstep=.2
pnum = 100   



iter = 0
gamma_c = 0.19
tol = 0.0001

# Create lambda-gamma mesh
#l = np.arange(lmin,lmax+lstep,lstep, float) # lambda array
#g = np.arange(gmin,gmax+gstep,gstep) # gamma array
l = np.linspace(lmin,lmax,41) # lambda array
g = np.linspace(gmin,gmax,41) # gamma array
Lx = len(l); Ly = len(g)
R0 = np.zeros( (Lx,Ly) )
R1 = np.zeros( (Lx,Ly) )
R2 = np.zeros( (Lx,Ly) )
I0 = np.zeros( (Lx,Ly) )
I1 = np.zeros( (Lx,Ly) )
I2 = np.zeros( (Lx,Ly) )
Lam, Gam = np.meshgrid(l,g)
print('l=\n',l)
print('g=\n',g)

ReEmid=2*eps+0.5*U
print('ReEmid=',ReEmid)
flag_imag = 0
tol = 0.0001

# Cubic root function
# -------------------
def cubicroot(U,K,L):
         coeff=[1,-U,-K,L]
         return np.roots(coeff) 

# Swapping function
# -------------------
def swap(var1,var2):
          tmp = var1
          var1 = var2 
          var2 = tmp
          return var1, var2




for i_l, lval in  enumerate(l): # lambda loop 
  tm=t-lval; tp=t+lval

  for i_g, gamma in enumerate(g): # gamma loop


    # Parameters
    # -----------
    Emid=2.0*(eps-gamma)+U/2.0
    igamma=complex(0,gamma)
    #epsp=eps+gamma
    #epsm=eps-gamma
    epsp=eps+igamma
    epsm=eps-igamma
    #S=epsm+epsm
    S=2.0*eps
    M=4*tp*tm
    Dsq=(epsp-epsm)**2
    K=Dsq+M
    L=Dsq*U

    roots=cubicroot(U,K,L)
    roots=S+U-roots
    #coeff=[1,-U,-K,L]
    #roots=S+U-np.roots(coeff)
    #global flag_imag

    print('type=',type(roots), 'cubic roots=',roots)
    print('roots=',roots)



    # Check if all the roots are real of not
    for x in roots:
        #if x.imag != 0:
        if abs(x.imag) >= 0.0001:
             flag_imag = 1 # Complex roots exist


    print('flag_imag =',flag_imag) 
  

    for i in range(2):
        # AbsE[i] = abs(E[i]) # Wrong syntax
        if (abs(roots[i])<1.e-6):
            roots[i]=0.0
 



    if abs(gamma) < gamma_c:


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



    # IMAGINARY PARTS

    if flag_imag == 0: # No complex roots, before hitting EP 
       if abs(roots[0].real - roots[1].real < tol): 
          E[0]=roots[0] 
          E[1]=roots[1]
          E[2]=roots[2]
       elif abs(roots[0].real - roots[2].real < tol):
          E[0]=roots[0]
          E[1]=roots[2]
          E[2]=roots[1]
       else:
          E[0]=roots[1]
          E[1]=roots[2]
          E[2]=roots[0]
    else:  # complex roots arise
         #if roots[1].imag >= 0.0: 
          E[0]=roots[0] # g
          E[1]=roots[1] # m
         #else:
          E[2]=roots[2] # orange


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

    if flagU != 0 and abs(gamma) >= gamma_c: 
          # locate conjugate pairs # swap ImE1 & ImE3
          if ImE2 < tol and ImE3 > tol:
             # Swap ImE3 & ImE2
             ImE3, ImE2 = swap(ImE3,ImE2)
          if abs(ImE3) < tol and ImE1 < tol: # ImE3 should be always -ve
             # Swap ImE1 & ImE3
             ImE1, ImE3 = swap(ImE1,ImE3)


     
   


    #print("lambda=",lval,"Eigenvalues=",E)
    R0[i_l][i_g] = ReE1 # Re E-
    R1[i_l][i_g] = ReE2 # Re E+
    R2[i_l][i_g] = ReE3 # Re E0

    I0[i_l][i_g] = ImE1 # Im E+
    I1[i_l][i_g] = ImE2 # Im E0
    I2[i_l][i_g] = ImE3 # Im E-

  

    iter += 1

print('R0 =\n',R0)
print('Types of Lambda, Gamma, R0 =',type(Lam), type(Gam), type(R0))
#exit()


# PLOTTING
# ***************************************************** 


# REAL PARTS PLOTTING
# =====================

# ---- Color map --------
#mycmap = plt.get_cmap('gist_earth')
mycmap1 = ('winter')
mycmap2 = ('autumn')
alpha1=0.4
alpha2=0.2
alpha3=0.3
# Get rid of colored axes planes
# First remove fill
ax1.xaxis.pane.fill = False
ax1.yaxis.pane.fill = False
ax1.zaxis.pane.fill = False
# Now set color to white (or whatever is "invisible")
ax1.xaxis.pane.set_edgecolor('w')
ax1.yaxis.pane.set_edgecolor('w')
ax1.zaxis.pane.set_edgecolor('w')
# --- 3D view setup -----
ax1.view_init(elev1,azim1) # elevation, azimuth
#
# Surface plots -------------->
#
if flag_surf == 1:
  surf1  = ax1.plot_surface(Gam, Lam, R0, color='b', rstride=1, cstride=1, linewidth=0, antialiased=False, alpha=alpha1)
  surf2  = ax1.plot_surface(Gam, Lam, R1, color='r', rstride=1, cstride=1, linewidth=0, antialiased=False, alpha=alpha2)
  surf3  = ax1.plot_surface(Gam, Lam, R2, color='c', rstride=1, cstride=1, linewidth=0, antialiased=False, alpha=alpha3)
else:
  surf1  = ax1.plot_wireframe(Gam, Lam, R0, color='b',  rstride=1, cstride=1, alpha=alpha1)
  surf2  = ax1.plot_wireframe(Gam, Lam, R1, color='r',  alpha=alpha2)
  surf3  = ax1.plot_wireframe(Gam, Lam, R2, color='c', rstride=1, cstride=1, alpha=alpha3)
# 
# Scatter plot ----------------->
for i_g,gval in enumerate(g):
     if gval > 0.1 and gval < 0.2:
     #if gval == 0.1: 
           scat = ax1.scatter(l, g[i_g], R0[:,i_g], c='b', marker='x', lw=2.0, zorder=4)
           scat = ax1.scatter(l, g[i_g], R1[:,i_g], c='r', marker='+', lw=2.0, zorder=4)
           scat = ax1.scatter(l, g[i_g], R2[:,i_g], c='c', marker='o', lw=2.0, zorder=4)
#
# Labeling etc 
ax1.tick_params(axis='both',which='major',width=2,length=10,labelsize=16)
#ax1.tick_params(axis='both',which='minor',width=2,length=5) 
#ax1.set_xlabel('$\lambda$', size=bsize, linespacing=3.4)
ax1.set_xlabel('$\lambda$', size=bsize)
ax1.xaxis.labelpad=10
ax1.set_ylabel('$\gamma$', size=bsize)
ax1.yaxis.labelpad=10
ax1.set_zlabel('Re $E$', size=bsize)
ax1.zaxis.labelpad=10
ax1.set_title('$t=${}, $U$={}, $\epsilon=${}'.format(t,U,eps), size=msize, y=0.95 )
#ax1.text2D(0.05, 0.95, 'E^+', transform=ax1.transAxes)
if flagU == 0: 
  #ax1.text(-1.98,-0.3,ReEmid-0.25, 'Re $E$ = {}'.format(ReEmid), color='orange', size=16, zorder=1)  
  ax1.text(1,-1.6,-1.45, '$E^-$', color='b', size=msize, zorder=1)
  ax1.text(1.6,-1.6,2.65, '$E^+$', color='r', size=msize)
  ax1.text(0.47,-1.7,-0.09, '$E^0$', color='c', size=msize, zorder=1)
else: 
  ax1.text(-1.7,1.6,1.25, '$E^-$', color='b', size=msize, zorder=1)
  ax1.text(1.7,-1.6,3, '$E^+$', color='r', size=msize)
  ax1.text(-2,0.9,4.8, '$E^0$', color='c', size=msize, zorder=1)
#ax1.legend(loc='upper right', prop=dict(size=18), fancybox=True, framealpha=0.8)
#xannote=0.0; yannote=1.0
bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.8)
#ax1.text(xannote, yannote, "$\lambda=${}".format(l), ha="center", va="center", size=msize, bbox=bbox_props)
#
# --- Set tick frequency -----
#fig1.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
#fig1.gca().yaxis.set_major_locator(plt.MultipleLocator(1))
# also
#ax1.xaxis.set_ticks(np.arange(lmin, lmax, 1.0))
#ax1.yaxis.set_ticks(np.arange(gmin, gmax, 1.0))
# also
ax1.xaxis.set_major_locator(plt.MultipleLocator(1))
ax1.yaxis.set_major_locator(plt.MultipleLocator(1))
#
# --- Format ticking numbers
#ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
#ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
   # Format ax.xaxis.set_ticks(np.arange(start, end, stepsize))
fig1.tight_layout()
ax1.set_rasterized(True)
fig1.savefig('3D_ReE_U{}.png'.format(U))
fig1.savefig('3D_ReE_U{}.eps'.format(U))
fig1.savefig('3D_ReE_U{}.pdf'.format(U))

L = Lx*Ly
col3 = np.arange(L)/L
cmap = plt.get_cmap('rainbow')
col = [cmap(i) for i in np.linspace(0, 100, L)]


# IMAGINARY PARTS PLOTTING
# ========================= 
mcol = 'xkcd:chocolate'
mcol = 'orange'
# ---- Color map --------
#mycmap = plt.get_cmap('gist_earth')
mycmap1 = ('summer')
mycmap2 = ('spring')
# --- 3D view setup -----
# Get rid of colored axes planes
# First remove fill
ax2.xaxis.pane.fill = False
ax2.yaxis.pane.fill = False
ax2.zaxis.pane.fill = False
# Now set color to white (or whatever is "invisible")
ax2.xaxis.pane.set_edgecolor('w')
ax2.yaxis.pane.set_edgecolor('w')
ax2.zaxis.pane.set_edgecolor('w')
ax2.view_init(elev2,azim2) # elevation, azimuth
for i_g,gval in enumerate(g):
     if gval >= 0.1 and gval < 0.2:
     #if gval == 0.1: 
           scat = ax2.scatter(l, g[i_g], I0[:,i_g], c='g', marker='x', lw=2.0, zorder=4)
           scat = ax2.scatter(l, g[i_g], I1[:,i_g], c='m', marker='+', lw=2.0, zorder=4)
           scat = ax2.scatter(l, g[i_g], I2[:,i_g], c='darkorange', marker='o', s=6, lw=2.0, zorder=4)
alpha1=0.2
alpha2=0.2
alpha3=0.3
# Surface plots --------------------->
#
if flag_surf == 1:
  surf1  = ax2.plot_surface(Gam, Lam, I0, color='g', rstride=1, cstride=1, linewidth=0, antialiased=False, alpha=alpha1)
  surf2  = ax2.plot_surface(Gam, Lam, I1, color='m', rstride=1, cstride=1, linewidth=0, antialiased=False, alpha=alpha2)
  surf3  = ax2.plot_surface(Gam, Lam, I2, color='orange', rstride=1, cstride=1, linewidth=0, antialiased=False, alpha=alpha3)
else:
  surf1  = ax2.plot_wireframe(Gam, Lam, I0, color='g', alpha=alpha1)
  surf2  = ax2.plot_wireframe(Gam, Lam, I1, color='m',  alpha=alpha2)
  surf3  = ax2.plot_wireframe(Gam, Lam, I2, color='orange',  alpha=alpha3)
#
ax2.tick_params(axis='both',which='major',width=2,length=10,labelsize=16)
#ax2.tick_params(axis='both',which='minor',width=2,length=5) 
ax2.set_xlabel('$\lambda$', size=bsize)
ax2.xaxis.labelpad=10
ax2.set_ylabel('$\gamma$', size=bsize)
ax2.yaxis.labelpad=10
ax2.set_zlabel('Im $E$', size=bsize)
ax2.zaxis.labelpad=10
ax2.set_title('$t=${}, $U$={}, $\epsilon=${}'.format(t,U,eps), size=msize, y=0.95)
if flagU == 0: 
  ax2.text(1.6, 1.633, -3.1, '$E^-$', color='g', size=msize, zorder=1)
  ax2.text(0.93,1.6,3.4, '$E^+$', color='m', size=msize)
  ax2.text(-1.66,-0.47,-0.78, '$E^0$', color=mcol, size=msize, zorder=1) # for U=0
else:
  #ax2.text(2,-2.6,-3, 'Im $E$ = 0', color=mcol, size=msize, zorder=4)
  ax2.text(1.6, 1.633, -0.22, '$E^-$', color='g', size=msize, zorder=1)
  ax2.text(0.93,1.6,3.9, '$E^+$', color='m', size=msize)
  ax2.text(-1.6,1.65,-3.2, '$E^0$', color=mcol, size=msize, zorder=1)
bbox_props = dict(boxstyle="round", fc="silver", ec="0.5", alpha=0.8)
ax2.xaxis.set_major_locator(plt.MultipleLocator(1))
ax2.yaxis.set_major_locator(plt.MultipleLocator(1))
fig2.tight_layout()
ax2.set_rasterized(True)
fig2.savefig('3D_ImE_U{}.png'.format(U))
fig2.savefig('3D_ImE_U{}.eps'.format(U))
fig2.savefig('3D_ImE_U{}.pdf'.format(U))



plt.show()

    



