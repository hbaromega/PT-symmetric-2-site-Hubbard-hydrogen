from matplotlib.widgets import Slider, Button, RadioButtons 
import numpy as np
import matplotlib.pyplot as plt
import cmath


t = 1.0
eps = 0.5
E = [0,0,0]
#E = np.zeros(3)
#E = np.zeros(3, dtype='complex_')
#E = (0,0,0)
gam_min = -1.0    
gam_max = 1.0   
gam_init = 0.1 
gam = gam_init
gamma_c = 0.19

U_min = -4    
U_max = 4   
U_init = 2 
U = U_init
tol = 0.0001 # tolerance factor

lam_list=np.linspace(-1.5, 1.5, 6001)

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.35)

# Slider layout
slider_Uax = plt.axes([0.25, 0.2, 0.65, 0.02])
slider_gamax = plt.axes([0.25, 0.1, 0.65, 0.02])
              # Syntax: plt.axes((left, bottom, width, height), facecolor='w')


# Parameters in the cubic eq.
S=2.0*eps
#Emid=2.0*(eps-gamma)+U/2.0
flag_imag = 0
count = 0

def swap(var1,var2):
          tmp = var1
          var1 = var2 
          var2 = tmp
          return var1, var2
      


def cubicroot(U,K,L):
         coeff=[1,-U,-K,L]
         return np.roots(coeff) 

def sol1(U,lam,gam):
      tm=t-lam; tp=t+lam
      M=4*tp*tm
      igamma=complex(0,gam)
      epsp=eps+igamma
      epsm=eps-igamma
      Dsq=(epsp-epsm)**2
      K=Dsq+M
      L=Dsq*U
      roots=cubicroot(U,K,L)
      roots=S+U-roots
      ImE1 = roots[0].imag
      ImE2 = roots[1].imag
      ImE3 = roots[2].imag
      #'''
      global flag_imag # Need to declare this since it's in a condn below

      for x in roots:
        #if x.imag != 0:
        if abs(x.imag) >= tol:
             flag_imag = 1 # Complex roots exist

      print('lambda =',lam, 'flag_imag =', flag_imag)


      #'''
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
      '''
      if abs(lam) >1.5 and roots[0].imag > 0.0: # Correcting anomaly 
          E[0]=roots[1] # g 
          E[1]=roots[0] # m

          flag_imag = 1
          #count
       '''
          

      ''' 
      if E[0].real > E[1].real: # If blue line comes over red line
          # Swap
        ctmp  = E[0] 
        E[0] = E[1]
        E[1] = ctmp
      '''




      ImE1 = E[0].imag
      ImE2 = E[1].imag
      ImE3 = E[2].imag

      '''
      if gam == gam_init:
          E[0]=roots[0]
          E[1]=roots[1]
          E[2]=roots[2]
          ImE1 = E[1].imag
          ImE2 = E[2].imag
          ImE3 = E[0].imag
      '''

      # Color convention
      # Im E1 --> green (E^-)
      # Im E2 ---> magenta (E^+) [always positive]
      # Im E3 ----> orange (E^0)

      # Make sure ImE2 always positive definite
      if ImE1 > 0 and ImE2 < 0 and abs(ImE3) < tol:
          print('lam,ImE1,ImE2=',lam,ImE1,ImE2)
          # Swap ImE1 & ImE2
          ImE1, ImE2 = swap(ImE1,ImE2)
      
      if ImE3 > 0 and ImE2 < 0 and abs(ImE1) < tol:
          print('lam,ImE3,ImE2=',lam,ImE3,ImE2)
          # Swap ImE3 & ImE2
          ImE3, ImE2 = swap(ImE3,ImE2)

      if abs(gam) >= gamma_c: 
          # locate conjugate pairs # swap ImE1 & ImE3
          if ImE2 < tol and ImE3 > tol:
             # Swap ImE3 & ImE2
             ImE3, ImE2 = swap(ImE3,ImE2)
          if abs(ImE3) < tol and ImE1 < tol: # ImE3 should be always -ve
             # Swap ImE1 & ImE3
             ImE1, ImE3 = swap(ImE1,ImE3)


      # Now between ImE1 and ImE3 
 
      #ReE1, ReE2, ReE3 = E.real
      #(E[0].real, E[1].real, E[2].real)
      #ImE1, ImE2, ImE3 = E.imag

      #'''
      return ImE1, ImE2, ImE3 
      #return E[0].imag, E[1].imag, E[2].imag

#Function to update solutions after modifying lambda
def update_sols(U,gam):
  res1=[]
  res2=[]
  res3=[]
  for lam in lam_list:
    a1, a2, a3 = sol1(U,lam,gam) 
    res1.append(a1)
    res2.append(a2)
    res3.append(a3)
  return res1,res2,res3


#Initialising plot with solutions for a_init
sols1,sols2,sols3 = update_sols(U_init,gam_init)
plot1, = ax.plot(lam_list, sols1, 'go', ms=4, label='$E^-$')
plot2, = ax.plot(lam_list, sols2, 'm.', ms=6, label='$E^+$')
plot3, = ax.plot(lam_list, sols3, 'orange', marker='.', ls='None', ms=1, label='$E^0$')
ax.set_xlabel('$\lambda$')
ax.set_ylabel('Im $E$')
ax.legend(loc='best')  
ax.set_title('Roots of cubic characteristic eq.: \n$X^3 − UX^2 − KX − L = 0$')
 
# Create a slider 
U_slider = Slider(ax=slider_Uax,label='$U$',valmin=U_min,valmax=U_max,valinit=U_init)
gam_slider = Slider(ax=slider_gamax,label='$\gamma$',valmin=gam_min,valmax=gam_max,valinit=gam_init)

# Update function
def update(val):
    #updating y data
    sols1,sols2,sols3 = update_sols(U_slider.val,gam_slider.val)
    plot1.set_ydata(sols1) 
    plot2.set_ydata(sols2) 
    plot3.set_ydata(sols3) 
   
    #updating y limits
    sols1_abs=np.abs(sols1)
    sols2_abs=np.abs(sols2)
    sols3_abs=np.abs(sols3)
    max_ylim=np.amax(np.concatenate((sols1_abs,sols2_abs,sols3_abs)))
    ax.set_ylim([-1.1*max_ylim,1.1*max_ylim])
    fig.canvas.draw_idle() # redraw the plot

# Execute when parameter gets updated 
U_slider.on_changed(update)
gam_slider.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
axcolor = 'lightgoldenrodyellow'
resetax = plt.axes([0.7, 0.02, 0.1, 0.04]) #  plt.axes(left, bottom, width, height)
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

def reset(event):
    U_slider.reset()
    gam_slider.reset()
button.on_clicked(reset)



plt.show()
