import math
import numpy as np
import matplotlib.pyplot as plt

wn = 250          #natural frequency in axis-1 direction at 30 degrees angle from main axis
Wn = 150          #natural frequency in axis-2 direction at -45 degrees angle from main axis
z = 0.012        #damping in axis-1 direction
Z = 0.01		 #damping in axis-2 direction
ky = 22600000
Ky = 2130000
#Real part of transfer function varying with chatter frequency wc
def realpart(wc): 
	r = wc/wn
	R = wc/Wn
	y = (1 - r**2)*(math.cos(30/math.pi))**2   # y/f real-part in axis-1 direction
	f = ky*((1-r**2)**2 + (2*z*r)**2)          # Y/F real-part in axis-2 direction
	Y = (1 - R**2)*(math.cos(45/math.pi))**2
	F = Ky*((1-R**2)**2 + (2*Z*R)**2)
	return (y/f + Y/F)

wc = np.arange(1., 300., 0.01)
plt.plot(wc, realpart(wc))
plt.xlabel('frequency(rad/sec)')
plt.ylabel('real part')
plt.show()