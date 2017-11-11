import math
import numpy as np
import matplotlib.pyplot as plt

wn = 250            #natural frequency
z = 0.012		    #damping ratio
ky = 2.26*10**8     
Kf = 10**9

def depthofcut(wc):
	r = wc/wn
	y = (1 - r**2)
	f = ky*((1-r**2)**2 + (2*z*r)**2)
	realpart = y/f
	depth = -1/(2*Kf*realpart)
	return depth
def spindlespeed(wc, k):
	r = wc/wn
	num = 60*r*wn*2*math.pi
	den = 2*math.pi*k + 3*math.pi + 2*math.atan(2*z*r/(r**2 - 1))
	n = num/den
	return n
	
wc = np.arange(0., 1000., 0.1)
S_0 = []    #k = 0       
S_1 = []    #k = 1
S_2 = []	#k = 2
S_3 = []	#k = 3
S_4 = []	#k = 4
D = []

for i in wc:
	S_0.append(spindlespeed(i, 0))
	S_1.append(spindlespeed(i, 1))
	S_2.append(spindlespeed(i, 2))
	S_3.append(spindlespeed(i, 3))
	S_4.append(spindlespeed(i, 4))

	D.append(depthofcut(i))

plt.plot(S_0, D, label = 'k = 0')
plt.plot(S_1, D, label = 'k = 1')
plt.plot(S_3, D, label = 'k =2 ')
plt.plot(S_4, D, label = 'k = 3')
plt.xlabel('spindle speed(rev/min)')
plt.ylabel('alim(mm)')
plt.legend(loc='best')
plt.ylim(0,0.1)
plt.xlim(0,15000)
plt.show()