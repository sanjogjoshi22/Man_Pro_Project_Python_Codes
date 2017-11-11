import math
import matplotlib.pyplot as plt
import numpy as np

kr = 0.212
st = 0                    #entry angle for up-milling
ex = np.arange(0, 180, 1) #exit angle 


def axx(ex):              #function for alpha in x-direction
    return 0.5 * (
        (math.cos(2 * ex * math.pi / 180) - 2 * kr * ex * math.pi / 180 + kr * math.sin(2 * ex * math.pi / 180)) - (
            math.cos(2 * st * math.pi / 180) - 2 * kr * st * math.pi / 180 + kr * math.sin(2 * st * math.pi / 180)))


def axy(ex):              #function for alpha in cross axis 
    return 0.5 * (
        (-math.sin(2 * ex * math.pi / 180) - 2 * ex * math.pi / 180 + kr * math.cos(2 * ex * math.pi / 180)) - (
            -math.sin(2 * st * math.pi / 180) - 2 * st * math.pi / 180 + kr * math.cos(2 * st * math.pi / 180)))


def ayx(ex):              #function for alpha in cross axis
    return 0.5 * (
        (-math.sin(2 * ex * math.pi / 180) + 2 * ex * math.pi / 180 + kr * math.cos(2 * ex * math.pi / 180)) - (
            -math.sin(2 * st * math.pi / 180) + 2 * st * math.pi / 180 + kr * math.cos(2 * st * math.pi / 180)))


def ayy(ex):              #function for alpha in y-direction
    return 0.5 * (
        (-math.cos(2 * ex * math.pi / 180) - 2 * kr * ex * math.pi / 180 - kr * math.sin(2 * ex * math.pi / 180)) - (
            -math.cos(2 * st * math.pi / 180) - 2 * kr * st * math.pi / 180 - kr * math.sin(2 * st * math.pi / 180)))


R = axx(60) * ayy(60) - axy(60) * ayx(60)   #taking 60 degrees exit angle
wc = np.arange(0.5, 1000, 0.5)
wn1 = 262.16                                #natural frequency in x-direction
z1 = 0.048                                  #damping ration in x-direction
wn2 = 285.53                                #natural frequency in x-direction
z2 = 0.021                                  #damping ration in x-direction
ky = 226000
kt = 796
N = 4                                       #number of teeths
k = 0                           


def phi(wc1, wn, z):                        #transfer function
    r = wc1 / wn
    phi1 = (1 - r ** 2) / (ky * ((1 - r ** 2) ** 2 + (2 * z * r) ** 2))
    phi2 = (-2 * z * r) / (ky * ((1 - r ** 2) ** 2 + (2 * z * r) ** 2))
    return complex(phi1, phi2)


def a0(wc1, wn1, z1, wn2, z2):              #function for a0
    return phi(wc1, wn1, z1) * phi(wc1, wn2, z2) * R


def a1(wc1, wn1, z1, wn2, z2):              #function  for a1
    return axx(60) * phi(wc1, wn1, z1) - ayy(60) * phi(wc1, wn2, z2)


def lm1(wc1, wn1, z1, wn2, z2):             #function for lambda taling positive sign in root
    lm1 = a1(wc1, wn1, z1, wn2, z2) + (a1(wc1, wn1, z1, wn2, z2) ** 2 - 4 * a0(wc1, wn1, z1, wn2, z2)) ** 0.5
    return -lm1 / (2 * a0(wc1, wn1, z1, wn2, z2))
	
def lm2(wc1, wn1, z1, wn2, z2):             #function for lambda taling negetive sign in root
    lm1 = a1(wc1, wn1, z1, wn2, z2) - (a1(wc1, wn1, z1, wn2, z2) ** 2 - 4 * a0(wc1, wn1, z1, wn2, z2)) ** 0.5
    return -lm1 / (2 * a0(wc1, wn1, z1, wn2, z2))


def alim1(wc1, wn1, z1, wn2, z2, N, kt):    #critical width for lm1
    sk = lm1(wc1, wn1, z1, wn2, z2).imag / lm1(wc1, wn1, z1, wn2, z2).real
    lmr = lm1(wc1, wn1, z1, wn2, z2).real
    return (-2 * math.pi * lmr * (1 + sk ** 2)) / (N * kt)

def alim2(wc1, wn1, z1, wn2, z2, N, kt):    #critical width for lm2
    sk = lm2(wc1, wn1, z1, wn2, z2).imag / lm2(wc1, wn1, z1, wn2, z2).real  #real part/imagionary part
    lmr = lm2(wc1, wn1, z1, wn2, z2).real
    return (-2 * math.pi * lmr * (1 + sk ** 2)) / (N * kt)


def n1(wc1, wn1, z1, wn2, z2, N, k):        #spindle speed for lm1
    n1 = 60 * wc1
    sk = lm1(wc1, wn1, z1, wn2, z2).imag / lm1(wc1, wn1, z1, wn2, z2).real  #real part/imagionary part
    n2 = N * (((2 * k) + 1) * math.pi - 2 * math.atan(sk))
    return n1 / n2
def n2(wc1, wn1, z1, wn2, z2, N, k):        #spindle speed for lm2
    n1 = 60 * wc1
    sk = lm2(wc1, wn1, z1, wn2, z2).imag / lm2(wc1, wn1, z1, wn2, z2).real
    n2 = N * (((2 * k) + 1) * math.pi - 2 * math.atan(sk))
    return n1 / n2


alim_r = []

n_r_0 = []
n_r_1 = []
n_r_2 = []
n_r_3 = []

for i in wc:
    alim_r.append(alim2(i, wn1, z1, wn2, z2, N, kt))
    n_r_0.append(n2(i, wn1, z1, wn2, z2, N, 0))
    n_r_1.append(n2(i, wn1, z1, wn2, z2, N, 1))
    n_r_2.append(n2(i, wn1, z1, wn2, z2, N, 2))
    n_r_3.append(n2(i, wn1, z1, wn2, z2, N, 3))
  	
#plt.plot(n_r_0, alim_r, n_r_1, alim_r, n_r_2, alim_r, n_r_3, alim_r) for line plot
plt.scatter( n_r_0,alim_r, c='b', s=5, label = 'k = 0')
plt.scatter( n_r_1,alim_r, c='r',s=5,label = 'k = 1')
plt.scatter( n_r_2,alim_r, c ='g',s=5,label = 'k = 2')
plt.scatter( n_r_3,alim_r, c ='k',s=5,label = 'k = 3')
plt.legend(loc = 'best')
plt.xlim(0,5000)
plt.ylim(0,12000)
plt.xlabel('spindlespeed(rev/min)')
plt.ylabel('alim(mm)')
plt.show()

