# -*- coding: utf-8 -*-
"""
Created on Mon May 22 17:22:07 2023

@author: carlo
"""

import matplotlib.pyplot as plt
from numpy import*
import numpy as np
import cmath as cmath

###### DEFINO CTES #######

G = 6.67e-11
mT = 5.9736*10**24 #masa de la Tierra
mL = 0.07349*10**24 #masa de la Luna
d = 3.844*10**8 #distancia Tierra-LUna
omega = 2.6617*10**(-6)
rT = 6.37816*10**6 #radio Tierra
rL = 1.7374*10**6 #radio Luna
mu = mL / mT
delta = G*mT/d**3 

###### reescalamos ######

def conv (v):

    v[0] = v[0] / d
    v[2] = v[2] / (mc * d)  #mc= masa que le damos al cohete
    v[3] = v[3] / (mc* (d ** 2))
    
def deconversion (v, h):

    h[0] = v[0] * d
    h[1] = v[1]
    h[2] = v[2] * mc * d
    h[3] = v[3] * mc * (d ** 2)
    h[4] = v[4]

    
# definimos funciones a emplear #


def f1(v):#derivada de ~r
    return v[2]

def f2(v):#derivada de phi
    return v[3] / (v[0]**2)

def f3(v):#derivada pr
    rprima = np.sqrt(1.0 + v[0]**2 - 2.0 * v[0] * np.cos(v[1] - omega * v[4]))
    return (v[3]**2) / (v[0]**3) - delta * (1.0 / (v[0]**2) + (mu / (rprima**3)) * (v[0] - np.cos(v[1] - omega * v[4])))

def f4(v):#derivada pphi
    rprima = np.sqrt(1.0 + v[0]**2 - 2.0 * v[0] * np.cos(v[1] - omega * v[4]))
    return -((delta * mu * v[0]) / (rprima**3)) * np.sin(v[1] - omega * v[4])

def hamiltoniano (v): #v no puede estar reescalado

    A1 = (v[2] ** 2) / (2 * mc)
    A2 = (v[3] ** 2) / (2 * mc * (v[0] ** 2))
    A3 = -1.0 * G * mc* mT / v[0]
    A4 = -1.0 * G * mc * mL /(np.sqrt(d ** 2 + v[0] ** 2 - 2 * v[0] * d * np.cos(v[1] - omega * v[4])))
    A5 = -1.0 * omega * v[3]
    return A1 + A2 + A3 + A4 + A5



###### BUCLE ######

#a t=0
phi_i = 25* np.pi/180 #consideramos theta=0, lanzamiento tangencial, n grados para que se choquen
pphi_i= 0
mc=1
r_i=rT
pr_i = mc*11200 #momento lineal (velocidad de escape de la Tierra 11.2 km/h)
t=0


f_1=np.zeros(5)
f_2=np.zeros(5)
f_3=np.zeros(5)
f_4= np.zeros(5)


f = np.zeros(5)
f[0] = r_i
f[1] = phi_i
f[2] = pr_i
f[3] = pphi_i
f[4] = t

#inicializamos el vect de de H y el valor ded H
vH=np.zeros(5) 
H=0
N= 10000
h = 60


datoscohete = open('datoscohete.txt', 'w')
hamilton = open ('hamiltoniano.txt', 'w')
conv(f) #realizamos el cambio de escala
fig,ax=plt.subplots(1,1)

for i in range(1, N+1):
     
    datoscohete.write(str(0.0) + ', ' + str(0.0) + '\n') #el centro de la Tierra para nosotros es el (0,0)
    x = f[0] * np.cos(f[1]) #posición en el eje x del cohete
    y = f[0] * np.sin(f[1]) #posición en el eje y del cohete
    datoscohete.write(str(x) + ', ' + str(y) + '\n')
    datoscohete.write(str(np.cos(omega * f[4])) + ', ' + str(np.sin(omega * f[4])) + '\n' + '\n')#posición de la luna dividida entre d (pq reescalamos)
    deconversion(f, vH)
    H= hamiltoniano(vH)
    hamilton.write(str(H) + ',')
    l=ax.plot(t,H, "*",color="mediumvioletred",ms=1, label="H respecto a T")
        
    k1 = np.zeros(5)
    k2 = np.zeros(5)
    k3 = np.zeros(5)
    k4 = np.zeros(5)
        
    k1[0] = h * f1(f)
    k1[1] = h * f2(f)
    k1[2] = h * f3(f)
    k1[3] = h * f4(f)
        
    for j in range(4):
        f_2[j] = f[j] + k1[j] / 2.0 
    f_2[4] = f[4] + h / 2.0
            
    k2[0] = h * f1(f_2)
    k2[1] = h * f2(f_2)
    k2[2] = h * f3(f_2)
    k2[3] = h * f4(f_2)
        
    for j in range(4):
        f_3[j] = f[j] + k2[j] / 2.0 
    f_3[4] = f[4] + h / 2.0

    k3[0] = h * f1(f_3)
    k3[1] = h * f2(f_3)
    k3[2] = h * f3(f_3)
    k3[3] = h * f4(f_3)
        
    for j in range(4):
        f_4[j] = f[j] + k3[j]
    f_4[4] = f[4] + h

    k4[0] = h * f1(f_4)
    k4[1] = h * f2(f_4)
    k4[2] = h * f3(f_4)
    k4[3] = h * f4(f_4)

    for j in range(4):
        f[j] = f[j] + (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6.0
    t=t+h
    f[4]=t 
            
datoscohete.close()
hamilton.close()
plt.legend()
plt.show()