# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 15:25:47 2022

@author: Mahathi
"""
import math
import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt

#constants
G = 6.67430*(10**(-11)) #Gravitational constant
Ms = 1.989*(10**30) #Mass of sun
Me = 5.97219*(10**24) #Mass of earth
day = 60*60*24 #in seconds
year = 60*60*24*365 #in seconds

#inital conditions
r0 = 149587871*1000 #meters #from distance from sun
v0 = 20*1000 #meters/sec #velocity at perihelion
ti = 4.7*year #time to start plotting in seconds
tf = 5*year #time to start plotting in seconds

#single calculations
k = G*Ms*Me
L = Me*v0*r0
c = (L**2)/ (k*Me)
E = 0.5*Me*(v0**2) - k/r0 
e = np.sqrt((E*2*c/k) + 1)
a = c/(1 - e**2)
T = np.sqrt(4*(np.pi**2)*(a**3)/(G*Ms))
w = 2*np.pi/T
psi_0 = np.arccos((1 - (r0/a))/e)
phi_0 = np.arccos(((c/r0) - 1)/e)
tp = (psi_0 - e*np.sin(psi_0))/w

#defining t
tstep = 1000
t = np.linspace(ti, tf, tstep)

#defining zeta
z = np.zeros(len(t))
for i in range(len(t)):
    z[i] = w*(t[i] + tp)

#finding psi

psi = np.zeros(len(t))

def f(zeta):
    return zeta - e*np.sin(zeta) - z[i]

def fd(zeta):
    return 1 - e*np.cos(zeta)

def newtonraphson(zeta, i):
    
        h = f(zeta)/fd(zeta)
        while abs(h) >= 0.0001:
            zeta = zeta - h
            h = f(zeta)/fd(zeta)

        psi[i] = zeta
    
for i in range(len(t)):    
    newtonraphson(psi_0, i)

#finding r from psi
r = np.zeros(len(t))
for i in range(len(t)):
    r[i] = a*(1-e*np.cos(psi[i]))

#finding phi from r 
phi = np.zeros(len(t))
for i in range(len(t)):
    if i!=0 and r[i] > r[i - 1]:
        phi[i] =  + np.arccos((c/r[i] - 1)/e) + phi_0
    else:
        phi[i] = - np.arccos((c/r[i] - 1)/e) + phi_0

#plotting
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
fig.set_facecolor("black")
ax.set_facecolor("black")
ax.set_rlabel_position(280)
ax.tick_params(axis="both",colors="White", labelsize = 10)
for spine in ax.spines.values():
    spine.set_edgecolor("White")
ax.plot(phi,r, color="turquoise", linewidth = 3)
fig.savefig("planetary_orbit.png")

fig, ax1= plt.subplots()
fig.set_facecolor("black")
ax1.set_facecolor("black")
ax1.set_ylabel("$\zeta$",color="white", size = 13)
ax1.set_xlabel("$\psi$",color="white", size = 13)
ax1.tick_params(axis="both",colors="White", labelsize = 10)
for spine in ax.spines.values():
    spine.set_edgecolor("White")
ax1.grid(color = "lightgrey", linestyle = "--")
ax1.plot(z, psi, color="turquoise", linewidth = 2)
fig.savefig("planetary_orbit2.png")


