# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 18:38:37 2019

@author: Miks
"""

import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns

class KiteModel:
    
    def __init__(self, model='point mass', mass=3.0, surface=6.):
        self.model=model
        self.surface=surface
        self.mass=mass
        self.va=None
        self.vaxz=None
        self.vaxy=None
        self.zk=None
        self.yk=None
        self.xk=None
        self.c=None
        self.AsAratio=0.5 #Apparently arbitrary, ok
        self.CS=0.6 # Found from KiteSim code
        self.alphatop=0
        self.betatop=0
        self.alphaside=0
        self.betaside=0
        self.CL=1.0 #Default Lift coefficient
        self.CD=0.2 #Default Drag coefficient
        self.alphapower=20.
        self.alphadepower=-20.
        self.alpharefL=[-180, -160, -90, -20, 0, 20, 40, 90, 160, 180]
        self.clref=[0, 0.5, 0, 0.1, 0.1, 1, 1, 0, -0.5, 0]
        self.alpharefD=[-180, -160, -90, -20, 0, 20, 90, 160, 180]
        self.cdref=[0.5, 0.5, 1.0, 0.2, 0.1, 0.2, 1, 0.5, 0.5]
        
    def getCL(self, alpha):
        self.CL=np.interp(alpha, self.alpharefL, self.clref)
    
    def getCD(self, alpha):
        self.CD=np.interp(alpha, self.alpharefD, self.cdref)
        
    def getForces(self, rho, steering, power):
        self.alphatop+=np.radians(self.alphadepower+power*(self.alphapower-self.alphadepower))
        self.getCL(self.alphatop)
        self.getCD(self.alphatop)
        self.L=0.5*rho*np.linalg.norm(self.va)*np.linalg.norm(self.va)*self.surface*self.CL*(np.cross(self.yk, self.va, axis=0)/np.linalg.norm(np.cross(self.yk, self.va, axis=0)))
        self.S=0.5*rho*np.linalg.norm(self.va)*np.linalg.norm(self.va)*self.surface*self.AsAratio*steering*self.CS*self.yk
        self.D=(0.5*rho*np.linalg.norm(self.va)*np.linalg.norm(self.va)*self.surface*self.CD*(self.va/np.linalg.norm(self.va)))*(1+0.6*steering)
        
#sns.set(context='poster', font_scale=1.1, style='whitegrid')
#sns.set_palette(sns.color_palette("hls", 4))
#kt=KiteModel()
#plt.figure()
#plt.plot(kt.alpharefL, kt.clref, '--', linewidth=3, marker='d', markersize=10, label=r'$C_{L}$')
#plt.plot(kt.alpharefD, kt.cdref, '--', linewidth=3, marker='d', markersize=10, label=r'$C_{D}$')
#plt.xlabel(r'$\alpha$, [deg]')
#plt.ylabel(r'$C_{L}, C_{D}$, [-]')
#plt.legend()
#plt.show()
#
#        
#        