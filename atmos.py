# -*- coding: utf-8 -*-
"""
Created on Sun Feb 03 18:34:45 2019

@author: Miks
"""
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns

class Atmosphere:
    
    def __init__(self, vwind=[[5.], [0], [0.]], dens=1.225):
        '''Class to set up the atmospheric model for simulation, including wind velocity, 
        its variations, density, etc.'''
        self.velocity=np.matrix(vwind)
        self.vwind=np.matrix(vwind)
        self.dens=dens
        self.v0=6.
        self.h0=5.
        self.hr=0.6
        
    def getWind(self, h):
        self.velocity=self.v0*np.log(h/self.hr)/np.log(self.h0/self.hr)*(self.vwind/np.linalg.norm(self.vwind))
        return self.velocity
        
