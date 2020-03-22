# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 09:47:09 2018

@author: Miks
"""

import math
import numpy as np
import scipy as sc
from visual import *

def mag(x):
    return math.sqrt(i**2 for i in x);
    
    
def matmul(a, b):
    out=np.empty((np.shape(a)[0], np.shape(b)[1]))
    a=np.matrix(a)
    b=np.matrix(b)
    for i in range(np.shape(b)[1]):
        for j in range(np.shape(a)[0]):
            out[j, i]=np.dot(a[j, :], b[:, i])
    return out
    

def part_positions(n, x, y, z=0):
    '''INITIALIZE PARTICLE POSITIONS, OUTPUTS: POSITION MATRIX FOR CALCULCATIONS, POSITION LIST FOR VISUALIZATION'''
    pos=np.array(np.ones((3,n)))
    pos[0]=np.linspace(0,x*(1/5.),n)     
    pos[1]=np.linspace(0, (0.5/5.)*y, n)   
    pos[2]=np.linspace(0, z, n)
    pos=np.matrix(pos)
    posd=[]
    for i in range(n):
        posd.append(list(np.array(np.transpose(pos[:, i])))[0])
            
def vectorize(matrixin):
    tmp=[]
    for i in range(np.shape(matrixin)[1]):
        tmp.append(vector(matrixin[0, i], matrixin[1, i], matrixin[2, i]))
    return tmp

def Tx(phi):
    return np.matrix([[1, 0, 0],
                      [0, np.cos(phi), np.sin(phi)],
                      [0, -np.sin(phi), np.cos(phi)]])
                      
def Ty(phi):
    return np.matrix([[np.cos(phi), 0, -np.sin(phi)],
                      [0, 1, 0],
                      [np.sin(phi), 0, np.cos(phi)]])

def Tz(phi):
    return np.matrix([[np.cos(phi), np.sin(phi), 0],
                      [-np.sin(phi), np.cos(phi), 0],
                      [0, 0, 1]])
    