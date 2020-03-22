# -*- coding: utf-8 -*-
"""
Created on Sun Sep 30 18:20:39 2018

@author: Miks
"""

from visual import *
import numpy as np
from time import time
from libraries import *



class disp:
    
    def __init__(self, scrw=1000, scrh=500, framerate=60.):
        self.scrw=scrw
        self.scrh=scrh
        self.framerate=framerate
        self.label='Simulation window'
        self.bgcolor=[255, 255, 255]
    def create(self, showaxis=True):
        blk=[0, 0, 0]
        self.window=display(title=self.label, background=self.bgcolor, width=self.scrw, height=self.scrh)
        if showaxis==True:
            xax = arrow(pos=(0, 0, 0), axis=(self.scrh*(2/5.), 0,0), shaftwidth=0.5, color=blk)
            yax = arrow(pos=(0, 0, 0), axis=(0,self.scrh*(2/5.),0), shaftwidth=0.5, color=blk)
            zax = arrow(pos=(0, 0, 0), axis=(0, 0, self.scrh*(2/5.)), shaftwidth=0.5, color=blk)
        print('-->Display instance created')
        
class ExParticleSystem:
    
   def __init__(self, nparticles=10, mass=0.25, x_last=200., y_last=40., z_last=0., g_acc=9.81, springc=9., dampc=15., framerate=60., vwind=10.):
        '''ParticleSystem assumes linear intial particle spacing between the origin and the 
          predefined position of the last particle.
            
          Intializing function, all predefined mass-spring-damper constants defined within here:
           - nparticles - amount of mass particles in the tether
           - masss - mass of a single particle
           - x_, y_, z_last - coordinates of the last particle of the tether
           - t_radius - radius of the tether (for vizualization purposes)
           - g - gravitational acceleration vector
           - pos - positional matrix of the particle system
           - Jx - positional Jacobian of the system
           - Jv - velocity Jacobian of the system
           - posd - intermediate positional matrix for visualization input
           - vel - velocity matrix of the system
           - fmag - force matrix of the system
           - acc - acceleration matrix of the system
           - framerate - framerate of the visualization
           - Rl, R - array containing initial lengths of the springs
           - incl_final - boolean, whether the position of last particle should be updated
           - vwind - wind velocity
         '''
        self.nparticles=nparticles
        self.x_last=x_last
        self.y_last=y_last
        self.z_last=z_last
        self.mass=mass
        self.t_radius=1.0
        self.g=np.matrix([[0], [g_acc], [0]])
        self.springc=springc
        self.dampc=dampc
        self.pos=np.array(np.ones((3,self.nparticles)))
        self.pos[0]=np.linspace(0,self.x_last*(4/5.),self.nparticles)     
        self.pos[1]=np.linspace(0, (2/5.)*self.y_last, self.nparticles)   
        self.pos[2]=np.linspace(0, self.z_last, self.nparticles)
        self.pos=np.matrix(self.pos)
        self.posd=[]
        for i in range(self.nparticles):
            self.posd.append(list(np.array(np.transpose(self.pos[:, i])))[0])
        self.teth=curve(pos=self.posd, radius=self.t_radius, color=(211, 211, 211))  
        self.vel=np.matrix(np.zeros((3,self.nparticles)))                                 
        self.fmag=np.matrix(np.zeros((3,self.nparticles)))                                
        self.acc=np.matrix(np.zeros((3,self.nparticles)))
        self.Rl=[]
        self.framerate=framerate
        self.dt=1/self.framerate
        for i in range(self.nparticles-1):
            self.Rl.append(np.abs((self.pos[:, i]-self.pos[:, i+1])))
        self.R=np.array(self.Rl)
        self.incl_final=False
        self.vwind=vwind
        print('-->Particle system initialized')
        
   def get_diff(self, i):
        '''Obtains the difference vectors for ith particle w.r.t. the adjacent particles'''
        if i==self.nparticles-1:
            self.dxba=self.pos[:, i]-self.pos[:, i-1]
            self.dvba=self.vel[:, i]-self.vel[:, i-1]
        else:
            self.dxab=self.pos[:, i]-self.pos[:, i+1]
            self.dxba=self.pos[:, i]-self.pos[:, i-1]
            self.dvab=self.vel[:, i]-self.vel[:, i+1]
            self.dvba=self.vel[:, i]-self.vel[:, i-1]
   
   def get_norm(self, i):
        '''Obtains the normal vectors between particles'''
        if i==self.nparticles-1:
            self.nba=self.dxba/np.linalg.norm(self.dxba)
        else:
            self.nab=self.dxab/np.linalg.norm(self.dxab)
            self.nba=self.dxba/np.linalg.norm(self.dxba)
   
   def get_forces(self, i):
        '''Obtains the implicit spring & damping forces'''
        if i==self.nparticles-1:
            self.Fsba=-self.springc*(np.linalg.norm(self.dxba)-np.linalg.norm(self.R[i-1]))*((np.linalg.norm(self.dxab)-np.linalg.norm(self.R[i-1]))>0)*self.nba
            self.Fdba=np.asscalar(-self.dampc*(np.dot(np.transpose(self.dvba), self.dxba)/np.linalg.norm(self.dxba)))*self.nba
        else:
            self.Fsab=-self.springc*(np.linalg.norm(self.dxab)-np.linalg.norm(self.R[i]))*((np.linalg.norm(self.dxab)-np.linalg.norm(self.R[i]))>0)*self.nab
            self.Fsba=-self.springc*(np.linalg.norm(self.dxba)-np.linalg.norm(self.R[i-1]))*((np.linalg.norm(self.dxab)-np.linalg.norm(self.R[i-1]))>0)*self.nba
            self.Fdab=np.asscalar(-self.dampc*(np.dot(np.transpose(self.dvab), self.dxab)/np.linalg.norm(self.dxab)))*self.nab
            self.Fdba=np.asscalar(-self.dampc*(np.dot(np.transpose(self.dvba), self.dxba)/np.linalg.norm(self.dxba)))*self.nba
   
   def updatePS(self):
        '''Update after each dt to obtain the new positions '''
        if self.incl_final==True:
            final=self.nparticles
        else:
            final=self.nparticles-1
        for i in range(1, final):
            self.get_diff(i)
            self.get_norm(i)
            self.get_forces(i)
            if i==self.nparticles-1:
                self.fmag[:, i]=-self.mass*self.g+self.Fsba+self.Fdba#+np.matrix([[0], [9.1], [10*sin(self.tr)]])
            else:
                self.fmag[:, i]=-self.mass*self.g+self.Fsab+self.Fsba+self.Fdab+self.Fdba
            self.acc[:, i]=self.fmag[:, i]/self.mass
            self.vel[:, i]+=self.dt*self.acc[:, i]
            self.pos[:, i]+=self.dt*self.vel[:, i]
            self.posd[i]=list(np.array(np.transpose(self.pos[:, i])))[0]
            self.teth.pos=self.posd
   def start_sim(self):
        self.tr=0
        running=True
        while running:
            self.t=time()
            self.tr+=self.dt
            self.updatePS()
            rate(self.framerate)
            

'''To-Dos: Finish implementing the Jacobian '''
screen=disp()
screen.create()
ps=ExParticleSystem()
ps.start_sim()