# -*- coding: utf-8 -*-
"""
Created on Sun Sep 30 18:20:39 2018

@author: Miks
"""
from visual import *
import numpy as np
import time
from libraries import *
from scipy.sparse.linalg import bicgstab
import matplotlib.pyplot as plt
from win32api import GetSystemMetrics
import seaborn as sns

class disp:
    
    def __init__(self, scrw=1000, scrh=500, framerate=60.):
        self.scrw=scrw
        self.scrh=scrh
        self.framerate=framerate
        self.label='Simulation window'
        self.bgcolor=[255, 255, 255]
        
    def create(self, showaxis=True):
        blk=[0, 0, 0]
        self.window=display(title=self.label, background=vector(self.bgcolor), width=GetSystemMetrics(0)-200, height=GetSystemMetrics(1), fullscreen=True)
        if showaxis==True:
            xax = arrow(pos=(0, 0, 0), axis=(self.scrh*(1/5.), 0,0), shaftwidth=0.5, color=blk)
            yax = arrow(pos=(0, 0, 0), axis=(0,self.scrh*(1/5.),0), shaftwidth=0.5, color=blk)
            zax = arrow(pos=(0, 0, 0), axis=(0, 0, self.scrh*(1/5.)), shaftwidth=0.5, color=blk)
        print('-->Display instance created')
            
class ImParticleSystem:
    
    def __init__(self, nparticles=26, x_last=300, y_last=0, z_last=0., 
                 g_acc=9.81, k0=614620, c0=40.,  framerate=60.,breeze=False, screen=None, kite=None):
        '''ParticleSystem assumes linear intial particle spacing between the origin and the 
          predefined position of the last particle.
            
          Intializing function, all predefined mass-spring-damper constants defined within here:
           - nparticles - amount of mass particles in the tether
           - masss - mass of a single particle
           - x_, y_, z_last - coordinates of the last particle of the tether
           - t_radius - radius of the tether (for vizualization purposes)
           - g - gravitational acceleration vector.
           - fmag - force matrix of the system
           - acc - acceleration matrix of the system
           - framerate - framerate of the visualization
           - Rl, R - array containing initial lengths of the springs
           - incl_final - boolean, whether the position of last particle should be updated
           - vwind - wind velocity
         '''
        self.nparticles=int(nparticles)
        self.x_last=x_last
        self.y_last=y_last
        self.z_last=z_last
        self.last=np.matrix([[x_last], [y_last], [z_last]])      
        self.rho=1.225
        self.rhot=960
        self.kite=kite
        self.t_radius=2e-03
        self.CDt=1.22
        self.screen=screen
        self.kpm=0.012
        self.mass=(self.kpm*np.linalg.norm(self.last))/(self.nparticles-1)
        self.g=np.matrix([[0], [-g_acc], [0]])
        self.young=4.891e10
        self.k0=k0
        self.c0=c0
        self.springc=(self.nparticles-1)*self.k0/np.linalg.norm(self.last)
        self.dampc=(self.nparticles-1)*self.c0/np.linalg.norm(self.last)
        self.dim=3*self.nparticles
        self.M=np.matrix(np.identity(self.dim)*self.mass)
        self.M[0:3, 0:3]=np.identity(3)*self.mass*0.5
        self.M[3*self.nparticles-3:3*self.nparticles, 3*self.nparticles-3:3*self.nparticles]=np.identity(3)*self.mass*0.5
        self.pos=np.array(np.zeros((3,self.nparticles)))
        self.pos[0]=np.linspace(0,self.x_last,self.nparticles)     
        self.pos[1]=np.linspace(0, self.y_last, self.nparticles)   
        self.pos[2]=np.linspace(0, self.z_last, self.nparticles)
        self.pos=np.matrix(self.pos)
        self.R=np.linalg.norm(np.abs((self.pos[:, 2]-self.pos[:, 1])))
        self.R1=self.R
        self.Ap=2*self.t_radius*(np.linalg.norm(self.last)/(self.nparticles-1))        
        self.Jx=np.matrix(np.zeros((self.nparticles*3, self.nparticles*3)))
        self.Jv=self.dampc*np.matrix(np.identity(self.nparticles*3))
        self.teth=curve(pos=vectorize(self.pos), display=self.screen, radius=0.4, color=vector(211, 211, 211))
        self.vel=np.matrix(np.zeros((3,self.nparticles)))                                 
        self.fmag=np.matrix(np.zeros((3,self.nparticles)))
        self.FDt=np.matrix(np.zeros((3,self.nparticles)))                        
        self.acc=np.matrix(np.zeros((3,self.nparticles)))
        self.fspring=np.matrix(np.zeros((3,self.nparticles)))
        self.framerate=framerate
        self.dt=1/self.framerate
        self.incl_final=False
        self.vwind=np.matrix([[0.], [0.], [0]])
        self.Fg=self.g*self.mass
        self.Xfactor=0.05
        self.orig=np.matrix([[0],
                             [0],
                             [0]])
        self.reelincr=0.001*self.R
        if np.linalg.norm(self.vwind)==0:
            self.breeze=False
        else:
            self.breeze=breeze
        print('-->Particle system initialized')
        self.reelout=False
        
    def getDiff(self, i):
        '''Obtains the difference vectors for ith particle w.r.t. the adjacent particles'''
        if i==0:
            self.dxab=self.pos[:, i]-self.pos[:, i+1]
            self.dvab=self.vel[:, i]-self.vel[:, i+1]
        elif i==self.nparticles-1:
            self.dxba=self.pos[:, i]-self.pos[:, i-1]
            self.dvba=self.vel[:, i]-self.vel[:, i-1]
        else:
            self.dxab=self.pos[:, i]-self.pos[:, i+1]
            self.dxba=self.pos[:, i]-self.pos[:, i-1]
            self.dvab=self.vel[:, i]-self.vel[:, i+1]
            self.dvba=self.vel[:, i]-self.vel[:, i-1]
               
    def getNorm(self, i):
        # type: (object) -> object
        '''Obtains the normal vectors between particles'''
        if i==0:
            self.nab=self.dxab/np.linalg.norm(self.dxab)
        elif i==self.nparticles-1:
            self.nba=self.dxba/np.linalg.norm(self.dxba)
        else:
            self.nab=self.dxab/np.linalg.norm(self.dxab)
            self.nba=self.dxba/np.linalg.norm(self.dxba)
   
    def getExForces(self, i, frame):
        '''Adds (revalues) the force magnitude values for a particle in the beginning of each loop,
        by addding the gravitational force + drag force'''
        self.vapp=self.vwind-self.vel[:, i]
        if frame>1:
            self.FDt[:, i]=0.5*self.rho*(np.linalg.norm(self.vapp)**2)*self.Ap*self.CDt*(self.vapp/np.linalg.norm(self.vapp))
        else:
            self.FDt[:, i]=0
        self.fmag[:, i]+=self.Fg+self.FDt[:, i]
        
    def getImForces(self, i):
        '''Obtains the implicit spring & damping forces'''
        if i==0:
            self.Fsab=-self.springc*(np.linalg.norm(self.dxab)-np.linalg.norm(self.R))*((np.linalg.norm(self.dxab)-np.linalg.norm(self.R))>0)*self.nab
            self.Fdab=np.asscalar(-self.dampc*(np.dot(np.transpose(self.dvab), self.dxab)/np.linalg.norm(self.dxab)))*self.nab
            self.fmag[:, i]+=self.Fsab+self.Fdab
            self.fspring[:, i]=self.Fsab+self.Fdab
        elif i==self.nparticles-1:
            self.Fsba=-self.springc*(np.linalg.norm(self.dxba)-np.linalg.norm(self.R))*((np.linalg.norm(self.dxab)-np.linalg.norm(self.R))>0)*self.nba
            self.Fdba=np.asscalar(-self.dampc*(np.dot(np.transpose(self.dvba), self.dxba)/np.linalg.norm(self.dxba)))*self.nba
            self.fmag[:, i]+=self.Fsba+self.Fdba
            self.fspring[:, i]=self.Fsba+self.Fdba
        elif i==1:
            self.Fsab=-self.springc*(np.linalg.norm(self.dxab)-np.linalg.norm(self.R1))*((np.linalg.norm(self.dxab)-np.linalg.norm(self.R1))>0)*self.nab
            self.Fsba=-self.springc*(np.linalg.norm(self.dxba)-np.linalg.norm(self.R1))*((np.linalg.norm(self.dxab)-np.linalg.norm(self.R1))>0)*self.nba
            self.Fdab=np.asscalar(-self.dampc*(np.dot(np.transpose(self.dvab), self.dxab)/np.linalg.norm(self.dxab)))*self.nab
            self.Fdba=np.asscalar(-self.dampc*(np.dot(np.transpose(self.dvba), self.dxba)/np.linalg.norm(self.dxba)))*self.nba
            self.fmag[:, i]+=self.Fsab+self.Fsba+self.Fdab+self.Fdba
            self.fspring[:, i]=self.Fsba+self.Fdba
        else:
            self.Fsab=-self.springc*(np.linalg.norm(self.dxab)-np.linalg.norm(self.R))*((np.linalg.norm(self.dxab)-np.linalg.norm(self.R))>0)*self.nab
            self.Fsba=-self.springc*(np.linalg.norm(self.dxba)-np.linalg.norm(self.R))*((np.linalg.norm(self.dxab)-np.linalg.norm(self.R))>0)*self.nba
            self.Fdab=np.asscalar(-self.dampc*(np.dot(np.transpose(self.dvab), self.dxab)/np.linalg.norm(self.dxab)))*self.nab
            self.Fdba=np.asscalar(-self.dampc*(np.dot(np.transpose(self.dvba), self.dxba)/np.linalg.norm(self.dxba)))*self.nba
            self.fmag[:, i]+=self.Fsab+self.Fsba+self.Fdab+self.Fdba
            self.fspring[:, i]=self.Fsba+self.Fdba

    def getJx(self, i):
        '''Obtains and fills in the positional Jacobian Jx for the ith particle, fills in the global one'''
        if i==0:
            self.Jxab=self.springc*(-np.identity(3)+((self.R)/(np.linalg.norm(self.dxab))) \
            *(np.identity(3)-(1/(np.linalg.norm(self.dxab)**2))*self.dxab*np.transpose(self.dxab)))
            self.Jx[3*i:3*i+3, 3*i:3*i+3]=self.Jxab
            self.Jx[3*i:i+3, 3*(i+1):3*(i+1)+3]=-self.Jxab    
        elif i==self.nparticles-1:
            self.Jxba=self.springc*(-np.identity(3)+((self.R)/(np.linalg.norm(self.dxba))) \
            *(np.identity(3)-(1/(np.linalg.norm(self.dxba)**2))*self.dxba*np.transpose(self.dxba)))
            self.Jx[3*i:3*i+3, 3*i:3*i+3]=self.Jxba
            self.Jx[3*i:3*i+3, 3*(i-1):3*(i-1)+3]=-self.Jxba
        else:
            self.Jxab=self.springc*(-np.identity(3)+((self.R)/(np.linalg.norm(self.dxab))) \
            *(np.identity(3)-(1/(np.linalg.norm(self.dxab)**2))*self.dxab*np.transpose(self.dxab)))
            self.Jxba=self.springc*(-np.identity(3)+((self.R)/(np.linalg.norm(self.dxba))) \
            *(np.identity(3)-(1/(np.linalg.norm(self.dxba)**2))*self.dxba*np.transpose(self.dxba)))
            self.Jx[3*i:3*i+3, 3*i:3*i+3]=self.Jxab+self.Jxba
            self.Jx[3*i:3*i+3, 3*(i+1):3*(i+1)+3]=-self.Jxab
            self.Jx[3*i:3*i+3, 3*(i-1):3*(i-1)+3]=-self.Jxba 
#            self.Jv[3*i:3*i+3, 3*i:3*i+3]=self.dampc*np.matrix(np.identity(3))
#            self.Jv[3*i:3*i+3, 3*(i+1):3*(i+1)+3]=self.dampc*np.matrix(np.identity(3))
#            self.Jv[3*i:3*i+3, 3*(i-1):3*(i-1)+3]=self.dampc*np.matrix(np.identity(3))
#
    
    def packStates(self):
        ''''Packs the velocity and force vectors to solve Adv=D system, calculates the A and D matrices'''
        self.vpack=self.vel.reshape((self.nparticles*3, 1), order='F')
        self.fpack=self.fmag.reshape((self.nparticles*3, 1), order='F')
        self.B=-(self.dt*self.Jv+self.dt*self.dt*self.Jx)
        self.A=np.identity(3*self.nparticles)+np.linalg.inv(self.M)*self.B
        self.D=self.dt*np.linalg.inv(self.M)*(self.fpack+self.dt*self.Jx*self.vpack)
        
    def solveSys(self):
        '''Solves for deltaV and deltaX by using the BICGStab solver from SciPy'''
        self.deltaV=bicgstab(self.A, self.D, tol=1e-10)[0].reshape((3, self.nparticles), order='F')
        self.deltaX=(self.vel+self.deltaV)*self.dt
        
    def updateStates(self):
        '''Updates the velocity and position matrices, feeds the new position to the 
        tether visualization instance'''

        if self.incl_final:  
            self.vel[:, 1:]+=self.deltaV[:, 1:]
            self.pos[:, 1:]+=self.deltaX[:, 1:]
        else:
            self.vel[:, 1:self.nparticles-1]+=self.deltaV[:, 1:self.nparticles-1]
            self.pos[:, 1:self.nparticles-1]+=self.deltaX[:, 1:self.nparticles-1]  
        self.teth.pos=vectorize(self.pos)   
        
    def updatePS(self, frame):
        '''Update after each dt to obtain the new positions '''
        if self.incl_final:
            final=self.nparticles
        else:
            final=self.nparticles-1
        for i in range(1, final):
            self.getDiff(i)
            self.getNorm(i)
            self.getExForces(i, frame)
            self.getImForces(i)
            self.getJx(i)
        self.packStates()
        self.solveSys()
        self.updateStates()
        
    def reelOut(self):
        self.pos[:, i]+=self.reelincr*(self.pos[:, i]/np.linalg.norm(self.pos[:, i]))
            #DO it only for the first particle, and then change the restitution length of the 1st spring to the new value
        self.M[0:3, 0:3]=(0.5*self.mass*(1+np.linalg.norm(self.pos[:, 1])/self.R))*(np.identity(3))
        if ((np.linalg.norm(self.pos[:, 1])/self.R)-1)>=self.Xfactor:
            self.addParticle()
             
    def addParticle(self):
        self.nparticles+=1
        self.pos=np.hstack((self.orig, self.pos))
        self.pos[:, 1]=(self.Xfactor/(self.Xfactor+1))*self.pos[:, 2]
        vproj=(np.dot(self.vel[:, 2], self.pos[:, 2])/(np.linalg.norm(self.pos[:, 2])*np.linalg.norm(self.pos[:, 2])))*self.pos[:, 2]
        vorth=self.vel[:, 2]-vproj
        self.vel=np.hstack((self.orig, self.vel))
        self.vel[:, 1]=vproj+(self.Xfactor/(self.Xfactor+1))*vorth        
        
        tmp1=np.zeros((3, np.shape(self.Jx)[1]))
        tmp2=np.zeros((self.nparticles, 3))
        self.Jx=np.vstack((tmp1, self.Jx))
        self.Jx=np.hstack((tmp2, self.Jx))
        
        self.M=np.vstack((tmp1, self.M))           
        self.M=np.hstack((tmp2, self.M))
        self.M[0:3, 0:3]=np.identity(3)*self.mass*0.5
        self.M[3:6, 3:6]=np.identity(3)*self.mass
        
        print 'PARTICLE ADDED!'
        
    def reelIn(self):
        pass
 
    def plotCatenary(self):
        sns.set(context='paper', font_scale=2.5, style='white')
        plt.plot(-150+np.array(self.pos[0])[0], np.array(self.pos[1])[0], 'k', label='Simulated tether', linewidth=5)
        self.shapepar=abs(np.array(self.fspring[0])[0][1])/(self.rhot*9.81*self.t_radius*self.t_radius*np.pi)       
        self.y=self.shapepar*np.cosh((-150+np.array(self.pos[0])[0])/self.shapepar)
        self.y-=self.y[0]        
        plt.plot(-150+np.array(self.pos[0])[0], self.y, 'r--', label='Catenary equation', linewidth=5)
        plt.xlabel('x [m]', fontsize=25)
        plt.ylabel('y [m]', fontsize=25)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.legend()
        plt.show()
        print "RMSE =", np.sqrt(np.sum((self.y-np.array(self.pos[1])[0])**2)/len((self.y-np.array(self.pos[1])[0])))
    def plotWind(self):
        sns.set(context='paper', font_scale=4.5, style='white')
        plt.plot(np.array(self.pos[0])[0], np.array(self.pos[1])[0], 'k', label='Simulated tether', linewidth=5, marker='o', markersize=12)
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
#        plt.legend()
        plt.show()
        
    def startSim(self):
        '''Method to start running the simulation within the particle system only'''
        self.tr=0      
        
        self.interval=1. 
        self.count=0
        running=True
        self.tmat=[]
        while running:
            self.fmag=np.matrix(np.zeros((3,self.nparticles)))
            self.t=time.time()
            self.tr+=self.dt
            self.tmat.append(self.tr)
            self.updatePS(self.tr/self.dt)
            self.max=np.max(abs(self.pos[1]))
            rate(self.framerate)
            if self.tr>10:
                self.plotCatenary()
                running=False
scr=disp()
scr.create()
ps=ImParticleSystem(screen=scr.window)
ps.startSim()