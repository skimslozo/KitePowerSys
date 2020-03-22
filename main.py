# -*- coding: utf-8 -*-  
"""
Created on Sun May 05 12:30:28 2019

@author: Miks
"""

from visual import *
from atmos import Atmosphere
from libraries import Tx, Ty, Tz
from kite import KiteModel
from ps_impl_euler_new import disp, ImParticleSystem
import time
import numpy as np
import wx
import pygame

class Simulation(Atmosphere, KiteModel, disp, ImParticleSystem):
    
    def __init__(self, showkiteaxes=True, showlift=True, showdrag=True, showwind=True, showinfo=True, keyboard=True, joystick=False):
        self.atm=Atmosphere()
        self.kite=KiteModel()
        self.screen=disp()
        self.screen.create()
        self.ps=ImParticleSystem(screen=self.screen.window)
        self.ps.M[3*self.ps.nparticles-3:3*self.ps.nparticles, 3*self.ps.nparticles-3:3*self.ps.nparticles]+=np.identity(3)*self.kite.mass
        self.showkiteaxes=showkiteaxes
        self.showlift=showlift
        self.showdrag=showdrag
        self.showinfo=showinfo
        self.showwind=showwind
        if self.showinfo:
            self.log=window(menus=False, title='Status', width=500, height=360)
            self.p=self.log.panel
            d=30
            self.count=0
            self.font = wx.Font(15, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, False)
            self.xinfo=wx.StaticText(self.p, pos=(1,0), label='')
            self.yinfo=wx.StaticText(self.p, pos=(1,1*d), label='')
            self.zinfo=wx.StaticText(self.p, pos=(1,2*d), label='')
            self.vappinfo=wx.StaticText(self.p, pos=(1,3*d), label='')
            self.vwindinfo=wx.StaticText(self.p, pos=(1,4*d), label='')
            self.liftinfo=wx.StaticText(self.p, pos=(1,5*d), label='') 
            self.draginfo=wx.StaticText(self.p, pos=(1,6*d), label='')
            self.alphainfo=wx.StaticText(self.p, pos=(1,7*d), label='')
            self.leninfo=wx.StaticText(self.p, pos=(1,8*d), label='')
            self.liftinfo.Font=self.font
            self.draginfo.Font=self.font
            self.xinfo.Font=self.font
            self.yinfo.Font=self.font
            self.zinfo.Font=self.font
            self.vwindinfo.Font=self.font
            self.vappinfo.Font=self.font
            self.alphainfo.Font=self.font
            self.leninfo.Font=self.font
        if keyboard==True:
            self.keyboard=True
            self.joystick=False
        if joystick==True:
            self.keyboard=False
            self.joystick=True
            pygame.init()
            pygame.joystick.init()
            self.joy=pygame.joystick.Joystick(0)
            self.joy.init()
        
    def createAxes(self):
        self.vaGlobal=self.atm.getWind(self.ps.pos[1, self.ps.nparticles-1])-self.ps.vel[:, self.ps.nparticles-1]
        self.kite.va=self.vaGlobal    
        self.kite.c=self.ps.pos[:, self.ps.nparticles-1]-self.ps.pos[:, self.ps.nparticles-2]
        self.kite.zk=-self.kite.c/np.linalg.norm(self.kite.c)
        self.kite.yk=np.matrix(np.cross(self.vaGlobal, self.kite.c, axis=0)/np.linalg.norm(np.cross(self.vaGlobal, self.kite.c, axis=0)))
        self.kite.xk=np.matrix(np.cross(self.kite.yk, self.kite.zk, axis=0))
        self.kite.getForces(rho=self.atm.dens, steering=0, power=0)
        if self.showkiteaxes:
            self.kzax=arrow(pos=self.ps.pos[:, self.ps.nparticles-1], axis=10*self.kite.zk, color=(0, 255, 0), maketrail=False, shaftwidth=0.7)
            self.kyax=arrow(pos=self.ps.pos[:, self.ps.nparticles-1], axis=10*self.kite.yk, color=(0, 255, 0), maketrail=False, shaftwidth=0.7)
            self.kxax=arrow(pos=self.ps.pos[:, self.ps.nparticles-1], axis=10*self.kite.xk, color=(0, 255, 0), maketrail=False, shaftwidth=0.7)
        if self.showlift:
            self.liftax=arrow(pos=self.ps.pos[:, self.ps.nparticles-1], axis=10*self.kite.L/np.linalg.norm(self.kite.D), color=(0, 0, 255), maketrail=False)
        if self.showdrag:
            self.dragax=arrow(pos=self.ps.pos[:, self.ps.nparticles-1], axis=10*self.kite.D/np.linalg.norm(self.kite.D), color=(255, 0, 0), maketrail=False)
        if self.showwind:
            self.vwindinit=np.linalg.norm(self.atm.velocity)
            self.windax=arrow(pos=[0, 0, 0], axis=15*self.atm.velocity/self.vwindinit, color=(51, 102, 0), maketrail=False)
            
    def updateKiteAxes(self):
        self.vaGlobal=self.atm.getWind(self.ps.pos[1, self.ps.nparticles-1])-self.ps.vel[:, self.ps.nparticles-1]
        self.kite.va=self.vaGlobal    
        self.kite.c=self.ps.pos[:, self.ps.nparticles-1]-self.ps.pos[:, self.ps.nparticles-2]
        self.kite.zk=-self.kite.c/np.linalg.norm(self.kite.c)
        self.kite.yk=np.matrix(np.cross(self.vaGlobal, self.kite.c, axis=0)/np.linalg.norm(np.cross(self.vaGlobal, self.kite.c, axis=0)))
        self.kite.xk=np.matrix(np.cross(self.kite.yk, self.kite.zk, axis=0))
        if self.showkiteaxes:
            self.kxax.pos=self.ps.pos[:, self.ps.nparticles-1]
            self.kxax.axis=10*self.kite.xk
            self.kyax.pos=self.ps.pos[:, self.ps.nparticles-1]
            self.kyax.axis=10*self.kite.yk
            self.kzax.pos=self.ps.pos[:, self.ps.nparticles-1]
            self.kzax.axis=10*self.kite.zk
        if self.showwind:
            self.windax.axis=10*self.atm.vwind/self.vwindinit
        
    def updateKiteAngles(self):
        self.kite.vaxz=self.vaGlobal-np.dot(np.transpose(self.vaGlobal), self.kite.yk)[0, 0]*self.kite.yk
        self.kite.vaxy=self.vaGlobal-np.dot(np.transpose(self.vaGlobal), self.kite.zk)[0, 0]*self.kite.zk
        self.kite.alphatop=np.arccos(np.dot(np.transpose(self.kite.vaxz), self.kite.xk)[0, 0]/(np.linalg.norm(self.kite.vaxz)*np.linalg.norm(self.kite.xk)))
        self.kite.betatop=np.arccos(np.dot(np.transpose(self.kite.vaxy), self.kite.xk)[0, 0]/(np.linalg.norm(self.kite.vaxy)*np.linalg.norm(self.kite.xk)))
    
    def getKeyboardInput(self):
        if self.keyboard:
            if self.screen.window.kb.keys:
                for i in range(self.screen.window.kb.keys):
                    self.key=self.screen.window.kb.getkey()
                    if self.key=='left':
                        self.steering=-1.0
                    elif self.key=='right':
                        self.steering=1.0
                    elif self.key=='up':
                        self.power=1.0
                    elif self.key=='down':
                        self.power=0.0
                    else:
                        self.steering=0.
                        self.power=0.5
                        if self.key=='d':
                            self.ps.reelOut()
                        if self.key=='a':
                            self.ps.reelIn()
            else:
                self.steering=0.
                self.power=0.5

    def getSteeringInput(self):
        if pygame.event.get()!=[] or self.joy.get_axis(0)!=0:
            return self.joy.get_axis(0)
        else:
            return 0
            
    def getPowerInput(self):
        if self.joystick:
            if pygame.event.get()!=[] or self.joy.get_axis(1)!=0:
                return self.joy.get_axis(1)
            else:
                return 0
        else:
            return 0
           
    def increaseWind(self):
        if self.joystick and self.joy.get_hat(0)!=(0, 0):
            self.atm.velocity[0, 0]+=-0.1*self.joy.get_hat(0)[1]
            self.atm.velocity[2, 0]+=-0.1*self.joy.get_hat(0)[0]
            
    def updateGlobeAngles(self):
        """Implement the angle definitions"""
        self.theta=np.arctan2(np.sqrt(self.ps.pos[0, self.ps.nparticles-1]**2+self.ps.pos[2, self.ps.nparticles-1]**2), self.ps.pos[1, self.ps.nparticles-1])
        self.phi=np.arctan2(self.ps.pos[0, self.ps.nparticles-1]**2, self.ps.pos[2, self.ps.nparticles-1]**2)
        self.increaseWind()
        
    def getKiteToGlobe(self, vector):
        """Set-up the transformation matrix, return the transformed vectors"""
        self.TGk=Tx(-np.pi/2)*Tz(-((np.pi/2)-self.phi))*Ty(self.theta)
        
    def addKiteForces(self):
        self.ps.fmag=np.matrix(np.zeros((3,self.ps.nparticles)))
        self.ps.fmag[:, self.ps.nparticles-1]+=self.kite.L+self.kite.D+self.kite.S
        if self.showlift:
            self.liftax.pos=self.ps.pos[:, self.ps.nparticles-1]
            self.liftax.axis=10*self.kite.L/np.linalg.norm(self.kite.L)
        if self.showdrag:
            self.dragax.pos=self.ps.pos[:, self.ps.nparticles-1]
            self.dragax.axis=10*self.kite.D/np.linalg.norm(self.kite.D)

    def displayInfo(self):
        if self.tr>self.count:
            self.count+=1
            self.liftinfo.LabelText='Lift: {:.2f} N'.format(np.linalg.norm(self.kite.L))
            self.draginfo.LabelText='Drag: {:.2f} N'.format(np.linalg.norm(self.kite.D))
            self.xinfo.LabelText='Xpos: {:.2f} m'.format(self.ps.pos[0, self.ps.nparticles-1])
            self.yinfo.LabelText='Ypos: {:.2f} m'.format(self.ps.pos[1, self.ps.nparticles-1])
            self.zinfo.LabelText='Zpos: {:.2f} m'.format(self.ps.pos[2, self.ps.nparticles-1])
            self.vappinfo.LabelText='Apparent velocity: {:.2f} m/s'.format(np.linalg.norm(self.vaGlobal))
            self.vwindinfo.LabelText='Wind velocity: {:.2f} m/s'.format(np.linalg.norm(self.atm.velocity))
            self.alphainfo.LabelText='Angle of attack: {:.2f} deg'.format(np.degrees(self.kite.alphatop))
            self.leninfo.LabelText='Length {:.2f} m'.format(np.linalg.norm(self.ps.pos[:, self.ps.nparticles-1]))
            
    def startSim(self):
        self.tr=0.
        self.running=True
        self.createAxes()
        while self.running:
            self.t=time.time() 
            self.tr+=self.ps.dt
            self.ps.fmag=np.matrix(np.zeros((3,self.ps.nparticles)))
            self.updateGlobeAngles()
            self.updateKiteAxes()
            self.updateKiteAngles()
            if self.keyboard:
                self.getKeyboardInput()
                self.kite.getForces(rho=self.atm.dens, power=self.power, steering=self.steering)
            if self.joystick:
                self.kite.getForces(rho=self.atm.dens, power=self.getPowerInput(), steering=self.getSteeringInput())
            self.addKiteForces()
            self.ps.updatePS(self.tr/self.ps.dt)
#            self.screen.window.kb.getkey() 
            if self.showinfo:
                self.displayInfo()
            rate(self.ps.framerate)


sim=Simulation()
sim.startSim()