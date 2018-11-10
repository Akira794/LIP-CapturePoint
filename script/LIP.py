#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np              # Numerical library
from scipy import *             # Load the scipy functions
from control.matlab import *    # Load the controls systems library
from matplotlib import pyplot as plt

Zc = 0.8
g = 9.81

def drange(begin, end, step):
    n = begin
    while n + step < end:
        yield n
        n += step

class LinearInvertedPendlum():

    T = 0
    foot_step_list = []
    cog_list_x = []
    cog_list_y =[]
    p_list_x =[]
    p_list_y =[]
    p_modi_list_x =[]
    p_modi_list_y =[]

    def __init__(self, t_sup, dt, a, b):

        self.__t_sup = t_sup
        self.__dt    = dt
        self.__a     = a
        self.__b     = b

        self.__Tc = np.sqrt(Zc/g)
        self.__D  = self.__a*pow((cosh(self.__t_sup/self.__Tc)-1),2) + self.__b*pow((sinh(self.__t_sup/self.__Tc)/self.__Tc),2)
        self.__S = sinh(0.8/self.__Tc)
        self.__C = cosh(0.8/self.__Tc)
        self.T = 0
        self.__x, self.__y = 0.0, 0.01
        self.__dx, self.__dy = 0, 0
        self.__xi, self.__yi  = self.__x, self.__y
        self.__dxi, self.__dyi = self.__dx, self.__dy
        self.__px, self.__py = 0, 0
        self.__pxa, self.__pya = self.__px, self.__py

    def SetFootStep(self):

        self.foot_step_list.append([0.0, 0.2, 0.0])
        self.foot_step_list.append([0.25, 0.2, 0.0])
        self.foot_step_list.append([0.25, 0.2, 0.0])
        self.foot_step_list.append([0.25, 0.2, 0.0])
        self.foot_step_list.append([0.0, 0.2, 0.0])
        self.foot_step_list.append([0.0, 0.0, 0.0])
        self.foot_step_list.append([0.0, 0.0, 0.0])

    def Integrate(self, count):
        self.__count = count
        for t in drange(0.0,self.__t_sup,self.__dt):
            self.x =  (self.__xi-self.__pxa)*cosh(t/self.__Tc) + self.__Tc*self.__dxi*sinh(t/self.__Tc) + self.__pxa
            self.y =  (self.__yi-self.__pya)*cosh(t/self.__Tc) + self.__Tc*self.__dyi*sinh(t/self.__Tc) + self.__pya
            self.dx = (self.__xi-self.__pxa)/self.__Tc*sinh(t/self.__Tc) + self.__dxi*cosh(t/self.__Tc)
            self.dy = (self.__yi-self.__pya)/self.__Tc*sinh(t/self.__Tc) + self.__dyi*cosh(t/self.__Tc)
            self.cog_list_x.append(self.x)
            self.cog_list_y.append(self.y)

        self.T += self.__t_sup
        self.__xi, self.__yi = self.x, self.y
        self.__dxi, self.__dyi = self.dx, self.dy

    def CalcLegLandingPos(self):
        self.__px = self.__px + self.foot_step_list[self.__count][0]
        self.__py = self.__py -(pow(-1, ((self.__count) + 1) )) * self.foot_step_list[self.__count][1]

        self.p_list_x.append(self.__px)
        self.p_list_y.append(self.__py)

    def CalcWalkFragment(self):
        self.__xb = self.foot_step_list[self.__count+1][0]/2
        self.__yb = pow(-1,self.__count+1)*self.foot_step_list[self.__count+1][1]/2;
        self.__vxb = ((cosh(self.__t_sup/self.__Tc)+1)/(self.__Tc*sinh(self.__t_sup/self.__Tc)))*self.__xb;
        self.__vyb = ((cosh(self.__t_sup/self.__Tc)-1)/(self.__Tc*sinh(self.__t_sup/self.__Tc)))*self.__yb;

    def CalcGoalState(self):
        self.__xd = self.__px + self.__xb
        self.__dxb = self.__vxb
        self.__yd = self.__py + self.__yb
        self.__dyb = self.__vyb

    def ModifyLandPos(self):
        self.__pxa = -self.__a*(self.__C-1)/self.__D*(self.__xd - self.__C*self.__xi - self.__Tc*self.__S*self.__dxi) - self.__b*self.__S/(self.__Tc*self.__D)*(self.__dxb - self.__S/self.__Tc*self.__xi - self.__C*self.__dxi)
        self.__pya = -self.__a*(self.__C-1)/self.__D*(self.__yd - self.__C*self.__yi - self.__Tc*self.__S*self.__dyi) - self.__b*self.__S/(self.__Tc*self.__D)*(self.__dyb - self.__S/self.__Tc*self.__yi - self.__C*self.__dyi)

        self.p_modi_list_x.append(self.__pxa)
        self.p_modi_list_y.append(self.__pya)

    def plot_gait_pattern_list(self):
        plt.plot(self.cog_list_x, self.cog_list_y, label="$cog$")
        plt.plot(self.p_modi_list_x,self.p_modi_list_y, 'r*', label="$p modi$")
        plt.legend()
        plt.show()

def main():
    LIP = LinearInvertedPendlum(0.8,0.02,10,1) #t_sup dt 評価関数 10 1
    LIP.SetFootStep()

    LIP.p_list_x.append(0.0)
    LIP.p_list_y.append(0.0)

    for n in range(6):
        LIP.Integrate(n)
        LIP.CalcLegLandingPos()
        LIP.CalcWalkFragment()
        LIP.CalcGoalState()
        LIP.ModifyLandPos()

    LIP.plot_gait_pattern_list()

if __name__ == '__main__':
    main()
