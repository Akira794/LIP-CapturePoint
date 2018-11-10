#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np              # Numerical library
from scipy import *             # Load the scipy functions
from control.matlab import *    # Load the controls systems library
from matplotlib import pyplot as plt

#ヒューマノイドロボット p108~p109 p111~-112

Zc = 0.8
g = 9.81

def drange(begin, end, step):
    n = begin
    while n + step < end:
        yield n
        n += step

class dynamics():

    T = 0
    E_lip_x = 0.0
    E_lip_y = 0.0
    cog_list_x = []
    cog_list_y =[]
    cog_list_t = []

    def __init__(self, t_sup, dt, init_x, init_dx):

        self.__t_sup = t_sup
        self.__dt    = dt
        self.__omega = np.sqrt(g/Zc)
        self.T = 0
        self.__x, self.__y = init_x, 0.01
        self.__dx, self.__dy = init_dx, 0
        self.__xi, self.__yi  = self.__x, self.__y
        self.__dxi, self.__dyi = self.__dx, self.__dy
        self.__px, self.__py = 0, 0
        self.__pxa, self.__pya = self.__px, self.__py

    def Integrate(self, count):
        self.__count = count
        for t in drange(0.0,self.__t_sup,self.__dt):
            self.x  = (self.__xi-self.__pxa)*cosh(self.__omega*t) + (1.0/self.__omega)*self.__dxi*sinh(self.__omega*t) + self.__pxa
            self.y  = (self.__yi-self.__pya)*cosh(self.__omega*t) + (1.0/self.__omega)*self.__dyi*sinh(self.__omega*t) + self.__pya
            self.dx = self.__omega*(self.__xi-self.__pxa)*sinh(self.__omega*t) + self.__dxi*cosh(self.__omega*t)
            self.dy = self.__omega*(self.__yi-self.__pya)*sinh(self.__omega*t) + self.__dyi*cosh(self.__omega*t)

            self.cog_list_x.append(self.x)
            self.cog_list_y.append(self.y)
            self.cog_list_t.append(t)

            self.E_lip_x = (1.0/2.0)*pow(self.__dxi,2) - (1.0/2.0)*pow(self.__omega,2)*pow(self.__xi - self.__pxa,2)
            self.E_lip_y = (1.0/2.0)*pow(self.__dyi,2) - (1.0/2.0)*pow(self.__omega,2)*pow(self.__yi - self.__pya,2)


        self.T += self.__t_sup
        self.__xi, self.__yi = self.x, self.y
        self.__dxi, self.__dyi = self.dx, self.dy


    def plot_gait_pattern_list(self):
        plt.xlim([0.0,0.8])
        plt.ylim([-0.2,0.2])
        plt.plot(self.cog_list_t, self.cog_list_x, label="$cog$")
        plt.legend()
        plt.grid()
        plt.show()

def main():
    EOM = dynamics(0.9,0.02, -0.2, 0.791)#-0.2, 0.791 # -0.151, 0.467
    for n in range(1):#6
        EOM.Integrate(n)

    print("E_lip_x: " + str(EOM.E_lip_x))
    if EOM.E_lip_x <= 0:
        print("x_apex = " + str(np.sqrt((-2*Zc*EOM.E_lip_x)/g)))
    else:
        print("dx_apex = " + str(np.sqrt(2*EOM.E_lip_x)))

    EOM.plot_gait_pattern_list()

if __name__ == '__main__':
    main()
