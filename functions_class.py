#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 19:03:01 2021

@author: vashisth
"""
# %%

import numpy as np
import math
from math import sqrt
import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from numba import jit
from ipywidgets.widgets.interaction import show_inline_matplotlib_plots
from IPython.display import clear_output
from ipywidgets import FloatProgress
from IPython.display import display


# %%
# define the particle class:
class particle:
    m=1.6726219e-27 #mass (kg)
    rho=0 #density (kg/m^3)
    P=0 #pressure (Pa)
    x=0 #location (m)
    v=0 #speed (m/s)
    rho_0=0 #reference density (kg/m^3)
    K0=0 #reference adiabat 
    alpha=0 #smoothing width (m)
    aindex=5/3.


# W (kernel function) and dW

def W(ri,alpha=1): 
    r=abs(ri/alpha)    
    y=0
    if (r<=1):
        y=21./2./np.pi*(1-r)**4*(1+4*r)
    return y/alpha**3

def dW(ri,alpha=1):
    r=abs(ri/alpha)
    y=0
    if (r<=1):
        y=210./np.pi*(r-1)**3*r
    return y/alpha**4

# rho, pressure, and gradient
def compute_rho(particles):
    for i in range (len(particles)):
        particles[i].rho=0
        for j in range(len(particles)):
            r=particles[i].x-particles[j].x
            r=(r[0]**2+r[1]**2+r[2]**2)**.5
            particles[i].rho+=W(r,particles[j].alpha)*particles[j].m

def compute_pressure(particles):
    for i in range (len(particles)):
        particles[i].P=particles[i].K0*particles[i].rho**particles[i].aindex

def compute_gradient(particles):
    grad=np.zeros((len(particles),len(particles),3))
    for i in range (len(particles)):
        for j in range (len(particles)):
            r=particles[i].x-particles[j].x
            r_val=(r[0]**2+r[1]**2+r[2]**2)**.5
            if (r_val>0):
                r/=r_val
                grad[i,j]=dW(r_val,particles[j].alpha)*r
    return grad            

# dv/ dt

def time_advance(particles,grad,dt=1,x_adv=1):
    global G,nu
    for i in range (len(particles)):
        for j in range(len(particles)):
            val=(particles[j].m*(particles[i].P/particles[i].rho**2+particles[j].P/particles[j].rho**2))*grad[j,i]
            val+=G*particles[j].rho*(particles[j].x-particles[i].x)
            particles[i].v+=val*dt/2.
        particles[i].v-=nu*particles[i].v*dt/2.
    if (x_adv==1):
        for i in range (len(particles)):
            particles[i].x+=particles[i].v*dt
            

def recenter(particles):
    x_ave=np.zeros(3)
    m=0
    for i in range (len(particles)):
        x_ave+=particles[i].x
        m+=particles[i].m
    x_ave/=m
    for i in range (len(particles)):
        particles[i].x-=x_ave

# min dt
def get_min_dt(particles):
    global length
    dt_min=1e10
    for i in range (len(particles)):
        vel=(particles[i].v[0]**2+particles[i].v[1]**2+particles[i].v[2]**2)**.5+1e-3
        dt=length/vel/100
        if (dt_min>dt):
            dt_min=dt
    return min(dt,.1)

# plot 3d and 2d 

def plot_3D(particles,plot_size=3,plot_axes=True,title=''):
    from mpl_toolkits.mplot3d import Axes3D
    n = len(particles)
    #coordinates
    data_x=np.zeros(n)
    data_y=np.zeros(n)
    data_z=np.zeros(n)
    data_m=np.zeros(n)
    for i in range(n):
        data_x[i]=cp[i].x[0]
        data_y[i]=cp[i].x[1]
        data_z[i]=cp[i].x[2]
        data_m[i]=cp[i].m
    #generate figure and axes
    fig = plt.figure(figsize=(plot_size,plot_size))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(data_x, data_y, data_z, c=data_m, linewidths=0,
               marker='o', s=8*plot_size, cmap=plt.cm.Wistia,alpha=.8) #s is the size of the plotted symbol 'o'
    #set autoscale parameters
    xc=(data_x.max()+data_x.min())/2.
    x_low=xc-(data_x.max()-data_x.min())/2.*1.1-1e-12
    x_high=xc+(data_x.max()-data_x.min())/2.*1.1+1e-12
    yc=(data_y.max()+data_y.min())/2.
    y_low=yc-(data_y.max()-data_y.min())/2.*1.1-1e-12
    y_high=yc+(data_y.max()-data_y.min())/2.*1.1+1e-12
    zc=(data_z.max()+data_z.min())/2.
    z_low=zc-(data_z.max()-data_z.min())/2.*1.1-1e-12
    z_high=zc+(data_z.max()-data_z.min())/2.*1.1+1e-12
    #set autoscale parameters
    ax.set_xlim(min(x_low,y_low,z_low),max(x_high,y_high,z_high))
    ax.set_ylim(min(x_low,y_low,z_low),max(x_high,y_high,z_high))
    ax.set_zlim(min(x_low,y_low,z_low),max(x_high,y_high,z_high))
    ax.set_box_aspect((1,1,1))
    if (plot_axes):#so we can switch the axis on or off
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_zlabel("z (m)")
        ax.grid(False) 
        ax.w_xaxis.pane.fill = False
        ax.w_yaxis.pane.fill = False
        ax.w_zaxis.pane.fill = False
    else:
        ax.set_axis_off()
    fig.set_facecolor('black')
    ax.set_facecolor('black')
    plt.suptitle(title)
    plt.show()

