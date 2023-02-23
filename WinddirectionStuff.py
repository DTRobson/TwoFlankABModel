# -*- coding: utf-8 -*-
"""
Created on Tue May  3 12:44:52 2022

@author: Robbo
"""

import numpy as np
import matplotlib.pyplot as plt



def GetVector(windangle_deg):
    rtheta = np.radians(windangle_deg)
    mat = np.array([[np.cos(rtheta), -np.sin(rtheta)],[(np.sin(rtheta)), np.cos(rtheta)]])
    unitvec = np.array([1,0])
    windirection = mat.dot(unitvec)
    windirection = windirection/(windirection.dot(windirection))
    return np.array(windirection)

def GetAngle(wind_direction):
    x, y = wind_direction
    if x == 0:
        if y < 0:
            theta = 3*np.pi/2
        if y > 0:
            theta = np.pi/2
    else:
        if  x > 0:
            theta = np.arctan(y/x)
        else:
            theta = np.pi + np.arctan(y/x)
            
    return theta
    

def Uniform(params):
    mintheta, maxtheta = params
    theta = np.random.uniform(mintheta, maxtheta)
    windirection = GetVector(theta)
    return windirection

def Gaussian(params):
    meantheta, thetasd = params
    theta = np.random.normal(meantheta, thetasd)
    windirection = GetVector(theta)
    return windirection

def BimodalGaussians(params):
    theta1, theta1sd, theta2, theta2sd, weight1 = params
    rand = np.random.uniform()
    if rand <= weight1:
        return Gaussian((theta1, theta1sd))
    else:
        return Gaussian((theta2, theta2sd))


def Fixed(params):
    theta = np.random.choice(params)
    windirection = GetVector(theta)
    return windirection

