# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 10:13:50 2023

@author: Robbo
"""
#!!!!!!!!!!!!!!!!!!!!!!
#A note on packages, you will need to use shapely 1.8.0 as the TFABM modules are currently incompatible with later versions

import numpy as np
import matplotlib.pyplot as plt
import ABModel as abm
import WinddirectionStuff as wind
from shapely.geometry import Polygon, Point, box
from GammaStuff import max_gamma as gamma_c
from multiprocessing import Pool, cpu_count
import datetime
from FluxAndAvalanching import Barchan, Rotate

def Uniform(params, num): #an example function to generate 
    return np.random.uniform(params[0], params[1], num)

                 
                 
####CAVE!!!
##########Note! In these modules, the primary wind direction is always assumed to be in the negative y-direction i.e. a bearing of 270 degrees
#Thus the central axes of dunes is aligned with the y-axis and the widths of dunes are their widths parallel to the x-axis

simwidth = 600 #set the maximum width of the simulated space
simlength = 600 #set the maximum length of the simulated space
fieldwidth = 200 #set the width of the actual dune swarm.  If dunes are to be injected they will have 
#x-coordinate in the range [(simwidth - fieldwidth)/2, (simwidth - fieldwidth)/2 + fieldwidth]
#i.e. the field in which dunes can be injected is situated in the middle laterally of the  
fieldlength = 200 #the length of the field, dunes migrating beyond this will be destroyed.
xs = [276., 335.] #the initial x-coordinates of dunes (remember this is lateral dimension i.e. widths of dunes not lengths!)
ys = [580., 570.] #the initial y-coordinates of dunes (remember the model assumes the 
                                                #primary wind causes migration in the negative y-direcion)
lws = [20., 15.4] #the initial left widths of dunes(i.e. starboard)
rws = [20., 15.4] #the initial right widths of dunes (i.e. port)
qsatinit = 100. #the initial value of qsat, used to calculate how many dunes will be injected in the first timestep
                # if injections are turned on
q0 = 0.2       #the ratio of q0 to qsat, can be adjusted during a run with model.q0ratio = q0new
qshift = 0.2   #the ratio of qshift to qsat, can be adjusted during a run with model.qshiftratio = qshiftnew
dt = 0.001    #the length of a timestep, the unit is not important so long as it matches the inverse time unit of qsat
                # e.g. if qsat is in m^2/year then dt must be in year, if qsat is in cm^2/sec then dt must be in seconds
                
######The remaining parameters are optional and will take default values if not included#####

collson = True #if False collisions will be turned off

inject = False #if True, a number of dunes will be injected into the swarm at every timestep.  This can be achieved manually instead
                # using swarm.add(lwnew, rwnew, xnew, ynew, newid = model.maxid + 1)
injectdist = Uniform #function for determining the sizes of newly injected dunes
                     #can be set to any function you want to define that takes a two input, params and num
                     # e.g Uniform takes two inputs e.g. params = (a,b,c,d), num =  2 to give
                     #an output (lw, rw) where lw is uniform in [a,b] and rw
                     # is uniform in [c,d]
injectparams = [15.3, 15.3, 15.3, 15.3]#the input parameter for the chosen injectdist
initdensity = 1e-5 #used to define the density of dunes just upwind of the simulated space which controls the rate at which
                    #they will be injected if inject == True
                
periodic = False #if True boundaries will be periodic

lambda1 = 1. #shape parameters for the dunes
lambda2 = 1.8
lambda3 = 1/6
alpha = 0.05
delta = 4.6

c = 50. #the wind speed-up governing the rate of migration of the dunes

w0 = 0. #appears in the migration rate vmig = c*qsat/(wtot + w0)

outfluxmode = 'Hersen' #sets the method for calculating outfluxes, Hersen uses the standard assumption that output
                        #is at saturate rate across horn width whorns = alpha*wtot + delta
                        # appears in Hersen, Pascal, and Stéphane Douady. "Collision of barchan dunes as a mechanism of size regulation." Geophysical Research Letters 32.21 (2005).
                        # Can be replaced with 'Duran' which uses the rule outlined in Durán, Orencio, et al. "Size distribution and structure of Barchan dune fields." Nonlinear Processes in Geophysics 18.4 (2011): 455-467.
a = 0.4 #the parameters if outfluxmode == 'Duran'
b = 1 #see Durán, Orencio, et al. "Size distribution and structure of Barchan dune fields." Nonlinear Processes in Geophysics 18.4 (2011): 455-467.
        #for details
        
plottinghornflux = False #if True, model.step will provide the shapely objects describing the flux streaming off of the horns of 
                            # the dunes

keep_coll_rec = False #if True, will save the inputs and outputs of collisions in model.collrecord
                        # each entry will have the order
                        #[(inx1, iny1), (inx2, iny2), (inlw1, inrw1), (inlw2, inrw2), outlws, outrws, outxs, outys])
                        
                        
#######################################################################################################################

#To create the model one then simply needs to run following the line:

s1 = abm.Swarm(simwidth, simlength, fieldwidth, fieldlength, xs, ys, lws, rws, qsatinit, q0, qshift, dt,
                  collson = collson, inject = inject, injectdist = injectdist, injectparams = injectparams,
                  initdensity = initdensity, periodic = periodic, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3,
                  alpha = alpha, delta = delta, c = c, w0 = w0, outfluxmode = outfluxmode, a = a, b = b,
                  plottinghornflux = plottinghornflux, keep_coll_rec = keep_coll_rec)

#Can extract the coordinates and widths as follows
print('xs=', s1.xs, ',ys=', s1.ys, ',total widths=', s1.lws + s1.rws, 'initially')
print()
#Progressing a simulation is achieved using command model.step which takes 4 arguments:
    #1) a value of qsat for that time interval
    #2) a wind direction for that time interval
    #3) a shapely polygon describinig what shape the flux field would take if no dunes were present
    #4) a boolean indicating wether to return geometries to make plotting easier

#A single iteration is shown below
qsat_iter = 100.
winddirection = [0,-1] #a vector in the negative y-direction, this should be your wind direction for most iterations
fluxfield = box(s1.minx, s1.miny, s1.maxx, s1.simlength) #the easiest choice, a box full filling the field
helpplotting = False

#we then run the iteration simply with 
s1.step(qsat_iter, winddirection, fluxfield, helpplotting)

#we could check to see that things have updated
print('xs=', s1.xs, ',ys=', s1.ys, ',total widths=', s1.lws + s1.rws, 'after 1 iteration')
print()
#if we wanted to run another iteration with the same parameters but this time we will want to plot the dunes afterwards
#we can do that by changing helpplotting to True

helpplotting = True

#and running another step but this time we will have some outputs

leftflanks, rightflanks, updatedfluxfield, lefthornfluxes, righthornfluxes = s1.step(qsat_iter,
                                                                                      winddirection, fluxfield, helpplotting)

#note that lefthornfluxes and righthornfluxes will be empty unless we set plottinghornflux to True
#if we want to change that during a run  we simply use s1.plottinghornflux = True
#leftflanks and rightflanks will be lists of shapely polygons while updatedfluxfield will either be a shapely polygon
#or a shapely MultiPolygon
#A simple way of plotting the swarm after the second iteration would be

plt.figure()
for l in leftflanks:
    plt.fill(*l.exterior.xy, color = 'blue')
for r in rightflanks:
    plt.fill(*r.exterior.xy, color = 'red')
  

plt.fill(*updatedfluxfield.exterior.xy, color = 'yellow')
plt.show()

#let's see what happens if we do that again with some horn fluxes shown
#first we need to turn on hornflux plotting

s1.plottinghornflux = True

#now we just do the same as we did before

leftflanks, rightflanks, updatedfluxfield, lefthornfluxes, righthornfluxes = s1.step(qsat_iter,
                                                                                      winddirection, fluxfield, helpplotting)
#and plot
plt.figure()
for l in leftflanks:
    plt.fill(*l.exterior.xy, color = 'blue')
for r in rightflanks:
    plt.fill(*r.exterior.xy, color = 'red')
for hl in lefthornfluxes:
    plt.fill(*hl.exterior.xy, color = 'green')
for hr in righthornfluxes:
    plt.fill(*hr.exterior.xy, color = 'brown')    


plt.fill(*updatedfluxfield.exterior.xy, color = 'yellow')
plt.show()

#if we want to manually add a dune then we can do that easily with s1.add(lwnew, rwnew, xnew, ynew, newid = swarm.maxid + 1) e.g

s1.add(12., 12., 330., 590, s1.maxid + 1)

#now let's run another iteration and plot again

leftflanks, rightflanks, updatedfluxfield, lefthornfluxes, righthornfluxes = s1.step(qsat_iter,
                                                                                      winddirection, fluxfield, helpplotting)

plt.figure()
for l in leftflanks:
    plt.fill(*l.exterior.xy, color = 'blue')
for r in rightflanks:
    plt.fill(*r.exterior.xy, color = 'red')
for hl in lefthornfluxes:
    plt.fill(*hl.exterior.xy, color = 'green')
for hr in righthornfluxes:
    plt.fill(*hr.exterior.xy, color = 'brown')    


plt.fill(*updatedfluxfield.exterior.xy, color = 'yellow')
plt.show()

#now let's say we want to do another iteration but with a slightly different wind direction for example a bearing of 300 degrees
# the WinddirectionStuff module contains lots of helpful functions including, GetVector so we don't need to know what 
#vector that corresponds to
print()
wind300deg = wind.GetVector(300)
print('300 degree vector=', wind300deg)

#let's see what the hornfluxes and fluxfield look like if we run an iteration with this winddirection

leftflanks, rightflanks, updatedfluxfield, lefthornfluxes, righthornfluxes = s1.step(qsat_iter,
                                                                                      wind300deg, fluxfield, helpplotting)

plt.figure()
plt.fill(*updatedfluxfield.exterior.xy, color = 'yellow')
for l in leftflanks:
    plt.fill(*l.exterior.xy, color = 'blue')
for r in rightflanks:
    plt.fill(*r.exterior.xy, color = 'red')
for hl in lefthornfluxes:
    if type(hl) == type(box(0,0,1,1)):
        plt.fill(*hl.exterior.xy, color = 'green')
    else:
        for geom in hl.geoms:
            plt.fill(*geom.exterior.xy, color = 'green')
for hr in righthornfluxes:
    if type(hr) == type(box(0,0,1,1)):
        plt.fill(*hr.exterior.xy, color = 'brown')
    else:
        for geom in hr.geoms:
            plt.fill(*geom.exterior.xy, color = 'green')  


plt.show()

#As a final demonstration, let's imagine we would like to run the simulation for half a year (500 iterations) and plot every
#tenth of a year i.e. every 100 iterations using a wind which is unidirectional

for i in range(500):
    leftflanks, rightflanks, updatedfluxfield, lefthornfluxes, righthornfluxes = s1.step(qsat_iter,
                                                                                      winddirection, fluxfield, helpplotting)
    print(i)
    if i%100 == 0:
        plt.figure()
        plt.fill(*updatedfluxfield.exterior.xy, color = 'yellow')
        for l in leftflanks:
            plt.fill(*l.exterior.xy, color = 'blue')
        for r in rightflanks:
            plt.fill(*r.exterior.xy, color = 'red')
        for hl in lefthornfluxes:
            if type(hl) == type(box(0,0,1,1)):
                plt.fill(*hl.exterior.xy, color = 'green')
            else:
                for geom in hl.geoms:
                    plt.fill(*geom.exterior.xy, color = 'green')
        for hr in righthornfluxes:
            if type(hr) == type(box(0,0,1,1)):
                plt.fill(*hr.exterior.xy, color = 'brown')
            else:
                for geom in hr.geoms:
                    plt.fill(*geom.exterior.xy, color = 'green')  

        plt.show()
    
#and plot the end state too
    
plt.figure()
plt.fill(*updatedfluxfield.exterior.xy, color = 'yellow')
for l in leftflanks:
    plt.fill(*l.exterior.xy, color = 'blue')
for r in rightflanks:
    plt.fill(*r.exterior.xy, color = 'red')
for hl in lefthornfluxes:
    if type(hl) == type(box(0,0,1,1)):
        plt.fill(*hl.exterior.xy, color = 'green')
    else:
        for geom in hl.geoms:
            plt.fill(*geom.exterior.xy, color = 'green')
for hr in righthornfluxes:
    if type(hr) == type(box(0,0,1,1)):
        plt.fill(*hr.exterior.xy, color = 'brown')
    else:
        for geom in hr.geoms:
            plt.fill(*geom.exterior.xy, color = 'green')  
            
plt.show()

#You will notice that during the simulation a collision took place.  If you wanted to look at the results of that collision in
#more detail, you could begin the simulation again but with s1.keep_coll_rec = True
#then the details of the collision would be stored in s1.collrecord






