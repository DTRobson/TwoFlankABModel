3# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 10:48:46 2022

@author: Robbo
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 11:38:01 2022

@author: Robbo
"""

from mesa import Agent, Model
from mesa.time import SimultaneousActivation
from mesa.space import ContinuousSpace
import numpy as np
import FluxAndAvalanching as bss
import CollisionStuffVersion2 as cs
from shapely.ops import unary_union
import matplotlib.pyplot as plt
import multiprocessing as mp
from GammaStuff import max_gamma as gamma_c
import WinddirectionStuff as wind


def EquivWidth(v, lambda2, lambda3, flankorbody):
    '''Function calculates the equivalent width of either a flank or full barchan with volume v'''
    if flankorbody == 'flank':
        w = (6*v/(lambda2*lambda3))**(1/3)
    if flankorbody == 'body':
        w = 2*(3*v/(lambda2*lambda3))**(1/3)
    return w
    


def VS(w, otherw, lambda1, lambda2, lambda3, onlytotal = False):
    '''Function calculates the volume of the body, horn and their total of a flank with width w
    where the other flank has width otherw'''
    if onlytotal == False:
        vb = 1/6 * lambda1 * lambda3 * w**2 * (w + otherw)/2
        vh = 1/6 * lambda3 * w**2 * (lambda2 * w - lambda1*(w + otherw)/2)
        return vb, vh, vb + vh
    else:
        return 1/6 * lambda2 * lambda3 * w**3


def COMS(x, y, wl, wr, lambda1, lambda2, lambda3, alpha = 0.05, delta = 4.6):
    '''Determines the COM of a full barchan and each of its flanks'''
    lbody = lambda1 * (wl + wr)/2
    xlb = x - wl/4
    ylb = y - 3*lbody/4
    xrb = x + wr/4
    yrb = ylb
    xlh = x - wl + (2*wl +alpha*wl/2 + delta/4)/4
    xrh = x + (2*wr + alpha*wr/2 + delta/4)/4
    ylh = y - lbody - (lambda2*wl - lambda1 * lbody)
    yrh = y - lbody - (lambda2*wr - lambda1*lbody)
    vrb, vrh, vrtot = VS(wr, wl, lambda1, lambda2, lambda3)
    vlb, vlh, vltot = VS(wl, wr, lambda1, lambda2, lambda3)
    xl = (xlb*vlb + xlh*vlh)/vltot
    yl = (ylb*vlb + ylh*vlh)/vltot
    xr = (xrb*vrb + xrh*vrh)/vrtot
    yr = (yrb*vrb + yrh*vrh)/vrtot
    x = (xl*vltot + xr*vrtot)/(vltot + vrtot)
    y = (yl*vltot + yr*vrtot)/(vltot + vrtot)
    return (x, y), (xl, yl), (xr, yr)
    
class Barchan(Agent):
    '''Creates new class representing a barchan dune, identified by the a unique_id, a left-width and a right-width.  The dune is situated in 
    a space belonging to the specified model'''
    
    def __init__(self, unique_id, model, lw, rw):
        self.unique_id = unique_id
        self.model = model
        self.lw = lw
        self.rw = rw
        self.leftinflux = 0
        self.rightinflux = 0
        self.tomove = (0, 0)
        self.leftoutflux = 0
        self.rightoutflux = 0
        self.tosplit = False
        self.lwnew = lw
        self.rwnew = rw
        self.destroy = False


        
    def update(self, xnew, ynew, lnew, rnew):
        self.model.desert.move_agent(self, (xnew, ynew))
        ind = np.where(self.model.dune_unique_ids == self.unique_id)
        self.model.xs[ind] = xnew
        self.model.ys[ind] = ynew
        self.model.lws[ind] = lnew
        self.model.rws[ind] =  rnew
        self.lw = lnew
        self.rw = rnew
        self.tomove = (0, 0)
        self.leftoutflux = 0
        self.rightoutflux = 0
        self.leftinflux = 0
        self.rightinflux = 0
        self.tosplit = False
        self.lwnew = lnew
        self.rwnew = rnew


        


    def split(self):
        

        
        lambda1 = self.model.lambda1
        lambda2 = self.model.lambda2
        lambda3 = self.model.lambda3
        alpha = self.model.alpha
        delta = self.model.delta
    
    
        x, y = self.pos
    
        lw = self.lw
        rw = self.rw
        
        
        lv = VS(lw, rw, lambda1, lambda2, lambda3, True)
        rv = VS(rw, lw, lambda1, lambda2, lambda3, True)
        
        (comx, comy), (comxl, comyl), (comxr, comyr) = COMS(x, y, lw, rw, lambda1, lambda2, lambda3, alpha, delta)
    
        which = self.tosplit
        
        if which == 'left':
 
            leftflanktotalwidth = EquivWidth(lv, lambda2, lambda3, 'body')
            newlw = leftflanktotalwidth/2
            newrw = leftflanktotalwidth/2
            self.lw = newlw
            self.rw = newrw
            print('flankvolume conservation', VS(newlw, newrw, lambda1, lambda2, lambda3, True)*2 - lv)
            xnewguess = x - lw + newlw
            ynewguess = y
            
            
            rightflanktotalwidth = EquivWidth(rv, lambda2, lambda3, 'body')
            daughterlw = rightflanktotalwidth/2
            daughterrw = daughterlw
            print(daughterlw/daughterrw)
            print('daughter volume conservation', VS(daughterlw, daughterrw, lambda1, lambda2, lambda3, True)*2 - rv)
            daughterxguess = x + rw - daughterrw
            daughteryguess = ynewguess - newrw * lambda2 - np.random.uniform(0, 0.0001)
            
            (guesscomx, guesscomy), blah1, blah2 = COMS(xnewguess, ynewguess, newlw, newrw, lambda1, lambda2, lambda3, alpha, delta)
            (daughterguesscomx, daughterguesscomy), blah1, blah2 = COMS(daughterxguess, daughteryguess, daughterlw, daughterrw,
                                                                        lambda1, lambda2, lambda2, alpha, delta)
            
            totalguesscomx = (guesscomx * lv + daughterguesscomx * rv)/(lv + rv)
            totalguesscomy = (guesscomy * lv + daughterguesscomy * rv)/(lv + rv)
            
            diffx = comx - totalguesscomx
            diffy = comy - totalguesscomy
            
            xnew = xnewguess + diffx
            ynew = ynewguess + diffy
            
            
            
            if xnew > self.model.minx and xnew < self.model.maxx and ynew < self.model.simlength and ynew > self.model.miny:
                self.model.add(newlw, newrw, xnew, ynew, self.model.maxid + 1)
            
            self.model.remove(self)
                
            
            
            
            daughterx = daughterxguess + diffx
            daughtery = daughteryguess + diffy
            
            self.model.add(daughterlw, daughterrw, daughterx, daughtery,
                           self.model.maxid + 1)
        
        else:
            rightflanktotalwidth = EquivWidth(rv, lambda2, lambda3, 'body')
            newlw = rightflanktotalwidth/2
            newrw = rightflanktotalwidth/2
            self.lw = newlw
            self.rw = newrw
            print(self.lw/self.rw)
            print('flankvolume conservation', VS(newlw, newrw, lambda1, lambda2, lambda3, True)*2 - rv)
            xnewguess = x + rw - newrw
            ynewguess = y
            
            
            leftflanktotalwidth = EquivWidth(lv, lambda2, lambda3, 'body')
            daughterlw = leftflanktotalwidth/2
            daughterrw = daughterlw
            print(daughterlw/daughterrw)
            print('daughter volume conservation', VS(daughterlw, daughterrw, lambda1, lambda2, lambda3, True) * 2 - lv)
            daughterxguess = x - lw + daughterlw
            daughteryguess = ynewguess - newlw * lambda2 - np.random.uniform(0, 0.0001)
            
            (guesscomx, guesscomy), blah1, blah2 = COMS(xnewguess, ynewguess, newlw, newrw, lambda1, lambda2, lambda3, alpha, delta)
            (daughterguesscomx, daughterguesscomy), blah1, blah2 = COMS(daughterxguess, daughteryguess, daughterlw, daughterrw,
                                                                        lambda1, lambda2, lambda2, alpha, delta)
            
            totalguesscomx = (guesscomx * rv + daughterguesscomx * lv)/(lv + rv)
            totalguesscomy = (guesscomy * rv + daughterguesscomy * lv)/(lv + rv)
            
            diffx = comx - totalguesscomx
            diffy = comy - totalguesscomy
            
            xnew = xnewguess + diffx
            ynew = ynewguess + diffy

            
            if xnew > self.model.minx and xnew < self.model.maxx and ynew < self.model.simlength and ynew > self.model.miny:
                self.model.add(newlw, newrw, xnew, ynew, self.model.maxid + 1)
            
            self.model.remove(self)
                

            
            
            daughterx = daughterxguess + diffx
            daughtery = daughteryguess + diffy
            
          
            
            self.model.add(daughterlw, daughterrw, daughterx, daughtery,
                           self.model.maxid + 1)
        
        self.model.numsplits += 1
        print(ynew, daughtery, 'ypositions')
        print(xnew, daughterx, 'xpositions')
        print(newlw, daughterlw, 'lws')
        print('split complete', self.model.lws/self.model.rws, lv + rv, 2* VS(self.lw, self.rw, lambda1, lambda2, lambda3, True) + 2*VS(daughterlw, daughterrw, lambda1, lambda2, lambda3, True))
        print('xs', self.model.xs)
        print('ys', self.model.ys)
        print('lws', self.model.lws)
        
            


            
            
    def step(self):

        
        lwold = self.lw
        rwold = self.rw
        lambda1 = self.model.lambda1
        lambda2 = self.model.lambda2
        lambda3 = self.model.lambda3
        
        lvold = VS(lwold, rwold, lambda1, lambda2, lambda3, True)
        rvold = VS(rwold, lwold, lambda1, lambda2, lambda3, True)

        lvnew = lvold + self.leftinflux - self.leftoutflux

        rvnew = rvold + self.rightinflux - self.rightoutflux
        
        if lvnew <= 0:
            self.lwnew = 0
        
        else:
            self.lwnew = EquivWidth(lvnew, lambda2, lambda3, 'flank')

        if rvnew <= 0:
            self.rwnew = 0
        else:
            self.rwnew = EquivWidth(rvnew, lambda2, lambda3, 'flank')


    def advance(self):
        
        if self.destroy == True:
            self.model.remove(self)
              
        else:
            
            lwold = self.lw
            rwold = self.rw
            lambda1 = self.model.lambda1
            lambda2 = self.model.lambda2
            lambda3 = self.model.lambda3
            alpha = self.model.alpha
            delta = self.model.delta
            
            minw = (delta/2)/(1 - alpha)
            
            x, y = self.pos
            xmig, ymig = x + self.tomove[0], y + self.tomove[1]
            
            xmove, ymove = self.tomove[0], self.tomove[1]
            

            
            if xmig < self.model.minx or xmig > self.model.maxx or ymig < self.model.miny:
                self.model.remove(self)
            
            else:
            
                lnew1 = self.lwnew
                rnew1 = self.rwnew
                wtot = lnew1 + rnew1
                (comx, comy), blah1, blah2 = COMS(xmig, ymig, lnew1, rnew1, lambda1, lambda2, lambda3, alpha, delta)
                
                
                if lnew1 < minw  or rnew1 < minw:
                    lnew = wtot/2
                    rnew = wtot/2
                    if rnew1 <= minw and lnew1 > minw:
                        xnew, ynew = xmig - lnew, ymig - lnew * self.model.lambda1 
                        (newcomx, newcomy), blah3, blah4 = COMS(xnew, ynew, lnew, rnew, lambda1, lambda2, lambda3, alpha, delta)
                        shiftx = comx - newcomx
                        shifty = comy - newcomy
                        xnew += shiftx
                        ynew += shifty
                        self.update(xnew, ynew, lnew, rnew)
                    elif lnew1 <= minw and rnew1 > minw:
                        xnew, ynew = xmig + lnew, ymig - lnew * self.model.lambda1
                        (newcomx, newcomy), blah3, blah4 = COMS(xnew, ynew, lnew, rnew, lambda1, lambda2, lambda3, alpha, delta)
                        shiftx = comx - newcomx
                        shifty = comy - newcomy
                        xnew += shiftx
                        ynew += shifty
                        self.update(xnew, ynew, lnew, rnew)
                    else:
                        self.model.remove(self)
                        
            
                else:
                    
                    maxasymmratio = gamma_c(min(lnew1, rnew1), self.model.alpha, self.model.delta,
                                            self.model.lambda2, self.model.qshiftratio)
                    
                    ynew = ymig
                    xnew = xmig
                    
                    lnew, rnew, xnew = bss.LateralShift(lnew1, rnew1, xmig, lambda1, lambda2, lambda3, self.model.qshift, self.model.dt, self.model.windangle)
                    self.update(xnew, ynew, lnew, rnew)
                    
                    if lnew1 > 0  and rnew1 > 0:
                        
                        if (lnew1/rnew1) > maxasymmratio:
                            self.tosplit = 'left'

                        elif (rnew1/lnew1) > maxasymmratio:
                            self.tosplit = 'right'
                        
                        if self.tosplit != False:
                            self.split()


                            
                    else:
                        self.model.remove(self)
                    


                    
                    

        
        
        
        
        
        

        
    
class Swarm(Model):
    '''Creates a new class representing a swarm of barchan dunes.  These dunes are situated in a continuous space which is a model attribute.  Other
    attributes are the positions, sizes, and id numbers of the dunes and model parameters including whether or not new dunes are injected at the upwind boundary and, if
    they are, the parameters of the injecting size distirbution, the shape parameters of the dunes, and the sand flux parameters.
    whether the space is periodic, the initial density of dunes,
    
    def __init__(self, simwidth, simlength, fieldwidth, fieldlength,
                 xs, ys, lws, rws, qsatinit,
                 q0, qshift, dt,
                 inject = False, periodic = False, initdensity = None,
                 injectparams = [2.3, 18, 2.3, 18],
                 lambda3 = 1/6, lambda1 = 1., lambda2 = 1.8, thetar = np.pi/6,
                 alpha = 0.05, delta = 4.6, collson = False, c = 1., w0 = 0., a = 0.45, b = 0.1, outfluxmode = 'Hersen'):
        
        

        
        self.xs = xs
        self.ys = ys
        self.lws = lws
        self.rws = rws
        self.simwidth = simwidth
        self.simlength = simlength + 10
        self.fieldwidth = fieldwidth
        self.fieldlength = fieldlength
        self.minx = (self.simwidth - self.fieldwidth)/2
        self.maxx = self.minx + self.fieldwidth 
        self.miny = self.simlength - self.fieldlength
        self.c = c
        self.w0 = w0
        self.a = a
        self.b = b
        self.outfluxmode = outfluxmode
        self.q0ratio = q0
        self.qshiftratio = qshift
        self.qsat = qsatinit
        self.q0 = q0*qsatinit
        self.qshift = qshift*qsatinit
        self.dt = dt
        self.lambda1 = lambda1
        self.lambda2 = lambda2
        self.lambda3 = lambda3
        self.thetar = thetar
        self.alpha = alpha
        self.delta = delta
        self.toinject = inject
        self.periodic = periodic
        self.desert = ContinuousSpace(self.simwidth + 50, self.simlength + 50, x_min = -50, y_min = -50, torus = periodic)
        self.collson = collson
        self.collrecord = []
        self.colltime = []
        self.windangle = 3*np.pi/2
        self.numsplits = 0
        if initdensity == None:
            self.initdensity = len(xs)/abs(max(ys) - min(ys))
        else:
            self.initdensity = initdensity
        self.injectparams = injectparams
        
        self.dunes = []
        self.dune_unique_ids = []
        self.schedule = SimultaneousActivation(self)

        for i in range(len(xs)):
            dune = Barchan(i, self, lws[i], rws[i])
            self.dunes.append(dune)
            self.dune_unique_ids.append(i)
            self.schedule.add(dune)
            self.desert.place_agent(dune, (xs[i], ys[i]))
        
        self.maxid = max(self.dune_unique_ids)
        self.dune_unique_ids = np.array(self.dune_unique_ids, dtype = int)
        
    def add(self, lw, rw, x, y, newid):
        if x > self.simwidth or x < 0 or y > self.simlength or y < 0 or lw <= 0.1  or rw <= 0.1:
            print('oob')
        else:
            dune = Barchan(newid, self, lw, rw)
            self.dunes =  np.append(self.dunes, dune)
            self.dune_unique_ids = np.append(self.dune_unique_ids, newid)
            self.schedule.add(dune)
            self.desert.place_agent(dune, (x, y))
            self.maxid = max(self.dune_unique_ids)
            self.lws = np.append(self.lws, lw)
            self.rws = np.append(self.rws, rw)
            self.xs = np.append(self.xs, x)
            self.ys = np.append(self.ys, y)

    
    def remove(self, dune):

        ind = np.where(self.dune_unique_ids == dune.unique_id)
        self.desert.remove_agent(dune)
        self.dunes = np.delete(self.dunes, ind)
        self.dune_unique_ids = np.delete(self.dune_unique_ids, ind)
        self.schedule.remove(dune)
        self.xs = np.delete(self.xs, ind)
        self.ys = np.delete(self.ys, ind)
        self.lws = np.delete(self.lws, ind)
        self.rws = np.delete(self.rws, ind)

        
    def injection(self, qsat):
        density = self.initdensity
        lwmin, lwmax, rwmin, rwmax = self.injectparams
        avemigdist = bss.UniformDistAverageMigrationRate(lwmin, lwmax, rwmin, rwmax, qsat, self.dt, c = self.c)
        injectnum = density * avemigdist
        
        print('inj', injectnum)
        intinjectnum = int(injectnum)
        pleftover = injectnum - intinjectnum
    
        rand = np.random.uniform()
        
        if rand <= pleftover:
            intinjectnum += 1

        for i in range(intinjectnum):
            lw = np.random.uniform(lwmin, lwmax)
            rw = np.random.uniform(rwmin, rwmax)
            x = np.random.uniform(self.minx, self.maxx)
            y = self.simlength - np.random.uniform(10, 10.01)
            newid = self.maxid + 1
            self.add(lw, rw, x, y, newid)
            
    def collisions(self, dunes, lefts, rights):
        
        alpha = self.alpha
        delta = self.delta
            
        minw = (delta/2)/(1 - alpha)
        
        pairs = cs.CollPairs(lefts, rights)
        cpairs = np.random.permutation(np.copy(pairs))
        for i in range(len(cpairs)):
            print('collision')
            c1, c2 = cpairs[i]
            d1 = dunes[c1]
            d2 = dunes[c2]
            
            if d1 in self.dunes and d2 in self.dunes:
                
                initlw1 = d1.lw
                initrw1 = d1.rw
                initlw2 = d2.lw
                initrw2 = d2.rw
                
                
                
                
                if any(np.array([initlw1, initrw1, initlw2, initrw2]) < minw):
                    print('yes', [initlw1, initrw1, initlw2, initrw2], minw)
                    if initlw1 < minw:
                        self.remove(d1)
                    elif initrw1 < minw:
                        self.remove(d1)
                    if initlw2 < minw:
                        self.remove(d2)
                    elif initrw2 < minw:
                        self.remove(d2)
                else:
                
                    initx1, inity1 = d1.pos
                    initx2, inity2 = d2.pos
    
    
                    print('inputs', [initlw1, initrw1, initlw2, initrw2, initx1, inity1, initx2, inity2])
                    
                    self.remove(d1)
                    self.remove(d2)
                    
                    l1 = lefts[c1]
                    r1 = rights[c1]
                    l2 = lefts[c2]
                    r2 = rights[c2]
                    
                    xs, ys, lws, rws = cs.OneRule(l1, l2, r1, r2, self.lambda1, 
                                                  self.lambda2, self.lambda3, self.alpha, self.delta, self.qshiftratio)
                    
                    self.collrecord.append([(initx1, inity1), (initx2, inity2), (initlw1, initrw1), (initlw2, initrw2), lws, rws, xs, ys])
                    
                    
                    for j in range(len(xs)):
                        lw = lws[j]
                        rw = rws[j]
                        x = xs[j]
                        y = ys[j]
                        newid = self.maxid + 1
                        self.add(lw, rw, x, y, newid)
                        

            
            
    def step(self, qsat, wd, fluxfield, plot):
        
        self.windangle = wind.GetAngle(wd)

        q0 = self.q0ratio*qsat
        qshift = self.qshiftratio*qsat
        self.q0 = q0
        self.qsat = qsat
        self.qshift = qshift
        dunes = self.dunes
        lws = self.lws
        rws = self.rws
        xs = self.xs
        ys = self.ys
        dt = self.dt
        
        plotlefts, plotrights, plotdunes, ff = bss.IterationCalculations(fluxfield, dunes, xs, ys, lws, rws, wd, qsat, q0, dt, w0 = self.w0 , c = self.c,
                          lambda1 = self.lambda1, lambda2 = self.lambda2, lambda3 = self.lambda3,
                          alpha = self.alpha, delta = self.delta, a = self.a, b = self.b, outfluxmode = self.outfluxmode)
        
        self.schedule.step()
        
        


        if self.collson == True:
            self.collisions(plotdunes, plotlefts, plotrights)
            

            
        if self.toinject == True:
            self.injection(qsat)
            
        if plot == True:
            return plotlefts, plotrights, ff
                    
        else:
            return None
            
        
        
        
        

        
        
        
        
        
        
