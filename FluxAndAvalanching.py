# -*- coding: utf-8 -*-
"""
Created on Tue May 31 13:00:39 2022

@author: Robbo
"""

from shapely.geometry import Polygon, MultiPolygon, box, MultiPoint, Point
import numpy as np
from shapely.ops import unary_union
import WinddirectionStuff as wind
import matplotlib.pyplot as plt
import multiprocessing as mp

BOUNDS = box(1,1,2,2)

def EquivWidth(v, lambda2, lambda3, flankorbody):
    '''Function calculates the equivalent width of either a flank or full barchan with volume v'''
    if v > 0:
        if flankorbody == 'flank':
            w = (6*v/(lambda2*lambda3))**(1/3)
        if flankorbody == 'body':
            w = 2*(3*v/(lambda2*lambda3))**(1/3)
        return w
    else:
        return 0
        

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
    

class simplemodel():
        def __init__(self, simw, siml):
            self.simwidth = simw
            self.simlength = siml
            

class simpledune():
    def __init__(self, simw, siml):
        self.leftinflux = 0
        self.rightinflux = 0
        self.rightoutflux = 0
        self.leftoutflux = 0
        self.lw = 0
        self.rw = 0
        self.model = simplemodel(simw, siml)

        

def Rotate(polygon, windirection):
    '''
    
    Rotate a polygon clockwise or anticlockwise to match the new winddireciton

    Parameters
    ----------
    polygon : shapely.geometry.Polygon
        The polygon to be rotated.
    windirection : tuple
        The new winddirection.
    clock : bool, optional
        Should the rotation be clockwise? The default is False.

    Returns
    -------
    rotated : shapely.geometry.Polygon
        The rotated Polygon.

    '''
    
    if len(windirection) == 2:
            windirection = np.append(windirection, 0)
        
    crossproduct = np.cross(windirection, (0, -1, 0))
        
    theta = np.arcsin(crossproduct)[-1]
        
    if np.sign(windirection[1]) == 1:
            theta = np.pi - theta
        
        
    
    mat = np.array([[np.cos(theta), - np.sin(theta)],
                        [np.sin(theta), np.cos(theta)]], dtype = np.float32)
    
    
    if type(polygon) == Polygon:
    

        coords = np.array(polygon.exterior.xy, dtype = np.float32)
        rotatedcoords = mat.dot(coords)
        rotated = Polygon([(rotatedcoords[0, i], rotatedcoords[1, i])
                           for i in range(len(rotatedcoords[0, :]))])
        

    elif type(polygon) == MultiPolygon:
        rotpols = []
        for pol in polygon.geoms:
            coords = np.array(pol.exterior.xy)
            rotatedcoords = mat.dot(coords)
            rot = Polygon([(rotatedcoords[0, i], rotatedcoords[1, i])
                           for i in range(len(rotatedcoords[0, :]))])
            rotpols.append(rot)
        rotated = unary_union(rotpols)
        
    elif type(polygon) == Point:
        rotated = polygon
        
    else:        
        rotated = polygon

    return rotated
                      
def RotateXY(x, y, wd):
    
    if len(wd) == 2:
            windirection = np.append(wd, 0)
        
    crossproduct = np.cross(windirection, (0, -1, 0))
        
    theta = np.arcsin(crossproduct)[-1]
        
    if np.sign(windirection[1]) == 1:
            theta = np.pi - theta
        
        
    
    mat = np.array([[np.cos(theta), - np.sin(theta)],
                        [np.sin(theta), np.cos(theta)]], dtype = np.float32)
    newx, newy = mat.dot(np.array([x, y]))
    return newx, newy
    

def Barchan(x, y, wl, wr, hornlwratio = 1.8,
            lwratio = 1., alpha = 0.05, delta = 4.6, lambda3 = 1/3):
    '''
    Function creates a shapely shapely.geometry.Polygon object
    for the left and right flanks of a barchan dune allowing for
    the total barchan to be asymmetric assuming the wind
    is in the negative y-direction.

    Parameters
    ----------
    x : float
        The x-coordinate of the toe of the dune.
    y : float
        The y-coordinate of the toe of the dune.
    wl : float
        The width of the left flank of the dune.
    wr : float
        The width of the right flank of the dune.
    hornlwratio : float, optional
        The length of a horn is given by the relevant flank
        width multiplied by this value. The default is 1.8.
    lwratio : float, optional
        The length of the body of a symmetric barchan is given
        by the total width multiplied by this value. The default is 1..
    alpha : float, optional
        See decsription of delta. The default is 0.05.
    delta : float, optional
        The width of a horn is given by delta/2 + alpha*flankwidth.
        The default is 4.6.

    Returns
    -------
    leftflank : shapely.geometry.Polygon
        Polygon of the left flank of the dune.
    rightflank : shapely.geometry.Polygon
        Polygon of the right flank of the dune.

    '''
    
    
    x = np.mean(x)
    y = np.mean(y)
    wl = np.mean(wl)
    wr = np.mean(wr)
    x = np.nan_to_num(x)
    y = np.nan_to_num(y)
    wl = np.nan_to_num(wl)
    wr = np.nan_to_num(wr)
    if type(x) == type(np.nan) or type(y) == type(np.nan) or type(wl) == type(np.nan) or type(wr) == type(np.nan):

        return Point((0,0)), Point((0,0))
    else:
        minw = (delta/2)/(1 - alpha)
        if wl <=0 or wr <= 0:
            return box(x - 1e-5,y, x + 1e-5, y-2e-5), box(x - 1e-5,y, x + 1e-5, y-2e-5)
          
        else:
        
            leftlength = lwratio * wl 
            lefthornwidth = alpha*wl + delta/2
            lefthornlength = wl*hornlwratio
            
            rightlength = lwratio * wr
            righthornwidth = alpha*wr + delta/2
            righthornlength = wr*hornlwratio
            
            #we assume that the length of the body of the asymmetric barchan
            #is equal to the average of the lengths of the bodies of symmetric
            #barchans with widths wl and wr
            bodylength = (leftlength + rightlength)/2
            
            pl1 = (x, y) #toe
            pl2 = (x - wl, y - bodylength) #leftmost extent
            pl3 = (x - wl + lefthornwidth/2, y - lefthornlength)# tip ofleft horn
            pl4 = (x, y - bodylength) #brink
            
            if wl >= minw:
                leftflank = Polygon([pl1, pl2, pl3, pl4])
            else:

                leftflank  = box(x-wl, y - lefthornlength, x, y)
            
            pr1 = (x, y)
            pr2 = (x + wr, y - bodylength) #rightmost extent
            pr3 = (x + wr - righthornwidth/2, y - righthornlength) # tip of right horn
            pr4 = (x, y -  bodylength)
            
        
            if wr >= minw:
                rightflank = Polygon([pr1, pr2, pr3, pr4])
            else:
                rightflank = box(x, y-righthornlength, x+wr, y)
            
            if leftflank.is_valid == True and rightflank.is_valid == True:
                return leftflank, rightflank
            else:
                vl = VS(wl, wr, lwratio, hornlwratio, lambda3, True)
                vr = VS(wr, wl, lwratio, hornlwratio, lambda3, True)
                v = vl + vr
                w = EquivWidth(v, hornlwratio, lambda3, 'body')
                leftflank = box(x-w/2, y-w, x, y)
                rightflank = box(x, y-w, x+w/2, y)
                return  leftflank, rightflank


def Shadow(fluxfield, leftflank, rightflank, windirection = [0, -1], shift = 1e-4):
    '''
    Subtract the left and right flanks and the shadow they create from the flux field

    Parameters
    ----------
    fluxfield : shapely.geometry.Polygon
        The current flux field.
    leftflank : shapely.geometry.Polygon
        The left flank of a barchan.
    rightflank : shapely.geometry.Polygon
        The right flank of a barchan.
    windirection : tuple, optional
        The wind direction. The default is [0, -1].
    shift : float, optional
        The distance to shift the barchan down so that it
        does not occupy its own shadow. The default is 0.1.

    Returns
    -------
    newfluxfield : shapely.geometry.Polygon
        The new flux field with the shadow of the relevant barchan removed.

    '''
    
    if type(leftflank) == Polygon and type(rightflank) == Polygon:
    
        barchan = unary_union([leftflank, rightflank])
        
        if type(barchan) == Polygon:
            xb, yb = barchan.exterior.xy
        else:
            print('babby')
            xb = np.empty(0, dtype = float)
            yb = np.empty(0, dtype = float)
            for geom in barchan.geoms:
                xg, yg = geom.exterior.xy
                xb = np.append(xb, xg)
                yb = np.append(yb, yg)

        
        if type(fluxfield) == Polygon:
            fieldx, fieldy = fluxfield.exterior.xy
        elif type(fluxfield) == MultiPolygon:
            fieldx = np.empty(0, dtype = object)
            fieldy = np.empty(0, dtype = object)
            for geom in fluxfield.geoms:
                x, y = geom.exterior.xy
                fieldx = np.append(fieldx, np.array(x))
                fieldy = np.append(fieldy, np.array(y))
            fieldx = np.ndarray.flatten(fieldx)
            fieldy = np.ndarray.flatten(fieldy)


        x = np.array(xb, dtype = np.float32)
        y = np.array(yb, dtype = np.float32)
        
        shiftedx = x + shift * windirection[0] * np.ones(len(x))
        shiftedy = y + shift * windirection[1] * np.ones(len(y))
        

        
        points = [(shiftedx[i], shiftedy[i]) for i in range(len(shiftedx))]
        
        points.append((np.min(shiftedx), np.min(fieldy)))
        points.append((np.max(shiftedx), np.min(fieldy)))
        
        mp = MultiPoint(points)
    
        
        shadow = mp.convex_hull
        
        newfluxfield = fluxfield.difference(shadow)

        
        return newfluxfield
        
    else:

        return fluxfield
    
    
        

def OverlapWidth(poly1, poly2):
    if poly1.intersects(poly2) == True:
        overlap = poly1.intersection(poly2)
        bs = overlap.bounds
        if len(bs) > 0:
            minx, miny, maxx, maxy = bs
            width = abs(maxx - minx)
            return width
        else:
            return 0
    else:
        return 0




def IndividualAmbientFlux(fluxfield, leftflank, rightflank, q0, dt, winddirection = [0, -1]):
    
    barchan = unary_union([leftflank, rightflank])
    leftinflux = 0
    rightinflux = 0
    if fluxfield.intersects(barchan):
        fluxfield = Shadow(fluxfield, leftflank, rightflank)
        leftoverlappingwidth = OverlapWidth(fluxfield, leftflank)
        leftinflux = leftoverlappingwidth * q0 * dt
        rightoverlappingwidth = OverlapWidth(fluxfield, rightflank)
        rightinflux = rightoverlappingwidth * q0 * dt

    


    return fluxfield, leftinflux, rightinflux

def PolyLengthWidth(poly):
    
    bs = poly.bounds
    if len(bs) > 0:
        minx, miny, maxx, maxy = bs
        length = abs(maxy - miny)
        width = abs(maxx - minx)
        return length, width
    else:
        return 0, 0
    
def MigrationRate(lw, rw, qsat, dt, w0 = 0, c = 1.):
    return c*qsat/(lw + rw + w0)*dt 

def UniformDistAverageMigrationRate(lwmin, lwmax, rwmin, rwmax, qsat, dt, w0 = 0, c = 1.):
    a = lwmin + rwmin
    b = lwmax + rwmax
    return c*qsat*np.log((b + w0)/(a + w0))/(b - a)*dt


def LateralShift(lw, rw, x, lambda1, lambda2, lambda3, qshift, dt, windangle):
    lvold = VS(lw, rw, lambda1, lambda2, lambda3, True)
    rvold = VS(rw, lw, lambda1, lambda2, lambda3, True)
    lx = x - lw
    rx = x + rw
    if lw != rw:
        bigw = max(lw, rw)
        littlew = min(lw, rw)
        dv = qshift*dt*lambda2*(bigw - littlew)*abs(np.sin(windangle))

    else:
        dv = 0.
        
    lvnew = lvold + np.sign(rw - lw) * dv
    rvnew = rvold - np.sign(rw - lw) * dv
    
    lwnew = EquivWidth(lvnew, lambda2, lambda3, 'flank')
    rwnew = EquivWidth(rvnew, lambda2, lambda3, 'flank')
    
    
    if lw < rw:
        xnew = lx + lwnew
    else:
        xnew = rx - rwnew


    return lwnew, rwnew, xnew


def HeavisideIsh(x, extents, posorneg):
    if posorneg == 'pos':
        outputs = np.where(x - extents >= 0, 1., 0.)
    elif posorneg == 'neg':
        outputs = np.where(extents - x >= 0, 1., 0.)
    return outputs
        

def HornPoly(x, y, hornw, miny):
    x = np.nan_to_num(x)
    y = np.nan_to_num(y)
    hornw = np.nan_to_num(hornw)
    return box(x-hornw/2, miny, x + hornw/2, y-1e-3)

def HornFlux(hpoly, dunes, leftflanks, rightflanks, qsat, dt):
    
    for i in range(len(dunes)):
        d = dunes[i]
        leftflank = leftflanks[i]
        rightflank = rightflanks[i]
        hpoly, lin, rin = IndividualAmbientFlux(hpoly, leftflank, rightflank, qsat, dt)
        d.leftinflux += lin
        d.rightinflux += rin
        d = None
        del d
    dunes = None
    del dunes
    return hpoly
    
        


def GenDunesEtc(dunes, xs, ys, lws, rws, wd, qsat, q0, dt, w0 = 0., c = 1.,
                          lambda1 = 1., lambda2 = 1.8, lambda3 = 1/6, alpha = 0.05, delta = 4.6):
    lefts = np.empty(0, dtype = object)
    rights = np.empty(0, dtype = object)
    rotys = np.empty(0, dtype = float)
    rotlefts = np.empty(0, dtype = object)
    rotrights = np.empty(0, dtype = object)
    rotlhxs = np.empty(0, dtype = float)
    rotlhys = np.empty(0, dtype = float)
    rotrhxs = np.empty(0, dtype = float)
    rotrhys = np.empty(0, dtype = float)
    lhws = np.empty(0, dtype = float)
    rhws = np.empty(0, dtype = float)

    
    for i in range(len(dunes)):
        x = xs[i]
        y = ys[i]
        lw = lws[i]
        rw = rws[i]
        left, right = Barchan(x, y, lw, rw, lambda2, lambda1, alpha, delta, lambda3)
        lefts = np.append(lefts, left)
        rights = np.append(rights, right)
        rotlefts= np.append(rotlefts, Rotate(left, wd))
        rotrights = np.append(rotrights, Rotate(right, wd))
        
        lhw = (alpha*lw + delta/2)
        lhx = x - lw + lhw/2
        lhy = y - lambda2*lw
        lhws = np.append(lhws, lhw)
        rotlhx, rotlhy = RotateXY(lhx, lhy, wd)
        rotlhxs = np.append(rotlhxs, rotlhx)
        rotlhys=np.append(rotlhys, rotlhy)
        
        rhw = (alpha*rw + delta/2)
        rhx = x + rw - rhw/2
        rhy = y - lambda2*rw
        rhws= np.append(rhws, rhw)
        rotrhx, rotrhy = RotateXY(rhx, rhy, wd)
        rotrhxs= np.append(rotrhxs, rotrhx)
        rotrhys = np.append(rotrhys, rotrhy)
        
        rotx, roty = RotateXY(x, y, wd)
        rotys = np.append(rotys, roty)
    dunes = None
    del dunes
    return lefts, rights, rotys, rotlefts, rotrights, rotlhxs, rotlhys, rotrhxs, rotrhys, lhws, rhws
        


    

def IterationCalculations(fluxfield, dunes, xs, ys, lws, rws, wd, qsat, q0, dt, w0 = 0., c = 1.,
                          lambda1 = 1., lambda2 = 1.8, lambda3 = 1/6, alpha = 0.05,
                          delta = 4.6, a = 0.45, b = 0.1, outfluxmode = 'Duran',
                          plotting = False):
    
    lvols = VS(lws, rws, lambda1, lambda2, lambda3, True)
    rvols = VS(rws, lws, lambda1, lambda2, lambda3, True)
    
    dunes = np.array(dunes, dtype = object)
    num = len(xs)
    
    rotatedflux = Rotate(fluxfield, wd)

    blah1, miny, blah2, blah3 = rotatedflux.bounds
    

    lefts, rights, rotys, rotlefts, rotrights, rotlhxs, rotlhys, rotrhxs, rotrhys, lhws, rhws = GenDunesEtc(dunes, xs, ys, lws, rws, wd,
                                                                                                            qsat, q0, dt,
                                                                                                            w0 = w0, c = c,
                          lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, alpha = alpha, delta = delta)

    orderedinds = [i for _, i in sorted(zip(-rotys, np.arange(num)))]
    
    lefthorns = []
    righthorns = []
    

    for ind in orderedinds:
        
        
        dune = dunes[ind]
        left = rotlefts[ind]
        right = rotrights[ind]
        lhx = rotlhxs[ind]
        lhy = rotlhys[ind]
        lhw = lhws[ind]
        rhx = rotrhxs[ind]
        rhy = rotrhys[ind]
        rhw = rhws[ind]    
        lhx = np.nan_to_num(lhx)
        lhy = np.nan_to_num(lhy)
        rhx = np.nan_to_num(rhx)
        rhy = np.nan_to_num(rhy)
        lv = lvols[ind]
        rv = rvols[ind]
        


        rotatedflux, leftinflux, rightinflux = IndividualAmbientFlux(rotatedflux, left, right, q0, dt)
        
        dune.leftinflux += leftinflux
        dune.rightinflux += rightinflux
        
        lefthorn = HornPoly(lhx, lhy, lhw, miny)
        
        righthorn = HornPoly(rhx, rhy, rhw, miny)
        


        if outfluxmode == 'Duran':
    
            qoutleft = (a*leftinflux + b * qsat) * lws[ind]/lhw
            qoutright = (a*rightinflux + b * qsat) * rws[ind]/rhw
            

        
        else:
            qoutleft = qsat
            qoutright = qsat
        
        
        qoutleft = min(qoutleft, (lv + leftinflux)/(lhw*dt))
        qoutright = min(qoutright, (rv + rightinflux)/(rhw*dt))
        
        
        dune.leftoutflux += qoutleft * lhw * dt
        dune.rightoutflux += qoutright*dt*rhw
        
        
        otherinds = np.array(orderedinds)[np.where(orderedinds != ind)]
        
        lefthorn = HornFlux(lefthorn, dunes[otherinds], rotlefts[otherinds], rotrights[otherinds], qoutleft, dt)
        righthorn = HornFlux(righthorn, dunes[otherinds], rotlefts[otherinds], rotrights[otherinds], qoutright, dt)
        
        

        if plotting == True:
            lefthorns.append(Rotate(lefthorn, (-wd[0], wd[1])))
            righthorns.append(Rotate(righthorn, (-wd[0], wd[1])))
            
        
        migrate = MigrationRate(lws[ind], rws[ind], qsat, dt, w0, c)

        
        dune.tomove = (migrate*wd[0], migrate*wd[1])
        dune = None
        del dune

        
    return xs, ys, lws, rws, lefts, rights, dunes, Rotate(rotatedflux, (-wd[0], wd[1])), lefthorns, righthorns
        
        
        



        
        
        
        
    
    
    
        
        

        
        
        
        
    
    
    
    
    
    
    
    
    
    
        


    
    
    
    
    
    
    
    
    

    
    

