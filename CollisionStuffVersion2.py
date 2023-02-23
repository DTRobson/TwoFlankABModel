# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 13:27:45 2022

@author: Robbo
"""

from shapely.geometry import Polygon, MultiPolygon, box, MultiPoint, Point
import matplotlib.pyplot as plt
import numpy as np
from shapely.ops import unary_union
import multiprocessing as mp
from GammaStuff import max_gamma as gamma_c
from FluxAndAvalanching import Barchan as barc


def XYFROMFLANK(flank_poly, leftorright):
    xs = np.array(flank_poly.exterior.xy[0],dtype = float)
    ys = np.array(flank_poly.exterior.xy[1],dtype = float)
    y = max(ys)
    if leftorright == 'left':
        x = max(xs)
    if leftorright == 'right':
        x = min(xs)
    return x, y
    
    
def EquivWidth(v, lambda2, lambda3, flankorbody):
    '''Function calculates the equivalent width of either a flank or full barchan with volume v'''
    if flankorbody == 'flank':
        w = (6*v/(lambda2*lambda3))**(1/3)
    if flankorbody == 'body':
        w = 2*(3*v/(lambda2*lambda3))**(1/3)
    return w
    

def VS(w, otherw, lambda1, lambda2, lambda3):
    '''Function calculates the volume of the body, horn and their total of a flank with width w
    where the other flank has width otherw'''
    vb = 1/6 * lambda1 * lambda3 * w**2 * (w + otherw)/2
    vh = 1/6 * lambda3 * w**2 * (lambda2 * w - lambda1*(w + otherw)/2)
    return vb, vh, vb + vh



def COMS(x, y, wl, wr, lambda1, lambda2, lambda3, alpha = 0.05, delta = 4.6, oneonly = None):
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
    if oneonly == None:
        return (x, y), (xl, yl), (xr, yr)
    elif oneonly == 'x':
        return x
    elif oneonly == 'y':
        return y


def PolyWidth(poly):
    minx, miny, maxx, maxy = poly.bounds
    w = abs(maxx - minx)
    return w

def CentroidIntersect(left1, right1, left2, right2):
    x1, y1 = XYFROMFLANK(left1, 'left')
    x2, y2 = XYFROMFLANK(left2, 'left')
    w1 = PolyWidth(left1) + PolyWidth(right1)
    w2 = PolyWidth(left2) + PolyWidth(right2)
    
    if w1 <= w2:
        cl = left1.centroid.buffer(0.01)
        cr = right1.centroid.buffer(0.01)
        cb = unary_union([left1, right1]).centroid.buffer(0.01)
        target = unary_union([left2, right2])
    else:
        cl = left2.centroid.buffer(0.01)
        cr = right2.centroid.buffer(0.01)
        cb = unary_union([left2, right2]).centroid.buffer(0.01)
        target = unary_union([left1, right1])
    
    collide = False
    
    if cl.intersects(target) == True or cr.intersects(target) == True or cb.intersects(target) == True:
        collide = True
    
    return collide
        
    


def CollPairs(lefts, rights):
    n = len(lefts)
    info = [(ind1, ind2) for ind1 in range(n) for ind2 in range(ind1 + 1, n) if CentroidIntersect(lefts[ind1], rights[ind1], lefts[ind2], rights[ind2]) == True]
    return info
    
def OneRule(left1, left2, right1, right2, lambda1, lambda2, lambda3, alpha, delta, qshiftratio):
    
   
    colliders = []
    cvols = []
    colxs = []
    colys = []
    others = []
    ovols = []
    ows = []
    oxs = []
    oys = []
    l1col = False
    r1col = False
    l2col = False
    r2col = False
        
    l1w = PolyWidth(left1)
    l2w = PolyWidth(left2)
    r1w = PolyWidth(right1)
    r2w = PolyWidth(right2)

    barc1 = unary_union([left1, right1])
    barc2 = unary_union([left2, right2])

    comx1 = barc1.centroid.x
    comy1 = barc1.centroid.y
    comxl1 = left1.centroid.x
    comyl1 = left1.centroid.y
    comxr1 = right1.centroid.x
    comyr1 = right1.centroid.y
    
    comx2 = barc2.centroid.x
    comy2 = barc2.centroid.y
    comxl2 = left2.centroid.x
    comyl2 = left2.centroid.y
    comxr2 = right2.centroid.x
    comyr2 = right2.centroid.y
    
    
    blah5, blah6, l1v = VS(l1w, r1w, lambda1, lambda2, lambda3)
    blah7, blah8, r1v = VS(r1w, l1w, lambda1, lambda2, lambda3)
    blah9, blah10, l2v = VS(l2w, r2w, lambda1, lambda2, lambda3)
    blah11, blah12, r2v = VS(r2w, l2w, lambda1, lambda2, lambda3)
    
    originalvtot = (l1v + r1v + l2v + r2v)
    
    totcomx = (comx1*(l1v + r1v) + comx2*(l2v + r2v))/originalvtot
    totcomy = (comy1*(l1v + r1v) + comy2*(l2v + r2v))/originalvtot
    

    if left1.intersects(left2):
        l1col = True
        l2col = True
    if left1.intersects(right2):

        l1col = True
        r2col = True
    if right1.intersects(left2):

        r1col =  True
        l2col = True
    if right1.intersects(right2):

        r2col = True
        r1col = True
    
    if l1col == True:

        colliders.append(left1)
        cvols.append(l1v)
        colxs.append(comxl1)
        colys.append(comyl1)
        
    else:
        others.append(left1)
        ovols.append(l1v)
        ows.append(l1w)
        oxs.append(comxl1)
        oys.append(comyl1)
    if l2col == True:

        colliders.append(left2)
        cvols.append(l2v)
        colxs.append(comxl2)
        colys.append(comyl2)
    else:
        others.append(left2)
        ovols.append(l2v)
        ows.append(l2w)
        oxs.append(comxl2)
        oys.append(comyl2)
    if r1col == True:

        colxs.append(comxr1)
        colys.append(comyr1)
        colliders.append(right2)
        cvols.append(r1v)
    else:
        others.append(right2)
        ovols.append(r1v)
        ows.append(r1w)
        oxs.append(comxr1)
        oys.append(comyr1)
    if r2col == True:

        colliders.append(right2)
        cvols.append(r2v)
        colxs.append(comxr2)
        colys.append(comyr2)
    else:
        others.append(right2)
        ovols.append(r2v)
        ows.append(r2w)
        oxs.append(comxr2)
        oys.append(comyr2)
        

    
    cvoltot = sum(cvols)   
    cwtot = EquivWidth(cvoltot, lambda2, lambda3, 'body')  
    
    toosmall = []
    tsvols = []
    tsxs = []
    tsys = []

    for i in range(len(ows)):
        v = ovols[i]
        w = EquivWidth(v, lambda2, lambda3, 'flank')
        p = others[i]
        x = oxs[i]
        y = oys[i]

        ratio = w/cwtot
        ratioinv = cwtot/w
        maxasymmratio = gamma_c(w, alpha, delta, lambda2, qshiftratio, lambda1)
        if  ratioinv < maxasymmratio and ratio < maxasymmratio:
            colliders.append(p)
            cvols.append(v)
            colxs.append(x)
            colys.append(y)
        else:

            toosmall.append(p)
            tsvols.append(v)
            tsxs.append(x)
            tsys.append(y)
            
    
    cvols = np.array(cvols)
            

    cvoltot2 = sum(cvols)

    colcomx = sum(np.array(colxs)*np.array(cvols))/cvoltot2
    colcomy = sum(np.array(colys)*np.array(cvols))/cvoltot2
    
    leftfrac = sum(cvols[np.where(colxs <= colcomx)])/cvoltot2
    rightfrac = sum(cvols[np.where(colxs > colcomx)])/cvoltot2


    collw = EquivWidth(cvoltot2*leftfrac, lambda2, lambda3, 'flank')
    colrw = EquivWidth(cvoltot2*rightfrac, lambda2, lambda3, 'flank')
    colw = (collw+colrw)/2
    coly = colcomy + lambda1*colw/2 - (lambda1**2 - (lambda1-lambda1)**2)*colw/(8*lambda1)
    tsvols = np.array(tsvols)
    tsws = EquivWidth(tsvols, lambda2, lambda3, 'body')
    tsws = np.array(tsws)/2
    txs = []
    tys = []
    

    for i in range(len(tsvols)):
        x = tsxs[i]
        thew = max(colrw, collw)
        y = coly - lambda2 * thew + np.random.normal(0, 0.00001)
        txs.append(x)
        tys.append(y)

    
    xs = txs
    ys = tys
    lws = tsws
    rws = tsws

    xs = np.append(xs,colcomx)
    ys = np.append(ys, coly)
    lws = np.append(lws, collw)
    rws = np.append(rws, colrw)
    tsvols = np.append(tsvols, cvoltot2)
    
    bs = [unary_union(barc(xs[i], ys[i], lws[i], rws[i], lambda2,
                                                    lambda1, alpha, delta, lambda3)) for i in range(len(xs))]


    tempcomx = sum(tsvols/originalvtot * np.array([b.centroid.x for b in bs]))
    tempcomy = sum(tsvols/originalvtot * np.array([b.centroid.y for b in bs]))
    
    bs = []

    diffx = totcomx - tempcomx
    diffy = totcomy - tempcomy
    
    xs += diffx
    ys += diffy
    

    print('Check conservation', originalvtot - sum(tsvols))
    return xs, ys, lws, rws
    
    
    
    
        
        
        
        
    
    
   
    
        
            
    
            
    
    
    
    

    


    



