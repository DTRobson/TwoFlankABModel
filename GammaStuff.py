# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 11:10:02 2022

@author: Robbo
"""

import numpy as np

def max_gamma(wl, alpha, delta, lambda2, qshiftratio):
    return 1 + alpha/(lambda2*qshiftratio - alpha) + delta/(2*(lambda2*qshiftratio - alpha)*wl)