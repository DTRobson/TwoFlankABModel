# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 11:10:02 2022

@author: Robbo
"""

import numpy as np

def max_gamma(wl, alpha, delta, lambda2, qshiftratio, lambda1):
    if 1 + alpha/(lambda2*qshiftratio - alpha) + delta/(2*(lambda2*qshiftratio - alpha)*wl) > 0:
        return min(1 + alpha/(lambda2*qshiftratio - alpha) + delta/(2*(lambda2*qshiftratio - alpha)*wl), 2*lambda2/lambda1 - 1)
    else:
        return 2*lambda2/lambda1 - 1
    
    
    