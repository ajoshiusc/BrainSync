#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 00:03:42 2017

@author: ajoshi
"""
import scipy as sp
from hrf import fast_fslgamma
import matplotlib.pyplot as plt


h = fast_fslgamma(sp.arange(22)*0.72)

plt.plot(h)
plt.show()
#plt.close()