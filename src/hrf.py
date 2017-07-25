#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 00:03:42 2017

@author: ajoshi
"""
import scipy as sp


def pdf_gamma(x, a, b):
    # Gamma distribution evaluated at x for parameters a and b
    nx = len(x)
    pdfx = sp.zeros([nx])
    ind = x > 0
    pdfx[ind] = (b**2) * (x[ind]**(a-1)) *\
        sp.exp(-b*x[ind]) / sp.special.gamma(a)
    return pdfx


def fast_fslgamma(t):

    meanlag = 6.0
    stddev = 3.0
    a = (meanlag/stddev)**2
    b = meanlag/(stddev**2)
    t = sp.arange(len(t))*0.72 + 1e-16
    h = pdf_gamma(t, a, b)
    return h

#            % standard SPM, FSL double gamma fn
#            dg_HRF_t= ((t.^5.*exp(-t))/gamma(6) ...
#                - (1.0/6)*(t.^15.*exp(-t))/gamma(16));

def double_gamma_hrf(t):
    dg_HRF_t= ((t**5*sp.exp(-t))/sp.special.gamma(6) -
                (1.0/6)*(t**15*sp.exp(-t))/sp.special.gamma(16))
    return dg_HRF_t
    