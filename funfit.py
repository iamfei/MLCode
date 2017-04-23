#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 23:36:24 2017

@author: iamfei
"""
import numpy as np

def gs1order(w0,w1,x):
    return w0+w1*x

def gs2order(w0,w1,w2,x):
    return w0+w1*x+w2*x**2

def gs5order(w0,w1,w2,w3,w4,w5,x):
    return w0+w1*x+w2*x**2+w3*x**3+w4*x**4+w5*x**5

def gsnorder(w,n,x):
    norder=np.dot(gsnxext(x,n),w)
    return norder.T

def gsnxext(x,n):
    num=len(x)
    xext=np.ones([num,n+1])
    for col in range(1,n+1):
        xext[:,col]=x**col
    return xext
    
def cf(x,y):
    return np.sum((x-y)**2)/2
