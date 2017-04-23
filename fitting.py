# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import time as tm
import math as ma
#import random as rd
import numpy as np
import matplotlib.pyplot as mp
import funfit as f1

xnum=10
t0=np.linspace(0,1,1000)
x=np.random.random(xnum)
y=np.zeros(xnum,)

for ind in range(0,xnum,):
    y[ind]=ma.sin(x[ind]*ma.pi*2)

mp.plot(x,y,'b*')

#==============================================================================
# Fitting with y=w0+w1*x
#==============================================================================
w0=0
w1=0
lb=.1

for iter in range(1,10000):
    yfit=f1.gs1order(w0,w1,x)
    cf1=f1.cf(y,yfit)
#    print '# Iter',str(iter),': ',str(cf1)
    
    w0=w0-lb*np.sum(yfit-y)
    w1=w1-lb*np.sum((yfit-y)*x)

t1=w0+w1*t0
mp.plot(x,y,'rx')
mp.plot(t0,t1)

#==============================================================================
# Fitting with y=w0+w1*x+w2*x^2
#==============================================================================
w0=0
w1=0
w2=0
lb=.1

for iter in range(1,5000):
    yfit=f1.gs2order(w0,w1,w2,x)
    cf2=f1.cf(y,yfit)
#    print '# Iter',str(iter),': ',str(cf2)
    
    w0=w0-lb*np.sum(yfit-y)
    w1=w1-lb*np.sum((yfit-y)*x)
    w2=w2-lb*np.sum((yfit-y)*x**2)
    
t2=w0+w1*t0+w2*t0**2
mp.plot(x,y,'rx')
mp.plot(t0,t2)

#==============================================================================
# Fitting with 5 degree polynomial
#==============================================================================
w0=0
w1=0
w2=0
w3=0
w4=0
w5=0
lb=.1

for iter in range(1,10000):
    yfit=f1.gs5order(w0,w1,w2,w3,w4,w5,x)
    cf5=f1.cf(y,yfit)
#    print '# Iter',str(iter),': ',str(cf5)
    
    w0=w0-lb*np.sum(yfit-y)
    w1=w1-lb*np.sum((yfit-y)*x)
    w2=w2-lb*np.sum((yfit-y)*x**2)
    w3=w3-lb*np.sum((yfit-y)*x**3)
    w4=w4-lb*np.sum((yfit-y)*x**4)
    w5=w5-lb*np.sum((yfit-y)*x**5)
    
t5=f1.gs5order(w0,w1,w2,w3,w4,w5,t0)
mp.plot(x,y,'rx')
mp.plot(t0,t5)

#==============================================================================
# Fitting with n degree polynomial
#==============================================================================
n=5
w=np.zeros([1,n+1])
lb=.05
xext=f1.gsnxext(x,n)
cfprev=f1.cf(y,f1.gsnorder(w.T,n,x))
gdsucc=1

start=tm.clock()
for iter in range(1,90000):
    yfit=f1.gsnorder(w.T,n,x)
    cf=f1.cf(y,yfit)
    
    if (cf-cfprev)<1e-2:
#        print '# Iter',str(iter),': ',str(cf)
        w=w-lb*np.dot((yfit-y),xext)
        cfprev=cf
    else:
        print 'Parameter is too large! Please choose a smaller one.'
        gdsucc=0
        break
    
end=tm.clock()

if gdsucc==1:        
    t1=np.dot(f1.gsnxext(t0,n),w.T)
    mp.plot(x,y,'ro')
    mp.plot(t0,t1)
    leg='M='+str(n)
    mp.legend(('Xn',leg))
    mp.title('Fitting with polynomial')
    print "Successful! Running time: %f s" % (end-start)
    print "> Cost Value: %f" % cf
    print "> Fitting with %d degree polynomial" % n
    print "> Parameter of the polynomial:",w
    print "> Absolute value of (y-yfit):",abs(y-yfit)
    
#==============================================================================
# Regularized fitting with n degree polynomial
#==============================================================================
n=12
w=np.zeros([1,n+1])
lb=.01
regpar=1
xext=f1.gsnxext(x,n)
cfprev=f1.cf(y,f1.gsnorder(w.T,n,x))
gdsucc=1

start=tm.clock()
for iter in range(1,90000):
    yfit=f1.gsnorder(w.T,n,x)
    cf=f1.cf(y,yfit)
    
    if (cf-cfprev)<1e-4:
#        print '# Iter',str(iter),': ',str(cf)
        w=w-lb*(np.dot((yfit-y),xext)+regpar*w)
        w[0]=w[0]-lb*(np.dot((yfit-y),xext)[0])
        cfprev=cf
    else:
        print 'Parameter is too large! Please choose a smaller one.'
        gdsucc=0
        break
    
end=tm.clock()

if gdsucc==1:        
    t1=np.dot(f1.gsnxext(t0,n),w.T)
    mp.plot(x,y,'ro')
    mp.plot(t0,t1)
    leg='M='+str(n)
    mp.legend(('Xn',leg))
    mp.title('Regularized fitting with polynomial (lambda=%0.4f)'%lb)
    print "Successful! Running time: %f s" % (end-start)
    print "> Cost Value: %f" % cf
    print "> Fitting with %d degree polynomial" % n
    print "> Parameter of the polynomial:",w
    print "> Absolute value of (y-yfit):",abs(y-yfit)
    