# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 11:12:21 2022

@author: jose.contreras
"""

from scipy.stats import skew
x = [55, 78, 65, 98, 97, 60, 67, 65, 83, 65, 5, 784, 874875]
N = len(x)

def Skew (x):
    xprom=0
    for i in range (0, N):
        xprom = xprom + x[i]
        
    xprom = xprom/N
    
    m3 = 0
    m2 = 0
    for i in range (0, N):
        m = (1/N)*(x[i]-xprom)**3
        m3 = m3 + m
        m1 = (1/N)*(x[i]-xprom)**2
        m2 = m2 + m1
    
        
    G = (((N*(N-1))**0.5)/(N-2))*(m3/(m2**(1.5)))
    return G

valor = Skew(x)
print('Skew = ', valor)

print(skew(x, bias=False))
