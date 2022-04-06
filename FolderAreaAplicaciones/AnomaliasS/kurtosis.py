# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 10:47:30 2022

@author: jose.contreras
"""


from scipy.stats import kurtosis
x = [55, 78, 65, 98, 97, 60, 67, 65, 83, 65, 5, 784, 874875]
N = len(x)

def Kurtosis (x):
    xprom=0
    for i in range (0, N):
        xprom = xprom + x[i]
        
    xprom = xprom/N
    
    m4 = 0
    m2 = 0
    for i in range (0, N):
        m = (1/N)*(x[i]-xprom)**4
        m4 = m4 + m
        m1 = (1/N)*(x[i]-xprom)**2
        m2 = m2 + m1
    
        
    K = m4/(m2**2)
    return K

valor = Kurtosis(x)
print('Kurtosis = ', valor)

print(kurtosis(x, fisher=False))
