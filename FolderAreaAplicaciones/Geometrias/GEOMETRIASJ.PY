import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
from   statsmodels.nonparametric.kernel_regression import KernelReg
from   scipy.fft import fft, fftfreq
from   scipy     import pi
import scipy.signal
import openpyxl
from   scipy.fftpack import fft
#===========================================================================
from   tkinter     import *
from   tkinter.ttk import *
from   pylab       import *
import math
from   statistics  import mean
from   six.moves   import tkinter as tk
import time 

PI = 3.1416
TYPE=1
n=140 #numero de puntos
#Medidas de geometría del balancín
A=55
C=48
I=48
P=57
H=96
G=41.5
R=18
CW=-1 #horario positivo antihorario negativo
#===========CALCULOS DE ANGULOS TEORICOS===========================#
K =(I**2 + (H-G)**2)**0.5
print(K)
#Calculo del angulo que da la mínima posición
if TYPE == 1: #Convencional
    PHIt = math.asin(I/K)
    THETAktopt = math.acos((K**2+(R+P)**2-C**2)/(2*K*(R+P)))
    if THETAktopt > PHIt:
        THETAtopt = PHIt + THETAktopt
    else:
        THETAtopt = PHIt - THETAktopt
elif TYPE == 2: #MARK II
    PHIt = math.atan(I/(H-G)) + PI
    THETAktopt=math.acos((K**2+(P+R)**2-C**2)/(2*K*(P-R)))
    THETAtopt=THETAktopt-(PHIt-PI)
else:  #Airbalance
    PHIt=PI - math.atan(I/(H-G))
    THETAktopt=math.acos((K**2+(P-R)**2-C**2)/(2*K*(P-R))) 

print(PHIt)
#crear el vector de tamaño n desde 0+THETAtopt hasta 2*PI+THETAtopt

THETAtt=np.arange(0, 2*PI, 2*PI/(n))
THETAtt1=abs(2*PI-THETAtt)


#====================================Calculos de ángulos=====================================================#
#inicio los vectores 
BETAt  = np.zeros(n)
PSIt   = np.zeros(n)
PSIbt  = np.zeros(n)
PSItt  = np.zeros(n)
ALPHAt = np.zeros(n)
Sprt   = np.zeros(n)
Jt     = np.zeros(n)
THETAk = np.zeros(n)

for i in range (0, n-1):
    
    if TYPE==1: #convencional
        THETAk[i]=THETAtt[i]+(CW*PHIt)
        PSIbt=math.acos((C**2 + K**2 - (P+R)**2)/(2*C*K))
        PSItt=math.acos((C**2 + K**2 - (P-R)**2)/(2*C*K))
        Jt[i]=(K**2 + R**2-(2*K*R*(math.cos(THETAk[i]))))**0.5
        PSIt[i]=math.acos((C**2 + Jt[i]**2 - P**2)/(2*C*Jt[i]))-math.asin((R/Jt[i])*math.sin(THETAk[i]))
        Sprt[i]=A*(PSIbt-PSIt[i])
    elif TYPE==2: #markII
        THETAk[i]=THETAtt[i]+(CW*PHIt)
        PSIbt=math.acos((C**2+K**2-(P-R)**2)/(2*C*K))
        PSItt=math.acos((C**2+K**2-(P+R)**2)/(2*C*K))
        Jt[i]=(K**2+R**2-2*K*R*math.cos(THETAk[i]))**0.5
        BETAt[i]=((C**2+P**2-Jt[i])/(2*C*P))
        PSIt[i]=math.asin(P*math.sin(BETAt[i])/Jt[i])-math.asin((R/Jt[i])*math.sin(THETAtt[i]-PHIt))
        Sprt[i]=-A*(PSIbt-PSIt[i])
    else:
        THETAk[i]=THETAtt[i]+(CW*PHIt)
        PSIbt=math.acos((C**2+K**2-(P-R)**2)/(2*C*K))
        PSItt=math.acos((C**2+K**2-(P+R)**2)/(2*C*K))
        Jt[i]=(K**2+R**2-2*K*R*math.cos(THETAk[i]))**0.5
        BETAt[i]=((C**2+P**2-Jt[i])/(2*C*P))
        PSIt[i]=math.asin(P*math.sin(BETAt[i])/Jt[i])+math.asin((R/Jt[i])*math.sin(THETAtt[i]-PHIt))
        Sprt[i]=-A*(PSIbt-PSIt[i])

print(Sprt)
plt.plot(Sprt)
plt.show()