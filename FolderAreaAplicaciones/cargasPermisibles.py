# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 17:27:26 2022

@author: jose.contreras
"""

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
import time 


import openpyxl
import math
import scipy.signal
from   scipy.stats     import skew
from   scipy.stats     import kurtosis
from   scipy.spatial   import distance
from   scipy.integrate import trapz
from   statistics      import mode

PI = 3.1416
TYPE=1
n=140 #numero de puntos
#Medidas de geometría del balancín
#convencional
A=222
C=186
I=124
P=135.8
H=209.5
G=94.1
R=31.5
# A=129
# C=111.072
# I=111
# P=132
# H=232
# G=96
# R=42
# A=336
# C=134.5
# I=130
# P=261.5
# H=303.5
# G=42
# R=47
#MarkII
# A=384
# C=306
# I=228
# P=234.38
# H=346
# G=132
# R=80.06
#AirBalance
# A=115
# C=48
# I=46.5
# P=114
# H=132
# G=18
# R=13.31
# A=336
# C=134.5
# I=130
# P=261.5
# H=303.5
# G=42
# R=47
CW=1 #horario positivo antihorario negativo
# =============================================================================
# #===========CALCULOS DE ANGULOS TEORICOS===========================#
# =============================================================================
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
    THETAktopt=math.acos((K**2+(P+R)**2-C**2)/(2*K*(P+R)))
    THETAtopt=THETAktopt-(PHIt-PI)
else:  #Airbalance
    PHIt=PI - math.atan(I/(H-G))
    THETAktopt=math.acos((K**2+(P-R)**2-C**2)/(2*K*(P+R))) 

print(PHIt)
#crear el vector de tamaño n desde 0+THETAtopt hasta 2*PI+THETAtopt

THETAtt=np.arange(0, 2*PI, 2*PI/(n))
THETAtt1=abs(2*PI-THETAtt)


# =============================================================================
# #===============================Calculos de ángulos====================#
# =============================================================================
#inicio los vectores 
BETAt  = np.zeros(n)
PSIt   = np.zeros(n)
PSIbt  = np.zeros(n)
PSItt  = np.zeros(n)
ALPHAt = np.zeros(n)
Sprt   = np.zeros(n)
Jt     = np.zeros(n)
THETAk = np.zeros(n)

for i in range (0, n):
    
    if TYPE==1: #convencional
        THETAk[i] = THETAtt[i]-(CW*PHIt)
        PSIbt     = math.acos((C**2 + K**2 - (P+R)**2)/(2*C*K))
        PSItt     = math.acos((C**2 + K**2 - (P-R)**2)/(2*C*K))
        Jt[i]     = (K**2 + R**2-(2*K*R*(math.cos(THETAk[i]))))**0.5
        BETAt[i]  = np.arccos((C**2+P**2-Jt[i]**2)/(2*C*P))
        PSIt[i]   = math.acos((C**2 + Jt[i]**2 - P**2)/(2*C*Jt[i]))-math.asin((R/Jt[i])*math.sin(THETAk[i]))
        ALPHAt[i] = BETAt[i] + PSIt[i] - (THETAtt[i] - PHIt)
        Sprt[i]   = A*(PSIbt-PSIt[i])
    elif TYPE==2: #markII
        THETAk[i]=THETAtt[i]-(CW*PHIt)
        PSIbt=math.acos((C**2+K**2-(P-R)**2)/(2*C*K))
        PSItt=math.acos((C**2+K**2-(P+R)**2)/(2*C*K))
        Jt[i]=(K**2+R**2-2*K*R*math.cos(THETAk[i]))**0.5
        BETAt[i]=np.arccos((C**2+P**2-Jt[i]**2)/(2*C*P))
        # BETAt[i]=((C**2+P**2-Jt[i])/(2*C*P))
        PSIt[i]=math.asin(P*math.sin(BETAt[i])/Jt[i]) - math.asin((R/Jt[i])*math.sin(THETAtt[i]-PHIt))
        ALPHAt[i] = -BETAt[i] - PSIt[i] + (THETAtt[i] - PHIt)
        Sprt[i]=-A*(PSIbt-PSIt[i])
    else:
        THETAk[i]=THETAtt[i]-(CW*PHIt)
        PSIbt=math.acos((C**2+K**2-(P-R)**2)/(2*C*K))
        PSItt=math.acos((C**2+K**2-(P+R)**2)/(2*C*K))
        Jt[i]=(K**2+R**2-2*K*R*np.cos(THETAk[i]))**0.5
        BETAt[i]=np.arccos((C**2+P**2-Jt[i]**2)/(2*C*P))
        # BETAt[i]=((C**2+P**2-Jt[i])/(2*C*P))
        PSIt[i]=math.asin(P*math.sin(BETAt[i])/Jt[i])+math.asin((R/Jt[i])*math.sin(THETAtt[i]-PHIt))
        ALPHAt[i] = BETAt[i] + PSIt[i] + (THETAtt[i] - PHIt)
        Sprt[i]=-A*(PSIbt-PSIt[i])

print(Jt)
plt.figure
# plt.plot(ALPHAt)
# plt.plot(Jt)
# plt.plot(THETAk)
# plt.plot(BETAt)
plt.plot(PSIt)
# plt.plot(Sprt)
# plt.plot(THETAtt)
plt.show()

# =============================================================================
#                          calculo cargas permisibles 
# =============================================================================
#definir
T_crank = 0
M_cb    = 0
Dcb     = 0
N_cb    = 0
W_cbn   = 0
N_cba   = 0
W_cba   = 0
GB      = 320000 #GearBox
SU      = 550
Beam    = 29800
T_cb    = np.zeros(len(THETAtt))
TF      = np.zeros(len(THETAtt))
PL      = np.zeros(len(THETAtt))

for i in range (0, len(THETAtt)):
    # T_cbmax = T_crank * (M_cb - Dcb) * (N_cb * W_cbn + N_cba * W_cba)
    T_cbmax = 705954.7
    T_cb[i] = T_cbmax * np.sin(THETAtt[i])
    TF[i] = (R * A * np.sin(ALPHAt[i])) / (C * np.sin (BETAt[i]))
    PL[i] = SU + ((GB + T_cb[i])/TF[i])
    if PL[i] > Beam:
        PL[i] = Beam

# =============================================================================
#                   Lectura de dinagrama
# =============================================================================

def Norma_Dyna_load(Load):# Normalizar dinagrama
    np = len(Load)   # Longitud muestras
    mL = min(Load)  # Min de la carga
    ML = max(Load)
    Load_Sn = (Load-mL)/(ML-mL)
    return Load_Sn

def Norma_Dyna_pos(Pos):# Normalizar dinagrama
    np = len(Pos)   # Longitud muestras
    mX = min(Pos)   # minimo de X
    MX = max(Pos)   # maximo de X
    Pos_Sn  = (Pos-mX)/(MX-mX) # X normalizado
    return Pos_Sn

class return_values:    
    def __init__(self, posicionS, cargaS, posicionF, cargaF, vector):
        self.posicionS = posicionS
        self.cargaS = cargaS
        self.posicionF = posicionF
        self.cargaF = cargaF
        self.vector = vector
    
def LeerExcel(path, sheet_name,header):   
    xl   = pd.ExcelFile(path)
    df   = xl.parse(sheet_name,header)  # sheet name
    vector = df.to_numpy()
    posicionS =  vector[:,0]
    cargaS = vector[:,1]
    posicionF =  vector[:,2]
    posicionF = Norma_Dyna_pos(posicionF)
    cargaF = vector[:,3]
    cargaF = Norma_Dyna_load(cargaF)
    t=return_values(posicionS, cargaS, posicionF, cargaF, vector)
    return  t

path = "DINAGRAMAS_VALIDACION.xls"
Average      = LeerExcel(path,7,2)

# =============================================================================
#                                   ploteo
# =============================================================================

plt.figure
plt.plot(TF)
# plt.xlim(0,100)
plt.show()

plt.figure
plt.plot(Sprt, PL/1000)
plt.plot(Average.posicionS,Average.cargaS)
plt.ylim (-0,29.8)
plt.show()
