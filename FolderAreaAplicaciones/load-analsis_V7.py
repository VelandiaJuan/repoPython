# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 11:22:24 2022

@author: jose.contreras
"""

import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt
import openpyxl
import scipy.signal
import xlsxwriter

def smooth(y,pts):# Signal Smoothing
    box     = np.ones(pts)/pts
    ysmooth = np.convolve(y,box,mode='same')
    
    return ysmooth

def Stroke_Pot1(Pos,Pot, Load, nlineas):           # Detection of stroke using VSD Power
    global PosX,PotX,LoadX, Idx, Idxn, vectorLV
    
    MPow     = Pot[0:nlineas // 2] # Inicia con la mitad de los datos
    max_pow  = np.max(MPow)          # Detección del maximo de la señal de Potencia
    min_pow  = np.min(MPow)
    Idx      = np.argmax(MPow)       # Indice del maximo de la señal

    vectorLV=np.zeros(nlineas)
    Auxz=0
    final=max_pow-(max_pow-min_pow)*0.3
    for i in range (Idx, nlineas):
        if Pot[i]>final:
            vectorLV[Auxz]=i
            Auxz+=1
    Aux5=0
    Aux6=0
    Aux7=0
    Idxn=0
    for i in range (0, nlineas-1):
        if vectorLV[i+1]-vectorLV[i] > 10:
            if Aux5==0:
                Idxj = int(vectorLV[i+1])
                Aux5=1
                ind = i+1
                    
    for i in range (ind+1, nlineas-1):
        if vectorLV[i+1]-vectorLV[i] > 10:
            if Aux6==0:
                Idxn = int(vectorLV[i+1])
                Aux6 = 1
     
    Idx = Idxj
     
    Delta_S     = Idxn - Idx                # numero de datos entre maximos (Stroke)
    PosX        = Pos_F[Idx:Idxn]
    PosX[Delta_S-1] = PosX[0]
    PotX        = Pot_F[Idx:Idxn]
    LoadX       = Load_F [Idx:Idxn]
    LoadX[Delta_S-1] = LoadX[0]

    
    return Idx,Idxn, Delta_S

def Savgol (X, caso):
    global y, ln
    k = 0
    aux_1 = 0
    aux_2 = 0
    aux_3 = 0
    aux_4 = 0
    aux_5 = 0
    
    if caso==1:
        for i in range (2, nlineas-2):
            if aux_1==0:
                y = np.zeros(nlineas-4)
                aux_1 = 1
            y[k] = (1/35)*(-3*X[i-2] + 12*X[i-1] + 17*X[i] + 12*X[i+1] - 3*X[i+2])
            k += 1
        
    if caso==2:
        for i in range (3, nlineas-3):
            if aux_2==0:
                y = np.zeros(nlineas-6)
                aux_2 = 1
            y[k] = (1/21)*(-2*X[i-3] + 3*X[i-2] + 6*X[i-1] + 7*X[i] + 6*X[i+1] + 3*X[i+2] - 2*X[i+3])
            k += 1
        
    if caso==3:
        for i in range (4, nlineas - 4):
            if aux_3==0:
                y = np.zeros(nlineas-8)
                aux_3 = 1
            y[k] = (1/231)*(-21*X[i-4] + 14*X[i-3] + 39*X[i-2] + 54*X[i-1] + 59*X[i] + 54*X[i+1] + 39*X[i+2] + 14*X[i+3] - 21*X[i+4])
            k +=1
        
    if caso==4:
        for i in range (3, nlineas-3):
            if aux_4==0:
                y = np.zeros(nlineas-6)
                aux_4 = 1
            y[k] = (1/231)*(5*X[i-3] - 30*X[i-2] + 75*X[i-1] + 131*X[i] + 75*X[i+1] + -30*X[i+2] + 5*X[i+3])
            k += 1
            
    if caso==5:
        for i in range (4, nlineas - 4):
            if aux_5==0:
                y = np.zeros(nlineas-8)
                aux_5 = 1
            y[k] = (1/429)*(15*X[i-4] - 55*X[i-3] + 30*X[i-2] + 135*X[i-1] + 179*X[i] + 135*X[i+1] + 30*X[i+2] - 55*X[i+3] + 15*X[i+4])
            k +=1
        
    ln = len(y)
    return
   

#==============================================================================
#path = "0Cargas_7SPM_22Feb.xls" # Datos del pozo average del TWM
#path = "3Cargas_7SPM_22Feb.xls"
#path = "2Cargas_7SPM_22Feb.xls"
#path = "1Carga_BP.xls"
#path = "3Cargas_7SPM_24Feb.xls"
#path = "DobleVel_3Cargas_7SPM_24Feb.xls"
path = 'Reverse_3Cargas_7SPM_25Feb.xls'
xl   = pd.ExcelFile(path)
df   = xl.parse(sheet_name=1,header=0)  # sheet name
Oil_array = df.to_numpy()
nlineas = (len(Oil_array)) # Numero de lineas de los vectores 0..Tmatriz
nl=nlineas
N_puntos=nl

Pos        = np.zeros(nlineas) # posicion
Load       = np.zeros(nlineas) # posicion
Load_F     = np.zeros(nlineas) # posicion
Pot        = np.zeros(nlineas)
Pot_F      = np.zeros(nlineas)
Current    = np.zeros(nlineas)

for i in range(nlineas):            # Estructura de datos en la hoja excel
    
        Pos[i]       = Oil_array[i,0] # carga de datos en las variables      
        Load[i]      = Oil_array[i,1] # correspondientes
        Pot[i]       = Oil_array[i,2]
        #Current[i]   = Oil_array[i,3]
#==============================================================================

Tipo = 1 # selección por posición
Load_F      = scipy.signal.savgol_filter(Load,25,3)
Pos_F       = scipy.signal.savgol_filter(Pos,25,3)
Pot_F       = scipy.signal.savgol_filter(Pot,25,3)
Current_F       = scipy.signal.savgol_filter(Current,25,3)
Current_F       = scipy.signal.savgol_filter(Current_F,25,3)


Savgol(Pot, 3)

plt.figure
#plt.plot(Pot, 'r')
plt.plot(y)
plt.show()

Stroke_Pot1(Pos, y, Load, ln)
nlineas = ln

plt.figure
plt.plot(PotX, 'r')
plt.twinx()
plt.plot(PosX)
plt.show()

plt.figure
plt.plot(PosX, LoadX)
plt.show()