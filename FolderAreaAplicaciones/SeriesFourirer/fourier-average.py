"""
           Created on Thu Jan 20 13:43:09 2022
@author: Alexander.Zambrano - BCPGroup R&D Department
"""
#==============================================================================
#                 Program name : fourier-average.py
#      Evaluar el FFT con datos del banco de pruebas de SRP
#                        y processing de la data
#==============================================================================
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.nonparametric.kernel_regression import KernelReg
from   scipy.fft import fft, fftfreq
from   scipy     import pi
import scipy.signal
import openpyxl
from   scipy.fftpack import fft
#==============================================================================
#                     Detection of Stroke Routines
#==============================================================================
def smooth(y,pts):
    box     = np.ones(pts)/pts
    ysmooth = np.convolve(y,box,mode='same')
    
    return ysmooth

def ZerStro(Prox,VectS,nlineas):
    j=0
    k=0
    for i in range(0,nlineas-2):

        if Prox[i]==1 and Prox[i+1] == 1 and Prox[i+2] == 0: 
           j=i+1
           print("**",j)    
           VectS[k]= j
           k=k+1
        else:
            continue

    return

def MSTROKE(Power,VSD,Pos,VectS):#Deteccion numero de puntos por stroke
    
    TP1    = int(VectS[1])
    TP2    = int(VectS[2])
    NDatos = int(TP2-TP1)
    Mit    = (NDatos//2)
    PosX   = np.zeros(NDatos)
    VSDX   = np.zeros(NDatos)
    PowerX = np.zeros(NDatos)
    k=0
    for i in range(TP1+Mit,TP2+Mit):
        PosX[k]  = Pos[i]
        VSDX[k]  = VSD[i]
        PowerX[k]= Power[i]
        k=k+1
               
    return NDatos,k,PosX,VSDX,PowerX


def Stroke(Pos,VCD,nlineas): # Detection of stroke and numbers of points
    MPos     = Pos[0:nlineas//2]
    min_pos  = np.min(MPos)    # Deteccion del minimo de la señal de Pos
    Idx      = np.argmin(MPos) # Posicion del minimo
    print("primer minimo",Idx," ",MPos[Idx])
    
    i=1 #Counter of zeros 
    n=0 # numbers of stroke points  
    k=0
    Data_Pos = np.zeros(nlineas)
    
    for j in range(Idx+1,(nlineas)):
        if Pos[j] == min_pos:
            print("++",j,Pos[Idx])
            n=0
            
        elif Pos[j]<=Pos[Idx]:
            Data_Pos[k]= Pos[j]
            k=k+1
            n=n+1
            print(n,k,Data_Pos[j])
            
    print("***",Idx,"***",j,k,min_pos,Data_Pos[k])    
        
    return j,k,min_pos,Data_Pos


#==============================================================================
#                   LEER DATOS EN EXCEL DEL BANCO DE PRUEBA SRP-SLACOL
#                         Power,Amp,RPM,Pos,veloc,Load 
#==============================================================================
#path = "E:/SUCKER ROD PUMP/power-amp-averg.xls"
path = "DATA_2CARGAS_28ENERO.xls"
xl   = pd.ExcelFile(path)
df   = xl.parse(sheet_name=0,header=0)  # sheet name
Oil_array = df.to_numpy()

nlineas = len(Oil_array)
#==============================================================================
Sample     = np.zeros(nlineas)
Power      = np.zeros(nlineas)
CT         = np.zeros(nlineas)
VSD        = np.zeros(nlineas)
Pos        = np.zeros(nlineas)
Vel        = np.zeros(nlineas)
Load       = np.zeros(nlineas)
Data_pos   = np.zeros(nlineas)   
Prox       = np.zeros(nlineas)
SPM        = np.zeros(nlineas)
VectS      = np.zeros(nlineas)

#==============================================================================
#          GENERACION DE VECTORES PRINCIPALES PARA N VARIABLES DEL EXCEL
#==============================================================================

for i in range(0,nlineas): # CFilas
        Sample[i]      = i*0.05
        Power[i]       = Oil_array[i,9]
        # Amp[i]       = Oil_array[i,1]
        VSD[i]         = Oil_array[i,6]
        CT[i]          = Oil_array[i,7]
        Pos[i]         = Oil_array[i,12]
        Prox[i]        = Oil_array[i,10]
        SPM[i]         = Oil_array[i,11]
#       Load[i]        = Oil_array[i,9]
#==============================================================================
sample_rate = 1024
N           = (2 - 0) * sample_rate
time_S      = np.linspace(0,0.5,nlineas)
time        = np.linspace(0,2,N)
freq1       = 60
magnitude1  = 25
freq2       = 270
magnitude2  = 2

#waveform1  = magnitude1 * np.sin (2 * pi * freq1 * time)
#waveform2  = magnitude2 * np.sin (2 * pi * freq2 * time)

noise       = np.random.normal (0,10,nlineas)

#time_data  = waveform1 + waveform2 + noise # noise added to signal
Amp = VSD
Pos_F       = scipy.signal.savgol_filter(Pos,5,3) # Filtro de la señal
VSD_F       = scipy.signal.savgol_filter(VSD,5,3)
plt.figure(1)
plt.plot   (Sample,Pos_F)
plt.plot   (Sample,VSD_F)
plt.title  ('Time Domain Signal')
plt.xlabel ('Time')
plt.ylabel ('Amplitude')
T =0.02 # sampling
f = 1/T # Hz

N=nlineas
print(N)
ZerStro(Prox,VectS,nlineas)
[NDatos,kk,PosX,VSDX,PowerX] = MSTROKE(Power,VSD,Pos,VectS)
Torq       = np.zeros(NDatos)
#Pos_XF    = scipy.signal.savgol_filter(PosX,11,3) # Filtro de la señal
Pos_XF     = smooth(PosX,30)
#Power_XF  = scipy.signal.savgol_filter(PowerX,11,3) # Filtro de la señal
Power_XF   = smooth(PowerX,30)
print("--------------------------------------------------------")
[j,k,min_pos,Data_pos]= Stroke(Pos,VSD,nlineas)

DPos = np.diff(Pos_XF)
DPos = scipy.signal.savgol_filter(DPos,11,3)
print("N-puntos:=",kk)
eff = 0.85 # Efficiencia motor
Torq = 84484*eff*Power_XF/SPM[2]
#=============================================================================
#                 Transform the signal in the frequency domain
#=============================================================================
x = np.linspace (0.0,f/2,NDatos)
freq_data  = fft(PosX)
yf         = NDatos * np.abs (freq_data)
xf = fftfreq(NDatos,T)[:NDatos//2]
plt.figure(2)
plt.plot(xf,2/NDatos *np.abs(yf[0:NDatos//2]))
plt.title( 'Frequency Domain Signal')
plt.xlabel('Frequency in Hz')
plt.ylabel('Amplitude')
#==============================================================================
plt.figure (3)
plt.plot(Pos_XF,'r',linewidth=2)
plt.twinx()
plt.plot(Power_XF,'b',linewidth=2)
plt.xlabel('numero de muestras')
plt.ylabel('Potencia-KW')
#===================================ZZ=========================================
ysmooth = smooth(DPos,3) # using convolution

plt.figure(4)
plt.plot(ysmooth)
plt.plot(smooth(DPos,30))
plt.xlabel(' Numero de muestras')
plt.ylabel(' Velocity [in/sec]')
#==============================================================================
xl = np.arange(0,NDatos-1,1)
kr = KernelReg(DPos,xl,'c') # regresion
plt.figure(5)
plt.title( 'Velocity Analysis')
plt.plot(xl,DPos,'+')
DPos_P,y_std = kr.fit(xl)
plt.plot(xl,DPos_P,'r')
plt.xlabel('numero de muestras')
plt.ylabel('Velocity fitting [in/sec]')
#==============================================================================
plt.figure(6)
plt.plot(Power_XF)
plt.twinx()
plt.plot(Torq)
plt.show()
#==============================================================================