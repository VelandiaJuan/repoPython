#==============================================================================
#                        DEPARTAMENTO DE R&D DE SLACOL
#                               FEBRERO -2022-
#==============================================================================
#                 ANALISIS DE LA GEOMETRIA DEL SRP - BANCO DE PRUEBA TIPO C
#==============================================================================
#                   Program name : SIERRACOL-GEOMETRY.py
#      Evaluar las señales principales con datos del banco de pruebas de SRP
#                        y processing de la data
#==============================================================================
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
from   statsmodels.nonparametric.kernel_regression import KernelReg
from   scipy.fft         import fft, fftfreq
from   scipy             import pi
from   scipy.ndimage     import uniform_filter1d
import scipy.signal
import openpyxl
from   scipy.fftpack import fft
#===========================================================================
from   tkinter     import *
from   tkinter.ttk import *
from   pylab       import *
from   numpy       import trapz
import math
from   statistics  import mean
from six.moves import tkinter as tk
import time 
#==============================================================================
#                     Detection of Stroke Routines
#==============================================================================
def smooth(y,pts):# Signal Smoothing
    box     = np.ones(pts)/pts
    ysmooth = np.convolve(y,box,mode='same')
    
    return ysmooth

def ZerStro(Prox,VectS,nlineas):# Stroke Detection from proximeter
    j=0
    k=0
    for i in range(0,nlineas-2):

        if Prox[i]==1 and Prox[i+1] == 1 and Prox[i+2] == 0: 
           j=i+1
           print("**",j)    
           VectS[k] = j
           k=k+1
        else:
           continue

    return

def MSTROKE(Torq_m,Power,VSD,Pos,VectS):# Best Stroke Detection using proximeter
    
    TP1     = int(VectS[1])
    TP2     = int(VectS[2])
    NDatos  = int(TP2-TP1)
    Mit     = (NDatos//2)
    N1      = TP1+Mit
    N2      = TP2+Mit
    PosX    = np.zeros(NDatos)
    VSDX    = np.zeros(NDatos)
    PowerX  = np.zeros(NDatos)
    Tsample = np.zeros(NDatos)
    TorqX   = np.zeros(NDatos)
    PuntosX = np.zeros(NDatos)
    PuntosX = range(N1,N2)
    k=0
    for i in range(TP1+Mit,TP2+Mit):
        PosX[k]    = Pos[i]
        VSDX[k]    = VSD[i]
        PowerX[k]  = Power[i]
        Tsample[k] = k*0.027 # One Stroke sampling 
        TorqX      = Torq_m[i]
        k=k+1
               
    return NDatos,k,PosX,VSDX,PowerX,TorqX,Tsample,PuntosX,N1,N2


def Stroke(Power,nlineas): # Detection of stroke using VSD Power
    MPow     = Power[0:nlineas // 2] # Inicia con la mitad de los datos
    max_pow  = np.max(MPow)        # Deteccion del maximo de la señal de Potencia
    Idx      = np.argmax(MPow)     # Indice del maximo de la señal

    Data_Power = np.zeros(nlineas)
    nmpow      = Power[Idx+1:nlineas] # segundo VEctor para buscar el max power
    max_pow_n  = np.max(nmpow)
    Idxn       = np.argmax(nmpow)
    Idxn       = Idxn+Idx
    Delta_S    = Idxn - Idx                # numero de datos entre maximos (Stroke)
    Data_Power = Power[Idx:Idxn]           # Caaptura del Power en un stroke
      
    return Idx,Idxn,max_pow,max_pow_n,Delta_S,Data_Power

def StrokeX(SX,nlineas):          # Detection of stroke using Pos-Geometrico
    DNX      = nlineas//4               # Inicia con la mitad de los datos
    MX       = SX[DNX : nlineas]
    minX     = np.min(MX)         # Deteccion del minimo de la señal de Pos
    IdX      = np.argmin(MX)      # Indice del minimo de la señal
    IdX      = IdX+DNX
# Buscando el segundo minimo de la posición - geometria
    nmpos      = SX[IdX+1:nlineas]
    min_X_n    = np.min(nmpos)
    IdXn       = np.argmin(nmpos)
    IdXn       = IdXn + IdX
    DXX        = IdXn - IdX
    XPos       = np.zeros(DXX)
    XPos       = SX[IdX:IdXn] # segundo VEctor para buscar el min POS-X
      
    return IdX,IdXn,minX,min_X_n,DXX,XPos
#=================================================================================================
def Nodal(casos,PWf_test,Q_actual,Pcas,Zperf,Zbomba,NDFD_test,Pyac,Pb,Q_test,Grad_P):
#=================================================================================================
    if casos == 1:
#		 print ("Datos provenientes de la prueba de pozo")
        
         IP    = Q_test/(Pyac-PWf_test)
         PIP   = PWf_test-(Zperf - Zbomba) * Grad_P
         Z_PIP = PIP/Grad_P
         SUM   = (PIP - Pcas) / Grad_P
         Nivel = Zbomba - SUM
         PWf   = PWf_test
         
        
    elif casos == 2: # DATOS PROVENIENTE DEL ECHOMETER
         
         print("Datos provenientes del Echometer")
         SUM = Zbomba - NDFD_test
         PWf = Pcas + (Grad_P * (Zperf - NDFD_test))
         PIP = PWf - (Grad_P * (Zperf - Zbomba))
         Z_PIP=PIP/Grad_P
         Nivel = NDFD_test
         IP= Q_test/(Pyac - PWf)
    else:
         print("Casos invalidos")	         
         PIP=Z_PIP=PWf=SUM=Nivel=IP = 0
    return PIP,Z_PIP,PWf,SUM,Nivel,IP
#=================================================================================================
 
#            RUTINA PARA GRAFICAR EL NODAL CONSIDERANDO EL TIPO DE YACIMIENTO SATURADO

#=================================================================================================
#              GRAFICO DEL NODAL YAC-SATURADO 
#=================================================================================================
def Graph_Nodal_Saturado(Pyac,Pb,PWf,Qp,Q_actual):
    Presion = np.zeros((299),dtype=float)
    Qo=np.zeros((299),dtype=float)
    M=0
          
    Presion[0]=Pyac
    Qmax = Qp/(1.0 -(0.2*(PWf/Pyac))-0.8*(PWf/Pyac)**2)
    dif=100.1
    Tyac=0   
    for kk in range(0,299):
                if (Q_actual > 0 and Q_actual < Qmax) and Presion[kk] >= 0:
                    
                        Qo[kk]=Qmax*(1.0 -(0.2*Presion[kk]/Pyac) - 0.8*(Presion[kk]/Pyac)**2)
                        dif = Q_actual - Qo[kk]
                        Presion[kk+1]=Presion[kk]-5
                        print (" Pr=", " ",Presion[kk]," ",dif)
                        
                        if dif < 0.1 and M == 0: # error de deteccion del PT
                            
                                kk=kk-1
                                PWf=Presion[kk]
                                Qn=Qo[kk]
                                M=1
                                dif=100.1
                                print ("Potencial","  ",PWf,"  ",Qn)
                        else:
                            continue
                                
                       
                else:
                    break #dif=100.1
                        
                
    print(" Se calculo todo el rango del nodal")
    return PWf,Qn,kk-1,Qo,Presion,Tyac 
#=================================================================================================
    
# RUTINA PARA GRAFICAR EL NODAL CONSIDERANDO EL TIPO DE YACIMIENTO SUB-SATURADO
    
#================================ANALISIS CASO YAC SUB-SATURADO ==================================  
  
def Graph_Nodal_Sub_Saturado(Pyac,Pb,PWf,Qp,Q_actual):
    Presion = np.zeros((299),dtype=float)
    Qo=np.zeros((299),dtype=float)
    M=0
    Tyac=1
    
    if PWf > Pb: # Zona liquido
        IPS= Qp/(Pyac-PWf)
    else:
        IPS=Q_actual/(Pyac-Pb+(Pb/1.8)*(1.0 -0.2*(PWf/Pyac)-0.8*(PWf/Pyac)**2))
        
    QbS=IPS*(Pyac-Pb)
    Qmax=QbS+(IPS*Pb)/1.8
    Qo_actual = IPS*(Pyac-PWf)
    if PWf < Pb:
        Qo_actual = QbS*(Qmax-QbS)*(1.0 -0.2*(PWf/Pb)-0.8*(PWf/Pb)**2)
    
    Presion[0]=Pyac
    Qo[0]=0
    dif=100.1
    M=0
    for kk in range(0,299):
            if Presion[kk] >= 0:
                
                    if Presion[kk] > Pb:
                        Qo[kk]=IPS*(Pyac-Presion[kk])
                        Presion[kk+1]=Presion[kk]-5
                        dif=Qo_actual-Qo[kk]
                        #print ("PWf>Pb"," ",Qo[kk]," ",Presion[kk], " ",dif," ",kk)
                    else:
                       
                        Qo[kk]=QbS+(Qmax-QbS)*(1.0 -(0.2*(Presion[kk]/Pb))-(0.8*(Presion[kk]/Pb)**2))
                        #print ("PWf<Pb"," ",Qo[kk]," ",Presion[kk], " ",dif," ",kk)
                        Presion[kk+1]=Presion[kk]-5
                        dif=Qo_actual-Qo[kk]
                    if dif < 0.1 and M==0: # Punto de trabajo del pozo
                        PWf=Presion[kk-1]
                        Qn = Qo[kk-1]
                        PT=kk-1
                        M=1
                        print ("PT",PWf,"  ",Qp," ",PT,"  ",dif," ",kk)
                        dif=100.1
                
                    else:
                        print ("Buscando el punto de trabajo")
                        #dif=100.1
                        #M=1
            else:
                print("Presion =",Presion[kk]," ",kk)
                break      
    print("Generación nodal correcta")
    
            
    return PWf,Qn,kk-1,Qo,Presion,Tyac 
    
#=================================================================================================
#                            Graph_nod al End Rutine 
#             Analisis de situacion de Pump-Off <br>
#         Se consideran las variables:Grad_P, NDFD y PIP <br>
#             Criterio es un 10% por encima del PIP <br>
#               Estatus Pump-Off VPOFF=0 , sino VPOff=1
#=================================================================================================

def Prueba_POff(Grad_P,NDFD,PIP):
    Pnivel = NDFD * Grad_P
    Delta_presion = Pnivel-PIP
    CPOff = (10 *PIP/100)+PIP
    if Delta_presion < CPOff:
        VPOFF=0
        print("Pozo en Pump-Off")
    else:
        VPOFF=1
    return CPOff,Delta_presion,Pnivel,VPOFF
    
def Prueba_variables(VPOff,BSW,BSW_cr,Corriente,Corriente_cr,PIP,PIP_cr,Vib,Vib_cr,HZ,HZ_cr):

    if Vib >=  0: # Existe la medida del sensor de fondo de vibración
        if (VPOff == 0) or (BSW > BSW_cr) or (Corriente > Corriente_cr) or (PIP > PIP_cr) or (Vib > Vib_cr) or (HZ > HZ_cr):
            print("No hay recomendacion sobre incremento de produccion")
            print ("paso 1",BSW,">",BSW_cr,PIP,">",PIP_cr,Vib,">",Vib_cr,Corriente,">",Corriente_cr,HZ,">",HZ_cr,"CPOFF=",VPOff)
            recomendacion=0
        else:
            print("Posibilidad de recomendacion sobre incremento de produccion")
            recomendacion=1    
    else:# analisis sin el sensor de vibración
        if (VPOff == 0) or (BSW > BSW_cr) or (Corriente > Corriente_cr) or (PIP > PIP_cr) or (HZ > HZ_cr):
            print("No hay recomendacion sobre incremento de produccion")
            recomendacion=0
        else:
            print("Posibilidad de recomendacion sobre incremento de produccion")
            recomendacion=1
    return recomendacion 
#=================================================================================================
#                 ANALISIS FACTIBILIDAD DE INCREMENTAR PRODUCCION    
#=================================================================================================

def Analisis_Produccion(Pyac,Pb,PWf,Q_test,HZ,HZ_ref,Q_actual,Delta_Hz,Tipo_Yac):
    Q_nuevo = Q_actual*(HZ+Delta_Hz)/HZ_ref
    HZ_n=((HZ+Delta_Hz)/HZ_ref)+HZ
    
    if Tipo_Yac == 1: # Sub-saturado
        [PWf,Qn,kk,Qo,Presion,Tyac]= Graph_Nodal_Sub_Saturado(Pyac,Pb,PWf,Q_nuevo,Q_nuevo)
    else:# Sub-saturado
        [PWf,Qn,kk,Qo,Presion,Tyac]= Graph_Nodal_Saturado(Pyac,Pb,PWf,Q_nuevo,Q_nuevo)
    return PWf,Qn,Q_nuevo,HZ_n,Qo,Presion,kk 

#=================================================================================================
#                             ACTIVACION DE BOTONES Y LECTURA POR GUIL
#================================================================================================= 
    
def Proc_Entry():# Procedencia del dato para analisis productividad del pozo
    #global DD
            
    def selec():
        global DD
        DD = IntVar()
        monitor.config(text="Opcion del dato {}".format(opcion.get()))
        DD= opcion.get()
        print("Data from user:=",DD)
        
        
    def reset_radio():
        opcion.set(None)
        monitor.config(text="")
        root.destroy()
        
        
    root = tk.Tk()
    root.title('Seleccion procedencia de datos')
    root.geometry('300x100')

    opcion = IntVar()
    
    rad1 = Radiobutton(root,text='IP', variable=opcion,value=1,command=selec).pack( anchor=W)
    rad2 = Radiobutton(root,text='Echometer', variable=opcion,value=2,command=selec).pack( anchor=W)
    rad3 = Radiobutton(root,text='SS', variable=opcion,value=3,command=selec).pack( anchor=W)
               
    #btn = Button(window, text = "Click Me", command=clicked)
    
    #rad1.grid(row=0, column=0)
    #rad2.grid(row=0, column=1)
    #rad3.grid(row=0, column=2)
    
# btn.grid(row=0, column=3)
# def clicked():
    monitor = Label(root)
    monitor.pack()
    boton = tk.Button(root,text='Continuar aplicacion REX',command=reset_radio)
    boton.pack()
    
    root.mainloop()
    
def velocity(N_puntos,X,tiempo):

#=============================================================================================    
#           Calculo de las velocidades por geometria y sensor de posición
#==============================================================================================
    global V_geometry,a_geo,vel_Xpr
    V_geometry    = np.zeros(N_puntos)
    a_geo         = np.zeros(N_puntos)
    Ft            = 0.0833 # conversion a ft
    m             = N_puntos
    Xm            = X[0]-X[m-1]
    DT            = tiempo[1]-tiempo[0]
    V_geometry[m-1] = (Xm/DT)*Ft 
    vel_Xpr[m-1]    = (Xpr[0]-Xpr[m-1])/DT
                            
    for i in range(0,N_puntos-1):
        V_geometry[i]= (X[i+1]-X[i])/(tiempo[i+1]-tiempo[i])
        vel_Xpr[i]   = (Xpr[i+1]-Xpr[i])/(tiempo[i+1]-tiempo[i]) # in/sec
    print(V_geometry)
    print(X)  
                 
    #a_geo[m-1]= (V_geometry[0]-V_geometry[m-1])/(tiempo[0]-tiempo[m-1]) 
                
    for i in range(0,N_puntos-1):
        a_geo[i]= (V_geometry[i+1]-V_geometry[i])/(tiempo[i+1]-tiempo[i])
        
    print(a_geo)             
 
   
#============================================================================================== 
#                          PROGRAMA PRINCIPAL ESTRUCTURA DE DATOS DEL SRP
#===============================================================================================
#                                     VERSION SWC-BULL
#===============================================================================================
Data_campo   = {"Q_Oil":27,"Q_Wat":60,"Q_Gas":40,"Q_actual": 87,"BSW": 68.9,"Pcas":127.9,"Ptub":50,"API":30,"Tsup":60,"NDFD":3000}
Data_fluidos = {"Visco_fluido":0.25,"Ggas": 0.9, "Gwat":1.05,"GOR":76,"SGg":1.05,"CpT1":2000,"CpT2":200,"T1":220,"T2":400}
Data_Comple  = {"Anchor_Depth":5035,"Zdatum": 5221,"Zperf":5221,"Zbomba":5115,"Zcas":5200,"Ztub":5115,"Dcas":4.95,"Dtub":1.995}
Data_Yac     = {"Pyac":1485,"Pest":975.79,"Pb":80,"PWf":437,"PIP":239,"Tyac":150,"IP":23}
Data_CILA    = {"Corriente": 89,"Voltaje":480,"V_BUS":600}
#===========================================================================================
#                                   PRUEBAS DE PRODUCCION
#===========================================================================================
Data_test = {"Q_test":8619,"BSW_test":59,"NDFD_test":1171.5,"PWf_test":601.6,
             "HZ_test":55,"THP_test":120,"IP_test":23,"Vib_test":0.4,"Amp_test":90}
#===========================================================================================
#                                   PARAMETROS DE CONTROL
#===========================================================================================
Data_control     = {"PIP_cr":700,"PIP_Optimo":90,"BSW_cr":98,"Corriente_cr":30,"HZ_cr":60}
#===========================================================================================
Data_SRP_ESP     = {"Manufacturer":"Lufkin","SUCKER-ROD":"C-57D-76-42","TIPO":"C"} #Especificaciones Sucker Rod"
Data_costos_elec = {"Pelec":0.05,"Ppde":8,"hpd":24,"dpm":30}
Data_Electrico   = {"HPm":9.2,"Amps":15.7,"RPM":1760,"Volt":480,"N_motor":0.91,"N_belts":0.97,
                    "Rated_RPM":1760,"N_phases":3,"Rated_Hp":9.2,"Rated_Amp":15.7,"Pf":0.90,"Fre":60}
Data_Rod         = {"D_Rod":1.230,"Rod_len":3854,"T_Rod_Weig":8945.25}
Data_Bomba       = {"D_embolo":1.25,"L_bomba":48,"Hol_bomba":0.009}
Data_SRP         = {"Ngb":0.97,"Rated_Torque":320,"GBR":43.6,"SPM":11,"DS":12.6,"SU":150,
                    "S":42,"StC":25.6,"TF":-47.48}
Data_SRP_T1      = {"Rod_t1":"D","L1":1200,"DT1":0.874,"WT1":2654.4,"WT1_U":2.212}
Data_SRP_T2      = {"Rod_t2":"D","L2":3875,"DT2":0.750,"WT2":6290.9,"WT2_U":1.6235}
Data_Dump        = {"Dump_Up":0.02,"Damp_Down":0.02}
Data_Geometricos     = {"A":55,"C":48,"P":57,"R":18,"I":48,"H":96,"G":41.5,"Cw":1.0}
Data_Counter_Weights = {" NCW":4,"NCAW":0,"WCW":50,"WACW":0,"Tcrank":567.25,
                        "MLACW":45.26,"DCW":45.26}
#===========================================================================================
Tyac         = Data_Yac.get("Tyac")
Tsup         = Data_campo.get("Tsup")
BSW          = Data_campo.get("BSW")
BSW_cr       = Data_control.get("BSW_cr")
Gwat         = Data_fluidos.get("Gwat")
Zbomba       = Data_Comple.get("Zbomba")
Zsensor      = Data_Comple.get("Zsensor")
Zperf        = Data_Comple.get("Zperf")
Q_actual     = Data_campo.get("Q_actual")
Pcas         = Data_campo.get("Pcas")
Ptub         = Data_campo.get("Ptub")
NDFD_test    = Data_test.get("NDFD_test")
NDFD         = Data_campo.get ("NDFD")
Pyac         = Data_Yac.get("Pyac")
Pb           = Data_Yac.get("Pb")
Q_test       = Data_test.get("Q_test")
PWf_test     = Data_test.get("PWf_test")
Dtub         = Data_Comple.get("Dtub")
Visco_fluido = Data_fluidos.get("Visco_fluido")
Voltaje      = Data_CILA.get("Voltaje")
Corriente    = Data_CILA.get("Corriente")
Corriente_cr = Data_control.get("Corriente_cr")
HZ_cr        = Data_control.get("HZ_cr")
HZ_test      = Data_test.get("HZ_test")
PIP_cr       = Data_control.get("PIP_cr")

#====================================================================================================
RPM       = Data_Electrico.get("RPM")          # Synchronous RPM
Volt      = Data_Electrico.get("Volt")
N_motor   = Data_Electrico.get("N_motor")# N_motor: Eficiencia del motor
N_belts   = Data_Electrico.get("N_belts")# N_belts  Eficiecia de las correas
Rated_RPM = Data_Electrico.get("Rated_RPM")
N_phases  = Data_Electrico.get("N_phases")
Rated_Hp  = Data_Electrico.get("Rated_Hp")
Rated_Amp = Data_Electrico.get("Rated_Amp")
Pf        = Data_Electrico.get("Pf")# Power Factor
Fre       = Data_Electrico.get("Fre")# Frecuency [hz]
HPm       = Data_Electrico.get("HPm")# Hp
#===================================================================================================
Pelec = Data_costos_elec.get("Pelec")          # precio electricidad $/Kwh
Ppde  = Data_costos_elec.get("Ppde")           # precio potencia demandada $/KVA
hpd   = Data_costos_elec.get("hpd")            # horas/diarias
dpm   = Data_costos_elec.get("dpm")            # dias por mes
#==============================================================================
Ngb  = Data_SRP.get("Ngb")
Rated_Torque= Data_SRP.get("Rated_Torque")# Rated_Torque [Kin-lb]
GBR  = Data_SRP.get("GBR")# Gear Box Ratio
SPM  = Data_SRP.get("SPM")
DS   = Data_SRP.get("DS") # 
SU   = Data_SRP.get("SU") #  Unbalanced 
S    = Data_SRP.get("S")  #  Max length Polished rod [in]
StC  = Data_SRP.get("StC")#  Structure Capacity [Klb]
TF   = Data_SRP.get("TF") #  Torque Factor @90
#==================================================================================================A = Data_Geometricos.get("A")
#                    Datos Geometricos del balancin del banco de prueba
#================================================================================================
A = Data_Geometricos.get("A")
C = Data_Geometricos.get("C")
P = Data_Geometricos.get("P")
R = Data_Geometricos.get("R")
I = Data_Geometricos.get("I")  # Es la distancia horizontal entre el cigüeñal y el cojinete central. (Inch)
H = Data_Geometricos.get("H")  # Altura vertical del cojinete central. (Inch)
G = Data_Geometricos.get("G")  # Altura vertical del cigüeñal. (Inch)
Cw = Data_Geometricos.get("Cw")# Sentido de giro [+1 o -1]
#=================================================================================================

NCW    = Data_Counter_Weights.get("NCW")# Number of counter weights on the crank
NCAW   = Data_Counter_Weights.get("NCAW")# Number of auxiliary weights on the crank
WCW    = Data_Counter_Weights.get("WCW")# weight of each counter weight [lb]
WACW   = Data_Counter_Weights.get("WACW")# weight of each auxiliary weight [lb]
Tcrank = Data_Counter_Weights.get("Tcrank")# Moments of the crank[Kin-lb]
MLACW  = Data_Counter_Weights.get("MLACW")# Maximum lever arm [In]
DCW    = Data_Counter_Weights.get("DCW")# Distance of counterweight from end of the crank
#=================================================================================================
#                          INICIALIZACION DE VARIABLES
#=================================================================================================
global Npuntos,tiempo
plt.clf() # Clear figures
#===============================================================================================
#                   LEER DATOS DE SUPERFICIE DEL BANCO DE PRUEBA SRP
#             Sampling, Potencia_KW, Amp,RPM, Xpr, Lpr, Xbomba, Lbomba 
#===============================================================================================
#==============================================================================================
path      = "DATA_2CARGAS_28ENERO.xls"
xl        = pd.ExcelFile(path)
df        = xl.parse(sheet_name=0,header=0)  # sheet name
Oil_array = df.to_numpy()
nlineas   = len(Oil_array) # numbers of data
#=============================================================================
Sample     = np.zeros(nlineas)
Power      = np.zeros(nlineas)
CT         = np.zeros(nlineas)
VSD        = np.zeros(nlineas)
Pos        = np.zeros(nlineas)
Vel        = np.zeros(nlineas)
Load       = np.zeros(nlineas)
Data_pos   = np.zeros(nlineas)   
Prox       = np.zeros(nlineas)
SPMB       = np.zeros(nlineas)
VectS      = np.zeros(nlineas)
Torq_m     = np.zeros(nlineas)

#=================================================================================================
#          RECUPERACION  DE VECTORES PRINCIPALES DEL EXCEL BANCO DE PRUEBAS
#=================================================================================================
for i in range(0,nlineas): # CFilas
        Sample[i]      = i*0.027         # SPM=11
        VSD[i]         = Oil_array[i,6]  # VSD Current
        CT[i]          = Oil_array[i,7]  # Motor input current
        Power[i]       = Oil_array[i,9]  # Motor Power
        Prox[i]        = Oil_array[i,10] # Proximeter
        SPMB[i]        = Oil_array[i,11] # SRP Velocity [SPM]
        Pos[i]         = Oil_array[i,12] # Jack Pump Position
        Torq_m[i]      = (5252*Power[i]*1.34)/RPM
#       Load[i]        = Oil_array[i,13]

SPM  = SPMB[2]            #  valor del SPM en la data excel de las pruebas

'''
# ============ Cierre de los datos de posicion & esfuerzo/carga================

Xpr[nlineas-1]=Xpr[0]
Lpr[nlineas-1]=Lpr[0]
X_bomba[nlineas-1]=X_bomba[0]
L_bomba[nlineas-1]=L_bomba[0]
                            
print(L_bomba)

#=================================================================================================
#                               GRAPHS OF THE CORE PARAMETERS
#====================================== RUTINAS DE SISTEMA SRP-PRO ===============================

plt.plot(Xpr,Lpr)
plt.show()
'''
#========================================================================
PI           = 3.1416
Rated_Slip   = 100*(RPM-Rated_RPM)/RPM # %
Rated_Torque = (39.3701 * 0.2248*(0.7457*HPm/(2*PI*Rated_RPM/60)))*1000 # Kin-lbf
#=================================================================================================
#                      Calculo del gravedad especifica crudo
#===============  =================================================================================
API = Data_campo.get("API")
#=================================================================================================
#                                   LECTURA DE DATOS
#=================================================================================================
GEo  = 141.5/(131.5 + API)
# ==================== ANALISIS DE LA PRODUCTIVIDAD DEL POZO =====================================


if    GEo < 0.63:
	alfa = 0.00097
	beta = 0.0000004
elif  GEo >= 0.95:
	alfa = 0.00066
	beta = -0.0000003
elif  GEo  >= 0.63 and GEo <= 0.78:
	alfa = 0.00097 + (GEo-0.63)*((0.00075-0.00097)/(0.78-0.63))
	beta = -0.0000004 + (GEo-0.63)*((0.0000004)/(0.78-0.63))
elif  GEo  < 0.85 and GEo >= 0.78: 
	alfa = 0.00075 + (GEo-0.75)*((0.00068-0.00075)/(0.85-0.78))
	beta = (GEo-0.75)*((0.0000001)/(0.85-0.78))

elif  GEo  < 0.95 and GEo >= 0.85:
	alfa = 0.00068 + (GEo-0.85)*((0.00066-0.00068)/(0.95-0.85))
	beta = 0.0000001 + (GEo-0.75)*((0.0000001 - 0.0000003)/(0.85-0.78))

else:
	Out_Range = "valor fuera de rango"
	print(Out_Range)
	

Grave_TF    = GEo -(alfa*(Tyac - Tsup)) + beta*(Tyac-Tsup)**2
Grave_F_BSW = ((100.0-BSW)/100.0)*Grave_TF + (BSW/100.0)*Gwat
Grad_P      = 0.434 * Grave_F_BSW
Dens_oil    = 62.4 * Grave_TF
Dens_fluido = 62.4 * Grave_F_BSW

print (GEo," ",alfa," ",beta," ",Grad_P," ",Grave_TF)

#=================================================================================================
#                        LECTURA DE LA PROCEDENCIA DE LA FUENTE DEL DATO
#=================================================================================================

Proc_Entry() # lectura procedencia del dato para analisis de la productividad
print("DD_Global=",DD)
casos=DD
#casos = int(input("Seleccion procedencia del dato:"))
print (casos)

#casos= 1 prueba de pozo dato=2 Echometer (nivel de fluidos), dato=3 SS

[PIP,Z_PIP,PWf,SUM,Nivel,IP] = Nodal(casos,PWf_test,Q_actual,Pcas,Zperf,Zbomba,NDFD_test,Pyac,Pb,Q_test,Grad_P)
print (PIP,Z_PIP,PWf,SUM,Nivel,IP,sep =' ')
#=================================================================================================
PWf_O=PWf


if Pyac > Pb :# Sub-saturado
    Tipo_Yac=1
    [PWf,Qn,kk,Qo,Presion,Tyac]= Graph_Nodal_Sub_Saturado(Pyac,Pb,PWf,Q_test,Q_actual)
else:
    Tipo_Yac=0 # Saturado
    [PWf,Qn,kk,Qo,Presion,Tyac]= Graph_Nodal_Saturado(Pyac,Pb,PWf,Q_test,Q_actual)
 
print (PWf,Qn,kk,Tipo_Yac,sep='  ')

Pres = np.zeros((kk+1),dtype=float)
QQ   = np.zeros((kk+1),dtype=float)

for i in range(0,kk+1):
    QQ[i] = Qo[i]
    Pres[i] = Presion[i]
    #print(i)
print(Pres) 
#==============================================================================
plt.figure(1) # IWOC
plt.plot(QQ,Pres,'r--',Qn,PWf,'g^',Q_actual,PWf_O,'ko')
plt.grid()
plt.xlabel('Barriles/D')
plt.ylabel('Presion[Psi]')
plt.title('Curva IPR')

#=================================================================================================
#                    ANALISIS DE POTENCIA, TORQUE Y ELECTRICO
#=================================================================================================
Amp = VSD
Pos_F       = scipy.signal.savgol_filter(Pos,5,3) # Filtro de la señal
VSD_F       = scipy.signal.savgol_filter(VSD,5,3)

plt.figure(2)
plt.plot   (Sample,Pos_F,'r',linewidth=2)
plt.twinx()
plt.plot   (Sample,VSD_F,'b',linewidth=2)
plt.title  ('Posicion vs Corriente variador')
plt.xlabel ('Time')
plt.ylabel ('Corriente-Amps')
#==============================================================================
T = 0.027 # sampling
f = 1/T # Hz
#==============================================================================
N=nlineas
print(N)
ZerStro(Prox,VectS,nlineas)
[NDatos,kk,PosX,VSDX,PowerX,TorqX,Tsample,PuntosX,N1,N2] = MSTROKE(Torq_m,Power,VSD,Pos,VectS)
Torq       = np.zeros(NDatos)
#Pos_XF    = scipy.signal.savgol_filter(PosX,11,3) # Filtro de la señal
Pos_XF     = smooth(PosX,20)
#Power_XF  = scipy.signal.savgol_filter(PowerX,11,3) # Filtro de la señal
Power_XF   = smooth(PowerX,30)
VSD_XF     = smooth(VSDX,20)
print("--------------------------------------------------------")
[Idx,Idxn,max_pow,max_pow_n,Delta_S,Data_Power]= Stroke(Power,nlineas) #seleccion stroke

#Data_Power = smooth(Data_Power,5)

DPos = np.diff(Pos_XF)
DPos = scipy.signal.savgol_filter(DPos,11,3)
print("N-puntos:=",kk)
eff  = 0.85 # Efficiencia motor
Torq = 84484*eff*Power_XF/SPM
#==============================================================================
pi      = 180 # Grados
Tm      = float(60/SPM) # periodo de un ciclo de bombeo dado los strokes/min
Tcycle  = Tm/SPM # tiempo de un stroke
Ndiv    = nlineas
Npt     = Tcycle/Ndiv #seg
Vel_ang = 2*PI*SPM/60 # rad/s
SPM_max = 0.7*(60000/S)**0.5
print(Ndiv,Tcycle,Npt,Vel_ang,SPM,SPM_max,nlineas)
print("==================  Angulos    =======================================")
Angulo   = np.arange(0,2*pi,(2*pi)/NDatos) # cada [n] grados/seg 
                                           # se analizan para 60 seg

Rango    = np.arange(len(Angulo))
N_puntos = len(Rango)
print(Angulo)
time.sleep(21)

Xpr     = Pos_XF
Nciclos = 6*N_puntos

print(N_puntos)
tiempo  = np.zeros(N_puntos)
Crank   = np.zeros(N_puntos)
Theta   = np.zeros(N_puntos)
Thet_k  = np.zeros(N_puntos)
JVector = np.zeros(N_puntos)
Chi_V   = np.zeros(N_puntos)
Beta    = np.zeros(N_puntos)
Rho_V   = np.zeros(N_puntos)
Psi_V   = np.zeros(N_puntos)
Alpha   = np.zeros(N_puntos)
Delta_V = np.zeros(N_puntos)
X       = np.zeros(N_puntos)
V       = np.zeros(N_puntos)
TF      = np.zeros(N_puntos)
C1      = np.zeros(N_puntos)
C2      = np.zeros(N_puntos)
C3      = np.zeros(N_puntos)
C4      = np.zeros(N_puntos)
Acc     = np.zeros(N_puntos)
aw      = np.zeros(N_puntos)
CCang   = np.zeros(N_puntos)
Vel_Head= np.zeros(N_puntos)
POT_Vsd = np.zeros(N_puntos)
F_rod   = np.zeros(N_puntos)
POW     = np.zeros(N_puntos)
vel_Xpr = np.zeros(N_puntos)
S_des_m = np.zeros(N_puntos)
FP_int  = np.zeros(N_puntos)
PEM     = np.zeros(N_puntos)
SX      = np.zeros(Nciclos)
RGen_m  = np.zeros(N_puntos) # Regeneracion motor
Con_m   = np.zeros(N_puntos) # Consumo motor
P_crank = np.zeros(N_puntos)
P_rod   = np.zeros(N_puntos)
Inv_V   = np.zeros(N_puntos)
#=================================================================================================
# Calculo del Pitch Diameter [ds]
ds = DS*GBR*SPM/RPM
#=================================================================================================
# Calculo de la distancia [K] entre cigueñal y el cojinete segun esp. LUFKIN
#=================================================================================================
K = (I**2 + ( H-G)**2)**0.5
# Calculo del angulo Phi en radianes
Phi  = math.asin(I/K)
Thet = 1.69
Theta_k= Thet -  Phi # Angulo entre R y K ,Rad
J = (K**2 + R**2-(2*K*R*(math.cos(Theta_k))))**0.5 # Desplazamiento J
Chi = math.acos((C**2 + J**2 - P**2)/(2*C*J)) # X
Rho = math.asin((R/J)*math.sin(Theta_k))*Cw   # angulo entre K y J [p]
Psi = Chi - Rho                               # Angulo entre C y K, Rad
Psi_b = math.acos((C**2 + K**2 - (P+R)**2)/(2*C*K)) # Angulo Psi mas bajo de la manibela
Psi_t = math.acos((C**2 + K**2 - (P-R)**2)/(2*C*K)) # Angulo Psi mas alta de la manibela
Delta = Psi_b - Psi # Ángulo de amplitud de movimiento de la viga viajera del lado de la cabeza de caballo
Delta_max = Psi_b - Psi_t # Ángulo máximo de amplitud de movimiento de la viga viajera
S         = A*Delta # Desplazamiento de la barra pulida in
S_max     = A*Delta_max
Paso      = S_max/60
vel_Xpr   = smooth(vel_Xpr,20)
#=================================================================================================
print("=======================  Vector Tiempo =====================================")


for i in range(NDatos):
    tiempo[i] = Tsample[i]
    Crank[i]  = 6*tiempo[i]*RPM*ds/(DS*GBR)

Ntt=len(tiempo)
print(Ntt)
#=================================================================================================
#                  CALCULO DE LOS ANGULOS POR GEOMETRIA DEL SRP
#=================================================================================================
j=0
for i in Angulo:
    Theta[j]   = i # angulo del crank 0-360
    Thet_k[j]  = Theta[j]*PI/pi-Phi
    JVector[j] = (K**2 + R**2-(2*K*R*(math.cos(Thet_k[j]))))**0.5          # Angle J rad
    Chi_V[j]   = math.acos((C**2 + JVector[j]**2 - P**2)/(2*C*JVector[j])) # X rad
    Beta[j]    = math.acos((C**2 + P**2 - JVector[j]**2)/(2*C*P)) # B rad
    Rho_V[j]   = math.asin((R/JVector[j])*math.sin(Thet_k[j]))*Cw
    Psi_V[j]   = Chi_V[j] - Rho_V[j]
    Alpha[j]   = Cw*(Beta[j] + Psi_V[j]) - Thet_k[j]
    Delta_V[j] = Psi_b - Psi_V[j]
    X[j]       = A*Delta_V[j] # Desplazamiento de la barra pulida [in]
    V[j]       = (A*R*math.sin(Alpha[j])/(C*math.sin(Beta[j])))*Vel_ang # Velocidad [in/seg]
    TF[j]      = A*R*math.sin(Alpha[j])/C*math.sin(Beta[j]) # Torque Factor [In]
    C1[j]      = (A*R*K*(Vel_ang)**2)/(P*C*math.sin(Beta[j])**2)
    C2[j]      = math.cos(Alpha[j])*math.sin(Psi_V[j])
    C3[j]      = R*math.cos(Beta[j])/C*math.sin(Beta[j])
    C4[j]      = math.sin(Alpha[j])*math.sin(Thet_k[j])
    Acc[j]     = C1[j]*(C2[j]+C3[j]*C4[j]) # acceleration in/seg**2
    j=j+1
   
#Max y Min de la posición
X_max = np.amax(X)
X_min = np.amin(X)
DXX   = diff(X) # derivada de X geometrica
#==============================================================================
#               Generacion continua de S = A*DElta_V
#==============================================================================
for jj in range(6):
    for i in range(NDatos):
        SX[jj*NDatos+i]= X[i]


[IdX,IdXn,minX,min_X_n,DXX,XPos]= StrokeX(SX,nlineas)

print("X_max :=","%5.2f" % (X_max)," In")
print("X_min :=","%5.2f" % (X_min)," In")

print("=========================== Vector Posición ===================================")
print(len(XPos))

#=================================================================================================
#            Graficos de los angulos y verificación
#=================================================================================================
plt.figure(3)
plt.plot(Theta,X)
#plot(tiempo,X)
plt.xlabel('Tiempo (seg)')
plt.ylabel('THeta - Angle[Grd]')
plt.title('Kynamatics for SRP - REX')
plt.grid(True)
#=================================================================================================
print("Phi:= ","% 5.2f" % (Phi)," Rad")
print("K:= ","% 5.2f" % (K)," in")
print("Theta_k :=","%5.2f" % (Theta_k)," Rad")
print("J :=","%5.2f" % (J)," Rad")
print("Chi :=","%5.2f" % (Chi)," Rad")
print("Psi :=","%5.2f" % (Psi)," Rad")
print("Psi_b :=","%5.2f" % (Psi_b)," Rad")
print("Psi_t :=","%5.2f" % (Psi_t)," Rad")
print("Delta :=","%5.2f" % (Delta)," Rad")
print("Delta_max :=","%5.2f" % (Delta_max)," Rad")
print("S :=","%5.2f" % (S)," In")
print("S_max :=","%5.2f" % (S_max)," In")

# MInima velocidad del balancin 1S/min
# Paso_des = Desplazamiento Maximo / 60sec 

Paso_Des = S_max/N_puntos
print(Paso_Des)

Desp = np.zeros(N_puntos)
for i in Rango:
    Desp[i] = i*Paso_Des
print(Desp)
#=================================================================================================
#                         Graficos de los angulos y verificacion
#=================================================================================================
plt.figure(4)
plot(Desp,Chi_V*pi/PI)
xlabel('Stroke (In)')
ylabel('Ange Chi (Rad)')
title('Kynamatics for SRP - REX')
grid(True)
#=================================================================================================
#                                   Calculo del angulo Theta
#=================================================================================================
fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(Tsample,X)
axs[0, 0].set_title('Polished Rod - Position')
axs[0, 1].plot(Theta,V, 'tab:orange')
axs[0, 1].set_title('Theta - Velocidad')
axs[1, 0].plot(Theta,TF, 'tab:green')
axs[1, 0].set_title('Theta -Torque Factor')
axs[1, 1].plot(Theta,Acc, 'tab:red')
axs[1, 1].set_title('Theta - Aceleración')
fig.tight_layout()

#=================================================================================================
# Inferencia de la fuerza [F_rod] en la head del SRP a partir de la potencia del  motor
#=================================================================================================
GG  = (((R**2)+(H**2)-(P**2))**0.5)+ A**2
P1  = WCW  # pesos del counterbalance [lb]
velocity(N_puntos,X,tiempo) # Calculo de las velocidades y aceleracion por geometria
vel_Xpr = smooth(vel_Xpr,20)
#Inv_V = 1/vel_Xpr
Inv_V = uniform_filter1d(1/vel_Xpr,size=150)
for i in range(N_puntos):
 
    aw = Vel_ang*tiempo[i]
    aa = C*(C-R*math.cos(aw))
    bb = C*(R*math.sin(aw)-H)
    cc =(aa**2+bb**2)**0.5
    ff = GG-R*C*math.cos(aw)-H*R*math.sin(aw)
    CCang[i] = math.acos(ff/cc)+math.atan(bb/aa)
    #   CCang[i]=Theta[i]*PI/pi
    A1= A*R*Vel_ang/C
    A2= C*math.sin(aw) - H*math.cos(aw)
    A3= C*math.sin(aw+CCang[i])
    A4= R*math.sin(aw+CCang[i]) - (C*math.sin(CCang[i])+ H*math.cos(CCang[i]))
    Vel_Head[i]= A1*(A2-A3)/A4 # in/seg
    conversion = 1/8.851
    P_crank[i] = ((P1/1000)*Vel_ang*R*math.cos(aw))*conversion     # Crank Power
    #F_rod[i]  = (1/A1)*((A4/(A2-A3))*(Power_XF[i]-P1*Vel_ang*R*math.cos(aw)))
    Power_M    =  Power_XF[i] # a lib-in/sec
    P_rod[i]   =  Power_M - P_crank[i]   # KW
    F_rod[i]   =  Inv_V[i]*P_rod[i]       # libras
print(CCang) 
#F_rod = uniform_filter1d(F_rod,size=150)
F_rod  = smooth(F_rod,30)
#==============================================================================
#                                 Derivada del angulo
#==============================================================================
D2Chi = diff(CCang)
D2Chi = A*D2Chi # velocidad en la PR
#==============================================================================
figure(6)
plot(tiempo,Vel_Head)# Velocidad en el rod
#plot(tiempo,V)# Crank Angule

#=================================================================================================
#==============================  ANALISIS DEL SLA DEL POZO =======================================
#=================================================================================================

print("====================== print(Angulo)==================================")
print(Angulo)
print (len(Angulo))
print(N_puntos)

plt.figure(7)
plot(tiempo,V_geometry)

plt.figure(8)
plot(tiempo,a_geo[0:N_puntos])

print(V_geometry)
print(tiempo)

#================================== DATOS DE LOS SENSORES FISICOS DEL POZO
print(nlineas)

figure(figsize=(10, 6), dpi=80)
plot(Tsample, Xpr, color="blue", linewidth=2.5, linestyle="-",label="posicion")
twinx() # crear dos ejes
plot(Tsample, Power_XF, color="red",  linewidth=2.5, linestyle="-",label="Carga-Barra")
plt.xlabel('Tiempo (seg)')
plt.ylabel('Power - KW')
plt.title('Power vs PR Position ')
plt.grid(True)
#==============================================================================
#========================= Potencia Instantanea [Kw] ==========================
#==============================================================================

for i in range(NDatos):
    
    #POW[i]= Lpr[i]*vel_Xpr[i]*1/6.6 #KW
    #S_des_m[i]= 100*(Rated_RPM - RPM_V[i])/Rated_RPM
    FP_int[i]= (Power_XF[i]*1000)/(1.73 *VSDX[i]*Volt)
    #PEM[i]= 1.73*VSDX[i]*Volt/1000 # KVA

#plot(tiempo,vel_Xpr)
#plot(tiempo,POW)
#show()

#==============================================================================
# =========================Potencia Promedio
#==============================================================================

Pow_mean = mean(Power_XF)  # Promedio de la potencia en la barra pulida

print("Potencia_BP :=","%5.2f" % (Pow_mean)," KW")

PEM=Power_XF
#=============================================================================
#      Relación potencia en la barra pulida BP y la corriente motor
#=============================================================================

plt.figure(10)
plot(tiempo,VSDX,color="blue",)
twinx() # crear dos ejes
#plot(tiempo,POW,color="black",)
#show()
plot(tiempo,S_des_m) # % de deslizamiento del motor

plt.figure(11)
plot(tiempo,FP_int)

#==============================================================================
#================ Potencia Electrica Promedio Consumida por el Motor ==========
#==============================================================================

FP=0

for i in range(N_puntos):
    if Power_XF[i] >= 0: # Motorización del motor
        FP = FP+Power_XF[i]
        
    else:
        continue
    
PEC_m =( FP/N_puntos+1) # porcentaje de Motorizacion
print("Potencia Electrica Consumida por el  Motor :=","%5.2f" % (PEC_m)," KW")
print("Potencia Electrica Consumida  por el Motor  :=","%5.2f" % (PEC_m * 1.3410)," Hp")

 
#==============================================================================
#================ Potencia Electrica Promedio Generada por el Motor ==========
#==============================================================================

FP=0

for i in range(N_puntos):
    if Power_XF[i]  < 0: # RE-generación del motor
        FP = FP+Power_XF[i]
        
    else:
        continue
    
PEG_g =(FP/N_puntos+1) # Porcentaje de Motorización
print("Potencia Electrica Generada Motor :=","%5.2f" % (PEG_g)," KW")
print("Potencia Electrica Generada Motor :=","%5.2f" % (PEG_g * 1.3410)," Hp")

#==============================================================================
#================= Potencia Electrica Aparente suministrada al motor
#==============================================================================
plt.figure(12)
plot(tiempo,PEM)

# Potencia Electrica Aparente Pico

PEA_p = max(Power_XF)
print("Potencia Electrica Aparente Motor:=","%5.2f" % (PEA_p)," KW")
Pot_N= PEC_m + PEG_g
print("Potencia Electrica Neta Motor:=","%5.2f" % (Pot_N)," KW") 
#==============================================================================   
#======================= EFiciencia de la unidad de bombeo ====================
#==============================================================================  
Ef_b = (Pow_mean/Pot_N)*100
print("Eficiencia de la unidad Bombeo:=","%5.2f" % (Ef_b)," %") 
PE_app= mean(Power_XF)
print("Potencia Electrica Aparente:=","%5.2f" % (PE_app)," KVA") 
#=============================================================================
#  =============== Factor de potencia promedio ===============================
#==============================================================================

Fp_prom = (PEC_m/PEA_p)*100 # %
print("Factor de potencia promedio:=","%5.2f" % (Fp_prom)," %")

#==============================================================================
#                        Factor de Carga Ciclica CLF
#==============================================================================
Clf=0
for i in range(N_puntos):
    Clf = Power_XF[i]**2 + Clf
    
CLF_AVE = (Clf/N_puntos+1)**0.5

    
P_motor_AVE = mean(Power_XF)
FL          = CLF_AVE/P_motor_AVE
print("Potencia promedio entegada por VSD:=","%5.2f" % (P_motor_AVE),"KVA")
print("Factor de carga ciclico :=","%5.2f" % (FL))

#==============================================================================
# =========================  Corriente Activa del motor =======================
#==============================================================================

Iactiva = Fp_prom*mean(VSDX)/100
print(Iactiva)
# ========================= Corriente REactiva motor ==========================

Ireactiva = ((1-(Fp_prom/100)**2)**0.5)*mean(VSDX)
print(Ireactiva)

#============================ IRMS ============================================
IR=0
for i in range(N_puntos):
    IR = VSDX[i]**2 + IR
    
IRMS = (IR/N_puntos+1)**0.5

print("IRMS:=","%5.2f" % (IRMS))


#=========================== Caraga del motor en % ============================
Cm = (IRMS/Rated_Amp)*100 # Mayor que 60% motor oversized
print("Carga Motor:=","%5.2f" % (Cm),"%")
#==============================================================================
# ========================= Energia consumida Motor ===========================
#==============================================================================
Econs = hpd*dpm*PEC_m
Egen  = hpd*dpm*PEG_g
CESC  = Econs*Pelec
CECC  = (Econs+Egen)*Pelec
CDPM  = PEA_p*Ppde
CTESC = CESC+CDPM
CTECC = CECC+CDPM

print("Energia consumida - motor:=", "%5.2f" % (Econs),"KWh")
print("Energia generada - motor:=", "%5.2f" % (Egen),"KWh")
print("Costo energia consumida - sin credito:=", "%5.2f" % (CESC),"$/mes")
print("Costo energia consumida - con credito:=", "%5.2f" % (CECC),"$/mes")
print("Costo demanda potencia pico - motor :=", "%5.2f"  % (CDPM),"$/mes")
print("Costo Total de energia - sin credito:=", "%5.2f"  % (CTESC),"$/mes")
print("Costo Total de eenrgia - con credito:=", "%5.2f"  % (CTECC),"$/mes")
#==============================================================================
#=============================================================================
#                 Transform the signal in the frequency domain
#=============================================================================
x          = np.linspace (0.0,f/2,NDatos)
freq_data  = fft(PosX)
yf         = NDatos * np.abs (freq_data)
xf         = fftfreq(NDatos,T)[:NDatos//2]

plt.figure(13)
plt.plot(xf,2/NDatos *np.abs(yf[0:NDatos//2]))
plt.title( 'Frequency Domain Signal')
plt.xlabel('Frequency in Hz')
plt.ylabel('Amplitude')
#==============================================================================
plt.figure (14)
plt.plot(Pos_XF,'r',linewidth=2)
plt.twinx()
plt.plot(Power_XF,'b',linewidth=2)
plt.xlabel('numero de muestras')
plt.ylabel('Potencia-KW')
#===================================ZZ=========================================
ysmooth = smooth(DPos,3) # using convolution

plt.figure(15)
plt.plot(ysmooth,'k')
plt.plot(smooth(DPos,30),'r')
plt.xlabel(' Numero de muestras')
plt.ylabel(' Velocity [in/sec]')
#==============================================================================
xl = np.arange(0,NDatos-1,1)
kr = KernelReg(DPos,xl,'c') # regresion

plt.figure(16)
plt.title( 'Velocity Analysis')
plt.plot(xl,DPos,'+')
DPos_P,y_std = kr.fit(xl)
plt.plot(xl,DPos_P,'r')
plt.xlabel('numero de muestras')
plt.ylabel('Velocity fitting [in/sec]')
#==============================================================================

#==============================================================================

#==============================================================================
# ============== Analisis del Torque en el sistema ============================
#==============================================================================
plt.figure(17)
plt.plot(Power_XF)
plt.twinx() # dos ejes
plt.plot(Torq)
#==============================================================================
plt.figure(18)
plt.plot(Tsample,X,'r')
plt.twinx()
plt.plot(Tsample,Pos_XF)
plt.title('Pos-Inclinometro vs Pos Geometria')
plt.xlabel('Tiempo en sec')
plt.ylabel('Polisher Rod - Position[in]')
#==============================================================================
plt.figure(19)
plt.plot(Pos)
plt.twinx()
plt.plot(PuntosX,PosX,'r',linewidth=2)
plt.title('Polisher Rod - Best Stroke Position[in')
plt.xlabel('numero de muestras')
plt.ylabel('Polisher Rod - Position[in]')
#==============================================================================
plt.figure(20)
plt.plot(PuntosX,VSD[N1:N2],'b',linewidth=3)
plt.twinx()
plt.plot(PuntosX,VSD_XF,'r',linewidth=2)
plt.title('Polisher Rod - Best Stroke')
plt.xlabel('numero de muestras')
plt.ylabel('Motor - Ampers[Amp]')

#==============================================================================
plt.figure(21)
plt.plot(Pos,'b',linewidth=3)
plt.twinx()
plt.plot(Prox,'r',linewidth=2)
plt.plot(PuntosX,Pos_XF,'k +',linewidth=2)
plt.title('Stroke Detection')
plt.xlabel('numero de muestras')
plt.ylabel('Motor - Power[KW]')

#==============================================================================
areapos     = trapz(Pos_XF,dx=4) # dx:= espaciado datos
areaX       = trapz(X,dx=4)
Certeza_Pos = (areapos - (areapos - areaX))/areapos
Certeza_Pos = round(Certeza_Pos,3)
print("Certeza_Pos",Certeza_Pos)
#==============================================================================
plt.figure(22)
plt.plot(PowerX,'r')
plt.plot(Power_XF)
#==============================================================================
plt.figure(23)
plt.plot(Tsample,P_rod)
plt.plot(Tsample,P_crank,'r')

#==============================================================================
plt.figure(24)
plt.plot(Tsample,Inv_V)
plt.plot(Tsample,F_rod,'r')
#==============================================================================
plt.figure(25)
plt.plot(X,F_rod)
plt.show()
#======================================ZZ======================================