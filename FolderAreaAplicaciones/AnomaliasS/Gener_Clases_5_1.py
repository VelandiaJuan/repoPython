# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 14:50:37 2022

@author: jose.contreras
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 15:42:26 2022

@author: jose.contreras
"""

import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt
import openpyxl
import math
import scipy.signal
from   scipy.stats     import skew
from   scipy.stats     import kurtosis
from   scipy.spatial   import distance
from   scipy.integrate import trapz
from   statistics      import mode
# =============================================================================
#               RUTINA PARA ENCONTRAR EL VECTOR ANOMALIAS
# =============================================================================
def vec_desc (Pos_F, Load_F):
    #==============================================================================
    #           Calculo de parametros estadisticos de los vectores
    #==============================================================================
    nlineas = len (Pos_F)
    Ske_P   = skew(Pos_F,bias=False)
    Ske_L   = skew(Load_F,bias=False)
    Kur_P   = kurtosis(Pos_F,fisher=False)        # Fisher Pearson correlation of skew
    Kur_L   = kurtosis(Load_F,fisher=False)
    Ske_P   = round(Ske_P,4)
    Ske_L   = round(Ske_L,4)
    Kur_P   = round(Kur_P,4)
    Kur_L   = round(Kur_L,4)
    Dist_PL = distance.euclidean(Pos_F,Load_F)
    Dist_PL = round(Dist_PL,2)
    Tpos    = range(0,nlineas,1)
    I_trap  = trapz(Load_F,Pos_F)
    I_trap  = round(I_trap,3)
    # =============================================================================
    # Calculo del centroide del dinagrama de superficie y de fondo
    # =============================================================================
    
    Xo = round(np.sum(Pos_F*I_trap)/I_trap,3)
    Yo = round(np.sum(Load_F*I_trap)/I_trap,3)
    Xc = round(np.sum(Pos_F)/nlineas,3)
    Yc = round(np.sum(Load_F)/nlineas,3)
    
    # =============================================================================
    #                         Calculo de la Media Aritmetica
    # =============================================================================
    Media_A_P = sum(Pos_F)/nlineas
    Media_A_L = sum(Load_F)/nlineas
    # =============================================================================
    #                 Calculo de la desviación respecto a la media
    # =============================================================================
    Des_M_P = abs(Pos_F - Media_A_P)
    Des_M_L = abs(Load_F - Media_A_L)
    Min_D_L = min(Des_M_L)
    Max_D_S = max(Des_M_P)
    # =============================================================================
    #                      Calculo de la varianza: Dispersion
    # =============================================================================
    Var_z_P = sum(Pos_F - Media_A_P)**2/nlineas
    Var_z_L = sum(Load_F - Media_A_L)**2/nlineas
    # =============================================================================
    #                          Desviacion tipica
    # =============================================================================
    Desv_T_P = (Var_z_P)**.5
    Desv_T_L = (Var_z_L)**.5
    # =============================================================================
    #            Calculo de la Covarianza: Sentido de la correlación
    # =============================================================================
    Co_var = sum((Pos_F-Media_A_P)*(Load_F-Media_A_L))/nlineas
    Co     = np.cov(Pos_F,Load_F)
    # =============================================================================
    #                  Calculo del mode del vector: Numero de repetición de valores
    # =============================================================================
    # Mode_P = mode(Pos_F)
    # Mode_L = mode(Load_F)
    # =============================================================================
    #  Calculo de la correlación entre Pos y Load en la bomba
    # =============================================================================
    CorrePL = np.corrcoef(Pos_F,Load_F)
    D_corre = CorrePL[0][1]
    # =============================================================================
    # Calculo de las distancias del centroide a los max y min de la carga
    # =============================================================================
    [ML,MP]  = Dist_M_Centro(Load_F,Pos_F) # Calculo de los indices de la max. carga
    [mL,mP]  = Dist_m_Centro(Load_F,Pos_F)
    [mPX,mLX]= Dist_m_Centro_X(Load_F,Pos_F)  # Distancia centro al minX
    [MPX,MLX]= Dist_M_Centro_X(Load_F,Pos_F)
    Dist_CML = ((Xc+MP)**2 + (Yc+ML)**2)**.5 # Distancia del centroide max- carga
    Dist_CmL = ((Xc+mP)**2 + (Yc+mL)**2)**.5 # Dist del centroide min carga
    Dist_CmX = ((Xc+mPX)**2 + (Yc+mLX)**2)**.5
    Dist_CMX = ((Xc+MPX)**2 + (Yc+MLX)**2)**.5
    
    # =============================================================================
    #  Calculo de normalizacion del Vector Descriptor del Dinagrama (VDD)
    # =============================================================================
    
    VDD = []
    VDD = [Xc,Yc,I_trap,Dist_PL,Kur_P,Kur_L,Ske_L,Ske_P,Co_var,Media_A_P,Media_A_L,D_corre,
           Dist_CML,Dist_CmL, Dist_CmX, Dist_CMX]
    VDD_n  = np.array(VDD)
    npp    = len(VDD_n)
    max_V  = max(VDD_n)
    min_V  = min(VDD_n)
    Norm_D = np.round((VDD_n-min_V)/(max_V-min_V),3)
    
    return VDD_n
# =============================================================================
#                 RUTINA PARA LECTURA DE DATOS DE HOJA DE EXCEL
# =============================================================================


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

# =============================================================================
#                                    rutinas
# =============================================================================
def Dist_M_Centro(Load_F,Pos_F):# Distancias del centroide al max  de la carga
    ML   = np.max(Load_F)
    ind1 = np.argmax(Load_F)
    MP   = Pos_F[ind1]
    return ML,MP

def Dist_m_Centro(Load_F,Pos_F):# Distancias del centroide al min de la carga
    mL   = np.min(Load_F)
    ind2 = np.argmin(Load_F)
    mP   = Pos_F[ind2]
    return mL,mP

def Dist_m_Centro_X(Load_F,Pos_F):# Distancias del centroide al min de la carga
    mPX   = np.min(Pos_F)
    ind3  = np.argmin(Pos_F)
    mLX   = Load_F[ind3]
    return mPX,mLX

def Dist_M_Centro_X(Load_F,Pos_F):# Distancias del centroide al min de la carga
    MPX   = np.max(Pos_F)
    ind4  = np.argmax(Pos_F)
    MLX   = Load_F[ind4]
    return MPX,MLX

def Norma_Dyna_pos(Pos):# Normalizar dinagrama

    np = len(Pos)   # Longitud muestras
    mX = min(Pos)   # minimo de X
    MX = max(Pos)   # maximo de X
    Pos_Sn  = (Pos-mX)/(MX-mX) # X normalizado
    return Pos_Sn

def Norma_Dyna_load(Load):# Normalizar dinagrama
    np = len(Load)   # Longitud muestras
    mL = min(Load)  # Min de la carga
    ML = max(Load)
    Load_Sn = (Load-mL)/(ML-mL)
    return Load_Sn

def Class_Dyna(matriz_n):# Calculo de la coordenada cartesiana de la clase 
    nclases = len(matriz_n)
    clases = np.zeros ((nclases, 2))
    np1  = len(matriz_n[0][:])
    nV1 = np1//2
    for i in range (0, nclases):
        clases[i][0]  = sum(matriz_n[i][0:nV1]/nV1)
        clases[i][1]  = sum(matriz_n[i][nV1:np1]/nV1)
        
    return clases

def Class_Dyna1(V_c):# Calculo de la coordenada cartesiana de la clase 
    np3  = len(V_c)   # Numero de clases
    nV1 = np3//2
    V1  = sum(V_c[0:nV1]/nV1)
    V2  = sum(V_c[nV1:np3]/nV1)
    return V1,V2

# =============================================================================
#              LECTURA DE DATOS DE HOJA DE EXCEL
# =============================================================================
path = "DINAGRAMAS_VALIDACION.xls"
Leak_sv      = LeerExcel(path,0,2)
Leak_tub     = LeerExcel(path,1,2)
Asphaltenes  = LeerExcel(path,2,2)
Tagging_down = LeerExcel(path,3,2)
Gunk_pump    = LeerExcel(path,4,2)
Deep_well    = LeerExcel(path,5,2)
Gas_inter    = LeerExcel(path,6,2)
Average      = LeerExcel(path,7,2)
Leak_TV      = LeerExcel(path,8,2)
Fluid_pound  = LeerExcel(path,9,2)
Flumping     = LeerExcel(path,10,2)

plt.figure
plt.plot(Leak_TV.posicionF, Leak_TV.cargaF)
plt.show()

Leak_svC = vec_desc (Leak_sv.posicionF, Leak_sv.cargaF)
Leak_tubC = vec_desc (Leak_tub.posicionF, Leak_tub.cargaF)
AsphaltenesC = vec_desc (Asphaltenes.posicionF, Asphaltenes.cargaF)
Tagging_downC = vec_desc (Tagging_down.posicionF, Tagging_down.cargaF)
Gunk_pumpC = vec_desc (Gunk_pump.posicionF, Gunk_pump.cargaF)
Deep_wellC = vec_desc (Deep_well.posicionF, Deep_well.cargaF)
Gas_interC = vec_desc (Gas_inter.posicionF, Gas_inter.cargaF)
AverageC = vec_desc (Average.posicionF, Average.cargaF)
Leak_TVC = vec_desc (Leak_TV.posicionF, Leak_TV.cargaF)
Fluid_poundC = vec_desc (Fluid_pound.posicionF, Fluid_pound.cargaF)
FlumpingC = vec_desc (Flumping.posicionF, Flumping.cargaF)
# =============================================================================
# 
# =============================================================================
path2 = "DINAGRAMAS_ANOMALIAS_SKYPE.xls"

AirLock                = LeerExcel(path2,0,0)
BrokenRod              = LeerExcel(path2,1,0)
GasExistence           = LeerExcel(path2,2,0)
InletValveDelayClosing = LeerExcel(path2,3,0)
InletValveLeakage      = LeerExcel(path2,4,0)
OutletSuddenUnloading  = LeerExcel(path2,5,0)
OutletValveLeakage     = LeerExcel(path2,6,0)
ParrafinWax            = LeerExcel(path2,7,0)
PistonSticking         = LeerExcel(path2,8,0)
PlungerFixedValveCollision = LeerExcel(path2,9,0)
PlungerGuideRingCollision = LeerExcel(path2,10,0)
ProperWork2            = LeerExcel(path2,11,0)
RodExcesiveVibration   = LeerExcel(path2,12,0)
SandProblem            = LeerExcel(path2,13,0)
SmallPlunger           = LeerExcel(path2,14,0)
ThickOil               = LeerExcel(path2,15,0)
ValveLeakage           = LeerExcel(path2,16,0)



plt.figure
plt.plot(PistonSticking.posicionF, PistonSticking.cargaF)
plt.show()

AirLockC = vec_desc (AirLock.posicionF, AirLock.cargaF)
BrokenRodC = vec_desc (BrokenRod.posicionF, BrokenRod.cargaF)
GasExistenceC = vec_desc (GasExistence.posicionF, GasExistence.cargaF)
InletValveDelayClosingC = vec_desc (InletValveDelayClosing.posicionF, InletValveDelayClosing.cargaF)
InletValveLeakageC = vec_desc (InletValveLeakage.posicionF, InletValveLeakage.cargaF)
OutletSuddenUnloadingC = vec_desc (OutletSuddenUnloading.posicionF, OutletSuddenUnloading.cargaF)
OutletValveLeakageC = vec_desc (OutletValveLeakage.posicionF, OutletValveLeakage.cargaF)
ParrafinWaxC = vec_desc (ParrafinWax.posicionF, ParrafinWax.cargaF)
PistonStickingC = vec_desc (PistonSticking.posicionF, PistonSticking.cargaF)
PlungerFixedValveCollisionC = vec_desc (PlungerFixedValveCollision.posicionF, PlungerFixedValveCollision.cargaF)
PlungerGuideRingCollisionC = vec_desc (PlungerGuideRingCollision.posicionF, PlungerGuideRingCollision.cargaF)
ProperWork2C = vec_desc (ProperWork2.posicionF, ProperWork2.cargaF)
RodExcesiveVibrationC = vec_desc (RodExcesiveVibration.posicionF, RodExcesiveVibration.cargaF)
SandProblemC = vec_desc (SandProblem.posicionF, SandProblem.cargaF)
SmallPlungerC = vec_desc (SmallPlunger.posicionF, SmallPlunger.cargaF)
ThickOilC = vec_desc (ThickOil.posicionF, ThickOil.cargaF)
ValveLeakageC = vec_desc (ValveLeakage.posicionF, ValveLeakage.cargaF)

# =============================================================================
# 
# =============================================================================
matriz=[]

matriz=[Leak_svC, Leak_tubC, AsphaltenesC, Tagging_downC, Gunk_pumpC, Deep_wellC, Gas_interC, AverageC, Leak_TVC, Fluid_poundC, FlumpingC]#, AirLockC, BrokenRodC, GasExistenceC, InletValveDelayClosingC, InletValveLeakageC, OutletSuddenUnloadingC, OutletValveLeakageC, ParrafinWaxC, PistonStickingC, PlungerFixedValveCollisionC, PlungerGuideRingCollisionC, ProperWork2C, RodExcesiveVibrationC, SandProblemC, SmallPlungerC, ThickOilC, ValveLeakageC]

# matriz=[Leak_svC, Leak_tubC, AsphaltenesC, Tagging_downC, Gunk_pumpC, Deep_wellC, Gas_interC, AverageC, Leak_TVC, Fluid_poundC, FlumpingC, AirLockC, BrokenRodC, GasExistenceC, InletValveDelayClosingC, InletValveLeakageC, OutletSuddenUnloadingC, OutletValveLeakageC, ParrafinWaxC, PistonStickingC, PlungerFixedValveCollisionC, PlungerGuideRingCollisionC, ProperWork2C, RodExcesiveVibrationC, SandProblemC, SmallPlungerC, ThickOilC, ValveLeakageC]

matriz1=['Leak_sv', 'Leak_tub', 'Asphaltenes', 'Tagging_down', 'Gunk_pump', 'Deep_well', 'Gas_inter', 'Average', 'Leak_TV', 'Fluid_pound', 'Flumping']#, 'AirLock', 'BrokenRod', 'GasExistence', 'InletValveDelayClosing', 'InletValveLeakage', 'OutletSuddenUnloading', 'OutletValveLeakage', 'ParrafinWax', 'PistonSticking', 'PlungerFixedValveCollision', 'PlungerGuideRingCollision', 'ProperWork2', 'RodExcesiveVibration', 'SandProblem', 'SmallPlunger', 'ThickOil', 'ValveLeakage']

# matriz1=['Leak_sv', 'Leak_tub', 'Asphaltenes', 'Tagging_down', 'Gunk_pump', 'Deep_well', 'Gas_inter', 'Average', 'Leak_TV', 'Fluid_pound', 'Flumping', 'AirLock', 'BrokenRod', 'GasExistence', 'InletValveDelayClosing', 'InletValveLeakage', 'OutletSuddenUnloading', 'OutletValveLeakage', 'ParrafinWax', 'PistonSticking', 'PlungerFixedValveCollision', 'PlungerGuideRingCollision', 'ProperWork2', 'RodExcesiveVibration', 'SandProblem', 'SmallPlunger', 'ThickOil', 'ValveLeakage']

matriz2=['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11']#, 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28']

# matriz2=['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28']

matriz3=['C1: Leak_sv', 'C2: Leak_tub', 'C3: Asphaltenes', 'C4: Tagging_down', 'C5: Gunk_pump', 'C6: Deep_well', 'C7: Gas_inter', 'C8: Average', 'C9: Leak_TV', 'C10: Fluid_pound']#, 'C11: Flumping', 'C12: AirLock', 'C13: BrokenRod', 'C14: GasExistence', 'C15: InletValveDelayClosing', 'C16: InletValveLeakage', 'C17: OutletSuddenUnloading', 'C18: OutletValveLeakage', 'C19: ParrafinWax', 'C20: PistonSticking', 'C21: PlungerFixedValveCollision', 'C22: PlungerGuideRingCollision', 'C23: ProperWork2', 'C24: RodExcesiveVibration', 'C25: SandProblem', 'C26: SmallPlunger', 'C27: ThickOil', 'C28: ValveLeakage']

# matriz3=['C1: Leak_sv', 'C2: Leak_tub', 'C3: Asphaltenes', 'C4: Tagging_down', 'C5: Gunk_pump', 'C6: Deep_well', 'C7: Gas_inter', 'C8: Average', 'C9: Leak_TV', 'C10: Fluid_pound', 'C11: Flumping', 'C12: AirLock', 'C13: BrokenRod', 'C14: GasExistence', 'C15: InletValveDelayClosing', 'C16: InletValveLeakage', 'C17: OutletSuddenUnloading', 'C18: OutletValveLeakage', 'C19: ParrafinWax', 'C20: PistonSticking', 'C21: PlungerFixedValveCollision', 'C22: PlungerGuideRingCollision', 'C23: ProperWork2', 'C24: RodExcesiveVibration', 'C25: SandProblem', 'C26: SmallPlunger', 'C27: ThickOil', 'C28: ValveLeakage']

matriz_n  = np.array(matriz)

clasesanomalias = Class_Dyna(matriz_n)

promedio = [np.mean(clasesanomalias[:,0]), np.mean(clasesanomalias[:, 1])]

# =============================================================================
#                        LECTURA DE UN NUEVO DINAGRAMA
# =============================================================================
path1 = 'DINAGRAMAS_VALIDACION_1.xls'
Tuberia_no_anclada = LeerExcel(path1, 12,2)
Tuberia_no_ancladaC = vec_desc (Tuberia_no_anclada.posicionF, Tuberia_no_anclada.cargaF)

[V1,V2] = Class_Dyna1(Tuberia_no_ancladaC)


# =============================================================================
#      CALCULO DEL PORENTAJE DE CERTEZA CON EL TERCER PUNTO MAS CERCANO
# =============================================================================
npuntos = len(clasesanomalias)
V_dist = np.zeros((len(clasesanomalias), 2))
V_dist_ord = np.zeros((len(clasesanomalias), 2))
certeza = np.zeros(len(clasesanomalias))

for i in range (0, npuntos):
    V_dist[i][0] = ((clasesanomalias[i][0] - V1)**2 + (clasesanomalias[i][1] - V2)**2)**0.5
    V_dist[i][1] = i


mindist = 99999
minind = 99999
for k in range (0, npuntos):
    for i in range (0, npuntos):
        if k == 0:
            if V_dist[i][0] < mindist:
                mindist = V_dist[i][0]
                minind = V_dist[i][1]
        if k != 0:
            if V_dist[i][0] >= V_dist_ord[k-1][0] and V_dist_ord[k-1][1] != V_dist[i][1]:
                if V_dist[i][0] < mindist:
                    mindist = V_dist[i][0]
                    minind = V_dist[i][1]
    V_dist_ord[k][0] = mindist
    V_dist_ord[k][1] = minind
    mindist = 99999
    minind = 99999
    

Aux10=0
Radio = V_dist_ord[3][0]
for i in range (0, 3):
    certeza[i] = (1-(V_dist_ord[i][0]/Radio))*100
    
for i in range (0, 3):
    if Aux10==0:
        if certeza[i] == 100:
            print('Se tiene un porcentaje de certeza del 100% de ser de la clase ', int(V_dist_ord[i][1]), matriz1[int(V_dist_ord[i][1])], '.')
            Aux10=1
for i in range (0, 3):
    if Aux10 != 1:
        if V_dist_ord[i][0]<Radio:
            print('Se tiene un porcentaje de certeza del ', round(certeza[i],1) ,'% de ser de la clase ', int(V_dist_ord[i][1])+1, matriz1[int(V_dist_ord[i][1])], '.')

    
# =============================================================================
#             CÁLCULO DE LA CERTEZA  CON EL PUNTO PROMEDIO
# =============================================================================
# npuntos = len(clasesanomalias)
# V_dist = np.zeros((len(clasesanomalias), 2))
# V_dist_ord = np.zeros((len(clasesanomalias), 2))
# certeza = np.zeros(len(clasesanomalias))

# for i in range (0, npuntos):
#     V_dist[i][0] = ((clasesanomalias[i][0] - V1)**2 + (clasesanomalias[i][1] - V2)**2)**0.5
#     V_dist[i][1] = i

# Radio = ((promedio[0]-V1)**2 + (promedio[1] - V2)**2)**0.5
# Aux10=0
# for i in range (0, len(clasesanomalias)):
#     certeza[i] = (1-(V_dist[i][0]/Radio))*100
    
# for i in range (0, len(clasesanomalias)):
#     if Aux10==0:
#         if certeza[i] == 100:
#             print('Se tiene un porcentaje de certeza del 100% de ser de la clase ', i, matriz1[i], '.')
#             Aux10=1
# for i in range (0, len(clasesanomalias)):
#     if Aux10 != 1:
#         if V_dist[i][0]<Radio:
#             print('Se tiene un porcentaje de certeza del ', round(certeza[i],1) ,'% de ser de la clase ', i, matriz1[i], '.')

# =============================================================================
#                  RADIO DE ATRACCIÓN PARA GRÁFICA 
# =============================================================================
num_segmentos = 100
rad = Radio
cx = V1
cy = V2

angulo = np.linspace(0, 2*np.pi, num_segmentos+1)
x = rad * np.cos(angulo) + cx
y = rad * np.sin(angulo) + cy
# =============================================================================
#                         PLOTEO DE LAS FIGURAS
# =============================================================================

plt.subplot (2,2,1)
plt.title ("Dinagrama de fondo actual normalizado")
plt.xlabel("Posición normalizada")
plt.ylabel("Carga normalizada")
plt.plot (Tuberia_no_anclada.posicionF,Tuberia_no_anclada.cargaF)
plt.text(0.5, -0.4, str(round(certeza[0],1))+' % '+ ' de certeza que es un pozo de clase '+ str(int(V_dist_ord[i][1])+1)+', '+ str(matriz1[int(V_dist_ord[0][1])]), fontsize=15, color='black')


plt.subplot (2,2,2)
for i in range (0,len (clasesanomalias)):
    plt.plot(float(clasesanomalias[i][0]), float(clasesanomalias[i][1]), marker="o", color="red")
for i, label in enumerate(matriz2):
    plt.annotate(label, (float(clasesanomalias[i][0])*1.005, float(clasesanomalias[i][1]*1.005)), fontsize=7)
plt.plot (promedio[0], promedio[1], marker="o", color="black")
plt.annotate('CC', (promedio[0]*1.005, promedio[1]*1.005), fontsize=7)
plt.plot (V1, V2, markersize= 7, marker=".", color="black")
plt.annotate('DP', (V1*1.005, V2*1.005), fontsize=7)
plt.plot(x, y, linestyle='--',color="blue")
plt.legend(matriz3, bbox_to_anchor=(1.05, 1.0), loc='upper left', fontsize=7)
plt.show()