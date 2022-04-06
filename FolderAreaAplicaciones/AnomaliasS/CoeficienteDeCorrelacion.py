import numpy as np

x = [55, 78, 65, 98, 97, 60, 67, 65, 83, 65, 5, 784, 874875]
y = [45, 68, 65, 95, 17, 60, 67, 85, 83, 65, 35, 78, 8775]
N = len(x)

Promx = 0
Promy = 0
Sum1 = 0
Sum2 = 0
Sum3 = 0

for i in range (0, N):
    Promx = x[i] + Promx
    Promy = y[i] + Promy

Promy = Promy / N
Promx = Promx / N

#Sumatorias 
for i in range (0, N):
    Sum1 = (x[i] - Promx) * (y[i] - Promy) + Sum1
    Sum2 = (x[i] - Promx)**2 + Sum2
    Sum3 = (y[i] - Promy)**2 + Sum3

CoeficienteCorrelacion = Sum1 / (((Sum2)**0.5)*((Sum3)**0.5))
R1 = np.corrcoef(x, y)
print ('Corrcoef = ', CoeficienteCorrelacion)
print (R1 [0][1])
