import numpy as np
from   scipy.spatial   import distance


x = [55, 78, 65, 98, 97, 60, 67, 65, 83, 65, 5, 784, 874875]
y = [45, 68, 65, 95, 17, 60, 67, 85, 83, 65, 35, 78, 8775]
N = len(x)

Eu = distance.euclidean(x,y)
Suma1 = 0
for i in range (0, N):
    Suma1 = (x[i]-y[i])**2 + Suma1

DistEuclidiana = Suma1**0.5
print(DistEuclidiana)
print (Eu)