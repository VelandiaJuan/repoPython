import numpy as np
from   scipy.integrate   import trapz


x = [55, 78, 65, 8, 97, 60, 67, 65, 83, 65, 5, 784, 87481]
y = [45, 68, 65, 95, 17, 60, 67, 85, 83, 65, 35, 78, 8775]
N = len(x)

trapz1 = 0
for i in range (0, N-1):
    trapz1 = (y[i]*(x[i+1]-x[i]) + ((y[i+1]-y[i])/2)*(x[i+1]-x[i])) + trapz1

Xd = np.trapz (y, x)
print ('a mano:' ,trapz1)
print ('funcion: ', Xd)