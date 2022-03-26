% programa
x = [1 -1 2 5 3 2];
k = 1;%componente continua de la señal
X = dftk(x',1);
% validacion de la transformada de fourier
fft(x)
%ecuacion con todos los coeficientes
dft(x')