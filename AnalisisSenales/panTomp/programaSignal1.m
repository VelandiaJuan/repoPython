% logica principal del programa 
clc
clear all 
close all
data = load('datos.mat');%Cargar los datos
data = data.data;
filtro1 = filtroPasaBaja(data);% aplicar filtro pasa baja a la senal
filtro2 = filtroPasaAlta(filtro1);%aplicar filtro pasa alta
derivador1 = derivador(filtro2);%aplicar derivador
elevador1 = elevarCuadrado(derivador1);
ventana = 22;
panTomp = integralVentana(elevador1,ventana);


% identificacion de los picos de la senal
tamano = size(panTomp);%identificar el tamano del vector
% crear un vector con el tamano del vector de la integral de ventana
x = linspace(0,tamano(2),tamano(2));

% identificacion de picos y su posicion
minPicIden = max(panTomp)/3;%Identifique si su valor supera la 3 parte del pico max
% funcion para encontrar los picos en la senal
[pks,locs] = findpeaks(panTomp,x,'MinPeakProminence',minPicIden);

% encontrar los ciclos en los cuales se repite la frecuencia en terminos de
% numero de muestras
meanCycle = mean(diff(locs));
fMuestreo = 0.0035;%s
frecuenciaPromedio = 60/(meanCycle*fMuestreo);

% diagnostico
if (frecuenciaPromedio <= 100)&(frecuenciaPromedio >= 60)
    mensaje = ['Frecuencia de :',num2str(frecuenciaPromedio),' es normal'];
elseif (frecuenciaPromedio < 60)
    mensaje = ['Frecuencia de :',num2str(frecuenciaPromedio),'  es baja'];
elseif (frecuenciaPromedio > 60)
    mensaje = ['Frecuencia de :',num2str(frecuenciaPromedio),'  es alta'];
else
    mensaje = 'Frecuencia de No IDentificada';
end
disp(mensaje)

% graficas
graficar(data,'Datos Capturados','Muestras','mV');
graficar(filtro1,'Señal Filtro PasaBaja','Muestras','mV');
graficar(filtro2,'Señal Filtro PasaAlta','Muestras','mV');
graficar(derivador1,'Señal con Derivador','Muestras','mV');
graficar(elevador1,'Señal Elevada al Cuadrado','Muestras','mV^2');
graficar(panTomp,'Señal con la Integral de Ventana','Muestras','Adimesional');
findpeaks(panTomp,'MinPeakProminence',minPicIden)