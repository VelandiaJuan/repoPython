function [Ysamp]=Sintetizador(xn,b)
%% LA MISI�N
% Dise�a una funci�n que tenga tres entradas, xn(se�al, en este caso ser� una canci�n (usa 'songname' para escribirla), y b frecuencia de muestreo final. 
% Obtendr�s una repetici�n de la canci�n con la nueva frecuencia de muestreo. 


% Primera Parte %
[y,Fs]=audioread(xn);  % Como tenemos un archivo de audio, puedes subir la canci�n que quieras, de aqu� sacamos las Fs y la se�al 
Fs2=b; % El muestreo final que se utilizar� 

t=1:length(y(:,2));  % Queremos trazar la se�al discretizada de un canal de la canci�n. 
figure(1);  
subplot(3,1,1)  % Nosotros graficamos usando stem
plot(t,y); 
xlabel(' Real Time ');
ylabel('Audio Signal ');


t=1:length(y(:,2));  % Queremos trazar la se�al discretizada de un canal de la canci�n. 
Td=t/Fs; % Teniendo nuestras Fs, podemos tener nuestro tiempo Discreto 
figure(1);  
subplot(3,1,2)  % Nosotros graficamos usando stem
stem(Td,y); 
xlabel(' Discrete Time ');
ylabel('Audio Signal Sampled by their frecuency');


%%%%%%%%%%%%% DIEZMADO E INTERPOLACI�N   %%%%%%%%  

[L,M]=rat(Fs2/Fs); % Usando RAT (que devuelve la aproximaci�n racional de la diferencia entre la frecuencia anterior y la nueva frecuencia) podemos calcular 
% los L y M (nuestros coeficientes de Decimaci�n e Interpolaci�n). �ESTOS SON MUY IMPORTANTES! 

%%% NUESTRO MAGN�FICO FILTRO (Filtro LOW PASS antialiasing con ventana KAISER)

f=(Fs/2)*min(1/L,1/M); % Definimos nuestro Fc (Recuerda que por nyquist nuestro Fc deber�a ser Fs/2, pero como tenemos los par�metros L y M, necesitamos calcular la 
% de frecuencia que podr�a pasar nuestro filtro) tratamos de establecer la banda del filtro
% a un 90 % y 110% de la frecuencia de corte.

MFK=designfilt('lowpassfir','PassbandFrequency',0.9*f,'StopbandFrequency',1.1*f,'PassbandRipple',5,'StopbandAttenuation',40, 'DesignMethod','kaiserwin','SampleRate',Fs);

% No usamos YULEWALK (muy pocho), creamos un filtro paso bajo antialiasing usando una ventana Kaiser. Nuestra ondulaci�n de la banda de paso fue de 5[dB] y una atenuaci�n de la banda de parada de 40 dB, finalmente, ajustamos la ganancia de la banda de paso a L.

h=L*tf(MFK);  % Convertimos este filtro en una funci�n de transferencia usando (tf), y multiplicado por L (A�adiendo los ceros deseados) 

Ysamp=upfirdn(y,h,L,M); % Aplicado el filtro FIR upsample, con, y(canci�n), h(funci�n de transferencia),L(Decimaci�n), M(Interpolaci�n)
delay = floor(((filtord(MFK)-1)/2-(L-1))/L); % Cuando usamos este Ysamp, hay un retraso introducido por el filtro, intentamos corregir este retraso
Ysamp=Ysamp(delay+1:end);
t_res = (0:(length(Ysamp)-1))/Fs2; % Tenemos una nueva hora de muestreo 


% GR�FICO DE LA SE�AL REMUESTREADA 
figure(1); 
subplot(3,1,3)
plot(t_res,Ysamp,'o'); 
legend('Resampled')
xlabel(' Discrete Time ');
ylabel('Audio Signal RESAMPLED ');