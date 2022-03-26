function [Ysamp]=Sintetizador(xn,b)
%% LA MISIÓN
% Diseña una función que tenga tres entradas, xn(señal, en este caso será una canción (usa 'songname' para escribirla), y b frecuencia de muestreo final. 
% Obtendrás una repetición de la canción con la nueva frecuencia de muestreo. 


% Primera Parte %
[y,Fs]=audioread(xn);  % Como tenemos un archivo de audio, puedes subir la canción que quieras, de aquí sacamos las Fs y la señal 
Fs2=b; % El muestreo final que se utilizará 

t=1:length(y(:,2));  % Queremos trazar la señal discretizada de un canal de la canción. 
figure(1);  
subplot(3,1,1)  % Nosotros graficamos usando stem
plot(t,y); 
xlabel(' Real Time ');
ylabel('Audio Signal ');


t=1:length(y(:,2));  % Queremos trazar la señal discretizada de un canal de la canción. 
Td=t/Fs; % Teniendo nuestras Fs, podemos tener nuestro tiempo Discreto 
figure(1);  
subplot(3,1,2)  % Nosotros graficamos usando stem
stem(Td,y); 
xlabel(' Discrete Time ');
ylabel('Audio Signal Sampled by their frecuency');


%%%%%%%%%%%%% DIEZMADO E INTERPOLACIÓN   %%%%%%%%  

[L,M]=rat(Fs2/Fs); % Usando RAT (que devuelve la aproximación racional de la diferencia entre la frecuencia anterior y la nueva frecuencia) podemos calcular 
% los L y M (nuestros coeficientes de Decimación e Interpolación). ¡ESTOS SON MUY IMPORTANTES! 

%%% NUESTRO MAGNÍFICO FILTRO (Filtro LOW PASS antialiasing con ventana KAISER)

f=(Fs/2)*min(1/L,1/M); % Definimos nuestro Fc (Recuerda que por nyquist nuestro Fc debería ser Fs/2, pero como tenemos los parámetros L y M, necesitamos calcular la 
% de frecuencia que podría pasar nuestro filtro) tratamos de establecer la banda del filtro
% a un 90 % y 110% de la frecuencia de corte.

MFK=designfilt('lowpassfir','PassbandFrequency',0.9*f,'StopbandFrequency',1.1*f,'PassbandRipple',5,'StopbandAttenuation',40, 'DesignMethod','kaiserwin','SampleRate',Fs);

% No usamos YULEWALK (muy pocho), creamos un filtro paso bajo antialiasing usando una ventana Kaiser. Nuestra ondulación de la banda de paso fue de 5[dB] y una atenuación de la banda de parada de 40 dB, finalmente, ajustamos la ganancia de la banda de paso a L.

h=L*tf(MFK);  % Convertimos este filtro en una función de transferencia usando (tf), y multiplicado por L (Añadiendo los ceros deseados) 

Ysamp=upfirdn(y,h,L,M); % Aplicado el filtro FIR upsample, con, y(canción), h(función de transferencia),L(Decimación), M(Interpolación)
delay = floor(((filtord(MFK)-1)/2-(L-1))/L); % Cuando usamos este Ysamp, hay un retraso introducido por el filtro, intentamos corregir este retraso
Ysamp=Ysamp(delay+1:end);
t_res = (0:(length(Ysamp)-1))/Fs2; % Tenemos una nueva hora de muestreo 


% GRÁFICO DE LA SEÑAL REMUESTREADA 
figure(1); 
subplot(3,1,3)
plot(t_res,Ysamp,'o'); 
legend('Resampled')
xlabel(' Discrete Time ');
ylabel('Audio Signal RESAMPLED ');