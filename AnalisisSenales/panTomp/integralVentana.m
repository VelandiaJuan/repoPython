function y = integralVentana(vector,ventana)
    salida = vector;
    tamano = round(ventana/2);
    len = length(vector);%tamano del vector
    
    for i = 1:len
        acomulador = 0;
        if(i > len -tamano)
            break
        end
        if(i > tamano)   
            for j = (i - tamano):(i+tamano)
                acomulador = acomulador + vector(j);
            end
            salida(i) = acomulador/(tamano+1);
        end 
    end
    y = salida;
end