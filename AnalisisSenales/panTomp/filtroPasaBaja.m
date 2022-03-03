function y = filtroPasaBaja(vector)
    salida = vector;
    len = length(vector);%tamano del vector
    for i = 1:len
        if i< 13
%             no realice accion 
        else
            salida(i) = 2*salida(i-1) - salida(i-2) + vector(i) - 2*vector(i-6) + vector(i-12) ;
        end
    end
    y = salida;
end