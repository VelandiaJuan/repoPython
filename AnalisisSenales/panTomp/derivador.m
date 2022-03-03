function y = derivador(vector)
    salida = vector;
    len = length(vector);%tamano del vector
    for i = 0:len
        if i< 5
%             no realice accion 
        else
            salida(i) = (2*vector(i)+vector(i-1)-vector(i-3) - (2*vector(i- 4)))/4;
        end
    end
    y = salida;
end