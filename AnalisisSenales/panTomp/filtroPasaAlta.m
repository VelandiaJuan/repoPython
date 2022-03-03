function y = filtroPasaAlta(vector)
    salida = vector;
    len = length(vector);%tamano del vector
    for i = 0:len
        if i< 33
%             no realice accion 
        else
            salida(i) = salida(i-1)-(vector(i)/32) + vector(i-16) - vector(i-17) + (vector(i-32)/32);
        end
    end
    y = salida;
end