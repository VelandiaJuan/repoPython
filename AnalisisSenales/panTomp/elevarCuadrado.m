function y = elevarCuadrado(vector)
    salida = vector;
    len = length(vector);%tamano del vector
    for i = 1:len
        salida(i) = vector(i)^2;
    end
    y = salida;
end