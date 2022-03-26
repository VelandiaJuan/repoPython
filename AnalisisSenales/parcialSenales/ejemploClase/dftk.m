function [X] = dftk(x,k)
    N = length(x);
    n = 0:N-1;
    exponenciales = exp(-2*pi*1i*k*n/N);
    X = exponenciales*x;
end