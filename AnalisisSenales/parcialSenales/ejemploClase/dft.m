function [X] = dft(x)
    N = length(x);
    n = 0:N-1;
    k = 0:N-1;
    kn = n'*k;
    X = exp(-2*pi*1i*kn/N)*x;
end