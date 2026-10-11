function B = shifted_legendre_fun_vector(N,eta,x)

syms k;

for i = 0:N-1
    sigma = ((-1)^(i+k)*gamma(i+k+1)*(x^k))/( ...
        gamma(i-k+1)*(eta^k)*(gamma(k+1)^2));

    A(i+1,1) = symsum(sigma,k,0,i);
end

B = A;

end
